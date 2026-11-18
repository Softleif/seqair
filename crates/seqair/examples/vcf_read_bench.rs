//! VCF/BCF read throughput per access pattern: seqair against rust-htslib and
//! noodles where they have the pattern. A measurement tool, not a test.
//!
//! `vcf_read_bench <lib> <pattern> <file> [threads] [reps]`, keys from the
//! environment: `INFO1=DP:int INFO2=CPG:flag FMT=DP SAMPLE=name` (types: int
//! ints float floats flag). Prints a TSV line with the best-of-`reps`
//! ns/record and a checksum that ties the libraries' reads together.
//!
//! Patterns (seqair): `count` `locus` `site` `info0` (INFO copied, never
//! indexed) `info1` `site_info1` `filter_info2` (two INFO keys on PASS
//! records) `sample` (`.sample(SAMPLE)`) `only_sample` `column` (GT of every
//! sample) `full` `all` (every value, decoded dynamically).
#![allow(
    clippy::unwrap_used,
    clippy::print_stdout,
    clippy::indexing_slicing,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::arithmetic_side_effects,
    reason = "benchmark tool"
)]

use rust_htslib::bcf::{self as hts, Read as _};
use seqair::vcf::{
    DynValue, Fields, FilterStatus, FormatGt, FormatInt, FullRecord, GenotypeRef, Info, InfoFlag,
    InfoFloat, InfoFloats, InfoInt, InfoInts, Locus, OnlySample, PlanContext, PlanError, RawRecord,
    Sample, Samples, Sections, Site, VcfFile, VcfPos, VcfReadError,
};
use std::hint::black_box;
use std::num::NonZeroUsize;
use std::time::Instant;

#[derive(Clone)]
enum IKey {
    Int(InfoInt),
    Ints(InfoInts),
    Float(InfoFloat),
    Floats(InfoFloats),
    Flag(InfoFlag),
}

fn ikey(file: &VcfFile, var: &str) -> Option<IKey> {
    let spec = std::env::var(var).ok()?;
    let (name, typ) = spec.split_once(':').unwrap();
    let h = file.header();
    Some(match typ {
        "int" => IKey::Int(h.info_key(name).unwrap()),
        "ints" => IKey::Ints(h.info_key(name).unwrap()),
        "float" => IKey::Float(h.info_key(name).unwrap()),
        "floats" => IKey::Floats(h.info_key(name).unwrap()),
        "flag" => IKey::Flag(h.info_key(name).unwrap()),
        t => panic!("type {t}"),
    })
}

fn touch_info(info: &Info, key: &IKey) -> Result<u64, VcfReadError> {
    Ok(match key {
        IKey::Int(k) => info.get(k)?.map_or(0, |v| v as u64),
        IKey::Ints(k) => info.get(k)?.map_or(0, |v| v.map(|x| x as u64).sum()),
        IKey::Float(k) => info.get(k)?.map_or(0, |v| u64::from(v.to_bits())),
        IKey::Floats(k) => info.get(k)?.map_or(0, |v| v.map(|x| u64::from(x.to_bits())).sum()),
        IKey::Flag(k) => u64::from(info.get(k)?),
    })
}

fn start(v: VcfPos) -> u64 {
    match v {
        VcfPos::Telomere => 0,
        VcfPos::At(p) => u64::from(p.as_u32()),
    }
}

fn pass(s: &Site) -> bool {
    matches!(s.filters().status(), FilterStatus::Pass)
}

fn gt_sum(g: Option<GenotypeRef<'_>>) -> u64 {
    g.map_or(0, |g| g.alleles().map(|a| a.allele.map_or(7, u64::from)).sum())
}

fn dyn_sum(v: DynValue<'_>) -> u64 {
    match v {
        DynValue::Int(v) => v.map(|x| x.map_or(1, |x| x as u64)).sum(),
        DynValue::Float(v) => v.map(|x| x.map_or(1, |x| u64::from(x.to_bits()))).sum(),
        DynValue::Flag => 1,
        DynValue::Str(s) => s.map_or(0, |s| s.len() as u64),
        DynValue::Gt(g) => gt_sum(g),
    }
}

/// Fixed fields only, nothing kept: the floor of a seqair read loop.
#[derive(Clone)]
struct Count;

impl Fields for Count {
    type Plan = ();

    fn plan(_: &PlanContext<'_>) -> Result<(), PlanError> {
        Ok(())
    }

    fn needs((): &()) -> Sections {
        Sections::FIXED
    }

    fn new((): &(), _: RawRecord<'_>) -> Result<Self, VcfReadError> {
        Ok(Self)
    }
}

fn seqair_run<F: Fields>(
    file: &VcfFile,
    threads: usize,
    sample: Option<&str>,
    mut each: impl FnMut(&F) -> Result<u64, VcfReadError>,
) -> (u64, u64) {
    let mut b = file.reader::<F>();
    if let Some(sample) = sample {
        b = b.sample(sample);
    }
    if let Some(n) = NonZeroUsize::new(threads) {
        b = b.threads(n);
    }
    let mut r = b.open().unwrap();
    let (mut n, mut x) = (0u64, 0u64);
    while let Some(rec) = r.next_record().unwrap() {
        n += 1;
        x = x.wrapping_add(each(rec).unwrap());
    }
    (n, x)
}

fn seqair(pattern: &str, path: &str, threads: usize) -> (u64, u64) {
    let file = VcfFile::open(path).unwrap();
    let h = file.header();
    let i1 = ikey(&file, "INFO1");
    let i2 = ikey(&file, "INFO2");
    let gt: Option<FormatGt> = h.format_key("GT").ok();
    let fmt: Option<FormatInt> = std::env::var("FMT").ok().map(|n| h.format_key(&n).unwrap());
    let sample_name = std::env::var("SAMPLE").ok();
    let one_sample = |s: seqair::vcf::SampleRef<'_>| -> Result<u64, VcfReadError> {
        let mut x = 0;
        if let Some(gt) = &gt {
            x += gt_sum(s.get(gt)?);
        }
        if let Some(f) = &fmt {
            x += s.get(f)?.map_or(0, |v| v as u64);
        }
        Ok(x)
    };
    let t = threads;
    match pattern {
        "count" => seqair_run::<Count>(&file, t, None, |_| Ok(1)),
        "locus" => seqair_run::<Locus>(&file, t, None, |l| Ok(start(l.start()) + start(l.end()))),
        "site" => seqair_run::<(Locus, Site)>(&file, t, None, |(l, s)| {
            Ok(start(l.start()) + u64::from(pass(s)))
        }),
        "info0" => seqair_run::<(Locus, Info)>(&file, t, None, |(l, _)| Ok(start(l.start()))),
        "info1" => seqair_run::<(Locus, Info)>(&file, t, None, |(l, i)| {
            Ok(start(l.start()) + touch_info(i, i1.as_ref().unwrap())?)
        }),
        "site_info1" => seqair_run::<(Locus, Site, Info)>(&file, t, None, |(l, s, i)| {
            Ok(start(l.start()) + u64::from(pass(s)) + touch_info(i, i1.as_ref().unwrap())?)
        }),
        "filter_info2" => seqair_run::<(Locus, Site, Info)>(&file, t, None, |(l, s, i)| {
            let mut x = start(l.start());
            if pass(s) {
                x += touch_info(i, i1.as_ref().unwrap())? + touch_info(i, i2.as_ref().unwrap())?;
            }
            Ok(x)
        }),
        "only_sample" => seqair_run::<(Locus, OnlySample)>(&file, t, None, |(l, s)| {
            let mut x = start(l.start());
            if let Some(gt) = &gt {
                x += gt_sum(s.get(gt)?);
            }
            if let Some(f) = &fmt {
                x += s.get(f)?.map_or(0, |v| v as u64);
            }
            Ok(x)
        }),
        "sample" => seqair_run::<(Locus, Sample)>(&file, t, sample_name.as_deref(), |(l, s)| {
            let mut x = start(l.start());
            if let Some(gt) = &gt {
                x += gt_sum(s.get(gt)?);
            }
            if let Some(f) = &fmt {
                x += s.get(f)?.map_or(0, |v| v as u64);
            }
            Ok(x)
        }),
        "column" => seqair_run::<(Locus, Samples)>(&file, t, None, |(_, s)| {
            Ok(s.column(gt.as_ref().unwrap())?.map_or(0, |c| c.map(gt_sum).sum()))
        }),
        "full" => seqair_run::<FullRecord>(&file, t, None, |r| Ok(start(r.locus().start()))),
        "all" => {
            let idx: Vec<_> = h.sample_indices().collect();
            seqair_run::<FullRecord>(&file, t, None, |r| {
                let mut x = start(r.locus().start()) + u64::from(pass(r.site()));
                for v in r.info().iter() {
                    x += dyn_sum(v?.2);
                }
                for &s in &idx {
                    x += one_sample(r.samples().sample(s))?;
                }
                Ok(x)
            })
        }
        p => panic!("pattern {p}"),
    }
}

fn info_spec(var: &str) -> Option<(Vec<u8>, String)> {
    let s = std::env::var(var).ok()?;
    let (n, t) = s.split_once(':').unwrap();
    Some((n.as_bytes().to_vec(), t.to_owned()))
}

fn hts_info(
    rec: &hts::Record,
    (name, typ): &(Vec<u8>, String),
    buf: &mut hts::record::Buffer,
) -> u64 {
    let mut info = rec.info_shared_buffer(name, &mut *buf);
    match typ.as_str() {
        "flag" => u64::from(info.flag().unwrap()),
        "int" | "ints" => info.integer().unwrap().map_or(0, |v| v.iter().map(|&x| x as u64).sum()),
        _ => info.float().unwrap().map_or(0, |v| v.iter().map(|x| u64::from(x.to_bits())).sum()),
    }
}

fn hts_gt_sum(g: &hts::record::Genotype) -> u64 {
    g.iter().map(|a| a.index().map_or(7, u64::from)).sum()
}

fn htslib(pattern: &str, path: &str, threads: usize) -> (u64, u64) {
    let mut r = hts::Reader::from_path(path).unwrap();
    if threads > 0 {
        r.set_threads(threads).unwrap();
    }
    let header = r.header().clone();
    let pass = header.name_to_id(b"PASS").unwrap();
    let i1 = info_spec("INFO1");
    let i2 = info_spec("INFO2");
    let fmt = std::env::var("FMT").ok();
    let sample = std::env::var("SAMPLE").ok().map(|s| header.sample_id(s.as_bytes()).unwrap());
    let mut rec = r.empty_record();
    let mut buf = hts::record::Buffer::new();
    let (mut n, mut x) = (0u64, 0u64);
    while let Some(res) = r.read(&mut rec) {
        res.unwrap();
        n += 1;
        x = x.wrapping_add(match pattern {
            "count" => 1,
            "locus" => rec.pos() as u64 + 1 + rec.end() as u64,
            "info1" => rec.pos() as u64 + 1 + hts_info(&rec, i1.as_ref().unwrap(), &mut buf),
            "filter_info2" => {
                let mut x = rec.pos() as u64 + 1;
                if rec.has_filter(&pass) {
                    x += hts_info(&rec, i1.as_ref().unwrap(), &mut buf);
                    x += hts_info(&rec, i2.as_ref().unwrap(), &mut buf);
                }
                x
            }
            "sample" => {
                let s = sample.unwrap();
                let mut x = rec.pos() as u64 + 1;
                x += hts_gt_sum(&rec.genotypes_shared_buffer(&mut buf).unwrap().get(s));
                if let Some(f) = &fmt {
                    // A record without the tag is an error in rust-htslib.
                    if let Ok(v) = rec.format_shared_buffer(f.as_bytes(), &mut buf).integer() {
                        x += v.get(s).and_then(|v| v.first()).map_or(0, |&v| v as u64);
                    }
                }
                x
            }
            "column" => {
                let g = rec.genotypes_shared_buffer(&mut buf).unwrap();
                (0..rec.sample_count() as usize).map(|s| hts_gt_sum(&g.get(s))).sum()
            }
            "all" => {
                rec.unpack();
                rec.pos() as u64
            }
            p => panic!("pattern {p}"),
        });
    }
    (n, x)
}

fn noodles(pattern: &str, path: &str, threads: usize) -> (u64, u64) {
    use noodles::vcf::variant::record::samples::Sample as _;
    use noodles::vcf::variant::record::{Filters as _, Record as _, Samples as _};
    let file = std::fs::File::open(path).unwrap();
    let bgzf: Box<dyn std::io::BufRead> = if let Some(n) = NonZeroUsize::new(threads) {
        Box::new(noodles::bgzf::io::MultithreadedReader::with_worker_count(n, file))
    } else {
        Box::new(noodles::bgzf::io::Reader::new(file))
    };
    let i1 = info_spec("INFO1").map(|s| String::from_utf8(s.0).unwrap());
    let i2 = info_spec("INFO2").map(|s| String::from_utf8(s.0).unwrap());
    let fmt = std::env::var("FMT").ok();
    let (mut n, mut x) = (0u64, 0u64);
    let touch = |v: Option<
        std::io::Result<Option<noodles::vcf::variant::record::info::field::Value<'_>>>,
    >| { u64::from(black_box(v).is_some_and(|v| v.unwrap().is_some())) };
    macro_rules! body {
        ($rec:expr, $header:expr, $end:expr) => {{
            let rec = &$rec;
            let header = &$header;
            let sample = std::env::var("SAMPLE")
                .ok()
                .map(|s| header.sample_names().get_index_of(&s).unwrap());
            match pattern {
                "count" => 1,
                "locus" => {
                    let s = rec.variant_start().map_or(0, |p| p.unwrap().get() as u64);
                    s + $end(rec, header)
                }
                "filter_info2" => {
                    let mut x = rec.variant_start().map_or(0, |p| p.unwrap().get() as u64);
                    if rec.filters().is_pass(header).unwrap() {
                        let info = rec.info();
                        x += touch(info.get(header, i1.as_deref().unwrap()));
                        x += touch(info.get(header, i2.as_deref().unwrap()));
                    }
                    x
                }
                "sample" => {
                    let samples = noodles::vcf::variant::Record::samples(rec).unwrap();
                    let s = samples.iter().nth(sample.unwrap()).unwrap();
                    let mut x =
                        u64::from(s.get(header, "GT").is_some_and(|v| v.unwrap().is_some()));
                    if let Some(f) = &fmt {
                        x += u64::from(s.get(header, f).is_some_and(|v| v.unwrap().is_some()));
                    }
                    x
                }
                p => panic!("pattern {p}"),
            }
        }};
    }
    if path.ends_with(".bcf") {
        let mut r = noodles::bcf::io::Reader::from(bgzf);
        let header = r.read_header().unwrap();
        if pattern == "all" {
            let mut rec = noodles::vcf::variant::RecordBuf::default();
            while r.read_record_buf(&header, &mut rec).unwrap() > 0 {
                n += 1;
                x = x.wrapping_add(rec.variant_start().map_or(0, |p| p.get() as u64));
            }
            return (n, x);
        }
        let mut rec = noodles::bcf::Record::default();
        while r.read_record(&mut rec).unwrap() > 0 {
            n += 1;
            let end =
                |r: &noodles::bcf::Record, _: &noodles::vcf::Header| r.end().unwrap().get() as u64;
            x = x.wrapping_add(body!(rec, header, end));
        }
    } else {
        let mut r = noodles::vcf::io::Reader::new(bgzf);
        let header = r.read_header().unwrap();
        if pattern == "all" {
            let mut rec = noodles::vcf::variant::RecordBuf::default();
            while r.read_record_buf(&header, &mut rec).unwrap() > 0 {
                n += 1;
                x = x.wrapping_add(rec.variant_start().map_or(0, |p| p.get() as u64));
            }
            return (n, x);
        }
        let mut rec = noodles::vcf::Record::default();
        while r.read_record(&mut rec).unwrap() > 0 {
            n += 1;
            let end = |r: &noodles::vcf::Record, h: &noodles::vcf::Header| {
                r.variant_end(h).unwrap().get() as u64
            };
            x = x.wrapping_add(body!(rec, header, end));
        }
    }
    (n, x)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let (lib, pattern, path) = (&args[1], &args[2], &args[3]);
    let threads: usize = args.get(4).map_or(0, |s| s.parse().unwrap());
    let reps: usize = args.get(5).map_or(3, |s| s.parse().unwrap());
    let mut best = f64::MAX;
    let mut out = (0, 0);
    for _ in 0..reps {
        let t = Instant::now();
        out = match lib.as_str() {
            "seqair" => seqair(pattern, path, threads),
            "htslib" => htslib(pattern, path, threads),
            "noodles" => noodles(pattern, path, threads),
            l => panic!("lib {l}"),
        };
        best = best.min(t.elapsed().as_secs_f64());
    }
    let (n, x) = black_box(out);
    println!(
        "{lib}\t{pattern}\t{}\t{threads}\t{n}\t{:.1}\t{:.2}\t{x}",
        path.rsplit('/').next().unwrap(),
        best * 1e9 / n as f64,
        n as f64 / best / 1e6
    );
}
