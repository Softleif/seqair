//! seqair's writer as the reader's generator (`r[vcf_values.equivalence]`):
//! random headers and records written as VCF text, BGZF VCF and BCF read
//! back, through `FullRecord`, as what was written, and the same through
//! every accessor in all three.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    reason = "test code"
)]

use hegel::prelude::*;
use seqair::vcf::record_encoder::{FilterFieldDef, FormatFieldDef, InfoFieldDef};
use seqair::vcf::{
    Allele, Alleles, Arr, ContigDef, ContigId, DynValue, FilterId, FilterStatus, Flag,
    FormatFloatOpts, FormatFloats, FormatGt, FormatInt, FormatString, FullRecord, Genotype, Gt,
    InfoFlag, InfoFloat, InfoFloats, InfoInt, InfoIntOpts, InfoInts, InfoString, Number, OptArr,
    OutputFormat, Scalar, Str, VcfFile, VcfHeader, VcfPos, Writer,
};
use seqair_types::{Pos1, SmallVec, SmolStr};

// ── What is written ─────────────────────────────────────────────────────

/// The header's optional fields, each `Some` when declared.
#[derive(Debug, Default)]
struct Keys {
    i1: Option<InfoInt>,
    f1: Option<InfoFloat>,
    ia: Option<InfoInts>,
    fa: Option<InfoFloats>,
    io: Option<InfoIntOpts>,
    fl: Option<InfoFlag>,
    s1: Option<InfoString>,
    gt: Option<FormatGt>,
    dp: Option<FormatInt>,
    af: Option<FormatFloats>,
    fo: Option<FormatFloatOpts>,
    st: Option<FormatString>,
}

struct Setup {
    header: VcfHeader,
    contigs: Vec<ContigId>,
    filters: Vec<FilterId>,
    keys: Keys,
    n_samples: usize,
}

#[derive(Debug, Clone)]
enum Filter {
    Pass,
    Missing,
    Failed(Vec<usize>),
}

#[derive(Debug, Clone, Default)]
struct Info {
    i1: Option<i32>,
    f1: Option<f32>,
    ia: Option<Vec<i32>>,
    fa: Option<Vec<f32>>,
    io: Option<Vec<Option<i32>>>,
    fl: bool,
    s1: Option<String>,
}

#[derive(Debug, Clone, Default)]
struct Format {
    gt: Option<Vec<Genotype>>,
    dp: Option<Vec<i32>>,
    af: Option<Vec<Vec<f32>>>,
    fo: Option<Vec<Vec<Option<f32>>>>,
    st: Option<Vec<String>>,
}

#[derive(Debug, Clone)]
struct Record {
    contig: usize,
    pos: u32,
    reference: String,
    alts: Vec<String>,
    qual: Option<f32>,
    filter: Filter,
    info: Info,
    format: Format,
}

fn setup(tc: &TestCase) -> Setup {
    let mut b = VcfHeader::builder();
    let n_contigs = tc.draw(gs::integers::<usize>().min_value(1).max_value(3));
    let contigs: Vec<ContigId> = ["chr1", "chr2", "chrX"][..n_contigs]
        .iter()
        .map(|c| b.register_contig(*c, ContigDef::new(Some(10_000_000))).unwrap())
        .collect();
    let mut b = b.filters();
    let n_filters = tc.draw(gs::integers::<usize>().max_value(3));
    let filters = [("q10", "Q10"), ("lowDP", "low depth"), ("x", "x")][..n_filters]
        .iter()
        .map(|(n, d)| b.register_filter(&FilterFieldDef::new(n, d)).unwrap())
        .collect();
    let mut b = b.infos();
    let mut keys = Keys::default();
    let on = || tc.draw(gs::booleans());
    let unknown = Number::Unknown;
    if on() {
        keys.i1 = Some(b.register_info(&InfoFieldDef::<Scalar<i32>>::new("I1", "i")).unwrap());
    }
    if on() {
        keys.f1 = Some(b.register_info(&InfoFieldDef::<Scalar<f32>>::new("F1", "f")).unwrap());
    }
    if on() {
        keys.ia =
            Some(b.register_info(&InfoFieldDef::<Arr<i32>>::new("IA", unknown, "ia")).unwrap());
    }
    if on() {
        keys.fa =
            Some(b.register_info(&InfoFieldDef::<Arr<f32>>::new("FA", unknown, "fa")).unwrap());
    }
    if on() {
        keys.io =
            Some(b.register_info(&InfoFieldDef::<OptArr<i32>>::new("IO", unknown, "io")).unwrap());
    }
    if on() {
        keys.fl = Some(b.register_info(&InfoFieldDef::<Flag>::new("FL", "fl")).unwrap());
    }
    if on() {
        keys.s1 =
            Some(b.register_info(&InfoFieldDef::<Str>::new("S1", Number::Count(1), "s")).unwrap());
    }
    let mut b = b.formats();
    if on() {
        keys.gt = Some(b.register_format(&FormatFieldDef::<Gt>::new("GT", "gt")).unwrap());
    }
    if on() {
        keys.dp = Some(b.register_format(&FormatFieldDef::<Scalar<i32>>::new("DP", "dp")).unwrap());
    }
    if on() {
        keys.af =
            Some(b.register_format(&FormatFieldDef::<Arr<f32>>::new("AF", unknown, "af")).unwrap());
    }
    if on() {
        keys.fo = Some(
            b.register_format(&FormatFieldDef::<OptArr<f32>>::new("FO", unknown, "fo")).unwrap(),
        );
    }
    if on() {
        keys.st = Some(
            b.register_format(&FormatFieldDef::<Str>::new("ST", Number::Count(1), "st")).unwrap(),
        );
    }
    let mut b = b.samples();
    let n_samples = tc.draw(gs::integers::<usize>().max_value(4));
    for i in 0..n_samples {
        b.add_sample(format!("S{i}")).unwrap();
    }
    Setup { header: b.build().unwrap(), contigs, filters, keys, n_samples }
}

/// Integers BCF can store (the int32 reserved band excluded).
fn int() -> impl PrintableGenerator<i32> {
    hegel::one_of!(
        gs::integers::<i32>().min_value(-200).max_value(200),
        gs::integers::<i32>().min_value(i32::MIN + 8)
    )
}

fn float() -> impl PrintableGenerator<f32> {
    gs::floats::<f32>().min_value(f32::MIN).max_value(f32::MAX)
}

fn word(tc: &TestCase) -> String {
    let n = tc.draw(gs::integers::<usize>().min_value(1).max_value(6));
    (0..n).map(|_| tc.draw(gs::sampled_from(vec!['a', 'Z', '0', '_', 'x']))).collect()
}

fn bases(tc: &TestCase, min: usize) -> String {
    let n = tc.draw(gs::integers::<usize>().min_value(min).max_value(4));
    (0..n).map(|_| tc.draw(gs::sampled_from(vec!['A', 'C', 'G', 'T']))).collect()
}

fn record(tc: &TestCase, s: &Setup) -> Record {
    let reference = bases(tc, 1);
    let n_alts = tc.draw(gs::integers::<usize>().max_value(3));
    let alts: Vec<String> = (0..n_alts)
        .map(|_| match tc.draw(gs::integers::<u8>().max_value(5)) {
            0 => "<DEL>".to_owned(),
            1 => "*".to_owned(),
            _ => bases(tc, 1),
        })
        .collect();
    let n_alleles = 1 + alts.len() as u16;
    let filter = match tc.draw(gs::integers::<u8>().max_value(2)) {
        0 => Filter::Pass,
        1 => Filter::Missing,
        _ if s.filters.is_empty() => Filter::Pass,
        _ => {
            let n = tc.draw(gs::integers::<usize>().min_value(1).max_value(s.filters.len()));
            Filter::Failed((0..n).collect())
        }
    };
    let maybe = || tc.draw(gs::booleans());
    let k = &s.keys;
    let info = Info {
        i1: if k.i1.is_some() && maybe() { Some(tc.draw(int())) } else { None },
        f1: if k.f1.is_some() && maybe() { Some(tc.draw(float())) } else { None },
        ia: if k.ia.is_some() && maybe() {
            Some(tc.draw(gs::vecs(int()).min_size(1).max_size(4)))
        } else {
            None
        },
        fa: if k.fa.is_some() && maybe() {
            Some(tc.draw(gs::vecs(float()).min_size(1).max_size(4)))
        } else {
            None
        },
        io: if k.io.is_some() && maybe() {
            Some(tc.draw(gs::vecs(gs::optional(int())).min_size(1).max_size(4)))
        } else {
            None
        },
        fl: k.fl.is_some() && maybe(),
        s1: if k.s1.is_some() && maybe() { Some(word(tc)) } else { None },
    };
    let n = s.n_samples;
    let format = if n == 0 {
        Format::default()
    } else {
        let ploidy = tc.draw(gs::integers::<usize>().min_value(1).max_value(3));
        let genotype = |tc: &TestCase| {
            let alleles: SmallVec<Option<u16>, 2> = (0..ploidy)
                .map(|_| tc.draw(gs::optional(gs::integers::<u16>().max_value(n_alleles - 1))))
                .collect();
            let phased: SmallVec<bool, 2> = (1..ploidy).map(|_| tc.draw(gs::booleans())).collect();
            Genotype { alleles, phased }
        };
        Format {
            gt: if k.gt.is_some() && maybe() {
                Some((0..n).map(|_| genotype(tc)).collect())
            } else {
                None
            },
            dp: if k.dp.is_some() && maybe() {
                Some((0..n).map(|_| tc.draw(int())).collect())
            } else {
                None
            },
            af: if k.af.is_some() && maybe() {
                Some((0..n).map(|_| tc.draw(gs::vecs(float()).max_size(3))).collect())
            } else {
                None
            },
            fo: if k.fo.is_some() && maybe() {
                Some((0..n).map(|_| tc.draw(gs::vecs(gs::optional(float())).max_size(3))).collect())
            } else {
                None
            },
            st: if k.st.is_some() && maybe() {
                Some(
                    (0..n)
                        .map(|_| if tc.draw(gs::booleans()) { word(tc) } else { String::new() })
                        .collect(),
                )
            } else {
                None
            },
        }
    };
    Record {
        contig: tc.draw(gs::integers::<usize>().max_value(s.contigs.len() - 1)),
        pos: tc.draw(gs::integers::<u32>().min_value(1).max_value(5_000_000)),
        reference,
        alts,
        qual: tc.draw(gs::optional(float())),
        filter,
        info,
        format,
    }
}

fn write(s: &Setup, records: &[Record], format: OutputFormat) -> Vec<u8> {
    let mut out = Vec::new();
    let mut w = Writer::new(&mut out, format).write_header(&s.header).unwrap();
    let k = &s.keys;
    for r in records {
        let alleles = Alleles::complex(
            SmolStr::from(&r.reference),
            r.alts.iter().map(SmolStr::from).collect(),
        );
        let enc = w
            .begin_record(&s.contigs[r.contig], Pos1::new(r.pos).unwrap(), &alleles, r.qual)
            .unwrap();
        let mut enc = match &r.filter {
            Filter::Pass => enc.filter_pass(),
            Filter::Missing => enc.no_filter(),
            Filter::Failed(ids) => enc.filter_fail(ids.iter().map(|&i| &s.filters[i])).unwrap(),
        };
        let i = &r.info;
        if let (Some(key), Some(v)) = (&k.i1, i.i1) {
            key.encode(&mut enc, v).unwrap();
        }
        if let (Some(key), Some(v)) = (&k.f1, i.f1) {
            key.encode(&mut enc, v).unwrap();
        }
        if let (Some(key), Some(v)) = (&k.ia, &i.ia) {
            key.encode(&mut enc, v).unwrap();
        }
        if let (Some(key), Some(v)) = (&k.fa, &i.fa) {
            key.encode(&mut enc, v).unwrap();
        }
        if let (Some(key), Some(v)) = (&k.io, &i.io) {
            key.encode(&mut enc, v).unwrap();
        }
        if let (Some(key), true) = (&k.fl, i.fl) {
            key.encode(&mut enc).unwrap();
        }
        if let (Some(key), Some(v)) = (&k.s1, &i.s1) {
            key.encode(&mut enc, v).unwrap();
        }
        if s.n_samples == 0 {
            enc.emit().unwrap();
            continue;
        }
        let mut enc = enc.begin_samples();
        let f = &r.format;
        if let (Some(key), Some(v)) = (&k.gt, &f.gt) {
            key.encode(&mut enc, v).unwrap();
        }
        if let (Some(key), Some(v)) = (&k.dp, &f.dp) {
            key.encode(&mut enc, v).unwrap();
        }
        if let (Some(key), Some(v)) = (&k.af, &f.af) {
            let per_sample: Vec<&[f32]> = v.iter().map(Vec::as_slice).collect();
            key.encode(&mut enc, &per_sample).unwrap();
        }
        if let (Some(key), Some(v)) = (&k.fo, &f.fo) {
            let per_sample: Vec<&[Option<f32>]> = v.iter().map(Vec::as_slice).collect();
            seqair::vcf::FormatEncoder::format_float_opts(&mut enc, key.id(), &per_sample).unwrap();
        }
        if let (Some(key), Some(v)) = (&k.st, &f.st) {
            let per_sample: Vec<&str> = v.iter().map(String::as_str).collect();
            key.encode(&mut enc, &per_sample).unwrap();
        }
        enc.emit().unwrap();
    }
    w.finish().unwrap();
    out
}

// ── What is read ────────────────────────────────────────────────────────

fn bits(v: impl IntoIterator<Item = f32>) -> Vec<u32> {
    v.into_iter().map(f32::to_bits).collect()
}

fn opt_bits(v: impl IntoIterator<Item = Option<f32>>) -> Vec<Option<u32>> {
    v.into_iter().map(|f| f.map(f32::to_bits)).collect()
}

/// A vector written with no elements, or one missing one, reads as absent.
fn missing<T>(v: &[Option<T>]) -> bool {
    v.is_empty() || (v.len() == 1 && v[0].is_none())
}

/// What a written float reads back as: itself from BCF; from VCF text,
/// what C's `%g` (6 significant digits, as the writer prints) parses to.
fn as_read(v: f32, text: bool) -> f32 {
    if text { format!("{:.5e}", f64::from(v)).parse().unwrap() } else { v }
}

/// What a written record reads back as from `format`.
fn as_read_record(w: &Record, format: OutputFormat) -> Record {
    let text = format != OutputFormat::Bcf;
    let f = |v: f32| as_read(v, text);
    let mut r = w.clone();
    r.qual = r.qual.map(f);
    r.info.f1 = r.info.f1.map(f);
    if let Some(v) = &mut r.info.fa {
        v.iter_mut().for_each(|x| *x = f(*x));
    }
    if let Some(v) = &mut r.format.af {
        v.iter_mut().flatten().for_each(|x| *x = f(*x));
    }
    if let Some(v) = &mut r.format.fo {
        v.iter_mut().flatten().flatten().for_each(|x| *x = f(*x));
    }
    r
}

fn check(s: &Setup, written: &Record, read: &FullRecord, format: OutputFormat) {
    let at = format!("{format:?}: {written:?}");
    let written = &as_read_record(written, format);
    let locus = read.locus();
    assert_eq!(locus.contig(), &s.contigs[written.contig], "{at}");
    let start = written.pos;
    let end = start + written.reference.len() as u32 - 1;
    assert_eq!(locus.start(), VcfPos::At(Pos1::new(start).unwrap()), "{at}");
    assert_eq!(locus.end(), VcfPos::At(Pos1::new(end).unwrap()), "{at}");

    let site = read.site();
    assert_eq!(site.ids().count(), 0, "{at}");
    assert_eq!(site.alleles().reference(), written.reference.as_bytes(), "{at}");
    let alts: Vec<Allele<'_>> = site.alleles().alts().collect();
    let expected: Vec<Allele<'_>> = written
        .alts
        .iter()
        .map(|a| match a.as_str() {
            "*" => Allele::Overlapping,
            "<DEL>" => Allele::Symbolic("<DEL>"),
            b => Allele::Bases(b.as_bytes()),
        })
        .collect();
    assert_eq!(alts, expected, "{at}");
    assert_eq!(site.qual().map(f32::to_bits), written.qual.map(f32::to_bits), "{at}");
    match (&written.filter, site.filters().status()) {
        (Filter::Pass, FilterStatus::Pass) | (Filter::Missing, FilterStatus::Missing) => {}
        (Filter::Failed(ids), FilterStatus::Failed(read)) => {
            let read: Vec<&FilterId> = read.collect();
            let expected: Vec<&FilterId> = ids.iter().map(|&i| &s.filters[i]).collect();
            assert_eq!(read, expected, "{at}");
            for f in &expected {
                assert!(site.filters().contains(f).unwrap(), "{at}");
            }
        }
        (w, r) => panic!("{at}: FILTER {w:?} read as {r:?}"),
    }

    let (k, w, info) = (&s.keys, &written.info, read.info());
    if let Some(key) = &k.i1 {
        assert_eq!(info.get(key).unwrap(), w.i1, "{at}");
    }
    if let Some(key) = &k.f1 {
        assert_eq!(info.get(key).unwrap().map(f32::to_bits), w.f1.map(f32::to_bits), "{at}");
    }
    if let Some(key) = &k.ia {
        assert_eq!(info.get(key).unwrap().map(Iterator::collect::<Vec<_>>), w.ia, "{at}");
    }
    if let Some(key) = &k.fa {
        assert_eq!(info.get(key).unwrap().map(bits), w.fa.clone().map(bits), "{at}");
    }
    if let Some(key) = &k.io {
        let expected = w.io.clone().filter(|v| !missing(v));
        assert_eq!(info.get(key).unwrap().map(Iterator::collect::<Vec<_>>), expected, "{at}");
    }
    if let Some(key) = &k.fl {
        assert_eq!(info.get(key).unwrap(), w.fl, "{at}");
    }
    if let Some(key) = &k.s1 {
        assert_eq!(info.get(key).unwrap(), w.s1.as_deref(), "{at}");
    }
    let present = [
        w.i1.is_some(),
        w.f1.is_some(),
        w.ia.is_some(),
        w.fa.is_some(),
        w.io.is_some(),
        w.fl,
        w.s1.is_some(),
    ];
    assert_eq!(info.iter().unwrap().count(), present.iter().filter(|&&p| p).count(), "{at}");

    let (f, samples) = (&written.format, read.samples());
    assert_eq!(samples.len(), s.n_samples, "{at}");
    for (i, idx) in s.header.sample_indices().enumerate() {
        let sample = samples.sample(idx);
        if let Some(key) = &k.gt {
            // A missing allele keeps its phase bit; a missing haploid call
            // reads as absent (`r[vcf_record.gt_first_phase+3]`).
            let expected =
                f.gt.as_ref()
                    .map(|g| &g[i])
                    .filter(|g| !(g.alleles.len() == 1 && g.alleles[0].is_none()))
                    .cloned();
            let expected = expected.as_ref();
            let got = sample.get(key).unwrap();
            assert_eq!(
                got.as_ref().map(|g| g.to_genotype()).as_ref(),
                expected,
                "{at}: sample {i} GT"
            );
            let column = samples
                .column(key)
                .unwrap()
                .map(|c| c.map(|g| g.map(|g| g.to_genotype())).collect::<Vec<_>>());
            assert_eq!(column.and_then(|c| c[i].clone()).as_ref(), expected, "{at}: GT column");
        }
        if let Some(key) = &k.dp {
            assert_eq!(
                sample.get(key).unwrap(),
                f.dp.as_ref().map(|v| v[i]),
                "{at}: sample {i} DP"
            );
        }
        if let Some(key) = &k.af {
            let expected = f.af.as_ref().map(|v| bits(v[i].clone())).filter(|v| !v.is_empty());
            assert_eq!(sample.get(key).unwrap().map(bits), expected, "{at}: sample {i} AF");
        }
        if let Some(key) = &k.fo {
            let expected =
                f.fo.as_ref().map(|v| v[i].clone()).filter(|v| !missing(v)).map(opt_bits);
            assert_eq!(sample.get(key).unwrap().map(opt_bits), expected, "{at}: sample {i} FO");
        }
        if let Some(key) = &k.st {
            let expected = f.st.as_ref().map(|v| v[i].as_str()).filter(|v| !v.is_empty());
            assert_eq!(sample.get(key).unwrap(), expected, "{at}: sample {i} ST");
        }
    }
    assert_eq!(read.header_id(), s.header.id());
}

/// Every value every accessor reads from a record, as text: floats by the
/// bits of `float(value)`, a genotype's phases only before an allele that
/// is present (see the shortcut in `check`) unless `exact`.
fn digest(
    header: &VcfHeader,
    rec: &FullRecord,
    exact: bool,
    float: impl Fn(f32) -> f32,
) -> Vec<String> {
    let bits = |v: Option<f32>| v.map(|v| float(v).to_bits());
    let mut out = Vec::new();
    let l = rec.locus();
    out.push(format!("{} {:?} {:?}", l.contig().name(), l.start(), l.end()));
    let site = rec.site();
    out.push(format!("{:?} {:?}", site.ids().collect::<Vec<_>>(), site.alleles().reference()));
    out.push(format!("{:?}", site.alleles().alts().collect::<Vec<_>>()));
    let filter = match site.filters().status() {
        FilterStatus::Pass => "PASS".to_owned(),
        FilterStatus::Missing => ".".to_owned(),
        FilterStatus::Failed(f) => f.map(FilterId::name).collect::<Vec<_>>().join(";"),
    };
    out.push(format!("{:?} {filter}", bits(site.qual())));
    let dyn_value = |v: &DynValue<'_>| match v {
        DynValue::Int(v) => format!("{:?}", v.clone().collect::<Vec<_>>()),
        DynValue::Float(v) => format!("{:?}", v.clone().map(bits).collect::<Vec<_>>()),
        DynValue::Flag => "flag".to_owned(),
        DynValue::Str(s) => format!("{s:?}"),
        DynValue::Gt(g) => format!(
            "{:?}",
            g.as_ref().map(|g| g
                .alleles()
                .enumerate()
                .map(|(i, a)| (
                    a.allele,
                    (exact || a.allele.is_some()) && (i > 0 || exact) && a.phased
                ))
                .collect::<Vec<_>>())
        ),
    };
    for field in rec.info().iter().unwrap() {
        let (name, _, v) = field.unwrap();
        out.push(format!("INFO {name} {}", dyn_value(&v)));
    }
    for idx in header.sample_indices() {
        for field in rec.samples().sample(idx).iter().unwrap() {
            let (name, _, v) = field.unwrap();
            out.push(format!("{} {name} {}", idx.index(), dyn_value(&v)));
        }
    }
    out
}

/// The writer's records read back from `format` with the keys of the
/// header the reader parsed.
fn read_back(s: &Setup, records: &[Record], format: OutputFormat) -> (Setup, Vec<FullRecord>) {
    let file = VcfFile::from_bytes(write(s, records, format), None).unwrap();
    let expected = match format {
        OutputFormat::Bcf => seqair::vcf::reader::VcfFormat::Bcf,
        OutputFormat::VcfGz => seqair::vcf::reader::VcfFormat::VcfGz,
        _ => seqair::vcf::reader::VcfFormat::Vcf,
    };
    assert_eq!(file.format(), expected);
    let read: Vec<FullRecord> =
        file.reader::<FullRecord>().open().unwrap().into_iter().collect::<Result<_, _>>().unwrap();
    // The reader's header is the writer's, re-parsed: keys resolve against it.
    let s = Setup {
        contigs: reresolve_contigs(&file, s),
        filters: reresolve_filters(&file, s),
        keys: rekey(&file, &s.keys),
        header: file.header().clone(),
        n_samples: s.n_samples,
    };
    (s, read)
}

// r[verify vcf_values.equivalence]
// r[verify vcf_values.access]
// r[verify vcf_values.dynamic.missing]
// r[verify vcf_values.bcf_sentinels]
// r[verify vcf_values.text_numbers]
// r[verify vcf_values.text_strings]
// r[verify vcf_values.gt]
// r[verify vcf_values.qual]
// r[verify vcf_fields.site]
// r[verify vcf_fields.info]
// r[verify vcf_fields.samples]
// r[verify vcf_reader.formats]
// r[verify vcf_reader.oracles]
/// Random headers and records, written with seqair's writer as VCF text,
/// BGZF VCF and BCF, read back through `FullRecord` as written (floats as
/// the writer's 6-digit `%g` reads back, in text): every value kind,
/// missing elements, flags, strings, genotypes of every ploidy and phase,
/// several samples, every FILTER state and missing QUAL. Every accessor
/// then reads the same from all three.
#[hegel::test]
fn every_format_written_by_the_writer_reads_back_as_written(tc: TestCase) {
    let s = setup(&tc);
    let n = tc.draw(gs::integers::<usize>().max_value(12));
    let mut records: Vec<Record> = (0..n).map(|_| record(&tc, &s)).collect();
    records.sort_by_key(|r| (r.contig, r.pos));
    let mut digests = Vec::new();
    for format in [OutputFormat::Bcf, OutputFormat::Vcf, OutputFormat::VcfGz] {
        let (read_setup, read) = read_back(&s, &records, format);
        assert_eq!(read.len(), records.len(), "{format:?}");
        for (w, r) in records.iter().zip(&read) {
            check(&read_setup, w, r, format);
        }
        // BCF holds the floats VCF text rounds: compare BCF's rounded.
        let rounded = |v: f32| as_read(v, format == OutputFormat::Bcf);
        let h = &read_setup.header;
        let exact: Vec<_> = read.iter().map(|r| digest(h, r, true, |v| v)).collect();
        let loose: Vec<_> = read.iter().map(|r| digest(h, r, false, rounded)).collect();
        digests.push((exact, loose));
    }
    let [(_, bcf), (vcf, vcf_loose), (gz, _)] = &digests[..] else { unreachable!() };
    assert_eq!(vcf, gz, "VCF text and BGZF VCF read differently");
    assert_eq!(vcf_loose, bcf, "VCF text and BCF read differently");
}

fn reresolve_contigs(file: &VcfFile, s: &Setup) -> Vec<ContigId> {
    s.contigs.iter().map(|c| file.header().contig(c.name()).unwrap()).collect()
}

fn reresolve_filters(file: &VcfFile, s: &Setup) -> Vec<FilterId> {
    s.filters.iter().map(|f| file.header().filter(f.name()).unwrap()).collect()
}

fn rekey(file: &VcfFile, k: &Keys) -> Keys {
    let h = file.header();
    Keys {
        i1: k.i1.as_ref().map(|_| h.info_key("I1").unwrap()),
        f1: k.f1.as_ref().map(|_| h.info_key("F1").unwrap()),
        ia: k.ia.as_ref().map(|_| h.info_key("IA").unwrap()),
        fa: k.fa.as_ref().map(|_| h.info_key("FA").unwrap()),
        io: k.io.as_ref().map(|_| h.info_key("IO").unwrap()),
        fl: k.fl.as_ref().map(|_| h.info_key("FL").unwrap()),
        s1: k.s1.as_ref().map(|_| h.info_key("S1").unwrap()),
        gt: k.gt.as_ref().map(|_| h.format_key("GT").unwrap()),
        dp: k.dp.as_ref().map(|_| h.format_key("DP").unwrap()),
        af: k.af.as_ref().map(|_| h.format_key("AF").unwrap()),
        fo: k.fo.as_ref().map(|_| h.format_key("FO").unwrap()),
        st: k.st.as_ref().map(|_| h.format_key("ST").unwrap()),
    }
}
