//! The VCF/BCF reader on real files and crafted records: every field of the
//! BCF fixtures against bcftools, readers on many threads, plan errors,
//! foreign keys, lazy index errors and record limits.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss,
    reason = "test code"
)]

use hegel::prelude::*;
use seqair::io::BgzfWriter;
use seqair::vcf::reader::{
    BcfValueError, FieldError, FilterStatus, KeyKind, MAX_LINE_BYTES, MAX_RECORD_BYTES, ReadOption,
    RecordPart, RecordSection,
};
use seqair::vcf::{
    Allele, DynValue, Fields, FormatAny, FormatGt, FullRecord, Info, InfoAny, InfoFlag, InfoInt,
    KeyError, Locus, OnlySample, OpenError, PlanContext, PlanError, RawRecord, RecordError, Sample,
    Samples, Sections, Site, VcfFile, VcfPos, VcfReadError,
};
use std::path::PathBuf;
use std::process::Command;

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/vcf").join(name)
}

const BCFTOOLS: &str = "/opt/homebrew/bin/bcftools";

/// `bcftools args…` stdout, or `None` when bcftools is not installed.
fn bcftools(args: &[&str]) -> Option<String> {
    let bin = if std::path::Path::new(BCFTOOLS).exists() { BCFTOOLS } else { "bcftools" };
    match Command::new(bin).args(args).output() {
        Ok(out) => {
            assert!(
                out.status.success(),
                "bcftools {args:?}: {}",
                String::from_utf8(out.stderr.clone()).unwrap_or_default()
            );
            Some(String::from_utf8(out.stdout).unwrap())
        }
        Err(_) => {
            eprintln!("bcftools not found, skipping");
            None
        }
    }
}

// ── bcftools as the oracle ──────────────────────────────────────────────

/// A value as text renders it, floats kept as numbers to compare with the
/// precision bcftools prints them at.
#[derive(Debug, Clone, PartialEq)]
enum Tok {
    Int(Option<i32>),
    Float(Option<f32>),
    Str(String),
    Flag,
}

fn tokens_of(value: &DynValue<'_>) -> Vec<Tok> {
    match value {
        DynValue::Int(v) => v.clone().map(Tok::Int).collect(),
        DynValue::Float(v) => v.clone().map(Tok::Float).collect(),
        DynValue::Flag => vec![Tok::Flag],
        DynValue::Str(s) => vec![Tok::Str(s.unwrap_or(".").to_owned())],
        DynValue::Gt(None) => vec![Tok::Str(".".to_owned())],
        DynValue::Gt(Some(gt)) => {
            let mut s = String::new();
            for (i, a) in gt.alleles().enumerate() {
                if i > 0 {
                    s.push(if a.phased { '|' } else { '/' });
                }
                match a.allele {
                    Some(n) => s.push_str(&n.to_string()),
                    None => s.push('.'),
                }
            }
            vec![Tok::Str(s)]
        }
    }
}

/// `text` (one INFO value or sample value) parsed as `like` says.
fn tokens_from_text(text: Option<&str>, like: &DynValue<'_>) -> Vec<Tok> {
    let Some(text) = text else { return vec![Tok::Flag] };
    let missing_vector = text == ".";
    let split = || text.split(',');
    match like {
        DynValue::Int(_) if missing_vector => vec![],
        DynValue::Float(_) if missing_vector => vec![],
        DynValue::Int(_) => split().map(|t| Tok::Int(t.parse().ok())).collect(),
        DynValue::Float(_) => split().map(|t| Tok::Float(t.parse().ok())).collect(),
        DynValue::Flag => vec![Tok::Flag],
        // VCF 4.4 prints an explicit leading phase; the BCF fixtures are older.
        DynValue::Gt(_) => vec![Tok::Str(text.trim_start_matches(['|', '/']).to_owned())],
        DynValue::Str(_) => vec![Tok::Str(text.to_owned())],
    }
}

/// Our tokens equal bcftools', floats to the 6 significant digits it prints;
/// a lone missing element is how text spells an empty vector.
fn same(ours: &[Tok], theirs: &[Tok]) -> bool {
    let ours: Vec<Tok> = match ours {
        [Tok::Int(None)] | [Tok::Float(None)] => vec![],
        o => o.to_vec(),
    };
    ours.len() == theirs.len()
        && ours.iter().zip(theirs).all(|pair| match pair {
            (Tok::Float(Some(a)), Tok::Float(Some(b))) => close(*a, *b),
            (a, b) => a == b,
        })
}

fn close(a: f32, b: f32) -> bool {
    a == b || (a - b).abs() <= 1e-5 * a.abs().max(b.abs()) || (a.is_nan() && b.is_nan())
}

fn render_alt(a: Allele<'_>) -> String {
    String::from_utf8(a.as_bytes().to_vec()).unwrap()
}

// r[verify vcf_reader.oracles+2]
// r[verify vcf_values.access+2]
// r[verify vcf_values.dynamic+2]
// r[verify vcf_values.qual]
// r[verify vcf_fields.full_record+2]
/// Every field of every record of a BCF fixture, read through `FullRecord`,
/// equals what `bcftools view` prints for it; the end equals `%END`.
fn check_against_bcftools(name: &str) {
    let path = fixture(name);
    let path_str = path.to_str().unwrap();
    let Some(text) = bcftools(&["view", "-H", path_str]) else { return };
    let Some(ends) = bcftools(&["query", "-f", "%END\n", path_str]) else { return };
    let file = VcfFile::open(&path).unwrap();
    // Every fixture has an index: CSI for BCF, TBI for VCF.gz.
    assert!(file.has_index());
    let header = file.header();
    let samples: Vec<_> = header.sample_indices().collect();
    let mut reader = file.reader::<FullRecord>().open().unwrap();
    let mut lines = text.lines().zip(ends.lines());
    let mut n = 0;
    while let Some(rec) = reader.next_record().unwrap() {
        let (line, end) = lines.next().expect("bcftools has as many records");
        let cols: Vec<&str> = line.split('\t').collect();
        let at = format!("{name} record {n}: {line:.120}");
        let locus = rec.locus();
        assert_eq!(locus.contig().name(), cols[0], "{at}");
        assert_eq!(
            locus.start(),
            VcfPos::At(cols[1].parse::<u32>().ok().and_then(seqair_types::Pos1::new).unwrap()),
            "{at}"
        );
        assert_eq!(
            locus.end(),
            VcfPos::At(end.parse::<u32>().ok().and_then(seqair_types::Pos1::new).unwrap()),
            "{at}"
        );

        let site = rec.site();
        let ids: Vec<&str> = site.ids().collect();
        assert_eq!(if ids.is_empty() { ".".to_owned() } else { ids.join(";") }, cols[2], "{at}");
        assert_eq!(site.alleles().reference(), cols[3].as_bytes(), "{at}");
        let alts: Vec<String> = site.alleles().alts().map(render_alt).collect();
        assert_eq!(if alts.is_empty() { ".".to_owned() } else { alts.join(",") }, cols[4], "{at}");
        match (site.qual(), cols[5]) {
            (None, ".") => {}
            (Some(q), t) => assert!(close(q, t.parse().unwrap()), "{at}: QUAL {q} vs {t}"),
            (q, t) => panic!("{at}: QUAL {q:?} vs {t}"),
        }
        let filter = match site.filters().status() {
            FilterStatus::Pass => "PASS".to_owned(),
            FilterStatus::Missing => ".".to_owned(),
            FilterStatus::Failed(failed) => {
                failed.map(|f| f.name().to_owned()).collect::<Vec<_>>().join(";")
            }
        };
        assert_eq!(filter, cols[6], "{at}");

        let info: Vec<_> = rec.info().iter().collect::<Result<_, _>>().unwrap();
        let theirs: Vec<(&str, Option<&str>)> = if cols[7] == "." {
            vec![]
        } else {
            cols[7]
                .split(';')
                .map(|kv| kv.split_once('=').map_or((kv, None), |(k, v)| (k, Some(v))))
                .collect()
        };
        assert_eq!(info.len(), theirs.len(), "{at}: INFO fields");
        for ((key, _, value), (their_key, their_value)) in info.iter().zip(&theirs) {
            assert_eq!(key, their_key, "{at}");
            let (ours, theirs) = (tokens_of(value), tokens_from_text(*their_value, value));
            assert!(same(&ours, &theirs), "{at}: INFO {key}: {ours:?} vs {theirs:?}");
        }

        if !samples.is_empty() {
            let keys: Vec<&str> = cols[8].split(':').collect();
            for (s, &sample) in samples.iter().enumerate() {
                let values: Vec<&str> = cols[9 + s].split(':').collect();
                let ours: Vec<_> =
                    rec.samples().sample(sample).iter().collect::<Result<_, _>>().unwrap();
                assert_eq!(ours.len(), keys.len(), "{at}: sample {s} FORMAT fields");
                for (i, (key, _, value)) in ours.iter().enumerate() {
                    assert_eq!(*key, keys[i], "{at}");
                    let theirs =
                        tokens_from_text(Some(values.get(i).copied().unwrap_or(".")), value);
                    let ours = tokens_of(value);
                    assert!(same(&ours, &theirs), "{at}: sample {s} {key}: {ours:?} vs {theirs:?}");
                }
            }
        }
        n += 1;
    }
    assert!(lines.next().is_none(), "{name}: bcftools has more records than were read");
    assert!(n > 0);
}

#[test]
fn rastair_bcf_reads_as_bcftools_prints_it() {
    check_against_bcftools("rastair.bcf");
}

#[test]
fn panel_bcf_reads_as_bcftools_prints_it() {
    check_against_bcftools("panel.bcf");
}

#[test]
fn rastair_vcf_gz_reads_as_bcftools_prints_it() {
    check_against_bcftools("rastair.vcf.gz");
}

#[test]
fn truth_vcf_gz_reads_as_bcftools_prints_it() {
    check_against_bcftools("truth.vcf.gz");
}

// ── noodles as a second implementation ──────────────────────────────────

fn noodles_strs<T: ToString>(v: &[Option<T>]) -> Vec<Tok> {
    let joined: Vec<String> =
        v.iter().map(|x| x.as_ref().map_or_else(|| ".".to_owned(), ToString::to_string)).collect();
    vec![Tok::Str(joined.join(","))]
}

fn noodles_info(v: Option<&noodles::vcf::variant::record_buf::info::field::Value>) -> Vec<Tok> {
    use noodles::vcf::variant::record_buf::info::field::{Value, value::Array};
    match v {
        None => vec![],
        Some(Value::Integer(i)) => vec![Tok::Int(Some(*i))],
        Some(Value::Float(f)) => vec![Tok::Float(Some(*f))],
        Some(Value::Flag) => vec![Tok::Flag],
        Some(Value::Character(c)) => vec![Tok::Str(c.to_string())],
        Some(Value::String(s)) => vec![Tok::Str(s.clone())],
        Some(Value::Array(Array::Integer(v))) => v.iter().copied().map(Tok::Int).collect(),
        Some(Value::Array(Array::Float(v))) => v.iter().copied().map(Tok::Float).collect(),
        Some(Value::Array(Array::Character(v))) => noodles_strs(v),
        Some(Value::Array(Array::String(v))) => noodles_strs(v),
    }
}

fn noodles_sample(
    v: Option<&noodles::vcf::variant::record_buf::samples::sample::Value>,
) -> Vec<Tok> {
    use noodles::vcf::variant::record::samples::series::value::genotype::Phasing;
    use noodles::vcf::variant::record_buf::samples::sample::{Value, value::Array};
    match v {
        None => vec![],
        Some(Value::Integer(i)) => vec![Tok::Int(Some(*i))],
        Some(Value::Float(f)) => vec![Tok::Float(Some(*f))],
        Some(Value::Character(c)) => vec![Tok::Str(c.to_string())],
        Some(Value::String(s)) => vec![Tok::Str(s.clone())],
        Some(Value::Genotype(g)) => {
            let mut s = String::new();
            for (i, a) in g.as_ref().iter().enumerate() {
                if i > 0 {
                    s.push(if a.phasing() == Phasing::Phased { '|' } else { '/' });
                }
                s.push_str(&a.position().map_or_else(|| ".".to_owned(), |p| p.to_string()));
            }
            vec![Tok::Str(s)]
        }
        Some(Value::Array(Array::Integer(v))) => v.iter().copied().map(Tok::Int).collect(),
        Some(Value::Array(Array::Float(v))) => v.iter().copied().map(Tok::Float).collect(),
        Some(Value::Array(Array::Character(v))) => noodles_strs(v),
        Some(Value::Array(Array::String(v))) => noodles_strs(v),
    }
}

/// A missing GT is `.` to us and no value to noodles.
fn same_or_missing(ours: &[Tok], theirs: &[Tok]) -> bool {
    same(ours, theirs) || (ours == [Tok::Str(".".to_owned())] && theirs.is_empty())
}

// r[verify vcf_reader.oracles+2]
/// Every field of every record of a VCF text fixture, read through
/// `FullRecord`, equals what noodles reads.
fn check_against_noodles(name: &str) {
    let path = fixture(name);
    let file = VcfFile::open(&path).unwrap();
    let samples: Vec<_> = file.header().sample_indices().collect();
    let bgzf = noodles::bgzf::io::Reader::new(std::fs::File::open(&path).unwrap());
    let mut theirs = noodles::vcf::io::Reader::new(std::io::BufReader::new(bgzf));
    let their_header = theirs.read_header().unwrap();
    let mut their_records = theirs.record_bufs(&their_header);
    let mut reader = file.reader::<FullRecord>().open().unwrap();
    let mut n = 0;
    while let Some(rec) = reader.next_record().unwrap() {
        let t = their_records.next().expect("noodles has as many records").unwrap();
        let at = format!("{name} record {n}");
        assert_eq!(rec.locus().contig().name(), t.reference_sequence_name(), "{at}");
        let pos = t.variant_start().map_or(0, |p| u32::try_from(usize::from(p)).unwrap());
        let start = match rec.locus().start() {
            VcfPos::At(p) => p.as_u32(),
            VcfPos::Telomere => 0,
        };
        assert_eq!(start, pos, "{at}");
        let site = rec.site();
        let ids: Vec<&str> = site.ids().collect();
        let their_ids: Vec<&str> = t.ids().as_ref().iter().map(String::as_str).collect();
        assert_eq!(ids, their_ids, "{at}");
        assert_eq!(site.alleles().reference(), t.reference_bases().as_bytes(), "{at}");
        let alts: Vec<String> = site.alleles().alts().map(render_alt).collect();
        assert_eq!(alts, t.alternate_bases().as_ref(), "{at}");
        match (site.qual(), t.quality_score()) {
            (None, None) => {}
            (Some(a), Some(b)) => assert!(close(a, b), "{at}: QUAL {a} vs {b}"),
            (a, b) => panic!("{at}: QUAL {a:?} vs {b:?}"),
        }
        let filters: Vec<String> = match site.filters().status() {
            FilterStatus::Pass => vec!["PASS".to_owned()],
            FilterStatus::Missing => vec![],
            FilterStatus::Failed(f) => f.map(|f| f.name().to_owned()).collect(),
        };
        let their_filters: Vec<String> = t.filters().as_ref().iter().cloned().collect();
        assert_eq!(filters, their_filters, "{at}");
        let info: Vec<_> = rec.info().iter().collect::<Result<_, _>>().unwrap();
        let their_info = t.info().as_ref();
        assert_eq!(info.len(), their_info.len(), "{at}: INFO fields");
        for ((key, _, value), (their_key, their_value)) in info.iter().zip(their_info) {
            assert_eq!(key, their_key, "{at}");
            let (ours, theirs) = (tokens_of(value), noodles_info(their_value.as_ref()));
            assert!(same(&ours, &theirs), "{at}: INFO {key}: {ours:?} vs {theirs:?}");
        }
        let keys: Vec<&String> = t.samples().keys().as_ref().iter().collect();
        for (sample, their_sample) in samples.iter().zip(t.samples().values()) {
            let ours: Vec<_> =
                rec.samples().sample(*sample).iter().collect::<Result<_, _>>().unwrap();
            assert_eq!(ours.len(), keys.len(), "{at}: FORMAT fields");
            for ((key, _, value), (their_key, their_value)) in
                ours.iter().zip(keys.iter().zip(their_sample.values()))
            {
                assert_eq!(key, their_key, "{at}");
                let (ours, theirs) = (tokens_of(value), noodles_sample(their_value.as_ref()));
                assert!(same_or_missing(&ours, &theirs), "{at}: {key}: {ours:?} vs {theirs:?}");
            }
        }
        n += 1;
    }
    assert!(their_records.next().is_none(), "{name}: noodles has more records");
    assert!(n > 0);
}

// rastair.vcf.gz is not compared: noodles refuses its `Number=M` FORMAT field.
#[test]
fn truth_vcf_gz_reads_as_noodles_reads_it() {
    check_against_noodles("truth.vcf.gz");
}

// ── Fixtures through every way of reading ──────────────────────────────

// r[verify vcf_reader.reader]
// r[verify vcf_reader.read_into+2]
// r[verify vcf_reader.next_record]
// r[verify vcf_reader.records]
// r[verify vcf_reader.passthrough_header]
// r[verify vcf_fields.owned+2]
#[test]
fn every_way_of_reading_yields_the_same_records() {
    let file = VcfFile::open(fixture("rastair.bcf")).unwrap();
    let owned: Vec<Locus> =
        file.reader::<Locus>().open().unwrap().records().collect::<Result<_, _>>().unwrap();
    assert_eq!(owned.len(), 400);
    let into: Vec<Locus> =
        file.reader::<Locus>().open().unwrap().into_iter().collect::<Result<_, _>>().unwrap();
    assert_eq!(owned, into);
    let mut reader = file.reader::<Locus>().open().unwrap();
    let mut rec = reader.read().unwrap().unwrap();
    let mut lent = vec![rec.clone()];
    while reader.read_into(&mut rec).unwrap() {
        lent.push(rec.clone());
    }
    assert_eq!(owned, lent);
    let last = rec.clone();
    assert!(!reader.read_into(&mut rec).unwrap(), "stays at the end");
    assert_eq!(rec, last, "the end leaves the record as it was");
    assert!(reader.read().unwrap().is_none());
    let mut reader = file.reader::<Locus>().open().unwrap();
    let mut lent_by_reader = Vec::new();
    while let Some(locus) = reader.next_record().unwrap() {
        lent_by_reader.push(locus.clone());
    }
    assert_eq!(owned, lent_by_reader);
    assert!(reader.next_record().unwrap().is_none(), "stays at the end");
    assert_eq!(reader.header().id(), file.header().id());
}

/// A record's contents as one comparable value.
fn digest(rec: &FullRecord) -> (Locus, Vec<u8>) {
    (rec.locus().clone(), rec.raw_bytes().to_vec())
}

// r[verify vcf_reader.fork_equivalence]
// r[verify vcf_reader.file.shared+2]
// r[verify vcf_reader.reader.threads+2]
/// Readers of one `VcfFile` (and of clones of it) on several threads at
/// once, some decompressing on a thread pool, read what one reader reads.
#[test]
fn readers_on_many_threads_read_the_same_records() {
    for name in ["rastair.bcf", "panel.bcf"] {
        let file = VcfFile::open(fixture(name)).unwrap();
        let expected: Vec<_> = file
            .reader::<FullRecord>()
            .open()
            .unwrap()
            .into_iter()
            .map(|r| digest(&r.unwrap()))
            .collect();
        let pool = std::sync::Arc::new(
            rayon_core::ThreadPoolBuilder::new().num_threads(2).build().unwrap(),
        );
        std::thread::scope(|s| {
            let handles: Vec<_> = (0..6)
                .map(|i| {
                    let file = file.clone();
                    let pool = std::sync::Arc::clone(&pool);
                    s.spawn(move || {
                        let two = std::num::NonZeroUsize::new(2).unwrap();
                        let builder = file.reader::<FullRecord>();
                        let builder = match i % 3 {
                            0 => builder,
                            1 => builder.threads(two),
                            _ => builder.threads_in(two, pool),
                        };
                        builder
                            .open()
                            .unwrap()
                            .into_iter()
                            .map(|r| digest(&r.unwrap()))
                            .collect::<Vec<_>>()
                    })
                })
                .collect();
            for h in handles {
                assert!(h.join().unwrap() == expected, "{name}: a thread read different records");
            }
        });
    }
}

#[test]
fn files_and_readers_cross_threads() {
    fn send_sync_clone<T: Send + Sync + Clone>() {}
    fn send_clone<T: Send + Clone>() {}
    fn send<T: Send>() {}
    send_sync_clone::<VcfFile>();
    send::<seqair::vcf::Reader<FullRecord>>();
    // r[verify vcf_fields.owned+2]
    send_clone::<FullRecord>();
}

// r[verify vcf_fields.sample+2]
// r[verify vcf_fields.samples+2]
// r[verify vcf_fields.only_sample]
#[test]
fn one_sample_reads_what_all_samples_read_for_it() {
    let file = VcfFile::open(fixture("panel.bcf")).unwrap();
    let gt: FormatGt = file.header().format_key("GT").unwrap();
    for name in file.header().samples().to_vec() {
        let idx = file.header().sample(&name).unwrap();
        let one = file.reader::<Sample>().sample(name.clone()).open().unwrap();
        let all = file.reader::<Samples>().open().unwrap();
        for (one, all) in one.into_iter().zip(all) {
            let (one, all) = (one.unwrap(), all.unwrap());
            assert_eq!(one.name(), name);
            let a = one.get(&gt).unwrap().map(|g| g.to_genotype());
            let b = all.sample(idx).get(&gt).unwrap().map(|g| g.to_genotype());
            let c = all.column(&gt).unwrap().unwrap().get(idx).unwrap().map(|g| g.to_genotype());
            assert_eq!(a, b);
            assert_eq!(a, c);
            assert_eq!(all.len(), 4);
        }
    }
    let file = VcfFile::open(fixture("rastair.bcf")).unwrap();
    let gt: FormatGt = file.header().format_key("GT").unwrap();
    let only = file.reader::<OnlySample>().open().unwrap();
    let sample = file.reader::<Sample>().open().unwrap();
    for (a, b) in only.into_iter().zip(sample) {
        let (a, b) = (a.unwrap(), b.unwrap());
        assert_eq!(
            a.get(&gt).unwrap().map(|g| g.to_genotype()),
            b.get(&gt).unwrap().map(|g| g.to_genotype())
        );
    }
}

// r[verify vcf_fields.tuples]
/// A tuple reads what its elements read alone, in any order, duplicates
/// included.
#[test]
fn tuples_read_what_their_elements_read() {
    let file = VcfFile::open(fixture("rastair.bcf")).unwrap();
    let dp: InfoInt = file.header().info_key("DP").unwrap();
    let loci: Vec<Locus> =
        file.reader::<Locus>().open().unwrap().into_iter().map(Result::unwrap).collect();
    let depths: Vec<Option<i32>> = file
        .reader::<Info>()
        .open()
        .unwrap()
        .into_iter()
        .map(|i| i.unwrap().get(&dp).unwrap())
        .collect();
    let mut reader = file.reader::<(Info, Locus, Info)>().open().unwrap();
    let mut n = 0;
    while let Some((a, locus, b)) = reader.next_record().unwrap() {
        assert_eq!(*locus, loci[n]);
        assert_eq!(a.get(&dp).unwrap(), depths[n]);
        assert_eq!(b.get(&dp).unwrap(), depths[n]);
        n += 1;
    }
    assert_eq!(n, loci.len());
    let needs = <(Locus, Site, Info) as Fields>::needs(&(file.header().id(), (), ()));
    assert_eq!(needs, Sections::SITE | Sections::INFO);
}

// ── Plans ───────────────────────────────────────────────────────────────

fn plan_error<F: Fields>(file: &VcfFile, sample: Option<&str>) -> PlanError {
    let builder = file.reader::<F>();
    let builder = match sample {
        Some(s) => builder.sample(s),
        None => builder,
    };
    match builder.open() {
        Err(VcfReadError::Plan(e)) => e,
        other => panic!("expected a plan error, got {other:?}"),
    }
}

// r[verify vcf_fields.plan]
// r[verify vcf_fields.unused_options+2]
// r[verify vcf_fields.only_sample]
// r[verify vcf_fields.sample+2]
#[test]
fn plans_fail_at_open() {
    let panel = VcfFile::open(fixture("panel.bcf")).unwrap();
    let rastair = VcfFile::open(fixture("rastair.bcf")).unwrap();
    assert_eq!(plan_error::<OnlySample>(&panel, None), PlanError::NotSingleSample { n: 4 });
    assert_eq!(plan_error::<Sample>(&panel, None), PlanError::SampleNotSelected { n: 4 });
    assert_eq!(
        plan_error::<Sample>(&panel, Some("nope")),
        PlanError::Key(KeyError::UnknownSample { name: "nope".into() })
    );
    assert_eq!(
        plan_error::<Locus>(&panel, Some("NA12878")),
        PlanError::UnusedOption { option: ReadOption::Sample }
    );
    assert_eq!(
        plan_error::<(Locus, OnlySample)>(&rastair, Some("x")),
        PlanError::UnusedOption { option: ReadOption::Sample }
    );
    // A tuple plans each element; one that consumes the option is enough.
    assert!(panel.reader::<(Locus, Sample, Info)>().sample("NA12878").open().is_ok());
    assert!(matches!(
        plan_error::<BadKey>(&rastair, None),
        PlanError::Key(KeyError::TypeMismatch { .. })
    ));
    // A sample by handle: this header's reads, another header's is refused.
    let na12878 = panel.header().sample("NA12878").unwrap();
    let mut reader = panel.reader::<Sample>().sample(na12878).open().unwrap();
    assert_eq!(reader.next_record().unwrap().unwrap().name(), "NA12878");
    let foreign = VcfFile::open(fixture("panel.bcf")).unwrap().header().sample("NA12878").unwrap();
    assert!(matches!(
        panel.reader::<Sample>().sample(foreign).open(),
        Err(VcfReadError::Plan(PlanError::ForeignKey { .. }))
    ));
    // A plan given by the caller skips `plan`, so an option is left unused.
    let plan = (na12878, "NA12878".into());
    assert!(panel.reader::<Sample>().with_plan(plan.clone()).open().is_ok());
    assert!(matches!(
        panel.reader::<Sample>().with_plan(plan).sample("NA12878").open(),
        Err(VcfReadError::Plan(PlanError::UnusedOption { option: ReadOption::Sample }))
    ));
}

/// A user type with runtime context and one sample's FORMAT: its plan is
/// built by the caller (`with_plan`), or by `plan` through the context.
#[derive(Clone, Debug)]
struct Depth(Option<bool>);

struct DepthPlan {
    sample: seqair::vcf::SampleIdx,
    gt: FormatGt,
    min_alt: u16,
}

impl Fields for Depth {
    type Plan = DepthPlan;
    fn plan(cx: &PlanContext<'_>) -> Result<DepthPlan, PlanError> {
        Ok(DepthPlan { sample: cx.selected_sample()?, gt: cx.format_key("GT")?, min_alt: 1 })
    }
    fn needs(plan: &DepthPlan) -> Sections {
        Sections::sample(plan.sample)
    }
    fn new(plan: &DepthPlan, raw: RawRecord<'_>) -> Result<Self, VcfReadError> {
        let gt = raw.sample(plan.sample)?.get(&plan.gt)?;
        Ok(Self(gt.map(|gt| gt.alleles().any(|a| a.allele >= Some(plan.min_alt)))))
    }
}

// r[verify vcf_fields.plan_context]
// r[verify vcf_fields.with_plan]
#[test]
fn a_user_type_plans_with_the_context_or_takes_a_plan() {
    let panel = VcfFile::open(fixture("panel.bcf")).unwrap();
    assert_eq!(plan_error::<Depth>(&panel, None), PlanError::SampleNotSelected { n: 4 });
    let read = |reader: seqair::vcf::Reader<Depth>| -> Vec<Option<bool>> {
        reader.into_iter().map(|d| d.unwrap().0).collect()
    };
    let planned = read(panel.reader::<Depth>().sample("NA12878").open().unwrap());
    let sample = panel.header().sample("NA12878").unwrap();
    let gt = panel.format_key("GT").unwrap();
    let given =
        read(panel.reader().with_plan(DepthPlan { sample, gt, min_alt: 1 }).open().unwrap());
    assert_eq!(planned, given);
    assert!(planned.contains(&Some(true)));
    let gt = panel.format_key("GT").unwrap();
    let strict =
        read(panel.reader().with_plan(DepthPlan { sample, gt, min_alt: 100 }).open().unwrap());
    assert!(strict.iter().all(|d| *d != Some(true)), "runtime context reaches `new`");
}

#[derive(Clone, Debug)]
struct BadKey;

impl Fields for BadKey {
    type Plan = InfoInt;
    fn plan(cx: &PlanContext<'_>) -> Result<InfoInt, PlanError> {
        Ok(cx.info_key("BQ")?) // declared Float
    }
    fn needs(_: &InfoInt) -> Sections {
        Sections::INFO
    }
    fn new(_: &InfoInt, _: RawRecord<'_>) -> Result<Self, VcfReadError> {
        Ok(Self)
    }
}

// ── Keys of another header ──────────────────────────────────────────────

// r[verify vcf_values.header_brand+2]
// r[verify vcf_fields.sample_idx]
#[test]
fn keys_and_samples_of_another_header_are_refused() {
    let rastair = VcfFile::open(fixture("rastair.bcf")).unwrap();
    let panel = VcfFile::open(fixture("panel.bcf")).unwrap();
    // Same file opened twice: same text, different headers.
    let again = VcfFile::open(fixture("panel.bcf")).unwrap();
    let foreign_gt: FormatGt = again.header().format_key("GT").unwrap();
    let foreign_sample = again.header().sample("NA12878").unwrap();
    let foreign_dp: InfoInt = rastair.header().info_key("DP").unwrap();
    let foreign_filter = rastair.header().filter("lowDp").unwrap();

    let mut reader = panel.reader::<(Site, Info, Samples)>().open().unwrap();
    let (site, info, samples) = reader.next_record().unwrap().unwrap();
    let is_foreign = |e: VcfReadError| matches!(e, VcfReadError::ForeignKey { .. });
    assert!(is_foreign(info.get(&foreign_dp).unwrap_err()));
    assert!(is_foreign(samples.column(&foreign_gt).unwrap_err()));
    let own_gt: FormatGt = panel.header().format_key("GT").unwrap();
    assert!(is_foreign(samples.sample(foreign_sample).get(&own_gt).unwrap_err()));
    assert!(is_foreign(samples.column(&own_gt).unwrap().unwrap().get(foreign_sample).unwrap_err()));
    assert!(is_foreign(site.filters().contains(&foreign_filter).unwrap_err()));
    // PASS is index 0 in every header.
    assert!(!site.filters().contains(&seqair::vcf::FilterId::PASS).unwrap());
}

// ── Your own record type ────────────────────────────────────────────────

/// The reading guide's `Truth` type. (The guide's example asserts a `CpG` SNV
/// exists, which this fixture does not have: its `CpG` sites have no ALT.)
#[derive(Clone, Debug)]
struct Truth(Option<(VcfPos, u8, u8, bool)>);

impl Fields for Truth {
    type Plan = InfoFlag;
    fn plan(cx: &PlanContext<'_>) -> Result<InfoFlag, PlanError> {
        Ok(cx.info_key("CPG")?)
    }
    fn needs(_: &InfoFlag) -> Sections {
        Sections::SITE | Sections::INFO
    }
    fn new(cpg: &InfoFlag, raw: RawRecord<'_>) -> Result<Self, VcfReadError> {
        let site = raw.site()?;
        let mut alts = site.alleles().alts();
        Ok(Self(
            match (site.filters().status(), site.alleles().reference(), alts.next(), alts.next()) {
                (FilterStatus::Pass, &[r], Some(Allele::Bases(&[a])), None) => {
                    Some((raw.locus()?.start(), r, a, raw.info()?.get(cpg)?))
                }
                _ => None,
            },
        ))
    }
}

// r[verify vcf_fields.trait+4]
// r[verify vcf_fields.user_types+2]
// r[verify vcf_fields.sections]
#[test]
fn a_user_type_reads_through_raw_records() {
    let file = VcfFile::open(fixture("rastair.bcf")).unwrap();
    let snvs: Vec<_> =
        file.reader::<Truth>().open().unwrap().into_iter().filter_map(|t| t.unwrap().0).collect();
    let pos = VcfPos::At(seqair_types::Pos1::new(10437).unwrap());
    assert_eq!(snvs, [(pos, b'C', b'T', false)]);

    // Viewing a section `needs` did not ask for is an error, not a silent read.
    #[derive(Clone, Debug)]
    struct Sneaky;
    impl Fields for Sneaky {
        type Plan = ();
        fn plan(_: &PlanContext<'_>) -> Result<(), PlanError> {
            Ok(())
        }
        fn needs(_: &()) -> Sections {
            Sections::SITE
        }
        fn new(_: &(), raw: RawRecord<'_>) -> Result<Self, VcfReadError> {
            raw.info().map(|_| Self)
        }
    }
    let err = file.reader::<Sneaky>().open().unwrap().next_record().unwrap_err();
    assert!(matches!(err, VcfReadError::NotRequested { section: RecordSection::Info }));
}

// ── Crafted BCF ─────────────────────────────────────────────────────────

const HEADER: &str = "##fileformat=VCFv4.3\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
##contig=<ID=chr1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n";

/// One BCF record, valid unless a test changes a part: chr1:11 A, PASS,
/// DP=5, GT 0/1.
#[derive(Clone)]
struct Rec {
    rid: i32,
    pos: i32,
    n_allele: u16,
    site: Vec<u8>,
    info: Vec<u8>,
    n_info: u16,
    indiv: Vec<u8>,
    n_fmt: u8,
    n_sample: u32,
    /// Overrides the record's length fields.
    lengths: Option<(u32, u32)>,
}

impl Default for Rec {
    fn default() -> Self {
        Self {
            rid: 0,
            pos: 10,
            n_allele: 1,
            site: vec![0x17, b'.', 0x17, b'A', 0x11, 0x00],
            info: vec![0x11, 0x01, 0x11, 5],
            n_info: 1,
            indiv: vec![0x11, 0x02, 0x21, 0x02, 0x04],
            n_fmt: 1,
            n_sample: 1,
            lengths: None,
        }
    }
}

impl Rec {
    fn bytes(&self) -> Vec<u8> {
        let mut shared = Vec::new();
        for v in [self.rid, self.pos, 1] {
            shared.extend_from_slice(&v.to_le_bytes());
        }
        shared.extend_from_slice(&0x7F80_0001u32.to_le_bytes());
        shared.extend_from_slice(
            &(u32::from(self.n_info) | (u32::from(self.n_allele) << 16)).to_le_bytes(),
        );
        shared.extend_from_slice(&(self.n_sample | (u32::from(self.n_fmt) << 24)).to_le_bytes());
        shared.extend_from_slice(&self.site);
        shared.extend_from_slice(&self.info);
        let (l_shared, l_indiv) =
            self.lengths.unwrap_or((shared.len() as u32, self.indiv.len() as u32));
        let mut out = Vec::new();
        out.extend_from_slice(&l_shared.to_le_bytes());
        out.extend_from_slice(&l_indiv.to_le_bytes());
        out.extend_from_slice(&shared);
        out.extend_from_slice(&self.indiv);
        out
    }
}

fn bcf(records: &[Rec]) -> Vec<u8> {
    bcf_with(HEADER, records)
}

fn bcf_with(header: &str, records: &[Rec]) -> Vec<u8> {
    let mut raw = b"BCF\x02\x02".to_vec();
    let text = format!("{header}\0");
    raw.extend_from_slice(&(text.len() as u32).to_le_bytes());
    raw.extend_from_slice(text.as_bytes());
    for r in records {
        raw.extend_from_slice(&r.bytes());
    }
    let mut w = BgzfWriter::new(Vec::new());
    w.write_all(&raw).unwrap();
    w.finish().unwrap()
}

fn first_error<F: Fields>(records: &[Rec]) -> VcfReadError {
    let file = VcfFile::from_bytes(bcf(records), None).unwrap();
    let mut reader = file.reader::<F>().open().unwrap();
    loop {
        match reader.next_record() {
            Ok(Some(_)) => {}
            Ok(None) => panic!("no error"),
            Err(e) => return e,
        }
    }
}

fn record_error(e: VcfReadError) -> (u64, RecordError) {
    match e {
        VcfReadError::Record { at, source } => {
            assert!(at.offset.is_some(), "a BGZF record error names its offset");
            assert_eq!(at.line, None, "BCF has no lines");
            (at.number.unwrap(), source)
        }
        e => panic!("expected a record error, got {e:?}"),
    }
}

#[test]
fn crafted_records_read_back() {
    let file = VcfFile::from_bytes(bcf(&[Rec::default(), Rec::default()]), None).unwrap();
    assert!(!file.has_index());
    let dp: InfoInt = file.header().info_key("DP").unwrap();
    let recs: Vec<FullRecord> =
        file.reader::<FullRecord>().open().unwrap().into_iter().collect::<Result<_, _>>().unwrap();
    assert_eq!(recs.len(), 2);
    assert_eq!(recs[0].info().get(&dp).unwrap(), Some(5));
    assert_eq!(recs[0].site().filters().status(), FilterStatus::Pass);
    assert_eq!(recs[0].locus().start(), VcfPos::At(seqair_types::Pos1::new(11).unwrap()));
}

// r[verify vcf_fields.site+2]
/// htslib writes a text ALT `.` as one `.` allele after REF: no ALT, as the
/// text reads; a `.` among ALTs keeps its place as `Missing`.
#[test]
fn a_bcf_alt_dot_is_no_alt() {
    let site = |alts: &[&[u8]]| {
        let mut site = vec![0x17, b'.', 0x17, b'A'];
        for alt in alts {
            site.push(0x07 | ((alt.len() as u8) << 4));
            site.extend_from_slice(alt);
        }
        site.extend_from_slice(&[0x11, 0x00]);
        Rec { n_allele: 1 + alts.len() as u16, site, ..Rec::default() }
    };
    let file =
        VcfFile::from_bytes(bcf(&[site(&[b"."]), site(&[b"T", b"."]), site(&[])]), None).unwrap();
    let sites: Vec<Site> =
        file.reader::<Site>().open().unwrap().into_iter().collect::<Result<_, _>>().unwrap();
    fn alts(s: &Site) -> Vec<Allele<'_>> {
        s.alleles().alts().collect()
    }
    assert_eq!(alts(&sites[0]), []);
    assert_eq!(sites[0].alleles().len(), 1);
    assert_eq!(alts(&sites[1]), [Allele::Bases(b"T"), Allele::Missing]);
    assert_eq!(sites[1].alleles().get(2), Some(Allele::Missing));
    assert_eq!(alts(&sites[2]), []);
    let text = text_records(&["chr1\t5\t.\tA\t.\t.\t.\t.", "chr1\t6\t.\tA\tT,.\t.\t.\t."]);
    assert_eq!(text[0].site().alleles().alts().count(), 0);
    let all: Vec<_> = text[1].site().alleles().iter().collect();
    assert_eq!(all, [Allele::Bases(b"A"), Allele::Bases(b"T"), Allele::Missing]);
    assert_eq!(text[1].site().alleles().get(0), Some(Allele::Bases(b"A")));
    assert_eq!(text[1].site().alleles().get(3), None);
}

// r[verify vcf_reader.limits+2]
// r[verify vcf_reader.errors+3]
// r[verify vcf_reader.field_table]
#[test]
fn malformed_records_are_typed_errors_with_their_place() {
    let ok = Rec::default();
    let (ordinal, e) =
        record_error(first_error::<Locus>(&[ok.clone(), Rec { rid: 3, ..ok.clone() }]));
    assert_eq!((ordinal, e), (2, RecordError::ContigIndexOutOfRange { index: 3 }));
    let (_, e) = record_error(first_error::<Locus>(&[Rec { pos: -2, ..ok.clone() }]));
    assert_eq!(e, RecordError::InvalidPosition { pos: -2 });
    let huge = Rec { lengths: Some((u32::MAX, 16)), ..ok.clone() };
    let (_, e) = record_error(first_error::<Locus>(&[huge]));
    assert!(matches!(e, RecordError::TooLarge { max: MAX_RECORD_BYTES, .. }));
    let short = Rec { lengths: Some((20, 0)), ..ok.clone() };
    assert!(matches!(
        record_error(first_error::<Locus>(&[short])).1,
        RecordError::SharedTooShort { len: 20 }
    ));
    // A record that claims more bytes than the file has.
    let cut = Rec { lengths: Some((1000, 0)), ..ok.clone() };
    assert!(matches!(record_error(first_error::<Locus>(&[cut])).1, RecordError::Truncated));
    // An allele string that runs past the shared section.
    let long_id = Rec { site: vec![0x17, b'.', 0x77, b'A'], info: vec![], n_info: 0, ..ok.clone() };
    assert!(matches!(record_error(first_error::<Site>(&[long_id])).1, RecordError::Value { .. }));
    // A FILTER index that names an INFO field.
    let bad_filter = Rec { site: vec![0x17, b'.', 0x17, b'A', 0x11, 0x01], ..ok.clone() };
    let e = record_error(first_error::<Site>(&[bad_filter])).1;
    assert_eq!(e, RecordError::UnknownKey { index: 1, kind: KeyKind::Filter });
    // A ploidy above 256.
    let mut gt = vec![0x11, 0x02, 0xF1, 0x12];
    gt.extend_from_slice(&300u16.to_le_bytes());
    gt.extend(std::iter::repeat_n(0x02, 300));
    let file = VcfFile::from_bytes(bcf(&[Rec { indiv: gt, ..ok }]), None).unwrap();
    let key: FormatGt = file.header().format_key("GT").unwrap();
    let mut reader = file.reader::<OnlySample>().open().unwrap();
    let sample = reader.next_record().unwrap().unwrap();
    let e = record_error(sample.get(&key).unwrap_err()).1;
    assert!(matches!(
        e,
        RecordError::Field { source: FieldError::Ploidy { ploidy: 300, max: 256 }, .. }
    ));
}

// r[verify vcf_reader.limits+2]
/// A BCF record that cannot be framed fails every later read with the same
/// error, never the end of input; a framed malformed record does not.
#[test]
fn a_record_that_cannot_be_framed_fails_every_later_read() {
    let ok = Rec::default();
    let huge = Rec { lengths: Some((u32::MAX, 16)), ..ok.clone() };
    let short = Rec { lengths: Some((20, 0)), ..ok.clone() };
    let cut = Rec { lengths: Some((1000, 0)), ..ok.clone() };
    for bad in [huge, short, cut] {
        let file = VcfFile::from_bytes(bcf(&[ok.clone(), bad, ok.clone()]), None).unwrap();
        let mut reader = file.reader::<Locus>().open().unwrap();
        assert!(reader.next_record().unwrap().is_some());
        let first = record_error(reader.next_record().unwrap_err());
        for _ in 0..3 {
            assert_eq!(record_error(reader.next_record().unwrap_err()), first);
        }
    }
    // A framed record with a contig the header lacks: the next one reads.
    let file = VcfFile::from_bytes(bcf(&[Rec { rid: 3, ..ok.clone() }, ok]), None).unwrap();
    let mut reader = file.reader::<Locus>().open().unwrap();
    assert!(reader.next_record().is_err());
    assert!(reader.next_record().unwrap().is_some());
}

// r[verify vcf_reader.limits+2]
/// A text line over the limit, read from a stream that never holds it, fails
/// every later read with the same error instead of reading the line's tail
/// as a record.
#[test]
fn a_text_line_over_the_limit_fails_every_later_read() {
    /// The header, a record whose INFO runs past the limit, and a record.
    struct Long(u64);
    impl std::io::Read for Long {
        fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
            let head = format!("{TEXT_HEADER}chr1\t5\t.\tA\tC\t.\t.\tS=");
            let tail = "\nchr1\t6\t.\tA\tC\t.\t.\t.\n";
            let (head_len, body) = (head.len() as u64, MAX_LINE_BYTES + 1);
            let at = self.0;
            let n = if at < head_len {
                let rest = &head.as_bytes()[at as usize..];
                let n = rest.len().min(buf.len());
                buf[..n].copy_from_slice(&rest[..n]);
                n
            } else if at < head_len + body {
                let n = ((head_len + body - at) as usize).min(buf.len());
                buf[..n].fill(b'x');
                n
            } else {
                let rest =
                    tail.as_bytes().get((at - head_len - body) as usize..).unwrap_or_default();
                let n = rest.len().min(buf.len());
                buf[..n].copy_from_slice(&rest[..n]);
                n
            };
            self.0 += n as u64;
            Ok(n)
        }
    }
    let mut reader = seqair::vcf::VcfStream::new(Long(0)).reader::<Locus>().open().unwrap();
    let text_error = |e: VcfReadError| match e {
        VcfReadError::Record { source, .. } => source,
        e => panic!("{e:?}"),
    };
    let first = text_error(reader.next_record().unwrap_err());
    assert!(matches!(first, RecordError::LineTooLong { .. }), "{first:?}");
    assert_eq!(text_error(reader.next_record().unwrap_err()), first);
}

// r[verify vcf_fields.only_sample]
// r[verify vcf_fields.lazy+4]
/// One sample's FORMAT data is copied whole, unchecked: malformed data (a
/// truncated column, bytes after the last field) fails on access with the
/// same typed error the multi-sample index gives, and the next record reads.
#[test]
fn a_malformed_single_sample_fails_on_access() {
    let trailing = Rec { indiv: vec![0x11, 0x02, 0x21, 0x02, 0x04, 0x99], ..Rec::default() };
    let cut = Rec { indiv: vec![0x11, 0x02, 0x21, 0x02], ..Rec::default() };
    let file = VcfFile::from_bytes(bcf(&[trailing, cut, Rec::default()]), None).unwrap();
    let gt: FormatGt = file.header().format_key("GT").unwrap();
    let mut only = file.reader::<OnlySample>().open().unwrap();
    let mut all = file.reader::<Samples>().open().unwrap();
    let sample = file.header().sample_indices().next().unwrap();
    let expected = [
        RecordError::TrailingBytes { part: RecordPart::Format, len: 1 },
        RecordError::Value {
            part: RecordPart::Format,
            source: BcfValueError::Truncated { needed: 3, available: 2 },
        },
    ];
    for e in expected {
        let one = only.next_record().unwrap().unwrap();
        assert_eq!(record_error(one.get(&gt).unwrap_err()).1, e);
        assert_eq!(record_error(one.has(&gt).unwrap_err()).1, e);
        let every = all.next_record().unwrap().unwrap();
        assert_eq!(record_error(every.sample(sample).get(&gt).unwrap_err()).1, e);
    }
    let one = only.next_record().unwrap().unwrap();
    assert!(one.get(&gt).unwrap().is_some());
}

// r[verify vcf_fields.sample+2]
// r[verify vcf_fields.lazy+4]
/// One sample of a multi-sample BCF is copied out of each column: malformed
/// data fails every access with the error the whole section's index gives,
/// not the read, and the next record reads.
#[test]
fn a_malformed_sample_of_many_fails_on_access() {
    let header = HEADER.replace("\tS1\n", "\tS1\tS2\n");
    let gt = |bytes: &[u8]| Rec { indiv: bytes.to_vec(), n_sample: 2, ..Rec::default() };
    let records = [
        gt(&[0x11, 0x02, 0x21, 0x02, 0x04, 0x02, 0x02, 0x99]),
        gt(&[0x11, 0x02, 0x21, 0x02, 0x04]),
        gt(&[0x11, 0x02, 0x21, 0x02, 0x04, 0x04, 0x04]),
    ];
    let file = VcfFile::from_bytes(bcf_with(&header, &records), None).unwrap();
    let key: FormatGt = file.header().format_key("GT").unwrap();
    let s2 = file.header().sample("S2").unwrap();
    let mut one = file.reader::<Sample>().sample("S2").open().unwrap();
    let mut all = file.reader::<Samples>().open().unwrap();
    for _ in 0..2 {
        let sample = one.next_record().unwrap().unwrap();
        let every = all.next_record().unwrap().unwrap();
        let expected = record_error(every.sample(s2).get(&key).unwrap_err()).1;
        assert!(
            matches!(expected, RecordError::TrailingBytes { .. } | RecordError::Value { .. }),
            "{expected:?}"
        );
        assert_eq!(record_error(sample.get(&key).unwrap_err()).1, expected);
        assert_eq!(record_error(sample.get(&key).unwrap_err()).1, expected, "every access");
        assert_eq!(record_error(sample.has(&key).unwrap_err()).1, expected);
        assert!(matches!(sample.iter().next(), Some(Err(_))));
    }
    let sample = one.next_record().unwrap().unwrap();
    assert!(sample.get(&key).unwrap().is_some());
    assert!(one.next_record().unwrap().is_none());
}

// r[verify vcf_fields.lazy+4]
// r[verify vcf_reader.field_table]
/// An INFO key the dictionary does not declare fails the lazily built index,
/// and every access returns the same error; the record's other parts read.
#[test]
fn a_failed_index_is_returned_on_every_access() {
    let bad = Rec { info: vec![0x11, 0x07, 0x11, 5], ..Rec::default() };
    let file = VcfFile::from_bytes(bcf(&[Rec::default(), bad, Rec::default()]), None).unwrap();
    let dp: InfoInt = file.header().info_key("DP").unwrap();
    let mut reader = file.reader::<(Locus, Info)>().open().unwrap();
    let mut values = Vec::new();
    while let Some((locus, info)) = reader.next_record().unwrap() {
        assert!(matches!(locus.start(), VcfPos::At(_)));
        let first = info.get(&dp).map_err(|e| e.to_string());
        let second = info.get(&dp).map_err(|e| e.to_string());
        assert_eq!(first, second);
        if let Err(e) = &first {
            // r[verify vcf_values.dynamic+2]
            let items: Vec<_> = info.iter().collect();
            assert!(matches!(items[..], [Err(_)]), "the index error is the one item");
            assert!(info.has(&dp).is_err());
            assert!(e.starts_with("malformed record 2 at chr1:11"), "{e}");
        }
        values.push(first.ok().flatten());
    }
    // The index buffers survive the failure: the next record reads.
    assert_eq!(values, [Some(5), None, Some(5)]);
}

// ── The lazy index under interleaved accesses ──────────────────────────

/// Three samples in VCF text: INFO flags, `key=`, a repeated and an
/// undeclared key; FORMAT with dropped trailing fields, a value that does
/// not decode, and a sample with more fields than FORMAT has keys (which
/// fails the index).
const TEXT_PANEL: &str = "##fileformat=VCFv4.3\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"a\">\n\
##INFO=<ID=DB,Number=0,Type=Flag,Description=\"b\">\n\
##INFO=<ID=NOTE,Number=1,Type=String,Description=\"n\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"q\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"a\">\n\
##contig=<ID=chr1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n\
chr1\t10\t.\tA\tC\t.\tPASS\tDP=5;DB;AF=0.5;NOTE=x\tGT:GQ:AD\t0/1:30:3,2\t1/1:.:0,5\t0/0\n\
chr1\t20\t.\tA\tC,G\t.\tPASS\tAF=0.1,0.2;DP=7;DP=8;XX=1;NOTE=\tGT:AD\t0|1:1,2,3\t.\t1/2:.\n\
chr1\t30\t.\tA\tC\t.\tPASS\t.\tGT:GQ\t0/1:x\t1/1:5\t0/0:6\n\
chr1\t40\t.\tA\tC\t.\tPASS\tDP=abc;DB\tGT\t0/1\t1/1:5\t0/0\n";

/// One access to a record's INFO or FORMAT data.
#[derive(Debug, Clone)]
enum Access {
    Get(usize),
    Has(usize),
    /// Start iterating INFO (held open across later accesses), or step it.
    Iter,
    Format(usize),
    SampleHas(usize, usize),
    SampleIter(usize),
    /// Read a column (held open across later accesses), or step it.
    Column(usize),
    Clone,
}

#[hegel::composite]
fn accesses(tc: &TestCase, n_info: usize, n_format: usize, n_samples: usize) -> Vec<Access> {
    let info = || gs::integers::<usize>().max_value(n_info - 1);
    let format = || gs::integers::<usize>().max_value(n_format - 1);
    let sample = || gs::integers::<usize>().max_value(n_samples - 1);
    let access = hegel::one_of!(
        info().map(Access::Get),
        info().map(Access::Has),
        gs::just(Access::Iter),
        format().map(Access::Format),
        gs::tuples2(sample(), format()).map(|(s, k)| Access::SampleHas(s, k)),
        sample().map(Access::SampleIter),
        format().map(Access::Column),
        gs::just(Access::Clone),
    );
    tc.draw(gs::vecs(access).max_size(12))
}

/// What one access reads, as text: values, errors and absence alike.
fn read(info: &Info, samples: &Samples, keys: &Keys, access: &Access) -> String {
    let s = |i: usize| samples.sample(keys.samples[i]);
    match *access {
        Access::Get(k) => format!("{:?}", info.get(&keys.info[k])),
        Access::Has(k) => format!("{:?}", info.has(&keys.info[k])),
        Access::Iter => format!("{:?}", info.iter().collect::<Vec<_>>()),
        Access::Format(k) => format!("{:?}", s(0).get(&keys.format[k])),
        Access::SampleHas(i, k) => format!("{:?}", s(i).has(&keys.format[k])),
        Access::SampleIter(i) => format!("{:?}", s(i).iter().collect::<Vec<_>>()),
        Access::Column(k) => {
            format!("{:?}", samples.column(&keys.format[k]).map(|c| c.map(Vec::from_iter)))
        }
        Access::Clone => String::new(),
    }
}

struct Keys {
    info: Vec<InfoAny>,
    format: Vec<FormatAny>,
    samples: Vec<seqair::vcf::SampleIdx>,
}

// r[verify vcf_fields.lazy+4]
/// Accesses interleaved on one record — lookups and `has` in any order, an
/// INFO iterator and a column held open across them, clones at any point —
/// never panic, and each reads what the same access reads on a fresh copy
/// of the record; a clone reads what its source reads.
#[hegel::test]
fn interleaved_accesses_read_what_a_fresh_record_reads(tc: TestCase) {
    let text = VcfFile::from_bytes(TEXT_PANEL.as_bytes().to_vec(), None).unwrap();
    let bcf = VcfFile::open(fixture("panel.bcf")).unwrap();
    let file = if tc.draw(gs::booleans()) { text } else { bcf };
    let header = file.header();
    let mut info: Vec<InfoAny> = header.infos().keys().map(|n| file.info_key(n).unwrap()).collect();
    info.extend(file.lenient().info_key::<seqair::vcf::Any>("XX"));
    let keys = Keys {
        info,
        format: header.formats().keys().map(|n| file.format_key(n).unwrap()).collect(),
        samples: header.sample_indices().collect(),
    };
    let records: Vec<(Info, Samples)> =
        file.reader::<(Info, Samples)>().open().unwrap().into_iter().map(Result::unwrap).collect();
    let (info, samples) = &records[tc.draw(gs::integers::<usize>().max_value(records.len() - 1))];
    let (pristine_info, pristine_samples) = (info.clone(), samples.clone());
    let (n_info, n_format, n_samples) = (keys.info.len(), keys.format.len(), keys.samples.len());
    let mut held_iter = None;
    let mut held_column = None;
    for access in tc.draw(accesses(n_info, n_format, n_samples)) {
        let (fresh_info, fresh_samples) = (pristine_info.clone(), pristine_samples.clone());
        match access {
            // Each step reads the index again while other accesses interleave.
            Access::Iter => match held_iter.as_mut() {
                None => held_iter = Some(info.iter()),
                Some(iter) => drop(iter.next()),
            },
            Access::Column(k) => match held_column.as_mut() {
                None => held_column = samples.column(&keys.format[k]).ok().flatten(),
                Some(column) => drop(column.next()),
            },
            Access::Clone => {
                let (copy_info, copy_samples) = (info.clone(), samples.clone());
                for k in 0..n_info {
                    let get = Access::Get(k);
                    let expected = read(&fresh_info, &fresh_samples, &keys, &get);
                    assert_eq!(read(&copy_info, &copy_samples, &keys, &get), expected);
                }
                let iter = Access::Iter;
                let expected = read(&fresh_info, &fresh_samples, &keys, &iter);
                assert_eq!(read(&copy_info, &copy_samples, &keys, &iter), expected);
            }
            _ => {}
        }
        let expected = read(&fresh_info, &fresh_samples, &keys, &access);
        assert_eq!(read(info, samples, &keys, &access), expected, "{access:?}");
    }
}

// r[verify vcf_reader.file.header_only]
/// A BCF whose header fills its first block opens with nothing readable
/// after that block: the open inflates the header's block and no other.
#[test]
fn opening_reads_no_block_past_the_header() {
    let bytes = std::fs::read(fixture("rastair.bcf")).unwrap();
    let first = usize::from(u16::from_le_bytes([bytes[16], bytes[17]])) + 1;
    let mut cut = bytes[..first].to_vec();
    cut.extend_from_slice(b"not a BGZF block");
    let file = VcfFile::from_bytes(cut, None).unwrap();
    let mut reader = file.reader::<Locus>().open().unwrap();
    assert!(reader.next_record().is_err(), "the records were cut off");
    // The whole file reads from the same start.
    let whole = VcfFile::from_bytes(bytes, None).unwrap();
    assert!(whole.reader::<Locus>().open().unwrap().next_record().unwrap().is_some());
}

// r[verify vcf_reader.formats]
#[test]
fn formats_are_detected_and_refused_by_content() {
    let gzip = {
        use std::io::Write;
        let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        e.write_all(HEADER.as_bytes()).unwrap();
        e.finish().unwrap()
    };
    assert!(matches!(
        VcfFile::from_bytes(gzip, None),
        Err(VcfReadError::Open { path: None, source: OpenError::NotBgzf })
    ));
    let mut bcf21 = bcf(&[]);
    let mut w = BgzfWriter::new(Vec::new());
    w.write_all(b"BCF\x02\x01").unwrap();
    bcf21.clone_from(&w.finish().unwrap());
    assert!(matches!(
        VcfFile::from_bytes(bcf21, None),
        Err(VcfReadError::Open { path: None, source: OpenError::Bcf21 })
    ));
    let plain = VcfFile::from_bytes(HEADER.as_bytes().to_vec(), None).unwrap();
    assert_eq!(plain.format(), seqair::vcf::reader::VcfFormat::Vcf);
    let mut w = BgzfWriter::new(Vec::new());
    w.write_all(HEADER.as_bytes()).unwrap();
    let gz = VcfFile::from_bytes(w.finish().unwrap(), None).unwrap();
    assert_eq!(gz.format(), seqair::vcf::reader::VcfFormat::VcfGz);
    assert_eq!(
        VcfFile::from_bytes(bcf(&[]), None).unwrap().format(),
        seqair::vcf::reader::VcfFormat::Bcf
    );
    assert!(matches!(
        VcfFile::from_bytes(b"hello".to_vec(), None),
        Err(VcfReadError::Open { path: None, source: OpenError::UnknownFormat })
    ));
    // A file with no records reads none.
    let empty = VcfFile::from_bytes(bcf(&[]), None).unwrap();
    assert!(empty.reader::<FullRecord>().open().unwrap().next_record().unwrap().is_none());
}

// r[verify vcf_reader.file+2]
#[test]
fn a_missing_index_is_fine_and_a_broken_one_fails_open() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("x.bcf");
    std::fs::write(&path, bcf(&[Rec::default()])).unwrap();
    let file = VcfFile::open(&path).unwrap();
    assert!(!file.has_index());
    assert_eq!(file.reader::<Locus>().open().unwrap().into_iter().count(), 1);
    std::fs::write(dir.path().join("x.bcf.csi"), b"not an index").unwrap();
    assert!(matches!(
        VcfFile::open(&path),
        Err(VcfReadError::Open { path: Some(p), source: OpenError::UnknownIndexFormat })
            if p == path
    ));
}

// r[verify vcf_reader.errors.open]
/// Every open error names the file, once in its whole chain: a header that
/// does not parse, a format refused, a file that is not there, an index
/// that does not parse; bytes in memory are "the input".
#[test]
fn open_errors_name_the_file_once() {
    let dir = tempfile::tempdir().unwrap();
    let chain = |e: VcfReadError| format!("{:#}", anyhow::Error::from(e));
    let once = |chain: &str, path: &std::path::Path| {
        let shown = path.display().to_string();
        assert_eq!(chain.matches(&shown).count(), 1, "{chain}");
    };

    let bad = dir.path().join("bad.vcf");
    let text = "##fileformat=VCFv4.3\n\
                ##INFO=<ID=X,Number=Q,Type=Integer,Description=\"x\">\n\
                #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    std::fs::write(&bad, text).unwrap();
    let e = VcfFile::open(&bad).unwrap_err();
    assert!(matches!(&e, VcfReadError::Open { source: OpenError::Header(_), .. }), "{e:?}");
    assert_eq!(e.to_string(), format!("could not open {}", bad.display()));
    let c = chain(e);
    assert!(c.contains("unknown Number=Q"), "{c}");
    once(&c, &bad);

    let hello = dir.path().join("hello.vcf");
    std::fs::write(&hello, "hello").unwrap();
    let c = chain(VcfFile::open(&hello).unwrap_err());
    assert!(c.ends_with("the input is neither VCF nor BCF"), "{c}");
    once(&c, &hello);

    let missing = dir.path().join("missing.vcf");
    let e = VcfFile::open(&missing).unwrap_err();
    assert!(matches!(&e, VcfReadError::Open { source: OpenError::Io(_), .. }), "{e:?}");
    once(&chain(e), &missing);

    let bcf_path = dir.path().join("x.bcf");
    std::fs::write(&bcf_path, bcf(&[Rec::default()])).unwrap();
    std::fs::write(dir.path().join("x.bcf.csi"), b"not an index").unwrap();
    let c = chain(VcfFile::open(&bcf_path).unwrap_err());
    assert!(c.ends_with("the index is neither TBI nor CSI"), "{c}");
    once(&c, &bcf_path);

    let e = VcfFile::from_bytes(b"hello".to_vec(), None).unwrap_err();
    assert_eq!(e.to_string(), "could not open the input");
}

// ── Crafted VCF text ────────────────────────────────────────────────────

const TEXT_HEADER: &str = "##fileformat=VCFv4.3\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##FILTER=<ID=q10,Description=\"q\">\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"a\">\n\
##INFO=<ID=DB,Number=0,Type=Flag,Description=\"b\">\n\
##INFO=<ID=S,Number=1,Type=String,Description=\"s\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"a\">\n\
##contig=<ID=chr1,length=1000>\n\
##contig=<ID=chr2,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n";

fn text_file(lines: &[&str]) -> VcfFile {
    let mut text = TEXT_HEADER.to_owned();
    for l in lines {
        text.push_str(l);
        text.push('\n');
    }
    VcfFile::from_bytes(text.into_bytes(), None).unwrap()
}

fn text_records(lines: &[&str]) -> Vec<FullRecord> {
    let file = text_file(lines);
    file.reader::<FullRecord>().open().unwrap().into_iter().collect::<Result<_, _>>().unwrap()
}

fn text_error(lines: &[&str]) -> (u64, Option<seqair::io::VirtualOffset>, RecordError) {
    let file = text_file(lines);
    let mut reader = file.reader::<FullRecord>().open().unwrap();
    loop {
        match reader.next_record() {
            Ok(Some(_)) => {}
            Ok(None) => panic!("no error"),
            Err(VcfReadError::Record { at, source }) => {
                return (at.number.unwrap(), at.offset, source);
            }
            Err(e) => panic!("expected a record error, got {e:?}"),
        }
    }
}

// r[verify vcf_reader.formats]
// r[verify vcf_fields.site+2]
// r[verify vcf_values.access+2]
/// The fixed columns, INFO and FORMAT of hand-written lines: IDs, symbolic
/// and missing ALTs, QUAL, every FILTER state, flags, strings, vectors.
#[test]
fn text_records_read_their_columns() {
    let recs = text_records(&[
        "chr1\t10\trs1;rs2\tA\tC,<DEL>\t29.5\tPASS\tDP=7;AF=0.5,0.25;DB;S=x%3By\tGT:DP:AD\t0/1:3:1,2\t1|1:4:0,4",
        "chr1\t0\t.\tN\t.\t.\t.\t.\tGT\t.\t./.",
        "chr2\t5\t.\tAC\tA\t-1\tq10\tDP=.\tGT:DP\t0/1:.\t1/1:9",
    ]);
    let at = |p| VcfPos::At(seqair_types::Pos1::new(p).unwrap());
    let [a, b, c] = &recs[..] else { panic!("three records") };
    // r[verify vcf_fields.locus.bases0]
    let pos0 = |p| seqair_types::Pos0::new(p).unwrap();
    assert_eq!(a.locus().bases0(), Some(pos0(9)..pos0(10)));
    assert_eq!(b.locus().bases0(), None, "the telomere only");
    assert_eq!(c.locus().bases0(), Some(pos0(4)..pos0(6)));
    assert_eq!(
        (a.locus().contig().name().as_str(), a.locus().start(), a.locus().end()),
        ("chr1", at(10), at(10))
    );
    assert_eq!(a.site().ids().collect::<Vec<_>>(), ["rs1", "rs2"]);
    assert_eq!(a.site().alleles().reference(), b"A");
    assert_eq!(
        a.site().alleles().alts().collect::<Vec<_>>(),
        [Allele::Bases(b"C"), Allele::Symbolic("<DEL>")]
    );
    assert_eq!(a.site().qual(), Some(29.5));
    assert_eq!(a.site().filters().status(), FilterStatus::Pass);
    assert_eq!(b.locus().start(), VcfPos::Telomere);
    assert_eq!(b.site().alleles().alts().count(), 0);
    assert_eq!((b.site().qual(), b.site().filters().status()), (None, FilterStatus::Missing));
    assert_eq!(c.locus().contig().name(), "chr2");
    assert_eq!((c.locus().start(), c.locus().end()), (at(5), at(6)));
    assert_eq!(c.site().qual(), Some(-1.0));
    let FilterStatus::Failed(failed) = c.site().filters().status() else { panic!("failed") };
    assert_eq!(failed.map(|f| f.name().to_owned()).collect::<Vec<_>>(), ["q10"]);

    let file = text_file(&[
        "chr1\t10\trs1;rs2\tA\tC,<DEL>\t29.5\tPASS\tDP=7;AF=0.5,0.25;DB;S=x%3By\tGT:DP:AD\t0/1:3:1,2\t1|1:4:0,4",
        "chr2\t5\t.\tAC\tA\t-1\tq10\tDP=.\tGT:DP\t0/1:.\t1/1:9",
    ]);
    let h = file.header();
    let dp: InfoInt = h.info_key("DP").unwrap();
    let af: seqair::vcf::InfoFloats = h.info_key("AF").unwrap();
    let db: InfoFlag = h.info_key("DB").unwrap();
    let st: seqair::vcf::InfoString = h.info_key("S").unwrap();
    let gt: FormatGt = h.format_key("GT").unwrap();
    let fdp: seqair::vcf::FormatInt = h.format_key("DP").unwrap();
    let ad: seqair::vcf::FormatInts = h.format_key("AD").unwrap();
    let (s1, s2) = (h.sample("S1").unwrap(), h.sample("S2").unwrap());
    let recs: Vec<FullRecord> =
        file.reader::<FullRecord>().open().unwrap().into_iter().collect::<Result<_, _>>().unwrap();
    let info = recs[0].info();
    assert_eq!(info.get(&dp).unwrap(), Some(7));
    assert_eq!(info.get(&af).unwrap().map(Iterator::collect::<Vec<_>>), Some(vec![0.5, 0.25]));
    assert!(info.get(&db).unwrap());
    assert_eq!(info.get(&st).unwrap(), Some("x%3By"), "percent-encoding is left alone");
    let samples = recs[0].samples();
    let gts: Vec<_> =
        samples.column(&gt).unwrap().unwrap().map(|g| g.unwrap().to_genotype()).collect();
    assert_eq!(
        gts,
        [seqair::vcf::Genotype::unphased(0, 1), seqair::vcf::Genotype::phased_diploid(1, 1)]
    );
    assert_eq!(samples.sample(s2).get(&fdp).unwrap(), Some(4));
    assert_eq!(
        samples.sample(s1).get(&ad).unwrap().map(Iterator::collect::<Vec<_>>),
        Some(vec![1, 2])
    );
    assert_eq!(recs[1].info().get(&dp).unwrap(), None);
    assert!(!recs[1].info().get(&db).unwrap());
    assert_eq!(recs[1].samples().sample(s1).get(&fdp).unwrap(), None);
    assert_eq!(recs[1].samples().sample(s2).get(&ad).unwrap().map(|v| v.count()), None);
}

// r[verify vcf_fields.samples+2]
// r[verify vcf_values.access+2]
// r[verify vcf_values.duplicates+2]
// r[verify vcf_reader.field_table]
// r[verify vcf_reader.undeclared.filter+2]
/// A sample may drop trailing FORMAT fields (they read as missing); a
/// repeated key reads its first value; an undeclared INFO or FORMAT key
/// reads as a String; a flag given a value is set; an integer key without a
/// value is absent; an undeclared FILTER fails the record without being
/// listed.
#[test]
fn text_edge_cases_read_as_htslib_reads_them() {
    let file = text_file(&[
        "chr1\t10\t.\tA\tC\t.\tPASS\tDP=1;XX=3;DP=2;DB=yes\tGT:DP:ZZ:AD\t0/1\t1/1:5:q:2,3",
        "chr1\t11\t.\tA\tC\t.\tnope;q10\tDP\tGT:GT\t0|1:1|1\t./.:0/0",
        "chr1\t12\t.\tA\tC\t.\tnope\t.\t.\t.\t.",
    ]);
    let h = file.header();
    let dp: InfoInt = h.info_key("DP").unwrap();
    let db: InfoFlag = h.info_key("DB").unwrap();
    let gt: FormatGt = h.format_key("GT").unwrap();
    let fdp: seqair::vcf::FormatInt = h.format_key("DP").unwrap();
    let ad: seqair::vcf::FormatInts = h.format_key("AD").unwrap();
    let (s1, s2) = (h.sample("S1").unwrap(), h.sample("S2").unwrap());
    let recs: Vec<FullRecord> =
        file.reader::<FullRecord>().open().unwrap().into_iter().collect::<Result<_, _>>().unwrap();
    let info = recs[0].info();
    assert_eq!(info.get(&dp).unwrap(), Some(1), "the first DP wins");
    assert!(info.get(&db).unwrap(), "a flag with a value is set");
    let fields: Vec<_> = info.iter().map(|f| f.unwrap()).collect();
    let names: Vec<&str> = fields.iter().map(|f| f.0).collect();
    assert_eq!(names, ["DP", "XX", "DB"], "the second DP is a duplicate");
    assert!(
        matches!(fields[1], (_, def, DynValue::Str(Some("3"))) if def.number == seqair::vcf::Number::Unknown)
    );
    let one = recs[0].samples().sample(s1);
    assert!(one.get(&gt).unwrap().is_some());
    assert_eq!(one.get(&fdp).unwrap(), None, "dropped trailing fields are missing");
    assert_eq!(one.get(&ad).unwrap().map(|v| v.count()), None);
    let fields: Vec<&str> = one.iter().map(|f| f.unwrap().0).collect();
    assert_eq!(fields, ["GT", "DP", "ZZ", "AD"]);
    let two = recs[0].samples().sample(s2);
    let zz = two.iter().map(|f| f.unwrap()).find(|f| f.0 == "ZZ").unwrap().2;
    assert!(matches!(zz, DynValue::Str(Some("q"))), "{zz:?}");
    assert_eq!(two.get(&ad).unwrap().map(Iterator::collect::<Vec<_>>), Some(vec![2, 3]));

    assert_eq!(recs[1].info().get(&dp).unwrap(), None, "a key without a value holds nothing");
    let FilterStatus::Failed(failed) = recs[1].site().filters().status() else { panic!() };
    assert_eq!(failed.map(|f| f.name().to_owned()).collect::<Vec<_>>(), ["q10"]);
    let g = recs[1].samples().sample(s1).get(&gt).unwrap().unwrap().to_genotype();
    assert_eq!(g, seqair::vcf::Genotype::phased_diploid(0, 1), "the first GT wins");
    let FilterStatus::Failed(failed) = recs[2].site().filters().status() else {
        panic!("an undeclared FILTER fails the record")
    };
    assert_eq!(failed.count(), 0);
    assert_eq!(recs[2].samples().sample(s1).iter().count(), 0, "FORMAT . has no fields");
}

/// Counts WARN events; nothing else about tracing matters here.
struct WarnCounter(std::sync::Arc<std::sync::atomic::AtomicUsize>);

impl tracing::Subscriber for WarnCounter {
    fn enabled(&self, _: &tracing::Metadata<'_>) -> bool {
        true
    }
    fn new_span(&self, _: &tracing::span::Attributes<'_>) -> tracing::span::Id {
        tracing::span::Id::from_u64(1)
    }
    fn record(&self, _: &tracing::span::Id, _: &tracing::span::Record<'_>) {}
    fn record_follows_from(&self, _: &tracing::span::Id, _: &tracing::span::Id) {}
    fn event(&self, event: &tracing::Event<'_>) {
        if *event.metadata().level() == tracing::Level::WARN {
            self.0.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        }
    }
    fn enter(&self, _: &tracing::span::Id) {}
    fn exit(&self, _: &tracing::span::Id) {}
}

/// The `warn!`s reading every record of `file` as `F` logs, after `visit`
/// sees each record.
fn warnings<F: Fields>(file: &VcfFile, mut visit: impl FnMut(&F)) -> usize {
    let warns = std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0));
    tracing::subscriber::with_default(WarnCounter(warns.clone()), || {
        let mut reader = file.reader::<F>().open().unwrap();
        while let Some(record) = reader.next_record().unwrap() {
            visit(record);
        }
    });
    warns.load(std::sync::atomic::Ordering::Relaxed)
}

// r[verify vcf_values.text_numbers+2]
// r[verify vcf_values.equivalence]
/// A text integer in BCF's reserved band reads as missing, as the BCF
/// bcftools writes from the same text does, with one warning per file; a
/// copy of the record writes it missing instead of failing.
#[test]
fn reserved_band_integers_in_text_read_as_bcftools_stores_them() {
    let lines = [
        "chr1\t10\t.\tA\tC\t.\tPASS\tDP=-2147483645\tGT:DP:AD\t0/1:-2147483648:1,-2147483641\t1/1:5:2,3",
        "chr1\t11\t.\tA\tC\t.\tPASS\tDP=-2147483640\tGT:DP:AD\t0/1:-2147483641:.,4\t1/1:5:2,3",
    ];
    let dir = tempfile::tempdir().unwrap();
    let (vcf, bcf) = (dir.path().join("in.vcf"), dir.path().join("in.bcf"));
    std::fs::write(&vcf, TEXT_HEADER.to_owned() + &lines.join("\n") + "\n").unwrap();
    let (vcf_str, bcf_str) = (vcf.to_str().unwrap(), bcf.to_str().unwrap());
    if bcftools(&["view", "-Ob", "-o", bcf_str, vcf_str]).is_none() {
        return;
    }
    let shown = |path: &PathBuf| -> Vec<String> {
        let file = VcfFile::open(path).unwrap();
        let mut reader = file.reader::<FullRecord>().open().unwrap();
        let mut out = Vec::new();
        while let Some(r) = reader.next_record().unwrap() {
            let samples = file.header().sample_indices().map(|i| r.samples().sample(i));
            out.push(format!("{:?} {:?}", r.info(), samples.collect::<Vec<_>>()));
        }
        out
    };
    let text = shown(&vcf);
    assert_eq!(text, shown(&bcf));
    assert!(text[0].contains("\"DP\": Int([])"), "{}", text[0]);
    assert!(text[1].contains("\"DP\": Int([Some(-2147483640)])"), "{}", text[1]);

    let file = VcfFile::open(&vcf).unwrap();
    let dp: InfoInt = file.header().info_key("DP").unwrap();
    let ad: seqair::vcf::FormatIntOpts = file.header().format_key("AD").unwrap();
    let s1 = file.header().sample("S1").unwrap();
    assert_eq!(
        warnings::<FullRecord>(&file, |r| {
            assert!(matches!(r.info().get(&dp).unwrap(), None | Some(-2147483640)));
            r.samples().sample(s1).get(&ad).unwrap().map(|v| v.count());
        }),
        1,
        "once per file"
    );
    let mut out = Vec::new();
    let mut writer = seqair::vcf::Writer::new(&mut out, seqair::vcf::OutputFormat::Bcf)
        .write_header(file.header())
        .unwrap();
    let mut reader = file.reader::<FullRecord>().open().unwrap();
    while let Some(record) = reader.next_record().unwrap() {
        writer.copy_record(record).unwrap().emit().unwrap();
    }
    writer.finish().unwrap();
    std::fs::write(&bcf, out).unwrap();
    assert_eq!(shown(&bcf), text, "the copy");
}

// r[verify vcf_reader.undeclared+5]
// r[verify vcf_reader.undeclared.filter+2]
/// Undeclared FILTER names warn once per file wherever `Site` reads them;
/// undeclared keys warn once more, but only where an index meets one.
#[test]
fn undeclared_filters_and_keys_warn_once_each_where_met() {
    let lines = [
        "chr1\t10\t.\tA\tC\t.\tnope\tDP=1;XX=3\tGT:ZZ\t0/1:a\t1/1:b",
        "chr1\t11\t.\tA\tC\t.\tnah;q10\tDP=2;YY=2\tGT:WW\t0/1:a\t1/1:b",
    ];
    let file = || text_file(&lines);
    assert_eq!(warnings::<Locus>(&file(), |_| {}), 0, "nothing reads FILTER or INFO");
    assert_eq!(warnings::<(Locus, Site)>(&file(), |_| {}), 1, "the FILTER names");
    let fresh = file();
    let dp: InfoInt = fresh.header().info_key("DP").unwrap();
    let early = warnings::<Info>(&fresh, |info| {
        info.get(&dp).unwrap();
    });
    assert_eq!(early, 0, "a lookup of DP stops before the undeclared keys");
    let all = warnings::<FullRecord>(&file(), |r| {
        r.info().iter().for_each(|f| assert!(f.is_ok()));
    });
    assert_eq!(all, 2, "one for the FILTER names, one for the keys");
}

// r[verify vcf_reader.contigs]
/// A CHROM the header does not declare is added after the header's
/// contigs, once, and keeps its id for every record of the file.
#[test]
fn undeclared_text_contigs_get_ids_after_the_headers() {
    let file = text_file(&[
        "chrUn\t1\t.\tA\t.\t.\t.\t.",
        "chrUn\t2\t.\tA\t.\t.\t.\t.",
        "chr1\t1\t.\tA\t.\t.\t.\t.",
        "chrM\t1\t.\tA\t.\t.\t.\t.",
        "chrUn\t3\t.\tA\t.\t.\t.\t.",
    ]);
    assert_eq!(file.contig("chrUn"), None, "no record named it yet");
    let recs: Vec<FullRecord> =
        file.reader::<FullRecord>().open().unwrap().into_iter().collect::<Result<_, _>>().unwrap();
    let ids: Vec<(&str, u32)> = recs
        .iter()
        .map(|r| (r.locus().contig().name().as_str(), r.locus().contig().tid()))
        .collect();
    assert_eq!(ids, [("chrUn", 2), ("chrUn", 2), ("chr1", 0), ("chrM", 3), ("chrUn", 2)]);
    // r[verify vcf_reader.contig_lookup]
    assert_eq!(file.contig("chrUn").map(|c| c.tid()), Some(2));
    assert_eq!(file.contig("chr1"), file.header().contig("chr1").ok());
    assert_eq!(file.contig("chrX"), None);
}

// r[verify vcf_record.gt_first_phase+3]
// r[verify vcf_values.gt+2]
/// Before VCF 4.4 the first allele's phase is unphased; from 4.4 a leading
/// indicator gives it, and without one it is inferred as htslib does.
#[test]
fn text_genotype_phase_follows_the_header_version() {
    let read = |version: &str, gt: &str| {
        let text = format!(
            "##fileformat=VCFv{version}\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
             ##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\
             chr1\t1\t.\tA\tC\t.\t.\t.\tGT\t{gt}\n"
        );
        let file = VcfFile::from_bytes(text.into_bytes(), None).unwrap();
        let key: FormatGt = file.header().format_key("GT").unwrap();
        let rec = file.reader::<OnlySample>().open().unwrap().read().unwrap().unwrap();
        rec.get(&key)
            .unwrap()
            .map(|g| g.alleles().map(|a| (a.allele, a.phased)).collect::<Vec<_>>())
    };
    assert_eq!(read("4.3", "0|1"), Some(vec![(Some(0), false), (Some(1), true)]));
    assert_eq!(read("4.4", "0|1"), Some(vec![(Some(0), true), (Some(1), true)]));
    assert_eq!(read("4.4", "/0|1"), Some(vec![(Some(0), false), (Some(1), true)]));
    assert_eq!(read("4.4", "1"), Some(vec![(Some(1), true)]));
    assert_eq!(read("4.5", ".|."), Some(vec![(None, true), (None, true)]));
    assert_eq!(read("4.5", "."), None);
}

// r[verify vcf_fields.sections]
// r[verify vcf_fields.sample+2]
/// `Sample` reads its column and no further: a malformed later column does
/// not fail it, while `Samples` reads (and refuses) them all.
#[test]
fn one_text_sample_scans_only_up_to_its_column() {
    let line = "chr1\t1\t.\tA\tC\t.\t.\t.\tGT\t0/1\t0/1:1:2:3:4:5";
    let file = text_file(&[line]);
    let gt: FormatGt = file.header().format_key("GT").unwrap();
    let mut one = file.reader::<Sample>().sample("S1").open().unwrap();
    assert!(one.next_record().unwrap().unwrap().get(&gt).unwrap().is_some());
    let mut all = file.reader::<Samples>().open().unwrap();
    let samples = all.next_record().unwrap().unwrap();
    let e = samples.column(&gt).unwrap_err();
    assert!(
        matches!(
            e,
            VcfReadError::Record { source: RecordError::SampleFields { sample: 1, .. }, .. }
        ),
        "{e:?}"
    );
    // Too few sample columns, for one sample and for all.
    let file = text_file(&["chr1\t1\t.\tA\tC\t.\t.\t.\tGT\t0/1"]);
    let e = file.reader::<Sample>().sample("S2").open().unwrap().read().unwrap_err();
    assert!(
        matches!(
            e,
            VcfReadError::Record { source: RecordError::SampleCount { found: 1, expected: 2 }, .. }
        ),
        "{e:?}"
    );
    let rec = file.reader::<Samples>().open().unwrap().read().unwrap().unwrap();
    assert!(rec.column(&gt).is_err());
}

// r[verify vcf_reader.errors+3]
// r[verify vcf_reader.limits+2]
/// Malformed lines are typed errors with the record's ordinal; in a BGZF
/// file, its virtual offset too. Blank lines are not records.
#[test]
fn malformed_text_lines_are_typed_errors_with_their_place() {
    let ok = "chr1\t1\t.\tA\tC\t.\t.\t.";
    let (ordinal, offset, e) = text_error(&[ok, "", ok, "chr1\t1\t.\tA"]);
    assert_eq!((ordinal, offset, e), (3, None, RecordError::TooFewColumns { found: 4 }));
    let (_, _, e) = text_error(&["chr1\tx\t.\tA\tC\t.\t.\t."]);
    assert!(matches!(e, RecordError::Text { part: seqair::vcf::reader::RecordPart::Pos, .. }));
    let (_, _, e) = text_error(&["chr1\t1\t.\tA\tC\tq\t.\t."]);
    assert!(matches!(e, RecordError::Text { part: seqair::vcf::reader::RecordPart::Qual, .. }));
    let (_, _, e) = text_error(&["chr1\t4294967296\t.\tA\tC\t.\t.\t."]);
    assert!(matches!(e, RecordError::InvalidPosition { .. }), "{e:?}");
    // Values are read on access: a malformed one fails there, with its
    // place; so does a vector where the header declares one value.
    let file = text_file(&[ok, "chr1\t1\t.\tA\tC\t.\t.\tDP=x", "chr1\t1\t.\tA\tC\t.\t.\tDP=1,2"]);
    let dp: InfoInt = file.header().info_key("DP").unwrap();
    let recs: Vec<Info> =
        file.reader::<Info>().open().unwrap().into_iter().collect::<Result<_, _>>().unwrap();
    let e = recs[1].get(&dp).unwrap_err();
    assert!(
        matches!(
            e,
            VcfReadError::Record {
                ref at,
                source: RecordError::Field { source: FieldError::Text(_), .. },
            } if at.number == Some(2)
        ),
        "{e:?}"
    );
    // The place says where in the file and the genome; the message leaves
    // the cause to the error's source, so a reporter prints it once.
    let place = "malformed record 2 (line 15) at chr1:1";
    assert_eq!(e.to_string(), place);
    let cause = std::error::Error::source(&e).unwrap();
    assert_eq!(cause.to_string(), "field DP");
    let value = std::error::Error::source(cause).unwrap().to_string();
    assert_eq!(value, "\"x\" is not an integer (an optional sign and decimal digits)");
    assert!(!format!("{:#}", anyhow::Error::from(e)).contains("field DP: field DP"));
    // A file on disk is named.
    let mut tmp = tempfile::NamedTempFile::new().unwrap();
    std::io::Write::write_all(
        &mut tmp,
        format!("{TEXT_HEADER}chr2\t7\t.\tA\tC\tq\t.\t.\n").as_bytes(),
    )
    .unwrap();
    let file = VcfFile::open(tmp.path()).unwrap();
    let e = file.reader::<Site>().open().unwrap().next_record().unwrap_err();
    let expected = format!("malformed record 1 (line 14) at chr2:7 in {}", tmp.path().display());
    assert_eq!(e.to_string(), expected);
    let e = recs[2].get(&dp).unwrap_err();
    assert!(
        matches!(
            e,
            VcfReadError::Record {
                source: RecordError::Field { source: FieldError::NotScalar { len: 2 }, .. },
                ..
            }
        ),
        "{e:?}"
    );

    // The same in BGZF: every error names where its record starts.
    let mut text = TEXT_HEADER.to_owned();
    text.push_str(&[ok, ok, "chr1\t1", ""].join("\n"));
    let mut w = BgzfWriter::new(Vec::new());
    w.write_all(text.as_bytes()).unwrap();
    let file = VcfFile::from_bytes(w.finish().unwrap(), None).unwrap();
    let mut reader = file.reader::<Locus>().open().unwrap();
    assert!(reader.next_record().unwrap().is_some());
    assert!(reader.next_record().unwrap().is_some());
    let e = reader.next_record().unwrap_err();
    let VcfReadError::Record { at, .. } = e else { panic!("{e:?}") };
    // Too few columns, but CHROM and POS parse: the place names them.
    let locus = at.locus.as_ref().map(|(c, p)| format!("{c}:{p}"));
    assert_eq!((at.number, at.line, locus.as_deref()), (Some(3), Some(16), Some("chr1:1")));
    let offset = at.offset.unwrap();
    let header_len = TEXT_HEADER.len() + 2 * (ok.len() + 1);
    assert_eq!(usize::from(offset.within_block()), header_len);
}

/// CRLF line ends, a last line without a newline, and lines longer than a
/// BGZF block or a file buffer, which are copied together.
#[test]
fn text_lines_are_framed_whatever_their_length_and_ending() {
    let long_info = format!("S={}", "x".repeat(300_000));
    let lines = [
        "chr1\t1\t.\tA\tC\t.\t.\t.\tGT\t0/1\t1/1\r".to_owned(),
        format!("chr1\t2\t.\tA\tC\t.\t.\t{long_info}\tGT\t0/1\t1/1"),
        "chr1\t3\t.\tA\tC\t.\t.\t.\tGT\t0/1\t1/1".to_owned(),
    ];
    let mut text = TEXT_HEADER.to_owned();
    text.push_str(&lines.join("\n")); // no final newline
    let dir = tempfile::tempdir().unwrap();
    let plain = dir.path().join("x.vcf");
    std::fs::write(&plain, &text).unwrap();
    let mut w = BgzfWriter::new(Vec::new());
    w.write_all(text.as_bytes()).unwrap();
    let gz = dir.path().join("x.vcf.gz");
    std::fs::write(&gz, w.finish().unwrap()).unwrap();
    for path in [plain, gz] {
        let file = VcfFile::open(&path).unwrap();
        let s: seqair::vcf::InfoString = file.header().info_key("S").unwrap();
        let gt: FormatGt = file.header().format_key("GT").unwrap();
        let s2 = file.header().sample("S2").unwrap();
        let recs: Vec<FullRecord> = file
            .reader::<FullRecord>()
            .open()
            .unwrap()
            .into_iter()
            .collect::<Result<_, _>>()
            .unwrap();
        assert_eq!(recs.len(), 3, "{path:?}");
        assert_eq!(recs[1].info().get(&s).unwrap().map(str::len), Some(300_000));
        for r in &recs {
            let g = r.samples().sample(s2).get(&gt).unwrap().unwrap().to_genotype();
            assert_eq!(g, seqair::vcf::Genotype::unphased(1, 1), "{path:?}: the CR is not a value");
        }
    }
}

/// riker's pattern (the reading guide's "Positions only"): `Locus` alone on
/// GIAB's truth set, the deletions spanning more than their first base.
#[test]
fn truth_loci_span_their_deletions() {
    let file = VcfFile::open(fixture("truth.vcf.gz")).unwrap();
    let mut reader = file.reader::<Locus>().open().unwrap();
    let mut longest = 0;
    while let Some(locus) = reader.next_record().unwrap() {
        if let (VcfPos::At(start), VcfPos::At(end)) = (locus.start(), locus.end()) {
            longest = longest.max(end.as_u32() - start.as_u32() + 1);
        }
    }
    assert!(longest > 1);
}

// r[verify vcf_reader.debug]
/// A record's `Debug` shows what it holds, not its file's tables; a locus
/// displays as `CHROM:POS`.
#[test]
fn records_debug_as_what_they_hold() {
    for name in ["rastair.bcf", "panel.bcf", "truth.vcf.gz"] {
        let file = VcfFile::open(fixture(name)).unwrap();
        let record = file.reader::<FullRecord>().open().unwrap().read().unwrap().unwrap();
        let debug = format!("{record:?}");
        // About the record's own size: its fields, not the header's.
        let limit = 6 * record.raw_bytes().len() + 300;
        assert!(debug.len() < limit, "{name}: {} bytes: {debug}", debug.len());
        for leak in ["string_map", "HeaderId", "contigs", "FileShared", "infos"] {
            assert!(!debug.contains(leak), "{name}: {debug}");
        }
        assert!(debug.contains(&format!("locus: {}", record.locus())), "{debug}");
        assert!(format!("{file:?}").len() < 200 + fixture(name).as_os_str().len(), "{file:?}");
    }
    let recs = text_records(&[
        "chr1\t5\trs1;rs2\tAC\tA,<DEL>\t3\tq10\tDP=4;AF=0.5,.;DB\tGT:DP\t0|1:7\t./.",
        "chr2\t0\t.\tN\t.\t.\tPASS\t.",
        "chr2\t1\t.\tN\t.\t.\tPASS;q10\t.",
    ]);
    // r[verify vcf_fields.site+2]
    assert_eq!(format!("{:?}", recs[2].site().filters()), r#"Failed(["PASS", "q10"])"#);
    assert!(recs[2].site().filters().contains(&seqair::vcf::FilterId::PASS).unwrap());
    assert_eq!(recs[0].locus().to_string(), "chr1:5-6");
    assert_eq!(recs[1].locus().to_string(), "chr2:0");
    let site = format!("{:?}", recs[0].site());
    assert_eq!(
        site,
        r#"Site { ids: ["rs1", "rs2"], alleles: [Bases("AC"), Bases("A"), Symbolic("<DEL>")], qual: Some(3.0), filters: Failed(["q10"]) }"#
    );
    assert_eq!(
        format!("{:?}", recs[0].info()),
        "Info {\"DP\": Int([Some(4)]), \"AF\": Float([Some(0.5), None]), \"DB\": Flag}"
    );
    assert_eq!(format!("{:?}", recs[0].samples()), r#"Samples { len: 2, fields: ["GT", "DP"] }"#);
}

// r[verify vcf_values.access+2]
/// One sample's view tells a field the record does not have from a `.`
/// and from a trailing field the sample dropped, without reading the other
/// samples (`HiPhase` phases a variant whose GQ is absent, not one whose GQ
/// is `.`).
#[test]
fn one_sample_tells_absent_from_missing() {
    let file = text_file(&[
        "chr1\t1\t.\tA\tC\t.\t.\tDP=3\tGT\t0/1\t0/1",
        "chr1\t2\t.\tA\tC\t.\t.\tDP=.\tGT:DP\t0/1:.\t0/1:4",
        "chr1\t3\t.\tA\tC\t.\t.\t.\tGT:DP\t0/1\t0/1:4",
        "chr1\t4\t.\tA\tC\t.\t.\t.\tGT:DP\t0/1:7\t0/1:4",
    ]);
    let dp: seqair::vcf::FormatInt = file.format_key("DP").unwrap();
    let any: seqair::vcf::FormatAny = file.format_key("DP").unwrap();
    let info_dp: InfoInt = file.info_key("DP").unwrap();
    let mut reader = file.reader::<(Info, Sample)>().sample("S1").open().unwrap();
    let mut seen = Vec::new();
    while let Some((info, sample)) = reader.next_record().unwrap() {
        let any = sample.get(&any).unwrap();
        seen.push((
            info.has(&info_dp).unwrap(),
            sample.has(&dp).unwrap(),
            sample.get(&dp).unwrap(),
            any.as_ref().map(DynValue::is_missing),
        ));
    }
    assert_eq!(
        seen,
        [
            (true, false, None, None),
            (true, true, None, Some(true)),
            (false, true, None, Some(true)),
            (false, true, Some(7), Some(false)),
        ]
    );
}

// ── A column reads as its samples ───────────────────────────────────────

/// The FORMAT fields of the column test, one per marker shape, in
/// dictionary order after `PASS`: `GT` is 1, `S` is 6.
const COLUMN_FIELDS: [(&str, &str, &str); 6] = [
    ("GT", "1", "String"),
    ("I1", "1", "Integer"),
    ("IN", ".", "Integer"),
    ("F1", "1", "Float"),
    ("FN", ".", "Float"),
    ("S", "1", "String"),
];

fn column_header(version: &str, n_samples: usize) -> String {
    use std::fmt::Write as _;
    let mut h = format!(
        "##fileformat=VCFv{version}\n##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
    );
    for (id, number, typ) in COLUMN_FIELDS {
        writeln!(h, "##FORMAT=<ID={id},Number={number},Type={typ},Description=\"x\">").unwrap();
    }
    h += "##contig=<ID=chr1,length=1000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for s in 0..n_samples {
        write!(h, "\tS{s}").unwrap();
    }
    h + "\n"
}

/// One BCF FORMAT column of `n` samples: any type, up to four elements a
/// sample, each element one a reader treats specially at its width (the
/// sentinels, the reserved band, GT values out of range, bad UTF-8).
#[hegel::composite]
fn bcf_column(tc: &TestCase, n: usize) -> Vec<u8> {
    let code = tc.draw(gs::sampled_from(vec![0u8, 1, 2, 3, 5, 7]));
    let per_sample = tc.draw(gs::integers::<usize>().max_value(4));
    let elements: Vec<Vec<u8>> = match code {
        1 => [0u8, 1, 2, 3, 4, 5, 0x7F, 0x80, 0x81, 0x82, 0xFB, 0xFF].map(|b| vec![b]).to_vec(),
        2 => [0u16, 2, 3, 4, 0x7FFF, 0x8000, 0x8001, 0x8002, 0xFFFF]
            .map(|v| v.to_le_bytes().to_vec())
            .to_vec(),
        3 => [0u32, 2, 3, 5, 131_074, 131_075, 131_076, 0x7FFF_FFFF, 0x8000_0000, 0x8000_0001]
            .into_iter()
            .chain([0x8000_0002, 0x8000_0007, 0x8000_0008, 0xFFFF_FFFF])
            .map(|v| v.to_le_bytes().to_vec())
            .collect(),
        5 => [0f32.to_bits(), 1.5f32.to_bits(), f32::NAN.to_bits(), 0x7F80_0001, 0x7F80_0002]
            .map(|v| v.to_le_bytes().to_vec())
            .to_vec(),
        7 => [b'a', b'.', b'0', b'/', 0, 0xC3, 0x28].map(|b| vec![b]).to_vec(),
        _ => vec![vec![]],
    };
    let mut bytes = vec![((per_sample as u8) << 4) | code];
    for _ in 0..per_sample * n {
        bytes.extend(tc.draw(gs::sampled_from(elements.clone())));
    }
    bytes
}

/// One VCF text sample: up to every field's cell, each a few tokens a
/// reader parses, or rejects, in some way.
#[hegel::composite]
fn text_sample(tc: &TestCase) -> String {
    let tokens = ["0", "1", "12", ".", "", "-1", "x", "65536", "2147483648", "-2147483648"]
        .into_iter()
        .chain(["0.5", "nan", "Inf", "a%3A"])
        .map(str::to_owned)
        .collect::<Vec<_>>();
    let n_cells = tc.draw(gs::integers::<usize>().min_value(1).max_value(COLUMN_FIELDS.len()));
    let cells: Vec<String> = (0..n_cells)
        .map(|_| {
            let mut cell = String::new();
            if tc.draw(gs::booleans()) {
                cell += tc.draw(gs::sampled_from(vec!["/", "|"]));
            }
            let n = tc.draw(gs::integers::<usize>().min_value(1).max_value(3));
            for i in 0..n {
                if i > 0 {
                    cell += tc.draw(gs::sampled_from(vec!["/", "|", ","]));
                }
                cell += &tc.draw(gs::sampled_from(tokens.clone()));
            }
            cell
        })
        .collect();
    cells.join(":")
}

/// `column(key)` reads what `sample(s).get(key)` reads for every sample,
/// and fails with the first sample's error when any sample's value is
/// malformed.
fn column_reads_as_each_sample<V: seqair::vcf::ReadKind>(
    samples: &Samples,
    name: &str,
    file: &VcfFile,
) where
    for<'a> V::Out<'a>: std::fmt::Debug,
{
    let key: seqair::vcf::FormatKey<V> = file.format_key(name).unwrap();
    let ids: Vec<_> = file.header().sample_indices().collect();
    let each: Result<Vec<String>, String> = ids
        .iter()
        .map(|&s| samples.sample(s).get(&key).map(|v| format!("{v:?}")).map_err(|e| e.to_string()))
        .collect();
    match (samples.column(&key), each) {
        (Ok(None), Ok(each)) => {
            assert!(!samples.has(&key).unwrap(), "{name}");
            let absent = format!("{:?}", V::absent());
            assert!(each.iter().all(|v| *v == absent), "{name}: {each:?}");
        }
        (Ok(Some(column)), Ok(each)) => {
            let by_get: Vec<String> =
                ids.iter().map(|&s| format!("{:?}", column.get(s).unwrap())).collect();
            assert_eq!(column.len(), ids.len(), "{name}");
            let read: Vec<String> = column.map(|v| format!("{v:?}")).collect();
            assert_eq!(read, each, "{name}");
            assert_eq!(by_get, each, "{name}");
        }
        (Err(e), Err(first)) => assert_eq!(e.to_string(), first, "{name}"),
        (column, each) => panic!("{name}: column {:?}, each sample {each:?}", column.map(|_| ())),
    }
}

// r[verify vcf_fields.samples+2]
// r[verify vcf_values.access+2]
/// A column reads, for every marker, what each sample's `get` reads, and
/// fails with the first malformed sample's error: BCF columns of every type
/// and integer width with sentinels and invalid values, and VCF text cells
/// of valid and malformed tokens, trailing fields dropped.
#[hegel::test]
fn a_column_reads_what_each_sample_reads(tc: TestCase) {
    let n = tc.draw(gs::integers::<usize>().min_value(1).max_value(4));
    let version = tc.draw(gs::sampled_from(vec!["4.3", "4.4"]));
    let header = column_header(version, n);
    let file = if tc.draw(gs::booleans()) {
        let mut indiv = Vec::new();
        for id in 1..=COLUMN_FIELDS.len() as u8 {
            indiv.extend([0x11, id]);
            indiv.extend(tc.draw(bcf_column(n)));
        }
        let rec = Rec {
            info: Vec::new(),
            n_info: 0,
            indiv,
            n_fmt: COLUMN_FIELDS.len() as u8,
            n_sample: n as u32,
            ..Rec::default()
        };
        VcfFile::from_bytes(bcf_with(&header, &[rec]), None).unwrap()
    } else {
        let keys: Vec<&str> = COLUMN_FIELDS.iter().map(|f| f.0).collect();
        let mut line = format!("chr1\t11\t.\tA\t.\t.\tPASS\t.\t{}", keys.join(":"));
        for _ in 0..n {
            line += "\t";
            line += &tc.draw(text_sample());
        }
        VcfFile::from_bytes(format!("{header}{line}\n").into_bytes(), None).unwrap()
    };
    let mut reader = file.reader::<Samples>().open().unwrap();
    let samples = reader.next_record().unwrap().unwrap();
    use seqair::vcf::{Any, Arr, Gt, OptArr, Scalar, Str};
    column_reads_as_each_sample::<Gt>(samples, "GT", &file);
    column_reads_as_each_sample::<Scalar<i32>>(samples, "I1", &file);
    column_reads_as_each_sample::<Arr<i32>>(samples, "IN", &file);
    column_reads_as_each_sample::<OptArr<i32>>(samples, "IN", &file);
    column_reads_as_each_sample::<OptArr<i32>>(samples, "I1", &file);
    column_reads_as_each_sample::<Scalar<f32>>(samples, "F1", &file);
    column_reads_as_each_sample::<Arr<f32>>(samples, "FN", &file);
    column_reads_as_each_sample::<OptArr<f32>>(samples, "FN", &file);
    column_reads_as_each_sample::<Str>(samples, "S", &file);
    column_reads_as_each_sample::<Str>(samples, "GT", &file);
    for (name, ..) in COLUMN_FIELDS {
        column_reads_as_each_sample::<Any>(samples, name, &file);
    }
}
