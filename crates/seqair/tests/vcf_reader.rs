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

use seqair::io::BgzfWriter;
use seqair::vcf::reader::{
    FieldError, FilterStatus, KeyKind, MAX_RECORD_BYTES, ReadOption, RecordSection,
};
use seqair::vcf::{
    Allele, DynValue, Fields, FormatGt, FullRecord, Info, InfoFlag, InfoInt, KeyError, Locus,
    OnlySample, PlanError, RawRecord, ReadOptions, RecordError, Sample, Samples, Sections, Site,
    VcfFile, VcfHeader, VcfPos, VcfReadError,
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
    match a {
        Allele::Bases(b) => String::from_utf8(b.to_vec()).unwrap(),
        Allele::Symbolic(s) | Allele::Breakend(s) => s.to_owned(),
        Allele::Overlapping => "*".to_owned(),
    }
}

// r[verify vcf_reader.oracles]
// r[verify vcf_values.access]
// r[verify vcf_values.dynamic]
// r[verify vcf_values.qual]
// r[verify vcf_fields.full_record]
/// Every field of every record of a BCF fixture, read through `FullRecord`,
/// equals what `bcftools view` prints for it; the end equals `%END`.
fn check_against_bcftools(name: &str) {
    let path = fixture(name);
    let path_str = path.to_str().unwrap();
    let Some(text) = bcftools(&["view", "-H", path_str]) else { return };
    let Some(ends) = bcftools(&["query", "-f", "%END\n", path_str]) else { return };
    let file = VcfFile::open(&path).unwrap();
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

        let info: Vec<_> = rec.info().iter().unwrap().collect::<Result<_, _>>().unwrap();
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
                    rec.samples().sample(sample).iter().unwrap().collect::<Result<_, _>>().unwrap();
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

// ── Fixtures through every way of reading ──────────────────────────────

// r[verify vcf_reader.reader]
// r[verify vcf_reader.read_into+2]
// r[verify vcf_reader.next_record]
// r[verify vcf_reader.records]
// r[verify vcf_reader.passthrough_header]
// r[verify vcf_fields.owned]
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
// r[verify vcf_reader.file.shared]
// r[verify vcf_reader.reader.threads]
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
                        let two = std::num::NonZeroU8::new(2).unwrap();
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
    fn send<T: Send>() {}
    send_sync_clone::<VcfFile>();
    send::<seqair::vcf::Reader<FullRecord>>();
    send_sync_clone::<FullRecord>();
}

// r[verify vcf_fields.sample]
// r[verify vcf_fields.samples]
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
// r[verify vcf_fields.unused_options]
// r[verify vcf_fields.only_sample]
// r[verify vcf_fields.sample]
#[test]
fn plans_fail_at_open() {
    let panel = VcfFile::open(fixture("panel.bcf")).unwrap();
    let rastair = VcfFile::open(fixture("rastair.bcf")).unwrap();
    assert_eq!(plan_error::<OnlySample>(&panel, None), PlanError::NotSingleSample { n: 4 });
    assert_eq!(plan_error::<Sample>(&panel, None), PlanError::SampleNotSelected { n: 4 });
    assert_eq!(
        plan_error::<Sample>(&panel, Some("nope")),
        PlanError::UnknownSample { name: "nope".into() }
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
}

#[derive(Clone, Debug)]
struct BadKey;

impl Fields for BadKey {
    type Plan = InfoInt;
    fn plan(header: &VcfHeader, _: &ReadOptions) -> Result<InfoInt, PlanError> {
        Ok(header.info_key("BQ")?) // declared Float
    }
    fn needs(_: &InfoInt) -> Sections {
        Sections::INFO
    }
    fn new(_: &InfoInt, _: RawRecord<'_>) -> Result<Self, VcfReadError> {
        Ok(Self)
    }
}

// ── Keys of another header ──────────────────────────────────────────────

// r[verify vcf_values.header_brand]
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
    fn plan(header: &VcfHeader, _: &ReadOptions) -> Result<InfoFlag, PlanError> {
        Ok(header.info_key("CPG")?)
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

// r[verify vcf_fields.trait+2]
// r[verify vcf_fields.user_types]
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
        fn plan(_: &VcfHeader, _: &ReadOptions) -> Result<(), PlanError> {
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
        shared.extend_from_slice(&(1 | (u32::from(self.n_fmt) << 24)).to_le_bytes());
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
    let mut raw = b"BCF\x02\x02".to_vec();
    let text = format!("{HEADER}\0");
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
        VcfReadError::Record { ordinal, offset, source } => {
            assert!(offset.is_some(), "a BGZF record error names its offset");
            (ordinal, source)
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

// r[verify vcf_reader.limits]
// r[verify vcf_reader.errors]
// r[verify vcf_reader.undeclared]
#[test]
fn malformed_records_are_typed_errors_with_their_place() {
    let ok = Rec::default();
    let (ordinal, e) =
        record_error(first_error::<Locus>(&[ok.clone(), Rec { rid: 3, ..ok.clone() }]));
    assert_eq!((ordinal, e), (1, RecordError::UnknownContig { index: 3 }));
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

// r[verify vcf_fields.lazy]
// r[verify vcf_reader.undeclared]
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
            assert!(info.iter().is_err());
            assert!(e.contains("record 1"), "{e}");
        }
        values.push(first.ok().flatten());
    }
    // The index buffers survive the failure: the next record reads.
    assert_eq!(values, [Some(5), None, Some(5)]);
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
    assert!(matches!(VcfFile::from_bytes(gzip, None), Err(VcfReadError::NotBgzf)));
    let mut bcf21 = bcf(&[]);
    let mut w = BgzfWriter::new(Vec::new());
    w.write_all(b"BCF\x02\x01").unwrap();
    bcf21.clone_from(&w.finish().unwrap());
    assert!(matches!(VcfFile::from_bytes(bcf21, None), Err(VcfReadError::Bcf21)));
    assert!(matches!(
        VcfFile::from_bytes(HEADER.as_bytes().to_vec(), None),
        Err(VcfReadError::Unsupported { .. })
    ));
    assert!(matches!(
        VcfFile::from_bytes(b"hello".to_vec(), None),
        Err(VcfReadError::UnknownFormat)
    ));
    // A file with no records reads none.
    let empty = VcfFile::from_bytes(bcf(&[]), None).unwrap();
    assert!(empty.reader::<FullRecord>().open().unwrap().next_record().unwrap().is_none());
}

// r[verify vcf_reader.file]
#[test]
fn a_missing_index_is_fine_and_a_broken_one_fails_open() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("x.bcf");
    std::fs::write(&path, bcf(&[Rec::default()])).unwrap();
    let file = VcfFile::open(&path).unwrap();
    assert!(!file.has_index());
    assert_eq!(file.reader::<Locus>().open().unwrap().into_iter().count(), 1);
    std::fs::write(dir.path().join("x.bcf.csi"), b"not an index").unwrap();
    assert!(matches!(VcfFile::open(&path), Err(VcfReadError::Index(_))));
}
