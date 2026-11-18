//! Region queries and partitions (`docs/spec/2-vcf-1-reader.md`, "Regions",
//! "Region queries", "Partitions") against a sequential read filtered in
//! the test and against `bcftools view -r`/`-t`, on files seqair's writer
//! makes (VCF.gz with seqair's or bcftools' CSI or TBI; BCF with CSI) and
//! on the fixtures.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    reason = "test code"
)]

#[path = "../src/pinned.rs"]
mod pinned;

use hegel::prelude::*;
use seqair::vcf::reader::ReadOption;
use seqair::vcf::{
    Alleles, ContigDef, ContigId, IndexKind, Info, InfoFieldDef, InfoInt, InfoString, Locus, Match,
    Number, OutputFormat, Partition, PlanError, Region, Scalar, Site, Str, VcfFile, VcfHeader,
    VcfPos, VcfReadError, Writer,
};
use seqair_types::{Pos0, Pos1, SmolStr};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::process::Command;

const BCFTOOLS: &str = "/opt/homebrew/bin/bcftools";

fn bcftools_bin() -> Option<&'static str> {
    let bin = if Path::new(BCFTOOLS).exists() { BCFTOOLS } else { "bcftools" };
    Command::new(bin).arg("--version").output().ok().map(|_| bin)
}

/// `bcftools args…` stdout; panics on failure.
fn bcftools(bin: &str, args: &[&str]) -> String {
    let out = Command::new(bin).args(args).output().unwrap();
    assert!(
        out.status.success(),
        "bcftools {args:?}: {}",
        String::from_utf8(out.stderr.clone()).unwrap_or_default()
    );
    String::from_utf8(out.stdout).unwrap()
}

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/vcf").join(name)
}

// ── Generated files ─────────────────────────────────────────────────────

const CONTIGS: [&str; 4] = ["chrA", "chrB", "chrC", "chrEmpty"];
const CONTIG_LEN: u32 = 1_000_000;

/// Who indexes a generated file.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Indexer {
    /// The writer's co-produced CSI.
    Seqair,
    /// The writer's co-produced TBI (`.vcf.gz` only).
    SeqairTbi,
    BcftoolsCsi,
    BcftoolsTbi,
}

#[derive(Debug, Clone)]
struct Rec {
    contig: usize,
    pos: u32,
    ref_len: usize,
    /// INFO END, 1-based inclusive, for a `<DEL>`.
    end: Option<u32>,
    pad: usize,
}

struct Generated {
    _dir: tempfile::TempDir,
    path: PathBuf,
}

#[hegel::composite]
fn records(tc: &TestCase, with_end: bool) -> Vec<Rec> {
    let n = tc.draw(gs::integers::<usize>().max_value(80));
    let mut recs: Vec<Rec> = (0..n)
        .map(|_| {
            let pos = tc.draw(gs::integers::<u32>().min_value(1).max_value(120_000));
            let end = (with_end && tc.draw(gs::integers::<u8>().max_value(5)) == 0)
                .then(|| pos + tc.draw(gs::integers::<u32>().max_value(40_000)));
            Rec {
                // chrEmpty is declared and never used.
                contig: tc.draw(gs::integers::<usize>().max_value(2)),
                pos,
                ref_len: tc.draw(gs::integers::<usize>().min_value(1).max_value(6)),
                end,
                pad: tc.draw(gs::sampled_from(vec![0, 10, 300, 3000])),
            }
        })
        .collect();
    recs.sort_by_key(|r| (r.contig, r.pos));
    recs
}

/// Write `recs` (record `i` has INFO `PAD` = `r{i}` and padding, so blocks
/// fill and every record is told apart) and index it with `indexer`.
fn generate(recs: &[Rec], format: OutputFormat, indexer: Indexer, bin: &str) -> Generated {
    let mut b = VcfHeader::builder();
    let contigs: Vec<ContigId> = CONTIGS
        .iter()
        .map(|c| b.register_contig(*c, ContigDef::new(Some(u64::from(CONTIG_LEN)))).unwrap())
        .collect();
    let mut b = b.filters().infos();
    let end: InfoInt = b.register_info(&InfoFieldDef::<Scalar<i32>>::new("END", "end")).unwrap();
    let pad: InfoString =
        b.register_info(&InfoFieldDef::<Str>::new("PAD", Number::Count(1), "pad")).unwrap();
    let header = b.build().unwrap();

    let mut out = Vec::new();
    let kind = if indexer == Indexer::SeqairTbi { IndexKind::Tbi } else { IndexKind::Csi };
    let mut w = Writer::new(&mut out, format).index_kind(kind).write_header(&header).unwrap();
    for (i, r) in recs.iter().enumerate() {
        let reference: String = "ACGTAC"[..r.ref_len].to_owned();
        let alt = if r.end.is_some() { "<DEL>" } else { "T" };
        let alleles =
            Alleles::complex(SmolStr::from(reference), [SmolStr::from(alt)].into_iter().collect());
        let enc =
            w.begin_record(&contigs[r.contig], Pos1::new(r.pos).unwrap(), &alleles, None).unwrap();
        let mut enc = enc.filter_pass();
        if let Some(e) = r.end {
            end.encode(&mut enc, e as i32).unwrap();
        }
        pad.encode(&mut enc, &format!("r{i}{}", "x".repeat(r.pad))).unwrap();
        enc.emit().unwrap();
    }
    let (_, index) = w.finish().unwrap();

    let dir = tempfile::tempdir().unwrap();
    let ext = if format == OutputFormat::Bcf { "bcf" } else { "vcf.gz" };
    let path = dir.path().join(format!("x.{ext}"));
    std::fs::write(&path, &out).unwrap();
    let p = path.to_str().unwrap();
    match indexer {
        Indexer::Seqair | Indexer::SeqairTbi => {
            let index = index.unwrap();
            index.write_to_path(format!("{p}.{}", index.suffix())).unwrap();
        }
        Indexer::BcftoolsCsi => {
            bcftools(bin, &["index", "-c", p]);
        }
        Indexer::BcftoolsTbi => {
            bcftools(bin, &["index", "-t", p]);
        }
    }
    Generated { _dir: dir, path }
}

/// A record as the tests compare it: contig, 0-based start, 0-based
/// exclusive end, and its INFO PAD tag (unique per record).
type Key = (SmolStr, i64, i64, String);

fn key(locus: &Locus, info: &Info, pad: &InfoString) -> Key {
    let at = |p: VcfPos| match p {
        VcfPos::At(p) => i64::from(p.as_u32()),
        VcfPos::Telomere => 0,
    };
    let tag = info.get(pad).unwrap().unwrap();
    let tag = tag.trim_end_matches('x').to_owned();
    (locus.contig().name().clone(), at(locus.start()) - 1, at(locus.end()), tag)
}

fn drain(reader: &mut seqair::vcf::Reader<(Locus, Info)>, pad: &InfoString) -> Vec<Key> {
    let mut out = Vec::new();
    while let Some((locus, info)) = reader.next_record().unwrap() {
        out.push(key(locus, info, pad));
    }
    out
}

fn sequential(file: &VcfFile, pad: &InfoString) -> Vec<Key> {
    drain(&mut file.reader::<(Locus, Info)>().open().unwrap(), pad)
}

/// The test's own reading of `r[vcf_reader.query.overlap+2]`: whether a
/// record `[start, end)` is selected by the 0-based half-open `[from, to)`.
fn selects(matching: Match, (start, end): (i64, i64), (from, to): (i64, i64)) -> bool {
    match matching {
        Match::StartsIn => from <= start && start < to,
        Match::Overlapping => from < to && start < to && end.max(start + 1) > from,
    }
}

/// A region as drawn: contig, 0-based `[from, to)` (`to` `None` for the
/// contig's end), and whether it is written 1-based.
#[derive(Debug, Clone)]
struct Drawn {
    contig: usize,
    from: u32,
    to: Option<u32>,
    one_based: bool,
}

impl Drawn {
    fn region(&self) -> Region {
        let c = CONTIGS[self.contig];
        let r = match (self.one_based, self.to) {
            (false, Some(to)) => {
                Region::new(c, Pos0::new(self.from).unwrap()..Pos0::new(to).unwrap())
            }
            (false, None) => Region::new(c, Pos0::new(self.from).unwrap()..),
            // [from, to) 0-based is from+1..=to 1-based; an empty one is
            // written half-open, as a..a.
            (true, Some(to)) if to == self.from => {
                Region::new(c, Pos1::new(self.from + 1).unwrap()..Pos1::new(to + 1).unwrap())
            }
            (true, Some(to)) => {
                Region::new(c, Pos1::new(self.from + 1).unwrap()..=Pos1::new(to).unwrap())
            }
            (true, None) => Region::new(c, Pos1::new(self.from + 1).unwrap()..),
        };
        r.unwrap()
    }

    fn span(&self) -> (i64, i64) {
        (i64::from(self.from), self.to.map_or(i64::MAX, i64::from))
    }

    /// The region as bcftools reads it: 1-based inclusive.
    fn bcftools(&self) -> String {
        format!("{}:{}-{}", CONTIGS[self.contig], self.from + 1, self.to.unwrap_or(CONTIG_LEN))
    }
}

#[hegel::composite]
fn drawn(tc: &TestCase) -> Drawn {
    let from = tc.draw(gs::integers::<u32>().max_value(170_000));
    let to = match tc.draw(gs::integers::<u8>().max_value(9)) {
        0 => None,
        1 => Some(from),
        _ => Some(from + tc.draw(gs::integers::<u32>().max_value(60_000))),
    };
    Drawn {
        contig: tc.draw(gs::integers::<usize>().max_value(3)),
        from,
        to,
        one_based: tc.draw(gs::booleans()),
    }
}

#[hegel::composite]
fn generated_file(tc: &TestCase) -> (Vec<Rec>, OutputFormat, Indexer) {
    let format = tc.draw(gs::sampled_from(vec![OutputFormat::VcfGz, OutputFormat::Bcf]));
    let indexers = match format {
        OutputFormat::Bcf => vec![Indexer::Seqair, Indexer::BcftoolsCsi],
        _ => vec![Indexer::Seqair, Indexer::SeqairTbi, Indexer::BcftoolsCsi, Indexer::BcftoolsTbi],
    };
    let indexer = tc.draw(gs::sampled_from(indexers));
    // Every index, the writer's own included, spans INFO END.
    let recs = tc.draw(records(true));
    (recs, format, indexer)
}

// r[verify vcf_reader.partitions.complete+2]
// r[verify vcf_reader.partitions+2]
// r[verify vcf_reader.partitions.ownership+2]
/// The partitions of a file, read one after another, are the file:
/// nothing missing, nothing twice, in order.
#[hegel::test(test_cases = pinned::cases(20))]
fn partitions_concatenate_to_the_file(tc: TestCase) {
    let Some(bin) = bcftools_bin() else { return };
    let (recs, format, indexer) = tc.draw(generated_file());
    let g = generate(&recs, format, indexer, bin);
    let file = VcfFile::open(&g.path).unwrap();
    assert!(file.has_index());
    let pad: InfoString = file.header().info_key("PAD").unwrap();
    let n = NonZeroUsize::new(tc.draw(gs::integers::<usize>().min_value(1).max_value(9))).unwrap();
    let partitions: Vec<Partition> = file.partitions(n).unwrap();
    let mut all = Vec::new();
    for p in &partitions {
        let got = drain(&mut file.reader::<(Locus, Info)>().partition(p).unwrap(), &pad);
        // Each partition's records start inside it.
        let (from, to) = (
            p.region().start().as_i64(),
            match p.region().end() {
                seqair::vcf::RegionEnd::ContigEnd => i64::MAX,
                seqair::vcf::RegionEnd::Before(e) => e.as_i64(),
            },
        );
        assert!(got.iter().all(|k| from <= k.1 && k.1 < to && k.0 == p.region().target().name()));
        all.extend(got);
    }
    assert_eq!(all, sequential(&file, &pad), "{partitions:?}");
    // Each indexed contig gets at least one partition, no unindexed one any.
    let used = recs.iter().map(|r| r.contig).collect::<std::collections::BTreeSet<_>>().len();
    let contigs = partitions
        .iter()
        .map(|p| p.region().target().name())
        .collect::<std::collections::BTreeSet<_>>();
    assert_eq!(contigs.len(), used);
}

// r[verify vcf_reader.query+2]
// r[verify vcf_reader.query.overlap+2]
// r[verify vcf_reader.query.many+2]
// r[verify vcf_reader.query.contigs+2]
// r[verify vcf_reader.region]
/// Random regions (0- and 1-based, open-ended, empty, on a contig without
/// records) yield what a sequential read filtered by the overlap rule
/// yields, under both match modes; `requery` the same on one reader;
/// `query_many` the union, each record once in file order; and bcftools
/// `view -r` (`Overlapping`) and `-t` (`StartsIn`) agree.
#[hegel::test(test_cases = pinned::cases(20))]
fn queries_equal_a_filtered_sequential_read(tc: TestCase) {
    let Some(bin) = bcftools_bin() else { return };
    let (recs, format, indexer) = tc.draw(generated_file());
    let g = generate(&recs, format, indexer, bin);
    let file = VcfFile::open(&g.path).unwrap();
    let pad: InfoString = file.header().info_key("PAD").unwrap();
    let all = sequential(&file, &pad);
    let regions = tc.draw(gs::vecs(drawn()).min_size(1).max_size(4));
    let matching = tc.draw(gs::sampled_from(vec![Match::Overlapping, Match::StartsIn]));
    let expect = |rs: &[Drawn]| -> Vec<Key> {
        all.iter()
            .filter(|k| {
                rs.iter()
                    .any(|r| k.0 == CONTIGS[r.contig] && selects(matching, (k.1, k.2), r.span()))
            })
            .cloned()
            .collect()
    };

    let builder = || file.reader::<(Locus, Info)>().matching(matching);
    let mut reader: Option<seqair::vcf::Reader<(Locus, Info)>> = None;
    let p = g.path.to_str().unwrap();
    for r in &regions {
        let want = expect(std::slice::from_ref(r));
        let got = drain(&mut builder().query(r.region()).unwrap(), &pad);
        assert_eq!(got, want, "{r:?} {matching:?}");
        // requery, on one reader across all regions.
        match &mut reader {
            Some(reader) => reader.requery(r.region()).unwrap(),
            None => reader = Some(builder().query(r.region()).unwrap()),
        }
        assert_eq!(drain(reader.as_mut().unwrap(), &pad), want, "requery {r:?}");
        // bcftools, on a non-empty region.
        if r.to != Some(r.from) {
            let flag = if matching == Match::Overlapping { "-r" } else { "-t" };
            let text = bcftools(bin, &["view", "-H", flag, &r.bcftools(), p]);
            let tags: Vec<String> = text
                .lines()
                .map(|l| {
                    let info = l.split('\t').nth(7).unwrap();
                    let pad = info.split(';').find_map(|e| e.strip_prefix("PAD=")).unwrap();
                    pad.trim_end_matches('x').to_owned()
                })
                .collect();
            let ours: Vec<String> = want.iter().map(|k| k.3.clone()).collect();
            assert_eq!(ours, tags, "bcftools view {flag} {}", r.bcftools());
        }
    }
    let union = drain(&mut builder().query_many(regions.iter().map(Drawn::region)).unwrap(), &pad);
    assert_eq!(union, expect(&regions), "query_many {regions:?}");
    let mut reader = reader.unwrap();
    reader.requery_many(regions.iter().map(Drawn::region)).unwrap();
    assert_eq!(drain(&mut reader, &pad), union, "requery_many");
}

// ── Fixtures ────────────────────────────────────────────────────────────

/// `(CHROM, POS, REF, ALT)` of each record.
fn site_keys(reader: &mut seqair::vcf::Reader<(Locus, Site)>) -> Vec<String> {
    let mut out = Vec::new();
    while let Some((locus, site)) = reader.next_record().unwrap() {
        let VcfPos::At(pos) = locus.start() else { panic!("telomere") };
        let alts: Vec<String> = site
            .alleles()
            .alts()
            .map(|a| String::from_utf8(a.as_bytes().to_vec()).unwrap())
            .collect();
        let alts = if alts.is_empty() { ".".to_owned() } else { alts.join(",") };
        let reference = String::from_utf8(site.alleles().reference().to_vec()).unwrap();
        out.push(format!("{}\t{pos}\t{reference}\t{alts}", locus.contig().name()));
    }
    out
}

fn bcftools_keys(bin: &str, flag: &str, region: &str, path: &Path) -> Vec<String> {
    bcftools(bin, &["view", "-H", flag, region, path.to_str().unwrap()])
        .lines()
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            format!("{}\t{}\t{}\t{}", c[0], c[1], c[3], c[4])
        })
        .collect()
}

// r[verify vcf_reader.query+2]
// r[verify vcf_reader.query.overlap+2]
// r[verify vcf_reader.oracles]
/// On every fixture (TBI and CSI), `query` equals `bcftools view -r` and
/// `StartsIn` equals `-t`, over regions inside, across and between records.
#[test]
fn fixture_queries_equal_bcftools() {
    let Some(bin) = bcftools_bin() else { return };
    let cases: &[(&str, &[&str])] = &[
        (
            "truth.vcf.gz",
            &[
                "chr1:600000-800000",
                "chr1:1-700000",
                "chr1:1000000-1500000",
                "chr2:1-30000",
                "chr2:390000-400000",
                "chr1:1482316-1482316",
                "chr3:1-100",
            ],
        ),
        ("rastair.vcf.gz", &["chr1:10009-10100", "chr2:10090-10500", "chr1:1-1000000"]),
        ("rastair.bcf", &["chr1:10009-10100", "chr2:10090-10500", "chr2:1-10000"]),
        ("panel.bcf", &["chr22:20000241-20010000", "chr22:20020000-20049792"]),
    ];
    for (name, regions) in cases {
        let path = fixture(name);
        let file = VcfFile::open(&path).unwrap();
        assert!(file.has_index(), "{name}");
        for region in *regions {
            for (matching, flag) in [(Match::Overlapping, "-r"), (Match::StartsIn, "-t")] {
                let r: Region = region.parse().unwrap();
                let mut q = file.reader::<(Locus, Site)>().matching(matching).query(&r).unwrap();
                assert_eq!(
                    site_keys(&mut q),
                    bcftools_keys(bin, flag, region, &path),
                    "{name} {flag} {region}"
                );
            }
        }
    }
}

// r[verify vcf_reader.query.contigs+2]
// r[verify vcf_reader.file+2]
/// A declared contig without records is an empty query; an unknown name is
/// `UnknownContig`; another file's `ContigId` is `ForeignKey`; a file
/// without an index reads front to back and fails only queries.
#[test]
fn query_errors() {
    let file = VcfFile::open(fixture("truth.vcf.gz")).unwrap();
    let mut q = file.reader::<Locus>().query(Region::contig("chr3")).unwrap();
    assert!(q.next_record().unwrap().is_none());
    assert!(matches!(
        file.reader::<Locus>().query(Region::contig("chrUn_nope")),
        Err(VcfReadError::UnknownContig { name }) if name == "chrUn_nope"
    ));
    let other = VcfFile::open(fixture("rastair.vcf.gz")).unwrap();
    let foreign = other.header().contig("chr1").unwrap();
    assert!(matches!(
        file.reader::<Locus>().query(Region::contig(&foreign)),
        Err(VcfReadError::ForeignKey { .. })
    ));
    // A ContigId of the file itself works.
    let own = file.header().contig("chr2").unwrap();
    let n = file.reader::<Locus>().query(Region::contig(&own)).unwrap().into_iter().count();
    assert_eq!(n, 497);

    let dir = tempfile::tempdir().unwrap();
    let bare = dir.path().join("bare.vcf.gz");
    std::fs::copy(fixture("truth.vcf.gz"), &bare).unwrap();
    let bare = VcfFile::open(&bare).unwrap();
    assert!(!bare.has_index());
    assert_eq!(bare.reader::<Locus>().open().unwrap().into_iter().count(), 1237);
    assert!(matches!(
        bare.reader::<Locus>().query(Region::contig("chr1")),
        Err(VcfReadError::NoIndex)
    ));
    assert!(matches!(bare.partitions(NonZeroUsize::MIN), Err(VcfReadError::NoIndex)));
    let mut reader = bare.reader::<Locus>().open().unwrap();
    assert!(matches!(reader.requery(Region::contig("chr1")), Err(VcfReadError::NoIndex)));
}

// r[verify vcf_reader.query.options]
#[test]
fn matching_is_a_query_option() {
    let file = VcfFile::open(fixture("rastair.bcf")).unwrap();
    let unused = |r: Result<_, VcfReadError>| {
        matches!(
            r,
            Err(VcfReadError::Plan(PlanError::UnusedOption { option: ReadOption::Matching }))
        )
    };
    assert!(unused(file.reader::<Locus>().matching(Match::StartsIn).open()));
    let p = file.partitions(NonZeroUsize::MIN).unwrap();
    assert!(unused(file.reader::<Locus>().matching(Match::StartsIn).partition(&p[0])));
}

// r[verify vcf_reader.file+2]
/// A VCF.gz with only a `.tbi` and one with only a `.csi` both open indexed
/// and query alike; `from_bytes` tells them apart by content.
#[test]
fn tbi_only_and_csi_only() {
    let Some(bin) = bcftools_bin() else { return };
    let dir = tempfile::tempdir().unwrap();
    let csi = dir.path().join("csi.vcf.gz");
    std::fs::copy(fixture("truth.vcf.gz"), &csi).unwrap();
    bcftools(bin, &["index", "-c", csi.to_str().unwrap()]);
    let tbi = fixture("truth.vcf.gz");
    let region: Region = "chr1:600000-800000".parse().unwrap();
    let mut results = Vec::new();
    for (path, index) in
        [(&tbi, format!("{}.tbi", tbi.display())), (&csi, format!("{}.csi", csi.display()))]
    {
        let file = VcfFile::open(path).unwrap();
        assert!(file.has_index(), "{}", path.display());
        results.push(
            file.reader::<Locus>()
                .query(&region)
                .unwrap()
                .into_iter()
                .collect::<Result<Vec<_>, _>>()
                .unwrap()
                .len(),
        );
        let bytes =
            VcfFile::from_bytes(std::fs::read(path).unwrap(), Some(&std::fs::read(index).unwrap()))
                .unwrap();
        assert!(bytes.has_index());
        results.push(bytes.reader::<Locus>().query(&region).unwrap().into_iter().count());
    }
    assert!(results[0] > 0 && results.iter().all(|&n| n == results[0]), "{results:?}");
}

// r[verify vcf_reader.fork_equivalence]
/// Threads querying one `VcfFile` at once — partitions and regions, each
/// on its own reader — read what one thread reads.
#[test]
fn concurrent_queries_equal_a_single_thread() {
    let file = VcfFile::open(fixture("truth.vcf.gz")).unwrap();
    let regions: Vec<Region> =
        ["chr1:600000-800000", "chr1:1-1500000", "chr2:100000-300000", "chr2"]
            .iter()
            .map(|r| r.parse().unwrap())
            .collect();
    let partitions = file.partitions(NonZeroUsize::new(6).unwrap()).unwrap();
    let read = |f: &VcfFile| -> Vec<Vec<Locus>> {
        let mut out: Vec<Vec<Locus>> = regions
            .iter()
            .map(|r| {
                f.reader::<Locus>().query(r).unwrap().into_iter().map(Result::unwrap).collect()
            })
            .collect();
        out.extend(partitions.iter().map(|p| {
            f.reader::<Locus>().partition(p).unwrap().into_iter().map(Result::unwrap).collect()
        }));
        out
    };
    let single = read(&file);
    assert_eq!(single[4..].iter().map(Vec::len).sum::<usize>(), 1237);
    std::thread::scope(|s| {
        let handles: Vec<_> = (0..8).map(|_| s.spawn(|| read(&file))).collect();
        for h in handles {
            assert_eq!(h.join().unwrap(), single);
        }
    });
}
