//! `IndexedCramReader::query_raw` (through `IndexedReader::query_raw`) on the
//! CRAM fixtures, against two oracles: what `samtools view -b <cram> <region>`
//! writes, and the records the stream decoder yields for the whole file,
//! filtered to the region by a CIGAR walk of their own.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    reason = "test code"
)]

#[path = "../src/pinned.rs"]
mod pinned;

use core::range::RangeInclusive;
use hegel::TestCase;
use hegel::generators as gs;
use seqair::bam::stream::BamStreamReader;
use seqair::bam::{RawBamRecord, RecordStore};
use seqair::cram::stream::{CramStreamError, CramStreamReader};
use seqair::reader::{IndexPath, IndexedReader, ReaderError, ResolveTid};
use seqair::sam::encode::RecordEncodeError;
use seqair_types::{IupacBase, Pos0};
use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::OnceLock;

fn data(name: &str) -> PathBuf {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data")).join(name)
}

fn run(cmd: &mut Command) -> Vec<u8> {
    let out = cmd.output().expect("tool on PATH");
    assert!(out.status.success(), "{cmd:?}: {}", out.stderr.escape_ascii());
    out.stdout
}

/// A CRAM, its reference and its `.crai`, with every record the stream
/// decoder makes of the whole file.
struct Fixture {
    cram: PathBuf,
    fasta: PathBuf,
    crai: PathBuf,
    streamed: Vec<RawBamRecord>,
}

impl Fixture {
    fn new(cram: PathBuf, fasta: PathBuf, crai: PathBuf) -> Self {
        let streamed = stream(&cram, &fasta);
        Self { cram, fasta, crai, streamed }
    }

    fn open(&self) -> IndexedReader {
        IndexedReader::open_with_reference_and_index(
            &self.cram,
            &self.fasta,
            &IndexPath::new(&self.crai),
        )
        .unwrap()
    }

    fn name(&self) -> String {
        self.cram.file_stem().unwrap().display().to_string()
    }
}

/// Every record of `cram`, named as `samtools view <cram>` names them.
fn stream(cram: &Path, fasta: &Path) -> Vec<RawBamRecord> {
    let mut reader = CramStreamReader::new(std::fs::File::open(cram).unwrap()).unwrap();
    let prefix = cram.file_name().unwrap().to_str().unwrap();
    let mut decoder = reader.decoder(fasta).unwrap().name_prefix(prefix);
    let mut records = Vec::new();
    while let Some(container) = reader.next_container().unwrap() {
        decoder.decode(&container, &mut records).unwrap();
    }
    records
}

/// Every CRAM fixture of `tests/cram_stream.rs`, the hts-specs files, which
/// come without a `.crai`, indexed into a temporary directory, and, last,
/// [`shapes`].
fn fixtures() -> &'static [Fixture] {
    static FIXTURES: OnceLock<(Vec<Fixture>, tempfile::TempDir, Generated)> = OnceLock::new();
    let (fixtures, _, _) = FIXTURES.get_or_init(|| {
        let dir = tempfile::tempdir().unwrap();
        let beside = |cram: PathBuf, fasta: PathBuf| {
            let crai = PathBuf::from(format!("{}.crai", cram.display()));
            Fixture::new(cram, fasta, crai)
        };
        let test_fa = data("test.fasta.gz");
        let fidelity_fa = data("cram-stream/fidelity.fa");
        let mut out = vec![
            beside(data("test.cram"), test_fa.clone()),
            beside(data("test_v30.cram"), test_fa.clone()),
            beside(data("test_gzip.cram"), test_fa),
            beside(data("unplaced_multiref.cram"), data("unplaced_multiref.fa")),
        ];
        for f in [
            "fidelity_v31",
            "fidelity_v30",
            "fidelity_slices",
            "fidelity_multiref",
            "fidelity_lossy_names",
            "fidelity_embed_ref2",
            "fidelity_no_ref",
            "onlyu",
        ] {
            out.push(beside(data(&format!("cram-stream/{f}.cram")), fidelity_fa.clone()));
        }
        for level in 1..=4 {
            let cram = data(&format!("hts-specs/cram/3.1/passed/level-{level}.cram"));
            let crai = dir.path().join(format!("level-{level}.cram.crai"));
            run(Command::new("samtools").arg("index").arg("-o").arg(&crai).arg(&cram));
            out.push(Fixture::new(cram, data("hts-specs/cram/ce.fa"), crai));
        }
        let shapes = generate_with(&["seqs_per_slice=5"], shapes);
        out.push(Fixture::new(shapes.cram.clone(), shapes.fasta.clone(), shapes.crai.clone()));
        let unmapped = out.last().unwrap().streamed.iter().filter(|r| r.flags().is_unmapped());
        assert_eq!((out.last().unwrap().streamed.len(), unmapped.count()), (72, 12));
        (out, dir, shapes)
    });
    fixtures
}

/// Last reference base a record covers, from its decoded CIGAR: its POS when
/// it is unmapped or spans no reference (htslib's `bam_endpos` minus one).
fn last_base(record: &RawBamRecord, pos: u32) -> u32 {
    if record.flags().is_unmapped() {
        return pos;
    }
    let rlen: u32 = record
        .cigar()
        .map(|op| op.to_bam_u32())
        .filter(|op| matches!(op & 0xF, 0 | 2 | 3 | 7 | 8))
        .map(|op| op >> 4)
        .sum();
    pos + rlen.saturating_sub(1)
}

/// The streamed records on `tid` that overlap `span`, in file order.
fn filtered(records: &[RawBamRecord], tid: u32, span: RangeInclusive<Pos0>) -> Vec<&[u8]> {
    let (start, last) = (span.start.as_u32(), span.last.as_u32());
    records
        .iter()
        .filter(|r| {
            let (Some(t), Some(pos)) = (r.tid(), r.pos()) else { return false };
            let pos = pos.as_u32();
            t.as_u32() == tid && pos <= last && last_base(r, pos) >= start
        })
        .map(RawBamRecord::as_payload)
        .collect()
}

/// `query_raw`'s records, read into one reused record, checked against the
/// same query's iterator.
fn query(reader: &mut IndexedReader, tid: u32, span: RangeInclusive<Pos0>) -> Vec<RawBamRecord> {
    let tid = tid.resolve_tid(reader.header()).unwrap();
    let mut q = reader.query_raw(tid, span).unwrap();
    let mut got = Vec::new();
    let mut record = RawBamRecord::default();
    while q.read_record(&mut record).unwrap() {
        got.push(record.clone());
    }
    let last = record.clone();
    assert!(!q.read_record(&mut record).unwrap(), "a spent query stays spent");
    assert!(record == last, "the end leaves the record as it was");
    drop(q);
    let iterated: Vec<RawBamRecord> =
        reader.query_raw(tid, span).unwrap().collect::<Result<_, _>>().unwrap();
    assert!(iterated == got, "read_record and the iterator differ");
    got
}

fn payloads(records: &[RawBamRecord]) -> Vec<&[u8]> {
    records.iter().map(RawBamRecord::as_payload).collect()
}

/// samtools' records for `region` of `cram`, as BAM payloads.
fn samtools_region(f: &Fixture, region: &str) -> Vec<Vec<u8>> {
    let out = run(Command::new("samtools")
        .args(["view", "--no-PG", "-b", "--reference"])
        .arg(&f.fasta)
        .arg("-X")
        .arg(&f.cram)
        .arg(&f.crai)
        .arg(region));
    let mut reader = BamStreamReader::new(&out[..]).unwrap();
    let mut payloads = Vec::new();
    loop {
        let mut payload = Vec::new();
        if !reader.read_payload_into(&mut payload).unwrap() {
            break;
        }
        payloads.push(payload);
    }
    payloads
}

fn region_string(name: &str, span: RangeInclusive<Pos0>) -> String {
    format!("{name}:{}-{}", span.start.as_u32() + 1, span.last.as_u32() + 1)
}

fn span(start: u32, last: u32) -> RangeInclusive<Pos0> {
    RangeInclusive { start: Pos0::new(start).unwrap(), last: Pos0::new(last).unwrap() }
}

/// First and last base of every streamed record on `tid`.
fn ends(f: &Fixture, tid: u32) -> Vec<(u32, u32)> {
    f.streamed
        .iter()
        .filter(|r| r.tid().is_some_and(|t| t.as_u32() == tid))
        .filter_map(|r| r.pos().map(|p| (p.as_u32(), last_base(r, p.as_u32()))))
        .collect()
}

// r[verify cram.query_raw]
// r[verify cram.query_raw.htslib]
// r[verify unified.query_raw+2]
/// On every fixture and every reference: the whole reference, and a base
/// at, before and after the first and last base of a few records, give the
/// records `samtools view -b` writes for that region, byte for byte.
#[test]
fn raw_regions_are_what_samtools_writes() {
    let mut regions = 0;
    for f in fixtures() {
        let mut reader = f.open();
        let header = reader.header().clone();
        for tid in 0..header.target_count() as u32 {
            let resolved = tid.resolve_tid(&header).unwrap();
            let name = header.target_name(resolved).unwrap().to_owned();
            let last = header.target_len(resolved).unwrap() as u32 - 1;
            let ends = ends(f, tid);
            let mut spans = vec![span(0, last)];
            let step = (ends.len() / 2).max(1);
            for &(first, end) in ends.iter().step_by(step) {
                for at in [first.saturating_sub(1), end + 1] {
                    let at = at.min(last);
                    spans.push(span(at, at));
                }
                spans.push(span(first, end.min(last)));
            }
            for s in spans {
                let ours: Vec<RawBamRecord> =
                    reader.query_raw(resolved, s).unwrap().collect::<Result<_, _>>().unwrap();
                let theirs = samtools_region(f, &region_string(&name, s));
                let ours = payloads(&ours);
                assert!(
                    ours.len() == theirs.len() && ours.iter().zip(&theirs).all(|(a, b)| a == b),
                    "{} {}: {} records, samtools {}",
                    f.name(),
                    region_string(&name, s),
                    ours.len(),
                    theirs.len()
                );
                regions += 1;
            }
        }
    }
    assert!(regions > 50, "{regions} regions");
}

/// A region of a contig `len` long: anchored a base or two off a streamed
/// record's first or last base, or anywhere, a few widths wide.
#[hegel::composite]
fn region(tc: &TestCase, len: u32, ends: Vec<(u32, u32)>) -> RangeInclusive<Pos0> {
    let last = len - 1;
    let width = tc.draw(gs::sampled_from(vec![0u32, 1, 150, 5_000, len]));
    let start = if !ends.is_empty() && tc.draw(gs::weighted_booleans(0.7)) {
        let (pos, end) = ends[tc.draw(gs::integers::<usize>().max_value(ends.len() - 1))];
        let anchor = if tc.draw(gs::booleans()) { pos } else { end };
        anchor.saturating_add_signed(tc.draw(gs::integers::<i32>().min_value(-2).max_value(2)))
    } else {
        tc.draw(gs::integers::<u32>().max_value(last))
    }
    .min(last);
    span(start, start.saturating_add(width).min(last))
}

// r[verify cram.query_raw]
// r[verify cram.query_raw.errors]
// r[verify cram.stream.slice_decode]
// r[verify unified.query_raw+2]
/// Queries in a row on one reader, some read only part way and dropped,
/// some after a `fetch_into` of the same span:
/// each query read whole gives exactly the streamed records that overlap
/// it, byte for byte and in file order, whatever came before it; and, at
/// random, what `samtools view -b` writes for it.
#[hegel::test(test_cases = pinned::cases(40))]
fn raw_regions_are_the_streamed_records(tc: TestCase) {
    let fixtures = fixtures();
    let pick = if tc.draw(gs::weighted_booleans(0.3)) {
        fixtures.len() - 1
    } else {
        tc.draw(gs::integers::<usize>().max_value(fixtures.len() - 1))
    };
    let f = &fixtures[pick];
    let mut reader = f.open();
    let header = reader.header().clone();
    for _ in 0..tc.draw(gs::integers::<usize>().min_value(1).max_value(4)) {
        let tid = tc.draw(gs::integers::<u32>().max_value(header.target_count() as u32 - 1));
        let resolved = tid.resolve_tid(&header).unwrap();
        let len = header.target_len(resolved).unwrap() as u32;
        let s = tc.draw(region(len, ends(f, tid)));
        let want = filtered(&f.streamed, tid, s);
        if tc.draw(gs::weighted_booleans(0.3)) {
            reader.fetch_into(resolved, s, &mut RecordStore::new()).unwrap();
        }
        if tc.draw(gs::weighted_booleans(0.3)) {
            let k = tc.draw(gs::integers::<usize>().max_value(want.len()));
            let mut q = reader.query_raw(resolved, s).unwrap();
            let mut record = RawBamRecord::default();
            for _ in 0..k {
                assert!(q.read_record(&mut record).unwrap());
            }
            continue;
        }
        let got = query(&mut reader, tid, s);
        let name = header.target_name(resolved).unwrap();
        assert!(
            payloads(&got) == want,
            "{} {} (1-based): {} records, the stream {}",
            f.name(),
            region_string(name, s),
            got.len(),
            want.len()
        );
        if tc.draw(gs::weighted_booleans(0.25)) {
            let theirs = samtools_region(f, &region_string(name, s));
            assert!(payloads(&got) == theirs, "{} {}: samtools", f.name(), region_string(name, s));
        }
        tc.event_value("records", got.len() as f64);
    }
}

// r[verify record_store.seq_lossless]
// r[verify cram.stream.iupac+1]
// r[verify cram.record.no_ref_fill+1]
/// `fetch_into` keeps each record's SEQ as `query_raw` of the same region
/// yields it, code for code: one decode rebuilds both. Every reference of
/// every fixture, whole (IUPAC codes, `no_ref`'s `=`, embedded references,
/// multi-reference slices, v3.0 and v3.1, the generated shapes).
#[test]
fn fetched_seq_is_the_raw_querys() {
    let mut ambiguous = 0;
    for f in fixtures() {
        let mut reader = f.open();
        let header = reader.header().clone();
        for tid in 0..header.target_count() as u32 {
            let resolved = tid.resolve_tid(&header).unwrap();
            let s = span(0, header.target_len(resolved).unwrap() as u32 - 1);
            let mut store = RecordStore::new();
            reader.fetch_into(resolved, s, &mut store).unwrap();
            let fetched: Vec<(Vec<u8>, Vec<IupacBase>)> = store
                .indices()
                .map(|i| {
                    let r = store.record(i).unwrap();
                    (r.qname().to_vec(), r.seq().to_vec())
                })
                .collect();
            let queried: Vec<(Vec<u8>, Vec<IupacBase>)> = query(&mut reader, tid, s)
                .iter()
                .map(|r| (r.qname().to_vec(), r.iupac_bases().collect()))
                .collect();
            let name = header.target_name(resolved).unwrap();
            assert!(fetched == queried, "{} {name}", f.name());
            ambiguous += fetched
                .iter()
                .flat_map(|(_, seq)| seq)
                .filter(|b| b.to_base() == seqair_types::Base::Unknown && **b != IupacBase::N)
                .count();
        }
    }
    assert!(ambiguous > 0, "the fixtures hold `=` or ambiguity codes");
}

// r[verify cram.query_raw.htslib]
/// A multi-reference slice decoded for one reference needs no other
/// reference: a FASTA with only `chr1` answers a `chr1` query of
/// `fidelity_multiref.cram` with the records the full FASTA gives.
#[test]
fn a_multi_reference_slice_needs_only_the_queried_reference() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = dir.path().join("chr1.fa");
    let chr1 =
        run(Command::new("samtools").arg("faidx").arg(data("cram-stream/fidelity.fa")).arg("chr1"));
    std::fs::write(&fasta, chr1).unwrap();
    run(Command::new("samtools").arg("faidx").arg(&fasta));
    let cram = data("cram-stream/fidelity_multiref.cram");
    let full = fixtures().iter().find(|f| f.cram == cram).unwrap();
    let only_chr1 =
        Fixture { cram: cram.clone(), fasta, crai: full.crai.clone(), streamed: Vec::new() };
    let whole = span(0, 299);
    let got = query(&mut only_chr1.open(), 0, whole);
    let want = query(&mut full.open(), 0, whole);
    assert!(!want.is_empty());
    assert!(got == want, "{} records with chr1 only, {} with both", got.len(), want.len());
}

/// A contig and reads for generated CRAMs: `n` 100 bp reads on `chr1`, the
/// one at `huge` (if any) 1.5 Mbp long instead, too large for a BAM record.
struct Generated {
    cram: PathBuf,
    fasta: PathBuf,
    crai: PathBuf,
    _dir: tempfile::TempDir,
}

const HUGE_LEN: u32 = 1_500_000;
const GENERATED_LEN: u32 = 1_600_000;

fn generate(n: u32, huge: Option<u32>, options: &[&str]) -> Generated {
    generate_with(options, |reference| {
        (0..n)
            .map(|i| {
                let pos = i * 20;
                let len = if huge == Some(i) { HUGE_LEN } else { 100 };
                let seq = &reference[pos as usize..(pos + len) as usize];
                format!("r{i}\t0\tchr1\t{}\t60\t{len}M\t*\t0\t0\t{seq}\t*", pos + 1)
            })
            .collect()
    })
}

/// Reads of the shapes whose reference span differs from their length: a
/// spliced `20M3000N20M`, a `10M500D10M`, a 2 kb read, and a pair whose
/// unmapped mate is placed at the mapped one's position; `r<i>` at
/// `50 * i`, 100 bp `M` reads between them.
fn shapes(reference: &str) -> Vec<String> {
    let seq = |from: u32, len: u32| &reference[from as usize..(from + len) as usize];
    let mut lines = Vec::new();
    for i in 0..60u32 {
        let (pos, p) = (i * 50, i * 50 + 1);
        let line = |flag: u32, cigar: &str, mate: &str, seq: String| {
            format!("r{i}\t{flag}\tchr1\t{p}\t60\t{cigar}\t{mate}\t0\t{seq}\t*")
        };
        match i % 5 {
            0 => lines.push(line(0, "100M", "*\t0", seq(pos, 100).into())),
            1 => lines.push(line(
                0,
                "20M3000N20M",
                "*\t0",
                [seq(pos, 20), seq(pos + 3_020, 20)].concat(),
            )),
            2 => lines.push(line(
                0,
                "10M500D10M",
                "*\t0",
                [seq(pos, 10), seq(pos + 510, 10)].concat(),
            )),
            3 => lines.push(line(0, "2000M", "*\t0", seq(pos, 2_000).into())),
            _ => {
                let mate = format!("=\t{p}");
                lines.push(line(73, "100M", &mate, seq(pos, 100).into()));
                lines.push(
                    line(133, "*", &mate, seq(pos + 7, 100).into())
                        .replacen("\t60\t*", "\t0\t*", 1),
                );
            }
        }
    }
    lines
}

/// A CRAM of the reads `lines` builds from the reference, written with
/// samtools' `options`.
fn generate_with(options: &[&str], lines: impl FnOnce(&str) -> Vec<String>) -> Generated {
    let dir = tempfile::tempdir().unwrap();
    let fasta = dir.path().join("ref.fa");
    let mut bases =
        String::with_capacity(GENERATED_LEN as usize + GENERATED_LEN as usize / 60 + 16);
    bases.push_str(">chr1\n");
    let mut x: u32 = 0x1234_5678;
    for i in 0..GENERATED_LEN {
        x ^= x << 13;
        x ^= x >> 17;
        x ^= x << 5;
        bases.push(['A', 'C', 'G', 'T'][(x % 4) as usize]);
        if i % 60 == 59 {
            bases.push('\n');
        }
    }
    bases.push('\n');
    std::fs::write(&fasta, &bases).unwrap();
    run(Command::new("samtools").arg("faidx").arg(&fasta));
    let reference: String = bases.lines().skip(1).collect();

    let mut sam = format!("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:{GENERATED_LEN}\n");
    for line in lines(&reference) {
        writeln!(sam, "{line}").unwrap();
    }
    let sam_path = dir.path().join("in.sam");
    std::fs::write(&sam_path, sam).unwrap();
    let cram = dir.path().join("generated.cram");
    let mut cmd = Command::new("samtools");
    cmd.args(["view", "--no-PG", "-C", "-T"]).arg(&fasta);
    for option in options {
        cmd.args(["--output-fmt-option", option]);
    }
    run(cmd.arg("-o").arg(&cram).arg(&sam_path));
    run(Command::new("samtools").arg("index").arg(&cram));
    let crai = PathBuf::from(format!("{}.crai", cram.display()));
    Generated { cram, fasta, crai, _dir: dir }
}

fn open_generated(g: &Generated) -> IndexedReader {
    IndexedReader::open_with_reference_and_index(&g.cram, &g.fasta, &IndexPath::new(&g.crai))
        .unwrap()
}

/// The read index a generated name `r<i>` carries.
fn generated_index(qname: &[u8]) -> u32 {
    std::str::from_utf8(qname).unwrap().trim_start_matches('r').parse().unwrap()
}

/// A CRAM of 60 reads, one of them too large to be a BAM record.
fn with_huge_read() -> &'static (Generated, u32) {
    static FILE: OnceLock<(Generated, u32)> = OnceLock::new();
    FILE.get_or_init(|| {
        let huge = 23;
        (generate(60, Some(huge), &["seqs_per_slice=10", "slices_per_container=2"]), huge)
    })
}

// r[verify cram.query_raw.errors]
/// Every record `fetch_into` keeps comes out of `query_raw` in its place: as
/// itself, or, for the read too large for BAM, as an error that carries its
/// record number and that the query continues past.
#[hegel::test(test_cases = pinned::cases(20))]
fn a_record_too_large_for_bam_is_an_error_in_its_place(tc: TestCase) {
    let (file, huge) = with_huge_read();
    let start = tc.draw(gs::integers::<u32>().max_value(60 * 20 + 200));
    let width = tc.draw(gs::sampled_from(vec![0u32, 1, 150, 1_000, GENERATED_LEN]));
    let s = span(start, start.saturating_add(width).min(GENERATED_LEN - 1));

    let mut reader = open_generated(file);
    let tid = "chr1".resolve_tid(reader.header()).unwrap();
    let mut store = RecordStore::new();
    reader.fetch_into(tid, s, &mut store).unwrap();
    let want: Vec<Option<u32>> = store
        .records()
        .map(|r| generated_index(r.qname(&store).unwrap()))
        .map(|i| (i != *huge).then_some(i))
        .collect();

    let mut query = reader.query_raw(tid, s).unwrap();
    let mut record = RawBamRecord::default();
    let mut got = Vec::new();
    loop {
        let before = record.clone();
        match query.read_record(&mut record) {
            Ok(true) => got.push(Some(generated_index(record.qname()))),
            Ok(false) => break,
            Err(e) => {
                let ReaderError::CramStream {
                    source:
                        CramStreamError::RecordEncode {
                            record: number,
                            source: RecordEncodeError::RecordTooLarge { .. },
                        },
                } = e
                else {
                    panic!("{e:?}");
                };
                assert_eq!(number, u64::from(*huge), "the error names the record");
                assert!(record == before, "an error changed the record");
                got.push(None);
            }
        }
    }
    assert_eq!(got, want, "{}..={}", s.start.as_u32(), s.last.as_u32());
    assert!(!query.read_record(&mut record).unwrap(), "spent, not poisoned");
    drop(query);

    let iterated: Vec<Option<u32>> = reader
        .query_raw(tid, s)
        .unwrap()
        .map(|r| r.ok().map(|r| generated_index(r.qname())))
        .collect();
    assert_eq!(iterated, want, "the iterator");
    tc.event_value("errors", want.iter().filter(|w| w.is_none()).count() as f64);
}

/// Byte offset of each data container in `cram`, from its index.
fn container_offsets(crai: &Path) -> Vec<u64> {
    let text = run(Command::new("gzip").arg("-dc").arg(crai));
    let mut offsets: Vec<u64> = std::str::from_utf8(&text)
        .unwrap()
        .lines()
        .map(|l| l.split('\t').nth(3).unwrap().parse().unwrap())
        .collect();
    offsets.dedup();
    offsets
}

// r[verify cram.query_raw.errors]
/// A CRAM cut off inside a later container, and one with a corrupt byte in a
/// later container: the query yields the records before it, then an error
/// that is not about one record, and from then on only `Poisoned`; the
/// iterator stops after the error.
#[test]
fn a_truncated_or_corrupt_cram_poisons_its_raw_query() {
    let file = generate(300, None, &["seqs_per_slice=20"]);
    let bytes = std::fs::read(&file.cram).unwrap();
    let offsets = container_offsets(&file.crai);
    assert!(offsets.len() > 5, "{offsets:?}");
    let broken = offsets[offsets.len() / 2] as usize;
    let next = offsets[offsets.len() / 2 + 1] as usize;
    let whole = span(0, GENERATED_LEN - 1);

    let truncated = bytes[..broken + (next - broken) / 2].to_vec();
    let mut corrupt = bytes.clone();
    corrupt[next - 6] ^= 0xFF;
    for (what, contents) in [("truncated", truncated), ("corrupt", corrupt)] {
        std::fs::write(&file.cram, &contents).unwrap();
        let mut reader = open_generated(&file);
        let tid = "chr1".resolve_tid(reader.header()).unwrap();
        let mut query = reader.query_raw(tid, whole).unwrap();
        let mut record = RawBamRecord::default();
        let mut read = 0;
        let err = loop {
            match query.read_record(&mut record) {
                Ok(true) => read += 1,
                Ok(false) => panic!("{what}: ended cleanly after {read} records"),
                Err(e) => break e,
            }
        };
        assert!(read > 0 && read < 300, "{what}: {read}");
        assert!(
            !matches!(
                err,
                ReaderError::CramStream {
                    source: CramStreamError::InvalidRecord { .. }
                        | CramStreamError::RecordEncode { .. }
                }
            ),
            "{what}: {err:?}"
        );
        if what == "truncated" {
            assert!(
                matches!(
                    err,
                    ReaderError::CramStream {
                        source: CramStreamError::TruncatedContainerAt { offset }
                    } if offset == broken as u64
                ),
                "{err:?}"
            );
        }
        for _ in 0..2 {
            let again = query.read_record(&mut record);
            assert!(
                matches!(again, Err(ReaderError::CramStream { source: CramStreamError::Poisoned })),
                "{what}: {again:?}"
            );
        }
        drop(query);

        let mut items = reader.query_raw(tid, whole).unwrap();
        let oks = items.by_ref().take_while(Result::is_ok).count();
        assert_eq!(oks, read, "{what}");
        assert!(items.next().is_none(), "{what}: the iterator ends at the error");
    }
    std::fs::write(&file.cram, &bytes).unwrap();
}

// r[verify cram.stream.container_header+2]
/// A container of 512 two-record slices has a container header longer than
/// a fixed 1 KiB read: `fetch_into` and `query_raw` both read every record.
#[test]
fn a_container_with_hundreds_of_slices_is_read_whole() {
    let file = generate(1_024, None, &["seqs_per_slice=2", "slices_per_container=512"]);
    let offsets = container_offsets(&file.crai);
    assert!(offsets.len() <= 2, "{} containers", offsets.len());
    let mut reader = open_generated(&file);
    let tid = "chr1".resolve_tid(reader.header()).unwrap();
    let whole = span(0, GENERATED_LEN - 1);
    let mut store = RecordStore::new();
    assert_eq!(reader.fetch_into(tid, whole, &mut store).unwrap(), 1_024);
    let raw: Vec<RawBamRecord> =
        reader.query_raw(tid, whole).unwrap().collect::<Result<_, _>>().unwrap();
    assert_eq!(raw.len(), 1_024);
}
