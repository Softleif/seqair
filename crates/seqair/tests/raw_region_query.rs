//! `query_raw` on the fixture BAM and a bgzipped SAM of it, against two
//! oracles: what `samtools view <file> <region>` prints, and the records the
//! streaming readers yield, filtered to the region by a CIGAR walk of their
//! own.
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
use seqair::bam::{BgzfStream, RawBamRecord};
use seqair::reader::{IndexedReader, ReaderError, ResolveTid};
use seqair::sam::SamStreamReader;
use seqair::sam::format::format_sam_record;
use seqair_types::Pos0;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::OnceLock;

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/data").join(name)
}

fn run(cmd: &mut Command) -> Vec<u8> {
    let out = cmd.output().unwrap();
    assert!(out.status.success(), "{cmd:?}: {}", out.stderr.escape_ascii());
    out.stdout
}

struct Fixture {
    bam: PathBuf,
    sam: PathBuf,
    bam_records: Vec<RawBamRecord>,
    sam_records: Vec<RawBamRecord>,
    _dir: tempfile::TempDir,
}

/// The fixture BAM, the same records as bgzipped, tabix-indexed SAM, and
/// every record of each, streamed.
fn fixture() -> &'static Fixture {
    static FIXTURE: OnceLock<Fixture> = OnceLock::new();
    FIXTURE.get_or_init(|| {
        let dir = tempfile::tempdir().unwrap();
        let bam = data("test.bam");
        let sam = dir.path().join("test.sam.gz");
        run(Command::new("samtools")
            .args(["view", "--no-PG", "-h", "--output-fmt", "SAM,level=6", "-o"])
            .arg(&sam)
            .arg(&bam));
        run(Command::new("tabix").args(["-p", "sam"]).arg(&sam));
        let open = |path: &Path| BgzfStream::builder(std::fs::File::open(path).unwrap()).build();
        let bam_records = BamStreamReader::from_bgzf(open(&bam))
            .unwrap()
            .records()
            .collect::<Result<_, _>>()
            .unwrap();
        let sam_records =
            SamStreamReader::new(open(&sam)).unwrap().records().collect::<Result<_, _>>().unwrap();
        Fixture { bam, sam, bam_records, sam_records, _dir: dir }
    })
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
fn streamed(records: &[RawBamRecord], tid: u32, span: RangeInclusive<Pos0>) -> Vec<&[u8]> {
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

/// A region of a contig `len` long: anchored a base or two off a streamed
/// record's first or last base, or anywhere, a few widths wide.
#[hegel::composite]
fn region(tc: &TestCase, len: u32, ends: Vec<(u32, u32)>) -> RangeInclusive<Pos0> {
    let last = len - 1;
    let width = tc.draw(gs::sampled_from(vec![0u32, 1, 150, 5_000]));
    let start = if !ends.is_empty() && tc.draw(gs::weighted_booleans(0.7)) {
        let (pos, end) = ends[tc.draw(gs::integers::<usize>().max_value(ends.len() - 1))];
        let anchor = if tc.draw(gs::booleans()) { pos } else { end };
        anchor.saturating_add_signed(tc.draw(gs::integers::<i32>().min_value(-2).max_value(2)))
    } else {
        tc.draw(gs::integers::<u32>().max_value(last))
    }
    .min(last);
    let start = Pos0::new(start).unwrap();
    RangeInclusive {
        start,
        last: Pos0::new(start.as_u32().saturating_add(width).min(last)).unwrap(),
    }
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

// r[verify bam.reader.query_raw]
// r[verify bam.reader.query_raw.errors]
// r[verify sam.reader.query_raw]
// r[verify unified.query_raw]
// r[verify bam.reader.resume]
/// A few regions in a row on one BAM and one SAM reader (so the resume
/// shortcut is in play): each query gives exactly the streamed records that
/// overlap it, byte for byte and in file order, and, formatted as SAM, the
/// lines `samtools view` prints for it.
#[hegel::test(test_cases = pinned::cases(12))]
fn raw_regions_are_the_streamed_records_and_what_samtools_prints(tc: TestCase) {
    let f = fixture();
    let mut bam = IndexedReader::open(&f.bam).unwrap();
    let mut sam = IndexedReader::open(&f.sam).unwrap();
    let header = bam.header().clone();
    let n_targets = header.target_count();

    for _ in 0..tc.draw(gs::integers::<usize>().min_value(1).max_value(4)) {
        let tid = tc.draw(gs::integers::<u32>().max_value(n_targets as u32 - 1));
        let resolved = tid.resolve_tid(&header).unwrap();
        let name = header.target_name(resolved).unwrap().to_owned();
        let len = header.target_len(resolved).unwrap() as u32;
        let ends: Vec<(u32, u32)> = f
            .bam_records
            .iter()
            .filter(|r| r.tid().is_some_and(|t| t.as_u32() == tid))
            .filter_map(|r| r.pos().map(|p| (p.as_u32(), last_base(r, p.as_u32()))))
            .collect();
        let span = tc.draw(region(len, ends));

        let from_bam = query(&mut bam, tid, span);
        assert!(
            payloads(&from_bam) == streamed(&f.bam_records, tid, span),
            "BAM {name}:{}-{} (0-based) differs from the stream",
            span.start.as_u32(),
            span.last.as_u32()
        );
        let from_sam = query(&mut sam, tid, span);
        assert!(
            payloads(&from_sam) == streamed(&f.sam_records, tid, span),
            "SAM {name}:{}-{} (0-based) differs from the stream",
            span.start.as_u32(),
            span.last.as_u32()
        );

        let region = format!("{name}:{}-{}", span.start.as_u32() + 1, span.last.as_u32() + 1);
        let want = run(Command::new("samtools").arg("view").arg(&f.bam).arg(&region));
        for (format, records) in [("BAM", &from_bam), ("SAM", &from_sam)] {
            let mut text = Vec::new();
            for record in records.iter() {
                format_sam_record(record, &header, &mut text).unwrap();
                text.push(b'\n');
            }
            assert!(
                text == want,
                "{format} {region} differs from samtools:\n{}\nwant\n{}",
                text.escape_ascii(),
                want.escape_ascii()
            );
        }
        tc.event_value("records", from_bam.len() as f64);
    }
}

// r[verify unified.query_raw]
#[test]
fn a_cram_raw_query_is_unsupported() {
    let mut reader =
        IndexedReader::open_with_reference(&data("test_v30.cram"), &data("test.fasta.gz")).unwrap();
    let tid = reader.header().tid("chr19").unwrap();
    let span = (Pos0::new(6_105_700).unwrap()..=Pos0::new(6_105_800).unwrap()).into();
    let err = reader.query_raw(tid, span).unwrap_err();
    assert!(matches!(err, ReaderError::RawQueryUnsupported), "{err:?}");
}
