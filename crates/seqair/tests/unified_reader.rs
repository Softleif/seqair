//! Tests for the unified reader (auto-detection of BAM vs SAM vs CRAM)
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    reason = "test code"
)]
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    reason = "test code with known small values"
)]

use seqair::bam::{Pos0, RecordStore, RejectUnmapped};
use seqair::reader::{IndexedReader, Readers};
use std::path::Path;
use std::process::{Command, Stdio};

fn test_bam_path() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.bam"))
}

fn test_fasta_path() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.fasta.gz"))
}

fn test_cram_path() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test_v30.cram"))
}

const CONTIGS: &[(&str, u64, u64)] = &[
    ("chr19", 6_103_076, 6_143_229),
    ("2kb_3_Unmodified", 1, 2_018),
    ("bacteriophage_lambda_CpG", 1, 48_502),
];

fn create_sam_gz(dir: &Path) -> std::path::PathBuf {
    let sam_gz = dir.join("test.sam.gz");
    let status = Command::new("samtools")
        .args(["view", "-h", "--output-fmt", "SAM,level=6", "-o"])
        .arg(&sam_gz)
        .arg(test_bam_path())
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .expect("samtools not found");
    assert!(status.success(), "samtools view failed");
    let status = Command::new("tabix")
        .args(["-p", "sam"])
        .arg(&sam_gz)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .expect("tabix not found");
    assert!(status.success(), "tabix indexing failed");
    sam_gz
}

// r[verify unified.detect_format]
// r[verify unified.reader_enum]
// r[verify unified.detect_index]
#[test]
fn detects_bam_format() {
    let reader = IndexedReader::open(test_bam_path()).expect("should detect BAM");
    assert!(matches!(reader, IndexedReader::Bam(_)));
}

// r[verify unified.detect_format]
// r[verify unified.reader_enum]
// r[verify unified.detect_index]
#[test]
fn detects_sam_format() {
    let dir = tempfile::tempdir().unwrap();
    let sam_gz = create_sam_gz(dir.path());
    let reader = IndexedReader::open(&sam_gz).expect("should detect SAM");
    assert!(matches!(reader, IndexedReader::Sam(_)));
}

// r[verify unified.detect_error]
#[test]
fn rejects_uncompressed_sam() {
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("test.sam");
    std::fs::write(&sam, b"@HD\tVN:1.6\n").unwrap();
    let err = IndexedReader::open(&sam).unwrap_err();
    let msg = err.to_string();
    assert!(msg.contains("uncompressed SAM"), "error should mention uncompressed SAM: {msg}");
}

// r[verify unified.detect_error]
#[test]
fn rejects_unknown_format() {
    let dir = tempfile::tempdir().unwrap();
    let garbage = dir.path().join("garbage.xyz");
    std::fs::write(&garbage, b"NOTAFORMAT12345678").unwrap();
    let err = IndexedReader::open(&garbage).unwrap_err();
    let msg = err.to_string();
    assert!(msg.contains("unrecognized"), "error should mention unrecognized: {msg}");
}

// r[verify unified.fetch_equivalence]
// r[verify unified.reader_api]
// r[verify index.shared_query]
// r[verify index.shared_types]
#[test]
fn bam_and_sam_produce_same_records() {
    let dir = tempfile::tempdir().unwrap();
    let sam_gz = create_sam_gz(dir.path());

    let mut bam_reader = IndexedReader::open(test_bam_path()).expect("bam");
    let mut sam_reader = IndexedReader::open(&sam_gz).expect("sam");

    for &(contig, start, end) in CONTIGS {
        let bam_tid = bam_reader.header().tid(contig).unwrap();
        let sam_tid = sam_reader.header().tid(contig).unwrap();

        let mut bam_store = RecordStore::new();
        bam_reader
            .fetch_into(
                bam_tid,
                Pos0::new(start as u32).unwrap(),
                Pos0::new(end as u32).unwrap(),
                &mut bam_store,
            )
            .unwrap();

        let mut sam_store = RecordStore::new();
        sam_reader
            .fetch_into(
                sam_tid,
                Pos0::new(start as u32).unwrap(),
                Pos0::new(end as u32).unwrap(),
                &mut sam_store,
            )
            .unwrap();

        assert_eq!(
            bam_store.len(),
            sam_store.len(),
            "{contig}: record count mismatch bam={} sam={}",
            bam_store.len(),
            sam_store.len()
        );

        for i in 0..bam_store.len() as u32 {
            assert_eq!(bam_store.record(i).pos, sam_store.record(i).pos, "{contig} rec {i}: pos");
            assert_eq!(
                bam_store.record(i).flags,
                sam_store.record(i).flags,
                "{contig} rec {i}: flags"
            );
            assert_eq!(bam_store.seq(i), sam_store.seq(i), "{contig} rec {i}: seq");
            assert_eq!(bam_store.qual(i), sam_store.qual(i), "{contig} rec {i}: qual");
        }
    }
}

// r[verify interval.overlap_test]
// r[verify interval.end_pos_inclusive]
// r[verify interval.inclusive_ends]
/// The one-base edges of a query are where the three readers used to disagree:
/// SAM tested a half-open `[start, end)` against an exclusive reading of
/// `end_pos`, so it dropped a record overlapping the region by exactly its last
/// base while BAM and CRAM kept it. Query each edge of a known record and
/// require all three formats to answer identically.
#[test]
fn every_format_agrees_on_the_edges_of_a_query() {
    let dir = tempfile::tempdir().unwrap();
    let sam_gz = create_sam_gz(dir.path());

    // Take a real record's span from the BAM, then probe its two edges.
    let mut bam_reader = IndexedReader::open(test_bam_path()).expect("bam");
    let tid = bam_reader.header().tid("chr19").unwrap();
    let mut store = RecordStore::new();
    bam_reader
        .fetch_into(tid, Pos0::new(6_103_076).unwrap(), Pos0::new(6_104_000).unwrap(), &mut store)
        .unwrap();
    let probe = store.record(0);
    let (first, last) = (probe.pos, probe.end_pos);
    assert!(last > first, "need a record spanning more than one base");
    let qname = store.qname(0).to_vec();

    /// A named way to open the same data in one of the three formats.
    type OpenReader<'a> = (&'a str, Box<dyn Fn() -> IndexedReader + 'a>);

    let readers: [OpenReader<'_>; 3] = [
        ("bam", Box::new(|| IndexedReader::open(test_bam_path()).expect("bam"))),
        (
            "sam",
            Box::new({
                let sam_gz = sam_gz.clone();
                move || IndexedReader::open(&sam_gz).expect("sam")
            }),
        ),
        (
            "cram",
            Box::new(|| {
                IndexedReader::open_with_reference(test_cram_path(), test_fasta_path())
                    .expect("cram")
            }),
        ),
    ];

    // A single-position query on the record's *last* base, and one on its first.
    for (label, edge) in [("last base", last), ("first base", first)] {
        let mut found = Vec::new();
        for (format, open) in &readers {
            let mut reader = open();
            let tid = reader.header().tid("chr19").unwrap();
            let mut store = RecordStore::new();
            reader.fetch_into(tid, edge, edge, &mut store).unwrap();
            let hit = (0..store.len() as u32).any(|i| store.qname(i) == qname.as_slice());
            assert!(
                hit,
                "{format}: the record covering its own {label} ({}) must overlap a query there",
                edge.as_u64()
            );
            found.push((*format, store.len()));
        }
        let (_, expected) = found[0];
        for &(format, len) in &found {
            assert_eq!(
                len, expected,
                "{format}: {label} query returned {len} records, bam returned {expected}"
            );
        }
    }

    // One base past the record's end, it must be gone from every format.
    let past = Pos0::new(last.as_u64() as u32 + 1).unwrap();
    for (format, open) in &readers {
        let mut reader = open();
        let tid = reader.header().tid("chr19").unwrap();
        let mut store = RecordStore::new();
        reader.fetch_into(tid, past, past, &mut store).unwrap();
        // Other records may cover that position; this one must not be reported
        // as covering it unless it genuinely does.
        for i in 0..store.len() as u32 {
            if store.qname(i) == qname.as_slice() {
                assert!(
                    store.record(i).end_pos >= past,
                    "{format}: record returned for a position it does not cover"
                );
            }
        }
    }
}

// r[verify unified.fork_bam]
// r[verify unified.fork_sam]
#[test]
fn fork_works_for_both_formats() {
    let dir = tempfile::tempdir().unwrap();
    let sam_gz = create_sam_gz(dir.path());

    let bam_reader = IndexedReader::open(test_bam_path()).expect("bam");
    let sam_reader = IndexedReader::open(&sam_gz).expect("sam");

    let mut bam_fork = bam_reader.fork().expect("bam fork");
    let mut sam_fork = sam_reader.fork().expect("sam fork");

    let tid = bam_fork.header().tid("chr19").unwrap();
    let mut store = RecordStore::new();
    bam_fork
        .fetch_into(tid, Pos0::new(6_105_700).unwrap(), Pos0::new(6_105_800).unwrap(), &mut store)
        .unwrap();
    assert!(!store.is_empty(), "bam fork should fetch records");

    let sam_tid = sam_fork.header().tid("chr19").unwrap();
    sam_fork
        .fetch_into(
            sam_tid,
            Pos0::new(6_105_700).unwrap(),
            Pos0::new(6_105_800).unwrap(),
            &mut store,
        )
        .unwrap();
    assert!(!store.is_empty(), "sam fork should fetch records");
}

// ── CRAM format detection ────────────────────────────────────────────

// r[verify unified.detect_format]
// r[verify unified.readers_backward_compat]
#[test]
fn indexed_reader_open_rejects_cram_without_fasta() {
    let err = IndexedReader::open(test_cram_path()).unwrap_err();
    let msg = err.to_string();
    assert!(
        msg.contains("CRAM") && msg.contains("reference"),
        "error should mention CRAM and reference: {msg}"
    );
}

// ── Readers struct tests ─────────────────────────────────────────────

// r[verify unified.readers_struct]
// r[verify unified.readers_open]
// r[verify unified.readers_accessors]
#[test]
fn readers_open_bam() {
    let readers =
        Readers::open_customized(test_bam_path(), test_fasta_path(), RejectUnmapped).unwrap();
    assert!(readers.header().target_count() > 0);
    assert!(matches!(readers.alignment(), &IndexedReader::Bam(_)));
}

// r[verify unified.readers_open]
#[test]
fn readers_open_sam() {
    let dir = tempfile::tempdir().unwrap();
    let sam_gz = create_sam_gz(dir.path());
    let readers = Readers::open_customized(&sam_gz, test_fasta_path(), RejectUnmapped).unwrap();
    assert!(matches!(readers.alignment(), &IndexedReader::Sam(_)));
}

// r[verify unified.readers_open]
// r[verify unified.detect_format]
// r[verify unified.fork_cram]
#[test]
fn readers_open_cram() {
    let readers =
        Readers::open_customized(test_cram_path(), test_fasta_path(), RejectUnmapped).unwrap();
    assert!(matches!(readers.alignment(), &IndexedReader::Cram(_)));
    assert!(readers.header().target_count() > 0);
}

// r[verify unified.readers_fork]
#[test]
fn readers_fork_cram() {
    let readers =
        Readers::open_customized(test_cram_path(), test_fasta_path(), RejectUnmapped).unwrap();
    let mut forked = readers.fork().unwrap();
    let tid = forked.header().tid("chr19").unwrap();
    let mut store = RecordStore::new();
    forked
        .fetch_into(tid, Pos0::new(6_105_700).unwrap(), Pos0::new(6_105_800).unwrap(), &mut store)
        .unwrap();
    assert!(!store.is_empty(), "cram fork should fetch records");
}

// r[verify unified.fetch_equivalence]
#[test]
fn bam_and_cram_produce_same_records() {
    let mut bam =
        Readers::open_customized(test_bam_path(), test_fasta_path(), RejectUnmapped).unwrap();
    let mut cram =
        Readers::open_customized(test_cram_path(), test_fasta_path(), RejectUnmapped).unwrap();

    for &(contig, start, end) in CONTIGS {
        let bam_tid = bam.header().tid(contig).unwrap();
        let cram_tid = cram.header().tid(contig).unwrap();

        let mut bam_store = RecordStore::new();
        bam.fetch_into(
            bam_tid,
            Pos0::new(start as u32).unwrap(),
            Pos0::new(end as u32).unwrap(),
            &mut bam_store,
        )
        .unwrap();

        let mut cram_store = RecordStore::new();
        cram.fetch_into(
            cram_tid,
            Pos0::new(start as u32).unwrap(),
            Pos0::new(end as u32).unwrap(),
            &mut cram_store,
        )
        .unwrap();

        assert_eq!(
            bam_store.len(),
            cram_store.len(),
            "{contig}: record count mismatch bam={} cram={}",
            bam_store.len(),
            cram_store.len()
        );

        for i in 0..bam_store.len() as u32 {
            assert_eq!(bam_store.record(i).pos, cram_store.record(i).pos, "{contig} rec {i}: pos");
            assert_eq!(
                bam_store.record(i).flags,
                cram_store.record(i).flags,
                "{contig} rec {i}: flags"
            );
            assert_eq!(
                bam_store.record(i).mapq,
                cram_store.record(i).mapq,
                "{contig} rec {i}: mapq"
            );
        }
    }
}

// r[verify unified.fetch_equivalence]
#[test]
fn all_three_formats_produce_same_records() {
    let dir = tempfile::tempdir().unwrap();
    let sam_gz = create_sam_gz(dir.path());

    let mut bam =
        Readers::open_customized(test_bam_path(), test_fasta_path(), RejectUnmapped).unwrap();
    let mut sam = Readers::open_customized(&sam_gz, test_fasta_path(), RejectUnmapped).unwrap();
    let mut cram =
        Readers::open_customized(test_cram_path(), test_fasta_path(), RejectUnmapped).unwrap();

    let contig = "chr19";
    let start = 6_105_700u64;
    let end = 6_106_200u64;

    let bam_tid = bam.header().tid(contig).unwrap();
    let sam_tid = sam.header().tid(contig).unwrap();
    let cram_tid = cram.header().tid(contig).unwrap();

    let mut bam_store = RecordStore::new();
    let mut sam_store = RecordStore::new();
    let mut cram_store = RecordStore::new();

    bam.fetch_into(
        bam_tid,
        Pos0::new(start as u32).unwrap(),
        Pos0::new(end as u32).unwrap(),
        &mut bam_store,
    )
    .unwrap();
    sam.fetch_into(
        sam_tid,
        Pos0::new(start as u32).unwrap(),
        Pos0::new(end as u32).unwrap(),
        &mut sam_store,
    )
    .unwrap();
    cram.fetch_into(
        cram_tid,
        Pos0::new(start as u32).unwrap(),
        Pos0::new(end as u32).unwrap(),
        &mut cram_store,
    )
    .unwrap();

    assert_eq!(bam_store.len(), sam_store.len(), "BAM vs SAM count");
    assert_eq!(bam_store.len(), cram_store.len(), "BAM vs CRAM count");

    for i in 0..bam_store.len() as u32 {
        let bp = bam_store.record(i).pos;
        let sp = sam_store.record(i).pos;
        let cp = cram_store.record(i).pos;
        assert_eq!(bp, sp, "BAM vs SAM pos at {i}");
        assert_eq!(bp, cp, "BAM vs CRAM pos at {i}");
    }
}
