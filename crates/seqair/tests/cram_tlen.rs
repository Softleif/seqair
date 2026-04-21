//! TLEN validation: decode htslib's 30 TLEN test files and compare
//! every record field (pos, flags, qname, TLEN) against the expected SAM.
//!
//! htslib's tlen/ directory contains CRAM/SAM pairs testing all combinations
//! of read pair orientations, start/end matching, and ordering.
//!
//! The CRAMs in this directory use attached (downstream) mates. TLEN is
//! validated both via the BAM path (SAM -> sorted BAM -> seqair) and the
//! CRAM path (pre-built CRAM -> seqair), covering both storage formats.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "test code"
)]
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    reason = "test code with known small values"
)]

use seqair::bam::{IndexedBamReader, Pos0, RecordStore};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn tlen_dir() -> PathBuf {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/htslib/tlen/")).to_path_buf()
}

/// Parse a SAM file into a list of (qname, flags, pos, tlen) tuples.
/// Only includes mapped records (flag & 0x4 == 0).
fn parse_expected_sam(sam_path: &Path) -> Vec<(Vec<u8>, u16, i64, i32)> {
    let content = std::fs::read_to_string(sam_path).expect("read SAM");
    let mut records = Vec::new();
    for line in content.lines() {
        if line.starts_with('@') || line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        assert!(fields.len() >= 11, "SAM line has too few fields: {line}");
        let qname = fields[0].as_bytes().to_vec();
        let flags: u16 = fields[1].parse().unwrap();
        if flags & 0x4 != 0 {
            continue;
        }
        let pos: i64 = fields[3].parse::<i64>().unwrap() - 1; // 1-based -> 0-based
        let tlen: i32 = fields[8].parse().unwrap();
        records.push((qname, flags, pos, tlen));
    }
    records
}

/// Convert SAM to sorted, indexed BAM.
fn sam_to_indexed_bam(dir: &Path, sam_path: &Path) -> PathBuf {
    let bam_path = dir.join("test.bam");
    let status = Command::new("samtools")
        .args(["sort", "-o"])
        .arg(&bam_path)
        .arg(sam_path)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .expect("samtools not found");
    assert!(status.success(), "samtools sort failed");
    let status = Command::new("samtools")
        .arg("index")
        .arg(&bam_path)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .expect("samtools index failed");
    assert!(status.success(), "samtools index failed");
    bam_path
}

/// Validate TLEN for a single test case via the BAM path.
fn assert_tlen_via_bam(name: &str) {
    let dir = tlen_dir();
    let sam_path = dir.join(format!("{name}.sam"));
    assert!(sam_path.exists(), "SAM not found: {}", sam_path.display());

    let mut expected = parse_expected_sam(&sam_path);
    assert!(!expected.is_empty(), "{name}: no mapped records");
    // Sort by (pos, qname, flags) to match BAM sort order
    expected.sort_by(|a, b| a.2.cmp(&b.2).then(a.0.cmp(&b.0)).then(a.1.cmp(&b.1)));

    let tmpdir = tempfile::tempdir().unwrap();
    let bam_path = sam_to_indexed_bam(tmpdir.path(), &sam_path);

    let mut reader = IndexedBamReader::open(&bam_path).expect("open BAM");
    let tid = reader.header().tid("ref").expect("contig 'ref'");
    let mut store = RecordStore::new();
    reader
        .fetch_into(tid, Pos0::new(0).unwrap(), Pos0::new(20).unwrap(), &mut store)
        .expect("fetch");

    assert_eq!(
        store.len(),
        expected.len(),
        "{name}: record count seqair={} expected={}",
        store.len(),
        expected.len()
    );

    // Extract seqair records and sort the same way
    let mut actual: Vec<_> = (0..store.len() as u32)
        .map(|idx| {
            let r = store.record(idx);
            (store.qname(idx).to_vec(), r.flags.raw(), r.pos.as_i64(), r.template_len)
        })
        .collect();
    actual.sort_by(|a, b| a.2.cmp(&b.2).then(a.0.cmp(&b.0)).then(a.1.cmp(&b.1)));

    for (i, ((sq, sf, sp, st), (eq, ef, ep, et))) in actual.iter().zip(&expected).enumerate() {
        let qn = String::from_utf8_lossy(eq);
        assert_eq!(sq.as_slice(), eq.as_slice(), "{name}[{i}]: qname");
        assert_eq!(*sf, *ef, "{name}[{i}] ({qn}): flags");
        assert_eq!(*sp, *ep, "{name}[{i}] ({qn}): pos");
        assert_eq!(*st, *et, "{name}[{i}] ({qn}): tlen seqair={st} expected={et}");
    }
}

/// Validate that the CRAM path preserves all fields including TLEN.
fn assert_cram_fields_with_tlen(name: &str) {
    let dir = tlen_dir();
    let cram_path = dir.join(format!("{name}.cram"));
    let sam_path = dir.join(format!("{name}.sam"));
    assert!(cram_path.exists(), "CRAM not found: {}", cram_path.display());

    let mut expected = parse_expected_sam(&sam_path);
    expected.sort_by(|a, b| a.2.cmp(&b.2).then(a.0.cmp(&b.0)).then(a.1.cmp(&b.1)));

    // Create minimal ref FASTA + index for Readers API
    let tmpdir = tempfile::tempdir().unwrap();
    let ref_fasta = tmpdir.path().join("ref.fa");
    let ref_fai = tmpdir.path().join("ref.fa.fai");
    std::fs::write(&ref_fasta, ">ref\nAAAAACCCCCGGGGGTTTTT\n").unwrap();
    std::fs::write(&ref_fai, "ref\t20\t5\t20\t21\n").unwrap();

    // Copy CRAM and create index
    let cram_copy = tmpdir.path().join(format!("{name}.cram"));
    std::fs::copy(&cram_path, &cram_copy).unwrap();
    let status = Command::new("samtools")
        .arg("index")
        .arg(&cram_copy)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .expect("samtools index");
    assert!(status.success(), "samtools index failed for {name}");

    let mut readers = seqair::reader::Readers::open(&cram_copy, &ref_fasta)
        .unwrap_or_else(|e| panic!("{name}: open failed: {e:#}"));
    let tid = readers.header().tid("ref").expect("contig 'ref'");
    let mut store = RecordStore::new();
    readers
        .fetch_into(tid, Pos0::new(0).unwrap(), Pos0::new(20).unwrap(), &mut store)
        .unwrap_or_else(|e| panic!("{name}: fetch failed: {e:#}"));

    assert_eq!(store.len(), expected.len(), "{name}: record count");

    let mut actual: Vec<_> = (0..store.len() as u32)
        .map(|idx| {
            let r = store.record(idx);
            (store.qname(idx).to_vec(), r.flags.raw(), r.pos.as_i64(), r.template_len)
        })
        .collect();
    actual.sort_by(|a, b| a.2.cmp(&b.2).then(a.0.cmp(&b.0)).then(a.1.cmp(&b.1)));

    for (i, ((sq, sf, sp, st), (eq, ef, ep, et))) in actual.iter().zip(&expected).enumerate() {
        let qn = String::from_utf8_lossy(eq);
        assert_eq!(sq.as_slice(), eq.as_slice(), "{name}[{i}]: qname");
        assert_eq!(*sf, *ef, "{name}[{i}] ({qn}): flags");
        assert_eq!(*sp, *ep, "{name}[{i}] ({qn}): pos");
        assert_eq!(*st, *et, "{name}[{i}] ({qn}): cram tlen seqair={st} expected={et}");
    }
}

// --- BAM path: full TLEN validation (30 tests) ---

#[test]
fn tlen_a7() {
    assert_tlen_via_bam("a7");
}
#[test]
fn tlen_a7b() {
    assert_tlen_via_bam("a7b");
}
#[test]
fn tlen_a8() {
    assert_tlen_via_bam("a8");
}
#[test]
fn tlen_a8b() {
    assert_tlen_via_bam("a8b");
}
#[test]
fn tlen_a9() {
    assert_tlen_via_bam("a9");
}
#[test]
fn tlen_a9b() {
    assert_tlen_via_bam("a9b");
}
#[test]
fn tlen_b7() {
    assert_tlen_via_bam("b7");
}
#[test]
fn tlen_b7b() {
    assert_tlen_via_bam("b7b");
}
#[test]
fn tlen_b8() {
    assert_tlen_via_bam("b8");
}
#[test]
fn tlen_b8b() {
    assert_tlen_via_bam("b8b");
}
#[test]
fn tlen_c7() {
    assert_tlen_via_bam("c7");
}
#[test]
fn tlen_c7b() {
    assert_tlen_via_bam("c7b");
}
#[test]
fn tlen_c8() {
    assert_tlen_via_bam("c8");
}
#[test]
fn tlen_c8b() {
    assert_tlen_via_bam("c8b");
}
#[test]
fn tlen_d7() {
    assert_tlen_via_bam("d7");
}
#[test]
fn tlen_d7b() {
    assert_tlen_via_bam("d7b");
}
#[test]
fn tlen_d4() {
    assert_tlen_via_bam("d4");
}
#[test]
fn tlen_d4b() {
    assert_tlen_via_bam("d4b");
}
#[test]
fn tlen_d4c() {
    assert_tlen_via_bam("d4c");
}
#[test]
fn tlen_d4d() {
    assert_tlen_via_bam("d4d");
}
#[test]
fn tlen_d4e() {
    assert_tlen_via_bam("d4e");
}
#[test]
fn tlen_d4f() {
    assert_tlen_via_bam("d4f");
}
#[test]
fn tlen_d5() {
    assert_tlen_via_bam("d5");
}
#[test]
fn tlen_d5b() {
    assert_tlen_via_bam("d5b");
}
#[test]
fn tlen_d5c() {
    assert_tlen_via_bam("d5c");
}
#[test]
fn tlen_d5d() {
    assert_tlen_via_bam("d5d");
}
#[test]
fn tlen_d5e() {
    assert_tlen_via_bam("d5e");
}
#[test]
fn tlen_d5f() {
    assert_tlen_via_bam("d5f");
}
#[test]
fn tlen_a4() {
    assert_tlen_via_bam("a4");
}
#[test]
fn tlen_a5() {
    assert_tlen_via_bam("a5");
}

// --- CRAM path: full TLEN validation via attached mate reconstruction ---

// r[verify cram.record.mate_tlen_reconstruction]
#[test]
fn cram_fields_a7() {
    assert_cram_fields_with_tlen("a7");
}
#[test]
fn cram_fields_b7() {
    assert_cram_fields_with_tlen("b7");
}
#[test]
fn cram_fields_c7() {
    assert_cram_fields_with_tlen("c7");
}
#[test]
fn cram_fields_d7() {
    assert_cram_fields_with_tlen("d7");
}
#[test]
fn cram_fields_d4() {
    assert_cram_fields_with_tlen("d4");
}
#[test]
fn cram_fields_a4() {
    assert_cram_fields_with_tlen("a4");
}
