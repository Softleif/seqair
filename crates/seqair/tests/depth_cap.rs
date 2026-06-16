//! Load-time depth-cap cross-validation against rust-htslib's `bam_plp_auto`.
//!
//! Two properties:
//! 1. **Transparency** — a cap larger than the deepest column produces a pileup
//!    byte-identical to htslib (the cap never fires, so nothing changes).
//! 2. **Bounding** — a cap smaller than the true depth holds every column to
//!    `<= cap`, while htslib (and an uncapped fetch) see the full depth.
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

use std::fmt::Write as _;
use std::num::NonZeroU32;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use rust_htslib::bam::{self, FetchDefinition, Read as _};
use seqair::bam::{DepthCap, IndexedBamReader, PileupEngine, Pos0, RecordStore};

/// Write `sam` to a temp dir, sort + index it into a BAM, return the BAM path.
fn write_and_index(dir: &Path, sam: &str) -> PathBuf {
    let sam_path = dir.join("in.sam");
    std::fs::write(&sam_path, sam).unwrap();
    let bam_path = dir.join("test.bam");
    let status = Command::new("samtools")
        .args(["sort", "-o"])
        .arg(&bam_path)
        .arg(&sam_path)
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
    assert!(status.success());
    bam_path
}

/// Build a SAM where `reads_per_start` reads begin at each 1-based start in
/// `1..=n_starts`, each spanning `read_len` reference bases (`<read_len>M`).
/// This produces a coverage ramp whose peak exceeds `reads_per_start`.
fn ramp_sam(
    contig: &str,
    contig_len: u32,
    n_starts: u32,
    reads_per_start: u32,
    read_len: u32,
) -> String {
    let mut sam = String::new();
    writeln!(sam, "@HD\tVN:1.6\tSO:coordinate").unwrap();
    writeln!(sam, "@SQ\tSN:{contig}\tLN:{contig_len}").unwrap();
    let seq: String = "A".repeat(read_len as usize);
    let qual: String = "I".repeat(read_len as usize);
    for start in 1..=n_starts {
        for r in 0..reads_per_start {
            writeln!(
                sam,
                "r{start}_{r}\t0\t{contig}\t{start}\t60\t{read_len}M\t*\t0\t0\t{seq}\t{qual}"
            )
            .unwrap();
        }
    }
    sam
}

/// (pos, depth) per column from rust-htslib's pileup.
fn htslib_depths(bam_path: &Path, contig: &[u8], start: u64, end: u64) -> Vec<(u32, u32)> {
    let mut reader = bam::IndexedReader::from_path(bam_path).expect("htslib open");
    let tid = reader.header().tid(contig).expect("tid");
    reader
        .fetch(FetchDefinition::Region(tid as i32, start as i64, end as i64))
        .expect("htslib fetch");
    let mut out = Vec::new();
    for p in reader.pileup() {
        let p = p.expect("htslib pileup");
        let pos = p.pos();
        if u64::from(pos) < start || u64::from(pos) >= end {
            continue;
        }
        out.push((pos, p.depth()));
    }
    out
}

/// (pos, depth) per column from seqair, fetching with an optional depth cap.
fn seqair_depths(
    bam_path: &Path,
    contig: &str,
    start: u32,
    end: u32,
    cap: Option<NonZeroU32>,
) -> Vec<(u32, u32)> {
    let mut reader = IndexedBamReader::open(bam_path).expect("seqair open");
    let tid = reader.header().tid(contig).expect("contig not found");
    let mut store = RecordStore::new();
    let mut customize = DepthCap::with_cap((), cap);
    reader
        .fetch_into_customized(
            tid,
            Pos0::new(start).unwrap(),
            Pos0::new(end).unwrap(),
            &mut store,
            &mut customize,
        )
        .expect("fetch");

    let mut engine = PileupEngine::new(store, Pos0::new(start).unwrap(), Pos0::new(end).unwrap());
    let mut out = Vec::new();
    while let Some(col) = engine.pileups() {
        out.push((u32::from(col.pos()), col.depth() as u32));
    }
    out
}

/// A cap larger than the deepest column must leave the pileup identical to
/// htslib — the cap never fires.
#[test]
fn cap_above_depth_matches_htslib() {
    let dir = tempfile::tempdir().unwrap();
    // Peak depth = 5 reads/start × up to 20 overlapping starts = 100.
    let bam = write_and_index(dir.path(), &ramp_sam("c1", 300, 50, 5, 20));

    let hts = htslib_depths(&bam, b"c1", 0, 300);
    let ours = seqair_depths(&bam, "c1", 0, 300, NonZeroU32::new(10_000));

    assert!(!hts.is_empty(), "expected covered columns");
    assert_eq!(ours, hts, "non-binding cap must reproduce htslib exactly");
    // Sanity: the data really is deeper than any cap we test below.
    assert!(hts.iter().map(|&(_, d)| d).max().unwrap() >= 50);
}

/// A cap below the true depth must bound every column to `<= cap`, while the
/// uncapped fetch (and htslib) still see the full depth.
#[test]
fn cap_below_depth_bounds_every_column() {
    let dir = tempfile::tempdir().unwrap();
    let bam = write_and_index(dir.path(), &ramp_sam("c1", 300, 50, 5, 20));

    let cap = 7;
    let capped = seqair_depths(&bam, "c1", 0, 300, NonZeroU32::new(cap));
    let uncapped = seqair_depths(&bam, "c1", 0, 300, None);
    let hts = htslib_depths(&bam, b"c1", 0, 300);

    assert_eq!(uncapped, hts, "uncapped fetch must match htslib");
    assert!(!capped.is_empty());
    for &(pos, depth) in &capped {
        assert!(depth <= cap, "pos {pos}: depth {depth} exceeds cap {cap}");
    }
    // The cap actually fired: somewhere the uncapped depth was above it.
    assert!(uncapped.iter().any(|&(_, d)| d > cap), "test data not deep enough to bind");
    // Capping never invents coverage — every kept column also exists uncapped.
    let uncapped_positions: std::collections::HashSet<u32> =
        uncapped.iter().map(|&(p, _)| p).collect();
    for &(pos, _) in &capped {
        assert!(uncapped_positions.contains(&pos), "capped produced phantom column {pos}");
    }
    // Greedy streaming capping keeps the earliest reads, so the tail of a
    // saturated ramp may be truncated. This only happens where depth exceeds
    // the cap, so the deeply-covered core is always retained.
    assert!(capped.iter().any(|&(_, d)| d == cap), "core columns should reach the cap");
}

/// Many reads stacked at a single position (the decoy/EBV hotspot shape):
/// the store must hold exactly `cap` reads, not the full pile.
#[test]
fn cap_bounds_single_position_pileup() {
    let dir = tempfile::tempdir().unwrap();
    // 500 reads all starting at position 1, each 30M.
    let bam = write_and_index(dir.path(), &ramp_sam("c1", 100, 1, 500, 30));

    let cap = 25;
    let capped = seqair_depths(&bam, "c1", 0, 100, NonZeroU32::new(cap));
    assert!(!capped.is_empty());
    for &(pos, depth) in &capped {
        assert!(depth <= cap, "pos {pos}: depth {depth} exceeds cap {cap}");
    }
    // The covered region (positions 0..30) is held at the cap, not 500.
    assert_eq!(capped.iter().map(|&(_, d)| d).max().unwrap(), cap);
}
