//! Comparison tests: Pileup engine.
//! Builds pileups over the same region using both htslib and seqair,
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::indexing_slicing)]

use rust_htslib::bam::{self, FetchDefinition, Read as _};
use seqair::bam::{Pos, Zero};
use std::path::Path;

fn test_bam_path() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.bam"))
}

const TEST_REGION: &str = "chr19";
const TEST_START: u64 = 6_105_700;
const TEST_END: u64 = 6_105_800;

struct HtsPileupColumn {
    pos: u32,
    depth: u32,
    /// (qpos, flags) for each alignment that has a qpos (non-deletion)
    alignments: Vec<(usize, u16)>,
}

fn fetch_htslib_pileup() -> Vec<HtsPileupColumn> {
    let bam_path = test_bam_path();
    let mut reader = bam::IndexedReader::from_path(bam_path).expect("htslib open");
    let tid = reader.header().tid(TEST_REGION.as_bytes()).expect("tid");
    reader
        .fetch(FetchDefinition::Region(tid as i32, TEST_START as i64, TEST_END as i64))
        .expect("htslib fetch");

    let mut columns = Vec::new();
    let pileup = reader.pileup();
    for p in pileup {
        let p = p.expect("htslib pileup");
        let pos = p.pos();
        if (pos as u64) < TEST_START || (pos as u64) > TEST_END {
            continue;
        }
        let mut alignments = Vec::new();
        for a in p.alignments() {
            if let Some(qpos) = a.qpos() {
                alignments.push((qpos, a.record().flags()));
            }
        }
        alignments.sort();
        columns.push(HtsPileupColumn { pos, depth: p.depth(), alignments });
    }
    columns
}

// r[verify pileup.htslib_compat]
// r[verify pileup.position_iteration]
#[test]
fn pileup_positions_match() {
    let hts_columns = fetch_htslib_pileup();
    assert!(!hts_columns.is_empty(), "htslib produced no pileup columns — check test data");

    let bam_path = test_bam_path();
    let mut reader = seqair::bam::IndexedBamReader::open(bam_path).expect("seqair open");
    let mut store = seqair::bam::RecordStore::new();
    let tid = reader.header().tid(TEST_REGION).expect("tid");
    reader
        .fetch_into(
            tid,
            Pos::<Zero>::new(TEST_START as u32),
            Pos::<Zero>::new(TEST_END as u32),
            &mut store,
        )
        .expect("seqair fetch");

    let engine = seqair::bam::PileupEngine::new(
        store,
        Pos::<Zero>::new(TEST_START as u32),
        Pos::<Zero>::new(TEST_END as u32),
    );
    let columns: Vec<_> = engine.collect();

    let hts_positions: Vec<u32> = hts_columns.iter().map(|c| c.pos).collect();
    let positions: Vec<u32> = columns.iter().map(|c| c.pos().get()).collect();

    assert_eq!(
        positions.len(),
        hts_positions.len(),
        "column count mismatch: seqair={} htslib={}",
        positions.len(),
        hts_positions.len()
    );

    for (i, (pos, hts_pos)) in positions.iter().zip(hts_positions.iter()).enumerate() {
        assert_eq!(*pos, *hts_pos, "position mismatch at column {i}");
    }
}

// r[verify pileup.column_contents]
// r[verify pileup.active_set]
#[test]
fn pileup_depth_matches() {
    let hts_columns = fetch_htslib_pileup();

    let bam_path = test_bam_path();
    let mut reader = seqair::bam::IndexedBamReader::open(bam_path).expect("seqair open");
    let mut store = seqair::bam::RecordStore::new();
    let tid = reader.header().tid(TEST_REGION).expect("tid");
    reader
        .fetch_into(
            tid,
            Pos::<Zero>::new(TEST_START as u32),
            Pos::<Zero>::new(TEST_END as u32),
            &mut store,
        )
        .expect("seqair fetch");

    let engine = seqair::bam::PileupEngine::new(
        store,
        Pos::<Zero>::new(TEST_START as u32),
        Pos::<Zero>::new(TEST_END as u32),
    );
    let columns: Vec<_> = engine.collect();

    for (i, (rio, hts)) in columns.iter().zip(hts_columns.iter()).enumerate() {
        // htslib depth includes deletions; alignments.len() only counts bases with qpos.
        // seqair depth should match the non-deletion count.
        assert!(
            hts.depth as usize >= hts.alignments.len(),
            "htslib depth < non-del count at position {} (column {i}): depth={} non_del={}",
            hts.pos,
            hts.depth,
            hts.alignments.len()
        );
        assert_eq!(
            rio.depth(),
            hts.alignments.len(),
            "depth mismatch at position {} (column {i}): seqair={} htslib={}",
            hts.pos,
            rio.depth(),
            hts.alignments.len()
        );
    }
}

// r[verify pileup.qpos]
// r[verify pileup.qpos_none]
// r[verify cigar.qpos_at]
// r[verify cigar.qpos_accuracy]
#[test]
fn pileup_qpos_matches() {
    let hts_columns = fetch_htslib_pileup();

    let bam_path = test_bam_path();
    let mut reader = seqair::bam::IndexedBamReader::open(bam_path).expect("seqair open");
    let mut store = seqair::bam::RecordStore::new();
    let tid = reader.header().tid(TEST_REGION).expect("tid");
    reader
        .fetch_into(
            tid,
            Pos::<Zero>::new(TEST_START as u32),
            Pos::<Zero>::new(TEST_END as u32),
            &mut store,
        )
        .expect("seqair fetch");

    let engine = seqair::bam::PileupEngine::new(
        store,
        Pos::<Zero>::new(TEST_START as u32),
        Pos::<Zero>::new(TEST_END as u32),
    );
    let columns: Vec<_> = engine.collect();

    for (col_idx, (rio, hts)) in columns.iter().zip(hts_columns.iter()).enumerate() {
        let mut alns: Vec<(usize, u16)> =
            rio.alignments().filter_map(|a| a.qpos().map(|q| (q, a.flags))).collect();
        alns.sort();

        assert_eq!(
            alns.len(),
            hts.alignments.len(),
            "alignment count mismatch at position {} (column {col_idx})",
            hts.pos,
        );

        for (aln_idx, (aln, hts_aln)) in alns.iter().zip(hts.alignments.iter()).enumerate() {
            assert_eq!(
                aln.0, hts_aln.0,
                "qpos mismatch at position {} alignment {aln_idx}: rio={} hts={}",
                hts.pos, aln.0, hts_aln.0
            );
            assert_eq!(
                aln.1, hts_aln.1,
                "flags mismatch at position {} alignment {aln_idx}",
                hts.pos,
            );
        }
    }
}
