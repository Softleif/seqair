//! Comparison tests: Pileup engine.
//! Builds pileups over the same region using both htslib and seqair,
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    dead_code
)]

use rust_htslib::bam::pileup::Indel;
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

struct HtsFullAlignment {
    qpos: Option<usize>,
    flags: u16,
    is_del: bool,
    is_refskip: bool,
    indel: Indel,
}

struct HtsFullColumn {
    pos: u32,
    depth: u32,
    alignments: Vec<HtsFullAlignment>,
}

fn fetch_htslib_pileup_region(region: &str, start: u64, end: u64) -> Vec<HtsFullColumn> {
    let bam_path = test_bam_path();
    let mut reader = bam::IndexedReader::from_path(bam_path).expect("htslib open");
    let tid = reader.header().tid(region.as_bytes()).expect("tid");
    reader
        .fetch(FetchDefinition::Region(tid as i32, start as i64, end as i64))
        .expect("htslib fetch");
    let mut columns = Vec::new();
    let pileup = reader.pileup();
    for p in pileup {
        let p = p.expect("htslib pileup");
        let pos = p.pos();
        if (pos as u64) < start || (pos as u64) > end {
            continue;
        }
        let alignments = p
            .alignments()
            .map(|a| HtsFullAlignment {
                qpos: a.qpos(),
                flags: a.record().flags(),
                is_del: a.is_del(),
                is_refskip: a.is_refskip(),
                indel: a.indel(),
            })
            .collect();
        columns.push(HtsFullColumn { pos, depth: p.depth(), alignments });
    }
    columns
}

fn fetch_seqair_pileup_region(
    region: &str,
    start: u64,
    end: u64,
) -> Vec<seqair::bam::PileupColumn> {
    let bam_path = test_bam_path();
    let mut reader = seqair::bam::IndexedBamReader::open(bam_path).expect("seqair open");
    let mut store = seqair::bam::RecordStore::new();
    let tid = reader.header().tid(region).expect("tid");
    reader
        .fetch_into(
            tid,
            Pos::<Zero>::new(start as u32).unwrap(),
            Pos::<Zero>::new(end as u32).unwrap(),
            &mut store,
        )
        .expect("seqair fetch");
    seqair::bam::PileupEngine::new(
        store,
        Pos::<Zero>::new(start as u32).unwrap(),
        Pos::<Zero>::new(end as u32).unwrap(),
    )
    .collect()
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
            Pos::<Zero>::new(TEST_START as u32).unwrap(),
            Pos::<Zero>::new(TEST_END as u32).unwrap(),
            &mut store,
        )
        .expect("seqair fetch");

    let engine = seqair::bam::PileupEngine::new(
        store,
        Pos::<Zero>::new(TEST_START as u32).unwrap(),
        Pos::<Zero>::new(TEST_END as u32).unwrap(),
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
// Note: this test validates match_depth() specifically (qpos-bearing alignments only).
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
            Pos::<Zero>::new(TEST_START as u32).unwrap(),
            Pos::<Zero>::new(TEST_END as u32).unwrap(),
            &mut store,
        )
        .expect("seqair fetch");

    let engine = seqair::bam::PileupEngine::new(
        store,
        Pos::<Zero>::new(TEST_START as u32).unwrap(),
        Pos::<Zero>::new(TEST_END as u32).unwrap(),
    );
    let columns: Vec<_> = engine.collect();

    for (i, (rio, hts)) in columns.iter().zip(hts_columns.iter()).enumerate() {
        // htslib depth includes deletions; alignments.len() only counts bases with qpos.
        // The non-deletion alignment count (match_depth) must match htslib's qpos-bearing count.
        assert!(
            hts.depth as usize >= hts.alignments.len(),
            "htslib depth < non-del count at position {} (column {i}): depth={} non_del={}",
            hts.pos,
            hts.depth,
            hts.alignments.len()
        );
        assert_eq!(
            rio.match_depth(),
            hts.alignments.len(),
            "match_depth mismatch at position {} (column {i}): seqair={} htslib={}",
            hts.pos,
            rio.match_depth(),
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
            Pos::<Zero>::new(TEST_START as u32).unwrap(),
            Pos::<Zero>::new(TEST_END as u32).unwrap(),
            &mut store,
        )
        .expect("seqair fetch");

    let engine = seqair::bam::PileupEngine::new(
        store,
        Pos::<Zero>::new(TEST_START as u32).unwrap(),
        Pos::<Zero>::new(TEST_END as u32).unwrap(),
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

// r[verify pileup_indel.depth_includes_all]
// r[verify pileup_indel.htslib_compat_update]
#[test]
fn total_depth_with_deletions_matches_htslib() {
    let hts_cols = fetch_htslib_pileup_region("bacteriophage_lambda_CpG", 73, 200);
    let seq_cols = fetch_seqair_pileup_region("bacteriophage_lambda_CpG", 73, 200);

    assert!(!hts_cols.is_empty(), "htslib produced no columns — check test data");
    let total_del_alns: usize =
        hts_cols.iter().flat_map(|c| c.alignments.iter()).filter(|a| a.is_del).count();
    assert!(total_del_alns > 0, "test region has no deletion alignments — test data problem");

    assert_eq!(hts_cols.len(), seq_cols.len(), "column count mismatch");
    for (i, (hts, seq)) in hts_cols.iter().zip(seq_cols.iter()).enumerate() {
        assert_eq!(hts.pos, seq.pos().get(), "pos mismatch at column {i}");
        assert_eq!(
            hts.depth as usize,
            seq.depth(),
            "total depth (including deletions) mismatch at pos {} (column {i}): htslib={} seqair={}",
            hts.pos,
            hts.depth,
            seq.depth()
        );
    }
}

// r[verify pileup_indel.deletions_included]
// r[verify pileup_indel.refskips_included]
// r[verify pileup_indel.accessors]
#[test]
fn deletion_ops_match_htslib() {
    let hts_cols = fetch_htslib_pileup_region("bacteriophage_lambda_CpG", 73, 200);
    let seq_cols = fetch_seqair_pileup_region("bacteriophage_lambda_CpG", 73, 200);

    assert_eq!(hts_cols.len(), seq_cols.len(), "column count mismatch");

    // Build position lookup for seqair columns to enable anchor→deletion cross-validation.
    let seq_by_pos: std::collections::HashMap<u32, &seqair::bam::PileupColumn> =
        seq_cols.iter().map(|c| (c.pos().get(), c)).collect();

    for (i, (hts, seq)) in hts_cols.iter().zip(seq_cols.iter()).enumerate() {
        let hts_del = hts.alignments.iter().filter(|a| a.is_del).count();
        let seq_del = seq.alignments().filter(|a| a.is_del()).count();
        assert_eq!(
            hts_del, seq_del,
            "deletion alignment count mismatch at pos {} (column {i}): htslib={hts_del} seqair={seq_del}",
            hts.pos,
        );

        // Verify seqair del_len > 0 for all deletion alignments.
        if seq_del > 0 {
            let seq_del_lens: Vec<u32> =
                seq.alignments().filter(|a| a.is_del()).map(|a| a.del_len()).collect();
            assert!(
                seq_del_lens.iter().all(|&len| len > 0),
                "deletion alignments must have nonzero del_len at pos {} (column {i})",
                hts.pos,
            );
        }

        // Cross-validate del_len against htslib anchor positions:
        // htslib reports Indel::Del(len) at the last M position before the deletion starts.
        // At position P+1, seqair must have at least as many Deletion{del_len: len} as htslib
        // anchor reads reporting Del(len) at pos P.
        let next_pos = hts.pos + 1;
        if let Some(next_col) = seq_by_pos.get(&next_pos) {
            // Group htslib Del(len) anchors by length.
            let mut hts_del_counts: std::collections::HashMap<u32, usize> =
                std::collections::HashMap::new();
            for hts_aln in &hts.alignments {
                if let Indel::Del(hts_del_len) = hts_aln.indel {
                    *hts_del_counts.entry(hts_del_len).or_insert(0) += 1;
                }
            }
            for (hts_del_len, hts_count) in &hts_del_counts {
                let seq_count = next_col
                    .alignments()
                    .filter(|a| a.is_del() && a.del_len() == *hts_del_len)
                    .count();
                assert!(
                    seq_count >= *hts_count,
                    "htslib reports {hts_count} Del({hts_del_len}) anchor(s) at pos {}, but seqair \
                     has only {seq_count} matching del_len at pos {next_pos} (column {i})",
                    hts.pos,
                );
            }
        }

        let hts_refskip = hts.alignments.iter().filter(|a| a.is_refskip).count();
        let seq_refskip = seq.alignments().filter(|a| a.is_refskip()).count();
        assert_eq!(
            hts_refskip, seq_refskip,
            "refskip alignment count mismatch at pos {} (column {i}): htslib={hts_refskip} seqair={seq_refskip}",
            hts.pos,
        );
    }
}

// r[verify pileup_indel.insertion_at_last_match]
// r[verify pileup_indel.insertion_len]
#[test]
fn insertion_ops_match_htslib() {
    let hts_cols = fetch_htslib_pileup_region("chr19", 6_110_698, 6_111_300);
    let seq_cols = fetch_seqair_pileup_region("chr19", 6_110_698, 6_111_300);

    let total_ins_alns: usize = hts_cols
        .iter()
        .flat_map(|c| c.alignments.iter())
        .filter(|a| matches!(a.indel, Indel::Ins(_)))
        .count();
    assert!(total_ins_alns > 0, "test region has no insertion alignments — test data problem");

    assert_eq!(hts_cols.len(), seq_cols.len(), "column count mismatch");
    for (i, (hts, seq)) in hts_cols.iter().zip(seq_cols.iter()).enumerate() {
        // Collect (qpos, insert_len) for insertion alignments from both engines.
        let mut hts_ins: Vec<(usize, u32)> =
            hts.alignments
                .iter()
                .filter_map(|a| {
                    if let Indel::Ins(len) = a.indel { a.qpos.map(|q| (q, len)) } else { None }
                })
                .collect();
        hts_ins.sort_unstable();

        let mut seq_ins: Vec<(usize, u32)> = seq
            .alignments()
            .filter(|a| a.insert_len() > 0)
            .filter_map(|a| a.qpos().map(|q| (q, a.insert_len())))
            .collect();
        seq_ins.sort_unstable();

        assert_eq!(
            hts_ins.len(),
            seq_ins.len(),
            "insertion count mismatch at pos {} (column {i}): htslib={} seqair={}",
            hts.pos,
            hts_ins.len(),
            seq_ins.len()
        );
        for (ins_idx, (hts_ins_aln, seq_ins_aln)) in hts_ins.iter().zip(seq_ins.iter()).enumerate()
        {
            assert_eq!(
                hts_ins_aln.0, seq_ins_aln.0,
                "insertion qpos mismatch at pos {} (column {i}, ins {ins_idx}): htslib={} seqair={}",
                hts.pos, hts_ins_aln.0, seq_ins_aln.0
            );
            assert_eq!(
                hts_ins_aln.1, seq_ins_aln.1,
                "insertion length mismatch at pos {} (column {i}, ins {ins_idx}): htslib={} seqair={}",
                hts.pos, hts_ins_aln.1, seq_ins_aln.1
            );
        }
    }
}
