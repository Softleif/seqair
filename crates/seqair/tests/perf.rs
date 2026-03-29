//! Tests for performance-related spec rules.
//! after refactoring) rather than measuring wall-clock time.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::indexing_slicing)]
mod helpers;

use helpers::{cigar_bytes, cigar_op, make_record, make_record_with_cigar};
use proptest::prelude::*;
use seqair::bam::cigar::{CigarMapping, CigarPosInfo};
use seqair::bam::{RecordStore, pileup::PileupEngine};

// ---- perf.reuse_alignment_vec ----
// Verified indirectly: if the engine reuses the vec internally, output must
// still be correct. This test ensures correctness is maintained after the
// optimization.

// r[verify perf.reuse_alignment_vec+2]
proptest! {
    #[test]
    fn reused_vec_depth_matches_sweep_line(
        records in prop::collection::vec((0i32..200, 20u32..=80), 2..=15),
    ) {
        // Sort by position to match fetch_into's coordinate-sorted guarantee.
        let mut sorted_records = records.clone();
        sorted_records.sort_by_key(|&(offset, _)| offset);

        let mut arena = RecordStore::new();
        // Sweep-line: +1 at read start, -1 at read end+1. Prefix sum = depth.
        const MAX_POS: usize = 301;
        let mut delta = vec![0i32; MAX_POS + 1];
        for &(offset, len) in &sorted_records {
            let start = offset as usize;
            let end_excl = (offset as usize + len as usize).min(MAX_POS);
            if start < MAX_POS {
                delta[start] += 1;
            }
            if end_excl <= MAX_POS {
                delta[end_excl] -= 1;
            }
            arena.push_raw(&make_record(0, offset, 99, 60, len)).unwrap();
        }

        let mut expected_depth = vec![0usize; MAX_POS];
        let mut running = 0i32;
        for (i, &d) in delta.iter().enumerate().take(MAX_POS) {
            running += d;
            expected_depth[i] = running.max(0) as usize;
        }

        let engine = PileupEngine::new(arena, 0, 300);
        for col in engine {
            let pos = col.pos() as usize;
            let exp = expected_depth.get(pos).copied().unwrap_or(0);
            prop_assert_eq!(col.depth(), exp, "depth wrong at pos {}", col.pos());
        }
    }
}

// ---- perf.no_sorted_indices ----
// The engine must produce correct output when records are already sorted
// (which they always are from fetch_into).

// r[verify perf.no_sorted_indices]
#[test]
fn already_sorted_records_produce_correct_pileup() {
    let mut arena = RecordStore::new();
    // Push in sorted order (as fetch_into guarantees)
    arena.push_raw(&make_record(0, 10, 99, 60, 30)).unwrap();
    arena.push_raw(&make_record(0, 20, 99, 60, 30)).unwrap();
    arena.push_raw(&make_record(0, 30, 99, 60, 30)).unwrap();

    let engine = PileupEngine::new(arena, 0, 70);
    let columns: Vec<_> = engine.collect();

    // pos 10-19: 1 read
    assert_eq!(columns.iter().find(|c| c.pos() == 15).unwrap().depth(), 1);
    // pos 20-29: 2 reads
    assert_eq!(columns.iter().find(|c| c.pos() == 25).unwrap().depth(), 2);
    // pos 30-39: 3 reads
    assert_eq!(columns.iter().find(|c| c.pos() == 35).unwrap().depth(), 3);
    // pos 40-49: 2 reads (first ended)
    assert_eq!(columns.iter().find(|c| c.pos() == 45).unwrap().depth(), 2);
}

// ---- perf.avoid_redundant_arena_get ----
// Correctness preserved when arena.get() is called fewer times.

// r[verify perf.avoid_redundant_arena_get+2]
#[test]
fn single_arena_get_per_record_entry_still_correct() {
    let mut arena = RecordStore::new();
    arena.push_raw(&make_record(0, 0, 99, 60, 50)).unwrap();
    arena.push_raw(&make_record(0, 10, 99 | 0x100, 40, 50)).unwrap(); // secondary flag

    let mut engine = PileupEngine::new(arena, 0, 59);
    engine.set_filter(|flags, _aux| flags & 0x100 == 0);
    let columns: Vec<_> = engine.collect();

    // Only the mapq=60 read should pass the filter
    for col in &columns {
        if col.pos() < 10 {
            assert_eq!(col.depth(), 1, "only high-mapq read at pos {}", col.pos());
        }
    }
    // At pos 10+, still only 1 because second read is filtered
    let col = columns.iter().find(|c| c.pos() == 25).unwrap();
    assert_eq!(col.depth(), 1);
}

// ---- perf.cigar_no_to_vec ----
// After removing .to_vec(), CigarMapping must still produce correct qpos.

// r[verify perf.cigar_no_to_vec]
#[test]
fn cigar_index_from_arena_slab_correct() {
    let mut arena = RecordStore::new();
    // 30M 5D 20M — tests that cigar bytes from arena slab work without copying
    let raw = make_record_with_cigar(
        0,
        100,
        99,
        60,
        &[cigar_op(30, 0), cigar_op(5, 2), cigar_op(20, 0)],
        50,
    );
    arena.push_raw(&raw).unwrap();

    let engine = PileupEngine::new(arena, 100, 154);
    let columns: Vec<_> = engine.collect();

    // Before deletion: qpos = pos - 100
    let col = columns.iter().find(|c| c.pos() == 110).unwrap();
    assert_eq!(col.alignments().next().unwrap().qpos(), Some(10));

    // Inside deletion: alignment present with Deletion op (no qpos)
    let del_col = columns.iter().find(|c| c.pos() == 132).unwrap();
    assert_eq!(del_col.depth(), 1);
    assert!(del_col.alignments().next().unwrap().is_del());

    // After deletion: qpos = pos - 100 - 5 (deletion consumes 5 ref but 0 query)
    let col = columns.iter().find(|c| c.pos() == 140).unwrap();
    assert_eq!(col.alignments().next().unwrap().qpos(), Some(35));
}

// ---- perf.precompute_matches_indels ----

// r[verify perf.precompute_matches_indels]
#[test]
fn precomputed_matches_indels_accessible_from_record() {
    let mut arena = RecordStore::new();
    // 30M 5I 15M = 45 matches, 5 indels
    let raw = make_record_with_cigar(
        0,
        0,
        99,
        60,
        &[cigar_op(30, 0), cigar_op(5, 1), cigar_op(15, 0)],
        50,
    );
    arena.push_raw(&raw).unwrap();

    let r = arena.record(0);
    assert_eq!(r.matching_bases, 45);
    assert_eq!(r.indel_bases, 5);
}

// r[verify perf.precompute_matches_indels]
proptest! {
    #[test]
    fn matches_indels_from_varied_cigar_ops(
        // Generate a sequence of (op_type, length) pairs from {M=0, I=1, D=2, S=4}.
        // S is only valid at the ends, so we place it as leading/trailing clips.
        leading_clip in 0u32..=10,
        trailing_clip in 0u32..=10,
        inner_ops in prop::collection::vec(
            (prop::sample::select(vec![0u8, 1u8, 2u8]), 1u32..=50u32),
            1..=8,
        ),
    ) {
        // Build CIGAR: optional leading S, inner ops {M,I,D}, optional trailing S.
        // Per r[cigar.matches_indels]: matching_bases = sum of M(0) lengths only,
        // indel_bases = sum of I(1) + D(2) lengths.
        let mut ops: Vec<u32> = Vec::new();
        let mut seq_len = 0u32;

        if leading_clip > 0 {
            ops.push(cigar_op(leading_clip, 4)); // S
            seq_len += leading_clip;
        }

        // Ensure at least one M op so the record has a ref-consuming op.
        let mut has_m = false;
        let mut exp_matches = 0u32;
        let mut exp_indels = 0u32;

        for &(op_type, len) in &inner_ops {
            ops.push(cigar_op(len, op_type));
            match op_type {
                0 => { exp_matches += len; seq_len += len; has_m = true; }
                1 => { exp_indels += len; seq_len += len; }
                2 => { exp_indels += len; } // D consumes ref, not query
                _ => {}
            }
        }

        // Guarantee at least one M op so the record is usable.
        if !has_m {
            ops.push(cigar_op(1, 0)); // M
            seq_len += 1;
            exp_matches += 1;
        }

        if trailing_clip > 0 {
            ops.push(cigar_op(trailing_clip, 4)); // S
            seq_len += trailing_clip;
        }

        prop_assume!(seq_len > 0);

        let raw = make_record_with_cigar(0, 0, 99, 60, &ops, seq_len);
        let mut arena = RecordStore::new();
        arena.push_raw(&raw).unwrap();
        let r = arena.record(0);

        prop_assert_eq!(r.matching_bases, exp_matches,
            "matching_bases mismatch for ops={:?}", inner_ops);
        prop_assert_eq!(r.indel_bases, exp_indels,
            "indel_bases mismatch for ops={:?}", inner_ops);
    }
}

// ---- perf.cigar_binary_search ----

// r[verify perf.cigar_binary_search]
#[test]
fn binary_search_correct_for_many_ops() {
    // Simulate RNA-seq: 30M 5000N 30M 3000N 40M (5 ops, triggers binary search)
    let ops = cigar_bytes(&[
        cigar_op(30, 0),
        cigar_op(5000, 3),
        cigar_op(30, 0),
        cigar_op(3000, 3),
        cigar_op(40, 0),
    ]);
    let mapping = CigarMapping::new(1000, &ops);

    // First M block: 1000-1029
    assert_eq!(mapping.pos_info_at(1000), Some(CigarPosInfo::Match { qpos: 0 }));
    assert_eq!(mapping.pos_info_at(1029), Some(CigarPosInfo::Match { qpos: 29 }));

    // N skip: 1030-6029 → RefSkip
    assert_eq!(mapping.pos_info_at(1030), Some(CigarPosInfo::RefSkip));
    assert_eq!(mapping.pos_info_at(5000), Some(CigarPosInfo::RefSkip));

    // Second M block: 6030-6059
    assert_eq!(mapping.pos_info_at(6030), Some(CigarPosInfo::Match { qpos: 30 }));
    assert_eq!(mapping.pos_info_at(6059), Some(CigarPosInfo::Match { qpos: 59 }));

    // Second N skip: 6060-9059 → RefSkip
    assert_eq!(mapping.pos_info_at(6060), Some(CigarPosInfo::RefSkip));

    // Third M block: 9060-9099
    assert_eq!(mapping.pos_info_at(9060), Some(CigarPosInfo::Match { qpos: 60 }));
    assert_eq!(mapping.pos_info_at(9099), Some(CigarPosInfo::Match { qpos: 99 }));
    assert_eq!(mapping.pos_info_at(9100), None);
}

// r[verify perf.cigar_binary_search]
proptest! {
    #[test]
    fn binary_search_correct_for_rna_seq_pattern(
        n_exons in 3usize..=8,
        exon_len in 20u32..=100,
        intron_len in 100u32..=5000,
    ) {
        // Build alternating M/N CIGAR (RNA-seq pattern, triggers bsearch path for ≥5 ops)
        let mut ops = Vec::new();
        let mut total_query = 0u32;
        for i in 0..n_exons {
            ops.push(cigar_op(exon_len, 0)); // M
            total_query += exon_len;
            if i < n_exons - 1 {
                ops.push(cigar_op(intron_len, 3)); // N
            }
        }
        let bytes = cigar_bytes(&ops);
        let mapping = CigarMapping::new(0, &bytes);

        // Compute total ref span
        let total_ref: i64 = ops.iter().map(|&op| {
            let len = i64::from(op >> 4);
            let op_type = (op & 0xF) as u8;
            if matches!(op_type, 0 | 2 | 3 | 7 | 8) { len } else { 0 }
        }).sum();

        // Check every position — qpos must be in range and monotonically increasing
        let mut last_qpos: Option<u32> = None;
        for pos in 0..total_ref {
            match mapping.pos_info_at(pos) {
                Some(CigarPosInfo::Match { qpos }) => {
                    prop_assert!((qpos as usize) < total_query as usize,
                        "qpos {qpos} out of range at ref {pos}");
                    if let Some(prev) = last_qpos {
                        prop_assert!(qpos > prev, "qpos not monotonic at ref {pos}");
                    }
                    last_qpos = Some(qpos);
                }
                Some(CigarPosInfo::RefSkip) => {}
                other => prop_assert!(false, "unexpected result at ref {pos}: {other:?}"),
            }
        }
        // Out of range
        prop_assert_eq!(mapping.pos_info_at(total_ref), None);
    }
}

// ---- perf.arena_capacity_hint ----

// r[verify perf.arena_capacity_hint+2]
#[test]
fn arena_with_capacity_avoids_realloc() {
    // Pre-size based on estimated compressed bytes
    let mut arena = RecordStore::with_byte_hint(50_000);

    for i in 0..100 {
        arena.push_raw(&make_record(0, i * 10, 99, 60, 50)).unwrap();
    }
    assert_eq!(arena.len(), 100);

    // All records accessible
    for i in 0..100u32 {
        let r = arena.record(i);
        assert_eq!(r.pos, i64::from(i) * 10);
    }
}
