//! Tests for PileupEngine: column iteration, filtering, max depth, edge cases.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::indexing_slicing)]
mod helpers;

use helpers::{cigar_op, make_record, make_record_with_cigar};
use proptest::prelude::*;
use seqair::bam::pileup::RefSeq;
use seqair::bam::{RecordStore, pileup::PileupEngine};
use seqair_types::Base;
use std::{cell::Cell, rc::Rc};

// ---- pileup.active_set + pileup.column_contents ----

// r[verify pileup.active_set]
// r[verify pileup.column_contents]
proptest! {
    #[test]
    fn depth_matches_sweep_line_count(
        records in prop::collection::vec((0i32..500, 10u32..=100), 1..=20),
    ) {
        let mut sorted = records.clone();
        sorted.sort_by_key(|&(offset, _)| offset);

        let mut arena = RecordStore::new();
        // Build a sweep-line depth array: for each read, +1 at start, -1 at end+1.
        // Prefix-sum gives expected depth at each position.
        const MAX_POS: usize = 601;
        let mut delta = vec![0i32; MAX_POS + 1];
        for &(offset, read_len) in &sorted {
            let start = offset as usize;
            let end_excl = (offset as usize + read_len as usize).min(MAX_POS);
            if start < MAX_POS {
                delta[start] += 1;
            }
            if end_excl <= MAX_POS {
                delta[end_excl] -= 1;
            }
            arena.push_raw(&make_record(0, offset, 99, 60, read_len)).unwrap();
        }

        // Compute prefix sum to get expected depth at each position.
        let mut expected_depth = vec![0usize; MAX_POS];
        let mut running = 0i32;
        for (i, &d) in delta.iter().enumerate().take(MAX_POS) {
            running += d;
            expected_depth[i] = running.max(0) as usize;
        }

        let engine = PileupEngine::new(arena, 0, 600);
        for col in engine {
            let pos = col.pos() as usize;
            let exp = expected_depth.get(pos).copied().unwrap_or(0);
            prop_assert_eq!(col.depth(), exp,
                "depth mismatch at pos {}", col.pos());
        }
    }
}

// ---- pileup.position_iteration ----

// r[verify pileup.position_iteration]
proptest! {
    #[test]
    fn columns_match_sweep_line_coverage(
        records in prop::collection::vec((0i32..500, 10u32..=100), 1..=20),
    ) {
        let mut sorted = records.clone();
        sorted.sort_by_key(|&(offset, _)| offset);

        let mut arena = RecordStore::new();
        // Build sweep-line depth array using event-based approach.
        const MAX_POS: usize = 601;
        let mut delta = vec![0i32; MAX_POS + 1];
        for &(offset, read_len) in &sorted {
            let start = offset as usize;
            let end_excl = (offset as usize + read_len as usize).min(MAX_POS);
            if start < MAX_POS {
                delta[start] += 1;
            }
            if end_excl <= MAX_POS {
                delta[end_excl] -= 1;
            }
            arena.push_raw(&make_record(0, offset, 99, 60, read_len)).unwrap();
        }

        // Prefix sum → expected depth per position.
        let mut expected_depth = vec![0i32; MAX_POS];
        let mut running = 0i32;
        for (i, &d) in delta.iter().enumerate().take(MAX_POS) {
            running += d;
            expected_depth[i] = running.max(0);
        }

        let engine = PileupEngine::new(arena, 0, 600);
        let columns: Vec<_> = engine.collect();
        let col_positions: std::collections::HashSet<i64> =
            columns.iter().map(|c| c.pos()).collect();

        for pos in 0..MAX_POS as i64 {
            let exp_covered = expected_depth.get(pos as usize).copied().unwrap_or(0) > 0;
            let actual_covered = col_positions.contains(&pos);
            prop_assert_eq!(actual_covered, exp_covered,
                "position {} coverage mismatch (sweep says covered={})", pos, exp_covered);
        }
    }
}

// ---- pileup.read_filter ----

// r[verify pileup.read_filter]
#[test]
fn filter_evaluated_once_per_record() {
    let mut arena = RecordStore::new();
    arena.push_raw(&make_record(0, 0, 99, 60, 20)).unwrap();

    let count = Rc::new(Cell::new(0usize));
    let count_clone = count.clone();
    let mut engine = PileupEngine::new(arena, 0, 19);
    engine.set_filter(move |_flags, _aux| {
        count_clone.set(count_clone.get() + 1);
        true
    });

    let columns: Vec<_> = engine.collect();
    assert_eq!(columns.len(), 20);
    assert_eq!(count.get(), 1, "filter called once per record, not per column");
}

// r[verify pileup.read_filter]
proptest! {
    #[test]
    fn filter_by_flags_excludes_correct_reads(
        pass_flags in prop::collection::vec(prop::bool::ANY, 1..=20),
    ) {
        let mut arena = RecordStore::new();
        for &pass in &pass_flags {
            let flags = if pass { 99 } else { 99 | 0x100 }; // 0x100 = secondary
            arena.push_raw(&make_record(0, 0, flags, 60, 10)).unwrap();
        }

        let expected = pass_flags.iter().filter(|&&p| p).count();
        let mut engine = PileupEngine::new(arena, 0, 0);
        engine.set_filter(move |flags, _aux| flags & 0x100 == 0);
        let columns: Vec<_> = engine.collect();

        if expected > 0 {
            prop_assert_eq!(columns[0].depth(), expected);
        } else {
            prop_assert!(columns.is_empty());
        }
    }
}

// ---- pileup.max_depth + pileup.max_depth_per_position ----

// r[verify pileup.max_depth]
proptest! {
    #[test]
    fn max_depth_never_exceeded(max in 1u32..=10, n in 1usize..=30, len in 10u32..=50) {
        let mut arena = RecordStore::new();
        for _ in 0..n {
            arena.push_raw(&make_record(0, 0, 99, 60, len)).unwrap();
        }

        let mut engine = PileupEngine::new(arena, 0, i64::from(len) - 1);
        engine.set_max_depth(max);
        for col in engine {
            prop_assert!(col.depth() <= max as usize);
        }
    }
}

// r[verify pileup.max_depth_per_position]
proptest! {
    #[test]
    fn reads_available_at_lower_coverage_despite_cap(extra in 3usize..=10) {
        let mut arena = RecordStore::new();
        for _ in 0..2 {
            arena.push_raw(&make_record(0, 0, 99, 60, 100)).unwrap();
        }
        for _ in 0..extra {
            arena.push_raw(&make_record(0, 0, 99, 60, 20)).unwrap();
        }

        let mut engine = PileupEngine::new(arena, 0, 99);
        engine.set_max_depth(3);
        let columns: Vec<_> = engine.collect();

        let col0 = columns.iter().find(|c| c.pos() == 0).unwrap();
        prop_assert_eq!(col0.depth(), 3); // capped

        let col50 = columns.iter().find(|c| c.pos() == 50).unwrap();
        prop_assert_eq!(col50.depth(), 2); // below cap
    }
}

// ---- pileup.qpos ----

// r[verify pileup.qpos]
proptest! {
    #[test]
    fn qpos_for_simple_cigar_is_offset_from_start(start in 0i32..100, len in 10u32..=50) {
        let mut arena = RecordStore::new();
        arena.push_raw(&make_record(0, start, 99, 60, len)).unwrap();

        let start_i64 = i64::from(start);
        let end_i64 = start_i64 + i64::from(len) - 1;
        let engine = PileupEngine::new(arena, start_i64, end_i64);

        for col in engine {
            for aln in col.alignments() {
                prop_assert_eq!(aln.qpos().unwrap(), (col.pos() - start_i64) as usize);
            }
        }
    }
}

// ---- pileup.empty_positions_skipped ----

// r[verify pileup.empty_positions_skipped]
proptest! {
    #[test]
    fn no_empty_columns(pos1 in 0i32..100, len1 in 5u32..=20, gap in 50i32..=500, len2 in 5u32..=20) {
        let pos2 = pos1 + i32::from(len1 as u16) + gap;
        let mut arena = RecordStore::new();
        arena.push_raw(&make_record(0, pos1, 99, 60, len1)).unwrap();
        arena.push_raw(&make_record(0, pos2, 99, 60, len2)).unwrap();

        let end = i64::from(pos2) + i64::from(len2) + 10;
        let engine = PileupEngine::new(arena, i64::from(pos1), end);
        for col in engine {
            prop_assert!(col.depth() > 0, "empty column at {}", col.pos());
        }
    }
}

// ---- pileup.unmapped_excluded ----

// r[verify pileup.unmapped_excluded]
#[test]
fn unmapped_reads_excluded_from_pileup() {
    let mut arena = RecordStore::new();
    arena.push_raw(&make_record(0, 100, 99, 60, 50)).unwrap();
    arena.push_raw(&make_record(0, 100, 0x4, 60, 50)).unwrap(); // unmapped
    arena.push_raw(&make_record(0, 100, 163, 60, 50)).unwrap();

    let engine = PileupEngine::new(arena, 100, 149);
    let columns: Vec<_> = engine.collect();
    assert_eq!(columns[0].depth(), 2, "unmapped should be excluded");
}

// ---- pileup.zero_refspan_reads ----

// r[verify pileup.zero_refspan_reads]
#[test]
fn zero_refspan_read_handled_gracefully() {
    let raw = make_record_with_cigar(0, 50, 99, 60, &[cigar_op(10, 4)], 10); // 10S
    let mut arena = RecordStore::new();
    arena.push_raw(&raw).unwrap();
    assert_eq!(arena.record(0).end_pos, 50);

    let engine = PileupEngine::new(arena, 48, 52);
    let columns: Vec<_> = engine.collect();
    // Pure soft-clip has no ref-consuming ops → no qpos → no column
    assert!(columns.is_empty() || columns.iter().all(|c| c.pos() != 50 || c.depth() == 0));
}

// ---- pileup.soft_clip_at_position ----

// r[verify pileup.soft_clip_at_position]
#[test]
fn leading_softclip_does_not_extend_ref_range() {
    // 5S 20M at pos 100 → ref 100-119, query 5-24
    let raw = make_record_with_cigar(0, 100, 99, 60, &[cigar_op(5, 4), cigar_op(20, 0)], 25);
    let mut arena = RecordStore::new();
    arena.push_raw(&raw).unwrap();

    let engine = PileupEngine::new(arena, 95, 125);
    let columns: Vec<_> = engine.collect();
    let positions: Vec<i64> = columns.iter().map(|c| c.pos()).collect();

    assert_eq!(positions.first(), Some(&100), "should start at first ref-consuming pos");
    assert_eq!(positions.last(), Some(&119));
    assert_eq!(columns[0].alignments().next().unwrap().qpos().unwrap(), 5);
}

// ---- bam.reader edge cases ----

// r[verify bam.reader.secondary_supplementary_included+2]
#[test]
fn secondary_and_supplementary_reads_in_arena() {
    let mut arena = RecordStore::new();
    arena.push_raw(&make_record(0, 100, 99, 60, 50)).unwrap();
    arena.push_raw(&make_record(0, 100, 0x100 | 99, 60, 50)).unwrap(); // secondary
    arena.push_raw(&make_record(0, 100, 0x800 | 99, 60, 50)).unwrap(); // supplementary

    assert_eq!(arena.len(), 3);

    let engine = PileupEngine::new(arena, 100, 149);
    let columns: Vec<_> = engine.collect();
    assert_eq!(columns[0].depth(), 3, "secondary+supplementary should appear by default");
}

// ---- bam.record.zero_refspan ----

// r[verify bam.record.zero_refspan]
#[test]
fn zero_refspan_endpos_equals_pos() {
    let raw = make_record_with_cigar(0, 100, 99, 60, &[cigar_op(10, 4)], 10); // 10S
    let mut arena = RecordStore::new();
    arena.push_raw(&raw).unwrap();
    let r = arena.record(0);
    assert_eq!(r.end_pos, 100, "zero-refspan should have end_pos == pos");
}

// ---- pileup.position_iteration (monotonic) ----

// r[verify pileup.position_iteration]
proptest! {
    #[test]
    fn column_positions_are_strictly_increasing(
        records in prop::collection::vec((0i32..200, 10u32..=50), 1..=15),
    ) {
        let mut sorted = records.clone();
        sorted.sort_by_key(|&(offset, _)| offset);

        let mut arena = RecordStore::new();
        for &(offset, read_len) in &sorted {
            arena.push_raw(&make_record(0, offset, 99, 60, read_len)).unwrap();
        }

        let engine = PileupEngine::new(arena, 0, 300);
        let mut prev_pos = -1i64;
        for col in engine {
            prop_assert!(col.pos() > prev_pos,
                "positions must be strictly increasing: {} not > {}", col.pos(), prev_pos);
            prev_pos = col.pos();
        }
    }
}

// ---- pileup.qpos (bounds check) ----

// r[verify pileup.qpos]
proptest! {
    #[test]
    fn qpos_within_read_length(
        start in 0i32..100,
        len in 10u32..=100,
    ) {
        let mut arena = RecordStore::new();
        // make_record uses a simple NxM CIGAR — no deletions, so all alignments have qpos.
        arena.push_raw(&make_record(0, start, 99, 60, len)).unwrap();

        let start_i64 = i64::from(start);
        let end_i64 = start_i64 + i64::from(len) - 1;
        let engine = PileupEngine::new(arena, start_i64, end_i64);

        for col in engine {
            for aln in col.alignments() {
                if let Some(qpos) = aln.qpos() {
                    prop_assert!(qpos < len as usize,
                        "qpos {} must be < seq_len {} at pos {}", qpos, len, col.pos());
                }
            }
        }
    }
}

// ---- pileup.column_contents (base and qual present) ----

// r[verify pileup.column_contents]
proptest! {
    #[test]
    fn all_alignments_have_valid_base_and_qual(
        n in 1usize..=10,
        len in 10u32..=50,
    ) {
        let mut arena = RecordStore::new();
        // make_record uses simple NxM CIGAR — all alignments are Match, so all have base/qual.
        for _ in 0..n {
            arena.push_raw(&make_record(0, 0, 99, 60, len)).unwrap();
        }

        let engine = PileupEngine::new(arena, 0, i64::from(len) - 1);
        for col in engine {
            prop_assert!(col.depth() > 0);
            for aln in col.alignments() {
                // Simple M CIGAR — only Match ops, so base and qual are always present.
                if let (Some(base), Some(qual)) = (aln.base(), aln.qual()) {
                    // Base must be a valid Base discriminant (A/C/G/T/Unknown)
                    let base_byte = base as u8;
                    prop_assert!(
                        matches!(base_byte, 65 | 67 | 71 | 84 | 78),
                        "invalid base byte {} at pos {}", base_byte, col.pos()
                    );
                    // Qual was set to 30 in make_record
                    prop_assert_eq!(qual, 30, "qual should be 30");
                }
            }
        }
    }
}

// ---- pileup.trailing_empty_termination ----

// r[verify pileup.trailing_empty_termination]
proptest! {
    #[test]
    fn trailing_empty_positions_terminate_early(
        start in 0i32..100,
        len in 10u32..=30,
        trailing in 100i64..=500,
    ) {
        let mut arena = RecordStore::new();
        arena.push_raw(&make_record(0, start, 99, 60, len)).unwrap();

        let start_i64 = i64::from(start);
        let region_end = start_i64 + i64::from(len) + trailing;
        let engine = PileupEngine::new(arena, start_i64, region_end);
        let columns: Vec<_> = engine.collect();

        // Should produce exactly `len` columns, not `len + trailing`
        prop_assert_eq!(columns.len(), len as usize,
            "should terminate after last record, not iterate {} trailing positions", trailing);
    }
}

// ---- reference_base propagation ----

// r[verify pileup.column_contents]
#[test]
fn reference_base_matches_ref_seq() {
    let ref_bases: Vec<Base> = vec![
        Base::A,
        Base::C,
        Base::G,
        Base::T,
        Base::A,
        Base::C,
        Base::G,
        Base::T,
        Base::A,
        Base::C,
    ];
    let ref_seq = RefSeq::new(Rc::from(ref_bases.as_slice()), 100);

    let mut arena = RecordStore::new();
    arena.push_raw(&make_record(0, 100, 99, 60, 10)).unwrap();

    let mut engine = PileupEngine::new(arena, 100, 109);
    engine.set_reference_seq(ref_seq);
    let columns: Vec<_> = engine.collect();

    assert_eq!(columns.len(), 10);
    for (i, col) in columns.iter().enumerate() {
        assert_eq!(
            col.reference_base(),
            ref_bases[i],
            "reference_base mismatch at pos {}",
            col.pos()
        );
    }
}

// r[verify pileup.column_contents]
#[test]
fn reference_base_unknown_without_ref_seq() {
    let mut arena = RecordStore::new();
    arena.push_raw(&make_record(0, 0, 99, 60, 5)).unwrap();

    let engine = PileupEngine::new(arena, 0, 4);
    for col in engine {
        assert_eq!(col.reference_base(), Base::Unknown, "should be Unknown when no ref_seq set");
    }
}

// ---- combined filter + dedup + max_depth ----

// r[verify pileup.read_filter]
// r[verify dedup.per_position]
// r[verify pileup.max_depth]
#[test]
fn filter_dedup_max_depth_combined() {
    use helpers::{BASE_A, make_named_record, pack_bases};

    let mut arena = RecordStore::new();
    let first: u16 = 0x41; // paired + first_in_template
    let second: u16 = 0x81; // paired + second_in_template

    let push = |arena: &mut RecordStore, name: &[u8], flags: u16| {
        let packed = vec![pack_bases(BASE_A, BASE_A); 10];
        let raw = make_named_record(name, 0, 0, flags, 60, 20, &packed);
        arena.push_raw(&raw).unwrap();
    };

    // 3 proper pairs (6 reads)
    for i in 0..3u8 {
        let name = [b'r', b'0' + i];
        push(&mut arena, &name, first);
        push(&mut arena, &name, second);
    }
    // 2 reads that will be filtered (secondary flag 0x100)
    arena.push_raw(&make_record(0, 0, 0x100 | 99, 60, 20)).unwrap();
    arena.push_raw(&make_record(0, 0, 0x100 | 99, 60, 20)).unwrap();

    let mut engine = PileupEngine::new(arena, 0, 19);
    engine.set_filter(|flags, _| flags & 0x100 == 0); // exclude secondary
    engine.set_dedup_overlapping();
    engine.set_max_depth(2);

    let columns: Vec<_> = engine.collect();
    let col = columns.first().unwrap();
    // 8 total → filter removes 2 secondary → 6 remain → dedup 3 pairs → 3 → max_depth caps to 2
    assert_eq!(col.depth(), 2, "filter→dedup→max_depth should yield 2");
}

// ---- take_store reuse ----

#[test]
fn take_store_reuse_across_regions() {
    let mut arena = RecordStore::new();
    arena.push_raw(&make_record(0, 0, 99, 60, 10)).unwrap();

    // Drive the engine to completion without consuming it via collect().
    let mut engine = PileupEngine::new(arena, 0, 9);
    let mut count = 0usize;
    while engine.next().is_some() {
        count += 1;
    }
    assert_eq!(count, 10);

    let mut store = engine.take_store().expect("store should be available");
    assert!(store.is_empty(), "store should be cleared after take");
    assert!(store.records_capacity() > 0, "capacity should be retained");

    store.push_raw(&make_record(0, 100, 99, 60, 5)).unwrap();
    let engine2 = PileupEngine::new(store, 100, 104);
    let columns2: Vec<_> = engine2.collect();
    assert_eq!(columns2.len(), 5);
    assert_eq!(columns2[0].pos(), 100);
}

#[test]
fn take_store_returns_none_after_second_take() {
    let mut arena = RecordStore::new();
    arena.push_raw(&make_record(0, 0, 99, 60, 5)).unwrap();

    let mut engine = PileupEngine::new(arena, 0, 4);
    while engine.next().is_some() {}

    let store = engine.take_store();
    assert!(store.is_some());

    let store2 = engine.take_store();
    assert!(store2.is_none(), "second take_store should return None");
}

// ---- region boundary behavior ----

// r[verify pileup.active_set]
proptest! {
    #[test]
    fn read_starting_before_region_contributes_within(
        read_start in 0i32..50,
        read_len in 60u32..=100,
        region_start_offset in 10i64..=40,
    ) {
        let mut arena = RecordStore::new();
        arena.push_raw(&make_record(0, read_start, 99, 60, read_len)).unwrap();

        let region_start = i64::from(read_start) + region_start_offset;
        let region_end = i64::from(read_start) + i64::from(read_len) - 1;
        let engine = PileupEngine::new(arena, region_start, region_end);
        let columns: Vec<_> = engine.collect();

        if !columns.is_empty() {
            prop_assert_eq!(columns[0].pos(), region_start,
                "first column should be at region_start");
        }
        for col in &columns {
            prop_assert!(col.pos() >= region_start, "no column before region_start");
            prop_assert!(col.pos() <= region_end, "no column after region_end");
        }
    }
}

// r[verify pileup.active_set]
#[test]
fn read_extending_past_region_end_is_truncated() {
    let mut arena = RecordStore::new();
    arena.push_raw(&make_record(0, 0, 99, 60, 100)).unwrap();

    // Region only covers positions 0-9, but read covers 0-99
    let engine = PileupEngine::new(arena, 0, 9);
    let columns: Vec<_> = engine.collect();
    assert_eq!(columns.len(), 10);
    assert_eq!(columns.last().unwrap().pos(), 9, "last column at region_end");
}

// ---- unmapped exclusion with complex reads ----

// r[verify pileup.unmapped_excluded]
proptest! {
    #[test]
    fn unmapped_reads_never_appear_in_columns(
        n_mapped in 1usize..=8,
        n_unmapped in 1usize..=5,
    ) {
        let mut arena = RecordStore::new();
        for _ in 0..n_mapped {
            arena.push_raw(&make_record(0, 0, 99, 60, 20)).unwrap();
        }
        for _ in 0..n_unmapped {
            arena.push_raw(&make_record(0, 0, 0x4, 60, 20)).unwrap(); // unmapped flag
        }

        let engine = PileupEngine::new(arena, 0, 19);
        for col in engine {
            prop_assert_eq!(col.depth(), n_mapped,
                "unmapped reads should be excluded, expected {} at pos {}", n_mapped, col.pos());
        }
    }
}

// ---- Complex CIGAR pileup tests using synthetic read strategy ----

use helpers::{arb_read, arb_read_set};

// r[verify pileup.active_set]
// r[verify pileup.qpos]
proptest! {
    #![proptest_config(ProptestConfig::with_cases(200))]

    /// Every column's match-depth must match the number of reads with M/=/X at that position.
    /// Deletions/RefSkips are now included in the pileup but not in covered_ref_positions.
    #[test]
    fn depth_matches_cigar_derived_coverage(reads in arb_read_set(15)) {
        let mut arena = RecordStore::new();
        for read in &reads {
            arena.push_raw(&read.raw).unwrap();
        }

        let region_start = reads.iter().map(|r| i64::from(r.pos)).min().unwrap();
        let region_end = reads.iter()
            .map(|r| i64::from(r.pos) + i64::from(r.ref_span) + 10)
            .max()
            .unwrap();

        let engine = PileupEngine::new(arena, region_start, region_end);
        for col in engine {
            // Only count M/=/X alignments (those with a qpos) — D/N are now included
            // in depth too, so compare match-depth against the oracle's M/=/X count.
            let match_depth = col.alignments().filter(|a| a.qpos().is_some()).count();
            let expected_match_depth = reads.iter()
                .filter(|r| r.covered_ref_positions().contains(&col.pos()))
                .count();
            prop_assert_eq!(match_depth, expected_match_depth,
                "match-depth mismatch at pos {} (expected from CIGAR analysis)", col.pos());
        }
    }

    /// qpos values from the pileup engine must match our independent CIGAR walk.
    #[test]
    fn qpos_matches_independent_cigar_walk(reads in arb_read_set(10)) {
        let mut arena = RecordStore::new();
        for read in &reads {
            arena.push_raw(&read.raw).unwrap();
        }

        let region_start = reads.iter().map(|r| i64::from(r.pos)).min().unwrap();
        let region_end = reads.iter()
            .map(|r| i64::from(r.pos) + i64::from(r.ref_span) + 10)
            .max()
            .unwrap();

        let engine = PileupEngine::new(arena, region_start, region_end);
        for col in engine {
            for aln in col.alignments() {
                let read = &reads[aln.record_idx() as usize];
                let expected_qpos = read.qpos_at(col.pos());
                prop_assert_eq!(aln.qpos(), expected_qpos,
                    "qpos mismatch at pos {} for read at pos {} with cigar {:?}",
                    col.pos(), read.pos, read.cigar_ops);
            }
        }
    }

    /// Deletions and ref-skips are now included in pileup columns.
    /// D/N positions have depth 1 with a Deletion or RefSkip op (no qpos).
    #[test]
    fn deletions_cause_absent_alignment(read in arb_read()) {
        let has_deletion = read.cigar_ops.iter().any(|&(_, op)| op == 2 || op == 3);
        if !has_deletion {
            return Ok(());
        }

        let mut arena = RecordStore::new();
        arena.push_raw(&read.raw).unwrap();

        let region_start = i64::from(read.pos);
        let region_end = region_start + i64::from(read.ref_span);
        let engine = PileupEngine::new(arena, region_start, region_end);

        let covered = read.covered_ref_positions();
        for col in engine {
            // Every ref position now yields depth 1 — M/=/X positions have qpos,
            // D/N positions have Deletion or RefSkip op without qpos.
            prop_assert_eq!(col.depth(), 1,
                "position {} should have depth 1 (all ref-consuming ops contribute)", col.pos());
            let aln = col.alignments().next().unwrap();
            if covered.contains(&col.pos()) {
                prop_assert!(aln.qpos().is_some(),
                    "position {} is M/=/X, should have qpos", col.pos());
            } else {
                prop_assert!(aln.is_del() || aln.is_refskip(),
                    "position {} is D/N, should be Deletion or RefSkip, got {:?}", col.pos(), aln.op);
                prop_assert!(aln.qpos().is_none(),
                    "position {} is D/N, should have no qpos", col.pos());
            }
        }
    }

    /// Insertions don't affect depth or ref positions — they only shift qpos.
    #[test]
    fn insertions_shift_qpos_but_not_depth(read in arb_read()) {
        let has_insertion = read.cigar_ops.iter().any(|&(_, op)| op == 1);
        if !has_insertion {
            return Ok(());
        }

        let mut arena = RecordStore::new();
        arena.push_raw(&read.raw).unwrap();

        let region_start = i64::from(read.pos);
        let region_end = region_start + i64::from(read.ref_span);
        let engine = PileupEngine::new(arena, region_start, region_end);

        for col in engine {
            let expected_qpos = read.qpos_at(col.pos());
            if let Some(expected) = expected_qpos {
                prop_assert_eq!(col.depth(), 1);
                let aln = col.alignments().next().unwrap();
                prop_assert_eq!(aln.qpos(), Some(expected),
                    "qpos wrong at pos {} — insertion should shift query offset", col.pos());
            }
        }
    }

    /// With multiple complex-CIGAR reads, no column should ever appear with
    /// zero depth (the engine skips empty positions).
    #[test]
    fn no_zero_depth_columns_with_complex_cigars(reads in arb_read_set(8)) {
        let mut arena = RecordStore::new();
        for read in &reads {
            arena.push_raw(&read.raw).unwrap();
        }

        let region_start = reads.iter().map(|r| i64::from(r.pos)).min().unwrap();
        let region_end = reads.iter()
            .map(|r| i64::from(r.pos) + i64::from(r.ref_span) + 100)
            .max()
            .unwrap();

        let engine = PileupEngine::new(arena, region_start, region_end);
        for col in engine {
            prop_assert!(col.depth() > 0,
                "empty column at pos {} with complex CIGARs", col.pos());
        }
    }

    /// qpos must always be within [0, seq_len) for every alignment that has one.
    /// Deletion and RefSkip alignments have no qpos — that is correct behavior.
    #[test]
    fn qpos_always_within_seq_bounds(reads in arb_read_set(10)) {
        let mut arena = RecordStore::new();
        for read in &reads {
            arena.push_raw(&read.raw).unwrap();
        }

        let region_start = reads.iter().map(|r| i64::from(r.pos)).min().unwrap();
        let region_end = reads.iter()
            .map(|r| i64::from(r.pos) + i64::from(r.ref_span) + 10)
            .max()
            .unwrap();

        let engine = PileupEngine::new(arena, region_start, region_end);
        for col in engine {
            for aln in col.alignments() {
                if let Some(qpos) = aln.qpos() {
                    prop_assert!(qpos < aln.seq_len as usize,
                        "qpos {} >= seq_len {} at pos {} for read with cigar {:?}",
                        qpos, aln.seq_len, col.pos(),
                        reads[aln.record_idx() as usize].cigar_ops);
                }
            }
        }
    }
}
