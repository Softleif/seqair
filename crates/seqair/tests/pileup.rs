//! Tests for PileupEngine: column iteration, filtering, max depth, edge cases.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::indexing_slicing)]
mod helpers;

use helpers::{cigar_op, make_record, make_record_with_cigar};
use proptest::prelude::*;
use seqair::bam::pileup::RefSeq;
use seqair::bam::{Pos, RecordStore, Zero, pileup::PileupEngine};
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

        let engine = PileupEngine::new(arena, Pos::<Zero>::new(0).unwrap(), Pos::<Zero>::new(600).unwrap());
        for col in engine {
            let pos = col.pos().as_usize();
            let exp = expected_depth.get(pos).copied().unwrap_or(0);
            prop_assert_eq!(col.depth(), exp,
                "depth mismatch at pos {}", col.pos().get());
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

        let engine = PileupEngine::new(arena, Pos::<Zero>::new(0).unwrap(), Pos::<Zero>::new(600).unwrap());
        let columns: Vec<_> = engine.collect();
        let col_positions: std::collections::HashSet<u32> =
            columns.iter().map(|c| c.pos().get()).collect();

        for pos in 0..MAX_POS as u32 {
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
    let mut engine =
        PileupEngine::new(arena, Pos::<Zero>::new(0).unwrap(), Pos::<Zero>::new(19).unwrap());
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
        let mut engine = PileupEngine::new(arena, Pos::<Zero>::new(0).unwrap(), Pos::<Zero>::new(0).unwrap());
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

        let mut engine = PileupEngine::new(arena, Pos::<Zero>::new(0).unwrap(), Pos::<Zero>::new(len - 1).unwrap());
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

        let mut engine = PileupEngine::new(arena, Pos::<Zero>::new(0).unwrap(), Pos::<Zero>::new(99).unwrap());
        engine.set_max_depth(3);
        let columns: Vec<_> = engine.collect();

        let col0 = columns.iter().find(|c| c.pos() == Pos::<Zero>::new(0).unwrap()).unwrap();
        prop_assert_eq!(col0.depth(), 3); // capped

        let col50 = columns.iter().find(|c| c.pos() == Pos::<Zero>::new(50).unwrap()).unwrap();
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

        let start_u32 = start as u32;
        let end_u32 = start_u32 + len - 1;
        let engine = PileupEngine::new(arena, Pos::<Zero>::new(start_u32).unwrap(), Pos::<Zero>::new(end_u32).unwrap());

        for col in engine {
            for aln in col.alignments() {
                prop_assert_eq!(aln.qpos().unwrap(), (col.pos().get() - start_u32) as usize);
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

        let end = pos2 as u32 + len2 + 10;
        let engine = PileupEngine::new(arena, Pos::<Zero>::new(pos1 as u32).unwrap(), Pos::<Zero>::new(end).unwrap());
        for col in engine {
            prop_assert!(col.depth() > 0, "empty column at {}", col.pos().get());
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

    let engine =
        PileupEngine::new(arena, Pos::<Zero>::new(100).unwrap(), Pos::<Zero>::new(149).unwrap());
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
    assert_eq!(arena.record(0).end_pos, Pos::<Zero>::new(50).unwrap());

    let engine =
        PileupEngine::new(arena, Pos::<Zero>::new(48).unwrap(), Pos::<Zero>::new(52).unwrap());
    let columns: Vec<_> = engine.collect();
    // Pure soft-clip has no ref-consuming ops → no qpos → no column
    assert!(
        columns.is_empty()
            || columns.iter().all(|c| c.pos() != Pos::<Zero>::new(50).unwrap() || c.depth() == 0)
    );
}

// ---- pileup.soft_clip_at_position ----

// r[verify pileup.soft_clip_at_position]
#[test]
fn leading_softclip_does_not_extend_ref_range() {
    // 5S 20M at pos 100 → ref 100-119, query 5-24
    let raw = make_record_with_cigar(0, 100, 99, 60, &[cigar_op(5, 4), cigar_op(20, 0)], 25);
    let mut arena = RecordStore::new();
    arena.push_raw(&raw).unwrap();

    let engine =
        PileupEngine::new(arena, Pos::<Zero>::new(95).unwrap(), Pos::<Zero>::new(125).unwrap());
    let columns: Vec<_> = engine.collect();
    let positions: Vec<u32> = columns.iter().map(|c| c.pos().get()).collect();

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

    let engine =
        PileupEngine::new(arena, Pos::<Zero>::new(100).unwrap(), Pos::<Zero>::new(149).unwrap());
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
    assert_eq!(
        r.end_pos,
        Pos::<Zero>::new(100).unwrap(),
        "zero-refspan should have end_pos == pos"
    );
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

        let engine = PileupEngine::new(arena, Pos::<Zero>::new(0).unwrap(), Pos::<Zero>::new(300).unwrap());
        let mut prev_pos: Option<u32> = None;
        for col in engine {
            let cur = col.pos().get();
            if let Some(prev) = prev_pos {
                prop_assert!(cur > prev,
                    "positions must be strictly increasing: {} not > {}", cur, prev);
            }
            prev_pos = Some(cur);
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

        let start_u32 = start as u32;
        let end_u32 = start_u32 + len - 1;
        let engine = PileupEngine::new(arena, Pos::<Zero>::new(start_u32).unwrap(), Pos::<Zero>::new(end_u32).unwrap());

        for col in engine {
            for aln in col.alignments() {
                if let Some(qpos) = aln.qpos() {
                    prop_assert!(qpos < len as usize,
                        "qpos {} must be < seq_len {} at pos {}", qpos, len, col.pos().get());
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

        let engine = PileupEngine::new(arena, Pos::<Zero>::new(0).unwrap(), Pos::<Zero>::new(len - 1).unwrap());
        for col in engine {
            prop_assert!(col.depth() > 0);
            for aln in col.alignments() {
                // Simple M CIGAR — only Match ops, so base and qual are always present.
                if let (Some(base), Some(qual)) = (aln.base(), aln.qual()) {
                    // Base must be a valid Base discriminant (A/C/G/T/Unknown)
                    let base_byte = base as u8;
                    prop_assert!(
                        matches!(base_byte, 65 | 67 | 71 | 84 | 78),
                        "invalid base byte {} at pos {}", base_byte, col.pos().get()
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
        trailing in 100u32..=500,
    ) {
        let mut arena = RecordStore::new();
        arena.push_raw(&make_record(0, start, 99, 60, len)).unwrap();

        let start_u32 = start as u32;
        let region_end = start_u32 + len + trailing;
        let engine = PileupEngine::new(arena, Pos::<Zero>::new(start_u32).unwrap(), Pos::<Zero>::new(region_end).unwrap());
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
    let ref_seq = RefSeq::new(Rc::from(ref_bases.as_slice()), Pos::<Zero>::new(100).unwrap());

    let mut arena = RecordStore::new();
    arena.push_raw(&make_record(0, 100, 99, 60, 10)).unwrap();

    let mut engine =
        PileupEngine::new(arena, Pos::<Zero>::new(100).unwrap(), Pos::<Zero>::new(109).unwrap());
    engine.set_reference_seq(ref_seq);
    let columns: Vec<_> = engine.collect();

    assert_eq!(columns.len(), 10);
    for (i, col) in columns.iter().enumerate() {
        assert_eq!(
            col.reference_base(),
            ref_bases[i],
            "reference_base mismatch at pos {}",
            col.pos().get()
        );
    }
}

// r[verify pileup.column_contents]
#[test]
fn reference_base_unknown_without_ref_seq() {
    let mut arena = RecordStore::new();
    arena.push_raw(&make_record(0, 0, 99, 60, 5)).unwrap();

    let engine =
        PileupEngine::new(arena, Pos::<Zero>::new(0).unwrap(), Pos::<Zero>::new(4).unwrap());
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

    let mut engine =
        PileupEngine::new(arena, Pos::<Zero>::new(0).unwrap(), Pos::<Zero>::new(19).unwrap());
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
    let mut engine =
        PileupEngine::new(arena, Pos::<Zero>::new(0).unwrap(), Pos::<Zero>::new(9).unwrap());
    let mut count = 0usize;
    while engine.next().is_some() {
        count += 1;
    }
    assert_eq!(count, 10);

    let mut store = engine.take_store().expect("store should be available");
    assert!(store.is_empty(), "store should be cleared after take");
    assert!(store.records_capacity() > 0, "capacity should be retained");

    store.push_raw(&make_record(0, 100, 99, 60, 5)).unwrap();
    let engine2 =
        PileupEngine::new(store, Pos::<Zero>::new(100).unwrap(), Pos::<Zero>::new(104).unwrap());
    let columns2: Vec<_> = engine2.collect();
    assert_eq!(columns2.len(), 5);
    assert_eq!(columns2[0].pos(), Pos::<Zero>::new(100).unwrap());
}

#[test]
fn take_store_returns_none_after_second_take() {
    let mut arena = RecordStore::new();
    arena.push_raw(&make_record(0, 0, 99, 60, 5)).unwrap();

    let mut engine =
        PileupEngine::new(arena, Pos::<Zero>::new(0).unwrap(), Pos::<Zero>::new(4).unwrap());
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
        region_start_offset in 10u32..=40,
    ) {
        let mut arena = RecordStore::new();
        arena.push_raw(&make_record(0, read_start, 99, 60, read_len)).unwrap();

        let region_start = read_start as u32 + region_start_offset;
        let region_end = read_start as u32 + read_len - 1;
        let engine = PileupEngine::new(arena, Pos::<Zero>::new(region_start).unwrap(), Pos::<Zero>::new(region_end).unwrap());
        let columns: Vec<_> = engine.collect();

        if !columns.is_empty() {
            prop_assert_eq!(columns[0].pos(), Pos::<Zero>::new(region_start).unwrap(),
                "first column should be at region_start");
        }
        for col in &columns {
            prop_assert!(col.pos() >= Pos::<Zero>::new(region_start).unwrap(), "no column before region_start");
            prop_assert!(col.pos() <= Pos::<Zero>::new(region_end).unwrap(), "no column after region_end");
        }
    }
}

// r[verify pileup.active_set]
#[test]
fn read_extending_past_region_end_is_truncated() {
    let mut arena = RecordStore::new();
    arena.push_raw(&make_record(0, 0, 99, 60, 100)).unwrap();

    // Region only covers positions 0-9, but read covers 0-99
    let engine =
        PileupEngine::new(arena, Pos::<Zero>::new(0).unwrap(), Pos::<Zero>::new(9).unwrap());
    let columns: Vec<_> = engine.collect();
    assert_eq!(columns.len(), 10);
    assert_eq!(
        columns.last().unwrap().pos(),
        Pos::<Zero>::new(9).unwrap(),
        "last column at region_end"
    );
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

        let engine = PileupEngine::new(arena, Pos::<Zero>::new(0).unwrap(), Pos::<Zero>::new(19).unwrap());
        for col in engine {
            prop_assert_eq!(col.depth(), n_mapped,
                "unmapped reads should be excluded, expected {} at pos {}", n_mapped, col.pos().get());
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

        let region_start = reads.iter().map(|r| r.pos as u32).min().unwrap();
        let region_end = reads.iter()
            .map(|r| r.pos as u32 + r.ref_span + 10)
            .max()
            .unwrap();

        let engine = PileupEngine::new(arena, Pos::<Zero>::new(region_start).unwrap(), Pos::<Zero>::new(region_end).unwrap());
        for col in engine {
            // Only count M/=/X alignments (those with a qpos) — D/N are now included
            // in depth too, so compare match-depth against the oracle's M/=/X count.
            let match_depth = col.alignments().filter(|a| a.qpos().is_some()).count();
            let expected_match_depth = reads.iter()
                .filter(|r| r.covered_ref_positions().contains(&col.pos().as_i64()))
                .count();
            prop_assert_eq!(match_depth, expected_match_depth,
                "match-depth mismatch at pos {} (expected from CIGAR analysis)", col.pos().get());
        }
    }

    /// qpos values from the pileup engine must match our independent CIGAR walk.
    #[test]
    fn qpos_matches_independent_cigar_walk(reads in arb_read_set(10)) {
        let mut arena = RecordStore::new();
        for read in &reads {
            arena.push_raw(&read.raw).unwrap();
        }

        let region_start = reads.iter().map(|r| r.pos as u32).min().unwrap();
        let region_end = reads.iter()
            .map(|r| r.pos as u32 + r.ref_span + 10)
            .max()
            .unwrap();

        let engine = PileupEngine::new(arena, Pos::<Zero>::new(region_start).unwrap(), Pos::<Zero>::new(region_end).unwrap());
        for col in engine {
            for aln in col.alignments() {
                let read = &reads[aln.record_idx() as usize];
                let expected_qpos = read.qpos_at(col.pos().as_i64());
                prop_assert_eq!(aln.qpos(), expected_qpos,
                    "qpos mismatch at pos {} for read at pos {} with cigar {:?}",
                    col.pos().get(), read.pos, read.cigar_ops);
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

        let region_start = read.pos as u32;
        let region_end = region_start + read.ref_span;
        let engine = PileupEngine::new(arena, Pos::<Zero>::new(region_start).unwrap(), Pos::<Zero>::new(region_end).unwrap());

        let covered = read.covered_ref_positions();
        for col in engine {
            // Every ref position now yields depth 1 — M/=/X positions have qpos,
            // D/N positions have Deletion or RefSkip op without qpos.
            prop_assert_eq!(col.depth(), 1,
                "position {} should have depth 1 (all ref-consuming ops contribute)", col.pos().get());
            let aln = col.alignments().next().unwrap();
            if covered.contains(&col.pos().as_i64()) {
                prop_assert!(aln.qpos().is_some(),
                    "position {} is M/=/X, should have qpos", col.pos().get());
            } else {
                prop_assert!(aln.is_del() || aln.is_refskip(),
                    "position {} is D/N, should be Deletion or RefSkip, got {:?}", col.pos().get(), aln.op);
                prop_assert!(aln.qpos().is_none(),
                    "position {} is D/N, should have no qpos", col.pos().get());
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

        let region_start = read.pos as u32;
        let region_end = region_start + read.ref_span;
        let engine = PileupEngine::new(arena, Pos::<Zero>::new(region_start).unwrap(), Pos::<Zero>::new(region_end).unwrap());

        for col in engine {
            let expected_qpos = read.qpos_at(col.pos().as_i64());
            if let Some(expected) = expected_qpos {
                prop_assert_eq!(col.depth(), 1);
                let aln = col.alignments().next().unwrap();
                prop_assert_eq!(aln.qpos(), Some(expected),
                    "qpos wrong at pos {} — insertion should shift query offset", col.pos().get());
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

        let region_start = reads.iter().map(|r| r.pos as u32).min().unwrap();
        let region_end = reads.iter()
            .map(|r| r.pos as u32 + r.ref_span + 100)
            .max()
            .unwrap();

        let engine = PileupEngine::new(arena, Pos::<Zero>::new(region_start).unwrap(), Pos::<Zero>::new(region_end).unwrap());
        for col in engine {
            prop_assert!(col.depth() > 0,
                "empty column at pos {} with complex CIGARs", col.pos().get());
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

        let region_start = reads.iter().map(|r| r.pos as u32).min().unwrap();
        let region_end = reads.iter()
            .map(|r| r.pos as u32 + r.ref_span + 10)
            .max()
            .unwrap();

        let engine = PileupEngine::new(arena, Pos::<Zero>::new(region_start).unwrap(), Pos::<Zero>::new(region_end).unwrap());
        for col in engine {
            for aln in col.alignments() {
                if let Some(qpos) = aln.qpos() {
                    prop_assert!(qpos < aln.seq_len as usize,
                        "qpos {} >= seq_len {} at pos {} for read with cigar {:?}",
                        qpos, aln.seq_len, col.pos().get(),
                        reads[aln.record_idx() as usize].cigar_ops);
                }
            }
        }
    }
}
