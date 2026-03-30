//! Tests for overlapping pair deduplication in the pileup engine.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::indexing_slicing)]
mod helpers;

use helpers::*;
use proptest::prelude::*;
use seqair::bam::{Pos, RecordStore, Zero, pileup::PileupEngine};

// Flags for first/second in pair
const FIRST: u16 = 0x41; // paired + first_in_template
const SECOND: u16 = 0x81; // paired + second_in_template

/// Push a named record with all-A sequence.
fn push_named(arena: &mut RecordStore, name: &[u8], pos: i32, flags: u16, seq_len: u32) {
    let packed_len = (seq_len as usize).div_ceil(2);
    let packed: Vec<u8> = vec![pack_bases(BASE_A, BASE_A); packed_len];
    let raw = make_named_record(name, 0, pos, flags, 60, seq_len, &packed);
    arena.push_raw(&raw).unwrap();
}

/// Push a named record with a specific base at position 0.
fn push_named_with_base(
    arena: &mut RecordStore,
    name: &[u8],
    pos: i32,
    flags: u16,
    seq_len: u32,
    first_base: u8,
) {
    let packed_len = (seq_len as usize).div_ceil(2);
    let mut packed = vec![pack_bases(BASE_A, BASE_A); packed_len];
    // Override first base
    packed[0] = pack_bases(first_base, packed[0] & 0x0F);
    let raw = make_named_record(name, 0, pos, flags, 60, seq_len, &packed);
    arena.push_raw(&raw).unwrap();
}

// ---- dedup.opt_in ----

// r[verify dedup.opt_in]
#[test]
fn dedup_disabled_by_default() {
    let mut arena = RecordStore::new();
    push_named(&mut arena, b"read1", 0, FIRST, 20);
    push_named(&mut arena, b"read1", 5, SECOND, 20);

    // Without dedup — both reads at overlapping positions
    let engine =
        PileupEngine::new(arena, Pos::<Zero>::new(0).unwrap(), Pos::<Zero>::new(24).unwrap());
    let columns: Vec<_> = engine.collect();
    let col = columns.iter().find(|c| c.pos() == Pos::<Zero>::new(10).unwrap()).unwrap();
    assert_eq!(col.depth(), 2, "without dedup, both mates should appear");
}

// r[verify dedup.opt_in]
#[test]
fn dedup_enabled_reduces_depth() {
    let mut arena = RecordStore::new();
    push_named(&mut arena, b"read1", 0, FIRST, 20);
    push_named(&mut arena, b"read1", 5, SECOND, 20);

    let mut engine =
        PileupEngine::new(arena, Pos::<Zero>::new(0).unwrap(), Pos::<Zero>::new(24).unwrap());
    engine.set_dedup_overlapping();
    let columns: Vec<_> = engine.collect();
    let col = columns.iter().find(|c| c.pos() == Pos::<Zero>::new(10).unwrap()).unwrap();
    assert_eq!(col.depth(), 1, "with dedup, only one mate at overlapping position");
}

// ---- dedup.mate_detection ----

// r[verify dedup.mate_detection]
#[test]
fn mates_detected_by_qname() {
    let mut arena = RecordStore::new();
    push_named(&mut arena, b"readA", 0, FIRST, 20);
    push_named(&mut arena, b"readA", 5, SECOND, 20);
    push_named(&mut arena, b"readB", 0, FIRST, 20);
    push_named(&mut arena, b"readB", 5, SECOND, 20);

    let mut engine =
        PileupEngine::new(arena, Pos::<Zero>::new(0).unwrap(), Pos::<Zero>::new(24).unwrap());
    engine.set_dedup_overlapping();
    let columns: Vec<_> = engine.collect();

    // At pos 10: 4 reads active, but 2 pairs → 2 after dedup
    let col = columns.iter().find(|c| c.pos() == Pos::<Zero>::new(10).unwrap()).unwrap();
    assert_eq!(col.depth(), 2, "two pairs should each keep one mate");
}

// r[verify dedup.mate_detection]
#[test]
fn unpaired_reads_unaffected() {
    let mut arena = RecordStore::new();
    push_named(&mut arena, b"solo1", 0, FIRST, 20);
    push_named(&mut arena, b"solo2", 0, FIRST, 20);
    push_named(&mut arena, b"solo3", 0, SECOND, 20);

    let mut engine =
        PileupEngine::new(arena, Pos::<Zero>::new(0).unwrap(), Pos::<Zero>::new(19).unwrap());
    engine.set_dedup_overlapping();
    let columns: Vec<_> = engine.collect();
    let col = columns.first().unwrap();
    assert_eq!(col.depth(), 3, "reads with unique names should all appear");
}

// ---- dedup.mate_pairs_only ----

// r[verify dedup.mate_pairs_only]
#[test]
fn third_record_with_same_name_not_deduped() {
    let mut arena = RecordStore::new();
    push_named(&mut arena, b"read1", 0, FIRST, 20);
    push_named(&mut arena, b"read1", 0, SECOND, 20);
    push_named(&mut arena, b"read1", 0, 0x800 | FIRST, 20); // supplementary

    let mut engine =
        PileupEngine::new(arena, Pos::<Zero>::new(0).unwrap(), Pos::<Zero>::new(19).unwrap());
    engine.set_dedup_overlapping();
    let columns: Vec<_> = engine.collect();
    let col = columns.first().unwrap();
    // First two are mates → deduplicated to 1. Third (supplementary) has no mate → kept.
    assert_eq!(col.depth(), 2, "supplementary with same name keeps one pair + supplementary");
}

// ---- dedup.resolution_same_base ----

// r[verify dedup.resolution_same_base]
#[test]
fn same_base_keeps_first_encountered() {
    let mut arena = RecordStore::new();
    // Both show A at pos 0
    push_named_with_base(&mut arena, b"read1", 0, FIRST, 10, BASE_A);
    push_named_with_base(&mut arena, b"read1", 0, SECOND, 10, BASE_A);

    let mut engine =
        PileupEngine::new(arena, Pos::<Zero>::new(0).unwrap(), Pos::<Zero>::new(0).unwrap());
    engine.set_dedup_overlapping();
    let columns: Vec<_> = engine.collect();
    assert_eq!(columns.len(), 1);
    assert_eq!(columns[0].depth(), 1);
    // The kept read should be the first-in-template
    let kept = columns[0].alignments().next().unwrap();
    assert!(kept.flags & 0x40 != 0);
}

// ---- dedup.resolution_different_base ----

// r[verify dedup.resolution_different_base]
#[test]
fn different_base_keeps_first_in_template() {
    let mut arena = RecordStore::new();
    // First-in-template shows A, second shows T
    push_named_with_base(&mut arena, b"read1", 0, FIRST, 10, BASE_A);
    push_named_with_base(&mut arena, b"read1", 0, SECOND, 10, BASE_T);

    let mut engine =
        PileupEngine::new(arena, Pos::<Zero>::new(0).unwrap(), Pos::<Zero>::new(0).unwrap());
    engine.set_dedup_overlapping();
    let columns: Vec<_> = engine.collect();
    assert_eq!(columns[0].depth(), 1);
    let kept = columns[0].alignments().next().unwrap();
    assert!(kept.flags & 0x40 != 0);
    assert_eq!(kept.base().unwrap() as u8, b'A');
}

// r[verify dedup.resolution_different_base]
#[test]
fn different_base_second_first_in_arena_order() {
    let mut arena = RecordStore::new();
    // Second-in-template pushed first into arena, shows T
    push_named_with_base(&mut arena, b"read1", 0, SECOND, 10, BASE_T);
    push_named_with_base(&mut arena, b"read1", 0, FIRST, 10, BASE_A);

    let mut engine =
        PileupEngine::new(arena, Pos::<Zero>::new(0).unwrap(), Pos::<Zero>::new(0).unwrap());
    engine.set_dedup_overlapping();
    let columns: Vec<_> = engine.collect();
    assert_eq!(columns[0].depth(), 1);
    // First-in-template should still be kept
    let kept = columns[0].alignments().next().unwrap();
    assert!(kept.flags & 0x40 != 0);
}

// ---- dedup.per_position ----

// r[verify dedup.per_position]
#[test]
fn both_mates_contribute_outside_overlap() {
    let mut arena = RecordStore::new();
    // mate1: pos 0-19, mate2: pos 10-29
    push_named(&mut arena, b"read1", 0, FIRST, 20);
    push_named(&mut arena, b"read1", 10, SECOND, 20);

    let mut engine =
        PileupEngine::new(arena, Pos::<Zero>::new(0).unwrap(), Pos::<Zero>::new(29).unwrap());
    engine.set_dedup_overlapping();
    let columns: Vec<_> = engine.collect();

    // pos 0-9: only mate1 → depth 1
    let col = columns.iter().find(|c| c.pos() == Pos::<Zero>::new(5).unwrap()).unwrap();
    assert_eq!(col.depth(), 1, "only mate1 at non-overlapping position");

    // pos 10-19: both mates overlap → depth 1 (deduped)
    let col = columns.iter().find(|c| c.pos() == Pos::<Zero>::new(15).unwrap()).unwrap();
    assert_eq!(col.depth(), 1, "deduped at overlapping position");

    // pos 20-29: only mate2 → depth 1
    let col = columns.iter().find(|c| c.pos() == Pos::<Zero>::new(25).unwrap()).unwrap();
    assert_eq!(col.depth(), 1, "only mate2 at non-overlapping position");

    // Total columns should be 30 (0-29, all covered)
    assert_eq!(columns.len(), 30);
}

// r[verify dedup.per_position]
proptest! {
    #[test]
    fn dedup_never_removes_both_mates(
        len in 20u32..=50,
        overlap in 1i32..=19,
    ) {
        let mut arena = RecordStore::new();
        let mate2_start = len as i32 - overlap;
        push_named(&mut arena, b"readX", 0, FIRST, len);
        push_named(&mut arena, b"readX", mate2_start, SECOND, len);

        let mut engine = PileupEngine::new(arena, Pos::<Zero>::new(0).unwrap(), Pos::<Zero>::new((mate2_start + len as i32 - 1) as u32).unwrap());
        engine.set_dedup_overlapping();
        let columns: Vec<_> = engine.collect();

        // Every position should have exactly depth 1 (never 0)
        for col in &columns {
            prop_assert_eq!(col.depth(), 1,
                "depth should be 1 at pos {}, not 0 (both mates removed)", col.pos().get());
        }
    }
}

// ---- dedup.filter_independent ----

// r[verify dedup.filter_independent]
#[test]
fn filtered_mate_does_not_trigger_dedup() {
    let mut arena = RecordStore::new();
    push_named(&mut arena, b"read1", 0, FIRST, 20);
    // Push second mate with duplicate flag so it gets filtered
    let packed = vec![pack_bases(BASE_A, BASE_A); 10];
    let raw = make_named_record(b"read1", 0, 0, SECOND | 0x400, 60, 20, &packed);
    arena.push_raw(&raw).unwrap();

    let mut engine =
        PileupEngine::new(arena, Pos::<Zero>::new(0).unwrap(), Pos::<Zero>::new(19).unwrap());
    engine.set_filter(|flags, _aux| flags & 0x400 == 0);
    engine.set_dedup_overlapping();
    let columns: Vec<_> = engine.collect();

    // Second mate filtered out by flag → first mate should be kept, not deduped
    let col = columns.first().unwrap();
    assert_eq!(col.depth(), 1, "first mate should survive when second is filtered");
}

// ---- dedup.max_depth_independent ----

// r[verify dedup.max_depth_independent]
#[test]
fn dedup_applied_before_max_depth() {
    let mut arena = RecordStore::new();
    // 3 pairs all overlapping at pos 0
    for i in 0..3 {
        let name = format!("read{i}");
        push_named(&mut arena, name.as_bytes(), 0, FIRST, 10);
        push_named(&mut arena, name.as_bytes(), 0, SECOND, 10);
    }

    let mut engine =
        PileupEngine::new(arena, Pos::<Zero>::new(0).unwrap(), Pos::<Zero>::new(9).unwrap());
    engine.set_dedup_overlapping();
    engine.set_max_depth(5);
    let columns: Vec<_> = engine.collect();

    // 6 reads → 3 after dedup → below max_depth of 5
    let col = columns.first().unwrap();
    assert_eq!(col.depth(), 3, "3 pairs deduped to 3, which is below max_depth 5");
}

// ---- property tests ----

// r[verify dedup.per_position]
// r[verify dedup.resolution_same_base]
proptest! {
    #[test]
    fn dedup_depth_never_exceeds_undeduplicated(
        n_pairs in 1usize..=10,
        overlap in 1i32..=20,
    ) {
        let len = 30u32;
        let build_store = |store: &mut RecordStore| {
            for i in 0..n_pairs {
                let name = format!("r{i}");
                let mate2_start = len as i32 - overlap;
                push_named(store, name.as_bytes(), 0, FIRST, len);
                push_named(store, name.as_bytes(), mate2_start, SECOND, len);
            }
        };

        let end = len + (len - overlap as u32) - 1;

        // Without dedup
        let mut arena_raw = RecordStore::new();
        build_store(&mut arena_raw);
        let engine_raw = PileupEngine::new(arena_raw, Pos::<Zero>::new(0).unwrap(), Pos::<Zero>::new(end).unwrap());
        let raw_columns: Vec<_> = engine_raw.collect();

        // With dedup
        let mut arena_dedup = RecordStore::new();
        build_store(&mut arena_dedup);
        let mut engine_dedup = PileupEngine::new(arena_dedup, Pos::<Zero>::new(0).unwrap(), Pos::<Zero>::new(end).unwrap());
        engine_dedup.set_dedup_overlapping();
        let dedup_columns: Vec<_> = engine_dedup.collect();

        // Same positions
        prop_assert_eq!(raw_columns.len(), dedup_columns.len());

        for (raw, dedup) in raw_columns.iter().zip(dedup_columns.iter()) {
            prop_assert_eq!(raw.pos(), dedup.pos());
            prop_assert!(dedup.depth() <= raw.depth(),
                "dedup depth {} exceeds raw depth {} at pos {}",
                dedup.depth(), raw.depth(), raw.pos().get());
            prop_assert!(dedup.depth() > 0,
                "dedup should not remove all reads at pos {}", dedup.pos().get());
        }
    }
}

// ---- cram.edge.unknown_read_names ----

// r[verify cram.edge.unknown_read_names]
#[test]
fn star_qname_records_not_paired() {
    let mut arena = RecordStore::new();
    push_named(&mut arena, b"*", 0, FIRST, 20);
    push_named(&mut arena, b"*", 0, SECOND, 20);

    let mut engine =
        PileupEngine::new(arena, Pos::<Zero>::new(0).unwrap(), Pos::<Zero>::new(19).unwrap());
    engine.set_dedup_overlapping();
    let columns: Vec<_> = engine.collect();
    let col = columns.first().unwrap();
    assert_eq!(col.depth(), 2, "records with qname '*' should not be treated as mates");
}

// r[verify cram.edge.unknown_read_names]
#[test]
fn empty_qname_records_not_paired() {
    let mut arena = RecordStore::new();
    push_named(&mut arena, b"", 0, FIRST, 20);
    push_named(&mut arena, b"", 0, SECOND, 20);

    let mut engine =
        PileupEngine::new(arena, Pos::<Zero>::new(0).unwrap(), Pos::<Zero>::new(19).unwrap());
    engine.set_dedup_overlapping();
    let columns: Vec<_> = engine.collect();
    let col = columns.first().unwrap();
    assert_eq!(col.depth(), 2, "records with empty qname should not be treated as mates");
}

// r[verify cram.edge.unknown_read_names]
proptest! {
    #[test]
    fn star_qnames_mixed_with_named_reads(n_named in 1usize..=5, n_star in 1usize..=5) {
        let mut arena = RecordStore::new();
        for i in 0..n_named {
            let name = format!("read{i}");
            push_named(&mut arena, name.as_bytes(), 0, FIRST, 20);
            push_named(&mut arena, name.as_bytes(), 0, SECOND, 20);
        }
        for _ in 0..n_star {
            push_named(&mut arena, b"*", 0, FIRST, 20);
        }

        let mut engine = PileupEngine::new(arena, Pos::<Zero>::new(0).unwrap(), Pos::<Zero>::new(19).unwrap());
        engine.set_dedup_overlapping();
        let columns: Vec<_> = engine.collect();
        let col = columns.first().unwrap();
        // n_named pairs dedup to n_named, plus n_star unpaired star records
        prop_assert_eq!(col.depth(), n_named + n_star,
            "expected {} named + {} star = {}", n_named, n_star, n_named + n_star);
    }
}
