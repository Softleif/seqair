//! The two preconditions `PileupEngine` relies on, and the type that carries them.
//!
//! Covers `record_store.pileup_input`.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]
#![allow(clippy::cast_possible_truncation, reason = "test code with known small values")]

use proptest::prelude::*;
use seqair::bam::cigar::{CigarOp, CigarOpType};
use seqair::bam::pileup::PileupEngine;
use seqair::bam::record_store::RecordStore;
use seqair_types::{BamFlags, Base, Pos0};

const PAIRED: u16 = 0x1;
const FIRST: u16 = 0x40;
const SECOND: u16 = 0x80;

fn push(store: &mut RecordStore, qname: &[u8], pos: u32, len: u32, mate_pos: i32, flags: u16) {
    let cigar = [CigarOp::new(CigarOpType::Match, len)];
    store
        .push_fields(
            Pos0::new(pos).unwrap(),
            Pos0::new(pos + len - 1).unwrap(),
            BamFlags::from(flags),
            60,
            len,
            0,
            qname,
            &cigar,
            &vec![Base::A; len as usize],
            &vec![30u8; len as usize],
            &[],
            0,
            0,
            mate_pos,
            0,
            &mut (),
        )
        .expect("well-formed record")
        .expect("no customizer");
}

/// Depth at every column the engine yields, keyed by position.
fn depths(input: seqair::bam::record_store::PileupInput, end: u32) -> Vec<(u64, usize)> {
    let mut engine = PileupEngine::new(input, Pos0::new(0).unwrap(), Pos0::new(end).unwrap());
    let mut out = Vec::new();
    while let Some(col) = engine.pileups() {
        out.push((col.pos().as_u64(), col.depth()));
    }
    out
}

/// Records pushed out of position order must still all reach the pileup.
///
/// The engine walks the store assuming ascending `pos` and evicts on
/// `end_pos`, so an out-of-order record is silently skipped — a column comes
/// back one read short, and nothing anywhere reports it. `prepare_for_pileup`
/// is what makes that unreachable: it sorts, and the engine takes only its
/// output.
#[test]
fn out_of_order_pushes_still_reach_the_pileup() {
    let mut store = RecordStore::new();
    push(&mut store, b"late", 40, 10, -1, 0);
    push(&mut store, b"early", 10, 10, -1, 0); // out of order
    push(&mut store, b"mid", 20, 10, -1, 0); // out of order

    let observed = depths(store.prepare_for_pileup().input, 60);
    let covered: Vec<u64> = observed.iter().map(|&(pos, _)| pos).collect();

    assert!(covered.contains(&10), "the second-pushed read must be pileuped: {observed:?}");
    assert!(covered.contains(&20), "the third-pushed read must be pileuped: {observed:?}");
    assert!(covered.contains(&40), "the first-pushed read must be pileuped: {observed:?}");
    assert_eq!(observed.len(), 30, "three 10bp reads, no overlap");
}

/// A store reaching the engine always has its mates linked.
///
/// `mate_idx()` returning `None` is indistinguishable from "this record has no
/// mate", so a consumer that deduplicates overlapping mates on an unlinked
/// store silently keeps both — and a test asserting dedup behaviour on one
/// passes while proving nothing.
#[test]
fn a_store_reaching_the_engine_has_its_mates_linked() {
    let mut store = RecordStore::new();
    push(&mut store, b"frag", 10, 20, 20, PAIRED | FIRST);
    push(&mut store, b"frag", 20, 20, 10, PAIRED | SECOND);

    let prepared = store.prepare_for_pileup();
    assert_eq!(prepared.stats.pairs, 1, "preparing must link the pair");

    let mut engine =
        PileupEngine::new(prepared.input, Pos0::new(0).unwrap(), Pos0::new(60).unwrap());
    let mut checked = 0;
    while let Some(col) = engine.pileups() {
        if col.pos().as_u64() != 25 {
            continue;
        }
        for view in col.alignments() {
            assert!(view.mate_idx().is_some(), "both mates are in the store and must be linked");
            assert!(view.in_mate_overlap(), "position 25 is inside the pair's overlap");
            checked += 1;
        }
    }
    assert_eq!(checked, 2, "both mates cover position 25");
}

/// A read the property tests place: position and length only.
#[derive(Debug, Clone, Copy)]
struct Placed {
    pos: u32,
    len: u32,
}

fn placed() -> impl Strategy<Value = Placed> {
    (0u32..200, 1u32..40).prop_map(|(pos, len)| Placed { pos, len })
}

/// Column depths for `reads` pushed in the given order, after preparation.
fn profile(reads: &[Placed]) -> Vec<(u64, usize)> {
    let mut store = RecordStore::new();
    for (i, r) in reads.iter().enumerate() {
        push(&mut store, format!("r{i}").as_bytes(), r.pos, r.len, -1, 0);
    }
    depths(store.prepare_for_pileup().input, 260)
}

proptest! {
    /// A pileup is a function of the record *set*, not of the order the records
    /// were pushed in. The engine only walks ascending positions, so before
    /// `PileupInput` this held exactly when the caller remembered to sort — and
    /// failed silently otherwise, dropping every record that arrived early.
    #[test]
    fn the_pileup_does_not_depend_on_push_order(
        reads in prop::collection::vec(placed(), 1..12),
        rotation in 0usize..12,
    ) {
        let mut sorted = reads.clone();
        sorted.sort_by_key(|r| r.pos);

        let mut rotated = reads.clone();
        let n = rotated.len();
        rotated.rotate_left(rotation % n);

        let expected = profile(&sorted);
        prop_assert_eq!(profile(&reads), expected.clone(), "as generated");
        prop_assert_eq!(profile(&rotated), expected, "rotated");
    }

    /// Total observations are conserved: every read contributes exactly its
    /// length in covered columns, whatever order it arrived in. This is the
    /// property the old behaviour broke outright — reads pushed early produced
    /// no columns at all.
    #[test]
    fn every_pushed_read_is_wholly_pileuped(reads in prop::collection::vec(placed(), 1..12)) {
        let expected: usize = reads.iter().map(|r| r.len as usize).sum();
        let observed: usize = profile(&reads).iter().map(|&(_, depth)| depth).sum();
        prop_assert_eq!(observed, expected);
    }

    /// Preparing a store whose properties already hold must be a no-op.
    ///
    /// `prepare_for_pileup` skips the sort when no push arrived out of order
    /// and skips linking when nothing invalidated it, so a store the caller
    /// already sorted and linked by hand takes a different path through it than
    /// a raw one. Both must reach the same pileup.
    #[test]
    fn preparing_an_already_prepared_store_is_a_no_op(
        reads in prop::collection::vec(placed(), 1..8),
    ) {
        let build = || {
            let mut store = RecordStore::new();
            for (i, r) in reads.iter().enumerate() {
                push(&mut store, format!("r{i}").as_bytes(), r.pos, r.len, -1, 0);
            }
            store
        };
        let mut warm = build();
        warm.sort_by_pos();
        let _ = warm.link_mates();

        prop_assert_eq!(
            depths(warm.prepare_for_pileup().input, 260),
            depths(build().prepare_for_pileup().input, 260)
        );
    }
}

/// `dedup` must establish the order it needs rather than document it.
///
/// It collapses *consecutive* equal records, so duplicates are only adjacent
/// once the store is sorted — which its doc comment used to require of the
/// caller ("Must be called after `sort_by_pos`"). An unsorted store silently
/// kept the duplicates: no error, just a record count that is wrong in the one
/// case the method exists for (overlapping BAM index chunks loading a record
/// twice).
#[test]
fn dedup_collapses_duplicates_that_arrive_out_of_order() {
    let mut store = RecordStore::new();
    push(&mut store, b"dup", 10, 10, -1, 0);
    push(&mut store, b"other", 40, 10, -1, 0); // separates the duplicates
    push(&mut store, b"dup", 10, 10, -1, 0);

    store.dedup();
    assert_eq!(store.len(), 2, "the two identical records must collapse to one");
}

proptest! {
    /// Whatever order distinct records and their duplicates arrive in, `dedup`
    /// leaves exactly the distinct ones.
    ///
    /// Positions are made distinct per record so that "duplicate" and
    /// "adjacent after sorting" coincide, which is the contract `dedup`
    /// documents; it is deliberately not a full-store uniqueness pass.
    #[test]
    fn dedup_leaves_exactly_the_distinct_records(
        count in 1usize..8,
        copies in prop::collection::vec(1usize..3, 1..8),
        rotation in 0usize..8,
    ) {
        let mut planned: Vec<usize> = Vec::new();
        for i in 0..count {
            let n = copies.get(i).copied().unwrap_or(1);
            planned.extend(std::iter::repeat_n(i, n));
        }
        let n = planned.len();
        planned.rotate_left(rotation % n);

        let mut store = RecordStore::new();
        for &i in &planned {
            // Distinct position per record, so duplicates and only duplicates
            // become adjacent once sorted.
            push(&mut store, format!("r{i}").as_bytes(), 10 + i as u32 * 20, 10, -1, 0);
        }
        store.dedup();
        prop_assert_eq!(store.len(), count);
    }
}
