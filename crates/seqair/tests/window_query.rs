//! The window query from outside the crate: real reads through the reader's
//! own push path, `pileup_with` and `records_overlapping` composed the way a
//! caller's slow path composes them, and nothing but the public API.
//!
//! Covers `record_store.window_query`, `pileup.records_overlapping` and
//! `unified.readers_pileup_store_mutation`.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, reason = "test code")]
#![allow(
    clippy::cast_possible_truncation,
    clippy::arithmetic_side_effects,
    reason = "test code with known small values"
)]

mod helpers;
use helpers::ri;
use seqair::bam::cigar::{CigarOp, CigarOpType};
use seqair::bam::record_store::{RecordIdx, RecordStore};
use seqair::reader::{DepthLimit, Readers, Segment, SegmentOptions};
use seqair_types::{Base, Pos0};
use std::collections::BTreeMap;
use std::num::NonZeroU32;
use std::path::Path;

fn open() -> Readers {
    Readers::open(
        Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.bam")),
        Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.fasta.gz")),
    )
    .unwrap()
}

fn segment(readers: &Readers) -> Segment {
    let opts = SegmentOptions::new(NonZeroU32::new(5_000).unwrap());
    readers
        .segments(("chr19", Pos0::new(6_103_500).unwrap(), Pos0::new(6_106_500).unwrap()), opts)
        .unwrap()
        .next()
        .expect("the fixture yields a segment")
}

/// The definition, over the store's public iterator: no search, no index.
fn brute_force(store: &RecordStore, start: Pos0, end: Pos0) -> Vec<RecordIdx> {
    if end < start {
        return Vec::new();
    }
    store
        .records()
        .enumerate()
        .filter(|(_, rec)| !rec.flags.is_unmapped() && rec.pos <= end && rec.end_pos >= start)
        .map(|(idx, _)| ri(u32::try_from(idx).unwrap()))
        .collect()
}

fn window(around: Pos0, flank: u32) -> (Pos0, Pos0) {
    (
        Pos0::new(around.as_u32().saturating_sub(flank)).unwrap(),
        Pos0::new(around.as_u32() + flank).unwrap(),
    )
}

// r[verify record_store.window_query]
// r[verify pileup.records_overlapping]
/// Real reads through the reader: at every sampled column the query from the
/// column matches the definition, every alignment the column reports is in
/// the answer, and every answer that covers the column resolves on it. After
/// the last column the engine answers the same.
#[test]
fn query_matches_brute_force_on_real_reads() {
    let mut readers = open();
    let segment = segment(&readers);
    let mut pileup = readers.pileup(&segment, DepthLimit::Unlimited).unwrap();

    let mut sampled = 0usize;
    let mut non_empty = 0usize;
    while let Some(col) = pileup.pileups() {
        if col.pos().as_u32() % 97 != 0 {
            continue;
        }
        let (start, end) = window(col.pos(), 150);
        let hits: Vec<RecordIdx> = col.records_overlapping(start, end).collect();
        assert_eq!(hits, brute_force(col.store(), start, end), "window around {}", col.pos());
        for view in col.alignments() {
            assert!(
                hits.binary_search(&view.record_idx()).is_ok(),
                "column {} reports record {} the window around it does not",
                col.pos(),
                view.record_idx()
            );
        }
        for &idx in &hits {
            let rec = col.store().record(idx);
            let covers = rec.pos <= col.pos() && rec.end_pos >= col.pos();
            assert_eq!(col.find_record(idx).is_some(), covers, "record {idx} at {}", col.pos());
        }
        sampled += 1;
        non_empty += usize::from(!hits.is_empty());
    }
    assert!(sampled > 5, "the fixture segment should give several sampled columns");
    assert!(non_empty > 0, "some window should hold reads");

    let (start, end) = window(segment.start(), 300);
    let hits: Vec<RecordIdx> = pileup.records_overlapping(start, end).collect();
    assert_eq!(hits, brute_force(pileup.store(), start, end));
    assert!(!hits.is_empty());
}

/// Soft-clip the first aligned base of every read whose CIGAR starts with
/// `M` of length ≥ 2, shifting `pos` right by one; returns each moved read's
/// original position, keyed by qname and `end_pos` — which the move leaves
/// alone, and which tells mates of one template apart.
fn shift_leading_base(store: &mut RecordStore) -> BTreeMap<(Vec<u8>, Pos0), Pos0> {
    struct Planned {
        idx: RecordIdx,
        key: (Vec<u8>, Pos0),
        old_pos: Pos0,
        new_pos: Pos0,
        cigar: Vec<CigarOp>,
    }
    let mut plan: Vec<Planned> = Vec::new();
    for idx in store.indices() {
        let rec = store.record(idx);
        if rec.flags.is_unmapped() {
            continue;
        }
        let Ok(cigar) = rec.cigar(store) else { continue };
        let Some(first) = cigar.first() else { continue };
        if first.op_type() != CigarOpType::Match || first.len() < 2 {
            continue;
        }
        let mut new = vec![
            CigarOp::new(CigarOpType::SoftClip, 1),
            CigarOp::new(CigarOpType::Match, first.len() - 1),
        ];
        new.extend_from_slice(cigar.get(1..).unwrap_or_default());
        let new_pos = Pos0::new(rec.pos.as_u32() + 1).unwrap();
        let key = (rec.qname(store).unwrap().to_vec(), rec.end_pos);
        plan.push(Planned { idx, key, old_pos: rec.pos, new_pos, cigar: new });
    }
    let mut moved = BTreeMap::new();
    for Planned { idx, key, old_pos, new_pos, cigar } in plan {
        store.set_alignment(idx, new_pos, &cigar).expect("query length preserved");
        moved.insert(key, old_pos);
    }
    moved
}

// r[verify unified.readers_pileup_store_mutation+3]
// r[verify record_store.window_query]
/// The two features composed: a hook that moves reads and reads the
/// reference, then the query over the re-sorted store. The indices name the
/// moved records where they now lie, the hook's reference is the one the
/// columns report, and the answer still matches the definition.
#[test]
fn query_after_a_realigning_hook_names_the_moved_reads() {
    let mut readers = open();
    let segment = segment(&readers);
    let (seg_start, seg_end) = (segment.start(), segment.end());

    let mut moved = BTreeMap::new();
    let mut reference: Vec<Base> = Vec::new();
    let mut pileup = readers
        .pileup_with(&segment, DepthLimit::Unlimited, |store, ref_seq| {
            moved = shift_leading_base(store);
            assert_eq!(ref_seq.start_pos(), seg_start, "the hook sees the segment's reference");
            reference = (seg_start.as_u32()..=seg_end.as_u32())
                .map(|p| ref_seq.base_at(Pos0::new(p).unwrap()))
                .collect();
        })
        .unwrap();
    assert!(!moved.is_empty(), "the fixture has reads with a leading M");

    let hits: Vec<RecordIdx> = pileup.records_overlapping(seg_start, seg_end).collect();
    assert_eq!(hits, brute_force(pileup.store(), seg_start, seg_end));
    let mut found_moved = 0usize;
    for &idx in &hits {
        let rec = pileup.store().record(idx);
        let key = (rec.qname(pileup.store()).unwrap().to_vec(), rec.end_pos);
        if let Some(&old_pos) = moved.get(&key) {
            assert_eq!(
                rec.pos.as_u32(),
                old_pos.as_u32() + 1,
                "record {idx} sits where the hook put it"
            );
            found_moved += 1;
        }
    }
    assert!(found_moved > 0, "the query should name moved reads");

    while let Some(col) = pileup.pileups() {
        let at = (col.pos().as_u32() - seg_start.as_u32()) as usize;
        assert_eq!(reference.get(at).copied(), Some(col.reference_base()), "at {}", col.pos());
    }
}
