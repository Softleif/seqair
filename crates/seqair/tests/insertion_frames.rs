//! The two query frames an insertion is reported in, pinned against each other.
//!
//! Covers `cigar.aligned_pairs.insertion_qpos` and
//! `pileup_indel.insertion_at_last_match`.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]
#![allow(clippy::arithmetic_side_effects, reason = "test code")]

mod helpers;
use hegel::prelude::*;
use helpers::ri;
use seqair::bam::aligned_pairs::AlignedPair;
use seqair::bam::cigar::{CigarOp, CigarOpType};
use seqair::bam::pileup::{PileupEngine, PileupOp};
use seqair::bam::record_store::RecordStore;
use seqair_types::{BamFlags, Base, Pos0};

/// `<lead>M <ins>I <trail>M` at `pos`, sequence all `A` except the inserted run.
fn store_with_insertion(pos: u32, lead: u32, ins: u32, trail: u32) -> RecordStore {
    let mut seq = vec![Base::A; lead as usize];
    seq.extend(vec![Base::C; ins as usize]);
    seq.extend(vec![Base::A; trail as usize]);
    let cigar = [
        CigarOp::new(CigarOpType::Match, lead),
        CigarOp::new(CigarOpType::Insertion, ins),
        CigarOp::new(CigarOpType::Match, trail),
    ];
    let mut store = RecordStore::new();
    store
        .push_fields(
            Pos0::new(pos).unwrap(),
            Pos0::new(pos + lead + trail - 1).unwrap(),
            BamFlags::from(0u16),
            60,
            lead + trail,
            ins,
            b"ins",
            &cigar,
            &seq,
            &vec![30u8; seq.len()],
            &[],
            0,
            0,
            -1,
            0,
            &mut (),
        )
        .unwrap()
        .unwrap();
    store
}

/// The aligned-pairs frame: the first inserted base.
fn first_inserted_from_pairs(store: &RecordStore) -> u32 {
    let rec = store.record(ri(0)).unwrap();
    let pairs = rec.aligned_pairs(store).expect("valid cigar");
    pairs
        .filter_map(|p| match p {
            AlignedPair::Insertion { first_inserted, .. } => Some(first_inserted.get()),
            _ => None,
        })
        .next()
        .expect("the read has an insertion")
}

/// The pileup frame: the last matched base *before* the insertion.
fn anchor_from_pileup(store: RecordStore, region_end: u32) -> u32 {
    let mut engine = PileupEngine::new(
        store.prepare_for_pileup().input,
        Pos0::new(0).unwrap(),
        Pos0::new(region_end).unwrap(),
    );
    let mut anchor = None;
    while anchor.is_none() {
        let Some(col) = engine.pileups() else {
            break;
        };
        for view in col.alignments() {
            if let PileupOp::Insertion { qpos, .. } = view.op() {
                anchor = Some(qpos.get());
                break;
            }
        }
    }
    anchor.expect("no insertion column")
}

/// One insertion, two frames, one base apart — always.
///
/// `AlignedPair::Insertion` reports the first *inserted* base;
/// `PileupOp::Insertion` reports the matched base that *precedes* the run,
/// because a pileup column has to sit on a reference position and an
/// insertion consumes none. Nothing pinned the relationship between them,
/// so a consumer moving between the two APIs had only prose to go on.
#[hegel::test]
fn the_pileup_anchor_is_one_before_the_first_inserted_base(tc: TestCase) {
    let pos = tc.draw(gs::integers::<u32>().max_value(49));
    let lead = tc.draw(gs::integers::<u32>().min_value(1).max_value(19));
    let ins = tc.draw(gs::integers::<u32>().min_value(1).max_value(7));
    let trail = tc.draw(gs::integers::<u32>().min_value(1).max_value(19));
    let store = store_with_insertion(pos, lead, ins, trail);
    let first_inserted = first_inserted_from_pairs(&store);
    let anchor = anchor_from_pileup(store, pos + lead + trail + 10);

    assert_eq!(first_inserted, anchor + 1, "the frames differ by exactly one base");
    assert_eq!(anchor, lead - 1, "the anchor is the last leading matched base");
}

/// The inserted run the pileup exposes starts at the aligned-pairs frame,
/// so a caller never has to reconstruct it from the anchor.
#[hegel::test]
fn inserted_bases_start_at_the_aligned_pairs_frame(tc: TestCase) {
    let lead = tc.draw(gs::integers::<u32>().min_value(1).max_value(19));
    let ins = tc.draw(gs::integers::<u32>().min_value(1).max_value(7));
    let trail = tc.draw(gs::integers::<u32>().min_value(1).max_value(19));
    let store = store_with_insertion(10, lead, ins, trail);
    let first_inserted = first_inserted_from_pairs(&store) as usize;
    let seq: Vec<Base> = store.record(ri(0)).unwrap().seq().to_vec();

    let mut engine = PileupEngine::new(
        store.prepare_for_pileup().input,
        Pos0::new(0).unwrap(),
        Pos0::new(10 + lead + trail + 10).unwrap(),
    );
    let mut checked = false;
    while let Some(col) = engine.pileups() {
        for view in col.alignments() {
            if !matches!(view.op(), PileupOp::Insertion { .. }) {
                continue;
            }
            let run = &seq[first_inserted..first_inserted + ins as usize];
            assert_eq!(view.inserted_bases(), run);
            checked = true;
        }
    }
    assert!(checked, "the fixture must produce an insertion column");
}
