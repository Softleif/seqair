//! Tests for pileup indel reporting: `PileupOp::Deletion`, `RefSkip`, Insertion.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    reason = "test code"
)]
#![allow(clippy::arithmetic_side_effects, reason = "test code")]
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    reason = "test code with known small values"
)]
mod helpers;

use helpers::{cigar_op, make_record, make_record_with_cigar};
use proptest::prelude::*;
use seqair::bam::{
    Pos0, RecordStore,
    pileup::{PileupEngine, PileupOp},
};

// ---- pileup_indel.deletions_included ----

// r[verify pileup_indel.deletions_included]
#[test]
fn deletion_positions_have_deletion_op() {
    // CIGAR: 10M 5D 10M at pos 100
    let raw = make_record_with_cigar(
        0,
        100,
        99,
        60,
        &[cigar_op(10, 0), cigar_op(5, 2), cigar_op(10, 0)],
        20,
    );
    let mut arena = RecordStore::new();
    arena.push_raw(&raw).unwrap();

    let mut engine = PileupEngine::new(arena, Pos0::new(100).unwrap(), Pos0::new(124).unwrap());
    let columns = helpers::collect_columns(&mut engine);

    // Positions 100-109: Match
    for col in columns
        .iter()
        .filter(|c| c.pos() >= Pos0::new(100).unwrap() && c.pos() <= Pos0::new(109).unwrap())
    {
        let aln = col.alignments().next().unwrap();
        assert!(
            matches!(aln.op, PileupOp::Match { .. }),
            "pos {} should be Match, got {:?}",
            col.pos(),
            aln.op
        );
    }

    // Positions 110-114: Deletion
    for col in columns
        .iter()
        .filter(|c| c.pos() >= Pos0::new(110).unwrap() && c.pos() <= Pos0::new(114).unwrap())
    {
        let aln = col.alignments().next().unwrap();
        assert!(aln.is_del(), "pos {} should be Deletion, got {:?}", col.pos(), aln.op);
        assert_eq!(aln.qpos(), None, "deletion should have no qpos");
        assert_eq!(aln.base(), None, "deletion should have no base");
        assert_eq!(aln.qual(), None, "deletion should have no qual");
    }

    // Positions 115-124: Match
    for col in columns
        .iter()
        .filter(|c| c.pos() >= Pos0::new(115).unwrap() && c.pos() <= Pos0::new(124).unwrap())
    {
        let aln = col.alignments().next().unwrap();
        assert!(
            matches!(aln.op, PileupOp::Match { .. }),
            "pos {} should be Match, got {:?}",
            col.pos(),
            aln.op
        );
    }

    // Total: 25 columns (10 + 5 + 10), not 20 (which was the old behavior)
    assert_eq!(columns.len(), 25, "should have 25 columns including deletion positions");
}

// ---- pileup_indel.refskips_included ----

// r[verify pileup_indel.refskips_included]
#[test]
fn refskip_positions_have_refskip_op() {
    // CIGAR: 10M 100N 10M at pos 0
    let raw = make_record_with_cigar(
        0,
        0,
        99,
        60,
        &[cigar_op(10, 0), cigar_op(100, 3), cigar_op(10, 0)],
        20,
    );
    let mut arena = RecordStore::new();
    arena.push_raw(&raw).unwrap();

    let mut engine = PileupEngine::new(arena, Pos0::new(0).unwrap(), Pos0::new(119).unwrap());
    let columns = helpers::collect_columns(&mut engine);

    // Should have 120 columns: 10 match + 100 refskip + 10 match
    assert_eq!(columns.len(), 120, "should include refskip positions");

    // Check a refskip position
    let col50 = columns.iter().find(|c| c.pos() == Pos0::new(50).unwrap()).unwrap();
    let aln = col50.alignments().next().unwrap();
    assert!(aln.is_refskip(), "pos 50 should be RefSkip");
    assert_eq!(aln.qpos(), None);
}

// ---- pileup_indel.depth_includes_all ----

// r[verify pileup_indel.depth_includes_all]
#[test]
fn depth_counts_deletions_and_refskips() {
    let mut arena = RecordStore::new();
    // Read 1: 20M (covers 0-19)
    arena.push_raw(&make_record(0, 0, 99, 60, 20)).unwrap();
    // Read 2: 5M 5D 10M at pos 0 (covers 0-19 with deletion at 5-9)
    arena
        .push_raw(&make_record_with_cigar(
            0,
            0,
            99,
            60,
            &[cigar_op(5, 0), cigar_op(5, 2), cigar_op(10, 0)],
            15,
        ))
        .unwrap();

    let mut engine = PileupEngine::new(arena, Pos0::new(0).unwrap(), Pos0::new(19).unwrap());
    let columns = helpers::collect_columns(&mut engine);

    // At deletion positions (5-9): depth should be 2 (one Match + one Deletion)
    let col7 = columns.iter().find(|c| c.pos() == Pos0::new(7).unwrap()).unwrap();
    assert_eq!(col7.depth(), 2, "depth should include deletion alignment");
}

// ---- pileup_indel.insertion_at_last_match ----

// r[verify pileup_indel.insertion_at_last_match]
#[test]
fn insertion_reported_at_last_match_before_insert() {
    // CIGAR: 10M 3I 10M at pos 0
    let raw = make_record_with_cigar(
        0,
        0,
        99,
        60,
        &[cigar_op(10, 0), cigar_op(3, 1), cigar_op(10, 0)],
        23,
    );
    let mut arena = RecordStore::new();
    arena.push_raw(&raw).unwrap();

    let mut engine = PileupEngine::new(arena, Pos0::new(0).unwrap(), Pos0::new(19).unwrap());
    let columns = helpers::collect_columns(&mut engine);
    assert_eq!(columns.len(), 20);

    // Position 9 (last base before insertion): should be Insertion with insert_len=3
    let col9 = columns.iter().find(|c| c.pos() == Pos0::new(9).unwrap()).unwrap();
    let aln = col9.alignments().next().unwrap();
    assert!(
        matches!(aln.op, PileupOp::Insertion { insert_len: 3, .. }),
        "pos 9 should be Insertion with len=3, got {:?}",
        aln.op
    );
    assert_eq!(aln.insert_len(), 3);
    assert_eq!(aln.qpos(), Some(9));

    // Position 8: regular Match (not the last before insertion)
    let col8 = columns.iter().find(|c| c.pos() == Pos0::new(8).unwrap()).unwrap();
    let aln8 = col8.alignments().next().unwrap();
    assert!(matches!(aln8.op, PileupOp::Match { .. }), "pos 8 should be Match, got {:?}", aln8.op);
    assert_eq!(aln8.insert_len(), 0);

    // Position 10 (first after insertion): Match with qpos = 10 + 3 = 13
    let col10 = columns.iter().find(|c| c.pos() == Pos0::new(10).unwrap()).unwrap();
    let aln10 = col10.alignments().next().unwrap();
    assert_eq!(aln10.qpos(), Some(13), "qpos after insertion should skip inserted bases");
}

// ---- pileup_indel.type_safety (compile-time) ----

// r[verify pileup_indel.type_safety]
#[test]
fn type_safety_deletion_has_no_base() {
    // This test verifies at compile time that you can't access base/qual/qpos
    // from a Deletion without matching. The fact that the convenience methods
    // return Option proves the API is safe.
    let op = PileupOp::Deletion { del_len: 1 };
    let _ = op; // Deletion carries del_len but no base/qual/qpos — that's the point
    // If someone tried: op.base — compile error. Must use match or convenience method.
}

// ---- pileup_indel.complex_indel ----

// r[verify pileup_indel.complex_indel]
#[test]
fn complex_indel_at_last_deletion_position() {
    // CIGAR: 10M 5D 3I 10M — insertion follows deletion → ComplexIndel at last D position
    let raw = make_record_with_cigar(
        0,
        0,
        99,
        60,
        &[cigar_op(10, 0), cigar_op(5, 2), cigar_op(3, 1), cigar_op(10, 0)],
        23,
    );
    let mut arena = RecordStore::new();
    arena.push_raw(&raw).unwrap();

    let mut engine = PileupEngine::new(arena, Pos0::new(0).unwrap(), Pos0::new(24).unwrap());
    let columns = helpers::collect_columns(&mut engine);

    // Position 9 (last match before D): should be plain Match, NOT Insertion
    // because the insertion follows the deletion, not this match
    let col9 = columns.iter().find(|c| c.pos() == Pos0::new(9).unwrap()).unwrap();
    let aln = col9.alignments().next().unwrap();
    assert!(
        matches!(aln.op, PileupOp::Match { .. }),
        "pos 9 should be Match (insertion is after deletion, not after this match), got {:?}",
        aln.op
    );

    // Interior deletion positions 10-13: plain Deletion
    for pos in 10..14u32 {
        let col = columns.iter().find(|c| c.pos() == Pos0::new(pos).unwrap()).unwrap();
        let aln = col.alignments().next().unwrap();
        assert!(
            matches!(aln.op, PileupOp::Deletion { del_len: 5 }),
            "pos {} should be plain Deletion(5), got {:?}",
            pos,
            aln.op
        );
    }

    // Last deletion position (14): ComplexIndel with del_len=5, insert_len=3
    let col14 = columns.iter().find(|c| c.pos() == Pos0::new(14).unwrap()).unwrap();
    let aln14 = col14.alignments().next().unwrap();
    assert!(
        matches!(aln14.op, PileupOp::ComplexIndel { del_len: 5, insert_len: 3, .. }),
        "pos 14 should be ComplexIndel(del=5,ins=3), got {:?}",
        aln14.op
    );
    assert!(aln14.is_del(), "ComplexIndel must report is_del() = true");
    assert_eq!(aln14.del_len(), 5);
    assert_eq!(aln14.insert_len(), 3);
    assert_eq!(aln14.qpos(), None, "ComplexIndel has no qpos");
    assert_eq!(aln14.base(), None, "ComplexIndel has no base");
    assert_eq!(aln14.qual(), None, "ComplexIndel has no qual");
}

// ---- insertion before deletion ----

// r[verify pileup_indel.insertion_at_last_match]
#[test]
fn insertion_before_deletion() {
    // CIGAR: 10M 3I 5D 10M at pos 0
    let raw = make_record_with_cigar(
        0,
        0,
        99,
        60,
        &[cigar_op(10, 0), cigar_op(3, 1), cigar_op(5, 2), cigar_op(10, 0)],
        23,
    );
    let mut arena = RecordStore::new();
    arena.push_raw(&raw).unwrap();

    let mut engine = PileupEngine::new(arena, Pos0::new(0).unwrap(), Pos0::new(24).unwrap());
    let columns = helpers::collect_columns(&mut engine);

    // pos 9: last M before I → should be Insertion with insert_len=3
    let col9 = columns.iter().find(|c| c.pos() == Pos0::new(9).unwrap()).unwrap();
    let aln = col9.alignments().next().unwrap();
    assert!(
        matches!(aln.op, PileupOp::Insertion { insert_len: 3, .. }),
        "pos 9 should be Insertion with len=3, got {:?}",
        aln.op
    );

    // pos 10-14: Deletion
    for pos in 10..15u32 {
        let col = columns.iter().find(|c| c.pos() == Pos0::new(pos).unwrap()).unwrap();
        let aln = col.alignments().next().unwrap();
        assert!(aln.is_del(), "pos {} should be Deletion, got {:?}", pos, aln.op);
    }

    // pos 15-24: Match (after deletion)
    let col15 = columns.iter().find(|c| c.pos() == Pos0::new(15).unwrap()).unwrap();
    let aln15 = col15.alignments().next().unwrap();
    assert!(
        matches!(aln15.op, PileupOp::Match { .. }),
        "pos 15 should be Match, got {:?}",
        aln15.op
    );
}

// r[verify pileup_indel.insertion_at_last_match]
// r[verify pileup_indel.complex_indel]
#[test]
fn insertion_with_anchor_and_complex_indel_after_deletion() {
    // CIGAR: 10M 2I 5D 3I 10M at pos 0
    // First I (len=2) has anchor at pos 9 → reported as Insertion
    // Second I (len=3) follows D → ComplexIndel at last D position (pos 14)
    let raw = make_record_with_cigar(
        0,
        0,
        99,
        60,
        &[cigar_op(10, 0), cigar_op(2, 1), cigar_op(5, 2), cigar_op(3, 1), cigar_op(10, 0)],
        25,
    );
    let mut arena = RecordStore::new();
    arena.push_raw(&raw).unwrap();

    let mut engine = PileupEngine::new(arena, Pos0::new(0).unwrap(), Pos0::new(24).unwrap());
    let columns = helpers::collect_columns(&mut engine);

    // pos 9: last M before first I → Insertion with insert_len=2
    let col9 = columns.iter().find(|c| c.pos() == Pos0::new(9).unwrap()).unwrap();
    let aln = col9.alignments().next().unwrap();
    assert!(
        matches!(aln.op, PileupOp::Insertion { insert_len: 2, .. }),
        "pos 9 should be Insertion with len=2 (anchored), got {:?}",
        aln.op
    );

    // pos 10-13: plain Deletion
    for pos in 10..14u32 {
        let col = columns.iter().find(|c| c.pos() == Pos0::new(pos).unwrap()).unwrap();
        let aln = col.alignments().next().unwrap();
        assert!(
            matches!(aln.op, PileupOp::Deletion { del_len: 5 }),
            "pos {} should be plain Deletion(5), got {:?}",
            pos,
            aln.op
        );
    }

    // pos 14: ComplexIndel (last D position, followed by I)
    let col14 = columns.iter().find(|c| c.pos() == Pos0::new(14).unwrap()).unwrap();
    let aln14 = col14.alignments().next().unwrap();
    assert!(
        matches!(aln14.op, PileupOp::ComplexIndel { del_len: 5, insert_len: 3, .. }),
        "pos 14 should be ComplexIndel(del=5,ins=3), got {:?}",
        aln14.op
    );

    // pos 15: Match — qpos accounts for both insertions: 10M + 2I + 5D + 3I = qpos 15
    let col15 = columns.iter().find(|c| c.pos() == Pos0::new(15).unwrap()).unwrap();
    let aln15 = col15.alignments().next().unwrap();
    assert!(
        matches!(aln15.op, PileupOp::Match { .. }),
        "pos 15 should be Match, got {:?}",
        aln15.op
    );
    // qpos at pos 15 should account for both insertions: 10 (first M) + 2 (first I) + 3 (second I) = 15
    assert_eq!(aln15.qpos(), Some(15), "qpos should skip both insertions' query bases");
}

// ---- minimal insertion ----

// r[verify pileup_indel.insertion_at_last_match]
#[test]
fn minimal_insertion_1m_1i_1m() {
    // CIGAR: 1M 1I 1M at pos 0
    let raw =
        make_record_with_cigar(0, 0, 99, 60, &[cigar_op(1, 0), cigar_op(1, 1), cigar_op(1, 0)], 3);
    let mut arena = RecordStore::new();
    arena.push_raw(&raw).unwrap();

    let mut engine = PileupEngine::new(arena, Pos0::new(0).unwrap(), Pos0::new(1).unwrap());
    let columns = helpers::collect_columns(&mut engine);
    assert_eq!(columns.len(), 2);

    // pos 0: single base that is BOTH first and last of its M block → Insertion
    let col0 = columns.iter().find(|c| c.pos() == Pos0::new(0).unwrap()).unwrap();
    let aln = col0.alignments().next().unwrap();
    assert!(
        matches!(aln.op, PileupOp::Insertion { insert_len: 1, .. }),
        "pos 0 should be Insertion with len=1, got {:?}",
        aln.op
    );
    assert_eq!(aln.qpos(), Some(0));

    // pos 1: Match after the insertion
    let col1 = columns.iter().find(|c| c.pos() == Pos0::new(1).unwrap()).unwrap();
    let aln1 = col1.alignments().next().unwrap();
    assert!(matches!(aln1.op, PileupOp::Match { .. }));
    assert_eq!(aln1.qpos(), Some(2)); // qpos 0 (M) + 1 (I) + 0 offset = 2
}

// ---- proptest: complex CIGARs with indels ----

use helpers::arb_read;

// r[verify pileup_indel.op_enum]
// r[verify pileup_indel.deletions_included]
proptest! {
    #![proptest_config(ProptestConfig::with_cases(200))]

    /// Every position in a read's ref span must produce a column.
    /// Deletions and ref-skips are now included, not skipped.
    #[test]
    fn every_ref_position_produces_column_for_single_read(read in arb_read()) {
        let mut arena = RecordStore::new();
        arena.push_raw(&read.raw).unwrap();

        let region_start = read.pos as u32;
        let region_end = region_start + read.ref_span - 1;
        let mut engine = PileupEngine::new(arena, Pos0::new(region_start).unwrap(), Pos0::new(region_end).unwrap());
        let columns = helpers::collect_columns(&mut engine);

        // With deletions/refskips included, every position in ref_span should have a column
        prop_assert_eq!(columns.len(), read.ref_span as usize,
            "should have one column per ref position (including D/N), got {} for ref_span {}",
            columns.len(), read.ref_span);
    }

    /// Match positions must have qpos, deletion/refskip positions must not.
    #[test]
    fn qpos_presence_matches_cigar_op_type(read in arb_read()) {
        let mut arena = RecordStore::new();
        arena.push_raw(&read.raw).unwrap();

        let region_start = read.pos as u32;
        let region_end = region_start + read.ref_span - 1;
        let mut engine = PileupEngine::new(arena, Pos0::new(region_start).unwrap(), Pos0::new(region_end).unwrap());

        let covered = read.covered_ref_positions();
        while let Some(col) = engine.pileups() {
            let aln = col.alignments().next().unwrap();
            if covered.contains(&col.pos().as_i64()) {
                prop_assert!(aln.qpos().is_some(),
                    "pos {} is covered (M/=/X) but has no qpos", col.pos());
            } else {
                prop_assert!(aln.qpos().is_none(),
                    "pos {} is not covered (D/N) but has qpos {:?}", col.pos(), aln.qpos());
                prop_assert!(aln.is_del() || aln.is_refskip(),
                    "pos {} is not covered but op is {:?}", col.pos(), aln.op);
            }
        }
    }
}

// r[verify pileup_indel.accessors]
// r[verify pileup_indel.op_enum]
proptest! {
    #![proptest_config(ProptestConfig::with_cases(200))]

    /// del_len() must equal the D CIGAR op length at deletion positions, and 0 elsewhere.
    #[test]
    fn del_len_matches_cigar_d_op_length(read in arb_read()) {
        let mut arena = RecordStore::new();
        arena.push_raw(&read.raw).unwrap();

        let region_start = read.pos as u32;
        let region_end = region_start + read.ref_span - 1;
        let mut engine = PileupEngine::new(arena, Pos0::new(region_start).unwrap(), Pos0::new(region_end).unwrap());

        while let Some(col) = engine.pileups() {
            let aln = col.alignments().next().unwrap();
            let rpos = col.pos().as_i64();
            match read.del_len_at(rpos) {
                Some(expected_len) => {
                    prop_assert!(aln.is_del(),
                        "pos {} should be Deletion but got {:?}", rpos, aln.op);
                    prop_assert_eq!(aln.del_len(), expected_len,
                        "pos {} del_len mismatch: expected {}, got {}", rpos, expected_len, aln.del_len());
                }
                None => {
                    prop_assert_eq!(aln.del_len(), 0,
                        "pos {} is not a deletion but del_len() returned {}", rpos, aln.del_len());
                }
            }
        }
    }

    /// Deletion ops must carry del_len > 0; non-deletion ops must have del_len == 0.
    #[test]
    fn del_len_nonzero_iff_deletion(read in arb_read()) {
        let mut arena = RecordStore::new();
        arena.push_raw(&read.raw).unwrap();

        let region_start = read.pos as u32;
        let region_end = region_start + read.ref_span - 1;
        let mut engine = PileupEngine::new(arena, Pos0::new(region_start).unwrap(), Pos0::new(region_end).unwrap());

        while let Some(col) = engine.pileups() {
            let aln = col.alignments().next().unwrap();
            if aln.is_del() {
                prop_assert!(aln.del_len() > 0,
                    "pos {}: Deletion op must have del_len > 0", col.pos());
            } else {
                prop_assert_eq!(aln.del_len(), 0,
                    "pos {}: non-Deletion op {:?} must have del_len == 0, got {}",
                    col.pos(), aln.op, aln.del_len());
            }
        }
    }

    /// All positions within the same D op must report the same del_len.
    #[test]
    fn del_len_consistent_across_deletion_span(read in arb_read()) {
        let mut arena = RecordStore::new();
        arena.push_raw(&read.raw).unwrap();

        let region_start = read.pos as u32;
        let region_end = region_start + read.ref_span - 1;
        let mut engine = PileupEngine::new(arena, Pos0::new(region_start).unwrap(), Pos0::new(region_end).unwrap());
        let columns = helpers::collect_columns(&mut engine);

        // Group consecutive deletion columns and verify they all have the same del_len
        let mut i = 0;
        while i < columns.len() {
            let aln = columns[i].alignments().next().unwrap();
            if aln.is_del() {
                let expected_del_len = aln.del_len();
                let start = i;
                while i < columns.len() {
                    let a = columns[i].alignments().next().unwrap();
                    if !a.is_del() { break; }
                    prop_assert_eq!(a.del_len(), expected_del_len,
                        "pos {}: del_len {} differs from first deletion pos {} del_len {}",
                        columns[i].pos(), a.del_len(),
                        columns[start].pos(), expected_del_len);
                    i += 1;
                }
                // The number of consecutive deletion columns must equal del_len
                let span = i - start;
                prop_assert_eq!(span as u32, expected_del_len,
                    "deletion span {} does not match del_len {}", span, expected_del_len);
            } else {
                i += 1;
            }
        }
    }
}
