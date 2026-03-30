//! Tests for CIGAR operations, CigarMapping, and pos_info_at computation.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::indexing_slicing)]
mod helpers;

use helpers::{cigar_bytes, cigar_op};
use proptest::prelude::*;
use seqair::bam::cigar::{CigarMapping, CigarPosInfo, calc_matches_indels};
use seqair::bam::{Pos, Zero};

// ---- cigar.matches_indels ----

// r[verify cigar.matches_indels]
#[test]
fn matches_indels_simple_match() {
    let ops = cigar_bytes(&[cigar_op(100, 0)]);
    let (matches, indels) = calc_matches_indels(&ops);
    assert_eq!(matches, 100);
    assert_eq!(indels, 0);
}

// r[verify cigar.matches_indels]
#[test]
fn matches_indels_with_insertions() {
    let ops = cigar_bytes(&[cigar_op(30, 0), cigar_op(5, 1), cigar_op(20, 0)]);
    let (matches, indels) = calc_matches_indels(&ops);
    assert_eq!(matches, 50);
    assert_eq!(indels, 5);
}

// r[verify cigar.matches_indels]
#[test]
fn matches_indels_with_deletions() {
    let ops = cigar_bytes(&[cigar_op(30, 0), cigar_op(5, 2), cigar_op(20, 0)]);
    let (matches, indels) = calc_matches_indels(&ops);
    assert_eq!(matches, 50);
    assert_eq!(indels, 5);
}

// r[verify cigar.matches_indels]
#[test]
fn matches_indels_soft_clip_ignored() {
    let ops = cigar_bytes(&[cigar_op(5, 4), cigar_op(90, 0), cigar_op(5, 4)]);
    let (matches, indels) = calc_matches_indels(&ops);
    assert_eq!(matches, 90);
    assert_eq!(indels, 0);
}

// ---- cigar.index + cigar.qpos_at ----

// r[verify cigar.index]
// r[verify cigar.qpos_at]
#[test]
fn cigar_mapping_simple_match() {
    let ops = cigar_bytes(&[cigar_op(100, 0)]);
    let mapping = CigarMapping::new(Pos::<Zero>::new(200).unwrap(), &ops);
    assert_eq!(
        mapping.pos_info_at(Pos::<Zero>::new(200).unwrap()),
        Some(CigarPosInfo::Match { qpos: 0 })
    );
    assert_eq!(
        mapping.pos_info_at(Pos::<Zero>::new(250).unwrap()),
        Some(CigarPosInfo::Match { qpos: 50 })
    );
    assert_eq!(
        mapping.pos_info_at(Pos::<Zero>::new(299).unwrap()),
        Some(CigarPosInfo::Match { qpos: 99 })
    );
    assert_eq!(mapping.pos_info_at(Pos::<Zero>::new(199).unwrap()), None);
    assert_eq!(mapping.pos_info_at(Pos::<Zero>::new(300).unwrap()), None);
}

// r[verify cigar.qpos_at]
#[test]
fn cigar_mapping_with_soft_clip() {
    let ops = cigar_bytes(&[cigar_op(5, 4), cigar_op(90, 0), cigar_op(5, 4)]);
    let mapping = CigarMapping::new(Pos::<Zero>::new(100).unwrap(), &ops);
    assert_eq!(
        mapping.pos_info_at(Pos::<Zero>::new(100).unwrap()),
        Some(CigarPosInfo::Match { qpos: 5 })
    );
    assert_eq!(
        mapping.pos_info_at(Pos::<Zero>::new(189).unwrap()),
        Some(CigarPosInfo::Match { qpos: 94 })
    );
    assert_eq!(mapping.pos_info_at(Pos::<Zero>::new(99).unwrap()), None);
    assert_eq!(mapping.pos_info_at(Pos::<Zero>::new(190).unwrap()), None);
}

// r[verify cigar.qpos_at]
#[test]
fn cigar_mapping_with_deletion() {
    let ops = cigar_bytes(&[cigar_op(30, 0), cigar_op(5, 2), cigar_op(20, 0)]);
    let mapping = CigarMapping::new(Pos::<Zero>::new(100).unwrap(), &ops);
    assert_eq!(
        mapping.pos_info_at(Pos::<Zero>::new(129).unwrap()),
        Some(CigarPosInfo::Match { qpos: 29 })
    );
    assert_eq!(mapping.pos_info_at(Pos::<Zero>::new(130).unwrap()), Some(CigarPosInfo::Deletion)); // inside deletion
    assert_eq!(mapping.pos_info_at(Pos::<Zero>::new(134).unwrap()), Some(CigarPosInfo::Deletion));
    assert_eq!(
        mapping.pos_info_at(Pos::<Zero>::new(135).unwrap()),
        Some(CigarPosInfo::Match { qpos: 30 })
    );
    // after deletion
}

// r[verify cigar.qpos_at]
#[test]
fn cigar_mapping_with_insertion() {
    let ops = cigar_bytes(&[cigar_op(30, 0), cigar_op(5, 1), cigar_op(20, 0)]);
    let mapping = CigarMapping::new(Pos::<Zero>::new(100).unwrap(), &ops);
    // pos 129 is the last base of the 30M block, and the next op is 5I — yields Insertion
    assert_eq!(
        mapping.pos_info_at(Pos::<Zero>::new(129).unwrap()),
        Some(CigarPosInfo::Insertion { qpos: 29, insert_len: 5 })
    );
    // pos 130 starts the 20M block; insertion skips query positions 30-34
    assert_eq!(
        mapping.pos_info_at(Pos::<Zero>::new(130).unwrap()),
        Some(CigarPosInfo::Match { qpos: 35 })
    );
}

// r[verify cigar.qpos_at]
#[test]
fn cigar_mapping_with_ref_skip() {
    let ops = cigar_bytes(&[cigar_op(30, 0), cigar_op(1000, 3), cigar_op(20, 0)]);
    let mapping = CigarMapping::new(Pos::<Zero>::new(100).unwrap(), &ops);
    assert_eq!(mapping.pos_info_at(Pos::<Zero>::new(130).unwrap()), Some(CigarPosInfo::RefSkip)); // inside N skip
    assert_eq!(
        mapping.pos_info_at(Pos::<Zero>::new(1130).unwrap()),
        Some(CigarPosInfo::Match { qpos: 30 })
    );
}

// r[verify cigar.qpos_at]
#[test]
fn cigar_mapping_seq_match_and_mismatch() {
    let ops = cigar_bytes(&[cigar_op(50, 7), cigar_op(10, 8), cigar_op(40, 7)]);
    let mapping = CigarMapping::new(Pos::<Zero>::new(0).unwrap(), &ops);
    assert_eq!(
        mapping.pos_info_at(Pos::<Zero>::new(50).unwrap()),
        Some(CigarPosInfo::Match { qpos: 50 })
    ); // X op
    assert_eq!(
        mapping.pos_info_at(Pos::<Zero>::new(60).unwrap()),
        Some(CigarPosInfo::Match { qpos: 60 })
    ); // back to =
    assert_eq!(mapping.pos_info_at(Pos::<Zero>::new(100).unwrap()), None);
}

// r[verify cigar.operations]
#[test]
fn cigar_mapping_hard_clip_ignored() {
    let ops = cigar_bytes(&[cigar_op(5, 5), cigar_op(90, 0), cigar_op(5, 5)]);
    let mapping = CigarMapping::new(Pos::<Zero>::new(100).unwrap(), &ops);
    assert_eq!(
        mapping.pos_info_at(Pos::<Zero>::new(100).unwrap()),
        Some(CigarPosInfo::Match { qpos: 0 })
    );
    assert_eq!(mapping.pos_info_at(Pos::<Zero>::new(190).unwrap()), None);
}

// r[verify io.named_constants]
#[test]
fn cigar_op_constants_match_sam_spec() {
    use seqair::bam::cigar::*;
    assert_eq!(CIGAR_M, 0);
    assert_eq!(CIGAR_I, 1);
    assert_eq!(CIGAR_D, 2);
    assert_eq!(CIGAR_N, 3);
    assert_eq!(CIGAR_S, 4);
    assert_eq!(CIGAR_H, 5);
    assert_eq!(CIGAR_P, 6);
    assert_eq!(CIGAR_EQ, 7);
    assert_eq!(CIGAR_X, 8);
}

// ---- CigarOpType enum ----

// r[verify io.typed_cigar_ops]
#[test]
fn cigar_op_type_roundtrip_all_valid_codes() {
    use seqair::bam::cigar::CigarOpType;

    for code in 0..=8u8 {
        let op = CigarOpType::from_bam(code);
        assert!(op.is_some(), "valid code {code} should parse");
    }
    // Invalid codes
    assert!(CigarOpType::from_bam(9).is_none());
    assert!(CigarOpType::from_bam(255).is_none());
}

// r[verify io.typed_cigar_ops]
#[test]
fn cigar_op_type_consumes_ref_and_query() {
    use seqair::bam::cigar::CigarOpType;

    // M, =, X consume both
    assert!(CigarOpType::Match.consumes_ref());
    assert!(CigarOpType::Match.consumes_query());
    assert!(CigarOpType::SeqMatch.consumes_ref());
    assert!(CigarOpType::SeqMatch.consumes_query());
    assert!(CigarOpType::SeqMismatch.consumes_ref());
    assert!(CigarOpType::SeqMismatch.consumes_query());

    // I, S consume query only
    assert!(!CigarOpType::Insertion.consumes_ref());
    assert!(CigarOpType::Insertion.consumes_query());
    assert!(!CigarOpType::SoftClip.consumes_ref());
    assert!(CigarOpType::SoftClip.consumes_query());

    // D, N consume ref only
    assert!(CigarOpType::Deletion.consumes_ref());
    assert!(!CigarOpType::Deletion.consumes_query());
    assert!(CigarOpType::RefSkip.consumes_ref());
    assert!(!CigarOpType::RefSkip.consumes_query());

    // H, P consume neither
    assert!(!CigarOpType::HardClip.consumes_ref());
    assert!(!CigarOpType::HardClip.consumes_query());
    assert!(!CigarOpType::Padding.consumes_ref());
    assert!(!CigarOpType::Padding.consumes_query());
}

// ---- property tests ----

const MATCH_OPS: [u8; 3] = [0, 7, 8];
const REF_ONLY_OPS: [u8; 2] = [2, 3];

fn arb_cigar() -> impl Strategy<Value = Vec<(u32, u8)>> {
    prop::option::of((1u32..=20, Just(4u8)))
        .prop_flat_map(|leading_clip| {
            let inner_ops = prop::collection::vec(
                prop_oneof![
                    (1u32..=200, prop::sample::select(&MATCH_OPS[..])),
                    (1u32..=20, prop::sample::select(&REF_ONLY_OPS[..])),
                    (1u32..=20, Just(1u8)), // I
                ],
                1..=8,
            );
            let trailing_clip = prop::option::of((1u32..=20, Just(4u8)));
            (Just(leading_clip), inner_ops, trailing_clip)
        })
        .prop_map(|(lead, inner, trail)| {
            let mut ops = Vec::new();
            if let Some(clip) = lead {
                ops.push(clip);
            }
            ops.extend(inner);
            if let Some(clip) = trail {
                ops.push(clip);
            }
            ops
        })
        .prop_filter("need at least one ref-consuming op", |ops| {
            ops.iter().any(|(_, op)| matches!(op, 0 | 2 | 3 | 7 | 8))
        })
}

// r[verify cigar.qpos_at]
proptest! {
    #[test]
    fn qpos_monotonically_increasing(ops in arb_cigar(), start in 0u32..1_000_000) {
        let packed: Vec<u32> = ops.iter().map(|&(len, op)| cigar_op(len, op)).collect();
        let bytes = cigar_bytes(&packed);
        let mapping = CigarMapping::new(Pos::<Zero>::new(start).unwrap(), &bytes);

        let ref_span: u32 = ops.iter().map(|&(len, op)| match op {
            0 | 2 | 3 | 7 | 8 => len,
            _ => 0,
        }).sum();

        let mut last_qpos: Option<u32> = None;
        for ref_pos in start..start + ref_span {
            let qpos = match mapping.pos_info_at(Pos::<Zero>::new(ref_pos).unwrap()) {
                Some(CigarPosInfo::Match { qpos }) => Some(qpos),
                Some(CigarPosInfo::Insertion { qpos, .. }) => Some(qpos),
                _ => None,
            };
            if let Some(q) = qpos {
                if let Some(prev) = last_qpos {
                    prop_assert!(q > prev, "not monotonic at ref {ref_pos}");
                }
                last_qpos = Some(q);
            }
        }
    }
}

// r[verify cigar.qpos_accuracy]
proptest! {
    #[test]
    fn pure_match_qpos_is_offset(len in 1u32..=500, start in 0u32..1_000_000) {
        let bytes = cigar_bytes(&[cigar_op(len, 0)]);
        let mapping = CigarMapping::new(Pos::<Zero>::new(start).unwrap(), &bytes);
        for offset in 0..len {
            prop_assert_eq!(
                mapping.pos_info_at(Pos::<Zero>::new(start + offset).unwrap()),
                Some(CigarPosInfo::Match { qpos: offset })
            );
        }
        prop_assert_eq!(mapping.pos_info_at(Pos::<Zero>::new(start + len).unwrap()), None);
    }
}

// r[verify cigar.qpos_at]
//
// Builds CIGARs from human-readable strings and derives expected Some/None
// from the string form — not from numeric op-code constants — so the test
// cannot be tautological.
proptest! {
    #[test]
    fn qpos_none_for_deletions_some_for_matches_from_strings(
        // Each part is (length, op_char) drawn from the SAM alphabet.
        // At least one ref-consuming op must be present so the CIGAR is valid.
        parts in prop::collection::vec(
            prop_oneof![
                (1u32..=200u32, Just('M')),
                (1u32..=200u32, Just('=')),
                (1u32..=200u32, Just('X')),
                (1u32..=20u32,  Just('D')),
                (1u32..=20u32,  Just('N')),
                (1u32..=20u32,  Just('I')),
            ],
            1..=8,
        ).prop_filter("need at least one ref-consuming op", |parts| {
            parts.iter().any(|(_, op)| matches!(op, 'M' | '=' | 'X' | 'D' | 'N'))
        }),
        start in 0u32..1_000_000,
    ) {
        // Build the byte representation from the human-readable string.
        let cigar_str: String = parts.iter().map(|(len, op)| format!("{len}{op}")).collect();
        let op_char_to_code = |c: char| -> u8 {
            match c {
                'M' => 0, 'I' => 1, 'D' => 2, 'N' => 3, 'S' => 4,
                'H' => 5, 'P' => 6, '=' => 7, 'X' => 8, _ => 0,
            }
        };
        let packed: Vec<u32> = parts.iter()
            .map(|&(len, op)| cigar_op(len, op_char_to_code(op)))
            .collect();
        let bytes = cigar_bytes(&packed);
        let mapping = CigarMapping::new(Pos::<Zero>::new(start).unwrap(), &bytes);

        // Walk the string-form parts to derive expectations independently.
        let mut ref_pos = start;
        for &(len, op) in &parts {
            match op {
                // M / = / X consume ref+query → every covered ref position must
                // return Some(Match) or Some(Insertion).
                'M' | '=' | 'X' => {
                    for i in 0..len {
                        let pos = ref_pos + i;
                        let result = mapping.pos_info_at(Pos::<Zero>::new(pos).unwrap());
                        prop_assert!(
                            matches!(
                                result,
                                Some(CigarPosInfo::Match { .. }) | Some(CigarPosInfo::Insertion { .. })
                            ),
                            "cigar={}: ref {} under {} op should be Match or Insertion, got {:?}",
                            cigar_str, pos, op, result
                        );
                    }
                    ref_pos += len;
                }
                // D consumes ref only → Deletion at every covered position.
                'D' => {
                    for i in 0..len {
                        let pos = ref_pos + i;
                        prop_assert_eq!(
                            mapping.pos_info_at(Pos::<Zero>::new(pos).unwrap()),
                            Some(CigarPosInfo::Deletion),
                            "cigar={}: ref {} under D op should be Deletion",
                            cigar_str, pos
                        );
                    }
                    ref_pos += len;
                }
                // N consumes ref only → RefSkip at every covered position.
                'N' => {
                    for i in 0..len {
                        let pos = ref_pos + i;
                        prop_assert_eq!(
                            mapping.pos_info_at(Pos::<Zero>::new(pos).unwrap()),
                            Some(CigarPosInfo::RefSkip),
                            "cigar={}: ref {} under N op should be RefSkip",
                            cigar_str, pos
                        );
                    }
                    ref_pos += len;
                }
                // I consumes query only — no ref positions to check.
                _ => {}
            }
        }
    }
}
