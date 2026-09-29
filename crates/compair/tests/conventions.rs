//! The GATK conventions the 104 published vectors cannot pin, pinned by hand.
//!
//! The vectors carry a **uniform** gap-continuation column (`+`, Q10) and
//! `read-ins-qual == read-del-qual` on every row, and their indel qualities are
//! all Q36-Q45. Three consequences, each measured by mutating the crate and
//! re-running `reference_reproduces_gatk`:
//!
//! - replacing `1 - (p_ins + p_del)` with `(1 - p_ins)(1 - p_del)` still
//!   passes all 104 at 1e-4 (the two differ by `p_ins * p_del`, ~6e-8 per row),
//! - swapping the insertion and deletion roles still passes all 104,
//! - **adding the deletion matrix to the final sum still passes all 104** — the
//!   excluded term is ~1e-4 of the total at Q40, just under the gate.
//!
//! So this file computes one 3 x 2 dynamic program by hand and asserts the
//! whole recurrence at once: the `1 / haplotype_len` start spread over every
//! haplotype column, the exclusion of the deletion matrix from the total, the
//! summed `match_to_match` with gap-open qualities raised to Q6 as GATK raises
//! them, the insertion/deletion role split, and
//! `indel_to_match = 1 - gap_continuation` read at the right read base.

use compair::{
    Band, BaseQuality, Haplotype, Read, StandardEmission, Strand, align_banded, align_banded_simd,
    align_full, align_strips, align_strips_simd,
};

/// `align_full` on `haplotype = ACT`, `read = AC`, written out cell by cell.
///
/// Every quality is distinct so that a transition taken from the wrong read
/// base, or an insertion quality used where a deletion quality belongs, moves
/// the answer.
fn expected_act_vs_ac(q: [u8; 2], ins: [u8; 2], del: [u8; 2], gcp: [u8; 2]) -> f64 {
    let p = |phred: u8| 10f64.powf(-f64::from(phred) / 10.0);
    let gap_open = |phred: u8| p(phred.max(6));
    let (e1, e2) = (p(q[0]), p(q[1]));
    let (pd1, g1) = (gap_open(del[0]), p(gcp[0]));
    let (pi2, pd2, g2) = (gap_open(ins[1]), gap_open(del[1]), p(gcp[1]));
    let mm2 = 1.0 - (pi2 + pd2);
    let (im1, im2) = (1.0 - g1, 1.0 - g2);
    let init = 1.0 / 3.0;

    // Row 1 emits read base `A` against `A`, `C`, `T`; its only predecessor is
    // the start row, which is `init` in the deletion matrix at every column.
    let m11 = (1.0 - e1) * init * im1;
    let m12 = (e1 / 3.0) * init * im1;
    let m13 = (e1 / 3.0) * init * im1;
    let d12 = m11 * pd1;

    // Row 2 emits `C`. `m21` is zero because column 0 is zero below row 0.
    let m22 = (1.0 - e2) * (m11 * mm2);
    let m23 = (e2 / 3.0) * (m12 * mm2 + d12 * im2);
    let (i21, i22, i23) = (m11 * pi2, m12 * pi2, m13 * pi2);

    // The deletion matrix of the last row is *not* summed: `d23 = m22 * pd2`
    // is a real, large value here and must not appear.
    i21 + m22 + i22 + m23 + i23
}

fn act_vs_ac(q: [u8; 2], ins: [u8; 2], del: [u8; 2], gcp: [u8; 2]) -> (Haplotype, Read) {
    let haplotype = Haplotype::from_ascii(b"ACT");
    let bases = Haplotype::from_ascii(b"AC").bases().to_vec();
    let read = Read::new(
        bases,
        &q.map(BaseQuality::from_byte),
        &ins.map(BaseQuality::from_byte),
        &del.map(BaseQuality::from_byte),
        &gcp.map(BaseQuality::from_byte),
        Strand::OT,
    )
    .expect("two bases and four two-entry quality tracks");
    (haplotype, read)
}

/// Base qualities, insertion qualities, deletion qualities and
/// gap-continuation qualities, two read bases each.
type QualitySet = ([u8; 2], [u8; 2], [u8; 2], [u8; 2]);

/// The parameter sets, and what each one is for.
const CASES: [QualitySet; 4] = [
    // p_ins + p_del = 0.26: this separates `1 - (p_ins + p_del)` from
    // `(1 - p_ins)(1 - p_del)` (0.7388 vs 0.7413), and p_del = 0.25 makes the
    // excluded `d23` a large share of the total.
    ([20, 25], [10, 20], [15, 6], [8, 15]),
    // Both gap opens at Q3, which count as Q6: `match_to_match` is 0.498, where
    // taken at face value the two would sum past one and leave it at zero.
    ([30, 12], [10, 3], [15, 3], [8, 15]),
    // An insertion quality of Q0, a literal `p_ins = 1`, counts as Q6 too.
    ([30, 12], [10, 0], [15, 20], [8, 15]),
    // GATK's own range, where every convention costs only a few 1e-8 per row --
    // which is why the vectors miss them.
    ([37, 22], [36, 45], [40, 39], [10, 10]),
];

#[test]
fn the_reference_recurrence_cell_by_cell() {
    for (index, (q, ins, del, gcp)) in CASES.into_iter().enumerate() {
        let (haplotype, read) = act_vs_ac(q, ins, del, gcp);
        let want = expected_act_vs_ac(q, ins, del, gcp).log10();
        let got = align_full(&haplotype, &read, &StandardEmission::default()).get();
        assert!((got - want).abs() < 1e-12, "case {index}: align_full {got}, hand-computed {want}");
    }
}

/// The banded kernels must use the same conventions, not merely agree with each
/// other. A band wider than the matrix visits every cell, so any difference
/// here is a difference in the recurrence rather than in the band.
#[test]
fn the_banded_kernels_use_the_same_recurrence() {
    let band = Band::new(64, 0).expect("wide enough for a 3 x 2 matrix");
    for (index, (q, ins, del, gcp)) in CASES.into_iter().enumerate() {
        let (haplotype, read) = act_vs_ac(q, ins, del, gcp);
        let want = expected_act_vs_ac(q, ins, del, gcp).log10();
        let scalar = align_banded(&haplotype, &read, &StandardEmission::default(), band).get();
        let simd = align_banded_simd(&haplotype, &read, &StandardEmission::default(), band).get();
        assert_eq!(scalar.to_bits(), simd.to_bits(), "case {index}: scalar {scalar} simd {simd}");
        assert!(
            (scalar - want).abs() < 1e-5,
            "case {index}: align_banded {scalar}, hand-computed {want}"
        );
        let scalar = align_strips(&haplotype, &read, &StandardEmission::default(), band).get();
        let simd = align_strips_simd(&haplotype, &read, &StandardEmission::default(), band).get();
        assert_eq!(scalar.to_bits(), simd.to_bits(), "case {index}: strips {scalar} simd {simd}");
        assert!(
            (scalar - want).abs() < 1e-5,
            "case {index}: align_strips {scalar}, hand-computed {want}"
        );
    }
}
