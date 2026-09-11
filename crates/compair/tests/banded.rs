#[allow(dead_code, reason = "each integration test uses a different part of this")]
mod support;

use compair::{
    Band, Base, BaseQuality, Betas, ConversionModel, Haplotype, Probability, Read,
    StandardEmission, Strand, TapsEmission, Workspace, align_banded, align_banded_simd, align_full,
    error_probability,
};
use proptest::prelude::*;
use support::{
    any_conversion, any_probability, arbitrary_case, derived_case, plausible_conversion,
};

proptest! {
    // The bit-parity gate is the cheapest of these and the one most worth
    // running hard, so it alone does not take the default 256. The other
    // properties live in their own blocks below: `proptest_config` is
    // block-wide, and a full-matrix oracle at 2048 cases was most of the
    // suite's wall-clock.
    #![proptest_config(ProptestConfig { cases: 2048, ..ProptestConfig::default() })]

    /// The C4 gate: the scalar and the eight-wide kernels agree to the bit, on
    /// inputs that include bands missing the alignment entirely, because parity
    /// has to hold there too. And a reused `Workspace` is the same kernel:
    /// nothing from the previous pair leaks into the next.
    #[test]
    fn simd_is_bit_identical_to_scalar(
        case in arbitrary_case(),
        conversion in any_conversion(),
        uniform in any_probability(),
    ) {
        let band = case.band();
        let taps = TapsEmission::new(conversion, Betas::PerSite(&case.betas));
        let uniform = TapsEmission::new(conversion, Betas::Uniform(uniform));
        let mut workspace = Workspace::new();
        for (name, banded, simd, reused) in [
            (
                "standard",
                align_banded(&case.haplotype, &case.read, &StandardEmission::default(), band),
                align_banded_simd(&case.haplotype, &case.read, &StandardEmission::default(), band),
                workspace.align_banded_simd(&case.haplotype, &case.read, &StandardEmission::default(), band),
            ),
            (
                "taps",
                align_banded(&case.haplotype, &case.read, &taps, band),
                align_banded_simd(&case.haplotype, &case.read, &taps, band),
                workspace.align_banded(&case.haplotype, &case.read, &taps, band),
            ),
            (
                "uniform",
                align_banded(&case.haplotype, &case.read, &uniform, band),
                align_banded_simd(&case.haplotype, &case.read, &uniform, band),
                workspace.align_banded_simd(&case.haplotype, &case.read, &uniform, band),
            ),
        ] {
            prop_assert_eq!(
                banded.get().to_bits(), simd.get().to_bits(),
                "{}: scalar {:?} vs simd {:?}", name, banded, simd
            );
            prop_assert_eq!(
                banded.get().to_bits(), reused.get().to_bits(),
                "{}: fresh {:?} vs reused workspace {:?}", name, banded, reused
            );
        }
    }
}

proptest! {
    /// The other C4 gate: where the optimal path is inside the band, the f32
    /// band agrees with the f64 full matrix. The reads here are cut out of the
    /// haplotype and given at most three edits, so a 48-column band contains the
    /// path by construction.
    ///
    /// `plausible_conversion` and not `any_conversion`, and the distinction is
    /// the premise rather than a tolerance: under `c = 0, f = 1` every cytosine
    /// converts, the read is unlikely against its own haplotype, and alignments
    /// elsewhere beat the seeded one -- so the band legitimately misses the
    /// optimum and this gate would be testing the band, not the kernel. The
    /// kernel against the same band, under *any* conversion model, is
    /// `the_f32_band_is_the_f64_recurrence_over_the_same_band`.
    #[test]
    fn banded_matches_the_reference_inside_the_band(
        case in derived_case(3),
        conversion in plausible_conversion(),
    ) {
        let band = case.band();
        let taps = TapsEmission::new(conversion, Betas::PerSite(&case.betas));
        for (name, banded, full) in [
            (
                "standard",
                align_banded(&case.haplotype, &case.read, &StandardEmission::default(), band).get(),
                align_full(&case.haplotype, &case.read, &StandardEmission::default()).get(),
            ),
            (
                "taps",
                align_banded(&case.haplotype, &case.read, &taps, band).get(),
                align_full(&case.haplotype, &case.read, &taps).get(),
            ),
        ] {
            prop_assert!(
                (banded - full).abs() < 1e-3,
                "{name}: banded {banded} against reference {full}"
            );
        }
    }

    /// Where the chemistry cannot act -- a haplotype of `A` and `T` only, so
    /// no cytosine on either strand -- `TapsEmission` is `StandardEmission`
    /// through the whole DP, not only per cell.
    ///
    /// This used to admit `C` as well, on the reading that the conversion rows
    /// were `CpG`-only. They are not: design §4.1 puts every unmethylated
    /// cytosine on the false-conversion row, so a non-`CpG` `C` on an OT read
    /// differs from a plain mismatch by `f`. See
    /// `a_non_cpg_cytosine_reads_t_at_the_false_conversion_rate`.
    #[test]
    fn taps_equals_standard_where_the_chemistry_cannot_act(
        bases in proptest::collection::vec(
            prop_oneof![Just(Base::A), Just(Base::T)],
            30..90,
        ),
        quals in proptest::collection::vec(2u8..=45, 30..90),
        strand in support::any_strand(),
        conversion in any_conversion(),
    ) {
        let haplotype = Haplotype::new(bases.clone());
        let read_bases: Vec<Base> = bases.iter().copied().take(quals.len().min(bases.len())).collect();
        let quals: Vec<BaseQuality> = quals.iter().take(read_bases.len()).map(|q| BaseQuality::from_byte(*q)).collect();
        let read = Read::uniform(read_bases, &quals, BaseQuality::from_byte(45), BaseQuality::from_byte(45), BaseQuality::from_byte(10), strand)
            .map_err(|_| TestCaseError::reject("read"))?;
        let betas = vec![Probability::ONE; haplotype.len()];
        let taps = TapsEmission::new(conversion, Betas::PerSite(&betas));
        prop_assert_eq!(
            align_full(&haplotype, &read, &taps).get().to_bits(),
            align_full(&haplotype, &read, &StandardEmission::default()).get().to_bits()
        );
    }
}

/// What happens when the optimal alignment is not inside the band: the score is
/// the best path that stays inside it, so it is finite and strictly lower.
/// It is never an error and never `-inf`, and a caller comparing two haplotypes
/// under the same band is still comparing like with like.
#[test]
fn a_band_that_misses_the_alignment_underestimates() {
    // Non-repetitive on purpose: the full matrix marginalises over every start
    // position, so on a tandem repeat it legitimately beats any single band.
    let haplotype = Haplotype::from_ascii(
        b"TCATTGGCTATCCTAACCCGACCCTAGGAGCGGTTGGCGTGTATGCCGTGAATTTTCTCATTTCCGCTAGACATAATCGTTCTGCCTATA",
    );
    let read = Read::uniform(
        haplotype.bases().get(45..80).expect("in range").to_vec(),
        &[BaseQuality::from_byte(35); 35],
        BaseQuality::from_byte(45),
        BaseQuality::from_byte(45),
        BaseQuality::from_byte(10),
        Strand::OT,
    )
    .expect("valid");

    let full = align_full(&haplotype, &read, &StandardEmission::default()).get();
    let on_target =
        align_banded(&haplotype, &read, &StandardEmission::default(), Band::anchored(45)).get();
    let missed =
        align_banded(&haplotype, &read, &StandardEmission::default(), Band::anchored(0)).get();

    assert!(
        (on_target - full).abs() < 1e-3,
        "an anchored band reproduces the reference: {on_target} against {full}"
    );
    assert!(
        missed.is_finite(),
        "a band that misses the alignment still returns a score, got {missed}"
    );
    assert!(
        missed < on_target - 1.0,
        "and it is strictly and visibly lower: {missed} against {on_target}"
    );
}

/// A property the model does **not** have, pinned so that nobody adds it as a
/// test and watches it fail.
///
/// The transitions at row `i` are derived from read base `i - 1`, so a gap is
/// charged the insertion or deletion quality of the base on one particular side
/// of it. Reverse-complementing puts the other base there, and the two scores
/// part company. With *uniform* indel qualities the effect all but vanishes
/// (~1e-4, from the boundary alone: the free start is encoded in the deletion
/// row and cannot begin in an insertion, while the free end is a plain sum over
/// the match and insertion rows); it is the per-base indel qualities that make
/// it large. So the OT/OB mirror is a property of the *emission*, tested exactly
/// in `tests/emission.rs`, and the DP's only dependence on strand is that one
/// call.
#[test]
#[allow(
    clippy::cast_possible_truncation,
    reason = "a deliberately tiny deterministic generator over 90-base sequences"
)]
fn the_dp_is_not_reverse_complement_symmetric() {
    let mut state = 0x2545_f491_4f6c_dd1d_u64;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let mut worst = 0.0f64;
    for _ in 0..200 {
        let length = 60 + (next() % 60) as usize;
        let bases: Vec<Base> = (0..length)
            .map(|_| match next() % 4 {
                0 => Base::A,
                1 => Base::C,
                2 => Base::G,
                _ => Base::T,
            })
            .collect();
        let start = (next() % 20) as usize;
        let read_length = 30 + (next() % 25) as usize;
        let end = (start + read_length).min(length);
        let Some(slice) = bases.get(start..end) else {
            continue;
        };
        let mut read_bases = slice.to_vec();
        // Mismatches matter: with a perfect substring one alignment dominates and
        // the asymmetry stays near the boundary term. Competing gapped paths are
        // what expose the per-base indel qualities.
        for _ in 0..(next() % 6) {
            let at = (next() as usize) % read_bases.len();
            if let Some(slot) = read_bases.get_mut(at) {
                *slot = match next() % 4 {
                    0 => Base::A,
                    1 => Base::C,
                    2 => Base::G,
                    _ => Base::T,
                };
            }
        }
        let quals: Vec<BaseQuality> =
            read_bases.iter().map(|_| BaseQuality::from_byte(6 + (next() % 34) as u8)).collect();
        let haplotype = Haplotype::new(bases);
        let insertion_quals: Vec<BaseQuality> =
            read_bases.iter().map(|_| BaseQuality::from_byte(30 + (next() % 20) as u8)).collect();
        let deletion_quals: Vec<BaseQuality> =
            read_bases.iter().map(|_| BaseQuality::from_byte(30 + (next() % 20) as u8)).collect();
        let gap_quals = vec![BaseQuality::from_byte(10); read_bases.len()];
        let Ok(read) = Read::new(
            read_bases.clone(),
            &quals,
            &insertion_quals,
            &deletion_quals,
            &gap_quals,
            Strand::OT,
        ) else {
            continue;
        };
        let mirrored_bases: Vec<Base> = read_bases.iter().rev().map(|b| b.inverse()).collect();
        let mirrored_quals: Vec<BaseQuality> = quals.iter().rev().copied().collect();
        let Ok(mirrored_read) = Read::new(
            mirrored_bases,
            &mirrored_quals,
            &reversed(read.insertion_quals()),
            &reversed(read.deletion_quals()),
            &reversed(read.gap_quals()),
            Strand::OB,
        ) else {
            continue;
        };
        let forward = align_full(&haplotype, &read, &StandardEmission::default()).get();
        let backward = align_full(
            &haplotype.reverse_complement(),
            &mirrored_read,
            &StandardEmission::default(),
        )
        .get();
        worst = worst.max((forward - backward).abs());
    }
    assert!(worst > 0.5, "the asymmetry is real and worth documenting; worst seen {worst}");
}

fn reversed(quals: &[BaseQuality]) -> Vec<BaseQuality> {
    quals.iter().rev().copied().collect()
}

/// The conversion model this project measured from its own spike-in controls.
#[test]
fn the_bundled_conversion_model_is_the_measured_one() {
    let model = ConversionModel::taps_default();
    assert!((*model.efficiency() - 0.96).abs() < 1e-12);
    assert!((*model.false_conversion() - 0.004).abs() < 1e-12);
    assert!((model.conversion_rate(Probability::ZERO) - 0.004).abs() < 1e-12);
    assert!((model.conversion_rate(Probability::ONE) - 0.96).abs() < 1e-12);
}

/// An independent `f64` full dynamic program restricted to the band's
/// documented predicate, sharing nothing with the kernel under test but the
/// emission. Written out rather than derived from `align_full` so that a
/// change to the kernel's index algebra has to disagree with something.
#[allow(
    clippy::indexing_slicing,
    reason = "every index is in 0..=h and every row is allocated with h + 1 entries"
)]
#[allow(clippy::cast_possible_wrap, reason = "test sequences are a few hundred bases")]
fn align_masked<E: compair::Emission>(
    haplotype: &Haplotype,
    read: &Read,
    emission: &E,
    band: Band,
) -> f64 {
    let (h, r) = (haplotype.len(), read.len());
    if h == 0 || r == 0 {
        return f64::NEG_INFINITY;
    }
    let half = i64::from(band.width() / 2);
    let offset = i64::from(band.offset());
    let live = |i: usize, j: usize| ((j as i64) - (i as i64) - offset).abs() <= half;
    let quality = |track: &[BaseQuality], index: usize| {
        error_probability(track.get(index).copied().unwrap_or(BaseQuality::from_byte(0)))
    };

    #[allow(clippy::cast_precision_loss, reason = "test haplotypes are short")]
    let init = 1.0 / h as f64;
    let mut prev_m = vec![0.0f64; h + 1];
    let mut prev_i = vec![0.0f64; h + 1];
    let mut prev_d: Vec<f64> = (0..=h).map(|j| if live(0, j) { init } else { 0.0 }).collect();
    let (mut cur_m, mut cur_i, mut cur_d) =
        (vec![0.0f64; h + 1], vec![0.0f64; h + 1], vec![0.0f64; h + 1]);

    for i in 1..=r {
        let Some(observation) = read.observation(i - 1) else { return f64::NEG_INFINITY };
        let p_ins = quality(read.insertion_quals(), i - 1);
        let p_del = quality(read.deletion_quals(), i - 1);
        let gap = quality(read.gap_quals(), i - 1);
        let m2m = 1.0 - (p_ins + p_del).min(1.0);
        let i2m = 1.0 - gap;
        cur_m[0] = 0.0;
        cur_i[0] = 0.0;
        cur_d[0] = 0.0;
        for j in 1..=h {
            if !live(i, j) {
                cur_m[j] = 0.0;
                cur_i[j] = 0.0;
                cur_d[j] = 0.0;
                continue;
            }
            let Some(site) = haplotype.site(j - 1) else { return f64::NEG_INFINITY };
            let prior = emission.match_probability(site, observation);
            cur_m[j] = prior * (prev_m[j - 1] * m2m + prev_i[j - 1] * i2m + prev_d[j - 1] * i2m);
            cur_i[j] = prev_m[j] * p_ins + prev_i[j] * gap;
            cur_d[j] = cur_m[j - 1] * p_del + cur_d[j - 1] * gap;
        }
        core::mem::swap(&mut prev_m, &mut cur_m);
        core::mem::swap(&mut prev_i, &mut cur_i);
        core::mem::swap(&mut prev_d, &mut cur_d);
    }
    let total: f64 = prev_m.iter().zip(prev_i.iter()).skip(1).map(|(m, i)| m + i).sum();
    if total <= 0.0 { f64::NEG_INFINITY } else { total.log10() }
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 1024, ..ProptestConfig::default() })]

    /// The band is exactly the strip `Band` documents: `|j - i - offset| <=
    /// width / 2` and nothing else. `arbitrary_case` reaches offsets outside
    /// the matrix on both sides, which is where the kernel's clamps live.
    #[test]
    fn the_band_is_the_documented_predicate(
        case in arbitrary_case(),
        width in 2u32..64,
        conversion in any_conversion(),
    ) {
        let band = Band::new(width, case.offset).map_err(|_| TestCaseError::reject("width"))?;
        let taps = TapsEmission::new(conversion, Betas::PerSite(&case.betas));
        for (name, got, want) in [
            (
                "standard",
                align_banded(&case.haplotype, &case.read, &StandardEmission::default(), band).get(),
                align_masked(&case.haplotype, &case.read, &StandardEmission::default(), band),
            ),
            (
                "taps",
                align_banded(&case.haplotype, &case.read, &taps, band).get(),
                align_masked(&case.haplotype, &case.read, &taps, band),
            ),
        ] {
            prop_assert_eq!(
                got.is_finite(), want.is_finite(),
                "{}: kernel {} masked oracle {}", name, got, want
            );
            // The finite/infinite agreement above is the geometry claim, and
            // it is unconditional -- it is what caught B1. The *value* is only
            // compared where f32 can hold it: `width` here runs to 64 against
            // reads as short as one base, so past the read length a diagonal
            // crosses every read row and the kernel underflows by a fraction of
            // the score, up to several log10
            // (`an_unbanded_f32_run_underflows_where_a_real_band_does_not`).
            if got.is_finite() && want > -30.0 {
                prop_assert!(
                    (got - want).abs() < 1e-3 * (1.0 + want.abs()),
                    "{}: kernel {} masked oracle {}", name, got, want
                );
            }
        }
    }

    /// Widening a band can only add paths, so it can only raise the score, and
    /// no band can beat the unbanded reference.
    ///
    /// The second half is unconditional: the `f32` kernel only ever loses
    /// mass, to rounding and to the subnormals it flushes. The first half
    /// holds where `f32` holds the score. A flushed cell is below `2^-126`
    /// of its diagonal's maximum, so below ~1e-38 in absolute terms, and
    /// cannot move a total of 1e-30 or more; below that a wider band's
    /// diagonals span more dynamic range and can flush mass a narrower band
    /// kept -- the fuzzer found a pair at `-133` where eight more columns
    /// cost 0.16 -- which is the regime `Band::MAX_WIDTH` documents.
    #[test]
    fn a_band_is_monotone_in_its_width_and_never_beats_the_reference(
        case in arbitrary_case(),
        width in 2u32..40,
    ) {
        let narrow = Band::new(width, case.offset).map_err(|_| TestCaseError::reject("w"))?;
        let wider = Band::new(width + 8, case.offset).map_err(|_| TestCaseError::reject("w"))?;
        let a = align_banded(&case.haplotype, &case.read, &StandardEmission::default(), narrow).get();
        let b = align_banded(&case.haplotype, &case.read, &StandardEmission::default(), wider).get();
        let full = align_full(&case.haplotype, &case.read, &StandardEmission::default()).get();
        prop_assert!(!a.is_nan() && !b.is_nan(), "a {} b {}", a, b);
        if b > -30.0 {
            prop_assert!(a <= b + 1e-4, "narrow {} beats wider {}", a, b);
        }
        prop_assert!(b <= full + 1e-4, "banded {} beats the reference {}", b, full);
    }

    /// The kernel, with the band held fixed: the `f32` band is the `f64`
    /// recurrence over that same band, under **any** conversion model however
    /// implausible. Separating this from the gate above is what says a
    /// disagreement there is the band missing a path and not the arithmetic.
    #[test]
    fn the_f32_band_is_the_f64_recurrence_over_the_same_band(
        case in derived_case(3),
        conversion in any_conversion(),
    ) {
        let band = case.band();
        let taps = TapsEmission::new(conversion, Betas::PerSite(&case.betas));
        let got = align_banded(&case.haplotype, &case.read, &taps, band).get();
        let want = align_masked(&case.haplotype, &case.read, &taps, band);
        prop_assert!(
            (got - want).abs() < 1e-4 * (1.0 + want.abs()),
            "f32 {} against the f64 recurrence over the same band {}", got, want
        );
    }

    /// The same, at every width: `the_band_is_the_documented_predicate` above
    /// only compares *values* where a random read scores above `-30`, which
    /// is short reads. Reads cut from the haplotype score near zero at any
    /// length, so this is where a 100-base read meets a 3-column band and
    /// the `f32` numbers still have to be the `f64` recurrence's.
    #[test]
    fn the_f32_band_is_the_f64_recurrence_at_every_width(
        case in derived_case(3),
        width in 2u32..64,
        conversion in any_conversion(),
    ) {
        let band = Band::new(width, case.offset).map_err(|_| TestCaseError::reject("w"))?;
        let taps = TapsEmission::new(conversion, Betas::PerSite(&case.betas));
        let got = align_banded_simd(&case.haplotype, &case.read, &taps, band).get();
        let want = align_masked(&case.haplotype, &case.read, &taps, band);
        prop_assert_eq!(got.is_finite(), want.is_finite(), "kernel {} oracle {}", got, want);
        if got.is_finite() {
            prop_assert!(
                (got - want).abs() < 1e-4 * (1.0 + want.abs()),
                "width {}: f32 {} against the f64 recurrence over the same band {}", width, got, want
            );
        }
    }

    /// Bit-parity is a property of the `Lane` trait, so it must hold at every
    /// width, not only at `DEFAULT_WIDTH` -- a band narrower than one vector is
    /// where the tail-lane masking is exercised.
    #[test]
    fn simd_is_bit_identical_to_scalar_at_every_width(
        case in arbitrary_case(),
        width in 2u32..64,
        conversion in any_conversion(),
    ) {
        let band = Band::new(width, case.offset).map_err(|_| TestCaseError::reject("w"))?;
        let taps = TapsEmission::new(conversion, Betas::PerSite(&case.betas));
        prop_assert_eq!(
            align_banded(&case.haplotype, &case.read, &taps, band).get().to_bits(),
            align_banded_simd(&case.haplotype, &case.read, &taps, band).get().to_bits()
        );
    }
}

/// Moving the haplotype and the band by the same amount leaves the alignment
/// alone; only the `1 / haplotype_len` start prior changes, and by exactly the
/// length ratio. This is what says the band's offset arithmetic and its
/// `lo`-history bookkeeping do not depend on where in the haplotype the read
/// sits.
#[test]
fn shifting_the_haplotype_only_moves_the_start_prior() {
    let bases = support::bases(
        "TCATTGGCTATCCTAACCCGACCCTAGGAGCGGTTGGCGTGTATGCCGTGAATTTTCTCATTTCCGCTAGACATAATCGTTCTGCCTATA",
    );
    let haplotype = Haplotype::new(bases.clone());
    let read = Read::uniform(
        bases.get(40..75).expect("in range").to_vec(),
        &[BaseQuality::from_byte(32); 35],
        BaseQuality::from_byte(45),
        BaseQuality::from_byte(45),
        BaseQuality::from_byte(10),
        Strand::OT,
    )
    .expect("valid");

    let here =
        align_banded(&haplotype, &read, &StandardEmission::default(), Band::anchored(40)).get();
    for pad in [1usize, 7, 30, 64] {
        let mut shifted_bases = vec![Base::A; pad];
        shifted_bases.extend_from_slice(&bases);
        let shifted = Haplotype::new(shifted_bases);
        #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap, reason = "pad <= 64")]
        let band = Band::anchored(40 + pad as i32);
        let there = align_banded(&shifted, &read, &StandardEmission::default(), band).get();
        #[allow(clippy::cast_precision_loss, reason = "lengths under 200")]
        let expected = here + (haplotype.len() as f64 / shifted.len() as f64).log10();
        assert!(
            (there - expected).abs() < 1e-4,
            "pad {pad}: {there} against the expected {expected} (unpadded {here})"
        );
    }
}

/// `f32` with a power-of-two renormalisation must carry a read far longer than
/// any it will see in production without underflowing or drifting from the
/// `f64` reference.
#[test]
fn a_long_read_neither_underflows_nor_drifts() {
    let mut state = 0x1234_5678_9abc_def0_u64;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let bases: Vec<Base> = (0..4000)
        .map(|_| match next() % 4 {
            0 => Base::A,
            1 => Base::C,
            2 => Base::G,
            _ => Base::T,
        })
        .collect();
    let haplotype = Haplotype::new(bases.clone());
    let read_bases = bases.get(1000..2500).expect("in range").to_vec();
    let length = read_bases.len();
    let read = Read::uniform(
        read_bases,
        &vec![BaseQuality::from_byte(30); length],
        BaseQuality::from_byte(45),
        BaseQuality::from_byte(45),
        BaseQuality::from_byte(10),
        Strand::OT,
    )
    .expect("valid");

    let full = align_full(&haplotype, &read, &StandardEmission::default()).get();
    let banded =
        align_banded(&haplotype, &read, &StandardEmission::default(), Band::anchored(1000)).get();
    let simd =
        align_banded_simd(&haplotype, &read, &StandardEmission::default(), Band::anchored(1000))
            .get();
    assert!(full.is_finite() && banded.is_finite(), "full {full}, banded {banded}");
    assert_eq!(banded.to_bits(), simd.to_bits(), "banded {banded}, simd {simd}");
    assert!((banded - full).abs() < 1e-2, "banded {banded} against the reference {full}");
}

/// Nothing in the degenerate corners panics, and every one of them is
/// `IMPOSSIBLE` rather than a `NaN` a caller would carry into a likelihood
/// ratio.
#[test]
fn the_degenerate_corners_are_impossible_not_a_panic() {
    let haplotype = Haplotype::from_ascii(b"ACGT");
    let empty_haplotype = Haplotype::new(Vec::new());
    let read = Read::uniform(
        vec![Base::A],
        &[BaseQuality::from_byte(30)],
        BaseQuality::from_byte(45),
        BaseQuality::from_byte(45),
        BaseQuality::from_byte(10),
        Strand::OT,
    )
    .expect("valid");
    let empty_read = Read::uniform(
        Vec::new(),
        &[],
        BaseQuality::from_byte(45),
        BaseQuality::from_byte(45),
        BaseQuality::from_byte(10),
        Strand::OT,
    )
    .expect("valid");

    for (name, got) in [
        ("full, no haplotype", align_full(&empty_haplotype, &read, &StandardEmission::default())),
        ("full, no read", align_full(&haplotype, &empty_read, &StandardEmission::default())),
        (
            "banded, no haplotype",
            align_banded(&empty_haplotype, &read, &StandardEmission::default(), Band::anchored(0)),
        ),
        (
            "banded, no read",
            align_banded(&haplotype, &empty_read, &StandardEmission::default(), Band::anchored(0)),
        ),
        (
            "simd, no haplotype",
            align_banded_simd(
                &empty_haplotype,
                &read,
                &StandardEmission::default(),
                Band::anchored(0),
            ),
        ),
    ] {
        assert_eq!(got, compair::Log10Likelihood::IMPOSSIBLE, "{name}: {got:?}");
    }

    // Extreme offsets are arithmetic on `i64`, so neither end wraps.
    for offset in [i32::MIN, i32::MIN + 1, -1, 0, 1, i32::MAX - 1, i32::MAX] {
        let got =
            align_banded(&haplotype, &read, &StandardEmission::default(), Band::anchored(offset))
                .get();
        assert!(!got.is_nan(), "offset {offset} gave {got}");
    }

    // A quality of zero means "this base is certainly wrong", which leaves no
    // transition out of a match and no path.
    let certain = Read::uniform(
        haplotype.bases().to_vec(),
        &[BaseQuality::from_byte(0); 4],
        BaseQuality::from_byte(0),
        BaseQuality::from_byte(0),
        BaseQuality::from_byte(0),
        Strand::OT,
    )
    .expect("valid");
    assert_eq!(
        align_full(&haplotype, &certain, &StandardEmission::default()),
        compair::Log10Likelihood::IMPOSSIBLE
    );
}

/// `Band::new` rejects a width under two and **truncates** every other odd one:
/// the band it builds is `2 * (width / 2) + 1` columns wide, so widths `2 * n`
/// and `2 * n + 1` are the same band. Pinned because `width()` reports the
/// truncated value rather than the argument.
#[test]
fn band_new_truncates_odd_widths() {
    assert!(Band::new(0, 0).is_err());
    assert!(Band::new(1, 0).is_err());
    for (asked, reported) in [(2u32, 2u32), (3, 2), (4, 4), (47, 46), (48, 48), (49, 48)] {
        let band = Band::new(asked, 0).expect("width at least two");
        assert_eq!(band.width(), reported, "Band::new({asked}).width()");
    }
    assert_eq!(Band::new(2, 0).ok(), Band::new(3, 0).ok());
    assert_eq!(Band::anchored(7).width(), Band::DEFAULT_WIDTH);
    assert_eq!(Band::anchored(7).offset(), 7);
    assert_eq!(Band::anchored(7), Band::new(Band::DEFAULT_WIDTH, 7).expect("the default is valid"));
}

/// The one place the `f32` kernel is not the `f64` recurrence, and the band is
/// what keeps it out of trouble.
///
/// Each anti-diagonal is renormalised so its largest cell lands in `[1, 2)`.
/// The diagonal holds one cell per read row it crosses, and those are prefix
/// alignments of different lengths, so their magnitudes span the alignment's
/// whole dynamic range. A band of half-width `w` crosses at most `w + 1` rows,
/// which bounds that span; a band wider than the read crosses all of them,
/// and everything more than ~1e-38 below the diagonal's maximum flushes to
/// zero in `f32`.
///
/// Measured here: at `DEFAULT_WIDTH` the `f32` kernel is bit-for-bit the `f64`
/// recurrence over the same band even at `log10 L = -75`; widen the band past
/// the read and it loses several log10. The loss is always an under-estimate,
/// never an over-estimate, so it cannot make a wrong haplotype win.
#[test]
fn an_unbanded_f32_run_underflows_where_a_real_band_does_not() {
    let mut state = 0x5151_2323_9999_0f0f_u64;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let draw = |next: &mut dyn FnMut() -> u64, n: usize| -> Vec<Base> {
        (0..n)
            .map(|_| match next() % 4 {
                0 => Base::A,
                1 => Base::C,
                2 => Base::G,
                _ => Base::T,
            })
            .collect()
    };

    let mut worst_banded = 0.0f64;
    let mut worst_unbanded = 0.0f64;
    for _ in 0..200 {
        let haplotype = Haplotype::new(draw(&mut next, 90));
        // Unrelated to the haplotype, so the score is deeply negative and the
        // dynamic range across read rows is large.
        let read = Read::uniform(
            draw(&mut next, 66),
            &[BaseQuality::from_byte(30); 66],
            BaseQuality::from_byte(45),
            BaseQuality::from_byte(45),
            BaseQuality::from_byte(10),
            Strand::OT,
        )
        .expect("valid");

        for (width, worst) in [(Band::DEFAULT_WIDTH, &mut worst_banded), (320, &mut worst_unbanded)]
        {
            let band = Band::new(width, 0).expect("width at least two");
            let got = align_banded(&haplotype, &read, &StandardEmission::default(), band).get();
            let want = align_masked(&haplotype, &read, &StandardEmission::default(), band);
            if !got.is_finite() || !want.is_finite() {
                continue;
            }
            assert!(got <= want + 1e-6, "width {width}: f32 {got} over the f64 recurrence {want}");
            *worst = worst.max((got - want).abs());
        }
    }
    assert!(
        worst_banded < 1e-5,
        "at DEFAULT_WIDTH the f32 kernel tracks the f64 recurrence to rounding, worst {worst_banded}"
    );
    assert!(
        worst_unbanded > 1.0,
        "a band wider than the read is expected to underflow visibly; worst seen {worst_unbanded}"
    );
}

/// A band wide enough to hold the whole matrix is not a constraint, so it must
/// give back the unbanded reference. Run on GATK's own 164 x 101 vectors, whose
/// scores are a few log10 and so well inside `f32`'s range -- see
/// `an_unbanded_f32_run_underflows_where_a_real_band_does_not` for what happens
/// when they are not.
#[test]
fn a_band_that_holds_the_whole_matrix_is_not_a_constraint() {
    let vectors = support::gatk_vectors().expect("the bundled test data parses");
    let mut worst = 0.0f64;
    for (index, vector) in vectors.iter().enumerate() {
        let span = vector.haplotype.len() + vector.read.len() + 4;
        let width = u32::try_from(2 * span).expect("a few hundred");
        let band = Band::new(width, 0).expect("width at least two");
        let banded =
            align_banded(&vector.haplotype, &vector.read, &StandardEmission::default(), band).get();
        let full = align_full(&vector.haplotype, &vector.read, &StandardEmission::default()).get();
        assert!(
            (banded - full).abs() < 1e-5,
            "vector {index}: unconstrained band {banded} against the reference {full}"
        );
        worst = worst.max((banded - full).abs());
    }
    assert!(worst < 1e-5, "worst {worst}");
}

/// The regression the proptest above found, kept as a deterministic case
/// because it needs a band of width two at a negative offset and proptest only
/// reaches that shape every few thousand draws.
///
/// The read's free start is the `1 / haplotype_len` in `D[0][j]`, and the
/// kernel used to place it at `(0, 0)` whenever its row index clamped to zero
/// -- which happens for every `offset < -(width / 2)`, where the band holds no
/// `(0, j)` at all. The result was a finite score for an alignment the band
/// excludes.
#[test]
fn the_free_start_is_never_granted_outside_the_band() {
    let bases = support::bases("ACGTACGTACGTACGTACGTACGTACGTACGTACGTA");
    let haplotype = Haplotype::new(bases.clone());
    for read_length in [1usize, 2, 5] {
        let read = Read::uniform(
            bases.get(0..read_length).expect("in range").to_vec(),
            &vec![BaseQuality::from_byte(30); read_length],
            BaseQuality::from_byte(45),
            BaseQuality::from_byte(45),
            BaseQuality::from_byte(10),
            Strand::OT,
        )
        .expect("valid");
        for (width, offset) in [(2u32, -2i32), (2, -10), (4, -3), (26, -15)] {
            let band = Band::new(width, offset).expect("width at least two");
            let got = align_banded(&haplotype, &read, &StandardEmission::default(), band);
            let simd = align_banded_simd(&haplotype, &read, &StandardEmission::default(), band);
            assert_eq!(
                got,
                compair::Log10Likelihood::IMPOSSIBLE,
                "width {width} at offset {offset} holds no (0, j): read {read_length} got {got:?}"
            );
            assert_eq!(got, simd);
        }
    }
}

/// The other half of the out-of-band contract, and the half the `Band` doc
/// comment used to get wrong: a band that contains no path at all scores
/// `IMPOSSIBLE`, not a small finite number.
///
/// `align_banded` grants the read's free start only at cells `(0, j)` the band
/// contains. Push the band off the left of the haplotype -- `offset` below
/// `-(width / 2)`, which is what an anchor from a read hanging off the start of
/// a short haplotype looks like -- and there is no such cell, so no path opens.
/// A caller comparing haplotypes under one band has to expect `-inf` here, and
/// not read it as a failure.
#[test]
fn a_band_with_no_start_cell_is_impossible() {
    let bases = support::bases(
        "TCATTGGCTATCCTAACCCGACCCTAGGAGCGGTTGGCGTGTATGCCGTGAATTTTCTCATTTCCGCTAGACATAATCGTTCTGCCTATA",
    );
    let haplotype = Haplotype::new(bases.clone());
    let read = Read::uniform(
        bases.get(0..40).expect("in range").to_vec(),
        &[BaseQuality::from_byte(32); 40],
        BaseQuality::from_byte(45),
        BaseQuality::from_byte(45),
        BaseQuality::from_byte(10),
        Strand::OT,
    )
    .expect("valid");

    let half = i32::try_from(Band::DEFAULT_WIDTH / 2).expect("24 fits");
    for offset in [-half, -half + 1, 0] {
        let got =
            align_banded(&haplotype, &read, &StandardEmission::default(), Band::anchored(offset))
                .get();
        assert!(got.is_finite(), "offset {offset} still holds a start cell and must score");
    }
    for offset in [-half - 1, -half - 10, -1_000, i32::MIN] {
        let got =
            align_banded(&haplotype, &read, &StandardEmission::default(), Band::anchored(offset));
        let simd = align_banded_simd(
            &haplotype,
            &read,
            &StandardEmission::default(),
            Band::anchored(offset),
        );
        assert_eq!(
            got,
            compair::Log10Likelihood::IMPOSSIBLE,
            "offset {offset} holds no (0, j) cell, got {got:?}"
        );
        assert_eq!(got, simd);
    }
}

/// Found by `fuzz_pair_hmm`: a read at Phred 0 whose last base is `N`, so
/// every match emission is exactly zero and the last row is reached through
/// the insertion state alone, against a haplotype that is mostly `N`.
///
/// On the diagonal where the free-start cell leaves the band, the only cell
/// left is the tail of a long deletion chain, `gap^n`, and the diagonal's
/// maximum collapses from `~1` to `~1e-37` in one step. That diagonal asks
/// the scaling for a shift of 123, the next one for 10, and a kernel that
/// composed the two shifts into one factor for the diagonal two back
/// computed `2^133`: infinity, then `0 * inf = NaN`, then `IMPOSSIBLE` for a
/// pair the reference scores at `-5.7`. The two lifts are applied one at a
/// time now, and this pins it at every width that holds a path.
#[test]
fn a_collapsing_diagonal_maximum_does_not_overflow_the_lift() {
    let mut sequence = b"gg0gg".to_vec();
    sequence.extend([255u8; 13]);
    sequence.extend(b"U");
    sequence.extend([255u8; 3]);
    sequence.extend(b"g;g(gggg");
    sequence.extend([255u8; 5]);
    sequence.extend(b"G");
    sequence.extend([255u8; 3]);
    sequence.extend([3u8]);
    sequence.extend(b"cc");
    sequence.extend([248u8]);
    sequence.extend([255u8; 9]);
    sequence.extend(b"GGG");
    let haplotype = Haplotype::from_ascii(&sequence);
    assert_eq!(haplotype.len(), 55);
    let read = Read::uniform(
        vec![Base::G, Base::G, Base::Unknown],
        &[BaseQuality::from_byte(0); 3],
        BaseQuality::from_byte(30),
        BaseQuality::from_byte(30),
        BaseQuality::from_byte(30),
        Strand::OB,
    )
    .expect("valid");
    let emission = StandardEmission::default();
    let full = align_full(&haplotype, &read, &emission).get();
    assert!((full - -5.6933).abs() < 1e-3, "the reference scores this pair at -5.69, got {full}");
    // The path ends at columns 41 and 42, so a band at offset 0 needs a
    // half-width of at least 39 to hold it at all; at 40 it holds the path
    // but not everything around it, and scores a little under the reference.
    let narrow = Band::new(80, 0).expect("valid width");
    assert!(align_banded(&haplotype, &read, &emission, narrow).get().is_finite());
    for width in [124u32, 200, Band::MAX_WIDTH] {
        let band = Band::new(width, 0).expect("valid width");
        let banded = align_banded(&haplotype, &read, &emission, band).get();
        let simd = align_banded_simd(&haplotype, &read, &emission, band).get();
        assert_eq!(banded.to_bits(), simd.to_bits());
        assert!(
            (banded - full).abs() < 1e-4,
            "width {width}: banded {banded} against the reference {full}"
        );
    }
}

/// `Band::MAX_WIDTH` is a bound on the allocation, and the error is typed.
///
/// The kernels allocate nine buffers of `width / 2 + O(1)` floats, so an
/// unbounded width is an unbounded allocation inside a safe function. The
/// rejection happens in `Band::new`, before any of it: this test never reaches
/// a kernel.
#[test]
fn band_new_rejects_a_width_it_would_have_to_allocate() {
    assert_eq!(Band::new(0, 0), Err(compair::Error::BandTooNarrow(0)));
    assert_eq!(Band::new(1, 0), Err(compair::Error::BandTooNarrow(1)));
    assert!(Band::new(Band::MAX_WIDTH, 0).is_ok());
    for width in [Band::MAX_WIDTH + 1, 1 << 20, u32::MAX] {
        assert_eq!(
            Band::new(width, 0),
            Err(compair::Error::BandTooWide { got: width, max: Band::MAX_WIDTH }),
            "width {width}"
        );
    }
    const { assert!(Band::DEFAULT_WIDTH <= Band::MAX_WIDTH) };
}
