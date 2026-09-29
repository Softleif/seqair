#[path = "../src/pinned.rs"]
mod pinned;
#[allow(dead_code, reason = "each integration test uses a different part of this")]
mod support;

use compair::{
    BATCH, Band, Base, BaseQuality, Betas, ConversionModel, Haplotype, MatchProbability,
    Probability, Read, StandardEmission, Strand, TapsEmission, Workspace, align_banded,
    align_banded_simd, align_candidates, align_full, align_strips, align_strips_simd,
    error_probability,
};
use hegel::TestCase;
use hegel::generators::{self as gs, Generator};
use support::{
    any_base_or_n, any_conversion, any_probability, arbitrary_case, bases, derived_case, length,
    plausible_conversion,
};

// The bit-parity gates are the cheapest of these and the ones most worth
// running hard, so they pin 2048 cases where the rest take fewer: a
// full-matrix oracle at 2048 cases was most of the suite's wall-clock.

/// The C4 gate: the scalar and the eight-wide kernels agree to the bit, on
/// inputs that include bands missing the alignment entirely, because parity
/// has to hold there too. And a reused `Workspace` is the same kernel:
/// nothing from the previous pair leaks into the next.
#[hegel::test(test_cases = pinned::cases(2048))]
fn simd_is_bit_identical_to_scalar(tc: TestCase) {
    let case = tc.draw(arbitrary_case());
    let conversion = tc.draw(any_conversion());
    let uniform = tc.draw(any_probability());
    let band = case.band();
    let taps = TapsEmission::new(conversion, Betas::PerSite(&case.betas));
    let uniform = TapsEmission::new(conversion, Betas::Uniform(uniform));
    let mut workspace = Workspace::new();
    for (name, banded, simd, reused) in [
        (
            "standard",
            align_banded(&case.haplotype, &case.read, &StandardEmission::default(), band),
            align_banded_simd(&case.haplotype, &case.read, &StandardEmission::default(), band),
            workspace.align_banded_simd(
                &case.haplotype,
                &case.read,
                &StandardEmission::default(),
                band,
            ),
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
        assert_eq!(
            banded.get().to_bits(),
            simd.get().to_bits(),
            "{}: scalar {:?} vs simd {:?}",
            name,
            banded,
            simd
        );
        assert_eq!(
            banded.get().to_bits(),
            reused.get().to_bits(),
            "{}: fresh {:?} vs reused workspace {:?}",
            name,
            banded,
            reused
        );
    }
}

/// The same gate for the strip kernel: its scalar and eight-lane
/// instances agree to the bit, fresh and through a reused `Workspace`,
/// and one workspace serves both traversals in any order.
#[hegel::test(test_cases = pinned::cases(2048))]
fn strips_simd_is_bit_identical_to_strips_scalar(tc: TestCase) {
    let case = tc.draw(arbitrary_case());
    let conversion = tc.draw(any_conversion());
    let uniform = tc.draw(any_probability());
    let band = case.band();
    let taps = TapsEmission::new(conversion, Betas::PerSite(&case.betas));
    let uniform = TapsEmission::new(conversion, Betas::Uniform(uniform));
    let mut workspace = Workspace::new();
    for (name, scalar, simd, reused) in [
        (
            "standard",
            align_strips(&case.haplotype, &case.read, &StandardEmission::default(), band),
            align_strips_simd(&case.haplotype, &case.read, &StandardEmission::default(), band),
            workspace.align_strips_simd(
                &case.haplotype,
                &case.read,
                &StandardEmission::default(),
                band,
            ),
        ),
        (
            "taps",
            align_strips(&case.haplotype, &case.read, &taps, band),
            align_strips_simd(&case.haplotype, &case.read, &taps, band),
            workspace.align_strips(&case.haplotype, &case.read, &taps, band),
        ),
        (
            "uniform",
            align_strips(&case.haplotype, &case.read, &uniform, band),
            align_strips_simd(&case.haplotype, &case.read, &uniform, band),
            {
                // A diagonal alignment in between: the two kernels share
                // the workspace's plan and must not share its state.
                workspace.align_banded_simd(&case.haplotype, &case.read, &taps, band);
                workspace.align_strips_simd(&case.haplotype, &case.read, &uniform, band)
            },
        ),
    ] {
        assert_eq!(
            scalar.get().to_bits(),
            simd.get().to_bits(),
            "{}: scalar {:?} vs simd {:?}",
            name,
            scalar,
            simd
        );
        assert_eq!(
            scalar.get().to_bits(),
            reused.get().to_bits(),
            "{}: fresh {:?} vs reused workspace {:?}",
            name,
            scalar,
            reused
        );
    }
}

/// The three eight-lane kernels at every SIMD level this CPU has, the
/// scalar `Fallback` included, against their scalar oracles: the diagonal
/// and strip kernels against their own scalar instances, and a ragged
/// batch against the scalar strip kernel one haplotype at a time. One
/// workspace serves every level in turn, so nothing a level leaves behind
/// may leak into the next.
#[hegel::test(test_cases = pinned::cases(512))]
fn every_simd_level_is_bit_identical_to_scalar(tc: TestCase) {
    let case = tc.draw(arbitrary_case());
    let conversion = tc.draw(any_conversion());
    let spread = tc.draw(gs::integers::<usize>().max_value(3));
    let band = case.band();
    let taps = TapsEmission::new(conversion, Betas::PerSite(&case.betas));
    let bases = case.haplotype.bases().to_vec();
    let batch: Vec<Haplotype> = (0..BATCH)
        .map(|k| Haplotype::new(bases[..bases.len().saturating_sub(k * spread)].to_vec()))
        .collect();
    let refs: Vec<&Haplotype> = batch.iter().collect();
    let banded = align_banded(&case.haplotype, &case.read, &taps, band).get().to_bits();
    let strips = align_strips(&case.haplotype, &case.read, &taps, band).get().to_bits();
    let one_at_a_time: Vec<u64> =
        refs.iter().map(|h| align_strips(h, &case.read, &taps, band).get().to_bits()).collect();
    let mut workspace = Workspace::new();
    let mut out = Vec::new();
    for (name, level) in support::levels() {
        let simd = workspace.align_banded_simd_at(level, &case.haplotype, &case.read, &taps, band);
        assert_eq!(banded, simd.get().to_bits(), "diagonal at {}", name);
        let simd = workspace.align_strips_simd_at(level, &case.haplotype, &case.read, &taps, band);
        assert_eq!(strips, simd.get().to_bits(), "strips at {}", name);
        workspace.align_batch_at(level, &refs, &case.read, &taps, band, &mut out);
        let batched: Vec<u64> = out.iter().map(|score| score.get().to_bits()).collect();
        assert_eq!(&one_at_a_time, &batched, "batch at {}", name);
    }
}

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
///
/// Nor does it hold for a haplotype that repeats itself: in a homopolymer
/// the read fits every offset equally well, most of the full matrix's mass
/// is outside any band, and the band is right to miss it. Those cases are
/// set aside rather than toleranced.
#[hegel::test]
fn banded_matches_the_reference_inside_the_band(tc: TestCase) {
    let case = tc.draw(derived_case(3));
    let conversion = tc.draw(plausible_conversion());
    tc.assume(!echoes_outside_the_band(&case));
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
        assert!((banded - full).abs() < 1e-3, "{name}: banded {banded} against reference {full}");
    }
}

/// Whether the read, placed without gaps at some offset the default band
/// leaves out, matches three in four of its bases there. A read cut from a
/// random haplotype matches about two in five anywhere but where it was
/// cut; one that matches this well elsewhere has a second alignment of real
/// weight, and the band holds only the first.
#[allow(
    clippy::cast_possible_wrap,
    reason = "test scaffolding over sequences of a few hundred bases"
)]
fn echoes_outside_the_band(case: &support::Case) -> bool {
    let (haplotype, read) = (case.haplotype.bases(), case.read.bases());
    // Three edits move the true diagonal by at most three columns.
    let reach = i64::from(Band::DEFAULT_WIDTH / 2) - 3;
    let seed = i64::from(case.offset);
    (-(read.len() as i64)..haplotype.len() as i64)
        .filter(|offset| (offset - seed).abs() > reach)
        .any(|offset| {
            let matches = read
                .iter()
                .enumerate()
                .filter(|&(index, &base)| {
                    usize::try_from(offset + index as i64)
                        .ok()
                        .and_then(|column| haplotype.get(column))
                        .is_some_and(|&hap| {
                            hap == base || hap == Base::Unknown || base == Base::Unknown
                        })
                })
                .count();
            4 * matches >= 3 * read.len()
        })
}

/// Where the chemistry cannot act -- a haplotype of `A` and `T` only, so
/// no cytosine on either strand -- `TapsEmission` is `StandardEmission`
/// through the whole DP, not only per cell.
///
/// This used to admit `C` as well, on the reading that the conversion rows
/// were `CpG`-only. They are not: the joint model puts every unmethylated
/// cytosine on the false-conversion row, so a non-`CpG` `C` on an OT read
/// differs from a plain mismatch by `f`. See
/// `a_non_cpg_cytosine_reads_t_at_the_false_conversion_rate`.
#[hegel::test]
fn taps_equals_standard_where_the_chemistry_cannot_act(tc: TestCase) {
    let len = tc.draw(length(30, 89));
    let bases = tc.draw(
        gs::vecs(gs::sampled_from(&[Base::A, Base::T]).print_as_debug())
            .min_size(len)
            .max_size(len),
    );
    let len = tc.draw(length(30, 89));
    let quals = tc.draw(
        gs::vecs(gs::integers::<u8>().min_value(2).max_value(45)).min_size(len).max_size(len),
    );
    let strand = tc.draw(support::any_strand());
    let conversion = tc.draw(any_conversion());
    let haplotype = Haplotype::new(bases.clone());
    let read_bases: Vec<Base> = bases.iter().copied().take(quals.len().min(bases.len())).collect();
    let quals: Vec<BaseQuality> =
        quals.iter().take(read_bases.len()).map(|q| BaseQuality::from_byte(*q)).collect();
    let read = Read::uniform(
        read_bases,
        &quals,
        BaseQuality::from_byte(45),
        BaseQuality::from_byte(45),
        BaseQuality::from_byte(10),
        strand,
    )
    .expect("read");
    let betas = vec![Probability::ONE; haplotype.len()];
    let taps = TapsEmission::new(conversion, Betas::PerSite(&betas));
    assert_eq!(
        align_full(&haplotype, &read, &taps).get().to_bits(),
        align_full(&haplotype, &read, &StandardEmission::default()).get().to_bits()
    );
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
///
/// Gap-open qualities below Q6 count as Q6, as GATK squashes them. Each row
/// is rescaled by a power of two, which is exact, so the oracle holds every
/// score `f64` can write down however far below zero it is.
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
    let gap_open =
        |track: &[BaseQuality], index: usize| quality(track, index).min(10f64.powf(-0.6));
    let mut log10_scale = 0.0f64;

    #[allow(clippy::cast_precision_loss, reason = "test haplotypes are short")]
    let init = 1.0 / h as f64;
    let mut prev_m = vec![0.0f64; h + 1];
    let mut prev_i = vec![0.0f64; h + 1];
    let mut prev_d: Vec<f64> = (0..=h).map(|j| if live(0, j) { init } else { 0.0 }).collect();
    let (mut cur_m, mut cur_i, mut cur_d) =
        (vec![0.0f64; h + 1], vec![0.0f64; h + 1], vec![0.0f64; h + 1]);

    for i in 1..=r {
        let Some(observation) = read.observation(i - 1) else { return f64::NEG_INFINITY };
        let p_ins = gap_open(read.insertion_quals(), i - 1);
        let p_del = gap_open(read.deletion_quals(), i - 1);
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
        let max = cur_m.iter().chain(&cur_i).chain(&cur_d).fold(0.0f64, |a, &b| a.max(b));
        if max > 0.0 {
            let exponent = max.log2().floor();
            let factor = 2f64.powf(-exponent);
            for cell in cur_m.iter_mut().chain(cur_i.iter_mut()).chain(cur_d.iter_mut()) {
                *cell *= factor;
            }
            log10_scale += exponent * core::f64::consts::LOG10_2;
        }
        core::mem::swap(&mut prev_m, &mut cur_m);
        core::mem::swap(&mut prev_i, &mut cur_i);
        core::mem::swap(&mut prev_d, &mut cur_d);
    }
    let total: f64 = prev_m.iter().zip(prev_i.iter()).skip(1).map(|(m, i)| m + i).sum();
    if total <= 0.0 { f64::NEG_INFINITY } else { total.log10() + log10_scale }
}

/// Every entry point, on one case: the `f64` banded recurrence, the scalar
/// and SIMD diagonal and strip kernels, the batch kernel (the haplotype eight
/// times) and the pairs kernel (the read eight times).
fn every_entry_point<E: compair::Emission>(
    case: &support::Case,
    emission: &E,
    band: Band,
) -> Vec<(&'static str, f64)> {
    let (haplotype, read) = (&case.haplotype, &case.read);
    let mut workspace = Workspace::new();
    let mut batch = Vec::new();
    workspace.align_candidates(&[haplotype; BATCH], read, emission, band, &mut batch);
    let mut pairs = Vec::new();
    workspace.align_reads(&[haplotype], &[(read, band); compair::PAIRS], emission, &mut pairs);
    let mut all = vec![
        ("f64", compair::align_banded_f64(haplotype, read, emission, band).get()),
        ("diagonals", align_banded(haplotype, read, emission, band).get()),
        ("diagonals/simd", align_banded_simd(haplotype, read, emission, band).get()),
        ("strips", align_strips(haplotype, read, emission, band).get()),
        ("strips/simd", align_strips_simd(haplotype, read, emission, band).get()),
    ];
    all.extend(batch.iter().map(|score| ("batch", score.get())));
    all.extend(pairs.iter().map(|score| ("pairs", score.get())));
    all
}

/// **Every score is the `f64` recurrence over its band**, at every quality a
/// read can carry and however far below zero the score is -- not only where
/// a random read happens to score above `-30`. The `f32` kernels cannot hold
/// every such score themselves: a cell more than `2^-126` below the row (or
/// diagonal) it is scaled by flushes to zero, which loses nothing measurable
/// while the total stays well above that, and everything once it does not.
/// So an entry point that finishes below the level where its flushes could
/// matter scores the pair again in `f64`, and that is what this holds.
#[hegel::test(test_cases = pinned::cases(1024))]
fn every_score_is_the_f64_recurrence_over_its_band(tc: TestCase) {
    let case = tc.draw(arbitrary_case());
    let width = tc.draw(gs::integers::<u32>().min_value(2).max_value(63));
    let conversion = tc.draw(any_conversion());
    let band = Band::new(width, case.offset).expect("width");
    let taps = TapsEmission::new(conversion, Betas::PerSite(&case.betas));
    for (emission, want, got) in [
        (
            "standard",
            align_masked(&case.haplotype, &case.read, &StandardEmission::default(), band),
            every_entry_point(&case, &StandardEmission::default(), band),
        ),
        (
            "taps",
            align_masked(&case.haplotype, &case.read, &taps, band),
            every_entry_point(&case, &taps, band),
        ),
    ] {
        for (kernel, got) in got {
            assert_eq!(
                got.is_finite(),
                want.is_finite(),
                "{emission}/{kernel}: {got} against the f64 recurrence {want}"
            );
            if want.is_finite() {
                assert!(
                    (got - want).abs() < 1e-4 * (1.0 + want.abs()),
                    "{emission}/{kernel}: {got} against the f64 recurrence {want}"
                );
            }
        }
    }
}

/// The strip kernels keep a row's largest cell at `2^STRIP_SCALE` (`2^96`),
/// so what their flushes can lose, and the floor below which their scores
/// are rescored, sits 96 bits lower than the diagonal kernel's. This pair
/// (hegel's, shrunk) has the strip kernel flush a cell that goes on to carry
/// part of the answer: its `f32` score is -69.268 against the recurrence's
/// -69.258, twelve decades below the strip floor (-57.4). A floor that took
/// the scale for more than it is -- here, by forty bits -- lets that score
/// through; this pins it.
#[test]
fn a_strip_flush_below_the_strip_floor_is_rescored() {
    let quals = |qs: &[u8]| qs.iter().map(|&q| BaseQuality::from_byte(q)).collect::<Vec<_>>();
    let gaps = quals(&[3, 2, 6, 0, 4, 86, 254, 95, 213, 254, 53, 254, 210, 254, 254, 29, 64]);
    let read = Read::new(
        bases("TNNANAGCGAANCCCAC"),
        &quals(&[0, 0, 0, 0, 0, 3, 61, 51, 85, 1, 0, 6, 237, 248, 214, 9, 56]),
        &gaps,
        &gaps,
        &gaps,
        Strand::OB,
    )
    .expect("a valid read");
    let case = support::Case {
        haplotype: Haplotype::new(bases("AANANAAACNAANCCAANAANAAAANGA")),
        read,
        offset: 19,
        betas: vec![Probability::ZERO; 28],
    };
    let band = Band::new(27, case.offset).expect("width");
    let emission = StandardEmission::default();
    let want = align_masked(&case.haplotype, &case.read, &emission, band);
    assert!((want - -69.2579).abs() < 1e-3, "{want}");
    for (kernel, got) in every_entry_point(&case, &emission, band) {
        assert!((got - want).abs() < 1e-4 * (1.0 + want.abs()), "{kernel}: {got} against {want}");
    }
}

/// A read whose every base is below Q2 says nothing about which haplotype it
/// came from: every base is equally likely at every column, so every
/// haplotype of one length scores the same, to the bit, through every
/// kernel. (At Q0 as a literal `eps = 1` a base that *matches* scores zero,
/// and such a read preferred whichever haplotype it matched least.)
#[hegel::test(test_cases = pinned::cases(256))]
fn a_read_with_no_information_scores_every_haplotype_alike(tc: TestCase) {
    let h = tc.draw(length(1, 60));
    let r = tc.draw(length(1, 60));
    let draw =
        |n: usize| -> Vec<Base> { (0..n).map(|_| tc.draw_silent(any_base_or_n())).collect() };
    let (one, other, bases) = (Haplotype::new(draw(h)), Haplotype::new(draw(h)), draw(r));
    let quals: Vec<BaseQuality> = (0..r)
        .map(|_| BaseQuality::from_byte(tc.draw(gs::integers::<u8>().max_value(1))))
        .collect();
    let gap_open = BaseQuality::from_byte(tc.draw(gs::integers::<u8>().min_value(6).max_value(60)));
    let read =
        Read::uniform(bases, &quals, gap_open, gap_open, BaseQuality::from_byte(10), Strand::OT)
            .expect("valid");
    let offset = tc.draw(gs::integers::<i32>().min_value(-10).max_value(30));
    let band = Band::new(tc.draw(gs::integers::<u32>().min_value(2).max_value(63)), offset)
        .expect("width");
    let emission = StandardEmission::default();
    let case = |haplotype: &Haplotype| support::Case {
        haplotype: haplotype.clone(),
        read: read.clone(),
        offset,
        betas: Vec::new(),
    };
    let (a, b) = (case(&one), case(&other));
    for ((kernel, x), (_, y)) in every_entry_point(&a, &emission, band)
        .into_iter()
        .zip(every_entry_point(&b, &emission, band))
    {
        assert_eq!(x.to_bits(), y.to_bits(), "{kernel}: {x} against {y}");
    }
    let (x, y) =
        (align_full(&one, &read, &emission).get(), align_full(&other, &read, &emission).get());
    assert!((x - y).abs() <= 1e-12 * (1.0 + x.abs()), "full: {x} against {y}");
}

/// A cell of GATK's recurrence is not a probability: the transitions out of
/// a match are taken from two rows, as are those out of a deletion, and
/// where the rows' qualities differ a cell hands on more than it holds. A
/// read whose gap-continuation quality alternates between Q254 and Q0 scores
/// far above zero against a homopolymer, and the strip kernel's cells grow
/// by about `2^20` every eight rows: at `2^115` that overflows `f32` (the
/// kernel's own total is NaN). The read's `Growth` says so before the kernel
/// runs, and every entry point returns the `f64` recurrence.
#[test]
fn a_read_whose_rows_hand_on_more_than_they_hold_is_scored_in_f64() {
    let q = BaseQuality::from_byte;
    let r = 64;
    let gaps: Vec<BaseQuality> = (0..r).map(|i| if i % 2 == 0 { q(254) } else { q(0) }).collect();
    let read =
        Read::new(vec![Base::A; r], &vec![q(60); r], &vec![q(6); r], &vec![q(6); r], &gaps, Strand::OT)
            .expect("a valid read");
    let case = support::Case {
        haplotype: Haplotype::new(vec![Base::A; 1000]),
        read,
        offset: 0,
        betas: vec![],
    };
    let band = Band::new(Band::MAX_WIDTH, case.offset).expect("width");
    let emission = StandardEmission::default();
    let want = align_masked(&case.haplotype, &case.read, &emission, band);
    assert!(want > 20.0, "{want}");
    for (kernel, got) in every_entry_point(&case, &emission, band) {
        assert!((got - want).abs() < 1e-9 * (1.0 + want.abs()), "{kernel}: {got} against {want}");
    }
}

/// A read that fits nowhere, with qualities a real run reports: eighty `C`s
/// at Q40 against a homopolymer of `A`s. Every path pays at least a Q40
/// mismatch or a Q60 gap per row, so eight rows take every cell below
/// `2^-126` of the row the kernels last scaled by, and the `f32` kernels
/// returned `-inf` for a pair the `f64` recurrence scores near `-440`.
#[test]
fn a_read_that_fits_nowhere_still_has_a_score() {
    let emission = StandardEmission::default();
    for h in [9, 17, 40] {
        let haplotype = Haplotype::new(vec![Base::A; h]);
        let read = Read::uniform(
            vec![Base::C; 80],
            &[BaseQuality::from_byte(40); 80],
            BaseQuality::from_byte(60),
            BaseQuality::from_byte(60),
            BaseQuality::from_byte(60),
            Strand::OT,
        )
        .expect("valid");
        let band = Band::new(Band::MAX_WIDTH, 0).expect("width");
        let want = align_full(&haplotype, &read, &emission).get();
        let case = support::Case {
            haplotype: haplotype.clone(),
            read: read.clone(),
            offset: 0,
            betas: vec![],
        };
        for (kernel, got) in every_entry_point(&case, &emission, band) {
            assert!(
                (got - want).abs() < 1e-4 * (1.0 + want.abs()),
                "haplotype {h}: {kernel} {got} against the full reference {want}"
            );
        }
    }
}

/// The band is exactly the strip `Band` documents: `|j - i - offset| <=
/// width / 2` and nothing else. `arbitrary_case` reaches offsets outside
/// the matrix on both sides, which is where the kernel's clamps live.
#[hegel::test(test_cases = pinned::cases(1024))]
fn the_band_is_the_documented_predicate(tc: TestCase) {
    let case = tc.draw(arbitrary_case());
    let width = tc.draw(gs::integers::<u32>().min_value(2).max_value(63));
    let conversion = tc.draw(any_conversion());
    let band = Band::new(width, case.offset).expect("width");
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
        (
            "strips/standard",
            align_strips(&case.haplotype, &case.read, &StandardEmission::default(), band).get(),
            align_masked(&case.haplotype, &case.read, &StandardEmission::default(), band),
        ),
        (
            "strips/taps",
            align_strips_simd(&case.haplotype, &case.read, &taps, band).get(),
            align_masked(&case.haplotype, &case.read, &taps, band),
        ),
    ] {
        assert_eq!(
            got.is_finite(),
            want.is_finite(),
            "{}: kernel {} masked oracle {}",
            name,
            got,
            want
        );
        // The finite/infinite agreement above is the geometry claim -- it is
        // what caught B1. The value holds at any score: where the `f32`
        // kernel's range runs out the entry point rescores in `f64`.
        if want.is_finite() {
            assert!(
                (got - want).abs() < 1e-4 * (1.0 + want.abs()),
                "{}: kernel {} masked oracle {}",
                name,
                got,
                want
            );
        }
    }
}

/// Widening a band can only add paths, so it can only raise the score, and
/// no band can beat the unbanded reference.
///
/// Both halves hold at any score. The `f32` kernels only ever lose mass,
/// to rounding and to the subnormals they flush; a wider band's diagonals
/// span more dynamic range and can flush mass a narrower band kept -- the
/// fuzzer found a pair at `-133` where eight more columns cost 0.16 -- but
/// that is below the floor where the entry points rescore in `f64`.
#[hegel::test(test_cases = pinned::cases(1024))]
fn a_band_is_monotone_in_its_width_and_never_beats_the_reference(tc: TestCase) {
    let case = tc.draw(arbitrary_case());
    let width = tc.draw(gs::integers::<u32>().min_value(2).max_value(39));
    let narrow = Band::new(width, case.offset).expect("w");
    let wider = Band::new(width + 8, case.offset).expect("w");
    let full = align_full(&case.haplotype, &case.read, &StandardEmission::default()).get();
    for (name, a, b) in [
        (
            "diagonals",
            align_banded(&case.haplotype, &case.read, &StandardEmission::default(), narrow).get(),
            align_banded(&case.haplotype, &case.read, &StandardEmission::default(), wider).get(),
        ),
        (
            "strips",
            align_strips_simd(&case.haplotype, &case.read, &StandardEmission::default(), narrow)
                .get(),
            align_strips_simd(&case.haplotype, &case.read, &StandardEmission::default(), wider)
                .get(),
        ),
    ] {
        assert!(!a.is_nan() && !b.is_nan(), "{}: a {} b {}", name, a, b);
        assert!(
            a <= b || a - b <= 1e-4 * (1.0 + b.abs()),
            "{}: narrow {} beats wider {}",
            name,
            a,
            b
        );
        assert!(
            b <= full || b - full <= 1e-4 * (1.0 + full.abs()),
            "{}: banded {} beats the reference {}",
            name,
            b,
            full
        );
    }
}

/// The kernel, with the band held fixed: the `f32` band is the `f64`
/// recurrence over that same band, under **any** conversion model however
/// implausible. Separating this from the gate above is what says a
/// disagreement there is the band missing a path and not the arithmetic.
#[hegel::test(test_cases = pinned::cases(1024))]
fn the_f32_band_is_the_f64_recurrence_over_the_same_band(tc: TestCase) {
    let case = tc.draw(derived_case(3));
    let conversion = tc.draw(any_conversion());
    let band = case.band();
    let taps = TapsEmission::new(conversion, Betas::PerSite(&case.betas));
    let want = align_masked(&case.haplotype, &case.read, &taps, band);
    for (name, got) in [
        ("diagonals", align_banded(&case.haplotype, &case.read, &taps, band).get()),
        ("strips", align_strips(&case.haplotype, &case.read, &taps, band).get()),
    ] {
        assert!(
            (got - want).abs() < 1e-4 * (1.0 + want.abs()),
            "{}: f32 {} against the f64 recurrence over the same band {}",
            name,
            got,
            want
        );
    }
}

/// The same, at every width, on reads cut from the haplotype: they score
/// near zero at any length, so this is where a 100-base read meets a
/// 3-column band and the `f32` kernels' own numbers -- not a rescue --
/// have to be the `f64` recurrence's.
#[hegel::test(test_cases = pinned::cases(1024))]
fn the_f32_band_is_the_f64_recurrence_at_every_width(tc: TestCase) {
    let case = tc.draw(derived_case(3));
    let width = tc.draw(gs::integers::<u32>().min_value(2).max_value(63));
    let conversion = tc.draw(any_conversion());
    let band = Band::new(width, case.offset).expect("w");
    let taps = TapsEmission::new(conversion, Betas::PerSite(&case.betas));
    let want = align_masked(&case.haplotype, &case.read, &taps, band);
    for (name, got) in [
        ("diagonals", align_banded_simd(&case.haplotype, &case.read, &taps, band).get()),
        ("strips", align_strips_simd(&case.haplotype, &case.read, &taps, band).get()),
    ] {
        assert_eq!(got.is_finite(), want.is_finite(), "{}: kernel {} oracle {}", name, got, want);
        if got.is_finite() {
            assert!(
                (got - want).abs() < 1e-4 * (1.0 + want.abs()),
                "{} at width {}: f32 {} against the f64 recurrence over the same band {}",
                name,
                width,
                got,
                want
            );
        }
    }
}

/// Bit-parity is a property of the `Lane` trait, so it must hold at every
/// width, not only at `DEFAULT_WIDTH` -- a band narrower than one vector is
/// where the tail-lane masking is exercised.
#[hegel::test(test_cases = pinned::cases(1024))]
fn simd_is_bit_identical_to_scalar_at_every_width(tc: TestCase) {
    let case = tc.draw(arbitrary_case());
    let width = tc.draw(gs::integers::<u32>().min_value(2).max_value(63));
    let conversion = tc.draw(any_conversion());
    let band = Band::new(width, case.offset).expect("w");
    let taps = TapsEmission::new(conversion, Betas::PerSite(&case.betas));
    assert_eq!(
        align_banded(&case.haplotype, &case.read, &taps, band).get().to_bits(),
        align_banded_simd(&case.haplotype, &case.read, &taps, band).get().to_bits()
    );
    assert_eq!(
        align_strips(&case.haplotype, &case.read, &taps, band).get().to_bits(),
        align_strips_simd(&case.haplotype, &case.read, &taps, band).get().to_bits()
    );
}

/// The two traversals compute the same `f32` recurrence over the same
/// band, so they agree wherever `f32` holds the score: to the rounding of
/// their different renormalisation points and of the final sum, which the
/// strip kernel takes in `f32` and the diagonal kernel in `f64`. Whether a
/// path exists at all does not depend on the traversal.
#[hegel::test(test_cases = pinned::cases(1024))]
fn strips_and_diagonals_are_the_same_recurrence(tc: TestCase) {
    let case = tc.draw(arbitrary_case());
    let width = tc.draw(gs::integers::<u32>().min_value(2).max_value(63));
    let conversion = tc.draw(any_conversion());
    let band = Band::new(width, case.offset).expect("w");
    let taps = TapsEmission::new(conversion, Betas::PerSite(&case.betas));
    let diagonals = align_banded_simd(&case.haplotype, &case.read, &taps, band).get();
    let strips = align_strips_simd(&case.haplotype, &case.read, &taps, band).get();
    assert_eq!(
        diagonals.is_finite(),
        strips.is_finite(),
        "diagonals {} strips {}",
        diagonals,
        strips
    );
    if diagonals.is_finite() {
        assert!(
            (diagonals - strips).abs() < 1e-4 * (1.0 + diagonals.abs()),
            "width {}: diagonals {} strips {}",
            width,
            diagonals,
            strips
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
        let strips = align_strips_simd(&shifted, &read, &StandardEmission::default(), band).get();
        assert!(
            (strips - expected).abs() < 1e-4,
            "strips, pad {pad}: {strips} against the expected {expected} (unpadded {here})"
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
    let band = Band::anchored(1000);
    for (name, scalar, simd) in [
        (
            "diagonals",
            align_banded(&haplotype, &read, &StandardEmission::default(), band).get(),
            align_banded_simd(&haplotype, &read, &StandardEmission::default(), band).get(),
        ),
        (
            "strips",
            align_strips(&haplotype, &read, &StandardEmission::default(), band).get(),
            align_strips_simd(&haplotype, &read, &StandardEmission::default(), band).get(),
        ),
    ] {
        assert!(full.is_finite() && scalar.is_finite(), "{name}: full {full}, banded {scalar}");
        assert_eq!(scalar.to_bits(), simd.to_bits(), "{name}: scalar {scalar}, simd {simd}");
        assert!((scalar - full).abs() < 1e-2, "{name}: {scalar} against the reference {full}");
    }
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
        (
            "strips, no haplotype",
            align_strips(&empty_haplotype, &read, &StandardEmission::default(), Band::anchored(0)),
        ),
        (
            "strips, no read",
            align_strips_simd(
                &haplotype,
                &empty_read,
                &StandardEmission::default(),
                Band::anchored(0),
            ),
        ),
    ] {
        assert_eq!(got, compair::Log10Likelihood::IMPOSSIBLE, "{name}: {got:?}");
    }

    // Extreme offsets are arithmetic on `i64`, so neither end wraps.
    for offset in [i32::MIN, i32::MIN + 1, -1, 0, 1, i32::MAX - 1, i32::MAX] {
        let band = Band::anchored(offset);
        let got = align_banded(&haplotype, &read, &StandardEmission::default(), band).get();
        assert!(!got.is_nan(), "offset {offset} gave {got}");
        let strips = align_strips_simd(&haplotype, &read, &StandardEmission::default(), band).get();
        assert!(!strips.is_nan(), "strips at offset {offset} gave {strips}");
        assert_eq!(got.is_finite(), strips.is_finite(), "offset {offset}: {got} vs {strips}");
        if got.is_finite() {
            assert!(
                (got - strips).abs() < 1e-4,
                "offset {offset}: diagonals {got} strips {strips}"
            );
        }
    }

    // A gap-continuation quality of zero is a gap that never closes, and the
    // read's free start is a gap: no path reaches row one.
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
/// zero in `f32`. On these pairs, which score near `-75`, the `f32` diagonal
/// kernel on its own lost several log10 at width 320.
///
/// That is below the floor where every entry point rescores a pair in `f64`,
/// so what a caller gets is the recurrence at both widths.
#[test]
fn a_band_wider_than_the_read_is_still_the_recurrence() {
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
        worst_unbanded < 1e-5,
        "a band wider than the read is rescored where f32 cannot hold it, worst {worst_unbanded}"
    );
}

/// The same pairs through the strip entry point. It renormalises every eight
/// *rows* rather than per anti-diagonal, so a wide band costs it less range
/// than the diagonal kernel -- but not none: a row's cells are prefix
/// alignments ending at every column of the band, and a read much longer
/// than its haplotype spreads them past `f32` too. These pairs score below
/// the floor, so what this holds is the entry point: the recurrence at both
/// widths, deep into the dynamic range.
#[test]
fn the_strip_kernel_tracks_the_recurrence_at_any_width() {
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

    let mut worst = 0.0f64;
    let mut worst_score = 0.0f64;
    for _ in 0..200 {
        let haplotype = Haplotype::new(draw(&mut next, 90));
        let read = Read::uniform(
            draw(&mut next, 66),
            &[BaseQuality::from_byte(30); 66],
            BaseQuality::from_byte(45),
            BaseQuality::from_byte(45),
            BaseQuality::from_byte(10),
            Strand::OT,
        )
        .expect("valid");
        for width in [Band::DEFAULT_WIDTH, 320] {
            let band = Band::new(width, 0).expect("width at least two");
            let got =
                align_strips_simd(&haplotype, &read, &StandardEmission::default(), band).get();
            let want = align_masked(&haplotype, &read, &StandardEmission::default(), band);
            assert_eq!(got.is_finite(), want.is_finite(), "width {width}: f32 {got} f64 {want}");
            if !want.is_finite() {
                continue;
            }
            assert!(got <= want + 1e-6, "width {width}: f32 {got} over the f64 recurrence {want}");
            worst = worst.max((got - want).abs());
            worst_score = worst_score.min(want);
        }
    }
    assert!(
        worst < 1e-5,
        "the strip kernel tracks the f64 recurrence at both widths, worst {worst}"
    );
    assert!(
        worst_score < -60.0,
        "the pairs reach deep into the dynamic range, worst {worst_score}"
    );
}

/// A band wide enough to hold the whole matrix is not a constraint, so it must
/// give back the unbanded reference. Run on GATK's own 164 x 101 vectors, whose
/// scores are a few log10 and so well inside `f32`'s range -- see
/// `a_band_wider_than_the_read_is_still_the_recurrence` for what happens
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
        let strips =
            align_strips_simd(&vector.haplotype, &vector.read, &StandardEmission::default(), band)
                .get();
        let full = align_full(&vector.haplotype, &vector.read, &StandardEmission::default()).get();
        assert!(
            (banded - full).abs() < 1e-5,
            "vector {index}: unconstrained band {banded} against the reference {full}"
        );
        assert!(
            (strips - full).abs() < 1e-5,
            "vector {index}: unconstrained strips {strips} against the reference {full}"
        );
        worst = worst.max((banded - full).abs()).max((strips - full).abs());
    }
    assert!(worst < 1e-5, "worst {worst}");
}

/// The regression the properties above found, kept as a deterministic case
/// because it needs a band of width two at a negative offset and the
/// generators only reach that shape every few thousand draws.
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
            let strips = align_strips_simd(&haplotype, &read, &StandardEmission::default(), band);
            assert_eq!(
                got,
                compair::Log10Likelihood::IMPOSSIBLE,
                "width {width} at offset {offset} holds no (0, j): read {read_length} got {got:?}"
            );
            assert_eq!(got, simd);
            assert_eq!(got, strips, "strips, width {width} at offset {offset}");
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
        let strips = align_strips_simd(
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
        assert_eq!(got, strips, "strips at offset {offset}");
    }
}

/// Found by `fuzz_pair_hmm`, when a read at Phred 0 meant `eps = 1`: a read
/// whose last base is `N`, so every match emission was exactly zero and the
/// last row was reached through the insertion state alone, against a
/// haplotype that is mostly `N`. With `eps` capped at 3/4 the emissions are a
/// quarter and the pair no longer collapses a diagonal; it stays as a
/// regression pin, and the rescue covers any pair where the lifts would not.
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
    assert!((full - -1.8201).abs() < 1e-3, "the reference scores this pair at -1.82, got {full}");
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
        // The strip kernel has no lift chain to overflow, but the pair is a
        // read of three bases against a haplotype of `N`, and worth a pin.
        let strips = align_strips(&haplotype, &read, &emission, band).get();
        assert_eq!(
            strips.to_bits(),
            align_strips_simd(&haplotype, &read, &emission, band).get().to_bits()
        );
        assert!(
            (strips - full).abs() < 1e-4,
            "width {width}: strips {strips} against the reference {full}"
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

/// The batched kernel's gate: scoring a haplotype in a batch of eight is
/// bit-identical to scoring it alone through the scalar strip kernel.
///
/// It is a different traversal -- row-wise, one haplotype per lane -- so
/// this is not parity by construction the way the two strip lanes are; it
/// asserts that the row-wise band, the free start, the flush to zero and
/// the renormalisation cadence all land on the strip kernel's numbers.
/// The one-lane instance is carried as well, because it and the eight-lane
/// one *are* parity by construction and a divergence there would be a lane
/// bug rather than a traversal bug.
#[hegel::test(test_cases = pinned::cases(1024))]
fn a_batch_is_bit_identical_to_one_alignment_at_a_time(tc: TestCase) {
    let case = tc.draw(arbitrary_case());
    let conversion = tc.draw(any_conversion());
    let uniform = tc.draw(any_probability());
    let spread = tc.draw(gs::integers::<usize>().max_value(3));
    let band = case.band();
    let taps = TapsEmission::new(conversion, Betas::PerSite(&case.betas));
    let uniform = TapsEmission::new(conversion, Betas::Uniform(uniform));
    // A ragged batch: each next haplotype is the same one with a few bases
    // cut off the end, so the lanes do not share a length.
    let bases = case.haplotype.bases().to_vec();
    let batch: Vec<Haplotype> = (0..BATCH)
        .map(|k| Haplotype::new(bases[..bases.len().saturating_sub(k * spread)].to_vec()))
        .collect();
    let refs: Vec<&Haplotype> = batch.iter().collect();
    let mut workspace = Workspace::new();
    for (name, emission) in [
        ("standard", &StandardEmission::default() as &dyn ErasedEmission),
        ("taps", &taps),
        ("uniform", &uniform),
    ] {
        let one_at_a_time: Vec<u64> =
            refs.iter().map(|h| emission.strips(h, &case.read, band).get().to_bits()).collect();
        let wide = emission.batch(&mut workspace, &refs, &case.read, band);
        let scalar = emission.batch_scalar(&mut workspace, &refs, &case.read, band);
        assert_eq!(
            &one_at_a_time,
            &wide.iter().map(|s| s.get().to_bits()).collect::<Vec<_>>(),
            "{}: strips vs eight-lane batch",
            name
        );
        assert_eq!(
            &one_at_a_time,
            &scalar.iter().map(|s| s.get().to_bits()).collect::<Vec<_>>(),
            "{}: strips vs one-lane batch",
            name
        );
    }
}

/// The three emissions the parity test sweeps, behind one object so the test
/// body is written once. `Emission` is generic, not a trait object, because
/// the kernels monomorphise on it; this is the test's own indirection.
trait ErasedEmission {
    fn strips(&self, haplotype: &Haplotype, read: &Read, band: Band) -> compair::Log10Likelihood;
    fn batch(
        &self,
        workspace: &mut Workspace,
        haplotypes: &[&Haplotype],
        read: &Read,
        band: Band,
    ) -> Vec<compair::Log10Likelihood>;
    fn batch_scalar(
        &self,
        workspace: &mut Workspace,
        haplotypes: &[&Haplotype],
        read: &Read,
        band: Band,
    ) -> Vec<compair::Log10Likelihood>;
}

impl<E: compair::Emission> ErasedEmission for E {
    fn strips(&self, haplotype: &Haplotype, read: &Read, band: Band) -> compair::Log10Likelihood {
        align_strips(haplotype, read, self, band)
    }
    fn batch(
        &self,
        workspace: &mut Workspace,
        haplotypes: &[&Haplotype],
        read: &Read,
        band: Band,
    ) -> Vec<compair::Log10Likelihood> {
        let mut out = Vec::new();
        workspace.align_batch(haplotypes, read, self, band, &mut out);
        out
    }
    fn batch_scalar(
        &self,
        workspace: &mut Workspace,
        haplotypes: &[&Haplotype],
        read: &Read,
        band: Band,
    ) -> Vec<compair::Log10Likelihood> {
        let mut out = Vec::new();
        workspace.align_batch_scalar(haplotypes, read, self, band, &mut out);
        out
    }
}

/// The bench's own shape -- a 150 bp read against eight 200 bp candidates --
/// through the batch, against one strip alignment each. The property test's
/// reads are under 70 bases, i.e. nine renormalisations; this one runs
/// nineteen and is the size the shadow path actually scores.
#[test]
fn the_bench_shape_batches_to_the_same_bits() {
    let bases: Vec<Base> = (0..200u32)
        .scan(0x9e37_79b9_7f4a_7c15u64, |state, _| {
            *state ^= *state << 13;
            *state ^= *state >> 7;
            *state ^= *state << 17;
            Some(match *state % 4 {
                0 => Base::A,
                1 => Base::C,
                2 => Base::G,
                _ => Base::T,
            })
        })
        .collect();
    let read_bases = bases[25..175].to_vec();
    let quals: Vec<BaseQuality> =
        (0..read_bases.len()).map(|i| BaseQuality::from_byte(25 + (i % 15) as u8)).collect();
    let read = Read::uniform(
        read_bases,
        &quals,
        BaseQuality::from_byte(45),
        BaseQuality::from_byte(45),
        BaseQuality::from_byte(10),
        Strand::OT,
    )
    .expect("the fixture read must build");
    let betas: Vec<Probability> = (0..200)
        .map(|i| Probability::new(f64::from(i % 11) / 10.0).unwrap_or(Probability::ZERO))
        .collect();
    let taps = TapsEmission::new(ConversionModel::taps_default(), Betas::PerSite(&betas));
    let band = Band::anchored(25);
    let batch: Vec<Haplotype> = (0..BATCH)
        .map(|variant| {
            let mut alternative = bases.clone();
            for edit in 0..variant {
                if let Some(base) = alternative.get_mut(35 + edit * 17) {
                    *base = base.inverse();
                }
            }
            Haplotype::new(alternative)
        })
        .collect();
    let refs: Vec<&Haplotype> = batch.iter().collect();
    let one_at_a_time: Vec<u64> =
        refs.iter().map(|h| align_strips(h, &read, &taps, band).get().to_bits()).collect();
    let batched: Vec<u64> =
        compair::align_batch(&refs, &read, &taps, band).iter().map(|s| s.get().to_bits()).collect();
    assert_eq!(one_at_a_time, batched);
    assert!(one_at_a_time.iter().all(|bits| f64::from_bits(*bits).is_finite()));
}

/// `align_candidates` picks a kernel per batch, so at nine haplotypes it uses
/// both of them for one call. Whichever it picks, the scores have to be the
/// ones the strip kernel gives -- that is what makes the dispatch a pure
/// performance decision rather than a second set of semantics.
///
/// Every group size from nothing to three full batches is checked, because the
/// interesting sizes are the boundaries: `BATCH_BREAK_EVEN - 1` and
/// `BATCH_BREAK_EVEN` take different paths, and `BATCH + 1` is the case that
/// was slower than `BATCH` before the dispatch existed.
#[test]
fn the_dispatch_never_changes_a_score() {
    let bases: Vec<Base> = (0..200u32)
        .scan(0x9e37_79b9_7f4a_7c15u64, |state, _| {
            *state ^= *state << 13;
            *state ^= *state >> 7;
            *state ^= *state << 17;
            Some(match *state % 4 {
                0 => Base::A,
                1 => Base::C,
                2 => Base::G,
                _ => Base::T,
            })
        })
        .collect();
    let read_bases = bases.get(25..175).map(<[Base]>::to_vec).expect("200 bases were built");
    let quals: Vec<BaseQuality> = (0..read_bases.len())
        .map(|index| BaseQuality::from_byte(25 + u8::try_from(index % 15).unwrap_or(0)))
        .collect();
    let read = Read::uniform(
        read_bases,
        &quals,
        BaseQuality::from_byte(45),
        BaseQuality::from_byte(45),
        BaseQuality::from_byte(10),
        Strand::OT,
    )
    .expect("the fixture read must build");
    let betas: Vec<Probability> = (0..200)
        .map(|index| {
            Probability::new(f64::from(u32::try_from(index % 11).unwrap_or(0)) / 10.0)
                .unwrap_or(Probability::ZERO)
        })
        .collect();
    let taps = TapsEmission::new(ConversionModel::taps_default(), Betas::PerSite(&betas));
    let band = Band::anchored(25);

    let haplotypes: Vec<Haplotype> = (0..3 * BATCH)
        .map(|variant| {
            let mut alternative = bases.clone();
            for edit in 0..variant {
                if let Some(base) = alternative.get_mut(35 + edit * 5) {
                    *base = base.inverse();
                }
            }
            Haplotype::new(alternative)
        })
        .collect();

    let mut workspace = Workspace::new();
    let mut out = Vec::new();
    for count in 0..=3 * BATCH {
        let refs: Vec<&Haplotype> = haplotypes.iter().take(count).collect();
        let one_at_a_time: Vec<u64> = refs
            .iter()
            .map(|haplotype| align_strips_simd(haplotype, &read, &taps, band).get().to_bits())
            .collect();

        workspace.align_candidates(&refs, &read, &taps, band, &mut out);
        let dispatched: Vec<u64> = out.iter().map(|score| score.get().to_bits()).collect();
        assert_eq!(dispatched, one_at_a_time, "{count} haplotypes through the workspace");

        let free: Vec<u64> = align_candidates(&refs, &read, &taps, band)
            .iter()
            .map(|score| score.get().to_bits())
            .collect();
        assert_eq!(free, one_at_a_time, "{count} haplotypes through the free function");
    }
}

/// Both entry points replace what is in `out` rather than appending to it.
/// Reusing the buffer is the reason it is passed in, and a version that
/// appended silently returned a previous call's scores.
#[test]
fn a_reused_output_buffer_holds_only_this_call() {
    let haplotype = Haplotype::from_ascii(b"ACGTACGTACGTACGTACGT");
    let read = Read::uniform(
        vec![Base::A, Base::C, Base::G, Base::T],
        &[BaseQuality::from_byte(30); 4],
        BaseQuality::from_byte(45),
        BaseQuality::from_byte(45),
        BaseQuality::from_byte(10),
        Strand::OT,
    )
    .expect("a four base read must build");
    let band = Band::anchored(0);
    let refs = [&haplotype, &haplotype, &haplotype];

    let mut workspace = Workspace::new();
    let mut out = Vec::new();
    for _ in 0..3 {
        workspace.align_candidates(&refs, &read, &StandardEmission::default(), band, &mut out);
        assert_eq!(out.len(), refs.len(), "align_candidates must not accumulate");
    }
    for _ in 0..3 {
        workspace.align_batch(&refs, &read, &StandardEmission::default(), band, &mut out);
        assert_eq!(out.len(), refs.len(), "align_batch must not accumulate");
    }
}
