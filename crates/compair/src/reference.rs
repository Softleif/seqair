use seqair_types::Base;

use crate::{
    banded::Band,
    emission::{Emission, MatchProbability, epsilon},
    haplotype::Haplotype,
    read::Read,
    scaling::{STRIP_SCALE, exp2_f64, normalising_shift_f64},
    transitions::{Growth, Transition},
    types::Log10Likelihood,
};

/// The full `(read + 1) x (haplotype + 1)` dynamic program in `f64`, with no
/// band.
///
/// This is the crate's oracle: it reproduces GATK's `PairHMM` semantics, which
/// are worth stating because two of them are easy to get wrong. The read may
/// begin anywhere on the haplotype, encoded as `D[0][j] = 1 / haplotype_len`
/// for every `j`, so the free start is paid for with one `indel_to_match`
/// factor. And the total is the final row of the match and insertion matrices
/// only: the deletion matrix is excluded, so a read may not end inside a gap.
///
/// Rows are renormalised by a power of two after each step, which is exact, so
/// the result is the unscaled recurrence's to the last bit and cannot underflow
/// on a long read against a haplotype it does not fit.
#[allow(
    clippy::indexing_slicing,
    reason = "every index is in 0..=h and every row is allocated with h + 1 entries"
)]
pub fn align_full<E: Emission>(
    haplotype: &Haplotype,
    read: &Read,
    emission: &E,
) -> Log10Likelihood {
    let (h, r) = (haplotype.len(), read.len());
    if h == 0 || r == 0 {
        return Log10Likelihood::IMPOSSIBLE;
    }
    let mut prev_m = vec![0.0f64; h + 1];
    let mut prev_i = vec![0.0f64; h + 1];
    #[allow(
        clippy::cast_precision_loss,
        reason = "a haplotype long enough to lose precision here does not exist"
    )]
    let mut prev_d = vec![1.0 / h as f64; h + 1];
    let mut exponent = 0i32;
    rescale_row(&mut prev_m, &mut prev_i, &mut prev_d, &mut exponent);

    let (mut cur_m, mut cur_i, mut cur_d) =
        (vec![0.0f64; h + 1], vec![0.0f64; h + 1], vec![0.0f64; h + 1]);

    for i in 1..=r {
        let Some(observation) = read.observation(i - 1) else {
            return Log10Likelihood::IMPOSSIBLE;
        };
        let Some(t) = read.transition(i - 1) else {
            return Log10Likelihood::IMPOSSIBLE;
        };
        cur_m[0] = 0.0;
        cur_i[0] = 0.0;
        cur_d[0] = 0.0;
        for j in 1..=h {
            let Some(site) = haplotype.site(j - 1) else {
                return Log10Likelihood::IMPOSSIBLE;
            };
            let prior = emission.match_probability(site, observation);
            {
                cur_m[j] = prior
                    * (prev_m[j - 1] * t.match_to_match
                        + prev_i[j - 1] * t.indel_to_match
                        + prev_d[j - 1] * t.indel_to_match);
                cur_i[j] = prev_m[j] * t.match_to_insertion + prev_i[j] * t.gap_continuation;
                cur_d[j] = cur_m[j - 1] * t.match_to_deletion + cur_d[j - 1] * t.gap_continuation;
            }
        }
        rescale_row(&mut cur_m, &mut cur_i, &mut cur_d, &mut exponent);
        core::mem::swap(&mut prev_m, &mut cur_m);
        core::mem::swap(&mut prev_i, &mut cur_i);
        core::mem::swap(&mut prev_d, &mut cur_d);
    }

    let total: f64 = prev_m.iter().zip(prev_i.iter()).skip(1).map(|(m, i)| m + i).sum();
    if total <= 0.0 {
        return Log10Likelihood::IMPOSSIBLE;
    }
    Log10Likelihood::new(total.log10() - f64::from(exponent) * core::f64::consts::LOG10_2)
}

fn rescale_row(m: &mut [f64], i: &mut [f64], d: &mut [f64], exponent: &mut i32) {
    let max = m
        .iter()
        .chain(i.iter())
        .chain(d.iter())
        .fold(0.0f64, |acc, v| if *v > acc { *v } else { acc });
    let shift = normalising_shift_f64(max);
    if shift == 0 {
        return;
    }
    let factor = exp2_f64(shift);
    for value in m.iter_mut().chain(i.iter_mut()).chain(d.iter_mut()) {
        *value *= factor;
    }
    *exponent += shift;
}

/// The `f64` recurrence over a band: [`align_full`] restricted to the cells
/// `|j - i - offset| <= width / 2`, which is what every banded kernel
/// computes in `f32`.
///
/// It is what an `f32` score falls back to below the level where `f32` can
/// vouch for it (see [`trusted`]), and it is exposed for callers that want the
/// banded answer without the `f32` kernels' range. It visits only the band's
/// cells, one row at a time, and rescales each row by a power of two, so it
/// holds any score `f64` can write down -- but it is scalar, several times
/// slower than the kernels.
#[allow(
    clippy::indexing_slicing,
    reason = "every index is in 0..=h and every row is allocated with h + 1 entries"
)]
pub fn align_banded_f64<E: Emission + ?Sized>(
    haplotype: &Haplotype,
    read: &Read,
    emission: &E,
    band: Band,
) -> Log10Likelihood {
    let (h, r) = (haplotype.len(), read.len());
    if h == 0 || r == 0 {
        return Log10Likelihood::IMPOSSIBLE;
    }
    // The columns of row `i` inside the band and the haplotype, `1..=h` for
    // every row but the start row's `0..=h`.
    let span = |i: usize, first: usize| -> Option<(usize, usize)> {
        let (i, last) = (i64::try_from(i).ok()?, i64::try_from(h).ok()?);
        let centre = i.checked_add(band.offset)?;
        let low = centre.checked_sub(band.half_width)?.max(i64::try_from(first).ok()?);
        let high = centre.checked_add(band.half_width)?.min(last);
        (low <= high).then_some((usize::try_from(low).ok()?, usize::try_from(high).ok()?))
    };

    let mut prev_m = vec![0.0f64; h + 1];
    let mut prev_i = vec![0.0f64; h + 1];
    let mut prev_d = vec![0.0f64; h + 1];
    let Some(mut prev_span) = span(0, 0) else {
        return Log10Likelihood::IMPOSSIBLE;
    };
    #[allow(
        clippy::cast_precision_loss,
        reason = "a haplotype long enough to lose precision here does not exist"
    )]
    let init = 1.0 / h as f64;
    prev_d[prev_span.0..=prev_span.1].fill(init);
    let mut exponent = 0i32;
    let (mut cur_m, mut cur_i, mut cur_d) =
        (vec![0.0f64; h + 1], vec![0.0f64; h + 1], vec![0.0f64; h + 1]);
    // What the `cur` buffers still hold from two rows up.
    let mut stale: Option<(usize, usize)> = None;
    // The emission's two halves, each once: per haplotype column the band
    // reaches, and per read row. The composition per cell is what
    // `MatchProbability` does.
    let Some((first, last)) = band.columns(h, r) else {
        return Log10Likelihood::IMPOSSIBLE;
    };
    // Each column's latent weight for every base a read can show, one track
    // per base, so a row reads its priors' weights as one contiguous slice.
    let strand = read.strand();
    let mut latent: [Vec<f64>; 5] = Default::default();
    for index in first..=last {
        let Some(site) = haplotype.site(index) else {
            return Log10Likelihood::IMPOSSIBLE;
        };
        let weights = emission.site_weights(site, strand);
        for (track, base) in latent.iter_mut().zip(OBSERVABLE) {
            track.push(weights.latent_weight(base));
        }
    }

    for i in 1..=r {
        let (Some(observation), Some(t)) = (read.observation(i - 1), read.transition(i - 1)) else {
            return Log10Likelihood::IMPOSSIBLE;
        };
        let eps = epsilon(emission, observation);
        if let Some((low, high)) = stale {
            for buffer in [&mut cur_m, &mut cur_i, &mut cur_d] {
                buffer[low..=high].fill(0.0);
            }
        }
        // A row with no cell in the band ends every path.
        let Some((low, high)) = span(i, 1) else {
            return Log10Likelihood::IMPOSSIBLE;
        };
        let track = latent.get(observation.base.known_index().unwrap_or(UNKNOWN_TRACK));
        let Some(weights) = (low - 1)
            .checked_sub(first)
            .zip((high - 1).checked_sub(first))
            .and_then(|(from, to)| track?.get(from..=to))
        else {
            return Log10Likelihood::IMPOSSIBLE;
        };
        let max = banded_row(
            [&prev_m[low - 1..=high], &prev_i[low - 1..=high], &prev_d[low - 1..=high]],
            [&mut cur_m[low..=high], &mut cur_i[low..=high], &mut cur_d[low..=high]],
            weights,
            eps,
            &t,
        );
        let shift = normalising_shift_f64(max);
        if shift != 0 {
            let factor = exp2_f64(shift);
            for buffer in [&mut cur_m, &mut cur_i, &mut cur_d] {
                for value in &mut buffer[low..=high] {
                    *value *= factor;
                }
            }
            exponent += shift;
        }
        core::mem::swap(&mut prev_m, &mut cur_m);
        core::mem::swap(&mut prev_i, &mut cur_i);
        core::mem::swap(&mut prev_d, &mut cur_d);
        stale = Some(prev_span);
        prev_span = (low, high);
    }

    let total: f64 = prev_m.iter().zip(prev_i.iter()).skip(1).map(|(m, i)| m + i).sum();
    if total <= 0.0 {
        return Log10Likelihood::IMPOSSIBLE;
    }
    Log10Likelihood::new(total.log10() - f64::from(exponent) * core::f64::consts::LOG10_2)
}

/// The bases a read can show, in [`Base::known_index`] order and then
/// [`Base::Unknown`].
const OBSERVABLE: [Base; 5] = [Base::A, Base::C, Base::G, Base::T, Base::Unknown];
const UNKNOWN_TRACK: usize = 4;

/// One row of [`align_banded_f64`] over the band's columns `low..=high`, and
/// the largest cell it stored.
///
/// `previous` is the row above over `low - 1..=high`, `current` this row over
/// `low..=high`, whose cell left of `low` is zero (the buffer holds nothing
/// outside its row's span). `weights` are the columns' latent weights for this
/// row's base. The operations are [`align_full`]'s, in its order; the
/// deletion chain is carried in registers rather than read back from the row.
#[inline(always)]
fn banded_row(
    previous: [&[f64]; 3],
    current: [&mut [f64]; 3],
    weights: &[f64],
    eps: f64,
    t: &Transition,
) -> f64 {
    let [prev_m, prev_i, prev_d] = previous;
    let [cur_m, cur_i, cur_d] = current;
    let (mut left_m, mut left_d) = (0.0f64, 0.0f64);
    // One running maximum per matrix, so no cell waits on another's compare.
    let mut max = [0.0f64; 3];
    let above = prev_m.windows(2).zip(prev_i.windows(2)).zip(prev_d);
    let cells = cur_m.iter_mut().zip(cur_i.iter_mut()).zip(cur_d.iter_mut());
    for ((((m_above, i_above), d_diagonal), ((m_cell, i_cell), d_cell)), weight) in
        above.zip(cells).zip(weights)
    {
        let (&[m_diagonal, m_up], &[i_diagonal, i_up]) = (m_above, i_above) else {
            continue;
        };
        let prior = weight * (1.0 - eps) + (1.0 - weight) * (eps / 3.0);
        let m = prior
            * (m_diagonal * t.match_to_match
                + i_diagonal * t.indel_to_match
                + d_diagonal * t.indel_to_match);
        let i = m_up * t.match_to_insertion + i_up * t.gap_continuation;
        let d = left_m * t.match_to_deletion + left_d * t.gap_continuation;
        (*m_cell, *i_cell, *d_cell) = (m, i, d);
        for (max, cell) in max.iter_mut().zip([m, i, d]) {
            if cell > *max {
                *max = cell;
            }
        }
        (left_m, left_d) = (m, d);
    }
    max.into_iter().fold(0.0, |max, cell| if cell > max { cell } else { max })
}

/// An `f32` kernel's score where `f32` can vouch for it, and the `f64`
/// recurrence over the same band where it cannot.
///
/// The strip kernels (and the batch, pairs and GPU kernels, which share their
/// sweep) flush every stored cell below `2^-126` to zero, and scale every
/// eighth row by a power of two so that its largest cell sits in
/// `[2^STRIP_SCALE, 2^(STRIP_SCALE + 1))`. What the flushes can lose, with
/// `spill` and `carry` the read's `Growth` and `total = 10^log10_total`:
///
/// - *The scale.* The free start puts `1 / h` into each row-0 column that
///   feeds row 1, `rho` in all. So a match or insertion cell of the
///   recurrence is at most `rho * total`, any cell at most `max(1, spill)`
///   times that, and the power of two `2^e` a renormalised row is stored at
///   has `2^-e <= max(1, spill) * rho * total * 2^-STRIP_SCALE`. Where the
///   shift was capped at 127, `e` only grew, from 0 or above: `2^-e <=
///   2^-127` then.
/// - *A drop.* A stored cell below `2^-126`, so below `2^(-126 - e)` in the
///   recurrence's units. A GPU that flushes subnormal intermediates also
///   drops the two products a cell sums, each below `2^-126`: at most three
///   drops per cell and matrix, in `r` rows of at most `columns` cells.
/// - *Where it goes.* The paths from a match or insertion cell to the total
///   sum to at most `total`, from a deletion cell to at most `carry * total`.
/// - *Rounding.* A value the kernel computes is at most `(1 + 2^-24)^n`
///   times the exact recurrence, `n` the roundings on the longest chain of
///   operations behind it: at most seven per path step, the `f32` transition
///   and prior among them, over at most `2r + columns` steps inside the band,
///   and `columns + 40` more for the total's sum and the `f64` arithmetic of
///   these bounds. That enters twice, in the row maximum a scale is taken
///   from and in what a drop hands on.
///
/// So the total the kernel returns is short of the recurrence's by less than
/// `3 * r * columns * (2 + carry) * total * (1 + 2^-24)^(2n) *
/// max(max(1, spill) * rho * 2^(-126 - STRIP_SCALE), 2^(-126 - 127))`, and a
/// total a million times that is right to within a millionth. With steady
/// qualities `total` and `carry` are one and `spill` below it; for 150 bases
/// against 290 in a 64-wide band the floor is near `-56`.
///
/// Below it the kernel may have flushed anything up to the whole answer --
/// eight rows of confident mismatches take every cell under `2^-126` at once
/// -- and the pair is scored again by [`align_banded_f64`]. So is a pair the
/// strip kernel could overflow on (see `scaling::STRIP_SCALE`), which no read
/// with steady qualities is. Both checks depend on the pair alone, and every
/// kernel returns the same `f32` score for a pair, so the fallback cannot make
/// two kernels disagree. The factors take a pass over the read's rows, so a
/// bound from its quality bytes (`Read::growth_bound`) is tried first; it
/// vouches only where they would.
///
/// Every CPU entry point applies it already. The GPU kernel's scores come
/// back without their pairs' inputs, so a caller of [`crate::gpu`] passes each
/// through this itself.
#[must_use]
pub fn trusted<E: Emission + ?Sized>(
    score: Log10Likelihood,
    haplotype: &Haplotype,
    read: &Read,
    emission: &E,
    band: Band,
) -> Log10Likelihood {
    let h = haplotype.len();
    let vouches = |growth: Growth| {
        strips_fit(h, band, growth)
            && score.get() >= trust_floor(h, read.len(), band, STRIP_SCALE, growth.log10_total, growth)
    };
    // The bound vouches only where the exact factors would: every check is
    // monotone in them.
    if read.growth_bound().is_some_and(vouches) || vouches(read.growth()) {
        score
    } else {
        align_banded_f64(haplotype, read, emission, band)
    }
}

/// [`trusted`] for a score from the diagonal kernel, which scales each
/// anti-diagonal's largest cell into `[1, 2)`: the scale is `2^0`. It takes
/// that scale from a whole anti-diagonal, whose largest cell can sit in a
/// later row than a cell it drops, so both the scale's bound and the drop's
/// reach can be `Growth::total`, and the floor takes it twice.
pub(crate) fn trusted_diagonal<E: Emission + ?Sized>(
    score: Log10Likelihood,
    haplotype: &Haplotype,
    read: &Read,
    emission: &E,
    band: Band,
) -> Log10Likelihood {
    let h = haplotype.len();
    let vouches = |growth: Growth| {
        score.get() >= trust_floor(h, read.len(), band, 0, 2.0 * growth.log10_total, growth)
    };
    if read.growth_bound().is_some_and(vouches) || vouches(read.growth()) {
        score
    } else {
        align_banded_f64(haplotype, read, emission, band)
    }
}

/// `(1 + 2^-24)^n` for the longest chain of roundings inside one strip
/// window, `n <= 7 * (16 + 1025)`: eight rows and a band's width of
/// deletions, seven roundings a step.
const WINDOW_ROUNDING: f64 = 1.001;

/// Whether no value the strip kernels compute for this read in this band can
/// overflow `f32`: the bound `scaling::STRIP_SCALE` proves, with this read's
/// [`Growth`].
#[allow(clippy::cast_precision_loss, reason = "a count of cells, far below 2^52")]
fn strips_fit(h: usize, band: Band, growth: Growth) -> bool {
    let columns = row_columns(band, h.saturating_add(1)) as f64;
    let cell = 3.0 * (1.0 + growth.spill);
    let total = (2.0 + growth.spill.min(1.0)) * columns;
    let bound = cell.max(total) * growth.window * WINDOW_ROUNDING * exp2_f64(STRIP_SCALE + 1);
    bound < f64::from(f32::MAX)
}

/// `log10` of the least total the flushes of an `f32` kernel that keeps its
/// rows at `2^scale` cannot move by more than a millionth; see [`trusted`].
/// `spread` is `log10` of how far a dropped cell's paths can grow,
/// `Growth::log10_total` for the strip kernels.
#[allow(clippy::cast_precision_loss, reason = "counts of cells, far below 2^52")]
fn trust_floor(h: usize, r: usize, band: Band, scale: i32, spread: f64, growth: Growth) -> f64 {
    let columns = row_columns(band, h) as f64;
    let rows = r as f64;
    let rho = free_starts(band, h) as f64 / h as f64;
    let roundings = 7.0 * (2.0 * rows + columns) + columns + 40.0;
    let rounding = roundings * f64::from(f32::EPSILON) * core::f64::consts::LOG10_E;
    let scale = (growth.spill.max(1.0) * rho * exp2_f64(-126 - scale)).max(exp2_f64(-253));
    6.0 + rounding + spread + (3.0 * rows * columns * (2.0 + growth.carry) * scale).log10()
}

/// The cells a row holds inside the band and `limit` columns.
fn row_columns(band: Band, limit: usize) -> usize {
    let width = usize::try_from(band.half_width).unwrap_or(usize::MAX).saturating_mul(2);
    width.saturating_add(1).min(limit)
}

/// The row-0 columns inside the band that feed row 1, `0..h`: the cells the
/// free start's `1 / h` goes into.
#[allow(clippy::cast_sign_loss, reason = "clamped at zero")]
fn free_starts(band: Band, h: usize) -> usize {
    let last = i64::try_from(h).unwrap_or(i64::MAX).saturating_sub(1);
    let low = band.offset.saturating_sub(band.half_width).max(0);
    let high = band.offset.saturating_add(band.half_width).min(last);
    high.saturating_sub(low).saturating_add(1).max(0) as usize
}

#[cfg(test)]
mod tests {
    use hegel::TestCase;
    use hegel::generators as gs;
    use seqair_types::{Base, BaseQuality, Strand};

    use super::{align_banded_f64, free_starts};
    use crate::{Band, Haplotype, Read, StandardEmission};

    /// No total of the recurrence exceeds what the free start puts in times
    /// the read's `Growth::total`, the bound every floor rests on: the `f64`
    /// recurrence is the independent side. Each quality track steps between
    /// two values drawn from a few far apart, at Q0 and Q1 the deletion runs
    /// hand on most, and the bases come from one or two letters at Q60, so a
    /// read fits its haplotype and nothing but the transitions holds the
    /// total down.
    #[hegel::test(test_cases = crate::pinned::cases(2048))]
    fn no_total_exceeds_the_free_start_times_the_growth(tc: TestCase) {
        let r = tc.draw(gs::integers::<usize>().min_value(1).max_value(48));
        let h = tc.draw(gs::integers::<usize>().min_value(1).max_value(200));
        let alphabet: &[Base] =
            if tc.draw(gs::booleans()) { &[Base::A] } else { &[Base::A, Base::C] };
        let letters = |n: usize| -> Vec<Base> {
            (0..n).map(|_| tc.draw_silent(gs::sampled_from(alphabet))).collect()
        };
        let (haplotype, bases) = (Haplotype::new(letters(h)), letters(r));
        let steps = [0u8, 1, 2, 6, 10, 20, 43, 254];
        let quals = |n: usize| -> Vec<BaseQuality> {
            let pair = [tc.draw(gs::sampled_from(&steps)), tc.draw(gs::sampled_from(&steps))];
            (0..n).map(|_| BaseQuality::from_byte(tc.draw_silent(gs::sampled_from(&pair)))).collect()
        };
        let (insertion, deletion, gap) = (quals(r), quals(r), quals(r));
        let read = Read::new(
            bases,
            &vec![BaseQuality::from_byte(60); r],
            &insertion,
            &deletion,
            &gap,
            Strand::OT,
        )
        .expect("a valid read");
        let width = tc.draw(gs::integers::<u32>().min_value(2).max_value(400));
        let offset = tc.draw(gs::integers::<i32>().min_value(-8).max_value(200));
        let band = Band::new(width, offset).expect("a legal width");
        let score = align_banded_f64(&haplotype, &read, &StandardEmission::default(), band).get();
        #[allow(clippy::cast_precision_loss, reason = "test sizes")]
        let start = free_starts(band, h) as f64 / h as f64;
        let bound = start.log10() + read.growth().log10_total;
        assert!(score <= bound + 1e-9, "{score} above the bound {bound}");
    }

    /// Steady qualities hand on exactly what they receive: every factor is
    /// one but `spill`, the deletion run's `match_to_deletion / indel_to_match`,
    /// and the scan's bound says so too.
    #[test]
    fn steady_qualities_do_not_grow() {
        let q = BaseQuality::from_byte;
        let read = Read::uniform(vec![Base::A; 40], &[q(30); 40], q(20), q(20), q(10), Strand::OT)
            .expect("a valid read");
        for growth in [read.growth(), read.growth_bound().expect("a steady gap quality")] {
            assert_eq!((growth.log10_total, growth.window, growth.carry), (0.0, 1.0, 1.0));
            assert!((growth.spill - 0.01 / 0.9).abs() < 1e-12, "{}", growth.spill);
        }
    }

    /// Where the gap-continuation quality holds still, the bound from the
    /// quality scans is at least the exact factors, every one of them: the
    /// entry points skip the exact pass whenever the bound vouches.
    #[hegel::test(test_cases = crate::pinned::cases(2048))]
    fn the_steady_bound_is_at_least_the_growth(tc: TestCase) {
        let r = tc.draw(gs::integers::<usize>().min_value(1).max_value(200));
        let steps = [0u8, 1, 2, 6, 10, 20, 30, 43, 60, 254];
        let quals = |n: usize| -> Vec<BaseQuality> {
            let pair = [tc.draw(gs::sampled_from(&steps)), tc.draw(gs::sampled_from(&steps))];
            (0..n).map(|_| BaseQuality::from_byte(tc.draw_silent(gs::sampled_from(&pair)))).collect()
        };
        let (insertion, deletion) = (quals(r), quals(r));
        let gap = BaseQuality::from_byte(tc.draw(gs::sampled_from(&steps)));
        let read = Read::new(
            vec![Base::A; r],
            &vec![BaseQuality::from_byte(30); r],
            &insertion,
            &deletion,
            &vec![gap; r],
            Strand::OT,
        )
        .expect("a valid read");
        let (exact, bound) = (read.growth(), read.growth_bound().expect("a steady gap quality"));
        let at_least = |bound: f64, exact: f64| bound >= exact * (1.0 - 1e-12);
        assert!(bound.log10_total >= exact.log10_total - 1e-12, "{bound:?} against {exact:?}");
        assert!(at_least(bound.window, exact.window), "{bound:?} against {exact:?}");
        assert!(at_least(bound.spill, exact.spill), "{bound:?} against {exact:?}");
        assert!(at_least(bound.carry, exact.carry), "{bound:?} against {exact:?}");
    }

    /// A read with a gap-continuation quality that changes has no bound
    /// from the scans; the entry points take the exact factors.
    #[test]
    fn a_changing_gap_quality_has_no_scan_bound() {
        let q = BaseQuality::from_byte;
        let read = Read::new(
            vec![Base::A; 3],
            &[q(30); 3],
            &[q(20); 3],
            &[q(20); 3],
            &[q(10), q(10), q(2)],
            Strand::OT,
        )
        .expect("a valid read");
        assert!(read.growth_bound().is_none());
    }
}
