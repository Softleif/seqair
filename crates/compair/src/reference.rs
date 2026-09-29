use seqair_types::Base;

use crate::{
    banded::Band,
    emission::{Emission, MatchProbability, epsilon},
    haplotype::Haplotype,
    read::Read,
    scaling::{STRIP_SCALE, exp2_f64, normalising_shift_f64},
    transitions::Transition,
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
/// sweep) scale every eighth row by a power of two so that its largest cell
/// sits in `[2^STRIP_SCALE, 2^(STRIP_SCALE + 1))`, `2^96`, and flush every
/// stored cell below `2^-126` to zero. With the transitions out of every
/// state a distribution and every emission at most one, each cell is a
/// probability, so the cell a scale is taken from is at most one, every scale
/// is at least `2^96`, and a flushed cell held less than `2^(-125 - 96)` in
/// absolute terms. The paths through it can reach the final row with no more
/// than that, and there are at most three cells per row per band column to
/// flush. So the total an `f32` kernel returns is short of the recurrence's
/// by less than `3 * (r + 1) * columns * 2^(-125 - 96)`, and a total a million
/// times that is right to within a millionth. (The diagonal kernel scales to
/// `[1, 2)` and checks against the same bound without the `2^-96`.)
///
/// Below that the kernel may have flushed anything up to the whole answer --
/// eight rows of confident mismatches take every cell under `2^-126` at once
/// -- and the pair is scored again by [`align_banded_f64`]. Every kernel
/// returns the same `f32` score for a pair and checks it against the same
/// floor, so the fallback cannot make two kernels disagree.
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
    trusted_at(STRIP_SCALE, score, haplotype, read, emission, band)
}

/// [`trusted`] for a score from the diagonal kernel, which scales each
/// anti-diagonal's largest cell into `[1, 2)` rather than to
/// `2^STRIP_SCALE`, so that a flush there lost less than `2^-125`.
pub(crate) fn trusted_diagonal<E: Emission + ?Sized>(
    score: Log10Likelihood,
    haplotype: &Haplotype,
    read: &Read,
    emission: &E,
    band: Band,
) -> Log10Likelihood {
    trusted_at(0, score, haplotype, read, emission, band)
}

/// [`trusted`] for a kernel that renormalises its rows to `2^scale`.
fn trusted_at<E: Emission + ?Sized>(
    scale: i32,
    score: Log10Likelihood,
    haplotype: &Haplotype,
    read: &Read,
    emission: &E,
    band: Band,
) -> Log10Likelihood {
    if score.get() >= trust_floor(haplotype.len(), read.len(), band, scale) {
        score
    } else {
        align_banded_f64(haplotype, read, emission, band)
    }
}

/// `log10` of the least total the flushes of an `f32` kernel that keeps its
/// rows at `2^scale` cannot move by more than a millionth:
/// `log10(3 * (r + 1) * columns * 2^(-125 - scale)) + 6`.
#[allow(clippy::cast_precision_loss, reason = "a count of cells, far below 2^52")]
fn trust_floor(h: usize, r: usize, band: Band, scale: i32) -> f64 {
    let band_columns = usize::try_from(band.half_width).unwrap_or(usize::MAX).saturating_mul(2);
    let columns = band_columns.saturating_add(1).min(h.saturating_add(1));
    let cells = 3.0 * (r as f64 + 1.0) * columns as f64;
    cells.log10() - (125.0 + f64::from(scale)) * core::f64::consts::LOG10_2 + 6.0
}
