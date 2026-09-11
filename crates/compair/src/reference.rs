use crate::{
    emission::Emission,
    haplotype::Haplotype,
    read::Read,
    scaling::{exp2_f64, normalising_shift_f64},
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
