//! Power-of-two rescaling.
//!
//! Every rescaling in this crate is a power of two so that it is exact in
//! binary floating point, which is what lets the scalar and SIMD kernels stay
//! bit-identical while their intermediate values are renormalised.
//!
//! The `f32` kernels shift a row (or anti-diagonal) so that its largest cell
//! sits at `2^STRIP_SCALE`, by at most 127, the largest power [`exp2_f32`]
//! represents. The banded `f64` recurrence shifts each row to
//! `2^F64_ROW_SCALE`, by at most 1023; `align_full` shifts into `[1, 2)`.

/// `2^exponent` exactly, or the nearest representable end of the range.
#[allow(
    clippy::cast_sign_loss,
    reason = "the branches above leave the biased exponent in 1..=2046"
)]
pub(crate) fn exp2_f64(exponent: i32) -> f64 {
    if exponent > 1023 {
        f64::INFINITY
    } else if exponent < -1022 {
        0.0
    } else {
        let bits = (u64::from(exponent.wrapping_add(1023) as u32)) << 52;
        f64::from_bits(bits)
    }
}

/// `2^exponent` exactly, or the nearest representable end of the range.
#[allow(clippy::cast_sign_loss, reason = "the branches above leave the biased exponent in 1..=254")]
pub(crate) fn exp2_f32(exponent: i32) -> f32 {
    if exponent > 127 {
        f32::INFINITY
    } else if exponent < -126 {
        0.0
    } else {
        let bits = (exponent.wrapping_add(127) as u32) << 23;
        f32::from_bits(bits)
    }
}

/// The shift that moves `max` into `[1, 2)`, or zero when there is nothing to
/// normalise: a non-positive, non-finite or subnormal maximum.
#[allow(clippy::cast_possible_wrap, reason = "the masked exponent field is 11 bits")]
pub(crate) fn normalising_shift_f64(max: f64) -> i32 {
    if max <= 0.0 || !max.is_finite() {
        return 0;
    }
    let exponent = ((max.to_bits() >> 52) & 0x7ff) as i32 - 1023;
    if exponent == -1023 { 0 } else { -exponent }
}

/// Where the banded `f64` recurrence (`align_banded_f64` and
/// `align_banded_f64_rows`) keeps a rescaled row: its largest cell in
/// `[2^F64_ROW_SCALE, 2^(F64_ROW_SCALE + 1))`, the start row included.
///
/// A value computed below `2^-1022` is rounded to a multiple of `2^-1074`,
/// which loses up to `2^-1075` of it, so what a row can lose is `2^-1075`
/// over its scale: the higher the scale, the further below its row's
/// largest cell a cell keeps every bit (`2^(-1022 - F64_ROW_SCALE)`), and the
/// lower the scores the recurrence is exact for (see
/// `crate::align_banded_f64`). At `[1, 2)` a cell 308 decades below its row's
/// largest was already rounding, and a read that switches paths after that
/// many decades scored wrong.
///
/// The bound. Every transition and every prior is at most one. So from a
/// row whose cells are all below `2^(S + 1)`, a match cell of the next row
/// is at most `3 * 2^(S + 1)` (its prior times three neighbours times a
/// transition each), an insertion cell `2 * 2^(S + 1)`, and a deletion cell,
/// a sum along the row of match cells each times a transition, at most
/// `columns * 3 * 2^(S + 1)`, with `columns <= 1025`: every value, the sums a
/// cell is made of among them, below `6150 * 2^S < 2^(S + 13)`. That holds
/// whatever the read's `Growth`, since it is one step from a row whose
/// largest cell is known. The total sums the last row's rescaled match and
/// insertion cells, at most `2 * 1025 * 2^(S + 1) < 2^(S + 13)`.
/// `align_banded_f64` also makes one row from a row whose deletions it has
/// not yet rescaled for (and throws that row away when they asked for
/// another scale): those deletions are at most `1025 * 2^(S + 1)`, so that
/// row's values stay below `2^(S + 24)`.
///
/// So `S + 24 <= 1023` keeps every value finite. 960 leaves the start row's
/// shift, `S - log2(1 / h) <= 960 + 64`, within the 1023 one factor holds,
/// and lets any other row fall `2^-1023` below the row above in one step
/// before its shift is capped: a fall that takes a prior or transition near
/// zero at every cell, which no quality byte gives (the smallest is about
/// `2^-85`). A capped row is only lifted less, never lowered, and the next
/// rows lift it the rest of the way.
pub(crate) const F64_ROW_SCALE: i32 = 960;

/// The shift that moves `max` into `[2^F64_ROW_SCALE, 2^(F64_ROW_SCALE +
/// 1))`, at most 1023 so that [`exp2_f64`] represents it, and zero when there
/// is nothing to scale: a non-positive or non-finite maximum. A subnormal
/// maximum counts as `2^-1023`, which lifts it exactly, if short of the
/// scale.
#[allow(clippy::cast_possible_wrap, reason = "the masked exponent field is 11 bits")]
pub(crate) fn row_shift_f64(max: f64) -> i32 {
    if max <= 0.0 || !max.is_finite() {
        return 0;
    }
    let exponent = ((max.to_bits() >> 52) & 0x7ff) as i32 - 1023;
    (F64_ROW_SCALE - exponent).min(1023)
}

/// Where the strip kernels (and the batch, pairs and GPU kernels, which share
/// their sweep) keep a renormalised row: its largest cell in
/// `[2^STRIP_SCALE, 2^(STRIP_SCALE + 1))`.
///
/// A flush loses less than `2^-126` of the scale, so the higher the scale the
/// less a flush can lose (see `reference::trusted`). 115 is the highest scale
/// at which no read with steady qualities can overflow `f32` in any band.
///
/// The bound. Let `window` and `spill` be the read's `Growth`, and row `R` a
/// renormalised row: every cell at most `2^(S + 1)`, in `columns` columns.
/// Split each path into a cell of rows `R + 1..=R + 8` at the last cell of
/// row `R` it visits. From there it takes a sequence of steps that, with the
/// cell it ends on, fixes which cell of row `R` it left: so the paths from
/// all of row `R` into one cell weigh no more than every step sequence from
/// each of the three kinds of cell, at most `window` apiece. Hence
///
/// - a match or insertion cell, and the sum a match cell takes its prior
///   times, is at most `3 * window * 2^(S + 1)`; a deletion cell at most
///   `spill` times that, and an insertion plus a deletion `1 + spill` times;
/// - the total, over the last row's match and insertion cells, of which a
///   path reaches one, is at most every cell of row `R` times `window`: with
///   a deletion at most `min(1, spill)` times its row's largest match, that is
///   `(2 + min(1, spill)) * columns * window * 2^(S + 1)`.
///
/// Rounding adds a factor below 1.001 inside a window. At `S = 115` with
/// steady qualities (`window` one) and the widest band (1025 columns), that is
/// at most `3 * 1025 * 1.001 * 2^116 < 2^127.6` and `3 * (1 + spill) * 2^116`
/// with `spill <= 257`: below `f32::MAX`. At 116 the total of a band that
/// wide could reach `2^128`. A read whose qualities step enough to break the
/// bound is scored in `f64` (`reference::trusted`); the check is `strips_fit`.
pub(crate) const STRIP_SCALE: i32 = 115;

/// The shift that moves a normal `max` into `[2^STRIP_SCALE, 2^(STRIP_SCALE +
/// 1))`, at most 127 so that [`exp2_f32`] represents it, and zero when there
/// is nothing to normalise: a non-positive, non-finite or subnormal maximum.
///
/// Capping the shift only ever leaves a row's exponent lower than asked,
/// never lower than it was (a shift up adds to a scale that was already
/// there), which is all the precision argument (`reference::trusted`) needs.
#[allow(clippy::cast_possible_wrap, reason = "the masked exponent field is 8 bits")]
pub(crate) fn strip_shift_f32(max: f32) -> i32 {
    if !max.is_normal() || max < 0.0 {
        return 0;
    }
    let exponent = ((max.to_bits() >> 23) & 0xff) as i32 - 127;
    (STRIP_SCALE - exponent).min(127)
}

#[cfg(test)]
mod tests {
    use super::{
        F64_ROW_SCALE, STRIP_SCALE, exp2_f32, exp2_f64, normalising_shift_f64, row_shift_f64,
        strip_shift_f32,
    };
    use hegel::TestCase;
    use hegel::generators as gs;

    /// `powi` is the independent oracle: it is a different computation
    /// that also happens to be exact for powers of two.
    #[hegel::test]
    fn exp2_f32_is_powi_over_the_whole_normal_range(tc: TestCase) {
        let exponent = tc.draw(gs::integers::<i32>().min_value(-126).max_value(127));
        assert_eq!(exp2_f32(exponent).to_bits(), 2f32.powi(exponent).to_bits());
    }

    #[hegel::test]
    fn exp2_f64_is_powi_over_the_whole_normal_range(tc: TestCase) {
        let exponent = tc.draw(gs::integers::<i32>().min_value(-1022).max_value(1023));
        assert_eq!(exp2_f64(exponent).to_bits(), 2f64.powi(exponent).to_bits());
    }

    /// Every normal maximum lands on the strip scale, `[2^S, 2^(S + 1))`, unless
    /// that takes more than the one factor `exp2_f32` has (a maximum below
    /// `2^(S - 127)`), where the shift stops at 127 and the maximum lands below
    /// it -- still a shift up, so a scale that was at least `2^S` stays so.
    /// Either way the shift is exact.
    #[hegel::test]
    fn a_normal_f32_maximum_is_moved_to_the_strip_scale(tc: TestCase) {
        let bits = tc.draw(
            gs::integers::<u32>()
                .min_value(f32::MIN_POSITIVE.to_bits())
                .max_value(f32::MAX.to_bits()),
        );
        let max = f32::from_bits(bits);
        let shift = strip_shift_f32(max);
        let scaled = max * exp2_f32(shift);
        let (low, high) = (exp2_f32(STRIP_SCALE), exp2_f32(STRIP_SCALE + 1));
        if shift == 127 {
            assert!(scaled < high, "{max} shifted by {shift} is {scaled}");
        } else {
            assert!((low..high).contains(&scaled), "{max} shifted by {shift} is {scaled}");
        }
        // Back in two factors: `2^-127` itself is not a normal `f32`.
        let back = scaled * exp2_f32(-shift / 2) * exp2_f32(-(shift - shift / 2));
        assert_eq!(back.to_bits(), max.to_bits());
    }

    /// A maximum with nothing to normalise asks for no shift.
    #[test]
    fn a_strip_maximum_that_is_not_normal_is_left_alone() {
        for max in [0.0, -1.0, f32::MIN_POSITIVE / 2.0, f32::INFINITY, f32::NAN] {
            assert_eq!(strip_shift_f32(max), 0, "{max}");
        }
    }

    #[hegel::test]
    fn a_normal_f64_maximum_is_moved_into_one_to_two(tc: TestCase) {
        let bits = tc.draw(
            gs::integers::<u64>()
                .min_value(f64::MIN_POSITIVE.to_bits())
                .max_value(f64::MAX.to_bits()),
        );
        let max = f64::from_bits(bits);
        let shift = normalising_shift_f64(max);
        let scaled = max * exp2_f64(shift);
        if shift == -1023 {
            assert!(max >= exp2_f64(1023));
            assert_eq!(scaled, 0.0);
        } else {
            assert!((1.0..2.0).contains(&scaled), "{max} shifted by {shift} is {scaled}");
            assert_eq!((scaled * exp2_f64(-shift)).to_bits(), max.to_bits());
        }
    }

    /// Every positive finite maximum is lifted exactly: onto the `f64` row
    /// scale, or where that takes more than 1023, by 1023 and below it. The
    /// largest possible value of a row (`2^(S + 24)`, see `F64_ROW_SCALE`)
    /// lands on the scale too.
    #[hegel::test]
    fn an_f64_maximum_is_moved_to_the_row_scale(tc: TestCase) {
        let bits = tc.draw(
            gs::integers::<u64>().min_value(1).max_value(exp2_f64(F64_ROW_SCALE + 24).to_bits()),
        );
        let max = f64::from_bits(bits);
        let shift = row_shift_f64(max);
        let scaled = max * exp2_f64(shift);
        assert_eq!((scaled * 2f64.powi(-shift)).to_bits(), max.to_bits(), "{max} by {shift}");
        let (low, high) = (exp2_f64(F64_ROW_SCALE), exp2_f64(F64_ROW_SCALE + 1));
        if shift == 1023 {
            assert!(scaled < high, "{max} shifted by {shift} is {scaled}");
        } else {
            assert!((low..high).contains(&scaled), "{max} shifted by {shift} is {scaled}");
        }
    }

    /// A subnormal maximum reports no shift: the kernel would carry it as
    /// is rather than scale it up.
    #[hegel::test]
    fn a_subnormal_f32_maximum_reports_no_shift(tc: TestCase) {
        let bits =
            tc.draw(gs::integers::<u32>().min_value(1).max_value(f32::MIN_POSITIVE.to_bits() - 1));
        let max = f32::from_bits(bits);
        assert!(max.is_subnormal());
        assert_eq!(strip_shift_f32(max), 0);
    }

    #[test]
    fn nothing_to_normalise_reports_no_shift() {
        for max in [0.0f32, -1.0, -0.0, f32::NEG_INFINITY, f32::INFINITY, f32::NAN] {
            assert_eq!(strip_shift_f32(max), 0, "{max}");
        }
        for max in [0.0f64, -1.0, -0.0, f64::NEG_INFINITY, f64::INFINITY, f64::NAN] {
            assert_eq!(normalising_shift_f64(max), 0, "{max}");
            assert_eq!(row_shift_f64(max), 0, "{max}");
        }
    }

    #[test]
    fn exp2_saturates_outside_the_normal_range() {
        assert_eq!(exp2_f32(128), f32::INFINITY);
        assert_eq!(exp2_f32(i32::MAX), f32::INFINITY);
        assert_eq!(exp2_f32(-127), 0.0);
        assert_eq!(exp2_f32(i32::MIN), 0.0);
        assert_eq!(exp2_f64(1024), f64::INFINITY);
        assert_eq!(exp2_f64(-1023), 0.0);
    }
}
