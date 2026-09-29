//! Power-of-two rescaling.
//!
//! Every rescaling in this crate is a power of two so that it is exact in
//! binary floating point, which is what lets the scalar and SIMD kernels stay
//! bit-identical while their intermediate values are renormalised.
//!
//! The two halves are asymmetric on purpose. [`normalising_shift_f32`] reports
//! the shift that moves a normal `max` into `[1, 2)`, which is a value in
//! `-127..=126`; [`exp2_f32`] represents `-126..=127` exactly and saturates
//! outside it. So a maximum at or above `2^127` asks for a shift of `-127`,
//! which `exp2_f32` rounds to zero. No diagonal of the DP gets there -- every
//! cell is a product of probabilities times a value in `[1, 2)` -- and a
//! subnormal maximum reports no shift at all, so it is never scaled *up* into
//! the normal range. Both are pinned by the tests below rather than handled,
//! because the kernel cannot reach them and a branch there would be dead.

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

/// The shift that moves `max` into `[1, 2)`, or zero when there is nothing to
/// normalise: a non-positive, non-finite or subnormal maximum.
#[allow(clippy::cast_possible_wrap, reason = "the masked exponent field is 8 bits")]
pub(crate) fn normalising_shift_f32(max: f32) -> i32 {
    if max <= 0.0 || !max.is_finite() {
        return 0;
    }
    let exponent = ((max.to_bits() >> 23) & 0xff) as i32 - 127;
    if exponent == -127 { 0 } else { -exponent }
}

/// Where the strip kernels (and the batch and pairs kernels, which share their
/// sweep) keep a renormalised row: its largest cell in
/// `[2^STRIP_SCALE, 2^(STRIP_SCALE + 1))`.
///
/// A flush loses less than `2^-126` of the scale, so the higher the scale the
/// less a flush can lose: at `2^96` it is `2^-222` of a row's largest cell
/// rather than `2^-126`, and the `f32` kernels vouch for 29 more decades of
/// scores (see `reference::trusted_strips`). Nothing can overflow: every cell
/// of the eight rows below a renormalised row is a sum over paths from that
/// row's `3 * columns` cells, each continuing with probability at most one,
/// so no cell exceeds `3 * 1025 * 2^97 < 2^109`, and a total sums fewer
/// than `2^11` of those.
pub(crate) const STRIP_SCALE: i32 = 96;

/// The shift that moves a normal `max` into `[2^STRIP_SCALE, 2^(STRIP_SCALE +
/// 1))`, at most 127 so that [`exp2_f32`] represents it, and zero when there
/// is nothing to normalise: a non-positive, non-finite or subnormal maximum.
///
/// Capping the shift only ever leaves a row's scale above `2^STRIP_SCALE`
/// times its true maximum (a shift up adds to a scale that was already
/// there), which is all the precision argument needs.
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
        STRIP_SCALE, exp2_f32, exp2_f64, normalising_shift_f32, normalising_shift_f64,
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

    /// The contract the kernels rely on: after the shift the maximum is
    /// in `[1, 2)`, and the shift is a lossless operation.
    #[hegel::test]
    fn a_normal_f32_maximum_is_moved_into_one_to_two(tc: TestCase) {
        let bits = tc.draw(
            gs::integers::<u32>()
                .min_value(f32::MIN_POSITIVE.to_bits())
                .max_value(f32::MAX.to_bits()),
        );
        let max = f32::from_bits(bits);
        let shift = normalising_shift_f32(max);
        let scaled = max * exp2_f32(shift);
        if shift == -127 {
            // Above `2^127` the factor saturates to zero; see the module
            // docs for why this is pinned and not handled.
            assert!(max >= exp2_f32(127));
            assert_eq!(scaled, 0.0);
        } else {
            assert!((1.0..2.0).contains(&scaled), "{max} shifted by {shift} is {scaled}");
            assert_eq!((scaled * exp2_f32(-shift)).to_bits(), max.to_bits());
        }
    }

    /// Every normal maximum lands on the strip scale, `[2^96, 2^97)`, unless
    /// that takes more than the one factor `exp2_f32` has (a maximum below
    /// `2^-31`), where the shift stops at 127 and the maximum lands below
    /// it -- still a shift up, so a scale that was at least `2^96` stays so.
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

    /// A subnormal maximum reports no shift: the kernel would carry it as
    /// is rather than scale it up.
    #[hegel::test]
    fn a_subnormal_f32_maximum_reports_no_shift(tc: TestCase) {
        let bits =
            tc.draw(gs::integers::<u32>().min_value(1).max_value(f32::MIN_POSITIVE.to_bits() - 1));
        let max = f32::from_bits(bits);
        assert!(max.is_subnormal());
        assert_eq!(normalising_shift_f32(max), 0);
    }

    #[test]
    fn nothing_to_normalise_reports_no_shift() {
        for max in [0.0f32, -1.0, -0.0, f32::NEG_INFINITY, f32::INFINITY, f32::NAN] {
            assert_eq!(normalising_shift_f32(max), 0, "{max}");
        }
        for max in [0.0f64, -1.0, -0.0, f64::NEG_INFINITY, f64::INFINITY, f64::NAN] {
            assert_eq!(normalising_shift_f64(max), 0, "{max}");
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
