/// `2^exponent` exactly, or the nearest representable end of the range.
///
/// Every rescaling in this crate is a power of two so that it is exact in
/// binary floating point, which is what lets the scalar and SIMD kernels stay
/// bit-identical while their intermediate values are renormalised.
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
/// normalise.
#[allow(clippy::cast_possible_wrap, reason = "the masked exponent field is 11 bits")]
pub(crate) fn normalising_shift_f64(max: f64) -> i32 {
    if max <= 0.0 || !max.is_finite() {
        return 0;
    }
    let exponent = ((max.to_bits() >> 52) & 0x7ff) as i32 - 1023;
    if exponent == -1023 { 0 } else { -exponent }
}

#[allow(clippy::cast_possible_wrap, reason = "the masked exponent field is 8 bits")]
pub(crate) fn normalising_shift_f32(max: f32) -> i32 {
    if max <= 0.0 || !max.is_finite() {
        return 0;
    }
    let exponent = ((max.to_bits() >> 23) & 0xff) as i32 - 127;
    if exponent == -127 { 0 } else { -exponent }
}
