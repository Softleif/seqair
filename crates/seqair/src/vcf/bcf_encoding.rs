//! Shared BCF typed value encoding constants and helpers.
//! Used by [`super::encoder`] and the unified writer in [`super::unified`].

use super::error::{ReservedMarker, VcfEncodeError};

// BCF type codes
pub const BCF_BT_NULL: u8 = 0;
pub const BCF_BT_INT8: u8 = 1;
pub const BCF_BT_INT16: u8 = 2;
pub const BCF_BT_INT32: u8 = 3;
pub const BCF_BT_FLOAT: u8 = 5;
pub const BCF_BT_CHAR: u8 = 7;

// Missing sentinels — r[bcf_writer.missing_sentinels]
pub const INT8_MISSING: u8 = 0x80;
pub const INT16_MISSING: u16 = 0x8000;
pub const INT32_MISSING: u32 = 0x80000000;
pub const FLOAT_MISSING: u32 = 0x7F800001;

// End-of-vector sentinels — r[bcf_writer.end_of_vector]
pub const INT8_END_OF_VECTOR: u8 = 0x81;
pub const INT16_END_OF_VECTOR: u16 = 0x8001;
pub const INT32_END_OF_VECTOR: u32 = 0x80000001;
pub const FLOAT_END_OF_VECTOR: u32 = 0x7F800002;

// Int ranges (sentinel-safe) — r[bcf_writer.smallest_int_type]
// The lower bounds are exactly `i8::MIN + 8` and `i16::MIN + 8`: the eight most
// negative values of each width are reserved, so this range excludes them.
pub const INT8_MIN: i32 = -120;
pub const INT8_MAX: i32 = 127;
pub const INT16_MIN: i32 = -32760;
pub const INT16_MAX: i32 = 32767;

// ── Reserved values ─────────────────────────────────────────────────────

// r[impl bcf_encoder.reserved_bands]
// r[impl bcf_encoder.reserved_width]
/// Last value of the int32 reserved band `0x80000000..=0x80000007`. It is the
/// only integer band a value can land in: `INT8_MIN` and `INT16_MIN` sit exactly
/// eight above their widths' minimums, so [`smallest_int_type`] promotes any
/// value in the int8 or int16 band to a wider type.
const INT32_RESERVED_LAST: i32 = i32::MIN + 7;
/// Number of reserved float bit patterns, `FLOAT_MISSING..=0x7F800007`.
const FLOAT_RESERVED_LEN: u32 = 7;

/// The band member at `offset` from the start of a reserved band.
fn marker_at(offset: u32) -> ReservedMarker {
    match offset {
        0 => ReservedMarker::Missing,
        1 => ReservedMarker::EndOfVector,
        _ => ReservedMarker::FutureUse,
    }
}

// r[impl bcf_encoder.reserved_bands]
/// The reserved band member `value` is encoded as, or `None` for ordinary data.
pub fn reserved_int(value: i32) -> Option<ReservedMarker> {
    (value <= INT32_RESERVED_LAST).then(|| marker_at(value.abs_diff(i32::MIN)))
}

// r[impl bcf_encoder.reserved_bands]
/// The reserved band member `value`'s bit pattern is, or `None` for ordinary
/// data. Compares exact bits, never `is_nan()`: quiet NaN is a valid Float.
pub fn reserved_float(value: f32) -> Option<ReservedMarker> {
    let offset = value.to_bits().wrapping_sub(FLOAT_MISSING);
    (offset < FLOAT_RESERVED_LEN).then(|| marker_at(offset))
}

// r[impl bcf_encoder.reserved_rejected]
// r[impl bcf_encoder.reserved_uniform]
/// Reject the first of `values` that BCF would read back as a reserved marker.
/// Called by the VCF text path too, so every output format accepts the same calls.
pub fn reject_reserved_ints(
    field: &str,
    values: impl IntoIterator<Item = i32>,
) -> Result<(), VcfEncodeError> {
    for value in values {
        if let Some(marker) = reserved_int(value) {
            return Err(VcfEncodeError::ReservedIntValue { field: field.into(), value, marker });
        }
    }
    Ok(())
}

// r[impl bcf_encoder.reserved_rejected]
// r[impl bcf_encoder.reserved_uniform]
/// The float twin of [`reject_reserved_ints`].
pub fn reject_reserved_floats(
    field: &str,
    values: impl IntoIterator<Item = f32>,
) -> Result<(), VcfEncodeError> {
    for value in values {
        if let Some(marker) = reserved_float(value) {
            return Err(VcfEncodeError::ReservedFloatValue {
                field: field.into(),
                bits: value.to_bits(),
                marker,
            });
        }
    }
    Ok(())
}

// r[impl bcf_writer.typed_values]
/// Write a BCF type byte: `(count << 4) | type_code`.
/// For count >= 15, emits overflow encoding with a typed integer count.
#[expect(
    clippy::cast_possible_truncation,
    reason = "each branch is guarded by a range check; else branch: count > u32::MAX requires >4 billion elements, exceeding addressable memory"
)]
pub fn encode_type_byte(buf: &mut Vec<u8>, count: usize, type_code: u8) {
    if count < 15 {
        buf.push(((count as u8) << 4) | type_code);
    } else {
        buf.push((15 << 4) | type_code);
        if count <= INT8_MAX as usize {
            buf.push((1 << 4) | BCF_BT_INT8);
            buf.push(count as u8);
        } else if count <= INT16_MAX as usize {
            buf.push((1 << 4) | BCF_BT_INT16);
            buf.extend_from_slice(&(count as u16).to_le_bytes());
        } else {
            buf.push((1 << 4) | BCF_BT_INT32);
            buf.extend_from_slice(&(count as u32).to_le_bytes());
        }
    }
}

// r[impl bcf_writer.string_encoding]
/// Encode a string as a BCF typed char vector (type code 7).
pub fn encode_typed_string(buf: &mut Vec<u8>, s: &[u8]) {
    encode_type_byte(buf, s.len(), BCF_BT_CHAR);
    buf.extend_from_slice(s);
}

/// Select the smallest BCF integer type that fits all values.
pub fn smallest_int_type(values: &[i32]) -> u8 {
    smallest_int_type_iter(values.iter().copied())
}

/// Like [`smallest_int_type`] but takes an iterator — avoids collecting into a Vec.
pub fn smallest_int_type_iter(values: impl Iterator<Item = i32>) -> u8 {
    let mut needs_16 = false;
    for v in values {
        if !(INT8_MIN..=INT8_MAX).contains(&v) {
            if !(INT16_MIN..=INT16_MAX).contains(&v) {
                return BCF_BT_INT32;
            }
            needs_16 = true;
        }
    }
    if needs_16 { BCF_BT_INT16 } else { BCF_BT_INT8 }
}

/// Encode a slice of i32 values as a BCF typed integer vector with smallest fitting type.
pub fn encode_typed_int_vec(buf: &mut Vec<u8>, values: &[i32]) {
    if values.is_empty() {
        encode_type_byte(buf, 0, BCF_BT_INT8);
        return;
    }
    let typ = smallest_int_type(values);
    encode_type_byte(buf, values.len(), typ);
    for &v in values {
        encode_int_as(buf, v, typ);
    }
}

/// Encode a dictionary index as a typed int (used for INFO/FORMAT/FILTER keys).
#[expect(
    clippy::cast_possible_truncation,
    reason = "each branch is guarded by a range check that ensures dict_idx fits in the target type"
)]
pub fn encode_typed_int_key(buf: &mut Vec<u8>, dict_idx: u32) {
    if dict_idx <= INT8_MAX as u32 {
        encode_type_byte(buf, 1, BCF_BT_INT8);
        buf.push(dict_idx as u8);
    } else if dict_idx <= INT16_MAX as u32 {
        encode_type_byte(buf, 1, BCF_BT_INT16);
        buf.extend_from_slice(&(dict_idx as u16).to_le_bytes());
    } else {
        encode_type_byte(buf, 1, BCF_BT_INT32);
        buf.extend_from_slice(&dict_idx.to_le_bytes());
    }
}

/// Write an i32 value using the specified BCF int type.
/// The caller must ensure `value` fits in the target type (via `smallest_int_type`).
#[expect(
    clippy::cast_possible_truncation,
    reason = "debug_assert above verifies value fits in the target type"
)]
pub fn encode_int_as(buf: &mut Vec<u8>, value: i32, typ: u8) {
    match typ {
        BCF_BT_INT8 => {
            debug_assert!((INT8_MIN..=INT8_MAX).contains(&value), "INT8 overflow: {value}");
            buf.push(value as u8);
        }
        BCF_BT_INT16 => {
            debug_assert!((INT16_MIN..=INT16_MAX).contains(&value), "INT16 overflow: {value}");
            buf.extend_from_slice(&(value as i16).to_le_bytes());
        }
        _ => buf.extend_from_slice(&value.to_le_bytes()),
    }
}

/// Write an integer value or the missing sentinel for the given type.
pub fn encode_int_value_or_missing(buf: &mut Vec<u8>, value: Option<i32>, typ: u8) {
    match value {
        Some(n) => encode_int_as(buf, n, typ),
        None => encode_int_missing(buf, typ),
    }
}

/// Write the missing sentinel for the given int type.
pub fn encode_int_missing(buf: &mut Vec<u8>, typ: u8) {
    match typ {
        BCF_BT_INT8 => buf.push(INT8_MISSING),
        BCF_BT_INT16 => buf.extend_from_slice(&INT16_MISSING.to_le_bytes()),
        _ => buf.extend_from_slice(&INT32_MISSING.to_le_bytes()),
    }
}

/// Write the end-of-vector sentinel for the given int type.
pub fn encode_int_eov(buf: &mut Vec<u8>, typ: u8) {
    match typ {
        BCF_BT_INT8 => buf.push(INT8_END_OF_VECTOR),
        BCF_BT_INT16 => buf.extend_from_slice(&INT16_END_OF_VECTOR.to_le_bytes()),
        _ => buf.extend_from_slice(&INT32_END_OF_VECTOR.to_le_bytes()),
    }
}

#[cfg(test)]
mod tests {
    #![allow(
        clippy::arithmetic_side_effects,
        clippy::indexing_slicing,
        reason = "test code with known small values"
    )]

    use super::*;
    use hegel::prelude::*;

    /// The band member of the bytes the encoder actually writes for `value`,
    /// read back as the unsigned pattern of the width it chose — the oracle
    /// `reserved_int` is checked against, sharing none of its reasoning about
    /// which widths are reachable.
    fn marker_of_emitted_bytes(value: i32) -> Option<ReservedMarker> {
        let typ = smallest_int_type(&[value]);
        let mut buf = Vec::new();
        encode_int_as(&mut buf, value, typ);
        let (bits, band_start) = match typ {
            BCF_BT_INT8 => (u32::from(buf[0]), u32::from(INT8_MISSING)),
            BCF_BT_INT16 => {
                (u32::from(u16::from_le_bytes([buf[0], buf[1]])), u32::from(INT16_MISSING))
            }
            _ => (u32::from_le_bytes([buf[0], buf[1], buf[2], buf[3]]), INT32_MISSING),
        };
        let offset = bits.wrapping_sub(band_start);
        (offset < 8).then(|| marker_at(offset))
    }

    // r[verify bcf_encoder.reserved_bands]
    // r[verify bcf_encoder.reserved_width]
    /// `reserved_int` looks at the int32 band only, on the argument that
    /// `smallest_int_type` never emits a value at a width whose band it lands
    /// in. The spec says that must be verified rather than assumed, so this
    /// checks it against the bytes actually written, with the draw weighted
    /// onto all three bands — a uniform `i32` reaches one about once in 2^29.
    #[hegel::test]
    fn reserved_int_matches_the_emitted_bytes(tc: TestCase) {
        let band = |min: i32| gs::integers::<i32>().min_value(min).max_value(min + 16);
        let value = tc.draw(hegel::one_of!(
            gs::integers::<i32>(),
            band(i32::MIN),
            band(i32::from(i16::MIN)),
            band(i32::from(i8::MIN)),
        ));
        assert_eq!(reserved_int(value), marker_of_emitted_bytes(value), "value {value}");
    }

    // r[verify bcf_encoder.reserved_bands]
    /// The bands start at the sentinels the encoder itself writes, and the
    /// float band has +inf just below it and quiet NaN above it.
    #[test]
    fn the_bands_start_at_the_sentinels_the_encoder_emits() {
        assert_eq!(reserved_int(i32::MIN), Some(ReservedMarker::Missing));
        assert_eq!(reserved_int(i32::MIN + 1), Some(ReservedMarker::EndOfVector));
        assert_eq!(reserved_int(i32::MIN + 7), Some(ReservedMarker::FutureUse));
        assert_eq!(reserved_int(i32::MIN + 8), None);
        assert_eq!(i32::MIN.to_le_bytes(), INT32_MISSING.to_le_bytes());

        assert_eq!(reserved_float(f32::from_bits(FLOAT_MISSING)), Some(ReservedMarker::Missing));
        assert_eq!(
            reserved_float(f32::from_bits(FLOAT_END_OF_VECTOR)),
            Some(ReservedMarker::EndOfVector)
        );
        assert_eq!(reserved_float(f32::from_bits(0x7F80_0007)), Some(ReservedMarker::FutureUse));
        assert_eq!(reserved_float(f32::from_bits(0x7F80_0008)), None);
        assert_eq!(reserved_float(f32::INFINITY), None);
        assert_eq!(reserved_float(f32::NAN), None);
        assert_eq!(reserved_float(f32::NEG_INFINITY), None);
    }
}
