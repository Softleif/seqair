//! Shared BCF typed value encoding constants and helpers.
//! Used by both [`super::bcf_writer`] (`VcfRecord` path) and [`super::encoder`] (direct path).

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
pub const INT8_MIN: i32 = -120;
pub const INT8_MAX: i32 = 127;
pub const INT16_MIN: i32 = -32760;
pub const INT16_MAX: i32 = 32767;

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
