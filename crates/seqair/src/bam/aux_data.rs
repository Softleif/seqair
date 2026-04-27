//! Mutable auxiliary tag storage for owned BAM records.
//!
//! [`AuxData`] stores tags in raw BAM binary format and provides methods for
//! reading, adding, replacing, and removing tags. Used by [`OwnedBamRecord`](super::owned_record::OwnedBamRecord).

use super::aux::{self, AuxValue};
use thiserror::Error;

// r[impl bam.owned_record.aux_data]
/// Mutable auxiliary tag storage backed by raw BAM bytes.
///
/// Tags are stored in standard BAM binary format (2-byte name + 1-byte type + value).
/// Read operations delegate to [`aux::find_tag`]. Mutation operations modify the
/// underlying byte buffer directly.
#[derive(Debug, Clone, Default)]
pub struct AuxData {
    data: Vec<u8>,
}

/// Errors from auxiliary tag operations.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum AuxDataError {
    /// Integer value does not fit in any BAM integer type (i32/u32 range).
    #[error("integer value {value} out of BAM aux tag range (i32/u32)")]
    IntegerOutOfRange { value: i64 },

    /// Array data length is not a multiple of the element size.
    #[error("array data length {len} is not a multiple of element size {elem_size}")]
    InvalidArrayLength { len: usize, elem_size: usize },
}

impl AuxData {
    /// Create empty auxiliary data.
    pub fn new() -> Self {
        Self { data: Vec::new() }
    }

    // r[impl bam.owned_record.aux_from_slab]
    // r[impl base_mod.passthrough] — MM/ML are never decoded on the read path; the raw
    //   BAM aux byte block is moved into the owned record without inspection.
    /// Wrap existing raw BAM auxiliary bytes. No parsing or validation is performed.
    pub fn from_bytes(data: Vec<u8>) -> Self {
        Self { data }
    }

    /// Raw bytes for serialization.
    pub fn as_bytes(&self) -> &[u8] {
        &self.data
    }

    /// Whether there are no tags.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Look up a tag by its 2-byte name.
    pub fn get(&self, tag: [u8; 2]) -> Option<AuxValue<'_>> {
        aux::find_tag(&self.data, tag)
    }

    // r[impl bam.owned_record.aux_uniqueness]
    // r[impl bam.owned_record.aux_replace_semantics]
    // r[impl base_mod.passthrough.bam] — Z-type set/get round-trips MM bytes verbatim
    // r[impl base_mod.passthrough.coexistence] — `remove(tag)` only touches the named tag
    /// Add or replace a Z-type (null-terminated string) tag.
    pub fn set_string(&mut self, tag: [u8; 2], value: &[u8]) {
        self.remove(tag);
        self.data.extend_from_slice(&tag);
        self.data.push(b'Z');
        self.data.extend_from_slice(value);
        self.data.push(0); // NUL terminator
    }

    /// Add or replace a char (A-type) tag.
    pub fn set_char(&mut self, tag: [u8; 2], value: u8) {
        self.remove(tag);
        self.data.extend_from_slice(&tag);
        self.data.push(b'A');
        self.data.push(value);
    }

    // r[impl bam.owned_record.aux_int_encoding]
    /// Add or replace an integer tag, auto-selecting the smallest BAM type.
    ///
    /// Unsigned types are preferred for non-negative values (matching htslib behavior):
    /// - Non-negative: C (u8), S (u16), I (u32)
    /// - Negative: c (i8), s (i16), i (i32)
    pub fn set_int(&mut self, tag: [u8; 2], value: i64) -> Result<(), AuxDataError> {
        // Validate range before modifying data to avoid orphaned bytes on error
        if value >= 0 {
            if u64::try_from(value).map_or(true, |v| v > u64::from(u32::MAX)) {
                return Err(AuxDataError::IntegerOutOfRange { value });
            }
        } else if value < i64::from(i32::MIN) {
            return Err(AuxDataError::IntegerOutOfRange { value });
        }

        self.remove(tag);
        self.data.extend_from_slice(&tag);

        if value >= 0 {
            #[allow(clippy::cast_sign_loss, reason = "validated non-negative above")]
            let v = value as u64;
            #[allow(clippy::cast_possible_truncation, reason = "validated max size")]
            if v <= u64::from(u8::MAX) {
                self.data.push(b'C');
                self.data.push(v as u8);
            } else if v <= u64::from(u16::MAX) {
                self.data.push(b'S');
                self.data.extend_from_slice(&(v as u16).to_le_bytes());
            } else {
                self.data.push(b'I');
                self.data.extend_from_slice(&(v as u32).to_le_bytes());
            }
        } else {
            // Negative values — range already validated above

            #[allow(clippy::cast_possible_truncation, reason = "validated max size")]
            if value >= i64::from(i8::MIN) {
                self.data.push(b'c');
                self.data.push(value as u8);
            } else if value >= i64::from(i16::MIN) {
                self.data.push(b's');
                self.data.extend_from_slice(&(value as i16).to_le_bytes());
            } else {
                self.data.push(b'i');
                self.data.extend_from_slice(&(value as i32).to_le_bytes());
            }
        }

        Ok(())
    }

    /// Add or replace a float (f-type) tag.
    pub fn set_float(&mut self, tag: [u8; 2], value: f32) {
        self.remove(tag);
        self.data.extend_from_slice(&tag);
        self.data.push(b'f');
        self.data.extend_from_slice(&value.to_le_bytes());
    }

    /// Add or replace a double-precision float (d-type) tag.
    pub fn set_double(&mut self, tag: [u8; 2], value: f64) {
        self.remove(tag);
        self.data.extend_from_slice(&tag);
        self.data.push(b'd');
        self.data.extend_from_slice(&value.to_le_bytes());
    }

    // r[impl bam.owned_record.aux_array_encoding]
    // r[impl base_mod.passthrough.bam] — B:C round-trips ML bytes verbatim
    /// Add or replace a B:C (unsigned byte array) tag.
    ///
    /// The BAM B-array element count is u32. In practice, array size is bounded by
    /// the 2 MiB record size limit, so values longer than `u32::MAX` cannot occur in
    /// valid BAM data. We validate the cast to be safe on 64-bit platforms.
    pub fn set_array_u8(&mut self, tag: [u8; 2], values: &[u8]) -> Result<(), AuxDataError> {
        let count = array_count(values.len(), 1)?;
        self.remove(tag);
        encode_array_header(&mut self.data, &tag, b'C', count);
        self.data.extend_from_slice(values);
        Ok(())
    }

    // r[impl bam.owned_record.aux_array_setters]
    /// Add or replace a B:s (signed 16-bit array) tag.
    pub fn set_array_i16(&mut self, tag: [u8; 2], values: &[u8]) -> Result<(), AuxDataError> {
        let count = array_count(values.len(), 2)?;
        self.remove(tag);
        encode_array_header(&mut self.data, &tag, b's', count);
        self.data.extend_from_slice(values);
        Ok(())
    }

    /// Add or replace a B:S (unsigned 16-bit array) tag.
    pub fn set_array_u16(&mut self, tag: [u8; 2], values: &[u8]) -> Result<(), AuxDataError> {
        let count = array_count(values.len(), 2)?;
        self.remove(tag);
        encode_array_header(&mut self.data, &tag, b'S', count);
        self.data.extend_from_slice(values);
        Ok(())
    }

    /// Add or replace a B:i (signed 32-bit array) tag.
    pub fn set_array_i32(&mut self, tag: [u8; 2], values: &[u8]) -> Result<(), AuxDataError> {
        let count = array_count(values.len(), 4)?;
        self.remove(tag);
        encode_array_header(&mut self.data, &tag, b'i', count);
        self.data.extend_from_slice(values);
        Ok(())
    }

    /// Add or replace a B:I (unsigned 32-bit array) tag.
    pub fn set_array_u32(&mut self, tag: [u8; 2], values: &[u8]) -> Result<(), AuxDataError> {
        let count = array_count(values.len(), 4)?;
        self.remove(tag);
        encode_array_header(&mut self.data, &tag, b'I', count);
        self.data.extend_from_slice(values);
        Ok(())
    }

    /// Add or replace a B:f (float array) tag.
    pub fn set_array_f32(&mut self, tag: [u8; 2], values: &[u8]) -> Result<(), AuxDataError> {
        let count = array_count(values.len(), 4)?;
        self.remove(tag);
        encode_array_header(&mut self.data, &tag, b'f', count);
        self.data.extend_from_slice(values);
        Ok(())
    }

    /// Remove a tag if present.
    pub fn remove(&mut self, tag: [u8; 2]) {
        if let Some(range) = find_tag_byte_range(&self.data, tag) {
            self.data.drain(range);
        }
    }
}

#[allow(
    clippy::arithmetic_side_effects,
    reason = "checked by wrapping_rem above; divisor is the compile-time-known element size"
)]
/// Validate array element count and compute it from byte length.
///
/// Returns `InvalidArrayLength` if `byte_len` is not a multiple of `elem_size`,
/// or `IntegerOutOfRange` if the count exceeds `u32::MAX`.
fn array_count(byte_len: usize, elem_size: usize) -> Result<u32, AuxDataError> {
    if byte_len.wrapping_rem(elem_size) != 0 {
        return Err(AuxDataError::InvalidArrayLength { len: byte_len, elem_size });
    }
    let count = byte_len / elem_size;
    u32::try_from(count).map_err(|_| AuxDataError::IntegerOutOfRange {
        #[allow(clippy::cast_possible_wrap, reason = "count fits in i64 on 64-bit")]
        value: count as i64,
    })
}

/// Append the BAM B-array header (tag name + `B` type + subtype + 4-byte count)
/// to the given buffer.
fn encode_array_header(buf: &mut Vec<u8>, tag: &[u8; 2], subtype: u8, count: u32) {
    buf.extend_from_slice(tag);
    buf.push(b'B');
    buf.push(subtype);
    buf.extend_from_slice(&count.to_le_bytes());
}

/// Find the byte range `[start..end)` of a complete tag entry (name + type + value)
/// in raw BAM auxiliary data.
fn find_tag_byte_range(data: &[u8], target: [u8; 2]) -> Option<std::ops::Range<usize>> {
    let mut pos = 0usize;
    while pos.saturating_add(3) <= data.len() {
        let tag_start = pos;
        let tag = [*data.get(pos)?, *data.get(pos.saturating_add(1))?];
        let typ = *data.get(pos.saturating_add(2))?;
        pos = pos.checked_add(3)?;

        // Advance past the value
        pos = advance_past_value(data, pos, typ)?;

        if tag == target {
            return Some(tag_start..pos);
        }
    }
    None
}

/// Advance `pos` past one auxiliary tag value of the given type.
/// Returns the new position, or `None` if data is truncated/malformed.
fn advance_past_value(data: &[u8], mut pos: usize, typ: u8) -> Option<usize> {
    match typ {
        b'A' | b'c' | b'C' => pos = pos.checked_add(1)?,
        b's' | b'S' => pos = pos.checked_add(2)?,
        b'i' | b'I' | b'f' => pos = pos.checked_add(4)?,
        b'd' => pos = pos.checked_add(8)?,
        b'Z' | b'H' => {
            while pos < data.len() && *data.get(pos)? != 0 {
                pos = pos.checked_add(1)?;
            }
            pos = pos.checked_add(1)?; // skip NUL
        }
        b'B' => {
            let elem_type = *data.get(pos)?;
            let count_end = pos.checked_add(5)?;
            if count_end > data.len() {
                return None;
            }
            let count = u32::from_le_bytes([
                *data.get(pos.checked_add(1)?)?,
                *data.get(pos.checked_add(2)?)?,
                *data.get(pos.checked_add(3)?)?,
                *data.get(pos.checked_add(4)?)?,
            ]) as usize;
            pos = count_end;
            let elem_size = match elem_type {
                b'c' | b'C' => 1,
                b's' | b'S' => 2,
                b'i' | b'I' | b'f' => 4,
                _ => return None,
            };
            pos = pos.checked_add(count.checked_mul(elem_size)?)?;
        }
        _ => return None,
    }
    Some(pos)
}

#[cfg(test)]
#[allow(clippy::arithmetic_side_effects, reason = "test arithmetic on known small values")]
mod tests {
    use super::*;
    use proptest::prelude::*;

    // r[verify bam.owned_record.aux_data]
    #[test]
    fn set_and_get_string() {
        let mut aux = AuxData::new();
        aux.set_string(*b"RG", b"group1");
        assert_eq!(aux.get(*b"RG"), Some(AuxValue::String(b"group1")));
    }

    #[test]
    fn set_and_get_float() {
        let mut aux = AuxData::new();
        aux.set_float(*b"XF", std::f32::consts::PI);
        assert_eq!(aux.get(*b"XF"), Some(AuxValue::Float(std::f32::consts::PI)));
    }

    // r[verify bam.owned_record.aux_array_encoding]
    #[test]
    fn set_and_get_array_u8() {
        let mut aux = AuxData::new();
        aux.set_array_u8(*b"ML", &[10, 20, 30]).unwrap();
        let val = aux.get(*b"ML");
        assert_eq!(val, Some(AuxValue::ArrayU8(&[10, 20, 30])));
    }

    // r[verify bam.owned_record.aux_uniqueness]
    #[test]
    fn set_replaces_existing() {
        let mut aux = AuxData::new();
        aux.set_string(*b"RG", b"old");
        aux.set_string(*b"RG", b"new");
        assert_eq!(aux.get(*b"RG"), Some(AuxValue::String(b"new")));
        // Should not have duplicate
        let tags: Vec<_> = aux::iter_tags(aux.as_bytes()).collect();
        assert_eq!(tags.len(), 1);
    }

    #[test]
    fn remove_tag() {
        let mut aux = AuxData::new();
        aux.set_string(*b"RG", b"grp");
        aux.set_int(*b"NM", 5).unwrap();
        aux.remove(*b"RG");
        assert_eq!(aux.get(*b"RG"), None);
        assert_eq!(aux.get(*b"NM"), Some(AuxValue::U8(5)));
    }

    #[test]
    fn remove_nonexistent_is_noop() {
        let mut aux = AuxData::new();
        aux.set_string(*b"RG", b"grp");
        aux.remove(*b"XX");
        assert_eq!(aux.get(*b"RG"), Some(AuxValue::String(b"grp")));
    }

    // r[verify bam.owned_record.aux_int_encoding]
    #[test]
    fn int_type_selection_unsigned() {
        let mut aux = AuxData::new();

        // u8 range
        aux.set_int(*b"X0", 0).unwrap();
        assert_eq!(aux.get(*b"X0"), Some(AuxValue::U8(0)));

        aux.set_int(*b"X1", 255).unwrap();
        assert_eq!(aux.get(*b"X1"), Some(AuxValue::U8(255)));

        // u16 range
        aux.set_int(*b"X2", 256).unwrap();
        assert_eq!(aux.get(*b"X2"), Some(AuxValue::U16(256)));

        aux.set_int(*b"X3", 65535).unwrap();
        assert_eq!(aux.get(*b"X3"), Some(AuxValue::U16(65535)));

        // u32 range
        aux.set_int(*b"X4", 65536).unwrap();
        assert_eq!(aux.get(*b"X4"), Some(AuxValue::U32(65536)));

        aux.set_int(*b"X5", i64::from(u32::MAX)).unwrap();
        assert_eq!(aux.get(*b"X5"), Some(AuxValue::U32(u32::MAX)));
    }

    #[test]
    fn int_type_selection_signed() {
        let mut aux = AuxData::new();

        // i8 range
        aux.set_int(*b"X0", -1).unwrap();
        assert_eq!(aux.get(*b"X0"), Some(AuxValue::I8(-1)));

        aux.set_int(*b"X1", -128).unwrap();
        assert_eq!(aux.get(*b"X1"), Some(AuxValue::I8(-128)));

        // i16 range
        aux.set_int(*b"X2", -129).unwrap();
        assert_eq!(aux.get(*b"X2"), Some(AuxValue::I16(-129)));

        aux.set_int(*b"X3", -32768).unwrap();
        assert_eq!(aux.get(*b"X3"), Some(AuxValue::I16(-32768)));

        // i32 range
        aux.set_int(*b"X4", -32769).unwrap();
        assert_eq!(aux.get(*b"X4"), Some(AuxValue::I32(-32769)));

        aux.set_int(*b"X5", i64::from(i32::MIN)).unwrap();
        assert_eq!(aux.get(*b"X5"), Some(AuxValue::I32(i32::MIN)));
    }

    #[test]
    fn int_out_of_range() {
        let mut aux = AuxData::new();
        assert!(aux.set_int(*b"X0", i64::from(u32::MAX) + 1).is_err());
        assert!(aux.set_int(*b"X1", i64::from(i32::MIN) - 1).is_err());
    }

    // r[verify bam.owned_record.aux_from_slab]
    #[test]
    fn from_bytes_preserves_data() {
        // Build raw aux bytes manually
        let mut raw = Vec::new();
        raw.extend_from_slice(b"NM");
        raw.push(b'C');
        raw.push(42);

        let aux = AuxData::from_bytes(raw.clone());
        assert_eq!(aux.as_bytes(), &raw);
        assert_eq!(aux.get(*b"NM"), Some(AuxValue::U8(42)));
    }

    #[test]
    fn multiple_tags_coexist() {
        let mut aux = AuxData::new();
        aux.set_string(*b"RG", b"grp");
        aux.set_int(*b"NM", 3).unwrap();
        aux.set_float(*b"XF", 1.5);
        aux.set_array_u8(*b"ML", &[10, 20]).unwrap();

        assert_eq!(aux.get(*b"RG"), Some(AuxValue::String(b"grp")));
        assert_eq!(aux.get(*b"NM"), Some(AuxValue::U8(3)));
        assert_eq!(aux.get(*b"XF"), Some(AuxValue::Float(1.5)));
        assert_eq!(aux.get(*b"ML"), Some(AuxValue::ArrayU8(&[10, 20])));
    }

    #[test]
    fn remove_middle_tag_preserves_others() {
        let mut aux = AuxData::new();
        aux.set_string(*b"AA", b"first");
        aux.set_int(*b"BB", 99).unwrap();
        aux.set_string(*b"CC", b"third");

        aux.remove(*b"BB");

        assert_eq!(aux.get(*b"AA"), Some(AuxValue::String(b"first")));
        assert_eq!(aux.get(*b"BB"), None);
        assert_eq!(aux.get(*b"CC"), Some(AuxValue::String(b"third")));
    }

    #[test]
    fn is_empty() {
        let mut aux = AuxData::new();
        assert!(aux.is_empty());
        aux.set_int(*b"NM", 0).unwrap();
        assert!(!aux.is_empty());
        aux.remove(*b"NM");
        assert!(aux.is_empty());
    }

    #[test]
    fn set_and_get_char() {
        let mut aux = AuxData::new();
        aux.set_char(*b"XA", b'Q');
        assert_eq!(aux.get(*b"XA"), Some(AuxValue::Char(b'Q')));
    }

    #[test]
    fn set_and_get_double() {
        let mut aux = AuxData::new();
        aux.set_double(*b"XD", std::f64::consts::E);
        assert_eq!(aux.get(*b"XD"), Some(AuxValue::Double(std::f64::consts::E)));
    }

    // r[verify bam.owned_record.aux_array_setters]
    #[test]
    fn set_and_get_array_i16() {
        let mut aux = AuxData::new();
        let raw = [(-1i16).to_le_bytes(), 42i16.to_le_bytes()].concat();
        aux.set_array_i16(*b"XA", &raw).unwrap();
        assert_eq!(aux.get(*b"XA"), Some(AuxValue::ArrayI16(&raw)));
    }

    #[test]
    fn set_and_get_array_u16() {
        let mut aux = AuxData::new();
        let raw = [0u16.to_le_bytes(), 65535u16.to_le_bytes()].concat();
        aux.set_array_u16(*b"XA", &raw).unwrap();
        assert_eq!(aux.get(*b"XA"), Some(AuxValue::ArrayU16(&raw)));
    }

    #[test]
    fn set_and_get_array_i32() {
        let mut aux = AuxData::new();
        let raw = [i32::MIN.to_le_bytes(), 42i32.to_le_bytes()].concat();
        aux.set_array_i32(*b"XA", &raw).unwrap();
        assert_eq!(aux.get(*b"XA"), Some(AuxValue::ArrayI32(&raw)));
    }

    #[test]
    fn set_and_get_array_u32() {
        let mut aux = AuxData::new();
        let raw = [0u32.to_le_bytes(), u32::MAX.to_le_bytes()].concat();
        aux.set_array_u32(*b"XA", &raw).unwrap();
        assert_eq!(aux.get(*b"XA"), Some(AuxValue::ArrayU32(&raw)));
    }

    #[test]
    fn set_and_get_array_f32() {
        let mut aux = AuxData::new();
        let raw = [1.0f32.to_le_bytes(), (-0.5f32).to_le_bytes()].concat();
        aux.set_array_f32(*b"XA", &raw).unwrap();
        assert_eq!(aux.get(*b"XA"), Some(AuxValue::ArrayFloat(&raw)));
    }

    #[test]
    fn array_misaligned_length_is_error() {
        let mut aux = AuxData::new();
        // i16 arrays need even byte length
        assert!(aux.set_array_i16(*b"XA", &[1, 2, 3]).is_err());
        // i32 arrays need length divisible by 4
        assert!(aux.set_array_i32(*b"XA", &[1, 2, 3, 4, 5]).is_err());
        // u8 arrays allow any length
        assert!(aux.set_array_u8(*b"XA", &[1, 2, 3]).is_ok());
    }

    #[test]
    fn all_array_types_roundtrip_with_other_tags() {
        let mut aux = AuxData::new();
        aux.set_string(*b"RG", b"sample");

        let i16_raw = [100i16.to_le_bytes()].concat();
        aux.set_array_i16(*b"X1", &i16_raw).unwrap();

        let u16_raw = [50000u16.to_le_bytes()].concat();
        aux.set_array_u16(*b"X2", &u16_raw).unwrap();

        let i32_raw = [(-1_000_000i32).to_le_bytes()].concat();
        aux.set_array_i32(*b"X3", &i32_raw).unwrap();

        let u32_raw = [3_000_000_000u32.to_le_bytes()].concat();
        aux.set_array_u32(*b"X4", &u32_raw).unwrap();

        let f32_raw = [std::f32::consts::PI.to_le_bytes()].concat();
        aux.set_array_f32(*b"X5", &f32_raw).unwrap();

        assert_eq!(aux.get(*b"RG"), Some(AuxValue::String(b"sample")));
        assert_eq!(aux.get(*b"X1"), Some(AuxValue::ArrayI16(&i16_raw)));
        assert_eq!(aux.get(*b"X2"), Some(AuxValue::ArrayU16(&u16_raw)));
        assert_eq!(aux.get(*b"X3"), Some(AuxValue::ArrayI32(&i32_raw)));
        assert_eq!(aux.get(*b"X4"), Some(AuxValue::ArrayU32(&u32_raw)));
        assert_eq!(aux.get(*b"X5"), Some(AuxValue::ArrayFloat(&f32_raw)));
    }

    #[test]
    fn char_replaces_existing() {
        let mut aux = AuxData::new();
        aux.set_char(*b"XA", b'A');
        aux.set_char(*b"XA", b'B');
        assert_eq!(aux.get(*b"XA"), Some(AuxValue::Char(b'B')));
        let tags: Vec<_> = aux::iter_tags(aux.as_bytes()).collect();
        assert_eq!(tags.len(), 1);
    }

    #[test]
    fn double_replaces_existing() {
        let mut aux = AuxData::new();
        aux.set_double(*b"XD", 1.0);
        aux.set_double(*b"XD", 2.0);
        assert_eq!(aux.get(*b"XD"), Some(AuxValue::Double(2.0)));
        let tags: Vec<_> = aux::iter_tags(aux.as_bytes()).collect();
        assert_eq!(tags.len(), 1);
    }

    // ═══════════════════════════════════════════════════════════════
    // Proptests: round-trip, uniqueness, and removal properties
    // ═══════════════════════════════════════════════════════════════

    proptest::proptest! {
        // r[verify bam.owned_record.aux_data]
        /// Parse random bytes → rebuild via `AuxData` → verify tags match.
        ///
        /// This is an independent oracle: the parser (`aux::iter_tags`) and the
        /// builder (`AuxData::set_*`) use separate code paths. Round-tripping
        /// catches both encoding bugs and parse/skip divergence.
        #[test]
        fn roundtrip_random_bytes(raw in proptest::collection::vec(0u8..=255, 0..=512)) {
            // Parse: extract all well-formed tags
            let parsed: Vec<_> = aux::iter_tags(&raw).collect();

            // Rebuild via AuxData
            let mut built = AuxData::new();
            for (tag, value) in &parsed {
                match value {
                    AuxValue::Char(v) => built.set_char(*tag, *v),
                    AuxValue::String(s) => built.set_string(*tag, s),
                    AuxValue::Hex(h) => built.set_string(*tag, h),
                    AuxValue::I8(v) => drop(built.set_int(*tag, i64::from(*v))),
                    AuxValue::U8(v) => drop(built.set_int(*tag, i64::from(*v))),
                    AuxValue::I16(v) => drop(built.set_int(*tag, i64::from(*v))),
                    AuxValue::U16(v) => drop(built.set_int(*tag, i64::from(*v))),
                    AuxValue::I32(v) => drop(built.set_int(*tag, i64::from(*v))),
                    AuxValue::U32(v) => drop(built.set_int(*tag, i64::from(*v))),
                    AuxValue::Float(v) => built.set_float(*tag, *v),
                    AuxValue::Double(v) => built.set_double(*tag, *v),
                    AuxValue::ArrayI8(a) => drop(built.set_array_i16(*tag, a)),
                    AuxValue::ArrayU8(a) => drop(built.set_array_u8(*tag, a)),
                    AuxValue::ArrayI16(a) => drop(built.set_array_i16(*tag, a)),
                    AuxValue::ArrayU16(a) => drop(built.set_array_u16(*tag, a)),
                    AuxValue::ArrayI32(a) => drop(built.set_array_i32(*tag, a)),
                    AuxValue::ArrayU32(a) => drop(built.set_array_u32(*tag, a)),
                    AuxValue::ArrayFloat(a) => drop(built.set_array_f32(*tag, a)),
                }
            }

            // Re-parse the built bytes
            let rebuilt: Vec<_> = aux::iter_tags(built.as_bytes()).collect();

            // Tag count must match (H→Z preserves count, ArrayI8→ArrayI16 may fail alignment)
            let non_i8_parsed = parsed.iter().filter(|(_, v)| !matches!(v, AuxValue::ArrayI8(_))).count();
            prop_assert_eq!(rebuilt.len(), non_i8_parsed,
                "tag count mismatch: parsed vs rebuilt");

            // Verify each tag round-trips (loose comparison — int types may change)
            for (tag, original) in &parsed {
                if matches!(original, AuxValue::Hex(_) | AuxValue::ArrayI8(_)) {
                    continue;
                }
                let rt = built.get(*tag);
                prop_assert!(rt.is_some(), "tag missing after rebuild");
                let rt = rt.unwrap();
                match original {
                    AuxValue::I8(v) => prop_assert_eq!(rt.as_i64(), Some(i64::from(*v))),
                    AuxValue::U8(v) => prop_assert_eq!(rt.as_i64(), Some(i64::from(*v))),
                    AuxValue::I16(v) => prop_assert_eq!(rt.as_i64(), Some(i64::from(*v))),
                    AuxValue::U16(v) => prop_assert_eq!(rt.as_i64(), Some(i64::from(*v))),
                    AuxValue::I32(v) => prop_assert_eq!(rt.as_i64(), Some(i64::from(*v))),
                    AuxValue::U32(v) => prop_assert_eq!(rt.as_i64(), Some(i64::from(*v))),
                    _ => prop_assert_eq!(&rt, original, "tag mismatch after rebuild"),
                }
            }
        }

        // r[verify bam.owned_record.aux_uniqueness]
        /// Setting the same tag name multiple times MUST NOT create duplicates.
        #[test]
        fn set_replaces_always_unique(
            tag_names in proptest::collection::vec(
                (b'A'..=b'Z', b'A'..=b'Z'), 1..=30
            ),
            values in proptest::collection::vec(1u8..=255, 0..=256),
        ) {
            let mut aux = AuxData::new();
            for &(hi, lo) in &tag_names {
                let tag = [hi, lo];
                aux.set_string(tag, &values);
            }
            // Deduplicated count must equal unique tag count
            let unique_count = tag_names.iter().collect::<std::collections::BTreeSet<_>>().len();
            let actual_count = aux::iter_tags(aux.as_bytes()).count();
            prop_assert_eq!(actual_count, unique_count,
                "duplicate tags after set_replace");
        }

        // r[verify bam.owned_record.aux_data]
        /// Sequential removal must leave `AuxData` empty.
        #[test]
        fn remove_all_tags_leaves_empty(
            pairs in proptest::collection::vec((b'A'..=b'Z', b'A'..=b'Z', 0u8..=255), 0..=20),
        ) {
            let mut aux = AuxData::new();
            let tags: Vec<[u8; 2]> = pairs.iter().map(|&(t0, t1, _v)| [t0, t1]).collect();
            for &(t0, t1, v) in &pairs {
                aux.set_int([t0, t1], i64::from(v)).ok();
            }
            for tag in &tags {
                aux.remove(*tag);
            }
            prop_assert!(aux.is_empty());
            prop_assert_eq!(aux::iter_tags(aux.as_bytes()).count(), 0);
        }

        // r[verify bam.owned_record.aux_int_encoding]
        /// Integer set_int → get round-trip preserves the value, even
        /// though the BAM type code may differ.
        #[test]
        fn int_set_get_roundtrip(
            v in (i64::from(i32::MIN)..=i64::from(u32::MAX)),
        ) {
            let mut aux = AuxData::new();
            aux.set_int(*b"XX", v).unwrap();
            let rt = aux.get(*b"XX");
            prop_assert!(rt.is_some());
            prop_assert_eq!(rt.unwrap().as_i64(), Some(v));
        }

        // r[verify bam.owned_record.aux_array_setters]
        /// B-array setters produce bytes that the parser reads back correctly.
        #[test]
        fn array_roundtrip_u8(
            values in proptest::collection::vec(0u8..=255, 0..=64),
        ) {
            let mut aux = AuxData::new();
            aux.set_array_u8(*b"XA", &values).unwrap();
            let rt = aux.get(*b"XA");
            prop_assert_eq!(rt, Some(AuxValue::ArrayU8(&values)));
        }

        #[test]
        fn array_roundtrip_i16(
            values in proptest::collection::vec(i16::MIN..=i16::MAX, 0..=32),
        ) {
            let raw: Vec<u8> = values.iter().flat_map(|v| v.to_le_bytes()).collect();
            let mut aux = AuxData::new();
            aux.set_array_i16(*b"XB", &raw).unwrap();
            let rt = aux.get(*b"XB");
            prop_assert_eq!(rt, Some(AuxValue::ArrayI16(&raw)));
        }

        #[test]
        fn array_roundtrip_u16(
            values in proptest::collection::vec(u16::MIN..=u16::MAX, 0..=32),
        ) {
            let raw: Vec<u8> = values.iter().flat_map(|v| v.to_le_bytes()).collect();
            let mut aux = AuxData::new();
            aux.set_array_u16(*b"XC", &raw).unwrap();
            let rt = aux.get(*b"XC");
            prop_assert_eq!(rt, Some(AuxValue::ArrayU16(&raw)));
        }

        #[test]
        fn array_roundtrip_i32(
            values in proptest::collection::vec(i32::MIN..=i32::MAX, 0..=16),
        ) {
            let raw: Vec<u8> = values.iter().flat_map(|v| v.to_le_bytes()).collect();
            let mut aux = AuxData::new();
            aux.set_array_i32(*b"XD", &raw).unwrap();
            let rt = aux.get(*b"XD");
            prop_assert_eq!(rt, Some(AuxValue::ArrayI32(&raw)));
        }

        #[test]
        fn array_roundtrip_u32(
            values in proptest::collection::vec(u32::MIN..=u32::MAX, 0..=16),
        ) {
            let raw: Vec<u8> = values.iter().flat_map(|v| v.to_le_bytes()).collect();
            let mut aux = AuxData::new();
            aux.set_array_u32(*b"XE", &raw).unwrap();
            let rt = aux.get(*b"XE");
            prop_assert_eq!(rt, Some(AuxValue::ArrayU32(&raw)));
        }

        #[test]
        fn array_roundtrip_f32(
            values in proptest::collection::vec(
                proptest::num::f32::ANY, 0..=16
            ),
        ) {
            let raw: Vec<u8> = values.iter().flat_map(|v| v.to_le_bytes()).collect();
            let mut aux = AuxData::new();
            aux.set_array_f32(*b"XF", &raw).unwrap();
            let rt = aux.get(*b"XF");
            prop_assert_eq!(rt, Some(AuxValue::ArrayFloat(&raw)));
        }

        // r[verify bam.record.aux_truncated]
        /// Garbage after a valid tag stops iteration after the valid tag —
        /// the valid tag is still collected and round-tripped correctly.
        #[test]
        fn valid_tag_before_garbage_roundtrips(
            garbage in proptest::collection::vec(0u8..=255, 0..=32),
        ) {
            let mut raw = Vec::new();
            // Build a valid NM:i:42 tag
            raw.extend_from_slice(b"NM");
            raw.push(b'i');
            raw.extend_from_slice(&42i32.to_le_bytes());
            // Append garbage
            raw.extend_from_slice(&garbage);

            let parsed: Vec<_> = aux::iter_tags(&raw).collect();
            // At minimum the valid NM tag should be parsed
            prop_assert!(!parsed.is_empty(), "valid tag before garbage not parsed");
            prop_assert_eq!(parsed[0].0, *b"NM");
            prop_assert_eq!(parsed[0].1.as_i64(), Some(42));
        }
    }
}
