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

    /// `set_char` received a byte outside the SAM A-type printable-ASCII grammar
    /// `[!-~]` (0x21..=0x7E).
    #[error("invalid A-type byte {value:#04x}: SAM grammar requires printable ASCII [!-~]")]
    InvalidCharByte { value: u8 },
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
    ///
    /// `value` MUST be in the SAM grammar `[!-~]` (printable ASCII,
    /// 0x21..=0x7E). Bytes outside that range return
    /// [`AuxDataError::InvalidCharByte`] and the underlying buffer is left
    /// untouched (no orphaned bytes — validation runs before any mutation).
    pub fn set_char(&mut self, tag: [u8; 2], value: u8) -> Result<(), AuxDataError> {
        if !matches!(value, 0x21..=0x7E) {
            return Err(AuxDataError::InvalidCharByte { value });
        }
        self.remove(tag);
        self.data.extend_from_slice(&tag);
        self.data.push(b'A');
        self.data.push(value);
        Ok(())
    }

    /// Add or replace an H-type (hex string) tag.
    ///
    /// The H type is encoded like Z (NUL-terminated), but with type byte `'H'`.
    /// `value` is written verbatim; the SAM grammar restricts it to even-length
    /// `[0-9A-F]` but the writer is lenient — the BAM wire format does not
    /// constrain the bytes, only the parser does.
    pub fn set_hex(&mut self, tag: [u8; 2], value: &[u8]) {
        self.remove(tag);
        self.data.extend_from_slice(&tag);
        self.data.push(b'H');
        self.data.extend_from_slice(value);
        self.data.push(0); // NUL terminator
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
    /// Add or replace a B:C (unsigned 8-bit array) tag.
    ///
    /// The BAM B-array element count is u32. In practice, array size is bounded by
    /// the 2 MiB record size limit, so values longer than `u32::MAX` cannot occur in
    /// valid BAM data. We validate the cast to be safe on 64-bit platforms.
    pub fn set_array_u8(&mut self, tag: [u8; 2], values: &[u8]) -> Result<(), AuxDataError> {
        let count = array_count_from_slice(values.len())?;
        self.remove(tag);
        encode_array_header(&mut self.data, &tag, b'C', count);
        self.data.extend_from_slice(values);
        Ok(())
    }

    // r[impl bam.owned_record.aux_array_setters]
    /// Add or replace a B:c (signed 8-bit array) tag.
    pub fn set_array_i8(&mut self, tag: [u8; 2], values: &[i8]) -> Result<(), AuxDataError> {
        let count = array_count_from_slice(values.len())?;
        self.remove(tag);
        encode_array_header(&mut self.data, &tag, b'c', count);
        // i8 → u8 is a sign-preserving 1:1 byte reinterpretation (LE-irrelevant: 1 byte).
        for &v in values {
            self.data.push(v.cast_unsigned());
        }
        Ok(())
    }

    /// Add or replace a B:s (signed 16-bit array) tag.
    pub fn set_array_i16(&mut self, tag: [u8; 2], values: &[i16]) -> Result<(), AuxDataError> {
        let count = array_count_from_slice(values.len())?;
        self.remove(tag);
        encode_array_header(&mut self.data, &tag, b's', count);
        for &v in values {
            self.data.extend_from_slice(&v.to_le_bytes());
        }
        Ok(())
    }

    /// Add or replace a B:S (unsigned 16-bit array) tag.
    pub fn set_array_u16(&mut self, tag: [u8; 2], values: &[u16]) -> Result<(), AuxDataError> {
        let count = array_count_from_slice(values.len())?;
        self.remove(tag);
        encode_array_header(&mut self.data, &tag, b'S', count);
        for &v in values {
            self.data.extend_from_slice(&v.to_le_bytes());
        }
        Ok(())
    }

    /// Add or replace a B:i (signed 32-bit array) tag.
    pub fn set_array_i32(&mut self, tag: [u8; 2], values: &[i32]) -> Result<(), AuxDataError> {
        let count = array_count_from_slice(values.len())?;
        self.remove(tag);
        encode_array_header(&mut self.data, &tag, b'i', count);
        for &v in values {
            self.data.extend_from_slice(&v.to_le_bytes());
        }
        Ok(())
    }

    /// Add or replace a B:I (unsigned 32-bit array) tag.
    pub fn set_array_u32(&mut self, tag: [u8; 2], values: &[u32]) -> Result<(), AuxDataError> {
        let count = array_count_from_slice(values.len())?;
        self.remove(tag);
        encode_array_header(&mut self.data, &tag, b'I', count);
        for &v in values {
            self.data.extend_from_slice(&v.to_le_bytes());
        }
        Ok(())
    }

    /// Add or replace a B:f (float array) tag.
    pub fn set_array_f32(&mut self, tag: [u8; 2], values: &[f32]) -> Result<(), AuxDataError> {
        let count = array_count_from_slice(values.len())?;
        self.remove(tag);
        encode_array_header(&mut self.data, &tag, b'f', count);
        for &v in values {
            self.data.extend_from_slice(&v.to_le_bytes());
        }
        Ok(())
    }

    /// Remove a tag if present.
    pub fn remove(&mut self, tag: [u8; 2]) {
        if let Some(range) = find_tag_byte_range(&self.data, tag) {
            self.data.drain(range);
        }
    }
}

/// Validate that a typed-slice element count fits in u32 (BAM B-array `count` is u32).
///
/// On 64-bit platforms, `usize` is wider than `u32`, so we must guard against
/// arrays whose element count exceeds the BAM wire-format limit.
fn array_count_from_slice(elem_count: usize) -> Result<u32, AuxDataError> {
    u32::try_from(elem_count).map_err(|_| AuxDataError::IntegerOutOfRange {
        #[allow(clippy::cast_possible_wrap, reason = "count fits in i64 on 64-bit")]
        value: elem_count as i64,
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
        aux.set_char(*b"XA", b'Q').unwrap();
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
    fn set_and_get_array_i8() {
        let mut aux = AuxData::new();
        aux.set_array_i8(*b"XA", &[-1i8, 0, 127, -128]).unwrap();
        let expected = [(-1i8).cast_unsigned(), 0, 127, (-128i8).cast_unsigned()];
        assert_eq!(aux.get(*b"XA"), Some(AuxValue::ArrayI8(&expected)));
    }

    #[test]
    fn set_and_get_array_i16() {
        let mut aux = AuxData::new();
        aux.set_array_i16(*b"XA", &[-1i16, 42]).unwrap();
        let expected = [(-1i16).to_le_bytes(), 42i16.to_le_bytes()].concat();
        assert_eq!(aux.get(*b"XA"), Some(AuxValue::ArrayI16(&expected)));
    }

    #[test]
    fn set_and_get_array_u16() {
        let mut aux = AuxData::new();
        aux.set_array_u16(*b"XA", &[0u16, 65535]).unwrap();
        let expected = [0u16.to_le_bytes(), 65535u16.to_le_bytes()].concat();
        assert_eq!(aux.get(*b"XA"), Some(AuxValue::ArrayU16(&expected)));
    }

    #[test]
    fn set_and_get_array_i32() {
        let mut aux = AuxData::new();
        aux.set_array_i32(*b"XA", &[i32::MIN, 42]).unwrap();
        let expected = [i32::MIN.to_le_bytes(), 42i32.to_le_bytes()].concat();
        assert_eq!(aux.get(*b"XA"), Some(AuxValue::ArrayI32(&expected)));
    }

    #[test]
    fn set_and_get_array_u32() {
        let mut aux = AuxData::new();
        aux.set_array_u32(*b"XA", &[0u32, u32::MAX]).unwrap();
        let expected = [0u32.to_le_bytes(), u32::MAX.to_le_bytes()].concat();
        assert_eq!(aux.get(*b"XA"), Some(AuxValue::ArrayU32(&expected)));
    }

    #[test]
    fn set_and_get_array_f32() {
        let mut aux = AuxData::new();
        aux.set_array_f32(*b"XA", &[1.0f32, -0.5]).unwrap();
        let expected = [1.0f32.to_le_bytes(), (-0.5f32).to_le_bytes()].concat();
        assert_eq!(aux.get(*b"XA"), Some(AuxValue::ArrayFloat(&expected)));
    }

    #[test]
    fn all_array_types_roundtrip_with_other_tags() {
        let mut aux = AuxData::new();
        aux.set_string(*b"RG", b"sample");

        aux.set_array_i8(*b"X0", &[-1i8, 0, 1]).unwrap();
        aux.set_array_i16(*b"X1", &[100i16]).unwrap();
        aux.set_array_u16(*b"X2", &[50000u16]).unwrap();
        aux.set_array_i32(*b"X3", &[-1_000_000i32]).unwrap();
        aux.set_array_u32(*b"X4", &[3_000_000_000u32]).unwrap();
        aux.set_array_f32(*b"X5", &[std::f32::consts::PI]).unwrap();

        let i8_raw = [(-1i8).cast_unsigned(), 0u8, 1u8];
        let i16_raw = 100i16.to_le_bytes();
        let u16_raw = 50000u16.to_le_bytes();
        let i32_raw = (-1_000_000i32).to_le_bytes();
        let u32_raw = 3_000_000_000u32.to_le_bytes();
        let f32_raw = std::f32::consts::PI.to_le_bytes();

        assert_eq!(aux.get(*b"RG"), Some(AuxValue::String(b"sample")));
        assert_eq!(aux.get(*b"X0"), Some(AuxValue::ArrayI8(&i8_raw)));
        assert_eq!(aux.get(*b"X1"), Some(AuxValue::ArrayI16(&i16_raw)));
        assert_eq!(aux.get(*b"X2"), Some(AuxValue::ArrayU16(&u16_raw)));
        assert_eq!(aux.get(*b"X3"), Some(AuxValue::ArrayI32(&i32_raw)));
        assert_eq!(aux.get(*b"X4"), Some(AuxValue::ArrayU32(&u32_raw)));
        assert_eq!(aux.get(*b"X5"), Some(AuxValue::ArrayFloat(&f32_raw)));
    }

    // r[verify bam.owned_record.aux_data]
    /// `set_char` rejects bytes outside `[!-~]` and does not mutate state.
    #[test]
    fn set_char_rejects_non_printable() {
        let mut aux = AuxData::new();
        aux.set_string(*b"RG", b"grp"); // unrelated tag

        // Pre-mutation snapshot.
        let before = aux.as_bytes().to_vec();

        // Below 0x21 (space, control bytes) → reject.
        assert!(matches!(
            aux.set_char(*b"XA", 0x20),
            Err(AuxDataError::InvalidCharByte { value: 0x20 })
        ));
        // 0x7F (DEL) → reject.
        assert!(matches!(
            aux.set_char(*b"XA", 0x7F),
            Err(AuxDataError::InvalidCharByte { value: 0x7F })
        ));
        // High-bit set → reject.
        assert!(matches!(
            aux.set_char(*b"XA", 0xC3),
            Err(AuxDataError::InvalidCharByte { value: 0xC3 })
        ));

        // Boundaries inclusive.
        aux.set_char(*b"XB", 0x21).unwrap();
        aux.set_char(*b"XC", 0x7E).unwrap();

        // No mutation from the rejected calls — only the boundary insertions added bytes.
        assert!(aux.as_bytes().starts_with(&before));
    }

    #[test]
    fn char_replaces_existing() {
        let mut aux = AuxData::new();
        aux.set_char(*b"XA", b'A').unwrap();
        aux.set_char(*b"XA", b'B').unwrap();
        assert_eq!(aux.get(*b"XA"), Some(AuxValue::Char(b'B')));
        let tags: Vec<_> = aux::iter_tags(aux.as_bytes()).collect();
        assert_eq!(tags.len(), 1);
    }

    // r[verify bam.owned_record.aux_data]
    /// H-type round-trip: bytes go in, parser returns them as `Hex`.
    #[test]
    fn set_and_get_hex() {
        let mut aux = AuxData::new();
        aux.set_hex(*b"XH", b"DEADBEEF");
        assert_eq!(aux.get(*b"XH"), Some(AuxValue::Hex(b"DEADBEEF")));
    }

    #[test]
    fn hex_replaces_existing() {
        let mut aux = AuxData::new();
        aux.set_hex(*b"XH", b"AA");
        aux.set_hex(*b"XH", b"FF");
        assert_eq!(aux.get(*b"XH"), Some(AuxValue::Hex(b"FF")));
        let tags: Vec<_> = aux::iter_tags(aux.as_bytes()).collect();
        assert_eq!(tags.len(), 1);
    }

    #[test]
    fn hex_does_not_collide_with_z() {
        // Hex and String share the same parser path (NUL-terminated) but
        // produce different `AuxValue` variants based on the type byte.
        let mut aux = AuxData::new();
        aux.set_hex(*b"XH", b"CAFE");
        aux.set_string(*b"XZ", b"CAFE");
        assert_eq!(aux.get(*b"XH"), Some(AuxValue::Hex(b"CAFE")));
        assert_eq!(aux.get(*b"XZ"), Some(AuxValue::String(b"CAFE")));
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
        ///
        /// `Char` tags whose byte falls outside the printable-ASCII grammar
        /// `[!-~]` are intentionally skipped — `set_char` rejects them by
        /// design, and that path is covered by `set_char_rejects_non_printable`.
        #[test]
        fn roundtrip_random_bytes(raw in proptest::collection::vec(0u8..=255, 0..=512)) {
            // Parse: extract all well-formed tags
            let parsed: Vec<_> = aux::iter_tags(&raw).collect();

            // Rebuild via AuxData using only the typed setters that match each variant.
            let mut built = AuxData::new();
            let mut skipped_char_tags = std::collections::BTreeSet::new();
            for (tag, value) in &parsed {
                match value {
                    AuxValue::Char(v) => {
                        // Random bytes are usually not printable ASCII; the writer
                        // rejects them. Track these so we don't expect them to round-trip.
                        if built.set_char(*tag, *v).is_err() {
                            skipped_char_tags.insert(*tag);
                        }
                    }
                    AuxValue::String(s) => built.set_string(*tag, s),
                    AuxValue::Hex(h) => built.set_hex(*tag, h),
                    AuxValue::I8(v) => drop(built.set_int(*tag, i64::from(*v))),
                    AuxValue::U8(v) => drop(built.set_int(*tag, i64::from(*v))),
                    AuxValue::I16(v) => drop(built.set_int(*tag, i64::from(*v))),
                    AuxValue::U16(v) => drop(built.set_int(*tag, i64::from(*v))),
                    AuxValue::I32(v) => drop(built.set_int(*tag, i64::from(*v))),
                    AuxValue::U32(v) => drop(built.set_int(*tag, i64::from(*v))),
                    AuxValue::Float(v) => built.set_float(*tag, *v),
                    AuxValue::Double(v) => built.set_double(*tag, *v),
                    AuxValue::ArrayI8(a) => {
                        let typed: Vec<i8> = a.iter().map(|&b| b.cast_signed()).collect();
                        drop(built.set_array_i8(*tag, &typed));
                    }
                    AuxValue::ArrayU8(a) => drop(built.set_array_u8(*tag, a)),
                    AuxValue::ArrayI16(a) => {
                        let typed: Vec<i16> = a.chunks_exact(2)
                            .map(|c| i16::from_le_bytes([c[0], c[1]]))
                            .collect();
                        drop(built.set_array_i16(*tag, &typed));
                    }
                    AuxValue::ArrayU16(a) => {
                        let typed: Vec<u16> = a.chunks_exact(2)
                            .map(|c| u16::from_le_bytes([c[0], c[1]]))
                            .collect();
                        drop(built.set_array_u16(*tag, &typed));
                    }
                    AuxValue::ArrayI32(a) => {
                        let typed: Vec<i32> = a.chunks_exact(4)
                            .map(|c| i32::from_le_bytes([c[0], c[1], c[2], c[3]]))
                            .collect();
                        drop(built.set_array_i32(*tag, &typed));
                    }
                    AuxValue::ArrayU32(a) => {
                        let typed: Vec<u32> = a.chunks_exact(4)
                            .map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]]))
                            .collect();
                        drop(built.set_array_u32(*tag, &typed));
                    }
                    AuxValue::ArrayFloat(a) => {
                        let typed: Vec<f32> = a.chunks_exact(4)
                            .map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]]))
                            .collect();
                        drop(built.set_array_f32(*tag, &typed));
                    }
                }
            }

            // Re-parse the built bytes
            let rebuilt: Vec<_> = aux::iter_tags(built.as_bytes()).collect();

            // Tag count must match (after dedup-by-name and after dropping skipped chars).
            let mut expected_tags: std::collections::BTreeSet<[u8; 2]> = parsed.iter()
                .map(|(t, _)| *t)
                .collect();
            for t in &skipped_char_tags { expected_tags.remove(t); }
            prop_assert_eq!(rebuilt.len(), expected_tags.len(),
                "tag count mismatch: parsed (deduped, minus rejected chars) vs rebuilt");

            // Verify each tag round-trips. For repeated-tag inputs, the LAST occurrence wins;
            // walk forward and let later writes overwrite the expected value.
            let mut expected: std::collections::BTreeMap<[u8; 2], &AuxValue<'_>> =
                std::collections::BTreeMap::new();
            for (tag, value) in &parsed {
                if matches!(value, AuxValue::Char(_)) && skipped_char_tags.contains(tag) {
                    expected.remove(tag);
                    continue;
                }
                expected.insert(*tag, value);
            }
            for (tag, original) in &expected {
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
                    _ => prop_assert_eq!(&rt, *original, "tag mismatch after rebuild"),
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
        fn array_roundtrip_i8(
            values in proptest::collection::vec(i8::MIN..=i8::MAX, 0..=64),
        ) {
            let mut aux = AuxData::new();
            aux.set_array_i8(*b"XA", &values).unwrap();
            let raw: Vec<u8> = values.iter().map(|&v| v.cast_unsigned()).collect();
            let rt = aux.get(*b"XA");
            prop_assert_eq!(rt, Some(AuxValue::ArrayI8(&raw)));
        }

        #[test]
        fn array_roundtrip_i16(
            values in proptest::collection::vec(i16::MIN..=i16::MAX, 0..=32),
        ) {
            let mut aux = AuxData::new();
            aux.set_array_i16(*b"XB", &values).unwrap();
            let raw: Vec<u8> = values.iter().flat_map(|v| v.to_le_bytes()).collect();
            let rt = aux.get(*b"XB");
            prop_assert_eq!(rt, Some(AuxValue::ArrayI16(&raw)));
        }

        #[test]
        fn array_roundtrip_u16(
            values in proptest::collection::vec(u16::MIN..=u16::MAX, 0..=32),
        ) {
            let mut aux = AuxData::new();
            aux.set_array_u16(*b"XC", &values).unwrap();
            let raw: Vec<u8> = values.iter().flat_map(|v| v.to_le_bytes()).collect();
            let rt = aux.get(*b"XC");
            prop_assert_eq!(rt, Some(AuxValue::ArrayU16(&raw)));
        }

        #[test]
        fn array_roundtrip_i32(
            values in proptest::collection::vec(i32::MIN..=i32::MAX, 0..=16),
        ) {
            let mut aux = AuxData::new();
            aux.set_array_i32(*b"XD", &values).unwrap();
            let raw: Vec<u8> = values.iter().flat_map(|v| v.to_le_bytes()).collect();
            let rt = aux.get(*b"XD");
            prop_assert_eq!(rt, Some(AuxValue::ArrayI32(&raw)));
        }

        #[test]
        fn array_roundtrip_u32(
            values in proptest::collection::vec(u32::MIN..=u32::MAX, 0..=16),
        ) {
            let mut aux = AuxData::new();
            aux.set_array_u32(*b"XE", &values).unwrap();
            let raw: Vec<u8> = values.iter().flat_map(|v| v.to_le_bytes()).collect();
            let rt = aux.get(*b"XE");
            prop_assert_eq!(rt, Some(AuxValue::ArrayU32(&raw)));
        }

        #[test]
        fn array_roundtrip_f32(
            values in proptest::collection::vec(
                proptest::num::f32::ANY, 0..=16
            ),
        ) {
            let mut aux = AuxData::new();
            aux.set_array_f32(*b"XF", &values).unwrap();
            let raw: Vec<u8> = values.iter().flat_map(|v| v.to_le_bytes()).collect();
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
