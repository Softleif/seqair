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
}

impl AuxData {
    /// Create empty auxiliary data.
    pub fn new() -> Self {
        Self { data: Vec::new() }
    }

    // r[impl bam.owned_record.aux_from_slab]
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
    /// Add or replace a Z-type (null-terminated string) tag.
    pub fn set_string(&mut self, tag: [u8; 2], value: &[u8]) {
        self.remove(tag);
        self.data.extend_from_slice(&tag);
        self.data.push(b'Z');
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

    // r[impl bam.owned_record.aux_array_encoding]
    /// Add or replace a B:C (unsigned byte array) tag.
    ///
    /// The BAM B-array element count is u32. In practice, array size is bounded by
    /// the 2 MiB record size limit, so values longer than u32::MAX cannot occur in
    /// valid BAM data. We validate the cast to be safe on 64-bit platforms.
    pub fn set_array_u8(&mut self, tag: [u8; 2], values: &[u8]) -> Result<(), AuxDataError> {
        let count = u32::try_from(values.len())
            .map_err(|_| AuxDataError::IntegerOutOfRange { value: values.len() as i64 })?;
        self.remove(tag);
        self.data.extend_from_slice(&tag);
        self.data.push(b'B');
        self.data.push(b'C');
        self.data.extend_from_slice(&count.to_le_bytes());
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
        aux.set_float(*b"XF", 3.14);
        assert_eq!(aux.get(*b"XF"), Some(AuxValue::Float(3.14)));
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
}
