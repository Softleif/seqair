//! BAM auxiliary tag parsing and [`FromAuxValue`] conversion.
//!
//! Parses raw BAM auxiliary data bytes into typed [`AuxValue`] values.
//! The [`Aux`] wrapper provides a fluent [`Aux::get`] API with automatic
//! type conversion via [`FromAuxValue`].

use seqair_types::SmolStr;

/// A single BAM auxiliary tag value.
#[derive(Debug, Clone, PartialEq)]
pub enum AuxValue<'a> {
    /// Printable character (`A` type code).
    Char(u8),
    /// Signed 8-bit integer (`c` type code).
    I8(i8),
    /// Unsigned 8-bit integer (`C` type code).
    U8(u8),
    /// Signed 16-bit integer (`s` type code).
    I16(i16),
    /// Unsigned 16-bit integer (`S` type code).
    U16(u16),
    /// Signed 32-bit integer (`i` type code).
    I32(i32),
    /// Unsigned 32-bit integer (`I` type code).
    U32(u32),
    /// Single-precision float (`f` type code).
    Float(f32),
    /// Double-precision float (`d` type code).
    Double(f64),
    /// Null-terminated string (`Z` type code), without the terminator.
    String(&'a [u8]),
    /// Hex string (`H` type code), without the null terminator.
    Hex(&'a [u8]),
    /// Array of signed 8-bit integers (`B`/`c` subtype). Raw little-endian bytes.
    ArrayI8(&'a [u8]),
    /// Array of unsigned 8-bit integers (`B`/`C` subtype). Raw bytes.
    ArrayU8(&'a [u8]),
    /// Array of signed 16-bit integers (`B`/`s` subtype). Raw little-endian bytes.
    ArrayI16(&'a [u8]),
    /// Array of unsigned 16-bit integers (`B`/`S` subtype). Raw little-endian bytes.
    ArrayU16(&'a [u8]),
    /// Array of signed 32-bit integers (`B`/`i` subtype). Raw little-endian bytes.
    ArrayI32(&'a [u8]),
    /// Array of unsigned 32-bit integers (`B`/`I` subtype). Raw little-endian bytes.
    ArrayU32(&'a [u8]),
    /// Array of single-precision floats (`B`/`f` subtype). Raw little-endian bytes.
    ArrayFloat(&'a [u8]),
}

impl<'a> AuxValue<'a> {
    /// If this is a `String` value, return the raw bytes (without null terminator).
    pub fn as_str(&self) -> Option<&'a [u8]> {
        match self {
            AuxValue::String(s) => Some(s),
            _ => None,
        }
    }

    /// Convert any integer variant to `i64`. Returns `None` for non-integer types.
    pub fn as_i64(&self) -> Option<i64> {
        match *self {
            AuxValue::I8(v) => Some(i64::from(v)),
            AuxValue::U8(v) => Some(i64::from(v)),
            AuxValue::I16(v) => Some(i64::from(v)),
            AuxValue::U16(v) => Some(i64::from(v)),
            AuxValue::I32(v) => Some(i64::from(v)),
            AuxValue::U32(v) => Some(i64::from(v)),
            _ => None,
        }
    }

    /// Human-readable type name for error messages (e.g. "i32", "Z").
    pub fn type_name(&self) -> &'static str {
        match self {
            AuxValue::Char(_) => "A",
            AuxValue::I8(_) => "c",
            AuxValue::U8(_) => "C",
            AuxValue::I16(_) => "s",
            AuxValue::U16(_) => "S",
            AuxValue::I32(_) => "i",
            AuxValue::U32(_) => "I",
            AuxValue::Float(_) => "f",
            AuxValue::Double(_) => "d",
            AuxValue::String(_) => "Z",
            AuxValue::Hex(_) => "H",
            AuxValue::ArrayI8(_) => "B:c",
            AuxValue::ArrayU8(_) => "B:C",
            AuxValue::ArrayI16(_) => "B:s",
            AuxValue::ArrayU16(_) => "B:S",
            AuxValue::ArrayI32(_) => "B:i",
            AuxValue::ArrayU32(_) => "B:I",
            AuxValue::ArrayFloat(_) => "B:f",
        }
    }
}

/// Find a tag in raw BAM auxiliary data and return its typed value.
///
/// Returns `None` if the tag is not found or the data is truncated.
pub fn find_tag<'a>(aux: &'a [u8], tag: [u8; 2]) -> Option<AuxValue<'a>> {
    for (t, value) in iter_tags(aux) {
        if t == tag {
            return Some(value);
        }
    }
    None
}

// r[impl bam.record.aux_get_error]
/// Error from [`Aux::get`].
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum GetAuxError {
    /// The requested tag is not present in the auxiliary data.
    #[error("tag not found: {tag:?}")]
    TagNotFound { tag: [u8; 2] },
    /// Tag name must be exactly 2 bytes (BAM format requirement).
    #[error("invalid tag name length: expected 2, got {len}")]
    InvalidTagName { len: usize },
    /// The tag exists but has a different BAM type than requested.
    #[error("type mismatch: expected {expected}, got {actual}")]
    TypeMismatch { expected: &'static str, actual: &'static str },
    /// A `Z`-type string contains bytes that are not valid UTF-8.
    #[error("Z-type string contains invalid UTF-8")]
    InvalidUtf8,
}

// r[impl bam.record.aux_from_aux_value]
/// Convert an [`AuxValue`] into a concrete Rust type.
///
/// Implementations perform type checking and (for integers) widening
/// conversions. Implemented for `i64`, `u64`, `f64`, `&str`, `&[u8]`,
/// `SmolStr`, `String`, `u8`, `u16`, `u32`, `i32`, `f32`, and `char`.
///
/// # Examples
///
/// ```
/// use seqair::bam::aux::{AuxValue, FromAuxValue};
///
/// let val = AuxValue::U8(42);
/// let wide: i64 = FromAuxValue::from_aux_value(val).unwrap();
/// assert_eq!(wide, 42);
/// ```
pub trait FromAuxValue<'a>: Sized {
    fn from_aux_value(value: AuxValue<'a>) -> Result<Self, GetAuxError>;
}

// r[impl bam.record.aux_widening]
impl<'a> FromAuxValue<'a> for i64 {
    fn from_aux_value(value: AuxValue<'a>) -> Result<Self, GetAuxError> {
        value.as_i64().ok_or_else(|| GetAuxError::TypeMismatch {
            expected: "integer",
            actual: value.type_name(),
        })
    }
}

impl<'a> FromAuxValue<'a> for i32 {
    fn from_aux_value(value: AuxValue<'a>) -> Result<Self, GetAuxError> {
        match value {
            AuxValue::I8(v) => Ok(i32::from(v)),
            AuxValue::U8(v) => Ok(i32::from(v)),
            AuxValue::I16(v) => Ok(i32::from(v)),
            AuxValue::U16(v) => Ok(i32::from(v)),
            AuxValue::I32(v) => Ok(v),
            AuxValue::U32(v) => i32::try_from(v)
                .map_err(|_| GetAuxError::TypeMismatch { expected: "i32", actual: "U32" }),
            ref other => {
                Err(GetAuxError::TypeMismatch { expected: "integer", actual: other.type_name() })
            }
        }
    }
}

impl<'a> FromAuxValue<'a> for u64 {
    fn from_aux_value(value: AuxValue<'a>) -> Result<Self, GetAuxError> {
        match value {
            AuxValue::U8(v) => Ok(u64::from(v)),
            AuxValue::U16(v) => Ok(u64::from(v)),
            AuxValue::U32(v) => Ok(u64::from(v)),
            ref other => Err(GetAuxError::TypeMismatch {
                expected: "unsigned integer",
                actual: other.type_name(),
            }),
        }
    }
}

impl<'a> FromAuxValue<'a> for u32 {
    fn from_aux_value(value: AuxValue<'a>) -> Result<Self, GetAuxError> {
        match value {
            AuxValue::U8(v) => Ok(u32::from(v)),
            AuxValue::U16(v) => Ok(u32::from(v)),
            AuxValue::U32(v) => Ok(v),
            ref other => Err(GetAuxError::TypeMismatch {
                expected: "unsigned integer",
                actual: other.type_name(),
            }),
        }
    }
}

impl<'a> FromAuxValue<'a> for u16 {
    fn from_aux_value(value: AuxValue<'a>) -> Result<Self, GetAuxError> {
        match value {
            AuxValue::U8(v) => Ok(u16::from(v)),
            AuxValue::U16(v) => Ok(v),
            ref other => Err(GetAuxError::TypeMismatch {
                expected: "unsigned integer",
                actual: other.type_name(),
            }),
        }
    }
}

impl<'a> FromAuxValue<'a> for u8 {
    fn from_aux_value(value: AuxValue<'a>) -> Result<Self, GetAuxError> {
        match value {
            AuxValue::U8(v) => Ok(v),
            ref other => {
                Err(GetAuxError::TypeMismatch { expected: "U8", actual: other.type_name() })
            }
        }
    }
}

impl<'a> FromAuxValue<'a> for f64 {
    fn from_aux_value(value: AuxValue<'a>) -> Result<Self, GetAuxError> {
        match value {
            AuxValue::Float(v) => Ok(f64::from(v)),
            AuxValue::Double(v) => Ok(v),
            ref other => Err(GetAuxError::TypeMismatch {
                expected: "float or double",
                actual: other.type_name(),
            }),
        }
    }
}

impl<'a> FromAuxValue<'a> for f32 {
    fn from_aux_value(value: AuxValue<'a>) -> Result<Self, GetAuxError> {
        match value {
            AuxValue::Float(v) => Ok(v),
            ref other => {
                Err(GetAuxError::TypeMismatch { expected: "float", actual: other.type_name() })
            }
        }
    }
}

impl<'a> FromAuxValue<'a> for char {
    fn from_aux_value(value: AuxValue<'a>) -> Result<Self, GetAuxError> {
        match value {
            AuxValue::Char(v) => Ok(v as char),
            ref other => {
                Err(GetAuxError::TypeMismatch { expected: "A", actual: other.type_name() })
            }
        }
    }
}

impl<'a> FromAuxValue<'a> for &'a [u8] {
    fn from_aux_value(value: AuxValue<'a>) -> Result<Self, GetAuxError> {
        match value {
            AuxValue::String(s) => Ok(s),
            ref other => {
                Err(GetAuxError::TypeMismatch { expected: "Z", actual: other.type_name() })
            }
        }
    }
}

impl<'a> FromAuxValue<'a> for &'a str {
    fn from_aux_value(value: AuxValue<'a>) -> Result<Self, GetAuxError> {
        match value {
            AuxValue::String(s) => std::str::from_utf8(s).map_err(|_| GetAuxError::InvalidUtf8),
            ref other => {
                Err(GetAuxError::TypeMismatch { expected: "Z", actual: other.type_name() })
            }
        }
    }
}

impl<'a> FromAuxValue<'a> for String {
    fn from_aux_value(value: AuxValue<'a>) -> Result<Self, GetAuxError> {
        let s: &str = FromAuxValue::from_aux_value(value)?;
        Ok(s.to_owned())
    }
}

impl<'a> FromAuxValue<'a> for SmolStr {
    fn from_aux_value(value: AuxValue<'a>) -> Result<Self, GetAuxError> {
        let s: &str = FromAuxValue::from_aux_value(value)?;
        Ok(SmolStr::from(s))
    }
}

// r[impl bam.record.aux_wrapper]
/// Borrowed view of a record's auxiliary tag bytes.
///
/// Provides a fluent [`get`](Aux::get) API with automatic type conversion
/// via [`FromAuxValue`]. Derefs to `[u8]` for backward compatibility.
///
/// # Examples
///
/// ```
/// use seqair::bam::aux::Aux;
///
/// let raw = vec![b'R', b'G', b'Z', b'g', b'r', b'p', b'1', 0];
/// let aux = Aux::new(&raw);
///
/// // Type inference from context:
/// let rg: &str = aux.get("RG").unwrap();
/// assert_eq!(rg, "grp1");
///
/// // Owned type:
/// let rg2: String = aux.get("RG").unwrap();
/// assert_eq!(rg2, "grp1");
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Aux<'a> {
    data: &'a [u8],
}

impl<'a> Aux<'a> {
    /// Wrap raw BAM auxiliary bytes.
    pub const fn new(data: &'a [u8]) -> Self {
        Self { data }
    }

    /// Raw bytes for pass-through / backward compatibility.
    pub fn as_bytes(&self) -> &'a [u8] {
        self.data
    }

    /// Iterate over all tags.
    pub fn iter(&self) -> AuxIter<'a> {
        iter_tags(self.data)
    }

    // r[impl bam.record.aux_get]
    /// Look up a tag by name and convert it to `T` via [`FromAuxValue`].
    ///
    /// The tag name is validated to exactly 2 bytes (BAM format requirement).
    /// String literals like `"RG"` are accepted.
    ///
    /// # Errors
    ///
    /// - [`GetAuxError::TagNotFound`] if the tag is not present.
    /// - [`GetAuxError::InvalidTagName`] if the tag name is not exactly 2 bytes.
    /// - [`GetAuxError::TypeMismatch`] if the tag has a different BAM type.
    /// - [`GetAuxError::InvalidUtf8`] if a `Z`-type string is not valid UTF-8.
    pub fn get<T: FromAuxValue<'a>>(&self, tag: impl AsRef<[u8]>) -> Result<T, GetAuxError> {
        let tag_bytes = tag.as_ref();
        let tag_arr: [u8; 2] = tag_bytes
            .try_into()
            .map_err(|_| GetAuxError::InvalidTagName { len: tag_bytes.len() })?;
        let value =
            find_tag(self.data, tag_arr).ok_or(GetAuxError::TagNotFound { tag: tag_arr })?;
        T::from_aux_value(value)
    }
}

impl<'a> std::ops::Deref for Aux<'a> {
    type Target = [u8];

    fn deref(&self) -> &Self::Target {
        self.data
    }
}

/// Iterator over all tags in raw BAM auxiliary data.
pub struct AuxIter<'a> {
    data: &'a [u8],
    pos: usize,
}

/// Iterate over all (tag, value) pairs in raw BAM auxiliary data.
///
/// Silently stops on truncated or malformed data.
pub fn iter_tags(aux: &[u8]) -> AuxIter<'_> {
    AuxIter { data: aux, pos: 0 }
}

impl<'a> Iterator for AuxIter<'a> {
    type Item = ([u8; 2], AuxValue<'a>);

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos.saturating_add(3) > self.data.len() {
            return None;
        }

        let data = self.data;
        let tag = [*data.get(self.pos)?, *data.get(self.pos.saturating_add(1))?];
        let typ = *data.get(self.pos.saturating_add(2))?;
        self.pos = self.pos.checked_add(3)?;

        let value = self.parse_value(typ)?;
        Some((tag, value))
    }
}

impl<'a> AuxIter<'a> {
    /// Returns `Some(value)` for a successfully parsed tag,
    /// or `None` for unknown type codes or malformed/truncated data.
    fn parse_value(&mut self, typ: u8) -> Option<AuxValue<'a>> {
        match typ {
            b'A' => {
                let v = *self.data.get(self.pos)?;
                self.pos = self.pos.checked_add(1)?;
                Some(AuxValue::Char(v))
            }
            b'c' => {
                let v = *self.data.get(self.pos)?;
                self.pos = self.pos.checked_add(1)?;
                Some(AuxValue::I8(v.cast_signed()))
            }
            b'C' => {
                let v = *self.data.get(self.pos)?;
                self.pos = self.pos.checked_add(1)?;
                Some(AuxValue::U8(v))
            }
            b's' => {
                let bytes = self.read_bytes::<2>()?;
                Some(AuxValue::I16(i16::from_le_bytes(bytes)))
            }
            b'S' => {
                let bytes = self.read_bytes::<2>()?;
                Some(AuxValue::U16(u16::from_le_bytes(bytes)))
            }
            b'i' => {
                let bytes = self.read_bytes::<4>()?;
                Some(AuxValue::I32(i32::from_le_bytes(bytes)))
            }
            b'I' => {
                let bytes = self.read_bytes::<4>()?;
                Some(AuxValue::U32(u32::from_le_bytes(bytes)))
            }
            b'f' => {
                let bytes = self.read_bytes::<4>()?;
                Some(AuxValue::Float(f32::from_le_bytes(bytes)))
            }
            b'd' => {
                let bytes = self.read_bytes::<8>()?;
                Some(AuxValue::Double(f64::from_le_bytes(bytes)))
            }
            b'Z' | b'H' => {
                let start = self.pos;
                while self.pos < self.data.len() && *self.data.get(self.pos)? != 0 {
                    self.pos = self.pos.checked_add(1)?;
                }
                let slice = self.data.get(start..self.pos)?;
                self.pos = self.pos.checked_add(1)?; // skip null terminator
                let v = if typ == b'Z' { AuxValue::String(slice) } else { AuxValue::Hex(slice) };
                Some(v)
            }
            // r[impl bam.record.aux_parse]
            b'B' => {
                let elem_type = *self.data.get(self.pos)?;
                let count_bytes =
                    self.data.get(self.pos.checked_add(1)?..self.pos.checked_add(5)?)?;
                let count = u32::from_le_bytes([
                    *count_bytes.first()?,
                    *count_bytes.get(1)?,
                    *count_bytes.get(2)?,
                    *count_bytes.get(3)?,
                ]) as usize;
                self.pos = self.pos.checked_add(5)?;

                let elem_size = match elem_type {
                    b'c' | b'C' => 1,
                    b's' | b'S' => 2,
                    b'i' | b'I' | b'f' => 4,
                    _ => return None,
                };

                let total = count.checked_mul(elem_size)?;
                let array_start = self.pos;
                let end = self.pos.checked_add(total)?;
                if end > self.data.len() {
                    return None;
                }
                self.pos = end;
                let array_data = self.data.get(array_start..end)?;
                let value = match elem_type {
                    b'c' => AuxValue::ArrayI8(array_data),
                    b'C' => AuxValue::ArrayU8(array_data),
                    b's' => AuxValue::ArrayI16(array_data),
                    b'S' => AuxValue::ArrayU16(array_data),
                    b'i' => AuxValue::ArrayI32(array_data),
                    b'I' => AuxValue::ArrayU32(array_data),
                    b'f' => AuxValue::ArrayFloat(array_data),
                    _ => return None,
                };
                Some(value)
            }
            _ => None,
        }
    }

    fn read_bytes<const N: usize>(&mut self) -> Option<[u8; N]> {
        let slice = self.data.get(self.pos..self.pos.checked_add(N)?)?;
        let mut arr = [0u8; N];
        arr.copy_from_slice(slice);
        self.pos = self.pos.checked_add(N)?;
        Some(arr)
    }
}

#[cfg(test)]
#[allow(clippy::arithmetic_side_effects, reason = "test arithmetic on known small values")]
mod tests {
    use super::*;

    fn build_aux(tags: &[(&[u8; 2], &[u8])]) -> Vec<u8> {
        let mut buf = Vec::new();
        for (tag, data) in tags {
            buf.extend_from_slice(&tag[..]);
            buf.extend_from_slice(data);
        }
        buf
    }

    fn z_tag(s: &[u8]) -> Vec<u8> {
        let mut v = vec![b'Z'];
        v.extend_from_slice(s);
        v.push(0);
        v
    }

    #[test]
    fn find_z_tag() {
        let value = z_tag(b"hello");
        let aux = build_aux(&[(b"RG", &value)]);
        let result = find_tag(&aux, *b"RG");
        assert_eq!(result, Some(AuxValue::String(b"hello")));
    }

    #[test]
    fn find_tag_not_present() {
        let value = z_tag(b"hello");
        let aux = build_aux(&[(b"RG", &value)]);
        assert_eq!(find_tag(&aux, *b"XY"), None);
    }

    #[test]
    fn find_tag_among_multiple() {
        let rg = z_tag(b"group1");
        let nm = [b'C', 42];
        let bc = z_tag(b"ACGT");
        let aux = build_aux(&[(b"RG", &rg), (b"NM", &nm), (b"BC", &bc)]);

        assert_eq!(find_tag(&aux, *b"RG"), Some(AuxValue::String(b"group1")));
        assert_eq!(find_tag(&aux, *b"NM"), Some(AuxValue::U8(42)));
        assert_eq!(find_tag(&aux, *b"BC"), Some(AuxValue::String(b"ACGT")));
    }

    #[test]
    fn char_type() {
        let aux = build_aux(&[(b"XA", b"AQ")]);
        assert_eq!(find_tag(&aux, *b"XA"), Some(AuxValue::Char(b'Q')));
    }

    #[test]
    fn i8_type() {
        let aux = build_aux(&[(b"X1", &[b'c', (-5i8) as u8])]);
        assert_eq!(find_tag(&aux, *b"X1"), Some(AuxValue::I8(-5)));
    }

    #[test]
    fn u8_type() {
        let aux = build_aux(&[(b"X1", &[b'C', 200])]);
        assert_eq!(find_tag(&aux, *b"X1"), Some(AuxValue::U8(200)));
    }

    #[test]
    fn i16_type() {
        let bytes = (-1234i16).to_le_bytes();
        let mut data = vec![b's'];
        data.extend_from_slice(&bytes);
        let aux = build_aux(&[(b"X1", &data)]);
        assert_eq!(find_tag(&aux, *b"X1"), Some(AuxValue::I16(-1234)));
    }

    #[test]
    fn u16_type() {
        let bytes = 60000u16.to_le_bytes();
        let mut data = vec![b'S'];
        data.extend_from_slice(&bytes);
        let aux = build_aux(&[(b"X1", &data)]);
        assert_eq!(find_tag(&aux, *b"X1"), Some(AuxValue::U16(60000)));
    }

    #[test]
    fn i32_type() {
        let bytes = (-100_000i32).to_le_bytes();
        let mut data = vec![b'i'];
        data.extend_from_slice(&bytes);
        let aux = build_aux(&[(b"X1", &data)]);
        assert_eq!(find_tag(&aux, *b"X1"), Some(AuxValue::I32(-100_000)));
    }

    #[test]
    fn u32_type() {
        let bytes = 3_000_000u32.to_le_bytes();
        let mut data = vec![b'I'];
        data.extend_from_slice(&bytes);
        let aux = build_aux(&[(b"X1", &data)]);
        assert_eq!(find_tag(&aux, *b"X1"), Some(AuxValue::U32(3_000_000)));
    }

    #[test]
    fn float_type() {
        let bytes = std::f32::consts::PI.to_le_bytes();
        let mut data = vec![b'f'];
        data.extend_from_slice(&bytes);
        let aux = build_aux(&[(b"X1", &data)]);
        assert_eq!(find_tag(&aux, *b"X1"), Some(AuxValue::Float(std::f32::consts::PI)));
    }

    #[test]
    fn double_type() {
        let bytes = std::f64::consts::E.to_le_bytes();
        let mut data = vec![b'd'];
        data.extend_from_slice(&bytes);
        let aux = build_aux(&[(b"X1", &data)]);
        assert_eq!(find_tag(&aux, *b"X1"), Some(AuxValue::Double(std::f64::consts::E)));
    }

    #[test]
    fn hex_type() {
        let mut data = vec![b'H'];
        data.extend_from_slice(b"DEADBEEF");
        data.push(0);
        let aux = build_aux(&[(b"X1", &data)]);
        assert_eq!(find_tag(&aux, *b"X1"), Some(AuxValue::Hex(b"DEADBEEF")));
    }

    #[test]
    fn empty_aux_data() {
        assert_eq!(find_tag(&[], *b"RG"), None);
    }

    #[test]
    fn truncated_header() {
        assert_eq!(find_tag(b"RG", *b"RG"), None);
    }

    #[test]
    fn truncated_value() {
        let aux = build_aux(&[(b"X1", &[b'i', 1, 2])]);
        assert_eq!(find_tag(&aux, *b"X1"), None);
    }

    #[test]
    fn array_tag_u32() {
        let mut arr = vec![b'B', b'I'];
        arr.extend_from_slice(&2u32.to_le_bytes());
        arr.extend_from_slice(&100u32.to_le_bytes());
        arr.extend_from_slice(&200u32.to_le_bytes());

        let expected: Vec<u8> =
            [100u32.to_le_bytes(), 200u32.to_le_bytes()].into_iter().flatten().collect();

        let target = z_tag(b"found");
        let aux = build_aux(&[(b"AR", &arr), (b"ZZ", &target)]);

        assert_eq!(find_tag(&aux, *b"AR"), Some(AuxValue::ArrayU32(&expected)));
        assert_eq!(find_tag(&aux, *b"ZZ"), Some(AuxValue::String(b"found")));
    }

    #[test]
    fn array_tag_i8() {
        let mut arr = vec![b'B', b'c'];
        arr.extend_from_slice(&3u32.to_le_bytes());
        arr.extend_from_slice(&[1u8, 0xfeu8, 0x7fu8]); // 1, -2, 127

        let aux = build_aux(&[(b"XB", &arr)]);
        assert_eq!(find_tag(&aux, *b"XB"), Some(AuxValue::ArrayI8(&[1, 0xfe, 0x7f])));
    }

    #[test]
    fn array_tag_float() {
        let v1 = 1.5f32.to_le_bytes();
        let v2 = (-0.5f32).to_le_bytes();
        let mut arr = vec![b'B', b'f'];
        arr.extend_from_slice(&2u32.to_le_bytes());
        arr.extend_from_slice(&v1);
        arr.extend_from_slice(&v2);

        let expected: Vec<u8> = [v1, v2].into_iter().flatten().collect();
        let aux = build_aux(&[(b"XF", &arr)]);
        assert_eq!(find_tag(&aux, *b"XF"), Some(AuxValue::ArrayFloat(&expected)));
    }

    #[test]
    fn iter_all_tags() {
        let rg = z_tag(b"grp");
        let nm = [b'C', 5];
        let aux = build_aux(&[(b"RG", &rg), (b"NM", &nm)]);

        let tags: Vec<_> = iter_tags(&aux).collect();
        assert_eq!(tags.len(), 2);
        assert_eq!(tags.first(), Some(&(*b"RG", AuxValue::String(b"grp" as &[u8]))));
        assert_eq!(tags.get(1), Some(&(*b"NM", AuxValue::U8(5))));
    }

    #[test]
    fn tag_at_end_of_data() {
        let aux = build_aux(&[(b"XA", &[b'C', 77])]);
        assert_eq!(find_tag(&aux, *b"XA"), Some(AuxValue::U8(77)));
    }

    #[test]
    fn as_str_on_string() {
        let v = AuxValue::String(b"test");
        assert_eq!(v.as_str(), Some(b"test" as &[u8]));
    }

    #[test]
    fn as_str_on_non_string() {
        let v = AuxValue::U8(42);
        assert_eq!(v.as_str(), None);
    }

    #[test]
    fn as_i64_conversions() {
        assert_eq!(AuxValue::I8(-1).as_i64(), Some(-1));
        assert_eq!(AuxValue::U8(255).as_i64(), Some(255));
        assert_eq!(AuxValue::I16(-1000).as_i64(), Some(-1000));
        assert_eq!(AuxValue::U16(65535).as_i64(), Some(65535));
        assert_eq!(AuxValue::I32(-100_000).as_i64(), Some(-100_000));
        assert_eq!(AuxValue::U32(4_000_000_000).as_i64(), Some(4_000_000_000));
        assert_eq!(AuxValue::Float(1.0).as_i64(), None);
        assert_eq!(AuxValue::String(b"x").as_i64(), None);
    }

    #[test]
    fn unknown_type_code_stops_iteration() {
        // Build: valid NM:C:5, then garbage "tag" with type 0x00, then valid RG:Z:grp
        let valid_after = z_tag(b"grp");
        let mut raw = Vec::new();
        // NM:C:5
        raw.extend_from_slice(b"NM");
        raw.push(b'C');
        raw.push(5);
        // Garbage: tag XX with unknown type 0x00
        raw.extend_from_slice(b"XX");
        raw.push(0x00);
        // RG:Z:grp
        raw.extend_from_slice(b"RG");
        raw.extend_from_slice(&valid_after);

        // Only NM should be found; RG is unreachable because unknown type code
        // stops iteration.
        assert_eq!(find_tag(&raw, *b"NM"), Some(AuxValue::U8(5)));
        assert_eq!(find_tag(&raw, *b"RG"), None);
    }

    #[test]
    fn array_with_unknown_element_type_stops_iteration() {
        let mut arr = vec![b'B', b'x']; // invalid element type
        arr.extend_from_slice(&1u32.to_le_bytes());
        arr.push(0);

        let target = z_tag(b"after");
        let aux = build_aux(&[(b"AR", &arr), (b"ZZ", &target)]);
        // Unknown B-array element type stops iteration; subsequent tags unreachable
        assert_eq!(find_tag(&aux, *b"ZZ"), None);
    }
}

// r[verify bam.record.aux_from_aux_value]
// r[verify bam.record.aux_get]
// r[verify bam.record.aux_widening]
#[cfg(test)]
mod prop_tests {
    use super::*;
    use proptest::prelude::*;

    /// Generate a valid 2-byte ASCII tag name.
    fn tag_name() -> impl Strategy<Value = [u8; 2]> {
        (prop::char::range('A', 'Z'), prop::char::range('A', 'Z'))
            .prop_map(|(a, b)| [a as u8, b as u8])
    }

    /// Encode an i64 into the smallest BAM integer type and return the raw type+value bytes.
    fn encode_int(value: i64) -> Vec<u8> {
        if value >= 0 {
            #[allow(clippy::cast_sign_loss, reason = "validated non-negative")]
            let v = value as u64;
            #[allow(clippy::cast_possible_truncation, reason = "validated max size")]
            if v <= u64::from(u8::MAX) {
                vec![b'C', v as u8]
            } else if v <= u64::from(u16::MAX) {
                let mut buf = vec![b'S'];
                buf.extend_from_slice(&(v as u16).to_le_bytes());
                buf
            } else {
                let mut buf = vec![b'I'];
                buf.extend_from_slice(&(v as u32).to_le_bytes());
                buf
            }
        } else {
            #[allow(clippy::cast_possible_truncation, reason = "validated max size")]
            if value >= i64::from(i8::MIN) {
                vec![b'c', value as u8]
            } else if value >= i64::from(i16::MIN) {
                let mut buf = vec![b's'];
                buf.extend_from_slice(&(value as i16).to_le_bytes());
                buf
            } else {
                let mut buf = vec![b'i'];
                buf.extend_from_slice(&(value as i32).to_le_bytes());
                buf
            }
        }
    }

    /// Build raw aux bytes from a list of (tag, raw_type_value) pairs.
    fn build_aux(tags: &[([u8; 2], &[u8])]) -> Vec<u8> {
        let mut buf = Vec::new();
        for (tag, data) in tags {
            buf.extend_from_slice(tag);
            buf.extend_from_slice(data);
        }
        buf
    }

    proptest! {
        // r[verify bam.record.aux_get]
        #[test]
        fn aux_get_roundtrips_ints(
            tag in tag_name(),
            value in i32::MIN as i64..=u32::MAX as i64,
        ) {
            let raw = encode_int(value);
            let aux_bytes = build_aux(&[(tag, &raw)]);
            let aux = Aux::new(&aux_bytes);
            let got: i64 = aux.get(&tag).unwrap();
            assert_eq!(got, value);
        }

        // r[verify bam.record.aux_widening]
        #[test]
        fn aux_get_widens_u8_to_u64(
            tag in tag_name(),
            value in prop::num::u8::ANY,
        ) {
            let raw = [b'C', value];
            let aux_bytes = build_aux(&[(tag, &raw)]);
            let aux = Aux::new(&aux_bytes);
            let got: u64 = aux.get(&tag).unwrap();
            assert_eq!(got, u64::from(value));
        }

        #[test]
        fn aux_get_widens_u16_to_u64(
            tag in tag_name(),
            value in prop::num::u16::ANY,
        ) {
            let mut raw = vec![b'S'];
            raw.extend_from_slice(&value.to_le_bytes());
            let aux_bytes = build_aux(&[(tag, &raw)]);
            let aux = Aux::new(&aux_bytes);
            let got: u64 = aux.get(&tag).unwrap();
            assert_eq!(got, u64::from(value));
        }

        // r[verify bam.record.aux_from_aux_value]
        #[test]
        fn aux_get_string_roundtrip(
            tag in tag_name(),
            value in "[a-zA-Z0-9]{1,20}",
        ) {
            let mut raw = vec![b'Z'];
            raw.extend_from_slice(value.as_bytes());
            raw.push(0);
            let aux_bytes = build_aux(&[(tag, &raw)]);
            let aux = Aux::new(&aux_bytes);

            let got: &str = aux.get(&tag).unwrap();
            assert_eq!(got, value);

            let owned: String = aux.get(&tag).unwrap();
            assert_eq!(owned, value);

            let smol: SmolStr = aux.get(&tag).unwrap();
            assert_eq!(smol, value);
        }

        #[test]
        fn tag_not_found_is_error(
            tag in tag_name(),
            other in tag_name(),
            value in prop::num::u8::ANY,
        ) {
            prop_assume!(tag != other);
            let raw = [b'C', value];
            let aux_bytes = build_aux(&[(tag, &raw)]);
            let aux = Aux::new(&aux_bytes);
            let err = aux.get::<u8>(&other).unwrap_err();
            assert!(matches!(err, GetAuxError::TagNotFound { .. }));
        }

        #[test]
        fn type_mismatch_is_error(
            tag in tag_name(),
            value in prop::num::u8::ANY,
        ) {
            // Tag is U8 (C type), request as &str (Z type)
            let raw = [b'C', value];
            let aux_bytes = build_aux(&[(tag, &raw)]);
            let aux = Aux::new(&aux_bytes);
            let err = aux.get::<&str>(&tag).unwrap_err();
            assert!(matches!(err, GetAuxError::TypeMismatch { .. }));
        }

        #[test]
        fn invalid_tag_name_length_is_error(
            name in prop::collection::vec(prop::num::u8::ANY, 0..=1)
                .prop_union(prop::collection::vec(prop::num::u8::ANY, 3..=10)),
        ) {
            prop_assume!(name.len() != 2);
            let aux = Aux::new(&[]);
            let err = aux.get::<u8>(&name).unwrap_err();
            assert!(matches!(err, GetAuxError::InvalidTagName { .. }));
        }
    }
}
