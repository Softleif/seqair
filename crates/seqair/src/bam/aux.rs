//! BAM auxiliary tag parsing.
//!
//! Parses raw BAM auxiliary data bytes into typed [`AuxValue`] values.

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
}

/// Find a tag in raw BAM auxiliary data and return its typed value.
///
/// Returns `None` if the tag is not found or the data is truncated.
///
/// B-type (array) tags are not yet supported and will return `None`, even when
/// the tag is present. The bytes are skipped correctly so that subsequent tags
/// in the same record remain reachable.
pub fn find_tag<'a>(aux: &'a [u8], tag: [u8; 2]) -> Option<AuxValue<'a>> {
    for (t, value) in iter_tags(aux) {
        if t == tag {
            return Some(value);
        }
    }
    None
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
        loop {
            let data = self.data;
            if self.pos + 3 > data.len() {
                return None;
            }

            let tag = [*data.get(self.pos)?, *data.get(self.pos + 1)?];
            let typ = *data.get(self.pos + 2)?;
            self.pos += 3;

            match self.parse_value(typ) {
                Some(Some(value)) => return Some((tag, value)),
                Some(None) => continue, // skippable tag (e.g. array)
                None => return None,    // malformed data
            }
        }
    }
}

impl<'a> AuxIter<'a> {
    /// Returns `Some(Some(value))` for a successfully parsed tag,
    /// `Some(None)` for a tag that was skipped (e.g. arrays),
    /// or `None` for malformed/truncated data.
    fn parse_value(&mut self, typ: u8) -> Option<Option<AuxValue<'a>>> {
        match typ {
            b'A' => {
                let v = *self.data.get(self.pos)?;
                self.pos += 1;
                Some(Some(AuxValue::Char(v)))
            }
            b'c' => {
                let v = *self.data.get(self.pos)?;
                self.pos += 1;
                Some(Some(AuxValue::I8(v as i8)))
            }
            b'C' => {
                let v = *self.data.get(self.pos)?;
                self.pos += 1;
                Some(Some(AuxValue::U8(v)))
            }
            b's' => {
                let bytes = self.read_bytes::<2>()?;
                Some(Some(AuxValue::I16(i16::from_le_bytes(bytes))))
            }
            b'S' => {
                let bytes = self.read_bytes::<2>()?;
                Some(Some(AuxValue::U16(u16::from_le_bytes(bytes))))
            }
            b'i' => {
                let bytes = self.read_bytes::<4>()?;
                Some(Some(AuxValue::I32(i32::from_le_bytes(bytes))))
            }
            b'I' => {
                let bytes = self.read_bytes::<4>()?;
                Some(Some(AuxValue::U32(u32::from_le_bytes(bytes))))
            }
            b'f' => {
                let bytes = self.read_bytes::<4>()?;
                Some(Some(AuxValue::Float(f32::from_le_bytes(bytes))))
            }
            b'd' => {
                let bytes = self.read_bytes::<8>()?;
                Some(Some(AuxValue::Double(f64::from_le_bytes(bytes))))
            }
            b'Z' | b'H' => {
                let start = self.pos;
                while self.pos < self.data.len() && *self.data.get(self.pos)? != 0 {
                    self.pos += 1;
                }
                let slice = self.data.get(start..self.pos)?;
                self.pos += 1; // skip null terminator
                let v = if typ == b'Z' { AuxValue::String(slice) } else { AuxValue::Hex(slice) };
                Some(Some(v))
            }
            b'B' => {
                let elem_type = *self.data.get(self.pos)?;
                let count_bytes = self.data.get(self.pos + 1..self.pos + 5)?;
                let count = u32::from_le_bytes([
                    *count_bytes.first()?,
                    *count_bytes.get(1)?,
                    *count_bytes.get(2)?,
                    *count_bytes.get(3)?,
                ]) as usize;
                self.pos += 5;

                let elem_size = match elem_type {
                    b'c' | b'C' => 1,
                    b's' | b'S' => 2,
                    b'i' | b'I' | b'f' => 4,
                    _ => return None,
                };

                let total = count.checked_mul(elem_size)?;
                let end = self.pos.checked_add(total)?;
                if end > self.data.len() {
                    return None;
                }
                self.pos = end;
                Some(None) // arrays not yet represented in AuxValue
            }
            _ => None,
        }
    }

    fn read_bytes<const N: usize>(&mut self) -> Option<[u8; N]> {
        let slice = self.data.get(self.pos..self.pos + N)?;
        let mut arr = [0u8; N];
        arr.copy_from_slice(slice);
        self.pos += N;
        Some(arr)
    }
}

#[cfg(test)]
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
    fn skip_array_tag() {
        let mut arr = vec![b'B', b'I'];
        arr.extend_from_slice(&2u32.to_le_bytes());
        arr.extend_from_slice(&100u32.to_le_bytes());
        arr.extend_from_slice(&200u32.to_le_bytes());

        let target = z_tag(b"found");
        let aux = build_aux(&[(b"AR", &arr), (b"ZZ", &target)]);

        assert_eq!(find_tag(&aux, *b"ZZ"), Some(AuxValue::String(b"found")));
        assert_eq!(find_tag(&aux, *b"AR"), None);
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
    fn array_with_unknown_element_type_stops_iteration() {
        let mut arr = vec![b'B', b'x']; // invalid element type
        arr.extend_from_slice(&1u32.to_le_bytes());
        arr.push(0);

        let target = z_tag(b"after");
        let aux = build_aux(&[(b"AR", &arr), (b"ZZ", &target)]);
        assert_eq!(find_tag(&aux, *b"ZZ"), None);
    }
}
