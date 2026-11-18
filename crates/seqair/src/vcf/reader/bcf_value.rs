//! BCF typed values ([VCF45] §6.3.3): a type descriptor byte (element count
//! in the high nibble, 15 meaning the count follows as a typed integer; type
//! code in the low nibble), then `count` elements of the type's width.
//!
//! Every count and length here comes from the file, so all of it is checked
//! arithmetic and `.get()` slicing; a value that runs past its bytes is
//! [`BcfValueError::Truncated`] (`r[io.fuzz.no_panic]`, `r[io.checked_arithmetic]`).

use core::fmt;
use core::marker::PhantomData;
use std::str::Utf8Error;

use super::values::GtAllele;
use crate::io::text_scan::find_byte;
use crate::vcf::bcf_encoding::{
    BCF_BT_CHAR, BCF_BT_FLOAT, BCF_BT_INT8, BCF_BT_INT16, BCF_BT_INT32, BCF_BT_NULL,
    FLOAT_END_OF_VECTOR, FLOAT_MISSING, INT8_END_OF_VECTOR, INT8_MISSING, INT16_END_OF_VECTOR,
    INT16_MISSING, INT32_END_OF_VECTOR, INT32_MISSING,
};

/// Malformed BCF typed values.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum BcfValueError {
    #[error("BCF typed value truncated: needs {needed} bytes, {available} available")]
    Truncated { needed: usize, available: usize },

    #[error("unknown BCF type code {code}")]
    UnknownType { code: u8 },

    #[error("expected a typed scalar integer, found type descriptor {descriptor:#04x}")]
    NotAScalarInt { descriptor: u8 },

    #[error("a typed scalar integer holds the missing or end-of-vector value")]
    MissingScalar,

    #[error("negative BCF vector length {count}")]
    NegativeCount { count: i32 },

    #[error("BCF vector of {count} {ty} values for {samples} samples overflows usize")]
    LengthOverflow { count: usize, ty: BcfType, samples: usize },

    #[error("expected integer values, found BCF type {found}")]
    NotInteger { found: BcfType },

    #[error("expected float values, found BCF type {found}")]
    NotFloat { found: BcfType },

    #[error("expected a character vector, found BCF type {found}")]
    NotChar { found: BcfType },

    #[error("BCF string is not UTF-8")]
    Utf8(#[from] Utf8Error),

    #[error("BCF genotype value {value} is negative or names an allele above 65535")]
    InvalidGenotype { value: i32 },
}

/// The type of a BCF typed value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BcfType {
    /// Type code 0: no elements, whatever the count says (a flag, or `.`).
    Missing,
    Int8,
    Int16,
    Int32,
    Float,
    Char,
}

impl BcfType {
    fn from_code(code: u8) -> Result<Self, BcfValueError> {
        Ok(match code {
            BCF_BT_NULL => Self::Missing,
            BCF_BT_INT8 => Self::Int8,
            BCF_BT_INT16 => Self::Int16,
            BCF_BT_INT32 => Self::Int32,
            BCF_BT_FLOAT => Self::Float,
            BCF_BT_CHAR => Self::Char,
            _ => return Err(BcfValueError::UnknownType { code }),
        })
    }

    /// Bytes per element.
    pub(crate) const fn width(self) -> usize {
        match self {
            Self::Missing => 0,
            Self::Int8 | Self::Char => 1,
            Self::Int16 => 2,
            Self::Int32 | Self::Float => 4,
        }
    }
}

impl fmt::Display for BcfType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Missing => "missing",
            Self::Int8 => "int8",
            Self::Int16 => "int16",
            Self::Int32 => "int32",
            Self::Float => "float",
            Self::Char => "char",
        })
    }
}

/// `bytes.get(start..start + len)`, or how many bytes it would have needed.
fn take(bytes: &[u8], start: usize, len: usize) -> Result<&[u8], BcfValueError> {
    let end = start
        .checked_add(len)
        .ok_or(BcfValueError::Truncated { needed: usize::MAX, available: bytes.len() })?;
    bytes.get(start..end).ok_or(BcfValueError::Truncated { needed: end, available: bytes.len() })
}

/// The type, element count and length of the descriptor at the start of `bytes`.
fn descriptor(bytes: &[u8]) -> Result<(BcfType, usize, usize), BcfValueError> {
    let &byte = bytes.first().ok_or(BcfValueError::Truncated { needed: 1, available: 0 })?;
    let ty = BcfType::from_code(byte & 0x0F)?;
    match byte >> 4 {
        15 => {
            let rest = bytes.get(1..).unwrap_or_default();
            let (count, used) = read_typed_int(rest).map_err(|e| shift_truncation(e, 1))?;
            let count =
                usize::try_from(count).map_err(|_| BcfValueError::NegativeCount { count })?;
            Ok((ty, count, used.saturating_add(1)))
        }
        n => Ok((ty, usize::from(n), 1)),
    }
}

/// A truncation inside a sub-slice starting at `by`, restated for the whole.
fn shift_truncation(e: BcfValueError, by: usize) -> BcfValueError {
    match e {
        BcfValueError::Truncated { needed, available } => BcfValueError::Truncated {
            needed: needed.saturating_add(by),
            available: available.saturating_add(by),
        },
        e => e,
    }
}

/// The typed scalar integer at the start of `bytes` (a dictionary key, an
/// overflow count) and its encoded length. Missing and end-of-vector are errors:
/// neither names anything.
pub(crate) fn read_typed_int(bytes: &[u8]) -> Result<(i32, usize), BcfValueError> {
    let &descriptor = bytes.first().ok_or(BcfValueError::Truncated { needed: 1, available: 0 })?;
    let ty = BcfType::from_code(descriptor & 0x0F)?;
    if descriptor >> 4 != 1 || !matches!(ty, BcfType::Int8 | BcfType::Int16 | BcfType::Int32) {
        return Err(BcfValueError::NotAScalarInt { descriptor });
    }
    let element = take(bytes, 1, ty.width())?;
    match i32::decode(ty, element) {
        Slot::Value(v) => Ok((v, ty.width().saturating_add(1))),
        Slot::Missing | Slot::End => Err(BcfValueError::MissingScalar),
    }
}

/// The encoded length of the typed value at the start of `bytes`.
pub(crate) fn skip_typed_value(bytes: &[u8]) -> Result<usize, BcfValueError> {
    TypedValue::parse(bytes).map(|(_, len)| len)
}

/// One typed value: its type, element count and the bytes of its elements.
#[derive(Debug, Clone, Copy)]
pub(crate) struct TypedValue<'a> {
    ty: BcfType,
    count: usize,
    payload: &'a [u8],
}

impl<'a> TypedValue<'a> {
    // r[impl vcf_values.bcf_sentinels]
    /// The typed value at the start of `bytes`, and its encoded length
    /// (descriptor and elements), so the caller can step to the next one.
    pub(crate) fn parse(bytes: &'a [u8]) -> Result<(Self, usize), BcfValueError> {
        let (ty, count, header) = descriptor(bytes)?;
        let len = payload_len(ty, count, 1)?;
        let payload = take(bytes, header, len)?;
        Ok((Self { ty, count, payload }, header.saturating_add(len)))
    }

    #[cfg(test)]
    pub(crate) fn ty(&self) -> BcfType {
        self.ty
    }

    /// The payload's length in bytes (the value after its descriptor).
    pub(crate) fn stored_len(&self) -> usize {
        self.payload.len()
    }

    /// The element count the descriptor declares (for [`BcfType::Missing`],
    /// with no bytes behind it).
    pub(crate) fn count(&self) -> usize {
        self.count
    }

    /// The elements as integers, at their stored width.
    pub(crate) fn ints(self) -> Result<BcfVec<'a, i32>, BcfValueError> {
        BcfVec::new(self)
    }

    /// The elements as floats.
    pub(crate) fn floats(self) -> Result<BcfVec<'a, f32>, BcfValueError> {
        BcfVec::new(self)
    }

    // r[impl vcf_values.gt+2]
    /// One sample's GT vector. `first_phase_is_read` is false for VCF 4.3 and
    /// earlier, whose first allele's phase bit is treated as unphased
    /// (`r[vcf_record.gt_first_phase+3]`).
    pub(crate) fn gt(self, first_phase_is_read: bool) -> Result<BcfGt<'a>, BcfValueError> {
        let alleles = self.ints()?;
        for value in alleles.clone().flatten() {
            gt_allele(value)?;
        }
        Ok(BcfGt { alleles, first: true, first_phase_is_read })
    }

    // r[impl vcf_values.text_strings]
    /// A character vector as a string, cut at its first NUL (padding).
    pub(crate) fn str(self) -> Result<&'a str, BcfValueError> {
        match self.ty {
            BcfType::Char => {
                let end = find_byte(self.payload, 0).unwrap_or(self.payload.len());
                Ok(std::str::from_utf8(self.payload.get(..end).unwrap_or_default())?)
            }
            BcfType::Missing => Ok(""),
            found => Err(BcfValueError::NotChar { found }),
        }
    }
}

fn payload_len(ty: BcfType, count: usize, samples: usize) -> Result<usize, BcfValueError> {
    count
        .checked_mul(ty.width())
        .and_then(|n| n.checked_mul(samples))
        .ok_or(BcfValueError::LengthOverflow { count, ty, samples })
}

/// A FORMAT field's values: one typed descriptor, then the same number of
/// elements for each sample, sample after sample.
#[derive(Debug, Clone, Copy)]
pub(crate) struct FormatColumn<'a> {
    ty: BcfType,
    /// Length of the type descriptor (with its overflow count).
    header: usize,
    per_sample: usize,
    n_samples: usize,
    stride: usize,
    payload: &'a [u8],
}

impl<'a> FormatColumn<'a> {
    /// The values of `n_samples` samples at the start of `bytes` (the bytes
    /// after the field's key), and their encoded length.
    pub(crate) fn parse(bytes: &'a [u8], n_samples: usize) -> Result<(Self, usize), BcfValueError> {
        let (ty, per_sample, header) = descriptor(bytes)?;
        let len = payload_len(ty, per_sample, n_samples)?;
        let payload = take(bytes, header, len)?;
        let stride = per_sample.saturating_mul(ty.width());
        let column = Self { ty, header, per_sample, n_samples, stride, payload };
        Ok((column, header.saturating_add(len)))
    }

    /// The length of the type descriptor at the start of the column's bytes.
    pub(crate) fn header_len(&self) -> usize {
        self.header
    }

    /// The bytes of sample `i`'s elements, `None` past the last sample.
    pub(crate) fn sample_bytes(&self, i: usize) -> Option<&'a [u8]> {
        if i >= self.n_samples {
            return None;
        }
        let start = i.checked_mul(self.stride)?;
        self.payload.get(start..start.checked_add(self.stride)?)
    }

    /// Sample `i`'s values, `None` past the last sample.
    pub(crate) fn sample(&self, i: usize) -> Option<TypedValue<'a>> {
        let payload = self.sample_bytes(i)?;
        Some(TypedValue { ty: self.ty, count: self.per_sample, payload })
    }
}

/// One stored element.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum Slot<T> {
    Value(T),
    Missing,
    End,
}

/// A number type BCF stores (`i32` at any int width, `f32`).
pub(crate) trait BcfNumber: Copy + fmt::Debug {
    fn check(ty: BcfType) -> Result<(), BcfValueError>;
    /// The element in `bytes` (at least `ty.width()` of them), [`Slot::End`]
    /// for too few.
    fn decode(ty: BcfType, bytes: &[u8]) -> Slot<Self>;
}

impl BcfNumber for i32 {
    fn check(ty: BcfType) -> Result<(), BcfValueError> {
        match ty {
            BcfType::Missing | BcfType::Int8 | BcfType::Int16 | BcfType::Int32 => Ok(()),
            found => Err(BcfValueError::NotInteger { found }),
        }
    }

    #[inline]
    fn decode(ty: BcfType, bytes: &[u8]) -> Slot<Self> {
        match ty {
            BcfType::Int8 => match bytes.first() {
                Some(&INT8_MISSING) => Slot::Missing,
                Some(&INT8_END_OF_VECTOR) => Slot::End,
                Some(&b) => Slot::Value(i32::from(i8::from_le_bytes([b]))),
                None => Slot::End,
            },
            BcfType::Int16 => match bytes.first_chunk::<2>().map(|b| u16::from_le_bytes(*b)) {
                Some(INT16_MISSING) => Slot::Missing,
                Some(INT16_END_OF_VECTOR) | None => Slot::End,
                Some(bits) => Slot::Value(i32::from(i16::from_le_bytes(bits.to_le_bytes()))),
            },
            BcfType::Int32 => match bytes.first_chunk::<4>().map(|b| u32::from_le_bytes(*b)) {
                Some(INT32_MISSING) => Slot::Missing,
                Some(INT32_END_OF_VECTOR) | None => Slot::End,
                Some(bits) => Slot::Value(i32::from_le_bytes(bits.to_le_bytes())),
            },
            BcfType::Missing | BcfType::Float | BcfType::Char => Slot::End,
        }
    }
}

impl BcfNumber for f32 {
    fn check(ty: BcfType) -> Result<(), BcfValueError> {
        match ty {
            BcfType::Missing | BcfType::Float => Ok(()),
            found => Err(BcfValueError::NotFloat { found }),
        }
    }

    /// By bit pattern: only the two sentinels are special, any other NaN is a value.
    #[inline]
    fn decode(ty: BcfType, bytes: &[u8]) -> Slot<Self> {
        if ty != BcfType::Float {
            return Slot::End;
        }
        match bytes.first_chunk::<4>().map(|b| u32::from_le_bytes(*b)) {
            Some(FLOAT_MISSING) => Slot::Missing,
            Some(FLOAT_END_OF_VECTOR) | None => Slot::End,
            Some(bits) => Slot::Value(f32::from_bits(bits)),
        }
    }
}

// r[impl vcf_values.bcf_sentinels]
/// A BCF number vector at its stored width, widened per element: missing is
/// `None`, and the first end-of-vector ends it (nothing after it is read).
#[derive(Debug, Clone)]
pub(crate) struct BcfVec<'a, T> {
    /// The elements not yet yielded, cut at the first end-of-vector.
    rest: &'a [u8],
    /// A numeric type, never [`BcfType::Missing`], so its width is nonzero.
    ty: BcfType,
    _t: PhantomData<T>,
}

impl<'a, T: BcfNumber> BcfVec<'a, T> {
    fn new(value: TypedValue<'a>) -> Result<Self, BcfValueError> {
        T::check(value.ty)?;
        if value.ty == BcfType::Missing {
            return Ok(Self { rest: &[], ty: BcfType::Int32, _t: PhantomData });
        }
        let width = value.ty.width();
        let len = value
            .payload
            .chunks_exact(width)
            .position(|e| matches!(T::decode(value.ty, e), Slot::End))
            .unwrap_or(value.count);
        let rest = value.payload.get(..len.saturating_mul(width)).unwrap_or(value.payload);
        Ok(Self { rest, ty: value.ty, _t: PhantomData })
    }

    /// Element `i` of those not yet yielded; `None` past the end.
    #[expect(
        clippy::option_option,
        reason = "in range, then the element as the iterator yields it"
    )]
    pub(crate) fn get(&self, i: usize) -> Option<Option<T>> {
        let width = self.ty.width();
        let start = i.checked_mul(width)?;
        match T::decode(self.ty, self.rest.get(start..start.checked_add(width)?)?) {
            Slot::Value(v) => Some(Some(v)),
            Slot::Missing => Some(None),
            Slot::End => None,
        }
    }
}

impl<T: BcfNumber> Iterator for BcfVec<'_, T> {
    type Item = Option<T>;

    #[inline]
    fn next(&mut self) -> Option<Option<T>> {
        let (element, rest) = self.rest.split_at_checked(self.ty.width())?;
        self.rest = rest;
        match T::decode(self.ty, element) {
            Slot::Value(v) => Some(Some(v)),
            Slot::Missing => Some(None),
            // `new` cut `rest` at the first end-of-vector.
            Slot::End => None,
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = self.rest.len().checked_div(self.ty.width()).unwrap_or(0);
        (len, Some(len))
    }
}

impl<T: BcfNumber> ExactSizeIterator for BcfVec<'_, T> {}

/// The allele of a GT element `(allele + 1) << 1 | phased`: `None` for `.`.
fn gt_allele(value: i32) -> Result<Option<u16>, BcfValueError> {
    if value < 0 {
        return Err(BcfValueError::InvalidGenotype { value });
    }
    match value >> 1 {
        0 => Ok(None),
        n => u16::try_from(n.saturating_sub(1))
            .map(Some)
            .map_err(|_| BcfValueError::InvalidGenotype { value }),
    }
}

// r[impl vcf_values.gt+2]
/// One sample's GT alleles, each with its phase bit; a ploidy below the
/// field's is padded with end-of-vector, which ends it. The int missing value
/// (htslib's fill for a sample without GT) is a missing, unphased allele.
#[derive(Debug, Clone)]
pub(crate) struct BcfGt<'a> {
    alleles: BcfVec<'a, i32>,
    first: bool,
    first_phase_is_read: bool,
}

impl Iterator for BcfGt<'_> {
    type Item = GtAllele;

    fn next(&mut self) -> Option<GtAllele> {
        let element = self.alleles.next()?;
        let first = core::mem::replace(&mut self.first, false);
        let Some(value) = element else {
            return Some(GtAllele { allele: None, phased: false });
        };
        let allele = gt_allele(value);
        debug_assert!(allele.is_ok(), "gt() checked every element");
        let phased = value & 1 == 1 && (self.first_phase_is_read || !first);
        Some(GtAllele { allele: allele.unwrap_or(None), phased })
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.alleles.size_hint()
    }
}

impl ExactSizeIterator for BcfGt<'_> {}

#[cfg(test)]
mod tests {
    #![allow(
        clippy::arithmetic_side_effects,
        clippy::cast_possible_truncation,
        reason = "test code with small known values"
    )]

    use super::*;
    use crate::vcf::bcf_encoding::{
        encode_int_eov, encode_int_value_or_missing, encode_type_byte, encode_typed_int_key,
        smallest_int_type_iter,
    };
    use hegel::prelude::*;

    fn ints(bytes: &[u8]) -> Vec<Option<i32>> {
        TypedValue::parse(bytes).unwrap().0.ints().unwrap().collect()
    }

    fn float_bits(bytes: &[u8]) -> Vec<Option<u32>> {
        let v = TypedValue::parse(bytes).unwrap().0.floats().unwrap();
        v.map(|f| f.map(f32::to_bits)).collect()
    }

    /// Integers BCF can store: everything above the int32 reserved band.
    fn storable() -> impl PrintableGenerator<i32> {
        gs::integers::<i32>().min_value(i32::MIN + 8)
    }

    /// `values` then `pad` end-of-vector values, at the smallest width that
    /// fits them or a wider one, then `trailing` bytes that belong to the next value.
    fn encode_ints(values: &[Option<i32>], pad: usize, min_ty: u8, trailing: &[u8]) -> Vec<u8> {
        let ty = smallest_int_type_iter(values.iter().flatten().copied()).max(min_ty);
        let mut buf = Vec::new();
        encode_type_byte(&mut buf, values.len() + pad, ty);
        for &v in values {
            encode_int_value_or_missing(&mut buf, v, ty);
        }
        for _ in 0..pad {
            encode_int_eov(&mut buf, ty);
        }
        buf.extend_from_slice(trailing);
        buf
    }

    // r[verify vcf_values.bcf_sentinels]
    /// The writer's encoding read back: every width, missing elements,
    /// end-of-vector padding, and the overflow count (16 elements and up).
    #[hegel::test]
    fn int_vectors_read_back_as_written(tc: TestCase) {
        let values = tc.draw(gs::vecs(gs::optional(storable())).max_size(40));
        let pad = tc.draw(gs::integers::<usize>().max_value(20));
        let min_ty = tc.draw(gs::sampled_from(vec![BCF_BT_INT8, BCF_BT_INT16, BCF_BT_INT32]));
        let trailing = tc.draw(gs::binary().max_size(4));
        let buf = encode_ints(&values, pad, min_ty, &trailing);

        let (value, len) = TypedValue::parse(&buf).unwrap();
        assert_eq!(len, buf.len() - trailing.len());
        assert_eq!(skip_typed_value(&buf).unwrap(), len);
        assert_eq!(value.count(), values.len() + pad);
        let v = value.ints().unwrap();
        assert_eq!(v.len(), values.len());
        for (i, &expected) in values.iter().enumerate() {
            assert_eq!(v.get(i), Some(expected), "element {i}");
        }
        assert_eq!(v.get(values.len()), None);
        assert_eq!(v.collect::<Vec<_>>(), values);
    }

    // r[verify vcf_values.bcf_sentinels]
    /// Any float bit pattern but the reserved band reads back bit for bit,
    /// NaNs included; the two sentinels become `None` and the end.
    #[hegel::test]
    fn float_vectors_read_back_bit_for_bit(tc: TestCase) {
        let bits = gs::integers::<u32>().filter(|b| !(FLOAT_MISSING..=0x7F80_0007).contains(b));
        let values = tc.draw(gs::vecs(gs::optional(bits)).max_size(40));
        let pad = tc.draw(gs::integers::<usize>().max_value(3));
        let mut buf = Vec::new();
        encode_type_byte(&mut buf, values.len() + pad, BCF_BT_FLOAT);
        for v in &values {
            buf.extend_from_slice(&v.unwrap_or(FLOAT_MISSING).to_le_bytes());
        }
        for _ in 0..pad {
            buf.extend_from_slice(&FLOAT_END_OF_VECTOR.to_le_bytes());
        }
        assert_eq!(float_bits(&buf), values);
    }

    // r[verify vcf_values.text_strings]
    #[hegel::test]
    fn strings_read_back_cut_at_the_first_nul(tc: TestCase) {
        let s = tc.draw(gs::text().max_size(40)).replace('\0', "");
        let padding = tc.draw(gs::vecs(gs::integers::<u8>()).max_size(5));
        let mut bytes = s.clone().into_bytes();
        if !padding.is_empty() {
            bytes.push(0);
            bytes.extend_from_slice(&padding);
        }
        let mut buf = Vec::new();
        encode_type_byte(&mut buf, bytes.len(), BCF_BT_CHAR);
        buf.extend_from_slice(&bytes);
        assert_eq!(TypedValue::parse(&buf).unwrap().0.str().unwrap(), s);
    }

    // r[verify vcf_values.bcf_sentinels]
    /// Per-sample vectors of different lengths, padded to the longest with
    /// end-of-vector, as `[VCF45]` §6.3.3 "Vectors of mixed length" lays them out.
    #[hegel::test]
    fn format_columns_read_back_per_sample(tc: TestCase) {
        let samples = tc.draw(gs::vecs(gs::vecs(gs::optional(storable())).max_size(5)).max_size(6));
        let per_sample = samples.iter().map(Vec::len).max().unwrap_or(0);
        let ty = smallest_int_type_iter(samples.iter().flatten().flatten().copied());
        let mut buf = Vec::new();
        encode_typed_int_key(&mut buf, 3);
        let key_len = buf.len();
        encode_type_byte(&mut buf, per_sample, ty);
        for s in &samples {
            for &v in s {
                encode_int_value_or_missing(&mut buf, v, ty);
            }
            for _ in s.len()..per_sample {
                encode_int_eov(&mut buf, ty);
            }
        }
        assert_eq!(read_typed_int(&buf).unwrap(), (3, key_len));
        let (column, len) = FormatColumn::parse(&buf[key_len..], samples.len()).unwrap();
        assert_eq!(key_len + len, buf.len());
        for (i, s) in samples.iter().enumerate() {
            let got: Vec<_> = column.sample(i).unwrap().ints().unwrap().collect();
            assert_eq!(&got, s, "sample {i}");
        }
        assert!(column.sample(samples.len()).is_none());
    }

    /// Arbitrary bytes never panic, and a parsed value never claims more
    /// bytes than it was given.
    // r[verify io.fuzz.no_panic]
    #[hegel::test]
    fn arbitrary_bytes_do_not_panic(tc: TestCase) {
        let bytes = tc.draw(gs::binary().max_size(64));
        let n_samples = tc.draw(gs::integers::<usize>().max_value(4));
        if let Ok((v, len)) = TypedValue::parse(&bytes) {
            assert!(len <= bytes.len());
            let _ = v.ints().map(Iterator::count);
            let _ = v.floats().map(Iterator::count);
            let _ = v.str();
        }
        if let Ok((c, len)) = FormatColumn::parse(&bytes, n_samples) {
            assert!(len <= bytes.len());
            for i in 0..=n_samples {
                let _ = c.sample(i).map(|s| s.ints().map(Iterator::count));
            }
        }
        let _ = read_typed_int(&bytes);
    }

    // r[verify vcf_values.bcf_sentinels]
    #[test]
    fn sentinels_at_every_width() {
        // int8: missing 0x80, end 0x81; 0x82 (future use) is the value -126.
        assert_eq!(ints(&[0x41, 0x01, 0x80, 0x82, 0x81]), [Some(1), None, Some(-126)]);
        // int16
        assert_eq!(ints(&[0x32, 0xFF, 0x7F, 0x00, 0x80, 0x01, 0x80]), [Some(32767), None]);
        // int32
        let mut buf = vec![0x33];
        for v in [i32::MAX as u32, INT32_MISSING, INT32_END_OF_VECTOR] {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        assert_eq!(ints(&buf), [Some(i32::MAX), None]);
        // Nothing after the first end-of-vector is read, not even a value.
        assert_eq!(ints(&[0x31, 0x81, 0x05, 0x06]), []);
        // float: quiet NaN and infinities are values, the two sentinels are not.
        let mut buf = vec![0x55];
        for v in [f32::NAN.to_bits(), f32::INFINITY.to_bits(), FLOAT_MISSING, 0x7F80_0003] {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        buf.extend_from_slice(&FLOAT_END_OF_VECTOR.to_le_bytes());
        assert_eq!(
            float_bits(&buf),
            [Some(f32::NAN.to_bits()), Some(f32::INFINITY.to_bits()), None, Some(0x7F80_0003)]
        );
    }

    #[test]
    fn missing_type_is_empty_whatever_its_count() {
        let (v, len) = TypedValue::parse(&[0x30, 0xAA]).unwrap();
        assert_eq!((v.ty(), v.count(), len), (BcfType::Missing, 3, 1));
        assert_eq!(v.ints().unwrap().len(), 0);
        assert_eq!(v.floats().unwrap().len(), 0);
        assert_eq!(v.str().unwrap(), "");
    }

    #[test]
    fn overflow_counts_at_every_width() {
        for count_bytes in [vec![0x11, 16], vec![0x12, 16, 0], vec![0x13, 16, 0, 0, 0]] {
            let mut buf = vec![0xF1];
            buf.extend_from_slice(&count_bytes);
            buf.extend(1..=16u8);
            let (v, len) = TypedValue::parse(&buf).unwrap();
            assert_eq!(len, buf.len());
            assert_eq!(
                v.ints().unwrap().flatten().collect::<Vec<_>>(),
                (1..=16).collect::<Vec<_>>()
            );
        }
    }

    #[test]
    fn malformed_values_are_typed_errors() {
        use BcfValueError as E;
        assert!(matches!(TypedValue::parse(&[]), Err(E::Truncated { needed: 1, available: 0 })));
        for code in [4, 6, 8, 15] {
            assert!(
                matches!(TypedValue::parse(&[0x10 | code]), Err(E::UnknownType { code: c }) if c == code)
            );
        }
        // Three int16 elements need six bytes.
        assert!(matches!(
            TypedValue::parse(&[0x32, 1, 2, 3]),
            Err(E::Truncated { needed: 7, available: 4 })
        ));
        // An overflow count that is negative, missing, a vector, a float, or cut off.
        assert!(matches!(
            TypedValue::parse(&[0xF1, 0x11, 0xFF]),
            Err(E::NegativeCount { count: -1 })
        ));
        assert!(matches!(TypedValue::parse(&[0xF1, 0x11, 0x80]), Err(E::MissingScalar)));
        assert!(matches!(
            TypedValue::parse(&[0xF1, 0x21, 1, 1]),
            Err(E::NotAScalarInt { descriptor: 0x21 })
        ));
        assert!(matches!(TypedValue::parse(&[0xF1, 0x15]), Err(E::NotAScalarInt { .. })));
        assert!(matches!(
            TypedValue::parse(&[0xF1, 0x13, 1]),
            Err(E::Truncated { needed: 6, available: 3 })
        ));
        // A huge count of wide elements overflows before anything is sliced.
        let mut huge = vec![0xF3, 0x13];
        huge.extend_from_slice(&i32::MAX.to_le_bytes());
        assert!(matches!(
            TypedValue::parse(&huge),
            Err(E::Truncated { .. } | E::LengthOverflow { .. })
        ));
        assert!(matches!(
            FormatColumn::parse(&huge, usize::MAX),
            Err(E::LengthOverflow { ty: BcfType::Int32, .. })
        ));
        // Asking for the wrong kind.
        let (int, _) = TypedValue::parse(&[0x11, 1]).unwrap();
        assert!(matches!(int.floats(), Err(E::NotFloat { found: BcfType::Int8 })));
        assert!(matches!(int.str(), Err(E::NotChar { found: BcfType::Int8 })));
        let (float, _) = TypedValue::parse(&[0x15, 0, 0, 0, 0]).unwrap();
        assert!(matches!(float.ints(), Err(E::NotInteger { found: BcfType::Float })));
        let (bad, _) = TypedValue::parse(&[0x27, 0xC3, 0x28]).unwrap();
        assert!(matches!(bad.str(), Err(E::Utf8(_))));
    }

    #[test]
    fn typed_int_keys_at_every_width() {
        for key in [0, 127, 128, 32767, 32768, 1 << 20] {
            let mut buf = Vec::new();
            encode_typed_int_key(&mut buf, key);
            assert_eq!(read_typed_int(&buf).unwrap(), (i32::try_from(key).unwrap(), buf.len()));
        }
        assert!(matches!(read_typed_int(&[0x13, 1]), Err(BcfValueError::Truncated { .. })));
    }

    // r[verify vcf_values.gt+2]
    #[test]
    fn malformed_bcf_genotypes_are_typed_errors() {
        fn gt(bytes: &[u8]) -> Result<BcfGt<'_>, BcfValueError> {
            TypedValue::parse(bytes)?.0.gt(true)
        }
        assert!(matches!(
            gt(&[0x21, 0x02, 0xFB]),
            Err(BcfValueError::InvalidGenotype { value: -5 })
        ));
        let mut too_big = vec![0x13];
        too_big.extend_from_slice(&((65537i32 << 1) | 1).to_le_bytes());
        assert!(matches!(gt(&too_big), Err(BcfValueError::InvalidGenotype { value: 131_075 })));
        assert!(matches!(gt(&[0x15, 0, 0, 0, 0]), Err(BcfValueError::NotInteger { .. })));
        // An empty GT has no alleles.
        assert_eq!(gt(&[0x01]).unwrap().len(), 0);
        assert_eq!(gt(&[0x11, 0x81]).unwrap().len(), 0);
    }
}
