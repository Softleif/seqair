//! VCF text values ([VCF45] §1.6): integers, floats and strings as written,
//! `.` for missing, vectors separated by `,`.

use core::marker::PhantomData;
use core::num::IntErrorKind;
use std::str::Utf8Error;

use super::values::GtAllele;
use crate::io::text_scan::find_byte;

/// Malformed VCF text values.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum TextValueError {
    #[error("not an integer (an optional sign and decimal digits)")]
    InvalidInteger,

    #[error("integer out of the 32-bit range")]
    IntegerOutOfRange,

    #[error("not a floating-point number")]
    InvalidFloat,

    #[error("string value is not UTF-8")]
    Utf8(#[from] Utf8Error),

    #[error("not a genotype (alleles as digits or `.`, separated by `/` or `|`)")]
    InvalidGenotype,

    #[error("genotype allele index above 65535")]
    GenotypeAlleleOutOfRange,
}

/// A number VCF text holds, `.` being missing.
pub(crate) trait TextNumber: Copy {
    fn parse_text(token: &[u8]) -> Result<Option<Self>, TextValueError>;
}

impl TextNumber for i32 {
    fn parse_text(token: &[u8]) -> Result<Option<Self>, TextValueError> {
        parse_i32(token)
    }
}

impl TextNumber for f32 {
    fn parse_text(token: &[u8]) -> Result<Option<Self>, TextValueError> {
        parse_f32(token)
    }
}

// r[impl vcf_values.text_numbers]
/// An optional sign and decimal digits, nothing else; `.` is missing.
pub(crate) fn parse_i32(token: &[u8]) -> Result<Option<i32>, TextValueError> {
    if token == b"." {
        return Ok(None);
    }
    // `i32::from_str` is exactly this grammar: `[+-]?[0-9]+`, no whitespace.
    let s = std::str::from_utf8(token).map_err(|_| TextValueError::InvalidInteger)?;
    s.parse::<i32>().map(Some).map_err(|e| match e.kind() {
        IntErrorKind::PosOverflow | IntErrorKind::NegOverflow => TextValueError::IntegerOutOfRange,
        _ => TextValueError::InvalidInteger,
    })
}

// r[impl vcf_values.text_numbers]
/// A float as `f32::from_str` reads it (correctly rounded), which also takes
/// VCF's `Inf`, `Infinity` and `NaN` in any case with an optional sign; `.`
/// is missing.
pub(crate) fn parse_f32(token: &[u8]) -> Result<Option<f32>, TextValueError> {
    if token == b"." {
        return Ok(None);
    }
    let s = std::str::from_utf8(token).map_err(|_| TextValueError::InvalidFloat)?;
    s.parse().map(Some).map_err(|_| TextValueError::InvalidFloat)
}

// r[impl vcf_values.text_strings]
/// A string value as stored: UTF-8 checked, percent-encoding left alone (as
/// htslib does, so VCF and BCF reads agree).
pub(crate) fn text_str(bytes: &[u8]) -> Result<&str, TextValueError> {
    Ok(std::str::from_utf8(bytes)?)
}

/// The token before the first `,` of `rest`, and what follows that comma.
#[inline]
fn split_comma(rest: &[u8]) -> (&[u8], Option<&[u8]>) {
    match find_byte(rest, b',') {
        Some(i) => (rest.get(..i).unwrap_or_default(), rest.get(i.saturating_add(1)..)),
        None => (rest, None),
    }
}

/// A `,`-separated number vector, missing elements (`.`) as `None`. Every
/// element is parsed once by [`TextVec::new`], so iterating cannot fail.
#[derive(Debug, Clone)]
pub(crate) struct TextVec<'a, T> {
    /// The elements not yet yielded; `None` once all are.
    rest: Option<&'a [u8]>,
    len: usize,
    _t: PhantomData<T>,
}

impl<'a, T: TextNumber> TextVec<'a, T> {
    /// No elements: an INFO key written without a value.
    pub(crate) fn empty() -> Self {
        Self { rest: None, len: 0, _t: PhantomData }
    }

    /// Checks every element of `value` (`""` is one empty, invalid element).
    pub(crate) fn new(value: &'a [u8]) -> Result<Self, TextValueError> {
        let mut len = 0usize;
        let mut rest = Some(value);
        while let Some(r) = rest {
            let (token, tail) = split_comma(r);
            T::parse_text(token)?;
            len = len.saturating_add(1);
            rest = tail;
        }
        Ok(Self { rest: Some(value), len, _t: PhantomData })
    }

    /// Element `i` of those not yet yielded; `None` past the end.
    #[expect(
        clippy::option_option,
        reason = "in range, then the element as the iterator yields it"
    )]
    pub(crate) fn get(&self, i: usize) -> Option<Option<T>> {
        self.clone().nth(i)
    }
}

impl<T: TextNumber> Iterator for TextVec<'_, T> {
    type Item = Option<T>;

    #[inline]
    fn next(&mut self) -> Option<Option<T>> {
        let (token, tail) = split_comma(self.rest?);
        self.rest = tail;
        self.len = self.len.saturating_sub(1);
        let v = T::parse_text(token);
        debug_assert!(v.is_ok(), "new() parsed every element");
        v.ok()
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.len, Some(self.len))
    }
}

impl<T: TextNumber> ExactSizeIterator for TextVec<'_, T> {}

/// The allele before the first `/` or `|` of `rest`, whether that separator
/// is `|`, and what follows it. Genotypes are a few bytes, so this is scalar.
#[inline]
fn split_allele(rest: &[u8]) -> (&[u8], Option<(bool, &[u8])>) {
    match rest.iter().position(|&b| b == b'/' || b == b'|') {
        Some(i) => (
            rest.get(..i).unwrap_or_default(),
            rest.get(i..).and_then(|sep| Some((*sep.first()? == b'|', sep.get(1..)?))),
        ),
        None => (rest, None),
    }
}

/// One allele: decimal digits, or `.` for missing.
fn parse_allele(token: &[u8]) -> Result<Option<u16>, TextValueError> {
    if token == b"." {
        return Ok(None);
    }
    if token.is_empty() || !token.iter().all(u8::is_ascii_digit) {
        return Err(TextValueError::InvalidGenotype);
    }
    let s = std::str::from_utf8(token).map_err(|_| TextValueError::InvalidGenotype)?;
    s.parse().map(Some).map_err(|_| TextValueError::GenotypeAlleleOutOfRange)
}

// r[impl vcf_values.gt]
/// A text GT (`0/1`, `0|1`, `.`, `./.`, `1`, `0/1|2`, VCF 4.4's `|0|1`),
/// allele by allele, checked once by [`TextGt::new`].
#[derive(Debug, Clone)]
pub(crate) struct TextGt<'a> {
    /// The alleles not yet yielded; `None` once all are.
    rest: Option<&'a [u8]>,
    len: usize,
    /// The phase of the next allele.
    phased: bool,
}

impl<'a> TextGt<'a> {
    // r[impl vcf_record.gt_first_phase+3]
    /// `gt` checked and its first allele's phase settled: with
    /// `first_phase_is_read` (VCF 4.4 and later) a leading `|` or `/` gives
    /// it; without one, htslib's rule infers it — a haploid allele is phased
    /// unless it is `.`, a longer genotype is phased iff every separator is
    /// `|` — so text and htslib-written BCF agree. Before 4.4 it is unphased.
    pub(crate) fn new(gt: &'a [u8], first_phase_is_read: bool) -> Result<Self, TextValueError> {
        let (prefix, body) = match gt.split_first() {
            Some((&b'|', body)) => (Some(true), body),
            Some((&b'/', body)) => (Some(false), body),
            _ => (None, gt),
        };
        let mut len = 0usize;
        let mut all_phased = true;
        let mut first_missing = false;
        let mut rest = Some(body);
        while let Some(r) = rest {
            let (token, sep) = split_allele(r);
            let allele = parse_allele(token)?;
            if len == 0 {
                first_missing = allele.is_none();
            }
            len = len.saturating_add(1);
            rest = sep.map(|(phased, tail)| {
                all_phased &= phased;
                tail
            });
        }
        let phased = first_phase_is_read
            && prefix.unwrap_or(if len == 1 { !first_missing } else { all_phased });
        Ok(Self { rest: Some(body), len, phased })
    }
}

impl Iterator for TextGt<'_> {
    type Item = GtAllele;

    fn next(&mut self) -> Option<GtAllele> {
        let (token, sep) = split_allele(self.rest?);
        let phased = self.phased;
        (self.phased, self.rest) = sep.map_or((false, None), |(p, tail)| (p, Some(tail)));
        self.len = self.len.saturating_sub(1);
        let allele = parse_allele(token);
        debug_assert!(allele.is_ok(), "new() checked every allele");
        Some(GtAllele { allele: allele.unwrap_or(None), phased })
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.len, Some(self.len))
    }
}

impl ExactSizeIterator for TextGt<'_> {}

#[cfg(test)]
mod tests {
    use super::*;
    use hegel::prelude::*;

    // r[verify vcf_values.text_numbers]
    /// Any `i32` std prints, with an optional `+` and leading zeros, reads back.
    #[hegel::test]
    fn integers_read_back_as_printed(tc: TestCase) {
        let v = tc.draw(gs::integers::<i32>());
        let zeros = "0".repeat(tc.draw(gs::integers::<usize>().max_value(3)));
        let plus = if v >= 0 && tc.draw(gs::booleans()) { "+" } else { "" };
        let digits = v.unsigned_abs().to_string();
        let sign = if v < 0 { "-" } else { plus };
        let text = format!("{sign}{zeros}{digits}");
        assert_eq!(parse_i32(text.as_bytes()).unwrap(), Some(v), "{text}");
    }

    // r[verify vcf_values.text_numbers]
    /// Any `f32` std prints, plainly or in exponent form, reads back bit for
    /// bit (NaN as some NaN: the text has no payload).
    #[hegel::test]
    fn floats_read_back_as_printed(tc: TestCase) {
        let v = f32::from_bits(tc.draw(gs::integers::<u32>()));
        let text = if tc.draw(gs::booleans()) { format!("{v}") } else { format!("{v:e}") };
        let got = parse_f32(text.as_bytes()).unwrap().unwrap();
        if v.is_nan() {
            assert!(got.is_nan(), "{text}");
        } else {
            assert_eq!(got.to_bits(), v.to_bits(), "{text}");
        }
    }

    // r[verify vcf_values.text_numbers]
    #[test]
    fn integers_are_strict() {
        assert_eq!(parse_i32(b".").unwrap(), None);
        assert_eq!(parse_i32(b"-0").unwrap(), Some(0));
        assert_eq!(parse_i32(b"-2147483648").unwrap(), Some(i32::MIN));
        for bad in
            [&b""[..], b"-", b"+", b" 1", b"1 ", b"1.0", b"1e3", b"0x1", b"..", b"\xff", b"1,2"]
        {
            assert!(matches!(parse_i32(bad), Err(TextValueError::InvalidInteger)), "{bad:?}");
        }
        for big in [&b"2147483648"[..], b"-2147483649", b"99999999999999999999"] {
            assert!(matches!(parse_i32(big), Err(TextValueError::IntegerOutOfRange)), "{big:?}");
        }
    }

    // r[verify vcf_values.text_numbers]
    #[test]
    fn floats_take_the_vcf_special_forms_in_any_case() {
        for (text, expected) in [
            (&b"Inf"[..], f32::INFINITY),
            (b"-inf", f32::NEG_INFINITY),
            (b"+INFINITY", f32::INFINITY),
            (b"-Infinity", f32::NEG_INFINITY),
            (b"1e3", 1000.0),
            (b".5", 0.5),
            (b"5.", 5.0),
            (b"1e100", f32::INFINITY),
        ] {
            assert_eq!(parse_f32(text).unwrap(), Some(expected), "{text:?}");
        }
        for nan in [&b"NaN"[..], b"nan", b"-NAN", b"+nan"] {
            assert!(parse_f32(nan).unwrap().unwrap().is_nan(), "{nan:?}");
        }
        assert_eq!(parse_f32(b".").unwrap(), None);
        for bad in [&b""[..], b"-", b" 1", b"inff", b"0x1p3", b"1,0", b"..", b"\xff"] {
            assert!(matches!(parse_f32(bad), Err(TextValueError::InvalidFloat)), "{bad:?}");
        }
    }

    // r[verify vcf_values.text_strings]
    #[test]
    fn strings_are_returned_as_stored() {
        assert_eq!(text_str(b"a%3Ab%2C").unwrap(), "a%3Ab%2C");
        assert_eq!(text_str("Ünïcode".as_bytes()).unwrap(), "Ünïcode");
        assert!(matches!(text_str(b"\xC3\x28"), Err(TextValueError::Utf8(_))));
    }

    fn join(values: &[Option<String>]) -> String {
        let tokens: Vec<&str> = values.iter().map(|v| v.as_deref().unwrap_or(".")).collect();
        tokens.join(",")
    }

    // r[verify vcf_values.text_numbers]
    #[hegel::test]
    fn int_vectors_read_back_as_joined(tc: TestCase) {
        let values =
            tc.draw(gs::vecs(gs::optional(gs::integers::<i32>())).min_size(1).max_size(40));
        let text = join(&values.iter().map(|v| v.map(|v| v.to_string())).collect::<Vec<_>>());
        let v = TextVec::<i32>::new(text.as_bytes()).unwrap();
        assert_eq!(v.len(), values.len());
        for (i, &expected) in values.iter().enumerate() {
            assert_eq!(v.get(i), Some(expected), "element {i}");
        }
        assert_eq!(v.get(values.len()), None);
        assert_eq!(v.collect::<Vec<_>>(), values);
    }

    #[test]
    fn a_bad_element_anywhere_fails_the_vector() {
        assert!(matches!(TextVec::<i32>::new(b"1,2,x"), Err(TextValueError::InvalidInteger)));
        assert!(matches!(TextVec::<i32>::new(b"1,,2"), Err(TextValueError::InvalidInteger)));
        assert!(matches!(TextVec::<i32>::new(b""), Err(TextValueError::InvalidInteger)));
        assert!(matches!(TextVec::<i32>::new(b"1,"), Err(TextValueError::InvalidInteger)));
        assert!(matches!(TextVec::<f32>::new(b"0.5,abc"), Err(TextValueError::InvalidFloat)));
        let v = TextVec::<f32>::new(b".,0.5,Inf").unwrap();
        assert_eq!(v.collect::<Vec<_>>(), [None, Some(0.5), Some(f32::INFINITY)]);
    }

    // r[verify vcf_values.gt]
    #[test]
    fn malformed_text_genotypes_are_typed_errors() {
        for bad in
            ["", "/", "|", "0/", "0//1", "/|0", "a", "-1", "+1", "0/1 ", "0\\1", "0,1", "0:1"]
        {
            assert!(
                matches!(TextGt::new(bad.as_bytes(), true), Err(TextValueError::InvalidGenotype)),
                "{bad:?}"
            );
        }
        assert!(matches!(
            TextGt::new(b"0/65536", true),
            Err(TextValueError::GenotypeAlleleOutOfRange)
        ));
        let gt: Vec<_> = TextGt::new(b"65535|00", true).unwrap().map(|a| a.allele).collect();
        assert_eq!(gt, [Some(65535), Some(0)]);
    }
}
