//! [`RawBamRecord`] → SAM text, byte for byte what htslib's `sam_format1`
//! prints (`r[sam.format]`), and the header as `samtools view -h --no-PG`
//! prints it (`r[sam.format.header]`).
//!
//! Float tags print as htslib prints them: `f` and `d` values with C's `%g`,
//! `B:f` elements with htslib's `kputd`. Both are reimplemented here, so no
//! libc call and no allocation beyond growing `out`.

use crate::bam::{
    BamHeader, RawBamRecord,
    aux::{AuxLayoutError, AuxValue, tag_end_in, whole_tag},
};
use crate::io::text_scan::find_byte;
use crate::reader::Tid;
use fearless_simd::{Level, dispatch, prelude::*};
use fearless_simd_macros::simd;
use std::fmt::{self, Write as _};
use std::num::NonZeroI32;

/// Which SAM field named a reference the header lacks.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReferenceField {
    /// RNAME, the record's own reference.
    Rname,
    /// RNEXT, the mate's reference.
    Rnext,
}

impl fmt::Display for ReferenceField {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Rname => "RNAME",
            Self::Rnext => "RNEXT",
        })
    }
}

/// Why a record or header could not be formatted as SAM.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum SamFormatError {
    // r[impl sam.format.unknown_reference]
    #[error("{field} reference id {} is not a target of the header ({target_count} targets)", .tid.as_u32())]
    UnknownReference { field: ReferenceField, tid: Tid, target_count: usize },
    // r[impl sam.format.aux_errors]
    #[error("malformed aux data")]
    Aux(#[from] AuxLayoutError),
    #[error("SAM header line {line} does not start with '@'")]
    MalformedHeaderLine { line: usize },
}

/// SAM's CIGAR letters by BAM op code, `?` for the reserved codes 10–15
/// (htslib's `BAM_CIGAR_STR "??????"`).
const CIGAR_CHARS: &[u8; 16] = b"MIDNSHP=XB??????";

// r[impl sam.format]
/// Append `record` to `out` as one SAM line, without a trailing newline,
/// exactly as htslib's `sam_format1` prints it after reading the record from a
/// BAM file with `header`.
///
/// # Errors
/// [`SamFormatError::UnknownReference`] for a reference id `header` lacks,
/// [`SamFormatError::Aux`] for a malformed aux block. `out` is then as it was.
pub fn format_sam_record(
    record: &RawBamRecord,
    header: &BamHeader,
    out: &mut Vec<u8>,
) -> Result<(), SamFormatError> {
    let start = out.len();
    let result = write_record(record, header, out);
    if result.is_err() {
        out.truncate(start);
    }
    result
}

fn target_name(
    header: &BamHeader,
    tid: Tid,
    field: ReferenceField,
) -> Result<&str, SamFormatError> {
    header.target_name(tid).ok_or(SamFormatError::UnknownReference {
        field,
        tid,
        target_count: header.target_count(),
    })
}

fn write_record(
    record: &RawBamRecord,
    header: &BamHeader,
    out: &mut Vec<u8>,
) -> Result<(), SamFormatError> {
    let mut int = itoa::Buffer::new();
    let mate = record.mate_fields();
    let rname = match record.tid() {
        Some(tid) => target_name(header, tid, ReferenceField::Rname)?,
        None => "*",
    };
    // r[impl sam.format.fixed_fields]
    let rnext = match mate.tid {
        None => "*",
        Some(tid) if Some(tid) == record.tid() => "=",
        Some(tid) => target_name(header, tid, ReferenceField::Rnext)?,
    };
    let one_based = |pos: Option<crate::bam::Pos0>| pos.map_or(0, |p| p.as_i64().saturating_add(1));

    // r[impl sam.format.qname]
    out.extend_from_slice(record.qname_field());
    out.push(b'\t');
    out.extend_from_slice(int.format(record.flags().raw()).as_bytes());
    out.push(b'\t');
    out.extend_from_slice(rname.as_bytes());
    out.push(b'\t');
    out.extend_from_slice(int.format(one_based(record.pos())).as_bytes());
    out.push(b'\t');
    out.extend_from_slice(int.format(record.mapq().as_byte()).as_bytes());
    out.push(b'\t');

    // r[impl sam.format.cigar]
    // r[impl sam.format.long_cigar]
    let long_cigar = record.long_cigar();
    let cigar = long_cigar.unwrap_or(record.stored_cigar_bytes());
    if cigar.is_empty() {
        out.push(b'*');
    }
    for op in cigar.as_chunks::<4>().0 {
        let op = u32::from_le_bytes(*op);
        out.extend_from_slice(int.format(op >> 4).as_bytes());
        out.push(CIGAR_CHARS.get((op & 0xF) as usize).copied().unwrap_or(b'?'));
    }
    out.push(b'\t');

    out.extend_from_slice(rnext.as_bytes());
    out.push(b'\t');
    out.extend_from_slice(int.format(one_based(mate.pos)).as_bytes());
    out.push(b'\t');
    out.extend_from_slice(int.format(mate.template_len.map_or(0, NonZeroI32::get)).as_bytes());
    out.push(b'\t');

    // r[impl sam.format.seq_qual]
    if record.seq_len() == 0 {
        out.extend_from_slice(b"*\t*");
    } else {
        record.decode_seq_into(out);
        out.push(b'\t');
        match record.qual() {
            None => out.push(b'*'),
            Some(qual) => {
                out.extend(qual.iter().map(|q| q.as_byte().wrapping_add(33)));
            }
        }
    }

    write_aux_at(Level::new(), record.aux().as_bytes(), long_cigar.is_some(), out, &mut int)?;
    Ok(())
}

// r[impl io.simd_portable]
fn write_aux_at(
    level: Level,
    aux: &[u8],
    skip_cg: bool,
    out: &mut Vec<u8>,
    int: &mut itoa::Buffer,
) -> Result<(), AuxLayoutError> {
    dispatch!(level, simd => write_aux_simd(simd, aux, skip_cg, out, int))
}

// r[impl sam.format.aux]
// r[impl sam.format.aux_walk]
/// Every tag as `\tTAG:TYPE:VALUE`, checked as `Aux::iter_checked` checks it,
/// in one kernel: each `Z`/`H` string's NUL is found once, without a dispatch.
#[simd]
fn write_aux_simd<S: Simd>(
    simd: S,
    aux: &[u8],
    mut skip_cg: bool,
    out: &mut Vec<u8>,
    int: &mut itoa::Buffer,
) -> Result<(), AuxLayoutError> {
    let mut at = 0usize;
    while at < aux.len() {
        let end = tag_end_in(simd, aux, at)?;
        let (tag, value) =
            aux.get(at..end).and_then(whole_tag).ok_or(AuxLayoutError::Truncated { offset: at })?;
        at = end;
        // The first CG tag is the one `bam_tag2cigar` moved into the CIGAR.
        if skip_cg && tag == *b"CG" {
            skip_cg = false;
            continue;
        }
        out.push(b'\t');
        out.extend_from_slice(&tag);
        out.push(b':');
        write_aux_value(&value, out, int);
    }
    Ok(())
}

fn write_aux_value(value: &AuxValue<'_>, out: &mut Vec<u8>, int: &mut itoa::Buffer) {
    let mut put_int = |out: &mut Vec<u8>, n: i64| {
        out.extend_from_slice(b"i:");
        out.extend_from_slice(int.format(n).as_bytes());
    };
    match *value {
        AuxValue::Char(c) => out.extend_from_slice(&[b'A', b':', c]),
        AuxValue::I8(n) => put_int(out, n.into()),
        AuxValue::U8(n) => put_int(out, n.into()),
        AuxValue::I16(n) => put_int(out, n.into()),
        AuxValue::U16(n) => put_int(out, n.into()),
        AuxValue::I32(n) => put_int(out, n.into()),
        AuxValue::U32(n) => put_int(out, n.into()),
        AuxValue::Float(x) => {
            out.extend_from_slice(b"f:");
            format_g(x.into(), out);
        }
        AuxValue::Double(x) => {
            out.extend_from_slice(b"d:");
            format_g(x, out);
        }
        AuxValue::String(s) => {
            out.extend_from_slice(b"Z:");
            out.extend_from_slice(s);
        }
        AuxValue::Hex(s) => {
            out.extend_from_slice(b"H:");
            out.extend_from_slice(s);
        }
        AuxValue::ArrayI8(a) => write_array(out, b'c', a, |out, [b]: [u8; 1]| {
            out.extend_from_slice(int.format(b.cast_signed()).as_bytes());
        }),
        AuxValue::ArrayU8(a) => write_array(out, b'C', a, |out, [b]: [u8; 1]| {
            out.extend_from_slice(int.format(b).as_bytes());
        }),
        AuxValue::ArrayI16(a) => write_array(out, b's', a, |out, b| {
            out.extend_from_slice(int.format(i16::from_le_bytes(b)).as_bytes());
        }),
        AuxValue::ArrayU16(a) => write_array(out, b'S', a, |out, b| {
            out.extend_from_slice(int.format(u16::from_le_bytes(b)).as_bytes());
        }),
        AuxValue::ArrayI32(a) => write_array(out, b'i', a, |out, b| {
            out.extend_from_slice(int.format(i32::from_le_bytes(b)).as_bytes());
        }),
        AuxValue::ArrayU32(a) => write_array(out, b'I', a, |out, b| {
            out.extend_from_slice(int.format(u32::from_le_bytes(b)).as_bytes());
        }),
        AuxValue::ArrayFloat(a) => write_array(out, b'f', a, |out, b| {
            format_kputd(f32::from_le_bytes(b).into(), out);
        }),
    }
}

/// `B:<subtype>` and `,<element>` per `N`-byte element of `bytes`.
fn write_array<const N: usize>(
    out: &mut Vec<u8>,
    subtype: u8,
    bytes: &[u8],
    mut element: impl FnMut(&mut Vec<u8>, [u8; N]),
) {
    out.extend_from_slice(&[b'B', b':', subtype]);
    for e in bytes.as_chunks::<N>().0 {
        out.push(b',');
        element(out, *e);
    }
}

/// A `fmt::Write` over a stack buffer, failing when it would overflow.
struct StackBuf<const N: usize> {
    buf: [u8; N],
    len: usize,
}

impl<const N: usize> fmt::Write for StackBuf<N> {
    fn write_str(&mut self, s: &str) -> fmt::Result {
        let end = self.len.checked_add(s.len()).ok_or(fmt::Error)?;
        self.buf.get_mut(self.len..end).ok_or(fmt::Error)?.copy_from_slice(s.as_bytes());
        self.len = end;
        Ok(())
    }
}

/// `value`, a finite non-negative number, rounded to six significant digits:
/// the digits and the decimal exponent of the first.
fn six_digits(value: f64) -> Option<([u8; 6], i32)> {
    // Longest output: "1.00000e-308" — 12 bytes.
    let mut text = StackBuf { buf: [0; 24], len: 0 };
    // Rust's float formatting is exact and rounds ties to even, as glibc's
    // printf does on the exact binary value.
    write!(text, "{value:.5e}").ok()?;
    let text = text.buf.get(..text.len)?;
    let (mantissa, exponent) = text.split_at(text.iter().position(|&b| b == b'e')?);
    let &[d0, b'.', d1, d2, d3, d4, d5] = mantissa else { return None };
    let exponent = std::str::from_utf8(exponent.get(1..)?).ok()?.parse().ok()?;
    Some(([d0, d1, d2, d3, d4, d5], exponent))
}

fn trim_trailing_zeros(digits: &[u8]) -> &[u8] {
    let len = digits.iter().rposition(|&b| b != b'0').map_or(0, |i| i.saturating_add(1));
    digits.get(..len).unwrap_or_default()
}

/// `int` and, when the fraction has a nonzero digit, `.` and the fraction
/// without trailing zeros.
fn put_decimal(out: &mut Vec<u8>, int: &[u8], frac: &[u8]) {
    out.extend_from_slice(int);
    let frac = trim_trailing_zeros(frac);
    if !frac.is_empty() {
        out.push(b'.');
        out.extend_from_slice(frac);
    }
}

// r[impl sam.format.float_g]
/// Append `value` as C's `printf("%g", value)` prints it: six significant
/// digits, fixed notation for decimal exponents −4 to 5 and `d.ddddde±XX`
/// otherwise, trailing zeros dropped. NaN prints `nan`, or `-nan` with its sign
/// bit set (glibc and musl; macOS's libc prints `nan` for both).
pub(crate) fn format_g(value: f64, out: &mut Vec<u8>) {
    if value.is_sign_negative() {
        out.push(b'-');
    }
    if value.is_nan() {
        out.extend_from_slice(b"nan");
        return;
    }
    let value = value.abs();
    if value.is_infinite() {
        out.extend_from_slice(b"inf");
        return;
    }
    let Some((digits, exp)) = six_digits(value) else {
        debug_assert!(false, "`{{:.5e}}` of a finite f64 is six digits and an exponent");
        return;
    };
    write_g_digits(&digits, exp, out);
}

/// `%g`'s layout of six rounded significant digits with decimal exponent
/// `exp`.
fn write_g_digits(digits: &[u8; 6], exp: i32, out: &mut Vec<u8>) {
    match exp {
        0..=5 => {
            let (int, frac) = digits.split_at((exp.unsigned_abs() as usize).saturating_add(1));
            put_decimal(out, int, frac);
        }
        -4..=-1 => {
            out.extend_from_slice(b"0.");
            out.extend(std::iter::repeat_n(b'0', (exp.unsigned_abs() as usize).saturating_sub(1)));
            out.extend_from_slice(trim_trailing_zeros(digits));
        }
        _ => {
            let (int, frac) = digits.split_at(1);
            put_decimal(out, int, frac);
            out.extend_from_slice(if exp < 0 { b"e-" } else { b"e+" });
            if exp.unsigned_abs() < 10 {
                out.push(b'0');
            }
            out.extend_from_slice(itoa::Buffer::new().format(exp.unsigned_abs()).as_bytes());
        }
    }
}

// r[impl sam.format.float_kputd]
/// Append `value` as htslib's `kputd` prints it — `%g` (see [`format_g`]),
/// except that magnitudes in `[0.0001, 999999]` are scaled to a six-digit
/// integer with round-half-even on the scaled `f64`, which can round
/// differently from `%g`'s exact rounding.
pub(crate) fn format_kputd(value: f64, out: &mut Vec<u8>) {
    if value == 0.0 {
        out.extend_from_slice(if value.is_sign_negative() { b"-0" } else { b"0" });
        return;
    }
    let mut d = value;
    if d < 0.0 {
        out.push(b'-');
        d = -d;
    }
    if !(0.0001..=999_999.0).contains(&d) {
        format_g(d, out);
        return;
    }
    // (scale, digits after the point): `d * scale` has six integer digits.
    let (scale, frac_digits) = match d {
        _ if d < 0.001 => (1e9, 9),
        _ if d < 0.01 => (1e8, 8),
        _ if d < 0.1 => (1e7, 7),
        _ if d < 1.0 => (1e6, 6),
        _ if d < 10.0 => (1e5, 5),
        _ if d < 100.0 => (1e4, 4),
        _ if d < 1000.0 => (1e3, 3),
        _ if d < 10000.0 => (1e2, 2),
        _ if d < 100_000.0 => (1e1, 1),
        _ => (1.0, 0),
    };
    let scaled = (d * scale).round_ties_even();
    debug_assert!((99_999.0..=1_000_000.0).contains(&scaled), "{d} scaled to {scaled}");
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "`scaled` is a whole number in 99_999..=1_000_000 by the band choice"
    )]
    let scaled = scaled as u32;
    let mut buf = itoa::Buffer::new();
    let digits = buf.format(scaled).as_bytes();
    match digits.len().checked_sub(frac_digits) {
        Some(int_len) if int_len > 0 => {
            let (int, frac) = digits.split_at(int_len);
            put_decimal(out, int, frac);
        }
        _ => {
            out.extend_from_slice(b"0.");
            out.extend(std::iter::repeat_n(b'0', frac_digits.saturating_sub(digits.len())));
            out.extend_from_slice(trim_trailing_zeros(digits));
        }
    }
}

// r[impl sam.format.header]
/// Append `header` as SAM header lines, as `samtools view -h --no-PG` prints a
/// BAM file's: the text up to its first NUL, a final newline if it lacks one,
/// and `@SQ` lines from the targets when the text has none.
///
/// # Errors
/// [`SamFormatError::MalformedHeaderLine`] for a line that does not start with
/// `@`, as htslib rejects such a header. `out` is then as it was.
pub fn format_sam_header(header: &BamHeader, out: &mut Vec<u8>) -> Result<(), SamFormatError> {
    let text = header.header_text().as_bytes();
    let text = text.get(..find_byte(text, 0).unwrap_or(text.len())).unwrap_or_default();
    let mut has_sq = false;
    let mut rest = text;
    let mut line_number = 0usize;
    while !rest.is_empty() {
        line_number = line_number.saturating_add(1);
        let (line, tail) = match find_byte(rest, b'\n') {
            Some(nl) => (rest.get(..nl), rest.get(nl.saturating_add(1)..)),
            None => (Some(rest), None),
        };
        let line = line.unwrap_or_default();
        if line.first() != Some(&b'@') {
            return Err(SamFormatError::MalformedHeaderLine { line: line_number });
        }
        has_sq |= line.starts_with(b"@SQ\t");
        rest = tail.unwrap_or_default();
    }
    out.extend_from_slice(text);
    if text.last().is_some_and(|&b| b != b'\n') {
        out.push(b'\n');
    }
    if !has_sq {
        let mut int = itoa::Buffer::new();
        for target in header.targets() {
            out.extend_from_slice(b"@SQ\tSN:");
            out.extend_from_slice(target.name.as_bytes());
            out.extend_from_slice(b"\tLN:");
            out.extend_from_slice(int.format(target.length).as_bytes());
            out.push(b'\n');
        }
    }
    Ok(())
}

#[cfg(test)]
mod kernel_tests {
    use super::*;
    use crate::bam::aux::kernel_tests::aux_block;
    use hegel::prelude::*;

    // r[verify sam.format.aux_walk]
    // r[verify io.simd_portable]
    #[hegel::test]
    fn write_aux_every_level_matches_fallback(tc: TestCase) {
        let aux = tc.draw(aux_block());
        let skip_cg = tc.draw(gs::booleans());
        let run = |level| {
            let mut out = b"prefix".to_vec();
            let result = write_aux_at(level, &aux, skip_cg, &mut out, &mut itoa::Buffer::new());
            (result, out)
        };
        let expected = run(Level::fallback());
        for level in crate::simd_levels::levels() {
            assert_eq!(run(level), expected, "{level:?}");
        }
    }
}

#[cfg(test)]
#[path = "../../tests/support/float_oracles.rs"]
mod float_oracles;

/// The float printers against their C originals: `format_g` against libc's
/// `snprintf("%g")`, `format_kputd` against htslib's `kputd`.
#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::indexing_slicing, reason = "test code")]
mod float_tests {
    use super::float_oracles::{arb_f32, arb_f64, expect_from_libc, htslib_kputd, libc_g};
    use super::{format_g, format_kputd};
    use hegel::TestCase;
    use hegel::generators as gs;

    fn ours(f: fn(f64, &mut Vec<u8>), value: f64) -> Vec<u8> {
        // Something already in `out` must stay there.
        let mut out = b"x".to_vec();
        f(value, &mut out);
        assert_eq!(out[0], b'x');
        out.split_off(1)
    }

    // r[verify sam.format.float_g]
    #[hegel::test]
    fn format_g_matches_libc_for_doubles(tc: TestCase) {
        let value = arb_f64(&tc, true);
        let libc = expect_from_libc(value, libc_g(value));
        assert_eq!(
            ours(format_g, value).escape_ascii().to_string(),
            libc.escape_ascii().to_string(),
            "{value:e}"
        );
    }

    // r[verify sam.format.float_g]
    #[hegel::test]
    fn format_g_matches_libc_for_floats(tc: TestCase) {
        let value = f64::from(arb_f32(&tc, true));
        let libc = expect_from_libc(value, libc_g(value));
        assert_eq!(
            ours(format_g, value).escape_ascii().to_string(),
            libc.escape_ascii().to_string(),
            "{value:e}"
        );
    }

    // r[verify sam.format.float_kputd]
    #[hegel::test]
    fn format_kputd_matches_htslib(tc: TestCase) {
        // `B:f` elements are f32, but kputd takes any double.
        let value = if tc.draw(gs::booleans()) {
            f64::from(arb_f32(&tc, true))
        } else {
            arb_f64(&tc, true)
        };
        let hts = expect_from_libc(value, htslib_kputd(value));
        assert_eq!(
            ours(format_kputd, value).escape_ascii().to_string(),
            hts.escape_ascii().to_string(),
            "{value:e}"
        );
    }

    // r[verify sam.format.float_g]
    // r[verify sam.format.float_kputd]
    #[test]
    fn float_printers_on_known_values() {
        let g = |v| String::from_utf8(ours(format_g, v)).unwrap();
        let k = |v| String::from_utf8(ours(format_kputd, v)).unwrap();
        for (value, expected) in [
            (0.0, "0"),
            (-0.0, "-0"),
            (1.0, "1"),
            (0.5, "0.5"),
            (1e-5, "1e-05"),
            (0.0001, "0.0001"),
            (123_456.0, "123456"),
            (1_234_567.0, "1.23457e+06"),
            (1_234_565.0, "1.23456e+06"),
            (9.999_995, "10"),
            (1e100, "1e+100"),
            (f64::INFINITY, "inf"),
            (f64::NEG_INFINITY, "-inf"),
            (f64::NAN, "nan"),
            (-f64::NAN, "-nan"),
            (f64::from_bits(1), "4.94066e-324"),
        ] {
            assert_eq!(g(value), expected, "%g of {value:e}");
        }
        // kputd agrees with %g except where its scaled rounding differs.
        // On doubles near a tie the scaled product rounds the other way (never on
        // an f32, which is all `B:f` holds: every f32 in kputd's range agrees).
        assert_eq!((g(1.000_005), k(1.000_005)), ("1.00001".into(), "1".into()));
        assert_eq!((g(1000.015), k(1000.015)), ("1000.01".into(), "1000.02".into()));
        assert_eq!((g(1.000_045), k(1.000_045)), ("1.00005".into(), "1.00004".into()));
        assert_eq!(k(999_999.0), "999999");
        assert_eq!(k(999_999.5), "1e+06");
        assert_eq!(k(-0.0), "-0");
        assert_eq!(k(-f64::NAN), "-nan");
        assert_eq!(k(f64::NEG_INFINITY), "-inf");
    }
}
