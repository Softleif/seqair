//! VCF text formatting helpers used by the unified writer.

/// Format a float as C's `%g` with 6 significant digits — no trailing zeros,
/// no trailing decimal point — which is what htslib and bcftools write into
/// VCF text, so seqair's output can be diffed against theirs.
///
/// `%g` picks the shorter of the two forms by the value's decimal exponent
/// `e`: scientific when `e < -4` or `e >= 6`, fixed otherwise. The exponent
/// itself is written C's way, signed and at least two digits.
///
/// Examples: 35.89775 → `35.8978`, 60.0 → `60`, 0.777778 → `0.777778`,
/// 0.0001 → `0.0001`, 0.00001 → `1e-05`, 999999 → `999999`,
/// 1234567 → `1.23457e+06`.
pub(crate) fn write_float_g(buf: &mut Vec<u8>, v: f32) -> Result<(), WriteError> {
    let v = f64::from(v);

    if v == 0.0 {
        // `-0.0 == 0.0`, but C prints `-0` and so must this: the sign is the
        // difference between a value that round-trips and one that does not.
        if v.is_sign_negative() {
            buf.push(b'-');
        }
        buf.push(b'0');
        return Ok(());
    }
    if !v.is_finite() {
        // VCF has no spelling for these; write what `Display` would and let
        // the caller's validation deal with it rather than silently dropping
        // the value.
        buf.extend_from_slice(v.to_string().as_bytes());
        return Ok(());
    }

    // `%g` chooses on the exponent of the value *as rounded to `PRECISION`
    // significant digits*, not of the value itself: `0.0001f32` is a hair
    // below `1e-4`, but it rounds to `1.00000e-04`, so C prints `0.0001` and
    // not `1e-04`. Rendering the scientific form first and reading the
    // exponent back off it is that rounding, done once.
    let mut tmp = [0u8; 32];
    let mut cursor = std::io::Cursor::new(&mut tmp[..]);
    std::io::Write::write_fmt(&mut cursor, format_args!("{v:.*e}", SIG_DECIMALS))
        .map_err(|_source| WriteError::FormattedFloatLongerThan32Chars)?;
    #[expect(
        clippy::cast_possible_truncation,
        reason = "the cursor cannot pass the 32-byte buffer it was built from"
    )]
    let end = cursor.position() as usize;
    let scientific = tmp.get(..end).ok_or(WriteError::FormattedFloatLongerThan32Chars)?;

    let split =
        scientific.iter().position(|&b| b == b'e').ok_or(WriteError::FailedToStripTrailingZeros)?;
    let (mantissa, exponent_text) =
        scientific.split_at_checked(split).ok_or(WriteError::FailedToStripTrailingZeros)?;
    let exponent_text = exponent_text.get(1..).ok_or(WriteError::FailedToStripTrailingZeros)?;
    let exponent: i32 = std::str::from_utf8(exponent_text)
        .map_err(|_source| WriteError::FailedToStripTrailingZeros)?
        .parse()
        .map_err(|_source| WriteError::FailedToStripTrailingZeros)?;

    if !(-4..PRECISION).contains(&exponent) {
        write_exponential(buf, mantissa, exponent);
        Ok(())
    } else {
        write_fixed(buf, v, exponent)
    }
}

/// Significant digits, as C's `%g` counts them, and the decimals that leaves
/// after the leading one in the scientific form.
const PRECISION: i32 = 6;
const SIG_DECIMALS: usize = 5;

/// `%f` with the decimals `%g` would use at this exponent, trailing zeros and
/// a bare decimal point stripped.
fn write_fixed(buf: &mut Vec<u8>, v: f64, exponent: i32) -> Result<(), WriteError> {
    #[expect(
        clippy::cast_sign_loss,
        reason = "the caller only takes this path for exponent in -4..6, so the difference is 1..=10"
    )]
    let decimals = PRECISION.saturating_sub(exponent).saturating_sub(1).max(0) as usize;

    // The widest form this path can produce is one integer digit plus a point
    // plus ten decimals.
    let mut tmp = [0u8; 32];
    let mut cursor = std::io::Cursor::new(&mut tmp[..]);
    std::io::Write::write_fmt(&mut cursor, format_args!("{v:.decimals$}"))
        .map_err(|_source| WriteError::FormattedFloatLongerThan32Chars)?;
    #[expect(
        clippy::cast_possible_truncation,
        reason = "the cursor cannot pass the 32-byte buffer it was built from"
    )]
    let end = cursor.position() as usize;
    let formatted = tmp.get(..end).ok_or(WriteError::FormattedFloatLongerThan32Chars)?;
    buf.extend_from_slice(strip_trailing_zeros(formatted));
    Ok(())
}

/// The scientific form, with the mantissa's trailing zeros stripped and the
/// exponent written C's way: a sign, and at least two digits.
fn write_exponential(buf: &mut Vec<u8>, mantissa: &[u8], exponent: i32) {
    buf.extend_from_slice(strip_trailing_zeros(mantissa));
    buf.push(b'e');
    buf.push(if exponent < 0 { b'-' } else { b'+' });
    let magnitude = exponent.unsigned_abs();
    if magnitude < 10 {
        buf.push(b'0');
    }
    let mut itoa = itoa::Buffer::new();
    buf.extend_from_slice(itoa.format(magnitude).as_bytes());
}

/// Drop trailing zeros after a decimal point, then the point itself. A number
/// with no point is returned untouched — `100` must not become `1`.
fn strip_trailing_zeros(formatted: &[u8]) -> &[u8] {
    if !formatted.contains(&b'.') {
        return formatted;
    }
    let mut trimmed = formatted;
    while let Some(stripped) = trimmed.strip_suffix(b"0") {
        trimmed = stripped;
    }
    trimmed.strip_suffix(b".").unwrap_or(trimmed)
}

/// Percent-encode special characters per VCF spec §1.0.2.
pub(crate) fn percent_encode_into(buf: &mut Vec<u8>, data: &[u8]) {
    for &b in data {
        match b {
            b':' => buf.extend_from_slice(b"%3A"),
            b';' => buf.extend_from_slice(b"%3B"),
            b'=' => buf.extend_from_slice(b"%3D"),
            b'%' => buf.extend_from_slice(b"%25"),
            b',' => buf.extend_from_slice(b"%2C"),
            b'\t' => buf.extend_from_slice(b"%09"),
            b'\n' => buf.extend_from_slice(b"%0A"),
            b'\r' => buf.extend_from_slice(b"%0D"),
            _ => buf.push(b),
        }
    }
}

#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum WriteError {
    #[error("failed to write formatted string")]
    WriteFormattedString { source: std::io::Error },
    #[error("Failed to strip trailing zeros from formatted float")]
    FailedToStripTrailingZeros,
    #[error("formatted float exceeds 32 bytes")]
    FormattedFloatLongerThan32Chars,
}

#[cfg(test)]
mod tests {
    use super::*;
    use hegel::prelude::*;

    fn formatted(v: f32) -> String {
        let mut buf = Vec::new();
        write_float_g(&mut buf, v).expect("every finite float must be writable");
        String::from_utf8(buf).expect("ASCII")
    }

    // r[verify vcf_writer.float_precision]
    /// Pinned against what `bcftools query` prints for the same values, which
    /// is C's `%g` with six significant digits.
    #[test]
    fn float_format_known_vectors() {
        assert_eq!(formatted(0.0), "0");
        assert_eq!(formatted(-0.0), "-0");
        assert_eq!(formatted(60.0), "60");
        assert_eq!(formatted(0.5), "0.5");
        assert_eq!(formatted(35.897_75), "35.8978");
        assert_eq!(formatted(0.777_778), "0.777778");
        // The boundaries `%g` switches on: `e < -4` or `e >= 6`.
        assert_eq!(formatted(0.0001), "0.0001");
        assert_eq!(formatted(0.000_123_456), "0.000123456");
        assert_eq!(formatted(0.00001), "1e-05");
        assert_eq!(formatted(6.103_515_6e-5), "6.10352e-05");
        assert_eq!(formatted(999_999.0), "999999");
        assert_eq!(formatted(1_234_567.0), "1.23457e+06");
        assert_eq!(formatted(1.234_567_9e8), "1.23457e+08");
        // Three-digit exponents, and the smallest and largest finite f32.
        assert_eq!(formatted(5.169_879e-26), "5.16988e-26");
        assert_eq!(formatted(f32::MIN_POSITIVE), "1.17549e-38");
        assert_eq!(formatted(f32::MAX), "3.40282e+38");
        assert_eq!(formatted(-0.00001), "-1e-05");
    }

    // r[verify vcf_writer.float_precision]
    /// Whatever notation it picks, the text must parse back to the same `f32`
    /// to within six significant digits — the round-trip the rule demands.
    #[hegel::test]
    fn float_format_round_trips_to_six_significant_digits(tc: TestCase) {
        let v = tc.draw(gs::floats::<f32>().min_value(f32::MIN).max_value(f32::MAX));
        let text = formatted(v);
        let back: f32 = text.parse().expect("the writer emits parseable text");
        if v == 0.0 {
            assert_eq!(back, 0.0);
            return;
        }
        let relative = f64::from((back - v).abs()) / f64::from(v.abs());
        assert!(
            relative <= 1e-5,
            "{v} -> {text:?} -> {back} is off by {relative}, more than six digits allows"
        );
    }

    #[test]
    fn percent_encoding_special_chars() {
        let mut buf = Vec::new();
        percent_encode_into(&mut buf, b"hello:world;key=val%done,x\ty\nz\rend");
        let result = String::from_utf8(buf).unwrap();
        assert_eq!(result, "hello%3Aworld%3Bkey%3Dval%25done%2Cx%09y%0Az%0Dend");
    }
}
