//! VCF text formatting helpers used by the unified writer.

/// Format a float like C's `%g` with 6 significant digits — no trailing zeros,
/// no trailing decimal point. This matches htslib/bcftools VCF text output.
///
/// Examples: 35.89775 → "35.8978", 60.0 → "60", 0.777778 → "0.777778"
///
/// Magnitudes too small for the fixed form to fit six significant digits in
/// 32 bytes (below about `1e-25`) fall back to scientific notation, the way
/// `%g` does — `5.16988e-26`. Those values are not quality scores, but they
/// do turn up in INFO fields carrying p-values, and refusing to write them
/// would lose the record.
#[expect(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::indexing_slicing,
    reason = "magnitude is log10 of f64 (range -308..=308), precision is a literal 6, cursor position is bounded by the 32-byte tmp buffer; all casts and slices are safe"
)]
pub(crate) fn write_float_g(buf: &mut Vec<u8>, v: f32) -> Result<(), WriteError> {
    let v = f64::from(v);
    let precision = 6usize;

    if v == 0.0 {
        buf.push(b'0');
        return Ok(());
    }

    let magnitude = v.abs().log10().floor() as i32;
    let decimal_places = (precision as i32).saturating_sub(magnitude).saturating_sub(1);
    let decimal_places = decimal_places.max(0) as usize;

    // Format into a small stack buffer
    let mut tmp = [0u8; 32];
    let mut cursor = std::io::Cursor::new(&mut tmp[..]);
    if std::io::Write::write_fmt(&mut cursor, format_args!("{v:.decimal_places$}")).is_err() {
        // The fixed form did not fit. Scientific notation always does: an
        // `f64` mantissa of six significant digits plus `e-308` is 13 bytes.
        return write_float_exponential(buf, v);
    }
    let end = cursor.position() as usize;

    // Strip trailing zeros after decimal point, then trailing dot
    debug_assert!(end <= tmp.len(), "cursor position must not exceed buffer size");
    let formatted = &tmp[..end];
    let trimmed = if formatted.contains(&b'.') {
        let without_zeros = formatted
            .strip_suffix(b"0")
            .map(|s| {
                let mut s = s;
                while let Some(stripped) = s.strip_suffix(b"0") {
                    s = stripped;
                }
                s
            })
            .unwrap_or(formatted);
        without_zeros.strip_suffix(b".").unwrap_or(without_zeros)
    } else {
        formatted
    };

    buf.extend_from_slice(trimmed);
    Ok(())
}

/// The scientific fallback of [`write_float_g`]: six significant digits, with
/// the mantissa's trailing zeros stripped the same way.
fn write_float_exponential(buf: &mut Vec<u8>, v: f64) -> Result<(), WriteError> {
    let mut tmp = [0u8; 32];
    let mut cursor = std::io::Cursor::new(&mut tmp[..]);
    std::io::Write::write_fmt(&mut cursor, format_args!("{v:.5e}"))
        .map_err(|_source| WriteError::FormattedFloatLongerThan32Chars)?;
    #[expect(
        clippy::cast_possible_truncation,
        reason = "the cursor cannot pass the 32-byte buffer it was built from"
    )]
    let end = cursor.position() as usize;
    let formatted = tmp.get(..end).ok_or(WriteError::FormattedFloatLongerThan32Chars)?;

    let Some(e) = formatted.iter().position(|&b| b == b'e') else {
        buf.extend_from_slice(formatted);
        return Ok(());
    };
    let (mantissa, exponent) =
        formatted.split_at_checked(e).ok_or(WriteError::FailedToStripTrailingZeros)?;
    let mut mantissa = mantissa;
    if mantissa.contains(&b'.') {
        while let Some(stripped) = mantissa.strip_suffix(b"0") {
            mantissa = stripped;
        }
        mantissa = mantissa.strip_suffix(b".").unwrap_or(mantissa);
    }
    buf.extend_from_slice(mantissa);
    buf.extend_from_slice(exponent);
    Ok(())
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
    #[test]
    fn float_format_known_vectors() {
        assert_eq!(formatted(0.0), "0");
        assert_eq!(formatted(60.0), "60");
        assert_eq!(formatted(35.897_75), "35.8978");
        assert_eq!(formatted(0.777_778), "0.777778");
        // Small enough that the fixed form no longer fits six significant
        // digits in the buffer, so the writer switches to `%g`'s other half.
        assert_eq!(formatted(5.169_879e-26), "5.16988e-26");
        assert_eq!(formatted(f32::MIN_POSITIVE), "1.17549e-38");
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
