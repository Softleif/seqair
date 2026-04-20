//! VCF text formatting helpers used by the unified writer.

/// Format a float like C's `%g` with 6 significant digits — no trailing zeros,
/// no trailing decimal point. This matches htslib/bcftools VCF text output.
///
/// Examples: 35.89775 → "35.8978", 60.0 → "60", 0.777778 → "0.777778"
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
    std::io::Write::write_fmt(&mut cursor, format_args!("{v:.decimal_places$}"))
        .map_err(|_source| WriteError::FormattedFloatLongerThan32Chars)?;
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

    #[test]
    fn percent_encoding_special_chars() {
        let mut buf = Vec::new();
        percent_encode_into(&mut buf, b"hello:world;key=val%done,x\ty\nz\rend");
        let result = String::from_utf8(buf).unwrap();
        assert_eq!(result, "hello%3Aworld%3Bkey%3Dval%25done%2Cx%09y%0Az%0Dend");
    }
}
