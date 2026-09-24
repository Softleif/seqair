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
    if write_float_g_exact(buf, v) {
        return Ok(());
    }
    write_float_g_fmt(buf, v)
}

/// `%g` through `core::fmt`: the fallback for what the integer path does not
/// cover (zero, non-finite, subnormal), and the oracle it is tested against.
/// `core::fmt`'s exact-precision printing falls back to bignum arithmetic
/// (Dragon4) for most inputs, and it runs twice here, which made this the
/// largest cost of writing VCF text.
fn write_float_g_fmt(buf: &mut Vec<u8>, v: f32) -> Result<(), WriteError> {
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
    // r[impl vcf_writer.float_non_finite]
    if !v.is_finite() {
        // VCF does have a spelling for these. §1.3 admits
        // `^[-+]?(INF|INFINITY|NAN)$`, case-insensitively, as a Float, and
        // §6.3.3 gives quiet NaN first-class status *alongside* the missing
        // sentinel rather than as a synonym for it. htslib hands non-finite
        // values to C's `%g`, which lowercases them; `Display` writes `NaN`,
        // and a case difference is exactly what this function exists to
        // avoid.
        buf.extend_from_slice(if v.is_nan() {
            b"nan".as_slice()
        } else if v.is_sign_negative() {
            b"-inf".as_slice()
        } else {
            b"inf".as_slice()
        });
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

/// `5^k` for `k` in `0..=MAX_POW5`.
const POW5: [u128; MAX_POW5 as usize + 1] = {
    let mut t = [1u128; MAX_POW5 as usize + 1];
    let mut i = 1;
    while i < t.len() {
        t[i] = t[i - 1] * 5;
        i += 1;
    }
    t
};
/// A 24-bit mantissa times `5^44` still fits a `u128` (2^24 · 5^44 < 2^127),
/// which covers every normal `f32`: the smallest, 1.18e-38, needs `k = 43`.
const MAX_POW5: u32 = 44;

/// `%g` of a normal, nonzero `f32` from its exact binary value, with integer
/// arithmetic only. Returns `false`, having written nothing, for zero,
/// non-finite and subnormal values, which take the `core::fmt` path.
///
/// `|v| = m · 2^q` exactly. With `E = floor(log10 |v|)`, the six significant
/// digits are `D = round_half_even(m · 2^q · 10^(5 - E))`, `10^5 <= D < 10^6`
/// — computed exactly, so it is the rounding C's `printf` does, not an
/// approximation of it.
// r[impl vcf_writer.float_precision]
#[allow(
    clippy::arithmetic_side_effects,
    reason = "exponents are within ±150 and digit counts within 0..=6; nothing here nears overflow"
)]
fn write_float_g_exact(buf: &mut Vec<u8>, v: f32) -> bool {
    let bits = v.to_bits();
    let biased = (bits >> 23) & 0xFF;
    if biased == 0 || biased == 0xFF {
        return false;
    }
    let m = u128::from((bits & 0x7F_FFFF) | 0x80_0000);
    #[expect(clippy::cast_possible_wrap, reason = "biased < 255")]
    let q = biased as i32 - 150;
    // floor(log2 |v|) = q + 23; times log10(2) ≈ 1233 / 4096 gives E or E - 1.
    let mut e = ((q + 23) * 1233) >> 12;
    let Some((mut d, round_up)) = scaled(m, q, 5 - e).and_then(|first| {
        if first.0 >= 1_000_000 {
            e = e.saturating_add(1);
            scaled(m, q, 5 - e)
        } else {
            Some(first)
        }
    }) else {
        return false;
    };
    debug_assert!((100_000..1_000_000).contains(&d), "d={d} for {v}");
    if round_up {
        d = d.saturating_add(1);
        if d == 1_000_000 {
            d = 100_000;
            e = e.saturating_add(1);
        }
    }

    let mut digits = [0u8; 6];
    let mut rest = d;
    for slot in digits.iter_mut().rev() {
        let digit = u8::try_from(rest % 10).unwrap_or_default();
        *slot = b'0'.saturating_add(digit);
        rest /= 10;
    }
    // Significant digits after trailing zeros are stripped: at least one.
    let sig = digits.iter().rposition(|&c| c != b'0').map_or(1, |i| i.saturating_add(1));
    let sig_digits = digits.get(..sig).unwrap_or_default();

    if bits >> 31 == 1 {
        buf.push(b'-');
    }
    if (-4..PRECISION).contains(&e) {
        #[expect(clippy::cast_sign_loss, reason = "e is in -4..6 here")]
        if e >= 0 {
            let int_len = e as usize + 1;
            buf.extend_from_slice(digits.get(..int_len).unwrap_or_default());
            if sig > int_len {
                buf.push(b'.');
                buf.extend_from_slice(sig_digits.get(int_len..).unwrap_or_default());
            }
        } else {
            buf.extend_from_slice(b"0.");
            let zeros = (-e - 1) as usize;
            buf.extend_from_slice(b"000".get(..zeros).unwrap_or_default());
            buf.extend_from_slice(sig_digits);
        }
    } else {
        let (lead, frac) = sig_digits.split_at(1);
        buf.extend_from_slice(lead);
        if !frac.is_empty() {
            buf.push(b'.');
            buf.extend_from_slice(frac);
        }
        buf.push(b'e');
        buf.push(if e < 0 { b'-' } else { b'+' });
        let magnitude = e.unsigned_abs();
        if magnitude < 10 {
            buf.push(b'0');
        }
        let mut itoa = itoa::Buffer::new();
        buf.extend_from_slice(itoa.format(magnitude).as_bytes());
    }
    true
}

/// `floor(m · 2^q · 10^k)` and whether rounding it half-to-even goes up, or
/// `None` outside the range the fixed-width arithmetic covers exactly.
fn scaled(m: u128, q: i32, k: i32) -> Option<(u64, bool)> {
    let (quotient, round_up) = if k >= 0 {
        // m · 5^k · 2^(q + k): multiply by the odd part, shift by the rest.
        let num = m.checked_mul(*POW5.get(k.unsigned_abs() as usize)?)?;
        let shift = q.checked_add(k)?;
        if shift >= 0 {
            (
                num.checked_shl(shift.unsigned_abs())
                    .filter(|n| n >> shift.unsigned_abs() == num)?,
                false,
            )
        } else {
            let s = shift.unsigned_abs();
            if s >= 128 {
                return Some((0, false));
            }
            let quotient = num >> s;
            let rem = num & ((1u128 << s).wrapping_sub(1));
            let half = 1u128 << s.saturating_sub(1);
            (quotient, rem > half || (rem == half && quotient & 1 == 1))
        }
    } else {
        // m · 2^q / 10^j with j = -k: exact integer division.
        let j = k.unsigned_abs();
        let pow10 = POW5.get(j as usize)?.checked_shl(j)?;
        let (num, den) = if q >= 0 {
            (m.checked_shl(q.unsigned_abs())?, pow10)
        } else {
            (m, pow10.checked_shl(q.unsigned_abs())?)
        };
        let quotient = num.checked_div(den)?;
        let rem = num.checked_rem(den)?;
        let twice = rem.checked_mul(2)?;
        (quotient, twice > den || (twice == den && quotient & 1 == 1))
    };
    Some((u64::try_from(quotient).ok()?, round_up))
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

    // r[verify vcf_writer.float_non_finite]
    /// The non-finite spellings, pinned against what `bcftools view` prints
    /// for the same values. `Display` would write `NaN`; C's `%g`, and so
    /// htslib, writes `nan`.
    ///
    /// NaN is deliberately *not* `.`: the missing sentinel is the signaling
    /// NaN `0x7F800001`, a quiet NaN is `0x7FC00000`, and htslib's own
    /// `bcf_float_is_missing` is an exact bit compare rather than `isnan`.
    /// Collapsing the two here would also put this arm at odds with the BCF
    /// arm, which writes the caller's bits through untouched.
    #[test]
    fn float_format_non_finite() {
        assert_eq!(formatted(f32::NAN), "nan");
        assert_eq!(formatted(-f32::NAN), "nan");
        assert_eq!(formatted(f32::INFINITY), "inf");
        assert_eq!(formatted(f32::NEG_INFINITY), "-inf");
    }

    // r[verify vcf_writer.float_non_finite]
    /// Whatever it emits for a non-finite value has to parse back to an equal
    /// value, the same round-trip the finite case owes — and `bcftools` uses
    /// `strtod`, which accepts exactly these three spellings.
    #[test]
    fn non_finite_text_parses_back() {
        assert!(formatted(f32::NAN).parse::<f32>().expect("`nan` parses").is_nan());
        assert_eq!(formatted(f32::INFINITY).parse::<f32>().expect("`inf` parses"), f32::INFINITY);
        assert_eq!(
            formatted(f32::NEG_INFINITY).parse::<f32>().expect("`-inf` parses"),
            f32::NEG_INFINITY
        );
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

    fn fmt_path(v: f32) -> String {
        let mut buf = Vec::new();
        write_float_g_fmt(&mut buf, v).expect("every finite float must be writable");
        String::from_utf8(buf).expect("ASCII")
    }

    /// C's own `%g`, which is what htslib falls back to and what
    /// `r[vcf_writer.float_precision]` names: an oracle that shares nothing
    /// with either Rust path. macOS's libc keeps the trailing zeros when a
    /// scientific-form value is an exact tie that rounds down (`1000005.0` →
    /// `1.00000e+06`, where glibc and C99 say `1e+06`), so the mantissa is
    /// normalized before comparing.
    fn c_printf_g(v: f32) -> String {
        unsafe extern "C" {
            fn snprintf(
                buf: *mut std::ffi::c_char,
                n: usize,
                fmt: *const std::ffi::c_char,
                ...
            ) -> i32;
        }
        let mut out = [0u8; 64];
        // SAFETY: the buffer is 64 bytes and `%g` of a double writes at most ~15.
        let n =
            unsafe { snprintf(out.as_mut_ptr().cast(), out.len(), c"%g".as_ptr(), f64::from(v)) };
        let n = usize::try_from(n).expect("snprintf succeeds");
        let text = String::from_utf8(out[..n].to_vec()).expect("ASCII");
        match text.split_once('e') {
            Some((mantissa, exp)) if mantissa.contains('.') => {
                format!("{}e{exp}", mantissa.trim_end_matches('0').trim_end_matches('.'))
            }
            _ => text,
        }
    }

    // r[verify vcf_writer.float_precision]
    #[hegel::test(test_cases = 2000)]
    fn exact_path_matches_c_printf_and_fmt(tc: TestCase) {
        let v = f32::from_bits(tc.draw(gs::integers::<u32>()));
        if !v.is_finite() {
            return;
        }
        let got = formatted(v);
        assert_eq!(got, fmt_path(v), "{v:e} ({:#x})", v.to_bits());
        assert_eq!(got, c_printf_g(v), "{v:e} ({:#x})", v.to_bits());
    }

    // r[verify vcf_writer.float_precision]
    /// Ties and near-ties at every decimal scale: the values rounding to six
    /// digits is most likely to get wrong, which uniform bit patterns rarely hit.
    #[test]
    fn exact_path_matches_c_printf_at_ties() {
        for exp in -40i32..=38 {
            for digits in [100_000u64, 123_455, 999_995, 999_999, 500_005, 250_000] {
                let v = (digits as f64 + 0.5) * 10f64.powi(exp - 5);
                #[expect(clippy::cast_possible_truncation, reason = "test inputs")]
                let v = v as f32;
                for w in [
                    v,
                    f32::from_bits(v.to_bits().wrapping_add(1)),
                    f32::from_bits(v.to_bits().wrapping_sub(1)),
                ] {
                    if w.is_finite() {
                        assert_eq!(formatted(w), c_printf_g(w), "{w:e}");
                    }
                }
            }
        }
        // Exactly representable ties: 0.5 steps at the sixth digit.
        for v in [1_000_000.5f32, 1_234_567.5, 999_999.5, 0.125, 1.000_005, 2.5e-5, 8_388_609.0] {
            assert_eq!(formatted(v), c_printf_g(v), "{v:e}");
        }
    }

    /// Every positive finite `f32` against the `core::fmt` path — about 2.1
    /// billion values, a few minutes across all cores. Run by hand:
    /// `cargo test --release -p seqair --lib exhaustive_f32 -- --ignored`.
    #[test]
    #[ignore = "exhaustive; minutes in release mode"]
    fn exhaustive_f32_matches_fmt_path() {
        let threads = std::thread::available_parallelism().map_or(8, |n| n.get());
        let threads = u32::try_from(threads).unwrap_or(8);
        let end = f32::INFINITY.to_bits();
        let step = end.div_ceil(threads);
        std::thread::scope(|scope| {
            for t in 0..threads {
                scope.spawn(move || {
                    let (mut a, mut b) = (Vec::new(), Vec::new());
                    for bits in t * step..((t + 1) * step).min(end) {
                        let v = f32::from_bits(bits);
                        a.clear();
                        b.clear();
                        write_float_g(&mut a, v).expect("finite");
                        write_float_g_fmt(&mut b, v).expect("finite");
                        assert_eq!(a, b, "{v:e} ({bits:#x})");
                    }
                });
            }
        });
    }

    #[test]
    fn percent_encoding_special_chars() {
        let mut buf = Vec::new();
        percent_encode_into(&mut buf, b"hello:world;key=val%done,x\ty\nz\rend");
        let result = String::from_utf8(buf).unwrap();
        assert_eq!(result, "hello%3Aworld%3Bkey%3Dval%25done%2Cx%09y%0Az%0Dend");
    }
}
