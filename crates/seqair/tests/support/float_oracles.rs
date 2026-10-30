//! Floats for the SAM float printers' tests: generators reaching every corner
//! of `f64` and `f32`, and the C printers they are checked against (libc's
//! `%g`, htslib's `kputd`) with macOS's departures from glibc mapped out.
//!
//! Shared by `src/sam/format.rs`'s unit tests and, through `raw_record_gen.rs`,
//! the integration tests; it names no crate, so it compiles in both.
#![allow(
    dead_code,
    clippy::cast_possible_truncation,
    clippy::indexing_slicing,
    clippy::unwrap_used,
    clippy::arithmetic_side_effects,
    reason = "test support; each including test uses a different part"
)]

use hegel::TestCase;
use hegel::generators as gs;
use rust_htslib::htslib;

fn n<T: hegel::generators::Integer>(lo: T, hi: T) -> gs::IntegerGenerator<T> {
    gs::integers::<T>().min_value(lo).max_value(hi)
}

fn special_f64() -> Vec<f64> {
    vec![
        0.0,
        -0.0,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NAN,
        -f64::NAN,
        f64::MIN_POSITIVE,
        f64::from_bits(1),
        f64::MAX,
        f64::MIN,
        f64::EPSILON,
        1e-5,
        0.0001,
        0.000_099_999_99,
        0.001,
        0.009_999_999,
        0.099_999_95,
        0.999_999_5,
        1.0,
        9.999_995,
        99_999.95,
        99_999.97,
        123_456.5,
        999_998.5,
        999_999.0,
        999_999.4,
        999_999.5,
        1e6,
        1_234_565.0,
        1e15,
        1e16,
        1e100,
        1e-100,
        0.5,
        0.25,
        2.5e-5,
        -123.456,
        f64::from(f32::MAX),
        f64::from(f32::MIN_POSITIVE),
        f64::from(f32::from_bits(1)),
    ]
}

/// A double from every corner: raw bit patterns, the special list, and
/// decimals with seven significant digits (near-ties for six-digit rounding).
pub fn arb_f64(tc: &TestCase, negative_nan: bool) -> f64 {
    let v = match tc.draw_silent(n(0u8, 3)) {
        0 => f64::from_bits(tc.draw_silent(gs::integers::<u64>())),
        1 => tc.draw_silent(gs::sampled_from(special_f64())),
        2 => {
            let mantissa = tc.draw_silent(n(-10_000_000i64, 10_000_000));
            let exp = tc.draw_silent(n(-15i32, 15));
            mantissa as f64 * 10f64.powi(exp)
        }
        _ => tc.draw_silent(gs::floats::<f64>()),
    };
    if v.is_nan() && !negative_nan { f64::from_bits(v.to_bits() & !(1 << 63)) } else { v }
}

/// A float from the same corners, in `f32`.
pub fn arb_f32(tc: &TestCase, negative_nan: bool) -> f32 {
    let v = match tc.draw_silent(n(0u8, 3)) {
        0 => f32::from_bits(tc.draw_silent(gs::integers::<u32>())),
        1 => tc.draw_silent(gs::sampled_from(special_f64())) as f32,
        2 => {
            let mantissa = tc.draw_silent(n(-10_000_000i64, 10_000_000));
            let exp = tc.draw_silent(n(-15i32, 15));
            (mantissa as f64 * 10f64.powi(exp)) as f32
        }
        _ => tc.draw_silent(gs::floats::<f32>()),
    };
    if v.is_nan() && !negative_nan { f32::from_bits(v.to_bits() & !(1 << 31)) } else { v }
}

pub fn libc_g(value: f64) -> Vec<u8> {
    let mut buf = [0u8; 64];
    // SAFETY: `buf` is writable for its length and the format takes one double.
    let n = unsafe { libc::snprintf(buf.as_mut_ptr().cast(), buf.len(), c"%g".as_ptr(), value) };
    buf[..usize::try_from(n).unwrap()].to_vec()
}

/// `%g` output with the trailing zeros of an exponent form's mantissa
/// removed, as the C standard requires. macOS's libc keeps them when the value
/// is an exact tie rounded down (`1000005` prints `1.00000e+06`, glibc `1e+06`).
pub fn strip_mantissa_zeros(g: &[u8]) -> Vec<u8> {
    let Some(e) = g.iter().position(|&b| b == b'e') else { return g.to_vec() };
    let (mantissa, exponent) = g.split_at(e);
    if !mantissa.contains(&b'.') {
        return g.to_vec();
    }
    let mut out = mantissa.to_vec();
    while out.last() == Some(&b'0') {
        out.pop();
    }
    if out.last() == Some(&b'.') {
        out.pop();
    }
    out.extend_from_slice(exponent);
    out
}

/// What `format_g` must print, given libc's output for `value`: the same,
/// except where macOS's libc departs from glibc and musl, which seqair
/// follows: a negative NaN prints `nan` there (`-nan` here), and the exact-tie
/// trailing zeros of [`strip_mantissa_zeros`].
pub fn expect_from_libc(value: f64, libc: Vec<u8>) -> Vec<u8> {
    if !cfg!(target_os = "macos") {
        return libc;
    }
    if value.is_nan() && value.is_sign_negative() {
        assert_eq!(libc, b"nan");
        return b"-nan".to_vec();
    }
    strip_mantissa_zeros(&libc)
}

pub fn kstring_bytes(ks: &mut htslib::kstring_t) -> Vec<u8> {
    let out = if ks.s.is_null() {
        Vec::new()
    } else {
        // SAFETY: htslib wrote `l` bytes at `s`.
        unsafe { std::slice::from_raw_parts(ks.s.cast::<u8>(), ks.l) }.to_vec()
    };
    // SAFETY: htslib allocates kstrings with malloc.
    unsafe { libc::free(ks.s.cast()) };
    ks.s = std::ptr::null_mut();
    out
}

pub fn htslib_kputd(value: f64) -> Vec<u8> {
    let mut ks = htslib::kstring_t { l: 0, m: 0, s: std::ptr::null_mut() };
    // SAFETY: an empty kstring is valid input; kputd grows it.
    let r = unsafe { htslib::kputd(value, &mut ks) };
    assert!(r >= 0, "kputd failed");
    kstring_bytes(&mut ks)
}
