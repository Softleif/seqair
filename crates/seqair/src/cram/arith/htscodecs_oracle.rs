//! htscodecs as the oracle for the arithmetic coder and tok3 over it: hegel
//! draws data and flags, htscodecs' own encoder (`arith_compress_to`,
//! `tok3_encode_names`, from the static htslib that the rust-htslib
//! dev-dependency links) compresses, and the production decoder and the
//! reference decoder must both return the original. Differential properties
//! feed all three decoders corrupted streams.
#![allow(
    clippy::unwrap_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "test code with bounded sizes"
)]

use std::ffi::{c_char, c_int, c_uint};

use hegel::prelude::*;

use super::{decode, reference};
use crate::cram::tok3;

// Pulls in hts-sys' static libhts, which contains htscodecs.
use rust_htslib as _;

/// htscodecs' `arith_dynamic.h` / `tokenise_name3.h`.
mod ffi {
    use std::ffi::{c_char, c_int, c_uint};

    unsafe extern "C" {
        pub(super) fn arith_compress_to(
            input: *mut u8,
            in_size: c_uint,
            out: *mut u8,
            out_size: *mut c_uint,
            order: c_int,
        ) -> *mut u8;

        pub(super) fn arith_uncompress_to(
            input: *mut u8,
            in_size: c_uint,
            out: *mut u8,
            out_size: *mut c_uint,
        ) -> *mut u8;

        pub(super) fn tok3_encode_names(
            blk: *mut c_char,
            len: c_int,
            level: c_int,
            use_arith: c_int,
            out_len: *mut c_int,
            last_start_p: *mut c_int,
        ) -> *mut u8;
    }
}

/// Take ownership of a buffer htscodecs `malloc`ed.
fn take_c_buffer(p: *mut u8, len: usize) -> Vec<u8> {
    assert!(!p.is_null(), "htscodecs returned NULL");
    // SAFETY: htscodecs returned `len` initialised bytes at `p`, allocated
    // with malloc, which nothing else frees.
    unsafe {
        let v = std::slice::from_raw_parts(p, len).to_vec();
        libc::free(p.cast());
        v
    }
}

/// `arith_compress_to` with `order` holding the flags (and the stripe count
/// in bits 8..16).
fn htscodecs_compress(data: &[u8], order: c_int) -> Vec<u8> {
    let mut input = data.to_vec();
    let mut out_size: c_uint = 0;
    // SAFETY: a null `out` makes htscodecs allocate the output and report
    // its size; `input` is live and `in_size` bytes long.
    let p = unsafe {
        ffi::arith_compress_to(
            input.as_mut_ptr(),
            c_uint::try_from(input.len()).unwrap(),
            std::ptr::null_mut(),
            &raw mut out_size,
            order,
        )
    };
    take_c_buffer(p, usize::try_from(out_size).unwrap())
}

/// `arith_uncompress_to` into a buffer of exactly `len` bytes, as tok3 and
/// the stripe path call it.
fn htscodecs_decompress(src: &[u8], len: usize) -> Option<Vec<u8>> {
    let mut input = src.to_vec();
    let mut out = vec![0u8; len];
    let mut out_size = c_uint::try_from(len).unwrap();
    // SAFETY: `out` holds `out_size` bytes, `input` holds `in_size`.
    let p = unsafe {
        ffi::arith_uncompress_to(
            input.as_mut_ptr(),
            c_uint::try_from(input.len()).unwrap(),
            out.as_mut_ptr(),
            &raw mut out_size,
        )
    };
    if p.is_null() {
        return None;
    }
    out.truncate(usize::try_from(out_size).unwrap());
    Some(out)
}

/// `tok3_encode_names` over newline-terminated names.
fn htscodecs_tokenise(names: &[Vec<u8>], level: c_int, use_arith: bool) -> Vec<u8> {
    let mut block: Vec<u8> = names.iter().flat_map(|n| n.iter().copied().chain(*b"\n")).collect();
    let mut out_len: c_int = 0;
    let mut last_start: c_int = 0;
    // SAFETY: `block` is live and `len` bytes long (htscodecs writes NULs
    // into it); the result is malloc'ed with `out_len` bytes.
    let p = unsafe {
        ffi::tok3_encode_names(
            block.as_mut_ptr().cast::<c_char>(),
            c_int::try_from(block.len()).unwrap(),
            level,
            c_int::from(use_arith),
            &raw mut out_len,
            &raw mut last_start,
        )
    };
    take_c_buffer(p, usize::try_from(out_len).unwrap())
}

const ORDER1: c_int = 0x01;
const EXT: c_int = 0x04;
const STRIPE: c_int = 0x08;
const NOSZ: c_int = 0x10;
const RLE: c_int = 0x40;
const PACK: c_int = 0x80;

/// Data in the shapes the coder meets: raw bytes, a small alphabet in runs,
/// long runs of a couple of symbols, and little-endian u32s.
#[hegel::composite]
fn arb_data(tc: &TestCase) -> Vec<u8> {
    match tc.draw_silent(gs::integers::<u8>().max_value(3)) {
        0 => tc.draw_silent(gs::binary().max_size(3000)),
        1 => {
            let k = tc.draw_silent(gs::integers::<u8>().min_value(1).max_value(24));
            let base = tc.draw_silent(gs::integers::<u8>());
            let runs = tc.draw_silent(
                gs::vecs(gs::tuples!(
                    gs::integers::<u8>().max_value(k - 1),
                    gs::integers::<usize>().max_value(12)
                ))
                .max_size(400),
            );
            runs.into_iter()
                .flat_map(|(s, n)| std::iter::repeat_n(base.wrapping_add(s), n + 1))
                .collect()
        }
        2 => {
            let runs = tc.draw_silent(
                gs::vecs(gs::tuples!(
                    gs::integers::<u8>().max_value(1),
                    gs::integers::<usize>().max_value(3000)
                ))
                .max_size(8),
            );
            runs.into_iter().flat_map(|(s, n)| std::iter::repeat_n(b'A' + s, n)).collect()
        }
        _ => {
            let max = tc.draw_silent(gs::integers::<u32>().min_value(1));
            let vals = tc.draw_silent(gs::vecs(gs::integers::<u32>().max_value(max)).max_size(800));
            vals.iter().flat_map(|v| v.to_le_bytes()).collect()
        }
    }
}

/// Any combination of the flags htscodecs' encoder round-trips, with a
/// stripe count of 1..=8. Left out: the reserved order bit (with RLE its
/// encoder writes order-1 data that its decoder reads as order-0), and CAT,
/// which the encoder picks by itself for data that does not shrink; asked
/// for, it can write a stream its own decoder rejects (bzip2 data under a
/// CAT flag with EXT, a truncated stripe, a PACK stream without its body).
#[hegel::composite]
fn arb_order(tc: &TestCase) -> c_int {
    let mut order = 0;
    for flag in [ORDER1, EXT, STRIPE, NOSZ, RLE, PACK] {
        if tc.draw_silent(gs::booleans()) {
            order |= flag;
        }
    }
    if order & STRIPE != 0 {
        order |= c_int::from(tc.draw_silent(gs::integers::<u8>().min_value(1).max_value(8))) << 8;
    }
    order
}

// r[verify cram.codec.arith+2]
// r[verify cram.codec.arith.wrapper]
// r[verify cram.codec.arith.order]
// r[verify cram.codec.arith.rle]
// r[verify cram.codec.arith.pack]
// r[verify cram.codec.arith.stripe]
// r[verify cram.codec.arith.ext]
#[hegel::test]
fn decodes_what_htscodecs_encodes(tc: TestCase) {
    let data = tc.draw(arb_data());
    let order = tc.draw(arb_order());
    let stream = htscodecs_compress(&data, order);
    // The oracle holds: htscodecs reads its own stream back.
    assert_eq!(htscodecs_decompress(&stream, data.len()).as_deref(), Some(data.as_slice()));
    // The size only matters to a NOSZ stream, which cannot say it itself.
    let decoded = decode(&stream, data.len())
        .unwrap_or_else(|e| panic!("order {order:#x}, {} bytes: {e}", data.len()));
    assert!(decoded == data, "order {order:#x}: decoded data differs");
    let decoded = reference::decode(&stream, data.len())
        .unwrap_or_else(|| panic!("order {order:#x}: the reference decoder rejects it"));
    assert!(decoded == data, "order {order:#x}: the reference decoder's data differs");
}

/// Whatever the production decoder accepts, htscodecs decodes to the same
/// bytes, and so does the reference decoder where it accepts it too (the
/// reference's module docs list where their verdicts differ). Corrupting
/// valid streams reaches the paths random bytes never get past.
///
/// htscodecs is not asked about STRIPE streams: it hands a substream every
/// byte to the end of the input, so a corrupt one can read on into the next
/// substream and fail (or succeed) where seqair, which stops at the
/// substream's own length, does the opposite. Valid streams never do that.
#[hegel::test(test_cases = 1000)]
fn agrees_with_htscodecs_on_corrupted_streams(tc: TestCase) {
    let data = tc.draw(arb_data());
    let order = tc.draw(arb_order());
    let mut stream = htscodecs_compress(&data, order);
    let flips = tc.draw(
        gs::vecs(gs::tuples!(gs::integers::<usize>(), gs::integers::<u8>().min_value(1)))
            .max_size(3),
    );
    for (at, x) in flips {
        let len = stream.len();
        stream[at % len] ^= x;
    }
    let Ok(ours) = decode(&stream, data.len()) else { return };
    if let Some(theirs) = reference::decode(&stream, data.len()) {
        assert!(ours == theirs, "order {order:#x}: the reference decoder disagrees");
    }
    if order & STRIPE == 0 {
        let theirs = htscodecs_decompress(&stream, ours.len());
        assert_eq!(theirs.as_deref(), Some(ours.as_slice()), "order {order:#x}");
    }
}

/// A read name: Illumina-like `PREFIX:lane:tile:x:y#idx/pair`, a counter
/// with a fixed prefix, or random printable ASCII.
#[hegel::composite]
fn arb_name(tc: &TestCase) -> Vec<u8> {
    let num = |max: u32| tc.draw_silent(gs::integers::<u32>().max_value(max));
    match tc.draw_silent(gs::integers::<u8>().max_value(2)) {
        0 => {
            let prefix =
                tc.draw_silent(gs::sampled_from(&["READ", "I17_08765", "SRR062634.1"][..]));
            let fields =
                format!("{prefix}:{}:{}:{}:{}", num(8), num(3000), num(99_999), num(99_999));
            let suffix = if tc.draw_silent(gs::booleans()) {
                format!("#{}/{}", num(9), num(1) + 1)
            } else {
                String::new()
            };
            format!("{fields}{suffix}").into_bytes()
        }
        1 => format!("m{:08}_{}", num(200), num(1_000_000)).into_bytes(),
        _ => tc.draw_silent(
            gs::vecs(gs::integers::<u8>().min_value(b' ').max_value(b'~')).max_size(40),
        ),
    }
}

// r[verify cram.codec.tok3_arith]
// r[verify cram.codec.tok3]
#[hegel::test(test_cases = 200)]
fn decodes_what_htscodecs_tokenises(tc: TestCase) {
    let names = tc.draw(gs::vecs(arb_name()).min_size(1).max_size(300));
    let level = tc.draw(gs::integers::<c_int>().min_value(1).max_value(9));
    let use_arith = tc.draw(gs::booleans());

    let encoded = htscodecs_tokenise(&names, level, use_arith);
    assert_eq!(encoded[8] != 0, use_arith, "tok3 header records the entropy coder");

    let expected: Vec<u8> = names.iter().flat_map(|n| n.iter().copied().chain([0])).collect();
    let decoded = tok3::decode(&encoded)
        .unwrap_or_else(|e| panic!("level {level}, use_arith {use_arith}: {e}"));
    assert!(decoded == expected, "level {level}, use_arith {use_arith}: names differ");
    let decoded = tok3::decode_with_arith_reference(&encoded)
        .unwrap_or_else(|e| panic!("reference, level {level}, use_arith {use_arith}: {e}"));
    assert!(decoded == expected, "reference, level {level}, use_arith {use_arith}: names differ");
}

/// tok3 over the production and the reference arithmetic decoder agree on
/// whatever both accept, corrupted streams included.
// r[verify cram.codec.tok3_arith]
#[hegel::test(test_cases = 500)]
fn tok3_agrees_with_the_reference_on_corrupted_blocks(tc: TestCase) {
    let names = tc.draw(gs::vecs(arb_name()).min_size(1).max_size(50));
    let level = tc.draw(gs::integers::<c_int>().min_value(1).max_value(9));
    let mut encoded = htscodecs_tokenise(&names, level, true);
    let flips = tc.draw(
        gs::vecs(gs::tuples!(gs::integers::<usize>(), gs::integers::<u8>().min_value(1)))
            .max_size(3),
    );
    for (at, x) in flips {
        // Past the header, which only sizes the output.
        let len = encoded.len() - 9;
        encoded[9 + at % len] ^= x;
    }
    if let (Ok(ours), Ok(theirs)) =
        (tok3::decode(&encoded), tok3::decode_with_arith_reference(&encoded))
    {
        assert!(ours == theirs, "tok3 over the two arith decoders disagrees");
    }
}
