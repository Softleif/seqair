#![no_main]

//! tok3 (CRAM method 8), differentially: the production decoder and the
//! spec-literal reference share no tok3 code, so over the same arithmetic
//! decoder they must agree on every input — the same names, or both an
//! error. Checked over the production arithmetic decoder and over the
//! reference one.

use libfuzzer_sys::fuzz_target;
use seqair::cram::{arith, tok3, tok3_reference};

/// Blocks claiming more names or output than these are skipped: the
/// reference materialises a regenerated type stream per position (up to 128
/// of the name count each) and is slow on large outputs. htscodecs caps its
/// fuzzing builds at 100 000 names.
const MAX_NAMES: u32 = 100_000;
const MAX_CLAIMED_OUTPUT: u32 = 1_000_000;

fn header_u32(data: &[u8], at: usize) -> Option<u32> {
    Some(u32::from_le_bytes(data.get(at..at + 4)?.try_into().ok()?))
}

fn assert_same(what: &str, production: Option<Vec<u8>>, reference: Option<Vec<u8>>) {
    match (&production, &reference) {
        (Some(p), Some(r)) => assert!(p == r, "{what}: decoders disagree on the names"),
        (None, None) => {}
        (Some(_), None) => panic!("{what}: production decodes, reference rejects"),
        (None, Some(_)) => panic!("{what}: reference decodes, production rejects"),
    }
}

fuzz_target!(|data: &[u8]| {
    if header_u32(data, 0).is_some_and(|n| n > MAX_CLAIMED_OUTPUT)
        || header_u32(data, 4).is_some_and(|n| n > MAX_NAMES)
    {
        return;
    }
    assert_same(
        "production arith",
        tok3::decode(data).ok(),
        tok3_reference::decode_with(data, |s| arith::decode(s, 0).ok()),
    );
    assert_same(
        "reference arith",
        tok3::decode_with_arith_reference(data).ok(),
        tok3_reference::decode(data),
    );
});
