#![no_main]

//! fqzcomp (CRAM method 7), differentially: the production decoder and the
//! spec-literal reference share no code, so they must agree on every input —
//! the same output, or both an error.

use libfuzzer_sys::fuzz_target;
use seqair::cram::{fqzcomp, fqzcomp_reference};

/// Streams claiming more output than this are skipped: a few input bytes can
/// legitimately expand to the full claimed size (duplicate records cost
/// almost nothing), and decoding hundreds of MiB per input starves the fuzzer.
const MAX_CLAIMED_OUTPUT: u64 = 1 << 20;

/// The uint7 output size at the start of the stream, if it has one.
fn claimed_output(data: &[u8]) -> Option<u64> {
    let mut n = 0u64;
    for &b in data.iter().take(5) {
        n = (n << 7) | u64::from(b & 0x7f);
        if b & 0x80 == 0 {
            return Some(n);
        }
    }
    None
}

fuzz_target!(|data: &[u8]| {
    if claimed_output(data).is_some_and(|n| n > MAX_CLAIMED_OUTPUT) {
        return;
    }
    let production = fqzcomp::decode(data);
    let reference = fqzcomp_reference::decode(data);
    match (&production, &reference) {
        (Ok(p), Some(r)) => assert!(p == r, "decoders disagree on the output"),
        (Err(_), None) => {}
        (Ok(_), None) => panic!("production decodes, reference rejects"),
        (Err(e), Some(_)) => panic!("reference decodes, production rejects: {e}"),
    }
});
