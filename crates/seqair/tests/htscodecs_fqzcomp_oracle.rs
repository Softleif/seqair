//! fqzcomp against htscodecs, the C implementation htslib writes CRAM 3.1 with.
//!
//! htscodecs is linked through `hts-sys` (the `rust-htslib` dev-dependency).
//! Generated records are compressed with its `fqz_compress` — with each of its
//! built-in strategies, and with generated parameter sets that reach what the
//! strategies never pick (several parameter blocks, custom tables, qtab,
//! selector tables) — and our decoder must give back the input.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    reason = "test code: indices and values bounded by the generators"
)]

#[path = "support/htscodecs_fqz.rs"]
mod htscodecs_fqz;

use hegel::TestCase;
use hegel::generators as gs;
use htscodecs_fqz::{FREVERSE, Records, builtin_stream, custom_stream, htscodecs_compress};
use seqair::cram::fqzcomp;

// ── Properties ───────────────────────────────────────────────────────

// r[verify cram.codec.fqzcomp]
// r[verify cram.codec.fqzcomp.params]
// r[verify cram.codec.fqzcomp.record]
// r[verify cram.codec.fqzcomp.context]
// r[verify cram.codec.fqzcomp.reverse]
// r[verify cram.codec.fqzcomp.models]
#[hegel::test(test_cases = 300)]
fn decodes_every_builtin_strategy(tc: TestCase) {
    let stream = builtin_stream(&tc);
    tc.assume(stream.is_some());
    let (records, stream) = stream.unwrap();
    assert_eq!(fqzcomp::decode(&stream).unwrap(), records.quals);
}

// r[verify cram.codec.fqzcomp]
// r[verify cram.codec.fqzcomp.params]
// r[verify cram.codec.fqzcomp.param_block]
// r[verify cram.codec.fqzcomp.array]
// r[verify cram.codec.fqzcomp.selector]
// r[verify cram.codec.fqzcomp.record]
// r[verify cram.codec.fqzcomp.context]
// r[verify cram.codec.fqzcomp.reverse]
#[hegel::test(test_cases = 300)]
fn decodes_generated_parameter_sets(tc: TestCase) {
    let stream = custom_stream(&tc);
    tc.assume(stream.is_some());
    let (records, stream) = stream.unwrap();
    assert_eq!(fqzcomp::decode(&stream).unwrap(), records.quals);
}

#[test]
fn record_lengths_past_two_bytes() {
    // 70 000 qualities need the third length byte's model.
    let quals: Vec<u8> = (0..70_000u32).map(|i| (i * 7 % 41) as u8).collect();
    let records = Records {
        lens: vec![5, 70_000 - 5],
        flags: vec![0, FREVERSE],
        quals,
        alphabet: (0..41).collect(),
    };
    for strat in 0..=4 {
        let compressed = htscodecs_compress(&records, 3, strat, None).unwrap();
        assert_eq!(fqzcomp::decode(&compressed).unwrap(), records.quals, "strategy {strat}");
    }
}

#[hegel::test(test_cases = 100)]
fn truncated_streams_are_errors(tc: TestCase) {
    let stream = builtin_stream(&tc);
    tc.assume(stream.is_some());
    let (_, compressed) = stream.unwrap();
    let cut = tc.draw(gs::integers::<usize>().max_value(compressed.len() - 1));
    assert!(fqzcomp::decode(&compressed[..cut]).is_err());
}

#[hegel::test(test_cases = 300)]
fn corrupted_streams_do_not_panic(tc: TestCase) {
    let stream = if tc.draw(gs::booleans()) { builtin_stream(&tc) } else { custom_stream(&tc) };
    tc.assume(stream.is_some());
    let (_, mut bytes) = stream.unwrap();
    // Mostly the parameter block, where each byte changes the layout.
    let flips: Vec<(usize, u8)> = tc.draw(
        gs::vecs(gs::tuples!(gs::integers::<usize>().max_value(80), gs::integers::<u8>()))
            .min_size(1)
            .max_size(4),
    );
    for (at, x) in flips {
        let at = at.min(bytes.len() - 1);
        bytes[at] ^= x;
    }
    if let Ok(out) = fqzcomp::decode(&bytes) {
        // Whatever it decodes, it decodes to the size the stream claims.
        let mut cur = bytes.as_slice();
        let mut claimed = 0usize;
        for _ in 0..5 {
            let b = cur[0];
            cur = &cur[1..];
            claimed = (claimed << 7) | usize::from(b & 0x7f);
            if b & 0x80 == 0 {
                break;
            }
        }
        assert_eq!(out.len(), claimed);
    }
}

#[hegel::test(test_cases = 500)]
fn arbitrary_bytes_do_not_panic(tc: TestCase) {
    let bytes = tc.draw(gs::binary().max_size(512));
    let _ = fqzcomp::decode(&bytes);
}
