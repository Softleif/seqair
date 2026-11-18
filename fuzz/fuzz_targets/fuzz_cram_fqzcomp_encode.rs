#![no_main]

//! The CRAM fqzcomp encoder (method 7): whatever it writes, under any of its
//! strategies, the decoder gives back the qualities.
//!
//! Input: a strategy byte, then records, each a length byte (0 makes an
//! empty record), a flags byte (bit 0 reverse,
//! bit 1 READ2, bit 2 multiplies the length by 8: long records reach the position
//! context's saturation) and that many quality bytes, cut to what is left;
//! a lone byte at the end is dropped. A record that copies its predecessor
//! exercises dedup.

use libfuzzer_sys::fuzz_target;
use seqair::cram::fqzcomp;
use seqair::cram::writer::fuzz::fqzcomp_encode;

fuzz_target!(|input: &[u8]| {
    let Some((&strategy, mut rest)) = input.split_first() else { return };
    let mut quals = Vec::new();
    let mut records = Vec::new();
    while let [len, flags, tail @ ..] = rest {
        let mut len = usize::from(*len);
        if flags & 4 != 0 {
            len *= 8;
        }
        let take = len.min(tail.len());
        let (q, after) = tail.split_at(take);
        quals.extend_from_slice(q);
        records.push((u32::try_from(take).unwrap_or(u32::MAX), flags & 1 != 0, flags & 2 != 0));
        rest = after;
    }
    let stream = fqzcomp_encode(&quals, &records, strategy).expect("the records add up");
    let decoded = fqzcomp::decode(&stream).expect("the decoder reads what the encoder wrote");
    assert!(decoded == quals, "strategy {}: round trip differs", strategy % 4);
});
