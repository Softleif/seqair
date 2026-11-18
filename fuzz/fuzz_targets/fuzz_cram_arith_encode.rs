#![no_main]

//! The CRAM arithmetic coder's encoder (method 6): whatever it writes, under
//! any combination of order, STRIPE, RLE and PACK, the decoder gives back
//! the input. The first byte picks the combination, the rest is the data.

use libfuzzer_sys::fuzz_target;
use seqair::cram::arith;
use seqair::cram::writer::fuzz::arith_encode;

fuzz_target!(|input: &[u8]| {
    let Some((&flags, data)) = input.split_first() else { return };
    let stream = arith_encode(data, flags).expect("a fuzz input fits the format's sizes");
    let decoded = arith::decode(&stream, 0).expect("the decoder reads what the encoder wrote");
    assert!(decoded == data, "flags {flags:#04x}: round trip differs");
});
