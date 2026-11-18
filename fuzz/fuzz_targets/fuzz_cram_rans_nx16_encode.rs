#![no_main]

//! The CRAM writer's rANS Nx16 encoder, differentially: the SIMD kernels
//! must write the scalar encoder's bytes (`r[cram_writer.rans_nx16.simd]`),
//! and the decoder must read them back, for every combination of order,
//! states, STRIPE, RLE and PACK the first byte names.

use libfuzzer_sys::fuzz_target;
use seqair::cram::rans_nx16::decode;
use seqair::cram::writer::rans_nx16_encode;

fuzz_target!(|input: &[u8]| {
    let Some((&method, data)) = input.split_first() else { return };
    // Order, X32, STRIPE, RLE, PACK: the bits `Nx16Flags::from_byte` reads.
    let method = method & 0xCD;
    let simd = rans_nx16_encode(data, method, false);
    let scalar = rans_nx16_encode(data, method, true);
    assert!(simd == scalar, "method {method:#04x}: SIMD and scalar encoders differ");
    if let Some(coded) = simd {
        let decoded = decode(&coded, data.len()).expect("the encoder's stream decodes");
        assert!(decoded == data, "method {method:#04x}: round trip differs");
    }
});
