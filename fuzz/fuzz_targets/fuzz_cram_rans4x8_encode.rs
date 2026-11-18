#![no_main]

//! The CRAM writer's rANS 4×8 encoder, differentially: the fast path (SIMD
//! kernels) must write the reference encoder's bytes
//! (`r[cram_writer.rans4x8.simd]`), and the decoder must read them back, at
//! the order the first byte's low bit names.

use libfuzzer_sys::fuzz_target;
use seqair::cram::rans::decode;
use seqair::cram::writer::rans4x8_encode;

fuzz_target!(|input: &[u8]| {
    let Some((&method, data)) = input.split_first() else { return };
    let order_1 = method & 1 != 0;
    let fast = rans4x8_encode(data, order_1, false);
    let reference = rans4x8_encode(data, order_1, true);
    assert!(fast == reference, "order {}: fast and reference encoders differ", u8::from(order_1));
    if let Some(coded) = fast {
        let decoded = decode(&coded).expect("the encoder's stream decodes");
        assert!(decoded == data, "order {}: round trip differs", u8::from(order_1));
    }
});
