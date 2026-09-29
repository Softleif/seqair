#![no_main]

use libfuzzer_sys::fuzz_target;
use seqair::cram::varint;

fuzz_target!(|data: &[u8]| {
    let _ = varint::decode_itf8(data);
    let _ = varint::decode_ltf8(data);

    // Also fuzz the cursor-advancing variants
    let mut cursor = data;
    let _ = varint::read_itf8_from(&mut cursor);

    let mut cursor = data;
    let _ = varint::read_ltf8_from(&mut cursor);
});
