#![no_main]

use libfuzzer_sys::fuzz_target;
use seqair::cram::block::parse_block;

fuzz_target!(|data: &[u8]| {
    let _ = parse_block(data);
});
