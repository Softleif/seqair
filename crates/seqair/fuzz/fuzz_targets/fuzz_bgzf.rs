#![no_main]

use libfuzzer_sys::fuzz_target;
use seqair::bam::bgzf::decode_bgzf_block;

fuzz_target!(|data: &[u8]| {
    let _ = decode_bgzf_block(data);
});
