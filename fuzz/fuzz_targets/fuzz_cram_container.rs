#![no_main]

use libfuzzer_sys::fuzz_target;
use seqair::cram::container::ContainerHeader;

fuzz_target!(|data: &[u8]| {
    let _ = ContainerHeader::parse(data);
});
