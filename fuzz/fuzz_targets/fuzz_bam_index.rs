#![no_main]

use libfuzzer_sys::fuzz_target;
use seqair::bam::index::BamIndex;

fuzz_target!(|data: &[u8]| {
    // Fuzz BAI index parsing from raw bytes (with magic prefix)
    let _ = BamIndex::from_bytes(data);
});
