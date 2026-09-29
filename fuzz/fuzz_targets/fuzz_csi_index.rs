#![no_main]

use libfuzzer_sys::fuzz_target;
use seqair::bam::csi_index::CsiIndex;

fuzz_target!(|data: &[u8]| {
    // Fuzz CSI index parsing from raw bytes (both uncompressed and BGZF-compressed)
    let _ = CsiIndex::from_bytes(data);
});
