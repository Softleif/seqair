#![no_main]

use libfuzzer_sys::fuzz_target;
use seqair::bam::header::BamHeader;

fuzz_target!(|data: &[u8]| {
    if let Ok(text) = std::str::from_utf8(data) {
        let _ = BamHeader::from_sam_text(text);
    }
});
