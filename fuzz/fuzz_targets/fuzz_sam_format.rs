#![no_main]

//! Arbitrary record payloads: whatever `RawBamRecord::from_payload` accepts,
//! `format_sam_record` must format or reject with a typed error — never panic,
//! and never leave a partial line behind.

use libfuzzer_sys::fuzz_target;
use seqair::bam::{BamHeader, RawBamRecord};
use seqair::sam::format::format_sam_record;
use std::sync::OnceLock;

fn header() -> &'static BamHeader {
    static HEADER: OnceLock<BamHeader> = OnceLock::new();
    HEADER.get_or_init(|| {
        BamHeader::from_sam_text("@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:2000\n")
            .expect("valid header")
    })
}

fuzz_target!(|data: &[u8]| {
    let Ok(record) = RawBamRecord::from_payload(data.to_vec(), header()) else { return };
    let mut out = b"prefix".to_vec();
    match format_sam_record(&record, header(), &mut out) {
        Ok(()) => assert!(out.starts_with(b"prefix")),
        Err(_) => assert_eq!(out, b"prefix", "a failed format must not leave bytes behind"),
    }
});
