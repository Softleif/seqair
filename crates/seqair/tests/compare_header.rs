//! Comparison tests: BAM header parsing.
//! These tests verify that seqair's header parser produces identical
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    reason = "test code"
)]
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    reason = "test code with known small values"
)]

use rust_htslib::bam::{self, Read as _};
use std::path::Path;

fn test_bam_path() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.bam"))
}

// r[verify bam.header.target_names]
// r[verify bam.header.target_len]
// r[verify bam.header.tid_lookup]
// r[verify bam.header.references]
// r[verify bam.header.magic]
// r[verify bam.header.text]
// r[verify bam.reader.open]
// r[verify bam.reader.header_access]
// r[verify bgzf.magic]
// r[verify bgzf.bsize]
// r[verify bgzf.decompression]
// r[verify bgzf.read_exact]
// r[verify bgzf.libdeflate]
#[test]
fn header_matches_htslib() {
    let bam_path = test_bam_path();

    // --- htslib reference ---
    let hts_reader = bam::IndexedReader::from_path(bam_path).expect("htslib open");
    let hts_header = hts_reader.header();
    let hts_target_count = hts_header.target_count();

    let hts_names: Vec<String> =
        hts_header.target_names().iter().map(|n| String::from_utf8_lossy(n).to_string()).collect();

    let hts_lens: Vec<u64> = (0..hts_target_count)
        .map(|tid| hts_header.target_len(tid).expect("htslib target_len"))
        .collect();

    // Verify tid lookup roundtrip with htslib
    for (tid, name) in hts_names.iter().enumerate() {
        let resolved = hts_header.tid(name.as_bytes()).expect("htslib tid lookup");
        assert_eq!(resolved, tid as u32, "htslib tid roundtrip for {name}");
    }

    // --- seqair ---
    let header = seqair::bam::BamHeader::from_bam_path(bam_path).expect("seqair open");

    // Target count
    assert_eq!(header.target_count(), hts_target_count as usize, "target count mismatch");

    // Target names match
    let names: Vec<&str> = header.target_names().collect();
    assert_eq!(names.len(), hts_names.len());
    for (i, (name, hts_name)) in names.iter().zip(hts_names.iter()).enumerate() {
        assert_eq!(*name, hts_name, "target name mismatch at tid {i}");
    }

    // Target lengths match
    for (tid, hts_len) in hts_lens.iter().enumerate() {
        let len = header
            .target_len(tid as u32)
            .unwrap_or_else(|| panic!("seqair target_len missing for tid {tid}"));
        assert_eq!(len, *hts_len, "target length mismatch at tid {tid}");
    }

    // tid lookup roundtrip
    for (tid, name) in hts_names.iter().enumerate() {
        let resolved =
            header.tid(name).unwrap_or_else(|| panic!("seqair tid lookup failed for {name}"));
        assert_eq!(resolved, tid as u32, "tid roundtrip mismatch for {name}");
    }
}
