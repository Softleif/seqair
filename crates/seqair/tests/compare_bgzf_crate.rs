//! Cross-validation: compress with the `bgzf` crate's Writer (which handles
//! and seqair's `RegionBuf`. Both must match the original byte-for-byte.
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
use bgzf::{CompressionLevel, Reader as BgzfReader, Writer as BgzfWriter};
use proptest::prelude::*;
use seqair::bam::{bgzf::VirtualOffset, index::Chunk, region_buf::RegionBuf};
use std::io::{Read, Write};

/// Compress data using `bgzf::Writer` (handles block splitting via Write trait).
fn bgzf_compress(data: &[u8], level: u8) -> Vec<u8> {
    let mut output = Vec::new();
    let mut writer = BgzfWriter::new(&mut output, CompressionLevel::new(level).unwrap());
    writer.write_all(data).unwrap();
    writer.finish().unwrap();
    output
}

/// Decompress using `bgzf::Reader`.
fn bgzf_decompress(compressed: &[u8]) -> Vec<u8> {
    let mut reader = BgzfReader::new(compressed);
    let mut output = Vec::new();
    let mut buf = [0u8; 4096];
    loop {
        match reader.read(&mut buf) {
            Ok(0) => break,
            Ok(n) => output.extend_from_slice(&buf[..n]),
            Err(_) => break,
        }
    }
    output
}

/// Decompress using seqair's `RegionBuf`.
/// `expected_len` is used to pre-allocate and read the exact amount.
fn seqair_decompress(compressed: &[u8], expected_len: usize) -> Vec<u8> {
    let mut cursor = std::io::Cursor::new(compressed.to_vec());
    let file_len = compressed.len() as u64;
    let chunks =
        vec![Chunk { begin: VirtualOffset::new(0, 0), end: VirtualOffset::new(file_len, 0) }];

    let mut buf = RegionBuf::load(&mut cursor, &chunks).unwrap();
    buf.seek_virtual(VirtualOffset::new(0, 0)).unwrap();

    let mut output = vec![0u8; expected_len];
    buf.read_exact_into(&mut output).unwrap();
    output
}

// r[verify bgzf.decompression]
// r[verify bgzf.crc32]
// r[verify region_buf.decompress]
#[test]
fn small_data_roundtrip() {
    let data = b"Hello, BGZF world! Testing cross-crate decompression compatibility.";
    let compressed = bgzf_compress(data, 6);

    let bgzf_out = bgzf_decompress(&compressed);
    let out = seqair_decompress(&compressed, data.len());

    assert_eq!(bgzf_out, data.to_vec(), "bgzf crate mismatch");
    assert_eq!(out, data.to_vec(), "seqair mismatch");
    assert_eq!(bgzf_out, out, "cross-crate mismatch");
}

#[test]
fn multi_block_roundtrip() {
    // 100KB — forces multiple BGZF blocks (each max ~64KB uncompressed)
    let data: Vec<u8> = (0..100_000).map(|i| (i % 256) as u8).collect();
    let compressed = bgzf_compress(&data, 3);

    let bgzf_out = bgzf_decompress(&compressed);
    let out = seqair_decompress(&compressed, data.len());

    assert_eq!(bgzf_out.len(), data.len(), "bgzf length");
    assert_eq!(out.len(), data.len(), "seqair length");
    assert_eq!(bgzf_out, data, "bgzf crate mismatch");
    assert_eq!(out, data, "seqair mismatch");
    assert_eq!(bgzf_out, out, "cross-crate mismatch");
}

#[test]
fn all_compression_levels() {
    let data: Vec<u8> = (0..50_000).map(|i| (i % 256) as u8).collect();
    for level in 1..=9 {
        let compressed = bgzf_compress(&data, level);

        let bgzf_out = bgzf_decompress(&compressed);
        let out = seqair_decompress(&compressed, data.len());

        assert_eq!(bgzf_out, data, "bgzf failed at level {level}");
        assert_eq!(out, data, "seqair failed at level {level}");
    }
}

#[test]
fn highly_repetitive_data() {
    // All zeros — maximally compressible
    let data = vec![0u8; 200_000];
    let compressed = bgzf_compress(&data, 6);

    let bgzf_out = bgzf_decompress(&compressed);
    let out = seqair_decompress(&compressed, data.len());

    assert_eq!(bgzf_out, data);
    assert_eq!(out, data);
}

#[test]
fn random_looking_data() {
    // Pseudorandom — minimally compressible
    let data: Vec<u8> = (0..50_000)
        .map(|i| ((i as u64).wrapping_mul(6364136223846793005).wrapping_add(1) >> 33) as u8)
        .collect();
    let compressed = bgzf_compress(&data, 6);

    let bgzf_out = bgzf_decompress(&compressed);
    let out = seqair_decompress(&compressed, data.len());

    assert_eq!(bgzf_out, data);
    assert_eq!(out, data);
}

proptest! {
    // this is slow!
    #![proptest_config(ProptestConfig { cases: 20, .. ProptestConfig::default() })]

    /// Both decompressors must produce identical output for arbitrary data
    /// compressed by the bgzf crate's Writer.
    // r[verify bgzf.decompression]
    // r[verify bgzf.crc32]
    // r[verify region_buf.decompress]
    #[test]
    fn proptest_both_decompressors_match(
        data in prop::collection::vec(0..255u8, 1..1_000_000),
        level in 1u8..10,
    ) {
        let compressed = bgzf_compress(&data, level);

        let bgzf_out = bgzf_decompress(&compressed);
        let out = seqair_decompress(&compressed, data.len());

        prop_assert_eq!(bgzf_out.len(), data.len());
        prop_assert_eq!(out.len(), data.len());
        prop_assert_eq!(&bgzf_out, &data);
        prop_assert_eq!(&out, &data);
        prop_assert_eq!(&bgzf_out, &out);
    }
}
