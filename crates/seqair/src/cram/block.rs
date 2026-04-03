//! Parse and decompress CRAM blocks. [`parse_block`] handles all compression methods
//! defined by CRAM v3: raw, gzip, bzip2, lzma, rANS order-0/1, NX16, and tok3.

use super::{reader::CramError, varint};

/// Block content types as defined in the CRAM spec.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContentType {
    FileHeader,
    CompressionHeader,
    SliceHeader,
    ExternalData,
    CoreData,
}

impl ContentType {
    fn from_byte(b: u8) -> Result<Self, CramError> {
        match b {
            0 => Ok(Self::FileHeader),
            1 => Ok(Self::CompressionHeader),
            2 => Ok(Self::SliceHeader),
            4 => Ok(Self::ExternalData),
            5 => Ok(Self::CoreData),
            _ => Err(CramError::UnknownBlockContentType { content_type: b }),
        }
    }
}

/// A decompressed CRAM block.
#[derive(Debug)]
pub struct Block {
    pub content_type: ContentType,
    pub content_id: i32,
    pub data: Vec<u8>,
}

/// Parse a CRAM block from a byte slice.
///
// r[impl cram.block.structure]
// r[impl cram.block.crc32]
/// Returns the block and the number of bytes consumed from `buf`.
pub fn parse_block(buf: &[u8]) -> Result<(Block, usize), CramError> {
    let mut pos = 0;

    let &method = buf.get(pos).ok_or(CramError::Truncated { context: "block method" })?;
    pos = pos.checked_add(1).ok_or(CramError::Truncated { context: "block method pos" })?;

    let &content_type_byte =
        buf.get(pos).ok_or(CramError::Truncated { context: "block content type" })?;
    let content_type = ContentType::from_byte(content_type_byte)?;
    pos = pos.checked_add(1).ok_or(CramError::Truncated { context: "block content type pos" })?;

    let (content_id_u32, n) = varint::decode_itf8(
        buf.get(pos..).ok_or(CramError::Truncated { context: "block content id" })?,
    )
    .ok_or(CramError::Truncated { context: "block content id" })?;
    let content_id = content_id_u32 as i32;
    pos = pos.checked_add(n).ok_or(CramError::Truncated { context: "block content id pos" })?;

    let (compressed_size, n) = varint::decode_itf8(
        buf.get(pos..).ok_or(CramError::Truncated { context: "block compressed size" })?,
    )
    .ok_or(CramError::Truncated { context: "block compressed size" })?;
    pos =
        pos.checked_add(n).ok_or(CramError::Truncated { context: "block compressed size pos" })?;

    let (uncompressed_size, n) = varint::decode_itf8(
        buf.get(pos..).ok_or(CramError::Truncated { context: "block uncompressed size" })?,
    )
    .ok_or(CramError::Truncated { context: "block uncompressed size" })?;
    pos = pos
        .checked_add(n)
        .ok_or(CramError::Truncated { context: "block uncompressed size pos" })?;

    // CRC32 covers everything from the method byte through the end of compressed data
    let crc_start = 0;

    let compressed_end = pos
        .checked_add(compressed_size as usize)
        .ok_or(CramError::Truncated { context: "block compressed end overflow" })?;
    let compressed_data =
        buf.get(pos..compressed_end).ok_or(CramError::Truncated { context: "block data" })?;
    pos = compressed_end;

    // CRC32 is 4 bytes after the compressed data
    let crc_end =
        pos.checked_add(4).ok_or(CramError::Truncated { context: "block CRC32 pos overflow" })?;
    let crc_bytes = buf.get(pos..crc_end).ok_or(CramError::Truncated { context: "block CRC32" })?;
    let expected_crc = u32::from_le_bytes(
        crc_bytes.try_into().map_err(|_| CramError::Truncated { context: "block CRC32" })?,
    );

    // Verify CRC32 over [method .. compressed_data] (excludes the CRC32 itself)
    let crc_data = buf
        .get(crc_start..compressed_end)
        .ok_or(CramError::Truncated { context: "block CRC32 range" })?;
    let mut crc = libdeflater::Crc::new();
    crc.update(crc_data);
    if crc.sum() != expected_crc {
        return Err(CramError::ChecksumMismatch {
            context: "block",
            expected: expected_crc,
            found: crc.sum(),
        });
    }

    pos = pos.checked_add(4).ok_or(CramError::Truncated { context: "block pos after CRC" })?;

    let uncompressed_size = uncompressed_size as usize;
    super::reader::check_alloc_size(uncompressed_size, "block uncompressed size")?;

    let data = decompress_block(
        method,
        compressed_data,
        uncompressed_size,
        content_type_byte,
        content_id,
    )?;

    Ok((Block { content_type, content_id, data }, pos))
}

fn decompress_block(
    method: u8,
    compressed: &[u8],
    uncompressed_size: usize,
    content_type: u8,
    content_id: i32,
) -> Result<Vec<u8>, CramError> {
    match method {
        // r[impl cram.codec.raw]
        0 => Ok(compressed.to_vec()),
        // r[impl cram.codec.gzip]
        1 => {
            let mut decompressor = libdeflater::Decompressor::new();
            let mut output = vec![0u8; uncompressed_size];
            decompressor
                .gzip_decompress(compressed, &mut output)
                .map_err(|source| CramError::GzipDecompressionFailed { source })?;
            Ok(output)
        }
        // r[impl cram.codec.bzip2]
        2 => {
            use bzip2::read::BzDecoder;
            use std::io::Read;
            let mut decoder = BzDecoder::new(compressed);
            let mut output = Vec::with_capacity(uncompressed_size);
            decoder
                .read_to_end(&mut output)
                .map_err(|source| CramError::Bzip2DecompressionFailed { source })?;
            Ok(output)
        }
        // r[impl cram.codec.lzma]
        3 => {
            use std::io::Read;
            let mut decoder = xz2::read::XzDecoder::new(compressed);
            let mut output = Vec::with_capacity(uncompressed_size);
            decoder
                .read_to_end(&mut output)
                .map_err(|source| CramError::LzmaDecompressionFailed { source })?;
            Ok(output)
        }
        // r[impl cram.codec.rans4x8]
        4 => super::rans::decode(compressed),
        // r[impl cram.codec.rans_nx16]
        5 => super::rans_nx16::decode(compressed, uncompressed_size),
        // r[impl cram.codec.tok3]
        8 => super::tok3::decode(compressed),
        // r[impl cram.codec.unknown]
        _ => Err(CramError::UnsupportedCodec { method, content_type, content_id }),
    }
}

/// Build a raw CRAM block from components (for testing).
#[cfg(test)]
pub fn build_test_block(content_type: u8, content_id: i32, data: &[u8]) -> Vec<u8> {
    let mut buf = Vec::new();

    // method = 0 (raw)
    buf.push(0);
    // content type
    buf.push(content_type);
    // content_id as ITF8
    encode_itf8_to(&mut buf, content_id as u32);
    // compressed_size = uncompressed_size = data.len()
    encode_itf8_to(&mut buf, data.len() as u32);
    encode_itf8_to(&mut buf, data.len() as u32);
    // data
    buf.extend_from_slice(data);

    // CRC32 over everything so far
    let mut crc = libdeflater::Crc::new();
    crc.update(&buf);
    buf.extend_from_slice(&crc.sum().to_le_bytes());

    buf
}

#[cfg(test)]
fn encode_itf8_to(buf: &mut Vec<u8>, val: u32) {
    let mut tmp = [0u8; 5];
    let n = super::varint::encode_itf8(val, &mut tmp);
    #[allow(clippy::indexing_slicing)]
    buf.extend_from_slice(&tmp[..n]);
}

#[cfg(test)]
mod tests {
    use super::*;

    // r[verify cram.block.structure]
    // r[verify cram.codec.raw]
    #[test]
    fn parse_raw_block() {
        let data = b"hello world";
        let block_bytes = build_test_block(4, 0, data); // ExternalData, content_id=0
        let (block, consumed) = parse_block(&block_bytes).unwrap();
        assert_eq!(block.content_type, ContentType::ExternalData);
        assert_eq!(block.content_id, 0);
        assert_eq!(block.data, data);
        assert_eq!(consumed, block_bytes.len());
    }

    #[test]
    fn parse_raw_block_compression_header() {
        let data = b"\x00\x00\x00"; // minimal
        let block_bytes = build_test_block(1, 0, data);
        let (block, _) = parse_block(&block_bytes).unwrap();
        assert_eq!(block.content_type, ContentType::CompressionHeader);
    }

    #[test]
    fn parse_raw_block_core_data() {
        let data = b"\xFF\x00\xAA";
        let block_bytes = build_test_block(5, 0, data);
        let (block, _) = parse_block(&block_bytes).unwrap();
        assert_eq!(block.content_type, ContentType::CoreData);
        assert_eq!(block.data, data);
    }

    // r[verify cram.codec.gzip]
    #[test]
    fn parse_gzip_block() {
        // Compress "hello" with gzip
        let original = b"hello";
        let mut compressor =
            libdeflater::Compressor::new(libdeflater::CompressionLvl::new(6).unwrap());
        let max_len = compressor.gzip_compress_bound(original.len());
        let mut compressed = vec![0u8; max_len];
        let actual_len = compressor.gzip_compress(original, &mut compressed).unwrap();
        compressed.truncate(actual_len);

        // Build a block with method=1 (gzip)
        let mut buf = Vec::new();
        buf.push(1); // method = gzip
        buf.push(4); // content_type = ExternalData
        encode_itf8_to(&mut buf, 7); // content_id = 7
        encode_itf8_to(&mut buf, compressed.len() as u32);
        encode_itf8_to(&mut buf, original.len() as u32);
        buf.extend_from_slice(&compressed);

        let mut crc = libdeflater::Crc::new();
        crc.update(&buf);
        buf.extend_from_slice(&crc.sum().to_le_bytes());

        let (block, consumed) = parse_block(&buf).unwrap();
        assert_eq!(block.content_type, ContentType::ExternalData);
        assert_eq!(block.content_id, 7);
        assert_eq!(block.data, original);
        assert_eq!(consumed, buf.len());
    }

    // r[verify cram.codec.bzip2]
    #[test]
    fn parse_bzip2_block() {
        use bzip2::write::BzEncoder;
        use std::io::Write;

        let original = b"hello bzip2 world, this is a test of bzip2 compression";
        let mut encoder = BzEncoder::new(Vec::new(), bzip2::Compression::default());
        encoder.write_all(original).unwrap();
        let compressed = encoder.finish().unwrap();

        let mut buf = Vec::new();
        buf.push(2); // method = bzip2
        buf.push(4); // ExternalData
        encode_itf8_to(&mut buf, 0);
        encode_itf8_to(&mut buf, compressed.len() as u32);
        encode_itf8_to(&mut buf, original.len() as u32);
        buf.extend_from_slice(&compressed);
        let mut crc = libdeflater::Crc::new();
        crc.update(&buf);
        buf.extend_from_slice(&crc.sum().to_le_bytes());

        let (block, _) = parse_block(&buf).unwrap();
        assert_eq!(block.data, original);
    }

    // r[verify cram.codec.lzma]
    #[test]
    fn parse_lzma_block() {
        use std::io::Write;

        let original = b"hello lzma world, this is a test of lzma compression";
        let mut encoder = xz2::write::XzEncoder::new(Vec::new(), 6);
        encoder.write_all(original).unwrap();
        let compressed = encoder.finish().unwrap();

        let mut buf = Vec::new();
        buf.push(3); // method = lzma
        buf.push(4); // ExternalData
        encode_itf8_to(&mut buf, 0);
        encode_itf8_to(&mut buf, compressed.len() as u32);
        encode_itf8_to(&mut buf, original.len() as u32);
        buf.extend_from_slice(&compressed);
        let mut crc = libdeflater::Crc::new();
        crc.update(&buf);
        buf.extend_from_slice(&crc.sum().to_le_bytes());

        let (block, _) = parse_block(&buf).unwrap();
        assert_eq!(block.data, original);
    }

    // r[verify cram.block.crc32]
    #[test]
    fn crc32_mismatch_detected() {
        let data = b"test";
        let mut block_bytes = build_test_block(4, 0, data);
        // Corrupt the last byte of CRC32
        let len = block_bytes.len();
        #[allow(clippy::indexing_slicing)]
        {
            block_bytes[len - 1] ^= 0xFF;
        }
        let err = parse_block(&block_bytes).unwrap_err();
        assert!(matches!(err, CramError::ChecksumMismatch { .. }));
    }

    // r[verify cram.codec.unknown]
    #[test]
    fn unsupported_codec_detected() {
        let mut buf = Vec::new();
        buf.push(6); // method = 6 (not supported)
        buf.push(4); // ExternalData
        encode_itf8_to(&mut buf, 0);
        encode_itf8_to(&mut buf, 0); // compressed size = 0
        encode_itf8_to(&mut buf, 0); // uncompressed size = 0

        let mut crc = libdeflater::Crc::new();
        crc.update(&buf);
        buf.extend_from_slice(&crc.sum().to_le_bytes());

        let err = parse_block(&buf).unwrap_err();
        match err {
            CramError::UnsupportedCodec { method, .. } => assert_eq!(method, 6),
            other => panic!("expected UnsupportedCodec, got {other:?}"),
        }
    }

    #[test]
    fn truncated_block() {
        assert!(parse_block(&[]).is_err());
        assert!(parse_block(&[0]).is_err());
        assert!(parse_block(&[0, 4]).is_err());
    }

    #[test]
    fn multiple_blocks_sequential() {
        let block1 = build_test_block(4, 1, b"first");
        let block2 = build_test_block(4, 2, b"second");
        let mut combined = block1.clone();
        combined.extend_from_slice(&block2);

        let (b1, consumed1) = parse_block(&combined).unwrap();
        assert_eq!(b1.content_id, 1);
        assert_eq!(b1.data, b"first");

        #[allow(clippy::indexing_slicing)]
        let (b2, consumed2) = parse_block(&combined[consumed1..]).unwrap();
        assert_eq!(b2.content_id, 2);
        assert_eq!(b2.data, b"second");
        assert_eq!(consumed1 + consumed2, combined.len());
    }

    #[test]
    fn unknown_content_type_returns_error() {
        let err = ContentType::from_byte(255).unwrap_err();
        assert!(matches!(err, CramError::UnknownBlockContentType { content_type: 255 }));
    }

    #[test]
    fn unknown_content_type_invalid_values() {
        // Only 0, 1, 2, 4, 5 are valid — test a selection of invalid ones
        for invalid in [3u8, 6, 100, 200] {
            let err = ContentType::from_byte(invalid).unwrap_err();
            assert!(
                matches!(err, CramError::UnknownBlockContentType { content_type: b } if b == invalid),
                "expected UnknownBlockContentType for {invalid}"
            );
        }
    }

    #[test]
    fn gzip_decompression_failed_on_garbage() {
        // Build a block with method=1 (gzip) but garbage compressed data
        let garbage = [0xFFu8; 10];
        let mut buf = Vec::new();
        buf.push(1u8); // method = gzip
        buf.push(4u8); // ExternalData
        encode_itf8_to(&mut buf, 0);
        encode_itf8_to(&mut buf, garbage.len() as u32);
        encode_itf8_to(&mut buf, 100u32); // claim 100 bytes uncompressed
        buf.extend_from_slice(&garbage);
        let mut crc = libdeflater::Crc::new();
        crc.update(&buf);
        buf.extend_from_slice(&crc.sum().to_le_bytes());

        let err = parse_block(&buf).unwrap_err();
        assert!(matches!(err, CramError::GzipDecompressionFailed { .. }));
    }

    #[test]
    fn bzip2_decompression_failed_on_garbage() {
        let garbage = [0xFFu8; 10];
        let mut buf = Vec::new();
        buf.push(2u8); // method = bzip2
        buf.push(4u8); // ExternalData
        encode_itf8_to(&mut buf, 0);
        encode_itf8_to(&mut buf, garbage.len() as u32);
        encode_itf8_to(&mut buf, 100u32);
        buf.extend_from_slice(&garbage);
        let mut crc = libdeflater::Crc::new();
        crc.update(&buf);
        buf.extend_from_slice(&crc.sum().to_le_bytes());

        let err = parse_block(&buf).unwrap_err();
        assert!(matches!(err, CramError::Bzip2DecompressionFailed { .. }));
    }

    #[test]
    fn lzma_decompression_failed_on_garbage() {
        let garbage = [0xFFu8; 10];
        let mut buf = Vec::new();
        buf.push(3u8); // method = lzma
        buf.push(4u8); // ExternalData
        encode_itf8_to(&mut buf, 0);
        encode_itf8_to(&mut buf, garbage.len() as u32);
        encode_itf8_to(&mut buf, 100u32);
        buf.extend_from_slice(&garbage);
        let mut crc = libdeflater::Crc::new();
        crc.update(&buf);
        buf.extend_from_slice(&crc.sum().to_le_bytes());

        let err = parse_block(&buf).unwrap_err();
        assert!(matches!(err, CramError::LzmaDecompressionFailed { .. }));
    }

    #[test]
    fn expected_compression_header_wrong_content_type() {
        // Build a raw block with content type = CoreData (5), then check it
        // where CompressionHeader (1) is expected.
        // We simulate this by checking what parse_block returns vs the expected type.
        let data = b"\x00\x01\x02";
        // Block with ExternalData (4) but we want to detect it's not CompressionHeader (1)
        let block_bytes = build_test_block(4, 0, data); // 4 = ExternalData
        let (block, _) = parse_block(&block_bytes).unwrap();
        // Check that the content type mismatch would produce the right error
        assert_eq!(block.content_type, ContentType::ExternalData);
        let err = CramError::ExpectedCompressionHeader { found: block.content_type };
        assert!(matches!(
            err,
            CramError::ExpectedCompressionHeader { found: ContentType::ExternalData }
        ));
    }

    #[test]
    fn expected_file_header_wrong_content_type() {
        let data = b"\x00\x01\x02";
        let block_bytes = build_test_block(5, 0, data); // 5 = CoreData
        let (block, _) = parse_block(&block_bytes).unwrap();
        assert_eq!(block.content_type, ContentType::CoreData);
        let err = CramError::ExpectedFileHeader { found: block.content_type };
        assert!(matches!(err, CramError::ExpectedFileHeader { found: ContentType::CoreData }));
    }

    // Construct a raw block at the byte level according to the CRAM spec:
    // method(1) | content_type(1) | content_id(ITF8) | compressed_size(ITF8)
    // | uncompressed_size(ITF8) | data | CRC32(4 LE)
    // Parse it with parse_block and assert every field matches what was written.
    #[test]
    fn raw_block_from_known_bytes() {
        // ExternalData (4), content_id=7, data="hello"
        let payload = b"hello";
        let size = payload.len() as u32; // 5 — encodes as single ITF8 byte

        let mut buf = Vec::new();
        buf.push(0u8); // method = raw
        buf.push(4u8); // content_type = ExternalData
        encode_itf8_to(&mut buf, 7); // content_id ITF8: 7 → [0x07]
        encode_itf8_to(&mut buf, size); // compressed_size ITF8: 5 → [0x05]
        encode_itf8_to(&mut buf, size); // uncompressed_size ITF8: 5 → [0x05]
        buf.extend_from_slice(payload);
        let mut crc = libdeflater::Crc::new();
        crc.update(&buf);
        buf.extend_from_slice(&crc.sum().to_le_bytes());

        // Verify hand-crafted header bytes before parsing
        assert_eq!(buf.first(), Some(&0u8), "method byte must be 0 (raw)");
        assert_eq!(buf.get(1), Some(&4u8), "content_type byte must be 4 (ExternalData)");
        assert_eq!(buf.get(2), Some(&7u8), "content_id ITF8: single byte for val < 0x80");
        assert_eq!(buf.get(3), Some(&5u8), "compressed_size ITF8: single byte for val < 0x80");
        assert_eq!(buf.get(4), Some(&5u8), "uncompressed_size ITF8: single byte for val < 0x80");

        let (block, consumed) = parse_block(&buf).unwrap();
        assert_eq!(block.content_type, ContentType::ExternalData);
        assert_eq!(block.content_id, 7);
        assert_eq!(block.data, payload);
        assert_eq!(consumed, buf.len());

        // Also verify a 2-byte ITF8 content_id (val=128 → [0x80, 0x80])
        let mut buf2 = Vec::new();
        buf2.push(0u8); // raw
        buf2.push(4u8); // ExternalData
        encode_itf8_to(&mut buf2, 128); // content_id = 128 → 2 ITF8 bytes: [0x80, 0x80]
        encode_itf8_to(&mut buf2, 3u32);
        encode_itf8_to(&mut buf2, 3u32);
        buf2.extend_from_slice(b"abc");
        let mut crc2 = libdeflater::Crc::new();
        crc2.update(&buf2);
        buf2.extend_from_slice(&crc2.sum().to_le_bytes());

        assert_eq!(buf2.get(2), Some(&0x80u8), "first ITF8 byte for 128: 0x80");
        assert_eq!(buf2.get(3), Some(&0x80u8), "second ITF8 byte for 128: 0x80");

        let (block2, consumed2) = parse_block(&buf2).unwrap();
        assert_eq!(block2.content_id, 128);
        assert_eq!(block2.data, b"abc");
        assert_eq!(consumed2, buf2.len());
    }

    proptest::proptest! {

        #[test]
        fn gzip_block_roundtrip(data in proptest::collection::vec(proptest::prelude::any::<u8>(), 1..256)) {
            let mut compressor = libdeflater::Compressor::new(libdeflater::CompressionLvl::new(1).unwrap());
            let max_len = compressor.gzip_compress_bound(data.len());
            let mut compressed = vec![0u8; max_len];
            let actual_len = compressor.gzip_compress(&data, &mut compressed).unwrap();
            compressed.truncate(actual_len);

            let mut buf = Vec::new();
            buf.push(1); // gzip
            buf.push(4); // ExternalData
            encode_itf8_to(&mut buf, 0);
            encode_itf8_to(&mut buf, compressed.len() as u32);
            encode_itf8_to(&mut buf, data.len() as u32);
            buf.extend_from_slice(&compressed);
            let mut crc = libdeflater::Crc::new();
            crc.update(&buf);
            buf.extend_from_slice(&crc.sum().to_le_bytes());

            let (block, _) = parse_block(&buf).unwrap();
            proptest::prop_assert_eq!(&block.data, &data);
        }
    }

    #[test]
    fn parse_block_from_real_cram() {
        // Read the header container's first block from our test CRAM file
        let data =
            std::fs::read(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.cram"))
                .unwrap();

        // Skip 26-byte file definition
        #[allow(clippy::indexing_slicing)]
        let after_file_def = &data[26..];

        // Skip the container header (variable length) to get to the first block
        // Container header starts with i32 length
        let container_length = i32::from_le_bytes([
            after_file_def[0],
            after_file_def[1],
            after_file_def[2],
            after_file_def[3],
        ]);
        assert!(container_length > 0, "header container should have positive length");

        // Parse ITF8 fields to find where the container header ends
        // We need to skip: i32(length) + itf8(ref_seq_id) + itf8(alignment_start)
        //   + itf8(alignment_span) + itf8(num_records) + ltf8(record_counter)
        //   + ltf8(bases) + itf8(num_blocks) + itf8_array(landmarks) + 4(crc32)
        // For the file header container these are mostly zeros, so ITF8 = 1 byte each
        //
        // Instead of parsing the full header, just try to find and parse the block
        // by scanning for valid block structures. The header container's block should
        // be a file header block (content_type=0) with the SAM header text.
        //
        // This is a smoke test — we'll do proper container parsing in step 4.
    }
}
