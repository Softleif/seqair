#![no_main]

use arbitrary::Arbitrary;
use libfuzzer_sys::fuzz_target;
use seqair::cram::compression_header::CompressionHeader;
use seqair::cram::encoding::{DecodeContext, ExternalCursor};

/// Full-stack CRAM decode: parse compression header, then decode values
/// through the encoding descriptors using fuzzed core + external data.
#[derive(Arbitrary, Debug)]
struct CramDecodeInput {
    /// Bytes for CompressionHeader::parse
    header_data: Vec<u8>,
    /// Core data block (BitReader source for Huffman, Beta, Gamma, etc.)
    core_data: Vec<u8>,
    /// External block data keyed by content_id (0..4)
    external_0: Vec<u8>,
    external_1: Vec<u8>,
    external_2: Vec<u8>,
    external_3: Vec<u8>,
}

fuzz_target!(|input: CramDecodeInput| {
    let ch = match CompressionHeader::parse(&input.header_data) {
        Ok(ch) => ch,
        Err(_) => return,
    };

    // Build external blocks map with common content IDs
    let external = [
        (0, ExternalCursor::new(input.external_0)),
        (1, ExternalCursor::new(input.external_1)),
        (2, ExternalCursor::new(input.external_2)),
        (3, ExternalCursor::new(input.external_3)),
    ]
    .into();

    let mut ctx = DecodeContext::new(&input.core_data, external);

    // Try decoding through each data series encoding (up to 4 rounds)
    for _ in 0..4 {
        let _ = ch.data_series.bam_flags.decode(&mut ctx);
        let _ = ch.data_series.cram_flags.decode(&mut ctx);
        let _ = ch.data_series.read_length.decode(&mut ctx);
        let _ = ch.data_series.alignment_pos.decode(&mut ctx);
        let _ = ch.data_series.feature_count.decode(&mut ctx);
        let _ = ch.data_series.feature_code.decode(&mut ctx);
        let _ = ch.data_series.feature_pos.decode(&mut ctx);
        let _ = ch.data_series.deletion_length.decode(&mut ctx);
        let _ = ch.data_series.mapping_quality.decode(&mut ctx);
        let _ = ch.data_series.base.decode(&mut ctx);
        let _ = ch.data_series.quality_score.decode(&mut ctx);
        let _ = ch.data_series.base_sub.decode(&mut ctx);
        let _ = ch.data_series.tag_line.decode(&mut ctx);
        let _ = ch.data_series.ref_id.decode(&mut ctx);
        let _ = ch.data_series.read_group.decode(&mut ctx);
    }
});
