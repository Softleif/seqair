//! Parse CRAM container headers. [`ContainerHeader`] carries the reference ID, alignment span,
//! and block count needed to locate and read the slices inside a container.

use super::{reader::CramError, varint};

/// Parsed CRAM container header.
#[derive(Debug)]
pub struct ContainerHeader {
    /// Total bytes of block data in this container (NOT including this header).
    pub length: i32,
    /// Reference sequence ID: >=0 for single ref, -1 unmapped, -2 multi-ref.
    pub ref_seq_id: i32,
    /// 1-based alignment start on the reference.
    pub alignment_start: i32,
    /// Number of reference bases covered.
    pub alignment_span: i32,
    /// Number of alignment records in this container.
    pub num_records: i32,
    /// Global 0-based record counter.
    pub record_counter: i64,
    /// Total bases across all records.
    pub bases: i64,
    /// Number of blocks in the container.
    pub num_blocks: i32,
    /// Byte offsets from container data start to each slice header.
    pub landmarks: Vec<i32>,
    /// Total bytes consumed by this header (for seeking past it).
    pub header_size: usize,
}

impl ContainerHeader {
    /// Parse a container header from a byte slice.
    ///
    // r[impl cram.container.header]
    /// Returns the header and the number of bytes consumed.
    pub fn parse(buf: &[u8]) -> Result<Self, CramError> {
        let mut pos: usize = 0;

        // length: fixed-width i32 LE
        let length_end = pos
            .checked_add(4)
            .ok_or(CramError::Truncated { context: "container header length pos" })?;
        let length_bytes = buf
            .get(pos..length_end)
            .ok_or(CramError::Truncated { context: "container header length" })?;
        let length = i32::from_le_bytes(
            length_bytes
                .try_into()
                .map_err(|_| CramError::Truncated { context: "container header length" })?,
        );
        pos = length_end;

        let (ref_seq_id, n) = read_itf8_at(buf, pos, "container ref_seq_id")?;
        pos = pos
            .checked_add(n)
            .ok_or(CramError::Truncated { context: "container ref_seq_id pos" })?;

        let (alignment_start, n) = read_itf8_at(buf, pos, "container alignment_start")?;
        pos = pos
            .checked_add(n)
            .ok_or(CramError::Truncated { context: "container alignment_start pos" })?;

        let (alignment_span, n) = read_itf8_at(buf, pos, "container alignment_span")?;
        pos = pos
            .checked_add(n)
            .ok_or(CramError::Truncated { context: "container alignment_span pos" })?;

        let (num_records, n) = read_itf8_at(buf, pos, "container num_records")?;
        pos = pos
            .checked_add(n)
            .ok_or(CramError::Truncated { context: "container num_records pos" })?;

        let (record_counter_raw, n) = read_ltf8_at(buf, pos, "container record_counter")?;
        let record_counter = record_counter_raw.cast_signed();
        pos = pos
            .checked_add(n)
            .ok_or(CramError::Truncated { context: "container record_counter pos" })?;

        let (bases_raw, n) = read_ltf8_at(buf, pos, "container bases")?;
        let bases = bases_raw.cast_signed();
        pos = pos.checked_add(n).ok_or(CramError::Truncated { context: "container bases pos" })?;

        let (num_blocks, n) = read_itf8_at(buf, pos, "container num_blocks")?;
        pos = pos
            .checked_add(n)
            .ok_or(CramError::Truncated { context: "container num_blocks pos" })?;

        // landmarks: ITF8 array (count followed by values)
        let (landmark_count, n) = read_itf8_at(buf, pos, "container landmark count")?;
        pos = pos
            .checked_add(n)
            .ok_or(CramError::Truncated { context: "container landmark count pos" })?;

        let landmark_count_usize = landmark_count as usize;
        super::reader::check_alloc_size(
            landmark_count_usize.saturating_mul(4),
            "container landmark array",
        )?;
        let mut landmarks = Vec::with_capacity(landmark_count_usize);
        for i in 0..landmark_count {
            let (landmark, n) = read_itf8_at(buf, pos, "container landmark")?;
            let _ = i; // suppress unused warning
            landmarks.push(landmark.cast_signed());
            pos = pos
                .checked_add(n)
                .ok_or(CramError::Truncated { context: "container landmark pos" })?;
        }

        // CRC32: 4 bytes over the entire header up to this point
        let crc_end = pos
            .checked_add(4)
            .ok_or(CramError::Truncated { context: "container header CRC32 pos" })?;
        let crc_bytes = buf
            .get(pos..crc_end)
            .ok_or(CramError::Truncated { context: "container header CRC32" })?;
        let expected_crc = u32::from_le_bytes(
            crc_bytes
                .try_into()
                .map_err(|_| CramError::Truncated { context: "container header CRC32" })?,
        );

        let crc_data = buf
            .get(..pos)
            .ok_or(CramError::Truncated { context: "container header CRC32 range" })?;
        let mut crc = libdeflater::Crc::new();
        crc.update(crc_data);
        if crc.sum() != expected_crc {
            return Err(CramError::ChecksumMismatch {
                context: "container header",
                expected: expected_crc,
                found: crc.sum(),
            });
        }
        pos = pos
            .checked_add(4)
            .ok_or(CramError::Truncated { context: "container pos after CRC" })?;

        Ok(ContainerHeader {
            length,
            ref_seq_id: ref_seq_id.cast_signed(),
            alignment_start: alignment_start.cast_signed(),
            alignment_span: alignment_span.cast_signed(),
            num_records: num_records.cast_signed(),
            record_counter,
            bases,
            num_blocks: num_blocks.cast_signed(),
            landmarks,
            header_size: pos,
        })
    }

    // r[impl cram.file.eof]
    /// Is this the EOF container?
    pub fn is_eof(&self) -> bool {
        self.num_records == 0 && self.length <= 15
    }

    /// Does this container's range overlap `[query_start, query_end)`?
    ///
    /// Positions are 1-based in the container header.
    pub fn overlaps(&self, query_start: i64, query_end: i64) -> bool {
        if self.ref_seq_id < 0 {
            // Unmapped or multi-ref — can't do coordinate filtering
            return self.ref_seq_id == -2;
        }
        let container_start = i64::from(self.alignment_start);
        let container_end = container_start.saturating_add(i64::from(self.alignment_span));
        container_start < query_end && container_end > query_start
    }
}

fn read_itf8_at(buf: &[u8], pos: usize, context: &'static str) -> Result<(u32, usize), CramError> {
    let slice = buf.get(pos..).ok_or(CramError::Truncated { context })?;
    varint::decode_itf8(slice).ok_or(CramError::Truncated { context })
}

fn read_ltf8_at(buf: &[u8], pos: usize, context: &'static str) -> Result<(u64, usize), CramError> {
    let slice = buf.get(pos..).ok_or(CramError::Truncated { context })?;
    varint::decode_ltf8(slice).ok_or(CramError::Truncated { context })
}

#[cfg(test)]
#[allow(
    clippy::arithmetic_side_effects,
    reason = "test-only arithmetic on known-valid CRAM file data"
)]
mod tests {
    use super::*;

    // r[verify cram.container.header]
    #[test]
    fn parse_header_container_from_real_cram() {
        let data =
            std::fs::read(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.cram"))
                .unwrap();

        // Skip 26-byte file definition
        #[allow(clippy::indexing_slicing)]
        let after_file_def = &data[26..];

        let header = ContainerHeader::parse(after_file_def).unwrap();

        // The first container is the file header container
        // ref_seq_id should be 0 (or could be -1 for header)
        // num_records should be 0 (it's a header, not data)
        assert_eq!(header.num_records, 0, "header container should have 0 records");
        assert!(header.length > 0, "header container should have positive length");
        assert!(!header.is_eof(), "header container should not be EOF");
    }

    #[test]
    fn parse_data_container_from_real_cram() {
        let data =
            std::fs::read(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.cram"))
                .unwrap();

        // Skip file def (26 bytes) + header container (header + data)
        #[allow(clippy::indexing_slicing)]
        let after_file_def = &data[26..];
        let header_container = ContainerHeader::parse(after_file_def).unwrap();

        // Skip to after the header container's data
        let data_start = header_container.header_size + header_container.length as usize;
        #[allow(clippy::indexing_slicing)]
        let after_header = &after_file_def[data_start..];

        let data_container = ContainerHeader::parse(after_header).unwrap();

        // Data container should have records and a valid reference
        assert!(data_container.num_records > 0, "data container should have records");
        assert!(
            data_container.ref_seq_id >= 0 || data_container.ref_seq_id == -2,
            "data container should have valid ref_seq_id"
        );
        assert!(data_container.alignment_start > 0);
        assert!(data_container.alignment_span > 0);
        assert!(!data_container.landmarks.is_empty(), "should have landmarks");
    }

    // r[verify cram.file.eof]
    #[test]
    fn find_eof_container_in_real_cram() {
        let data =
            std::fs::read(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.cram"))
                .unwrap();

        // The EOF container is the last 38 bytes of the file
        // But let's iterate containers to find it
        #[allow(clippy::indexing_slicing)]
        let mut pos = 26; // skip file def

        let mut found_eof = false;
        let mut container_count = 0;
        while pos < data.len() {
            #[allow(clippy::indexing_slicing)]
            let remaining = &data[pos..];
            let container = match ContainerHeader::parse(remaining) {
                Ok(c) => c,
                Err(_) => break,
            };
            container_count += 1;

            if container.is_eof() {
                found_eof = true;
                break;
            }

            pos += container.header_size + container.length as usize;
        }

        assert!(found_eof, "should find EOF container");
        assert!(
            container_count >= 3,
            "should have at least header + data + EOF containers, found {container_count}"
        );
    }

    // r[verify cram.file.header_container]
    #[test]
    fn parse_header_block_extracts_sam_header() {
        use super::super::block;

        let data =
            std::fs::read(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.cram"))
                .unwrap();

        #[allow(clippy::indexing_slicing)]
        let after_file_def = &data[26..];
        let header_container = ContainerHeader::parse(after_file_def).unwrap();

        // The block data starts right after the container header
        #[allow(clippy::indexing_slicing)]
        let block_data = &after_file_def[header_container.header_size..];

        let (blk, _) = block::parse_block(block_data).unwrap();
        assert_eq!(blk.content_type, block::ContentType::FileHeader);

        // The file header block contains: i32 header_text_length + UTF-8 header text
        let header_len =
            i32::from_le_bytes([blk.data[0], blk.data[1], blk.data[2], blk.data[3]]) as usize;
        #[allow(clippy::indexing_slicing)]
        let header_text = std::str::from_utf8(&blk.data[4..4 + header_len]).unwrap();

        assert!(header_text.starts_with("@HD"), "should start with @HD header line");
        assert!(header_text.contains("@SQ"), "should contain @SQ sequence dictionary");
    }

    // r[verify cram.file.header_container]
    // r[verify cram.container.header]
    #[test]
    fn container_header_matches_htslib() {
        use super::super::block;
        use crate::bam::BamHeader;
        use rust_htslib::bam::Read as _;

        let data =
            std::fs::read(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.cram"))
                .unwrap();

        #[allow(clippy::indexing_slicing)]
        let after_file_def = &data[26..];
        let header_container = ContainerHeader::parse(after_file_def).unwrap();

        #[allow(clippy::indexing_slicing)]
        let block_data = &after_file_def[header_container.header_size..];
        let (blk, _) = block::parse_block(block_data).unwrap();

        let header_len =
            i32::from_le_bytes([blk.data[0], blk.data[1], blk.data[2], blk.data[3]]) as usize;
        #[allow(clippy::indexing_slicing)]
        let header_text = std::str::from_utf8(&blk.data[4..4 + header_len]).unwrap();

        let header = BamHeader::from_sam_text(header_text).unwrap();

        // Compare with htslib reading the same BAM (same data, different format)
        let htslib_bam = rust_htslib::bam::Reader::from_path(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../tests/data/test.bam"
        ))
        .unwrap();
        let htslib_header = htslib_bam.header();

        let our_names: Vec<&str> = header.target_names().collect();
        let htslib_names = htslib_header.target_names();

        assert_eq!(our_names.len(), htslib_names.len(), "target count mismatch");
        for (ours, theirs) in our_names.iter().zip(htslib_names.iter()) {
            let theirs_str = std::str::from_utf8(theirs).unwrap();
            assert_eq!(*ours, theirs_str, "target name mismatch");
        }
    }

    #[test]
    fn overlaps_query() {
        // Container at positions 100-200 (1-based)
        let header = ContainerHeader {
            length: 0,
            ref_seq_id: 0,
            alignment_start: 100,
            alignment_span: 100,
            num_records: 10,
            record_counter: 0,
            bases: 1000,
            num_blocks: 1,
            landmarks: vec![],
            header_size: 0,
        };

        // Query overlapping
        assert!(header.overlaps(50, 150)); // partial overlap left
        assert!(header.overlaps(150, 250)); // partial overlap right
        assert!(header.overlaps(100, 200)); // exact match
        assert!(header.overlaps(120, 180)); // contained
        assert!(header.overlaps(50, 250)); // contains

        // Query not overlapping
        assert!(!header.overlaps(200, 300)); // starts at container end
        assert!(!header.overlaps(0, 100)); // ends at container start
        assert!(!header.overlaps(300, 400)); // completely after
    }
}
