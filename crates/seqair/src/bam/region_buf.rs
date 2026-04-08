//! Bulk-read BGZF buffer for high-latency I/O (cluster/NFS storage).
//!
//! Reads all compressed bytes for a set of BAM index chunks into memory
//! in one large I/O operation, then decompresses from RAM.

use tracing::instrument;

use super::{
    bgzf::{self, BgzfError, VirtualOffset},
    index::Chunk,
};
use std::io::{Read, Seek, SeekFrom};

pub const PROFILE_TARGET: &str = "seqair::profile";

const BGZF_HEADER_SIZE: usize = 18;
const BGZF_FOOTER_SIZE: usize = 8;
const BGZF_MAGIC: [u8; 4] = [0x1f, 0x8b, 0x08, 0x04];
const MAX_BLOCK_SIZE: usize = 65536;

/// Pre-computed byte range covering one or more merged index chunks.
struct MergedRange {
    file_start: u64,
    file_end: u64,
}

/// Maps a file offset range to a position in the buffer.
#[derive(Debug, Clone)]
struct RangeMapping {
    /// Start file offset of this range.
    file_start: u64,
    /// End file offset (exclusive) of this range.
    file_end: u64,
    /// Offset in the data buffer where this range's bytes begin.
    buf_start: usize,
}

/// In-memory BGZF decompressor for a pre-fetched region.
///
/// Created by [`RegionBuf::load`], which bulk-reads all compressed bytes
/// for a region's index chunks into memory. All subsequent decompression
/// happens from RAM with zero network I/O.
pub struct RegionBuf {
    /// Raw compressed bytes covering all merged chunk ranges (concatenated).
    data: Vec<u8>,
    /// Maps file offset ranges to buffer positions for disjoint ranges.
    ranges: Vec<RangeMapping>,
    /// Current read position within `data`.
    cursor: usize,
    /// Decompressed block buffer (reused across blocks).
    buf: Vec<u8>,
    /// Current position within `buf`.
    buf_pos: usize,
    /// File offset of the current BGZF block.
    block_offset: u64,
    eof: bool,
    decompressor: libdeflater::Decompressor,
    blocks_decompressed: u32,
    decompressed_bytes: u64,
}

impl std::fmt::Debug for RegionBuf {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("RegionBuf")
            .field("data_len", &self.data.len())
            .field("ranges", &self.ranges.len())
            .field("cursor", &self.cursor)
            .field("block_offset", &self.block_offset)
            .field("eof", &self.eof)
            .finish()
    }
}

impl RegionBuf {
    // r[impl region_buf.load]
    // r[impl region_buf.empty]
    // r[related region_buf.merge_chunks]
    /// Bulk-read all compressed bytes for the given chunks into memory.
    ///
    /// Merges overlapping/adjacent chunks to minimize I/O, then performs
    /// one large read per merged range.
    #[instrument(level = "trace", skip_all, fields(input_chunks = chunks.len()))]
    pub fn load<R: Read + Seek>(reader: &mut R, chunks: &[Chunk]) -> Result<Self, BgzfError> {
        if chunks.is_empty() {
            return Ok(Self {
                data: Vec::new(),
                ranges: Vec::new(),
                cursor: 0,
                buf: Vec::new(),
                buf_pos: 0,
                block_offset: 0,
                eof: true,
                decompressor: libdeflater::Decompressor::new(),
                blocks_decompressed: 0,
                decompressed_bytes: 0,
            });
        }

        let merged = merge_chunks(chunks);
        let total_bytes: usize = merged
            .iter()
            .map(|r| r.file_end.saturating_sub(r.file_start) as usize)
            .fold(0usize, usize::saturating_add);

        // Reject obviously corrupt chunk ranges that would cause OOM
        const MAX_REGION_BYTES: usize = 256 * 1024 * 1024; // 256 MiB
        if total_bytes > MAX_REGION_BYTES {
            return Err(BgzfError::RecordTooLarge { block_size: total_bytes });
        }

        let load_start = std::time::Instant::now();
        let mut data = Vec::with_capacity(total_bytes);
        let mut range_map = Vec::with_capacity(merged.len());
        let mut max_range_us: u64 = 0;

        for range in &merged {
            let range_start = std::time::Instant::now();
            reader.seek(SeekFrom::Start(range.file_start)).map_err(|_| BgzfError::SeekFailed)?;

            let len = range.file_end.saturating_sub(range.file_start) as usize;
            if len > MAX_REGION_BYTES {
                return Err(BgzfError::RecordTooLarge { block_size: len });
            }
            let buf_start = data.len();
            data.resize(buf_start.saturating_add(len), 0);
            // Read as much as available — the last range may extend past EOF
            #[allow(
                clippy::indexing_slicing,
                reason = "buf_start = pre-resize len, within bounds after resize"
            )]
            let actually_read = read_all(reader, &mut data[buf_start..]);
            let actual_file_end = range.file_start.wrapping_add(actually_read as u64);
            data.truncate(buf_start.wrapping_add(actually_read));

            max_range_us = max_range_us.max(range_start.elapsed().as_micros() as u64);

            range_map.push(RangeMapping {
                file_start: range.file_start,
                file_end: actual_file_end,
                buf_start,
            });
        }

        let first_offset = merged.first().map_or(0, |r| r.file_start);
        let last_offset = merged.last().map_or(0, |r| r.file_start);
        let file_span = last_offset.saturating_sub(first_offset);

        tracing::debug!(
            target: PROFILE_TARGET,
            input_chunks = chunks.len(),
            merged_ranges = merged.len(),
            total_read_bytes = data.len(),
            file_span,
            first_offset,
            last_offset,
            max_range_us,
            elapsed_us = load_start.elapsed().as_micros() as u64,
            "region_buf::load",
        );

        let first_file_start = range_map.first().map(|r| r.file_start).unwrap_or(0);

        Ok(RegionBuf {
            data,
            ranges: range_map,
            cursor: 0,
            buf: Vec::with_capacity(MAX_BLOCK_SIZE),
            buf_pos: 0,
            block_offset: first_file_start,
            eof: false,
            decompressor: libdeflater::Decompressor::new(),
            blocks_decompressed: 0,
            decompressed_bytes: 0,
        })
    }

    // r[impl region_buf.seek_virtual]
    /// Seek to a virtual offset within the pre-loaded data.
    #[instrument(level = "trace", skip(self))]
    pub fn seek_virtual(&mut self, voff: VirtualOffset) -> Result<(), BgzfError> {
        let block_off = voff.block_offset();
        let within = voff.within_block() as usize;

        let cursor_pos = self.file_offset_to_cursor(block_off)?;

        self.cursor = cursor_pos;
        self.block_offset = block_off;
        self.buf.clear();
        self.buf_pos = 0;
        self.eof = false;

        if within > 0 {
            self.read_block()?;
            self.buf_pos = within.min(self.buf.len());
        }

        Ok(())
    }

    /// Translate a file offset to a buffer cursor position using the range map.
    fn file_offset_to_cursor(&self, file_off: u64) -> Result<usize, BgzfError> {
        for rm in &self.ranges {
            if file_off >= rm.file_start && file_off < rm.file_end {
                return Ok(rm
                    .buf_start
                    .wrapping_add(file_off.wrapping_sub(rm.file_start) as usize));
            }
        }
        Err(BgzfError::VirtualOffsetOutOfRange { offset: file_off })
    }

    /// Translate a buffer cursor position back to a file offset.
    fn cursor_to_file_offset(&self) -> u64 {
        for rm in &self.ranges {
            let buf_end =
                rm.buf_start.wrapping_add(rm.file_end.wrapping_sub(rm.file_start) as usize);
            if self.cursor >= rm.buf_start && self.cursor < buf_end {
                return rm.file_start.wrapping_add(self.cursor.wrapping_sub(rm.buf_start) as u64);
            }
        }
        // Past all ranges — return best estimate from the last range
        if let Some(last) = self.ranges.last() {
            let buf_end =
                last.buf_start.wrapping_add(last.file_end.wrapping_sub(last.file_start) as usize);
            last.file_end.wrapping_add(self.cursor.saturating_sub(buf_end) as u64)
        } else {
            self.cursor as u64
        }
    }

    // r[impl region_buf.virtual_offset]
    pub fn virtual_offset(&self) -> VirtualOffset {
        VirtualOffset::new(self.block_offset, self.buf_pos as u16)
    }

    // r[impl region_buf.decompress]
    // r[impl region_buf.fast_header]
    fn read_block(&mut self) -> Result<bool, BgzfError> {
        self.block_offset = self.cursor_to_file_offset();

        if self.cursor.wrapping_add(BGZF_HEADER_SIZE) > self.data.len() {
            self.eof = true;
            self.buf.clear();
            self.buf_pos = 0;
            return Ok(false);
        }

        // header is always BGZF_HEADER_SIZE=18 bytes (bounds checked above).
        debug_assert!(
            self.cursor.wrapping_add(BGZF_HEADER_SIZE) <= self.data.len(),
            "header overrun: {} > {}",
            self.cursor.wrapping_add(BGZF_HEADER_SIZE),
            self.data.len()
        );
        #[allow(
            clippy::indexing_slicing,
            reason = "header length = BGZF_HEADER_SIZE = 18; all indices < 18"
        )]
        let header = &self.data[self.cursor..self.cursor.wrapping_add(BGZF_HEADER_SIZE)];

        #[allow(clippy::indexing_slicing, reason = "header.len() = 18; indices < 18")]
        if header[..4] != BGZF_MAGIC {
            return Err(BgzfError::InvalidMagic);
        }

        #[allow(clippy::indexing_slicing, reason = "header.len() = 18; indices 10, 11 < 18")]
        let xlen = u16::from_le_bytes([header[10], header[11]]) as usize;

        // Fast path: standard BGZF header
        #[allow(clippy::indexing_slicing, reason = "header.len() = 18; indices 12-17 < 18")]
        let bsize = if xlen == 6
            && header[12] == b'B'
            && header[13] == b'C'
            && header[14] == 2
            && header[15] == 0
        {
            u16::from_le_bytes([header[16], header[17]])
        } else {
            let extra_start = self.cursor.wrapping_add(12);
            let extra_end = extra_start.wrapping_add(xlen);
            if extra_end > self.data.len() {
                return Err(BgzfError::TruncatedBlock);
            }
            debug_assert!(
                extra_end <= self.data.len(),
                "extra field overrun: {extra_end} > {}",
                self.data.len()
            );
            #[allow(clippy::indexing_slicing, reason = "extra_end ≤ data.len() checked above")]
            bgzf::find_bsize(&self.data[extra_start..extra_end]).ok_or(BgzfError::MissingBsize)?
        };

        let total_block_size = (bsize as usize).wrapping_add(1);
        let block_end = self.cursor.wrapping_add(total_block_size);

        if block_end > self.data.len() {
            // Partial block at end of loaded data — treat as EOF
            self.eof = true;
            self.buf.clear();
            self.buf_pos = 0;
            return Ok(false);
        }

        let data_start = self.cursor.saturating_add(12).saturating_add(xlen);
        if data_start > block_end {
            return Err(BgzfError::BlockSizeTooSmall { bsize });
        }
        let remaining = self.data.get(data_start..block_end).ok_or(BgzfError::TruncatedBlock)?;

        if remaining.len() < BGZF_FOOTER_SIZE {
            return Err(BgzfError::TruncatedBlock);
        }

        let footer_start = remaining.len().wrapping_sub(BGZF_FOOTER_SIZE);
        // footer_start + 8 = remaining.len() ≤ remaining.len() is trivially true.
        debug_assert!(
            footer_start.wrapping_add(8) <= remaining.len(),
            "footer overrun: {} > {}",
            footer_start.wrapping_add(8),
            remaining.len()
        );
        #[allow(clippy::indexing_slicing, reason = "footer_start + 8 = remaining.len()")]
        let crc32_bytes: [u8; 4] = remaining[footer_start..footer_start.wrapping_add(4)]
            .try_into()
            .map_err(|_| BgzfError::TruncatedBlock)?;
        let expected_crc = u32::from_le_bytes(crc32_bytes);

        #[allow(clippy::indexing_slicing, reason = "footer_start + 8 = remaining.len()")]
        let isize_bytes: [u8; 4] = remaining
            [footer_start.wrapping_add(4)..footer_start.wrapping_add(8)]
            .try_into()
            .map_err(|_| BgzfError::TruncatedBlock)?;
        let uncompressed_size = u32::from_le_bytes(isize_bytes) as usize;

        if uncompressed_size > MAX_BLOCK_SIZE {
            return Err(BgzfError::UncompressedSizeTooLarge { isize_value: uncompressed_size });
        }

        self.cursor = block_end;

        if uncompressed_size == 0 {
            self.eof = true;
            self.buf.clear();
            self.buf_pos = 0;
            return Ok(false);
        }

        debug_assert!(
            footer_start <= remaining.len(),
            "deflate slice overrun: {footer_start} > {}",
            remaining.len()
        );
        #[allow(clippy::indexing_slicing, reason = "footer_start ≤ remaining.len()")]
        let deflate_data = &remaining[..footer_start];

        // Safety: all bytes written by deflate_decompress before any read.
        unsafe { bgzf::resize_uninit(&mut self.buf, uncompressed_size) };
        let actual = self
            .decompressor
            .deflate_decompress(deflate_data, &mut self.buf)
            .map_err(|source| BgzfError::DecompressionFailed { source })?;
        self.buf.truncate(actual);

        // Verify CRC32
        let mut crc = libdeflater::Crc::new();
        crc.update(&self.buf);
        if crc.sum() != expected_crc {
            return Err(BgzfError::ChecksumMismatch { expected: expected_crc, found: crc.sum() });
        }

        self.buf_pos = 0;
        self.blocks_decompressed = self.blocks_decompressed.wrapping_add(1);
        self.decompressed_bytes = self.decompressed_bytes.wrapping_add(actual as u64);
        Ok(true)
    }

    // r[impl region_buf.read_exact]
    #[inline]
    pub fn read_exact_into(&mut self, out: &mut [u8]) -> Result<(), BgzfError> {
        let mut written = 0;
        while written < out.len() {
            if self.buf_pos >= self.buf.len() && !self.read_block()? {
                return Err(BgzfError::UnexpectedEof);
            }
            let avail = self.buf.len().wrapping_sub(self.buf_pos);
            let need = out.len().wrapping_sub(written);
            let n = avail.min(need);
            let dst =
                out.get_mut(written..written.wrapping_add(n)).ok_or(BgzfError::TruncatedBlock)?;
            let src = self
                .buf
                .get(self.buf_pos..self.buf_pos.wrapping_add(n))
                .ok_or(BgzfError::TruncatedBlock)?;
            dst.copy_from_slice(src);
            self.buf_pos = self.buf_pos.wrapping_add(n);
            written = written.wrapping_add(n);
        }
        Ok(())
    }

    #[inline]
    pub fn read_byte(&mut self) -> Result<u8, BgzfError> {
        if self.buf_pos >= self.buf.len() && !self.read_block()? {
            return Err(BgzfError::UnexpectedEof);
        }
        let b = self.buf.get(self.buf_pos).copied().ok_or(BgzfError::TruncatedBlock)?;
        self.buf_pos = self.buf_pos.wrapping_add(1);
        Ok(b)
    }

    pub fn read_u32(&mut self) -> Result<u32, BgzfError> {
        let mut buf = [0u8; 4];
        self.read_exact_into(&mut buf)?;
        Ok(u32::from_le_bytes(buf))
    }

    /// Read a complete BAM record body (block_size bytes, not including the
    /// 4-byte length prefix which this method reads itself) and return a slice
    /// of its bytes.
    ///
    /// Fast path: when the entire record body lies within the current
    /// decompressed BGZF block, returns a zero-copy `&[u8]` directly from the
    /// internal buffer, avoiding any allocation or copy.
    ///
    /// Slow path: when the record spans a block boundary (rare — BGZF blocks
    /// are 64 KB and typical BAM records are ≪64 KB), `scratch` is resized and
    /// filled, and a slice of `scratch` is returned instead.
    ///
    /// The returned slice is valid for the lifetime of whichever buffer it
    /// points into, expressed here as `'a` covering both `self` and `scratch`.
    pub fn read_record<'a>(&'a mut self, scratch: &'a mut Vec<u8>) -> Result<&'a [u8], BgzfError> {
        // Ensure the decompressed buffer has data to read the 4-byte length from.
        if self.buf_pos >= self.buf.len() && !self.read_block()? {
            return Err(BgzfError::UnexpectedEof);
        }

        // Fast-path u32 read: all 4 bytes in the current block.
        let block_size = if self.buf_pos.wrapping_add(4) <= self.buf.len() {
            let bytes = self
                .buf
                .get(self.buf_pos..self.buf_pos.wrapping_add(4))
                .ok_or(BgzfError::TruncatedBlock)?;
            let val = u32::from_le_bytes(bytes.try_into().map_err(|_| BgzfError::TruncatedBlock)?)
                as usize;
            self.buf_pos = self.buf_pos.wrapping_add(4);
            val
        } else {
            let mut len_buf = [0u8; 4];
            self.read_exact_into(&mut len_buf)?;
            u32::from_le_bytes(len_buf) as usize
        };

        // r[impl bam.record.max_size]
        const MAX_RECORD_SIZE: usize = 2 * 1024 * 1024; // 2 MiB
        if block_size > MAX_RECORD_SIZE {
            return Err(BgzfError::RecordTooLarge { block_size });
        }

        // Fast path: the entire record body is already in the decompressed buffer.
        if self.buf_pos.wrapping_add(block_size) <= self.buf.len() {
            let slice = self
                .buf
                .get(self.buf_pos..self.buf_pos.wrapping_add(block_size))
                .ok_or(BgzfError::TruncatedBlock)?;
            self.buf_pos = self.buf_pos.wrapping_add(block_size);
            return Ok(slice);
        }

        // Slow path: record spans a block boundary — copy into scratch.
        scratch.clear();
        // Safety: read_exact_into overwrites every byte before any are read.
        unsafe { super::bgzf::resize_uninit(scratch, block_size) };
        self.read_exact_into(scratch)?;
        Ok(scratch)
    }
}

// r[impl region_buf.drop_no_panic]
impl Drop for RegionBuf {
    fn drop(&mut self) {
        if self.blocks_decompressed > 0 {
            let max_gap = self
                .ranges
                .windows(2)
                .map(|w| {
                    w.get(1)
                        .map_or(0, |r| r.file_start)
                        .saturating_sub(w.first().map_or(0, |r| r.file_end))
                })
                .max()
                .unwrap_or(0);

            tracing::debug!(
                target: PROFILE_TARGET,
                blocks = self.blocks_decompressed,
                compressed_bytes = self.data.len(),
                decompressed_bytes = self.decompressed_bytes,
                ranges = self.ranges.len(),
                max_gap_bytes = max_gap,
                buf_capacity = self.buf.capacity(),
                "region_buf summary",
            );
        }
    }
}

fn read_all<R: Read>(reader: &mut R, buf: &mut [u8]) -> usize {
    let mut total = 0;
    while total < buf.len() {
        debug_assert!(
            total < buf.len(),
            "read_all index out of bounds: total={total}, len={}",
            buf.len()
        );
        #[allow(clippy::indexing_slicing, reason = "total < buf.len() by loop condition")]
        match reader.read(&mut buf[total..]) {
            Ok(0) => break,
            Ok(n) => total = total.wrapping_add(n),
            Err(_) => break,
        }
    }
    total
}

// r[impl region_buf.merge_chunks]
/// Merge chunks into non-overlapping byte ranges sorted by file offset.
fn merge_chunks(chunks: &[Chunk]) -> Vec<MergedRange> {
    let mut offsets: Vec<(u64, u64)> = chunks
        .iter()
        .map(|c| {
            let start = c.begin.block_offset();
            // end's block_offset points to the block containing the last byte;
            // extend by max block size to ensure we capture the full final block.
            let end = c.end.block_offset().saturating_add(MAX_BLOCK_SIZE as u64);
            (start, end)
        })
        .collect();

    offsets.sort_unstable();

    let mut merged: Vec<MergedRange> = Vec::with_capacity(offsets.len());
    for (start, end) in offsets {
        if let Some(last) = merged.last_mut()
            && start <= last.file_end
        {
            last.file_end = last.file_end.max(end);
            continue;
        }
        merged.push(MergedRange { file_start: start, file_end: end });
    }

    merged
}

#[cfg(test)]
#[allow(clippy::arithmetic_side_effects, reason = "test arithmetic is not safety-critical")]
mod tests {
    use super::*;

    // r[verify region_buf.drop_no_panic]
    #[test]
    fn drop_does_not_panic_with_overlapping_ranges() {
        // Simulate a RegionBuf with overlapping ranges where file_start < prev file_end
        // (which would cause subtraction underflow without saturating_sub).
        let buf = RegionBuf {
            data: vec![0; 100],
            ranges: vec![
                RangeMapping { file_start: 100, file_end: 300, buf_start: 0 },
                RangeMapping { file_start: 200, file_end: 400, buf_start: 50 },
            ],
            cursor: 0,
            buf: Vec::new(),
            buf_pos: 0,
            block_offset: 0,
            eof: false,
            blocks_decompressed: 1,
            decompressed_bytes: 100,
            decompressor: libdeflater::Decompressor::new(),
        };
        // If Drop panics, the test will fail.
        drop(buf);
    }

    // r[verify bam.record.max_size]
    #[test]
    fn read_record_rejects_huge_block_size() {
        // Build a BGZF file with a single block whose decompressed content
        // starts with a u32 block_size = 0xFFFF_FFFF (4 GB), which should be rejected.
        let mut payload = Vec::new();
        payload.extend_from_slice(&u32::MAX.to_le_bytes()); // block_size = 4GB
        payload.extend_from_slice(&[0u8; 28]); // some padding

        let block = make_bgzf_block(&payload);
        let mut file = Vec::new();
        file.extend_from_slice(&block);
        file.extend_from_slice(&make_bgzf_eof());

        let offsets = [0u64];
        let chunks = vec![Chunk {
            begin: VirtualOffset::new(offsets[0], 0),
            end: VirtualOffset::new(offsets[0] + 1, 0),
        }];

        let mut cursor = std::io::Cursor::new(file);
        let mut buf = RegionBuf::load(&mut cursor, &chunks).unwrap();
        buf.seek_virtual(VirtualOffset::new(0, 0)).unwrap();

        let mut scratch = Vec::new();
        let result = buf.read_record(&mut scratch);
        assert!(result.is_err(), "read_record should reject block_size > 2MB");
        let err = result.unwrap_err();
        assert!(
            matches!(err, BgzfError::RecordTooLarge { .. }),
            "expected RecordTooLarge, got {err:?}"
        );
    }

    #[test]
    fn merge_overlapping_chunks() {
        let chunks = vec![
            Chunk { begin: VirtualOffset::new(100, 0), end: VirtualOffset::new(200, 0) },
            Chunk { begin: VirtualOffset::new(150, 0), end: VirtualOffset::new(300, 0) },
        ];
        let ranges = merge_chunks(&chunks);
        assert_eq!(ranges.len(), 1);
        assert_eq!(ranges[0].file_start, 100);
    }

    #[test]
    fn keep_disjoint_chunks() {
        let chunks = vec![
            Chunk { begin: VirtualOffset::new(100, 0), end: VirtualOffset::new(200, 0) },
            Chunk { begin: VirtualOffset::new(500_000, 0), end: VirtualOffset::new(600_000, 0) },
        ];
        let ranges = merge_chunks(&chunks);
        assert_eq!(ranges.len(), 2);
    }

    #[test]
    fn empty_chunks_is_eof() {
        let mut cursor = std::io::Cursor::new(vec![]);
        let buf = RegionBuf::load(&mut cursor, &[]).unwrap();
        assert!(buf.eof);
    }

    // --- Disjoint range regression tests ---

    /// Helper: build fake file data of `len` bytes (just sequential bytes).
    fn fake_file(len: usize) -> Vec<u8> {
        (0..len).map(|i| (i % 256) as u8).collect()
    }

    #[test]
    fn single_range_seek_uses_simple_offset() {
        let file_data = fake_file(1000);
        let mut cursor = std::io::Cursor::new(file_data);

        let chunks =
            vec![Chunk { begin: VirtualOffset::new(100, 0), end: VirtualOffset::new(200, 0) }];

        let mut buf = RegionBuf::load(&mut cursor, &chunks).unwrap();

        // Seek to file offset 100 → cursor should be 0 (start of buffer)
        buf.seek_virtual(VirtualOffset::new(100, 0)).unwrap();
        assert_eq!(buf.cursor, 0);

        // Seek to file offset 150 → cursor should be 50
        buf.seek_virtual(VirtualOffset::new(150, 0)).unwrap();
        assert_eq!(buf.cursor, 50);
    }

    #[test]
    fn disjoint_ranges_data_loaded_correctly() {
        // With MAX_BLOCK_SIZE=65536, chunks within 65536 of each other merge.
        // Use offsets > 65536 apart.
        let chunks = vec![
            Chunk { begin: VirtualOffset::new(100, 0), end: VirtualOffset::new(200, 0) },
            Chunk { begin: VirtualOffset::new(100_000, 0), end: VirtualOffset::new(100_100, 0) },
        ];
        let ranges = merge_chunks(&chunks);
        assert_eq!(ranges.len(), 2, "chunks should be disjoint");

        // Now load from a file big enough to contain both ranges
        let big_file = fake_file(200_000);
        let mut cursor = std::io::Cursor::new(big_file.clone());
        let buf = RegionBuf::load(&mut cursor, &chunks).unwrap();

        // Buffer should contain data from BOTH ranges concatenated.
        // Range1: file[100..100+65536+200] (with 64KB extension)
        // Range2: file[100000..100000+65536+100100]
        // Total: both ranges' bytes.
        assert!(!buf.data.is_empty());

        // Verify first range's data starts with file[100]
        assert_eq!(buf.data[0], big_file[100]);
        assert_eq!(buf.data[1], big_file[101]);
    }

    #[test]
    fn disjoint_ranges_seek_to_second_range_is_correct() {
        // This is the core regression test.
        // Two disjoint ranges: [100, 300) and [100000, 100200).
        // After load, seeking to file offset 100000 must land at the
        // correct position in the buffer (start of second range's data),
        // NOT at offset 100000 - 100 = 99900 (which would be garbage).

        let big_file = fake_file(200_000);
        let mut cursor = std::io::Cursor::new(big_file.clone());

        let chunks = vec![
            Chunk { begin: VirtualOffset::new(100, 0), end: VirtualOffset::new(200, 0) },
            Chunk { begin: VirtualOffset::new(100_000, 0), end: VirtualOffset::new(100_100, 0) },
        ];
        let ranges = merge_chunks(&chunks);
        assert_eq!(ranges.len(), 2, "should have 2 disjoint ranges");

        let range1_len = (ranges[0].file_end - ranges[0].file_start) as usize;

        let mut buf = RegionBuf::load(&mut cursor, &chunks).unwrap();

        // Seek to start of first range
        buf.seek_virtual(VirtualOffset::new(100, 0)).unwrap();
        assert_eq!(buf.cursor, 0, "seek to range1 start");

        // The data at cursor=0 should be file[100]
        assert_eq!(buf.data[buf.cursor], big_file[100]);

        // Seek to start of second range — this is the bug
        buf.seek_virtual(VirtualOffset::new(100_000, 0)).unwrap();

        // cursor should point to the start of range2's data in the buffer,
        // which is at offset range1_len (after range1's data).
        assert_eq!(
            buf.cursor, range1_len,
            "seek to range2 should land at buffer offset {range1_len}, got {}",
            buf.cursor
        );

        // And the data there should be file[100000]
        assert_eq!(buf.data[buf.cursor], big_file[100_000]);
    }

    #[test]
    fn disjoint_ranges_seek_within_second_range() {
        let big_file = fake_file(200_000);
        let mut cursor = std::io::Cursor::new(big_file.clone());

        let chunks = vec![
            Chunk { begin: VirtualOffset::new(100, 0), end: VirtualOffset::new(200, 0) },
            Chunk { begin: VirtualOffset::new(100_000, 0), end: VirtualOffset::new(100_100, 0) },
        ];
        let ranges = merge_chunks(&chunks);
        let range1_len = (ranges[0].file_end - ranges[0].file_start) as usize;

        let mut buf = RegionBuf::load(&mut cursor, &chunks).unwrap();

        // Seek to file offset 100050 (50 bytes into second range)
        buf.seek_virtual(VirtualOffset::new(100_050, 0)).unwrap();
        assert_eq!(buf.cursor, range1_len + 50, "seek 50 bytes into range2");
        assert_eq!(buf.data[buf.cursor], big_file[100_050]);
    }

    #[test]
    fn disjoint_ranges_seek_before_loaded_region_fails() {
        let big_file = fake_file(200_000);
        let mut cursor = std::io::Cursor::new(big_file);

        let chunks =
            vec![Chunk { begin: VirtualOffset::new(1000, 0), end: VirtualOffset::new(2000, 0) }];

        let mut buf = RegionBuf::load(&mut cursor, &chunks).unwrap();

        // Seeking to file offset 500 (before the loaded range) must fail
        let result = buf.seek_virtual(VirtualOffset::new(500, 0));
        assert!(result.is_err(), "seek before loaded region should fail");
    }

    #[test]
    fn disjoint_ranges_seek_in_gap_between_ranges_fails() {
        let big_file = fake_file(200_000);
        let mut cursor = std::io::Cursor::new(big_file);

        let chunks = vec![
            Chunk { begin: VirtualOffset::new(100, 0), end: VirtualOffset::new(200, 0) },
            Chunk { begin: VirtualOffset::new(100_000, 0), end: VirtualOffset::new(100_100, 0) },
        ];

        let ranges = merge_chunks(&chunks);
        // Range1 ends at 200 + 65536 = 65736. Range2 starts at 100000.
        // The gap is [65736, 100000). Pick an offset in that gap.
        let gap_offset = ranges[0].file_end + 1;
        assert!(
            gap_offset < ranges[1].file_start,
            "gap_offset {gap_offset} should be before range2 start {}",
            ranges[1].file_start
        );

        let mut buf = RegionBuf::load(&mut cursor, &chunks).unwrap();

        let result = buf.seek_virtual(VirtualOffset::new(gap_offset, 0));
        assert!(result.is_err(), "seek into gap between ranges should fail");
    }

    // r[verify region_buf.decompress]
    #[test]
    fn read_block_rejects_corrupt_isize_in_footer() {
        // Build a valid BGZF block, then corrupt the ISIZE footer field
        // to claim a huge uncompressed size (e.g. ~4 GB). This must be
        // rejected with UncompressedSizeTooLarge, not cause an OOM.
        let payload = vec![0u8; 32];
        let mut block = make_bgzf_block(&payload);

        // The ISIZE field is the last 4 bytes of the BGZF block.
        // Overwrite it with a value larger than MAX_BLOCK_SIZE (65536).
        let isize_offset = block.len() - 4;
        let corrupt_isize: u32 = 0xF0F0_F0F0; // ~4 GB
        block[isize_offset..isize_offset + 4].copy_from_slice(&corrupt_isize.to_le_bytes());

        let mut file = Vec::new();
        file.extend_from_slice(&block);
        file.extend_from_slice(&make_bgzf_eof());

        let chunks = vec![Chunk { begin: VirtualOffset::new(0, 0), end: VirtualOffset::new(1, 0) }];

        let mut cursor = std::io::Cursor::new(file);
        let mut buf = RegionBuf::load(&mut cursor, &chunks).unwrap();
        buf.seek_virtual(VirtualOffset::new(0, 0)).unwrap();

        // read_block is called internally by seek_virtual (via within=0 path)
        // or by read_exact_into. Force a read_block via read_exact_into.
        let mut out = [0u8; 1];
        let result = buf.read_exact_into(&mut out);
        assert!(result.is_err(), "should reject corrupt ISIZE");
        let err = result.unwrap_err();
        assert!(
            matches!(err, BgzfError::UncompressedSizeTooLarge { .. }),
            "expected UncompressedSizeTooLarge, got {err:?}"
        );
    }

    // --- BGZF round-trip proptests ---

    /// Build a single BGZF block from uncompressed data. Returns the full
    /// block bytes (header + compressed payload + footer).
    fn make_bgzf_block(data: &[u8]) -> Vec<u8> {
        let mut compressor =
            libdeflater::Compressor::new(libdeflater::CompressionLvl::new(1).unwrap());
        let bound = compressor.deflate_compress_bound(data.len());
        let mut compressed = vec![0u8; bound];
        let compressed_len =
            compressor.deflate_compress(data, &mut compressed).expect("compression");
        compressed.truncate(compressed_len);

        let mut crc = libdeflater::Crc::new();
        crc.update(data);

        // BGZF header (18 bytes)
        let bsize = (18 + compressed_len + 8 - 1) as u16; // total block size - 1
        let mut block = Vec::with_capacity(18 + compressed_len + 8);
        // gzip magic + DEFLATE + FEXTRA
        block.extend_from_slice(&[0x1f, 0x8b, 0x08, 0x04]);
        block.extend_from_slice(&[0; 4]); // MTIME
        block.push(0); // XFL
        block.push(0xff); // OS
        block.extend_from_slice(&6u16.to_le_bytes()); // XLEN = 6
        block.extend_from_slice(&[b'B', b'C', 2, 0]); // BC subfield, SLEN=2
        block.extend_from_slice(&bsize.to_le_bytes()); // BSIZE

        // Compressed data
        block.extend_from_slice(&compressed);

        // Footer: CRC32 + ISIZE
        block.extend_from_slice(&crc.sum().to_le_bytes());
        block.extend_from_slice(&(data.len() as u32).to_le_bytes());

        block
    }

    /// Build a BGZF EOF marker block.
    fn make_bgzf_eof() -> Vec<u8> {
        vec![
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        ]
    }

    /// Build a fake BGZF file with N blocks, each containing `block_data[i]`.
    /// Returns (file_bytes, block_offsets) where block_offsets[i] is the
    /// file offset of block i.
    fn make_bgzf_file(blocks: &[Vec<u8>]) -> (Vec<u8>, Vec<u64>) {
        let mut file = Vec::new();
        let mut offsets = Vec::with_capacity(blocks.len());

        for block_data in blocks {
            offsets.push(file.len() as u64);
            file.extend_from_slice(&make_bgzf_block(block_data));
        }
        // EOF block
        file.extend_from_slice(&make_bgzf_eof());

        (file, offsets)
    }

    use proptest::prelude::*;

    proptest! {
        /// Single contiguous range: load all blocks, read them sequentially,
        /// verify decompressed content matches original.
        #[test]
        fn proptest_single_range_roundtrip(
            n_blocks in 1usize..8,
            block_size in 10usize..500,
            seed in 0u8..255,
        ) {
            // Build blocks with deterministic content
            let blocks: Vec<Vec<u8>> = (0..n_blocks)
                .map(|i| {
                    (0..block_size)
                        .map(|j| seed.wrapping_add(i as u8).wrapping_add(j as u8))
                        .collect()
                })
                .collect();

            let (file, offsets) = make_bgzf_file(&blocks);
            let last_offset = *offsets.last().unwrap();

            // One chunk covering all blocks
            let chunks = vec![Chunk {
                begin: VirtualOffset::new(offsets[0], 0),
                end: VirtualOffset::new(last_offset + 1, 0),
            }];

            let mut cursor = std::io::Cursor::new(file);
            let mut buf = RegionBuf::load(&mut cursor, &chunks).unwrap();

            // Seek to the first block and read all data
            buf.seek_virtual(VirtualOffset::new(offsets[0], 0)).unwrap();

            let total_bytes: usize = blocks.iter().map(|b| b.len()).sum();
            let mut output = vec![0u8; total_bytes];
            buf.read_exact_into(&mut output).unwrap();

            // Verify
            let expected: Vec<u8> = blocks.iter().flatten().copied().collect();
            prop_assert_eq!(output, expected);
        }

        /// Disjoint ranges: create blocks in two groups separated by padding,
        /// load both groups, seek to each and verify content.
        #[test]
        fn proptest_disjoint_ranges_roundtrip(
            n_blocks_a in 1usize..4,
            n_blocks_b in 1usize..4,
            block_size in 10usize..300,
            padding in 100_000usize..200_000,
            seed in 0u8..255,
        ) {
            let blocks_a: Vec<Vec<u8>> = (0..n_blocks_a)
                .map(|i| {
                    (0..block_size)
                        .map(|j| seed.wrapping_add(i as u8).wrapping_add(j as u8))
                        .collect()
                })
                .collect();
            let blocks_b: Vec<Vec<u8>> = (0..n_blocks_b)
                .map(|i| {
                    (0..block_size)
                        .map(|j| seed.wrapping_add(100).wrapping_add(i as u8).wrapping_add(j as u8))
                        .collect()
                })
                .collect();

            // Build group A
            let (mut file, offsets_a) = make_bgzf_file(&blocks_a);

            // Add padding to create a gap > MAX_BLOCK_SIZE so chunks are disjoint
            let pad_start = file.len();
            file.resize(pad_start + padding, 0);

            // Build group B at the padded offset
            let group_b_start = file.len() as u64;
            let mut offsets_b = Vec::new();
            for block_data in &blocks_b {
                offsets_b.push(file.len() as u64);
                file.extend_from_slice(&make_bgzf_block(block_data));
            }
            file.extend_from_slice(&make_bgzf_eof());

            // Two disjoint chunks
            let last_a = *offsets_a.last().unwrap();
            let last_b = *offsets_b.last().unwrap();
            let chunk_a = Chunk {
                begin: VirtualOffset::new(offsets_a[0], 0),
                end: VirtualOffset::new(last_a + 1, 0),
            };
            let chunk_b = Chunk {
                begin: VirtualOffset::new(group_b_start, 0),
                end: VirtualOffset::new(last_b + 1, 0),
            };

            let ranges = merge_chunks(&[chunk_a, chunk_b]);
            prop_assert_eq!(ranges.len(), 2);

            let mut cursor = std::io::Cursor::new(file);
            let mut buf = RegionBuf::load(&mut cursor, &[chunk_a, chunk_b]).unwrap();

            // Read group A
            buf.seek_virtual(VirtualOffset::new(offsets_a[0], 0)).unwrap();
            let total_a: usize = blocks_a.iter().map(|b| b.len()).sum();
            let mut out_a = vec![0u8; total_a];
            buf.read_exact_into(&mut out_a).unwrap();
            let expected_a: Vec<u8> = blocks_a.iter().flatten().copied().collect();
            prop_assert_eq!(out_a, expected_a, "group A content mismatch");

            // Read group B
            buf.seek_virtual(VirtualOffset::new(group_b_start, 0)).unwrap();
            let total_b: usize = blocks_b.iter().map(|b| b.len()).sum();
            let mut out_b = vec![0u8; total_b];
            buf.read_exact_into(&mut out_b).unwrap();
            let expected_b: Vec<u8> = blocks_b.iter().flatten().copied().collect();
            prop_assert_eq!(out_b, expected_b, "group B content mismatch");
        }

        /// Seek to mid-block positions (within_block > 0) should work correctly.
        #[test]
        fn proptest_within_block_seek(
            block_size in 20usize..500,
            within in 1usize..19, // at most block_size-1, capped at 19 for simplicity
            seed in 0u8..255,
        ) {
            let data: Vec<u8> = (0..block_size)
                .map(|j| seed.wrapping_add(j as u8))
                .collect();

            let (file, offsets) = make_bgzf_file(std::slice::from_ref(&data));

            let chunks = vec![Chunk {
                begin: VirtualOffset::new(offsets[0], 0),
                end: VirtualOffset::new(offsets[0] + 1, 0),
            }];

            let mut cursor = std::io::Cursor::new(file);
            let mut buf = RegionBuf::load(&mut cursor, &chunks).unwrap();

            let within_clamped = within.min(block_size - 1);
            buf.seek_virtual(VirtualOffset::new(offsets[0], within_clamped as u16)).unwrap();

            let remaining = block_size - within_clamped;
            let mut output = vec![0u8; remaining];
            buf.read_exact_into(&mut output).unwrap();

            prop_assert_eq!(output, data[within_clamped..].to_vec());
        }

        /// CRC32 mismatch detection: corrupt a byte in the compressed data
        /// and verify decompression or CRC check fails.
        #[test]
        fn proptest_crc32_detects_corruption(
            block_size in 20usize..200,
            seed in 0u8..255,
        ) {
            let data: Vec<u8> = (0..block_size)
                .map(|j| seed.wrapping_add(j as u8))
                .collect();

            let (mut file, offsets) = make_bgzf_file(&[data]);

            // Corrupt a byte in the compressed payload (after the 18-byte header)
            let corrupt_pos = offsets[0] as usize + 18;
            if corrupt_pos < file.len() - 8 {
                file[corrupt_pos] ^= 0xFF;

                let chunks = vec![Chunk {
                    begin: VirtualOffset::new(offsets[0], 0),
                    end: VirtualOffset::new(offsets[0] + 1, 0),
                }];

                let mut cursor = std::io::Cursor::new(file);
                let mut buf = RegionBuf::load(&mut cursor, &chunks).unwrap();
                buf.seek_virtual(VirtualOffset::new(offsets[0], 0)).unwrap();

                let mut output = vec![0u8; block_size];
                let result = buf.read_exact_into(&mut output);
                // Should fail with either DecompressionFailed or ChecksumMismatch
                prop_assert!(result.is_err(), "corrupted block should fail");
            }
        }
    }

    // --- read_record tests ---

    /// Encode a BAM-style length-prefixed record into a byte vec:
    /// 4 bytes LE u32 `body.len()` followed by `body`.
    fn make_length_prefixed(body: &[u8]) -> Vec<u8> {
        let mut out = Vec::with_capacity(4 + body.len());
        out.extend_from_slice(&(body.len() as u32).to_le_bytes());
        out.extend_from_slice(body);
        out
    }

    /// Fast path: length prefix AND body both fit inside one BGZF block.
    /// `read_record` must return a slice pointing directly into `RegionBuf::buf`
    /// (zero-copy), and scratch must remain empty.
    #[test]
    fn read_record_single_block_fast_path() {
        let body: Vec<u8> = (0u8..64).collect();
        let block_payload = make_length_prefixed(&body);

        let (file, offsets) = make_bgzf_file(&[block_payload]);
        let chunks = vec![Chunk {
            begin: VirtualOffset::new(offsets[0], 0),
            end: VirtualOffset::new(offsets[0] + 1, 0),
        }];

        let mut cursor = std::io::Cursor::new(file);
        let mut region = RegionBuf::load(&mut cursor, &chunks).unwrap();
        region.seek_virtual(VirtualOffset::new(offsets[0], 0)).unwrap();

        let mut scratch: Vec<u8> = Vec::new();
        let result = region.read_record(&mut scratch).unwrap();

        assert_eq!(result, body.as_slice());
        // scratch untouched — fast path never fills it
        assert!(scratch.is_empty());
    }

    /// Slow path: the record body straddles a BGZF block boundary.
    /// Block 1 ends with the 4-byte length prefix; block 2 holds the body.
    /// `read_record` must fall back to the scratch buffer and still return
    /// the correct bytes.
    #[test]
    fn read_record_cross_block_slow_path() {
        let body: Vec<u8> = (0u8..32).map(|b| b.wrapping_mul(3)).collect();
        let len_bytes = (body.len() as u32).to_le_bytes();

        // Block 1: 60 filler bytes then the 4-byte length prefix (total 64 bytes).
        let mut block1_data = vec![0xffu8; 60];
        block1_data.extend_from_slice(&len_bytes);

        // Block 2: just the record body.
        let block2_data = body.clone();

        let (file, offsets) = make_bgzf_file(&[block1_data, block2_data]);
        let chunks = vec![Chunk {
            begin: VirtualOffset::new(offsets[0], 0),
            end: VirtualOffset::new(offsets[1] + 1, 0),
        }];

        let mut cursor = std::io::Cursor::new(file);
        let mut region = RegionBuf::load(&mut cursor, &chunks).unwrap();
        // Seek past the filler to where the length prefix starts (within_block = 60).
        region.seek_virtual(VirtualOffset::new(offsets[0], 60)).unwrap();

        let mut scratch: Vec<u8> = Vec::new();
        let result = region.read_record(&mut scratch).unwrap();

        assert_eq!(result, body.as_slice());
        // Slow path: scratch was used.
        assert_eq!(scratch.as_slice(), body.as_slice());
    }

    /// When the length prefix itself straddles a block boundary, `read_record`
    /// must still decode it correctly and return the right body.
    #[test]
    fn read_record_length_prefix_straddles_block_boundary() {
        let body: Vec<u8> = (0u8..16).collect();
        let len_bytes = (body.len() as u32).to_le_bytes();

        // Block 1: 62 filler bytes then the first 2 bytes of the 4-byte prefix.
        let mut block1_data = vec![0xaau8; 62];
        block1_data.extend_from_slice(&len_bytes[..2]);

        // Block 2: remaining 2 bytes of prefix then the body.
        let mut block2_data = Vec::new();
        block2_data.extend_from_slice(&len_bytes[2..]);
        block2_data.extend_from_slice(&body);

        let (file, offsets) = make_bgzf_file(&[block1_data, block2_data]);
        let chunks = vec![Chunk {
            begin: VirtualOffset::new(offsets[0], 0),
            end: VirtualOffset::new(offsets[1] + 1, 0),
        }];

        let mut cursor = std::io::Cursor::new(file);
        let mut region = RegionBuf::load(&mut cursor, &chunks).unwrap();
        // Seek to where the first 2 bytes of the prefix begin (within_block = 62).
        region.seek_virtual(VirtualOffset::new(offsets[0], 62)).unwrap();

        let mut scratch: Vec<u8> = Vec::new();
        let result = region.read_record(&mut scratch).unwrap();

        assert_eq!(result, body.as_slice());
    }

    /// Multiple sequential records, all within a single BGZF block.
    /// Each call to `read_record` must advance the position and return
    /// the correct body independently.
    #[test]
    fn read_record_multiple_sequential_records_single_block() {
        let records: Vec<Vec<u8>> =
            vec![(0u8..8).collect(), (10u8..26).collect(), vec![0xdeu8, 0xad, 0xbe, 0xef]];

        let mut block_payload = Vec::new();
        for rec in &records {
            block_payload.extend_from_slice(&make_length_prefixed(rec));
        }

        let (file, offsets) = make_bgzf_file(&[block_payload]);
        let chunks = vec![Chunk {
            begin: VirtualOffset::new(offsets[0], 0),
            end: VirtualOffset::new(offsets[0] + 1, 0),
        }];

        let mut cursor = std::io::Cursor::new(file);
        let mut region = RegionBuf::load(&mut cursor, &chunks).unwrap();
        region.seek_virtual(VirtualOffset::new(offsets[0], 0)).unwrap();

        let mut scratch = Vec::new();
        for expected in &records {
            let got = region.read_record(&mut scratch).unwrap();
            assert_eq!(got, expected.as_slice());
        }
    }

    /// Proptest: varied record counts and body sizes, all packed into a single
    /// BGZF block. Every record must round-trip correctly.
    #[test]
    fn proptest_read_record_roundtrip() {
        use proptest::prelude::*;

        proptest!(|(
            n_records in 1usize..10,
            body_size in 4usize..200,
            seed in 0u8..255,
        )| {
            let records: Vec<Vec<u8>> = (0..n_records)
                .map(|i| {
                    (0..body_size)
                        .map(|j| seed.wrapping_add(i as u8).wrapping_add(j as u8))
                        .collect()
                })
                .collect();

            let mut block_payload = Vec::new();
            for rec in &records {
                block_payload.extend_from_slice(&make_length_prefixed(rec));
            }

            let (file, offsets) = make_bgzf_file(&[block_payload]);
            let chunks = vec![Chunk {
                begin: VirtualOffset::new(offsets[0], 0),
                end: VirtualOffset::new(offsets[0] + 1, 0),
            }];

            let mut cursor = std::io::Cursor::new(file);
            let mut region = RegionBuf::load(&mut cursor, &chunks).unwrap();
            region.seek_virtual(VirtualOffset::new(offsets[0], 0)).unwrap();

            let mut scratch = Vec::new();
            for (i, expected) in records.iter().enumerate() {
                let got = region.read_record(&mut scratch)
                    .unwrap_or_else(|e| panic!("record {i} failed: {e}"));
                prop_assert_eq!(got, expected.as_slice());
            }
        });
    }

    /// Zero-length body: a record with block_size=0 should be read back as
    /// an empty slice without error.
    #[test]
    fn read_record_zero_length_body() {
        let block_payload = make_length_prefixed(&[]);

        let (file, offsets) = make_bgzf_file(&[block_payload]);
        let chunks = vec![Chunk {
            begin: VirtualOffset::new(offsets[0], 0),
            end: VirtualOffset::new(offsets[0] + 1, 0),
        }];

        let mut cursor = std::io::Cursor::new(file);
        let mut region = RegionBuf::load(&mut cursor, &chunks).unwrap();
        region.seek_virtual(VirtualOffset::new(offsets[0], 0)).unwrap();

        let mut scratch = Vec::new();
        let result = region.read_record(&mut scratch).unwrap();
        assert_eq!(result, &[] as &[u8]);
    }

    /// A truncated stream (length prefix present but body cut short) must
    /// return an error rather than silently returning wrong data.
    #[test]
    fn read_record_truncated_body_returns_error() {
        let body_len: u32 = 32;
        let mut block_payload = Vec::new();
        block_payload.extend_from_slice(&body_len.to_le_bytes());
        // Only write 10 bytes of the promised 32.
        block_payload.extend_from_slice(&[0u8; 10]);

        let (file, offsets) = make_bgzf_file(&[block_payload]);
        let chunks = vec![Chunk {
            begin: VirtualOffset::new(offsets[0], 0),
            end: VirtualOffset::new(offsets[0] + 1, 0),
        }];

        let mut cursor = std::io::Cursor::new(file);
        let mut region = RegionBuf::load(&mut cursor, &chunks).unwrap();
        region.seek_virtual(VirtualOffset::new(offsets[0], 0)).unwrap();

        let mut scratch = Vec::new();
        assert!(region.read_record(&mut scratch).is_err());
    }
}
