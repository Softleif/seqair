//! BGZF block writer. [`BgzfWriter`] compresses data into independent BGZF blocks
//! via libdeflate, tracking virtual offsets for index co-production.

use super::bgzf::{BgzfError, VirtualOffset};
use std::io::Write;
use tracing::warn;

/// Maximum uncompressed data per BGZF block (64 KiB).
const MAX_UNCOMPRESSED_SIZE: usize = 65536;

/// BGZF header template: gzip magic + DEFLATE + FEXTRA, then BC subfield.
/// Bytes 16-17 (BSIZE) are filled per block.
const BGZF_HEADER: [u8; 18] = [
    0x1f, 0x8b, // gzip magic
    0x08, // CM = DEFLATE
    0x04, // FLG = FEXTRA
    0x00, 0x00, 0x00, 0x00, // MTIME
    0x00, // XFL
    0xff, // OS = unknown
    0x06, 0x00, // XLEN = 6
    0x42, 0x43, // SI1, SI2 = 'B', 'C'
    0x02, 0x00, // SLEN = 2
    0x00, 0x00, // BSIZE placeholder (filled per block)
];

/// Standard 28-byte BGZF EOF marker (empty gzip member with ISIZE=0).
const EOF_BLOCK: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, // BSIZE = 27 (total size 28 - 1)
    0x03, 0x00, // empty DEFLATE block
    0x00, 0x00, 0x00, 0x00, // CRC32 = 0
    0x00, 0x00, 0x00, 0x00, // ISIZE = 0
];

// r[impl bgzf.writer]
// r[impl bgzf.writer.buffer]
// r[impl bgzf.writer.compression]
// r[impl bgzf.writer.virtual_offset]
/// BGZF block writer that compresses data into independent gzip blocks.
///
/// Accumulates up to 64 KB of uncompressed data, then compresses and emits a
/// BGZF block. Tracks virtual offsets for index co-production.
pub struct BgzfWriter<W: Write> {
    inner: Option<W>,
    /// Uncompressed data buffer (up to MAX_UNCOMPRESSED_SIZE).
    buf: Vec<u8>,
    /// Reusable buffer for compressed output.
    compressed_buf: Vec<u8>,
    compressor: libdeflater::Compressor,
    /// Compressed file offset of the current (not yet flushed) block.
    block_offset: u64,
}

impl<W: Write> BgzfWriter<W> {
    /// Create a new BGZF writer with default compression level (6).
    pub fn new(inner: W) -> Self {
        Self::with_compression_level(inner, 6)
    }

    /// Create a new BGZF writer with the specified compression level (0-12).
    pub fn with_compression_level(inner: W, level: i32) -> Self {
        Self {
            inner: Some(inner),
            buf: Vec::with_capacity(MAX_UNCOMPRESSED_SIZE),
            compressed_buf: Vec::with_capacity(MAX_UNCOMPRESSED_SIZE),
            compressor: libdeflater::Compressor::new(
                libdeflater::CompressionLvl::new(level).unwrap_or_default(),
            ),
            block_offset: 0,
        }
    }

    /// Get a mutable reference to the inner writer.
    fn writer(&mut self) -> Result<&mut W, BgzfError> {
        self.inner.as_mut().ok_or(BgzfError::AlreadyFinished)
    }

    // r[impl bgzf.writer.virtual_offset]
    /// Current write position as a virtual offset.
    ///
    /// Upper 48 bits = compressed offset of the current block.
    /// Lower 16 bits = bytes written into the current (unflushed) block.
    pub fn virtual_offset(&self) -> VirtualOffset {
        // Buffer can be at most MAX_UNCOMPRESSED_SIZE (65536) which doesn't fit u16.
        // When buffer is exactly full, flush_block hasn't run yet — cap at 65535.
        // This is safe: a full buffer will be flushed on the next write_all/flush_if_needed.
        let within = self.buf.len().min(u16::MAX as usize) as u16;
        VirtualOffset::new(self.block_offset, within)
    }

    // r[impl bgzf.writer.flush_if_needed]
    /// Flush the current block if `upcoming_bytes` would exceed the 64 KB limit.
    ///
    /// Call this before writing a record to keep it from spanning block boundaries,
    /// improving seek granularity for index-based random access.
    pub fn flush_if_needed(&mut self, upcoming_bytes: usize) -> Result<(), BgzfError> {
        if self.buf.len().saturating_add(upcoming_bytes) > MAX_UNCOMPRESSED_SIZE {
            self.flush_block()?;
        }
        Ok(())
    }

    /// Write data into the BGZF stream. Flushes blocks as needed.
    pub fn write_all(&mut self, data: &[u8]) -> Result<(), BgzfError> {
        let mut remaining = data;
        while !remaining.is_empty() {
            let space = MAX_UNCOMPRESSED_SIZE.saturating_sub(self.buf.len());
            if space == 0 {
                self.flush_block()?;
                continue;
            }
            let take = remaining.len().min(space);
            // Safe: take <= remaining.len()
            #[allow(clippy::indexing_slicing)]
            {
                self.buf.extend_from_slice(&remaining[..take]);
                remaining = &remaining[take..];
            }
        }
        Ok(())
    }

    /// Compress and emit the current buffer as a BGZF block.
    fn flush_block(&mut self) -> Result<(), BgzfError> {
        if self.buf.is_empty() {
            return Ok(());
        }

        // Compress the data
        let max_compressed = self.compressor.deflate_compress_bound(self.buf.len());
        self.compressed_buf.clear();
        self.compressed_buf.resize(max_compressed, 0);
        let compressed_len = self
            .compressor
            .deflate_compress(&self.buf, &mut self.compressed_buf)
            .map_err(|source| BgzfError::CompressionFailed { source })?;
        self.compressed_buf.truncate(compressed_len);

        // Compute CRC32
        let mut crc = libdeflater::Crc::new();
        crc.update(&self.buf);
        let crc32 = crc.sum();

        // BSIZE = total block size - 1
        // total = header(18) + compressed_len + footer(8)
        let total_block_size = 18usize
            .checked_add(compressed_len)
            .and_then(|n| n.checked_add(8))
            .ok_or(BgzfError::CorruptHeader)?;
        let bsize = total_block_size.checked_sub(1).ok_or(BgzfError::CorruptHeader)?;

        // Write header with BSIZE filled in
        let mut header = BGZF_HEADER;
        let bsize_bytes = u16::try_from(bsize).map_err(|_| BgzfError::CorruptHeader)?.to_le_bytes();
        // Bytes 16-17 are the BSIZE field
        #[allow(clippy::indexing_slicing)]
        {
            header[16] = bsize_bytes[0];
            header[17] = bsize_bytes[1];
        }

        let isize_val = u32::try_from(self.buf.len()).map_err(|_| BgzfError::CorruptHeader)?;

        // Borrow inner writer directly to avoid conflicting borrow with self.compressed_buf/buf
        let w = self.inner.as_mut().ok_or(BgzfError::AlreadyFinished)?;
        w.write_all(&header).map_err(|source| BgzfError::WriteFailed { source })?;
        w.write_all(&self.compressed_buf).map_err(|source| BgzfError::WriteFailed { source })?;
        w.write_all(&crc32.to_le_bytes()).map_err(|source| BgzfError::WriteFailed { source })?;
        w.write_all(&isize_val.to_le_bytes())
            .map_err(|source| BgzfError::WriteFailed { source })?;

        // Advance block_offset by the total compressed block size
        self.block_offset = self
            .block_offset
            .checked_add(total_block_size as u64)
            .ok_or(BgzfError::CorruptHeader)?;
        self.buf.clear();

        Ok(())
    }

    // r[impl bgzf.writer.eof_marker]
    // r[impl bgzf.writer.finish]
    /// Flush remaining data, write the EOF marker, and return the inner writer.
    pub fn finish(mut self) -> Result<W, BgzfError> {
        self.flush_block()?;
        let w = self.writer()?;
        w.write_all(&EOF_BLOCK).map_err(|source| BgzfError::WriteFailed { source })?;
        w.flush().map_err(|source| BgzfError::WriteFailed { source })?;
        // Take the inner writer so Drop doesn't try to flush again
        self.inner.take().ok_or(BgzfError::AlreadyFinished)
    }
}

impl<W: Write> Drop for BgzfWriter<W> {
    fn drop(&mut self) {
        // Best-effort flush on drop — log errors since we can't propagate from Drop.
        if !self.buf.is_empty()
            && self.inner.is_some()
            && let Err(e) = self.flush_block()
        {
            warn!("BgzfWriter::drop: failed to flush {} buffered bytes: {e}", self.buf.len());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam::bgzf::BgzfReader;
    use std::io::Cursor;

    fn write_and_finish(data: &[u8]) -> Vec<u8> {
        let mut output = Vec::new();
        let mut writer = BgzfWriter::new(&mut output);
        writer.write_all(data).unwrap();
        writer.finish().unwrap();
        output
    }

    fn read_all(compressed: &[u8]) -> Vec<u8> {
        let mut reader = BgzfReader::from_reader(Cursor::new(compressed));
        let mut result = Vec::new();
        reader.read_to_end(&mut result).unwrap();
        result
    }

    // r[verify bgzf.writer]
    // r[verify bgzf.writer.buffer]
    #[test]
    fn round_trip_small() {
        let data = b"Hello, BGZF world!";
        let output = write_and_finish(data);
        assert_eq!(read_all(&output), data);
    }

    // r[verify bgzf.writer.buffer]
    #[test]
    fn round_trip_exact_block_boundary() {
        let data: Vec<u8> = (0..MAX_UNCOMPRESSED_SIZE).map(|i| (i & 0xFF) as u8).collect();
        let output = write_and_finish(&data);
        assert_eq!(read_all(&output), data);
    }

    // r[verify bgzf.writer.buffer]
    #[test]
    fn round_trip_spanning_multiple_blocks() {
        let data: Vec<u8> = (0..200_000).map(|i| (i % 251) as u8).collect();
        let output = write_and_finish(&data);
        assert_eq!(read_all(&output), data);
    }

    // r[verify bgzf.writer.eof_marker]
    #[test]
    fn eof_marker_present() {
        let mut output = Vec::new();
        let writer = BgzfWriter::new(&mut output);
        writer.finish().unwrap();
        assert!(output.len() >= 28);
        #[allow(clippy::indexing_slicing)]
        {
            assert_eq!(&output[output.len() - 28..], &EOF_BLOCK);
        }
    }

    // r[verify bgzf.writer.virtual_offset]
    #[test]
    fn virtual_offset_tracking() {
        let mut output = Vec::new();
        let mut writer = BgzfWriter::new(&mut output);

        let voff0 = writer.virtual_offset();
        assert_eq!(voff0.block_offset(), 0);
        assert_eq!(voff0.within_block(), 0);

        writer.write_all(b"test data").unwrap();
        let voff1 = writer.virtual_offset();
        assert_eq!(voff1.block_offset(), 0);
        assert_eq!(voff1.within_block(), 9);

        // Fill the block and write one more byte to force a flush
        let fill = vec![0u8; MAX_UNCOMPRESSED_SIZE - 9];
        writer.write_all(&fill).unwrap();
        // Buffer is full (64KB) but not yet flushed — flush happens on next write
        writer.write_all(&[0x42]).unwrap();
        let voff2 = writer.virtual_offset();
        assert!(voff2.block_offset() > 0, "should have advanced to a new block");
        assert_eq!(voff2.within_block(), 1, "one byte in the new block");

        writer.finish().unwrap();
    }

    // r[verify bgzf.writer.flush_if_needed]
    #[test]
    fn flush_if_needed_triggers_flush() {
        let mut output = Vec::new();
        let mut writer = BgzfWriter::new(&mut output);

        writer.write_all(&vec![0u8; 60_000]).unwrap();
        assert_eq!(writer.virtual_offset().block_offset(), 0);

        // 60K + 10K > 64K — should flush
        writer.flush_if_needed(10_000).unwrap();
        assert!(writer.virtual_offset().block_offset() > 0);
        assert_eq!(writer.virtual_offset().within_block(), 0);

        writer.finish().unwrap();
    }

    // r[verify bgzf.writer.flush_if_needed]
    #[test]
    fn flush_if_needed_no_flush_when_fits() {
        let mut output = Vec::new();
        let mut writer = BgzfWriter::new(&mut output);

        writer.write_all(&vec![0u8; 1000]).unwrap();
        let before = writer.virtual_offset();

        writer.flush_if_needed(100).unwrap();
        let after = writer.virtual_offset();
        assert_eq!(before, after, "should not have flushed");

        writer.finish().unwrap();
    }

    // r[verify bgzf.writer.compression]
    #[test]
    fn compression_level_configurable() {
        let data = b"compression test data repeated many times for better ratio ";
        let repeated: Vec<u8> = data.iter().copied().cycle().take(10_000).collect();

        let mut out1 = Vec::new();
        let mut w1 = BgzfWriter::with_compression_level(&mut out1, 1);
        w1.write_all(&repeated).unwrap();
        w1.finish().unwrap();

        let mut out9 = Vec::new();
        let mut w9 = BgzfWriter::with_compression_level(&mut out9, 9);
        w9.write_all(&repeated).unwrap();
        w9.finish().unwrap();

        assert!(out9.len() <= out1.len());
        assert_eq!(read_all(&out1), repeated);
        assert_eq!(read_all(&out9), repeated);
    }

    // r[verify bgzf.writer]
    #[test]
    fn empty_write_produces_only_eof() {
        let mut output = Vec::new();
        let writer = BgzfWriter::new(&mut output);
        writer.finish().unwrap();
        assert_eq!(output.len(), 28);
    }
}
