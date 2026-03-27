//! Read and seek BGZF files block-by-block. [`BgzfReader`] decompresses blocks via libdeflate
//! with CRC32 verification; [`VirtualOffset`] encodes block + intra-block positions.

use std::{
    io::{self, BufReader, Read, Seek, SeekFrom},
    path::Path,
};
use tracing::instrument;

// r[impl bgzf.virtual_offset]
/// Virtual file offset in a BGZF file.
///
/// Upper 48 bits: compressed block offset. Lower 16 bits: uncompressed offset
/// within the block.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct VirtualOffset(pub u64);

impl VirtualOffset {
    pub fn new(block_offset: u64, within_block: u16) -> Self {
        Self((block_offset << 16) | u64::from(within_block))
    }

    pub fn block_offset(self) -> u64 {
        self.0 >> 16
    }

    pub fn within_block(self) -> u16 {
        (self.0 & 0xFFFF) as u16
    }
}

// r[impl io.errors]
#[derive(Debug, thiserror::Error)]
pub enum BgzfError {
    #[error("I/O error opening {path}")]
    Open { path: std::path::PathBuf, source: std::io::Error },

    #[error("seek failed")]
    SeekFailed,

    #[error("invalid BGZF magic bytes")]
    InvalidMagic,

    #[error("BGZF extra subfields do not contain a BSIZE (BC) entry")]
    MissingBsize,

    #[error("BGZF BSIZE ({bsize}) is too small to account for the block header")]
    BlockSizeTooSmall { bsize: u16 },

    #[error("BGZF block data is truncated")]
    TruncatedBlock,

    #[error("BGZF decompression failed")]
    DecompressionFailed { source: libdeflater::DecompressionError },

    #[error("BGZF CRC32 mismatch: expected {expected:#010x}, got {found:#010x}")]
    ChecksumMismatch { expected: u32, found: u32 },

    #[error("unexpected EOF in BGZF stream")]
    UnexpectedEof,

    #[error("virtual offset {offset:#x} is not within any loaded BGZF range")]
    VirtualOffsetOutOfRange { offset: u64 },

    #[error("BGZF ISIZE ({isize_value}) exceeds maximum block size (65536)")]
    UncompressedSizeTooLarge { isize_value: usize },

    #[error("BAM record block_size ({block_size}) exceeds maximum ({max_size})")]
    RecordTooLarge { block_size: usize, max_size: usize },
}

// r[impl bgzf.libdeflate]
/// BGZF block reader that decompresses gzip blocks from a BAM file.
///
/// Tracks virtual file offsets for index-based seeking. Uses `libdeflater`
/// for decompression.
pub struct BgzfReader {
    inner: BufReader<std::fs::File>,
    buf: Vec<u8>,
    buf_pos: usize,
    // r[impl bgzf.block_offset_tracking]
    /// Tracked manually: updated on seek and after each block header read.
    block_offset: u64,
    eof: bool,
    // Reusable buffers to avoid per-block allocations
    compressed_buf: Vec<u8>,
    decompressor: libdeflater::Decompressor,
}

impl std::fmt::Debug for BgzfReader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BgzfReader")
            .field("block_offset", &self.block_offset)
            .field("buf_pos", &self.buf_pos)
            .field("eof", &self.eof)
            .finish()
    }
}

/// Full BGZF header: 12-byte gzip header + 6-byte BC subfield = 18 bytes.
const BGZF_HEADER_SIZE: usize = 18;
const BGZF_MAGIC: [u8; 4] = [0x1f, 0x8b, 0x08, 0x04]; // gzip magic + DEFLATE + FEXTRA
const BGZF_FOOTER_SIZE: usize = 8; // CRC32 (4) + ISIZE (4)
const MAX_BLOCK_SIZE: usize = 65536;

/// Resize a Vec<u8> without zero-filling.
///
/// # Safety
/// Caller must ensure all bytes in `0..new_len` are written before being read.
/// Safe because `u8` has no invalid bit patterns and `reserve_exact` ensures capacity.
#[inline(always)]
pub(crate) unsafe fn resize_uninit(buf: &mut Vec<u8>, new_len: usize) {
    buf.clear();
    buf.reserve_exact(new_len);
    debug_assert!(
        new_len <= buf.capacity(),
        "reserve_exact didn't provide enough capacity: need {new_len}, got {}",
        buf.capacity()
    );
    // Safety: u8 has no invalid bit patterns; caller guarantees all bytes written before read.
    unsafe { buf.set_len(new_len) };
}

impl BgzfReader {
    #[instrument(level = "debug", fields(path = %path.display()), err)]
    pub fn open(path: &Path) -> Result<Self, BgzfError> {
        let file = std::fs::File::open(path)
            .map_err(|source| BgzfError::Open { path: path.to_path_buf(), source })?;
        Ok(Self {
            inner: BufReader::with_capacity(128 * 1024, file),
            buf: Vec::with_capacity(MAX_BLOCK_SIZE),
            buf_pos: 0,
            block_offset: 0,
            eof: false,
            compressed_buf: Vec::with_capacity(MAX_BLOCK_SIZE),
            decompressor: libdeflater::Decompressor::new(),
        })
    }

    /// Current virtual offset (block start + position within decompressed data).
    pub fn virtual_offset(&self) -> VirtualOffset {
        VirtualOffset::new(self.block_offset, self.buf_pos as u16)
    }

    // r[impl bgzf.seek]
    #[instrument(level = "debug", skip(self), fields(voff = %voff.0), err)]
    pub fn seek_virtual(&mut self, voff: VirtualOffset) -> Result<(), BgzfError> {
        let block_off = voff.block_offset();
        let within = voff.within_block() as usize;

        self.inner.seek(SeekFrom::Start(block_off)).map_err(|_| BgzfError::SeekFailed)?;
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

    // r[impl bgzf.magic]
    // r[impl bgzf.bsize]
    // r[impl bgzf.decompression]
    // r[impl bgzf.eof]
    // r[impl bgzf.crc32]
    // r[impl bgzf.block_offset_tracking]
    // r[impl bgzf.resize_uninit]
    // r[impl bgzf.fast_header]
    fn read_block(&mut self) -> Result<bool, BgzfError> {
        // Read the full 18-byte BGZF header (gzip header + BC subfield).
        // Standard BAM files always have XLEN=6 with BC as the only subfield,
        // so we read 18 bytes and fast-path parse if the layout matches.
        let mut header = [0u8; BGZF_HEADER_SIZE];
        match self.inner.read_exact(&mut header) {
            Ok(()) => {}
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => {
                self.eof = true;
                self.buf.clear();
                self.buf_pos = 0;
                return Ok(false);
            }
            Err(_) => {
                return Err(BgzfError::TruncatedBlock);
            }
        }

        if header[..4] != BGZF_MAGIC {
            return Err(BgzfError::InvalidMagic);
        }

        let xlen = u16::from_le_bytes([header[10], header[11]]) as usize;

        // Fast path: standard BGZF header has XLEN=6 with BC subfield at bytes 12-17.
        // This covers all standard BAM files and avoids extra_buf allocation entirely.
        let (bsize, _extra_bytes_remaining) = if xlen == 6
            && header[12] == b'B'
            && header[13] == b'C'
            && header[14] == 2
            && header[15] == 0
        {
            let bsize = u16::from_le_bytes([header[16], header[17]]);
            (bsize, 0)
        } else {
            // Slow path: non-standard header layout, need to read remaining extra fields
            // We already consumed 6 bytes of extra data (bytes 12-17 in header).
            // If xlen > 6, read the rest. Then search all extra data for BC.
            let already_read = 6usize.min(xlen);
            let remaining = xlen - already_read;

            let mut extra_data = Vec::with_capacity(xlen);
            extra_data.extend_from_slice(header.get(12..12 + already_read).unwrap_or(&[]));

            if remaining > 0 {
                let start = extra_data.len();
                extra_data.resize(start + remaining, 0);
                #[allow(
                    clippy::indexing_slicing,
                    reason = "start = pre-resize len, within bounds after resize"
                )]
                self.inner
                    .read_exact(&mut extra_data[start..])
                    .map_err(|_| BgzfError::TruncatedBlock)?;
            }

            let bsize = find_bsize(&extra_data).ok_or(BgzfError::MissingBsize)?;
            (bsize, remaining)
        };

        // Total block size = BSIZE + 1; remaining = total - header - extra
        let total_block_size = bsize as usize + 1;
        let remaining_data = total_block_size
            .checked_sub(12 + xlen) // gzip header (12) + extra fields (xlen)
            .ok_or(BgzfError::BlockSizeTooSmall { bsize })?;

        // Safety: all bytes will be written by read_exact before any read.
        unsafe { resize_uninit(&mut self.compressed_buf, remaining_data) };
        self.inner.read_exact(&mut self.compressed_buf).map_err(|_| BgzfError::TruncatedBlock)?;

        if self.compressed_buf.len() < BGZF_FOOTER_SIZE {
            return Err(BgzfError::TruncatedBlock);
        }

        // Extract CRC32 and ISIZE from footer
        let footer_start = self.compressed_buf.len() - BGZF_FOOTER_SIZE;
        let crc32_bytes: [u8; 4] = self
            .compressed_buf
            .get(footer_start..footer_start + 4)
            .and_then(|s| s.try_into().ok())
            .ok_or(BgzfError::TruncatedBlock)?;
        let expected_crc = u32::from_le_bytes(crc32_bytes);

        let isize_bytes: [u8; 4] = self
            .compressed_buf
            .get(footer_start + 4..footer_start + 8)
            .and_then(|s| s.try_into().ok())
            .ok_or(BgzfError::TruncatedBlock)?;
        let uncompressed_size = u32::from_le_bytes(isize_bytes) as usize;

        // r[impl bgzf.max_block_size]
        if uncompressed_size > MAX_BLOCK_SIZE {
            return Err(BgzfError::UncompressedSizeTooLarge { isize_value: uncompressed_size });
        }

        // EOF marker: zero-length uncompressed data
        if uncompressed_size == 0 {
            self.eof = true;
            self.buf.clear();
            self.buf_pos = 0;
            return Ok(false);
        }

        let deflate_data =
            self.compressed_buf.get(..footer_start).ok_or(BgzfError::TruncatedBlock)?;

        // Safety: all bytes will be written by deflate_decompress before any read.
        unsafe { resize_uninit(&mut self.buf, uncompressed_size) };
        let actual = self
            .decompressor
            .deflate_decompress(deflate_data, &mut self.buf)
            .map_err(|source| BgzfError::DecompressionFailed { source })?;

        self.buf.truncate(actual);

        // Verify CRC32 of decompressed data
        let mut crc = libdeflater::Crc::new();
        crc.update(&self.buf);
        if crc.sum() != expected_crc {
            return Err(BgzfError::ChecksumMismatch { expected: expected_crc, found: crc.sum() });
        }

        // block_offset was set by the caller (advance_block_offset or seek_virtual)
        // before read_block was called; it stays pointing at this block's start.
        self.buf_pos = 0;

        Ok(true)
    }

    #[inline]
    pub fn read_byte(&mut self) -> Result<u8, BgzfError> {
        if self.buf_pos >= self.buf.len() {
            self.advance_block_offset();
            if !self.read_block()? {
                return Err(BgzfError::UnexpectedEof);
            }
        }
        let b = self.buf.get(self.buf_pos).copied().ok_or(BgzfError::TruncatedBlock)?;
        self.buf_pos += 1;
        Ok(b)
    }

    // r[impl bgzf.read_exact]
    #[inline]
    pub fn read_exact_into(&mut self, out: &mut [u8]) -> Result<(), BgzfError> {
        let mut written = 0;
        while written < out.len() {
            if self.buf_pos >= self.buf.len() {
                // About to read next block — advance block_offset to where we are now
                self.advance_block_offset();
                if !self.read_block()? {
                    return Err(BgzfError::UnexpectedEof);
                }
            }
            let avail = self.buf.len() - self.buf_pos;
            let need = out.len() - written;
            let n = avail.min(need);
            // n ≤ avail = buf.len() - buf_pos and n ≤ need = out.len() - written
            let dst = out.get_mut(written..written + n).ok_or(BgzfError::TruncatedBlock)?;
            let src =
                self.buf.get(self.buf_pos..self.buf_pos + n).ok_or(BgzfError::TruncatedBlock)?;
            dst.copy_from_slice(src);
            self.buf_pos += n;
            written += n;
        }
        Ok(())
    }

    // r[impl bgzf.read_partial]
    pub fn read_up_to(&mut self, out: &mut [u8]) -> Result<usize, BgzfError> {
        if self.eof && self.buf_pos >= self.buf.len() {
            return Ok(0);
        }
        if self.buf_pos >= self.buf.len() {
            self.advance_block_offset();
            if !self.read_block()? {
                return Ok(0);
            }
        }
        let avail = self.buf.len() - self.buf_pos;
        let n = avail.min(out.len());
        // n ≤ avail = buf.len() - buf_pos and n ≤ out.len()
        let dst = out.get_mut(..n).ok_or(BgzfError::TruncatedBlock)?;
        let src = self.buf.get(self.buf_pos..self.buf_pos + n).ok_or(BgzfError::TruncatedBlock)?;
        dst.copy_from_slice(src);
        self.buf_pos += n;
        Ok(n)
    }

    /// Read a little-endian i32.
    pub fn read_i32(&mut self) -> Result<i32, BgzfError> {
        let mut buf = [0u8; 4];
        self.read_exact_into(&mut buf)?;
        Ok(i32::from_le_bytes(buf))
    }

    /// Read a little-endian u32.
    pub fn read_u32(&mut self) -> Result<u32, BgzfError> {
        let mut buf = [0u8; 4];
        self.read_exact_into(&mut buf)?;
        Ok(u32::from_le_bytes(buf))
    }

    /// Advance block_offset to the stream position where the next block starts.
    /// Called before reading a new block in sequential mode.
    fn advance_block_offset(&mut self) {
        // We know block_offset points to the start of the current block.
        // We can't know the compressed block size without the header,
        // so we query stream_position only when transitioning between blocks.
        if let Ok(pos) = self.inner.stream_position() {
            self.block_offset = pos;
        }
    }
}

/// Find the BSIZE subfield (ID=BC) in the BGZF extra data.
pub(crate) fn find_bsize(extra: &[u8]) -> Option<u16> {
    let mut pos = 0;
    while pos + 4 <= extra.len() {
        let &[si1, si2, slen_lo, slen_hi, ..] = extra.get(pos..)? else { break };
        let slen = u16::from_le_bytes([slen_lo, slen_hi]) as usize;
        if si1 == b'B' && si2 == b'C' && slen == 2 {
            let bsize_bytes: [u8; 2] = extra.get(pos + 4..pos + 6)?.try_into().ok()?;
            return Some(u16::from_le_bytes(bsize_bytes));
        }
        pos += 4 + slen;
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn virtual_offset_roundtrip() {
        let vo = VirtualOffset::new(12345, 678);
        assert_eq!(vo.block_offset(), 12345);
        assert_eq!(vo.within_block(), 678);
    }

    #[test]
    fn virtual_offset_ordering() {
        let a = VirtualOffset::new(100, 50);
        let b = VirtualOffset::new(100, 60);
        let c = VirtualOffset::new(200, 0);
        assert!(a < b);
        assert!(b < c);
    }

    #[test]
    fn find_bsize_valid() {
        // BC subfield with SLEN=2, BSIZE=1234
        let mut extra = vec![b'B', b'C', 2, 0];
        extra.extend_from_slice(&1234u16.to_le_bytes());
        assert_eq!(find_bsize(&extra), Some(1234));
    }

    #[test]
    fn find_bsize_with_other_subfields() {
        // XX subfield first (SLEN=3, 3 bytes data), then BC
        let mut extra = vec![b'X', b'X', 3, 0, 0, 0, 0];
        extra.extend_from_slice(&[b'B', b'C', 2, 0]);
        extra.extend_from_slice(&999u16.to_le_bytes());
        assert_eq!(find_bsize(&extra), Some(999));
    }

    #[test]
    fn find_bsize_missing() {
        let extra = vec![b'X', b'X', 2, 0, 0, 0];
        assert_eq!(find_bsize(&extra), None);
    }

    // r[verify bgzf.max_block_size]
    #[test]
    fn rejects_uncompressed_size_exceeding_max_block_size() {
        use std::io::Write;

        // Build a BGZF block with a valid header but ISIZE claiming > 65536.
        let data = b"hello";
        let mut compressor =
            libdeflater::Compressor::new(libdeflater::CompressionLvl::new(1).unwrap());
        let bound = compressor.deflate_compress_bound(data.len());
        let mut compressed = vec![0u8; bound];
        let compressed_len = compressor.deflate_compress(data, &mut compressed).unwrap();
        compressed.truncate(compressed_len);

        let mut crc = libdeflater::Crc::new();
        crc.update(data);

        let bsize = (18 + compressed_len + 8 - 1) as u16;
        let mut block = Vec::new();
        block.extend_from_slice(&[0x1f, 0x8b, 0x08, 0x04]);
        block.extend_from_slice(&[0; 4]); // MTIME
        block.push(0); // XFL
        block.push(0xff); // OS
        block.extend_from_slice(&6u16.to_le_bytes()); // XLEN = 6
        block.extend_from_slice(&[b'B', b'C', 2, 0]); // BC subfield
        block.extend_from_slice(&bsize.to_le_bytes()); // BSIZE
        block.extend_from_slice(&compressed);
        block.extend_from_slice(&crc.sum().to_le_bytes());
        // Tamper: write ISIZE = 65537 (> MAX_BLOCK_SIZE)
        block.extend_from_slice(&65537u32.to_le_bytes());

        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("bad_isize.bgzf");
        let mut f = std::fs::File::create(&path).unwrap();
        f.write_all(&block).unwrap();
        drop(f);

        let mut reader = BgzfReader::open(&path).unwrap();
        let result = reader.read_block();
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(
            matches!(err, BgzfError::UncompressedSizeTooLarge { isize_value: 65537 }),
            "expected UncompressedSizeTooLarge, got {err:?}"
        );
    }

    #[test]
    fn fast_path_header_standard_layout() {
        // Standard BGZF header: XLEN=6, BC at fixed offset
        let header: [u8; BGZF_HEADER_SIZE] = [
            0x1f, 0x8b, 0x08, 0x04, // magic
            0x00, 0x00, 0x00, 0x00, // MTIME
            0x00, 0xff, // XFL, OS
            0x06, 0x00, // XLEN = 6
            b'B', b'C', 0x02, 0x00, // BC subfield, SLEN=2
            0x1b, 0x00, // BSIZE = 27
        ];
        // Verify the fast path conditions
        let xlen = u16::from_le_bytes([header[10], header[11]]) as usize;
        assert_eq!(xlen, 6);
        assert_eq!(header[12], b'B');
        assert_eq!(header[13], b'C');
        assert_eq!(header[14], 2);
        assert_eq!(header[15], 0);
        let bsize = u16::from_le_bytes([header[16], header[17]]);
        assert_eq!(bsize, 27);
    }
}
