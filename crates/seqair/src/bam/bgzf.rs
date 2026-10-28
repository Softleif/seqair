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
#[non_exhaustive]
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

    #[error("BGZF header arithmetic overflow (corrupt data)")]
    CorruptHeader,

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

    #[error("BAM record block_size ({block_size}) exceeds maximum (2 MiB)")]
    RecordTooLarge { block_size: usize },

    #[error(
        "region requires {total_bytes} bytes of compressed data, exceeding the {max_bytes} byte limit"
    )]
    RegionTooLarge { total_bytes: usize, max_bytes: usize },

    // r[impl bgzf.writer]
    #[error("BGZF compression failed")]
    CompressionFailed { source: libdeflater::CompressionError },

    #[error("BGZF write failed")]
    WriteFailed { source: std::io::Error },

    // r[impl bgzf.writer.finish]
    #[error("BgzfWriter already finished")]
    AlreadyFinished,

    // r[impl bgzf.writer.parallel]
    #[error("could not start a BGZF compression thread")]
    ThreadSpawn { source: std::io::Error },

    // r[impl bgzf.writer.parallel]
    #[error("a BGZF compression worker stopped before finishing its block")]
    CompressionWorkerLost,

    // Last, so no earlier variant's discriminant moves.
    // r[impl bgzf.writer.block_size]
    #[error("BGZF block of {size} bytes does not fit BGZF's 65536-byte limit")]
    BlockTooLarge { size: usize },

    // r[impl bgzf.read_errors]
    #[error("reading BGZF input failed")]
    ReadFailed { source: std::io::Error },

    // r[impl bgzf.decompression+2]
    #[error("BGZF block inflated to {found} bytes, but its ISIZE footer says {expected}")]
    SizeMismatch { expected: u32, found: usize },
}

// r[impl bgzf.libdeflate]
/// BGZF block reader that decompresses gzip blocks from a BAM file.
///
/// Tracks virtual file offsets for index-based seeking. Uses `libdeflater`
/// for decompression.
pub struct BgzfReader<R: Read + Seek> {
    inner: BufReader<R>,
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

impl<R: Read + Seek> std::fmt::Debug for BgzfReader<R> {
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

impl BgzfReader<std::fs::File> {
    #[instrument(level = "debug", fields(path = %path.display()))]
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
}

impl<R: Read + Seek> BgzfReader<R> {
    /// Create a BGZF reader from any `Read + Seek` source.
    pub fn from_reader(inner: R) -> Self {
        Self {
            inner: BufReader::with_capacity(128 * 1024, inner),
            buf: Vec::with_capacity(MAX_BLOCK_SIZE),
            buf_pos: 0,
            block_offset: 0,
            eof: false,
            compressed_buf: Vec::with_capacity(MAX_BLOCK_SIZE),
            decompressor: libdeflater::Decompressor::new(),
        }
    }

    /// Read all remaining data into a Vec.
    pub fn read_to_end(&mut self, out: &mut Vec<u8>) -> Result<(), BgzfError> {
        let mut tmp = [0u8; 8192];
        loop {
            let n = self.read_up_to(&mut tmp)?;
            if n == 0 {
                break;
            }
            out.extend_from_slice(tmp.get(..n).ok_or(BgzfError::TruncatedBlock)?);
        }
        Ok(())
    }
}

#[cfg(feature = "fuzz")]
impl BgzfReader<std::io::Cursor<Vec<u8>>> {
    pub fn from_cursor(data: Vec<u8>) -> Self {
        Self {
            inner: BufReader::with_capacity(128 * 1024, std::io::Cursor::new(data)),
            buf: Vec::with_capacity(MAX_BLOCK_SIZE),
            buf_pos: 0,
            block_offset: 0,
            eof: false,
            compressed_buf: Vec::with_capacity(MAX_BLOCK_SIZE),
            decompressor: libdeflater::Decompressor::new(),
        }
    }
}

impl<R: Read + Seek> BgzfReader<R> {
    /// Current virtual offset (block start + position within decompressed data).
    pub fn virtual_offset(&self) -> VirtualOffset {
        debug_assert!(
            self.buf_pos <= self.buf.len(),
            "buf_pos {} out of bounds for buf len {}",
            self.buf_pos,
            self.buf.len()
        );
        debug_assert!(
            self.buf_pos <= u16::MAX as usize,
            "buf_pos {} exceeds maximum representable within-block offset",
            self.buf_pos
        );
        #[allow(clippy::cast_possible_truncation, reason = "position in small buffer")]
        VirtualOffset::new(self.block_offset, self.buf_pos as u16)
    }

    // r[impl bgzf.seek]
    #[instrument(level = "debug", skip(self), fields(voff = %voff.0))]
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

    // r[impl bgzf.eof+2]
    // r[impl bgzf.empty_block]
    // r[impl bgzf.block_offset_tracking]
    // r[impl bgzf.resize_uninit]
    /// Load the next block with data into `buf`. `false` at a clean end of
    /// input; empty blocks before it are skipped.
    fn read_block(&mut self) -> Result<bool, BgzfError> {
        loop {
            self.compressed_buf.clear();
            match read_frame(&mut self.inner, &mut self.compressed_buf)? {
                Frame::EndOfInput => {
                    self.eof = true;
                    self.buf.clear();
                    self.buf_pos = 0;
                    return Ok(false);
                }
                // The next block starts where this one ended.
                Frame::Empty => self.advance_block_offset(),
                Frame::Data(block) => {
                    // Safety: `inflate_into` writes every byte or fails, and on
                    // failure `buf` is cleared before anything can read it.
                    unsafe { resize_uninit(&mut self.buf, block.isize()) };
                    let inflated = inflate_into(
                        &mut self.decompressor,
                        &self.compressed_buf,
                        block,
                        &mut self.buf,
                    );
                    self.buf_pos = 0;
                    if let Err(e) = inflated {
                        self.buf.clear();
                        return Err(e);
                    }
                    // block_offset was set by the caller (advance_block_offset or
                    // seek_virtual) before read_block; it names this block's start.
                    return Ok(true);
                }
            }
        }
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
        self.buf_pos = self.buf_pos.checked_add(1).ok_or(BgzfError::TruncatedBlock)?;
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
            let avail =
                self.buf.len().checked_sub(self.buf_pos).ok_or(BgzfError::TruncatedBlock)?;
            let need = out.len().checked_sub(written).ok_or(BgzfError::TruncatedBlock)?;
            let n = avail.min(need);
            // n ≤ avail = buf.len() - buf_pos and n ≤ need = out.len() - written
            let dst = out
                .get_mut(written..written.checked_add(n).ok_or(BgzfError::TruncatedBlock)?)
                .ok_or(BgzfError::TruncatedBlock)?;
            let src = self
                .buf
                .get(self.buf_pos..self.buf_pos.checked_add(n).ok_or(BgzfError::TruncatedBlock)?)
                .ok_or(BgzfError::TruncatedBlock)?;
            dst.copy_from_slice(src);
            self.buf_pos = self.buf_pos.checked_add(n).ok_or(BgzfError::TruncatedBlock)?;
            written = written.checked_add(n).ok_or(BgzfError::TruncatedBlock)?;
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
        let avail = self.buf.len().checked_sub(self.buf_pos).ok_or(BgzfError::TruncatedBlock)?;
        let n = avail.min(out.len());
        // n ≤ avail = buf.len() - buf_pos and n ≤ out.len()
        let dst = out.get_mut(..n).ok_or(BgzfError::TruncatedBlock)?;
        let src = self
            .buf
            .get(self.buf_pos..self.buf_pos.checked_add(n).ok_or(BgzfError::TruncatedBlock)?)
            .ok_or(BgzfError::TruncatedBlock)?;
        dst.copy_from_slice(src);
        self.buf_pos = self.buf_pos.checked_add(n).ok_or(BgzfError::TruncatedBlock)?;
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

    /// Advance `block_offset` to the stream position where the next block starts.
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

/// Decode a single BGZF block from raw bytes (header + compressed data + footer).
///
/// Returns the decompressed data. The input must start at the gzip magic bytes.
#[cfg(feature = "fuzz")]
pub fn decode_bgzf_block(data: &[u8]) -> Result<Vec<u8>, BgzfError> {
    let mut input = data;
    let mut payload = Vec::new();
    match read_frame(&mut input, &mut payload)? {
        Frame::EndOfInput => Err(BgzfError::TruncatedBlock),
        Frame::Empty => Ok(Vec::new()),
        Frame::Data(block) => {
            let mut out = vec![0u8; block.isize()];
            inflate_into(&mut libdeflater::Decompressor::new(), &payload, block, &mut out)?;
            Ok(out)
        }
    }
}

/// What [`read_frame`] found at the reader's position.
#[derive(Debug)]
pub(crate) enum Frame {
    /// Zero bytes before the next header: the input ended at a block boundary.
    EndOfInput,
    /// A block with `ISIZE = 0` (an EOF marker, or an empty block mid-stream).
    Empty,
    /// A block with data; its DEFLATE payload was appended to the caller's buffer.
    Data(DataBlock),
}

/// Where one block's DEFLATE payload sits in a buffer, and what its footer
/// says it inflates to.
#[derive(Debug, Clone, Copy)]
pub(crate) struct DataBlock {
    /// Offset of the payload in the buffer [`read_frame`] appended to.
    start: usize,
    /// Payload length in bytes (footer excluded).
    len: usize,
    crc: u32,
    /// `ISIZE`, in `1..=MAX_BLOCK_SIZE`.
    isize: u32,
}

impl DataBlock {
    /// The uncompressed size the footer promises, in `1..=65536`.
    pub(crate) fn isize(self) -> usize {
        self.isize as usize
    }
}

/// Fill `buf` from `inner`, retrying `Interrupted`; the number of bytes read,
/// short only at the end of input.
fn read_full(inner: &mut impl Read, buf: &mut [u8]) -> Result<usize, BgzfError> {
    let mut filled = 0;
    while let Some(rest) = buf.get_mut(filled..)
        && !rest.is_empty()
    {
        match inner.read(rest) {
            Ok(0) => break,
            Ok(n) => filled = filled.saturating_add(n),
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(source) => return Err(BgzfError::ReadFailed { source }),
        }
    }
    Ok(filled)
}

/// `read_exact`, with an early end of input reported as truncation.
fn read_exact_or_truncated(inner: &mut impl Read, buf: &mut [u8]) -> Result<(), BgzfError> {
    if read_full(inner, buf)? == buf.len() { Ok(()) } else { Err(BgzfError::TruncatedBlock) }
}

// r[impl bgzf.magic]
// r[impl bgzf.bsize]
// r[impl bgzf.fast_header]
// r[impl bgzf.max_block_size]
// r[impl bgzf.truncated_header]
// r[impl bgzf.read_errors]
/// Read one block's framing from `inner`, appending its DEFLATE payload (for a
/// block with data) to `payload`.
///
/// The block parser of [`BgzfReader`]: an end of input at a block boundary is [`Frame::EndOfInput`]; anywhere else it is
/// [`BgzfError::TruncatedBlock`]. On error `payload` may hold a partial
/// payload past its original length.
pub(crate) fn read_frame(inner: &mut impl Read, payload: &mut Vec<u8>) -> Result<Frame, BgzfError> {
    let mut header = [0u8; BGZF_HEADER_SIZE];
    match read_full(inner, &mut header)? {
        0 => return Ok(Frame::EndOfInput),
        BGZF_HEADER_SIZE => {}
        _ => return Err(BgzfError::TruncatedBlock),
    }
    let [m0, m1, m2, m3, _, _, _, _, _, _, x0, x1, si1, si2, sl0, sl1, b0, b1] = header;
    if [m0, m1, m2, m3] != BGZF_MAGIC {
        return Err(BgzfError::InvalidMagic);
    }
    let xlen = usize::from(u16::from_le_bytes([x0, x1]));

    // Every BAM writer we know of puts BC alone in a 6-byte extra field, which
    // the 18-byte header already holds.
    let bsize = if xlen == 6 && [si1, si2, sl0, sl1] == [b'B', b'C', 2, 0] {
        u16::from_le_bytes([b0, b1])
    } else {
        let mut extra = vec![0u8; xlen];
        let held = xlen.min(6);
        let (head, tail) = extra.split_at_mut(held);
        head.copy_from_slice(header.get(12..12usize.saturating_add(held)).unwrap_or_default());
        read_exact_or_truncated(inner, tail)?;
        find_bsize(&extra).ok_or(BgzfError::MissingBsize)?
    };

    // BSIZE + 1 is the whole block; 12 + XLEN of it is the header.
    let rest = (usize::from(bsize).saturating_add(1))
        .checked_sub(12usize.saturating_add(xlen))
        .ok_or(BgzfError::BlockSizeTooSmall { bsize })?;
    let len = rest.checked_sub(BGZF_FOOTER_SIZE).ok_or(BgzfError::TruncatedBlock)?;

    let start = payload.len();
    payload.resize(start.saturating_add(len), 0);
    read_exact_or_truncated(inner, payload.get_mut(start..).unwrap_or_default())?;
    let mut footer = [0u8; BGZF_FOOTER_SIZE];
    read_exact_or_truncated(inner, &mut footer)?;
    let [c0, c1, c2, c3, i0, i1, i2, i3] = footer;
    let crc = u32::from_le_bytes([c0, c1, c2, c3]);
    let isize = u32::from_le_bytes([i0, i1, i2, i3]);

    if isize as usize > MAX_BLOCK_SIZE {
        payload.truncate(start);
        return Err(BgzfError::UncompressedSizeTooLarge { isize_value: isize as usize });
    }
    if isize == 0 {
        payload.truncate(start);
        return Ok(Frame::Empty);
    }
    Ok(Frame::Data(DataBlock { start, len, crc, isize }))
}

// r[impl bgzf.decompression+2]
// r[impl bgzf.crc32]
// r[impl bgzf.libdeflate]
/// Inflate `block`, whose payload is in `payloads` where [`read_frame`] put
/// it, into `out` — exactly `block.isize()` bytes — and check its CRC32.
pub(crate) fn inflate_into(
    decompressor: &mut libdeflater::Decompressor,
    payloads: &[u8],
    block: DataBlock,
    out: &mut [u8],
) -> Result<(), BgzfError> {
    debug_assert_eq!(out.len(), block.isize(), "inflate_into needs an ISIZE-sized buffer");
    let deflate = payloads
        .get(block.start..block.start.saturating_add(block.len))
        .ok_or(BgzfError::TruncatedBlock)?;
    let actual = decompressor
        .deflate_decompress(deflate, out)
        .map_err(|source| BgzfError::DecompressionFailed { source })?;
    if actual != out.len() {
        return Err(BgzfError::SizeMismatch { expected: block.isize, found: actual });
    }
    let mut crc = libdeflater::Crc::new();
    crc.update(out);
    if crc.sum() != block.crc {
        return Err(BgzfError::ChecksumMismatch { expected: block.crc, found: crc.sum() });
    }
    Ok(())
}

/// Find the BSIZE subfield (ID=BC) in the BGZF extra data.
pub(crate) fn find_bsize(extra: &[u8]) -> Option<u16> {
    let mut pos: usize = 0;
    while pos.saturating_add(4) <= extra.len() {
        let &[si1, si2, slen_lo, slen_hi, ..] = extra.get(pos..)? else { break };
        let slen = u16::from_le_bytes([slen_lo, slen_hi]) as usize;
        if si1 == b'B' && si2 == b'C' && slen == 2 {
            let bsize_bytes: [u8; 2] =
                extra.get(pos.checked_add(4)?..pos.checked_add(6)?)?.try_into().ok()?;
            return Some(u16::from_le_bytes(bsize_bytes));
        }
        pos = pos.checked_add(4)?.checked_add(slen)?;
    }
    None
}

#[cfg(test)]
#[allow(clippy::arithmetic_side_effects, reason = "tests")]
#[allow(clippy::cast_possible_truncation, reason = "tests")]
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
