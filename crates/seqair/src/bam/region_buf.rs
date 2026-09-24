//! Streaming windowed BGZF buffer for high-latency I/O (cluster/NFS storage).
//!
//! Reads compressed bytes for a set of BAM index chunks in a bounded sliding
//! window (a few large sequential reads), refilling forward on demand as
//! decompression advances. Peak resident memory is the window budget, not the
//! whole region — see `WINDOW_BUDGET`.

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

/// Padding added past each chunk's end block offset when computing the byte
/// range to load. Covers the final BGZF block plus a little extra for the
/// record that straddles the end boundary.
pub(super) const CHUNK_END_PAD: usize = MAX_BLOCK_SIZE;

/// Target resident compressed bytes per [`RegionBuf`] window.
///
/// Each refill reads up to this much in one large sequential read, so memory
/// stays bounded regardless of how large the queried region is. The effective
/// budget is floored at one max BGZF block ([`MAX_BLOCK_SIZE`]) so any single
/// block always fits.
pub(super) const WINDOW_BUDGET: usize = 64 * 1024 * 1024; // 64 MiB

/// How many decompressed BGZF blocks a [`BlockCache`] keeps.
///
/// A query starts at the linear-index minimum of its first 16 kb window, so
/// neighbouring small queries re-read the same blocks: at ~30× short-read
/// coverage a 16 kb window is ~24 blocks, and a read with a long deletion or
/// skip can pull a query's start back further. Measured on 1 bp queries every
/// 100 bp over 10 Mb of 30× WGS, 32 blocks still re-inflated 1.4× the blocks
/// touched (the LRU thrashes on a query longer than the cache); 64 re-inflated
/// none. That is ≤ 4 MiB per reader, allocated only as blocks arrive.
pub(crate) const BLOCK_CACHE_BLOCKS: usize = 64;

// r[impl region_buf.block_cache]
/// Decompressed BGZF blocks kept across the queries of one reader, keyed by
/// the block's compressed file offset.
///
/// A [`RegionBuf`] borrows the cache for one query. Each block it decompresses
/// is handed back when it moves past the block, and a block that is already
/// here is taken instead of decompressed again. Buffers move between the
/// cache and the `RegionBuf` by swapping `Vec`s, so a cold scan pays no copy
/// for it, only the bookkeeping.
///
/// Only blocks that decompressed and passed their CRC check get in, so a hit
/// skips both. Forks start with an empty cache of their own.
#[derive(Debug)]
pub(crate) struct BlockCache {
    entries: Vec<CachedBlock>,
    /// Most entries kept; [`BLOCK_CACHE_BLOCKS`] outside tests.
    capacity: usize,
    /// Recency clock: bumped on every `put`, stamped on the entry.
    clock: u64,
    /// An empty buffer to hand out while the cache is still filling.
    spare: Vec<u8>,
}

#[derive(Debug)]
struct CachedBlock {
    /// File offset of the block's first compressed byte.
    offset: u64,
    /// Total compressed length (`BSIZE + 1`), checked on lookup.
    block_len: usize,
    data: Vec<u8>,
    last_used: u64,
}

impl BlockCache {
    pub(crate) fn new() -> Self {
        Self::with_capacity(BLOCK_CACHE_BLOCKS)
    }

    fn with_capacity(capacity: usize) -> Self {
        BlockCache { entries: Vec::new(), capacity, clock: 0, spare: Vec::new() }
    }

    /// Compressed length of the cached block at `offset`.
    fn block_len(&self, offset: u64) -> Option<usize> {
        self.entries.iter().find(|e| e.offset == offset).map(|e| e.block_len)
    }

    /// Remove and return the decompressed block at `offset`, if cached.
    fn take(&mut self, offset: u64, block_len: usize) -> Option<Vec<u8>> {
        let i = self.entries.iter().position(|e| e.offset == offset)?;
        let entry = self.entries.swap_remove(i);
        (entry.block_len == block_len).then_some(entry.data)
    }

    /// Store `data` as the block at `offset`, returning an empty buffer to
    /// decompress the next block into (the evicted entry's, when full).
    fn put(&mut self, offset: u64, block_len: usize, data: Vec<u8>) -> Vec<u8> {
        self.clock = self.clock.wrapping_add(1);
        let entry = CachedBlock { offset, block_len, data, last_used: self.clock };
        // Replace a stale copy of the same block, else the least recently
        // used entry once full.
        let slot = match self.entries.iter().position(|e| e.offset == offset) {
            Some(i) => self.entries.get_mut(i),
            None if self.entries.len() < self.capacity => None,
            None => self.entries.iter_mut().min_by_key(|e| e.last_used),
        };
        let mut recycled = match slot {
            Some(slot) => std::mem::replace(slot, entry).data,
            None => {
                self.entries.push(entry);
                std::mem::take(&mut self.spare)
            }
        };
        recycled.clear();
        recycled
    }

    /// Keep an unused buffer for the next `put` that doesn't evict.
    fn recycle(&mut self, mut buf: Vec<u8>) {
        if self.spare.capacity() < buf.capacity() {
            buf.clear();
            self.spare = buf;
        }
    }
}

/// Pre-computed byte range covering one or more merged index chunks.
struct MergedRange {
    file_start: u64,
    file_end: u64,
}

/// Streaming BGZF decompressor over a region's index chunks.
///
/// Created by [`RegionBuf::new`], which borrows a seekable reader and computes
/// the merged byte ranges to stream, but reads nothing eagerly. Compressed
/// bytes are pulled into a bounded `WINDOW_BUDGET`-sized window on demand and
/// decompressed block-by-block.
pub struct RegionBuf<'r, R: Read + Seek> {
    /// Reader kept alive to refill the window. Bulk, unbuffered.
    reader: &'r mut R,
    /// Merged byte ranges to stream, sorted by file offset (the read plan).
    ranges: Vec<MergedRange>,
    /// Index of the range currently being streamed.
    range_idx: usize,
    /// Resident compressed bytes for the current position within `ranges[range_idx]`.
    window: Vec<u8>,
    /// File offset of `window[0]`.
    ///
    /// Invariant: `window_file_start + cursor` is the true file offset of
    /// `window[cursor]`, and the window lies within `ranges[range_idx]`.
    window_file_start: u64,
    /// Read position within `window`.
    cursor: usize,
    /// Authoritative file length, captured once; refills never read past it.
    file_size: u64,
    /// Per-window compressed-byte budget (floored at `MAX_BLOCK_SIZE`).
    budget: usize,
    /// Decompressed block buffer (reused across blocks).
    buf: Vec<u8>,
    /// Current position within `buf`.
    buf_pos: usize,
    /// True file offset of the current BGZF block.
    block_offset: u64,
    eof: bool,
    decompressor: libdeflater::Decompressor,
    blocks_decompressed: u32,
    decompressed_bytes: u64,
    /// Blocks shared with the reader's other queries, when it lent one.
    cache: Option<&'r mut BlockCache>,
    /// `(file offset, compressed length)` of the verified block in `buf`,
    /// while it is one the cache should get back.
    buf_block: Option<(u64, usize)>,
}

impl<R: Read + Seek> std::fmt::Debug for RegionBuf<'_, R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("RegionBuf")
            .field("ranges", &self.ranges.len())
            .field("range_idx", &self.range_idx)
            .field("window_len", &self.window.len())
            .field("window_file_start", &self.window_file_start)
            .field("cursor", &self.cursor)
            .field("block_offset", &self.block_offset)
            .field("eof", &self.eof)
            .finish()
    }
}

impl<'r, R: Read + Seek> RegionBuf<'r, R> {
    // r[impl region_buf.new]
    // r[impl region_buf.empty]
    // r[related region_buf.merge_chunks]
    /// Plan a streaming read over the given chunks.
    ///
    /// Merges overlapping/adjacent chunks into the ranges to stream, but reads
    /// nothing — the first window is pulled lazily on the first decode (or
    /// [`seek_virtual`](Self::seek_virtual) with a non-zero intra-block offset).
    pub fn new(reader: &'r mut R, chunks: &[Chunk]) -> Result<Self, BgzfError> {
        Self::with_budget(reader, chunks, WINDOW_BUDGET)
    }

    // r[impl region_buf.window_budget]
    /// Like [`new`](Self::new) but with an explicit window budget. Exposed for
    /// tests that force many refills with a tiny budget; the budget is floored
    /// at one max BGZF block so a single block always fits.
    #[doc(hidden)]
    pub fn with_budget(
        reader: &'r mut R,
        chunks: &[Chunk],
        budget: usize,
    ) -> Result<Self, BgzfError> {
        Self::build(reader, chunks, budget, None)
    }

    /// Like [`new`](Self::new), keeping decompressed blocks in `cache` so the
    /// reader's next query can reuse them.
    pub(crate) fn with_cache(
        reader: &'r mut R,
        chunks: &[Chunk],
        cache: &'r mut BlockCache,
    ) -> Result<Self, BgzfError> {
        Self::build(reader, chunks, WINDOW_BUDGET, Some(cache))
    }

    /// [`with_cache`](Self::with_cache) with an explicit window budget, for tests.
    #[cfg(test)]
    pub(crate) fn with_cache_and_budget(
        reader: &'r mut R,
        chunks: &[Chunk],
        cache: &'r mut BlockCache,
        budget: usize,
    ) -> Result<Self, BgzfError> {
        Self::build(reader, chunks, budget, Some(cache))
    }

    fn build(
        reader: &'r mut R,
        chunks: &[Chunk],
        budget: usize,
        cache: Option<&'r mut BlockCache>,
    ) -> Result<Self, BgzfError> {
        let budget = budget.max(MAX_BLOCK_SIZE);
        let ranges = merge_chunks(chunks);
        let eof = ranges.is_empty();
        let window_file_start = ranges.first().map(|r| r.file_start).unwrap_or(0);
        // r[impl io.fuzz.alloc_limits]
        // Capture the real file length so refills never over-allocate when a
        // corrupt index points a chunk's end far past EOF.
        let file_size = reader.seek(SeekFrom::End(0)).map_err(|_| BgzfError::SeekFailed)?;
        Ok(RegionBuf {
            reader,
            ranges,
            range_idx: 0,
            window: Vec::new(),
            window_file_start,
            cursor: 0,
            file_size,
            budget,
            buf: Vec::with_capacity(MAX_BLOCK_SIZE),
            buf_pos: 0,
            block_offset: window_file_start,
            eof,
            decompressor: libdeflater::Decompressor::new(),
            blocks_decompressed: 0,
            decompressed_bytes: 0,
            cache,
            buf_block: None,
        })
    }

    /// Ensure at least `need` compressed bytes are resident at/after `cursor`
    /// within the current range, refilling from the reader as necessary.
    /// Returns the bytes available at/after `cursor` (may be `< need` at a
    /// range tail or EOF — the caller then advances the range or stops).
    fn ensure_available(&mut self, need: usize) -> Result<usize, BgzfError> {
        let avail = self.window.len().wrapping_sub(self.cursor);
        if avail >= need {
            return Ok(avail);
        }

        // Drop the consumed prefix, keeping `window_file_start + cursor` exact.
        if self.cursor > 0 {
            let keep = self.window.len().wrapping_sub(self.cursor);
            self.window.copy_within(self.cursor.., 0);
            self.window.truncate(keep);
            self.window_file_start = self.window_file_start.wrapping_add(self.cursor as u64);
            self.cursor = 0;
        }

        let range_end = match self.ranges.get(self.range_idx) {
            Some(range) => range.file_end.min(self.file_size),
            None => return Ok(self.window.len()),
        };

        while self.window.len() < need {
            let next_file = self.window_file_start.wrapping_add(self.window.len() as u64);
            let range_remaining = range_end.saturating_sub(next_file);
            if range_remaining == 0 {
                break;
            }
            // Grow to at least `need`, but prefer a full budget-sized read so
            // refills stay rare. Never read past the range (clamped to EOF).
            let want = need
                .saturating_sub(self.window.len())
                .max(self.budget.saturating_sub(self.window.len()));
            #[expect(
                clippy::cast_possible_truncation,
                reason = "min(want, range_remaining) ≤ want ≤ budget+block, fits usize on 64-bit"
            )]
            let read_len = (want as u64).min(range_remaining) as usize;
            if read_len == 0 {
                break;
            }

            self.reader.seek(SeekFrom::Start(next_file)).map_err(|_| BgzfError::SeekFailed)?;
            let old_len = self.window.len();
            // Safety: read_all overwrites [old_len..] before any read; we
            // truncate to the bytes actually read immediately after.
            unsafe { bgzf::resize_uninit(&mut self.window, old_len.wrapping_add(read_len)) };
            #[allow(clippy::indexing_slicing, reason = "old_len ≤ window.len() after resize above")]
            let got = read_all(self.reader, &mut self.window[old_len..]);
            self.window.truncate(old_len.wrapping_add(got));
            if got < read_len {
                // Reader delivered fewer bytes than the (EOF-clamped) range
                // promised — treat as the end of this range.
                break;
            }
        }

        Ok(self.window.len().wrapping_sub(self.cursor))
    }

    /// Step to the next merged range, resetting the window to its start.
    /// Returns `false` when no further range exists.
    fn advance_range(&mut self) -> bool {
        self.range_idx = self.range_idx.wrapping_add(1);
        let Some(range) = self.ranges.get(self.range_idx) else {
            return false;
        };
        self.window.clear();
        self.window_file_start = range.file_start;
        self.cursor = 0;
        true
    }

    /// Empty `buf`, handing the block it holds back to the cache first.
    fn release_block(&mut self) {
        release_block(self.cache.as_deref_mut(), &mut self.buf_block, &mut self.buf);
    }

    // r[impl region_buf.seek_virtual]
    /// Seek to a virtual offset within the planned ranges.
    ///
    /// Repositions the window if the target lies outside the resident bytes.
    /// Errors if the offset falls in a gap between ranges or outside all of
    /// them (validated against the *planned* ranges, which may extend up to
    /// `CHUNK_END_PAD` past the last record).
    pub fn seek_virtual(&mut self, voff: VirtualOffset) -> Result<(), BgzfError> {
        let block_off = voff.block_offset();
        let within = voff.within_block() as usize;

        let idx = self
            .ranges
            .iter()
            .position(|r| block_off >= r.file_start && block_off < r.file_end)
            .ok_or(BgzfError::VirtualOffsetOutOfRange { offset: block_off })?;

        let win_end = self.window_file_start.wrapping_add(self.window.len() as u64);
        if idx == self.range_idx && block_off >= self.window_file_start && block_off < win_end {
            #[expect(
                clippy::cast_possible_truncation,
                reason = "block_off - window_file_start < window.len() ≤ budget, fits usize"
            )]
            let new_cursor = block_off.wrapping_sub(self.window_file_start) as usize;
            self.cursor = new_cursor;
        } else {
            self.range_idx = idx;
            self.window.clear();
            self.window_file_start = block_off;
            self.cursor = 0;
        }

        self.block_offset = block_off;
        self.release_block();
        self.buf_pos = 0;
        self.eof = false;

        if within > 0 {
            self.read_block()?;
            self.buf_pos = within.min(self.buf.len());
        }

        Ok(())
    }

    // r[impl region_buf.virtual_offset+2]
    pub fn virtual_offset(&self) -> VirtualOffset {
        match u16::try_from(self.buf_pos) {
            Ok(within) => VirtualOffset::new(self.block_offset, within),
            // The end of a full 64 KiB block: the same position is the next
            // block's first byte, and `window_file_start + cursor` is where the
            // next block starts once this one has been read.
            Err(_) => {
                VirtualOffset::new(self.window_file_start.wrapping_add(self.cursor as u64), 0)
            }
        }
    }

    // r[impl region_buf.decompress]
    // r[impl region_buf.fast_header]
    fn read_block(&mut self) -> Result<bool, BgzfError> {
        loop {
            if self.take_cached_block() {
                return Ok(true);
            }

            // Ensure the 18-byte header is resident (refilling if needed).
            let avail = self.ensure_available(BGZF_HEADER_SIZE)?;
            if avail < BGZF_HEADER_SIZE {
                if self.advance_range() {
                    continue;
                }
                self.eof = true;
                self.release_block();
                self.buf_pos = 0;
                return Ok(false);
            }

            // The window may have compacted in `ensure_available`; `cursor` now
            // points at the block start. Record its TRUE file offset here — any
            // later compaction preserves `window_file_start + cursor` for this
            // block, so this value stays correct through the full-block ensure.
            self.block_offset = self.window_file_start.wrapping_add(self.cursor as u64);

            // Copy the header out so later `ensure_available` calls (which take
            // `&mut self`) don't conflict with a borrow of `window`.
            #[allow(clippy::indexing_slicing, reason = "avail ≥ 18 checked above")]
            let header: [u8; BGZF_HEADER_SIZE] = self.window
                [self.cursor..self.cursor.wrapping_add(BGZF_HEADER_SIZE)]
                .try_into()
                .map_err(|_| BgzfError::TruncatedBlock)?;

            if header[..4] != BGZF_MAGIC {
                return Err(BgzfError::InvalidMagic);
            }

            let xlen = u16::from_le_bytes([header[10], header[11]]) as usize;

            // Fast path: standard BGZF header carries BSIZE in the BC subfield.
            let bsize = if xlen == 6
                && header[12] == b'B'
                && header[13] == b'C'
                && header[14] == 2
                && header[15] == 0
            {
                u16::from_le_bytes([header[16], header[17]])
            } else {
                // Non-standard extra field: ensure it's resident, then scan it.
                let extra_len = 12usize.wrapping_add(xlen);
                let avail = self.ensure_available(extra_len)?;
                if avail < extra_len {
                    return Err(BgzfError::TruncatedBlock);
                }
                let extra_start = self.cursor.wrapping_add(12);
                let extra_end = self.cursor.wrapping_add(extra_len);
                #[allow(clippy::indexing_slicing, reason = "avail ≥ extra_len checked above")]
                bgzf::find_bsize(&self.window[extra_start..extra_end])
                    .ok_or(BgzfError::MissingBsize)?
            };

            let total_block_size = (bsize as usize).wrapping_add(1);

            // Ensure the whole block is resident.
            let avail = self.ensure_available(total_block_size)?;
            if avail < total_block_size {
                // Partial block at the range tail — advance, or EOF.
                if self.advance_range() {
                    continue;
                }
                self.eof = true;
                self.release_block();
                self.buf_pos = 0;
                return Ok(false);
            }

            let block_end = self.cursor.wrapping_add(total_block_size);
            let data_start = self.cursor.wrapping_add(12).wrapping_add(xlen);
            if data_start > block_end {
                return Err(BgzfError::BlockSizeTooSmall { bsize });
            }
            let remaining =
                self.window.get(data_start..block_end).ok_or(BgzfError::TruncatedBlock)?;

            if remaining.len() < BGZF_FOOTER_SIZE {
                return Err(BgzfError::TruncatedBlock);
            }

            let footer_start = remaining.len().wrapping_sub(BGZF_FOOTER_SIZE);
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
                self.release_block();
                self.buf_pos = 0;
                return Ok(false);
            }

            // The block we are leaving goes back to the cache before `buf` is
            // overwritten.
            let block_file_offset = self.block_offset;
            release_block(self.cache.as_deref_mut(), &mut self.buf_block, &mut self.buf);

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
                return Err(BgzfError::ChecksumMismatch {
                    expected: expected_crc,
                    found: crc.sum(),
                });
            }

            self.buf_pos = 0;
            if self.cache.is_some() {
                self.buf_block = Some((block_file_offset, total_block_size));
            }
            self.blocks_decompressed = self.blocks_decompressed.wrapping_add(1);
            self.decompressed_bytes = self.decompressed_bytes.wrapping_add(actual as u64);
            return Ok(true);
        }
    }

    // r[impl region_buf.block_cache]
    /// Make the cached block at the cursor current, without reading or
    /// decompressing anything. `false` when there is no cache, the block isn't
    /// in it, or the block would not lie wholly inside the current range —
    /// the uncached path then reads it, and decides as it always has.
    fn take_cached_block(&mut self) -> bool {
        let Some(cache) = self.cache.as_deref_mut() else {
            return false;
        };
        let Some(range) = self.ranges.get(self.range_idx) else {
            return false;
        };
        let offset = self.window_file_start.wrapping_add(self.cursor as u64);
        let Some(block_len) = cache.block_len(offset) else {
            return false;
        };
        let range_end = range.file_end.min(self.file_size);
        if offset.saturating_add(block_len as u64) > range_end {
            return false;
        }
        release_block(Some(&mut *cache), &mut self.buf_block, &mut self.buf);
        let Some(data) = cache.take(offset, block_len) else {
            return false;
        };
        cache.recycle(std::mem::replace(&mut self.buf, data));
        self.buf_block = Some((offset, block_len));
        self.buf_pos = 0;
        self.block_offset = offset;

        // Step past the block's compressed bytes, which may not be resident:
        // `window_file_start + cursor` stays the true file offset either way.
        let next = self.cursor.saturating_add(block_len);
        if next <= self.window.len() {
            self.cursor = next;
        } else {
            self.window.clear();
            self.window_file_start = offset.wrapping_add(block_len as u64);
            self.cursor = 0;
        }
        true
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
            // Bounds are checked here rather than with `get().ok_or(..)`: an
            // eagerly built error is dropped on every call that succeeds.
            let (Some(dst), Some(src)) = (
                out.get_mut(written..written.wrapping_add(n)),
                self.buf.get(self.buf_pos..self.buf_pos.wrapping_add(n)),
            ) else {
                return Err(BgzfError::TruncatedBlock);
            };
            dst.copy_from_slice(src);
            self.buf_pos = self.buf_pos.wrapping_add(n);
            written = written.wrapping_add(n);
        }
        Ok(())
    }

    /// The unread rest of the current decompressed block and the virtual
    /// offset of its first byte, loading the next block first when this one
    /// is exhausted. Empty only once every planned range is exhausted. Pair
    /// with [`Self::consume`].
    #[inline]
    pub fn fill_buf(&mut self) -> Result<(VirtualOffset, &[u8]), BgzfError> {
        if self.buf_pos >= self.buf.len() && !self.read_block()? {
            return Ok((self.virtual_offset(), &[]));
        }
        Ok((self.virtual_offset(), self.buf.get(self.buf_pos..).unwrap_or_default()))
    }

    /// Mark `n` bytes of what [`Self::fill_buf`] returned as read.
    #[inline]
    pub fn consume(&mut self, n: usize) {
        self.buf_pos = self.buf_pos.saturating_add(n).min(self.buf.len());
    }

    #[inline]
    pub fn read_byte(&mut self) -> Result<u8, BgzfError> {
        if self.buf_pos >= self.buf.len() && !self.read_block()? {
            return Err(BgzfError::UnexpectedEof);
        }
        let Some(&b) = self.buf.get(self.buf_pos) else {
            return Err(BgzfError::TruncatedBlock);
        };
        self.buf_pos = self.buf_pos.wrapping_add(1);
        Ok(b)
    }

    pub fn read_u32(&mut self) -> Result<u32, BgzfError> {
        let mut buf = [0u8; 4];
        self.read_exact_into(&mut buf)?;
        Ok(u32::from_le_bytes(buf))
    }

    /// Read a complete BAM record body (`block_size` bytes, not including the
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
    ///
    /// # Returns
    ///
    /// `Ok(None)` when every planned byte range is exhausted *at a record
    /// boundary* — there is no next record and there never will be. That is a
    /// different fact from `Err(UnexpectedEof)`, which means the data ran out
    /// partway through a record and so the file is truncated. Returning both as
    /// the same error let a caller mistake the first for "refill and retry" and
    /// spin forever; see `r[region_buf.record_boundary_eof]`.
    // r[impl region_buf.record_boundary_eof]
    pub fn read_record<'a>(
        &'a mut self,
        scratch: &'a mut Vec<u8>,
    ) -> Result<Option<&'a [u8]>, BgzfError> {
        // Ensure the decompressed buffer has data to read the 4-byte length from.
        if self.buf_pos >= self.buf.len() && !self.read_block()? {
            return Ok(None);
        }

        // Fast-path u32 read: all 4 bytes in the current block. No eagerly
        // built `BgzfError` on this path: dropping the unused one on every
        // record showed up in profiles.
        let head = self.buf.get(self.buf_pos..).and_then(<[u8]>::first_chunk::<4>).copied();
        let block_size = if let Some(bytes) = head {
            self.buf_pos = self.buf_pos.wrapping_add(4);
            u32::from_le_bytes(bytes) as usize
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
        let body_end = self.buf_pos.wrapping_add(block_size);
        if body_end <= self.buf.len() {
            debug_assert!(
                self.buf_pos <= body_end,
                "buf_pos {} > body_end {body_end}",
                self.buf_pos
            );
            #[allow(
                clippy::indexing_slicing,
                reason = "buf_pos ≤ body_end ≤ buf.len() checked above"
            )]
            let slice = &self.buf[self.buf_pos..body_end];
            self.buf_pos = body_end;
            return Ok(Some(slice));
        }

        // Slow path: record spans a block boundary — copy into scratch.
        scratch.clear();
        // Safety: read_exact_into overwrites every byte before any are read.
        unsafe { super::bgzf::resize_uninit(scratch, block_size) };
        self.read_exact_into(scratch)?;
        Ok(Some(scratch))
    }
}

// r[impl region_buf.drop_no_panic]
impl<R: Read + Seek> Drop for RegionBuf<'_, R> {
    fn drop(&mut self) {
        self.release_block();
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
                decompressed_bytes = self.decompressed_bytes,
                ranges = self.ranges.len(),
                max_gap_bytes = max_gap,
                window_capacity = self.window.capacity(),
                buf_capacity = self.buf.capacity(),
                "region_buf summary",
            );
        }
    }
}

/// Empty `buf`, handing the verified block it holds back to `cache` first.
///
/// Free-standing so `read_block` can call it while it still borrows the window.
fn release_block(
    cache: Option<&mut BlockCache>,
    buf_block: &mut Option<(u64, usize)>,
    buf: &mut Vec<u8>,
) {
    if let Some((offset, block_len)) = buf_block.take()
        && let Some(cache) = cache
    {
        let data = std::mem::take(buf);
        *buf = cache.put(offset, block_len, data);
    } else {
        buf.clear();
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

/// Compute the total compressed bytes a [`RegionBuf`] would stream for `chunks`.
///
/// This is the merged byte span of the chunks' index ranges — what byte-aware
/// segmentation budgets against (see `estimate_region_bytes`). It reflects the
/// planned read size, independent of the smaller resident window.
#[expect(
    clippy::cast_possible_truncation,
    reason = "BAM files are < 2^63 bytes; merged span fits in usize on 64-bit"
)]
pub(crate) fn merged_byte_size(chunks: &[Chunk]) -> usize {
    merge_chunks(chunks)
        .iter()
        .map(|r| {
            debug_assert!(
                r.file_end > r.file_start,
                "merged range has non-positive size: file_start={}, file_end={}",
                r.file_start,
                r.file_end
            );
            r.file_end.saturating_sub(r.file_start) as usize
        })
        .fold(0usize, usize::saturating_add)
}

// r[impl region_buf.merge_chunks]
/// Merge chunks into non-overlapping byte ranges sorted by file offset.
fn merge_chunks(chunks: &[Chunk]) -> Vec<MergedRange> {
    let mut offsets: Vec<(u64, u64)> = chunks
        .iter()
        .filter_map(|c| {
            let start = c.begin.block_offset();
            // end's block_offset points to the block containing the last byte;
            // extend by CHUNK_END_PAD to capture the full final block plus any
            // record that may straddle a sub-chunk boundary during batching.
            let end = c.end.block_offset().saturating_add(CHUNK_END_PAD as u64);
            // Skip degenerate chunks (e.g. from corrupt index data where begin > end).
            (start < end).then_some((start, end))
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
#[allow(
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    reason = "test code with known small values"
)]
mod tests {
    use super::*;

    // r[verify region_buf.drop_no_panic]
    #[test]
    fn drop_does_not_panic_with_overlapping_ranges() {
        // Simulate a RegionBuf with overlapping ranges where file_start < prev file_end
        // (which would cause subtraction underflow without saturating_sub).
        let mut cursor = std::io::Cursor::new(vec![0u8; 100]);
        let buf = RegionBuf {
            reader: &mut cursor,
            ranges: vec![
                MergedRange { file_start: 100, file_end: 300 },
                MergedRange { file_start: 200, file_end: 400 },
            ],
            range_idx: 0,
            window: vec![0; 100],
            window_file_start: 0,
            cursor: 0,
            file_size: 100,
            budget: WINDOW_BUDGET,
            buf: Vec::new(),
            buf_pos: 0,
            block_offset: 0,
            eof: false,
            blocks_decompressed: 1,
            decompressed_bytes: 100,
            decompressor: libdeflater::Decompressor::new(),
            cache: None,
            buf_block: None,
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
        let mut buf = RegionBuf::new(&mut cursor, &chunks).unwrap();
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
            Chunk { begin: VirtualOffset::new(200_000, 0), end: VirtualOffset::new(300_000, 0) },
        ];
        let ranges = merge_chunks(&chunks);
        assert_eq!(ranges.len(), 2);
    }

    #[test]
    fn empty_chunks_is_eof() {
        let mut cursor = std::io::Cursor::new(vec![]);
        let buf = RegionBuf::new(&mut cursor, &[]).unwrap();
        assert!(buf.eof);
    }

    // --- Disjoint range regression tests ---

    /// Helper: build fake file data of `len` bytes (just sequential bytes).
    fn fake_file(len: usize) -> Vec<u8> {
        (0..len).map(|i| (i % 256) as u8).collect()
    }

    /// The true file offset of the byte at `cursor`.
    fn file_pos<R: Read + Seek>(buf: &RegionBuf<'_, R>) -> u64 {
        buf.window_file_start + buf.cursor as u64
    }

    #[test]
    fn single_range_seek_tracks_file_offset() {
        let file_data = fake_file(1000);
        let mut cursor = std::io::Cursor::new(file_data);

        let chunks =
            vec![Chunk { begin: VirtualOffset::new(100, 0), end: VirtualOffset::new(200, 0) }];

        let mut buf = RegionBuf::new(&mut cursor, &chunks).unwrap();

        buf.seek_virtual(VirtualOffset::new(100, 0)).unwrap();
        assert_eq!(file_pos(&buf), 100);

        buf.seek_virtual(VirtualOffset::new(150, 0)).unwrap();
        assert_eq!(file_pos(&buf), 150);
    }

    #[test]
    fn disjoint_ranges_data_streamed_correctly() {
        // With CHUNK_END_PAD = MAX_BLOCK_SIZE = 65536, chunks more than 65536
        // bytes apart stay disjoint after padding.
        let chunks = vec![
            Chunk { begin: VirtualOffset::new(100, 0), end: VirtualOffset::new(200, 0) },
            Chunk { begin: VirtualOffset::new(200_000, 0), end: VirtualOffset::new(200_100, 0) },
        ];
        let ranges = merge_chunks(&chunks);
        assert_eq!(ranges.len(), 2, "chunks should be disjoint");

        let big_file = fake_file(400_000);
        let mut cursor = std::io::Cursor::new(big_file.clone());
        let mut buf = RegionBuf::new(&mut cursor, &chunks).unwrap();

        // Position at the first range and pull a window — it must start at file[100].
        buf.seek_virtual(VirtualOffset::new(100, 0)).unwrap();
        buf.ensure_available(2).unwrap();
        assert_eq!(buf.window[buf.cursor], big_file[100]);
        assert_eq!(buf.window[buf.cursor + 1], big_file[101]);
    }

    #[test]
    fn disjoint_ranges_seek_to_second_range_is_correct() {
        // The streaming analog of the old concat regression: seeking to the
        // second range must reposition the window there (range_idx advances)
        // and expose file[200000], not garbage from a miscomputed offset.
        let big_file = fake_file(400_000);
        let mut cursor = std::io::Cursor::new(big_file.clone());

        let chunks = vec![
            Chunk { begin: VirtualOffset::new(100, 0), end: VirtualOffset::new(200, 0) },
            Chunk { begin: VirtualOffset::new(200_000, 0), end: VirtualOffset::new(200_100, 0) },
        ];
        assert_eq!(merge_chunks(&chunks).len(), 2, "should have 2 disjoint ranges");

        let mut buf = RegionBuf::new(&mut cursor, &chunks).unwrap();

        buf.seek_virtual(VirtualOffset::new(100, 0)).unwrap();
        buf.ensure_available(1).unwrap();
        assert_eq!(buf.range_idx, 0, "seek to range1");
        assert_eq!(file_pos(&buf), 100);
        assert_eq!(buf.window[buf.cursor], big_file[100]);

        buf.seek_virtual(VirtualOffset::new(200_000, 0)).unwrap();
        buf.ensure_available(1).unwrap();
        assert_eq!(buf.range_idx, 1, "seek to range2 advances range_idx");
        assert_eq!(file_pos(&buf), 200_000);
        assert_eq!(buf.window[buf.cursor], big_file[200_000]);
    }

    #[test]
    fn disjoint_ranges_seek_within_second_range() {
        let big_file = fake_file(400_000);
        let mut cursor = std::io::Cursor::new(big_file.clone());

        let chunks = vec![
            Chunk { begin: VirtualOffset::new(100, 0), end: VirtualOffset::new(200, 0) },
            Chunk { begin: VirtualOffset::new(200_000, 0), end: VirtualOffset::new(200_100, 0) },
        ];

        let mut buf = RegionBuf::new(&mut cursor, &chunks).unwrap();

        // Seek to file offset 200050 (50 bytes into second range)
        buf.seek_virtual(VirtualOffset::new(200_050, 0)).unwrap();
        buf.ensure_available(1).unwrap();
        assert_eq!(buf.range_idx, 1);
        assert_eq!(file_pos(&buf), 200_050);
        assert_eq!(buf.window[buf.cursor], big_file[200_050]);
    }

    #[test]
    fn disjoint_ranges_seek_before_loaded_region_fails() {
        let big_file = fake_file(400_000);
        let mut cursor = std::io::Cursor::new(big_file);

        let chunks =
            vec![Chunk { begin: VirtualOffset::new(1000, 0), end: VirtualOffset::new(2000, 0) }];

        let mut buf = RegionBuf::new(&mut cursor, &chunks).unwrap();

        // Seeking to file offset 500 (before the loaded range) must fail
        let result = buf.seek_virtual(VirtualOffset::new(500, 0));
        assert!(result.is_err(), "seek before loaded region should fail");
    }

    #[test]
    fn disjoint_ranges_seek_in_gap_between_ranges_fails() {
        let big_file = fake_file(400_000);
        let mut cursor = std::io::Cursor::new(big_file);

        let chunks = vec![
            Chunk { begin: VirtualOffset::new(100, 0), end: VirtualOffset::new(200, 0) },
            Chunk { begin: VirtualOffset::new(200_000, 0), end: VirtualOffset::new(200_100, 0) },
        ];

        let ranges = merge_chunks(&chunks);
        // Range1 ends at 200 + CHUNK_END_PAD. Range2 starts at 200_000.
        // The gap is between them. Pick an offset in that gap.
        let gap_offset = ranges[0].file_end + 1;
        assert!(
            gap_offset < ranges[1].file_start,
            "gap_offset {gap_offset} should be before range2 start {}",
            ranges[1].file_start
        );

        let mut buf = RegionBuf::new(&mut cursor, &chunks).unwrap();

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
        let mut buf = RegionBuf::new(&mut cursor, &chunks).unwrap();
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

    // --- BGZF round-trip properties ---

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
    /// Returns (`file_bytes`, `block_offsets`) where `block_offsets`[i] is the
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

    use hegel::prelude::*;

    /// Single contiguous range: load all blocks, read them sequentially,
    /// verify decompressed content matches original.
    #[hegel::test]
    fn single_range_roundtrip(tc: TestCase) {
        let n_blocks = tc.draw(gs::integers::<usize>().min_value(1).max_value(7));
        let block_size = tc.draw(gs::integers::<usize>().min_value(10).max_value(499));
        let seed = tc.draw(gs::integers::<u8>());
        // Build blocks with deterministic content
        let blocks: Vec<Vec<u8>> = (0..n_blocks)
            .map(|i| {
                (0..block_size).map(|j| seed.wrapping_add(i as u8).wrapping_add(j as u8)).collect()
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
        let mut buf = RegionBuf::new(&mut cursor, &chunks).unwrap();

        // Seek to the first block and read all data
        buf.seek_virtual(VirtualOffset::new(offsets[0], 0)).unwrap();

        let total_bytes: usize = blocks.iter().map(|b| b.len()).sum();
        let mut output = vec![0u8; total_bytes];
        buf.read_exact_into(&mut output).unwrap();

        // Verify
        let expected: Vec<u8> = blocks.iter().flatten().copied().collect();
        assert_eq!(output, expected);
    }

    /// Disjoint ranges: create blocks in two groups separated by padding,
    /// load both groups, seek to each and verify content.
    #[hegel::test]
    fn disjoint_ranges_roundtrip(tc: TestCase) {
        let group = || gs::integers::<usize>().min_value(1).max_value(3);
        let n_blocks_a = tc.draw(group());
        let n_blocks_b = tc.draw(group());
        let block_size = tc.draw(gs::integers::<usize>().min_value(10).max_value(299));
        let padding = tc.draw(gs::integers::<usize>().min_value(100_000).max_value(199_999));
        let seed = tc.draw(gs::integers::<u8>());
        let blocks_a: Vec<Vec<u8>> = (0..n_blocks_a)
            .map(|i| {
                (0..block_size).map(|j| seed.wrapping_add(i as u8).wrapping_add(j as u8)).collect()
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

        // Add padding to create a gap > CHUNK_END_PAD so chunks are disjoint
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
        assert_eq!(ranges.len(), 2);

        let mut cursor = std::io::Cursor::new(file);
        let mut buf = RegionBuf::new(&mut cursor, &[chunk_a, chunk_b]).unwrap();

        // Read group A
        buf.seek_virtual(VirtualOffset::new(offsets_a[0], 0)).unwrap();
        let total_a: usize = blocks_a.iter().map(|b| b.len()).sum();
        let mut out_a = vec![0u8; total_a];
        buf.read_exact_into(&mut out_a).unwrap();
        let expected_a: Vec<u8> = blocks_a.iter().flatten().copied().collect();
        assert_eq!(out_a, expected_a, "group A content mismatch");

        // Read group B
        buf.seek_virtual(VirtualOffset::new(group_b_start, 0)).unwrap();
        let total_b: usize = blocks_b.iter().map(|b| b.len()).sum();
        let mut out_b = vec![0u8; total_b];
        buf.read_exact_into(&mut out_b).unwrap();
        let expected_b: Vec<u8> = blocks_b.iter().flatten().copied().collect();
        assert_eq!(out_b, expected_b, "group B content mismatch");
    }

    /// Seek to mid-block positions (`within_block > 0`) should work correctly.
    #[hegel::test]
    fn within_block_seek(tc: TestCase) {
        let block_size = tc.draw(gs::integers::<usize>().min_value(20).max_value(499));
        // at most block_size-1, capped at 18 for simplicity
        let within = tc.draw(gs::integers::<usize>().min_value(1).max_value(18));
        let seed = tc.draw(gs::integers::<u8>());
        let data: Vec<u8> = (0..block_size).map(|j| seed.wrapping_add(j as u8)).collect();

        let (file, offsets) = make_bgzf_file(std::slice::from_ref(&data));

        let chunks = vec![Chunk {
            begin: VirtualOffset::new(offsets[0], 0),
            end: VirtualOffset::new(offsets[0] + 1, 0),
        }];

        let mut cursor = std::io::Cursor::new(file);
        let mut buf = RegionBuf::new(&mut cursor, &chunks).unwrap();

        let within_clamped = within.min(block_size - 1);
        buf.seek_virtual(VirtualOffset::new(offsets[0], within_clamped as u16)).unwrap();

        let remaining = block_size - within_clamped;
        let mut output = vec![0u8; remaining];
        buf.read_exact_into(&mut output).unwrap();

        assert_eq!(output, data[within_clamped..].to_vec());
    }

    /// CRC32 mismatch detection: corrupt a byte in the compressed data
    /// and verify decompression or CRC check fails.
    #[hegel::test]
    fn crc32_detects_corruption(tc: TestCase) {
        let block_size = tc.draw(gs::integers::<usize>().min_value(20).max_value(199));
        let seed = tc.draw(gs::integers::<u8>());
        let data: Vec<u8> = (0..block_size).map(|j| seed.wrapping_add(j as u8)).collect();

        let (mut file, offsets) = make_bgzf_file(std::slice::from_ref(&data));

        // Corrupt a byte in the compressed payload (after the 18-byte header)
        let corrupt_pos = offsets[0] as usize + 18;
        tc.assume(corrupt_pos < file.len() - 8);
        file[corrupt_pos] ^= 0xFF;

        let chunks = vec![Chunk {
            begin: VirtualOffset::new(offsets[0], 0),
            end: VirtualOffset::new(offsets[0] + 1, 0),
        }];

        let mut cursor = std::io::Cursor::new(file);
        let mut buf = RegionBuf::new(&mut cursor, &chunks).unwrap();
        buf.seek_virtual(VirtualOffset::new(offsets[0], 0)).unwrap();

        let mut output = vec![0u8; block_size];
        let err = buf
            .read_exact_into(&mut output)
            .expect_err("a corrupted deflate payload must not decode silently");
        // The flipped byte either breaks the deflate stream or survives it and
        // fails the block's CRC32. Anything else — a short read, a truncation
        // error — would mean the corruption was noticed for the wrong reason.
        assert!(
            matches!(
                err,
                BgzfError::DecompressionFailed { .. } | BgzfError::ChecksumMismatch { .. }
            ),
            "corrupt payload reported as {err:?}"
        );
        // Whatever it decoded, it must not be the original data.
        assert_ne!(output, data, "a corrupted block decoded to the original bytes");
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
        let mut region = RegionBuf::new(&mut cursor, &chunks).unwrap();
        region.seek_virtual(VirtualOffset::new(offsets[0], 0)).unwrap();

        let mut scratch: Vec<u8> = Vec::new();
        let result = region.read_record(&mut scratch).unwrap();

        assert_eq!(result, Some(body.as_slice()));
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
        let mut region = RegionBuf::new(&mut cursor, &chunks).unwrap();
        // Seek past the filler to where the length prefix starts (within_block = 60).
        region.seek_virtual(VirtualOffset::new(offsets[0], 60)).unwrap();

        let mut scratch: Vec<u8> = Vec::new();
        let result = region.read_record(&mut scratch).unwrap();

        assert_eq!(result, Some(body.as_slice()));
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
        let mut region = RegionBuf::new(&mut cursor, &chunks).unwrap();
        // Seek to where the first 2 bytes of the prefix begin (within_block = 62).
        region.seek_virtual(VirtualOffset::new(offsets[0], 62)).unwrap();

        let mut scratch: Vec<u8> = Vec::new();
        let result = region.read_record(&mut scratch).unwrap();

        assert_eq!(result, Some(body.as_slice()));
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
        let mut region = RegionBuf::new(&mut cursor, &chunks).unwrap();
        region.seek_virtual(VirtualOffset::new(offsets[0], 0)).unwrap();

        let mut scratch = Vec::new();
        for expected in &records {
            let got = region.read_record(&mut scratch).unwrap();
            assert_eq!(got, Some(expected.as_slice()));
        }
    }

    /// Varied record counts and body sizes, all packed into a single
    /// BGZF block. Every record must round-trip correctly.
    #[hegel::test]
    fn read_record_roundtrip(tc: TestCase) {
        let n_records = tc.draw(gs::integers::<usize>().min_value(1).max_value(9));
        let body_size = tc.draw(gs::integers::<usize>().min_value(4).max_value(199));
        let seed = tc.draw(gs::integers::<u8>());
        {
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
            let mut region = RegionBuf::new(&mut cursor, &chunks).unwrap();
            region.seek_virtual(VirtualOffset::new(offsets[0], 0)).unwrap();

            let mut scratch = Vec::new();
            for (i, expected) in records.iter().enumerate() {
                let got = region
                    .read_record(&mut scratch)
                    .unwrap_or_else(|e| panic!("record {i} failed: {e}"));
                assert_eq!(got, Some(expected.as_slice()));
            }
        }
    }

    /// Zero-length body: a record with `block_size=0` should be read back as
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
        let mut region = RegionBuf::new(&mut cursor, &chunks).unwrap();
        region.seek_virtual(VirtualOffset::new(offsets[0], 0)).unwrap();

        let mut scratch = Vec::new();
        let result = region.read_record(&mut scratch).unwrap();
        assert_eq!(result, Some(&[] as &[u8]));
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
        let mut region = RegionBuf::new(&mut cursor, &chunks).unwrap();
        region.seek_virtual(VirtualOffset::new(offsets[0], 0)).unwrap();

        let mut scratch = Vec::new();
        assert!(region.read_record(&mut scratch).is_err());
    }

    // --- RegionTooLarge / merged_byte_size tests ---

    #[test]
    fn merged_byte_size_single_chunk() {
        let chunks =
            vec![Chunk { begin: VirtualOffset::new(1000, 0), end: VirtualOffset::new(2000, 0) }];
        let size = merged_byte_size(&chunks);
        // end block_offset=2000, extended by CHUNK_END_PAD
        // So range is [1000, 2000+CHUNK_END_PAD)
        assert_eq!(size, 1000 + CHUNK_END_PAD);
    }

    #[test]
    fn merged_byte_size_overlapping_chunks_merge() {
        // Two chunks whose block ranges overlap after CHUNK_END_PAD extension.
        // Place them far enough apart to be disjoint with the old MAX_BLOCK_SIZE
        // padding, but they now overlap with CHUNK_END_PAD.
        let chunks = vec![
            Chunk { begin: VirtualOffset::new(1000, 0), end: VirtualOffset::new(2000, 0) },
            Chunk { begin: VirtualOffset::new(50_000, 0), end: VirtualOffset::new(60_000, 0) },
        ];
        let size = merged_byte_size(&chunks);
        // Range 1: [1000, 2000+CHUNK_END_PAD)
        // Range 2: [50000, 60000+CHUNK_END_PAD)
        // 50000 < 2000+CHUNK_END_PAD, so they merge to [1000, 60000+CHUNK_END_PAD)
        assert_eq!(size, 59000 + CHUNK_END_PAD);
    }

    #[test]
    fn merged_byte_size_disjoint_chunks() {
        // Place chunks far enough apart that they remain disjoint even with CHUNK_END_PAD.
        let far = 10_000_000u64;
        let chunks = vec![
            Chunk { begin: VirtualOffset::new(1000, 0), end: VirtualOffset::new(2000, 0) },
            Chunk { begin: VirtualOffset::new(far, 0), end: VirtualOffset::new(far + 1000, 0) },
        ];
        let size = merged_byte_size(&chunks);
        // Disjoint: sum of two ranges
        let range1 = 1000 + CHUNK_END_PAD;
        let range2 = 1000 + CHUNK_END_PAD;
        assert_eq!(size, range1 + range2);
    }

    // r[verify region_buf.window_budget]
    /// `new` is lazy: planning a huge region reads nothing up front. The window
    /// only fills (bounded) when decoding actually advances.
    #[test]
    fn new_is_lazy_for_oversized_region() {
        let far_end = (WINDOW_BUDGET as u64) * 8 + 1_000_000;
        let chunks =
            vec![Chunk { begin: VirtualOffset::new(0, 0), end: VirtualOffset::new(far_end, 0) }];

        let mut cursor = std::io::Cursor::new(vec![0u8; 100]);
        let buf = RegionBuf::new(&mut cursor, &chunks).expect("plan oversized region");
        assert!(buf.window.is_empty(), "new must not eagerly read the region");
    }

    // r[verify io.fuzz.alloc_limits]
    /// A corrupt index entry whose virtual offset points far past EOF must not
    /// trigger a huge allocation: refills are bounded by the window budget and
    /// clamped to the real file size.
    #[test]
    fn refill_caps_allocation_when_range_exceeds_file_size() {
        // Virtual offset block_offsets near the 48-bit max — would request
        // ~280 TiB without the cap.
        let huge = 0xffff_0801_0000u64;
        let chunks =
            vec![Chunk { begin: VirtualOffset::new(0, 0), end: VirtualOffset::new(huge, 0) }];

        let mut cursor = std::io::Cursor::new(vec![0u8; 64]);
        let mut buf = RegionBuf::new(&mut cursor, &chunks).expect("must not crash");
        buf.seek_virtual(VirtualOffset::new(0, 0)).unwrap();
        // Force a refill; the window must stay bounded by the 64-byte file.
        let avail = buf.ensure_available(BGZF_HEADER_SIZE).unwrap();
        assert!(buf.window.len() <= 64, "window bounded by file size, got {}", buf.window.len());
        assert!(avail <= 64);
    }

    /// Seeking backward to an offset already consumed (now behind the window)
    /// must reposition and re-read it — the window is not append-only.
    #[test]
    fn seek_virtual_backward_within_range_rereads() {
        let blocks: Vec<Vec<u8>> =
            (0..4u8).map(|i| (0..200u16).map(|j| i.wrapping_add(j as u8)).collect()).collect();
        let (file, offsets) = make_bgzf_file(&blocks);
        let last = *offsets.last().unwrap();
        let chunks = vec![Chunk {
            begin: VirtualOffset::new(offsets[0], 0),
            end: VirtualOffset::new(last + 1, 0),
        }];

        let mut cursor = std::io::Cursor::new(file);
        // Tiny budget so reading forward compacts the early blocks out of the window.
        let mut buf = RegionBuf::with_budget(&mut cursor, &chunks, MAX_BLOCK_SIZE).unwrap();

        // Read the last block to advance the window well past block 0.
        buf.seek_virtual(VirtualOffset::new(last, 0)).unwrap();
        let mut tail = vec![0u8; blocks[3].len()];
        buf.read_exact_into(&mut tail).unwrap();
        assert_eq!(tail, blocks[3]);

        // Seek back to block 0 (now behind window_file_start) and re-read it.
        buf.seek_virtual(VirtualOffset::new(offsets[0], 0)).unwrap();
        assert_eq!(buf.range_idx, 0);
        let mut head = vec![0u8; blocks[0].len()];
        buf.read_exact_into(&mut head).unwrap();
        assert_eq!(head, blocks[0], "backward seek must re-read the earlier block");
    }

    // --- Window-budget parity against a trivial oracle ---

    /// Budgets that force many refills (≥ one block floor), a few blocks per
    /// window, and a single huge window (the old "load everything" behavior).
    const PARITY_BUDGETS: [usize; 3] = [MAX_BLOCK_SIZE, 200_000, usize::MAX / 2];

    /// Decode every block of a single-range region via `read_exact_into` and
    /// return the concatenated decompressed bytes.
    fn stream_all_bytes(file: &[u8], chunks: &[Chunk], budget: usize, total: usize) -> Vec<u8> {
        let mut cursor = std::io::Cursor::new(file.to_vec());
        let mut buf = RegionBuf::with_budget(&mut cursor, chunks, budget).unwrap();
        buf.seek_virtual(chunks[0].begin).unwrap();
        let mut out = vec![0u8; total];
        buf.read_exact_into(&mut out).unwrap();
        out
    }

    /// Decoded output (and per-block virtual offsets) must be byte-identical
    /// no matter the window budget — windowing only changes *when* bytes are
    /// read, never *what* is produced. The oracle is the original payloads.
    #[hegel::test]
    fn window_budget_parity(tc: TestCase) {
        let n_blocks = tc.draw(gs::integers::<usize>().min_value(1).max_value(11));
        let block_size = tc.draw(gs::integers::<usize>().min_value(1).max_value(3_999));
        let seed = tc.draw(gs::integers::<u8>());
        let blocks: Vec<Vec<u8>> = (0..n_blocks)
            .map(|i| {
                (0..block_size)
                    .map(|j| seed.wrapping_add(i as u8).wrapping_mul(31).wrapping_add(j as u8))
                    .collect()
            })
            .collect();
        let (file, offsets) = make_bgzf_file(&blocks);
        let last = *offsets.last().unwrap();
        let chunks = vec![Chunk {
            begin: VirtualOffset::new(offsets[0], 0),
            end: VirtualOffset::new(last + 1, 0),
        }];
        let oracle: Vec<u8> = blocks.iter().flatten().copied().collect();

        for budget in PARITY_BUDGETS {
            let got = stream_all_bytes(&file, &chunks, budget, oracle.len());
            assert_eq!(&got, &oracle, "budget {} mismatch", budget);
        }

        // Per-block virtual offset (= true file offset) must be stable across
        // budgets: read block-by-block and compare the recorded block_offset.
        let mut cursor = std::io::Cursor::new(file.clone());
        let mut buf = RegionBuf::with_budget(&mut cursor, &chunks, MAX_BLOCK_SIZE).unwrap();
        buf.seek_virtual(VirtualOffset::new(offsets[0], 0)).unwrap();
        for &expected_off in &offsets {
            assert!(buf.read_block().unwrap());
            assert_eq!(buf.block_offset, expected_off, "block_offset drift across refill");
        }
    }

    /// Disjoint ranges with a tiny budget must still decode both groups
    /// correctly (each range refills independently).
    #[hegel::test]
    fn disjoint_budget_parity(tc: TestCase) {
        let group = || gs::integers::<usize>().min_value(1).max_value(3);
        let n_a = tc.draw(group());
        let n_b = tc.draw(group());
        let block_size = tc.draw(gs::integers::<usize>().min_value(10).max_value(299));
        let padding = tc.draw(gs::integers::<usize>().min_value(100_000).max_value(199_999));
        let seed = tc.draw(gs::integers::<u8>());
        let mk = |base: u8, n: usize| -> Vec<Vec<u8>> {
            (0..n)
                .map(|i| {
                    (0..block_size)
                        .map(|j| {
                            seed.wrapping_add(base).wrapping_add(i as u8).wrapping_add(j as u8)
                        })
                        .collect()
                })
                .collect()
        };
        let blocks_a = mk(0, n_a);
        let blocks_b = mk(100, n_b);

        let (mut file, offsets_a) = make_bgzf_file(&blocks_a);
        let pad_start = file.len();
        file.resize(pad_start + padding, 0);
        let group_b_start = file.len() as u64;
        let mut offsets_b = Vec::new();
        for block in &blocks_b {
            offsets_b.push(file.len() as u64);
            file.extend_from_slice(&make_bgzf_block(block));
        }
        file.extend_from_slice(&make_bgzf_eof());

        let chunk_a = Chunk {
            begin: VirtualOffset::new(offsets_a[0], 0),
            end: VirtualOffset::new(*offsets_a.last().unwrap() + 1, 0),
        };
        let chunk_b = Chunk {
            begin: VirtualOffset::new(group_b_start, 0),
            end: VirtualOffset::new(*offsets_b.last().unwrap() + 1, 0),
        };
        assert_eq!(merge_chunks(&[chunk_a, chunk_b]).len(), 2);

        let mut cursor = std::io::Cursor::new(file);
        let mut buf =
            RegionBuf::with_budget(&mut cursor, &[chunk_a, chunk_b], MAX_BLOCK_SIZE).unwrap();

        buf.seek_virtual(VirtualOffset::new(offsets_a[0], 0)).unwrap();
        let total_a: usize = blocks_a.iter().map(|b| b.len()).sum();
        let mut out_a = vec![0u8; total_a];
        buf.read_exact_into(&mut out_a).unwrap();
        let exp_a: Vec<u8> = blocks_a.iter().flatten().copied().collect();
        assert_eq!(out_a, exp_a, "group A");

        buf.seek_virtual(VirtualOffset::new(group_b_start, 0)).unwrap();
        let total_b: usize = blocks_b.iter().map(|b| b.len()).sum();
        let mut out_b = vec![0u8; total_b];
        buf.read_exact_into(&mut out_b).unwrap();
        let exp_b: Vec<u8> = blocks_b.iter().flatten().copied().collect();
        assert_eq!(out_b, exp_b, "group B");
    }

    // r[verify region_buf.virtual_offset+2]
    /// Reading a full 64 KiB block to its end leaves the cursor at the next
    /// block's first byte — the offset `BgzfWriter` gives the same position —
    /// with or without a cache, and whether or not the next block is resident.
    #[test]
    fn end_of_full_block_is_next_block_start() {
        let full: Vec<u8> = (0..MAX_BLOCK_SIZE).map(|i| (i % 251) as u8).collect();
        let (file, offsets) = make_bgzf_file(&[full.clone(), vec![7u8; 100]]);
        let chunks = [Chunk {
            begin: VirtualOffset::new(offsets[0], 0),
            end: VirtualOffset::new(offsets[1], 50),
        }];
        let mut cache = BlockCache::new();
        for pass in 0..3 {
            let mut cursor = std::io::Cursor::new(file.clone());
            let mut buf = if pass == 0 {
                RegionBuf::new(&mut cursor, &chunks).unwrap()
            } else {
                RegionBuf::with_cache(&mut cursor, &chunks, &mut cache).unwrap()
            };
            buf.seek_virtual(chunks[0].begin).unwrap();
            let mut out = vec![0u8; full.len()];
            buf.read_exact_into(&mut out).unwrap();
            assert_eq!(out, full);
            assert_eq!(buf.virtual_offset(), VirtualOffset::new(offsets[1], 0), "pass {pass}");
        }
    }

    // --- Block cache: cached queries against uncached ones ---

    /// A virtual offset strictly inside one of `blocks`.
    #[hegel::composite]
    fn arb_point(tc: &TestCase, lens: Vec<usize>) -> (usize, usize) {
        let block = tc.draw(gs::integers::<usize>().max_value(lens.len() - 1));
        let within = tc.draw(gs::integers::<usize>().max_value(lens[block] - 1));
        (block, within)
    }

    /// Read `chunks` the way `BamQuery` walks them: seek to each begin, take
    /// bytes until the cursor reaches the end. Returns every `(virtual offset,
    /// bytes)` piece handed out.
    fn read_chunks<R: Read + Seek>(
        buf: &mut RegionBuf<'_, R>,
        chunks: &[Chunk],
    ) -> Vec<(VirtualOffset, Vec<u8>)> {
        let mut out = Vec::new();
        for c in chunks {
            buf.seek_virtual(c.begin).unwrap();
            loop {
                let (voff, data) = buf.fill_buf().unwrap();
                if data.is_empty() || voff >= c.end {
                    break;
                }
                let n = if voff.block_offset() == c.end.block_offset() {
                    usize::from(c.end.within_block() - voff.within_block())
                } else {
                    data.len()
                };
                let n = n.min(data.len());
                out.push((voff, data[..n].to_vec()));
                buf.consume(n);
            }
        }
        out
    }

    // r[verify region_buf.block_cache]
    /// A sequence of queries sharing one cache — repeats, backward jumps,
    /// disjoint ranges, a cache small enough to evict — hands out exactly the
    /// pieces a fresh uncached buffer hands out for each query, and those
    /// pieces are the original payloads.
    #[hegel::test]
    fn cached_queries_match_uncached(tc: TestCase) {
        let n_blocks = tc.draw(gs::integers::<usize>().min_value(1).max_value(24));
        let seed = tc.draw(gs::integers::<u64>());
        let mut state = seed | 1;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        // Incompressible payloads, some large, so that blocks far apart in the
        // file land in disjoint planned ranges.
        let blocks: Vec<Vec<u8>> = (0..n_blocks)
            .map(|_| {
                let len = if next() % 4 == 0 { 20_000 } else { 1 + (next() % 3_000) as usize };
                (0..len).map(|_| next() as u8).collect()
            })
            .collect();
        let lens: Vec<usize> = blocks.iter().map(Vec::len).collect();
        let (file, offsets) = make_bgzf_file(&blocks);
        let voff = |(b, w): (usize, usize)| VirtualOffset::new(offsets[b], w as u16);

        let budget = PARITY_BUDGETS[tc.draw(gs::integers::<usize>().max_value(2))];
        let capacity = tc.draw(gs::integers::<usize>().min_value(1).max_value(8));
        let mut cache = BlockCache::with_capacity(capacity);
        let mut cached_file = std::io::Cursor::new(file.clone());

        let n_queries = tc.draw(gs::integers::<usize>().min_value(1).max_value(10));
        for _ in 0..n_queries {
            // Sorted, non-overlapping chunks: pair up consecutive sorted points.
            let n_chunks = tc.draw(gs::integers::<usize>().min_value(1).max_value(3));
            let mut points: Vec<(usize, usize)> =
                (0..2 * n_chunks).map(|_| tc.draw(arb_point(lens.clone()))).collect();
            points.sort_unstable();
            let pairs = points.as_chunks::<2>().0;
            let chunks: Vec<Chunk> =
                pairs.iter().map(|&[a, b]| Chunk { begin: voff(a), end: voff(b) }).collect();

            let mut uncached_file = std::io::Cursor::new(file.clone());
            let mut plain = RegionBuf::with_budget(&mut uncached_file, &chunks, budget).unwrap();
            let expected = read_chunks(&mut plain, &chunks);
            drop(plain);

            let mut buf =
                RegionBuf::with_cache_and_budget(&mut cached_file, &chunks, &mut cache, budget)
                    .unwrap();
            let got = read_chunks(&mut buf, &chunks);
            drop(buf);
            assert_eq!(got, expected, "chunks {chunks:?}");

            let mut oracle: Vec<u8> = Vec::new();
            for &[(b0, w0), (b1, w1)] in pairs {
                for (b, block) in blocks.iter().enumerate().take(b1 + 1).skip(b0) {
                    let from = if b == b0 { w0 } else { 0 };
                    let to = if b == b1 { w1 } else { block.len() };
                    oracle.extend_from_slice(&block[from..to]);
                }
            }
            let flat: Vec<u8> = got.into_iter().flat_map(|(_, bytes)| bytes).collect();
            assert_eq!(flat, oracle, "payload bytes for {chunks:?}");
            assert!(cache.entries.len() <= capacity);
        }
    }
}
