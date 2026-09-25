//! BGZF block writer. [`BgzfWriter`] compresses data into independent BGZF blocks
//! via libdeflate, tracking virtual offsets for index co-production.

use super::bgzf::{BgzfError, VirtualOffset};
use crate::io::IndexBuilder;
use std::collections::VecDeque;
use std::io::Write;
use std::panic::AssertUnwindSafe;
use std::sync::{Arc, Mutex, mpsc};
use std::thread;
use tracing::warn;

// r[impl bgzf.writer.block_size]
/// Maximum uncompressed data per BGZF block: htslib's `BGZF_BLOCK_SIZE`.
///
/// Not 64 KiB. A block, header and footer included, must fit in 64 KiB, and
/// 64 KiB of data that does not compress (level 0 always, random bytes at any
/// level) stores as 64 KiB plus DEFLATE's stored-block overhead — no room for
/// the 26 bytes of gzip framing. At `0xff00` the worst case fits.
const MAX_UNCOMPRESSED_SIZE: usize = 0xff00;

/// Gzip member header (with the BC subfield) and footer (CRC32 + ISIZE) sizes.
const HEADER_LEN: usize = 18;
const FOOTER_LEN: usize = 8;

/// Largest whole block, framing included: BSIZE is a u16 holding size − 1.
const MAX_BLOCK_LEN: usize = 1 << 16;

/// libdeflate's worst case for `n` input bytes, as
/// `libdeflate_deflate_compress_bound` computes it: all stored blocks, each of
/// at least 5000 bytes and with 5 bytes of framing.
#[allow(clippy::arithmetic_side_effects, reason = "const-evaluated; overflow fails the build")]
const fn deflate_bound(n: usize) -> usize {
    let blocks = n.div_ceil(5000);
    n + 5 * if blocks == 0 { 1 } else { blocks }
}

// r[impl bgzf.writer.block_size]
// Whatever the data and level, a full buffer compresses into a block that fits.
const _: () = assert!(
    HEADER_LEN + deflate_bound(MAX_UNCOMPRESSED_SIZE) + FOOTER_LEN <= MAX_BLOCK_LEN,
    "a full block's worst case must fit in 64 KiB"
);

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

/// Length of a buffer that holds any block `compressor` can produce from a
/// full buffer: header, worst-case payload, footer.
fn max_block_len(compressor: &mut libdeflater::Compressor) -> usize {
    HEADER_LEN
        .saturating_add(compressor.deflate_compress_bound(MAX_UNCOMPRESSED_SIZE))
        .saturating_add(FOOTER_LEN)
}

// r[impl bgzf.writer.single_write]
/// Compress `data` into one complete BGZF block at the front of `block` —
/// the payload straight after the header, then the footer — and return the
/// block's length. `block` must be at least [`max_block_len`] long.
fn compress_block(
    compressor: &mut libdeflater::Compressor,
    data: &[u8],
    block: &mut [u8],
) -> Result<usize, BgzfError> {
    let payload_end = block.len().saturating_sub(FOOTER_LEN);
    let payload = block.get_mut(HEADER_LEN..payload_end).ok_or(BgzfError::CorruptHeader)?;
    let compressed_len = compressor
        .deflate_compress(data, payload)
        .map_err(|source| BgzfError::CompressionFailed { source })?;

    let mut crc = libdeflater::Crc::new();
    crc.update(data);
    let isize_val = u32::try_from(data.len()).map_err(|_| BgzfError::CorruptHeader)?;

    // total = header(18) + compressed_len + footer(8); BSIZE = total - 1
    let footer_start = HEADER_LEN.saturating_add(compressed_len);
    let total_block_size = footer_start.saturating_add(FOOTER_LEN);
    let bsize = total_block_size.saturating_sub(1);
    // r[impl bgzf.writer.block_size]
    let bsize_bytes = u16::try_from(bsize)
        .map_err(|_| BgzfError::BlockTooLarge { size: total_block_size })?
        .to_le_bytes();

    let header = block.get_mut(..HEADER_LEN).ok_or(BgzfError::CorruptHeader)?;
    header.copy_from_slice(&BGZF_HEADER);
    // Bytes 16-17 are the BSIZE field
    #[allow(clippy::indexing_slicing, reason = "fixed-size header with known offsets")]
    {
        header[16] = bsize_bytes[0];
        header[17] = bsize_bytes[1];
    }
    let footer = block.get_mut(footer_start..total_block_size).ok_or(BgzfError::CorruptHeader)?;
    let (crc_bytes, isize_bytes) = footer.split_at_mut(4);
    crc_bytes.copy_from_slice(&crc.sum().to_le_bytes());
    isize_bytes.copy_from_slice(&isize_val.to_le_bytes());
    Ok(total_block_size)
}

// r[impl bgzf.writer]
// r[impl bgzf.writer.buffer]
// r[impl bgzf.writer.compression]
// r[impl bgzf.writer.virtual_offset]
/// BGZF block writer that compresses data into independent gzip blocks.
///
/// Accumulates up to 65280 bytes of uncompressed data, then compresses and emits a
/// BGZF block. Tracks virtual offsets for index co-production.
pub struct BgzfWriter<W: Write> {
    inner: Option<W>,
    /// Uncompressed data buffer (up to `MAX_UNCOMPRESSED_SIZE`).
    buf: Vec<u8>,
    // r[impl bgzf.writer.single_write]
    /// One whole block — header, DEFLATE payload, footer — assembled in place.
    /// Sized once for the largest payload a full buffer can compress to and
    /// never shrunk, so no block pays for zero-filling it.
    block: Vec<u8>,
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
        let mut compressor = libdeflater::Compressor::new(
            libdeflater::CompressionLvl::new(level).unwrap_or_default(),
        );
        let block_len = max_block_len(&mut compressor);
        Self {
            inner: Some(inner),
            buf: Vec::with_capacity(MAX_UNCOMPRESSED_SIZE),
            block: vec![0; block_len],
            compressor,
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
        // `write_all` flushes the moment the buffer reaches MAX_UNCOMPRESSED_SIZE, so
        // the length is always representable. Clamping instead would name the block's
        // last byte — a position one short of the true one, and in the middle of the
        // record that just ended.
        let within = u16::try_from(self.buf.len()).unwrap_or_else(|_| {
            warn!(
                "BgzfWriter: {} buffered bytes exceed a BGZF block; virtual offset clamped",
                self.buf.len()
            );
            u16::MAX
        });
        VirtualOffset::new(self.block_offset, within)
    }

    // r[impl bgzf.writer.flush_if_needed]
    /// Flush the current block if `upcoming_bytes` would exceed the block limit.
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
        loop {
            let space = MAX_UNCOMPRESSED_SIZE.saturating_sub(self.buf.len());
            let take = remaining.len().min(space);
            #[allow(clippy::indexing_slicing, reason = "take <= remaining.len()")]
            {
                self.buf.extend_from_slice(&remaining[..take]);
                remaining = &remaining[take..];
            }
            // r[impl bgzf.writer.buffer]
            // Flush the instant the buffer is full, before returning to the caller:
            // a within-block offset of 65536 does not exist, and the position after
            // that byte is the *next* block's (offset, 0). Leaving the buffer full
            // would make `virtual_offset()` name a byte inside the record just
            // written. Block boundaries are unchanged — the next write would have
            // flushed here anyway.
            if self.buf.len() >= MAX_UNCOMPRESSED_SIZE {
                self.flush_block()?;
            }
            if remaining.is_empty() {
                return Ok(());
            }
        }
    }

    // r[impl bgzf.writer.single_write]
    /// Compress and emit the current buffer as a BGZF block.
    fn flush_block(&mut self) -> Result<(), BgzfError> {
        if self.buf.is_empty() {
            return Ok(());
        }
        let total_block_size = compress_block(&mut self.compressor, &self.buf, &mut self.block)?;
        let block = self.block.get(..total_block_size).ok_or(BgzfError::CorruptHeader)?;
        let w = self.inner.as_mut().ok_or(BgzfError::AlreadyFinished)?;
        w.write_all(block).map_err(|source| BgzfError::WriteFailed { source })?;

        // Advance block_offset by the total compressed block size
        self.block_offset = self
            .block_offset
            .checked_add(total_block_size as u64)
            .ok_or(BgzfError::CorruptHeader)?;
        self.buf.clear();

        Ok(())
    }

    /// The inner writer back, provided nothing has been written to it yet.
    fn into_unwritten_inner(mut self) -> Result<W, BgzfError> {
        if !self.buf.is_empty() || self.block_offset != 0 {
            return Err(BgzfError::AlreadyFinished);
        }
        self.inner.take().ok_or(BgzfError::AlreadyFinished)
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

// ── Parallel compression ────────────────────────────────────────────────

/// A block handed to a compression worker, with the buffer to compress it into.
struct Job {
    index: u64,
    data: Vec<u8>,
    block: Vec<u8>,
}

/// A worker's answer: both buffers come back so the writer can reuse them.
struct Done {
    index: u64,
    data: Vec<u8>,
    block: Vec<u8>,
    len: Result<usize, BgzfError>,
}

// r[impl bgzf.writer.parallel]
// r[impl bgzf.writer.parallel.identical_output]
/// BGZF writer that compresses blocks on `threads` worker threads and writes
/// them to the inner stream in order, from the calling thread.
///
/// Blocks are cut exactly where [`BgzfWriter`] cuts them — the cuts depend on
/// uncompressed sizes only — and every worker compresses at the same level, so
/// the output is byte-identical to the serial writer's.
///
/// A block's file offset is only known once every block before it has been
/// compressed, so this writer hands out *index offsets* instead of virtual
/// offsets: `(block number << 16) | within`. They order exactly as the real
/// ones do, which is all [`IndexBuilder`](crate::io::IndexBuilder) relies on
/// while it builds; [`finish`](Self::finish) translates them once every block
/// is on disk (see `r[bgzf.writer.parallel.index_offsets]`).
pub(crate) struct ParallelBgzfWriter<W: Write> {
    inner: Option<W>,
    buf: Vec<u8>,
    // The channels and thread handles are wrapped in `AssertUnwindSafe` so the
    // writers keep the auto traits they had before the pool existed; nothing a
    // panic could leave half-updated lives in them.
    jobs: AssertUnwindSafe<Option<mpsc::SyncSender<Job>>>,
    /// Behind a `Mutex` only so the writer stays `Sync`, as the serial one is;
    /// just the owning thread ever locks it.
    done: AssertUnwindSafe<Mutex<mpsc::Receiver<Done>>>,
    workers: AssertUnwindSafe<Vec<thread::JoinHandle<()>>>,
    /// Compressed blocks that arrived before an earlier one; slot `i` holds
    /// block `next_write + i`.
    ready: VecDeque<Option<(Vec<u8>, usize)>>,
    /// Number of the block `buf` will become.
    next_submit: u64,
    /// Number of the next block to write to `inner`.
    next_write: u64,
    max_in_flight: u64,
    spare_data: Vec<Vec<u8>>,
    spare_blocks: Vec<Vec<u8>>,
    block_len: usize,
    /// File offset after the last block written.
    written: u64,
    // r[impl bgzf.writer.parallel.index_offsets]
    /// File offset of every block written so far, by block number.
    block_starts: Vec<u64>,
    /// Set once a block failed to compress or write. Its place in the stream
    /// is lost for good, so nothing may be written — or waited for — after it.
    failed: bool,
}

impl<W: Write> ParallelBgzfWriter<W> {
    /// Start `threads` (≥ 1) compression workers at `level`.
    pub(crate) fn new(inner: W, level: i32, threads: usize) -> Result<Self, BgzfError> {
        let threads = threads.max(1);
        let lvl = libdeflater::CompressionLvl::new(level).unwrap_or_default();
        let block_len = max_block_len(&mut libdeflater::Compressor::new(lvl));
        // Enough blocks queued that no worker waits while the writer thread
        // is busy writing, without holding more than a few MiB.
        let max_in_flight = u64::try_from(threads.saturating_mul(4)).unwrap_or(u64::MAX);
        let (job_tx, job_rx) = mpsc::sync_channel::<Job>(threads.saturating_mul(4));
        let (done_tx, done_rx) = mpsc::channel::<Done>();
        let job_rx = Arc::new(Mutex::new(job_rx));
        let mut workers = Vec::with_capacity(threads);
        for i in 0..threads {
            let jobs = Arc::clone(&job_rx);
            let done = done_tx.clone();
            let handle = thread::Builder::new()
                .name(format!("seqair-bgzf-{i}"))
                .spawn(move || compress_worker(lvl, &jobs, &done))
                .map_err(|source| BgzfError::ThreadSpawn { source })?;
            workers.push(handle);
        }
        Ok(Self {
            inner: Some(inner),
            buf: Vec::with_capacity(MAX_UNCOMPRESSED_SIZE),
            jobs: AssertUnwindSafe(Some(job_tx)),
            done: AssertUnwindSafe(Mutex::new(done_rx)),
            workers: AssertUnwindSafe(workers),
            ready: VecDeque::new(),
            next_submit: 0,
            next_write: 0,
            max_in_flight,
            spare_data: Vec::new(),
            spare_blocks: Vec::new(),
            block_len,
            written: 0,
            block_starts: Vec::new(),
            failed: false,
        })
    }

    /// The position after the last byte written, as an index offset: the
    /// current block's *number* in the upper 48 bits.
    pub(crate) fn index_offset(&self) -> VirtualOffset {
        // Same invariant as `BgzfWriter::virtual_offset`: the buffer is flushed
        // the moment it fills, so its length fits a within-block offset.
        let within = u16::try_from(self.buf.len()).unwrap_or(u16::MAX);
        debug_assert!(self.buf.len() < MAX_UNCOMPRESSED_SIZE, "buffer is flushed when full");
        VirtualOffset::new(self.next_submit, within)
    }

    pub(crate) fn flush_if_needed(&mut self, upcoming_bytes: usize) -> Result<(), BgzfError> {
        if self.buf.len().saturating_add(upcoming_bytes) > MAX_UNCOMPRESSED_SIZE {
            self.flush_block()?;
        }
        Ok(())
    }

    /// Same block cutting as [`BgzfWriter::write_all`].
    pub(crate) fn write_all(&mut self, data: &[u8]) -> Result<(), BgzfError> {
        let mut remaining = data;
        loop {
            let space = MAX_UNCOMPRESSED_SIZE.saturating_sub(self.buf.len());
            let (now, later) = remaining.split_at(remaining.len().min(space));
            self.buf.extend_from_slice(now);
            remaining = later;
            if self.buf.len() >= MAX_UNCOMPRESSED_SIZE {
                self.flush_block()?;
            }
            if remaining.is_empty() {
                return Ok(());
            }
        }
    }

    /// Hand the buffer to a worker, then write whatever has come back in order.
    fn flush_block(&mut self) -> Result<(), BgzfError> {
        if self.buf.is_empty() {
            return Ok(());
        }
        if self.failed {
            return Err(BgzfError::CompressionWorkerLost);
        }
        let fresh =
            self.spare_data.pop().unwrap_or_else(|| Vec::with_capacity(MAX_UNCOMPRESSED_SIZE));
        let data = std::mem::replace(&mut self.buf, fresh);
        let block = self.spare_blocks.pop().unwrap_or_else(|| vec![0; self.block_len]);
        let jobs = self.jobs.as_ref().ok_or(BgzfError::AlreadyFinished)?;
        jobs.send(Job { index: self.next_submit, data, block })
            .map_err(|_| BgzfError::CompressionWorkerLost)?;
        self.next_submit = self.next_submit.saturating_add(1);
        self.pump(false)
    }

    /// Collect finished blocks and write the in-order prefix. Waits while more
    /// than `max_in_flight` blocks are outstanding, or — with `all` — until
    /// every submitted block is written.
    fn pump(&mut self, all: bool) -> Result<(), BgzfError> {
        let result = self.pump_inner(all);
        if result.is_err() {
            self.failed = true;
        }
        result
    }

    fn pump_inner(&mut self, all: bool) -> Result<(), BgzfError> {
        if self.failed {
            return Err(BgzfError::CompressionWorkerLost);
        }
        loop {
            let outstanding = self.next_submit.saturating_sub(self.next_write);
            let wait = if all { outstanding > 0 } else { outstanding >= self.max_in_flight };
            let rx = self.done.get_mut().map_err(|_| BgzfError::CompressionWorkerLost)?;
            let done = if wait {
                rx.recv().map_err(|_| BgzfError::CompressionWorkerLost)?
            } else {
                match rx.try_recv() {
                    Ok(done) => done,
                    Err(mpsc::TryRecvError::Empty) => return Ok(()),
                    Err(mpsc::TryRecvError::Disconnected) => {
                        return Err(BgzfError::CompressionWorkerLost);
                    }
                }
            };
            self.accept(done)?;
        }
    }

    fn accept(&mut self, done: Done) -> Result<(), BgzfError> {
        let Done { index, mut data, block, len } = done;
        let len = len?;
        data.clear();
        self.spare_data.push(data);
        let slot = usize::try_from(index.saturating_sub(self.next_write))
            .map_err(|_| BgzfError::CorruptHeader)?;
        if self.ready.len() <= slot {
            self.ready.resize_with(slot.saturating_add(1), || None);
        }
        let entry = self.ready.get_mut(slot).ok_or(BgzfError::CorruptHeader)?;
        *entry = Some((block, len));

        while let Some(Some(_)) = self.ready.front() {
            let Some(Some((block, len))) = self.ready.pop_front() else { break };
            let bytes = block.get(..len).ok_or(BgzfError::CorruptHeader)?;
            let w = self.inner.as_mut().ok_or(BgzfError::AlreadyFinished)?;
            w.write_all(bytes).map_err(|source| BgzfError::WriteFailed { source })?;
            self.block_starts.push(self.written);
            self.written = self.written.checked_add(len as u64).ok_or(BgzfError::CorruptHeader)?;
            self.next_write = self.next_write.saturating_add(1);
            self.spare_blocks.push(block);
        }
        Ok(())
    }

    // r[impl bgzf.writer.parallel.index_offsets]
    /// The virtual offset an index offset stands for, once its block is written.
    fn resolve(&self, offset: VirtualOffset) -> Option<VirtualOffset> {
        let block = usize::try_from(offset.block_offset()).ok()?;
        let start = match self.block_starts.get(block) {
            Some(&start) => start,
            // The position just past the last block written: where the next
            // block (or the EOF marker) starts.
            None if block == self.block_starts.len() => self.written,
            None => return None,
        };
        Some(VirtualOffset::new(start, offset.within_block()))
    }

    // r[impl bgzf.writer.eof_marker]
    // r[impl bgzf.writer.finish]
    // r[impl bgzf.writer.parallel.index_offsets]
    /// Flush and write every block, translate `index`'s offsets to file
    /// offsets, write the EOF marker and return the inner writer.
    pub(crate) fn finish(mut self, index: Option<&mut IndexBuilder>) -> Result<W, BgzfError> {
        self.flush_block()?;
        self.pump(true)?;
        if let Some(index) = index {
            index.map_offsets(|offset| {
                let resolved = self.resolve(offset);
                debug_assert!(resolved.is_some(), "{offset:?} names a block not yet written");
                resolved.unwrap_or(offset)
            });
        }
        let mut w = self.inner.take().ok_or(BgzfError::AlreadyFinished)?;
        w.write_all(&EOF_BLOCK).map_err(|source| BgzfError::WriteFailed { source })?;
        w.flush().map_err(|source| BgzfError::WriteFailed { source })?;
        self.shutdown();
        Ok(w)
    }

    /// The inner writer back, provided nothing has been written to it yet.
    fn into_unwritten_inner(mut self) -> Result<W, BgzfError> {
        if !self.buf.is_empty() || self.next_submit != 0 {
            return Err(BgzfError::AlreadyFinished);
        }
        let inner = self.inner.take().ok_or(BgzfError::AlreadyFinished)?;
        self.shutdown();
        Ok(inner)
    }

    /// Close the job queue and wait for the workers to exit.
    fn shutdown(&mut self) {
        *self.jobs = None;
        for worker in self.workers.drain(..) {
            if worker.join().is_err() {
                warn!("BGZF compression worker panicked");
            }
        }
    }
}

impl<W: Write> Drop for ParallelBgzfWriter<W> {
    fn drop(&mut self) {
        // Best-effort, like `BgzfWriter`: write out what was buffered and queued.
        if self.inner.is_some()
            && !self.failed
            && let Err(e) = self.flush_block().and_then(|()| self.pump(true))
        {
            warn!("ParallelBgzfWriter::drop: failed to flush buffered blocks: {e}");
        }
        self.shutdown();
    }
}

fn compress_worker(
    level: libdeflater::CompressionLvl,
    jobs: &Mutex<mpsc::Receiver<Job>>,
    done: &mpsc::Sender<Done>,
) {
    let mut compressor = libdeflater::Compressor::new(level);
    loop {
        let job = match jobs.lock() {
            Ok(rx) => rx.recv(),
            Err(_) => return,
        };
        let Ok(Job { index, data, mut block }) = job else { return };
        // A panic here must not leave the writer waiting for this block forever.
        let len = std::panic::catch_unwind(AssertUnwindSafe(|| {
            compress_block(&mut compressor, &data, &mut block)
        }))
        .unwrap_or(Err(BgzfError::CompressionWorkerLost));
        if done.send(Done { index, data, block, len }).is_err() {
            return;
        }
    }
}

// ── Serial or parallel ─────────────────────────────────────────────────

/// The BGZF stream behind the BAM and VCF/BCF writers: one thread, or a pool.
///
/// Offsets for index co-production come from [`index_offset`](Self::index_offset)
/// and are only final after [`finish`](Self::finish) has translated them.
pub(crate) enum BgzfSink<W: Write> {
    Serial(BgzfWriter<W>),
    Parallel(ParallelBgzfWriter<W>),
}

impl<W: Write> BgzfSink<W> {
    /// Serial for `threads == 0`, else a pool of `threads` workers.
    pub(crate) fn new(inner: W, level: i32, threads: usize) -> Result<Self, BgzfError> {
        if threads == 0 {
            Ok(Self::Serial(BgzfWriter::with_compression_level(inner, level)))
        } else {
            ParallelBgzfWriter::new(inner, level, threads).map(Self::Parallel)
        }
    }

    /// The same stream with `threads` workers instead, before anything is written.
    pub(crate) fn with_threads(self, level: i32, threads: usize) -> Result<Self, BgzfError> {
        let inner = match self {
            Self::Serial(w) => w.into_unwritten_inner()?,
            Self::Parallel(w) => w.into_unwritten_inner()?,
        };
        Self::new(inner, level, threads)
    }

    #[expect(
        clippy::same_name_method,
        reason = "inherent method is the concrete impl; BgzfWrite trait delegates to it for dyn dispatch"
    )]
    pub(crate) fn index_offset(&self) -> VirtualOffset {
        match self {
            Self::Serial(w) => w.virtual_offset(),
            Self::Parallel(w) => w.index_offset(),
        }
    }

    #[expect(
        clippy::same_name_method,
        reason = "inherent method is the concrete impl; BgzfWrite trait delegates to it for dyn dispatch"
    )]
    pub(crate) fn flush_if_needed(&mut self, upcoming_bytes: usize) -> Result<(), BgzfError> {
        match self {
            Self::Serial(w) => w.flush_if_needed(upcoming_bytes),
            Self::Parallel(w) => w.flush_if_needed(upcoming_bytes),
        }
    }

    #[expect(
        clippy::same_name_method,
        reason = "inherent method is the concrete impl; BgzfWrite trait delegates to it for dyn dispatch"
    )]
    pub(crate) fn write_all(&mut self, data: &[u8]) -> Result<(), BgzfError> {
        match self {
            Self::Serial(w) => w.write_all(data),
            Self::Parallel(w) => w.write_all(data),
        }
    }

    /// Finish the stream. `index`, built from [`index_offset`](Self::index_offset)s
    /// and already [`finish`](IndexBuilder::finish)ed, has its offsets translated
    /// to file offsets.
    pub(crate) fn finish(self, index: Option<&mut IndexBuilder>) -> Result<W, BgzfError> {
        match self {
            Self::Serial(w) => w.finish(),
            Self::Parallel(w) => w.finish(index),
        }
    }
}

#[allow(clippy::cast_possible_truncation, reason = "tests")]
#[allow(clippy::indexing_slicing, reason = "tests")]
#[allow(clippy::arithmetic_side_effects, reason = "tests")]
#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam::bgzf::BgzfReader;
    use hegel::prelude::*;
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

    /// Counts `write` calls and keeps the bytes.
    #[derive(Default)]
    struct CountingSink {
        writes: usize,
        bytes: Vec<u8>,
    }

    impl Write for CountingSink {
        fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
            self.writes += 1;
            self.bytes.extend_from_slice(buf);
            Ok(buf.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    // r[verify bgzf.writer.single_write]
    /// An unbuffered sink sees one `write` per block (plus one for the EOF
    /// marker), and the stream still decodes.
    #[test]
    fn one_write_per_block() {
        let data: Vec<u8> = (0..200_000).map(|i| (i % 251) as u8).collect();
        let mut writer = BgzfWriter::new(CountingSink::default());
        writer.write_all(&data).unwrap();
        let sink = writer.finish().unwrap();
        let blocks = data.len().div_ceil(MAX_UNCOMPRESSED_SIZE);
        assert_eq!(sink.writes, blocks + 1, "one write per block plus the EOF marker");
        assert_eq!(read_all(&sink.bytes), data);
    }

    // r[verify bgzf.writer.eof_marker]
    #[test]
    fn eof_marker_present() {
        let mut output = Vec::new();
        let writer = BgzfWriter::new(&mut output);
        writer.finish().unwrap();
        assert!(output.len() >= 28);
        assert_eq!(&output[output.len() - 28..], &EOF_BLOCK);
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

        // Fill the block exactly; the writer flushes at once (see below)
        let fill = vec![0u8; MAX_UNCOMPRESSED_SIZE - 9];
        writer.write_all(&fill).unwrap();
        writer.write_all(&[0x42]).unwrap();
        let voff2 = writer.virtual_offset();
        assert!(voff2.block_offset() > 0, "should have advanced to a new block");
        assert_eq!(voff2.within_block(), 1, "one byte in the new block");

        writer.finish().unwrap();
    }

    // r[verify bgzf.writer.buffer]
    // r[verify bgzf.writer.virtual_offset]
    #[test]
    fn exactly_full_buffer_flushes_before_returning() {
        let mut output = Vec::new();
        let mut writer = BgzfWriter::new(&mut output);

        writer.write_all(&vec![0u8; MAX_UNCOMPRESSED_SIZE]).unwrap();

        let voff = writer.virtual_offset();
        assert!(voff.block_offset() > 0, "an exactly-full buffer must have been flushed");
        assert_eq!(
            voff.within_block(),
            0,
            "the position after a block's last byte is the next block's (offset, 0), \
             never 65535 — that names a byte inside the record just written"
        );

        writer.finish().unwrap();
    }

    // r[verify bgzf.writer.virtual_offset]
    /// Every offset the writer hands out must resolve to the uncompressed byte
    /// count at that moment. This is what an index co-produced during writing
    /// records as a record boundary.
    #[test]
    fn virtual_offsets_resolve_to_the_uncompressed_position() {
        // Sizes chosen so one write lands exactly on the block boundary.
        let sizes =
            [40_000usize, MAX_UNCOMPRESSED_SIZE - 40_000, 1, 70_000, 300, MAX_UNCOMPRESSED_SIZE];
        let mut output = Vec::new();
        let mut writer = BgzfWriter::new(&mut output);

        let mut observed = Vec::new();
        let mut written = 0usize;
        let mut payload = Vec::new();
        for (i, &n) in sizes.iter().enumerate() {
            let chunk: Vec<u8> = (0..n).map(|j| (j.wrapping_add(i) & 0xFF) as u8).collect();
            writer.write_all(&chunk).unwrap();
            payload.extend_from_slice(&chunk);
            written += n;
            observed.push((writer.virtual_offset(), written));
        }
        writer.finish().unwrap();

        let block_starts = uncompressed_offsets_of_blocks(&output);
        for (voff, expected) in observed {
            let base = block_starts
                .get(&voff.block_offset())
                .copied()
                .unwrap_or_else(|| panic!("{voff:?} names no BGZF block in the output"));
            assert_eq!(
                base + usize::from(voff.within_block()),
                expected,
                "{voff:?} resolves to the wrong uncompressed position"
            );
        }
        assert_eq!(read_all(&output), payload);
    }

    /// Map each block's compressed file offset to its first uncompressed byte.
    fn uncompressed_offsets_of_blocks(data: &[u8]) -> std::collections::HashMap<u64, usize> {
        let mut map = std::collections::HashMap::new();
        let mut pos = 0usize;
        let mut uncompressed = 0usize;
        while pos + 18 <= data.len() {
            let bsize = u16::from_le_bytes([data[pos + 16], data[pos + 17]]) as usize + 1;
            let isize_val =
                u32::from_le_bytes(data[pos + bsize - 4..pos + bsize].try_into().unwrap()) as usize;
            map.insert(pos as u64, uncompressed);
            uncompressed += isize_val;
            pos += bsize;
        }
        map
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

    // r[verify bgzf.writer.block_size]
    /// The compile-time bound models libdeflate's own: at no level can the
    /// library's worst case for a full buffer exceed it, or overflow a block.
    #[test]
    fn libdeflates_bound_for_a_full_block_fits_at_every_level() {
        let modelled = HEADER_LEN + deflate_bound(MAX_UNCOMPRESSED_SIZE) + FOOTER_LEN;
        for level in 0..=12 {
            let lvl = libdeflater::CompressionLvl::new(level).unwrap();
            let actual = max_block_len(&mut libdeflater::Compressor::new(lvl));
            assert!(actual <= modelled, "level {level}: libdeflate bound {actual} > {modelled}");
            assert!(actual <= MAX_BLOCK_LEN, "level {level}: {actual}");
        }
    }

    // r[verify bgzf.writer.block_size]
    /// A block that cannot fit says so. 64 KiB stored at level 0 is two stored
    /// DEFLATE blocks (65535 + 1 bytes, 5 bytes of framing each): 65546 bytes
    /// of payload, 65572 with the gzip header and footer.
    #[test]
    fn a_block_that_cannot_fit_is_block_too_large() {
        let mut compressor =
            libdeflater::Compressor::new(libdeflater::CompressionLvl::new(0).unwrap());
        let data = vec![0u8; 65536];
        let mut block =
            vec![0; HEADER_LEN + compressor.deflate_compress_bound(data.len()) + FOOTER_LEN];
        let err = compress_block(&mut compressor, &data, &mut block).unwrap_err();
        assert!(matches!(err, BgzfError::BlockTooLarge { size: 65572 }), "{err:?}");
    }

    // r[verify bgzf.writer.block_size]
    /// Data that does not compress still makes valid blocks: level 0 stores
    /// every block, and noise does not shrink at level 6.
    #[test]
    fn incompressible_blocks_fit() {
        let mut x: u64 = 0x9e37_79b9_7f4a_7c15;
        let data: Vec<u8> = (0..300_000)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                (x >> 24) as u8
            })
            .collect();
        for level in [0, 6] {
            let mut out = Vec::new();
            let mut writer = BgzfWriter::with_compression_level(&mut out, level);
            writer.write_all(&data).unwrap();
            writer.finish().unwrap();
            assert_eq!(read_all(&out), data, "level {level}");
            for &start in uncompressed_offsets_of_blocks(&out).keys() {
                let bsize =
                    u16::from_le_bytes([out[start as usize + 16], out[start as usize + 17]]);
                assert!(usize::from(bsize) < 65536);
            }
        }
    }

    /// One step of a write stream: a record-sized write, preceded (as the BAM
    /// and BCF writers do) by `flush_if_needed` when `keep_whole`.
    #[derive(Debug, Clone)]
    struct Step {
        len: usize,
        fill: u8,
        keep_whole: bool,
    }

    #[hegel::composite]
    fn arb_step_inner(tc: &TestCase) -> Step {
        // Mostly record-sized, now and then more than a whole block.
        let len = if tc.draw(gs::integers::<u8>().max_value(9)) == 0 {
            tc.draw(gs::integers::<usize>().min_value(60_000).max_value(140_000))
        } else {
            tc.draw(gs::integers::<usize>().max_value(3_000))
        };
        Step { len, fill: tc.draw(gs::integers::<u8>()), keep_whole: tc.draw(gs::booleans()) }
    }

    fn arb_step() -> impl PrintableGenerator<Step> {
        arb_step_inner().print_as_debug()
    }

    /// Bytes for a step: half a repeated byte, half noise, so blocks neither
    /// vanish nor refuse to compress.
    fn step_bytes(step: &Step, seed: &mut u64) -> Vec<u8> {
        (0..step.len)
            .map(|i| {
                if i % 2 == 0 {
                    step.fill
                } else {
                    *seed ^= *seed << 13;
                    *seed ^= *seed >> 7;
                    *seed ^= *seed << 17;
                    (*seed >> 40) as u8 & 0x3f
                }
            })
            .collect()
    }

    // r[verify bgzf.writer.parallel]
    // r[verify bgzf.writer.parallel.identical_output]
    // r[verify bgzf.writer.parallel.index_offsets]
    /// The parallel writer produces the serial writer's bytes, and every index
    /// offset it handed out resolves to the virtual offset the serial writer
    /// reported at the same point of the stream.
    #[hegel::test(test_cases = 60)]
    fn parallel_writer_matches_serial(tc: TestCase) {
        let steps = tc.draw(gs::vecs(arb_step()).max_size(60));
        let threads = tc.draw(gs::integers::<usize>().min_value(1).max_value(4));
        let level = tc.draw(gs::sampled_from(&[0, 1, 6]));

        let mut serial = BgzfWriter::with_compression_level(Vec::new(), level);
        let mut parallel = ParallelBgzfWriter::new(Vec::new(), level, threads).unwrap();
        let mut expected = vec![serial.virtual_offset()];
        let mut logical = vec![parallel.index_offset()];
        let mut seed = 0x9e37_79b9_7f4a_7c15u64;
        for step in &steps {
            let bytes = step_bytes(step, &mut seed);
            if step.keep_whole {
                serial.flush_if_needed(bytes.len()).unwrap();
                parallel.flush_if_needed(bytes.len()).unwrap();
            }
            serial.write_all(&bytes).unwrap();
            parallel.write_all(&bytes).unwrap();
            expected.push(serial.virtual_offset());
            logical.push(parallel.index_offset());
        }

        parallel.flush_block().unwrap();
        parallel.pump(true).unwrap();
        let resolved: Vec<VirtualOffset> =
            logical.iter().map(|&v| parallel.resolve(v).expect("block written")).collect();
        assert_eq!(resolved, expected, "index offsets resolve to the serial virtual offsets");

        let serial_bytes = serial.finish().unwrap();
        let parallel_bytes = parallel.finish(None).unwrap();
        assert_eq!(parallel_bytes, serial_bytes, "parallel output is byte-identical");
    }

    // r[verify bgzf.writer.parallel]
    /// Dropping the parallel writer without `finish` still writes every block.
    #[test]
    fn parallel_writer_drop_flushes() {
        let data: Vec<u8> = (0..300_000).map(|i| (i % 251) as u8).collect();
        let mut serial = BgzfWriter::new(Vec::new());
        serial.write_all(&data).unwrap();
        let mut expected = serial.finish().unwrap();
        expected.truncate(expected.len() - EOF_BLOCK.len());

        let mut sink = Vec::new();
        {
            let mut parallel = ParallelBgzfWriter::new(&mut sink, 6, 3).unwrap();
            parallel.write_all(&data).unwrap();
        }
        assert_eq!(sink, expected);
    }

    /// A sink whose every write fails.
    struct BrokenSink;

    impl Write for BrokenSink {
        fn write(&mut self, _: &[u8]) -> std::io::Result<usize> {
            Err(std::io::Error::from(std::io::ErrorKind::BrokenPipe))
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    // r[verify bgzf.writer.parallel]
    /// A failed block stops the stream for good: later calls fail instead of
    /// waiting for a block that will never be written, and drop returns.
    #[test]
    fn parallel_writer_fails_fast_after_an_error() {
        let data = vec![7u8; 1_000_000];
        let mut parallel = ParallelBgzfWriter::new(BrokenSink, 6, 2).unwrap();
        let first = parallel.write_all(&data).and_then(|()| parallel.flush_block());
        let err = first.and_then(|()| parallel.pump(true));
        assert!(matches!(err, Err(BgzfError::WriteFailed { .. })), "{err:?}");
        assert!(parallel.write_all(&data).is_err(), "writes after a failure fail");
        drop(parallel);
    }
}
