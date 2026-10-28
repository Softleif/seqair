//! Front-to-back BGZF decoding from any [`Read`], with optional read-ahead
//! inflation on worker threads. See `r[bgzf.stream]`.

use super::{BgzfError, DataBlock, Frame, inflate_into, read_frame, resize_uninit};
use std::{
    collections::VecDeque,
    io::{self, BufRead, Read},
    sync::{Arc, Mutex, mpsc},
    thread,
};
use tracing::warn;

/// Blocks inflated per batch: about 512 KiB of output, enough that the channel
/// round trip is noise next to the inflation.
const BLOCKS_PER_BATCH: usize = 8;

/// What to do when the input ends without the BGZF EOF marker block.
///
/// See `r[bgzf.stream.eof_marker]`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum EofMarker {
    /// Fail with [`BgzfError::MissingEofMarker`]: the input was probably cut
    /// short. The default.
    #[default]
    Required,
    /// Accept the end of input at a block boundary as the end of the stream.
    Optional,
}

/// Configures a [`BgzfStream`].
#[derive(Debug)]
pub struct BgzfStreamBuilder<R> {
    inner: R,
    threads: usize,
    eof_marker: EofMarker,
}

impl<R: Read> BgzfStreamBuilder<R> {
    /// Inflate on `threads` worker threads; `0` (the default) inflates on the
    /// thread that reads.
    #[must_use]
    pub fn threads(mut self, threads: usize) -> Self {
        self.threads = threads;
        self
    }

    /// How to treat input that ends without an EOF marker block.
    #[must_use]
    pub fn eof_marker(mut self, eof_marker: EofMarker) -> Self {
        self.eof_marker = eof_marker;
        self
    }

    /// Start the stream (and its workers, if any). Reads nothing yet.
    ///
    /// # Errors
    /// [`BgzfError::ThreadSpawn`] if a worker thread cannot be started.
    pub fn build(self) -> Result<BgzfStream<R>, BgzfError> {
        let pool = match self.threads {
            0 => None,
            n => Some(Pool::start(n)?),
        };
        Ok(BgzfStream {
            inner: self.inner,
            current: Vec::new(),
            pos: 0,
            pending: VecDeque::new(),
            depth: self.threads.saturating_mul(2).max(1),
            input: Input::Open { last_was_empty: false },
            deferred: None,
            failed: false,
            eof_marker: self.eof_marker,
            decompressor: libdeflater::Decompressor::new(),
            spare: Vec::new(),
            pool,
        })
    }
}

/// How far the compressed input has been read.
#[derive(Debug, Clone, Copy)]
enum Input {
    /// More blocks may follow. `last_was_empty` tells whether the block read
    /// last was an empty one, which is what an EOF marker is.
    Open { last_was_empty: bool },
    /// Every block has been read (or reading failed; the failure is queued).
    Done,
}

/// What a batch yields, in stream order: the bytes of every block before the
/// first that failed, then that failure (or a read error that followed the
/// batch's blocks).
struct Batch {
    bytes: Vec<u8>,
    error: Option<BgzfError>,
}

impl Batch {
    fn failed(error: BgzfError) -> Self {
        Self { bytes: Vec::new(), error: Some(error) }
    }
}

/// A batch's output and the compressed buffer it came from, which goes back to
/// the stream for reuse.
struct Inflated {
    batch: Batch,
    compressed: Vec<u8>,
}

/// One batch, in stream order.
enum Pending {
    Ready(Batch),
    Waiting(mpsc::Receiver<Inflated>),
}

/// A batch of consecutive blocks for a worker.
struct Job {
    compressed: Vec<u8>,
    blocks: Vec<DataBlock>,
    out: Vec<u8>,
    reply: mpsc::SyncSender<Inflated>,
}

/// The inflation workers. Dropping closes the queue and joins them.
struct Pool {
    jobs: Option<mpsc::SyncSender<Job>>,
    workers: Vec<thread::JoinHandle<()>>,
}

impl Pool {
    fn start(threads: usize) -> Result<Self, BgzfError> {
        // The stream never has more than `2 * threads` batches out, so the
        // queue never makes it wait.
        let (jobs, queue) = mpsc::sync_channel::<Job>(threads.saturating_mul(2));
        let queue = Arc::new(Mutex::new(queue));
        let mut pool = Self { jobs: Some(jobs), workers: Vec::with_capacity(threads) };
        for i in 0..threads {
            let queue = Arc::clone(&queue);
            let worker = thread::Builder::new()
                .name(format!("seqair-bgzf-inflate-{i}"))
                .spawn(move || inflate_worker(&queue))
                // Dropping `pool` joins the workers already started.
                .map_err(|source| BgzfError::ThreadSpawn { source })?;
            pool.workers.push(worker);
        }
        Ok(pool)
    }
}

impl Drop for Pool {
    fn drop(&mut self) {
        // Closing the queue ends each worker after the batch it holds.
        drop(self.jobs.take());
        for worker in self.workers.drain(..) {
            if worker.join().is_err() {
                warn!("a BGZF inflation worker panicked");
            }
        }
    }
}

fn inflate_worker(queue: &Mutex<mpsc::Receiver<Job>>) {
    let mut decompressor = libdeflater::Decompressor::new();
    loop {
        // Hold the lock only while waiting for a job, not while inflating.
        let job = match queue.lock() {
            Ok(queue) => queue.recv(),
            Err(_) => return,
        };
        let Ok(Job { compressed, blocks, out, reply }) = job else { return };
        let batch = inflate_batch(&mut decompressor, &compressed, &blocks, out);
        // A failed send means the stream was dropped and wants nothing more.
        if reply.send(Inflated { batch, compressed }).is_err() {
            return;
        }
    }
}

/// Inflate `blocks`, in order, into `out` (replacing its contents). A block
/// that fails ends the batch: the blocks before it are kept.
fn inflate_batch(
    decompressor: &mut libdeflater::Decompressor,
    compressed: &[u8],
    blocks: &[DataBlock],
    mut out: Vec<u8>,
) -> Batch {
    let total = blocks.iter().map(|b| b.isize()).fold(0usize, usize::saturating_add);
    // Safety: every byte is written by `inflate_into` below, which covers
    // `0..total` block by block; on failure `out` is truncated to the blocks
    // written before return.
    unsafe { resize_uninit(&mut out, total) };
    let mut at = 0usize;
    for &block in blocks {
        let end = at.saturating_add(block.isize());
        let result = match out.get_mut(at..end) {
            Some(dst) => inflate_into(decompressor, compressed, block, dst),
            None => Err(BgzfError::TruncatedBlock),
        };
        if let Err(e) = result {
            out.truncate(at);
            return Batch { bytes: out, error: Some(e) };
        }
        at = end;
    }
    Batch { bytes: out, error: None }
}

// r[impl bgzf.stream]
/// A BGZF stream decoded front to back, from any [`Read`] — no `Seek` needed.
///
/// Built for whole-file reads (a BAM piped from an aligner, `stdin`):
/// [`BgzfReader`](super::BgzfReader) serves index-driven random access instead.
/// Empty blocks mid-stream are skipped (`r[bgzf.empty_block]`), input cut
/// short anywhere is an error (`r[bgzf.truncated_header]`,
/// `r[bgzf.stream.eof_marker]`), and after an error the stream stays failed
/// (`r[bgzf.stream.poison]`).
///
/// ```no_run
/// use seqair::bam::bgzf::{BgzfStream, EofMarker};
/// use std::io::BufRead;
///
/// let stdin = std::io::stdin().lock();
/// let mut stream = BgzfStream::builder(stdin).threads(4).build()?;
/// let mut line = Vec::new();
/// stream.read_until(b'\n', &mut line)?;
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub struct BgzfStream<R> {
    inner: R,
    /// Decompressed bytes of the batch being read, and the position in it.
    current: Vec<u8>,
    pos: usize,
    /// Batches read from the input but not yet returned, in stream order.
    pending: VecDeque<Pending>,
    /// Most batches `pending` may hold.
    depth: usize,
    input: Input,
    /// The error that ended the current batch, raised once its bytes are read.
    deferred: Option<BgzfError>,
    /// Set by the first error; every later read is [`BgzfError::Poisoned`].
    failed: bool,
    eof_marker: EofMarker,
    decompressor: libdeflater::Decompressor,
    /// Buffers to reuse for compressed batches.
    spare: Vec<Vec<u8>>,
    // Last, so the queue (and every reply channel) closes before the workers
    // are joined.
    pool: Option<Pool>,
}

impl<R> std::fmt::Debug for BgzfStream<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BgzfStream")
            .field("buffered", &self.current.len().saturating_sub(self.pos))
            .field("pending", &self.pending.len())
            .field("input", &self.input)
            .field("failed", &self.failed)
            .field("eof_marker", &self.eof_marker)
            .field("threads", &self.pool.as_ref().map_or(0, |p| p.workers.len()))
            .finish_non_exhaustive()
    }
}

impl<R: Read> BgzfStream<R> {
    /// A stream that inflates on the calling thread and requires the EOF marker.
    pub fn new(inner: R) -> Self {
        Self {
            inner,
            current: Vec::new(),
            pos: 0,
            pending: VecDeque::new(),
            depth: 1,
            input: Input::Open { last_was_empty: false },
            deferred: None,
            failed: false,
            eof_marker: EofMarker::Required,
            decompressor: libdeflater::Decompressor::new(),
            spare: Vec::new(),
            pool: None,
        }
    }

    /// Configure threads and EOF-marker handling.
    pub fn builder(inner: R) -> BgzfStreamBuilder<R> {
        BgzfStreamBuilder { inner, threads: 0, eof_marker: EofMarker::Required }
    }

    // r[impl bgzf.stream.api]
    /// The rest of the current block batch, refilling first if it is used up;
    /// empty only at the end of the stream. Pair with [`Self::advance`].
    ///
    /// # Errors
    /// Any [`BgzfError`] from reading or inflating; then
    /// [`BgzfError::Poisoned`] on every later call.
    pub fn fill(&mut self) -> Result<&[u8], BgzfError> {
        if self.pos >= self.current.len() && !self.refill()? {
            return Ok(&[]);
        }
        Ok(self.current.get(self.pos..).unwrap_or_default())
    }

    /// Mark `n` bytes of what [`Self::fill`] returned as read.
    pub fn advance(&mut self, n: usize) {
        self.pos = self.pos.saturating_add(n).min(self.current.len());
    }

    /// Read up to `out.len()` bytes; `0` only at the end of the stream.
    ///
    /// # Errors
    /// As [`Self::fill`].
    pub fn read_up_to(&mut self, out: &mut [u8]) -> Result<usize, BgzfError> {
        let available = self.fill()?;
        let n = available.len().min(out.len());
        if let (Some(dst), Some(src)) = (out.get_mut(..n), available.get(..n)) {
            dst.copy_from_slice(src);
        }
        self.advance(n);
        Ok(n)
    }

    /// Fill `out` exactly.
    ///
    /// # Errors
    /// [`BgzfError::UnexpectedEof`] if the stream ends first (what was read is
    /// lost), or as [`Self::fill`].
    pub fn read_exact_into(&mut self, out: &mut [u8]) -> Result<(), BgzfError> {
        let mut done = 0;
        while let Some(rest) = out.get_mut(done..)
            && !rest.is_empty()
        {
            match self.read_up_to(rest)? {
                0 => return Err(BgzfError::UnexpectedEof),
                n => done = done.saturating_add(n),
            }
        }
        Ok(())
    }

    /// Whether the stream has ended cleanly: everything was read and the
    /// input's end was accepted. Does not read.
    pub fn is_finished(&self) -> bool {
        !self.failed
            && self.deferred.is_none()
            && self.pos >= self.current.len()
            && self.pending.is_empty()
            && matches!(self.input, Input::Done)
    }

    /// Make the next non-empty batch current. `false` at the end.
    fn refill(&mut self) -> Result<bool, BgzfError> {
        if self.failed {
            return Err(BgzfError::Poisoned);
        }
        let result = self.refill_inner();
        if result.is_err() {
            self.failed = true;
            self.deferred = None;
            self.pending.clear();
            self.current.clear();
            self.pos = 0;
        }
        result
    }

    fn refill_inner(&mut self) -> Result<bool, BgzfError> {
        loop {
            if let Some(e) = self.deferred.take() {
                return Err(e);
            }
            self.top_up();
            let Some(next) = self.pending.pop_front() else { return Ok(false) };
            let Batch { bytes, error } = match next {
                Pending::Ready(batch) => batch,
                Pending::Waiting(reply) => {
                    let Inflated { batch, compressed } =
                        reply.recv().map_err(|_| BgzfError::DecompressionWorkerLost)?;
                    self.recycle(compressed);
                    batch
                }
            };
            self.deferred = error;
            // The batch we leave behind becomes the next one's output buffer.
            let done = std::mem::replace(&mut self.current, bytes);
            self.recycle(done);
            self.pos = 0;
            if !self.current.is_empty() {
                return Ok(true);
            }
        }
    }

    // r[impl bgzf.eof+2]
    // r[impl bgzf.empty_block]
    // r[impl bgzf.stream.eof_marker]
    // r[impl bgzf.stream.threads]
    /// Read batches from the input until `depth` are pending or the input is
    /// done. A read error, or a missing EOF marker, is queued after the batch
    /// it follows, so it reaches the consumer in stream order.
    fn top_up(&mut self) {
        while self.pending.len() < self.depth {
            let Input::Open { mut last_was_empty } = self.input else { return };
            let mut compressed = self.spare.pop().unwrap_or_default();
            compressed.clear();
            let mut blocks = Vec::with_capacity(BLOCKS_PER_BATCH);
            let mut tail = None;
            while blocks.len() < BLOCKS_PER_BATCH {
                match read_frame(&mut self.inner, &mut compressed, &mut self.decompressor) {
                    Ok(Frame::Data(block)) => {
                        blocks.push(block);
                        last_was_empty = false;
                    }
                    Ok(Frame::Empty) => last_was_empty = true,
                    Ok(Frame::EndOfInput) => {
                        self.input = Input::Done;
                        if !last_was_empty && self.eof_marker == EofMarker::Required {
                            tail = Some(BgzfError::MissingEofMarker);
                        }
                        break;
                    }
                    Err(e) => {
                        self.input = Input::Done;
                        tail = Some(e);
                        break;
                    }
                }
            }
            if let Input::Open { .. } = self.input {
                self.input = Input::Open { last_was_empty };
            }
            if !blocks.is_empty() {
                let batch = self.dispatch(compressed, blocks);
                self.pending.push_back(batch);
            }
            if let Some(e) = tail {
                self.pending.push_back(Pending::Ready(Batch::failed(e)));
            }
        }
    }

    /// Keep `buf` for reuse, up to what the batches in flight can need.
    fn recycle(&mut self, buf: Vec<u8>) {
        if self.spare.len() < self.depth.saturating_mul(2).saturating_add(2) {
            self.spare.push(buf);
        }
    }

    /// Inflate a batch now, or hand it to a worker.
    fn dispatch(&mut self, compressed: Vec<u8>, blocks: Vec<DataBlock>) -> Pending {
        let out = self.spare.pop().unwrap_or_default();
        let Some(jobs) = self.pool.as_ref().and_then(|p| p.jobs.as_ref()) else {
            let batch = inflate_batch(&mut self.decompressor, &compressed, &blocks, out);
            self.recycle(compressed);
            return Pending::Ready(batch);
        };
        let (reply, receiver) = mpsc::sync_channel(1);
        match jobs.send(Job { compressed, blocks, out, reply }) {
            Ok(()) => Pending::Waiting(receiver),
            // Every worker is gone.
            Err(_) => Pending::Ready(Batch::failed(BgzfError::DecompressionWorkerLost)),
        }
    }
}

/// `BgzfError` as an `io::Error` of kind `InvalidData`, recoverable by downcast.
fn to_io(e: BgzfError) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, e)
}

impl<R: Read> Read for BgzfStream<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        self.read_up_to(buf).map_err(to_io)
    }
}

impl<R: Read> BufRead for BgzfStream<R> {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        self.fill().map_err(to_io)
    }

    fn consume(&mut self, amt: usize) {
        self.advance(amt);
    }
}
