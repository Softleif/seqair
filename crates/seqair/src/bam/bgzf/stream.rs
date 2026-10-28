//! Front-to-back BGZF decoding from any [`Read`], with optional read-ahead
//! inflation on worker threads. See `r[bgzf.stream]`.

use super::{BgzfError, DataBlock, Frame, inflate_into, read_frame, resize_uninit};
use rayon_core::ThreadPool;
use std::{
    cell::RefCell,
    collections::VecDeque,
    io::{self, BufRead, Read},
    num::NonZeroU8,
    panic::{self, AssertUnwindSafe},
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
        mpsc,
    },
};
use tracing::warn;

/// Blocks inflated per batch: about 1 MiB of output. A rayon `spawn` plus the
/// reply round trip costs about 5% of throughput at 8 blocks (M4, 4 threads,
/// a 51 MB BAM: 59.8 ms against 57.0 ms at 16), so batches are this large.
const BLOCKS_PER_BATCH: usize = 16;

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
    workers: Workers,
    eof_marker: EofMarker,
}

/// Where a stream's batches are inflated.
#[derive(Debug, Clone)]
pub(crate) enum Workers {
    /// On the thread that reads.
    Inline,
    /// On a rayon pool — the reading thread's (else the global one) when
    /// `pool` is `None` — with at most `n` batches in flight.
    Threads { n: NonZeroU8, pool: Option<Arc<ThreadPool>> },
}

impl<R: Read> BgzfStreamBuilder<R> {
    // r[impl bgzf.stream.current_pool]
    /// Inflate up to `n` batches at once on the rayon pool the reading thread
    /// runs in — the one a caller's `pool.install(..)` or parallel iterator
    /// put it in — or, read from any other thread, rayon's global pool. A
    /// rayon program shares its threads without naming its pool. Without this
    /// (or [`Self::threads_in`]) the stream inflates on the reading thread.
    #[must_use]
    pub fn threads(mut self, n: NonZeroU8) -> Self {
        self.workers = Workers::Threads { n, pool: None };
        self
    }

    // r[impl bgzf.stream.shared_pool]
    /// Inflate up to `n` batches at once on `pool`, a rayon pool shared with
    /// other work (rastair's pileup workers, say).
    #[must_use]
    pub fn threads_in(mut self, n: NonZeroU8, pool: impl Into<Arc<ThreadPool>>) -> Self {
        self.workers = Workers::Threads { n, pool: Some(pool.into()) };
        self
    }

    /// Where to inflate, as chosen on another builder that forwards here.
    #[must_use]
    pub(crate) fn workers(mut self, workers: Workers) -> Self {
        self.workers = workers;
        self
    }

    /// How to treat input that ends without an EOF marker block.
    #[must_use]
    pub fn eof_marker(mut self, eof_marker: EofMarker) -> Self {
        self.eof_marker = eof_marker;
        self
    }

    /// Configure the stream. Reads nothing yet, and starts no threads: a
    /// rayon pool starts its own (the global one on first use).
    pub fn build(self) -> BgzfStream<R> {
        let (pool, depth) = match self.workers {
            Workers::Inline => (None, 1),
            Workers::Threads { n, pool } => {
                let pool = pool.map_or(Pool::Current, Pool::Given);
                (
                    Some(pool),
                    std::num::NonZeroUsize::from(n).get().saturating_mul(BATCHES_PER_THREAD),
                )
            }
        };
        BgzfStream::with(self.inner, self.eof_marker, pool, depth)
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

/// The rayon pool a stream inflates on.
#[derive(Debug, Clone)]
enum Pool {
    /// One the stream owns or was given.
    Given(Arc<ThreadPool>),
    /// The reading thread's pool, else the global one, decided per batch.
    Current,
}

impl Pool {
    /// Whether the calling thread is one of this pool's.
    fn is_current(&self) -> bool {
        match self {
            Self::Given(pool) => pool.current_thread_index().is_some(),
            Self::Current => rayon_core::current_thread_index().is_some(),
        }
    }

    fn spawn(&self, job: impl FnOnce() + Send + 'static) {
        match self {
            Self::Given(pool) => pool.spawn(job),
            Self::Current => rayon_core::spawn(job),
        }
    }
}

/// Batches in flight per thread asked for: one inflating, one inflated and
/// waiting for the reader, so no thread waits on the reader to hand it work.
const BATCHES_PER_THREAD: usize = 2;

thread_local! {
    /// Each pool thread's decompressor, kept between batches.
    static DECOMPRESSOR: RefCell<Option<libdeflater::Decompressor>> = const { RefCell::new(None) };
}

/// Inflate a batch on a pool thread and send it back. A panic (there should
/// be none) drops `reply` unanswered, which the stream reports as
/// [`BgzfError::DecompressionWorkerLost`]; rayon would otherwise abort.
fn inflate_job(
    compressed: Vec<u8>,
    blocks: Vec<DataBlock>,
    out: Vec<u8>,
    reply: mpsc::SyncSender<Inflated>,
    dropped: &AtomicBool,
) {
    if dropped.load(Ordering::Relaxed) {
        // The stream is gone; nobody will read this batch.
        return;
    }
    let inflated = panic::catch_unwind(AssertUnwindSafe(|| {
        DECOMPRESSOR.with_borrow_mut(|d| {
            let d = d.get_or_insert_with(libdeflater::Decompressor::new);
            inflate_batch(d, &compressed, &blocks, out)
        })
    }));
    match inflated {
        // A failed send means the stream was dropped meanwhile.
        Ok(batch) => if reply.send(Inflated { batch, compressed }).is_err() {},
        Err(_) => warn!("a BGZF inflation job panicked"),
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
/// use seqair::bam::bgzf::BgzfStream;
/// use std::{io::BufRead, num::NonZeroU8};
///
/// let stdin = std::io::stdin().lock();
/// let four = NonZeroU8::new(4).expect("4 is not 0");
/// let mut stream = BgzfStream::builder(stdin).threads(four).build();
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
    /// Batches read so far.
    batches: u64,
    /// Buffers to reuse for compressed batches.
    spare: Vec<Vec<u8>>,
    /// Where batches are inflated; `None` inline.
    pool: Option<Pool>,
    /// Set when the stream is dropped, so queued jobs skip their batch.
    dropped: Arc<AtomicBool>,
}

impl<R> Drop for BgzfStream<R> {
    fn drop(&mut self) {
        self.dropped.store(true, Ordering::Relaxed);
    }
}

impl<R> std::fmt::Debug for BgzfStream<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BgzfStream")
            .field("buffered", &self.current.len().saturating_sub(self.pos))
            .field("pending", &self.pending.len())
            .field("input", &self.input)
            .field("failed", &self.failed)
            .field("eof_marker", &self.eof_marker)
            .field("in_flight", &self.depth)
            .finish_non_exhaustive()
    }
}

impl<R: Read> BgzfStream<R> {
    /// A stream that inflates on the calling thread and requires the EOF marker.
    pub fn new(inner: R) -> Self {
        Self::with(inner, EofMarker::Required, None, 1)
    }

    fn with(inner: R, eof_marker: EofMarker, pool: Option<Pool>, depth: usize) -> Self {
        Self {
            inner,
            current: Vec::new(),
            pos: 0,
            pending: VecDeque::new(),
            depth,
            input: Input::Open { last_was_empty: false },
            deferred: None,
            failed: false,
            eof_marker,
            decompressor: libdeflater::Decompressor::new(),
            batches: 0,
            spare: Vec::new(),
            pool,
            dropped: Arc::new(AtomicBool::new(false)),
        }
    }

    /// Configure threads (or a shared pool) and EOF-marker handling.
    pub fn builder(inner: R) -> BgzfStreamBuilder<R> {
        BgzfStreamBuilder { inner, workers: Workers::Inline, eof_marker: EofMarker::Required }
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
                    let Inflated { batch, compressed } = self.wait(&reply)?;
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

    // r[impl bgzf.stream.shared_pool]
    /// The reply for a batch. On one of the pool's own threads, blocking would
    /// take that thread from the pool the batch waits for (with one thread,
    /// forever), so it runs the pool's queued jobs instead until the reply is
    /// in, as rayon's own `join` does.
    fn wait(&self, reply: &mpsc::Receiver<Inflated>) -> Result<Inflated, BgzfError> {
        let lost = |_| BgzfError::DecompressionWorkerLost;
        let on_pool = self.pool.as_ref().is_some_and(Pool::is_current);
        if !on_pool {
            return reply.recv().map_err(lost);
        }
        loop {
            match reply.try_recv() {
                Ok(inflated) => return Ok(inflated),
                Err(mpsc::TryRecvError::Disconnected) => {
                    return Err(BgzfError::DecompressionWorkerLost);
                }
                Err(mpsc::TryRecvError::Empty) => match rayon_core::yield_now() {
                    Some(rayon_core::Yield::Executed) => {}
                    // Our batch is on another thread: wait for it there.
                    _ => return reply.recv().map_err(lost),
                },
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
        // At most two new batches per refill: the queue still grows to `depth`
        // one batch per batch consumed, but a slow producer (a pipe) is never
        // waited on for the whole read-ahead before the first bytes return.
        let mut added = 0u8;
        while self.pending.len() < self.depth && added < 2 {
            let Input::Open { mut last_was_empty } = self.input else { return };
            added = added.saturating_add(1);
            let mut compressed = self.spare.pop().unwrap_or_default();
            compressed.clear();
            // The first batch is one block, so the header (and the format a
            // caller sniffs from it) arrives as soon as that block has.
            let size = if self.batches == 0 { 1 } else { BLOCKS_PER_BATCH };
            self.batches = self.batches.saturating_add(1);
            let mut blocks = Vec::with_capacity(size);
            let mut tail = None;
            while blocks.len() < size {
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
        let Some(pool) = &self.pool else {
            let batch = inflate_batch(&mut self.decompressor, &compressed, &blocks, out);
            self.recycle(compressed);
            return Pending::Ready(batch);
        };
        let (reply, receiver) = mpsc::sync_channel(1);
        let dropped = Arc::clone(&self.dropped);
        pool.spawn(move || inflate_job(compressed, blocks, out, reply, &dropped));
        Pending::Waiting(receiver)
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
