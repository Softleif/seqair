//! [`AlignmentStream`]: read BAM, SAM (plain, bgzipped or gzipped) or CRAM
//! front to back from any reader, deciding the format from its first bytes.
//! See `r[unified.stream.detect+4]`.

use crate::bam::{
    BamHeader, RawBamRecord,
    bgzf::{BgzfError, BgzfStream, EofMarker, Pool, Workers},
    header::BAM_MAGIC,
    stream::{BamStreamError, BamStreamReader},
};
use crate::cram::{
    slice::Plan,
    stream::{CramContainer, CramContainerDecoder, CramStreamError, CramStreamReader},
};
use crate::reader::{Complete, Partial, RecordFields, RecordKind};
use crate::sam::{SamStreamError, SamStreamReader};
use seqair_types::SmolStr;
use std::collections::VecDeque;
use std::io::{self, BufRead, BufReader, Cursor, Read};
use std::marker::PhantomData;
use std::panic::{self, AssertUnwindSafe};
use std::path::{Path, PathBuf};
use std::sync::{
    Arc,
    atomic::{AtomicBool, Ordering},
    mpsc,
};
use tracing::{trace, warn};

/// The format an [`AlignmentStream`] detected.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum StreamFormat {
    /// BGZF-compressed BAM.
    Bam,
    /// SAM text in BGZF blocks (`bgzip`).
    BgzfSam,
    /// SAM text in ordinary gzip members.
    GzipSam,
    /// Uncompressed SAM text (empty input included).
    Sam,
    /// CRAM v3.0 or v3.1.
    Cram,
}

/// Why an [`AlignmentStream`] could not be opened or read.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum AlignmentStreamError {
    #[error("reading the input failed")]
    Io { source: io::Error },

    #[error("BGZF error")]
    Bgzf {
        #[from]
        source: BgzfError,
    },

    #[error("BAM stream error")]
    Bam {
        #[from]
        source: BamStreamError,
    },

    #[error("SAM stream error")]
    Sam {
        #[from]
        source: SamStreamError,
    },

    #[error("CRAM stream error")]
    Cram {
        #[from]
        source: CramStreamError,
    },

    /// CRAM is decoded against a reference unless the fields read leave out
    /// the sequence and `MD`/`NM`; give one with
    /// [`AlignmentStreamBuilder::reference`].
    #[error(
        "the input is CRAM, whose sequence needs a FASTA reference (AlignmentStreamBuilder::reference)"
    )]
    CramNeedsReference,

    /// BAM in ordinary gzip rather than BGZF. htslib reads it; recompress with
    /// `samtools view -b` to stream it here.
    #[error("the input is BAM compressed with gzip rather than BGZF; recompress it with samtools")]
    GzipBamNotSupported,

    #[error("the input is uncompressed BAM, which is not supported")]
    UncompressedBamNotSupported,
}

type Source = std::io::Chain<Cursor<Vec<u8>>, Box<dyn Read + Send>>;

/// The input of a CRAM [`AlignmentStream`] hands over with
/// [`AlignmentStream::into_cram`]: the bytes it read to detect the format,
/// then the rest.
pub struct CramInput(Source);

impl Read for CramInput {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        self.0.read(buf)
    }
}

impl std::fmt::Debug for CramInput {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("CramInput").finish_non_exhaustive()
    }
}

enum Inner {
    Bam(Box<BamStreamReader<Source>>),
    Sam(Box<SamStreamReader<Box<dyn BufRead + Send>>>),
    Cram(Box<CramRecords>),
}

/// A CRAM stream: containers read on the reading thread, decoded there or,
/// with threads, on a rayon pool (`r[unified.stream.cram_threads]`).
struct CramRecords {
    reader: CramStreamReader<CramInput>,
    /// Decodes inline, and is what workers' decoders are forked from.
    decoder: CramContainerDecoder,
    /// The current container's records, `..len` of them valid, those before
    /// `next` handed out. Those past `len` are buffers kept for reuse.
    records: Vec<RawBamRecord>,
    len: usize,
    next: usize,
    /// `None` decodes on the reading thread.
    threads: Option<CramThreads>,
}

impl CramRecords {
    fn new(
        reader: CramStreamReader<CramInput>,
        decoder: CramContainerDecoder,
        workers: Workers,
    ) -> Self {
        let threads = workers.pool().map(|(pool, threads)| CramThreads::new(pool, threads));
        Self { reader, decoder, records: Vec::new(), len: 0, next: 0, threads }
    }

    fn read_record(&mut self, record: &mut RawBamRecord) -> Result<bool, CramStreamError> {
        loop {
            if self.next < self.len
                && let Some(slot) = self.records.get_mut(self.next)
            {
                // A swap: the caller's old record becomes a buffer the next
                // container is decoded into, and nothing is copied.
                std::mem::swap(record, slot);
                self.next = self.next.saturating_add(1);
                return Ok(true);
            }
            self.len = 0;
            self.next = 0;
            let Some(threads) = &mut self.threads else {
                let Some(container) = self.reader.next_container()? else { return Ok(false) };
                self.len = self.decoder.decode_reusing(&container, &mut self.records, 0)?;
                continue;
            };
            let Some(Decoded { records, result }) =
                threads.next(&mut self.reader, &mut self.decoder)
            else {
                return Ok(false);
            };
            let done = std::mem::replace(&mut self.records, records);
            threads.recycle(done);
            self.len = result?;
        }
    }

    /// Whether an error ended the stream and was returned.
    fn is_failed(&self) -> bool {
        self.reader.is_failed() && self.threads.as_ref().is_none_or(|t| t.pending.is_empty())
    }

    /// Whether every record read from the input has been handed out.
    fn is_drained(&self) -> bool {
        self.next >= self.len && self.threads.as_ref().is_none_or(|t| t.pending.is_empty())
    }
}

/// A container's records, `result` how many of them it holds (or why it
/// failed), and the buffer they are in.
struct Decoded {
    records: Vec<RawBamRecord>,
    result: Result<usize, CramStreamError>,
}

/// One container, in file order.
enum Pending {
    /// Read, waiting for a decoder.
    Queued(CramContainer),
    /// On a worker.
    Waiting {
        container: u64,
        reply: mpsc::Receiver<Decoded>,
    },
    Ready(Decoded),
}

/// Containers decoded in flight per thread asked for: one decoding, one read
/// and waiting for a decoder (or decoded and waiting for the consumer).
const CONTAINERS_PER_THREAD: usize = 2;

// r[impl unified.stream.cram_threads]
/// Containers decoded on a rayon pool by at most `threads` decoders, so at
/// most `threads` at once, with up to twice as many read ahead.
struct CramThreads {
    pool: Pool,
    threads: usize,
    pending: VecDeque<Pending>,
    /// Decoders forked so far and not lost to a panic.
    decoders: usize,
    /// Decoders free for the next container.
    idle: Vec<CramContainerDecoder>,
    /// Workers hand their decoder back here as soon as they are done, before
    /// their records, so a container decoded early frees its decoder at once.
    returned: mpsc::Receiver<CramContainerDecoder>,
    give_back: mpsc::Sender<CramContainerDecoder>,
    /// Record buffers to decode the next containers into.
    spare: Vec<Vec<RawBamRecord>>,
    /// Set when the stream is dropped, so queued jobs skip their container.
    dropped: Arc<AtomicBool>,
}

impl Drop for CramThreads {
    fn drop(&mut self) {
        self.dropped.store(true, Ordering::Relaxed);
    }
}

impl CramThreads {
    fn new(pool: Pool, threads: usize) -> Self {
        let (give_back, returned) = mpsc::channel();
        Self {
            pool,
            threads,
            pending: VecDeque::new(),
            decoders: 0,
            idle: Vec::new(),
            returned,
            give_back,
            spare: Vec::new(),
            dropped: Arc::new(AtomicBool::new(false)),
        }
    }

    /// The next container's records, in file order; `None` at the end.
    fn next(
        &mut self,
        reader: &mut CramStreamReader<CramInput>,
        template: &mut CramContainerDecoder,
    ) -> Option<Decoded> {
        self.top_up(reader);
        self.dispatch(template);
        Some(match self.pending.pop_front()? {
            Pending::Ready(decoded) => decoded,
            Pending::Waiting { container, reply } => match self.pool.wait(&reply) {
                Ok(decoded) => decoded,
                Err(_) => {
                    // Its decoder went down with it.
                    self.decoders = self.decoders.saturating_sub(1);
                    Decoded {
                        records: Vec::new(),
                        result: Err(CramStreamError::WorkerLost { container }),
                    }
                }
            },
            // Every earlier container is consumed, so every decoder is
            // free and `dispatch` gave it one, unless forking failed.
            Pending::Queued(container) => {
                let mut records = self.spare.pop().unwrap_or_default();
                let result = template.decode_reusing(&container, &mut records, 0);
                Decoded { records, result }
            }
        })
    }

    /// Read containers until `CONTAINERS_PER_THREAD` per thread are pending
    /// or the input is done, at most two per call so a slow producer (a pipe)
    /// is not waited on for the whole read-ahead before the first records
    /// return. A read error is queued after the containers before it, so it
    /// reaches the consumer in file order.
    fn top_up(&mut self, reader: &mut CramStreamReader<CramInput>) {
        let depth = self.threads.saturating_mul(CONTAINERS_PER_THREAD);
        for _ in 0..2 {
            // A failed reader's error is queued already; read again (for
            // `Poisoned`) only once it has been returned.
            if self.pending.len() >= depth || (reader.is_failed() && !self.pending.is_empty()) {
                return;
            }
            self.pending.push_back(match reader.next_container() {
                Ok(Some(container)) => Pending::Queued(container),
                // Ended: later calls to `next_container` keep saying so.
                Ok(None) => return,
                Err(e) => Pending::Ready(Decoded { records: Vec::new(), result: Err(e) }),
            });
        }
    }

    /// Hand queued containers, in file order, to workers while a decoder is
    /// free or another may be forked.
    fn dispatch(&mut self, template: &CramContainerDecoder) {
        self.idle.extend(self.returned.try_iter());
        for slot in &mut self.pending {
            let Pending::Queued(queued) = slot else { continue };
            let number = queued.number();
            let decoder = match self.idle.pop() {
                Some(decoder) => decoder,
                None if self.decoders < self.threads => match template.fork() {
                    Ok(decoder) => {
                        self.decoders = self.decoders.saturating_add(1);
                        decoder
                    }
                    Err(error) => {
                        warn!(%error, "cannot fork a CRAM decoder; decoding on fewer threads");
                        return;
                    }
                },
                None => return,
            };
            let (reply, receiver) = mpsc::sync_channel(1);
            let waiting = Pending::Waiting { container: number, reply: receiver };
            // Always `Queued`, matched above.
            if let Pending::Queued(container) = std::mem::replace(slot, waiting) {
                let records = self.spare.pop().unwrap_or_default();
                let give_back = self.give_back.clone();
                let dropped = Arc::clone(&self.dropped);
                self.pool.spawn(move || {
                    decode_job(decoder, &container, records, &give_back, &reply, &dropped);
                });
            }
        }
    }

    /// Keep `records` for reuse, up to what the containers in flight can need.
    fn recycle(&mut self, records: Vec<RawBamRecord>) {
        if self.spare.len() < self.threads.saturating_mul(CONTAINERS_PER_THREAD) {
            self.spare.push(records);
        }
    }
}

/// Decode a container on a pool thread, hand the decoder back, then send the
/// records. A panic (there should be none) drops `reply` unanswered, which
/// the stream reports as [`CramStreamError::WorkerLost`]; rayon would
/// otherwise abort.
fn decode_job(
    mut decoder: CramContainerDecoder,
    container: &CramContainer,
    mut records: Vec<RawBamRecord>,
    give_back: &mpsc::Sender<CramContainerDecoder>,
    reply: &mpsc::SyncSender<Decoded>,
    dropped: &AtomicBool,
) {
    if dropped.load(Ordering::Relaxed) {
        // The stream is gone; nobody will read these records.
        return;
    }
    let decoded = panic::catch_unwind(AssertUnwindSafe(move || {
        let result = decoder.decode_reusing(container, &mut records, 0);
        (decoder, Decoded { records, result })
    }));
    match decoded {
        Ok((decoder, decoded)) => {
            if give_back.send(decoder).is_err() || reply.send(decoded).is_err() {
                trace!("the CRAM stream was dropped while a container decoded");
            }
        }
        Err(_) => warn!(container = container.number(), "a CRAM decoding job panicked"),
    }
}

/// Configures an [`AlignmentStream`]. The thread options apply to the BGZF
/// formats and CRAM, the EOF-marker option to the BGZF formats.
#[derive(Debug)]
pub struct AlignmentStreamBuilder<R, K: RecordKind = Complete> {
    inner: R,
    workers: Workers,
    eof_marker: EofMarker,
    reference: Option<PathBuf>,
    name_prefix: Option<SmolStr>,
    fields: RecordFields,
    kind: PhantomData<K>,
}

impl<R: Read + Send + 'static> AlignmentStreamBuilder<R, Complete> {
    // r[impl unified.stream.fields]
    /// Read only `fields`: the stream yields
    /// [`PartialRecord`](crate::reader::PartialRecord)s, which refuse the
    /// rest, and may skip decoding it (`r[unified.record_fields]`).
    pub fn fields(self, fields: RecordFields) -> AlignmentStreamBuilder<R, Partial> {
        let Self { inner, workers, eof_marker, reference, name_prefix, .. } = self;
        AlignmentStreamBuilder {
            inner,
            workers,
            eof_marker,
            reference,
            name_prefix,
            fields,
            kind: PhantomData,
        }
    }
}

impl<R: Read + Send + 'static, K: RecordKind> AlignmentStreamBuilder<R, K> {
    /// Inflate up to `n` BGZF batches, or decode up to `n` CRAM containers,
    /// at once on the reading thread's rayon pool, else the global one
    /// (`r[bgzf.stream.current_pool]`, `r[unified.stream.cram_threads]`).
    #[must_use]
    pub fn threads(mut self, n: std::num::NonZeroU8) -> Self {
        self.workers = Workers::Threads { n, pool: None };
        self
    }

    /// Inflate up to `n` BGZF batches, or decode up to `n` CRAM containers,
    /// at once on `pool` (`r[bgzf.stream.shared_pool]`).
    #[must_use]
    pub fn threads_in(
        mut self,
        n: std::num::NonZeroU8,
        pool: impl Into<std::sync::Arc<rayon_core::ThreadPool>>,
    ) -> Self {
        self.workers = Workers::Threads { n, pool: Some(pool.into()) };
        self
    }

    /// The FASTA (with its `.fai`, and `.gzi` if bgzipped) CRAM input is
    /// decoded against. Other formats ignore it, and CRAM read with fields
    /// that leave out the sequence and `MD`/`NM` needs none
    /// (`r[cram.stream.fields.reference_optional]`).
    #[must_use]
    pub fn reference(mut self, fasta: &Path) -> Self {
        self.reference = Some(fasta.to_path_buf());
        self
    }

    /// The prefix of the names CRAM input without stored names (`lossy_names`)
    /// gets, `<prefix>:<record number>` (`r[cram.record.read_name+2]`);
    /// `-` by default, as htslib names stdin. [`AlignmentStream::open_path`]
    /// sets the file's basename, as `samtools view <file>` does.
    #[must_use]
    pub fn name_prefix(mut self, prefix: impl Into<SmolStr>) -> Self {
        self.name_prefix = Some(prefix.into());
        self
    }

    /// How BGZF input without an EOF marker is treated (`r[bgzf.stream.eof_marker]`).
    #[must_use]
    pub fn eof_marker(mut self, eof_marker: EofMarker) -> Self {
        self.eof_marker = eof_marker;
        self
    }

    /// Detect the format and read the header.
    ///
    /// # Errors
    /// Reading or decompression errors, a header error, or an unsupported
    /// format.
    pub fn build(self) -> Result<AlignmentStream<K>, AlignmentStreamError> {
        AlignmentStream::open(self)
    }
}

// r[impl unified.stream.records+2]
/// Any BAM, SAM or CRAM input, read front to back as [`RawBamRecord`]s, or
/// as [`PartialRecord`](crate::reader::PartialRecord)s when built with
/// [`AlignmentStreamBuilder::fields`].
///
/// ```no_run
/// use seqair::bam::RawBamRecord;
/// use seqair::reader::AlignmentStream;
///
/// let two = std::num::NonZeroU8::new(2).expect("2 is not 0");
/// let mut stream = AlignmentStream::builder(std::io::stdin()).threads(two).build()?;
/// let mut record = RawBamRecord::default();
/// while stream.read_record(&mut record)? {
///     // ...
/// }
/// // or, one allocation per record:
/// for record in AlignmentStream::from_reader(std::io::stdin())?.records() {
///     let record = record?;
/// }
/// # Ok::<(), seqair::reader::AlignmentStreamError>(())
/// ```
pub struct AlignmentStream<K: RecordKind = Complete> {
    inner: Inner,
    format: StreamFormat,
    fields: RecordFields,
    kind: PhantomData<K>,
}

impl<K: RecordKind> std::fmt::Debug for AlignmentStream<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("AlignmentStream").field("format", &self.format).finish_non_exhaustive()
    }
}

/// gzip magic, and the flag byte's FEXTRA bit.
const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];
const FEXTRA: u8 = 0x04;

impl AlignmentStream {
    /// Detect the format of `inner` with default options.
    ///
    /// # Errors
    /// As [`AlignmentStreamBuilder::build`].
    pub fn from_reader(inner: impl Read + Send + 'static) -> Result<Self, AlignmentStreamError> {
        Self::builder(inner).build()
    }

    /// Configure threads and EOF-marker handling before detecting the format.
    pub fn builder<R: Read + Send + 'static>(inner: R) -> AlignmentStreamBuilder<R> {
        AlignmentStreamBuilder {
            inner,
            workers: Workers::Inline,
            eof_marker: EofMarker::Required,
            reference: None,
            name_prefix: None,
            fields: RecordFields::all(),
            kind: PhantomData,
        }
    }

    /// A builder for the file at `path`, which names the records of CRAM
    /// input without stored names after its basename, as htslib does.
    ///
    /// # Errors
    /// The file cannot be opened.
    pub fn open_path(
        path: &Path,
    ) -> Result<AlignmentStreamBuilder<std::fs::File>, AlignmentStreamError> {
        let file =
            std::fs::File::open(path).map_err(|source| AlignmentStreamError::Io { source })?;
        let builder = Self::builder(file);
        Ok(match path.file_name().and_then(|n| n.to_str()) {
            Some(name) => builder.name_prefix(name),
            None => builder,
        })
    }
}

impl<K: RecordKind> AlignmentStream<K> {
    fn new(inner: Inner, format: StreamFormat, fields: RecordFields) -> Self {
        Self { inner, format, fields, kind: PhantomData }
    }

    // r[impl unified.stream.detect+4]
    fn open<R: Read + Send + 'static>(
        options: AlignmentStreamBuilder<R, K>,
    ) -> Result<Self, AlignmentStreamError> {
        let fields = options.fields;
        let mut inner: Box<dyn Read + Send> = Box::new(options.inner);
        let mut prefix = Vec::with_capacity(18);
        read_up_to(&mut inner, &mut prefix, 18)?;

        if prefix.starts_with(&GZIP_MAGIC) {
            if is_bgzf(&mut inner, &mut prefix)? {
                let source = Cursor::new(prefix).chain(inner);
                let mut bgzf = BgzfStream::builder(source)
                    .workers(options.workers)
                    .eof_marker(options.eof_marker)
                    .build();
                let magic = peek_magic(&mut bgzf)?;
                if magic == BAM_MAGIC {
                    let reader = BamStreamReader::after_magic(bgzf)?;
                    return Ok(Self::new(Inner::Bam(Box::new(reader)), StreamFormat::Bam, fields));
                }
                let text: Box<dyn BufRead + Send> = Box::new(Cursor::new(magic).chain(bgzf));
                return Self::sam(text, StreamFormat::BgzfSam, fields);
            }
            let source = Cursor::new(prefix).chain(inner);
            let mut gzip = BufReader::new(flate2::read::MultiGzDecoder::new(source));
            let mut magic = Vec::with_capacity(4);
            read_up_to(&mut gzip, &mut magic, 4)?;
            if magic == BAM_MAGIC {
                return Err(AlignmentStreamError::GzipBamNotSupported);
            }
            let text: Box<dyn BufRead + Send> = Box::new(Cursor::new(magic).chain(gzip));
            return Self::sam(text, StreamFormat::GzipSam, fields);
        }
        if prefix.starts_with(b"CRAM") {
            // r[impl cram.stream.fields.reference_optional]
            // Before the header: the decoder generates MD/NM by default.
            if options.reference.is_none() && Plan::new(&fields, true).needs_reference() {
                return Err(AlignmentStreamError::CramNeedsReference);
            }
            let reader = CramStreamReader::new(CramInput(Cursor::new(prefix).chain(inner)))?;
            let mut decoder = match options.reference {
                Some(fasta) => reader.decoder(&fasta)?,
                None => reader.decoder_without_reference(),
            }
            .with_fields(fields.clone());
            if let Some(prefix) = options.name_prefix {
                decoder = decoder.name_prefix(prefix);
            }
            let records = CramRecords::new(reader, decoder, options.workers);
            return Ok(Self::new(Inner::Cram(Box::new(records)), StreamFormat::Cram, fields));
        }
        if prefix.starts_with(&BAM_MAGIC) {
            return Err(AlignmentStreamError::UncompressedBamNotSupported);
        }
        let text: Box<dyn BufRead + Send> =
            Box::new(BufReader::new(Cursor::new(prefix).chain(inner)));
        Self::sam(text, StreamFormat::Sam, fields)
    }

    fn sam(
        text: Box<dyn BufRead + Send>,
        format: StreamFormat,
        fields: RecordFields,
    ) -> Result<Self, AlignmentStreamError> {
        let reader = SamStreamReader::new(text)?.with_aux_tags(fields.aux().clone());
        Ok(Self::new(Inner::Sam(Box::new(reader)), format, fields))
    }

    /// The fields read: [`RecordFields::all`] unless built with
    /// [`AlignmentStreamBuilder::fields`].
    pub fn fields(&self) -> &RecordFields {
        &self.fields
    }

    /// The format detected.
    pub fn format(&self) -> StreamFormat {
        self.format
    }

    /// The header, with its original text.
    pub fn header(&self) -> &BamHeader {
        match &self.inner {
            Inner::Bam(r) => r.header(),
            Inner::Sam(r) => r.header(),
            Inner::Cram(r) => r.reader.header(),
        }
    }

    /// Read the next record into `record`. `Ok(false)` at the end.
    ///
    /// # Errors
    /// The underlying reader's error; see [`BamStreamReader::read_record`] and
    /// [`SamStreamReader::read_record`] for which ones leave the stream usable.
    /// For CRAM, an error reading a container ends the stream
    /// ([`CramStreamReader::next_container`]), but a container that fails to
    /// decode (a reference sequence the FASTA lacks, a corrupt slice) is one
    /// error for the whole container: none of its records are returned, and
    /// the next call continues with the following container. Decoding on
    /// threads changes none of this; a worker that panics is one
    /// [`CramStreamError::WorkerLost`] for its container.
    pub fn read_record(&mut self, record: &mut K::Record) -> Result<bool, AlignmentStreamError> {
        let record = K::raw_mut(record, &self.fields);
        Ok(match &mut self.inner {
            Inner::Bam(r) => r.read_record(record)?,
            Inner::Sam(r) => r.read_record(record)?,
            Inner::Cram(r) => r.read_record(record)?,
        })
    }

    // r[impl unified.stream.into_cram+2]
    /// For CRAM, the reader and decoder this stream decodes with, set up as
    /// [`AlignmentStreamBuilder::build`] set them up (reference, name prefix),
    /// so the caller can decode containers on threads of its own
    /// ([`CramContainerDecoder::fork`]). The reader continues with the
    /// containers this stream has not read.
    ///
    /// # Errors
    /// The stream itself, unchanged, if it is not CRAM or holds records it
    /// has not handed out yet: with threads, any container read ahead, so a
    /// threaded stream hands over only before its first read or after its
    /// end.
    pub fn into_cram(self) -> Result<(CramStreamReader<CramInput>, CramContainerDecoder<K>), Self> {
        match self.inner {
            Inner::Cram(cram) if cram.is_drained() => {
                let CramRecords { reader, decoder, .. } = *cram;
                Ok((reader, decoder.retype()))
            }
            inner => Err(Self::new(inner, self.format, self.fields)),
        }
    }

    // r[impl unified.stream.iter]
    /// Iterate over the remaining records, one allocation each, as
    /// `Result<RawBamRecord, AlignmentStreamError>`. An error the stream can
    /// continue past is an item and iteration goes on: an invalid record, or
    /// for CRAM a container that fails to decode, which skips that
    /// container's records (see [`Self::read_record`]). An error that ends the
    /// stream is the last item. For a reused buffer, loop over
    /// [`Self::read_record`] instead.
    pub fn records(&mut self) -> Records<'_, K> {
        Records { stream: self, done: false }
    }

    /// Whether an error ended the stream.
    fn is_failed(&self) -> bool {
        match &self.inner {
            Inner::Bam(r) => r.is_failed(),
            Inner::Sam(r) => r.is_failed(),
            Inner::Cram(r) => r.is_failed(),
        }
    }
}

/// Iterator over an [`AlignmentStream`]'s records. See
/// [`AlignmentStream::records`].
#[derive(Debug)]
pub struct Records<'a, K: RecordKind = Complete> {
    stream: &'a mut AlignmentStream<K>,
    /// `None` was returned; the readers do not promise to keep returning it.
    done: bool,
}

impl<K: RecordKind> Iterator for Records<'_, K> {
    type Item = Result<K::Record, AlignmentStreamError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done || self.stream.is_failed() {
            return None;
        }
        let mut record = K::Record::default();
        match self.stream.read_record(&mut record) {
            Ok(true) => Some(Ok(record)),
            Ok(false) => {
                self.done = true;
                None
            }
            Err(e) => Some(Err(e)),
        }
    }
}

impl<K: RecordKind> std::iter::FusedIterator for Records<'_, K> {}

/// Append bytes from `inner` to `buf` until it holds `len` or the input ends.
fn read_up_to(
    inner: &mut impl Read,
    buf: &mut Vec<u8>,
    len: usize,
) -> Result<(), AlignmentStreamError> {
    let mut chunk = [0u8; 64];
    while buf.len() < len {
        let want = len.saturating_sub(buf.len()).min(chunk.len());
        match inner.read(chunk.get_mut(..want).unwrap_or_default()) {
            Ok(0) => break,
            Ok(n) => buf.extend_from_slice(chunk.get(..n).unwrap_or_default()),
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(source) => return Err(AlignmentStreamError::Io { source }),
        }
    }
    Ok(())
}

/// Whether the gzip header at the start of `prefix` is a BGZF block header
/// (`FEXTRA` with a `BC` subfield), reading the rest of a non-standard extra
/// field into `prefix` if needed.
fn is_bgzf(inner: &mut impl Read, prefix: &mut Vec<u8>) -> Result<bool, AlignmentStreamError> {
    let (Some(&flags), Some(&x0), Some(&x1)) = (prefix.get(3), prefix.get(10), prefix.get(11))
    else {
        return Ok(false);
    };
    if flags & FEXTRA == 0 {
        return Ok(false);
    }
    let extra_end = 12usize.saturating_add(usize::from(u16::from_le_bytes([x0, x1])));
    read_up_to(inner, prefix, extra_end)?;
    Ok(prefix.get(12..extra_end).and_then(crate::bam::bgzf::find_bsize).is_some())
}

/// The first (up to) four decompressed bytes, consumed from `bgzf`.
fn peek_magic<R: Read>(bgzf: &mut BgzfStream<R>) -> Result<Vec<u8>, BgzfError> {
    let mut magic = Vec::with_capacity(4);
    while magic.len() < 4 {
        let available = bgzf.fill()?;
        if available.is_empty() {
            break;
        }
        let n = available.len().min(4usize.saturating_sub(magic.len()));
        magic.extend_from_slice(available.get(..n).unwrap_or_default());
        bgzf.advance(n);
    }
    Ok(magic)
}
