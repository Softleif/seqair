//! [`Reader`]: records read front to back into a [`Fields`] type
//! (`r[vcf_reader.reader]`, `r[vcf_reader.read_into]`).

use super::error::{RecordAt, RecordError, VcfReadError};
use super::fields::{Fields, PlanError, ReadOptions, SampleSelection, Sections};
use super::file::{DataStart, FileShared, Input, VcfFile};
use super::raw::{BcfRecord, Fixed, RawRecord, Repr};
use super::site::SiteSpans;
use crate::bam::bgzf::{BgzfError, BgzfStream, Workers};
use crate::vcf::header::VcfHeader;
use rayon_core::ThreadPool;
use seqair_types::SmolStr;
use std::io::{Seek, SeekFrom};
use std::marker::PhantomData;
use std::num::NonZeroU8;
use std::sync::Arc;

/// Most bytes a BCF record's `l_shared + l_indiv` may claim
/// (`r[vcf_reader.limits]`).
pub const MAX_RECORD_BYTES: u64 = 256 << 20;

/// Configures a [`Reader`]; [`open`](Self::open) plans `F` and opens it.
#[must_use]
pub struct ReaderBuilder<'f, F> {
    file: &'f VcfFile,
    workers: Workers,
    sample: Option<SmolStr>,
    _f: PhantomData<fn() -> F>,
}

impl<F> std::fmt::Debug for ReaderBuilder<'_, F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ReaderBuilder")
            .field("workers", &self.workers)
            .field("sample", &self.sample)
            .finish_non_exhaustive()
    }
}

impl<'f, F: Fields> ReaderBuilder<'f, F> {
    pub(crate) fn new(file: &'f VcfFile) -> Self {
        Self { file, workers: Workers::Inline, sample: None, _f: PhantomData }
    }

    // r[impl vcf_reader.reader.threads]
    /// Decompress with up to `n` batches in flight on the rayon pool the
    /// reading thread runs in (else the global one).
    pub fn threads(mut self, n: NonZeroU8) -> Self {
        self.workers = Workers::Threads { n, pool: None };
        self
    }

    // r[impl vcf_reader.reader.threads]
    /// Decompress with up to `n` batches in flight on `pool`, so many
    /// readers can share a bounded set of threads.
    pub fn threads_in(mut self, n: NonZeroU8, pool: impl Into<Arc<ThreadPool>>) -> Self {
        self.workers = Workers::Threads { n, pool: Some(pool.into()) };
        self
    }

    /// The sample a [`Sample`](super::Sample) component reads.
    pub fn sample(mut self, name: impl Into<SmolStr>) -> Self {
        self.sample = Some(name.into());
        self
    }

    // r[impl vcf_fields.plan]
    /// Plan `F` against the header and open a reader at the first record:
    /// one file open, no header or index parse.
    ///
    /// # Errors
    /// [`VcfReadError::Plan`] if `F` cannot read this file, or an I/O error.
    pub fn open(self) -> Result<Reader<F>, VcfReadError> {
        let shared = &self.file.0;
        let opts = ReadOptions::new(self.sample);
        let plan = F::plan(&shared.header, &opts)?;
        if let Some(option) = opts.unused() {
            return Err(PlanError::UnusedOption { option }.into());
        }
        let sections = F::needs(&plan);
        let (start, skip) = match shared.data_start {
            DataStart::At(offset) => (offset.block_offset(), u64::from(offset.within_block())),
            DataStart::After(n) => (0, n),
        };
        let mut input = shared.open_input()?;
        input.seek(SeekFrom::Start(start)).map_err(|source| BgzfError::ReadFailed { source })?;
        // r[impl vcf_reader.reuse_io]
        // r[impl vcf_reader.fork_equivalence]
        let mut stream =
            BgzfStream::builder(input).workers(self.workers).start_offset(start).build();
        skip_bytes(&mut stream, skip).map_err(|e| e.into_read_error(RecordAt::default()))?;
        let current = F::empty(&plan);
        let core = Core {
            file: Arc::clone(shared),
            plan,
            sections,
            stream,
            scratch: Vec::new(),
            spans: SiteSpans::default(),
            ordinal: 0,
        };
        Ok(Reader { core, current })
    }
}

/// A failure while framing a record: a cut-short input is a record error.
enum Framing {
    Truncated,
    Bgzf(BgzfError),
}

impl From<BgzfError> for Framing {
    fn from(e: BgzfError) -> Self {
        match e {
            BgzfError::UnexpectedEof => Self::Truncated,
            e => Self::Bgzf(e),
        }
    }
}

impl Framing {
    fn into_read_error(self, at: RecordAt) -> VcfReadError {
        match self {
            Self::Truncated => at.error(RecordError::Truncated),
            Self::Bgzf(e) => e.into(),
        }
    }
}

fn skip_bytes(stream: &mut BgzfStream<Input>, mut n: u64) -> Result<(), Framing> {
    while n > 0 {
        let available = stream.fill()?.len();
        if available == 0 {
            return Err(Framing::Truncated);
        }
        let step = usize::try_from(n).map_or(available, |n| n.min(available));
        stream.advance(step);
        n = n.saturating_sub(step as u64);
    }
    Ok(())
}

/// A reader's I/O state and plan; everything but the record it lends out.
struct Core<P> {
    file: Arc<FileShared>,
    plan: P,
    sections: Sections,
    stream: BgzfStream<Input>,
    /// A record that spans the stream's batches, copied together.
    scratch: Vec<u8>,
    spans: SiteSpans,
    ordinal: u64,
}

impl<P> Core<P> {
    // r[impl vcf_reader.read_into]
    // r[impl vcf_fields.sections]
    /// Fill `rec` with the next record; `false` at the end.
    fn read_into<F: Fields<Plan = P>>(&mut self, rec: &mut F) -> Result<bool, VcfReadError> {
        if self.stream.fill()?.is_empty() {
            return Ok(false);
        }
        let at = RecordAt { ordinal: self.ordinal, offset: self.stream.virtual_offset() };
        self.ordinal = self.ordinal.saturating_add(1);
        let Self { file, plan, sections, stream, scratch, spans, .. } = self;

        // The record header straight from the batch, or, cut by the batch's
        // end, read (and kept) here.
        let (head, head_read) = match stream.fill()?.first_chunk::<8>() {
            Some(head) => (*head, false),
            None => {
                let mut head = [0u8; 8];
                stream
                    .read_exact_into(&mut head)
                    .map_err(|e| Framing::from(e).into_read_error(at))?;
                (head, true)
            }
        };
        let [s0, s1, s2, s3, i0, i1, i2, i3] = head;
        let l_shared = u32::from_le_bytes([s0, s1, s2, s3]);
        let l_indiv = u32::from_le_bytes([i0, i1, i2, i3]);
        // r[impl vcf_reader.limits]
        let size = u64::from(l_shared).saturating_add(u64::from(l_indiv));
        if size > MAX_RECORD_BYTES {
            return Err(at.error(RecordError::TooLarge { size, max: MAX_RECORD_BYTES }));
        }
        if l_shared < 24 {
            return Err(at.error(RecordError::SharedTooShort { len: l_shared }));
        }
        let shared_end = 8usize.saturating_add(l_shared as usize);
        let end = shared_end.saturating_add(l_indiv as usize);
        let wants_indiv = sections.samples() != SampleSelection::None;

        // The whole record in the current batch: read it in place.
        if !head_read && let Some(record) = stream.fill()?.get(..end) {
            let result = fill(file, plan, *sections, spans, at, record, shared_end, rec);
            stream.advance(end);
            return result.map(|()| true);
        }
        // Otherwise copy it together, skipping samples nobody reads.
        let framing = |e: Framing| e.into_read_error(at);
        scratch.clear();
        let mut copied = 0;
        if head_read {
            scratch.extend_from_slice(&head);
            copied = head.len();
        }
        scratch.resize(if wants_indiv { end } else { shared_end }, 0);
        stream
            .read_exact_into(scratch.get_mut(copied..).unwrap_or_default())
            .map_err(|e| framing(e.into()))?;
        if !wants_indiv {
            skip_bytes(stream, u64::from(l_indiv)).map_err(framing)?;
        }
        fill(file, plan, *sections, spans, at, scratch, shared_end, rec).map(|()| true)
    }
}

/// Fill `rec` from `record`, a BCF record whose shared section ends at
/// `shared_end` (its individual section, if read, follows).
#[expect(clippy::too_many_arguments, reason = "the reader's fields, split to borrow apart")]
fn fill<F: Fields>(
    file: &Arc<FileShared>,
    plan: &F::Plan,
    sections: Sections,
    spans: &mut SiteSpans,
    at: RecordAt,
    record: &[u8],
    shared_end: usize,
    rec: &mut F,
) -> Result<(), VcfReadError> {
    let shared = record.get(8..shared_end).unwrap_or_default();
    let indiv = record.get(shared_end..).unwrap_or_default();
    let fixed =
        Fixed::parse(shared).ok_or_else(|| at.error(RecordError::SharedTooShort { len: 0 }))?;
    // r[impl vcf_reader.contigs]
    let contig = file
        .contigs
        .by_rid(fixed.rid)
        .ok_or_else(|| at.error(RecordError::UnknownContig { index: fixed.rid }))?;
    if sections.site() || sections.info() {
        *spans = SiteSpans::find(shared, fixed.n_allele).map_err(|e| at.error(e))?;
    }
    let bcf = BcfRecord { bytes: record, shared, indiv, fixed, spans };
    let raw = RawRecord { file, contig, at, sections, repr: Repr::Bcf(bcf) };
    rec.fill(plan, raw)
}

/// Records read front to back into `F`.
///
/// Read with [`next_record`](Self::next_record) (lends one record at a time,
/// allocation-free) or [`read_into`](Self::read_into) (fills a record of
/// yours); [`records`](Self::records) and `into_iter()` hand out owned
/// records, allocating each.
pub struct Reader<F: Fields> {
    core: Core<F::Plan>,
    current: F,
}

impl<F: Fields> std::fmt::Debug for Reader<F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Reader")
            .field("records_read", &self.core.ordinal)
            .field("sections", &self.core.sections)
            .finish_non_exhaustive()
    }
}

impl<F: Fields> Reader<F> {
    // r[impl vcf_reader.passthrough_header]
    /// The file's header, as [`Writer::write_header`](crate::vcf::Writer::write_header) takes it.
    pub fn header(&self) -> &VcfHeader {
        &self.core.file.header
    }

    /// A record to [`read_into`](Self::read_into).
    pub fn empty_record(&self) -> F {
        F::empty(&self.core.plan)
    }

    // r[impl vcf_reader.read_into]
    /// Fill `rec` with the next record, reusing its buffers; `false` at the
    /// end of input.
    ///
    /// # Errors
    /// [`VcfReadError::Record`] for a malformed record (with its ordinal and
    /// virtual offset), or a decompression error.
    pub fn read_into(&mut self, rec: &mut F) -> Result<bool, VcfReadError> {
        self.core.read_into(rec)
    }

    // r[impl vcf_reader.next_record]
    /// The next record, lent until the next call; `None` at the end.
    ///
    /// # Errors
    /// As [`Self::read_into`].
    pub fn next_record(&mut self) -> Result<Option<&F>, VcfReadError> {
        Ok(self.core.read_into(&mut self.current)?.then_some(&self.current))
    }

    // r[impl vcf_reader.records]
    /// The records as owned values, for `peekable()`, channels and
    /// `collect`. Allocates per record; a hot loop wants
    /// [`next_record`](Self::next_record) or [`read_into`](Self::read_into).
    /// Ends after the first error.
    pub fn records(&mut self) -> Records<'_, F> {
        Records { reader: self, done: false }
    }
}

/// Owned records of a borrowed [`Reader`] ([`Reader::records`]).
#[derive(Debug)]
pub struct Records<'r, F: Fields> {
    reader: &'r mut Reader<F>,
    done: bool,
}

fn next_owned<F: Fields>(
    core: &mut Core<F::Plan>,
    done: &mut bool,
) -> Option<Result<F, VcfReadError>> {
    if *done {
        return None;
    }
    let mut rec = F::empty(&core.plan);
    match core.read_into(&mut rec) {
        Ok(true) => Some(Ok(rec)),
        Ok(false) => {
            *done = true;
            None
        }
        Err(e) => {
            *done = true;
            Some(Err(e))
        }
    }
}

impl<F: Fields> Iterator for Records<'_, F> {
    type Item = Result<F, VcfReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        next_owned(&mut self.reader.core, &mut self.done)
    }
}

/// Owned records of an owned [`Reader`] (`reader.into_iter()`), for
/// `par_bridge` and channels.
#[derive(Debug)]
pub struct IntoRecords<F: Fields> {
    reader: Reader<F>,
    done: bool,
}

impl<F: Fields> Iterator for IntoRecords<F> {
    type Item = Result<F, VcfReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        next_owned(&mut self.reader.core, &mut self.done)
    }
}

// r[impl vcf_reader.records]
impl<F: Fields> IntoIterator for Reader<F> {
    type Item = Result<F, VcfReadError>;
    type IntoIter = IntoRecords<F>;

    fn into_iter(self) -> IntoRecords<F> {
        IntoRecords { reader: self, done: false }
    }
}
