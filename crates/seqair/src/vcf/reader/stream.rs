//! [`Reader`]: records read front to back into a [`Fields`] type
//! (`r[vcf_reader.reader]`, `r[vcf_reader.read_into+2]`).

use super::entries::{FormatIndex, LazyIndex};
use super::error::{RecordAt, RecordError, RecordPart, VcfFormat, VcfReadError};
use super::fields::{
    Fields, PlanContext, PlanError, ReadOption, ReadOptions, SampleSel, SampleSelection, Sections,
};
use super::file::{DataStart, FileShared, Input, VcfFile};
use super::keys::Keys;
use super::query::{self, Partition, Query, RegionInput, Verdict};
use super::raw::{BcfRecord, Col, Fixed, RawRecord, Repr, TextCols, TextRecord, parse_pos};
use super::region::{Match, Region};
use super::rlen::{RlenColumns, text_rlen};
use super::site::{SiteParts, bcf_info_start};
use crate::bam::bgzf::{BgzfError, BgzfStream, VirtualOffset, Workers};
use crate::io::text_scan::find_byte;
use crate::vcf::header::VcfHeader;
use crate::vcf::record_encoder::ContigId;
use rayon_core::ThreadPool;
use std::io::{BufRead, Seek, SeekFrom};
use std::num::{NonZeroU8, NonZeroUsize};
use std::sync::Arc;

/// Most bytes a BCF record's `l_shared + l_indiv` may claim
/// (`r[vcf_reader.limits]`).
pub const MAX_RECORD_BYTES: u64 = 256 << 20;

/// Most bytes a VCF text line may hold, its newline not counted
/// (`r[vcf_reader.limits]`).
pub const MAX_LINE_BYTES: u64 = 256 << 20;

/// Configures a [`Reader`]; [`open`](Self::open) plans `F` and opens it.
#[must_use]
pub struct ReaderBuilder<'f, F: Fields> {
    file: &'f VcfFile,
    workers: Workers,
    sample: Option<SampleSel>,
    plan: Option<F::Plan>,
    matching: Option<Match>,
}

impl<F: Fields> std::fmt::Debug for ReaderBuilder<'_, F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ReaderBuilder")
            .field("workers", &self.workers)
            .field("sample", &self.sample)
            .field("with_plan", &self.plan.is_some())
            .field("matching", &self.matching)
            .finish_non_exhaustive()
    }
}

impl<'f, F: Fields> ReaderBuilder<'f, F> {
    pub(crate) fn new(file: &'f VcfFile) -> Self {
        Self { file, workers: Workers::Inline, sample: None, plan: None, matching: None }
    }

    // r[impl vcf_reader.reader.threads+2]
    /// Decompress with up to `n` (at most 255) batches in flight on the
    /// rayon pool the reading thread runs in (else the global one). Queries
    /// and partitions decompress on the calling thread.
    pub fn threads(mut self, n: NonZeroUsize) -> Self {
        self.workers = Workers::Threads { n: batches(n), pool: None };
        self
    }

    // r[impl vcf_reader.reader.threads+2]
    /// Decompress with up to `n` (at most 255) batches in flight on `pool`,
    /// so many readers can share a bounded set of threads.
    pub fn threads_in(mut self, n: NonZeroUsize, pool: impl Into<Arc<ThreadPool>>) -> Self {
        self.workers = Workers::Threads { n: batches(n), pool: Some(pool.into()) };
        self
    }

    // r[impl vcf_reader.query.overlap+2]
    /// How a [`query`](Self::query) selects records: [`Match::Overlapping`]
    /// (the default) or [`Match::StartsIn`].
    pub fn matching(mut self, matching: Match) -> Self {
        self.matching = Some(matching);
        self
    }

    // r[impl vcf_fields.sample+2]
    /// The sample a [`Sample`](super::Sample) component reads (and
    /// [`PlanContext::selected_sample`](super::PlanContext::selected_sample)
    /// gives a type of yours): a name, or a `SampleIdx` of the file's header.
    pub fn sample(mut self, sample: impl Into<SampleSel>) -> Self {
        self.sample = Some(sample.into());
        self
    }

    // r[impl vcf_fields.with_plan]
    /// Read with `plan` instead of planning `F` against the header: how a
    /// type of yours takes runtime context (a length band, a region), which
    /// lives in its `Plan` beside the keys you resolved with
    /// [`VcfFile::info_key`] and friends. Options the plan would have read
    /// ([`sample`](Self::sample)) then fail the open as unused.
    pub fn with_plan(mut self, plan: F::Plan) -> Self {
        self.plan = Some(plan);
        self
    }

    // r[impl vcf_fields.plan]
    /// Plan `F` against the header and open a reader at the first record:
    /// one file open, no header or index parse.
    ///
    /// # Errors
    /// [`VcfReadError::Plan`] if `F` cannot read this file (or
    /// [`matching`](Self::matching) was set), or an I/O error.
    pub fn open(self) -> Result<Reader<F>, VcfReadError> {
        self.start(Start::Sequential)
    }

    // r[impl vcf_reader.query+2]
    /// Plan `F` like [`open`](Self::open) and read only the records of
    /// `region` (a [`Region`], a `&`[`RegionString`](seqair_types::RegionString)),
    /// in file order.
    ///
    /// # Errors
    /// As [`open`](Self::open); [`VcfReadError::NoIndex`] without an index,
    /// [`VcfReadError::UnknownContig`] for a contig the file does not know,
    /// [`VcfReadError::ForeignKey`] for a contig of another file.
    pub fn query(self, region: impl Into<Region>) -> Result<Reader<F>, VcfReadError> {
        self.query_many([region])
    }

    // r[impl vcf_reader.query.many+2]
    /// Like [`query`](Self::query), over the union of `regions`: each record
    /// at most once, in file order.
    ///
    /// # Errors
    /// As [`query`](Self::query).
    pub fn query_many<I>(self, regions: I) -> Result<Reader<F>, VcfReadError>
    where
        I: IntoIterator,
        I::Item: Into<Region>,
    {
        let matching = self.matching.unwrap_or_default();
        self.start(Start::Query(regions.into_iter().map(Into::into).collect(), matching))
    }

    // r[impl vcf_reader.partitions.ownership+2]
    /// Read the records that start in `partition`
    /// ([`VcfFile::partitions`]): every record of the file is in exactly one.
    ///
    /// # Errors
    /// As [`query`](Self::query), and [`PlanError::UnusedOption`] if
    /// [`matching`](Self::matching) was set: a partition is always read
    /// [`Match::StartsIn`].
    pub fn partition(self, partition: &Partition) -> Result<Reader<F>, VcfReadError> {
        self.start(Start::Partition(partition.region().clone()))
    }

    pub(crate) fn start(self, start: Start) -> Result<Reader<F>, VcfReadError> {
        let shared = &self.file.0;
        let opts = ReadOptions::new(self.sample);
        let plan = match self.plan {
            Some(plan) => plan,
            None => F::plan(&PlanContext::new(shared, &opts))?,
        };
        if let Some(option) = opts.unused() {
            return Err(PlanError::UnusedOption { option }.into());
        }
        let sections = F::needs(&plan);
        // Line numbers are known reading VCF text from its start, not inside
        // a query's chunks.
        let mut line = shared.is_text().then_some(shared.header_lines);
        let (input, query) = match start {
            // r[impl vcf_reader.query.options]
            Start::Sequential | Start::Stream(_) | Start::Partition(_)
                if self.matching.is_some() =>
            {
                return Err(PlanError::UnusedOption { option: ReadOption::Matching }.into());
            }
            Start::Sequential => (open_sequential(shared, self.workers)?, None),
            Start::Stream(input) => (input, None),
            Start::Query(regions, matching) => {
                line = None;
                query_input(shared, regions, matching)?
            }
            Start::Partition(region) => {
                line = None;
                query_input(shared, vec![region], Match::StartsIn)?
            }
        };
        let core = Core {
            file: Arc::clone(shared),
            plan,
            sections,
            input,
            query,
            scratch: Vec::new(),
            site: SiteParts::default(),
            info_index: LazyIndex::default(),
            format_index: LazyIndex::default(),
            text: TextState::default(),
            ordinal: 0,
            line,
        };
        Ok(Reader { core, current: None })
    }
}

/// Where a reader starts.
pub(crate) enum Start {
    Sequential,
    /// A [`VcfStream`](super::VcfStream)'s input, past its header.
    Stream(Framed),
    Query(Vec<Region>, Match),
    Partition(Region),
}

/// A query's input and filter.
fn query_input(
    shared: &FileShared,
    regions: Vec<Region>,
    matching: Match,
) -> Result<(Framed, Option<Query>), VcfReadError> {
    let (query, chunks) = query::plan(shared, regions, matching)?;
    let input = RegionInput::new(shared.open_input()?, chunks)?;
    Ok((Framed::Region(Box::new(input)), Some(query)))
}

/// `n` batches in flight, as many as a `BgzfStream` takes.
pub(crate) fn batches(n: NonZeroUsize) -> NonZeroU8 {
    NonZeroU8::try_from(n).unwrap_or(NonZeroU8::MAX)
}

/// A file's records from the first, front to back.
fn open_sequential(shared: &FileShared, workers: Workers) -> Result<Framed, VcfReadError> {
    let (start, skip) = match shared.data_start {
        DataStart::At(offset) => (offset.block_offset(), u64::from(offset.within_block())),
        DataStart::After(n) => (0, n),
    };
    let mut input = shared.open_input()?;
    let seek = |input: &mut Input, at: u64| {
        input.seek(SeekFrom::Start(at)).map_err(|source| BgzfError::ReadFailed { source })
    };
    if shared.format == VcfFormat::Vcf {
        seek(&mut input, start.saturating_add(skip))?;
        return Ok(Framed::Plain(input));
    }
    seek(&mut input, start)?;
    // r[impl vcf_reader.reuse_io+2]
    // r[impl vcf_reader.fork_equivalence]
    let stream = BgzfStream::builder(input).workers(workers).start_offset(start).build();
    let mut input = Framed::Bgzf(Box::new(stream));
    skip_bytes(&mut input, skip).map_err(|e| e.into_read_error(RecordAt::default(), shared))?;
    Ok(input)
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
    fn into_read_error(self, at: RecordAt, file: &FileShared) -> VcfReadError {
        match self {
            Self::Truncated => at.error(file, RecordError::Truncated),
            Self::Bgzf(e) => e.into(),
        }
    }
}

/// A reader's decompressed input: a BGZF stream, a plain text file, or a
/// query's chunks.
pub(crate) enum Framed {
    Bgzf(Box<BgzfStream<Input>>),
    Plain(Input),
    Region(Box<RegionInput>),
}

impl Framed {
    fn fill(&mut self) -> Result<&[u8], BgzfError> {
        match self {
            Self::Bgzf(s) => s.fill(),
            Self::Plain(p) => p.fill_buf().map_err(|source| BgzfError::ReadFailed { source }),
            Self::Region(r) => r.fill(),
        }
    }

    fn advance(&mut self, n: usize) {
        match self {
            Self::Bgzf(s) => s.advance(n),
            Self::Plain(p) => p.consume(n),
            Self::Region(r) => r.advance(n),
        }
    }

    /// Whether another record follows: input is left, or, for a query,
    /// left in its chunks.
    fn has_record(&mut self) -> Result<bool, BgzfError> {
        if let Self::Region(r) = self {
            return r.has_record();
        }
        Ok(!self.fill()?.is_empty())
    }

    /// The query group of the chunk being read.
    fn group(&self) -> Option<usize> {
        match self {
            Self::Region(r) => r.group(),
            _ => None,
        }
    }

    fn read_exact_into(&mut self, out: &mut [u8]) -> Result<(), BgzfError> {
        let mut done = 0usize;
        while let Some(rest) = out.get_mut(done..)
            && !rest.is_empty()
        {
            let available = self.fill()?;
            let n = available.len().min(rest.len());
            if n == 0 {
                return Err(BgzfError::UnexpectedEof);
            }
            if let (Some(dst), Some(src)) = (rest.get_mut(..n), available.get(..n)) {
                dst.copy_from_slice(src);
            }
            self.advance(n);
            done = done.saturating_add(n);
        }
        Ok(())
    }

    /// The virtual offset of the next byte; `None` for plain text.
    #[inline]
    fn virtual_offset(&mut self) -> Option<VirtualOffset> {
        match self {
            Self::Bgzf(s) => s.next_virtual_offset(),
            Self::Plain(_) => None,
            Self::Region(r) => Some(r.virtual_offset()),
        }
    }
}

fn skip_bytes(input: &mut Framed, mut n: u64) -> Result<(), Framing> {
    while n > 0 {
        let available = input.fill()?.len();
        if available == 0 {
            return Err(Framing::Truncated);
        }
        let step = usize::try_from(n).map_or(available, |n| n.min(available));
        input.advance(step);
        n = n.saturating_sub(step as u64);
    }
    Ok(())
}

/// What a text reader keeps between lines: the columns of the current one,
/// and the last CHROM with its contig, so a contig is looked up only when
/// CHROM changes (`r[vcf_reader.contigs]`).
#[derive(Debug, Default)]
struct TextState {
    cols: TextCols,
    chrom: Vec<u8>,
    contig: Option<ContigId>,
}

impl TextState {
    /// Find `line`'s columns, its contig and its position (into `at`).
    fn read(
        &mut self,
        line: &[u8],
        file: &FileShared,
        at: &mut RecordAt,
    ) -> Result<(), RecordError> {
        self.cols.find(line)?;
        let chrom = self.cols.col(line, Col::Chrom);
        if self.contig.is_none() || self.chrom != chrom {
            let name = std::str::from_utf8(chrom)
                .map_err(|source| RecordError::Utf8 { part: RecordPart::Chrom, source })?;
            // r[impl vcf_reader.contigs]
            self.contig = Some(file.contigs.by_name(&file.header, name));
            self.chrom.clear();
            self.chrom.extend_from_slice(chrom);
        }
        at.contig = self.contig.as_ref().map(ContigId::tid);
        at.pos = parse_pos(self.cols.col(line, Col::Pos)).map(|pos1| pos1.saturating_sub(1));
        Ok(())
    }
}

/// A reader's I/O state and plan; everything but the record it lends out.
struct Core<P> {
    file: Arc<FileShared>,
    plan: P,
    sections: Sections,
    input: Framed,
    /// What a query yields of the records it reads; `None` reading front to back.
    query: Option<Query>,
    /// A record that spans the stream's batches, copied together.
    scratch: Vec<u8>,
    /// The current record's site columns, decoded.
    site: SiteParts,
    /// The current record's INFO and FORMAT indexes, for [`RawRecord`]'s views.
    info_index: LazyIndex,
    format_index: LazyIndex<FormatIndex>,
    text: TextState,
    ordinal: u64,
    /// The line last framed, for VCF text read from its start.
    line: Option<u64>,
}

/// What reading one record came to.
enum Step<R> {
    Took(R),
    /// Not in the query.
    Skip,
    /// Past the query's ranges on the group's contig.
    Stop(usize),
    End,
}

/// What a query does with a record: `None` to go on decoding it.
fn query_step<R>(
    query: Option<&Query>,
    group: Option<usize>,
    tid: u32,
    pos0: i64,
    rlen: impl FnOnce() -> i64,
) -> Option<Step<R>> {
    let (query, group) = (query?, group?);
    match query.verdict(group, tid, pos0, rlen) {
        Verdict::Take => None,
        Verdict::Skip => Some(Step::Skip),
        Verdict::Stop => Some(Step::Stop(group)),
    }
}

/// Hand `raw` to the caller's `with`, once.
fn take<P, R, W>(
    with: &mut Option<W>,
    plan: &P,
    raw: RawRecord<'_>,
) -> Result<Step<R>, VcfReadError>
where
    W: FnOnce(&P, RawRecord<'_>) -> Result<R, VcfReadError>,
{
    match with.take() {
        Some(with) => with(plan, raw).map(Step::Took),
        None => Ok(Step::Skip),
    }
}

/// Where a line is: in the input's current batch (its length, and the
/// bytes to step past it), or copied into the scratch buffer.
enum Line {
    InPlace { len: usize, consumed: usize },
    Copied,
}

impl<P> Core<P> {
    // r[impl vcf_reader.read_into+2]
    // r[impl vcf_fields.sections]
    /// Hand the next record to `with`; `None` at the end.
    fn next<R>(
        &mut self,
        with: impl FnOnce(&P, RawRecord<'_>) -> Result<R, VcfReadError>,
    ) -> Result<Option<R>, VcfReadError> {
        let mut with = Some(with);
        loop {
            self.info_index.reset();
            self.format_index.reset();
            let step = if self.file.is_text() {
                self.next_text(&mut with)?
            } else {
                self.next_bcf(&mut with)?
            };
            match step {
                Step::Took(r) => return Ok(Some(r)),
                Step::End => return Ok(None),
                Step::Skip => {}
                Step::Stop(group) => {
                    if let Framed::Region(r) = &mut self.input {
                        r.skip_group(group);
                    }
                }
            }
        }
    }

    /// Frame the next non-empty line; `None` at the end.
    fn frame_line(&mut self) -> Result<Option<(Line, RecordAt)>, VcfReadError> {
        loop {
            if !self.input.has_record()? {
                return Ok(None);
            }
            self.line = self.line.map(|l| l.saturating_add(1));
            let at = RecordAt {
                ordinal: self.ordinal,
                line: self.line,
                offset: self.input.virtual_offset(),
                ..RecordAt::default()
            };
            let buf = self.input.fill()?;
            let line = match find_byte(buf, b'\n') {
                Some(len) if len > 0 && buf.get(..len) != Some(b"\r") => {
                    Line::InPlace { len, consumed: len.saturating_add(1) }
                }
                Some(len) => {
                    // An empty line holds no record.
                    self.input.advance(len.saturating_add(1));
                    continue;
                }
                None => {
                    self.copy_line(at)?;
                    if self.scratch.is_empty() || self.scratch == b"\r" {
                        continue;
                    }
                    Line::Copied
                }
            };
            self.ordinal = self.ordinal.saturating_add(1);
            return Ok(Some((line, at)));
        }
    }

    // r[impl vcf_reader.limits]
    /// Copy a line that runs past the current batch into `scratch`, up to
    /// [`MAX_LINE_BYTES`]; the input's end ends a last line too.
    fn copy_line(&mut self, at: RecordAt) -> Result<(), VcfReadError> {
        self.scratch.clear();
        loop {
            let buf = self.input.fill()?;
            if buf.is_empty() {
                return Ok(());
            }
            let (take, done) = match find_byte(buf, b'\n') {
                Some(i) => (i, true),
                None => (buf.len(), false),
            };
            let len = self.scratch.len().saturating_add(take) as u64;
            if len > MAX_LINE_BYTES {
                let e = RecordError::LineTooLong { len, max: MAX_LINE_BYTES };
                return Err(at.error(&self.file, e));
            }
            self.scratch.extend_from_slice(buf.get(..take).unwrap_or_default());
            self.input.advance(if done { take.saturating_add(1) } else { take });
            if done {
                return Ok(());
            }
        }
    }

    fn next_text<R, W>(&mut self, with: &mut Option<W>) -> Result<Step<R>, VcfReadError>
    where
        W: FnOnce(&P, RawRecord<'_>) -> Result<R, VcfReadError>,
    {
        let Some((framed, mut at)) = self.frame_line()? else { return Ok(Step::End) };
        let group = self.input.group();
        let Self {
            file,
            plan,
            sections,
            input,
            query,
            scratch,
            site,
            info_index,
            format_index,
            text,
            ..
        } = self;
        let (line, consumed) = match framed {
            Line::InPlace { len, consumed } => {
                (input.fill()?.get(..len).unwrap_or_default(), Some(consumed))
            }
            Line::Copied => (scratch.as_slice(), None),
        };
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        let result = (|| {
            text.read(line, file, &mut at).map_err(|e| at.error(file, e))?;
            if let (Some(q), Some(contig)) = (query.as_ref(), text.contig.as_ref()) {
                let cols = &text.cols;
                let pos0 = at.pos.ok_or_else(|| {
                    at.error(
                        file,
                        RecordError::Text {
                            part: RecordPart::Pos,
                            source: super::text_value::TextValueError::invalid_integer(
                                cols.col(line, Col::Pos),
                            ),
                        },
                    )
                })?;
                let rlen = || {
                    let columns = RlenColumns {
                        reference: cols.col(line, Col::Ref),
                        alt: cols.col(line, Col::Alt),
                        info: cols.col(line, Col::Info),
                        format: cols.rest(line),
                    };
                    text_rlen(&columns, file, pos0)
                };
                if let Some(step) = query_step(Some(q), group, contig.tid(), pos0, rlen) {
                    return Ok(step);
                }
            }
            if sections.site() {
                site.decode_text(line, &text.cols, file).map_err(|e| at.error(file, e))?;
            }
            let contig =
                text.contig.as_ref().ok_or_else(|| at.error(file, RecordError::Truncated))?;
            let raw = RawRecord {
                file,
                contig,
                at,
                sections: *sections,
                site: sections.site().then_some(&*site),
                info_index,
                format_index,
                repr: Repr::Text(TextRecord { line, cols: &text.cols }),
            };
            take(with, plan, raw)
        })();
        if let Some(n) = consumed {
            input.advance(n);
        }
        result
    }

    fn next_bcf<R, W>(&mut self, with: &mut Option<W>) -> Result<Step<R>, VcfReadError>
    where
        W: FnOnce(&P, RawRecord<'_>) -> Result<R, VcfReadError>,
    {
        if !self.input.has_record()? {
            return Ok(Step::End);
        }
        let at = RecordAt {
            ordinal: self.ordinal,
            offset: self.input.virtual_offset(),
            ..RecordAt::default()
        };
        self.ordinal = self.ordinal.saturating_add(1);
        let group = self.input.group();
        let Self {
            file,
            plan,
            sections,
            input,
            query,
            scratch,
            site,
            info_index,
            format_index,
            ..
        } = self;
        let parts = Parts {
            file,
            sections: *sections,
            site,
            info_index,
            format_index,
            at,
            query: query.as_ref(),
            group,
        };

        // The record header straight from the batch, or, cut by the batch's
        // end, read (and kept) here.
        let (head, head_read) = match input.fill()?.first_chunk::<8>() {
            Some(head) => (*head, false),
            None => {
                let mut head = [0u8; 8];
                input
                    .read_exact_into(&mut head)
                    .map_err(|e| Framing::from(e).into_read_error(at, file))?;
                (head, true)
            }
        };
        let [s0, s1, s2, s3, i0, i1, i2, i3] = head;
        let l_shared = u32::from_le_bytes([s0, s1, s2, s3]);
        let l_indiv = u32::from_le_bytes([i0, i1, i2, i3]);
        // r[impl vcf_reader.limits]
        let size = u64::from(l_shared).saturating_add(u64::from(l_indiv));
        if size > MAX_RECORD_BYTES {
            return Err(at.error(file, RecordError::TooLarge { size, max: MAX_RECORD_BYTES }));
        }
        if l_shared < 24 {
            return Err(at.error(file, RecordError::SharedTooShort { len: l_shared }));
        }
        let shared_end = 8usize.saturating_add(l_shared as usize);
        let end = shared_end.saturating_add(l_indiv as usize);
        let wants_indiv = sections.samples() != SampleSelection::None;

        // The whole record in the current batch: read it in place.
        if !head_read && let Some(record) = input.fill()?.get(..end) {
            let result = parts.hand_out(plan, record, shared_end, with);
            input.advance(end);
            return result;
        }
        // Otherwise copy it together, skipping samples nobody reads.
        let framing = |e: Framing| e.into_read_error(at, file);
        scratch.clear();
        let mut copied = 0;
        if head_read {
            scratch.extend_from_slice(&head);
            copied = head.len();
        }
        scratch.resize(if wants_indiv { end } else { shared_end }, 0);
        input
            .read_exact_into(scratch.get_mut(copied..).unwrap_or_default())
            .map_err(|e| framing(e.into()))?;
        if !wants_indiv {
            skip_bytes(input, u64::from(l_indiv)).map_err(framing)?;
        }
        parts.hand_out(plan, scratch, shared_end, with)
    }
}

/// The parts of a [`Core`] a BCF record is decoded with, borrowed apart
/// from its input.
struct Parts<'c> {
    file: &'c Arc<FileShared>,
    sections: Sections,
    site: &'c mut SiteParts,
    info_index: &'c LazyIndex,
    format_index: &'c LazyIndex<FormatIndex>,
    at: RecordAt,
    query: Option<&'c Query>,
    group: Option<usize>,
}

impl Parts<'_> {
    /// Hand `with` the BCF record `record`, whose shared section ends at
    /// `shared_end` (its individual section, if read, follows).
    fn hand_out<P, R, W>(
        self,
        plan: &P,
        record: &[u8],
        shared_end: usize,
        with: &mut Option<W>,
    ) -> Result<Step<R>, VcfReadError>
    where
        W: FnOnce(&P, RawRecord<'_>) -> Result<R, VcfReadError>,
    {
        let Self { file, sections, site, info_index, format_index, mut at, query, group } = self;
        let shared = record.get(8..shared_end).unwrap_or_default();
        let indiv = record.get(shared_end..).unwrap_or_default();
        let fixed = Fixed::parse(shared)
            .ok_or_else(|| at.error(file, RecordError::SharedTooShort { len: 0 }))?;
        at.pos = Some(i64::from(fixed.pos));
        // r[impl vcf_reader.contigs]
        let contig = file.contigs.by_rid(fixed.rid).ok_or_else(|| {
            at.error(file, RecordError::ContigIndexOutOfRange { index: fixed.rid })
        })?;
        at.contig = Some(contig.tid());
        let rlen = || i64::from(fixed.rlen);
        if let Some(step) = query_step(query, group, contig.tid(), i64::from(fixed.pos), rlen) {
            return Ok(step);
        }
        let info_start = if sections.site() {
            site.decode_bcf(shared, fixed.n_allele, fixed.qual, file)
        } else if sections.info() {
            bcf_info_start(shared, fixed.n_allele)
        } else {
            Ok(shared.len())
        };
        let info_start = info_start.map_err(|e| at.error(file, e))?;
        let bcf = BcfRecord { bytes: record, shared, indiv, fixed, info_start };
        let site = sections.site().then_some(&*site);
        let raw = RawRecord {
            file,
            contig,
            at,
            sections,
            site,
            info_index,
            format_index,
            repr: Repr::Bcf(bcf),
        };
        take(with, plan, raw)
    }
}

/// Records read front to back into `F`.
///
/// Read with [`next_record`](Self::next_record) (lends one record at a time,
/// allocation-free), or own them: [`read`](Self::read) hands out a record
/// and [`read_into`](Self::read_into) overwrites one, reusing its buffers;
/// [`records`](Self::records) and `into_iter()` hand out owned records,
/// allocating each.
pub struct Reader<F: Fields> {
    core: Core<F::Plan>,
    /// The record [`next_record`](Self::next_record) lends, once there is one.
    current: Option<F>,
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

    // r[impl vcf_values.reader_keys]
    /// Keys for reading this reader's records, resolved as its file
    /// resolves them: for a [`VcfStream`](super::VcfStream), whose file
    /// exists only inside its reader.
    pub fn keys(&self) -> Keys<'_> {
        Keys(&self.core.file)
    }

    // r[impl vcf_reader.query+2]
    /// Read the records of `region` next, as [`ReaderBuilder::query`] would,
    /// reusing this reader's file handle and block cache. The match mode
    /// is the one the reader was opened with ([`Match::StartsIn`] for a
    /// partition reader).
    ///
    /// # Errors
    /// As [`ReaderBuilder::query`]; the reader is unchanged on error.
    pub fn requery(&mut self, region: impl Into<Region>) -> Result<(), VcfReadError> {
        self.requery_many([region])
    }

    // r[impl vcf_reader.query.many+2]
    /// [`requery`](Self::requery) over the union of `regions`, as
    /// [`ReaderBuilder::query_many`].
    ///
    /// # Errors
    /// As [`ReaderBuilder::query`].
    pub fn requery_many<I>(&mut self, regions: I) -> Result<(), VcfReadError>
    where
        I: IntoIterator,
        I::Item: Into<Region>,
    {
        let core = &mut self.core;
        let matching = core.query.as_ref().map_or_else(Match::default, Query::matching);
        let (query, chunks) =
            query::plan(&core.file, regions.into_iter().map(Into::into), matching)?;
        match &mut core.input {
            Framed::Region(input) => input.replan(chunks),
            _ => {
                let input = RegionInput::new(core.file.open_input()?, chunks)?;
                core.input = Framed::Region(Box::new(input));
            }
        }
        core.query = Some(query);
        core.line = None;
        Ok(())
    }

    // r[impl vcf_reader.read_into+2]
    /// The next record, owned; `None` at the end of input. Keep reading into
    /// it with [`read_into`](Self::read_into).
    ///
    /// # Errors
    /// [`VcfReadError::Record`] for a malformed record (with its ordinal and
    /// virtual offset), or a decompression error.
    pub fn read(&mut self) -> Result<Option<F>, VcfReadError> {
        self.core.next(|plan, raw| F::new(plan, raw))
    }

    // r[impl vcf_reader.read_into+2]
    /// Overwrite `rec` with the next record, reusing its buffers; `false`
    /// at the end of input, leaving `rec` as it was.
    ///
    /// # Errors
    /// As [`Self::read`]; `rec` may then hold parts of the failed record.
    pub fn read_into(&mut self, rec: &mut F) -> Result<bool, VcfReadError> {
        Ok(self.core.next(|plan, raw| rec.fill(plan, raw))?.is_some())
    }

    // r[impl vcf_reader.next_record]
    /// The next record, lent until the next call; `None` at the end.
    ///
    /// # Errors
    /// As [`Self::read`].
    pub fn next_record(&mut self) -> Result<Option<&F>, VcfReadError> {
        let Self { core, current } = self;
        let read = match current {
            Some(rec) => core.next(|plan, raw| rec.fill(plan, raw))?.is_some(),
            None => {
                *current = core.next(|plan, raw| F::new(plan, raw))?;
                current.is_some()
            }
        };
        Ok(if read { current.as_ref() } else { None })
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
    match core.next(|plan, raw| F::new(plan, raw)) {
        Ok(Some(rec)) => Some(Ok(rec)),
        Ok(None) => {
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
