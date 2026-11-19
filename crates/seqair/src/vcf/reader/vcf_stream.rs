//! [`VcfStream`]: VCF or BCF from input that cannot seek — a pipe, standard
//! input (`r[vcf_reader.stream]`).

use super::error::{VcfFormat, VcfReadError};
use super::fields::{Fields, SampleSel};
use super::file::{
    Container, DataStart, FILE_BUFFER, HeaderStart, Input, Source, VcfFile, container, format_of,
    read_bcf_header, read_failed, read_text_header,
};
use super::stream::{Framed, Reader, ReaderBuilder, Start, batches};
use crate::bam::bgzf::{BgzfStream, Workers};
use rayon_core::ThreadPool;
use std::io::{BufReader, Cursor, Read};
use std::num::NonZeroUsize;
use std::sync::Arc;

// r[impl vcf_reader.stream]
/// VCF, BGZF VCF or BCF read once, front to back, from input that cannot
/// seek. Its one [`reader`](Self::reader) yields the same [`Reader`] a
/// [`VcfFile`] does; a stream has no index, so it has no queries.
///
/// ```
/// use seqair::vcf::{Locus, VcfStream};
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// # let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/vcf/rastair.bcf");
/// let input = std::fs::File::open(path)?; // or a pipe, or `VcfStream::stdin()`
/// let mut reader = VcfStream::new(input).reader::<Locus>().open()?;
/// let mut n = 0;
/// while reader.next_record()?.is_some() {
///     n += 1;
/// }
/// assert_eq!(n, 400);
/// # Ok(()) }
/// ```
pub struct VcfStream {
    input: Box<dyn Read + Send>,
}

impl std::fmt::Debug for VcfStream {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("VcfStream").finish_non_exhaustive()
    }
}

impl VcfStream {
    /// Records from `input`; nothing is read until the reader opens.
    pub fn new(input: impl Read + Send + 'static) -> Self {
        Self { input: Box::new(input) }
    }

    /// Records from standard input.
    pub fn stdin() -> Self {
        Self::new(std::io::stdin())
    }

    /// The stream's reader of records into `F`, to configure and
    /// [`open`](StreamReaderBuilder::open).
    pub fn reader<F: Fields>(self) -> StreamReaderBuilder<F> {
        StreamReaderBuilder { input: self.input, workers: None, sample: None, plan: None }
    }
}

/// Configures a [`VcfStream`]'s reader: as a [`ReaderBuilder`], without
/// queries.
#[must_use]
pub struct StreamReaderBuilder<F: Fields> {
    input: Box<dyn Read + Send>,
    workers: Option<(NonZeroUsize, Option<Arc<ThreadPool>>)>,
    sample: Option<SampleSel>,
    plan: Option<F::Plan>,
}

impl<F: Fields> std::fmt::Debug for StreamReaderBuilder<F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("StreamReaderBuilder").field("sample", &self.sample).finish_non_exhaustive()
    }
}

impl<F: Fields> StreamReaderBuilder<F> {
    /// As [`ReaderBuilder::threads`].
    pub fn threads(mut self, n: NonZeroUsize) -> Self {
        self.workers = Some((n, None));
        self
    }

    /// As [`ReaderBuilder::threads_in`].
    pub fn threads_in(mut self, n: NonZeroUsize, pool: impl Into<Arc<ThreadPool>>) -> Self {
        self.workers = Some((n, Some(pool.into())));
        self
    }

    /// As [`ReaderBuilder::sample`].
    pub fn sample(mut self, sample: impl Into<SampleSel>) -> Self {
        self.sample = Some(sample.into());
        self
    }

    /// As [`ReaderBuilder::with_plan`].
    pub fn with_plan(mut self, plan: F::Plan) -> Self {
        self.plan = Some(plan);
        self
    }

    /// Detect the format, read the header and plan `F` against it.
    ///
    /// # Errors
    /// As [`VcfFile::open`] and [`ReaderBuilder::open`].
    pub fn open(self) -> Result<Reader<F>, VcfReadError> {
        let workers = match &self.workers {
            Some((n, pool)) => Workers::Threads { n: batches(*n), pool: pool.clone() },
            None => Workers::Inline,
        };
        let (start, input) = read_stream_header(self.input, workers)?;
        let file = VcfFile::new(Source::Stream, start, None)?;
        let mut builder = ReaderBuilder::<F>::new(&file);
        if let Some(sample) = self.sample {
            builder = builder.sample(sample);
        }
        if let Some(plan) = self.plan {
            builder = builder.with_plan(plan);
        }
        builder.start(Start::Stream(input))
    }
}

/// Detect `input`'s format and read its header without seeking: the format,
/// header and data start, and the input positioned at the first record.
fn read_stream_header(
    mut input: Box<dyn Read + Send>,
    workers: Workers,
) -> Result<(HeaderStart, Framed), VcfReadError> {
    let mut prefix = Vec::with_capacity(64);
    (&mut input).take(64).read_to_end(&mut prefix).map_err(read_failed)?;
    let container = container(&prefix);
    let plain_format = format_of(&prefix, false);
    let chained: Box<dyn Read + Send> = Box::new(Cursor::new(prefix).chain(input));
    let mut input = Input::Stream(BufReader::with_capacity(FILE_BUFFER, chained));
    match container {
        Container::Gzip => Err(VcfReadError::NotBgzf),
        Container::Plain => {
            let format = plain_format?;
            let (header, len, lines) = read_text_header(&mut input)?;
            Ok(((format, header, DataStart::After(len), lines), Framed::Plain(input)))
        }
        Container::Bgzf => {
            let mut stream = BgzfStream::builder(input).workers(workers).build();
            let format = format_of(stream.fill()?, true)?;
            let (header, len, lines) = if format == VcfFormat::Bcf {
                let (header, len) = read_bcf_header(&mut stream)?;
                (header, len, 0)
            } else {
                read_text_header(&mut stream)?
            };
            Ok(((format, header, DataStart::After(len), lines), Framed::Bgzf(Box::new(stream))))
        }
    }
}
