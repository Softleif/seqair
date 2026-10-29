//! [`AlignmentStream`]: read BAM, SAM (plain, bgzipped or gzipped) or CRAM
//! front to back from any reader, deciding the format from its first bytes.
//! See `r[unified.stream.detect+2]`.

use crate::bam::{
    BamHeader, RawBamRecord,
    bgzf::{BgzfError, BgzfStream, EofMarker, Workers},
    header::BAM_MAGIC,
    stream::{BamStreamError, BamStreamReader},
};
use crate::cram::stream::{CramContainerDecoder, CramStreamError, CramStreamReader};
use crate::sam::{SamStreamError, SamStreamReader};
use std::collections::VecDeque;
use std::io::{self, BufRead, BufReader, Cursor, Read};
use std::path::PathBuf;

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

    /// CRAM is decoded against a reference; give one with
    /// [`AlignmentStreamBuilder::reference`].
    #[error("the input is CRAM, which needs a FASTA reference (AlignmentStreamBuilder::reference)")]
    CramNeedsReference,

    /// BAM in ordinary gzip rather than BGZF. htslib reads it; recompress with
    /// `samtools view -b` to stream it here.
    #[error("the input is BAM compressed with gzip rather than BGZF; recompress it with samtools")]
    GzipBamNotSupported,

    #[error("the input is uncompressed BAM, which is not supported")]
    UncompressedBamNotSupported,
}

type Source = std::io::Chain<Cursor<Vec<u8>>, Box<dyn Read + Send>>;

enum Inner {
    Bam(Box<BamStreamReader<Source>>),
    Sam(Box<SamStreamReader<Box<dyn BufRead + Send>>>),
    Cram(Box<CramRecords>),
}

/// A CRAM stream decoded one container at a time on the reading thread.
struct CramRecords {
    reader: CramStreamReader<Source>,
    decoder: CramContainerDecoder,
    pending: VecDeque<RawBamRecord>,
    scratch: Vec<RawBamRecord>,
}

impl CramRecords {
    fn read_record(&mut self, record: &mut RawBamRecord) -> Result<bool, CramStreamError> {
        while self.pending.is_empty() {
            let Some(container) = self.reader.next_container()? else { return Ok(false) };
            self.scratch.clear();
            self.decoder.decode(&container, &mut self.scratch)?;
            self.pending.extend(self.scratch.drain(..));
        }
        match self.pending.pop_front() {
            Some(next) => {
                *record = next;
                Ok(true)
            }
            None => Ok(false),
        }
    }
}

/// Configures an [`AlignmentStream`]. The options apply to the BGZF formats.
#[derive(Debug)]
pub struct AlignmentStreamBuilder<R> {
    inner: R,
    workers: Workers,
    eof_marker: EofMarker,
    reference: Option<PathBuf>,
}

impl<R: Read + Send + 'static> AlignmentStreamBuilder<R> {
    /// Inflate up to `n` BGZF batches at once on the reading thread's rayon
    /// pool, else the global one (`r[bgzf.stream.current_pool]`).
    #[must_use]
    pub fn threads(mut self, n: std::num::NonZeroU8) -> Self {
        self.workers = Workers::Threads { n, pool: None };
        self
    }

    /// Inflate up to `n` BGZF batches at once on `pool`
    /// (`r[bgzf.stream.shared_pool]`).
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
    /// decoded against. Other formats ignore it.
    #[must_use]
    pub fn reference(mut self, fasta: impl Into<PathBuf>) -> Self {
        self.reference = Some(fasta.into());
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
    pub fn build(self) -> Result<AlignmentStream, AlignmentStreamError> {
        AlignmentStream::open(self)
    }
}

// r[impl unified.stream.records]
/// Any BAM or SAM input, read front to back as [`RawBamRecord`]s.
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
/// # Ok::<(), seqair::reader::AlignmentStreamError>(())
/// ```
pub struct AlignmentStream {
    inner: Inner,
    format: StreamFormat,
}

impl std::fmt::Debug for AlignmentStream {
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
        }
    }

    // r[impl unified.stream.detect+2]
    fn open<R: Read + Send + 'static>(
        options: AlignmentStreamBuilder<R>,
    ) -> Result<Self, AlignmentStreamError> {
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
                    return Ok(Self {
                        inner: Inner::Bam(Box::new(reader)),
                        format: StreamFormat::Bam,
                    });
                }
                let text: Box<dyn BufRead + Send> = Box::new(Cursor::new(magic).chain(bgzf));
                return Self::sam(text, StreamFormat::BgzfSam);
            }
            let source = Cursor::new(prefix).chain(inner);
            let mut gzip = BufReader::new(flate2::read::MultiGzDecoder::new(source));
            let mut magic = Vec::with_capacity(4);
            read_up_to(&mut gzip, &mut magic, 4)?;
            if magic == BAM_MAGIC {
                return Err(AlignmentStreamError::GzipBamNotSupported);
            }
            let text: Box<dyn BufRead + Send> = Box::new(Cursor::new(magic).chain(gzip));
            return Self::sam(text, StreamFormat::GzipSam);
        }
        if prefix.starts_with(b"CRAM") {
            let fasta = options.reference.ok_or(AlignmentStreamError::CramNeedsReference)?;
            let reader = CramStreamReader::new(Cursor::new(prefix).chain(inner))?;
            let decoder = reader.decoder(&fasta)?;
            return Ok(Self {
                inner: Inner::Cram(Box::new(CramRecords {
                    reader,
                    decoder,
                    pending: VecDeque::new(),
                    scratch: Vec::new(),
                })),
                format: StreamFormat::Cram,
            });
        }
        if prefix.starts_with(&BAM_MAGIC) {
            return Err(AlignmentStreamError::UncompressedBamNotSupported);
        }
        let text: Box<dyn BufRead + Send> =
            Box::new(BufReader::new(Cursor::new(prefix).chain(inner)));
        Self::sam(text, StreamFormat::Sam)
    }

    fn sam(
        text: Box<dyn BufRead + Send>,
        format: StreamFormat,
    ) -> Result<Self, AlignmentStreamError> {
        let reader = SamStreamReader::new(text)?;
        Ok(Self { inner: Inner::Sam(Box::new(reader)), format })
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
    pub fn read_record(&mut self, record: &mut RawBamRecord) -> Result<bool, AlignmentStreamError> {
        Ok(match &mut self.inner {
            Inner::Bam(r) => r.read_record(record)?,
            Inner::Sam(r) => r.read_record(record)?,
            Inner::Cram(r) => r.read_record(record)?,
        })
    }
}

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
