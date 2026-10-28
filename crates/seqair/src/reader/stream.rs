//! [`AlignmentStream`]: read BAM or SAM (plain, bgzipped or gzipped) front to
//! back from any reader, deciding the format from its first bytes.
//! See `r[unified.stream.detect]`.

use crate::bam::{
    BamHeader, RawBamRecord,
    bgzf::{BgzfError, BgzfStream, EofMarker},
    header::BAM_MAGIC,
    stream::{BamStreamError, BamStreamReader},
};
use crate::sam::{SamStreamError, SamStreamReader};
use std::io::{self, BufRead, BufReader, Cursor, Read};

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

    #[error("the input is CRAM, which can only be read with an index (IndexedReader)")]
    CramNotSupported,

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
}

/// Configures an [`AlignmentStream`]. The options apply to the BGZF formats.
#[derive(Debug)]
pub struct AlignmentStreamBuilder<R> {
    inner: R,
    threads: usize,
    eof_marker: EofMarker,
}

impl<R: Read + Send + 'static> AlignmentStreamBuilder<R> {
    /// Inflate BGZF blocks on `threads` worker threads (`r[bgzf.stream.threads]`).
    #[must_use]
    pub fn threads(mut self, threads: usize) -> Self {
        self.threads = threads;
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
/// let mut stream = AlignmentStream::builder(std::io::stdin()).threads(2).build()?;
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
        AlignmentStreamBuilder { inner, threads: 0, eof_marker: EofMarker::Required }
    }

    // r[impl unified.stream.detect]
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
                    .threads(options.threads)
                    .eof_marker(options.eof_marker)
                    .build()?;
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
            return Err(AlignmentStreamError::CramNotSupported);
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
