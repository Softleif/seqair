//! Read a whole SAM front to back, without an index. See `r[sam.stream]`.

use super::encode::{SamRecordEncoder, SamRecordError};
use crate::bam::{BamHeader, BamHeaderError, RawBamRecord};
use crate::io::text_scan;
use std::io::{self, BufRead};

/// The longest line [`SamStreamReader`] buffers (`r[sam.stream.lines]`).
/// Text of a 2 MiB record (`r[bam.record.max_size]`) fits: its densest
/// encoding, a `B:c` array, is five text bytes per stored byte.
pub const MAX_LINE_LEN: usize = 16 * 1024 * 1024;

/// Why reading a SAM stream failed.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum SamStreamError {
    /// The underlying reader failed — including a truncated gzip or BGZF
    /// stream, whose error is the `io::Error`'s source.
    #[error("I/O error reading SAM")]
    Io {
        #[from]
        source: io::Error,
    },

    #[error("SAM header line {line} is not valid UTF-8")]
    HeaderNotUtf8 { line: u64, source: std::str::Utf8Error },

    #[error("SAM header error")]
    Header {
        #[from]
        source: BamHeaderError,
    },

    /// A line that is not a valid alignment. The stream stays usable.
    #[error("SAM line {line} is invalid")]
    InvalidRecord { line: u64, source: SamRecordError },

    #[error("SAM line {line} is longer than {MAX_LINE_LEN} bytes")]
    LineTooLong { line: u64 },

    /// An earlier I/O error or over-long line ended the stream.
    #[error("SAM stream already failed; an earlier read returned the error")]
    Poisoned,
}

// r[impl sam.stream]
/// A SAM read front to back from any [`BufRead`] — plain text, a
/// [`BgzfStream`](crate::bam::BgzfStream) for bgzipped SAM, a
/// [`BufReader`](std::io::BufReader) over `flate2`'s `MultiGzDecoder` for
/// plain gzip — yielding every alignment line as a [`RawBamRecord`], encoded
/// as htslib encodes it (`r[sam.encode]`).
///
/// ```no_run
/// use seqair::bam::RawBamRecord;
/// use seqair::sam::SamStreamReader;
///
/// let mut reader = SamStreamReader::new(std::io::stdin().lock())?;
/// let mut record = RawBamRecord::default();
/// while reader.read_record(&mut record)? {
///     println!("{}", String::from_utf8_lossy(record.qname()));
/// }
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug)]
pub struct SamStreamReader<R> {
    inner: R,
    header: BamHeader,
    encoder: SamRecordEncoder,
    line: Vec<u8>,
    /// Lines read so far, header lines included.
    line_no: u64,
    failed: bool,
}

impl<R: BufRead> SamStreamReader<R> {
    // r[impl sam.stream.header]
    /// Read the header — the `@` lines at the start of `inner` — leaving
    /// `inner` at the first alignment line.
    ///
    /// # Errors
    /// An I/O error, a header line that is not UTF-8, or a header
    /// [`BamHeader::from_sam_text_allow_no_targets`] rejects.
    pub fn new(mut inner: R) -> Result<Self, SamStreamError> {
        let mut text = String::new();
        let mut line = Vec::new();
        let mut line_no = 0u64;
        while peek(&mut inner)? == Some(b'@') {
            line_no = line_no.saturating_add(1);
            read_line(&mut inner, &mut line, line_no)?;
            let s = std::str::from_utf8(&line)
                .map_err(|source| SamStreamError::HeaderNotUtf8 { line: line_no, source })?;
            text.push_str(s);
            text.push('\n');
        }
        let header = BamHeader::from_sam_text_allow_no_targets(&text)?;
        Ok(Self { inner, header, encoder: SamRecordEncoder::new(), line, line_no, failed: false })
    }

    /// The header read when the stream was opened.
    pub fn header(&self) -> &BamHeader {
        &self.header
    }

    // r[impl sam.stream.lines]
    // r[impl sam.stream.errors]
    /// Read the next line into `record`, reusing its buffer. `Ok(false)` at
    /// the end of the input, which leaves `record` unchanged.
    ///
    /// # Errors
    /// [`SamStreamError::InvalidRecord`] for a line that does not encode —
    /// the stream stays usable. An I/O error or [`SamStreamError::LineTooLong`]
    /// ends the stream: later calls return [`SamStreamError::Poisoned`].
    pub fn read_record(&mut self, record: &mut RawBamRecord) -> Result<bool, SamStreamError> {
        if self.failed {
            return Err(SamStreamError::Poisoned);
        }
        let line_no = self.line_no.saturating_add(1);
        match read_line(&mut self.inner, &mut self.line, line_no) {
            Ok(true) => {}
            Ok(false) => return Ok(false),
            Err(e) => {
                self.failed = true;
                return Err(e);
            }
        }
        self.line_no = line_no;
        self.encoder
            .encode(&self.line, &self.header, record)
            .map_err(|source| SamStreamError::InvalidRecord { line: line_no, source })?;
        Ok(true)
    }

    /// Iterate over the remaining records, one allocation each. For a reused
    /// buffer, loop over [`Self::read_record`] instead.
    pub fn records(&mut self) -> Records<'_, R> {
        Records { reader: self }
    }

    /// Give up the reader, keeping the header.
    pub fn into_header(self) -> BamHeader {
        self.header
    }
}

/// The next byte of `inner` without consuming it; `None` at the end.
fn peek(inner: &mut impl BufRead) -> io::Result<Option<u8>> {
    loop {
        match inner.fill_buf() {
            Ok(buf) => return Ok(buf.first().copied()),
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
}

/// Read one line into `buf` (cleared first) without its `\n` or `\r\n`.
/// `Ok(false)` when `inner` is already at its end.
fn read_line(
    inner: &mut impl BufRead,
    buf: &mut Vec<u8>,
    line: u64,
) -> Result<bool, SamStreamError> {
    buf.clear();
    let mut any = false;
    loop {
        let chunk = match inner.fill_buf() {
            Ok(chunk) => chunk,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
            Err(e) => return Err(e.into()),
        };
        if chunk.is_empty() {
            break;
        }
        any = true;
        let newline = text_scan::find_byte(chunk, b'\n');
        let take = newline.unwrap_or(chunk.len());
        if buf.len().saturating_add(take) > MAX_LINE_LEN {
            return Err(SamStreamError::LineTooLong { line });
        }
        buf.extend_from_slice(chunk.get(..take).unwrap_or_default());
        inner.consume(if newline.is_some() { take.saturating_add(1) } else { take });
        if newline.is_some() {
            break;
        }
    }
    if buf.last() == Some(&b'\r') {
        buf.pop();
    }
    Ok(any)
}

/// Iterator over a [`SamStreamReader`]'s records. See
/// [`SamStreamReader::records`].
#[derive(Debug)]
pub struct Records<'a, R> {
    reader: &'a mut SamStreamReader<R>,
}

impl<R: BufRead> Iterator for Records<'_, R> {
    type Item = Result<RawBamRecord, SamStreamError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.reader.failed {
            return None;
        }
        let mut record = RawBamRecord::default();
        match self.reader.read_record(&mut record) {
            Ok(true) => Some(Ok(record)),
            Ok(false) => None,
            Err(e) => Some(Err(e)),
        }
    }
}
