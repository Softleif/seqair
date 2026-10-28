//! Read a whole BAM front to back, without an index and without seeking.
//! See `r[bam.stream]`.

use super::{
    BamHeader, BamHeaderError, DecodeError, RawBamRecord,
    bgzf::{BgzfError, BgzfStream},
    raw_record::Refill,
    record::MAX_RECORD_SIZE,
};
use std::io::Read;

/// Why reading a BAM stream failed.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum BamStreamError {
    #[error("BGZF error")]
    Bgzf {
        #[from]
        source: BgzfError,
    },

    #[error("BAM header error")]
    Header {
        #[from]
        source: BamHeaderError,
    },

    /// A record's `block_size` is below the 32 fixed bytes or above the
    /// record size limit (`r[bam.record.max_size]`).
    #[error("BAM record {record} declares block_size {size}, outside 32..={MAX_RECORD_SIZE}")]
    RecordSizeOutOfRange { record: u64, size: i32 },

    /// The input ended inside a record (its `block_size` or its payload).
    #[error("BAM input ends inside record {record}")]
    TruncatedRecord { record: u64 },

    /// A record's bytes are framed correctly but are not a valid record. The
    /// stream stays aligned, so reading on yields the next record.
    #[error("BAM record {record} is invalid")]
    InvalidRecord { record: u64, source: DecodeError },

    /// An earlier error left the stream misaligned; nothing more can be read.
    #[error("BAM stream already failed; an earlier read returned the error")]
    Poisoned,
}

// r[impl bam.stream]
/// A BAM file read front to back from any [`Read`] — `stdin`, a pipe, a
/// socket — with no index and no `Seek` bound.
///
/// Yields every record, unplaced ones included, as a [`RawBamRecord`] in the
/// order the file stores them.
///
/// ```no_run
/// use seqair::bam::{RawBamRecord, stream::BamStreamReader};
///
/// let mut reader = BamStreamReader::new(std::io::stdin().lock())?;
/// let mut record = RawBamRecord::default();
/// while reader.read_record(&mut record)? {
///     println!("{}", String::from_utf8_lossy(record.qname()));
/// }
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug)]
pub struct BamStreamReader<R> {
    stream: BgzfStream<R>,
    header: BamHeader,
    /// Records framed so far, for error messages (1-based numbering).
    records: u64,
    failed: bool,
}

impl<R: Read> BamStreamReader<R> {
    /// Read the header from `inner` (BGZF-compressed BAM), inflating on the
    /// calling thread and requiring the EOF marker.
    ///
    /// # Errors
    /// A BGZF or header error.
    pub fn new(inner: R) -> Result<Self, BamStreamError> {
        Self::from_bgzf(BgzfStream::new(inner))
    }

    /// Read the header from a configured [`BgzfStream`] (threads, EOF-marker
    /// handling) positioned at the start of a BAM file.
    ///
    /// # Errors
    /// A BGZF or header error.
    pub fn from_bgzf(mut stream: BgzfStream<R>) -> Result<Self, BamStreamError> {
        let header = BamHeader::parse_inner(&mut stream)?;
        Ok(Self { stream, header, records: 0, failed: false })
    }

    /// Read the header from a stream whose `BAM\1` magic the caller has
    /// already read and checked (format sniffing).
    pub(crate) fn after_magic(mut stream: BgzfStream<R>) -> Result<Self, BamStreamError> {
        let header = BamHeader::parse_after_magic(&mut stream)?;
        Ok(Self { stream, header, records: 0, failed: false })
    }

    /// The header read when the stream was opened.
    pub fn header(&self) -> &BamHeader {
        &self.header
    }

    // r[impl bam.stream.framing]
    /// Read the next record into `record`, reusing its buffer. `Ok(false)` at
    /// the end of the stream, which leaves `record` unchanged.
    ///
    /// # Errors
    /// [`BamStreamError::InvalidRecord`] for a record that fails validation —
    /// the stream stays usable. Anything else ends the stream: later calls
    /// return [`BamStreamError::Poisoned`] (or the BGZF layer's own poison).
    pub fn read_record(&mut self, record: &mut RawBamRecord) -> Result<bool, BamStreamError> {
        if self.failed {
            return Err(BamStreamError::Poisoned);
        }
        let result = self.read_record_inner(record);
        if let Err(e) = &result
            && !matches!(e, BamStreamError::InvalidRecord { .. })
        {
            self.failed = true;
        }
        result
    }

    fn read_record_inner(&mut self, record: &mut RawBamRecord) -> Result<bool, BamStreamError> {
        let Some(len) = self.next_frame()? else { return Ok(false) };
        let (number, stream) = (self.records, &mut self.stream);
        record
            .refill(&self.header, |payload| {
                payload.resize(len, 0);
                stream.read_exact_into(payload).map_err(|e| truncated(number, e))
            })
            .map_err(|e| match e {
                Refill::Fill(e) => e,
                Refill::Invalid(source) => BamStreamError::InvalidRecord { record: number, source },
            })?;
        Ok(true)
    }

    // r[impl bam.stream.payloads]
    /// Append the next record's payload (its bytes after `block_size`) to
    /// `out` *without* validating it, for callers that validate elsewhere —
    /// typically [`RawBamRecord::from_payload`] on worker threads. `Ok(false)`
    /// at the end of the stream. Framing is checked as by [`Self::read_record`].
    ///
    /// # Errors
    /// As [`Self::read_record`], except that nothing here is an
    /// [`BamStreamError::InvalidRecord`]: every error ends the stream.
    pub fn read_payload_into(&mut self, out: &mut Vec<u8>) -> Result<bool, BamStreamError> {
        if self.failed {
            return Err(BamStreamError::Poisoned);
        }
        let start = out.len();
        let result = self.next_frame().and_then(|frame| {
            let Some(len) = frame else { return Ok(false) };
            out.resize(start.saturating_add(len), 0);
            let payload = out.get_mut(start..).unwrap_or_default();
            self.stream.read_exact_into(payload).map_err(|e| truncated(self.records, e))?;
            Ok(true)
        });
        if result.is_err() {
            out.truncate(start);
            self.failed = true;
        }
        result
    }

    /// Read the next `block_size`: `None` at a clean end, else the payload
    /// length, checked against the limits.
    fn next_frame(&mut self) -> Result<Option<usize>, BamStreamError> {
        // Zero bytes before a record is the end; anything else must be whole.
        if self.stream.fill()?.is_empty() {
            return Ok(None);
        }
        self.records = self.records.saturating_add(1);
        let number = self.records;
        let mut size = [0u8; 4];
        self.stream.read_exact_into(&mut size).map_err(|e| truncated(number, e))?;
        let size = i32::from_le_bytes(size);
        usize::try_from(size)
            .ok()
            .filter(|n| (32..=MAX_RECORD_SIZE).contains(n))
            .map(Some)
            .ok_or(BamStreamError::RecordSizeOutOfRange { record: number, size })
    }

    /// Iterate over the remaining records, one allocation each. For a
    /// reused buffer, loop over [`Self::read_record`] instead.
    pub fn records(&mut self) -> Records<'_, R> {
        Records { reader: self }
    }

    /// Give up the reader, keeping the header.
    pub fn into_header(self) -> BamHeader {
        self.header
    }
}

/// An end of input inside record `record` is truncation.
fn truncated(record: u64, e: BgzfError) -> BamStreamError {
    match e {
        BgzfError::UnexpectedEof => BamStreamError::TruncatedRecord { record },
        e => BamStreamError::Bgzf { source: e },
    }
}

/// Iterator over a [`BamStreamReader`]'s records. See
/// [`BamStreamReader::records`].
#[derive(Debug)]
pub struct Records<'a, R> {
    reader: &'a mut BamStreamReader<R>,
}

impl<R: Read> Iterator for Records<'_, R> {
    type Item = Result<RawBamRecord, BamStreamError>;

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
