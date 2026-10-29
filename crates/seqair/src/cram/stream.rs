//! Read a whole CRAM front to back, without an index and without seeking.
//! [`CramStreamReader`] reads containers off the input; a
//! [`CramContainerDecoder`] per thread turns them into [`RawBamRecord`]s.
//! See `r[cram.stream]`.

use super::{
    block,
    compression_header::CompressionHeader,
    container::{ContainerHeader, RefSeqId},
    reader::{self, CramError},
    slice::{self, Keep, SliceHeader, SliceReference, SliceScratch, TidWindow},
    varint,
};
use crate::bam::{BamHeader, DecodeError, RawBamRecord};
use crate::fasta::{FastaError, IndexedFastaReader};
use crate::reader::Tid;
use crate::sam::encode::{RecordEncodeError, encode_alignment};
use core::range::RangeInclusive;
use seqair_types::{IupacBase, Pos0, SmolStr};
use std::io::{self, Read};
use std::path::Path;
use std::sync::Arc;
use tracing::warn;

/// Why reading or decoding a CRAM stream failed.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum CramStreamError {
    #[error("reading the CRAM input failed")]
    Io { source: io::Error },

    #[error(transparent)]
    Cram {
        #[from]
        source: CramError,
    },

    /// The input ended inside the file definition.
    #[error("CRAM input ends inside the 26-byte file definition")]
    TruncatedFileDefinition,

    /// The input ended inside a container's header or body. Container 0 is
    /// the file header container.
    #[error("CRAM input ends inside container {container}")]
    TruncatedContainer { container: u64 },

    /// A decoded record is not a valid BAM record for the header.
    #[error("CRAM record {record} is invalid")]
    InvalidRecord { record: u64, source: DecodeError },

    /// A decoded record has a field BAM cannot hold.
    #[error("CRAM record {record} cannot be encoded as BAM")]
    RecordEncode { record: u64, source: RecordEncodeError },

    /// An earlier error left the stream misaligned; nothing more can be read.
    #[error("CRAM stream already failed; an earlier read returned the error")]
    Poisoned,
}

/// The parsed file header, shared by a reader and its decoders.
#[derive(Debug, Clone)]
struct StreamShared {
    header: BamHeader,
    /// `@RG` ids in header order, as CRAM's `RG` series numbers them.
    read_group_ids: Vec<SmolStr>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum State {
    Open,
    Ended,
    Failed,
}

// r[impl cram.stream]
/// A CRAM v3.0/v3.1 file read front to back from any [`Read`] — `stdin`, a
/// pipe — with no index and no `Seek` bound.
///
/// [`Self::next_container`] is cheap; decoding is the CPU cost, and can run on
/// other threads: hand each [`CramContainer`] to a [`CramContainerDecoder`]
/// (one per thread, from [`Self::decoder`] and [`CramContainerDecoder::fork`]).
///
/// ```no_run
/// use seqair::cram::stream::CramStreamReader;
///
/// let mut reader = CramStreamReader::new(std::io::stdin().lock())?;
/// let mut decoder = reader.decoder("ref.fa".as_ref())?;
/// let mut records = Vec::new();
/// while let Some(container) = reader.next_container()? {
///     records.clear();
///     decoder.decode(&container, &mut records)?;
///     for record in &records {
///         println!("{}", record.qname().escape_ascii());
///     }
/// }
/// # Ok::<(), seqair::cram::stream::CramStreamError>(())
/// ```
#[derive(Debug)]
pub struct CramStreamReader<R> {
    inner: R,
    shared: Arc<StreamShared>,
    /// Containers read so far, the header container included.
    containers: u64,
    state: State,
}

// r[impl cram.stream.container]
/// One data container, read whole: its header and the bytes of its blocks.
/// Owned and `Send`, for decoding on any thread.
#[derive(Debug)]
pub struct CramContainer {
    number: u64,
    header: ContainerHeader,
    data: Vec<u8>,
}

impl CramContainer {
    /// The container's 0-based number in the file; the file header container
    /// is 0.
    pub fn number(&self) -> u64 {
        self.number
    }

    /// How many records the container header says it holds.
    pub fn num_records(&self) -> u64 {
        u64::try_from(self.header.num_records).unwrap_or(0)
    }
}

impl<R: Read> CramStreamReader<R> {
    /// Read the file definition and the file header container.
    ///
    /// # Errors
    /// Not CRAM v3, a truncated or corrupt header container, or a header the
    /// indexed reader would also reject.
    pub fn new(mut inner: R) -> Result<Self, CramStreamError> {
        let mut file_def = [0u8; 26];
        inner.read_exact(&mut file_def).map_err(|e| match e.kind() {
            io::ErrorKind::UnexpectedEof => CramStreamError::TruncatedFileDefinition,
            _ => CramStreamError::Io { source: e },
        })?;
        let [m0, m1, m2, m3, major, minor, ..] = file_def;
        // r[impl cram.file.magic]
        if [m0, m1, m2, m3] != *b"CRAM" {
            return Err(CramError::InvalidMagic { found: [m0, m1, m2, m3] }.into());
        }
        // r[impl cram.scope.versions]
        if major != 3 {
            return Err(CramError::UnsupportedVersion { major, minor }.into());
        }

        // r[impl cram.file.header_container]
        let (_, data) = read_container(&mut inner, 0)?
            .ok_or(CramStreamError::TruncatedContainer { container: 0 })?;
        let header = reader::header_from_container_data(&data)?;
        // r[impl unified.sort_order]
        header.validate_sort_order().map_err(CramError::from)?;
        let read_group_ids = reader::parse_read_group_ids(header.header_text());
        Ok(Self {
            inner,
            shared: Arc::new(StreamShared { header, read_group_ids }),
            containers: 1,
            state: State::Open,
        })
    }

    /// The header read when the stream was opened.
    pub fn header(&self) -> &BamHeader {
        &self.shared.header
    }

    /// A decoder for this file's containers, reconstructing sequences against
    /// the FASTA at `fasta` (it needs a `.fai`, and a `.gzi` if bgzipped).
    ///
    /// # Errors
    /// The FASTA or its index cannot be opened.
    pub fn decoder(&self, fasta: &Path) -> Result<CramContainerDecoder, CramStreamError> {
        let fasta = IndexedFastaReader::open(fasta).map_err(CramError::from)?;
        Ok(CramContainerDecoder::new(Arc::clone(&self.shared), fasta))
    }

    // r[impl cram.stream.framing]
    /// Read the next container that holds records. `Ok(None)` at the end of
    /// the stream.
    ///
    /// # Errors
    /// Input that ends inside a container, a corrupt container header, or an
    /// I/O error. Each ends the stream: later calls return
    /// [`CramStreamError::Poisoned`].
    pub fn next_container(&mut self) -> Result<Option<CramContainer>, CramStreamError> {
        match self.state {
            State::Ended => return Ok(None),
            State::Failed => return Err(CramStreamError::Poisoned),
            State::Open => {}
        }
        loop {
            let number = self.containers;
            let read = read_container(&mut self.inner, number);
            let Some((header, data)) = read.inspect_err(|_| self.state = State::Failed)? else {
                // r[impl cram.file.eof]
                warn!(containers = number, "CRAM stream ends without an EOF container");
                self.state = State::Ended;
                return Ok(None);
            };
            self.containers = number.saturating_add(1);
            if header.is_eof() {
                self.state = State::Ended;
                return Ok(None);
            }
            if header.num_records > 0 {
                return Ok(Some(CramContainer { number, header, data }));
            }
        }
    }

    /// Give up the reader, keeping the header.
    pub fn into_header(self) -> BamHeader {
        Arc::unwrap_or_clone(self.shared).header
    }
}

/// Copies every byte read through it into `bytes`, so a container header
/// read field by field can be CRC-checked whole.
struct Recording<'a, R> {
    inner: &'a mut R,
    bytes: &'a mut Vec<u8>,
}

impl<R: Read> Read for Recording<'_, R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n = self.inner.read(buf)?;
        self.bytes.extend_from_slice(buf.get(..n).unwrap_or_default());
        Ok(n)
    }
}

// r[impl cram.stream.container_header]
/// Read one container (header and body) from `inner`. `None` when the input
/// ends before its first byte.
fn read_container<R: Read>(
    inner: &mut R,
    number: u64,
) -> Result<Option<(ContainerHeader, Vec<u8>)>, CramStreamError> {
    let truncated = |e: io::Error| match e.kind() {
        io::ErrorKind::UnexpectedEof => CramStreamError::TruncatedContainer { container: number },
        _ => CramStreamError::Io { source: e },
    };

    let mut bytes = Vec::with_capacity(64);
    // The first byte tells a clean end from a cut one.
    let mut first = [0u8; 1];
    loop {
        match inner.read(&mut first) {
            Ok(0) => return Ok(None),
            Ok(_) => break,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(source) => return Err(CramStreamError::Io { source }),
        }
    }
    bytes.extend_from_slice(&first);
    {
        let mut r = Recording { inner: &mut *inner, bytes: &mut bytes };
        let mut rest_of_length = [0u8; 3];
        r.read_exact(&mut rest_of_length).map_err(truncated)?;
        // reference id, alignment start, span, record count
        for _ in 0..4 {
            varint::read_itf8(&mut r).map_err(truncated)?;
        }
        // record counter, bases
        for _ in 0..2 {
            varint::read_ltf8(&mut r).map_err(truncated)?;
        }
        varint::read_itf8(&mut r).map_err(truncated)?; // block count
        let landmarks = varint::read_itf8(&mut r).map_err(truncated)?;
        reader::check_alloc_size(
            usize::try_from(landmarks).unwrap_or(usize::MAX).saturating_mul(4),
            "container landmark array",
        )?;
        for _ in 0..landmarks {
            varint::read_itf8(&mut r).map_err(truncated)?;
        }
        let mut crc = [0u8; 4];
        r.read_exact(&mut crc).map_err(truncated)?;
    }
    let header = ContainerHeader::parse(&bytes)?;

    let len = usize::try_from(header.length)
        .map_err(|_| CramError::InvalidLength { value: header.length })?;
    reader::check_alloc_size(len, "container data")?;
    // Grown as the bytes arrive, so a corrupt length cannot allocate ahead.
    let mut data = Vec::new();
    let got = inner
        .take(u64::try_from(len).unwrap_or(u64::MAX))
        .read_to_end(&mut data)
        .map_err(|source| CramStreamError::Io { source })?;
    if got != len {
        return Err(CramStreamError::TruncatedContainer { container: number });
    }
    Ok(Some((header, data)))
}

// r[impl cram.stream.decoder]
/// Decodes [`CramContainer`]s of one file into [`RawBamRecord`]s. Keep one
/// per thread ([`Self::fork`]); its buffers and codec tables are reused from
/// container to container.
pub struct CramContainerDecoder {
    shared: Arc<StreamShared>,
    fasta: IndexedFastaReader,
    scratch: SliceScratch<IupacBase>,
    ref_buf: Vec<u8>,
    windows: Vec<TidWindow<u8>>,
    decode_md: bool,
}

impl std::fmt::Debug for CramContainerDecoder {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("CramContainerDecoder").finish_non_exhaustive()
    }
}

impl CramContainerDecoder {
    fn new(shared: Arc<StreamShared>, fasta: IndexedFastaReader) -> Self {
        Self {
            shared,
            fasta,
            scratch: SliceScratch::default(),
            ref_buf: Vec::new(),
            windows: Vec::new(),
            decode_md: true,
        }
    }

    /// A decoder for another thread, sharing the header and the FASTA index.
    ///
    /// # Errors
    /// The FASTA cannot be reopened.
    pub fn fork(&self) -> Result<Self, CramStreamError> {
        let fasta = self.fasta.fork().map_err(CramError::from)?;
        Ok(Self { decode_md: self.decode_md, ..Self::new(Arc::clone(&self.shared), fasta) })
    }

    // r[impl cram.stream.md_nm.option]
    /// Whether mapped records that store no `MD`/`NM` get them generated, as
    /// htslib does by default (`r[cram.stream.md_nm]`); `false` is
    /// `samtools view --input-fmt-option decode_md=0`. On by default.
    #[must_use]
    pub fn decode_md(mut self, on: bool) -> Self {
        self.decode_md = on;
        self
    }

    /// The file header records are validated against.
    pub fn header(&self) -> &BamHeader {
        &self.shared.header
    }

    // r[impl cram.stream.all_records]
    /// Append every record of `container` to `out`, in file order.
    ///
    /// # Errors
    /// A decoding error, a missing or mismatching reference, or a record that
    /// is not valid BAM. `out` is then as it was.
    pub fn decode(
        &mut self,
        container: &CramContainer,
        out: &mut Vec<RawBamRecord>,
    ) -> Result<(), CramStreamError> {
        let start = out.len();
        let result = self.decode_inner(container, out);
        if result.is_err() {
            out.truncate(start);
        }
        result
    }

    fn decode_inner(
        &mut self,
        container: &CramContainer,
        out: &mut Vec<RawBamRecord>,
    ) -> Result<(), CramStreamError> {
        let data = container.data.as_slice();
        let (comp_block, _) = block::parse_block_with(data, &mut self.scratch.codec)?;
        if comp_block.content_type != block::ContentType::CompressionHeader {
            return Err(
                CramError::ExpectedCompressionHeader { found: comp_block.content_type }.into()
            );
        }
        let ch = CompressionHeader::parse(&comp_block.data)?;
        let header = &self.shared.header;
        let read_group_ids = self.shared.read_group_ids.as_slice();

        for &landmark in &container.header.landmarks {
            let offset = usize::try_from(landmark)
                .map_err(|_| CramError::InvalidLength { value: landmark })?;
            let slice_data =
                data.get(offset..).ok_or(CramError::Truncated { context: "slice offset" })?;
            let (sh_block, _) = block::parse_block_with(slice_data, &mut self.scratch.codec)?;
            if sh_block.content_type != block::ContentType::SliceHeader {
                return Err(CramError::ExpectedSliceHeader { found: sh_block.content_type }.into());
            }
            let sh = SliceHeader::parse(&sh_block.data)?;
            let slice_ref = sh.ref_seq_id.resolve(header).map_err(CramError::from)?;

            // r[impl cram.stream.reference_window]
            let reference = match slice_ref {
                RefSeqId::Ref(tid) => {
                    let begin =
                        u64::try_from(sh.alignment_start.max(1)).unwrap_or(1).saturating_sub(1);
                    let span = u64::try_from(sh.alignment_span).unwrap_or(0);
                    let end = begin.saturating_add(span).min(header.target_len(tid).unwrap_or(0));
                    let found = fetch(&mut self.fasta, header, tid, begin, end, &mut self.ref_buf)?;
                    SliceReference::Window {
                        bases: found.then_some(self.ref_buf.as_slice()),
                        start: begin.cast_signed(),
                    }
                }
                RefSeqId::Unmapped => SliceReference::Window { bases: Some(&[]), start: 0 },
                RefSeqId::Multi => {
                    let spans = slice::multi_ref_spans(
                        &ch,
                        data,
                        offset,
                        header,
                        read_group_ids,
                        &mut self.scratch,
                    )?;
                    self.windows.clear();
                    for (tid, span) in spans {
                        let begin = u64::try_from(span.start).unwrap_or(0);
                        let end = u64::try_from(span.end)
                            .unwrap_or(0)
                            .min(header.target_len(tid).unwrap_or(0));
                        if !fetch(&mut self.fasta, header, tid, begin, end, &mut self.ref_buf)? {
                            return Err(CramError::MissingReference {
                                contig: header.target_name(tid).unwrap_or("?").into(),
                            }
                            .into());
                        }
                        let bases = self.ref_buf.to_ascii_uppercase();
                        self.windows.push(TidWindow { tid, start: begin.cast_signed(), bases });
                    }
                    SliceReference::PerTid(&self.windows)
                }
            };

            let decoded = slice::decode_slice(
                &ch,
                data,
                offset,
                reference,
                header,
                read_group_ids,
                Keep::All,
                self.decode_md,
                &mut self.scratch,
            )?;

            // r[impl cram.stream.raw_records]
            // r[impl cram.stream.fidelity+2]
            let first = u64::try_from(sh.record_counter).unwrap_or(0);
            out.reserve(decoded.alignments().len());
            for (i, fields) in (0u64..).zip(decoded.alignments()) {
                let record = first.saturating_add(i);
                let mut payload = Vec::new();
                encode_alignment(&mut payload, &fields)
                    .map_err(|source| CramStreamError::RecordEncode { record, source })?;
                let raw = RawBamRecord::from_payload(payload, header)
                    .map_err(|source| CramStreamError::InvalidRecord { record, source })?;
                out.push(raw);
            }
        }
        Ok(())
    }
}

/// Fetch `begin..end` (0-based, half-open) of `tid` from the FASTA into
/// `buf`. `false` when the FASTA lacks the contig; an empty span fetches
/// nothing and succeeds.
fn fetch(
    fasta: &mut IndexedFastaReader,
    header: &BamHeader,
    tid: Tid,
    begin: u64,
    end: u64,
    buf: &mut Vec<u8>,
) -> Result<bool, CramError> {
    buf.clear();
    if begin >= end {
        return Ok(true);
    }
    let name = header.target_name(tid).ok_or(CramError::UnknownTid { tid: tid.as_u32() })?;
    let position = |value: u64| {
        Pos0::try_from(value).map_err(|_| CramError::InvalidPosition {
            value: i64::try_from(value).unwrap_or(i64::MAX),
        })
    };
    let span = RangeInclusive { start: position(begin)?, last: position(end.saturating_sub(1))? };
    // r[impl cram.edge.missing_reference+2]
    match fasta.fetch_seq_into(name, span, buf) {
        Ok(()) => Ok(true),
        Err(FastaError::SequenceNotFound { .. }) => Ok(false),
        Err(e) => Err(e.into()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn containers_and_decoders_are_send() {
        fn send<T: Send>() {}
        send::<CramContainer>();
        send::<CramContainerDecoder>();
    }
}
