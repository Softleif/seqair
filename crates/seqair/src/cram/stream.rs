//! Read a whole CRAM front to back, without an index and without seeking.
//! [`CramStreamReader`] reads containers off the input; a
//! [`CramContainerDecoder`] per thread turns them into [`RawBamRecord`]s.
//! See `r[cram.stream+2]`.

use super::{
    block,
    compression_header::CompressionHeader,
    container::{ContainerHeader, RefSeqId},
    reader::{self, CramError},
    slice::{
        self, DecodedSlice, Keep, RefBases, SliceHeader, SliceReference, SliceScratch, Synthesis,
        TidWindow,
    },
    varint,
};
use crate::bam::{BamHeader, DecodeError, RawBamRecord};
use crate::fasta::{FastaError, IndexedFastaReader};
use crate::reader::Tid;
use crate::sam::encode::{RecordEncodeError, encode_alignment};
use core::range::RangeInclusive;
use seqair_types::{Pos0, SmolStr};
use std::io::{self, Read, Seek, SeekFrom};
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

    /// The file ends inside the container a region query's index entry
    /// points at, at this byte offset (`r[cram.query_raw]`).
    #[error("CRAM file ends inside the container at byte {offset}")]
    TruncatedContainerAt { offset: u64 },

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

// r[impl cram.stream+2]
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
    pub(crate) header: ContainerHeader,
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
        let truncated = || CramStreamError::TruncatedContainer { container: 0 };
        let mut data = Vec::new();
        read_container(&mut inner, &mut data, truncated)?.ok_or_else(truncated)?;
        let header = reader::header_from_container_data(&data)?;
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
        Ok(CramContainerDecoder { fasta, slices: SliceDecoder::new(Arc::clone(&self.shared)) })
    }

    // r[impl cram.stream.framing+2]
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
            let mut data = Vec::new();
            let read = read_container(&mut self.inner, &mut data, || {
                CramStreamError::TruncatedContainer { container: number }
            });
            let Some(header) = read.inspect_err(|_| self.state = State::Failed)? else {
                // r[impl cram.file.eof+2]
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

// r[impl cram.query_raw]
// r[impl cram.stream.container_header+2]
/// Read the container at byte `offset` of `file` as the stream reads one,
/// its body into `data`. `None` for the EOF container.
pub(crate) fn read_container_at<R: Read + Seek>(
    file: &mut R,
    offset: u64,
    data: &mut Vec<u8>,
) -> Result<Option<ContainerHeader>, CramStreamError> {
    let truncated = || CramStreamError::TruncatedContainerAt { offset };
    file.seek(SeekFrom::Start(offset)).map_err(|source| CramStreamError::Io { source })?;
    // One read for the header, which `read_container` takes a field at a
    // time: a few hundred bytes cover all but containers of many slices,
    // whose rest comes from the file.
    let mut prefix = [0u8; 512];
    let n = loop {
        match file.read(&mut prefix) {
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            read => break read.map_err(|source| CramStreamError::Io { source })?,
        }
    };
    let mut input = prefix.get(..n).unwrap_or_default().chain(file);
    let header = read_container(&mut input, data, truncated)?.ok_or_else(truncated)?;
    Ok((!header.is_eof()).then_some(header))
}

// r[impl cram.stream.container_header+2]
/// Read one container from `inner`: its header, and its body into `data`
/// (cleared first). `None` when the input ends before its first byte;
/// `truncated` is the error for input that ends inside it.
pub(crate) fn read_container<R: Read>(
    inner: &mut R,
    data: &mut Vec<u8>,
    truncated: impl Fn() -> CramStreamError,
) -> Result<Option<ContainerHeader>, CramStreamError> {
    let truncated_by = |e: io::Error| match e.kind() {
        io::ErrorKind::UnexpectedEof => truncated(),
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
        r.read_exact(&mut rest_of_length).map_err(&truncated_by)?;
        // reference id, alignment start, span, record count
        for _ in 0..4 {
            varint::read_itf8(&mut r).map_err(&truncated_by)?;
        }
        // record counter, bases
        for _ in 0..2 {
            varint::read_ltf8(&mut r).map_err(&truncated_by)?;
        }
        varint::read_itf8(&mut r).map_err(&truncated_by)?; // block count
        let landmarks = varint::read_itf8(&mut r).map_err(&truncated_by)?.cast_signed();
        let count = usize::try_from(landmarks)
            .map_err(|_| CramError::InvalidLength { value: landmarks })?;
        reader::check_alloc_size(count.saturating_mul(4), "container landmark array")?;
        for _ in 0..count {
            varint::read_itf8(&mut r).map_err(&truncated_by)?;
        }
        let mut crc = [0u8; 4];
        r.read_exact(&mut crc).map_err(&truncated_by)?;
    }
    let header = ContainerHeader::parse(&bytes)?;

    let len = usize::try_from(header.length)
        .map_err(|_| CramError::InvalidLength { value: header.length })?;
    reader::check_alloc_size(len, "container data")?;
    // Grown as the bytes arrive, so a corrupt length cannot allocate ahead.
    data.clear();
    let got = inner
        .take(u64::try_from(len).unwrap_or(u64::MAX))
        .read_to_end(data)
        .map_err(|source| CramStreamError::Io { source })?;
    if got != len {
        return Err(truncated());
    }
    Ok(Some(header))
}

// r[impl cram.stream.decoder]
/// Decodes [`CramContainer`]s of one file into [`RawBamRecord`]s. Keep one
/// per thread ([`Self::fork`]); its buffers and codec tables are reused from
/// container to container.
pub struct CramContainerDecoder {
    fasta: IndexedFastaReader,
    slices: SliceDecoder,
}

impl std::fmt::Debug for CramContainerDecoder {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("CramContainerDecoder").finish_non_exhaustive()
    }
}

impl CramContainerDecoder {
    /// A decoder for another thread, sharing the header and the FASTA index.
    ///
    /// # Errors
    /// The FASTA cannot be reopened.
    pub fn fork(&self) -> Result<Self, CramStreamError> {
        let fasta = self.fasta.fork().map_err(CramError::from)?;
        let slices = SliceDecoder {
            decode_md: self.slices.decode_md,
            name_prefix: self.slices.name_prefix.clone(),
            ..SliceDecoder::new(Arc::clone(&self.slices.shared))
        };
        Ok(Self { fasta, slices })
    }

    // r[impl cram.stream.md_nm.option]
    /// Whether mapped records that store no `MD`/`NM` get them generated, as
    /// htslib does by default (`r[cram.stream.md_nm]`); `false` is
    /// `samtools view --input-fmt-option decode_md=0`. On by default.
    #[must_use]
    pub fn decode_md(mut self, on: bool) -> Self {
        self.slices.decode_md = on;
        self
    }

    // r[impl cram.record.read_name+2]
    /// The prefix of the names generated for records a `lossy_names` file
    /// stores none for, `<prefix>:<record number>`. htslib uses the basename
    /// of the path it opened (`-` for stdin, the default here); samtools'
    /// `--input-fmt-option name_prefix=` sets it.
    #[must_use]
    pub fn name_prefix(mut self, prefix: impl Into<SmolStr>) -> Self {
        self.slices.name_prefix = prefix.into();
        self
    }

    /// The file header records are validated against.
    pub fn header(&self) -> &BamHeader {
        self.slices.header()
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
        let ch = self.slices.compression_header(data)?;
        for &landmark in &container.header.landmarks {
            let offset = usize::try_from(landmark)
                .map_err(|_| CramError::InvalidLength { value: landmark })?;
            let decoded =
                self.slices.decode_slice(&mut self.fasta, &ch, data, offset, Keep::All)?;
            let header = self.slices.header();

            // r[impl cram.stream.raw_records]
            // r[impl cram.stream.fidelity+2]
            out.reserve(decoded.records.alignments().len());
            for (i, fields) in decoded.records.alignments().enumerate() {
                let record = decoded.number(i);
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

/// One slice's records, decoded by [`SliceDecoder::decode_slice`].
#[derive(Debug, Default)]
pub(crate) struct StreamSlice {
    pub(crate) records: DecodedSlice,
    /// The slice's `record_counter`: the number of its first record.
    first: u64,
    /// Each decoded record's index among all of the slice's records.
    kept: Vec<u32>,
}

impl StreamSlice {
    /// The 0-based number in the file of decoded record `i`, for errors.
    pub(crate) fn number(&self, i: usize) -> u64 {
        let index = self.kept.get(i).copied().unwrap_or_default();
        self.first.saturating_add(u64::from(index))
    }
}

// r[impl cram.stream.slice_decode]
/// The state decoding a slice needs besides the FASTA, which the caller
/// lends: a [`CramContainerDecoder`] owns its own, the indexed reader's raw
/// query lends the reader's.
pub(crate) struct SliceDecoder {
    shared: Arc<StreamShared>,
    scratch: SliceScratch,
    ref_buf: Vec<u8>,
    windows: Vec<TidWindow>,
    decode_md: bool,
    name_prefix: SmolStr,
}

impl SliceDecoder {
    fn new(shared: Arc<StreamShared>) -> Self {
        Self {
            shared,
            scratch: SliceScratch::default(),
            ref_buf: Vec::new(),
            windows: Vec::new(),
            decode_md: true,
            name_prefix: SmolStr::new_static("-"),
        }
    }

    /// A decoder for the indexed reader's raw queries, with htslib's
    /// defaults: MD/NM generated, names prefixed `name_prefix`.
    pub(crate) fn for_query(
        header: BamHeader,
        read_group_ids: Vec<SmolStr>,
        name_prefix: SmolStr,
    ) -> Self {
        Self { name_prefix, ..Self::new(Arc::new(StreamShared { header, read_group_ids })) }
    }

    pub(crate) fn header(&self) -> &BamHeader {
        &self.shared.header
    }

    /// The compression header, the first block of a container's `data`.
    pub(crate) fn compression_header(
        &mut self,
        data: &[u8],
    ) -> Result<CompressionHeader, CramError> {
        let (comp_block, _) = block::parse_block_with(data, &mut self.scratch.codec)?;
        if comp_block.content_type != block::ContentType::CompressionHeader {
            return Err(CramError::ExpectedCompressionHeader { found: comp_block.content_type });
        }
        CompressionHeader::parse(&comp_block.data)
    }

    /// Decode the records `keep` names of the slice at `offset` in a
    /// container's `data`, against the reference `fasta` holds for it.
    pub(crate) fn decode_slice<F: Read + Seek>(
        &mut self,
        fasta: &mut IndexedFastaReader<F>,
        ch: &CompressionHeader,
        data: &[u8],
        offset: usize,
        keep: Keep,
    ) -> Result<StreamSlice, CramError> {
        let header = &self.shared.header;
        let read_group_ids = self.shared.read_group_ids.as_slice();
        let slice_data =
            data.get(offset..).ok_or(CramError::Truncated { context: "slice offset" })?;
        let (sh_block, _) = block::parse_block_with(slice_data, &mut self.scratch.codec)?;
        if sh_block.content_type != block::ContentType::SliceHeader {
            return Err(CramError::ExpectedSliceHeader { found: sh_block.content_type });
        }
        let sh = SliceHeader::parse(&sh_block.data)?;
        let slice_ref = sh.ref_seq_id.resolve(header)?;

        // r[impl cram.stream.reference_window+2]
        // A file written without a reference (`RR` false) is decoded
        // without one (`r[cram.record.no_ref_fill+1]`): nothing to fetch.
        let no_ref = !ch.preservation.reference_required;
        let reference = match slice_ref {
            RefSeqId::Unmapped => SliceReference::Window { bases: Some(&[]), start: 0 },
            _ if no_ref => SliceReference::Window { bases: Some(&[]), start: 0 },
            RefSeqId::Ref(tid) => {
                let begin = u64::try_from(sh.alignment_start.max(1)).unwrap_or(1).saturating_sub(1);
                let span = u64::try_from(sh.alignment_span).unwrap_or(0);
                let end = begin.saturating_add(span).min(header.target_len(tid).unwrap_or(0));
                let found = fetch(fasta, header, tid, begin, end, &mut self.ref_buf)?;
                SliceReference::Window {
                    bases: found.then_some(self.ref_buf.as_slice()),
                    start: begin.cast_signed(),
                }
            }
            RefSeqId::Multi => {
                let spans = slice::multi_ref_spans(
                    ch,
                    data,
                    offset,
                    header,
                    read_group_ids,
                    &mut self.scratch,
                )?;
                self.windows.clear();
                // r[impl cram.query_raw.htslib]
                // Only the references whose records are kept need theirs.
                for (tid, span) in spans.into_iter().filter(|(tid, _)| keep.keeps(Some(*tid))) {
                    let begin = u64::try_from(span.start).unwrap_or(0);
                    let end = u64::try_from(span.end)
                        .unwrap_or(0)
                        .min(header.target_len(tid).unwrap_or(0));
                    if !fetch(fasta, header, tid, begin, end, &mut self.ref_buf)? {
                        return Err(CramError::MissingReference {
                            contig: header.target_name(tid).unwrap_or("?").into(),
                        });
                    }
                    let bases = RefBases::new(&self.ref_buf);
                    self.windows.push(TidWindow { tid, start: begin.cast_signed(), bases });
                }
                SliceReference::PerTid(&self.windows)
            }
        };

        let (records, indices) = slice::decode_slice_indexed(
            ch,
            data,
            offset,
            reference,
            header,
            read_group_ids,
            keep,
            Synthesis { md_nm: self.decode_md, name_prefix: Some(&self.name_prefix) },
            &mut self.scratch,
        )?;
        debug_assert_eq!(indices.len(), records.alignments().len());
        let first = u64::try_from(sh.record_counter).unwrap_or(0);
        Ok(StreamSlice { records, first, kept: indices })
    }
}

/// Fetch `begin..end` (0-based, half-open) of `tid` from the FASTA into
/// `buf`. `false` when the FASTA lacks the contig; an empty span fetches
/// nothing and succeeds.
fn fetch<F: Read + Seek>(
    fasta: &mut IndexedFastaReader<F>,
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
        // r[impl cram.stream.reference_window+2]
        Err(FastaError::RegionOutOfBounds { seq_len, .. }) => Err(CramError::ReferenceTooShort {
            contig: name.into(),
            header_len: header.target_len(tid).unwrap_or(0),
            fasta_len: seq_len,
        }),
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
