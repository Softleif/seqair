//! Open and query bgzf-compressed SAM files. [`IndexedSamReader`] uses a tabix index to fetch
//! records for a region into a [`RecordStore`], with the same API as [`crate::bam::reader::IndexedBamReader`].

use super::encode::{self, SamRecordEncoder, SamRecordError};
use crate::bam::{
    BaiError, BamHeader, BamHeaderError, BamIndex, MateFields, RawBamRecord,
    bgzf::VirtualOffset,
    cigar::{self, CigarOp},
    flags::BamFlags,
    index::Chunk,
    record::DecodeError,
    record_store::RecordStore,
    region_buf::{BlockCache, RegionBuf},
};
use crate::io::text_scan;
use crate::reader::{ResolveTid, Tid};
use core::range::RangeInclusive;
use seqair_types::{Base, MappingQuality, Pos0, Pos1};
use std::num::NonZeroI32;
use std::{
    convert::Infallible,
    fs::File,
    io::{Read, Seek},
    path::{Path, PathBuf},
    sync::Arc,
};
use tracing::instrument;

use super::super::bam::record_store::CustomizeRecordStore;
use crate::reader::FetchedRecord;

#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum SamError {
    #[error("I/O error opening {path}")]
    Open { path: PathBuf, source: std::io::Error },

    #[error("SAM header error")]
    Header {
        #[from]
        source: BamHeaderError,
    },

    #[error("tabix index error")]
    Index {
        #[from]
        source: BaiError,
    },

    #[error("tabix index not found for {sam_path} (tried .tbi, .bai and .csi)")]
    IndexNotFound { sam_path: PathBuf },

    // r[impl sam.index.csi+2]
    #[error(
        "found CSI index at {path} but CSI indexes are not yet supported. \
         Re-index with `tabix -p sam` to create a .tbi index instead.",
        path = path.display()
    )]
    CsiNotSupported { path: PathBuf },

    #[error(transparent)]
    MalformedRecord {
        #[from]
        source: SamRecordError,
    },

    #[error("BGZF error")]
    Bgzf {
        #[from]
        source: crate::bam::BgzfError,
    },

    #[error(
        "plain (uncompressed) SAM cannot be indexed. \
         Compress with `bgzip {path}` then index with `tabix -p sam {path}.gz`.",
        path = path.display()
    )]
    UncompressedSam { path: PathBuf },

    #[error("record decode error")]
    RecordDecode {
        #[from]
        source: DecodeError,
    },

    /// An earlier error ended this raw query (`r[sam.reader.query_raw]`).
    #[error("SAM query already failed; an earlier read returned the error")]
    Poisoned,
}

pub struct SamShared {
    index: BamIndex,
    header: BamHeader,
    sam_path: PathBuf,
}

impl SamShared {
    pub fn index(&self) -> &BamIndex {
        &self.index
    }
}

pub struct IndexedSamReader<R = File> {
    bulk_reader: R,
    shared: Arc<SamShared>,
    /// Decompressed blocks kept between queries; see `r[region_buf.block_cache]`.
    block_cache: BlockCache,
}

impl<R> std::fmt::Debug for IndexedSamReader<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("IndexedSamReader").field("sam_path", &self.shared.sam_path).finish()
    }
}

impl IndexedSamReader<File> {
    // r[impl sam.reader.open]
    /// Open a bgzf-compressed SAM and the index beside it (`r[index.locate]`).
    #[instrument(level = "debug", fields(path = %path.display()))]
    pub fn open(path: &Path) -> Result<Self, SamError> {
        Self::open_inner(path, None)
    }

    // r[impl index.explicit]
    /// Open a bgzf-compressed SAM with the tabix index or BAI at `index`,
    /// wherever it is; nothing beside the SAM is looked at.
    #[instrument(level = "debug", fields(path = %path.display(), index = %index.as_ref().display()))]
    pub fn open_with_index(
        path: &Path,
        index: &crate::reader::IndexPath,
    ) -> Result<Self, SamError> {
        Self::open_inner(path, Some(index.as_ref()))
    }

    fn open_inner(path: &Path, index: Option<&Path>) -> Result<Self, SamError> {
        // r[impl sam.reader.reject_uncompressed]
        {
            let mut f = File::open(path)
                .map_err(|source| SamError::Open { path: path.to_path_buf(), source })?;
            let mut magic = [0u8; 1];
            if f.read_exact(&mut magic).is_ok() && magic[0] == b'@' {
                return Err(SamError::UncompressedSam { path: path.to_path_buf() });
            }
        }

        // Parse SAM header from BGZF stream
        let mut bgzf = crate::bam::bgzf::BgzfReader::open(path)?;
        let header_text = read_sam_header(&mut bgzf)?;
        let header = BamHeader::from_sam_text(&header_text)?;
        // r[impl unified.sort_order]
        header.validate_sort_order()?;

        let index = match index {
            Some(index) => open_index(index, &header)?,
            None => find_tabix_index(path, &header)?,
        };

        let bulk_file = File::open(path)
            .map_err(|source| SamError::Open { path: path.to_path_buf(), source })?;

        Ok(IndexedSamReader {
            bulk_reader: bulk_file,
            shared: Arc::new(SamShared { index, header, sam_path: path.to_path_buf() }),
            block_cache: BlockCache::new(),
        })
    }

    pub fn fork(&self) -> Result<Self, SamError> {
        let bulk_file = File::open(&self.shared.sam_path)
            .map_err(|source| SamError::Open { path: self.shared.sam_path.clone(), source })?;

        Ok(IndexedSamReader {
            bulk_reader: bulk_file,
            shared: Arc::clone(&self.shared),
            block_cache: BlockCache::new(),
        })
    }
}

#[cfg(feature = "fuzz")]
impl IndexedSamReader<std::io::Cursor<Vec<u8>>> {
    /// Build an in-memory BGZF SAM reader from raw bytes (BGZF-compressed SAM + tabix index).
    pub fn from_bytes(sam_data: Vec<u8>, tbi_data: &[u8]) -> Result<Self, SamError> {
        let mut bgzf = crate::bam::bgzf::BgzfReader::from_cursor(sam_data.clone());
        let header_text = read_sam_header(&mut bgzf)?;
        let header = BamHeader::from_sam_text(&header_text)?;
        header.validate_sort_order()?;
        let index = BamIndex::from_tabix_bytes(tbi_data, &header)?;

        Ok(IndexedSamReader {
            bulk_reader: std::io::Cursor::new(sam_data),
            shared: Arc::new(SamShared { index, header, sam_path: PathBuf::from("<fuzz>") }),
            block_cache: BlockCache::new(),
        })
    }

    /// Build an in-memory plain (uncompressed) SAM reader for fuzzing.
    /// Parses header from the text, then linearly scans all records into a `RecordStore`
    /// on each `fetch_plain_into` call.
    pub fn from_plain_bytes(sam_data: Vec<u8>) -> Result<Self, SamError> {
        let text = String::from_utf8(sam_data.clone())
            .map_err(|source| SamRecordError::HeaderNotUtf8 { source })?;

        // Extract header lines (starting with @)
        let mut header_text = String::new();
        for line in text.lines() {
            if line.starts_with('@') {
                header_text.push_str(line);
                header_text.push('\n');
            } else {
                break;
            }
        }

        let header = BamHeader::from_sam_text(&header_text)?;

        // No index needed — we'll do linear scans
        // Create a dummy empty index
        let index = BamIndex::empty();

        Ok(IndexedSamReader {
            bulk_reader: std::io::Cursor::new(sam_data),
            shared: Arc::new(SamShared { index, header, sam_path: PathBuf::from("<fuzz-plain>") }),
            block_cache: BlockCache::new(),
        })
    }

    /// Linear scan of plain (uncompressed) SAM text, parsing all records that overlap [start, end).
    pub fn fetch_plain_into(
        &mut self,
        tid: Tid,
        span: RangeInclusive<Pos0>,
        store: &mut RecordStore,
    ) -> Result<usize, SamError> {
        store.clear();
        // r[impl interval.empty_span]
        if span.is_empty() {
            return Ok(0);
        }

        self.bulk_reader.set_position(0);
        let mut all_data = Vec::new();
        std::io::Read::read_to_end(&mut self.bulk_reader, &mut all_data)
            .map_err(|source| SamError::Open { path: PathBuf::from("<fuzz-plain>"), source })?;

        let start_i64 = span.start.as_i64();
        let end_i64 = span.last.as_i64();

        let mut cigar_buf = Vec::with_capacity(256);
        let mut bases_buf = Vec::with_capacity(256);
        let mut qual_buf = Vec::with_capacity(256);
        let mut aux_buf = Vec::with_capacity(256);

        for line in all_data.split(|&b| b == b'\n') {
            if line.is_empty() || line.first() == Some(&b'@') {
                continue;
            }
            // Strip \r
            let line = if line.last() == Some(&b'\r') {
                line.get(..line.len().saturating_sub(1)).unwrap_or(line)
            } else {
                line
            };
            if line.is_empty() {
                continue;
            }
            parse_sam_line(
                line,
                &self.shared.header,
                tid,
                start_i64,
                end_i64,
                store,
                &mut cigar_buf,
                &mut bases_buf,
                &mut qual_buf,
                &mut aux_buf,
                &mut (),
            )?;
        }

        Ok(store.len())
    }
}

impl<R: Read + Seek> IndexedSamReader<R> {
    pub fn shared(&self) -> &Arc<SamShared> {
        &self.shared
    }

    pub fn header(&self) -> &BamHeader {
        &self.shared.header
    }

    /// Estimate the compressed bytes a `[start, end]` region query would load,
    /// for byte-aware segmentation. See
    /// [`IndexedBamReader::estimate_region_bytes`](crate::bam::IndexedBamReader::estimate_region_bytes).
    pub fn estimate_region_bytes(&self, tid: Tid, span: RangeInclusive<Pos0>) -> u64 {
        let chunks = self.shared.index.query(tid.as_u32(), span.start, span.last);
        crate::bam::region_buf::merged_byte_size(&chunks) as u64
    }

    // r[impl sam.reader.fetch_into]
    // r[impl sam.reader.sorted_order]
    #[instrument(level = "debug", skip(self, store))]
    // r[impl sam.perf.bulk_read]
    pub fn fetch_into(
        &mut self,
        tid: Tid,
        span: RangeInclusive<Pos0>,
        store: &mut RecordStore,
    ) -> Result<usize, SamError> {
        self.fetch_into_customized(tid, span, store, &mut ()).map(|c| c.kept)
    }

    // r[impl unified.fetch_into_customized]
    pub fn fetch_into_customized<E: CustomizeRecordStore>(
        &mut self,
        tid: Tid,
        span: RangeInclusive<Pos0>,
        store: &mut RecordStore<E::Extra>,
        customize: &mut E,
    ) -> Result<crate::reader::FetchCounts, SamError> {
        store.clear();

        // r[impl interval.empty_span]
        if span.is_empty() {
            return Ok(crate::reader::FetchCounts::default());
        }
        let chunks = self.shared.index.query(tid.as_u32(), span.start, span.last);
        if chunks.is_empty() {
            return Ok(crate::reader::FetchCounts::default());
        }

        let region = RegionBuf::with_cache(&mut self.bulk_reader, &chunks, &mut self.block_cache)?;
        let mut lines = ChunkLines::new(region, chunks);

        let start_i64 = span.start.as_i64();
        let end_i64 = span.last.as_i64();
        let header = &self.shared.header;

        // Scratch buffers reused across records
        let mut cigar_buf = Vec::with_capacity(256);
        let mut seq_buf = Vec::with_capacity(256);
        let mut qual_buf = Vec::with_capacity(256);
        let mut aux_buf = Vec::with_capacity(256);
        let mut counts = crate::reader::FetchCounts::default();

        let None::<Infallible> = lines.next_line(|line| {
            if let Some(record) = parse_sam_line(
                line,
                header,
                tid,
                start_i64,
                end_i64,
                store,
                &mut cigar_buf,
                &mut seq_buf,
                &mut qual_buf,
                &mut aux_buf,
                customize,
            )? {
                counts.add(record);
            }
            Ok(None)
        })?;

        Ok(counts)
    }

    // r[impl sam.reader.query_raw]
    /// The lines overlapping `span` on `tid`, encoded as `SamStreamReader`
    /// encodes them: the selection and order of [`Self::fetch_into`].
    pub fn query_raw(
        &mut self,
        tid: Tid,
        span: RangeInclusive<Pos0>,
    ) -> Result<SamRawQuery<'_, R>, SamError> {
        // r[impl interval.empty_span]
        let chunks = if span.is_empty() {
            Vec::new()
        } else {
            self.shared.index.query(tid.as_u32(), span.start, span.last)
        };
        let region = RegionBuf::with_cache(&mut self.bulk_reader, &chunks, &mut self.block_cache)?;
        Ok(SamRawQuery {
            lines: ChunkLines::new(region, chunks),
            shared: Arc::clone(&self.shared),
            encoder: SamRecordEncoder::new(),
            scratch: RawBamRecord::default(),
            cigar: Vec::new(),
            tid,
            span,
            failed: false,
        })
    }
}

/// The records of one region of a SAM as [`RawBamRecord`]s; see
/// [`IndexedSamReader::query_raw`]. Borrows the reader until dropped.
pub struct SamRawQuery<'r, R: Read + Seek> {
    lines: ChunkLines<'r, R>,
    shared: Arc<SamShared>,
    encoder: SamRecordEncoder,
    /// Each line is encoded here, and swapped out once it is kept, so a
    /// query that ends leaves the caller's record as it was.
    scratch: RawBamRecord,
    /// For [`select_line`] on a line the encoder rejects.
    cigar: Vec<CigarOp>,
    tid: Tid,
    span: RangeInclusive<Pos0>,
    failed: bool,
}

impl<R: Read + Seek> std::fmt::Debug for SamRawQuery<'_, R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SamRawQuery")
            .field("sam_path", &self.shared.sam_path)
            .field("tid", &self.tid)
            .field("failed", &self.failed)
            .finish_non_exhaustive()
    }
}

impl<R: Read + Seek> SamRawQuery<'_, R> {
    pub fn header(&self) -> &BamHeader {
        &self.shared.header
    }

    /// Read the next record into `record`, reusing its buffer. `Ok(false)` at
    /// the end of the region, which leaves `record` unchanged.
    ///
    /// # Errors
    /// [`SamError::MalformedRecord`] for a kept line that does not encode —
    /// the query stays usable. Anything else ends the query: later calls
    /// return [`SamError::Poisoned`].
    pub fn read_record(&mut self, record: &mut RawBamRecord) -> Result<bool, SamError> {
        if self.failed {
            return Err(SamError::Poisoned);
        }
        let Self { lines, shared, encoder, scratch, cigar, tid, span, .. } = self;
        let header = &shared.header;
        let result = lines.next_line(|line| {
            if line_reference(line, header)? != Some(*tid) {
                return Ok(None);
            }
            match encoder.encode(line, header, scratch) {
                Ok(()) => Ok(overlaps(scratch, *span).then_some(())),
                // Without a record there is no overlap to test: the line is
                // an error only where `fetch_into` selects it.
                Err(e) => {
                    let (fields, n_fields) = text_scan::split_tabs(line);
                    let fields = fields.get(..n_fields).unwrap_or_default();
                    let (start, end) = (span.start.as_i64(), span.last.as_i64());
                    match select_line(fields, header, *tid, start, end, cigar)? {
                        Some(_) => Err(e.into()),
                        None => Ok(None),
                    }
                }
            }
        });
        match result {
            Ok(Some(())) => {
                std::mem::swap(record, &mut self.scratch);
                Ok(true)
            }
            Ok(None) => Ok(false),
            Err(e) => {
                if !matches!(e, SamError::MalformedRecord { .. }) {
                    self.failed = true;
                }
                Err(e)
            }
        }
    }
}

impl<R: Read + Seek> Iterator for SamRawQuery<'_, R> {
    type Item = Result<RawBamRecord, SamError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.failed {
            return None;
        }
        let mut record = RawBamRecord::default();
        match self.read_record(&mut record) {
            Ok(true) => Some(Ok(record)),
            Ok(false) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

/// The reference a line's RNAME names, checked as [`parse_sam_line`] checks
/// it: `None` for `*` or a name the header lacks, which the query skips.
fn line_reference(line: &[u8], header: &BamHeader) -> Result<Option<Tid>, SamError> {
    let (fields, n_fields) = text_scan::split_tabs(line);
    if n_fields < 11 {
        return Err(SamRecordError::TooFewFields { found: n_fields }.into());
    }
    let [_, _, rname, ..] = fields;
    if rname == b"*" {
        return Ok(None);
    }
    let rname = std::str::from_utf8(rname)
        .map_err(|_| SamRecordError::InvalidRname { value: rname.into() })?;
    Ok(rname.resolve_tid(header).ok())
}

/// `r[sam.reader.overlap_filter+4]` on an encoded record: an unmapped one
/// sits at its POS whatever its CIGAR.
fn overlaps(record: &RawBamRecord, span: RangeInclusive<Pos0>) -> bool {
    let Some(pos) = record.pos() else { return false };
    let last = if record.flags().is_unmapped() { pos } else { record.end_pos().unwrap_or(pos) };
    pos <= span.last && last >= span.start
}

/// The alignment lines of a query's index chunks, handed out one at a time
/// and resumable between them.
struct ChunkLines<'r, R: Read + Seek> {
    region: RegionBuf<'r, R>,
    chunks: Vec<Chunk>,
    next_chunk: usize,
    /// End of the chunk being read; `None` between chunks.
    chunk_end: Option<VirtualOffset>,
    /// A line spanning BGZF blocks, assembled; a line inside one block is
    /// handed out straight from the block.
    line: Vec<u8>,
}

impl<'r, R: Read + Seek> ChunkLines<'r, R> {
    fn new(region: RegionBuf<'r, R>, chunks: Vec<Chunk>) -> Self {
        Self { region, chunks, next_chunk: 0, chunk_end: None, line: Vec::with_capacity(1024) }
    }

    // r[impl sam.perf.text_parsing]
    // r[impl sam.edge.line_spanning_blocks]
    /// Hand each alignment line to `on_line` until it returns `Some`, and
    /// return that; `None` once every chunk is read.
    fn next_line<T>(
        &mut self,
        mut on_line: impl FnMut(&[u8]) -> Result<Option<T>, SamError>,
    ) -> Result<Option<T>, SamError> {
        loop {
            let Some(end) = self.chunk_end else {
                let Some(chunk) = self.chunks.get(self.next_chunk) else { return Ok(None) };
                self.next_chunk = self.next_chunk.saturating_add(1);
                self.region.seek_virtual(chunk.begin)?;
                self.line.clear();
                self.chunk_end = Some(chunk.end);
                continue;
            };
            // `fill_buf` may load the next block, so the chunk end is checked
            // again against where the returned bytes start: the end of one
            // block and the start of the next are one place.
            let block = if self.region.virtual_offset() >= end {
                None
            } else {
                let (voff, block) = self.region.fill_buf()?;
                (!block.is_empty() && voff < end).then_some((voff, block))
            };
            let Some((voff, block)) = block else {
                // The chunk is read; a last line may lack its newline.
                self.chunk_end = None;
                let found = alignment_line(&self.line).map(&mut on_line).transpose();
                self.line.clear();
                match found?.flatten() {
                    Some(t) => return Ok(Some(t)),
                    None => continue,
                }
            };
            let limit = if voff.block_offset() == end.block_offset() {
                usize::from(end.within_block().saturating_sub(voff.within_block()))
            } else {
                block.len()
            };
            let window = block.get(..limit).unwrap_or(block);
            let (consumed, found) = match text_scan::find_byte(window, b'\n') {
                Some(nl) => {
                    let rest = window.get(..nl).unwrap_or_default();
                    let found = if self.line.is_empty() {
                        alignment_line(rest).map(&mut on_line).transpose()
                    } else {
                        self.line.extend_from_slice(rest);
                        let found = alignment_line(&self.line).map(&mut on_line).transpose();
                        self.line.clear();
                        found
                    };
                    (nl.saturating_add(1), found)
                }
                None => {
                    self.line.extend_from_slice(window);
                    (window.len(), Ok(None))
                }
            };
            self.region.consume(consumed);
            if let Some(t) = found?.flatten() {
                return Ok(Some(t));
            }
        }
    }
}

// r[impl sam.edge.empty_lines]
/// `line` without a Windows `\r`, unless it is empty or a header line.
fn alignment_line(line: &[u8]) -> Option<&[u8]> {
    let line = line.strip_suffix(b"\r").unwrap_or(line);
    (!line.is_empty() && line.first() != Some(&b'@')).then_some(line)
}

/// Where a line [`parse_sam_line`] keeps sits, and what it parsed to find out.
struct SelectedLine {
    tid: Tid,
    pos: Pos0,
    end_pos: Pos0,
    flags: BamFlags,
    mapq: MappingQuality,
}

/// The fields `fetch_into` parses to select a line, in its order and with
/// its errors: `None` for a line it skips (RNAME `*`, unknown or another
/// reference, POS 0, or outside `start..=end`). Leaves the CIGAR in
/// `cigar_buf`.
fn select_line(
    fields: &[&[u8]],
    header: &BamHeader,
    tid_filter: Tid,
    start: i64,
    end: i64,
    cigar_buf: &mut Vec<CigarOp>,
) -> Result<Option<SelectedLine>, SamError> {
    if fields.len() < 11 {
        return Err(SamRecordError::TooFewFields { found: fields.len() }.into());
    }

    // Field 2: FLAG
    let flag_field = fields.get(1).copied().unwrap_or(b"0");
    let flags = BamFlags::from(
        parse_u16(flag_field)
            .ok_or_else(|| SamRecordError::InvalidFlag { value: flag_field.into() })?,
    );

    // r[impl sam.edge.rname_star]
    let rname = fields.get(2).copied().unwrap_or(b"*");
    if rname == b"*" {
        return Ok(None);
    }

    let rname_str = std::str::from_utf8(rname)
        .map_err(|_| SamRecordError::InvalidRname { value: rname.into() })?;
    // Unknown contig: skip.
    let Ok(rec_tid) = rname_str.resolve_tid(header) else { return Ok(None) };

    if rec_tid != tid_filter {
        return Ok(None);
    }

    // r[impl sam.record.coordinate_conversion]
    let pos_field = fields.get(3).copied().unwrap_or(b"0");
    let pos_1based = parse_i64(pos_field)
        .ok_or_else(|| SamRecordError::InvalidPos { value: pos_field.into() })?;
    // SAM POS is 1-based; POS=0 means unmapped/unplaced — skip it (mirrors the
    // BAM reader skipping pos = -1) rather than erroring the whole query. A
    // negative or out-of-range POS is still a malformed record and errors.
    if pos_1based == 0 {
        return Ok(None);
    }
    let pos = Pos1::try_from(pos_1based)
        .map(|p| p.to_zero_based())
        .map_err(|_| SamRecordError::InvalidPos { value: pos_field.into() })?;

    // Field 5: MAPQ
    let mapq_field = fields.get(4).copied().unwrap_or(b"0");
    let mapq = parse_u8(mapq_field)
        .map(MappingQuality::from_byte)
        .ok_or_else(|| SamRecordError::InvalidMapq { value: mapq_field.into() })?;

    // Field 6: CIGAR
    let cigar_str = fields.get(5).copied().unwrap_or(b"*");
    cigar_buf.clear();
    let cigar_available = parse_cigar(cigar_str, cigar_buf)?;

    // r[impl record_store.end_pos_htslib]
    // An unmapped record sits at `pos` whatever CIGAR it carries. Otherwise
    // the CIGAR's end; compute_end_pos returns None only if the ref span
    // overflows u32 range, so fall back to pos (zero-span) on malformed input.
    let end_pos = if cigar_available && !flags.is_unmapped() {
        cigar::compute_end_pos(pos, cigar_buf).unwrap_or(pos)
    } else {
        pos
    };

    // r[impl sam.reader.overlap_filter+4]
    // r[impl sam.reader.overlap_inclusive]
    // r[impl interval.overlap_test]
    if pos.as_i64() > end || end_pos.as_i64() < start {
        return Ok(None);
    }
    Ok(Some(SelectedLine { tid: rec_tid, pos, end_pos, flags, mapq }))
}

// r[impl sam.record.parse]
#[expect(
    clippy::too_many_arguments,
    reason = "SAM line parsing needs header, customize, region, store, and per-record output parameters"
)]
fn parse_sam_line<E: CustomizeRecordStore>(
    line: &[u8],
    header: &BamHeader,
    tid_filter: Tid,
    start: i64,
    end: i64,
    store: &mut RecordStore<E::Extra>,
    cigar_buf: &mut Vec<CigarOp>,
    seq_buf: &mut Vec<u8>,
    qual_buf: &mut Vec<u8>,
    aux_buf: &mut Vec<u8>,
    customize: &mut E,
) -> Result<Option<FetchedRecord>, SamError> {
    let (fields, n_fields) = text_scan::split_tabs(line);
    let fields = fields.get(..n_fields).unwrap_or_default();
    let Some(SelectedLine { tid: rec_tid, pos, end_pos, flags, mapq }) =
        select_line(fields, header, tid_filter, start, end, cigar_buf)?
    else {
        return Ok(None);
    };

    // Field 1: QNAME
    let qname = fields.first().copied().unwrap_or(b"*");

    // r[impl sam.record.seq_decode]
    // r[impl sam.edge.missing_seq]
    let seq_field = fields.get(9).copied().unwrap_or(b"*");
    seq_buf.clear();
    if seq_field != b"*" {
        seq_buf.extend_from_slice(seq_field);
    }
    let bases = Base::convert_ascii_in_place_as_slice(seq_buf);

    // r[impl sam.record.qual_decode]
    // r[impl sam.edge.missing_qual]
    let qual_field = fields.get(10).copied().unwrap_or(b"*");
    qual_buf.clear();
    if qual_field == b"*" {
        qual_buf.resize(bases.len(), 0xFF);
    } else {
        qual_buf.extend(qual_field.iter().map(|b| b.wrapping_sub(33)));
    }

    // Optional fields (12+): aux tags
    aux_buf.clear();
    if fields.len() > 11 {
        // The 12th element contains everything after the 11th TAB (because of splitn(12,...))
        let aux_text = fields.get(11).copied().unwrap_or(b"");
        // r[impl sam.record.aux_tags+2]
        // r[impl sam.record.aux_parse_strict+2]
        // r[impl sam.record.aux_int_range+3]
        encode::parse_aux_fields(aux_text, aux_buf)?;
    }

    // Field 7: RNEXT (mate reference name; * = unavailable, = = same as RNAME)
    let rnext_field = fields.get(6).copied().unwrap_or(b"*");
    // r[impl sam.record.mate_fields]
    let next_ref_id = if rnext_field == b"*" {
        None
    } else if rnext_field == b"=" {
        Some(rec_tid)
    } else {
        let rnext_str = std::str::from_utf8(rnext_field)
            .map_err(|_| SamRecordError::InvalidRname { value: rnext_field.into() })?;
        let next_ref_id = rnext_str.resolve_tid(header).ok();
        if next_ref_id.is_none() {
            encode::warn_unknown_mate_reference(rnext_field);
        }
        next_ref_id
    };

    // Field 8: PNEXT (1-based mate position, 0 = unavailable)
    let pnext_field = fields.get(7).copied().unwrap_or(b"0");
    let pnext = parse_i64(pnext_field)
        .ok_or_else(|| SamRecordError::InvalidPnext { value: pnext_field.into() })?;
    let next_pos = match pnext {
        0 => None,
        _ => Some(
            Pos1::try_from(pnext)
                .map_err(|_| SamRecordError::InvalidPnext { value: pnext_field.into() })?
                .to_zero_based(),
        ),
    };

    // Field 9: TLEN
    let tlen_field = fields.get(8).copied().unwrap_or(b"0");
    let template_len = parse_i64(tlen_field)
        .and_then(|t| i32::try_from(t).ok())
        .ok_or_else(|| SamRecordError::InvalidTlen { value: tlen_field.into() })?;
    let template_len = NonZeroI32::new(template_len);

    let (matching_bases, indel_bases) = cigar::calc_matches_indels(cigar_buf);

    let kept = store
        .push_fields(
            pos,
            end_pos,
            flags,
            mapq,
            matching_bases,
            indel_bases,
            qname,
            cigar_buf,
            bases,
            qual_buf,
            aux_buf,
            Some(rec_tid),
            MateFields { tid: next_ref_id, pos: next_pos, template_len },
            customize,
        )?
        .is_some();

    Ok(Some(FetchedRecord { flags, started: pos.as_i64() >= start, kept }))
}

// r[impl sam.record.cigar_parse+2]
// r[impl sam.edge.long_cigar]
/// Parse CIGAR into `buf`; `false` for `*`.
fn parse_cigar(cigar_str: &[u8], buf: &mut Vec<CigarOp>) -> Result<bool, SamRecordError> {
    encode::parse_cigar(cigar_str, buf)?;
    Ok(cigar_str != b"*")
}

// r[impl sam.header.parse]
// r[impl sam.header.sq_required]
// r[impl sam.header.preserve_full_text]
// r[impl sam.edge.line_spanning_blocks]
fn read_sam_header<R: std::io::Read + std::io::Seek>(
    bgzf: &mut crate::bam::bgzf::BgzfReader<R>,
) -> Result<String, SamError> {
    let mut header = Vec::with_capacity(4096);
    let mut line = Vec::with_capacity(256);

    while let Ok(b) = bgzf.read_byte() {
        if b == b'\n' {
            if line.first() == Some(&b'@') {
                header.extend_from_slice(&line);
                header.push(b'\n');
                line.clear();
            } else {
                // First non-header line — stop
                break;
            }
        } else {
            line.push(b);
        }
    }

    String::from_utf8(header).map_err(|source| SamRecordError::HeaderNotUtf8 { source }.into())
}

// r[impl sam.index.tabix]
// r[impl unified.detect_index+2]
// r[impl unified.detect_error]
// r[impl sam.index.locate+3]
fn find_tabix_index(sam_path: &Path, header: &BamHeader) -> Result<BamIndex, SamError> {
    let index = crate::utils::find_index(sam_path, &["tbi", "bai", "csi"])
        .ok_or_else(|| SamError::IndexNotFound { sam_path: sam_path.to_path_buf() })?;
    // r[depends sam.index.csi+2]
    if index.extension().is_some_and(|ext| ext == "csi") {
        return Err(SamError::CsiNotSupported { path: index });
    }
    open_index(&index, header)
}

/// Open a tabix index or a BAI, told apart by content: a BAI starts with its
/// magic, a tabix index is BGZF-compressed. A CSI is recognised and refused.
fn open_index(path: &Path, header: &BamHeader) -> Result<BamIndex, SamError> {
    let mut magic = [0u8; 4];
    let read = File::open(path).and_then(|mut f| f.read_exact(&mut magic));
    read.map_err(|source| BaiError::Read { path: path.to_path_buf(), source })?;
    if magic == *b"BAI\x01" {
        return Ok(BamIndex::from_path(path)?);
    }
    match BamIndex::from_tabix_path(path, header) {
        // r[depends sam.index.csi+2]
        Err(BaiError::InvalidTabixMagic { found }) if found == *b"CSI\x01" => {
            Err(SamError::CsiNotSupported { path: path.to_path_buf() })
        }
        result => Ok(result?),
    }
}

// --- Integer parsing from byte slices (no String allocation) ---

// r[impl sam.record.parse_int_nonempty]
fn parse_u8(bytes: &[u8]) -> Option<u8> {
    if bytes.is_empty() {
        return None;
    }
    let mut val = 0u16;
    for &b in bytes {
        if !b.is_ascii_digit() {
            return None;
        }
        val = val.checked_mul(10)?.checked_add(u16::from(b.checked_sub(b'0')?))?;
    }
    u8::try_from(val).ok()
}

fn parse_u16(bytes: &[u8]) -> Option<u16> {
    if bytes.is_empty() {
        return None;
    }
    let mut val = 0u32;
    for &b in bytes {
        if !b.is_ascii_digit() {
            return None;
        }
        val = val.checked_mul(10)?.checked_add(u32::from(b.checked_sub(b'0')?))?;
    }
    u16::try_from(val).ok()
}

fn parse_i64(bytes: &[u8]) -> Option<i64> {
    if bytes.is_empty() {
        return None;
    }
    let (negative, digits) =
        if bytes.first() == Some(&b'-') { (true, bytes.get(1..)?) } else { (false, bytes) };
    let mut val = 0i64;
    for &b in digits {
        if !b.is_ascii_digit() {
            return None;
        }
        val = val.checked_mul(10)?.checked_add(i64::from(b.checked_sub(b'0')?))?;
    }
    if negative { Some(val.checked_neg()?) } else { Some(val) }
}

#[cfg(test)]
#[allow(clippy::arithmetic_side_effects, reason = "test arithmetic is not safety-critical")]
mod tests {
    use super::*;
    use crate::bam::record_store::RecordStore;
    use crate::bam::test_util::tid;
    use std::io::Write;

    /// The aux bytes the shared parser produces for one `TAG:TYPE:VALUE` field.
    fn aux_bytes(field: &str) -> Result<Vec<u8>, SamError> {
        let mut buf = Vec::new();
        encode::parse_aux_fields(field.as_bytes(), &mut buf)?;
        Ok(buf)
    }

    // r[verify sam.record.aux_tags+2]
    /// A `B` array's subtype says whether its elements are signed. Reading a
    /// signed one through the unsigned type of the same width refuses every
    /// negative element, which is most of what a signed array is for.
    #[test]
    fn signed_b_arrays_keep_their_negative_elements() {
        assert_eq!(
            aux_bytes("Xa:B:c,-1,-128,127").unwrap(),
            b"XaBc\x03\x00\x00\x00\xff\x80\x7f".to_vec()
        );
        assert_eq!(
            aux_bytes("Xa:B:s,-1,-32768").unwrap(),
            b"XaBs\x02\x00\x00\x00\xff\xff\x00\x80".to_vec()
        );
        assert_eq!(
            aux_bytes("Xa:B:i,-1").unwrap(),
            b"XaBi\x01\x00\x00\x00\xff\xff\xff\xff".to_vec()
        );
        // The unsigned subtypes keep their full positive range.
        assert_eq!(aux_bytes("Xa:B:C,255").unwrap(), b"XaBC\x01\x00\x00\x00\xff".to_vec());
        assert_eq!(
            aux_bytes("Xa:B:I,4294967295").unwrap(),
            b"XaBI\x01\x00\x00\x00\xff\xff\xff\xff".to_vec()
        );
    }

    // r[verify sam.record.aux_tags+2]
    /// Out of range for the subtype, and a subtype that is not one at all: the
    /// element count is written before the elements, so neither may be skipped.
    #[test]
    fn b_arrays_reject_what_their_subtype_cannot_hold() {
        for field in ["Xa:B:c,128", "Xa:B:C,-1", "Xa:B:s,32768", "Xa:B:I,-1"] {
            let err = aux_bytes(field).unwrap_err();
            assert!(
                matches!(
                    err,
                    SamError::MalformedRecord { source: SamRecordError::InvalidAuxValue { .. } }
                ),
                "{field}: {err:?}"
            );
        }
        let err = aux_bytes("Xa:B:q,1").unwrap_err();
        assert!(
            matches!(
                err,
                SamError::MalformedRecord {
                    source: SamRecordError::InvalidAuxSubtype { subtype: b'q', .. }
                }
            ),
            "{err:?}"
        );
    }

    fn make_header() -> BamHeader {
        BamHeader::from_sam_text("@SQ\tSN:chr1\tLN:1000\n").expect("failed to build test header")
    }

    fn make_store_and_bufs() -> (RecordStore, Vec<CigarOp>, Vec<u8>, Vec<u8>, Vec<u8>) {
        (RecordStore::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new())
    }

    fn call_parse(line: &[u8], header: &BamHeader) -> Result<Option<bool>, SamError> {
        let (mut store, mut cigar_buf, mut bases_buf, mut qual_buf, mut aux_buf) =
            make_store_and_bufs();
        parse_sam_line(
            line,
            header,
            tid(0),
            0,
            10000,
            &mut store,
            &mut cigar_buf,
            &mut bases_buf,
            &mut qual_buf,
            &mut aux_buf,
            &mut (),
        )
        .map(|outcome| outcome.map(|record| record.kept))
    }

    #[test]
    fn error_too_few_fields() {
        let header = make_header();
        let line = b"read1\t0\tchr1";
        let err = call_parse(line, &header).expect_err("expected TooFewFields error");
        assert!(
            matches!(err, SamError::MalformedRecord { source: SamRecordError::TooFewFields { found } } if found < 11),
            "expected TooFewFields, got: {err}"
        );
    }

    #[test]
    fn error_invalid_flag() {
        let header = make_header();
        let line = b"r\tNOTANUM\tchr1\t100\t60\t10M\t*\t0\t0\tACGT\t~~~~";
        let err = call_parse(line, &header).expect_err("expected InvalidFlag error");
        assert!(
            matches!(err, SamError::MalformedRecord { source: SamRecordError::InvalidFlag { .. } }),
            "expected InvalidFlag, got: {err}"
        );
    }

    #[test]
    fn error_invalid_rname_not_utf8() {
        let header = make_header();
        // Build a line with invalid UTF-8 bytes in the RNAME field (field[2])
        let mut line = b"r\t0\t".to_vec();
        line.extend_from_slice(b"\xFF\xFE"); // invalid UTF-8
        line.extend_from_slice(b"\t100\t60\t10M\t*\t0\t0\tACGT\t~~~~");
        let err = call_parse(&line, &header).expect_err("expected InvalidRname error");
        assert!(
            matches!(
                err,
                SamError::MalformedRecord { source: SamRecordError::InvalidRname { .. } }
            ),
            "expected InvalidRname, got: {err}"
        );
    }

    // r[verify sam.record.mate_fields]
    /// RNEXT names resolve against the header (an unknown one is no mate
    /// reference), PNEXT 0 is no mate position.
    #[test]
    fn mate_fields_resolve_against_the_header() {
        let header =
            BamHeader::from_sam_text("@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:1000\n").unwrap();
        let tid = |n: u32| n.resolve_tid(&header).ok();
        for (rnext, pnext, want_ref, want_pos) in [
            ("*", "0", None, None),
            ("=", "50", tid(0), Pos0::new(49)),
            ("chr2", "1", tid(1), Pos0::new(0)),
            ("chr9", "7", None, Pos0::new(6)),
        ] {
            let line = format!("r\t1\tchr1\t100\t60\t4M\t{rnext}\t{pnext}\t0\tACGT\t~~~~");
            let (mut store, mut cigar_buf, mut bases_buf, mut qual_buf, mut aux_buf) =
                make_store_and_bufs();
            parse_sam_line(
                line.as_bytes(),
                &header,
                crate::bam::test_util::tid(0),
                0,
                10000,
                &mut store,
                &mut cigar_buf,
                &mut bases_buf,
                &mut qual_buf,
                &mut aux_buf,
                &mut (),
            )
            .unwrap()
            .unwrap();
            let rec = store.record(crate::bam::RecordIdx::ZERO).unwrap();
            assert_eq!((rec.mate_fields.tid, rec.mate_fields.pos), (want_ref, want_pos), "{line}");
        }
    }

    // r[verify sam.record.mate_fields]
    /// PNEXT and TLEN past `i32`, and a negative PNEXT, are malformed, not
    /// truncated into range.
    #[test]
    fn error_mate_fields_beyond_i32() {
        let header = make_header();
        for (line, want_pnext) in [
            (&b"r\t1\tchr1\t100\t60\t4M\t=\t2147483649\t0\tACGT\t~~~~"[..], true),
            (&b"r\t1\tchr1\t100\t60\t4M\t=\t-5\t0\tACGT\t~~~~"[..], true),
            (&b"r\t1\tchr1\t100\t60\t4M\t=\t200\t2147483648\tACGT\t~~~~"[..], false),
            (&b"r\t1\tchr1\t100\t60\t4M\t=\t200\t-2147483649\tACGT\t~~~~"[..], false),
        ] {
            let err = call_parse(line, &header).expect_err("mate field out of range");
            let SamError::MalformedRecord { source } = err else { panic!("{err}") };
            assert!(
                if want_pnext {
                    matches!(source, SamRecordError::InvalidPnext { .. })
                } else {
                    matches!(source, SamRecordError::InvalidTlen { .. })
                },
                "{source}"
            );
        }
    }

    #[test]
    fn error_invalid_pos() {
        let header = make_header();
        let line = b"r\t0\tchr1\tNOTANUM\t60\t10M\t*\t0\t0\tACGT\t~~~~";
        let err = call_parse(line, &header).expect_err("expected InvalidPos error");
        assert!(
            matches!(err, SamError::MalformedRecord { source: SamRecordError::InvalidPos { .. } }),
            "expected InvalidPos, got: {err}"
        );
    }

    // r[verify sam.edge.rname_star]
    /// A record on a known contig but with POS=0 (the SAM marker for "no
    /// position", i.e. unmapped/unplaced) must be skipped, not error the query.
    /// Mirrors the BAM reader skipping pos = -1.
    #[test]
    fn pos_zero_is_skipped_not_errored() {
        let header = make_header();
        let line = b"r\t4\tchr1\t0\t0\t*\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII";
        let result = call_parse(line, &header).expect("POS=0 must not error");
        assert_eq!(result, None, "POS=0 record should be skipped");
    }

    /// An unplaced read (RNAME=*) is skipped before POS is even parsed.
    #[test]
    fn rname_star_is_skipped() {
        let header = make_header();
        let line = b"r\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII";
        let result = call_parse(line, &header).expect("RNAME=* must not error");
        assert_eq!(result, None, "RNAME=* record should be skipped");
    }

    // r[verify types.mapping_quality.field_type]
    #[test]
    fn mapq_255_is_unavailable() {
        let header = make_header();
        let (mut store, mut cigar_buf, mut bases_buf, mut qual_buf, mut aux_buf) =
            make_store_and_bufs();
        let line = b"r\t0\tchr1\t100\t255\t4M\t*\t0\t0\tACGT\t~~~~";
        parse_sam_line(
            line,
            &header,
            tid(0),
            0,
            10000,
            &mut store,
            &mut cigar_buf,
            &mut bases_buf,
            &mut qual_buf,
            &mut aux_buf,
            &mut (),
        )
        .expect("MAPQ 255 is valid SAM");
        let rec = store.records().next().expect("record kept");
        assert_eq!(rec.mapq, MappingQuality::UNAVAILABLE);
    }

    #[test]
    fn error_invalid_mapq() {
        let header = make_header();
        // 999 overflows u8 (max 255)
        let line = b"r\t0\tchr1\t100\t999\t10M\t*\t0\t0\tACGT\t~~~~";
        let err = call_parse(line, &header).expect_err("expected InvalidMapq error");
        assert!(
            matches!(err, SamError::MalformedRecord { source: SamRecordError::InvalidMapq { .. } }),
            "expected InvalidMapq, got: {err}"
        );
    }

    #[test]
    fn error_invalid_cigar_length() {
        let header = make_header();
        // Length 4294967296 exceeds u32::MAX, causing parse_u32 to return None
        let line = b"r\t0\tchr1\t100\t60\t4294967296M\t*\t0\t0\tACGT\t~~~~";
        let err = call_parse(line, &header).expect_err("expected InvalidCigarLength error");
        assert!(
            matches!(
                err,
                SamError::MalformedRecord { source: SamRecordError::InvalidCigarLength { .. } }
            ),
            "expected InvalidCigarLength, got: {err}"
        );
    }

    #[test]
    fn error_unknown_cigar_op() {
        let header = make_header();
        // "10Z" — Z is not a valid CIGAR operation
        let line = b"r\t0\tchr1\t100\t60\t10Z\t*\t0\t0\tACGT\t~~~~";
        let err = call_parse(line, &header).expect_err("expected UnknownCigarOp error");
        assert!(
            matches!(
                err,
                SamError::MalformedRecord { source: SamRecordError::UnknownCigarOp { op: b'Z' } }
            ),
            "expected UnknownCigarOp {{ op: 'Z' }}, got: {err}"
        );
    }

    // TODO: HeaderNotUtf8 requires crafting a BGZF stream with invalid UTF-8 header bytes —
    // it is raised by read_sam_header() which reads from a BgzfReader, not from parse_sam_line.

    // r[verify sam.reader.reject_uncompressed]
    #[test]
    fn reject_uncompressed_sam() {
        let dir = tempfile::tempdir().expect("failed to create temp dir");
        let sam_path = dir.path().join("test.sam");
        {
            let mut f = File::create(&sam_path).expect("failed to create temp file");
            f.write_all(b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n").expect("failed to write");
        }

        let err = IndexedSamReader::open(&sam_path).expect_err("should reject uncompressed SAM");
        assert!(
            matches!(err, SamError::UncompressedSam { .. }),
            "expected UncompressedSam, got: {err}"
        );
        let msg = err.to_string();
        assert!(msg.contains("bgzip"), "error message should mention bgzip: {msg}");
        assert!(msg.contains("test.sam"), "error message should contain file name: {msg}");
    }

    fn call_parse_with_region(
        line: &[u8],
        header: &BamHeader,
        start: i64,
        end: i64,
    ) -> Result<(Option<bool>, RecordStore), SamError> {
        let (mut store, mut cigar_buf, mut bases_buf, mut qual_buf, mut aux_buf) =
            make_store_and_bufs();
        let result = parse_sam_line(
            line,
            header,
            tid(0),
            start,
            end,
            &mut store,
            &mut cigar_buf,
            &mut bases_buf,
            &mut qual_buf,
            &mut aux_buf,
            &mut (),
        )?;
        Ok((result.map(|record| record.kept), store))
    }

    // r[verify sam.reader.overlap_inclusive]
    /// A record starting exactly at the query end still overlaps: both ends of
    /// the query are inclusive.
    #[test]
    fn overlap_filter_record_starting_at_query_end() {
        let header = make_header();
        // pos=100 (1-based) -> 0-based 99; 4M covers 99..=102.
        let line = b"r\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\t~~~~";
        let (result, store) =
            call_parse_with_region(line, &header, 0, 99).expect("parse should succeed");
        assert!(result.is_some(), "record starting at the inclusive query end overlaps");
        assert_eq!(store.len(), 1);
    }

    // r[verify sam.reader.overlap_inclusive]
    /// A record whose *last* covered base is the query start overlaps by one
    /// base. The old half-open test dropped exactly this record.
    #[test]
    fn overlap_filter_record_ending_at_query_start() {
        let header = make_header();
        // 4M at 0-based 99 covers 99..=102, so end_pos == 102.
        let line = b"r\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\t~~~~";
        let (result, store) =
            call_parse_with_region(line, &header, 102, 200).expect("parse should succeed");
        assert!(result.is_some(), "record covering the query start overlaps");
        assert_eq!(store.len(), 1);
    }

    // r[verify sam.reader.overlap_inclusive]
    #[test]
    fn overlap_filter_rejects_records_outside_the_region() {
        let header = make_header();
        let line = b"r\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\t~~~~";
        // Region ends one before the record starts.
        let (before, store) =
            call_parse_with_region(line, &header, 0, 98).expect("parse should succeed");
        assert!(before.is_none(), "record starting past the region must be filtered");
        assert_eq!(store.len(), 0);
        // Region starts one past the record's last covered base.
        let (after, store) =
            call_parse_with_region(line, &header, 103, 200).expect("parse should succeed");
        assert!(after.is_none(), "record ending before the region must be filtered");
        assert_eq!(store.len(), 0);
    }

    // r[verify sam.reader.overlap_inclusive]
    #[test]
    fn overlap_filter_overlapping() {
        let header = make_header();
        let line = b"r\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\t~~~~";
        let (result, store) =
            call_parse_with_region(line, &header, 100, 103).expect("parse should succeed");
        assert!(result.is_some(), "record should overlap region [100, 103]");
        assert_eq!(store.len(), 1);
    }

    // r[verify sam.record.aux_parse_strict+2]
    #[test]
    fn error_malformed_aux_int_tag() {
        let header = make_header();
        let line = b"r\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\t~~~~\tNH:i:NOTANUMBER";
        let err = call_parse(line, &header).expect_err("expected InvalidAuxValue error");
        assert!(
            matches!(
                err,
                SamError::MalformedRecord { source: SamRecordError::InvalidAuxValue { .. } }
            ),
            "expected InvalidAuxValue, got: {err}"
        );
    }

    // r[verify sam.record.parse_int_nonempty]
    #[test]
    fn parse_u8_rejects_empty() {
        assert_eq!(parse_u8(b""), None, "parse_u8 should reject empty input");
    }

    // r[verify sam.record.parse_int_nonempty]
    #[test]
    fn parse_u16_rejects_empty() {
        assert_eq!(parse_u16(b""), None, "parse_u16 should reject empty input");
    }

    // r[verify sam.record.aux_int_range+3]
    #[test]
    fn aux_ints_outside_every_bam_type_are_rejected() {
        for field in ["NH:i:4294967296", "NH:i:-2147483649", "NH:i:9223372036854775808"] {
            let err = aux_bytes(field).unwrap_err();
            assert!(
                matches!(
                    err,
                    SamError::MalformedRecord { source: SamRecordError::AuxIntOutOfRange { .. } }
                ),
                "{field}: {err:?}"
            );
        }
        assert_eq!(aux_bytes("NH:i:4294967295").unwrap(), b"NHI\xff\xff\xff\xff");
        assert_eq!(aux_bytes("NH:i:-2147483648").unwrap(), b"NHi\x00\x00\x00\x80");
    }

    // r[verify sam.record.aux_tags+2]
    /// htslib's width choice: unsigned first, signed only for a leading `-`.
    #[test]
    fn aux_ints_are_stored_unsigned_first() {
        assert_eq!(aux_bytes("NH:i:1").unwrap(), b"NHC\x01");
        assert_eq!(aux_bytes("NH:i:-1").unwrap(), b"NHc\xff");
        assert_eq!(aux_bytes("NH:i:40000").unwrap(), b"NHS\x40\x9c");
    }

    // r[verify sam.record.aux_parse_strict+2]
    /// Fields the old parser skipped silently now fail the record.
    #[test]
    fn malformed_aux_fields_fail_the_record() {
        let header = make_header();
        for aux in ["NM", "NM:i", "NMxi:1", "NM:q:1", "NM:A:"] {
            let line = format!("r\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\t~~~~\t{aux}");
            let err = call_parse(line.as_bytes(), &header).expect_err(aux);
            assert!(matches!(err, SamError::MalformedRecord { .. }), "{aux}: {err:?}");
        }
    }

    // r[verify sam.record.aux_int_range+3]
    #[test]
    fn aux_tag_with_huge_int_returns_error() {
        let header = make_header();
        // 9999999999999 exceeds u32::MAX
        let line = b"r\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\t~~~~\tNH:i:9999999999999";
        let err = call_parse(line, &header).expect_err("expected AuxIntOutOfRange error");
        assert!(
            matches!(
                err,
                SamError::MalformedRecord { source: SamRecordError::AuxIntOutOfRange { .. } }
            ),
            "expected AuxIntOutOfRange, got: {err}"
        );
    }

    // r[depends sam.index.csi+2]
    #[test]
    fn csi_index_returns_clear_error() {
        let dir = tempfile::tempdir().expect("failed to create temp dir");
        let sam_path = dir.path().join("test.sam.gz");

        // Create a fake .csi file next to the SAM
        let csi_path = dir.path().join("test.sam.gz.csi");
        File::create(&csi_path).expect("failed to create .csi file");

        let err = find_tabix_index(&sam_path, &make_header()).expect_err("should reject CSI index");
        assert!(
            matches!(err, SamError::CsiNotSupported { .. }),
            "expected CsiNotSupported, got: {err}"
        );
        let msg = err.to_string();
        assert!(msg.contains("CSI"), "error message should mention CSI: {msg}");
        assert!(msg.contains("tabix"), "error message should mention tabix: {msg}");
    }
}
