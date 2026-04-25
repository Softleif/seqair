//! Open and query CRAM files. [`IndexedCramReader`] handles CRAM v3.0/v3.1 including
//! multi-ref slices, embedded references, coordinate clamping, and per-slice MD5 verification.
//! Call [`IndexedCramReader::fork`] for cheap per-thread readers.

use super::{
    block,
    compression_header::CompressionHeader,
    container::ContainerHeader,
    index::{self, CramIndex, CramIndexError},
    slice,
};
use crate::bam::record::DecodeError;
use crate::bam::record_store::RecordStore;
use crate::bam::{BamHeader, BamHeaderError, BgzfError};
use crate::fasta::{FastaError, IndexedFastaReader};
use seqair_types::{Base, Pos0, Pos1, SmolStr};
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use tracing::instrument;

#[cfg(feature = "fuzz")]
use std::io::Cursor;

#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum CramError {
    #[error("I/O error opening {path}")]
    Open { path: PathBuf, source: std::io::Error },

    #[error("invalid CRAM magic (expected 'CRAM', found {found:?})")]
    InvalidMagic { found: [u8; 4] },

    #[error("unsupported CRAM version {major}.{minor}")]
    UnsupportedVersion { major: u8, minor: u8 },

    #[error("CRAM header error")]
    Header {
        #[from]
        source: BamHeaderError,
    },

    #[error("CRAI index not found for {cram_path} (expected .cram.crai or .crai)")]
    IndexNotFound { cram_path: PathBuf },

    #[error(transparent)]
    IndexParse {
        #[from]
        source: CramIndexError,
    },

    #[error("CRC32 mismatch in {context}: expected {expected:#010x}, got {found:#010x}")]
    ChecksumMismatch { context: &'static str, expected: u32, found: u32 },

    #[error(
        "unsupported compression codec {method} in block (content_type={content_type}, \
             content_id={content_id}). Convert to BAM with `samtools view -b`."
    )]
    UnsupportedCodec { method: u8, content_type: u8, content_id: i32 },

    #[error("unsupported encoding type {encoding_id}")]
    UnsupportedEncoding { encoding_id: i32 },

    #[error("truncated CRAM data: {context}")]
    Truncated { context: &'static str },

    #[error("FASTA reference error")]
    Fasta {
        #[from]
        source: FastaError,
    },

    #[error("reference MD5 mismatch for {contig}:{start}-{end}")]
    ReferenceMd5Mismatch { contig: SmolStr, start: u64, end: u64 },

    #[error(
        "missing reference for {contig} and no embedded reference in slice. \
             Consider setting REF_PATH/REF_CACHE for automatic reference download."
    )]
    MissingReference { contig: SmolStr },

    #[error("multi-ref slices not yet supported (Phase 2)")]
    MultiRefNotSupported,

    #[error("external block not found for content_id {content_id}")]
    ExternalBlockNotFound { content_id: i32 },

    // ── Block decoding ───────────────────────────────────────────────────────
    #[error("unknown block content type: {content_type}")]
    UnknownBlockContentType { content_type: u8 },

    #[error("gzip decompression failed in block")]
    GzipDecompressionFailed { source: libdeflater::DecompressionError },

    #[error("bzip2 decompression failed in block")]
    Bzip2DecompressionFailed { source: std::io::Error },

    #[error("lzma decompression failed in block")]
    LzmaDecompressionFailed { source: std::io::Error },

    // ── Block type mismatch ──────────────────────────────────────────────────
    #[error("expected compression header block, got {found:?}")]
    ExpectedCompressionHeader { found: block::ContentType },

    #[error("expected file header block, got {found:?}")]
    ExpectedFileHeader { found: block::ContentType },

    #[error("expected slice header block, got {found:?}")]
    ExpectedSliceHeader { found: block::ContentType },

    // ── Reader-level ─────────────────────────────────────────────────────────
    #[error("unknown target ID {tid}")]
    UnknownTid { tid: u32 },

    #[error("invalid CRAM alignment position {value}: out of valid range")]
    InvalidPosition { value: i64 },

    #[error("file header text is not valid UTF-8")]
    HeaderNotUtf8 { source: std::str::Utf8Error },

    // ── Slice-level ──────────────────────────────────────────────────────────
    #[error("no core data block in slice")]
    MissingCoreDataBlock,

    #[error("unknown CRAM feature code {feature_code:#04x} ('{}')", *feature_code as char)]
    UnknownFeatureCode { feature_code: u8 },

    // ── Encoding-level ───────────────────────────────────────────────────────
    #[error("Huffman alphabet size ({alphabet_size}) != bit_lengths size ({bit_lengths_size})")]
    HuffmanSizeMismatch { alphabet_size: usize, bit_lengths_size: usize },

    #[error("external byte array encoding for content_id={content_id} requires explicit length")]
    ExternalByteArrayNeedsLength { content_id: i32 },

    // ── rANS ─────────────────────────────────────────────────────────────────
    #[error("invalid rANS 4x8 order: {order}")]
    InvalidRansOrder { order: u8 },

    #[error("rANS Nx16 stripe chunk_count must be > 0")]
    RansStripeZeroChunks,

    #[error("rANS Nx16 bit-pack symbol_count must be > 0")]
    RansBitPackZeroSymbols,

    #[error("rANS Nx16 bit-pack symbol count must be <= 16, got {symbol_count}")]
    RansBitPackTooManySymbols { symbol_count: usize },

    // ── tok3 ─────────────────────────────────────────────────────────────────
    #[error("invalid tok3 token type: {token_type}")]
    InvalidTok3TokenType { token_type: u8 },

    #[error("tok3 adaptive arithmetic coder not supported")]
    Tok3ArithmeticCoderUnsupported,

    #[error("tok3 dup position {dup_pos} out of range")]
    Tok3DupPositionOutOfRange { dup_pos: usize },

    #[error("tok3 Delta token requires Digits predecessor, got token discriminant {found}")]
    Tok3DeltaRequiresDigits { found: u8 },

    #[error("tok3 Delta0 token requires PaddedDigits predecessor, got token discriminant {found}")]
    Tok3Delta0RequiresPaddedDigits { found: u8 },

    #[error("tok3 distance {distance} exceeds name index {name_index}")]
    Tok3DistanceExceedsIndex { distance: usize, name_index: usize },

    #[error("tok3 dup reference index {index} out of range")]
    Tok3DupRefOutOfRange { index: usize },

    #[error("tok3 name_count {count} exceeds limit {limit}")]
    Tok3NameCountExceedsLimit { count: usize, limit: usize },

    // ── Codec safety ─────────────────────────────────────────────────────────
    #[error("uint7 overflow: more than 5 continuation bytes")]
    Uint7Overflow,

    #[error("frequency normalization overflow: sum {sum} cannot be represented in u32")]
    FrequencyNormalizationOverflow { sum: u64 },

    #[error("invalid length from ITF8: {value} is negative")]
    InvalidLength { value: i32 },

    #[error("invalid BAM flags value: {value:#06x}")]
    InvalidBamFlags { value: i32, source: std::num::TryFromIntError },

    #[error("BGZF error")]
    Bgzf {
        #[from]
        source: BgzfError,
    },

    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("record decode error")]
    RecordDecode {
        #[from]
        source: DecodeError,
    },

    #[error(
        "claimed size {size} exceeds allocation limit ({limit}) in {context} — data may be corrupt"
    )]
    AllocationTooLarge { size: usize, limit: usize, context: &'static str },
}

/// Maximum allocation size (256 MiB) for a single decompressed block or parsed structure.
/// Real CRAM data stays well under this; exceeding it indicates corrupt/malicious input.
pub const MAX_ALLOC_SIZE: usize = 256 * 1024 * 1024;

/// Check that a size is within the allocation limit.
pub(crate) fn check_alloc_size(size: usize, context: &'static str) -> Result<(), CramError> {
    if size > MAX_ALLOC_SIZE {
        return Err(CramError::AllocationTooLarge { size, limit: MAX_ALLOC_SIZE, context });
    }
    Ok(())
}

pub struct CramShared {
    pub index: CramIndex,
    pub header: BamHeader,
    pub cram_path: PathBuf,
    pub fasta_path: PathBuf,
}

pub struct IndexedCramReader<R: Read + Seek = File> {
    file: R,
    fasta: IndexedFastaReader<R>,
    shared: Arc<CramShared>,
    // Scratch buffers reused across fetch_into calls
    container_buf: Vec<u8>,
    cigar_buf: Vec<u8>,
    bases_buf: Vec<Base>,
    qual_buf: Vec<u8>,
    aux_buf: Vec<u8>,
    ref_seq_buf: Vec<u8>,
}

impl<R: Read + Seek> std::fmt::Debug for IndexedCramReader<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("IndexedCramReader").field("path", &self.shared.cram_path).finish()
    }
}

impl IndexedCramReader<File> {
    /// Open a CRAM file with a FASTA reference for sequence reconstruction.
    #[instrument(level = "debug", fields(cram = %cram_path.display(), fasta = %fasta_path.display()))]
    pub fn open(cram_path: &Path, fasta_path: &Path) -> Result<Self, CramError> {
        let mut file = File::open(cram_path)
            .map_err(|source| CramError::Open { path: cram_path.to_path_buf(), source })?;

        // Read and verify file definition (26 bytes)
        let mut file_def = [0u8; 26];
        file.read_exact(&mut file_def)?;

        // r[impl cram.file.magic]
        if file_def[..4] != *b"CRAM" {
            return Err(CramError::InvalidMagic {
                found: [file_def[0], file_def[1], file_def[2], file_def[3]],
            });
        }

        // r[impl cram.scope.versions]
        let major = file_def[4];
        let minor = file_def[5];
        if major != 3 {
            return Err(CramError::UnsupportedVersion { major, minor });
        }

        // Read header container
        let header = read_header_container(&mut file)?;
        // r[impl unified.sort_order]
        header.validate_sort_order()?;

        // Find and parse CRAI index
        let crai_path = index::find_crai_path(cram_path)?;
        let cram_index = CramIndex::from_path(&crai_path)?;

        // r[impl cram.scope.reference_required]
        // Open FASTA reader
        let fasta = IndexedFastaReader::open(fasta_path)?;

        Ok(IndexedCramReader {
            file,
            fasta,
            shared: Arc::new(CramShared {
                index: cram_index,
                header,
                cram_path: cram_path.to_path_buf(),
                fasta_path: fasta_path.to_path_buf(),
            }),
            container_buf: Vec::new(),
            cigar_buf: Vec::new(),
            bases_buf: Vec::new(),
            qual_buf: Vec::new(),
            aux_buf: Vec::new(),
            ref_seq_buf: Vec::new(),
        })
    }

    pub fn fork(&self) -> Result<Self, CramError> {
        let file = File::open(&self.shared.cram_path)
            .map_err(|source| CramError::Open { path: self.shared.cram_path.clone(), source })?;
        let fasta = self.fasta.fork()?;
        Ok(IndexedCramReader {
            file,
            fasta,
            shared: Arc::clone(&self.shared),
            container_buf: Vec::new(),
            cigar_buf: Vec::new(),
            bases_buf: Vec::new(),
            qual_buf: Vec::new(),
            aux_buf: Vec::new(),
            ref_seq_buf: Vec::new(),
        })
    }
}

#[cfg(feature = "fuzz")]
impl IndexedCramReader<Cursor<Vec<u8>>> {
    pub fn from_bytes(
        cram_data: Vec<u8>,
        crai_data: &[u8],
        fasta_gz_data: Vec<u8>,
        fai_contents: &str,
        gzi_data: &[u8],
    ) -> Result<Self, CramError> {
        // Parse the 26-byte file definition
        let mut file_def = [0u8; 26];
        file_def.copy_from_slice(
            cram_data.get(..26).ok_or(CramError::Truncated { context: "CRAM file definition" })?,
        );

        // r[impl cram.file.magic]
        if file_def[..4] != *b"CRAM" {
            return Err(CramError::InvalidMagic {
                found: [file_def[0], file_def[1], file_def[2], file_def[3]],
            });
        }

        // r[impl cram.scope.versions]
        let major = file_def[4];
        let minor = file_def[5];
        if major != 3 {
            return Err(CramError::UnsupportedVersion { major, minor });
        }

        let mut cursor = Cursor::new(cram_data.clone());
        cursor.seek(SeekFrom::Start(26))?;
        let header = read_header_container(&mut cursor)?;
        header.validate_sort_order()?;

        let crai_text =
            std::str::from_utf8(crai_data).map_err(|source| CramError::HeaderNotUtf8 { source })?;
        let cram_index = CramIndex::from_text(crai_text)?;

        let fasta = IndexedFastaReader::from_bytes(fasta_gz_data, fai_contents, gzi_data)?;

        Ok(IndexedCramReader {
            file: Cursor::new(cram_data),
            fasta,
            shared: Arc::new(CramShared {
                index: cram_index,
                header,
                cram_path: PathBuf::from("<fuzz>"),
                fasta_path: PathBuf::from("<fuzz>"),
            }),
            container_buf: Vec::new(),
            cigar_buf: Vec::new(),
            bases_buf: Vec::new(),
            qual_buf: Vec::new(),
            aux_buf: Vec::new(),
            ref_seq_buf: Vec::new(),
        })
    }
}

impl<R: Read + Seek> IndexedCramReader<R> {
    pub fn header(&self) -> &BamHeader {
        &self.shared.header
    }

    // r[impl region_buf.not_cram]
    // r[impl cram.container.region_skip]
    // r[impl cram.perf.slice_granularity]
    // r[impl cram.perf.reference_caching]
    // r[impl cram.perf.codec_overhead]
    /// Fetch records overlapping `[start, end)` (0-based) for the given tid.
    pub fn fetch_into(
        &mut self,
        tid: u32,
        start: Pos0,
        end: Pos0,
        store: &mut RecordStore,
    ) -> Result<usize, CramError> {
        self.fetch_into_customized(tid, start, end, store, &mut ()).map(|c| c.kept)
    }

    // r[impl unified.fetch_into_customized]
    // r[impl cram.fetch_into_filtered.push_time]
    /// Customized variant: each record that passes the reader's built-in
    /// overlap/tid/unmapped checks is pushed into the store, with
    /// `customize.keep_record` consulted at push time. Rejection triggers
    /// the same zero-waste rollback used by BAM/SAM — the slabs are
    /// truncated back to their pre-push lengths. The returned
    /// [`FetchCounts`](crate::reader::FetchCounts) reports `fetched`
    /// (produced by the reader) vs `kept` (survived the filter).
    pub fn fetch_into_customized<E: super::super::bam::record_store::CustomizeRecordStore>(
        &mut self,
        tid: u32,
        start: Pos0,
        end: Pos0,
        store: &mut RecordStore,
        customize: &mut E,
    ) -> Result<crate::reader::FetchCounts, CramError> {
        store.clear();

        let start_u64 = start.as_u64();
        let end_u64 = end.as_u64();

        #[expect(
            clippy::cast_possible_wrap,
            reason = "tid bounded by BAM header limits (MAX_REFERENCES = 1M), fits i32"
        )]
        let entries = self.shared.index.query(tid as i32, start_u64, end_u64);
        if entries.is_empty() {
            return Ok(crate::reader::FetchCounts::default());
        }

        let mut fetched_total = 0usize;
        let mut kept_total = 0usize;

        // Group entries by container_offset (multiple slices may be in same container)
        let mut container_offsets: Vec<u64> = entries.iter().map(|e| e.container_offset).collect();
        container_offsets.sort_unstable();
        container_offsets.dedup();

        // Get reference name for FASTA lookup
        let ref_name =
            self.shared.header.target_name(tid).ok_or(CramError::UnknownTid { tid })?.to_string();

        for &container_offset in &container_offsets {
            self.file.seek(SeekFrom::Start(container_offset))?;

            // Read container header
            let mut header_buf = [0u8; 1024]; // container headers are typically < 100 bytes
            let bytes_read = self.file.read(&mut header_buf)?;
            let container_header = ContainerHeader::parse(
                header_buf
                    .get(..bytes_read)
                    .ok_or(CramError::Truncated { context: "container header buf" })?,
            )?;

            if container_header.is_eof() {
                continue;
            }

            // Read container data; non-negative checked, so try_from is infallible here.
            let data_len = usize::try_from(container_header.length)
                .map_err(|_| CramError::Truncated { context: "container negative length" })?;
            check_alloc_size(data_len, "container data")?;
            self.container_buf.clear();
            self.container_buf.resize(data_len, 0);
            self.file.seek(SeekFrom::Start(
                container_offset
                    .checked_add(container_header.header_size as u64)
                    .ok_or(CramError::Truncated { context: "container seek offset overflow" })?,
            ))?;
            self.file.read_exact(&mut self.container_buf)?;

            // Parse compression header (first block)
            let (comp_block, _) = block::parse_block(&self.container_buf)?;
            if comp_block.content_type != block::ContentType::CompressionHeader {
                return Err(CramError::ExpectedCompressionHeader {
                    found: comp_block.content_type,
                });
            }
            let ch = CompressionHeader::parse(&comp_block.data)?;

            // Fetch reference sequence for the query tid's range.
            // For multi-ref containers (ref_seq_id=-2), the container's alignment
            // range covers multiple references — use the CRAI entry's range for our tid instead.
            let (ref_start, ref_end_clamped) = if container_header.ref_seq_id == -2 {
                // Multi-ref: find the CRAI entry for this container + our tid
                let crai_entry = entries.iter().find(|e| e.container_offset == container_offset);
                match crai_entry {
                    Some(e) if e.alignment_span > 0 => {
                        let s = Pos1::try_from(e.alignment_start.max(1))
                            .map_err(|_| CramError::InvalidPosition { value: e.alignment_start })?
                            .to_zero_based()
                            .as_u64();
                        let e_end =
                            s.checked_add(e.alignment_span as u64).ok_or(CramError::Truncated {
                                context: "crai alignment end overflow",
                            })?;
                        let ref_len = self.shared.header.target_len(tid).unwrap_or(0);
                        (s, e_end.min(ref_len))
                    }
                    _ => {
                        // No CRAI span info — fetch the full contig
                        let ref_len = self.shared.header.target_len(tid).unwrap_or(0);
                        (start_u64.min(ref_len), end_u64.min(ref_len))
                    }
                }
            } else {
                let ref_start = Pos1::try_from(container_header.alignment_start.max(1))
                    .map_err(|_| CramError::InvalidPosition {
                        value: i64::from(container_header.alignment_start),
                    })?
                    .to_zero_based()
                    .as_u64();
                let ref_end = ref_start
                    .checked_add(container_header.alignment_span as u64)
                    .ok_or(CramError::Truncated { context: "container alignment end overflow" })?;
                let ref_len = self.shared.header.target_len(tid).unwrap_or(0);
                (ref_start, ref_end.min(ref_len))
            };

            self.ref_seq_buf.clear();
            if ref_start < ref_end_clamped {
                // r[impl cram.edge.missing_reference]
                self.fasta
                    .fetch_seq_into(
                        &ref_name,
                        Pos0::try_from(ref_start)
                            .map_err(|_| CramError::InvalidPosition {
                            #[expect(
                                clippy::cast_possible_wrap,
                                reason = "error reporting only; value may exceed i64::MAX but is only used for diagnostics"
                            )]
                            value: ref_start as i64,
                        })?,
                        Pos0::try_from(ref_end_clamped).unwrap_or(Pos0::max_value()),
                        &mut self.ref_seq_buf,
                    )
                    .map_err(|e| match &e {
                        FastaError::SequenceNotFound { .. } => {
                            CramError::MissingReference { contig: SmolStr::new(&ref_name) }
                        }
                        _ => CramError::from(e),
                    })?;
            }

            // Decode each slice that belongs to this container and overlaps our query
            for &landmark in &container_header.landmarks {
                let slice_offset = usize::try_from(landmark)
                    .map_err(|_| CramError::InvalidLength { value: landmark })?;

                // r[impl cram.edge.coordinate_clamp]
                let (slice_fetched, slice_kept) = slice::decode_slice(
                    &ch,
                    &self.container_buf,
                    slice_offset,
                    &self.ref_seq_buf,
                    ref_start.cast_signed(),
                    &self.shared.header,
                    tid,
                    start,
                    end,
                    store,
                    &mut self.cigar_buf,
                    &mut self.bases_buf,
                    &mut self.qual_buf,
                    &mut self.aux_buf,
                    customize,
                )?;
                fetched_total = fetched_total.wrapping_add(slice_fetched);
                kept_total = kept_total.wrapping_add(slice_kept);
            }
        }

        Ok(crate::reader::FetchCounts { fetched: fetched_total, kept: kept_total })
    }
}

// r[impl cram.file.header_container]
/// Read the file header container and extract the SAM header text.
fn read_header_container<R: Read + Seek>(file: &mut R) -> Result<BamHeader, CramError> {
    // Read enough for the container header (usually < 100 bytes)
    let mut buf = [0u8; 4096];
    let pos = file.stream_position()?;
    let bytes_read = file.read(&mut buf)?;

    let container_header = ContainerHeader::parse(
        buf.get(..bytes_read).ok_or(CramError::Truncated { context: "header container buf" })?,
    )?;

    // Read the container data; non-negative checked above so the cast is safe.
    let data_len = usize::try_from(container_header.length)
        .map_err(|_| CramError::Truncated { context: "header container negative length" })?;
    check_alloc_size(data_len, "header container data")?;
    let mut data = vec![0u8; data_len];
    file.seek(SeekFrom::Start(
        pos.checked_add(container_header.header_size as u64)
            .ok_or(CramError::Truncated { context: "header container seek offset overflow" })?,
    ))?;
    file.read_exact(&mut data)?;

    // Parse the file header block
    let (blk, _) = block::parse_block(&data)?;
    if blk.content_type != block::ContentType::FileHeader {
        return Err(CramError::ExpectedFileHeader { found: blk.content_type });
    }

    // Block data: i32 header_text_length + UTF-8 text
    if blk.data.len() < 4 {
        return Err(CramError::Truncated { context: "file header block data" });
    }
    let text_len_i32 = i32::from_le_bytes(
        blk.data
            .get(..4)
            .ok_or(CramError::Truncated { context: "file header block data" })?
            .try_into()
            .map_err(|_| CramError::Truncated { context: "file header block data" })?,
    );
    let text_len = usize::try_from(text_len_i32)
        .map_err(|_| CramError::InvalidLength { value: text_len_i32 })?;
    let text_end = 4usize
        .checked_add(text_len)
        .ok_or(CramError::Truncated { context: "file header text length overflow" })?;
    let text_bytes =
        blk.data.get(4..text_end).ok_or(CramError::Truncated { context: "file header text" })?;
    let header_text =
        std::str::from_utf8(text_bytes).map_err(|source| CramError::HeaderNotUtf8 { source })?;

    Ok(BamHeader::from_sam_text(header_text)?)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cram_path() -> &'static Path {
        Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test_v30.cram"))
    }

    fn fasta_path() -> &'static Path {
        Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.fasta.gz"))
    }

    #[test]
    fn open_cram_reader() {
        let reader = IndexedCramReader::open(cram_path(), fasta_path()).unwrap();
        assert!(reader.header().target_count() > 0);
    }

    #[test]
    fn fork_cram_reader() {
        let reader = IndexedCramReader::open(cram_path(), fasta_path()).unwrap();
        let forked = reader.fork().unwrap();
        assert_eq!(reader.header().target_count(), forked.header().target_count());
    }

    #[test]
    fn fetch_into_returns_records() {
        let mut reader = IndexedCramReader::open(cram_path(), fasta_path()).unwrap();
        let mut store = RecordStore::new();

        // Fetch from first reference
        let tid = 0;
        let count =
            reader.fetch_into(tid, Pos0::new(0).unwrap(), Pos0::max_value(), &mut store).unwrap();
        assert!(count > 0, "should fetch records from tid={tid}");
    }

    // r[verify cram.edge.missing_reference]
    #[test]
    fn missing_reference_gives_helpful_error() {
        let err = CramError::MissingReference { contig: SmolStr::new("chr19") };
        let msg = format!("{err}");
        assert!(msg.contains("chr19"), "error should mention the contig name");
        assert!(
            msg.contains("REF_PATH") || msg.contains("REF_CACHE"),
            "error should mention REF_PATH/REF_CACHE as alternatives, got: {msg}"
        );
    }

    #[test]
    fn unknown_tid_error_variant() {
        // UnknownTid is returned in fetch_into when target_name(tid) returns None.
        let err = CramError::UnknownTid { tid: 9999 };
        assert!(matches!(err, CramError::UnknownTid { tid: 9999 }));
        let msg = format!("{err}");
        assert!(msg.contains("9999"), "error message should contain the tid: {msg}");
    }

    #[test]
    fn header_not_utf8_from_invalid_bytes() {
        // read_header_container calls std::str::from_utf8 on the header text bytes and
        // returns CramError::HeaderNotUtf8 on failure. We test the construction path
        // directly since read_header_container is private.
        // Build the invalid byte at runtime to avoid the invalid_from_utf8 lint.
        let invalid_byte = 0xFFu8;
        let invalid_utf8 = std::slice::from_ref(&invalid_byte);
        let err = std::str::from_utf8(invalid_utf8).unwrap_err();
        let cram_err = CramError::HeaderNotUtf8 { source: err };
        assert!(matches!(cram_err, CramError::HeaderNotUtf8 { .. }));
        let msg = format!("{cram_err}");
        assert!(
            msg.contains("UTF-8") || msg.contains("utf-8") || msg.contains("valid"),
            "message: {msg}"
        );
    }
}
