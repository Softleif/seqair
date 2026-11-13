//! Open and query CRAM files. [`IndexedCramReader`] handles CRAM v3.0/v3.1 including
//! multi-ref slices, embedded references, coordinate clamping, and per-slice MD5 verification.
//! Call [`IndexedCramReader::fork`] for cheap per-thread readers.

use super::{
    super::bam::record_store::CustomizeRecordStore,
    block,
    compression_header::CompressionHeader,
    container::{ContainerHeader, RefSeqId},
    index::{self, CramIndex, CramIndexError},
    slice::{self, DecodedSlice, Keep, SliceReference, SliceScratch},
    slice_cache::{SliceCache, SliceCacheBudget, SliceCacheStats, SliceKey},
    stream::{self, CramStreamError, SliceDecoder, StreamSlice},
};
use crate::bam::RawBamRecord;
use crate::bam::raw_record::Refill;
use crate::bam::record::DecodeError;
use crate::bam::record_store::RecordStore;
use crate::bam::{BamHeader, BamHeaderError, BgzfError};
use crate::fasta::{FastaError, IndexedFastaReader};
use crate::reader::{RecordFields, Tid};
use crate::sam::encode::encode_alignment;
use core::range::RangeInclusive;
use seqair_types::{Pos0, Pos1, SmallVec, SmolStr};
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{Read, Seek};
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

    #[error(
        "CRAI index not found for {cram_path} (expected {cram_path}.crai or the name with .crai for its extension)"
    )]
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

    // r[impl cram.stream.fields.reference_optional]
    /// The reader was opened without a reference, and a slice of `contig`
    /// needs one: its records' sequence was asked for, or a fetch, which
    /// always reads it.
    #[error("a slice of {contig} needs the reference, and the CRAM reader was opened without one")]
    NoReference { contig: SmolStr },

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
    Bzip2DecompressionFailed { source: bz2ngrs::Error },

    #[error("lzma decompression failed in block")]
    LzmaDecompressionFailed { source: std::io::Error },

    #[error(
        "block with compression method {method} decoded to {found} bytes, but its header says {expected}"
    )]
    BlockSizeMismatch { method: u8, expected: usize, found: usize },

    #[error(
        "block with compression method {method} decodes past the {expected} bytes its header says"
    )]
    BlockDecodesPastSize { method: u8, expected: usize },

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

    /// The FASTA's sequence is shorter than the header's `@SQ` `LN` and the
    /// slice needs bases past its end: the FASTA is not the file's reference.
    #[error(
        "reference {contig} in the FASTA is {fasta_len} bases, the CRAM header says {header_len}: wrong reference?"
    )]
    ReferenceTooShort { contig: SmolStr, header_len: u64, fasta_len: u64 },

    /// A record's `RG` names no `@RG` line of the header.
    #[error("CRAM record read group {read_group} is not one of the header's {count}")]
    ReadGroupOutOfRange { read_group: i32, count: usize },

    /// A record whose features rebuild a different number of bases than
    /// its qualities.
    #[error("CRAM record has {bases} bases but {qualities} qualities")]
    QualityLengthMismatch { bases: usize, qualities: usize },

    #[error("unknown CRAM feature code {feature_code:#04x} ('{}')", *feature_code as char)]
    UnknownFeatureCode { feature_code: u8 },

    // ── Encoding-level ───────────────────────────────────────────────────────
    #[error("Huffman alphabet size ({alphabet_size}) != bit_lengths size ({bit_lengths_size})")]
    HuffmanSizeMismatch { alphabet_size: usize, bit_lengths_size: usize },

    #[error("external byte array encoding for content_id={content_id} requires explicit length")]
    ExternalByteArrayNeedsLength { content_id: i32 },

    #[error("BETA encoding bit width {bits} exceeds 32")]
    InvalidBetaBits { bits: u32 },

    // ── rANS ─────────────────────────────────────────────────────────────────
    #[error("invalid rANS 4x8 order: {order}")]
    InvalidRansOrder { order: u8 },

    #[error("invalid rANS Nx16 bits: {bits} (must be ≤ 12)")]
    InvalidRansNx16Bits { bits: u32 },

    #[error(
        "rANS 4x8 header compressed_size ({advertised}) does not match remaining payload \
         length ({actual}) — file truncated or header corrupt"
    )]
    Rans4x8CompressedSizeMismatch { advertised: u32, actual: usize },

    #[error("rANS Nx16 stripe chunk_count must be > 0")]
    RansStripeZeroChunks,

    #[error("rANS Nx16 STRIPE substreams nest deeper than {limit} levels")]
    RansStripeTooDeep { limit: u8 },

    #[error("rANS Nx16 bit-pack symbol_count must be > 0")]
    RansBitPackZeroSymbols,

    #[error("rANS Nx16 bit-pack symbol count must be <= 16, got {symbol_count}")]
    RansBitPackTooManySymbols { symbol_count: usize },

    #[error(
        "malformed rANS alphabet/frequency-table run: starts at sym {start} with len {len} \
         (would extend past 255)"
    )]
    MalformedAlphabetRun { start: u8, len: u8 },

    // ── Arithmetic coder ─────────────────────────────────────────────────────
    #[error("arith: corrupt range-coded data in {context}")]
    ArithCorruptData { context: &'static str },

    #[error("arith stripe stream count must be > 0")]
    ArithStripeZeroStreams,

    #[error("arith stripe streams nested deeper than {limit} levels")]
    ArithStripeTooDeep { limit: u8 },

    #[error("arith decoded {actual} bytes where the stream needs {expected}")]
    ArithLengthMismatch { expected: usize, actual: usize },

    #[error("arith packed length {packed} exceeds the uncompressed length {unpacked}")]
    ArithPackedLengthTooLarge { packed: usize, unpacked: usize },

    #[error("arith EXT data is not bzip2 (unknown magic number)")]
    ArithExtUnknownCodec,

    // ── tok3 ─────────────────────────────────────────────────────────────────
    #[error("invalid tok3 token type: {token_type}")]
    InvalidTok3TokenType { token_type: u8 },

    #[error("tok3 stream copies position {position} type {token_type}, which was never set")]
    Tok3DupStreamUnset { position: usize, token_type: u8 },

    #[error("tok3 block has more than {limit} token positions")]
    Tok3TooManyPositions { limit: usize },

    #[error(
        "tok3 DELTA token requires a DIGITS value in the previous name, got token kind {found} \
         (0 none, 1 text, 2 digits, 3 zero-padded digits, 4 empty)"
    )]
    Tok3DeltaRequiresDigits { found: u8 },

    #[error(
        "tok3 DELTA0 token requires a DIGITS0 value in the previous name, got token kind {found} \
         (0 none, 1 text, 2 digits, 3 zero-padded digits, 4 empty)"
    )]
    Tok3Delta0RequiresPaddedDigits { found: u8 },

    #[error("tok3 MATCH at token position {position} has no value in the previous name to copy")]
    Tok3MatchWithoutValue { position: usize },

    #[error("tok3 names exceed the block's uncompressed length of {limit} bytes")]
    Tok3OutputOverflow { limit: usize },

    #[error("tok3 name_count {count} exceeds the uncompressed length {length}")]
    Tok3NameCountExceedsLength { count: usize, length: usize },

    #[error("tok3 distance {distance} exceeds name index {name_index}")]
    Tok3DistanceExceedsIndex { distance: usize, name_index: usize },

    #[error("tok3 dup reference index {index} out of range")]
    Tok3DupRefOutOfRange { index: usize },

    #[error("tok3 name_count {count} exceeds limit {limit}")]
    Tok3NameCountExceedsLimit { count: usize, limit: usize },

    // ── fqzcomp ──────────────────────────────────────────────────────────────
    #[error("unsupported fqzcomp format version {version} (expected 5)")]
    FqzcompVersion { version: u8 },

    #[error("fqzcomp stream declares zero parameter blocks")]
    FqzcompNoParams,

    #[error("fqzcomp parameter block {block} uses the selector, but max_sel is 0")]
    FqzcompSelectorWithoutRange { block: u8 },

    #[error("malformed fqzcomp {table} table")]
    FqzcompMalformedTable { table: &'static str },

    #[error("fqzcomp selector picks parameter block {index}, but there are only {nparam}")]
    FqzcompParamOutOfRange { index: u16, nparam: usize },

    #[error("fqzcomp record length {len} is zero or exceeds the {remaining} output bytes left")]
    FqzcompRecordLength { len: u32, remaining: usize },

    #[error("fqzcomp duplicate record of length {len} with only {available} bytes decoded")]
    FqzcompDuplicateTooLong { len: u32, available: usize },

    #[error("fqzcomp range-coded data is corrupt or truncated")]
    FqzcompRangeCoder,

    // ── Codec safety ─────────────────────────────────────────────────────────
    #[error("uint7 overflow: more than 5 continuation bytes")]
    Uint7Overflow,

    #[error("frequency normalization overflow: sum {sum} cannot be represented in u32")]
    FrequencyNormalizationOverflow { sum: u64 },

    #[error("invalid length from ITF8: {value} is negative")]
    InvalidLength { value: i32 },

    /// A container or slice header has a reference id below `-2`; only `-1`
    /// (unmapped) and `-2` (multiple references) are negative ids.
    #[error("invalid reference sequence id {value} (expected >= 0, -1 unmapped or -2 multi-ref)")]
    InvalidRefSeqId { value: i32 },

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

    #[error("huffman {context} count {count} reaches the limit of {limit}")]
    HuffmanCountExceedsLimit { count: u32, limit: u32, context: &'static str },
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

/// The largest output one codec call may produce. `MAX_ALLOC_SIZE` in normal
/// builds. Under `cfg(fuzzing)` (cargo-fuzz) it is 4 MiB: a few bytes of a
/// valid rANS, arith or fqzcomp stream can expand to the full limit, and
/// filling 256 MiB per input starves the fuzzer without reaching new code.
/// htscodecs caps at 100 000 bytes in its fuzzing builds; 4 MiB keeps the
/// seeds' real blocks (up to ~2.3 MB) decodable.
#[cfg(not(fuzzing))]
pub(crate) const MAX_CODEC_OUTPUT: usize = MAX_ALLOC_SIZE;
#[cfg(fuzzing)]
pub(crate) const MAX_CODEC_OUTPUT: usize = 4 * 1024 * 1024;

// r[impl io.fuzz.codec_output_cap]
/// Check a codec's claimed output size against [`MAX_CODEC_OUTPUT`].
pub(crate) fn check_codec_output(size: usize, context: &'static str) -> Result<(), CramError> {
    if size > MAX_CODEC_OUTPUT {
        return Err(CramError::AllocationTooLarge { size, limit: MAX_CODEC_OUTPUT, context });
    }
    Ok(())
}

pub struct CramShared {
    pub index: CramIndex,
    pub header: BamHeader,
    /// Pre-parsed `@RG ID:` values, indexed by RG order in the SAM header.
    /// CRAM records reference read groups by index into this list, so the
    /// alternative — re-scanning `header.header_text()` per record inside
    /// `decode_record` — was the dominant cost in `fetch_into_customized`
    /// profiles. Computed once at `open()` time.
    ///
    /// TODO(perf): this lives on `CramShared` because BAM/SAM keep RG in
    /// raw aux bytes and don't need the indexed lookup. If a future
    /// reader (e.g. SAM-from-CRAM converter) starts needing it too,
    /// promote the parsed list to `BamHeader` itself so every format
    /// shares one canonical cache instead of each opening parser
    /// re-scanning the header text.
    pub read_group_ids: Vec<SmolStr>,
    pub cram_path: PathBuf,
    /// `None` for a reader opened without a reference
    /// (`r[cram.stream.fields.reference_optional]`).
    pub fasta_path: Option<PathBuf>,
    /// Decoded slices, shared by the reader and its forks.
    slice_cache: SliceCache,
    /// What every fetch and query of the reader and its forks decodes. Fixed
    /// at open, since the cache's slices hold only these fields.
    pub(crate) fields: RecordFields,
}

/// Parse `@RG ID:` values from SAM header text in declaration order.
pub(crate) fn parse_read_group_ids(header_text: &str) -> Vec<SmolStr> {
    let mut ids = Vec::new();
    for line in header_text.lines() {
        if !line.starts_with("@RG") {
            continue;
        }
        for field in line.split('\t') {
            if let Some(id) = field.strip_prefix("ID:") {
                ids.push(SmolStr::new(id));
                break;
            }
        }
    }
    ids
}

pub struct IndexedCramReader<R: Read + Seek = File> {
    file: R,
    /// `None` when opened without a reference: then only slices that need
    /// none decode.
    fasta: Option<IndexedFastaReader<R>>,
    shared: Arc<CramShared>,
    // Scratch buffers reused across fetch_into calls, so the second record
    // onward pays no allocation.
    container_buf: Vec<u8>,
    ref_seq_buf: Vec<u8>,
    scratch: SliceScratch,
    /// The decoder of [`Self::query_raw`], made by the first one.
    raw: Option<SliceDecoder>,
}

impl<R: Read + Seek> std::fmt::Debug for IndexedCramReader<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("IndexedCramReader").field("path", &self.shared.cram_path).finish()
    }
}

impl IndexedCramReader<File> {
    /// Open a CRAM file with a FASTA reference for sequence reconstruction,
    /// and the `.crai` beside it (`r[index.locate]`).
    #[instrument(level = "debug", fields(cram = %cram_path.display(), fasta = %fasta_path.display()))]
    pub fn open(cram_path: &Path, fasta_path: &Path) -> Result<Self, CramError> {
        Self::open_inner(cram_path, Some(fasta_path), None, RecordFields::all())
    }

    // r[impl index.explicit]
    /// Open a CRAM with the `.crai` at `index`, wherever it is; nothing beside
    /// the CRAM is looked at.
    #[instrument(level = "debug", fields(cram = %cram_path.display(), index = %index.as_ref().display()))]
    pub fn open_with_index(
        cram_path: &Path,
        fasta_path: &Path,
        index: &crate::reader::IndexPath,
    ) -> Result<Self, CramError> {
        Self::open_inner(cram_path, Some(fasta_path), Some(index.as_ref()), RecordFields::all())
    }

    /// Open with the reference at `fasta_path` and the index at `index` if
    /// given, else the one beside the CRAM, decoding only `fields`. Without
    /// a reference, a query decodes only the slices that need none
    /// (`r[cram.stream.fields.reference_optional]`); the caller checks that
    /// `fields` need none.
    pub(crate) fn open_inner(
        cram_path: &Path,
        fasta_path: Option<&Path>,
        index: Option<&Path>,
        fields: RecordFields,
    ) -> Result<Self, CramError> {
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
        let crai_path = match index {
            Some(index) => index.to_path_buf(),
            None => index::find_crai_path(cram_path)?,
        };
        let cram_index = CramIndex::from_path(&crai_path)?;

        // r[impl cram.scope.reference_required]
        // Open FASTA reader
        let fasta = fasta_path.map(IndexedFastaReader::open).transpose()?;

        let read_group_ids = parse_read_group_ids(header.header_text());

        Ok(IndexedCramReader {
            file,
            fasta,
            shared: Arc::new(CramShared {
                index: cram_index,
                header,
                read_group_ids,
                cram_path: cram_path.to_path_buf(),
                fasta_path: fasta_path.map(Path::to_path_buf),
                slice_cache: SliceCache::new(SliceCacheBudget::Auto),
                fields,
            }),
            container_buf: Vec::new(),
            ref_seq_buf: Vec::new(),
            scratch: SliceScratch::default(),
            raw: None,
        })
    }

    // r[impl unified.record_fields.indexed]
    /// Decode against a fork of `fasta` if opened without a reference, as the
    /// pileup's fetches, which read the sequence, need one.
    pub(crate) fn ensure_reference(&mut self, fasta: &IndexedFastaReader) -> Result<(), CramError> {
        if self.fasta.is_none() {
            self.fasta = Some(fasta.fork()?);
        }
        Ok(())
    }

    pub fn fork(&self) -> Result<Self, CramError> {
        let file = File::open(&self.shared.cram_path)
            .map_err(|source| CramError::Open { path: self.shared.cram_path.clone(), source })?;
        let fasta = self.fasta.as_ref().map(IndexedFastaReader::fork).transpose()?;
        self.shared.slice_cache.handle_opened();
        Ok(IndexedCramReader {
            file,
            fasta,
            shared: Arc::clone(&self.shared),
            container_buf: Vec::new(),
            ref_seq_buf: Vec::new(),
            scratch: SliceScratch::default(),
            raw: None,
        })
    }
}

impl<R: Read + Seek> Drop for IndexedCramReader<R> {
    fn drop(&mut self) {
        self.shared.slice_cache.handle_closed();
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
        cursor.seek(std::io::SeekFrom::Start(26))?;
        let header = read_header_container(&mut cursor)?;
        header.validate_sort_order()?;

        let crai_text =
            std::str::from_utf8(crai_data).map_err(|source| CramError::HeaderNotUtf8 { source })?;
        let cram_index = CramIndex::from_text(crai_text)?;

        let fasta = Some(IndexedFastaReader::from_bytes(fasta_gz_data, fai_contents, gzi_data)?);

        let read_group_ids = parse_read_group_ids(header.header_text());

        Ok(IndexedCramReader {
            file: Cursor::new(cram_data),
            fasta,
            shared: Arc::new(CramShared {
                index: cram_index,
                header,
                read_group_ids,
                cram_path: PathBuf::from("<fuzz>"),
                fasta_path: Some(PathBuf::from("<fuzz>")),
                slice_cache: SliceCache::new(SliceCacheBudget::Auto),
                fields: RecordFields::all(),
            }),
            container_buf: Vec::new(),
            ref_seq_buf: Vec::new(),
            scratch: SliceScratch::default(),
            raw: None,
        })
    }
}

impl<R: Read + Seek> IndexedCramReader<R> {
    pub fn header(&self) -> &BamHeader {
        &self.shared.header
    }

    /// What this reader decodes (`r[unified.record_fields]`).
    pub fn fields(&self) -> &RecordFields {
        &self.shared.fields
    }

    // r[impl cram.slice_cache.budget+3]
    /// Set the budget of the decoded-slice cache this reader shares with
    /// its forks (default [`SliceCacheBudget::Auto`]). It bounds the slices
    /// the cache keeps, not the ones fetches in flight hold.
    pub fn set_slice_cache_budget(&self, budget: SliceCacheBudget) {
        self.shared.slice_cache.set_budget(budget);
    }

    /// What the decoded-slice cache this reader shares with its forks has
    /// done so far, for sizing its budget.
    pub fn slice_cache_stats(&self) -> SliceCacheStats {
        self.shared.slice_cache.stats()
    }

    /// Estimate the compressed bytes of the records a query for `span` on
    /// `tid` returns, from the index alone (`r[cram.index.region_bytes]`).
    pub fn estimate_region_bytes(&self, tid: Tid, span: RangeInclusive<Pos0>) -> u64 {
        self.shared.index.region_bytes(tid, span.start.as_u64(), span.last.as_u64())
    }

    // r[impl region_buf.not_cram]
    // r[impl cram.container.region_skip]
    // r[impl cram.perf.slice_granularity]
    // r[impl cram.perf.reference_caching]
    // r[impl cram.perf.codec_overhead]
    /// Fetch records overlapping the closed interval `span` (0-based) on `tid`.
    pub fn fetch_into(
        &mut self,
        tid: Tid,
        span: RangeInclusive<Pos0>,
        store: &mut RecordStore,
    ) -> Result<usize, CramError> {
        self.fetch_into_customized(tid, span, store, &mut ()).map(|c| c.kept)
    }

    // r[impl unified.fetch_into_customized]
    /// Customized variant: each record that passes the reader's built-in
    /// overlap/tid/unmapped checks is pushed into the store, with
    /// `customize.filter` consulted at push time. Rejection triggers
    /// the same zero-waste rollback used by BAM/SAM — the slabs are
    /// truncated back to their pre-push lengths. The returned
    /// [`FetchCounts`](crate::reader::FetchCounts) reports `fetched`
    /// (produced by the reader) vs `kept` (survived the filter).
    pub fn fetch_into_customized<E: CustomizeRecordStore>(
        &mut self,
        tid: Tid,
        span: RangeInclusive<Pos0>,
        store: &mut RecordStore<E::Extra>,
        customize: &mut E,
    ) -> Result<crate::reader::FetchCounts, CramError> {
        store.clear();
        // r[impl interval.empty_span]
        if span.is_empty() {
            return Ok(crate::reader::FetchCounts::default());
        }

        let start_u64 = span.start.as_u64();
        let end_u64 = span.last.as_u64();

        // Its own handle on the shared state, so decoding below can borrow
        // the reader mutably while the index entries and the reference name
        // are borrowed from it.
        let shared = Arc::clone(&self.shared);
        let wanted = wanted_slices(&shared.index, tid, start_u64, end_u64);
        if wanted.is_empty() {
            return Ok(crate::reader::FetchCounts::default());
        }

        let mut counts = crate::reader::FetchCounts::default();

        let ref_name: &str =
            shared.header.target_name(tid).ok_or(CramError::UnknownTid { tid: tid.as_u32() })?;
        shared.slice_cache.record_fetch(wanted.values().map(|slices| slices.len()).sum());

        // r[impl cram.slice_cache.shared]
        // Each wanted slice comes from the cache shared with this reader's
        // forks; only a miss reads its container, parses the compression
        // header and fetches the reference.
        for (&container_offset, slice_offsets) in &wanted {
            // Read and parsed only once a slice needs decoding; `Some(None)`
            // is the EOF container, which holds no slices.
            let mut container: Option<Option<LoadedContainer>> = None;
            for &slice_offset in slice_offsets {
                let key = SliceKey { container_offset, slice_offset, tid };
                let decoded = shared.slice_cache.get_or_decode(key, || {
                    let loaded = match &mut container {
                        Some(loaded) => loaded,
                        None => container.insert(self.load_container(
                            container_offset,
                            tid,
                            ref_name,
                        )?),
                    };
                    match loaded {
                        Some(loaded) => self.decode_slice_in(loaded, slice_offset, tid),
                        None => Ok(DecodedSlice::default()),
                    }
                })?;
                decoded.copy_into(span.start, span.last, store, customize, &mut counts)?;
            }
        }

        Ok(counts)
    }
}

impl<R: Read + Seek> IndexedCramReader<R> {
    // r[impl cram.query_raw]
    /// The records on `tid` overlapping the closed interval `span`
    /// (0-based), in file order, each as htslib reconstructs it: what
    /// `samtools view -b <cram> <region>` writes, MD/NM generated and every
    /// IUPAC code kept. The slices are the ones [`Self::fetch_into`] decodes,
    /// decoded by the stream's decoder.
    ///
    /// ```no_run
    /// use seqair::bam::RawBamRecord;
    /// use seqair::cram::reader::IndexedCramReader;
    /// use seqair_types::Pos0;
    ///
    /// let mut reader = IndexedCramReader::open("in.cram".as_ref(), "ref.fa".as_ref())?;
    /// let tid = reader.header().tid("chr1").expect("chr1 in header");
    /// let span = (Pos0::new(1_000).expect("pos")..=Pos0::new(2_000).expect("pos")).into();
    /// let mut query = reader.query_raw(tid, span);
    /// let mut record = RawBamRecord::default();
    /// while query.read_record(&mut record)? {
    ///     println!("{}", record.qname().escape_ascii());
    /// }
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn query_raw(&mut self, tid: Tid, span: RangeInclusive<Pos0>) -> CramRawQuery<'_, R> {
        let slices: Vec<(u64, u64)> = if span.is_empty() {
            Vec::new()
        } else {
            wanted_slices(&self.shared.index, tid, span.start.as_u64(), span.last.as_u64())
                .into_iter()
                .flat_map(|(container, slices)| slices.into_iter().map(move |s| (container, s)))
                .collect()
        };
        let Self { file, fasta, shared, raw, .. } = self;
        let decoder = raw.get_or_insert_with(|| {
            let prefix = shared.cram_path.file_name().and_then(|n| n.to_str()).unwrap_or("-");
            SliceDecoder::for_query(
                shared.header.clone(),
                shared.read_group_ids.clone(),
                SmolStr::new(prefix),
                shared.fields.clone(),
            )
        });
        CramRawQuery {
            file,
            fasta: fasta.as_mut(),
            decoder,
            slices: slices.into_iter(),
            container: None,
            data: Vec::new(),
            slice: StreamSlice::default(),
            next: 0,
            tid,
            span,
            scratch: RawBamRecord::default(),
            failed: false,
        }
    }
}

// r[impl cram.index.crai_per_slice]
/// The slices the CRAI lists for `start..=end` (0-based) on `tid`: their
/// containers' byte offsets, ascending, each with its slices' offsets,
/// ascending — the file's order, as landmarks ascend.
fn wanted_slices(
    index: &CramIndex,
    tid: Tid,
    start: u64,
    end: u64,
) -> BTreeMap<u64, SmallVec<u64, 4>> {
    let mut wanted: BTreeMap<u64, SmallVec<u64, 4>> = BTreeMap::new();
    for e in &index.query(tid, start, end) {
        let slices = wanted.entry(e.container_offset).or_default();
        if !slices.contains(&e.slice_offset) {
            slices.push(e.slice_offset);
        }
    }
    for slices in wanted.values_mut() {
        slices.sort_unstable();
    }
    wanted
}

/// The records of one CRAM region as [`RawBamRecord`]s; see
/// [`IndexedCramReader::query_raw`]. Borrows the reader until dropped.
pub struct CramRawQuery<'r, R: Read + Seek> {
    file: &'r mut R,
    fasta: Option<&'r mut IndexedFastaReader<R>>,
    decoder: &'r mut SliceDecoder,
    /// `(container offset, slice offset)` of the slices still to decode.
    slices: std::vec::IntoIter<(u64, u64)>,
    /// The container the last slice came from, kept for its next slices.
    container: Option<Box<RawContainer>>,
    /// The bytes of `container`, reused from container to container.
    data: Vec<u8>,
    slice: StreamSlice,
    /// Index in `slice` of the next record to look at.
    next: usize,
    tid: Tid,
    span: RangeInclusive<Pos0>,
    /// Each record is encoded here and swapped out once it validates, so an
    /// error leaves the caller's record as it was.
    scratch: RawBamRecord,
    failed: bool,
}

/// A container a raw query read: its header, bytes and compression header,
/// or nothing for the EOF container.
struct RawContainer {
    offset: u64,
    body: Option<(ContainerHeader, CompressionHeader)>,
}

impl<R: Read + Seek> std::fmt::Debug for CramRawQuery<'_, R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("CramRawQuery")
            .field("slices_left", &self.slices.len())
            .field("failed", &self.failed)
            .finish_non_exhaustive()
    }
}

impl<R: Read + Seek> CramRawQuery<'_, R> {
    pub fn header(&self) -> &BamHeader {
        self.decoder.header()
    }

    // r[impl cram.query_raw.errors]
    /// Read the next record into `record`, reusing its buffer. `Ok(false)` at
    /// the end of the region and any error leave `record` unchanged.
    ///
    /// # Errors
    /// [`CramStreamError::RecordEncode`] or [`CramStreamError::InvalidRecord`]
    /// for a record BAM cannot hold or that fails validation — the query
    /// stays usable. Anything else ends the query: later calls return
    /// [`CramStreamError::Poisoned`].
    pub fn read_record(&mut self, record: &mut RawBamRecord) -> Result<bool, CramStreamError> {
        if self.failed {
            return Err(CramStreamError::Poisoned);
        }
        match self.advance() {
            Ok(true) => {
                std::mem::swap(record, &mut self.scratch);
                Ok(true)
            }
            Ok(false) => Ok(false),
            Err(e) => {
                if !matches!(
                    e,
                    CramStreamError::RecordEncode { .. } | CramStreamError::InvalidRecord { .. }
                ) {
                    self.failed = true;
                }
                Err(e)
            }
        }
    }

    /// Encode the next overlapping record into `scratch`, decoding slices
    /// until there is one. `false` when none is left.
    fn advance(&mut self) -> Result<bool, CramStreamError> {
        loop {
            let found =
                self.slice.records.next_overlapping(self.next, self.span.start, self.span.last);
            if let Some((i, fields)) = found {
                self.next = i.saturating_add(1);
                let record = self.slice.number(i);
                let header = self.decoder.header();
                self.scratch.refill(header, |payload| encode_alignment(payload, &fields)).map_err(
                    |e| match e {
                        Refill::Fill(source) => CramStreamError::RecordEncode { record, source },
                        Refill::Invalid(source) => {
                            CramStreamError::InvalidRecord { record, source }
                        }
                    },
                )?;
                return Ok(true);
            }
            let Some((container_offset, slice_offset)) = self.slices.next() else {
                return Ok(false);
            };
            self.slice = StreamSlice::default();
            self.next = 0;
            self.decode(container_offset, slice_offset)?;
        }
    }

    /// Decode the slice at `slice_offset` of the container at
    /// `container_offset`, reading the container unless the last slice was
    /// in it.
    fn decode(&mut self, container_offset: u64, slice_offset: u64) -> Result<(), CramStreamError> {
        let container = match &mut self.container {
            Some(c) if c.offset == container_offset => c,
            slot => {
                let read =
                    stream::read_container_at(&mut *self.file, container_offset, &mut self.data);
                let body = match read? {
                    Some(header) => {
                        let compression = self.decoder.compression_header(&self.data)?;
                        Some((header, compression))
                    }
                    None => None,
                };
                slot.insert(Box::new(RawContainer { offset: container_offset, body }))
            }
        };
        let Some((header, compression)) = &container.body else { return Ok(()) };
        // r[impl cram.index.crai_per_slice]
        if !header.landmarks.iter().any(|&l| u64::try_from(l).is_ok_and(|l| l == slice_offset)) {
            return Ok(());
        }
        let offset = usize::try_from(slice_offset)
            .map_err(|_| CramError::Truncated { context: "slice offset" })?;
        self.slice = self.decoder.decode_slice(
            self.fasta.as_deref_mut(),
            compression,
            &self.data,
            offset,
            Keep::Tid(self.tid),
        )?;
        Ok(())
    }
}

impl<R: Read + Seek> Iterator for CramRawQuery<'_, R> {
    type Item = Result<RawBamRecord, CramStreamError>;

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

/// A container read and parsed for decoding its slices: its header (for
/// the landmarks), its compression header, and where the reference in the
/// reader's `ref_seq_buf` starts.
struct LoadedContainer {
    header: ContainerHeader,
    compression: CompressionHeader,
    ref_start: u64,
    /// The FASTA lacks the contig; only a slice that needs it fails.
    reference_missing: bool,
}

impl<R: Read + Seek> IndexedCramReader<R> {
    /// Read the container at `container_offset` into `container_buf`, parse
    /// its compression header, and fetch the reference its slices decode
    /// against for `tid` into `ref_seq_buf`. `None` for the EOF container.
    fn load_container(
        &mut self,
        container_offset: u64,
        tid: Tid,
        ref_name: &str,
    ) -> Result<Option<LoadedContainer>, CramError> {
        // r[impl cram.stream.container_header+2]
        let Some(container_header) =
            stream::read_container_at(&mut self.file, container_offset, &mut self.container_buf)
                .map_err(container_read_error)?
        else {
            return Ok(None);
        };

        // Parse compression header (first block)
        let (comp_block, _) =
            block::parse_block_with(&self.container_buf, &mut self.scratch.codec)?;
        if comp_block.content_type != block::ContentType::CompressionHeader {
            return Err(CramError::ExpectedCompressionHeader { found: comp_block.content_type });
        }
        let ch = CompressionHeader::parse(&comp_block.data)?;

        // Fetch reference sequence for the query tid's range.
        // For multi-ref containers (`RefSeqId::Multi`), the container's alignment
        // range covers multiple references — use the CRAI entry's range for our tid instead.
        let (ref_start, ref_end_clamped) = if container_header.ref_seq_id == RefSeqId::Multi {
            // r[impl cram.slice.multi_ref_reference_window+2]
            // Multi-ref: a container holds one CRAI entry per slice per
            // reference, so there is not *a* entry for this container and
            // this tid — there are as many as it has slices touching the
            // tid. The window has to span all of them. Taking only the
            // first stops the reference short of whatever a later slice
            // reaches, and the tail of those reads decodes as `N` with no
            // error at all, because running past the end of the fetched
            // reference is only a warning (r[`cram.slice.ref_bounds_warning`]).
            // They come from the whole index, not the query's entries, so
            // the window — and the decoded slice — is the same whichever
            // query decodes it.
            let ref_len = self.shared.header.target_len(tid).unwrap_or(0);
            let tid_entries = self.shared.index.entries_for(tid);
            let in_container =
                || tid_entries.iter().filter(|e| e.container_offset == container_offset);
            // An entry without an extent says nothing about where its
            // reads lie, so only the whole reference is sure to hold them.
            let span = if in_container().any(|e| e.alignment_span <= 0) {
                Some((0, ref_len))
            } else {
                in_container().try_fold(None::<(u64, u64)>, |window, e| {
                    let s = Pos1::try_from(e.alignment_start.max(1))
                        .map_err(|_| CramError::InvalidPosition { value: e.alignment_start })?
                        .to_zero_based()
                        .as_u64();
                    let e_end = s
                        .checked_add(e.alignment_span.unsigned_abs())
                        .ok_or(CramError::Truncated { context: "crai alignment end overflow" })?;
                    Ok::<_, CramError>(Some(match window {
                        Some((lo, hi)) => (lo.min(s), hi.max(e_end)),
                        None => (s, e_end),
                    }))
                })?
            };
            match span {
                Some((s, e_end)) => (s, e_end.min(ref_len)),
                // No entry for this reference: there is nothing to decode.
                None => (0, 0),
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
        // A contig the FASTA lacks is only an error for a slice that
        // needs it; one that embeds its reference does not.
        let mut reference_missing = false;
        if let Some(fasta) = self.fasta.as_mut()
            && ref_start < ref_end_clamped
        {
            // `ref_end_clamped` is one past the last base the slices
            // reach, and it is above `ref_start`, so the last base is
            // `ref_end_clamped - 1`. Fetching that closed span means a
            // container reaching the last representable position gets its
            // last base — the half-open `end` there is `i32::MAX + 1`,
            // which is not a `Pos0`, and clamping it silently dropped one.
            let position = |value: u64| {
                Pos0::try_from(value).map_err(|_| CramError::InvalidPosition {
                    value: i64::try_from(value).unwrap_or(i64::MAX),
                })
            };
            let ref_span = RangeInclusive {
                start: position(ref_start)?,
                last: position(ref_end_clamped.saturating_sub(1))?,
            };
            // r[impl cram.edge.missing_reference+2]
            match fasta.fetch_seq_into(ref_name, ref_span, &mut self.ref_seq_buf) {
                Ok(()) => {}
                Err(FastaError::SequenceNotFound { .. }) => reference_missing = true,
                Err(e) => return Err(e.into()),
            }
        }
        Ok(Some(LoadedContainer {
            header: container_header,
            compression: ch,
            ref_start,
            reference_missing,
        }))
    }

    /// Decode the slice at `slice_offset` of a loaded container for `tid`.
    /// An offset that is not one of the container's landmarks holds no
    /// slice, and decodes to none.
    fn decode_slice_in(
        &mut self,
        container: &LoadedContainer,
        slice_offset: u64,
        tid: Tid,
    ) -> Result<DecodedSlice, CramError> {
        // CRAI's `slice_offset` matches the container's landmark value
        // (both are byte offsets from the start of the container's data
        // block to the slice header).
        // r[impl cram.index.crai_per_slice]
        let is_landmark = container
            .header
            .landmarks
            .iter()
            .any(|&landmark| u64::try_from(landmark).is_ok_and(|l| l == slice_offset));
        if !is_landmark {
            return Ok(DecodedSlice::default());
        }
        let slice_offset = usize::try_from(slice_offset)
            .map_err(|_| CramError::Truncated { context: "slice offset" })?;
        let bases = (!container.reference_missing).then_some(self.ref_seq_buf.as_slice());
        let reference = match self.fasta {
            Some(_) => SliceReference::Window { bases, start: container.ref_start.cast_signed() },
            // r[impl cram.stream.fields.reference_optional]
            None => SliceReference::None,
        };
        // r[impl unified.record_fields.indexed]
        // The pileup reads the sequence and qualities whatever was asked;
        // only the tags are cut down.
        let fields = self.shared.fields.clone().with_sequence().with_quality();
        // r[impl cram.edge.coordinate_clamp]
        slice::decode_slice(
            &container.compression,
            &self.container_buf,
            slice_offset,
            reference,
            &self.shared.header,
            &self.shared.read_group_ids,
            Keep::Tid(tid),
            &fields,
            // r[impl cram.record.read_name+2]
            slice::Synthesis {
                md_nm: false,
                name_prefix: self.shared.cram_path.file_name().and_then(|n| n.to_str()),
            },
            &mut self.scratch,
        )
    }
}

// r[impl cram.file.header_container]
/// Read the file header container and extract the SAM header text.
// r[impl cram.stream.container_header+2]
fn read_header_container<R: Read>(file: &mut R) -> Result<BamHeader, CramError> {
    let mut data = Vec::new();
    stream::read_container(file, &mut data, || CramStreamError::TruncatedContainer {
        container: 0,
    })
    .map_err(container_read_error)?
    .ok_or(CramError::Truncated { context: "header container" })?;
    header_from_container_data(&data)
}

/// A container read's error as the indexed reader reports it.
fn container_read_error(e: CramStreamError) -> CramError {
    match e {
        CramStreamError::Cram { source } => source,
        CramStreamError::Io { source } => CramError::Io(source),
        _ => CramError::Truncated { context: "container" },
    }
}

/// The SAM header held by the file header container whose body is `data`.
pub(crate) fn header_from_container_data(data: &[u8]) -> Result<BamHeader, CramError> {
    // Parse the file header block. Cold: once per open, before any reader (and
    // its scratch) exists; `parse_block` uses a fresh scratch.
    let (blk, _) = block::parse_block(data)?;
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
    use crate::reader::ResolveTid;
    use std::sync::atomic::Ordering;

    fn cram_path() -> &'static Path {
        Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test_v30.cram"))
    }

    fn fasta_path() -> &'static Path {
        Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.fasta.gz"))
    }

    /// The whole of the first reference of `reader`.
    fn whole_first_reference<R: Read + Seek>(
        reader: &IndexedCramReader<R>,
    ) -> (Tid, RangeInclusive<Pos0>) {
        let tid = crate::bam::test_util::tid(0);
        let len = reader.header().target_len(tid).unwrap();
        let last = Pos0::new(u32::try_from(len.saturating_sub(1)).unwrap()).unwrap();
        (tid, RangeInclusive { start: Pos0::ZERO, last })
    }

    // r[verify cram.stream.fields.dropped_blocks]
    // r[verify unified.record_fields.fetch_into]
    /// A fetch of a reader opened for no tags decompresses no tag block, and
    /// stores no tag but nothing else less.
    #[test]
    fn a_fetch_without_tags_decompresses_no_tag_block() {
        use crate::reader::AuxTags;
        let cram = Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.cram"));
        let tags = super::stream::tests::SeriesBlocks::of(cram).tags;
        let fetch = |fields: RecordFields| {
            let mut reader =
                IndexedCramReader::open_inner(cram, Some(fasta_path()), None, fields).unwrap();
            let (tid, span) = whole_first_reference(&reader);
            let mut store = RecordStore::new();
            reader.fetch_into(tid, span, &mut store).unwrap();
            (store, reader.scratch.codec.decompressed.clone())
        };
        let (full, all) = fetch(RecordFields::all());
        assert!(all.iter().any(|id| tags.contains(id)), "tags are read at all");
        // Sequence and qualities are kept whatever the fields say.
        let (store, decompressed) = fetch(RecordFields::core().with_aux(AuxTags::None));
        assert!(!decompressed.iter().any(|id| tags.contains(id)), "a tag block decompressed");
        assert_eq!(store.len(), full.len());
        for (i, j) in store.indices().zip(full.indices()) {
            let (got, want) = (store.record(i).unwrap(), full.record(j).unwrap());
            assert!(got.aux().is_empty());
            assert_eq!(got.seq(), want.seq());
            assert_eq!(got.qual(), want.qual());
            assert_eq!(got.cigar(), want.cigar());
        }
    }

    // r[verify cram.stream.fields.reference_optional]
    /// Opened without a reference for fields that need none, a raw query
    /// decodes; a fetch, which reads the sequence, fails with its own error.
    #[test]
    fn without_a_reference_queries_decode_and_fetches_fail() {
        let cram = cram_path();
        let mut with = IndexedCramReader::open(cram, fasta_path()).unwrap();
        let mut without =
            IndexedCramReader::open_inner(cram, None, None, RecordFields::core()).unwrap();
        let (tid, span) = whole_first_reference(&with);
        let full: Vec<RawBamRecord> = with.query_raw(tid, span).map(Result::unwrap).collect();
        let core: Vec<RawBamRecord> = without.query_raw(tid, span).map(Result::unwrap).collect();
        assert!(!full.is_empty());
        assert_eq!(core.len(), full.len());
        for (c, f) in core.iter().zip(&full) {
            assert_eq!(c.cigar().collect::<Vec<_>>(), f.cigar().collect::<Vec<_>>());
            assert_eq!((c.pos(), c.flags(), c.seq_len()), (f.pos(), f.flags(), f.seq_len()));
        }
        let err = without.fetch_into(tid, span, &mut RecordStore::new()).unwrap_err();
        assert!(matches!(err, CramError::NoReference { .. }), "{err:?}");
    }

    // r[verify cram.query_raw.errors]
    // r[verify cram.stream.slice_decode]
    /// The number a raw query gives a record (for its errors) is its place
    /// among the file's records, also in multi-reference slices, where the
    /// query keeps only some of a slice's records.
    #[test]
    fn a_raw_query_numbers_records_by_their_place_in_the_file() {
        let data = Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/cram-stream"));
        let (cram, fasta) = (data.join("fidelity_multiref.cram"), data.join("fidelity.fa"));
        let mut stream = super::stream::CramStreamReader::new(File::open(&cram).unwrap()).unwrap();
        let mut decoder = stream.decoder(&fasta).unwrap().name_prefix("fidelity_multiref.cram");
        let mut streamed = Vec::new();
        let mut multi_ref = false;
        while let Some(container) = stream.next_container().unwrap() {
            decoder.decode(&container, &mut streamed).unwrap();
            multi_ref |= container.header.ref_seq_id == RefSeqId::Multi;
        }
        assert!(multi_ref, "the fixture has multi-reference slices");

        let mut reader = IndexedCramReader::open(&cram, &fasta).unwrap();
        let mut checked = 0;
        for name in ["chr1", "chr2"] {
            let tid = name.resolve_tid(reader.header()).unwrap();
            let len = reader.header().target_len(tid).unwrap();
            let span = RangeInclusive {
                start: Pos0::new(0).unwrap(),
                last: Pos0::new(u32::try_from(len - 1).unwrap()).unwrap(),
            };
            let mut query = reader.query_raw(tid, span);
            let mut record = RawBamRecord::default();
            while query.read_record(&mut record).unwrap() {
                let number = query.slice.number(query.next - 1);
                let at = usize::try_from(number).unwrap();
                assert!(streamed[at] == record, "{name}: record {number}");
                checked += 1;
            }
        }
        assert!(checked > 10, "{checked}");
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
        let tid = 0u32.resolve_tid(reader.header()).unwrap();
        let count =
            reader.fetch_into(tid, (Pos0::new(0).unwrap()..=Pos0::MAX).into(), &mut store).unwrap();
        assert!(count > 0, "should fetch records from tid={tid:?}");
    }

    // r[verify cram.slice_cache.shared]
    // r[verify cram.slice_cache.budget+3]
    #[test]
    fn forks_share_decoded_slices() {
        let mut reader = IndexedCramReader::open(cram_path(), fasta_path()).unwrap();
        let mut fork = reader.fork().unwrap();
        let mut store = RecordStore::new();
        let whole = || (Pos0::ZERO..=Pos0::MAX).into();
        let tid0 = 0u32.resolve_tid(reader.header()).unwrap();

        let count = reader.fetch_into(tid0, whole(), &mut store).unwrap();
        let decoded = reader.slice_cache_stats().decoded;
        assert!(decoded > 1, "tid 0 spans several slices");

        assert_eq!(fork.fetch_into(tid0, whole(), &mut store).unwrap(), count);
        assert_eq!(reader.fetch_into(tid0, whole(), &mut store).unwrap(), count);
        assert_eq!(reader.slice_cache_stats().decoded, decoded, "a cached slice was decoded again");

        fork.set_slice_cache_budget(SliceCacheBudget::Bytes(0));
        assert_eq!(reader.slice_cache_stats().budget_bytes, 0, "the budget is shared too");
        assert_eq!(reader.fetch_into(tid0, whole(), &mut store).unwrap(), count);
        assert_eq!(
            reader.slice_cache_stats().decoded,
            decoded * 2,
            "budget 0 still cached a slice"
        );
    }

    // r[verify cram.slice_cache.budget+3]
    #[test]
    fn forks_count_as_handles_until_dropped() {
        let reader = IndexedCramReader::open(cram_path(), fasta_path()).unwrap();
        let handles = || reader.shared.slice_cache.handles.load(Ordering::Relaxed);
        assert_eq!(handles(), 1);
        let forks: Vec<_> = (0..3).map(|_| reader.fork().unwrap()).collect();
        assert_eq!(handles(), 4);
        drop(forks);
        assert_eq!(handles(), 1);
    }

    // r[verify cram.edge.unmapped_reads+2]
    /// A multi-ref slice can interleave a truly unplaced read (RNAME=*, AP=0)
    /// with mapped reads. Decoding that slice for a positional query must not
    /// error on the unplaced read's `alignment_pos = 0` — it has no coordinate
    /// and is mapped to 0-based 0. Mirrors the BAM reader skipping pos = -1.
    #[test]
    fn fetch_into_tolerates_unplaced_read_in_multiref_slice() {
        let cram = Path::new(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../tests/data/unplaced_multiref.cram"
        ));
        let fasta = Path::new(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../tests/data/unplaced_multiref.fa"
        ));
        let mut reader = IndexedCramReader::open(cram, fasta).unwrap();
        let mut store = RecordStore::new();

        // Querying chr1 decodes the multi-ref slice that also holds the unplaced
        // read; before the fix this errored with InvalidPosition { value: 0 }.
        let chr1 = reader.header().tid("chr1").unwrap();
        let count = reader
            .fetch_into(chr1, (Pos0::new(0).unwrap()..=Pos0::MAX).into(), &mut store)
            .expect("query must not error on an unplaced read in a multi-ref slice");
        // Only the chr1 read: the unplaced one is on no reference.
        assert_eq!(count, 1);
        let first = store.indices().next().unwrap();
        assert_eq!(store.record(first).unwrap().qname(), b"r1");
    }

    // r[verify cram.edge.seq_unknown+3]
    /// A mapped read stored with SEQ and QUAL `*` (CRAM flag 0x8) comes back
    /// without bases or qualities, as htslib gives it, and the records after
    /// it decode normally.
    #[test]
    fn mapped_read_without_sequence_has_none() {
        let cram = Path::new(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../tests/data/cram-stream/fidelity_v31.cram"
        ));
        let fasta = Path::new(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../tests/data/cram-stream/fidelity.fa"
        ));
        let mut reader = IndexedCramReader::open(cram, fasta).unwrap();
        let mut store = RecordStore::new();
        let chr1 = reader.header().tid("chr1").unwrap();
        reader.fetch_into(chr1, (Pos0::ZERO..=Pos0::MAX).into(), &mut store).unwrap();
        let record = |name: &[u8]| {
            store.indices().map(|i| store.record(i).unwrap()).find(|r| r.qname() == name).unwrap()
        };
        let sec = record(b"secStar");
        assert_eq!(sec.seq_len, 0);
        assert_eq!(record(b"qualStar").seq_len, 20);
    }

    // r[verify cram.edge.missing_reference+2]
    /// `test.cram` does not embed its reference, so a FASTA without its
    /// contigs fails the fetch — per slice, once decoding needs the bases.
    #[test]
    fn fasta_without_the_contig_fails_a_slice_that_needs_it() {
        let other_fasta = Path::new(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../tests/data/hts-specs/cram/ce.fa"
        ));
        let mut reader = IndexedCramReader::open(cram_path(), other_fasta).unwrap();
        let mut store = RecordStore::new();
        let tid0 = 0u32.resolve_tid(reader.header()).unwrap();
        let err = reader
            .fetch_into(tid0, (Pos0::new(0).unwrap()..=Pos0::MAX).into(), &mut store)
            .unwrap_err();
        assert!(matches!(err, CramError::MissingReference { .. }), "{err:?}");
    }

    // r[verify cram.edge.missing_reference+2]
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
