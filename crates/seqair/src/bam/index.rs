//! Parse BAI indexes and query them for regions. [`BamIndex::query`] returns merged chunks for a
//! region; [`BamIndex::query_split`] separates nearby (levels 3–5) from distant (levels 0–2) chunks.

use super::bgzf::{BgzfError, VirtualOffset};
use super::csi_index::{CsiIndex, csi_pseudo_bin};
use super::header::BamHeader;
use crate::reader::Tid;
use seqair_types::Pos0;
use std::path::Path;
use tracing::{instrument, warn};

const PSEUDO_BIN: u32 = 37450;
const BAI_MIN_SHIFT: u32 = 14;
const BAI_DEPTH: u32 = 5;

// r[impl io.errors]
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum BaiError {
    #[error("I/O error reading index {path}")]
    Read { path: std::path::PathBuf, source: std::io::Error },

    #[error("invalid BAI magic bytes (expected BAI\\x01)")]
    InvalidMagic,

    #[error("BAI/tabix data is truncated")]
    Truncated,

    #[error("BGZF error decompressing tabix file")]
    Bgzf {
        #[from]
        source: BgzfError,
    },

    #[error("invalid tabix magic bytes: expected TBI\\x01, got {found:?}")]
    InvalidTabixMagic { found: [u8; 4] },

    #[error("tabix header is truncated: {reason}")]
    TruncatedTabixHeader { reason: &'static str },

    #[error("BAI/tabix field has negative count: {value}")]
    NegativeCount { value: i32 },

    #[error("BAI/tabix field `{field}` count {value} exceeds limit {limit}")]
    FieldTooLarge { field: &'static str, value: usize, limit: usize },

    #[error("tabix header names {names} sequences but indexes {n_ref}")]
    TabixNameCount { names: usize, n_ref: usize },
}

#[derive(Debug)]
pub struct BamIndex {
    references: Vec<RefIndex>,
    stats: Vec<ReferenceStats>,
    unplaced_unmapped: Option<u64>,
}

#[derive(Debug)]
struct RefIndex {
    /// Sorted by `bin_id` (stable, so duplicate ids in a malformed file keep
    /// file order), so a query looks up its handful of candidate bins by
    /// binary search instead of testing every bin of the reference against
    /// the candidate list — which, per 10 kb tile on a human chromosome, cost
    /// thousands of bins × ~8 candidates on every query.
    bins: Vec<Bin>,
    linear_index: Vec<VirtualOffset>,
}

impl RefIndex {
    fn new(mut bins: Vec<Bin>, linear_index: Vec<VirtualOffset>) -> Self {
        bins.sort_by_key(|b| b.bin_id);
        RefIndex { bins, linear_index }
    }

    /// The bins whose id is in `candidates`.
    fn candidate_bins<'a>(&'a self, candidates: &'a [u32]) -> impl Iterator<Item = &'a Bin> + 'a {
        candidates.iter().flat_map(move |&id| {
            let lo = self.bins.partition_point(|b| b.bin_id < id);
            self.bins.get(lo..).unwrap_or(&[]).iter().take_while(move |b| b.bin_id == id)
        })
    }
}

#[derive(Debug)]
struct Bin {
    bin_id: u32,
    chunks: Vec<Chunk>,
}

/// A chunk is a pair of virtual offsets defining a contiguous range in the BAM.
#[derive(Debug, Clone, Copy)]
pub struct Chunk {
    pub begin: VirtualOffset,
    pub end: VirtualOffset,
}

// r[impl index.reference_stats]
/// What an index records about one reference, read without touching the
/// alignment file: the counts `samtools idxstats` prints, and how much of the
/// reference holds records.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReferenceStats {
    /// From the reference's pseudo-bin; `None` when the index has none for it
    /// (an index written without them, or a reference with no records).
    pub counts: Option<MappedCounts>,
    /// Finest-level bins with a chunk (a malformed index repeating a bin id
    /// counts it twice). A record spanning more than one
    /// leaf sits in a coarser bin and marks none, so this is a lower bound on
    /// the bins records touch.
    pub covered_leaf_bins: u64,
    /// Bases one leaf bin spans: `2^min_shift`, 16,384 for a BAI.
    pub leaf_bin_len: u64,
}

impl ReferenceStats {
    /// `covered_leaf_bins` in bases. The last bin of a reference may reach past
    /// its end, so cap this at the reference length when dividing by it.
    pub fn covered_bases(&self) -> u64 {
        self.covered_leaf_bins.saturating_mul(self.leaf_bin_len)
    }

    /// Read from one reference's bins, pseudo-bin included, for a binning
    /// scheme of `depth` levels below the root whose leaves span
    /// `2^min_shift` bases.
    pub(super) fn from_bins<'a>(
        bins: impl IntoIterator<Item = (u32, &'a [Chunk])>,
        min_shift: u32,
        depth: u32,
    ) -> Self {
        let pseudo_bin = csi_pseudo_bin(depth);
        // Leaves are the last level: from `((8^depth) - 1) / 7` up to the bin limit.
        let first_leaf =
            1u64.checked_shl(depth.saturating_mul(3)).map_or(u64::MAX, |n| n.saturating_sub(1) / 7);
        let mut counts = None;
        let mut covered_leaf_bins = 0u64;
        for (bin_id, chunks) in bins {
            if bin_id == pseudo_bin {
                // SAM spec §5.2: the second chunk holds (n_mapped, n_unmapped)
                // in place of its two virtual offsets.
                counts =
                    chunks.get(1).map(|c| MappedCounts { mapped: c.begin.0, unmapped: c.end.0 });
            } else if u64::from(bin_id) >= first_leaf
                && bin_id < pseudo_bin.saturating_sub(1)
                && !chunks.is_empty()
            {
                covered_leaf_bins = covered_leaf_bins.saturating_add(1);
            }
        }
        ReferenceStats {
            counts,
            covered_leaf_bins,
            leaf_bin_len: 1u64.checked_shl(min_shift).unwrap_or(u64::MAX),
        }
    }
}

/// A reference's record counts as its index's pseudo-bin keeps them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MappedCounts {
    /// Records without flag 0x4.
    pub mapped: u64,
    /// Records with flag 0x4 placed on this reference (at their mate's
    /// position). Unplaced ones (`RNAME *`) belong to no reference.
    pub unmapped: u64,
}

/// A chunk annotated with its source bin ID.
#[derive(Debug, Clone, Copy)]
pub struct AnnotatedChunk {
    pub chunk: Chunk,
    pub bin_id: u32,
}

// r[impl bam.index.bin0_separate+2]
// r[impl bam.index.chunk_separation+2]
/// Result of a BAM index query, with distant higher-level bin chunks separated
/// from nearby leaf-level chunks.
#[derive(Debug)]
pub struct QueryChunks {
    /// Chunks from bins at levels 3–5 (≤1 Mbp span, clustered near the query).
    /// These are loaded via `RegionBuf` per region.
    pub nearby: Vec<Chunk>,
    /// Chunks from bins at levels 0–2 (≥8 Mbp span, scattered across the file).
    /// These should be cached and reused across regions.
    pub distant: Vec<Chunk>,
}

/// Bin level threshold: bins at this level or lower (covering ≥8 Mbp)
/// are classified as "distant" by `query_split`.
const CACHE_LEVEL_THRESHOLD: u8 = 2;

impl BamIndex {
    // r[impl bam.index.bai_magic]
    // r[impl bam.index.bai_parse]
    #[instrument(level = "debug", fields(path = %path.display()))]
    pub fn from_path(path: &Path) -> Result<Self, BaiError> {
        let data = std::fs::read(path)
            .map_err(|source| BaiError::Read { path: path.to_path_buf(), source })?;

        if data.len() < 8 {
            return Err(BaiError::InvalidMagic);
        }
        if !data.starts_with(b"BAI\x01") {
            return Err(BaiError::InvalidMagic);
        }

        let mut pos = 4;
        parse_refs(&data, &mut pos)
    }

    /// Parse BAI index data from raw bytes (includes magic prefix).
    #[cfg(feature = "fuzz")]
    pub fn from_bytes(data: &[u8]) -> Result<Self, BaiError> {
        if data.len() < 8 {
            return Err(BaiError::InvalidMagic);
        }
        if !data.starts_with(b"BAI\x01") {
            return Err(BaiError::InvalidMagic);
        }
        let mut pos = 4;
        parse_refs(data, &mut pos)
    }

    /// Parse a tabix index from raw (already decompressed) bytes, its
    /// references put in `header`'s order (`r[index.edge.name_mismatch+2]`).
    #[cfg(feature = "fuzz")]
    pub fn from_tabix_bytes(data: &[u8], header: &BamHeader) -> Result<Self, BaiError> {
        Self::from_tabix_data(data, header)
    }

    /// Create an empty index with no references (for plain SAM fuzzing).
    #[cfg(feature = "fuzz")]
    pub fn empty() -> Self {
        BamIndex { references: Vec::new(), stats: Vec::new(), unplaced_unmapped: None }
    }

    // r[impl tabix.bai_reuse]
    // r[impl tabix.index_data]
    // r[impl tabix.compression]
    // r[impl index.shared_query]
    // r[impl index.shared_types]
    /// Parse a tabix (.tbi) index file.
    ///
    /// Tabix files are BGZF-compressed. After decompression and the tabix-specific
    /// header, the bin/chunk/linear-index data is identical to BAI format.
    #[instrument(level = "debug", skip(header), fields(path = %path.display()))]
    pub fn from_tabix_path(path: &Path, header: &BamHeader) -> Result<Self, BaiError> {
        let compressed = std::fs::read(path)
            .map_err(|source| BaiError::Read { path: path.to_path_buf(), source })?;
        let data = decompress_bgzf_file(&compressed)?;
        Self::from_tabix_data(&data, header)
    }

    fn from_tabix_data(data: &[u8], header: &BamHeader) -> Result<Self, BaiError> {
        if data.len() < 4 {
            return Err(BaiError::TruncatedTabixHeader { reason: "file too short for magic" });
        }
        // r[impl tabix.magic]
        if !data.starts_with(b"TBI\x01") {
            let mut found = [0u8; 4];
            if let Some(src) = data.get(..4) {
                found.copy_from_slice(src);
            }
            return Err(BaiError::InvalidTabixMagic { found });
        }

        let mut pos = 4;
        // Read tabix header fields (skip over them — we reuse BAI parsing)
        // r[impl tabix.header]
        let n_ref = read_count(data, &mut pos)?;
        let _format = read_i32(data, &mut pos)?;
        let _col_seq = read_i32(data, &mut pos)?;
        let _col_beg = read_i32(data, &mut pos)?;
        let _col_end = read_i32(data, &mut pos)?;
        let _meta = read_i32(data, &mut pos)?;
        let _skip = read_i32(data, &mut pos)?;
        let l_nm = read_count(data, &mut pos)?;

        let names = data
            .get(pos..pos.wrapping_add(l_nm))
            .ok_or(BaiError::TruncatedTabixHeader { reason: "names extend past end of data" })?;
        pos = pos.wrapping_add(l_nm);
        let names: Vec<&[u8]> = names.split(|&b| b == 0).filter(|n| !n.is_empty()).collect();
        if names.len() != n_ref {
            return Err(BaiError::TabixNameCount { names: names.len(), n_ref });
        }

        // The remaining data is identical to BAI format (bins, chunks, linear index)
        // but we already read n_ref from the tabix header
        let index = parse_refs_with_count(data, &mut pos, n_ref)?;
        Ok(index.in_header_order(&names, header))
    }

    // r[impl index.edge.name_mismatch+2]
    /// Re-key references numbered as `names` lists them (tabix numbers only
    /// the sequences that have records, in file order) by `header`'s tids;
    /// a header sequence the index lacks gets no bins.
    fn in_header_order(self, names: &[&[u8]], header: &BamHeader) -> Self {
        let n = header.target_count();
        let mut references: Vec<RefIndex> =
            std::iter::repeat_with(|| RefIndex::new(Vec::new(), Vec::new())).take(n).collect();
        let mut stats = vec![ReferenceStats::from_bins([], BAI_MIN_SHIFT, BAI_DEPTH); n];
        for ((name, reference), stat) in names.iter().zip(self.references).zip(self.stats) {
            let tid = std::str::from_utf8(name).ok().and_then(|n| header.tid(n));
            match tid.map(Tid::as_usize) {
                Some(tid) => {
                    if let Some(slot) = references.get_mut(tid) {
                        *slot = reference;
                    }
                    if let Some(slot) = stats.get_mut(tid) {
                        *slot = stat;
                    }
                }
                // `tabix -p sam` indexes unplaced lines under RNAME `*`.
                None if *name == b"*" => {}
                None => warn!(
                    name = %name.escape_ascii(),
                    "tabix index names a sequence the SAM header lacks; the index may be stale"
                ),
            }
        }
        BamIndex { references, stats, unplaced_unmapped: self.unplaced_unmapped }
    }

    // r[impl index.reference_stats]
    /// One entry per reference, in header (tid) order.
    pub fn reference_stats(&self) -> &[ReferenceStats] {
        &self.stats
    }

    /// Records with no reference (`RNAME *`, idxstats' `*` line); `None` when
    /// the index does not record them.
    pub fn unplaced_unmapped(&self) -> Option<u64> {
        self.unplaced_unmapped
    }

    // r[impl bam.index.region_query]
    /// Query the index for chunks overlapping a region [start, end] (0-based inclusive).
    // r[impl tabix.query]
    pub fn query(&self, tid: u32, start: Pos0, end: Pos0) -> Vec<Chunk> {
        let result = self.query_split(tid, start, end);
        let mut all = result.nearby;
        all.extend(result.distant);
        all.sort_by_key(|c| c.begin);
        merge_overlapping_chunks(&mut all);
        all
    }

    /// Query the index, returning chunks annotated with their source bin ID.
    /// For diagnostics only — use `query_split` for production code.
    pub fn query_annotated(&self, tid: u32, start: Pos0, end: Pos0) -> Vec<AnnotatedChunk> {
        let Some(ref_idx) = self.references.get(tid as usize) else {
            return Vec::new();
        };
        let start_u64 = start.as_u64();
        let end_u64 = end.as_u64();
        let candidate_bins = reg2bins(start_u64, end_u64.wrapping_add(1));
        let linear_min = linear_index_min(ref_idx, start_u64);
        let mut result = Vec::new();
        for bin in ref_idx.candidate_bins(&candidate_bins) {
            for chunk in &bin.chunks {
                if chunk.end > linear_min {
                    result.push(AnnotatedChunk { chunk: *chunk, bin_id: bin.bin_id });
                }
            }
        }
        result.sort_by_key(|c| c.chunk.begin);
        result
    }

    // r[impl bam.index.bin0_required]
    // r[impl bam.index.bin0_separate+2]
    // r[impl bam.index.chunk_separation+2]
    // r[impl index.edge.no_records]
    /// Query the index, separating distant (level 0–2) from nearby (level 3–5) chunks.
    pub fn query_split(&self, tid: u32, start: Pos0, end: Pos0) -> QueryChunks {
        let Some(ref_idx) = self.references.get(tid as usize) else {
            return QueryChunks { nearby: Vec::new(), distant: Vec::new() };
        };

        let start_u64 = start.as_u64();
        let end_u64 = end.as_u64();
        let candidate_bins = reg2bins(start_u64, end_u64.wrapping_add(1)); // reg2bins uses half-open

        // Linear index: minimum virtual offset for reads starting at or after `start`
        let linear_min = linear_index_min(ref_idx, start_u64);

        let mut nearby = Vec::new();
        let mut distant = Vec::new();
        for bin in ref_idx.candidate_bins(&candidate_bins) {
            let target = if bin_level(bin.bin_id) <= CACHE_LEVEL_THRESHOLD {
                &mut distant
            } else {
                &mut nearby
            };
            for chunk in &bin.chunks {
                // Skip chunks entirely before the linear index minimum
                if chunk.end > linear_min {
                    target.push(*chunk);
                }
            }
        }

        nearby.sort_by_key(|c| c.begin);
        merge_overlapping_chunks(&mut nearby);
        distant.sort_by_key(|c| c.begin);
        merge_overlapping_chunks(&mut distant);
        QueryChunks { nearby, distant }
    }
}

// r[impl csi.unified_enum]
// r[impl csi.query_interface]
/// Unified index type: dispatches queries transparently to BAI or CSI.
#[derive(Debug)]
pub enum AlignmentIndex {
    Bai(BamIndex),
    Csi(CsiIndex),
}

impl AlignmentIndex {
    // r[impl index.reference_stats]
    /// One entry per reference, in header (tid) order.
    pub fn reference_stats(&self) -> &[ReferenceStats] {
        match self {
            Self::Bai(idx) => idx.reference_stats(),
            Self::Csi(idx) => idx.reference_stats(),
        }
    }

    /// Records with no reference (`RNAME *`, idxstats' `*` line); `None` when
    /// the index does not record them.
    pub fn unplaced_unmapped(&self) -> Option<u64> {
        match self {
            Self::Bai(idx) => idx.unplaced_unmapped(),
            Self::Csi(idx) => idx.unplaced_unmapped(),
        }
    }

    /// Query the index for chunks overlapping a region [start, end] (0-based inclusive).
    pub fn query(&self, tid: u32, start: Pos0, end: Pos0) -> Vec<Chunk> {
        match self {
            Self::Bai(idx) => idx.query(tid, start, end),
            Self::Csi(idx) => idx.query(tid, start, end),
        }
    }

    /// Query the index, returning chunks annotated with their source bin ID.
    /// For diagnostics only — use `query_split` for production code.
    pub fn query_annotated(&self, tid: u32, start: Pos0, end: Pos0) -> Vec<AnnotatedChunk> {
        match self {
            Self::Bai(idx) => idx.query_annotated(tid, start, end),
            Self::Csi(idx) => idx.query_annotated(tid, start, end),
        }
    }

    /// Query the index, separating distant from nearby chunks.
    pub fn query_split(&self, tid: u32, start: Pos0, end: Pos0) -> QueryChunks {
        match self {
            Self::Bai(idx) => idx.query_split(tid, start, end),
            Self::Csi(idx) => idx.query_split(tid, start, end),
        }
    }
}

impl From<BamIndex> for AlignmentIndex {
    fn from(idx: BamIndex) -> Self {
        Self::Bai(idx)
    }
}

impl From<CsiIndex> for AlignmentIndex {
    fn from(idx: CsiIndex) -> Self {
        Self::Csi(idx)
    }
}

/// Merge sorted chunks whose byte ranges overlap or are adjacent.
///
/// Chunks from different bins can cover overlapping file regions. Without
/// merging, `fetch_into` would read records in the overlap zone twice.
/// Input must be sorted by `begin`; output preserves that invariant.
pub(super) fn merge_overlapping_chunks(chunks: &mut Vec<Chunk>) {
    if chunks.len() <= 1 {
        return;
    }
    let mut write = 0;
    for read in 1..chunks.len() {
        debug_assert!(read < chunks.len(), "read index OOB: read={read}, len={}", chunks.len());
        debug_assert!(write < read, "write must lag behind read: write={write}, read={read}");
        #[allow(clippy::indexing_slicing, reason = "read < chunks.len() and write < read")]
        if chunks[read].begin <= chunks[write].end {
            // Overlapping or adjacent — extend the current merged chunk
            if chunks[read].end > chunks[write].end {
                chunks[write].end = chunks[read].end;
            }
        } else {
            write = write.wrapping_add(1);
            chunks[write] = chunks[read];
        }
    }
    chunks.truncate(write.wrapping_add(1));
}

/// Parse the reference index data (shared between BAI and tabix).
fn parse_refs(data: &[u8], pos: &mut usize) -> Result<BamIndex, BaiError> {
    let n_ref = read_count(data, pos)?;
    parse_refs_with_count(data, pos, n_ref)
}

// r[impl io.fuzz.alloc_limits]
// r[impl csi.alloc_limits]
// Practical upper bounds for index fields. These exceed any real-world genome
// but catch corrupt/malicious files before they cause multi-GB allocations.
const MAX_INDEX_REFS: usize = 100_000; // human: ~3000 contigs
const MAX_BINS_PER_REF: usize = 100_000; // BAI depth=5: max 37451
const MAX_CHUNKS_PER_BIN: usize = 1_000_000;
const MAX_LINEAR_INDEX: usize = 500_000; // 3 Gbp / 16KiB ≈ 183K

fn check_index_limit(value: usize, limit: usize, field: &'static str) -> Result<(), BaiError> {
    if value > limit {
        return Err(BaiError::FieldTooLarge { field, value, limit });
    }
    Ok(())
}

fn parse_refs_with_count(data: &[u8], pos: &mut usize, n_ref: usize) -> Result<BamIndex, BaiError> {
    check_index_limit(n_ref, MAX_INDEX_REFS, "n_ref")?;
    let mut references = Vec::with_capacity(n_ref);
    let mut stats = Vec::with_capacity(n_ref);

    for _ in 0..n_ref {
        let n_bin = read_count(data, pos)?;
        check_index_limit(n_bin, MAX_BINS_PER_REF, "n_bin")?;
        let mut bins = Vec::with_capacity(n_bin);

        for _ in 0..n_bin {
            let bin_id = read_u32(data, pos)?;
            let n_chunk = read_count(data, pos)?;
            check_index_limit(n_chunk, MAX_CHUNKS_PER_BIN, "n_chunk")?;
            let mut chunks = Vec::with_capacity(n_chunk);

            for _ in 0..n_chunk {
                let begin = VirtualOffset(read_u64(data, pos)?);
                let end = VirtualOffset(read_u64(data, pos)?);
                chunks.push(Chunk { begin, end });
            }
            bins.push(Bin { bin_id, chunks });
        }

        stats.push(ReferenceStats::from_bins(
            bins.iter().map(|b| (b.bin_id, b.chunks.as_slice())),
            BAI_MIN_SHIFT,
            BAI_DEPTH,
        ));
        // r[impl tabix.pseudo_bin]
        // r[impl index.pseudo_bin_dropped]
        // Metadata only, now in `stats`: dropped here, as the CSI parser does,
        // so no query can reach it.
        bins.retain(|b| b.bin_id != PSEUDO_BIN);

        let n_intv = read_count(data, pos)?;
        check_index_limit(n_intv, MAX_LINEAR_INDEX, "n_intv")?;
        let mut linear_index = Vec::with_capacity(n_intv);
        for _ in 0..n_intv {
            linear_index.push(VirtualOffset(read_u64(data, pos)?));
        }

        references.push(RefIndex::new(bins, linear_index));
    }

    let unplaced_unmapped = read_n_no_coor(data, pos, read_u64)?;
    Ok(BamIndex { references, stats, unplaced_unmapped })
}

// r[impl index.unplaced_unmapped]
/// The optional `n_no_coor` after the last reference: records with no
/// reference at all (`RNAME *`). `None` when the file ends without it.
pub(super) fn read_n_no_coor<E>(
    data: &[u8],
    pos: &mut usize,
    read_u64: impl Fn(&[u8], &mut usize) -> Result<u64, E>,
) -> Result<Option<u64>, E> {
    if data.len().saturating_sub(*pos) >= 8 { read_u64(data, pos).map(Some) } else { Ok(None) }
}

/// Decompress all BGZF blocks from a BGZF-compressed file into a single buffer.
pub(super) fn decompress_bgzf_file(compressed: &[u8]) -> Result<Vec<u8>, BaiError> {
    let mut result = Vec::with_capacity(compressed.len().wrapping_mul(3));
    let mut offset = 0;
    let mut decompressor = libdeflater::Decompressor::new();

    while offset < compressed.len() {
        // Check for BGZF header
        let remaining = compressed.len().wrapping_sub(offset);
        if remaining < 18 {
            break; // EOF or truncated block
        }

        debug_assert!(remaining >= 18, "BGZF block too short: remaining={remaining}");
        #[expect(clippy::indexing_slicing, reason = "remaining >= 18 checked above")]
        let block = &compressed[offset..];

        // Verify gzip magic
        if block.get(..2) != Some(&[0x1f, 0x8b]) {
            return Err(BaiError::Bgzf { source: BgzfError::InvalidMagic });
        }

        // BSIZE from the BGZF extra field at offset 16-17 (after standard 12-byte gzip
        // header + 6-byte extra field header)
        #[expect(clippy::indexing_slicing, reason = "remaining >= 18 checked above")]
        let bsize = (u16::from_le_bytes([block[16], block[17]]) as usize).wrapping_add(1);

        if offset.wrapping_add(bsize) > compressed.len() {
            return Err(BaiError::Bgzf { source: BgzfError::TruncatedBlock });
        }

        // r[impl bam.index.bgzf_block_validation]
        // ISIZE (uncompressed size) from the last 4 bytes of the block
        if bsize < 18 {
            #[expect(clippy::cast_possible_truncation, reason = "bsize < 18")]
            return Err(BaiError::Bgzf {
                source: BgzfError::BlockSizeTooSmall { bsize: bsize as u16 },
            });
        }
        debug_assert!(bsize >= 18, "bsize should be >= 18 at this point: bsize={bsize}");
        #[expect(clippy::indexing_slicing, reason = "bsize >= 18 checked above")]
        let isize_bytes = &block[bsize.wrapping_sub(4)..bsize];
        let isize = u32::from_le_bytes(
            isize_bytes.try_into().unwrap_or_else(|_| unreachable!("bsize >= 18")),
        ) as usize;

        if isize == 0 {
            offset = offset.wrapping_add(bsize);
            continue; // EOF block
        }

        // Compressed data starts at offset 18, ends at bsize - 8 (before CRC32 + ISIZE)
        if bsize < 26 {
            #[expect(clippy::cast_possible_truncation, reason = "bsize < 26")]
            return Err(BaiError::Bgzf {
                source: BgzfError::BlockSizeTooSmall { bsize: bsize as u16 },
            });
        }
        debug_assert!(bsize >= 26, "bsize should be >= 26 at this point: bsize={bsize}");
        #[allow(clippy::indexing_slicing, reason = "bsize >= 26 checked above")]
        let cdata = &block[18..bsize.wrapping_sub(8)];

        let start = result.len();
        result.resize(start.wrapping_add(isize), 0);
        debug_assert!(
            start.wrapping_add(isize) <= result.len(),
            "result overrun: {} > {}",
            start.wrapping_add(isize),
            result.len()
        );
        #[allow(clippy::indexing_slicing, reason = "just resized to start + isize")]
        decompressor
            .deflate_decompress(cdata, &mut result[start..start.wrapping_add(isize)])
            .map_err(|source| BaiError::Bgzf {
                source: BgzfError::DecompressionFailed { source },
            })?;

        offset = offset.wrapping_add(bsize);
    }

    Ok(result)
}

// r[impl bam.index.bin_calculation]
/// BAI bin level (0 = root/512 Mbp, 5 = leaf/16 Kbp).
pub fn bin_level(bin_id: u32) -> u8 {
    match bin_id {
        0 => 0,
        1..=8 => 1,
        9..=72 => 2,
        73..=584 => 3,
        585..=4680 => 4,
        4681..=37449 => 5,
        _ => 6, // pseudo-bin or invalid
    }
}

/// Genomic span covered by one bin at a given level.
pub fn bin_span(level: u8) -> u64 {
    match level {
        0 => 1 << 29, // 512 Mbp
        1 => 1 << 26, // 64 Mbp
        2 => 1 << 23, // 8 Mbp
        3 => 1 << 20, // 1 Mbp
        4 => 1 << 17, // 128 Kbp
        5 => 1 << 14, // 16 Kbp
        _ => 0,
    }
}

fn linear_index_min(ref_idx: &RefIndex, start: u64) -> VirtualOffset {
    let window = (start / 16384) as usize;
    ref_idx.linear_index.get(window).copied().unwrap_or(VirtualOffset(0))
}

/// Compute all BAI bins that overlap the interval [beg, end) (0-based, half-open).
fn reg2bins(beg: u64, end: u64) -> Vec<u32> {
    let mut bins = Vec::with_capacity(32);
    bins.push(0); // bin 0 always included

    let levels: [(u32, u32); 5] = [
        (1, 26),    // level 1: 64 Mbp bins
        (9, 23),    // level 2: 8 Mbp bins
        (73, 20),   // level 3: 1 Mbp bins
        (585, 17),  // level 4: 128 Kbp bins
        (4681, 14), // level 5: 16 Kbp bins (leaf)
    ];

    for &(offset, shift) in &levels {
        #[expect(
            clippy::cast_possible_truncation,
            reason = "BAM positions are capped at 2^29; shifted by ≥14 gives values ≤ 2^15, fits in u32"
        )]
        let mut k = offset.wrapping_add((beg >> shift) as u32);
        #[expect(
            clippy::cast_possible_truncation,
            reason = "BAM positions are capped at 2^29; shifted by ≥14 gives values ≤ 2^15, fits in u32"
        )]
        let end_k = offset.wrapping_add((end.saturating_sub(1) >> shift) as u32);
        while k <= end_k {
            bins.push(k);
            k = k.wrapping_add(1);
        }
    }

    bins
}

fn read_i32(data: &[u8], pos: &mut usize) -> Result<i32, BaiError> {
    let bytes: [u8; 4] = data
        .get(*pos..pos.wrapping_add(4))
        .and_then(|s| s.try_into().ok())
        .ok_or(BaiError::Truncated)?;
    *pos = pos.wrapping_add(4);
    Ok(i32::from_le_bytes(bytes))
}

// r[impl bam.index.negative_count]
/// An `int32_t` count: a negative one is an error of its own, not a huge
/// `usize` that some later limit happens to catch.
fn read_count(data: &[u8], pos: &mut usize) -> Result<usize, BaiError> {
    let raw = read_i32(data, pos)?;
    usize::try_from(raw).map_err(|_| BaiError::NegativeCount { value: raw })
}

fn read_u32(data: &[u8], pos: &mut usize) -> Result<u32, BaiError> {
    let bytes: [u8; 4] = data
        .get(*pos..pos.wrapping_add(4))
        .and_then(|s| s.try_into().ok())
        .ok_or(BaiError::Truncated)?;
    *pos = pos.wrapping_add(4);
    Ok(u32::from_le_bytes(bytes))
}

fn read_u64(data: &[u8], pos: &mut usize) -> Result<u64, BaiError> {
    let bytes: [u8; 8] = data
        .get(*pos..pos.wrapping_add(8))
        .and_then(|s| s.try_into().ok())
        .ok_or(BaiError::Truncated)?;
    *pos = pos.wrapping_add(8);
    Ok(u64::from_le_bytes(bytes))
}

#[cfg(test)]
#[allow(clippy::arithmetic_side_effects, reason = "tests")]
#[allow(clippy::cast_possible_truncation, reason = "tests")]
mod tests {
    use super::*;
    use seqair_types::Pos0;

    #[test]
    fn reg2bins_includes_bin0() {
        let bins = reg2bins(0, 100);
        assert!(bins.contains(&0));
    }

    #[test]
    fn reg2bins_small_region() {
        let bins = reg2bins(100, 200);
        assert!(bins.contains(&0));
        // Should include leaf-level bin for positions 100-200
        assert!(bins.contains(&(4681 + (100 >> 14) as u32)));
    }

    // --- query_split tests ---

    /// Build a minimal `BamIndex` with specific bins populated.
    fn index_with_bins(bins: Vec<(u32, Vec<Chunk>)>) -> BamIndex {
        BamIndex {
            references: vec![RefIndex::new(
                bins.into_iter().map(|(bin_id, chunks)| Bin { bin_id, chunks }).collect(),
                vec![VirtualOffset(0); 64],
            )],
            stats: Vec::new(),
            unplaced_unmapped: None,
        }
    }

    fn chunk(begin: u64, end: u64) -> Chunk {
        Chunk { begin: VirtualOffset(begin), end: VirtualOffset(end) }
    }

    // r[verify bam.index.bin0_separate+2]
    // r[verify bam.index.chunk_separation+2]
    #[test]
    fn query_split_separates_bin0() {
        let idx = index_with_bins(vec![
            (0, vec![chunk(9000, 9500)]),  // bin 0: distant
            (4681, vec![chunk(100, 500)]), // leaf bin: main data
        ]);

        let result = idx.query_split(0, Pos0::new(0).unwrap(), Pos0::new(100).unwrap());

        assert_eq!(result.distant.len(), 1, "bin 0 should have 1 chunk");
        assert_eq!(result.distant[0].begin.0, 9000);
        assert!(!result.nearby.is_empty(), "regular should have chunks");
        assert!(
            result.nearby.iter().all(|c| c.begin.0 != 9000),
            "bin 0 chunk should not appear in regular"
        );
    }

    // r[verify bam.index.bin0_required]
    #[test]
    fn query_split_union_equals_query() {
        let idx = index_with_bins(vec![
            (0, vec![chunk(9000, 9500)]),
            (4681, vec![chunk(100, 500)]),
            (4682, vec![chunk(500, 900)]),
        ]);

        let split = idx.query_split(0, Pos0::new(0).unwrap(), Pos0::new(200).unwrap());
        let flat = idx.query(0, Pos0::new(0).unwrap(), Pos0::new(200).unwrap());

        let mut combined: Vec<u64> =
            split.nearby.iter().chain(&split.distant).map(|c| c.begin.0).collect();
        combined.sort();

        let mut flat_begins: Vec<u64> = flat.iter().map(|c| c.begin.0).collect();
        flat_begins.sort();

        assert_eq!(combined, flat_begins, "split union must equal flat query");
    }

    #[test]
    fn query_split_no_bin0_in_index() {
        let idx = index_with_bins(vec![(4681, vec![chunk(100, 500)])]);

        let result = idx.query_split(0, Pos0::new(0).unwrap(), Pos0::new(100).unwrap());
        assert!(result.distant.is_empty(), "no bin 0 in index → empty bin0");
        assert!(!result.nearby.is_empty());
    }

    // r[verify index.edge.no_records]
    #[test]
    fn query_split_empty_reference() {
        let idx = BamIndex { references: Vec::new(), stats: Vec::new(), unplaced_unmapped: None };
        let result = idx.query_split(0, Pos0::new(0).unwrap(), Pos0::new(100).unwrap());
        assert!(result.nearby.is_empty());
        assert!(result.distant.is_empty());
    }

    // --- Bin 0 assignment: reads at 64 MB boundaries ---

    #[test]
    fn reg2bin_at_64mb_boundary_returns_bin0() {
        // A read spanning the 64 MB boundary (2^26 = 67_108_864)
        // beg in one level-1 bin, end in another → must be bin 0
        let boundary: u64 = 1 << 26; // 67_108_864
        let beg = boundary - 10;
        let end = boundary + 10;

        let bins = reg2bins(beg, end);
        assert!(bins.contains(&0), "read spanning 64 MB boundary must include bin 0");

        // Verify level-1 bins differ for beg and end
        let beg_l1 = 1 + (beg >> 26) as u32;
        let end_l1 = 1 + ((end - 1) >> 26) as u32;
        assert_ne!(beg_l1, end_l1, "beg and end should be in different level-1 bins");
    }

    #[test]
    fn reg2bin_not_at_boundary_still_has_bin0_in_candidates() {
        // Even a small region not near a boundary includes bin 0 in the
        // candidate list — but there just won't be any bin 0 chunks in the index
        // for that region (the linear index filters them).
        let bins = reg2bins(1000, 1100);
        assert!(bins.contains(&0));
    }

    #[test]
    fn from_path_invalid_magic_returns_error() {
        use std::io::Write;
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("bad.bai");
        let mut f = std::fs::File::create(&path).unwrap();
        f.write_all(b"NOTBAI\x01\x00extra").unwrap();
        let err = BamIndex::from_path(&path).unwrap_err();
        assert!(matches!(err, BaiError::InvalidMagic));
    }

    #[test]
    fn from_path_too_short_returns_error() {
        use std::io::Write;
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("short.bai");
        let mut f = std::fs::File::create(&path).unwrap();
        f.write_all(b"BAI").unwrap(); // only 3 bytes, < 8
        let err = BamIndex::from_path(&path).unwrap_err();
        assert!(matches!(err, BaiError::InvalidMagic));
    }

    // r[verify index.pseudo_bin_dropped]
    // r[verify tabix.pseudo_bin]
    #[test]
    fn the_pseudo_bin_is_read_into_stats_and_dropped() {
        let mut bai = b"BAI\x01".to_vec();
        bai.extend_from_slice(&1i32.to_le_bytes()); // n_ref
        bai.extend_from_slice(&2i32.to_le_bytes()); // n_bin
        for (bin, chunks) in
            [(4681u32, [(10u64, 20u64)].as_slice()), (PSEUDO_BIN, &[(10, 20), (7, 3)])]
        {
            bai.extend_from_slice(&bin.to_le_bytes());
            bai.extend_from_slice(&i32::try_from(chunks.len()).unwrap().to_le_bytes());
            for (begin, end) in chunks {
                bai.extend_from_slice(&begin.to_le_bytes());
                bai.extend_from_slice(&end.to_le_bytes());
            }
        }
        bai.extend_from_slice(&0i32.to_le_bytes()); // n_intv
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("pseudo.bai");
        std::fs::write(&path, &bai).unwrap();

        let index = BamIndex::from_path(&path).unwrap();
        assert_eq!(
            index.reference_stats(),
            [ReferenceStats {
                counts: Some(MappedCounts { mapped: 7, unmapped: 3 }),
                covered_leaf_bins: 1,
                leaf_bin_len: 16_384,
            }]
        );
        let ids: Vec<u32> = index.references[0].bins.iter().map(|b| b.bin_id).collect();
        assert_eq!(ids, [4681]);
    }

    // r[verify bam.index.negative_count]
    /// One reference with one bin of one chunk and no linear index, with the
    /// count at `field` (0 `n_ref`, 1 `n_bin`, 2 `n_chunk`, 3 `n_intv`) set to -1.
    #[test]
    fn a_negative_count_is_reported_as_negative() {
        for field in 0..4 {
            let count = |i: usize, n: i32| if i == field { -1 } else { n };
            let mut bai = b"BAI\x01".to_vec();
            bai.extend_from_slice(&count(0, 1).to_le_bytes());
            bai.extend_from_slice(&count(1, 1).to_le_bytes());
            bai.extend_from_slice(&4681u32.to_le_bytes());
            bai.extend_from_slice(&count(2, 1).to_le_bytes());
            bai.extend_from_slice(&[0; 16]);
            bai.extend_from_slice(&count(3, 0).to_le_bytes());
            let dir = tempfile::tempdir().unwrap();
            let path = dir.path().join("negative.bai");
            std::fs::write(&path, &bai).unwrap();
            let err = BamIndex::from_path(&path).unwrap_err();
            assert!(matches!(err, BaiError::NegativeCount { value: -1 }), "field {field}: {err:?}");
        }
    }

    // r[verify tabix.magic]
    #[test]
    fn from_tabix_path_wrong_magic_returns_invalid_tabix_magic() {
        use std::io::Write;

        // Craft a minimal valid BGZF block whose decompressed content is
        // [0xDE, 0xAD, 0xBE, 0xEF] — valid BGZF structure, wrong tabix magic.
        //
        // BGZF block layout (35 bytes total):
        //   [0..2]   gzip magic: 0x1f 0x8b
        //   [2]      CM=8 (deflate)
        //   [3]      FLG=0x04 (FEXTRA)
        //   [4..8]   MTIME=0
        //   [8]      XFL=0
        //   [9]      OS=0xff
        //   [10..12] XLEN=6
        //   [12..14] SI1='B' SI2='C'
        //   [14..16] SLEN=2
        //   [16..18] BSIZE-1 = 34 (stored value; actual block size = 35)
        //   [18..27] deflate stored block: BFINAL=1 BTYPE=00
        //            + LEN(4) + NLEN(~4) + 4 payload bytes
        //   [27..31] CRC32 (not verified by decompress_bgzf_file)
        //   [31..35] ISIZE=4
        //
        // Followed by the standard 28-byte BGZF EOF block.
        #[rustfmt::skip]
        let bgzf_block: &[u8] = &[
            // gzip header
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff,
            // extra field: XLEN=6, BC subfield, SLEN=2, BSIZE-1=34
            0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x22, 0x00,
            // deflate stored block (BFINAL=1, BTYPE=00): LEN=4, NLEN=~4, data
            0x01, 0x04, 0x00, 0xfb, 0xff, 0xde, 0xad, 0xbe, 0xef,
            // CRC32 (unchecked)
            0x00, 0x00, 0x00, 0x00,
            // ISIZE=4
            0x04, 0x00, 0x00, 0x00,
        ];
        // Standard BGZF EOF block (SAM spec section 4.1)
        #[rustfmt::skip]
        let bgzf_eof: &[u8] = &[
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff,
            0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00,
            0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        ];

        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("bad_magic.tbi");
        let mut f = std::fs::File::create(&path).unwrap();
        f.write_all(bgzf_block).unwrap();
        f.write_all(bgzf_eof).unwrap();
        drop(f);

        let err = BamIndex::from_tabix_path(&path, &BamHeader::synthetic(1, 1000)).unwrap_err();
        match err {
            BaiError::InvalidTabixMagic { found } => {
                assert_eq!(found, [0xde, 0xad, 0xbe, 0xef]);
            }
            other => panic!("expected InvalidTabixMagic, got {other:?}"),
        }
    }

    // r[verify bam.index.bgzf_block_validation]
    #[test]
    fn decompress_bgzf_rejects_bsize_too_small() {
        // Build a malformed BGZF block with bsize = 17 (< 18 minimum).
        // Valid gzip magic, but BSIZE field claiming a very small block.
        let mut block = vec![0x1f, 0x8b, 0x08, 0x04]; // gzip magic
        block.extend_from_slice(&[0; 4]); // MTIME
        block.push(0); // XFL
        block.push(0xff); // OS
        block.extend_from_slice(&6u16.to_le_bytes()); // XLEN = 6
        block.extend_from_slice(&[b'B', b'C', 2, 0]); // BC subfield
        // BSIZE = 16 means total block size = 17, which is < 18
        block.extend_from_slice(&16u16.to_le_bytes());
        // Pad to make the block at least 18 bytes for the outer remaining check
        block.resize(18, 0);

        let result = decompress_bgzf_file(&block);
        assert!(result.is_err(), "should reject bsize < 18");
        let err = result.unwrap_err();
        assert!(
            matches!(err, BaiError::Bgzf { source: BgzfError::BlockSizeTooSmall { .. } }),
            "expected BlockSizeTooSmall, got {err:?}"
        );
    }

    /// Chunks from different nearby bins (e.g. level 3 and level 5) can have
    /// overlapping byte ranges in the BAM file. `query_split` must merge these
    /// so that `fetch_into` never reads the same file region twice.
    #[test]
    fn query_split_merges_overlapping_nearby_chunks() {
        // bin 73 (level 3) and bin 4685 (level 5) both match region [65500, 65700].
        // Their chunks overlap in file space: bin 73 ends at 800, bin 4685 starts at 700.
        let idx = index_with_bins(vec![
            (73, vec![chunk(500, 800)]),    // level 3 nearby
            (4685, vec![chunk(700, 1000)]), // level 5 nearby — overlaps with bin 73
        ]);

        let result = idx.query_split(0, Pos0::new(65500).unwrap(), Pos0::new(65700).unwrap());

        // After merging, there should be a single chunk [500, 1000]
        assert_eq!(
            result.nearby.len(),
            1,
            "overlapping nearby chunks must be merged, got: {:?}",
            result.nearby.iter().map(|c| (c.begin.0, c.end.0)).collect::<Vec<_>>()
        );
        assert_eq!(result.nearby[0].begin.0, 500);
        assert_eq!(result.nearby[0].end.0, 1000);
    }

    #[test]
    fn query_split_merges_adjacent_nearby_chunks() {
        // Two chunks that are exactly adjacent (no gap, no overlap)
        let idx =
            index_with_bins(vec![(73, vec![chunk(500, 700)]), (4685, vec![chunk(700, 1000)])]);

        let result = idx.query_split(0, Pos0::new(65500).unwrap(), Pos0::new(65700).unwrap());

        assert_eq!(
            result.nearby.len(),
            1,
            "adjacent chunks must be merged, got: {:?}",
            result.nearby.iter().map(|c| (c.begin.0, c.end.0)).collect::<Vec<_>>()
        );
    }

    #[test]
    fn query_split_keeps_non_overlapping_nearby_chunks_separate() {
        let idx =
            index_with_bins(vec![(73, vec![chunk(500, 700)]), (4685, vec![chunk(900, 1200)])]);

        let result = idx.query_split(0, Pos0::new(65500).unwrap(), Pos0::new(65700).unwrap());

        assert_eq!(
            result.nearby.len(),
            2,
            "non-overlapping chunks must stay separate, got: {:?}",
            result.nearby.iter().map(|c| (c.begin.0, c.end.0)).collect::<Vec<_>>()
        );
    }

    // --- merge_overlapping_chunks unit tests ---

    fn merge(pairs: &[(u64, u64)]) -> Vec<(u64, u64)> {
        let mut chunks: Vec<Chunk> = pairs.iter().map(|&(b, e)| chunk(b, e)).collect();
        chunks.sort_by_key(|c| c.begin);
        merge_overlapping_chunks(&mut chunks);
        chunks.iter().map(|c| (c.begin.0, c.end.0)).collect()
    }

    #[test]
    fn merge_empty() {
        assert_eq!(merge(&[]), Vec::<(u64, u64)>::new());
    }

    #[test]
    fn merge_single() {
        assert_eq!(merge(&[(10, 20)]), vec![(10, 20)]);
    }

    #[test]
    fn merge_subset_chunk() {
        // B is entirely inside A
        assert_eq!(merge(&[(10, 100), (30, 50)]), vec![(10, 100)]);
    }

    #[test]
    fn merge_chain_of_three() {
        // A overlaps B, B overlaps C → all merge into one
        assert_eq!(merge(&[(10, 30), (20, 50), (40, 70)]), vec![(10, 70)]);
    }

    #[test]
    fn merge_mixed_overlap_and_gap() {
        assert_eq!(
            merge(&[(10, 30), (20, 50), (100, 200), (150, 250)]),
            vec![(10, 50), (100, 250)]
        );
    }

    #[test]
    fn merge_identical_chunks() {
        assert_eq!(merge(&[(10, 20), (10, 20), (10, 20)]), vec![(10, 20)]);
    }

    // --- merge output is non-overlapping and covers the same range ---

    use hegel::prelude::*;

    /// A chunk `(begin, end)` with `begin < end`, spanning up to 5000 units.
    #[hegel::composite]
    fn arb_chunk(tc: &TestCase) -> (u64, u64) {
        let begin = tc.draw(gs::integers::<u64>().max_value(9_999));
        let end = tc.draw(gs::integers::<u64>().min_value(begin + 1).max_value(begin + 4_999));
        (begin, end)
    }

    #[hegel::test]
    fn merge_output_is_sorted_and_non_overlapping(tc: TestCase) {
        let chunks = tc.draw(gs::vecs(arb_chunk()).max_size(19));
        let merged = merge(&chunks);
        // Sorted
        for w in merged.windows(2) {
            assert!(w[0].0 < w[1].0, "not sorted: {merged:?}");
        }
        // Non-overlapping: each chunk's begin is strictly after the previous end
        for w in merged.windows(2) {
            assert!(w[1].0 > w[0].1, "overlapping after merge: {merged:?}");
        }
    }

    /// The merged list must be the input's union — the same bytes, no more.
    ///
    /// The oracle marks every byte each input chunk covers in a flat array and
    /// reads the runs back off it, which shares no arithmetic with the in-place
    /// two-pointer sweep under test. Checking only that the input is *covered*
    /// (as this test used to) passes for a `merge` that returns one chunk
    /// spanning everything and reads the whole file back.
    #[hegel::test]
    fn merge_is_exactly_the_union_of_its_input(tc: TestCase) {
        let chunks = tc.draw(gs::vecs(arb_chunk()).min_size(1).max_size(19));

        // `arb_chunk` bounds begin < 10_000 and end < begin + 5_000.
        const SPAN: usize = 15_000;
        let mut covered = vec![false; SPAN];
        for &(b, e) in &chunks {
            for byte in covered.get_mut(b as usize..e as usize).expect("bounded by arb_chunk") {
                *byte = true;
            }
        }
        let mut expected: Vec<(u64, u64)> = Vec::new();
        let mut run_start: Option<usize> = None;
        for i in 0..=SPAN {
            match (run_start, covered.get(i).copied().unwrap_or(false)) {
                (None, true) => run_start = Some(i),
                (Some(start), false) => {
                    expected.push((start as u64, i as u64));
                    run_start = None;
                }
                _ => {}
            }
        }

        assert_eq!(merge(&chunks), expected);
    }

    /// A bin id near a region's candidates: any level-0..5 bin over the
    /// first 4 Mbp, or the pseudo-bin.
    #[hegel::composite]
    fn arb_bin_id(tc: &TestCase) -> u32 {
        let pick = tc.draw(gs::integers::<u8>().max_value(6));
        let (offset, shift): (u32, u32) = match pick {
            0 => return 0,
            1 => (1, 26),
            2 => (9, 23),
            3 => (73, 20),
            4 => (585, 17),
            5 => (4681, 14),
            _ => return PSEUDO_BIN,
        };
        offset + tc.draw(gs::integers::<u32>().max_value((1 << 22) >> shift))
    }

    // r[verify bam.index.chunk_separation+2]
    /// The sorted-bin lookup returns what testing every bin of the reference
    /// against the candidate list returns — bins in any file order, duplicate
    /// ids included.
    #[hegel::test]
    fn query_split_matches_full_bin_scan(tc: TestCase) {
        let bins: Vec<(u32, Vec<Chunk>)> = tc
            .draw(
                gs::vecs(gs::tuples!(arb_bin_id(), gs::vecs(arb_chunk()).max_size(3))).max_size(40),
            )
            .into_iter()
            .map(|(id, cs)| (id, cs.into_iter().map(|(b, e)| chunk(b, e)).collect()))
            .collect();
        let start = tc.draw(gs::integers::<u32>().max_value(1 << 22));
        let end = tc.draw(gs::integers::<u32>().min_value(start).max_value(start + 200_000));

        let candidates = reg2bins(u64::from(start), u64::from(end) + 1);
        let (mut nearby, mut distant) = (Vec::new(), Vec::new());
        for (id, chunks) in &bins {
            if *id == PSEUDO_BIN || !candidates.contains(id) {
                continue;
            }
            let target =
                if bin_level(*id) <= CACHE_LEVEL_THRESHOLD { &mut distant } else { &mut nearby };
            target.extend(chunks.iter().filter(|c| c.end > VirtualOffset(0)).copied());
        }
        let key = |v: &mut Vec<Chunk>| {
            v.sort_by_key(|c| c.begin);
            merge_overlapping_chunks(v);
            v.iter().map(|c| (c.begin.0, c.end.0)).collect::<Vec<_>>()
        };

        let got = index_with_bins(bins).query_split(
            0,
            Pos0::new(start).unwrap(),
            Pos0::new(end).unwrap(),
        );
        let (mut got_near, mut got_far) = (got.nearby, got.distant);
        assert_eq!(key(&mut got_near), key(&mut nearby));
        assert_eq!(key(&mut got_far), key(&mut distant));
    }
}
