//! Parse CSI (Coordinate-Sorted Index) files and query them for regions.
//! CSI generalizes BAI with configurable `min_shift`/`depth` and per-bin `loffset`
//! instead of a separate linear index.

use super::{
    bgzf::{BgzfError, VirtualOffset},
    index::{
        AnnotatedChunk, BaiError, Chunk, QueryChunks, decompress_bgzf_file,
        merge_overlapping_chunks,
    },
};
use seqair_types::Pos0;
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use tracing::instrument;

// r[impl csi.alloc_limits]
const MAX_INDEX_REFS: usize = 100_000;
const MAX_BINS_PER_REF: usize = 100_000;
const MAX_CHUNKS_PER_BIN: usize = 1_000_000;

// r[impl csi.header_bounds]
// depth * 3 is used as a shift amount on u64 values. depth=17 with min_shift=14
// gives shift=65 which overflows. depth=16 with min_shift=1 gives shift=49, safe.
// No real genome needs depth > 10 (depth 10, min_shift 14 → 2^44 = 17.6 Tbp).
const MAX_CSI_DEPTH: u32 = 16;

// `min_shift + depth * 3` is used as a shift amount on u64 values. Rust masks
// shift amounts by the type width, so `x >> 64` silently becomes `x >> 0`.
// We require `min_shift + depth * 3 < 64`, i.e. `min_shift <= 63 - depth*3`.
// This is the only constraint that protects `csi_reg2bins` and `csi_min_offset`
// from silently wrong bin calculations on a malformed file.
const MAX_CSI_SHIFT_SUM: u32 = 63;

// r[impl csi.errors]
#[derive(Debug, thiserror::Error)]
pub enum CsiError {
    #[error("I/O error reading CSI index {path}")]
    Read { path: PathBuf, source: std::io::Error },

    #[error("invalid CSI magic bytes (expected CSI\\x01)")]
    InvalidMagic,

    #[error("CSI data is truncated")]
    Truncated,

    #[error("CSI field has negative count: {value}")]
    NegativeCount { value: i32 },

    #[error("CSI field `{field}` count {value} exceeds limit {limit}")]
    FieldTooLarge { field: &'static str, value: usize, limit: usize },

    #[error("CSI depth {depth} exceeds maximum {max}")]
    DepthTooLarge { depth: u32, max: u32 },

    #[error(
        "CSI min_shift {min_shift} with depth {depth} produces shift sum {sum} which would \
         overflow the u64 shift width (must satisfy min_shift + depth*3 <= 63)"
    )]
    MinShiftTooLarge { min_shift: u32, depth: u32, sum: u64 },

    #[error("BGZF error decompressing CSI file")]
    Bgzf {
        #[from]
        source: BgzfError,
    },
}

fn check_limit(value: usize, limit: usize, field: &'static str) -> Result<(), CsiError> {
    if value > limit {
        return Err(CsiError::FieldTooLarge { field, value, limit });
    }
    Ok(())
}

fn read_i32(data: &[u8], pos: &mut usize) -> Result<i32, CsiError> {
    let bytes: [u8; 4] = data
        .get(*pos..pos.wrapping_add(4))
        .and_then(|s| s.try_into().ok())
        .ok_or(CsiError::Truncated)?;
    *pos = pos.wrapping_add(4);
    Ok(i32::from_le_bytes(bytes))
}

fn read_u32(data: &[u8], pos: &mut usize) -> Result<u32, CsiError> {
    let bytes: [u8; 4] = data
        .get(*pos..pos.wrapping_add(4))
        .and_then(|s| s.try_into().ok())
        .ok_or(CsiError::Truncated)?;
    *pos = pos.wrapping_add(4);
    Ok(u32::from_le_bytes(bytes))
}

fn read_u64(data: &[u8], pos: &mut usize) -> Result<u64, CsiError> {
    let bytes: [u8; 8] = data
        .get(*pos..pos.wrapping_add(8))
        .and_then(|s| s.try_into().ok())
        .ok_or(CsiError::Truncated)?;
    *pos = pos.wrapping_add(8);
    Ok(u64::from_le_bytes(bytes))
}

/// Read an i32 count field, rejecting negative values before casting to usize.
fn read_count(data: &[u8], pos: &mut usize) -> Result<usize, CsiError> {
    let raw = read_i32(data, pos)?;
    if raw < 0 {
        return Err(CsiError::NegativeCount { value: raw });
    }
    Ok(raw as usize)
}

// r[impl csi.header]
/// A CSI index with configurable binning parameters.
#[derive(Debug)]
pub struct CsiIndex {
    min_shift: u32,
    depth: u32,
    aux: Vec<u8>,
    references: Vec<CsiRefIndex>,
}

#[derive(Debug)]
struct CsiRefIndex {
    bins: Vec<CsiBin>,
    /// Indexes into `bins` keyed by `bin_id`. Built once at parse time so that
    /// `csi_min_offset` does not have to do a linear scan per leaf-to-root
    /// step (relevant for dense indexes on large chromosomes).
    by_bin_id: HashMap<u32, usize>,
}

impl CsiRefIndex {
    fn new(bins: Vec<CsiBin>) -> Self {
        let by_bin_id = bins.iter().enumerate().map(|(i, b)| (b.bin_id, i)).collect();
        CsiRefIndex { bins, by_bin_id }
    }
}

#[derive(Debug)]
struct CsiBin {
    bin_id: u32,
    loffset: VirtualOffset,
    chunks: Vec<Chunk>,
}

impl CsiIndex {
    // r[impl csi.magic]
    // r[impl csi.header]
    // r[impl csi.aux_data]
    // r[impl csi.n_ref]
    // r[impl csi.bins]
    #[instrument(level = "debug", fields(path = %path.display()))]
    pub fn from_path(path: &Path) -> Result<Self, CsiError> {
        let raw = std::fs::read(path)
            .map_err(|source| CsiError::Read { path: path.to_path_buf(), source })?;

        // CSI files may be BGZF-compressed (samtools always compresses them).
        // Detect by checking for gzip magic (0x1f 0x8b) vs CSI magic (CSI\x01).
        let data = if raw.starts_with(&[0x1f, 0x8b]) {
            decompress_bgzf_file(&raw).map_err(|e| match e {
                BaiError::Bgzf { source } => CsiError::Bgzf { source },
                _ => CsiError::Truncated,
            })?
        } else {
            raw
        };

        Self::from_bytes_inner(&data)
    }

    /// Parse CSI index data from raw bytes.
    #[cfg(feature = "fuzz")]
    pub fn from_bytes(data: &[u8]) -> Result<Self, CsiError> {
        Self::from_bytes_inner(data)
    }

    fn from_bytes_inner(data: &[u8]) -> Result<Self, CsiError> {
        if data.len() < 4 {
            return Err(CsiError::InvalidMagic);
        }
        if !data.starts_with(b"CSI\x01") {
            return Err(CsiError::InvalidMagic);
        }

        let mut pos = 4;
        let min_shift = read_i32(data, &mut pos)?;
        if min_shift < 0 {
            return Err(CsiError::NegativeCount { value: min_shift });
        }
        let min_shift = min_shift as u32;

        let depth = read_i32(data, &mut pos)?;
        if depth < 0 {
            return Err(CsiError::NegativeCount { value: depth });
        }
        let depth = depth as u32;
        // r[impl csi.header_bounds]
        if depth > MAX_CSI_DEPTH {
            return Err(CsiError::DepthTooLarge { depth, max: MAX_CSI_DEPTH });
        }
        // r[impl csi.header_bounds]
        let shift_sum = u64::from(min_shift).saturating_add(u64::from(depth).saturating_mul(3));
        if shift_sum > u64::from(MAX_CSI_SHIFT_SUM) {
            return Err(CsiError::MinShiftTooLarge { min_shift, depth, sum: shift_sum });
        }

        // r[impl csi.aux_data]
        let l_aux = read_i32(data, &mut pos)?;
        if l_aux < 0 {
            return Err(CsiError::NegativeCount { value: l_aux });
        }
        let l_aux = l_aux as usize;
        if data.len() < pos.wrapping_add(l_aux) {
            return Err(CsiError::Truncated);
        }
        let aux = data.get(pos..pos.wrapping_add(l_aux)).ok_or(CsiError::Truncated)?.to_vec();
        pos = pos.wrapping_add(l_aux);

        // r[impl csi.n_ref]
        // r[impl csi.alloc_limits]
        let n_ref = read_count(data, &mut pos)?;
        check_limit(n_ref, MAX_INDEX_REFS, "n_ref")?;

        let pseudo_bin_id = csi_pseudo_bin(depth);
        let mut references = Vec::with_capacity(n_ref);

        for _ in 0..n_ref {
            let n_bin = read_count(data, &mut pos)?;
            check_limit(n_bin, MAX_BINS_PER_REF, "n_bin")?;
            let mut bins = Vec::with_capacity(n_bin);

            for _ in 0..n_bin {
                let bin_id = read_u32(data, &mut pos)?;
                // r[impl csi.no_linear_index]
                // r[impl csi.loffset]
                let loffset = VirtualOffset(read_u64(data, &mut pos)?);
                let n_chunk = read_count(data, &mut pos)?;
                check_limit(n_chunk, MAX_CHUNKS_PER_BIN, "n_chunk")?;
                let mut chunks = Vec::with_capacity(n_chunk);

                for _ in 0..n_chunk {
                    let begin = VirtualOffset(read_u64(data, &mut pos)?);
                    let end = VirtualOffset(read_u64(data, &mut pos)?);
                    chunks.push(Chunk { begin, end });
                }
                // r[impl csi.pseudo_bin]
                // Skip pseudo-bins during storage — they're metadata only
                if bin_id != pseudo_bin_id {
                    bins.push(CsiBin { bin_id, loffset, chunks });
                }
            }

            references.push(CsiRefIndex::new(bins));
        }

        Ok(CsiIndex { min_shift, depth, aux, references })
    }

    /// Binning parameters.
    pub fn min_shift(&self) -> u32 {
        self.min_shift
    }

    /// Binning depth.
    pub fn depth(&self) -> u32 {
        self.depth
    }

    /// Auxiliary data (tabix metadata if present).
    pub fn aux(&self) -> &[u8] {
        &self.aux
    }

    // r[impl csi.query]
    // r[impl csi.query_interface]
    /// Query the index for chunks overlapping a region [start, end] (0-based inclusive).
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
        let candidate_bins =
            csi_reg2bins(start_u64, end_u64.wrapping_add(1), self.min_shift, self.depth);
        let min_offset = csi_min_offset(ref_idx, start_u64, self.min_shift, self.depth);

        let mut result = Vec::new();
        for bin in &ref_idx.bins {
            if candidate_bins.contains(&bin.bin_id) {
                for chunk in &bin.chunks {
                    if chunk.end > min_offset {
                        result.push(AnnotatedChunk { chunk: *chunk, bin_id: bin.bin_id });
                    }
                }
            }
        }
        result.sort_by_key(|c| c.chunk.begin);
        result
    }

    // r[impl csi.query_split]
    /// Query the index, separating distant from nearby chunks.
    ///
    /// The split threshold is `depth - 3`: bins at levels 0 through `depth - 3`
    /// are "distant" (cached per-chromosome), levels `depth - 2` through `depth`
    /// are "nearby" (loaded per-query).
    pub fn query_split(&self, tid: u32, start: Pos0, end: Pos0) -> QueryChunks {
        let Some(ref_idx) = self.references.get(tid as usize) else {
            return QueryChunks { nearby: Vec::new(), distant: Vec::new() };
        };

        let start_u64 = start.as_u64();
        let end_u64 = end.as_u64();
        // reg2bins uses half-open interval
        let candidate_bins =
            csi_reg2bins(start_u64, end_u64.wrapping_add(1), self.min_shift, self.depth);

        // CSI uses per-bin loffset instead of a separate linear index
        let min_offset = csi_min_offset(ref_idx, start_u64, self.min_shift, self.depth);

        let cache_threshold = self.depth.saturating_sub(3);

        let mut nearby = Vec::new();
        let mut distant = Vec::new();
        for bin in &ref_idx.bins {
            if candidate_bins.contains(&bin.bin_id) {
                let level = csi_bin_level(bin.bin_id, self.depth);
                let target = if level <= cache_threshold { &mut distant } else { &mut nearby };
                for chunk in &bin.chunks {
                    if chunk.end > min_offset {
                        target.push(*chunk);
                    }
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

// r[impl csi.reg2bins]
/// Compute all CSI bins overlapping the interval [beg, end) (0-based, half-open),
/// parameterized by `min_shift` and `depth`.
#[expect(
    clippy::cast_possible_truncation,
    reason = "bin values are bounded by the CSI spec, fits in u32"
)]
fn csi_reg2bins(beg: u64, end: u64, min_shift: u32, depth: u32) -> Vec<u32> {
    let end = end.saturating_sub(1);
    let mut bins = Vec::with_capacity(32);
    let mut s = min_shift.wrapping_add(depth.wrapping_mul(3));
    let mut t = 0u64;

    for l in 0..=depth {
        let b = t.wrapping_add(beg >> s) as u32;
        let e = t.wrapping_add(end >> s) as u32;
        let mut i = b;
        while i <= e {
            bins.push(i);
            i = i.wrapping_add(1);
        }
        s = s.wrapping_sub(3);
        t = t.wrapping_add(1u64 << (l.wrapping_mul(3)));
    }

    bins
}

/// Compute the minimum virtual offset for a CSI query starting at `beg`.
///
/// Walks from the leaf-level bin containing `beg` up to the root, returning
/// the `loffset` of the first bin found. This replaces BAI's linear index.
fn csi_min_offset(ref_idx: &CsiRefIndex, beg: u64, min_shift: u32, depth: u32) -> VirtualOffset {
    // Leaf-level offset: ((1 << depth*3) - 1) / 7
    let leaf_offset = ((1u64 << depth.wrapping_mul(3)).wrapping_sub(1)).checked_div(7).unwrap_or(0);
    #[expect(
        clippy::cast_possible_truncation,
        reason = "bin values are bounded by the CSI spec, fits in u32"
    )]
    let mut bin = leaf_offset.wrapping_add(beg >> min_shift) as u32;

    // Walk from leaf to root looking for an existing bin with loffset.
    // Uses the parse-time HashMap to avoid an O(n_bins) scan per step.
    for _ in 0..=depth {
        if let Some(&i) = ref_idx.by_bin_id.get(&bin)
            && let Some(b) = ref_idx.bins.get(i)
        {
            return b.loffset;
        }
        // Parent bin: (bin - 1) / 8
        if bin == 0 {
            break;
        }
        bin = bin.wrapping_sub(1) / 8;
    }

    VirtualOffset(0)
}

/// Compute the level of a CSI bin given the binning depth.
fn csi_bin_level(bin_id: u32, depth: u32) -> u32 {
    let bin_id = u64::from(bin_id);
    let mut t = 0u64;
    for l in 0..=depth {
        let n_bins_at_level = 1u64 << (l.wrapping_mul(3));
        if bin_id < t.wrapping_add(n_bins_at_level) {
            return l;
        }
        t = t.wrapping_add(n_bins_at_level);
    }
    depth.wrapping_add(1) // pseudo-bin or invalid
}

// r[impl csi.pseudo_bin]
/// Pseudo-bin ID for the given depth (one past the max valid bin).
#[expect(clippy::cast_possible_truncation, reason = "pseudo-bin is bounded by depth, fits in u32")]
fn csi_pseudo_bin(depth: u32) -> u32 {
    let bits = depth.saturating_add(1).saturating_mul(3);
    let max_id = ((1u64 << bits).saturating_sub(1)).checked_div(7).unwrap_or(0);
    max_id.saturating_add(1) as u32
}

#[cfg(test)]
#[allow(clippy::arithmetic_side_effects, reason = "tests")]
#[allow(clippy::cast_possible_truncation, reason = "tests")]
#[allow(clippy::cast_possible_wrap, reason = "tests")]
#[allow(clippy::type_complexity, reason = "tests")]
mod tests {
    use super::*;
    use crate::vcf::index_builder::reg2bin;

    // --- csi_reg2bins tests ---

    #[test]
    fn csi_reg2bins_bai_compat_includes_bin0() {
        let bins = csi_reg2bins(0, 100, 14, 5);
        assert!(bins.contains(&0), "must include root bin");
    }

    #[test]
    fn csi_reg2bins_bai_compat_small_region() {
        let bins = csi_reg2bins(100, 200, 14, 5);
        // Should include leaf bin for positions 100-199
        let leaf_offset = ((1u64 << 15) - 1) / 7; // 4681
        let expected_leaf = leaf_offset as u32 + (100 >> 14) as u32;
        assert!(bins.contains(&expected_leaf));
    }

    #[test]
    fn csi_reg2bins_matches_bai_hardcoded() {
        // BAI's hardcoded reg2bins should match CSI's parameterized version
        // at min_shift=14, depth=5
        let bai_bins = {
            let mut bins = Vec::with_capacity(32);
            bins.push(0);
            let levels: [(u32, u32); 5] = [(1, 26), (9, 23), (73, 20), (585, 17), (4681, 14)];
            for &(offset, shift) in &levels {
                let mut k = offset + (100u64 >> shift) as u32;
                let end_k = offset + ((200u64 - 1) >> shift) as u32;
                while k <= end_k {
                    bins.push(k);
                    k += 1;
                }
            }
            bins
        };
        let csi_bins = csi_reg2bins(100, 200, 14, 5);
        assert_eq!(bai_bins, csi_bins, "CSI with BAI params must match BAI");
    }

    // --- csi_bin_level tests ---

    #[test]
    fn bin_level_bai_compat() {
        assert_eq!(csi_bin_level(0, 5), 0);
        assert_eq!(csi_bin_level(1, 5), 1);
        assert_eq!(csi_bin_level(8, 5), 1);
        assert_eq!(csi_bin_level(9, 5), 2);
        assert_eq!(csi_bin_level(72, 5), 2);
        assert_eq!(csi_bin_level(73, 5), 3);
        assert_eq!(csi_bin_level(584, 5), 3);
        assert_eq!(csi_bin_level(585, 5), 4);
        assert_eq!(csi_bin_level(4680, 5), 4);
        assert_eq!(csi_bin_level(4681, 5), 5);
        assert_eq!(csi_bin_level(37448, 5), 5); // last valid level-5 bin
        // bin 37449 is bin_limit (past all valid bins), 37450 is pseudo-bin
        assert_eq!(csi_bin_level(37449, 5), 6);
        assert_eq!(csi_bin_level(37450, 5), 6);
    }

    #[test]
    fn pseudo_bin_bai_compat() {
        assert_eq!(csi_pseudo_bin(5), 37450);
    }

    // --- reg2bin consistency ---

    #[test]
    fn csi_reg2bins_contains_reg2bin_result() {
        // The single bin from reg2bin must always be in the reg2bins list
        let cases = [(0u64, 100u64), (16384, 32768), (0, 512_000_000), (100_000, 200_000)];
        for (beg, end) in cases {
            let single = reg2bin(beg, end, 14, 5);
            let all = csi_reg2bins(beg, end, 14, 5);
            assert!(all.contains(&single), "reg2bin({beg},{end})={single} not in reg2bins={all:?}");
        }
    }

    // --- Proptest: CSI reg2bins matches BAI reg2bins at BAI-compatible params ---

    use proptest::prelude::*;

    proptest! {
        #[test]
        fn csi_reg2bins_matches_bai_for_all_positions(
            beg in 0u64..536_870_000u64,
            len in 1u64..100_000u64,
        ) {
            let end = beg.saturating_add(len).min(536_870_912);
            if beg >= end { return Ok(()); }

            let csi_bins = csi_reg2bins(beg, end, 14, 5);

            // Reproduce BAI's hardcoded algorithm
            let mut bai_bins = Vec::with_capacity(32);
            bai_bins.push(0);
            let levels: [(u32, u32); 5] = [(1, 26), (9, 23), (73, 20), (585, 17), (4681, 14)];
            for &(offset, shift) in &levels {
                let mut k = offset + (beg >> shift) as u32;
                let end_k = offset + ((end - 1) >> shift) as u32;
                while k <= end_k {
                    bai_bins.push(k);
                    k += 1;
                }
            }

            prop_assert_eq!(csi_bins, bai_bins, "CSI != BAI for [{}, {})", beg, end);
        }

        #[test]
        fn reg2bin_always_in_reg2bins(
            beg in 0u64..536_870_000u64,
            len in 1u64..100_000u64,
        ) {
            let end = beg.saturating_add(len).min(536_870_912);
            if beg >= end { return Ok(()); }

            let single = reg2bin(beg, end, 14, 5);
            let all = csi_reg2bins(beg, end, 14, 5);
            prop_assert!(
                all.contains(&single),
                "reg2bin({},{})={} not in reg2bins={:?}", beg, end, single, all
            );
        }
    }

    // --- CSI parsing tests ---

    /// Build a minimal valid CSI file in memory.
    fn build_csi_bytes(
        min_shift: i32,
        depth: i32,
        aux: &[u8],
        refs: &[Vec<(u32, u64, Vec<(u64, u64)>)>],
    ) -> Vec<u8> {
        let mut buf = Vec::new();
        buf.extend_from_slice(b"CSI\x01");
        buf.extend_from_slice(&min_shift.to_le_bytes());
        buf.extend_from_slice(&depth.to_le_bytes());
        buf.extend_from_slice(&(aux.len() as i32).to_le_bytes());
        buf.extend_from_slice(aux);
        buf.extend_from_slice(&(refs.len() as i32).to_le_bytes());
        for ref_bins in refs {
            buf.extend_from_slice(&(ref_bins.len() as i32).to_le_bytes());
            for (bin_id, loffset, chunks) in ref_bins {
                buf.extend_from_slice(&bin_id.to_le_bytes());
                buf.extend_from_slice(&loffset.to_le_bytes());
                buf.extend_from_slice(&(chunks.len() as i32).to_le_bytes());
                for (beg, end) in chunks {
                    buf.extend_from_slice(&beg.to_le_bytes());
                    buf.extend_from_slice(&end.to_le_bytes());
                }
            }
        }
        buf
    }

    #[test]
    fn parse_minimal_csi() {
        let data = build_csi_bytes(14, 5, &[], &[]);
        let idx = CsiIndex::from_bytes_inner(&data).unwrap();
        assert_eq!(idx.min_shift, 14);
        assert_eq!(idx.depth, 5);
        assert!(idx.aux.is_empty());
        assert!(idx.references.is_empty());
    }

    #[test]
    fn parse_csi_with_bins() {
        let data = build_csi_bytes(14, 5, &[], &[vec![(4681, 100, vec![(200, 300)])]]);
        let idx = CsiIndex::from_bytes_inner(&data).unwrap();
        assert_eq!(idx.references.len(), 1);
        assert_eq!(idx.references[0].bins.len(), 1);
        assert_eq!(idx.references[0].bins[0].bin_id, 4681);
        assert_eq!(idx.references[0].bins[0].loffset.0, 100);
        assert_eq!(idx.references[0].bins[0].chunks.len(), 1);
    }

    #[test]
    fn parse_csi_with_aux() {
        let aux = vec![1, 2, 3, 4, 5];
        let data = build_csi_bytes(14, 5, &aux, &[]);
        let idx = CsiIndex::from_bytes_inner(&data).unwrap();
        assert_eq!(idx.aux, aux);
    }

    #[test]
    fn parse_csi_skips_pseudo_bin() {
        // Pseudo-bin for depth=5 is 37450
        let data = build_csi_bytes(
            14,
            5,
            &[],
            &[vec![
                (4681, 100, vec![(200, 300)]),
                (37450, 0, vec![(0, 999), (10, 20)]), // pseudo-bin
            ]],
        );
        let idx = CsiIndex::from_bytes_inner(&data).unwrap();
        assert_eq!(idx.references[0].bins.len(), 1, "pseudo-bin should be filtered");
        assert_eq!(idx.references[0].bins[0].bin_id, 4681);
    }

    #[test]
    fn invalid_csi_magic() {
        let data = b"BAI\x01\x00\x00\x00\x00";
        let err = CsiIndex::from_bytes_inner(data).unwrap_err();
        assert!(matches!(err, CsiError::InvalidMagic));
    }

    #[test]
    fn truncated_csi() {
        let data = b"CSI\x01\x0e"; // too short for header
        let err = CsiIndex::from_bytes_inner(data).unwrap_err();
        assert!(matches!(err, CsiError::Truncated));
    }

    // r[verify csi.header_bounds]
    #[test]
    fn depth_22_panicked_before_fix() {
        // depth=22 causes 1u64 << 66 in csi_pseudo_bin — must return error, not panic
        let data = build_csi_bytes(14, 22, &[], &[]);
        let err = CsiIndex::from_bytes_inner(&data).unwrap_err();
        assert!(matches!(err, CsiError::DepthTooLarge { depth: 22, .. }));
    }

    // r[verify csi.header_bounds]
    #[test]
    fn depth_17_rejected() {
        // depth=17 with min_shift=14 → shift = 14+51 = 65, overflows u64 shift
        let data = build_csi_bytes(14, 17, &[], &[]);
        let err = CsiIndex::from_bytes_inner(&data).unwrap_err();
        assert!(matches!(err, CsiError::DepthTooLarge { .. }));
    }

    // r[verify csi.header_bounds]
    #[test]
    fn min_shift_upper_bound_rejected() {
        // min_shift + depth*3 must remain < 64 to avoid Rust's shift-amount
        // masking on u64 (e.g. `beg >> 64` evaluates as `beg >> 0`).
        // min_shift=64 with depth=0 already overflows.
        let data = build_csi_bytes(64, 0, &[], &[]);
        let err = CsiIndex::from_bytes_inner(&data).unwrap_err();
        assert!(matches!(err, CsiError::MinShiftTooLarge { .. }), "got {err:?}");
    }

    // r[verify csi.header_bounds]
    #[test]
    fn min_shift_combined_with_depth_rejected() {
        // min_shift=50, depth=5 → s = 50 + 15 = 65, overflows u64 shift
        let data = build_csi_bytes(50, 5, &[], &[]);
        let err = CsiIndex::from_bytes_inner(&data).unwrap_err();
        assert!(matches!(err, CsiError::MinShiftTooLarge { .. }), "got {err:?}");
    }

    // r[verify csi.header_bounds]
    #[test]
    fn min_shift_at_boundary_accepted() {
        // min_shift=15, depth=16 → s = 15 + 48 = 63, just under 64. Accepted.
        let data = build_csi_bytes(15, 16, &[], &[]);
        let idx = CsiIndex::from_bytes_inner(&data).unwrap();
        assert_eq!(idx.min_shift(), 15);
        assert_eq!(idx.depth(), 16);
    }

    // r[verify csi.header_bounds]
    #[test]
    fn depth_16_accepted() {
        // depth=16 is the maximum allowed (min_shift+depth*3 = 14+48 = 62, safe)
        let data = build_csi_bytes(14, 16, &[], &[]);
        let idx = CsiIndex::from_bytes_inner(&data).unwrap();
        assert_eq!(idx.depth(), 16);
    }

    // r[verify csi.alloc_limits]
    #[test]
    fn negative_n_ref_rejected() {
        // n_ref = -5 should produce NegativeCount, not a giant allocation
        let mut data = Vec::new();
        data.extend_from_slice(b"CSI\x01");
        data.extend_from_slice(&14i32.to_le_bytes()); // min_shift
        data.extend_from_slice(&5i32.to_le_bytes()); // depth
        data.extend_from_slice(&0i32.to_le_bytes()); // l_aux
        data.extend_from_slice(&(-5i32).to_le_bytes()); // n_ref = -5
        let err = CsiIndex::from_bytes_inner(&data).unwrap_err();
        assert!(matches!(err, CsiError::NegativeCount { value: -5 }), "got: {err:?}");
    }

    // r[verify csi.alloc_limits]
    #[test]
    fn negative_n_bin_rejected() {
        // First ref has n_bin = -1
        let mut data = Vec::new();
        data.extend_from_slice(b"CSI\x01");
        data.extend_from_slice(&14i32.to_le_bytes());
        data.extend_from_slice(&5i32.to_le_bytes());
        data.extend_from_slice(&0i32.to_le_bytes());
        data.extend_from_slice(&1i32.to_le_bytes()); // n_ref = 1
        data.extend_from_slice(&(-1i32).to_le_bytes()); // n_bin = -1
        let err = CsiIndex::from_bytes_inner(&data).unwrap_err();
        assert!(matches!(err, CsiError::NegativeCount { value: -1 }), "got: {err:?}");
    }

    // --- Query tests ---

    fn make_csi_index(bins: Vec<(u32, u64, Vec<(u64, u64)>)>) -> CsiIndex {
        let bins: Vec<CsiBin> = bins
            .into_iter()
            .map(|(bin_id, loffset, chunks)| CsiBin {
                bin_id,
                loffset: VirtualOffset(loffset),
                chunks: chunks
                    .into_iter()
                    .map(|(b, e)| Chunk { begin: VirtualOffset(b), end: VirtualOffset(e) })
                    .collect(),
            })
            .collect();
        CsiIndex {
            min_shift: 14,
            depth: 5,
            aux: Vec::new(),
            references: vec![CsiRefIndex::new(bins)],
        }
    }

    #[test]
    fn query_returns_matching_chunks() {
        let idx = make_csi_index(vec![(4681, 0, vec![(100, 500)])]);
        let chunks = idx.query(0, Pos0::new(0).unwrap(), Pos0::new(100).unwrap());
        assert_eq!(chunks.len(), 1);
        assert_eq!(chunks[0].begin.0, 100);
    }

    #[test]
    fn query_split_separates_by_level() {
        let idx = make_csi_index(vec![(0, 0, vec![(9000, 9500)]), (4681, 0, vec![(100, 500)])]);
        let result = idx.query_split(0, Pos0::new(0).unwrap(), Pos0::new(100).unwrap());
        assert_eq!(result.distant.len(), 1);
        assert_eq!(result.distant[0].begin.0, 9000);
        assert!(!result.nearby.is_empty());
    }

    #[test]
    fn query_empty_tid() {
        let idx = CsiIndex { min_shift: 14, depth: 5, aux: Vec::new(), references: Vec::new() };
        let chunks = idx.query(0, Pos0::new(0).unwrap(), Pos0::new(100).unwrap());
        assert!(chunks.is_empty());
    }

    /// `csi_min_offset` walks from leaf to root. With a dense bin set spanning
    /// many sibling bins the lookup must still return the correct `loffset`
    /// for the queried position. This exercises the bin-id index used in
    /// place of the prior linear scan.
    #[test]
    fn csi_min_offset_dense_bins() {
        // Build a ref with many leaf-level bins, each with a distinct loffset.
        // depth=5, min_shift=14 ⇒ leaf offset = 4681. Insert leaf bins
        // 4681..4781 with loffset = (bin_id - 4681) * 1000.
        let bins: Vec<(u32, u64, Vec<(u64, u64)>)> = (0..100u32)
            .map(|i| {
                let bin_id = 4681 + i;
                let loff = u64::from(i) * 1000;
                (bin_id, loff, vec![(loff, loff + 500)])
            })
            .collect();
        let idx = make_csi_index(bins);
        let ref_idx = &idx.references[0];

        // The leaf bin for position p (with min_shift=14) is 4681 + (p >> 14).
        for i in 0..100u32 {
            let pos = u64::from(i) << 14; // start of the i-th 16 KiB window
            let off = csi_min_offset(ref_idx, pos, 14, 5);
            assert_eq!(off.0, u64::from(i) * 1000, "dense-bin lookup wrong at i={i}, pos={pos}");
        }
    }

    #[test]
    fn loffset_filters_chunks() {
        // bin 4681 has loffset=400, so chunks ending before 400 should be skipped
        let idx = make_csi_index(vec![(0, 0, vec![(50, 200)]), (4681, 400, vec![(400, 800)])]);
        let result = idx.query_split(0, Pos0::new(0).unwrap(), Pos0::new(100).unwrap());
        // csi_min_offset returns 400 from the leaf bin, so the level-0 chunk
        // (50, 200) with end=200 < 400 is filtered out
        assert!(
            result.distant.is_empty() || result.distant.iter().all(|c| c.end.0 >= 400),
            "chunks ending before loffset should be filtered"
        );
    }

    #[test]
    fn query_split_union_equals_query() {
        let idx = make_csi_index(vec![
            (0, 0, vec![(9000, 9500)]),
            (4681, 0, vec![(100, 500)]),
            (4682, 0, vec![(500, 900)]),
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

    // --- Non-default min_shift/depth tests ---

    #[test]
    fn csi_reg2bins_different_params() {
        // With min_shift=15, depth=6: larger coordinate space
        let bins = csi_reg2bins(0, 100, 15, 6);
        assert!(bins.contains(&0), "must include root bin");
        // Leaf-level offset for depth=6: ((1<<18)-1)/7 = 37449
        let leaf_offset = ((1u64 << 18) - 1) / 7;
        let expected_leaf = leaf_offset as u32;
        assert!(bins.contains(&expected_leaf), "must include leaf bin");
    }

    #[test]
    fn bin_level_different_depth() {
        // depth=3: levels 0-3
        assert_eq!(csi_bin_level(0, 3), 0);
        assert_eq!(csi_bin_level(1, 3), 1);
        assert_eq!(csi_bin_level(9, 3), 2);
        assert_eq!(csi_bin_level(73, 3), 3);
        assert_eq!(csi_bin_level(584, 3), 3);
    }
}
