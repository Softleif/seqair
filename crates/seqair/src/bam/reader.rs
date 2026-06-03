//! Open and query BAM files. [`IndexedBamReader`] parses the BAI index, then fetches records
//! for a region into a [`RecordStore`]. Call [`IndexedBamReader::fork`] to get a cheap
//! per-thread reader that shares the parsed index and header via [`Arc<BamShared>`].

use crate::reader::FetchCounts;

use super::{
    bgzf::{BgzfError, BgzfReader, VirtualOffset},
    csi_index::CsiIndex,
    flags::BamFlags,
    header::{BamHeader, BamHeaderError},
    index::{AlignmentIndex, BaiError, BamIndex, Chunk},
    record::{DecodeError, compute_end_pos_from_raw},
    record_store::{CustomizeRecordStore, RecordStore},
    region_buf::{self, RegionBuf},
};
use seqair_types::{Pos0, SmolStr};
use std::{
    fs::File,
    io::{Read, Seek},
    path::{Path, PathBuf},
    sync::Arc,
};
use tracing::instrument;

// r[impl io.errors]
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum BamError {
    #[error("I/O error opening {path}")]
    Open { path: PathBuf, source: std::io::Error },

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

    #[error("BAM index error")]
    Index {
        #[from]
        source: BaiError,
    },

    #[error("CSI index error")]
    CsiIndex {
        #[from]
        source: super::csi_index::CsiError,
    },

    #[error("truncated BAM record at virtual offset {offset:#x}")]
    TruncatedRecord { offset: u64 },

    #[error("contig not found: {name}")]
    ContigNotFound { name: SmolStr },

    #[error("region {contig}:{start}-{end} is out of bounds (contig length: {contig_len})")]
    RegionOutOfBounds { contig: SmolStr, start: u64, end: u64, contig_len: u64 },

    #[error("index not found for {bam_path} (checked .csi and .bai)")]
    IndexNotFound { bam_path: PathBuf },

    #[error("record decode error")]
    RecordDecode {
        #[from]
        source: DecodeError,
    },

    // r[impl bam.reader.coordinate_overflow]
    #[error("coordinate overflow: tid value {value} exceeds {max}")]
    TidOverflow { value: u64, max: u64 },

    #[error("invalid BAM position value {value}: negative positions are reserved")]
    InvalidPosition { value: i32 },
}

// r[impl bam.reader.coordinate_overflow]
fn validate_tid(tid: u32) -> Result<i32, BamError> {
    i32::try_from(tid)
        .map_err(|_| BamError::TidOverflow { value: u64::from(tid), max: i32::MAX as u64 })
}

// r[impl bam.reader.shared_state]
pub struct BamShared {
    index: AlignmentIndex,
    header: BamHeader,
    pub bam_path: PathBuf,
}

impl BamShared {
    pub fn index(&self) -> &AlignmentIndex {
        &self.index
    }
}

// r[impl bam.reader.open]
// r[impl bam.reader.header_access]
pub struct IndexedBamReader<R: Read + Seek = File> {
    /// Separate reader handle for bulk region reads (unbuffered — `RegionBuf` does
    /// large sequential reads that don't benefit from `BufReader`).
    bulk_reader: R,
    shared: Arc<BamShared>,
}

impl<R: Read + Seek> std::fmt::Debug for IndexedBamReader<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("IndexedBamReader").field("bam_path", &self.shared.bam_path).finish()
    }
}

impl IndexedBamReader<File> {
    #[instrument(level = "debug", fields(path = %path.display()))]
    pub fn open(path: &Path) -> Result<Self, BamError> {
        let mut bgzf = BgzfReader::open(path)?;
        let header = BamHeader::parse(&mut bgzf)?;
        // r[impl unified.sort_order]
        header.validate_sort_order()?;
        let index = find_and_open_index(path)?;

        let bulk_file = File::open(path)
            .map_err(|source| BamError::Open { path: path.to_path_buf(), source })?;

        Ok(IndexedBamReader {
            bulk_reader: bulk_file,
            shared: Arc::new(BamShared { index, header, bam_path: path.to_path_buf() }),
        })
    }

    // r[impl bam.reader.fork]
    // r[impl bam.reader.fork_independence]
    // r[impl bam.reader.fork_equivalence]
    // r[impl bam.reader.fork_concurrent]
    #[instrument(level = "debug", skip(self), fields(path = %self.shared.bam_path.display()))]
    pub fn fork(&self) -> Result<Self, BamError> {
        let bulk_file = File::open(&self.shared.bam_path)
            .map_err(|source| BamError::Open { path: self.shared.bam_path.clone(), source })?;

        Ok(IndexedBamReader { bulk_reader: bulk_file, shared: Arc::clone(&self.shared) })
    }
}

#[cfg(feature = "fuzz")]
impl IndexedBamReader<std::io::Cursor<Vec<u8>>> {
    pub fn from_bytes(bam_data: Vec<u8>, bai_data: &[u8]) -> Result<Self, BamError> {
        let mut bgzf = BgzfReader::from_cursor(bam_data.clone());
        let header = BamHeader::parse(&mut bgzf)?;
        header.validate_sort_order()?;
        let index = AlignmentIndex::Bai(BamIndex::from_bytes(bai_data)?);

        Ok(IndexedBamReader {
            bulk_reader: std::io::Cursor::new(bam_data),
            shared: Arc::new(BamShared { index, header, bam_path: PathBuf::from("<fuzz>") }),
        })
    }
}

impl<R: Read + Seek> IndexedBamReader<R> {
    // r[impl bam.reader.fork_arc_identity]
    pub fn shared(&self) -> &Arc<BamShared> {
        &self.shared
    }

    pub fn header(&self) -> &BamHeader {
        &self.shared.header
    }

    // r[impl bam.reader.fetch_into+2]
    // r[impl bam.reader.overlap_filter]
    // r[impl bam.reader.sorted_order+2]
    // r[impl bam.reader.secondary_supplementary_included+2]
    // r[impl region_buf.fetch_into+2]
    // r[impl region_buf.no_bin0]
    #[instrument(level = "debug", skip(self, store), fields(tid, start, end))]
    pub fn fetch_into(
        &mut self,
        tid: u32,
        start: Pos0,
        end: Pos0,
        store: &mut RecordStore,
    ) -> Result<usize, BamError> {
        let mut query = self.query(tid, start, end)?;
        store.clear();

        let mut kept = 0usize;
        query.for_each_result::<_, BamError>(|raw| {
            if store.push_raw(raw, &mut ())?.is_some() {
                kept = kept.saturating_add(1);
            }
            Ok(())
        })?;

        Ok(kept)
    }

    // r[impl unified.fetch_into_customized]
    /// Customized variant: each record that passes the reader's built-in
    /// checks is pushed via [`RecordStore::push_raw`], which forwards the
    /// `customize` value's `filter` to decide retention. Rejected
    /// records roll back their slab writes. The returned [`FetchCounts`]
    /// reports `fetched` (produced by the reader) vs `kept` (survived the
    /// filter).
    pub fn fetch_into_customized<E: CustomizeRecordStore>(
        &mut self,
        tid: u32,
        start: Pos0,
        end: Pos0,
        store: &mut RecordStore<E::Extra>,
        customize: &mut E,
    ) -> Result<FetchCounts, BamError> {
        store.clear();

        let mut query = self.query(tid, start, end)?;
        let mut kept_count: usize = 0;
        query.for_each_result::<_, BamError>(|raw| {
            if store.push_raw(raw, customize)?.is_some() {
                kept_count = kept_count.saturating_add(1);
            }
            Ok(())
        })?;

        let qc = query.counts();
        Ok(FetchCounts { fetched: qc.fetched, kept: kept_count })
    }

    /// Create a cursor over raw BAM records in a genomic region.
    ///
    /// Returns a [`BamQuery`] that yields raw `&[u8]` record bytes (including the
    /// 32-byte BAM fixed header) for records overlapping `[start, end)` on `tid`.
    /// Records are pre-filtered by tid and position range — the caller sees only
    /// records that fall within the requested region.
    ///
    /// Unlike [`fetch_into_customized`](Self::fetch_into_customized), this
    /// does not involve [`RecordStore`] or [`CustomizeRecordStore`]. The caller
    /// receives raw bytes and can parse, filter, or store them however they want.
    ///
    /// # Example
    ///
    /// ```ignore
    /// let mut query = reader.query(tid, start, end)?;
    /// query.for_each(|raw| {
    ///     let ref_id = i32::from_le_bytes(raw[0..4].try_into().unwrap());
    ///     // ... process raw bytes ...
    /// })?;
    /// let counts = query.counts();
    /// ```
    pub fn query(&mut self, tid: u32, start: Pos0, end: Pos0) -> Result<BamQuery<'_, R>, BamError> {
        let chunks = self.shared.index.query(tid, start, end);
        let tid_i32 = validate_tid(tid)?;
        let batches = if chunks.is_empty() {
            Vec::new()
        } else {
            partition_chunks(&chunks, region_buf::MAX_REGION_BYTES)
        };

        Ok(BamQuery {
            reader: &mut self.bulk_reader,
            batches,
            batch_idx: 0,
            region: None,
            chunk_idx: 0,
            chunk_end: VirtualOffset(0),
            scratch: Vec::new(),
            tid_i32,
            start,
            end,
            accepted: 0,
            skipped_tid: 0,
            skipped_out_of_range: 0,
        })
    }
}

// ── BamQuery ────────────────────────────────────────────────────────────

/// Streaming cursor over raw BAM records within a genomic region.
///
/// Created by [`IndexedBamReader::query`]. Call [`for_each`](Self::for_each) to
/// process pre-filtered (tid + position) raw record bytes with zero-copy
/// closure-based iteration. Use [`for_each_result`](Self::for_each_result) if
/// the closure can fail.
///
/// The cursor borrows the reader's file handle exclusively. Drop it to
/// release the borrow.
pub struct BamQuery<'r, R: Read + Seek> {
    reader: &'r mut R,
    batches: Vec<Vec<Chunk>>,
    batch_idx: usize,
    region: Option<RegionBuf>,
    chunk_idx: usize,
    chunk_end: VirtualOffset,
    scratch: Vec<u8>,
    tid_i32: i32,
    start: Pos0,
    end: Pos0,
    accepted: u32,
    skipped_tid: u32,
    skipped_out_of_range: u32,
}

impl<R: Read + Seek> std::fmt::Debug for BamQuery<'_, R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BamQuery")
            .field("batches", &self.batches.len())
            .field("batch_idx", &self.batch_idx)
            .field("chunk_idx", &self.chunk_idx)
            .field("chunk_end", &self.chunk_end)
            .field("accepted", &self.accepted)
            .field("skipped_tid", &self.skipped_tid)
            .field("skipped_out_of_range", &self.skipped_out_of_range)
            .finish()
    }
}

/// Snapshot of the cursor's internal counters.
///
/// Returned by [`BamQuery::counts`]. Unlike [`FetchCounts`] there is no
/// `kept` field — the cursor does not filter beyond position and tid.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BamQueryCounts {
    /// Records that passed the built-in tid and position checks.
    pub fetched: usize,
    /// Records skipped because their reference ID didn't match.
    pub skipped_tid: usize,
    /// Records skipped because they fell outside the requested range.
    pub skipped_out_of_range: usize,
}

impl<'r, R: Read + Seek> BamQuery<'r, R> {
    /// Iterate over raw BAM records in the region, calling `f` for each.
    ///
    /// Zero-copy: the `&[u8]` passed to `f` borrows from the cursor's
    /// internal decompressed buffer (or scratch buffer for cross-block records)
    /// and is valid only within the closure call.
    ///
    /// Records are already filtered by `tid` and genomic position `[start, end)`.
    /// I/O and decode errors are returned immediately (the iteration stops).
    ///
    /// Returns the final [`BamQueryCounts`] on success.
    pub fn for_each<F>(&mut self, mut f: F) -> Result<BamQueryCounts, BamError>
    where
        F: FnMut(&[u8]),
    {
        self.run_loop(|raw| {
            f(raw);
            Ok(())
        })
    }

    /// Like [`for_each`](Self::for_each), but the closure can return an error
    /// that propagates out of the iteration.
    ///
    /// This is used internally by [`fetch_into_customized`] and is also useful
    /// for callers that want to short-circuit on their own error conditions.
    pub fn for_each_result<F, E>(&mut self, f: F) -> Result<BamQueryCounts, E>
    where
        F: FnMut(&[u8]) -> Result<(), E>,
        E: From<BamError>,
    {
        self.run_loop(f)
    }

    fn run_loop<F, E>(&mut self, mut f: F) -> Result<BamQueryCounts, E>
    where
        F: FnMut(&[u8]) -> Result<(), E>,
        E: From<BamError>,
    {
        loop {
            if self.region.is_none() {
                if self.batch_idx >= self.batches.len() {
                    break;
                }
                debug_assert!(self.batch_idx < self.batches.len());
                #[allow(clippy::indexing_slicing, reason = "bounds checked above")]
                let batch = &self.batches[self.batch_idx];
                self.batch_idx = self.batch_idx.saturating_add(1);
                debug_assert!(!batch.is_empty());

                let mut region =
                    RegionBuf::load(self.reader, batch).map_err(|e| E::from(BamError::from(e)))?;
                #[allow(clippy::indexing_slicing, reason = "non-empty assert above")]
                let first_chunk = batch[0];
                region.seek_virtual(first_chunk.begin).map_err(|e| E::from(BamError::from(e)))?;
                self.chunk_end = first_chunk.end;
                self.chunk_idx = 0;
                self.region = Some(region);
            }

            let region = self.region.as_mut().expect("region guaranteed Some above");

            if region.virtual_offset() >= self.chunk_end {
                self.chunk_idx = self.chunk_idx.saturating_add(1);
                debug_assert!(self.batch_idx >= 1, "batch_idx incremented on load");
                #[allow(
                    clippy::indexing_slicing,
                    reason = "batch_idx >= 1 after first load; saturating_sub(1) safe"
                )]
                let batch = &self.batches[self.batch_idx.saturating_sub(1)];
                if self.chunk_idx >= batch.len() {
                    self.region = None;
                    continue;
                }
                #[allow(clippy::indexing_slicing, reason = "chunk_idx < batch.len() checked above")]
                let chunk = &batch[self.chunk_idx];
                region.seek_virtual(chunk.begin).map_err(|e| E::from(BamError::from(e)))?;
                self.chunk_end = chunk.end;
                continue;
            }

            // r[impl bam.reader.propagate_errors]
            let current_voff = region.virtual_offset();
            let raw = match region.read_record(&mut self.scratch) {
                Ok(s) => s,
                Err(BgzfError::UnexpectedEof) => continue,
                Err(e) => return Err(E::from(BamError::from(e))),
            };

            if raw.len() < 32 {
                return Err(E::from(BamError::TruncatedRecord { offset: current_voff.0 }));
            }

            // r[impl bam.reader.unmapped_skipped]
            debug_assert!(raw.len() >= 32, "raw record too short: {}", raw.len());
            #[allow(clippy::indexing_slicing, reason = "raw.len() >= 32 checked above")]
            let rec_tid = i32::from_le_bytes([raw[0], raw[1], raw[2], raw[3]]);
            #[allow(clippy::indexing_slicing, reason = "raw.len() >= 32 checked above")]
            let rec_pos_raw = i32::from_le_bytes([raw[4], raw[5], raw[6], raw[7]]);
            let rec_pos = Pos0::try_from(rec_pos_raw)
                .map_err(|_| E::from(BamError::InvalidPosition { value: rec_pos_raw }))?;
            #[allow(clippy::indexing_slicing, reason = "raw.len() >= 32 checked above")]
            let rec_flags = BamFlags::from(u16::from_le_bytes([raw[14], raw[15]]));

            if rec_tid != self.tid_i32 {
                self.skipped_tid = self.skipped_tid.saturating_add(1);
                continue;
            }

            // r[impl record_store.end_pos_htslib]
            let rec_end = if rec_flags.is_unmapped() {
                rec_pos
            } else {
                compute_end_pos_from_raw(raw).unwrap_or(rec_pos)
            };
            if rec_pos > self.end || rec_end < self.start {
                self.skipped_out_of_range = self.skipped_out_of_range.saturating_add(1);
                continue;
            }

            self.accepted = self.accepted.saturating_add(1);
            f(raw)?;
        }

        Ok(self.counts())
    }

    /// Return the current counter snapshot.
    ///
    /// Can be called at any point during iteration — mid-stream, after
    /// exhaustion, or after an error. Useful for diagnostics.
    pub fn counts(&self) -> BamQueryCounts {
        BamQueryCounts {
            fetched: self.accepted as usize,
            skipped_tid: self.skipped_tid as usize,
            skipped_out_of_range: self.skipped_out_of_range as usize,
        }
    }
}

/// Partition chunks into batches where each batch's merged compressed byte
/// range fits within `max_bytes`.
///
/// Chunks are added greedily in order. When adding the next chunk would push
/// the batch over the limit, a new batch is started. A single chunk that
/// exceeds `max_bytes` on its own gets its own batch — `RegionBuf::load`
/// handles oversized ranges by warning and allocating the needed memory.
fn partition_chunks(chunks: &[Chunk], max_bytes: usize) -> Vec<Vec<Chunk>> {
    if chunks.is_empty() {
        return Vec::new();
    }

    // Fast path: if everything fits, return a single batch (avoids recomputing).
    if region_buf::merged_byte_size(chunks) <= max_bytes {
        return vec![chunks.to_vec()];
    }

    let mut batches: Vec<Vec<Chunk>> = Vec::new();
    let mut current_batch: Vec<Chunk> = Vec::new();

    for chunk in chunks {
        if current_batch.is_empty() {
            current_batch.push(*chunk);
            continue;
        }

        // Tentatively add this chunk and check if we still fit.
        current_batch.push(*chunk);
        if region_buf::merged_byte_size(&current_batch) <= max_bytes {
            continue;
        }

        // Doesn't fit — remove it and start a new batch.
        current_batch.pop();
        batches.push(std::mem::take(&mut current_batch));
        current_batch.push(*chunk);
    }

    if !current_batch.is_empty() {
        batches.push(current_batch);
    }

    tracing::info!(
        batches = batches.len(),
        total_chunks = chunks.len(),
        "region split into multiple batches due to size"
    );

    batches
}

// r[impl unified.detect_index]
// r[impl unified.detect_error]
// r[impl csi.detect]
fn find_and_open_index(bam_path: &Path) -> Result<AlignmentIndex, BamError> {
    // CSI preferred: try {file}.csi first, then {file_without_ext}.csi (per htslib convention)
    let csi_path = bam_path.with_extension("bam.csi");
    if csi_path.exists() {
        return Ok(AlignmentIndex::Csi(CsiIndex::from_path(&csi_path)?));
    }

    let csi_path2 = bam_path.with_extension("csi");
    if csi_path2.exists() {
        return Ok(AlignmentIndex::Csi(CsiIndex::from_path(&csi_path2)?));
    }

    // Fall back to BAI: try {file}.bai first, then {file_without_ext}.bai
    let bai_path = bam_path.with_extension("bam.bai");
    if bai_path.exists() {
        return Ok(AlignmentIndex::Bai(
            BamIndex::from_path(&bai_path).map_err(|source| BamError::Index { source })?,
        ));
    }

    let bai_path2 = bam_path.with_extension("bai");
    if bai_path2.exists() {
        return Ok(AlignmentIndex::Bai(
            BamIndex::from_path(&bai_path2).map_err(|source| BamError::Index { source })?,
        ));
    }

    Err(BamError::IndexNotFound { bam_path: bam_path.to_path_buf() })
}

#[cfg(test)]
#[allow(clippy::arithmetic_side_effects, reason = "test code with controlled values")]
mod tests {
    use super::*;

    // r[verify bam.reader.coordinate_overflow]
    #[test]
    fn fetch_into_rejects_tid_exceeding_i32_max() {
        // tid that exceeds i32::MAX should return TidOverflow
        let tid: u32 = u32::try_from(i32::MAX as u64 + 1).unwrap(); // 2_147_483_648
        let result = validate_tid(tid);
        assert!(result.is_err(), "tid > i32::MAX must error");
        let err = result.unwrap_err();
        assert!(
            matches!(err, BamError::TidOverflow { .. }),
            "expected TidOverflow for tid, got: {err}"
        );
    }

    // r[verify bam.reader.coordinate_overflow]
    #[test]
    fn fetch_into_accepts_max_valid_tid() {
        let result = validate_tid(i32::MAX as u32);
        assert!(result.is_ok(), "max valid tid must succeed");
    }

    use super::super::bgzf::VirtualOffset;

    #[test]
    fn partition_chunks_empty() {
        let result = partition_chunks(&[], 1024);
        assert!(result.is_empty());
    }

    #[test]
    fn partition_chunks_single_batch_when_small() {
        let chunks = vec![
            Chunk { begin: VirtualOffset::new(100, 0), end: VirtualOffset::new(200, 0) },
            Chunk { begin: VirtualOffset::new(300, 0), end: VirtualOffset::new(400, 0) },
        ];
        let result = partition_chunks(&chunks, region_buf::MAX_REGION_BYTES);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].len(), 2);
    }

    #[test]
    fn partition_chunks_splits_large_ranges() {
        // Create two disjoint chunks that individually fit but together exceed max_bytes.
        // Chunks must be spaced > CHUNK_END_PAD apart to stay disjoint after padding.
        let max_bytes = 10_000_000; // 10 MiB
        let gap = (region_buf::CHUNK_END_PAD as u64) + 1_000_000; // > CHUNK_END_PAD
        let span = 5_000_000u64; // each chunk spans 5 MiB

        let chunks = vec![
            Chunk { begin: VirtualOffset::new(0, 0), end: VirtualOffset::new(span, 0) },
            Chunk {
                begin: VirtualOffset::new(span + gap, 0),
                end: VirtualOffset::new(span + gap + span, 0),
            },
        ];

        // Together these exceed max_bytes (two disjoint ranges of ~7 MiB each)
        let total = region_buf::merged_byte_size(&chunks);
        assert!(total > max_bytes, "test setup: total {total} must exceed {max_bytes}");

        let result = partition_chunks(&chunks, max_bytes);
        assert_eq!(result.len(), 2, "should split into 2 batches");

        // Each batch individually fits
        for batch in &result {
            let size = region_buf::merged_byte_size(batch);
            assert!(size <= max_bytes, "batch size {size} exceeds {max_bytes}");
        }
    }

    // r[verify bam.reader.chunk_batching]
    /// Reproduces the 122 GB BAM scenario: many chunks that merge to >256 MiB.
    /// Verifies that `partition_chunks` splits them into batches that each fit.
    #[test]
    fn partition_chunks_122gb_bam_scenario() {
        // Simulate chunk layout from the real BAM:
        // - One huge chunk spanning ~131 MiB (bins 13421-13422)
        // - Many smaller chunks from higher-level bins that overlap and extend the range
        let mut chunks = vec![
            // Large level-5 bin chunk: 131 MiB compressed span
            Chunk {
                begin: VirtualOffset::new(5_971_912_384, 12584),
                end: VirtualOffset::new(6_141_438_390, 36748),
            },
        ];

        // Add many smaller chunks from higher-level bins scattered in and around
        for offset in (6_141_500_000..6_252_000_000u64).step_by(15_000) {
            chunks.push(Chunk {
                begin: VirtualOffset::new(offset, 0),
                end: VirtualOffset::new(offset, 50000),
            });
        }
        // Sort by begin (as query() does)
        chunks.sort_by_key(|c| c.begin);

        let total = region_buf::merged_byte_size(&chunks);
        assert!(
            total > region_buf::MAX_REGION_BYTES,
            "test setup: total {total} must exceed MAX_REGION_BYTES"
        );

        let batches = partition_chunks(&chunks, region_buf::MAX_REGION_BYTES);
        assert!(batches.len() >= 2, "should need at least 2 batches, got {}", batches.len());

        // Every batch must fit within the limit
        for (i, batch) in batches.iter().enumerate() {
            let size = region_buf::merged_byte_size(batch);
            assert!(
                size <= region_buf::MAX_REGION_BYTES,
                "batch {i} has size {size} exceeding limit {}",
                region_buf::MAX_REGION_BYTES
            );
        }

        // All chunks are preserved (no loss, no duplication)
        let total_chunks: usize = batches.iter().map(|b| b.len()).sum();
        assert_eq!(total_chunks, chunks.len());
    }

    proptest::proptest! {
        /// Any set of chunks must be partitioned such that every batch fits
        /// within the limit, and no chunks are lost or duplicated.
        #[test]
        fn proptest_partition_preserves_all_chunks(
            n_chunks in 1usize..20,
            seed in 0u64..10_000,
        ) {
            let max_bytes = 500_000;
            // Space chunks far enough apart to stay disjoint after CHUNK_END_PAD.
            let spacing = 200_000u64;
            let chunks: Vec<Chunk> = (0..n_chunks)
                .map(|i| {
                    let base = seed + (i as u64) * spacing;
                    Chunk {
                        begin: VirtualOffset::new(base, 0),
                        end: VirtualOffset::new(base + 100_000, 0),
                    }
                })
                .collect();

            let batches = partition_chunks(&chunks, max_bytes);

            // All chunks preserved (no loss, no duplication from splitting)
            let total: usize = batches.iter().map(|b| b.len()).sum();
            proptest::prop_assert_eq!(total, chunks.len());

            // Each batch fits (or is a single oversized chunk)
            for batch in &batches {
                if batch.len() > 1 {
                    let size = region_buf::merged_byte_size(batch);
                    proptest::prop_assert!(
                        size <= max_bytes,
                        "batch size {size} > {max_bytes}"
                    );
                }
            }
        }
    }
}
