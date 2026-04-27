//! Open and query BAM files. [`IndexedBamReader`] parses the BAI index, then fetches records
//! for a region into a [`RecordStore`]. Call [`IndexedBamReader::fork`] to get a cheap
//! per-thread reader that shares the parsed index and header via [`Arc<BamShared>`].

use crate::reader::FetchCounts;

use super::{
    bgzf::{BgzfError, BgzfReader},
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
        self.fetch_into_customized(tid, start, end, store, &mut ()).map(|c| c.kept)
    }

    // r[impl unified.fetch_into_customized]
    /// Customized variant: each record that passes the reader's built-in
    /// checks is pushed via [`RecordStore::push_raw`], which forwards the
    /// `customize` value's `keep_record` to decide retention. Rejected
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

        let chunks = self.shared.index.query(tid, start, end);
        if chunks.is_empty() {
            return Ok(FetchCounts::default());
        }

        let tid_i32 = validate_tid(tid)?;

        let mut skipped_tid: u32 = 0;
        let mut skipped_unmapped: u32 = 0;
        let mut skipped_out_of_range: u32 = 0;
        let mut accepted: u32 = 0;
        let mut kept_count: usize = 0;

        // Scratch buffer for the rare case where a record body straddles a BGZF
        // block boundary; zero-copy slice from RegionBuf::buf is used otherwise.
        let mut scratch: Vec<u8> = Vec::new();

        // r[impl bam.reader.chunk_batching]
        // Partition chunks into batches that each fit within MAX_REGION_BYTES.
        // For typical regions this produces a single batch (no overhead).
        // For very large BAM files where BAI bins span >256 MiB of compressed
        // data, this splits the work into multiple RegionBuf loads.
        let batches = partition_chunks(&chunks, region_buf::MAX_REGION_BYTES);

        for batch in &batches {
            let mut region = RegionBuf::load(&mut self.bulk_reader, batch)?;

            for chunk in batch {
                region.seek_virtual(chunk.begin)?;

                loop {
                    let current_voff = region.virtual_offset();
                    if current_voff >= chunk.end {
                        break;
                    }

                    // r[impl bam.reader.propagate_errors]
                    let raw = match region.read_record(&mut scratch) {
                        Ok(s) => s,
                        Err(BgzfError::UnexpectedEof) => break,
                        Err(e) => return Err(BamError::from(e)),
                    };

                    if raw.len() < 32 {
                        return Err(BamError::TruncatedRecord { offset: current_voff.0 });
                    }

                    // r[impl bam.reader.unmapped_skipped+2]
                    debug_assert!(raw.len() >= 32, "raw record too short: {}", raw.len());
                    #[allow(clippy::indexing_slicing, reason = "raw.len() >= 32 checked above")]
                    let rec_tid = i32::from_le_bytes([raw[0], raw[1], raw[2], raw[3]]);
                    #[allow(clippy::indexing_slicing, reason = "raw.len() >= 32 checked above")]
                    let rec_pos_raw = i32::from_le_bytes([raw[4], raw[5], raw[6], raw[7]]);
                    let rec_pos = Pos0::try_from(rec_pos_raw)
                        .map_err(|_| BamError::InvalidPosition { value: rec_pos_raw })?;
                    #[allow(clippy::indexing_slicing, reason = "raw.len() >= 32 checked above")]
                    let rec_flags = BamFlags::from(u16::from_le_bytes([raw[14], raw[15]]));

                    if rec_tid != tid_i32 {
                        skipped_tid = skipped_tid.saturating_add(1);
                        continue;
                    }

                    if rec_flags.is_unmapped() {
                        skipped_unmapped = skipped_unmapped.saturating_add(1);
                        continue;
                    }

                    let rec_end = compute_end_pos_from_raw(raw).unwrap_or(rec_pos);
                    if rec_pos > end || rec_end < start {
                        skipped_out_of_range = skipped_out_of_range.saturating_add(1);
                        continue;
                    }

                    accepted = accepted.saturating_add(1);
                    if store.push_raw(raw, customize)?.is_some() {
                        kept_count = kept_count.saturating_add(1);
                    }
                }
            }
        }

        Ok(FetchCounts { fetched: accepted as usize, kept: kept_count })
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
