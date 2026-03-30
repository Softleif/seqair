//! Open and query BAM files. [`IndexedBamReader`] parses the BAI index, then fetches records
//! for a region into a [`RecordStore`]. Call [`IndexedBamReader::fork`] to get a cheap
//! per-thread reader that shares the parsed index and header via [`Arc<BamShared>`].

use super::{
    bgzf::{BgzfError, BgzfReader},
    flags::FLAG_UNMAPPED,
    header::{BamHeader, BamHeaderError},
    index::{BaiError, BamIndex},
    record::DecodeError,
    record::compute_end_pos_from_raw,
    record_store::RecordStore,
    region_buf::RegionBuf,
};
use seqair_types::{Pos, SmolStr, Zero};
use std::{
    fs::File,
    path::{Path, PathBuf},
    sync::Arc,
};
use tracing::instrument;

// r[impl io.errors]
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

    #[error("truncated BAM record at virtual offset {offset:#x}")]
    TruncatedRecord { offset: u64 },

    #[error("contig not found: {name}")]
    ContigNotFound { name: SmolStr },

    #[error("region {contig}:{start}-{end} is out of bounds (contig length: {contig_len})")]
    RegionOutOfBounds { contig: SmolStr, start: u64, end: u64, contig_len: u64 },

    #[error("BAI index not found for {bam_path}")]
    IndexNotFound { bam_path: PathBuf },

    #[error("record decode error: {source}")]
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
    index: BamIndex,
    header: BamHeader,
    pub bam_path: PathBuf,
}

impl BamShared {
    pub fn index(&self) -> &BamIndex {
        &self.index
    }
}

// r[impl bam.reader.open]
// r[impl bam.reader.header_access]
pub struct IndexedBamReader {
    /// Separate file handle for bulk region reads (unbuffered — RegionBuf does
    /// large sequential reads that don't benefit from BufReader).
    bulk_reader: File,
    shared: Arc<BamShared>,
    /// Per-tid cache for records from distant bins (levels 0–2, ≥8 Mbp).
    chunk_cache: ChunkCache,
}

impl std::fmt::Debug for IndexedBamReader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("IndexedBamReader").field("bam_path", &self.shared.bam_path).finish()
    }
}

impl IndexedBamReader {
    #[instrument(level = "debug", fields(path = %path.display()), err)]
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
            chunk_cache: ChunkCache::default(),
        })
    }

    // r[impl bam.reader.fork]
    // r[impl bam.reader.fork_independence]
    // r[impl bam.reader.fork_equivalence]
    // r[impl bam.reader.fork_concurrent]
    #[instrument(level = "debug", skip(self), fields(path = %self.shared.bam_path.display()), err)]
    pub fn fork(&self) -> Result<Self, BamError> {
        let bulk_file = File::open(&self.shared.bam_path)
            .map_err(|source| BamError::Open { path: self.shared.bam_path.clone(), source })?;

        Ok(IndexedBamReader {
            bulk_reader: bulk_file,
            shared: Arc::clone(&self.shared),
            chunk_cache: ChunkCache::default(),
        })
    }

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
    // r[impl bam.index.chunk_cache]
    #[instrument(level = "debug", skip(self, store), fields(tid, start, end), err)]
    pub fn fetch_into(
        &mut self,
        tid: u32,
        start: Pos<Zero>,
        end: Pos<Zero>,
        store: &mut RecordStore,
    ) -> Result<usize, BamError> {
        store.clear();

        let query = self.shared.index.query_split(tid, start, end);
        if query.nearby.is_empty() && query.distant.is_empty() {
            return Ok(0);
        }

        // Lazily load ALL bin 0 records for this tid (once per chromosome).
        if !query.distant.is_empty() && self.chunk_cache.tid != Some(tid) {
            let all_distant = self.shared.index.distant_chunks(tid);
            self.chunk_cache.load(&mut self.bulk_reader, tid, &all_distant)?;
        }

        let tid_i32 = validate_tid(tid)?;

        let mut skipped_tid: u32 = 0;
        let mut skipped_unmapped: u32 = 0;
        let mut skipped_out_of_range: u32 = 0;
        let mut accepted: u32 = 0;

        if !query.nearby.is_empty() {
            let mut region = RegionBuf::load(&mut self.bulk_reader, &query.nearby)?;

            // Scratch buffer for the rare case where a record body straddles a BGZF
            // block boundary; zero-copy slice from RegionBuf::buf is used otherwise.
            let mut scratch: Vec<u8> = Vec::new();

            for chunk in &query.nearby {
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
                    let rec_pos = Pos::<Zero>::try_from_i64(i64::from(rec_pos_raw))
                        .ok_or(BamError::InvalidPosition { value: rec_pos_raw })?;
                    #[allow(clippy::indexing_slicing, reason = "raw.len() >= 32 checked above")]
                    let rec_flags = u16::from_le_bytes([raw[14], raw[15]]);

                    if rec_tid != tid_i32 {
                        skipped_tid += 1;
                        continue;
                    }

                    if rec_flags & FLAG_UNMAPPED != 0 {
                        skipped_unmapped += 1;
                        continue;
                    }

                    let rec_end = compute_end_pos_from_raw(raw).unwrap_or(rec_pos);
                    if rec_pos > end || rec_end < start {
                        skipped_out_of_range += 1;
                        continue;
                    }

                    accepted += 1;
                    store.push_raw(raw).map_err(|source| BamError::RecordDecode { source })?;
                }
            }
        }

        // Inject matching records from the distant-bin cache.
        let cache_injected = self.chunk_cache.inject_overlapping(tid, start, end, store)?;
        accepted += cache_injected;

        tracing::debug!(
            target: super::region_buf::PROFILE_TARGET,
            accepted,
            cache_injected,
            skipped_tid,
            skipped_unmapped,
            skipped_out_of_range,
            "fetch_into",
        );

        tracing::debug!(
            target: super::region_buf::PROFILE_TARGET,
            records = store.len(),
            records_cap = store.records_capacity(),
            names_bytes = store.names_capacity(),
            bases_bytes = store.bases_capacity(),
            data_bytes = store.data_capacity(),
            "record_store",
        );

        Ok(store.len())
    }
}

// r[impl bam.index.chunk_cache]
// r[related bam.index.chunk_separation+2]
/// Caches records from distant bins (levels 0–2, covering ≥8 Mbp) for a
/// single reference sequence.
///
/// These bins contain reads near higher-level bin boundaries. Their chunks are
/// scattered across the BAM file (often GB away from the query's main data).
/// Loading them once per tid and scanning per-query is much cheaper than
/// seeking to their distant file offsets on every region fetch.
#[derive(Debug, Default)]
pub(crate) struct ChunkCache {
    /// The reference id these records belong to, or `None` if empty.
    tid: Option<u32>,
    /// Raw BAM record bodies (the bytes after the 4-byte block_size prefix).
    /// Stored as raw bytes so they can be `push_raw`'d into any RecordStore.
    records: Vec<Vec<u8>>,
}

impl ChunkCache {
    /// Inject all cached records that overlap [start, end] and match `tid`
    /// into the given store.
    pub fn inject_overlapping(
        &self,
        tid: u32,
        start: Pos<Zero>,
        end: Pos<Zero>,
        store: &mut RecordStore,
    ) -> Result<u32, BamError> {
        if self.tid != Some(tid) {
            return Ok(0);
        }
        let mut count = 0u32;
        for raw in &self.records {
            if raw.len() < 32 {
                continue;
            }
            debug_assert!(raw.len() >= 32, "cached record too short: {}", raw.len());
            #[allow(clippy::indexing_slicing, reason = "raw.len() >= 32 checked above")]
            let rec_tid = i32::from_le_bytes([raw[0], raw[1], raw[2], raw[3]]);
            // Safety: tid was validated by validate_tid before calling this method
            if rec_tid != tid as i32 {
                continue;
            }
            #[allow(clippy::indexing_slicing, reason = "raw.len() >= 32 checked above")]
            let rec_flags = u16::from_le_bytes([raw[14], raw[15]]);
            if rec_flags & FLAG_UNMAPPED != 0 {
                continue;
            }
            #[allow(clippy::indexing_slicing, reason = "raw.len() >= 32 checked above")]
            let rec_pos_raw = i32::from_le_bytes([raw[4], raw[5], raw[6], raw[7]]);
            let rec_pos = Pos::<Zero>::try_from_i64(i64::from(rec_pos_raw))
                .ok_or(BamError::InvalidPosition { value: rec_pos_raw })?;
            let rec_end = compute_end_pos_from_raw(raw).unwrap_or(rec_pos);
            if rec_pos > end || rec_end < start {
                continue;
            }
            store.push_raw(raw).map_err(|source| BamError::RecordDecode { source })?;
            count += 1;
        }
        Ok(count)
    }

    /// Load bin 0 chunks for a given tid. Replaces any previously cached data.
    pub fn load<R: std::io::Read + std::io::Seek>(
        &mut self,
        reader: &mut R,
        tid: u32,
        chunks: &[super::index::Chunk],
    ) -> Result<(), BamError> {
        self.tid = Some(tid);
        self.records.clear();

        if chunks.is_empty() {
            return Ok(());
        }

        let mut region = RegionBuf::load(reader, chunks)?;
        let mut scratch: Vec<u8> = Vec::new();

        for chunk_ref in chunks {
            if region.seek_virtual(chunk_ref.begin).is_err() {
                continue;
            }
            loop {
                let current_voff = region.virtual_offset();
                if current_voff >= chunk_ref.end {
                    break;
                }
                // r[impl bam.reader.propagate_errors]
                let raw = match region.read_record(&mut scratch) {
                    Ok(s) => s,
                    Err(BgzfError::UnexpectedEof) => break,
                    Err(e) => return Err(BamError::from(e)),
                };
                self.records.push(raw.to_vec());
            }
        }

        tracing::debug!(
            target: super::region_buf::PROFILE_TARGET,
            tid,
            cached_records = self.records.len(),
            "chunk_cache::load",
        );

        Ok(())
    }
}

// r[impl unified.detect_index]
// r[impl unified.detect_error]
fn find_and_open_index(bam_path: &Path) -> Result<BamIndex, BamError> {
    let bai_path = bam_path.with_extension("bam.bai");
    if bai_path.exists() {
        return BamIndex::from_path(&bai_path).map_err(|source| BamError::Index { source });
    }

    let mut bai_path2 = bam_path.to_path_buf();
    let mut name = bai_path2.file_name().unwrap_or_default().to_os_string();
    name.push(".bai");
    bai_path2.set_file_name(name);
    if bai_path2.exists() {
        return BamIndex::from_path(&bai_path2).map_err(|source| BamError::Index { source });
    }

    Err(BamError::IndexNotFound { bam_path: bam_path.to_path_buf() })
}

#[cfg(test)]
mod tests {
    use super::*;

    // r[verify bam.index.chunk_cache]
    #[test]
    fn chunk_cache_inject_filters_by_region() {
        // Construct minimal raw BAM records (32 bytes minimum).
        // Layout: tid(4) pos(4) name_len(1) mapq(1) bin(2) n_cigar(2) flags(2) seq_len(4) ...
        fn make_raw_record(tid: i32, pos: i32, flags: u16) -> Vec<u8> {
            let mut raw = vec![0u8; 36];
            // tid at offset 0
            raw[..4].copy_from_slice(&tid.to_le_bytes());
            // pos at offset 4
            raw[4..8].copy_from_slice(&pos.to_le_bytes());
            // name_len at offset 8 (must be >= 1 for qname)
            raw[8] = 1;
            // flags at offset 14
            raw[14..16].copy_from_slice(&flags.to_le_bytes());
            // n_cigar_op at offset 12 = 0 (no cigar → end_pos = pos)
            raw
        }

        let cache = ChunkCache {
            tid: Some(0),
            records: vec![
                make_raw_record(0, 100, 0),   // overlaps [50, 200]
                make_raw_record(0, 300, 0),   // outside [50, 200]
                make_raw_record(0, 150, 0x4), // unmapped, should be skipped
                make_raw_record(1, 100, 0),   // wrong tid
            ],
        };

        let mut store = RecordStore::new();
        let count = cache
            .inject_overlapping(
                0,
                Pos::<Zero>::new(50).unwrap(),
                Pos::<Zero>::new(200).unwrap(),
                &mut store,
            )
            .unwrap();

        assert_eq!(count, 1, "only the record at pos=100 on tid=0 should match");
        assert_eq!(store.len(), 1);
    }

    #[test]
    fn chunk_cache_empty_for_wrong_tid() {
        let cache = ChunkCache { tid: Some(5), records: vec![] };
        let mut store = RecordStore::new();
        let count = cache
            .inject_overlapping(
                3,
                Pos::<Zero>::new(0).unwrap(),
                Pos::<Zero>::new(1000).unwrap(),
                &mut store,
            )
            .unwrap();
        assert_eq!(count, 0);
    }

    #[test]
    fn chunk_cache_empty_when_uninitialized() {
        let cache = ChunkCache::default();
        let mut store = RecordStore::new();
        let count = cache
            .inject_overlapping(
                0,
                Pos::<Zero>::new(0).unwrap(),
                Pos::<Zero>::new(1000).unwrap(),
                &mut store,
            )
            .unwrap();
        assert_eq!(count, 0);
    }

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
}
