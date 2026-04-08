//! Open and query BAM files. [`IndexedBamReader`] parses the BAI index, then fetches records
//! for a region into a [`RecordStore`]. Call [`IndexedBamReader::fork`] to get a cheap
//! per-thread reader that shares the parsed index and header via [`Arc<BamShared>`].

use super::{
    bgzf::{BgzfError, BgzfReader},
    flags::FLAG_UNMAPPED,
    header::{BamHeader, BamHeaderError},
    index::{BaiError, BamIndex},
    record::{DecodeError, compute_end_pos_from_raw},
    record_store::RecordStore,
    region_buf::RegionBuf,
};
use seqair_types::{Pos, SmolStr, Zero};
use std::{
    fs::File,
    io::{Read, Seek},
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
        let index = BamIndex::from_bytes(bai_data)?;

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
        start: Pos<Zero>,
        end: Pos<Zero>,
        store: &mut RecordStore,
    ) -> Result<usize, BamError> {
        store.clear();

        let chunks = self.shared.index.query(tid, start, end);
        if chunks.is_empty() {
            return Ok(0);
        }

        let tid_i32 = validate_tid(tid)?;

        let mut skipped_tid: u32 = 0;
        let mut skipped_unmapped: u32 = 0;
        let mut skipped_out_of_range: u32 = 0;
        let mut accepted: u32 = 0;

        {
            let mut region = RegionBuf::load(&mut self.bulk_reader, &chunks)?;

            // Scratch buffer for the rare case where a record body straddles a BGZF
            // block boundary; zero-copy slice from RegionBuf::buf is used otherwise.
            let mut scratch: Vec<u8> = Vec::new();

            for chunk in &chunks {
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
                        skipped_tid = skipped_tid.saturating_add(1);
                        continue;
                    }

                    if rec_flags & FLAG_UNMAPPED != 0 {
                        skipped_unmapped = skipped_unmapped.saturating_add(1);
                        continue;
                    }

                    let rec_end = compute_end_pos_from_raw(raw).unwrap_or(rec_pos);
                    if rec_pos > end || rec_end < start {
                        skipped_out_of_range = skipped_out_of_range.saturating_add(1);
                        continue;
                    }

                    accepted = accepted.saturating_add(1);
                    store.push_raw(raw)?;
                }
            }
        }

        tracing::debug!(
            target: super::region_buf::PROFILE_TARGET,
            accepted,
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
}
