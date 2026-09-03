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

    /// Estimate the **compressed** bytes a `[start, end]` region query would
    /// load into a [`RegionBuf`] — the merged byte span of its index chunks.
    ///
    /// This is what byte-aware segmentation budgets against. It's only as fine
    /// as the index (~16 kb leaf bins), so a sub-bin range still reports its
    /// whole leaf bin's size.
    pub fn estimate_region_bytes(&self, tid: u32, start: Pos0, end: Pos0) -> u64 {
        let chunks = self.shared.index.query(tid, start, end);
        region_buf::merged_byte_size(&chunks) as u64
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
    /// 32-byte BAM fixed header) for records overlapping the inclusive region
    /// `[start, end]` on `tid` (r[`interval.overlap_test`]).
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

        // One streaming RegionBuf spans all chunks: the sliding window bounds
        // peak memory, so there's no need to pre-partition into batches.
        let mut region = RegionBuf::new(&mut self.bulk_reader, &chunks)?;
        let chunk_end = match chunks.first() {
            Some(first) => {
                region.seek_virtual(first.begin)?;
                first.end
            }
            None => VirtualOffset(0),
        };

        Ok(BamQuery {
            region,
            chunks,
            chunk_idx: 0,
            chunk_end,
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
    /// Streaming buffer spanning all of the query's chunks; owns the reader borrow.
    region: RegionBuf<'r, R>,
    /// All index chunks for the query, in begin order.
    chunks: Vec<Chunk>,
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
            .field("chunks", &self.chunks.len())
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

    /// Drive the cursor: a flat `loop` + `continue` state machine over two
    /// phases, resumable so `for_each` can stop early. State persists in
    /// `self.{region, chunk_idx, chunk_end}` between iterations.
    ///
    /// 1. **Advance** — the cursor passed the current chunk's `end`: step to
    ///    the next chunk and seek the streaming buffer to its `begin` (skipping
    ///    any bytes between chunks). The window refills lazily as needed.
    /// 2. **Record** — read one record from the streaming buffer, filter by tid
    ///    then position, and hand accepted records to `f`.
    ///
    /// The phases are mutually exclusive per iteration; each ends in `continue`
    /// (or `break`) so the next iteration re-evaluates from the top.
    fn run_loop<F, E>(&mut self, mut f: F) -> Result<BamQueryCounts, E>
    where
        F: FnMut(&[u8]) -> Result<(), E>,
        E: From<BamError>,
    {
        loop {
            if self.chunk_idx >= self.chunks.len() {
                break; // all chunks consumed → done
            }

            // ── Phase 1: cursor reached this chunk's end → advance ───────────
            if self.region.virtual_offset() >= self.chunk_end {
                self.chunk_idx = self.chunk_idx.saturating_add(1);
                let Some(chunk) = self.chunks.get(self.chunk_idx) else {
                    break; // no more chunks
                };
                let chunk = *chunk;
                self.region.seek_virtual(chunk.begin).map_err(|e| E::from(BamError::from(e)))?;
                self.chunk_end = chunk.end;
                continue;
            }

            // ── Phase 2: read one record from the streaming buffer ───────────
            // r[impl bam.reader.propagate_errors]
            let current_voff = self.region.virtual_offset();
            let raw = match self.region.read_record(&mut self.scratch) {
                Ok(s) => s,
                // Buffer exhausted before chunk_end (e.g. a record straddling the
                // loaded range's tail): loop back so phase 2 advances the chunk.
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

            // Skip reads on other references before touching position. Unplaced
            // unmapped reads carry tid = -1 / pos = -1 and can appear at chunk
            // boundaries (the unmapped block trails the last contig's records),
            // so checking tid first also keeps us from parsing their reserved
            // pos = -1, which would otherwise error the whole query.
            if rec_tid != self.tid_i32 {
                self.skipped_tid = self.skipped_tid.saturating_add(1);
                continue;
            }

            #[allow(clippy::indexing_slicing, reason = "raw.len() >= 32 checked above")]
            let rec_pos_raw = i32::from_le_bytes([raw[4], raw[5], raw[6], raw[7]]);
            // A record claiming this reference but with pos = -1 can't overlap
            // any region; skip it rather than erroring out the whole query.
            if rec_pos_raw < 0 {
                self.skipped_out_of_range = self.skipped_out_of_range.saturating_add(1);
                continue;
            }
            let rec_pos = Pos0::try_from(rec_pos_raw)
                .map_err(|_| E::from(BamError::InvalidPosition { value: rec_pos_raw }))?;
            #[allow(clippy::indexing_slicing, reason = "raw.len() >= 32 checked above")]
            let rec_flags = BamFlags::from(u16::from_le_bytes([raw[14], raw[15]]));

            // r[impl record_store.end_pos_htslib]
            let rec_end = if rec_flags.is_unmapped() {
                rec_pos
            } else {
                compute_end_pos_from_raw(raw).unwrap_or(rec_pos)
            };
            // r[impl interval.overlap_test]
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

    /// Frame a record with its 4-byte `block_size` prefix, as it appears on the
    /// BAM wire (what [`RegionBuf::read_record`] expects to find).
    fn framed_record(rec: &super::super::owned_record::OwnedBamRecord) -> Vec<u8> {
        let mut body = Vec::new();
        rec.to_bam_bytes(&mut body).unwrap();
        let mut out = i32::try_from(body.len()).unwrap().to_le_bytes().to_vec();
        out.extend_from_slice(&body);
        out
    }

    // r[verify bam.reader.unmapped_skipped]
    /// When a query's chunk range extends past the queried contig's reads into
    /// the trailing block of unplaced reads (tid = -1 / pos = -1), those records
    /// must be skipped, not error the whole query. This reproduces the `GRCh38`
    /// ALT/decoy crash where the last contigs' BAI chunk end abuts the unmapped
    /// block at the end of the file.
    ///
    /// We build `BamQuery` directly (rather than via a real index) because a
    /// correctly-built index keeps `chunk_end` tight; the bug only surfaces when
    /// the scanned range reaches into the unplaced block, which we model with an
    /// explicit chunk covering both records.
    #[test]
    fn query_skips_trailing_unplaced_reads() {
        use super::super::owned_record::OwnedBamRecord;
        use crate::io::BgzfWriter;
        use seqair_types::Base;

        let mapped = OwnedBamRecord::builder(0, Some(Pos0::new(100).unwrap()), b"mapped".to_vec())
            .seq(vec![Base::A, Base::C, Base::G])
            .build()
            .unwrap();
        // Fully unplaced read: ref_id = -1, pos = -1 (BAM wire), unmapped flag.
        let unplaced = OwnedBamRecord::builder(-1, None, b"unplaced".to_vec())
            .flags(BamFlags::from(0x4))
            .seq(vec![Base::A, Base::C, Base::G])
            .build()
            .unwrap();

        // One BGZF block holding [mapped(tid=0), unplaced(tid=-1/pos=-1)].
        let mut payload = framed_record(&mapped);
        payload.extend_from_slice(&framed_record(&unplaced));
        let payload_len = u16::try_from(payload.len()).expect("payload fits one block");

        let mut bgzf = BgzfWriter::new(Vec::new());
        bgzf.write_all(&payload).unwrap();
        let bam_bytes = bgzf.finish().unwrap();

        // A single chunk whose end sits past the mapped record, so the scan
        // reads on into the unplaced record — the production index quirk.
        let batch = vec![Chunk {
            begin: VirtualOffset::new(0, 0),
            end: VirtualOffset::new(0, payload_len),
        }];

        let mut reader = std::io::Cursor::new(bam_bytes);
        let mut region = RegionBuf::new(&mut reader, &batch).unwrap();
        region.seek_virtual(batch[0].begin).unwrap();
        let mut query = BamQuery {
            region,
            chunks: batch,
            chunk_idx: 0,
            chunk_end: VirtualOffset::new(0, payload_len),
            scratch: Vec::new(),
            tid_i32: 0,
            start: Pos0::new(0).unwrap(),
            end: Pos0::new(1000).unwrap(),
            accepted: 0,
            skipped_tid: 0,
            skipped_out_of_range: 0,
        };

        let mut seen = 0usize;
        let counts = query
            .for_each(|_raw| seen += 1)
            .expect("query must not error on trailing unplaced reads");

        assert_eq!(seen, 1, "only the mapped record should be yielded");
        assert_eq!(counts.fetched, 1);
        assert_eq!(counts.skipped_tid, 1, "the unplaced read must be skipped by tid");
    }

    use super::super::bgzf::VirtualOffset;
}
