//! Parse CRAI indexes to locate slices within a CRAM file. [`CramIndex`] maps reference IDs
//! to [`CraiEntry`] lists with byte offsets and alignment spans for region queries.

use super::reader::CramError;
use std::io::Read;
use std::path::{Path, PathBuf};

#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum CramIndexError {
    #[error("failed to decompress CRAI index")]
    DecompressionFailed,

    #[error("expected 6 TAB-separated fields, got {found}")]
    InvalidFieldCount { found: usize },

    #[error("invalid CRAI field {name}: {value:?}")]
    InvalidField { name: &'static str, value: String },
}

/// A single entry in a CRAI index.
#[derive(Debug, Clone)]
pub struct CraiEntry {
    pub ref_id: i32,
    pub alignment_start: i64,
    pub alignment_span: i64,
    pub container_offset: u64,
    pub slice_offset: u64,
    pub slice_size: u64,
}

/// Parsed CRAM index for region-based random access.
#[derive(Debug)]
pub struct CramIndex {
    /// Sorted by reference, then by 0-based start.
    entries: Vec<CraiEntry>,
    /// Per entry, the furthest 0-based last position it or an earlier entry
    /// of its reference reaches — `u64::MAX` from the first entry of unknown
    /// extent on. Non-decreasing within a reference, so a query finds the
    /// first entry that can reach it by binary search.
    reach: Vec<u64>,
}

impl CraiEntry {
    /// The entry's 0-based start. CRAI stores a **1-based** one; `0` only
    /// marks unmapped entries (see [`Self::placed_last`]).
    fn start0(&self) -> u64 {
        self.alignment_start.unsigned_abs().saturating_sub(1)
    }

    /// The 0-based last position the entry reaches: `None` for an unmapped
    /// entry, `u64::MAX` for one of unknown extent.
    fn placed_last(&self) -> Option<u64> {
        // r[impl cram.index.unmapped]
        if self.alignment_start == 0 && self.alignment_span == 0 {
            return None;
        }
        // r[impl cram.index.zero_span+2]
        if self.alignment_span == 0 {
            return Some(u64::MAX);
        }
        Some(self.start0().saturating_add(self.alignment_span.unsigned_abs()).saturating_sub(1))
    }
}

impl CramIndex {
    /// Parse a `.crai` index file (gzip-compressed TSV).
    // r[impl cram.index.parse]
    pub fn from_path(path: &Path) -> Result<Self, CramError> {
        let compressed = std::fs::read(path)
            .map_err(|source| CramError::Open { path: path.to_path_buf(), source })?;

        let mut decoder = flate2::read::GzDecoder::new(&compressed[..]);
        let mut text = String::new();
        decoder
            .read_to_string(&mut text)
            .map_err(|_| CramError::from(CramIndexError::DecompressionFailed))?;

        Self::from_text_lines(&text)
    }

    #[cfg(feature = "fuzz")]
    /// Parse a CRAI index from uncompressed TSV text (skips gzip decompression).
    pub fn from_text(text: &str) -> Result<Self, CramError> {
        Self::from_text_lines(text)
    }

    fn from_text_lines(text: &str) -> Result<Self, CramError> {
        let entries = text
            .lines()
            .filter(|line| !line.is_empty())
            .map(parse_crai_line)
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self::from_entries(entries))
    }

    fn from_entries(mut entries: Vec<CraiEntry>) -> Self {
        entries.sort_by_key(|e| (e.ref_id, e.start0()));
        let mut reach = Vec::with_capacity(entries.len());
        let mut previous: Option<(i32, u64)> = None;
        for e in &entries {
            let before = previous.filter(|&(ref_id, _)| ref_id == e.ref_id).map(|(_, r)| r);
            let furthest = before.max(e.placed_last()).unwrap_or(0);
            reach.push(furthest);
            previous = Some((e.ref_id, furthest));
        }
        CramIndex { entries, reach }
    }

    /// Find all index entries whose range overlaps the 0-based inclusive region
    /// `[start, end]` for the given reference ID. Returns entries sorted by
    /// start.
    ///
    /// CRAI stores a **1-based** `alignment_start`, so it is converted here
    /// before the comparison — the two conventions are not interchangeable even
    /// in a symmetric overlap test, and reading the entry one position too high
    /// prunes the container holding a narrow query's only records.
    // r[impl cram.index.query+2]
    // r[impl cram.index.unmapped]
    // r[impl interval.overlap_test]
    pub fn query(&self, tid: i32, start: u64, end: u64) -> Vec<&CraiEntry> {
        let (lo, hi) = self.tid_range(tid);
        let (Some(entries), Some(reach)) = (self.entries.get(lo..hi), self.reach.get(lo..hi))
        else {
            return Vec::new();
        };
        // Past the query end nothing starts in it; before the first entry
        // whose reach gets to the query start nothing reaches into it.
        let last = entries.partition_point(|e| e.start0() <= end);
        let first = reach.get(..last).map_or(0, |reach| reach.partition_point(|&r| r < start));
        entries
            .get(first..last)
            .unwrap_or_default()
            .iter()
            .filter(|e| e.placed_last().is_some_and(|reach| reach >= start))
            .collect()
    }

    // r[impl cram.index.region_bytes]
    /// Estimate the compressed bytes of the records a query for `[start, end]`
    /// (0-based, inclusive) on `tid` returns: each selected slice's size,
    /// prorated by the share of its span the query covers.
    pub fn region_bytes(&self, tid: i32, start: u64, end: u64) -> u64 {
        self.query(tid, start, end)
            .into_iter()
            .map(|e| {
                let Some(last) = e.placed_last().filter(|&last| last != u64::MAX) else {
                    return e.slice_size;
                };
                let first = e.start0();
                let covered = last.min(end).saturating_sub(first.max(start)).saturating_add(1);
                let span = last.saturating_sub(first).saturating_add(1);
                // Both fit in u64, so the product fits in u128; span >= 1.
                let bytes = u128::from(e.slice_size)
                    .saturating_mul(u128::from(covered))
                    .checked_div(u128::from(span))
                    .unwrap_or(0);
                u64::try_from(bytes).unwrap_or(u64::MAX)
            })
            .fold(0, u64::saturating_add)
    }

    /// The positions in `entries` of reference `tid`'s entries.
    fn tid_range(&self, tid: i32) -> (usize, usize) {
        let lo = self.entries.partition_point(|e| e.ref_id < tid);
        let hi = self.entries.partition_point(|e| e.ref_id <= tid);
        (lo, hi)
    }

    /// Every entry for reference `tid`, placed or not, sorted by
    /// `alignment_start`.
    pub fn entries_for(&self, tid: i32) -> &[CraiEntry] {
        let (lo, hi) = self.tid_range(tid);
        self.entries.get(lo..hi).unwrap_or_default()
    }

    /// Get all entries (for debugging/testing).
    pub fn entries(&self) -> &[CraiEntry] {
        &self.entries
    }
}

/// Find the `.crai` index file for a CRAM file.
pub fn find_crai_path(cram_path: &Path) -> Result<PathBuf, CramError> {
    // Try <file>.cram.crai first, then <file>.crai
    let mut crai = cram_path.to_path_buf();
    crai.set_extension("cram.crai");
    if crai.exists() {
        return Ok(crai);
    }

    // Try appending .crai to the full path
    let with_crai = PathBuf::from(format!("{}.crai", cram_path.display()));
    if with_crai.exists() {
        return Ok(with_crai);
    }

    Err(CramError::IndexNotFound { cram_path: cram_path.to_path_buf() })
}

fn parse_crai_line(line: &str) -> Result<CraiEntry, CramError> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 6 {
        return Err(CramIndexError::InvalidFieldCount { found: fields.len() }.into());
    }

    let parse_field = |idx: usize, name: &'static str| -> Result<i64, CramError> {
        debug_assert!(idx < fields.len(), "field index out of bounds: {idx} >= {}", fields.len());
        #[allow(clippy::indexing_slicing, reason = "bounds checked above")]
        fields[idx].parse::<i64>().map_err(|_| {
            CramIndexError::InvalidField { name, value: fields[idx].to_string() }.into()
        })
    };

    #[expect(
        clippy::cast_possible_truncation,
        reason = "CRAI ref_id is an ITF8 value bounded by i32 range; parsed as i64 for uniformity"
    )]
    let ref_id = parse_field(0, "ref_id")? as i32;
    Ok(CraiEntry {
        ref_id,
        alignment_start: parse_field(1, "alignment_start")?,
        alignment_span: parse_field(2, "alignment_span")?,
        container_offset: parse_field(3, "container_offset")? as u64,
        slice_offset: parse_field(4, "slice_offset")? as u64,
        slice_size: parse_field(5, "slice_size")? as u64,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use hegel::prelude::*;
    use tempfile::tempdir;

    // r[verify cram.index.parse]
    #[test]
    fn parse_real_crai() {
        let crai_path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.cram.crai");
        let index = CramIndex::from_path(Path::new(crai_path)).unwrap();
        let entries = index.entries();

        assert!(!entries.is_empty(), "CRAI should have entries");

        // All entries should have valid container offsets
        for entry in entries {
            assert!(entry.container_offset > 0, "container offset should be > 0");
        }
    }

    // r[verify cram.index.query+2]
    #[test]
    fn query_crai_for_known_region() {
        let crai_path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.cram.crai");
        let index = CramIndex::from_path(Path::new(crai_path)).unwrap();

        // Find a non-unmapped entry to query
        let first = index.entries().iter().find(|e| e.ref_id >= 0).unwrap();
        let tid = first.ref_id;

        // Query a broad region that should overlap
        let results = index.query(tid, 0, u64::MAX);
        assert!(!results.is_empty(), "should find entries for tid={tid}");
    }

    // r[verify cram.index.unmapped]
    #[test]
    fn query_crai_non_overlapping() {
        let crai_path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.cram.crai");
        let index = CramIndex::from_path(Path::new(crai_path)).unwrap();

        // Query for a tid that doesn't exist
        let results = index.query(9999, 0, u64::MAX);
        assert!(results.is_empty(), "should find no entries for non-existent tid");
    }

    #[test]
    fn find_crai_path_for_test_cram() {
        let cram_path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.cram");
        let found = find_crai_path(Path::new(cram_path)).unwrap();
        assert!(found.exists());
    }

    #[test]
    fn find_crai_path_missing() {
        let result = find_crai_path(Path::new("/nonexistent/file.cram"));
        assert!(matches!(result, Err(CramError::IndexNotFound { .. })));
    }

    // r[verify cram.index.zero_span+2]
    #[test]
    fn query_with_zero_span_entries() {
        // CRAI entries with span=0 occur when samtools writes CRAM with
        // embedded references or when the span is unknown. These entries
        // must still be returned by queries that overlap their start position.
        let index = CramIndex::from_entries(vec![
            CraiEntry {
                ref_id: 0,
                alignment_start: 1000,
                alignment_span: 0,
                container_offset: 100,
                slice_offset: 0,
                slice_size: 500,
            },
            CraiEntry {
                ref_id: 0,
                alignment_start: 2000,
                alignment_span: 500,
                container_offset: 600,
                slice_offset: 0,
                slice_size: 500,
            },
        ]);

        // Query [500, 1500) — entry_start=1000 is within query, should match
        let results = index.query(0, 500, 1500);
        assert!(!results.is_empty(), "span=0 entry within range should match");

        // Query [1500, 2500) — entry_start=1000 is before query start.
        // With span=0, entry_end=1000 which is before query start, so this
        // would be missed. But the slice actually covers records well past 1000.
        // span=0 means "unknown extent" — the entry MUST be included.
        let results = index.query(0, 1500, 2500);
        assert_eq!(results.len(), 2, "span=0 entry should be included when span is unknown");

        // Query [3000, 4000) should match nothing (span=0 entry at 1000 < 4000
        // but it's still included, however the span=500 entry at 2000 ends at 2500 < 3000)
        // Actually span=0 entry has start=1000 < end=4000, so it IS included.
        // Only entries past query end are excluded.
        let results = index.query(0, 3000, 4000);
        assert_eq!(results.len(), 1, "span=0 entry should still match (start < query_end)");
    }

    #[hegel::composite]
    fn arb_entry(tc: &TestCase) -> CraiEntry {
        CraiEntry {
            ref_id: tc.draw_silent(gs::integers::<i32>().min_value(-1).max_value(2)),
            // 0 with span 0 is an unmapped entry; span 0 alone is unknown extent.
            alignment_start: tc.draw_silent(gs::integers::<i64>().min_value(0).max_value(60)),
            alignment_span: tc.draw_silent(gs::sampled_from(&[0, 1, 2, 5, 20])),
            container_offset: tc.draw_silent(gs::integers::<u64>()),
            slice_offset: tc.draw_silent(gs::integers::<u64>()),
            slice_size: 1,
        }
    }

    // r[verify cram.index.query+2]
    // r[verify cram.index.zero_span+2]
    // r[verify cram.index.unmapped]
    /// The binary-searched query returns what testing every entry of the
    /// reference for overlap returns, in the same order.
    #[hegel::test]
    #[allow(
        clippy::arithmetic_side_effects,
        clippy::cast_sign_loss,
        reason = "generated starts and spans are small and non-negative"
    )]
    fn query_matches_a_scan_of_every_entry(tc: TestCase) {
        let index =
            CramIndex::from_entries(tc.draw(gs::vecs(arb_entry()).max_size(40).print_as_debug()));
        let tid = tc.draw(gs::integers::<i32>().min_value(-1).max_value(3));
        let start = tc.draw(gs::integers::<u64>().max_value(80));
        let end = tc.draw(gs::integers::<u64>().min_value(start).max_value(80));

        let scanned: Vec<&CraiEntry> = index
            .entries()
            .iter()
            .filter(|e| {
                // 1-based start and span, as the CRAI stores them.
                let unmapped = e.alignment_start == 0 && e.alignment_span == 0;
                let first = (e.alignment_start as u64).saturating_sub(1);
                let reaches = e.alignment_span == 0
                    || (first + e.alignment_span as u64).saturating_sub(1) >= start;
                e.ref_id == tid && !unmapped && first <= end && reaches
            })
            .collect();
        let queried = index.query(tid, start, end);
        assert_eq!(queried.len(), scanned.len());
        assert!(queried.iter().zip(&scanned).all(|(a, b)| std::ptr::eq(*a, *b)));
    }

    // r[verify cram.index.region_bytes]
    /// The estimate never shrinks as the query grows, and a query over a
    /// whole reference counts every slice on it once, whole.
    #[hegel::test]
    fn region_bytes_grows_with_the_query_to_every_slice(tc: TestCase) {
        let index = CramIndex::from_entries(
            tc.draw(gs::vecs(arb_sized_entry()).max_size(40).print_as_debug()),
        );
        let tid = tc.draw(gs::integers::<i32>().min_value(0).max_value(2));
        let start = tc.draw(gs::integers::<u64>().max_value(80));
        let end = tc.draw(gs::integers::<u64>().min_value(start).max_value(80));
        let further = tc.draw(gs::integers::<u64>().min_value(end).max_value(80));

        assert!(index.region_bytes(tid, start, end) <= index.region_bytes(tid, start, further));
        let whole: u64 = index
            .entries_for(tid)
            .iter()
            .filter(|e| e.alignment_start != 0 || e.alignment_span != 0)
            .map(|e| e.slice_size)
            .sum();
        assert_eq!(index.region_bytes(tid, 0, u64::MAX), whole);
    }

    #[hegel::composite]
    fn arb_sized_entry(tc: &TestCase) -> CraiEntry {
        CraiEntry {
            slice_size: tc.draw_silent(gs::integers::<u64>().max_value(1 << 20)),
            ..tc.draw_silent(arb_entry())
        }
    }

    // r[verify cram.index.region_bytes]
    #[test]
    fn region_bytes_prorates_a_slice_by_the_share_covered() {
        let entry = |start, span, slice_size| CraiEntry {
            ref_id: 0,
            alignment_start: start,
            alignment_span: span,
            container_offset: 0,
            slice_offset: 0,
            slice_size,
        };
        // 1-based 101..=200 and 201..=300, then one of unknown extent at 1001.
        let index = CramIndex::from_entries(vec![
            entry(101, 100, 1000),
            entry(201, 100, 3000),
            entry(1001, 0, 50),
        ]);
        assert_eq!(index.region_bytes(0, 100, 124), 250, "a quarter of the first");
        assert_eq!(index.region_bytes(0, 150, 249), 500 + 1500, "half of each");
        assert_eq!(index.region_bytes(0, 0, 999), 4000);
        assert_eq!(index.region_bytes(0, 5000, 6000), 50, "unknown extent counts whole");
        assert_eq!(index.region_bytes(1, 0, 6000), 0);
    }

    #[test]
    fn parse_crai_line_valid() {
        let entry = parse_crai_line("0\t100\t500\t1234\t0\t5678").unwrap();
        assert_eq!(entry.ref_id, 0);
        assert_eq!(entry.alignment_start, 100);
        assert_eq!(entry.alignment_span, 500);
        assert_eq!(entry.container_offset, 1234);
        assert_eq!(entry.slice_offset, 0);
        assert_eq!(entry.slice_size, 5678);
    }

    #[test]
    fn crai_entries_match_cram_containers() {
        // Verify that CRAI entries point to valid container offsets in the CRAM file
        let crai_path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.cram.crai");
        let cram_data =
            std::fs::read(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.cram"))
                .unwrap();
        let index = CramIndex::from_path(Path::new(crai_path)).unwrap();

        for entry in index.entries() {
            if entry.ref_id < 0 {
                continue; // skip unmapped
            }
            // Verify the container offset points to a valid container header
            #[allow(
                clippy::cast_possible_truncation,
                reason = "test runs on 64-bit; offset verified < file size"
            )]
            let offset = entry.container_offset as usize;
            assert!(offset < cram_data.len(), "container offset out of bounds");
            #[allow(
                clippy::indexing_slicing,
                reason = "offset is verified < cram_data.len() above"
            )]
            let container =
                super::super::container::ContainerHeader::parse(&cram_data[offset..]).unwrap();
            assert!(container.num_records > 0, "CRAI entry points to empty container");
        }
    }

    #[test]
    fn parse_crai_line_invalid_field_count() {
        let err = parse_crai_line("too\tfew\tfields").unwrap_err();
        assert!(
            matches!(
                err,
                CramError::IndexParse { source: CramIndexError::InvalidFieldCount { found: 3 } }
            ),
            "expected InvalidFieldCount {{ found: 3 }}, got {err:?}"
        );
    }

    #[test]
    fn parse_crai_line_invalid_field_value() {
        let err = parse_crai_line("not_a_number\t0\t0\t0\t0\t0").unwrap_err();
        assert!(
            matches!(
                &err,
                CramError::IndexParse {
                    source: CramIndexError::InvalidField { name: "ref_id", .. }
                }
            ),
            "expected InvalidField {{ name: \"ref_id\", .. }}, got {err:?}"
        );
    }

    #[test]
    fn from_path_decompression_failed() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("seqair_test_bad.crai");
        std::fs::write(&path, b"this is not valid gzip data").unwrap();
        let err = CramIndex::from_path(&path).unwrap_err();
        assert!(
            matches!(err, CramError::IndexParse { source: CramIndexError::DecompressionFailed }),
            "expected DecompressionFailed, got {err:?}"
        );
        std::fs::remove_file(&path).ok();
    }
}
