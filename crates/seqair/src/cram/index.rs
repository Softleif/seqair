//! Parse CRAI indexes to locate slices within a CRAM file. [`CramIndex`] maps reference IDs
//! to [`CraiEntry`] lists with byte offsets and alignment spans for region queries.

use super::reader::CramError;
use std::io::Read;
use std::path::{Path, PathBuf};

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
    entries: Vec<CraiEntry>,
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

        let mut entries = Vec::new();
        for line in text.lines() {
            if line.is_empty() {
                continue;
            }
            let entry = parse_crai_line(line)?;
            entries.push(entry);
        }

        // Sort by (ref_id, alignment_start) for efficient querying
        entries.sort_by(|a, b| {
            a.ref_id.cmp(&b.ref_id).then(a.alignment_start.cmp(&b.alignment_start))
        });

        Ok(CramIndex { entries })
    }

    #[cfg(feature = "fuzz")]
    /// Parse a CRAI index from uncompressed TSV text (skips gzip decompression).
    pub fn from_text(text: &str) -> Result<Self, CramError> {
        let mut entries = Vec::new();
        for line in text.lines() {
            if line.is_empty() {
                continue;
            }
            let entry = parse_crai_line(line)?;
            entries.push(entry);
        }

        entries.sort_by(|a, b| {
            a.ref_id.cmp(&b.ref_id).then(a.alignment_start.cmp(&b.alignment_start))
        });

        Ok(CramIndex { entries })
    }

    /// Find all index entries whose range overlaps `[start, end)` for the given
    /// reference ID. Returns entries sorted by `alignment_start`.
    ///
    /// Positions are 0-based half-open (matching the query convention).
    /// CRAI entries use 1-based positions internally, but the comparison
    /// handles both correctly since we just check overlap.
    // r[impl cram.index.query]
    // r[impl cram.index.unmapped]
    pub fn query(&self, tid: i32, start: u64, end: u64) -> Vec<&CraiEntry> {
        self.entries
            .iter()
            .filter(|e| {
                if e.ref_id != tid {
                    return false;
                }
                // r[impl cram.index.unmapped]
                if e.alignment_start == 0 && e.alignment_span == 0 {
                    return false;
                }
                let entry_start = e.alignment_start as u64;
                // r[impl cram.index.zero_span]
                if e.alignment_span == 0 {
                    return entry_start < end;
                }
                let entry_end = entry_start.saturating_add(e.alignment_span as u64);
                entry_start < end && entry_end > start
            })
            .collect()
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

    // r[verify cram.index.query]
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

    // r[verify cram.index.zero_span]
    #[test]
    fn query_with_zero_span_entries() {
        // CRAI entries with span=0 occur when samtools writes CRAM with
        // embedded references or when the span is unknown. These entries
        // must still be returned by queries that overlap their start position.
        let index = CramIndex {
            entries: vec![
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
            ],
        };

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
