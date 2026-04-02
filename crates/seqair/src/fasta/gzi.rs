//! Parse GZI indexes to seek within bgzf-compressed FASTA files. [`GziIndex`] maps uncompressed
//! byte offsets to compressed block positions, enabling random access without full decompression.

use std::path::Path;
use tracing::instrument;

// r[impl io.errors]
#[derive(Debug, thiserror::Error)]
pub enum GziError {
    #[error("I/O error reading GZI index {path}")]
    Read { path: std::path::PathBuf, source: std::io::Error },

    #[error("GZI index {path} is too short ({actual} bytes, need at least 8 for the entry count)")]
    TooShort { path: std::path::PathBuf, actual: usize },

    #[error("GZI index {path}: entry count {count} causes size overflow")]
    SizeOverflow { path: std::path::PathBuf, count: u64 },

    #[error("GZI index {path}: expected {expected} bytes for {count} entries, got {actual}")]
    WrongSize { path: std::path::PathBuf, count: u64, expected: u64, actual: usize },

    #[error("GZI index {path}: entry {index} is truncated")]
    TruncatedEntry { path: std::path::PathBuf, index: usize },

    #[error(
        "GZI index {path}: entries are not sorted at entry {index} (uncompressed offset {offset} <= previous {prev_offset})"
    )]
    UnsortedEntries { path: std::path::PathBuf, index: usize, offset: u64, prev_offset: u64 },

    #[error(
        "GZI within-block offset {within_block} exceeds BGZF max block size (65535) — GZI index may be corrupt or missing entries"
    )]
    WithinBlockOverflow { within_block: u64 },

    #[error("internal GZI error: binary search returned out-of-bounds index {index}")]
    InternalBoundsError { index: usize },
}

/// An entry mapping a compressed block offset to an uncompressed byte offset.
#[derive(Debug, Clone, Copy)]
struct GziEntry {
    compressed_offset: u64,
    uncompressed_offset: u64,
}

// r[impl fasta.gzi.parse]
#[derive(Debug, Clone)]
pub struct GziIndex {
    entries: Vec<GziEntry>,
}

impl GziIndex {
    #[instrument(level = "debug")]
    pub fn from_file(path: &Path) -> Result<Self, GziError> {
        let data = std::fs::read(path)
            .map_err(|source| GziError::Read { path: path.to_path_buf(), source })?;
        Self::parse(&data, path)
    }

    fn parse(data: &[u8], path: &Path) -> Result<Self, GziError> {
        if data.len() < 8 {
            return Err(GziError::TooShort { path: path.to_path_buf(), actual: data.len() });
        }

        let count =
            u64::from_le_bytes(data.get(..8).and_then(|s| s.try_into().ok()).ok_or_else(|| {
                GziError::TooShort { path: path.to_path_buf(), actual: data.len() }
            })?);

        let expected_size = count
            .checked_mul(16)
            .and_then(|n| n.checked_add(8))
            .ok_or_else(|| GziError::SizeOverflow { path: path.to_path_buf(), count })?;

        if data.len() as u64 != expected_size {
            return Err(GziError::WrongSize {
                path: path.to_path_buf(),
                count,
                expected: expected_size,
                actual: data.len(),
            });
        }

        let mut entries: Vec<GziEntry> = Vec::with_capacity(count as usize);
        for i in 0..count as usize {
            let base = 8 + i * 16;
            let compressed_offset = u64::from_le_bytes(
                data.get(base..base + 8).and_then(|s| s.try_into().ok()).ok_or_else(|| {
                    GziError::TruncatedEntry { path: path.to_path_buf(), index: i }
                })?,
            );
            let uncompressed_offset = u64::from_le_bytes(
                data.get(base + 8..base + 16).and_then(|s| s.try_into().ok()).ok_or_else(|| {
                    GziError::TruncatedEntry { path: path.to_path_buf(), index: i }
                })?,
            );
            // r[impl fasta.gzi.sorted]
            if let Some(prev) = entries.last()
                && uncompressed_offset <= prev.uncompressed_offset
            {
                return Err(GziError::UnsortedEntries {
                    path: path.to_path_buf(),
                    index: i,
                    offset: uncompressed_offset,
                    prev_offset: prev.uncompressed_offset,
                });
            }

            entries.push(GziEntry { compressed_offset, uncompressed_offset });
        }

        Ok(GziIndex { entries })
    }

    #[doc(hidden)]
    pub fn parse_test(data: &[u8]) -> Self {
        Self::parse(data, Path::new("<test>")).expect("test GZI data should be valid")
    }

    /// Parse GZI index from raw bytes.
    #[cfg(feature = "fuzz")]
    pub fn from_bytes(data: &[u8]) -> Result<Self, GziError> {
        Self::parse(data, Path::new("<fuzz>"))
    }

    // r[impl fasta.gzi.translate]
    pub fn translate(&self, uncompressed_offset: u64) -> Result<BlockLocation, GziError> {
        // The GZI entries are sorted by uncompressed_offset. We want the last entry
        // whose uncompressed_offset <= target. If no entries or all entries are after
        // the target, the block starts at file position 0, uncompressed offset 0
        // (the implicit first block that GZI omits).
        let (compressed_offset, within_block) = match self
            .entries
            .binary_search_by_key(&uncompressed_offset, |e| e.uncompressed_offset)
        {
            Ok(i) => {
                let e = self.entries.get(i).ok_or(GziError::InternalBoundsError { index: i })?;
                (e.compressed_offset, 0u64)
            }
            Err(0) => (0, uncompressed_offset),
            Err(i) => {
                let e = self
                    .entries
                    .get(i - 1)
                    .ok_or(GziError::InternalBoundsError { index: i - 1 })?;
                (e.compressed_offset, uncompressed_offset - e.uncompressed_offset)
            }
        };

        // r[impl fasta.gzi.within_block_bounds]
        let within_block_offset = u16::try_from(within_block)
            .map_err(|_| GziError::WithinBlockOverflow { within_block })?;

        Ok(BlockLocation { compressed_offset, within_block_offset })
    }
}

#[derive(Debug, Clone, Copy)]
pub struct BlockLocation {
    /// File offset of the BGZF block containing the target byte.
    pub compressed_offset: u64,
    /// Offset within the decompressed block data. BGZF blocks decompress to at
    /// most 65,536 bytes, so this always fits in `u16`.
    pub within_block_offset: u16,
}

#[cfg(test)]
#[allow(clippy::arithmetic_side_effects, reason = "test arithmetic is not safety-critical")]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn test_path() -> PathBuf {
        PathBuf::from("test.gzi")
    }

    fn make_gzi_data(entries: &[(u64, u64)]) -> Vec<u8> {
        let mut data = Vec::new();
        data.extend_from_slice(&(entries.len() as u64).to_le_bytes());
        for &(compressed, uncompressed) in entries {
            data.extend_from_slice(&compressed.to_le_bytes());
            data.extend_from_slice(&uncompressed.to_le_bytes());
        }
        data
    }

    // r[verify fasta.gzi.parse]
    #[test]
    fn parse_empty_gzi() {
        let data = make_gzi_data(&[]);
        let gzi = GziIndex::parse(&data, &test_path()).unwrap();
        assert!(gzi.entries.is_empty());
    }

    // r[verify fasta.gzi.parse]
    #[test]
    fn parse_gzi_with_entries() {
        let data = make_gzi_data(&[(65536, 65280), (131072, 130560)]);
        let gzi = GziIndex::parse(&data, &test_path()).unwrap();
        assert_eq!(gzi.entries.len(), 2);
        assert_eq!(gzi.entries[0].compressed_offset, 65536);
        assert_eq!(gzi.entries[0].uncompressed_offset, 65280);
    }

    // r[verify fasta.gzi.parse]
    #[test]
    fn parse_truncated_gzi() {
        let err = GziIndex::parse(&[0; 4], &test_path()).unwrap_err();
        assert!(matches!(err, GziError::TooShort { .. }));
    }

    // r[verify fasta.gzi.parse]
    #[test]
    fn parse_wrong_size_gzi() {
        // Says 2 entries but only has data for 1
        let mut data = make_gzi_data(&[(100, 200)]);
        // Override count to 2
        data[..8].copy_from_slice(&2u64.to_le_bytes());
        let err = GziIndex::parse(&data, &test_path()).unwrap_err();
        assert!(matches!(err, GziError::WrongSize { .. }));
    }

    #[test]
    fn translate_before_first_entry() {
        let data = make_gzi_data(&[(65536, 65280)]);
        let gzi = GziIndex::parse(&data, &test_path()).unwrap();

        let loc = gzi.translate(100).unwrap();
        assert_eq!(loc.compressed_offset, 0);
        assert_eq!(loc.within_block_offset, 100);
    }

    #[test]
    fn translate_exact_match() {
        let data = make_gzi_data(&[(65536, 65280), (131072, 130560)]);
        let gzi = GziIndex::parse(&data, &test_path()).unwrap();

        let loc = gzi.translate(65280).unwrap();
        assert_eq!(loc.compressed_offset, 65536);
        assert_eq!(loc.within_block_offset, 0);
    }

    #[test]
    fn translate_between_entries() {
        let data = make_gzi_data(&[(65536, 65280), (131072, 130560)]);
        let gzi = GziIndex::parse(&data, &test_path()).unwrap();

        let loc = gzi.translate(70000).unwrap();
        assert_eq!(loc.compressed_offset, 65536);
        assert_eq!(loc.within_block_offset, (70000 - 65280) as u16);
    }

    #[test]
    fn translate_after_last_entry() {
        let data = make_gzi_data(&[(65536, 65280)]);
        let gzi = GziIndex::parse(&data, &test_path()).unwrap();

        let loc = gzi.translate(100000).unwrap();
        assert_eq!(loc.compressed_offset, 65536);
        assert_eq!(loc.within_block_offset, (100000 - 65280) as u16);
    }

    #[test]
    fn translate_empty_gzi() {
        let data = make_gzi_data(&[]);
        let gzi = GziIndex::parse(&data, &test_path()).unwrap();

        let loc = gzi.translate(12345).unwrap();
        assert_eq!(loc.compressed_offset, 0);
        assert_eq!(loc.within_block_offset, 12345);
    }

    // r[verify fasta.gzi.within_block_bounds]
    #[test]
    fn translate_rejects_within_block_exceeding_u16() {
        // Single entry at uncompressed offset 0 — any target > 65535 would
        // need a within-block offset that doesn't fit in u16
        let data = make_gzi_data(&[]);
        let gzi = GziIndex::parse(&data, &test_path()).unwrap();

        // 65535 should be fine
        let loc = gzi.translate(65535).unwrap();
        assert_eq!(loc.within_block_offset, 65535);

        // 65536 should fail — exceeds max BGZF uncompressed block size
        let err = gzi.translate(65536).unwrap_err();
        assert!(matches!(err, GziError::WithinBlockOverflow { .. }));
    }

    // r[verify fasta.gzi.sorted]
    #[test]
    fn reject_unsorted_entries() {
        // Second entry has smaller uncompressed_offset than first
        let data = make_gzi_data(&[(100, 200), (300, 100)]);
        let err = GziIndex::parse(&data, &test_path()).unwrap_err();
        assert!(
            matches!(err, GziError::UnsortedEntries { .. }),
            "expected UnsortedEntries, got {err:?}"
        );
    }

    // r[verify fasta.gzi.sorted]
    #[test]
    fn reject_duplicate_uncompressed_offsets() {
        let data = make_gzi_data(&[(100, 200), (300, 200)]);
        let err = GziIndex::parse(&data, &test_path()).unwrap_err();
        assert!(matches!(err, GziError::UnsortedEntries { .. }));
    }

    // r[verify fasta.gzi.parse]
    #[test]
    fn reject_overflow_count() {
        // Craft a GZI header claiming u64::MAX entries
        let mut data = vec![0u8; 8];
        data[..8].copy_from_slice(&u64::MAX.to_le_bytes());
        let err = GziIndex::parse(&data, &test_path()).unwrap_err();
        assert!(matches!(err, GziError::SizeOverflow { .. }), "expected SizeOverflow, got {err:?}");
    }
}
