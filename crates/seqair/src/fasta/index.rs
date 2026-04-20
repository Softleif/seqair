//! Parse FAI indexes to locate sequences within a FASTA file. [`FastaIndex`] maps sequence
//! names to [`FaiEntry`] offsets and line widths; used by [`IndexedFastaReader`](super::IndexedFastaReader)
//! to seek directly to any subsequence.

use rustc_hash::FxHashMap;
use seqair_types::SmolStr;
use std::path::Path;
use tracing::instrument;

// r[impl io.errors]
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum FaiError {
    #[error("I/O error reading FAI index {path}")]
    Read { path: std::path::PathBuf, source: std::io::Error },

    #[error("FAI index {path} is empty")]
    Empty { path: std::path::PathBuf },

    #[error("FAI index {path} line {line_number}: {kind}")]
    InvalidEntry { path: std::path::PathBuf, line_number: usize, kind: FaiEntryError },
}

#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum FaiEntryError {
    #[error("expected 5 tab-separated fields, found {found}")]
    TooFewFields { found: usize },

    #[error("expected 5 tab-separated fields, found {found}")]
    TooManyFields { found: usize },

    #[error("invalid {field} value {raw_value:?}: not a valid integer")]
    InvalidField { field: &'static str, raw_value: String },

    #[error("linebases is 0 (would cause division by zero) for sequence {name:?}")]
    ZeroLinebases { name: SmolStr },

    #[error("linewidth ({linewidth}) < linebases ({linebases}) for sequence {name:?}")]
    LinewidthTooSmall { name: SmolStr, linebases: u64, linewidth: u64 },

    #[error("sequence length is 0 for sequence {name:?}")]
    ZeroLength { name: SmolStr },

    #[error("duplicate sequence name {name:?}")]
    DuplicateName { name: SmolStr },
}

/// A single entry from a `.fai` index file.
#[derive(Debug, Clone)]
pub struct FaiEntry {
    pub name: SmolStr,
    pub length: u64,
    pub offset: u64,
    pub linebases: u64,
    pub linewidth: u64,
}

impl FaiEntry {
    // r[impl fasta.index.offset_calculation]
    pub fn byte_offset(&self, pos: u64) -> u64 {
        let line =
            pos.checked_div(self.linebases).expect("linebases is non-zero, enforced by FAI parser");
        let col =
            pos.checked_rem(self.linebases).expect("linebases is non-zero, enforced by FAI parser");
        self.offset.wrapping_add(line.wrapping_mul(self.linewidth)).wrapping_add(col)
    }
}

// r[impl fasta.index.parse]
// r[impl fasta.index.name_lookup]
#[derive(Debug, Clone)]
pub struct FastaIndex {
    entries: Vec<FaiEntry>,
    name_to_idx: FxHashMap<SmolStr, usize>,
}

impl FastaIndex {
    #[instrument(level = "debug", fields(path = %path.display()))]
    pub fn from_file(path: &Path) -> Result<Self, FaiError> {
        let contents = std::fs::read_to_string(path)
            .map_err(|source| FaiError::Read { path: path.to_path_buf(), source })?;
        Self::parse(&contents, path)
    }

    /// Parse FAI index from string contents.
    #[cfg(feature = "fuzz")]
    pub fn from_contents(contents: &str) -> Result<Self, FaiError> {
        Self::parse(contents, Path::new("<fuzz>"))
    }

    fn parse(contents: &str, path: &Path) -> Result<Self, FaiError> {
        let mut entries = Vec::new();
        let mut name_to_idx = FxHashMap::default();

        for (line_num, line) in contents.lines().enumerate() {
            if line.is_empty() {
                continue;
            }

            // r[impl fasta.index.fields]
            let mut fields = line.splitn(6, '\t');
            let (Some(name_str), Some(len_s), Some(off_s), Some(lb_s), Some(lw_s)) =
                (fields.next(), fields.next(), fields.next(), fields.next(), fields.next())
            else {
                let found = line.split('\t').count();
                return Err(FaiError::InvalidEntry {
                    path: path.to_path_buf(),
                    line_number: line_num.wrapping_add(1),
                    kind: FaiEntryError::TooFewFields { found },
                });
            };
            if fields.next().is_some() {
                let found = line.split('\t').count();
                return Err(FaiError::InvalidEntry {
                    path: path.to_path_buf(),
                    line_number: line_num.wrapping_add(1),
                    kind: FaiEntryError::TooManyFields { found },
                });
            }

            let name = SmolStr::new(name_str);
            let length = parse_fai_field(len_s, "length", line_num, path)?;
            let offset = parse_fai_field(off_s, "offset", line_num, path)?;
            let linebases = parse_fai_field(lb_s, "linebases", line_num, path)?;
            let linewidth = parse_fai_field(lw_s, "linewidth", line_num, path)?;

            // r[impl fasta.index.validation]
            if linebases == 0 {
                return Err(FaiError::InvalidEntry {
                    path: path.to_path_buf(),
                    line_number: line_num.wrapping_add(1),
                    kind: FaiEntryError::ZeroLinebases { name },
                });
            }
            if linewidth < linebases {
                return Err(FaiError::InvalidEntry {
                    path: path.to_path_buf(),
                    line_number: line_num.wrapping_add(1),
                    kind: FaiEntryError::LinewidthTooSmall { name, linebases, linewidth },
                });
            }
            if length == 0 {
                return Err(FaiError::InvalidEntry {
                    path: path.to_path_buf(),
                    line_number: line_num.wrapping_add(1),
                    kind: FaiEntryError::ZeroLength { name },
                });
            }

            if name_to_idx.contains_key(&name) {
                return Err(FaiError::InvalidEntry {
                    path: path.to_path_buf(),
                    line_number: line_num.wrapping_add(1),
                    kind: FaiEntryError::DuplicateName { name },
                });
            }

            name_to_idx.insert(name.clone(), entries.len());
            entries.push(FaiEntry { name, length, offset, linebases, linewidth });
        }

        if entries.is_empty() {
            return Err(FaiError::Empty { path: path.to_path_buf() });
        }

        Ok(FastaIndex { entries, name_to_idx })
    }

    pub fn get(&self, name: &str) -> Option<&FaiEntry> {
        let &idx = self.name_to_idx.get(name)?;
        self.entries.get(idx)
    }

    pub fn sequence_names(&self) -> Vec<SmolStr> {
        self.entries.iter().map(|e| e.name.clone()).collect()
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

fn parse_fai_field(
    s: &str,
    field_name: &'static str,
    line_num: usize,
    path: &Path,
) -> Result<u64, FaiError> {
    s.parse().map_err(|_| FaiError::InvalidEntry {
        path: path.to_path_buf(),
        line_number: line_num.wrapping_add(1),
        kind: FaiEntryError::InvalidField { field: field_name, raw_value: s.to_string() },
    })
}

#[cfg(test)]
#[allow(clippy::arithmetic_side_effects, reason = "test arithmetic is not safety-critical")]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn test_path() -> PathBuf {
        PathBuf::from("test.fai")
    }

    // r[verify fasta.index.parse]
    // r[verify fasta.index.name_lookup]
    #[test]
    fn parse_simple_index() {
        let contents = "seq1\t100\t6\t50\t51\nseq2\t200\t160\t80\t81\n";
        let idx = FastaIndex::parse(contents, &test_path()).unwrap();

        assert_eq!(idx.len(), 2);
        let e1 = idx.get("seq1").unwrap();
        assert_eq!(e1.length, 100);
        assert_eq!(e1.offset, 6);
        assert_eq!(e1.linebases, 50);
        assert_eq!(e1.linewidth, 51);

        let e2 = idx.get("seq2").unwrap();
        assert_eq!(e2.length, 200);
        assert_eq!(e2.offset, 160);
    }

    // r[verify fasta.index.name_lookup]
    #[test]
    fn lookup_missing_sequence() {
        let contents = "seq1\t100\t6\t50\t51\n";
        let idx = FastaIndex::parse(contents, &test_path()).unwrap();
        assert!(idx.get("nonexistent").is_none());
    }

    // r[verify fasta.index.fields]
    #[test]
    fn reject_too_few_fields() {
        let contents = "seq1\t100\t6\t50\n"; // only 4 fields
        let err = FastaIndex::parse(contents, &test_path()).unwrap_err();
        match err {
            FaiError::InvalidEntry {
                kind: FaiEntryError::TooFewFields { found },
                line_number,
                ..
            } => {
                assert_eq!(found, 4, "should report 4 fields found");
                assert_eq!(line_number, 1);
            }
            other => panic!("expected TooFewFields, got {other:?}"),
        }
    }

    // r[verify fasta.index.fields]
    #[test]
    fn reject_too_many_fields() {
        let contents = "seq1\t100\t6\t50\t51\textra\n"; // 6 fields
        let err = FastaIndex::parse(contents, &test_path()).unwrap_err();
        match err {
            FaiError::InvalidEntry {
                kind: FaiEntryError::TooManyFields { found },
                line_number,
                ..
            } => {
                assert_eq!(found, 6, "should report 6 fields found");
                assert_eq!(line_number, 1);
            }
            other => panic!("expected TooManyFields, got {other:?}"),
        }
    }

    // r[verify fasta.index.fields]
    #[test]
    fn reject_invalid_field_value() {
        let contents = "chr1\tNOTANUM\t0\t50\t51\n"; // non-integer length
        let err = FastaIndex::parse(contents, &test_path()).unwrap_err();
        match err {
            FaiError::InvalidEntry {
                kind: FaiEntryError::InvalidField { field, raw_value },
                line_number,
                ..
            } => {
                assert_eq!(field, "length");
                assert_eq!(raw_value, "NOTANUM");
                assert_eq!(line_number, 1);
            }
            other => panic!("expected InvalidField, got {other:?}"),
        }
    }

    // r[verify fasta.index.fields]
    #[test]
    fn reject_invalid_field_via_file() {
        use std::io::Write;
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("bad.fai");
        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(f, "chr1\t1000\tBAD_OFFSET\t50\t51").unwrap();
        drop(f);

        let err = FastaIndex::from_file(&path).unwrap_err();
        match err {
            FaiError::InvalidEntry {
                kind: FaiEntryError::InvalidField { field, raw_value },
                ..
            } => {
                assert_eq!(field, "offset");
                assert_eq!(raw_value, "BAD_OFFSET");
            }
            other => panic!("expected InvalidField, got {other:?}"),
        }
    }

    // r[verify fasta.index.validation]
    #[test]
    fn reject_zero_linebases() {
        let contents = "seq1\t100\t6\t0\t1\n";
        let err = FastaIndex::parse(contents, &test_path()).unwrap_err();
        assert!(
            matches!(err, FaiError::InvalidEntry { kind: FaiEntryError::ZeroLinebases { .. }, .. }),
            "expected ZeroLinebases, got {err:?}"
        );
    }

    // r[verify fasta.index.validation]
    #[test]
    fn reject_linewidth_less_than_linebases() {
        let contents = "seq1\t100\t6\t50\t49\n";
        let err = FastaIndex::parse(contents, &test_path()).unwrap_err();
        assert!(
            matches!(
                err,
                FaiError::InvalidEntry { kind: FaiEntryError::LinewidthTooSmall { .. }, .. }
            ),
            "expected LinewidthTooSmall, got {err:?}"
        );
    }

    // r[verify fasta.index.validation]
    #[test]
    fn reject_zero_length() {
        let contents = "seq1\t0\t6\t50\t51\n";
        let err = FastaIndex::parse(contents, &test_path()).unwrap_err();
        assert!(
            matches!(err, FaiError::InvalidEntry { kind: FaiEntryError::ZeroLength { .. }, .. }),
            "expected ZeroLength, got {err:?}"
        );
    }

    // r[verify fasta.index.name_lookup]
    #[test]
    fn reject_duplicate_names() {
        let contents = "seq1\t100\t6\t50\t51\nseq1\t200\t160\t80\t81\n";
        let err = FastaIndex::parse(contents, &test_path()).unwrap_err();
        assert!(
            matches!(err, FaiError::InvalidEntry { kind: FaiEntryError::DuplicateName { .. }, .. }),
            "expected DuplicateName, got {err:?}"
        );
    }

    #[test]
    fn reject_empty_index() {
        let err = FastaIndex::parse("", &test_path()).unwrap_err();
        assert!(matches!(err, FaiError::Empty { .. }), "expected Empty, got {err:?}");
    }

    #[test]
    fn byte_offset_single_line() {
        // Sequence fits on one line: 4 bases, linebases=4, linewidth=5
        let entry = FaiEntry { name: "s".into(), length: 4, offset: 6, linebases: 4, linewidth: 5 };
        // pos 0 → offset + 0 = 6
        assert_eq!(entry.byte_offset(0), 6);
        // pos 3 → offset + 3 = 9
        assert_eq!(entry.byte_offset(3), 9);
    }

    #[test]
    fn byte_offset_multi_line() {
        // 100 bases, 50 per line, linewidth 51 (50 bases + \n)
        let entry =
            FaiEntry { name: "s".into(), length: 100, offset: 10, linebases: 50, linewidth: 51 };
        // pos 0 → 10
        assert_eq!(entry.byte_offset(0), 10);
        // pos 49 → 10 + 49 = 59 (last base of first line)
        assert_eq!(entry.byte_offset(49), 59);
        // pos 50 → 10 + 1*51 + 0 = 61 (first base of second line, skipping \n)
        assert_eq!(entry.byte_offset(50), 61);
        // pos 51 → 10 + 1*51 + 1 = 62
        assert_eq!(entry.byte_offset(51), 62);
    }

    #[test]
    fn byte_offset_crlf() {
        // Windows-style line endings: linewidth = linebases + 2
        let entry = FaiEntry {
            name: "s".into(),
            length: 100,
            offset: 10,
            linebases: 50,
            linewidth: 52, // 50 bases + \r\n
        };
        // pos 50 → 10 + 1*52 + 0 = 62 (skips \r\n)
        assert_eq!(entry.byte_offset(50), 62);
    }

    #[test]
    fn parse_real_fai() {
        let contents = "chr19\t61431566\t7\t50\t51\n\
                         2kb_3_Unmodified\t2018\t62660223\t2018\t2019\n\
                         bacteriophage_lambda_CpG\t48502\t62662268\t48502\t48503\n";
        let idx = FastaIndex::parse(contents, &test_path()).unwrap();
        assert_eq!(idx.len(), 3);

        let chr19 = idx.get("chr19").unwrap();
        assert_eq!(chr19.length, 61_431_566);
        assert_eq!(chr19.offset, 7);

        let lambda = idx.get("bacteriophage_lambda_CpG").unwrap();
        assert_eq!(lambda.length, 48502);
    }
}
