//! Fetch reference sequences by name and coordinate range. [`IndexedFastaReader`] uses a
//! [`FastaIndex`](super::FastaIndex) for seeking and an optional [`GziIndex`](super::GziIndex)
//! for bgzf-compressed references. Returns raw bytes; call [`Base::convert_ascii_in_place`](seqair_types::Base::convert_ascii_in_place) to convert.

use super::{
    gzi::{GziError, GziIndex},
    index::{FaiError, FastaIndex},
};
use crate::bam::bgzf::{BgzfError, BgzfReader};
use seqair_types::{Pos0, SmolStr};
use std::{
    fs::File,
    io::{Read, Seek, SeekFrom},
    path::{Path, PathBuf},
    sync::Arc,
};
use tracing::instrument;

#[cfg(feature = "fuzz")]
use std::io::Cursor;

// r[impl fasta.errors]
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum FastaError {
    #[error(
        "FAI index not found at {expected_path}\n\nCreate it with: samtools faidx {fasta_path}"
    )]
    FaiIndexNotFound { fasta_path: PathBuf, expected_path: PathBuf },

    #[error("FAI index error")]
    Fai {
        #[from]
        source: FaiError,
    },

    #[error(
        "GZI index not found at {expected_path}\n\nCreate it with: samtools faidx {fasta_path}"
    )]
    GziIndexNotFound { fasta_path: PathBuf, expected_path: PathBuf },

    #[error("GZI index error")]
    Gzi {
        #[from]
        source: GziError,
    },

    #[error("I/O error reading FASTA {path}")]
    Read { path: PathBuf, source: std::io::Error },

    #[error("BGZF error reading FASTA")]
    Bgzf {
        #[from]
        source: BgzfError,
    },

    #[error("sequence not found: {name:?}{}", format_available_sequences(.available))]
    SequenceNotFound { name: SmolStr, available: Vec<SmolStr> },

    #[error("FASTA region {name:?}:{start}-{end} is out of bounds (sequence length: {seq_len})")]
    RegionOutOfBounds { name: SmolStr, start: u64, end: u64, seq_len: u64 },

    #[error("BGZF FASTA missing GZI index (internal error)")]
    MissingGzi,

    #[error("no FASTA reference attached; reopen with `Readers::open` to enable reference lookups")]
    NotConfigured,
}

fn format_available_sequences(available: &[SmolStr]) -> String {
    if available.is_empty() || available.len() > 20 {
        String::new()
    } else {
        format!("\n\nAvailable sequences: {}", available.join(", "))
    }
}

// r[impl fasta.fork.shared_state]
struct FastaShared {
    index: FastaIndex,
    gzi: Option<GziIndex>,
    fasta_path: PathBuf,
    is_bgzf: bool,
}

enum FileHandle<R: Read + Seek> {
    Plain(R),
    Bgzf(BgzfReader<R>),
}

/// Indexed FASTA reader for random-access subsequence fetching.
///
/// Supports plain FASTA (`.fa` / `.fasta`) and BGZF-compressed FASTA
/// (`.fa.gz`). Both require an `.fai` index; BGZF additionally requires a
/// `.gzi` index. Create both with `samtools faidx`.
///
/// # Example
///
/// ```no_run
/// use seqair::fasta::IndexedFastaReader;
/// use seqair_types::{Base, Pos0};
/// use std::path::Path;
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let mut reader = IndexedFastaReader::open(Path::new("reference.fa"))?;
///
/// // Fetch raw ASCII bytes for a region (half-open: start inclusive, end exclusive)
/// let seq = reader.fetch_seq("chr1", Pos0::new(1000).unwrap(), Pos0::new(2000).unwrap())?;
///
/// // Convert to typed Base values (A/C/G/T/Unknown) for downstream use
/// let bases = Base::from_ascii_vec(seq);
/// assert_eq!(bases.len(), 1000);
/// # Ok(())
/// # }
/// ```
// r[impl fasta.fetch.coordinates]
pub struct IndexedFastaReader<R: Read + Seek = File> {
    handle: FileHandle<R>,
    shared: Arc<FastaShared>,
    /// Reusable buffer for raw bytes read from plain FASTA (includes newlines).
    raw_buf: Vec<u8>,
}

impl<R: Read + Seek> std::fmt::Debug for IndexedFastaReader<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("IndexedFastaReader")
            .field("fasta_path", &self.shared.fasta_path)
            .field("is_bgzf", &self.shared.is_bgzf)
            .finish()
    }
}

impl IndexedFastaReader<File> {
    #[instrument(level = "debug", fields(path = %path.display()))]
    pub fn open(path: &Path) -> Result<Self, FastaError> {
        let fai_path = fai_path_for(path);
        // r[impl fasta.index.missing]
        if !fai_path.exists() {
            return Err(FastaError::FaiIndexNotFound {
                fasta_path: path.to_path_buf(),
                expected_path: fai_path,
            });
        }
        let index = FastaIndex::from_file(&fai_path)?;

        // r[impl fasta.bgzf.detect]
        let is_bgzf = detect_bgzf(path)?;

        // r[impl fasta.bgzf.gzi_required]
        let gzi = if is_bgzf {
            let gzi_path = gzi_path_for(path);
            // r[impl fasta.bgzf.gzi_missing]
            if !gzi_path.exists() {
                return Err(FastaError::GziIndexNotFound {
                    fasta_path: path.to_path_buf(),
                    expected_path: gzi_path,
                });
            }
            Some(GziIndex::from_file(&gzi_path)?)
        } else {
            None
        };

        let handle = if is_bgzf {
            FileHandle::Bgzf(BgzfReader::open(path)?)
        } else {
            let file = File::open(path)
                .map_err(|source| FastaError::Read { path: path.to_path_buf(), source })?;
            FileHandle::Plain(file)
        };

        Ok(IndexedFastaReader {
            handle,
            shared: Arc::new(FastaShared { index, gzi, fasta_path: path.to_path_buf(), is_bgzf }),
            raw_buf: Vec::with_capacity(64 * 1024),
        })
    }

    // r[impl fasta.fork.operation]
    // r[impl fasta.fork.independence]
    #[instrument(level = "debug", skip(self), fields(path = %self.shared.fasta_path.display()))]
    pub fn fork(&self) -> Result<Self, FastaError> {
        let handle = if self.shared.is_bgzf {
            FileHandle::Bgzf(BgzfReader::open(&self.shared.fasta_path)?)
        } else {
            let file = File::open(&self.shared.fasta_path).map_err(|source| FastaError::Read {
                path: self.shared.fasta_path.clone(),
                source,
            })?;
            FileHandle::Plain(file)
        };

        Ok(IndexedFastaReader {
            handle,
            shared: Arc::clone(&self.shared),
            raw_buf: Vec::with_capacity(64 * 1024),
        })
    }
}

#[cfg(feature = "fuzz")]
impl IndexedFastaReader<Cursor<Vec<u8>>> {
    pub fn from_bytes(
        fasta_gz_data: Vec<u8>,
        fai_contents: &str,
        gzi_data: &[u8],
    ) -> Result<Self, FastaError> {
        let index = FastaIndex::from_contents(fai_contents)?;
        let gzi = GziIndex::from_bytes(gzi_data)?;
        let bgzf = BgzfReader::from_cursor(fasta_gz_data);
        Ok(IndexedFastaReader {
            handle: FileHandle::Bgzf(bgzf),
            shared: Arc::new(FastaShared {
                index,
                gzi: Some(gzi),
                fasta_path: PathBuf::from("<fuzz>"),
                is_bgzf: true,
            }),
            raw_buf: Vec::with_capacity(64 * 1024),
        })
    }
}

impl<R: Read + Seek> IndexedFastaReader<R> {
    pub fn shares_index_with(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.shared, &other.shared)
    }

    pub fn index(&self) -> &FastaIndex {
        &self.shared.index
    }

    // r[impl fasta.fetch.raw_bytes]
    // r[impl fasta.fetch.buffer_reuse]
    #[instrument(level = "debug", skip(self), fields(name, start, stop))]
    pub fn fetch_seq(
        &mut self,
        name: &str,
        start: Pos0,
        stop: Pos0,
    ) -> Result<Vec<u8>, FastaError> {
        let mut buf = Vec::new();
        self.fetch_seq_into(name, start, stop, &mut buf)?;
        Ok(buf)
    }

    // r[impl fasta.fetch.raw_bytes]
    #[instrument(level = "debug", skip(self, out), fields(name, start, stop))]
    pub fn fetch_seq_into(
        &mut self,
        name: &str,
        start: Pos0,
        stop: Pos0,
        out: &mut Vec<u8>,
    ) -> Result<(), FastaError> {
        self.fetch_seq_into_u64(name, start.as_u64(), stop.as_u64(), out)
    }

    /// Same as [`fetch_seq_into`] but accepts raw `u64` bounds so callers
    /// can express `end + 1 == i32::MAX + 1` (one past the last representable
    /// `Pos0`). Used by the pileup path to handle a segment whose inclusive
    /// `end` is `Pos0::max_value()` without losing the last base.
    ///
    /// [`fetch_seq_into`]: Self::fetch_seq_into
    pub fn fetch_seq_into_u64(
        &mut self,
        name: &str,
        start: u64,
        stop: u64,
        out: &mut Vec<u8>,
    ) -> Result<(), FastaError> {
        out.clear();

        // r[impl fasta.fetch.unknown_sequence]
        let entry = self
            .shared
            .index
            .get(name)
            .ok_or_else(|| FastaError::SequenceNotFound {
                name: SmolStr::new(name),
                available: self.shared.index.sequence_names(),
            })?
            .clone();

        // r[impl fasta.fetch.bounds_check]
        if start >= stop {
            return Err(FastaError::RegionOutOfBounds {
                name: SmolStr::new(name),
                start,
                end: stop,
                seq_len: entry.length,
            });
        }
        if stop > entry.length {
            return Err(FastaError::RegionOutOfBounds {
                name: SmolStr::new(name),
                start,
                end: stop,
                seq_len: entry.length,
            });
        }

        #[expect(
            clippy::cast_possible_truncation,
            reason = "FASTA sequences fit in usize on supported platforms; 32-bit addressable memory bounds the read size"
        )]
        let num_bases =
            stop.checked_sub(start).expect("stop > start is guaranteed by bounds check above")
                as usize;
        let start_byte = entry.byte_offset(start);
        let last_base_byte = entry.byte_offset(
            (num_bases as u64)
                .checked_sub(1)
                .and_then(|n| start.checked_add(n))
                .expect("num_bases >= 1 is guaranteed by stop > start"),
        );
        let end_byte = last_base_byte
            .checked_add(1)
            .expect("byte_offset is within the file, so adding 1 cannot overflow u64");
        #[expect(
            clippy::cast_possible_truncation,
            reason = "raw_len is bounded by the FASTA region size which fits in usize on supported platforms"
        )]
        let raw_len = end_byte
            .checked_sub(start_byte)
            .expect("end_byte >= start_byte for any valid FAI entry")
            as usize;

        // r[impl fasta.plain.read]
        // r[impl fasta.plain.read_size]
        // r[impl fasta.bgzf.sequential_read]
        // r[impl fasta.bgzf.decompress]
        self.raw_buf.clear();
        self.raw_buf.resize(raw_len, 0);

        match &mut self.handle {
            FileHandle::Plain(file) => {
                file.seek(SeekFrom::Start(start_byte)).map_err(|source| FastaError::Read {
                    path: self.shared.fasta_path.clone(),
                    source,
                })?;
                file.read_exact(&mut self.raw_buf).map_err(|source| FastaError::Read {
                    path: self.shared.fasta_path.clone(),
                    source,
                })?;
            }
            FileHandle::Bgzf(bgzf) => {
                let gzi = self.shared.gzi.as_ref().ok_or(FastaError::MissingGzi)?;

                let block_loc = gzi.translate(start_byte)?;
                let voff = crate::bam::bgzf::VirtualOffset::new(
                    block_loc.compressed_offset,
                    block_loc.within_block_offset,
                );

                bgzf.seek_virtual(voff)?;
                bgzf.read_exact_into(&mut self.raw_buf)?;
            }
        }

        // r[impl fasta.fetch.newline_stripping]
        // r[impl fasta.fetch.uppercase]
        out.reserve(num_bases);
        for &b in &self.raw_buf {
            if b != b'\n' && b != b'\r' {
                out.push(b.to_ascii_uppercase());
            }
        }

        Ok(())
    }
}

// r[impl fasta.index.location]
fn fai_path_for(fasta_path: &Path) -> PathBuf {
    let mut p = fasta_path.as_os_str().to_owned();
    p.push(".fai");
    PathBuf::from(p)
}

fn gzi_path_for(fasta_path: &Path) -> PathBuf {
    let mut p = fasta_path.as_os_str().to_owned();
    p.push(".gzi");
    PathBuf::from(p)
}

// r[impl fasta.bgzf.detect]
fn detect_bgzf(path: &Path) -> Result<bool, FastaError> {
    let mut file =
        File::open(path).map_err(|source| FastaError::Read { path: path.to_path_buf(), source })?;

    // Full BGZF header is 18 bytes: 10-byte gzip header + 2 XLEN + 6 BC subfield.
    // Check gzip magic (1f 8b), DEFLATE method (08), FEXTRA flag (04),
    // then verify XLEN and BC subfield presence.
    let mut header = [0u8; 18];
    match file.read_exact(&mut header) {
        Ok(()) => {}
        Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => return Ok(false),
        Err(source) => return Err(FastaError::Read { path: path.to_path_buf(), source }),
    }

    // Bytes 0-3: gzip magic + method + FEXTRA flag
    if header.get(..4) != Some(&[0x1f, 0x8b, 0x08, 0x04]) {
        return Ok(false);
    }

    // Bytes 12-13: BC subfield identifier
    Ok(header.get(12..14) == Some(b"BC"))
}

#[cfg(test)]
#[allow(clippy::arithmetic_side_effects, reason = "test arithmetic is not safety-critical")]
mod tests {
    use super::*;
    use seqair_types::Pos0;
    use std::io::Write;
    use tempfile::TempDir;

    fn p(v: u32) -> Pos0 {
        Pos0::new(v).unwrap()
    }

    fn make_plain_fasta(dir: &TempDir) -> (PathBuf, PathBuf) {
        let fasta_path = dir.path().join("test.fa");
        let fai_path = dir.path().join("test.fa.fai");

        let mut f = File::create(&fasta_path).unwrap();
        // Two sequences: seq1 (8 bases, 4 per line) and seq2 (4 bases, single line)
        f.write_all(b">seq1\nACGT\nTGCA\n>seq2\nGGCC\n").unwrap();

        let mut idx = File::create(&fai_path).unwrap();
        // >seq1\n = 6 bytes, ACGT\n = 5, TGCA\n = 5, >seq2\n = 6
        // seq1: 8 bases, offset=6, 4 bases/line, 5 bytes/line
        // seq2: 4 bases, offset=22, 4 bases/line, 5 bytes/line
        idx.write_all(b"seq1\t8\t6\t4\t5\nseq2\t4\t22\t4\t5\n").unwrap();

        (fasta_path, fai_path)
    }

    // r[verify fasta.index.location]
    #[test]
    fn open_plain_fasta() {
        let dir = TempDir::new().unwrap();
        let (fasta_path, _) = make_plain_fasta(&dir);

        let reader = IndexedFastaReader::open(&fasta_path).unwrap();
        assert!(!reader.shared.is_bgzf);
        assert_eq!(reader.index().len(), 2);
    }

    // r[verify fasta.fetch.raw_bytes]
    // r[verify fasta.plain.read_size]
    #[test]
    fn fetch_full_sequence() {
        let dir = TempDir::new().unwrap();
        let (fasta_path, _) = make_plain_fasta(&dir);

        let mut reader = IndexedFastaReader::open(&fasta_path).unwrap();
        let seq = reader.fetch_seq("seq1", p(0), p(8)).unwrap();
        assert_eq!(seq, b"ACGTTGCA");
    }

    #[test]
    fn fetch_subsequence() {
        let dir = TempDir::new().unwrap();
        let (fasta_path, _) = make_plain_fasta(&dir);

        let mut reader = IndexedFastaReader::open(&fasta_path).unwrap();
        // Fetch positions 2-6 (0-based half-open), crossing the line boundary
        let seq = reader.fetch_seq("seq1", p(2), p(6)).unwrap();
        assert_eq!(seq, b"GTTG");
    }

    #[test]
    fn fetch_second_sequence() {
        let dir = TempDir::new().unwrap();
        let (fasta_path, _) = make_plain_fasta(&dir);

        let mut reader = IndexedFastaReader::open(&fasta_path).unwrap();
        let seq = reader.fetch_seq("seq2", p(0), p(4)).unwrap();
        assert_eq!(seq, b"GGCC");
    }

    #[test]
    fn fetch_single_base() {
        let dir = TempDir::new().unwrap();
        let (fasta_path, _) = make_plain_fasta(&dir);

        let mut reader = IndexedFastaReader::open(&fasta_path).unwrap();
        let seq = reader.fetch_seq("seq1", p(0), p(1)).unwrap();
        assert_eq!(seq, b"A");
    }

    #[test]
    fn fetch_uppercases_lowercase() {
        let dir = TempDir::new().unwrap();
        let fasta_path = dir.path().join("lower.fa");
        let fai_path = dir.path().join("lower.fa.fai");

        let mut f = File::create(&fasta_path).unwrap();
        f.write_all(b">seq1\nacgtNn\n").unwrap();
        let mut idx = File::create(&fai_path).unwrap();
        idx.write_all(b"seq1\t6\t6\t6\t7\n").unwrap();

        let mut reader = IndexedFastaReader::open(&fasta_path).unwrap();
        let seq = reader.fetch_seq("seq1", p(0), p(6)).unwrap();
        assert_eq!(seq, b"ACGTNN");
    }

    // r[verify fasta.fetch.bounds_check]
    // r[verify fasta.errors]
    #[test]
    fn fetch_out_of_bounds() {
        let dir = TempDir::new().unwrap();
        let (fasta_path, _) = make_plain_fasta(&dir);

        let mut reader = IndexedFastaReader::open(&fasta_path).unwrap();
        let err = reader.fetch_seq("seq1", p(0), p(100)).unwrap_err();
        assert!(matches!(err, FastaError::RegionOutOfBounds { .. }));
    }

    // r[verify fasta.fetch.unknown_sequence]
    // r[verify fasta.errors]
    #[test]
    fn fetch_unknown_sequence() {
        let dir = TempDir::new().unwrap();
        let (fasta_path, _) = make_plain_fasta(&dir);

        let mut reader = IndexedFastaReader::open(&fasta_path).unwrap();
        let err = reader.fetch_seq("nonexistent", p(0), p(4)).unwrap_err();
        match &err {
            FastaError::SequenceNotFound { name, available } => {
                assert_eq!(name, "nonexistent");
                assert!(available.contains(&SmolStr::new("seq1")));
                assert!(available.contains(&SmolStr::new("seq2")));
            }
            _ => panic!("expected SequenceNotFound, got {err:?}"),
        }
    }

    // r[verify fasta.index.missing]
    // r[verify fasta.index.location]
    // r[verify fasta.errors]
    #[test]
    fn missing_fai_error() {
        let dir = TempDir::new().unwrap();
        let fasta_path = dir.path().join("noidx.fa");
        File::create(&fasta_path).unwrap();

        let err = IndexedFastaReader::open(&fasta_path).unwrap_err();
        match &err {
            FastaError::FaiIndexNotFound { fasta_path: fp, expected_path } => {
                assert_eq!(fp, &fasta_path);
                assert!(expected_path.to_str().unwrap().ends_with(".fai"));
            }
            _ => panic!("expected FaiIndexNotFound, got {err:?}"),
        }
    }

    // r[impl fasta.fork.equivalence]
    // r[verify fasta.fork.operation]
    // r[verify fasta.fork.shared_state]
    #[test]
    fn fork_produces_identical_results() {
        let dir = TempDir::new().unwrap();
        let (fasta_path, _) = make_plain_fasta(&dir);

        let mut reader = IndexedFastaReader::open(&fasta_path).unwrap();
        let mut forked = reader.fork().unwrap();

        let seq_orig = reader.fetch_seq("seq1", p(0), p(8)).unwrap();
        let seq_fork = forked.fetch_seq("seq1", p(0), p(8)).unwrap();
        assert_eq!(seq_orig, seq_fork);
    }

    // r[verify fasta.fork.shared_state]
    #[test]
    fn fork_shares_arc() {
        let dir = TempDir::new().unwrap();
        let (fasta_path, _) = make_plain_fasta(&dir);

        let reader = IndexedFastaReader::open(&fasta_path).unwrap();
        let forked = reader.fork().unwrap();

        assert!(reader.shares_index_with(&forked));
    }

    // r[verify fasta.fork.independence]
    #[test]
    fn fork_independence() {
        let dir = TempDir::new().unwrap();
        let (fasta_path, _) = make_plain_fasta(&dir);

        let mut reader = IndexedFastaReader::open(&fasta_path).unwrap();
        let mut forked = reader.fork().unwrap();

        // Interleaved reads on different sequences
        let s1 = reader.fetch_seq("seq1", p(0), p(4)).unwrap();
        let s2 = forked.fetch_seq("seq2", p(0), p(4)).unwrap();
        let s3 = reader.fetch_seq("seq2", p(0), p(4)).unwrap();
        let s4 = forked.fetch_seq("seq1", p(0), p(4)).unwrap();

        assert_eq!(s1, b"ACGT");
        assert_eq!(s2, b"GGCC");
        assert_eq!(s3, b"GGCC");
        assert_eq!(s4, b"ACGT");
    }

    #[test]
    fn fetch_into_reuses_buffer() {
        let dir = TempDir::new().unwrap();
        let (fasta_path, _) = make_plain_fasta(&dir);

        let mut reader = IndexedFastaReader::open(&fasta_path).unwrap();
        let mut buf = Vec::new();

        reader.fetch_seq_into("seq1", p(0), p(4), &mut buf).unwrap();
        assert_eq!(buf, b"ACGT");

        // Reuse the same buffer
        reader.fetch_seq_into("seq2", p(0), p(4), &mut buf).unwrap();
        assert_eq!(buf, b"GGCC");
    }

    // r[verify fasta.plain.read_size]
    #[test]
    fn multiline_cross_boundary() {
        let dir = TempDir::new().unwrap();
        let fasta_path = dir.path().join("multi.fa");
        let fai_path = dir.path().join("multi.fa.fai");

        let mut f = File::create(&fasta_path).unwrap();
        // 3 bases per line, 10 bases total: ACG\nTAC\nGTA\nC\n
        f.write_all(b">seq1\nACG\nTAC\nGTA\nC\n").unwrap();
        let mut idx = File::create(&fai_path).unwrap();
        idx.write_all(b"seq1\t10\t6\t3\t4\n").unwrap();

        let mut reader = IndexedFastaReader::open(&fasta_path).unwrap();

        // Full sequence
        assert_eq!(reader.fetch_seq("seq1", p(0), p(10)).unwrap(), b"ACGTACGTAC");

        // Cross first boundary: pos 2-4
        assert_eq!(reader.fetch_seq("seq1", p(2), p(5)).unwrap(), b"GTA");

        // Last base
        assert_eq!(reader.fetch_seq("seq1", p(9), p(10)).unwrap(), b"C");
    }

    // r[verify fasta.fetch.bounds_check]
    #[test]
    fn fetch_empty_range_errors() {
        let dir = TempDir::new().unwrap();
        let (fasta_path, _) = make_plain_fasta(&dir);

        let mut reader = IndexedFastaReader::open(&fasta_path).unwrap();
        // start == stop → empty range, should error
        let err = reader.fetch_seq("seq1", p(5), p(5)).unwrap_err();
        assert!(matches!(err, FastaError::RegionOutOfBounds { .. }));
    }

    // r[verify fasta.fetch.bounds_check]
    #[test]
    fn fetch_reversed_range_errors() {
        let dir = TempDir::new().unwrap();
        let (fasta_path, _) = make_plain_fasta(&dir);

        let mut reader = IndexedFastaReader::open(&fasta_path).unwrap();
        // start > stop
        let err = reader.fetch_seq("seq1", p(6), p(3)).unwrap_err();
        assert!(matches!(err, FastaError::RegionOutOfBounds { .. }));
    }

    #[test]
    fn linewidth_equals_linebases_no_newlines() {
        // Degenerate case: entire sequence on one line, no trailing newline
        // before EOF (linewidth == linebases means zero newline bytes).
        let dir = TempDir::new().unwrap();
        let fasta_path = dir.path().join("noline.fa");
        let fai_path = dir.path().join("noline.fa.fai");

        let mut f = File::create(&fasta_path).unwrap();
        // Write without newlines between bases (only after header)
        f.write_all(b">seq1\nACGTACGT\n").unwrap();
        let mut idx = File::create(&fai_path).unwrap();
        // linewidth == linebases (8 == 8): the entire sequence fits on one
        // "line" with no newline mid-sequence. linewidth == 9 includes the
        // trailing \n, but since there's only one line, linebases == linewidth
        // effectively means "one line, no mid-sequence newlines". Actually for
        // a single-line sequence, samtools produces linebases==8, linewidth==9.
        // But linewidth==linebases is valid when a sequence exactly fills its
        // lines with no trailing newline.
        idx.write_all(b"seq1\t8\t6\t8\t8\n").unwrap();

        let mut reader = IndexedFastaReader::open(&fasta_path).unwrap();
        let seq = reader.fetch_seq("seq1", p(0), p(8)).unwrap();
        assert_eq!(seq, b"ACGTACGT");

        // Sub-range
        let seq = reader.fetch_seq("seq1", p(2), p(6)).unwrap();
        assert_eq!(seq, b"GTAC");
    }

    // r[verify fasta.plain.read_size]
    #[test]
    fn last_line_shorter_than_linebases() {
        // Sequence with a short last line — the FAI offset formula must still
        // work because byte_offset only does integer division on full lines.
        let dir = TempDir::new().unwrap();
        let fasta_path = dir.path().join("short_last.fa");
        let fai_path = dir.path().join("short_last.fa.fai");

        let mut f = File::create(&fasta_path).unwrap();
        // 7 bases, 4 per line: ACGT\nTGC\n (last line has 3 bases, not 4)
        f.write_all(b">seq1\nACGT\nTGC\n").unwrap();
        let mut idx = File::create(&fai_path).unwrap();
        idx.write_all(b"seq1\t7\t6\t4\t5\n").unwrap();

        let mut reader = IndexedFastaReader::open(&fasta_path).unwrap();

        // Full sequence
        assert_eq!(reader.fetch_seq("seq1", p(0), p(7)).unwrap(), b"ACGTTGC");

        // Just the last line
        assert_eq!(reader.fetch_seq("seq1", p(4), p(7)).unwrap(), b"TGC");

        // Last single base
        assert_eq!(reader.fetch_seq("seq1", p(6), p(7)).unwrap(), b"C");

        // Cross boundary into short last line
        assert_eq!(reader.fetch_seq("seq1", p(3), p(6)).unwrap(), b"TTG");
    }

    // r[verify fasta.fetch.raw_bytes]
    #[test]
    fn iupac_codes_preserved_as_raw_bytes() {
        let dir = TempDir::new().unwrap();
        let fasta_path = dir.path().join("iupac.fa");
        let fai_path = dir.path().join("iupac.fa.fai");

        let mut f = File::create(&fasta_path).unwrap();
        // Sequence with IUPAC ambiguity codes — these must NOT be converted to N
        f.write_all(b">seq1\nACMRWSYKVHDBN\n").unwrap();

        let mut idx = File::create(&fai_path).unwrap();
        idx.write_all(b"seq1\t13\t6\t13\t14\n").unwrap();

        let mut reader = IndexedFastaReader::open(&fasta_path).unwrap();
        let seq = reader.fetch_seq("seq1", p(0), p(13)).unwrap();
        assert_eq!(seq, b"ACMRWSYKVHDBN", "IUPAC codes must be preserved as raw bytes");

        // Lowercase IUPAC must be uppercased but not converted
        let dir2 = TempDir::new().unwrap();
        let fasta_path2 = dir2.path().join("iupac_lower.fa");
        let fai_path2 = dir2.path().join("iupac_lower.fa.fai");

        let mut f2 = File::create(&fasta_path2).unwrap();
        f2.write_all(b">seq1\nacmrwsykvhdbn\n").unwrap();
        let mut idx2 = File::create(&fai_path2).unwrap();
        idx2.write_all(b"seq1\t13\t6\t13\t14\n").unwrap();

        let mut reader2 = IndexedFastaReader::open(&fasta_path2).unwrap();
        let seq2 = reader2.fetch_seq("seq1", p(0), p(13)).unwrap();
        assert_eq!(
            seq2, b"ACMRWSYKVHDBN",
            "lowercase IUPAC must be uppercased, not converted to N"
        );
    }
}
