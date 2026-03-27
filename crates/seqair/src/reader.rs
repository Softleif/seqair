//! Format-agnostic entry point. Use [`Readers::open`] to open an alignment file (BAM/SAM/CRAM)
//! together with a FASTA reference, then call [`Readers::pileup`] to iterate columns or
//! [`Readers::fetch_base_seq`] for reference sequence. [`IndexedReader`] is the bare alignment
//! handle when no FASTA is needed. Both types are forkable for multi-threaded use.

use crate::{
    bam::{
        BamError, BamHeader, IndexedBamReader, bgzf::BgzfReader, pileup::PileupEngine,
        record_store::RecordStore,
    },
    cram::reader::{CramError, IndexedCramReader},
    fasta::{FastaError, IndexedFastaReader},
    sam::reader::{IndexedSamReader, SamError},
};
use seqair_types::Base;
use std::path::{Path, PathBuf};
use std::rc::Rc;
use tracing::instrument;

// r[impl io.non_exhaustive_enums]
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum FormatDetectionError {
    #[error(
        "CRAM format detected but no FASTA reference provided. \
         Use Readers::open(alignment_path, fasta_path) instead, \
         or convert to BAM with `samtools view -b`"
    )]
    CramRequiresFasta,

    #[error("file too short to determine format ({len} bytes)")]
    FileTooShort { path: PathBuf, len: u64 },

    #[error(
        "uncompressed SAM detected. \
         Compress with `bgzip {path}` then index with `tabix -p sam {path}.gz`"
    )]
    UncompressedSam { path: PathBuf },

    #[error(
        "could not decompress first BGZF block from {path}. \
         If this is a plain gzip file, use `bgzip` instead of `gzip`"
    )]
    NotBgzf { path: PathBuf, source: crate::bam::bgzf::BgzfError },

    #[error(
        "BGZF-compressed file {path} does not contain BAM or SAM data \
         (first decompressed bytes: {magic:?}). \
         Supported formats: BAM (.bam), bgzf-compressed SAM (.sam.gz), CRAM (.cram)."
    )]
    UnrecognizedBgzfContent { path: PathBuf, magic: [u8; 4] },

    #[error(
        "unrecognized file format for {path} (magic bytes: {magic:02x?}). \
         Supported formats: BAM (.bam), bgzf-compressed SAM (.sam.gz), CRAM (.cram)."
    )]
    UnrecognizedFormat { path: PathBuf, magic: [u8; 4] },
}

#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum ReaderError {
    #[error("BAM reader error")]
    Bam {
        #[from]
        source: BamError,
    },

    #[error("SAM reader error")]
    Sam {
        #[from]
        source: SamError,
    },

    #[error("CRAM reader error")]
    Cram {
        #[from]
        source: CramError,
    },

    #[error("I/O error reading {path}")]
    Open { path: PathBuf, source: std::io::Error },

    #[error(transparent)]
    Format {
        #[from]
        source: FormatDetectionError,
    },

    #[error("failed to open FASTA reference")]
    FastaOpen { source: FastaError },

    #[error("failed to fork FASTA reader")]
    FastaFork { source: FastaError },

    #[error("region coordinate {value} overflows i64")]
    RegionOverflow { value: u64 },
}

/// Format-agnostic indexed alignment reader.
///
/// Auto-detects BAM, bgzf-compressed SAM, or CRAM by inspecting magic bytes.
/// All formats populate the same `RecordStore` via `fetch_into()`.
// r[impl unified.reader_enum]
#[non_exhaustive]
pub enum IndexedReader {
    Bam(IndexedBamReader),
    Sam(IndexedSamReader),
    Cram(Box<IndexedCramReader>),
}

impl std::fmt::Debug for IndexedReader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Bam(r) => r.fmt(f),
            Self::Sam(r) => r.fmt(f),
            Self::Cram(r) => r.fmt(f),
        }
    }
}

// r[impl unified.reader_api]
impl IndexedReader {
    /// Open a BAM or bgzf-compressed SAM file, auto-detecting the format.
    ///
    /// CRAM files are detected but require a FASTA reference — use
    /// [`Readers::open`] instead. This method returns an error for CRAM
    /// with a message directing the user to provide a reference.
    // r[impl unified.readers_backward_compat]
    #[instrument(level = "debug", fields(path = %path.display()), err)]
    pub fn open(path: &Path) -> Result<Self, ReaderError> {
        match detect_format(path)? {
            Format::Bam => {
                let reader = IndexedBamReader::open(path)?;
                Ok(IndexedReader::Bam(reader))
            }
            Format::Sam => {
                let reader = IndexedSamReader::open(path)?;
                Ok(IndexedReader::Sam(reader))
            }
            Format::Cram => Err(FormatDetectionError::CramRequiresFasta.into()),
        }
    }

    /// Open any format, with a FASTA path for CRAM support.
    fn open_with_fasta(path: &Path, fasta_path: &Path) -> Result<Self, ReaderError> {
        match detect_format(path)? {
            Format::Bam => {
                let reader = IndexedBamReader::open(path)?;
                Ok(IndexedReader::Bam(reader))
            }
            Format::Sam => {
                let reader = IndexedSamReader::open(path)?;
                Ok(IndexedReader::Sam(reader))
            }
            Format::Cram => {
                let reader = IndexedCramReader::open(path, fasta_path)?;
                Ok(IndexedReader::Cram(Box::new(reader)))
            }
        }
    }

    pub fn header(&self) -> &BamHeader {
        match self {
            Self::Bam(r) => r.header(),
            Self::Sam(r) => r.header(),
            Self::Cram(r) => r.header(),
        }
    }

    // r[impl unified.fetch_equivalence]
    pub fn fetch_into(
        &mut self,
        tid: u32,
        start: u64,
        end: u64,
        store: &mut RecordStore,
    ) -> Result<usize, ReaderError> {
        match self {
            Self::Bam(r) => r.fetch_into(tid, start, end, store).map_err(ReaderError::from),
            Self::Sam(r) => r.fetch_into(tid, start, end, store).map_err(ReaderError::from),
            Self::Cram(r) => r.fetch_into(tid, start, end, store).map_err(ReaderError::from),
        }
    }

    // r[impl unified.fork_bam]
    // r[impl unified.fork_sam]
    // r[impl unified.fork_cram]
    pub fn fork(&self) -> Result<Self, ReaderError> {
        match self {
            Self::Bam(r) => Ok(Self::Bam(r.fork()?)),
            Self::Sam(r) => Ok(Self::Sam(r.fork()?)),
            Self::Cram(r) => Ok(Self::Cram(Box::new(r.fork()?))),
        }
    }
}

/// Alignment + reference reader bundle.
///
/// Bundles an [`IndexedReader`] (BAM/SAM/CRAM) with an [`IndexedFastaReader`]
/// so that CRAM has access to the reference it needs and all formats have
/// uniform open/fork/fetch semantics.
// r[impl unified.readers_struct]
pub struct Readers {
    alignment: IndexedReader,
    fasta: IndexedFastaReader,
    store: RecordStore,
    fasta_buf: Vec<u8>,
}

impl std::fmt::Debug for Readers {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Readers").field("alignment", &self.alignment).finish()
    }
}

impl Readers {
    /// Open an alignment file (BAM/SAM/CRAM) and a FASTA reference.
    ///
    /// Auto-detects the alignment format. For CRAM, the FASTA path is passed
    /// to the CRAM reader for sequence reconstruction.
    // r[impl unified.readers_open]
    #[instrument(level = "debug", fields(alignment = %alignment_path.display(), fasta = %fasta_path.display()), err)]
    pub fn open(alignment_path: &Path, fasta_path: &Path) -> Result<Self, ReaderError> {
        let fasta = IndexedFastaReader::open(fasta_path)
            .map_err(|source| ReaderError::FastaOpen { source })?;
        let alignment = IndexedReader::open_with_fasta(alignment_path, fasta_path)?;
        Ok(Readers { alignment, fasta, store: RecordStore::default(), fasta_buf: Vec::new() })
    }

    /// Fork both the alignment reader and the FASTA reader.
    // r[impl unified.readers_fork]
    pub fn fork(&self) -> Result<Self, ReaderError> {
        let alignment = self.alignment.fork()?;
        let fasta = self.fasta.fork().map_err(|source| ReaderError::FastaFork { source })?;
        Ok(Readers { alignment, fasta, store: RecordStore::default(), fasta_buf: Vec::new() })
    }

    // r[impl unified.readers_accessors]
    pub fn header(&self) -> &BamHeader {
        self.alignment.header()
    }

    pub fn fetch_into(
        &mut self,
        tid: u32,
        start: u64,
        end: u64,
        store: &mut RecordStore,
    ) -> Result<usize, ReaderError> {
        self.alignment.fetch_into(tid, start, end, store)
    }

    /// Fetch records for a region and return a [`PileupEngine`] ready for iteration.
    ///
    /// Uses an internal [`RecordStore`] whose capacity is retained across calls.
    /// After iterating the engine, call [`recover_store`](Self::recover_store) to
    /// return the store for reuse. If not called, the next `pileup()` call
    /// allocates a fresh store (small perf hit, not a correctness issue).
    pub fn pileup(&mut self, tid: u32, start: u64, end: u64) -> Result<PileupEngine, ReaderError> {
        let start_i64 =
            i64::try_from(start).map_err(|_| ReaderError::RegionOverflow { value: start })?;
        let end_i64 = i64::try_from(end).map_err(|_| ReaderError::RegionOverflow { value: end })?;

        self.alignment.fetch_into(tid, start, end, &mut self.store)?;

        let store = std::mem::take(&mut self.store);
        Ok(PileupEngine::new(store, start_i64, end_i64))
    }

    /// Recover the [`RecordStore`] from a consumed [`PileupEngine`] for reuse.
    ///
    /// Call this after iteration is complete. The store retains its allocated
    /// capacity, avoiding ~39 MB of re-allocation on the next `pileup()` call.
    pub fn recover_store(&mut self, engine: &mut PileupEngine) {
        if let Some(store) = engine.take_store() {
            self.store = store;
        }
    }

    pub fn fasta(&self) -> &IndexedFastaReader {
        &self.fasta
    }

    pub fn fasta_mut(&mut self) -> &mut IndexedFastaReader {
        &mut self.fasta
    }

    // r[impl fasta.fetch.buffer_reuse]
    /// Fetch a reference sequence region and return it as `Rc<[Base]>`.
    ///
    /// Uses an internal buffer whose capacity is retained across calls,
    /// avoiding per-segment allocation. The returned `Rc<[Base]>` owns its
    /// own copy — the internal buffer is reused on the next call.
    pub fn fetch_base_seq(
        &mut self,
        name: &str,
        start: u64,
        stop: u64,
    ) -> Result<Rc<[Base]>, FastaError> {
        self.fasta.fetch_seq_into(name, start, stop, &mut self.fasta_buf)?;
        // Take the buffer so from_ascii_vec can reinterpret it in-place (safe),
        // then create the Rc (which copies). The buffer capacity is lost but
        // fasta_buf re-grows on the next call.
        let buf = std::mem::take(&mut self.fasta_buf);
        let bases = Base::from_ascii_vec(buf);
        Ok(Rc::from(bases))
    }

    pub fn alignment(&self) -> &IndexedReader {
        &self.alignment
    }

    pub fn alignment_mut(&mut self) -> &mut IndexedReader {
        &mut self.alignment
    }
}

#[derive(Debug, Clone, Copy)]
enum Format {
    Bam,
    Sam,
    Cram,
}

// r[impl unified.detect_format]
fn detect_format(path: &Path) -> Result<Format, ReaderError> {
    let mut magic = [0u8; 4];
    let file_len = std::fs::metadata(path)
        .map_err(|source| ReaderError::Open { path: path.to_path_buf(), source })?
        .len();

    if file_len < 4 {
        return Err(
            FormatDetectionError::FileTooShort { path: path.to_path_buf(), len: file_len }.into()
        );
    }

    let mut file = std::fs::File::open(path)
        .map_err(|source| ReaderError::Open { path: path.to_path_buf(), source })?;

    use std::io::Read;
    file.read_exact(&mut magic)
        .map_err(|source| ReaderError::Open { path: path.to_path_buf(), source })?;

    // CRAM: magic bytes "CRAM"
    if magic == *b"CRAM" {
        return Ok(Format::Cram);
    }

    // Uncompressed SAM: starts with '@'
    if magic[0] == b'@' {
        return Err(FormatDetectionError::UncompressedSam { path: path.to_path_buf() }.into());
    }

    // BGZF: gzip magic 1f 8b
    if magic[0] == 0x1f && magic[1] == 0x8b {
        // Decompress first block to check BAM magic vs SAM text
        let mut bgzf = BgzfReader::open(path)
            .map_err(|source| FormatDetectionError::NotBgzf { path: path.to_path_buf(), source })?;

        let mut first_bytes = [0u8; 4];
        bgzf.read_exact_into(&mut first_bytes)
            .map_err(|source| FormatDetectionError::NotBgzf { path: path.to_path_buf(), source })?;

        if first_bytes == *b"BAM\x01" {
            return Ok(Format::Bam);
        }

        if first_bytes[0] == b'@' {
            return Ok(Format::Sam);
        }

        return Err(FormatDetectionError::UnrecognizedBgzfContent {
            path: path.to_path_buf(),
            magic: first_bytes,
        }
        .into());
    }

    Err(FormatDetectionError::UnrecognizedFormat { path: path.to_path_buf(), magic }.into())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_bgzf(content: &[u8]) -> Vec<u8> {
        use bgzf::{CompressionLevel, Writer as BgzfWriter};
        let mut output = Vec::new();
        let mut writer =
            BgzfWriter::new(&mut output, CompressionLevel::new(1).expect("valid level"));
        writer.write_all(content).expect("write");
        writer.finish().expect("finish");
        output
    }

    // r[verify unified.detect_format]
    #[test]
    fn file_too_short() {
        let mut f = NamedTempFile::new().expect("tempfile");
        f.write_all(b"BA").expect("write");
        f.flush().expect("flush");
        let err = detect_format(f.path()).unwrap_err();
        assert!(
            matches!(
                err,
                ReaderError::Format { source: FormatDetectionError::FileTooShort { len: 2, .. } }
            ),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn file_empty() {
        let f = NamedTempFile::new().expect("tempfile");
        let err = detect_format(f.path()).unwrap_err();
        assert!(
            matches!(
                err,
                ReaderError::Format { source: FormatDetectionError::FileTooShort { len: 0, .. } }
            ),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn uncompressed_sam() {
        let mut f = NamedTempFile::new().expect("tempfile");
        f.write_all(b"@HD\tVN:1.6\n").expect("write");
        f.flush().expect("flush");
        let err = detect_format(f.path()).unwrap_err();
        assert!(
            matches!(
                err,
                ReaderError::Format { source: FormatDetectionError::UncompressedSam { .. } }
            ),
            "unexpected error: {err}"
        );
    }

    // TODO: `FormatDetectionError::NotBgzf` is only triggered when `BgzfReader::open`
    // fails with an I/O error. Because `detect_format` already opens the file once before
    // calling `BgzfReader::open`, this path is not reachable via content manipulation —
    // it would require a race condition (file deleted or permissions revoked between the
    // two opens). Portable unit testing of this variant is not feasible.
    //
    #[test]
    fn not_valid_bgzf_returns_not_bgzf() {
        // gzip magic but invalid BGZF header (byte 3 is 0x00, not 0x04 FEXTRA flag).
        // BgzfReader::open succeeds (pure file open); read_exact_into fails → NotBgzf.
        let mut f = NamedTempFile::new().expect("tempfile");
        f.write_all(&[0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00]).expect("write");
        f.flush().expect("flush");
        let err = detect_format(f.path()).unwrap_err();
        assert!(
            matches!(err, ReaderError::Format { source: FormatDetectionError::NotBgzf { .. } }),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn short_bgzf_content_returns_not_bgzf() {
        // A valid BGZF block with only 2 bytes of decompressed content — too short
        // to read 4 bytes, triggering NotBgzf.
        let compressed = write_bgzf(b"AB");
        let mut f = NamedTempFile::new().expect("tempfile");
        f.write_all(&compressed).expect("write");
        f.flush().expect("flush");
        let err = detect_format(f.path()).unwrap_err();
        assert!(
            matches!(err, ReaderError::Format { source: FormatDetectionError::NotBgzf { .. } }),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn unrecognized_bgzf_content() {
        // A valid BGZF block whose decompressed content starts with neither
        // "BAM\x01" nor "@" — triggers UnrecognizedBgzfContent.
        let compressed = write_bgzf(b"\xDE\xAD\xBE\xEF and some more data");
        let mut f = NamedTempFile::new().expect("tempfile");
        f.write_all(&compressed).expect("write");
        f.flush().expect("flush");
        let err = detect_format(f.path()).unwrap_err();
        assert!(
            matches!(
                err,
                ReaderError::Format {
                    source: FormatDetectionError::UnrecognizedBgzfContent {
                        magic: [0xDE, 0xAD, 0xBE, 0xEF],
                        ..
                    }
                }
            ),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn bgzf_compressed_sam_detected() {
        // Sanity-check the Sam fast-path: "@"-prefixed BGZF content → Format::Sam
        let compressed = write_bgzf(b"@HD\tVN:1.6\n");
        let mut f = NamedTempFile::new().expect("tempfile");
        f.write_all(&compressed).expect("write");
        f.flush().expect("flush");
        // detect_format should succeed with Format::Sam (no error)
        let result = detect_format(f.path());
        assert!(matches!(result, Ok(Format::Sam)), "expected Sam, got {result:?}");
    }

    #[test]
    fn unrecognized_format() {
        // 4+ raw bytes that are not gzip magic, not "@", not "CRAM"
        let mut f = NamedTempFile::new().expect("tempfile");
        f.write_all(&[0xDE, 0xAD, 0xBE, 0xEF]).expect("write");
        f.flush().expect("flush");
        let err = detect_format(f.path()).unwrap_err();
        assert!(
            matches!(
                err,
                ReaderError::Format {
                    source: FormatDetectionError::UnrecognizedFormat {
                        magic: [0xDE, 0xAD, 0xBE, 0xEF],
                        ..
                    }
                }
            ),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn cram_requires_fasta() {
        // A file starting with the CRAM magic bytes — IndexedReader::open (no FASTA)
        // must return CramRequiresFasta.
        let mut f = NamedTempFile::new().expect("tempfile");
        f.write_all(b"CRAM\x03\x00").expect("write");
        f.flush().expect("flush");
        let err = IndexedReader::open(f.path()).unwrap_err();
        assert!(
            matches!(err, ReaderError::Format { source: FormatDetectionError::CramRequiresFasta }),
            "unexpected error: {err}"
        );
    }
}
