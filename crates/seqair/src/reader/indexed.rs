use super::ReaderError;
use super::formats::{self, Format};
use crate::{
    bam::{
        BamHeader, IndexedBamReader,
        record_store::{CustomizeRecordStore, RecordStore},
    },
    cram::reader::IndexedCramReader,
    sam::reader::IndexedSamReader,
};
use seqair_types::Pos0;
use std::{
    io::{Read, Seek},
    path::Path,
};
use tracing::instrument;

// r[impl unified.fetch_counts]
/// Return value of the customize-aware `fetch_into_customized` methods on each
/// reader. `kept <= fetched` always holds.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct FetchCounts {
    /// Records that passed the reader's built-in overlap/unmapped checks,
    /// before the user's pre-filter. This is the count you'd get from a
    /// filter-free fetch.
    pub fetched: usize,
    /// Records that also passed the user's pre-filter and remain in the store.
    /// Equal to `fetched` when no filter is installed.
    pub kept: usize,
}

/// Format-agnostic indexed alignment reader, generic over the I/O backend.
///
/// `R = File` for production (file-based I/O); `R = Cursor<Vec<u8>>` for fuzzing
/// (in-memory, no file I/O).
// r[impl unified.reader_enum]
#[non_exhaustive]
pub enum IndexedReader<R: Read + Seek = std::fs::File> {
    Bam(IndexedBamReader<R>),
    Sam(IndexedSamReader<R>),
    Cram(Box<IndexedCramReader<R>>),
}

/// Cursor-backed type alias for in-memory fuzzing.
#[cfg(feature = "fuzz")]
pub type CursorReader = IndexedReader<std::io::Cursor<Vec<u8>>>;

impl<R: Read + Seek> std::fmt::Debug for IndexedReader<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Bam(r) => r.fmt(f),
            Self::Sam(r) => r.fmt(f),
            Self::Cram(r) => r.fmt(f),
        }
    }
}

impl<R: Read + Seek> IndexedReader<R> {
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
        start: Pos0,
        end: Pos0,
        store: &mut RecordStore,
    ) -> Result<usize, ReaderError> {
        self.fetch_into_customized(tid, start, end, store, &mut ()).map(|c| c.kept)
    }

    // r[impl unified.fetch_into_customized]
    /// Customized variant of [`Self::fetch_into`]. For each record that would
    /// normally enter the store, `customize.filter` is invoked; if it
    /// returns `false`, the push is rolled back (zero slab waste, see
    /// [`RecordStore::push_raw`]). The returned [`FetchCounts`]
    /// distinguishes records the reader produced (`fetched`) from those
    /// the filter retained (`kept`). Pass `&mut ()` for no filtering.
    pub fn fetch_into_customized<E: CustomizeRecordStore>(
        &mut self,
        tid: u32,
        start: Pos0,
        end: Pos0,
        store: &mut RecordStore<E::Extra>,
        customize: &mut E,
    ) -> Result<FetchCounts, ReaderError> {
        match self {
            Self::Bam(r) => r
                .fetch_into_customized(tid, start, end, store, customize)
                .map_err(ReaderError::from),
            Self::Sam(r) => r
                .fetch_into_customized(tid, start, end, store, customize)
                .map_err(ReaderError::from),
            Self::Cram(r) => r
                .fetch_into_customized(tid, start, end, store, customize)
                .map_err(ReaderError::from),
        }
    }
}

// r[impl unified.reader_api]
impl IndexedReader<std::fs::File> {
    /// Open a BAM, bgzf-compressed SAM, or CRAM file, auto-detecting the format.
    ///
    /// For CRAM, opening without a FASTA always succeeds; missing-reference
    /// errors only surface at fetch time when a slice actually needs an
    /// external reference (`r[cram.fasta.optional]`).
    // r[impl unified.readers_backward_compat]
    #[instrument(level = "debug", fields(path = %path.display()))]
    pub fn open(path: &Path) -> Result<Self, ReaderError> {
        Self::open_with_optional_fasta(path, None)
    }

    /// Open any format, with a FASTA path for CRAM support.
    pub(crate) fn open_with_fasta(path: &Path, fasta_path: &Path) -> Result<Self, ReaderError> {
        Self::open_with_optional_fasta(path, Some(fasta_path))
    }

    fn open_with_optional_fasta(
        path: &Path,
        fasta_path: Option<&Path>,
    ) -> Result<Self, ReaderError> {
        match formats::detect(path)? {
            Format::Bam => {
                let reader = IndexedBamReader::open(path)?;
                Ok(IndexedReader::Bam(reader))
            }
            Format::Sam => {
                let reader = IndexedSamReader::open(path)?;
                Ok(IndexedReader::Sam(reader))
            }
            Format::Cram => {
                let reader = match fasta_path {
                    Some(fp) => IndexedCramReader::open(path, fp)?,
                    None => IndexedCramReader::open_without_reference(path)?,
                };
                Ok(IndexedReader::Cram(Box::new(reader)))
            }
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write as _;
    use tempfile::NamedTempFile;

    // r[verify cram.fasta.optional]
    /// `IndexedReader::open` MUST accept CRAM files without a FASTA. The
    /// reader fails with a downstream CRAM error (here, missing CRAI for our
    /// stub file) rather than with `CramRequiresFasta`.
    #[test]
    fn cram_open_without_fasta_does_not_short_circuit() {
        let mut f = NamedTempFile::new().expect("tempfile");
        f.write_all(b"CRAM\x03\x00").expect("write");
        f.flush().expect("flush");
        let err = IndexedReader::open(f.path()).unwrap_err();
        assert!(
            !matches!(&err, ReaderError::Format { .. }),
            "open should reach the CRAM reader, not error during format detection: {err}"
        );
    }
}
