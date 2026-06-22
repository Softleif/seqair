use super::ReaderError;
use super::formats::{self, Format, FormatDetectionError};
use super::readers::Readers;
use crate::{
    bam::{
        BamHeader, IndexedBamReader,
        record_store::{CustomizeRecordStore, RecordStore},
    },
    cram::reader::IndexedCramReader,
    fasta::IndexedFastaReader,
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

    /// Estimate the compressed bytes a `[start, end]` region query would load,
    /// for byte-aware segmentation. Returns `None` for CRAM, whose slice-based
    /// reader bounds memory differently (there is no `RegionBuf` bulk-load), so
    /// callers should skip byte-budget subdivision for it.
    pub fn estimate_region_bytes(&self, tid: u32, start: Pos0, end: Pos0) -> Option<u64> {
        match self {
            Self::Bam(r) => Some(r.estimate_region_bytes(tid, start, end)),
            Self::Sam(r) => Some(r.estimate_region_bytes(tid, start, end)),
            Self::Cram(_) => None,
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
    /// Open a BAM or bgzf-compressed SAM file, auto-detecting the format.
    ///
    /// CRAM files are detected but require a reference. To stream CRAM
    /// records without the pileup engine, use
    /// [`open_with_reference`](Self::open_with_reference); for reference-aware
    /// pileup, use [`crate::Readers::open`]. This method returns an error for
    /// CRAM with a message directing the user to provide a reference.
    // r[impl unified.readers_backward_compat]
    #[instrument(level = "debug", fields(path = %path.display()))]
    pub fn open(path: &Path) -> Result<Self, ReaderError> {
        match formats::detect(path)? {
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

    /// Open any alignment format, supplying a FASTA reference for CRAM.
    ///
    /// This is the **record-iteration** handle for all three formats: it does
    /// not build a pileup engine and allocates no reference-fetch buffer. For
    /// BAM/SAM the reference is not consulted (records carry their own bases);
    /// for CRAM it is required to reconstruct read bases during decoding. Use
    /// this when you want to stream records from any format without paying for
    /// pileup machinery. For reference-aware pileup, see [`crate::Readers`], or
    /// call [`into_pileup`](Self::into_pileup) to upgrade afterwards.
    ///
    /// ```no_run
    /// use seqair::reader::IndexedReader;
    /// use seqair::bam::record_store::RecordStore;
    /// use seqair_types::Pos0;
    /// use std::path::Path;
    ///
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// // Records only, any format (here CRAM), no pileup cost.
    /// let mut reader =
    ///     IndexedReader::open_with_reference(Path::new("sample.cram"), Path::new("ref.fa"))?;
    /// let tid = reader.header().tid("chr1").expect("chr1 in header");
    /// let mut store = RecordStore::default();
    /// let (start, end) = (Pos0::new(1_000).unwrap(), Pos0::new(2_000).unwrap());
    /// reader.fetch_into(tid, start, end, &mut store)?;
    /// # Ok(())
    /// # }
    /// ```
    // r[impl unified.open_with_reference]
    pub fn open_with_reference(path: &Path, reference_path: &Path) -> Result<Self, ReaderError> {
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
                let reader = IndexedCramReader::open(path, reference_path)?;
                Ok(IndexedReader::Cram(Box::new(reader)))
            }
        }
    }

    /// Upgrade this record reader into a reference-aware [`Readers`] for pileup.
    ///
    /// The alignment handle is reused as-is; `fasta` supplies the per-column
    /// reference sequence consumed by the pileup engine. Note that the CRAM
    /// decode reference (passed at open time) and this analysis reference are
    /// independent handles — this method does not check they point at the same
    /// file, matching the existing two-handle design.
    ///
    /// ```no_run
    /// use seqair::reader::IndexedReader;
    /// use seqair::fasta::IndexedFastaReader;
    /// use std::path::Path;
    ///
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// // Opened a records-only reader earlier; now decide we want pileup.
    /// let reader =
    ///     IndexedReader::open_with_reference(Path::new("sample.cram"), Path::new("ref.fa"))?;
    /// let readers = reader.into_pileup(IndexedFastaReader::open(Path::new("ref.fa"))?);
    /// // `readers` is a full reference-aware `Readers` ready for `segments()`/`pileup()`.
    /// # let _ = readers;
    /// # Ok(())
    /// # }
    /// ```
    // r[impl unified.into_pileup]
    pub fn into_pileup(self, fasta: IndexedFastaReader) -> Readers<()> {
        Readers::from_parts(self, fasta)
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
