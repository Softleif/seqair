use super::ReaderError;
use super::formats::{self, Format, FormatDetectionError};
use crate::{
    bam::{BamHeader, IndexedBamReader, record_store::RecordStore},
    cram::reader::IndexedCramReader,
    sam::reader::IndexedSamReader,
};
use seqair_types::{Pos, Zero};
use std::{
    io::{Read, Seek},
    path::Path,
};
use tracing::instrument;

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
        start: Pos<Zero>,
        end: Pos<Zero>,
        store: &mut RecordStore,
    ) -> Result<usize, ReaderError> {
        match self {
            Self::Bam(r) => r.fetch_into(tid, start, end, store).map_err(ReaderError::from),
            Self::Sam(r) => r.fetch_into(tid, start, end, store).map_err(ReaderError::from),
            Self::Cram(r) => r.fetch_into(tid, start, end, store).map_err(ReaderError::from),
        }
    }
}

// r[impl unified.reader_api]
impl IndexedReader<std::fs::File> {
    /// Open a BAM or bgzf-compressed SAM file, auto-detecting the format.
    ///
    /// CRAM files are detected but require a FASTA reference — use
    /// [`crate::Readers::open`] instead. This method returns an error for CRAM
    /// with a message directing the user to provide a reference.
    // r[impl unified.readers_backward_compat]
    #[instrument(level = "debug", fields(path = %path.display()), err)]
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

    /// Open any format, with a FASTA path for CRAM support.
    pub(crate) fn open_with_fasta(path: &Path, fasta_path: &Path) -> Result<Self, ReaderError> {
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
                let reader = IndexedCramReader::open(path, fasta_path)?;
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
