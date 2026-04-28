//! Format-agnostic entry point. Use [`Readers::open`] to open an alignment file (BAM/SAM/CRAM)
//! together with a FASTA reference, then call [`Readers::pileup`] to iterate columns or
//! [`Readers::fetch_base_seq`] for reference sequence. [`IndexedReader`] is the bare alignment
//! handle when no FASTA is needed. Both types are forkable for multi-threaded use.

use crate::{bam::BamError, cram::reader::CramError, fasta::FastaError, sam::reader::SamError};
use seqair_types::SmolStr;
use std::path::PathBuf;

mod formats;
mod indexed;
mod readers;
mod resolve;
mod segment;

#[cfg(feature = "fuzz")]
mod fuzz;

pub use formats::FormatDetectionError;
pub use indexed::{FetchCounts, IndexedReader};
pub use readers::Readers;
pub use resolve::{ResolveTid, Tid, TidError};
pub use segment::{IntoSegmentTarget, Segment, SegmentOptions, SegmentOptionsError, Segments};

#[cfg(feature = "fuzz")]
pub use fuzz::FuzzReaders;

// r[impl io.non_exhaustive_enums]
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

    #[error("failed to fetch reference sequence for {contig}:{start}-{end}")]
    FastaFetch { contig: SmolStr, start: u64, end: u64, source: FastaError },

    #[error("could not resolve target id")]
    Tid {
        #[from]
        source: resolve::TidError,
    },

    #[error("contig '{name}' has zero length; cannot resolve a region")]
    EmptyContig { name: SmolStr },

    #[error("region end {end} exceeds i32::MAX; not representable as Pos0")]
    RegionEndTooLarge { end: u64 },

    #[error("region start {start} > end {end} on contig '{contig}'")]
    RegionStartAfterEnd { contig: SmolStr, start: u64, end: u64 },

    #[error(
        "segment for contig '{contig}' (tid {expected_tid}) does not match this Readers' header"
    )]
    SegmentHeaderMismatch { contig: SmolStr, expected_tid: u32 },

    #[error("BAM header target count {count} does not fit in u32")]
    HeaderTargetCountOverflow { count: usize },
}
