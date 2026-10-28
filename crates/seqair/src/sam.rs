//! SAM reading. [`reader::IndexedSamReader`] fetches records of a
//! bgzf-compressed, indexed SAM into a [`crate::bam::RecordStore`];
//! [`SamRecordEncoder`] turns one SAM line
//! into a BAM record as htslib does, and [`format::format_sam_record`] turns
//! one back into SAM text as htslib does.

// r[impl io.minimal_public_api]
pub mod encode;
pub mod format;
pub mod reader;

pub use encode::{SamRecordEncoder, SamRecordError, parse_cigar};
pub use reader::SamError;
