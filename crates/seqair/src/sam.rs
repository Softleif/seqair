//! bgzf-compressed SAM reading. Use [`reader::IndexedSamReader`] to fetch records into a
//! [`crate::bam::RecordStore`] via a tabix index.

// r[impl io.minimal_public_api]
pub mod reader;

pub use reader::{SamError, SamRecordError};
