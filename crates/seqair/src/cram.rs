//! CRAM v3.0/v3.1 reading. Use [`reader::IndexedCramReader`] to fetch records into a
//! [`crate::bam::RecordStore`]; the sub-modules handle the compression codec stack underneath.

// r[impl unified.minimal_public_api]
pub mod bitstream;
pub mod block;
pub mod compression_header;
pub mod container;
pub mod encoding;
pub mod index;
pub mod rans;
pub mod rans_nx16;
pub mod reader;
pub mod slice;
pub mod tok3;
pub mod varint;

pub use index::CramIndexError;
pub use reader::CramError;
