//! CRAM v3.0/v3.1 reading. Use [`reader::IndexedCramReader`] to fetch records into a
//! [`crate::bam::RecordStore`]; the sub-modules handle the compression codec stack underneath.

// r[impl io.minimal_public_api]
pub mod arith;
pub mod bitstream;
pub mod block;
pub(crate) mod codec_io;
pub mod compression_header;
pub mod container;
pub mod encoding;
pub mod fqzcomp;
/// A plain, spec-literal fqzcomp decoder that tests and fuzz targets check
/// [`fqzcomp`] against.
#[cfg(any(test, feature = "fuzz"))]
#[doc(hidden)]
pub mod fqzcomp_reference;
pub mod index;
pub(crate) mod range_coder;
pub mod rans;
pub mod rans_nx16;
// The per-ISA Nx16 kernels became one `fearless_simd` kernel in `rans_nx16`.
// These were public (with only crate-private items) in 0.3.0; keep the paths
// so 0.3.1 stays semver-compatible.
#[doc(hidden)]
pub mod rans_nx16_avx2 {}
#[doc(hidden)]
pub mod rans_nx16_neon {}
pub mod reader;
pub mod slice;
pub mod tok3;
/// A plain, spec-literal tok3 decoder that tests and fuzz targets check
/// [`tok3`] against.
#[cfg(any(test, feature = "fuzz"))]
#[doc(hidden)]
pub mod tok3_reference;
pub mod varint;

pub use index::CramIndexError;
pub use reader::CramError;
