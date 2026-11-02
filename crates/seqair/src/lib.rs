//! Pure-Rust BAM/SAM/CRAM/FASTA reader, pileup engine, and VCF/BCF writer.
//!
//! seqair provides indexed random access to alignment files (BAM, SAM, CRAM),
//! reference sequences (FASTA), and a column-based pileup engine.
//! See [`Readers`] for examples.
//!
//! It also includes VCF/BCF writers with type-safe allele representation and
//! single-pass index co-production.
//! See [`vcf`] for examples.

// BAM, CRAM and BCF are little-endian on disk, and the readers and writers
// view CIGARs, qualities and integer arrays in place as their native types.
#[cfg(not(target_endian = "little"))]
compile_error!("seqair supports little-endian targets only");

pub mod bam;
pub mod cram;
pub mod fasta;
pub mod io;
#[cfg(test)]
mod pinned;
pub mod reader;
pub mod sam;
#[cfg(test)]
mod simd_levels;
pub(crate) mod utils;
pub mod vcf;

pub use reader::Readers;
