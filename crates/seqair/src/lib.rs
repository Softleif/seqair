//! Pure-Rust BAM/SAM/CRAM/FASTA reader, pileup engine, and VCF/BCF writer.
//!
//! seqair provides indexed random access to alignment files (BAM, SAM, CRAM),
//! reference sequences (FASTA), and a column-based pileup engine.
//! See [`Readers`] for examples.
//!
//! It also includes VCF/BCF writers with type-safe allele representation and
//! single-pass index co-production.
//! See [`vcf`] for examples.

pub mod bam;
pub mod cram;
pub mod fasta;
pub mod reader;
pub mod sam;
pub(crate) mod utils;
pub mod vcf;

pub use reader::Readers;
