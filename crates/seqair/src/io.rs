//! Shared I/O primitives used by both BAM and VCF writers.
//!
//! This module provides BGZF compression and index-building types that are
//! format-agnostic. Use these when working with [`crate::bam::BamWriter`] or
//! [`crate::vcf::Writer`] output.

mod index_builder;
pub(crate) mod text_scan;

// BGZF compression layer
pub use crate::bam::bgzf::{BgzfError, VirtualOffset};
pub use crate::bam::bgzf_writer::BgzfWriter;

// Single-pass index co-production
pub use index_builder::{IndexBuilder, IndexChunk, IndexError, reg2bin};
