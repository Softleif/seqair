//! VCF/BCF writing with type-safe alleles and zero-allocation encoding.
//!
//! # Two encoding paths
//!
//! 1. **[`VcfRecord`] path** — construct records via [`VcfRecordBuilder`], pass to
//!    [`VcfWriter`](writer::VcfWriter) or [`BcfWriter`](bcf_writer::BcfWriter).
//!    Simple and correct; suitable for most use cases.
//!
//! 2. **[`encoder`] path** — pre-resolve [`ScalarInfoHandle`](encoder::ScalarInfoHandle),
//!    [`GtFormatHandle`](encoder::GtFormatHandle), etc. at setup time, then encode
//!    directly into BCF buffers with zero allocations per record. Use this when
//!    writing millions of records in a hot loop (e.g., per-base variant calling).
//!
//! Both paths produce identical BCF output (verified by proptest equivalence tests).
//!
//! # Output formats
//!
//! Use [`open_writer`] with [`OutputFormat`] to create a writer for any format:
//! - `.vcf` → plain text VCF
//! - `.vcf.gz` → BGZF-compressed VCF with automatic TBI index co-production
//! - `.bcf` → BCF binary with automatic CSI index co-production

pub mod alleles;
pub(crate) mod bcf_encoding;
pub mod bcf_writer;
pub mod encoder;
pub mod error;
pub mod header;
pub mod index_builder;
pub mod record;
pub mod writer;

pub use alleles::Alleles;
pub use error::{AllelesError, VcfEncodeError, VcfError, VcfHeaderError};
pub use header::{
    ContigDef, FilterDef, FormatDef, InfoDef, Number, ValueType, VcfHeader, VcfHeaderBuilder,
};
pub use record::{
    Filters, Genotype, InfoFields, InfoValue, SampleFields, SampleValue, VcfRecord,
    VcfRecordBuilder,
};

use std::io::Write;
use std::sync::Arc;

/// Trait for VCF/BCF writers.
pub trait VcfWrite {
    /// Write the file header. Must be called exactly once before any record.
    fn write_header(&mut self) -> Result<(), VcfError>;
    /// Write a single record.
    fn write_record(&mut self, record: &VcfRecord) -> Result<(), VcfError>;
    /// Finalize the writer. Returns the index builder if one was co-produced.
    fn finish(self: Box<Self>) -> Result<Option<index_builder::IndexBuilder>, VcfError>;
}

impl<W: Write> VcfWrite for writer::VcfWriter<W> {
    fn write_header(&mut self) -> Result<(), VcfError> {
        self.write_header()
    }
    fn write_record(&mut self, record: &VcfRecord) -> Result<(), VcfError> {
        self.write_record(record)
    }
    fn finish(self: Box<Self>) -> Result<Option<index_builder::IndexBuilder>, VcfError> {
        (*self).finish()
    }
}

impl<W: Write> VcfWrite for bcf_writer::BcfWriter<W> {
    fn write_header(&mut self) -> Result<(), VcfError> {
        self.write_header()
    }
    fn write_record(&mut self, record: &VcfRecord) -> Result<(), VcfError> {
        self.write_vcf_record(record)
    }
    fn finish(self: Box<Self>) -> Result<Option<index_builder::IndexBuilder>, VcfError> {
        (*self).finish()
    }
}

/// Open a VCF/BCF writer based on file extension.
///
/// - `.vcf` → plain text VCF
/// - `.vcf.gz` → BGZF-compressed VCF with TBI index
/// - `.bcf` → BCF binary with CSI index
pub fn open_writer<W: Write + 'static>(
    inner: W,
    header: Arc<VcfHeader>,
    format: OutputFormat,
) -> Box<dyn VcfWrite> {
    match format {
        OutputFormat::Vcf => Box::new(writer::VcfWriter::new(inner, header)),
        OutputFormat::VcfGz => Box::new(writer::VcfWriter::bgzf(inner, header)),
        OutputFormat::Bcf => Box::new(bcf_writer::BcfWriter::new(inner, header, true)),
    }
}

/// Output format for [`open_writer`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputFormat {
    /// Plain text VCF.
    Vcf,
    /// BGZF-compressed VCF (.vcf.gz) with TBI index.
    VcfGz,
    /// BCF binary with CSI index.
    Bcf,
}

impl OutputFormat {
    // r[impl vcf_writer.output_formats]
    /// Detect format from a file path extension.
    /// Returns an error for unrecognized extensions.
    pub fn from_path(path: &std::path::Path) -> Result<Self, VcfError> {
        let name = path.to_str().unwrap_or("");
        if name.ends_with(".bcf") {
            Ok(Self::Bcf)
        } else if name.ends_with(".vcf.gz") || name.ends_with(".gz") {
            Ok(Self::VcfGz)
        } else if name.ends_with(".vcf") {
            Ok(Self::Vcf)
        } else {
            Err(VcfError::UnrecognizedFormat { path: name.to_string() })
        }
    }
}
