//! VCF/BCF writing. Construct a [`VcfHeader`] via its builder, create
//! [`VcfRecord`]s with type-safe [`Alleles`], and write to VCF text or BCF binary.
//!
//! Use [`open_writer`] to create a writer based on file extension, or construct
//! [`VcfWriter`](writer::VcfWriter) / [`BcfWriter`](bcf_writer::BcfWriter) directly.

pub mod alleles;
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
        self.write_record(record)
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
    /// Detect format from a file path extension.
    pub fn from_path(path: &std::path::Path) -> Self {
        let name = path.to_str().unwrap_or("");
        if name.ends_with(".bcf") {
            Self::Bcf
        } else if name.ends_with(".vcf.gz") || name.ends_with(".gz") {
            Self::VcfGz
        } else {
            Self::Vcf
        }
    }
}
