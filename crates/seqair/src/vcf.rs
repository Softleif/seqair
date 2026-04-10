//! VCF/BCF writing with type-safe alleles and zero-allocation encoding.
//!
//! Use [`Writer`] with [`OutputFormat`] for all VCF/BCF output:
//! - `.vcf` → plain text VCF
//! - `.vcf.gz` → BGZF-compressed VCF with automatic TBI index co-production
//! - `.bcf` → BCF binary with automatic CSI index co-production
//!
//! Records are encoded through a typestate [`RecordEncoder`] that enforces
//! the correct field ordering at compile time.

pub mod alleles;
pub(crate) mod bcf_encoding;
pub(crate) mod encoder;
pub mod error;
pub mod header;
pub mod index_builder;
pub mod record;
pub mod record_encoder;
pub mod unified;
pub(crate) mod writer;

pub use alleles::Alleles;
pub use error::{AllelesError, VcfEncodeError, VcfError, VcfHeaderError};
pub use header::{
    ContigDef, FilterDef, FormatDef, InfoDef, Number, ValueType, VcfHeader, VcfHeaderBuilder,
};
pub use record::Genotype;
pub use record_encoder::{
    Arr, ContigId, EncodeFormat, EncodeInfo, FieldDescription, FieldId, FilterFieldDef, FilterId,
    Flag, FormatEncoder, FormatFieldDef, FormatFloat, FormatGt, FormatInt, FormatKey, Gt,
    InfoEncoder, InfoFieldDef, InfoFlag, InfoFloat, InfoFloats, InfoInt, InfoIntOpts, InfoInts,
    InfoKey, InfoString, OptArr, Scalar, Str,
};
pub use unified::{Begun, Filtered, Ready, RecordEncoder, Unstarted, WithSamples, Writer};

/// Output format for [`Writer`].
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
