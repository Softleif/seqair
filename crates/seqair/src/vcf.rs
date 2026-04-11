//! VCF/BCF writing with type-safe alleles and zero-allocation encoding.
//!
//! Use [`Writer`] with [`OutputFormat`] for all VCF/BCF output:
//! - `.vcf` → plain text VCF
//! - `.vcf.gz` → BGZF-compressed VCF with automatic TBI index co-production
//! - `.bcf` → BCF binary with automatic CSI index co-production
//!
//! Records are encoded through a typestate [`RecordEncoder`] that enforces
//! the correct field ordering at compile time.
//!
//! A single [`Writer`] type handles all output formats (VCF, VCF.gz, BCF).
//! Records are encoded through a typestate chain that is enforced at compile time:
//! `Begun` → `Filtered` → `WithSamples` → `emit()`.
//!
//! # Writing VCF records
//!
//! ```
//! use seqair_types::{Base, Pos1};
//! use seqair::vcf::{
//!     Alleles, ContigDef, Genotype, Number, OutputFormat, ValueType, VcfHeader, Writer, FormatGt, FormatInt, InfoInt,
//!     record_encoder::{FormatFieldDef, Gt, InfoFieldDef, Scalar},
//! };
//! use std::sync::Arc;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! // 1. Build header — typed field handles are resolved once at setup, not per record.
//! // In production use VcfHeaderBuilder::from_bam_header() to copy contig info.
//! let mut builder = VcfHeader::builder();
//! let chr1 = builder.register_contig("chr1", ContigDef { length: Some(248_956_422) })?;
//! let dp_info: InfoInt = builder.register_info(
//!     &InfoFieldDef::<Scalar<i32>>::new("DP", Number::Count(1), ValueType::Integer, "Total read depth")
//! )?;
//! let gt_fmt: FormatGt = builder.register_format(
//!     &FormatFieldDef::<Gt>::new("GT", Number::Count(1), ValueType::String, "Genotype")
//! )?;
//! let dp_fmt: FormatInt = builder.register_format(
//!     &FormatFieldDef::<Scalar<i32>>::new("DP", Number::Count(1), ValueType::Integer, "Sample depth")
//! )?;
//! let header = Arc::new(builder.add_sample("sample1")?.build()?);
//!
//! // 2. Write to an in-memory buffer (or any `impl Write`)
//! let mut output = Vec::new();
//! let writer = Writer::new(&mut output, OutputFormat::Vcf);
//! let mut writer = writer.write_header(&header)?;
//!
//! // 3. Encode one record — the typestate chain is enforced at compile time
//! let alleles = Alleles::snv(Base::A, Base::T)?;
//! let enc = writer.begin_record(&chr1, Pos1::new(12345).unwrap(), &alleles, Some(30.0))?;
//! let mut enc = enc.filter_pass();    // Begun → Filtered
//! dp_info.encode(&mut enc, 50);
//! let mut enc = enc.begin_samples(); // Filtered → WithSamples
//! gt_fmt.encode(&mut enc, &[Genotype::unphased(0, 1)]);
//! dp_fmt.encode(&mut enc, &[45]);
//! enc.emit()?;
//!
//! writer.finish()?;
//! # let vcf = String::from_utf8(output)?;
//! # assert!(vcf.contains("chr1\t12345\t.\tA\tT\t30\tPASS\tDP=50"));
//! # Ok(())
//! # }
//! ```
//!
//! # Writing BCF (binary VCF)
//!
//! Switch [`OutputFormat::Vcf`] to [`OutputFormat::Bcf`] — the encoding API is
//! identical. Pre-resolved handles perform direct BCF encoding with no per-record
//! allocations or string dictionary lookups.
//!
//! ```
//! use seqair_types::{Base, Pos1};
//! use seqair::vcf::{
//!     Alleles, ContigDef, Genotype, Number, OutputFormat, ValueType, VcfHeader, Writer, FormatGt, InfoInt,
//!     record_encoder::{FormatFieldDef, Gt, InfoFieldDef, Scalar}
//! };
//! use std::sync::Arc;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let mut builder = VcfHeader::builder();
//! let contig = builder.register_contig("chr1", ContigDef { length: Some(1000) })?;
//! let dp: InfoInt = builder.register_info(
//!     &InfoFieldDef::<Scalar<i32>>::new("DP", Number::Count(1), ValueType::Integer, "Depth")
//! )?;
//! let gt: FormatGt = builder.register_format(
//!     &FormatFieldDef::<Gt>::new("GT", Number::Count(1), ValueType::String, "Genotype")
//! )?;
//! let header = Arc::new(builder.add_sample("sample1")?.build()?);
//!
//! let mut buf = Vec::new();
//! let writer = Writer::new(&mut buf, OutputFormat::Bcf);
//! let mut writer = writer.write_header(&header)?;
//!
//! let alleles = Alleles::snv(Base::A, Base::T)?;
//! let enc = writer.begin_record(&contig, Pos1::new(100).unwrap(), &alleles, Some(30.0))?;
//! let mut enc = enc.filter_pass();
//! dp.encode(&mut enc, 50);
//! let mut enc = enc.begin_samples();
//! gt.encode(&mut enc, &[Genotype::unphased(0, 1)]);
//! enc.emit()?;
//!
//! writer.finish()?;
//! # assert!(!buf.is_empty());
//! # Ok(())
//! # }
//! ```

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
