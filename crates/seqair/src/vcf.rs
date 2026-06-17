//! VCF/BCF writing with type-safe alleles and zero-allocation encoding.
//!
//! Use [`Writer`] with [`OutputFormat`] for all VCF/BCF output:
//! - `.vcf` → plain text VCF (no index)
//! - `.vcf.gz` → BGZF-compressed VCF; [`Writer::finish`] returns a CSI [`CoordinateIndex`]
//! - `.bcf` → BCF binary; [`Writer::finish`] returns a CSI [`CoordinateIndex`]
//!
//! For compressed output, [`Writer::finish`] hands back a [`CoordinateIndex`]
//! alongside the inner writer; write it next to the data file (e.g.
//! `calls.bcf.csi`) with [`CoordinateIndex::write`] or
//! [`write_to_path`](CoordinateIndex::write_to_path).
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
//! let mut builder = builder.infos();
//! let dp_info: InfoInt = builder.register_info(
//!     &InfoFieldDef::<Scalar<i32>>::new("DP", Number::Count(1), ValueType::Integer, "Total read depth")
//! )?;
//! let mut builder = builder.formats();
//! let gt_fmt: FormatGt = builder.register_format(
//!     &FormatFieldDef::<Gt>::new("GT", Number::Count(1), ValueType::String, "Genotype")
//! )?;
//! let dp_fmt: FormatInt = builder.register_format(
//!     &FormatFieldDef::<Scalar<i32>>::new("DP", Number::Count(1), ValueType::Integer, "Sample depth")
//! )?;
//! let mut builder = builder.samples();
//! builder.add_sample("sample1")?;
//! let header = Arc::new(builder.build()?);
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
//! gt_fmt.encode(&mut enc, &[Genotype::unphased(0, 1)])?;
//! dp_fmt.encode(&mut enc, &[45])?;
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
//! let mut builder = builder.infos();
//! let dp: InfoInt = builder.register_info(
//!     &InfoFieldDef::<Scalar<i32>>::new("DP", Number::Count(1), ValueType::Integer, "Depth")
//! )?;
//! let mut builder = builder.formats();
//! let gt: FormatGt = builder.register_format(
//!     &FormatFieldDef::<Gt>::new("GT", Number::Count(1), ValueType::String, "Genotype")
//! )?;
//! let mut builder = builder.samples();
//! builder.add_sample("sample1")?;
//! let header = Arc::new(builder.build()?);
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
//! gt.encode(&mut enc, &[Genotype::unphased(0, 1)])?;
//! enc.emit()?;
//!
//! // For compressed output, finish() returns a CoordinateIndex. Write it next
//! // to the data file as `<path>.csi`; here we serialize to an in-memory buffer.
//! let (_, index) = writer.finish()?;
//! let mut csi = Vec::new();
//! index.expect("BCF output is indexed").write(&mut csi)?;
//! # assert!(!buf.is_empty());
//! # assert!(!csi.is_empty());
//! # Ok(())
//! # }
//! ```
//!
//! # Collecting field defs and resolved keys in a struct
//!
//! Co-locate each field's header definition with the resolved, typed key your
//! application encodes through. The typed key (`InfoInt`, `InfoFloats`, …)
//! enforces the value type at the call site, so there is no need for a custom
//! encoding trait.
//!
//! ```rust
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! use seqair_types::{Base, Pos1};
//! use seqair::vcf;
//!
//! // header metadata for the field, defined once
//! const DP_DEF: vcf::InfoFieldDef<vcf::Scalar<i32>> = vcf::InfoFieldDef::new(
//!     "DP", vcf::Number::Count(1), vcf::ValueType::Integer, "Combined depth"
//! );
//!
//! // all resolved keys your app will encode through
//! struct MyInfoFields {
//!   depth: vcf::InfoInt,
//! }
//!
//! // default header setup
//! let mut builder = vcf::VcfHeader::builder();
//! let contig =
//!     builder.register_contig("chr1", vcf::ContigDef { length: Some(1000) })?;
//! let mut builder = builder.infos();
//! // register our field, keeping the typed key it resolves to
//! let depth_key = builder.register_info(&DP_DEF)?;
//! let header = builder.build()?;
//! let my_fields = MyInfoFields { depth: depth_key };
//!
//! let mut buf = Vec::new();
//! let writer = vcf::Writer::new(&mut buf, vcf::OutputFormat::Bcf);
//!
//! // start by writing header
//! let mut writer = writer.write_header(&header)?;
//!
//! // new record
//! let mut enc = writer.begin_record(&contig,
//!     Pos1::new(1).ok_or("invalid position")?,
//!     &vcf::Alleles::snv(Base::A, Base::T).unwrap(),
//!     None
//! )?;
//! // PASS
//! let mut enc = enc.filter_pass();
//! // add the info field through its typed key (only i32 is accepted)
//! my_fields.depth.encode(&mut enc, 30);
//! // write record
//! enc.emit()?;
//!
//! writer.finish()?;
//! # assert!(!buf.is_empty());
//! # Ok(()) }
//! ```

pub mod alleles;
pub(crate) mod bcf_encoding;
pub(crate) mod encoder;
pub mod error;
pub mod header;
pub mod record;
pub mod record_encoder;
pub mod unified;
pub(crate) mod writer;

pub use alleles::Alleles;
pub use error::{AllelesError, VcfEncodeError, VcfError, VcfHeaderError};
pub use header::{
    ContigDef, Contigs, FilterDef, Filters, FormatDef, Formats, FromBamHeader, InfoDef, Infos,
    Number, Samples, ValueType, VcfHeader, VcfHeaderBuilder,
};
pub use record::Genotype;
pub use record_encoder::{
    Arr, ContigId, FieldDescription, FieldId, FilterFieldDef, FilterId, Flag, FormatEncoder,
    FormatFieldDef, FormatFloat, FormatFloats, FormatGt, FormatInt, FormatInts, FormatKey,
    FormatString, Gt, InfoEncoder, InfoFieldDef, InfoFlag, InfoFloat, InfoFloats, InfoInt,
    InfoIntOpts, InfoInts, InfoKey, InfoString, OptArr, Scalar, Str,
};
pub use unified::{
    Begun, CoordinateIndex, Filtered, Ready, RecordEncoder, Unstarted, WithSamples, Writer,
};

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
            Err(VcfError::UnrecognizedFormat { path: path.to_path_buf() })
        }
    }
}
