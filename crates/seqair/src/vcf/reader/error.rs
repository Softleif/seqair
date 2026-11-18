//! What reading a VCF or BCF file can fail with (`r[vcf_reader.errors]`).

use super::bcf_value::BcfValueError;
use super::fields::PlanError;
use super::text_value::TextValueError;
use super::values::MissingElement;
use crate::bam::bgzf::{BgzfError, VirtualOffset};
use crate::bam::csi_index::CsiError;
use crate::vcf::error::VcfHeaderError;
use seqair_types::SmolStr;
use std::path::PathBuf;
use std::str::Utf8Error;

// r[impl vcf_reader.errors]
/// Why opening or reading a VCF/BCF file failed.
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum VcfReadError {
    #[error("could not open {path}")]
    Open { path: PathBuf, source: std::io::Error },

    #[error(transparent)]
    Bgzf(#[from] BgzfError),

    #[error(transparent)]
    Header(#[from] VcfHeaderError),

    #[error("the BCF header text is not UTF-8")]
    HeaderNotUtf8 { source: Utf8Error },

    #[error(transparent)]
    Index(#[from] CsiError),

    #[error(transparent)]
    Plan(#[from] PlanError),

    // r[impl vcf_reader.formats]
    #[error("the input is gzip-compressed but not BGZF; recompress it with `bgzip`")]
    NotBgzf,

    #[error("BCF 2.1 is not supported; convert the file with `bcftools view -Ob`")]
    Bcf21,

    #[error("uncompressed BCF is not supported; compress the file with `bcftools view -Ob`")]
    UncompressedBcf,

    #[error("the input is neither VCF nor BCF")]
    UnknownFormat,

    #[error("reading {format} is not supported yet")]
    Unsupported { format: VcfFormat },

    /// A record is malformed: `ordinal` counts records from 0, `offset` is
    /// where it starts in a BGZF file.
    #[error("record {ordinal}{}: {source}", at(*.offset))]
    Record { ordinal: u64, offset: Option<VirtualOffset>, source: RecordError },

    // r[impl vcf_values.header_brand]
    #[error("{name} was resolved from another header than this record's")]
    ForeignKey { name: SmolStr },

    #[error("Fields::needs did not ask for the {section} of records")]
    NotRequested { section: RecordSection },
}

fn at(offset: Option<VirtualOffset>) -> String {
    offset.map_or_else(String::new, |o| format!(" (at {}:{})", o.block_offset(), o.within_block()))
}

/// A file format the reader detects.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VcfFormat {
    /// Plain-text VCF.
    Vcf,
    /// BGZF-compressed VCF.
    VcfGz,
    /// BGZF-compressed BCF 2.2.
    Bcf,
}

impl std::fmt::Display for VcfFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Vcf => "VCF",
            Self::VcfGz => "BGZF-compressed VCF",
            Self::Bcf => "BCF",
        })
    }
}

/// A part of a record a [`Fields`](super::Fields) type asks for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RecordSection {
    /// ID, REF, ALT, QUAL and FILTER.
    Site,
    Info,
    /// The FORMAT values of one sample (by its column) or all of them.
    Sample {
        index: usize,
    },
    Samples,
}

impl std::fmt::Display for RecordSection {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Site => f.write_str("site columns"),
            Self::Info => f.write_str("INFO column"),
            Self::Sample { index } => write!(f, "sample {index}"),
            Self::Samples => f.write_str("samples"),
        }
    }
}

/// A part of a BCF record, to say where a value was malformed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RecordPart {
    Id,
    Alleles,
    Filter,
    Info,
    Format,
}

impl std::fmt::Display for RecordPart {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Id => "ID",
            Self::Alleles => "alleles",
            Self::Filter => "FILTER",
            Self::Info => "INFO",
            Self::Format => "FORMAT",
        })
    }
}

/// Which dictionary a record's key should be in.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KeyKind {
    Filter,
    Info,
    Format,
}

impl std::fmt::Display for KeyKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Filter => "FILTER",
            Self::Info => "INFO",
            Self::Format => "FORMAT",
        })
    }
}

/// What is wrong with one record. `Clone`, so a lazily built index that
/// failed returns the same error on every access (`r[vcf_fields.lazy]`).
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum RecordError {
    // r[impl vcf_reader.limits]
    #[error("the record is {size} bytes, more than the {max} a record may hold")]
    TooLarge { size: u64, max: u64 },

    #[error("the shared section is {len} bytes, shorter than its 24 fixed bytes")]
    SharedTooShort { len: u32 },

    #[error("the input ends inside the record")]
    Truncated,

    #[error("contig index {index} is not in the header's contig dictionary")]
    UnknownContig { index: i32 },

    #[error("POS {pos} (0-based) is before the telomere")]
    InvalidPosition { pos: i32 },

    #[error("POS {pos} (0-based) with rlen {rlen} ends outside 1..=2^31-1")]
    InvalidSpan { pos: i32, rlen: i32 },

    #[error("the record has {found} samples, the header {expected}")]
    SampleCount { found: u32, expected: usize },

    // r[impl vcf_reader.undeclared]
    #[error("dictionary index {index} names no {kind} field of the header")]
    UnknownKey { index: i32, kind: KeyKind },

    #[error("malformed {part}")]
    Value { part: RecordPart, source: BcfValueError },

    #[error("{part} is not UTF-8")]
    Utf8 { part: RecordPart, source: Utf8Error },

    #[error("{len} bytes are left over after the {part}")]
    TrailingBytes { part: RecordPart, len: usize },

    #[error("field {name}")]
    Field { name: SmolStr, source: FieldError },
}

/// What is wrong with one field's value.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum FieldError {
    #[error(transparent)]
    MissingElement(#[from] MissingElement),

    #[error("holds {len} values where its declaration allows one")]
    NotScalar { len: usize },

    // r[impl vcf_reader.limits]
    #[error("genotype ploidy {ploidy} is above {max}")]
    Ploidy { ploidy: usize, max: usize },

    #[error(transparent)]
    Bcf(#[from] BcfValueError),

    #[error(transparent)]
    Text(#[from] TextValueError),
}

/// Where a record is: its ordinal and, in a BGZF file, its virtual offset.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub(crate) struct RecordAt {
    pub(crate) ordinal: u64,
    pub(crate) offset: Option<VirtualOffset>,
}

impl RecordAt {
    pub(crate) fn error(self, source: RecordError) -> VcfReadError {
        VcfReadError::Record { ordinal: self.ordinal, offset: self.offset, source }
    }
}
