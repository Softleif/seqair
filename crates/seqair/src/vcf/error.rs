//! VCF/BCF error types. Typed fields only, no `String` payloads.

use crate::{bam::bgzf::BgzfError, vcf::writer::WriteError};
use seqair_types::SmolStr;

// r[impl vcf_header.builder]
// r[impl vcf_header.no_duplicates]
#[derive(Debug, thiserror::Error)]
pub enum VcfHeaderError {
    #[error("duplicate contig: {name}")]
    DuplicateContig { name: SmolStr },

    #[error("duplicate INFO field: {id}")]
    DuplicateInfo { id: SmolStr },

    #[error("duplicate FORMAT field: {id}")]
    DuplicateFormat { id: SmolStr },

    #[error("duplicate FILTER: {id}")]
    DuplicateFilter { id: SmolStr },

    #[error("duplicate sample: {name}")]
    DuplicateSample { name: SmolStr },

    #[error("contig not declared in header: {name}")]
    MissingContig { name: SmolStr },

    #[error("INFO field not declared in header: {id}")]
    MissingInfo { id: SmolStr },

    #[error("FORMAT field not declared in header: {id}")]
    MissingFormat { id: SmolStr },

    #[error("FILTER not declared in header: {id}")]
    MissingFilter { id: SmolStr },

    #[error("Flag INFO field {id} must have Number=0")]
    FlagNumberMismatch { id: SmolStr },

    #[error("FORMAT field {id} must not be Flag type")]
    FormatFlagNotAllowed { id: SmolStr },

    // r[impl vcf_record.sample_count]
    #[error("sample count mismatch: header declares {expected}, record has {actual}")]
    SampleCountMismatch { expected: usize, actual: usize },

    // r[impl vcf_record.format_gt_first]
    #[error("GT must be the first FORMAT key, but found at index {index}")]
    GtNotFirst { index: usize },
}

// r[impl vcf_record.alleles_typed]
#[derive(Debug, thiserror::Error)]
pub enum AllelesError {
    #[error("SNV alt base must differ from ref")]
    SnvAltEqualsRef,

    #[error("SNV alt bases must not be empty")]
    SnvEmpty,

    #[error("SNV alt bases contain duplicates")]
    SnvDuplicateAlt,

    #[error("insertion sequence must not be empty")]
    InsertionEmpty,

    #[error("deletion sequence must not be empty")]
    DeletionEmpty,
}

#[derive(Debug, thiserror::Error)]
pub enum VcfEncodeError {
    #[error("value type mismatch for field {field}: expected {expected}, got {got}")]
    TypeMismatch { field: SmolStr, expected: SmolStr, got: SmolStr },

    #[error("integer overflow for field {field}: value {value}")]
    IntegerOverflow { field: SmolStr, value: i64 },
}

#[derive(Debug, thiserror::Error)]
pub enum VcfError {
    #[error(transparent)]
    Bgzf(#[from] BgzfError),

    #[error(transparent)]
    Header(#[from] VcfHeaderError),

    #[error(transparent)]
    Encode(#[from] VcfEncodeError),

    #[error(transparent)]
    Alleles(#[from] AllelesError),

    #[error(transparent)]
    Io(#[from] std::io::Error),

    #[error(transparent)]
    Index(#[from] super::index_builder::IndexError),

    #[error("write_header() must be called before write_record()")]
    HeaderNotWritten,

    #[error("header text too large for BCF (exceeds u32::MAX)")]
    HeaderTooLarge,

    #[error("BCF record too large: {section} is {size} bytes (exceeds u32::MAX)")]
    RecordTooLarge { section: &'static str, size: usize },

    // r[impl bcf_encoder.checked_casts]
    #[error("value overflow: {field} value {value} exceeds {target_type} range")]
    ValueOverflow { field: &'static str, value: u64, target_type: &'static str },

    // r[impl vcf_writer.output_formats]
    #[error("unrecognized output format for path: {path}")]
    UnrecognizedFormat { path: String },

    #[error("failed to write field {field}")]
    FailedToWriteFormattedString { field: SmolStr, source: WriteError },
}
