//! VCF/BCF error types. Typed fields only, no `String` payloads.

use crate::{io::BgzfError, vcf::writer::WriteError};
use seqair_types::SmolStr;
use std::path::PathBuf;

// r[impl vcf_header.builder]
// r[impl vcf_header.no_duplicates]
#[non_exhaustive]
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

    #[error("too many contigs (exceeds u32::MAX)")]
    TooManyContigs,

    #[error(
        "too many header fields: the FILTER/INFO/FORMAT string dictionary is limited to {max} \
         entries (so the per-record duplicate-field tracker can use a u64 bitset keyed by \
         dictionary index)"
    )]
    TooManyFields { max: usize },
}

// r[impl vcf_record.alleles_typed]
#[non_exhaustive]
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

// r[impl bcf_encoder.reserved_bands]
/// Which member of a BCF reserved band a value collided with.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReservedMarker {
    /// The MISSING marker — a reader renders the element as `.`.
    Missing,
    /// The `END_OF_VECTOR` marker — a reader truncates the row at this element.
    EndOfVector,
    /// Reserved by BCF2 for future use; not currently interpreted, but not data either.
    FutureUse,
}

impl std::fmt::Display for ReservedMarker {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Missing => "MISSING",
            Self::EndOfVector => "END_OF_VECTOR",
            Self::FutureUse => "reserved (future use)",
        })
    }
}

#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum VcfEncodeError {
    #[error("value type mismatch for field {field}: expected {expected}, got {got}")]
    TypeMismatch { field: SmolStr, expected: SmolStr, got: SmolStr },

    #[error("integer overflow for field {field}: value {value}")]
    IntegerOverflow { field: SmolStr, value: i64 },

    // r[impl bcf_encoder.reserved_rejected]
    #[error(
        "field {field}: BCF would read the integer {value} back as its {marker} marker, not as data"
    )]
    ReservedIntValue { field: SmolStr, value: i32, marker: ReservedMarker },

    // r[impl bcf_encoder.reserved_rejected]
    #[error(
        "field {field}: BCF would read the float bit pattern {bits:#010X} back as its {marker} \
         marker, not as data"
    )]
    ReservedFloatValue { field: SmolStr, bits: u32, marker: ReservedMarker },
}

#[non_exhaustive]
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
    Index(#[from] crate::io::IndexError),

    #[error("write_header() must be called before write_record()")]
    HeaderNotWritten,

    #[error("FORMAT fields were written but begin_samples() was not called")]
    FormatDataWithoutSamples,

    #[error("FORMAT field has {got} values but header declares {expected} samples")]
    SampleCountMismatch { expected: usize, got: usize },

    #[error(
        "mixed ploidy in GT field: first sample has ploidy {first_ploidy} but sample {mismatch_index} has ploidy {mismatch_ploidy}"
    )]
    MixedPloidy { first_ploidy: usize, mismatch_index: usize, mismatch_ploidy: usize },

    #[error("header text too large for BCF (exceeds u32::MAX)")]
    HeaderTooLarge,

    #[error("BCF record too large: {section} is {size} bytes (exceeds u32::MAX)")]
    RecordTooLarge { section: &'static str, size: usize },

    // r[impl bcf_encoder.checked_casts]
    #[error("value overflow: {field} value {value} exceeds {target_type} range")]
    ValueOverflow { field: &'static str, value: u64, target_type: &'static str },

    // r[impl vcf_writer.output_formats]
    #[error("unrecognized output format for path: {}", path.display())]
    UnrecognizedFormat { path: PathBuf },

    #[error("failed to write float field")]
    FailedToWriteFormattedString { source: WriteError },
}
