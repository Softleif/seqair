//! VCF/BCF error types. Typed fields only, no `String` payloads.

use crate::{
    io::BgzfError,
    vcf::header::{Number, ValueType, VcfVersion},
    vcf::writer::WriteError,
};
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

    // r[impl vcf_header.file_format.single_version+4]
    #[error("field {id}: Number={number} needs {needs}, but the header declares {version}")]
    NumberNeedsNewerVersion { id: SmolStr, number: Number, version: VcfVersion, needs: VcfVersion },

    #[error("FORMAT field {id} must not be Flag type")]
    FormatFlagNotAllowed { id: SmolStr },

    // r[impl vcf_record.sample_count]
    #[error("sample count mismatch: header declares {expected}, record has {actual}")]
    SampleCountMismatch { expected: usize, actual: usize },

    // r[impl vcf_record.format_gt_first]
    #[error("GT must be the first FORMAT key, but found at index {index}")]
    GtNotFirst { index: usize },

    #[error("too many contigs: a header may declare at most {max}")]
    TooManyContigs { max: usize },

    #[error("too many samples: a header may declare at most {max}")]
    TooManySamples { max: usize },

    #[error(
        "too many header fields: the FILTER/INFO/FORMAT string dictionary is limited to {max} \
         entries"
    )]
    TooManyFields { max: usize },

    #[error(transparent)]
    Key(#[from] KeyError),

    // r[impl vcf_header.parse.limits]
    #[error("header text is {len} bytes, more than the {max} a header may hold")]
    TextTooLarge { len: usize, max: usize },

    // r[impl vcf_header.parse.version]
    #[error("the first header line must be ##fileformat=VCFvX.Y")]
    MissingFileFormat,

    #[error("malformed ##fileformat value {value}: expected VCFvMAJOR.MINOR")]
    MalformedVersion { value: SmolStr },

    #[error("{version} is not a VCF 4 version")]
    UnsupportedVersion { version: VcfVersion },

    // r[impl vcf_header.parse]
    #[error("header line {line}: {reason}")]
    MalformedLine { line: usize, reason: MetaLineError },

    #[error("header line {line}: ##{key} definition has no ID")]
    MissingId { line: usize, key: SmolStr },

    // r[impl vcf_header.parse.defs]
    #[error("field {id}: unknown Number={value}")]
    UnknownNumber { id: SmolStr, value: SmolStr },

    #[error("field {id}: unknown Type={value}")]
    UnknownType { id: SmolStr, value: SmolStr },

    #[error("contig {name}: length={value} is not a non-negative integer")]
    InvalidContigLength { name: SmolStr, value: SmolStr },

    // r[impl vcf_header.parse.dictionary]
    #[error("{id}: IDX={value} is not a dictionary index")]
    InvalidIdx { id: SmolStr, value: SmolStr },

    #[error("IDX={idx} is claimed by both {first} and {second}")]
    ConflictingIdx { idx: u32, first: SmolStr, second: SmolStr },

    #[error("{id} is given IDX={first} and IDX={second}")]
    InconsistentIdx { id: SmolStr, first: u32, second: u32 },

    #[error("the header has no #CHROM line")]
    MissingHeaderLine,

    #[error("#CHROM line: column {column} is missing or misnamed")]
    MalformedHeaderLine { column: usize },

    #[error("#CHROM line: sample column {column} is empty")]
    EmptySampleName { column: usize },

    #[error("header line {line} follows the #CHROM line")]
    TrailingData { line: usize },
}

/// Why a structured `##key=<…>` line could not be parsed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum MetaLineError {
    #[error("an attribute has no key")]
    MissingKey,
    #[error("an attribute has no '='")]
    MissingEquals,
    #[error("a quoted value has no closing '\"'")]
    UnterminatedQuote,
    #[error("a bracketed value has no closing ']'")]
    UnterminatedBracket,
    #[error("the line has no closing '>'")]
    Unterminated,
    #[error("a value is followed by something other than ',' or '>'")]
    JunkAfterValue,
}

/// Which section of a header or record a key names a field of.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
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

// r[impl vcf_values.keys+2]
// r[impl vcf_reader.errors+2]
/// Why a key could not be resolved from a header.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum KeyError {
    #[error("{kind} field {name} is not declared in the header")]
    UnknownField { name: SmolStr, kind: KeyKind },

    #[error("FILTER {name} is not declared in the header")]
    UnknownFilter { name: SmolStr },

    #[error("contig {name} is not declared in the header")]
    UnknownContig { name: SmolStr },

    #[error("sample {name} is not in the header")]
    UnknownSample { name: SmolStr },

    #[error("field {name} is declared Type={declared}, not {requested}")]
    TypeMismatch { name: SmolStr, declared: ValueType, requested: ValueType },

    #[error("field {name} is declared Number={declared}, which the requested key cannot hold")]
    NumberMismatch { name: SmolStr, declared: Number },

    #[error("field {name} is not the FORMAT GT field, the only one a genotype key names")]
    NotGenotype { name: SmolStr },
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

// r[impl record_encoder.record_id]
/// Why a record ID is not one (VCF §1.6.1: "no whitespace or semicolons permitted").
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RecordIdProblem {
    /// An empty identifier, which would read back as no identifier at all.
    Empty,
    /// `.`, the missing value: a reader cannot tell it from no ID.
    Missing,
    /// It contains whitespace.
    Whitespace,
    /// It contains `;`, the separator between identifiers.
    Semicolon,
}

impl std::fmt::Display for RecordIdProblem {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Empty => "it is empty",
            Self::Missing => "`.` is the missing value",
            Self::Whitespace => "it contains whitespace",
            Self::Semicolon => "it contains `;`, the separator",
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

    #[error("{name} was resolved from another header than the one being written")]
    ForeignKey { name: SmolStr },

    // r[impl record_encoder.count_check+3]
    #[error(
        "field {field} is declared Number={number}, so it takes {expected} values, not {found}"
    )]
    CountMismatch { field: SmolStr, number: Number, expected: usize, found: usize },

    #[error("write_header() must be called before write_record()")]
    HeaderNotWritten,

    #[error("FORMAT fields were written but begin_samples() was not called")]
    FormatDataWithoutSamples,

    #[error("FORMAT field has {got} values but header declares {expected} samples")]
    SampleCountMismatch { expected: usize, got: usize },

    // r[impl record_encoder.index_kind]
    #[error("BCF is indexed by CSI only; a TBI index was asked for")]
    TbiForBcf,

    // r[impl record_encoder.record_id]
    #[error("record ID {id:?} is not an identifier: {problem}")]
    InvalidRecordId { id: SmolStr, problem: RecordIdProblem },

    // r[impl record_encoder.record_id]
    #[error("record ID {id} is given twice")]
    DuplicateRecordId { id: SmolStr },

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

    // r[impl vcf_passthrough.copy]
    #[error("the writer's header is not the output header of the field map")]
    ForeignFieldMap,

    // r[impl vcf_passthrough.copy]
    #[error("the record was not read with the input header of the field map")]
    ForeignRecord,

    // r[impl vcf_passthrough.copy]
    #[error("a record on {contig} at the telomere (POS 0) cannot be written")]
    TelomereRecord { contig: SmolStr },

    // r[impl vcf_passthrough.undeclared]
    #[error("{kind} {name} is declared by neither the input nor the output header")]
    UndeclaredField { kind: KeyKind, name: SmolStr },

    // r[impl vcf_passthrough.field_map.contigs]
    #[error("contig {name} is not declared by the output header")]
    UndeclaredContig { name: SmolStr },

    #[error(transparent)]
    FieldMap(#[from] FieldMapError),

    /// A copied record's value could not be read.
    #[error(transparent)]
    Read(Box<crate::vcf::reader::VcfReadError>),
}

impl From<crate::vcf::reader::VcfReadError> for VcfError {
    fn from(e: crate::vcf::reader::VcfReadError) -> Self {
        Self::Read(Box::new(e))
    }
}

// r[impl vcf_passthrough.field_map]
/// Why an input header's fields cannot be mapped to an output header's.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum FieldMapError {
    // r[impl vcf_passthrough.field_map.compatible]
    #[error(
        "{kind} field {name} is Type={input_type},Number={input_number} in the input but \
         Type={output_type},Number={output_number} in the output"
    )]
    Incompatible {
        kind: KeyKind,
        name: SmolStr,
        input_type: ValueType,
        input_number: Number,
        output_type: ValueType,
        output_number: Number,
    },

    // r[impl vcf_passthrough.field_map.contigs]
    #[error("contig {name} of the input is not declared by the output header")]
    MissingContig { name: SmolStr },

    // r[impl vcf_passthrough.field_map.samples]
    #[error("sample column {index} is {input:?} in the input but {output:?} in the output")]
    SampleMismatch { index: usize, input: Option<SmolStr>, output: Option<SmolStr> },
}
