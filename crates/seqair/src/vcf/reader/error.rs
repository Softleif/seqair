//! What reading a VCF or BCF file can fail with (`r[vcf_reader.errors+2]`).

use super::bcf_value::BcfValueError;
use super::fields::PlanError;
use super::file::FileShared;
use super::locus::VcfPos;
use super::text_value::TextValueError;
use super::values::MissingElement;
use crate::bam::bgzf::{BgzfError, VirtualOffset};
use crate::bam::csi_index::CsiError;
pub use crate::vcf::error::KeyKind;
use crate::vcf::error::VcfHeaderError;
use seqair_types::SmolStr;
use std::path::{Path, PathBuf};
use std::str::Utf8Error;

// r[impl vcf_reader.errors+2]
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

    #[error("the header text is not UTF-8")]
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

    // r[impl vcf_reader.errors+2]
    /// A record is malformed: `at` says where it is, `source` what is wrong.
    #[error("malformed record {at}")]
    Record { at: Box<RecordPlace>, source: RecordError },

    // r[impl vcf_values.header_brand+2]
    #[error("{name} was resolved from another header than this record's")]
    ForeignKey { name: SmolStr },

    #[error("Fields::needs did not ask for the {section} of records")]
    NotRequested { section: RecordSection },
}

/// Where a malformed record is, as far as the reader knows.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RecordPlace {
    /// The record's number in the stream, from 1.
    pub number: u64,
    /// VCF text read from its start: the record's line, from 1.
    pub line: Option<u64>,
    /// BGZF input: where the record starts.
    pub offset: Option<VirtualOffset>,
    /// The file's path, unless it was read from memory.
    pub path: Option<PathBuf>,
    /// The record's CHROM and POS, once they decoded.
    pub locus: Option<(SmolStr, VcfPos)>,
}

impl std::fmt::Display for RecordPlace {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.number)?;
        if let Some(line) = self.line {
            write!(f, " (line {line})")?;
        }
        if let Some((chrom, pos)) = &self.locus {
            write!(f, " at {chrom}:{pos}")?;
        }
        if let Some(path) = &self.path {
            write!(f, " in {}", path.display())?;
        }
        if let Some(o) = self.offset {
            write!(f, " (virtual offset {}:{})", o.block_offset(), o.within_block())?;
        }
        Ok(())
    }
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

/// A part of a record, to say where a value was malformed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RecordPart {
    Chrom,
    Pos,
    Id,
    Alleles,
    Qual,
    Filter,
    Info,
    Format,
}

impl std::fmt::Display for RecordPart {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Chrom => "CHROM",
            Self::Pos => "POS",
            Self::Id => "ID",
            Self::Alleles => "alleles",
            Self::Qual => "QUAL",
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
    ContigIndexOutOfRange { index: i32 },

    #[error("POS {pos} (0-based) is before the telomere or past 2^31-2")]
    InvalidPosition { pos: i64 },

    #[error("POS {pos} (0-based) with rlen {rlen} ends outside 0..=2^31-1")]
    InvalidSpan { pos: i64, rlen: i64 },

    #[error("the record has {found} samples, the header {expected}")]
    SampleCount { found: u32, expected: usize },

    // r[impl vcf_reader.limits]
    #[error("the line is {len} bytes or more, more than the {max} a line may hold")]
    LineTooLong { len: u64, max: u64 },

    #[error("the line has {found} TAB-separated columns, fewer than VCF's 8")]
    TooFewColumns { found: usize },

    // r[impl vcf_reader.limits]
    #[error("the record has {n} alleles, more than BCF's 16-bit count holds")]
    TooManyAlleles { n: usize },

    #[error("sample {sample} has {found} fields or more; FORMAT names {keys}")]
    SampleFields { sample: usize, found: usize, keys: usize },

    #[error("malformed {part}")]
    Text { part: RecordPart, source: TextValueError },

    // r[impl vcf_reader.field_table]
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

/// Where a record is, as cheaply as a reader can say it per record: its
/// ordinal, its line (VCF text from the start), its virtual offset (BGZF),
/// and its contig and 0-based position once they decoded. The names are
/// looked up only for an error ([`Self::error`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub(crate) struct RecordAt {
    pub(crate) ordinal: u64,
    pub(crate) line: Option<u64>,
    pub(crate) offset: Option<VirtualOffset>,
    pub(crate) contig: Option<u32>,
    /// 0-based; −1 is the telomere.
    pub(crate) pos: Option<i64>,
}

impl RecordAt {
    pub(crate) fn error(self, file: &FileShared, source: RecordError) -> VcfReadError {
        let locus = self.contig.and_then(|tid| file.contigs.name(tid)).and_then(|chrom| {
            let pos = VcfPos::from_one_based(self.pos?.saturating_add(1))?;
            Some((chrom, pos))
        });
        let at = RecordPlace {
            number: self.ordinal.saturating_add(1),
            line: self.line,
            offset: self.offset,
            path: file.path().map(Path::to_path_buf),
            locus,
        };
        VcfReadError::Record { at: Box::new(at), source }
    }
}
