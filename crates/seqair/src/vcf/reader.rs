// r[impl vcf_fields.guide]
#![doc = include_str!("reading.md")]
//!
//! # Layout
//!
//! [`VcfFile`] opens a file once; [`VcfFile::reader`] plans a [`Fields`] type
//! against its header and opens a [`Reader`]. A reader hands each record to
//! `F::new` or `F::fill` as a [`RawRecord`], whose views ([`SiteRef`], [`InfoRef`],
//! [`SampleRef`], [`SamplesRef`]) decode values through the typed keys of
//! [`ReadKind`]. The components ([`Locus`], [`Site`], [`Info`],
//! [`OnlySample`], [`Sample`], [`Samples`], [`FullRecord`]) own their bytes.
//! The value decoders underneath know no header: BCF typed values
//! (`bcf_value`), VCF text values (`text_value`), and the views both share
//! ([`OptValues`], [`Values`], [`GenotypeRef`]).

mod access;
pub(crate) mod bcf_value;
mod entries;
mod error;
mod fields;
mod file;
mod full;
mod index;
mod info;
mod keys;
mod locus;
mod raw;
mod rlen;
mod samples;
mod site;
mod stream;
mod text_keys;
pub(crate) mod text_value;
pub(crate) mod values;

pub use access::{DecodeCx, DynValue, ReadKind, Stored};
pub use bcf_value::{BcfType, BcfValueError};
pub use error::{
    FieldError, KeyKind, RecordError, RecordPart, RecordPlace, RecordSection, VcfFormat,
    VcfReadError,
};
pub use fields::{
    Fields, PlanContext, PlanError, ReadOption, SampleSel, SampleSelection, Sections,
};
pub use file::VcfFile;
pub use full::FullRecord;
pub use info::{Info, InfoRef};
pub use locus::{Locus, VcfPos};
pub use raw::RawRecord;
pub use samples::{Column, OnlySample, Sample, SampleRef, Samples, SamplesRef};
pub use site::{Allele, AllelesRef, FailedFilters, FilterStatus, FiltersRef, Site, SiteRef};
pub use stream::{IntoRecords, MAX_LINE_BYTES, MAX_RECORD_BYTES, Reader, ReaderBuilder, Records};
pub use text_value::{Offending, TextValueError};
pub use values::{GenotypeRef, GtAllele, MissingElement, OptValues, Values};
