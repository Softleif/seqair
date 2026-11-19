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
//!
//! # How a field lookup flows
//!
//! [`Info`] and the sample components own their section's bytes and its
//! index, a `SectionIndex` (INFO) or `FormatIndex` (FORMAT) in a `RefCell`
//! that `fill` resets in place; a [`RawRecord`]'s views use two the reader
//! keeps. A view ([`InfoRef`], [`SampleRef`], [`SamplesRef`]) borrows the
//! bytes and the index and fixes the layout once (`InfoSource` /
//! `FormatSource`: BCF or VCF text). [`InfoRef::get`] then
//!
//! 1. checks the key's header brand (`info::check_brand`);
//! 2. asks the index for the key's `DictId` (`SectionIndex::find`,
//!    `index.rs`): a hit is one table read; a miss scans on from where the
//!    last access stopped (`entries.rs`) and records how far it got
//!    (`Progress`: still scanning, done, or failed with the error every
//!    later access returns). FORMAT is scanned whole on first access
//!    (`FormatIndex::finish`), then looked up the same way;
//! 3. takes the value's bytes (`InfoSource::value`, a sample's cell for
//!    FORMAT), checks them as the key's marker reads them (`Valid::check`,
//!    the only step that fails), and reads the checked value
//!    (`ReadKind::read`, `access.rs`). [`SamplesRef::column`] checks every
//!    sample first (a BCF column at once), so its [`Column`] reads each
//!    without a check.
//!
//! No borrow of the index outlives the method that takes it (iteration and
//! a VCF text [`Column`] borrow again per item), so the `RefCell` never
//! sees two at once (`r[vcf_fields.lazy+4]`).
//!
//! # Where a key's type lives
//!
//! | when | what is decided | where |
//! |---|---|---|
//! | compile time | the result's shape, `V::Out<'a>` | the marker `V` of `InfoKey<V>` / `FormatKey<V>` |
//! | once per key (`file.info_key::<V>(name)`) | that `V` reads the declared `Type` and `Number`; the field's `DictId`, `Type`, whether it is GT, and its header | `check_kind`, kept in the key's `FieldId` |
//! | per record | only that the bytes are well formed | `V::check`, failing with a [`FieldError`]; `V::read` then cannot fail |
//!
//! An access never consults the header: [`Any`] decodes by the `Type` its
//! key carries. Generic iteration has no key, and reads each field's
//! definition from the file's field table (`keys.rs`) by id, which also
//! knows the id of `GT`.

mod access;
pub(crate) mod bcf_value;
mod entries;
mod error;
mod fields;
mod file;
mod full;
mod index;
pub(crate) use index::DictId;
pub(crate) use samples::DynColumn;
mod info;
mod keys;
mod locus;
mod query;
mod raw;
mod region;
mod region_index;
mod reserved;
pub(crate) mod rlen;
mod samples;
mod site;
mod stream;
mod text_keys;
pub(crate) mod text_value;
pub(crate) mod values;
mod vcf_stream;

pub use access::{Any, DecodeCx, DynValue, FormatAny, InfoAny, ReadKind, Stored};
pub use bcf_value::{BcfType, BcfValueError};
pub use error::{
    FieldError, KeyKind, OpenError, RecordError, RecordPart, RecordPlace, RecordSection, VcfFormat,
    VcfReadError,
};
pub use fields::{
    Fields, PlanContext, PlanError, ReadOption, SampleSel, SampleSelection, Sections,
};
pub use file::VcfFile;
pub use full::FullRecord;
pub use info::{Info, InfoRef};
pub use keys::{Keys, LenientKeys};
pub use locus::{Locus, VcfPos};
pub use query::Partition;
pub use raw::RawRecord;
pub use region::{Coordinate, Match, Region, RegionContig, RegionEnd, RegionError};
pub use samples::{Column, OnlySample, Sample, SampleRef, Samples, SamplesRef};
pub use site::{Allele, AllelesRef, FailedFilters, FilterStatus, FiltersRef, Site, SiteRef};
pub use stream::{IntoRecords, MAX_LINE_BYTES, MAX_RECORD_BYTES, Reader, ReaderBuilder, Records};
pub use text_value::{Offending, TextValueError};
pub use values::{GenotypeRef, GtAllele, MissingElement, OptValues, Values};
pub use vcf_stream::{StreamReaderBuilder, VcfStream};
