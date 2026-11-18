//! What a reader decodes: the [`Fields`] trait, the [`Sections`] a type
//! needs, and the options it plans with (`docs/spec/2-vcf-2-fields.md`).

use super::error::VcfReadError;
use super::file::FileShared;
use super::raw::RawRecord;
use crate::vcf::error::KeyError;
use crate::vcf::header::{SampleIdx, VcfHeader};
use crate::vcf::record_encoder::{FieldKind, FilterId, FormatKey, InfoKey};
use seqair_types::SmolStr;
use std::cell::Cell;

// r[impl vcf_fields.trait+3]
/// A type a reader reads records into. It names what is decoded: a field
/// the type does not hold is never read.
///
/// Implemented by every component ([`Locus`](super::Locus),
/// [`Site`](super::Site), [`Info`](super::Info),
/// [`OnlySample`](super::OnlySample), [`Sample`](super::Sample),
/// [`Samples`](super::Samples)), by [`FullRecord`](super::FullRecord), by
/// tuples of up to eight of them, and by any type of yours.
///
/// There is no error type of its own: a domain outcome that is not a
/// malformed record (a call without GT, a site to skip) belongs in the
/// value, as `struct Call(Result<VariantCall, CallError>)` or
/// `Truth(Option<TruthSnv>)`. A record that makes 0..n values of yours is
/// read as a tuple of components and turned into them afterwards.
pub trait Fields: Sized + Send + Sync + Clone {
    /// Everything that depends only on the header: resolved keys, the
    /// selected sample. A plan of your own with runtime context goes to
    /// [`ReaderBuilder::with_plan`](super::ReaderBuilder::with_plan).
    type Plan: Send + Sync;

    /// Resolve and check everything against the file, once per reader.
    fn plan(cx: &PlanContext<'_>) -> Result<Self::Plan, PlanError>;

    /// The parts of a record [`fill`](Self::fill) reads.
    fn needs(plan: &Self::Plan) -> Sections;

    /// The record `raw`, decoded.
    ///
    /// # Errors
    /// What decoding `raw` fails with.
    fn new(plan: &Self::Plan, raw: RawRecord<'_>) -> Result<Self, VcfReadError>;

    /// Overwrite `self` with the record `raw`. Override it to reuse
    /// `self`'s buffers; the default builds a new value.
    ///
    /// # Errors
    /// As [`new`](Self::new).
    fn fill(&mut self, plan: &Self::Plan, raw: RawRecord<'_>) -> Result<(), VcfReadError> {
        *self = Self::new(plan, raw)?;
        Ok(())
    }
}

// r[impl vcf_fields.sample+2]
/// The sample a reader's builder selects with
/// [`.sample(sel)`](super::ReaderBuilder::sample): by name, or by a
/// [`SampleIdx`] of the file's header.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SampleSel {
    Name(SmolStr),
    Idx(SampleIdx),
}

impl From<&str> for SampleSel {
    fn from(name: &str) -> Self {
        Self::Name(name.into())
    }
}

impl From<&String> for SampleSel {
    fn from(name: &String) -> Self {
        Self::Name(name.into())
    }
}

impl From<String> for SampleSel {
    fn from(name: String) -> Self {
        Self::Name(name.into())
    }
}

impl From<SmolStr> for SampleSel {
    fn from(name: SmolStr) -> Self {
        Self::Name(name)
    }
}

impl From<SampleIdx> for SampleSel {
    fn from(sample: SampleIdx) -> Self {
        Self::Idx(sample)
    }
}

/// An option a reader's builder can set.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadOption {
    /// `.sample(name)`.
    Sample,
}

impl std::fmt::Display for ReadOption {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Sample => "sample",
        })
    }
}

/// Why a [`Fields`] type cannot read a file.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum PlanError {
    #[error(transparent)]
    Key(#[from] KeyError),

    #[error("OnlySample reads files with exactly one sample; this one has {n}")]
    NotSingleSample { n: usize },

    #[error("the file has {n} samples: choose one with .sample(name)")]
    SampleNotSelected { n: usize },

    // r[impl vcf_fields.unused_options+2]
    #[error("the {option} option was set, but nothing read uses it")]
    UnusedOption { option: ReadOption },

    #[error("{name} is of another header than this file's")]
    ForeignKey { name: SmolStr },
}

// r[impl vcf_fields.unused_options+2]
/// Options a component consumes when it plans, as the reader's builder set
/// them. Reading an option marks it used; one left unread fails the open
/// with [`PlanError::UnusedOption`].
#[derive(Debug, Default)]
pub(crate) struct ReadOptions {
    sample: Option<SampleSel>,
    sample_read: Cell<bool>,
}

impl ReadOptions {
    pub(crate) fn new(sample: Option<SampleSel>) -> Self {
        Self { sample, sample_read: Cell::new(false) }
    }

    /// The sample `.sample(sel)` selected.
    fn sample(&self) -> Option<&SampleSel> {
        self.sample_read.set(true);
        self.sample.as_ref()
    }

    /// The first option set that nothing read.
    pub(crate) fn unused(&self) -> Option<ReadOption> {
        (self.sample.is_some() && !self.sample_read.get()).then_some(ReadOption::Sample)
    }
}

// r[impl vcf_fields.plan_context]
/// What [`Fields::plan`] resolves against: the file's header and field
/// table, and the options the reader's builder set.
#[derive(Debug)]
pub struct PlanContext<'a> {
    file: &'a FileShared,
    opts: &'a ReadOptions,
}

/// `Ok(None)` for a key the file does not have, every other error kept.
pub(crate) fn optional<K>(key: Result<K, KeyError>) -> Result<Option<K>, KeyError> {
    match key {
        Ok(key) => Ok(Some(key)),
        Err(KeyError::UnknownField { .. }) => Ok(None),
        Err(e) => Err(e),
    }
}

impl<'a> PlanContext<'a> {
    pub(crate) fn new(file: &'a FileShared, opts: &'a ReadOptions) -> Self {
        Self { file, opts }
    }

    pub fn header(&self) -> &'a VcfHeader {
        &self.file.header
    }

    /// An INFO key, as [`VcfFile::info_key`](super::VcfFile::info_key)
    /// resolves it.
    ///
    /// # Errors
    /// As `VcfFile::info_key`.
    pub fn info_key<V: FieldKind>(&self, name: &str) -> Result<InfoKey<V>, KeyError> {
        self.file.info_key(name)
    }

    /// A FORMAT key, as [`VcfFile::format_key`](super::VcfFile::format_key)
    /// resolves it.
    ///
    /// # Errors
    /// As `VcfFile::format_key`.
    pub fn format_key<V: FieldKind>(&self, name: &str) -> Result<FormatKey<V>, KeyError> {
        self.file.format_key(name)
    }

    /// An INFO key, or `None` if the file does not have it.
    ///
    /// # Errors
    /// As [`Self::info_key`], but for `UnknownField`.
    pub fn info_key_opt<V: FieldKind>(&self, name: &str) -> Result<Option<InfoKey<V>>, KeyError> {
        optional(self.info_key(name))
    }

    /// A FORMAT key, or `None` if the file does not have it.
    ///
    /// # Errors
    /// As [`Self::format_key`], but for `UnknownField`.
    pub fn format_key_opt<V: FieldKind>(
        &self,
        name: &str,
    ) -> Result<Option<FormatKey<V>>, KeyError> {
        optional(self.format_key(name))
    }

    /// A FILTER the header declares.
    ///
    /// # Errors
    /// [`KeyError::UnknownFilter`].
    pub fn filter(&self, name: &str) -> Result<FilterId, KeyError> {
        self.file.header.filter(name)
    }

    // r[impl vcf_fields.sample+2]
    /// The sample the builder's `.sample(sel)` selected, or the only one:
    /// what [`Sample`](super::Sample) reads.
    ///
    /// # Errors
    /// [`KeyError::UnknownSample`] for a name the header lacks,
    /// [`PlanError::ForeignKey`] for a `SampleIdx` of another header, and
    /// [`PlanError::SampleNotSelected`] when nothing was selected and there
    /// is not exactly one sample.
    pub fn selected_sample(&self) -> Result<SampleIdx, PlanError> {
        let header = &self.file.header;
        match self.opts.sample() {
            Some(SampleSel::Name(name)) => Ok(header.sample(name)?),
            Some(&SampleSel::Idx(sample)) if sample.header == header.id() => Ok(sample),
            Some(SampleSel::Idx(sample)) => {
                Err(PlanError::ForeignKey { name: format!("sample {}", sample.index()).into() })
            }
            None => {
                let mut samples = header.sample_indices();
                match (samples.len(), samples.next()) {
                    (1, Some(sample)) => Ok(sample),
                    (n, _) => Err(PlanError::SampleNotSelected { n }),
                }
            }
        }
    }
}

/// Which samples' FORMAT values a type reads.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SampleSelection {
    #[default]
    None,
    One(SampleIdx),
    All,
}

// r[impl vcf_fields.sections]
/// The parts of a record a [`Fields`] type reads, combined with `|`:
/// `Sections::SITE | Sections::INFO`. The fixed fields (contig, position,
/// reference length) are always read; whatever is not named is skipped
/// without being copied.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Sections {
    site: bool,
    info: bool,
    samples: SampleSelection,
}

impl Sections {
    /// The fixed fields only.
    pub const FIXED: Self = Self { site: false, info: false, samples: SampleSelection::None };
    /// ID, REF, ALT, QUAL and FILTER.
    pub const SITE: Self = Self { site: true, ..Self::FIXED };
    pub const INFO: Self = Self { info: true, ..Self::FIXED };
    /// Every sample's FORMAT values.
    pub const SAMPLES: Self = Self { samples: SampleSelection::All, ..Self::FIXED };

    /// One sample's FORMAT values.
    pub const fn sample(sample: SampleIdx) -> Self {
        Self { samples: SampleSelection::One(sample), ..Self::FIXED }
    }

    /// Everything either reads; two different single samples read all.
    #[must_use]
    pub fn union(self, other: Self) -> Self {
        use SampleSelection::{All, None, One};
        let samples = match (self.samples, other.samples) {
            (None, s) | (s, None) => s,
            (One(a), One(b)) if a == b => One(a),
            _ => All,
        };
        Self { site: self.site || other.site, info: self.info || other.info, samples }
    }

    pub fn site(&self) -> bool {
        self.site
    }

    pub fn info(&self) -> bool {
        self.info
    }

    pub fn samples(&self) -> SampleSelection {
        self.samples
    }

    /// Whether sample `index`'s values are read.
    pub(crate) fn has_sample(&self, index: usize) -> bool {
        match self.samples {
            SampleSelection::None => false,
            SampleSelection::One(s) => s.index() == index,
            SampleSelection::All => true,
        }
    }
}

impl std::ops::BitOr for Sections {
    type Output = Self;

    fn bitor(self, other: Self) -> Self {
        self.union(other)
    }
}

impl std::ops::BitOrAssign for Sections {
    fn bitor_assign(&mut self, other: Self) {
        *self = self.union(other);
    }
}

// r[impl vcf_fields.tuples]
macro_rules! tuple_fields {
    ($($t:ident $i:tt),+) => {
        impl<$($t: Fields),+> Fields for ($($t,)+) {
            type Plan = ($($t::Plan,)+);

            fn plan(cx: &PlanContext<'_>) -> Result<Self::Plan, PlanError> {
                Ok(($($t::plan(cx)?,)+))
            }

            fn needs(plan: &Self::Plan) -> Sections {
                Sections::FIXED $(| $t::needs(&plan.$i))+
            }

            fn new(plan: &Self::Plan, raw: RawRecord<'_>) -> Result<Self, VcfReadError> {
                Ok(($($t::new(&plan.$i, raw)?,)+))
            }

            fn fill(&mut self, plan: &Self::Plan, raw: RawRecord<'_>) -> Result<(), VcfReadError> {
                $(self.$i.fill(&plan.$i, raw)?;)+
                Ok(())
            }
        }
    };
}

tuple_fields!(A 0);
tuple_fields!(A 0, B 1);
tuple_fields!(A 0, B 1, C 2);
tuple_fields!(A 0, B 1, C 2, D 3);
tuple_fields!(A 0, B 1, C 2, D 3, E 4);
tuple_fields!(A 0, B 1, C 2, D 3, E 4, F 5);
tuple_fields!(A 0, B 1, C 2, D 3, E 4, F 5, G 6);
tuple_fields!(A 0, B 1, C 2, D 3, E 4, F 5, G 6, H 7);
