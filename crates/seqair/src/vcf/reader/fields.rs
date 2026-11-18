//! What a reader decodes: the [`Fields`] trait, the [`Sections`] a type
//! needs, and the options it plans with (`docs/spec/2-vcf-2-fields.md`).

use super::error::VcfReadError;
use super::raw::RawRecord;
use crate::vcf::error::KeyError;
use crate::vcf::header::{SampleIdx, VcfHeader};
use seqair_types::SmolStr;
use std::cell::Cell;

// r[impl vcf_fields.trait+2]
/// A type a reader reads records into. It names what is decoded: a field
/// the type does not hold is never read.
///
/// Implemented by every component ([`Locus`](super::Locus),
/// [`Site`](super::Site), [`Info`](super::Info),
/// [`OnlySample`](super::OnlySample), [`Sample`](super::Sample),
/// [`Samples`](super::Samples)), by [`FullRecord`](super::FullRecord), by
/// tuples of up to eight of them, and by any type of yours.
pub trait Fields: Sized + Send + Sync + Clone {
    /// Everything that depends only on the header: resolved keys, the
    /// selected sample.
    type Plan: Send + Sync;

    /// Resolve and check everything against the header, once per reader.
    fn plan(header: &VcfHeader, opts: &ReadOptions) -> Result<Self::Plan, PlanError>;

    /// The parts of a record [`fill`](Self::fill) reads.
    fn needs(plan: &Self::Plan) -> Sections;

    /// A value to fill.
    fn empty(plan: &Self::Plan) -> Self;

    /// Replace `self` with the record `raw`, reusing `self`'s buffers.
    fn fill(&mut self, plan: &Self::Plan, raw: RawRecord<'_>) -> Result<(), VcfReadError>;
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

    #[error("sample {name} is not in the header")]
    UnknownSample { name: SmolStr },

    // r[impl vcf_fields.unused_options]
    #[error("the {option} option was set, but nothing read uses it")]
    UnusedOption { option: ReadOption },
}

// r[impl vcf_fields.unused_options]
/// Options a component consumes when it plans, as the reader's builder set
/// them. Reading an option marks it used; one left unread fails the open
/// with [`PlanError::UnusedOption`].
#[derive(Debug, Default)]
pub struct ReadOptions {
    sample: Option<SmolStr>,
    sample_read: Cell<bool>,
}

impl ReadOptions {
    pub(crate) fn new(sample: Option<SmolStr>) -> Self {
        Self { sample, sample_read: Cell::new(false) }
    }

    /// The sample `.sample(name)` named.
    pub fn sample(&self) -> Option<&str> {
        self.sample_read.set(true);
        self.sample.as_deref()
    }

    /// The first option set that nothing read.
    pub(crate) fn unused(&self) -> Option<ReadOption> {
        (self.sample.is_some() && !self.sample_read.get()).then_some(ReadOption::Sample)
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

            fn plan(header: &VcfHeader, opts: &ReadOptions) -> Result<Self::Plan, PlanError> {
                Ok(($($t::plan(header, opts)?,)+))
            }

            fn needs(plan: &Self::Plan) -> Sections {
                Sections::FIXED $(| $t::needs(&plan.$i))+
            }

            fn empty(plan: &Self::Plan) -> Self {
                ($($t::empty(&plan.$i),)+)
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
