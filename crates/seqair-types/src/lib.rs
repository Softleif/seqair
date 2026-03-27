//! This crate provides common types used throughout Seqair.
#![deny(missing_docs)]
#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

mod base;
mod phred;
mod probability;
mod region_string;
mod rms;
mod strand;

pub use smol_str::{self, SmolStr};

#[cfg(debug_assertions)]
/// A simple replacement for smallvec in debug mode.
pub mod smallvec {
    pub use std::vec as smallvec;
    pub use std::vec as smallvec_inline;
    /// A simple replacement for smallvec in debug mode.
    pub type SmallVec<T, const _N: usize> = Vec<T>;
}
#[cfg(not(debug_assertions))]
pub use smallvec;

pub use smallvec::SmallVec;

pub use {
    base::{Base, BaseError},
    phred::Phred,
    probability::{Probability, ProbabilityError},
    region_string::{RegionString, RegionStringError},
    rms::{RmsAccumulator, RootMeanSquare, RootMeanSquareExt},
    strand::{Strand, StrandFromRecord, strand_from_flags},
};
