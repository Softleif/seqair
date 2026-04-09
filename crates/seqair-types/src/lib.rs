//! This crate provides common types used throughout Seqair.
#![deny(missing_docs)]
#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

pub mod bam_flags;
mod base;
mod phred;
pub mod pos;
mod probability;
mod region_string;
mod rms;
mod strand;

pub use smol_str::{self, SmolStr};

pub use smallvec;
pub use smallvec::SmallVec;

pub use pos::{Offset, One, Pos, Pos0, Pos1, Zero};

pub use {
    bam_flags::BamFlags,
    base::{Base, BaseError},
    phred::Phred,
    probability::{Probability, ProbabilityError},
    region_string::{RegionString, RegionStringError},
    rms::{RmsAccumulator, RootMeanSquare, RootMeanSquareExt},
    strand::{Strand, StrandFromRecord, strand_from_flags},
};
