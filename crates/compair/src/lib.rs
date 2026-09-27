//! A banded affine-gap pair-HMM with a pluggable emission model.
//!
//! The dynamic program is GATK's: match, insertion and deletion states, with
//! the transitions out of a match derived from per-base insertion and deletion
//! qualities and the transitions out of an indel from a gap-continuation
//! penalty. What is pluggable is the match state's emission, and the reason is
//! conversion chemistry: under TAPS a methylated cytosine reads as `T`, so
//! scoring a `T` against a `C` as a plain mismatch loses a real alignment while
//! scoring it as a free match loses every real `C>T` variant. `TapsEmission`
//! scores it as a probability at the site's methylation level instead.
//!
//! **If you are scoring one read against its candidate haplotypes, which is
//! what a variant caller does, call [`Workspace::align_candidates`].** It picks
//! between the kernels below for you. The rest of this list is what it picks
//! between, and what to reach for when the shape of the work is different.
//!
//! Implementations of the same recurrence, in two families.
//!
//! One pair at a time, eight *cells* of it per vector:
//!
//! - [`align_full`] in `f64` over the whole matrix, the reference,
//! - [`align_banded`] in `f32` over a diagonal band, one anti-diagonal at a
//!   time, and [`align_banded_simd`], the same eight lanes at a time,
//! - [`align_strips`] over the same band, one read row at a time along the
//!   haplotype, and [`align_strips_simd`], eight rows at a time -- the
//!   traversal of Intel's Genomics Kernel Library, and the faster of the two.
//!
//! Eight *alignments* per vector, one haplotype per lane:
//!
//! - [`align_batch`], which scores a whole batch of candidate haplotypes
//!   against one read in lockstep. It computes all [`BATCH`] lanes whatever the
//!   caller asked for, so it wins on a full batch (~1.3x over the strip kernel)
//!   and loses badly on a short one (~0.2x at a single haplotype).
//!   [`BATCH_BREAK_EVEN`] is where the two meet.
//!
//! Each banded pair is one generic function over a lane type, so its scalar
//! and SIMD kernels are bit-identical by construction rather than by
//! agreement, and the batch kernel is bit-identical to the strip kernel for the
//! same reason -- which is what lets [`Workspace::align_candidates`] switch
//! between them without changing an answer. A [`Workspace`] keeps every
//! kernel's buffers between calls, which makes an alignment allocation-free.

mod banded;
mod batch;
mod emission;
mod error;
#[cfg(feature = "gpu")]
pub mod gpu;
mod haplotype;
mod read;
mod reference;
mod scaling;
mod simd;
mod strips;
mod transitions;

pub use banded::{Band, Workspace, align_banded, align_banded_simd};
pub use batch::{BATCH, BATCH_BREAK_EVEN, align_batch, align_candidates};
pub use emission::{
    Betas, ConversionModel, Emission, MatchProbability, SiteWeights, StandardEmission, TapsEmission,
};
pub use error::Error;
pub use haplotype::Haplotype;
pub use read::Read;
pub use reference::align_full;
pub use seqair_types::{Base, BaseQuality, Probability, QPos, Strand};
pub use simd::simd_level;
pub use strips::{align_strips, align_strips_simd};
pub use types::{CpgRole, HapPos, HapSite, Log10Likelihood, Observation, error_probability};

mod types;
