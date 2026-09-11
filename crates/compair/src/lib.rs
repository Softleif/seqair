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
//! Five implementations of the same recurrence:
//!
//! - [`align_full`] in `f64` over the whole matrix, the reference,
//! - [`align_banded`] in `f32` over a diagonal band, one anti-diagonal at a
//!   time, and [`align_banded_simd`], the same eight lanes at a time,
//! - [`align_strips`] over the same band, one read row at a time along the
//!   haplotype, and [`align_strips_simd`], eight rows at a time -- the
//!   traversal of Intel's Genomics Kernel Library, and the faster of the two.
//!
//! Each banded pair is one generic function over a lane type, so its scalar
//! and SIMD kernels are bit-identical by construction rather than by
//! agreement. A [`Workspace`] keeps every kernel's buffers between calls,
//! which makes an alignment allocation-free.

mod banded;
mod emission;
mod error;
mod haplotype;
mod read;
mod reference;
mod scaling;
mod strips;
mod transitions;

pub use banded::{Band, Workspace, align_banded, align_banded_simd};
pub use emission::{Betas, ConversionModel, Emission, SiteWeights, StandardEmission, TapsEmission};
pub use error::Error;
pub use haplotype::Haplotype;
pub use read::Read;
pub use reference::align_full;
pub use seqair_types::{Base, BaseQuality, Probability, QPos, Strand};
pub use strips::{align_strips, align_strips_simd};
pub use types::{CpgRole, HapPos, HapSite, Log10Likelihood, Observation, error_probability};

mod types;
