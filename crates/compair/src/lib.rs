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
//! Three implementations of the same recurrence:
//!
//! - [`align_full`] in `f64` over the whole matrix, the reference,
//! - [`align_banded`] in `f32` over a diagonal band,
//! - [`align_banded_simd`], the same band eight lanes at a time.
//!
//! The two banded kernels are one generic function over a lane type, so they
//! are bit-identical by construction rather than by agreement.

mod banded;
mod emission;
mod error;
mod haplotype;
mod read;
mod reference;
mod scaling;
mod transitions;

pub use banded::{Band, align_banded, align_banded_simd};
pub use emission::{Betas, ConversionModel, Emission, SiteWeights, StandardEmission, TapsEmission};
pub use error::Error;
pub use haplotype::Haplotype;
pub use read::Read;
pub use reference::align_full;
pub use seqair_types::{Base, BaseQuality, Probability, QPos, Strand};
pub use types::{CpgRole, HapPos, HapSite, Log10Likelihood, Observation, error_probability};

mod types;
