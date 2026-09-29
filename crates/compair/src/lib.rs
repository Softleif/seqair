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
//! Scoring two reads against a reference and an alternate haplotype:
//!
//! ```
//! use compair::{Band, Base, BaseQuality, Haplotype, Read, StandardEmission, Strand, Workspace};
//!
//! let reference = Haplotype::from_ascii(b"GACAATTACATAACATACACATCAGCACAAAACTTGTTGG");
//! let alternate = Haplotype::from_ascii(b"GACAATTACATAACATACACATCATCACAAAACTTGTTGG");
//! let read = |seq: &[u8]| {
//!     Read::uniform(
//!         Base::from_ascii_vec(seq.to_vec()),
//!         &vec![BaseQuality::from_byte(30); seq.len()], // base qualities
//!         BaseQuality::from_byte(45),                   // insertion gap open
//!         BaseQuality::from_byte(45),                   // deletion gap open
//!         BaseQuality::from_byte(10),                   // gap continuation
//!         Strand::OT,
//!     )
//! };
//! let (from_ref, from_alt) = (read(b"TACATAACATACACATCAGCACAA")?, read(b"TACATAACATACACATCATCACAA")?);
//! // Each read with a band centred on where its first base aligns.
//! let reads = [(&from_ref, Band::anchored(6)), (&from_alt, Band::anchored(6))];
//!
//! let mut workspace = Workspace::new(); // reuse it: scoring then allocates nothing
//! let mut scores = Vec::new();
//! workspace.align_reads(&[&reference, &alternate], &reads, &StandardEmission::default(), &mut scores);
//!
//! // Read-major log10 likelihoods: [ref read vs ref, vs alt, alt read vs ref, vs alt].
//! let log10: Vec<f64> = scores.iter().map(|score| score.get()).collect();
//! assert!(matches!(log10[..], [a, b, c, d] if a > b && d > c));
//! # Ok::<(), compair::Error>(())
//! ```
//!
//! `examples/` has runnable versions: `quickstart` (per-read evidence and
//! genotype likelihoods), `taps` (the conversion-aware emission against the
//! plain one), `candidates` (many haplotypes, one read at a time) and `gpu`.
//!
//! **Which entry point:**
//!
//! - many reads against a *few* haplotypes -- a variant caller's shadow
//!   scoring, every read of a locus against the reference and each candidate
//!   allele -- is [`Workspace::align_reads`], which packs the reads-by-haplotypes
//!   product eight pairs to a vector;
//! - one read, or reads arriving one at a time, against *many* candidate
//!   haplotypes is [`Workspace::candidates`] once and then [`Candidates::align`]
//!   per read, which prepares each haplotype's column tracks once per strand
//!   and runs full groups through the batch kernel; [`Workspace::align_candidates`]
//!   is the same dispatch without the preparation, for a single read.
//!
//! Each picks between the kernels below for you, and none of them changes a
//! score: every path is bit-identical to the strip kernel. Measured on rastair's
//! shape (3 haplotypes, 128 reads) `align_reads` is 1.5x `Candidates` on a
//! 3950X and 1.2x on an M4 Pro; on the 10s dataset, whose groups run to 24
//! haplotypes, `Candidates` is 1.13x `align_reads` on the 3950X and 0.96x on
//! the M4. The rest of this list is what they pick between, and
//! what to reach for when the shape of the work is different.
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
//! - [`align_pairs`], which scores eight *unrelated* pairs in lockstep: each
//!   lane its own read, haplotype, strand and band offset, the half-width
//!   shared. It is the batch kernel's traversal with the read and the band
//!   made per lane, so a group fills from any list of pairs -- in particular
//!   from many reads against two or three haplotypes, where the batch kernel
//!   would run three lanes of eight. [`align_reads`] packs a reads-by-haplotypes
//!   product into it and sends a short last group through the strip kernel.
//!
//! Each banded pair is one generic function over a lane type, so its scalar
//! and SIMD kernels are bit-identical by construction rather than by
//! agreement, and the batch kernel is bit-identical to the strip kernel for the
//! same reason, as is every lane of the pairs kernel -- which is what lets
//! [`Workspace::align_candidates`] and [`Workspace::align_reads`] switch
//! between kernels without changing an answer. A [`Workspace`] keeps every
//! kernel's buffers between calls, which makes an alignment allocation-free.
//!
//! **Precision.** Every entry point returns the `f64` recurrence over the
//! band, at any score. The kernels compute in `f32`, renormalising by powers
//! of two, and flush stored cells below `2^-126`. The strip, batch, pairs and
//! GPU kernels keep a renormalised row's largest cell at `2^96`, so a flush
//! there provably cannot move a total above a floor near `-56` (for 150 bases
//! in the default band) by a millionth; the diagonal kernel keeps it in
//! `[1, 2)`, and its floor is near `-27`. A pair that finishes below its
//! kernel's floor is scored again by [`align_banded_f64`]. [`trusted`] is
//! that check, for the GPU kernel's scores. Two model rules keep the proof's
//! premise -- every cell a probability -- true for every quality a [`Read`]
//! accepts: a base error probability is capped at [`MAX_EPSILON`], and
//! gap-open qualities below Q6 count as Q6, as GATK raises them.

mod banded;
mod batch;
mod emission;
mod error;
#[cfg(feature = "gpu")]
pub mod gpu;
mod haplotype;
mod lanes;
mod pairs;
#[cfg(test)]
mod pinned;
mod prepared;
mod read;
mod reference;
mod scaling;
mod simd;
mod strips;
mod transitions;

pub use banded::{Band, Workspace, align_banded, align_banded_simd};
pub use batch::{BATCH, BATCH_BREAK_EVEN, align_batch, align_candidates};
pub use emission::{
    Betas, ConversionModel, Emission, MAX_EPSILON, MatchProbability, SiteWeights, StandardEmission,
    TapsEmission,
};
pub use error::Error;
pub use haplotype::Haplotype;
pub use pairs::{PAIRS, PAIRS_BREAK_EVEN, Pair, align_pairs, align_reads};
pub use prepared::Candidates;
pub use read::Read;
pub use reference::{align_banded_f64, align_full, trusted};
pub use seqair_types::{Base, BaseQuality, Probability, QPos, Strand};
pub use simd::simd_level;
pub use strips::{align_strips, align_strips_simd};
pub use types::{CpgRole, HapPos, HapSite, Log10Likelihood, Observation, error_probability};

mod types;
