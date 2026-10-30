//! A banded affine-gap pair-HMM with a pluggable emission model.
//!
//! The dynamic program is GATK's: match, insertion and deletion states, with
//! the transitions out of a match derived from per-base insertion and deletion
//! qualities and the transitions out of an indel from a gap-continuation
//! penalty. It departs from GATK's only at the two ends of a read, which it
//! scores alike ([below](#where-compair-differs-from-gatk)). What is pluggable is the match state's emission, and the reason is
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
//! - [`Workspace::align_batch`], which scores a whole batch of candidate haplotypes
//!   against one read in lockstep. It computes all [`BATCH`] lanes whatever the
//!   caller asked for, so it wins on a full batch (~1.3x over the strip kernel)
//!   and loses badly on a short one (~0.2x at a single haplotype).
//!   [`BATCH_BREAK_EVEN`] is where the two meet.
//! - [`Workspace::align_pairs`], which scores eight *unrelated* pairs in lockstep: each
//!   lane its own read, haplotype, strand and band offset, the half-width
//!   shared. It is the batch kernel's traversal with the read and the band
//!   made per lane, so a group fills from any list of pairs -- in particular
//!   from many reads against two or three haplotypes, where the batch kernel
//!   would run three lanes of eight. [`Workspace::align_reads`] packs a reads-by-haplotypes
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
//! band, to a millionth, down to a floor near `-600` (for 150 bases; see
//! [`align_banded_f64`], which keeps its rows at `2^960`). The kernels compute in `f32`, renormalising by powers
//! of two, and flush stored cells below `2^-126`. The strip, batch, pairs and
//! GPU kernels keep a renormalised row's largest cell at `2^115`, and the
//! diagonal kernel an anti-diagonal's, so a flush provably cannot move a
//! total above a floor near `-62` (for 150 bases in a 64-wide band) by a
//! millionth. A pair that finishes below the floor is scored again by
//! [`align_banded_f64`], and so is a read whose qualities could overflow the
//! kernels. The GPU's scores carry the same check with them (`gpu::GpuScore`).
//! The proof
//! needs every emission and the transitions out of every state to be at most
//! one: a base error
//! probability is capped at [`MAX_EPSILON`], and gap-open qualities below Q6
//! count as Q6, as GATK raises them. It does not need a cell to be a
//! probability, which it is not where a read's gap qualities change from
//! row to row: GATK's recurrence takes a cell's transitions from two rows.
//!
//! # Where compair differs from GATK
//!
//! At the two ends of a read, and nowhere else.
//!
//! **What GATK does.** Its pair-HMM starts every read in the deletion row
//! (`D[0][j] = 1 / |H|`, the match and insertion rows zero) and sums the last
//! row's match *and* insertion cells. So a read's first base must be matched,
//! but its last base may be an insertion. GATK's own write-up of the model
//! (David Benjamin, *Pair HMM probabilistic realignment in HaplotypeCaller and
//! Mutect*, `docs/pair_hmm.pdf` in the GATK repository, footnote 4) notes that
//! an alignment cannot start with an insertion though it may end with one,
//! and that "this limitation is unnecessary and could be fixed".
//!
//! **Why it matters.** A base at either end of a read that disagrees with the
//! haplotype is a mismatch, `eps / 3`; at the last base it can also be a
//! one-base insertion, `match_to_insertion`. At Q40 bases and a Q44 gap open
//! that is `3.3e-5` against `4.0e-5`, so a mismatch at a read's last base is
//! about 2.2 times as likely as the same mismatch at its first base: `0.34`
//! in log10 per read, against the haplotype that lacks the base. For a base
//! the conversion-aware emission scores as converted the gap is smaller
//! (`0.015` at a conversion rate of `1.1e-3`). Nothing in the chemistry
//! favours one end in reference orientation: the pair-HMM sees a forward
//! read's 3' end on the right and a reverse read's on the left. rastair found
//! this because a scenario and its reverse complement, every read on the
//! other strand, were genotyped with different likelihoods: in the mirror the
//! reads that had ended on a variant began on it.
//!
//! **What compair does instead.**
//!
//! - *A read may begin with an insertion*, opened at its first base's
//!   insertion quality, as it may end with one. The free start goes to the
//!   first match at `indel_to_match` per column, as in GATK, and to an
//!   insertion of the first base at `match_to_insertion`. This is the begin
//!   state of the pair-HMM of Durbin et al. (*Biological Sequence Analysis*,
//!   ch. 4), which may enter any state, made to mirror GATK's end.
//! - *An insertion lies between two haplotype bases*, at the end as at the
//!   start: the insertion cell of the last haplotype column is left out of
//!   the total, as no read can insert before the first column. So a read must
//!   lie within its haplotype, as GATK assumes of its haplotypes; a read that
//!   runs past the haplotype's end, which GATK scores as an insertion of the
//!   overhang, is impossible here, as a read that starts before it already is
//!   in both.
//!
//! With these the model is mirror-symmetric: scoring a read and a haplotype
//! gives what scoring their reverse complements gives, wherever the read's gap
//! qualities are the same at every base. Where they vary along the read (as a
//! tandem-repeat gap model makes them) it is not and cannot be, since GATK's
//! recurrence charges a gap's opening and closing to one of the two rows it
//! joins.
//!
//! **When the numbers differ from GATK's.** Every read picks up the paths that
//! begin with an insertion, about `match_to_insertion` of its likelihood
//! (`1.7e-5` in log10 at Q44) when its first base fits the haplotype. GATK's
//! 104 published vectors (Q40 insertion qualities) move by `4e-5` to
//! `1.1e-4`, and agree with GATK to `1e-5` once the paths that begin with an
//! insertion are added to its numbers, which `tests/gatk_vectors.rs` does
//! independently of the kernels. The difference is large only
//! where a read's first base disagrees with the haplotype, which now scores
//! as a mismatch at the last base always did, and where a read runs to the
//! haplotype's last base and its last base disagrees, which now scores as a
//! mismatch at the first base always did.
//!
//! **Why GATK keeps it**, as far as we can tell: the effect is confined to
//! the last one or two bases of a read (a two-base insertion costs another
//! `gap_continuation`, `0.1`), GATK caps each read's likelihood with a global
//! mismapping rate and clips low-quality read tails, and its recurrence is
//! reproduced bit for bit by the Java code, Intel's vectorised native kernels,
//! FPGA and GPU ports, all of which a change would have to follow.

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
mod rescue;
mod scaling;
mod simd;
mod strips;
mod transitions;

pub use banded::{Band, Workspace, align_banded, align_banded_simd};
pub use batch::{BATCH, BATCH_BREAK_EVEN};
pub use emission::{
    Betas, ConversionModel, Emission, MAX_EPSILON, MatchProbability, SiteWeights, StandardEmission,
    TapsEmission,
};
pub use error::Error;
pub use haplotype::Haplotype;
pub use pairs::{PAIRS, PAIRS_BREAK_EVEN, Pair, ScoreMatrix};
pub use prepared::Candidates;
pub use read::{PcrIndelModel, Read};
pub use reference::{align_banded_f64_rows, align_full};
pub use rescue::{align_banded_f64, align_banded_f64_at, poison_align_banded_f64_scratch};
pub use seqair_types::{Base, BaseQuality, Probability, QPos, Strand};
pub use simd::simd_level;
pub use strips::{align_strips, align_strips_simd};
pub use types::{CpgRole, HapPos, HapSite, Log10Likelihood, Observation, error_probability};

mod types;
