# compair

Banded affine-gap pair-HMM for scoring a read against candidate haplotypes,
with the emission model as a trait so conversion chemistry (TAPS, 5-base)
is a probability rather than a mismatch.
Read-likelihood engine for [rastair](https://github.com/bsblabludwig/rastair).

> [!WARNING]
> **This is an experimental project!**
> The kernel reproduces GATK's numbers on GATK's published test vectors,
> but it has not yet been used outside rastair.

```rust
use compair::{Band, Base, BaseQuality, Betas, ConversionModel, Haplotype, Probability, Read, Strand, TapsEmission, Workspace};

let reference = Haplotype::from_ascii(b"GACAATTACATAACATACACATCAGCACAAAACTTGTTGG");
let alternate = Haplotype::from_ascii(b"GACAATTACATAACATACACATCATCACAAAACTTGTTGG");
let read = Read::uniform(
    Base::from_ascii_vec(b"TACATAACATACACATCATCACAA".to_vec()),
    &quals,                        // one BaseQuality per base
    BaseQuality::from_byte(45),    // gap open, insertion side
    BaseQuality::from_byte(45),    // gap open, deletion side
    BaseQuality::from_byte(10),    // gap continuation
    Strand::OT,
)?;
let emission = TapsEmission::new(ConversionModel::taps_default(), Betas::Uniform(Probability::new(0.5)?));
let mut workspace = Workspace::new();
let mut scores = Vec::new();
// Every read (each with a band around where it aligns) against every
// haplotype; log10 likelihoods, read-major.
workspace.align_reads(&[&reference, &alternate], &[(&read, Band::anchored(6))], &emission, &mut scores);
```

Runnable examples, `cargo run -p compair --release --example <name>`:

- `quickstart`: reads against a reference and an alternate haplotype, turned
  into per-read evidence and diploid genotype likelihoods
- `taps`: why the conversion-aware emission matters, on reads from a
  methylated `CpG`
- `candidates`: many candidate haplotypes (repeat lengths), one read at a time
- `gpu` (`--features gpu`): one launch for thousands of pairs, checked
  against the CPU

The other examples are measurement tools for maintainers: `batchfill` and
`pairsfill` regenerate the tables behind `BATCH_BREAK_EVEN` and
`PAIRS_BREAK_EVEN`, and the rest are profiling loops and sweeps whose header
comments say what they measure.

`align_full` is the `f64` reference over the whole matrix; `align_banded` and
`align_banded_simd` are the same `f32` recurrence over a diagonal band, one
generic function over a lane type, so the two are bit-identical by
construction. `align_strips` and `align_strips_simd` are the same band
through a different traversal, the one Intel's Genomics Kernel Library uses:
one read row per lane, swept along the haplotype, so the per-row inputs stay
in registers and the neighbours are the previous two steps' vectors shifted
by a lane. It is the faster of the two on every target measured so far and
its precision does not depend on the band width, since it renormalises every
eight rows rather than per anti-diagonal.

Two kernels put eight *alignments* in a vector instead of eight cells of one:
`align_batch` scores eight haplotypes against one read, and `align_pairs`
eight unrelated read-haplotype pairs, each with its own band offset. Both are
bit-identical to the strip kernel, lane by lane, so the two entry points that
choose between the kernels cannot change a score by choosing:
**`Workspace::align_reads` for many reads against a few haplotypes** (a
variant caller's shadow scoring: every read of a locus against the reference
and each allele) and `Workspace::align_candidates` for one read against many
haplotypes.

A `Workspace` keeps the kernels' buffers between calls, so
scoring a read against its candidate haplotypes allocates nothing; the free
functions of the same names build a fresh one per call. A `Read` computes
every `10^(-Q/10)` it needs at construction, so an alignment does no
transcendental arithmetic at all. `Workspace::candidates` prepares a set of
haplotypes under one emission for scoring read after read against, deriving
each haplotype's per-column terms once per strand rather than once per pair. `StandardEmission` is GATK's plain emission; `TapsEmission`
scores an OT `T` over a haplotype `C` (and an OB `A` over a `G`) at the site's
methylation level, reading CpG context off the haplotype's own sequence so a
de-novo CpG needs no special case.

## Highlights

- GATK `PairHMM` semantics: match, insertion and deletion states, transitions
  from per-base insertion, deletion and gap-continuation qualities, free start
  and end on the haplotype, `N` matches everything
- The emission is a trait: `StandardEmission`, `TapsEmission` with a
  `ConversionModel` and per-site or uniform `Betas`, and an artifact floor
- A stated band contract: a path outside the band scores lower, never wrong;
  no path at all is `Log10Likelihood::IMPOSSIBLE`, and `Band::MAX_WIDTH` bounds
  the allocation
- Every score is the `f64` recurrence over its band, at any score: the strip
  kernels keep a renormalised row at `2^96`, so their flushes provably cannot
  move a total above `log10(9 * r * columns * rho) - 60.8` by a millionth,
  `rho` the share of the haplotype's columns the band starts in (`- 31.9` for
  the diagonal kernel, which scales to one), and a pair that finishes below
  that is rescored in `f64` (`align_banded_f64`). On rastair's chr12 pairs
  that is 2.7% of pairs. A read whose gap qualities change along it can grow
  a cell past one, as GATK's recurrence can; the floor and an overflow check
  take that from the read's qualities
- Qualities mean what a base call can mean: a base error probability is
  capped at 3/4 (Q0 and Q1 carry no information rather than "certainly
  wrong"), and gap-open qualities below Q6 count as Q6, as GATK raises them
- 8-wide `f32` SIMD via `fearless_simd`, picked at run time (NEON; SSE2
  through AVX-512), with no build flags and every SIMD kernel bit-identical
  to its scalar instance. Two lane layouts: eight cells of one alignment per
  vector (strip and diagonal kernels), or eight alignments per vector (batch
  and pairs kernels), which the entry points choose between. The default width is 46
  and not 48 because a diagonal of `width / 2 + 1` cells fills three vectors
  of eight exactly at 46 and spills one cell into a fourth at 48
- A GPU kernel behind the `gpu` feature (wgpu: Metal, Vulkan, DX12), one
  thread per (read, haplotype) pair, for callers with thousands of pairs per
  launch. It has its own oracle rather than bit-parity by construction:
  measured bit-identical to the strip kernel on an Apple M4 Pro and an AMD RX
  5700 XT except for pairs scoring below log10 -40, which a GPU that flushes
  subnormal intermediates can score lower.
- Types shared with `seqair-types`: `Base`, `Strand`, `BaseQuality`,
  `Probability`, `QPos`

## Testing

```sh
cargo test -p compair
cargo bench -p compair --bench align
```

`tests/gatk_vectors.rs` runs GATK's 104 published test vectors through the
reference (agreement to 1e-4 in log10, measured 6e-6) and through both
banded traversals; `tests/conventions.rs` pins the conventions those vectors
cannot, with a hand-computed 3 × 2 matrix. Property tests compare the banded
kernels against the reference and an independently masked `f64` DP at every
width, with `N` on both sides, check each traversal's scalar and SIMD
kernels bit for bit, fresh and through a reused `Workspace`, and check the
two traversals against each other. `tests/read.rs` pins what the constructors
reject and that lowercase is a base, not an `N`. The power-of-two scaling the
bit-parity rests on has its own property tests in `src/scaling.rs`.

`tests/rust_bio.rs` scores the same pairs through
[rust-bio](https://github.com/rust-bio/rust-bio)'s `PairHMM`, an
implementation that shares no code with this one. The two conventions map
onto each other exactly once rust-bio's `x` is the haplotype and its free
start carries the `1 / h` prior; the one term that does not map, rust-bio
summing the deletion matrix at the free end, is bounded in closed form and
asserted as a rigorous upper bound. Below that bound the agreement is
empirical, because rust-bio's three-term log-sum shortcut mis-sorts its
arguments and drops gap mass on reads that need an indel (the test documents
the measurement). `benches/rust_bio.rs` times it on the same fixture:
rust-bio takes ~900 µs where `align_full` takes ~85 µs and the SIMD band
~10 µs.

A fuzz target, `fuzz_pair_hmm` in `crates/seqair/fuzz`, drives arbitrary
input through all three kernels with the crate's documented invariants as the
oracle; see the workspace `CLAUDE.md` for how to run it.

## What's adapted from existing projects

- The recurrence and its conventions are GATK HaplotypeCaller's `PairHMM`
  (`LOGLESS_CACHING`): the `1 / hap_len` free start, the clamped
  match-to-match transition, and the final sum without the deletion matrix.
  Written from the published algorithm and checked against its numbers; no
  GATK code.
- `tests/data/pairhmm-testdata.txt` is GATK4's own test resource
  (Broad Institute, BSD-3-Clause).

## What's new

- The emission trait, and the conversion-aware emission with per-site
  methylation levels
- The banded `f32` kernel with per-anti-diagonal renormalisation, and the
  lane abstraction that makes the SIMD and scalar kernels one function
- The strip kernel's traversal is GKL's (Intel, MIT licence), written from
  the published description; its band handling, per-strip renormalisation
  and lane abstraction are this crate's
- The band as a contract rather than a heuristic

## License

MIT/Apache-2.0
