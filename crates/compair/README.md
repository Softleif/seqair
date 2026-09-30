# compair

A fast pair-HMM for scoring sequencing reads against candidate haplotypes.
It computes the same likelihoods as GATK's `PairHMM`, but lets you plug in
the emission model, so that conversion chemistry (TAPS, 5-base) is scored as
a probability instead of a mismatch.
It is the read-likelihood engine of [rastair](https://github.com/bsbludwig/rastair).

> [!WARNING]
> **This is an experimental project!**
> It reproduces GATK's numbers on GATK's published test vectors,
> but it has not yet been used outside rastair.

## Usage

```rust
use compair::{Band, Base, BaseQuality, Haplotype, Read, StandardEmission, Strand, Workspace};

let reference = Haplotype::from_ascii(b"GACAATTACATAACATACACATCAGCACAAAACTTGTTGG");
let alternate = Haplotype::from_ascii(b"GACAATTACATAACATACACATCATCACAAAACTTGTTGG");
let read = Read::uniform(
    Base::from_ascii_vec(b"TACATAACATACACATCATCACAA".to_vec()),
    &vec![BaseQuality::from_byte(30); 24], // base qualities
    BaseQuality::from_byte(45),            // insertion gap open
    BaseQuality::from_byte(45),            // deletion gap open
    BaseQuality::from_byte(10),            // gap continuation
    Strand::OT,
)?;

let mut workspace = Workspace::new(); // reuse it: scoring then allocates nothing
let mut scores = Vec::new();
workspace.align_reads(
    &[&reference, &alternate],
    &[(&read, Band::anchored(6))], // band around where the read's first base aligns
    &StandardEmission::default(),
    &mut scores,
);
// scores: log10 likelihoods, read-major; here the read favours the alternate.
```

Pick the entry point by the shape of your work:

- **Many reads, few haplotypes** (a variant caller scoring every read at a
  locus against the reference and each allele): `Workspace::align_reads`.
- **One read at a time, many haplotypes**: `Workspace::candidates`, then
  `Candidates::align` per read.

Both choose the fastest kernel for you, and the choice never changes a score.
Swap `StandardEmission` for `TapsEmission` to score TAPS reads at each
site's methylation level.

Runnable examples (`cargo run -p compair --release --example <name>`):

- `quickstart`: per-read evidence and diploid genotype likelihoods
- `taps`: why the conversion-aware emission matters
- `candidates`: many candidate haplotypes, one read at a time
- `gpu` (`--features gpu`): thousands of pairs in one launch, checked
  against the CPU

The remaining examples are maintainer tools; their header comments say what
they measure.

## Features

- GATK `PairHMM` semantics, checked against GATK's 104 test vectors and
  against rust-bio as an independent implementation
- Pluggable emission model (`Emission` trait)
- Banded, with a clear contract: a path outside the band scores lower,
  never wrong
- Exact: every score is the `f64` recurrence over its band. The kernels run
  in `f32` and rescore in `f64` the rare pair where that could matter
- SIMD on x86-64 and aarch64, selected at run time, no build flags needed.
  SIMD and scalar results are bit-identical
- Optional GPU backend (`gpu` feature, via wgpu)

How it all works, including the precision argument and the kernel trade-offs,
is in the [API docs](https://docs.rs/compair).

## Development

```sh
cargo nextest run -p compair
cargo bench -p compair --bench align
```

There is also a fuzz target, `fuzz_pair_hmm`; see the workspace `CLAUDE.md`
for how to run it.

## Credits

- The recurrence and its conventions are GATK HaplotypeCaller's `PairHMM`,
  written from the published algorithm; no GATK code.
  `tests/data/pairhmm-testdata.txt` is GATK4's own test resource
  (Broad Institute, BSD-3-Clause).
- The strip kernel's traversal is that of Intel's Genomics Kernel Library
  (MIT), written from its published description.

## License

MIT/Apache-2.0
