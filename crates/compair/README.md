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
use compair::{Band, Base, BaseQuality, Betas, ConversionModel, Haplotype, Probability, Read, Strand, TapsEmission, align_banded_simd};

let haplotype = Haplotype::from_ascii(b"ACGTTAGCATCGGATCC...");
let read = Read::uniform(
    Base::from_ascii_vec(b"TTAGCATCGG...".to_vec()),
    &quals,                        // one BaseQuality per base
    BaseQuality::from_byte(45),    // gap open, insertion side
    BaseQuality::from_byte(45),    // gap open, deletion side
    BaseQuality::from_byte(10),    // gap continuation
    Strand::OT,
)?;
let emission = TapsEmission::new(ConversionModel::taps_default(), Betas::Uniform(Probability::new(0.5)?));
let log10 = align_banded_simd(&haplotype, &read, &emission, Band::anchored(3));
```

`align_full` is the `f64` reference over the whole matrix; `align_banded` and
`align_banded_simd` are the same `f32` recurrence over a diagonal band, one
generic function over a lane type, so the two are bit-identical by
construction. `StandardEmission` is GATK's plain emission; `TapsEmission`
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
- 8-wide `f32` SIMD via `wide`, ~17 µs per 150 bp read × 48-wide band on an
  Apple M4 Pro, bit-identical to the scalar kernel
- Types shared with `seqair-types`: `Base`, `Strand`, `BaseQuality`,
  `Probability`, `QPos`

## Testing

```sh
cargo test -p compair
cargo bench -p compair --bench align
```

`tests/gatk_vectors.rs` runs GATK's 104 published test vectors through the
reference (agreement to 1e-4 in log10, measured 6e-6) and through the band;
`tests/conventions.rs` pins the conventions those vectors cannot, with a
hand-computed 3 × 2 matrix. Property tests compare the banded kernels against
the reference and an independently masked `f64` DP, and check the scalar and
SIMD kernels bit for bit.

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
- The band as a contract rather than a heuristic

## License

MIT/Apache-2.0
