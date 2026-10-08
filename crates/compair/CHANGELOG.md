# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Breaking

- **GPU scores carry their trust floor; `compair::trusted` is private.** `Handle::collect` (and
  the emulator) return `GpuScore`s instead of `Log10Likelihood`s: `trusted()` is the score where
  a CPU entry point would keep it, `or_rescore(haplotype, read, emission, band)` scores the pair
  again in `f64` where it would not (~1% of pairs), and `untrusted()` is the kernel's own number,
  for kernel comparisons. The floor is worked out when the pair is pushed, so a caller no longer
  has to remember to pass every score through `trusted` — forgetting it gave silently wrong scores
  below log10 ~ -40.
- **The free `align_reads`, `align_pairs`, `align_batch` and `align_candidates` are gone.** Each
  built a fresh `Workspace` per call and duplicated a `Workspace` method one-for-one; keep a
  `Workspace` (`Workspace::new().align_reads(..)` for a one-off). The per-pair free functions
  stay.

### Added

- **Per-position methylation levels on a `Haplotype`** (`with_betas`, `betas`), kept through
  `with_insertion`, `with_deletion` and `reverse_complement`, and read by the new
  `Betas::OfHaplotype(fallback)`. Haplotypes with different indels applied no longer line up
  with one `Betas::PerSite` slice; this is what lets a TAPS caller score them in one call with
  the column's own methylation estimate, so a deletion spelled across a methylated `CpG` is not
  mistaken for a cheaper conversion. `HapSite` gains `beta` and loses `Eq`.
- **`PcrIndelModel`**, GATK's PCR indel gap model: the gap-open quality falls inside tandem
  repeats (units of 1 to `max_unit` bases) as `gap_open - exp(repeats / (rate * ln 2)) + 1`,
  rounded and kept between `floor` and `gap_open`. `Default` is GATK's conservative setting (Q45,
  floor Q10, rate 3, units up to 8, 20 repetitions). `gap_open_qual(repeats)` and
  `fill_gap_open(bases, track)` compute it, and **`Read::with_pcr_indel_model(bases, base_quals,
  &model, gap_qual, strand)`** builds a read from it, so a caller passes no per-base gap tracks.
- **`Workspace::align_reads_into` and `ScoreMatrix`**: `align_reads` into a buffer that knows its
  haplotype count, with `row(read)`, `rows()`, `reads()`, `haplotypes()` and `as_flat()`, so a
  caller cannot slice the read-major scores with the wrong stride. `Workspace::align_reads` stays
  for callers that want the flat buffer.
- **`Haplotype::with_insertion(anchor, &inserted)` and `Haplotype::with_deletion(anchor, len)`**:
  the haplotype with an allele applied after the base at `anchor`, `None` when it does not fit,
  with `CpG` roles recomputed from the new sequence.
- **`Log10Likelihood::as_f32`**, the documented narrowing for callers that store scores as `f32`.
- **`GpuPairs::push_reads(haplotypes, reads, emission)`**: every read against every haplotype of
  its strand in one call, read-major like `Workspace::align_reads`, returning where the scores
  land in the launch. Replaces the push loop every GPU caller wrote by hand.

### Changed

- **A `Read` keeps its four quality tracks in one allocation**, and `Read::uniform` and
  `Read::with_pcr_indel_model` write the indel tracks straight into it: building a read allocates
  four times (bases, qualities, error probabilities, transitions) instead of seven, plus three
  temporary tracks for `uniform`. Scores are unchanged. When a read has more than one problem,
  `Read::new` now reports a wrong-length track before a missing quality in an earlier one.
- **`Workspace::align_reads` and `align_reads_into` take owned or borrowed haplotypes and
  reads** (`&[H]` and `&[(R, Band)]` with `H: Borrow<Haplotype>`, `R: Borrow<Read>`), and
  `align_candidates` and `align_batch` take owned or borrowed haplotypes, so a caller's
  `Vec<(Read, Band)>` goes in as it is. Existing `&[&T]` call sites compile unchanged.
- The bit-parity oracles (`Workspace::align_batch_scalar`, `align_pairs_scalar`,
  `Candidates::align_strips`, `align_strips_simd`) are hidden from the docs.
