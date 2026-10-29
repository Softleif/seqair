# seqair

A bioinformatics toolkit written in Rust:

- [`seqair`](crates/seqair) is a BAM/SAM/CRAM/FASTA reader with a pileup engine and BCF/BAM writing.
  It is the I/O backend for [rastair](https://github.com/bsbludwig/rastair).
- [`compair`](crates/compair) is a banded pair-HMM calculating read likelihoods with a conversion-aware emission
- [`seqair-types`](crates/seqair-types) contains core types shared by both, e.g. `Base`, `Strand`, `Phred`, `Probability`, `RegionString`

> [!WARNING]
> **This is still an experimental library!**
> Only a few brave projects have started using it,
> and breaking changes still happen.

## Why?

We wrote Seqair and its sibling crates to improve performance and type-safety in Rastair.
The idea was to experiment with a different structure compared to htslib
aiming for improved memory locality and I/O.
We also wanted to provide idiomtic Rust APIs
that make it easy to do the performant and correct thing.
Both are expressed in the

Did we succeed? Fair question.
In October 2026, when switching Rastair over from rust-htslib,
its performance improved by about 4x on CPU-bound work.
Along the way, we added lots of little tweaks and features for our use-cases
that would not have been as easy in a different design,
which makes us happy but is probably not as reproducible.

Please note that most of the code was written with Claude Code.
We want to be clear and upfront about this, because this is a big library
that would not exist otherwise (or be even as close to feature complete).
The library architecture and most API design is very deliberate and human-made.
Aside from code review, we validate correctness with comparison tests against
[noodles](https://github.com/zaeleus/noodles) and [rust-htslib](https://github.com/rust-bio/rust-htslib),
including property tests and fuzzing.

You can read more about Seqair's design in [this blog post](https://deterministic.space/seqair.html).

## License

MIT/Apache-2.0
