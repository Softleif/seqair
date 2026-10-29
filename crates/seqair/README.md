# Seqair

Seqair is a library featuring BAM/SAM/CRAM/FASTA readers, a pileup engine, and BCF/BAM writing tools.

[Documentation](https://docs.rs/seqair)

It was written around improved memory locality and I/O:

- It is based around a slab-based `RecordStore`
- `RegionBuf`, the built-in readers, and `ChunkCache` use bulk I/O, aiming for few but large reads instead of small recurring ones (better on network storage)
- Seqair readers can be `fork()`ed for cheap per-thread readers sharing index and headers

It supports reading SAM, BAM, as well as CRAM v3.0/v3.1,
with SIMD-opimtized custom implementations (for both ARM64 NEON and X86 with AVX2 or AVX512).
In the same style we also provide SAM/BAM writers
as well as a VCF/BCF writer,
which uses pre-resolved typed handles to write directly into output buffers.

For variant callers like Rastair,
we also added a way to attach custom data while reads are being loaded
so you can access them in the pileup engine.

You can read the introduction blog post [here](https://deterministic.space/seqair.html).

## What's new vs. adapted

The `RecordStore` design with its customizable column-based memory layout,
and the type-safe VCF/BCF writer
are new to Seqair.
The heavy use of newtypes (e.g. `Pos0`/`Pos1`, `Phred`, `BaseQuality`, etc.)
are our versions of the data types described in the hts-specs encoded in the type-system.
While this was design by us, most of the code was written with Claude Code
which enabled use to reach the current feature-rich state early.
Aside from code review, we validate correctness with comparison tests against
[noodles](https://github.com/zaeleus/noodles) and [rust-htslib](https://github.com/rust-bio/rust-htslib)/htslib,
including property tests and fuzzing.

We also based a lot of our implemenation on the fantastic work from other people and projects:

- BAM, SAM, CRAM, VCF/BCF, FASTA, BGZF, BAI/CSI/TBI, FAI/GZI: implementations following the [hts-specs](https://github.com/samtools/hts-specs)
- CRAM codecs (rANS 4x8, rANS NX16, tok3, Huffman, ITF8/LTF8): implemented from the CRAM spec with test vectors taken from noodles
- Index builder (`hts_idx_push` state machine): ported from htslib's single-pass index co-production algorithm
- CRAM substitution matrix defaults: validated against htslib
- Pileup compatibility: output conventions match htslib's pileup
- VCF number formatting: matches htslib output

## License

MIT/Apache-2.0
