# seqair

Pure-Rust BAM/SAM/CRAM/FASTA reader with a pileup engine and BCF/BAM writing.
I/O backend for [rastair](https://github.com/bsblabludwig/rastair).

> [!WARNING]
> **This is an experimental project!**
> While it seems that it works as expected,
> it has not undergone the extensive real-world testing needed to be confident.

This repo has two crates:

- [`seqair`](crates/seqair): Indexed readers, pileup engine, VCF/BCF/BAM writer
- [`seqair-types`](crates/seqair-types): Core types: `Base`, `Strand`, `Phred`, `Probability`, `RegionString`

## Highlights

- Spec-driven development using [tracey](https://tracey.bearcove.eu/)
- Slab-based `RecordStore` — zero per-record heap allocation
- `RegionBuf` bulk I/O — one large read per region, then decompress from memory (good for NFS/Lustre)
- CRAM v3.0/v3.1 (rANS, tok3, bzip2, lzma, multi-ref slices, embedded refs, MD5 verification)
- Overlapping read-pair dedup built into the pileup engine
- `fork()` for cheap per-thread readers sharing index and header via `Arc`
- `ChunkCache` for distant BAI bins — loads wide-spanning chunks once per chromosome

## Testing

```sh
cargo test
```

Comparison tests in `crates/seqair/tests/` validate against `rust-htslib` (and sometimes `noodles`).
Property tests use `proptest`.

## Why?

We wanted a Rust-native BAM/CRAM/FASTA/VCF stack tailored to [rastair](https://www.rastair.com/)'s specific needs
(pileup-centric, cluster-friendly I/O, zero per-record allocation)
rather than wrapping htslib's C API.
The [SAM/BAM/CRAM](https://samtools.github.io/hts-specs/) and [VCF/BCF](https://samtools.github.io/hts-specs/) specifications were our primary references throughout.

Most of the code was written with [Claude Code](https://claude.ai/code).
Aside from code review, we validate correctness
with comparison tests against [noodles](https://github.com/zaeleus/noodles) and [rust-htslib](https://github.com/rust-bio/rust-htslib)/htslib,
plus property tests and fuzzing.

### What's adapted from existing projects

- BAM, SAM, CRAM, VCF/BCF, FASTA, BGZF, BAI/CSI/TBI, FAI/GZI\*: implementations of the samtools/hts-specs
- CRAM codecs (rANS 4x8, rANS NX16, tok3, Huffman, ITF8/LTF8): implemented from the CRAM spec with test vectors taken from noodles
- Index builder (`hts_idx_push` state machine): ported from htslib's single-pass index co-production algorithm
- CRAM substitution matrix defaults: validated against htslib
- Pileup compatibility: output conventions match htslib's pileup
- VCF number formatting: matches htslib output

### What's new in seqair

- RecordStore: slab-based BAM storage packing variable-length into contiguous buffers; zero per-record heap allocation
- RegionBuf: bulk-reads all compressed bytes for a genomic region in one I/O call, then decompresses from memory; designed for high-latency filesystems (NFS/Lustre)
- ChunkCache: loads wide-spanning BAI bins (L0–L2) once per chromosome per thread instead of re-reading them per query
- `Pos<S>`: compile-time zero-based vs one-based positions instead of manual `+1`/`-1`
- Type-safe `Alleles`: `Reference`/`Snv`/`Insertion`/`Deletion`/`Complex` enum enforcing VCF structural invariants at construction
- Zero-alloc BCF encoder: pre-resolved typed handles (`ScalarInfoHandle<T>`, `GtFormatHandle`, etc.) write directly into BCF buffers in hot loops
- Forkable readers: shared index and header, `fork()` gives a fresh file handle per thread with no lock contention

## License

MIT/Apache-2.0
