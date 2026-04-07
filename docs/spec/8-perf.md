# Performance

The pileup engine is the innermost loop of Rastair's processing pipeline: for a typical 30× whole-genome BAM, it processes ~3 billion reference positions, each with ~30 aligned reads. These rules specify performance-sensitive behaviors to avoid unnecessary overhead in this hot path.

> **Sources:** Seqair-specific performance rules with no upstream spec counterpart. Implementation strategies are guided by profiling results and cross-checked against [htslib] behaviour. See [References](./99-references.md).

## Allocation avoidance

r[perf.reuse_alignment_vec+2]
The pileup engine allocates a `Vec<PileupAlignment>` per column using `Vec::with_capacity(active.len())`. Reusing a buffer via clone was measured to be slower than fresh allocation due to the copy cost of entries. The allocator efficiently reuses recently-freed blocks of the same size.

r[perf.avoid_redundant_arena_get+2]
When a record enters the pileup active set, the engine MUST retrieve the record data once and cache what it needs (cigar, flags, strand, mapq, etc.) in the `ActiveRecord` — NOT perform repeated store lookups for the same record in the hot loop.

r[perf.no_sorted_indices]
Since BAM records arrive in coordinate-sorted order from `fetch_into`, the pileup engine MUST NOT sort record indices. It MUST iterate records in arena order directly.

r[perf.cigar_no_to_vec]
Building a `CigarIndex` MUST NOT clone the CIGAR bytes via `.to_vec()`. The CIGAR bytes live in the arena slab for the lifetime of the region; the `CigarIndex` MUST borrow or reference them by offset, not own a copy.

r[perf.precompute_matches_indels]
Matches and indels counts MUST be computed once per record during decode and stored in `BamRecord`, NOT recomputed from CIGAR at every pileup position.

r[perf.cigar_binary_search]
For records with more than 4 CIGAR operations, `qpos_at` SHOULD use binary search on `ref_start` instead of linear scan. For 1–4 ops, linear scan is acceptable.

r[perf.arena_capacity_hint+2]
The RecordStore MUST support a capacity hint so that the first region's allocation can be pre-sized based on an estimate (e.g., from the BAI index chunk sizes), avoiding repeated reallocation during the first `fetch_into`. This is provided by `RecordStore::with_byte_hint()`.
