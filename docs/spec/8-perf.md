# Performance

The pileup engine is the innermost loop of Rastair's processing pipeline: for a typical 30× whole-genome BAM, it processes ~3 billion reference positions, each with ~30 aligned reads. These rules specify performance-sensitive behaviors to avoid unnecessary overhead in this hot path.

> **Sources:** Seqair-specific performance rules with no upstream spec counterpart. Implementation strategies are guided by profiling results and cross-checked against [htslib] behaviour. See [References](./99-references.md).

## Allocation avoidance

r[perf.reuse_alignment_vec+3]
The pileup engine MUST hold one `Vec<PileupAlignment>` on the engine and reuse it for every column: `clear`, one `reserve` for the whole column, then a straight write per entry into the spare capacity. It MUST NOT allocate a fresh `Vec` per column, and it MUST NOT `push` per entry.

The earlier version of this rule said the opposite — allocate per column, because "reusing a buffer via clone was measured to be slower". That measured the wrong thing: it compared a fresh allocation against *cloning* into a reused buffer, which pays the copy the clone implies, rather than against clearing and writing in place. Against the in-place form the reuse wins twice, once on the allocator and once in the loop: `Vec::push` cannot hoist its capacity compare, because it may reallocate, so it repeats compare, store and length update per read per column — 1.5 % of a variant caller's worker CPU spent rediscovering a bound the loop already knows, namely at most one entry per active record.

`internal_buf_retains_capacity` pins the reuse: it is the observable half, and an engine that went back to allocating per column would fail it.

r[perf.avoid_redundant_arena_get+2]
When a record enters the pileup active set, the engine MUST retrieve the record data once and cache what it needs (cigar, flags, strand, mapq, etc.) in the `ActiveRecord` — NOT perform repeated store lookups for the same record in the hot loop.

r[perf.no_sorted_indices]
Since BAM records arrive in coordinate-sorted order from `fetch_into`, the pileup engine MUST NOT sort record indices. It MUST iterate records in arena order directly.

r[perf.cigar_no_to_vec]
Building a `CigarMapping` MUST NOT clone the CIGAR ops via `.to_vec()`. The typed CIGAR ops live in the arena slab for the lifetime of the region (see `r[record_store.slim_record.field_getters]`); `CigarMapping::new` MUST take a borrowed `&[CigarOp]` and either pre-extract a small `CompactOp` array (Complex path) or compute summary state without copying the ops (Linear path).

r[perf.precompute_matches_indels]
Matches and indels counts MUST be computed once per record during decode and stored in `SlimRecord`, NOT recomputed from CIGAR at every pileup position.

r[perf.cigar_binary_search]
For records with more than 4 CIGAR operations, `qpos_at` SHOULD use binary search on `ref_start` instead of linear scan. For 1–4 ops, linear scan is acceptable.

r[perf.arena_capacity_hint+2]
The RecordStore MUST support a capacity hint so that the first region's allocation can be pre-sized based on an estimate (e.g., from the BAI index chunk sizes), avoiding repeated reallocation during the first `fetch_into`. This is provided by `RecordStore::with_byte_hint()`.
