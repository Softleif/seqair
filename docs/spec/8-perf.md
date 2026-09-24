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

r[perf.cigar_cursor]
The pileup engine MUST resolve each active read's column from an incremental cursor over the read's CIGAR — htslib `resolve_cigar2`'s shape: the current reference-consuming op with its reference and query start — advanced forward as the pileup moves right, NOT from a stateless per-column lookup. The indel anchored at the end of an op (the insertion or deletion that follows it, skipping `P`) MUST be resolved once when the cursor steps onto the op, not re-derived at every column. A column inside a match op with nothing anchored MUST cost one compare and one add, whatever the CIGAR. The cursor reads the ops from the store's CIGAR slab by index; it MUST NOT copy them. Its answers MUST equal `CigarMapping::pos_info_at` and `CigarMapping::deletion_after_at` for every non-decreasing sequence of positions.

The stateless lookup paid on every column what only changes at op boundaries: a scan (or binary search) for the covering op, then a second scan for the deletion after it, which together were a fifth of the engine's time; its per-read `SmallVec` of pre-computed ops was also most of the 128-byte active record the eviction pass moves (92 bytes without it).

r[perf.arena_capacity_hint+2]
The RecordStore MUST support a capacity hint so that the first region's allocation can be pre-sized based on an estimate (e.g., from the BAI index chunk sizes), avoiding repeated reallocation during the first `fetch_into`. This is provided by `RecordStore::with_byte_hint()`.
