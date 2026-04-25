# Pileup Engine

A pileup is a column-by-column view of all reads aligned to a reference region. At each reference position, the pileup shows which reads are present and what base each read has at that position. This is the fundamental data structure for variant calling: to determine whether position X has a mutation, you look at the pileup column at X and count how many reads show the reference base vs. an alternate base.

The pileup engine takes a set of BAM records (stored in a `RecordStore`) and iterates over reference positions from left to right, maintaining an "active set" of reads that overlap the current position. At each position, it yields a `PileupColumn` containing the aligned reads with their query positions (the position within each read that corresponds to this reference position).

> **Sources:** The pileup algorithm is seqair-specific design, verified for correctness against [htslib] `bam_plp_auto`. FLAG bits referenced here (0x4 unmapped, 0x10 reverse strand) are defined in [SAM1] §1.4. CIGAR consume semantics (consumes query/reference) are from [SAM1] §1.4. See [References](./99-references.md).

## Column iteration

r[pileup.position_iteration]
The engine MUST iterate over every reference position from the region start to the region end (inclusive), yielding a `PileupColumn` at each position where at least one read is active.

r[pileup.active_set]
At each position, the active set contains all records whose alignment span `[pos, end_pos]` includes the current position. Records MUST enter the active set when the current position reaches their `pos` and MUST be evicted when the current position exceeds their `end_pos`.

> r[pileup.column_contents]
> Each `PileupColumn` MUST provide:
>
> - the reference position
> - the number of active alignments (depth)
> - access to each alignment with its query position (`qpos`) and a `RecordRef`

## Filtering

r[pileup.read_filter]
The engine MUST support a per-read filter function that is evaluated once per record when it enters the active set (NOT re-evaluated at each position). Records that fail the filter MUST be excluded from all subsequent columns.

r[pileup.max_depth]
The engine MUST support a maximum depth setting. When more records pass the filter than `max_depth` allows, excess records MUST be dropped (preferring records already in the active set).

r[pileup.max_depth_per_position]
Max depth MUST be enforced per-position, not globally at arena-load time. A read that is rejected at one high-coverage position may still be included at an adjacent lower-coverage position. This matches htslib's behavior where `MAXCNT` is checked against the current column's depth.

## Query position

The query position (`qpos`) maps a reference position back into a read's coordinate system. For example, if a read starts at reference position 100 and has CIGAR `50M`, then at reference position 125, `qpos` = 25. For complex CIGARs with insertions or deletions, the mapping requires walking the CIGAR operations (see `cigar.md`).

r[pileup.qpos]
For each alignment in a column, the engine MUST provide `qpos`: the position within the read's query sequence that aligns to the current reference position. This MUST be computed using the `CigarIndex`.

r[pileup.qpos_none]
When the current reference position falls within a deletion or reference skip in the read's CIGAR, that alignment MUST be excluded from the column (no qpos exists).

## Edge cases

r[pileup.empty_positions_skipped]
Positions with zero active reads MUST be skipped (no column yielded). When there is a gap between reads, the engine MUST jump to the next read's start position rather than iterating through empty positions one by one.

r[pileup.trailing_empty_termination]
When the active set is empty and all records in the store have been consumed (no more records ahead), the engine MUST terminate iteration immediately rather than continuing through remaining empty positions to the region end. Callers that need zero-depth columns for trailing positions MUST generate them externally.

r[pileup.unmapped_excluded]
Reads with the unmapped flag (0x4) MUST be excluded from the pileup. htslib's `bam_plp_push` rejects unmapped reads before they enter the buffer.

r[pileup.zero_refspan_reads]
Reads with zero reference-consuming CIGAR operations (pure soft-clip, insertion-only) have `end_pos == pos` and MUST appear in exactly one pileup column at their `pos`. They MUST NOT be skipped or cause off-by-one errors.

r[pileup.soft_clip_at_position]
When a reference position falls within the soft-clipped portion of a read's alignment, the read is NOT active at that position (soft clips do not consume reference). The read only becomes active at its first reference-consuming position.

## Per-record extras in pileup

r[pileup.extras.generic_param]
`PileupEngine` MUST accept a type parameter `U` (default `()`) matching its `RecordStore<U>`. When `U` is `()`, the engine MUST behave identically to the non-generic version — column construction and all non-extras-specific APIs MUST be unchanged.

r[pileup.extras.constructor_accepts_any_u]
`PileupEngine::new` MUST accept `RecordStore<U>` for any `U`. Callers MUST compute extras on the store via `RecordStore::apply_customize` with a [`CustomizeRecordStore`](../../crates/seqair/src/bam/record_store.rs) value before constructing the engine. There is no engine-level customize method — the constructor is the only entry point.

r[pileup.lending_iterator]
`PileupEngine<U>` MUST expose iteration via a single lending method `pileups(&mut self) -> Option<PileupColumn<'_, U>>`. The returned `PileupColumn<'store, U>` MUST borrow the engine's store for the duration of its use. `PileupEngine` MUST NOT implement `Iterator` — the store borrow held by each column would be incompatible with `Iterator::next`'s `&mut self` contract. Callers MUST use a `while let Some(col) = engine.pileups() { ... }` loop.

r[pileup.alignment_view]
`PileupColumn<'store, U>::alignments(&self)` MUST yield `AlignmentView<'_, 'store, U>` items. `AlignmentView` MUST deref to `PileupAlignment` so existing field and method access on alignments is unchanged. `AlignmentView` MUST additionally provide `extra(&self) -> &'store U`, `qname(&self) -> &'store [u8]`, and `aux(&self) -> &'store [u8]` methods that read from the borrowed store without additional arguments. `PileupColumn` MUST also expose `raw_alignments()` yielding `&PileupAlignment` for callers that do not need store access.

r[pileup.extras.recover_store]
`Readers::recover_store` MUST accept `PileupEngine<U>` for any `U`. It MUST strip extras via `RecordStore::strip_extras` before storing the recovered `RecordStore<()>` for reuse.

## Compatibility

r[pileup.htslib_compat]
For any BAM file and region, the pileup engine MUST produce identical columns to htslib's `bam_plp_auto`-based pileup: same positions, same depth at each position, and same `qpos` values per alignment. This is verified by comparison tests.
