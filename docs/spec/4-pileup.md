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

Read filtering is the responsibility of `CustomizeRecordStore::keep_record` (see [Record Store](./3-record_store.md)) at fetch time. By the time records reach the pileup engine they have already been filtered, so `PileupEngine` MUST NOT expose a separate read-filter API. The previous engine-level `set_filter` callback was removed because it duplicated `keep_record` while allocating a boxed closure and forcing pileup to be `!Send`; push-time filtering is strictly better (zero slab waste, no `'static` bound, sees the whole `SlimRecord` rather than just `(flags, aux)`).

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

r[pileup.empty_seq_unknown_base]
Mapped records with `SEQ=*` (`seq_len == 0`) MUST appear in every pileup column they overlap, with `base = Base::Unknown` and `qual = BaseQuality::UNAVAILABLE` at every M/=/X/I CIGAR op. htslib's pileup includes such records — `bam_seqi` reads beyond the empty SEQ buffer and the `0xFF` qual sentinel decodes to `N`. Silently dropping them would produce a different `depth()` and `match_depth()` than htslib for any BAM containing secondary alignments without sequence (a common SAM/BAM convention). Insertions on empty-SEQ records use the same Unknown/UNAVAILABLE pair for the anchor base.

r[pileup.zero_refspan_reads]
Reads with zero reference-consuming CIGAR operations (pure soft-clip, insertion-only) have `end_pos == pos` and MUST appear in exactly one pileup column at their `pos`. They MUST NOT be skipped or cause off-by-one errors.

r[pileup.soft_clip_at_position]
When a reference position falls within the soft-clipped portion of a read's alignment, the read is NOT active at that position (soft clips do not consume reference). The read only becomes active at its first reference-consuming position.

## Soft-clip overhang mode (opt-in)

By default soft clips are invisible to the pileup (r[pileup.soft_clip_at_position]). The overhang mode opts into emitting the clipped fringe bases immediately adjacent to each alignment, so a caller can recover evidence the aligner trimmed — for example a TAPS methylation base (a `T` over a CpG `C`) clipped because it looked like a mismatch at the read end. The mode is off by default to preserve htslib parity.

r[pileup.soft_clip_overhang.default_off]
The soft-clip overhang MUST default to 0. With overhang 0 the engine MUST behave exactly as r[pileup.soft_clip_at_position] specifies — soft-clipped positions yield no alignment — and the column stream MUST remain byte-identical to the htslib-parity behavior (r[pileup.htslib_compat]).

r[pileup.soft_clip_overhang.activation]
When the overhang is `n > 0`, a record MUST enter the active set up to `n` reference positions before its `pos` and MUST remain active up to `n` positions after its `end_pos`, so the columns flanking its alignment can carry its clipped fringe bases. The activation order MUST remain valid for the store's `pos`-sorted records (`pos - n` is monotonic in `pos`).

r[pileup.soft_clip_overhang.emit]
At a flanking position, a record whose soft clip on that side extends far enough MUST be emitted as a `PileupOp::SoftClip { qpos, base, qual }`, where `qpos` is `CigarMapping::soft_clip_qpos_at` (r[cigar.soft_clip_qpos]). The `base`, `qual`, and `qpos` accessors MUST resolve `SoftClip` identically to `Match`, and `is_soft_clip()` MUST distinguish it. A record with no soft clip on the flanking side MUST contribute nothing there.

r[pileup.soft_clip_overhang.additive]
Enabling the overhang MUST be purely additive: for the same records, the sequence of non-`SoftClip` column entries — positions, per-alignment `(record, op, qpos)`, and the depth excluding soft clips — MUST be identical to an overhang-0 run. The mode only adds `SoftClip` entries; it MUST NOT alter or remove any aligned/deletion/refskip entry, and MUST NOT emit columns whose only contents would be zero alignments.

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
`Readers::pileup` MUST return a guard type (`PileupGuard<'_, E::Extra>`) that derefs to `PileupEngine<E::Extra>` and whose `Drop` impl moves the underlying `RecordStore` back into the originating `Readers`, retaining its allocated capacity for the next pileup call. Recovery MUST happen on every drop path — end of scope, `?`-propagated error mid-iteration, `break` out of the loop — so callers do not need an explicit recover step. There MUST NOT be a separate `Readers::recover_store` method on the public API; if extras-stripping is needed it happens inside the guard's drop path.

## Mate links in columns

The pileup engine does not deduplicate overlapping mates itself (see
[Overlapping Pair Deduplication](./4-pileup-dedup.md) for why that has to happen
after per-position quality filtering). What it MUST do is hand the consumer the
mate link computed once per store by
[`RecordStore::link_mates`](./3-record_store.md), so the consumer never has to
hash or compare a qname per column.

r[pileup.mate_link_cache]
When a record enters the active set, the engine MUST cache its mate index, its
qname hash, and its mate-overlap interval, so that column construction reads
flat fields only and never looks the mate up in the store. `PileupAlignment`
MUST expose `mate_idx() -> Option<u32>` (the store index of the mate, `None`
when unlinked), `in_mate_overlap() -> bool` (true when this column's position
lies inside the mate-overlap interval), and `qname_hash() -> u64`. When the
soft-clip overhang is enabled, the cached interval MUST be widened by the
overhang at both ends, so a projected fringe base is reported as being inside
the overlap: it is the same molecule as its mate's aligned base at that
position, which is the only thing a consumer's dedup asks.

r[pileup.column_record_order]
The alignments of a column MUST be ordered by ascending `record_idx`, since the
active set preserves store order and the store is position-sorted. Consumers MAY
rely on this to look up another alignment of the same column by record index.

r[pileup.column_find_record]
`PileupColumn::find_record(record_idx)` MUST return the `AlignmentView` of that
record in this column, or `None` when the record has no entry here. It MUST run
in `O(log depth)` by binary search over
[`pileup.column_record_order`](#pileupcolumn_record_order). This is how a
consumer reaches the other mate of an overlapping pair while walking the column
exactly once.

## Compatibility

r[pileup.htslib_compat]
For any BAM file and region, the pileup engine MUST produce identical columns to htslib's `bam_plp_auto`-based pileup: same positions, same depth at each position, and same `qpos` values per alignment. This is verified by comparison tests.
