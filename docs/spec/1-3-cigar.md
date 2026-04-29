# CIGAR Operations

A CIGAR string describes how a read aligns to the reference genome. It's a sequence of operations like "50M2I48M" meaning "50 bases match, 2 bases inserted in the read, 48 bases match." Each operation has a type and a length.

In BAM's binary format, CIGAR operations are packed as u32 values: the lower 4 bits encode the operation type (0–8) and the upper 28 bits encode the length. A 150bp read with a simple alignment has just one operation (`150M`), while reads spanning splice junctions or containing indels can have 10+ operations.

The key concept is which operations **consume** the reference and/or query (read) sequences:

- **Consumes reference**: the operation advances the position on the reference genome
- **Consumes query**: the operation advances the position within the read

This determines how to translate between reference positions and read positions — the core question the pileup engine asks thousands of times per second.

> **Sources:** [SAM1] §1.4 "The alignment section: mandatory fields" — CIGAR operations table with op codes, BAM numeric encodings, and "consumes query"/"consumes reference" flags. §4.2 "The BAM format" — packed u32 CIGAR encoding (`op_len << 4 | op`). See [References](./99-references.md).

## Operations

> _[SAM1] §1.4 "The alignment section: mandatory fields" — CIGAR operations table (M=0, I=1, D=2, N=3, S=4, H=5, P=6, ==7, X=8)_

> r[cigar.operations]
> The following CIGAR operations MUST be supported:
>
> - M (0): alignment match (consumes query and reference)
> - I (1): insertion (consumes query only)
> - D (2): deletion (consumes reference only)
> - N (3): reference skip (consumes reference only)
> - S (4): soft clip (consumes query only)
> - H (5): hard clip (consumes neither)
> - P (6): padding (consumes neither)
> - = (7): sequence match (consumes query and reference)
> - X (8): sequence mismatch (consumes query and reference)

## Matches and indels

These aggregate statistics are used downstream for read quality metrics (e.g. fraction of matching bases, indel rate).

r[cigar.matches_indels]
`calc_matches_indels(cigar)` MUST compute the total number of matching bases (sum of M operation lengths) and the total number of indel bases (sum of I and D operation lengths) from a packed CIGAR array.

## CIGAR index

For the pileup engine, the critical operation is: "at reference position X, what is the corresponding position in this read?" This requires walking the CIGAR to accumulate reference and query offsets. The `CigarIndex` precomputes this to avoid re-walking for every position.

r[cigar.index]
A `CigarIndex` MUST be constructable from a record's CIGAR and start position. It precomputes cumulative reference and query offsets per operation.

r[cigar.qpos_at]
`CigarIndex::qpos_at(ref_pos)` MUST return `Some(query_pos)` when the reference position falls within a query-consuming operation (M, =, X, I at the boundary), or `None` when the position falls within a deletion, reference skip, or outside the alignment.

r[cigar.qpos_bounds]
`qpos_at` MUST return `None` for positions outside the alignment's reference span. This applies to all mapping variants, including the linear fast-path. Positions before the alignment start or at/after the alignment end MUST return `None`.

r[cigar.compact_op_position_invariant]
`CompactOp::ref_start` is stored as `i32`. Since BAM positions are defined as `i32` values (max 2^31 - 1), this is sufficient. The construction of `CompactOp` MUST assert (via `debug_assert!`) that the computed reference position fits in `i32`.

r[cigar.qpos_accuracy]
`qpos_at` MUST produce identical results to htslib's `Alignment::qpos()` for every reference position within the alignment span. This is verified by comparison tests.

## Aligned pairs

The `AlignedPairs` iterator walks a record's CIGAR operations once per record and yields typed `AlignedPair` variants describing the aligned query/reference positions. This is analogous to htslib's `aligned_pairs` / `aligned_pairs_full`, but uses a typed enum instead of `(Option<i64>, Option<i64>)` tuples.

> **Runnable example:** `cargo run --example aligned_pairs_walk -- --input X.bam --reference X.fa --region chrZ:A-B` exercises all three layers (bare walk, `with_read`, `with_reference`) against a real BAM/FASTA pair. It also has stand-alone showcase functions for indel collection and `MatchKind`-based mismatch counting.

r[cigar.aligned_pairs.types]
The `AlignedPair` enum MUST have the following variants:

- `Match { qpos: u32, rpos: Pos0, kind: MatchKind }` — M/=/X op, one yield per position; `kind` distinguishes the source CIGAR op (see `r[cigar.aligned_pairs.match_kind]`)
- `Insertion { qpos: u32, insert_len: u32 }` — I op, one yield per op (summary form)
- `Deletion { rpos: Pos0, del_len: u32 }` — D op, one yield per op (summary form)
- `RefSkip { rpos: Pos0, skip_len: u32 }` — N op, one yield per op (summary form)
- `SoftClip { qpos: u32, len: u32 }` — S op, one yield per op; yield controlled by options
- `Padding { len: u32 }` — P op, one yield per op; yield controlled by options
- `Unknown { code: u8, len: u32 }` — reserved op code (9–15); yield controlled by options
  Hard clips (H) MUST never be yielded. The total in-memory size of `AlignedPair` MUST be ≤ 16 bytes; this is enforced by a compile-time `const _` assert.

r[cigar.aligned_pairs.match_kind]
The `MatchKind` enum MUST have three variants — `Match` (M, BAM op 0), `SeqMatch` (=, BAM op 7), and `SeqMismatch` (X, BAM op 8) — and MUST be carried inside every `AlignedPair::Match`. The qpos/rpos sequence yielded by the iterator MUST be identical regardless of `kind`; the field is informational only. This preserves the per-op distinction that htslib's `aligned_pairs_full` discards.

r[cigar.aligned_pairs.zero_length_ops]
The iterator MUST skip zero-length CIGAR ops uniformly. No `len: 0` summary variant (`Insertion`, `Deletion`, `RefSkip`, `SoftClip`, `Padding`, `Unknown`) is ever yielded, and zero-length M/=/X ops do not enter expansion. This mirrors htslib's CIGAR walker, whose per-op for-loops simply do not iterate at length 0.

r[cigar.aligned_pairs.iterator]
The `AlignedPairs` iterator MUST maintain `qpos` (query position) and `rpos` (reference position) cursors and advance them according to CIGAR semantics: M/=/X advances both, I/S advance qpos only, D/N advance rpos only, H/P advance neither. Match ops MUST expand to per-position yields; Insertion/Deletion/RefSkip/SoftClip MUST yield once per op (summary form). The iterator MUST be infallible once constructed — all validation happens at construction time when the CIGAR slice is obtained.

r[cigar.aligned_pairs.default_mode]
Default mode (`AlignedPairs::new`) MUST yield only `Match`, `Insertion`, `Deletion`, and `RefSkip` variants. `SoftClip`, `Padding`, and `Unknown` MUST be silently skipped.

r[cigar.aligned_pairs.options]
Builder methods `.with_soft_clips()` and `.full()` MUST control which additional variants are yielded:

- `.with_soft_clips()` adds `SoftClip` to the default set
- `.full()` adds `SoftClip`, `Padding`, and `Unknown` to the default set

These MUST consume `self` and return `Self`. Calling either more than once is permitted and idempotent — the second call has no additional effect because the underlying option flags are already set. The methods are also available on the layered `AlignedPairsWithRead` and `AlignedPairsWithRef` iterators, where they forward to the inner iterator.

r[cigar.aligned_pairs.qpos_semantics]
`qpos` MUST be 0-based within the full read sequence including soft clips. The first match after a leading soft clip MUST have `qpos` equal to the soft-clip length. This matches htslib behavior where `aligned_pairs_full()` reports `Some(clip_start)` for the first match after clips.

r[cigar.aligned_pairs.hard_clips]
Hard clips (H) MUST NOT advance `qpos` and MUST NOT yield any variant. This is consistent with htslib, which also excludes hard clips from `aligned_pairs_full()`.

r[cigar.aligned_pairs.insertion_qpos]
For `Insertion { qpos, insert_len }`, `qpos` MUST be the position of the _first_ inserted base (the query position before the I op advances).

r[cigar.aligned_pairs.deletion_rpos]
For `Deletion { rpos, del_len }`, `rpos` MUST be the reference position where the deletion begins (before the D op advances rpos).

r[cigar.aligned_pairs.slim_record]
`SlimRecord::aligned_pairs(&self, store: &RecordStore<U>)` MUST return `Result<AlignedPairs, RecordAccessError>` by first reading the CIGAR slice from the store. Iteration over the returned `AlignedPairs` is infallible.

r[cigar.aligned_pairs.owned_record]
`OwnedBamRecord::aligned_pairs(&self) -> Result<AlignedPairs, AlignedPairsError>` MUST validate before iteration starts: if `pos = None` (unmapped) AND the CIGAR is non-empty, it MUST return `AlignedPairsError::UnmappedWithCigar { cigar_ops }`. For mapped records (`pos = Some(_)`) and fully-unmapped records (`pos = None` + empty CIGAR), iteration is infallible. The reason: walking a non-empty CIGAR for an unmapped record without a base reference position would produce nonsense `rpos` values anchored at `Pos0::ZERO`; refusing is safer than silently corrupting output.

**Design note (non-verifiable):** `AlignedPairs` and `CigarMapping::pos_info_at` are complementary, not redundant. `AlignedPairs` iterates per-CIGAR-op (record-walking); `pos_info_at` looks up by reference position (column-at-a-time). `PileupEngine` continues to use `CigarMapping::pos_info_at` — do not refactor it onto `AlignedPairs`. (No `r[…]` because there is no per-call behavior to verify; this is a structural caveat.)

r[cigar.aligned_pairs.htslib_equivalence]
The expanded form (one `AlignedPair` per consumed base, as would be produced by `htslib::aligned_pairs_full`) MUST produce identical `(qpos, rpos)` pairs when compared against htslib for test fixtures. This is verified by property tests gated on the `rust-htslib` dev-dependency.

r[cigar.aligned_pairs.position_monotonicity]
Both `qpos` and `rpos` across yielded `AlignedPair` variants MUST be monotone non-decreasing. This is verified by property tests.

## Layered iterator views

The base [`AlignedPairs`] iterator yields position-only events. Two adapter layers attach per-event read and reference data on demand. The layering exists for pay-for-what-you-use cost on hot paths: bare iteration has no slab lookups; `with_read` adds two slab references; `with_reference` adds one window borrow. Each layer's yield type is a distinct enum, so the compiler enforces the presence of the relevant data.

r[cigar.aligned_pairs.with_read.types]
`AlignedPairWithRead` MUST have one variant per `AlignedPair` variant. `Match` MUST gain `query: Base` and `qual: BaseQuality` (the read base and Phred score at `qpos`). `Insertion` and `SoftClip` MUST gain `query: &[Base]` and `qual: &[BaseQuality]` slices borrowed from the record's seq/qual slabs (the inserted/clipped run, not the whole record). `Deletion`, `RefSkip`, `Padding`, `Unknown` MUST pass through unchanged (no read data applies).

r[cigar.aligned_pairs.with_read.iterator]
`AlignedPairs::with_read(seq, qual) -> Result<AlignedPairsWithRead<'cigar, 'read>, AlignedPairsError>` MUST validate length invariants up front: the iterator's total query-consuming CIGAR length MUST equal `seq.len()` (otherwise `CigarSeqLengthMismatch`), and `qual.len()` MUST be either zero (BAM missing-qual sentinel) or equal to `seq.len()` (otherwise `SeqQualLengthMismatch`). Once construction succeeds, every yielded `Match`/`Insertion`/`SoftClip` index is provably in-bounds; iteration is infallible. Lifetimes `'cigar` and `'read` are independent so callers can compose borrows from different sources. The number of yielded events MUST be identical to the underlying `AlignedPairs` (the layer enriches; it does not filter or split events). `with_soft_clips()` and `full()` MUST be available on `AlignedPairsWithRead`, forwarding to the inner iterator.

r[cigar.aligned_pairs.with_read.slim_record]
`SlimRecord::aligned_pairs_with_read(store) -> Result<AlignedPairsWithRead, AlignedPairsError>` MUST be a one-shot helper equivalent to `self.aligned_pairs(store)?.with_read(self.seq(store)?, self.qual(store)?)`. It exists so callers do not have to fan out three slab accesses by hand. Errors from the slab fetch surface as `AlignedPairsError::Access`; mismatched lengths surface as `CigarSeqLengthMismatch` / `SeqQualLengthMismatch`.

r[cigar.aligned_pairs.with_reference.types]
`AlignedPairWithRef` MUST mirror `AlignedPairWithRead`'s variants. `Match` MUST gain `ref_base: Option<Base>` (`None` if `rpos` is outside the loaded `RefSeq` window; `Some(Base::Unknown)` for in-window N's). `Deletion` MUST gain `ref_bases: Option<&[Base]>` (`None` if any position in the deletion span is outside the window — partial overlap does not pad). All other variants pass through unchanged.

r[cigar.aligned_pairs.with_reference.iterator]
`AlignedPairsWithRead::with_reference(ref_seq) -> AlignedPairsWithRef<'cigar, 'read, 'ref_seq>` MUST reuse the existing `RefSeq` type from the pileup engine — no new wrapper, no copy. All three lifetimes are independent. `with_soft_clips()` and `full()` MUST also be available here, forwarding through both layers.

## Matches-only adapters

For methylation and SNV calling, callers usually only care about `Match` events — indels and clips are noise. Each layer provides a `matches_only()` adapter that filters to Match events and yields a flat struct (no enum dispatch).

r[cigar.aligned_pairs.matches_only.types]
Three flat structs MUST be provided, one per layer:

- `MatchPosition { qpos: u32, rpos: Pos0, kind: MatchKind }` — bare layer, equivalent to pysam's `aligned_pairs(matches_only=True)` plus the `kind` distinction htslib loses.
- `MatchedBase { qpos, rpos, kind, query: Base, qual: BaseQuality }` — with-read layer.
- `MatchedRef { qpos, rpos, kind, query, qual, ref_base: Option<Base> }` — with-reference layer; the most common shape for methylation and per-record SNV calling.

None of these carry lifetimes — `Base`, `BaseQuality`, `MatchKind`, and `Pos0` are all `Copy`, so the iterators yield owned values.

r[cigar.aligned_pairs.matches_only.bare]
`AlignedPairs::matches_only() -> MatchesOnly` MUST yield only `Match` events, dropping indels and clips silently. Field values MUST equal the corresponding `AlignedPair::Match` fields one-to-one.

r[cigar.aligned_pairs.matches_only.with_read]
`AlignedPairsWithRead::matches_only() -> MatchedBases` MUST yield only `Match` events with read base/qual attached.

r[cigar.aligned_pairs.matches_only.with_reference]
`AlignedPairsWithRef::matches_only() -> MatchedRefs` MUST yield only `Match` events with read base/qual AND `ref_base` attached.

## Iterator invariants

r[cigar.aligned_pairs.iterator_traits]
The bare and layered iterators (`AlignedPairs`, `AlignedPairsWithRead`, `AlignedPairsWithRef`) MUST implement `Clone`, `ExactSizeIterator`, and `FusedIterator`. `size_hint()` MUST return `(n, Some(n))` where `n` is the exact number of remaining events (the iterator can compute this from the unprocessed CIGAR plus the current expansion state).

The matches-only adapters (`MatchesOnly`, `MatchedBases`, `MatchedRefs`) MUST implement `Clone` and `FusedIterator`. Because they filter, `size_hint()` MUST return `(0, Some(upper))` where `upper` is the inner iterator's upper bound — they're not exact-sized.

r[cigar.aligned_pairs.with_read.owned_record]
`OwnedBamRecord::aligned_pairs_with_read()` MUST be a one-shot helper symmetric with `SlimRecord::aligned_pairs_with_read`, equivalent to `self.aligned_pairs()?.with_read(&self.seq, &self.qual)`. It validates the same length invariants and returns `AlignedPairsError` on failure.

## NM and MD recompute

The `NM` (edit distance) and `MD` (mismatch positions) tags are derivable from a record's CIGAR walk against a reference. `AlignedPairsWithRef` exposes them as consuming methods so callers can compute on demand without re-walking.

r[cigar.aligned_pairs.nm.compute]
`AlignedPairsWithRef::nm() -> u32` MUST consume the iterator and return the SAM-spec edit distance: `M`-op mismatches (decided by `query != ref_base`) + `=` always 0 + `X` always 1 per position + `I` length + `D` length. `N` and clips contribute 0. The function is best-effort: positions where the reference is outside the loaded `RefSeq` window MUST be skipped (count 0) rather than erroring. Callers needing strict computation MUST verify `RefSeq` coverage before calling. Comparison semantics: `Base::Unknown` (`N`) is treated as a mismatch against any base, including another `N`.

r[cigar.aligned_pairs.md.compute]
`AlignedPairsWithRef::md() -> Result<Vec<u8>, NmMdError>` MUST consume the iterator and return the SAM-spec MD tag bytes: alternating digit runs (counts of consecutive matching positions) and either single ref bases (mismatches) or `^` followed by ref bases (deletions). The output MUST start AND end with a digit (possibly `0`). Two adjacent deletions MUST be separated by a `0` to remain unambiguous. Insertions and soft clips MUST NOT appear. Reference skips (`N`) MUST NOT contribute to MD (matches `samtools calmd` behavior). Returns `NmMdError::MissingReference` if any required reference position is outside the loaded `RefSeq` window — MD encodes ref bases verbatim and partial coverage would corrupt the output.

r[cigar.aligned_pairs.nm_md.consistency]
For the same record + reference + iterator configuration, `nm()` and `md()` MUST agree: the computed `NM` MUST equal `(MD letters outside ^) + (MD bases inside ^) + (CIGAR insertion lengths)`. This invariant is verified by a property test.
