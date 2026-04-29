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

r[cigar.aligned_pairs.types]
The `AlignedPair` enum MUST have the following variants:

- `Match { qpos: u32, rpos: Pos0 }` — M/=/X op, one yield per position
- `Insertion { qpos: u32, insert_len: u32 }` — I op, one yield per op (summary form)
- `Deletion { rpos: Pos0, del_len: u32 }` — D op, one yield per op (summary form)
- `RefSkip { rpos: Pos0, skip_len: u32 }` — N op, one yield per op (summary form)
- `SoftClip { qpos: u32, len: u32 }` — S op, one yield per op; yield controlled by options
- `Padding` — P op; yield controlled by options
- `Unknown { code: u8, len: u32 }` — reserved op code (9–15); yield controlled by options
  Hard clips (H) MUST never be yielded.

r[cigar.aligned_pairs.iterator]
The `AlignedPairs` iterator MUST maintain `qpos` (query position) and `rpos` (reference position) cursors and advance them according to CIGAR semantics: M/=/X advances both, I/S advance qpos only, D/N advance rpos only, H/P advance neither. Match ops MUST expand to per-position yields; Insertion/Deletion/RefSkip/SoftClip MUST yield once per op (summary form). The iterator MUST be infallible once constructed — all validation happens at construction time when the CIGAR slice is obtained.

r[cigar.aligned_pairs.default_mode]
Default mode (`AlignedPairs::new`) MUST yield only `Match`, `Insertion`, `Deletion`, and `RefSkip` variants. `SoftClip`, `Padding`, and `Unknown` MUST be silently skipped.

r[cigar.aligned_pairs.options]
Builder methods `.with_soft_clips()` and `.full()` MUST control which additional variants are yielded:

- `.with_soft_clips()` adds `SoftClip` to the default set
- `.full()` adds `SoftClip`, `Padding`, and `Unknown` to the default set
  These MUST be set-once methods that consume and return `Self`.

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
`OwnedBamRecord::aligned_pairs(&self)` MUST return `AlignedPairs` directly from the owned CIGAR without an intermediate `Result`.

r[cigar.aligned_pairs.complementary]
`AlignedPairs` and `CigarMapping::pos_info_at` are complementary, not redundant. `AlignedPairs` iterates per-CIGAR-op (record-walking); `pos_info_at` looks up by reference position (column-at-a-time). `PileupEngine` MUST continue to use `CigarMapping::pos_info_at` — do not refactor it onto `AlignedPairs`.

r[cigar.aligned_pairs.htslib_equivalence]
The expanded form (one `AlignedPair` per consumed base, as would be produced by `htslib::aligned_pairs_full`) MUST produce identical `(qpos, rpos)` pairs when compared against htslib for test fixtures. This is verified by property tests gated on the `rust-htslib` dev-dependency.

r[cigar.aligned_pairs.position_monotonicity]
Both `qpos` and `rpos` across yielded `AlignedPair` variants MUST be monotone non-decreasing. This is verified by property tests.
