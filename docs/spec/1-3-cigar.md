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
