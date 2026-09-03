# Genomic Position Types

## Motivation

Genomic positions appear throughout the codebase in multiple integer types (`i64`, `i32`, `u64`, `u32`, `NonZeroU32`) and two coordinate systems (0-based for BAM/BED/internal, 1-based for SAM/VCF/CRAM/user-facing). Manual `- 1` / `+ 1` conversions between systems are scattered across ~15 sites, creating two classes of bugs:

1. **Off-by-one errors**: A bare `pos - 1` doesn't communicate whether the subtraction is coordinate conversion or interval arithmetic. A future refactor could accidentally remove or duplicate it.
2. **Type confusion**: Position values flow through `u64` (public API) → `i64` (records) → `i32` (CRAM headers, CompactOp) with unchecked casts at each boundary. Errors like `.max(1) as u64 - 1` are fragile and opaque.

A phantom-typed newtype `Pos<S>` solves both: the type system distinguishes coordinate systems at compile time, and a single internal integer type (`u32`) eliminates gratuitous casts.

### Why `u32` capped at `i32::MAX`?

The SAM spec caps reference sequence length (`@SQ LN`) at `[1, 2^31-1]` (~2.1 billion) [SAM1 §1.3], and the BAM binary `l_ref` field is `uint32_t` with the same `< 2^31` constraint [SAM1 §4.2]. BAM binary `pos` is `int32_t`, giving a valid position range of `[0, 2^31-1]`. All current SAM/BAM positions therefore fit in `u32` (max ~4.3 billion) with room to spare, and capping at `i32::MAX` makes `as_i32()` infallible. Using `u32`:

- Halves position storage (4 bytes vs 8), improving cache density in hot structs (`SlimRecord`, `PileupAlignment`, `CompactOp`)
- Eliminates the i64↔i32 casts currently scattered through cigar.rs
- `Offset` (distance between positions) uses `i64` for safe intermediate arithmetic

> **Note on large genomes and future widening:** The CSI index format [CSI] uses `int64_t` for positions to support sequences larger than `2^29` bases. Enormous total genome sizes do exist — the fern _Tmesipteris oblanceolata_ has a ~160 Gbp genome ([Hidalgo et al., iScience 2024](<https://www.cell.com/iscience/fulltext/S2589-0042(24)01111-8>)) — but that size is spread across many chromosomes. However, individual chromosomes exceeding `2^31-1` already exist in practice: users have reported single chromosomes of ~2.36 Gbp that BAM simply cannot represent (`[E::bam_write1] Positional data is too large for BAM format`). This is an active unresolved issue in the SAM/BAM spec ([samtools/hts-specs#655](https://github.com/samtools/hts-specs/issues/655)); htslib already uses 64-bit coordinates internally for SAM, and CRAM 4.0 (still in draft) is planned to support them, but standard BAM and CRAM remain limited to `int32_t`. The proposed `PN`/`PO` header tags for stitching split chromosomes were never standardized. For seqair, `u32` is correct today because BAM enforces the `2^31-1` cap at the format level — but if support for SAM-only large-chromosome files or future CRAM 4.0 is ever added, `Pos<S>` would need to widen to `u64`.

> **Sources:** Coordinate conventions follow [SAM1] §1.4 (1-based POS), §4.2 (0-based BAM binary POS). Reference length limit from [SAM1] §1.3 (`@SQ LN` range `[1, 2^31-1]`) and §4.2 (`l_ref: uint32_t < 2^31`). CSI position width from [CSI]. See [References](./99-references.md).

## Position type

r[pos.type]
`Pos<S>` MUST be a `#[repr(transparent)]` newtype over `u32`, parameterized by a coordinate system marker `S`. The phantom type parameter MUST be zero-sized (`PhantomData<S>`), ensuring `Pos<S>` has identical layout and ABI to `u32`.

r[pos.size]
`size_of::<Pos<S>>()` MUST equal `size_of::<u32>()` (4 bytes) for all `S`.

r[pos.systems]
Two coordinate system markers MUST be defined:

- `Zero`: 0-based coordinates (BAM binary, BED, internal engine)
- `One`: 1-based coordinates (SAM text, VCF, CRAM, user-facing)

r[pos.incompatible]
`Pos<Zero>` and `Pos<One>` MUST be distinct types. The compiler MUST reject comparison, assignment, or arithmetic between positions of different systems.

## Construction

r[pos.zero_new]
`Pos::<Zero>::new(u32) -> Option<Self>` MUST return `None` if the value exceeds `i32::MAX`. All values in `0..=i32::MAX` are valid 0-based positions.

r[pos.one_new]
`Pos::<One>::new(u32) -> Option<Self>` MUST return `None` if the value is 0 or exceeds `i32::MAX`. The value 0 is not a valid 1-based position.

r[pos.try_from]
Fallible `TryFrom` trait impls MUST be provided for common integer types, all returning `Result<Self, PosOverflow>`:

- `TryFrom<i32>`: rejects negative values (for `Zero`) or values < 1 (for `One`).
- `TryFrom<u32>`: delegates to `new()`.
- `TryFrom<i64>`: rejects negative values (for `Zero`) or values < 1 (for `One`), and values exceeding `i32::MAX`.
- `TryFrom<u64>`: rejects values exceeding `i32::MAX` (and 0 for `One`).

These and `new()` are the ONLY entry points for external integer values. Raw `as` casts from integers to `Pos` MUST NOT be possible outside the module.

## Conversion

r[pos.to_one_based]
`Pos<Zero>::to_one_based()` MUST return `Result<Pos<One>, PosOverflow>` by adding 1 to the raw value. This is fallible: a 0-based position of `i32::MAX` would require a 1-based value of `i32::MAX + 1`, which exceeds the valid range for `Pos<One>`.

r[pos.to_zero_based]
`Pos<One>::to_zero_based()` MUST return `Pos<Zero>` by subtracting 1 from the raw value. This is infallible because the minimum 1-based value is 1, yielding 0-based 0.

r[pos.explicit_conversion]
Coordinate system conversion MUST only occur through `to_one_based()` and `to_zero_based()`. There MUST NOT be implicit `From`/`Into` conversions between systems.

## Accessors

r[pos.as_i32]
`pos.as_i32() -> i32` MUST return the raw value as `i32`.

r[pos.as_usize]
`pos.as_usize() -> usize` MUST return the raw value as `usize` for array indexing.

r[pos.as_i64]
`pos.as_i64() -> i64` MUST return the raw value as `i64` for wider arithmetic.

## Arithmetic

r[pos.offset]
`Offset` MUST be a `#[repr(transparent)]` newtype over `i64` representing a signed distance between positions.

r[pos.sub_pos]
`Pos<S> - Pos<S> = Offset`: subtracting two positions of the same system MUST produce an `Offset`. Subtracting positions of different systems MUST be a compile error.

r[pos.add_offset]
`Pos<S>::checked_add_offset(Offset) -> Option<Pos<S>>`: adding an offset to a position MUST produce a position in the same system. MUST return `None` if the result is negative or exceeds `i32::MAX`.

r[pos.sub_offset]
`Pos<S>::checked_sub_offset(Offset) -> Option<Pos<S>>`: subtracting an offset from a position MUST produce a position in the same system. Same bounds checking as `checked_add_offset`.

r[pos.no_add_pos]
Adding two positions (`Pos<S> + Pos<S>`) MUST be a compile error. Adding positions is not a meaningful operation.

## Reference intervals

Every reader and the pileup engine work on the same kind of interval, and the
one number worth being pedantic about is where it stops. Half-open ends and
inclusive ends both appear in the file formats this library reads (htslib's
`bam_endpos` is exclusive, CRAM's slice spans are inclusive), so seqair picks
one internally and converts at the format boundary.

r[interval.inclusive_ends]
A reference interval in seqair's API is 0-based with **both ends inclusive**:
`[start, end]` covers `end - start + 1` positions. This applies to the range
passed to `fetch_into`/`fetch_into_customized`, to `Segment`, and to the pileup
engine's region. A half-open interval MUST NOT be used at these boundaries;
where an external format needs one, the conversion happens in the format's
own reader.

r[interval.end_pos_inclusive]
A record's `end_pos` is the **last reference position it covers**
(`pos + ref_len - 1`), not one past it. A record with a CIGAR that consumes no
reference has `end_pos == pos`, matching the single position htslib reports as
`end = start + 1` in its exclusive convention.

r[interval.overlap_test]
A record overlaps `[start, end]` iff `record.pos <= end && record.end_pos >= start`.
Every reader MUST use exactly this test, so that the same query returns the same
records from BAM, SAM, and CRAM. A `Range<Pos0>` built from these values —
[`RecordStore::mate_overlap`](./3-record_store.md), for instance — is half-open
and therefore ends at `end_pos + 1`; that conversion belongs at the point the
`Range` is constructed.

## Traits

r[pos.derives]
`Pos<S>` MUST derive or implement: `Copy`, `Clone`, `PartialEq`, `Eq`, `PartialOrd`, `Ord`, `Hash`, `Debug`, `Display`.

r[pos.must_use]
All accessor methods and conversion methods MUST be `#[must_use]`.

## Niche optimization

r[pos.niche]
`Pos<S>` SHOULD use a niche-optimized inner type so that `Option<Pos<S>>` has the same size as `Pos<S>` (4 bytes). This is acceptable because the SAM/BAM spec caps reference sequence length at `2^31-1` (~2.1B) [SAM1 §1.3, §4.2], well below `u32::MAX` (~4.3B). **Status:** not yet implemented — currently `Option<Pos<S>>` is 8 bytes.
