# Genomic Position Types

## Motivation

Genomic positions appear throughout the codebase in multiple integer types (`i64`, `i32`, `u64`, `u32`, `NonZeroU32`) and two coordinate systems (0-based for BAM/BED/internal, 1-based for SAM/VCF/CRAM/user-facing). Manual `- 1` / `+ 1` conversions between systems are scattered across ~15 sites, creating two classes of bugs:

1. **Off-by-one errors**: A bare `pos - 1` doesn't communicate whether the subtraction is coordinate conversion or interval arithmetic. A future refactor could accidentally remove or duplicate it.
2. **Type confusion**: Position values flow through `u64` (public API) → `i64` (records) → `i32` (CRAM headers, CompactOp) with unchecked casts at each boundary. Errors like `.max(1) as u64 - 1` are fragile and opaque.

A phantom-typed newtype `Pos<S>` solves both: the type system distinguishes coordinate systems at compile time, and a single internal integer type (`u32`) eliminates gratuitous casts.

### Why `u32` and not `i64`?

The SAM spec caps reference sequence length (`@SQ LN`) at `[1, 2^31-1]` (~2.1 billion) [SAM1 §1.3], and the BAM binary `l_ref` field is `uint32_t` with the same `< 2^31` constraint [SAM1 §4.2]. BAM binary `pos` is `int32_t`, giving a valid position range of `[0, 2^31-1]`. All current SAM/BAM positions therefore fit in `u32` (max ~4.3 billion) with room to spare. Using `u32`:
- Halves position storage (4 bytes vs 8), improving cache density in hot structs (`SlimRecord`, `PileupAlignment`, `CompactOp`)
- Eliminates the i64↔i32 casts currently scattered through cigar.rs
- `Offset` (distance between positions) uses `i64` for safe intermediate arithmetic

> **Note:** The CSI index format [CSI] uses `int64_t` for positions to support sequences larger than `2^29` bases. If seqair ever needs CSI-indexed files with sequences > `2^31-1` bases, `Pos<S>` would need to widen to `u64`. For now, SAM/BAM's `2^31-1` cap makes `u32` sufficient.

> **Sources:** Coordinate conventions follow [SAM1] §1.4 (1-based POS), §4.2 (0-based BAM binary POS). Reference length limit from [SAM1] §1.3 (`@SQ LN` range `[1, 2^31-1]`) and §4.2 (`l_ref: uint32_t < 2^31`). CSI position width from [CSI]. See [references.md](99-references.md).

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
`Pos::<Zero>::new(u32)` MUST be an infallible constructor accepting any `u32` value. All `u32` values are valid 0-based positions.

r[pos.one_new]
`Pos::<One>::new(u32)` MUST reject 0. In debug builds, it MUST panic. The value 0 is not a valid 1-based position.

r[pos.try_from]
Fallible constructors MUST be provided for wider integer types:
- `try_from_i64(i64) -> Option<Self>`: rejects negative values (for `Zero`) or values < 1 (for `One`), and values exceeding `u32::MAX`.
- `try_from_u64(u64) -> Option<Self>`: rejects values exceeding `u32::MAX`.
- `try_from_i32(i32) -> Option<Self>` (for `One`): rejects values < 1.

These are the ONLY entry points for external integer values. Raw `as` casts from integers to `Pos` MUST NOT be possible outside the module.

## Conversion

r[pos.to_one_based]
`Pos<Zero>::to_one_based()` MUST return `Pos<One>` by adding 1 to the raw value. This is infallible because the maximum valid BAM position (`i32::MAX` ≈ 2.1B) plus 1 fits in `u32` (max ≈ 4.3B).

r[pos.to_zero_based]
`Pos<One>::to_zero_based()` MUST return `Pos<Zero>` by subtracting 1 from the raw value. This is infallible because the minimum 1-based value is 1, yielding 0-based 0.

r[pos.explicit_conversion]
Coordinate system conversion MUST only occur through `to_one_based()` and `to_zero_based()`. There MUST NOT be implicit `From`/`Into` conversions between systems.

## Accessors

r[pos.get]
`pos.get() -> u32` MUST return the raw value in the position's native coordinate system.

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
`Pos<S> + Offset = Pos<S>`: adding an offset to a position MUST produce a position in the same system. The operation MUST `debug_assert!` that the result is non-negative and fits in `u32`.

r[pos.sub_offset]
`Pos<S> - Offset = Pos<S>`: subtracting an offset from a position MUST produce a position in the same system. Same bounds checking as `add_offset`.

r[pos.no_add_pos]
Adding two positions (`Pos<S> + Pos<S>`) MUST be a compile error. Adding positions is not a meaningful operation.

## Traits

r[pos.derives]
`Pos<S>` MUST derive or implement: `Copy`, `Clone`, `PartialEq`, `Eq`, `PartialOrd`, `Ord`, `Hash`, `Debug`, `Display`.

r[pos.must_use]
All accessor methods and conversion methods MUST be `#[must_use]`.

## Niche optimization

r[pos.niche]
`Pos<S>` SHOULD use a niche-optimized inner type (e.g., via the `nonmax` crate's `NonMaxU32`) so that `Option<Pos<S>>` has the same size as `Pos<S>` (4 bytes). The value `u32::MAX` is reserved as the niche and MUST NOT be a valid position. This is acceptable because the SAM/BAM spec caps reference sequence length at `2^31-1` (~2.1B) [SAM1 §1.3, §4.2], well below `u32::MAX` (~4.3B).
