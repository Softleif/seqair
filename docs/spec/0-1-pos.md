# Genomic Position Types

## Motivation

Genomic positions appear throughout the codebase in multiple integer types (`i64`, `i32`, `u64`, `u32`, `NonZeroU32`) and two coordinate systems (0-based for BAM/BED/internal, 1-based for SAM/VCF/CRAM/user-facing). Manual `- 1` / `+ 1` conversions between systems are scattered across ~15 sites, creating two classes of bugs:

1. **Off-by-one errors**: A bare `pos - 1` doesn't communicate whether the subtraction is coordinate conversion or interval arithmetic. A future refactor could accidentally remove or duplicate it.
2. **Type confusion**: Position values flow through `u64` (public API) → `i64` (records) → `i32` (CRAM headers, CompactOp) with unchecked casts at each boundary. Errors like `.max(1) as u64 - 1` are fragile and opaque.

Two newtypes, `Pos0` and `Pos1`, solve both: the type system distinguishes coordinate systems at compile time, and a single internal integer width (32 bits) eliminates gratuitous casts.

### Why `u32` capped at `i32::MAX`?

The SAM spec caps reference sequence length (`@SQ LN`) at `[1, 2^31-1]` (~2.1 billion) [SAM1 §1.3], and the BAM binary `l_ref` field is `uint32_t` with the same `< 2^31` constraint [SAM1 §4.2]. BAM binary `pos` is `int32_t`, giving a valid position range of `[0, 2^31-1]`. All current SAM/BAM positions therefore fit in `u32` (max ~4.3 billion) with room to spare, and capping at `i32::MAX` makes `as_i32()` infallible. Using `u32`:

- Halves position storage (4 bytes vs 8), improving cache density in hot structs (`SlimRecord`, `PileupAlignment`, `CompactOp`)
- Eliminates the i64↔i32 casts currently scattered through cigar.rs
- `Offset` (distance between positions) uses `i64` for safe intermediate arithmetic

> **Note on large genomes and future widening:** The CSI index format [CSI] uses `int64_t` for positions to support sequences larger than `2^29` bases. Enormous total genome sizes do exist — the fern _Tmesipteris oblanceolata_ has a ~160 Gbp genome ([Hidalgo et al., iScience 2024](<https://www.cell.com/iscience/fulltext/S2589-0042(24)01111-8>)) — but that size is spread across many chromosomes. However, individual chromosomes exceeding `2^31-1` already exist in practice: users have reported single chromosomes of ~2.36 Gbp that BAM simply cannot represent (`[E::bam_write1] Positional data is too large for BAM format`). This is an active unresolved issue in the SAM/BAM spec ([samtools/hts-specs#655](https://github.com/samtools/hts-specs/issues/655)); htslib already uses 64-bit coordinates internally for SAM, and CRAM 4.0 (still in draft) is planned to support them, but standard BAM and CRAM remain limited to `int32_t`. The proposed `PN`/`PO` header tags for stitching split chromosomes were never standardized. For seqair, 32 bits are correct today because BAM enforces the `2^31-1` cap at the format level — but if support for SAM-only large-chromosome files or future CRAM 4.0 is ever added, `Pos0`/`Pos1` would need to widen to 64.

> **Sources:** Coordinate conventions follow [SAM1] §1.4 (1-based POS), §4.2 (0-based BAM binary POS). Reference length limit from [SAM1] §1.3 (`@SQ LN` range `[1, 2^31-1]`) and §4.2 (`l_ref: uint32_t < 2^31`). CSI position width from [CSI]. See [References](./99-references.md).

## Position type

r[pos.type]
`Pos0` and `Pos1` MUST be distinct `#[repr(transparent)]` newtypes over
`NonZeroU32`. `Pos0` MUST store the position **plus one**; `Pos1` MUST store
the value itself. A `Pos0` and the `Pos1` of the same genomic base therefore
have the same representation, which is what makes `to_one_based` a
reinterpretation plus a bound check and `to_zero_based` a bare
reinterpretation. Adding one is monotone, so the derived ordering on the
stored number is the ordering of the positions.

r[pos.size]
`size_of::<Pos0>()` and `size_of::<Pos1>()` MUST each equal `size_of::<u32>()`
(4 bytes).

r[pos.systems]
Two position types MUST be defined:

- `Pos0`: 0-based coordinates (BAM binary, BED, internal engine)
- `Pos1`: 1-based coordinates (SAM text, VCF, CRAM, user-facing)

r[pos.incompatible]
`Pos0` and `Pos1` MUST be distinct types. The compiler MUST reject comparison, assignment, or arithmetic between positions of different systems.

## Construction

r[pos.two_types]
Every bound-aware operation — `new`, `MIN`/`MAX`, the `TryFrom` impls, the
checked and saturating offset arithmetic, `Debug`, and the serde impls — MUST
exist on both types under the same name with the same signature, and each MUST
use its own type's floor (`0` for `Pos0`, `1` for `Pos1`). Neither may inherit
the other's. (That is how `checked_add_offset` shipped accepting
`Pos1(1) + (-1)`: its check was "negative", which is the 0-based floor, and
every test of it used `Pos0`.) The two implementations are written out in full
rather than generated from a shared trait or macro, so the pairing is a
property the tests establish, not one the reader has to take on faith: the
unit tests MUST state each of these properties once and run it against both
types.

r[pos.min_max]
`Pos0::MIN` MUST be `Pos0(0)` and `Pos1::MIN` MUST be `Pos1(1)`; `Pos0::MAX`
and `Pos1::MAX` MUST both be the last representable position, `i32::MAX`.
`Pos0::ZERO` is `Pos0::MIN` under the name it is usually wanted by.

r[pos.zero_new]
`Pos0::new(u32) -> Option<Self>` MUST return `None` if the value exceeds `i32::MAX`. All values in `0..=i32::MAX` are valid 0-based positions.

r[pos.one_new]
`Pos1::new(u32) -> Option<Self>` MUST return `None` if the value is 0 or exceeds `i32::MAX`. The value 0 is not a valid 1-based position — and it is not a runtime check that could be forgotten: a `Pos1` stores a `NonZeroU32` holding the value itself, so the representation cannot hold a 1-based zero.

Both MUST be `const fn`, since `MIN`, `MAX` and `Pos0::ZERO` are built from them.

r[pos.try_from]
Fallible `TryFrom` trait impls MUST be provided on each type for common integer types, all returning `Result<Self, PosOverflow>`:

- `TryFrom<i32>`: rejects negative values (for `Pos0`) or values < 1 (for `Pos1`).
- `TryFrom<u32>`: delegates to `new()`.
- `TryFrom<i64>`: rejects negative values (for `Pos0`) or values < 1 (for `Pos1`), and values exceeding `i32::MAX`.
- `TryFrom<u64>`: rejects values exceeding `i32::MAX` (and 0 for `Pos1`).

Each is `new` behind a widening or narrowing that cannot itself invent a
position. These and `new()` are the ONLY entry points for external integer
values. Raw `as` casts from integers to a position MUST NOT be possible outside
the module.

## Conversion

r[pos.to_one_based]
`Pos0::to_one_based()` MUST return `Result<Pos1, PosOverflow>` naming the base
one further along in 1-based counting. This is fallible: a 0-based position of
`i32::MAX` would require a 1-based value of `i32::MAX + 1`, which exceeds the
valid range for `Pos1`. Because the two types share a representation
(r[`pos.type`]) it MUST be implemented as a reinterpretation of the stored
number plus that bound check, not as arithmetic on the logical value.

r[pos.to_zero_based]
`Pos1::to_zero_based()` MUST return the `Pos0` of the same base. This is
infallible because the minimum 1-based value is 1, yielding 0-based 0, and it
MUST likewise be a bare reinterpretation of the stored number.

r[pos.explicit_conversion]
Coordinate system conversion MUST only occur through `to_one_based()` and
`to_zero_based()`. There MUST NOT be `From`/`Into` conversions between `Pos0`
and `Pos1` in either direction, even though the two share a representation —
the reinterpretation is an implementation detail of those two named methods.

## Accessors

r[pos.as_i32]
`pos.as_i32() -> i32` MUST return the raw value as `i32`.

r[pos.as_usize]
`pos.as_usize() -> usize` MUST return the raw value as `usize` for array indexing.

r[pos.as_i64]
`pos.as_i64() -> i64` MUST return the raw value as `i64` for wider arithmetic.

r[pos.as_u32]
`pos.as_u32() -> u32` MUST return the logical value — for `Pos0`, the stored
number minus one. It is infallible, and so are `as_i32`, `as_i64`, `as_u64` and
`as_usize`, which MUST all report the same logical value. No accessor may ever
expose the stored `NonZeroU32`.

r[pos.no_deref]
`Pos0` and `Pos1` MUST NOT implement `Deref<Target = u32>`. A `Deref` to the
representation defeats the newtype: `*pos` yields a bare integer that flows into
any `u32` parameter, and autoderef silently exposes every `u32` method
(`pos.saturating_sub(n)` resolves and compiles). Callers wanting the number MUST
say so through a named accessor. This mirrors the same prohibition on `Aux` in
[BAM records](./2-bam-3-2-record.md).

## Arithmetic

r[pos.offset]
`Offset` MUST be a `#[repr(transparent)]` newtype over `i64` representing a signed distance between positions.

r[pos.sub_pos]
`Pos0 - Pos0 = Offset` and `Pos1 - Pos1 = Offset`: subtracting two positions of the same system MUST produce an `Offset` holding the difference of the logical values. Subtracting positions of different systems MUST be a compile error.

r[pos.add_offset]
`checked_add_offset(Offset) -> Option<Self>` MUST exist on both types: adding
an offset to a position produces a position in the same system, and MUST return
`None` if the result is below that type's `MIN` or exceeds `i32::MAX`. The
floor is the type's own, not zero: `Pos1(1) + (-1)` is `None`, never a 1-based
zero.

r[pos.sub_offset]
`checked_sub_offset(Offset) -> Option<Self>` MUST exist on both types:
subtracting an offset from a position produces a position in the same system.
Same bounds checking as `checked_add_offset`; `Offset::new(i64::MIN)` has no
negation, so subtracting it is `None` for every position.

r[pos.no_add_pos]
Adding two positions (`Pos0 + Pos0`, `Pos1 + Pos1`) MUST be a compile error. Adding positions is not a meaningful operation.

r[pos.saturating_offset]
`saturating_add_offset(Offset) -> Self` and `saturating_sub_offset(Offset) -> Self`
MUST exist on both types and MUST clamp instead of failing: a result below
`MIN` becomes `MIN`, and a result above `i32::MAX` becomes `MAX`. Wherever the
checked form succeeds, the saturating form MUST return the same position; it
may only differ where the checked form is `None`, and there it MUST land on the
bound that was crossed.

These exist because padding a window and trimming an overlap are saturating by
nature — two bases before the start of a contig is the start of the contig, not
an error — and expressing that with `checked_*` forces every such call site to
invent the same `unwrap_or` fallback. A caller that needs to *know* the bound
was hit MUST use the checked form; these MUST NOT be used to paper over a
position that should have been validated.

## Serialization

r[pos.serde]
Under the `serde` feature, `Pos0` and `Pos1` MUST each implement `Serialize`
and `Deserialize`. Both MUST use the identical wire format: the bare logical
value as a `u32`, with no wrapper and no coordinate-system tag — which system a
field is in is the field's type, not something the document restates. In
particular `Pos0` MUST serialize its logical value, never the stored
`value + 1`. Deserialization MUST re-check the bound on the way in and reject
exactly what that type's `new` rejects — above `i32::MAX` for either, `0` for
`Pos1` — so a round trip cannot introduce a position the type says is
impossible, and the error MUST name the offending value and the range it was
expected in.

r[pos.serde_optional]
The serde impls MUST live behind a `serde` feature, so a caller that only
reads BAM can drop the dependency with `default-features = false`. The feature
MUST be **on by default**: the impls on every other value type existed
unconditionally through 0.2, and removing them from a default build would
break callers that never asked for a feature at all.

## Reference intervals

Every reader and the pileup engine work on the same kind of interval, and the
one number worth being pedantic about is where it stops. Half-open ends and
inclusive ends both appear in the file formats this library reads (htslib's
`bam_endpos` is exclusive, CRAM's slice spans are inclusive), so seqair picks
one internally and converts at the format boundary.

r[interval.inclusive_ends]
A reference interval in seqair's API is 0-based with **both ends inclusive**:
`[start, last]` covers `last - start + 1` positions. This applies to the span
passed to `fetch_into`/`fetch_into_customized`, to `Segment`, to the pileup
engine's region, and to the FASTA fetch. A half-open interval MUST NOT be used
at these boundaries; where an external format needs one (htslib's `bam_endpos`,
a BED line), the conversion happens in the format's own reader or writer.

r[interval.span_type]
An interval that crosses an API boundary MUST be carried by one value of
`core::range::RangeInclusive<Pos0>`, not by two loose positions. Every query
(`fetch_into`, `fetch_into_customized`, `estimate_region_bytes`,
`IndexedBamReader::query`, `PileupEngine::new`, the `(resolver, span)` segment
target) and every reference fetch (`fetch_seq`, `fetch_seq_into`,
`fetch_base_seq`) takes one, and `Segment::span()` / `core_span()` produce one,
so a tile is handed to a fetch whole and the reference for a region is fetched
with the very value that queried it. There MUST NOT be a half-open
`Range<Pos0>` anywhere in the public API, and a loose-position form MUST NOT be
offered alongside the span form: as long as `fetch_seq(name, start, stop)`
exists it is what gets called, and the type carries nothing.

One convention buys two things. There is no `± 1` at the boundary for a caller
to get wrong — the historical bug class here is a half-open reference buffer
handed to code that treats the segment's `end` as inclusive, one base short and
silent. And a closed span can name the last representable position, where the
half-open `end = i32::MAX + 1` is not a `Pos0` and every reader had grown its
own `u64` side door to reach it.

`core::range::RangeInclusive` is preferred over `std::ops::RangeInclusive`: it
is `Copy`, 8 bytes rather than 12 (no exhaustion flag), does not implement
`Iterator` so a position span cannot be walked by accident, and names its
inclusive end `last`, which is what it is. `(a..=b).into()` converts from the
`..=` literal.

r[interval.empty_span]
A span whose `last` is before its `start` is empty, as it is for Rust's own
ranges, and every query over it MUST return no records and no error. The
overlap test alone does not give this: `record.pos <= last && record.end_pos
>= start` holds for a record that covers the whole gap between the two ends,
which is how a reversed span used to fetch exactly the reads spanning it.

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

## Query positions

r[qpos.type]
`QPos` MUST be a `#[repr(transparent)]` newtype over `u32` representing a
0-based offset into a read's sequence.

r[qpos.not_a_pos]
`QPos`, `Pos0` and `Pos1` MUST be distinct, non-interconvertible types. A
`Pos0`/`Pos1` names a location on the reference; a `QPos` indexes a read. The
compiler MUST reject passing one where the other is expected.

> **Why this exists.** Code that resolves reference bases from a segment computes
> `pos - segment_start`. Handed a query offset instead of a genomic position that
> subtraction underflows — and where it was written with `saturating_sub`, it
> clamps to 0 and silently returns bases from the start of the segment: real
> sequence from the wrong locus, no error. This was a live bug in rastair's
> deletion alleles, where 82.6 % of multi-base REF alleles did not match the
> reference. Both integers were `usize`, so nothing could catch it.
>
> **Boundary of the guarantee.** The newtype stops mix-ups only where a
> signature names `QPos` or `Pos`; the raw accessors are deliberate escape
> hatches (slice indexing requires `usize`) and `Pos0::try_from(qpos.get())`
> still compiles. Consumers MUST type their own position-carrying parameters
> with `QPos`/`Pos0` for the guarantee to extend past the seqair API boundary.
> Neither `QPos` nor `Pos0`/`Pos1` implements `From`/`Into` for integer types —
> extraction is always the explicit, named accessor.

r[qpos.new]
`QPos::new(u32) -> Self` MUST be infallible. Unlike a reference position, no
format caps a query offset below `u32::MAX`; it is bounded by the read length.

r[qpos.get]
`QPos::get() -> u32` and `QPos::as_usize() -> usize` MUST return the raw offset.

r[qpos.checked]
`QPos::checked_add(u32)` and `QPos::checked_sub(u32)` MUST return `Option<Self>`.

r[qpos.saturating]
`QPos::saturating_add(u32) -> Self` MUST saturate rather than wrap, so a
malformed CIGAR cannot walk a query cursor back to the start of the read.

## Traits

r[pos.derives]
`Pos0` and `Pos1` MUST each derive or implement: `Copy`, `Clone`, `PartialEq`,
`Eq`, `PartialOrd`, `Ord`, `Hash`, `Debug`, `Display`. `Ord` MUST order by the
logical value, which the derive on the stored `NonZeroU32` already gives
(r[`pos.type`]). `Debug` MUST print the logical value under the type's own
name — `Pos0(7)`, `Pos1(7)` — never the stored one, and `Display` the bare
number.

r[pos.must_use]
All accessor methods and conversion methods MUST be `#[must_use]`, on both types.

## Niche optimization

r[pos.niche]
`Option<Pos0>` and `Option<Pos1>` MUST be 4 bytes, the same as the position
itself. Both types store a `NonZeroU32` (r[`pos.type`]), whose niche the
`Option` discriminant occupies: `Pos1` holds the value, which is never zero
anyway, and `Pos0` holds the value plus one. That is affordable because the
SAM/BAM spec caps reference sequence length at `2^31-1` (~2.1B) [SAM1 §1.3,
§4.2], so `value + 1` cannot overflow `u32` (~4.3B). This matters wherever an
absent position is a normal state — `FilterRawFields::end_pos`,
`OwnedBamRecord::pos` and `next_pos` are all `Option<Pos0>`.
