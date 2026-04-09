# BAM Flags Type

> **Sources:** SAM specification [SAM1] §1.4 (FLAG field), §4.2 (BAM binary FLAG). See [References](./99-references.md).

## Motivation

BAM flags are a 16-bit bitfield with 12 defined bits. Raw `u16` fields invite magic-number comparisons (`flags & 0x4 != 0`) scattered across the codebase, which are hard to read and easy to get wrong. A strongly-typed newtype prevents accidental mixing with other `u16` values and provides self-documenting predicate methods.

## Type definition

r[flags.type]
`BamFlags` MUST be a `#[repr(transparent)]` newtype over `u16`.

r[flags.size]
`size_of::<BamFlags>()` MUST equal `size_of::<u16>()` (2 bytes).

r[flags.copy]
`BamFlags` MUST implement `Copy`, `Clone`, `PartialEq`, `Eq`, `Hash`, and `Debug`.

## Construction

r[flags.from]
`BamFlags` MUST implement `From<u16>` and `Into<u16>` (via `From<BamFlags> for u16`). These are the canonical conversions — all 16 bits are valid storage even if only 12 are currently assigned by the SAM spec.

r[flags.raw]
`BamFlags::raw(self) -> u16` MUST return the underlying value unchanged (convenience alias for `u16::from`).

r[flags.roundtrip]
For all `v: u16`, `BamFlags::from(v).raw() == v` MUST hold.

r[flags.empty]
`BamFlags::empty()` MUST return `BamFlags(0)` (no flags set).

## Named constants

> r[flags.constants]
> Named constants MUST be provided for all 12 SAM-defined flag bits:
>
> | Constant                  | Value    | SAM bit | Meaning                      |
> | ------------------------- | -------- | ------- | ---------------------------- |
> | `FLAG_PAIRED`             | `0x0001` | 0x1     | Read is paired in sequencing |
> | `FLAG_PROPER_PAIR`        | `0x0002` | 0x2     | Read mapped in a proper pair |
> | `FLAG_UNMAPPED`           | `0x0004` | 0x4     | Read is unmapped             |
> | `FLAG_MATE_UNMAPPED`      | `0x0008` | 0x8     | Mate is unmapped             |
> | `FLAG_REVERSE`            | `0x0010` | 0x10    | Read on reverse strand       |
> | `FLAG_MATE_REVERSE`       | `0x0020` | 0x20    | Mate on reverse strand       |
> | `FLAG_FIRST_IN_TEMPLATE`  | `0x0040` | 0x40    | First in template            |
> | `FLAG_SECOND_IN_TEMPLATE` | `0x0080` | 0x80    | Last in template             |
> | `FLAG_SECONDARY`          | `0x0100` | 0x100   | Secondary alignment          |
> | `FLAG_FAILED_QC`          | `0x0200` | 0x200   | Failed quality checks        |
> | `FLAG_DUPLICATE`          | `0x0400` | 0x400   | PCR or optical duplicate     |
> | `FLAG_SUPPLEMENTARY`      | `0x0800` | 0x800   | Supplementary alignment      |

## Predicates

r[flags.predicates]
`BamFlags` MUST provide `is_*` predicate methods for every named constant. Each predicate MUST return `self.0 & CONSTANT != 0`.

r[flags.is_set]
`BamFlags::is_set(self, flag: u16) -> bool` MUST test an arbitrary flag bit, enabling downstream code to check vendor-defined or future bits without needing a named method.

## Mutators

Mutator and predicate methods accept raw `u16` constants (not `BamFlags`) as the flag argument. This is intentional: named constants like `FLAG_UNMAPPED` are `u16`, and the SAM spec may add new bits or vendors may use currently-undefined bits. Accepting `u16` keeps the API ergonomic (`flags.is_set(FLAG_UNMAPPED)`) without requiring a wrapping step for every constant.

r[flags.set]
`BamFlags::set(&mut self, flag: u16)` MUST set the given bit(s): `self.0 |= flag`.

r[flags.unset]
`BamFlags::unset(&mut self, flag: u16)` MUST clear the given bit(s): `self.0 &= !flag`.

r[flags.with]
`BamFlags::with(self, flag: u16) -> Self` MUST return a copy with the given bit(s) set (builder-style, non-mutating).

r[flags.without]
`BamFlags::without(self, flag: u16) -> Self` MUST return a copy with the given bit(s) cleared.

## Bitwise operations

r[flags.bitor]
`BamFlags` MUST implement `BitOr<BamFlags>` returning `BamFlags`, so flags can be combined: `BamFlags::from(FLAG_PAIRED) | BamFlags::from(FLAG_REVERSE)`.

r[flags.bitor_assign]
`BamFlags` MUST implement `BitOrAssign<BamFlags>` (`|=`).

r[flags.bitand]
`BamFlags` MUST implement `BitAnd<BamFlags>` returning `BamFlags`, so flags can be intersected.

r[flags.bitand_assign]
`BamFlags` MUST implement `BitAndAssign<BamFlags>` (`&=`).

## Display

r[flags.display]
`Display` for `BamFlags` MUST output the decimal `u16` value (matching SAM text FLAG field format).

r[flags.format_traits]
`BamFlags` MUST implement `LowerHex`, `UpperHex`, `Octal`, and `Binary`, all delegating to the underlying `u16`. This enables standard format strings (`{:x}`, `{:#06x}`, `{:b}`, etc.) to work identically to raw `u16`.

## Field usage

r[flags.field_type]
All BAM flag fields in the codebase (`BamRecord`, `SlimRecord`, `OwnedBamRecord`, `PileupAlignment`, `ActiveRecord`, `ParsedHeader`, etc.) MUST use `BamFlags` instead of raw `u16`. Exception: test helpers that construct raw BAM wire-format bytes (e.g. `make_test_record`) and comparison structs that mirror external library types (`HtsRecord`, `SyntheticRead`) MAY use `u16` at the serialization/interop boundary.

r[flags.filter_signature]
The pileup `RecordFilter` closure MUST accept `BamFlags` (not `u16`) as its first argument.
