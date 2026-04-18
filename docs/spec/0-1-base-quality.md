# Base Quality Type

> **Sources:** SAM specification [SAM1] §1.4 (QUAL field, Phred+33 text encoding), §4.2.3 (BAM binary QUAL layout and 0xFF "unavailable" sentinel). See [References](./99-references.md).

## Motivation

BAM and SAM represent per-base quality as a single byte per base. The `0xFF` value is reserved as a whole-record "quality unavailable" sentinel (SAM text `*`, BAM binary array of all `0xFF`). Treating quality as raw `u8` means every consumer must remember to special-case `0xFF` — a constraint that is easy to forget and silently gives wrong answers (e.g. a threshold check `q >= 20` silently admits unavailable-qual reads as maximum quality).

`BaseQuality` is a `#[repr(transparent)]` wrapper that preserves the wire byte but forces consumers to acknowledge the sentinel at access time. This is distinct from `Phred` (which is an `f64` for error-probability math); `BaseQuality` is the wire-level per-base type.

Parsing is zero-cost: `&[u8]` from any BAM/SAM/CRAM decode path is reinterpreted as `&[BaseQuality]` via a pointer cast, relying on `#[repr(transparent)]` over `u8`. There is no scan, copy, or conditional at parse time.

## Type definition

r[types.base_quality.type]
`BaseQuality` MUST be a `#[repr(transparent)]` newtype over `u8`.

r[types.base_quality.size]
`size_of::<BaseQuality>()` MUST equal `size_of::<u8>()` (1 byte), and `align_of::<BaseQuality>()` MUST equal `align_of::<u8>()`. This is a prerequisite for the zero-copy slice cast.

r[types.base_quality.copy]
`BaseQuality` MUST implement `Copy`, `Clone`, `PartialEq`, `Eq`, `Hash`, and `Debug`.

r[types.base_quality.no_ord]
`BaseQuality` MUST NOT implement `Ord` or `PartialOrd`. The `0xFF` sentinel has no meaningful ordering against real Phred values, and a derived `Ord` would silently place "unavailable" above all valid qualities. Consumers that need ordering MUST extract the Phred via `get()` and decide explicitly.

## Unavailable sentinel

r[types.base_quality.unavailable]
`BaseQuality::UNAVAILABLE` MUST be the constant `BaseQuality(0xFF)`, matching the SAM/BAM wire sentinel for "quality unavailable" per [SAM1] §4.2.3.

## Construction

r[types.base_quality.from_byte]
`BaseQuality::from_byte(b: u8) -> Self` MUST construct a value from a raw wire byte without validation. All 256 byte values are valid (including `0xFF` = unavailable).

## Access

r[types.base_quality.get]
`BaseQuality::get(self) -> Option<u8>` MUST return `None` when the wire byte is `0xFF` (unavailable) and `Some(q)` for any other byte. This is the primary accessor and forces consumers to handle the absent case.

r[types.base_quality.as_byte]
`BaseQuality::as_byte(self) -> u8` MUST return the raw wire byte unchanged, including the `0xFF` sentinel. This is an escape hatch for serialization and interop at wire boundaries; consumers that want the Phred integer MUST use `get()`.

## Zero-copy slice cast

r[types.base_quality.slice_from_bytes]
`BaseQuality::slice_from_bytes(bytes: &[u8]) -> &[BaseQuality]` MUST return a view of the input bytes with identical length and lifetime. The implementation is a pointer cast sound by `r[types.base_quality.type]` and `r[types.base_quality.size]`. No allocation, no copy.

r[types.base_quality.slice_to_bytes]
`BaseQuality::slice_to_bytes(quals: &[BaseQuality]) -> &[u8]` MUST return the inverse view, preserving length and lifetime.

## Field usage

r[types.base_quality.field_type]
Public qual APIs — `RecordStore::qual`, `BamRecord::qual`, `OwnedBamRecord::qual`/`set_qual`, and `PileupOp::Match { qual }` / `PileupOp::Insertion { qual }` — MUST expose quality as `&[BaseQuality]` (or `BaseQuality` for per-base positions). Internal slab storage and reader scratch buffers MAY remain `Vec<u8>` / `&mut Vec<u8>`; the cast happens at the public API boundary.

## Round-trip

r[types.base_quality.roundtrip]
For all `v: u8`, `BaseQuality::from_byte(v).as_byte() == v` MUST hold. For all `&[u8]` slices, `slice_to_bytes(slice_from_bytes(s))` MUST produce a slice pointing to the same bytes with the same length.
