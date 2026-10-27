# Mapping Quality Type

> **Sources:** SAM specification [SAM1] §1.4 (MAPQ field: `−10 log10 Pr{mapping position is wrong}`, rounded; "A value 255 indicates that the mapping quality is not available"), §4.2 (BAM binary `mapq` byte). See [References](./99-references.md).

## Motivation

SAM, BAM and CRAM store a record's mapping quality as one byte, and `255` means "not available". Kept as a raw `u8`, the sentinel reads as the most confident mapping there is: a mean over reads, a weight, or a Phred-to-probability conversion silently treats an aligner that writes no MAPQ as if it were certain of every placement.

`MappingQuality` is a `#[repr(transparent)]` wrapper over the wire byte that makes consumers handle the sentinel when they read the number. It mirrors `BaseQuality` (`r[types.base_quality.type]`).

## Type definition

r[types.mapping_quality.type]
`MappingQuality` MUST be a `#[repr(transparent)]` newtype over `u8`, with `size_of` and `align_of` equal to `u8`'s.

r[types.mapping_quality.copy]
`MappingQuality` MUST implement `Copy`, `Clone`, `PartialEq`, `Eq`, `Hash`, and `Debug`.

r[types.mapping_quality.no_ord]
`MappingQuality` MUST NOT implement `Ord` or `PartialOrd`. A derived order would place "unavailable" above every real MAPQ, so a minimum-MAPQ filter written as a comparison would decide what happens to unavailable reads without saying so. A filter MUST extract the value via `get()` and handle `None` explicitly; the SAM convention (samtools `view -q`, `mpileup -q`) is that 255 passes a threshold.

## Unavailable sentinel

r[types.mapping_quality.unavailable]
`MappingQuality::UNAVAILABLE` MUST be the constant `MappingQuality(255)`, the "not available" value of [SAM1] §1.4.

## Construction and access

r[types.mapping_quality.from_byte]
`MappingQuality::from_byte(b: u8) -> Self` MUST construct a value from a raw wire byte without validation. All 256 byte values are valid.

r[types.mapping_quality.get]
`MappingQuality::get(self) -> Option<u8>` MUST return `None` for `255` and `Some(q)` for any other byte. This is the primary accessor.

r[types.mapping_quality.as_byte]
`MappingQuality::as_byte(self) -> u8` MUST return the raw wire byte unchanged, including `255`. It is for serialization and interop; consumers that want the number MUST use `get()`.

r[types.mapping_quality.roundtrip]
For all `v: u8`, `MappingQuality::from_byte(v).as_byte() == v` MUST hold.

## Field usage

r[types.mapping_quality.field_type]
Every public MAPQ field and accessor — `SlimRecord::mapq`, `FilterRawFields::mapq`, `PileupAlignment::mapq`, `OwnedBamRecord::mapq` and its builder, and the `mapq` argument of `RecordStore::push_fields` — MUST be a `MappingQuality`. The SAM, BAM and CRAM decoders wrap the decoded byte unchanged, so a 255 in any input format reaches the consumer as `UNAVAILABLE`; the BAM encoder writes `as_byte()`.
