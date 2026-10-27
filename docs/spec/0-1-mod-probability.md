# Modification Probability Type

> **Sources:** SAM specification [SAM1] §1.7 (base modifications; the `ML:B:C` tag: "the byte x encodes the probability range [x/256, (x+1)/256)"). See [References](./99-references.md).

## Motivation

The `ML` tag stores one byte per modification call. As a raw `u8` it is easily read as a probability on a 0–255 scale (`x/255`) or as a quality; it is the index of a 1/256-wide bucket. `ModProbability` is a `#[repr(transparent)]` wrapper that names the byte and converts to a probability in one place, at the bucket midpoint. It mirrors `BaseQuality` and `MappingQuality`, but has no sentinel value.

## Type definition

r[types.mod_probability.type]
`ModProbability` MUST be a `#[repr(transparent)]` newtype over `u8`, with `size_of` and `align_of` equal to `u8`'s.

r[types.mod_probability.copy]
`ModProbability` MUST implement `Copy`, `Clone`, `PartialEq`, `Eq`, `Hash`, and `Debug`.

r[types.mod_probability.ord]
`ModProbability` MUST implement `Ord` and `PartialOrd` by the byte. All 256 values are real buckets, so a larger byte is a higher probability and a byte comparison is a probability comparison.

## Construction and access

r[types.mod_probability.from_byte]
`ModProbability::from_byte(b: u8) -> Self` MUST construct a value from a raw `ML` byte without validation. All 256 byte values are valid.

r[types.mod_probability.as_byte]
`ModProbability::as_byte(self) -> u8` MUST return the raw `ML` byte unchanged.

r[types.mod_probability.roundtrip]
For all `v: u8`, `ModProbability::from_byte(v).as_byte() == v` MUST hold.

r[types.mod_probability.midpoint]
`ModProbability::midpoint(self) -> Probability` MUST return `(x + 0.5) / 256`, the centre of the range `[x/256, (x+1)/256)` of [SAM1] §1.7. The result is strictly between 0 and 1.

## Field usage

r[types.mod_probability.field_type]
`Modification::probability` MUST be a `ModProbability`; the `ML` decoder wraps each byte unchanged.
