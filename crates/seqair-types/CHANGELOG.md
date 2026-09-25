# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Changed

- New dependencies `fearless_simd` and `fearless_simd_macros`.

### Performance

- The ASCII→`Base` kernel behind `Base::from_ascii_vec`, `convert_ascii_in_place` and
  `convert_ascii_in_place_as_slice` is one portable `fearless_simd` kernel instead of AVX2, SSSE3
  and NEON copies with their `unsafe` blocks. It picks the widest level the CPU has at run time, so
  AVX-512 machines get 64-byte vectors, and handles a ragged tail as one overlapping vector: 0.70–0.78×
  the time of 0.3.0 on a Ryzen 3950X.

## v0.3.0 (2026-09-22)

### Breaking

- **`Pos<S>` is gone; `Pos0` and `Pos1` are two structs.** The generic form shared code between the
  systems without sharing the floor (see Fixed). Both are `#[repr(transparent)]` over `NonZeroU32`:
  `Pos1` holds its value, so the 1-based zero is not a check but a value that cannot exist; `Pos0`
  holds `value + 1`, so a `Pos0` and the `Pos1` of the same base share one bit pattern and the two
  conversions are reinterpretations. `Option<Pos0>` and `Option<Pos1>` are now 4 bytes. The marker
  types `Zero` and `One` and the `Pos` name are removed; code that spelled `Pos<Zero>` writes `Pos0`.
- `Pos0::max_value()` / `Pos1::max_value()` → the constants `Pos0::MAX` / `Pos1::MAX`, alongside
  new `MIN` (`Pos0(0)`, `Pos1(1)`).
- serde support now lives behind a `serde` feature, **on by default**. `default-features = false`
  drops serde, `smol_str/serde` and `smallvec/serde` for callers who only parse BAM. One caveat,
  invisible to `cargo semver-checks`: 0.2 had no `default` feature, so `default-features = false`
  written alongside `features = ["hts-compat"]` was a no-op. On 0.3 that same line removes every
  `Serialize`/`Deserialize` impl; such callers need `features = ["hts-compat", "serde"]`.

### Fixed

- **`Pos1::checked_add_offset` could produce a 1-based zero.** `Pos1::new(1).checked_add_offset(-1)`
  returned `Some(Pos1(0))`, the value `Pos1::new` exists to refuse; `checked_sub_offset` likewise.
  The bound check was "negative", which is the 0-based floor, and every test of the arithmetic used
  `Pos0`, where the floor and "negative" coincide. Each type now checks against its own floor, and
  every property test is stated once and run against both types.

### Added

- `saturating_add_offset` / `saturating_sub_offset` on both types: clamped to the type's range (`0`
  for `Pos0`, `1` for `Pos1`) instead of returning `None`. Padding a window and trimming an overlap
  are saturating by nature — two bases before the start of a contig is the start of the contig — and
  the checked forms made every such call site invent the same `unwrap_or`.
- `Serialize`/`Deserialize` for `Pos0` and `Pos1`. The wire format is the bare value as a `u32`;
  deserialization re-checks what the constructor checks, so a document cannot produce a `Pos1` of
  `0`, and the error names the offending value and the range it was expected in.

## v0.2.0 (2026-09-14)

### Breaking

- `strand_from_flags(flags)` removed → `Strand::from(flags)`.
- `Pos<S>` no longer implements `Deref<Target = u32>`. The deref defeated the newtype — `*pos` flowed
  into any `u32` parameter and autoderef exposed every `u32` method. Use `Pos::as_u32()`.
- MSRV raised to 1.98.1 (algebraic float ops in `Probability` → `Phred` and `RmsAccumulator`).

### Added

- `QPos`: a `repr(transparent)` 0-based offset into a read's sequence, deliberately not
  interconvertible with `Pos` — confusing the two silently resolves reference bases from the wrong
  locus. Carried by seqair's `PileupOp`, `AlignedPair` and `CigarPosInfo` families.
- `SumOfSquares`: a bare sum whose count is supplied at `finish(count)`, for callers keeping several
  sums over the same values (per-strand, per-allele) where `RmsAccumulator`'s own count is
  duplication. Bit-for-bit identical to `RmsAccumulator`. Merging is `AddAssign`, not `Add`.
- `RootMeanSquare::from_sum_of_squares(sum, count)` for a sum accumulated elsewhere.
- `Pos::as_u32()`, replacing the removed deref.

## v0.1.0 (2026-05-08)

First ever release. Changelog entries will start in next one.
