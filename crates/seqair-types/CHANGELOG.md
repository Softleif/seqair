# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
