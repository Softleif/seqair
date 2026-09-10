# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.2.0 (Unreleased)

### Breaking

- Removed `strand_from_flags(bam_flags)`. This is now `Strand::from(bam_flags)`.
- `Pos<S>` no longer implements `Deref<Target = u32>`. The deref defeated the newtype: `*pos` yielded a
  bare integer that flowed into any `u32` parameter, and autoderef silently exposed every `u32` method
  (`pos.saturating_sub(n)` resolved and compiled). Use the new named accessors — `Pos::as_u32()` — to say
  so explicitly.

### Added

- `QPos`: a `repr(transparent)` newtype for a 0-based offset into a read's sequence (query coordinates).
  Deliberately not interconvertible with `Pos` — a `Pos` names a place on the reference, a `QPos` indexes a
  read, and confusing the two silently resolves reference bases from the wrong locus. Construction is
  infallible (no format caps a query offset; it is bounded by the read length). Carried by seqair's
  `PileupOp`, `AlignedPair`, and `CigarPosInfo` families.
- `Pos::as_u32()` accessor, replacing the removed deref.
- `SumOfSquares`: a bare sum of squares whose count is supplied at `finish(count)`, for callers keeping
  several sums over the same values — per-strand and per-allele tallies, say — where `RmsAccumulator`'s
  own `count` is duplication. rastair holds nine such sums per allele and updates them once per read per
  column; eight of the nine counter increments were copies of a number it already had, measuring 5.2 % of
  its worker CPU. It squares through the same algebraic ops and routes to the same constructor
  `RmsAccumulator::finish` uses, so the two agree bit for bit (pinned by a test on the bits, not a
  tolerance — these values are printed into files people diff). Merging is `AddAssign` and deliberately
  not `Add`, so `add(value)` cannot be confused with adding another sum.
- `RootMeanSquare::from_sum_of_squares(sum, count)`, for a sum accumulated elsewhere. Previously the only
  route was a one-element `RmsAccumulator` fed the mean square — exact, but a riddle.
- `Probability` → `Phred` conversions and `RmsAccumulator` arithmetic use algebraic float ops
  (`algebraic_mul` / `algebraic_add` / `algebraic_div`); MSRV raised to 1.98.1 for them.

## v0.1.0 (2026-05-08)

First ever release. Changelog entries will start in next one.
