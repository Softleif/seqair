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

## v0.1.0 (2026-05-08)

First ever release. Changelog entries will start in next one.
