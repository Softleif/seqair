# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.2.0 (Unreleased)

Highlights: a streaming-window rewrite of the BAM region reader, a unified filter/pileup API, and VCF encoder ergonomics.

### Breaking

- **`RegionBuf` streams a bounded sliding window** instead of loading whole regions.
  It now borrows the reader, so it gained a lifetime and a required generic type parameter (`RegionBuf<'a, _>`),
  no longer implements `UnwindSafe`, and `RegionBuf::load` was removed (use `ensure_available` / `advance_range`).
- **`Readers::pileup` takes an extra parameter**, the optional reference is now threaded through.
- **`PileupAlignment`** does not expose `strand` anymore. Use `Strand::from(rec.flags)` (or your own logic) instead.
- **Unmapped-read filtering unified on `filter_raw`.**
  Removed `IndexedBamReader::keep_unmapped` / `keeps_unmapped`.
  `FilterRawFields` is now type-level enums: public fields `raw_cigar_bytes`, `cigar_ops`, `packed_seq`, `bases` removed;
  `end_pos` is now `Option<Pos0>`;
  new `RejectUnmapped` customizer.
- **VCF encoder traits `EncodeInfo` and `EncodeFormat` removed**.
  `FormatEncoder` gained three required methods (no defaults): `format_ints`, `format_floats`, `format_string`.
- **Trait `TargetInfoAccess` removed** from `seqair::bam::header`.
- **`ReaderError`** gained a new variant distinguishing "region end past contig" from "exceeds `i32::MAX`.
- **`VcfHeaderError::TooManyFields` changed from a unit variant to data-carrying.**
- **Query offsets are now the `QPos` newtype, not bare `u32`/`usize`.**
  Changed: `AlignedPair::{Match, Insertion, SoftClip}.qpos`,
  `MatchPosition` / `MatchedBase` / `MatchedRef` / `AlignedPairWithRead::*` / `AlignedPairWithRef::*`,
  `CigarPosInfo::{Match, Insertion}.qpos`, `CigarMapping::soft_clip_qpos_at` (now returns `Option<QPos>`),
  `PileupOp::{Match, Insertion, SoftClip}.qpos`, and `PileupAlignment::qpos()` (now returns `Option<QPos>`).
  `QPos` and `Pos` are deliberately non-interconvertible: a read-local offset handed to code
  expecting a genomic position previously compiled and silently resolved bases from the wrong locus.
- **`BaseModState::mod_at_qpos` / `is_unmodified` take `QPos`** instead of `usize`, so a genomic
  position can no longer be passed where a read offset is expected.

### Added

- `engine.pileups_into(&mut buf)`, `AlignmentView::inserted_bases` / `inserted_quals`, and anchor-addressed `PileupAlignment::indel_after`.
- `engine.set_soft_clip_overhang(n)` allows getting `n` soft-clipped bases at the fringes of alignments
- `SlimRecord::end_pos_htslib()` for htslib-compatible unmapped read spans.
- Byte-aware segment planning to bound per-segment memory.
- VCF FORMAT parity: `FormatInts`, `FormatString`, array-of-floats; percent-encode FORMAT string values;
  O(1) duplicate-field detection with htslib-compatible in-place overwrite.
- `Writer::finish` returns a self-serializing `CoordinateIndex`.

### Fixed

- CRAM: missing rANS Nx16 validation and allocation-size validation.
- Region queries no longer error on unplaced reads (`pos = -1` / `AP = 0`).
- Pileup: always set depth cap; stable `retain` replaces `swap_remove` in eviction.
- Empty-sample BCF missing-value handling.

### Performance

- Pileup scratch buffers pooled across regions.
- `CompactOp` shrunk 16 → 12 bytes (len + op type packed into one `u32`).
- `CigarSlice` enum collapsed to `&[CigarOp]` via zero-cost transmute.

## v0.1.0 (2026-05-08)

First ever release. Changelog entries will start in next one.
