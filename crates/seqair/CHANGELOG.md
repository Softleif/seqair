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
- **`Readers::pileup` now takes a `DepthLimit`** (it previously took just the `Segment`).
  Supplying reference bases or rewriting the store first moved to the new `pileup_with_reference` /
  `pileup_with` methods instead.
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
- `Readers::pileup_with(mutator)` runs a caller-supplied mutation on the freshly fetched `RecordStore`
  (then re-sorts by position) — the hook for in-place local realignment via `RecordStore::set_alignment`;
  `Readers::pileup_with_reference(&RefSeq)` drives the engine from a reference the caller already holds
  instead of re-reading the FASTA per segment, rejecting one that doesn't cover the segment
  (`ReaderError::SuppliedReferenceTooSmall`).
- Record-store mate linking for overlap-dedup consumers: `RecordStore::link_mates()` pairs a template's
  primary alignments once per store (qname-hash table, verified by qname bytes and reciprocal mate
  positions; secondary/supplementary and nameless reads never link). The pileup exposes
  `PileupAlignment::mate_idx()` / `in_mate_overlap()` / `qname_hash()`, `PileupColumn::find_record()`,
  and `PileupColumn::pair_indel()` — which recovers a fragment's indel from its linked mate in the same
  column when the kept read carries none (own indel always wins). The cached mate overlap is widened by
  the soft-clip overhang, so rescued soft-clipped fringe bases still count as overlapping.
- Indexed FASTQ references: the FAI parser accepts the sixth `qual_offset` column and `FaiEntry` exposes
  `is_fastq()` / `qual_byte_offset()`. Fetch is byte-correct on wrapped records whose quality lines may
  begin with `@`/`+`. htslib-validated.

### Fixed

- CRAM: missing rANS Nx16 validation and allocation-size validation.
- Region queries no longer error on unplaced reads (`pos = -1` / `AP = 0`).
- Pileup: always set depth cap; stable `retain` replaces `swap_remove` in eviction.
- Empty-sample BCF missing-value handling.
- SAM and CRAM region queries could silently lose records overlapping the query by its last base
  (SAM tested half-open `[start, end)` against an exclusive `end_pos`; CRAI compared a 1-based
  `alignment_start` to a 0-based query). Readers now share one interval convention — 0-based, both
  ends inclusive — and `DepthCap` counts a read in the same columns across BAM, SAM, and CRAM.
- FASTA fetches are bounded by data that can exist: an untrusted `.fai` can no longer drive a
  multi-PiB allocation (line width capped to what a writer can produce) or serve bytes from a
  wrapped file offset (`IndexOffsetOverflow` instead of wrapping); a hand-built `linebases == 0`
  entry errors instead of panicking.
- FASTA errors carry the offending path and distinguish plain gzip (`PlainGzipUnsupported`) from a
  missing `.gzi` (`GziIndexNotFound`); a FASTQ entry whose `qual_offset` precedes its sequence block
  is rejected (`QualBeforeSequence`).

### Performance

- Pileup scratch buffers pooled across regions.
- `CompactOp` shrunk 16 → 12 bytes (len + op type packed into one `u32`).
- `CigarSlice` enum collapsed to `&[CigarOp]` via zero-cost transmute.
- Plain-FASTA fetches use one positional `pread` per span instead of `seek` + `read` (~28% faster on
  the short-slice pileup pattern); plain-FASTA forks share a single handle since `pread` carries no
  cursor (BGZF forks still open their own).
- Fetched spans are stripped and uppercased on a vectorized path (memchr-style newline scan + SIMD
  range-masked ASCII uppercase).

## v0.1.0 (2026-05-08)

First ever release. Changelog entries will start in next one.
