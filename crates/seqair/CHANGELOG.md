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
- **`Readers::pileup(segment, depth)` returns a `Pileup` plan, not a started pileup** — finish it with
  `.run()`. It previously took just the `Segment` and returned the guard directly.
  `pileup_with` and `pileup_with_reference` are gone: they are `.mutate(f)` and `.with_reference(r)` on the
  plan, in any order and in any combination. The combination — a caller that both holds the region's bases
  and rewrites the store against them — had no spelling before, though the private implementation already
  supported it.
- **Record-store indices are the `RecordIdx` newtype, not bare `u32`.**
  Changed: `RecordStore::{record, try_record, qname, cigar, seq, seq_at, qual, aux, extra, extra_mut,
  set_alignment, set_template_len, set_mate_info, mate_overlap}`, the `Option<RecordIdx>` returned by
  `push_raw`/`push_fields`, `SlimRecord::mate_idx()`, `PileupAlignment::{record_idx, mate_idx}`,
  `PileupColumn::{find_record, position_of}`, every `records_overlapping`, and
  `BamWriter::write_store_record`. `u32::MAX` is not a representable index, which makes
  `Option<RecordIdx>` free and retires the `u32::MAX` "no mate" sentinel on both mate fields.
  New: `RecordStore::try_record` (`None` past the store's end — `record` panics there, documented) and
  `RecordStore::indices()` in place of `0..store.len() as u32`.
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
- **`PileupEngine::new` takes a `PileupInput`, not a `RecordStore`.**
  The engine assumed ascending position order and linked mates and checked neither; both failures were
  silent (records after an out-of-order push were never pileuped; an unlinked store reports
  `mate_idx() == None` *and* `in_mate_overlap() == false`, which is what a read with no mate looks like,
  so overlap dedup quietly did nothing). `PileupInput` has no public constructor —
  `RecordStore::prepare_for_pileup()` is the only thing that mints one, returning it alongside
  `MateLinkStats`. Preparing is idempotent: the sort is skipped when nothing arrived out of order and the
  linking when nothing invalidated it, and `set_alignment` retracts the ordering property.
- **`PileupEngine::take_store` renamed to `reclaim_allocation`.**
  It returns an *empty* store that keeps its slab capacity for the next region; the old name promised
  the records back.
- **`PileupEngine::set_max_depth` takes a `NonZeroU32`.**
  It stored whatever `u32` it was given, so `set_max_depth(0)` meant "emit no alignments anywhere" —
  which a caller spelling unlimited as `0` asks for by accident. Unset still means no cap.
- **`AlignedPair::Insertion.qpos` renamed to `first_inserted`** (and likewise on `AlignedPairWithRead`
  and `AlignedPairWithRef`). An insertion is reported in two frames one base apart —
  `AlignedPair` points at the first inserted base, `PileupOp` at the matched base *preceding* the run —
  and both fields being `qpos` meant a pattern copied between them compiled and read one base off.
  Both frames are `QPos`, so only the name can catch it.
- **`PileupColumn::pair_indel` and `PairIndel` removed**, replaced by `PileupColumn::mate_of`.
  `pair_indel` baked in "the view's own indel wins and the mate is never consulted", which is wrong for
  any consumer whose filters can reject an indel the read does carry — the rule belongs with the filters.
- **VCF headers always declare `VCFv4.5`.**
  `VcfHeaderBuilder::file_format()` is removed and `VcfHeader::file_format()` returns `&'static str`;
  the version is the new `VcfHeader::FILE_FORMAT` constant. A cardinality is versioned (`Number=M` is
  4.5 and nothing earlier), so a header free to declare an older version could promise a grammar it then
  violates, undetectably. NB: `Number::BaseModification` is spec-legal under this header but noodles
  rejects `Number=M` at any declared version — `Number::Unknown` states the same cardinality readably.

### Added

- `PileupColumn::position_of(record_idx)` / `alignment_at(index)` — `find_record` split at its
  result, for a consumer whose own per-column scratch is addressed by position.
- `engine.pileups_into(&mut buf)`, `AlignmentView::inserted_bases` / `inserted_quals`, and anchor-addressed `PileupAlignment::indel_after`.
- `engine.set_soft_clip_overhang(n)` allows getting `n` soft-clipped bases at the fringes of alignments
- `SlimRecord::end_pos_htslib()` for htslib-compatible unmapped read spans.
- Byte-aware segment planning to bound per-segment memory.
- VCF FORMAT parity: `FormatInts`, `FormatString`, array-of-floats; percent-encode FORMAT string values;
  O(1) duplicate-field detection with htslib-compatible in-place overwrite.
- `Writer::finish` returns a self-serializing `CoordinateIndex`.
- `Pileup::mutate(mutator)` runs a caller-supplied mutation on the freshly fetched `RecordStore`
  (then re-sorts by position) — the hook for in-place local realignment via `RecordStore::set_alignment`.
  The mutator also receives the segment's `RefSeq`, the very one the engine later reports through
  `PileupColumn::reference_base`, so a hook that rescores or normalises alignments against the reference
  never loads its own copy; it covers the segment, not the reads, so read past its ends with
  `RefSeq::try_base_at`. The hook runs after the reference fetch, so it does not run when that fetch fails.
  `Pileup::with_reference(RefSeq)` drives the engine from a reference the caller already holds
  instead of re-reading the FASTA per segment, rejecting one that doesn't cover the segment
  (`ReaderError::SuppliedReferenceTooSmall`) — and a hook set alongside it reads that reference, not a
  second fetch.
- Window query over a prepared store: `PileupInput::records_overlapping(start, end)`, and the same on
  `PileupEngine` (between columns) and `PileupColumn` (while holding one), yield the indices of the mapped
  records whose alignment overlaps an inclusive span, ascending — a binary search over a running maximum
  of `end_pos` that `prepare_for_pileup` builds, so one long read only costs the windows it overlaps.
  Not offered on a bare `RecordStore`, which cannot prove it is in position order. `PileupInput::store()`
  gives read-only access to the store the indices address.
- Record-store mate linking for overlap-dedup consumers: `RecordStore::link_mates()` pairs a template's
  primary alignments once per store (qname-hash table, verified by qname bytes and reciprocal mate
  positions; secondary/supplementary and nameless reads never link). The pileup exposes
  `PileupAlignment::mate_idx()` / `in_mate_overlap()`, `AlignmentView::qname_hash()`,
  `PileupColumn::find_record()`, and `PileupColumn::mate_of()` — the linked mate of a view when that
  mate is also in this column, so a consumer can apply its own filters before deciding what the
  fragment says. The cached mate overlap is widened by the soft-clip overhang, so rescued
  soft-clipped fringe bases still count as overlapping.
- Indexed FASTQ references: the FAI parser accepts the sixth `qual_offset` column and `FaiEntry` exposes
  `is_fastq()` / `qual_byte_offset()`. Fetch is byte-correct on wrapped records whose quality lines may
  begin with `@`/`+`. htslib-validated.
- `PileupColumn::mate_of(view)`: the linked mate of an alignment when it is also in this column. The
  relation is symmetric and irreflexive, so a pairwise rule gives the same answer from either mate, and
  the mate is reachable whatever either read carries.
- `pileup_tiled` bench and `examples/tiled_pileup`, which measure the many-small-queries access pattern a
  variant caller actually has (per-query BAI lookup, `RegionBuf` setup, store sort and mate link) rather
  than one large region. `tools/samply_hot.py` symbolicates a `--save-only` samply profile from its
  `--unstable-presymbolicate` sidecar.

### Fixed

- CRAM: missing rANS Nx16 validation and allocation-size validation.
- Region queries no longer error on unplaced reads (`pos = -1` / `AP = 0`).
- Pileup: always set depth cap; stable `retain` replaces `swap_remove` in eviction.
- Region queries stop at the query end instead of reading past it. A record starting past the end proves
  nothing can follow in a coordinate-sorted file, but the loop counted it as merely out of range and kept
  inflating BGZF blocks from the ancestor bins BAI hands back, whose chunk lists run far beyond the
  region. On `tests/data/test.bam` a 1 kb query examined 2529 records to keep 313; it now examines 313.
- `RecordStore::dedup` sorts first instead of documenting that the caller must. `dedup` collapses
  *consecutive* equal records, so an unsorted store silently kept duplicates — exactly the case the
  method exists for (overlapping BAM index chunks loading one record twice).
- BCF genotypes: the first allele of a GT carries its phase bit. It was hardcoded unphased, which VCF 4.3
  and earlier say to ignore, but 4.4+ readers honour it — htslib rendered a fully phased `0|1` as `/0|1`.
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

- Column eviction uses `Vec::retain` instead of a manual swap-compaction: one `copy_nonoverlapping` per
  survivor instead of three moves of a 120-byte `ActiveRecord`, and no work at all when nothing was
  evicted (the common case between adjacent columns). 12.6 % of rastair's worker CPU before; 1.06x
  end-to-end there, byte-identical output.
- `PileupAlignment` no longer caches the qname hash — eight bytes written once per read per column
  (638 M times for 20 Mb of a 26x chromosome) to answer a question `mate_idx` already answers.
  `AlignmentView::qname_hash()` resolves it through the record instead. With `#[inline]` on
  `CigarMapping::deletion_after_at` (it returns `None` immediately for the ~96 % clips-match-clips reads,
  but was an out-of-line call costing 5.9 % of the read+pileup path) and `Base::known_index` as a table
  lookup rather than a data-dependent branch: column phase 3.80 s → 3.25 s on NA12878 chr12, 20 Mb.
- Region queries decompress about what htslib's iterator does: 45.8 GB → 22.7 GB on chr12 in 10 kb tiles
  (a single pass is 11.6 GB), from stopping at the query end.
- Column entries are written into one reserved block instead of pushed one at a time.
  `Vec::push` cannot hoist its capacity compare, and the loop already knows it emits at
  most one entry per active record. ~1 % of rastair's wall time on a 26x human chromosome.
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
