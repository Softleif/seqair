# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Fixed

- **A float below about `1e-25` could not be written to VCF at all.** `write_float_g` sized its
  decimal places from the value's magnitude and formatted into a 32-byte buffer, so six significant
  digits of `5.169879e-26` overflowed it and the write failed with `FormattedFloatLongerThan32Chars`,
  taking the whole record with it. Those magnitudes do not occur in QUAL, but they do in a p-value
  carried in INFO. The writer now falls back to scientific notation there, which is what `%g` and
  htslib emit. Values that already fit are byte-identical to before.

### Internal

- Property-based tests moved from `proptest` to [hegel](https://hegel.dev/) (`hegeltest`). Case
  counts live in `hegel.toml` at the workspace root: 256 locally, 1000 on CI, and a `thorough`
  profile at 10000. No public API change; `proptest` is gone from the dev-dependencies.

## v0.2.0 (2026-09-14)

A streaming-window BAM region reader, newtyped record and query indices, a unified
pileup API, and VCF/BCF encoder ergonomics.

### Breaking

**Pileup API**

- `Readers::pileup(segment, depth)` returns a `Pileup` plan; finish it with `.run()`.
  `pileup_with` / `pileup_with_reference` → `.mutate(f)` / `.with_reference(r)` on the plan,
  combinable in any order.
- `PileupEngine::new` takes a `PileupInput`, minted only by `RecordStore::prepare_for_pileup()`
  (which also returns `MateLinkStats`). The engine used to silently assume position order and
  linked mates.
- `PileupEngine::take_store` → `reclaim_allocation` — it returns an _empty_ store that keeps its
  slab capacity, not the records.
- `PileupEngine::set_max_depth` takes `NonZeroU32`; `0` used to mean "emit nothing".
- `PileupAlignment::strand` removed — use `Strand::from(rec.flags)`.
- `PileupColumn::pair_indel` and `PairIndel` removed → `PileupColumn::mate_of`, so your own
  filters decide what the fragment says.

**Record store**

- Record indices are the `RecordIdx` newtype (`bam::record_idx`), not bare `u32`. `u32::MAX` is
  unrepresentable, so `Option<RecordIdx>` is free and the `u32::MAX` "no mate" sentinel is gone.
- Records are read through `RecordStore::record(idx) -> Option<RecordRef<'_>>`, which does not panic.
  The by-index readers (`try_record`, `qname(idx)`, `cigar(idx)`, `seq(idx)`, `seq_at`, `qual(idx)`,
  `aux(idx)`, `extra(idx)`, `mate_overlap(idx)`) are methods on the handle instead — `rec.qname()`,
  `rec.cigar()`, `rec.base_at(qpos)`, … — with the record's own fields via `Deref`. The handle borrows
  the store, so `clear`, pushes, `sort_by_pos` and `dedup` are rejected while one is alive.
  `set_alignment` / `write_store_record` now report `NoSuchRecord`; `extra_mut` returns `Option`.
  Use `RecordStore::indices()` instead of `0..store.len() as u32`.

**Coordinates**

- Query offsets are the `QPos` newtype, deliberately not interconvertible with `Pos`. Affects
  `AlignedPair`, `MatchPosition` / `MatchedBase` / `MatchedRef`, `CigarPosInfo`, `PileupOp`,
  `PileupAlignment::qpos()`, `CigarMapping::soft_clip_qpos_at`, and
  `BaseModState::mod_at_qpos` / `is_unmodified`.
- `AlignedPair::Insertion.qpos` → `first_inserted` (also on `AlignedPairWithRead` / `WithRef`):
  it names the _first inserted_ base, while `PileupOp::Insertion.qpos` names the matched base
  _before_ the run.

**Readers and writers**

- `RegionBuf` streams a bounded sliding window instead of loading whole regions: it borrows the
  reader (so `RegionBuf<'a, _>`), no longer implements `UnwindSafe`, and `load` is gone
  (`ensure_available` / `advance_range`).
- Unmapped-read filtering unified on `filter_raw`: `IndexedBamReader::keep_unmapped` / `keeps_unmapped`
  removed, `FilterRawFields` is type-level enums (public `raw_cigar_bytes`, `cigar_ops`, `packed_seq`,
  `bases` fields gone; `end_pos` is `Option<Pos0>`), new `RejectUnmapped` customizer.
- `TargetInfoAccess` removed from `seqair::bam::header`.
- `ReaderError` gained a variant separating "region end past contig" from "exceeds `i32::MAX`".

**VCF/BCF**

- VCF headers always declare `VCFv4.5`: `VcfHeaderBuilder::file_format()` removed,
  `VcfHeader::file_format()` returns `&'static str`, version is `VcfHeader::FILE_FORMAT`.
  A cardinality is versioned, so a header free to declare an older version could promise a grammar
  it then violates.
- `EncodeInfo` and `EncodeFormat` removed; `FormatEncoder` gained three required methods
  (`format_ints`, `format_floats`, `format_string`).
- `VcfHeaderError::TooManyFields` now carries data.

### Added

- **Reference hooks on the pileup plan**: `Pileup::mutate(f)` rewrites the fetched `RecordStore`
  before pileup (local realignment via `set_alignment`) and receives the segment's `RefSeq` — the
  same one `PileupColumn::reference_base` reports. `with_reference(RefSeq)` drives the engine from a
  reference you already hold; `reference_covers_reads()` widens the fetch to the span the records
  actually cover.
- **Mate linking**: `RecordStore::link_mates()` pairs a template's primary alignments once per store.
  Exposed as `PileupAlignment::mate_idx()` / `in_mate_overlap()`, `AlignmentView::qname_hash()`,
  `PileupColumn::find_record()` / `mate_of()`.
- **Window queries** over a prepared store: `records_overlapping(start, end)` on `PileupInput`,
  `PileupEngine` and `PileupColumn` — a binary search over a running `end_pos` maximum, so one long
  read only costs the windows it overlaps. Plus `PileupInput::store()`.
- `PileupColumn::position_of(record_idx)` / `alignment_at(index)`;
  `AlignmentView::inserted_bases` / `inserted_quals`; `PileupAlignment::indel_after`;
  `engine.set_soft_clip_overhang(n)` to keep `n` soft-clipped bases at alignment fringes.
- Indexed FASTQ references: the FAI parser accepts the sixth `qual_offset` column;
  `FaiEntry::is_fastq()` / `qual_byte_offset()`. htslib-validated.
- VCF FORMAT parity: `FormatInts`, `FormatString`, array-of-floats, percent-encoded string values,
  htslib-compatible duplicate-field overwrite. `Writer::finish` returns a self-serializing
  `CoordinateIndex`.
- Byte-aware segment planning to bound per-segment memory.

### Fixed

- **`BgzfWriter::virtual_offset()` could name a byte _inside_ a record** when a write filled a block
  exactly, corrupting every co-produced CSI/BAI/TBI offset that landed there — htslib rejected such
  files outright. ~0.4 % of blocks fill exactly. Compressed bytes are unchanged; only the recorded
  offsets differ.
- **Region queries stop at the query end** instead of inflating BGZF blocks far past it: a 1 kb query
  on `tests/data/test.bam` examined 2529 records to keep 313; now 313.
- **SAM and CRAM region queries lost records overlapping by the query's last base** (SAM tested a
  half-open range against an exclusive `end_pos`; CRAI compared a 1-based start to a 0-based query).
  All readers now share one convention — 0-based, both ends inclusive.
- **`RecordStore::dedup` sorts first.** It collapses _consecutive_ equal records, so an unsorted store
  silently kept the duplicates the method exists to remove.
- **BCF genotypes carry the phase bit on the first allele.** Hardcoded unphased, so htslib rendered a
  fully phased `0|1` as `/0|1`.
- `SlimRecord::end_pos` is htslib-compatible on unmapped reads (`end_pos == pos`).
- Region queries no longer error on unplaced reads (`pos = -1` / `AP = 0`).
- FASTA: an untrusted `.fai` can no longer drive a multi-PiB allocation, serve bytes from a wrapped
  offset, or panic on `linebases == 0`; errors carry the path and distinguish `PlainGzipUnsupported`
  from `GziIndexNotFound`.
- CRAM: missing rANS Nx16 and allocation-size validation.
- Empty-sample BCF missing-value handling; pileup depth cap always set.

### Performance

- Region queries decompress about what htslib's iterator does: 45.8 GB → 22.7 GB on chr12 in 10 kb
  tiles, from stopping at the query end.
- Pileup columns: `Vec::retain` eviction, entries written into one reserved block, pooled scratch
  buffers, no cached qname hash, `Base::known_index` as a table lookup. Column phase 3.80 s → 3.25 s
  on NA12878 chr12 (20 Mb); 1.06x end-to-end in rastair, byte-identical output.
- Plain-FASTA fetches use one positional `pread` per span (~28 % faster on short slices) and forks
  share a handle; fetched spans are stripped and uppercased on a vectorized path.
- `CompactOp` 16 → 12 bytes; `CigarSlice` collapsed to `&[CigarOp]`.
- New `pileup_tiled` bench and `examples/tiled_pileup` measure the many-small-queries pattern a
  variant caller actually has.

## v0.1.0 (2026-05-08)

First ever release. Changelog entries will start in next one.
