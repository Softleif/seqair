# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Added

- **CRAM byte data series may use the BETA encoding.** `FC`, `BA`, `QS`, `BS` and the values of
  `BYTE_ARRAY_LEN` accept `ByteEncoding::Beta`, as htslib does; the byte is the low 8 bits of
  `raw - offset`, matching htslib. BETA widths above 32 bits are now rejected with
  `CramError::InvalidBetaBits` for integer series too, as htslib rejects them.

## v0.3.1 (2026-09-25)

A performance release. Every hand-written `core::arch` kernel is now one portable `fearless_simd`
kernel, and the readers, pileup engine, CRAM codecs and writers were profiled and tightened around
it. Nothing public is removed or changed (`cargo semver-checks` against 0.3.0 is clean), and the
additions are opt-in. In rastair (chr12, byte-identical output) end-to-end cycles dropped 12 % on
a Ryzen 3950X and 1 % on an Apple M4, measured before the last rounds of CRAM, query and pileup
work.

### Fixed

- **BGZF writing failed on data that does not compress.** A block, gzip framing included, must fit
  in 64 KiB, and the writer filled blocks with 65,536 bytes of data. Data that does not compress is
  *stored*, a little larger than it went in, so every full block at `compression_level(0)` — any
  level-0 BAM, BCF or VCF.gz over 64 KiB — and any full block of random-looking bytes at other
  levels failed with `BgzfError::CorruptHeader`. Blocks now hold at most 65,280 bytes, htslib's
  `BGZF_BLOCK_SIZE`, which fits the worst case; a compile-time check keeps it that way. **Block
  boundaries move, so output is no longer byte-identical to 0.3.0's** (still valid, ~0.4 % more
  blocks).
- **Reading a block filled to exactly 64 KiB could read past the query's chunk.** Once such a block
  was read to its end, `RegionBuf` named the position as the block's own first byte — the
  within-block offset 65,536 wrapped to 0 — so the chunk-end check compared a position one block
  behind and kept going; a debug build panicked. 0.3.0's own writer produces such blocks. The
  position is now the next block's `(offset, 0)`, as the writer names it.
- **Bgzipped SAM queries took a BGZF error inside a chunk for the end of data** and returned short
  without an error. The error is returned now. They also no longer read one byte past a chunk end
  that falls on a block boundary.

### Added

- `BamWriterBuilder::compression_threads(n)` and `vcf::Writer::compression_threads(n)` compress
  BGZF blocks on `n` worker threads (`0`, the default, compresses on the calling thread), like
  htslib's `hts_set_threads`. The calling thread keeps serializing and writes blocks in order, so
  the data and the co-produced BAI/CSI are byte-identical to the single-threaded writer's. Writing
  10 Mb of 30x chr12 as BAM at level 6 with 3 threads: 4.0 → 1.34 s wall on the M4.
- `vcf::Writer::compression_level(level)` for `VcfGz` and `Bcf` (default 6), matching
  `BamWriterBuilder::compression_level`. Plain VCF ignores both settings.
- `BgzfError::ThreadSpawn`, `CompressionWorkerLost` and `BlockTooLarge`.
- `bam::seq::encode_bases_into`, and `RegionBuf::fill_buf` / `consume` for reading a region's
  decompressed bytes a block at a time.

### Changed

- **Plain VCF output is buffered** (128 KiB) instead of written with one `write_all` per record, so
  records reach the sink in large writes. `finish()` flushes; a writer dropped without it flushes
  best-effort and logs a failure with `warn!`.
- A BGZF block that cannot fit is reported as `BgzfError::BlockTooLarge { size }`, not
  `CorruptHeader`.
- New dependencies `fearless_simd` and `fearless_simd_macros`. The kernels pick the widest level
  the CPU has at run time, so AVX-512 machines now get 64-byte vectors, and the `unsafe` blocks of
  the old AVX2/SSSE3/NEON copies are gone. `cram::rans_nx16_avx2` and `cram::rans_nx16_neon` are
  now empty and hidden; they only ever held crate-private items.

### Performance

- **SIMD kernels** (3950X, time vs 0.3.0): ASCII→Base 0.70–0.78×, 4-bit sequence decode
  0.69–0.76×, plain-FASTA fetch 0.80–0.85×, the rANS Nx16 32-state order-0 kernel 0.29×.
  `OwnedBamRecord::to_bam_bytes` encodes the sequence in place (1.9× fewer cycles serialize-only).
  MM/ML base-modification resolution walks a SIMD bitmask of the target base: 5.7–7.6× faster on
  sparse calls, 1.8–2.4× on dense ones.
- **CRAM**, chr20:10–30 Mb (813k reads) end to end: CRAM 3.1 2.13 → 1.15 s on the 3950X (M4 1.08 →
  0.54 s), CRAM 3.0 1.86 → 1.34 s (M4 1.07 → 0.69 s) — before the last round of codec work, which
  took off several percent more. The rANS 4x8 and Nx16 decoders, PACK and RLE now use htscodecs'
  table layouts and loop shapes (rANS 4x8 order-0 331 → 513 MB/s, order-1 292 → 338; PACK up to
  2.1×). Reference runs are copied instead of converted base by base, codec buffers are reused, tag
  encodings are resolved once per slice, and external blocks are ordered by first use.
- **Small region queries**: a reader keeps a 64-block cache of decompressed BGZF blocks and its
  per-query setup between queries, starts a query where the previous one first found a record,
  and grows its reads from one block instead of the whole remaining range. BAI and CSI look up
  candidate bins instead of scanning every bin. chr12, M4, per query: 1 bp every 100 bp 277 →
  6.7 µs, 1 bp every 1 kb 293 → 43 µs, 1 kb tiles 334 → 53 µs; 1 Mb tiles unchanged. Bgzipped SAM
  shares the cache (1 kb tiles 1.30 → 0.50 ms) and scans lines and fields with SIMD (2.7× fewer
  cycles per fetch).
- **Pileup**: each read's CIGAR is walked with an incremental cursor (as htslib's `resolve_cigar2`),
  D and N columns take the fast path, columns where nothing expires skip eviction, and a column
  stops at `max_depth` instead of being built and truncated. Spliced synthetic pileup 2.81 →
  2.35 s; ~1000-deep columns capped at 50 1.85 → 0.29 s; a rastair-shaped TAPS pileup −16.5 %
  cycles on the M4.
- **Writers**: plain VCF of 2.76 M rastair-shaped records 5.2 → 0.83 s on the M4, from the
  buffering above and an O(1) FORMAT duplicate check. That is on top of `%g` floats printed with
  integer math instead of `core::fmt`'s Dragon4 fallback, which had already cut VCF text to 3.8×
  fewer cycles and VCF.gz to 2.1× fewer on the 3950X. Each BGZF block is assembled in place and
  written with one call.
- **Decoding**: BAM read-name and Z/H aux terminators are found with a vector scan (−2.8 % cycles
  decode-only), and BAM and CRAM errors are built only on failure instead of on every call (−8 %
  per small BAM query, −9 % on CRAM 3.1).

## v0.3.0 (2026-09-22)

### Breaking

Every genomic interval in the API is now one closed span, `core::range::RangeInclusive<Pos0>`,
instead of two loose positions (see `r[interval.span_type]`). The reference for a region is fetched
with the very value that queried it, so there is no `+ 1` at any boundary for a caller to get wrong,
and a span ending on `Pos0::MAX` names its last base where a half-open `end` could not.

- `fetch_into(tid, start, end, store)` → `fetch_into(tid, span, store)`, on `Readers`,
  `IndexedReader` and every format reader; same for `fetch_into_customized`,
  `estimate_region_bytes`, `IndexedBamReader::query` and `PileupEngine::new`. These were already
  inclusive on both ends: `(start..=end).into()` is the whole migration.
- **`IndexedFastaReader::fetch_seq` / `fetch_seq_into` and `Readers::fetch_base_seq` were half-open
  and are now closed.** `fetch_seq(name, p(0), p(4))` becomes `fetch_seq(name, (p(0)..=p(3)).into())`.
  `fetch_seq_into_u64`, the side door for reaching the last representable base, is gone — the
  closed span reaches it. `FastaError::RegionOutOfBounds` reports `last` (inclusive) instead of
  `end`, and fires for `start > last` or `last >= seq_len`.
- The segment target `(resolver, start, end)` is `(resolver, span)`.
- `Segment::end()` → `Segment::last()`: it returned the last covered position under a name that
  reads as one past it. `Segment::core_range()` → `core_span()`, now a `core::range::RangeInclusive`
  (`.start` / `.last`; `Copy`, 8 bytes, not an iterator). New `Segment::span()`.
- `Pos0::max_value()` → `Pos0::MAX`.
- MSRV 1.92.0 → 1.98.1. It was already effectively 1.98.1 via `seqair-types`; `core::range` needs
  1.96.
- `seqair-types` is depended on with `default-features = false` and a forwarding `serde` feature was
  added. seqair serializes nothing itself, so its users no longer compile serde for nothing; enable
  `seqair/serde` to serialize seqair-types values.

### Fixed

- **A reversed query returned the reads spanning it.** `fetch_into` with `start > end` handed the
  overlap test an interval that names no positions, and the test — `pos <= end && end_pos >= start`
  — is satisfied by exactly the records covering the whole gap. Every reader now treats a span with
  `last < start` as what it is, empty, and returns nothing. Found by a property test written for the
  new span API; the behaviour predates it.

- **A query on a truncated BAM never returned.** Not slowly — never. `IndexedBamReader::fetch_into`
  on a file truncated anywhere past its header, with its index left intact, spun forever with no
  panic, no allocation and no error, which is why 26 fuzz targets never saw it: libFuzzer could
  only ever have reported it as a timeout. `RegionBuf::read_record` reported two different facts as
  the same `BgzfError::UnexpectedEof` — "the window ran out mid-chunk, refill" and "the file ended
  here" — and the query loop read every one as the first, looping back on the assumption that the
  chunk step would happen on the next turn. At the end of the planned ranges the cursor stops
  advancing, so it never did. `read_record` now returns `Ok(None)` when the ranges are exhausted at
  a record boundary and keeps `UnexpectedEof` for running out partway through a record, and the
  chunk step is explicit, so every turn of the loop consumes either a record or a chunk. A
  truncated file now yields the records before the cut, as htslib does with one.
- **CRAM: a multi-reference container could decode reads as `N` with no error.** Such a container
  holds one index entry per slice per reference, and the reference window was taken from the first
  entry whose container offset matched. When a later slice reached further along the reference, the
  window stopped short and the tail of those reads came back as `N` — silently, since running past
  the fetched reference is only a warning. It is the union of every matching entry now. htslib turns
  multi-reference slices on by itself once a container would hold few records per reference, so an
  ordinary file with short contigs reaches this; two slices one base apart are enough.
- **CRAM: `embed_ref=2` files could not be opened at all.** The slice-header MD5 was checked against
  the FASTA whether or not the slice carried its own reference. Under `embed_ref=2` htslib embeds a
  *consensus* computed from the reads and digests that, so it matches the external reference only
  where the reads happen to agree with it — every such file with low coverage or a real difference
  failed with `ReferenceMd5Mismatch`. A slice with an embedded reference is now exempt.
- **CRAM: a no-reference file mis-decoded any read with an insertion.** The `Q` and `q` features
  carry quality and nothing else, but both were decoded as an anchoring reference match, adding a
  base and an `M` operation. In no-reference mode htslib emits one `Q` per *inserted* base next to
  the `I` feature, so those reads came back one base too long and were refused with
  `QualLenMismatch`. It needs no option to reach: htslib drops into no-reference mode by itself when
  `embed_ref` meets `multi_seq_per_slice`.
- **A CSI could silently drop records from a region query.** A bin's `loffset` was written as
  that bin's own first chunk offset. htslib derives it from the linear index instead — the first
  record at or after the start of the bin's leftmost leaf window — and the two differ whenever a
  record in another bin begins earlier inside that window, which a record straddling a window
  boundary does routinely. The bin's own chunk is then the larger of the two, and that is the
  dangerous direction: a reader takes `min_off` from `loffset` and discards every chunk ending
  before it, so `tabix` and `bcftools view -r` both returned short, with exit status 0. The record
  was never lost from the file and a query for its exact position still found it; only a query
  whose range started earlier missed it. Smallest case: a record at 0-based 16384 with a second
  straddling the next window, queried from position 1.
- **The VCF/BCF writer could not index a contig over 512 Mbp.** Its index was built with
  `min_shift=14, depth=5` whatever the header said, so bins ran out at 2^29 — a record above that
  was written to the file, pushed to the index, and then unreachable, with `bcftools view -r`
  returning nothing while the same query over `bcftools index -c`'s CSI returned it. Nothing
  errored. The depth is now derived from the longest reference, as htslib's
  `hts_adjust_csi_settings` does. This is the case CSI exists for; TBI and BAI cannot express it.
  Indexes that already fitted the depth-5 scheme are byte-identical, since the search floors at 5.
- **`OwnedBamRecord::to_bam_bytes` accepted a record whose quality no longer matched its sequence.**
  `set_seq` validates the new sequence against the CIGAR and deliberately leaves `qual` alone, so
  `set_seq` with no following `set_qual` was a reachable state where the two disagreed. It
  serialized: `l_seq` governs how many quality bytes a reader consumes, so the record did not
  bounce, it was misread from the sequence onwards. It now raises the same `SeqQualLengthMismatch`
  the builder and `set_qual` raise. An empty `qual` remains legal at any length.
- **Non-finite floats in VCF text now match htslib's spelling.** `write_float_g` sent NaN and the
  infinities through `Display`, writing `NaN` where C's `%g` — and so htslib — writes `nan`. VCF 4.3
  §1.3 admits `INF`/`INFINITY`/`NAN` case-insensitively as Float values and §6.3.3 gives quiet NaN
  first-class status, distinct from the missing sentinel, so these are values rather than errors and
  are deliberately *not* written as `.`: that would turn a value into a missing value and put the
  text output at odds with the BCF output, which writes the caller's bits through unchanged.

- **The SAM reader rejected every negative element of a signed `B` array.** `B:c`, `B:s` and
  `B:i` are the signed subtypes, but all three were converted through the unsigned type of the
  same width, so `XX:B:c,-1` failed the whole record with `InvalidAuxValue`. Each subtype is now
  converted through the type it names. An unrecognised subtype is also an error rather than a
  `B` tag whose element count promises more bytes than follow it.
- **VCF float text now matches `bcftools` exactly.** `write_float_g` claimed to write C's `%g`
  with six significant digits — the format htslib emits — but only ever produced the fixed form,
  so a value outside `1e-4 .. 1e6` came out as `0.0000610352` or `1234567` where bcftools writes
  `6.10352e-05` and `1.23457e+06`. It now switches forms on the value's exponent *after* rounding
  to six significant digits (so `0.0001` stays `0.0001` and does not become `1e-04`), writes the
  exponent C's way with a sign and at least two digits, and writes negative zero as `-0`. Values
  in `1e-4 .. 1e6` are unchanged, which is every quality score and most everything else.
- **A float below about `1e-25` could not be written to VCF at all.** `write_float_g` sized its
  decimal places from the value's magnitude and formatted into a 32-byte buffer, so six significant
  digits of `5.169879e-26` overflowed it and the write failed with `FormattedFloatLongerThan32Chars`,
  taking the whole record with it. Those magnitudes do not occur in QUAL, but they do in a p-value
  carried in INFO. The writer now falls back to scientific notation there, which is what `%g` and
  htslib emit. Values that already fit are byte-identical to before.

### Added

- `IndexBuilder::csi` builds a CSI whose depth covers a given longest-reference length, and
  `IndexBuilder::csi_depth_for` exposes that calculation on its own. `IndexBuilder::BAI_DEPTH` names
  the 5 that BAI and TBI are fixed at and that CSI now treats as a floor.
- `seqair::bam::DecodeError` is re-exported. `RecordStore::set_alignment`, `push_raw` and
  `push_fields` are public and return it, but it was only reachable through a private module, so
  a caller outside the crate could not match on the error.

### Changed

- `RegionBuf::read_record` returns `Result<Option<&[u8]>, BgzfError>`. `Ok(None)` means the planned
  ranges are exhausted at a record boundary — there is no next record; `Err(UnexpectedEof)` now
  means only that a record was cut in half. Callers that treated the old single error as "retry"
  must handle `Ok(None)` as an end, not a hiccup; see the truncation fix above for why.

### Internal

- Property-based tests moved from `proptest` to [hegel](https://hegel.dev/) (`hegeltest`). Case
  counts live in `hegel.toml` at the workspace root: 256 locally, 1000 on CI, and a `thorough`
  profile at 10000. No public API change; `proptest` is gone from the dev-dependencies.
- `BgzfWriter`'s virtual offsets, the two `CustomizeRecordStore` filter hooks, and the
  `Pos`/`QPos` type walls all gained property tests. The BGZF ones were checked against the
  offset bug fixed in the previous release: reverting it makes them fail on a two-write case.
- Six cross-implementation comparisons that ran on fixed fixtures now run on generated inputs:
  the pileup against htslib's `bam_plp_auto` (and now comparing per-alignment `qpos`, not only
  depth); the write → index → query round-trip against both `samtools view` and the generated
  records; CRAM against BAM decodes of the same reads; multi-sample FORMAT arrays against
  bcftools; SAM through the reader and `write_store_record` and back out through samtools; and
  records built with `OwnedBamRecord::builder` through `BamWriter` and back.
- Generated coverage now reaches the corners that were fixture-only, and found the six bugs listed
  above: the CRAM writer-option matrix (version × `embed_ref` × `seqs_per_slice` ×
  `slices_per_container`, with multi-reference slices and span-0 index entries asserted *per case*
  rather than counted across the run), `OwnedBamRecord` under sequences of mutations refereed by
  samtools, reader-level filtering compared field by field against the unfiltered fetch, extras
  through the pileup, INFO types and missing-value integer arrays against bcftools, CSI output
  including contigs past 2^29, mate fields refereed by `samtools fixmate`, and which typed error a
  corrupt or truncated file actually surfaces.
- Spec traceability: 13 dangling rule ids resolved to zero — six were renames, seven rules the code
  was annotated against but nobody had written — one stale reference fixed, and 28 rules in the
  `bcf_writer.*` and `csi.*` namespaces annotated against round-trips that already pinned them.
  Five rules were deliberately left unverified rather than annotated on tests that would not fail if
  they broke. Three rules whose prose had drifted from the code (`record_store.customize.trait`,
  `perf.reuse_alignment_vec`, `unified.segment_options`) were retold to match it.
- The six bare `compile_fail` blocks in the VCF writer each gained a passing twin that uses every
  name the rejection uses, so a rename fails the guard rather than silently making the rejection
  free; each was mutated into its legal form to confirm it rejects for the reason it claims.
- The nightly workflow (renamed `fuzz.yml` → `nightly.yml`) grew from one job to three, because
  there are three axes and they do not substitute for one another: `fuzz` as before, `thorough` at
  the 50000-case profile, and `deep-oracles` with `HEGEL_TEST_CASES=500`. The last is the only
  thing that reaches the 83 tests which pin their own case count — the subprocess round-trips
  against samtools, bcftools, htslib and noodles — and it found the CSI `loffset` bug above on its
  first outing. Both numbers are measured, not guessed, and the curves are recorded in `hegel.toml`
  and the workflow comments.

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
