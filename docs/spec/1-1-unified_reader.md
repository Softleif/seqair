# Unified Reader

seqair reads BAM, bgzf-compressed SAM, and CRAM files through a format-agnostic reader interface that auto-detects the input format and dispatches to the appropriate parser, while presenting a uniform API to the rest of the codebase.

> **Sources:** The unified reader design is seqair-specific. Format magic bytes come from [SAM1] §4.2 (BAM: `BAM\1`), [SAM1] §4.1 (BGZF: `1f 8b`), and [CRAM3] §6 (CRAM: `CRAM`). Sort order detection uses the `SO` tag from [SAM1] §1.3 "The header section" (`@HD` line). See [References](./99-references.md).

## Goals

1. **Backward compatibility**: existing code that calls `IndexedBamReader::open()` continues to work unchanged.
2. **Auto-detection**: `IndexedReader::open(path)` inspects the file to determine BAM vs SAM vs CRAM and opens the appropriate reader.
3. **Uniform downstream API**: all three formats populate the same `RecordStore` via `fetch_into()`. The pileup engine, call pipeline, and all downstream code are format-agnostic.
4. **Forkable**: all reader types support `fork()` for thread-safe parallel processing with shared immutable state.

## Format detection

> r[unified.detect_format]
> `IndexedReader::open(path)` MUST auto-detect the file format by inspecting magic bytes:
>
> - Bytes `1f 8b` (gzip magic) → verify BGZF structure (extra field with `BC` subfield). If not BGZF, reject with error: "file is gzip-compressed but not BGZF; use `bgzip` instead of `gzip`." If BGZF, decompress the first block:
> - Starts with `BAM\x01` → BAM format
> - Starts with `@` (0x40) → bgzf-compressed SAM
> - Otherwise → reject with error naming supported formats
> - Bytes `43 52 41 4d` (`CRAM`) → CRAM format
> - Byte `@` (0x40) at position 0 → uncompressed SAM. Reject with error: "uncompressed SAM files cannot be indexed; compress with `bgzip file.sam` then index with `samtools index file.sam.gz`."
> - Otherwise → return an error listing supported formats (BAM, bgzf-compressed SAM, CRAM)

> r[unified.detect_index]
> After format detection, the reader MUST locate the appropriate index file:
>
> - BAM: `.bai`, `.bam.bai`, or `.csi` (CSI supports references > 512 Mbp that BAI cannot)
> - SAM: `.tbi` or `.csi` (tabix or CSI index)
> - CRAM: `.crai` (CRAM index)

r[unified.detect_error]
If the format is detected but no matching index is found, the error MUST name the expected index extension and suggest the tool to create it (`samtools index` for BAM/CRAM, `samtools index` or `tabix` for SAM).

## Reader enum

r[unified.reader_enum]
The unified reader MUST be an enum dispatching to format-specific readers, not a trait object. This avoids dynamic dispatch overhead in the hot path and keeps the type concrete for `fork()`.

> r[unified.reader_api]
> The unified reader MUST expose:
>
> - `header() -> &BamHeader` — all formats produce the same header type (target names, lengths, tid lookup). The SAM and CRAM parsers convert their native header representations to `BamHeader`.
> - `fetch_into(tid, start, end, store) -> Result<usize>` — populates a `RecordStore` with records overlapping the region, identically to the BAM path.
> - `fork() -> Result<Self>` — creates a lightweight copy sharing immutable state.
> - `shared() -> &Arc<_>` — access to shared state for `Arc::ptr_eq` testing.

r[unified.fetch_equivalence]
For a given BAM file and its SAM/CRAM representations of the same data, `fetch_into` for the same region MUST produce records with equivalent logical content: same positions, flags, sequences, qualities, and CIGAR operations. Aux tags MUST contain the same set of tag names and values, but tag ordering and integer type codes (e.g., `c` vs `i` for small values) MAY differ between formats because BAM writers choose specific integer widths that SAM text cannot preserve.

## RecordStore integration

> r[unified.record_store_push]
> Currently `RecordStore::push_raw()` takes raw BAM bytes. For SAM and CRAM, records arrive as parsed fields rather than BAM binary.
>
> The RecordStore MUST provide `push_fields()` (or equivalent) that accepts pre-parsed record fields:
>
> - pos, end_pos, flags, mapq, seq_len, matching_bases, indel_bases (fixed fields)
> - qname bytes
> - CIGAR as packed BAM-format u32 ops (SAM/CRAM parsers convert to this representation)
> - sequence as `&[Base]` (already the enum type stored in the bases slab — no conversion needed)
> - quality bytes
> - aux tag bytes in BAM binary format (SAM/CRAM parsers serialize tags to this format)

This avoids a wasteful round-trip through BAM binary encoding.

r[unified.push_fields_equivalence]
`push_fields()` and `push_raw()` MUST produce identical `SlimRecord` fixed fields and identical name/bases/cigar/qual slab contents for the same logical record. Aux tag slab contents MAY differ in integer type codes (a SAM-derived `i:42` may serialize as `c` while the original BAM used `C`), since the original BAM writer's type choice is not recoverable from SAM text. Tests MUST verify fixed-field and non-aux-slab equivalence by pushing the same record both ways.

## Shared header type

r[unified.header_from_sam_text]
SAM headers are text (the `@HD`, `@SQ`, `@RG`, `@PG`, `@CO` lines). CRAM containers embed SAM header text. `BamHeader` already stores `header_text: String` and parses `@SQ` lines for targets.

`BamHeader` MUST gain a `from_sam_text(text: &str)` constructor that parses SAM header text directly, without requiring BGZF/BAM framing. This is used by both the SAM reader (header lines at start of file) and the CRAM reader (SAM header block in file header container).

## Sort order validation

r[unified.sort_order]
Indexed random access assumes coordinate-sorted data. If the `@HD` header line contains `SO:unsorted` or `SO:queryname`, the reader MUST return an error explaining that indexed region queries require coordinate-sorted input. `SO:coordinate` or absent `SO` (common in older files) MUST be accepted.

## Fork semantics per format

r[unified.fork_bam]
BAM fork: shares `Arc<BamShared>` (index + header), opens fresh `File` handle.

r[unified.fork_sam]
SAM fork: shares `Arc` holding parsed tabix/CSI index + header, opens fresh `File` handle for BGZF reading.

r[unified.fork_cram]
CRAM fork: shares `Arc` holding CRAI index + header, opens fresh `File` handle for container reading. When the source reader carries a FASTA, each fork MUST also have its own FASTA reader (via `fasta_reader.fork()`) for thread-safe reference lookups. When the source reader was opened without a FASTA (`r[cram.fasta.optional]`), the fork MUST also have no FASTA. Reference caching (`r[cram.perf.reference_caching]`) MUST be per-fork, not shared.

## Readers: alignment + reference bundle

r[unified.readers_struct]
The `Readers` struct bundles an `IndexedReader` (alignment) with an `IndexedFastaReader` (reference) in a single type. This eliminates the need for separate FASTA path passing — CRAM can access the reference it needs for sequence reconstruction, and all formats have uniform open/fork semantics. `Readers` is parameterized on `E: CustomizeRecordStore = ()` so callers can attach a per-record customize value at open time (see `r[record_store.customize.trait]`); the customize value's `keep_record` runs at fetch time and `compute` runs by `pileup()` after every fetch.

r[unified.readers_open]
`Readers::open(alignment_path, fasta_path)` MUST auto-detect the alignment format (via `r[unified.detect_format]`) and open the appropriate reader. The `fasta_path` parameter MUST accept both a `&Path` (FASTA attached) and `Option<&Path>` (`None` for no FASTA, `Some(path)` for FASTA) — concretely, the parameter type is `impl Into<Option<&Path>>`, so `Readers::open(bam, &ref_fa)`, `Readers::open(bam, Some(&ref_fa))`, and `Readers::open(bam, None)` MUST all compile and work. With a FASTA, the path is passed to `IndexedCramReader::open()` for CRAM sequence reconstruction; for BAM/SAM the FASTA reader is opened but not used during decoding. With no FASTA, the alignment reader is opened reference-free: BAM and SAM never need a reference for record decoding, and CRAM follows `r[cram.fasta.optional]`. The returned `Readers` then carries `fasta = None`, so `pileup()` MUST skip the FASTA fetch and the resulting `PileupColumn::reference_base()` MUST be `Base::Unknown` for every column; callers that later need bases for a CRAM slice without an embedded reference will see `CramError::MissingReference` propagated through `ReaderError::Cram`. `Readers::open` is only available when `E = ()`; for non-trivial customizers use `open_customized`.

r[unified.readers_open_customized]
`Readers::<E>::open_customized(alignment_path, fasta_path, customize)` MUST open the readers the same way as `Readers::open` (including the `impl Into<Option<&Path>>` polymorphism for the FASTA argument) but carry the user-supplied customize value along for `pileup()` to invoke on each fetched region. The customize value is stored by ownership inside the `Readers` struct.

r[unified.readers_fork]
`Readers::fork()` MUST fork both the alignment reader and the FASTA reader, returning a new `Readers` with independent I/O handles but shared immutable state (indices, headers). The CRAM fork gets its own FASTA reader via `IndexedFastaReader::fork()`. When `E ≠ ()`, `fork()` MUST clone the customize value into the new `Readers` (the `Clone` bound on `CustomizeRecordStore` guarantees this is cheap).

r[unified.tid.newtype]
A validated `Tid(u32)` newtype MUST wrap BAM target ids so downstream code cannot pass an unvalidated `u32`. Construction is only through the `ResolveTid::resolve_tid(header)` method, which MUST be implemented for `u32` (range-check), `&str`/`String`/`SmolStr` (by-name lookup), and `Tid` (passthrough). `ResolveTid::resolve_tid` MUST return a typed `TidError` on failure (unknown contig name or out-of-range `u32`).

## Segments and pileup

The pileup pipeline is intentionally split into a *planning* step (`segments()`)
and an *execution* step (`pileup()`). The planning step produces a `Segment` —
a small, pre-resolved tile of the genome. The only way to drive a `pileup()`
is to hand it a `Segment`, which makes it impossible to silently fetch a whole
chromosome in one call without thinking about tile size.

> r[unified.segment_struct]
> The crate MUST expose a `Segment` type representing one bounded tile of a
> genomic region. `Segment` MUST:
>
> - Be constructible only by `Readers::segments` (i.e. no public constructor).
> - Carry a validated `Tid`, the contig name (`SmolStr`), an inclusive
>   `[start, end]` range as `Pos0`, an `overlap_start` and `overlap_end`
>   (`u32`, in bases), and the contig's last valid 0-based position
>   (`contig_last_pos: Pos0`).
> - Implement `Clone + Debug + PartialEq + Eq + Hash` and be `Send`. It MUST
>   NOT carry FASTA bases or any per-fetch buffer (so segments can be cheaply
>   collected, shipped between threads, or stored).
> - Expose accessors: `tid() -> Tid`, `contig() -> &SmolStr`,
>   `start() -> Pos0`, `end() -> Pos0` (inclusive), `len() -> u32`
>   (= `end - start + 1`), `overlap_start() -> u32`, `overlap_end() -> u32`,
>   `contig_last_pos() -> Pos0`, `starts_at_contig_start() -> bool` (true
>   iff `start == Pos0::ZERO`), `ends_at_contig_end() -> bool` (true iff
>   `end == contig_last_pos`). The "contig" predicates check the contig
>   boundary, not the requested target range — for a sub-range query the
>   user's first/last segments do not necessarily satisfy them.

> r[unified.segment_overlap]
> Overlap fields on `Segment` MUST be expressed as a number of bases shared
> with neighboring segments and default to 0 when no overlap is requested.
>
> Overlap rules apply to the **requested target range**, not the contig:
> if a caller asks for `("chr1", 1000, 2000)` the iterator MUST NOT extend a
> segment to position 900 just because there are bases there — the user did
> not ask for them.
>
> - The first segment of every requested target range MUST have
>   `overlap_start == 0`.
> - The last segment of every requested target range MUST have
>   `overlap_end == 0`.
> - Internal segments (segments where the previous/next core is fully inside
>   the requested range) MUST carry `overlap_start == overlap_end ==
>   requested_overlap`.
> - When a tile's expanded right edge would extend past the requested
>   range's `end`, the actual `overlap_end` MUST equal the clamped distance
>   `tile_end - core_end`, which is `< requested_overlap`. Same on the left
>   for `overlap_start`.
>
> `Segment` MUST expose `core_range() -> RangeInclusive<Pos0>` returning the
> inclusive sub-range of `[start, end]` excluding both overlaps — this is
> the part of the segment a downstream tool should treat as "owned" by this
> tile when deduplicating against neighbors.

> r[unified.segment_options]
> `SegmentOptions` MUST expose:
>
> - `SegmentOptions::new(max_len: NonZeroU32) -> Self` — sets `overlap = 0`.
> - `with_overlap(self, overlap: u32) -> Result<Self, SegmentOptionsError>` —
>   returns `Err(SegmentOptionsError::OverlapTooLarge { max_len, overlap })`
>   when `overlap >= max_len.get()` (would produce zero forward progress).
> - `max_len() -> NonZeroU32`, `overlap() -> u32` accessors.
>
> `SegmentOptions` MUST NOT implement `Default` — every caller is required to
> commit to a `max_len`, since the absence of a sensible universal default is
> the whole point of forcing a planning step.

> r[unified.into_segment_target]
> `IntoSegmentTarget` MUST be a **sealed** trait (external crates cannot
> implement it) that resolves a user-supplied target against the header
> into one or more contig ranges. Implementations MUST be provided for at
> least:
>
> - `&str` / `String` / `SmolStr` / `Tid` / `u32` — entire contig, half-open
>   `[0, contig_len)` mapped to inclusive `[Pos0::ZERO, contig_last_pos]`.
> - `&RegionString` — parsed `chrom:start-end`, with missing start defaulting
>   to position 1 (0-based 0) and missing end to the contig's last base.
> - `(R, Pos0, Pos0)` where `R: ResolveTid` — explicit inclusive range
>   `[start, end]` against the resolved tid.
> - `()` — every contig in header order whose length is non-zero.
>
> An empty contig (length 0) MUST return `ReaderError::EmptyContig` for
> single-target conversions and MUST be silently skipped for `()` (whole-genome).
> Resolution failures (unknown contig, out-of-range tid, start > end,
> end > contig_last_pos) MUST surface as typed `ReaderError` variants — never
> by silently clamping.

> r[unified.readers_segments]
> `Readers::segments(target: impl IntoSegmentTarget, opts: SegmentOptions)`
> MUST return `Result<impl Iterator<Item = Segment>, ReaderError>`. The
> iterator MUST:
>
> - Cover every base of every input range exactly once across `core_range()`
>   of consecutive tiles — i.e. the union of `core_range()` of all yielded
>   segments equals the union of input ranges, with no overlap and no gap.
> - Produce tile **cores** (`Segment::core_range()`) of length `<= max_len`,
>   never splitting a tile across contigs. Each tile's full `[start, end]`
>   is its core expanded by `overlap` bases on each side, clipped to the
>   requested target range; internal tiles therefore have
>   `len() == max_len + 2 * overlap`, while edge tiles are shorter.
> - Set `overlap_start` / `overlap_end` per `r[unified.segment_overlap]`. The
>   actual `[start, end]` of every tile is its core range expanded by
>   `overlap` bases on each side, clamped to the **requested target range**
>   (not to the contig). At the requested range's edges this clamping
>   reduces the corresponding overlap field to 0.
> - When the remaining contig after a tile's core would be `<= overlap` bases
>   long, the next tile's core MAY be shorter than `max_len - 2*overlap`; in
>   the limiting case the last tile of a contig is one tile covering the
>   remainder. The iterator MUST NOT yield zero-length tiles.
> - Yield tiles in `(tid, start)` order. For multi-contig targets, all tiles
>   for one contig MUST be emitted contiguously before moving on.
> - Borrow `&self` only — building the iterator MUST NOT mutably borrow
>   `Readers`, so callers can build the plan once and re-acquire `&mut self`
>   for each `pileup(&segment)` call.

> r[unified.readers_pileup]
> `Readers::pileup(segment: &Segment) -> Result<PileupEngine<E::Extra>, ReaderError>`
> is the only entry point for driving a pileup. It MUST:
>
> 1. **Validate header consistency.** Look up `segment.contig()` against the
>    current header. If `header.tid(segment.contig().as_str())` does not
>    return the same numeric tid stored in `segment.tid()`, return
>    `ReaderError::SegmentHeaderMismatch { contig, expected_tid }`.
>    Additionally, if `header.target_len(tid) - 1` does not equal
>    `segment.contig_last_pos().as_u64()`, return
>    `ReaderError::SegmentContigLengthMismatch { contig, segment_last_pos,
>    header_last_pos }`. Together these catch the foot-gun where a
>    `Segment` built against one [`Readers`] is fed to a different
>    [`Readers`] whose header doesn't match — both contig-order and
>    contig-length differences across reference panels are surfaced as
>    typed errors instead of silent zero-record fetches.
> 2. Call `fetch_into_customized` with `segment.tid().as_u32()`,
>    `segment.start()`, `segment.end()` to load records into the internal
>    `RecordStore<E::Extra>`, passing the customize value through. `compute`
>    and `keep_record` both run inline during push: rejected records are
>    rolled back with zero slab waste and kept records already have their
>    extras populated.
> 3. Fetch the reference sequence for `[segment.start(), segment.end()]`
>    inclusive via the FASTA reader, using the contig name carried by the
>    segment (no header lookup). The fetch MUST use a u64-bounded path so
>    `segment.end() == Pos0::max_value()` does not silently truncate the
>    last base.
> 4. Take the store (now `RecordStore<E::Extra>` with populated extras) and
>    construct a `PileupEngine<E::Extra>` with the fetched reference
>    sequence pre-attached via `set_reference_seq`.
>
> No separate `apply_customize` pass is needed. The returned engine is
> ready for `pileups()` iteration with no further configuration required.
> There MUST NOT be a `pileup(tid, start, end)` overload — callers wanting
> a one-shot region build a `Segment` via `Readers::segments`.

r[unified.fetch_into_customized]
Each format reader (`IndexedBamReader`, `IndexedSamReader`, `IndexedCramReader`, and the format-agnostic `IndexedReader`) MUST expose `fetch_into_customized<E: CustomizeRecordStore>(tid, start, end, store, customize) -> Result<FetchCounts>` in addition to `fetch_into`. The reader MUST forward `customize` to `RecordStore::push_raw`/`push_fields` so its `keep_record` runs at push time. `FetchCounts { fetched, kept }` reports records produced by the reader's built-in overlap/unmapped checks (`fetched`) vs those that also passed `keep_record` (`kept`). Existing `fetch_into` MUST remain a thin wrapper passing `&mut ()` (whose default `keep_record` returns `true`) so its signature and behavior are unchanged.

r[unified.fetch_counts]
`FetchCounts` MUST be `#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]` with `fetched: usize` and `kept: usize` fields. `kept <= fetched` MUST always hold. The struct is re-exported from `crate::reader::FetchCounts` for callers that need to report filter statistics.

> r[unified.readers_accessors]
> `Readers` MUST expose:
>
> - `header() -> &BamHeader` — delegates to the alignment reader's header.
> - `segments(target, opts) -> Result<impl Iterator<Item = Segment>>` — see `r[unified.readers_segments]`. The only way to obtain a `Segment`.
> - `pileup(&Segment) -> Result<PileupEngine<E::Extra>>` — see `r[unified.readers_pileup]`.
> - `fetch_into(tid, start, end, store) -> Result<usize>` — delegates to the alignment reader. Always loads into a `RecordStore<()>` (for custom extras, use `pileup` directly — extras are populated inline at push time).
> - `fasta() -> &IndexedFastaReader` and `fasta_mut() -> &mut IndexedFastaReader` — direct access for callers that need reference sequences independently of the alignment reader (e.g., the call pipeline's segment fetching).
> - `alignment() -> &IndexedReader` and `alignment_mut() -> &mut IndexedReader` — direct access when needed.
> - `customize() -> &E` / `customize_mut() -> &mut E` — direct access to inspect or reset the customize value's state between regions.

r[unified.readers_backward_compat]
`IndexedReader::open(path)` MUST continue to work for BAM, SAM, and CRAM files without a FASTA path. For CRAM the reader follows `r[cram.fasta.optional]`: opening always succeeds, and missing-reference errors only surface at fetch time when a slice actually needs an external reference.

## API surface

See `r[io.non_exhaustive_enums]` and `r[io.minimal_public_api]` in `general.md` for the general rules. For this module, `IndexedReader`, `FormatDetectionError`, and `ReaderError` are the primary enums subject to `r[io.non_exhaustive_enums]`.
