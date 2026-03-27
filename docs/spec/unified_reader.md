# Unified Reader

seqair reads BAM, bgzf-compressed SAM, and CRAM files through a format-agnostic reader interface that auto-detects the input format and dispatches to the appropriate parser, while presenting a uniform API to the rest of the codebase.

> **Sources:** The unified reader design is seqair-specific. Format magic bytes come from [SAM1] §4.2 (BAM: `BAM\1`), [SAM1] §4.1 (BGZF: `1f 8b`), and [CRAM3] §6 (CRAM: `CRAM`). Sort order detection uses the `SO` tag from [SAM1] §1.3 "The header section" (`@HD` line). See [references.md](references.md).

## Goals

1. **Backward compatibility**: existing code that calls `IndexedBamReader::open()` continues to work unchanged.
2. **Auto-detection**: `IndexedReader::open(path)` inspects the file to determine BAM vs SAM vs CRAM and opens the appropriate reader.
3. **Uniform downstream API**: all three formats populate the same `RecordStore` via `fetch_into()`. The pileup engine, call pipeline, and all downstream code are format-agnostic.
4. **Forkable**: all reader types support `fork()` for thread-safe parallel processing with shared immutable state.

## Format detection

r[unified.detect_format]
`IndexedReader::open(path)` MUST auto-detect the file format by inspecting magic bytes:
- Bytes `1f 8b` (gzip magic) → verify BGZF structure (extra field with `BC` subfield). If not BGZF, reject with error: "file is gzip-compressed but not BGZF; use `bgzip` instead of `gzip`." If BGZF, decompress the first block:
  - Starts with `BAM\x01` → BAM format
  - Starts with `@` (0x40) → bgzf-compressed SAM
  - Otherwise → reject with error naming supported formats
- Bytes `43 52 41 4d` (`CRAM`) → CRAM format
- Byte `@` (0x40) at position 0 → uncompressed SAM. Reject with error: "uncompressed SAM files cannot be indexed; compress with `bgzip file.sam` then index with `samtools index file.sam.gz`."
- Otherwise → return an error listing supported formats (BAM, bgzf-compressed SAM, CRAM)

r[unified.detect_index]
After format detection, the reader MUST locate the appropriate index file:
- BAM: `.bai`, `.bam.bai`, or `.csi` (CSI supports references > 512 Mbp that BAI cannot)
- SAM: `.tbi` or `.csi` (tabix or CSI index)
- CRAM: `.crai` (CRAM index)

r[unified.detect_error]
If the format is detected but no matching index is found, the error MUST name the expected index extension and suggest the tool to create it (`samtools index` for BAM/CRAM, `samtools index` or `tabix` for SAM).

## Reader enum

r[unified.reader_enum]
The unified reader MUST be an enum dispatching to format-specific readers, not a trait object. This avoids dynamic dispatch overhead in the hot path and keeps the type concrete for `fork()`.

```
enum IndexedReader {
    Bam(IndexedBamReader),
    Sam(IndexedSamReader),
    Cram(IndexedCramReader),
}
```

r[unified.reader_api]
The unified reader MUST expose:
- `header() -> &BamHeader` — all formats produce the same header type (target names, lengths, tid lookup). The SAM and CRAM parsers convert their native header representations to `BamHeader`.
- `fetch_into(tid, start, end, store) -> Result<usize>` — populates a `RecordStore` with records overlapping the region, identically to the BAM path.
- `fork() -> Result<Self>` — creates a lightweight copy sharing immutable state.
- `shared() -> &Arc<_>` — access to shared state for `Arc::ptr_eq` testing.

r[unified.fetch_equivalence]
For a given BAM file and its SAM/CRAM representations of the same data, `fetch_into` for the same region MUST produce records with equivalent logical content: same positions, flags, sequences, qualities, and CIGAR operations. Aux tags MUST contain the same set of tag names and values, but tag ordering and integer type codes (e.g., `c` vs `i` for small values) MAY differ between formats because BAM writers choose specific integer widths that SAM text cannot preserve.

## RecordStore integration

r[unified.record_store_push]
Currently `RecordStore::push_raw()` takes raw BAM bytes. For SAM and CRAM, records arrive as parsed fields rather than BAM binary.

The RecordStore MUST provide `push_fields()` (or equivalent) that accepts pre-parsed record fields:
- pos, end_pos, flags, mapq, seq_len, matching_bases, indel_bases (fixed fields)
- qname bytes
- CIGAR as packed BAM-format u32 ops (SAM/CRAM parsers convert to this representation)
- sequence as `&[Base]` (already the enum type stored in the bases slab — no conversion needed)
- quality bytes
- aux tag bytes in BAM binary format (SAM/CRAM parsers serialize tags to this format)

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
BAM fork: shares `Arc<BamShared>` (index + header), opens fresh `File` handle. (Already implemented.)

r[unified.fork_sam]
SAM fork: shares `Arc` holding parsed tabix/CSI index + header, opens fresh `File` handle for BGZF reading.

r[unified.fork_cram]
CRAM fork: shares `Arc` holding CRAI index + header, opens fresh `File` handle for container reading. Each fork MUST have its own FASTA reader (via `fasta_reader.fork()`) for thread-safe reference lookups. Reference caching (`r[cram.perf.reference_caching]`) MUST be per-fork, not shared.

## Readers: alignment + reference bundle

r[unified.readers_struct]
The `Readers` struct bundles an `IndexedReader` (alignment) with an `IndexedFastaReader` (reference) in a single type. This eliminates the need for separate FASTA path passing — CRAM can access the reference it needs for sequence reconstruction, and all formats have uniform open/fork semantics.

```
pub struct Readers {
    alignment: IndexedReader,
    fasta: IndexedFastaReader,
}
```

r[unified.readers_open]
`Readers::open(alignment_path, fasta_path)` MUST auto-detect the alignment format (via `r[unified.detect_format]`), open the appropriate reader, and open the FASTA reader. For CRAM, the fasta_path is passed to `IndexedCramReader::open()` for sequence reconstruction. For BAM/SAM, the FASTA reader is opened but not used by the alignment reader.

r[unified.readers_fork]
`Readers::fork()` MUST fork both the alignment reader and the FASTA reader, returning a new `Readers` with independent I/O handles but shared immutable state (indices, headers). The CRAM fork gets its own FASTA reader via `IndexedFastaReader::fork()`.

r[unified.readers_accessors]
`Readers` MUST expose:
- `header() -> &BamHeader` — delegates to the alignment reader's header.
- `fetch_into(tid, start, end, store) -> Result<usize>` — delegates to the alignment reader.
- `fasta() -> &IndexedFastaReader` and `fasta_mut() -> &mut IndexedFastaReader` — direct access for callers that need reference sequences independently of the alignment reader (e.g., the call pipeline's segment fetching).
- `alignment() -> &IndexedReader` and `alignment_mut() -> &mut IndexedReader` — direct access when needed.

r[unified.readers_backward_compat]
`IndexedReader::open(path)` MUST continue to work for BAM and SAM files without a FASTA path. CRAM detection in `IndexedReader::open()` MUST return an error explaining that CRAM requires a reference and suggesting `Readers::open()` instead. This preserves backward compatibility for code that only needs BAM/SAM.

## API surface

See `r[io.non_exhaustive_enums]` and `r[io.minimal_public_api]` in `general.md` for the general rules. For this module, `IndexedReader`, `FormatDetectionError`, and `ReaderError` are the primary enums subject to `r[io.non_exhaustive_enums]`.
