# Tabix and CSI Index

Tabix (`.tbi`) and CSI (`.csi`) are generic indexes for bgzf-compressed, coordinate-sorted, TAB-delimited files. They enable region-based random access via BGZF virtual offsets. seqair needs these for bgzf-compressed SAM files.

> **Sources:** [TABIX] — full TBI format specification (magic, header, bins, chunks, linear index, n_no_coor, reg2bin/reg2bins). [CSI] — full CSI format specification (magic, min_shift, depth, l_aux, per-bin loffset, parameterised reg2bin/reg2bins). [SAM1] §5.3 — reg2bin/reg2bins C code for BAI-compatible binning parameters (min_shift=14, depth=5). See [References](./99-references.md).

## Context

BAM uses BAI indexes, which are a specialized form of the binning scheme. Tabix and CSI generalize this:

- **Tabix** (`.tbi`): fixed parameters (min_shift=14, depth=5), identical binning to BAI. Includes column configuration for the text format.
- **CSI** (`.csi`): configurable min_shift and depth, supporting references longer than 512 Mbp.

Both use the same bin/chunk/virtual-offset query mechanism as BAI.

**Implementation strategy**: Tabix with format=SAM uses identical binning parameters to BAI (min_shift=14, depth=5). The index data after the tabix-specific header is structurally identical to BAI (bins, chunks, linear index). The initial implementation SHOULD reuse the existing `BamIndex` parsing and query code by reading past the tabix header and delegating to the BAI parser. CSI support (parameterized binning, per-bin `loff`, no linear index) requires a separate parser and is deferred.

## File compression

r[tabix.compression]
Tabix (`.tbi`) files are BGZF-compressed. The parser MUST decompress via BGZF before parsing the binary format. Reading raw bytes without decompression produces garbage.

r[csi.compression]
CSI (`.csi`) files are BGZF-compressed. Same decompression requirement as tabix.

## Tabix format

> _[TABIX] — TBI format table: magic, n_ref, format, col_seq, col_beg, col_end, meta, skip, l_nm, names; bin/chunk/interval lists; n_no_coor_

r[tabix.magic]
Tabix files MUST start with magic bytes `TBI\x01` (4 bytes), after BGZF decompression.

> r[tabix.header]
> Tabix files MUST continue with these fields after magic bytes:
>
> - `n_ref` (i32): number of reference sequences
> - `format` (i32): 0=generic, 1=SAM, 2=VCF. For SAM files, the format code is 1 (not 0).
> - `col_seq` (i32): column for sequence name (1-based; 3 for SAM)
> - `col_beg` (i32): column for begin position (1-based; 4 for SAM)
> - `col_end` (i32): column for end position (0 = computed from begin)
> - `meta` (i32): comment/header prefix character (64 = `@` for SAM)
> - `skip` (i32): number of header lines to skip (0 for SAM; headers detected by `meta` prefix)
> - `l_nm` (i32): length of concatenated sequence names
> - `names` (bytes): null-terminated sequence names, concatenated

Note: when `format=1` (SAM), htslib hard-codes the column configuration internally. The stored col_seq/col_beg/col_end values may be ignored by some implementations. The reader SHOULD NOT reject a tabix file based on stored column values when format=1.

r[tabix.bai_reuse]
Since tabix uses identical binning parameters to BAI (min_shift=14, depth=5), the index data after the tabix header MUST be parsed by reusing `BamIndex::from_reader()` (or equivalent). This avoids duplicating the bin/chunk/linear-index parsing code. The tabix header is read and validated first, then the remaining bytes are handed to the BAI parser.

> r[tabix.index_data]
> After the header, the index data MUST follow the same structure as BAI:
>
> - For each reference:
>   - `n_bin` (i32): number of bins
>   - For each bin: `bin_id` (u32), `n_chunk` (i32), chunks as virtual offset pairs
>   - `n_intv` (i32): number of linear index entries
>   - Linear index: virtual offsets at 16 kbp intervals
>
> After all reference data, an optional 8-byte `n_no_coor` field (u64) MAY be present (number of unmapped reads with no coordinate). Its presence is detected by whether 8 more bytes remain.

r[tabix.pseudo_bin]
Bin ID 37450 is a pseudo-bin that stores metadata (virtual offset range of all records for a reference). Its chunks MUST NOT be included in query results. The reader MUST recognize and skip this bin during region queries.

> r[tabix.query]
> Region queries use the same `reg2bins` algorithm as BAI (with min_shift=14, depth=5). The bin offset for level `l` (counting from root, l=0=root) is `((1 << 3*l) - 1) / 7`, and the bin ID for position `p` at level `l` is:
>
> `((1 << 3*l) - 1) / 7 + (p >> (min_shift + 3*(depth - l)))`

For BAI/tabix with min_shift=14, depth=5:

- Level 0 (root): offset=0, shift=29. Bin 0 covers the entire reference.
- Level 5 (leaf): offset=4681, shift=14. Each bin covers 16 kbp.

## Shared query infrastructure

r[index.shared_query]
The existing `BamIndex::query()` returns `Vec<Chunk>` with begin/end virtual offsets. The tabix and CSI parsers MUST produce the same `Vec<Chunk>` output so that `RegionBuf::load()` can be reused unchanged. The binning parameters differ but the output format is identical.

r[index.shared_types]
`Chunk` and `VirtualOffset` types from the BAM index module MUST be reused by tabix/CSI indexes.

## Edge cases

r[index.edge.large_refs]
References longer than 512 Mbp require CSI (BAI/tabix with min_shift=14, depth=5 cannot represent them). The reader MUST support CSI indexes for BAM files too (not just SAM), since some BAM files from plant/polyploid genomes use CSI.

r[index.edge.no_records]
A reference with no records has no bins and no linear index entries. Queries for such references MUST return an empty chunk list (not an error).

r[index.edge.name_mismatch]
Tabix stores sequence names in the index header. These MUST match the names in the SAM header's `@SQ` lines. A mismatch SHOULD produce a warning (the index may be stale). CSI does NOT store sequence names — this check applies only to tabix.
