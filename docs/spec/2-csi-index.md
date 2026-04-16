# CSI Index (Coordinate-Sorted Index)

CSI is a generalization of the BAI index format that supports larger coordinate spaces and configurable bin sizes. While BAI uses fixed parameters (min_shift=14, depth=5) limiting coordinates to 2^29 (~512 Mbp), CSI stores these parameters in the index header, enabling support for arbitrarily large genomes.

CSI is required for organisms with chromosomes exceeding 512 Mbp (e.g., wheat, axolotl) and for positions larger than ~2 billion bases. The same `reg2bin`/`reg2bins` algorithms are used, but parameterized by `min_shift` and `depth` read from the file header.

> **Sources:** [CSI] CSIv1 specification — magic, header, parameterized binning. [SAM1] section 5 — BAI binning scheme (CSI generalizes this). See [References](./99-references.md).

## File format

> _[CSI] — binary layout_

r[csi.magic]
A CSI file MUST begin with the magic bytes `CSI\1` (0x43, 0x53, 0x49, 0x01). The reader MUST reject files that do not match.

r[csi.header]
After the magic, the header contains three `int32_t` fields: `min_shift` (number of bits for the minimum interval, default 14), `depth` (depth of the binning index, default 5), and `l_aux` (length of auxiliary data in bytes). The reader MUST read all three to parameterize the binning scheme.

r[csi.aux_data]
After the header fields, `l_aux` bytes of auxiliary data follow. For tabix-style CSI files, this contains the tabix metadata (format, col_seq, col_beg, col_end, meta character, skip lines, and the sequence name dictionary). The reader MUST parse the auxiliary data when present.

r[csi.aux_tabix]
When `l_aux > 0`, the auxiliary data block uses the same layout as a tabix header (see `r[tabix.header]`): `int32_t format`, `int32_t col_seq`, `int32_t col_beg`, `int32_t col_end`, `int32_t meta` (meta character), `int32_t skip` (header lines to skip), `int32_t l_nm` (length of concatenated sequence names), followed by `l_nm` bytes of NUL-terminated sequence names. This metadata enables CSI files to serve as tabix indexes for VCF/BED/GFF files.

r[csi.n_ref]
After the auxiliary data, an `int32_t` gives the number of reference sequences (`n_ref`). This MUST match the number of references in the associated BAM/VCF/etc. file's header.

## Per-reference index

r[csi.bins]
For each reference sequence, the index stores `n_bin` (int32_t) followed by that many bin entries. Each bin entry contains:
- `bin` (uint32_t): the bin number
- `loffset` (uint64_t): virtual file offset of the first overlapping record (replaces BAI's separate linear index)
- `n_chunk` (int32_t): number of chunks
- For each chunk: `chunk_beg` (uint64_t) and `chunk_end` (uint64_t) — virtual file offsets

r[csi.no_linear_index]
CSI does NOT have a separate linear index array (unlike BAI which stores one `ioffset` per 16 KiB window per reference). Instead, linear offset information is embedded per-bin via the `loffset` field. This simplifies the format at the cost of slightly less granular skip optimization.

r[csi.loffset]
Each bin's `loffset` records the virtual file offset of the first record overlapping the bin's genomic range. During a query, the reader MUST use `loffset` to skip chunks whose `chunk_end` is before the relevant `loffset`, matching the role of BAI's linear index.

## Binning scheme

r[csi.reg2bin]
The bin for a region `[beg, end)` MUST be computed using the parameterized formula from the CSI spec:
```
reg2bin(beg, end, min_shift, depth):
    s = min_shift; t = ((1 << depth*3) - 1) / 7
    for l = depth down to 1:
        if beg >> s == (end-1) >> s: return t + (beg >> s)
        s += 3; t -= 1 << l*3
    return 0
```
With BAI-compatible parameters (min_shift=14, depth=5), this produces identical bin numbers to BAI's `reg2bin`.

r[csi.reg2bins]
The list of bins overlapping a region `[beg, end)` MUST be computed using the parameterized `reg2bins` from the CSI spec. The implementation MUST iterate all levels from 0 to `depth`, computing candidate bins at each level.

r[csi.bin_limit]
The maximum valid bin number (`bin_limit`) is `((1 << (depth+1)*3) - 1) / 7`. For BAI-compatible parameters (depth=5) this is 37,449. The maximum addressable position is `2^(min_shift + depth*3)` — for BAI-compatible parameters this is 2^29 = 536,870,912. Bin numbers outside `[0, bin_limit)` MUST be treated as errors (except for pseudo-bins, see below).

## Pseudo-bins

r[csi.pseudo_bin]
CSI files MAY contain metadata pseudo-bins for each reference sequence, using bin number `bin_limit(min_shift, depth) + 1`. These pseudo-bins contain the same metadata as BAI pseudo-bins (mapped/unmapped counts, file span). The reader MUST recognize and handle pseudo-bins without treating them as regular bins.

## Unmapped reads

r[csi.n_no_coor]
After all per-reference indices, the file MAY contain an optional `uint64_t` giving the count of unmapped unplaced reads (RNAME `*`). The reader MUST handle both the presence and absence of this field (check for EOF).

## Region queries

r[csi.query]
Given a region (tid, start, end), the reader MUST:
1. Compute the list of overlapping bins using `reg2bins(start, end, min_shift, depth)`
2. For each overlapping bin that exists in the index, collect all chunks
3. Use `loffset` values to skip chunks that cannot contain relevant records
4. Return the merged list of chunks as virtual offset pairs

r[csi.query_interface]
The CSI reader MUST expose the same query interface as the BAI reader. Callers MUST NOT need to know whether the underlying index is BAI or CSI. The `IndexedBamReader` (and `IndexedReader`) MUST transparently support both formats.

## Unified index type

r[csi.unified_enum]
The implementation MUST provide a unified index type (e.g., `enum AlignmentIndex { Bai(BaiIndex), Csi(CsiIndex) }`) that dispatches queries transparently. All query methods (`query`, `query_split`) MUST be available on the unified type. The `IndexedBamReader` MUST accept either format without the caller specifying which.

r[csi.query_split]
The CSI reader SHOULD support `query_split` (separating distant from nearby chunks), parameterized by `depth`. The split threshold SHOULD be `depth - 3` — bins at levels 0 through `depth - 3` are "distant" (suitable for caching per-chromosome), bins at levels `depth - 2` through `depth` are "nearby" (loaded per-query). With BAI-compatible parameters (depth=5), this produces the same L0-L2 / L3-L5 split as BAI. This optimization is significant for large BAM files on network storage, where re-reading distant chunks (spanning gigabytes) for every region query wastes I/O.

## Format detection

r[csi.detect]
When opening an alignment file with an index, the reader MUST check for both `.csi` and `.bai` extensions. Following htslib convention, CSI is preferred: try `{file}.csi` first, then `{file_without_ext}.csi`, then `{file}.bai`, then `{file_without_ext}.bai`. The magic bytes MUST be verified regardless of file extension.

## Index writing

r[csi.write]
The `IndexBuilder` MUST support writing CSI indexes in addition to BAI. The builder MUST accept `min_shift` and `depth` parameters at construction time. When the maximum position in the data exceeds 2^29, the builder SHOULD automatically select CSI output (with appropriate `min_shift`/`depth`) instead of BAI.

r[csi.write_format]
The CSI writer MUST emit the full binary format: magic, `min_shift`, `depth`, `l_aux`, auxiliary data (if any), `n_ref`, then per-reference bins with `loffset` and chunks, and the optional `n_no_coor` trailer.

r[csi.write_loffset]
When building a CSI index, the builder MUST track the `loffset` for each bin — the minimum virtual offset of any record whose alignment start falls within the bin's genomic range. This is computed during the single-pass index construction.

r[csi.write_tabix_aux]
When building a CSI index for tabix (VCF/BED/GFF), the builder MUST write the tabix metadata into the `aux` block. The `format`, `col_seq`, `col_beg`, `col_end`, `meta`, `skip`, and sequence name dictionary MUST be serialized in the standard tabix header layout.

## Alloc limits

r[csi.alloc_limits]
The CSI reader MUST enforce the same allocation limits as the BAI reader (see `r[io.fuzz.alloc_limits]`): `n_ref` MUST NOT exceed 100,000, `n_bin` per reference MUST NOT exceed 100,000, `n_chunk` per bin MUST NOT exceed 1,000,000. These limits prevent denial-of-service from malformed index files.

## Error handling

r[csi.errors]
All CSI-specific errors MUST be part of a typed error enum (extending or paralleling `BaiError`). Errors MUST contain context (file path, field values) sufficient to diagnose the problem. The implementation MUST NOT use `panic!`, `unwrap()`, or string-based errors.
