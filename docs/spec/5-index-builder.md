# Index Builder

> **Sources:** [SAM1] §5 "Indexing BAM" — binning algorithm, reg2bin/reg2bins; §5.2 "The BAI index format" — bin/chunk/linear index binary format. [TABIX] — TBI format (magic, column config, same bin/chunk/linear data as BAI). [CSI] — CSI format (parameterised binning, per-bin loffset). htslib `hts.c` `hts_idx_push` / `hts_idx_finish` — single-pass incremental algorithm. See [references.md](99-references.md).

r[index_builder.single_pass]
The index MUST be built incrementally during writing in a single pass. The builder receives (tid, beg, end, virtual_offset) after each record is written. No second pass over the compressed file is required.

r[index_builder.sort_validation]
The builder MUST validate coordinate sort order: tid MUST be non-decreasing, and pos (beg) MUST be non-decreasing within the same tid. Violations MUST produce a typed error.

r[index_builder.offset_convention]
The virtual offset passed to the builder is the position AFTER writing the record (the end of the record). The builder uses the previous call's offset as the start of the current record, following htslib's convention.

r[index_builder.binning]
The builder MUST assign each record to the smallest bin that fully contains its `[beg, end)` interval, using the parameterised `reg2bin(beg, end, min_shift, depth)` formula. Default parameters are min_shift=14, depth=5 (matching BAI/TBI).

r[index_builder.chunk_flush]
Chunks MUST be emitted when the bin changes or the reference sequence changes. A chunk is the virtual offset range `[save_off, last_off)` covering a contiguous run of records in the same bin.

r[index_builder.linear_index]
The builder MUST maintain a linear index: for each `min_shift`-sized window (16 KiB for min_shift=14), store the minimum virtual offset of any record overlapping that window. Only the first record to touch a window slot sets its value.

r[index_builder.linear_backfill]
After all records, unset linear index slots MUST be backfilled right-to-left: each unset slot receives the value of the next valid slot to its right.

r[index_builder.pseudo_bin]
Each reference sequence MUST include a pseudo-bin (ID = `((1 << (depth*3 + 3)) - 1) / 7`, which is 37450 for depth=5) with two chunks: the first stores the virtual offset range `[off_beg, off_end)` of all records for that reference; the second stores mapped and unmapped counts as u64 pairs packed into the chunk's begin/end fields.

r[index_builder.tbi_empty_refs]
References with no records MUST be omitted from TBI output. Only references that received at least one `push()` call are included in n_ref, the sequence names section, and the per-reference bin/chunk/linear data. This matches bcftools behavior.

r[index_builder.tbi_format]
TBI output MUST be BGZF-compressed with magic `TBI\x01`, followed by: n_ref (only refs with data), format (2 for VCF), col_seq (1), col_beg (2), col_end (0), meta (`#` = 35), skip (0), sequence names (only for refs with data); then BAI-compatible bin/chunk/linear data per reference.

r[index_builder.csi_format]
CSI output MUST be BGZF-compressed with magic `CSI\x01`, followed by: min_shift, depth, l_aux, aux data; then per-reference bin data with per-bin loffset (computed from linear index), no separate linear index array.
