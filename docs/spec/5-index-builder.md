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
Each reference sequence MUST include a pseudo-bin (ID = `((1 << (depth+1)*3) - 1) / 7 + 1`, which is 37450 for depth=5) with two chunks: the first stores the virtual offset range `[off_beg, off_end)` of all records for that reference; the second stores mapped and unmapped counts as u64 pairs packed into the chunk's begin/end fields.

r[index_builder.tbi_empty_refs]
References with no records MUST be omitted from TBI output. Only references that received at least one `push()` call are included in n_ref, the sequence names section, and the per-reference bin/chunk/linear data. This matches bcftools behavior.

r[index_builder.tbi_format]
TBI output MUST be BGZF-compressed with magic `TBI\x01`, followed by: n_ref (only refs with data), format (2 for VCF), col_seq (1), col_beg (2), col_end (0), meta (`#` = 35), skip (0), sequence names (only for refs with data); then BAI-compatible bin/chunk/linear data per reference.

r[index_builder.csi_format]
CSI output MUST be BGZF-compressed with magic `CSI\x01`, followed by: min_shift, depth, l_aux, aux data; then per-reference bin data with per-bin loffset (computed from linear index), no separate linear index array.

## BAI output

The IndexBuilder already accumulates the same bin/chunk/linear data structures that BAI requires — the internal representation is format-agnostic. BAI, TBI, and CSI differ only in their serialization: BAI is uncompressed with a simpler header; TBI adds BGZF compression and column config; CSI replaces the linear index with per-bin loffsets.

Adding `write_bai()` to the existing IndexBuilder enables BAM index co-production (see `r[bam_writer.index_coproduction]` in [bam-writer.md](6-bam-writer.md)) without duplicating the single-pass state machine. This is the same builder that already produces TBI for VCF and CSI for BCF.

> *[SAM1] §5.2 "The BAI index format for BAM files" — magic, n_ref, bins, chunks, linear index*

r[index_builder.bai_format]
BAI output MUST be uncompressed (not BGZF-wrapped) with magic `BAI\x01` (0x42, 0x41, 0x49, 0x01), followed by: n_ref (i32 LE, always equal to the total number of reference sequences in the BAM header), then per-reference: n_bin (i32 LE), for each bin: bin_id (u32 LE) + n_chunk (i32 LE) + chunks (begin/end as u64 LE virtual offsets), then n_intv (i32 LE) + linear index entries (u64 LE virtual offsets). This is the same bin/chunk/linear data as TBI but without BGZF compression, without the TBI column config header, and with a different magic.

r[index_builder.bai_all_refs]
Unlike TBI (`r[index_builder.tbi_empty_refs]`), BAI output MUST include an entry for every reference sequence in the header, even those with no records. References with no records MUST be written with n_bin=0 and n_intv=0. This matches samtools index behavior and is required because BAI uses tid as a positional index — the reader looks up the Nth entry by seeking past the first N entries, so gaps would cause tid/entry misalignment.

r[index_builder.bai_constructor]
`IndexBuilder::bai(n_refs)` MUST construct a builder configured for BAI output with min_shift=14, depth=5. The `header_end_offset` parameter is not needed for BAI (it is only used for TBI/CSI to record where the header ends). The builder MUST pre-allocate `n_refs` empty `RefIndexBuilder` entries so that `write_bai()` can emit all references, including those that received no records.

r[index_builder.bai_write]
`write_bai<W: Write>(&self, writer: W, n_refs: usize)` MUST serialize the accumulated index data into BAI format. The `n_refs` parameter specifies the total number of reference sequences in the BAM header — this may exceed the number of `RefIndexBuilder` entries if the builder received records for only a subset of references. For reference indices `0..n_refs`, the method MUST emit the builder's data if available, or n_bin=0 and n_intv=0 if not (per `r[index_builder.bai_all_refs]`). Bins within each reference MUST be written in ascending bin_id order (using `BTreeMap` iteration, matching htslib's deterministic output). The method MUST validate that `finish()` was called before writing, returning an error otherwise.
