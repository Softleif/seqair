# BAM Writer

The BAM writer serializes alignment records into the BAM binary format: a BGZF-compressed file containing a header followed by a stream of binary records. This is the write counterpart to the BAM reader ([BAM Reader](./2-bam-1-reader.md)) and consumes `BamRecord` values from the owned record type ([BAM Record Builder](./6-bam-record-builder.md)).

BAM writing is needed in two primary contexts:

1. **BAM rewriting** — the TAPS methylation pipeline reads a BAM file, annotates each record with methylation information (MM/ML or XR/XG/XM tags), optionally rewrites the sequence (T->C, A->G to reverse methylation evidence), and writes a new BAM file. This is a record-by-record transformation where each output record is a modified copy of an input record.

2. **Post-realignment output** — after in-memory realignment (modifying CIGAR and position), modified records are written to a new BAM file. This may also involve sorting, though initial support can require pre-sorted input.

> **Sources:** [SAM1] §4 "The BAM Format Specification" — file structure (magic, header, records); §4.2 "The BAM format" — binary record layout, bin_mq_nl / flag_nc packing, l_text / n_ref / l_name / l_ref header fields; §4.1 "The BGZF compression format" — block compression (implemented by `BgzfWriter`, see [BZGF](./1-2-bgzf.md)). See [References](./99-references.md).

## Background

### Current state

The BAM rewrite pipeline in rastair (`src/bam.rs`, ~2150 lines) uses `rust_htslib::bam::Writer` for output. The writer is configured with a header template (cloned from the input BAM, with an added PG record), compression level (fastest), background threads (3), and output target (file path or stdout). Records are modified in memory using `rust_htslib::bam::Record` methods (`set_seq`, `push_aux`) and written one at a time.

The pipeline is parallelized: rayon workers each process a genomic region, producing modified records. An ordered channel ensures output records maintain coordinate order. A single writer thread consumes the channel and writes to the BAM file.

### Design goals

- **API compatibility**: method names and semantics should be familiar to users of rust-htslib's writer
- **Single-threaded first**: start with single-threaded BGZF compression (using `BgzfWriter` from [BZGF](./1-2-bgzf.md)). Multi-threaded compression is a future optimization
- **Streaming**: records are written one at a time, not buffered in memory. The writer handles BGZF block boundaries transparently
- **Ordered output**: the writer preserves insertion order — records appear in the output file in the order `write()` is called

## Writer

> _[SAM1] §4 — BAM file = magic + header + records, all BGZF-compressed_

r[bam_writer.create_from_path]
`BamWriter::from_path(path, header, build_index)` MUST create a BAM file at the given path (truncating any existing file), write the BAM magic and serialized header, and return a writer ready to accept records. The file MUST be BGZF-compressed using `BgzfWriter`. When `build_index` is true, an `IndexBuilder` configured for BAI MUST be created with `n_refs = header.n_targets()` (see `r[bam_writer.index_coproduction]`).

r[bam_writer.header_eager]
The BAM header (magic + header text + reference sequences) is written eagerly during construction (`from_path` / `from_stdout`). There is no separate `write_header()` method — unlike the VCF writer (`r[vcf_writer.header_first]`), which defers header writing. This is because BAM always requires a header (there is no "headerless BAM"), and writing it eagerly simplifies the state machine (no "header not yet written" error path on `write()`).

r[bam_writer.create_from_stdout]
`BamWriter::from_stdout(header)` MUST write BGZF-compressed BAM to stdout. Index co-production is not supported for stdout output (there is no sidecar path for the `.bai` file). The caller is responsible for determining whether stdout is the intended target (e.g. by checking for `-` or `/dev/stdout` in a CLI argument).

r[bam_writer.write_record]
`write(record: &BamRecord)` MUST serialize the record into the writer's reusable buffer via `BamRecord::to_bam_bytes()` (see `r[bam.owned_record.to_bam_bytes]`), prepend the 4-byte `block_size` (i32 little-endian, equal to the serialized byte count), and write the result to the BGZF stream. If index co-production is enabled, the writer MUST then push the record to the IndexBuilder (see `r[bam_writer.index_record_dispatch]`).

r[bam_writer.insertion_order]
Records MUST appear in the output file in the exact order that `write()` is called. The writer MUST NOT reorder, buffer, or batch records.

r[bam_writer.finish]
`finish()` MUST flush any buffered BGZF data. If index co-production is enabled, it MUST call `index.finish(final_virtual_offset)` after flushing all record data (see `r[bam_writer.index_finish]`). It MUST then write the BGZF EOF marker block (as specified by `r[bgzf.writer.eof_marker]`). The method MUST return `Result<(W, Option<IndexBuilder>), BamWriteError>` — the inner writer `W` (matching the VCF writer's `r[vcf_writer.finish]` pattern) and the finished IndexBuilder if co-production was enabled. Dropping the writer without calling `finish()` SHOULD flush on a best-effort basis (see `r[bgzf.writer.finish]`).

r[bam_writer.compression_level]
The compression level MUST be configurable as a constructor parameter (not a post-construction setter), since the header is written during construction and its BGZF blocks use whatever level is set at that point. The default SHOULD be level 6 (matching htslib). Level 1 (fastest) is commonly used in pipelines where downstream tools re-compress. The constructor signatures are: `from_path(path, header, build_index)` uses the default level; `from_path_with_level(path, header, build_index, level)` accepts an explicit level.

## BGZF block management

BAI virtual offsets encode both the compressed block offset and the position within the uncompressed block. When a record starts near the end of one block and continues into the next, the index still points to the correct position — but flushing proactively keeps most records block-aligned, which means random access typically decompresses a single block rather than two.

r[bam_writer.flush_before_record]
Before writing each record, the writer MUST call `BgzfWriter::flush_if_needed(record_size)` where `record_size` includes the 4-byte `block_size` prefix plus the serialized record bytes. This ensures the record starts at a fresh BGZF block when the current block cannot fit it, improving index seek granularity.

r[bam_writer.record_spanning_blocks]
Records larger than the BGZF uncompressed block size (64 KiB) MUST be written correctly across multiple blocks. The writer MUST NOT reject records that span block boundaries — `BgzfWriter::write_all` handles multi-block writes transparently. In practice, most BAM records are 200–500 bytes (150bp read + CIGAR + qual + a few aux tags), so block-spanning is rare.

r[bam_writer.record_size_limit]
The BAM `block_size` field is an i32, so a single record's serialized size (excluding the 4-byte prefix) MUST NOT exceed `i32::MAX` (2,147,483,647 bytes). Additionally, for round-trip compatibility with seqair's BAM reader (which enforces `r[bam.record.max_size]`), the writer SHOULD reject records exceeding 2 MiB and MUST reject records exceeding `i32::MAX`. The 2 MiB soft limit matches htslib's default and covers all practical BAM records; records approaching this size almost certainly indicate a bug.

## Index co-production

Without co-production, generating a BAI requires a second pass over the entire BAM file (`samtools index`). For whole-genome BAM rewriting (~50 GB output), this adds minutes of wall time and doubles I/O. Co-producing the index during writing eliminates the second pass entirely — the IndexBuilder accumulates bin/chunk/linear data as records stream through, using virtual offsets captured from the BgzfWriter.

This follows the same pattern established for VCF/BCF writing, where TBI/CSI indexes are co-produced during writing (see [Index Builder](./5-index-builder.md)).

> _[SAM1] §5.2 — BAI index format. The index builder (`r[index_builder.single_pass]`) supports BAI output via `r[index_builder.bai_format]`._

r[bam_writer.index_coproduction]
The writer MUST support co-producing a BAI index during writing. When enabled (via the `build_index` constructor parameter), the writer MUST hold an `IndexBuilder` configured for BAI via `IndexBuilder::bai(header.n_targets())` (see `r[index_builder.bai_constructor]`). After each record is written to the BGZF stream, the writer MUST dispatch the record to the IndexBuilder according to `r[bam_writer.index_record_dispatch]`.

> r[bam_writer.index_record_dispatch]
> When index co-production is enabled, the writer MUST classify each record into one of three cases and handle it accordingly:
>
> 1. **Mapped** (flags & 0x4 == 0, ref_id >= 0): push to the IndexBuilder with `beg` = 0-based pos, `end` = 0-based exclusive end position (pos + reference-consuming CIGAR ops), `virtual_offset` = the BGZF virtual offset after the record was written.
> 2. **Placed unmapped** (flags & 0x4 != 0, ref_id >= 0): push to the IndexBuilder with `beg = end = pos`. These are reads that are unmapped themselves but placed at a position — typically their mate's position in coordinate-sorted BAM files. They must be indexed so that region queries find them (e.g., `samtools view chr1:1000-2000` returns placed-unmapped reads in that range).
> 3. **Fully unmapped** (ref_id == -1): MUST NOT be pushed to the IndexBuilder. These reads have no genomic position and are not retrievable by region query. They appear at the end of coordinate-sorted BAM files.
>
> A mapped record with ref_id == -1 is structurally invalid (mapped flag set but no reference assigned). The writer MUST return a `BamWriteError` for such records when index co-production is enabled, since the IndexBuilder cannot index a record with no valid tid.

r[bam_writer.index_finish]
When `finish()` is called and index co-production is enabled, the writer MUST call `index.finish(final_virtual_offset)` after flushing all record data and before writing the BGZF EOF block. The finished IndexBuilder is returned via `finish()` (see `r[bam_writer.finish]`). The caller is responsible for writing the `.bai` file via `index.write_bai()`. The writer itself MUST NOT write the index file — the caller controls the output path. This separation is important for atomic-rename workflows (write index to a temp file, then rename) and because stdout output has no sidecar path.

r[bam_writer.index_sort_order]
When index co-production is enabled, the writer relies on the IndexBuilder's sort validation (`r[index_builder.sort_validation]`). If records are written out of coordinate order, the IndexBuilder MUST return an error, which the writer MUST propagate as a `BamWriteError`. The writer itself does not enforce sort order — it is the caller's responsibility to write records in coordinate-sorted order when indexing is desired. In rastair's pipeline, the ordered channel from rayon workers already guarantees this.

## Header serialization

> _[SAM1] §4.2 — BAM header: magic bytes 0x42 0x41 0x4D 0x01 (`BAM\1`), l_text (i32 LE), header text, n_ref (i32 LE), then per-reference: l_name (i32 LE) + name + NUL + l_ref (i32 LE)_

r[bam_writer.magic]
The writer MUST begin the BAM file with the 4-byte magic `0x42 0x41 0x4D 0x01` (ASCII `BAM` followed by byte value 1). This is the same magic validated by `r[bam.header.magic]` in [BAM Header](./2-bam-3-1-header.md).

r[bam_writer.header_text]
After the magic, the writer MUST write the header text length as an i32 (little-endian) followed by the header text bytes. The header text MUST NOT include a NUL terminator (the length field is sufficient per [SAM1] §4.2). The header text contains the `@HD`, `@SQ`, `@RG`, `@PG`, and `@CO` lines from the `BamHeader`.

r[bam_writer.header_references]
After the header text, the writer MUST write the number of reference sequences as an i32 (little-endian), followed by each reference's name length as an i32 (little-endian, including NUL byte), name bytes (NUL-terminated), and sequence length as an i32 (little-endian). These MUST match the `BamHeader`'s target list.

## Header construction

The following rules extend `BamHeader` (defined in [BAM Header](./2-bam-3-1-header.md)) with construction and mutation capabilities needed for BAM writing.

r[bam_writer.header_from_template]
`BamHeader` MUST support creating a copy from an existing header that clones all reference sequences, read groups, and other header lines. The caller can then add PG records or modify lines before passing it to the writer.

r[bam_writer.header_add_pg]
`BamHeader` MUST support adding `@PG` header records with fields: `ID` (required), `PN` (program name), `VN` (version), `CL` (command line), `DS` (description), and `PP` (previous program ID). The `PP` field links PG records into a provenance chain per [SAM1] §1.3 — when adding a new PG record, `PP` SHOULD be set to the `ID` of the last existing PG record so that tools like Picard's `ValidateSamFile` can trace the processing history. This is used by downstream tools (e.g. rastair) to record processing operations in the output BAM header.

r[bam_writer.header_text_generation]
`BamHeader` MUST support generating SAM-format header text from its structured data. The generated text MUST include `@HD` (with `VN` and `SO` sort order if known), `@SQ` lines for all references, and any `@RG`, `@PG`, and `@CO` lines that have been added. The `@HD SO:coordinate` tag MUST be preserved from the source header when using `from_template`, and MUST be settable when building from scratch.

## Error handling

r[bam_writer.error_type]
Write errors MUST be represented by a dedicated `BamWriteError` enum with typed variants (per `r[io.errors.typed_variants]`).

Examples:

- I/O errors from the underlying stream
- BGZF compression failures
- Record serialization errors (e.g. qname too long, sequence/quality length mismatch, CIGAR op count overflow, pos/ref_id overflow)
- Index errors (unsorted input, invalid record for indexing — from `IndexError` via `#[from]`)
- Writer poisoned (subsequent writes after a previous failure)

r[bam_writer.error_poisoning]
If writing a record fails (I/O error, compression error), the writer MUST return the error immediately and enter a poisoned state. All subsequent `write()` calls on a poisoned writer MUST return a `Poisoned` error without attempting I/O. The caller is responsible for discarding the partial output file. This is necessary because a multi-block record may have partially flushed BGZF blocks before the error, making "discard the block" infeasible.

## Performance

r[bam_writer.reuse_buffers]
The writer MUST reuse a single `Vec<u8>` serialization buffer across `write()` calls. The buffer is cleared (not deallocated) before each record to avoid per-record allocation. `BamRecord::to_bam_bytes()` appends into this buffer (see `r[bam.owned_record.to_bam_bytes]`).

r[bam_writer.multithreaded_compression]
The writer SHOULD support multi-threaded BGZF compression as a future optimization. A `set_threads(n)` method SHOULD configure background compression workers, where `n=0` (default) means synchronous compression in the calling thread. This matches htslib's `hts_set_threads` pattern and is important for write-heavy pipelines like BAM rewriting at whole-genome scale.

## Testing

r[bam_writer.test_roundtrip]
A BAM file written by `BamWriter` MUST be readable by seqair's BAM reader, by `samtools view`, and by `noodles::bam::Reader`. All record fields (pos, flags, mapq, cigar, seq, qual, aux tags, qname) MUST round-trip exactly. Noodles validation is important because it is a pure-Rust implementation with independent parsing logic — bugs that both seqair and samtools miss due to shared htslib heritage may be caught by noodles.

r[bam_writer.test_samtools_quickcheck]
Output BAM files MUST pass `samtools quickcheck` (validates BGZF structure, EOF block, and basic header integrity).

r[bam_writer.test_header_roundtrip]
Headers written by `BamWriter` MUST be readable by seqair's BAM reader. Reference sequence names and lengths MUST match. PG records (including PP chain) and sort order MUST be preserved.

r[bam_writer.test_empty_file]
Writing a BAM file with a header but zero records MUST produce a valid BAM file (magic + header + EOF block) that `samtools view` accepts.

r[bam_writer.test_large_record]
Records near the 2 MiB size limit (see `r[bam_writer.record_size_limit]`) MUST serialize and round-trip correctly. Records exceeding the limit MUST be rejected with a clear error.

r[bam_writer.test_cross_validation]
For the initial release, a comparison test MUST write the same records with both seqair's `BamWriter` and `rust_htslib::bam::Writer`, then compare the outputs field-by-field (record contents must match; BGZF block boundaries may differ).

r[bam_writer.test_index_roundtrip]
A BAM file written with index co-production MUST produce a BAI that passes `samtools quickcheck` on the BAM+BAI pair. The BAI MUST be parseable by seqair's `BamIndex::from_path()` and MUST return correct chunks for region queries — verified by comparing query results against `samtools view -c` counts for the same regions.

r[bam_writer.test_index_unsorted]
Writing records out of coordinate order with index co-production enabled MUST return an error from the IndexBuilder. The error MUST be propagated through `write()` as a `BamWriteError`.

r[bam_writer.test_fuzz]
A fuzz target SHOULD construct random `BamRecord` values (with arbitrary field combinations including edge-case lengths, empty CIGARs, maximum-length qnames, and varied aux tags), serialize them via `BamWriter`, and read them back via seqair's BAM reader. This catches encoder/decoder divergence. While the writer receives structured input (not raw untrusted bytes), a read-modify-write pipeline operating on a malicious source BAM can produce adversarial field combinations that the writer must handle without panicking.
