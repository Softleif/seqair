# BAM Writer

The BAM writer serializes alignment records into the BAM binary format: a BGZF-compressed file containing a header followed by a stream of binary records. This is the write counterpart to the BAM reader ([bam_reader.md](2-bam-1-reader.md)) and consumes `BamRecord` values from the owned record type ([bam_record_builder.md](6-bam-record-builder.md)).

BAM writing is needed in two primary contexts:

1. **BAM rewriting** — the TAPS methylation pipeline reads a BAM file, annotates each record with methylation information (MM/ML or XR/XG/XM tags), optionally rewrites the sequence (T->C, A->G to reverse methylation evidence), and writes a new BAM file. This is a record-by-record transformation where each output record is a modified copy of an input record.

2. **Post-realignment output** — after in-memory realignment (modifying CIGAR and position), modified records are written to a new BAM file. This may also involve sorting, though initial support can require pre-sorted input.

> **Sources:** [SAM1] §4 "The BAM Format Specification" — file structure (magic, header, records); §4.2 "The BAM format" — binary record layout, bin_mq_nl / flag_nc packing, l_text / n_ref / l_name / l_ref header fields; §4.1 "The BGZF compression format" — block compression (implemented by `BgzfWriter`, see [bgzf.md](1-2-bgzf.md)). See [references.md](99-references.md).

## Background

### Current state

The BAM rewrite pipeline in rastair (`src/bam.rs`, ~2150 lines) uses `rust_htslib::bam::Writer` for output. The writer is configured with a header template (cloned from the input BAM, with an added PG record), compression level (fastest), background threads (3), and output target (file path or stdout). Records are modified in memory using `rust_htslib::bam::Record` methods (`set_seq`, `push_aux`) and written one at a time.

The pipeline is parallelized: rayon workers each process a genomic region, producing modified records. An ordered channel ensures output records maintain coordinate order. A single writer thread consumes the channel and writes to the BAM file.

### Design goals

- **API compatibility**: method names and semantics should be familiar to users of rust-htslib's writer
- **Single-threaded first**: start with single-threaded BGZF compression (using `BgzfWriter` from [bgzf.md](1-2-bgzf.md)). Multi-threaded compression is a future optimization
- **Streaming**: records are written one at a time, not buffered in memory. The writer handles BGZF block boundaries transparently
- **Ordered output**: the writer preserves insertion order — records appear in the output file in the order `write()` is called

## Writer

> *[SAM1] §4 — BAM file = magic + header + records, all BGZF-compressed*

r[bam_writer.create_from_path]
`BamWriter::from_path(path, header)` MUST create a BAM file at the given path (truncating any existing file), write the BAM magic and serialized header, and return a writer ready to accept records. The file MUST be BGZF-compressed using `BgzfWriter`.

r[bam_writer.create_from_stdout]
`BamWriter::from_stdout(header)` MUST write BGZF-compressed BAM to stdout. The caller is responsible for determining whether stdout is the intended target (e.g. by checking for `-` or `/dev/stdout` in a CLI argument).

r[bam_writer.write_record]
`write(record: &BamRecord)` MUST serialize the record via `BamRecord::to_bam_bytes()` (see `r[bam.owned_record.to_bam_bytes]` in [bam_record_builder.md](6-bam-record-builder.md)), prepend the 4-byte `block_size` (i32 little-endian, equal to the serialized byte count), and write the result to the BGZF stream.

r[bam_writer.insertion_order]
Records MUST appear in the output file in the exact order that `write()` is called. The writer MUST NOT reorder, buffer, or batch records.

r[bam_writer.finish]
`finish()` MUST flush any buffered BGZF data and write the BGZF EOF marker block (as specified by `r[bgzf.writer.eof_marker]` in [bgzf.md](1-2-bgzf.md)), then close the output stream. Dropping the writer without calling `finish()` SHOULD flush on a best-effort basis (see `r[bgzf.writer.finish]`).

r[bam_writer.compression_level]
The writer MUST support configurable compression levels via `set_compression_level()`. The default SHOULD be level 6 (matching htslib). Level 1 (fastest) is commonly used in pipelines where downstream tools re-compress.

## Header serialization

> *[SAM1] §4.2 — BAM header: magic bytes 0x42 0x41 0x4D 0x01 (`BAM\1`), l_text (i32 LE), header text, n_ref (i32 LE), then per-reference: l_name (i32 LE) + name + NUL + l_ref (i32 LE)*

r[bam_writer.magic]
The writer MUST begin the BAM file with the 4-byte magic `0x42 0x41 0x4D 0x01` (ASCII `BAM` followed by byte value 1). This is the same magic validated by `r[bam.header.magic]` in [header.md](2-bam-3-1-header.md).

r[bam_writer.header_text]
After the magic, the writer MUST write the header text length as an i32 (little-endian) followed by the header text bytes. The header text MUST NOT include a NUL terminator (the length field is sufficient per [SAM1] §4.2). The header text contains the `@HD`, `@SQ`, `@RG`, `@PG`, and `@CO` lines from the `BamHeader`.

r[bam_writer.header_references]
After the header text, the writer MUST write the number of reference sequences as an i32 (little-endian), followed by each reference's name length as an i32 (little-endian, including NUL byte), name bytes (NUL-terminated), and sequence length as an i32 (little-endian). These MUST match the `BamHeader`'s target list.

## Header construction

The following rules extend `BamHeader` (defined in [header.md](2-bam-3-1-header.md)) with construction and mutation capabilities needed for BAM writing.

r[bam_writer.header_from_template]
`BamHeader` MUST support creating a copy from an existing header that clones all reference sequences, read groups, and other header lines. The caller can then add PG records or modify lines before passing it to the writer.

r[bam_writer.header_add_pg]
`BamHeader` MUST support adding `@PG` header records with fields: `ID` (required), `PN` (program name), `VN` (version), `CL` (command line), `DS` (description). This is used by downstream tools (e.g. rastair) to record processing operations in the output BAM header.

r[bam_writer.header_text_generation]
`BamHeader` MUST support generating SAM-format header text from its structured data. The generated text MUST include `@HD` (with `VN` and `SO` sort order if known), `@SQ` lines for all references, and any `@RG`, `@PG`, and `@CO` lines that have been added. The `@HD SO:coordinate` tag MUST be preserved from the source header when using `from_template`, and MUST be settable when building from scratch.

## Error handling

r[bam_writer.error_type]
Write errors MUST be represented by a dedicated `BamWriteError` enum with typed variants (per `r[io.errors.typed_variants]` in [general.md](0-general.md)):
- I/O errors from the underlying stream
- BGZF compression failures
- Record serialization errors (e.g. qname too long, sequence/quality length mismatch, CIGAR op count overflow)

r[bam_writer.error_poisoning]
If writing a record fails (I/O error, compression error), the writer MUST return the error immediately and enter a poisoned state. All subsequent `write()` calls on a poisoned writer MUST return an error without attempting I/O. The caller is responsible for discarding the partial output file. This is necessary because a multi-block record may have partially flushed BGZF blocks before the error, making "discard the block" infeasible.

## Performance

r[bam_writer.reuse_buffers]
The writer MUST reuse serialization buffers across `write()` calls. A single `Vec<u8>` for record serialization (cleared and reused each call) avoids per-record allocation.

r[bam_writer.multithreaded_compression]
The writer SHOULD support multi-threaded BGZF compression as a future optimization. A `set_threads(n)` method SHOULD configure background compression workers, where `n=0` (default) means synchronous compression in the calling thread. This matches htslib's `hts_set_threads` pattern and is important for write-heavy pipelines like BAM rewriting at whole-genome scale.

## Testing

r[bam_writer.test_roundtrip]
A BAM file written by `BamWriter` MUST be readable by seqair's BAM reader and by `samtools view`. All record fields (pos, flags, mapq, cigar, seq, qual, aux tags, qname) MUST round-trip exactly.

r[bam_writer.test_samtools_quickcheck]
Output BAM files MUST pass `samtools quickcheck` (validates BGZF structure, EOF block, and basic header integrity).

r[bam_writer.test_header_roundtrip]
Headers written by `BamWriter` MUST be readable by seqair's BAM reader. Reference sequence names and lengths MUST match. PG records and sort order MUST be preserved.

r[bam_writer.test_empty_file]
Writing a BAM file with a header but zero records MUST produce a valid BAM file (magic + header + EOF block) that `samtools view` accepts.

r[bam_writer.test_large_record]
Records near the BAM size limit (~2 MiB) MUST serialize and round-trip correctly. Records exceeding the limit MUST be rejected with a clear error.

r[bam_writer.test_cross_validation]
For the initial release, a comparison test MUST write the same records with both seqair's `BamWriter` and `rust_htslib::bam::Writer`, then compare the outputs field-by-field (record contents must match; BGZF block boundaries may differ).
