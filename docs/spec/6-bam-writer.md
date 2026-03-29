# BAM Writer

The BAM writer serializes alignment records into the BAM binary format: a BGZF-compressed file containing a header followed by a stream of binary records. This is the write counterpart to the BAM reader ([bam_reader.md](2-bam-1-reader.md)) and consumes `BamRecord` values from the owned record type ([bam_record_builder.md](6-bam-record-builder.md)).

BAM writing is needed in two primary contexts:

1. **BAM rewriting** — the TAPS methylation pipeline reads a BAM file, annotates each record with methylation information (MM/ML or XR/XG/XM tags), optionally rewrites the sequence (T->C, A->G to reverse methylation evidence), and writes a new BAM file. This is a record-by-record transformation where each output record is a modified copy of an input record.

2. **Post-realignment output** — after in-memory realignment (modifying CIGAR and position), modified records are written to a new BAM file. This may also involve sorting, though initial support can require pre-sorted input.

> **Sources:** [SAM1] §4 "The BAM Format Specification" — file structure (magic, header, records); §4.2 "The BAM format" — binary record layout; §4.1 "The BGZF compression format" — block compression (implemented by `BgzfWriter`, see [bgzf.md](1-2-bgzf.md)). See [references.md](99-references.md).

## Background

### Current state

The BAM rewrite pipeline in rastair (`src/bam.rs`, ~2150 lines) uses `rust_htslib::bam::Writer` for output. The writer is configured with a header template (cloned from the input BAM, with an added PG record), compression level (fastest), background threads (3), and output target (file path or stdout). Records are modified in memory using `rust_htslib::bam::Record` methods (`set_seq`, `push_aux`) and written one at a time.

The pipeline is parallelized: rayon workers each process a genomic region, producing modified records. An ordered channel ensures output records maintain coordinate order. A single writer thread consumes the channel and writes to the BAM file.

### Design goals

- **API compatibility**: method names and semantics should be familiar to users of rust-htslib's writer
- **Single-threaded first**: start with single-threaded BGZF compression (using `BgzfWriter` from [bgzf.md](1-2-bgzf.md)). Multi-threaded compression is a future optimization
- **Streaming**: records are written one at a time, not buffered in memory. The writer handles BGZF block boundaries transparently
- **Header construction**: support both cloning from an existing `BamHeader` and building from scratch

## Writer

> *[SAM1] §4 — BAM file = magic + header + records, all BGZF-compressed*

r[bam_writer.create_from_path]
`BamWriter::from_path(path, header)` MUST create a BAM file at the given path, write the BAM magic and serialized header, and return a writer ready to accept records. The file MUST be BGZF-compressed using `BgzfWriter`.

r[bam_writer.create_from_stdout]
`BamWriter::from_stdout(header)` MUST write BGZF-compressed BAM to stdout. This is used when the output file is `-` or `/dev/stdout`.

r[bam_writer.write_record]
`write(record: &BamRecord)` MUST serialize the record using `BamRecord::to_bam_bytes()` (see [bam_record_builder.md](6-bam-record-builder.md)) and write the resulting bytes to the BGZF stream. The 4-byte block_size prefix (total record bytes minus 4) MUST be prepended.

r[bam_writer.finish]
`finish()` MUST flush any buffered data in the BGZF writer, write the BGZF EOF marker block, and close the output stream. Dropping the writer without calling `finish()` SHOULD flush on a best-effort basis (see [bgzf.md](1-2-bgzf.md) `bgzf.writer.finish`).

r[bam_writer.compression_level]
The writer MUST support configurable compression levels via `set_compression_level()`. The default SHOULD be level 6 (matching htslib). Level 1 (fastest) is commonly used in pipelines where downstream tools re-compress.

## Header serialization

> *[SAM1] §4.2 — BAM header: magic `BAM\1`, l_text (i32), header text, n_ref (i32), then per-reference: l_name (i32) + name + NUL + l_ref (i32)*

r[bam_writer.magic]
The writer MUST begin the BAM file with the 4-byte magic `BAM\1`.

r[bam_writer.header_text]
After the magic, the writer MUST write the header text length (i32, little-endian) followed by the header text bytes. The header text contains the `@HD`, `@SQ`, `@RG`, `@PG`, and `@CO` lines from the `BamHeader`.

r[bam_writer.header_references]
After the header text, the writer MUST write the number of reference sequences (i32), followed by each reference's name length (i32, including NUL), name bytes (NUL-terminated), and sequence length (i32). These MUST match the `BamHeader`'s target list.

## Header construction

> *Building headers for output BAM files requires cloning the input header and adding program records.*

r[bam_writer.header_from_template]
`BamHeader::from_template(other)` MUST create a new header that clones all reference sequences, read groups, and other header lines from the source. The caller can then add PG records or modify lines before passing it to the writer.

r[bam_writer.header_add_pg]
`BamHeader` MUST support adding `@PG` header records with fields: `ID` (required), `PN` (program name), `VN` (version), `CL` (command line), `DS` (description). This is used by rastair to record the rewrite operation in the output BAM header.

r[bam_writer.header_text_generation]
`BamHeader` MUST support generating the SAM-format header text (`@HD`, `@SQ`, `@RG`, `@PG`, `@CO` lines) from its structured data. This text is what gets written in the BAM header section.

## Error handling

r[bam_writer.error_type]
Write errors MUST be represented by a dedicated `BamWriteError` enum with typed variants:
- I/O errors from the underlying stream
- BGZF compression failures
- Record serialization errors (e.g. qname too long, sequence/quality length mismatch)

r[bam_writer.no_silent_truncation]
If writing a record fails (I/O error, compression error), the writer MUST return the error immediately. Partial records MUST NOT be written — the BGZF block containing the failed record MUST be discarded, not flushed.

## Performance

r[bam_writer.reuse_buffers]
The writer MUST reuse serialization buffers across `write()` calls. A single `Vec<u8>` for record serialization (cleared and reused each call) avoids per-record allocation.

r[bam_writer.multithreaded_compression]
The writer SHOULD support multi-threaded BGZF compression as a future optimization. The API MUST support a `set_threads(n)` method that configures background compression workers. With `n=0` (default), compression happens synchronously in the calling thread. This matches htslib's `hts_set_threads` pattern and is critical for write-heavy pipelines like BAM rewriting at whole-genome scale.

## Testing

r[bam_writer.test_roundtrip]
A BAM file written by `BamWriter` MUST be readable by seqair's `IndexedBamReader` (after indexing) and by `samtools view`. All record fields (pos, flags, mapq, cigar, seq, qual, aux tags, qname) MUST round-trip exactly.

r[bam_writer.test_samtools_quickcheck]
Output BAM files MUST pass `samtools quickcheck` (validates BGZF structure, EOF block, and basic header integrity).

r[bam_writer.test_header_roundtrip]
Headers written by `BamWriter` MUST be readable by seqair's BAM reader. Reference sequence names and lengths MUST match. PG records MUST be preserved.

r[bam_writer.test_empty_file]
Writing a BAM file with a header but zero records MUST produce a valid BAM file (magic + header + EOF block) that `samtools view` accepts.

r[bam_writer.test_large_record]
Records near the BAM size limit (~2 MiB) MUST serialize and round-trip correctly. Records exceeding the limit MUST be rejected with a clear error.

r[bam_writer.test_cross_validation]
For the initial release, a comparison test MUST write the same records with both seqair's `BamWriter` and `rust_htslib::bam::Writer`, then compare the outputs field-by-field (record contents must match; BGZF block boundaries may differ).
