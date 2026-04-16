# FASTQ Reader

The FASTQ reader provides sequential access to FASTQ files, producing unmapped records for each entry. FASTQ is the standard format for raw sequencing reads before alignment. Seqair supports FASTQ for feature parity with htslib.

FASTQ has no formal specification. The format was invented by Jim Mullikin at the Sanger Institute. The authoritative reference is Cock et al., "The Sanger FASTQ file format for sequences with quality scores" (Nucleic Acids Research, 2010). The quality encoding follows the Sanger/Phred+33 convention used by all modern Illumina instruments.

> **Sources:** Cock et al. 2010, NAR 38(6):1767-1771 (Sanger FASTQ definition). Illumina CASAVA output format documentation. htslib FASTQ implementation (`sam.c`, `test/fastq/`). See [References](./99-references.md).

## Record format

r[fastq.record.structure]
A FASTQ record consists of exactly four lines:
1. Header line: begins with `@`, followed by the read name and optional description (separated by whitespace)
2. Sequence line: raw base letters (A, C, G, T, N, and IUPAC ambiguity codes)
3. Separator line: begins with `+`, optionally followed by the same read name
4. Quality line: ASCII-encoded quality scores, same length as the sequence line

r[fastq.record.name]
The read name is the text between `@` and the first whitespace character (space or tab). The reader MUST preserve the full name without truncation. Everything after the first whitespace on the header line is the comment/description and MUST be available separately.

r[fastq.record.sequence]
The sequence line MUST contain only valid IUPAC nucleotide characters (A, C, G, T, U, R, Y, S, W, K, M, B, D, H, V, N) in upper or lower case. The reader MUST uppercase all bases. Unrecognized characters MUST produce an error.

r[fastq.record.quality]
The quality line MUST have exactly the same number of characters as the sequence line. Each character represents a Phred quality score encoded as ASCII value minus 33 (Sanger/Phred+33 encoding). Valid quality characters range from `!` (ASCII 33, Phred 0) to `~` (ASCII 126, Phred 93).

r[fastq.record.quality_mismatch]
If the quality line length differs from the sequence line length, the reader MUST produce an error identifying the record name and the two lengths.

## Multi-line FASTQ

r[fastq.multiline]
Some FASTQ files use multiple lines for the sequence and/or quality data (e.g., wrapping at 80 characters). The reader MUST support multi-line records by concatenating sequence lines until the `+` separator is encountered, then concatenating quality lines until the quality length equals the sequence length. This matches htslib behavior.

## Output representation

r[fastq.output.record]
Each FASTQ record MUST be represented as a `FastqRecord` struct containing:
- `name`: read name (`Vec<u8>`)
- `comment`: optional description/comment after the name (`Vec<u8>`)
- `seq`: sequence bases (`Vec<u8>`, uppercase ASCII)
- `qual`: quality scores (`Vec<u8>`, raw Phred values, NOT ASCII+33)

The `FastqRecord` type is lightweight and does not carry BAM-specific fields (ref_id, pos, cigar, flags). Conversion to `OwnedBamRecord` MAY be provided as a separate method for callers that need BAM-compatible records.

r[fastq.output.to_bam]
`FastqRecord::to_owned_bam_record()` MUST produce an `OwnedBamRecord` with:
- `ref_id = -1`, `pos = -1` (unmapped)
- `flags` with UNMAPPED (0x4) set, plus READ1/READ2 if detected from the name
- `mapq = 0`, `cigar` empty
- `qname` from the FASTQ read name (with `/1` or `/2` suffix stripped if present)
- `seq` converted to `Vec<Base>`
- `qual` from the raw Phred scores

> **Design note:** An alternative API would yield pairs `(R1, R2)` directly from an interleaved reader. This can be provided as a higher-level wrapper over the individual-record iterator if needed.

## Paired-end support

### Read number detection

r[fastq.paired.readnum]
The reader MUST detect read number (1 or 2) from `/1` or `/2` suffix on the read name. When detected, `to_owned_bam_record()` MUST set the appropriate BAM flags (READ1 = 0x40, READ2 = 0x80, PAIRED = 0x1).

### Interleaved files

r[fastq.paired.interleaved]
In interleaved FASTQ, read 1 and read 2 of each pair alternate: R1, R2, R1, R2, ... The reader MUST support an interleaved mode that validates consecutive records have matching read names (after stripping `/1` and `/2` suffixes). Mismatched names MUST produce an error.

### Split files

r[fastq.paired.split]
Paired-end data is commonly stored in two separate files (R1.fq and R2.fq) with records in the same order. The reader MUST support opening two files simultaneously and iterating them in lockstep. Mismatched record names between files MUST produce an error.

## CASAVA header format

CASAVA is Illumina's pipeline output format. Headers look like: `@instrument:run:flowcell:lane:tile:x:y read:filtered:control:barcode`. The space-separated second field encodes read number (1 or 2), filter status (Y = failed QC, N = passed), control bits, and index/barcode sequence.

r[fastq.casava.parse]
When CASAVA parsing is enabled (opt-in), the reader MUST:
- Set read number from the `read` field (1 or 2)
- Set the FILTERED flag (0x200) when `filtered` is `Y`
- Store the barcode in the `BC` aux tag (or a configurable tag name)

r[fastq.casava.filter]
When CASAVA filtering is enabled, the reader MUST skip records where the filter field is `Y` (failed Illumina chastity filter). This matches htslib's `fastq_casava` option behavior.

## UMI / barcode support

r[fastq.umi]
When UMI parsing is enabled (opt-in, configurable tag name, default `RX`), the reader MUST extract UMI sequences from the read name (last colon-delimited field in Illumina format) and make them available on the `FastqRecord`. When converting to `OwnedBamRecord`, UMIs MUST be stored in the specified aux tag.

## Aux tags in comments

r[fastq.output.aux]
When aux tag parsing is enabled (opt-in), the reader MUST parse SAM-style aux tags from the comment field. Tags are tab-separated, formatted as `TAG:TYPE:VALUE` (e.g., `RG:Z:sample1`). This matches htslib's `fastq_aux` option.

## Compression

r[fastq.compression.detect]
The reader MUST auto-detect compression from magic bytes, not file extensions. BGZF is detected by the gzip magic (0x1f 0x8b) plus the BC subfield; plain gzip is detected by gzip magic without BC.

r[fastq.compression.bgzf]
BGZF-compressed FASTQ files MUST be read using the existing `BgzfReader`. This reuses the proven BGZF decompression path.

r[fastq.compression.gzip]
Plain gzip-compressed FASTQ files MUST be supported by wrapping the input in a gzip decompressor (using the `BgzfReader` in compatibility mode or a separate `flate2`/`miniz_oxide` decoder). Implementation detail TBD.

r[fastq.compression.plain]
Uncompressed FASTQ files MUST be supported without requiring any compression.

## Streaming access

r[fastq.access.streaming]
The FASTQ reader provides streaming (sequential) access only. FASTQ files are not indexed — there is no random access by position or read name. The reader MUST iterate records in file order.

r[fastq.access.no_index]
The FASTQ reader MUST NOT require an index file. Unlike BAM/CRAM, FASTQ is always read sequentially.

r[fastq.access.fqi_future]
htslib defines an FQI (FASTQ index) format for random access into BGZF-compressed FASTQ, analogous to FAI for FASTA. FQI support MAY be added as a follow-up if random-access FASTQ queries prove useful.

## Error handling

r[fastq.errors]
All FASTQ-specific errors MUST contain context (file path, line number, record name) sufficient to diagnose the problem. The implementation MUST NOT use `panic!`, `unwrap()`, or string-based errors. Error variants MUST include at minimum: `InvalidHeader` (line not starting with `@`), `InvalidSeparator` (line not starting with `+`), `QualityLengthMismatch`, `InvalidBase`, `UnexpectedEof`.
