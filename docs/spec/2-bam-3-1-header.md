# BAM Header

A BAM file begins with a header that describes the reference genome and records metadata about how the alignments were produced. The header contains two parts: a SAM-format text section (with `@HD`, `@SQ`, `@RG`, `@PG` lines) and a binary reference sequence table. The reference table maps each contig (chromosome) to a numeric **target ID** (tid) used throughout the rest of the file — every aligned read stores the tid of its reference contig rather than the name string.

> **Sources:** [SAM1] §4.2 "The BAM format" (BAM magic, l_text, n_ref, reference sequence table); §1.3 "The header section" (SAM header line types @HD/@SQ/@RG/@PG). See [references.md](references.md).

## Parsing

> *[SAM1] §4.2 "The BAM format" — magic, l_text, text, n_ref, l_name, name, l_ref fields*

r[bam.header.magic]
A BAM file MUST begin with the magic bytes `BAM\1` (0x42, 0x41, 0x4d, 0x01) after BGZF decompression. The reader MUST reject files that do not start with this magic.

r[bam.header.text]
After the magic, a 4-byte little-endian `l_text` gives the length of the SAM header text. The reader MUST parse this text to extract `@HD`, `@SQ`, `@RG`, and `@PG` header lines.

r[bam.header.non_negative_lengths]
Header length fields (`l_text`, `n_ref`, `l_name`) are stored as `i32` in the BAM format. The reader MUST validate that these values are non-negative before using them as sizes. Negative values from corrupt files MUST return a `BamHeaderError` rather than causing a massive allocation via `i32 as usize` wrapping.

r[bam.header.references]
After the header text, a 4-byte little-endian `n_ref` gives the number of reference sequences. Each reference is encoded as a 4-byte `l_name` (including NUL), the name bytes, and a 4-byte `l_ref` (sequence length). The reader MUST parse all reference sequences.

## Contig lookup

The rest of the BAM file and its index use tids (integer IDs) to refer to contigs. The caller typically has a contig name (e.g. "chr19") from user input or a BED file and needs the tid for index queries. Conversely, when reporting results, the tid needs to be mapped back to a name.

r[bam.header.tid_lookup]
The header MUST support O(1) lookup from contig name to target ID (tid) and O(1) lookup from tid to contig name and length.

r[bam.header.target_names]
The header MUST provide access to the list of all target names in tid order.

r[bam.header.target_len]
The header MUST provide the length of each target sequence by tid.
