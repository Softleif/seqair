# BAM Record

Each aligned read in a BAM file is stored as a **record**: a binary structure encoding the read's position on the reference genome, its DNA sequence, quality scores, alignment details (CIGAR), and optional auxiliary tags. Records are variable-length — a 150 bp read with few tags is ~300 bytes, while a long read with many tags can be kilobytes.

The record's binary layout starts with 32 fixed bytes (position, flags, lengths, etc.), followed by variable-length fields in order: read name (qname), CIGAR operations, 4-bit packed sequence, quality scores, and auxiliary tags.

> **Sources:** [SAM1] §4.2 "The BAM format" — binary record layout (refID, pos, l_read_name, mapq, bin, n_cigar_op, flag, l_seq, cigar, seq, qual, auxiliary); §4.2.4 "SEQ and QUAL encoding" — 4-bit sequence encoding; §4.2.5 "Auxiliary data encoding" — tag types A/c/C/s/S/i/I/f/Z/H/B; §1.4 "The alignment section: mandatory fields" — FLAG bit definitions. See [References](./99-references.md).

## Decoding

> _[SAM1] §4.2 "The BAM format" — block_size, fixed 32-byte header, variable-length field order_

r[bam.record.decode]
A BAM record MUST be decoded from its binary representation: 32 fixed bytes followed by variable-length qname, CIGAR, sequence, quality, and auxiliary data. The 4-byte block_size prefix is read by the caller.

r[bam.record.max_size]
BAM records have a 2 MiB size cap. The reader MUST reject records whose block_size exceeds this limit rather than allocating an unbounded buffer from untrusted input.

r[bam.record.checked_offsets]
Variable-length field offset calculations (var_start, cigar_end, seq_end, qual_end) MUST use checked arithmetic to prevent wrapping on malformed input. Overflow MUST return a decode error.

r[bam.record.fields]
The decoder MUST extract: tid (i32), pos (i32→i64), mapq (u8), flags (u16), seq_len (u32), n_cigar_ops (u16), qname (NUL-terminated, stored without NUL), packed sequence (4-bit), quality scores, CIGAR operations (packed u32), and raw auxiliary data.

r[bam.record.end_pos]
The decoder MUST precompute the reference end position from the CIGAR. Reference-consuming operations are: M(0), D(2), N(3), =(7), X(8). The end position is `pos + ref_consumed - 1` (inclusive).

r[bam.record.zero_refspan]
When a read has zero reference-consuming CIGAR operations (e.g., pure soft-clip or insertion-only), the reference span is zero. In this case `end_pos` MUST equal `pos` (matching htslib's `bam_cigar2rlen` which returns 0, placing the read at exactly one position). The implementation MUST NOT compute `pos - 1` for such reads.

## Sequence encoding

> _[SAM1] §4.2.4 "SEQ and QUAL encoding" — 4-bit nibble encoding, `=ACMGRSVTWYHKDBN` lookup table_

BAM stores DNA sequences in a compact 4-bit encoding: two bases per byte, high nibble first. This halves storage compared to ASCII but requires decoding before use.

r[bam.record.seq_4bit]
BAM encodes sequences in 4-bit nibbles: two bases per byte (high nibble first). The standard lookup table maps 0→`=`, 1→`A`, 2→`C`, 4→`G`, 8→`T`, 15→`N`, and other values to IUPAC ambiguity codes.

r[bam.record.seq_at]
The record MUST support decoding a single base at a given read position.

r[bam.record.seq_at_simd+2]
The decoder MAY decode the full 4-bit sequence into `Base` enum values at construction time using SIMD-accelerated bulk decoding (see `base_decode.md`). Single-base access via `seq_at(idx, pos)` returns the pre-decoded `Base` value directly.

## Flag access

> _[SAM1] §1.4 "The alignment section: mandatory fields" — FLAG bit table (0x4 unmapped, 0x10 reverse strand, 0x40 first in template, 0x80 second in template)_

Each record has a 16-bit flags field encoding properties like strand orientation, pairing status, and mapping quality. These are checked frequently during pileup construction and filtering.

r[bam.record.flag_reverse]
`is_reverse()` MUST return true when the 0x10 bit is set.

r[bam.record.flag_first]
`is_first_in_template()` MUST return true when the 0x40 bit is set.

r[bam.record.flag_second]
`is_second_in_template()` MUST return true when the 0x80 bit is set.

r[bam.record.flag_unmapped]
`is_unmapped()` MUST return true when the 0x4 bit is set.

## Auxiliary tags

> _[SAM1] §4.2.5 "Auxiliary data encoding" — tag byte layout, type codes A/c/C/s/S/i/I/f/Z/H/B_

BAM records can carry optional key-value tags (e.g. `XR:Z:CT` for bismark strand, `MD:Z:...` for mismatch string). Tags are stored as a flat byte array at the end of the record, each prefixed by a 2-byte name and a 1-byte type code.

r[bam.record.aux_parse]
The record MUST support looking up auxiliary tags by their 2-byte name. Tag types A, c, C, s, S, i, I, f, d, Z, H, and B (typed array) MUST be supported.

r[bam.record.aux_array_unsupported]
B-type (array) auxiliary tags are not yet parsed into `AuxValue`. `find_tag` returns `None` for array tags. The tag bytes MUST be correctly skipped so that all subsequent tags in the record remain accessible.

r[bam.record.raw_aux]
The record MUST provide access to raw auxiliary data bytes for efficient filtering without full tag parsing.
