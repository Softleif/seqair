# SAM Reader

A bgzf-compressed, indexed SAM reader for region-based random access. Plain (uncompressed) SAM files cannot be indexed and are NOT supported — users must compress with `bgzip` first.

> **Sources:** [SAM1] §1 "The SAM Format Specification" — SAM text format, mandatory fields, FLAG bits, CIGAR string, header lines (@HD/@SQ/@RG/@PG); §1.3 "The header section"; §1.4 "The alignment section: mandatory fields". [TABIX] — TBI index format and SAM column configuration. [SAM1] §4.1 — BGZF decompression (shared with BAM). The coordinate conversion (1-based → 0-based) and aux-to-BAM-binary serialization are seqair-specific. See [References](./99-references.md).

## Context

SAM is the text representation of the same data BAM encodes in binary. A bgzf-compressed SAM file (`.sam.gz`) wraps SAM text in BGZF blocks, enabling random access via tabix (`.tbi`) or CSI (`.csi`) indexes. htslib reads these transparently; seqair needs to support them for backward compatibility after removing the htslib reading path.

The key difference from BAM: records are text lines that must be parsed field-by-field, and coordinates are 1-based (BAM is 0-based). The BGZF and index layers are shared with BAM.

## Opening

r[sam.reader.open]
The SAM reader MUST open a bgzf-compressed SAM file, parse its header, and locate a tabix (`.tbi`) or CSI (`.csi`) index file. If neither index exists, return an error suggesting `samtools index` or `tabix -p sam`.

r[sam.reader.reject_uncompressed]
If the file does not start with BGZF magic (`1f 8b` with `BC` subfield), the reader MUST reject it with an error explaining that plain SAM files must be bgzf-compressed for indexed access (`bgzip file.sam`).

## Header parsing

> r[sam.header.parse]
> The SAM header consists of all lines starting with `@` at the beginning of the file. The reader MUST:
>
> 1. Read lines from the BGZF stream until a line does NOT start with `@`.
> 2. Record the BGZF virtual offset of the first alignment line (for seeking past the header).
> 3. Parse `@SQ` lines to extract target names (`SN` tag) and lengths (`LN` tag).
> 4. Construct a `BamHeader` from the parsed text via `BamHeader::from_sam_text()`.
>
> Lines starting with `@` that appear after the first alignment line are NOT header continuation — they are malformed alignment lines and MUST be treated as parse errors.

r[sam.header.sq_required]
If no `@SQ` lines are present, the reader MUST return an error. Without `@SQ` lines, target name→tid mapping is impossible and `fetch_into` cannot work.

r[sam.header.preserve_full_text]
The full header text (all `@` lines joined with newlines) MUST be stored in `BamHeader::header_text` for downstream use (e.g., constructing htslib `Header` for BAM writing).

## Alignment line parsing

> _[SAM1] §1.4 "The alignment section: mandatory fields" — 11 mandatory fields: QNAME, FLAG, RNAME, POS (1-based), MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL; QUAL Phred+33 encoding_

> r[sam.record.parse]
> Each alignment line MUST be split on TAB characters into the 11 mandatory fields plus optional tag fields. The parser MUST handle:
>
> | Field | SAM type        | Conversion to internal representation                                      |
> | ----- | --------------- | -------------------------------------------------------------------------- | ---- |
> | QNAME | string          | stored as bytes in name slab                                               |
> | FLAG  | integer         | u16 (identical to BAM flags)                                               |
> | RNAME | string          | converted to tid via header's `name_to_tid` map                            |
> | POS   | 1-based integer | subtract 1 for 0-based internal pos (i64)                                  |
> | MAPQ  | integer         | u8 (255 means unavailable)                                                 |
> | CIGAR | string          | parsed to `Vec<CigarOp>` (typed BAM-on-disk packed layout, see `r[bam.owned_record.cigar_op]`) |
> | RNEXT | string          | `=` means same as RNAME, `*` means unavailable — not stored in RecordStore |
> | PNEXT | integer         | not stored in RecordStore (mate info)                                      |
> | TLEN  | integer         | not stored in RecordStore (template length)                                |
> | SEQ   | string          | each ASCII char decoded to `Base` enum                                     |
> | QUAL  | string          | each char minus 33 for Phred score; `*` means all 0xFF                     |
>
> The parser MUST validate each field and produce clear error messages for malformed records: non-integer POS, FLAG > 65535, negative MAPQ, etc. Real-world SAM files from non-standard tools may contain hand-edited or corrupted records.

r[sam.record.coordinate_conversion]
SAM POS is 1-based; the internal representation is 0-based (matching BAM). The parser MUST subtract 1 from POS. A POS of 0 in SAM means unmapped — after subtraction this becomes -1 (as i64). Unmapped records are filtered by FLAG 0x4 before reaching this point, so -1 positions should not appear in the RecordStore.

r[sam.record.cigar_parse]
The CIGAR string MUST be parsed into a `Vec<CigarOp>`. Each `CigarOp` wraps the BAM packed `u32` (`(length << 4) | op_code`), with op_codes M=0, I=1, D=2, N=3, S=4, H=5, P=6, =7, X=8. Unknown op characters MUST produce a typed parse error (this is stricter than BAM ingest, which tolerates reserved op codes via `CigarOpType::Unknown` per `r[io.typed_cigar_ops]`, because in SAM the lexical character is the only signal). The string `*` means CIGAR unavailable (0 operations, and end_pos = pos).

r[sam.record.seq_decode]
SEQ characters MUST be decoded to `Base` values: A→A, C→C, G→G, T→T, all others (N, IUPAC codes)→Unknown. The string `*` means sequence unavailable (length 0).

Note: the SAM character `=` means "identical to the reference base at this position." This is an intentional lossy conversion to `Unknown` — matching BAM's behavior, which also stores `=` as a distinct 4-bit code but does not resolve it against the reference. Some tools (Illumina DRAGEN) use `=` extensively. If exact reference-matching base recovery is needed in the future, the `=` case would need reference lookup during parsing.

r[sam.record.qual_decode]
QUAL characters MUST have 33 subtracted (Phred+33 encoding). The string `*` means quality unavailable; all quality bytes MUST be set to 0xFF (255). QUAL length MUST equal SEQ length when both are present.

> r[sam.record.aux_tags]
> Optional fields after the 11th column follow the format `TAG:TYPE:VALUE`. They MUST be serialized to BAM binary aux format for storage in the aux slab:
>
> - `A` (char): `[tag0, tag1, 'A', char]`
> - `i` (integer): `[tag0, tag1, type, bytes...]` where type is the smallest fitting BAM integer type:
>   - [-128, 127] → `c` (int8)
>   - [128, 255] → `C` (uint8)
>   - [-32768, -129] or [256, 32767] → `s` (int16)
>   - [32768, 65535] → `S` (uint16)
>   - [-2147483648, -32769] or [65536, 2147483647] → `i` (int32)
>   - [2147483648, 4294967295] → `I` (uint32)
> - `f` (float): `[tag0, tag1, 'f', le_f32_bytes]`
> - `Z` (string): `[tag0, tag1, 'Z', string_bytes, 0x00]`
> - `H` (hex): `[tag0, tag1, 'H', hex_bytes, 0x00]`
> - `B` (array): `[tag0, tag1, 'B', subtype, count_le32, values...]`

Note: the integer type selected here may differ from the type used by the original BAM writer, since SAM text loses the specific width information. This is an inherent limitation — see `r[unified.push_fields_equivalence]`.

r[sam.record.aux_parse_strict]
Malformed aux tag values MUST return an error, not silently default to zero. For integer (`i`) typed tags, if the value cannot be parsed as an integer, the reader MUST return a `SamRecordError::InvalidAuxValue` error.

r[sam.record.aux_int_range+2]
SAM aux integer values (type `i`) are serialized into the smallest BAM integer type that fits, using the ranges defined in `r[sam.record.aux_tags]`. Values that exceed the u32 range (i.e., `val > u32::MAX` or `val < i32::MIN` when interpreted as i64) cannot be represented in any BAM integer type and MUST return a `SamRecordError::AuxIntOutOfRange` error instead of silently truncating.

## Region fetching

> r[sam.reader.fetch_into]
> `fetch_into(tid, start, end, store)` MUST:
>
> 1. Query the tabix/CSI index for BGZF virtual offset ranges overlapping the region.
> 2. Use `RegionBuf` (or equivalent bulk read) to load compressed bytes.
> 3. Decompress BGZF blocks to text.
> 4. Split text into lines (on `\n`).
> 5. Parse each alignment line.
> 6. Filter: skip lines where RNAME's tid doesn't match, or where the record doesn't overlap `[start, end]`.
> 7. Skip unmapped reads (FLAG 0x4).
> 8. Push passing records into the RecordStore.

r[sam.reader.overlap_filter+2]
The overlap filter uses half-open intervals (0-based). `end_pos` MUST be computed from the parsed CIGAR. When CIGAR is `*` (unavailable), `end_pos = pos` (point record).

r[sam.reader.overlap_halfopen]
A record overlaps `[start, end)` iff `pos < end && end_pos > start`. Records where `pos == end` or `end_pos == start` do NOT overlap and MUST be filtered out.

r[sam.reader.sorted_order]
Records MUST be added to the store in coordinate-sorted order. Since bgzf-compressed indexed SAM files are coordinate-sorted by construction (indexing requires it), the natural parse order is correct.

## Tabix / CSI index

> _[TABIX] — TBI magic `TBI\x01`, header fields, n_ref, bin/chunk/linear-index structure, format=1 for SAM_
> _[CSI] — CSI magic `CSI\x01`, min_shift, depth, l_aux, per-bin loffset_

r[sam.index.tabix]
The reader MUST support tabix (`.tbi`) indexes. Tabix is a generic index for bgzf-compressed, TAB-delimited files. For SAM, the format code is `1` (not `0`/generic). The column configuration (seq_col=3, beg_col=4, end_col=0) is implicit for format=1 and may be ignored by some implementations; the reader SHOULD NOT reject a tabix file based on stored column values when format=1 is present.

Tabix files are BGZF-compressed and MUST be decompressed before parsing.

r[sam.index.csi+2]
The reader SHOULD also support CSI (`.csi`) indexes, which generalize the binning scheme with configurable `min_shift` and `depth` parameters. CSI supports references longer than 512 Mbp (BAI's limit). CSI files are BGZF-compressed and MUST be decompressed before parsing.

r[sam.index.locate+2]
Index file lookup order: `<path>.tbi`, `<path>.bai`. If neither exists, return a clear error. CSI (`.csi`) lookup should be added when CSI support is implemented.

## Edge cases

r[sam.edge.missing_seq]
When SEQ is `*`, the record has no sequence data. `seq_len` MUST be 0 and no bases are pushed to the bases slab. This is valid per the SAM spec for secondary alignments that omit sequence.

r[sam.edge.missing_qual]
When QUAL is `*`, all quality bytes MUST be 0xFF. This matches BAM's representation of unavailable quality.

r[sam.edge.long_cigar]
BAM encodes CIGARs with >65535 operations using a `CG:B:I` aux tag and a placeholder CIGAR. SAM has no such limitation — the CIGAR string can be arbitrarily long. The parser MUST handle CIGARs of any length.

r[sam.edge.empty_lines]
Blank lines (empty or whitespace-only) between alignment records MUST be silently skipped.

r[sam.edge.rname_star]
RNAME of `*` means unmapped. Combined with FLAG 0x4, these records MUST be filtered out by `fetch_into` (same as BAM's unmapped filtering). RNAME `*` without FLAG 0x4 is malformed but SHOULD be treated as unmapped.

r[sam.edge.line_spanning_blocks]
A single SAM alignment line may span multiple BGZF blocks. The reader MUST buffer partial lines across block boundaries and only parse complete lines (terminated by `\n`). The parser MUST also strip `\r` at line boundaries to handle `\r\n` (Windows-style) line endings that may appear in files from mixed-platform pipelines.

## Integer parsing

r[sam.record.parse_int_nonempty]
Integer parsing functions (`parse_u8`, `parse_u16`, `parse_u32`) MUST reject empty input by returning `None`. An empty byte slice is not a valid integer representation and MUST NOT silently produce zero.

## Performance considerations

r[sam.perf.bulk_read]
Like BAM, the SAM reader SHOULD use `RegionBuf`-style bulk I/O: read all compressed bytes for the region in one large I/O operation, then decompress and parse from memory. This is critical for NFS/Lustre performance.

Text parsing is inherently slower than BAM binary decoding. The SAM reader is expected to be 2-5× slower than BAM for the same data due to text parsing overhead and larger compressed size. This is acceptable since SAM is not the recommended input format. Implementation guidance:

- Use `memchr` for fast TAB and newline scanning.
- Avoid allocations during parsing (reuse line buffers).
- Parse integers without going through `String` (direct byte-to-int conversion).
- Decode SEQ characters to Base with a lookup table (not per-char branching).
