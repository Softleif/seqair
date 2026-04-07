# CRAM Reader

A CRAM reader for region-based random access. CRAM is a reference-based, column-oriented compression format that achieves significantly smaller file sizes than BAM (typically 40-60% of BAM size) at the cost of decoding complexity and reference dependency.

> **Sources:** [CRAM3] throughout — file structure, containers, slices, blocks, variable-length integers, bit stream, compression header, data series encodings, record decode order, read features, CRAI index. [CRAMcodecs] — rANS 4×8 (§2), rANS Nx16 (§3), arithmetic coder (§4), tok3 (§5), fqzcomp (§6). See [References](./99-references.md).

## Context

CRAM stores alignment data as differences against a reference genome rather than full sequences. Where BAM stores `ACGTACGT` for a read, CRAM stores "matches reference at this position, except substitution at offset 5". This reference-based approach, combined with column-oriented storage and sophisticated entropy coding, makes CRAM the most compact alignment format.

Rastair needs CRAM support because many large-scale sequencing projects store data exclusively in CRAM format.

## Goals and scope

r[cram.scope.versions]
The reader MUST support CRAM v3.0 and v3.1. CRAM v2.0 and v2.1 MAY be supported but are not required — they are rarely encountered in practice. The magic bytes `CRAM` followed by major version 3 MUST be accepted; major version 2 SHOULD produce a warning; other versions MUST be rejected.

The reader is read-only. CRAM writing is out of scope — rastair writes BAM (via htslib) for output.

r[cram.scope.reference_required]
CRAM decoding requires access to the reference genome (FASTA). The reader MUST accept a FASTA reader (the one being built on the other branch) for sequence lookups. Embedded reference slices (where the reference is stored in the CRAM container itself) MUST also be supported as a fallback.

Note: htslib supports automatic reference download via `REF_PATH`/`REF_CACHE` environment variables and the EBI reference server. This is out of scope for seqair, but when the FASTA reader is missing or the required contig is absent, the error message SHOULD mention `REF_PATH`/`REF_CACHE` as an alternative for users familiar with the htslib workflow.

## File structure

> _[CRAM3] §6 "File definition" — CRAM magic, major/minor version, file ID_

### File definition

r[cram.file.magic]
The file MUST start with `CRAM` (bytes `43 52 41 4d`), followed by major version (1 byte) and minor version (1 byte), followed by a 20-byte file ID.

> r[cram.file.header_container]
> The first container after the file definition MUST be the file header container. It contains a single block (content_type=FILE_HEADER) whose data is:
>
> - `header_text_length` (i32, fixed-width little-endian): byte length of the SAM header text
> - `header_text` (UTF-8 bytes): the SAM header text (`@HD`, `@SQ`, `@RG`, etc.)

The header text is parsed via `BamHeader::from_sam_text()` to populate target names, lengths, and tid mappings.

**Note:** The i32 length prefix is part of the block's _decompressed data_, not the block header itself.

r[cram.file.eof]
The file MUST end with an EOF container. The EOF container produced by samtools has `length=15` (not 0) and `num_records=0`. The reader MUST detect the EOF container as any container where `num_records=0` and `length <= 15`. The reader SHOULD verify the EOF marker but MUST NOT fail if it's missing (truncated files should produce a warning, not a hard error).

### Containers

> _[CRAM3] §7 "Container header structure" — length (fixed i32), itf8 fields, landmarks, CRC32_

> r[cram.container.header]
> Each container header contains:
>
> - `length` (i32, **fixed-width** little-endian): total byte length of all blocks in the container. This is the _only_ fixed-width field in the header; all subsequent fields use variable-length encoding.
> - `reference_sequence_id` (itf8): reference ID, or -1 (unmapped), or -2 (multi-ref)
> - `alignment_start` (itf8): 1-based start position on reference
> - `alignment_span` (itf8): number of reference bases covered
> - `num_records` (itf8): number of alignment records
> - `record_counter` (ltf8): global 0-based record index (v3.0+)
> - `bases` (ltf8): total bases in container (v3.0+)
> - `num_blocks` (itf8): block count
> - `landmarks` (itf8 array): count (itf8) followed by that many itf8 values; byte offsets from container data start to each slice header
> - `crc32` (4 bytes, **fixed-width** little-endian): CRC32 over all preceding header bytes (from `length` through end of `landmarks`)

r[cram.container.region_skip]
During `fetch_into`, the reader MUST skip containers whose `[alignment_start, alignment_start + alignment_span)` range does not overlap the query region. The CRAI index provides byte offsets for this.

### Slices

> _[CRAM3] §8.4 "Slice header block" — slice header fields, multi-ref slices, embedded reference, reference_md5_

> r[cram.slice.header]
> Each slice header contains:
>
> - `reference_sequence_id` (itf8): ref ID for this slice, -1 unmapped, -2 multi-ref
> - `alignment_start`, `alignment_span` (itf8): genomic range covered
> - `num_records` (itf8): records in this slice
> - `record_counter` (ltf8): global record counter
> - `num_blocks` (itf8): data blocks in this slice
> - `block_content_ids` (itf8 array): content IDs of the blocks
> - `embedded_reference` (itf8): -1 if no embedded reference, else block content ID
> - `reference_md5` (16 bytes): MD5 of the reference segment (all zeros if N/A)
> - optional tags (BAM-format aux bytes)

r[cram.slice.multi_ref]
When `reference_sequence_id == -2`, the slice contains records from multiple references. The `RI` (reference ID) data series MUST be decoded per-record to determine which reference each record aligns to. Multi-ref slices are common in coordinate-sorted CRAM files at chromosome boundaries and for reads with unmapped mates.

r[cram.slice.embedded_ref]
When `embedded_reference >= 0`, the corresponding external data block contains the reference bases for this slice. The reader MUST use these instead of the FASTA reader. This enables CRAM files to be self-contained (produced by `samtools view --output-fmt-option embed_ref` or cramtools).

### Blocks

> _[CRAM3] §8 "Block structure" — method, content_type, content_id, compressed_size, uncompressed_size, data, CRC32 (v3.0+)_

> r[cram.block.structure]
> Each block contains:
>
> - `method` (byte): compression codec (see codecs section)
> - `content_type` (byte): 0=FILE_HEADER, 1=COMPRESSION_HEADER, 2=SLICE_HEADER, 4=EXTERNAL, 5=CORE
> - `content_id` (itf8): identifies which data series this block holds
> - `compressed_size`, `uncompressed_size` (itf8)
> - `data` (bytes): compressed payload
> - `crc32` (4 bytes): CRC32 of the block (v3.0+)

r[cram.block.crc32]
CRC32 MUST be verified on every block in v3.0+. The CRC32 covers all bytes from the start of the block (method byte) through the end of the compressed data, excluding the 4 CRC32 bytes themselves. Mismatches MUST produce an error with the block's content_type and content_id for diagnostics.

## Variable-length integers

> _[CRAM3] §2.1 "Logical data types" — ITF8 (up to 32 bits) and LTF8 (up to 64 bits) variable-length integer encodings_

> r[cram.itf8]
> ITF8 (up to 32 bits) MUST be decoded as:
>
> - `0xxxxxxx`: 1 byte, value = bits 6..0 (7 bits)
> - `10xxxxxx + 1 byte`: 2 bytes, value = `(b0 & 0x3F) << 8 | b1` (14 bits)
> - `110xxxxx + 2 bytes`: 3 bytes, value = `(b0 & 0x1F) << 16 | b1 << 8 | b2` (21 bits)
> - `1110xxxx + 3 bytes`: 4 bytes, value = `(b0 & 0x0F) << 24 | b1 << 16 | b2 << 8 | b3` (28 bits)
> - `1111xxxx + 4 bytes`: 5 bytes, value = `(b0 & 0x0F) << 28 | b1 << 20 | b2 << 12 | b3 << 4 | (b4 & 0x0F)` (32 bits)

Note: in the 5-byte form, only the low 4 bits of byte 4 are used. The high 4 bits of byte 4 are unused/ignored.

The decoded value is an unsigned u32. For fields that can be negative (e.g., `reference_sequence_id` where -1=unmapped, -2=multi-ref), the u32 MUST be reinterpreted as signed i32 via two's complement (`value as i32`).

r[cram.ltf8]
LTF8 (up to 64 bits) MUST be decoded with the same leading-bit pattern extended to 9 bytes, with `0xFF` prefix indicating 8 continuation bytes.

## Bit stream reading

> _[CRAM3] §2.2 "Reading and writing bits in a bit stream" — MSB-first bit ordering_

> r[cram.bitstream]
> The core data block (content_type=5) is read as a bit stream. Bits MUST be read MSB-first (most significant bit first) within each byte:
>
> - Bit offset 0 is the most significant bit (0x80)
> - Bit offset 7 is the least significant bit (0x01)
> - Reading N bits shifts them out from the current position, crossing byte boundaries as needed

Encodings that read from the core (Huffman, Beta, Subexp, Gamma) all use this bit stream. Encodings that read from external blocks (External, ByteArrayStop) use byte-level access to the specific external block identified by `content_id`.

The core bit stream is shared across all records in a slice — fields are interleaved sequentially. Reading fields in the wrong order corrupts all subsequent data.

## Compression header

> _[CRAM3] §8.3 "Compression header block" — preservation map, data series encoding map, tag encoding map_

The compression header is a block (content_type=1) at the start of each container. It defines how all data in the container's slices is encoded.

### Preservation map

> _[CRAM3] §8.3 "Preservation map" — RN, AP, RR, SM (substitution matrix), TD (tag dictionary)_

> r[cram.compression.preservation]
> The preservation map MUST be decoded as a series of key-value pairs:
>
> - `RN` (bool): whether read names are preserved (if false, names are generated)
> - `AP` (bool): whether alignment positions are delta-coded (true) or absolute
> - `RR` (bool): whether the reference is required for decoding
> - `SM` (5 bytes): substitution matrix (see `r[cram.compression.substitution_matrix]`)
> - `TD` (bytes): tag dictionary — null-separated lists of 3-byte entries (2-char tag name + 1-char BAM type)

> r[cram.compression.substitution_matrix]
> The substitution matrix is 5 bytes, one per reference base in order: A, C, G, T, N. Each byte encodes the ordering of substitution alternatives using 2-bit codes packed MSB-first:
>
> `byte = (alt0 << 6) | (alt1 << 4) | (alt2 << 2) | alt3`
>
> where each 2-bit value maps to a base: 0=A, 1=C, 2=G, 3=T. When decoding a substitution (BS data series value = code), the read base is `\alt[code]` from this ordering.

Example: for reference base C, the possible substitutions are A, G, T (not C itself). If the substitution byte is `0b_00_10_11_01` (0x2D):

- code 0 → bits `00` → A
- code 1 → bits `10` → G
- code 2 → bits `11` → T
- code 3 → bits `01` → C (unused for ref=C, but the slot exists)

For non-N reference bases, only codes 0-2 are meaningful (3 alternatives, excluding the reference base itself). Code 3 exists in the byte but is never emitted by conforming CRAM writers — the decoder MAY treat it as an error or map it to the remaining unused base. For reference base N, all 4 codes are meaningful (A, C, G, T).

### Data series encoding map

> _[CRAM3] §8.3 "Data series encodings" — two-byte data series codes and their types (BF, CF, RI, RL, AP, RG, RN, MF, NS, NP, TS, NF, FN, FC, FP, MQ, BA, QS, BS, IN, SC, DL, RS, PD, HC, BB, QQ, TL)_

> r[cram.compression.ds_encodings]
> The data series encoding map MUST be decoded as key-value pairs where:
>
> - Key: 2-byte data series code
> - Value: encoding descriptor (encoding ID + parameters)

Complete data series list:

| Code | Name              | Type       | Description                                        |
| ---- | ----------------- | ---------- | -------------------------------------------------- |
| `BF` | BAM flags         | int        | BAM FLAG field (u16)                               |
| `CF` | CRAM flags        | int        | CRAM-specific flags                                |
| `RI` | Ref ID            | int        | Reference sequence ID (multi-ref slices only)      |
| `RL` | Read length       | int        | Number of bases                                    |
| `AP` | Alignment pos     | int        | Position (delta or absolute per preservation map)  |
| `RG` | Read group        | int        | Index into @RG header entries                      |
| `RN` | Read name         | byte array | QNAME                                              |
| `MF` | Mate flags        | int        | 0x1=mate reverse, 0x2=mate unmapped                |
| `NS` | Next segment ref  | int        | Mate reference ID                                  |
| `NP` | Next mate pos     | int        | Mate position                                      |
| `TS` | Template size     | int        | Insert size                                        |
| `NF` | Next fragment     | int        | Record-skip to next fragment (attached mates)      |
| `FN` | Feature count     | int        | Number of read features                            |
| `FC` | Feature code      | byte       | Feature type character                             |
| `FP` | Feature position  | int        | Position within read (delta from previous feature) |
| `MQ` | Mapping quality   | int        | MAPQ                                               |
| `BA` | Base              | byte       | Single base                                        |
| `QS` | Quality score     | byte       | Single quality value                               |
| `BS` | Substitution code | byte       | Index into substitution matrix                     |
| `IN` | Insertion         | byte array | Inserted bases                                     |
| `SC` | Soft clip         | byte array | Soft-clipped bases                                 |
| `DL` | Deletion length   | int        | Number of deleted reference bases                  |
| `RS` | Ref skip          | int        | Reference skip length (intron)                     |
| `PD` | Padding           | int        | Padding length                                     |
| `HC` | Hard clip         | int        | Hard clip length                                   |
| `BB` | Bases block       | byte array | Stretch of bases                                   |
| `QQ` | Quality block     | byte array | Stretch of quality scores                          |
| `TL` | Tag line          | int        | Index into tag dictionary                          |

### Tag encoding map

> _[CRAM3] §8.3 "Tag encodings" — ITF8 key as `(char1 << 16) + (char2 << 8) + type_char`, BAM-format tag values_

r[cram.compression.tag_encodings]
The tag encoding map uses ITF8 keys computed as `(char1 << 16) + (char2 << 8) + type_char`. Each key maps to an encoding descriptor for that tag type. Tag values are stored as raw BAM-format bytes (without the 3-byte tag+type prefix).

Base modification tags (MM/ML) use `Z` and `B:C` types respectively and can be very large on long reads (thousands of bytes). The decoder MUST handle arbitrary-length tag values.

## Encoding types

> _[CRAM3] §13 "Encodings" — NULL (0), EXTERNAL (1), HUFFMAN (3), BYTE_ARRAY_LEN (4), BYTE_ARRAY_STOP (5), BETA (6), SUBEXP (7), GAMMA (9)_

r[cram.encoding.null]
Encoding ID 0 (NULL): no data. Used for absent data series.

r[cram.encoding.external]
Encoding ID 1 (EXTERNAL): reads raw bytes from an external data block identified by `content_id`. The most common encoding — used for quality scores, bases, tags, etc.

r[cram.encoding.huffman]
Encoding ID 3 (HUFFMAN): canonical Huffman coding. Parameters: alphabet (itf8 array) and bit lengths (itf8 array). Reads from the core data bit stream.

Canonical Huffman construction: sort by (bit_length, symbol_value), assign codes starting from 0 with incrementing and left-shifting when bit length increases. A single-symbol Huffman code has bit length 0 and consumes 0 bits from the stream — this is common for fields where all records have the same value (e.g., all flags identical).

r[cram.encoding.byte_array_len]
Encoding ID 4 (BYTE_ARRAY_LEN): a length encoding followed by a value encoding. The length is decoded first (from any encoding), then that many bytes are decoded using the value encoding.

r[cram.encoding.byte_array_stop]
Encoding ID 5 (BYTE_ARRAY_STOP): reads bytes from an external block until a stop byte is encountered. Used for read names (stop=0x00) and variable-length strings.

r[cram.encoding.beta]
Encoding ID 6 (BETA): fixed-width integer. Parameters: offset and number of bits. Reads from core bit stream.

> r[cram.encoding.subexp]
> Encoding ID 7 (SUBEXP): sub-exponential code. Parameters: `offset` (itf8) and `K` (itf8). Reads from core bit stream. Decode procedure:
>
> 1. Read unary: count consecutive `1` bits until a `0` bit is read. Call this count `n`.
> 2. If `n == 0`: read `K` bits from the stream as the value.
> 3. If `n > 0`: read `n + K - 1` bits from the stream as `remainder`. The value is `(1 << (n + K - 1)) - (1 << K) + remainder`.
> 4. The final decoded value is `value - offset`.

> r[cram.encoding.gamma]
> Encoding ID 9 (GAMMA): Elias gamma code with offset. Parameters: `offset` (itf8). Reads from core bit stream. Decode procedure:
>
> 1. Read unary: count consecutive `0` bits until a `1` bit is read. Call this count `n`.
> 2. If `n == 0`: the value is 0.
> 3. If `n > 0`: read `n` more bits from the stream as `suffix`. The value is `(1 << n) | suffix - 1`.
> 4. The final decoded value is `value - offset`.

## Compression codecs

> _[CRAM3] §14 "External compression methods" — RAW (0), GZIP (1), BZIP2 (2), LZMA (3), rANS4x8 (4), rANS Nx16 (5), arithmetic coder (6), fqzcomp (7), tok3 (8)_
> _[CRAMcodecs] §2 "rANS 4x8" — frequency table RLE, order-0/order-1 modes, chunk-based order-1 interleaving_
> _[CRAMcodecs] §3 "rANS Nx16" — flags byte (ORDER, N32, STRIPE, NoSize, CAT, RLE, PACK)_
> _[CRAMcodecs] §5 "Name tokenisation codec" — tok3_

r[cram.codec.raw]
Method 0 (RAW): uncompressed data.

r[cram.codec.gzip]
Method 1 (GZIP): deflate compression. MUST be supported — this is the most common codec.

r[cram.codec.bzip2]
Method 2 (BZIP2): bzip2 compression. MUST be supported for v3.0 compatibility.

r[cram.codec.lzma]
Method 3 (LZMA): LZMA compression. MUST be supported for v3.0 compatibility.

r[cram.codec.rans4x8]
Method 4 (rANS 4×8): rANS entropy coder with 4-way interleaved states and 8-bit renormalization (L=2^23). Supports order-0 (context-free) and order-1 (previous-symbol context) modes. MUST be supported — this is the default codec in CRAM v3.0 writers including samtools.

The frequency table header uses a compact run-length encoding: after reading a symbol and its frequency, if the next symbol is `prev_sym + 1`, a run length byte follows, then that many consecutive symbols' frequencies. The terminator is a NUL (0x00) byte. Order-1 nests this same structure: the outer loop reads context symbols (with the same run-length scheme), and for each context reads an order-0 frequency table.

Order-1 interleaving is chunk-based, not round-robin: the output is split into 4 equal segments, each decoded by one state sequentially. Remainder bytes (when `len % 4 != 0`) are handled by state 3.

r[cram.codec.rans_nx16]
Method 5 (rANS Nx16): v3.1 codec. rANS with N-way interleaving (N=4 or 32), 16-bit renormalization (L=2^15), and optional transforms controlled by an 8-bit flags byte: ORDER(1), N32(4), STRIPE(8), NoSize(16), CAT(32), RLE(64), PACK(128). Transforms are applied in order: decode entropy → unRLE → unPack. MUST be supported for v3.1 files — samtools uses this for most data blocks in v3.1 output.

r[cram.codec.arith]
Method 6 (arithmetic coder): v3.1 adaptive arithmetic coder. SHOULD be supported.

r[cram.codec.fqzcomp]
Method 7 (fqzcomp): v3.1 quality-score-specific compressor. SHOULD be supported.

r[cram.codec.tok3]
Method 8 (tok3): v3.1 read-name tokeniser. MUST be supported for v3.1 files — samtools uses tok3 for read name blocks by default in v3.1 output. Without tok3 support, v3.1 CRAM files produced by `samtools view -C` cannot be read.

r[cram.codec.unknown]
Unknown codec methods MUST produce a clear error naming the method ID and suggesting conversion to BAM (`samtools view -b`).

## Record decoding

> _[CRAM3] §10 "Record structure" — CRAM record layout and field decode order_

### Decoding order

> _[CRAM3] §10.1 "CRAM record" — strict field decode order (BF, CF, RI, RL, AP, RG, QS, RN, mate info, TL, FN/FC/FP, MQ, per-base QS)_

> r[cram.record.decode_order]
> Records within a slice MUST be decoded in a strict field order, because the core bit stream is sequential. Reading fields out of order corrupts all subsequent data. The order is:
>
> 1. `BF` (BAM flags)
> 2. `CF` (CRAM flags)
> 3. `RI` (reference ID — only if slice is multi-ref, i.e., ref_seq_id == -2)
> 4. `RL` (read length)
> 5. `AP` (alignment position — delta or absolute per preservation map)
> 6. `RG` (read group)
> 7. `QS` (quality scores — only if CRAM flags bit 0x1 is clear, meaning quality stored per-record not per-base)
> 8. `RN` (read name — only if preservation map `RN=true`)
> 9. Mate information:
>    - If CRAM flags bit 0x2 set (detached): `MF`, `NS`, `NP`, `TS`
>    - If CRAM flags bit 0x4 set (mate downstream): `NF`
> 10. `TL` (tag line index)
> 11. `FN` (feature count), then for each feature: `FC`, `FP`, feature-specific data
> 12. `MQ` (mapping quality)
> 13. `QS` per-base (only if CRAM flags bit 0x1 is set — quality as array)

**Note:** Getting this order wrong is the most common implementation bug in CRAM readers.

### Per-record data series

> r[cram.record.flags]
> `BF` (BAM flags) and `CF` (CRAM flags) MUST be decoded per record. BAM flags map directly to the u16 flags in `SlimRecord`. CRAM flags are separate and encode:
>
> - Bit 0x1: quality scores stored as array (vs per-record)
> - Bit 0x2: detached (mate info stored inline)
> - Bit 0x4: has mate downstream in same slice
> - Bit 0x8: sequence not stored (unknown)

r[cram.record.position]
`AP` (alignment position): if the preservation map `AP=true` (delta mode), each record's position is a delta from the previous record's position within the slice. The reader MUST maintain a running position accumulator and convert to absolute 0-based coordinates (CRAM positions are 1-based). Delta-coded positions can be negative (valid for supplementary alignments in coordinate-sorted CRAM).

r[cram.record.read_length]
`RL` (read length): the number of bases in the read. Used to allocate sequence and quality arrays. Long reads (PacBio/ONT) can be 10,000-1,000,000+ bases — buffer pre-allocation MUST NOT assume short-read sizes.

r[cram.record.read_group]
`RG` (read group): an index into the `@RG` entries in the header. Stored as itf8. The reader MUST convert to a `RG:Z:<id>` aux tag with the read group ID string.

r[cram.record.read_name]
`RN` (read name): if preservation map `RN=true`, names are stored. If `RN=false`, names are auto-generated. Rastair uses qnames for mate-pair dedup, so `RN=false` CRAM files MUST produce a warning that dedup will not work correctly. Generated names SHOULD be deterministic (e.g., `cram_<record_counter>`).

### Read feature decoding

> _[CRAM3] §10.6 "Mapped reads" and §10.6.1 "Read feature codes" — feature codes X/I/i/D/S/B/b/Q/q/N/H/P and their data series_

> r[cram.record.features]
> For each record, `FN` gives the number of read features. For each feature:
>
> 1. Decode `FC` (feature code byte)
> 2. Decode `FP` (position within read, delta-coded from the previous feature's position within the same record — first feature's FP is absolute within the read, starting at 1)
> 3. Decode feature-specific data:
>
> | FC  | Char | Data series | Description                                                |
> | --- | ---- | ----------- | ---------------------------------------------------------- |
> | `X` | 0x58 | BS          | Substitution: `base = substitution_matrix[ref_base][code]` |
> | `I` | 0x49 | IN          | Insertion: multi-base insertion                            |
> | `i` | 0x69 | BA          | Single-base insertion                                      |
> | `D` | 0x44 | DL          | Deletion: skip DL reference bases                          |
> | `S` | 0x53 | SC          | Soft clip: bases not aligned to reference                  |
> | `B` | 0x42 | BA + QS     | Single base replacement + quality                          |
> | `b` | 0x62 | BB + QQ     | Stretch of bases + qualities                               |
> | `Q` | 0x51 | QS          | Quality score at position                                  |
> | `q` | 0x71 | QQ          | Stretch of quality scores                                  |
> | `N` | 0x4e | RS          | Reference skip (intron): skip RS reference bases           |
> | `H` | 0x48 | HC          | Hard clip: record HC bases clipped                         |
> | `P` | 0x50 | PD          | Padding: PD padding bases                                  |

### Sequence reconstruction

> _[CRAM3] §10.6 "Mapped reads" — reference-based sequence reconstruction from read features; §11 "Reference sequences" — reference MD5 verification_

> r[cram.record.sequence]
> Read sequences are NOT stored directly. They MUST be reconstructed from the reference and read features:
>
> 1. Fetch the reference sequence for `[alignment_start, alignment_start + alignment_span)` from the FASTA reader (or embedded reference block).
> 2. Initialize the read's bases array with `RL` entries.
> 3. Walk the reference and read in parallel, applying features in order:
>    - Between features: copy reference bases to read bases (these are the matching segments).
>    - `X` (substitution): look up `substitution_matrix[ref_base][BS_code]` for the read base.
>    - `I`/`i` (insertion): insert bases from `IN`/`BA` into the read (these don't consume reference).
>    - `D` (deletion): skip `DL` reference bases (these don't consume read).
>    - `S` (soft clip): insert bases from `SC` into the read (no reference consumption).
>    - `B` (base+quality): replace base and quality at the feature position.
>    - `N` (ref skip): skip `RS` reference bases.
>    - `H` (hard clip): record length only (no bases in read or reference).
> 4. Convert reconstructed bases to `Base` enum values for the bases slab.

r[cram.record.quality]
Quality scores are decoded from `QS` (per-base) or `QQ` (stretch) data series, depending on CRAM flags bit 0x1. If CRAM flags bit 0x8 is set (sequence unknown), quality is unavailable (all 0xFF). Quality length MUST equal read length.

r[cram.feature.quality_score]
When `quality_as_array` is false (CRAM flags bit 0x1 unset), per-position quality overrides are encoded as `Q` (0x51) read features. Each `Q` feature carries a single quality score from the `QS` data series at a specific read position. The decoder MUST apply these quality scores to `qual_buf` at the feature's read position (1-based `FP`, converted to 0-based index). Positions without a `Q` feature retain the default quality (0xFF). This is distinct from the bulk `QS` decode path used when `quality_as_array` is true.

r[cram.feature.quality_block]
The `q` (0x71) read feature carries a stretch of quality scores from the `QQ` data series, starting at the feature's read position. The decoder MUST write the decoded quality bytes into `qual_buf` starting at the feature's 0-based read position. The number of quality scores in the block is determined by the `QQ` encoding. Like `Q`, this feature only appears when `quality_as_array` is false.

> r[cram.record.cigar_reconstruction]
> CRAM does not store CIGAR explicitly. The CIGAR MUST be reconstructed from the read features:
>
> - Stretches of matching reference → M ops (or =/X if distinguishing is needed, but M is standard)
> - `X` (substitution) → part of M ops (substitutions are still alignment matches)
> - `I`/`i` → I ops
> - `D` → D ops
> - `S` → S ops
> - `N` → N ops
> - `H` → H ops
> - `P` → P ops

The reconstructed CIGAR MUST be packed into BAM u32 format for the data slab.

### Mate information

> _[CRAM3] §10.4 "Mate records" — detached (CF bit 0x2) vs attached (NF) mate encoding_

r[cram.record.mate_detached]
When CRAM flags bit 0x2 is set (detached), mate information is stored inline: `MF` (mate flags), `NS` (next ref ID), `NP` (next mate position), `TS` (template size). These fields are NOT needed for RecordStore (which doesn't store mate info), but the decoder MUST read and discard them to advance the data stream correctly.

r[cram.record.mate_attached]
When bit 0x2 is clear, the mate is "attached" — `NF` (next fragment distance) gives the record-skip to the mate within the same slice. The decoder MUST read `NF` but does not need to resolve the linkage for RecordStore purposes.

### Auxiliary tags

> _[CRAM3] §10.5 "Auxiliary tags" — TL tag line index, tag dictionary, tag encoding map_

r[cram.record.aux_tags]
`TL` (tag line index) selects which tag combination this record has from the tag dictionary in the preservation map. For each tag in the combination, the tag encoding map provides the encoding for its value. Tag values MUST be decoded and serialized to BAM binary aux format for storage in the data slab.

r[cram.record.rg_tag]
The `RG` data series is separate from aux tags. If a read group is present, the reader MUST emit an `RG:Z:<id>` aux tag using the read group ID from the header's `@RG` entries. This tag MUST be included in the aux slab alongside dictionary-decoded tags.

## CRAM index (.crai)

> _[CRAM3] §12 "Indexing" — CRAI format: gzip-compressed TSV, 6 fields per line (ref_id, alignment_start, alignment_span, container_offset, slice_offset, slice_size)_

> r[cram.index.parse]
> The `.crai` file is gzip-compressed TAB-delimited text. Each line has 6 fields:
>
> 1. Reference sequence ID (int)
> 2. Alignment start (int, **1-based** matching the container header convention, 0 for unmapped)
> 3. Alignment span (int, 0 for unmapped)
> 4. Container byte offset in file (int, 0-based byte offset from start of file)
> 5. Slice byte offset relative to container data start (int)
> 6. Slice size in bytes (int)

r[cram.index.query]
Region queries MUST find all index entries where `[start, start+span)` overlaps the query region for the given reference ID. The reader then seeks to the container byte offset, skips to the slice offset, and decodes the slice.

r[cram.index.zero_span]
CRAI entries with `alignment_span == 0` (but `alignment_start > 0`) indicate unknown extent — this occurs when samtools writes CRAM with embedded references or when the span is not computed. These entries MUST be included in query results whenever `entry.alignment_start < query_end`, since the slice may contain records at any position past the start. Treating span=0 as "zero-width" would silently drop all records in that slice.

r[cram.index.multi_ref_slices]
Multi-ref slices produce multiple index entries (one per reference they span). A query for one reference may hit a multi-ref slice that also contains records from other references. The reader MUST filter records by reference ID after decoding.

r[cram.index.unmapped]
Unmapped records (ref ID = -1) have start=0, span=0. They are only returned if explicitly queried. For `fetch_into`, unmapped records MUST be skipped (same as BAM).

## Edge cases

r[cram.edge.reference_mismatch]
If the reference MD5 in the slice header does not match the MD5 of the FASTA sequence for that region, the reader MUST return an error. The MD5 is computed over the uppercase reference sequence for the exact range `[alignment_start, alignment_start + alignment_span)`, with no newlines or other formatting. If the range extends beyond the reference length, this is a file/reference mismatch error.

r[cram.edge.missing_reference]
If the FASTA reader cannot provide the reference sequence for a slice's region (e.g., contig not in FASTA), the reader MUST return an error unless the slice has an embedded reference. The error SHOULD mention `REF_PATH`/`REF_CACHE` as alternatives.

r[cram.edge.unknown_read_names]
When `RN=false` in the preservation map, read names are not available. The reader MUST log a warning (once per slice). Records without read names are stored with empty qnames. The pileup engine's overlapping-pair dedup MUST skip records with empty or `*` qnames, since mate matching requires real qnames.

r[cram.edge.seq_unknown]
When CRAM flags bit 0x8 is set, the read's sequence is unknown (`*` in SAM). `seq_len` MUST be 0 and no bases are pushed to the bases slab.

r[cram.edge.empty_slice]
Slices with 0 records MUST be silently skipped.

r[cram.edge.position_overflow]
Delta-coded positions can produce negative values if records are not strictly sorted. The decoder MUST handle this gracefully (negative deltas are valid for supplementary alignments within a coordinate-sorted CRAM).

r[cram.edge.unmapped_reads]
Unmapped reads (flag 0x4) may appear in the last containers of the file, in slices with `reference_sequence_id = -1`. Their bases are stored in the `BA` data series (no reference reconstruction). `fetch_into` MUST skip them, same as BAM.

r[cram.edge.long_reads]
PacBio and ONT CRAM files contain reads of 10,000-1,000,000+ bases. Read feature lists can have thousands of entries, and quality/sequence arrays are correspondingly large. Buffer pre-allocation MUST scale with `RL`, not assume short-read sizes.

r[cram.edge.coordinate_clamp]
Query coordinates passed to `fetch_into` are u64, but internal CRAM positions are i32/i64. When converting `u64::MAX` to i64 for overlap filtering, the cast wraps to -1, which silently filters out all records. The reader MUST clamp coordinates: `end.min(i64::MAX as u64) as i64`.

r[cram.slice.ref_bounds_warning+2]
During sequence reconstruction, if the reference sequence is shorter than the positions accessed (i.e., `ref_offset + ref_pos >= reference_seq.len()`), the reader MUST log a warning on the first occurrence per reconstruction and fall back to `b'N'`. This may indicate a reference/CRAM mismatch but can legitimately occur near chromosome ends. The warning MUST be emitted at most once per record to avoid flooding logs.

r[cram.edge.rans_sym_overflow]
The rANS 4x8 frequency table uses run-length encoding where a symbol counter `sym: u8` is incremented for each run element. When `sym == 255`, `sym += 1` overflows. The reader MUST use wrapping arithmetic (`wrapping_add(1)`) for this counter — the overflow is harmless because the loop terminates before using the overflowed value.

r[cram.codec.rans_sym_bounded+2]
Symbol variables in rANS frequency table run-length reading MUST use `u8` type, making overflow past 255 impossible by construction. After a run-length loop, `sym` may have wrapped to 0 via `wrapping_add(1)` when the previous value was 255. The code MUST NOT use the post-loop `sym` value as a table index if it has wrapped to 0. The run-length loop MUST break when `sym` would exceed 255.

r[cram.codec.state_step_safety]
The rANS `state_step` function computes `f * (s >> bits) + (s & mask) - g`. For valid frequency tables produced by conforming CRAM writers, `g <= f * (s >> bits) + (s & mask)` always holds, so the subtraction does not underflow. The implementation MUST use `wrapping_sub` for robustness against malformed data in release mode, with a `debug_assert!` to catch violations during development.

r[cram.codec.normalize_checked]
Frequency normalization (`normalize_frequencies`) sums all 256 frequency entries and doubles the sum in a loop. Both the sum and the doubling MUST use checked arithmetic to detect overflow from malformed frequency tables. Overflow MUST produce an error, not silent wraparound.

r[cram.codec.uint7_bounded]
The uint7 variable-length integer decoding loop MUST be bounded to at most 5 iterations. A u32 can hold at most 32 bits; 5 iterations of 7-bit groups produce 35 bits, which is the maximum meaningful input. More than 5 continuation bytes indicate malformed data and MUST produce an error.

r[cram.slice.validated_lengths]
Length fields decoded from ITF8 (e.g., `num_content_ids`, `num_blocks`, `alignment_span`) may be negative when interpreted as i32. Before using such values as `usize` for allocation or iteration, the reader MUST validate they are non-negative via `i32::try_from` or equivalent. Negative values MUST produce an error, not wrap to huge `usize` values causing OOM.

r[cram.tok3.dz_len_reader]
The `TokenReader::get()` and `get_mut()` methods MUST return the same reader for every `TokenType` variant. In particular, `DZLen` MUST map to `dz_len_reader` in both methods. A mismatch causes the dup-copy path (which uses `get()`) to read from the wrong stream.

r[cram.tok3.name_count_limit]
The `name_count` field in tok3 headers comes from untrusted data. Allocation based on `name_count` MUST be bounded to a reasonable limit (e.g., 10,000,000 names). The per-name token vector MUST be grown dynamically rather than pre-allocated to a fixed size, since the number of token positions varies per name.

## Performance considerations

CRAM decoding is expected to be 2-4× slower than BAM due to codec complexity and reference reconstruction. The I/O savings (40-60% smaller files) often compensate on network filesystems.

r[cram.perf.slice_granularity]
CRAM random access is slice-granular, not record-granular. A region query may decompress an entire slice (typically 10,000 records) to extract a few records at the edges. This is inherently coarser than BAM's record-level seeking.

r[cram.perf.reference_caching]
Reference sequence lookups happen per-slice (to reconstruct all records in the slice). The reader SHOULD cache the most recently used reference region to avoid repeated FASTA lookups for consecutive slices on the same contig. This cache MUST be per-fork, not shared across threads.

r[cram.perf.codec_overhead]
rANS and arithmetic decoders have higher per-byte CPU cost than zlib. The reader SHOULD pre-allocate decode buffers and reuse them across blocks within a slice.
