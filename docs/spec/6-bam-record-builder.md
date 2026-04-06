# Owned BAM Record

When a BAM record needs to be constructed from scratch (for writing) or modified in memory (for realignment), the read-only slab-based `SlimRecord`/`RecordStore` architecture is insufficient — it packs variable-length data into shared contiguous buffers with no mechanism for per-record mutation. An owned, mutable record type bridges this gap.

The `BamRecord` is an owned BAM record with individually-allocated fields that supports construction, modification, serialization to BAM binary format, and round-trip conversion to/from `RecordStore` entries. It is the unit of exchange between the read path (RecordStore/pileup) and the write path (BamWriter), and enables in-memory record manipulation for realignment workflows.

> **Sources:** [SAM1] §4.2 "The BAM format" — binary record layout (block_size, refID, pos, bin_mq_nl, flag_nc, l_seq, next_refID, next_pos, tlen, qname, cigar, seq, qual, aux); §1.4 — FLAG bit definitions; §4.2.4 — SEQ 4-bit encoding; §4.2.5 — auxiliary tag encoding. CIGAR operation semantics from [SAM1] §1.4 "CIGAR". See [references.md](99-references.md).

## Background

### Why a separate type from SlimRecord

`SlimRecord` (see [record_store.md](3-record_store.md)) stores fixed fields inline and references variable-length data by offset into shared slabs. This design eliminates per-record allocation and is optimal for the read-only pileup hot path, but it makes mutation impossible: changing one record's CIGAR would require shifting all subsequent records' data in the slab.

`BamRecord` owns its data in separate `Vec`s, trading cache density for mutability. The typical use case is extracting a small number of records (e.g. reads at a realignment candidate site) from the store, modifying them, and either writing them to a BAM file or re-inserting them into a fresh `RecordStore` for re-pileup.

### Realignment context

In TAPS methylation sequencing, C-to-T conversions caused by methylation are indistinguishable from true C>T SNPs at the sequence level. Realignment against local haplotypes can resolve ambiguous alignments at these positions by considering the methylation context. This requires modifying the CIGAR (new alignment path), potentially the position (if the optimal alignment starts elsewhere), and recomputing alignment-derived tags (NM, MD).

The owned record type makes these modifications safe and explicit: callers get a mutable `BamRecord`, change what they need, and the record recomputes derived fields (end position, BAM bin) on serialization.

## Owned record structure

> *[SAM1] §4.2 — fixed 32-byte header fields, variable-length qname/cigar/seq/qual/aux*

r[bam.owned_record.fields]
`BamRecord` MUST store the following fields, all owned and mutable:
- `ref_id: i32` — reference sequence index (-1 for unmapped)
- `pos: i64` — 0-based leftmost mapping position (-1 for unmapped)
- `mapq: u8` — mapping quality
- `flags: u16` — SAM flag bits
- `next_ref_id: i32` — mate's reference sequence index
- `next_pos: i64` — mate's 0-based position
- `template_len: i32` — observed template length
- `qname: Vec<u8>` — query name (without NUL terminator)
- `cigar: Vec<CigarOp>` — typed CIGAR operations
- `seq: Vec<Base>` — sequence as `Base` enum values
- `qual: Vec<u8>` — Phred quality scores
- `aux: AuxData` — auxiliary tags (see aux data section)

r[bam.owned_record.cigar_op]
CIGAR operations MUST be represented as a typed `CigarOp` struct with an `op: CigarOpType` field (using the existing `CigarOpType` enum defined in `r[io.typed_cigar_ops]` in [general.md](0-general.md)) and a `len: u32` field. Conversion to/from the BAM packed u32 format (`len << 4 | op_code`) MUST be provided.

## Construction

> *[SAM1] §4.2 — l_read_name (u8, max 255 including NUL), n_cigar_op (u16, max 65535), l_seq (i32)*

r[bam.owned_record.builder]
`BamRecord` MUST provide a builder API for constructing records from individual fields. The builder MUST require at minimum `ref_id`, `pos`, and `qname`. Fields `cigar`, `seq`, and `qual` MUST default to empty (representing an unmapped read). Optional fields (`mapq`, `flags`, `next_ref_id`, `next_pos`, `template_len`, `aux`) MUST default to zero or -1 for reference IDs. The `flags` field MUST default to 0; callers are responsible for setting appropriate flag bits.

r[bam.owned_record.qname_limit]
The qname (without NUL terminator) MUST NOT exceed 254 bytes. The `l_read_name` field in the BAM 32-byte header is a u8 that stores the length including NUL, so the maximum stored length is 255. The builder and `to_bam_bytes` MUST reject qnames exceeding this limit with a typed error.

r[bam.owned_record.cigar_count_limit]
The CIGAR operation count MUST NOT exceed 65535. The `n_cigar_op` field in the BAM header is a u16. The builder, `set_alignment`, and `to_bam_bytes` MUST reject CIGAR vectors exceeding this limit with a typed error.

r[bam.owned_record.seq_length_limit]
The sequence length MUST NOT exceed `i32::MAX` (2,147,483,647). The `l_seq` field in the BAM header is an i32. The builder, `set_seq`, and `to_bam_bytes` MUST reject sequences exceeding this limit with a typed error.

r[bam.owned_record.from_record_store]
`RecordStore` MUST provide a `to_owned_record(idx: u32) -> BamRecord` method that extracts a complete owned copy of a record from the slabs. CIGAR bytes MUST be converted to `Vec<CigarOp>`, sequence MUST be cloned from the bases slab, and auxiliary data MUST be copied from the data slab.

## Modification

The following rules support in-memory record manipulation for realignment and BAM rewriting. In TAPS methylation analysis, the BAM rewrite pipeline modifies SEQ (T->C or A->G to reverse methylation evidence) and adds auxiliary tags (MM/ML for SAM 4.5 modification tags, or XR/XG/XM for legacy Bismark format). Realignment additionally requires changing the CIGAR and possibly the position.

r[bam.owned_record.set_alignment]
`BamRecord` MUST support replacing the alignment by setting a new `pos` and `cigar` simultaneously. For mapped reads (flags & 0x4 == 0), the implementation MUST validate that the CIGAR's query-consuming length equals `seq.len()`, returning an error on mismatch. For unmapped reads (flags & 0x4 != 0), a zero-operation CIGAR with non-empty `seq` MUST be allowed.

r[bam.owned_record.set_seq]
`BamRecord` MUST support replacing the sequence. For mapped reads, the new sequence length MUST match the CIGAR's query-consuming length. A `seq_mut() -> &mut [Base]` accessor MUST also be provided for in-place base modification (e.g. T->C rewriting in the TAPS BAM rewrite pipeline).

r[bam.owned_record.set_qual]
`BamRecord` MUST support replacing quality scores. The new array length MUST equal `seq.len()`.

r[bam.owned_record.aligned_pairs]
`BamRecord` MUST provide an `aligned_pairs()` iterator yielding `(Option<usize>, Option<i64>)` tuples: `(query_pos, ref_pos)`. The iterator handles each CIGAR operation as follows:

- **M/=/X** (alignment match, sequence match, mismatch): both values are `Some`, one tuple per base.
- **I** (insertion to reference): `ref_pos = None`, one tuple per inserted base.
- **D** (deletion from reference): `query_pos = None`, one tuple per deleted base.
- **N** (reference skip / intron): `query_pos = None`, one tuple per skipped reference base. N consumes reference but not query, like D. This is important for RNA-seq BAMs where introns appear as N operations spanning thousands of reference bases.
- **S** (soft clip): MUST be skipped entirely. Soft-clipped bases are present in the sequence but not aligned to the reference.
- **H** (hard clip): MUST be skipped entirely. Hard-clipped bases are not present in the sequence.
- **P** (padding): MUST be skipped entirely. Padding is used for padded-reference alignment and consumes neither query nor reference.

For unmapped reads (empty CIGAR), the iterator MUST be empty.

r[bam.owned_record.end_pos]
`BamRecord` MUST provide `end_pos() -> i64` computing the rightmost reference position from `pos` + CIGAR reference-consuming operations. This uses the same logic as `r[bam.record.end_pos]` in [record.md](2-bam-3-2-record.md), but MUST be recomputed from the current CIGAR state (not cached from a stale value).

r[bam.owned_record.bin]
`BamRecord` MUST provide `bin() -> u16` computing the BAM bin value from `pos` and `end_pos` using the `reg2bin` algorithm from [SAM1] §5.3. This method is specific to the BAI binning scheme (min_shift=14, depth=5), where the maximum non-pseudo bin ID is 37449 — well within u16 range. CSI indexes with larger depths can produce bin IDs exceeding u16; those are handled by the IndexBuilder directly and do not flow through `BamRecord::bin()`. In the serialized 32-byte header, this value is packed into the `bin_mq_nl` field as `bin << 16 | mapq << 8 | l_read_name` — the serialization logic in `to_bam_bytes` handles this packing.

## Serialization

> *[SAM1] §4.2 — BAM record binary layout: block_size (i32), 32-byte fixed header, variable fields*

r[bam.owned_record.to_bam_bytes]
`BamRecord` MUST serialize to BAM binary format by appending into a caller-provided `&mut Vec<u8>` (the method appends; clearing the buffer is the caller's responsibility — see `r[bam_writer.reuse_buffers]` in [bam-writer.md](6-bam-writer.md)). The layout MUST follow [SAM1] §4.2: 32-byte fixed fields (refID, pos, bin_mq_nl, flag_nc, l_seq, next_refID, next_pos, tlen), NUL-terminated qname, packed CIGAR (u32 per op), 4-bit packed sequence (via `seq::encode_seq` from [seq_codec.md](2-bam-3-4-seq_codec.md)), quality scores, and raw auxiliary bytes. The method MUST NOT include the 4-byte `block_size` prefix — that is the caller's (writer's) responsibility, since the caller needs to know the total byte count to write the prefix.

r[bam.owned_record.to_bam_field_overflow]
The BAM binary format stores `pos` and `next_pos` as i32, and `ref_id` and `next_ref_id` as i32. The `BamRecord` stores `pos` and `next_pos` as i64 for internal flexibility (e.g., arithmetic without overflow risk). `to_bam_bytes` MUST validate that `pos` is in `[-1, i32::MAX]` and `next_pos` is in `[-1, i32::MAX]`, returning a typed error on overflow. The value -1 is the canonical sentinel for unmapped/unavailable. This validation catches bugs where CIGAR-based position calculations produce values exceeding the BAM coordinate space.

r[bam.owned_record.to_bam_bin_field]
The BAM `bin` value in the serialized `bin_mq_nl` field MUST be computed from the current `pos` and `end_pos` at serialization time, not stored as a persistent field. This ensures the bin is always consistent with the current alignment, even after modification.

## Re-insertion into RecordStore

For realignment workflows, modified records need to be fed back into the pileup engine for re-calling. This requires converting the owned record back into the slab-based format.

r[bam.owned_record.push_owned]
`RecordStore` MUST provide a `push_owned(record: &BamRecord) -> Result<u32>` method that inserts an owned record into the slabs. This MUST encode the sequence from `Vec<Base>` into the bases slab, pack CIGAR ops into u32 format for the data slab, and copy qname/qual/aux into their respective slabs. The returned index MUST be valid for subsequent `record()`, `seq()`, `qual()`, etc. lookups.

r[bam.owned_record.roundtrip]
A record extracted via `to_owned_record(idx)` and re-inserted via `push_owned()` (without modification) MUST produce a store entry with field values identical to the original. This MUST be verified by tests comparing all accessible fields.

## Auxiliary tag data

> *[SAM1] §4.2.5 — tag byte layout: 2-byte tag name, 1-byte type code, variable-length value. Tag types: A (char), c/C/s/S/i/I (integers), f (float), d (double), Z (string), H (hex string), B (typed array with subtype byte + 4-byte count + values).*

r[bam.owned_record.aux_data]
`AuxData` MUST store auxiliary tags in BAM binary format (packed bytes). Tag lookup via `get()` MUST delegate to the existing `aux::find_tag` implementation (see `r[bam.record.aux_parse]` in [record.md](2-bam-3-2-record.md)). Mutation methods MUST be provided:
- `set_string(tag, value)` — add or replace a Z-type string tag
- `set_int(tag, value: i64)` — add or replace an integer tag
- `set_float(tag, value)` — add or replace an f-type float tag
- `set_array_u8(tag, values)` — add or replace a B:C array tag
- `remove(tag)` — remove a tag if present
- `as_bytes() -> &[u8]` — raw bytes for serialization

r[bam.owned_record.aux_uniqueness]
Tag names MUST be unique within a record (per [SAM1] §1.5). `set_*` methods MUST enforce this: if a tag with the given name already exists, it MUST be replaced, not duplicated. `get()` MUST be unambiguous.

r[bam.owned_record.aux_int_encoding]
`set_int(tag, value: i64)` MUST auto-select the smallest BAM integer type that fits the value. For values where both signed and unsigned types could apply (e.g. 42 fits in both i8 and u8), the unsigned type MUST be preferred — this matches htslib behavior and produces more compact encodings for the common case of non-negative values. The selection order is:

1. Non-negative values: `C` (u8, 0..=255), `S` (u16, 256..=65535), `I` (u32, 65536..=2^32-1)
2. Negative values: `c` (i8, -128..=-1), `s` (i16, -32768..=-129), `i` (i32, -2^31..=-32769)

Values outside the union of i32 and u32 ranges MUST return a typed error — BAM auxiliary tags have no 64-bit integer type per [SAM1] §4.2.5.

r[bam.owned_record.aux_array_encoding]
`set_array_u8(tag, values)` MUST encode the tag in BAM B-type array format: 2-byte tag name, type byte `B`, subtype byte `C`, 4-byte little-endian element count, followed by the raw u8 values. This is used for the ML (modification likelihood) tag in SAM 4.5 methylation annotations.

r[bam.owned_record.aux_replace_semantics]
When `set_*` is called for a tag that already exists, the tag MUST be removed and re-appended with the new value. Implementations MAY optimize the case where the new encoded byte length (including tag name and type byte) matches the old length by replacing in-place, but this is not required.

r[bam.owned_record.aux_from_slab]
When extracting a record from `RecordStore`, the auxiliary data MUST be copied verbatim from the data slab into `AuxData`. No parsing or re-encoding is needed — the slab already stores tags in BAM binary format.

## Flag convenience methods

`BamRecord` MUST implement the same flag interface as `SlimRecord`, as specified by `r[io.typed_flags]` in [general.md](0-general.md) and the flag rules in [record.md](2-bam-3-2-record.md) (`r[bam.record.flag_reverse]`, `r[bam.record.flag_first]`, `r[bam.record.flag_second]`, `r[bam.record.flag_unmapped]`).

r[bam.owned_record.flag_methods]
In addition to the read-only flag predicates, `BamRecord` MUST provide `set_flags(flags: u16)` for bulk flag updates.

## Testing

r[bam.owned_record.test_roundtrip_store]
A record loaded via `push_raw` and extracted via `to_owned_record` MUST have identical field values to the original BAM record. This MUST be tested with records containing: simple CIGARs (match-only), complex CIGARs (insertions, deletions, soft clips), multiple auxiliary tags of different types, reverse-strand reads, and paired reads.

r[bam.owned_record.test_roundtrip_bytes]
A `BamRecord` serialized via `to_bam_bytes` (with a prepended block_size) and decoded via `push_raw` MUST produce a store entry with identical field values. This verifies the encoder/decoder symmetry.

r[bam.owned_record.test_aligned_pairs]
`aligned_pairs()` output MUST match `rust_htslib::bam::ext::BamRecordExtensions::aligned_pairs_full()` for the same record. Test cases MUST include: match-only CIGAR, CIGAR with insertions (one tuple per inserted base), CIGAR with deletions (one tuple per deleted base), CIGAR with soft clips (skipped), CIGAR with mixed operations, and unmapped records (empty iterator).

r[bam.owned_record.test_modification]
After `set_alignment()` with a new CIGAR, `end_pos()` and `bin()` MUST reflect the new alignment. `aligned_pairs()` MUST yield positions consistent with the new CIGAR. The old alignment's values MUST NOT leak through.
