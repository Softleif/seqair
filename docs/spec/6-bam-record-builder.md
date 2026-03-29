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

r[owned_record.fields]
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

r[owned_record.cigar_op]
CIGAR operations MUST be represented as a typed `CigarOp` struct with an `op: CigarOpType` field (using the existing `CigarOpType` enum from [cigar.md](1-3-cigar.md)) and a `len: u32` field. Conversion to/from the BAM packed u32 format (`len << 4 | op_code`) MUST be provided.

## Construction

r[owned_record.builder]
`BamRecord` MUST provide a builder API for constructing records from individual fields. The builder MUST require at minimum `ref_id`, `pos`, `qname`, `cigar`, `seq`, and `qual`. Optional fields (`mapq`, `flags`, `next_ref_id`, `next_pos`, `template_len`, `aux`) MUST default to sensible values (0, 4 for unmapped flags, -1 for unset reference IDs).

r[owned_record.from_record_store]
`RecordStore` MUST provide a `to_owned_record(idx: u32) -> BamRecord` method that extracts a complete owned copy of a record from the slabs. CIGAR bytes MUST be converted to `Vec<CigarOp>`, sequence MUST be cloned from the bases slab, and auxiliary data MUST be copied from the data slab.

r[owned_record.from_record_store_range]
For bulk extraction (e.g. all reads overlapping a realignment window), the store SHOULD provide a method to extract multiple records efficiently, reusing intermediate buffers.

## Modification

> *The following rules support in-memory record manipulation for realignment and BAM rewriting. In TAPS methylation analysis, the BAM rewrite pipeline modifies SEQ (T->C or A->G to reverse methylation evidence) and adds auxiliary tags (MM/ML for SAM 4.5 modification tags, or XR/XG/XM for legacy Bismark format). Realignment additionally requires changing the CIGAR and possibly the position.*

r[owned_record.set_alignment]
`BamRecord` MUST support replacing the alignment by setting a new `pos` and `cigar` simultaneously. The implementation MUST validate that the CIGAR's query-consuming length equals `seq.len()`, returning an error on mismatch.

r[owned_record.set_seq]
`BamRecord` MUST support replacing the sequence. The new sequence length MUST match the CIGAR's query-consuming length. A `seq_mut() -> &mut [Base]` accessor MUST also be provided for in-place base modification (e.g. T->C rewriting in the TAPS BAM rewrite pipeline).

r[owned_record.set_qual]
`BamRecord` MUST support replacing quality scores. The new array length MUST equal `seq.len()`.

r[owned_record.aligned_pairs]
`BamRecord` MUST provide an `aligned_pairs()` iterator yielding `(Option<usize>, Option<i64>)` tuples: `(query_pos, ref_pos)`. For matches/mismatches, both are `Some`. For insertions, `ref_pos` is `None`. For deletions, `query_pos` is `None`. This replaces `rust_htslib::bam::ext::BamRecordExtensions::aligned_pairs_full()` which the BAM rewrite pipeline (`src/bam.rs`) uses to map read positions to reference positions for methylation tag generation.

r[owned_record.end_pos]
`BamRecord` MUST provide `end_pos() -> i64` computing the rightmost reference position from `pos` + CIGAR reference-consuming operations. This MUST be recomputed from the current CIGAR (not cached from a stale value).

r[owned_record.bin]
`BamRecord` MUST provide `bin() -> u16` computing the BAM bin value from `pos` and `end_pos` using the `reg2bin` algorithm from [SAM1] §5.3.

## Serialization

> *[SAM1] §4.2 — BAM record binary layout*

r[owned_record.to_bam_bytes]
`BamRecord` MUST serialize to BAM binary format into a caller-provided `&mut Vec<u8>`. The layout MUST follow [SAM1] §4.2: 4-byte block_size, 32-byte fixed fields, NUL-terminated qname, packed CIGAR (u32 per op), 4-bit packed sequence (via `seq::encode_seq`), quality scores, and raw auxiliary bytes.

r[owned_record.to_bam_bin_field]
The BAM `bin` field in the serialized 32-byte header MUST be computed from the current `pos` and `end_pos`, not stored as a field. This ensures the bin is always consistent with the current alignment, even after modification.

## Re-insertion into RecordStore

> *For realignment workflows, modified records need to be fed back into the pileup engine for re-calling. This requires converting the owned record back into the slab-based format.*

r[owned_record.push_owned]
`RecordStore` MUST provide a `push_owned(record: &BamRecord) -> Result<u32>` method that inserts an owned record into the slabs. This MUST encode the sequence from `Vec<Base>` into the bases slab, pack CIGAR ops into u32 format for the data slab, and copy qname/qual/aux into their respective slabs. The returned index MUST be valid for subsequent `record()`, `seq()`, `qual()`, etc. lookups.

r[owned_record.roundtrip]
A record extracted via `to_owned_record(idx)` and re-inserted via `push_owned()` (without modification) MUST produce a store entry with field values identical to the original. This MUST be verified by tests comparing all accessible fields.

## Auxiliary tag data

> *[SAM1] §4.2.5 — tag byte layout: 2-byte tag name, 1-byte type code, variable-length value*

r[owned_record.aux_data]
`AuxData` MUST store auxiliary tags in BAM binary format (packed bytes). It MUST support:
- `get(tag: [u8; 2]) -> Option<AuxValue>` — read a tag by name (delegating to existing `aux::find_tag`)
- `set_string(tag, value)` — add or replace a Z-type string tag
- `set_int(tag, value)` — add or replace an integer tag (auto-selecting smallest type code: c/C/s/S/i/I)
- `set_float(tag, value)` — add or replace an f-type float tag
- `set_array_u8(tag, values)` — add or replace a B:C array tag
- `remove(tag)` — remove a tag if present
- `as_bytes() -> &[u8]` — raw bytes for serialization

r[owned_record.aux_replace_semantics]
When `set_*` is called for a tag that already exists: if the new value has the same byte length as the old value, the replacement SHOULD be done in-place. Otherwise, the tag MUST be removed and re-appended. This avoids unnecessary reallocations for the common case of updating MM/ML tags with similar-length values.

r[owned_record.aux_from_slab]
When extracting a record from `RecordStore`, the auxiliary data MUST be copied verbatim from the data slab into `AuxData`. No parsing or re-encoding is needed — the slab already stores tags in BAM binary format.

## Flag convenience methods

> *[SAM1] §1.4 — FLAG bit definitions*

r[owned_record.flag_methods]
`BamRecord` MUST provide the same flag query methods as `SlimRecord`: `is_reverse()`, `is_first_in_template()`, `is_second_in_template()`, `is_unmapped()`. It MUST also provide `set_flags(flags: u16)` for bulk flag updates and individual flag setters where commonly needed (e.g. `set_reverse(bool)`).

## Testing

r[owned_record.test_roundtrip_store]
A record loaded via `push_raw` and extracted via `to_owned_record` MUST have identical field values to the original BAM record. This MUST be tested with records containing: simple CIGARs (match-only), complex CIGARs (insertions, deletions, soft clips), multiple auxiliary tags of different types, reverse-strand reads, and paired reads.

r[owned_record.test_roundtrip_bytes]
A `BamRecord` serialized via `to_bam_bytes` and decoded via `push_raw` MUST produce a store entry with identical field values. This verifies the encoder/decoder symmetry.

r[owned_record.test_aligned_pairs]
`aligned_pairs()` output MUST match `rust_htslib::bam::ext::BamRecordExtensions::aligned_pairs_full()` for the same record. Test cases MUST include: match-only CIGAR, CIGAR with insertions, CIGAR with deletions, CIGAR with soft clips, and CIGAR with mixed operations.

r[owned_record.test_modification]
After `set_alignment()` with a new CIGAR, `end_pos()` and `bin()` MUST reflect the new alignment. `aligned_pairs()` MUST yield positions consistent with the new CIGAR. The old alignment's values MUST NOT leak through.
