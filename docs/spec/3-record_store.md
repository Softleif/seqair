# Record Store

> **Sources:** The slab-based RecordStore design is seqair-specific. The data it stores is defined by [SAM1] §4.2 (BAM record fields, CIGAR packed u32 format, quality Phred encoding, aux tag binary format). See [References](./99-references.md).

## Background

When loading a region from any alignment format (BAM, SAM.gz, or CRAM), the reader decodes hundreds to thousands of records. Each record has fixed-size fields (position, flags, quality) and variable-length fields (read name, CIGAR, sequence, quality scores, auxiliary tags). The original design wrapped each record in `Rc<BamRecord>` with 5 separate heap allocations for the variable-length fields, plus the Rc allocation itself — 6 allocations per record.

The record store replaces this with a **slab-based design**: all variable-length data is packed into contiguous buffers (slabs), and fixed fields are stored in a compact struct that references the slabs by offset. This eliminates all per-record heap allocation.

## Layout

The store contains four vectors:

- **Record table** (`Vec<SlimRecord>`): compact fixed-size structs with alignment fields and offsets into the slabs.
- **Name slab** (`Vec<u8>`): all read names (qnames) packed sequentially. Separated because qnames are only accessed during dedup mate detection — a linear scan of a compact, contiguous buffer.
- **Bases slab** (`Vec<Base>`): decoded base sequences stored as `Base` enum values (A/C/G/T/Unknown). For BAM, decoded from 4-bit nibbles via SIMD at push time. For SAM/CRAM, provided directly as `&[Base]` by the reader.
- **Data slab** (`Vec<u8>`): packed per-record as `[cigar | qual | aux]` — CIGAR in BAM packed u32 format, quality as raw Phred bytes, aux tags in BAM binary format.

```
Record table:  [SlimRecord₀] [SlimRecord₁] [SlimRecord₂] ...
Name slab:     [qname₀|qname₁|qname₂|...]
Bases slab:    [bases₀|bases₁|bases₂|...]
Data slab:     [cigar₀|qual₀|aux₀|cigar₁|qual₁|aux₁|...]
```

## Capacity estimation

r[record_store.capacity]
On construction, the store MUST pre-allocate all three vectors with estimated capacities to avoid reallocation during `push_raw`. Estimates SHOULD be based on the number of compressed bytes loaded for the region: a typical BAM compression ratio of ~3:1, an average record size of ~400 bytes, ~25 bytes per read name, and ~375 bytes per record of non-name data. The store MUST accept a byte-count hint for this purpose.

## Decoding into slabs

r[record_store.push_raw+2]
`push_raw(raw_bytes)` MUST decode a BAM record's fixed fields and append its variable-length data directly into the slabs, without allocating intermediate `Box<[u8]>` or `Vec<u8>` per field. The 4-bit packed sequence MUST be decoded (via SIMD when available) to `Base` enum values in the bases slab.

r[record_store.checked_offsets]
Slab offsets (name_off, bases_off, data_off) MUST be checked for u32 overflow before storing in both `push_raw` and `push_fields`. If any slab exceeds `u32::MAX` bytes, the store MUST return an error rather than silently truncating the offset.

r[record_store.push_fields]
`push_fields(...)` MUST accept pre-parsed record fields for SAM and CRAM readers: pos, end_pos, flags, mapq, matching_bases, indel_bases, qname bytes, CIGAR as packed BAM-format u32 ops, sequence as `&[Base]`, quality bytes, and aux tag bytes in BAM binary format. This avoids encoding to BAM binary only to immediately decode it again. SAM and CRAM parsers convert their native representations to these types before pushing.

## Field access

r[record_store.field_access]
The store MUST provide methods to access each variable-length field for a given record index: `qname(idx)`, `cigar(idx)`, `seq(idx)`, `seq_at(idx, pos)`, `qual(idx)`, `aux(idx)`. These return borrowed slices into the slabs.

## Region lifecycle

r[record_store.clear+2]
Clearing the store for a new region MUST retain the allocated capacity of all four vectors, so that subsequent regions reuse the same memory without reallocation.

r[record_store.region_scoped]
All records in the store MUST remain valid and accessible until the store is cleared.

## Integration

r[record_store.no_rc]
The pileup engine MUST NOT use `Rc` for record sharing. `PileupAlignment` carries pre-extracted flat fields (base, qual, mapq, flags, etc.), and `ActiveRecord` stores a record index into the store. The store is either owned by or borrowed by the engine for the duration of iteration.
