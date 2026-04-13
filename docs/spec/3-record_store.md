# Record Store

> **Sources:** The slab-based RecordStore design is seqair-specific. The data it stores is defined by [SAM1] §4.2 (BAM record fields, CIGAR packed u32 format, quality Phred encoding, aux tag binary format). See [References](./99-references.md).

## Background

When loading a region from any alignment format (BAM, SAM.gz, or CRAM), the reader decodes hundreds to thousands of records. Each record has fixed-size fields (position, flags, quality) and variable-length fields (read name, CIGAR, sequence, quality scores, auxiliary tags). The original design wrapped each record in `Rc<BamRecord>` with 5 separate heap allocations for the variable-length fields, plus the Rc allocation itself — 6 allocations per record.

The record store replaces this with a **slab-based design**: all variable-length data is packed into contiguous buffers (slabs), and fixed fields are stored in a compact struct that references the slabs by offset. This eliminates all per-record heap allocation.

## Layout

The store contains five vectors:

- **Record table** (`Vec<SlimRecord>`): compact fixed-size structs with alignment fields and offsets into the slabs.
- **Name slab** (`Vec<u8>`): all read names (qnames) packed sequentially. Separated because qnames are only accessed during dedup mate detection — a linear scan of a compact, contiguous buffer.
- **Bases slab** (`Vec<Base>`): decoded base sequences stored as `Base` enum values (A/C/G/T/Unknown). For BAM, decoded from 4-bit nibbles via SIMD at push time. For SAM/CRAM, provided directly as `&[Base]` by the reader.
- **Cigar slab** (`Vec<u8>`): CIGAR ops in BAM packed u32 format. Separated from the data slab because CIGAR is the field that changes during local realignment — isolating it allows append-only mutation without copying qual/aux.
- **Data slab** (`Vec<u8>`): packed per-record as `[qual | aux]` — quality as raw Phred bytes, aux tags in BAM binary format.

```
Record table:  [SlimRecord₀] [SlimRecord₁] [SlimRecord₂] ...
Name slab:     [qname₀|qname₁|qname₂|...]
Bases slab:    [bases₀|bases₁|bases₂|...]
Cigar slab:    [cigar₀|cigar₁|cigar₂|...]
Data slab:     [qual₀|aux₀|qual₁|aux₁|...]
```

## Record fields

r[record_store.slim_record_fields]
Each `SlimRecord` MUST store the following fixed fields: `tid` (reference ID), `pos`, `end_pos`, `flags`, `mapq`, `n_cigar_ops`, `seq_len`, `matching_bases`, `indel_bases`, `next_pos` (mate position), `template_len` (signed insert size). The `tid`, `next_pos`, and `template_len` fields are required for BAM write-back after local realignment.

## Capacity estimation

r[record_store.capacity]
On construction, the store MUST pre-allocate all five vectors with estimated capacities to avoid reallocation during `push_raw`. Estimates SHOULD be based on the number of compressed bytes loaded for the region: a typical BAM compression ratio of ~3:1, an average record size of ~400 bytes, ~25 bytes per read name, ~20 bytes per record of CIGAR data, and ~355 bytes per record of non-name/non-cigar data. The store MUST accept a byte-count hint for this purpose.

## Decoding into slabs

r[record_store.push_raw+2]
`push_raw(raw_bytes)` MUST decode a BAM record's fixed fields and append its variable-length data directly into the slabs, without allocating intermediate `Box<[u8]>` or `Vec<u8>` per field. The 4-bit packed sequence MUST be decoded (via SIMD when available) to `Base` enum values in the bases slab.

r[record_store.checked_offsets]
Slab offsets (name_off, bases_off, cigar_off, data_off) MUST be checked for u32 overflow before storing in both `push_raw` and `push_fields`. If any slab exceeds `u32::MAX` bytes, the store MUST return an error rather than silently truncating the offset.

r[record_store.push_fields]
`push_fields(...)` MUST accept pre-parsed record fields for SAM and CRAM readers: pos, end_pos, flags, mapq, matching_bases, indel_bases, qname bytes, CIGAR as packed BAM-format u32 ops, sequence as `&[Base]`, quality bytes, and aux tag bytes in BAM binary format. This avoids encoding to BAM binary only to immediately decode it again. SAM and CRAM parsers convert their native representations to these types before pushing.

## Field access

r[record_store.field_access]
The store MUST provide methods to access each variable-length field for a given record index: `qname(idx)`, `cigar(idx)`, `seq(idx)`, `seq_at(idx, pos)`, `qual(idx)`, `aux(idx)`. These return borrowed slices into the slabs.

## Region lifecycle

r[record_store.clear+2]
Clearing the store for a new region MUST retain the allocated capacity of all five vectors, so that subsequent regions reuse the same memory without reallocation.

r[record_store.region_scoped]
All records in the store MUST remain valid and accessible until the store is cleared.

## Alignment mutation

r[record_store.set_alignment]
`set_alignment(idx, new_pos, new_cigar_packed)` MUST replace a record's alignment by appending the new CIGAR to the end of the cigar slab and updating the record's `cigar_off`, `n_cigar_ops`, `pos`, `end_pos`, `matching_bases`, and `indel_bases`. The old CIGAR bytes MUST be left in place as dead data (append-only mutation). The sequence, quality, aux, and name slabs MUST NOT be modified. After calling `set_alignment` on one or more records, the caller MUST call `sort_by_pos()` before using the store for pileup iteration.

r[record_store.set_alignment.validation]
`set_alignment` MUST validate that the new CIGAR's query-consuming length equals the record's `seq_len` and that the CIGAR op count fits in `u16`. All validation MUST complete before any slab mutation. If validation fails, the method MUST return an error without modifying the store.

r[record_store.set_alignment.dead_data]
Each `set_alignment` call appends new CIGAR bytes to the cigar slab; old bytes become dead data. For a typical realignment pass where ~10% of records are realigned once, the dead-data overhead is ~10% of cigar slab size. Iterative algorithms that call `set_alignment` multiple times per record accumulate proportionally more dead data. Dead data is reclaimed when `clear()` is called (the cigar slab retains capacity but resets its length). No compaction is provided between `clear()` calls.

## Integration

r[record_store.no_rc]
The pileup engine MUST NOT use `Rc` for record sharing. `PileupAlignment` carries pre-extracted flat fields (base, qual, mapq, flags, etc.), and `ActiveRecord` stores a record index into the store. The store is either owned by or borrowed by the engine for the duration of iteration.
