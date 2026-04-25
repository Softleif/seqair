# Record Store

> **Sources:** The slab-based RecordStore design is seqair-specific. The data it stores is defined by [SAM1] §4.2 (BAM record fields, CIGAR packed u32 format, quality Phred encoding, aux tag binary format). See [References](./99-references.md).

## Background

When loading a region from any alignment format (BAM, SAM.gz, or CRAM), the reader decodes hundreds to thousands of records. Each record has fixed-size fields (position, flags, quality) and variable-length fields (read name, CIGAR, sequence, quality scores, auxiliary tags). The original design wrapped each record in `Rc<BamRecord>` with 5 separate heap allocations for the variable-length fields, plus the Rc allocation itself — 6 allocations per record.

The record store replaces this with a **slab-based design**: all variable-length data is packed into contiguous buffers (slabs), and fixed fields are stored in a compact struct that references the slabs by offset. This eliminates all per-record heap allocation.

## Layout

The store contains six vectors:

- **Record table** (`Vec<SlimRecord>`): compact fixed-size structs with alignment fields and offsets into the slabs.
- **Name slab** (`Vec<u8>`): all read names (qnames) packed sequentially. Separated because qnames are only accessed during dedup mate detection — a linear scan of a compact, contiguous buffer.
- **Bases slab** (`Vec<Base>`): decoded base sequences stored as `Base` enum values (A/C/G/T/Unknown). For BAM, decoded from 4-bit nibbles via SIMD at push time. For SAM/CRAM, provided directly as `&[Base]` by the reader.
- **Cigar slab** (`Vec<u8>`): CIGAR ops in BAM packed u32 format. Separated because CIGAR is the field that changes during local realignment — isolating it allows append-only mutation without copying other per-record data.
- **Qual slab** (`Vec<u8>`): per-record quality as raw Phred bytes. Separated from aux so that the pileup hot path (which reads qual but rarely aux) scans a dense, contiguous buffer without aux bytes polluting the cache.
- **Aux slab** (`Vec<u8>`): per-record auxiliary tag bytes in BAM binary format.

```
Record table:  [SlimRecord₀] [SlimRecord₁] [SlimRecord₂] ...
Name slab:     [qname₀|qname₁|qname₂|...]
Bases slab:    [bases₀|bases₁|bases₂|...]
Cigar slab:    [cigar₀|cigar₁|cigar₂|...]
Qual slab:     [qual₀|qual₁|qual₂|...]
Aux slab:      [aux₀|aux₁|aux₂|...]
```

## Record fields

r[record_store.slim_record_fields]
Each `SlimRecord` MUST store the following fixed fields: `tid` (reference ID), `pos`, `end_pos`, `flags`, `mapq`, `n_cigar_ops`, `seq_len`, `matching_bases`, `indel_bases`, `next_ref_id` (mate reference ID), `next_pos` (mate position), `template_len` (signed insert size). The `tid`, `next_ref_id`, `next_pos`, and `template_len` fields are required for BAM write-back after local realignment.

## Capacity estimation

r[record_store.capacity]
On construction, the store MUST pre-allocate all six vectors with estimated capacities to avoid reallocation during `push_raw`. Estimates SHOULD be based on the number of compressed bytes loaded for the region: a typical BAM compression ratio of ~3:1, an average record size of ~400 bytes, ~25 bytes per read name, ~20 bytes per record of CIGAR data, ~150 bytes per record of quality data, and the remaining bytes budgeted for aux. The store MUST accept a byte-count hint for this purpose.

## Decoding into slabs

r[record_store.push_raw+2]
`push_raw(raw_bytes)` MUST decode a BAM record's fixed fields and append its variable-length data directly into the slabs, without allocating intermediate `Box<[u8]>` or `Vec<u8>` per field. The 4-bit packed sequence MUST be decoded (via SIMD when available) to `Base` enum values in the bases slab.

r[record_store.checked_offsets]
Slab offsets (name_off, bases_off, cigar_off, qual_off, aux_off) MUST be checked for u32 overflow before storing in both `push_raw` and `push_fields`. If any slab exceeds `u32::MAX` bytes, the store MUST return an error rather than silently truncating the offset.

r[record_store.push_fields]
`push_fields(...)` MUST accept pre-parsed record fields for SAM and CRAM readers: pos, end_pos, flags, mapq, matching_bases, indel_bases, qname bytes, CIGAR as packed BAM-format u32 ops, sequence as `&[Base]`, quality bytes, and aux tag bytes in BAM binary format. This avoids encoding to BAM binary only to immediately decode it again. SAM and CRAM parsers convert their native representations to these types before pushing.

## Field access

r[record_store.field_access]
The store MUST provide methods to access each variable-length field for a given record index: `qname(idx)`, `cigar(idx)`, `seq(idx)`, `seq_at(idx, pos)`, `qual(idx)`, `aux(idx)`. These return borrowed slices into the slabs.

## Region lifecycle

r[record_store.clear+2]
Clearing the store for a new region MUST retain the allocated capacity of all six vectors, so that subsequent regions reuse the same memory without reallocation.

r[record_store.region_scoped]
All records in the store MUST remain valid and accessible until the store is cleared.

## Alignment mutation

r[record_store.set_alignment]
`set_alignment(idx, new_pos, new_cigar_packed)` MUST replace a record's alignment by appending the new CIGAR to the end of the cigar slab and updating the record's `cigar_off`, `n_cigar_ops`, `pos`, `end_pos`, `matching_bases`, and `indel_bases`. The old CIGAR bytes MUST be left in place as dead data (append-only mutation). The sequence, quality, aux, and name slabs MUST NOT be modified. After calling `set_alignment` on one or more records, the caller MUST call `sort_by_pos()` before using the store for pileup iteration.

r[record_store.set_alignment.validation]
`set_alignment` MUST validate that the new CIGAR's query-consuming length equals the record's `seq_len` and that the CIGAR op count fits in `u16`. All validation MUST complete before any slab mutation. If validation fails, the method MUST return an error without modifying the store.

r[record_store.set_alignment.dead_data]
Each `set_alignment` call appends new CIGAR bytes to the cigar slab; old bytes become dead data. For a typical realignment pass where ~10% of records are realigned once, the dead-data overhead is ~10% of cigar slab size. Iterative algorithms that call `set_alignment` multiple times per record accumulate proportionally more dead data. Dead data is reclaimed when `clear()` is called (the cigar slab retains capacity but resets its length). No compaction is provided between `clear()` calls.

## Per-record extras and filtering

r[record_store.extras.generic_param]
`RecordStore` MUST accept a type parameter `U` (default `()`) for per-record user data. A `Vec<U>` slab MUST grow alongside the existing slabs. When `U` is `()`, `Vec<()>` is a ZST vector with no heap allocation for elements — the extras slab MUST be zero-cost.

r[record_store.extras.push_unit]
`push_raw` and `push_fields` MUST be available only on `RecordStore<()>`. They MUST push `()` to the extras slab for each record added. Readers always produce `RecordStore<()>`.

r[record_store.customize.trait]
The customize trait MUST be defined as `trait CustomizeRecordStore: Clone { type Extra; fn keep_record(&mut self, rec: &SlimRecord, store: &RecordStore<()>) -> bool { true } fn compute(&mut self, rec: &SlimRecord, store: &RecordStore<()>) -> Self::Extra; }`. The `Clone` bound allows `Readers::fork` to duplicate the customize value into the forked reader. The default `keep_record` returns `true` so simple extras-only customizers do not need to override it. A blanket `impl CustomizeRecordStore for ()` MUST exist with `type Extra = ()` and `compute` returning `()` so the default no-customize case costs zero at runtime (`Vec<()>` is a ZST vector). `RecordStore` MUST NOT expose closure-based filter or extras APIs — the trait is the only API, so customize values are always reusable and clone-forwardable. Both `keep_record` and `compute` receive `&SlimRecord` so callers can use `SlimRecord::seq/qual/cigar/aux` getters directly.

r[record_store.customize.apply]
`RecordStore<()>` MUST provide `apply_customize<E: CustomizeRecordStore>(self, customize: &mut E) -> RecordStore<E::Extra>` that computes one `E::Extra` per record by calling `customize.compute(rec, &self)`. The customize value is reusable: `Readers` holds one and runs it on every region fetch, and may carry state (counters, caches) that persist across regions.

r[record_store.pre_filter.rollback]
`RecordStore<()>::push_raw` and `RecordStore<()>::push_fields` MUST accept `&mut E` (`E: CustomizeRecordStore`) as their last argument and return `Result<Option<u32>, DecodeError>`. After parsing and writing the record's slab data they MUST call `customize.keep_record(rec, &self)`; if it returns `false`, the record just pushed MUST be rolled back so that every slab (names/bases/cigar/qual/aux/records/extras) is byte-identical to its pre-push state. Rollback uses the just-pushed `SlimRecord`'s `*_off` offsets as pre-push slab lengths, since push is append-only. The returned `Option<u32>` is `Some(idx)` when kept, `None` when rejected — rejected records MUST NOT consume an index. Passing `&mut ()` (whose default `keep_record` is `true`) is the no-filter form. A proptest verifies that pushing a mixed accept/reject sequence produces byte-identical state to pushing only the accepted inputs with no filter.

r[record_store.pre_filter.retain]
`RecordStore<U>` MUST provide `retain<F: FnMut(&SlimRecord) -> bool>(&mut self, predicate: F)` that removes records where the predicate returns `false`. Unlike the push-time filter, `retain` runs post-push: slab data for dropped records is left in place as dead bytes until `clear()`. This is a general-purpose post-filter for callers that need to drop records after other mutations (e.g., `set_alignment`); the reader-level `fetch_into_customized` on all three formats (BAM, SAM, CRAM) filters at push time with zero waste and does not use `retain`.

r[record_store.slim_record.field_getters]
`SlimRecord` MUST provide getter methods `seq(&store) -> Result<&[Base]>`, `qual(&store) -> Result<&[BaseQuality]>`, `cigar(&store) -> Result<&[u8]>` (packed BAM CIGAR), and `aux(&store) -> Result<&[u8]>` that take `&RecordStore<U>` for any `U` and return slab-backed slices using the record's offset fields. They MUST validate offsets and return a typed `RecordAccessError` on overflow rather than panic. `extra(&store) -> Result<&U>` follows the same shape but reads the extras slab. These getters are the canonical way for `CustomizeRecordStore` impls to access a record's variable-length data without going through `RecordStore::record(idx)` first.

r[record_store.extras.access]
`RecordStore<U>` MUST provide `extra(idx) -> &U` and `extra_mut(idx) -> &mut U` to access the per-record extra by record index. Access MUST go through the record's `extras_idx` field (indirection into the extras slab), so that extras remain valid after `sort_by_pos` or `dedup` reorder/remove records.

r[record_store.extras.clear]
`clear()` on `RecordStore<U>` MUST also clear the extras slab, retaining its allocated capacity.

r[record_store.extras.strip]
`RecordStore<U>` MUST provide `strip_extras(self) -> RecordStore<()>` that discards the extras slab while preserving all other slab capacity. This is used by `Readers::recover_store` to recycle a `RecordStore<U>` back to `RecordStore<()>`.

r[record_store.extras.sort_dedup_generic]
`sort_by_pos` and `dedup` MUST be available on `RecordStore<U>` for any `U`. Each `SlimRecord` carries an `extras_idx` field that indexes into the extras slab, so reordering or removing records does not invalidate the extras mapping. Dead extras entries from `dedup` are left in place (minor waste, same as dead slab data for names, cigar, etc.).

## Integration

r[record_store.no_rc]
The pileup engine MUST NOT use `Rc` for record sharing. `PileupAlignment` carries pre-extracted flat fields (base, qual, mapq, flags, etc.), and `ActiveRecord` stores a record index into the store. The store is either owned by or borrowed by the engine for the duration of iteration.
