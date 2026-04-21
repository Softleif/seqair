# Record Encoder

> **Sources:** The record encoder provides a unified, format-agnostic API for writing VCF/BCF records with compile-time state enforcement. A single `Writer` type handles all output formats (VCF, VCF.gz, BCF). Field definitions provide a single source of truth for header construction, key resolution, and documentation generation. See also [VCF Header](./5-vcf-header.md) for string dictionary assignment.

## Unified Writer

r[record_encoder.writer]
A single `Writer<W, S>` type MUST support all output formats (VCF text, BGZF-compressed VCF, BCF binary). The output format MUST be selected at construction time via `OutputFormat`. There MUST NOT be separate `BcfWriter` / `VcfWriter` types in the public API.

r[record_encoder.writer_new]
`Writer::new(inner: W, format: OutputFormat)` MUST create a writer in the `Unstarted` state. No header is required at construction time. Index co-production MUST be automatic for compressed formats (`VcfGz` → TBI, `Bcf` → CSI).

r[record_encoder.writer_typestate]
The writer MUST use a typestate pattern: `Writer<W, Unstarted>` only exposes `write_header()`, which consumes the writer and returns `Writer<W, Ready>`. `Writer<W, Ready>` exposes `begin_record()` and `finish()`.

r[record_encoder.write_header]
`write_header(self, header: &VcfHeader)` MUST consume the `Unstarted` writer, write the header to the output, set up the index builder, and return a `Ready` writer. The header MUST NOT be stored — all field keys are pre-resolved at registration time.

r[record_encoder.finish]
`finish(self)` on `Writer<W, Ready>` MUST flush all buffered data, write the BGZF EOF marker (for compressed formats), finalize the index builder, and return `Result<Option<IndexBuilder>>`.

## Field Definitions

r[record_encoder.field_def]
Field definitions MUST carry all metadata needed for header construction, key resolution, and documentation: name (`&'static str`), Number, ValueType, and description (`&'static str`). Field definitions MUST be parameterized by a value type marker to enable type-safe key resolution.

> r[record_encoder.field_def_types]
> Three definition types MUST be provided:
>
> - `InfoFieldDef<V>` — INFO field definition, produces `InfoKey<V>` on registration
> - `FormatFieldDef<V>` — FORMAT field definition, produces `FormatKey<V>` on registration
> - `FilterFieldDef` — FILTER definition, produces `FilterId` on registration

r[record_encoder.field_def_const]
Field definitions MUST be constructible as `const` values so that downstream crates can define them at compile time.

r[record_encoder.field_description]
A `FieldDescription` trait MUST provide read access to name, number, value type, and description. Both `InfoFieldDef<V>` and `FormatFieldDef<V>` MUST implement this trait (via a type-erased form) to enable documentation generation from heterogeneous collections of field definitions.

## Field Registration

r[record_encoder.register]
`VcfHeaderBuilder` MUST provide `register_info`, `register_format`, and `register_filter` methods that accept a field definition, add the corresponding header entry, and return a resolved typed key. Each method is only available in its corresponding typestate phase (`register_filter` in the Filters phase, `register_info` in the Infos phase, `register_format` in the Formats phase). This combines header construction and key resolution in a single step, eliminating the possibility of declaring a header field without resolving its key or vice versa, and enforces BCF string dictionary ordering at compile time.

r[record_encoder.register_contig]
`VcfHeaderBuilder` MUST provide a `register_contig` method (available in the Contigs phase) that adds a contig to the header and returns a `ContigId` carrying both the integer tid and the contig name.

## Field Identifiers and Keys

r[record_encoder.field_id]
`FieldId` MUST carry both a BCF dictionary index (`u32`) and a string name (`SmolStr`). The BCF encoder uses the dictionary index; the VCF text encoder uses the name. Both are resolved at registration time, not per-record.

r[record_encoder.contig_id]
`ContigId` MUST carry both the integer tid (`u32`) and the contig name (`SmolStr`).

r[record_encoder.filter_id]
`FilterId` MUST carry both the BCF dictionary index and the filter name. A `FilterId::PASS` constant MUST be provided with dictionary index 0 and name `"PASS"`.

### Typed Keys

r[record_encoder.typed_keys]
Typed keys MUST encode the VCF value type at the Rust type level via uninhabited marker enums. This prevents passing a float value where an integer key is expected, or using an INFO key in a FORMAT context.

> r[record_encoder.key_types]
> The following key types MUST be provided:
>
> - `InfoKey<Scalar<i32>>` — single integer INFO value
> - `InfoKey<Scalar<f32>>` — single float INFO value
> - `InfoKey<Arr<i32>>` — integer array INFO value
> - `InfoKey<Arr<f32>>` — float array INFO value
> - `InfoKey<Flag>` — flag INFO value (no data)
> - `InfoKey<Str>` — string INFO value
> - `InfoKey<OptArr<i32>>` — integer array with optional (missing) elements
> - `FormatKey<Gt>` — genotype FORMAT value
> - `FormatKey<Scalar<i32>>` — single integer FORMAT value
> - `FormatKey<Scalar<f32>>` — single float FORMAT value

r[record_encoder.key_encode]
Each typed key MUST provide an `encode` method that delegates to the corresponding encoder trait method. INFO keys MUST accept `&mut impl InfoEncoder` and a single value. FORMAT keys MUST accept `&mut impl FormatEncoder` and a slice of values (one per sample):

- `FormatKey<Gt>::encode(&self, enc: &mut impl FormatEncoder, gts: &[Genotype])`
- `FormatKey<Scalar<T>>::encode(&self, enc: &mut impl FormatEncoder, values: &[T])`

## Typestate Record Encoder

r[record_encoder.typestate]
The record encoder MUST use a typestate pattern to enforce the correct calling sequence at compile time. A single `RecordEncoder<'a, S>` type MUST handle both BCF and VCF encoding via internal enum dispatch. The encoder MUST be parameterized by a state type that restricts which methods are available.

> r[record_encoder.typestate_states]
> Three record states MUST be provided:
>
> - `Begun` — record has been started (fixed fields written). Only filter methods are available.
> - `Filtered` — filters have been written. INFO methods and state transition methods are available.
> - `WithSamples` — sample count has been declared. FORMAT methods and `emit()` are available.

r[record_encoder.typestate_transitions]
State transitions MUST consume the encoder by value and return the encoder in the new state, preventing use of the old state:

- `writer.begin_record(contig, pos, alleles, qual)` → `RecordEncoder<Begun>`
- `filter_pass(self)` / `filter_fail(self)` / `no_filter(self)` on `Begun` → `Filtered`
- `begin_samples(self, n)` on `Filtered` → `WithSamples`
- `emit(self)` on `Filtered` (no samples) or `WithSamples` → consumed (borrow released)

r[record_encoder.typestate_must_use]
The encoder types MUST be marked `#[must_use]` to warn at compile time if a record is silently discarded without calling `emit()`.

r[record_encoder.typestate_w_erased]
The `RecordEncoder` type MUST NOT expose the writer's `W: Write` type parameter. The VCF text path MUST use `&mut dyn Write` internally so that the encoder type is `RecordEncoder<'a, S>` with no extra generic parameters.

r[record_encoder.begin_record]
`begin_record()` MUST be a method on `Writer<W, Ready>`. It MUST accept a `ContigId`, 1-based position, `Alleles`, and optional quality. It MUST clear all per-record state and write/buffer the fixed fields. Overflow checks on position and allele counts MUST use checked conversions and return typed errors.

r[record_encoder.filters]
Exactly one filter method MUST be called per record: `filter_pass()` for PASS, `filter_fail(&[&FilterId])` for one or more failed filters, or `no_filter()` for not-applied. All three MUST consume the `Begun` state and return `Filtered`.

## InfoEncoder Trait

r[record_encoder.info_encoder]
An `InfoEncoder` trait MUST be provided for encoding INFO fields. It MUST be object-safe (all methods use `&mut self` and concrete parameter types). It MUST be implemented by `RecordEncoder<'_, Filtered>`.

r[record_encoder.info_methods]
INFO methods (`info_int`, `info_float`, `info_ints`, `info_floats`, `info_flag`, `info_string`, `info_int_opts`) MUST accept a `&FieldId` and the appropriate value. These methods MUST be infallible — they write to in-memory buffers which cannot fail.

r[record_encoder.info_state_queries]
`InfoEncoder` MUST provide `n_allele()` and `n_alt()` methods returning the number of alleles and alternate alleles for the current record.

r[record_encoder.info_dedup]
If the same INFO field (identified by `FieldId`) is encoded more than once within a single record, the encoder SHOULD overwrite the previously-written value. The final output MUST contain at most one entry per distinct INFO field. The replacement field MUST appear after all non-replaced fields (i.e. it moves to the end of the INFO column). The deduplication tracker SHOULD reuse its allocation across records.

## FormatEncoder Trait

r[record_encoder.format_encoder]
A `FormatEncoder` trait MUST be provided for encoding FORMAT fields. It MUST be object-safe. It MUST be implemented by `RecordEncoder<'_, WithSamples>`.

r[record_encoder.format_methods]
FORMAT methods (`format_gt`, `format_int`, `format_float`) MUST accept a `&FieldId` and a **slice of values** — one element per sample. The slice length MUST equal the sample count declared by `begin_samples()`. These methods MUST be infallible. For single-sample records, callers pass a 1-element slice.

r[record_encoder.format_state_queries]
`FormatEncoder` MUST provide `n_allele()` and `n_alt()` methods returning the number of alleles and alternate alleles for the current record.

r[record_encoder.format_dedup]
If the same FORMAT field (identified by `FieldId`) is encoded more than once within a single record, the encoder SHOULD overwrite the previously-written value. The final output MUST contain at most one entry per distinct FORMAT field. The replacement field MUST appear after all non-replaced fields (i.e. it moves to the end of the FORMAT column). The deduplication tracker SHOULD reuse its allocation across records.

## Emit

r[record_encoder.emit]
`emit()` MUST finalize the record and write it to the output. For BCF, this patches the fixed header and flushes to BGZF. For VCF text, this writes accumulated FORMAT fields and flushes the line buffer. `emit()` MUST consume the encoder by value, releasing the writer borrow. `emit()` is the only encoding method that performs I/O and MAY return an error.

r[record_encoder.emit_no_samples]
`emit()` MUST be available on both `Filtered` (for records without samples) and `WithSamples` (for records with samples).

## Custom Type Encoding

r[record_encoder.encode_info_trait]
An `EncodeInfo` trait MUST be provided with an associated `Key` type and an `encode_info(&self, enc: &mut dyn InfoEncoder, key: &Self::Key)` method. This allows domain types to encapsulate their VCF encoding logic, including multi-field expansion (e.g., a strand-specific type producing separate OT and OB fields).

r[record_encoder.encode_format_trait]
An `EncodeFormat` trait MUST be provided with the same pattern as `EncodeInfo`, using `&mut dyn FormatEncoder`. Implementations that produce no output for certain values (e.g., unknown methylation status) simply do not call any encoder methods, and the FORMAT key does not appear.

r[record_encoder.encode_dyn]
`EncodeInfo` and `EncodeFormat` MUST use `&mut dyn InfoEncoder` and `&mut dyn FormatEncoder` respectively (not `&mut impl`) so that implementations do not need to be generic over the encoder type.

## Format-Specific Encoding

r[record_encoder.bcf_encoding]
The BCF arm of the encoder MUST write INFO fields to `shared_buf` and FORMAT fields to `indiv_buf` following all rules from the [BCF Writer](./5-bcf-writer.md) spec: typed value encoding, smallest int type selection, missing sentinels, field-major FORMAT layout, GT binary encoding.

r[record_encoder.vcf_encoding]
The VCF arm of the encoder MUST write tab-delimited text following all rules from the [VCF Text Format](./5-vcf-writer.md) spec: semicolon-separated INFO, colon-separated FORMAT keys and values, percent-encoding, float precision.

r[record_encoder.buffer_reuse]
All internal buffers MUST be reused across records. After warmup, encoding MUST NOT allocate.

## Equivalence

r[record_encoder.vcf_bcf_equivalence]
For the same logical record, encoding as BCF and as VCF text MUST produce output that, when read back, yields identical field values (within floating-point formatting precision).
