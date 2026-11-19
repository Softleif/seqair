# Record Encoder

> **Sources:** The record encoder provides a unified, format-agnostic API for writing VCF/BCF records with compile-time state enforcement. A single `Writer` type handles all output formats (VCF, VCF.gz, BCF). Field definitions provide a single source of truth for header construction, key resolution, and documentation generation. See also [VCF Header](./5-vcf-header.md) for string dictionary assignment.

## Unified Writer

r[record_encoder.writer]
A single `Writer<W, S>` type MUST support all output formats (VCF text, BGZF-compressed VCF, BCF binary). The output format MUST be selected at construction time via `OutputFormat`. There MUST NOT be separate `BcfWriter` / `VcfWriter` types in the public API.

r[record_encoder.writer_new+2]
`Writer::new(inner: W, format: OutputFormat)` MUST create a writer in the `Unstarted` state. No header is required at construction time. Index co-production MUST be automatic for compressed formats: a CSI for `Bcf`, and for `VcfGz` the index `index_kind` selects, CSI by default (`r[record_encoder.index_kind]`).

r[record_encoder.index_kind]
`Writer<W, Unstarted>::index_kind(kind: IndexKind)` MUST select the coordinate index a `VcfGz` writer co-produces — `IndexKind::Csi` (the default) or `IndexKind::Tbi` — and return the writer, still `Unstarted`, like `compression_level`. CSI is the default because it addresses every contig length (`r[index_builder.csi_depth]`) while TBI stops at 2^29; TBI is there for tools that read nothing else. BCF has only CSI (htslib refuses to build a TBI for BCF), so `write_header` on a `Bcf` writer asked for `Tbi` MUST fail with `VcfError::TbiForBcf` before writing anything. Plain VCF is not indexed and MUST ignore the setting. A TBI builder is fixed at `min_shift` 14, depth 5; a record that ends past 2^29 MUST fail its `emit` with `IndexError::BeyondIndexRange` (`r[index_builder.range_check]`), as htslib fails it, rather than be indexed where no query finds it. `IndexKind::suffix()` (a `const fn`: `"csi"` or `"tbi"`) names the index file, `<data path>.<suffix>`.

r[record_encoder.compression_threads]
`Writer<W, Unstarted>::compression_threads(n)` MUST switch the `VcfGz`/`Bcf` BGZF stream to `n` compression worker threads (`n = 0`, the default, compresses on the calling thread) and return the writer, still `Unstarted`. It cannot fail: the threads start in `write_header`, which reports a failure to start them — as `BamWriterBuilder::build` does for BAM. Plain VCF is not compressed and MUST ignore it. Output and coordinate index MUST be byte-identical to the single-threaded writer's (see `r[bgzf.writer.parallel]`).

r[record_encoder.compression_level]
`Writer<W, Unstarted>::compression_level(level)` MUST set the `VcfGz`/`Bcf` BGZF compression level (a `CompressionLevel`, `r[bgzf.writer.compression_level]`; default 6, htslib's default) and return the writer, still `Unstarted`, taking effect from the header on, like `BamWriterBuilder::compression_level`. It combines with `compression_threads` in either order. Plain VCF MUST ignore it.

r[record_encoder.writer_typestate]
The writer MUST use a typestate pattern: `Writer<W, Unstarted>` only exposes `write_header()`, which consumes the writer and returns `Writer<W, Ready>`. `Writer<W, Ready>` exposes `begin_record()` and `finish()`.

r[record_encoder.write_header+2]
`write_header(self, header: &VcfHeader)` MUST consume the `Unstarted` writer, write the header to the output, set up the index builder, and return a `Ready` writer. Encoding MUST NOT consult the header — all field keys are pre-resolved at registration time; the writer keeps one copy of it only to build a passthrough copy's field map (`r[vcf_passthrough.implicit_map+2]`).

r[record_encoder.finish+2]
`finish(self)` on `Writer<W, Ready>` MUST flush all buffered data, write the BGZF EOF marker (for compressed formats), finalize the index builder, and return `Result<(W, Option<CoordinateIndex>), VcfError>`: `Some` for `VcfGz` and `Bcf`, `None` for plain VCF. A `CoordinateIndex` knows its `kind()` and so its `suffix()`, and `write`/`write_to_path` serialize it in that format: TBI, CSI with tabix metadata for VCF.gz, plain CSI for BCF.

## Field Definitions

r[record_encoder.field_def+5]
Field definitions MUST carry all metadata needed for header construction, key resolution, and documentation: name (`&'static str`), Number, and description (`&'static str`). Field definitions MUST be parameterized by a value type marker to enable type-safe key resolution, and their ValueType MUST come from that marker (`FieldKind::VALUE_TYPE`, `r[vcf_values.field_kind]`), never from a separate argument that could disagree with it. The `Number` MUST come from the marker too wherever the marker fixes it: `Scalar<T>` and `Gt` declare `Number=1` and `Flag` declares `Number=0`, and their definitions are built with `new(name, description)`, which takes no `Number`. `Arr<T>`, `OptArr<T>` and `Str` definitions are built with `with_number(name, number, description)`; neither constructor exists for the other shape, so a scalar with another `Number` does not compile. The constructor MUST be selected by sealed marker traits (`InfoFixedNumber`, `FormatFixedNumber`, `ChosenNumber`) on one generic `impl` each, so the marker is inferred from a typed binding (`const DP: InfoFieldDef<Scalar<i32>> = InfoFieldDef::new("DP", "Depth")`); one inherent `new` per marker makes every untyped call ambiguous (E0034) and forced callers to name the marker twice. A string's cardinality is not fixed by its marker: real files declare String fields `Number=1` and `Number=.` (ANN, CSQ, GIAB's `platformnames`). A `Str` value is one string, or a list through `encode_list` (`r[record_encoder.string_list]`); `encode`'s single string is not count-checked (`r[record_encoder.count_check+3]`). A value per ALT is an `Arr` with `Number=A`, even at a biallelic site. The fields are private and read through `FieldDescription`, so a `Number` cannot be changed after construction. There is no INFO `Gt` and no FORMAT `Flag` constructor; registering a FORMAT `Gt` definition under any name but `GT` MUST fail with `KeyError::NotGenotype`.

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

r[record_encoder.field_id+2]
`FieldId` MUST carry the field's BCF dictionary id (a `DictId`, 16 bits: the dictionary holds at most `MAX_DICT_ENTRIES` = 2^16 ids, `r[vcf_header.dict_capacity+2]`, so the index is narrowed once, where the key is resolved), its name (`SmolStr`), its declared `Type` (`String` for a field no header declares) and `Number`, whether it is the FORMAT `GT` field, and its header's `HeaderId`. The BCF encoder uses the dictionary id; the VCF text encoder uses the name; a reader decodes an `Any` value by the declared `Type` and the `GT` flag (`r[vcf_values.dynamic+2]`). All of it is resolved at registration or key resolution, not per record: an access never looks the field up again.

r[record_encoder.contig_id]
`ContigId` MUST carry both the integer tid (`u32`) and the contig name (`SmolStr`).

r[record_encoder.filter_id]
`FilterId` MUST carry both the BCF dictionary index and the filter name. A `FilterId::PASS` constant MUST be provided with dictionary index 0 and name `"PASS"`.

### Typed Keys

r[record_encoder.typed_keys]
Typed keys MUST encode the VCF value type at the Rust type level via uninhabited marker enums. This prevents passing a float value where an integer key is expected, or using an INFO key in a FORMAT context.

> r[record_encoder.key_types+2]
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
> - `FormatKey<Scalar<i32>>` — single integer FORMAT value, a sample of which may be missing (`encode_opt`)
> - `FormatKey<Scalar<f32>>` — single float FORMAT value, a sample of which may be missing (`encode_opt`)
> - `FormatKey<Arr<i32>>` — multi-valued integer FORMAT value (one array per sample)
> - `FormatKey<Arr<f32>>` — multi-valued float FORMAT value (one array per sample)
> - `FormatKey<OptArr<f32>>` — multi-valued float FORMAT value with optional (missing) elements, encoded for a single sample with `encode_single_sample`; a `None` element is the BCF float missing value and `.` in VCF text
> - `FormatKey<Str>` — string FORMAT value (one string per sample)

r[record_encoder.key_encode+2]
Each typed key MUST provide an `encode` method that delegates to the corresponding encoder trait method. INFO keys MUST accept `&mut impl InfoEncoder` and a single value. FORMAT keys MUST accept `&mut impl FormatEncoder` and a slice of per-sample values (one element, or one inner slice for array keys, per sample):

- `FormatKey<Gt>::encode(&self, enc: &mut impl FormatEncoder, gts: &[Genotype])`
- `FormatKey<Scalar<T>>::encode(&self, enc: &mut impl FormatEncoder, values: &[T])`, and `encode_opt(&self, enc, values: &[Option<T>])` (`r[record_encoder.format_scalar_opt]`)
- `FormatKey<Arr<T>>::encode(&self, enc: &mut impl FormatEncoder, per_sample: &[&[T]])`
- `FormatKey<Str>::encode(&self, enc: &mut impl FormatEncoder, values: &[&str])`

## Typestate Record Encoder

r[record_encoder.typestate]
The record encoder MUST use a typestate pattern to enforce the correct calling sequence at compile time. A single `RecordEncoder<'a, S>` type MUST handle both BCF and VCF encoding via internal enum dispatch. The encoder MUST be parameterized by a state type that restricts which methods are available.

> r[record_encoder.typestate_states+2]
> Four record states MUST be provided:
>
> - `Begun` — record has been started (fixed fields written). Only filter methods are available.
> - `Filtered` — filters have been written. INFO methods and state transition methods are available.
> - `WithSamples` — sample count has been declared. FORMAT methods and `emit()` are available.
> - `Edit` — a copied record (`RecordEdit`, `r[vcf_passthrough.edit+2]`). INFO and FORMAT methods in any order, `remove_info`/`remove_format`, and `emit()` are available.

r[record_encoder.typestate_transitions+4]
State transitions MUST consume the encoder by value and return the encoder in the new state, preventing use of the old state:

- `writer.begin_record(contig, pos, alleles, qual)` → `RecordEncoder<Begun>`
- `ids(self, ids)` on `Begun` → `Result<Begun, VcfError>`, optional, since the ID column is checked (`r[record_encoder.record_id]`)
- `filter_pass(self)` / `no_filter(self)` on `Begun` → `Filtered`; `filter_fail(self, filters)` → `Result<Filtered, VcfError>`, since its handles are checked (`r[record_encoder.foreign_keys]`)
- `begin_samples(self)` on `Filtered` → `WithSamples` (the sample count is the header's)
- `writer.copy_record(record)` / `copy_record_with(map, record)` → `Result<RecordEncoder<Edit>, VcfError>` (`r[vcf_passthrough.copy]`), the whole record written and open for edits
- `emit(self)` on `Filtered` (no samples), `WithSamples` or `Edit` → consumed (borrow released)

r[record_encoder.typestate_must_use]
The encoder types MUST be marked `#[must_use]` to warn at compile time if a record is silently discarded without calling `emit()`.

r[record_encoder.typestate_w_erased]
The `RecordEncoder` type MUST NOT expose the writer's `W: Write` type parameter. The VCF text path MUST use `&mut dyn Write` internally so that the encoder type is `RecordEncoder<'a, S>` with no extra generic parameters.

r[record_encoder.begin_record]
`begin_record()` MUST be a method on `Writer<W, Ready>`. It MUST accept a `ContigId`, 1-based position, `Alleles`, and optional quality. It MUST clear all per-record state and write/buffer the fixed fields. Overflow checks on position and allele counts MUST use checked conversions and return typed errors.

r[record_encoder.record_id]
A record's ID column is `.` unless `ids(self, ids: impl IntoIterator<Item = &str>)` is called on the `Begun` encoder, which writes the given identifiers `;`-separated, in order, in VCF text and BCF alike; an empty iterator writes `.`. Calling it again replaces the column. Each identifier MUST be checked against VCF's ID grammar (§1.6.1: "no whitespace or semicolons permitted, duplicate values not allowed") before anything is written (the encoder is consumed, as by `filter_fail`, so a rejected record is never emitted): an empty identifier, one that contains whitespace or `;`, or one that is `.` (the missing value, which a reader cannot tell from no ID) MUST fail with `VcfError::InvalidRecordId { id, problem: RecordIdProblem }`, and an identifier given twice MUST fail with `VcfError::DuplicateRecordId { id }` (the spec's "duplicate values not allowed"; that no identifier appears in two records is the caller's to keep, since the writer holds one record at a time). It lives on `Begun` because ID precedes FILTER in the record, so the encoder writes it in column order.

r[record_encoder.filters+3]
Exactly one filter method MUST be called per record: `filter_pass()` for PASS, `filter_fail(filters)` for one or more failed filters, or `no_filter()` for not-applied. All three MUST consume the `Begun` state; `filter_pass` and `no_filter` return `Filtered`, `filter_fail` returns `Result<Filtered, VcfError>` and fails with `ForeignKey` for a handle from another header (`r[record_encoder.foreign_keys]`) and with `EmptyFilterList` for no filters: VCF text would write an empty column, BCF an empty vector that reads as `.`, so the two formats would disagree and neither says the record failed.

## InfoEncoder Trait

r[record_encoder.info_encoder+2]
An `InfoEncoder` trait MUST be provided for encoding INFO fields. It MUST be implemented by `RecordEncoder<'_, Filtered>` and `RecordEncoder<'_, Edit>`. Its methods take `&mut self` and concrete parameter types; the typed [`InfoKey`] encode wrappers dispatch through it.

r[record_encoder.info_methods+3]
INFO methods (`info_int`, `info_float`, `info_ints`, `info_floats`, `info_flag`, `info_string`, `info_strings`, `info_int_opts`, `info_float_opts`) MUST accept a `&FieldId` and the appropriate value. They perform no I/O — they write to in-memory buffers, which cannot fail. All of them nevertheless return `Result<(), VcfError>`, because they validate before writing: the key's header (`r[record_encoder.foreign_keys]`), the value count (`r[record_encoder.count_check+3]`) and reserved numbers (`r[record_encoder.reserved_rejection+2]`).

r[record_encoder.string_list]
`InfoKey<Str>::encode_list(&mut enc, &[&str])` MUST write a list of strings as one value: each element percent-encoded (so an element may hold `,`), joined with `,`, `.` for an empty list, the element count checked against the declared `Number` (`r[record_encoder.count_check+3]`). `encode` writes one string and percent-encodes its commas, so a caller-joined list (`vcf_0,vcf_1`) would come out as `vcf_0%2Cvcf_1`, as aardvark's `SOURCES` did.

r[record_encoder.reserved_rejection+2]
Every entry point that accepts a caller-supplied number destined for a BCF numeric field MUST reject the reserved bit patterns of `r[bcf_encoder.reserved_bands]` before writing any bytes, and MUST report `VcfEncodeError::ReservedIntValue` / `VcfEncodeError::ReservedFloatValue` wrapped in `VcfError`. The affected entry points are the INFO methods `info_int`, `info_ints`, `info_int_opts`, `info_float`, `info_floats` and `info_float_opts`; the FORMAT methods `format_int`, `format_int_opt`, `format_ints`, `format_int_opts`, `format_float`, `format_float_opt`, `format_floats` and `format_float_opts` (the `_opt` forms for their concrete elements); and the `qual` argument of `begin_record`. `format_gt` needs no test: its BCF encoding is `(allele_index + 1) << 1 | phased` over a `u16` allele index, which is never negative, and the FILTER methods only ever write non-negative dictionary indices. Rejection MUST leave the encoder's buffers as they were, so a caller that handles the error can continue with the same record.

r[record_encoder.count_check+3]
Encoding an `Arr` or `OptArr` value (not a `Str`, which is one string whatever its `Number`) MUST check its length (missing elements included) against the field's declared `Number`: `A` expects the record's number of ALT alleles, `R` its number of alleles, `Count(n)` exactly `n`, and `G` the number of genotypes for the sample's own ploidy, `binomial(n_allele + ploidy - 1, ploidy)`, when the record's GT was already encoded (samples may differ in ploidy, `r[record_encoder.mixed_ploidy]`, so each sample is checked against its own GT; a sample whose GT has no alleles is not checked); before GT, and for INFO, `G` is not checked. `.`, `M`, `P` and the local numbers `LA`/`LR`/`LG` are not checked. FORMAT checks apply per sample, and a sample with no values is missing (`.`), which every `Number` allows. A mismatch MUST return `VcfError::CountMismatch { field, number, expected, found }` before any byte of the field is written, so the record stays valid and can still be emitted. htslib does not check this; a wrong count silently produces a file other tools misread (a reader indexing `Number=R` values by allele reads past the end, or the wrong allele).

r[record_encoder.foreign_keys]
Every encode entry point MUST compare the handle's `HeaderId` with the header the writer wrote and fail with `VcfError::ForeignKey { name }` before writing anything: the contig of `begin_record`, every `FilterId` of `filter_fail`, and the `FieldId` of every INFO and FORMAT method (GT included). `FilterId::PASS` is index 0 in every header and MUST be accepted by every writer. A key of another header would otherwise write a valid-looking record naming whatever field holds the same dictionary index. A later header-merging feature MAY let a header derived from another, with its indices preserved, accept its ancestor's keys; until then the comparison is exact.

r[record_encoder.info_state_queries]
`InfoEncoder` MUST provide `n_allele()` and `n_alt()` methods returning the number of alleles and alternate alleles for the current record.

r[record_encoder.info_dedup+2]
If the same INFO field (identified by `FieldId`) is encoded more than once within a single record, the encoder SHOULD overwrite the previously-written value. The final output MUST contain at most one entry per distinct INFO field. The replacement MUST keep the field's original column position (overwrite in place), matching htslib's `bcf_update_info_*` semantics. Duplicate **detection** MUST be O(1) per field (a `FieldSet` keyed by `dict_idx`, `r[vcf_header.field_set]`); the linear offset bookkeeping is only consulted when a duplicate is actually present, so the common write-each-field-once path does no scan. The deduplication tracker SHOULD reuse its allocations across records.

## FormatEncoder Trait

r[record_encoder.format_encoder+2]
A `FormatEncoder` trait MUST be provided for encoding FORMAT fields. It MUST be implemented by `RecordEncoder<'_, WithSamples>` and `RecordEncoder<'_, Edit>`.

r[record_encoder.format_methods+3]
FORMAT methods (`format_gt`, `format_int`, `format_int_opt`, `format_float`, `format_float_opt`, `format_ints`, `format_int_opts`, `format_floats`, `format_float_opts`, `format_string`) MUST accept a `&FieldId` and a **slice of per-sample values** — one element per sample for scalar/string/GT methods, or one inner slice per sample (`&[&[T]]`) for the multi-valued array methods. The outer slice length MUST equal the sample count the header declares: the typed `FormatKey` wrappers, the only way to call these methods from outside the crate (`FieldId` is crate-private), check it and fail with `VcfError::SampleCountMismatch`; the methods themselves only `debug_assert!` it. They return `Result<(), VcfError>` for what they do check: the key's header, the FORMAT field limit, value counts and reserved values. For single-sample records, callers pass a 1-element slice. For the multi-valued array methods, the BCF encoding stores a single column width equal to the longest inner slice; a sample with fewer values MUST be padded to the column width with the end-of-vector sentinel for the chosen type. A sample with **no** values (empty slice / empty string) MUST be encoded as the BCF missing value (so it renders as `.`, matching the VCF text path) — the per-type integer/float missing sentinel followed by end-of-vector padding, or the literal `.` char NUL-padded for strings. For multi-valued integers the BCF integer type MUST be chosen to fit every value across all samples (the column stores one type), not just the first sample.

r[record_encoder.format_scalar_opt]
`FormatKey<Scalar<i32>>` and `FormatKey<Scalar<f32>>` MUST provide `encode_opt(&self, enc, values: &[Option<T>])`, one element per sample (else `VcfError::SampleCountMismatch`), where `None` is that sample's missing value: `.` in VCF text, and in BCF the missing value of the width the column is stored in — `0x80`, `0x8000` or `0x80000000` for the integer type chosen to fit every concrete value (int8 when every sample is missing, as htslib chooses), `0x7F800001` for a float. Concrete values are checked for reserved patterns as `encode` checks them (`r[record_encoder.reserved_rejection+2]`); `None` never is.

r[record_encoder.mixed_ploidy]
Samples of one record MAY differ in ploidy (a haploid chrX call next to a diploid one, a triploid sample in a diploid cohort): `FormatKey<Gt>::encode` MUST accept them and MUST NOT fail for it. VCF text writes each sample's own genotype (`r[vcf_writer.genotype_serialization+2]`); BCF stores the GT column at the largest ploidy and pads each shorter sample with the end-of-vector value of the column's integer type (`r[bcf_writer.gt_encoding]`), which is how htslib stores a ploidy mix and how it reads it back. A genotype with no alleles is the missing call `.`, a missing haploid allele in BCF.

r[record_encoder.format_single_sample]
The multi-valued array FORMAT keys (`FormatKey<Arr<i32>>`, `FormatKey<Arr<f32>>`) MUST provide an `encode_single_sample` convenience that takes one sample's values as a flat `&[T]` (rather than `&[&[T]]`) and encodes them for a single-sample record. It MUST return [`VcfError::SampleCountMismatch`] if the header declares anything other than exactly one sample.

r[record_encoder.format_state_queries]
`FormatEncoder` MUST provide `n_allele()` and `n_alt()` methods returning the number of alleles and alternate alleles for the current record, and `n_samples()` returning the sample count the header declares — the same count the per-sample slice length is validated against in `r[record_encoder.format_methods+3]`.

r[record_encoder.format_dedup+2]
If the same FORMAT field (identified by `FieldId`) is encoded more than once within a single record, the encoder SHOULD overwrite the previously-written value. The final output MUST contain at most one entry per distinct FORMAT field. The replacement MUST keep the field's original column position (overwrite in place), matching htslib's `bcf_update_format_*` semantics — for VCF text this means rewriting the corresponding colon-delimited subfield of every sample. The deduplication tracker SHOULD reuse its allocations across records. Detecting a duplicate SHOULD cost O(1) on the common path where each field is written once — a `FieldSet` keyed by `dict_idx` (`r[vcf_header.field_set]`), as `FieldTracker` keeps for INFO — rather than a name comparison against every key already written; the text path's per-record FORMAT key list SHOULD clone the header's names rather than rebuild them from `&str`.

r[record_encoder.format_count_limit]
A record MUST NOT hold more than 255 distinct FORMAT fields: BCF stores the count in 8 bits (`n_fmt`), and htslib rejects more on read in VCF text as well (`MAX_N_FMT`). Since the dictionary holds far more than 255 ids (`r[vcf_header.dict_capacity+2]`), the header alone no longer bounds this. Encoding a 256th distinct FORMAT field into one record MUST return `VcfError::ValueOverflow` (field `n_fmt`, target `u8`), in VCF text and BCF alike, before writing any of its bytes, so the record stays valid and can still be emitted. Rewriting a field already present (`r[record_encoder.format_dedup+2]`) is always allowed.

## Emit

r[record_encoder.emit]
`emit()` MUST finalize the record and write it to the output. For BCF, this patches the fixed header and flushes to BGZF. For VCF text, this writes accumulated FORMAT fields and flushes the line buffer. `emit()` MUST consume the encoder by value, releasing the writer borrow. `emit()` is the only encoding method that performs I/O and MAY return an error.

r[record_encoder.rlen]
The record's reference length — BCF's `rlen` field and the end the co-produced index files the record under — MUST be computed at `emit()` by the same function the reader uses to compute a VCF text record's end (`r[vcf_values.rlen+3]`, htslib's `get_rlen`), from the record as encoded: the REF length, extended by INFO `END`, by the `SVLEN` of `<DEL>`/`<DUP>`/`<INV>`/`<CNV>` alleles and by FORMAT `LEN` of a `<*>`/`<NON_REF>` allele, each only where the header declares it `Type=Integer`. An INFO or FORMAT field encoded twice counts with the value the record keeps (`r[record_encoder.info_dedup+2]`). `Alleles::rlen()` alone is the REF length and is wrong for structural variants and gVCF blocks: their index entry would end at the anchor base, and a region query that overlaps only the rest of the variant would not find it. Only these three fields can change the length, so a record that encodes none of them takes `Alleles::rlen()` without running the function. An `rlen` past `i32::MAX` MUST fail with `VcfError::ValueOverflow`. The BCF `rlen` and the index entry MUST be the same value `bcftools` writes and indexes for the same record (`bcftools view -Ob`, `bcftools index`).

r[record_encoder.emit_no_samples]
`emit()` MUST be available on both `Filtered` (for records without samples) and `WithSamples` (for records with samples).

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
