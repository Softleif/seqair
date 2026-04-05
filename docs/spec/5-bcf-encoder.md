# BCF Direct Encoder

> **Sources:** [BCF2] — typed value encoding, record layout, field-major FORMAT encoding, GT encoding. See [5-bcf-writer.md](5-bcf-writer.md) for the underlying BCF format rules. This spec extends the BCF writer with a zero-allocation direct-encode API.

## Motivation

The `BcfWriter::write_record(&VcfRecord)` path constructs an intermediate `VcfRecord` with heap-allocated `InfoFields`, `SmolStr` keys, and `InfoValue` enum boxing. For high-throughput callers (millions of records per chromosome), the direct encoder bypasses this intermediate entirely: pre-resolved typed handles write values straight into the BCF buffers.

## Handles

r[bcf_encoder.handles]
Field handles MUST be pre-resolved from the header at setup time via `resolve_*` methods on `BcfWriter`. Each handle stores the BCF dictionary index as an opaque `u32`. Resolution MUST validate that the field ID exists in the header and that the requested Number/Type matches the header declaration.

r[bcf_encoder.handle_types]
Handle types MUST encode the VCF Number semantics at the type level:
- `ScalarInfoHandle<T>` / `ScalarFormatHandle<T>` — `Number::Count(1)`, encodes exactly 1 value
- `FlagInfoHandle` — `Number::Count(0)`, encodes no value (presence = set)
- `PerAltInfoHandle<T>` / `PerAltFormatHandle<T>` — `Number::A`, encodes `n_alt` values
- `PerAlleleInfoHandle<T>` / `PerAlleleFormatHandle<T>` — `Number::R`, encodes `n_allele` values
- `GtFormatHandle` — genotype encoding with `(allele+1)<<1|phased` scheme
- `ContigHandle` — pre-resolved chromosome tid
- `FilterHandle` — pre-resolved filter dictionary index; PASS is always index 0

r[bcf_encoder.handle_value_type]
INFO and FORMAT handles MUST be generic over a value type `T: BcfValue`. The `BcfValue` trait defines how a Rust type maps to BCF typed encoding (type code, byte encoding, missing sentinel). Standard impls: `i32` (auto-selects smallest int type), `f32` (IEEE 754 LE). Domain types (`RootMeanSquare`, `Phred`) implement `BcfValue` by encoding as `f32`.

r[bcf_encoder.handle_encode]
Each handle type MUST provide an `encode(&self, enc: &mut BcfRecordEncoder, ...)` method. The method signature varies by handle kind:
- Scalar: `encode(&self, enc, value: T)`
- Flag: `encode(&self, enc)` (no value argument)
- PerAlt: `encode(&self, enc, values: &[T])` — debug-asserts `values.len() == enc.n_alt()`
- PerAllele: `encode(&self, enc, values: &[T])` — debug-asserts `values.len() == enc.n_allele()`
- GT: `encode(&self, enc, gt: &Genotype)`
- Filter: `encode(&self, enc)`
- Contig: used as argument to `Alleles::begin_record()`

## BcfValue trait

r[bcf_encoder.bcf_value]
The `BcfValue` trait MUST define: `bcf_type_code() -> u8` (BCF type code for this Rust type), `encode_bcf(&self, buf: &mut Vec<u8>)` (write value bytes), `encode_missing(buf: &mut Vec<u8>)` (write missing sentinel), `encode_end_of_vector(buf: &mut Vec<u8>)` (write EOV sentinel).

r[bcf_encoder.bcf_value_int]
Scalar `i32` values MUST select the smallest BCF integer type that fits, matching `r[bcf_writer.smallest_int_type]`. For arrays, the type is determined by scanning all concrete (non-missing) values first. Missing values within integer arrays MUST use the per-type sentinel (int8=0x80, int16=0x8000, int32=0x80000000) matching the selected type, not a fixed i32::MIN.

r[bcf_encoder.bcf_value_float]
`f32` values MUST be encoded as IEEE 754 single-precision LE bytes. Missing sentinel `0x7F800001` and EOV sentinel `0x7F800002` MUST be written as raw bytes (never through float arithmetic) per `r[bcf_writer.missing_sentinels]`.

## BcfRecordEncoder

r[bcf_encoder.encoder]
`BcfRecordEncoder` MUST borrow the `BcfWriter`'s `shared_buf` and `indiv_buf` buffers. It MUST NOT allocate during the encode-emit cycle.

r[bcf_encoder.begin_record]
`Alleles::begin_record(&self, enc, contig, pos, qual)` MUST write the 24-byte fixed header (with placeholder n_info/n_fmt), ID as `.`, and REF/ALT allele strings using zero-alloc `write_ref_into`/`write_alts_into`. It MUST set `enc.n_allele` and `enc.n_alt` for downstream validation.

r[bcf_encoder.emit]
`emit()` MUST patch the n_info, n_fmt, n_sample fields in the 24-byte fixed header, then flush the record to BGZF (with `flush_if_needed`), then push to the index builder if present.

r[bcf_encoder.info_counting]
Each INFO handle `encode()` call MUST increment the encoder's `n_info` counter. Each FORMAT handle `encode()` call MUST increment `n_fmt`.

r[bcf_encoder.format_field_major]
FORMAT fields MUST be encoded in field-major order per `r[bcf_writer.indiv_field_major]`. For single-sample records (the common case), each format encode call writes the key + type descriptor + 1 value.
