# BCF Encoding Primitives

> **Sources:** [BCF2] — typed value encoding, record layout, field-major FORMAT encoding, GT encoding. See [BCF Writer](./5-bcf-writer.md) for the underlying BCF format rules. These primitives are used internally by the unified [`RecordEncoder`](./5-record-encoder.md) for BCF output. Also see [References](./99-references.md).

## BcfValue trait

r[bcf_encoder.bcf_value]
The `BcfValue` trait MUST define: `bcf_type_code() -> u8` (BCF type code for this Rust type), `encode_bcf(&self, buf: &mut Vec<u8>)` (write value bytes), `encode_missing(buf: &mut Vec<u8>)` (write missing sentinel), `encode_end_of_vector(buf: &mut Vec<u8>)` (write EOV sentinel).

r[bcf_encoder.bcf_value_int]
Scalar `i32` values MUST select the smallest BCF integer type that fits, matching `r[bcf_writer.smallest_int_type]`. For arrays, the type is determined by scanning all concrete (non-missing) values first. Missing values within integer arrays MUST use the per-type sentinel (int8=0x80, int16=0x8000, int32=0x80000000) matching the selected type, not a fixed i32::MIN.

r[bcf_encoder.bcf_value_float]
`f32` values MUST be encoded as IEEE 754 single-precision LE bytes. Missing sentinel `0x7F800001` and EOV sentinel `0x7F800002` MUST be written as raw bytes (never through float arithmetic) per `r[bcf_writer.missing_sentinels]`.

## BCF Record Encoding

r[bcf_encoder.begin_record]
Beginning a BCF record MUST write the 24-byte fixed header (with placeholder n_info/n_fmt), ID as `.`, and REF/ALT allele strings using zero-alloc `write_ref_into`/`write_alts_into`. It MUST set `n_allele` and `n_alt` for downstream validation.

r[bcf_encoder.checked_casts]
All conversions from `usize`/`u32` to narrower integer types (`i32`, `u16`, `u8`) for user-controlled values (contig tid, position, allele count, sample count) MUST use checked conversion (`try_from`) and return a typed error on overflow rather than silently truncating via `as`.

r[bcf_encoder.emit]
BCF `emit()` MUST patch the n_info, n_fmt, n_sample fields in the 24-byte fixed header, then flush the record to BGZF (with `flush_if_needed`), then push to the index builder if present.

r[bcf_encoder.info_counting]
Each INFO field encoded MUST increment the encoder's `n_info` counter. Each FORMAT field encoded MUST increment `n_fmt`.

r[bcf_encoder.format_field_major]
FORMAT fields MUST be encoded in field-major order per `r[bcf_writer.indiv_field_major]`. For single-sample records (the common case), each format encode call writes the key + type descriptor + 1 value.

## Reserved values

> _[BCF2] — "In total, eight values are reserved for future use: 0x80--0x87, 0x8000--0x8007, 0x80000000--0x80000007" for integers, and 0x7F800001--0x7F800007 for floats. See also `r[bcf_writer.smallest_int_type]`, `r[bcf_writer.missing_sentinels]` and `r[bcf_writer.end_of_vector]`._

r[bcf_encoder.reserved_bands]
Each BCF numeric width has a band of bit patterns that carry meaning to the reader instead of data. For integers the band is the eight most negative values of the width: `MIN` is MISSING, `MIN + 1` is END_OF_VECTOR, and `MIN + 2 ..= MIN + 7` are reserved for future use — `0x80..=0x87` for int8, `0x8000..=0x8007` for int16, `0x80000000..=0x80000007` for int32. For floats the band is the bit patterns `0x7F800001 ..= 0x7F800007` (MISSING, END_OF_VECTOR, then five reserved). Positive infinity (`0x7F800000`) and quiet NaN (`0x7FC00000`) are NOT in the band: [BCF2] gives both first-class status as ordinary Float values, so a reserved-value test MUST compare exact bit patterns and MUST NOT use `is_nan()`.

r[bcf_encoder.reserved_rejected]
A caller-supplied value whose encoding would land in the reserved band MUST be rejected with a typed error instead of written. Writing END_OF_VECTOR as data makes a conforming reader truncate the row at that element, and writing MISSING makes it drop the value — both are silent data loss that no later read-back can detect. The error MUST name the field, the offending value (its bit pattern, for a float), and which band member it collided with.

r[bcf_encoder.reserved_width]
A value is reserved only for the width it is actually emitted at, so the test MUST use the width `r[bcf_writer.smallest_int_type]` selects, not the width of the caller's Rust type: `-127` is END_OF_VECTOR as an int8 but ordinary data as an int16, and `0x81` as an `i32` is the perfectly ordinary value 129. Because the selection bounds `-120` and `-32760` are exactly `INT8_MIN + 8` and `INT16_MIN + 8`, no value can ever be emitted at int8 or int16 width with a reserved pattern — `-128` is promoted to int16 and written `0xFF80`. The test therefore reduces, for every integer entry point and for arrays as well as scalars, to the width-independent `value <= i32::MIN + 7`; that reduction MUST be verified against the bytes the encoder actually emits rather than assumed.

r[bcf_encoder.reserved_uniform]
The rejection MUST apply to VCF text output as well as BCF, even though VCF text could render the value losslessly. The unified writer's contract is that the same encoding calls produce the same accepted records in every output format; a value that only survives as text would be lost the moment the file is converted to BCF, which is what htslib does internally for all output.
