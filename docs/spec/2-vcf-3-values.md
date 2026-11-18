# VCF/BCF Values: typed access to fields

Fields are read through the same typed keys the writer uses (`InfoKey<V>`, `FormatKey<V>`, with markers `Scalar<T>`, `Arr<T>`, `OptArr<T>`, `Flag`, `Str`, `Gt`), resolved once from the header and checked against its declarations. On-disk missing and end-of-vector values become `Option` at the boundary; no sentinel reaches the caller.

> **Sources:** [VCF45] §1.6.1 (fixed fields, END/SVLEN/LEN), §1.6.2 (genotype fields), §6.3 (BCF typed values, missing and end-of-vector values). [htslib] `vcf.c` `get_rlen`. See [References](./99-references.md).

## Keys

r[vcf_values.keys]
`VcfHeader` MUST resolve keys by name: `info_key::<V>(name) -> Result<InfoKey<V>, KeyError>`, `format_key::<V>(name)`, `filter(name) -> Result<FilterId, _>`, `contig(name) -> Result<ContigId, _>`. Resolution MUST check the marker against the header's declaration through a `FieldKind` trait on each marker: `UnknownField { name }` for an undeclared name, `TypeMismatch { name, declared, requested }` when the `ValueType` disagrees, `NumberMismatch { name, declared }` when the cardinality cannot hold (a `Scalar` needs `Number=1`, a `Flag` needs `Type=Flag`, `Gt` needs the FORMAT `GT` field).

r[vcf_values.field_kind]
`FieldKind` MUST be the single source of a marker's `ValueType` and accepted `Number`s, used by reader key resolution and by the writer's `InfoFieldDef`/`FormatFieldDef`, so a definition whose marker disagrees with its declared `ValueType` (`InfoFieldDef::<Scalar<i32>>` with `ValueType::Float`) can no longer be written.

r[vcf_values.header_brand]
Every key and `SampleIdx` MUST carry its header's `HeaderId` (`r[vcf_header.brand]`), and every access MUST compare it with the record's, in release builds too, failing with `ForeignKey`. Without it, a key from another file reads whatever field holds the same dictionary index. Moving records between headers goes through an explicit mapping, never a portable key.

## Access

> r[vcf_values.access]
> `get(&key)` MUST return `Result<V::Out<'_>, VcfReadError>`, with `Err` only for malformed data, and:
>
> | marker | `Out` | absent or `.` |
> |---|---|---|
> | `Scalar<i32>`, `Scalar<f32>` | `Option<T>` | `None` |
> | `Arr<T>` | `Option<Values<'_, T>>` | `None`; a missing element inside is `Err(MissingElement)` |
> | `OptArr<T>` | `Option<OptValues<'_, T>>` yielding `Option<T>` | `None` |
> | `Flag` | `bool` | `false` |
> | `Str` | `Option<&str>` | `None` |
> | `Gt` | `Option<GenotypeRef<'_>>` | `None` |
>
> `Values` and `OptValues` MUST be `ExactSizeIterator`s with `get(i)` over the vector as stored, whatever its integer width. A field absent from the record and a field present as `.` are both `None`: VCF lets a sample drop trailing FORMAT fields, so callers cannot rely on the difference.

r[vcf_values.bcf_sentinels]
BCF values MUST be decoded at their stored width (int8, int16, int32, float) and widened. The per-width missing value (`0x80`, `0x8000`, `0x8000_0000`; float bits `0x7F80_0001`) MUST become `None`, and the end-of-vector value (`0x81`, `0x8001`, `0x8000_0001`; float bits `0x7F80_0002`) MUST end the vector, never yielded. Floats MUST be compared by bit pattern: any other NaN is a value, `Some(NaN)`.

r[vcf_values.text_numbers]
VCF text integers MUST parse as an optional sign and decimal digits into `i32`; anything else, or a value out of range, is a typed error naming the field. Floats MUST parse with `f32::from_str` (correctly rounded) and accept the VCF forms `Inf`, `Infinity`, `NaN` in any case with an optional sign. `.` is missing.

r[vcf_values.text_strings]
String values MUST be returned as stored, after UTF-8 validation (`from_utf8`, typed error); percent-encoding (`%3A`, …) MUST NOT be decoded, matching htslib, so VCF and BCF reads agree. BCF strings MUST be trimmed at their first NUL (padding).

r[vcf_values.gt]
`GenotypeRef` MUST yield each allele as `Option<u16>` (`None` for `.`) with its phase, read the first allele's phase bit per `r[vcf_record.gt_first_phase]`, handle a ploidy that varies per sample (BCF pads with end-of-vector), and convert to the writer's owned `Genotype`.

r[vcf_values.qual]
`qual()` MUST return `Option<f32>`: `None` for `.` or BCF's missing bit pattern, the stored `f32` otherwise, negative, huge and NaN values included. It MUST NOT return `Phred`: that type is `f64`, documented non-negative, and its `as_int` clamps to 0..=99, while QUAL values reach the thousands. BCF to BCF stays bit-exact.

r[vcf_values.duplicates]
When a record holds the same INFO or FORMAT key twice, the first occurrence MUST win, as htslib's lookups do, and the reader MUST log one `warn!` per file. (The writer keeps the last write of a key, `r[record_encoder.info_dedup+2]`; that is about writing, not reading.)

r[vcf_values.dynamic]
`Info::iter()` and the sample views MUST offer untyped iteration over the declared fields present, yielding the field's header definition and a `DynValue` (`Int(Values<i32>)`, `Float(Values<f32>)`, `Flag`, `Str(&str)`, with missing elements as `Option`), for generic tools such as `convert` and `view`.

## Reference span

r[vcf_values.rlen]
A record's end MUST be `start + rlen − 1`. For BCF, `rlen` is the stored field. For VCF text, `rlen` MUST be computed as htslib's `get_rlen` computes it for the file's version — the REF length, extended by INFO END, by the SVLEN of symbolic structural-variant alleles, and by FORMAT LEN of `<*>` — because indexes were built from that value and a reader that disagrees misses records. A differential test MUST compare the computed end with the `rlen` bcftools writes when converting the same text to BCF, on generated records covering every input.

## Per-record index

r[vcf_values.index]
A section's index MUST be one list of present fields in record order (`DictId`, kind, byte span) and an epoch-stamped table with one `u32` per dictionary id (`epoch << 16 | slot`), so a lookup is one table read, an epoch compare and one list read; an absent key is an epoch mismatch. Advancing the epoch replaces clearing; on wrap the table is zeroed. There is one representation for VCF and BCF and for every record size: on real files a linear scan beats the table by at most ~5 ns per record below ~16 keys, inside the noise, and a second code path is not worth it.

r[vcf_values.text_resolution]
VCF text keys MUST be mapped to `DictId`s through a table built once per file from the header: open addressing at ≤ 50% load, each slot holding the key's first 8 bytes, its length and its id, keys longer than 8 bytes verified against one flat name buffer. The hash MUST be one 8-byte load (plus the last 8 bytes for longer keys). A perfect hash measured no faster than this table on real headers and is not used.

## Equivalence

r[vcf_values.equivalence]
For the same logical record, reading its VCF text and its BCF MUST give identical values through every accessor, floats compared by value after the writer's round trip. A property test MUST generate records with the writer, write them as VCF and as BCF, and compare everything read back with what was written and between the two formats (extending `r[record_encoder.vcf_bcf_equivalence]` to reading).
