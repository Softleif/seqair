# VCF/BCF Values: typed access to fields

Fields are read through the same typed keys the writer uses (`InfoKey<V>`, `FormatKey<V>`, with markers `Scalar<T>`, `Arr<T>`, `OptArr<T>`, `Flag`, `Str`, `Gt`), resolved once from the header and checked against its declarations. On-disk missing and end-of-vector values become `Option` at the boundary; no sentinel reaches the caller.

> **Sources:** [VCF45] §1.6.1 (fixed fields, END/SVLEN/LEN), §1.6.2 (genotype fields), §6.3 (BCF typed values, missing and end-of-vector values). [htslib] `vcf.c` `get_rlen`. See [References](./99-references.md).

## Keys

r[vcf_values.keys+2]
`VcfHeader` MUST resolve keys by name: `info_key::<V>(name) -> Result<InfoKey<V>, KeyError>`, `format_key::<V>(name)`, `filter(name) -> Result<FilterId, _>`, `contig(name) -> Result<ContigId, _>`. A header resolves only what it declares: it is what a writer writes, and the keys it hands out are valid for writing. Resolution MUST check the marker against the declaration through a `FieldKind` trait on each marker: `UnknownField { name, kind }` for an undeclared name, naming the section (`KeyKind`), `TypeMismatch { name, declared, requested }` when the `ValueType` disagrees, `NumberMismatch { name, declared }` when the cardinality cannot hold (a `Scalar` needs `Number=1`, a `Flag` needs `Type=Flag`, `Gt` needs the FORMAT `GT` field).

r[vcf_values.file_keys+2]
Readers resolve keys through their file: `file.info_key::<V>(name)` and `file.format_key::<V>(name)` (and the same on `PlanContext`) MUST resolve exactly what the header declares, as `r[vcf_values.keys+2]` does, for VCF text and BCF alike. The same file in both formats then resolves the same keys: HiPhase decided whether a file carries tandem-repeat annotations by whether `TRID` resolved, and with an undeclared `TRID` resolving in text but not in BCF it classified every text record as a tandem repeat. A key a header leaves out needs a lenient key (`r[vcf_values.lenient_keys]`); a record that uses one is still reported (`r[vcf_reader.undeclared+4]`).

r[vcf_values.lenient_keys]
`file.lenient()` (and `PlanContext::lenient()`) MUST return a `LenientKeys` view with the same four methods, resolving against the file's field table (`r[vcf_reader.field_table]`):

| name | result |
|---|---|
| declared in the header | the header's key (`r[vcf_values.keys+2]`) |
| reserved by VCF 4.5, undeclared | a key checked against the reserved `Type`/`Number` |
| other undeclared, `V` is `Str` or `Any` | a key read as `Type=String, Number=.` |
| other undeclared, any other `V` | `UnknownField { name, kind }` |

Resolving an undeclared name logs a `debug!` naming the `Type` and `Number` it is read with; a record that uses one is warned of (`r[vcf_reader.undeclared+4]`). The result MUST NOT depend on the format: in BCF the name gets an id past the dictionary, which no record can name, so the key resolves and finds nothing. A key the header does not declare MUST carry the file's own id rather than the header's, so the file's readers accept it and every writer refuses it as foreign (`r[vcf_values.header_brand+2]`): a writer would write a dictionary index its header lacks.

r[vcf_values.optional_keys+2]
`info_key_opt`/`format_key_opt`, on the file, on `PlanContext` and on `LenientKeys`, MUST return `Result<Option<Key<V>>, KeyError>`: `Ok(None)` exactly where the plain form fails with `UnknownField`, every other error kept. An optional annotation (GIAB has no CPG, many headers no SVTYPE) is then one `?` and one `Option`. There MUST NOT be a `get(Option<&Key>)`: the `Option` belongs to the plan, not to every record.

r[vcf_values.field_kind]
`FieldKind` MUST be the single source of a marker's `ValueType` and accepted `Number`s, used by reader key resolution and by the writer's `InfoFieldDef`/`FormatFieldDef`, so a definition whose marker disagrees with its declared `ValueType` (`InfoFieldDef::<Scalar<i32>>` with `ValueType::Float`) can no longer be written.

r[vcf_values.header_brand+2]
Every key and `SampleIdx` MUST carry its header's `HeaderId` (`r[vcf_header.brand]`), and every access MUST compare it with the record's, in release builds too, failing with `ForeignKey`. A key of a file's field table that its header does not declare carries a second id, the file's own, which that file's records accept as well. Without it, a key from another file reads whatever field holds the same dictionary index. Moving records between headers goes through an explicit mapping, never a portable key.

r[vcf_values.reader_keys]
`Reader::keys()` MUST return a `Keys` view resolving as `file.info_key`, `format_key`, their `_opt` forms and `lenient()` do, for the reader's file: a `VcfStream` reader has no `VcfFile` at hand, and without it a caller re-implements `_opt` by matching `UnknownField`.

## Access

> r[vcf_values.access+2]
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
> | `Any` | `Option<DynValue<'_>>` | absent: `None`; `.`: `Some(v)` with `v.is_missing()` |
>
> `Values` and `OptValues` MUST be `ExactSizeIterator`s with `get(i)` over the vector as stored, whatever its integer width; `Values::to_array::<N>()` MUST give `Some([T; N])` exactly when it holds `N` elements. `OptArr<T>` reads a field of any `Number` (it is the lenient key: `Number=.` GQ, a scalar declared as a vector); `Arr<T>` fails with `MissingElement` on an inner `.`. `Any` resolves for any declared type and reads the field as its declaration says (`r[vcf_values.dynamic+2]`).
>
> `get` folds a field absent from the record and a field present as `.` into the marker's absent value; `has(&key) -> Result<bool, VcfReadError>`, on `Info` and every sample view and for every marker, MUST tell them apart: for INFO, whether the record holds the key; for FORMAT, whether the record's FORMAT column names it, whatever the sample's value (a sample may drop trailing fields, which reads as `.`). So absent is `!has`, and `.` is `has` with an absent `get`. `has` on one sample's view MUST need only that sample's section: HiPhase phases a variant whose GQ is absent but not one whose GQ is `.`, and asking for every sample to learn it costs a column decode.

r[vcf_values.bcf_sentinels]
BCF values MUST be decoded at their stored width (int8, int16, int32, float) and widened. The per-width missing value (`0x80`, `0x8000`, `0x8000_0000`; float bits `0x7F80_0001`) MUST become `None`, and the end-of-vector value (`0x81`, `0x8001`, `0x8000_0001`; float bits `0x7F80_0002`) MUST end the vector, never yielded. Floats MUST be compared by bit pattern: any other NaN is a value, `Some(NaN)`.

r[vcf_values.text_numbers]
VCF text integers MUST parse as an optional sign and decimal digits into `i32`; anything else, or a value out of range, is a typed error naming the field. Floats MUST parse with `f32::from_str` (correctly rounded) and accept the VCF forms `Inf`, `Infinity`, `NaN` in any case with an optional sign. `.` is missing.

r[vcf_values.text_strings]
String values MUST be returned as stored, after UTF-8 validation (`from_utf8`, typed error); percent-encoding (`%3A`, …) MUST NOT be decoded, matching htslib, so VCF and BCF reads agree. BCF strings MUST be trimmed at their first NUL (padding).

r[vcf_values.gt+2]
`GenotypeRef` MUST yield each allele as `Option<u16>` (`None` for `.`) with its phase, offer `get(i)` and `ploidy()`, read the first allele's phase bit per `r[vcf_record.gt_first_phase]`, handle a ploidy that varies per sample (BCF pads with end-of-vector), and convert to the writer's owned `Genotype`.

r[vcf_values.gt_phased]
`GenotypeRef::is_phased()` MUST say whether the genotype is phased whatever the file's version, by htslib's inference for a first allele without a phase indicator: a haploid genotype is phased unless its allele is `.`, any other is phased iff every allele after the first is. Before VCF 4.4 the first allele's own phase bit is always unphased (`r[vcf_record.gt_first_phase]`), so code ported from htslib that compares the first allele's phase with the others' reads `0|1` as unphased in VCF 4.3 (rastair's `convert` did); `is_phased` is what such code means.

r[vcf_values.qual]
`qual()` MUST return `Option<f32>`: `None` for `.` or BCF's missing bit pattern, the stored `f32` otherwise, negative, huge and NaN values included. It MUST NOT return `Phred`: that type is `f64`, documented non-negative, and its `as_int` clamps to 0..=99, while QUAL values reach the thousands. BCF to BCF stays bit-exact.

r[vcf_values.duplicates+2]
When a record holds the same INFO or FORMAT key twice, the first occurrence MUST win, as htslib's lookups do, and the reader MUST log one `warn!` per file when an index meets a repeated key (a lookup that stops early, `r[vcf_values.index+2]`, may not reach it). (The writer keeps the last write of a key, `r[record_encoder.info_dedup+2]`; that is about writing, not reading.)

r[vcf_values.dynamic+2]
`Info::iter()` and the sample views MUST offer untyped iteration over the fields present, yielding each field's name, its definition in the field table (the header's, a reserved one, or the undeclared `String` one) and a `DynValue` (`Int(OptValues<i32>)`, `Float(OptValues<f32>)`, `Flag`, `Str(Option<&str>)`, `Gt(Option<GenotypeRef>)`), for generic tools such as `convert` and `view`. The iterator MUST be one layer of `Result`: a section that does not index is its first item's error. `DynValue::is_missing()` MUST say whether the value is `.`, whatever its variant.

r[vcf_values.dynamic.missing]
A number vector missing as a whole — `.` in VCF text, an empty vector or a lone missing value in BCF — MUST read as a `DynValue` with no elements. VCF text cannot spell an empty vector, so this is the only reading on which the formats agree (`r[vcf_values.equivalence]`).

## Reference span

r[vcf_values.rlen+2]
A record's end MUST be `start + rlen − 1`. For BCF, `rlen` is the stored field. For VCF text, `rlen` MUST be computed as htslib's `get_rlen` computes it for the file's version — the REF length, extended by INFO END, by the SVLEN of symbolic structural-variant alleles, and by FORMAT LEN of `<*>` — because indexes were built from that value and a reader that disagrees misses records. END, SVLEN and LEN count when the field table (`r[vcf_reader.field_table]`) gives them `Type=Integer`: declared so, or not declared at all, as tabix reads them (htslib's text parser reads an undeclared one as a String and ignores it). A differential test MUST compare the computed end with the `rlen` bcftools writes when converting the same text to BCF, on generated records covering every input.

r[vcf_values.rlen.end_before_pos]
An INFO END before POS (1-based `END < POS`; `END = POS` is a one-base span) MUST NOT shorten the span — `get_rlen` ignores it, so the REF length stands — and MUST log a single `warn!` per file, with the file's path and the first such record's `CHROM:POS` and END, when the reader first computes such a record's span (a VCF text record's span is computed lazily: for its `Locus`, or a query that needs its end). bcftools warns of the same record (`INFO/END=… is smaller than POS`); seqair used to take it silently. Of a duplicated END the first is the record's (`r[vcf_values.duplicates+2]`): it alone sets the span and is checked, as tabix's `tbx_parse1` — which indexes VCF text — reads it. (htslib's `get_rlen` is not consistent here: when the header has an SVLEN id it keeps overwriting END until it meets SVLEN, so bcftools converting `END=50;END=30` to BCF stores the span to 30, and without one it stops at the first.) A BCF record's `rlen` is stored, so nothing is computed or warned of.

## Per-record index

r[vcf_values.index+2]
A section's index MUST be one list of present fields in record order (`DictId`, kind, byte span) and an epoch-stamped table with one `u32` per dictionary id (`epoch << 16 | slot`), so a lookup is one table read, an epoch compare and one list read; an absent key is an epoch mismatch. Advancing the epoch replaces clearing; on wrap the table is zeroed. There is one representation for VCF and BCF and for every record size: on real files a linear scan beats the table by at most ~5 ns per record below ~16 keys, inside the noise, and a second code path is not worth it.

The INFO index MUST be built incrementally: a lookup of one key (`get`, `has`) indexes the section no further than that key's first entry, resuming where the last lookup stopped; iteration, and a lookup of a key the section lacks, index the rest. On wide INFO columns (gnomAD's hundreds of keys) a lookup of an early key then skips most of the work. Malformed data fails only the accesses that reach it: a lookup of a key whose first entry precedes it reads that entry, whatever the order of accesses, and every other access returns the error. A property test MUST check lookups and whole-index accesses, in any order, against the index built whole.

r[vcf_values.text_resolution+2]
VCF text keys MUST be mapped to `DictId`s through a table built once per file from the header and the reserved keys it does not declare: open addressing at ≤ 50% load, each slot holding the key's first 8 bytes, its length and its id, keys longer than 8 bytes verified against one flat name buffer. The hash MUST be one 8-byte load (plus the last 8 bytes for longer keys). A perfect hash measured no faster than this table on real headers and is not used. Keys appended to the field table later are found in a second, locked map consulted only when the first misses, so a file whose keys are all declared never takes the lock.

## Equivalence

r[vcf_values.equivalence]
For the same logical record, reading its VCF text and its BCF MUST give identical values through every accessor, floats compared by value after the writer's round trip. A property test MUST generate records with the writer, write them as VCF and as BCF, and compare everything read back with what was written and between the two formats (extending `r[record_encoder.vcf_bcf_equivalence]` to reading).
