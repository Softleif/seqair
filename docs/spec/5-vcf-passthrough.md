# VCF Passthrough Writing

> **Sources:** Copying the records a [reader](./2-vcf-1-reader.md) yields into a [`Writer`](./5-record-encoder.md), optionally editing fields on the way: VCF↔BCF conversion (rastair `convert`), and tools that rewrite a few FORMAT fields of every record and keep the rest (HiPhase: GT, PS, PF). htslib does this with `bcf_translate` and `bcf_update_*`; this spec does it through a map checked once and the record encoder's own rules. See [VCF Header](./5-vcf-header.md) for deriving the output header (`r[vcf_header.to_builder]`, `r[vcf_header.remove]`).

## Field map

r[vcf_passthrough.field_map]
`FieldMap::new(input: &VcfFile, output: &VcfHeader) -> Result<FieldMap, FieldMapError>` MUST build, once, the translation from the input file's ids to the output header's: every INFO, FORMAT and FILTER id of the input's dictionary (and, for VCF text, its reserved keys, `r[vcf_reader.field_table]`) to the output field of the same name in the same section, and every input contig to the output contig of the same name. It MUST remember both headers' ids, so a copy can check that the record and the writer are the map's (`r[vcf_passthrough.copy]`). Matching is by name: the output's dictionary indices need not be the input's (a derived header appends, removes and may renumber), and the copy writes the output's.

r[vcf_passthrough.field_map.compatible]
A field the input defines — in its header, or as a VCF 4.5 reserved key of a VCF text input — that the output declares in the same section MUST have the same `Type` and the same `Number` in both, else `FieldMapError::Incompatible { kind, name, input_type, input_number, output_type, output_number }` from `FieldMap::new`. Values are decoded by the input's definition and encoded by the output's: a changed `Type` cannot be encoded, and a changed `Number` would make the copied values break the definition the output declares. The check runs once, so no record fails for it.

r[vcf_passthrough.field_map.removed]
An INFO, FORMAT or FILTER definition the input header declares and the output header does not is removed: every copy MUST leave that field (or filter) out, as `bcftools annotate -x` does. Removing a definition from a derived header is how a caller says so (`r[vcf_header.remove]`), and that header is then the place to `register` a replacement under the same name with a new definition. A record whose failed filters are all removed MUST be written `PASS`, as `bcftools annotate -x FILTER/…` writes it: the filters it failed no longer apply, and `.` would say that none were applied.

r[vcf_passthrough.field_map.contigs]
Every contig the input header declares MUST be declared by the output header under the same name, else `FieldMapError::MissingContig { name }` from `FieldMap::new`: a record cannot be written without its contig, and failing at setup is better than failing at the first record on that contig after writing every one before it. A contig a VCF text input's records name without a declaration (`r[vcf_reader.contigs]`) MUST be looked up by name when a record uses it; if the output does not declare it either, the copy MUST fail with `VcfError::UndeclaredContig { name }`.

r[vcf_passthrough.field_map.samples]
The output header MUST have the input's samples, by name and in order, else `FieldMapError::SampleMismatch { index, input, output }`, naming the first column where they differ (`None` past the end of either list). A copy writes every sample's values in column order; subsetting or reordering samples is not offered.

r[vcf_passthrough.undeclared]
A VCF text input can hold keys and filters its header does not declare (`r[vcf_reader.undeclared+3]`, `r[vcf_reader.undeclared.filter]`). The output header decides what happens to them, by name: one it declares in the same section is copied (a VCF 4.5 reserved key checked as `r[vcf_passthrough.field_map.compatible]` says; any other undeclared key is read as a `String`, so the output's `Type` MUST be `String` or `Character`, else the copy fails with `FieldMapError::Incompatible`); one it does not declare cannot be written, since BCF needs a dictionary index and VCF text a definition, and the copy MUST fail with `VcfError::UndeclaredField { kind, name }`. `FieldMap::drop_undeclared(self) -> FieldMap` MUST instead leave such keys and filters out of every copy, logging one `warn!` per map naming the first one dropped. htslib keeps them by adding a `Number=1,Type=String` definition to its header when it meets them; a streaming writer has written its header before the first record, so it cannot. Failing is the default because dropping is data loss that nothing else would report; a caller who knows its input either declares the keys in the output header (the htslib outcome, chosen up front) or asks for them to be dropped.

## Copying a record

r[vcf_passthrough.copy]
`Writer<W, Ready>::copy_record(&mut self, map: &FieldMap, record: &FullRecord) -> Result<RecordEdit<'_>, VcfError>` MUST begin a record of the writer's output that holds everything `record` holds, translated through `map`. It MUST fail with `VcfError::ForeignFieldMap` unless the writer's header is `map`'s output, and with `VcfError::ForeignRecord` unless `record` was read with `map`'s input header, before writing anything. It writes, in the order of the record: CHROM (the mapped contig), POS, the IDs, REF and ALT, QUAL, FILTER (`PASS`, `.`, or the mapped failed filters in record order, `r[vcf_passthrough.field_map.removed]`), every INFO field in record order, and every FORMAT field in record order for every sample. A record at the telomere (POS 0) MUST fail with `VcfError::TelomereRecord { contig }`: the writer takes a 1-based position. A copy MUST be the same whatever the two formats are: VCF text, BGZF VCF or BCF in, any of them out.

r[vcf_passthrough.copy.encoder]
Values MUST be re-encoded by the record encoder from their `DynValue` (`r[vcf_values.dynamic+2]`), not copied as bytes, so the output follows the writer's rules whatever the input's format: the smallest integer type, missing and end-of-vector values, string padding, GT encoding, the record's `rlen` from END/SVLEN/LEN (`r[record_encoder.rlen]`), and the co-produced index. This is also what bcftools writes: `bcftools view` parses VCF text into binary and formats it again (so `0.10` comes out `0.1`), and a BCF that htslib wrote already holds the smallest types.

r[vcf_passthrough.copy.values]
A copied value MUST read back as the input's: a number vector missing as a whole (VCF `.`, BCF's lone missing value or empty vector) MUST be written as one missing element (`.`), and missing elements stay missing; a FORMAT value a sample does not hold (VCF text lets a sample drop trailing fields) MUST be written missing, as htslib fills it when it parses the line; a missing string MUST be written `.`, a missing genotype `.`, and a flag present. Copied values MUST NOT be checked against their `Number` (`r[record_encoder.count_check+3]`): they are the input's, which the reader accepted and bcftools copies as they are. Strings MUST be written as the reader returns them, without percent-encoding (`r[vcf_values.text_strings]` returns them undecoded, so encoding `%` or `,` again would change them), except for the bytes that would break a VCF text line's structure — TAB, LF, CR, and `;` in an INFO value or `:` in a FORMAT value — which MUST be percent-encoded.

r[vcf_passthrough.copy.equivalence]
Copying MUST be checked against bcftools, in every direction (VCF text and BCF in; VCF, BGZF VCF and BCF out): `bcftools view` of the copy MUST print the records `bcftools view` prints for the input, on generated records that cover every value type, missing values and elements, dropped trailing FORMAT fields, mixed ploidy, IDs, filters and symbolic alleles, and noodles MUST read every copied record.

## Editing a copied record

r[vcf_passthrough.edit]
`RecordEdit<'a>` MUST be the record encoder in an `Edit` state (`RecordEncoder<'a, Edit>`), which implements both `InfoEncoder` and `FormatEncoder`, so every typed key edits a copy the way it encodes a new record: `key.encode(&mut edit, value)`. A key MUST be the output header's (`r[record_encoder.foreign_keys]`) and its value is checked as any encode checks it (count, reserved values, sample count, the FORMAT field limit). Encoding a field the copy holds MUST replace its value where it stands (`r[record_encoder.info_dedup+2]`, `r[record_encoder.format_dedup+2]`, htslib's `bcf_update_*`); a field it does not hold is appended to its section. INFO and FORMAT edits may come in any order: the copy has written both sections, and they are kept apart until `emit`.

r[vcf_passthrough.edit.remove]
`RecordEdit::remove_info(&mut self, key: &InfoKey<V>) -> Result<bool, VcfError>` and `remove_format(&mut self, key: &FormatKey<V>)` MUST remove the field from the record, keeping the other fields in order, and return whether it was there; a key of another header fails with `ForeignKey`. Removing INFO END or SVLEN or FORMAT LEN MUST take it out of the record's `rlen` too, and removing GT MUST forget the ploidies `Number=G` counts are checked against.

r[vcf_passthrough.edit.emit]
`RecordEdit::emit(self) -> Result<(), VcfError>` MUST write the record as `RecordEncoder::emit` does (`r[record_encoder.emit]`). A `RecordEdit` dropped without `emit` writes nothing.

r[vcf_passthrough.edit.alleles]
A copy MUST NOT offer to change REF or ALT: every `Number=A`, `R` and `G` field and every GT would then name the wrong alleles. A caller who needs other alleles writes a new record with `begin_record`.

## Encoder additions

r[vcf_passthrough.opt_arrays]
Arrays with missing elements MUST be encodable in every section, since a copy holds them wherever the input did: `InfoKey<OptArr<f32>>::encode(enc, &[Option<f32>])` through `InfoEncoder::info_float_opts`, `FormatKey<OptArr<i32>>::encode(enc, &[&[Option<i32>]])` through `FormatEncoder::format_int_opts`, and `FormatKey<OptArr<f32>>::encode(enc, &[&[Option<f32>]])` for every sample (beside `encode_single_sample`). A `None` element is `.` in VCF text and the missing value in BCF — of the integer type chosen to fit every concrete value of the column, as `r[record_encoder.format_scalar_opt]` chooses it; samples shorter than the column are padded with end-of-vector values and a sample with no values is missing, as `r[record_encoder.format_methods+2]` says for arrays without missing elements. Counts and reserved values are checked as for the other arrays.
