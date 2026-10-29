# Base-typed Sequence Decoding

> **Sources:** [SAM1] §4.2.4 "SEQ and QUAL encoding" — 4-bit nibble-to-base mapping (`=ACMGRSVTWYHKDBN`). The `Base` enum and where the conversion to it happens are seqair-specific design choices; the upstream spec defines only the nibble codes. See [References](./99-references.md).

## Background

Records keep SEQ as `IupacBase`, every code the file stores (`r[record_store.seq_lossless]`), so a record read from the store writes back unchanged. Calling wants the five-state `Base` (A, C, G, T, Unknown), so the pileup converts each base when it builds a column entry (`r[pileup.base_from_iupac]`). Both types are one byte, and the conversion is one load from a 16-entry table.

## Storage in RecordStore

r[base_decode.slab+1]
The `RecordStore` MUST store a record's bases in a dedicated slab (`Vec<IupacBase>`) separate from the cigar, qual, and aux slabs. `seq` and `base_at` MUST return the stored codes.

## PileupAlignment

r[base_decode.alignment+1]
`PileupAlignment.base` MUST be of type `Base`, not `u8` or `IupacBase`, made by `IupacBase::to_base` when the column entry is built (`r[pileup.base_from_iupac]`). Consumers of a column never convert.

## ASCII-to-Base batch conversion

Reference sequences (FASTA) and SAM text sequences arrive as ASCII bytes, not 4-bit packed nibbles. A reference becomes `Base` once per load; SAM text becomes `IupacBase` (`r[sam.record.seq_decode+1]`).

r[base_decode.ascii_batch+2]
`Base::from_ascii_vec(Vec<u8>) -> Vec<Base>` MUST convert a vector of ASCII bytes to `Base` values in-place (reusing the allocation). Maps A/a→A, C/c→C, G/g→G, T/t→T, U/u→T (RNA, as htslib's SAM parser does since 1.22), all other bytes→Unknown. `Base::from(u8)` uses the same mapping. The function MUST NOT allocate a new vector — it operates on the input's buffer and transmutes the result. This is safe because `Base` is `repr(u8)` and the function only writes valid `Base` discriminant bytes (65, 67, 71, 78, 84).

r[base_decode.ascii_simd+2]
The ASCII batch converter MUST be a native-width SIMD kernel per `r[io.simd_portable]`, following the same dispatch pattern as the BAM nibble decoder (`r[seq.decode_simd+2]`). This is the required path for all bulk u8→Base conversions: FASTA reference sequences and any other ASCII-encoded base data that is read as `Base`.

r[base_decode.ascii_scalar_equivalence]
The SIMD ASCII converter MUST produce identical output to applying `Base::from(u8)` element-wise. This MUST be verified by property-based tests covering arbitrary byte values and lengths including SIMD boundary lengths (0, 1, 15, 16, 31, 32, 33, 63, 64, 65, 128).

## FromStr validation

r[types.base.from_str_validation+2]
`FromStr for Base` MUST accept only a single character A, C, G, T, or N (case-insensitive) after trimming whitespace. Input that is more than one *byte* after trimming MUST return `Err(BaseError::MultipleChars)` — that includes a single non-ASCII character, whose leading byte would name nothing useful in an `InvalidBase`. Any single byte that is not A/C/G/T/N MUST return `Err(BaseError::InvalidBase(byte))` where `byte` is that byte. Empty input MUST return `Err(BaseError::Empty)`.

## Error display

r[types.base.error_display]
`BaseError` variants MUST display the actual invalid value.
`InvalidBase(byte)` SHOULD format as `"Invalid base: 0x{byte:02x}"`.
`Empty` SHOULD format as `"Empty"`.
