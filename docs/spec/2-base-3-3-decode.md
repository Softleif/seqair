# Base-typed Sequence Decoding

> **Sources:** [SAM1] §4.2.4 "SEQ and QUAL encoding" — 4-bit nibble-to-base mapping (`=ACMGRSVTWYHKDBN`). The `Base` enum, SIMD acceleration, and direct-to-`Base` decode table are seqair-specific design choices; the upstream spec defines only the nibble codes. See [References](./99-references.md).

## Background

BAM sequences are stored as 4-bit nibbles (two bases per byte). The standard decode table maps nibbles to ASCII characters (`=ACMGRSVTWYHKDBN`). Downstream, the ASCII bytes must be converted to `Base` enum values (A, C, G, T, Unknown) at every pileup position — billions of times for whole-genome data.

By decoding directly into `Base` values at record construction time, we:

1. Eliminate the per-position `Base::from(u8)` conversion (a LUT lookup per read per position)
2. Get compile-time type safety — a `Base` is guaranteed to be one of 5 valid values, not an arbitrary byte
3. Collapse IUPAC ambiguity codes (M, R, W, S, Y, K, etc.) and `=` to `Base::Unknown` at decode time, instead of at each use site

## Decode table

r[base_decode.table]
A 16-entry lookup table MUST map 4-bit BAM nibble codes directly to `Base` discriminant values: 1→A(65), 2→C(67), 4→G(71), 8→T(84), and all other codes (0, 3, 5, 6, 7, 9–14, 15) → Unknown(78).

## Table invariant

r[base_decode.table_invariant]
Every entry in `DECODE_BASE_TYPED` (16 entries) and `DECODE_PAIR_TYPED` (256 × 2 = 512 bytes) MUST be a valid `Base` discriminant: A(65), C(67), G(71), T(84), or Unknown(78). This MUST be verified by an exhaustive test covering all entries.

## Decoding function

r[base_decode.decode]
`decode_bases(packed_seq, seq_len)` MUST decode a 4-bit packed BAM sequence into a `Vec<Base>`. The function SHOULD use the SIMD-accelerated decode path (SSSE3/NEON) with a Base-valued lookup table for identical throughput to the ASCII decoder.

## Storage in RecordStore

r[base_decode.slab]
The `RecordStore` MUST store decoded bases in a dedicated slab (`Vec<Base>`) separate from the data slab. The data slab continues to hold cigar, qual, and aux data. The `seq_at` method MUST return `Base` instead of `u8`.

## PileupAlignment

r[base_decode.alignment]
`PileupAlignment.base` MUST be of type `Base`, not `u8`. This eliminates the `Base::from(u8)` conversion in the metrics accumulation hot loop.

## ASCII-to-Base batch conversion

Reference sequences (FASTA) and SAM text sequences arrive as ASCII bytes, not 4-bit packed nibbles. These need the same one-time batch conversion to `Base` that BAM sequences get via `decode_bases`, but using an ASCII-aware path instead of the nibble decode table.

r[base_decode.ascii_batch]
`Base::from_ascii_vec(Vec<u8>) -> Vec<Base>` MUST convert a vector of ASCII bytes to `Base` values in-place (reusing the allocation). Maps A/a→A, C/c→C, G/g→G, T/t→T, all other bytes→Unknown. The function MUST NOT allocate a new vector — it operates on the input's buffer and transmutes the result. This is safe because `Base` is `repr(u8)` and the function only writes valid `Base` discriminant bytes (65, 67, 71, 78, 84).

r[base_decode.ascii_simd]
The ASCII batch converter MUST use SIMD acceleration (SSSE3 on x86_64, NEON on aarch64) with a scalar fallback, following the same dispatch pattern as `decode_bases`. This is the required path for all bulk u8→Base conversions outside BAM 4-bit decoding: FASTA reference sequences, SAM text sequences, and any other ASCII-encoded base data.

r[base_decode.ascii_scalar_equivalence]
The SIMD ASCII converter MUST produce identical output to applying `Base::from(u8)` element-wise. This MUST be verified by property-based tests covering arbitrary byte values and lengths including SIMD boundary lengths (0, 1, 15, 16, 31, 32, 33, 63, 64, 65, 128).

## FromStr validation

r[types.base.from_str_validation+2]
`FromStr for Base` MUST accept only a single character A, C, G, T, or N (case-insensitive) after trimming whitespace. Multi-character input (after trimming) MUST return `Err(BaseError::MultipleChars)`. Any single character that is not A/C/G/T/N MUST return `Err(BaseError::InvalidBase(byte))` where `byte` is the first byte of the trimmed input. Empty input MUST return `Err(BaseError::Empty)`.

## Error display

r[types.base.error_display]
`BaseError` variants MUST display the actual invalid value.
`InvalidBase(byte)` SHOULD format as `"Invalid base: 0x{byte:02x}"`.
`Empty` SHOULD format as `"Empty"`.
