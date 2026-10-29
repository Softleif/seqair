# IUPAC Base Type

> **Sources:** SAM specification [SAM1] §4.2.3 (BAM `seq`: 4-bit codes for `=ACMGRSVTWYHKDBN`); htslib `seq_nt16_str` and `seq_nt16_table`. See [References](./99-references.md).

## Motivation

`Base` holds A, C, G, T and one `Unknown`, which is right for calling and wrong for a record that passes through unchanged: BAM can store sixteen codes, and collapsing `R` or `=` to `N` rewrites the file. `IupacBase` is the lossless base. It is the type of a record's sequence; the pileup and calling paths keep using `Base`.

## Type definition

r[types.iupac_base.type]
`IupacBase` MUST be a `#[repr(u8)]` enum with one variant per BAM 4-bit code, each variant's discriminant being that code: `Eq` (`=`) 0, `A` 1, `C` 2, `M` 3, `G` 4, `R` 5, `S` 6, `V` 7, `T` 8, `W` 9, `Y` 10, `H` 11, `K` 12, `D` 13, `B` 14, `N` 15. `size_of` MUST be 1.

## Codes

r[types.iupac_base.code]
`code()` MUST return the discriminant. `from_code(u8)` MUST return `Some` for 0..=15 and `None` for any larger value, and `from_code(b.code())` MUST be `Some(b)`. `from_nibble(u8)` MUST use the low four bits only, so it is total.

## ASCII

r[types.iupac_base.from_ascii]
`from_ascii(u8)` MUST map a byte as htslib's `seq_nt16_table` does: the sixteen letters of `=ACMGRSVTWYHKDBN` in either case, `U` and `u` as `T`, the colour-space digits `0`..`3` as `A`, `C`, `G`, `T`, and every other byte as `N`.

r[types.iupac_base.to_ascii]
`to_ascii()` MUST return the upper-case letter of `=ACMGRSVTWYHKDBN` at the code's index. For every byte `x`, `from_ascii(x).to_ascii()` is a letter that `from_ascii` maps back to the same variant.

r[types.iupac_base.slices]
`from_ascii_slice` and `to_ascii_slice` MUST equal the element-wise `from_ascii` and `to_ascii`. They are the entry points for a vectorized implementation, which MUST agree with the scalar one on every input.

## Conversions

r[types.iupac_base.complement]
`complement()` MUST be the IUPAC complement: A and T, C and G, M and K, R and Y, V and B, H and D swap; W, S, N and `=` map to themselves. It MUST be an involution.

r[types.iupac_base.base_conversions]
`From<Base>` MUST map A, C, G, T to the same letter and `Base::Unknown` to `N`. `to_base()` MUST be the lossy inverse: A, C, G, T to `Base`, every other variant to `Base::Unknown`. `to_base(from(b)) == b` for every `Base`.
