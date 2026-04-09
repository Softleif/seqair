# Seqair-Types — Phred Quality Score

> **Sources:** Phred quality scores are defined in [SAM1] section 1.4 (QUAL field, Phred+33 encoding). The `Phred` newtype wraps the mathematical definition: `Q = -10 log10(P)`. See [References](./99-references.md).

## Construction

r[types.phred.non_negative]
Phred scores MUST be non-negative. `from_phred` MUST accept only `u8` values (0..=255), making negative inputs unrepresentable at the type level.

## Integer conversion

r[types.phred.as_int_clamp]
`as_int` MUST clamp the result to the range [0, 99]. Values less than or equal to zero (including negative results from floating-point arithmetic) MUST return 0. NaN MUST return 0.
