# Sequence Codec

BAM encodes DNA sequences in a compact 4-bit format: two bases per byte, high nibble first. A 150bp read takes 75 bytes instead of 150. The 16 possible nibble values map to the 4 standard bases (A=1, C=2, G=4, T=8), N (15), and IUPAC ambiguity codes for the remaining values.

Decoding these sequences is on the hot path: during pileup construction, the engine needs the base at each read's query position for every reference position. With ~30 reads at ~3 billion positions, that's ~90 billion base lookups for a whole-genome BAM. Encoding is needed when rewriting BAM files with modified sequences.

> **Sources:** [SAM1] §4.2.4 "SEQ and QUAL encoding" — 4-bit encoding table `=ACMGRSVTWYHKDBN → [0,15]`, high-nibble-first packing. The SIMD decode path and pair-table optimisation are seqair-specific. See [References](./99-references.md).

## Decoding (4-bit → ASCII)

r[seq.decode_scalar]
The scalar decoder MUST use a 16-entry lookup table mapping 4-bit codes to ASCII bases (`=ACMGRSVTWYHKDBN`). It MUST support decoding full sequences (all bases) and single bases at arbitrary positions.

r[seq.decode_pair_table]
For bulk decoding, a 256-entry pair table (`DECODE_PAIR[byte] → [high_base, low_base]`) SHOULD be used to decode two bases per byte lookup, avoiding repeated shift-and-mask operations.

r[seq.decode_simd+2]
Bulk sequence decoding MUST be one native-width SIMD kernel (`r[io.simd_portable]`) that looks each nibble up with a byte shuffle against the 16-entry table (`swizzle_dyn_within_blocks`: `pshufb` on x86, `tbl` on aarch64) and interleaves the high- and low-nibble results back into read order, decoding `2 × LEN` bases from `LEN` packed bytes per iteration. The remaining whole bytes MUST go through the pair table, and an odd final base through the high nibble of the last byte. The ASCII (`=ACMGRSVTWYHKDBN`) and `Base`-valued decoders MUST be this one kernel with a different table.

r[seq.decode_dispatch+2]
Runtime dispatch MUST go through `dispatch!(Level::new(), ..)`, which selects the widest level the CPU supports: AVX-512, AVX2, SSE4.2 or SSE2 on x86_64, NEON on aarch64, and the scalar `Fallback` elsewhere.

## Encoding (ASCII → 4-bit)

r[seq.encode_scalar]
The scalar encoder MUST use a 256-entry lookup table mapping ASCII bytes to 4-bit codes. Unknown bases MUST map to 15 (`N`). Two bases MUST be packed per byte (high nibble first).

## Correctness

r[seq.simd_scalar_equivalence]
The SIMD decoder MUST produce byte-identical output to the scalar decoder for all inputs. This MUST be verified by property-based tests covering arbitrary sequence lengths and byte patterns, including SIMD boundary lengths (0, 1, 15, 16, 31, 32, 33, 63, 64, 65, 128).
