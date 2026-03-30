All 4 agents are done. Here's the consolidated report.

---

## Proptest Audit: Tautological vs. Valid Tests

### Summary

Across 20 files and **~50 proptests**, I found **~22 tautological or weak tests** and **~28 valid tests**.

---

### Tautological Tests (reimplement production logic or test trivial identities)

| File | Test | Issue |
|---|---|---|
| **probability.rs** | `proptest_inverted_stays_in_range` | Checks `[0,1]` bounds — trivially true for `1-x` on `[0,1]`. Should check `p + p.inverted() ≈ 1.0` |
| **probability.rs** | `proptest_probability_phred_roundtrip` | Test reimplements `10^(-Q/10)` — same formula as production |
| **phred.rs** | `proptest_from_phred_roundtrip` | Tests IEEE 754 integer exactness, not code logic |
| **phred.rs** | `proptest_as_int_always_in_range` | Redundant duplicate of `proptest_from_phred_produces_valid` |
| **phred.rs** | `proptest_phred_probability_roundtrip` | Same formula both sides — algebraic inverse |
| **base.rs** | `proptest_from_u8_roundtrip` | Tests repr consistency, doesn't verify case-folding |
| **base.rs** | `proptest_from_str_roundtrip` | Byte-value proxy, not enum identity |
| **region_string.rs** | `proptest_roundtrip_random_string` | Pure `Display(parse(s)) == s` roundtrip, no oracle |
| **region_string.rs** | `proptest_roundtrip_region_string` | Field checks valid, but `to_string()` roundtrip is redundant |
| **varint.rs** | `itf8_roundtrip` | Encode/decode are inverses in same file |
| **varint.rs** | `itf8_read_matches_decode` | Same encoder; only proves two decoders agree |
| **varint.rs** | `itf8_cursor_matches_decode` | Same + relies on zero-fill of unused buffer |
| **varint.rs** | `ltf8_read_matches_decode` | Test encoder copies same shift/mask formulas as production |
| **encoding.rs** | `beta_roundtrip` | `bits` always 8, `offset` always 0 — tests identity read |
| **block.rs** | `raw_block_roundtrip` | `build_test_block` mirrors `parse_block`'s layout |
| **bitstream.rs** | `roundtrip_byte_values` | "8 MSB-first bits = the byte" — definitionally true |
| **bitstream.rs** | `roundtrip_u16_values` | Same — oracle is the definition of big-endian |
| **bgzf.rs** | `virtual_offset_roundtrips` | Pack/unpack roundtrip, doesn't pin bit layout to BAI spec |
| **bgzf.rs** | `ordering_matches_file_order` | Re-derives unsigned integer ordering by hand |
| **seq_codec.rs** | `encode_decode_roundtrip` | Roundtrip only — doesn't verify nibble values per BAM spec |
| **fasta.rs** | `byte_offset_*` (4 tests) | All derive expected values from the same `offset + (pos/linebases)*linewidth + (pos%linebases)` formula |
| **pileup.rs** | `depth_equals_covering_records` | Oracle uses same `pos + len - 1` as production |
| **pileup.rs** | `every_covered_position_yields_a_column` | Same interval predicate as engine |
| **perf.rs** | `reused_alignment_vec_produces_correct_depth` | Duplicate of pileup pattern |
| **perf.rs** | `precomputed_matches_indels_consistent_with_calc` | Sums same op codes as production; misses `=`/`X` ops |
| **sam_proptests.rs** | `cigar_roundtrip_preserves_end_pos` | Helper functions `cigar_query_len`/`cigar_ref_len` reimplement production logic |
| **cigar.rs** | `none_at_deletions_some_at_matches` | Test loop reuses same op-code constants |

### Valid Tests (independent oracles or genuine invariants)

| File | Test | Why it's good |
|---|---|---|
| **rms.rs** | All 6 proptests | Mathematical invariants (constant=value, non-negative, RMS>=mean, RMS<=max) |
| **probability.rs** | `valid_range_accepted`, `above/below_rejected`, `from_str_roundtrip` | Boundary/rejection tests with independent oracles |
| **phred.rs** | `proptest_from_phred_produces_valid` | Tests clamping contract |
| **base.rs** | `from_str_rejects_multi_char`, `from_str_never_panics` | Error/panic safety |
| **base.rs** | `from_ascii_vec_equivalence/all_byte_values` | SIMD vs scalar differential test |
| **varint.rs** | `ltf8_small_values_match_itf8` | Cross-format spec invariant |
| **encoding.rs** | `external_byte_roundtrip` | Cursor advancement (weak but valid) |
| **compression_header.rs** | `substitution_matrix_code_in_bounds` | Safety invariant (weak but valid) |
| **block.rs** | `gzip_block_roundtrip` | Independent codec (libdeflater) as oracle |
| **bitstream.rs** | `split_read_equals_combined` | Non-trivial: position state consistency across partial reads |
| **seq_codec.rs** | `simd_matches_scalar_arbitrary` | Genuine differential test |
| **cigar.rs** | `qpos_monotonically_increasing`, `pure_match_qpos_is_offset` | Real structural invariants |
| **sam_proptests.rs** | `qual_decode_roundtrip` | Phred+33 standard as oracle |
| **compare_bgzf_crate.rs** | `proptest_both_decompressors_match` | Third-party library as oracle — strongest test |
| **region_buf.rs** | All 5 proptests | Independent encoder (libdeflater) + CRC corruption test |
| **fasta.rs** | `plain_fasta_roundtrip`, `bgzf_random_regions_match_htslib`, `fork_matches_original`, `gzi_translate_before_first_entry`, `gzi_translate_rejects_overflow` | htslib oracle, safety checks |
| **pileup.rs** | `filter_by_flags`, `max_depth`, `qpos_for_simple_cigar`, `no_empty_columns`, `reads_available_at_lower_coverage_despite_cap` | Behavioral properties |
| **perf.rs** | `binary_search_matches_linear_for_many_ops` | Monotonicity + intron/exon invariants |
| **dedup.rs** | Both tests | Independent reference (non-deduped engine) |

---

### Top Priority Fixes

The most impactful tautological tests to fix, ranked by risk:

1. **pileup.rs `depth_equals_covering_records`** + **perf.rs duplicate** — Use hand-verified fixture data or an independent interval-tree oracle instead of the same `pos + len - 1` formula
2. **fasta.rs `byte_offset_*`** (4 tests) — Validate against actual file bytes: seek to `byte_offset(pos)` in a real `.fa` and assert the character matches `bases[pos]`
3. **sam_proptests.rs `cigar_roundtrip_preserves_end_pos`** — Replace in-test CIGAR parsers with known-correct test vectors from the SAM spec
4. **varint.rs** (3 ITF8 roundtrips) — Assert expected wire bytes per the CRAM spec table, not just encode/decode consistency
5. **perf.rs `precomputed_matches_indels`** — Generate from full `{M,I,D,N,S,H,=,X}` op set; this test currently **can't detect** if `=`/`X` ops are mishandled
6. **probability.rs `proptest_inverted_stays_in_range`** — Trivial fix: assert `*p + *p.inverted() ≈ 1.0`

Want me to fix any of these?
