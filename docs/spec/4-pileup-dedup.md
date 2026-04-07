# Overlapping Pair Deduplication

> **Status: Not implemented in seqair.** This spec documents the correct design
> for overlap dedup based on lessons learned integrating with htslib. The
> implementation was removed from seqair's pileup engine because correct dedup
> requires running **after** per-position quality filtering (see
> [Quality filtering order](#quality-filtering-order) below). Consumers like
> rastair handle dedup themselves using a name-based collector after
> quality/masking filtering, matching htslib's approach.

In paired-end sequencing, two reads (mates) are generated from opposite ends of the same DNA fragment. When the fragment is shorter than twice the read length, the mates overlap — both reads cover the same genomic positions in the middle. Without deduplication, these overlapping positions would be counted twice in the pileup, inflating the apparent coverage and biasing variant allele frequencies.

> **Sources:** Mate detection relies on QNAME matching as defined in [SAM1] §1.4. The first-in-template / second-in-template distinction uses FLAG bits 0x40 and 0x80 defined in [SAM1] §1.4. See [References](./99-references.md).

For example, with 150bp reads from a 200bp fragment, positions 50–150 are covered by both mates. At each of these positions, the pileup would show two observations from the same molecule. Deduplication removes one mate at each overlapping position, keeping the observation count honest.

## Mate detection

Mates are identified by sharing the same read name (qname).

r[dedup.mate_detection]
Mates MUST be detected by matching qnames within the same pileup column. Records with no mate at the current position are unaffected.

r[dedup.mate_pairs_only]
Only the first two records sharing a qname are treated as mates. If three or more records share the same qname (supplementary alignments, etc.), only the first pair is linked; additional records are treated as unpaired for dedup purposes.

## Resolution rule

When both mates are present at a position, one must be removed. Since both mates originate from the same DNA molecule, either could serve as the observation. The resolution rule must be **deterministic regardless of read iteration order** — different pileup engines (htslib, seqair) may iterate reads in different physical file orders.

r[dedup.resolution]
When both mates have a base at a position: keep first-in-template (FLAG 0x40), remove second-in-template (FLAG 0x80). When both have the same template flag (rare, e.g. supplementary alignments), keep the first-encountered. This rule applies whether the mates agree or disagree on the base call.

> **Design note:** An earlier version used different rules for same-base
> (keep first-encountered) vs different-base (keep first-in-template). The
> "keep first-encountered" rule is iteration-order-dependent, causing different
> pileup engines to produce different results for the same data. Using
> first-in-template consistently eliminates this source of non-determinism.

r[dedup.no_base_vs_base]
When one mate has a base at the position and the other does not (deletion or ref-skip), the mate with the base MUST be kept.

## Scope

r[dedup.per_position]
Deduplication MUST be applied per-position, not globally. A read that is removed at one overlapping position MUST still appear at positions where its mate is not active. Both mates contribute normally outside their overlap region.

## Quality filtering order

r[dedup.after_quality_filter]
Deduplication MUST be applied **after** per-position quality and masking filters. A read that fails quality filtering at a specific position MUST NOT participate in dedup at that position.

> **Rationale:** htslib applies quality/masking filters before overlap dedup.
> If dedup runs first, a mate pair where one read has low base quality at a
> position gets deduped (removing the good mate), then quality filtering also
> removes the bad mate — both are lost. With filter-then-dedup, the low-quality
> mate is filtered out first, the pair is never matched, and the good mate
> survives. This was the root cause of depth discrepancies between seqair's
> engine-level dedup and htslib's post-filter dedup.
>
> This constraint means overlap dedup cannot be implemented inside the pileup
> engine, which operates before per-position quality filtering. It must be
> done by the consumer after extracting and filtering the column's alignments.

r[dedup.max_depth_independent]
Deduplication SHOULD be applied before max_depth truncation. The depth after dedup may be below max_depth even if the raw depth exceeded it.
