# Overlapping Pair Deduplication

In paired-end sequencing, two reads (mates) are generated from opposite ends of the same DNA fragment. When the fragment is shorter than twice the read length, the mates overlap — both reads cover the same genomic positions in the middle. Without deduplication, these overlapping positions would be counted twice in the pileup, inflating the apparent coverage and biasing variant allele frequencies.

> **Sources:** Seqair-specific design with no upstream spec counterpart. Mate detection relies on QNAME matching as defined in [SAM1] §1.4. The first-in-template / second-in-template distinction uses FLAG bits 0x40 and 0x80 defined in [SAM1] §1.4. See [references.md](references.md).

For example, with 150bp reads from a 200bp fragment, positions 50–150 are covered by both mates. At each of these positions, the pileup would show two observations from the same molecule. Deduplication removes one mate at each overlapping position, keeping the observation count honest.

## Opt-in behavior

r[dedup.opt_in]
Overlapping pair deduplication MUST be opt-in via `set_dedup_overlapping()`. When not enabled, the pileup engine MUST yield all alignments without deduplication.

## Mate detection

Mates are identified by sharing the same read name (qname). The engine pre-computes a mate lookup table when deduplication is enabled, so the per-position dedup check is a fast hash lookup rather than a string comparison.

r[dedup.mate_detection]
Mates are detected by matching qnames. When `set_dedup_overlapping()` is called, the engine MUST build a mate index that maps each record to its mate's index (if both are in the arena). Records with no mate in the arena are unaffected.

r[dedup.mate_pairs_only]
Only the first two records sharing a qname are treated as mates. If three or more records share the same qname (supplementary alignments, etc.), only the first pair is linked; additional records are treated as unpaired for dedup purposes.

## Resolution rule

When both mates are present at a position, one must be removed. The rule depends on whether the mates agree on the base call:

r[dedup.resolution_same_base]
When both mates show the same base at a position, the first encountered (by arena order) is kept and the other is removed.

r[dedup.resolution_different_base]
When mates show different bases at a position, the first-in-template read is kept and the second-in-template read is removed. This preserves the primary observation.

## Scope

r[dedup.per_position]
Deduplication MUST be applied per-position, not globally. A read that is removed at one overlapping position MUST still appear at positions where its mate is not active. Both mates contribute normally outside their overlap region.

## Non-interference

r[dedup.filter_independent]
Deduplication MUST be applied after the read filter and after qpos computation. A read that was excluded by the filter MUST NOT participate in dedup. A read that has no qpos at a position (deletion/skip) MUST NOT participate in dedup at that position.

r[dedup.max_depth_independent]
Deduplication MUST be applied before max_depth truncation. The depth after dedup may be below max_depth even if the raw depth exceeded it.
