# Pileup Indel Reporting

## Motivation

A variant caller inspects each pileup column and asks: "what does each read show at this reference position?" For SNPs, the answer is a base (A/C/G/T). But for indels, the answer is richer:

- **Deletion**: the read's CIGAR has a D operation spanning this position. The read is aligned across this position but has no query base here — the reference base was deleted. A variant caller counts these as evidence for a deletion allele.
- **Insertion**: after a matched position, the read's CIGAR has an I operation. The read has extra bases not present in the reference, inserted between two reference positions. A variant caller needs to know where insertions occur and how long they are to call insertion alleles.
- **Reference skip (intron)**: the read's CIGAR has an N operation. This represents a spliced alignment (RNA-seq intron), not a variant. The variant caller must distinguish this from a deletion to avoid false calls.

Without indel reporting, the pileup engine silently drops reads at deletion positions (they have no `qpos`) and insertions are completely invisible. This makes indel calling impossible from the pileup alone — the caller would have to re-parse every read's CIGAR independently, defeating the purpose of the pileup abstraction.

### Design choice: type-safe enum vs. htslib flags

htslib's `bam_pileup1_t` uses flat fields: `is_del`, `is_refskip`, `indel`, and `qpos` are all available on every alignment regardless of state. This means a caller can accidentally read `qpos` from a deletion alignment (where it's meaningless) or forget to check `is_del` before using `base`.

We use a Rust enum (`PileupOp`) where `qpos`, `base`, and `qual` are only accessible in the `Match` and `Insertion` variants. A `Deletion` alignment cannot produce a `base` value without an explicit match. This makes the calling code at the use site visually obvious about which cases it handles, and the compiler rejects code that forgets to handle deletions.

> **Sources:** Seqair-specific design. Inspired by [htslib] `bam_pileup1_t` which provides `is_del`, `is_refskip`, and `indel` fields, but redesigned as a Rust enum for compile-time safety. CIGAR operation semantics from [SAM1] §1.4. See [references.md](99-references.md).

## PileupOp enum

> *[SAM1] §1.4 — CIGAR operations: M/=/X consume ref+query, D consumes ref only, N is ref-skip, I consumes query only*

r[pileup_indel.op_enum]
Each alignment in a pileup column MUST carry a `PileupOp` that describes the CIGAR state at that reference position:

- **`Match`**: the read has a base aligned at this position (M, =, or X CIGAR op). Carries `qpos`, `base`, and `qual`.
- **`Insertion`**: same as Match, but an insertion of `insert_len` query bases follows before the next reference position. Carries `qpos`, `base`, `qual`, and `insert_len`.
- **`Deletion`**: the read spans this position via a deletion (D op). No query base exists. Does not carry `qpos`, `base`, or `qual`.
- **`RefSkip`**: the read spans this position via a reference skip (N op, e.g., intron). No query base exists. Does not carry `qpos`, `base`, or `qual`.

r[pileup_indel.type_safety]
`qpos`, `base`, and `qual` MUST only be accessible when the op is `Match` or `Insertion`. The type system (enum variants) MUST enforce this at compile time. Callers MUST NOT be able to read a `base` value from a `Deletion` or `RefSkip` alignment without an explicit match/conversion.

## Inclusion in columns

r[pileup_indel.deletions_included]
Reads with a deletion at the current reference position MUST be included in the pileup column with `PileupOp::Deletion`. They MUST NOT be silently excluded. This supersedes the original `r[pileup.qpos_none]` exclusion rule.

r[pileup_indel.refskips_included]
Reads with a reference skip (N op) at the current position MUST be included in the column with `PileupOp::RefSkip`.

r[pileup_indel.depth_includes_all]
`PileupColumn::depth()` MUST count all alignments including deletions and ref-skips. This matches htslib's depth counting behavior.

## Insertion detection

> *Insertions in CIGAR (I op) consume query but not reference. They occur between two reference positions.*

r[pileup_indel.insertion_at_last_match]
An insertion MUST be reported at the last reference-consuming, query-consuming position (M/=/X) before the I op in the CIGAR. The `insert_len` field gives the total length of inserted query bases. The caller can access the inserted bases via the read's sequence at `qpos + 1 .. qpos + 1 + insert_len`.

r[pileup_indel.insertion_len]
If multiple consecutive I ops appear (which is invalid per SAM spec but may occur), their lengths MUST be summed into a single `insert_len`.

r[pileup_indel.no_orphan_insertions]
An insertion that is not preceded by a M/=/X op within the same CIGAR (e.g., an insertion after a deletion: `D I M`) MUST NOT be reported as an `Insertion` op. The insertion's query bases are still consumed, but the insertion event is not surfaced because there is no anchor match position to attach it to.

## Convenience accessors

r[pileup_indel.accessors]
`PileupAlignment` MUST provide convenience methods:
- `qpos() -> Option<usize>`: returns `Some(qpos)` for Match/Insertion, `None` for Deletion/RefSkip.
- `base() -> Option<Base>`: returns `Some(base)` for Match/Insertion, `None` for Deletion/RefSkip.
- `qual() -> Option<u8>`: returns `Some(qual)` for Match/Insertion, `None` for Deletion/RefSkip.
- `is_del() -> bool`: true for Deletion.
- `is_refskip() -> bool`: true for RefSkip.
- `insert_len() -> u32`: returns the insertion length for Insertion, 0 for all other ops.
- `op() -> &PileupOp`: returns a reference to the op enum.

## Deduplication interaction

r[pileup_indel.dedup_with_deletions]
Overlapping-pair deduplication MUST handle Deletion, RefSkip, and Insertion alignments. Resolution rules in priority order:

1. When one mate has a base (Match/Insertion) and the other does not (Deletion/RefSkip), the mate with a base MUST be kept.
2. When both have bases and show different bases, the first-in-template read MUST be kept.
3. When both have the same base but one carries an Insertion op, the Insertion MUST be kept (it carries strictly more information for indel calling).
4. When both have the same base and same insertion status, the first encountered (by arena order) MUST be kept.
5. When neither has a base (both Deletion/RefSkip), the first encountered MUST be kept.

## Backward compatibility

r[pileup_indel.htslib_compat_update]
With deletions and ref-skips included in columns, the engine MUST produce depth values matching htslib's `bam_plp_auto` at all positions, including those within deletions and ref-skips. This strengthens the existing `r[pileup.htslib_compat]` rule.
