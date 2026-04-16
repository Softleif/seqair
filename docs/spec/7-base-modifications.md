# Base Modifications (MM/ML Tags)

Base modifications encode epigenetic marks (e.g. 5mC methylation, 6mA, 5hmC) on aligned reads. The SAM 4.5 specification defines two auxiliary tags for this:

- **MM** (Z-type string): modification type and positions, e.g. `C+m,0,2,5;` means "cytosine with 5mC modification at the 1st, 4th, and 10th unmodified C in the query sequence"
- **ML** (B:C array): modification likelihoods (0–255 scale, 255 = highest confidence), one per modification call in MM

These tags are critical for bisulfite/TAPS methylation pipelines and nanopore base modification calling. seqair's TAPS pipeline already writes MM/ML tags during BAM rewriting (see `r[bam.owned_record.aux_array_encoding]`).

> **Sources:** MM/ML tag format is defined in [SAM1] §1.7 (Base modifications). htslib's `hts_base_mod_state` provides the reference parser. The modification type codes (m = 5mC, h = 5hmC, etc.) are defined in the SAM spec's base modification table and registered ChEBI ontology codes.

## Current support: raw pass-through

r[base_mod.passthrough]
MM and ML tags MUST be preserved as raw auxiliary bytes through all I/O paths (BAM read, BAM write, CRAM decode). No parsing or validation of the MM string or ML array contents is performed during pass-through. This ensures that files with base modification annotations are not silently corrupted.

r[base_mod.passthrough.bam]
BAM read → RecordStore → OwnedBamRecord → BAM write MUST preserve MM (Z-type) and ML (B:C array) tags byte-for-byte.

r[base_mod.passthrough.cram]
CRAM decode → RecordStore MUST preserve MM and ML tags. CRAM encodes these via the tag encoding map (`r[cram.compression.tag_encodings]`); the decoder MUST handle arbitrary-length Z and B:C values without truncation.

r[base_mod.passthrough.coexistence]
MM/ML tags MUST coexist with other auxiliary tags (NM, RG, MD, XS, etc.) without interference. Adding, removing, or updating other tags via `AuxData` methods MUST NOT corrupt MM/ML data.

## Future: structured modification API

The following rules describe the planned structured API for base modifications. They are not yet implemented.

r[base_mod.parse_mm]
A `parse_mm(mm_string: &[u8], canonical_base: Base) -> Result<Vec<ModificationCall>>` function MUST parse the MM tag string into structured modification calls. Each call specifies:
- The canonical base (C, G, A, T)
- The modification type code (single char like `m` for 5mC, `h` for 5hmC, or ChEBI code)
- The strand (`+` for top strand, `-` for bottom strand)
- A list of position deltas (skip counts relative to unmodified bases of the canonical type in the query sequence)

The parser MUST handle:
- Multiple modification types per base (e.g. `C+m,0,2;C+h,1,3;`)
- Implicit vs explicit mode (trailing `.` means unlisted positions are unmodified; without `.`, they are unknown)
- ChEBI numeric codes (e.g. `C+27551,0,2;` for 5mC by ChEBI ID)

r[base_mod.resolve_positions]
Given parsed MM deltas and the query sequence, a `resolve_positions(calls: &[ModificationCall], seq: &[Base]) -> Vec<(usize, ModType, u8)>` function MUST convert relative skip-counts into absolute query positions. For each modification call:
1. Iterate through the query sequence
2. Count only bases matching the canonical base
3. At each skip-count boundary, record the absolute position

The output is a list of `(query_position, modification_type, probability)` tuples, where probability comes from the corresponding ML array entry.

r[base_mod.pileup_integration]
The pileup engine SHOULD support exposing base modifications at each column position. For each alignment in a pileup column that has MM/ML tags, the engine SHOULD provide:
- Whether the base at the current `qpos` has a modification call
- The modification type and probability

This enables methylation-aware pileup visualization and variant calling. The exact API (e.g. extending `PileupOp` or a separate accessor) is TBD.

r[base_mod.validation]
When parsing MM tags, the implementation MUST validate:
- The canonical base character is one of A, C, G, T (case-insensitive)
- The strand is `+` or `-`
- Position deltas are non-negative integers
- The total number of modification calls across all MM entries equals the ML array length

Validation failures MUST return typed errors, not panics.

r[base_mod.reverse_complement]
For reverse-strand reads (flag 0x10), modification positions in MM refer to the reverse-complemented query sequence. The implementation MUST account for strand when resolving positions: reverse-strand reads have their sequence stored in reverse-complement form in BAM, so the MM positions index into the stored (already-complemented) sequence directly.
