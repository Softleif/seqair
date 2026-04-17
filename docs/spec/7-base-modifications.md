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

## Planned: structured modification API

The following rules describe the planned structured API for base modifications. They are not yet implemented.

### Design choices

**Type-based API, not free functions.** A `BaseModState` struct holds the parsed MM + ML data for a single record and provides query methods. This mirrors seqair's "types are the primary abstraction" philosophy and parallels `CigarMapping` (pre-extracted CIGAR data) and `PileupAlignment` (pre-extracted alignment fields).

**Lazy parsing with caching.** Modification data is NOT parsed when a record enters the pileup active set. Instead, `BaseModState` is constructed on first access (e.g. when a caller asks for modification data at a column position). The parsed state is cached for subsequent queries on the same record. This avoids overhead when modifications aren't needed (the common case for non-methylation workflows).

**Dual coordinate access.** Modifications can be queried in query-space (matching the stored BAM sequence, zero-cost) or reference-space (mapped through the CIGAR). The pileup integration uses query-space internally since it already has `qpos`.

### BaseModState

r[base_mod.state]
`BaseModState` MUST be constructable from a record's MM string and ML array. It parses all modification entries in the MM string (not filtered by canonical base) and pairs them with the corresponding ML probabilities. The struct is immutable after construction.

r[base_mod.parse_mm]
The MM parser MUST handle the full SAM 4.5 modification string syntax:
- Multiple canonical bases and modification types (e.g. `C+m,0,2;C+h,1,3;A+a,0;`)
- Mode markers after the modification code: `.` (explicit, unlisted positions are definitively unmodified), `?` (ambiguous, unlisted positions are explicitly marked as unknown), or absent (implicit, unlisted positions have unknown status). The API MUST preserve this three-way distinction; `?` and absent produce identical `is_unmodified` results but are retained separately in the parsed mode field.
- ChEBI numeric codes (e.g. `C+27551,0,2;` for 5mC by ChEBI ID)
- Combined modification codes in a single entry (e.g. `C+mh,0,2;`), where the ML array carries one value per (position, mod-code) pair in the order the codes appear. For an entry with N combined codes and K position deltas, the entry consumes `N * K` values from the ML array.
- Both `+` (top/forward strand) and `-` (bottom/reverse strand) modification strands

r[base_mod.resolve_positions]
`BaseModState` MUST resolve the relative skip-counts in MM into absolute query-sequence positions at construction time. For each modification entry:
1. Iterate through the query sequence (provided at construction)
2. Count only bases matching the canonical base for that entry
3. At each skip-count boundary, record the absolute query position and pair it with the corresponding ML probability

The resolved positions are stored internally as a sorted list for efficient lookup.

### Query API

r[base_mod.query_qpos]
`BaseModState` MUST provide `mod_at_qpos(qpos: usize) -> Option<&[Modification]>` to look up modifications at a given query-sequence position. Returns `None` if no modification is called at that position. Multiple modification types at the same position (e.g. 5mC and 5hmC on the same C) MUST all be returned.

r[base_mod.query_refpos]
`BaseModState` MUST provide `mod_at_ref_pos(ref_pos: i64, cigar: &CigarMapping) -> Option<&[Modification]>` to look up modifications at a reference position. This maps the reference position to a query position via the CIGAR, then delegates to `mod_at_qpos`. Returns `None` if the reference position falls in a deletion/ref-skip (no qpos) or has no modification.

r[base_mod.modification_struct]
Each `Modification` MUST carry:
- `mod_type`: the modification type (single-char code or ChEBI number)
- `probability`: the ML value (u8, 0–255)
- `canonical_base`: which base is modified (A, C, G, T)
- `strand`: `+` or `-`

r[base_mod.implicit_explicit]
`BaseModState` MUST track the mode marker per entry as one of `Implicit` (no marker), `Unmodified` (`.`), or `Ambiguous` (`?`). A `is_unmodified(qpos: usize, canonical_base: Base) -> Option<bool>` accessor MUST return:
- `Some(true)` if the position's canonical base has an `Unmodified`-mode entry and the position is not listed (definitively unmodified)
- `Some(false)` if the position has a modification call
- `None` if the position has no mod call and all entries for its canonical base are `Implicit` or `Ambiguous` (unknown status)

### Pileup integration

r[base_mod.pileup_integration]
*Deferred to a follow-up milestone.* The pileup engine will eventually expose modifications via a separate accessor on `PileupColumn` (NOT by extending `PileupOp`, which has a ≤16 byte size constraint):

```
column.modification_at(alignment_idx: usize) -> Option<&[Modification]>
```

This will return the modification(s) at the current column's reference position for the given alignment, using `mod_at_qpos` internally, with `BaseModState` constructed lazily and cached per active record. The initial `BaseModState` implementation is usable directly by callers (e.g. rastair) without pileup integration.

### Reverse complement

r[base_mod.reverse_complement]
Per SAM §1.7, MM position lists are always relative to the **original, unreversed** read sequence (5′→3′ of the sequenced molecule), regardless of alignment orientation. BAM stores the sequence reverse-complemented for reverse-strand alignments, so MM positions do NOT directly index into the stored BAM sequence for such reads.

`BaseModState` MUST take an `is_reverse` flag at construction and resolve MM positions into stored-BAM-sequence coordinates:
- If `is_reverse == false`: count occurrences of the canonical base in the stored sequence from index 0 upward. The resolved `qpos` is the stored-sequence index of the matched base.
- If `is_reverse == true`: count occurrences of the **complement** of the canonical base (A↔T, C↔G) in the stored sequence from the last index downward. The resolved `qpos` is the stored-sequence index of the matched base. This is equivalent to: reverse-complement the stored sequence to reconstruct the original, walk it 5′→3′ counting the canonical base, and map the resulting original-index `i` back to stored-index `seq_len - 1 - i`.

After resolution, all public queries (`mod_at_qpos`, `mod_at_ref_pos`) operate in stored-BAM-sequence coordinates. `mod_at_ref_pos` delegates to the CIGAR mapping (which also operates in stored-sequence coordinates).

Callers who need to reason about the original biological strand (e.g. "this is a methylated C on the bottom strand") can combine the `strand` field from `Modification` with the read's reverse flag.

### Validation

r[base_mod.validation]
When parsing MM tags, the implementation MUST validate:
- The canonical base character is one of A, C, G, T (case-insensitive)
- The strand is `+` or `-`
- Position deltas are non-negative integers
- The total number of modification calls across all MM entries equals the ML array length

Validation failures MUST return typed errors, not panics.
