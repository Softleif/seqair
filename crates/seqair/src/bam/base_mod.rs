//! Parse SAM §1.7 base-modification tags (`MM`, `ML`) into a queryable
//! [`BaseModState`]. Positions are resolved once at construction into
//! stored-BAM-sequence coordinates; queries are zero-alloc lookups.
//!
//! MM positions in SAM are always relative to the **original, unreversed** read
//! sequence. For reverse-strand alignments BAM stores the reverse-complemented
//! sequence, so `BaseModState::parse` takes an `is_reverse` flag and resolves
//! positions into stored-sequence coordinates up front — subsequent queries
//! never re-read the strand.

use seqair_types::Base;
use thiserror::Error;

use super::cigar::{CigarMapping, CigarPosInfo};

// r[impl base_mod.modification_struct]
/// Modification type code: a single-character SAM code (`m`, `h`, `a`, …) or a
/// numeric `ChEBI` ontology id.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ModType {
    /// Single-character modification code (e.g. `m` for 5mC, `h` for 5hmC).
    Code(u8),
    /// Numeric `ChEBI` id (e.g. 27551 for 5mC).
    ChEBI(u32),
}

/// Strand marker from the MM entry header.
///
/// This is the strand annotated in the MM tag itself, independent of whether
/// the read itself is aligned on the reverse strand.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ModStrand {
    Plus,
    Minus,
}

/// Mode marker following the modification code(s) in an MM entry.
///
/// Determines how to treat positions not explicitly listed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ModMode {
    /// No marker. Unlisted positions have unknown modification status.
    Implicit,
    /// `.` marker. Unlisted positions are definitively unmodified.
    Unmodified,
    /// `?` marker. Unlisted positions are explicitly unknown.
    Ambiguous,
}

/// A single modification call on a single base.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Modification {
    pub mod_type: ModType,
    pub probability: u8,
    pub canonical_base: Base,
    pub strand: ModStrand,
}

/// Errors produced by [`BaseModState::parse`].
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum BaseModError {
    #[error("invalid canonical base in MM tag: {0:#x}")]
    InvalidCanonicalBase(u8),
    #[error("canonical base 'N' is not supported for modification resolution")]
    UnsupportedUnknownBase,
    #[error("invalid strand byte in MM tag: {0:#x}")]
    InvalidStrand(u8),
    #[error("missing modification code in MM entry")]
    MissingModCode,
    #[error("invalid modification code byte: {0:#x}")]
    InvalidModCode(u8),
    #[error("invalid ChEBI numeric id in MM tag")]
    InvalidChebi,
    #[error("invalid position delta in MM tag")]
    InvalidDelta,
    #[error("MM entry missing terminating ';'")]
    MissingTerminator,
    #[error("ML array length {ml} does not match total MM call count {mm}")]
    MlLengthMismatch { ml: usize, mm: usize },
    #[error(
        "MM delta {delta} exceeds remaining count of canonical base {base:?} in sequence \
         (is_reverse={is_reverse})"
    )]
    DeltaOutOfRange { delta: u32, base: Base, is_reverse: bool },
}

// r[impl base_mod.state]
/// Resolved modification calls for a single record, stored in
/// stored-BAM-sequence coordinates.
#[derive(Debug, Clone, Default)]
pub struct BaseModState {
    /// Parallel to `mods`: sorted ascending, ties grouped (same-qpos entries
    /// are adjacent in insertion order, which preserves combined-code order).
    qpos: Vec<u32>,
    mods: Vec<Modification>,
    /// For each canonical base (A,C,G,T → index 0..=3): is there any entry
    /// whose mode is `Unmodified`? Used by `is_unmodified` to decide whether
    /// an unlisted position can be called definitively unmodified.
    any_unmodified_mode: [bool; 4],
}

impl BaseModState {
    // r[impl base_mod.parse_mm]
    // r[impl base_mod.resolve_positions]
    // r[impl base_mod.reverse_complement]
    // r[impl base_mod.validation]
    /// Parse MM/ML tag values and resolve positions against `seq`.
    ///
    /// `seq` is the sequence as stored in BAM (already reverse-complemented
    /// for reverse-strand reads). `is_reverse` is the read's reverse-strand
    /// flag (`flag & 0x10 != 0`).
    pub fn parse(
        mm: &[u8],
        ml: &[u8],
        seq: &[Base],
        is_reverse: bool,
    ) -> Result<Self, BaseModError> {
        let mut state = BaseModState::default();
        // Collect (qpos, Modification) in insertion order so that
        // combined-code entries stay adjacent at the same qpos; we sort at
        // the end by qpos with a stable sort.
        let mut pending: Vec<(u32, Modification)> = Vec::new();
        let mut ml_cursor: usize = 0;

        // Tolerate a trailing NUL (Z-type strings are sometimes handed in with
        // the NUL included) and an optional trailing empty segment.
        let mm = mm.strip_suffix(&[0]).unwrap_or(mm);

        let mut rest = mm;
        while !rest.is_empty() {
            let (entry, next) = split_entry(rest)?;
            rest = next;
            if entry.is_empty() {
                // Stray `;` — skip.
                continue;
            }
            parse_entry(entry, seq, is_reverse, ml, &mut ml_cursor, &mut pending, &mut state)?;
        }

        if ml_cursor != ml.len() {
            return Err(BaseModError::MlLengthMismatch { ml: ml.len(), mm: ml_cursor });
        }

        // Stable sort by qpos preserves same-qpos insertion order.
        pending.sort_by_key(|&(qp, _)| qp);
        state.qpos.reserve(pending.len());
        state.mods.reserve(pending.len());
        for (qp, m) in pending {
            state.qpos.push(qp);
            state.mods.push(m);
        }
        Ok(state)
    }

    /// Number of resolved modification calls.
    pub fn len(&self) -> usize {
        self.mods.len()
    }

    pub fn is_empty(&self) -> bool {
        self.mods.is_empty()
    }

    // r[impl base_mod.query_qpos]
    /// All modifications called at the given stored-sequence position, or
    /// `None` if no call is present.
    pub fn mod_at_qpos(&self, qpos: usize) -> Option<&[Modification]> {
        let target = u32::try_from(qpos).ok()?;
        let lo = self.qpos.partition_point(|&q| q < target);
        let hi = self.qpos.partition_point(|&q| q <= target);
        if lo == hi {
            return None;
        }
        self.mods.get(lo..hi)
    }

    // r[impl base_mod.query_refpos]
    /// Modifications at a reference position, mapped through `cigar`. Returns
    /// `None` if the reference position is in a deletion/ref-skip or has no
    /// modification call.
    pub fn mod_at_ref_pos(
        &self,
        ref_pos: seqair_types::Pos0,
        cigar: &CigarMapping,
    ) -> Option<&[Modification]> {
        match cigar.pos_info_at(ref_pos)? {
            CigarPosInfo::Match { qpos } | CigarPosInfo::Insertion { qpos, .. } => {
                self.mod_at_qpos(qpos as usize)
            }
            CigarPosInfo::Deletion { .. }
            | CigarPosInfo::RefSkip
            | CigarPosInfo::ComplexIndel { .. } => None,
        }
    }

    // r[impl base_mod.implicit_explicit]
    /// Whether `qpos` is definitively unmodified, definitively modified, or
    /// unknown. See [`base_mod.implicit_explicit`] in the spec.
    pub fn is_unmodified(&self, qpos: usize, canonical_base: Base) -> Option<bool> {
        if self.mod_at_qpos(qpos).is_some() {
            return Some(false);
        }
        let idx = canonical_base.known_index()?;
        #[allow(clippy::indexing_slicing, reason = "known_index() returns Some(0..=3) for ACGT")]
        let has_unmodified_mode = self.any_unmodified_mode[idx];
        if has_unmodified_mode { Some(true) } else { None }
    }
}

// ---- Parser ----

/// Split off one `;`-terminated entry from the head of `mm`.
fn split_entry(mm: &[u8]) -> Result<(&[u8], &[u8]), BaseModError> {
    let sc = mm.iter().position(|&b| b == b';').ok_or(BaseModError::MissingTerminator)?;
    let (entry, rest) = mm.split_at(sc);
    // Skip the `;`.
    let rest = rest.get(1..).unwrap_or(&[]);
    Ok((entry, rest))
}

/// Parse a single MM entry and push resolved modifications into `pending`.
#[allow(clippy::too_many_arguments, reason = "parser scratch state")]
fn parse_entry(
    entry: &[u8],
    seq: &[Base],
    is_reverse: bool,
    ml: &[u8],
    ml_cursor: &mut usize,
    pending: &mut Vec<(u32, Modification)>,
    state: &mut BaseModState,
) -> Result<(), BaseModError> {
    // [BASE][STRAND][MOD_CODES][MODE]?(,DELTA)*
    let canonical = parse_canonical_base(*entry.first().ok_or(BaseModError::MissingModCode)?)?;
    let strand = parse_strand(*entry.get(1).ok_or(BaseModError::MissingModCode)?)?;
    let body = entry.get(2..).unwrap_or(&[]);

    let (mod_types, mode, deltas_bytes) = parse_body(body)?;

    // Parse the delta list (zero or more comma-separated unsigned integers).
    let deltas = parse_deltas(deltas_bytes)?;

    // ML consumes `deltas.len() * mod_types.len()` values for this entry.
    let total = deltas
        .len()
        .checked_mul(mod_types.len())
        .ok_or(BaseModError::MlLengthMismatch { ml: ml.len(), mm: usize::MAX })?;
    let ml_end = ml_cursor
        .checked_add(total)
        .ok_or(BaseModError::MlLengthMismatch { ml: ml.len(), mm: usize::MAX })?;
    if ml_end > ml.len() {
        return Err(BaseModError::MlLengthMismatch { ml: ml.len(), mm: ml_end });
    }
    #[allow(clippy::indexing_slicing, reason = "bounds checked immediately above")]
    let ml_slice = &ml[*ml_cursor..ml_end];
    *ml_cursor = ml_end;

    // Resolve positions and emit Modifications.
    resolve_and_emit(
        canonical, strand, &mod_types, mode, &deltas, ml_slice, seq, is_reverse, pending, state,
    )?;
    Ok(())
}

fn parse_canonical_base(byte: u8) -> Result<Base, BaseModError> {
    match byte | 0x20 {
        b'a' => Ok(Base::A),
        b'c' => Ok(Base::C),
        b'g' => Ok(Base::G),
        b't' => Ok(Base::T),
        b'n' => Err(BaseModError::UnsupportedUnknownBase),
        _ => Err(BaseModError::InvalidCanonicalBase(byte)),
    }
}

fn parse_strand(byte: u8) -> Result<ModStrand, BaseModError> {
    match byte {
        b'+' => Ok(ModStrand::Plus),
        b'-' => Ok(ModStrand::Minus),
        other => Err(BaseModError::InvalidStrand(other)),
    }
}

/// Parse the body of an entry after `[BASE][STRAND]`: mod codes, optional
/// mode marker, and the `,`-prefixed delta list.
///
/// Returns `(mod_types, mode, remaining_bytes)` where remaining bytes include
/// the leading comma (or are empty if no deltas follow).
fn parse_body(body: &[u8]) -> Result<(Vec<ModType>, ModMode, &[u8]), BaseModError> {
    let first = *body.first().ok_or(BaseModError::MissingModCode)?;
    let (mod_types, after_codes) = if first.is_ascii_digit() {
        // ChEBI numeric id: read digits.
        let end = body.iter().position(|b| !b.is_ascii_digit()).unwrap_or(body.len());
        #[allow(clippy::indexing_slicing, reason = "end <= body.len() by construction")]
        let digits = &body[..end];
        let s = std::str::from_utf8(digits).map_err(|_| BaseModError::InvalidChebi)?;
        let id: u32 = s.parse().map_err(|_| BaseModError::InvalidChebi)?;
        (vec![ModType::ChEBI(id)], body.get(end..).unwrap_or(&[]))
    } else {
        // One or more single-character codes.
        let mut codes = Vec::new();
        let mut i = 0;
        while let Some(&b) = body.get(i) {
            if matches!(b, b',' | b';' | b'.' | b'?') {
                break;
            }
            if !b.is_ascii_alphabetic() {
                return Err(BaseModError::InvalidModCode(b));
            }
            codes.push(ModType::Code(b));
            i = i.saturating_add(1);
        }
        if codes.is_empty() {
            return Err(BaseModError::MissingModCode);
        }
        (codes, body.get(i..).unwrap_or(&[]))
    };

    let (mode, remaining) = match after_codes.first() {
        Some(b'.') => (ModMode::Unmodified, after_codes.get(1..).unwrap_or(&[])),
        Some(b'?') => (ModMode::Ambiguous, after_codes.get(1..).unwrap_or(&[])),
        _ => (ModMode::Implicit, after_codes),
    };
    Ok((mod_types, mode, remaining))
}

fn parse_deltas(bytes: &[u8]) -> Result<Vec<u32>, BaseModError> {
    let mut out = Vec::new();
    let mut rest = bytes;
    while let Some((&first, tail)) = rest.split_first() {
        if first != b',' {
            return Err(BaseModError::InvalidDelta);
        }
        let end = tail.iter().position(|b| !b.is_ascii_digit()).unwrap_or(tail.len());
        if end == 0 {
            return Err(BaseModError::InvalidDelta);
        }
        #[allow(clippy::indexing_slicing, reason = "end > 0 and end <= tail.len()")]
        let digits = &tail[..end];
        let s = std::str::from_utf8(digits).map_err(|_| BaseModError::InvalidDelta)?;
        let n: u32 = s.parse().map_err(|_| BaseModError::InvalidDelta)?;
        out.push(n);
        rest = tail.get(end..).unwrap_or(&[]);
    }
    Ok(out)
}

#[allow(clippy::too_many_arguments, reason = "per-entry resolution is a local helper")]
fn resolve_and_emit(
    canonical: Base,
    strand: ModStrand,
    mod_types: &[ModType],
    mode: ModMode,
    deltas: &[u32],
    ml_slice: &[u8],
    seq: &[Base],
    is_reverse: bool,
    pending: &mut Vec<(u32, Modification)>,
    state: &mut BaseModState,
) -> Result<(), BaseModError> {
    // Track the Unmodified flag per canonical base for is_unmodified().
    if let Some(idx) = canonical.known_index()
        && mode == ModMode::Unmodified
    {
        #[allow(clippy::indexing_slicing, reason = "known_index() returns Some(0..=3)")]
        {
            state.any_unmodified_mode[idx] = true;
        }
    }

    if deltas.is_empty() {
        return Ok(());
    }

    // Target base to match while iterating stored sequence.
    let target = if is_reverse { canonical.inverse() } else { canonical };

    // Iterate stored indices in the correct order.
    //
    // Forward: i = 0, 1, 2, ..., seq.len()-1
    // Reverse: i = seq.len()-1, ..., 0  (so the i-th match in the original
    //          sequence corresponds to the stored index we hit).
    let seq_len = seq.len();
    let mut delta_idx: usize = 0;
    let mut remaining: u32 = *deltas.first().ok_or(BaseModError::InvalidDelta)?;
    let num_mods = mod_types.len();

    let mut stored_idx: usize = if is_reverse { seq_len.saturating_sub(1) } else { 0 };
    let mut steps: usize = 0;
    while steps < seq_len && delta_idx < deltas.len() {
        #[allow(
            clippy::indexing_slicing,
            reason = "stored_idx < seq_len by the while guard / decrement logic"
        )]
        let b = seq[stored_idx];
        if b == target {
            if remaining == 0 {
                #[allow(
                    clippy::cast_possible_truncation,
                    reason = "stored_idx < seq_len and seq.len() fits in u32 for practical BAM records"
                )]
                let qp = stored_idx as u32;
                // Push one Modification per mod_type (combined codes interleaved in ML).
                let base_offset = delta_idx
                    .checked_mul(num_mods)
                    .ok_or(BaseModError::MlLengthMismatch { ml: ml_slice.len(), mm: usize::MAX })?;
                for (k, mt) in mod_types.iter().enumerate() {
                    let ml_idx =
                        base_offset.checked_add(k).ok_or(BaseModError::MlLengthMismatch {
                            ml: ml_slice.len(),
                            mm: usize::MAX,
                        })?;
                    let prob = *ml_slice.get(ml_idx).ok_or(BaseModError::MlLengthMismatch {
                        ml: ml_slice.len(),
                        mm: usize::MAX,
                    })?;
                    pending.push((
                        qp,
                        Modification {
                            mod_type: *mt,
                            probability: prob,
                            canonical_base: canonical,
                            strand,
                        },
                    ));
                }
                // Signal match/mode for is_unmodified bookkeeping unaffected;
                // modifications themselves are recorded as mod calls.
                let _ = mode;

                delta_idx = delta_idx.saturating_add(1);
                if delta_idx >= deltas.len() {
                    break;
                }
                #[allow(clippy::indexing_slicing, reason = "delta_idx < deltas.len()")]
                {
                    remaining = deltas[delta_idx];
                }
            } else {
                remaining = remaining.saturating_sub(1);
            }
        }
        steps = steps.saturating_add(1);
        if is_reverse {
            if stored_idx == 0 {
                break;
            }
            stored_idx = stored_idx.saturating_sub(1);
        } else {
            stored_idx = stored_idx.saturating_add(1);
        }
    }

    if delta_idx < deltas.len() {
        #[allow(clippy::indexing_slicing, reason = "delta_idx < deltas.len() by the guard")]
        return Err(BaseModError::DeltaOutOfRange {
            delta: deltas[delta_idx],
            base: canonical,
            is_reverse,
        });
    }
    Ok(())
}

// TODO(pileup-integration): once `PileupColumn::modification_at` is added
// (see spec `r[base_mod.pileup_integration]`), the pileup engine needs a
// lazy per-active-record cache of `BaseModState`. The natural shape is a
// side-table parallel to `active` indexed by alignment_idx, to keep
// `ActiveRecord` lean (hot loop only touches `active_end_pos`).

#[cfg(test)]
#[allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    reason = "test code"
)]
mod tests {
    use super::*;
    use seqair_types::Base::{A, C, G, T};

    fn seq(bases: &[Base]) -> Vec<Base> {
        bases.to_vec()
    }

    // r[verify base_mod.parse_mm]
    // r[verify base_mod.state]
    // r[verify base_mod.modification_struct]
    // r[verify base_mod.query_qpos]
    #[test]
    fn parses_simple_cpg_entry() {
        // Sequence: ACGT ACGC GACG T  (positions of C: 1, 5, 7, 10)
        let s = seq(&[A, C, G, T, A, C, G, C, G, A, C, G, T]);
        let state = BaseModState::parse(b"C+m,0,2;", &[200, 180], &s, false).unwrap();
        assert_eq!(state.len(), 2);
        // First C is at qpos 1; after skipping 2 C's (5, 7), next is at 10.
        let m1 = state.mod_at_qpos(1).unwrap();
        assert_eq!(m1.len(), 1);
        assert_eq!(m1[0].probability, 200);
        assert!(matches!(m1[0].mod_type, ModType::Code(b'm')));
        assert_eq!(m1[0].canonical_base, C);
        assert!(matches!(m1[0].strand, ModStrand::Plus));

        let m2 = state.mod_at_qpos(10).unwrap();
        assert_eq!(m2[0].probability, 180);
    }

    // r[verify base_mod.resolve_positions]
    #[test]
    fn delta_zero_is_first_occurrence() {
        let s = seq(&[A, C, C, C]);
        let state = BaseModState::parse(b"C+m,0;", &[128], &s, false).unwrap();
        assert_eq!(state.mod_at_qpos(1).unwrap()[0].probability, 128);
        assert!(state.mod_at_qpos(0).is_none());
        assert!(state.mod_at_qpos(2).is_none());
    }

    // r[verify base_mod.reverse_complement]
    #[test]
    fn reverse_strand_uses_complement_counting_from_end() {
        // Original 5'→3': ACGT ACG (C at original pos 1 and 5)
        // Stored (reverse-complemented): CGT ACGT (G at stored 0 and 4)
        // Original: A(0) C(1) G(2) T(3) A(4) C(5) G(6)
        // Reverse complement: complement each base, then reverse.
        // complement: T G C A T G C
        // reverse:    C G T A C G T
        // So stored is: C G T A C G T (length 7)
        // G occurrences in stored (from end): pos 5, pos 1
        // The original 0th C (original pos 1) should map to stored pos 5 (the
        // last-from-start, first-from-end G). The original 1st C (original
        // pos 5) should map to stored pos 1.
        let stored = seq(&[C, G, T, A, C, G, T]);
        let state = BaseModState::parse(b"C+m,0,0;", &[200, 210], &stored, true).unwrap();
        let m_first = state.mod_at_qpos(5).unwrap();
        assert_eq!(m_first[0].probability, 200);
        let m_second = state.mod_at_qpos(1).unwrap();
        assert_eq!(m_second[0].probability, 210);
    }

    // r[verify base_mod.parse_mm]
    #[test]
    fn parses_multiple_entries() {
        let s = seq(&[A, C, G, C, G, A, C, G]);
        let state = BaseModState::parse(b"C+m,0;C+h,1;", &[200, 150], &s, false).unwrap();
        // m on the 0th C (qpos 1); h on the 1st C (qpos 3) after skipping 1.
        assert!(matches!(state.mod_at_qpos(1).unwrap()[0].mod_type, ModType::Code(b'm')));
        assert!(matches!(state.mod_at_qpos(3).unwrap()[0].mod_type, ModType::Code(b'h')));
    }

    // r[verify base_mod.parse_mm]
    #[test]
    fn combined_codes_interleave_ml() {
        // `C+mh,0,2;` with ML = [m0, h0, m1, h1]
        // 4 C's at qpos 1, 3, 5, 7. delta 0 → C@1; delta 2 → skip 2 more (3,5) → C@7.
        let s = seq(&[A, C, A, C, A, C, A, C]);
        let state = BaseModState::parse(b"C+mh,0,2;", &[200, 100, 220, 120], &s, false).unwrap();
        let at1 = state.mod_at_qpos(1).unwrap();
        assert_eq!(at1.len(), 2);
        assert!(matches!(at1[0].mod_type, ModType::Code(b'm')));
        assert_eq!(at1[0].probability, 200);
        assert!(matches!(at1[1].mod_type, ModType::Code(b'h')));
        assert_eq!(at1[1].probability, 100);
        let at7 = state.mod_at_qpos(7).unwrap();
        assert_eq!(at7.len(), 2);
        assert_eq!(at7[0].probability, 220);
        assert_eq!(at7[1].probability, 120);
    }

    // r[verify base_mod.parse_mm]
    #[test]
    fn parses_chebi_code() {
        let s = seq(&[A, C, G, T]);
        let state = BaseModState::parse(b"C+27551,0;", &[200], &s, false).unwrap();
        assert!(matches!(state.mod_at_qpos(1).unwrap()[0].mod_type, ModType::ChEBI(27551)));
    }

    // r[verify base_mod.implicit_explicit]
    #[test]
    fn is_unmodified_three_modes() {
        let s = seq(&[A, C, G, C, G, C]);

        // Implicit mode: unlisted positions are unknown.
        let implicit = BaseModState::parse(b"C+m,0;", &[200], &s, false).unwrap();
        assert_eq!(implicit.is_unmodified(1, C), Some(false));
        assert_eq!(implicit.is_unmodified(3, C), None);

        // Explicit `.`: unlisted positions are definitively unmodified.
        let explicit = BaseModState::parse(b"C+m.,0;", &[200], &s, false).unwrap();
        assert_eq!(explicit.is_unmodified(1, C), Some(false));
        assert_eq!(explicit.is_unmodified(3, C), Some(true));

        // Ambiguous `?`: same observable behavior as Implicit.
        let ambiguous = BaseModState::parse(b"C+m?,0;", &[200], &s, false).unwrap();
        assert_eq!(ambiguous.is_unmodified(3, C), None);
    }

    // r[verify base_mod.query_refpos]
    #[test]
    fn query_refpos_via_cigar() {
        use super::super::cigar::CigarMapping;
        use seqair_types::Pos0;

        // Read starts at ref 1000, 10M. C's at qpos 1, 3 → ref 1001, 1003.
        let s = seq(&[A, C, G, C, G, A, C, G, T, A]);
        let state = BaseModState::parse(b"C+m,0,0;", &[200, 180], &s, false).unwrap();

        // Build a simple 10M CIGAR via bytes (BAM-packed: len<<4 | op, M=0).
        let cigar_op: u32 = 10u32 << 4;
        let cigar_bytes = cigar_op.to_le_bytes().to_vec();
        let cm = CigarMapping::new(Pos0::new(1000).unwrap(), &cigar_bytes).unwrap();

        let m = state.mod_at_ref_pos(Pos0::new(1001).unwrap(), &cm).unwrap();
        assert_eq!(m[0].probability, 200);
        let m2 = state.mod_at_ref_pos(Pos0::new(1003).unwrap(), &cm).unwrap();
        assert_eq!(m2[0].probability, 180);
        assert!(state.mod_at_ref_pos(Pos0::new(1002).unwrap(), &cm).is_none());
    }

    // r[verify base_mod.validation]
    #[test]
    fn validation_ml_length_mismatch() {
        let s = seq(&[A, C, G]);
        let err = BaseModState::parse(b"C+m,0;", &[200, 180], &s, false).unwrap_err();
        assert!(matches!(err, BaseModError::MlLengthMismatch { .. }), "got {err:?}");
    }

    // r[verify base_mod.validation]
    #[test]
    fn validation_missing_terminator() {
        let s = seq(&[A, C, G]);
        let err = BaseModState::parse(b"C+m,0", &[200], &s, false).unwrap_err();
        assert!(matches!(err, BaseModError::MissingTerminator), "got {err:?}");
    }

    // r[verify base_mod.validation]
    #[test]
    fn validation_delta_out_of_range() {
        let s = seq(&[A, C]); // only 1 C
        let err = BaseModState::parse(b"C+m,5;", &[200], &s, false).unwrap_err();
        assert!(matches!(err, BaseModError::DeltaOutOfRange { .. }), "got {err:?}");
    }

    // r[verify base_mod.validation]
    #[test]
    fn validation_bad_canonical() {
        let s = seq(&[A, C, G]);
        let err = BaseModState::parse(b"Z+m,0;", &[200], &s, false).unwrap_err();
        assert!(matches!(err, BaseModError::InvalidCanonicalBase(_)), "got {err:?}");
    }

    // r[verify base_mod.validation]
    #[test]
    fn validation_bad_strand() {
        let s = seq(&[A, C, G]);
        let err = BaseModState::parse(b"C*m,0;", &[200], &s, false).unwrap_err();
        assert!(matches!(err, BaseModError::InvalidStrand(_)), "got {err:?}");
    }

    #[test]
    fn empty_mm_empty_ml_ok() {
        let s = seq(&[A, C, G]);
        let state = BaseModState::parse(b"", &[], &s, false).unwrap();
        assert!(state.is_empty());
        assert!(state.mod_at_qpos(1).is_none());
        assert_eq!(state.is_unmodified(1, C), None);
    }

    #[test]
    fn trailing_nul_tolerated() {
        let s = seq(&[A, C, G]);
        let state = BaseModState::parse(b"C+m,0;\0", &[200], &s, false).unwrap();
        assert_eq!(state.len(), 1);
    }

    #[test]
    fn strand_minus_parses() {
        let s = seq(&[A, C, G]);
        let state = BaseModState::parse(b"C-m,0;", &[200], &s, false).unwrap();
        assert!(matches!(state.mod_at_qpos(1).unwrap()[0].strand, ModStrand::Minus));
    }

    #[test]
    fn lowercase_canonical_base() {
        let s = seq(&[A, C, G]);
        let state = BaseModState::parse(b"c+m,0;", &[200], &s, false).unwrap();
        assert_eq!(state.mod_at_qpos(1).unwrap()[0].canonical_base, C);
    }

    #[test]
    fn two_entries_same_canonical_sort_by_qpos() {
        // Second entry's first hit is earlier than first entry's first hit.
        let s = seq(&[A, C, G, C, G, C]);
        // C+m,2 → skip 2 C's → C@5. C+h,0 → C@1.
        let state = BaseModState::parse(b"C+m,2;C+h,0;", &[200, 100], &s, false).unwrap();
        assert!(matches!(state.mod_at_qpos(1).unwrap()[0].mod_type, ModType::Code(b'h')));
        assert!(matches!(state.mod_at_qpos(5).unwrap()[0].mod_type, ModType::Code(b'm')));
    }
}
