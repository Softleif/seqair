//! Walk a record's CIGAR operations, yielding typed [`AlignedPair`] variants.
//!
//! This is the record-walking complement to [`CigarMapping::pos_info_at`](super::cigar::CigarMapping::pos_info_at),
//! which looks up query positions by reference position (column-at-a-time for the pileup engine).
//! `AlignedPairs` iterates per-CIGAR-op — useful for recomputing NM/MD tags, per-position
//! annotations, motif inference, and any per-record walk that needs both `qpos` and `rpos`.
//!
//! Default mode yields only position-bearing ops (`Match`, `Insertion`, `Deletion`, `RefSkip`).
//! Opt in to `SoftClip` via [`.with_soft_clips()`](AlignedPairs::with_soft_clips) or
//! all ops via [`.full()`](AlignedPairs::full).
//!
//! # Layered views
//!
//! [`AlignedPairs::with_read`] attaches the read's seq/qual to the events;
//! `.with_reference(&ref_seq)` further attaches reference base lookups.
//! See [`super::aligned_pairs_view`] for the layered API and the
//! `aligned_pairs_walk` example crate-side for an end-to-end walk over a real
//! BAM/FASTA pair (`cargo run --example aligned_pairs_walk -- --help`).
//!
//! # Quick reference
//!
//! | Need | API |
//! |---|---|
//! | Positions only (cheapest) | `slim.aligned_pairs(store)?` |
//! | Positions + read base/qual | `slim.aligned_pairs_with_read(store)?` |
//! | Positions + read + ref base | `…aligned_pairs_with_read(store)?.with_reference(&ref_seq)` |
//! | M-vs-=-vs-X distinction | `AlignedPair::Match { kind, .. }` (`MatchKind` enum) |
//! | Soft clips visible | `.with_soft_clips()` on any layer |
//! | Padding/Unknown visible | `.full()` on any layer |

use super::cigar::{CigarOp, CigarOpType};
use super::record_store::{RecordAccessError, RecordStore, SlimRecord};
use seqair_types::Pos0;
use thiserror::Error;

// ── Errors ─────────────────────────────────────────────────────────────────

/// Errors raised by `aligned_pairs*` constructors.
///
/// Iteration itself is infallible — once the iterator is constructed it does
/// not return errors. All validation happens up front.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum AlignedPairsError {
    /// Failed to read CIGAR / seq / qual bytes from the [`RecordStore`].
    #[error("record access error")]
    Access {
        #[from]
        source: RecordAccessError,
    },

    /// An unmapped record (`pos = None`) carries CIGAR operations.
    /// Per [SAM1] §1.4, unmapped reads must have an empty CIGAR. Walking
    /// CIGAR ops without a base reference position would produce nonsense
    /// `rpos` values, so we refuse rather than silently anchor to position 0.
    #[error(
        "unmapped record (pos = None) has {cigar_ops} CIGAR op(s); only an empty CIGAR is valid \
         for unmapped records — check flags and CIGAR consistency before constructing the walk"
    )]
    UnmappedWithCigar { cigar_ops: usize },

    /// CIGAR's query-consuming length disagrees with `seq.len()`. The walk
    /// would slice past the end of the sequence and silently substitute
    /// [`Base::Unknown`](seqair_types::Base::Unknown), which is
    /// indistinguishable from a legitimate `N` in the read.
    #[error(
        "CIGAR query-consuming length {cigar_qlen} != seq length {seq_len} — record's seq slab \
         and CIGAR are out of sync"
    )]
    CigarSeqLengthMismatch { cigar_qlen: u64, seq_len: usize },

    /// Quality slab length disagrees with `seq.len()`. BAM allows missing
    /// qual (filled with `0xFF`), in which case `qual.len() == 0` is
    /// accepted; any other length is a structural mismatch.
    #[error(
        "qual length {qual_len} != seq length {seq_len} — qual slab must be either empty (BAM \
         missing-qual sentinel) or match seq length exactly"
    )]
    SeqQualLengthMismatch { seq_len: usize, qual_len: usize },
}

// ── MatchKind ──────────────────────────────────────────────────────────────

// r[impl cigar.aligned_pairs.match_kind]
/// Distinguishes the three aligner-emitted "matched-base" CIGAR ops.
///
/// htslib's `aligned_pairs_full` collapses M/=/X into a single tuple shape,
/// losing the per-op distinction. [`AlignedPair::Match`] preserves it via this
/// field so callers that care (e.g. NM/MD recompute, variant calling against
/// `=`/`X` CIGARs) don't have to re-walk the CIGAR.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MatchKind {
    /// `M` (BAM op 0) — alignment match, ambiguous between sequence match and
    /// mismatch. The aligner did not commit either way.
    Match,
    /// `=` (BAM op 7) — explicit sequence match.
    SeqMatch,
    /// `X` (BAM op 8) — explicit sequence mismatch.
    SeqMismatch,
}

// ── AlignedPair ────────────────────────────────────────────────────────────

// r[impl cigar.aligned_pairs.types]
/// A typed aligned pair — one event from walking a record's CIGAR.
///
/// Unlike htslib's `(Option<i64>, Option<i64>)`, this enum uses typed variants:
/// `Match` carries both positions plus a [`MatchKind`] tag, `Insertion`
/// carries `qpos`+len, `Deletion` carries `rpos`+len, `RefSkip` carries
/// `rpos`+len. Match ops expand per-position; I/D/N/S ops yield once per op
/// (summary form).
///
/// # Example
/// ```ignore
/// for pair in record.aligned_pairs(store)? {
///     match pair {
///         AlignedPair::Match { qpos, rpos, kind } => { /* M, =, or X */ }
///         AlignedPair::Insertion { qpos, insert_len } => { /* I */ }
///         AlignedPair::Deletion { rpos, del_len } => { /* D */ }
///         AlignedPair::RefSkip { rpos, skip_len } => { /* N */ }
///         _ => {}
///     }
/// }
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignedPair {
    /// M / = / X op — read base aligned to a reference base.
    /// Yields one variant per consumed position. Use `kind` to distinguish
    /// `M` (ambiguous) from `=` (explicit match) and `X` (explicit mismatch).
    Match { qpos: u32, rpos: Pos0, kind: MatchKind },

    /// I op — `insert_len` query bases follow `qpos`, no reference span.
    /// `qpos` is the position of the first inserted base.
    Insertion { qpos: u32, insert_len: u32 },

    /// D op — `del_len` reference bases starting at `rpos`, no query span.
    Deletion { rpos: Pos0, del_len: u32 },

    /// N op — reference skip (e.g. intron).
    RefSkip { rpos: Pos0, skip_len: u32 },

    /// S op — soft clip. `qpos` is the start of the clipped run.
    /// Hidden by default; opt in via [`.with_soft_clips()`](AlignedPairs::with_soft_clips) or [`.full()`](AlignedPairs::full).
    SoftClip { qpos: u32, len: u32 },

    /// P op — padding. `len` is the op length (typically small — most aligners
    /// emit `1P` or `0P`). Hidden by default; opt in via [`.full()`](AlignedPairs::full).
    Padding { len: u32 },

    /// Reserved op code (9..=15) we tolerate on read. `code` is the raw 4-bit op.
    /// Yielded by [`.full()`](AlignedPairs::full) only.
    Unknown { code: u8, len: u32 },
}

// Compile-time size guard: enum stays small enough for a hot-path iterator.
// Match { u32, u32, MatchKind(u8) } = 9 bytes payload, padded to 12. Plus a
// 1-byte discriminant rounded up by alignment → 16 bytes total.
const _: () = assert!(
    std::mem::size_of::<AlignedPair>() <= 16,
    "AlignedPair grew past 16 bytes — revisit before accepting the hit"
);

// ── AlignedPairsOptions ────────────────────────────────────────────────────

// r[impl cigar.aligned_pairs.options]
/// Controls which [`AlignedPair`] variants the iterator yields.
#[derive(Debug, Clone, Copy, Default)]
pub struct AlignedPairsOptions {
    /// Yield [`AlignedPair::SoftClip`] variants. Default `false`.
    pub soft_clips: bool,
    /// Yield [`AlignedPair::Padding`] and [`AlignedPair::Unknown`].
    /// Default `false` — most consumers are confused by P ops.
    pub padding_and_unknown: bool,
}

// ── AlignedPairs iterator ─────────────────────────────────────────────────

// r[impl cigar.aligned_pairs.iterator]
// r[impl cigar.aligned_pairs.default_mode]
/// Iterator over a record's CIGAR operations, yielding typed [`AlignedPair`]s.
///
/// # Default mode
///
/// Only `Match`, `Insertion`, `Deletion`, and `RefSkip` are yielded.
/// Soft clips, padding, and unknown ops are silently skipped.
///
/// # Builder methods
///
/// - `.with_soft_clips()` — also yield `SoftClip` variants
/// - `.full()` — yield everything including `Padding` and `Unknown`
///
/// # Example
///
/// ```ignore
/// for pair in record.aligned_pairs(store)?.full() {
///     match pair {
///         AlignedPair::Match { qpos, rpos } => { /* ... */ }
///         AlignedPair::Insertion { qpos, insert_len } => { /* ... */ }
///         // ...
///     }
/// }
/// ```
#[derive(Debug)]
pub struct AlignedPairs<'a> {
    /// Remaining CIGAR ops to process.
    ops: &'a [CigarOp],
    /// Current 0-based query position (includes soft clips).
    qpos: u32,
    /// Current 0-based reference position.
    rpos: Pos0,
    /// Which variants to yield.
    options: AlignedPairsOptions,
    /// Multi-base expansion state for M/=/X ops.
    expanding: ExpandingState,
}

#[derive(Debug, Clone, Copy)]
enum ExpandingState {
    None,
    /// Yield one Match per remaining base. `kind` is fixed for the run
    /// (all expanded yields of an M op are M, not interleaved with =/X).
    PerBase {
        remaining: u32,
        kind: MatchKind,
    },
}

impl<'a> AlignedPairs<'a> {
    // r[impl cigar.aligned_pairs.default_mode]
    /// Create a new iterator over `cigar` starting at reference position `rec_pos`.
    ///
    /// Default mode yields `Match`, `Insertion`, `Deletion`, `RefSkip` only.
    pub fn new(rec_pos: Pos0, cigar: &'a [CigarOp]) -> Self {
        Self {
            ops: cigar,
            qpos: 0,
            rpos: rec_pos,
            options: AlignedPairsOptions::default(),
            expanding: ExpandingState::None,
        }
    }

    // r[impl cigar.aligned_pairs.options]
    /// Also yield [`AlignedPair::SoftClip`] variants.
    pub fn with_soft_clips(mut self) -> Self {
        self.options.soft_clips = true;
        self
    }

    // r[impl cigar.aligned_pairs.options]
    /// Yield all variants including [`AlignedPair::SoftClip`],
    /// [`AlignedPair::Padding`], and [`AlignedPair::Unknown`].
    pub fn full(mut self) -> Self {
        self.options.soft_clips = true;
        self.options.padding_and_unknown = true;
        self
    }

    // r[impl cigar.aligned_pairs.with_read.iterator]
    /// Attach the read's sequence and quality slabs, producing an
    /// [`AlignedPairsWithRead`](super::aligned_pairs_view::AlignedPairsWithRead)
    /// that yields events enriched with the per-position read base and qual,
    /// and pre-sliced inserted/clipped runs.
    ///
    /// Validates up front:
    /// - the iterator's remaining query-consuming CIGAR length fits inside
    ///   `seq[qpos..]` — otherwise returns
    ///   [`AlignedPairsError::CigarSeqLengthMismatch`];
    /// - `qual.len()` equals `seq.len()` or is zero (BAM missing-qual
    ///   sentinel) — otherwise returns
    ///   [`AlignedPairsError::SeqQualLengthMismatch`].
    ///
    /// Once construction succeeds, iteration is infallible — every yielded
    /// `Match`/`Insertion`/`SoftClip` index is provably in-bounds.
    ///
    /// Lifetimes are split: `'cigar` ties to the CIGAR slice, `'read` ties to
    /// the seq/qual slabs. They may differ, allowing callers to compose
    /// borrows from different stores or buffers.
    ///
    /// # Example
    /// ```ignore
    /// for ev in slim.aligned_pairs(&store)?.with_read(seq, qual)? {
    ///     // ...
    /// }
    /// ```
    pub fn with_read<'read>(
        self,
        seq: &'read [seqair_types::Base],
        qual: &'read [seqair_types::BaseQuality],
    ) -> Result<super::aligned_pairs_view::AlignedPairsWithRead<'a, 'read>, AlignedPairsError> {
        // Compute the remaining query-consuming CIGAR length and compare to
        // `seq[qpos..]`. We work in u64 so a runaway 32-bit CIGAR can't
        // wrap; real-world records are far below this.
        let cigar_remaining_qlen: u64 =
            self.ops.iter().filter(|op| op.consumes_query()).map(|op| u64::from(op.len())).sum();
        let qpos_consumed = u64::from(self.qpos);
        let total_qlen = cigar_remaining_qlen.saturating_add(qpos_consumed);
        // total_qlen represents the total query-consuming length the iterator
        // expects across the *full* CIGAR (already-consumed + remaining). The
        // seq slab must cover all of it. Strict equality lets us also catch
        // the inverse mismatch (seq longer than CIGAR claims).
        if total_qlen != seq.len() as u64 {
            return Err(AlignedPairsError::CigarSeqLengthMismatch {
                cigar_qlen: total_qlen,
                seq_len: seq.len(),
            });
        }
        if !qual.is_empty() && qual.len() != seq.len() {
            return Err(AlignedPairsError::SeqQualLengthMismatch {
                seq_len: seq.len(),
                qual_len: qual.len(),
            });
        }
        Ok(super::aligned_pairs_view::AlignedPairsWithRead::new(self, seq, qual))
    }
}

impl Iterator for AlignedPairs<'_> {
    type Item = AlignedPair;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // ── Expand multi-base ops (M/=/X) ──
            match std::mem::replace(&mut self.expanding, ExpandingState::None) {
                ExpandingState::PerBase { remaining, kind } if remaining > 0 => {
                    self.expanding = ExpandingState::PerBase { remaining: remaining - 1, kind };
                    let qpos = self.qpos;
                    let rpos = self.rpos;
                    self.qpos = self.qpos.saturating_add(1);
                    self.rpos = advance_rpos(self.rpos, 1);
                    return Some(AlignedPair::Match { qpos, rpos, kind });
                }
                ExpandingState::PerBase { .. } => {
                    // Remaining was 0 — already replaced with None, loop to next op.
                    continue;
                }
                ExpandingState::None => {}
            }

            // ── Next CIGAR op ──
            let (op, rest) = self.ops.split_first()?;
            self.ops = rest;

            let len = op.len();
            let op_type = op.op_type();
            // r[impl cigar.aligned_pairs.zero_length_ops]
            // Zero-length ops are degenerate — htslib's CIGAR walker skips them
            // (its for-loops simply don't iterate). We mirror that for every op
            // kind so callers never see a `len: 0` summary.
            // HardClip/Padding/Unknown still need to be considered for advancing
            // qpos/rpos but those don't move when len == 0 anyway.
            if len == 0 && !matches!(op_type, CigarOpType::HardClip) {
                continue;
            }
            match op_type {
                CigarOpType::Match => {
                    self.expanding =
                        ExpandingState::PerBase { remaining: len, kind: MatchKind::Match };
                    continue;
                }
                CigarOpType::SeqMatch => {
                    self.expanding =
                        ExpandingState::PerBase { remaining: len, kind: MatchKind::SeqMatch };
                    continue;
                }
                CigarOpType::SeqMismatch => {
                    self.expanding =
                        ExpandingState::PerBase { remaining: len, kind: MatchKind::SeqMismatch };
                    continue;
                }

                // r[impl cigar.aligned_pairs.insertion_qpos]
                CigarOpType::Insertion => {
                    let qpos = self.qpos;
                    self.qpos = self.qpos.saturating_add(len);
                    return Some(AlignedPair::Insertion { qpos, insert_len: len });
                }

                // r[impl cigar.aligned_pairs.deletion_rpos]
                CigarOpType::Deletion => {
                    let rpos = self.rpos;
                    self.rpos = advance_rpos(self.rpos, len);
                    return Some(AlignedPair::Deletion { rpos, del_len: len });
                }

                CigarOpType::RefSkip => {
                    let rpos = self.rpos;
                    self.rpos = advance_rpos(self.rpos, len);
                    return Some(AlignedPair::RefSkip { rpos, skip_len: len });
                }

                CigarOpType::SoftClip => {
                    if self.options.soft_clips {
                        let qpos = self.qpos;
                        self.qpos = self.qpos.saturating_add(len);
                        return Some(AlignedPair::SoftClip { qpos, len });
                    }
                    // r[impl cigar.aligned_pairs.qpos_semantics]
                    // Even when hidden, soft clips still advance qpos.
                    self.qpos = self.qpos.saturating_add(len);
                }

                // r[impl cigar.aligned_pairs.hard_clips]
                CigarOpType::HardClip => {
                    // Hard clips are never yielded and don't advance anything.
                }

                CigarOpType::Padding => {
                    if self.options.padding_and_unknown {
                        return Some(AlignedPair::Padding { len });
                    }
                }

                CigarOpType::Unknown(code) => {
                    if self.options.padding_and_unknown {
                        return Some(AlignedPair::Unknown { code, len });
                    }
                }
            }
        }
    }
}

/// Advance a `Pos0` by `n` positions, saturating at `i32::MAX`.
/// Used internally by the iterator; overflow only happens on degenerate BAM
/// (positions would already exceed BAM's i32 range) and saturation is the
/// safe fallback.
#[inline]
fn advance_rpos(rpos: Pos0, n: u32) -> Pos0 {
    let val = rpos.as_u64().saturating_add(u64::from(n));
    let clamped = val.min(i32::MAX as u64) as u32;
    // clamped ≤ i32::MAX, so Pos0::new can't fail
    Pos0::new(clamped).unwrap_or(rpos)
}

// ── SlimRecord integration ────────────────────────────────────────────────

// r[impl cigar.aligned_pairs.slim_record]
impl SlimRecord {
    /// Walk this record's CIGAR, yielding typed [`AlignedPair`]s.
    ///
    /// Reads the CIGAR slice from the store at construction time.
    /// Once the slice is obtained, iteration is infallible.
    ///
    /// Default mode yields `Match`, `Insertion`, `Deletion`, `RefSkip`.
    /// Chain `.with_soft_clips()` or `.full()` on the result to widen.
    ///
    /// # Example
    /// ```ignore
    /// for pair in slim_record.aligned_pairs(&store)? {
    ///     if let AlignedPair::Match { qpos, rpos, kind } = pair {
    ///         // ...
    ///     }
    /// }
    /// ```
    pub fn aligned_pairs<'store, U>(
        &self,
        store: &'store RecordStore<U>,
    ) -> Result<AlignedPairs<'store>, RecordAccessError> {
        Ok(AlignedPairs::new(self.pos, self.cigar(store)?))
    }

    // r[impl cigar.aligned_pairs.with_read.slim_record]
    /// Walk this record's CIGAR with read seq/qual attached. One-shot helper:
    /// equivalent to `self.aligned_pairs(store)?.with_read(self.seq(store)?, self.qual(store)?)`
    /// but only borrows `store` once.
    ///
    /// Chain `.with_reference(&ref_seq)` to add reference-base lookup.
    ///
    /// # Example
    /// ```ignore
    /// for ev in slim.aligned_pairs_with_read(&store)? {
    ///     if let AlignedPairWithRead::Match { qpos, query, qual, .. } = ev {
    ///         // ...
    ///     }
    /// }
    /// ```
    pub fn aligned_pairs_with_read<'store, U>(
        &self,
        store: &'store RecordStore<U>,
    ) -> Result<super::aligned_pairs_view::AlignedPairsWithRead<'store, 'store>, AlignedPairsError>
    {
        let cigar = self.cigar(store)?;
        let seq = self.seq(store)?;
        let qual = self.qual(store)?;
        // Records pushed via `push_raw` always have cigar.qlen == seq.len() ==
        // qual.len(), so this validation should always succeed for store-resident
        // records — but the typed error makes the "should" provable.
        AlignedPairs::new(self.pos, cigar).with_read(seq, qual)
    }
}

// ── OwnedBamRecord integration ────────────────────────────────────────────

// r[impl cigar.aligned_pairs.owned_record]
impl super::owned_record::OwnedBamRecord {
    /// Walk this record's CIGAR, yielding typed [`AlignedPair`]s.
    ///
    /// Returns [`AlignedPairsError::UnmappedWithCigar`] if `pos = None` but
    /// the CIGAR is non-empty — walking that CIGAR without a base reference
    /// position would produce nonsense `rpos` values. For mapped records
    /// (`pos = Some(_)`) and fully-unmapped records (`pos = None` + empty
    /// CIGAR), iteration over the returned `AlignedPairs` is infallible.
    ///
    /// # Example
    /// ```ignore
    /// for pair in record.aligned_pairs()? {
    ///     // ...
    /// }
    /// ```
    pub fn aligned_pairs(&self) -> Result<AlignedPairs<'_>, AlignedPairsError> {
        match self.pos {
            Some(p) => Ok(AlignedPairs::new(p, &self.cigar)),
            None if self.cigar.is_empty() => Ok(AlignedPairs::new(Pos0::ZERO, &self.cigar)),
            None => Err(AlignedPairsError::UnmappedWithCigar { cigar_ops: self.cigar.len() }),
        }
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────

#[cfg(test)]
#[allow(clippy::arithmetic_side_effects, reason = "test arithmetic on known small values")]
mod tests {
    use super::super::cigar::CigarOp as Op;
    use super::*;

    // Helper: construct a CigarOp
    fn m(len: u32) -> Op {
        Op::new(CigarOpType::Match, len)
    }
    fn ins(len: u32) -> Op {
        Op::new(CigarOpType::Insertion, len)
    }
    fn del(len: u32) -> Op {
        Op::new(CigarOpType::Deletion, len)
    }
    fn skip(len: u32) -> Op {
        Op::new(CigarOpType::RefSkip, len)
    }
    fn soft(len: u32) -> Op {
        Op::new(CigarOpType::SoftClip, len)
    }
    fn hard(len: u32) -> Op {
        Op::new(CigarOpType::HardClip, len)
    }
    fn pad(len: u32) -> Op {
        Op::new(CigarOpType::Padding, len)
    }
    fn unk(code: u8, len: u32) -> Op {
        Op::new(CigarOpType::Unknown(code), len)
    }

    fn p0(v: u32) -> Pos0 {
        Pos0::new(v).unwrap()
    }

    fn match_m(qpos: u32, rpos: Pos0) -> AlignedPair {
        AlignedPair::Match { qpos, rpos, kind: MatchKind::Match }
    }

    fn match_eq(qpos: u32, rpos: Pos0) -> AlignedPair {
        AlignedPair::Match { qpos, rpos, kind: MatchKind::SeqMatch }
    }

    fn match_x(qpos: u32, rpos: Pos0) -> AlignedPair {
        AlignedPair::Match { qpos, rpos, kind: MatchKind::SeqMismatch }
    }

    fn seq_match(len: u32) -> Op {
        Op::new(CigarOpType::SeqMatch, len)
    }

    fn seq_mismatch(len: u32) -> Op {
        Op::new(CigarOpType::SeqMismatch, len)
    }

    // ── Basic unit tests ──────────────────────────────────────────────────

    // r[verify cigar.aligned_pairs.default_mode]
    #[test]
    fn simple_match() {
        let cigar = [m(5)];
        let pairs: Vec<_> = AlignedPairs::new(p0(100), &cigar).collect();
        assert_eq!(pairs.len(), 5);
        assert_eq!(pairs[0], match_m(0, p0(100)));
        assert_eq!(pairs[4], match_m(4, p0(104)));
    }

    // r[verify cigar.aligned_pairs.match_kind]
    #[test]
    fn match_kind_is_preserved_per_op() {
        // 2M + 2= + 2X — each op type yields its own MatchKind
        let cigar = [m(2), seq_match(2), seq_mismatch(2)];
        let pairs: Vec<_> = AlignedPairs::new(p0(0), &cigar).collect();
        assert_eq!(pairs.len(), 6);
        assert_eq!(pairs[0], match_m(0, p0(0)));
        assert_eq!(pairs[1], match_m(1, p0(1)));
        assert_eq!(pairs[2], match_eq(2, p0(2)));
        assert_eq!(pairs[3], match_eq(3, p0(3)));
        assert_eq!(pairs[4], match_x(4, p0(4)));
        assert_eq!(pairs[5], match_x(5, p0(5)));
    }

    // r[verify cigar.aligned_pairs.match_kind]
    #[test]
    fn match_kind_does_not_change_qpos_rpos() {
        // Walk a single CIGAR mixing all three kinds. The (qpos, rpos)
        // sequence under kind variation must equal the same sequence with all
        // ops collapsed to plain `M`. If the iterator's position cursor ever
        // drifted based on `kind`, this would fire.
        //
        // Earlier version of this test compared three separate CIGARs (5M,
        // 5=, 5X) and would pass even if MatchKind handling were completely
        // broken — included this as a regression note.
        let mixed = [m(2), seq_match(3), seq_mismatch(2), m(1)];
        let collapsed = [m(8)];

        let strip_kind = |pairs: Vec<AlignedPair>| -> Vec<(u32, u32)> {
            pairs
                .iter()
                .filter_map(|p| match p {
                    AlignedPair::Match { qpos, rpos, .. } => Some((*qpos, **rpos)),
                    _ => None,
                })
                .collect()
        };
        let mixed_positions = strip_kind(AlignedPairs::new(p0(50), &mixed).collect());
        let collapsed_positions = strip_kind(AlignedPairs::new(p0(50), &collapsed).collect());
        assert_eq!(mixed_positions, collapsed_positions, "kind must not affect qpos/rpos");
        // Sanity check: the mixed walk really did expose all three kinds.
        let kinds: Vec<MatchKind> = AlignedPairs::new(p0(50), &mixed)
            .filter_map(|p| match p {
                AlignedPair::Match { kind, .. } => Some(kind),
                _ => None,
            })
            .collect();
        assert!(kinds.contains(&MatchKind::Match), "mixed walk must have M");
        assert!(kinds.contains(&MatchKind::SeqMatch), "mixed walk must have =");
        assert!(kinds.contains(&MatchKind::SeqMismatch), "mixed walk must have X");
    }

    // r[verify cigar.aligned_pairs.insertion_qpos]
    #[test]
    fn with_insertion() {
        let cigar = [m(2), ins(1), m(2)];
        let pairs: Vec<_> = AlignedPairs::new(p0(100), &cigar).collect();
        assert_eq!(pairs.len(), 5, "2M + 1I(summary) + 2M = 5 items");
        assert_eq!(pairs[0], match_m(0, p0(100)));
        assert_eq!(pairs[1], match_m(1, p0(101)));
        // Insertion summary: qpos = 2 (first inserted base), len = 1
        assert_eq!(pairs[2], AlignedPair::Insertion { qpos: 2, insert_len: 1 });
        assert_eq!(pairs[3], match_m(3, p0(102)));
        assert_eq!(pairs[4], match_m(4, p0(103)));
    }

    // r[verify cigar.aligned_pairs.deletion_rpos]
    #[test]
    fn with_deletion() {
        let cigar = [m(2), del(3), m(2)];
        let pairs: Vec<_> = AlignedPairs::new(p0(100), &cigar).collect();
        assert_eq!(pairs.len(), 5, "2M + 1D(summary) + 2M = 5 items");
        assert_eq!(pairs[0], match_m(0, p0(100)));
        assert_eq!(pairs[1], match_m(1, p0(101)));
        // Deletion summary: rpos = 102 (start of deletion), del_len = 3
        assert_eq!(pairs[2], AlignedPair::Deletion { rpos: p0(102), del_len: 3 });
        assert_eq!(pairs[3], match_m(2, p0(105)));
        assert_eq!(pairs[4], match_m(3, p0(106)));
    }

    #[test]
    fn with_refskip() {
        let cigar = [m(2), skip(5), m(3)];
        let pairs: Vec<_> = AlignedPairs::new(p0(100), &cigar).collect();
        assert_eq!(pairs.len(), 6, "2M + 1N(summary) + 3M = 6 items");
        assert_eq!(pairs[2], AlignedPair::RefSkip { rpos: p0(102), skip_len: 5 });
        assert_eq!(pairs[3], match_m(2, p0(107)));
    }

    // r[verify cigar.aligned_pairs.qpos_semantics]
    #[test]
    fn qpos_includes_soft_clips() {
        // 5S + 3M: first match should have qpos = 5
        let cigar = [soft(5), m(3)];
        let pairs: Vec<_> = AlignedPairs::new(p0(0), &cigar).collect();
        assert_eq!(pairs.len(), 3);
        assert_eq!(pairs[0], match_m(5, p0(0)));
        assert_eq!(pairs[2], match_m(7, p0(2)));
    }

    // r[verify cigar.aligned_pairs.hard_clips]
    #[test]
    fn hard_clips_ignored() {
        let cigar = [hard(3), m(2), hard(1)];
        let pairs: Vec<_> = AlignedPairs::new(p0(0), &cigar).collect();
        assert_eq!(pairs.len(), 2);
        // qpos starts at 0 — hard clips don't advance it
        assert_eq!(pairs[0], match_m(0, p0(0)));
        assert_eq!(pairs[1], match_m(1, p0(1)));
    }

    // r[verify cigar.aligned_pairs.default_mode]
    #[test]
    fn soft_clips_hidden_by_default() {
        let cigar = [soft(3), m(2)];
        let pairs: Vec<_> = AlignedPairs::new(p0(0), &cigar).collect();
        // Only Match variants, no SoftClip
        assert_eq!(pairs.len(), 2);
        assert!(pairs.iter().all(|p| matches!(p, AlignedPair::Match { .. })));
    }

    // r[verify cigar.aligned_pairs.options]
    #[test]
    fn with_soft_clips_enables_them() {
        let cigar = [soft(3), m(2)];
        let pairs: Vec<_> = AlignedPairs::new(p0(0), &cigar).with_soft_clips().collect();
        assert_eq!(pairs.len(), 3, "SoftClip(summary) + 2M");
        assert_eq!(pairs[0], AlignedPair::SoftClip { qpos: 0, len: 3 });
        assert_eq!(pairs[1], match_m(3, p0(0)));
    }

    // r[verify cigar.aligned_pairs.options]
    #[test]
    fn full_yields_everything() {
        let cigar = [soft(1), m(1), pad(2), ins(1), unk(9, 2)];
        let pairs: Vec<_> = AlignedPairs::new(p0(0), &cigar).full().collect();
        // 5 ops: SoftClip + Match + Padding + Insertion + Unknown
        assert_eq!(pairs.len(), 5);
        assert!(matches!(pairs[0], AlignedPair::SoftClip { .. }));
        assert!(matches!(pairs[1], AlignedPair::Match { .. }));
        assert_eq!(pairs[2], AlignedPair::Padding { len: 2 });
        assert!(matches!(pairs[3], AlignedPair::Insertion { .. }));
        assert!(matches!(pairs[4], AlignedPair::Unknown { code: 9, len: 2 }));
    }

    #[test]
    fn full_yields_padding_and_unknown() {
        // Simpler: just pad and unknown
        let cigar = [m(1), pad(3), unk(13, 4), m(1)];
        let pairs: Vec<_> = AlignedPairs::new(p0(0), &cigar).full().collect();
        assert_eq!(pairs.len(), 4);
        assert!(matches!(pairs[0], AlignedPair::Match { .. }));
        assert_eq!(pairs[1], AlignedPair::Padding { len: 3 });
        assert!(matches!(pairs[2], AlignedPair::Unknown { code: 13, len: 4 }));
        assert!(matches!(pairs[3], AlignedPair::Match { .. }));
    }

    #[test]
    fn default_skips_padding_and_unknown() {
        let cigar = [m(1), pad(1), unk(13, 4), m(1)];
        let pairs: Vec<_> = AlignedPairs::new(p0(0), &cigar).collect();
        assert_eq!(pairs.len(), 2);
        assert!(matches!(pairs[0], AlignedPair::Match { .. }));
        assert!(matches!(pairs[1], AlignedPair::Match { .. }));
    }

    // r[verify cigar.aligned_pairs.default_mode]
    #[test]
    fn empty_cigar() {
        let pairs: Vec<_> = AlignedPairs::new(p0(0), &[]).collect();
        assert!(pairs.is_empty());
    }

    // r[verify cigar.aligned_pairs.zero_length_ops]
    #[test]
    fn zero_length_match_skipped() {
        let cigar = [m(0), m(2)];
        let pairs: Vec<_> = AlignedPairs::new(p0(0), &cigar).collect();
        assert_eq!(pairs.len(), 2);
        assert_eq!(pairs[0], match_m(0, p0(0)));
    }

    // r[verify cigar.aligned_pairs.zero_length_ops]
    #[test]
    fn zero_length_indel_skipped() {
        // Zero-length I/D/N/S must not yield summary variants with len=0
        let cigar = [m(1), ins(0), del(0), skip(0), Op::new(CigarOpType::SoftClip, 0), m(1)];
        let pairs: Vec<_> = AlignedPairs::new(p0(0), &cigar).with_soft_clips().collect();
        assert_eq!(pairs.len(), 2, "only the two Match positions remain");
        assert_eq!(pairs[0], match_m(0, p0(0)));
        assert_eq!(pairs[1], match_m(1, p0(1)));
    }

    // r[verify cigar.aligned_pairs.zero_length_ops]
    #[test]
    fn zero_length_padding_and_unknown_skipped() {
        let cigar = [m(1), pad(0), unk(9, 0), m(1)];
        let pairs: Vec<_> = AlignedPairs::new(p0(0), &cigar).full().collect();
        assert_eq!(pairs.len(), 2);
        assert!(matches!(pairs[0], AlignedPair::Match { .. }));
        assert!(matches!(pairs[1], AlignedPair::Match { .. }));
    }

    // ── D-I adjacency ─────────────────────────────────────────────────────

    // r[verify cigar.aligned_pairs.iterator]
    #[test]
    fn deletion_followed_by_insertion_yields_two_events() {
        // D-I: should yield Deletion then Insertion (two separate events)
        // Unlike the pileup engine's ComplexIndel which collapses them
        let cigar = [del(2), ins(3)];
        let pairs: Vec<_> = AlignedPairs::new(p0(100), &cigar).collect();
        assert_eq!(pairs.len(), 2);
        assert_eq!(pairs[0], AlignedPair::Deletion { rpos: p0(100), del_len: 2 });
        assert_eq!(pairs[1], AlignedPair::Insertion { qpos: 0, insert_len: 3 });
    }

    // ── Full CIGAR walk ───────────────────────────────────────────────────

    #[test]
    fn complex_cigar_walk() {
        // 2S + 3M + 1I + 2M + 1D + 1N + 2M
        let cigar = [soft(2), m(3), ins(1), m(2), del(1), skip(1), m(2)];
        let pairs: Vec<_> = AlignedPairs::new(p0(50), &cigar).collect();

        // Expected: 3M + I(summary) + 2M + D(summary) + N(summary) + 2M = 10 items
        let expected = vec![
            match_m(2, p0(50)), // after 2S
            match_m(3, p0(51)),
            match_m(4, p0(52)),
            AlignedPair::Insertion { qpos: 5, insert_len: 1 },
            match_m(6, p0(53)),
            match_m(7, p0(54)),
            AlignedPair::Deletion { rpos: p0(55), del_len: 1 },
            AlignedPair::RefSkip { rpos: p0(56), skip_len: 1 },
            match_m(8, p0(57)),
            match_m(9, p0(58)),
        ];
        assert_eq!(pairs, expected);
    }

    // ── Property: position monotonicity ───────────────────────────────────

    // r[verify cigar.aligned_pairs.position_monotonicity]
    #[test]
    fn qpos_is_monotone_nondecreasing() {
        // 2S3M1I2M2D1M
        let cigar = [soft(2), m(3), ins(1), m(2), del(2), m(1)];
        let pairs: Vec<_> = AlignedPairs::new(p0(0), &cigar).with_soft_clips().collect();

        // Just verify the qpos values in order
        let qposes: Vec<u32> = pairs
            .iter()
            .filter_map(|p| match p {
                AlignedPair::Match { qpos, .. } => Some(*qpos),
                AlignedPair::Insertion { qpos, .. } => Some(*qpos),
                AlignedPair::SoftClip { qpos, .. } => Some(*qpos),
                _ => None,
            })
            .collect();
        assert!(!qposes.is_empty());
        assert!(qposes.windows(2).all(|w| w[0] <= w[1]));
    }

    // r[verify cigar.aligned_pairs.position_monotonicity]
    #[test]
    fn rpos_is_monotone_nondecreasing() {
        let cigar = [m(3), del(2), skip(1), m(2)];
        let pairs: Vec<_> = AlignedPairs::new(p0(100), &cigar).collect();
        let rposes: Vec<u32> = pairs
            .iter()
            .filter_map(|p| match p {
                AlignedPair::Match { rpos, .. } => Some(**rpos),
                AlignedPair::Deletion { rpos, .. } => Some(**rpos),
                AlignedPair::RefSkip { rpos, .. } => Some(**rpos),
                _ => None,
            })
            .collect();
        assert!(!rposes.is_empty());
        assert!(rposes.windows(2).all(|w| w[0] <= w[1]));
    }

    // ── Property test oracle ──────────────────────────────────────────────

    /// Independent "oracle" implementation that walks the CIGAR the same way
    /// `AlignedPairs` does, but without the state-machine pattern.
    /// Used in property tests to verify equivalence.
    fn oracle_walk(
        cigar: &[CigarOp],
        rec_pos: Pos0,
        options: AlignedPairsOptions,
    ) -> Vec<AlignedPair> {
        let mut pairs = Vec::new();
        let mut qpos = 0u32;
        let mut rpos = rec_pos;

        for op in cigar {
            let len = op.len();
            // Skip zero-length ops uniformly (htslib behavior).
            if len == 0 && !matches!(op.op_type(), CigarOpType::HardClip) {
                continue;
            }
            let match_kind = match op.op_type() {
                CigarOpType::Match => Some(MatchKind::Match),
                CigarOpType::SeqMatch => Some(MatchKind::SeqMatch),
                CigarOpType::SeqMismatch => Some(MatchKind::SeqMismatch),
                _ => None,
            };
            match op.op_type() {
                CigarOpType::Match | CigarOpType::SeqMatch | CigarOpType::SeqMismatch => {
                    let kind = match_kind.expect("set above for M/=/X");
                    for i in 0..len {
                        pairs.push(AlignedPair::Match {
                            qpos: qpos.saturating_add(i),
                            rpos: advance_rpos(rpos, i),
                            kind,
                        });
                    }
                    qpos = qpos.saturating_add(len);
                    rpos = advance_rpos(rpos, len);
                }
                CigarOpType::Insertion => {
                    pairs.push(AlignedPair::Insertion { qpos, insert_len: len });
                    qpos = qpos.saturating_add(len);
                }
                CigarOpType::Deletion => {
                    pairs.push(AlignedPair::Deletion { rpos, del_len: len });
                    rpos = advance_rpos(rpos, len);
                }
                CigarOpType::RefSkip => {
                    pairs.push(AlignedPair::RefSkip { rpos, skip_len: len });
                    rpos = advance_rpos(rpos, len);
                }
                CigarOpType::SoftClip => {
                    if options.soft_clips {
                        pairs.push(AlignedPair::SoftClip { qpos, len });
                    }
                    qpos = qpos.saturating_add(len);
                }
                CigarOpType::HardClip => {
                    // never yielded, don't advance
                }
                CigarOpType::Padding => {
                    if options.padding_and_unknown {
                        pairs.push(AlignedPair::Padding { len });
                    }
                }
                CigarOpType::Unknown(code) => {
                    if options.padding_and_unknown {
                        pairs.push(AlignedPair::Unknown { code, len });
                    }
                }
            }
        }

        pairs
    }

    #[test]
    fn oracle_matches_iterator_default() {
        let cigars: &[&[CigarOp]] = &[
            &[m(5)],
            &[m(2), ins(1), m(2)],
            &[m(2), del(3), m(2)],
            &[soft(2), m(3)],
            &[m(2), skip(5), m(3)],
            &[del(2), ins(3)],
            &[hard(3), m(2)],
            &[m(1), pad(2), m(1)],
            &[],
            &[soft(1), m(1), ins(1), del(1), skip(1), m(1)],
            &[m(1), seq_match(2), seq_mismatch(1), m(1)],
        ];

        for cigar in cigars {
            let actual: Vec<_> = AlignedPairs::new(p0(100), cigar).collect();
            let expected = oracle_walk(cigar, p0(100), AlignedPairsOptions::default());
            assert_eq!(actual, expected, "mismatch for CIGAR: {cigar:?}");
        }
    }

    #[test]
    fn oracle_matches_iterator_full() {
        let cigar = &[soft(1), m(1), pad(2), ins(1), del(1), unk(9, 2), m(1)];
        let actual: Vec<_> = AlignedPairs::new(p0(0), cigar).full().collect();
        let mut opts = AlignedPairsOptions::default();
        opts.soft_clips = true;
        opts.padding_and_unknown = true;
        let expected = oracle_walk(cigar, p0(0), opts);
        assert_eq!(actual, expected);
    }

    // ── Proptest: random CIGAR walk matches oracle ────────────────────────

    mod proptests {
        use super::*;
        use proptest::prelude::*;

        fn arb_cigar_op() -> impl Strategy<Value = CigarOp> {
            // Generate valid CIGAR ops within reasonable bounds.
            // Includes len=0 to exercise the zero-length-skip behavior.
            (0u8..=14u8, 0u32..=50u32).prop_map(|(code, len)| {
                // Map codes 9..=14 to Unknown variant (reserved codes)
                let op_type = CigarOpType::from_bam(code);
                CigarOp::new(op_type, len)
            })
        }

        fn arb_cigar() -> impl Strategy<Value = Vec<CigarOp>> {
            proptest::collection::vec(arb_cigar_op(), 0..20)
        }

        fn arb_pos0() -> impl Strategy<Value = Pos0> {
            (0u32..=100_000u32).prop_map(|v| Pos0::new(v).unwrap())
        }

        // r[verify cigar.aligned_pairs.default_mode]
        proptest! {
            #[test]
            fn matches_oracle_default(
                cigar in arb_cigar(),
                pos in arb_pos0(),
            ) {
                let actual: Vec<_> = AlignedPairs::new(pos, &cigar).collect();
                let expected = oracle_walk(&cigar, pos, AlignedPairsOptions::default());
                prop_assert_eq!(actual, expected);
            }
        }

        proptest! {
            #[test]
            fn matches_oracle_full(
                cigar in arb_cigar(),
                pos in arb_pos0(),
            ) {
                let actual: Vec<_> = AlignedPairs::new(pos, &cigar).full().collect();
                let mut opts = AlignedPairsOptions::default();
                opts.soft_clips = true;
                opts.padding_and_unknown = true;
                let expected = oracle_walk(&cigar, pos, opts);
                prop_assert_eq!(actual, expected);
            }
        }

        // r[verify cigar.aligned_pairs.position_monotonicity]
        proptest! {
            #[test]
            fn ref_pos_monotone_nondecreasing(
                cigar in arb_cigar(),
            ) {
                let pairs: Vec<_> = AlignedPairs::new(Pos0::ZERO, &cigar).collect();
                let rposes: Vec<u32> = pairs.iter().filter_map(|p| match p {
                    AlignedPair::Match { rpos, .. } => Some(**rpos),
                    AlignedPair::Deletion { rpos, .. } => Some(**rpos),
                    AlignedPair::RefSkip { rpos, .. } => Some(**rpos),
                    _ => None,
                }).collect();
                prop_assert!(rposes.windows(2).all(|w| w[0] <= w[1]));
            }
        }

        proptest! {
            #[test]
            fn qpos_monotone_nondecreasing(
                cigar in arb_cigar(),
            ) {
                let pairs: Vec<_> = AlignedPairs::new(Pos0::ZERO, &cigar).full().collect();
                let qposes: Vec<u32> = pairs.iter().filter_map(|p| match p {
                    AlignedPair::Match { qpos, .. } => Some(*qpos),
                    AlignedPair::Insertion { qpos, .. } => Some(*qpos),
                    AlignedPair::SoftClip { qpos, .. } => Some(*qpos),
                    _ => None,
                }).collect();
                prop_assert!(qposes.windows(2).all(|w| w[0] <= w[1]));
            }
        }
    }

    // ── htslib comparison tests ──────────────────────────────────────────

    // r[verify cigar.aligned_pairs.htslib_equivalence]
    #[cfg(test)]
    mod htslib_tests {
        use super::*;

        /// Expand the actual `AlignedPairs` iterator into per-base
        /// `(Option<qpos>, Option<rpos>)` tuples — the shape rust-htslib's
        /// `aligned_pairs_full()` yields. Closes the loop: this test exercises
        /// the real iterator instead of comparing two seqair-internal oracles.
        ///
        /// Soft clips are enabled via `.with_soft_clips()` to match htslib;
        /// padding is intentionally excluded because rust-htslib's
        /// `aligned_pairs_full` panics on `Cigar::Pad`.
        fn expand_iterator(cigar: &[CigarOp], rec_pos: Pos0) -> Vec<(Option<u32>, Option<u32>)> {
            let mut result = Vec::new();
            for pair in AlignedPairs::new(rec_pos, cigar).with_soft_clips() {
                match pair {
                    AlignedPair::Match { qpos, rpos, .. } => {
                        result.push((Some(qpos), Some(*rpos)));
                    }
                    AlignedPair::Insertion { qpos, insert_len } => {
                        for i in 0..insert_len {
                            result.push((Some(qpos.saturating_add(i)), None));
                        }
                    }
                    AlignedPair::Deletion { rpos, del_len } => {
                        for i in 0..del_len {
                            result.push((None, Some(advance_rpos(rpos, i).as_u64() as u32)));
                        }
                    }
                    AlignedPair::RefSkip { rpos, skip_len } => {
                        for i in 0..skip_len {
                            result.push((None, Some(advance_rpos(rpos, i).as_u64() as u32)));
                        }
                    }
                    AlignedPair::SoftClip { qpos, len } => {
                        for i in 0..len {
                            result.push((Some(qpos.saturating_add(i)), None));
                        }
                    }
                    AlignedPair::Padding { .. } | AlignedPair::Unknown { .. } => {
                        // Default mode hides these; with_soft_clips alone
                        // doesn't enable them — included here for completeness.
                    }
                }
            }
            result
        }

        /// Drive `rust_htslib::bam::ext::BamRecordExtensions::aligned_pairs_full()`
        /// over the same CIGAR + pos and return its yields in the same shape
        /// `expand_iterator` produces.
        fn htslib_pairs(cigar_str: &str, pos: i64) -> Vec<(Option<u32>, Option<u32>)> {
            use rust_htslib::bam::ext::BamRecordExtensions;
            use rust_htslib::bam::record::Record;

            let cigar_string = build_htslib_cigar(cigar_str);
            let mut rec = Record::new();
            rec.set_pos(pos);
            rec.set_cigar(Some(&cigar_string));

            #[expect(
                clippy::cast_sign_loss,
                clippy::cast_possible_truncation,
                reason = "test fixtures: positions are non-negative and ≤ i32::MAX"
            )]
            let pairs = rec
                .aligned_pairs_full()
                .map(|arr| (arr[0].map(|v| v as u32), arr[1].map(|v| v as u32)))
                .collect();
            pairs
        }

        // r[verify cigar.aligned_pairs.htslib_equivalence]
        #[test]
        fn matches_htslib_aligned_pairs_full() {
            // (CIGAR string, pos). Cigar::Pad excluded — rust-htslib panics.
            // = and X are exercised here too — htslib collapses them into the
            // same (Some, Some) shape we produce.
            //
            // Hard-clip cases (`H`) explicitly verify htslib's "exclude
            // entirely, do not advance qpos" semantics — covered by spec
            // rule cigar.aligned_pairs.hard_clips.
            //
            // Adjacent indel cases (D-I, I-D, D-N, N-D) verify that the
            // summary form's separation of events matches htslib's expanded
            // form when both are flattened.
            let test_cases: &[(&str, i64)] = &[
                // ── Single-op CIGARs ──
                ("5M", 100),
                ("3=", 100),
                ("3X", 100),
                // ── Match with one indel ──
                ("2M1I2M", 100),
                ("2M3D2M", 100),
                ("3M2N4M", 100),
                // ── Soft clips ──
                ("5S3M", 0),
                ("3S5M2S", 0),
                ("3M5S", 0),
                // ── Hard clips: must NOT advance qpos and MUST be excluded ──
                ("3H5M", 100),
                ("5M2H", 100),
                ("3H5M2H", 100),
                ("3H2S5M", 100),
                // ── Hard + soft + match combinations ──
                ("2H3S5M2S1H", 100),
                // ── Adjacent indels ──
                ("2M1D1I2M", 100),
                ("2M1I2D3M", 100),
                ("3M1I1D2M", 100),
                ("2M1D2N1M", 100),
                // ── = and X interleaved ──
                ("3=2X1=", 50),
                ("2S3=1X2=2S", 0),
                ("1=1X1=1X1M", 50),
                // ── Long single ops to stress saturating arithmetic ──
                ("100M", 1_000_000),
                ("50M50D50M", 0),
            ];

            for &(cigar_str, pos) in test_cases {
                let our_pairs =
                    expand_iterator(&parse_cigar_string(cigar_str), Pos0::new(pos as u32).unwrap());
                let hts_pairs = htslib_pairs(cigar_str, pos);
                assert_eq!(
                    our_pairs, hts_pairs,
                    "CIGAR '{cigar_str}' at pos {pos}: seqair AlignedPairs (expanded) \
                     differs from rust-htslib::aligned_pairs_full"
                );
            }
        }

        // r[verify cigar.aligned_pairs.htslib_equivalence]
        #[test]
        fn match_kind_does_not_change_htslib_shape() {
            // Critical invariant: adding `kind: MatchKind` to AlignedPair::Match
            // must NOT change the (qpos, rpos) sequence — htslib parity is
            // about positions, not op kinds. M/=/X with the same lengths must
            // produce byte-identical expansions.
            let cigars = ["5M", "5=", "5X"];
            let pos = 100;

            let mut expansions = Vec::new();
            for s in cigars {
                expansions.push(expand_iterator(&parse_cigar_string(s), Pos0::new(pos).unwrap()));
            }
            assert_eq!(expansions[0], expansions[1], "5M and 5= must expand identically");
            assert_eq!(expansions[1], expansions[2], "5= and 5X must expand identically");
        }

        /// Build a rust-htslib CigarString from a CIGAR string like "5M2I3D".
        fn build_htslib_cigar(s: &str) -> rust_htslib::bam::record::CigarString {
            use rust_htslib::bam::record::{Cigar, CigarString};
            let mut cigars = Vec::new();
            let mut num = 0u32;
            for b in s.bytes() {
                match b {
                    b'0'..=b'9' => {
                        num = num.saturating_mul(10).saturating_add(u32::from(b - b'0'));
                    }
                    b'M' => {
                        cigars.push(Cigar::Match(num));
                        num = 0;
                    }
                    b'I' => {
                        cigars.push(Cigar::Ins(num));
                        num = 0;
                    }
                    b'D' => {
                        cigars.push(Cigar::Del(num));
                        num = 0;
                    }
                    b'N' => {
                        cigars.push(Cigar::RefSkip(num));
                        num = 0;
                    }
                    b'S' => {
                        cigars.push(Cigar::SoftClip(num));
                        num = 0;
                    }
                    b'H' => {
                        cigars.push(Cigar::HardClip(num));
                        num = 0;
                    }
                    b'P' => {
                        cigars.push(Cigar::Pad(num));
                        num = 0;
                    }
                    b'=' => {
                        cigars.push(Cigar::Equal(num));
                        num = 0;
                    }
                    b'X' => {
                        cigars.push(Cigar::Diff(num));
                        num = 0;
                    }
                    _ => {}
                }
            }
            CigarString(cigars)
        }

        // r[verify cigar.aligned_pairs.htslib_equivalence]
        mod htslib_proptests {
            use super::*;
            use proptest::prelude::*;

            /// Generate a CIGAR fragment as `(char, len)`. Excludes `P` because
            /// rust-htslib's `aligned_pairs_full` panics on `Cigar::Pad`. Excludes
            /// op codes 9..=14 because rust-htslib has no representation for them.
            /// Restricts to non-empty (1..=20) lengths because zero-length op skip
            /// is htslib-equivalent (its for-loop simply doesn't iterate) but our
            /// expansion explicitly drops len=0 events; testing both produces the
            /// same empty contribution but the `assert_eq!` on `Vec<...>` would
            /// pass trivially. Keeping len > 0 stresses the actual walking logic.
            fn arb_op_char_len() -> impl Strategy<Value = (char, u32)> {
                let chars = prop_oneof![
                    Just('M'),
                    Just('I'),
                    Just('D'),
                    Just('N'),
                    Just('S'),
                    Just('H'),
                    Just('='),
                    Just('X'),
                ];
                (chars, 1u32..=20u32)
            }

            fn arb_cigar_string() -> impl Strategy<Value = String> {
                proptest::collection::vec(arb_op_char_len(), 1..=10).prop_map(|ops| {
                    let mut s = String::new();
                    for (c, n) in ops {
                        s.push_str(&n.to_string());
                        s.push(c);
                    }
                    s
                })
            }

            proptest! {
                /// For random valid CIGARs (M/I/D/N/S/H/=/X with non-zero
                /// lengths) at random positions, seqair's expanded
                /// `AlignedPairs` output MUST match rust-htslib's
                /// `aligned_pairs_full()` exactly.
                ///
                /// This is the strongest possible parity test we can build:
                /// rust-htslib is the htslib-binding ground truth, and the
                /// proptest exercises the full CIGAR-walk state machine on
                /// inputs we never hand-wrote.
                #[test]
                fn matches_htslib_random_cigars(
                    cigar_str in arb_cigar_string(),
                    pos in 0i64..=100_000i64,
                ) {
                    let our_pairs = expand_iterator(
                        &parse_cigar_string(&cigar_str),
                        Pos0::new(pos as u32).unwrap(),
                    );
                    let hts_pairs = htslib_pairs(&cigar_str, pos);
                    prop_assert_eq!(
                        our_pairs,
                        hts_pairs,
                        "CIGAR '{}' at pos {} diverges from rust-htslib",
                        cigar_str,
                        pos,
                    );
                }
            }
        }

        /// Parse a CIGAR string like "5M2I3D" into our CigarOp vec.
        fn parse_cigar_string(s: &str) -> Vec<CigarOp> {
            let mut ops = Vec::new();
            let mut num = 0u32;
            for b in s.bytes() {
                match b {
                    b'0'..=b'9' => {
                        num = num.saturating_mul(10).saturating_add(u32::from(b - b'0'));
                    }
                    b'M' => {
                        ops.push(Op::new(CigarOpType::Match, num));
                        num = 0;
                    }
                    b'I' => {
                        ops.push(Op::new(CigarOpType::Insertion, num));
                        num = 0;
                    }
                    b'D' => {
                        ops.push(Op::new(CigarOpType::Deletion, num));
                        num = 0;
                    }
                    b'N' => {
                        ops.push(Op::new(CigarOpType::RefSkip, num));
                        num = 0;
                    }
                    b'S' => {
                        ops.push(Op::new(CigarOpType::SoftClip, num));
                        num = 0;
                    }
                    b'H' => {
                        ops.push(Op::new(CigarOpType::HardClip, num));
                        num = 0;
                    }
                    b'P' => {
                        ops.push(Op::new(CigarOpType::Padding, num));
                        num = 0;
                    }
                    b'=' => {
                        ops.push(Op::new(CigarOpType::SeqMatch, num));
                        num = 0;
                    }
                    b'X' => {
                        ops.push(Op::new(CigarOpType::SeqMismatch, num));
                        num = 0;
                    }
                    _ => {}
                }
            }
            ops
        }
    }

    // ── SlimRecord integration via push_raw → store → SlimRecord round-trip ──

    mod slim_record_tests {
        use super::super::super::owned_record::OwnedBamRecord;
        use super::super::super::record_store::RecordStore;
        use super::*;
        use seqair_types::{BamFlags, Base, BaseQuality};

        /// Build a record via OwnedBamRecord, serialize to BAM bytes, push into
        /// a fresh RecordStore via push_raw.
        fn store_with_record(pos: u32, cigar: &[CigarOp], seq_len: usize) -> RecordStore<()> {
            let rec = OwnedBamRecord::builder(0, Some(Pos0::new(pos).unwrap()), b"r".to_vec())
                .flags(BamFlags::empty())
                .cigar(cigar.to_vec())
                .seq(vec![Base::A; seq_len])
                .qual(vec![BaseQuality::from_byte(30); seq_len])
                .build()
                .unwrap();
            let mut buf = Vec::new();
            rec.to_bam_bytes(&mut buf).unwrap();

            let mut store = RecordStore::<()>::new();
            let _ = store.push_raw(&buf, &mut ()).unwrap();
            store
        }

        // r[verify cigar.aligned_pairs.slim_record]
        #[test]
        fn slim_record_round_trips_simple_match() {
            let cigar = [m(5)];
            let store = store_with_record(100, &cigar, 5);
            let rec = store.record(0);
            let pairs: Vec<_> = rec.aligned_pairs(&store).unwrap().collect();
            assert_eq!(pairs.len(), 5);
            assert_eq!(pairs[0], match_m(0, p0(100)));
            assert_eq!(pairs[4], match_m(4, p0(104)));
        }

        // r[verify cigar.aligned_pairs.slim_record]
        #[test]
        fn slim_record_matches_owned_record_walk() {
            // Build once, push into a store, then walk via both entry points.
            // The AlignedPair sequences must be identical.
            let cigar = vec![m(2), ins(1), m(2), del(1), m(2)];
            let owned = OwnedBamRecord::builder(0, Some(Pos0::new(200).unwrap()), b"r".to_vec())
                .flags(BamFlags::empty())
                .cigar(cigar.clone())
                .seq(vec![Base::A; 7])
                .qual(vec![BaseQuality::from_byte(30); 7])
                .build()
                .unwrap();
            let mut buf = Vec::new();
            owned.to_bam_bytes(&mut buf).unwrap();

            let mut store = RecordStore::<()>::new();
            let _ = store.push_raw(&buf, &mut ()).unwrap();

            let slim = store.record(0);
            let slim_pairs: Vec<_> = slim.aligned_pairs(&store).unwrap().collect();
            let owned_pairs: Vec<_> = owned.aligned_pairs().unwrap().collect();
            assert_eq!(slim_pairs, owned_pairs);
        }
    }
}
