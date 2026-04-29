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

use super::cigar::{CigarOp, CigarOpType};
use super::record_store::{RecordAccessError, RecordStore, SlimRecord};
use seqair_types::Pos0;

// ── AlignedPair ────────────────────────────────────────────────────────────

// r[impl cigar.aligned_pairs.types]
/// A typed aligned pair — one event from walking a record's CIGAR.
///
/// Unlike htslib's `(Option<i64>, Option<i64>)`, this enum uses typed variants:
/// `Match` carries both positions, `Insertion` carries `qpos`+len, `Deletion`
/// carries `rpos`+len, `RefSkip` carries `rpos`+len. Match ops expand
/// per-position; I/D/N/S ops yield once per op (summary form).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignedPair {
    /// M / = / X op — read base aligned to a reference base.
    /// Yields one variant per consumed position.
    Match { qpos: u32, rpos: Pos0 },

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

    /// P op — padding. Hidden by default; opt in via [`.full()`](AlignedPairs::full).
    Padding,

    /// Reserved op code (9..=15) we tolerate on read. `code` is the raw 4-bit op.
    /// Yielded by [`.full()`](AlignedPairs::full) only.
    Unknown { code: u8, len: u32 },
}

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
    /// Yield one Match per remaining base.
    PerBase {
        remaining: u32,
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
}

impl Iterator for AlignedPairs<'_> {
    type Item = AlignedPair;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // ── Expand multi-base ops (M/=/X) ──
            match std::mem::replace(&mut self.expanding, ExpandingState::None) {
                ExpandingState::PerBase { remaining } if remaining > 0 => {
                    self.expanding = ExpandingState::PerBase { remaining: remaining - 1 };
                    let qpos = self.qpos;
                    let rpos = self.rpos;
                    self.qpos = self.qpos.saturating_add(1);
                    self.rpos = advance_rpos(self.rpos, 1);
                    return Some(AlignedPair::Match { qpos, rpos });
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
            match op.op_type() {
                CigarOpType::Match | CigarOpType::SeqMatch | CigarOpType::SeqMismatch => {
                    if len == 0 {
                        continue;
                    }
                    // Enter expansion state; loop back to yield first position.
                    self.expanding = ExpandingState::PerBase { remaining: len };
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
                        return Some(AlignedPair::Padding);
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
    pub fn aligned_pairs<'store, U>(
        &self,
        store: &'store RecordStore<U>,
    ) -> Result<AlignedPairs<'store>, RecordAccessError> {
        Ok(AlignedPairs::new(self.pos, self.cigar(store)?))
    }
}

// ── OwnedBamRecord integration ────────────────────────────────────────────

// r[impl cigar.aligned_pairs.owned_record]
impl super::owned_record::OwnedBamRecord {
    /// Walk this record's CIGAR, yielding typed [`AlignedPair`]s.
    ///
    /// Uses the owned CIGAR directly — no store access needed.
    /// For unmapped records (pos = -1), `rpos` starts at 0.
    pub fn aligned_pairs(&self) -> AlignedPairs<'_> {
        let rpos = if self.pos >= 0 {
            Pos0::new(self.pos as u32).unwrap_or(Pos0::ZERO)
        } else {
            Pos0::ZERO
        };
        AlignedPairs::new(rpos, &self.cigar)
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
    fn pad() -> Op {
        Op::new(CigarOpType::Padding, 0)
    }
    fn unk(code: u8, len: u32) -> Op {
        Op::new(CigarOpType::Unknown(code), len)
    }

    fn p0(v: u32) -> Pos0 {
        Pos0::new(v).unwrap()
    }

    // ── Basic unit tests ──────────────────────────────────────────────────

    // r[verify cigar.aligned_pairs.default_mode]
    #[test]
    fn simple_match() {
        let cigar = [m(5)];
        let pairs: Vec<_> = AlignedPairs::new(p0(100), &cigar).collect();
        assert_eq!(pairs.len(), 5);
        assert_eq!(pairs[0], AlignedPair::Match { qpos: 0, rpos: p0(100) });
        assert_eq!(pairs[4], AlignedPair::Match { qpos: 4, rpos: p0(104) });
    }

    // r[verify cigar.aligned_pairs.insertion_qpos]
    #[test]
    fn with_insertion() {
        let cigar = [m(2), ins(1), m(2)];
        let pairs: Vec<_> = AlignedPairs::new(p0(100), &cigar).collect();
        assert_eq!(pairs.len(), 5, "2M + 1I(summary) + 2M = 5 items");
        assert_eq!(pairs[0], AlignedPair::Match { qpos: 0, rpos: p0(100) });
        assert_eq!(pairs[1], AlignedPair::Match { qpos: 1, rpos: p0(101) });
        // Insertion summary: qpos = 2 (first inserted base), len = 1
        assert_eq!(pairs[2], AlignedPair::Insertion { qpos: 2, insert_len: 1 });
        assert_eq!(pairs[3], AlignedPair::Match { qpos: 3, rpos: p0(102) });
        assert_eq!(pairs[4], AlignedPair::Match { qpos: 4, rpos: p0(103) });
    }

    // r[verify cigar.aligned_pairs.deletion_rpos]
    #[test]
    fn with_deletion() {
        let cigar = [m(2), del(3), m(2)];
        let pairs: Vec<_> = AlignedPairs::new(p0(100), &cigar).collect();
        assert_eq!(pairs.len(), 5, "2M + 1D(summary) + 2M = 5 items");
        assert_eq!(pairs[0], AlignedPair::Match { qpos: 0, rpos: p0(100) });
        assert_eq!(pairs[1], AlignedPair::Match { qpos: 1, rpos: p0(101) });
        // Deletion summary: rpos = 102 (start of deletion), del_len = 3
        assert_eq!(pairs[2], AlignedPair::Deletion { rpos: p0(102), del_len: 3 });
        assert_eq!(pairs[3], AlignedPair::Match { qpos: 2, rpos: p0(105) });
        assert_eq!(pairs[4], AlignedPair::Match { qpos: 3, rpos: p0(106) });
    }

    #[test]
    fn with_refskip() {
        let cigar = [m(2), skip(5), m(3)];
        let pairs: Vec<_> = AlignedPairs::new(p0(100), &cigar).collect();
        assert_eq!(pairs.len(), 6, "2M + 1N(summary) + 3M = 6 items");
        assert_eq!(pairs[2], AlignedPair::RefSkip { rpos: p0(102), skip_len: 5 });
        assert_eq!(pairs[3], AlignedPair::Match { qpos: 2, rpos: p0(107) });
    }

    // r[verify cigar.aligned_pairs.qpos_semantics]
    #[test]
    fn qpos_includes_soft_clips() {
        // 5S + 3M: first match should have qpos = 5
        let cigar = [soft(5), m(3)];
        let pairs: Vec<_> = AlignedPairs::new(p0(0), &cigar).collect();
        assert_eq!(pairs.len(), 3);
        assert_eq!(pairs[0], AlignedPair::Match { qpos: 5, rpos: p0(0) });
        assert_eq!(pairs[2], AlignedPair::Match { qpos: 7, rpos: p0(2) });
    }

    // r[verify cigar.aligned_pairs.hard_clips]
    #[test]
    fn hard_clips_ignored() {
        let cigar = [hard(3), m(2), hard(1)];
        let pairs: Vec<_> = AlignedPairs::new(p0(0), &cigar).collect();
        assert_eq!(pairs.len(), 2);
        // qpos starts at 0 — hard clips don't advance it
        assert_eq!(pairs[0], AlignedPair::Match { qpos: 0, rpos: p0(0) });
        assert_eq!(pairs[1], AlignedPair::Match { qpos: 1, rpos: p0(1) });
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
        assert_eq!(pairs[1], AlignedPair::Match { qpos: 3, rpos: p0(0) });
    }

    // r[verify cigar.aligned_pairs.options]
    #[test]
    fn full_yields_everything() {
        let cigar = [soft(1), m(1), pad(), ins(1), unk(9, 2)];
        let pairs: Vec<_> = AlignedPairs::new(p0(0), &cigar).full().collect();
        // 5 ops: SoftClip + Match + Padding + Insertion + Unknown
        assert_eq!(pairs.len(), 5);
        assert!(matches!(pairs[0], AlignedPair::SoftClip { .. }));
        assert!(matches!(pairs[1], AlignedPair::Match { .. }));
        assert!(matches!(pairs[2], AlignedPair::Padding));
        assert!(matches!(pairs[3], AlignedPair::Insertion { .. }));
        assert!(matches!(pairs[4], AlignedPair::Unknown { code: 9, len: 2 }));
    }

    #[test]
    fn full_yields_padding_and_unknown() {
        // Simpler: just pad and unknown
        let cigar = [m(1), pad(), unk(13, 4), m(1)];
        let pairs: Vec<_> = AlignedPairs::new(p0(0), &cigar).full().collect();
        assert_eq!(pairs.len(), 4);
        assert!(matches!(pairs[0], AlignedPair::Match { .. }));
        assert!(matches!(pairs[1], AlignedPair::Padding));
        assert!(matches!(pairs[2], AlignedPair::Unknown { code: 13, len: 4 }));
        assert!(matches!(pairs[3], AlignedPair::Match { .. }));
    }

    #[test]
    fn default_skips_padding_and_unknown() {
        let cigar = [m(1), pad(), unk(13, 4), m(1)];
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

    #[test]
    fn zero_length_op_skipped() {
        // Degenerate: a 0-length M op should be skipped
        let cigar = [m(0), m(2)];
        let pairs: Vec<_> = AlignedPairs::new(p0(0), &cigar).collect();
        assert_eq!(pairs.len(), 2);
        assert_eq!(pairs[0], AlignedPair::Match { qpos: 0, rpos: p0(0) });
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
            AlignedPair::Match { qpos: 2, rpos: p0(50) }, // after 2S
            AlignedPair::Match { qpos: 3, rpos: p0(51) },
            AlignedPair::Match { qpos: 4, rpos: p0(52) },
            AlignedPair::Insertion { qpos: 5, insert_len: 1 },
            AlignedPair::Match { qpos: 6, rpos: p0(53) },
            AlignedPair::Match { qpos: 7, rpos: p0(54) },
            AlignedPair::Deletion { rpos: p0(55), del_len: 1 },
            AlignedPair::RefSkip { rpos: p0(56), skip_len: 1 },
            AlignedPair::Match { qpos: 8, rpos: p0(57) },
            AlignedPair::Match { qpos: 9, rpos: p0(58) },
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
            match op.op_type() {
                CigarOpType::Match | CigarOpType::SeqMatch | CigarOpType::SeqMismatch => {
                    for i in 0..len {
                        pairs.push(AlignedPair::Match {
                            qpos: qpos.saturating_add(i),
                            rpos: advance_rpos(rpos, i),
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
                        pairs.push(AlignedPair::Padding);
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
            &[m(1), pad(), m(1)],
            &[],
            &[soft(1), m(1), ins(1), del(1), skip(1), m(1)],
        ];

        for cigar in cigars {
            let actual: Vec<_> = AlignedPairs::new(p0(100), cigar).collect();
            let expected = oracle_walk(cigar, p0(100), AlignedPairsOptions::default());
            assert_eq!(actual, expected, "mismatch for CIGAR: {cigar:?}");
        }
    }

    #[test]
    fn oracle_matches_iterator_full() {
        let cigar = &[soft(1), m(1), pad(), ins(1), del(1), unk(9, 2), m(1)];
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
            // Generate valid CIGAR ops within reasonable bounds
            (0u8..=14u8, 1u32..=50u32).prop_map(|(code, len)| {
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

        /// Walk a CIGAR in expanded form (one yield per consumed base, like
        /// htslib's `aligned_pairs_full`) and return `(Option<qpos>, Option<rpos>)` tuples.
        fn expanded_pairs(cigar: &[CigarOp], rec_pos: Pos0) -> Vec<(Option<u32>, Option<Pos0>)> {
            let mut result = Vec::new();
            let mut qpos = 0u32;
            let mut rpos = rec_pos;

            for op in cigar {
                let len = op.len();
                match op.op_type() {
                    CigarOpType::Match | CigarOpType::SeqMatch | CigarOpType::SeqMismatch => {
                        for i in 0..len {
                            result.push((Some(qpos + i), Some(advance_rpos(rpos, i))));
                        }
                        qpos = qpos.saturating_add(len);
                        rpos = advance_rpos(rpos, len);
                    }
                    CigarOpType::Insertion => {
                        for i in 0..len {
                            result.push((Some(qpos + i), None));
                        }
                        qpos = qpos.saturating_add(len);
                    }
                    CigarOpType::Deletion | CigarOpType::RefSkip => {
                        for i in 0..len {
                            result.push((None, Some(advance_rpos(rpos, i))));
                        }
                        rpos = advance_rpos(rpos, len);
                    }
                    CigarOpType::SoftClip => {
                        for i in 0..len {
                            result.push((Some(qpos + i), None));
                        }
                        qpos = qpos.saturating_add(len);
                    }
                    CigarOpType::HardClip => {
                        // htslib excludes hard clips entirely
                    }
                    CigarOpType::Padding => {
                        // htslib includes padding as (None, None) pairs if requested
                    }
                    CigarOpType::Unknown(_) => {
                        // htslib treats unknown ops like padding
                    }
                }
            }

            result
        }

        /// Compare our expanded output against htslib's `aligned_pairs_full`.
        /// This test requires the `rust-htslib` dev-dependency.
        #[test]
        fn matches_htslib_aligned_pairs_full() {
            // Build a set of CIGARs and verify expanded output matches
            // what htslib would produce for similarly-configured records.
            //
            // We verify equivalence by constructing htslib records with the
            // same CIGAR and comparing position by position.

            use rust_htslib::bam::ext::BamRecordExtensions;
            use rust_htslib::bam::record::Record;

            let test_cases: &[(&str, i64)] = &[
                // (CIGAR string, pos)
                ("5M", 100),
                ("2M1I2M", 100),
                ("2M3D2M", 100),
                ("5S3M", 0),
                ("3S5M2S", 0),
                ("3M2N4M", 100),
                ("2M1D1I2M", 100),
                // htslib panics on Cigar::Pad in aligned_pairs_full, so we skip P ops here
                ("2M1I2D3M", 100),
            ];

            for &(cigar_str, pos) in test_cases {
                // Build htslib record
                let cigar_string = build_htslib_cigar(cigar_str);
                let mut rec = Record::new();
                rec.set_pos(pos);
                rec.set_cigar(Some(&cigar_string));

                // Get htslib's aligned_pairs_full
                let hts_pairs: Vec<_> = rec
                    .aligned_pairs_full()
                    .map(|arr| (arr[0].map(|v| v as u32), arr[1].map(|v| v as u32)))
                    .collect();

                // Get our expanded pairs
                let cigar_ops: Vec<CigarOp> = {
                    // Parse the same CIGAR string into our ops
                    parse_cigar_string(cigar_str)
                };
                let rec_pos = Pos0::new(pos as u32).unwrap();
                let our_pairs: Vec<_> = expanded_pairs(&cigar_ops, rec_pos)
                    .into_iter()
                    .map(|(q, r)| (q, r.map(|p| *p)))
                    .collect();

                assert_eq!(
                    our_pairs, hts_pairs,
                    "CIGAR '{cigar_str}' at pos {pos}: seqair expanded pairs differ from htslib"
                );
            }
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
}
