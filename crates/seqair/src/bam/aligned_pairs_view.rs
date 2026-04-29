//! Layered iterator views over [`AlignedPairs`] that bundle per-event read
//! and reference data.
//!
//! The base [`AlignedPairs`](super::aligned_pairs::AlignedPairs) iterator yields
//! position-only events (`(qpos, rpos)` plus op type). Most callers want
//! richer information at each event — the bases and qualities the read placed
//! at this position, and the reference base it aligns against.
//!
//! Two adapters layer on top:
//!
//! - [`AlignedPairs::with_read`] — attaches the read's `seq` and `qual` slabs.
//!   `Match` yields gain `query: Base` and `qual: BaseQuality`. `Insertion` and
//!   `SoftClip` yields gain pre-sliced `query` / `qual` slices (the bases and
//!   quals of just that op's run — no off-by-one risk for callers).
//! - [`AlignedPairsWithRead::with_reference`] — further attaches a [`RefSeq`].
//!   `Match` yields gain `ref_base: Option<Base>` (`None` outside the loaded
//!   window). `Deletion` yields gain `ref_bases: Option<&[Base]>` (the reference
//!   bases that were deleted).
//!
//! ## Why a separate layer instead of one rich iterator?
//!
//! Pay-for-what-you-use: the base iterator stays cheap on the hot path. The
//! read layer adds two slab references; the reference layer adds one more.
//! Each layer's yield type is a distinct enum so the compiler enforces that
//! a caller using `with_reference` actually has reference data.
//!
//! ## Use case mapping
//!
//! - **Methylation calling** wants `with_read` + `with_reference`. The Match
//!   variant alone is sufficient; methylation logic is a one-line `match` arm.
//! - **SNV calling per record** wants the same — Match's `query` vs `ref_base`
//!   is the per-position evidence.
//! - **Indel calling** wants `with_read` only — the Insertion variant has
//!   pre-sliced bases and quals; Deletion has the length.
//! - **Hot-loop CIGAR walks** (e.g. NM/MD recompute) want the bare
//!   [`AlignedPairs`] and look up bases themselves only when needed.
//!
//! ## Example: methylation
//!
//! ```ignore
//! use seqair::bam::{AlignedPairWithRef, MatchKind};
//! use seqair_types::Base;
//!
//! for ev in slim.aligned_pairs_with_read(&store)?.with_reference(&ref_seq) {
//!     if let AlignedPairWithRef::Match {
//!         rpos,
//!         query: Base::T,
//!         ref_base: Some(Base::C),
//!         qual,
//!         ..
//!     } = ev
//!     {
//!         // C → T conversion: candidate methylation evidence at rpos.
//!         if qual.get().unwrap_or(0) >= 20 {
//!             // record evidence ...
//!         }
//!     }
//! }
//! ```

use super::aligned_pairs::{AlignedPair, AlignedPairs, MatchKind};
use super::pileup::RefSeq;
use seqair_types::{Base, BaseQuality, Pos0};

// ── AlignedPairWithRead ────────────────────────────────────────────────────

// r[impl cigar.aligned_pairs.with_read.types]
/// One event from [`AlignedPairsWithRead`] — an [`AlignedPair`] enriched with
/// the read's bases and qualities at the event's range.
///
/// Variants mirror [`AlignedPair`] one-to-one. `Match` positions carry a
/// single `query`/`qual`. `Insertion` and `SoftClip` carry slices of the
/// inserted/clipped run. `Deletion`/`RefSkip`/`Padding`/`Unknown` have no
/// read data, so they pass through unchanged.
///
/// The lifetime `'read` ties to the seq/qual slabs (typically a
/// [`RecordStore`](super::record_store::RecordStore)).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignedPairWithRead<'read> {
    /// M / = / X — read base aligned to a reference base.
    Match {
        qpos: u32,
        rpos: Pos0,
        kind: MatchKind,
        /// The read base at `qpos`. Always in bounds — `with_read`
        /// constructor enforces seq/CIGAR length consistency.
        query: Base,
        /// Phred quality at `qpos`. [`BaseQuality::UNAVAILABLE`] if `qual`
        /// was passed empty (BAM missing-qual sentinel).
        qual: BaseQuality,
    },
    /// I — `query` and `qual` are the inserted bases and their Phred scores.
    Insertion { qpos: u32, query: &'read [Base], qual: &'read [BaseQuality] },
    /// D — deletion from the reference. No read data.
    Deletion { rpos: Pos0, del_len: u32 },
    /// N — reference skip (e.g. intron). No read data.
    RefSkip { rpos: Pos0, skip_len: u32 },
    /// S — soft clip. `query` and `qual` are the clipped run.
    SoftClip { qpos: u32, query: &'read [Base], qual: &'read [BaseQuality] },
    /// P — padding. No read data.
    Padding { len: u32 },
    /// Reserved op code (9..=15). No read data.
    Unknown { code: u8, len: u32 },
}

// ── AlignedPairsWithRead iterator ──────────────────────────────────────────

// r[impl cigar.aligned_pairs.with_read.iterator]
/// Iterator over [`AlignedPairWithRead`] — produced by
/// [`AlignedPairs::with_read`](super::aligned_pairs::AlignedPairs::with_read).
///
/// Carries the read's `seq` and `qual` slabs (borrowed from a
/// [`RecordStore`](super::record_store::RecordStore)) and slices into them as
/// it walks the CIGAR.
///
/// Lifetimes: `'cigar` ties to the CIGAR slice; `'read` ties to the seq/qual
/// slabs. They are independent so callers can compose borrows from different
/// sources (e.g. CIGAR from a realignment buffer + seq/qual from the store).
#[derive(Debug)]
pub struct AlignedPairsWithRead<'cigar, 'read> {
    inner: AlignedPairs<'cigar>,
    seq: &'read [Base],
    qual: &'read [BaseQuality],
}

impl<'cigar, 'read> AlignedPairsWithRead<'cigar, 'read> {
    pub(super) fn new(
        inner: AlignedPairs<'cigar>,
        seq: &'read [Base],
        qual: &'read [BaseQuality],
    ) -> Self {
        Self { inner, seq, qual }
    }

    /// Forward `with_soft_clips()` to the inner iterator. Soft clips are
    /// hidden by default; opt in here just like on bare [`AlignedPairs`].
    pub fn with_soft_clips(mut self) -> Self {
        self.inner = self.inner.with_soft_clips();
        self
    }

    /// Forward `full()` to the inner iterator (yields `SoftClip`, `Padding`,
    /// and `Unknown` in addition to the default set).
    pub fn full(mut self) -> Self {
        self.inner = self.inner.full();
        self
    }

    /// Attach a reference window. Yields [`AlignedPairWithRef`] events, where
    /// Match gains `ref_base: Option<Base>` and Deletion gains
    /// `ref_bases: Option<&[Base]>`. `None` means the position is outside the
    /// loaded window (`RefSeq` only covers the queried region).
    pub fn with_reference<'ref_seq>(
        self,
        ref_seq: &'ref_seq RefSeq,
    ) -> AlignedPairsWithRef<'cigar, 'read, 'ref_seq> {
        AlignedPairsWithRef { inner: self, ref_seq }
    }

    fn attach_read(&self, pair: AlignedPair) -> AlignedPairWithRead<'read> {
        match pair {
            AlignedPair::Match { qpos, rpos, kind } => {
                let q = qpos as usize;
                AlignedPairWithRead::Match {
                    qpos,
                    rpos,
                    kind,
                    query: self.seq.get(q).copied().unwrap_or(Base::Unknown),
                    qual: self.qual.get(q).copied().unwrap_or(BaseQuality::UNAVAILABLE),
                }
            }
            AlignedPair::Insertion { qpos, insert_len } => {
                let (q, q_end) = range_for(qpos, insert_len);
                AlignedPairWithRead::Insertion {
                    qpos,
                    query: self.seq.get(q..q_end).unwrap_or(&[]),
                    qual: self.qual.get(q..q_end).unwrap_or(&[]),
                }
            }
            AlignedPair::Deletion { rpos, del_len } => {
                AlignedPairWithRead::Deletion { rpos, del_len }
            }
            AlignedPair::RefSkip { rpos, skip_len } => {
                AlignedPairWithRead::RefSkip { rpos, skip_len }
            }
            AlignedPair::SoftClip { qpos, len } => {
                let (q, q_end) = range_for(qpos, len);
                AlignedPairWithRead::SoftClip {
                    qpos,
                    query: self.seq.get(q..q_end).unwrap_or(&[]),
                    qual: self.qual.get(q..q_end).unwrap_or(&[]),
                }
            }
            AlignedPair::Padding { len } => AlignedPairWithRead::Padding { len },
            AlignedPair::Unknown { code, len } => AlignedPairWithRead::Unknown { code, len },
        }
    }
}

impl<'read> Iterator for AlignedPairsWithRead<'_, 'read> {
    type Item = AlignedPairWithRead<'read>;

    fn next(&mut self) -> Option<Self::Item> {
        let pair = self.inner.next()?;
        Some(self.attach_read(pair))
    }
}

#[inline]
fn range_for(qpos: u32, len: u32) -> (usize, usize) {
    let start = qpos as usize;
    let end = start.saturating_add(len as usize);
    (start, end)
}

// ── AlignedPairWithRef ─────────────────────────────────────────────────────

// r[impl cigar.aligned_pairs.with_reference.types]
/// One event from [`AlignedPairsWithRef`] — an [`AlignedPairWithRead`]
/// enriched with reference base lookups.
///
/// Match positions gain `ref_base: Option<Base>`; Deletion gains
/// `ref_bases: Option<&[Base]>`. `None` means the corresponding reference
/// position is outside the loaded [`RefSeq`] window.
///
/// Lifetimes: `'read` is the seq/qual slab borrow; `'ref_seq` is the
/// reference-window borrow. They are independent.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignedPairWithRef<'read, 'ref_seq> {
    /// M / = / X — query base, ref base, kind.
    Match {
        qpos: u32,
        rpos: Pos0,
        kind: MatchKind,
        query: Base,
        qual: BaseQuality,
        /// Reference base at `rpos`. `None` if `rpos` is outside the loaded
        /// `RefSeq` window. Inside the window, an `N` in the FASTA shows up
        /// as `Some(Base::Unknown)`.
        ref_base: Option<Base>,
    },
    /// I — inserted query bases and quals (no reference span).
    Insertion { qpos: u32, query: &'read [Base], qual: &'read [BaseQuality] },
    /// D — deleted reference bases. `ref_bases` is `None` if any position in
    /// the deletion span falls outside the loaded `RefSeq` window.
    Deletion { rpos: Pos0, del_len: u32, ref_bases: Option<&'ref_seq [Base]> },
    /// N — reference skip.
    RefSkip { rpos: Pos0, skip_len: u32 },
    /// S — soft clip (read-only; reference doesn't participate).
    SoftClip { qpos: u32, query: &'read [Base], qual: &'read [BaseQuality] },
    /// P — padding.
    Padding { len: u32 },
    /// Reserved op code.
    Unknown { code: u8, len: u32 },
}

// ── AlignedPairsWithRef iterator ───────────────────────────────────────────

// r[impl cigar.aligned_pairs.with_reference.iterator]
/// Iterator over [`AlignedPairWithRef`] — produced by
/// [`AlignedPairsWithRead::with_reference`].
///
/// Reuses the [`RefSeq`] already loaded by the reader (no copy or new window).
/// Borrows the slab references from `with_read` plus a borrow of the
/// reference window for the duration of iteration.
#[derive(Debug)]
pub struct AlignedPairsWithRef<'cigar, 'read, 'ref_seq> {
    inner: AlignedPairsWithRead<'cigar, 'read>,
    ref_seq: &'ref_seq RefSeq,
}

impl<'cigar, 'read, 'ref_seq> AlignedPairsWithRef<'cigar, 'read, 'ref_seq> {
    /// Forward `with_soft_clips()` to the inner iterator.
    pub fn with_soft_clips(mut self) -> Self {
        self.inner = self.inner.with_soft_clips();
        self
    }

    /// Forward `full()` to the inner iterator.
    pub fn full(mut self) -> Self {
        self.inner = self.inner.full();
        self
    }
}

impl<'read, 'ref_seq> Iterator for AlignedPairsWithRef<'_, 'read, 'ref_seq> {
    type Item = AlignedPairWithRef<'read, 'ref_seq>;

    fn next(&mut self) -> Option<Self::Item> {
        let pair = self.inner.next()?;
        Some(match pair {
            AlignedPairWithRead::Match { qpos, rpos, kind, query, qual } => {
                AlignedPairWithRef::Match {
                    qpos,
                    rpos,
                    kind,
                    query,
                    qual,
                    ref_base: self.ref_seq.try_base_at(rpos),
                }
            }
            AlignedPairWithRead::Insertion { qpos, query, qual } => {
                AlignedPairWithRef::Insertion { qpos, query, qual }
            }
            AlignedPairWithRead::Deletion { rpos, del_len } => AlignedPairWithRef::Deletion {
                rpos,
                del_len,
                ref_bases: self.ref_seq.range(rpos, del_len),
            },
            AlignedPairWithRead::RefSkip { rpos, skip_len } => {
                AlignedPairWithRef::RefSkip { rpos, skip_len }
            }
            AlignedPairWithRead::SoftClip { qpos, query, qual } => {
                AlignedPairWithRef::SoftClip { qpos, query, qual }
            }
            AlignedPairWithRead::Padding { len } => AlignedPairWithRef::Padding { len },
            AlignedPairWithRead::Unknown { code, len } => AlignedPairWithRef::Unknown { code, len },
        })
    }
}

// ── Tests ──────────────────────────────────────────────────────────────────

#[cfg(test)]
#[allow(clippy::arithmetic_side_effects, reason = "test arithmetic on small values")]
mod tests {
    use super::super::cigar::{CigarOp, CigarOpType};
    use super::super::owned_record::OwnedBamRecord;
    use super::super::record_store::RecordStore;
    use super::*;
    use seqair_types::{BamFlags, Base, BaseQuality, Pos0};
    use std::rc::Rc;

    fn op(t: CigarOpType, len: u32) -> CigarOp {
        CigarOp::new(t, len)
    }

    fn p0(v: u32) -> Pos0 {
        Pos0::new(v).unwrap()
    }

    /// Build a single-record store via OwnedBamRecord -> raw BAM -> push_raw.
    fn make_store(
        pos: Pos0,
        cigar: Vec<CigarOp>,
        seq: Vec<Base>,
        qual: Vec<BaseQuality>,
    ) -> RecordStore<()> {
        let owned = OwnedBamRecord::builder(0, Some(pos), b"r".to_vec())
            .flags(BamFlags::empty())
            .cigar(cigar)
            .seq(seq)
            .qual(qual)
            .build()
            .unwrap();
        let mut buf = Vec::new();
        owned.to_bam_bytes(&mut buf).unwrap();
        let mut store = RecordStore::<()>::new();
        let _ = store.push_raw(&buf, &mut ()).unwrap();
        store
    }

    // ── with_read: each variant ──────────────────────────────────────────

    // r[verify cigar.aligned_pairs.with_read.iterator]
    #[test]
    fn with_read_yields_query_base_and_qual_for_match() {
        let cigar = vec![op(CigarOpType::Match, 3)];
        let seq = vec![Base::A, Base::C, Base::G];
        let qual = vec![
            BaseQuality::from_byte(30),
            BaseQuality::from_byte(31),
            BaseQuality::from_byte(32),
        ];
        let store = make_store(p0(100), cigar, seq, qual);
        let rec = store.record(0);

        let events: Vec<_> = rec.aligned_pairs_with_read(&store).unwrap().collect();
        assert_eq!(events.len(), 3);
        assert_eq!(
            events[0],
            AlignedPairWithRead::Match {
                qpos: 0,
                rpos: p0(100),
                kind: MatchKind::Match,
                query: Base::A,
                qual: BaseQuality::from_byte(30),
            }
        );
        assert_eq!(
            events[2],
            AlignedPairWithRead::Match {
                qpos: 2,
                rpos: p0(102),
                kind: MatchKind::Match,
                query: Base::G,
                qual: BaseQuality::from_byte(32),
            }
        );
    }

    // r[verify cigar.aligned_pairs.with_read.iterator]
    #[test]
    fn with_read_slices_inserted_bases() {
        // 2M + 3I + 2M, seq = [A,C, T,G,T, A,A]
        let cigar = vec![
            op(CigarOpType::Match, 2),
            op(CigarOpType::Insertion, 3),
            op(CigarOpType::Match, 2),
        ];
        let seq = vec![Base::A, Base::C, Base::T, Base::G, Base::T, Base::A, Base::A];
        let qual: Vec<_> = (0..7).map(|i| BaseQuality::from_byte(20 + i)).collect();
        let store = make_store(p0(50), cigar, seq, qual);
        let rec = store.record(0);

        let events: Vec<_> = rec.aligned_pairs_with_read(&store).unwrap().collect();
        // 2 Match + 1 Insertion (summary) + 2 Match = 5 events
        assert_eq!(events.len(), 5);
        match events[2] {
            AlignedPairWithRead::Insertion { qpos, query, qual } => {
                assert_eq!(qpos, 2);
                assert_eq!(query, &[Base::T, Base::G, Base::T]);
                assert_eq!(qual.len(), 3);
                assert_eq!(qual[0].as_byte(), 22);
                assert_eq!(qual[2].as_byte(), 24);
            }
            other => panic!("expected Insertion, got {other:?}"),
        }
    }

    // r[verify cigar.aligned_pairs.with_read.iterator]
    #[test]
    fn with_read_slices_soft_clip_bases() {
        // 2S + 3M, seq = [N,N, A,C,G]
        let cigar = vec![op(CigarOpType::SoftClip, 2), op(CigarOpType::Match, 3)];
        let seq = vec![Base::Unknown, Base::Unknown, Base::A, Base::C, Base::G];
        let qual: Vec<_> = (0..5).map(|_| BaseQuality::from_byte(20)).collect();
        let store = make_store(p0(0), cigar, seq, qual);
        let rec = store.record(0);

        let events: Vec<_> =
            rec.aligned_pairs_with_read(&store).unwrap().with_soft_clips().collect();
        match events[0] {
            AlignedPairWithRead::SoftClip { qpos, query, qual } => {
                assert_eq!(qpos, 0);
                assert_eq!(query, &[Base::Unknown, Base::Unknown]);
                assert_eq!(qual.len(), 2);
            }
            other => panic!("expected SoftClip, got {other:?}"),
        }
    }

    // r[verify cigar.aligned_pairs.with_read.iterator]
    #[test]
    fn with_read_passes_deletion_through_unchanged() {
        let cigar = vec![
            op(CigarOpType::Match, 2),
            op(CigarOpType::Deletion, 3),
            op(CigarOpType::Match, 2),
        ];
        let seq = vec![Base::A; 4];
        let qual = vec![BaseQuality::from_byte(30); 4];
        let store = make_store(p0(100), cigar, seq, qual);
        let rec = store.record(0);

        let events: Vec<_> = rec.aligned_pairs_with_read(&store).unwrap().collect();
        // 2M + D(summary) + 2M = 5 events
        assert_eq!(events.len(), 5);
        assert_eq!(events[2], AlignedPairWithRead::Deletion { rpos: p0(102), del_len: 3 });
    }

    // ── with_reference: each variant ─────────────────────────────────────

    fn make_ref_seq(start: u32, bases: &[Base]) -> RefSeq {
        RefSeq::new(Rc::from(bases.to_vec()), p0(start))
    }

    // r[verify cigar.aligned_pairs.with_reference.iterator]
    #[test]
    fn with_reference_attaches_ref_base_to_match() {
        // Read 3M at ref pos 100, seq = [A,T,G]
        // Reference window starts at 100, contains [A,C,G,T]
        let cigar = vec![op(CigarOpType::Match, 3)];
        let seq = vec![Base::A, Base::T, Base::G];
        let qual = vec![BaseQuality::from_byte(30); 3];
        let store = make_store(p0(100), cigar, seq, qual);
        let rec = store.record(0);
        let ref_seq = make_ref_seq(100, &[Base::A, Base::C, Base::G, Base::T]);

        let events: Vec<_> =
            rec.aligned_pairs_with_read(&store).unwrap().with_reference(&ref_seq).collect();
        // qpos 0 / rpos 100: query=A, ref=A (match)
        // qpos 1 / rpos 101: query=T, ref=C (mismatch — variant evidence)
        // qpos 2 / rpos 102: query=G, ref=G (match)
        assert!(matches!(
            events[0],
            AlignedPairWithRef::Match { ref_base: Some(Base::A), query: Base::A, .. }
        ));
        assert!(matches!(
            events[1],
            AlignedPairWithRef::Match { ref_base: Some(Base::C), query: Base::T, .. }
        ));
        assert!(matches!(
            events[2],
            AlignedPairWithRef::Match { ref_base: Some(Base::G), query: Base::G, .. }
        ));
    }

    // r[verify cigar.aligned_pairs.with_reference.iterator]
    #[test]
    fn with_reference_returns_none_outside_loaded_window() {
        // Reference window starts at 200 — read at 100 is fully outside.
        let cigar = vec![op(CigarOpType::Match, 2)];
        let seq = vec![Base::A; 2];
        let qual = vec![BaseQuality::from_byte(30); 2];
        let store = make_store(p0(100), cigar, seq, qual);
        let rec = store.record(0);
        let ref_seq = make_ref_seq(200, &[Base::A, Base::C, Base::G]);

        let events: Vec<_> =
            rec.aligned_pairs_with_read(&store).unwrap().with_reference(&ref_seq).collect();
        for ev in &events {
            assert!(
                matches!(ev, AlignedPairWithRef::Match { ref_base: None, .. }),
                "expected ref_base=None, got {ev:?}"
            );
        }
    }

    // r[verify cigar.aligned_pairs.with_reference.iterator]
    #[test]
    fn with_reference_slices_ref_bases_for_deletion() {
        // 2M + 3D + 1M at pos 100, ref window covers 100..110
        let cigar = vec![
            op(CigarOpType::Match, 2),
            op(CigarOpType::Deletion, 3),
            op(CigarOpType::Match, 1),
        ];
        let seq = vec![Base::A; 3];
        let qual = vec![BaseQuality::from_byte(30); 3];
        let store = make_store(p0(100), cigar, seq, qual);
        let rec = store.record(0);
        // ref bases at positions 100..110 = [A,C,G,T,A,C,G,T,A,C]
        let ref_bases: Vec<Base> = b"ACGTACGTAC".iter().map(|&b| Base::from(b)).collect();
        let ref_seq = make_ref_seq(100, &ref_bases);

        let events: Vec<_> =
            rec.aligned_pairs_with_read(&store).unwrap().with_reference(&ref_seq).collect();
        // events: 2 Match + 1 Deletion + 1 Match = 4
        assert_eq!(events.len(), 4);
        match events[2] {
            AlignedPairWithRef::Deletion { rpos, del_len, ref_bases: Some(bases) } => {
                assert_eq!(rpos, p0(102));
                assert_eq!(del_len, 3);
                // ref pos 102, 103, 104 = G,T,A
                assert_eq!(bases, &[Base::G, Base::T, Base::A]);
            }
            other => panic!("expected Deletion with ref_bases, got {other:?}"),
        }
    }

    // r[verify cigar.aligned_pairs.with_reference.iterator]
    #[test]
    fn with_reference_returns_none_for_partial_deletion_overlap() {
        // 2M + 3D + 1M at pos 100, but ref window only covers 100..103
        // Deletion span 102..105 → 102 is in, 103/104 out → ref_bases=None
        let cigar = vec![
            op(CigarOpType::Match, 2),
            op(CigarOpType::Deletion, 3),
            op(CigarOpType::Match, 1),
        ];
        let seq = vec![Base::A; 3];
        let qual = vec![BaseQuality::from_byte(30); 3];
        let store = make_store(p0(100), cigar, seq, qual);
        let rec = store.record(0);
        let ref_seq = make_ref_seq(100, &[Base::A, Base::C, Base::G]);

        let events: Vec<_> =
            rec.aligned_pairs_with_read(&store).unwrap().with_reference(&ref_seq).collect();
        match events[2] {
            AlignedPairWithRef::Deletion { ref_bases: None, .. } => {}
            other => panic!("expected Deletion with ref_bases=None, got {other:?}"),
        }
    }

    // r[verify cigar.aligned_pairs.with_reference.iterator]
    #[test]
    fn with_reference_then_with_soft_clips_chains_correctly() {
        // 2S + 3M at pos 100. Without `.with_soft_clips()` the SoftClip event
        // is skipped; with it, we get 1 SoftClip + 3 Match. Verifies the
        // forwarder on AlignedPairsWithRef (not just AlignedPairsWithRead).
        let cigar = vec![op(CigarOpType::SoftClip, 2), op(CigarOpType::Match, 3)];
        let seq = vec![Base::Unknown, Base::Unknown, Base::A, Base::C, Base::G];
        let qual = vec![BaseQuality::from_byte(20); 5];
        let store = make_store(p0(100), cigar, seq, qual);
        let rec = store.record(0);
        let ref_seq = make_ref_seq(100, &[Base::A, Base::C, Base::G]);

        // Default: SoftClip hidden — 3 events.
        let default_events: Vec<_> =
            rec.aligned_pairs_with_read(&store).unwrap().with_reference(&ref_seq).collect();
        assert_eq!(default_events.len(), 3);
        assert!(default_events.iter().all(|e| matches!(e, AlignedPairWithRef::Match { .. })));

        // With soft clips chained AFTER with_reference — must propagate through.
        let chained_events: Vec<_> = rec
            .aligned_pairs_with_read(&store)
            .unwrap()
            .with_reference(&ref_seq)
            .with_soft_clips()
            .collect();
        assert_eq!(chained_events.len(), 4, "1 SoftClip + 3 Match");
        assert!(matches!(chained_events[0], AlignedPairWithRef::SoftClip { .. }));
        assert!(matches!(chained_events[1], AlignedPairWithRef::Match { .. }));
        assert!(matches!(chained_events[3], AlignedPairWithRef::Match { .. }));
    }

    // r[verify cigar.aligned_pairs.with_reference.iterator]
    #[test]
    fn with_reference_full_yields_padding_and_unknown() {
        // Confirm `.full()` propagates through both layers — Padding/Unknown
        // events surface even in the with_reference iterator.
        let cigar = vec![
            op(CigarOpType::Match, 1),
            op(CigarOpType::Padding, 2),
            op(CigarOpType::Unknown(9), 3),
            op(CigarOpType::Match, 1),
        ];
        let seq = vec![Base::A, Base::C];
        let qual = vec![BaseQuality::from_byte(30); 2];
        let store = make_store(p0(50), cigar, seq, qual);
        let rec = store.record(0);
        let ref_seq = make_ref_seq(50, &[Base::A, Base::C]);

        let events: Vec<_> =
            rec.aligned_pairs_with_read(&store).unwrap().with_reference(&ref_seq).full().collect();
        assert_eq!(events.len(), 4);
        assert!(matches!(events[0], AlignedPairWithRef::Match { .. }));
        assert_eq!(events[1], AlignedPairWithRef::Padding { len: 2 });
        assert!(matches!(events[2], AlignedPairWithRef::Unknown { code: 9, len: 3 }));
        assert!(matches!(events[3], AlignedPairWithRef::Match { .. }));
    }

    // ── Validation errors ─────────────────────────────────────────────────

    #[test]
    fn with_read_rejects_seq_shorter_than_cigar_query_len() {
        // CIGAR claims 3M (3 query bases) but seq has only 2 bases. The
        // earlier silent Base::Unknown substitution would mask this; the
        // validation now surfaces a typed error instead.
        let cigar = [op(CigarOpType::Match, 3)];
        let seq = vec![Base::A, Base::C];
        let qual = vec![BaseQuality::from_byte(30); 2];
        let result = AlignedPairs::new(p0(0), &cigar).with_read(&seq, &qual);
        match result {
            Err(super::super::aligned_pairs::AlignedPairsError::CigarSeqLengthMismatch {
                cigar_qlen: 3,
                seq_len: 2,
            }) => {}
            other => panic!("expected CigarSeqLengthMismatch{{3,2}}, got {other:?}"),
        }
    }

    #[test]
    fn with_read_rejects_seq_longer_than_cigar_query_len() {
        // Inverse mismatch: CIGAR 2M but seq 3 bases.
        let cigar = [op(CigarOpType::Match, 2)];
        let seq = vec![Base::A, Base::C, Base::G];
        let qual = vec![BaseQuality::from_byte(30); 3];
        let result = AlignedPairs::new(p0(0), &cigar).with_read(&seq, &qual);
        assert!(matches!(
            result,
            Err(super::super::aligned_pairs::AlignedPairsError::CigarSeqLengthMismatch {
                cigar_qlen: 2,
                seq_len: 3,
            })
        ));
    }

    #[test]
    fn with_read_rejects_qual_length_mismatch() {
        // qual.len() must equal seq.len() OR be zero (BAM missing-qual sentinel).
        let cigar = [op(CigarOpType::Match, 3)];
        let seq = vec![Base::A; 3];
        let qual = vec![BaseQuality::from_byte(30); 2]; // length 2 ≠ 3
        let result = AlignedPairs::new(p0(0), &cigar).with_read(&seq, &qual);
        assert!(matches!(
            result,
            Err(super::super::aligned_pairs::AlignedPairsError::SeqQualLengthMismatch {
                seq_len: 3,
                qual_len: 2,
            })
        ));
    }

    #[test]
    fn with_read_accepts_empty_qual_as_missing_sentinel() {
        // qual.len() == 0 represents BAM's "qual missing" encoding (0xFF bytes).
        // Should succeed; iteration uses BaseQuality::UNAVAILABLE.
        let cigar = [op(CigarOpType::Match, 3)];
        let seq = vec![Base::A; 3];
        let qual: Vec<BaseQuality> = vec![];
        let it = AlignedPairs::new(p0(0), &cigar).with_read(&seq, &qual).unwrap();
        let events: Vec<_> = it.collect();
        for ev in events {
            match ev {
                AlignedPairWithRead::Match { qual, .. } => {
                    assert_eq!(qual, BaseQuality::UNAVAILABLE);
                }
                _ => panic!("expected only Match events for 3M, got {ev:?}"),
            }
        }
    }

    // ── End-to-end: methylation calling pattern ──────────────────────────

    // r[verify cigar.aligned_pairs.with_reference.iterator]
    #[test]
    fn methylation_pattern_finds_c_to_t_conversion() {
        // Simulated TAPS scenario: ref has C at pos 100, read shows T → methylation.
        let cigar = vec![op(CigarOpType::Match, 3)];
        let seq = vec![Base::T, Base::A, Base::G];
        let qual = vec![BaseQuality::from_byte(35); 3];
        let store = make_store(p0(100), cigar, seq, qual);
        let rec = store.record(0);
        let ref_seq = make_ref_seq(100, &[Base::C, Base::A, Base::G]);

        let mut conversions = 0;
        for ev in rec.aligned_pairs_with_read(&store).unwrap().with_reference(&ref_seq) {
            if let AlignedPairWithRef::Match {
                ref_base: Some(Base::C), query: Base::T, qual, ..
            } = ev
                && qual.get().unwrap_or(0) >= 20
            {
                conversions += 1;
            }
        }
        assert_eq!(conversions, 1, "expected exactly one C→T conversion at pos 100");
    }

    // ── Property test: with_read pass-through ─────────────────────────────

    mod proptests {
        use super::super::super::aligned_pairs::AlignedPairs;
        use super::*;
        use proptest::prelude::*;

        fn arb_cigar_op() -> impl Strategy<Value = CigarOp> {
            // Match the bare-iterator proptest: full op-code range (0..=14
            // covers M/I/D/N/S/H/P/=/X plus reserved Unknown variants) and
            // include zero-length ops (0..=10) so zero-length skip behavior is
            // exercised through the with_read layer too.
            (0u8..=14u8, 0u32..=10u32).prop_map(|(code, len)| {
                let t = CigarOpType::from_bam(code);
                CigarOp::new(t, len)
            })
        }

        proptest! {
            /// `with_read` must yield exactly the same number of events as the
            /// underlying iterator, in the same order — it only enriches.
            #[test]
            fn with_read_preserves_event_sequence(
                ops in proptest::collection::vec(arb_cigar_op(), 0..15),
                start in 0u32..=10_000u32,
            ) {
                // Build a synthetic seq/qual sized to the query-consuming length.
                let qlen: u32 = ops
                    .iter()
                    .filter(|o| matches!(
                        o.op_type(),
                        CigarOpType::Match
                            | CigarOpType::Insertion
                            | CigarOpType::SoftClip
                            | CigarOpType::SeqMatch
                            | CigarOpType::SeqMismatch
                    ))
                    .map(|o| o.len())
                    .sum();
                let seq: Vec<Base> = (0..qlen).map(|_| Base::A).collect();
                let qual: Vec<BaseQuality> = (0..qlen).map(|_| BaseQuality::from_byte(30)).collect();

                let pos = Pos0::new(start).unwrap();
                let bare: Vec<_> = AlignedPairs::new(pos, &ops).with_soft_clips().collect();
                let rich_iter = AlignedPairs::new(pos, &ops)
                    .with_soft_clips()
                    .with_read(&seq, &qual);
                let rich = match rich_iter {
                    Ok(it) => it.collect::<Vec<_>>(),
                    Err(e) => {
                        // We sized seq/qual to qlen exactly, so validation
                        // should always pass. If it ever fails, surface the
                        // diagnostic loudly rather than silently skipping.
                        prop_assert!(false, "with_read validation failed unexpectedly: {e}");
                        return Ok(());
                    }
                };

                prop_assert_eq!(bare.len(), rich.len());
                for (b, r) in bare.iter().zip(rich.iter()) {
                    let consistent = match (b, r) {
                        (
                            AlignedPair::Match { qpos: bq, rpos: br, kind: bk },
                            AlignedPairWithRead::Match { qpos: rq, rpos: rr, kind: rk, .. },
                        ) => bq == rq && br == rr && bk == rk,
                        (
                            AlignedPair::Insertion { qpos: bq, insert_len: bl },
                            AlignedPairWithRead::Insertion { qpos: rq, query, .. },
                        ) => bq == rq && (*bl) as usize == query.len(),
                        (
                            AlignedPair::Deletion { rpos: br, del_len: bl },
                            AlignedPairWithRead::Deletion { rpos: rr, del_len: rl },
                        ) => br == rr && bl == rl,
                        (
                            AlignedPair::RefSkip { rpos: br, skip_len: bl },
                            AlignedPairWithRead::RefSkip { rpos: rr, skip_len: rl },
                        ) => br == rr && bl == rl,
                        (
                            AlignedPair::SoftClip { qpos: bq, len: bl },
                            AlignedPairWithRead::SoftClip { qpos: rq, query, .. },
                        ) => bq == rq && (*bl) as usize == query.len(),
                        (
                            AlignedPair::Padding { len: bl },
                            AlignedPairWithRead::Padding { len: rl },
                        ) => bl == rl,
                        (
                            AlignedPair::Unknown { code: bc, len: bl },
                            AlignedPairWithRead::Unknown { code: rc, len: rl },
                        ) => bc == rc && bl == rl,
                        _ => false,
                    };
                    prop_assert!(consistent, "variant mismatch: bare={:?} rich={:?}", b, r);
                }
            }
        }
    }
}
