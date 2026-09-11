//! Addressing a record in a [`RecordStore`]: the index that travels, and the
//! handle that can be read.
//!
//! The two are separate because they answer different questions. A
//! [`RecordIdx`] is a plain number — it is what a [`SlimRecord`] stores for its
//! mate, what a [`PileupAlignment`] reports, and what a window query yields, so
//! it must be `Copy`, storable, and free of borrows. A [`RecordRef`] is an
//! index that has been *checked against a particular store and borrowed from
//! it*, so reading through it cannot fail and cannot go stale.
//!
//! [`PileupAlignment`]: crate::bam::pileup::PileupAlignment

use super::{
    aux::Aux,
    cigar::CigarOp,
    record_store::{RecordStore, SlimRecord},
};
use seqair_types::{Base, BaseQuality, Pos0, QPos};
use std::{num::NonZeroU32, ops::Range};

// r[impl record_store.record_idx]
/// An index into a [`RecordStore`]'s records.
///
/// Store indices used to travel as bare `u32`, the same type as query offsets,
/// slab offsets, lengths and depths — none of which is a record index. This is
/// the type that says which of those a value is.
///
/// `u32::MAX` is deliberately not representable. A store cannot hold that many
/// records (a push mints indices with `u32::try_from(len)`, and every record
/// occupies at least one byte of the `u32`-addressed name slab), so nothing is
/// lost — and the gap left behind is a niche, which makes `Option<RecordIdx>`
/// four bytes wide, the same as the bare index. That is what lets
/// [`SlimRecord`] and [`PileupAlignment`] spell "no linked mate" as `None`
/// rather than as a sentinel a caller could mistake for a real index.
///
/// An index means something only against the store that minted it, and only
/// until that store is cleared or refilled. Resolve one with
/// [`RecordStore::record`], which answers `None` past the store's end and
/// otherwise hands back a [`RecordRef`] that cannot go stale.
///
/// [`PileupAlignment`]: crate::bam::pileup::PileupAlignment
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct RecordIdx(NonZeroU32);

impl RecordIdx {
    /// The first record of a store.
    pub const ZERO: Self = Self(NonZeroU32::MIN);

    // r[impl record_store.record_idx]
    /// The index `idx`, or `None` for `u32::MAX` — the one value a record
    /// index can never take (see the type docs).
    ///
    /// The stored form is `idx + 1`, which is order-preserving, so the derived
    /// [`Ord`] sorts by index. The binary searches in
    /// [`PileupColumn::position_of`] and in the window query depend on that.
    ///
    /// [`PileupColumn::position_of`]: crate::bam::pileup::PileupColumn::position_of
    #[inline]
    #[must_use]
    pub const fn new(idx: u32) -> Option<Self> {
        match NonZeroU32::new(idx.wrapping_add(1)) {
            Some(bumped) => Some(Self(bumped)),
            None => None,
        }
    }

    /// The index of the record in slot `n`. The bridge from a `usize` loop
    /// counter or a `Vec` length, in place of an `as u32` truncation.
    #[inline]
    #[must_use]
    pub fn from_usize(n: usize) -> Option<Self> {
        Self::new(u32::try_from(n).ok()?)
    }

    /// The raw index.
    #[inline]
    #[must_use]
    pub const fn get(self) -> u32 {
        self.0.get().wrapping_sub(1)
    }

    /// The raw index, for reaching into the store's slabs.
    #[inline]
    #[must_use]
    pub const fn as_usize(self) -> usize {
        self.get() as usize
    }
}

impl std::fmt::Debug for RecordIdx {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "RecordIdx({})", self.get())
    }
}

impl std::fmt::Display for RecordIdx {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.get())
    }
}

// r[impl record_store.record_ref]
/// A record, resolved against the store that holds it.
///
/// [`RecordStore::record`] is the only way to make one, and it checks the index
/// once. Everything reachable from the handle — the record's fields, its slab
/// slices, its mate — is then infallible, because the two things that could
/// have gone wrong have both been ruled out:
///
/// * **Out of range.** Checked at construction, and the index cannot change
///   afterwards.
/// * **Stale.** The handle borrows the store, so the store cannot be cleared,
///   refilled, sorted, deduplicated or pushed to while the handle exists — the
///   borrow checker rejects it. An index kept across regions used to resolve
///   silently against whatever record now sits at that slot; a handle cannot
///   outlive the records it was checked against.
///
/// It is also paired with *one* store, not merely with some store: the slab
/// accessors are methods on the handle rather than functions taking a store, so
/// there is no call that could supply a different one. That is the failure
/// [`SlimRecord::qname`] and its siblings return a [`RecordAccessError`] for;
/// through a handle it cannot be spelled.
///
/// ```
/// # use seqair::bam::{RecordIdx, RecordStore};
/// # fn demo(store: &RecordStore) -> Option<()> {
/// let rec = store.record(RecordIdx::ZERO)?;
/// let (pos, name) = (rec.pos, rec.qname());   // derefs to `SlimRecord`
/// let mate_start = rec.mate().map(|m| m.pos); // links resolve through the handle
/// # let _ = (pos, name, mate_start);
/// # Some(())
/// # }
/// ```
///
/// The store cannot be disturbed while a handle is out — the error code is
/// pinned, so these fail for the borrow and not for a typo:
///
/// ```compile_fail,E0502
/// # use seqair::bam::{RecordIdx, RecordStore};
/// # fn demo(store: &mut RecordStore) -> Option<()> {
/// let rec = store.record(RecordIdx::ZERO)?;
/// store.clear();               // cannot borrow `*store` as mutable
/// let _ = rec.qname();
/// # Some(())
/// # }
/// ```
///
/// ```compile_fail,E0502
/// # use seqair::bam::{RecordIdx, RecordStore};
/// # fn demo(store: &mut RecordStore, raw: &[u8]) -> Option<()> {
/// let rec = store.record(RecordIdx::ZERO)?;
/// store.push_raw(raw, &mut ()).ok()?;   // a refill would move the slabs
/// let _ = rec.seq();
/// # Some(())
/// # }
/// ```
///
/// ```compile_fail,E0502
/// # use seqair::bam::{RecordIdx, RecordStore};
/// # fn demo(store: &mut RecordStore) -> Option<()> {
/// let rec = store.record(RecordIdx::ZERO)?;
/// store.sort_by_pos();         // a reorder would leave the index pointing elsewhere
/// let _ = rec.pos;
/// # Some(())
/// # }
/// ```
///
/// And a handle cannot outlive the store it was taken from:
///
/// ```compile_fail,E0515
/// # use seqair::bam::{RecordIdx, RecordRef, RecordStore};
/// fn escape() -> Option<RecordRef<'static>> {
///     let store = RecordStore::new();
///     store.record(RecordIdx::ZERO)    // cannot return a value referencing a local
/// }
/// ```
///
/// [`RecordAccessError`]: crate::bam::record_store::RecordAccessError
pub struct RecordRef<'store, U = ()> {
    store: &'store RecordStore<U>,
    slim: &'store SlimRecord,
    idx: RecordIdx,
}

// Hand-written: `derive` would demand `U: Clone`, and the handle is three
// references' worth of `Copy` data whatever `U` is.
impl<U> Clone for RecordRef<'_, U> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<U> Copy for RecordRef<'_, U> {}

impl<U> std::fmt::Debug for RecordRef<'_, U> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("RecordRef")
            .field("idx", &self.idx)
            .field("pos", &self.slim.pos)
            .finish_non_exhaustive()
    }
}

/// The record's own fields (`pos`, `end_pos`, `flags`, `mapq`, …) read straight
/// off the handle.
impl<U> std::ops::Deref for RecordRef<'_, U> {
    type Target = SlimRecord;

    fn deref(&self) -> &Self::Target {
        self.slim
    }
}

impl<'store, U> RecordRef<'store, U> {
    /// Pair a record with the store it was found in. Crate-private: the check
    /// that `idx` addresses `slim` inside `store` happens in
    /// [`RecordStore::record`], and this is what records that it did.
    pub(crate) fn new(
        store: &'store RecordStore<U>,
        slim: &'store SlimRecord,
        idx: RecordIdx,
    ) -> Self {
        Self { store, slim, idx }
    }

    /// Where this record sits in the store — the plain index, for storing or
    /// comparing against one a [`PileupAlignment`] reported.
    ///
    /// [`PileupAlignment`]: crate::bam::pileup::PileupAlignment
    #[must_use]
    pub fn idx(&self) -> RecordIdx {
        self.idx
    }

    /// The record itself, borrowed for as long as the store is.
    ///
    /// [`Deref`](std::ops::Deref) already exposes its fields; this is for when
    /// the shorter borrow of `&self` is not enough.
    #[must_use]
    pub fn slim(&self) -> &'store SlimRecord {
        self.slim
    }

    /// The store this record lives in.
    #[must_use]
    pub fn store(&self) -> &'store RecordStore<U> {
        self.store
    }

    // r[impl record_store.field_access]
    /// The read's QNAME bytes.
    #[must_use]
    pub fn qname(&self) -> &'store [u8] {
        self.store.qname_of(self.slim)
    }

    // r[impl record_store.field_access]
    /// The read's CIGAR operations.
    #[must_use]
    pub fn cigar(&self) -> &'store [CigarOp] {
        self.store.cigar_of(self.slim)
    }

    // r[impl record_store.field_access]
    /// The read's decoded sequence.
    pub fn seq(&self) -> &'store [Base] {
        self.store.seq_of(self.slim)
    }

    // r[impl bam.record.seq_at]
    /// The base at `qpos`, or [`Base::Unknown`] past the end of the read —
    /// which is also what a record with `SEQ=*` reports at every offset.
    pub fn base_at(&self, qpos: QPos) -> Base {
        self.store.base_at_of(self.slim, qpos)
    }

    // r[impl types.base_quality.field_type]
    /// The read's per-base quality scores.
    #[must_use]
    pub fn qual(&self) -> &'store [BaseQuality] {
        self.store.qual_of(self.slim)
    }

    // r[impl bam.record.raw_aux]
    /// The read's auxiliary tag bytes, in BAM wire format.
    #[must_use]
    pub fn aux(&self) -> &'store [u8] {
        self.store.aux_of(self.slim)
    }

    // r[impl bam.record.aux_wrapper]
    /// The auxiliary tags, as a queryable [`Aux`] wrapper.
    #[must_use]
    pub fn aux_tags(&self) -> Aux<'store> {
        Aux::new(self.aux())
    }

    // r[impl record_store.extras.access]
    /// The per-record extra computed by the store's customize value.
    #[must_use]
    pub fn extra(&self) -> &'store U {
        self.store.extra_of(self.slim)
    }

    // r[impl record_store.link_mates+2]
    /// This record's linked mate, when [`RecordStore::link_mates`] paired them.
    ///
    /// The link is a store index, so this is where it should be followed: the
    /// mate is resolved against the same store, through the same check.
    #[must_use]
    pub fn mate(&self) -> Option<RecordRef<'store, U>> {
        self.store.record(self.slim.mate_idx()?)
    }

    // r[impl record_store.mate_overlap]
    /// The half-open reference interval both mates cover, or `None` when this
    /// record has no linked mate. Empty when the two alignments are disjoint.
    #[must_use]
    pub fn mate_overlap(&self) -> Option<Range<Pos0>> {
        let mate = self.mate()?;
        let start = self.slim.pos.max(mate.pos);
        let last = self.slim.end_pos.min(mate.end_pos);
        if last < start {
            return Some(start..start);
        }
        let end =
            last.checked_add_offset(seqair_types::Offset::new(1)).unwrap_or_else(Pos0::max_value);
        Some(start..end)
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]
mod tests {
    use super::*;
    use crate::bam::record_store::RecordStore;
    use crate::bam::test_util::ri;
    use seqair_types::BamFlags;

    // r[verify record_store.record_idx]
    #[test]
    fn record_idx_round_trips_every_representable_index() {
        for raw in [0u32, 1, 2, 255, 65_535, u32::MAX - 1] {
            let idx = RecordIdx::new(raw).expect("representable");
            assert_eq!(idx.get(), raw);
            assert_eq!(idx.as_usize(), raw as usize);
            assert_eq!(RecordIdx::from_usize(raw as usize), Some(idx));
        }
        assert_eq!(RecordIdx::ZERO.get(), 0);
    }

    // r[verify record_store.record_idx]
    /// The one excluded value, and the niche it buys: the `Option` must cost
    /// nothing over the index itself, which is what lets `mate_idx` drop its
    /// sentinel.
    #[test]
    fn u32_max_is_not_an_index_and_that_is_the_niche() {
        assert_eq!(RecordIdx::new(u32::MAX), None);
        assert_eq!(RecordIdx::from_usize(u32::MAX as usize), None);
        assert_eq!(RecordIdx::from_usize(usize::MAX), None);
        assert_eq!(size_of::<Option<RecordIdx>>(), size_of::<RecordIdx>());
        assert_eq!(size_of::<RecordIdx>(), size_of::<u32>());
    }

    // r[verify record_store.record_idx]
    /// Sorting by the stored form must be sorting by the index — the column
    /// binary searches in `PileupColumn::position_of` rely on it, and the
    /// stored form is not the index.
    #[test]
    fn ordering_follows_the_index_not_the_representation() {
        let raws = [0u32, 1, 7, 300, 70_000, u32::MAX - 2, u32::MAX - 1];
        for pair in raws.windows(2) {
            let (lo, hi) = (pair[0], pair[1]);
            assert!(ri(lo) < ri(hi), "{lo} must sort before {hi}");
        }
        let mut shuffled: Vec<RecordIdx> = [9u32, 2, 40, 0, 7].into_iter().map(ri).collect();
        shuffled.sort_unstable();
        assert_eq!(shuffled.iter().map(|i| i.get()).collect::<Vec<_>>(), [0, 2, 7, 9, 40]);
    }

    /// Two mates of one template, linked, so a `RecordRef` has somewhere to go.
    fn paired_store() -> RecordStore {
        let mut store = RecordStore::new();
        let bases = [Base::A, Base::C, Base::G, Base::T];
        let quals = [30u8, 31, 32, 33];
        for (pos, end, mate_pos) in [(100u32, 103u32, 102i32), (102, 105, 100)] {
            store
                .push_fields(
                    Pos0::new(pos).unwrap(),
                    Pos0::new(end).unwrap(),
                    BamFlags::from(0x1u16),
                    60,
                    4,
                    0,
                    b"frag",
                    &[CigarOp::new(crate::bam::cigar::CigarOpType::Match, 4)],
                    &bases,
                    &quals,
                    b"XAi\x07\x00\x00\x00",
                    0,
                    0,
                    mate_pos,
                    0,
                    &mut (),
                )
                .unwrap()
                .unwrap();
        }
        let stats = store.link_mates();
        assert_eq!(stats.pairs, 1, "the fixture must be one linked pair");
        store
    }

    // r[verify record_store.record_ref]
    /// One check, then every slab is readable off the handle — and the record's
    /// own fields come through `Deref` without a second lookup.
    #[test]
    fn a_handle_reads_every_slab_of_its_own_record() {
        let store = paired_store();
        let rec = store.record(ri(0)).expect("record 0 exists");

        assert_eq!(rec.idx(), ri(0));
        assert_eq!(rec.pos, Pos0::new(100).unwrap(), "fields arrive through Deref");
        assert_eq!(rec.qname(), b"frag");
        assert_eq!(rec.seq(), [Base::A, Base::C, Base::G, Base::T]);
        assert_eq!(rec.qual().len(), 4);
        assert_eq!(rec.cigar().len(), 1);
        assert_eq!(rec.aux(), b"XAi\x07\x00\x00\x00");
        assert_eq!(
            rec.aux_tags().get::<i32>(b"XA").ok(),
            Some(7),
            "the Aux wrapper reads the same bytes"
        );
        assert_eq!(rec.base_at(QPos::new(2)), Base::G);
        assert_eq!(rec.base_at(QPos::new(99)), Base::Unknown, "past the read, not the next record");
        assert!(std::ptr::eq(rec.store(), &store));
    }

    // r[verify record_store.record_ref]
    /// A mate link is a store index, so it is followed through the same check —
    /// and the hop lands on the other half of the pair, not on index 1 of
    /// whatever store happens to be at hand.
    #[test]
    fn a_handle_follows_its_mate_link() {
        let store = paired_store();
        let first = store.record(ri(0)).expect("record 0 exists");
        let mate = first.mate().expect("the pair is linked");

        assert_eq!(mate.idx(), ri(1));
        assert_eq!(mate.pos, Pos0::new(102).unwrap());
        assert_eq!(mate.mate().map(|back| back.idx()), Some(ri(0)), "the link is symmetric");
        assert_eq!(first.mate_overlap(), Some(Pos0::new(102).unwrap()..Pos0::new(104).unwrap()));
        assert_eq!(first.mate_overlap(), mate.mate_overlap(), "both mates see one overlap");
    }

    // r[verify record_store.record_ref]
    /// A record with no linked mate has no hop to make, and no overlap.
    #[test]
    fn an_unlinked_record_has_no_mate_and_no_overlap() {
        let mut store = paired_store();
        store.sort_by_pos();
        let rec = store.record(ri(0)).expect("record 0 exists");
        assert!(rec.mate().is_none(), "sort_by_pos clears the links");
        assert!(rec.mate_overlap().is_none());
    }
}
