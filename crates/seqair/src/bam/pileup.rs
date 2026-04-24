//! Iterate over pileup columns. [`PileupEngine`] yields [`PileupColumn`]s with pre-extracted
//! flat fields per active read. Supports read filtering and per-position max-depth.

use seqair_types::{BamFlags, Base, BaseQuality, Offset, Pos0, Strand, strand_from_flags};
// Rc is used only for RefSeq (reference sequence), not for BAM records.
// PileupEngine is intentionally !Send due to RecordFilter: Box<dyn Fn(...)>.
use std::rc::Rc;

use crate::utils::TraceErr;

use super::{
    cigar::{CigarMapping, CigarPosInfo},
    record_store::RecordStore,
};

pub struct RefSeq {
    bases: Rc<[Base]>,
    start_pos: Pos0,
}

impl std::fmt::Debug for RefSeq {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("RefSeq")
            .field("start_pos", &self.start_pos)
            .field("len", &self.bases.len())
            .finish()
    }
}

impl RefSeq {
    pub fn new(bases: Rc<[Base]>, start_pos: Pos0) -> Self {
        Self { bases, start_pos }
    }

    pub fn base_at(&self, pos: Pos0) -> Base {
        let Some(offset) = pos.as_i64().checked_sub(self.start_pos.as_i64()) else {
            return Base::Unknown;
        };
        if offset < 0 {
            return Base::Unknown;
        }
        #[expect(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            reason = "offset is non-negative (checked above) and bounded by reference sequence length; safe to cast on supported platforms"
        )]
        self.bases.get(offset as usize).copied().unwrap_or(Base::Unknown)
    }
}

/// Filter function receiving (flags, `aux_bytes`) for the candidate record.
// r[impl flags.filter_signature]
type RecordFilter = Box<dyn Fn(BamFlags, &[u8]) -> bool>;

// r[impl pileup.active_set]
// r[impl pileup.zero_refspan_reads]
// r[impl pileup.soft_clip_at_position]
// r[impl perf.no_sorted_indices]
// r[impl perf.avoid_redundant_arena_get+2]
// r[impl perf.cigar_no_to_vec]
// r[impl perf.reuse_alignment_vec+2]
// r[impl pileup.extras.generic_param]
pub struct PileupEngine<U = ()> {
    store: RecordStore<U>,
    current_pos: Pos0,
    region_end: Pos0,
    next_entry: usize,
    /// Hot field: checked every column during retain. Stored separately so the
    /// retain loop strides 4 bytes instead of the full `ActiveRecord` size (~144 bytes).
    active_end_pos: Vec<Pos0>,
    /// Cold fields: only accessed for records that survive retain.
    active: Vec<ActiveRecord>,
    max_depth: Option<u32>,
    filter: Option<RecordFilter>,
    ref_seq: Option<RefSeq>,
    // Profiling counters — u32 is sufficient for single-region pileups (max ~250M positions).
    // For whole-genome streaming, these saturate at u32::MAX which is acceptable for debug output.
    columns_produced: u32,
    max_active_depth: u32,
}

#[derive(Debug)]
struct ActiveRecord {
    record_idx: u32,
    cigar: CigarMapping,
    // Cached from SlimRecord to avoid store lookups in the hot loop
    flags: BamFlags,
    strand: Strand,
    mapq: u8,
    seq_len: u32,
    matching_bases: u32,
    indel_bases: u32,
}

// r[impl pileup.column_contents]
// r[impl pileup.htslib_compat]
// r[impl pileup.lending_iterator]
/// A single pileup column, borrowing the record store for the duration of its use.
///
/// Returned by [`PileupEngine::pileups`]. Valid until the next call to `pileups`
/// on the same engine. Access per-record extras or slab data (qname, aux) via
/// [`alignments`](Self::alignments), which yields [`AlignmentView`] wrappers, or
/// via [`store`](Self::store) directly.
pub struct PileupColumn<'store, U = ()> {
    pos: Pos0,
    reference_base: Base,
    alignments: Vec<PileupAlignment>,
    store: &'store RecordStore<U>,
}

impl<U> std::fmt::Debug for PileupColumn<'_, U> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PileupColumn")
            .field("pos", &self.pos)
            .field("reference_base", &self.reference_base)
            .field("depth", &self.alignments.len())
            .finish_non_exhaustive()
    }
}

impl<'store, U> PileupColumn<'store, U> {
    #[must_use]
    pub fn pos(&self) -> Pos0 {
        self.pos
    }

    // r[impl pileup_indel.depth_includes_all]
    #[must_use]
    pub fn depth(&self) -> usize {
        self.alignments.len()
    }

    /// Iterate alignments with store access for per-record extras, qname, and aux.
    ///
    /// Each yielded [`AlignmentView`] derefs to [`PileupAlignment`] for the flat
    /// per-position fields, and exposes [`extra`](AlignmentView::extra),
    /// [`qname`](AlignmentView::qname), and [`aux`](AlignmentView::aux) for slab data.
    pub fn alignments(&self) -> impl Iterator<Item = AlignmentView<'_, 'store, U>> + '_ {
        let store = self.store;
        self.alignments.iter().map(move |aln| AlignmentView { aln, store })
    }

    /// Iterate the raw alignments without store access (equivalent to the old API).
    pub fn raw_alignments(&self) -> impl Iterator<Item = &PileupAlignment> + '_ {
        self.alignments.iter()
    }

    pub fn reference_base(&self) -> Base {
        self.reference_base
    }

    /// Borrow the record store for custom slab access (e.g., record fields by index).
    pub fn store(&self) -> &'store RecordStore<U> {
        self.store
    }

    /// Count of alignments with a query base at this position.
    ///
    /// Unlike [`depth`](Self::depth), deletions and ref-skips are not counted.
    /// Use `match_depth` when you need the number of reads that actually cover
    /// the position with a base (e.g., for allele frequency calculations).
    #[must_use]
    pub fn match_depth(&self) -> usize {
        self.alignments.iter().filter(|a| a.qpos().is_some()).count()
    }
}

/// A view over a single alignment in a column, with access to the record store.
///
/// Derefs to [`PileupAlignment`] so existing field and method access works
/// unchanged. Adds [`extra`](Self::extra), [`qname`](Self::qname), and
/// [`aux`](Self::aux) for data stored in the record store's slabs.
// r[impl pileup.alignment_view]
pub struct AlignmentView<'a, 'store, U> {
    aln: &'a PileupAlignment,
    store: &'store RecordStore<U>,
}

impl<U> std::fmt::Debug for AlignmentView<'_, '_, U> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("AlignmentView").field("aln", &self.aln).finish_non_exhaustive()
    }
}

impl<'a, 'store, U> AlignmentView<'a, 'store, U> {
    /// The per-record extra, as computed by the extras provider / `with_extras`.
    pub fn extra(&self) -> &'store U {
        self.store.extra(self.aln.record_idx())
    }

    /// The read's QNAME bytes in the store's name slab.
    pub fn qname(&self) -> &'store [u8] {
        self.store.qname(self.aln.record_idx())
    }

    /// The raw BAM aux bytes for this record.
    pub fn aux(&self) -> &'store [u8] {
        self.store.aux(self.aln.record_idx())
    }

    /// The underlying [`PileupAlignment`] — use when you want to call
    /// [`Clone::clone`] or otherwise escape the deref coercion.
    pub fn alignment(&self) -> &'a PileupAlignment {
        self.aln
    }

    /// The store this view references.
    pub fn store(&self) -> &'store RecordStore<U> {
        self.store
    }
}

impl<U> std::ops::Deref for AlignmentView<'_, '_, U> {
    type Target = PileupAlignment;
    fn deref(&self) -> &PileupAlignment {
        self.aln
    }
}

/// What a read shows at a specific reference position in the pileup.
///
/// `base`, `qual`, and `qpos` are only available for `Match` and `Insertion`
/// variants — the type system prevents reading a base from a deletion.
///
/// Use the convenience accessors [`PileupAlignment::base`], [`PileupAlignment::qual`],
/// and [`PileupAlignment::qpos`] when you only need one field, or match exhaustively
/// when different ops need different handling:
///
/// ```
/// use seqair::bam::pileup::PileupOp;
///
/// fn summarize(op: &PileupOp) -> &'static str {
///     match op {
///         PileupOp::Match { .. }        => "match",
///         PileupOp::Insertion { .. }    => "insertion",
///         PileupOp::Deletion { .. }     => "deletion",
///         PileupOp::ComplexIndel { .. } => "complex-indel",
///         PileupOp::RefSkip             => "ref-skip",
///     }
/// }
///
/// assert_eq!(summarize(&PileupOp::Match { qpos: 10, base: seqair_types::Base::A, qual: seqair_types::BaseQuality::from_byte(30) }), "match");
/// assert_eq!(summarize(&PileupOp::Deletion { del_len: 3 }), "deletion");
/// assert_eq!(summarize(&PileupOp::ComplexIndel { del_len: 3, insert_len: 2, is_refskip: false }), "complex-indel");
/// ```
// r[impl pileup_indel.op_enum]
// r[impl pileup_indel.type_safety]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PileupOp {
    /// Read has a base aligned at this position (M, =, or X CIGAR op).
    // r[impl types.base_quality.field_type]
    Match { qpos: u32, base: Base, qual: BaseQuality },
    /// Read has a base aligned at this position AND an insertion of
    /// `insert_len` query bases follows before the next reference position.
    /// Access inserted bases via the read's sequence at `qpos + 1 .. qpos + 1 + insert_len`.
    // r[impl types.base_quality.field_type]
    Insertion { qpos: u32, base: Base, qual: BaseQuality, insert_len: u32 },
    /// Read has a deletion spanning this position (D CIGAR op). `del_len` is the total length
    /// of the D CIGAR op — how many reference bases are deleted. No query base.
    Deletion { del_len: u32 },
    /// Deletion or ref-skip at this position with a following insertion
    /// (complex indel, e.g. D→I or N→I in CIGAR). Only emitted at the last
    /// position of the D/N op. No query base exists at this position.
    /// `is_refskip` distinguishes N→I (true) from D→I (false), matching
    /// htslib where `is_del=true, is_refskip=true/false, indel>0` coexist.
    // r[impl pileup_indel.complex_indel]
    ComplexIndel { del_len: u32, insert_len: u32, is_refskip: bool },
    /// Read has a reference skip at this position (N CIGAR op, e.g. intron). No query base.
    RefSkip,
}

const _: () = assert!(std::mem::size_of::<PileupOp>() <= 12, "PileupOp grew unexpectedly large");

// r[impl pileup.qpos]
// r[impl base_decode.alignment]
// r[impl pileup_indel.op_enum]
#[derive(Clone, Debug)]
pub struct PileupAlignment {
    record_idx: u32,
    pub op: PileupOp,
    pub mapq: u8,
    pub flags: BamFlags,
    pub strand: Strand,
    pub seq_len: u32,
    pub matching_bases: u32,
    pub indel_bases: u32,
}

// r[impl pileup_indel.accessors]
impl PileupAlignment {
    #[must_use]
    pub fn record_idx(&self) -> u32 {
        self.record_idx
    }

    #[must_use]
    pub fn op(&self) -> &PileupOp {
        &self.op
    }

    // r[impl pileup.qpos_none]
    #[must_use]
    pub fn qpos(&self) -> Option<usize> {
        match self.op {
            PileupOp::Match { qpos, .. } | PileupOp::Insertion { qpos, .. } => Some(qpos as usize),
            PileupOp::Deletion { .. } | PileupOp::ComplexIndel { .. } | PileupOp::RefSkip => None,
        }
    }

    #[must_use]
    pub fn base(&self) -> Option<Base> {
        match self.op {
            PileupOp::Match { base, .. } | PileupOp::Insertion { base, .. } => Some(base),
            PileupOp::Deletion { .. } | PileupOp::ComplexIndel { .. } | PileupOp::RefSkip => None,
        }
    }

    #[must_use]
    pub fn qual(&self) -> Option<BaseQuality> {
        match self.op {
            PileupOp::Match { qual, .. } | PileupOp::Insertion { qual, .. } => Some(qual),
            PileupOp::Deletion { .. } | PileupOp::ComplexIndel { .. } | PileupOp::RefSkip => None,
        }
    }

    // r[impl pileup_indel.accessors]
    #[must_use]
    pub fn is_del(&self) -> bool {
        matches!(self.op, PileupOp::Deletion { .. } | PileupOp::ComplexIndel { .. })
    }

    #[must_use]
    pub fn is_refskip(&self) -> bool {
        matches!(self.op, PileupOp::RefSkip | PileupOp::ComplexIndel { is_refskip: true, .. })
    }

    #[must_use]
    pub fn insert_len(&self) -> u32 {
        match self.op {
            PileupOp::Insertion { insert_len, .. } | PileupOp::ComplexIndel { insert_len, .. } => {
                insert_len
            }
            _ => 0,
        }
    }

    /// Returns the deletion length for a `Deletion` or `ComplexIndel` op, or 0 for all other ops.
    /// All positions within the same D CIGAR op report the same `del_len` (the total
    /// D op length, not the remaining bases in the deletion).
    // r[impl pileup_indel.accessors]
    #[must_use]
    pub fn del_len(&self) -> u32 {
        match self.op {
            PileupOp::Deletion { del_len } | PileupOp::ComplexIndel { del_len, .. } => del_len,
            _ => 0,
        }
    }
}

impl<U> PileupEngine<U> {
    /// Create a pileup engine that owns the record store.
    pub fn new(store: RecordStore<U>, region_start: Pos0, region_end: Pos0) -> Self {
        PileupEngine {
            store,
            current_pos: region_start,
            region_end,
            next_entry: 0,
            active_end_pos: Vec::new(),
            active: Vec::new(),
            max_depth: None,
            filter: None,
            ref_seq: None,
            columns_produced: 0,
            max_active_depth: 0,
        }
    }

    pub fn set_reference_seq(&mut self, ref_seq: RefSeq) {
        self.ref_seq = Some(ref_seq);
    }

    // r[impl pileup.max_depth]
    // r[impl pileup.max_depth_per_position]
    pub fn set_max_depth(&mut self, max: u32) {
        self.max_depth = Some(max);
    }

    // r[impl pileup.read_filter]
    /// Set a filter that receives `(flags, aux_bytes)` for each record entering the active set.
    pub fn set_filter(&mut self, f: impl Fn(BamFlags, &[u8]) -> bool + 'static) {
        self.filter = Some(Box::new(f));
    }

    pub fn remaining_positions(&self) -> usize {
        let diff = self
            .region_end
            .as_i64()
            .checked_sub(self.current_pos.as_i64())
            .and_then(|d| d.checked_add(1))
            .unwrap_or(0);
        #[expect(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            reason = "diff.max(0) is non-negative; bounded by region size which fits in usize on supported platforms"
        )]
        let r = diff.max(0) as usize;
        r
    }

    /// Borrow the underlying `RecordStore` for qname lookups during iteration.
    pub fn store(&self) -> &RecordStore<U> {
        &self.store
    }

    /// Take the `RecordStore` out for reuse. Returns `None` if already taken.
    ///
    /// Call this after iteration is complete. The returned store retains its
    /// allocated capacity but is cleared.
    pub fn take_store(&mut self) -> Option<RecordStore<U>> {
        if self.store.is_empty() && self.store.records_capacity() == 0 {
            return None;
        }
        let mut store = self.store.take_contents();
        store.clear();
        Some(store)
    }
}

impl<U> std::fmt::Debug for PileupEngine<U> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PileupEngine")
            .field("current_pos", &self.current_pos)
            .field("region_end", &self.region_end)
            .field("next_entry", &self.next_entry)
            .field("active_count", &self.active.len())
            .field("max_depth", &self.max_depth)
            .field("columns_produced", &self.columns_produced)
            .finish_non_exhaustive()
    }
}

// r[impl pileup.position_iteration]
impl<U> PileupEngine<U> {
    /// Core iteration logic used by [`pileups`](Self::pileups).
    ///
    /// Returns the per-column data (pos, reference base, alignments) so the caller
    /// can attach a store reference to build a [`PileupColumn`].
    fn advance(&mut self) -> Option<(Pos0, Base, Vec<PileupAlignment>)> {
        loop {
            if self.current_pos > self.region_end {
                return None;
            }

            let pos = self.current_pos;
            // current_pos <= region_end (checked above), and region_end < u32::MAX (niche),
            // so current_pos + 1 always fits.
            #[allow(
                clippy::unwrap_in_result,
                reason = "infallible: current_pos <= region_end < u32::MAX - 1"
            )]
            {
                self.current_pos = self
                    .current_pos
                    .checked_add_offset(Offset::new(1))
                    .trace_err("BUG: current_pos + 1 overflowed despite being <= region_end")?;
            }

            // Evict expired records. Iterate the compact end_pos vec (4-byte stride)
            // and swap-remove from both vecs in lockstep.
            {
                let mut i = 0;
                while i < self.active_end_pos.len() {
                    debug_assert!(i < self.active.len());
                    #[allow(
                        clippy::indexing_slicing,
                        reason = "i < active_end_pos.len() == active.len()"
                    )]
                    if self.active_end_pos[i] < pos {
                        self.active_end_pos.swap_remove(i);
                        self.active.swap_remove(i);
                    } else {
                        i = i.checked_add(1).trace_err("active set size exceeded usize::MAX")?;
                    }
                }
            }

            while self.next_entry < self.store.len() {
                #[expect(
                    clippy::cast_possible_truncation,
                    reason = "RecordStore capacity is bounded by SlabOverflow (u32); debug_assert enforces invariant"
                )]
                let idx = self.next_entry as u32;
                debug_assert_eq!(idx as usize, self.next_entry, "next_entry exceeds u32::MAX");

                let rec = self.store.record(idx);
                if rec.pos > pos {
                    break;
                }
                self.next_entry =
                    self.next_entry.checked_add(1).trace_err("next_entry exceeded usize::MAX")?;

                if rec.end_pos < pos {
                    continue;
                }

                // r[impl pileup.unmapped_excluded]
                if rec.flags.is_unmapped() {
                    continue;
                }

                if let Some(filter) = &self.filter
                    && !filter(rec.flags, self.store.aux(idx))
                {
                    continue;
                }

                let cigar = CigarMapping::new(rec.pos, self.store.cigar(idx))
                    .trace_err("failed to generate cigar mapping")?;

                self.active_end_pos.push(rec.end_pos);
                self.active.push(ActiveRecord {
                    record_idx: idx,
                    cigar,
                    flags: rec.flags,
                    strand: strand_from_flags(rec.flags),
                    mapq: rec.mapq,
                    seq_len: rec.seq_len,
                    matching_bases: rec.matching_bases,
                    indel_bases: rec.indel_bases,
                });
            }

            // r[impl pileup.empty_positions_skipped]
            if self.active.is_empty() {
                // r[impl pileup.trailing_empty_termination]
                if self.next_entry >= self.store.len() {
                    return None;
                }
                #[expect(
                    clippy::cast_possible_truncation,
                    reason = "RecordStore capacity is bounded by SlabOverflow (u32); debug_assert enforces invariant"
                )]
                let next_entry_u32 = self.next_entry as u32;
                debug_assert_eq!(
                    next_entry_u32 as usize, self.next_entry,
                    "next_entry exceeds u32::MAX"
                );
                let next_pos = self.store.record(next_entry_u32).pos;
                if next_pos > self.current_pos {
                    self.current_pos = next_pos;
                }
                continue;
            }

            // r[impl pileup_indel.op_enum]
            // r[impl pileup_indel.deletions_included]
            // r[impl pileup_indel.refskips_included]
            // r[related record_store.field_access]
            // r[impl perf.reuse_alignment_vec+2]
            let mut alignments = Vec::with_capacity(self.active.len());
            for active in &self.active {
                let Some(info) = active.cigar.pos_info_at(pos) else { continue };

                let op = match info {
                    CigarPosInfo::Match { qpos } => {
                        let qual = self.store.qual(active.record_idx);
                        let Some(&q) = qual.get(qpos as usize) else { continue };
                        PileupOp::Match {
                            qpos,
                            base: self.store.seq_at(active.record_idx, qpos as usize),
                            qual: q,
                        }
                    }
                    CigarPosInfo::Insertion { qpos, insert_len } => {
                        let qual = self.store.qual(active.record_idx);
                        let Some(&q) = qual.get(qpos as usize) else { continue };
                        PileupOp::Insertion {
                            qpos,
                            base: self.store.seq_at(active.record_idx, qpos as usize),
                            qual: q,
                            insert_len,
                        }
                    }
                    CigarPosInfo::Deletion { del_len } => PileupOp::Deletion { del_len },
                    // r[impl pileup_indel.complex_indel]
                    CigarPosInfo::ComplexIndel { del_len, insert_len, is_refskip } => {
                        PileupOp::ComplexIndel { del_len, insert_len, is_refskip }
                    }
                    CigarPosInfo::RefSkip => PileupOp::RefSkip,
                };

                alignments.push(PileupAlignment {
                    op,
                    mapq: active.mapq,
                    flags: active.flags,
                    strand: active.strand,
                    seq_len: active.seq_len,
                    matching_bases: active.matching_bases,
                    indel_bases: active.indel_bases,
                    record_idx: active.record_idx,
                });
            }

            // r[impl pileup.max_depth_per_position]
            if let Some(max) = self.max_depth {
                alignments.truncate(max as usize);
            }

            if !alignments.is_empty() {
                self.columns_produced = self.columns_produced.saturating_add(1);
                #[expect(
                    clippy::cast_possible_truncation,
                    reason = "depth is bounded by max_depth (u32) or typical pileup sizes; saturates at u32::MAX for profiling"
                )]
                let depth_u32 = alignments.len() as u32;
                self.max_active_depth = self.max_active_depth.max(depth_u32);
                let reference_base =
                    self.ref_seq.as_ref().map_or(Base::Unknown, |r| r.base_at(pos));
                return Some((pos, reference_base, alignments));
            }
        }
    }

    // r[impl pileup.lending_iterator]
    /// Advance to the next pileup column.
    ///
    /// Returns `Some(PileupColumn<'_, U>)` borrowing the record store. Call
    /// this in a `while let` loop. The column is valid until the next call to
    /// `pileups` on the same engine.
    ///
    /// This is a lending iterator — the returned column holds a borrow of the
    /// engine's store, so it cannot be collected into a `Vec` or held across
    /// subsequent `pileups` calls. Extract primitive data (pos, depth, etc.)
    /// if you need to retain it.
    pub fn pileups(&mut self) -> Option<PileupColumn<'_, U>> {
        let (pos, reference_base, alignments) = self.advance()?;
        Some(PileupColumn { pos, reference_base, alignments, store: &self.store })
    }

    /// Remaining positions in the current region — lower-bound estimate for
    /// pre-allocation of result vectors.
    pub fn size_hint(&self) -> (usize, Option<usize>) {
        (0, Some(self.remaining_positions()))
    }
}

impl<U> Drop for PileupEngine<U> {
    fn drop(&mut self) {
        if self.columns_produced > 0 {
            tracing::debug!(
                target: super::region_buf::PROFILE_TARGET,
                columns = self.columns_produced,
                max_depth = self.max_active_depth,
                active_cap = self.active.capacity().saturating_add(self.active_end_pos.capacity()),
                "pileup_engine",
            );
        }
    }
}

#[cfg(test)]
#[allow(clippy::arithmetic_side_effects, reason = "test arithmetic on known small values")]
mod tests {
    use super::*;

    #[test]
    fn ref_seq_base_at_within_range() {
        let ref_seq =
            RefSeq::new(Rc::from([Base::A, Base::C, Base::G, Base::T]), Pos0::new(100).unwrap());
        assert_eq!(ref_seq.base_at(Pos0::new(100).unwrap()), Base::A);
        assert_eq!(ref_seq.base_at(Pos0::new(101).unwrap()), Base::C);
        assert_eq!(ref_seq.base_at(Pos0::new(102).unwrap()), Base::G);
        assert_eq!(ref_seq.base_at(Pos0::new(103).unwrap()), Base::T);
    }

    #[test]
    fn ref_seq_base_at_before_start() {
        let ref_seq = RefSeq::new(Rc::from([Base::A, Base::C]), Pos0::new(100).unwrap());
        assert_eq!(ref_seq.base_at(Pos0::new(99).unwrap()), Base::Unknown);
        assert_eq!(ref_seq.base_at(Pos0::new(0).unwrap()), Base::Unknown);
    }

    #[test]
    fn ref_seq_base_at_after_end() {
        let ref_seq = RefSeq::new(Rc::from([Base::A, Base::C]), Pos0::new(100).unwrap());
        assert_eq!(ref_seq.base_at(Pos0::new(102).unwrap()), Base::Unknown);
        assert_eq!(ref_seq.base_at(Pos0::new(1000).unwrap()), Base::Unknown);
    }

    #[test]
    fn ref_seq_base_at_empty() {
        let ref_seq = RefSeq::new(Rc::from([]), Pos0::new(100).unwrap());
        assert_eq!(ref_seq.base_at(Pos0::new(100).unwrap()), Base::Unknown);
    }
}
