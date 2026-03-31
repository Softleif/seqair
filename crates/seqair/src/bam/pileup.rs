//! Iterate over pileup columns. [`PileupEngine`] yields [`PileupColumn`]s with pre-extracted
//! flat fields per active read. Supports read filtering and per-position max-depth.

use seqair_types::{Base, Offset, Pos, Strand, Zero, strand_from_flags};
// Rc is used only for RefSeq (reference sequence), not for BAM records.
// PileupEngine is intentionally !Send due to RecordFilter: Box<dyn Fn(...)>.
use std::rc::Rc;

use super::{
    cigar::{CigarMapping, CigarPosInfo},
    flags::FLAG_UNMAPPED,
    record_store::RecordStore,
};

pub struct RefSeq {
    bases: Rc<[Base]>,
    start_pos: Pos<Zero>,
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
    pub fn new(bases: Rc<[Base]>, start_pos: Pos<Zero>) -> Self {
        Self { bases, start_pos }
    }

    pub fn base_at(&self, pos: Pos<Zero>) -> Base {
        let offset = pos.as_i64() - self.start_pos.as_i64();
        if offset < 0 {
            return Base::Unknown;
        }
        self.bases.get(offset as usize).copied().unwrap_or(Base::Unknown)
    }
}

/// Filter function receiving (flags, aux_bytes) for the candidate record.
type RecordFilter = Box<dyn Fn(u16, &[u8]) -> bool>;

// r[impl pileup.active_set]
// r[impl pileup.zero_refspan_reads]
// r[impl pileup.soft_clip_at_position]
// r[impl perf.no_sorted_indices]
// r[impl perf.avoid_redundant_arena_get+2]
// r[impl perf.cigar_no_to_vec]
// r[impl perf.reuse_alignment_vec+2]
pub struct PileupEngine {
    store: RecordStore,
    current_pos: Pos<Zero>,
    region_end: Pos<Zero>,
    next_entry: usize,
    /// Hot field: checked every column during retain. Stored separately so the
    /// retain loop strides 4 bytes instead of the full ActiveRecord size (~144 bytes).
    active_end_pos: Vec<Pos<Zero>>,
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
    flags: u16,
    strand: Strand,
    mapq: u8,
    seq_len: u32,
    matching_bases: u32,
    indel_bases: u32,
}

// r[impl pileup.column_contents]
// r[impl pileup.htslib_compat]
#[derive(Debug)]
pub struct PileupColumn {
    pos: Pos<Zero>,
    reference_base: Base,
    alignments: Vec<PileupAlignment>,
}

impl PileupColumn {
    #[must_use]
    pub fn pos(&self) -> Pos<Zero> {
        self.pos
    }

    #[must_use]
    pub fn depth(&self) -> usize {
        self.alignments.len()
    }

    pub fn alignments(&self) -> impl Iterator<Item = &PileupAlignment> {
        self.alignments.iter()
    }

    pub fn reference_base(&self) -> Base {
        self.reference_base
    }

    /// Count of alignments with a base at this position (excludes deletions and ref-skips).
    #[must_use]
    pub fn match_depth(&self) -> usize {
        self.alignments.iter().filter(|a| a.qpos().is_some()).count()
    }
}

/// What a read shows at a specific reference position in the pileup.
///
/// `base`, `qual`, and `qpos` are only available for `Match` and `Insertion`
/// variants — the type system prevents reading a base from a deletion.
// r[impl pileup_indel.op_enum]
// r[impl pileup_indel.type_safety]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PileupOp {
    /// Read has a base aligned at this position (M, =, or X CIGAR op).
    Match { qpos: u32, base: Base, qual: u8 },
    /// Read has a base aligned at this position AND an insertion of
    /// `insert_len` query bases follows before the next reference position.
    /// Access inserted bases via the read's sequence at `qpos + 1 .. qpos + 1 + insert_len`.
    Insertion { qpos: u32, base: Base, qual: u8, insert_len: u32 },
    /// Read has a deletion spanning this position (D CIGAR op). `del_len` is the total length
    /// of the D CIGAR op — how many reference bases are deleted. No query base.
    Deletion { del_len: u32 },
    /// Read has a reference skip at this position (N CIGAR op, e.g. intron). No query base.
    RefSkip,
}

const _: () = assert!(std::mem::size_of::<PileupOp>() <= 16, "PileupOp grew unexpectedly large");

// r[impl pileup.qpos]
// r[impl base_decode.alignment]
// r[impl pileup_indel.op_enum]
#[derive(Clone, Debug)]
pub struct PileupAlignment {
    record_idx: u32,
    pub op: PileupOp,
    pub mapq: u8,
    pub flags: u16,
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

    #[must_use]
    pub fn qpos(&self) -> Option<usize> {
        match self.op {
            PileupOp::Match { qpos, .. } | PileupOp::Insertion { qpos, .. } => Some(qpos as usize),
            PileupOp::Deletion { .. } | PileupOp::RefSkip => None,
        }
    }

    #[must_use]
    pub fn base(&self) -> Option<Base> {
        match self.op {
            PileupOp::Match { base, .. } | PileupOp::Insertion { base, .. } => Some(base),
            PileupOp::Deletion { .. } | PileupOp::RefSkip => None,
        }
    }

    #[must_use]
    pub fn qual(&self) -> Option<u8> {
        match self.op {
            PileupOp::Match { qual, .. } | PileupOp::Insertion { qual, .. } => Some(qual),
            PileupOp::Deletion { .. } | PileupOp::RefSkip => None,
        }
    }

    #[must_use]
    pub fn is_del(&self) -> bool {
        matches!(self.op, PileupOp::Deletion { .. })
    }

    #[must_use]
    pub fn is_refskip(&self) -> bool {
        matches!(self.op, PileupOp::RefSkip)
    }

    #[must_use]
    pub fn insert_len(&self) -> u32 {
        match self.op {
            PileupOp::Insertion { insert_len, .. } => insert_len,
            _ => 0,
        }
    }

    /// Returns the deletion length for a `Deletion` op, or 0 for all other ops.
    /// All positions within the same D CIGAR op report the same `del_len` (the total
    /// D op length, not the remaining bases in the deletion).
    // r[impl pileup_indel.accessors]
    #[must_use]
    pub fn del_len(&self) -> u32 {
        match self.op {
            PileupOp::Deletion { del_len } => del_len,
            _ => 0,
        }
    }
}

impl PileupEngine {
    /// Create a pileup engine that owns the record store.
    pub fn new(store: RecordStore, region_start: Pos<Zero>, region_end: Pos<Zero>) -> Self {
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
    pub fn set_filter(&mut self, f: impl Fn(u16, &[u8]) -> bool + 'static) {
        self.filter = Some(Box::new(f));
    }

    pub fn remaining_positions(&self) -> usize {
        let diff = self.region_end.as_i64() - self.current_pos.as_i64() + 1;
        diff.max(0) as usize
    }

    /// Borrow the underlying `RecordStore` for qname lookups during iteration.
    pub fn store(&self) -> &RecordStore {
        &self.store
    }

    /// Take the `RecordStore` out for reuse. Returns `None` if already taken.
    ///
    /// Call this after iteration is complete. The returned store retains its
    /// allocated capacity but is cleared.
    pub fn take_store(&mut self) -> Option<RecordStore> {
        if self.store.is_empty() && self.store.records_capacity() == 0 {
            return None;
        }
        let mut store = std::mem::take(&mut self.store);
        store.clear();
        Some(store)
    }
}

impl std::fmt::Debug for PileupEngine {
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
impl Iterator for PileupEngine {
    type Item = PileupColumn;

    fn size_hint(&self) -> (usize, Option<usize>) {
        (0, Some(self.remaining_positions()))
    }

    fn next(&mut self) -> Option<PileupColumn> {
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
                    .expect("BUG: current_pos + 1 overflowed despite being <= region_end");
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
                        i += 1;
                    }
                }
            }

            while self.next_entry < self.store.len() {
                let idx = self.next_entry as u32;
                debug_assert_eq!(idx as usize, self.next_entry, "next_entry exceeds u32::MAX");

                let rec = self.store.record(idx);
                if rec.pos > pos {
                    break;
                }
                self.next_entry += 1;

                if rec.end_pos < pos {
                    continue;
                }

                // r[impl pileup.unmapped_excluded]
                if rec.flags & FLAG_UNMAPPED != 0 {
                    continue;
                }

                if let Some(filter) = &self.filter
                    && !filter(rec.flags, self.store.aux(idx))
                {
                    continue;
                }

                let cigar = CigarMapping::new(rec.pos, self.store.cigar(idx));

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
                self.max_active_depth = self.max_active_depth.max(alignments.len() as u32);
                let reference_base =
                    self.ref_seq.as_ref().map_or(Base::Unknown, |r| r.base_at(pos));
                return Some(PileupColumn { pos, reference_base, alignments });
            }
        }
    }
}

impl Drop for PileupEngine {
    fn drop(&mut self) {
        if self.columns_produced > 0 {
            tracing::debug!(
                target: super::region_buf::PROFILE_TARGET,
                columns = self.columns_produced,
                max_depth = self.max_active_depth,
                active_cap = self.active.capacity() + self.active_end_pos.capacity(),
                "pileup_engine",
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ref_seq_base_at_within_range() {
        let ref_seq = RefSeq::new(
            Rc::from([Base::A, Base::C, Base::G, Base::T]),
            Pos::<Zero>::new(100).unwrap(),
        );
        assert_eq!(ref_seq.base_at(Pos::<Zero>::new(100).unwrap()), Base::A);
        assert_eq!(ref_seq.base_at(Pos::<Zero>::new(101).unwrap()), Base::C);
        assert_eq!(ref_seq.base_at(Pos::<Zero>::new(102).unwrap()), Base::G);
        assert_eq!(ref_seq.base_at(Pos::<Zero>::new(103).unwrap()), Base::T);
    }

    #[test]
    fn ref_seq_base_at_before_start() {
        let ref_seq = RefSeq::new(Rc::from([Base::A, Base::C]), Pos::<Zero>::new(100).unwrap());
        assert_eq!(ref_seq.base_at(Pos::<Zero>::new(99).unwrap()), Base::Unknown);
        assert_eq!(ref_seq.base_at(Pos::<Zero>::new(0).unwrap()), Base::Unknown);
    }

    #[test]
    fn ref_seq_base_at_after_end() {
        let ref_seq = RefSeq::new(Rc::from([Base::A, Base::C]), Pos::<Zero>::new(100).unwrap());
        assert_eq!(ref_seq.base_at(Pos::<Zero>::new(102).unwrap()), Base::Unknown);
        assert_eq!(ref_seq.base_at(Pos::<Zero>::new(1000).unwrap()), Base::Unknown);
    }

    #[test]
    fn ref_seq_base_at_empty() {
        let ref_seq = RefSeq::new(Rc::from([]), Pos::<Zero>::new(100).unwrap());
        assert_eq!(ref_seq.base_at(Pos::<Zero>::new(100).unwrap()), Base::Unknown);
    }
}
