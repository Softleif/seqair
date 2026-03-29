//! Iterate over pileup columns. [`PileupEngine`] yields [`PileupColumn`]s with pre-extracted
//! flat fields per active read. Supports read filtering, per-position max-depth, and
//! opt-in overlapping-pair deduplication via [`PileupEngine::set_dedup_overlapping`].

use rustc_hash::FxHashMap;
use seqair_types::{Base, Strand, strand_from_flags};
// Rc is used only for RefSeq (reference sequence), not for BAM records.
// PileupEngine is intentionally !Send due to RecordFilter: Box<dyn Fn(...)>.
use std::rc::Rc;
use tracing::instrument;

use super::{
    cigar::CigarMapping,
    flags::{FLAG_SECOND_IN_TEMPLATE, FLAG_UNMAPPED},
    record_store::RecordStore,
};

pub struct RefSeq {
    bases: Rc<[Base]>,
    start_pos: i64,
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
    pub fn new(bases: Rc<[Base]>, start_pos: i64) -> Self {
        Self { bases, start_pos }
    }

    pub fn base_at(&self, pos: i64) -> Base {
        let offset = pos.wrapping_sub(self.start_pos);
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
    current_pos: i64,
    region_end: i64,
    next_entry: usize,
    /// Hot field: checked every column during retain. Stored separately so the
    /// retain loop strides 8 bytes instead of the full ActiveRecord size (~144 bytes).
    active_end_pos: Vec<i64>,
    /// Cold fields: only accessed for records that survive retain.
    active: Vec<ActiveRecord>,
    max_depth: Option<u32>,
    filter: Option<RecordFilter>,
    ref_seq: Option<RefSeq>,
    /// Map from record index to its mate's record index (if paired and both in store).
    mate_of: Option<Vec<Option<u32>>>,
    /// Reusable scratch buffers for dedup (avoids per-position allocation).
    dedup_to_remove: Vec<usize>,
    dedup_seen: FxHashMap<u32, usize>,
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
    pos: i64,
    reference_base: Base,
    alignments: Vec<PileupAlignment>,
}

impl PileupColumn {
    pub fn pos(&self) -> i64 {
        self.pos
    }

    pub fn depth(&self) -> usize {
        self.alignments.len()
    }

    pub fn alignments(&self) -> impl Iterator<Item = &PileupAlignment> {
        self.alignments.iter()
    }

    pub fn reference_base(&self) -> Base {
        self.reference_base
    }
}

// r[impl pileup.qpos]
// r[impl base_decode.alignment]
#[derive(Clone, Debug)]
pub struct PileupAlignment {
    record_idx: u32,
    qpos: usize,
    /// Pre-extracted fields for the hot loop — no store lookup needed.
    pub base: Base,
    pub qual: u8,
    pub mapq: u8,
    pub flags: u16,
    pub strand: Strand,
    pub seq_len: u32,
    pub matching_bases: u32,
    pub indel_bases: u32,
}

impl PileupAlignment {
    pub fn qpos(&self) -> usize {
        self.qpos
    }

    pub fn record_idx(&self) -> u32 {
        self.record_idx
    }
}

impl PileupEngine {
    /// Create a pileup engine that owns the record store.
    pub fn new(store: RecordStore, region_start: i64, region_end: i64) -> Self {
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
            mate_of: None,
            dedup_to_remove: Vec::new(),
            dedup_seen: FxHashMap::default(),
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

    // r[impl dedup.opt_in]
    // r[impl dedup.mate_detection]
    // r[impl dedup.mate_pairs_only]
    #[instrument(level = "debug", skip_all, fields(store_len = self.store.len()))]
    pub fn set_dedup_overlapping(&mut self) {
        use std::collections::hash_map::Entry;
        let n = self.store.len();
        let mut name_to_idx: FxHashMap<&[u8], u32> =
            FxHashMap::with_capacity_and_hasher(n / 2, Default::default());
        let mut mate_of = vec![None; n];

        for i in 0..n {
            let idx = i as u32;
            let qname = self.store.qname(idx);
            // r[impl cram.edge.unknown_read_names]
            if qname.is_empty() || qname == b"*" {
                continue;
            }
            match name_to_idx.entry(qname) {
                Entry::Vacant(e) => {
                    e.insert(idx);
                }
                Entry::Occupied(e) => {
                    let mate_idx = *e.get();
                    debug_assert!(
                        (mate_idx as usize) < mate_of.len(),
                        "mate_idx out of bounds: {} >= {}",
                        mate_idx,
                        mate_of.len()
                    );
                    debug_assert!(i < mate_of.len(), "i out of bounds: {i} >= {}", mate_of.len());
                    #[allow(clippy::indexing_slicing, reason = "both indices < n == mate_of.len()")]
                    if mate_of[mate_idx as usize].is_none() {
                        mate_of[i] = Some(mate_idx);
                        mate_of[mate_idx as usize] = Some(idx);
                    }
                }
            }
        }

        self.mate_of = Some(mate_of);
    }

    pub fn remaining_positions(&self) -> usize {
        (self.region_end - self.current_pos + 1).max(0) as usize
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
            .field("dedup_enabled", &self.mate_of.is_some())
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
            self.current_pos += 1;

            // Evict expired records. Iterate the compact end_pos vec (8-byte stride)
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

            // r[impl pileup.qpos_none]
            // r[related record_store.field_access]
            // r[impl perf.reuse_alignment_vec+2]
            let mut alignments = Vec::with_capacity(self.active.len());
            for active in &self.active {
                if let Some(qpos) = active.cigar.qpos_at(pos) {
                    let qual = self.store.qual(active.record_idx);
                    let Some(&q) = qual.get(qpos) else { continue };
                    alignments.push(PileupAlignment {
                        base: self.store.seq_at(active.record_idx, qpos),
                        qual: q,
                        mapq: active.mapq,
                        flags: active.flags,
                        strand: active.strand,
                        seq_len: active.seq_len,
                        matching_bases: active.matching_bases,
                        indel_bases: active.indel_bases,
                        record_idx: active.record_idx,
                        qpos,
                    });
                }
            }

            // r[impl dedup.per_position]
            // r[impl dedup.filter_independent]
            // r[impl dedup.max_depth_independent]
            if let Some(ref mate_of) = self.mate_of {
                dedup_overlapping_pairs(
                    &mut alignments,
                    mate_of,
                    &mut self.dedup_to_remove,
                    &mut self.dedup_seen,
                );
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
                dedup = self.mate_of.is_some(),
                dedup_remove_cap = self.dedup_to_remove.capacity(),
                dedup_seen_cap = self.dedup_seen.capacity(),
                "pileup_engine",
            );
        }
    }
}

// r[impl dedup.resolution_same_base]
// r[impl dedup.resolution_different_base]
fn dedup_overlapping_pairs(
    alignments: &mut Vec<PileupAlignment>,
    mate_of: &[Option<u32>],
    to_remove: &mut Vec<usize>,
    seen_record_idx: &mut FxHashMap<u32, usize>,
) {
    if alignments.len() <= 1 {
        return;
    }

    to_remove.clear();
    seen_record_idx.clear();

    for aln_idx in 0..alignments.len() {
        debug_assert!(
            aln_idx < alignments.len(),
            "aln_idx out of bounds: {aln_idx} >= {}",
            alignments.len()
        );
        #[allow(clippy::indexing_slicing, reason = "aln_idx < alignments.len() by loop bound")]
        let rec_idx = alignments[aln_idx].record_idx;

        #[allow(
            clippy::indexing_slicing,
            reason = "aln_idx/mate_aln_idx < alignments.len() by construction"
        )]
        if let Some(Some(mate_idx)) = mate_of.get(rec_idx as usize)
            && let Some(&mate_aln_idx) = seen_record_idx.get(mate_idx)
        {
            debug_assert!(
                mate_aln_idx < alignments.len(),
                "mate_aln_idx out of bounds: {mate_aln_idx} >= {}",
                alignments.len()
            );
            let this_base = alignments[aln_idx].base;
            let mate_base = alignments[mate_aln_idx].base;

            if this_base == mate_base || alignments[aln_idx].flags & FLAG_SECOND_IN_TEMPLATE != 0 {
                to_remove.push(aln_idx);
            } else {
                to_remove.push(mate_aln_idx);
            }
            continue;
        }

        seen_record_idx.insert(rec_idx, aln_idx);
    }

    // Reverse-order swap_remove is correct because:
    // - to_remove is sorted+deduped, so indices are unique and ascending
    // - swap_remove(i) moves the last element to position i
    // - iterating in reverse means we only swap_remove from positions >= current,
    //   so previously-removed indices are never invalidated
    to_remove.sort_unstable();
    to_remove.dedup();
    for &idx in to_remove.iter().rev() {
        alignments.swap_remove(idx);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ref_seq_base_at_within_range() {
        let ref_seq = RefSeq::new(Rc::from([Base::A, Base::C, Base::G, Base::T]), 100);
        assert_eq!(ref_seq.base_at(100), Base::A);
        assert_eq!(ref_seq.base_at(101), Base::C);
        assert_eq!(ref_seq.base_at(102), Base::G);
        assert_eq!(ref_seq.base_at(103), Base::T);
    }

    #[test]
    fn ref_seq_base_at_before_start() {
        let ref_seq = RefSeq::new(Rc::from([Base::A, Base::C]), 100);
        assert_eq!(ref_seq.base_at(99), Base::Unknown);
        assert_eq!(ref_seq.base_at(0), Base::Unknown);
        assert_eq!(ref_seq.base_at(-1), Base::Unknown);
    }

    #[test]
    fn ref_seq_base_at_after_end() {
        let ref_seq = RefSeq::new(Rc::from([Base::A, Base::C]), 100);
        assert_eq!(ref_seq.base_at(102), Base::Unknown);
        assert_eq!(ref_seq.base_at(1000), Base::Unknown);
    }

    #[test]
    fn ref_seq_base_at_empty() {
        let ref_seq = RefSeq::new(Rc::from([]), 100);
        assert_eq!(ref_seq.base_at(100), Base::Unknown);
    }
}
