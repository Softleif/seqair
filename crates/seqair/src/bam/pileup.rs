//! Iterate over pileup columns. [`PileupEngine`] yields [`PileupColumn`]s with pre-extracted
//! flat fields per active read. Supports per-position max-depth.
//!
//! Read filtering happens at fetch time via
//! [`CustomizeRecordStore::filter`](crate::bam::record_store::CustomizeRecordStore::filter),
//! not in this engine. Records that should not pollute the pileup never enter
//! the store in the first place.

use seqair_types::{BamFlags, Base, BaseQuality, Offset, Pos0, Strand, strand_from_flags};
// Rc is used only for RefSeq (reference sequence), not for BAM records.
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

    /// 0-based reference position of the first base in the loaded window.
    #[inline]
    pub fn start_pos(&self) -> Pos0 {
        self.start_pos
    }

    /// Number of bases loaded in this window.
    #[inline]
    pub fn len(&self) -> usize {
        self.bases.len()
    }

    /// `true` if no bases are loaded.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.bases.is_empty()
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

    /// Look up a single base at `pos`. Returns `None` if `pos` falls outside
    /// the loaded window. Use this when callers need to distinguish "outside
    /// the loaded reference window" from "inside the window but `Base::Unknown`
    /// (e.g. an `N` in the FASTA)" — [`base_at`](Self::base_at) collapses both
    /// into `Base::Unknown`.
    pub fn try_base_at(&self, pos: Pos0) -> Option<Base> {
        let offset = pos.as_i64().checked_sub(self.start_pos.as_i64())?;
        if offset < 0 {
            return None;
        }
        #[expect(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            reason = "offset is non-negative (checked above)"
        )]
        self.bases.get(offset as usize).copied()
    }

    /// Borrow `len` reference bases starting at `pos`. Returns `None` if any
    /// position in `pos..pos+len` falls outside the loaded window.
    ///
    /// Unlike [`base_at`](Self::base_at), partial overlaps return `None`
    /// rather than padding with `Base::Unknown`.
    pub fn range(&self, pos: Pos0, len: u32) -> Option<&[Base]> {
        let offset = pos.as_i64().checked_sub(self.start_pos.as_i64())?;
        if offset < 0 {
            return None;
        }
        #[expect(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            reason = "offset is non-negative (checked above)"
        )]
        let start = offset as usize;
        let end = start.checked_add(len as usize)?;
        self.bases.get(start..end)
    }
}

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

// r[impl pileup.empty_seq_unknown_base]
// Resolve the base + quality at `qpos` for a record on a Match/Insertion op.
// Records with `SEQ=*` (`seq_len == 0`) have an empty bases/qual slab slice;
// htslib's pileup keeps them in every covered column with base=N and the
// missing-quality sentinel. We mirror that here by returning Unknown +
// UNAVAILABLE so depth() and match_depth() stay in parity with htslib.
fn base_qual_at<U>(
    store: &RecordStore<U>,
    active: &ActiveRecord,
    qpos: u32,
) -> (Base, BaseQuality) {
    if active.seq_len == 0 {
        return (Base::Unknown, BaseQuality::UNAVAILABLE);
    }
    let qual = store.qual(active.record_idx);
    let q = qual.get(qpos as usize).copied().unwrap_or(BaseQuality::UNAVAILABLE);
    (store.seq_at(active.record_idx, qpos as usize), q)
}

// r[impl pileup.summary_base_counts]
/// Counts of canonical and unknown bases observed in a pileup column.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct BaseCounts {
    /// Count of `A` bases.
    pub a: usize,
    /// Count of `C` bases.
    pub c: usize,
    /// Count of `G` bases.
    pub g: usize,
    /// Count of `T` bases.
    pub t: usize,
    /// Count of unknown bases.
    pub unknown: usize,
}

impl BaseCounts {
    /// Add one observed base to the counts.
    pub fn add_base(&mut self, base: Base) {
        match base {
            Base::A => self.a = self.a.saturating_add(1),
            Base::C => self.c = self.c.saturating_add(1),
            Base::G => self.g = self.g.saturating_add(1),
            Base::T => self.t = self.t.saturating_add(1),
            Base::Unknown => self.unknown = self.unknown.saturating_add(1),
        }
    }
}

/// Counts of raw [`PileupOp`] variants observed in a pileup column.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PileupOpCounts {
    /// Count of [`PileupOp::Match`] operations.
    pub match_ops: usize,
    /// Count of [`PileupOp::Insertion`] operations.
    pub insertion_ops: usize,
    /// Count of [`PileupOp::Deletion`] operations.
    pub deletion_ops: usize,
    /// Count of [`PileupOp::ComplexIndel`] operations.
    pub complex_indel_ops: usize,
    /// Count of [`PileupOp::RefSkip`] operations.
    pub refskip_ops: usize,
}

// r[impl pileup.summary]
// r[impl pileup.summary_indel_counts]
/// Common counts for a single emitted pileup column.
///
/// Summary counts describe the same post-max-depth alignments exposed by
/// [`PileupColumn::raw_alignments`]. They do not apply downstream read filters
/// or base-quality thresholds.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PileupSummary {
    /// Total number of emitted alignments.
    pub depth: usize,
    /// Number of emitted alignments with a query base.
    pub match_depth: usize,
    /// Number of emitted alignments with an insertion signal.
    pub insertion_count: usize,
    /// Number of emitted alignments with a deletion signal.
    pub deletion_count: usize,
    /// Number of emitted alignments with a reference-skip signal.
    pub ref_skip_count: usize,
    /// Counts of raw pileup operation variants.
    pub op_counts: PileupOpCounts,
    /// Counts of bases from operations that carry a query base.
    pub base_counts: BaseCounts,
}

impl PileupSummary {
    /// Add one emitted alignment to the summary.
    pub fn observe(&mut self, alignment: &PileupAlignment) {
        self.depth = self.depth.saturating_add(1);
        match alignment.op {
            PileupOp::Match { base, .. } => {
                self.match_depth = self.match_depth.saturating_add(1);
                self.op_counts.match_ops = self.op_counts.match_ops.saturating_add(1);
                self.base_counts.add_base(base);
            }
            PileupOp::Insertion { base, .. } => {
                self.match_depth = self.match_depth.saturating_add(1);
                self.insertion_count = self.insertion_count.saturating_add(1);
                self.op_counts.insertion_ops = self.op_counts.insertion_ops.saturating_add(1);
                self.base_counts.add_base(base);
            }
            PileupOp::Deletion { .. } => {
                self.deletion_count = self.deletion_count.saturating_add(1);
                self.op_counts.deletion_ops = self.op_counts.deletion_ops.saturating_add(1);
            }
            PileupOp::ComplexIndel { insert_len, is_refskip, .. } => {
                self.op_counts.complex_indel_ops =
                    self.op_counts.complex_indel_ops.saturating_add(1);
                if insert_len > 0 {
                    self.insertion_count = self.insertion_count.saturating_add(1);
                }
                if is_refskip {
                    self.ref_skip_count = self.ref_skip_count.saturating_add(1);
                } else {
                    self.deletion_count = self.deletion_count.saturating_add(1);
                }
            }
            PileupOp::RefSkip => {
                self.ref_skip_count = self.ref_skip_count.saturating_add(1);
                self.op_counts.refskip_ops = self.op_counts.refskip_ops.saturating_add(1);
            }
        }
    }
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
    summary: PileupSummary,
    store: &'store RecordStore<U>,
}

impl<U> std::fmt::Debug for PileupColumn<'_, U> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PileupColumn")
            .field("pos", &self.pos)
            .field("reference_base", &self.reference_base)
            .field("depth", &self.summary.depth)
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
        self.summary.depth
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

    /// Common counts for this emitted column.
    // r[impl pileup.summary]
    #[must_use]
    pub fn summary(&self) -> &PileupSummary {
        &self.summary
    }

    /// Count of alignments with a query base at this position.
    ///
    /// Unlike [`depth`](Self::depth), deletions and ref-skips are not counted.
    /// Use `match_depth` when you need the number of reads that actually cover
    /// the position with a base (e.g., for allele frequency calculations).
    #[must_use]
    pub fn match_depth(&self) -> usize {
        self.summary.match_depth
    }

    /// Count of alignments with an insertion signal at this position.
    #[must_use]
    pub fn insertion_count(&self) -> usize {
        self.summary.insertion_count
    }

    /// Count of alignments with a deletion signal at this position.
    #[must_use]
    pub fn deletion_count(&self) -> usize {
        self.summary.deletion_count
    }

    /// Count of alignments with a reference-skip signal at this position.
    #[must_use]
    pub fn ref_skip_count(&self) -> usize {
        self.summary.ref_skip_count
    }

    /// Base counts for operations with a query base at this position.
    #[must_use]
    pub fn base_counts(&self) -> BaseCounts {
        self.summary.base_counts
    }

    /// Raw pileup operation counts for this position.
    #[must_use]
    pub fn op_counts(&self) -> PileupOpCounts {
        self.summary.op_counts
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
    /// The per-record extra, as computed by the customize value's `compute` method.
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

/// Metadata supplied when a custom pileup accumulator finishes a column.
pub struct PileupColumnContext<'store, U> {
    pos: Pos0,
    reference_base: Base,
    depth: usize,
    store: &'store RecordStore<U>,
}

impl<'store, U> PileupColumnContext<'store, U> {
    /// Reference position for the column.
    #[must_use]
    pub fn pos(&self) -> Pos0 {
        self.pos
    }

    /// Reference base for the column, or [`Base::Unknown`] when unavailable.
    pub fn reference_base(&self) -> Base {
        self.reference_base
    }

    /// Number of alignments observed by the accumulator.
    #[must_use]
    pub fn depth(&self) -> usize {
        self.depth
    }

    /// Borrow the record store for per-record extras, names, and auxiliary data.
    #[must_use]
    pub fn store(&self) -> &'store RecordStore<U> {
        self.store
    }
}

// r[impl pileup.custom_aggregation]
/// Accumulates caller-defined output while the pileup engine builds one column.
///
/// This is an opt-in alternative to [`PileupEngine::pileups`] for callers that
/// want to compute per-column summaries without materializing a public
/// [`PileupColumn`]. The accumulator observes the same post-max-depth
/// alignments that `pileups()` would expose.
pub trait PileupColumnAccumulator<U> {
    /// Owned output produced for each emitted column.
    type Output;

    /// Start a new emitted pileup column.
    fn begin_column(&mut self, pos: Pos0, reference_base: Base);

    /// Optional reservation hint for accumulators that store per-alignment data.
    fn reserve_column(&mut self, _active_upper_bound: usize) {}

    /// Observe one emitted alignment.
    fn observe_alignment(&mut self, alignment: PileupAlignment, store: &RecordStore<U>);

    /// Finish the current emitted column.
    ///
    /// Returning `None` skips user-visible output for this column; iteration will
    /// continue with the next column.
    fn finish_column(&mut self, context: PileupColumnContext<'_, U>) -> Option<Self::Output>;
}

enum AccumulatorAdvance<T> {
    End,
    Skip,
    Output(T),
}

struct MaterializedColumn {
    pos: Pos0,
    reference_base: Base,
    alignments: Vec<PileupAlignment>,
    summary: PileupSummary,
}

#[derive(Default)]
struct MaterializeColumnAccumulator {
    pos: Option<Pos0>,
    reference_base: Base,
    alignments: Vec<PileupAlignment>,
    summary: PileupSummary,
}

impl<U> PileupColumnAccumulator<U> for MaterializeColumnAccumulator {
    type Output = MaterializedColumn;

    fn begin_column(&mut self, pos: Pos0, reference_base: Base) {
        self.pos = Some(pos);
        self.reference_base = reference_base;
        self.alignments.clear();
        self.summary = PileupSummary::default();
    }

    fn reserve_column(&mut self, active_upper_bound: usize) {
        self.alignments.reserve(active_upper_bound);
    }

    fn observe_alignment(&mut self, alignment: PileupAlignment, _store: &RecordStore<U>) {
        self.summary.observe(&alignment);
        self.alignments.push(alignment);
    }

    fn finish_column(&mut self, _context: PileupColumnContext<'_, U>) -> Option<Self::Output> {
        let pos = self.pos.take()?;
        Some(MaterializedColumn {
            pos,
            reference_base: self.reference_base,
            alignments: std::mem::take(&mut self.alignments),
            summary: self.summary,
        })
    }
}

fn build_alignment<U>(
    store: &RecordStore<U>,
    active: &ActiveRecord,
    pos: Pos0,
) -> Option<PileupAlignment> {
    let info = active.cigar.pos_info_at(pos)?;

    let op = match info {
        CigarPosInfo::Match { qpos } => {
            let (base, qual) = base_qual_at(store, active, qpos);
            PileupOp::Match { qpos, base, qual }
        }
        CigarPosInfo::Insertion { qpos, insert_len } => {
            let (base, qual) = base_qual_at(store, active, qpos);
            PileupOp::Insertion { qpos, base, qual, insert_len }
        }
        CigarPosInfo::Deletion { del_len } => PileupOp::Deletion { del_len },
        // r[impl pileup_indel.complex_indel]
        CigarPosInfo::ComplexIndel { del_len, insert_len, is_refskip } => {
            PileupOp::ComplexIndel { del_len, insert_len, is_refskip }
        }
        CigarPosInfo::RefSkip => PileupOp::RefSkip,
    };

    Some(PileupAlignment {
        op,
        mapq: active.mapq,
        flags: active.flags,
        strand: active.strand,
        seq_len: active.seq_len,
        matching_bases: active.matching_bases,
        indel_bases: active.indel_bases,
        record_idx: active.record_idx,
    })
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

// r[impl pileup.extras.recover_store]
/// RAII guard around a [`PileupEngine`] that returns the underlying
/// [`RecordStore`] to a caller-provided slot when dropped.
///
/// Returned by [`Readers::pileup`](crate::reader::Readers::pileup). Derefs
/// to the inner [`PileupEngine`], so all engine methods (`pileups`,
/// `set_max_depth`, `store`, …) are available transparently. The point of
/// the guard is that buffer recovery is automatic: drop the guard (end of
/// scope, `break` out of the loop, `?`-propagated error mid-iteration)
/// and the populated-then-cleared store is moved back into the slot,
/// retaining its allocated capacity for the next pileup. Callers no
/// longer need an explicit `recover_store` step.
///
/// # Footgun: don't drain the engine via `Deref`
///
/// Because the guard `Deref`s to [`PileupEngine`], the engine's
/// [`PileupEngine::take_store`] is reachable as `guard.take_store()`.
/// Calling it leaves the engine's store empty, and then the guard's
/// `Drop` finds nothing to recover, so the next call to
/// [`Readers::pileup`](crate::reader::Readers::pileup) allocates a fresh
/// ~39 MB store. **Do not call `take_store` on a guard.** If you need
/// the underlying engine without recovery, use [`Self::into_inner`] —
/// that path is documented to disable recovery and at least makes the
/// intent explicit.
pub struct PileupGuard<'a, U = ()> {
    /// `Some` for the lifetime of the guard, then taken to `None` by
    /// [`Self::into_inner`] (so the `Drop` impl knows to skip recovery).
    /// Using an `Option` avoids the `unsafe` `ManuallyDrop` + `ptr::read`
    /// dance that an owning field would require.
    engine: Option<PileupEngine<U>>,
    slot: &'a mut RecordStore<U>,
}

impl<'a, U> PileupGuard<'a, U> {
    /// Build a guard from a populated engine and the slot to recover into.
    /// Used by [`Readers::pileup`] to wire up automatic store recovery;
    /// `pub(crate)` because external callers should always go through
    /// `Readers::pileup` to obtain one.
    pub(crate) fn new(engine: PileupEngine<U>, slot: &'a mut RecordStore<U>) -> Self {
        Self { engine: Some(engine), slot }
    }

    /// Consume the guard and return the inner engine, **without**
    /// recovering the store. Rare escape hatch — prefer letting the
    /// guard drop normally.
    ///
    /// Note that the originating `Readers`' internal store slot is
    /// already empty at this point (it was moved into the engine when
    /// the guard was constructed), and `into_inner` does not put
    /// anything back. The next call to
    /// [`Readers::pileup`](crate::reader::Readers::pileup) on that
    /// `Readers` will therefore allocate a fresh `RecordStore`. If you
    /// want the engine for inspection without disabling recovery,
    /// inspect through the guard via [`Deref`](std::ops::Deref) and
    /// drop the guard normally instead.
    pub fn into_inner(mut self) -> PileupEngine<U> {
        self.engine.take().expect("engine present until into_inner is called")
    }

    fn engine_ref(&self) -> &PileupEngine<U> {
        self.engine.as_ref().expect("engine present until guard is dropped or into_inner is called")
    }

    fn engine_mut(&mut self) -> &mut PileupEngine<U> {
        self.engine.as_mut().expect("engine present until guard is dropped or into_inner is called")
    }
}

impl<U> std::ops::Deref for PileupGuard<'_, U> {
    type Target = PileupEngine<U>;
    fn deref(&self) -> &PileupEngine<U> {
        self.engine_ref()
    }
}

impl<U> std::ops::DerefMut for PileupGuard<'_, U> {
    fn deref_mut(&mut self) -> &mut PileupEngine<U> {
        self.engine_mut()
    }
}

impl<U> std::fmt::Debug for PileupGuard<'_, U> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PileupGuard").field("engine", &self.engine).finish_non_exhaustive()
    }
}

impl<U> Drop for PileupGuard<'_, U> {
    fn drop(&mut self) {
        // Recover the store into the caller's slot. The engine is `None`
        // only after [`PileupGuard::into_inner`], in which case recovery
        // is intentionally skipped. `take_store` further returns `None`
        // if the store was already drained through `Deref` (the
        // documented footgun) or was empty with zero capacity from the
        // start. In any of those cases, leave the slot untouched.
        if let Some(mut engine) = self.engine.take()
            && let Some(store) = engine.take_store()
        {
            *self.slot = store;
        }
    }
}

// r[impl pileup.position_iteration]
impl<U> PileupEngine<U> {
    /// Core iteration logic used by [`pileups`](Self::pileups) and
    /// [`pileup_with`](Self::pileup_with).
    fn advance_with_accumulator<A>(&mut self, accumulator: &mut A) -> AccumulatorAdvance<A::Output>
    where
        A: PileupColumnAccumulator<U>,
    {
        loop {
            if self.current_pos > self.region_end {
                return AccumulatorAdvance::End;
            }

            let pos = self.current_pos;
            // current_pos <= region_end (checked above), and region_end < u32::MAX (niche),
            // so current_pos + 1 always fits.
            #[allow(
                clippy::unwrap_in_result,
                reason = "infallible: current_pos <= region_end < u32::MAX - 1"
            )]
            {
                let Some(next_pos) = self
                    .current_pos
                    .checked_add_offset(Offset::new(1))
                    .trace_err("BUG: current_pos + 1 overflowed despite being <= region_end")
                else {
                    return AccumulatorAdvance::End;
                };
                self.current_pos = next_pos;
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
                        let Some(next_i) =
                            i.checked_add(1).trace_err("active set size exceeded usize::MAX")
                        else {
                            return AccumulatorAdvance::End;
                        };
                        i = next_i;
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
                let Some(next_entry) =
                    self.next_entry.checked_add(1).trace_err("next_entry exceeded usize::MAX")
                else {
                    return AccumulatorAdvance::End;
                };
                self.next_entry = next_entry;

                if rec.end_pos < pos {
                    continue;
                }

                // r[impl pileup.unmapped_excluded]
                if rec.flags.is_unmapped() {
                    continue;
                }

                let Some(cigar) = CigarMapping::new(rec.pos, self.store.cigar(idx))
                    .trace_err("failed to generate cigar mapping")
                else {
                    return AccumulatorAdvance::End;
                };

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
                    return AccumulatorAdvance::End;
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
            // r[impl pileup.custom_aggregation_max_depth]
            let max_depth = self.max_depth.map(|max| usize::try_from(max).unwrap_or(usize::MAX));
            let mut emitted_depth = 0usize;
            let mut started = false;
            let reference_base = self.ref_seq.as_ref().map_or(Base::Unknown, |r| r.base_at(pos));

            for active in &self.active {
                if max_depth.is_some_and(|max| emitted_depth >= max) {
                    break;
                }

                let Some(alignment) = build_alignment(&self.store, active, pos) else { continue };

                if !started {
                    accumulator.begin_column(pos, reference_base);
                    accumulator.reserve_column(self.active.len());
                    started = true;
                }

                emitted_depth = emitted_depth.saturating_add(1);
                accumulator.observe_alignment(alignment, &self.store);
            }

            if emitted_depth > 0 {
                self.columns_produced = self.columns_produced.saturating_add(1);
                let depth_u32 = u32::try_from(emitted_depth).unwrap_or(u32::MAX);
                self.max_active_depth = self.max_active_depth.max(depth_u32);
                let context = PileupColumnContext {
                    pos,
                    reference_base,
                    depth: emitted_depth,
                    store: &self.store,
                };
                return match accumulator.finish_column(context) {
                    Some(output) => AccumulatorAdvance::Output(output),
                    None => AccumulatorAdvance::Skip,
                };
            }
        }
    }

    // r[impl pileup.custom_aggregation]
    /// Advance to the next accumulator-produced pileup output.
    ///
    /// The accumulator observes the same post-max-depth alignments that
    /// [`pileups`](Self::pileups) would expose, but it can summarize them
    /// without materializing a [`PileupColumn`]. If the accumulator returns
    /// `None` from [`PileupColumnAccumulator::finish_column`], iteration skips
    /// that column and continues to the next one.
    pub fn pileup_with<A>(&mut self, accumulator: &mut A) -> Option<A::Output>
    where
        A: PileupColumnAccumulator<U>,
    {
        loop {
            match self.advance_with_accumulator(accumulator) {
                AccumulatorAdvance::End => return None,
                AccumulatorAdvance::Skip => continue,
                AccumulatorAdvance::Output(output) => return Some(output),
            }
        }
    }

    // r[impl pileup.lending_iterator]
    // r[impl pileup.custom_aggregation_existing_api_unchanged]
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
        let mut accumulator = MaterializeColumnAccumulator::default();
        let materialized = self.pileup_with(&mut accumulator)?;
        Some(PileupColumn {
            pos: materialized.pos,
            reference_base: materialized.reference_base,
            alignments: materialized.alignments,
            summary: materialized.summary,
            store: &self.store,
        })
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
