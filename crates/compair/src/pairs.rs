//! One alignment per lane, any alignment: eight independent (read,
//! haplotype, band) pairs scored in lockstep, row by row.
//!
//! The batch kernel puts eight *haplotypes* against one read in a vector,
//! which is the right axis for a caller with many candidates per read. A
//! variant caller's shadow scoring has it the other way round: two to five
//! haplotypes and tens to hundreds of reads, each read with its own band.
//! Batching over haplotypes fills three lanes of eight there, below where the
//! batch kernel wins, so the whole workload went through the strip kernel one
//! pair at a time. This kernel lets lane `k` hold any pair at all, so the
//! lanes fill from the product of reads and haplotypes and a group is short
//! only at the very end of it.
//!
//! It is the batch kernel with three things made per lane that were per
//! batch:
//!
//! - **The band's offset.** Lane `k`'s columns are re-indexed as
//!   `c = j + delta_k`, with `delta_k = offset - offset_k` against the largest
//!   offset in the group. A cell is in lane `k`'s band when `|j - i -
//!   offset_k| <= half_width`, which is `|c - i - offset| <= half_width` for
//!   every lane: one band in `c`, so the column loop is still the band and has
//!   no per-lane range. The plan writes lane `k`'s haplotype column `j` at
//!   `c = j + delta_k`; a lane's haplotype then lives on `delta_k + 1 ..=
//!   delta_k + h_k`, and the half-width has to be the same for every lane,
//!   which is what a caller with one band width per locus has anyway.
//! - **The read.** Each of the eight row tracks holds one entry per lane, so
//!   a row is one vector load per track rather than a splat -- eight loads
//!   per row, against a column loop of ~47 steps. Rows
//!   past a lane's own read have all-zero tracks and compute to exactly zero,
//!   as the strip kernel's rows past the read do.
//! - **The last row.** The score is `m + i` summed along lane `k`'s own last
//!   row, so the sum is masked per lane, in the strip kernel's form `total +
//!   masked(m + i)`, and a lane's renormalisation stops after its last row:
//!   the strip kernel never renormalises below the read, and a shift after
//!   the total was taken would be applied to an exponent the total never saw.
//!
//! Clipping to the haplotype is two lane masks in principle, `j >= 1` and
//! `j <= h_k`. Only the second is in the column loop, as the batch kernel's
//! `past_end`. The first is needed only on row 0, the free start: a cell at
//! `j <= 0` reads its diagonal from `j - 1 < 0`, its insertion from the cell
//! above and its deletion from the cell to its left, and all three are zero
//! when row 0 is zero at `j < 0` -- so everything left of the haplotype
//! computes to exactly zero on every row without a mask, by induction along
//! the row and down the column.
//!
//! Two things differ from the batch kernel that are not about lanes:
//!
//! - **Blocks, not tracks.** A row's eight tracks, a column's five and a
//!   cell's three matrices are each one block of lane-windows rather than one
//!   interleaved `Vec` per track. The kernel reads the same vectors either way;
//!   what changes is that the column loop holds two base pointers instead of
//!   eight, which x86-64 has the registers for. With eight, the batch kernel's
//!   loop reloads the spilled ones every step: this loop is 50 instructions a
//!   step on AVX2 where the batch kernel's is 69.
//! - **Derived once per call.** The rows depend only on the read and the
//!   columns only on the haplotype and the strand, and `align_reads` puts
//!   every read against every haplotype. So the plan derives each read's rows
//!   and each haplotype's columns per strand once per call, keyed by address,
//!   and copies them into every lane that holds them; see `Recent`.
//!
//! The arithmetic, the order of the operations, the flush to zero, the
//! renormalisation cadence and the free start are the strip kernel's, cell for
//! cell, so every lane is bit-identical to [`align_strips`] on its own pair.
//! That is the gate; the parity tests hold every lane of mixed groups to it.
//!
//! As with the batch kernel, a lane is not free: a group of `n < PAIRS` pairs
//! computes all eight lanes, and the loop runs to the longest read in the
//! group and over the widest span of shifted columns. [`PAIRS_BREAK_EVEN`] is
//! where a short group stops being worth it.
//!
//! [`align_strips`]: crate::align_strips

use fearless_simd::Level;

use crate::{
    Strand,
    banded::{
        Band, CODE_NO_CONVERSION, CODE_NO_PLAIN_MATCH, ColumnLanes, LANE_MAX, Lane, RowLanes,
        TransitionLanes, Window, Workspace, code, prior,
    },
    emission::Emission,
    haplotype::Haplotype,
    read::Read,
    scaling::{exp2_f32, normalising_shift_f32},
    types::Log10Likelihood,
};

/// Pairs per group: the lane count.
pub const PAIRS: usize = LANE_MAX;

/// How many pairs a group needs before it beats scoring them one at a time
/// through the strip kernel.
///
/// Measured by `examples/pairsfill.rs` on the shadow fixture (150 bp reads,
/// each at its own offset, width 48), pairs kernel against strip kernel on the
/// same SIMD level, per alignment:
///
/// | pairs | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
/// |---|---|---|---|---|---|---|---|---|
/// | M4 Pro (NEON) | 0.19x | 0.38x | 0.55x | 0.70x | 0.93x | **0.96x** | 1.18x | 1.32x |
/// | 3950X (AVX2) | 0.30x | 0.57x | 0.83x | **0.98x** | 1.19x | 1.39x | 1.59x | 1.76x |
///
/// The machines disagree as they do for [`BATCH_BREAK_EVEN`], and for the same
/// reason six is the number: it costs the M4 ~4% on a group of six and the
/// 3950X ~1.4x on groups of four and five. `align_reads` only ever has one
/// short group per call -- the end of the product -- so on a locus of tens of
/// reads the choice moves a few per cent of one group.
///
/// [`BATCH_BREAK_EVEN`]: crate::BATCH_BREAK_EVEN
pub const PAIRS_BREAK_EVEN: usize = 6;

// The same bounds `BATCH_BREAK_EVEN` has, for the same reasons.
const _: () = assert!(
    PAIRS_BREAK_EVEN > PAIRS / 2 && PAIRS_BREAK_EVEN <= PAIRS,
    "the break-even has to be reachable and above half fill"
);

/// One alignment to score: a read against a haplotype, through a band.
#[derive(Debug, Clone, Copy)]
pub struct Pair<'a> {
    pub haplotype: &'a Haplotype,
    pub read: &'a Read,
    pub band: Band,
}

impl<'a> Pair<'a> {
    #[must_use]
    pub const fn new(haplotype: &'a Haplotype, read: &'a Read, band: Band) -> Self {
        Self { haplotype, read, band }
    }
}

/// One group scored, plus how many lanes that kernel fills; see
/// `BatchKernel` for why a function pointer.
#[derive(Clone, Copy)]
struct PairsKernel {
    lanes: usize,
    run: fn(&PairsPlan, &mut PairsBuffer) -> [Log10Likelihood; PAIRS],
}

impl PairsKernel {
    const fn over<L: Lane<Token = ()>>() -> Self {
        Self {
            lanes: if L::LANES < PAIRS { L::LANES } else { PAIRS },
            run: |plan, buffer| pairs_kernel::<L>((), plan, buffer),
        }
    }
}

/// The pairs kernel at the best SIMD level this CPU has.
const SIMD_PAIRS_KERNEL: PairsKernel = PairsKernel { lanes: PAIRS, run: pairs_kernel_simd };

fn pairs_kernel_simd(plan: &PairsPlan, buffer: &mut PairsBuffer) -> [Log10Likelihood; PAIRS] {
    crate::simd::pairs_kernel_at(Level::new(), plan, buffer)
}

/// Rows between two renormalisations: the strip kernel's constant.
const STRIP_ROWS: usize = LANE_MAX;

/// One read row of every lane: each track's entry for lane `k` at index `k`.
///
/// The tracks of a row are one block rather than eight interleaved `Vec`s so
/// that a row is one bounds check and one base pointer, in the fill and in
/// the kernel alike. Row 0 is the free start and holds zeros, and so does
/// every row past a lane's own read.
#[derive(Debug, Default, Clone, Copy)]
struct RowBlock {
    base: Window,
    /// `(1 - eps) - eps / 3`, as in `banded::RowTracks`.
    spread: Window,
    mismatched: Window,
    match_to_match: Window,
    match_to_insertion: Window,
    match_to_deletion: Window,
    indel_to_match: Window,
    gap_continuation: Window,
}

/// One shifted column of every lane: lane `k`'s haplotype column `c -
/// delta_k`, one-based, so `c = delta_k` is its free-start column. A block
/// for the same reason as [`RowBlock`], and it matters more here: the column
/// loop held eight track pointers, more than x86-64 has registers to spare,
/// and reloaded the spilled ones every step.
#[derive(Debug, Clone, Copy)]
struct ColumnBlock {
    base: Window,
    converted: Window,
    plain: Window,
    rate: Window,
    unconverted: Window,
}

impl ColumnBlock {
    /// What a column no lane's band reaches holds: never read, but finite.
    const UNVISITED: Self = Self {
        base: [CODE_NO_PLAIN_MATCH; PAIRS],
        converted: [CODE_NO_CONVERSION; PAIRS],
        plain: [CODE_NO_PLAIN_MATCH; PAIRS],
        rate: [0.0; PAIRS],
        unconverted: [0.0; PAIRS],
    };
}

/// The three matrices at one shifted column of one read row, every lane.
#[derive(Debug, Default, Clone, Copy)]
struct Cells {
    m: Window,
    i: Window,
    d: Window,
}

/// One read row of the three matrices; column `c` is `cells[c]`, and `c =
/// width + 1` is the pad the last row's `up` reads.
#[derive(Debug, Default)]
pub(crate) struct PairsBuffer {
    cells: Vec<Cells>,
}

/// One read row's entry in every track, for one lane: what a [`RowBlock`]
/// holds at index `k`.
#[derive(Debug, Default, Clone, Copy)]
struct RowEntry {
    base: f32,
    spread: f32,
    mismatched: f32,
    match_to_match: f32,
    match_to_insertion: f32,
    match_to_deletion: f32,
    indel_to_match: f32,
    gap_continuation: f32,
}

impl RowEntry {
    /// Read base `index`'s row, exactly as `banded::Plan::fill` derives it.
    #[allow(clippy::cast_possible_truncation, reason = "the f32 narrowing is the point")]
    #[inline]
    fn new<E: Emission>(read: &Read, emission: &E, index: usize) -> Option<Self> {
        let observation = read.observation(index)?;
        let eps = emission.epsilon(observation);
        let t = read.transition(index)?;
        Some(Self {
            base: code(observation.base),
            spread: ((1.0 - eps) - eps / 3.0) as f32,
            mismatched: (eps / 3.0) as f32,
            match_to_match: t.match_to_match as f32,
            match_to_insertion: t.match_to_insertion as f32,
            match_to_deletion: t.match_to_deletion as f32,
            indel_to_match: t.indel_to_match as f32,
            gap_continuation: t.gap_continuation as f32,
        })
    }

    /// Into lane `lane` of `row`. Inlined, and the lane taken modulo
    /// [`PAIRS`], so that the eight stores of a row carry no bounds checks:
    /// out of line, with them, this was 8% of a batch-shaped group's cycles on
    /// the 3950X.
    #[inline(always)]
    fn put(self, row: &mut RowBlock, lane: usize) -> Option<()> {
        debug_assert!(lane < PAIRS);
        let lane = lane % PAIRS;
        *row.base.get_mut(lane)? = self.base;
        *row.spread.get_mut(lane)? = self.spread;
        *row.mismatched.get_mut(lane)? = self.mismatched;
        *row.match_to_match.get_mut(lane)? = self.match_to_match;
        *row.match_to_insertion.get_mut(lane)? = self.match_to_insertion;
        *row.match_to_deletion.get_mut(lane)? = self.match_to_deletion;
        *row.indel_to_match.get_mut(lane)? = self.indel_to_match;
        *row.gap_continuation.get_mut(lane)? = self.gap_continuation;
        Some(())
    }
}

/// One haplotype column's entry in every track, for one lane: what a
/// [`ColumnBlock`] holds at index `k`.
#[derive(Debug, Clone, Copy)]
struct ColumnEntry {
    base: f32,
    converted: f32,
    plain: f32,
    rate: f32,
    unconverted: f32,
}

impl ColumnEntry {
    /// Haplotype site `index` as a read on `strand` sees it, exactly as
    /// `banded::Plan::fill` derives it; the track a site does not use keeps
    /// the sentinel an unvisited column has.
    #[allow(clippy::cast_possible_truncation, reason = "the f32 narrowing is the point")]
    #[inline]
    fn new<E: Emission>(
        haplotype: &Haplotype,
        emission: &E,
        strand: Strand,
        index: usize,
    ) -> Option<Self> {
        let weights = emission.site_weights(haplotype.site(index)?, strand);
        let base = code(weights.base);
        Some(match weights.converted {
            Some(converted) => Self {
                base,
                converted: code(converted),
                plain: CODE_NO_PLAIN_MATCH,
                rate: weights.rate as f32,
                unconverted: (1.0 - weights.rate) as f32,
            },
            None => Self {
                base,
                converted: CODE_NO_CONVERSION,
                plain: base,
                rate: 0.0,
                unconverted: 0.0,
            },
        })
    }

    /// Into lane `lane` of `column`. Inlined, and the lane taken modulo
    /// [`PAIRS`], so that the eight stores of a row carry no bounds checks:
    /// out of line, with them, this was 8% of a batch-shaped group's cycles on
    /// the 3950X.
    #[inline(always)]
    fn put(self, column: &mut ColumnBlock, lane: usize) -> Option<()> {
        debug_assert!(lane < PAIRS);
        let lane = lane % PAIRS;
        *column.base.get_mut(lane)? = self.base;
        *column.converted.get_mut(lane)? = self.converted;
        *column.plain.get_mut(lane)? = self.plain;
        *column.rate.get_mut(lane)? = self.rate;
        *column.unconverted.get_mut(lane)? = self.unconverted;
        Some(())
    }
}

/// How many reads' rows a call keeps. `align_reads` packs a read's pairs
/// next to each other, so a read is reused within a group or two of where it
/// was derived.
const RECENT_READS: usize = 16;

/// How many haplotypes' columns a call keeps, per strand: the whole of a
/// locus's candidates up to 32 of them. A cache smaller than the set it cycles
/// through misses on every lookup -- 24 haplotypes through 16 slots did, on
/// the 10s dataset's largest group, and derived every column twice.
const RECENT_HAPLOTYPES: usize = 64;

/// The tables a call has derived, keyed by the address of what they were
/// derived from, so that a read scored against three haplotypes has its rows
/// derived once rather than three times, and a haplotype its columns once per
/// strand rather than once per read.
///
/// An address is a sound key only while the thing it points at is borrowed,
/// which the pairs are for exactly one call -- so [`Recent::forget`] runs at
/// the start of every call, and a key never outlives the borrow it came from.
/// The emission is the call's too, and is not part of the key for the same
/// reason.
#[derive(Debug)]
struct Recent<K, T, const N: usize> {
    keys: Vec<Option<K>>,
    tables: Vec<Vec<T>>,
    /// The slot the next miss evicts, once every slot is in use.
    next: usize,
}

impl<K, T, const N: usize> Default for Recent<K, T, N> {
    fn default() -> Self {
        Self { keys: Vec::new(), tables: Vec::new(), next: 0 }
    }
}

impl<K: Copy + PartialEq, T, const N: usize> Recent<K, T, N> {
    fn forget(&mut self) {
        self.keys.clear();
        self.next = 0;
    }

    /// The table for `key`, derived by `derive` into a cleared table if this
    /// call has not derived it yet. A table whose derivation fails is
    /// forgotten rather than kept half-filled.
    #[inline]
    fn get_or_derive(
        &mut self,
        key: K,
        derive: impl FnOnce(&mut Vec<T>) -> Option<()>,
    ) -> Option<&[T]> {
        if let Some(slot) = self.keys.iter().position(|known| *known == Some(key)) {
            return self.tables.get(slot).map(Vec::as_slice);
        }
        let slot = if self.keys.len() < N {
            self.keys.push(None);
            self.keys.len() - 1
        } else {
            let slot = self.next;
            self.next = (slot + 1) % N;
            slot
        };
        if self.tables.len() <= slot {
            self.tables.resize_with(slot + 1, Vec::new);
        }
        let table = self.tables.get_mut(slot)?;
        table.clear();
        let known = self.keys.get_mut(slot)?;
        *known = None;
        derive(table)?;
        *known = Some(key);
        Some(table)
    }
}

/// A haplotype whose band reaches only a small part of it is cheaper to derive
/// column by column than to derive whole and cache: past this many times the
/// lane's band span, a lane derives its own columns.
const WHOLE_HAPLOTYPE_LIMIT: usize = 4;

/// Everything one group hoists out of its loops.
#[derive(Debug, Default)]
pub(crate) struct PairsPlan {
    rows: Vec<RowBlock>,
    columns: Vec<ColumnBlock>,
    /// `1 / h` per lane, the read's free start; zero for a lane that cannot
    /// align, which then scores [`Log10Likelihood::IMPOSSIBLE`].
    init: Window,
    /// `delta_k` per lane: row 0 is live from here.
    front: Window,
    /// `delta_k + h_k + 1` per lane, so `c < past_end[k]` is `j <= h_k`.
    past_end: Window,
    /// Each lane's read length; zero for a lane that cannot align.
    lengths: [usize; PAIRS],
    /// The longest read in the group.
    rows_len: usize,
    /// The last shifted column any lane has.
    width: usize,
    /// The group's band in shifted columns: the common half-width about the
    /// largest live offset.
    offset: i64,
    half_width: i64,
    /// Each read's rows, per call.
    read_rows: Recent<usize, RowEntry, RECENT_READS>,
    /// Each haplotype's columns, per strand, per call.
    haplotype_columns: Recent<(usize, Strand), ColumnEntry, RECENT_HAPLOTYPES>,
}

impl PairsPlan {
    /// Starts a call: nothing derived for an earlier call's pairs is valid.
    fn forget(&mut self) {
        self.read_rows.forget();
        self.haplotype_columns.forget();
    }

    /// `None` if the lanes do not share a half-width, which the grouping rules
    /// out, or if a read or haplotype cannot be addressed, which their
    /// constructors rule out.
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_precision_loss,
        reason = "the f32 narrowing is the point of this kernel, and a span is a few hundred columns"
    )]
    fn fill<E: Emission>(&mut self, group: &[Pair<'_>], emission: &E) -> Option<()> {
        let group = group.get(..group.len().min(PAIRS))?;
        // A lane whose band misses its haplotype scores impossible whatever
        // the kernel does, so it is left empty and does not get a say in the
        // group's offset: an arbitrary offset there would stretch every other
        // lane's shifted columns without bound.
        let mut live = [false; PAIRS];
        let mut offset = i64::MIN;
        let mut half_width = None;
        let mut rows_len = 0;
        for (slot, pair) in live.iter_mut().zip(group) {
            if half_width.is_some_and(|w| w != pair.band.half_width) {
                debug_assert!(false, "a group shares one half-width; `route` cuts it there");
                return None;
            }
            half_width = Some(pair.band.half_width);
            let (h, r) = (pair.haplotype.len(), pair.read.len());
            if pair.band.columns(h, r).is_none() {
                continue;
            }
            *slot = true;
            offset = offset.max(pair.band.offset);
            rows_len = rows_len.max(r);
        }
        let half_width = half_width?;

        // Every live offset is within `r + w` below and `h + w` above the
        // haplotype, which `Band::columns` just checked, so a shift is at most
        // a read and a haplotype and two half-widths.
        let mut width = 0;
        let mut deltas = [0usize; PAIRS];
        for ((delta, pair), &alive) in deltas.iter_mut().zip(group).zip(&live) {
            if alive {
                *delta = usize::try_from(offset.checked_sub(pair.band.offset)?).ok()?;
                width = width.max(*delta + pair.haplotype.len());
            }
        }

        let Self {
            rows,
            columns,
            init,
            front,
            past_end,
            lengths,
            read_rows,
            haplotype_columns,
            ..
        } = self;
        rows.clear();
        rows.resize(rows_len + 1, RowBlock::default());
        columns.clear();
        columns.resize(width + 2, ColumnBlock::UNVISITED);
        *init = [0.0; PAIRS];
        *front = [0.0; PAIRS];
        *past_end = [0.0; PAIRS];
        *lengths = [0; PAIRS];

        for (lane, ((pair, &alive), &delta)) in group.iter().zip(&live).zip(&deltas).enumerate() {
            if !alive {
                continue;
            }
            let (haplotype, read) = (pair.haplotype, pair.read);
            let (h, r) = (haplotype.len(), read.len());
            *init.get_mut(lane)? = 1.0 / h as f32;
            *front.get_mut(lane)? = delta as f32;
            *past_end.get_mut(lane)? = (delta + h + 1) as f32;
            *lengths.get_mut(lane)? = r;

            let entries = read_rows.get_or_derive(core::ptr::from_ref(read).addr(), |table| {
                for index in 0..r {
                    table.push(RowEntry::new(read, emission, index)?);
                }
                Some(())
            })?;
            for (entry, row) in entries.iter().zip(rows.get_mut(1..=r)?) {
                entry.put(row, lane)?;
            }

            // Only the columns the lane's band can reach, as `BatchPlan::fill`
            // does and for the same reason.
            let strand = read.strand();
            let (lo, hi) = pair.band.columns(h, r)?;
            let blocks = columns.get_mut(lo + 1 + delta..=hi + 1 + delta)?;
            if h > WHOLE_HAPLOTYPE_LIMIT * (hi - lo + 1) {
                for (index, column) in (lo..=hi).zip(blocks) {
                    ColumnEntry::new(haplotype, emission, strand, index)?.put(column, lane)?;
                }
            } else {
                let key = (core::ptr::from_ref(haplotype).addr(), strand);
                let entries = haplotype_columns.get_or_derive(key, |table| {
                    for index in 0..h {
                        table.push(ColumnEntry::new(haplotype, emission, strand, index)?);
                    }
                    Some(())
                })?;
                for (entry, column) in entries.get(lo..=hi)?.iter().zip(blocks) {
                    entry.put(column, lane)?;
                }
            }
        }

        // A group with no live lane has no rows, and its offset is never read.
        self.rows_len = rows_len;
        self.width = width;
        self.offset = if offset == i64::MIN { 0 } else { offset };
        self.half_width = half_width;
        Some(())
    }
}

impl Workspace {
    /// Every read against every haplotype, through whichever kernel is faster
    /// for the number of pairs. **This is the entry point for many reads
    /// against a few haplotypes**, which is a variant caller's shadow scoring;
    /// [`Workspace::align_candidates`] is the one for one read against many.
    ///
    /// `reads` carries each read's own band. The scores replace `out`'s
    /// contents, read-major: read `r` against haplotype `h` is
    /// `out[r * haplotypes.len() + h]`.
    ///
    /// The pairs are packed [`PAIRS`] at a time into the pairs kernel, in that
    /// order. A group that ends short of [`PAIRS_BREAK_EVEN`] -- the last one,
    /// or one cut short where the band width changes, since a group shares
    /// one -- goes through the strip kernel one pair at a time instead. The
    /// two kernels are bit-identical, so the split cannot change a score.
    pub fn align_reads<E: Emission>(
        &mut self,
        haplotypes: &[&Haplotype],
        reads: &[(&Read, Band)],
        emission: &E,
        out: &mut Vec<Log10Likelihood>,
    ) {
        let pairs = reads.iter().flat_map(|&(read, band)| {
            haplotypes.iter().map(move |&haplotype| Pair { haplotype, read, band })
        });
        out.clear();
        self.route(pairs, emission, out, PAIRS_BREAK_EVEN, SIMD_PAIRS_KERNEL.lanes, |this| {
            (SIMD_PAIRS_KERNEL.run)(&this.pairs_plan, &mut this.pairs_rows)
        });
    }

    /// Up to [`PAIRS`] pairs at a time, one pair per lane, whatever the fill.
    ///
    /// The scores replace `out`'s contents, in the pairs' order. Consecutive
    /// pairs share a group only while they share a band width. Prefer
    /// [`Workspace::align_reads`], which is this kernel where it wins and the
    /// strip kernel where it does not.
    pub fn align_pairs<E: Emission>(
        &mut self,
        pairs: &[Pair<'_>],
        emission: &E,
        out: &mut Vec<Log10Likelihood>,
    ) {
        out.clear();
        self.route(pairs.iter().copied(), emission, out, 0, SIMD_PAIRS_KERNEL.lanes, |this| {
            (SIMD_PAIRS_KERNEL.run)(&this.pairs_plan, &mut this.pairs_rows)
        });
    }

    /// [`Workspace::align_pairs`] at a given `fearless_simd` level, so the
    /// parity tests can hold every level the CPU has to the scalar kernel.
    #[doc(hidden)]
    pub fn align_pairs_at<E: Emission>(
        &mut self,
        level: Level,
        pairs: &[Pair<'_>],
        emission: &E,
        out: &mut Vec<Log10Likelihood>,
    ) {
        out.clear();
        self.route(pairs.iter().copied(), emission, out, 0, PAIRS, |this| {
            crate::simd::pairs_kernel_at(level, &this.pairs_plan, &mut this.pairs_rows)
        });
    }

    /// [`Workspace::align_pairs`] with one lane, the bit-parity oracle.
    pub fn align_pairs_scalar<E: Emission>(
        &mut self,
        pairs: &[Pair<'_>],
        emission: &E,
        out: &mut Vec<Log10Likelihood>,
    ) {
        let kernel = PairsKernel::over::<f32>();
        out.clear();
        self.route(pairs.iter().copied(), emission, out, 0, kernel.lanes, |this| {
            (kernel.run)(&this.pairs_plan, &mut this.pairs_rows)
        });
    }

    /// Cuts `pairs` into groups of at most `lanes` that share a half-width,
    /// and scores each group through `run` when it holds at least
    /// `break_even` pairs, or one pair at a time through the strip kernel
    /// otherwise. Appends to `out`.
    fn route<'a, E: Emission>(
        &mut self,
        pairs: impl Iterator<Item = Pair<'a>>,
        emission: &E,
        out: &mut Vec<Log10Likelihood>,
        break_even: usize,
        lanes: usize,
        mut run: impl FnMut(&mut Self) -> [Log10Likelihood; PAIRS],
    ) {
        self.pairs_plan.forget();
        let mut pairs = pairs.peekable();
        let Some(&first) = pairs.peek() else {
            return;
        };
        let lanes = lanes.clamp(1, PAIRS);
        let mut group = [first; PAIRS];
        let mut len = 0;
        for pair in pairs {
            let full = len == lanes;
            let other_width = group
                .first()
                .is_some_and(|head| len > 0 && head.band.half_width != pair.band.half_width);
            if full || other_width {
                self.group(
                    group.get(..len).unwrap_or_default(),
                    emission,
                    out,
                    break_even,
                    &mut run,
                );
                len = 0;
            }
            if let Some(slot) = group.get_mut(len) {
                *slot = pair;
                len += 1;
            }
        }
        self.group(group.get(..len).unwrap_or_default(), emission, out, break_even, &mut run);
    }

    /// One group, appended to `out`.
    fn group<E: Emission>(
        &mut self,
        group: &[Pair<'_>],
        emission: &E,
        out: &mut Vec<Log10Likelihood>,
        break_even: usize,
        run: &mut impl FnMut(&mut Self) -> [Log10Likelihood; PAIRS],
    ) {
        if group.is_empty() {
            return;
        }
        if group.len() < break_even {
            for pair in group {
                out.push(self.align_strips_simd(pair.haplotype, pair.read, emission, pair.band));
            }
            return;
        }
        let scores = if self.pairs_plan.fill(group, emission).is_some() {
            run(self)
        } else {
            [Log10Likelihood::IMPOSSIBLE; PAIRS]
        };
        out.extend(scores.iter().take(group.len()).copied());
    }
}

/// Up to [`PAIRS`] pairs at a time, one pair per lane.
pub fn align_pairs<E: Emission>(pairs: &[Pair<'_>], emission: &E) -> Vec<Log10Likelihood> {
    let mut out = Vec::with_capacity(pairs.len());
    Workspace::new().align_pairs(pairs, emission, &mut out);
    out
}

/// Every read against every haplotype, read-major, through whichever kernel
/// is faster for the number of pairs.
///
/// Allocates a [`Workspace`] per call; keep one and call
/// [`Workspace::align_reads`] when scoring more than one locus.
pub fn align_reads<E: Emission>(
    haplotypes: &[&Haplotype],
    reads: &[(&Read, Band)],
    emission: &E,
) -> Vec<Log10Likelihood> {
    let mut out = Vec::with_capacity(haplotypes.len() * reads.len());
    Workspace::new().align_reads(haplotypes, reads, emission, &mut out);
    out
}

/// The previous column's cells and the diagonal carried across a column step.
struct Carry<L> {
    left_m: L,
    left_d: L,
    diag_m: L,
    diag_indel: L,
}

/// Subnormals to zero, as in the strip kernel and for the same reason.
#[inline(always)]
fn flush<L: Lane>(value: L) -> L {
    L::masked_out(value.below(L::splat(value.token(), f32::MIN_POSITIVE)), value)
}

/// One row at a time along the shifted columns, every lane a different pair.
#[allow(
    clippy::too_many_lines,
    reason = "one traversal, and splitting it would hide the band algebra it exists to get right"
)]
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss,
    reason = "the f32 narrowing is the point of this kernel, and every count here is a few hundred"
)]
#[inline(always)]
pub(crate) fn pairs_kernel<L: Lane>(
    token: L::Token,
    plan: &PairsPlan,
    buffer: &mut PairsBuffer,
) -> [Log10Likelihood; PAIRS] {
    let impossible = [Log10Likelihood::IMPOSSIBLE; PAIRS];
    let (width, read_len) = (plan.width, plan.rows_len);
    if width == 0 || read_len == 0 {
        return impossible;
    }
    let span = width + 2;
    buffer.cells.clear();
    buffer.cells.resize(span, Cells::default());

    // Slices cut to one length, not the `Vec`s: see `batch_kernel`.
    let (Some(cells), Some(columns), Some(rows)) =
        (buffer.cells.get_mut(..span), plan.columns.get(..span), plan.rows.get(..=read_len))
    else {
        return impossible;
    };

    let zero = L::splat(token, 0.0);
    let one = L::splat(token, 1.0);
    let past_end = L::load(token, &plan.past_end);
    let init = L::load(token, &plan.init);
    let (o, w) = (plan.offset, plan.half_width);

    // Row 0, the free start: `1 / h_k` in the deletion matrix at every column
    // of the band that lane `k`'s haplotype reaches, which starts at its own
    // column 0, `delta_k`. Left of that, row 0 has to be zero: see the module
    // docs for why that alone keeps every row clear there.
    let front = L::load(token, &plan.front);
    let mut crossing = zero;
    for (column, slot) in cells.iter_mut().enumerate().take(width + 1) {
        if (column as i64 - o).abs() > w {
            continue;
        }
        let at = L::splat(token, column as f32);
        let cell = L::masked_out(at.below(front), L::masked(at.below(past_end), init));
        cell.store(&mut slot.d);
        crossing = crossing.vmax(cell);
    }

    let mut exponent = [0i32; PAIRS];
    let mut total = zero;
    let mut scratch = [0.0f32; LANE_MAX];
    let lanes = L::LANES.min(PAIRS);

    for (row, tracks) in rows.iter().enumerate().skip(1) {
        let first = (row as i64 + o - w).max(1);
        let last = (row as i64 + o + w).min(width as i64);
        if (row - 1) % STRIP_ROWS == 0 {
            crossing.store(&mut scratch);
            let mut lift = [1.0f32; LANE_MAX];
            let mut any = false;
            for lane in 0..lanes {
                // Only while the lane is still inside its read: the strip
                // kernel never renormalises below the last row, and the total
                // taken there is on the scale of the exponent at that row.
                let inside = plan.lengths.get(lane).is_some_and(|&r| row <= r);
                let shift = if inside {
                    scratch.get(lane).copied().map_or(0, normalising_shift_f32)
                } else {
                    0
                };
                any |= shift != 0;
                if let (Some(slot), Some(sum)) = (lift.get_mut(lane), exponent.get_mut(lane)) {
                    *slot = exp2_f32(shift);
                    *sum += shift;
                }
            }
            // Only the cells this row reads, `first - 1..=last`: left of them
            // nothing is read again, since the band only moves right, and
            // right of them every cell is still zero. Lifting the whole row
            // was ~9% of the batch kernel's cycles on the 3950X.
            let live = if first <= last {
                cells.get_mut(first as usize - 1..=last as usize)
            } else {
                None
            };
            if any && let Some(live) = live {
                let lift = L::load(token, &lift);
                for cell in live {
                    (L::load(token, &cell.m) * lift).store(&mut cell.m);
                    (L::load(token, &cell.i) * lift).store(&mut cell.i);
                    (L::load(token, &cell.d) * lift).store(&mut cell.d);
                }
            }
        }

        let lanes_row = RowLanes::new(
            L::load(token, &tracks.base),
            L::load(token, &tracks.spread),
            L::load(token, &tracks.mismatched),
        );
        let t = TransitionLanes {
            match_to_match: L::load(token, &tracks.match_to_match),
            match_to_insertion: L::load(token, &tracks.match_to_insertion),
            match_to_deletion: L::load(token, &tracks.match_to_deletion),
            indel_to_match: L::load(token, &tracks.indel_to_match),
            gap_continuation: L::load(token, &tracks.gap_continuation),
        };

        // The lanes whose read ends on this row, as a mask, built from the
        // integer lengths so that no `f32` compare of a row number is
        // involved; `None` on the rows where no lane ends.
        let summing = if plan.lengths.iter().take(lanes).any(|&r| r == row) {
            let mut ends = [0.0f32; LANE_MAX];
            for (slot, &r) in ends.iter_mut().zip(&plan.lengths) {
                *slot = if r == row { 1.0 } else { 0.0 };
            }
            Some(L::load(token, &ends).equals(one))
        } else {
            None
        };

        let mut running = zero;
        if first <= last {
            let (first, last) = (first as usize, last as usize);
            let (Some(before), Some(band), Some(band_columns)) =
                (cells.get(first - 1), cells.get(first..=last), columns.get(first..=last))
            else {
                return impossible;
            };
            debug_assert_eq!(band.len(), band_columns.len());
            let mut carry = Carry {
                left_m: zero,
                left_d: zero,
                diag_m: L::load(token, &before.m),
                diag_indel: L::load(token, &before.i) + L::load(token, &before.d),
            };

            // Carried and incremented rather than converted per step; see
            // `batch_kernel`.
            let mut column_lane = L::splat(token, first as f32);

            let Some(band) = cells.get_mut(first..=last) else {
                return impossible;
            };
            for (cell, column) in band.iter_mut().zip(band_columns) {
                let (up_m, up_i, up_d) =
                    (L::load(token, &cell.m), L::load(token, &cell.i), L::load(token, &cell.d));
                let prior_v = prior::<L>(
                    ColumnLanes {
                        base: L::load(token, &column.base),
                        converted: L::load(token, &column.converted),
                        plain: L::load(token, &column.plain),
                        rate: L::load(token, &column.rate),
                        unconverted_rate: L::load(token, &column.unconverted),
                    },
                    lanes_row,
                );

                let m = prior_v
                    * (carry.diag_m * t.match_to_match + carry.diag_indel * t.indel_to_match);
                let i = up_m * t.match_to_insertion + up_i * t.gap_continuation;
                let d = carry.left_m * t.match_to_deletion + carry.left_d * t.gap_continuation;
                let (m, i, d) = (flush(m), flush(i), flush(d));
                // Past a lane's haplotype, as in `batch_kernel`: `m` and `d`
                // masked, `i` zero by construction.
                let keep = column_lane.below(past_end);
                column_lane = column_lane + one;
                let (m, d) = (L::masked(keep, m), L::masked(keep, d));

                m.store(&mut cell.m);
                i.store(&mut cell.i);
                d.store(&mut cell.d);

                running = running.vmax(m).vmax(i).vmax(d);
                if let Some(ends) = summing {
                    // The strip kernel's form, parentheses and mask included:
                    // `f32` addition is not associative, and a lane that does
                    // not end here adds an exact zero.
                    total = total + L::masked(ends, m + i);
                }
                carry.left_m = m;
                carry.left_d = d;
                carry.diag_m = up_m;
                carry.diag_indel = up_i + up_d;
            }

            // The two cells the next row reads that this one did not write;
            // see `batch_kernel`.
            for at in [first - 1, last + 1] {
                if let Some(cell) = cells.get_mut(at) {
                    *cell = Cells::default();
                }
            }
        }
        if row % STRIP_ROWS == 0 {
            crossing = running;
        }
    }

    total.store(&mut scratch);
    let mut out = impossible;
    for lane in 0..lanes {
        let (Some(sum), Some(shift)) = (scratch.get(lane).copied(), exponent.get(lane).copied())
        else {
            continue;
        };
        if sum > 0.0
            && let Some(slot) = out.get_mut(lane)
        {
            *slot = Log10Likelihood::new(
                f64::from(sum).log10() - f64::from(shift) * core::f64::consts::LOG10_2,
            );
        }
    }
    out
}
