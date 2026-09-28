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
//! Two more things, neither of them about lanes:
//!
//! - **Blocks, not tracks.** A column's five tracks and a cell's three
//!   matrices are each one block of lane-windows, so the column loop holds
//!   two base pointers; the row sweep that reads them is `lanes`, which the
//!   batch kernel runs too.
//! - **Derived once per call, gathered per group.** The rows depend only on
//!   the read and the columns only on the haplotype and the strand, and
//!   `align_reads` puts every read against every haplotype. So the plan
//!   derives each read's rows and each haplotype's columns per strand once
//!   per call, keyed by address (see `Recent`), as one window per row or
//!   column. A group only names each lane's tables: the kernel gathers the
//!   eight lanes' windows of a column or a row and transposes them into
//!   vectors, [`Lane::transpose`], rather than the plan scattering them into
//!   lane-interleaved tracks one `f32` at a time -- which was 17% of
//!   `align_reads` on the 3950X. The columns are gathered once per group
//!   into blocks, since every row of the band reads them; a row is
//!   transposed where the kernel reaches it, and never stored.
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
        Band, CODE_NO_CONVERSION, CODE_NO_PLAIN_MATCH, LANE_MAX, Lane, Window, Workspace, code,
    },
    emission::Emission,
    haplotype::Haplotype,
    lanes::{Cells, ColumnBlock, LaneRows, LanesView, RowEntry, lanes_kernel},
    read::Read,
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
/// | M4 Pro (NEON) | 0.22x | 0.44x | 0.65x | 0.84x | **1.05x** | 1.26x | 1.46x | 1.66x |
/// | 3950X (AVX2) | | | 0.83x | **1.03x** | 1.28x | 1.52x | | 2.03x |
///
/// The M4 crosses at five, the 3950X at four, and five is the number that is
/// safe on both -- it is also the lowest the bound below allows. Four would put
/// the M4 into the pairs kernel at 0.84x. `align_reads` only ever has one short
/// group per call -- the end of the product -- so on a locus of tens of reads
/// the choice moves a few per cent of one group; on rastair's real inputs
/// (notes §16.6) five against six measured within noise.
pub const PAIRS_BREAK_EVEN: usize = 5;

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

/// What one group's kernel writes: a read row of the three matrices, where
/// column `c` is `cells[c]` and `c = width + 1` is the pad the last row's `up`
/// reads, and the group's column blocks, gathered from its lanes' tables.
#[derive(Debug, Default)]
pub(crate) struct PairsBuffer {
    cells: Vec<Cells>,
    /// Shifted column `c` is `columns[c - first]` for the `first` the kernel
    /// gathered from; kept at its longest length, and only the gathered range
    /// is read.
    columns: Vec<ColumnBlock>,
}

/// Read base `index`'s row, as `banded::RowTracks::fill` derives it, in
/// [`RowEntry`]'s order.
#[allow(clippy::cast_possible_truncation, reason = "the f32 narrowing is the point")]
#[inline]
fn row_entry<E: Emission>(read: &Read, emission: &E, index: usize) -> Option<RowEntry> {
    let observation = read.observation(index)?;
    let eps = emission.epsilon(observation);
    let t = read.transition(index)?;
    Some([
        code(observation.base),
        ((1.0 - eps) - eps / 3.0) as f32,
        (eps / 3.0) as f32,
        t.match_to_match as f32,
        t.match_to_insertion as f32,
        t.match_to_deletion as f32,
        t.indel_to_match as f32,
        t.gap_continuation as f32,
    ])
}

/// One haplotype column's entry in every track, for one lane, in
/// [`ColumnBlock`]'s order and padded to a window: base, converted, plain,
/// rate, unconverted. The track a site does not use keeps the sentinel an
/// unvisited column has.
#[allow(clippy::cast_possible_truncation, reason = "the f32 narrowing is the point")]
#[inline]
fn column_entry<E: Emission>(
    haplotype: &Haplotype,
    emission: &E,
    strand: Strand,
    index: usize,
) -> Option<Window> {
    let weights = emission.site_weights(haplotype.site(index)?, strand);
    let base = code(weights.base);
    Some(match weights.converted {
        Some(converted) => [
            base,
            code(converted),
            CODE_NO_PLAIN_MATCH,
            weights.rate as f32,
            (1.0 - weights.rate) as f32,
            0.0,
            0.0,
            0.0,
        ],
        None => [base, CODE_NO_CONVERSION, base, 0.0, 0.0, 0.0, 0.0, 0.0],
    })
}

/// A column outside a lane's table: [`ColumnBlock::UNVISITED`]'s lane.
const NO_COLUMN: Window =
    [CODE_NO_PLAIN_MATCH, CODE_NO_CONVERSION, CODE_NO_PLAIN_MATCH, 0.0, 0.0, 0.0, 0.0, 0.0];

/// How many reads' rows a call keeps. `align_reads` packs a read's pairs
/// next to each other, so a read is reused within a group or two of where it
/// was derived.
const RECENT_READS: usize = 16;

/// How many haplotypes' columns a call keeps, per strand: the whole of a
/// locus's candidates up to 32 of them. A cache smaller than the set it cycles
/// through misses on every lookup -- 24 haplotypes through 16 slots did, on
/// the 10s dataset's largest group, and derived every column twice.
const RECENT_HAPLOTYPES: usize = 64;

// A group pins one slot per lane, so every cache needs a slot to spare, and a
// pin is a bit of a `u64`.
const _: () = assert!(RECENT_READS > PAIRS && RECENT_HAPLOTYPES > PAIRS);
const _: () = assert!(RECENT_READS <= 64 && RECENT_HAPLOTYPES <= 64);

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
///
/// A group's lanes refer to their tables by slot until its kernel has run, so
/// a slot handed out since [`Recent::unpin`] is never evicted: a group of
/// eight new reads would otherwise evict the table its first lane hit.
#[derive(Debug)]
struct Recent<K, T, const N: usize> {
    keys: Vec<Option<K>>,
    tables: Vec<Vec<T>>,
    /// The slot the next miss evicts, once every slot is in use.
    next: usize,
    /// The slots the current group uses, as bits.
    pinned: u64,
}

impl<K, T, const N: usize> Default for Recent<K, T, N> {
    fn default() -> Self {
        Self { keys: Vec::new(), tables: Vec::new(), next: 0, pinned: 0 }
    }
}

impl<K: Copy + PartialEq, T, const N: usize> Recent<K, T, N> {
    fn forget(&mut self) {
        self.keys.clear();
        self.next = 0;
        self.pinned = 0;
    }

    /// A new group: the previous one's slots may be evicted again.
    fn unpin(&mut self) {
        self.pinned = 0;
    }

    /// The slot holding the table for `key`, derived by `derive` into a
    /// cleared table if this call has not derived it yet, and pinned until
    /// the next [`Recent::unpin`]. A table whose derivation fails is
    /// forgotten rather than kept half-filled.
    #[inline]
    fn get_or_derive(
        &mut self,
        key: K,
        derive: impl FnOnce(&mut Vec<T>) -> Option<()>,
    ) -> Option<usize> {
        if let Some(slot) = self.keys.iter().position(|known| *known == Some(key)) {
            self.pinned |= 1 << slot;
            return Some(slot);
        }
        let slot = if self.keys.len() < N {
            self.keys.push(None);
            self.keys.len() - 1
        } else {
            // At most `PAIRS` slots are pinned and `N > PAIRS`, so this ends
            // within `PAIRS + 1` steps.
            while self.pinned & (1 << self.next) != 0 {
                self.next = (self.next + 1) % N;
            }
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
        self.pinned |= 1 << slot;
        Some(slot)
    }

    /// The table in `slot`; empty for a slot never filled.
    fn table(&self, slot: usize) -> &[T] {
        self.tables.get(slot).map_or(&[], Vec::as_slice)
    }
}

/// A haplotype whose band reaches only a small part of it is cheaper to derive
/// column by column than to derive whole and cache: past this many times the
/// lane's band span, a lane derives its own columns.
const WHOLE_HAPLOTYPE_LIMIT: usize = 4;

/// Where a lane's column entries come from.
#[derive(Debug, Default, Clone, Copy)]
enum ColumnTable {
    /// A lane with no pair.
    #[default]
    None,
    /// A haplotype's whole table, in [`PairsPlan::haplotype_columns`]'s slot.
    Recent(usize),
    /// Just the lane's band's columns, in its [`PairsPlan::own_columns`].
    Own,
}

/// Everything one group hoists out of its loops. The rows and columns are
/// not copied into the group: each lane names the tables it reads, and the
/// kernel gathers eight lanes' entries at a time into vectors.
#[derive(Debug, Default)]
pub(crate) struct PairsPlan {
    /// Each lane's rows, as a slot of `read_rows`; `None` for a lane that
    /// cannot align.
    rows_from: [Option<usize>; PAIRS],
    /// Each lane's columns, and the shifted column of the table's first
    /// entry.
    columns_from: [(ColumnTable, usize); PAIRS],
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
    read_rows: Recent<usize, Window, RECENT_READS>,
    /// Each haplotype's columns, per strand, per call.
    haplotype_columns: Recent<(usize, Strand), Window, RECENT_HAPLOTYPES>,
    /// The band's columns of a lane whose haplotype is too long to derive
    /// whole; see [`WHOLE_HAPLOTYPE_LIMIT`].
    own_columns: [Vec<Window>; PAIRS],
}

impl PairsPlan {
    /// Starts a call: nothing derived for an earlier call's pairs is valid.
    fn forget(&mut self) {
        self.read_rows.forget();
        self.haplotype_columns.forget();
    }

    /// Lane `lane`'s row table: its read's rows `1..=r` at `0..r`.
    fn rows(&self, lane: usize) -> &[Window] {
        match self.rows_from.get(lane).copied().flatten() {
            Some(slot) => self.read_rows.table(slot),
            None => &[],
        }
    }

    /// Lane `lane`'s column table and the shifted column of its first entry.
    fn columns(&self, lane: usize) -> (&[Window], usize) {
        match self.columns_from.get(lane).copied() {
            Some((ColumnTable::Recent(slot), first)) => (self.haplotype_columns.table(slot), first),
            Some((ColumnTable::Own, first)) => {
                (self.own_columns.get(lane).map_or(&[], Vec::as_slice), first)
            }
            Some((ColumnTable::None, _)) | None => (&[], 0),
        }
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
            rows_from,
            columns_from,
            init,
            front,
            past_end,
            lengths,
            read_rows,
            haplotype_columns,
            own_columns,
            ..
        } = self;
        *rows_from = [None; PAIRS];
        *columns_from = [(ColumnTable::None, 0); PAIRS];
        *init = [0.0; PAIRS];
        *front = [0.0; PAIRS];
        *past_end = [0.0; PAIRS];
        *lengths = [0; PAIRS];
        read_rows.unpin();
        haplotype_columns.unpin();

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

            *rows_from.get_mut(lane)? =
                Some(read_rows.get_or_derive(core::ptr::from_ref(read).addr(), |table| {
                    for index in 0..r {
                        table.push(row_entry(read, emission, index)?);
                    }
                    Some(())
                })?);

            // A lane's haplotype column `j`, zero-based, is shifted column
            // `j + 1 + delta`. A whole table serves every band; a band that
            // reaches only a sliver of a long haplotype derives just that.
            // The kernel reads a lane's columns past its band only on rows
            // past its read, whose zero rows compute zero from any finite
            // entry, so the two cannot score differently.
            let strand = read.strand();
            let (lo, hi) = pair.band.columns(h, r)?;
            *columns_from.get_mut(lane)? = if h > WHOLE_HAPLOTYPE_LIMIT * (hi - lo + 1) {
                let own = own_columns.get_mut(lane)?;
                own.clear();
                for index in lo..=hi {
                    own.push(column_entry(haplotype, emission, strand, index)?);
                }
                (ColumnTable::Own, lo + 1 + delta)
            } else {
                let key = (core::ptr::from_ref(haplotype).addr(), strand);
                let slot = haplotype_columns.get_or_derive(key, |table| {
                    for index in 0..h {
                        table.push(column_entry(haplotype, emission, strand, index)?);
                    }
                    Some(())
                })?;
                (ColumnTable::Recent(slot), 1 + delta)
            };
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

/// Gathers the group's column blocks from its lanes' tables, then sweeps it.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss,
    reason = "every count here is a few hundred"
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
    let (o, w) = (plan.offset, plan.half_width);
    let PairsBuffer { cells, columns } = buffer;

    // Each lane's tables, looked up once.
    let mut row_tables: [&[RowEntry]; PAIRS] = [&[]; PAIRS];
    let mut column_tables: [(&[Window], usize); PAIRS] = [(&[], 0); PAIRS];
    for (lane, (rows, columns)) in row_tables.iter_mut().zip(&mut column_tables).enumerate() {
        *rows = plan.rows(lane);
        *columns = plan.columns(lane);
    }

    // The columns any row's band reaches, `first..=last` over rows
    // `1..=read_len`, gathered eight lanes at a time: lane `k`'s entry for
    // shifted column `c` is its table's `c - first_k`, or the unvisited
    // sentinel outside it.
    let gathered_from = (1 + o - w).max(1);
    let gathered_to = (read_len as i64 + o + w).min(width as i64);
    let (origin, gathered) = if gathered_from <= gathered_to {
        (gathered_from as usize, (gathered_to - gathered_from + 1) as usize)
    } else {
        (1, 0)
    };
    if columns.len() < gathered {
        columns.resize(gathered, ColumnBlock::UNVISITED);
    }
    let Some(columns) = columns.get_mut(..gathered) else {
        return impossible;
    };
    for (column, block) in (origin..).zip(columns.iter_mut()) {
        let mut entries = [&NO_COLUMN; PAIRS];
        for (entry, &(table, first)) in entries.iter_mut().zip(&column_tables) {
            if let Some(found) = table.get(column.wrapping_sub(first)) {
                *entry = found;
            }
        }
        let [base, converted, plain, rate, unconverted, ..] = L::transpose(token, entries);
        base.store(&mut block.base);
        converted.store(&mut block.converted);
        plain.store(&mut block.plain);
        rate.store(&mut block.rate);
        unconverted.store(&mut block.unconverted);
    }

    let view = LanesView {
        columns,
        origin,
        init: &plan.init,
        front: &plan.front,
        past_end: &plan.past_end,
        lengths: &plan.lengths,
        rows_len: read_len,
        width,
        offset: o,
        half_width: w,
    };
    lanes_kernel::<L, _>(token, view, &LaneRows(row_tables), cells)
}
