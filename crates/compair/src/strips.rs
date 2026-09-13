//! The strip traversal: one read row per lane, swept along the haplotype.
//!
//! The diagonal kernel in `banded` walks whole anti-diagonals, so every chunk
//! of eight cells reloads its three neighbours and its eight per-row inputs
//! from memory. This kernel is the traversal of Intel's Genomics Kernel
//! Library: the read is cut into strips of one row per lane and each strip is
//! swept along the haplotype, lane `l` holding cell `(r0 + l, d - l)` at step
//! `d` -- an eight-cell anti-diagonal segment that moves one column per step.
//! What that buys:
//!
//! - everything per row, the read base, the emission's two terms and the five
//!   transitions, is loaded once per strip and stays in registers for the
//!   whole sweep;
//! - the three neighbours of a cell are the previous two steps' vectors, still
//!   in registers. `left` is the previous step's own lane. `up` is the previous
//!   step shifted up by one lane, with the row above the strip entering at
//!   lane 0. `diag` is the step before that, shifted the same way -- which is
//!   exactly the previous step's `up`, so it costs nothing. A step shifts three
//!   vectors and reads three scalars from a row buffer where a diagonal chunk
//!   loaded seven vectors;
//! - the haplotype columns are still one contiguous load per track, because
//!   lane `l` reads column `d - l` and the plan holds the columns reversed.
//!
//! What it costs: the band's edges cross a strip diagonally, so a sweep begins
//! and ends with steps that are only partly inside the band and are masked,
//! and the eight rows of a strip are eight rows of dynamic range between two
//! renormalisations rather than a diagonal's one.
//!
//! Renormalisation happens between strips, on the row that crosses from one to
//! the next: the row buffer is multiplied by the power of two that puts that
//! row's maximum into `[1, 2)` before the next strip reads it, and since the
//! row buffer is the only way any value enters a strip, that rescales
//! everything the strip will compute. It is a power of two, so it is exact,
//! and it is applied every eight rows whatever the lane count -- the scalar
//! instance's strip is one row, and renormalising every row instead would put
//! the flush-to-zero threshold at a different place. That, and every operation
//! on a cell being lanewise, is what makes the scalar and the eight-lane
//! instances bit-identical by construction, as the diagonal kernels are.

use wide::f32x8;

use crate::{
    banded::{
        Band, COLUMN_FRONT, ColumnLanes, LANE_MAX, Lane, Plan, RowLanes, Shape, TransitionLanes,
        View, Window, Workspace, prior, reset,
    },
    emission::Emission,
    haplotype::Haplotype,
    read::Read,
    scaling::{exp2_f32, normalising_shift_f32},
    types::Log10Likelihood,
};

/// Rows between two renormalisations, whatever the lane count; see the module
/// docs for why it is not the lane count.
const STRIP_ROWS: usize = LANE_MAX;

/// Entries before column 0 of the row buffer. Lane `l` stores column `d - l`,
/// so the widest lane stores `LANE_MAX - 1` columns behind the step and the
/// first steps of a sweep put that below zero. Those cells are masked; the pad
/// only has to exist.
const FRONT: usize = LANE_MAX - 1;

/// The three matrices of one read row, indexed by `FRONT + column`.
///
/// Before a strip it holds the row above the strip; the strip's last lane
/// overwrites it, column by column, into the row below. One buffer serves both
/// because a step reads column `d` before it stores column `d - LANES + 1`,
/// and because every column the next strip reads has either been overwritten
/// or was never non-zero: the columns a strip stores run from `LANES - 1`
/// before its first step to `LANES - 1` before its last, and that range covers
/// the live range of the row it replaces on both sides.
#[derive(Debug, Default)]
pub(crate) struct RowBuffer {
    m: Vec<f32>,
    i: Vec<f32>,
    d: Vec<f32>,
}

impl RowBuffer {
    fn reset(&mut self, len: usize) {
        reset(&mut self.m, len, 0.0);
        reset(&mut self.i, len, 0.0);
        reset(&mut self.d, len, 0.0);
    }
}

impl Workspace {
    /// Scalar `f32`, one row at a time along the haplotype, in this workspace.
    pub fn align_strips<E: Emission>(
        &mut self,
        haplotype: &Haplotype,
        read: &Read,
        emission: &E,
        band: Band,
    ) -> Log10Likelihood {
        self.strips::<f32, E>(haplotype, read, emission, band)
    }

    /// Eight read rows at a time along the haplotype, in this workspace.
    pub fn align_strips_simd<E: Emission>(
        &mut self,
        haplotype: &Haplotype,
        read: &Read,
        emission: &E,
        band: Band,
    ) -> Log10Likelihood {
        self.strips::<f32x8, E>(haplotype, read, emission, band)
    }

    fn strips<L: Lane, E: Emission>(
        &mut self,
        haplotype: &Haplotype,
        read: &Read,
        emission: &E,
        band: Band,
    ) -> Log10Likelihood {
        let (h, r) = (haplotype.len(), read.len());
        if h == 0 || r == 0 {
            return Log10Likelihood::IMPOSSIBLE;
        }
        if self.plan.fill(haplotype, read, emission).is_none() {
            return Log10Likelihood::IMPOSSIBLE;
        }
        strip_kernel::<L>(&self.plan, &mut self.rows, Shape { haplotype: h, read: r }, band)
    }
}

/// Scalar `f32`, one row at a time along the haplotype.
///
/// The same recurrence and band as [`align_banded`], traversed by row rather
/// than by anti-diagonal; see the `strips` module for what that changes.
///
/// [`align_banded`]: crate::align_banded
pub fn align_strips<E: Emission>(
    haplotype: &Haplotype,
    read: &Read,
    emission: &E,
    band: Band,
) -> Log10Likelihood {
    Workspace::new().align_strips(haplotype, read, emission, band)
}

/// Eight read rows at a time along the haplotype.
///
/// Bit-identical to [`align_strips`], for the same reason [`align_banded_simd`]
/// is to [`align_banded`].
///
/// [`align_banded`]: crate::align_banded
/// [`align_banded_simd`]: crate::align_banded_simd
pub fn align_strips_simd<E: Emission>(
    haplotype: &Haplotype,
    read: &Read,
    emission: &E,
    band: Band,
) -> Log10Likelihood {
    Workspace::new().align_strips_simd(haplotype, read, emission, band)
}

/// Where one strip's sweep begins and ends, and where every lane is live.
///
/// Lane `l` is live at step `d` when cell `(r0 + l, d - l)` is inside the band
/// and inside the matrix's columns. Rows past the end of the read are *not*
/// excluded: their plan tracks are all zero, so they compute to exactly zero
/// without a mask, and nothing below them reads them.
struct Sweep {
    /// The first and last step with any live lane.
    first: usize,
    last: usize,
    /// The steps with every lane live, `full_first..=full_last`; empty, with
    /// `full_first > full_last`, when some lane is never live.
    full_first: usize,
    full_last: usize,
    /// Per lane, the first live step and one past the last, as the numbers a
    /// mask compares the step against.
    lane_first: Window,
    lane_past: Window,
}

impl Sweep {
    /// `None` when no lane of the strip is live.
    #[allow(
        clippy::cast_possible_wrap,
        clippy::cast_precision_loss,
        clippy::cast_sign_loss,
        clippy::cast_possible_truncation,
        reason = "steps are bounded by the haplotype length plus the lane count"
    )]
    fn new(shape: Shape, band: Band, r0: usize, lanes: usize) -> Option<Self> {
        let (h, o, w) = (shape.haplotype as i64, band.offset, band.half_width);
        let (mut first, mut last) = (i64::MAX, i64::MIN);
        let (mut full_first, mut full_last) = (i64::MIN, i64::MAX);
        let mut lane_first = [f32::INFINITY; LANE_MAX];
        let mut lane_past = [f32::NEG_INFINITY; LANE_MAX];
        for lane in 0..lanes.min(LANE_MAX) {
            let i = (r0 + lane) as i64;
            let column_first = (i + o - w).max(1);
            let column_last = (i + o + w).min(h);
            if column_first > column_last {
                full_first = i64::MAX;
                continue;
            }
            let (lo, hi) = (column_first + lane as i64, column_last + lane as i64);
            first = first.min(lo);
            last = last.max(hi);
            full_first = full_first.max(lo);
            full_last = full_last.min(hi);
            if let Some(slot) = lane_first.get_mut(lane) {
                *slot = lo as f32;
            }
            if let Some(slot) = lane_past.get_mut(lane) {
                *slot = (hi + 1) as f32;
            }
        }
        if first > last {
            return None;
        }
        let (full_first, full_last) = if full_first > full_last {
            // Empty, with both ends inside the sweep so the arithmetic below
            // stays in `usize`.
            (last + 1, last)
        } else {
            (full_first, full_last)
        };
        Some(Self {
            first: first as usize,
            last: last as usize,
            full_first: full_first as usize,
            full_last: full_last as usize,
            lane_first,
            lane_past,
        })
    }
}

/// What a strip holds in registers for its whole sweep.
struct Strip<'a, L> {
    row: RowLanes<L>,
    transitions: TransitionLanes<L>,
    /// The five column tracks, viewed from the origin the sweep's *last* step
    /// reads: the reversed index falls as the step rises.
    base: View<'a>,
    converted: View<'a>,
    plain: View<'a>,
    rate: View<'a>,
    unconverted: View<'a>,
    /// `last - first`: the window a step reads is at `reach - at`.
    reach: usize,
}

/// The row buffer's slice for one strip, from `FRONT + first - LANE_MAX + 1`,
/// exactly `reach + LANE_MAX` entries long so that every step's offsets have
/// been proven in bounds once; the same reasoning as `View` in `banded`.
///
/// A step at `at = d - first` reads column `d`, at offset `at + LANE_MAX - 1`,
/// and stores column `d - LANES + 1`, at offset `at + LANE_MAX - LANES`. Both
/// are below `reach + LANE_MAX` whenever `at <= reach`.
struct Buffers<'a> {
    m: &'a mut [f32],
    i: &'a mut [f32],
    d: &'a mut [f32],
}

impl<'a> Buffers<'a> {
    /// `None` unless each buffer holds `from..from + reach + LANE_MAX`.
    fn new(rows: &'a mut RowBuffer, from: usize, reach: usize) -> Option<Self> {
        let span = reach.checked_add(LANE_MAX)?;
        let end = from.checked_add(span)?;
        Some(Self {
            m: rows.m.get_mut(from..end)?,
            i: rows.i.get_mut(from..end)?,
            d: rows.d.get_mut(from..end)?,
        })
    }

    /// The row above the strip at the column step `at` reads.
    ///
    /// # Safety
    ///
    /// `at <= reach`, the reach `new` was given.
    #[inline(always)]
    unsafe fn above(&self, at: usize) -> (f32, f32, f32) {
        let offset = at + LANE_MAX - 1;
        debug_assert!(offset < self.m.len());
        // SAFETY: `new` proved every buffer holds `reach + LANE_MAX` entries
        // and the caller keeps `at <= reach`.
        unsafe {
            (
                *self.m.as_ptr().add(offset),
                *self.i.as_ptr().add(offset),
                *self.d.as_ptr().add(offset),
            )
        }
    }

    /// The last lane's cells, into the column step `at` stores.
    ///
    /// # Safety
    ///
    /// `at <= reach`, the reach `new` was given.
    #[inline(always)]
    unsafe fn store<L: Lane>(&mut self, at: usize, cells: (L, L, L)) {
        let offset = at + LANE_MAX - L::LANES;
        debug_assert!(offset < self.m.len());
        // SAFETY: as for `above`, through a unique borrow.
        unsafe {
            *self.m.as_mut_ptr().add(offset) = cells.0.last();
            *self.i.as_mut_ptr().add(offset) = cells.1.last();
            *self.d.as_mut_ptr().add(offset) = cells.2.last();
        }
    }
}

/// What every step feeds: the running maximum the next renormalisation
/// reads, and the read's last row where this strip holds it.
struct Totals<L> {
    /// The lanewise maximum of every step's cells; only its last lane, the
    /// row that crosses into the next strip, is read.
    running: L,
    /// `m + i` of the read's last row, accumulated in that row's lane.
    total: L,
    /// That lane, as a mask, in the strip that holds the last row.
    summed: Option<L>,
}

impl<L: Lane> Totals<L> {
    #[inline(always)]
    fn absorb(&mut self, (m, i, d): (L, L, L)) {
        self.running = self.running.vmax(m).vmax(i).vmax(d);
        if let Some(keep) = self.summed {
            self.total = self.total + keep.select(m + i, L::splat(0.0));
        }
    }
}

/// The vectors that carry from one step to the next.
struct State<L> {
    /// The previous step: every lane's `left`.
    m: L,
    i: L,
    d: L,
    /// The previous step's `up`, which is this step's `diag`. The insertion
    /// and deletion diagonals enter a match through the same transition, so
    /// they are carried as one sum: a lane shift is linear, and adding before
    /// it is adding after it.
    m_diag: L,
    indel_diag: L,
}

/// Subnormals to zero, for the reason `Sink::store` in `banded` gives.
#[inline(always)]
fn flush<L: Lane>(value: L) -> L {
    value.below(L::splat(f32::MIN_POSITIVE)).select(L::splat(0.0), value)
}

/// One step of the sweep: the cells `(r0 + l, first + at - l)` for every lane.
///
/// `mask` is `None` where every lane is live and the lanewise mask otherwise.
/// `None` only past the sweep's reach, which the kernel's loops never are.
#[inline(always)]
fn step<L: Lane>(
    strip: &Strip<'_, L>,
    buffers: &mut Buffers<'_>,
    state: &mut State<L>,
    at: usize,
    mask: Option<L>,
) -> Option<(L, L, L)> {
    // The one check a step makes: `at <= reach` is what every unchecked
    // access below rests on.
    let window = strip.reach.checked_sub(at)?;
    // SAFETY: the buffers were built with this strip's reach and `at <=
    // reach` was just checked.
    let (above_m, above_i, above_d) = unsafe { buffers.above(at) };
    let up_m = state.m.shift_in(above_m);
    let up_i = state.i.shift_in(above_i);
    let up_indel = (state.i + state.d).shift_in(above_i + above_d);

    // SAFETY: every view was built with `reach + LANE_MAX` entries and
    // `window <= reach`.
    let load = |view: View<'_>| unsafe { L::load(view.window(window)) };
    let prior_v = prior::<L>(
        ColumnLanes {
            base: load(strip.base),
            converted: load(strip.converted),
            plain: load(strip.plain),
            rate: load(strip.rate),
            unconverted_rate: load(strip.unconverted),
        },
        strip.row,
    );
    let t = strip.transitions;

    let m = prior_v * (state.m_diag * t.match_to_match + state.indel_diag * t.indel_to_match);
    let i = up_m * t.match_to_insertion + up_i * t.gap_continuation;
    let d = state.m * t.match_to_deletion + state.d * t.gap_continuation;
    let (m, i, d) = (flush(m), flush(i), flush(d));
    let (m, i, d) = match mask {
        Some(keep) => {
            let zero = L::splat(0.0);
            (keep.select(m, zero), keep.select(i, zero), keep.select(d, zero))
        }
        None => (m, i, d),
    };

    // SAFETY: as for `above`.
    unsafe { buffers.store::<L>(at, (m, i, d)) };

    *state = State { m, i, d, m_diag: up_m, indel_diag: up_indel };
    Some((m, i, d))
}

/// The lanewise mask for step `d`: every lane whose live steps include it.
#[inline(always)]
fn live_lanes<L: Lane>(lane_first: L, lane_past: L, d: usize) -> L {
    #[allow(clippy::cast_precision_loss, reason = "a step is a few hundred")]
    let now = L::splat(d as f32);
    lane_first.below(now + L::splat(1.0)).both(now.below(lane_past))
}

/// One phase of a sweep: the steps `from..=to`, masked or not, two at a time
/// so that the state alternates between two sets of registers rather than
/// being copied at the end of every step. `false` past the sweep's reach,
/// which the kernel's phases never are.
#[inline(always)]
#[allow(clippy::too_many_arguments, reason = "the sweep's registers, passed as such")]
fn phase<L: Lane, const MASKED: bool>(
    strip: &Strip<'_, L>,
    buffers: &mut Buffers<'_>,
    state: &mut State<L>,
    totals: &mut Totals<L>,
    edges: (L, L),
    first: usize,
    from: usize,
    to: usize,
) -> bool {
    let mask = |d: usize| MASKED.then(|| live_lanes(edges.0, edges.1, d));
    let mut d = from;
    while d < to {
        let Some(cells) = step(strip, buffers, state, d - first, mask(d)) else {
            return false;
        };
        totals.absorb(cells);
        let Some(cells) = step(strip, buffers, state, d + 1 - first, mask(d + 1)) else {
            return false;
        };
        totals.absorb(cells);
        d += 2;
    }
    if d == to {
        let Some(cells) = step(strip, buffers, state, d - first, mask(d)) else {
            return false;
        };
        totals.absorb(cells);
    }
    true
}

/// Strips of one row per lane, each swept along the haplotype; see the
/// module docs.
#[allow(
    clippy::too_many_lines,
    reason = "one traversal; splitting it would hide the index algebra it exists to get right"
)]
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_possible_wrap,
    reason = "the f32 narrowing is the point of this kernel, and every count here is a few hundred"
)]
pub(crate) fn strip_kernel<L: Lane>(
    plan: &Plan,
    rows: &mut RowBuffer,
    shape: Shape,
    band: Band,
) -> Log10Likelihood {
    let Shape { haplotype: h, read: r } = shape;
    // Columns `0..=h` behind the front pad, and one window past `h` for the
    // lanes that run past the haplotype.
    rows.reset(FRONT + h + 1 + LANE_MAX);

    // The read's free start: row 0 holds `1 / h` in the deletion matrix at
    // every column the band contains.
    let init = 1.0 / h as f32;
    let mut crossing_max = 0.0f32;
    for column in 0..=h {
        if (column as i64 - band.offset).abs() <= band.half_width
            && let Some(cell) = rows.d.get_mut(FRONT + column)
        {
            *cell = init;
            crossing_max = init;
        }
    }

    let mut exponent = 0i32;
    let zero = L::splat(0.0);
    let offsets = L::offsets();
    let mut total = zero;

    let strips = r.div_ceil(L::LANES);
    for s in 0..strips {
        let r0 = 1 + s * L::LANES;
        // The row above this strip, `r0 - 1`, is where the scale changes
        // every `STRIP_ROWS` rows; `crossing_max` is that row's maximum.
        let shift = if (r0 - 1) % STRIP_ROWS == 0 {
            let shift = normalising_shift_f32(crossing_max);
            exponent += shift;
            shift
        } else {
            0
        };

        let Some(sweep) = Sweep::new(shape, band, r0, L::LANES) else {
            // No cell of these rows is inside the band. The row buffer is
            // left alone: it is all zero, because a strip with no live lane
            // has every strip above it dead too.
            continue;
        };
        let reach = sweep.last - sweep.first;
        let span = reach + LANE_MAX;

        let row_lanes = |track: &[f32]| -> Option<L> {
            Some(L::load(track.get(r0..r0.checked_add(LANE_MAX)?)?.first_chunk()?))
        };
        let tracks = &plan.rows;
        let (Some(row_base), Some(spread), Some(mismatched)) =
            (row_lanes(&tracks.base), row_lanes(&tracks.spread), row_lanes(&tracks.mismatched))
        else {
            return Log10Likelihood::IMPOSSIBLE;
        };
        let (Some(match_to_match), Some(match_to_insertion), Some(match_to_deletion)) = (
            row_lanes(&tracks.match_to_match),
            row_lanes(&tracks.match_to_insertion),
            row_lanes(&tracks.match_to_deletion),
        ) else {
            return Log10Likelihood::IMPOSSIBLE;
        };
        let (Some(indel_to_match), Some(gap_continuation)) =
            (row_lanes(&tracks.indel_to_match), row_lanes(&tracks.gap_continuation))
        else {
            return Log10Likelihood::IMPOSSIBLE;
        };

        // Lane `l` at step `d` reads column `d - l`, reversed at
        // `COLUMN_FRONT + h - d + l`: the last step reads from the lowest
        // origin, and `d <= h + LANES - 1` keeps it at or past the front pad.
        let Some(column_origin) = (COLUMN_FRONT + h).checked_sub(sweep.last) else {
            return Log10Likelihood::IMPOSSIBLE;
        };
        let columns = &plan.columns;
        let (Some(base), Some(converted), Some(plain), Some(rate), Some(unconverted)) = (
            View::new(&columns.base, column_origin, span),
            View::new(&columns.converted, column_origin, span),
            View::new(&columns.plain, column_origin, span),
            View::new(&columns.rate, column_origin, span),
            View::new(&columns.unconverted, column_origin, span),
        ) else {
            return Log10Likelihood::IMPOSSIBLE;
        };
        let strip = Strip {
            row: RowLanes::new(row_base, spread, mismatched),
            transitions: TransitionLanes {
                match_to_match,
                match_to_insertion,
                match_to_deletion,
                indel_to_match,
                gap_continuation,
            },
            base,
            converted,
            plain,
            rate,
            unconverted,
            reach,
        };

        // The buffer slice this sweep touches, from the column the first step
        // stores to the column the last step reads.
        let Some(mut buffers) = Buffers::new(rows, sweep.first, reach) else {
            return Log10Likelihood::IMPOSSIBLE;
        };

        // Bring the row above onto this strip's scale. The strip reads
        // columns `first - 1..=last`, at offsets `LANE_MAX - 2..`; whatever
        // lies outside is zero or never read again.
        if shift != 0 {
            let lift = exp2_f32(shift);
            for buffer in [&mut *buffers.m, &mut *buffers.i, &mut *buffers.d] {
                for cell in buffer.iter_mut().skip(LANE_MAX - 2) {
                    *cell *= lift;
                }
            }
        }

        // Before the first step, `diag` is the row above at column `first -
        // 1` in lane 0 and zero elsewhere: every other lane's diagonal
        // neighbour is a cell of this strip before its first live step.
        let Some(prime) = LANE_MAX.checked_sub(2) else {
            return Log10Likelihood::IMPOSSIBLE;
        };
        let (Some(&above_m), Some(&above_i), Some(&above_d)) =
            (buffers.m.get(prime), buffers.i.get(prime), buffers.d.get(prime))
        else {
            return Log10Likelihood::IMPOSSIBLE;
        };
        let mut state = State {
            m: zero,
            i: zero,
            d: zero,
            m_diag: zero.shift_in(above_m),
            indel_diag: zero.shift_in(above_i + above_d),
        };

        // The last strip holds the read's last row in one lane; its match
        // and insertion cells are the total.
        let mut totals = Totals {
            running: zero,
            total,
            summed: (r0 + L::LANES > r).then(|| offsets.equals(L::splat((r - r0) as f32))),
        };
        let edges = (L::load(&sweep.lane_first), L::load(&sweep.lane_past));

        // Three phases: the leading edge of the band, where lanes come live
        // one by one; the middle, with nothing to mask; and the trailing
        // edge, where they go dead again. Each phase's range is empty, with
        // `from > to`, where the sweep has no such steps.
        let Sweep { first, last, full_first, full_last, .. } = sweep;
        let leading = (first, full_first.min(last + 1).saturating_sub(1));
        let trailing = ((full_last + 1).max(first), last);
        let mut ran = phase::<L, true>(
            &strip,
            &mut buffers,
            &mut state,
            &mut totals,
            edges,
            first,
            leading.0,
            leading.1,
        );
        ran = ran
            && phase::<L, false>(
                &strip,
                &mut buffers,
                &mut state,
                &mut totals,
                edges,
                first,
                full_first,
                full_last,
            );
        ran = ran
            && phase::<L, true>(
                &strip,
                &mut buffers,
                &mut state,
                &mut totals,
                edges,
                first,
                trailing.0,
                trailing.1,
            );
        if !ran {
            return Log10Likelihood::IMPOSSIBLE;
        }
        crossing_max = totals.running.last();
        total = totals.total;

        // Column 0 belongs to the free start alone: every row below row 0 is
        // zero there. The eight-lane strip's last lane stores it, as a
        // masked zero, on the way past; the scalar strip's stores begin at
        // its first live column, which is 1, and would leave row 0's `init`
        // for the next row's diagonal prime to read.
        if s == 0 {
            for buffer in [&mut *buffers.m, &mut *buffers.i, &mut *buffers.d] {
                if let Some(cell) = FRONT.checked_sub(sweep.first).and_then(|at| buffer.get_mut(at))
                {
                    *cell = 0.0;
                }
            }
        }
    }

    let total = f64::from(total.horizontal_sum());
    if total > 0.0 {
        Log10Likelihood::new(total.log10() - f64::from(exponent) * core::f64::consts::LOG10_2)
    } else {
        Log10Likelihood::IMPOSSIBLE
    }
}
