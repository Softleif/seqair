//! One alignment per lane: the row sweep the batch and the pairs kernels
//! share.
//!
//! The batch kernel puts eight haplotypes against one read, the pairs kernel
//! eight unrelated pairs. Below their plans they are one kernel: lane `k` is
//! an alignment whose band, in the shifted columns `c = j + delta_k`, is the
//! same `|c - i - offset| <= half_width` for every lane, so the column loop is
//! the band and the three neighbours of a cell are the same lane of the
//! previous column and of the row above. What differs is only where a row
//! comes from -- one read's row as splats for a batch, eight reads' rows
//! transposed for a group of pairs -- and that is [`Rows`]. A batch is the
//! case `delta_k = 0` and every lane's read the same; see `batch` and `pairs`
//! for what each plan proves about its lanes.
//!
//! The layout is blocks of lane-windows: a column's five tracks are one
//! [`ColumnBlock`] and a cell's three matrices one [`Cells`], so the column
//! loop holds two base pointers. With one interleaved `Vec` per track it held
//! eight, which x86-64 cannot keep in registers, and reloaded the spilled
//! ones every step (notes §13.1).
//!
//! The column loop is `steps::<L, MASKED, SUMMING>`, four instances, so that
//! each stretch of a row pays only for what its columns need: the past-end
//! mask only past the group's shortest haplotype, and the sum only on the
//! rows where some lane's read ends (notes §14.5(e)).
//!
//! The arithmetic, the order of the operations, the flush to zero, the
//! renormalisation cadence and the free start are the strip kernel's, cell for
//! cell, so every lane is bit-identical to [`align_strips`] on its own pair.
//!
//! [`align_strips`]: crate::align_strips

use crate::{
    banded::{
        CODE_NO_CONVERSION, CODE_NO_PLAIN_MATCH, ColumnLanes, LANE_MAX, Lane, RowLanes, RowTracks,
        TransitionLanes, Window, prior,
    },
    scaling::{exp2_f32, normalising_shift_f32},
    types::Log10Likelihood,
};

/// Rows between two renormalisations: the strip kernel's constant, and it has
/// to be, or the kernels would not round the same way.
const STRIP_ROWS: usize = LANE_MAX;

/// One shifted column of every lane, a block so that the column loop holds one
/// base pointer for all five tracks.
#[derive(Debug, Clone, Copy)]
pub(crate) struct ColumnBlock {
    pub(crate) base: Window,
    pub(crate) converted: Window,
    pub(crate) plain: Window,
    pub(crate) rate: Window,
    pub(crate) unconverted: Window,
}

impl ColumnBlock {
    /// What a column no lane's band reaches holds: never read, but finite.
    pub(crate) const UNVISITED: Self = Self {
        base: [CODE_NO_PLAIN_MATCH; LANE_MAX],
        converted: [CODE_NO_CONVERSION; LANE_MAX],
        plain: [CODE_NO_PLAIN_MATCH; LANE_MAX],
        rate: [0.0; LANE_MAX],
        unconverted: [0.0; LANE_MAX],
    };
}

/// The three matrices at one shifted column of one read row, every lane.
#[derive(Debug, Default, Clone, Copy)]
pub(crate) struct Cells {
    pub(crate) m: Window,
    pub(crate) i: Window,
    pub(crate) d: Window,
}

/// One read row's entry in every track, for one lane, in this order: the
/// base's code, `(1 - eps) - eps / 3` as in `banded::RowTracks`, `eps / 3`,
/// then the five transitions `match_to_match`, `match_to_insertion`,
/// `match_to_deletion`, `indel_to_match` and `gap_continuation`.
pub(crate) type RowEntry = Window;

/// A row past a lane's read, and every row of a lane with no pair: zeros,
/// which compute to zero.
pub(crate) const NO_ROW: RowEntry = [0.0; LANE_MAX];

/// Where a kernel's read rows come from.
pub(crate) trait Rows<L: Lane> {
    /// Row `row`, one-based, of every lane, in [`RowEntry`]'s track order;
    /// `None` only past what the plan was built for, which the kernel's row
    /// count rules out.
    fn row(&self, token: L::Token, row: usize) -> Option<[L; LANE_MAX]>;
}

/// One read's rows for every lane, as splats: the batch kernel's.
#[derive(Clone, Copy)]
pub(crate) struct SharedRows<'a>(pub(crate) &'a RowTracks);

impl<L: Lane> Rows<L> for SharedRows<'_> {
    #[inline(always)]
    fn row(&self, token: L::Token, row: usize) -> Option<[L; LANE_MAX]> {
        let RowTracks {
            base,
            spread,
            mismatched,
            match_to_match,
            match_to_insertion,
            match_to_deletion,
            indel_to_match,
            gap_continuation,
        } = self.0;
        Some([
            splat_at(token, base, row)?,
            splat_at(token, spread, row)?,
            splat_at(token, mismatched, row)?,
            splat_at(token, match_to_match, row)?,
            splat_at(token, match_to_insertion, row)?,
            splat_at(token, match_to_deletion, row)?,
            splat_at(token, indel_to_match, row)?,
            splat_at(token, gap_continuation, row)?,
        ])
    }
}

/// `track[row]` in every lane. A function, not a closure: a closure inside a
/// `#[simd]` body can lose the level's target features.
#[inline(always)]
fn splat_at<L: Lane>(token: L::Token, track: &[f32], row: usize) -> Option<L> {
    track.get(row).map(|&value| L::splat(token, value))
}

/// Each lane's own read, one [`RowEntry`] per row, transposed eight lanes at a
/// time: the pairs kernel's. A lane's table is its rows `1..=r` at `0..r`, and
/// a lane with no pair has an empty one.
#[derive(Clone, Copy)]
pub(crate) struct LaneRows<'a>(pub(crate) [&'a [RowEntry]; LANE_MAX]);

impl<L: Lane> Rows<L> for LaneRows<'_> {
    #[inline(always)]
    fn row(&self, token: L::Token, row: usize) -> Option<[L; LANE_MAX]> {
        let mut entries = [&NO_ROW; LANE_MAX];
        for (entry, table) in entries.iter_mut().zip(&self.0) {
            if let Some(found) = table.get(row.wrapping_sub(1)) {
                *entry = found;
            }
        }
        Some(L::transpose(token, entries))
    }
}

/// What one group's sweep reads besides its rows.
#[derive(Clone, Copy)]
pub(crate) struct LanesView<'a> {
    /// Shifted column `c` is `columns[c - origin]`, for every `c` a row's band
    /// reaches.
    pub(crate) columns: &'a [ColumnBlock],
    pub(crate) origin: usize,
    /// `1 / h` per lane, the read's free start; zero for a lane that cannot
    /// align.
    pub(crate) init: &'a Window,
    /// `delta_k` per lane: row 0 is live from here.
    pub(crate) front: &'a Window,
    /// `delta_k + h_k + 1` per lane, so `c < past_end[k]` is `j <= h_k`; zero
    /// for a lane with no haplotype.
    pub(crate) past_end: &'a Window,
    /// Each lane's read length; zero for a lane that cannot align.
    pub(crate) lengths: &'a [usize; LANE_MAX],
    /// The longest read.
    pub(crate) rows_len: usize,
    /// The last shifted column any lane has.
    pub(crate) width: usize,
    /// The band in shifted columns.
    pub(crate) offset: i64,
    pub(crate) half_width: i64,
}

/// The previous column's cells and the diagonal carried across a column step.
struct Carry<L> {
    left_m: L,
    left_d: L,
    diag_m: L,
    diag_indel: L,
}

/// What a row's column steps carry from one to the next, and accumulate.
struct Sweep<L> {
    carry: Carry<L>,
    running: L,
    total: L,
}

/// What a row's column steps read and never change.
#[derive(Clone, Copy)]
struct RowConstants<L: Lane> {
    lanes: RowLanes<L>,
    t: TransitionLanes<L>,
    past_end: L,
    /// The lanes whose read ends on this row, for the `SUMMING` steps.
    ends: L::Mask,
}

/// Subnormals to zero, as in the strip kernel and for the same reason.
#[inline(always)]
fn flush<L: Lane>(value: L) -> L {
    L::masked_out(value.below(L::splat(value.token(), f32::MIN_POSITIVE)), value)
}

/// The column steps over `cells` and `columns`, which are cut to the same
/// columns, the first of which is `from`.
///
/// `MASKED` zeroes the lanes whose haplotype ended before the column, which
/// only the columns past the group's shortest haplotype need, and `SUMMING`
/// adds the cells of the lanes whose read ends on this row into the total.
/// Both are constants so that each of the four loops pays only for what its
/// columns need: the mask is a compare, an add and two `and`s per step, and a
/// runtime `summing` was a compare and a branch per step on x86 and, on NEON,
/// a select and a round trip of the total through the stack.
#[inline(always)]
#[allow(clippy::cast_precision_loss, reason = "a column is a few hundred")]
fn steps<L: Lane, const MASKED: bool, const SUMMING: bool>(
    token: L::Token,
    row: RowConstants<L>,
    sweep: &mut Sweep<L>,
    cells: &mut [Cells],
    columns: &[ColumnBlock],
    from: usize,
) {
    debug_assert_eq!(cells.len(), columns.len());
    let t = row.t;
    // The column index as a lane, carried and incremented rather than
    // converted per step. `column as f32` is a `vcvtsi2ss` -- two uops, a
    // false dependency on the destination register and ~5 cycles -- plus a
    // broadcast; one `vaddps` replaces both, and integers this small are
    // exact in `f32`, so the value is the same bit for bit.
    let one = L::splat(token, 1.0);
    let mut column_lane = L::splat(token, from as f32);
    let Sweep { carry, running, total } = sweep;

    for (cell, column) in cells.iter_mut().zip(columns) {
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
            row.lanes,
        );

        let m = prior_v * (carry.diag_m * t.match_to_match + carry.diag_indel * t.indel_to_match);
        let i = up_m * t.match_to_insertion + up_i * t.gap_continuation;
        let d = carry.left_m * t.match_to_deletion + carry.left_d * t.gap_continuation;
        let (m, i, d) = (flush(m), flush(i), flush(d));
        // `i` needs no mask: past a lane's haplotype its `m` is masked to
        // zero at every row, and `i` reads only the cell above --
        // `up_m * match_to_insertion + up_i * gap_continuation` -- so a dead
        // lane's `i` starts at zero and stays there. `d` does need one: it
        // reads the cell to its *left*, which at the column just past the
        // haplotype is still live.
        let (m, d) = if MASKED {
            let keep = column_lane.below(row.past_end);
            column_lane = column_lane + one;
            (L::masked(keep, m), L::masked(keep, d))
        } else {
            (m, d)
        };

        m.store(&mut cell.m);
        i.store(&mut cell.i);
        d.store(&mut cell.d);

        *running = running.vmax(m).vmax(i).vmax(d);
        if SUMMING {
            // The strip kernel's form, parentheses and mask included: `f32`
            // addition is not associative, and a lane that does not end here
            // adds an exact zero.
            *total = *total + L::masked(row.ends, m + i);
        }
        carry.left_m = m;
        carry.left_d = d;
        carry.diag_m = up_m;
        carry.diag_indel = up_i + up_d;
    }
}

/// One row at a time along the shifted columns, every lane its own alignment.
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
pub(crate) fn lanes_kernel<L: Lane, R: Rows<L>>(
    token: L::Token,
    view: LanesView<'_>,
    rows: &R,
    cells: &mut Vec<Cells>,
) -> [Log10Likelihood; LANE_MAX] {
    let impossible = [Log10Likelihood::IMPOSSIBLE; LANE_MAX];
    let (width, read_len) = (view.width, view.rows_len);
    if width == 0 || read_len == 0 {
        return impossible;
    }
    let (o, w) = (view.offset, view.half_width);
    let span = width + 2;
    cells.clear();
    cells.resize(span, Cells::default());
    // A slice, not the `Vec`: on the SIMD lane this kernel is inlined into a
    // `#[simd]` closure, where the buffer arrives as a field of the closure
    // rather than as a `noalias` argument, and LLVM then has to assume every
    // store into a row may have moved the `Vec` and reloads its header every
    // column. Locals cannot move.
    let Some(cells) = cells.get_mut(..span) else {
        return impossible;
    };
    let (columns, origin) = (view.columns, view.origin);

    let zero = L::splat(token, 0.0);
    let one = L::splat(token, 1.0);
    let past_end = L::load(token, view.past_end);
    let init = L::load(token, view.init);

    // Row 0, the free start: `1 / h_k` in the deletion matrix at every column
    // of the band that lane `k`'s haplotype reaches, which starts at its own
    // column 0, `delta_k`. Left of that, row 0 has to be zero: see `pairs`
    // for why that alone keeps every row clear there.
    let front = L::load(token, view.front);
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

    // The last column at which every lane with a haplotype is still live. A
    // lane without one needs no mask: its columns are sentinels and zero
    // rates and its free start is zero, so every one of its cells is a
    // product with a zero whether masked or not.
    let shortest = view
        .past_end
        .iter()
        .filter(|&&past| past > 0.0)
        .fold(f32::INFINITY, |shortest, &past| shortest.min(past - 1.0));
    let shortest = if shortest.is_finite() { shortest as usize } else { width };

    let mut exponent = [0i32; LANE_MAX];
    let mut total = zero;
    let mut scratch = [0.0f32; LANE_MAX];
    let lanes = L::LANES.min(LANE_MAX);

    for row in 1..=read_len {
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
                let inside = view.lengths.get(lane).is_some_and(|&r| row <= r);
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

        let Some(
            [
                base,
                spread,
                mismatched,
                match_to_match,
                match_to_insertion,
                match_to_deletion,
                indel_to_match,
                gap_continuation,
            ],
        ) = rows.row(token, row)
        else {
            return impossible;
        };

        // The lanes whose read ends on this row, as a mask, built from the
        // integer lengths so that no `f32` compare of a row number is
        // involved; no mask on the rows where no lane ends.
        let summing = view.lengths.iter().take(lanes).any(|&r| r == row);
        let mut ends = [0.0f32; LANE_MAX];
        if summing {
            for (slot, &r) in ends.iter_mut().zip(view.lengths) {
                *slot = if r == row { 1.0 } else { 0.0 };
            }
        }
        let constants = RowConstants {
            lanes: RowLanes::new(base, spread, mismatched),
            t: TransitionLanes {
                match_to_match,
                match_to_insertion,
                match_to_deletion,
                indel_to_match,
                gap_continuation,
            },
            past_end,
            ends: L::load(token, &ends).equals(one),
        };

        let mut running = zero;
        if first <= last {
            let (first, last) = (first as usize, last as usize);
            let Some(before) = cells.get(first - 1) else {
                return impossible;
            };
            let carry = Carry {
                left_m: zero,
                left_d: zero,
                diag_m: L::load(token, &before.m),
                diag_indel: L::load(token, &before.i) + L::load(token, &before.d),
            };
            let mut sweep = Sweep { carry, running: zero, total };

            // Every lane is live up to the shortest haplotype, and the rest
            // of the row is masked; either part may be empty.
            let alive = last.min(shortest).max(first - 1);
            let (Some(band), Some(band_columns)) =
                (cells.get_mut(first..=last), columns.get(first - origin..=last - origin))
            else {
                return impossible;
            };
            let (live_cells, masked_cells) = band.split_at_mut(alive + 1 - first);
            let (live_columns, masked_columns) = band_columns.split_at(alive + 1 - first);
            if summing {
                steps::<L, false, true>(
                    token,
                    constants,
                    &mut sweep,
                    live_cells,
                    live_columns,
                    first,
                );
                steps::<L, true, true>(
                    token,
                    constants,
                    &mut sweep,
                    masked_cells,
                    masked_columns,
                    alive + 1,
                );
            } else {
                steps::<L, false, false>(
                    token,
                    constants,
                    &mut sweep,
                    live_cells,
                    live_columns,
                    first,
                );
                steps::<L, true, false>(
                    token,
                    constants,
                    &mut sweep,
                    masked_cells,
                    masked_columns,
                    alive + 1,
                );
            }
            running = sweep.running;
            total = sweep.total;

            // The two cells the next row reads that this one did not write:
            // its own `first - 1`, where the band was clamped at column 1, and
            // `last + 1`, where the band still has room to move right.
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
