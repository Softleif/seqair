//! One alignment per lane: a batch of haplotypes scored against one read in
//! lockstep, row by row.
//!
//! The two banded kernels put eight *cells of one alignment* in a vector, so
//! every lane of every step belongs to the same pair -- and pays for it: the
//! neighbours are lane shifts, the band's edges cross a vector diagonally and
//! have to be masked lane by lane, and ~23 % of the cells a strip computes are
//! outside the band. This kernel puts eight *alignments* in a vector instead.
//! Lane `k` holds candidate haplotype `k` at the same `(row, column)` as every
//! other lane, so:
//!
//! - the three neighbours are the same lane of the previous column and of the
//!   row buffer -- no shifts at all;
//! - the band is a property of the row, not of the lane, so the column loop
//!   *is* the band and there is no per-lane edge mask. What stays is one
//!   lanewise mask for haplotypes shorter than the longest in the batch;
//! - the read's row is the same for every lane, so its base, its two emission
//!   terms and its five transitions are splats held for the whole row;
//! - the haplotype columns are interleaved across the batch, so lane `k`'s
//!   column `j` is element `j * BATCH + k` and a column is one contiguous
//!   load, as in the strip kernel.
//!
//! The arithmetic, the order of the operations, the flush to zero, the
//! power-of-two renormalisation every [`STRIP_ROWS`] rows and the free start
//! are the strip kernel's, cell for cell. So a batch is bit-identical to
//! running [`align_strips`] once per haplotype, which is the gate this kernel
//! is worth nothing without.
//!
//! What it does *not* do is make a lane free. A batch of `n < BATCH`
//! haplotypes computes the full eight lanes and throws `8 - n` away, so the
//! per-alignment cost falls only with the fill. See the crate's notes for the
//! step-count algebra.
//!
//! [`align_strips`]: crate::align_strips
//! [`STRIP_ROWS`]: crate::batch::STRIP_ROWS

use fearless_simd::Level;

use crate::{
    banded::{
        Band, CODE_NO_CONVERSION, CODE_NO_PLAIN_MATCH, ColumnLanes, LANE_MAX, Lane, RowLanes,
        TransitionLanes, Window, Workspace, code, prior, reset,
    },
    emission::Emission,
    haplotype::Haplotype,
    read::Read,
    scaling::{exp2_f32, normalising_shift_f32},
    types::Log10Likelihood,
};

/// Alignments per batch: the lane count, which is what the whole kernel is.
pub const BATCH: usize = LANE_MAX;

/// How many haplotypes a batch needs before it beats scoring them one at a
/// time through the strip kernel.
///
/// A batch computes all [`BATCH`] lanes whatever the caller asked for, so its
/// cost is nearly flat in the group size and the per-alignment cost falls only
/// with the fill. Measured against the strip kernel *on the same lane*, over a
/// 200 bp haplotype through the default band, at both read lengths (the two
/// agreed to within 0.02x at every point):
///
/// | haplotypes | 1 | 2 | 3 | 4 | 5 | **6** | 7 | 8 |
/// |---|---|---|---|---|---|---|---|---|
/// | M4 Pro (NEON) | 0.18x | 0.36x | 0.53x | 0.70x | 0.86x | **1.03x** | 1.18x | 1.33x |
/// | 3950X (AVX2) | 0.28x | 0.54x | 0.79x | **1.02x** | 1.25x | 1.48x | 1.69x | 1.90x |
///
/// **The two machines disagree -- the M4 crosses at six, the 3950X at four --
/// and six is the number that is safe on both.** Four would put the M4 into the
/// batch kernel at half fill, where its strip kernel is 1.4x faster; six leaves
/// the 3950X's half-full groups about 1.15x on the table. The number is the
/// ratio between the two kernels' per-cell costs and nothing else, so a
/// per-target constant is defensible -- but not on two machines. See §9.5 of
/// `docs/notes/benchmarking.md`.
///
/// `examples/batchfill.rs` regenerates the table, racing the two kernels on
/// the same lane, so the crossover it reports is the cost of an empty lane and
/// not the gap between two lanes. (The table was measured on the hand-written
/// intrinsics lane that `fearless_simd` replaced; §11 of the notes has the
/// replacement at parity on both kernels.)
///
/// Whatever a machine says, the number is bounded below: the batch kernel
/// cannot be worth running below `BATCH / 2` lanes, because at that fill it is
/// doing more than twice the work for the same answers.
pub const BATCH_BREAK_EVEN: usize = 6;

// A break-even above `BATCH` would mean the batch kernel never runs, and one at
// or below `BATCH / 2` would mean it runs where it is doing more than twice the
// work for the same answers.
const _: () = assert!(
    BATCH_BREAK_EVEN > BATCH / 2 && BATCH_BREAK_EVEN <= BATCH,
    "the break-even has to be reachable and above half fill"
);

/// One batch of haplotypes scored, plus how many lanes that kernel fills --
/// which is how many haplotypes a group may hold, so it travels with it.
///
/// A function pointer rather than a type parameter because the eight-lane
/// kernel is not one type: it is an instance per SIMD level, chosen at run
/// time by `dispatch!` inside the function this points to. One indirect call
/// per batch is nothing, and every caller passes a constant.
#[derive(Clone, Copy)]
pub(crate) struct BatchKernel {
    lanes: usize,
    run: fn(&BatchPlan, &mut BatchBuffer, usize, Band) -> [Log10Likelihood; BATCH],
}

impl BatchKernel {
    /// [`batch_kernel`] over one lane type, for the paths that can name it.
    const fn over<L: Lane<Token = ()>>() -> Self {
        Self {
            lanes: if L::LANES < BATCH { L::LANES } else { BATCH },
            run: |plan, buffer, read_len, band| batch_kernel::<L>((), plan, buffer, read_len, band),
        }
    }
}

/// The batch kernel at the best SIMD level this CPU has.
const SIMD_BATCH_KERNEL: BatchKernel = BatchKernel { lanes: BATCH, run: batch_kernel_simd };

/// [`SIMD_BATCH_KERNEL`]'s function, in [`BatchKernel`]'s shape.
fn batch_kernel_simd(
    plan: &BatchPlan,
    buffer: &mut BatchBuffer,
    read_len: usize,
    band: Band,
) -> [Log10Likelihood; BATCH] {
    crate::simd::batch_kernel_at(Level::new(), plan, buffer, read_len, band)
}

/// Rows between two renormalisations. The strip kernel's constant, and it has
/// to be, or the two would not round the same way.
const STRIP_ROWS: usize = LANE_MAX;

/// Per read row, exactly the strip kernel's row tracks; every lane reads the
/// same entry, so they are scalars here rather than vectors.
#[derive(Debug, Default)]
struct Rows {
    base: Vec<f32>,
    /// `(1 - eps) - eps / 3`, as in `banded::RowTracks`.
    spread: Vec<f32>,
    mismatched: Vec<f32>,
    match_to_match: Vec<f32>,
    match_to_insertion: Vec<f32>,
    match_to_deletion: Vec<f32>,
    indel_to_match: Vec<f32>,
    gap_continuation: Vec<f32>,
}

/// Per haplotype column, interleaved across the batch: index `j * BATCH + k`
/// is lane `k`'s column `j`, one-based, so `j = 0` is the free-start column
/// and `j = width + 1` the pad the last row's `up` reads.
#[derive(Debug, Default)]
struct Columns {
    base: Vec<f32>,
    converted: Vec<f32>,
    plain: Vec<f32>,
    rate: Vec<f32>,
    unconverted: Vec<f32>,
}

/// The three matrices of one read row, interleaved across the batch.
#[derive(Debug, Default)]
pub(crate) struct BatchBuffer {
    m: Vec<f32>,
    i: Vec<f32>,
    d: Vec<f32>,
}

/// Everything one batch hoists out of its loops.
#[derive(Debug, Default)]
pub(crate) struct BatchPlan {
    rows: Rows,
    columns: Columns,
    /// `1 / h` per lane, the read's free start; zero for an absent or empty
    /// haplotype, which then scores [`Log10Likelihood::IMPOSSIBLE`].
    init: Window,
    /// `h + 1` per lane, so `column < lengths[k]` is `column <= h`.
    past_end: Window,
    /// The longest haplotype in the batch.
    width: usize,
}

impl BatchPlan {
    /// `None` only if the read or a haplotype cannot be addressed, which their
    /// constructors rule out.
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_precision_loss,
        reason = "the f32 narrowing is the point of this kernel, and a haplotype is a few hundred bases"
    )]
    fn fill<E: Emission>(
        &mut self,
        haplotypes: &[&Haplotype],
        read: &Read,
        emission: &E,
        band: Band,
    ) -> Option<()> {
        let r = read.len();
        self.width = haplotypes.iter().map(|h| h.len()).max().unwrap_or(0);
        let Self { rows, columns, init, past_end, width } = self;
        let width = *width;

        let Rows {
            base,
            spread,
            mismatched,
            match_to_match,
            match_to_insertion,
            match_to_deletion,
            indel_to_match,
            gap_continuation,
        } = rows;
        for track in [
            &mut *base,
            spread,
            mismatched,
            match_to_match,
            match_to_insertion,
            match_to_deletion,
            indel_to_match,
            gap_continuation,
        ] {
            reset(track, r + 1, 0.0);
        }
        for index in 0..r {
            let observation = read.observation(index)?;
            let eps = emission.epsilon(observation);
            let t = read.transition(index)?;
            let row = index + 1;
            *base.get_mut(row)? = code(observation.base);
            *spread.get_mut(row)? = ((1.0 - eps) - eps / 3.0) as f32;
            *mismatched.get_mut(row)? = (eps / 3.0) as f32;
            *match_to_match.get_mut(row)? = t.match_to_match as f32;
            *match_to_insertion.get_mut(row)? = t.match_to_insertion as f32;
            *match_to_deletion.get_mut(row)? = t.match_to_deletion as f32;
            *indel_to_match.get_mut(row)? = t.indel_to_match as f32;
            *gap_continuation.get_mut(row)? = t.gap_continuation as f32;
        }

        let span = (width + 2) * BATCH;
        let Columns { base, converted, plain, rate, unconverted } = columns;
        reset(base, span, CODE_NO_PLAIN_MATCH);
        reset(converted, span, CODE_NO_CONVERSION);
        reset(plain, span, CODE_NO_PLAIN_MATCH);
        reset(rate, span, 0.0);
        reset(unconverted, span, 0.0);
        *init = [0.0; BATCH];
        *past_end = [0.0; BATCH];

        let strand = read.strand();
        for (lane, haplotype) in haplotypes.iter().enumerate().take(BATCH) {
            let h = haplotype.len();
            if h == 0 {
                continue;
            }
            *init.get_mut(lane)? = 1.0 / h as f32;
            *past_end.get_mut(lane)? = (h + 1) as f32;
            // Only the columns the band can reach. Everything else keeps the
            // sentinels `reset` just wrote, which is what an unvisited column
            // holds anyway -- the kernel's column loop is the band, so it never
            // reads them. Deriving them cost ~29% of a call at 200 bp
            // haplotypes and grew without bound with haplotype length, because
            // the band's span is `read + width` however long the haplotype is.
            let Some((lo, hi)) = band.columns(h, r) else { continue };
            for index in lo..=hi {
                let weights = emission.site_weights(haplotype.site(index)?, strand);
                let site_base = code(weights.base);
                let at = (index + 1) * BATCH + lane;
                *base.get_mut(at)? = site_base;
                match weights.converted {
                    Some(converted_base) => {
                        *converted.get_mut(at)? = code(converted_base);
                        *rate.get_mut(at)? = weights.rate as f32;
                        *unconverted.get_mut(at)? = (1.0 - weights.rate) as f32;
                    }
                    None => *plain.get_mut(at)? = site_base,
                }
            }
        }
        Some(())
    }
}

impl Workspace {
    /// One read against its candidate haplotypes, through whichever kernel is
    /// faster for the number of them. **This is the entry point a caller
    /// wants**; [`Workspace::align_batch`] and [`Workspace::align_strips_simd`]
    /// are the two kernels it chooses between.
    ///
    /// The scores replace `out`'s contents, in the haplotypes' order.
    ///
    /// The choice is per batch, not per call, and it has to be: a caller with
    /// nine haplotypes gets one full batch and a remainder of one, and running
    /// the remainder through eight lanes to score a single haplotype is what
    /// made nine haplotypes *slower* than eight before this existed. Full
    /// batches go through the batch kernel, a remainder shorter than
    /// [`BATCH_BREAK_EVEN`] goes through the strip kernel one haplotype at a
    /// time, and because the two kernels are bit-identical the split cannot
    /// change a score.
    pub fn align_candidates<E: Emission>(
        &mut self,
        haplotypes: &[&Haplotype],
        read: &Read,
        emission: &E,
        band: Band,
        out: &mut Vec<Log10Likelihood>,
    ) {
        out.clear();
        for group in haplotypes.chunks(BATCH) {
            if group.len() >= BATCH_BREAK_EVEN {
                self.batch_group::<E>(group, read, emission, band, out, SIMD_BATCH_KERNEL);
            } else {
                for haplotype in group {
                    out.push(self.align_strips_simd(haplotype, read, emission, band));
                }
            }
        }
    }

    /// Up to [`BATCH`] haplotypes against one read, one haplotype per lane.
    ///
    /// The scores replace `out`'s contents, in the haplotypes' order. More than
    /// [`BATCH`] haplotypes are scored in successive batches; a batch short of
    /// [`BATCH`] still computes every lane, so the cost per alignment falls
    /// only with the fill. Prefer [`Workspace::align_candidates`], which is
    /// this kernel where it wins and the strip kernel where it does not.
    pub fn align_batch<E: Emission>(
        &mut self,
        haplotypes: &[&Haplotype],
        read: &Read,
        emission: &E,
        band: Band,
        out: &mut Vec<Log10Likelihood>,
    ) {
        self.batch::<E>(haplotypes, read, emission, band, out, SIMD_BATCH_KERNEL);
    }

    /// [`Workspace::align_batch`] at a given `fearless_simd` level, so the
    /// parity tests can hold every level the CPU has to the scalar kernel.
    #[doc(hidden)]
    pub fn align_batch_at<E: Emission>(
        &mut self,
        level: Level,
        haplotypes: &[&Haplotype],
        read: &Read,
        emission: &E,
        band: Band,
        out: &mut Vec<Log10Likelihood>,
    ) {
        out.clear();
        for group in haplotypes.chunks(BATCH) {
            let scores =
                if read.is_empty() || self.batch_plan.fill(group, read, emission, band).is_none() {
                    [Log10Likelihood::IMPOSSIBLE; BATCH]
                } else {
                    crate::simd::batch_kernel_at(
                        level,
                        &self.batch_plan,
                        &mut self.batch_rows,
                        read.len(),
                        band,
                    )
                };
            out.extend(scores.iter().take(group.len()).copied());
        }
    }

    /// [`Workspace::align_batch`] with one lane, the bit-parity oracle.
    pub fn align_batch_scalar<E: Emission>(
        &mut self,
        haplotypes: &[&Haplotype],
        read: &Read,
        emission: &E,
        band: Band,
        out: &mut Vec<Log10Likelihood>,
    ) {
        self.batch::<E>(haplotypes, read, emission, band, out, BatchKernel::over::<f32>());
    }

    fn batch<E: Emission>(
        &mut self,
        haplotypes: &[&Haplotype],
        read: &Read,
        emission: &E,
        band: Band,
        out: &mut Vec<Log10Likelihood>,
        kernel: BatchKernel,
    ) {
        out.clear();
        for group in haplotypes.chunks(kernel.lanes) {
            self.batch_group::<E>(group, read, emission, band, out, kernel);
        }
    }

    /// One batch of at most [`BATCH`] haplotypes, appended to `out`.
    fn batch_group<E: Emission>(
        &mut self,
        group: &[&Haplotype],
        read: &Read,
        emission: &E,
        band: Band,
        out: &mut Vec<Log10Likelihood>,
        kernel: BatchKernel,
    ) {
        let scores =
            if read.is_empty() || self.batch_plan.fill(group, read, emission, band).is_none() {
                [Log10Likelihood::IMPOSSIBLE; BATCH]
            } else {
                (kernel.run)(&self.batch_plan, &mut self.batch_rows, read.len(), band)
            };
        out.extend(scores.iter().take(group.len()).copied());
    }
}

/// Up to [`BATCH`] haplotypes against one read, one haplotype per lane.
pub fn align_batch<E: Emission>(
    haplotypes: &[&Haplotype],
    read: &Read,
    emission: &E,
    band: Band,
) -> Vec<Log10Likelihood> {
    let mut out = Vec::with_capacity(haplotypes.len());
    Workspace::new().align_batch(haplotypes, read, emission, band, &mut out);
    out
}

/// One read against its candidate haplotypes, through whichever kernel is
/// faster for the number of them.
///
/// Allocates a [`Workspace`] per call; keep one and call
/// [`Workspace::align_candidates`] when scoring more than one read.
pub fn align_candidates<E: Emission>(
    haplotypes: &[&Haplotype],
    read: &Read,
    emission: &E,
    band: Band,
) -> Vec<Log10Likelihood> {
    let mut out = Vec::with_capacity(haplotypes.len());
    Workspace::new().align_candidates(haplotypes, read, emission, band, &mut out);
    out
}

/// The three matrices of the previous read row, plus the diagonal carried
/// across a column step.
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

/// One row at a time along the haplotype, every lane a different haplotype.
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
pub(crate) fn batch_kernel<L: Lane>(
    token: L::Token,
    plan: &BatchPlan,
    buffer: &mut BatchBuffer,
    read_len: usize,
    band: Band,
) -> [Log10Likelihood; BATCH] {
    let impossible = [Log10Likelihood::IMPOSSIBLE; BATCH];
    let width = plan.width;
    if width == 0 {
        return impossible;
    }
    let span = (width + 2) * BATCH;
    for track in [&mut buffer.m, &mut buffer.i, &mut buffer.d] {
        reset(track, span, 0.0);
    }

    let zero = L::splat(token, 0.0);
    let past_end = L::load(token, &plan.past_end);
    let init = L::load(token, &plan.init);
    let (o, w) = (band.offset, band.half_width);

    // Row 0, the read's free start: `1 / h` in the deletion matrix at every
    // column of the band that the lane's own haplotype reaches.
    let mut crossing = zero;
    for column in 0..=width {
        if (column as i64 - o).abs() > w {
            continue;
        }
        let live = L::splat(token, column as f32).below(past_end);
        let cell = L::masked(live, init);
        let at = column * BATCH;
        let Some(slot) = buffer.d.get_mut(at..at + LANE_MAX) else {
            return impossible;
        };
        let Some(slot) = slot.first_chunk_mut::<LANE_MAX>() else {
            return impossible;
        };
        cell.store(slot);
        crossing = crossing.vmax(cell);
    }

    let mut exponent = [0i32; BATCH];
    let mut total = zero;
    let mut scratch = [0.0f32; LANE_MAX];

    for row in 1..=read_len {
        if (row - 1) % STRIP_ROWS == 0 {
            crossing.store(&mut scratch);
            let mut lift = [0.0f32; LANE_MAX];
            let mut any = false;
            for lane in 0..L::LANES.min(BATCH) {
                let shift = scratch.get(lane).copied().map_or(0, normalising_shift_f32);
                any |= shift != 0;
                if let (Some(slot), Some(sum)) = (lift.get_mut(lane), exponent.get_mut(lane)) {
                    *slot = exp2_f32(shift);
                    *sum += shift;
                }
            }
            if any {
                let lift = L::load(token, &lift);
                for track in [&mut buffer.m, &mut buffer.i, &mut buffer.d] {
                    for chunk in track.as_chunks_mut::<LANE_MAX>().0 {
                        (L::load(token, chunk) * lift).store(chunk);
                    }
                }
            }
        }

        let Some(row_base) = plan.rows.base.get(row).copied() else {
            return impossible;
        };
        let (Some(spread), Some(mismatched)) =
            (plan.rows.spread.get(row).copied(), plan.rows.mismatched.get(row).copied())
        else {
            return impossible;
        };
        let lanes = RowLanes::new(
            L::splat(token, row_base),
            L::splat(token, spread),
            L::splat(token, mismatched),
        );
        let (
            Some(match_to_match),
            Some(match_to_insertion),
            Some(match_to_deletion),
            Some(indel_to_match),
            Some(gap_continuation),
        ) = (
            plan.rows.match_to_match.get(row).copied(),
            plan.rows.match_to_insertion.get(row).copied(),
            plan.rows.match_to_deletion.get(row).copied(),
            plan.rows.indel_to_match.get(row).copied(),
            plan.rows.gap_continuation.get(row).copied(),
        )
        else {
            return impossible;
        };
        let t = TransitionLanes {
            match_to_match: L::splat(token, match_to_match),
            match_to_insertion: L::splat(token, match_to_insertion),
            match_to_deletion: L::splat(token, match_to_deletion),
            indel_to_match: L::splat(token, indel_to_match),
            gap_continuation: L::splat(token, gap_continuation),
        };

        let first = (row as i64 + o - w).max(1);
        let last = (row as i64 + o + w).min(width as i64);
        let summing = row == read_len;
        let mut running = zero;
        if first <= last {
            let (first, last) = (first as usize, last as usize);
            // The row's own left neighbour at its first column is outside the
            // band, and the diagonal there is the previous row's `first - 1`.
            let mut carry = Carry { left_m: zero, left_d: zero, diag_m: zero, diag_indel: zero };
            let at = (first - 1) * BATCH;
            let (Some(m), Some(i), Some(d)) = (
                buffer.m.get(at..at + LANE_MAX).and_then(<[f32]>::first_chunk::<LANE_MAX>),
                buffer.i.get(at..at + LANE_MAX).and_then(<[f32]>::first_chunk::<LANE_MAX>),
                buffer.d.get(at..at + LANE_MAX).and_then(<[f32]>::first_chunk::<LANE_MAX>),
            ) else {
                return impossible;
            };
            carry.diag_m = L::load(token, m);
            carry.diag_indel = L::load(token, i) + L::load(token, d);

            // The column index as a lane, carried and incremented rather than
            // converted per step. `column as f32` is a `vcvtsi2ss` -- two uops,
            // a false dependency on the destination register and ~5 cycles --
            // plus a broadcast; one `vaddps` replaces both, and integers this
            // small are exact in `f32`, so the value is the same bit for bit.
            let one = L::splat(token, 1.0);
            let mut column_lane = L::splat(token, first as f32);

            for column in first..=last {
                let at = column * BATCH;
                let (Some(up_m), Some(up_i), Some(up_d)) = (
                    buffer.m.get(at..at + LANE_MAX).and_then(<[f32]>::first_chunk::<LANE_MAX>),
                    buffer.i.get(at..at + LANE_MAX).and_then(<[f32]>::first_chunk::<LANE_MAX>),
                    buffer.d.get(at..at + LANE_MAX).and_then(<[f32]>::first_chunk::<LANE_MAX>),
                ) else {
                    return impossible;
                };
                let (up_m, up_i, up_d) =
                    (L::load(token, up_m), L::load(token, up_i), L::load(token, up_d));

                let (Some(base), Some(converted), Some(plain), Some(rate), Some(unconverted)) = (
                    plan.columns.base.get(at..at + LANE_MAX).and_then(<[f32]>::first_chunk),
                    plan.columns.converted.get(at..at + LANE_MAX).and_then(<[f32]>::first_chunk),
                    plan.columns.plain.get(at..at + LANE_MAX).and_then(<[f32]>::first_chunk),
                    plan.columns.rate.get(at..at + LANE_MAX).and_then(<[f32]>::first_chunk),
                    plan.columns.unconverted.get(at..at + LANE_MAX).and_then(<[f32]>::first_chunk),
                ) else {
                    return impossible;
                };
                let prior_v = prior::<L>(
                    ColumnLanes {
                        base: L::load(token, base),
                        converted: L::load(token, converted),
                        plain: L::load(token, plain),
                        rate: L::load(token, rate),
                        unconverted_rate: L::load(token, unconverted),
                    },
                    lanes,
                );

                let m = prior_v
                    * (carry.diag_m * t.match_to_match + carry.diag_indel * t.indel_to_match);
                let i = up_m * t.match_to_insertion + up_i * t.gap_continuation;
                let d = carry.left_m * t.match_to_deletion + carry.left_d * t.gap_continuation;
                let (m, i, d) = (flush(m), flush(i), flush(d));
                let keep = column_lane.below(past_end);
                column_lane = column_lane + one;
                // `i` needs no mask: past a lane's haplotype its `m` is masked
                // to zero at every row, and `i` reads only the cell above --
                // `up_m * match_to_insertion + up_i * gap_continuation` -- so a
                // dead lane's `i` starts at zero and stays there. `d` does need
                // one: it reads the cell to its *left*, which at the column
                // just past the haplotype is still live.
                let (m, d) = (L::masked(keep, m), L::masked(keep, d));

                let (Some(mm), Some(ii), Some(dd)) = (
                    buffer
                        .m
                        .get_mut(at..at + LANE_MAX)
                        .and_then(<[f32]>::first_chunk_mut::<LANE_MAX>),
                    buffer
                        .i
                        .get_mut(at..at + LANE_MAX)
                        .and_then(<[f32]>::first_chunk_mut::<LANE_MAX>),
                    buffer
                        .d
                        .get_mut(at..at + LANE_MAX)
                        .and_then(<[f32]>::first_chunk_mut::<LANE_MAX>),
                ) else {
                    return impossible;
                };
                m.store(mm);
                i.store(ii);
                d.store(dd);

                running = running.vmax(m).vmax(i).vmax(d);
                if summing {
                    // `(m + i)` and not `total + m + i`: the strip kernel adds
                    // the pair before the accumulator, and `f32` addition is
                    // not associative, so the parentheses are the parity.
                    total = total + (m + i);
                }
                carry.left_m = m;
                carry.left_d = d;
                carry.diag_m = up_m;
                carry.diag_indel = up_i + up_d;
            }

            // The two cells the next row reads that this one did not write:
            // its own `first - 1`, where the band was clamped at column 1, and
            // `last + 1`, where the band still has room to move right.
            for at in [(first - 1) * BATCH, (last + 1) * BATCH] {
                for track in [&mut buffer.m, &mut buffer.i, &mut buffer.d] {
                    if let Some(slot) = track.get_mut(at..at + LANE_MAX) {
                        slot.fill(0.0);
                    }
                }
            }
        }
        if row % STRIP_ROWS == 0 {
            crossing = running;
        }
    }

    total.store(&mut scratch);
    let mut out = impossible;
    for lane in 0..L::LANES.min(BATCH) {
        let Some(sum) = scratch.get(lane).copied() else {
            continue;
        };
        let Some(shift) = exponent.get(lane).copied() else {
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
