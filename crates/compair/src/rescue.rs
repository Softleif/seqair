//! [`align_banded_f64`]: the `f64` recurrence over a band, two rows at a time.
//!
//! Bit for bit the row-at-a-time recurrence (`align_banded_f64_rows`, which
//! it replaced and is pinned to): every cell is the same IEEE operations on
//! the same values in the same order, with no fused multiply-add, and every
//! row is rescaled by the same power of two.
//!
//! In a row, the match and insertion cells read only the row above, so they
//! are vectors; the deletions are a chain along the row, a multiply and an
//! add per cell whose latency no reordering may shorten. Two rows share one
//! sweep so that their chains overlap: [`prepass`] makes row `a`'s match and
//! insertion cells, and [`sweep`] runs row `a`'s chain and the whole of row
//! `b` a vector behind it.
//!
//! What stands between two rows is the rescale: row `b` reads row `a`
//! multiplied by the power of two that row `a`'s largest cell asks for, and
//! that is not known until row `a`'s chain ends. Row `a`'s largest match or
//! insertion cell is known after the prepass, and its deletions almost never
//! outgrow it, so the sweep takes that scale and the caller checks it
//! against the chain's maximum afterwards; where the deletions did set the
//! scale, row `b` is made again from row `a` at the right one. The rescale
//! itself is applied as a row is read, which is the same multiply the
//! reference stores.
//!
//! A row is never multiplied in place, so each buffer holds one row unscaled
//! and is zero outside its span, as the reference's are.

use fearless_simd::{Level, Select, Simd, SimdBase, dispatch};
use fearless_simd_macros::simd;
use seqair_types::Base;

use crate::{
    banded::Band,
    emission::{Emission, epsilon},
    haplotype::Haplotype,
    read::Read,
    scaling::{exp2_f64, row_shift_f64},
    transitions::Transition,
    types::Log10Likelihood,
};

/// The `f64` recurrence over a band: [`crate::align_full`] restricted to the
/// cells `|j - i - offset| <= width / 2`, which is what every banded kernel
/// computes in `f32`.
///
/// It is what an `f32` score falls back to below the level where `f32` can
/// vouch for it (see `reference::trusted`), and it is exposed for callers that
/// want the banded answer without the `f32` kernels' range. It visits only
/// the band's cells and rescales each row by a power of two, so that its
/// largest cell sits at `2^960` (`scaling::F64_ROW_SCALE`, which also shows
/// nothing overflows).
///
/// **Precision.** A value below `2^-1022` rounds to a multiple of `2^-1074`,
/// losing up to `2^-1075` of it -- in either direction -- and that is the
/// only loss beyond `f64`'s relative rounding, which stays below `1e-12` of
/// the total. With `rho` the free start's mass, and `total`, `spill` and
/// `carry` the read's `Growth` (all one or below with steady qualities; see
/// `reference::trusted` for the same argument about the `f32` kernels):
///
/// - a row's largest cell of the recurrence is at most
///   `max(1, spill) * rho * total`, so the row's scale `2^e >= 2^960 / that`,
///   and a rounding in it is at most `2^(-1075 - e)` of the recurrence's
///   units -- also where a row's lift was capped, which only follows a fall;
/// - at most sixteen roundings per cell (its three values, the products
///   they are sums of, and the rescale), in `r + 1` rows of `columns` cells;
/// - each hands on at most `max(1, carry) * total`.
///
/// So the total is within `16 * (r + 1) * columns * max(1, carry) * max(1,
/// spill) * rho * total * 2^(-1075 - 960)` of the recurrence's, and a total a
/// million times that is right to within a millionth: for 150 bases in a
/// 48-wide band with steady qualities, any score above about `-601`. Below
/// that the answer can be off by up to the bound either way: a total far
/// under it can come out as zero ([`Log10Likelihood::IMPOSSIBLE`]), or above
/// what it should be by up to the bound. At the old scale, `[1, 2)`, the same
/// held only above about `-310`, and a read that switched paths after 78
/// decades more scored -430 where the recurrence gives -352.
pub fn align_banded_f64<E: Emission + ?Sized>(
    haplotype: &Haplotype,
    read: &Read,
    emission: &E,
    band: Band,
) -> Log10Likelihood {
    align_banded_f64_at(Level::new(), haplotype, read, emission, band)
}

/// [`align_banded_f64`] at a given SIMD level, for tests that run every one.
#[doc(hidden)]
pub fn align_banded_f64_at<E: Emission + ?Sized>(
    level: Level,
    haplotype: &Haplotype,
    read: &Read,
    emission: &E,
    band: Band,
) -> Log10Likelihood {
    SCRATCH.with(|scratch| match scratch.try_borrow_mut() {
        Ok(mut scratch) => score(level, haplotype, read, emission, band, &mut scratch),
        // An emission that scores a pair from inside `site_weights`.
        Err(_) => score(level, haplotype, read, emission, band, &mut Scratch::default()),
    })
}

/// Every buffer [`align_banded_f64`] uses, kept per thread so that a call
/// allocates nothing once the thread has scored a pair as long.
#[derive(Default)]
struct Scratch {
    /// Each column's latent weight for every base a read can show
    /// (`hits`), and one minus it (`misses`), one track per base, indexed
    /// like [`Rows`]. Only the columns the band reaches are written, and a
    /// track only if the read shows its base.
    hits: [Vec<f64>; 5],
    misses: [Vec<f64>; 5],
    plans: Vec<RowPlan>,
    above: Rows,
    below: Rows,
    a: Rows,
}

std::thread_local! {
    static SCRATCH: core::cell::RefCell<Scratch> = core::cell::RefCell::default();
}

fn score<E: Emission + ?Sized>(
    level: Level,
    haplotype: &Haplotype,
    read: &Read,
    emission: &E,
    band: Band,
    scratch: &mut Scratch,
) -> Log10Likelihood {
    let (h, r) = (haplotype.len(), read.len());
    if h == 0 || r == 0 {
        return Log10Likelihood::IMPOSSIBLE;
    }
    // The columns of row `i` inside the band and the haplotype, `1..=h` for
    // every row but the start row's `0..=h`, as buffer indices.
    let span = |i: usize, first: usize| -> Option<(usize, usize)> {
        let (i, last) = (i64::try_from(i).ok()?, i64::try_from(h).ok()?);
        let centre = i.checked_add(band.offset)?;
        let low = centre.checked_sub(band.half_width)?.max(i64::try_from(first).ok()?);
        let high = centre.checked_add(band.half_width)?.min(last);
        (low <= high)
            .then_some((usize::try_from(low).ok()? + PAD, usize::try_from(high).ok()? + PAD))
    };
    let Some((first, last)) = band.columns(h, r) else {
        return Log10Likelihood::IMPOSSIBLE;
    };
    let width = h + 1 + 2 * PAD;

    // Every row's plan, up front: a row the reference gives up on ends the
    // pair wherever it is.
    let Scratch { hits, misses, plans, above, below, a } = scratch;
    plans.clear();
    let mut shown = [false; 5];
    for i in 1..=r {
        let (Some(observation), Some(t), Some((low, high))) =
            (read.observation(i - 1), read.transition(i - 1), span(i, 1))
        else {
            return Log10Likelihood::IMPOSSIBLE;
        };
        // The weights the reference slices out must exist.
        if low - PAD - 1 < first || high - PAD - 1 > last {
            return Log10Likelihood::IMPOSSIBLE;
        }
        let eps = epsilon(emission, observation);
        let track = observation.base.known_index().unwrap_or(UNKNOWN_TRACK);
        if let Some(shown) = shown.get_mut(track) {
            *shown = true;
        }
        plans.push(RowPlan { low, high, track, hit: 1.0 - eps, miss: eps / 3.0, t });
    }
    let Some(start) = span(0, 0) else {
        return Log10Likelihood::IMPOSSIBLE;
    };

    let strand = read.strand();
    for track in hits.iter_mut().chain(misses.iter_mut()) {
        // Past what the band reaches a track is only read into lanes that
        // are masked off, so what is left there does not matter.
        track.resize(width, 0.0);
    }
    for index in first..=last {
        let Some(site) = haplotype.site(index) else {
            return Log10Likelihood::IMPOSSIBLE;
        };
        let weights = emission.site_weights(site, strand);
        for (((hit, miss), base), shown) in
            hits.iter_mut().zip(misses.iter_mut()).zip(OBSERVABLE).zip(shown)
        {
            if !shown {
                continue;
            }
            let weight = weights.latent_weight(base);
            if let (Some(hit), Some(miss)) =
                (hit.get_mut(index + 1 + PAD), miss.get_mut(index + 1 + PAD))
            {
                (*hit, *miss) = (weight, 1.0 - weight);
            }
        }
    }
    for rows in [&mut *above, &mut *below, &mut *a] {
        rows.prepare(width);
    }
    // The start row: zero everywhere the first row reads it, but for the
    // free start in its deletions.
    for buffer in [&mut above.m, &mut above.i, &mut above.d] {
        if let Some(cells) = buffer.get_mut(..(start.1 + PAD).min(width)) {
            cells.fill(0.0);
        }
    }
    #[allow(
        clippy::cast_precision_loss,
        reason = "a haplotype long enough to lose precision here does not exist"
    )]
    let init = 1.0 / h as f64;
    if let Some(cells) = above.d.get_mut(start.0..=start.1) {
        cells.fill(init);
    }
    let tracks = Tracks { hits, misses };
    // The start row is stored unscaled and read at the row scale, as every
    // row is.
    let start = (start, row_shift_f64(init));
    dispatch!(level, simd => kernel_simd(simd, plans, tracks, start, [above, below, a]))
}

/// Zero columns either side of `0..=h` in every buffer, so a vector may read
/// and write past the band's edge without a check: a row's loads reach two
/// vectors past its span. Twice the widest vector's lane count.
const PAD: usize = 16;

/// The bases a read can show, in [`Base::known_index`] order and then
/// [`Base::Unknown`].
const OBSERVABLE: [Base; 5] = [Base::A, Base::C, Base::G, Base::T, Base::Unknown];
const UNKNOWN_TRACK: usize = 4;

/// The larger of two cells, neither a NaN.
#[inline(always)]
fn larger(a: f64, b: f64) -> f64 {
    if b > a { b } else { a }
}

/// One row's three matrices, unscaled, indexed by column + [`PAD`].
///
/// What the next row reads of it is its span, the match and insertion cells
/// one column right of it -- zero: every writer stores a zero vector past
/// its last block -- and, where the band starts at column 1, column 0, which
/// no row but the start row writes. The rest holds whatever an earlier row
/// or an earlier call left, and reaches only lanes that are masked off; the
/// buffers are never cleared.
#[derive(Default)]
struct Rows {
    m: Vec<f64>,
    i: Vec<f64>,
    d: Vec<f64>,
}

impl Rows {
    /// `width` cells in every matrix, column 0 zero.
    fn prepare(&mut self, width: usize) {
        for buffer in [&mut self.m, &mut self.i, &mut self.d] {
            buffer.resize(width, 0.0);
            if let Some(cell) = buffer.get_mut(PAD) {
                *cell = 0.0;
            }
        }
    }
}

/// Fills every buffer of this thread's [`align_banded_f64`] with NaN, for
/// tests that prove no call reads what an earlier one left.
#[doc(hidden)]
pub fn poison_align_banded_f64_scratch() {
    SCRATCH.with(|scratch| {
        if let Ok(mut scratch) = scratch.try_borrow_mut() {
            let Scratch { hits, misses, plans: _, above, below, a } = &mut *scratch;
            let rows = [above, below, a]
                .into_iter()
                .flat_map(|rows| [&mut rows.m, &mut rows.i, &mut rows.d]);
            for buffer in hits.iter_mut().chain(misses.iter_mut()).chain(rows) {
                buffer.clear();
                buffer.resize(4096, f64::NAN);
            }
        }
    });
}

/// Each column's latent weight for every base (`hits`) and one minus it
/// (`misses`), indexed like [`Rows`].
#[derive(Clone, Copy)]
struct Tracks<'a> {
    hits: &'a [Vec<f64>; 5],
    misses: &'a [Vec<f64>; 5],
}

/// One read row: its span as buffer indices, its base's weight track, its
/// emission terms and its transitions.
#[derive(Clone, Copy)]
struct RowPlan {
    low: usize,
    high: usize,
    track: usize,
    hit: f64,
    miss: f64,
    t: Transition,
}

type V<S> = <S as Simd>::f64s;

/// A row's transitions and emission terms as vectors.
#[derive(Clone, Copy)]
struct Terms<S: Simd> {
    hit: V<S>,
    miss: V<S>,
    match_to_match: V<S>,
    indel_to_match: V<S>,
    match_to_insertion: V<S>,
    gap_continuation: V<S>,
}

impl<S: Simd> Terms<S> {
    #[inline(always)]
    fn new(simd: S, plan: &RowPlan) -> Self {
        let t = &plan.t;
        Self {
            hit: V::<S>::splat(simd, plan.hit),
            miss: V::<S>::splat(simd, plan.miss),
            match_to_match: V::<S>::splat(simd, t.match_to_match),
            indel_to_match: V::<S>::splat(simd, t.indel_to_match),
            match_to_insertion: V::<S>::splat(simd, t.match_to_insertion),
            gap_continuation: V::<S>::splat(simd, t.gap_continuation),
        }
    }

    /// The match and insertion cells of `LANES` columns from their
    /// neighbours in the row above (already scaled) and their weights: the
    /// reference's operations, in its order.
    #[inline(always)]
    #[allow(clippy::too_many_arguments, reason = "the five neighbours and two weights of a cell")]
    fn cells(
        &self,
        m_diagonal: V<S>,
        i_diagonal: V<S>,
        d_diagonal: V<S>,
        m_up: V<S>,
        i_up: V<S>,
        hit: V<S>,
        miss: V<S>,
    ) -> (V<S>, V<S>) {
        let prior = hit * self.hit + miss * self.miss;
        let m = prior
            * (m_diagonal * self.match_to_match
                + i_diagonal * self.indel_to_match
                + d_diagonal * self.indel_to_match);
        let i = m_up * self.match_to_insertion + i_up * self.gap_continuation;
        (m, i)
    }
}

/// Lanes `0..count` of `value`, zero in the rest.
#[inline(always)]
fn first_lanes<S: Simd>(simd: S, count: usize, value: V<S>) -> V<S> {
    let lanes = lane_indices(simd);
    #[allow(clippy::cast_precision_loss, reason = "a lane count")]
    let count = V::<S>::splat(simd, count as f64);
    lanes.simd_lt(count).select(value, V::<S>::splat(simd, 0.0))
}

#[simd]
fn kernel_simd<S: Simd>(
    simd: S,
    plans: &[RowPlan],
    tracks: Tracks<'_>,
    start: ((usize, usize), i32),
    rows: [&mut Rows; 3],
) -> Log10Likelihood {
    kernel(simd, plans, tracks, start, rows)
}

/// The rows two at a time: [`prepass`] for the first row's match and
/// insertion cells, then [`sweep`] for its deletion chain and the whole
/// second row.
#[inline(always)]
fn kernel<S: Simd>(
    simd: S,
    plans: &[RowPlan],
    tracks: Tracks<'_>,
    (start, start_shift): ((usize, usize), i32),
    [above, below, a]: [&mut Rows; 3],
) -> Log10Likelihood {
    // `above` holds the start row.
    let (mut above, mut below, mut a) = (above, below, a);
    let (mut above_span, mut above_scale) = (start, exp2_f64(start_shift));
    let mut exponent = start_shift;

    let mut pairs = plans.chunks(2);
    for pair in &mut pairs {
        let (pa, pb) = match pair {
            [pa, pb] => (pa, Some(pb)),
            [pa] => (pa, None),
            _ => break,
        };
        let max_a = prepass(simd, above, above_scale, pa, tracks, a);
        // Nothing reads `above` again but at column 0, where only the
        // start row has a cell.
        if let Some(cell) = above.d.get_mut(PAD) {
            *cell = 0.0;
        }
        let guess = row_shift_f64(max_a);
        let Some(pb) = pb else {
            // The last row, on its own.
            let max_d = chain(a, pa);
            let shift = row_shift_f64(larger(max_a, max_d));
            exponent += shift;
            (above_span, above_scale) = ((pa.low, pa.high), exp2_f64(shift));
            core::mem::swap(&mut above, &mut a);
            break;
        };
        let fused = if pb.low == pa.low + 1 {
            let (max_da, max_b) = sweep(simd, a, exp2_f64(guess), pa, pb, tracks, below);
            let shift_a = row_shift_f64(larger(max_a, max_da));
            (shift_a == guess).then_some((shift_a, max_b))
        } else {
            None
        };
        // Otherwise row `a` first, then row `b` from it: where the band
        // starts at column 1, or where row `a`'s deletions set its scale.
        let (shift_a, max_b) = if let Some(fused) = fused {
            fused
        } else {
            let shift_a = row_shift_f64(larger(max_a, chain(a, pa)));
            let max_b = prepass(simd, a, exp2_f64(shift_a), pb, tracks, below);
            (shift_a, larger(max_b, chain(below, pb)))
        };
        let shift_b = row_shift_f64(max_b);
        exponent += shift_a + shift_b;
        core::mem::swap(&mut above, &mut below);
        (above_span, above_scale) = ((pb.low, pb.high), exp2_f64(shift_b));
    }

    let (m, i) =
        (above.m.get(above_span.0..=above_span.1), above.i.get(above_span.0..=above_span.1));
    let total: f64 = m
        .into_iter()
        .flatten()
        .zip(i.into_iter().flatten())
        .map(|(m, i)| m * above_scale + i * above_scale)
        .sum();
    if total <= 0.0 {
        return Log10Likelihood::IMPOSSIBLE;
    }
    Log10Likelihood::new(total.log10() - f64::from(exponent) * core::f64::consts::LOG10_2)
}

/// `LANES` cells from a slice of exactly that many.
#[inline(always)]
fn load<S: Simd>(simd: S, cells: &[f64]) -> V<S> {
    V::<S>::from_slice(simd, cells)
}

/// A row's match and insertion cells from the row above, which is scaled by
/// `scale` on load, and the largest of them. The deletions are [`chain`]'s.
/// Writes zeros past the row's span, up to a whole vector.
///
/// Every load starts a vector where the row above stored one (its span
/// starts a column left of this row's): a load straddling two stores still
/// in flight cannot be forwarded and waits for both to reach the cache. The
/// cells straight above are the diagonal ones slid by a lane.
#[inline(always)]
fn prepass<S: Simd>(
    simd: S,
    above: &Rows,
    scale: f64,
    plan: &RowPlan,
    tracks: Tracks<'_>,
    out: &mut Rows,
) -> f64 {
    let lanes = V::<S>::LEN;
    let (low, high) = (plan.low, plan.high);
    let cells = high + 1 - low;
    let (full, tail) = (cells / lanes, cells % lanes);
    let n = cells.div_ceil(lanes) * lanes;
    // One vector more on the diagonal, for the last block's cells above.
    let (Some(m_d), Some(i_d), Some(d_d)) = (
        above.m.get(low - 1..low - 1 + n + lanes),
        above.i.get(low - 1..low - 1 + n + lanes),
        above.d.get(low - 1..low - 1 + n),
    ) else {
        return 0.0;
    };
    let (Some(hits), Some(misses)) = (
        tracks.hits.get(plan.track).and_then(|track| track.get(low..low + n)),
        tracks.misses.get(plan.track).and_then(|track| track.get(low..low + n)),
    ) else {
        return 0.0;
    };
    let (Some(out_m), Some(out_i)) =
        (out.m.get_mut(low..low + n + lanes), out.i.get_mut(low..low + n + lanes))
    else {
        return 0.0;
    };
    let (out_m, past_m) = out_m.split_at_mut(n);
    let (out_i, past_i) = out_i.split_at_mut(n);
    let zero = V::<S>::splat(simd, 0.0);
    zero.store_slice(past_m);
    zero.store_slice(past_i);
    let (m_first, m_d) = m_d.split_at(lanes);
    let (i_first, i_d) = i_d.split_at(lanes);
    let scale = V::<S>::splat(simd, scale);
    let mut block = PrepassBlock {
        simd,
        terms: Terms::new(simd, plan),
        scale,
        m_diagonal: load(simd, m_first) * scale,
        i_diagonal: load(simd, i_first) * scale,
        max_m: V::<S>::splat(simd, 0.0),
        max_i: V::<S>::splat(simd, 0.0),
    };
    let next = m_d.chunks_exact(lanes).zip(i_d.chunks_exact(lanes));
    let weights = hits.chunks_exact(lanes).zip(misses.chunks_exact(lanes));
    let outputs = out_m.chunks_exact_mut(lanes).zip(out_i.chunks_exact_mut(lanes));
    let mut blocks = next.zip(d_d.chunks_exact(lanes)).zip(weights).zip(outputs);
    for (((next, d_d), weights), outputs) in blocks.by_ref().take(full) {
        block.step(next, d_d, weights, outputs, None);
    }
    if let Some((((next, d_d), weights), outputs)) = blocks.next() {
        block.step(next, d_d, weights, outputs, Some(tail));
    }
    block.max_m.max(block.max_i).reduce_max()
}

type Pair<'a> = (&'a [f64], &'a [f64]);
type PairMut<'a> = (&'a mut [f64], &'a mut [f64]);

/// [`prepass`]'s loop state: the block's diagonal neighbours, scaled.
struct PrepassBlock<S: Simd> {
    simd: S,
    terms: Terms<S>,
    scale: V<S>,
    m_diagonal: V<S>,
    i_diagonal: V<S>,
    max_m: V<S>,
    max_i: V<S>,
}

impl<S: Simd> PrepassBlock<S> {
    /// One vector of cells, given the next block's diagonal neighbours;
    /// with `live`, only its first `live` lanes are in the span.
    #[inline(always)]
    fn step(
        &mut self,
        (m_next, i_next): Pair<'_>,
        d_d: &[f64],
        (hit, miss): Pair<'_>,
        (out_m, out_i): PairMut<'_>,
        live: Option<usize>,
    ) {
        let simd = self.simd;
        let (m_next, i_next) = (load(simd, m_next) * self.scale, load(simd, i_next) * self.scale);
        let (mut m, mut i) = self.terms.cells(
            self.m_diagonal,
            self.i_diagonal,
            load(simd, d_d) * self.scale,
            self.m_diagonal.slide::<1>(m_next),
            self.i_diagonal.slide::<1>(i_next),
            load(simd, hit),
            load(simd, miss),
        );
        (self.m_diagonal, self.i_diagonal) = (m_next, i_next);
        if let Some(live) = live {
            (m, i) = (first_lanes(simd, live, m), first_lanes(simd, live, i));
        }
        m.store_slice(out_m);
        i.store_slice(out_i);
        (self.max_m, self.max_i) = (self.max_m.max(m), self.max_i.max(i));
    }
}

/// `0, 1, .., LANES - 1`.
#[inline(always)]
fn lane_indices<S: Simd>(simd: S) -> V<S> {
    #[allow(clippy::cast_precision_loss, reason = "a lane index")]
    V::<S>::from_fn(simd, |lane| lane as f64)
}

/// A row's deletion chain from its match cells, stored, and the largest.
#[inline(always)]
fn chain(rows: &mut Rows, plan: &RowPlan) -> f64 {
    let t = &plan.t;
    let (mut d, mut max) = (0.0f64, 0.0f64);
    let (Some(m), Some(out)) =
        (rows.m.get(plan.low..plan.high), rows.d.get_mut(plan.low..=plan.high))
    else {
        return 0.0;
    };
    for (m_left, cell) in core::iter::once(&0.0).chain(m).zip(out) {
        d = m_left * t.match_to_deletion + d * t.gap_continuation;
        *cell = d;
        max = larger(max, d);
    }
    max
}

/// Row `a`'s deletion chain and row `b` (the row below it), swept along the
/// band together, so that the two deletion chains -- each a multiply and an
/// add per cell, whose latency is the whole cost of a row -- overlap. Row
/// `b`'s span must start one column right of row `a`'s.
///
/// Row `a`'s match and insertion cells are [`prepass`]'s, and `scale` is the
/// scale their maximum asks for. Row `b` reads row `a` at that scale, which is
/// the reference's only if row `a`'s deletions do not ask for another: the
/// caller checks, and redoes row `b` if they did. Returns the largest of row
/// `a`'s deletions, and the largest cell of row `b`.
///
/// A block of row `b`'s cells at columns `c..c + LANES` reads row `a`'s
/// deletions at `c - 1..c + LANES - 1`, which is the block of the chain the
/// same step makes. Loads line up with [`prepass`]'s stores, as its own do
/// with this one's.
#[inline(always)]
fn sweep<S: Simd>(
    simd: S,
    a: &Rows,
    scale: f64,
    pa: &RowPlan,
    pb: &RowPlan,
    tracks: Tracks<'_>,
    out: &mut Rows,
) -> (f64, f64) {
    let lanes = V::<S>::LEN;
    let (low, high) = (pb.low, pb.high);
    debug_assert_eq!(low, pa.low + 1);
    let cells = high + 1 - low;
    let (full, tail) = (cells / lanes, cells % lanes);
    let n = cells.div_ceil(lanes) * lanes;
    let (Some(m_d), Some(i_d)) =
        (a.m.get(low - 1..low - 1 + n + lanes), a.i.get(low - 1..low - 1 + n + lanes))
    else {
        return (0.0, 0.0);
    };
    let (Some(hits), Some(misses)) = (
        tracks.hits.get(pb.track).and_then(|track| track.get(low..low + n)),
        tracks.misses.get(pb.track).and_then(|track| track.get(low..low + n)),
    ) else {
        return (0.0, 0.0);
    };
    let (Some(out_m), Some(out_i), Some(out_d)) = (
        out.m.get_mut(low..low + n + lanes),
        out.i.get_mut(low..low + n + lanes),
        out.d.get_mut(low..low + n),
    ) else {
        return (0.0, 0.0);
    };
    let (out_m, past_m) = out_m.split_at_mut(n);
    let (out_i, past_i) = out_i.split_at_mut(n);
    let zero = V::<S>::splat(simd, 0.0);
    zero.store_slice(past_m);
    zero.store_slice(past_i);
    let (m_first, m_d) = m_d.split_at(lanes);
    let (i_first, i_d) = i_d.split_at(lanes);
    let scale = V::<S>::splat(simd, scale);
    let m_first = load(simd, m_first);
    let mut block = SweepBlock {
        simd,
        terms: Terms::new(simd, pb),
        ta: pa.t,
        tb: pb.t,
        scale,
        m_raw: m_first,
        m_diagonal: m_first * scale,
        i_diagonal: load(simd, i_first) * scale,
        m_a: 0.0,
        d_a: 0.0,
        m_b: 0.0,
        d_b: 0.0,
        max_da: zero,
        max_m: zero,
        max_i: zero,
        max_db: zero,
    };
    let next = m_d.chunks_exact(lanes).zip(i_d.chunks_exact(lanes));
    let weights = hits.chunks_exact(lanes).zip(misses.chunks_exact(lanes));
    let outputs = out_m
        .chunks_exact_mut(lanes)
        .zip(out_i.chunks_exact_mut(lanes))
        .zip(out_d.chunks_exact_mut(lanes));
    let mut blocks = next.zip(weights).zip(outputs);
    for ((next, weights), outputs) in blocks.by_ref().take(full) {
        block.step(next, weights, outputs, None);
    }
    if let Some(((next, weights), outputs)) = blocks.next() {
        // Past row `a`'s span its chain is zero, past row `b`'s everything.
        let live_a = pa.high + 2 - (low + full * lanes);
        block.step(next, weights, outputs, Some((live_a, tail)));
    }
    let mut max_da = block.max_da.reduce_max();
    if tail == 0 && pa.high == high {
        // The chain's last cell, which no block of row `b` reads.
        let m_left = a.m.get(high - 1).copied().unwrap_or(0.0);
        max_da =
            larger(max_da, m_left * pa.t.match_to_deletion + block.d_a * pa.t.gap_continuation);
    }
    (max_da, block.max_m.max(block.max_i).max(block.max_db).reduce_max())
}

type Triple<'a> = ((&'a mut [f64], &'a mut [f64]), &'a mut [f64]);

/// [`sweep`]'s loop state: the block's diagonal neighbours in row `a`, as
/// stored and scaled, each chain's last match cell and deletion, and the
/// running maxima.
struct SweepBlock<S: Simd> {
    simd: S,
    terms: Terms<S>,
    ta: Transition,
    tb: Transition,
    scale: V<S>,
    m_raw: V<S>,
    m_diagonal: V<S>,
    i_diagonal: V<S>,
    m_a: f64,
    d_a: f64,
    m_b: f64,
    d_b: f64,
    max_da: V<S>,
    max_m: V<S>,
    max_i: V<S>,
    max_db: V<S>,
}

impl<S: Simd> SweepBlock<S> {
    /// One vector of row `b`'s cells and the block of row `a`'s chain they
    /// read, given the next block's diagonal neighbours; with `live`, only
    /// the first lanes of each are in its row's span.
    #[inline(always)]
    fn step(
        &mut self,
        (m_next, i_next): Pair<'_>,
        (hit, miss): Pair<'_>,
        ((out_m, out_i), out_d): Triple<'_>,
        live: Option<(usize, usize)>,
    ) {
        let simd = self.simd;
        let (ta, tb) = (&self.ta, &self.tb);
        // Row `a`'s deletions at `c - 1..c + LANES - 1`.
        let mut d_diagonal = V::<S>::splat(simd, 0.0);
        for (cell, &m_left) in d_diagonal
            .as_mut_slice()
            .iter_mut()
            .zip(core::iter::once(&self.m_a).chain(self.m_raw.as_slice()))
        {
            self.d_a = m_left * ta.match_to_deletion + self.d_a * ta.gap_continuation;
            *cell = self.d_a;
        }
        self.m_a = self.m_raw.as_slice().last().copied().unwrap_or(0.0);
        if let Some((live_a, _)) = live {
            d_diagonal = first_lanes(simd, live_a, d_diagonal);
        }
        self.max_da = self.max_da.max(d_diagonal);
        let m_raw = load(simd, m_next);
        let (m_next, i_next) = (m_raw * self.scale, load(simd, i_next) * self.scale);
        let (mut m, mut i) = self.terms.cells(
            self.m_diagonal,
            self.i_diagonal,
            d_diagonal * self.scale,
            self.m_diagonal.slide::<1>(m_next),
            self.i_diagonal.slide::<1>(i_next),
            load(simd, hit),
            load(simd, miss),
        );
        (self.m_raw, self.m_diagonal, self.i_diagonal) = (m_raw, m_next, i_next);
        if let Some((_, live_b)) = live {
            (m, i) = (first_lanes(simd, live_b, m), first_lanes(simd, live_b, i));
        }
        m.store_slice(out_m);
        i.store_slice(out_i);
        (self.max_m, self.max_i) = (self.max_m.max(m), self.max_i.max(i));
        let mut d = V::<S>::splat(simd, 0.0);
        for (cell, &m_left) in
            d.as_mut_slice().iter_mut().zip(core::iter::once(&self.m_b).chain(m.as_slice()))
        {
            self.d_b = m_left * tb.match_to_deletion + self.d_b * tb.gap_continuation;
            *cell = self.d_b;
        }
        self.m_b = m.as_slice().last().copied().unwrap_or(0.0);
        if let Some((_, live_b)) = live {
            d = first_lanes(simd, live_b, d);
        }
        d.store_slice(out_d);
        self.max_db = self.max_db.max(d);
    }
}
