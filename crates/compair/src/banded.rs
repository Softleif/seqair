use wide::{Select, f32x8};

use crate::{
    emission::Emission,
    error::Error,
    haplotype::Haplotype,
    read::Read,
    scaling::{exp2_f32, exp2_f64, normalising_shift_f32},
    transitions::transitions,
    types::Log10Likelihood,
};
use seqair_types::Base;

/// The diagonal strip of the matrix the banded kernels visit.
///
/// A cell `(i, j)` is inside the band exactly when `|j - i - offset| <= width /
/// 2`, so `offset` is the haplotype position the read's first base is expected
/// to sit near. `width` is rounded **down** to an even number, and the strip is
/// then `2 * (width / 2) + 1` columns wide.
///
/// Two things about a band that does not hold the optimal alignment:
///
/// - if some path still crosses the strip, the score is that path's: finite,
///   lower than the unbanded one, and never an error;
/// - if none does, the score is [`Log10Likelihood::IMPOSSIBLE`]. The read's
///   free start lives on the cells `(0, j)`, so an `offset` below `-(width /
///   2)` puts every one of them outside the strip and no path opens at all.
///
/// Either way the strip is the only thing either kernel reads, so two
/// haplotypes scored under one band are compared like with like.
///
/// [`Log10Likelihood::IMPOSSIBLE`]: crate::Log10Likelihood::IMPOSSIBLE
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(
    clippy::cast_lossless,
    clippy::cast_possible_wrap,
    reason = "`anchored` is const, so it widens with `as` rather than `from`"
)]
pub struct Band {
    half_width: i64,
    offset: i64,
}

impl Band {
    /// The width §6.4 of the design budgets for: 150 bp of read against 48
    /// columns is ~7,200 cells.
    pub const DEFAULT_WIDTH: u32 = 48;

    /// The widest band either kernel will build, for two reasons that are the
    /// same reason.
    ///
    /// The banded kernels allocate nine buffers of `width / 2 + O(1)` floats,
    /// so an unbounded `width` is an unbounded allocation inside a safe
    /// function — `u32::MAX` would ask for ~77 GB.
    ///
    /// And a wide band is not a better answer, it is a worse one. The `f32`
    /// kernels renormalise per **anti-diagonal**, and a diagonal holds one cell
    /// per read row it crosses — prefix alignments of that many lengths, whose
    /// magnitudes span the whole alignment's dynamic range. The band is what
    /// bounds that span: at `DEFAULT_WIDTH` the `f32` kernel reproduces the
    /// `f64` recurrence over the same band to rounding even at
    /// `log10 L = -200`, while a band wider than the read loses several log10
    /// to underflow. Use [`align_full`] when the band is not wanted, not a
    /// band wide enough to disable itself.
    ///
    /// [`align_full`]: crate::align_full
    pub const MAX_WIDTH: u32 = 1024;

    pub fn new(width: u32, offset: i32) -> Result<Self, Error> {
        if width < 2 {
            return Err(Error::BandTooNarrow(width));
        }
        if width > Self::MAX_WIDTH {
            return Err(Error::BandTooWide { got: width, max: Self::MAX_WIDTH });
        }
        Ok(Self { half_width: i64::from(width / 2), offset: i64::from(offset) })
    }

    /// `DEFAULT_WIDTH` columns centred on the given haplotype offset.
    #[must_use]
    pub const fn anchored(offset: i32) -> Self {
        Self { half_width: (Self::DEFAULT_WIDTH / 2) as i64, offset: offset as i64 }
    }

    #[must_use]
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        clippy::cast_lossless,
        reason = "both fields were built from a u32 and an i32 in `new`"
    )]
    pub const fn width(self) -> u32 {
        (self.half_width * 2) as u32
    }

    #[must_use]
    #[allow(clippy::cast_possible_truncation, reason = "the offset was built from an i32 in `new`")]
    pub const fn offset(self) -> i32 {
        self.offset as i32
    }
}

/// One scalar or one vector of the DP's cells.
///
/// The banded kernel is written once against this trait, so the scalar and the
/// eight-wide kernels perform the same operations in the same order on the same
/// values; bit-parity is then a property of the code, not a test result that
/// could quietly stop holding.
///
/// There is deliberately no fused multiply-add here. `f32::mul_add` fuses on
/// every target, while `wide`'s vector one falls back to a separate multiply
/// and add wherever the vector unit has none -- so a `mul_add` on this trait
/// would round identically on aarch64 and differently on a baseline x86-64,
/// which is the one kind of parity failure the trait exists to rule out.
trait Lane: Copy + core::ops::Add<Output = Self> + core::ops::Mul<Output = Self> {
    const LANES: usize;
    fn splat(value: f32) -> Self;
    fn load(source: &Window) -> Self;
    fn store(self, destination: &mut Window);
    fn vmax(self, other: Self) -> Self;
    fn horizontal_max(self) -> f32;
    /// Every bit set where the two lanes are equal, no bit set where they are
    /// not: a mask for [`Lane::select`], not a number.
    fn equals(self, other: Self) -> Self;
    /// Bitwise or of two masks.
    fn either(self, other: Self) -> Self;
    /// Every bit set where `self` is strictly below `other`.
    fn below(self, other: Self) -> Self;
    /// `0.0, 1.0, ...` up to [`Lane::LANES`], for turning a live-cell count
    /// into a mask.
    fn offsets() -> Self;
    /// `self` is a mask from [`Lane::equals`] or [`Lane::either`].
    fn select(self, if_true: Self, if_false: Self) -> Self;
}

/// What a lanewise comparison reports for "equal". As a number it is a `NaN`;
/// it is only ever consumed bitwise.
const MASK_SET: f32 = f32::from_bits(u32::MAX);

impl Lane for f32 {
    const LANES: usize = 1;
    #[inline]
    fn splat(value: f32) -> Self {
        value
    }
    #[inline]
    fn load(source: &Window) -> Self {
        source[0]
    }
    #[inline]
    fn store(self, destination: &mut Window) {
        destination[0] = self;
    }
    #[inline]
    fn vmax(self, other: Self) -> Self {
        if self > other { self } else { other }
    }
    #[inline]
    fn horizontal_max(self) -> f32 {
        self
    }
    #[inline]
    fn equals(self, other: Self) -> Self {
        if self == other { MASK_SET } else { 0.0 }
    }
    #[inline]
    fn either(self, other: Self) -> Self {
        Self::from_bits(self.to_bits() | other.to_bits())
    }
    #[inline]
    fn below(self, other: Self) -> Self {
        if self < other { MASK_SET } else { 0.0 }
    }
    #[inline]
    fn offsets() -> Self {
        0.0
    }
    #[inline]
    fn select(self, if_true: Self, if_false: Self) -> Self {
        if self.to_bits() == 0 { if_false } else { if_true }
    }
}

impl Lane for f32x8 {
    const LANES: usize = 8;
    #[inline]
    fn splat(value: f32) -> Self {
        Self::splat(value)
    }
    #[inline]
    fn load(source: &Window) -> Self {
        Self::from(*source)
    }
    #[inline]
    fn store(self, destination: &mut Window) {
        *destination = self.to_array();
    }
    #[inline]
    fn vmax(self, other: Self) -> Self {
        self.max(other)
    }
    #[inline]
    fn horizontal_max(self) -> f32 {
        self.to_array().into_iter().fold(f32::NEG_INFINITY, f32::max)
    }
    #[inline]
    fn equals(self, other: Self) -> Self {
        self.simd_eq(other)
    }
    #[inline]
    fn either(self, other: Self) -> Self {
        self | other
    }
    #[inline]
    fn below(self, other: Self) -> Self {
        self.simd_lt(other)
    }
    #[inline]
    fn offsets() -> Self {
        Self::from([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
    }
    #[inline]
    fn select(self, if_true: Self, if_false: Self) -> Self {
        Select::select(self, if_true, if_false)
    }
}

/// Scalar `f32` over the band.
pub fn align_banded<E: Emission>(
    haplotype: &Haplotype,
    read: &Read,
    emission: &E,
    band: Band,
) -> Log10Likelihood {
    banded_kernel::<f32, E>(haplotype, read, emission, band)
}

/// Eight `f32` lanes of one anti-diagonal at a time.
///
/// Bit-identical to [`align_banded`]; see the `Lane` trait for why.
pub fn align_banded_simd<E: Emission>(
    haplotype: &Haplotype,
    read: &Read,
    emission: &E,
    band: Band,
) -> Log10Likelihood {
    banded_kernel::<f32x8, E>(haplotype, read, emission, band)
}

/// A base as a small integer in `f32`, so a whole lane of them is one vector
/// load and one lanewise compare.
///
/// An unknown base has its own code so that it never compares equal to a real
/// base, and two further codes -- impossible as a read base, so never equal to
/// anything -- stand for "this column does not convert" and "this column has no
/// plain match row".
const CODE_N: f32 = 4.0;
const CODE_NO_CONVERSION: f32 = 6.0;
const CODE_NO_PLAIN_MATCH: f32 = 7.0;

/// The code of each of [`Base::KNOWN`], in that order.
const CODE_KNOWN: [f32; 4] = [0.0, 1.0, 2.0, 3.0];

#[inline]
fn code(base: Base) -> f32 {
    base.known_index().and_then(|index| CODE_KNOWN.get(index).copied()).unwrap_or(CODE_N)
}

/// Everything the inner loop reads, hoisted out of it and split into tracks so
/// that a chunk of lanes loads each one contiguously.
///
/// The emission is a function of `(site, read base, eps)` and the kernel visits
/// every cell, so evaluating the `f64` trait per cell costs `h * r` calls for
/// `h + r` distinct inputs. [`Emission::site_weights`] and [`Emission::epsilon`]
/// are the split that makes the per-cell work a handful of lanewise selects.
///
/// The row tracks are indexed by read row, which is the axis an anti-diagonal
/// runs along. The column tracks are **reversed** -- index `h - j` holds the
/// column a cell at haplotype position `j` reads -- because `j` runs down an
/// anti-diagonal while `i` runs up, and reversing is what makes the haplotype a
/// second ascending contiguous load rather than a gather. Index `h` is the
/// `j == 0` boundary column, which has no site; the kernel zeroes that cell
/// after computing it, so the padding there only has to be finite.
///
/// Every track carries a vector of slack past its live range so that a load at
/// the last live index is still in bounds.
struct Plan {
    row_base: Vec<f32>,
    row_matched: Vec<f32>,
    row_mismatched: Vec<f32>,
    match_to_match: Vec<f32>,
    match_to_insertion: Vec<f32>,
    match_to_deletion: Vec<f32>,
    indel_to_match: Vec<f32>,
    gap_continuation: Vec<f32>,
    column_base: Vec<f32>,
    column_converted: Vec<f32>,
    column_plain: Vec<f32>,
    column_rate: Vec<f32>,
    column_unconverted: Vec<f32>,
}

impl Plan {
    #[allow(
        clippy::cast_possible_truncation,
        reason = "the f32 narrowing is the point of this kernel"
    )]
    fn build<E: Emission>(haplotype: &Haplotype, read: &Read, emission: &E) -> Option<Self> {
        let (h, r) = (haplotype.len(), read.len());
        let rows = r + 1 + LANE_MAX;
        let columns = h + 1 + LANE_MAX;
        let mut plan = Self {
            row_base: vec![0.0; rows],
            row_matched: vec![0.0; rows],
            row_mismatched: vec![0.0; rows],
            match_to_match: vec![0.0; rows],
            match_to_insertion: vec![0.0; rows],
            match_to_deletion: vec![0.0; rows],
            indel_to_match: vec![0.0; rows],
            gap_continuation: vec![0.0; rows],
            column_base: vec![CODE_NO_PLAIN_MATCH; columns],
            column_converted: vec![CODE_NO_CONVERSION; columns],
            column_plain: vec![CODE_NO_PLAIN_MATCH; columns],
            column_rate: vec![0.0; columns],
            column_unconverted: vec![0.0; columns],
        };
        let transitions = transitions(read);
        for index in 0..r {
            let observation = read.observation(index)?;
            let eps = emission.epsilon(observation);
            *plan.row_base.get_mut(index + 1)? = code(observation.base);
            *plan.row_matched.get_mut(index + 1)? = (1.0 - eps) as f32;
            *plan.row_mismatched.get_mut(index + 1)? = (eps / 3.0) as f32;
            let t = transitions.get(index + 1)?;
            *plan.match_to_match.get_mut(index + 1)? = t.match_to_match as f32;
            *plan.match_to_insertion.get_mut(index + 1)? = t.match_to_insertion as f32;
            *plan.match_to_deletion.get_mut(index + 1)? = t.match_to_deletion as f32;
            *plan.indel_to_match.get_mut(index + 1)? = t.indel_to_match as f32;
            *plan.gap_continuation.get_mut(index + 1)? = t.gap_continuation as f32;
        }
        let strand = read.strand();
        for index in 0..h {
            let weights = emission.site_weights(haplotype.site(index)?, strand);
            let base = code(weights.base);
            let reversed = h - 1 - index;
            *plan.column_base.get_mut(reversed)? = base;
            match weights.converted {
                Some(converted) => {
                    *plan.column_converted.get_mut(reversed)? = code(converted);
                    *plan.column_rate.get_mut(reversed)? = weights.rate as f32;
                    *plan.column_unconverted.get_mut(reversed)? = (1.0 - weights.rate) as f32;
                }
                None => *plan.column_plain.get_mut(reversed)? = base,
            }
        }
        Some(plan)
    }
}

/// One haplotype column per lane, in the reversed order the diagonal reads.
#[derive(Clone, Copy)]
struct ColumnLanes<L> {
    base: L,
    converted: L,
    plain: L,
    rate: L,
    unconverted_rate: L,
}

/// One read row per lane.
#[derive(Clone, Copy)]
struct RowLanes<L> {
    base: L,
    matched: L,
    mismatched: L,
}

impl<L: Lane> ColumnLanes<L> {
    #[inline]
    fn load(plan: &Plan, index: usize) -> Self {
        Self {
            base: L::load(window(&plan.column_base, index)),
            converted: L::load(window(&plan.column_converted, index)),
            plain: L::load(window(&plan.column_plain, index)),
            rate: L::load(window(&plan.column_rate, index)),
            unconverted_rate: L::load(window(&plan.column_unconverted, index)),
        }
    }
}

impl<L: Lane> RowLanes<L> {
    #[inline]
    fn load(plan: &Plan, index: usize) -> Self {
        Self {
            base: L::load(window(&plan.row_base, index)),
            matched: L::load(window(&plan.row_matched, index)),
            mismatched: L::load(window(&plan.row_mismatched, index)),
        }
    }
}

/// `SiteWeights::probability` in `f32`, one lane per cell.
/// `emission_table_matches_the_trait` asserts the two agree.
///
/// Selects, not branches: whether a column converts is a property of the
/// haplotype, so a branch on it mispredicts once per cytosine -- and the vector
/// kernel has no branch to take at all. `matched * weight + offset` covers all
/// five rows of the emission because `weight = 1, offset = 0` reproduces a
/// plain match and `weight = 0` a plain mismatch, exactly, `matched` and
/// `offset` being non-negative.
#[inline]
fn prior<L: Lane>(column: ColumnLanes<L>, row: RowLanes<L>) -> L {
    let unknown = row.base.equals(L::splat(CODE_N)).either(column.base.equals(L::splat(CODE_N)));
    let plain = unknown.either(row.base.equals(column.plain));
    let weight = plain.select(
        L::splat(1.0),
        row.base.equals(column.converted).select(
            column.rate,
            row.base.equals(column.base).select(column.unconverted_rate, L::splat(0.0)),
        ),
    );
    row.matched * weight + plain.select(L::splat(0.0), row.mismatched)
}

/// The widest `Lane`, so the per-row tracks can be padded once.
const LANE_MAX: usize = 8;

/// One vector's worth of a track, whatever the `Lane`.
type Window = [f32; LANE_MAX];

/// Every buffer the kernel reads carries [`LANE_MAX`] slots of slack past the
/// largest index it can reach, so the zero fallback here is unreachable.
///
/// It is a fixed-size window rather than a slice because the fallback has to be
/// a *pointer*: hand `Lane::load` an `Option<&[f32]>` and the loaded vector
/// becomes the merge of two control-flow paths, which the register allocator
/// resolves by spilling a `q` pair to the stack -- once for each of the twenty
/// loads an eight-wide chunk takes.
#[inline]
fn window(track: &[f32], index: usize) -> &Window {
    static ZEROS: Window = [0.0; LANE_MAX];
    track.get(index..).and_then(|rest| rest.first_chunk::<LANE_MAX>()).unwrap_or(&ZEROS)
}

/// [`window`] for a store. There is no fallback: a store past the end of a DP
/// buffer would be a lost cell, not a harmless zero, so it is skipped instead.
#[inline]
fn window_mut(track: &mut [f32], index: usize) -> Option<&mut Window> {
    track.get_mut(index..).and_then(|rest| rest.first_chunk_mut::<LANE_MAX>())
}

const PAD: usize = 2;

/// The anti-diagonals a cell can reach back to, and so the length of the ring.
const RING: usize = 3;

/// `M`, `I` and `D`, and where each one's [`RING`] slots start in the kernel's
/// single allocation.
const MATRICES: usize = 3;
const MATCH: usize = 0;
const INSERTION: usize = RING;
const DELETION: usize = 2 * RING;

fn ceil_div2(value: i64) -> i64 {
    (value + 1).div_euclid(2)
}

/// Cells on one anti-diagonal `k = i + j` depend only on diagonals `k - 1` and
/// `k - 2`, so a whole diagonal is computed in parallel. Within a diagonal the
/// row index `i` runs over a contiguous range, which is what makes the read
/// qualities and the three neighbour buffers contiguous loads rather than
/// gathers; only the haplotype bases run backwards, and those reach the kernel
/// through the scalar emission call anyway.
#[allow(
    clippy::too_many_lines,
    reason = "one traversal; splitting it would hide the index algebra it exists to get right"
)]
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap,
    reason = "the f32 narrowing is the point of this kernel, and every index cast is clamped above"
)]
#[allow(
    clippy::indexing_slicing,
    reason = "every literal index here is a ring-buffer slot in 0..9 or a lane in 0..8; the DP's own buffers are reached through `get`"
)]
fn banded_kernel<L: Lane, E: Emission>(
    haplotype: &Haplotype,
    read: &Read,
    emission: &E,
    band: Band,
) -> Log10Likelihood {
    let (h, r) = (haplotype.len(), read.len());
    if h == 0 || r == 0 {
        return Log10Likelihood::IMPOSSIBLE;
    }
    let Some(plan) = Plan::build(haplotype, read, emission) else {
        return Log10Likelihood::IMPOSSIBLE;
    };
    let (hi64, ri64) = (h as i64, r as i64);
    let (w, o) = (band.half_width, band.offset);

    // A diagonal holds `hi - lo + 1` cells with `lo, hi` in `0..=r` and
    // `hi - lo <= half_width`, so this bound is tight at both ends.
    let max_len = (w.min(ri64) as usize) + 2;
    let capacity = 2 * PAD + max_len + 4 * L::LANES + 8;
    // Three matrices over three anti-diagonals, in one allocation rather than
    // nine. They are a few hundred bytes each and C5 scores one read against up
    // to eight haplotypes, so nine `malloc`/`free` pairs per alignment were
    // most of what a short one spent outside the loop.
    let mut cells = vec![0.0f32; MATRICES * RING * capacity];
    let mut ring = cells.chunks_exact_mut(capacity);
    let slots: [&mut [f32]; MATRICES * RING] =
        core::array::from_fn(|_| ring.next().unwrap_or_default());
    let mut lo_history = [0usize; RING];

    let init = 1.0f32 / h as f32;
    let mut init_scaled = init;
    let mut exponent = 0i32;
    let mut previous_shift = 0i32;
    let mut accumulator = 0.0f64;
    let mut reference_exponent: Option<i32> = None;

    for k in 0..=(r + h) {
        let ki = k as i64;
        // The band's own lower bound, before it is clamped into `0..=r` to be a
        // usable row index. The clamp can push `lo` *below* the band, and the
        // read's free start must not be granted on a cell the band excludes.
        let lo_in_band = ceil_div2(ki - o - w).max(ki - hi64).max(0);
        let lo = lo_in_band.min(ri64.min(ki)) as usize;
        let hi_signed = (ki - o + w).div_euclid(2).min(ri64).min(ki);
        let len = if hi_signed < lo as i64 { 0 } else { (hi_signed - lo as i64 + 1) as usize }
            .min(max_len);

        let cur = k % 3;
        let prev1 = (k + 2) % 3;
        let prev2 = (k + 1) % 3;
        let back1 = lo - lo_history[prev1].min(lo);
        let back2 = lo - lo_history[prev2].min(lo);

        // The cells this diagonal owns. A diagonal that reaches `j == 0` has a
        // topmost cell with no haplotype column behind it, and that cell is not
        // one of them: masking it out here is what used to be a `set` to zero
        // after the fact.
        let live = if len > 0 && lo + len - 1 == k { len - 1 } else { len };

        // Cell `(i, j)` on this diagonal reads haplotype column `j - 1`, which
        // the plan holds reversed at `h - j = h - k + i`. `i` ascends with the
        // lane, so this index does too: the haplotype is a contiguous load and
        // not a gather.
        let column_origin = (hi64 + lo as i64 - ki).max(0) as usize;

        let lift = L::splat(exp2_f32(previous_shift));
        let zero = L::splat(0.0);
        let offsets = L::offsets();
        let mut running = zero;
        let mut chunk = 0usize;
        while chunk < len {
            let diag = PAD + chunk + back2 - 1;
            let up = PAD + chunk + back1 - 1;
            let left = PAD + chunk + back1;
            let m_diag = L::load(window(slots[MATCH + prev2], diag));
            let i_diag = L::load(window(slots[INSERTION + prev2], diag));
            let d_diag = L::load(window(slots[DELETION + prev2], diag));
            let m_up = L::load(window(slots[MATCH + prev1], up));
            let i_up = L::load(window(slots[INSERTION + prev1], up));
            let m_left = L::load(window(slots[MATCH + prev1], left));
            let d_left = L::load(window(slots[DELETION + prev1], left));

            // `i = lo + chunk + lane` runs upwards over the lanes, so every
            // per-row track is one contiguous load.
            let row = lo + chunk;
            let prior_v = prior::<L>(
                ColumnLanes::load(&plan, column_origin + chunk),
                RowLanes::load(&plan, row),
            );
            let m2m_v = L::load(window(&plan.match_to_match, row));
            let m2i_v = L::load(window(&plan.match_to_insertion, row));
            let m2d_v = L::load(window(&plan.match_to_deletion, row));
            let i2m_v = L::load(window(&plan.indel_to_match, row));
            let gap_v = L::load(window(&plan.gap_continuation, row));

            let new_m = prior_v * (m_diag * m2m_v + i_diag * i2m_v + d_diag * i2m_v) * lift;
            let new_i = m_up * m2i_v + i_up * gap_v;
            let new_d = m_left * m2d_v + d_left * gap_v;

            #[allow(
                clippy::cast_precision_loss,
                reason = "a diagonal holds at most `half_width + 1` cells"
            )]
            let keep = offsets.below(L::splat(live.saturating_sub(chunk) as f32));
            let (new_m, new_i, new_d) =
                (keep.select(new_m, zero), keep.select(new_i, zero), keep.select(new_d, zero));
            running = running.vmax(new_m).vmax(new_i).vmax(new_d);

            let slot = PAD + chunk;
            for (matrix, value) in [(MATCH, new_m), (INSERTION, new_i), (DELETION, new_d)] {
                if let Some(destination) = window_mut(slots[matrix + cur], slot) {
                    value.store(destination);
                }
            }
            chunk += L::LANES;
        }

        // Above the live range the buffer has to read as zero for the two
        // diagonals that will load across it: they reach one vector past their
        // own live range at an offset of up to `PAD`, and `len` grows by at
        // most one per diagonal. A constant three vectors covers that, and
        // being constant is the point -- the variable-length clear this
        // replaces was a `memset` call per matrix per diagonal.
        for step in 0..3 {
            let slot = PAD + chunk + step * L::LANES;
            for matrix in [MATCH, INSERTION, DELETION] {
                if let Some(destination) = window_mut(slots[matrix + cur], slot) {
                    zero.store(destination);
                }
            }
        }

        // The read's free start. It replaces a cell that computed to exactly
        // zero in all three matrices -- read row 0 carries no base, so its
        // prior is zero, and its five transitions are zero too -- which is also
        // why the running maximum above could be taken before this write.
        if len > 0 && lo_in_band == 0 {
            set(slots[DELETION + cur], PAD, init_scaled);
            running = running.vmax(L::splat(init_scaled));
        }

        let shift = normalising_shift_f32(running.horizontal_max());
        if shift != 0 {
            let factor = L::splat(exp2_f32(shift));
            let mut slot = PAD;
            while slot < PAD + chunk {
                for matrix in [MATCH, INSERTION, DELETION] {
                    if let Some(destination) = window_mut(slots[matrix + cur], slot) {
                        (L::load(destination) * factor).store(destination);
                    }
                }
                slot += L::LANES;
            }
            exponent += shift;
            if lo_in_band == 0 {
                init_scaled *= exp2_f32(shift);
            }
        }
        previous_shift = shift;

        if len > 0 && lo + len - 1 == r {
            let slot = PAD + (r - lo);
            let value =
                f64::from(get(slots[MATCH + cur], slot) + get(slots[INSERTION + cur], slot));
            let base = *reference_exponent.get_or_insert(exponent);
            accumulator += value * exp2_f64(base - exponent);
        }

        lo_history[cur] = lo;
    }

    match reference_exponent {
        Some(base) if accumulator > 0.0 => {
            Log10Likelihood::new(accumulator.log10() - f64::from(base) * core::f64::consts::LOG10_2)
        }
        _ => Log10Likelihood::IMPOSSIBLE,
    }
}

#[inline]
fn set(buffer: &mut [f32], index: usize, value: f32) {
    if let Some(slot) = buffer.get_mut(index) {
        *slot = value;
    }
}

#[inline]
fn get(buffer: &[f32], index: usize) -> f32 {
    buffer.get(index).copied().unwrap_or(0.0)
}

#[cfg(test)]
mod tests {
    use super::{ColumnLanes, Plan, RowLanes, prior};
    use crate::{
        BaseQuality, Betas, ConversionModel, Emission, Haplotype, Probability, Read,
        StandardEmission, Strand, TapsEmission,
    };

    /// The kernel's hoisted tables are a re-encoding of the `Emission` trait,
    /// not a second model. This asserts it cell by cell over the whole matrix:
    /// `prior` in the scalar lane is `match_probability` narrowed to `f32`, to
    /// within the rounding that narrowing the factors first can cost.
    ///
    /// Four ulp rather than one, and only for the conversion rows: the trait
    /// computes `(1 - eps) * rate + eps / 3` in `f64` and rounds once, while the
    /// kernel rounds `1 - eps`, `rate` and `eps / 3` first and then multiplies
    /// and adds them unfused -- unfused because a fused vector multiply-add is
    /// not available on every target and bit-parity has to be. Every other row
    /// is exact.
    #[test]
    fn emission_table_matches_the_trait() {
        let haplotype = Haplotype::from_ascii(b"ACGTNCGGCGATCGATTACGNNCGATCGGGCCATCGAT");
        let betas: Vec<Probability> = (0..haplotype.len())
            .map(|index| {
                Probability::new(f64::from(u32::try_from(index % 7).unwrap_or(0)) / 6.0)
                    .unwrap_or(Probability::ZERO)
            })
            .collect();

        for strand in [Strand::OT, Strand::OB] {
            let bases = haplotype.bases().to_vec();
            let quals: Vec<BaseQuality> = (0..bases.len())
                .map(|index| BaseQuality::from_byte(u8::try_from(index % 60).unwrap_or(30)))
                .collect();
            let read = Read::uniform(
                bases,
                &quals,
                BaseQuality::from_byte(45),
                BaseQuality::from_byte(45),
                BaseQuality::from_byte(10),
                strand,
            )
            .expect("valid");

            let floor = Probability::new(0.01).expect("in [0, 1]");
            let standard = StandardEmission::default().with_artifact_floor(floor);
            let taps = TapsEmission::new(ConversionModel::taps_default(), Betas::PerSite(&betas))
                .with_artifact_floor(floor);

            for (name, exact) in [
                ("standard", check(&haplotype, &read, &standard)),
                ("taps", check(&haplotype, &read, &taps)),
            ] {
                assert!(exact, "{name} on {strand:?}: the table left the trait");
            }
        }
    }

    fn check<E: Emission>(haplotype: &Haplotype, read: &Read, emission: &E) -> bool {
        let plan = Plan::build(haplotype, read, emission).expect("the fixture is well formed");
        let Some(last) = haplotype.len().checked_sub(1) else {
            return false;
        };
        for column_index in 0..haplotype.len() {
            for row_index in 1..=read.len() {
                let (Some(site), Some(observation)) =
                    (haplotype.site(column_index), read.observation(row_index - 1))
                else {
                    return false;
                };
                #[allow(
                    clippy::cast_possible_truncation,
                    reason = "the narrowing is what is under test"
                )]
                let want = emission.match_probability(site, observation) as f32;
                let got = prior::<f32>(
                    ColumnLanes::load(&plan, last - column_index),
                    RowLanes::load(&plan, row_index),
                );
                let slack = 4.0 * f32::EPSILON * want.abs();
                if (got - want).abs() > slack {
                    return false;
                }
            }
        }
        true
    }
}
