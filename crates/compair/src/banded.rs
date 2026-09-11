use wide::{Select, f32x8};

use crate::{
    emission::Emission,
    error::Error,
    haplotype::Haplotype,
    read::Read,
    scaling::{exp2_f32, exp2_f64, normalising_shift_f32},
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
    pub(crate) half_width: i64,
    pub(crate) offset: i64,
}

impl Band {
    /// The width §6.4 of the design budgets for: 150 bp of read against ~48
    /// columns is ~7,000 cells.
    ///
    /// It is 46 and not 48 because of how the SIMD kernel fills its lanes. An
    /// anti-diagonal holds at most `width / 2 + 1` cells, and the kernel
    /// computes them eight at a time: at a half-width of 23 that is exactly
    /// three vectors, while 24 spills one cell into a fourth and costs ~8% for
    /// one more column on each side. Widths of the form `16n - 2` fill every
    /// vector; widths of the form `16n` waste one per diagonal.
    pub const DEFAULT_WIDTH: u32 = 46;

    /// The widest band either kernel will build, for two reasons that are the
    /// same reason.
    ///
    /// The banded kernels allocate nine buffers of `width / 2 + O(1)` floats,
    /// so an unbounded `width` is an unbounded allocation inside a safe
    /// function -- `u32::MAX` would ask for ~77 GB.
    ///
    /// And a wide band is not a better answer, it is a worse one. The `f32`
    /// kernels renormalise per **anti-diagonal**, and a diagonal holds one cell
    /// per read row it crosses -- prefix alignments of that many lengths, whose
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

// `anchored` builds a band without going through `new`, so the default has to
// satisfy `new`'s checks by construction.
const _: () = assert!(Band::DEFAULT_WIDTH >= 2 && Band::DEFAULT_WIDTH <= Band::MAX_WIDTH);

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
pub(crate) trait Lane:
    Copy + core::ops::Add<Output = Self> + core::ops::Mul<Output = Self>
{
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
    /// Bitwise and of two masks.
    fn both(self, other: Self) -> Self;
    /// `self` is a mask from [`Lane::equals`] or [`Lane::either`].
    fn select(self, if_true: Self, if_false: Self) -> Self;
    /// Every lane moved up by one: `first` enters at lane 0 and the last lane
    /// leaves. On one lane, `first`.
    fn shift_in(self, first: f32) -> Self;
    /// The lane [`Lane::shift_in`] would drop.
    fn last(self) -> f32;
    /// The sum of the lanes. Only ever called with at most one non-zero lane,
    /// so the order of the additions is not observable.
    fn horizontal_sum(self) -> f32;
}

/// What a lanewise comparison reports for "equal". As a number it is a `NaN`;
/// it is only ever consumed bitwise.
pub(crate) const MASK_SET: f32 = f32::from_bits(u32::MAX);

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
    fn both(self, other: Self) -> Self {
        Self::from_bits(self.to_bits() & other.to_bits())
    }
    #[inline]
    fn select(self, if_true: Self, if_false: Self) -> Self {
        if self.to_bits() == 0 { if_false } else { if_true }
    }
    #[inline]
    fn shift_in(self, first: f32) -> Self {
        first
    }
    #[inline]
    fn last(self) -> f32 {
        self
    }
    #[inline]
    fn horizontal_sum(self) -> f32 {
        self
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
        // A tree, not a fold: three dependent `fmaxnm` instead of seven, on
        // the one reduction that sits on every diagonal's critical path.
        let [a, b, c, d, e, f, g, h] = self.to_array();
        a.max(b).max(c.max(d)).max(e.max(f).max(g.max(h)))
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
    fn both(self, other: Self) -> Self {
        self & other
    }
    #[inline]
    fn select(self, if_true: Self, if_false: Self) -> Self {
        Select::select(self, if_true, if_false)
    }
    #[inline]
    fn shift_in(self, first: f32) -> Self {
        // Written as a permutation of the array so that it holds by
        // construction; the compiler recognises the pattern and emits lane
        // permutes, see `strips_hot_loop_has_no_stack_traffic` in the
        // profiling notes.
        let [a, b, c, d, e, f, g, _] = self.to_array();
        Self::from([first, a, b, c, d, e, f, g])
    }
    #[inline]
    fn last(self) -> f32 {
        let [_, _, _, _, _, _, _, h] = self.to_array();
        h
    }
    #[inline]
    fn horizontal_sum(self) -> f32 {
        self.reduce_add()
    }
}

/// The buffers a banded alignment works in, kept between calls.
///
/// One alignment needs fifteen small buffers: the hoisted per-row and
/// per-column tracks of the [`Plan`] and the ring of anti-diagonals. A caller
/// scores one read against several haplotypes and many reads in a row, so
/// allocating them per call was a measurable share of a short alignment.
/// Reusing a workspace makes an alignment allocation-free; the free functions
/// [`align_banded`] and [`align_banded_simd`] build a fresh one each time.
#[derive(Debug, Default)]
pub struct Workspace {
    pub(crate) plan: Plan,
    ring: Ring,
    pub(crate) rows: crate::strips::RowBuffer,
}

impl Workspace {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Scalar `f32` over the band, in this workspace.
    pub fn align_banded<E: Emission>(
        &mut self,
        haplotype: &Haplotype,
        read: &Read,
        emission: &E,
        band: Band,
    ) -> Log10Likelihood {
        self.align::<f32, E>(haplotype, read, emission, band)
    }

    /// Eight `f32` lanes of one anti-diagonal at a time, in this workspace.
    pub fn align_banded_simd<E: Emission>(
        &mut self,
        haplotype: &Haplotype,
        read: &Read,
        emission: &E,
        band: Band,
    ) -> Log10Likelihood {
        self.align::<f32x8, E>(haplotype, read, emission, band)
    }

    fn align<L: Lane, E: Emission>(
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
        // The emission is folded into the plan, so from here on the kernel is
        // generic over the lane type only: one hot loop per lane, not one per
        // emission model.
        banded_kernel::<L>(&self.plan, &mut self.ring, Shape { haplotype: h, read: r }, band)
    }
}

/// Scalar `f32` over the band.
pub fn align_banded<E: Emission>(
    haplotype: &Haplotype,
    read: &Read,
    emission: &E,
    band: Band,
) -> Log10Likelihood {
    Workspace::new().align_banded(haplotype, read, emission, band)
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
    Workspace::new().align_banded_simd(haplotype, read, emission, band)
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
/// runs along. The column tracks are **reversed** -- index `COLUMN_FRONT + h -
/// 1 - j` holds the column a cell at haplotype position `j` reads -- because
/// `j` runs down an anti-diagonal while `i` runs up, and reversing is what
/// makes the haplotype a second ascending contiguous load rather than a
/// gather. Index `COLUMN_FRONT + h` is the `j == 0` boundary column, which has
/// no site; the kernel masks that cell to zero, so the padding there only has
/// to be finite.
///
/// Every track carries [`TRACK_SLACK`] entries past its live range so that
/// the last window a diagonal loads is still in bounds, and the column tracks
/// carry [`COLUMN_FRONT`] entries before it: a lane that has run past the
/// end of the haplotype loads from there, masked, and the window it loads
/// still has to exist.
#[derive(Debug, Default)]
pub(crate) struct Plan {
    pub(crate) rows: RowTracks,
    pub(crate) columns: ColumnTracks,
}

/// One entry per read row, row 0 being the base-less start row.
#[derive(Debug, Default)]
pub(crate) struct RowTracks {
    pub(crate) base: Vec<f32>,
    pub(crate) matched: Vec<f32>,
    pub(crate) mismatched: Vec<f32>,
    pub(crate) match_to_match: Vec<f32>,
    pub(crate) match_to_insertion: Vec<f32>,
    pub(crate) match_to_deletion: Vec<f32>,
    pub(crate) indel_to_match: Vec<f32>,
    pub(crate) gap_continuation: Vec<f32>,
}

/// One entry per haplotype column, reversed; see [`Plan`].
#[derive(Debug, Default)]
pub(crate) struct ColumnTracks {
    pub(crate) base: Vec<f32>,
    pub(crate) converted: Vec<f32>,
    pub(crate) plain: Vec<f32>,
    pub(crate) rate: Vec<f32>,
    pub(crate) unconverted: Vec<f32>,
}

/// `len` copies of `fill`, reusing the allocation.
pub(crate) fn reset(track: &mut Vec<f32>, len: usize, fill: f32) {
    track.clear();
    track.resize(len, fill);
}

impl Plan {
    /// Rebuilds every track for this pair, reusing the allocations. `None`
    /// only if the read or haplotype cannot be addressed, which their
    /// constructors rule out.
    #[allow(
        clippy::cast_possible_truncation,
        reason = "the f32 narrowing is the point of this kernel"
    )]
    pub(crate) fn fill<E: Emission>(
        &mut self,
        haplotype: &Haplotype,
        read: &Read,
        emission: &E,
    ) -> Option<()> {
        let (h, r) = (haplotype.len(), read.len());
        let rows = r + 1 + TRACK_SLACK;
        let columns = COLUMN_FRONT + h + 1 + TRACK_SLACK;

        let RowTracks {
            base,
            matched,
            mismatched,
            match_to_match,
            match_to_insertion,
            match_to_deletion,
            indel_to_match,
            gap_continuation,
        } = &mut self.rows;
        for track in [
            &mut *base,
            matched,
            mismatched,
            match_to_match,
            match_to_insertion,
            match_to_deletion,
            indel_to_match,
            gap_continuation,
        ] {
            reset(track, rows, 0.0);
        }
        for index in 0..r {
            let observation = read.observation(index)?;
            let eps = emission.epsilon(observation);
            let t = read.transition(index)?;
            let row = index + 1;
            *base.get_mut(row)? = code(observation.base);
            *matched.get_mut(row)? = (1.0 - eps) as f32;
            *mismatched.get_mut(row)? = (eps / 3.0) as f32;
            *match_to_match.get_mut(row)? = t.match_to_match as f32;
            *match_to_insertion.get_mut(row)? = t.match_to_insertion as f32;
            *match_to_deletion.get_mut(row)? = t.match_to_deletion as f32;
            *indel_to_match.get_mut(row)? = t.indel_to_match as f32;
            *gap_continuation.get_mut(row)? = t.gap_continuation as f32;
        }

        let ColumnTracks { base, converted, plain, rate, unconverted } = &mut self.columns;
        reset(base, columns, CODE_NO_PLAIN_MATCH);
        reset(converted, columns, CODE_NO_CONVERSION);
        reset(plain, columns, CODE_NO_PLAIN_MATCH);
        reset(rate, columns, 0.0);
        reset(unconverted, columns, 0.0);
        let strand = read.strand();
        for index in 0..h {
            let weights = emission.site_weights(haplotype.site(index)?, strand);
            let site_base = code(weights.base);
            let reversed = COLUMN_FRONT + h - 1 - index;
            *base.get_mut(reversed)? = site_base;
            match weights.converted {
                Some(converted_base) => {
                    *converted.get_mut(reversed)? = code(converted_base);
                    *rate.get_mut(reversed)? = weights.rate as f32;
                    *unconverted.get_mut(reversed)? = (1.0 - weights.rate) as f32;
                }
                None => *plain.get_mut(reversed)? = site_base,
            }
        }
        Some(())
    }
}

/// One haplotype column per lane, in the reversed order the diagonal reads.
#[derive(Clone, Copy)]
pub(crate) struct ColumnLanes<L> {
    pub(crate) base: L,
    pub(crate) converted: L,
    pub(crate) plain: L,
    pub(crate) rate: L,
    pub(crate) unconverted_rate: L,
}

/// One read row per lane: what the emission needs.
#[derive(Clone, Copy)]
pub(crate) struct RowLanes<L> {
    base: L,
    matched: L,
    mismatched: L,
    /// `base` is `N`: a mask, computed here so that a kernel whose rows are
    /// fixed for many cells computes it once.
    unknown: L,
}

impl<L: Lane> RowLanes<L> {
    #[inline(always)]
    pub(crate) fn new(base: L, matched: L, mismatched: L) -> Self {
        Self { base, matched, mismatched, unknown: base.equals(L::splat(CODE_N)) }
    }
}

/// One read row per lane: what the transitions need.
#[derive(Clone, Copy)]
pub(crate) struct TransitionLanes<L> {
    pub(crate) match_to_match: L,
    pub(crate) match_to_insertion: L,
    pub(crate) match_to_deletion: L,
    pub(crate) indel_to_match: L,
    pub(crate) gap_continuation: L,
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
pub(crate) fn prior<L: Lane>(column: ColumnLanes<L>, row: RowLanes<L>) -> L {
    let unknown = row.unknown.either(column.base.equals(L::splat(CODE_N)));
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
pub(crate) const LANE_MAX: usize = 8;

/// One vector's worth of a track, whatever the `Lane`.
pub(crate) type Window = [f32; LANE_MAX];

/// Past the last live entry of a plan track: a full chunk of the widest lane
/// beyond the last cell, and the window that chunk's last lane reads.
pub(crate) const TRACK_SLACK: usize = 2 * LANE_MAX;

/// Before the first live entry of a column track: one window of the widest
/// lane, for a lane that reads up to `LANE_MAX - 1` columns past the end of
/// the haplotype.
pub(crate) const COLUMN_FRONT: usize = LANE_MAX;

/// Slack below slot zero, so that `diag` and `up` at chunk zero of a diagonal
/// whose predecessor started one row later still address the buffer.
const PAD: usize = 2;

/// The anti-diagonals a cell can reach back to, and so the length of the ring.
const RING: usize = 3;

/// `M`, `I` and `D`.
const MATRICES: usize = 3;

/// Three matrices over three anti-diagonals, in one allocation rather than
/// nine, laid out phase-major: the three matrices of one diagonal are one
/// contiguous region, so a diagonal borrows its own region mutably and the two
/// before it immutably.
#[derive(Debug, Default)]
struct Ring {
    cells: Vec<f32>,
    capacity: usize,
}

/// The three matrices of one finished anti-diagonal.
#[derive(Clone, Copy)]
struct Phase<'a> {
    m: &'a [f32],
    i: &'a [f32],
    d: &'a [f32],
}

/// The three matrices of the anti-diagonal being computed.
struct PhaseMut<'a> {
    m: &'a mut [f32],
    i: &'a mut [f32],
    d: &'a mut [f32],
}

impl<'a> Phase<'a> {
    fn split(region: &'a [f32], capacity: usize) -> Option<Self> {
        let (m, rest) = region.split_at_checked(capacity)?;
        let (i, d) = rest.split_at_checked(capacity)?;
        Some(Self { m, i, d })
    }
}

impl<'a> PhaseMut<'a> {
    fn split(region: &'a mut [f32], capacity: usize) -> Option<Self> {
        let (m, rest) = region.split_at_mut_checked(capacity)?;
        let (i, d) = rest.split_at_mut_checked(capacity)?;
        Some(Self { m, i, d })
    }
}

impl Ring {
    /// Zeroes the ring at a new per-diagonal capacity, reusing the allocation.
    fn reset(&mut self, capacity: usize) {
        self.capacity = capacity;
        reset(&mut self.cells, MATRICES * RING * capacity, 0.0);
    }

    /// The diagonal at ring phase `cur`, to write, and the two before it.
    ///
    /// `reset` sizes `cells` to exactly `MATRICES * RING * capacity`, so the
    /// `None` is unreachable; it is typed rather than asserted so that the
    /// kernel returns `IMPOSSIBLE` instead of panicking should that ever stop
    /// being true.
    #[inline(always)]
    fn phases(&mut self, cur: usize) -> Option<(PhaseMut<'_>, Phase<'_>, Phase<'_>)> {
        let region = MATRICES * self.capacity;
        let (first, rest) = self.cells.split_at_mut_checked(region)?;
        let (second, third) = rest.split_at_mut_checked(region)?;
        let (current, previous, before) = match cur {
            0 => (first, third, second),
            1 => (second, first, third),
            _ => (third, second, first),
        };
        Some((
            PhaseMut::split(current, self.capacity)?,
            Phase::split(previous, self.capacity)?,
            Phase::split(before, self.capacity)?,
        ))
    }
}

/// A track from one origin, exactly `span` entries long, so that every window
/// a diagonal will load from it has been proven in bounds once.
///
/// This is why the kernel has an `unsafe` block at all. With a checked
/// `get(..).first_chunk()` per window, the compiler emitted, for each of the
/// twenty loads a chunk takes, a compare, a conditional compare and a select
/// on a pointer and a length it had spilled to the stack -- two thirds of the
/// inner loop's instructions were bounds checks on loads that cannot fail.
/// Here the check moves to [`View::new`], once per track per diagonal, and
/// [`Sources::step`] checks the one offset it uses against the shared span.
#[derive(Clone, Copy)]
pub(crate) struct View<'a> {
    track: &'a [f32],
}

impl<'a> View<'a> {
    /// `None` unless `track[origin..origin + span]` exists.
    #[inline]
    pub(crate) fn new(track: &'a [f32], origin: usize, span: usize) -> Option<Self> {
        Some(Self { track: track.get(origin..origin.checked_add(span)?)? })
    }

    /// The window at `offset`.
    ///
    /// # Safety
    ///
    /// `offset + LANE_MAX <= span`, the length `new` was given.
    #[inline(always)]
    pub(crate) unsafe fn window(self, offset: usize) -> &'a Window {
        debug_assert!(offset + LANE_MAX <= self.track.len());
        // SAFETY: `new` proved `track` holds `span` initialised `f32`s and the
        // caller keeps `offset + LANE_MAX` within it; `[f32; LANE_MAX]` has
        // the alignment of `f32` and no invalid bit patterns.
        unsafe { &*self.track.as_ptr().add(offset).cast::<Window>() }
    }
}

/// [`View`] for the diagonal being written.
struct ViewMut<'a> {
    track: &'a mut [f32],
}

impl<'a> ViewMut<'a> {
    #[inline]
    fn new(track: &'a mut [f32], origin: usize, span: usize) -> Option<Self> {
        Some(Self { track: track.get_mut(origin..origin.checked_add(span)?)? })
    }

    /// # Safety
    ///
    /// `offset + LANE_MAX <= span`, the length `new` was given.
    #[inline(always)]
    unsafe fn window(&mut self, offset: usize) -> &mut Window {
        debug_assert!(offset + LANE_MAX <= self.track.len());
        // SAFETY: as for `View::window`, through a unique borrow.
        unsafe { &mut *self.track.as_mut_ptr().add(offset).cast::<Window>() }
    }
}

#[inline]
fn ceil_div2(value: i64) -> i64 {
    (value + 1).div_euclid(2)
}

/// The two dimensions of the matrix.
#[derive(Clone, Copy)]
pub(crate) struct Shape {
    pub(crate) haplotype: usize,
    pub(crate) read: usize,
}

/// Where one anti-diagonal reads its neighbours from.
#[derive(Clone, Copy)]
struct Diagonal {
    /// The lowest read row on this diagonal.
    lo: usize,
    /// How many slots later than this diagonal the previous one started.
    back1: usize,
    /// The same for the one before it.
    back2: usize,
    /// The plan's reversed column index of the cell at `lo`.
    column_origin: usize,
    /// The chunk offsets this diagonal visits are `0..extent` in steps of the
    /// lane count; `extent` is the live length rounded up to a whole chunk.
    extent: usize,
}

/// Everything one anti-diagonal loads, each track viewed from the origin its
/// chunk offset zero reads, all with one span.
struct Sources<'a> {
    span: usize,
    m_diag: View<'a>,
    i_diag: View<'a>,
    d_diag: View<'a>,
    m_up: View<'a>,
    i_up: View<'a>,
    m_left: View<'a>,
    d_left: View<'a>,
    row_base: View<'a>,
    row_matched: View<'a>,
    row_mismatched: View<'a>,
    match_to_match: View<'a>,
    match_to_insertion: View<'a>,
    match_to_deletion: View<'a>,
    indel_to_match: View<'a>,
    gap_continuation: View<'a>,
    column_base: View<'a>,
    column_converted: View<'a>,
    column_plain: View<'a>,
    column_rate: View<'a>,
    column_unconverted: View<'a>,
}

impl<'a> Sources<'a> {
    /// `None` if any track is too short for this diagonal, which the plan's
    /// and the ring's slack rule out.
    ///
    /// `inline(always)`, like `Ring::phases` and `Sink::new`: left to itself
    /// the compiler kept this out of line, wrote the twenty views to the
    /// stack and read them back inside the chunk loop, which was most of the
    /// per-diagonal overhead.
    #[inline(always)]
    fn new(plan: &'a Plan, prev1: Phase<'a>, prev2: Phase<'a>, diagonal: Diagonal) -> Option<Self> {
        let Diagonal { lo, back1, back2, column_origin, extent } = diagonal;
        let span = extent + LANE_MAX;
        let diag = PAD + back2 - 1;
        let up = PAD + back1 - 1;
        let left = PAD + back1;
        let rows = &plan.rows;
        let columns = &plan.columns;
        Some(Self {
            span,
            m_diag: View::new(prev2.m, diag, span)?,
            i_diag: View::new(prev2.i, diag, span)?,
            d_diag: View::new(prev2.d, diag, span)?,
            m_up: View::new(prev1.m, up, span)?,
            i_up: View::new(prev1.i, up, span)?,
            m_left: View::new(prev1.m, left, span)?,
            d_left: View::new(prev1.d, left, span)?,
            row_base: View::new(&rows.base, lo, span)?,
            row_matched: View::new(&rows.matched, lo, span)?,
            row_mismatched: View::new(&rows.mismatched, lo, span)?,
            match_to_match: View::new(&rows.match_to_match, lo, span)?,
            match_to_insertion: View::new(&rows.match_to_insertion, lo, span)?,
            match_to_deletion: View::new(&rows.match_to_deletion, lo, span)?,
            indel_to_match: View::new(&rows.indel_to_match, lo, span)?,
            gap_continuation: View::new(&rows.gap_continuation, lo, span)?,
            column_base: View::new(&columns.base, column_origin, span)?,
            column_converted: View::new(&columns.converted, column_origin, span)?,
            column_plain: View::new(&columns.plain, column_origin, span)?,
            column_rate: View::new(&columns.rate, column_origin, span)?,
            column_unconverted: View::new(&columns.unconverted, column_origin, span)?,
        })
    }

    /// One vector of cells at chunk offset `chunk`: the recurrence, with
    /// nothing masked yet. `None` past the span, which the kernel's chunk loop
    /// never reaches.
    ///
    /// `lift1` brings the previous diagonal onto this diagonal's scale and
    /// `lift_before` brings the one before it onto the previous diagonal's;
    /// see `banded_kernel` for why the scaling is applied on the way in
    /// rather than in place. They are two factors and not their product on
    /// purpose: each is at most `2^126`, but a diagonal whose maximum
    /// collapses -- the free-start cell leaving the band with only a long
    /// deletion chain behind -- can ask for a shift of over a hundred, and
    /// the product of two such shifts is not an `f32`. Applied one at a
    /// time, every intermediate is a value some diagonal already held.
    #[inline(always)]
    fn step<L: Lane>(&self, chunk: usize, lift1: L, lift_before: L) -> Option<(L, L, L)> {
        if chunk.checked_add(LANE_MAX)? > self.span {
            return None;
        }
        // SAFETY: every view was built with `self.span` entries and
        // `chunk + LANE_MAX <= self.span` was just checked.
        let load = |view: View<'a>| unsafe { L::load(view.window(chunk)) };

        let m_diag = load(self.m_diag);
        let i_diag = load(self.i_diag);
        let d_diag = load(self.d_diag);
        let m_up = load(self.m_up);
        let i_up = load(self.i_up);
        let m_left = load(self.m_left);
        let d_left = load(self.d_left);

        // `i = lo + chunk + lane` runs upwards over the lanes, so every per-row
        // track is one contiguous load, and the reversed column tracks are too.
        let prior_v = prior::<L>(
            ColumnLanes {
                base: load(self.column_base),
                converted: load(self.column_converted),
                plain: load(self.column_plain),
                rate: load(self.column_rate),
                unconverted_rate: load(self.column_unconverted),
            },
            RowLanes::new(load(self.row_base), load(self.row_matched), load(self.row_mismatched)),
        );
        let t = TransitionLanes {
            match_to_match: load(self.match_to_match),
            match_to_insertion: load(self.match_to_insertion),
            match_to_deletion: load(self.match_to_deletion),
            indel_to_match: load(self.indel_to_match),
            gap_continuation: load(self.gap_continuation),
        };

        let new_m = prior_v
            * (m_diag * t.match_to_match + i_diag * t.indel_to_match + d_diag * t.indel_to_match)
            * lift_before
            * lift1;
        let new_i = (m_up * t.match_to_insertion + i_up * t.gap_continuation) * lift1;
        let new_d = (m_left * t.match_to_deletion + d_left * t.gap_continuation) * lift1;
        Some((new_m, new_i, new_d))
    }
}

/// The three matrices of the diagonal being written, from slot `PAD`, with
/// room for every chunk plus the three vectors of zeros cleared past them.
struct Sink<'a> {
    span: usize,
    m: ViewMut<'a>,
    i: ViewMut<'a>,
    d: ViewMut<'a>,
}

impl<'a> Sink<'a> {
    #[inline(always)]
    fn new<L: Lane>(cur: PhaseMut<'a>, diagonal: Diagonal) -> Option<Self> {
        let span = diagonal.extent + 3 * L::LANES + LANE_MAX;
        Some(Self {
            span,
            m: ViewMut::new(cur.m, PAD, span)?,
            i: ViewMut::new(cur.i, PAD, span)?,
            d: ViewMut::new(cur.d, PAD, span)?,
        })
    }

    /// Stores one vector of each matrix at `offset`, subnormals flushed to
    /// zero; `false` past the span, which the kernel never reaches.
    ///
    /// The flush is what makes "the band under-estimates, never
    /// over-estimates" a property of the kernel rather than an observation:
    /// a cell more than `2^-126` below its diagonal's maximum would be stored
    /// with fewer mantissa bits than a normal `f32` and round either way, and
    /// its share of the total is far below anything a caller could act on.
    /// Zero only ever removes mass. It also makes the kernel's numbers the
    /// same on a target that flushes denormals in hardware and one that does
    /// not.
    #[inline(always)]
    fn store<L: Lane>(&mut self, offset: usize, cells: (L, L, L)) -> bool {
        let Some(end) = offset.checked_add(LANE_MAX) else { return false };
        if end > self.span {
            return false;
        }
        let (m, i, d) = cells;
        let (tiny, zero) = (L::splat(f32::MIN_POSITIVE), L::splat(0.0));
        let (m, i, d) = (
            m.below(tiny).select(zero, m),
            i.below(tiny).select(zero, i),
            d.below(tiny).select(zero, d),
        );
        // SAFETY: every view was built with `self.span` entries and
        // `offset + LANE_MAX <= self.span` was just checked.
        unsafe {
            m.store(self.m.window(offset));
            i.store(self.i.window(offset));
            d.store(self.d.window(offset));
        }
        true
    }
}

/// Cells on one anti-diagonal `k = i + j` depend only on diagonals `k - 1` and
/// `k - 2`, so a whole diagonal is computed in parallel. Within a diagonal the
/// row index `i` runs over a contiguous range, which is what makes the read
/// tracks and the three neighbour buffers contiguous loads rather than
/// gathers; the haplotype runs the other way, and the plan holds it reversed
/// so that it is a contiguous load too.
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
fn banded_kernel<L: Lane>(
    plan: &Plan,
    ring: &mut Ring,
    shape: Shape,
    band: Band,
) -> Log10Likelihood {
    let Shape { haplotype: h, read: r } = shape;
    let (hi64, ri64) = (h as i64, r as i64);
    let (w, o) = (band.half_width, band.offset);

    // A diagonal holds `hi - lo + 1` cells with `lo, hi` in `0..=r` and
    // `hi - lo <= half_width`, so this bound is tight at both ends.
    let max_len = (w.min(ri64) as usize) + 2;
    // Room for the widest diagonal rounded up to whole chunks, the three
    // vectors of zeros cleared past it, and the window the last lane reads.
    let capacity = 2 * PAD + max_len + 4 * L::LANES + LANE_MAX;
    ring.reset(capacity);
    let mut lo_history = [0usize; RING];

    // Renormalisation is lazy. Each diagonal is computed and stored in one
    // scale, `2^exponent`, chosen before it is computed as the previous
    // diagonal's scale plus the shift that would have put *that* diagonal's
    // maximum into `[1, 2)`. The two diagonals a cell reads are then brought
    // onto the current scale by a power of two on the way in -- exact in
    // binary floating point -- rather than rewritten in place, so nothing
    // between one diagonal and the next waits on a store: no rescale pass,
    // and no branch on whether one is needed.
    let init = 1.0f32 / h as f32;
    let mut init_scaled = init;
    let mut exponent = 0i32;
    let (mut shift1, mut shift2) = (0i32, 0i32);
    let mut accumulator = 0.0f64;
    let mut reference_exponent: Option<i32> = None;

    let zero = L::splat(0.0);
    let offsets = L::offsets();

    for k in 0..=(r + h) {
        let ki = k as i64;
        // The band's own lower bound, before it is clamped into `0..=r` to be a
        // usable row index. The clamp can push `lo` *below* the band, and the
        // read's free start must not be granted on a cell the band excludes.
        let lo_in_band = ceil_div2(ki - o - w).max(ki - hi64).max(0);
        let lo = lo_in_band.min(ri64.min(ki)) as usize;
        let hi_signed = (ki - o + w).div_euclid(2).min(ri64).min(ki);
        let len = if hi_signed < lo as i64 { 0 } else { (hi_signed - lo as i64 + 1) as usize };
        debug_assert!(len <= max_len, "diagonal {k} holds {len} cells, over the bound {max_len}");
        let len = len.min(max_len);

        let phase = k % RING;
        let Some((cur, prev1, prev2)) = ring.phases(phase) else {
            return Log10Likelihood::IMPOSSIBLE;
        };
        // `lo` grows by at most one per diagonal, so these are 0..=1 and
        // 0..=2, which is what `PAD` is sized for.
        let back1 = lo - lo_history[(k + RING - 1) % RING].min(lo);
        let back2 = lo - lo_history[(k + RING - 2) % RING].min(lo);
        debug_assert!(back1 <= 1 && back2 <= PAD);

        // The cells this diagonal owns. A diagonal that reaches `j == 0` has a
        // topmost cell with no haplotype column behind it, and that cell is not
        // one of them: it is masked to zero rather than skipped, so that the
        // slot it occupies reads as zero for the diagonals that load across it.
        let live = if len > 0 && lo + len - 1 == k { len - 1 } else { len };

        // Cell `(i, j)` on this diagonal reads haplotype column `j - 1`, which
        // the plan holds reversed at `COLUMN_FRONT + h - j`, with `h - j = h -
        // k + i`. `i` ascends with the lane, so this index does too.
        let column_origin = COLUMN_FRONT + (hi64 + lo as i64 - ki).max(0) as usize;

        let extent = len.div_ceil(L::LANES) * L::LANES;
        let diagonal = Diagonal { lo, back1, back2, column_origin, extent };
        let Some(sources) = Sources::new(plan, prev1, prev2, diagonal) else {
            return Log10Likelihood::IMPOSSIBLE;
        };
        let Some(mut sink) = Sink::new::<L>(cur, diagonal) else {
            return Log10Likelihood::IMPOSSIBLE;
        };
        // This diagonal's scale, and the lifts onto it from the two before.
        exponent += shift1;
        init_scaled *= exp2_f32(shift1);
        let lift1 = L::splat(exp2_f32(shift1));
        let lift_before = L::splat(exp2_f32(shift2));
        // One running maximum per matrix: three short dependency chains
        // rather than one three times as long.
        let (mut running_m, mut running_i, mut running_d) = (zero, zero, zero);

        // Full vectors first, with nothing to mask.
        let mut chunk = 0usize;
        while chunk + L::LANES <= live {
            let Some(cells) = sources.step::<L>(chunk, lift1, lift_before) else {
                return Log10Likelihood::IMPOSSIBLE;
            };
            running_m = running_m.vmax(cells.0);
            running_i = running_i.vmax(cells.1);
            running_d = running_d.vmax(cells.2);
            if !sink.store(chunk, cells) {
                return Log10Likelihood::IMPOSSIBLE;
            }
            chunk += L::LANES;
        }
        // Then the one partial vector, its lanes past `live` zeroed. That
        // includes the `j == 0` cell, and it is what the scalar lane goes
        // through for that cell alone.
        if chunk < len {
            let Some((m, i, d)) = sources.step::<L>(chunk, lift1, lift_before) else {
                return Log10Likelihood::IMPOSSIBLE;
            };
            #[allow(
                clippy::cast_precision_loss,
                reason = "a diagonal holds at most `half_width + 1` cells"
            )]
            let keep = offsets.below(L::splat((live - chunk) as f32));
            let cells = (keep.select(m, zero), keep.select(i, zero), keep.select(d, zero));
            running_m = running_m.vmax(cells.0);
            running_i = running_i.vmax(cells.1);
            running_d = running_d.vmax(cells.2);
            if !sink.store(chunk, cells) {
                return Log10Likelihood::IMPOSSIBLE;
            }
            chunk += L::LANES;
        }
        debug_assert_eq!(chunk, extent);

        // Above the live range the buffer has to read as zero for the two
        // diagonals that will load across it: they reach one vector past their
        // own live range at an offset of up to `PAD`, and `len` grows by at
        // most one per diagonal. A constant three vectors covers that, and
        // being constant is the point -- the variable-length clear this
        // replaces was a `memset` call per matrix per diagonal.
        for vector in 0..3 {
            if !sink.store(extent + vector * L::LANES, (zero, zero, zero)) {
                return Log10Likelihood::IMPOSSIBLE;
            }
        }

        // The read's free start. It replaces a cell that computed to exactly
        // zero in all three matrices -- read row 0 carries no base, so its
        // prior is zero, and its five transitions are zero too -- which is also
        // why the running maximum above could be taken before this write.
        if len > 0 && lo_in_band == 0 {
            if let Some(start) = sink.d.track.first_mut() {
                *start = init_scaled;
            }
            running_d = running_d.vmax(L::splat(init_scaled));
        }

        // The shift this diagonal's maximum asks for is applied to the next
        // one's scale, not to this one's cells.
        let running = running_m.vmax(running_i).vmax(running_d);
        shift2 = shift1;
        shift1 = normalising_shift_f32(running.horizontal_max());

        if len > 0 && lo + len - 1 == r {
            let slot = r - lo;
            let value = f64::from(get(sink.m.track, slot) + get(sink.i.track, slot));
            let base = *reference_exponent.get_or_insert(exponent);
            accumulator += value * exp2_f64(base - exponent);
        }

        lo_history[phase] = lo;
    }

    match reference_exponent {
        Some(base) if accumulator > 0.0 => {
            Log10Likelihood::new(accumulator.log10() - f64::from(base) * core::f64::consts::LOG10_2)
        }
        _ => Log10Likelihood::IMPOSSIBLE,
    }
}

#[inline]
fn get(buffer: &[f32], index: usize) -> f32 {
    buffer.get(index).copied().unwrap_or(0.0)
}

#[cfg(test)]
mod tests {
    use super::{ColumnLanes, Lane, MASK_SET, Plan, RowLanes, prior};
    use crate::{
        BaseQuality, Betas, ConversionModel, Emission, Haplotype, Probability, Read,
        StandardEmission, Strand, TapsEmission,
    };
    use proptest::prelude::*;
    use wide::f32x8;

    proptest! {
        /// The vector lane operations the strip kernel adds, against the
        /// scalar lane applied lane by lane: `shift_in` moves every lane up
        /// by one and drops the one `last` reports, and a mask `both` another
        /// mask is the lanewise and.
        #[test]
        fn vector_lane_ops_are_the_scalar_ops_lanewise(
            values in proptest::array::uniform8(-1e30f32..1e30),
            first in -1e30f32..1e30,
            left in proptest::array::uniform8(any::<bool>()),
            right in proptest::array::uniform8(any::<bool>()),
        ) {
            let vector = f32x8::from(values);
            let shifted = vector.shift_in(first).to_array();
            let [v0, v1, v2, v3, v4, v5, v6, v7] = values;
            prop_assert_eq!(shifted, [first, v0, v1, v2, v3, v4, v5, v6]);
            prop_assert_eq!(vector.last().to_bits(), v7.to_bits());
            prop_assert_eq!(<f32 as Lane>::shift_in(v0, first).to_bits(), first.to_bits());
            prop_assert_eq!(<f32 as Lane>::last(v0).to_bits(), v0.to_bits());

            let mask = |bits: [bool; 8]| f32x8::from(bits.map(|b| if b { MASK_SET } else { 0.0 }));
            let got = mask(left).both(mask(right)).to_array();
            for lane in 0..8 {
                let want = <f32 as Lane>::both(
                    if left[lane] { MASK_SET } else { 0.0 },
                    if right[lane] { MASK_SET } else { 0.0 },
                );
                prop_assert_eq!(got[lane].to_bits(), want.to_bits(), "lane {}", lane);
            }
        }

        /// `horizontal_sum` with one live lane is that lane exactly, which is
        /// the only way the kernel uses it.
        #[test]
        fn horizontal_sum_of_one_live_lane_is_that_lane(
            value in -1e30f32..1e30,
            lane in 0usize..8,
        ) {
            let mut lanes = [0.0f32; 8];
            *lanes.get_mut(lane).expect("lane in 0..8") = value;
            prop_assert_eq!(f32x8::from(lanes).horizontal_sum().to_bits(), value.to_bits());
        }
    }

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

    /// The scalar lanes of one cell, read straight off the plan's tracks.
    fn lanes_at(
        plan: &Plan,
        column: usize,
        row: usize,
    ) -> Option<(ColumnLanes<f32>, RowLanes<f32>)> {
        let at = |track: &[f32], index: usize| track.get(index).copied();
        let (columns, rows) = (&plan.columns, &plan.rows);
        Some((
            ColumnLanes {
                base: at(&columns.base, column)?,
                converted: at(&columns.converted, column)?,
                plain: at(&columns.plain, column)?,
                rate: at(&columns.rate, column)?,
                unconverted_rate: at(&columns.unconverted, column)?,
            },
            RowLanes::new(
                at(&rows.base, row)?,
                at(&rows.matched, row)?,
                at(&rows.mismatched, row)?,
            ),
        ))
    }

    fn check<E: Emission>(haplotype: &Haplotype, read: &Read, emission: &E) -> bool {
        let mut plan = Plan::default();
        plan.fill(haplotype, read, emission).expect("the fixture is well formed");
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
                let Some((column, row)) =
                    lanes_at(&plan, super::COLUMN_FRONT + last - column_index, row_index)
                else {
                    return false;
                };
                let got = prior::<f32>(column, row);
                let slack = 4.0 * f32::EPSILON * want.abs();
                if (got - want).abs() > slack {
                    return false;
                }
            }
        }
        true
    }

    /// A workspace reused across pairs of different sizes gives the same
    /// numbers as a fresh one: nothing from a previous alignment leaks.
    #[test]
    fn a_reused_workspace_is_a_fresh_workspace() {
        let long = Haplotype::from_ascii(b"ACGTTAGCATCGGATCCGATTACAGGCATTACGGATCCAGT");
        let short = Haplotype::from_ascii(b"TTAGCAT");
        let make = |bases: &[u8]| {
            Read::uniform(
                Haplotype::from_ascii(bases).bases().to_vec(),
                &vec![BaseQuality::from_byte(30); bases.len()],
                BaseQuality::from_byte(45),
                BaseQuality::from_byte(45),
                BaseQuality::from_byte(10),
                Strand::OT,
            )
            .expect("valid")
        };
        let reads = [make(b"TAGCATCGGATCCGATTACAGG"), make(b"AGCAT"), make(b"GGATCCAGTAAAA")];
        let emission = StandardEmission::default();
        let mut workspace = super::Workspace::new();
        for _ in 0..2 {
            for haplotype in [&long, &short, &long] {
                for read in &reads {
                    for band in [super::Band::anchored(2), super::Band::anchored(20)] {
                        let fresh = super::align_banded_simd(haplotype, read, &emission, band);
                        let reused = workspace.align_banded_simd(haplotype, read, &emission, band);
                        assert_eq!(fresh.get().to_bits(), reused.get().to_bits());
                        let fresh = super::align_banded(haplotype, read, &emission, band);
                        let reused = workspace.align_banded(haplotype, read, &emission, band);
                        assert_eq!(fresh.get().to_bits(), reused.get().to_bits());
                        let fresh = crate::align_strips_simd(haplotype, read, &emission, band);
                        let reused = workspace.align_strips_simd(haplotype, read, &emission, band);
                        assert_eq!(fresh.get().to_bits(), reused.get().to_bits());
                        let fresh = crate::align_strips(haplotype, read, &emission, band);
                        let reused = workspace.align_strips(haplotype, read, &emission, band);
                        assert_eq!(fresh.get().to_bits(), reused.get().to_bits());
                    }
                }
            }
        }
    }
}
