use std::{
    fmt,
    ops::{self, Deref},
};

/// The root mean square (RMS) of a set of values.
///
/// RMS is a statistical measure of the magnitude of a varying quantity.
///
/// # Examples
///
/// ```rust
/// # use seqair_types::RootMeanSquare;
/// let data = [1, 2, 3, 4, 5];
/// // Explicitly construct from an iterator
/// let rms = RootMeanSquare::from_iter(data);
/// // our using `collect`
/// let rms: RootMeanSquare = data.into_iter().collect();
/// // you can use the value as a float
/// assert_eq!(rms.round(), 3.0);
/// ```
#[derive(Clone, Copy, Default, serde::Serialize, serde::Deserialize)]
#[must_use]
pub struct RootMeanSquare(f64);

impl<T: Into<f64>> FromIterator<T> for RootMeanSquare {
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        let mut acc = RmsAccumulator::new();
        for x in iter {
            acc.add(x.into());
        }
        acc.finish()
    }
}

/// Incremental accumulator for computing RMS in a single pass.
///
/// Useful when multiple RMS values need to be computed from the same data
/// without iterating multiple times.
#[derive(Clone, Copy, Default, Debug)]
pub struct RmsAccumulator {
    sum_of_squares: f64,
    count: u32,
}

impl RmsAccumulator {
    /// Creates a new empty accumulator.
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds a value to the accumulator.
    #[inline]
    pub fn add(&mut self, x: f64) {
        self.add_squared(x.algebraic_mul(x));
    }

    /// Adds a pre-computed squared value to the accumulator.
    ///
    /// Use this when the square of `x` is already available to avoid redundant multiplication.
    #[inline]
    pub fn add_squared(&mut self, x_sq: f64) {
        self.sum_of_squares = self.sum_of_squares.algebraic_add(x_sq);
        self.count = self.count.saturating_add(1);
    }

    /// Computes the final RMS from accumulated values.
    pub fn finish(self) -> RootMeanSquare {
        RootMeanSquare::from_sum_of_squares(self.sum_of_squares, self.count)
    }
}

impl RootMeanSquare {
    /// Build from a sum of squares accumulated elsewhere, and how many values
    /// went into it.
    ///
    /// [`RmsAccumulator`] carries its own `count`, which is the right default
    /// for one accumulator but pure duplication for a caller keeping several
    /// that necessarily share a count — per-strand and per-allele tallies of the
    /// same reads, say. Such a caller can keep bare `f64` sums and come here at
    /// the end, instead of incrementing one counter per sum per value.
    ///
    /// A `count` of zero is not an error: it is an empty set, whose RMS is zero.
    #[must_use]
    #[inline]
    pub fn from_sum_of_squares(sum_of_squares: f64, count: u32) -> Self {
        if count == 0 {
            return Self(0.0);
        }
        let average = sum_of_squares.algebraic_div(f64::from(count));
        Self(average.sqrt())
    }
}

/// A running sum of squared values, for an RMS whose count the caller keeps.
///
/// [`RmsAccumulator`] stores a count beside its sum, which is the right default
/// when there is one accumulator. It becomes duplication as soon as a caller
/// keeps *several* sums over the same values — a per-strand and a per-allele
/// tally of one set of reads, say — because then one count describes many sums,
/// and incrementing a private copy of it per sum per value is measurable work.
/// Keep the count once and the sums bare, and finish with [`Self::finish`].
///
/// # Examples
///
/// The RMS of 3 and 4 is `sqrt((9 + 16) / 2)`:
///
/// ```
/// use seqair_types::SumOfSquares;
///
/// let mut sum = SumOfSquares::new();
/// sum.add(3.0);
/// sum.add(4.0);
/// assert_eq!(*sum.finish(2), 12.5_f64.sqrt());
/// ```
///
/// The case it exists for — several sums, one count:
///
/// ```
/// use seqair_types::SumOfSquares;
///
/// // (base quality, mapping quality) for three reads at one position
/// let reads = [(30.0, 60.0), (35.0, 60.0), (25.0, 12.0)];
/// let (mut baseq, mut mapq) = (SumOfSquares::new(), SumOfSquares::new());
/// for (b, m) in reads {
///     baseq.add(b);
///     mapq.add(m);
/// }
///
/// let depth = u32::try_from(reads.len()).expect("three fits");
/// assert_eq!(*baseq.finish(depth), (((900.0 + 1225.0 + 625.0) / 3.0) as f64).sqrt());
/// assert_eq!(*mapq.finish(depth), (((3600.0 + 3600.0 + 144.0) / 3.0) as f64).sqrt());
/// ```
///
/// It agrees with [`RmsAccumulator`] *exactly*, not approximately, so the two
/// can be mixed and one can replace the other without moving any output:
///
/// ```
/// use seqair_types::{RmsAccumulator, SumOfSquares};
///
/// let values = [37.0, 41.0, 12.0];
/// let mut accumulator = RmsAccumulator::new();
/// let mut sum = SumOfSquares::new();
/// for v in values {
///     accumulator.add(v);
///     sum.add(v);
/// }
/// assert_eq!(accumulator.finish().to_bits(), sum.finish(3).to_bits());
/// ```
///
/// A count of zero is an empty set, not an error, and its RMS is zero:
///
/// ```
/// # use seqair_types::SumOfSquares;
/// assert_eq!(*SumOfSquares::new().finish(0), 0.0);
/// ```
///
/// Partial sums combine, so the values can be split across passes or threads:
///
/// ```
/// use seqair_types::SumOfSquares;
///
/// let mut whole = SumOfSquares::new();
/// for v in [1.0, 2.0, 3.0, 4.0] {
///     whole.add(v);
/// }
///
/// let mut left = SumOfSquares::new();
/// let mut right = SumOfSquares::new();
/// left.add(1.0);
/// left.add(2.0);
/// right.add(3.0);
/// right.add(4.0);
/// left += right;
///
/// assert_eq!(whole.finish(4).to_bits(), left.finish(4).to_bits());
/// ```
/// Every method here is `#[inline]` on purpose. The point of the type is to be
/// updated in a consumer's innermost loop, and a consumer without LTO — which
/// is the default release profile — can only inline across the crate boundary
/// what is marked. Un-marking one of these turns a paired `addpd` back into a
/// function call.
#[derive(Clone, Copy, Default, Debug, PartialEq, PartialOrd)]
pub struct SumOfSquares(f64);

impl SumOfSquares {
    /// An empty sum.
    #[must_use]
    #[inline]
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds a value, squaring it.
    #[inline]
    pub fn add(&mut self, x: f64) {
        self.add_squared(x.algebraic_mul(x));
    }

    /// Adds an already-squared value.
    ///
    /// Use this when the square is available anyway — a caller feeding the same
    /// `x²` to several sums should square it once.
    #[inline]
    pub fn add_squared(&mut self, x_sq: f64) {
        self.0 = self.0.algebraic_add(x_sq);
    }

    /// The RMS of the `count` values added.
    ///
    /// `count` is the caller's to track; passing a different one than was added
    /// computes the mean over that many values, which is occasionally what you
    /// want and otherwise a bug this type cannot catch for you.
    #[inline]
    pub fn finish(self, count: u32) -> RootMeanSquare {
        RootMeanSquare::from_sum_of_squares(self.0, count)
    }

    /// The accumulated sum of squares itself.
    #[must_use]
    #[inline]
    pub fn get(self) -> f64 {
        self.0
    }
}

impl From<f64> for SumOfSquares {
    /// Resumes from a sum accumulated elsewhere.
    #[inline]
    fn from(sum_of_squares: f64) -> Self {
        Self(sum_of_squares)
    }
}

// Merging is `+=` and deliberately not `+`: an inherent `add(value)` and an
// `Add::add(other_sum)` would be two different meanings of the same word on one
// type, which is how a caller ends up squaring a sum by accident.
impl ops::AddAssign for SumOfSquares {
    #[inline]
    fn add_assign(&mut self, other: Self) {
        self.add_squared(other.0);
    }
}

impl Deref for RootMeanSquare {
    type Target = f64;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl ops::Mul<f64> for RootMeanSquare {
    type Output = f64;

    fn mul(self, rhs: f64) -> Self::Output {
        self.0 * rhs
    }
}

impl ops::Mul<RootMeanSquare> for f64 {
    type Output = f64;

    fn mul(self, rhs: RootMeanSquare) -> Self::Output {
        self * rhs.0
    }
}

#[cfg_attr(coverage_nightly, coverage(off))]
impl fmt::Debug for RootMeanSquare {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "RMS({:.2})", self.0)
    }
}

#[cfg_attr(coverage_nightly, coverage(off))]
impl fmt::Display for RootMeanSquare {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

/// Extension trait for iterators that adds an `.rms()` method.
///
/// This trait provides a convenient way to calculate the root mean square
/// of an iterator's items without explicitly calling `collect()`.
///
/// # Examples
///
/// ```rust
/// # use seqair_types::{RootMeanSquare, RootMeanSquareExt};
/// let data = [1, 2, 3, 4, 5];
/// let rms = data.into_iter().rms();
/// assert_eq!(rms.round(), 3.0);
/// ```
pub trait RootMeanSquareExt: Iterator {
    /// Calculates the root mean square of the iterator's items.
    ///
    /// This is equivalent to collecting into a `RootMeanSquare`.
    fn rms<T>(self) -> RootMeanSquare
    where
        Self: Sized + Iterator<Item = T>,
        T: Into<f64>,
    {
        RootMeanSquare::from_iter(self)
    }
}

impl<I: Iterator> RootMeanSquareExt for I {}

#[cfg(test)]
mod tests {
    use super::*;

    /// The whole promise of [`SumOfSquares`] is that dropping the count changes
    /// nothing about the answer. Pinned on the bits, not a tolerance: these
    /// values are printed into files that get diffed.
    #[test]
    fn sum_of_squares_matches_the_accumulator_bit_for_bit() {
        let cases: &[&[f64]] = &[
            &[],
            &[0.0],
            &[37.0],
            &[37.0, 41.0, 12.0],
            &[60.0; 1000],
            &[1e-8, 1e8, 3.5, 0.25],
            &[0.0, 0.0, 40.0],
        ];
        for values in cases {
            let mut accumulator = RmsAccumulator::new();
            let mut sum = SumOfSquares::new();
            for &v in *values {
                accumulator.add(v);
                sum.add(v);
            }
            let count = u32::try_from(values.len()).expect("test cases are small");
            assert_eq!(
                accumulator.finish().to_bits(),
                sum.finish(count).to_bits(),
                "diverged on {values:?}"
            );
        }
    }

    /// Splitting the values across two sums and merging must land in the same
    /// place as one pass, or the `AddAssign` the docs advertise is a trap.
    #[test]
    fn partial_sums_merge_exactly() {
        let mut whole = SumOfSquares::new();
        let (mut left, mut right) = (SumOfSquares::new(), SumOfSquares::new());
        for (i, v) in [3.0_f64, 1.5, 9.0, 0.25, 7.0].into_iter().enumerate() {
            whole.add(v);
            if i % 2 == 0 { left.add(v) } else { right.add(v) }
        }
        left += right;
        assert_eq!(whole.finish(5).to_bits(), left.finish(5).to_bits());
    }

    /// The accumulator now *is* this constructor, so the two can only agree —
    /// but that is the property callers rely on when they mix the styles, and
    /// it is one refactor away from silently not holding.
    #[test]
    fn from_sum_of_squares_agrees_with_the_accumulator() {
        for values in [&[][..], &[0.0], &[37.0], &[37.0, 41.0, 12.0], &[1e-8, 1e8, 3.5]] {
            let mut acc = RmsAccumulator::new();
            let mut sum = 0.0f64;
            for &v in values {
                acc.add(v);
                sum = sum.algebraic_add(v.algebraic_mul(v));
            }
            let count = u32::try_from(values.len()).expect("fits");
            assert_eq!(
                acc.finish().to_bits(),
                RootMeanSquare::from_sum_of_squares(sum, count).to_bits(),
                "diverged on {values:?}"
            );
        }
    }
    use proptest::{collection::vec, prelude::*};

    #[test]
    fn test_empty_accumulator_gives_zero() {
        let acc = RmsAccumulator::new();
        let rms = acc.finish();
        assert_eq!(*rms, 0.0);
    }

    proptest! {
        #[test]
        fn test_rms_constant_value(value: u8) {
            let data = vec![value; 100];
            let rms = RootMeanSquare::from_iter(data);
            prop_assert_eq!(*rms, f64::from(value));
        }

        #[test]
        fn test_rms_accumulator_identical_values(value in 0u8..=255, count in 1usize..200) {
            let mut acc = RmsAccumulator::new();
            for _ in 0..count {
                acc.add(f64::from(value));
            }
            let rms = acc.finish();
            // RMS of N identical values v is v
            let diff = (*rms - f64::from(value)).abs();
            prop_assert!(diff < 1e-10, "RMS of {} identical {value} = {}, expected {value}", count, *rms);
        }

        #[test]
        fn test_rms_never_negative(data: Vec<u8>) {
            let rms = RootMeanSquare::from_iter(data);
            prop_assert!(*rms >= 0.0);
        }

        #[test]
        fn test_rms_zero_iff_all_zeros(data: Vec<u8>) {
            let rms = RootMeanSquare::from_iter(data.clone());
            if data.iter().all(|&x| x == 0) {
                prop_assert_eq!(rms.0, 0.0);
            } else if !data.is_empty() {
                prop_assert!(*rms > 0.0);
            }
        }

        #[test]
        fn test_rms_greater_than_or_equal_to_mean(data in vec(any::<u8>(), 1..300)) {
            let mean = data.iter().map(|&x| f64::from(x)).sum::<f64>() / data.len() as f64;
            let rms = RootMeanSquare::from_iter(data);
            prop_assert!(*rms >= mean);
        }

        #[test]
        fn test_rms_less_than_or_equal_to_max(data in vec(any::<u8>(), 1..300)) {
            let max = f64::from(*data.iter().max().unwrap());
            let rms = RootMeanSquare::from_iter(data);
            prop_assert!(*rms <= max);
        }
    }
}
