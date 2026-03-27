use crate::Probability;
use std::{fmt, ops::Deref};

/// Phred-scaled quality score
///
/// # Definition
///
/// Source: [Wikipedia](https://en.wikipedia.org/wiki/Phred_quality_score)
///
/// > Phred quality scores `Q` are logarithmically related to the base-calling error probabilities `P` and defined as
/// > `Q = -10 log_10 P`.
/// >
/// > This relation can also be written as
/// > `P = 10^(-Q / 10)`.
/// >
/// > For example, if Phred assigns a quality score of 30 to a base, the chances that this base is called incorrectly are 1 in 1000.
/// >
/// > | Phred Quality Score | Probability of incorrect base call | Base call accuracy |
/// > | --- | --- | --- |
/// > | 10 | 1 in 10 | 90% |
/// > | 20 | 1 in 100 | 99% |
/// > | 30 | 1 in 1000 | 99.9% |
/// > | 40 | 1 in 10,000 | 99.99% |
/// > | 50 | 1 in 100,000 | 99.999% |
/// > | 60 | 1 in 1,000,000 | 99.9999% |
/// >
/// > The phred quality score is the negative ratio of the error probability to the reference level of `P=1` expressed in [Decibel (dB)](https://en.wikipedia.org/wiki/20_log_rule "20 log rule").
#[derive(Clone, Copy, serde::Serialize, serde::Deserialize)]
#[must_use]
pub struct Phred(f64);

impl Phred {
    /// Get the Phred quality score as an integer, clamped to [0, 99].
    // r[impl types.phred.as_int_clamp]
    pub fn as_int(self) -> i32 {
        let phred = self.0;
        // NaN and values <= 0 clamp to 0
        if !(phred > 0.0) {
            return 0;
        }
        if phred >= 99.0 {
            return 99;
        }
        phred.round() as i32
    }

    /// Create a Phred quality score from an integer.
    // r[impl types.phred.non_negative]
    pub fn from_phred(phred: u8) -> Self {
        Self(f64::from(phred))
    }
}

impl From<Probability> for Phred {
    fn from(probability: Probability) -> Self {
        let phred = -10.0 * probability.log10();
        Self(phred)
    }
}

impl Deref for Phred {
    type Target = f64;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

#[cfg_attr(coverage_nightly, coverage(off))]
impl fmt::Debug for Phred {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Phred({:.2})", self.0)
    }
}

#[cfg_attr(coverage_nightly, coverage(off))]
impl fmt::Display for Phred {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_wikipedia_examples() {
        assert_eq!(10., *Phred::from(Probability::new_panicky(0.1)));
        assert_eq!(20., *Phred::from(Probability::new_panicky(0.01)));
        assert_eq!(30., *Phred::from(Probability::new_panicky(0.001)));
        assert_eq!(40., *Phred::from(Probability::new_panicky(0.0001)));
        assert_eq!(50., *Phred::from(Probability::new_panicky(0.00001)));
        assert_eq!(60., *Phred::from(Probability::new_panicky(0.000001)));
    }

    #[test]
    fn test_phred_zero() {
        // A probability of 1 (100% certainty) should yield a Phred score of 0
        assert_eq!(0., *Phred::from(Probability::new_panicky(1.)));
    }

    // r[verify types.phred.non_negative]
    #[test]
    fn from_phred_accepts_u8() {
        // from_phred takes u8, so negative values are unrepresentable
        let p = Phred::from_phred(0);
        assert_eq!(p.as_int(), 0);
        let p = Phred::from_phred(42);
        assert_eq!(p.as_int(), 42);
        let p = Phred::from_phred(255);
        assert_eq!(p.as_int(), 99); // clamped
    }

    // r[verify types.phred.as_int_clamp]
    #[test]
    fn as_int_clamps_low_to_zero() {
        // A probability of 1.0 yields Phred 0.0; as_int must return 0
        let p = Phred::from(Probability::new_panicky(1.0));
        assert_eq!(p.as_int(), 0);
    }

    // r[verify types.phred.as_int_clamp]
    #[test]
    fn as_int_clamps_high_to_99() {
        let p = Phred::from_phred(99);
        assert_eq!(p.as_int(), 99);
        let p = Phred::from_phred(100);
        assert_eq!(p.as_int(), 99);
    }

    // r[verify types.phred.as_int_clamp]
    #[test]
    fn as_int_handles_nan() {
        // Construct a Phred with NaN inner value (via From<Probability> edge case
        // is impossible since Probability rejects NaN, so we test the raw struct)
        let p = Phred(f64::NAN);
        assert_eq!(p.as_int(), 0);
    }
}
