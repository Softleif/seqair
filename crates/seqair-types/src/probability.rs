use std::{num::ParseFloatError, ops::Deref, str::FromStr};

/// A probability value guaranteed to be in the range [0, 1]
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
#[must_use]
pub struct Probability(f64);

/// Errors that can occur when creating or parsing a [`Probability`].
#[derive(Debug, thiserror::Error)]
pub enum ProbabilityError {
    /// The value is outside the valid range [0, 1]
    #[error("Probability value `{value}` is outside the valid range [0, 1]")]
    OutOfRange {
        /// The invalid value
        value: f64,
    },
    /// Failed to parse a float from a string
    #[error("Failed to parse probability float")]
    ParseFloat(#[from] ParseFloatError),
}

impl Probability {
    /// Creates a new Probability, returning an error if the value is outside [0, 1]
    pub const fn new(value: f64) -> Result<Self, ProbabilityError> {
        if value < 0.0 || value > 1.0 || !value.is_finite() {
            Err(ProbabilityError::OutOfRange { value })
        } else {
            Ok(Self(value))
        }
    }

    /// Creates a new Probability, panicking if the value is outside [0, 1]
    ///
    /// Can be used to build Probability in const contexts.
    #[expect(clippy::panic, reason = "panicking is intentional here")]
    pub const fn new_panicky(value: f64) -> Self {
        match Self::new(value) {
            Ok(p) => p,
            Err(_e) => panic!("invalid probability"),
        }
    }

    /// Returns the inverted probability (1 - p)
    pub const fn inverted(self) -> Self {
        Self(1.0 - self.0)
    }

    /// The probability value 0.0
    pub const ZERO: Self = Self(0.0);
    /// The probability value 1.0
    pub const ONE: Self = Self(1.0);
}

impl std::fmt::Display for Probability {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if let Some(precision) = f.precision() {
            // If we received a precision, we use it.
            write!(f, "{1:.*}", precision, self.0)
        } else {
            // Otherwise we default to 2.
            write!(f, "{:.2}", self.0)
        }
    }
}

impl FromStr for Probability {
    type Err = ProbabilityError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let v = s.parse::<f64>()?;
        Self::new(v)
    }
}

impl Deref for Probability {
    type Target = f64;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl TryFrom<f64> for Probability {
    type Error = ProbabilityError;

    fn try_from(value: f64) -> Result<Self, Self::Error> {
        Self::new(value)
    }
}

impl Default for Probability {
    fn default() -> Self {
        Self(0.0)
    }
}

impl serde::Serialize for Probability {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        self.0.serialize(serializer)
    }
}

impl<'de> serde::Deserialize<'de> for Probability {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let value = f64::deserialize(deserializer)?;
        Self::new(value).map_err(serde::de::Error::custom)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fails_outside_range() {
        assert!(Probability::new(-0.1).is_err());
        assert!(Probability::new(1.1).is_err());
        assert!(Probability::new(f64::NAN).is_err());
        assert!(Probability::new(f64::INFINITY).is_err());
        assert!(Probability::new(f64::NEG_INFINITY).is_err());
    }

    // r[verify types.probability.r23_from_str_typed_error]
    #[test]
    fn from_str_returns_typed_error() {
        // Valid parse
        let p: Probability = "0.5".parse().expect("valid probability");
        assert_eq!(*p, 0.5);

        // ParseFloat variant for non-numeric input
        let err = "not_a_number".parse::<Probability>().unwrap_err();
        assert!(matches!(err, ProbabilityError::ParseFloat(_)));

        // OutOfRange variant for out-of-range value
        let err = "1.5".parse::<Probability>().unwrap_err();
        assert!(matches!(err, ProbabilityError::OutOfRange { value } if value == 1.5));
    }

    // r[verify types.probability.r23_new_returns_out_of_range]
    #[test]
    fn new_returns_out_of_range_variant() {
        let err = Probability::new(2.0).unwrap_err();
        assert!(matches!(err, ProbabilityError::OutOfRange { value } if value == 2.0));
    }

    // r[verify types.probability.r24_error_reexport]
    #[test]
    fn error_accessible_from_crate_root() {
        let err: crate::ProbabilityError = ProbabilityError::OutOfRange { value: 2.0 };
        assert!(matches!(err, crate::ProbabilityError::OutOfRange { .. }));
    }
}
