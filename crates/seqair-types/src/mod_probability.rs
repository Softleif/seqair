//! A base-modification call's probability as stored in the `ML` tag.
//!
//! See `docs/spec/0-1-mod-probability.md` for the full contract.

use crate::Probability;
use std::fmt;

/// One entry of the `ML` tag (`[SAM1] §1.7`): the byte `x` says the modification
/// probability lies in `[x/256, (x+1)/256)`. It is a bucket, not a probability, so
/// `x/255` and `x/256` are both wrong; use [`ModProbability::midpoint`] for a number.
///
/// Ordered by the byte: a larger byte is a higher probability bucket.
// r[impl types.mod_probability.type]
// r[impl types.mod_probability.copy]
// r[impl types.mod_probability.ord]
#[repr(transparent)]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ModProbability(u8);

const _: () = assert!(
    std::mem::size_of::<ModProbability>() == std::mem::size_of::<u8>(),
    "ModProbability must be 1 byte",
);
const _: () = assert!(
    std::mem::align_of::<ModProbability>() == std::mem::align_of::<u8>(),
    "ModProbability must have u8 alignment",
);

impl ModProbability {
    /// Wrap a raw `ML` byte. All 256 values are valid.
    // r[impl types.mod_probability.from_byte]
    // r[impl types.mod_probability.roundtrip]
    #[inline]
    pub const fn from_byte(b: u8) -> Self {
        Self(b)
    }

    /// The raw `ML` byte, for serialization.
    // r[impl types.mod_probability.as_byte]
    #[inline]
    pub const fn as_byte(self) -> u8 {
        self.0
    }

    /// Midpoint of the bucket, `(x + 0.5) / 256`: never exactly 0 or 1.
    // r[impl types.mod_probability.midpoint]
    pub fn midpoint(self) -> Probability {
        // (x + 0.5) / 256 is in (0, 1) for every byte, so this cannot panic.
        Probability::new_panicky((f64::from(self.0) + 0.5) / 256.0)
    }
}

#[cfg_attr(coverage_nightly, coverage(off))]
impl fmt::Debug for ModProbability {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "ModProbability({})", self.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // r[verify types.mod_probability.roundtrip]
    // r[verify types.mod_probability.from_byte]
    // r[verify types.mod_probability.as_byte]
    #[test]
    fn from_byte_as_byte_roundtrip_all_values() {
        for v in 0u8..=255 {
            assert_eq!(ModProbability::from_byte(v).as_byte(), v);
        }
    }

    // r[verify types.mod_probability.midpoint]
    #[test]
    fn midpoint_is_the_bucket_centre_and_strictly_inside_zero_one() {
        assert_eq!(*ModProbability::from_byte(0).midpoint(), 0.5 / 256.0);
        assert_eq!(*ModProbability::from_byte(128).midpoint(), 128.5 / 256.0);
        assert_eq!(*ModProbability::from_byte(255).midpoint(), 255.5 / 256.0);
        for v in 0u8..=255 {
            let p = *ModProbability::from_byte(v).midpoint();
            assert!(p > 0.0 && p < 1.0);
        }
    }

    // r[verify types.mod_probability.ord]
    // r[verify types.mod_probability.copy]
    #[test]
    fn order_follows_the_byte_and_the_midpoint() {
        let (lo, hi) = (ModProbability::from_byte(127), ModProbability::from_byte(128));
        assert!(lo < hi);
        assert!(lo.midpoint() < hi.midpoint());
        assert_eq!(format!("{lo:?}"), "ModProbability(127)");
    }
}
