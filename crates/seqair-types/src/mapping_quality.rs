//! Per-record mapping quality wrapper with explicit "unavailable" handling.
//!
//! See `docs/spec/0-1-mapping-quality.md` for the full contract.

use std::fmt;

/// A record's mapping quality (MAPQ) as stored on the wire (SAM, BAM, CRAM).
///
/// `255` is reserved by `[SAM1] §1.4` as "mapping quality not available".
/// Access the number via [`MappingQuality::get`], which returns `None` for
/// the sentinel. Intentionally does not implement `Ord`/`PartialOrd`, so a
/// minimum-MAPQ filter has to say what it does with an unavailable MAPQ —
/// see `r[types.mapping_quality.no_ord]`.
// r[impl types.mapping_quality.type]
// r[impl types.mapping_quality.copy]
// r[impl types.mapping_quality.no_ord]
#[repr(transparent)]
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct MappingQuality(u8);

const _: () = assert!(
    std::mem::size_of::<MappingQuality>() == std::mem::size_of::<u8>(),
    "MappingQuality must be 1 byte",
);
const _: () = assert!(
    std::mem::align_of::<MappingQuality>() == std::mem::align_of::<u8>(),
    "MappingQuality must have u8 alignment",
);

impl MappingQuality {
    /// The "mapping quality not available" sentinel from SAM/BAM/CRAM.
    // r[impl types.mapping_quality.unavailable]
    pub const UNAVAILABLE: Self = Self(255);

    /// Wrap a raw wire byte. All 256 values are valid.
    // r[impl types.mapping_quality.from_byte]
    // r[impl types.mapping_quality.roundtrip]
    #[inline]
    pub const fn from_byte(b: u8) -> Self {
        Self(b)
    }

    /// The MAPQ, or `None` if the wire byte signals "unavailable".
    // r[impl types.mapping_quality.get]
    #[inline]
    pub const fn get(self) -> Option<u8> {
        match self.0 {
            255 => None,
            q => Some(q),
        }
    }

    /// Raw wire byte, including the `255` sentinel. Escape hatch for
    /// serialization; prefer [`MappingQuality::get`] in consumer code.
    // r[impl types.mapping_quality.as_byte]
    #[inline]
    pub const fn as_byte(self) -> u8 {
        self.0
    }
}

#[cfg_attr(coverage_nightly, coverage(off))]
impl fmt::Debug for MappingQuality {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.get() {
            Some(q) => write!(f, "MappingQuality({q})"),
            None => f.write_str("MappingQuality(unavailable)"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // r[verify types.mapping_quality.get]
    // r[verify types.mapping_quality.unavailable]
    #[test]
    fn get_returns_none_only_for_the_sentinel() {
        assert_eq!(MappingQuality::UNAVAILABLE.get(), None);
        assert_eq!(MappingQuality::UNAVAILABLE.as_byte(), 255);
        for v in 0u8..255 {
            assert_eq!(MappingQuality::from_byte(v).get(), Some(v));
        }
    }

    // r[verify types.mapping_quality.roundtrip]
    // r[verify types.mapping_quality.from_byte]
    // r[verify types.mapping_quality.as_byte]
    #[test]
    fn from_byte_as_byte_roundtrip_all_values() {
        for v in 0u8..=255 {
            assert_eq!(MappingQuality::from_byte(v).as_byte(), v);
        }
    }

    // r[verify types.mapping_quality.copy]
    #[test]
    fn derives_copy_eq_hash_and_debug_names_unavailable() {
        use std::collections::HashSet;
        let a = MappingQuality::from_byte(60);
        let b = a;
        assert_eq!(a, b);
        assert!(HashSet::from([a]).contains(&b));
        assert_eq!(format!("{a:?}"), "MappingQuality(60)");
        assert_eq!(format!("{:?}", MappingQuality::UNAVAILABLE), "MappingQuality(unavailable)");
    }
}
