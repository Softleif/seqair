//! Per-base quality score wrapper with explicit "unavailable" handling.
//!
//! See `docs/spec/0-1-base-quality.md` for the full contract.

use std::fmt;

/// Per-base quality score as stored on the wire (BAM binary, SAM decoded).
///
/// This is distinct from [`crate::Phred`], which is an `f64` for probability
/// math. `BaseQuality` is the wire-level byte, and `0xFF` is reserved by
/// `[SAM1] §4.2.3` as the per-record "quality unavailable" sentinel.
///
/// Access the Phred integer via [`BaseQuality::get`], which returns `None`
/// for the sentinel. Intentionally does not implement `Ord`/`PartialOrd` —
/// see `r[types.base_quality.no_ord]`.
// r[impl types.base_quality.type]
// r[impl types.base_quality.copy]
// r[impl types.base_quality.no_ord]
#[repr(transparent)]
#[derive(Clone, Copy, PartialEq, Eq, Hash, bytemuck::Pod, bytemuck::Zeroable)]
pub struct BaseQuality(u8);

// r[impl types.base_quality.size]
const _: () = assert!(
    std::mem::size_of::<BaseQuality>() == std::mem::size_of::<u8>(),
    "BaseQuality must be 1 byte for the zero-copy slice cast",
);
const _: () = assert!(
    std::mem::align_of::<BaseQuality>() == std::mem::align_of::<u8>(),
    "BaseQuality must have u8 alignment for the zero-copy slice cast",
);

impl BaseQuality {
    /// The per-record "quality unavailable" sentinel from SAM/BAM.
    // r[impl types.base_quality.unavailable]
    pub const UNAVAILABLE: Self = Self(0xFF);

    /// Wrap a raw wire byte. All 256 values are valid.
    // r[impl types.base_quality.from_byte]
    // r[impl types.base_quality.roundtrip]
    #[inline]
    pub const fn from_byte(b: u8) -> Self {
        Self(b)
    }

    /// The Phred integer, or `None` if the wire byte signals "unavailable".
    // r[impl types.base_quality.get]
    #[inline]
    pub const fn get(self) -> Option<u8> {
        match self.0 {
            0xFF => None,
            q => Some(q),
        }
    }

    /// Raw wire byte, including the `0xFF` sentinel. Escape hatch for
    /// serialization; prefer [`BaseQuality::get`] in consumer code.
    // r[impl types.base_quality.as_byte]
    #[inline]
    pub const fn as_byte(self) -> u8 {
        self.0
    }

    /// Zero-copy view of a qual byte slab as `&[BaseQuality]`.
    ///
    /// Sound by `r[types.base_quality.type]` + `r[types.base_quality.size]`:
    /// `#[repr(transparent)]` over `u8` guarantees identical layout and
    /// alignment, and every byte value is a valid `BaseQuality` — encoded
    /// in the `Pod` derive on the type.
    // r[impl types.base_quality.slice_from_bytes]
    #[inline]
    pub fn slice_from_bytes(bytes: &[u8]) -> &[BaseQuality] {
        bytemuck::cast_slice(bytes)
    }

    /// Inverse of [`BaseQuality::slice_from_bytes`]: view `&[BaseQuality]` as
    /// raw bytes, preserving length and lifetime.
    // r[impl types.base_quality.slice_to_bytes]
    #[inline]
    pub fn slice_to_bytes(quals: &[BaseQuality]) -> &[u8] {
        bytemuck::cast_slice(quals)
    }
}

#[cfg_attr(coverage_nightly, coverage(off))]
impl fmt::Debug for BaseQuality {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.get() {
            Some(q) => write!(f, "BaseQuality({q})"),
            None => f.write_str("BaseQuality(unavailable)"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // r[verify types.base_quality.get]
    #[test]
    fn get_returns_none_for_sentinel() {
        assert_eq!(BaseQuality::UNAVAILABLE.get(), None);
        assert_eq!(BaseQuality::from_byte(0xFF).get(), None);
    }

    // r[verify types.base_quality.get]
    #[test]
    fn get_returns_some_for_phred_values() {
        assert_eq!(BaseQuality::from_byte(0).get(), Some(0));
        assert_eq!(BaseQuality::from_byte(30).get(), Some(30));
        assert_eq!(BaseQuality::from_byte(93).get(), Some(93));
        // 0xFE is the largest valid (non-sentinel) byte. Not a realistic Phred,
        // but BaseQuality is defined on the wire byte, not on Phred range.
        assert_eq!(BaseQuality::from_byte(0xFE).get(), Some(0xFE));
    }

    // r[verify types.base_quality.unavailable]
    #[test]
    fn unavailable_constant_matches_wire_sentinel() {
        assert_eq!(BaseQuality::UNAVAILABLE.as_byte(), 0xFF);
    }

    // r[verify types.base_quality.roundtrip]
    #[test]
    fn from_byte_as_byte_roundtrip_all_values() {
        for v in 0u8..=255 {
            assert_eq!(BaseQuality::from_byte(v).as_byte(), v, "roundtrip failed for {v}");
        }
    }

    // r[verify types.base_quality.slice_from_bytes]
    // r[verify types.base_quality.slice_to_bytes]
    #[test]
    fn slice_cast_preserves_length_and_bytes() {
        let bytes: Vec<u8> = (0u8..=200).collect();
        let quals = BaseQuality::slice_from_bytes(&bytes);
        assert_eq!(quals.len(), bytes.len());
        for (i, q) in quals.iter().enumerate() {
            #[allow(clippy::indexing_slicing, reason = "test; i < bytes.len() by construction")]
            let expected = bytes[i];
            assert_eq!(q.as_byte(), expected);
        }
        let back = BaseQuality::slice_to_bytes(quals);
        assert_eq!(back, bytes.as_slice());
        // Pointer identity: the cast must not have copied.
        assert_eq!(back.as_ptr(), bytes.as_ptr());
    }

    // r[verify types.base_quality.slice_from_bytes]
    #[test]
    fn slice_cast_handles_empty() {
        let empty: &[u8] = &[];
        let quals = BaseQuality::slice_from_bytes(empty);
        assert!(quals.is_empty());
    }

    // r[verify types.base_quality.slice_from_bytes]
    #[test]
    fn slice_cast_preserves_sentinel_and_mixed() {
        let bytes = [30u8, 0xFF, 25, 0xFF, 0, 93];
        let quals = BaseQuality::slice_from_bytes(&bytes);
        assert_eq!(quals[0].get(), Some(30));
        assert_eq!(quals[1].get(), None);
        assert_eq!(quals[2].get(), Some(25));
        assert_eq!(quals[3].get(), None);
        assert_eq!(quals[4].get(), Some(0));
        assert_eq!(quals[5].get(), Some(93));
    }

    // r[verify types.base_quality.copy]
    #[test]
    fn derives_copy_eq_hash() {
        use std::collections::HashSet;
        let a = BaseQuality::from_byte(30);
        let b = a; // Copy
        assert_eq!(a, b);
        let mut set = HashSet::new();
        set.insert(a);
        assert!(set.contains(&b));
    }

    #[test]
    fn debug_format_distinguishes_unavailable() {
        assert_eq!(format!("{:?}", BaseQuality::from_byte(30)), "BaseQuality(30)");
        assert_eq!(format!("{:?}", BaseQuality::UNAVAILABLE), "BaseQuality(unavailable)");
    }
}
