//! Genomic position newtype parameterized by coordinate system.
//!
//! `Pos<Zero>` is 0-based (BAM, BED, internal engine). `Pos<One>` is 1-based
//! (SAM, VCF, CRAM, user-facing). The type system prevents mixing coordinate
//! systems at compile time. `Offset` represents a signed distance between positions.
//!
//! # Niche optimization
//!
//! The internal storage is `NonMaxU32`, so `u32::MAX` is reserved as the niche.
//! This makes `Option<Pos<S>>` the same size as `Pos<S>` (4 bytes).
//! The maximum valid position value is `u32::MAX - 1`.

use std::fmt;
use std::marker::PhantomData;
use std::ops::Sub;

use nonmax::NonMaxU32;

// r[impl pos.systems]
/// 0-based coordinate system (BAM binary, BED, internal engine).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Zero;

// r[impl pos.systems]
/// 1-based coordinate system (SAM text, VCF, CRAM, user-facing).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct One;

// r[impl pos.type]
// r[impl pos.size]
// r[impl pos.incompatible]
// r[impl pos.derives]
// r[impl pos.niche]
/// A genomic position parameterized by coordinate system.
///
/// Zero runtime overhead: identical layout to `u32` via `#[repr(transparent)]`.
/// The phantom type parameter prevents mixing 0-based and 1-based positions
/// at compile time.
///
/// `u32::MAX` is reserved as the niche so that `Option<Pos<S>>` costs no extra
/// space. The maximum valid position value is `u32::MAX - 1`.
///
/// # Construction
///
/// ```
/// use seqair_types::pos::{Pos, Zero, One};
///
/// let bam_pos = Pos::<Zero>::new(100).unwrap();  // 0-based position 100
/// let sam_pos = Pos::<One>::new(101).unwrap();    // 1-based position 101
/// assert_eq!(bam_pos, sam_pos.to_zero_based()); // same genomic location
/// ```
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(transparent)]
pub struct Pos<S> {
    value: NonMaxU32,
    _system: PhantomData<S>,
}

// r[impl pos.offset]
/// Signed distance between two positions.
///
/// `Pos - Pos = Offset` and `Pos + Offset = Pos`. You cannot add two positions.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(transparent)]
pub struct Offset(i64);

// ---- Pos<Zero> construction ----

impl Pos<Zero> {
    // r[impl pos.zero_new]
    /// Create a 0-based position from a `u32`. Returns `None` if `value == u32::MAX`
    /// (reserved as the `Option<Pos>` niche).
    #[inline]
    pub fn new(value: u32) -> Option<Self> {
        NonMaxU32::new(value).map(|v| Self { value: v, _system: PhantomData })
    }

    // r[impl pos.try_from]
    /// Create a 0-based position from an `i64`. Returns `None` if negative,
    /// > `u32::MAX - 1`, or exactly `u32::MAX`.
    #[inline]
    pub fn try_from_i64(value: i64) -> Option<Self> {
        let v = u32::try_from(value).ok()?;
        NonMaxU32::new(v).map(|v| Self { value: v, _system: PhantomData })
    }

    // r[impl pos.try_from]
    /// Create a 0-based position from a `u64`. Returns `None` if > `u32::MAX - 1`
    /// or exactly `u32::MAX`.
    #[inline]
    pub fn try_from_u64(value: u64) -> Option<Self> {
        let v = u32::try_from(value).ok()?;
        NonMaxU32::new(v).map(|v| Self { value: v, _system: PhantomData })
    }

    // r[impl pos.to_one_based]
    // r[impl pos.explicit_conversion]
    /// Convert to 1-based. Returns `None` if the result would be `u32::MAX` (reserved niche).
    ///
    /// Only fails when `self` holds `u32::MAX - 1` (the maximum valid `Pos<Zero>`,
    /// since `u32::MAX` itself is the niche). Adding 1 would produce `u32::MAX`,
    /// which is reserved. This value is far beyond any real chromosome length.
    #[inline]
    pub fn to_one_based(self) -> Option<Pos<One>> {
        let new_val = self.value.get().checked_add(1)?;
        NonMaxU32::new(new_val).map(|v| Pos { value: v, _system: PhantomData })
    }

    // r[impl pos.try_from]
    /// Create a 0-based position from an `i32`. Returns `None` if negative.
    #[inline]
    pub fn try_from_i32(value: i32) -> Option<Self> {
        if value < 0 {
            return None;
        }
        Self::new(value as u32)
    }

    /// Maximum valid 0-based position (u32::MAX - 1).
    /// Used as "end of contig" sentinel in queries.
    pub fn max_value() -> Self {
        // u32::MAX - 1 is always valid for NonMaxU32
        Self::new(u32::MAX - 1).expect("BUG: u32::MAX - 1 is a valid NonMaxU32 value")
    }
}

// ---- Pos<One> construction ----

impl Pos<One> {
    // r[impl pos.one_new]
    /// Create a 1-based position from a `u32`. Returns `None` if value is 0 or `u32::MAX`.
    #[inline]
    pub fn new(value: u32) -> Option<Self> {
        if value == 0 {
            return None;
        }
        NonMaxU32::new(value).map(|v| Self { value: v, _system: PhantomData })
    }

    // r[impl pos.try_from]
    /// Create a 1-based position from an `i64`. Returns `None` if < 1, > `u32::MAX - 1`,
    /// or exactly `u32::MAX`.
    #[inline]
    pub fn try_from_i64(value: i64) -> Option<Self> {
        if value < 1 {
            return None;
        }
        u32::try_from(value).ok().and_then(Self::new)
    }

    // r[impl pos.try_from]
    /// Create a 1-based position from an `i32`. Returns `None` if < 1.
    #[inline]
    pub fn try_from_i32(value: i32) -> Option<Self> {
        if value < 1 {
            return None;
        }
        Self::new(value as u32)
    }

    // r[impl pos.to_zero_based]
    // r[impl pos.explicit_conversion]
    /// Convert to 0-based. Infallible: 1-based 1 → 0-based 0.
    ///
    /// The result of subtracting 1 from a value in 1..u32::MAX is always in
    /// 0..(u32::MAX-1), which is always a valid NonMaxU32.
    #[inline]
    pub fn to_zero_based(self) -> Pos<Zero> {
        // 1-based min is 1, so result is >= 0 and < u32::MAX — NonMaxU32 invariant always holds.
        NonMaxU32::new(self.value.get() - 1)
            .map(|v| Pos { value: v, _system: PhantomData })
            .expect("BUG: to_zero_based produced u32::MAX (impossible for valid Pos<One>)")
    }

    /// Maximum valid 1-based position (u32::MAX - 1).
    pub fn max_value() -> Self {
        // u32::MAX - 1 is non-zero and < u32::MAX, always valid
        Self::new(u32::MAX - 1).expect("BUG: u32::MAX - 1 is a valid Pos<One> value")
    }
}

// ---- Common methods (both systems) ----

impl<S> Pos<S> {
    // r[impl pos.get]
    // r[impl pos.must_use]
    /// Raw u32 value in the position's native coordinate system.
    #[inline]
    #[must_use]
    pub fn get(self) -> u32 {
        self.value.get()
    }

    // r[impl pos.as_usize]
    /// Convenience for indexing: returns the raw value as usize.
    #[inline]
    #[must_use]
    pub fn as_usize(self) -> usize {
        self.value.get() as usize
    }

    // r[impl pos.as_i64]
    /// Convenience for wider arithmetic: returns the raw value as i64.
    #[inline]
    #[must_use]
    pub fn as_i64(self) -> i64 {
        self.value.get() as i64
    }

    /// Convenience for 64-bit arithmetic: returns the raw value as u64.
    #[inline]
    #[must_use]
    pub fn as_u64(self) -> u64 {
        u64::from(self.value.get())
    }

    // r[impl pos.add_offset]
    /// Checked position + offset. Returns `None` if result is negative, exceeds u32 range, or hits the niche.
    #[inline]
    pub fn checked_add_offset(self, offset: Offset) -> Option<Self> {
        let result = self.value.get() as i64 + offset.0;
        let result_u32 = u32::try_from(result).ok()?;
        NonMaxU32::new(result_u32).map(|v| Self { value: v, _system: PhantomData })
    }

    // r[impl pos.sub_offset]
    /// Checked position - offset. Returns `None` if result is negative, exceeds u32 range, or hits the niche.
    #[inline]
    pub fn checked_sub_offset(self, offset: Offset) -> Option<Self> {
        self.checked_add_offset(Offset(-offset.0))
    }
}

// ---- Arithmetic ----

// r[impl pos.sub_pos]
// r[impl pos.no_add_pos]
// Pos - Pos = Offset (same system only); Add<Pos> is deliberately not implemented.
impl<S> Sub for Pos<S> {
    type Output = Offset;
    #[inline]
    fn sub(self, rhs: Self) -> Offset {
        Offset(self.value.get() as i64 - rhs.value.get() as i64)
    }
}

impl Offset {
    /// Construct an `Offset` from a signed 64-bit integer.
    #[inline]
    pub const fn new(value: i64) -> Self {
        Self(value)
    }

    /// Absolute value as usize (for lengths, capacities).
    #[inline]
    #[must_use]
    pub const fn abs_usize(self) -> usize {
        self.0.unsigned_abs() as usize
    }

    /// Raw i64 value.
    #[inline]
    #[must_use]
    pub const fn get(self) -> i64 {
        self.0
    }

    /// Checked addition. Returns `None` on overflow.
    #[inline]
    pub fn checked_add(self, rhs: Self) -> Option<Self> {
        self.0.checked_add(rhs.0).map(Offset)
    }

    /// Checked subtraction. Returns `None` on overflow.
    #[inline]
    pub fn checked_sub(self, rhs: Self) -> Option<Self> {
        self.0.checked_sub(rhs.0).map(Offset)
    }
}

// ---- From impls ----

impl<S> From<Pos<S>> for u32 {
    fn from(pos: Pos<S>) -> u32 {
        pos.value.get()
    }
}

impl<S> From<Pos<S>> for u64 {
    fn from(pos: Pos<S>) -> u64 {
        u64::from(pos.value.get())
    }
}

impl<S> From<Pos<S>> for i64 {
    fn from(pos: Pos<S>) -> i64 {
        i64::from(pos.value.get())
    }
}

impl<S> From<Pos<S>> for usize {
    fn from(pos: Pos<S>) -> usize {
        pos.value.get() as usize
    }
}

// ---- Display / Debug ----

impl fmt::Debug for Pos<Zero> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Pos0({})", self.value.get())
    }
}

impl fmt::Debug for Pos<One> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Pos1({})", self.value.get())
    }
}

impl<S> fmt::Display for Pos<S> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.value)
    }
}

impl fmt::Debug for Offset {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Offset({})", self.0)
    }
}

impl fmt::Display for Offset {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

#[cfg(test)]
impl Pos<Zero> {
    /// Test-only convenience: panics if value is u32::MAX.
    pub fn at(value: u32) -> Self {
        Self::new(value).expect("test position must not be u32::MAX")
    }
}

#[cfg(test)]
impl Pos<One> {
    /// Test-only convenience: panics if value is 0 or u32::MAX.
    pub fn at(value: u32) -> Self {
        Self::new(value).expect("test position must be 1..u32::MAX")
    }
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;

    // r[verify pos.to_one_based]
    // r[verify pos.to_zero_based]
    #[test]
    fn zero_based_roundtrip() {
        let z = Pos::<Zero>::new(100).unwrap();
        let o = z.to_one_based().unwrap();
        assert_eq!(o.get(), 101);
        assert_eq!(o.to_zero_based(), z);
    }

    // r[verify pos.to_zero_based]
    // r[verify pos.to_one_based]
    #[test]
    fn one_based_roundtrip() {
        let o = Pos::<One>::new(1).unwrap();
        let z = o.to_zero_based();
        assert_eq!(z.get(), 0);
        assert_eq!(z.to_one_based().unwrap(), o);
    }

    // r[verify pos.sub_pos]
    // r[verify pos.add_offset]
    #[test]
    fn pos_minus_pos_is_offset() {
        let a = Pos::<Zero>::new(100).unwrap();
        let b = Pos::<Zero>::new(50).unwrap();
        let off = a - b;
        assert_eq!(off.get(), 50);
        assert_eq!(b.checked_add_offset(off).unwrap(), a);
    }

    // r[verify pos.sub_pos]
    #[test]
    fn pos_minus_pos_negative_offset() {
        let a = Pos::<Zero>::new(10).unwrap();
        let b = Pos::<Zero>::new(50).unwrap();
        let off = a - b;
        assert_eq!(off.get(), -40);
    }

    // r[verify pos.add_offset]
    #[test]
    fn pos_plus_offset() {
        let p = Pos::<Zero>::new(10).unwrap();
        let q = p.checked_add_offset(Offset::new(5)).unwrap();
        assert_eq!(q.get(), 15);
    }

    // r[verify pos.sub_offset]
    #[test]
    fn pos_minus_offset() {
        let p = Pos::<Zero>::new(10).unwrap();
        let q = p.checked_sub_offset(Offset::new(3)).unwrap();
        assert_eq!(q.get(), 7);
    }

    // r[verify pos.try_from]
    #[test]
    fn try_from_i64_rejects_negative() {
        assert!(Pos::<Zero>::try_from_i64(-1).is_none());
        assert!(Pos::<One>::try_from_i64(0).is_none());
        assert!(Pos::<One>::try_from_i64(-1).is_none());
    }

    // r[verify pos.try_from]
    // r[verify pos.niche]
    #[test]
    fn try_from_i64_rejects_u32_max() {
        assert!(Pos::<Zero>::try_from_i64(u32::MAX as i64).is_none());
        assert!(Pos::<One>::try_from_i64(u32::MAX as i64).is_none());
    }

    // r[verify pos.try_from]
    #[test]
    fn try_from_i64_accepts_valid() {
        assert_eq!(Pos::<Zero>::try_from_i64(0).unwrap().get(), 0);
        assert_eq!(Pos::<One>::try_from_i64(1).unwrap().get(), 1);
        assert_eq!(Pos::<Zero>::try_from_i64(100).unwrap().get(), 100);
    }

    // r[verify pos.try_from]
    #[test]
    fn try_from_u64_rejects_overflow() {
        assert!(Pos::<Zero>::try_from_u64(u64::from(u32::MAX) + 1).is_none());
    }

    // r[verify pos.try_from]
    // r[verify pos.niche]
    #[test]
    fn try_from_u64_rejects_u32_max() {
        assert!(Pos::<Zero>::try_from_u64(u32::MAX as u64).is_none());
    }

    // r[verify pos.derives]
    #[test]
    fn ordering() {
        let a = Pos::<Zero>::new(10).unwrap();
        let b = Pos::<Zero>::new(20).unwrap();
        assert!(a < b);
        assert!(b > a);
    }

    // r[verify pos.as_usize]
    #[test]
    fn as_usize() {
        let p = Pos::<Zero>::new(42).unwrap();
        assert_eq!(p.as_usize(), 42);
    }

    // r[verify pos.as_i64]
    #[test]
    fn as_i64() {
        let p = Pos::<Zero>::new(100).unwrap();
        assert_eq!(p.as_i64(), 100);
    }

    // r[verify pos.size]
    #[test]
    fn size_is_u32() {
        assert_eq!(std::mem::size_of::<Pos<Zero>>(), std::mem::size_of::<u32>());
        assert_eq!(std::mem::size_of::<Pos<One>>(), std::mem::size_of::<u32>());
    }

    // r[verify pos.niche]
    #[test]
    fn option_pos_same_size_as_pos() {
        assert_eq!(std::mem::size_of::<Option<Pos<Zero>>>(), std::mem::size_of::<Pos<Zero>>(),);
        assert_eq!(std::mem::size_of::<Option<Pos<Zero>>>(), 4);
    }

    // r[verify pos.derives]
    #[test]
    fn debug_shows_coordinate_system() {
        let z = Pos::<Zero>::new(42).unwrap();
        let o = Pos::<One>::new(43).unwrap();
        assert_eq!(format!("{z:?}"), "Pos0(42)");
        assert_eq!(format!("{o:?}"), "Pos1(43)");
    }

    // r[verify pos.niche]
    #[test]
    fn max_value_zero() {
        let m = Pos::<Zero>::max_value();
        assert_eq!(m.get(), u32::MAX - 1);
    }

    // r[verify pos.niche]
    #[test]
    fn max_value_one() {
        let m = Pos::<One>::max_value();
        assert_eq!(m.get(), u32::MAX - 1);
    }

    // r[verify pos.niche]
    #[test]
    fn zero_based_new_rejects_niche() {
        assert!(Pos::<Zero>::new(u32::MAX).is_none(), "u32::MAX should be rejected (niche)");
    }

    // r[verify pos.niche]
    #[test]
    fn zero_based_new_accepts_max_minus_one() {
        assert!(Pos::<Zero>::new(u32::MAX - 1).is_some(), "u32::MAX - 1 should be accepted");
    }

    // r[verify pos.one_new]
    #[test]
    fn one_based_new_rejects_zero() {
        assert!(Pos::<One>::new(0).is_none(), "0 should be rejected for 1-based");
    }

    // r[verify pos.niche]
    #[test]
    fn one_based_new_rejects_niche() {
        assert!(Pos::<One>::new(u32::MAX).is_none(), "u32::MAX should be rejected (niche)");
    }

    // r[verify pos.to_one_based]
    #[test]
    fn to_one_based_returns_none_at_boundary() {
        // u32::MAX - 1 is the max valid Pos<Zero>. to_one_based would give u32::MAX (niche).
        let p = Pos::<Zero>::new(u32::MAX - 1).unwrap();
        assert!(p.to_one_based().is_none(), "to_one_based at niche boundary should return None");
    }

    // r[verify pos.add_offset]
    #[test]
    fn checked_add_offset_rejects_overflow() {
        // p = u32::MAX - 2. p + 1 = u32::MAX - 1 (valid). p + 2 = u32::MAX (niche, invalid).
        let p = Pos::<Zero>::new(u32::MAX - 2).unwrap();
        assert!(p.checked_add_offset(Offset::new(1)).is_some());
        assert!(p.checked_add_offset(Offset::new(2)).is_none(), "result would be u32::MAX (niche)");
    }

    // r[verify pos.add_offset]
    #[test]
    fn checked_add_offset_rejects_negative_result() {
        let p = Pos::<Zero>::new(5).unwrap();
        assert!(p.checked_add_offset(Offset::new(-10)).is_none(), "result would be negative");
    }

    // r[verify pos.sub_offset]
    #[test]
    fn checked_sub_offset_works() {
        let p = Pos::<Zero>::new(10).unwrap();
        assert_eq!(p.checked_sub_offset(Offset::new(3)).unwrap().get(), 7);
        assert!(p.checked_sub_offset(Offset::new(11)).is_none());
    }

    proptest! {
        // r[verify pos.to_one_based]
        // r[verify pos.to_zero_based]
        #[test]
        fn roundtrip_zero_to_one_and_back(v in 0u32..=u32::MAX - 2) {
            // u32::MAX - 1 maps to u32::MAX which is the niche, so to_one_based returns None.
            // Limit to u32::MAX - 2 for the roundtrip test.
            let z = Pos::<Zero>::new(v).unwrap();
            let o = z.to_one_based().unwrap();
            prop_assert_eq!(o.to_zero_based(), z);
        }

        // r[verify pos.to_zero_based]
        // r[verify pos.to_one_based]
        #[test]
        fn roundtrip_one_to_zero_and_back(v in 1u32..=u32::MAX - 1) {
            let o = Pos::<One>::new(v).unwrap();
            let z = o.to_zero_based();
            prop_assert_eq!(z.to_one_based().unwrap(), o);
        }

        // r[verify pos.add_offset]
        // r[verify pos.sub_offset]
        #[test]
        fn pos_plus_minus_offset_roundtrip(
            v in 0u32..1_000_000,
            off in -500_000i64..=500_000,
        ) {
            let result = v as i64 + off;
            if result >= 0 && result < u32::MAX as i64 {
                let p = Pos::<Zero>::new(v).unwrap();
                let q = p.checked_add_offset(Offset::new(off)).unwrap();
                let r = q.checked_sub_offset(Offset::new(off)).unwrap();
                prop_assert_eq!(r, p);
            }
        }

        // r[verify pos.sub_pos]
        #[test]
        fn pos_sub_pos_is_offset(a in 0u32..1_000_000, b in 0u32..1_000_000) {
            let pa = Pos::<Zero>::new(a).unwrap();
            let pb = Pos::<Zero>::new(b).unwrap();
            let off = pa - pb;
            prop_assert_eq!(off.get(), a as i64 - b as i64);
        }

        // r[verify pos.niche]
        #[test]
        fn option_pos_niche_always_4_bytes(v in 0u32..u32::MAX - 1) {
            let p = Pos::<Zero>::new(v);
            assert_eq!(std::mem::size_of_val(&p), 4);
        }

        // r[verify pos.niche]
        #[test]
        fn new_never_accepts_u32_max(v in u32::MAX..=u32::MAX) {
            prop_assert!(Pos::<Zero>::new(v).is_none());
            prop_assert!(Pos::<One>::new(v).is_none());
        }

        // r[verify pos.add_offset]
        #[test]
        fn checked_add_never_produces_u32_max(
            v in 0u32..u32::MAX - 1,
            off in 0i64..=100,
        ) {
            let p = Pos::<Zero>::new(v).unwrap();
            if let Some(result) = p.checked_add_offset(Offset::new(off)) {
                prop_assert_ne!(result.get(), u32::MAX, "checked_add must never produce u32::MAX");
            }
        }
    }
}
