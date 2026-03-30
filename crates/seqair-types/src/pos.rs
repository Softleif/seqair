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
use std::ops::{Add, AddAssign, Sub, SubAssign};

use nonmax::NonMaxU32;

/// 0-based coordinate system (BAM binary, BED, internal engine).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Zero;

/// 1-based coordinate system (SAM text, VCF, CRAM, user-facing).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct One;

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
/// let bam_pos = Pos::<Zero>::new(100);       // 0-based position 100
/// let sam_pos = Pos::<One>::new(101);         // 1-based position 101
/// assert_eq!(bam_pos, sam_pos.to_zero_based()); // same genomic location
/// ```
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(transparent)]
pub struct Pos<S> {
    value: NonMaxU32,
    _system: PhantomData<S>,
}

/// Signed distance between two positions.
///
/// `Pos - Pos = Offset` and `Pos + Offset = Pos`. You cannot add two positions.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(transparent)]
pub struct Offset(pub i64);

// ---- Pos<Zero> construction ----

impl Pos<Zero> {
    /// Create a 0-based position from a `u32`.
    ///
    /// Panics in debug builds if `value == u32::MAX` (reserved niche).
    /// BAM positions are capped at i32::MAX (~2.1B), well below the limit.
    #[inline]
    pub fn new(value: u32) -> Self {
        debug_assert!(value != u32::MAX, "position u32::MAX is reserved as niche");
        // Safety: debug_assert above ensures value != u32::MAX. In practice BAM
        // positions never approach this limit.
        Self { value: unsafe { NonMaxU32::new_unchecked(value) }, _system: PhantomData }
    }

    /// Create a 0-based position from an `i64`. Returns `None` if negative,
    /// > `u32::MAX - 1`, or exactly `u32::MAX`.
    #[inline]
    pub fn try_from_i64(value: i64) -> Option<Self> {
        let v = u32::try_from(value).ok()?;
        if v == u32::MAX {
            return None;
        }
        Some(Self::new(v))
    }

    /// Create a 0-based position from a `u64`. Returns `None` if > `u32::MAX - 1`
    /// or exactly `u32::MAX`.
    #[inline]
    pub fn try_from_u64(value: u64) -> Option<Self> {
        let v = u32::try_from(value).ok()?;
        if v == u32::MAX {
            return None;
        }
        Some(Self::new(v))
    }

    /// Convert to 1-based. Infallible: 0-based 0 → 1-based 1.
    ///
    /// BAM positions cap at i32::MAX (~2.1B), so +1 stays far below `u32::MAX - 1`.
    #[inline]
    pub fn to_one_based(self) -> Pos<One> {
        let new_val = self.value.get() + 1;
        debug_assert!(new_val != u32::MAX, "to_one_based would produce reserved niche value");
        // Safety: BAM positions are at most i32::MAX; +1 gives at most ~2.1B+1,
        // far below u32::MAX.
        Pos { value: unsafe { NonMaxU32::new_unchecked(new_val) }, _system: PhantomData }
    }
}

// ---- Pos<One> construction ----

impl Pos<One> {
    /// Create a 1-based position. Panics in debug mode if value is 0 or `u32::MAX`.
    #[inline]
    pub fn new(value: u32) -> Self {
        debug_assert!(value > 0, "1-based position must be >= 1");
        debug_assert!(value != u32::MAX, "position u32::MAX is reserved as niche");
        // Safety: debug_asserts above ensure value is in 1..u32::MAX.
        Self { value: unsafe { NonMaxU32::new_unchecked(value) }, _system: PhantomData }
    }

    /// Create a 1-based position from a `u32`. Returns `None` if value is 0 or `u32::MAX`.
    #[inline]
    pub fn try_new(value: u32) -> Option<Self> {
        if value == 0 || value == u32::MAX {
            return None;
        }
        Some(Self::new(value))
    }

    /// Create a 1-based position from an `i64`. Returns `None` if < 1, > `u32::MAX - 1`,
    /// or exactly `u32::MAX`.
    #[inline]
    pub fn try_from_i64(value: i64) -> Option<Self> {
        if value < 1 {
            return None;
        }
        u32::try_from(value).ok().and_then(Self::try_new)
    }

    /// Create a 1-based position from an `i32`. Returns `None` if < 1.
    #[inline]
    pub fn try_from_i32(value: i32) -> Option<Self> {
        if value < 1 {
            return None;
        }
        Self::try_new(value as u32)
    }

    /// Convert to 0-based. Infallible: 1-based 1 → 0-based 0.
    ///
    /// The result of subtracting 1 from a value in 1..u32::MAX is always in
    /// 0..(u32::MAX-1), which is always a valid NonMaxU32.
    #[inline]
    pub fn to_zero_based(self) -> Pos<Zero> {
        let new_val = self.value.get() - 1;
        // Safety: self.value >= 1 (enforced by construction), so new_val >= 0
        // and new_val <= u32::MAX - 2, which is always != u32::MAX.
        Pos { value: unsafe { NonMaxU32::new_unchecked(new_val) }, _system: PhantomData }
    }
}

// ---- Common methods (both systems) ----

impl<S> Pos<S> {
    /// Raw u32 value in the position's native coordinate system.
    #[inline]
    #[must_use]
    pub fn get(self) -> u32 {
        self.value.get()
    }

    /// Convenience for indexing: returns the raw value as usize.
    #[inline]
    #[must_use]
    pub fn as_usize(self) -> usize {
        self.value.get() as usize
    }

    /// Convenience for wider arithmetic: returns the raw value as i64.
    #[inline]
    #[must_use]
    pub fn as_i64(self) -> i64 {
        self.value.get() as i64
    }
}

// ---- Arithmetic ----

// Pos + Offset = Pos (same system)
impl<S> Add<Offset> for Pos<S> {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Offset) -> Self {
        let result = self.value.get() as i64 + rhs.0;
        debug_assert!(result >= 0 && result < u32::MAX as i64, "position overflow");
        // Safety: debug_assert ensures result is in 0..(u32::MAX-1).
        Pos { value: unsafe { NonMaxU32::new_unchecked(result as u32) }, _system: PhantomData }
    }
}

impl<S: Copy> AddAssign<Offset> for Pos<S> {
    #[inline]
    fn add_assign(&mut self, rhs: Offset) {
        *self = *self + rhs;
    }
}

// Pos - Offset = Pos (same system)
impl<S> Sub<Offset> for Pos<S> {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Offset) -> Self {
        self + Offset(-rhs.0)
    }
}

impl<S: Copy> SubAssign<Offset> for Pos<S> {
    #[inline]
    fn sub_assign(&mut self, rhs: Offset) {
        *self = *self - rhs;
    }
}

// Pos - Pos = Offset (same system only)
impl<S> Sub for Pos<S> {
    type Output = Offset;
    #[inline]
    fn sub(self, rhs: Self) -> Offset {
        Offset(self.value.get() as i64 - rhs.value.get() as i64)
    }
}

// Offset arithmetic
impl Add for Offset {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Offset(self.0 + rhs.0)
    }
}

impl Sub for Offset {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Offset(self.0 - rhs.0)
    }
}

impl Offset {
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
}

// ---- Display / Debug ----

impl<S> fmt::Debug for Pos<S> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Pos({})", self.value)
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
mod tests {
    use super::*;

    #[test]
    fn zero_based_roundtrip() {
        let z = Pos::<Zero>::new(100);
        let o = z.to_one_based();
        assert_eq!(o.get(), 101);
        assert_eq!(o.to_zero_based(), z);
    }

    #[test]
    fn one_based_roundtrip() {
        let o = Pos::<One>::new(1);
        let z = o.to_zero_based();
        assert_eq!(z.get(), 0);
        assert_eq!(z.to_one_based(), o);
    }

    #[test]
    fn pos_minus_pos_is_offset() {
        let a = Pos::<Zero>::new(100);
        let b = Pos::<Zero>::new(50);
        let off = a - b;
        assert_eq!(off.get(), 50);
        assert_eq!(b + off, a);
    }

    #[test]
    fn pos_minus_pos_negative_offset() {
        let a = Pos::<Zero>::new(10);
        let b = Pos::<Zero>::new(50);
        let off = a - b;
        assert_eq!(off.get(), -40);
    }

    #[test]
    fn pos_plus_offset() {
        let p = Pos::<Zero>::new(10);
        let q = p + Offset(5);
        assert_eq!(q.get(), 15);
    }

    #[test]
    fn pos_minus_offset() {
        let p = Pos::<Zero>::new(10);
        let q = p - Offset(3);
        assert_eq!(q.get(), 7);
    }

    #[test]
    fn try_from_i64_rejects_negative() {
        assert!(Pos::<Zero>::try_from_i64(-1).is_none());
        assert!(Pos::<One>::try_from_i64(0).is_none());
        assert!(Pos::<One>::try_from_i64(-1).is_none());
    }

    #[test]
    fn try_from_i64_rejects_u32_max() {
        assert!(Pos::<Zero>::try_from_i64(u32::MAX as i64).is_none());
        assert!(Pos::<One>::try_from_i64(u32::MAX as i64).is_none());
    }

    #[test]
    fn try_from_i64_accepts_valid() {
        assert_eq!(Pos::<Zero>::try_from_i64(0).unwrap().get(), 0);
        assert_eq!(Pos::<One>::try_from_i64(1).unwrap().get(), 1);
        assert_eq!(Pos::<Zero>::try_from_i64(100).unwrap().get(), 100);
    }

    #[test]
    fn try_from_u64_rejects_overflow() {
        assert!(Pos::<Zero>::try_from_u64(u64::from(u32::MAX) + 1).is_none());
    }

    #[test]
    fn try_from_u64_rejects_u32_max() {
        assert!(Pos::<Zero>::try_from_u64(u32::MAX as u64).is_none());
    }

    #[test]
    fn ordering() {
        let a = Pos::<Zero>::new(10);
        let b = Pos::<Zero>::new(20);
        assert!(a < b);
        assert!(b > a);
    }

    #[test]
    fn as_usize() {
        let p = Pos::<Zero>::new(42);
        assert_eq!(p.as_usize(), 42);
    }

    #[test]
    fn as_i64() {
        let p = Pos::<Zero>::new(100);
        assert_eq!(p.as_i64(), 100);
    }

    #[test]
    fn size_is_u32() {
        assert_eq!(std::mem::size_of::<Pos<Zero>>(), std::mem::size_of::<u32>());
        assert_eq!(std::mem::size_of::<Pos<One>>(), std::mem::size_of::<u32>());
    }

    #[test]
    fn option_pos_same_size_as_pos() {
        assert_eq!(std::mem::size_of::<Option<Pos<Zero>>>(), std::mem::size_of::<Pos<Zero>>(),);
        assert_eq!(std::mem::size_of::<Option<Pos<Zero>>>(), 4);
    }
}
