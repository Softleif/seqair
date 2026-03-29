//! Genomic position newtype parameterized by coordinate system.
//!
//! `Pos<Zero>` is 0-based (BAM, BED, internal engine). `Pos<One>` is 1-based
//! (SAM, VCF, CRAM, user-facing). The type system prevents mixing coordinate
//! systems at compile time. `Offset` represents a signed distance between positions.

use std::fmt;
use std::marker::PhantomData;
use std::ops::{Add, AddAssign, Sub, SubAssign};

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
    value: u32,
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
    /// Create a 0-based position from a `u32`. All u32 values are valid.
    #[inline]
    pub const fn new(value: u32) -> Self {
        Self { value, _system: PhantomData }
    }

    /// Create a 0-based position from an `i64`. Returns `None` if negative or > u32::MAX.
    #[inline]
    pub fn try_from_i64(value: i64) -> Option<Self> {
        u32::try_from(value).ok().map(Self::new)
    }

    /// Create a 0-based position from a `u64`. Returns `None` if > u32::MAX.
    #[inline]
    pub fn try_from_u64(value: u64) -> Option<Self> {
        u32::try_from(value).ok().map(Self::new)
    }

    /// Convert to 1-based. Infallible: 0-based 0 → 1-based 1.
    ///
    /// This is safe because the maximum 0-based BAM position (i32::MAX ≈ 2.1B)
    /// plus 1 still fits in u32 (max ≈ 4.3B).
    #[inline]
    pub const fn to_one_based(self) -> Pos<One> {
        Pos { value: self.value + 1, _system: PhantomData }
    }
}

// ---- Pos<One> construction ----

impl Pos<One> {
    /// Create a 1-based position. Panics in debug mode if value is 0.
    #[inline]
    pub const fn new(value: u32) -> Self {
        debug_assert!(value > 0, "1-based position must be >= 1");
        Self { value, _system: PhantomData }
    }

    /// Create a 1-based position from a `u32`. Returns `None` if value is 0.
    #[inline]
    pub fn try_new(value: u32) -> Option<Self> {
        if value > 0 { Some(Self { value, _system: PhantomData }) } else { None }
    }

    /// Create a 1-based position from an `i64`. Returns `None` if < 1 or > u32::MAX.
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
        Some(Self::new(value as u32))
    }

    /// Convert to 0-based. Infallible: 1-based 1 → 0-based 0.
    #[inline]
    pub const fn to_zero_based(self) -> Pos<Zero> {
        Pos { value: self.value - 1, _system: PhantomData }
    }
}

// ---- Common methods (both systems) ----

impl<S> Pos<S> {
    /// Raw u32 value in the position's native coordinate system.
    #[inline]
    #[must_use]
    pub const fn get(self) -> u32 {
        self.value
    }

    /// Convenience for indexing: returns the raw value as usize.
    #[inline]
    #[must_use]
    pub const fn as_usize(self) -> usize {
        self.value as usize
    }

    /// Convenience for wider arithmetic: returns the raw value as i64.
    #[inline]
    #[must_use]
    pub const fn as_i64(self) -> i64 {
        self.value as i64
    }
}

// ---- Arithmetic ----

// Pos + Offset = Pos (same system)
impl<S> Add<Offset> for Pos<S> {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Offset) -> Self {
        let result = self.value as i64 + rhs.0;
        debug_assert!(result >= 0 && result <= u32::MAX as i64, "position overflow");
        Pos { value: result as u32, _system: PhantomData }
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
        Offset(self.value as i64 - rhs.value as i64)
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
}
