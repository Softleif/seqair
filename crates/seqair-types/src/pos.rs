//! Genomic position newtype parameterized by coordinate system.
//!
//! `Pos<Zero>` is 0-based (BAM, BED, internal engine). `Pos<One>` is 1-based
//! (SAM, VCF, CRAM, user-facing). The type system prevents mixing coordinate
//! systems at compile time. `Offset` represents a signed distance between positions.
//!
//! # Range and niche optimization
//!
//! Valid positions are `0..=i32::MAX` (for `Pos<Zero>`) or `1..=i32::MAX`
//! (for `Pos<One>`). This cap matches BAM/CRAM/BCF which store positions as
//! `i32`, and makes `as_i32()` infallible.

use std::fmt;
use std::marker::PhantomData;
use std::ops::{Deref, Sub};

/// Maximum valid raw value for any `Pos`: `i32::MAX` (2,147,483,647).
///
/// We cap at `i32::MAX` rather than `u32::MAX - 1` because BAM, CRAM, and BCF
/// all store positions as `i32`.  This lets `as_i32()` be infallible.
const POS_MAX: u32 = i32::MAX as u32;

// r[impl pos.systems]
/// 0-based coordinate system (BAM binary, BED, internal engine).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Zero;

// r[impl pos.systems]
/// 1-based coordinate system (SAM text, VCF, CRAM, user-facing).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct One;

/// 0-based position (BAM, BED). Alias for `Pos<Zero>`.
pub type Pos0 = Pos<Zero>;
/// 1-based position (SAM, VCF). Alias for `Pos<One>`.
pub type Pos1 = Pos<One>;

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
/// Valid range is `0..=i32::MAX` for `Pos<Zero>` and `1..=i32::MAX` for
/// `Pos<One>`. `u32::MAX` is reserved as the niche so `Option<Pos<S>>` is
/// 4 bytes.
///
/// # Construction
///
/// ```
/// use seqair_types::pos::{Pos0, Pos1};
///
/// let bam_pos = Pos0::new(100).unwrap();   // 0-based position 100
/// let sam_pos = Pos1::new(101).unwrap();   // 1-based position 101
/// assert_eq!(bam_pos, sam_pos.to_zero_based()); // same genomic location
///
/// // From i32 (BAM wire format):
/// let from_bam = Pos0::try_from(42i32).unwrap();
/// assert_eq!(from_bam.as_i32(), 42);
/// ```
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(transparent)]
pub struct Pos<S> {
    value: u32,
    _system: PhantomData<S>,
}

// r[impl pos.offset]
/// Signed distance between two positions.
///
/// `Pos - Pos = Offset` and `Pos + Offset = Pos`. You cannot add two positions.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(transparent)]
pub struct Offset(i64);

/// Error returned when a value cannot be converted to a [`Pos`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PosOverflow;

impl fmt::Display for PosOverflow {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("position value out of valid range")
    }
}

impl std::error::Error for PosOverflow {}

// ---- Pos<Zero> construction ----

impl Pos<Zero> {
    // r[impl pos.zero_new]
    /// Create a 0-based position from a `u32`. Returns `None` if `value > i32::MAX`.
    #[inline]
    pub const fn new(value: u32) -> Option<Self> {
        if value > POS_MAX {
            return None;
        }
        Some(Self { value, _system: PhantomData })
    }

    // r[impl pos.to_one_based]
    // r[impl pos.explicit_conversion]
    /// Convert to 1-based. Fails only at `i32::MAX` (0-based) where the
    /// 1-based result would exceed `i32::MAX`.
    #[inline]
    #[must_use]
    pub const fn to_one_based(self) -> Result<Pos<One>, PosOverflow> {
        let Some(new_val) = self.value.checked_add(1) else {
            // impossible by construction
            return Err(PosOverflow);
        };
        match Pos::<One>::new(new_val) {
            Some(v) => Ok(v),
            None => Err(PosOverflow),
        }
    }

    /// Maximum valid 0-based position (`i32::MAX`).
    /// Used as "end of contig" sentinel in queries.
    pub const fn max_value() -> Self {
        match Self::new(POS_MAX) {
            Some(v) => v,
            None => unreachable!(),
        }
    }
}

// ---- Pos<One> construction ----

impl Pos<One> {
    // r[impl pos.one_new]
    /// Create a 1-based position from a `u32`. Returns `None` if value is 0 or > `i32::MAX`.
    #[inline]
    pub const fn new(value: u32) -> Option<Self> {
        if value == 0 || value > POS_MAX {
            return None;
        }
        Some(Self { value, _system: PhantomData })
    }

    // r[impl pos.to_zero_based]
    // r[impl pos.explicit_conversion]
    /// Convert to 0-based. Infallible: 1-based values are in `1..=i32::MAX`,
    /// so subtracting 1 gives `0..=i32::MAX - 1`, always valid.
    #[inline]
    #[must_use]
    pub const fn to_zero_based(self) -> Pos<Zero> {
        let Some(new_val) = self.value.checked_sub(1) else {
            // always >0 by construction
            unreachable!()
        };
        Pos { value: new_val, _system: PhantomData }
    }

    /// Maximum valid 1-based position (`i32::MAX`).
    pub const fn max_value() -> Self {
        match Self::new(POS_MAX) {
            Some(v) => v,
            None => unreachable!(),
        }
    }
}

// ---- Common methods (both systems) ----

impl<S> Pos<S> {
    // r[impl pos.as_i32]
    // r[impl pos.must_use]
    /// Raw value as `i32`. Infallible because all valid positions are `<= i32::MAX`.
    #[inline]
    #[must_use]
    #[expect(clippy::cast_possible_wrap, reason = "value ≤ i32::MAX by construction")]
    pub const fn as_i32(self) -> i32 {
        self.value as i32
    }

    // r[impl pos.as_usize]
    /// Convenience for indexing: returns the raw value as usize.
    #[inline]
    #[must_use]
    pub const fn as_usize(self) -> usize {
        self.value as usize
    }

    // r[impl pos.as_i64]
    /// Convenience for wider arithmetic: returns the raw value as i64.
    #[inline]
    #[must_use]
    pub const fn as_i64(self) -> i64 {
        self.value as i64
    }

    /// Convenience for 64-bit arithmetic: returns the raw value as u64.
    #[inline]
    #[must_use]
    pub const fn as_u64(self) -> u64 {
        self.value as u64
    }

    // r[impl pos.add_offset]
    /// Checked position + offset. Returns `None` if result is negative or > `i32::MAX`.
    #[inline]
    #[must_use]
    #[expect(
        clippy::cast_sign_loss,
        clippy::cast_possible_truncation,
        reason = "result is checked to be in 0..=i32::MAX"
    )]
    pub fn checked_add_offset(self, offset: Offset) -> Option<Self> {
        let result = i64::from(self.value).wrapping_add(offset.0);
        if result < 0 || result > i64::from(POS_MAX) {
            return None;
        }
        // result is in 0..=i32::MAX
        Some(Self { value: result as u32, _system: PhantomData })
    }

    // r[impl pos.sub_offset]
    /// Checked position - offset. Returns `None` if result is negative or > `i32::MAX`.
    #[inline]
    #[must_use]
    pub fn checked_sub_offset(self, offset: Offset) -> Option<Self> {
        let negated = Offset(offset.0.checked_neg()?);
        self.checked_add_offset(negated)
    }
}

impl<S> Deref for Pos<S> {
    type Target = u32;
    #[inline]
    fn deref(&self) -> &u32 {
        &self.value
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
        Offset(i64::from(self.value).wrapping_sub(i64::from(rhs.value)))
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
    #[expect(
        clippy::cast_possible_truncation,
        reason = "Offset wraps i64; unsigned_abs() ≤ i64::MAX which fits in usize on 64-bit; on 32-bit platforms the calling code guarantees values bounded by BAM/FASTA region sizes"
    )]
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
    pub const fn checked_add(self, rhs: Self) -> Option<Self> {
        match self.0.checked_add(rhs.0) {
            Some(v) => Some(Offset(v)),
            None => None,
        }
    }

    /// Checked subtraction. Returns `None` on overflow.
    #[inline]
    pub const fn checked_sub(self, rhs: Self) -> Option<Self> {
        match self.0.checked_sub(rhs.0) {
            Some(v) => Some(Offset(v)),
            None => None,
        }
    }
}

// ---- From / TryFrom impls ----

impl<S> From<Pos<S>> for u32 {
    #[inline]
    fn from(pos: Pos<S>) -> u32 {
        pos.value
    }
}

impl<S> From<Pos<S>> for i32 {
    #[inline]
    #[expect(clippy::cast_possible_wrap, reason = "value ≤ i32::MAX by construction")]
    fn from(pos: Pos<S>) -> i32 {
        pos.value as i32
    }
}

impl<S> From<Pos<S>> for u64 {
    #[inline]
    fn from(pos: Pos<S>) -> u64 {
        u64::from(pos.value)
    }
}

impl<S> From<Pos<S>> for i64 {
    #[inline]
    fn from(pos: Pos<S>) -> i64 {
        i64::from(pos.value)
    }
}

impl<S> From<Pos<S>> for usize {
    #[inline]
    fn from(pos: Pos<S>) -> usize {
        pos.value as usize
    }
}

// r[impl pos.try_from]
impl TryFrom<i32> for Pos<Zero> {
    type Error = PosOverflow;
    /// Create a 0-based position from an `i32`. Fails if negative.
    #[inline]
    fn try_from(value: i32) -> Result<Self, PosOverflow> {
        if value < 0 {
            return Err(PosOverflow);
        }
        // 0..=i32::MAX always valid.
        Ok(Self { value: value as u32, _system: PhantomData })
    }
}

// r[impl pos.try_from]
impl TryFrom<i32> for Pos<One> {
    type Error = PosOverflow;
    /// Create a 1-based position from an `i32`. Fails if < 1.
    #[inline]
    fn try_from(value: i32) -> Result<Self, PosOverflow> {
        if value < 1 {
            return Err(PosOverflow);
        }
        Ok(Self { value: value as u32, _system: PhantomData })
    }
}

// r[impl pos.try_from]
impl TryFrom<u32> for Pos<Zero> {
    type Error = PosOverflow;
    /// Create a 0-based position from a `u32`. Fails if > `i32::MAX`.
    #[inline]
    fn try_from(value: u32) -> Result<Self, PosOverflow> {
        Self::new(value).ok_or(PosOverflow)
    }
}

// r[impl pos.try_from]
impl TryFrom<u32> for Pos<One> {
    type Error = PosOverflow;
    /// Create a 1-based position from a `u32`. Fails if 0 or > `i32::MAX`.
    #[inline]
    fn try_from(value: u32) -> Result<Self, PosOverflow> {
        Self::new(value).ok_or(PosOverflow)
    }
}

// r[impl pos.try_from]
impl TryFrom<i64> for Pos<Zero> {
    type Error = PosOverflow;
    /// Create a 0-based position from an `i64`. Fails if negative or > `i32::MAX`.
    #[inline]
    fn try_from(value: i64) -> Result<Self, PosOverflow> {
        let v = i32::try_from(value).map_err(|_| PosOverflow)?;
        Self::try_from(v)
    }
}

// r[impl pos.try_from]
impl TryFrom<i64> for Pos<One> {
    type Error = PosOverflow;
    /// Create a 1-based position from an `i64`. Fails if < 1 or > `i32::MAX`.
    #[inline]
    fn try_from(value: i64) -> Result<Self, PosOverflow> {
        let v = i32::try_from(value).map_err(|_| PosOverflow)?;
        Self::try_from(v)
    }
}

// r[impl pos.try_from]
impl TryFrom<u64> for Pos<Zero> {
    type Error = PosOverflow;
    /// Create a 0-based position from a `u64`. Fails if > `i32::MAX`.
    #[inline]
    fn try_from(value: u64) -> Result<Self, PosOverflow> {
        let v = i32::try_from(value).map_err(|_| PosOverflow)?;
        Self::try_from(v)
    }
}

// r[impl pos.try_from]
impl TryFrom<u64> for Pos<One> {
    type Error = PosOverflow;
    /// Create a 1-based position from a `u64`. Fails if 0 or > `i32::MAX`.
    #[inline]
    fn try_from(value: u64) -> Result<Self, PosOverflow> {
        let v = i32::try_from(value).map_err(|_| PosOverflow)?;
        Self::try_from(v)
    }
}

// ---- Display / Debug ----

impl fmt::Debug for Pos<Zero> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Pos0({})", self.value)
    }
}

impl fmt::Debug for Pos<One> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Pos1({})", self.value)
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
    /// Test-only convenience: panics if value is out of range.
    pub fn at(value: u32) -> Self {
        Self::new(value).expect("test position out of range")
    }
}

#[cfg(test)]
impl Pos<One> {
    /// Test-only convenience: panics if value is out of range.
    pub fn at(value: u32) -> Self {
        Self::new(value).expect("test position out of range")
    }
}

#[cfg(test)]
#[allow(clippy::arithmetic_side_effects, reason = "tests")]
mod tests {
    use proptest::prelude::*;

    use super::*;

    const I32_MAX_U32: u32 = i32::MAX as u32;

    // r[verify pos.to_one_based]
    // r[verify pos.to_zero_based]
    #[test]
    fn zero_based_roundtrip() {
        let z = Pos0::new(100).unwrap();
        let o = z.to_one_based().unwrap();
        assert_eq!(*o, 101);
        assert_eq!(o.to_zero_based(), z);
    }

    // r[verify pos.to_zero_based]
    // r[verify pos.to_one_based]
    #[test]
    fn one_based_roundtrip() {
        let o = Pos1::new(1).unwrap();
        let z = o.to_zero_based();
        assert_eq!(*z, 0);
        assert_eq!(z.to_one_based().unwrap(), o);
    }

    // r[verify pos.sub_pos]
    // r[verify pos.add_offset]
    #[test]
    fn pos_minus_pos_is_offset() {
        let a = Pos0::new(100).unwrap();
        let b = Pos0::new(50).unwrap();
        let off = a - b;
        assert_eq!(off.get(), 50);
        assert_eq!(b.checked_add_offset(off).unwrap(), a);
    }

    // r[verify pos.sub_pos]
    #[test]
    fn pos_minus_pos_negative_offset() {
        let a = Pos0::new(10).unwrap();
        let b = Pos0::new(50).unwrap();
        let off = a - b;
        assert_eq!(off.get(), -40);
    }

    // r[verify pos.add_offset]
    #[test]
    fn pos_plus_offset() {
        let p = Pos0::new(10).unwrap();
        let q = p.checked_add_offset(Offset::new(5)).unwrap();
        assert_eq!(*q, 15);
    }

    // r[verify pos.sub_offset]
    #[test]
    fn pos_minus_offset() {
        let p = Pos0::new(10).unwrap();
        let q = p.checked_sub_offset(Offset::new(3)).unwrap();
        assert_eq!(*q, 7);
    }

    // r[verify pos.try_from]
    #[test]
    fn try_from_i32_rejects_negative() {
        assert!(Pos0::try_from(-1i32).is_err());
        assert!(Pos1::try_from(0i32).is_err());
        assert!(Pos1::try_from(-1i32).is_err());
    }

    // r[verify pos.try_from]
    #[test]
    fn try_from_i32_accepts_valid() {
        assert_eq!(*Pos0::try_from(0i32).unwrap(), 0);
        assert_eq!(*Pos1::try_from(1i32).unwrap(), 1);
        assert_eq!(*Pos0::try_from(100i32).unwrap(), 100);
    }

    // r[verify pos.try_from]
    #[test]
    fn try_from_i64_rejects_negative() {
        assert!(Pos0::try_from(-1i64).is_err());
        assert!(Pos1::try_from(0i64).is_err());
        assert!(Pos1::try_from(-1i64).is_err());
    }

    // r[verify pos.try_from]
    #[test]
    fn try_from_i64_rejects_overflow() {
        assert!(Pos0::try_from(i64::from(i32::MAX) + 1).is_err());
        assert!(Pos1::try_from(i64::from(i32::MAX) + 1).is_err());
        assert!(Pos0::try_from(i64::from(u32::MAX)).is_err());
    }

    // r[verify pos.try_from]
    #[test]
    fn try_from_i64_accepts_valid() {
        assert_eq!(*Pos0::try_from(0i64).unwrap(), 0);
        assert_eq!(*Pos1::try_from(1i64).unwrap(), 1);
        assert_eq!(*Pos0::try_from(100i64).unwrap(), 100);
    }

    // r[verify pos.try_from]
    #[test]
    fn try_from_u64_rejects_overflow() {
        assert!(Pos0::try_from(u64::from(u32::MAX) + 1).is_err());
        assert!(Pos0::try_from(u64::from(u32::MAX)).is_err());
        assert!(Pos0::try_from(i32::MAX as u64 + 1).is_err());
        assert!(Pos1::try_from(0u64).is_err(), "Pos1 rejects 0");
        assert!(Pos1::try_from(i32::MAX as u64 + 1).is_err());
        assert_eq!(*Pos1::try_from(1u64).unwrap(), 1);
    }

    // r[verify pos.try_from]
    #[test]
    fn try_from_u32() {
        assert_eq!(*Pos0::try_from(0u32).unwrap(), 0);
        assert_eq!(*Pos0::try_from(i32::MAX as u32).unwrap(), i32::MAX as u32);
        assert!(Pos0::try_from(i32::MAX as u32 + 1).is_err());
        assert!(Pos1::try_from(0u32).is_err());
        assert_eq!(*Pos1::try_from(1u32).unwrap(), 1);
        assert!(Pos1::try_from(i32::MAX as u32 + 1).is_err());
    }

    #[test]
    fn as_i32_roundtrips() {
        let p = Pos0::new(42).unwrap();
        assert_eq!(p.as_i32(), 42);
        assert_eq!(Pos0::try_from(p.as_i32()).unwrap(), p);

        let max = Pos0::max_value();
        assert_eq!(max.as_i32(), i32::MAX);
    }

    // r[verify pos.derives]
    #[test]
    fn ordering() {
        let a = Pos0::new(10).unwrap();
        let b = Pos0::new(20).unwrap();
        assert!(a < b);
        assert!(b > a);
    }

    // r[verify pos.as_usize]
    #[test]
    fn as_usize() {
        let p = Pos0::new(42).unwrap();
        assert_eq!(p.as_usize(), 42);
    }

    // r[verify pos.as_i64]
    #[test]
    fn as_i64() {
        let p = Pos0::new(100).unwrap();
        assert_eq!(p.as_i64(), 100);
    }

    // r[verify pos.size]
    #[test]
    fn size_is_u32() {
        assert_eq!(std::mem::size_of::<Pos0>(), std::mem::size_of::<u32>());
        assert_eq!(std::mem::size_of::<Pos1>(), std::mem::size_of::<u32>());
    }

    // r[verify pos.derives]
    #[test]
    fn debug_shows_coordinate_system() {
        let z = Pos0::new(42).unwrap();
        let o = Pos1::new(43).unwrap();
        assert_eq!(format!("{z:?}"), "Pos0(42)");
        assert_eq!(format!("{o:?}"), "Pos1(43)");
    }

    // r[verify pos.zero_new]
    #[test]
    fn max_value_zero() {
        let m = Pos0::max_value();
        assert_eq!(*m, I32_MAX_U32);
    }

    // r[verify pos.one_new]
    #[test]
    fn max_value_one() {
        let m = Pos1::max_value();
        assert_eq!(*m, I32_MAX_U32);
    }

    #[test]
    fn new_rejects_above_i32_max() {
        assert!(Pos0::new(I32_MAX_U32 + 1).is_none());
        assert!(Pos1::new(I32_MAX_U32 + 1).is_none());
        assert!(Pos0::new(u32::MAX).is_none());
        assert!(Pos1::new(u32::MAX).is_none());
    }

    #[test]
    fn new_accepts_i32_max() {
        assert!(Pos0::new(I32_MAX_U32).is_some());
        assert!(Pos1::new(I32_MAX_U32).is_some());
    }

    // r[verify pos.one_new]
    #[test]
    fn one_based_new_rejects_zero() {
        assert!(Pos1::new(0).is_none(), "0 should be rejected for 1-based");
    }

    // r[verify pos.to_one_based]
    #[test]
    fn to_one_based_returns_err_at_max() {
        let p = Pos0::max_value();
        assert!(p.to_one_based().is_err(), "i32::MAX + 1 would exceed i32::MAX");
    }

    // r[verify pos.add_offset]
    #[test]
    fn checked_add_offset_rejects_overflow() {
        let p = Pos0::new(I32_MAX_U32 - 1).unwrap();
        assert!(p.checked_add_offset(Offset::new(1)).is_some());
        assert!(p.checked_add_offset(Offset::new(2)).is_none(), "result would exceed i32::MAX");
    }

    // r[verify pos.add_offset]
    #[test]
    fn checked_add_offset_rejects_negative_result() {
        let p = Pos0::new(5).unwrap();
        assert!(p.checked_add_offset(Offset::new(-10)).is_none(), "result would be negative");
    }

    // r[verify pos.sub_offset]
    #[test]
    fn checked_sub_offset_works() {
        let p = Pos0::new(10).unwrap();
        assert_eq!(*p.checked_sub_offset(Offset::new(3)).unwrap(), 7);
        assert!(p.checked_sub_offset(Offset::new(11)).is_none());
    }

    // r[verify pos.try_from]
    #[test]
    fn from_i32_infallible_roundtrip() {
        let p = Pos0::try_from(12345i32).unwrap();
        let back: i32 = p.into();
        assert_eq!(back, 12345);
    }

    proptest! {
        // r[verify pos.to_one_based]
        // r[verify pos.to_zero_based]
        #[test]
        fn roundtrip_zero_to_one_and_back(v in 0u32..I32_MAX_U32) {
            // i32::MAX itself can't round-trip (to_one_based would overflow).
            let z = Pos0::new(v).unwrap();
            let o = z.to_one_based().unwrap();
            prop_assert_eq!(o.to_zero_based(), z);
        }

        // r[verify pos.to_zero_based]
        // r[verify pos.to_one_based]
        #[test]
        fn roundtrip_one_to_zero_and_back(v in 1u32..=I32_MAX_U32) {
            let o = Pos1::new(v).unwrap();
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
            let result = i64::from(v) + off;
            if result >= 0 && result <= i64::from(I32_MAX_U32) {
                let p = Pos0::new(v).unwrap();
                let q = p.checked_add_offset(Offset::new(off)).unwrap();
                let r = q.checked_sub_offset(Offset::new(off)).unwrap();
                prop_assert_eq!(r, p);
            }
        }

        // r[verify pos.sub_pos]
        #[test]
        fn pos_sub_pos_is_offset(a in 0u32..1_000_000, b in 0u32..1_000_000) {
            let pa = Pos0::new(a).unwrap();
            let pb = Pos0::new(b).unwrap();
            let off = pa - pb;
            prop_assert_eq!(off.get(), i64::from(a) - i64::from(b));
        }

        #[test]
        fn new_never_accepts_above_i32_max(v in (I32_MAX_U32 + 1)..=u32::MAX) {
            prop_assert!(Pos0::new(v).is_none());
            prop_assert!(Pos1::new(v).is_none());
        }

        // r[verify pos.add_offset]
        #[test]
        fn checked_add_never_exceeds_i32_max(
            v in 0u32..=I32_MAX_U32,
            off in 0i64..=100,
        ) {
            if let Some(p) = Pos0::new(v)
                && let Some(result) = p.checked_add_offset(Offset::new(off))
            {
                prop_assert!(*result <= I32_MAX_U32, "checked_add must not exceed i32::MAX");
            }
        }

        // r[verify pos.as_i32]
        #[test]
        fn as_i32_always_valid(v in 0u32..=I32_MAX_U32) {
            let p = Pos0::new(v).unwrap();
            let i = p.as_i32();
            prop_assert!(i >= 0);
            prop_assert_eq!(i as u32, v);
        }
    }
}
