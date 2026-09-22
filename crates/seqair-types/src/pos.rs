//! Genomic positions, in the two coordinate systems the formats use.
//!
//! [`Pos0`] is 0-based (BAM, BED, internal engine). [`Pos1`] is 1-based
//! (SAM, VCF, CRAM, user-facing). They are two separate types carrying the
//! same API rather than one type parameterized by a marker: mixing them is a
//! compile error, and neither can inherit the other's floor. [`Offset`]
//! represents a signed distance between positions.
//!
//! # Range and niche optimization
//!
//! Valid positions are `0..=i32::MAX` (for [`Pos0`]) or `1..=i32::MAX`
//! (for [`Pos1`]). This cap matches BAM/CRAM/BCF which store positions as
//! `i32`, and makes `as_i32()` infallible. Both types store a
//! [`NonZeroU32`], so `Option<Pos0>` and `Option<Pos1>` are four bytes.

use std::fmt;
use std::num::NonZeroU32;
use std::ops::Sub;

/// Maximum valid raw value for either position type: `i32::MAX`
/// (2,147,483,647).
///
/// We cap at `i32::MAX` rather than `u32::MAX - 1` because BAM, CRAM, and BCF
/// all store positions as `i32`.  This lets `as_i32()` be infallible.
const POS_MAX: u32 = i32::MAX as u32;

/// `POS_MAX` as the stored `NonZeroU32`: the representation of [`Pos1::MAX`],
/// and one less than the representation of [`Pos0::MAX`].
const POS_MAX_NZ: NonZeroU32 = NonZeroU32::new(POS_MAX).expect("i32::MAX is not zero");

/// The serde "expected" text depends on `POS_MAX` being spelled out, since a
/// `&'static str` cannot interpolate a constant.
const _: () = assert!(POS_MAX == 2_147_483_647);

// r[impl pos.type]
// r[impl pos.size]
// r[impl pos.systems]
// r[impl pos.incompatible]
// r[impl pos.derives]
// r[impl pos.niche]
/// A 0-based position: BAM binary, BED, and everything internal to the engine.
///
/// Valid values are `0..=i32::MAX`.
///
/// # Representation
///
/// The field holds the position **plus one**, so that the niche of
/// [`NonZeroU32`] makes `Option<Pos0>` four bytes rather than eight. Storing
/// `value + 1` rather than `value` is what lets a `Pos0` and the [`Pos1`] of
/// the same genomic base have the identical bit pattern: [`Pos1`] stores its
/// value as it is, and the 1-based value of a base *is* its 0-based value plus
/// one. [`Self::to_one_based`] is therefore a reinterpretation plus the
/// `i32::MAX` bound check, and [`Pos1::to_zero_based`] a bare reinterpretation.
/// Adding one is monotone, so the derived `Ord` on the stored number is the
/// ordering of the positions.
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
///
/// // Each type knows its own floor: there is no 1-based zero.
/// assert_eq!(Pos1::new(0), None);
/// assert_eq!(Pos1::MIN, Pos1::new(1).unwrap());
/// assert_eq!(Pos0::MIN, Pos0::ZERO);
/// ```
///
/// # What the compiler must reject
///
/// These are the guarantees the two separate types exist for, and the only way
/// to test a compile error is to fail to compile. Note what that does and does
/// not buy: rustdoc checks that the block *fails*, but it does not check the
/// error code — `compile_fail,E0999` passes just as happily. So a block goes
/// vacuous the moment its example stops compiling for some unrelated reason,
/// a renamed constructor being the obvious one.
///
/// The twin below is the guard against that: it uses every name these blocks
/// use, and it has to compile. If `Pos0::new` or `QPos::new` is renamed, it
/// fails rather than the rejections silently becoming free.
///
/// ```
/// use seqair_types::pos::{Pos0, Pos1, QPos};
/// let zero = Pos0::new(100).unwrap();
/// let one = Pos1::new(100).unwrap();
/// let offset = QPos::new(5);
/// fn takes_a_position(_: Pos0) {}
/// fn takes_an_offset(_: QPos) {}
/// takes_a_position(one.to_zero_based());
/// takes_an_offset(offset);
/// assert_eq!(zero.as_u32(), 100);
/// ```
///
/// A 0-based and a 1-based position are different types
/// (r[`pos.incompatible`]):
///
/// ```compile_fail,E0308
/// use seqair_types::pos::{Pos0, Pos1};
/// let zero = Pos0::new(100).unwrap();
/// let one: Pos1 = zero;
/// ```
///
/// and they do not compare (r[`pos.incompatible`]):
///
/// ```compile_fail,E0308
/// use seqair_types::pos::{Pos0, Pos1};
/// let _ = Pos0::new(100).unwrap() == Pos1::new(100).unwrap();
/// ```
///
/// Conversion is only ever explicit (r[`pos.explicit_conversion`]):
///
/// ```compile_fail,E0277
/// use seqair_types::pos::{Pos0, Pos1};
/// let one: Pos1 = Pos0::new(100).unwrap().into();
/// ```
///
/// Two positions do not add — the sum of two locations is not a location
/// (r[`pos.no_add_pos`]):
///
/// ```compile_fail,E0369
/// use seqair_types::pos::Pos0;
/// let _ = Pos0::new(100).unwrap() + Pos0::new(1).unwrap();
/// ```
///
/// and a position is not an integer (r[`pos.explicit_conversion`]):
///
/// ```compile_fail,E0308
/// use seqair_types::pos::Pos0;
/// let n: u32 = Pos0::new(100).unwrap();
/// ```
///
/// A query offset indexes a read, not the reference, so it is not a position
/// either way round (r[`qpos.not_a_pos`]):
///
/// ```compile_fail,E0308
/// use seqair_types::pos::{Pos0, QPos};
/// fn takes_a_position(_: Pos0) {}
/// takes_a_position(QPos::new(5));
/// ```
///
/// ```compile_fail,E0308
/// use seqair_types::pos::{Pos0, QPos};
/// fn takes_an_offset(_: QPos) {}
/// takes_an_offset(Pos0::new(5).unwrap());
/// ```
// r[verify pos.incompatible]
// r[verify pos.explicit_conversion]
// r[verify pos.no_add_pos]
// r[verify qpos.not_a_pos]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(transparent)]
pub struct Pos0(NonZeroU32);

// r[impl pos.type]
// r[impl pos.size]
// r[impl pos.systems]
// r[impl pos.incompatible]
// r[impl pos.derives]
// r[impl pos.niche]
/// A 1-based position: SAM text, VCF, CRAM, and everything user-facing.
///
/// Valid values are `1..=i32::MAX`. The field holds the value as it is; the
/// [`NonZeroU32`] is both the niche that makes `Option<Pos1>` four bytes and
/// the reason no runtime check can ever produce the 1-based zero. See
/// [`Pos0`]'s representation note for why the two types share a bit pattern.
///
/// ```
/// use seqair_types::pos::Pos1;
///
/// assert_eq!(Pos1::new(0), None);
/// assert_eq!(Pos1::MIN.as_u32(), 1);
/// assert_eq!(Pos1::new(1).unwrap().to_zero_based().as_u32(), 0);
/// ```
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(transparent)]
pub struct Pos1(NonZeroU32);

// r[impl pos.offset]
/// Signed distance between two positions.
///
/// `Pos - Pos = Offset` and `Pos + Offset = Pos`. You cannot add two positions.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(transparent)]
pub struct Offset(i64);

/// Error returned when a value cannot be converted to a [`Pos0`] or [`Pos1`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PosOverflow;

impl fmt::Display for PosOverflow {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("position value out of valid range")
    }
}

impl std::error::Error for PosOverflow {}

// ---- Pos0: construction, accessors, arithmetic ----

impl Pos0 {
    /// Zero position (0-based). Valid and commonly used, so provided as a constant.
    pub const ZERO: Self = Self::MIN;

    // r[impl pos.min_max]
    // r[impl pos.two_types]
    /// The first position this type admits, `Pos0(0)`.
    pub const MIN: Self = Self(NonZeroU32::MIN);

    // r[impl pos.min_max]
    // r[impl pos.two_types]
    /// The last representable position, `i32::MAX`.
    /// Used as the "end of contig" sentinel in queries.
    pub const MAX: Self = Self(NonZeroU32::MIN.saturating_add(POS_MAX));

    // r[impl pos.zero_new]
    // r[impl pos.two_types]
    /// Create a 0-based position from a `u32`, rejecting anything above
    /// `i32::MAX`. Every value in `0..=i32::MAX` is a 0-based position.
    #[inline]
    pub const fn new(value: u32) -> Option<Self> {
        if value > POS_MAX {
            return None;
        }
        // The stored number is `value + 1`, built as "one, plus the value" so
        // it is a `NonZeroU32` by construction. Nothing saturates: `value` is
        // at most `i32::MAX`.
        Some(Self(NonZeroU32::MIN.saturating_add(value)))
    }

    // r[impl pos.as_u32]
    // r[impl pos.must_use]
    /// The 0-based value.
    #[inline]
    #[must_use]
    pub const fn as_u32(self) -> u32 {
        // The stored number is the value plus one and therefore at least one;
        // `wrapping_sub` cannot wrap. `-` would trip `arithmetic_side_effects`.
        self.0.get().wrapping_sub(1)
    }

    // r[impl pos.as_i32]
    // r[impl pos.must_use]
    /// Raw value as `i32`. Infallible because all valid positions are `<= i32::MAX`.
    #[inline]
    #[must_use]
    #[expect(clippy::cast_possible_wrap, reason = "value ≤ i32::MAX by construction")]
    pub const fn as_i32(self) -> i32 {
        self.as_u32() as i32
    }

    // r[impl pos.as_usize]
    // r[impl pos.must_use]
    /// Convenience for indexing: returns the raw value as usize.
    #[inline]
    #[must_use]
    pub const fn as_usize(self) -> usize {
        self.as_u32() as usize
    }

    // r[impl pos.as_i64]
    // r[impl pos.must_use]
    /// Convenience for wider arithmetic: returns the raw value as i64.
    #[inline]
    #[must_use]
    pub const fn as_i64(self) -> i64 {
        self.as_u32() as i64
    }

    // r[impl pos.must_use]
    /// Convenience for 64-bit arithmetic: returns the raw value as u64.
    #[inline]
    #[must_use]
    pub const fn as_u64(self) -> u64 {
        self.as_u32() as u64
    }

    // r[impl pos.to_one_based]
    // r[impl pos.explicit_conversion]
    /// Convert to 1-based. Fails only at `i32::MAX` (0-based) where the
    /// 1-based result would exceed `i32::MAX`.
    ///
    /// A reinterpretation: a `Pos0` already stores `value + 1`, which is the
    /// 1-based value. Only the bound differs between the two types.
    #[inline]
    pub const fn to_one_based(self) -> Result<Pos1, PosOverflow> {
        if self.0.get() > POS_MAX {
            return Err(PosOverflow);
        }
        Ok(Pos1(self.0))
    }

    // r[impl pos.add_offset]
    // r[impl pos.two_types]
    /// Checked `position + offset`. `None` if the result leaves this type's
    /// range: below [`Pos0::MIN`] — so below zero — or above `i32::MAX`.
    #[inline]
    #[must_use]
    pub fn checked_add_offset(self, offset: Offset) -> Option<Self> {
        let result = i64::from(self.as_u32()).checked_add(offset.0)?;
        u32::try_from(result).ok().and_then(Self::new)
    }

    // r[impl pos.sub_offset]
    // r[impl pos.two_types]
    /// Checked `position - offset`. Same bounds as [`Self::checked_add_offset`].
    ///
    /// `Offset::new(i64::MIN)` has no negation, so subtracting it is `None`
    /// whatever the position.
    #[inline]
    #[must_use]
    pub fn checked_sub_offset(self, offset: Offset) -> Option<Self> {
        self.checked_add_offset(Offset(offset.0.checked_neg()?))
    }

    // r[impl pos.saturating_offset]
    // r[impl pos.two_types]
    /// `position + offset`, clamped to this type's range instead of failing:
    /// a result below [`Pos0::MIN`] becomes `MIN`, one above `i32::MAX`
    /// becomes [`Pos0::MAX`].
    ///
    /// For padding and trimming, where running off the start of a contig means
    /// the start of the contig. A caller that needs to *know* the bound was hit
    /// wants [`Self::checked_add_offset`].
    #[inline]
    #[must_use]
    #[expect(
        clippy::cast_sign_loss,
        clippy::cast_possible_truncation,
        reason = "clamped into 0..=POS_MAX before the cast"
    )]
    pub fn saturating_add_offset(self, offset: Offset) -> Self {
        let result = i64::from(self.as_u32()).saturating_add(offset.0);
        let clamped = result.clamp(0, i64::from(POS_MAX));
        // Same construction as `new`, minus the check `clamp` made redundant.
        Self(NonZeroU32::MIN.saturating_add(clamped as u32))
    }

    // r[impl pos.saturating_offset]
    // r[impl pos.two_types]
    /// `position - offset`, clamped to this type's range instead of failing.
    ///
    /// `Offset::new(i64::MIN)` has no negation; subtracting it saturates
    /// to [`Pos0::MAX`], since the most negative offset really does run off
    /// the far end.
    #[inline]
    #[must_use]
    pub fn saturating_sub_offset(self, offset: Offset) -> Self {
        self.saturating_add_offset(Offset(offset.0.saturating_neg()))
    }
}

// ---- Pos1: construction, accessors, arithmetic ----

impl Pos1 {
    // r[impl pos.min_max]
    // r[impl pos.two_types]
    /// The first position this type admits, `Pos1(1)`.
    pub const MIN: Self = Self(NonZeroU32::MIN);

    // r[impl pos.min_max]
    // r[impl pos.two_types]
    /// The last representable position, `i32::MAX`.
    /// Used as the "end of contig" sentinel in queries.
    pub const MAX: Self = Self(POS_MAX_NZ);

    // r[impl pos.one_new]
    // r[impl pos.two_types]
    /// Create a 1-based position from a `u32`, rejecting `0` and anything
    /// above `i32::MAX`. The zero is not a runtime check that could be
    /// forgotten: the representation cannot hold it.
    #[inline]
    pub const fn new(value: u32) -> Option<Self> {
        if value > POS_MAX {
            return None;
        }
        match NonZeroU32::new(value) {
            Some(stored) => Some(Self(stored)),
            None => None,
        }
    }

    // r[impl pos.as_u32]
    // r[impl pos.must_use]
    /// The 1-based value. Infallible; this is the representation.
    #[inline]
    #[must_use]
    pub const fn as_u32(self) -> u32 {
        self.0.get()
    }

    // r[impl pos.as_i32]
    // r[impl pos.must_use]
    /// Raw value as `i32`. Infallible because all valid positions are `<= i32::MAX`.
    #[inline]
    #[must_use]
    #[expect(clippy::cast_possible_wrap, reason = "value ≤ i32::MAX by construction")]
    pub const fn as_i32(self) -> i32 {
        self.as_u32() as i32
    }

    // r[impl pos.as_usize]
    // r[impl pos.must_use]
    /// Convenience for indexing: returns the raw value as usize.
    #[inline]
    #[must_use]
    pub const fn as_usize(self) -> usize {
        self.as_u32() as usize
    }

    // r[impl pos.as_i64]
    // r[impl pos.must_use]
    /// Convenience for wider arithmetic: returns the raw value as i64.
    #[inline]
    #[must_use]
    pub const fn as_i64(self) -> i64 {
        self.as_u32() as i64
    }

    // r[impl pos.must_use]
    /// Convenience for 64-bit arithmetic: returns the raw value as u64.
    #[inline]
    #[must_use]
    pub const fn as_u64(self) -> u64 {
        self.as_u32() as u64
    }

    // r[impl pos.to_zero_based]
    // r[impl pos.explicit_conversion]
    /// Convert to 0-based. Infallible: 1-based values are in `1..=i32::MAX`,
    /// so subtracting 1 gives `0..=i32::MAX - 1`, always valid.
    ///
    /// A reinterpretation: a `Pos0` stores its value plus one, which is
    /// exactly the number a `Pos1` stores for the same base.
    #[inline]
    #[must_use]
    pub const fn to_zero_based(self) -> Pos0 {
        Pos0(self.0)
    }

    // r[impl pos.add_offset]
    // r[impl pos.two_types]
    /// Checked `position + offset`. `None` if the result leaves this type's
    /// range: below [`Pos1::MIN`] — so `Pos1(1) + (-1)` is `None`, not a
    /// 1-based zero — or above `i32::MAX`.
    #[inline]
    #[must_use]
    pub fn checked_add_offset(self, offset: Offset) -> Option<Self> {
        let result = i64::from(self.as_u32()).checked_add(offset.0)?;
        u32::try_from(result).ok().and_then(Self::new)
    }

    // r[impl pos.sub_offset]
    // r[impl pos.two_types]
    /// Checked `position - offset`. Same bounds as [`Self::checked_add_offset`].
    ///
    /// `Offset::new(i64::MIN)` has no negation, so subtracting it is `None`
    /// whatever the position.
    #[inline]
    #[must_use]
    pub fn checked_sub_offset(self, offset: Offset) -> Option<Self> {
        self.checked_add_offset(Offset(offset.0.checked_neg()?))
    }

    // r[impl pos.saturating_offset]
    // r[impl pos.two_types]
    /// `position + offset`, clamped to this type's range instead of failing:
    /// a result below [`Pos1::MIN`] becomes `MIN`, one above `i32::MAX`
    /// becomes [`Pos1::MAX`].
    ///
    /// For padding and trimming, where running off the start of a contig means
    /// the start of the contig. A caller that needs to *know* the bound was hit
    /// wants [`Self::checked_add_offset`].
    #[inline]
    #[must_use]
    #[expect(
        clippy::cast_sign_loss,
        clippy::cast_possible_truncation,
        reason = "clamped into 0..=POS_MAX - 1 before the cast"
    )]
    pub fn saturating_add_offset(self, offset: Offset) -> Self {
        // Clamp the distance above `MIN`, then build "one, plus that": a
        // `NonZeroU32` by construction, with no `Option` to resolve.
        let result = i64::from(self.as_u32()).saturating_add(offset.0);
        let above_min = result.saturating_sub(1).clamp(0, i64::from(POS_MAX).saturating_sub(1));
        Self(NonZeroU32::MIN.saturating_add(above_min as u32))
    }

    // r[impl pos.saturating_offset]
    // r[impl pos.two_types]
    /// `position - offset`, clamped to this type's range instead of failing.
    ///
    /// `Offset::new(i64::MIN)` has no negation; subtracting it saturates
    /// to [`Pos1::MAX`], since the most negative offset really does run off
    /// the far end.
    #[inline]
    #[must_use]
    pub fn saturating_sub_offset(self, offset: Offset) -> Self {
        self.saturating_add_offset(Offset(offset.0.saturating_neg()))
    }
}

// ---- Arithmetic ----

// r[impl pos.sub_pos]
// r[impl pos.no_add_pos]
// Pos0 - Pos0 = Offset; Add<Pos0> is deliberately not implemented.
impl Sub for Pos0 {
    type Output = Offset;
    #[inline]
    fn sub(self, rhs: Self) -> Offset {
        Offset(i64::from(self.as_u32()).wrapping_sub(i64::from(rhs.as_u32())))
    }
}

// r[impl pos.sub_pos]
// r[impl pos.no_add_pos]
// Pos1 - Pos1 = Offset; Add<Pos1> is deliberately not implemented.
impl Sub for Pos1 {
    type Output = Offset;
    #[inline]
    fn sub(self, rhs: Self) -> Offset {
        Offset(i64::from(self.as_u32()).wrapping_sub(i64::from(rhs.as_u32())))
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

impl From<Pos0> for u32 {
    #[inline]
    fn from(pos: Pos0) -> u32 {
        pos.as_u32()
    }
}

impl From<Pos0> for i32 {
    #[inline]
    fn from(pos: Pos0) -> i32 {
        pos.as_i32()
    }
}

impl From<Pos0> for u64 {
    #[inline]
    fn from(pos: Pos0) -> u64 {
        pos.as_u64()
    }
}

impl From<Pos0> for i64 {
    #[inline]
    fn from(pos: Pos0) -> i64 {
        pos.as_i64()
    }
}

impl From<Pos0> for usize {
    #[inline]
    fn from(pos: Pos0) -> usize {
        pos.as_usize()
    }
}

impl From<Pos1> for u32 {
    #[inline]
    fn from(pos: Pos1) -> u32 {
        pos.as_u32()
    }
}

impl From<Pos1> for i32 {
    #[inline]
    fn from(pos: Pos1) -> i32 {
        pos.as_i32()
    }
}

impl From<Pos1> for u64 {
    #[inline]
    fn from(pos: Pos1) -> u64 {
        pos.as_u64()
    }
}

impl From<Pos1> for i64 {
    #[inline]
    fn from(pos: Pos1) -> i64 {
        pos.as_i64()
    }
}

impl From<Pos1> for usize {
    #[inline]
    fn from(pos: Pos1) -> usize {
        pos.as_usize()
    }
}

// Every integer entry point below is `new` behind a widening or narrowing that
// cannot itself invent a position: a value that does not fit `u32` is out of
// range, and one that does is handed to that type's own check.

// r[impl pos.try_from]
impl TryFrom<u32> for Pos0 {
    type Error = PosOverflow;
    /// Fails above `i32::MAX`.
    #[inline]
    fn try_from(value: u32) -> Result<Self, PosOverflow> {
        Self::new(value).ok_or(PosOverflow)
    }
}

// r[impl pos.try_from]
impl TryFrom<i32> for Pos0 {
    type Error = PosOverflow;
    /// Fails below the floor: negative values are not 0-based positions.
    #[inline]
    fn try_from(value: i32) -> Result<Self, PosOverflow> {
        u32::try_from(value).ok().and_then(Self::new).ok_or(PosOverflow)
    }
}

// r[impl pos.try_from]
impl TryFrom<i64> for Pos0 {
    type Error = PosOverflow;
    /// Fails below zero or above `i32::MAX`.
    #[inline]
    fn try_from(value: i64) -> Result<Self, PosOverflow> {
        u32::try_from(value).ok().and_then(Self::new).ok_or(PosOverflow)
    }
}

// r[impl pos.try_from]
impl TryFrom<u64> for Pos0 {
    type Error = PosOverflow;
    /// Fails above `i32::MAX`.
    #[inline]
    fn try_from(value: u64) -> Result<Self, PosOverflow> {
        u32::try_from(value).ok().and_then(Self::new).ok_or(PosOverflow)
    }
}

// r[impl pos.try_from]
impl TryFrom<u32> for Pos1 {
    type Error = PosOverflow;
    /// Fails at `0` and above `i32::MAX`.
    #[inline]
    fn try_from(value: u32) -> Result<Self, PosOverflow> {
        Self::new(value).ok_or(PosOverflow)
    }
}

// r[impl pos.try_from]
impl TryFrom<i32> for Pos1 {
    type Error = PosOverflow;
    /// Fails below the floor: `0` and every negative value.
    #[inline]
    fn try_from(value: i32) -> Result<Self, PosOverflow> {
        u32::try_from(value).ok().and_then(Self::new).ok_or(PosOverflow)
    }
}

// r[impl pos.try_from]
impl TryFrom<i64> for Pos1 {
    type Error = PosOverflow;
    /// Fails below `1` or above `i32::MAX`.
    #[inline]
    fn try_from(value: i64) -> Result<Self, PosOverflow> {
        u32::try_from(value).ok().and_then(Self::new).ok_or(PosOverflow)
    }
}

// r[impl pos.try_from]
impl TryFrom<u64> for Pos1 {
    type Error = PosOverflow;
    /// Fails at `0` and above `i32::MAX`.
    #[inline]
    fn try_from(value: u64) -> Result<Self, PosOverflow> {
        u32::try_from(value).ok().and_then(Self::new).ok_or(PosOverflow)
    }
}

// ---- Query positions ----

// r[impl qpos.type]
// r[impl qpos.not_a_pos]
/// A 0-based offset into a read's sequence (query coordinates).
///
/// Deliberately **not** a [`Pos0`] or [`Pos1`]. Those name a place on the
/// reference; a `QPos` indexes a read. They are different spaces, and code
/// that resolves reference bases from a segment or window needs the former —
/// handed the latter it computes a wrong answer that still type-checks.
/// Keeping the two apart at the type level is the entire purpose of this
/// newtype.
///
/// Unlike a reference position there is no `i32::MAX` cap: a query offset is
/// bounded by the read length, and no format constrains it further.
/// Construction is therefore infallible.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
#[repr(transparent)]
pub struct QPos(u32);

impl QPos {
    /// The first base of a read.
    pub const ZERO: Self = Self(0);

    // r[impl qpos.new]
    /// Create a query offset.
    #[inline]
    #[must_use]
    pub const fn new(value: u32) -> Self {
        Self(value)
    }

    // r[impl qpos.get]
    /// The raw offset.
    #[inline]
    #[must_use]
    pub const fn get(self) -> u32 {
        self.0
    }

    // r[impl qpos.get]
    /// Convenience for indexing into a read's sequence or qualities.
    #[inline]
    #[must_use]
    pub const fn as_usize(self) -> usize {
        self.0 as usize
    }

    // r[impl qpos.checked]
    /// Checked offset + delta.
    #[inline]
    #[must_use]
    pub const fn checked_add(self, delta: u32) -> Option<Self> {
        match self.0.checked_add(delta) {
            Some(v) => Some(Self(v)),
            None => None,
        }
    }

    // r[impl qpos.checked]
    /// Checked offset - delta.
    #[inline]
    #[must_use]
    pub const fn checked_sub(self, delta: u32) -> Option<Self> {
        match self.0.checked_sub(delta) {
            Some(v) => Some(Self(v)),
            None => None,
        }
    }

    // r[impl qpos.saturating]
    /// Saturating offset + delta, for advancing a query cursor across CIGAR
    /// ops where a malformed record must not wrap.
    #[inline]
    #[must_use]
    pub const fn saturating_add(self, delta: u32) -> Self {
        Self(self.0.saturating_add(delta))
    }
}

impl fmt::Debug for QPos {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "QPos({})", self.0)
    }
}

impl fmt::Display for QPos {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

// ---- Display / Debug ----

// r[impl pos.derives]
// r[impl pos.two_types]
impl fmt::Debug for Pos0 {
    /// The logical value, never the stored one.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Pos0({})", self.as_u32())
    }
}

// r[impl pos.derives]
// r[impl pos.two_types]
impl fmt::Debug for Pos1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Pos1({})", self.as_u32())
    }
}

// r[impl pos.derives]
impl fmt::Display for Pos0 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_u32())
    }
}

// r[impl pos.derives]
impl fmt::Display for Pos1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_u32())
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

// ---- Serde ----

/// The "expected" half of a serde range error. Spelled out rather than built
/// from `POS_MAX`, which a `const _` above pins to the same number.
#[cfg(feature = "serde")]
const EXPECTED_POS0: &str = "a 0-based position in 0..=2147483647";

/// The 1-based counterpart of [`EXPECTED_POS0`].
#[cfg(feature = "serde")]
const EXPECTED_POS1: &str = "a 1-based position in 1..=2147483647";

// r[impl pos.serde]
// r[impl pos.two_types]
/// The wire format is the bare `u32`. Which coordinate system a field is in is
/// the field's type, not something the document restates.
#[cfg(feature = "serde")]
impl serde::Serialize for Pos0 {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_u32(self.as_u32())
    }
}

// r[impl pos.serde]
// r[impl pos.two_types]
/// Identical wire format to [`Pos0`]: the bare number.
#[cfg(feature = "serde")]
impl serde::Serialize for Pos1 {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_u32(self.as_u32())
    }
}

// r[impl pos.serde]
// r[impl pos.two_types]
/// The bound is re-checked on the way in, so a hand-edited or foreign document
/// cannot produce a position above `i32::MAX`.
#[cfg(feature = "serde")]
impl<'de> serde::Deserialize<'de> for Pos0 {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let value = u32::deserialize(deserializer)?;
        Self::new(value).ok_or_else(|| {
            serde::de::Error::invalid_value(
                serde::de::Unexpected::Unsigned(u64::from(value)),
                &EXPECTED_POS0,
            )
        })
    }
}

// r[impl pos.serde]
// r[impl pos.two_types]
/// The bound is re-checked on the way in, so a hand-edited or foreign document
/// cannot produce a `Pos1` of `0`.
#[cfg(feature = "serde")]
impl<'de> serde::Deserialize<'de> for Pos1 {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let value = u32::deserialize(deserializer)?;
        Self::new(value).ok_or_else(|| {
            serde::de::Error::invalid_value(
                serde::de::Unexpected::Unsigned(u64::from(value)),
                &EXPECTED_POS1,
            )
        })
    }
}

#[cfg(test)]
impl Pos0 {
    /// Test-only convenience: panics if value is out of range.
    pub fn at(value: u32) -> Self {
        Self::new(value).expect("test position out of range")
    }
}

#[cfg(test)]
impl Pos1 {
    /// Test-only convenience: panics if value is out of range.
    pub fn at(value: u32) -> Self {
        Self::new(value).expect("test position out of range")
    }
}

#[cfg(test)]
#[allow(clippy::arithmetic_side_effects, reason = "tests")]
mod tests {
    use hegel::prelude::*;

    use super::*;

    const I32_MAX_U32: u32 = i32::MAX as u32;

    // r[verify pos.to_one_based]
    // r[verify pos.to_zero_based]
    #[test]
    fn zero_based_roundtrip() {
        let z = Pos0::new(100).unwrap();
        let o = z.to_one_based().unwrap();
        assert_eq!(o.as_u32(), 101);
        assert_eq!(o.to_zero_based(), z);
    }

    // r[verify pos.to_zero_based]
    // r[verify pos.to_one_based]
    #[test]
    fn one_based_roundtrip() {
        let o = Pos1::new(1).unwrap();
        let z = o.to_zero_based();
        assert_eq!(z.as_u32(), 0);
        assert_eq!(z.to_one_based().unwrap(), o);
    }

    // r[verify pos.type]
    /// The one non-obvious thing about the representation: a `Pos0` and the
    /// `Pos1` of the same base store the same number, which is what makes
    /// both conversions a reinterpretation.
    #[test]
    fn the_same_base_has_the_same_representation() {
        let z = Pos0::new(100).unwrap();
        let o = Pos1::new(101).unwrap();
        assert_eq!(z.0, o.0);
        assert_eq!(z.to_one_based().unwrap().0, z.0);
        assert_eq!(o.to_zero_based().0, o.0);
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

    // r[verify pos.sub_pos]
    #[test]
    fn one_based_pos_minus_pos_is_offset() {
        let a = Pos1::new(100).unwrap();
        let b = Pos1::new(50).unwrap();
        assert_eq!((a - b).get(), 50);
        assert_eq!((b - a).get(), -50);
    }

    // r[verify pos.add_offset]
    #[test]
    fn pos_plus_offset() {
        let p = Pos0::new(10).unwrap();
        let q = p.checked_add_offset(Offset::new(5)).unwrap();
        assert_eq!(q.as_u32(), 15);
    }

    // r[verify pos.sub_offset]
    #[test]
    fn pos_minus_offset() {
        let p = Pos0::new(10).unwrap();
        let q = p.checked_sub_offset(Offset::new(3)).unwrap();
        assert_eq!(q.as_u32(), 7);
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
        assert_eq!(Pos0::try_from(0i32).unwrap().as_u32(), 0);
        assert_eq!(Pos1::try_from(1i32).unwrap().as_u32(), 1);
        assert_eq!(Pos0::try_from(100i32).unwrap().as_u32(), 100);
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
        assert_eq!(Pos0::try_from(0i64).unwrap().as_u32(), 0);
        assert_eq!(Pos1::try_from(1i64).unwrap().as_u32(), 1);
        assert_eq!(Pos0::try_from(100i64).unwrap().as_u32(), 100);
    }

    // r[verify pos.try_from]
    #[test]
    fn try_from_u64_rejects_overflow() {
        assert!(Pos0::try_from(u64::from(u32::MAX) + 1).is_err());
        assert!(Pos0::try_from(u64::from(u32::MAX)).is_err());
        assert!(Pos0::try_from(i32::MAX as u64 + 1).is_err());
        assert!(Pos1::try_from(0u64).is_err(), "Pos1 rejects 0");
        assert!(Pos1::try_from(i32::MAX as u64 + 1).is_err());
        assert_eq!(Pos1::try_from(1u64).unwrap().as_u32(), 1);
    }

    // r[verify pos.try_from]
    #[test]
    fn try_from_u32() {
        assert_eq!(Pos0::try_from(0u32).unwrap().as_u32(), 0);
        assert_eq!(Pos0::try_from(i32::MAX as u32).unwrap().as_u32(), i32::MAX as u32);
        assert!(Pos0::try_from(i32::MAX as u32 + 1).is_err());
        assert!(Pos1::try_from(0u32).is_err());
        assert_eq!(Pos1::try_from(1u32).unwrap().as_u32(), 1);
        assert!(Pos1::try_from(i32::MAX as u32 + 1).is_err());
    }

    #[test]
    fn as_i32_roundtrips() {
        let p = Pos0::new(42).unwrap();
        assert_eq!(p.as_i32(), 42);
        assert_eq!(Pos0::try_from(p.as_i32()).unwrap(), p);

        let max = Pos0::MAX;
        assert_eq!(max.as_i32(), i32::MAX);
        assert_eq!(Pos1::MAX.as_i32(), i32::MAX);
    }

    // r[verify pos.derives]
    #[test]
    fn ordering() {
        let a = Pos0::new(10).unwrap();
        let b = Pos0::new(20).unwrap();
        assert!(a < b);
        assert!(b > a);
        // Storing `value + 1` is monotone, so the derived `Ord` on the stored
        // number is the ordering of the positions — including at the floor.
        assert!(Pos0::ZERO < a);
        assert!(Pos1::MIN < Pos1::new(2).unwrap());
    }

    // r[verify pos.as_usize]
    #[test]
    fn as_usize() {
        let p = Pos0::new(42).unwrap();
        assert_eq!(p.as_usize(), 42);
        assert_eq!(Pos1::new(42).unwrap().as_usize(), 42);
    }

    // r[verify pos.as_i64]
    #[test]
    fn as_i64() {
        let p = Pos0::new(100).unwrap();
        assert_eq!(p.as_i64(), 100);
        assert_eq!(Pos1::new(100).unwrap().as_i64(), 100);
        assert_eq!(p.as_u64(), 100);
    }

    // r[verify pos.size]
    // r[verify pos.niche]
    #[test]
    fn size_is_u32() {
        assert_eq!(size_of::<Pos0>(), size_of::<u32>());
        assert_eq!(size_of::<Pos1>(), size_of::<u32>());
        // The `NonZeroU32` niche: an absent position costs nothing.
        assert_eq!(size_of::<Option<Pos0>>(), 4);
        assert_eq!(size_of::<Option<Pos1>>(), 4);
    }

    // r[verify pos.derives]
    #[test]
    fn debug_shows_coordinate_system() {
        let z = Pos0::new(42).unwrap();
        let o = Pos1::new(43).unwrap();
        assert_eq!(format!("{z:?}"), "Pos0(42)");
        assert_eq!(format!("{o:?}"), "Pos1(43)");
        // The logical value, not the stored one: these two share a
        // representation and must still print differently.
        assert_eq!(format!("{:?}", Pos0::ZERO), "Pos0(0)");
        assert_eq!(format!("{:?}", Pos1::MIN), "Pos1(1)");
        assert_eq!(format!("{z}"), "42");
        assert_eq!(format!("{o}"), "43");
    }

    // r[verify pos.min_max]
    #[test]
    fn min_and_max_are_the_types_bounds() {
        assert_eq!(Pos0::MIN, Pos0::ZERO);
        assert_eq!(Pos0::MIN.as_u32(), 0);
        assert_eq!(Pos1::MIN.as_u32(), 1);
        assert_eq!(Pos0::MAX.as_u32(), I32_MAX_U32);
        assert_eq!(Pos1::MAX.as_u32(), I32_MAX_U32);
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
        let p = Pos0::MAX;
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
        assert_eq!(p.checked_sub_offset(Offset::new(3)).unwrap().as_u32(), 7);
        assert!(p.checked_sub_offset(Offset::new(11)).is_none());
    }

    // r[verify pos.try_from]
    #[test]
    fn from_i32_infallible_roundtrip() {
        let p = Pos0::try_from(12345i32).unwrap();
        let back: i32 = p.into();
        assert_eq!(back, 12345);
        let q = Pos1::try_from(12345i32).unwrap();
        let back: i32 = q.into();
        assert_eq!(back, 12345);
    }

    // r[verify pos.to_one_based]
    // r[verify pos.to_zero_based]
    #[hegel::test]
    fn roundtrip_zero_to_one_and_back(tc: TestCase) {
        // i32::MAX itself can't round-trip (to_one_based would overflow).
        let v = tc.draw(gs::integers::<u32>().max_value(I32_MAX_U32 - 1));
        let z = Pos0::new(v).unwrap();
        let o = z.to_one_based().unwrap();
        assert_eq!(o.to_zero_based(), z);
    }

    // r[verify pos.to_zero_based]
    // r[verify pos.to_one_based]
    #[hegel::test]
    fn roundtrip_one_to_zero_and_back(tc: TestCase) {
        let v = tc.draw(gs::integers::<u32>().min_value(1).max_value(I32_MAX_U32));
        let o = Pos1::new(v).unwrap();
        let z = o.to_zero_based();
        assert_eq!(z.to_one_based().unwrap(), o);
    }

    // r[verify pos.to_one_based]
    // r[verify pos.to_zero_based]
    // r[verify pos.type]
    #[hegel::test]
    /// Both conversions are reinterpretations of the stored number, so what is
    /// worth checking is that they still agree with the `± 1` they replaced,
    /// and that the two types really do share a representation. A stored zero
    /// is not among the things that can go wrong: `NonZeroU32` has no such
    /// value, so no reinterpretation can produce one.
    fn conversion_agrees_with_the_arithmetic(tc: TestCase) {
        let v = tc.draw(gs::integers::<u32>().max_value(I32_MAX_U32));

        let zero = Pos0::new(v).expect("v <= i32::MAX");
        match zero.to_one_based() {
            Ok(one) => {
                assert_eq!(i64::from(one.as_u32()), i64::from(v) + 1);
                assert_eq!(one.0, zero.0, "same base, same representation");
                assert_eq!(one.to_zero_based(), zero);
            }
            Err(PosOverflow) => assert_eq!(v, I32_MAX_U32, "only the last position has no twin"),
        }

        if let Some(one) = Pos1::new(v) {
            let back = one.to_zero_based();
            assert_eq!(i64::from(back.as_u32()), i64::from(v) - 1);
            assert_eq!(back.0, one.0, "same base, same representation");
            assert_eq!(back.to_one_based(), Ok(one));
        } else {
            assert_eq!(v, 0, "only zero is not a 1-based position");
        }
    }

    // r[verify pos.add_offset]
    // r[verify pos.sub_offset]
    #[hegel::test]
    fn pos_plus_minus_offset_roundtrip(tc: TestCase) {
        let v = tc.draw(gs::integers::<u32>().max_value(999_999));
        let off = tc.draw(gs::integers::<i64>().min_value(-500_000).max_value(500_000));
        let result = i64::from(v) + off;
        if result >= 0 && result <= i64::from(I32_MAX_U32) {
            let p = Pos0::new(v).unwrap();
            let q = p.checked_add_offset(Offset::new(off)).unwrap();
            let r = q.checked_sub_offset(Offset::new(off)).unwrap();
            assert_eq!(r, p);
        }
    }

    // r[verify pos.sub_pos]
    #[hegel::test]
    fn pos_sub_pos_is_offset(tc: TestCase) {
        let coord = || gs::integers::<u32>().max_value(999_999);
        let a = tc.draw(coord());
        let b = tc.draw(coord());
        let pa = Pos0::new(a).unwrap();
        let pb = Pos0::new(b).unwrap();
        let off = pa - pb;
        assert_eq!(off.get(), i64::from(a) - i64::from(b));
    }

    #[hegel::test]
    fn new_never_accepts_above_i32_max(tc: TestCase) {
        let v = tc.draw(gs::integers::<u32>().min_value(I32_MAX_U32 + 1));
        assert!(Pos0::new(v).is_none());
        assert!(Pos1::new(v).is_none());
    }

    // r[verify pos.add_offset]
    #[hegel::test]
    /// `checked_add_offset` returns `Some` exactly when the sum is a position,
    /// and then it is the sum. Drawing the offset on both sides of zero, and
    /// far enough past `i32::MAX` to leave the range, covers both answers —
    /// the `if let` this used to be could not tell a correct `None` from a
    /// missing position.
    fn checked_add_is_exact_and_bounded(tc: TestCase) {
        let v = tc.draw(gs::integers::<u32>().max_value(I32_MAX_U32));
        let off = tc.draw(gs::integers::<i64>().min_value(-1000).max_value(1000));
        let p = Pos0::new(v).expect("v <= i32::MAX");
        let sum = i64::from(v) + off;
        let representable = (0..=i64::from(I32_MAX_U32)).contains(&sum);

        match p.checked_add_offset(Offset::new(off)) {
            Some(result) => {
                assert!(representable, "{v} + {off} = {sum} is not a position, but got {result:?}");
                assert_eq!(i64::from(result.as_u32()), sum);
                assert!(result.as_u32() <= I32_MAX_U32, "checked_add must not exceed i32::MAX");
            }
            None => assert!(!representable, "{v} + {off} = {sum} is a position, but got None"),
        }
    }

    // r[verify pos.as_i32]
    #[hegel::test]
    fn as_i32_always_valid(tc: TestCase) {
        let v = tc.draw(gs::integers::<u32>().max_value(I32_MAX_U32));
        let p = Pos0::new(v).unwrap();
        let i = p.as_i32();
        assert!(i >= 0);
        assert_eq!(i as u32, v);
    }

    // ---- The two types, one property at a time ----
    //
    // Each property is stated once for `P: Position` and run against both
    // types, so nothing can hold for `Pos0` and quietly not for `Pos1`.
    // That is how the 1-based floor went unenforced in `checked_add_offset`:
    // every test of it used `Pos0`, where the floor and "negative" coincide.

    /// The two position types, as one thing a property can be stated about.
    ///
    /// Test-only, and deliberately so: the library has no trait over `Pos0`
    /// and `Pos1` — writing the two out in full is what keeps each one's floor
    /// its own. The risk that buys is a property checked for one type and not
    /// the other, which is exactly what this trait puts back, in the only
    /// place where a shared abstraction cannot leak into a signature.
    trait Position: Copy + Ord + fmt::Debug + Sized {
        /// The first value this type admits.
        const FLOOR: u32;
        /// How the type names itself, for assertion messages.
        const NAME: &'static str;
        /// The type's own `MIN`.
        const FIRST: Self;
        /// The type's own `MAX`.
        const LAST: Self;

        fn build(value: u32) -> Option<Self>;
        fn raw(self) -> u32;
        fn raw_i64(self) -> i64;
        /// The `From<Self> for u32` impl, which must agree with `raw`.
        fn into_raw(self) -> u32;
        fn add_offset(self, offset: Offset) -> Option<Self>;
        fn sub_offset(self, offset: Offset) -> Option<Self>;
        fn sat_add_offset(self, offset: Offset) -> Self;
        fn sat_sub_offset(self, offset: Offset) -> Self;
        fn from_u32(value: u32) -> Result<Self, PosOverflow>;
        fn from_i32(value: i32) -> Result<Self, PosOverflow>;
        fn from_i64(value: i64) -> Result<Self, PosOverflow>;
        fn from_u64(value: u64) -> Result<Self, PosOverflow>;
    }

    impl Position for Pos0 {
        const FLOOR: u32 = 0;
        const NAME: &'static str = "Pos0";
        const FIRST: Self = Pos0::MIN;
        const LAST: Self = Pos0::MAX;

        fn build(value: u32) -> Option<Self> {
            Pos0::new(value)
        }
        fn raw(self) -> u32 {
            self.as_u32()
        }
        fn raw_i64(self) -> i64 {
            self.as_i64()
        }
        fn into_raw(self) -> u32 {
            u32::from(self)
        }
        fn add_offset(self, offset: Offset) -> Option<Self> {
            self.checked_add_offset(offset)
        }
        fn sub_offset(self, offset: Offset) -> Option<Self> {
            self.checked_sub_offset(offset)
        }
        fn sat_add_offset(self, offset: Offset) -> Self {
            self.saturating_add_offset(offset)
        }
        fn sat_sub_offset(self, offset: Offset) -> Self {
            self.saturating_sub_offset(offset)
        }
        fn from_u32(value: u32) -> Result<Self, PosOverflow> {
            Self::try_from(value)
        }
        fn from_i32(value: i32) -> Result<Self, PosOverflow> {
            Self::try_from(value)
        }
        fn from_i64(value: i64) -> Result<Self, PosOverflow> {
            Self::try_from(value)
        }
        fn from_u64(value: u64) -> Result<Self, PosOverflow> {
            Self::try_from(value)
        }
    }

    impl Position for Pos1 {
        const FLOOR: u32 = 1;
        const NAME: &'static str = "Pos1";
        const FIRST: Self = Pos1::MIN;
        const LAST: Self = Pos1::MAX;

        fn build(value: u32) -> Option<Self> {
            Pos1::new(value)
        }
        fn raw(self) -> u32 {
            self.as_u32()
        }
        fn raw_i64(self) -> i64 {
            self.as_i64()
        }
        fn into_raw(self) -> u32 {
            u32::from(self)
        }
        fn add_offset(self, offset: Offset) -> Option<Self> {
            self.checked_add_offset(offset)
        }
        fn sub_offset(self, offset: Offset) -> Option<Self> {
            self.checked_sub_offset(offset)
        }
        fn sat_add_offset(self, offset: Offset) -> Self {
            self.saturating_add_offset(offset)
        }
        fn sat_sub_offset(self, offset: Offset) -> Self {
            self.saturating_sub_offset(offset)
        }
        fn from_u32(value: u32) -> Result<Self, PosOverflow> {
            Self::try_from(value)
        }
        fn from_i32(value: i32) -> Result<Self, PosOverflow> {
            Self::try_from(value)
        }
        fn from_i64(value: i64) -> Result<Self, PosOverflow> {
            Self::try_from(value)
        }
        fn from_u64(value: u64) -> Result<Self, PosOverflow> {
            Self::try_from(value)
        }
    }

    /// A position anywhere in `P`'s range, biased towards the two ends so the
    /// boundary cases are drawn often rather than once in four billion.
    fn any_position<P: Position>(tc: &TestCase) -> P {
        let near_min = P::FLOOR.saturating_add(3);
        let near_max = POS_MAX.saturating_sub(3);
        let value = match tc.draw(gs::integers::<u8>().max_value(2)) {
            0 => tc.draw(gs::integers::<u32>().min_value(P::FLOOR).max_value(near_min)),
            1 => tc.draw(gs::integers::<u32>().min_value(near_max).max_value(POS_MAX)),
            _ => tc.draw(gs::integers::<u32>().min_value(P::FLOOR).max_value(POS_MAX)),
        };
        P::build(value).expect("drawn inside the type's range")
    }

    /// An offset that is small enough to land next to a boundary from
    /// `any_position` about half the time, and large enough to overshoot
    /// the whole range the rest of the time.
    fn any_offset(tc: &TestCase) -> Offset {
        let magnitude = if tc.draw(gs::booleans()) {
            tc.draw(gs::integers::<i64>().min_value(-8).max_value(8))
        } else {
            tc.draw(gs::integers::<i64>())
        };
        Offset::new(magnitude)
    }

    fn in_range<P: Position>(value: i64) -> bool {
        (i64::from(P::FLOOR)..=i64::from(POS_MAX)).contains(&value)
    }

    // r[verify pos.zero_new]
    // r[verify pos.one_new]
    // r[verify pos.two_types]
    fn new_admits_exactly_the_types_range<P: Position>(value: u32) {
        let pos = P::build(value);
        assert_eq!(pos.is_some(), in_range::<P>(i64::from(value)), "{value} for {}", P::NAME);
        if let Some(pos) = pos {
            assert_eq!(pos.raw(), value);
            assert_eq!(pos.into_raw(), value);
        }
    }

    #[hegel::test]
    fn new_admits_exactly_the_types_range_in_both_types(tc: TestCase) {
        let value = tc.draw(gs::integers::<u32>());
        new_admits_exactly_the_types_range::<Pos0>(value);
        new_admits_exactly_the_types_range::<Pos1>(value);
    }

    // r[verify pos.try_from]
    /// Every integer entry point agrees with `new` on the same number.
    fn try_from_agrees_with_new<P: Position>(value: i64) {
        let expected = u32::try_from(value).ok().and_then(P::build);
        assert_eq!(P::from_i64(value).ok(), expected, "i64 {value} for {}", P::NAME);
        if let Ok(v) = i32::try_from(value) {
            assert_eq!(P::from_i32(v).ok(), expected, "i32 {value} for {}", P::NAME);
        }
        if let Ok(v) = u32::try_from(value) {
            assert_eq!(P::from_u32(v).ok(), expected, "u32 {value} for {}", P::NAME);
        }
        if let Ok(v) = u64::try_from(value) {
            assert_eq!(P::from_u64(v).ok(), expected, "u64 {value} for {}", P::NAME);
        }
    }

    #[hegel::test]
    fn try_from_agrees_with_new_in_both_types(tc: TestCase) {
        // Around the two floors and the shared ceiling, plus anywhere at all.
        let value = match tc.draw(gs::integers::<u8>().max_value(2)) {
            0 => tc.draw(gs::integers::<i64>().min_value(-3).max_value(3)),
            1 => tc.draw(
                gs::integers::<i64>()
                    .min_value(i64::from(POS_MAX) - 3)
                    .max_value(i64::from(POS_MAX) + 3),
            ),
            _ => tc.draw(gs::integers::<i64>()),
        };
        try_from_agrees_with_new::<Pos0>(value);
        try_from_agrees_with_new::<Pos1>(value);
    }

    // r[verify pos.add_offset]
    // r[verify pos.sub_offset]
    /// `checked_add_offset` is plain integer arithmetic followed by the type's
    /// own range check — `Some` exactly when the sum is a position of that
    /// type, and then it is the sum.
    fn checked_offset_is_arithmetic_then_a_range_check<P: Position>(pos: P, offset: Offset) {
        let sum = pos.raw_i64().checked_add(offset.get());
        match pos.add_offset(offset) {
            Some(moved) => {
                assert_eq!(Some(moved.raw_i64()), sum, "{pos:?} + {offset:?}");
                assert!(in_range::<P>(moved.raw_i64()), "{moved:?} left {}'s range", P::NAME);
            }
            None => assert!(
                !sum.is_some_and(in_range::<P>),
                "{pos:?} + {offset:?} = {sum:?} is a position, but got None"
            ),
        }

        // Subtracting is adding the negation, where the negation exists.
        let via_add = offset.get().checked_neg().and_then(|n| pos.add_offset(Offset::new(n)));
        assert_eq!(pos.sub_offset(offset), via_add);
    }

    #[hegel::test]
    fn checked_offset_is_arithmetic_then_a_range_check_in_both_types(tc: TestCase) {
        let offset = any_offset(&tc);
        checked_offset_is_arithmetic_then_a_range_check(any_position::<Pos0>(&tc), offset);
        checked_offset_is_arithmetic_then_a_range_check(any_position::<Pos1>(&tc), offset);
    }

    // r[verify pos.add_offset]
    // r[verify pos.two_types]
    /// The regression the generic property exists for: a 1-based position
    /// at its floor, moved down by one, is not a position.
    #[test]
    fn checked_offset_respects_the_one_based_floor() {
        let first = Pos1::new(1).unwrap();
        assert_eq!(first.checked_add_offset(Offset::new(-1)), None);
        assert_eq!(first.checked_sub_offset(Offset::new(1)), None);
        assert_eq!(first.checked_add_offset(Offset::new(0)), Some(first));
        // The 0-based floor is where "negative" starts, so the same move is fine.
        assert_eq!(Pos0::new(1).unwrap().checked_sub_offset(Offset::new(1)), Some(Pos0::ZERO));
    }

    // r[verify pos.saturating_offset]
    /// Wherever the checked form succeeds the saturating form agrees with it;
    /// where it fails, the saturating form lands on the bound that was crossed.
    fn saturating_agrees_with_checked_or_lands_on_the_bound<P: Position>(pos: P, offset: Offset) {
        let saturated = pos.sat_add_offset(offset);
        match pos.add_offset(offset) {
            Some(checked) => assert_eq!(saturated, checked, "{pos:?} + {offset:?}"),
            None if offset.get() < 0 => assert_eq!(saturated, P::FIRST, "{pos:?} + {offset:?}"),
            None => assert_eq!(saturated, P::LAST, "{pos:?} + {offset:?}"),
        }
        assert!(in_range::<P>(saturated.raw_i64()));

        // The same for subtraction, whose one extra case is the offset with
        // no negation: subtracting it runs off the far end.
        let saturated = pos.sat_sub_offset(offset);
        match pos.sub_offset(offset) {
            Some(checked) => assert_eq!(saturated, checked, "{pos:?} - {offset:?}"),
            None if offset.get() > 0 => assert_eq!(saturated, P::FIRST, "{pos:?} - {offset:?}"),
            None => assert_eq!(saturated, P::LAST, "{pos:?} - {offset:?}"),
        }
    }

    #[hegel::test]
    fn saturating_agrees_with_checked_or_lands_on_the_bound_in_both_types(tc: TestCase) {
        let offset = any_offset(&tc);
        saturating_agrees_with_checked_or_lands_on_the_bound(any_position::<Pos0>(&tc), offset);
        saturating_agrees_with_checked_or_lands_on_the_bound(any_position::<Pos1>(&tc), offset);
    }

    // r[verify pos.saturating_offset]
    /// Moving further never lands earlier.
    fn saturating_is_monotone_in_the_offset<P: Position>(pos: P, a: Offset, b: Offset) {
        let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
        assert!(pos.sat_add_offset(lo) <= pos.sat_add_offset(hi));
    }

    #[hegel::test]
    fn saturating_is_monotone_in_the_offset_in_both_types(tc: TestCase) {
        let (a, b) = (any_offset(&tc), any_offset(&tc));
        saturating_is_monotone_in_the_offset(any_position::<Pos0>(&tc), a, b);
        saturating_is_monotone_in_the_offset(any_position::<Pos1>(&tc), a, b);
    }

    // r[verify pos.saturating_offset]
    // r[verify pos.min_max]
    #[test]
    fn saturating_offset_clamps_at_each_types_own_floor() {
        let zero = Pos0::new(2).unwrap();
        let one = Pos1::new(2).unwrap();
        assert_eq!(zero.saturating_sub_offset(Offset::new(10)), Pos0::ZERO);
        assert_eq!(one.saturating_sub_offset(Offset::new(10)), Pos1::MIN);
        assert_eq!(zero.saturating_add_offset(Offset::new(i64::MAX)), Pos0::MAX);
        assert_eq!(one.saturating_add_offset(Offset::new(i64::MAX)), Pos1::MAX);
        // `i64::MIN` has no negation: subtracting it saturates upwards, where
        // the checked form gives up.
        assert_eq!(zero.saturating_sub_offset(Offset::new(i64::MIN)), Pos0::MAX);
        assert_eq!(zero.checked_sub_offset(Offset::new(i64::MIN)), None);
        assert_eq!(one.saturating_sub_offset(Offset::new(i64::MIN)), Pos1::MAX);
        assert_eq!(one.checked_sub_offset(Offset::new(i64::MIN)), None);
    }

    #[cfg(feature = "serde")]
    mod serde_tests {
        use serde::Serialize;
        use serde::de::DeserializeOwned;

        use super::*;

        // r[verify pos.serde]
        #[test]
        fn a_position_serializes_as_its_bare_number() {
            assert_eq!(serde_json::to_string(&Pos0::new(41).unwrap()).unwrap(), "41");
            assert_eq!(serde_json::to_string(&Pos1::new(42).unwrap()).unwrap(), "42");
            // The two share a representation; the wire format is the logical
            // value, so the same base serializes differently in each system.
            assert_eq!(serde_json::to_string(&Pos0::ZERO).unwrap(), "0");
            assert_eq!(serde_json::to_string(&Pos1::MIN).unwrap(), "1");
        }

        // r[verify pos.serde]
        #[test]
        fn deserialization_rejects_what_the_constructor_rejects() {
            assert!(serde_json::from_str::<Pos1>("0").is_err(), "1-based has no zero");
            assert!(serde_json::from_str::<Pos0>("0").is_ok());
            let too_big = u64::from(POS_MAX) + 1;
            assert!(serde_json::from_str::<Pos0>(&too_big.to_string()).is_err());
            assert!(serde_json::from_str::<Pos1>(&too_big.to_string()).is_err());

            let err = serde_json::from_str::<Pos1>("0").unwrap_err().to_string();
            assert!(err.contains("1-based position in 1..=2147483647"), "{err}");
            let err = serde_json::from_str::<Pos0>(&too_big.to_string()).unwrap_err().to_string();
            assert!(err.contains("0-based position in 0..=2147483647"), "{err}");
        }

        // r[verify pos.serde]
        /// The wire format is the number and nothing else, and deserialization
        /// admits exactly what `new` admits — the two checks cannot drift.
        fn serde_matches_construction<P: Position + Serialize + DeserializeOwned>(value: u32) {
            let parsed = serde_json::from_str::<P>(&value.to_string()).ok();
            assert_eq!(parsed, P::build(value), "{value} for {}", P::NAME);
            let Some(pos) = parsed else { return };
            let json = serde_json::to_string(&pos).expect("a number always serializes");
            assert_eq!(json, value.to_string(), "no wrapper, no system tag");
            assert_eq!(serde_json::from_str::<P>(&json).ok(), Some(pos));
        }

        #[hegel::test]
        fn serde_admits_exactly_what_new_does(tc: TestCase) {
            let value = tc.draw(gs::integers::<u32>());
            serde_matches_construction::<Pos0>(value);
            serde_matches_construction::<Pos1>(value);
        }

        // r[verify pos.serde]
        /// `"5"`, `5.0` and `-1` are not positions; none is coerced.
        #[hegel::test]
        fn serde_rejects_anything_but_a_bare_unsigned_integer(tc: TestCase) {
            let value = tc.draw(gs::integers::<u32>().max_value(POS_MAX));
            for text in [format!("\"{value}\""), format!("{value}.0"), format!("-{}", value + 1)] {
                assert!(
                    serde_json::from_str::<Pos0>(&text).is_err(),
                    "{text} parsed as a position"
                );
                assert!(
                    serde_json::from_str::<Pos1>(&text).is_err(),
                    "{text} parsed as a position"
                );
            }
        }
    }

    // ---- QPos ----

    // r[verify pos.as_u32]
    #[test]
    fn as_u32_returns_raw_value() {
        assert_eq!(Pos0::new(I32_MAX_U32).unwrap().as_u32(), I32_MAX_U32);
        assert_eq!(Pos0::ZERO.as_u32(), 0);
        assert_eq!(Pos1::new(I32_MAX_U32).unwrap().as_u32(), I32_MAX_U32);
        assert_eq!(Pos1::MIN.as_u32(), 1);
    }

    // r[verify qpos.type]
    #[test]
    fn qpos_is_transparent_u32() {
        assert_eq!(size_of::<QPos>(), 4);
        assert_eq!(align_of::<QPos>(), 4);
    }

    // r[verify qpos.new]
    // r[verify qpos.get]
    #[test]
    fn qpos_roundtrip() {
        assert_eq!(QPos::new(0).get(), 0);
        assert_eq!(QPos::new(42).as_usize(), 42);
        assert_eq!(QPos::default(), QPos::ZERO);
    }

    // r[verify qpos.new]
    // Contrast with a reference position: no format caps a query offset, so
    // construction is infallible even above `i32::MAX` — the cap lives in the
    // read length, which the iterator layers check, not the type.
    #[test]
    fn qpos_accepts_above_i32_max() {
        assert_eq!(QPos::new(u32::MAX).get(), u32::MAX);
    }

    // r[verify qpos.checked]
    #[test]
    fn qpos_checked_arith() {
        assert_eq!(QPos::new(10).checked_add(5), Some(QPos::new(15)));
        assert_eq!(QPos::new(10).checked_sub(3), Some(QPos::new(7)));
        assert_eq!(QPos::new(10).checked_sub(11), None);
        assert_eq!(QPos::new(u32::MAX).checked_add(1), None);
    }

    // r[verify qpos.saturating]
    #[test]
    fn qpos_saturates_at_u32_max() {
        assert_eq!(QPos::new(u32::MAX).saturating_add(1), QPos::new(u32::MAX));
        assert_eq!(QPos::new(u32::MAX - 1).saturating_add(100), QPos::new(u32::MAX));
        assert_eq!(QPos::new(5).saturating_add(5), QPos::new(10), "no saturation below the cap");
    }
}
