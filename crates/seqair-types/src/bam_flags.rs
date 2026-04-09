//! Strongly-typed BAM flags that prevent accidentally mixing raw `u16` values.
//!
//! Named constants (e.g. [`FLAG_UNMAPPED`], [`FLAG_REVERSE`]) and predicate methods
//! (e.g. [`BamFlags::is_unmapped`]) replace scattered magic-number comparisons.

use self::consts::*;
use std::fmt;
use std::ops::{BitAnd, BitAndAssign, BitOr, BitOrAssign};

// r[impl flags.type]
// r[impl flags.copy]
/// Strongly-typed BAM flags newtype over `u16`.
///
/// Wraps the 16-bit SAM/BAM FLAG field, providing named predicates and mutators
/// instead of raw bitwise operations.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
#[repr(transparent)]
pub struct BamFlags(u16);

// r[impl flags.size]
const _: () = assert!(
    std::mem::size_of::<BamFlags>() == std::mem::size_of::<u16>(),
    "BamFlags must be exactly 2 bytes"
);

impl fmt::Debug for BamFlags {
    #[cfg_attr(coverage_nightly, coverage(off))]
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "BamFlags(0x{:04x})", self.0)
    }
}

// r[impl flags.display]
impl fmt::Display for BamFlags {
    #[cfg_attr(coverage_nightly, coverage(off))]
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(&self.0, f)
    }
}

// r[impl flags.format_traits]
impl fmt::LowerHex for BamFlags {
    #[cfg_attr(coverage_nightly, coverage(off))]
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::LowerHex::fmt(&self.0, f)
    }
}

impl fmt::UpperHex for BamFlags {
    #[cfg_attr(coverage_nightly, coverage(off))]
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::UpperHex::fmt(&self.0, f)
    }
}

impl fmt::Octal for BamFlags {
    #[cfg_attr(coverage_nightly, coverage(off))]
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Octal::fmt(&self.0, f)
    }
}

impl fmt::Binary for BamFlags {
    #[cfg_attr(coverage_nightly, coverage(off))]
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Binary::fmt(&self.0, f)
    }
}

impl BamFlags {
    // r[impl flags.empty]
    /// Create flags with no bits set.
    pub const fn empty() -> Self {
        Self(0)
    }

    // r[impl flags.raw]
    // r[impl flags.roundtrip]
    /// Return the underlying `u16` value.
    pub const fn raw(self) -> u16 {
        self.0
    }

    // r[impl flags.is_set]
    /// Test whether an arbitrary flag bit is set.
    #[inline]
    pub const fn is_set(self, flag: u16) -> bool {
        self.0 & flag != 0
    }

    // r[impl flags.predicates]
    /// Read is paired in sequencing.
    pub const fn is_paired(self) -> bool {
        self.is_set(FLAG_PAIRED)
    }

    /// Read is mapped in a proper pair.
    pub const fn is_proper_pair(self) -> bool {
        self.is_set(FLAG_PROPER_PAIR)
    }

    /// Read is unmapped.
    pub const fn is_unmapped(self) -> bool {
        self.is_set(FLAG_UNMAPPED)
    }

    /// Mate is unmapped.
    pub const fn is_mate_unmapped(self) -> bool {
        self.is_set(FLAG_MATE_UNMAPPED)
    }

    /// Read is mapped to the reverse strand.
    pub const fn is_reverse(self) -> bool {
        self.is_set(FLAG_REVERSE)
    }

    /// Mate is mapped to the reverse strand.
    pub const fn is_mate_reverse(self) -> bool {
        self.is_set(FLAG_MATE_REVERSE)
    }

    /// Read is the first in a template (pair).
    pub const fn is_first_in_template(self) -> bool {
        self.is_set(FLAG_FIRST_IN_TEMPLATE)
    }

    /// Read is the second (last) in a template (pair).
    pub const fn is_second_in_template(self) -> bool {
        self.is_set(FLAG_SECOND_IN_TEMPLATE)
    }

    /// Alignment is secondary (not primary).
    pub const fn is_secondary(self) -> bool {
        self.is_set(FLAG_SECONDARY)
    }

    /// Read fails platform/vendor quality checks.
    pub const fn is_failed_qc(self) -> bool {
        self.is_set(FLAG_FAILED_QC)
    }

    /// Read is a PCR or optical duplicate.
    pub const fn is_duplicate(self) -> bool {
        self.is_set(FLAG_DUPLICATE)
    }

    /// Alignment is supplementary.
    pub const fn is_supplementary(self) -> bool {
        self.is_set(FLAG_SUPPLEMENTARY)
    }

    // r[impl flags.set]
    /// Set the given flag bit(s).
    pub const fn set(&mut self, flag: u16) {
        self.0 |= flag;
    }

    // r[impl flags.unset]
    /// Clear the given flag bit(s).
    pub const fn unset(&mut self, flag: u16) {
        self.0 &= !flag;
    }

    // r[impl flags.with]
    /// Return a copy with the given flag bit(s) set.
    #[must_use]
    pub const fn with(self, flag: u16) -> Self {
        Self(self.0 | flag)
    }

    // r[impl flags.without]
    /// Return a copy with the given flag bit(s) cleared.
    #[must_use]
    pub const fn without(self, flag: u16) -> Self {
        Self(self.0 & !flag)
    }
}

// r[impl flags.bitor]
impl BitOr for BamFlags {
    type Output = Self;
    fn bitor(self, rhs: Self) -> Self {
        Self(self.0 | rhs.0)
    }
}

// r[impl flags.bitor_assign]
impl BitOrAssign for BamFlags {
    fn bitor_assign(&mut self, rhs: Self) {
        self.0 |= rhs.0;
    }
}

// r[impl flags.bitand]
impl BitAnd for BamFlags {
    type Output = Self;
    fn bitand(self, rhs: Self) -> Self {
        Self(self.0 & rhs.0)
    }
}

// r[impl flags.bitand_assign]
impl BitAndAssign for BamFlags {
    fn bitand_assign(&mut self, rhs: Self) {
        self.0 &= rhs.0;
    }
}

// r[impl flags.from]
impl From<u16> for BamFlags {
    fn from(raw: u16) -> Self {
        Self(raw)
    }
}

impl From<BamFlags> for u16 {
    fn from(flags: BamFlags) -> Self {
        flags.0
    }
}

// r[impl flags.constants]
/// SAM/BAM FLAG field bit constants
pub mod consts {
    /// Read is paired in sequencing.
    pub const FLAG_PAIRED: u16 = 0x1;
    /// Read is mapped in a proper pair.
    pub const FLAG_PROPER_PAIR: u16 = 0x2;
    /// Read is unmapped.
    pub const FLAG_UNMAPPED: u16 = 0x4;
    /// Mate is unmapped.
    pub const FLAG_MATE_UNMAPPED: u16 = 0x8;
    /// Read is mapped to the reverse strand.
    pub const FLAG_REVERSE: u16 = 0x10;
    /// Mate is mapped to the reverse strand.
    pub const FLAG_MATE_REVERSE: u16 = 0x20;
    /// Read is the first in a template (pair).
    pub const FLAG_FIRST_IN_TEMPLATE: u16 = 0x40;
    /// Read is the second (last) in a template (pair).
    pub const FLAG_SECOND_IN_TEMPLATE: u16 = 0x80;
    /// Alignment is secondary (not primary).
    pub const FLAG_SECONDARY: u16 = 0x100;
    /// Read fails platform/vendor quality checks.
    pub const FLAG_FAILED_QC: u16 = 0x200;
    /// Read is a PCR or optical duplicate.
    pub const FLAG_DUPLICATE: u16 = 0x400;
    /// Alignment is supplementary.
    pub const FLAG_SUPPLEMENTARY: u16 = 0x800;

    /// All 12 SAM-defined flag constants in bit order, for exhaustive testing.
    #[doc(hidden)]
    pub const ALL_FLAGS: [u16; 12] = [
        FLAG_PAIRED,
        FLAG_PROPER_PAIR,
        FLAG_UNMAPPED,
        FLAG_MATE_UNMAPPED,
        FLAG_REVERSE,
        FLAG_MATE_REVERSE,
        FLAG_FIRST_IN_TEMPLATE,
        FLAG_SECOND_IN_TEMPLATE,
        FLAG_SECONDARY,
        FLAG_FAILED_QC,
        FLAG_DUPLICATE,
        FLAG_SUPPLEMENTARY,
    ];
}

#[cfg(test)]
#[allow(clippy::type_complexity, reason = "tests")]
mod tests {
    use super::*;

    // r[verify flags.from]
    // r[verify flags.raw]
    // r[verify flags.roundtrip]
    #[test]
    fn roundtrip_all_u16() {
        for v in [0u16, 1, 0x4, 0xFFF, 0x8000, u16::MAX] {
            assert_eq!(BamFlags::from(v).raw(), v);
        }
    }

    // r[verify flags.empty]
    #[test]
    fn empty_is_zero() {
        let f = BamFlags::empty();
        assert_eq!(f.raw(), 0);
        assert!(!f.is_paired());
        assert!(!f.is_unmapped());
        assert!(!f.is_reverse());
    }

    // r[verify flags.type]
    // r[verify flags.size]
    #[test]
    fn size_is_two_bytes() {
        assert_eq!(std::mem::size_of::<BamFlags>(), 2);
    }

    // r[verify flags.copy]
    #[test]
    fn copy_semantics() {
        let a = BamFlags::from(0x63);
        let b = a; // Copy
        assert_eq!(a, b); // both still usable (PartialEq)

        let mut set = std::collections::HashSet::new();
        set.insert(a); // Hash
        assert!(set.contains(&b));
    }

    // r[verify flags.predicates]
    #[test]
    fn individual_predicates() {
        let predicates: Vec<(u16, fn(BamFlags) -> bool)> = vec![
            (FLAG_PAIRED, BamFlags::is_paired),
            (FLAG_PROPER_PAIR, BamFlags::is_proper_pair),
            (FLAG_UNMAPPED, BamFlags::is_unmapped),
            (FLAG_MATE_UNMAPPED, BamFlags::is_mate_unmapped),
            (FLAG_REVERSE, BamFlags::is_reverse),
            (FLAG_MATE_REVERSE, BamFlags::is_mate_reverse),
            (FLAG_FIRST_IN_TEMPLATE, BamFlags::is_first_in_template),
            (FLAG_SECOND_IN_TEMPLATE, BamFlags::is_second_in_template),
            (FLAG_SECONDARY, BamFlags::is_secondary),
            (FLAG_FAILED_QC, BamFlags::is_failed_qc),
            (FLAG_DUPLICATE, BamFlags::is_duplicate),
            (FLAG_SUPPLEMENTARY, BamFlags::is_supplementary),
        ];

        for (flag, pred) in &predicates {
            let f = BamFlags::from(*flag);
            assert!(pred(f), "predicate should be true for flag 0x{flag:04x}");

            // Other predicates should be false (each flag is a single bit)
            for (other_flag, other_pred) in &predicates {
                if other_flag != flag {
                    assert!(
                        !other_pred(f),
                        "predicate for 0x{other_flag:04x} should be false when only 0x{flag:04x} is set"
                    );
                }
            }
        }
    }

    // r[verify flags.constants]
    #[test]
    fn constants_are_single_bits() {
        for flag in ALL_FLAGS {
            assert_eq!(flag.count_ones(), 1, "FLAG 0x{flag:04x} is not a single bit");
        }
    }

    #[test]
    fn constants_are_unique() {
        for (i, a) in ALL_FLAGS.iter().enumerate() {
            for b in &ALL_FLAGS[i + 1..] {
                assert_ne!(a, b, "duplicate flag constant");
            }
        }
    }

    #[test]
    fn all_defined_flags_within_12_bits() {
        for flag in ALL_FLAGS {
            assert!(flag <= 0x800, "FLAG 0x{flag:04x} exceeds 12-bit range");
        }
    }

    // r[verify flags.set]
    // r[verify flags.unset]
    #[test]
    fn set_and_unset() {
        let mut f = BamFlags::empty();
        assert!(!f.is_unmapped());

        f.set(FLAG_UNMAPPED);
        assert!(f.is_unmapped());

        f.unset(FLAG_UNMAPPED);
        assert!(!f.is_unmapped());
    }

    #[test]
    fn set_preserves_existing() {
        let mut f = BamFlags::from(FLAG_PAIRED);
        f.set(FLAG_REVERSE);
        assert!(f.is_paired());
        assert!(f.is_reverse());
    }

    #[test]
    fn unset_preserves_other_bits() {
        let mut f = BamFlags::from(FLAG_PAIRED | FLAG_REVERSE);
        f.unset(FLAG_PAIRED);
        assert!(!f.is_paired());
        assert!(f.is_reverse());
    }

    // r[verify flags.with]
    // r[verify flags.without]
    #[test]
    fn with_and_without() {
        let f = BamFlags::empty().with(FLAG_UNMAPPED).with(FLAG_REVERSE);
        assert!(f.is_unmapped());
        assert!(f.is_reverse());

        let g = f.without(FLAG_UNMAPPED);
        assert!(!g.is_unmapped());
        assert!(g.is_reverse());
        // Original unchanged (Copy semantics)
        assert!(f.is_unmapped());
    }

    // r[verify flags.bitor]
    #[test]
    fn bitor_combines_flags() {
        let a = BamFlags::from(FLAG_PAIRED);
        let b = BamFlags::from(FLAG_REVERSE);
        let c = a | b;
        assert!(c.is_paired());
        assert!(c.is_reverse());
        assert_eq!(c.raw(), FLAG_PAIRED | FLAG_REVERSE);
    }

    // r[verify flags.bitor_assign]
    #[test]
    fn bitor_assign() {
        let mut f = BamFlags::from(FLAG_PAIRED);
        f |= BamFlags::from(FLAG_REVERSE);
        assert!(f.is_paired());
        assert!(f.is_reverse());
    }

    // r[verify flags.bitand]
    #[test]
    fn bitand_intersects_flags() {
        let a = BamFlags::from(FLAG_PAIRED | FLAG_REVERSE);
        let b = BamFlags::from(FLAG_PAIRED | FLAG_UNMAPPED);
        let c = a & b;
        assert!(c.is_paired());
        assert!(!c.is_reverse());
        assert!(!c.is_unmapped());
    }

    // r[verify flags.bitand_assign]
    #[test]
    fn bitand_assign() {
        let mut f = BamFlags::from(FLAG_PAIRED | FLAG_REVERSE);
        f &= BamFlags::from(FLAG_PAIRED | FLAG_UNMAPPED);
        assert!(f.is_paired());
        assert!(!f.is_reverse());
    }

    // r[verify flags.display]
    #[test]
    fn display_is_decimal() {
        assert_eq!(format!("{}", BamFlags::from(99)), "99");
        assert_eq!(format!("{}", BamFlags::from(0)), "0");
        assert_eq!(format!("{}", BamFlags::from(0xFFF)), "4095");
    }

    #[test]
    fn display_respects_width_and_fill() {
        assert_eq!(format!("{:>6}", BamFlags::from(99)), "    99");
        assert_eq!(format!("{:0>6}", BamFlags::from(99)), "000099");
    }

    // r[verify flags.format_traits]
    #[test]
    fn format_lower_hex() {
        assert_eq!(format!("{:x}", BamFlags::from(0xFFF)), "fff");
        assert_eq!(format!("{:#06x}", BamFlags::from(0x63)), "0x0063");
    }

    #[test]
    fn format_upper_hex() {
        assert_eq!(format!("{:X}", BamFlags::from(0xFFF)), "FFF");
        assert_eq!(format!("{:#06X}", BamFlags::from(0x63)), "0x0063");
    }

    #[test]
    fn format_binary() {
        assert_eq!(format!("{:b}", BamFlags::from(0x5)), "101");
        assert_eq!(format!("{:#018b}", BamFlags::from(0xFFF)), "0b0000111111111111");
    }

    #[test]
    fn format_octal() {
        assert_eq!(format!("{:o}", BamFlags::from(0x1FF)), "777");
        assert_eq!(format!("{:#o}", BamFlags::from(8)), "0o10");
    }

    #[test]
    fn debug_is_hex() {
        assert_eq!(format!("{:?}", BamFlags::from(0x63)), "BamFlags(0x0063)");
    }

    // r[verify flags.is_set]
    #[test]
    fn is_set_arbitrary_bit() {
        let f = BamFlags::from(0x8000);
        assert!(f.is_set(0x8000));
        assert!(!f.is_set(0x4000));
    }

    #[test]
    fn typical_bam_flag_combinations() {
        // 0x63 = 0x01|0x02|0x20|0x40 = paired + proper_pair + mate_reverse + first_in_template
        let r1 = BamFlags::from(0x63);
        assert!(r1.is_paired());
        assert!(r1.is_proper_pair());
        assert!(r1.is_mate_reverse());
        assert!(r1.is_first_in_template());
        assert!(!r1.is_unmapped());
        assert!(!r1.is_reverse());

        // Secondary alignment
        let sec = BamFlags::from(FLAG_SECONDARY);
        assert!(sec.is_secondary());
        assert!(!sec.is_supplementary());

        // Supplementary alignment
        let sup = BamFlags::from(FLAG_SUPPLEMENTARY);
        assert!(sup.is_supplementary());
        assert!(!sup.is_secondary());
    }

    #[test]
    fn set_multiple_bits_at_once() {
        let mut f = BamFlags::empty();
        f.set(FLAG_PAIRED | FLAG_PROPER_PAIR | FLAG_REVERSE);
        assert!(f.is_paired());
        assert!(f.is_proper_pair());
        assert!(f.is_reverse());
        assert!(!f.is_unmapped());
    }

    #[test]
    fn unset_multiple_bits_at_once() {
        let mut f = BamFlags::from(FLAG_PAIRED | FLAG_PROPER_PAIR | FLAG_REVERSE);
        f.unset(FLAG_PAIRED | FLAG_REVERSE);
        assert!(!f.is_paired());
        assert!(f.is_proper_pair());
        assert!(!f.is_reverse());
    }

    #[test]
    fn with_is_idempotent() {
        let f = BamFlags::from(FLAG_PAIRED);
        let g = f.with(FLAG_PAIRED);
        assert_eq!(f, g);
    }

    #[test]
    fn without_on_absent_bit_is_noop() {
        let f = BamFlags::from(FLAG_PAIRED);
        let g = f.without(FLAG_REVERSE);
        assert_eq!(f, g);
    }

    mod proptests {
        use super::*;
        use proptest::prelude::*;

        // r[verify flags.from]
        // r[verify flags.raw]
        // r[verify flags.roundtrip]
        proptest! {
            #[test]
            fn roundtrip(v: u16) {
                prop_assert_eq!(BamFlags::from(v).raw(), v);
            }

            /// `with` and `set` are two code paths to the same result.
            #[test]
            fn with_equals_set(v: u16, bit: u16) {
                let f = BamFlags::from(v);
                let via_with = f.with(bit);
                let mut via_set = f;
                via_set.set(bit);
                prop_assert_eq!(via_with, via_set);
            }

            /// `without` and `unset` are two code paths to the same result.
            #[test]
            fn without_equals_unset(v: u16, bit: u16) {
                let f = BamFlags::from(v);
                let via_without = f.without(bit);
                let mut via_unset = f;
                via_unset.unset(bit);
                prop_assert_eq!(via_without, via_unset);
            }

            /// set then unset clears only the target bit, leaves others unchanged.
            #[test]
            fn set_then_unset_clears_only_target(v: u16, bit: u16) {
                let mut f = BamFlags::from(v);
                f.set(bit);
                f.unset(bit);
                prop_assert_eq!(f.raw(), v & !bit);
            }

            /// with(bit).without(bit) on a value where bit was clear restores original.
            #[test]
            fn with_without_roundtrip_when_bit_was_clear(v: u16, bit_idx in 0u16..16) {
                let bit = 1u16 << bit_idx;
                let cleared = BamFlags::from(v & !bit);
                let toggled = cleared.with(bit).without(bit);
                prop_assert_eq!(toggled, cleared);
            }

            /// Display roundtrips through decimal parse.
            #[test]
            fn display_roundtrips_through_parse(v: u16) {
                let s = format!("{}", BamFlags::from(v));
                let parsed = s.parse::<u16>();
                prop_assert!(parsed.is_ok(), "failed to parse {:?}", s);
                prop_assert_eq!(parsed.unwrap(), v);
            }

            /// Undefined bits above 0x800 don't trigger any named predicate.
            #[test]
            fn undefined_bits_dont_trigger_predicates(extra_bits in 0u16..=0xF000) {
                let f = BamFlags::from(extra_bits & 0xF000);
                prop_assert!(!f.is_paired());
                prop_assert!(!f.is_proper_pair());
                prop_assert!(!f.is_unmapped());
                prop_assert!(!f.is_mate_unmapped());
                prop_assert!(!f.is_reverse());
                prop_assert!(!f.is_mate_reverse());
                prop_assert!(!f.is_first_in_template());
                prop_assert!(!f.is_second_in_template());
                prop_assert!(!f.is_secondary());
                prop_assert!(!f.is_failed_qc());
                prop_assert!(!f.is_duplicate());
                prop_assert!(!f.is_supplementary());
            }
        }
    }
}
