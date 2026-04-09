use crate::bam_flags::{BamFlags, consts::*};
use std::fmt;

/// Original top or bottom strand of a read
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[must_use]
pub enum Strand {
    /// Original top
    OT,
    /// Original bottom
    OB,
    /// Unknown
    Unknown,
}

impl Strand {
    /// Return `Some(self)` if strand is known, otherwise `None`
    pub fn ok(self) -> Option<Self> {
        match self {
            Strand::OT | Strand::OB => Some(self),
            Strand::Unknown => None,
        }
    }

    /// As symbol `+` (OT) or `-` (OB), or `.` (unknown)
    pub fn as_symbol(&self) -> &'static str {
        match self {
            Strand::OT => "+",
            Strand::OB => "-",
            Strand::Unknown => ".",
        }
    }
}

impl AsRef<str> for Strand {
    fn as_ref(&self) -> &str {
        match self {
            Strand::OT => "OT",
            Strand::OB => "OB",
            Strand::Unknown => "NA",
        }
    }
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_ref())
    }
}

/// Determine strand from BAM flags.
///
/// # Flags used
///
/// |   Flag | Decimal | Meaning                              |
/// | -----: | ------: | ------------------------------------ |
/// | `0x10` |      16 | Read is mapped to the reverse strand |
/// | `0x20` |      32 | Mate is mapped to the reverse strand |
/// | `0x40` |      64 | Read is first in pair                |
/// | `0x80` |     128 | Read is second in pair               |
#[allow(clippy::collapsible_else_if, reason = "clearer")]
pub fn strand_from_flags(flags: BamFlags) -> Strand {
    if !flags.is_set(FLAG_PAIRED) {
        // Unpaired read: strand determined solely by alignment direction
        if flags.is_set(FLAG_REVERSE) { Strand::OB } else { Strand::OT }
    } else if flags.is_set(FLAG_FIRST_IN_TEMPLATE) {
        // First in pair
        if flags.is_set(FLAG_REVERSE) { Strand::OB } else { Strand::OT }
    } else if flags.is_set(FLAG_SECOND_IN_TEMPLATE) {
        // Last in pair
        if flags.is_set(FLAG_MATE_REVERSE) { Strand::OB } else { Strand::OT }
    } else {
        Strand::Unknown
    }
}

/// Extension trait to get strand information from a BCF record
pub trait StrandFromRecord {
    /// Create a `Strand` from a BCF record
    fn strand(&self) -> Strand;
}

#[cfg(feature = "hts-compat")]
mod hts {
    use super::*;
    use rust_htslib::bam::Record;

    impl StrandFromRecord for Record {
        fn strand(&self) -> Strand {
            strand_from_flags(BamFlags::from(self.flags()))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn flags(raw: u16) -> BamFlags {
        BamFlags::from(raw)
    }

    #[test]
    fn test_various_records() {
        assert_eq!(strand_from_flags(flags(0x1 | 0x2 | 0x40 | 0x10)), Strand::OB); // First in pair, reverse strand
        assert_eq!(strand_from_flags(flags(0x1 | 0x2 | 0x80 | 0x20)), Strand::OB); // Second in pair, reverse strand
        assert_eq!(strand_from_flags(flags(0x1 | 0x2 | 0x40 | 0x20)), Strand::OT); // First in pair, mate reverse strand
        assert_eq!(strand_from_flags(flags(0x1 | 0x2 | 0x80 | 0x10)), Strand::OT); // Second in pair, mate reverse strand
        assert_eq!(strand_from_flags(flags(0x1 | 0x2 | 0x40)), Strand::OT); // First in pair, forward strand
        assert_eq!(strand_from_flags(flags(0x1 | 0x2 | 0x80)), Strand::OT); // Second in pair, forward strand
        assert_eq!(strand_from_flags(flags(0x00)), Strand::OT); // No flags set, ie top strand
        assert_eq!(strand_from_flags(flags(0x10)), Strand::OB); // No pairing flags, but read reverse strand -> OB
        assert_eq!(strand_from_flags(flags(0x01)), Strand::Unknown); // Paired but no first/second information
    }

    #[test]
    fn test_unpaired_mode() {
        assert_eq!(strand_from_flags(flags(0x00)), Strand::OT); // Single-end, forward
        assert_eq!(strand_from_flags(flags(0x10)), Strand::OB); // Single-end, reverse
        assert_eq!(strand_from_flags(flags(0x40 | 0x10)), Strand::OB); // Pair flags ignored in unpaired mode
        assert_eq!(strand_from_flags(flags(0x80 | 0x20)), Strand::OT); // Pair/mate flags ignored, only 0x10 matters
    }

    #[test]
    fn test_unpaired_mode_ignores_paired_flags() {
        assert_eq!(strand_from_flags(flags(0x01 | 0x02 | 0x20 | 0x40)), Strand::OT);
        assert_eq!(strand_from_flags(flags(0x01 | 0x02 | 0x20 | 0x40 | 0x10)), Strand::OB);
    }
}
