use crate::bam_flags::{BamFlags, consts::*};
use std::fmt;

/// Original top or bottom strand of a read.
///
/// A read whose strand cannot be told is `Option<Strand>::None`, so code that
/// needs a strand gets one.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[must_use]
pub enum Strand {
    /// Original top
    OT,
    /// Original bottom
    OB,
}

impl Strand {
    /// As symbol `+` (OT) or `-` (OB)
    pub fn as_symbol(&self) -> &'static str {
        match self {
            Strand::OT => "+",
            Strand::OB => "-",
        }
    }

    /// The strand BAM flags imply; `None` for a paired read flagged neither
    /// first nor second in its template.
    ///
    /// # Flags used
    ///
    /// |   Flag | Decimal | Meaning                              |
    /// | -----: | ------: | ------------------------------------ |
    /// | `0x10` |      16 | Read is mapped to the reverse strand |
    /// | `0x20` |      32 | Mate is mapped to the reverse strand |
    /// | `0x40` |      64 | Read is first in pair                |
    /// | `0x80` |     128 | Read is second in pair               |
    pub fn from_flags(flags: BamFlags) -> Option<Self> {
        let reverse = if !flags.is_set(FLAG_PAIRED) || flags.is_set(FLAG_FIRST_IN_TEMPLATE) {
            flags.is_set(FLAG_REVERSE)
        } else if flags.is_set(FLAG_SECOND_IN_TEMPLATE) {
            flags.is_set(FLAG_MATE_REVERSE)
        } else {
            return None;
        };
        Some(if reverse { Strand::OB } else { Strand::OT })
    }
}

impl AsRef<str> for Strand {
    fn as_ref(&self) -> &str {
        match self {
            Strand::OT => "OT",
            Strand::OB => "OB",
        }
    }
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_ref())
    }
}

/// Extension trait to get strand information from a BCF record
pub trait StrandFromRecord {
    /// The strand the record's flags imply, see [`Strand::from_flags`]
    fn strand(&self) -> Option<Strand>;
}

#[cfg(feature = "hts-compat")]
mod hts {
    use super::*;
    use rust_htslib::bam::Record;

    impl StrandFromRecord for Record {
        fn strand(&self) -> Option<Strand> {
            Strand::from_flags(BamFlags::from(self.flags()))
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
        assert_eq!(Strand::from_flags(flags(0x1 | 0x2 | 0x40 | 0x10)), Some(Strand::OB)); // First in pair, reverse strand
        assert_eq!(Strand::from_flags(flags(0x1 | 0x2 | 0x80 | 0x20)), Some(Strand::OB)); // Second in pair, reverse strand
        assert_eq!(Strand::from_flags(flags(0x1 | 0x2 | 0x40 | 0x20)), Some(Strand::OT)); // First in pair, mate reverse strand
        assert_eq!(Strand::from_flags(flags(0x1 | 0x2 | 0x80 | 0x10)), Some(Strand::OT)); // Second in pair, mate reverse strand
        assert_eq!(Strand::from_flags(flags(0x1 | 0x2 | 0x40)), Some(Strand::OT)); // First in pair, forward strand
        assert_eq!(Strand::from_flags(flags(0x1 | 0x2 | 0x80)), Some(Strand::OT)); // Second in pair, forward strand
        assert_eq!(Strand::from_flags(flags(0x00)), Some(Strand::OT)); // No flags set, ie top strand
        assert_eq!(Strand::from_flags(flags(0x10)), Some(Strand::OB)); // No pairing flags, but read reverse strand -> OB
        assert_eq!(Strand::from_flags(flags(0x01)), None); // Paired but no first/second information
    }

    #[test]
    fn test_unpaired_mode() {
        assert_eq!(Strand::from_flags(flags(0x00)), Some(Strand::OT)); // Single-end, forward
        assert_eq!(Strand::from_flags(flags(0x10)), Some(Strand::OB)); // Single-end, reverse
        assert_eq!(Strand::from_flags(flags(0x40 | 0x10)), Some(Strand::OB)); // Pair flags ignored in unpaired mode
        assert_eq!(Strand::from_flags(flags(0x80 | 0x20)), Some(Strand::OT)); // Pair/mate flags ignored, only 0x10 matters
    }

    #[test]
    fn test_unpaired_mode_ignores_paired_flags() {
        assert_eq!(Strand::from_flags(flags(0x01 | 0x02 | 0x20 | 0x40)), Some(Strand::OT));
        assert_eq!(Strand::from_flags(flags(0x01 | 0x02 | 0x20 | 0x40 | 0x10)), Some(Strand::OB));
    }
}
