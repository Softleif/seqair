//! [`Region`]: a contig and a range on it, for queries (`r[vcf_reader.region]`).

use crate::vcf::record_encoder::ContigId;
use seqair_types::{Pos0, Pos1, RegionString, RegionStringError, SmolStr};
use std::fmt;
use std::ops::{Bound, RangeBounds};

/// The contig a [`Region`] is on: by name, or by a file's [`ContigId`].
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum RegionContig {
    Name(SmolStr),
    Id(ContigId),
}

impl From<&str> for RegionContig {
    fn from(name: &str) -> Self {
        Self::Name(SmolStr::from(name))
    }
}

impl From<SmolStr> for RegionContig {
    fn from(name: SmolStr) -> Self {
        Self::Name(name)
    }
}

impl From<&ContigId> for RegionContig {
    fn from(id: &ContigId) -> Self {
        Self::Id(id.clone())
    }
}

impl From<ContigId> for RegionContig {
    fn from(id: ContigId) -> Self {
        Self::Id(id)
    }
}

impl RegionContig {
    /// The contig's name.
    pub fn name(&self) -> &str {
        match self {
            Self::Name(name) => name,
            Self::Id(id) => id.name(),
        }
    }
}

/// Where a [`Region`] ends.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RegionEnd {
    /// At the contig's last base.
    ContigEnd,
    /// Before this 0-based position (exclusive).
    Before(Pos0),
}

/// Why a [`Region`] could not be built.
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum RegionError {
    #[error("the region starts at {start} (0-based) but ends before {end}")]
    Reversed { start: Pos0, end: Pos0 },

    #[error("the region starts after the last representable position")]
    StartOverflow,

    #[error(transparent)]
    Parse(#[from] RegionStringError),
}

mod sealed {
    pub trait Sealed: Copy {
        /// The 0-based position of this coordinate's base.
        fn pos0(self) -> seqair_types::Pos0;
    }
}

/// A position type a [`Region`] range can be written in: [`Pos0`] (0-based,
/// `a..b` half-open) or [`Pos1`] (1-based, `a..=b` inclusive, as VCF and
/// region strings count). Sealed.
pub trait Coordinate: sealed::Sealed {}

impl sealed::Sealed for Pos0 {
    fn pos0(self) -> Pos0 {
        self
    }
}
impl Coordinate for Pos0 {}

impl sealed::Sealed for Pos1 {
    fn pos0(self) -> Pos0 {
        self.to_zero_based()
    }
}
impl Coordinate for Pos1 {}

// r[impl vcf_reader.region]
/// A contig and a range on it: a 0-based start and an end that is either the
/// contig's end or an exclusive 0-based position. Empty when they meet.
///
/// ```
/// use seqair::vcf::{Region, RegionEnd};
/// use seqair_types::{Pos0, Pos1};
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// // The same 100 bases, 1-based inclusive and 0-based half-open:
/// let one = Region::new("chr1", Pos1::try_from(101u32)?..=Pos1::try_from(200u32)?)?;
/// let zero = Region::new("chr1", Pos0::try_from(100u32)?..Pos0::try_from(200u32)?)?;
/// assert_eq!(one, zero);
/// assert_eq!(one.end(), RegionEnd::Before(Pos0::try_from(200u32)?));
/// // Open-ended, a whole contig, a region string:
/// assert_eq!(Region::new("chr1", Pos1::try_from(101u32)?..)?.end(), RegionEnd::ContigEnd);
/// assert_eq!(Region::contig("chr2"), "chr2".parse()?);
/// assert_eq!("chr1:101-200".parse::<Region>()?, one);
/// // Reversed is an error; `a..a` is empty.
/// assert!(Region::new("chr1", Pos0::try_from(5u32)?..Pos0::try_from(4u32)?).is_err());
/// assert!(Region::new("chr1", Pos0::ZERO..Pos0::ZERO)?.is_empty());
/// # Ok(()) }
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Region {
    contig: RegionContig,
    start: Pos0,
    end: RegionEnd,
}

impl Region {
    /// All of `contig`.
    pub fn contig(contig: impl Into<RegionContig>) -> Self {
        Self { contig: contig.into(), start: Pos0::ZERO, end: RegionEnd::ContigEnd }
    }

    /// `range` on `contig`, in 0-based ([`Pos0`]) or 1-based ([`Pos1`])
    /// coordinates. An unbounded start is the contig's first base, an
    /// unbounded end its last.
    ///
    /// # Errors
    /// [`RegionError::Reversed`] if the range starts after it ends;
    /// [`RegionError::StartOverflow`] for an excluded start at the last position.
    pub fn new<P: Coordinate>(
        contig: impl Into<RegionContig>,
        range: impl RangeBounds<P>,
    ) -> Result<Self, RegionError> {
        let after = |p: P| p.pos0().as_u32().checked_add(1).and_then(Pos0::new);
        let start = match range.start_bound() {
            Bound::Included(&p) => p.pos0(),
            Bound::Excluded(&p) => after(p).ok_or(RegionError::StartOverflow)?,
            Bound::Unbounded => Pos0::ZERO,
        };
        let end = match range.end_bound() {
            Bound::Excluded(&p) => RegionEnd::Before(p.pos0()),
            // Through the last representable position is through any contig's end.
            Bound::Included(&p) => after(p).map_or(RegionEnd::ContigEnd, RegionEnd::Before),
            Bound::Unbounded => RegionEnd::ContigEnd,
        };
        if let RegionEnd::Before(end) = end
            && end < start
        {
            return Err(RegionError::Reversed { start, end });
        }
        Ok(Self { contig: contig.into(), start, end })
    }

    /// The contig the region is on.
    pub fn target(&self) -> &RegionContig {
        &self.contig
    }

    /// The first base, 0-based.
    pub fn start(&self) -> Pos0 {
        self.start
    }

    pub fn end(&self) -> RegionEnd {
        self.end
    }

    /// Whether the region holds no base.
    pub fn is_empty(&self) -> bool {
        self.end == RegionEnd::Before(self.start)
    }

    /// The range as 0-based half-open `i64`s, the end `i64::MAX` at a contig's end.
    pub(crate) fn span(&self) -> (i64, i64) {
        let end = match self.end {
            RegionEnd::ContigEnd => i64::MAX,
            RegionEnd::Before(end) => end.as_i64(),
        };
        (self.start.as_i64(), end)
    }
}

/// `chr1:100-200` (1-based, inclusive), `chr1:100` (to the contig's end),
/// `chr1`. A region string whose start is after its end, which
/// [`RegionString`]'s parser rejects, is an empty region at its start.
impl From<&RegionString> for Region {
    fn from(r: &RegionString) -> Self {
        let start = r.start.map_or(Pos0::ZERO, Pos1::to_zero_based);
        let end = match r.end {
            None => RegionEnd::ContigEnd,
            Some(end) => RegionEnd::Before(Pos0::new(end.as_u32()).unwrap_or(Pos0::MAX).max(start)),
        };
        Self { contig: RegionContig::Name(r.chromosome.clone()), start, end }
    }
}

impl From<&Region> for Region {
    fn from(r: &Region) -> Self {
        r.clone()
    }
}

impl std::str::FromStr for Region {
    type Err = RegionError;

    /// A region string, as [`RegionString`] parses it.
    fn from_str(s: &str) -> Result<Self, RegionError> {
        Ok(Self::from(&s.parse::<RegionString>()?))
    }
}

/// `chr1:101-200`, 1-based inclusive as a region string reads.
impl fmt::Display for Region {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.contig.name())?;
        let start = self.start.as_u64().saturating_add(1);
        match self.end {
            RegionEnd::ContigEnd if self.start == Pos0::ZERO => Ok(()),
            RegionEnd::ContigEnd => write!(f, ":{start}"),
            RegionEnd::Before(end) => write!(f, ":{start}-{}", end.as_u64()),
        }
    }
}

/// How a query selects records (`r[vcf_reader.query.overlap+2]`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Match {
    /// Records whose reference span overlaps the region (`bcftools view -r`).
    #[default]
    Overlapping,
    /// Records whose start lies in the region (`bcftools view -t`): each
    /// record is in at most one of a set of disjoint regions.
    StartsIn,
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, clippy::arithmetic_side_effects, reason = "test code")]
    use super::*;
    use hegel::prelude::*;

    // r[verify vcf_reader.region]
    /// A 1-based inclusive range is the 0-based half-open range one lower at
    /// its start, and both agree with the region string that names it.
    #[hegel::test]
    fn one_and_zero_based_ranges_agree(tc: TestCase) {
        let a = tc.draw(gs::integers::<u32>().min_value(1).max_value(i32::MAX as u32));
        let b = tc.draw(gs::integers::<u32>().min_value(1).max_value(i32::MAX as u32));
        let one = Region::new("c", Pos1::new(a).unwrap()..=Pos1::new(b).unwrap());
        let zero = Region::new("c", Pos0::new(a - 1).unwrap()..Pos0::new(b).unwrap());
        if a > b + 1 {
            assert!(matches!(one, Err(RegionError::Reversed { .. })));
            assert!(matches!(zero, Err(RegionError::Reversed { .. })));
            return;
        }
        let (one, zero) = (one.unwrap(), zero.unwrap());
        assert_eq!(one, zero);
        assert_eq!(one.is_empty(), a == b + 1);
        assert_eq!(
            (one.start().as_u32(), one.end()),
            (a - 1, RegionEnd::Before(Pos0::new(b).unwrap()))
        );
        if a <= b {
            assert_eq!(format!("c:{a}-{b}").parse::<Region>().unwrap(), one);
            assert_eq!(one.to_string().parse::<Region>().unwrap(), one);
        }
    }

    #[test]
    fn bounds() {
        let p0 = |v| Pos0::new(v).unwrap();
        assert_eq!(Region::new::<Pos0>("c", ..).unwrap(), Region::contig("c"));
        assert_eq!(Region::new("c", ..=Pos0::MAX).unwrap().end(), RegionEnd::ContigEnd);
        assert_eq!(Region::new("c", ..=Pos1::MAX).unwrap().end(), RegionEnd::Before(Pos0::MAX));
        assert_eq!(
            Region::new("c", (Bound::Excluded(p0(4)), Bound::Unbounded)).unwrap().start(),
            p0(5)
        );
        assert!(matches!(
            Region::new("c", (Bound::Excluded(Pos0::MAX), Bound::Unbounded)),
            Err(RegionError::StartOverflow)
        ));
        assert_eq!(
            "c:7".parse::<Region>().unwrap(),
            Region::new("c", Pos1::new(7).unwrap()..).unwrap()
        );
        assert_eq!(Region::contig("c").to_string(), "c");
        assert_eq!("c:7".parse::<Region>().unwrap().to_string(), "c:7");
    }
}
