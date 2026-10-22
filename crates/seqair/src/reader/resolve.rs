//! Resolve user-supplied contig identifiers (`u32` tid, contig name, or [`Tid`])
//! against a [`BamHeader`](crate::bam::BamHeader).

use crate::bam::BamHeader;
use seqair_types::SmolStr;
use std::num::NonZeroU32;

// r[impl unified.tid.newtype]
/// A validated BAM target id — an index into the header's target list.
///
/// Constructed through [`ResolveTid::resolve_tid`], which checks the id
/// against the header, or crate-internally from a record's target field
/// ([`Tid::from_wire`]), which trusts the file. `Clone` and `Copy` so callers
/// can reuse the value without re-validating.
///
/// Stored as the id plus one, so `Option<Tid>` — a record's "unmapped" or
/// "no mate reference" — costs nothing over the `u32`, like [`Pos0`].
///
/// Pattern types would drop the plus one: on rustc 1.96,
/// `Option<pattern_type!(u32 is 0..=0x7FFF_FFFF)>` is four bytes, the values
/// above `i32::MAX` serving as the niche, and `new`/`as_u32` would lose their
/// `±1`. Worth switching once `pattern_types` is stable, or at least no
/// longer an `internal_features` lint with only `transmute` to build a value.
///
/// [`Pos0`]: seqair_types::Pos0
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Tid(NonZeroU32);

impl Tid {
    /// `None` for an id BAM cannot write (above `i32::MAX`).
    pub(crate) const fn new(index: u32) -> Option<Self> {
        if index > i32::MAX as u32 {
            return None;
        }
        match NonZeroU32::new(index.wrapping_add(1)) {
            Some(n) => Some(Self(n)),
            None => None,
        }
    }

    /// A target id as a BAM record stores it: `None` for `-1` (unmapped, no
    /// mate reference) and for any other negative value.
    pub(crate) fn from_wire(raw: i32) -> Option<Self> {
        u32::try_from(raw).ok().and_then(Self::new)
    }

    /// The BAM wire form of an optional target id: `-1` for `None`.
    pub(crate) fn to_wire(id: Option<Self>) -> i32 {
        // `new` rejects ids above `i32::MAX`, so this cannot wrap.
        id.map_or(-1, |t| t.as_u32().cast_signed())
    }

    /// The raw target id as a `u32`.
    #[must_use]
    pub fn as_u32(self) -> u32 {
        self.0.get().wrapping_sub(1)
    }

    /// The raw target id as a `usize`, for indexing into per-contig vectors.
    ///
    /// Always lossless on supported platforms: BAM tids fit in a `u32`, and
    /// `usize` is at least 32 bits everywhere we run.
    #[must_use]
    pub fn as_usize(self) -> usize {
        self.as_u32() as usize
    }
}

impl From<Tid> for u32 {
    fn from(t: Tid) -> u32 {
        t.as_u32()
    }
}

// r[impl unified.tid.newtype]
/// Error returned by [`ResolveTid::resolve_tid`].
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum TidError {
    #[error("contig '{name}' not found in header")]
    UnknownContig { name: SmolStr },
    #[error("tid {tid} out of range (header has {n_targets} targets)")]
    TidOutOfRange { tid: u32, n_targets: u32 },
}

// r[impl unified.tid.newtype]
/// Convert a user-facing contig identifier into a validated [`Tid`].
///
/// Implementations exist for `u32` (range-check against the header), `&str` /
/// `String` / `SmolStr` (lookup by name), and [`Tid`] (passthrough).
pub trait ResolveTid {
    /// Validate and resolve this value against `header`.
    fn resolve_tid(&self, header: &BamHeader) -> Result<Tid, TidError>;
}

impl ResolveTid for Tid {
    fn resolve_tid(&self, _header: &BamHeader) -> Result<Tid, TidError> {
        Ok(*self)
    }
}

impl ResolveTid for u32 {
    fn resolve_tid(&self, header: &BamHeader) -> Result<Tid, TidError> {
        let n = u32::try_from(header.target_count()).unwrap_or(u32::MAX);
        match Tid::new(*self) {
            Some(tid) if *self < n => Ok(tid),
            _ => Err(TidError::TidOutOfRange { tid: *self, n_targets: n }),
        }
    }
}

impl ResolveTid for str {
    fn resolve_tid(&self, header: &BamHeader) -> Result<Tid, TidError> {
        match header.tid(self).and_then(Tid::new) {
            Some(tid) => Ok(tid),
            None => Err(TidError::UnknownContig { name: self.into() }),
        }
    }
}

impl ResolveTid for &str {
    fn resolve_tid(&self, header: &BamHeader) -> Result<Tid, TidError> {
        (*self).resolve_tid(header)
    }
}

impl ResolveTid for String {
    fn resolve_tid(&self, header: &BamHeader) -> Result<Tid, TidError> {
        self.as_str().resolve_tid(header)
    }
}

impl ResolveTid for SmolStr {
    fn resolve_tid(&self, header: &BamHeader) -> Result<Tid, TidError> {
        self.as_str().resolve_tid(header)
    }
}

// `Option<Tid>` is a record field: the niche is what keeps `SlimRecord` at 96 bytes.
const _: () = assert!(std::mem::size_of::<Option<Tid>>() == 4);

#[cfg(test)]
mod tests {
    use super::Tid;
    use hegel::prelude::*;

    /// CRAI entries sort by `Option<Tid>`; that has to be the wire order, with
    /// the unplaced `-1` first.
    #[hegel::test]
    fn option_tid_orders_as_the_wire_ids(tc: TestCase) {
        let a = tc.draw(gs::integers::<i32>().min_value(-1));
        let b = tc.draw(gs::integers::<i32>().min_value(-1));
        assert_eq!(Tid::from_wire(a).cmp(&Tid::from_wire(b)), a.cmp(&b));
    }

    // r[verify unified.tid.newtype]
    #[test]
    fn wire_round_trip_and_bounds() {
        assert_eq!(Tid::from_wire(-1), None);
        assert_eq!(Tid::from_wire(-2), None);
        assert_eq!(Tid::to_wire(None), -1);
        for raw in [0, 1, 7, i32::MAX] {
            let tid = Tid::from_wire(raw).expect("non-negative ids are tids");
            assert_eq!(tid.as_u32(), raw.cast_unsigned());
            assert_eq!(Tid::to_wire(Some(tid)), raw);
        }
        assert_eq!(Tid::new(u32::MAX), None);
        assert_eq!(Tid::new(i32::MAX as u32 + 1), None);
    }
}
