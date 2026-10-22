//! Resolve user-supplied contig identifiers (`u32` tid, contig name, or [`Tid`])
//! against a [`BamHeader`](crate::bam::BamHeader).

use crate::bam::BamHeader;
use seqair_types::SmolStr;
use std::num::NonZeroU32;

// r[impl unified.tid.newtype+2]
/// A validated BAM target id — an index into the header's target list.
///
/// Every `Tid` was checked against a header: by [`ResolveTid::resolve_tid`],
/// or crate-internally when a reader decodes a record's reference field
/// ([`Tid::from_wire`]), which rejects an id the header does not have. `Clone`
/// and `Copy` so callers can reuse the value without re-validating.
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
    const fn from_index(index: u32) -> Option<Self> {
        if index > i32::MAX as u32 {
            return None;
        }
        match NonZeroU32::new(index.wrapping_add(1)) {
            Some(n) => Some(Self(n)),
            None => None,
        }
    }

    /// Target `index` of `header`, if it has one.
    pub(crate) fn of_header(index: u32, header: &BamHeader) -> Option<Self> {
        usize::try_from(index)
            .is_ok_and(|i| i < header.target_count())
            .then(|| Self::from_index(index))
            .flatten()
    }

    /// A reference id as a record stores it: `None` for `-1` (unmapped, no
    /// mate reference), a `Tid` for one of `header`'s targets, and an error
    /// for anything else — as htslib's `sam_read1` rejects it.
    pub(crate) fn from_wire(raw: i32, header: &BamHeader) -> Result<Option<Self>, TidError> {
        if raw == -1 {
            return Ok(None);
        }
        let index = u32::try_from(raw).map_err(|_| TidError::Negative { value: raw })?;
        index.resolve_tid(header).map(Some)
    }

    /// The BAM wire form of an optional target id: `-1` for `None`.
    pub(crate) fn to_wire(id: Option<Self>) -> i32 {
        // `from_index` rejects ids above `i32::MAX`, so this cannot wrap.
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

/// The target index, as samtools prints it.
impl std::fmt::Display for Tid {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.as_u32().fmt(f)
    }
}

impl From<Tid> for u32 {
    fn from(t: Tid) -> u32 {
        t.as_u32()
    }
}

// r[impl unified.tid.newtype+2]
/// Error returned by [`ResolveTid::resolve_tid`].
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum TidError {
    #[error("contig '{name}' not found in header")]
    UnknownContig { name: SmolStr },
    #[error("tid {tid} out of range (header has {n_targets} targets)")]
    TidOutOfRange { tid: u32, n_targets: u32 },
    /// A record's reference id below `-1`, the only negative one that means
    /// something ("none").
    #[error("reference id {value} is negative; only -1 marks an absent reference")]
    Negative { value: i32 },
}

// r[impl unified.tid.newtype+2]
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
        Tid::of_header(*self, header).ok_or_else(|| TidError::TidOutOfRange {
            tid: *self,
            n_targets: u32::try_from(header.target_count()).unwrap_or(u32::MAX),
        })
    }
}

impl ResolveTid for str {
    fn resolve_tid(&self, header: &BamHeader) -> Result<Tid, TidError> {
        match header.tid(self) {
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
    use super::{Tid, TidError};
    use crate::bam::BamHeader;
    use hegel::prelude::*;

    fn header(n: usize) -> BamHeader {
        let text = (0..n).map(|i| format!("@SQ\tSN:c{i}\tLN:10\n")).collect::<Vec<_>>().concat();
        BamHeader::from_sam_text(&text).unwrap()
    }

    // r[verify unified.tid.newtype+2]
    /// A wire id is `None` for `-1`, a `Tid` below the header's count, and an
    /// error otherwise; a `Tid` writes back the id it was read from.
    #[hegel::test]
    fn a_wire_id_is_none_a_target_or_an_error(tc: TestCase) {
        // A header without references does not parse; one is the fewest.
        let n = tc.draw(gs::integers::<usize>().min_value(1).max_value(8));
        let raw = tc.draw(hegel::one_of!(
            gs::integers::<i32>().min_value(-3).max_value(10),
            gs::integers::<i32>(),
        ));
        let got = Tid::from_wire(raw, &header(n));
        match usize::try_from(raw) {
            Ok(i) if i < n => {
                let tid = got.unwrap().unwrap();
                assert_eq!(Tid::to_wire(Some(tid)), raw);
            }
            _ if raw == -1 => assert!(matches!(got, Ok(None))),
            _ if raw < 0 => {
                assert!(matches!(got, Err(TidError::Negative { value }) if value == raw));
            }
            _ => assert!(
                matches!(got, Err(TidError::TidOutOfRange { tid, .. }) if i64::from(tid) == i64::from(raw))
            ),
        }
    }

    /// `Tid` orders as the index it holds.
    #[hegel::test]
    fn tids_order_as_their_indices(tc: TestCase) {
        let h = header(8);
        let a = tc.draw(gs::integers::<u32>().max_value(7));
        let b = tc.draw(gs::integers::<u32>().max_value(7));
        assert_eq!(Tid::of_header(a, &h).cmp(&Tid::of_header(b, &h)), a.cmp(&b));
    }

    #[test]
    fn ids_bam_cannot_write_are_not_tids() {
        assert_eq!(Tid::from_index(u32::MAX), None);
        assert_eq!(Tid::from_index(i32::MAX.cast_unsigned() + 1), None);
        assert_eq!(
            Tid::to_wire(Some(Tid::from_index(i32::MAX.cast_unsigned()).unwrap())),
            i32::MAX
        );
    }
}
