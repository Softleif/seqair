//! Resolve user-supplied contig identifiers (`u32` tid, contig name, or [`Tid`])
//! against a [`BamHeader`](crate::bam::BamHeader).

use crate::bam::BamHeader;
use seqair_types::SmolStr;

// r[impl unified.tid.newtype]
/// A validated BAM target id — an index into the header's target list.
///
/// Constructed only through [`ResolveTid::resolve_tid`]: the wrapper guarantees
/// the underlying `u32` is in range for the header it was resolved against.
/// `Clone` and `Copy` so callers can reuse the value without re-validating.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Tid(u32);

impl Tid {
    /// The raw target id as a `u32`.
    #[must_use]
    pub fn as_u32(self) -> u32 {
        self.0
    }

    /// The raw target id as a `usize`, for indexing into per-contig vectors.
    ///
    /// Always lossless on supported platforms: BAM tids fit in a `u32`, and
    /// `usize` is at least 32 bits everywhere we run.
    #[must_use]
    pub fn as_usize(self) -> usize {
        self.0 as usize
    }
}

impl From<Tid> for u32 {
    fn from(t: Tid) -> u32 {
        t.0
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
        if *self >= n {
            return Err(TidError::TidOutOfRange { tid: *self, n_targets: n });
        }
        Ok(Tid(*self))
    }
}

impl ResolveTid for str {
    fn resolve_tid(&self, header: &BamHeader) -> Result<Tid, TidError> {
        match header.tid(self) {
            Some(tid) => Ok(Tid(tid)),
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
