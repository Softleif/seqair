// `IntoSegmentTarget` is sealed via `private::IntoSegmentTargetSealed`, which
// uses the crate-private `ResolvedRange`. Rust's reachability lint flags this
// as "leaking" the private type, but the seal is the whole point — external
// crates cannot implement the trait, so they cannot ever construct or observe
// a `ResolvedRange`. Allow the warning at module scope.
#![allow(
    private_interfaces,
    reason = "ResolvedRange is intentionally crate-private; the IntoSegmentTarget trait is sealed"
)]

//! Tile-based planning for [`Readers::pileup`](super::Readers::pileup).
//!
//! `Readers::pileup` only ever takes a [`Segment`]. The only way to obtain one
//! is to call [`Readers::segments`](super::Readers::segments) with a target and
//! [`SegmentOptions`]. This shape makes "pile up the entire chromosome in one
//! call" a deliberate choice — pick a `max_len`, iterate the segments, drive
//! `pileup` on each.
//!
//! See `r[unified.segment_struct]`, `r[unified.segment_overlap]`,
//! `r[unified.segment_options]`, `r[unified.into_segment_target]`,
//! `r[unified.readers_segments]`, `r[unified.readers_pileup]`.

use crate::bam::BamHeader;
use seqair_types::{Pos0, RegionString, SmolStr};
use std::num::NonZeroU32;

use super::{
    ReaderError,
    resolve::{ResolveTid, Tid},
};

/// Errors returned by [`Segment::new`] when invariants are violated.
///
/// The constructor is `pub(crate)`, so external callers never hit these —
/// they're surfaced for the iterator's tests and any future internal caller.
#[derive(Debug, thiserror::Error, PartialEq, Eq)]
#[non_exhaustive]
pub enum SegmentInvariantError {
    #[error("segment start {start} > end {end}")]
    StartAfterEnd { start: u64, end: u64 },
    #[error("segment end {end} exceeds contig_last_pos {contig_last_pos}")]
    EndPastContig { end: u64, contig_last_pos: u64 },
    #[error(
        "overlap_start ({overlap_start}) + overlap_end ({overlap_end}) \
         must be < tile length ({len})"
    )]
    OverlapTooLarge { overlap_start: u32, overlap_end: u32, len: u32 },
}

// r[impl unified.segment_struct]
/// One bounded tile of a genomic region, ready for [`Readers::pileup`].
///
/// Carries a pre-resolved [`Tid`], the contig name, an inclusive
/// `[start, end]` range, and explicit overlap with neighboring tiles.
/// `Segment` has no public constructor — obtain instances from
/// [`Readers::segments`](super::Readers::segments).
///
/// `Segment` is intentionally lightweight (no FASTA cache, no buffers) so
/// callers can collect or send segments freely between threads.
///
/// [`Readers::pileup`]: super::Readers::pileup
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Segment {
    tid: Tid,
    contig: SmolStr,
    /// Inclusive 0-based start of this tile.
    start: Pos0,
    /// Inclusive 0-based end of this tile.
    end: Pos0,
    /// Number of bases covered by this tile (= `end - start + 1`). Computed
    /// once at construction so accessors don't have to re-derive it.
    len: u32,
    overlap_start: u32,
    overlap_end: u32,
    /// Inclusive 0-based last position of the contig this tile belongs to.
    contig_last_pos: Pos0,
}

impl Segment {
    /// Construct a segment. `pub(crate)` so only the iterator and tests in
    /// this crate can build them.
    ///
    /// Returns [`SegmentInvariantError`] when:
    ///
    /// * `start > end`
    /// * `end > contig_last_pos`
    /// * `overlap_start + overlap_end >= len()` — would make `core_range()`
    ///   cover the entire tile, breaking neighbor-dedupe downstream.
    pub(crate) fn new(
        tid: Tid,
        contig: SmolStr,
        start: Pos0,
        end: Pos0,
        overlap_start: u32,
        overlap_end: u32,
        contig_last_pos: Pos0,
    ) -> Result<Self, SegmentInvariantError> {
        if start > end {
            return Err(SegmentInvariantError::StartAfterEnd {
                start: start.as_u64(),
                end: end.as_u64(),
            });
        }
        if end > contig_last_pos {
            return Err(SegmentInvariantError::EndPastContig {
                end: end.as_u64(),
                contig_last_pos: contig_last_pos.as_u64(),
            });
        }
        // start <= end <= i32::MAX (both are Pos0), so end - start + 1 fits
        // in u32 exactly: max value is i32::MAX + 1 = 2_147_483_648, and that
        // round-trips through u32 since u32::MAX = 4_294_967_295. Use checked
        // arithmetic so a future Pos0 representation widening can't quietly
        // wrap.
        let span_u64 =
            end.as_u64().checked_sub(start.as_u64()).and_then(|d| d.checked_add(1)).ok_or(
                SegmentInvariantError::StartAfterEnd { start: start.as_u64(), end: end.as_u64() },
            )?;
        let len = u32::try_from(span_u64).map_err(|_| SegmentInvariantError::EndPastContig {
            end: end.as_u64(),
            contig_last_pos: contig_last_pos.as_u64(),
        })?;
        let overlap_total = u64::from(overlap_start).saturating_add(u64::from(overlap_end));
        if overlap_total >= u64::from(len) {
            return Err(SegmentInvariantError::OverlapTooLarge { overlap_start, overlap_end, len });
        }
        Ok(Self { tid, contig, start, end, len, overlap_start, overlap_end, contig_last_pos })
    }

    /// The validated target id of the contig this tile lies on.
    #[must_use]
    pub fn tid(&self) -> Tid {
        self.tid
    }

    /// The contig name, as carried in the BAM header.
    #[must_use]
    pub fn contig(&self) -> &SmolStr {
        &self.contig
    }

    /// Inclusive 0-based start of the tile (first pileup position).
    #[must_use]
    pub fn start(&self) -> Pos0 {
        self.start
    }

    /// Inclusive 0-based end of the tile (last pileup position).
    #[must_use]
    pub fn end(&self) -> Pos0 {
        self.end
    }

    /// Number of bases covered by this tile (= `end - start + 1`). Always >= 1.
    #[must_use]
    pub fn len(&self) -> u32 {
        self.len
    }

    /// Always `false`. A `Segment` is built only by the iterator, which
    /// never yields zero-length tiles (`r[unified.readers_segments]`). Kept
    /// to satisfy `clippy::len_without_is_empty`.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        false
    }

    /// Bases shared with the previous tile of the same contig. `0` for the
    /// first tile of a contig.
    #[must_use]
    pub fn overlap_start(&self) -> u32 {
        self.overlap_start
    }

    /// Bases shared with the next tile of the same contig. `0` for the last
    /// tile of a contig.
    #[must_use]
    pub fn overlap_end(&self) -> u32 {
        self.overlap_end
    }

    /// Inclusive 0-based last position of the contig this tile lives on.
    #[must_use]
    pub fn contig_last_pos(&self) -> Pos0 {
        self.contig_last_pos
    }

    /// True when this tile's `start` is exactly the contig's first base
    /// (`Pos0::ZERO`).
    ///
    /// Note: this checks the contig boundary, not the requested target
    /// range — for a sub-range query like `("chr1", 1000, 2000)` the first
    /// yielded segment will have `start == 1000` and this method returns
    /// `false`. If you want "is this the first segment of the user's
    /// query?", inspect the iterator directly.
    #[must_use]
    pub fn starts_at_contig_start(&self) -> bool {
        self.start == Pos0::ZERO
    }

    /// True when this tile's `end` is exactly the contig's last base
    /// (`contig_last_pos`). See [`starts_at_contig_start`] for the analogous
    /// caveat about sub-range queries.
    ///
    /// [`starts_at_contig_start`]: Self::starts_at_contig_start
    #[must_use]
    pub fn ends_at_contig_end(&self) -> bool {
        self.end == self.contig_last_pos
    }

    // r[impl unified.segment_overlap]
    /// The inclusive sub-range "owned" by this tile — `[start, end]` shrunk
    /// by `overlap_start` on the left and `overlap_end` on the right.
    ///
    /// Downstream tools that emit per-position output should restrict their
    /// emission to this range to avoid double-counting positions that
    /// neighboring tiles also cover.
    #[must_use]
    pub fn core_range(&self) -> std::ops::RangeInclusive<Pos0> {
        // overlap_start <= len - 1 (validated by SegmentOptions / iterator),
        // so checked_add_offset succeeds.
        let core_start = self
            .start
            .checked_add_offset(seqair_types::Offset::new(i64::from(self.overlap_start)))
            .unwrap_or(self.end);
        let core_end = self
            .end
            .checked_sub_offset(seqair_types::Offset::new(i64::from(self.overlap_end)))
            .unwrap_or(self.start);
        // core_start <= core_end always holds because the iterator guarantees
        // overlap_start + overlap_end < len.
        let (core_start, core_end) =
            if core_start <= core_end { (core_start, core_end) } else { (self.start, self.end) };
        core_start..=core_end
    }
}

// ── SegmentOptions ────────────────────────────────────────────────────────

// r[impl unified.segment_options]
/// Tile-size policy for [`Readers::segments`](super::Readers::segments).
///
/// There is no `Default` impl — every caller must commit to a `max_len`,
/// since the whole point of a segmenter is that there is no universal sane
/// default for "how much of a contig should one pileup load at once".
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct SegmentOptions {
    max_len: NonZeroU32,
    overlap: u32,
}

impl SegmentOptions {
    /// Plain options: tiles up to `max_len` bases, no overlap.
    ///
    /// **Picking `max_len`.** Each tile drives one BAM fetch and one FASTA
    /// fetch. A single tile occupies roughly `~40 MB` of `RecordStore`
    /// state for typical 30× coverage plus the reference bases for the
    /// tile. As a rule of thumb:
    ///
    /// * **50–500 kb** is a good default for whole-genome work — small
    ///   enough to bound peak memory, large enough that BAI bin lookup
    ///   overhead is amortized.
    /// * **5–50 kb** if you're parallelizing across many forks and want
    ///   tiles to be cheap to ship.
    /// * **>1 Mb** if you're processing very low-depth data and want to
    ///   minimize fetch round-trips. Be aware that reference bases plus
    ///   record store grow linearly with tile length.
    ///
    /// There is no universally correct value, which is why this method
    /// requires a `NonZeroU32` rather than offering a default.
    #[must_use]
    pub const fn new(max_len: NonZeroU32) -> Self {
        Self { max_len, overlap: 0 }
    }

    /// Set the per-edge overlap (bases shared with each neighboring tile).
    ///
    /// Returns `Err(SegmentOptionsError::OverlapTooLarge)` if `overlap >=
    /// max_len.get()`. Allowing that would let the iterator make zero forward
    /// progress: tiles would advance by `max_len - overlap == 0` bases.
    pub const fn with_overlap(self, overlap: u32) -> Result<Self, SegmentOptionsError> {
        if overlap >= self.max_len.get() {
            return Err(SegmentOptionsError::OverlapTooLarge {
                max_len: self.max_len.get(),
                overlap,
            });
        }
        Ok(Self { max_len: self.max_len, overlap })
    }

    /// Maximum tile length, in bases.
    #[must_use]
    pub const fn max_len(&self) -> NonZeroU32 {
        self.max_len
    }

    /// Per-edge overlap with neighboring tiles, in bases.
    #[must_use]
    pub const fn overlap(&self) -> u32 {
        self.overlap
    }
}

#[derive(Debug, thiserror::Error, PartialEq, Eq)]
#[non_exhaustive]
pub enum SegmentOptionsError {
    #[error("segment overlap {overlap} must be < max_len {max_len}")]
    OverlapTooLarge { max_len: u32, overlap: u32 },
}

// ── IntoSegmentTarget ─────────────────────────────────────────────────────

/// One concrete inclusive range on one contig, ready for tiling.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ResolvedRange {
    pub tid: Tid,
    pub contig: SmolStr,
    pub start: Pos0,
    pub end: Pos0,
    pub contig_last_pos: Pos0,
}

mod private {
    use super::{BamHeader, ReaderError, ResolvedRange};
    /// Sealed inner trait: callers can pass values whose types implement
    /// [`super::IntoSegmentTarget`], but cannot implement the trait themselves
    /// (the inner trait is in a private module).
    #[allow(
        private_interfaces,
        reason = "method signature uses crate-private ResolvedRange; trait is sealed"
    )]
    pub trait IntoSegmentTargetSealed {
        fn resolve_target(self, header: &BamHeader) -> Result<Vec<ResolvedRange>, ReaderError>;
    }
}

// r[impl unified.into_segment_target]
/// Anything that resolves into one or more contig ranges to be tiled.
///
/// Implementations are provided for every "I want a pileup over X" variant
/// that exists today: `&str`, `String`, `SmolStr`, `Tid`, `u32` (whole
/// contig); `&RegionString` and `(impl ResolveTid, Pos0, Pos0)` (explicit
/// range); `()` (whole genome). See `r[unified.into_segment_target]` for
/// the contract.
///
/// This trait is *sealed*: external crates cannot implement it. The intent
/// is that this is the closed list of targets seqair understands, not an
/// extension point.
pub trait IntoSegmentTarget: private::IntoSegmentTargetSealed {}
impl<T: private::IntoSegmentTargetSealed> IntoSegmentTarget for T {}

/// Helper used internally by [`super::Readers::segments`] to drive trait
/// dispatch through the sealed inner trait.
pub(crate) fn resolve_target<T: IntoSegmentTarget>(
    target: T,
    header: &BamHeader,
) -> Result<Vec<ResolvedRange>, ReaderError> {
    private::IntoSegmentTargetSealed::resolve_target(target, header)
}

/// Cap a `usize` target count to `u32`, returning a typed error if the cap
/// would silently lose information. BAM allows up to `i32::MAX` references,
/// so this overflow is theoretical, but we surface it as a real error rather
/// than saturating to `u32::MAX`.
fn target_count_u32(header: &BamHeader) -> Result<u32, ReaderError> {
    let count = header.target_count();
    u32::try_from(count).map_err(|_| ReaderError::HeaderTargetCountOverflow { count })
}

/// Helper: build a `ResolvedRange` for the contig identified by `tid`,
/// reusing `name` instead of looking it up in the header. Used when the
/// caller already has the contig name as a [`SmolStr`] (e.g. `&SmolStr`,
/// `&RegionString`) so we skip a redundant `target_name` + `String::to_owned`
/// + `Into<SmolStr>` round-trip.
fn whole_contig_with_name(
    header: &BamHeader,
    tid: Tid,
    name: SmolStr,
) -> Result<ResolvedRange, ReaderError> {
    let len = header.target_len(tid.as_u32()).unwrap_or(0);
    if len == 0 {
        return Err(ReaderError::EmptyContig { name });
    }
    let last = len.checked_sub(1).expect("len > 0 checked above");
    let last_u32 = u32::try_from(last).map_err(|_| ReaderError::RegionEndTooLarge { end: last })?;
    let contig_last_pos =
        Pos0::new(last_u32).ok_or(ReaderError::RegionEndTooLarge { end: last })?;
    Ok(ResolvedRange {
        tid,
        contig: name,
        start: Pos0::ZERO,
        end: contig_last_pos,
        contig_last_pos,
    })
}

/// Helper: build a `ResolvedRange` covering the entire contig identified by
/// `tid`. Errors if the contig has zero length, the tid is out of range, or
/// the header advertises more targets than fit in a `u32`.
fn whole_contig(header: &BamHeader, tid: Tid) -> Result<ResolvedRange, ReaderError> {
    let n_targets = target_count_u32(header)?;
    let name = header
        .target_name(tid.as_u32())
        .ok_or(super::resolve::TidError::TidOutOfRange { tid: tid.as_u32(), n_targets })?
        .to_owned();
    whole_contig_with_name(header, tid, name.into())
}

impl private::IntoSegmentTargetSealed for Tid {
    fn resolve_target(self, header: &BamHeader) -> Result<Vec<ResolvedRange>, ReaderError> {
        Ok(vec![whole_contig(header, self)?])
    }
}

impl private::IntoSegmentTargetSealed for u32 {
    fn resolve_target(self, header: &BamHeader) -> Result<Vec<ResolvedRange>, ReaderError> {
        let tid = self.resolve_tid(header)?;
        Ok(vec![whole_contig(header, tid)?])
    }
}

impl private::IntoSegmentTargetSealed for &str {
    fn resolve_target(self, header: &BamHeader) -> Result<Vec<ResolvedRange>, ReaderError> {
        let tid = self.resolve_tid(header)?;
        Ok(vec![whole_contig(header, tid)?])
    }
}

// Ergonomic forwarder: `&String` and `&SmolStr` exist because the trait is
// sealed, so we can't take `impl AsRef<str>` and let users coerce. They both
// route to a `whole_contig` path; `&SmolStr` short-circuits the header
// `target_name` lookup since the caller already has the name owned.
impl private::IntoSegmentTargetSealed for &String {
    fn resolve_target(self, header: &BamHeader) -> Result<Vec<ResolvedRange>, ReaderError> {
        private::IntoSegmentTargetSealed::resolve_target(self.as_str(), header)
    }
}

impl private::IntoSegmentTargetSealed for &SmolStr {
    fn resolve_target(self, header: &BamHeader) -> Result<Vec<ResolvedRange>, ReaderError> {
        let tid = self.as_str().resolve_tid(header)?;
        Ok(vec![whole_contig_with_name(header, tid, self.clone())?])
    }
}

impl private::IntoSegmentTargetSealed for &RegionString {
    fn resolve_target(self, header: &BamHeader) -> Result<Vec<ResolvedRange>, ReaderError> {
        let tid = self.chromosome.as_str().resolve_tid(header)?;
        // Reuse the parsed chromosome SmolStr instead of re-fetching it from
        // the header.
        let mut range = whole_contig_with_name(header, tid, self.chromosome.clone())?;
        if let Some(p1) = self.start {
            range.start = p1.to_zero_based();
        }
        if let Some(p1) = self.end {
            let end0 = p1.to_zero_based();
            if end0 > range.contig_last_pos {
                return Err(ReaderError::RegionEndTooLarge { end: end0.as_u64() });
            }
            range.end = end0;
        }
        if range.start > range.end {
            return Err(ReaderError::RegionStartAfterEnd {
                contig: range.contig.clone(),
                start: range.start.as_u64(),
                end: range.end.as_u64(),
            });
        }
        Ok(vec![range])
    }
}

impl<R: ResolveTid> private::IntoSegmentTargetSealed for (R, Pos0, Pos0) {
    fn resolve_target(self, header: &BamHeader) -> Result<Vec<ResolvedRange>, ReaderError> {
        let (resolver, start, end) = self;
        let tid = resolver.resolve_tid(header)?;
        let mut range = whole_contig(header, tid)?;
        if start > end {
            return Err(ReaderError::RegionStartAfterEnd {
                contig: range.contig.clone(),
                start: start.as_u64(),
                end: end.as_u64(),
            });
        }
        if end > range.contig_last_pos {
            return Err(ReaderError::RegionEndTooLarge { end: end.as_u64() });
        }
        range.start = start;
        range.end = end;
        Ok(vec![range])
    }
}

/// Whole-genome scan: every contig with non-zero length, in header order.
impl private::IntoSegmentTargetSealed for () {
    fn resolve_target(self, header: &BamHeader) -> Result<Vec<ResolvedRange>, ReaderError> {
        let n = target_count_u32(header)?;
        let mut out = Vec::with_capacity(header.target_count());
        for tid_u32 in 0..n {
            let tid = tid_u32.resolve_tid(header)?;
            match whole_contig(header, tid) {
                Ok(r) => out.push(r),
                Err(ReaderError::EmptyContig { .. }) => continue,
                Err(e) => return Err(e),
            }
        }
        Ok(out)
    }
}

// ── Segments iterator ─────────────────────────────────────────────────────

// r[impl unified.readers_segments]
/// Iterator that walks one or more `ResolvedRange`s and yields `Segment`s
/// per [`SegmentOptions`].
///
/// Built by [`Readers::segments`](super::Readers::segments). Encapsulates the
/// slightly fiddly tile arithmetic so callers can't get it wrong:
///
/// * The union of `core_range()` over consecutive tiles equals the input
///   range exactly (no overlap, no gap).
/// * Internal tiles are exactly `max_len` bases long; the last tile of a
///   range may be shorter.
/// * `overlap_start == 0` for the first tile of each range; `overlap_end ==
///   0` for the last tile of each range; internal tiles carry the requested
///   overlap on both edges.
pub struct Segments {
    ranges: std::vec::IntoIter<ResolvedRange>,
    current: Option<ResolvedRange>,
    /// 0-based inclusive start of the *next core* to emit within the current
    /// range. Advanced by `max_len - 2*overlap` after each emit (less for the
    /// first tile of a range, where there is no left overlap to subtract).
    next_core_start: Pos0,
    opts: SegmentOptions,
}

impl Segments {
    pub(crate) fn new(ranges: Vec<ResolvedRange>, opts: SegmentOptions) -> Self {
        let mut iter = ranges.into_iter();
        let current = iter.next();
        let next_core_start = current.as_ref().map_or(Pos0::ZERO, |r| r.start);
        Self { ranges: iter, current, next_core_start, opts }
    }

    fn advance_to_next_range(&mut self) {
        self.current = self.ranges.next();
        if let Some(r) = self.current.as_ref() {
            self.next_core_start = r.start;
        }
    }
}

impl Iterator for Segments {
    type Item = Segment;

    fn next(&mut self) -> Option<Segment> {
        loop {
            let range = self.current.as_ref()?;
            // Skip exhausted ranges.
            if self.next_core_start > range.end {
                self.advance_to_next_range();
                continue;
            }

            let max_len = self.opts.max_len.get();
            let overlap = self.opts.overlap;

            let range_start_u64 = range.start.as_u64();
            let range_end_u64 = range.end.as_u64();
            let core_start_u64 = self.next_core_start.as_u64();
            // Cores partition [range.start, range.end] in `max_len`-sized
            // chunks. The last chunk may be shorter if the range doesn't
            // divide evenly. core_end is inclusive.
            let core_end_u64 = core_start_u64
                .saturating_add(u64::from(max_len))
                .saturating_sub(1)
                .min(range_end_u64);
            // Tile = core ± overlap, clamped to [range.start, range.end].
            // The `overlap_start`/`overlap_end` fields record the *actual*
            // amount of overlap after clamping (0 at range edges).
            let tile_start_u64 =
                core_start_u64.saturating_sub(u64::from(overlap)).max(range_start_u64);
            let tile_end_u64 = core_end_u64.saturating_add(u64::from(overlap)).min(range_end_u64);
            let overlap_start =
                u32::try_from(core_start_u64.saturating_sub(tile_start_u64)).unwrap_or(u32::MAX);
            let overlap_end =
                u32::try_from(tile_end_u64.saturating_sub(core_end_u64)).unwrap_or(u32::MAX);

            let tile_start = pos0_from_u64(tile_start_u64).unwrap_or(range.start);
            let tile_end = pos0_from_u64(tile_end_u64).unwrap_or(range.end);

            // The iterator's arithmetic is verified by unit + property tests
            // (see `segments_match_oracle`, `cores_tile_input_exactly`), so any
            // construction failure here is an internal bug, not user input.
            let seg = Segment::new(
                range.tid,
                range.contig.clone(),
                tile_start,
                tile_end,
                overlap_start,
                overlap_end,
                range.contig_last_pos,
            )
            .expect("iterator arithmetic produces valid Segment invariants");

            // Advance to the next range explicitly when this tile reached the
            // last base of the current range. We can't rely on the
            // `next_core_start > range.end` check at the top of the loop
            // because `range.end == Pos0::max_value()` would saturate
            // `core_end + 1` back to `Pos0::max_value()` (Pos0 caps at
            // i32::MAX), keeping the loop alive forever.
            if core_end_u64 >= range_end_u64 {
                self.advance_to_next_range();
            } else {
                // In this branch `core_end_u64 < range_end_u64 <= i32::MAX`,
                // so `core_end_u64 + 1` always fits in `Pos0`. The `None`
                // arm is unreachable today; if a future change ever breaks
                // that invariant we end the iterator cleanly instead of
                // looping forever or panicking.
                let next = core_end_u64.saturating_add(1);
                match pos0_from_u64(next) {
                    Some(p) => self.next_core_start = p,
                    None => {
                        self.current = None;
                        self.ranges = Vec::new().into_iter();
                    }
                }
            }

            return Some(seg);
        }
    }
}

fn pos0_from_u64(v: u64) -> Option<Pos0> {
    let v32 = u32::try_from(v).ok()?;
    Pos0::new(v32)
}

#[cfg(test)]
#[allow(
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    clippy::format_push_string,
    clippy::comparison_chain,
    clippy::comparison_to_empty,
    reason = "test code with bounded values"
)]
mod tests {
    use super::private::IntoSegmentTargetSealed as _;
    use super::*;
    use proptest::prelude::*;

    fn p(n: u32) -> Pos0 {
        Pos0::new(n).expect("test position")
    }

    fn tid(n: u32) -> Tid {
        // Build a Tid via a tiny header so we don't reach into private constructors.
        // Round-trip through ResolveTid.
        let mut header_text = String::from("@HD\tVN:1.6\n");
        for i in 0..=n {
            header_text.push_str(&format!("@SQ\tSN:c{i}\tLN:1000\n"));
        }
        let header = crate::bam::BamHeader::from_sam_text(&header_text).unwrap();
        n.resolve_tid(&header).unwrap()
    }

    fn header_with_contigs(contigs: &[(&str, u32)]) -> crate::bam::BamHeader {
        let mut text = String::from("@HD\tVN:1.6\n");
        for (name, len) in contigs {
            text.push_str(&format!("@SQ\tSN:{name}\tLN:{len}\n"));
        }
        crate::bam::BamHeader::from_sam_text(&text).unwrap()
    }

    /// Independent tile oracle: given an inclusive range on a contig, produce
    /// the expected sequence of `(tile_start, tile_end, overlap_start, overlap_end)`
    /// using *different* arithmetic from the iterator (cumulative core ranges,
    /// then expand). This gives us a real second source of truth.
    fn oracle_tiles(
        range_start: u32,
        range_end: u32,
        max_len: u32,
        overlap: u32,
    ) -> Vec<(u32, u32, u32, u32)> {
        assert!(range_start <= range_end);
        assert!(max_len >= 1);
        assert!(overlap < max_len);

        // Step 1: produce non-overlapping cores covering [range_start, range_end].
        // Each core has length `max_len`, except the last which may be shorter.
        let mut cores: Vec<(u32, u32)> = Vec::new();
        let mut s = range_start;
        loop {
            let len_minus_one = max_len - 1;
            let e = s.saturating_add(len_minus_one).min(range_end);
            cores.push((s, e));
            if e == range_end {
                break;
            }
            s = e + 1;
        }

        // Step 2: expand each core by `overlap` on each side, then clip to
        // [range_start, range_end]. The actual overlap fields are the
        // *measured* extension on each side after clipping (not the
        // requested `overlap`), so cores at the very start or end of the
        // range get a smaller `overlap_start` / `overlap_end`.
        let mut out = Vec::with_capacity(cores.len());
        for (cs, ce) in cores {
            let want_left = u64::from(cs).saturating_sub(u64::from(overlap));
            let want_right = u64::from(ce).saturating_add(u64::from(overlap));
            let tile_start =
                u32::try_from(want_left.max(u64::from(range_start))).expect("range fits in u32");
            let tile_end =
                u32::try_from(want_right.min(u64::from(range_end))).expect("range fits in u32");
            let overlap_start = cs - tile_start;
            let overlap_end = tile_end - ce;
            out.push((tile_start, tile_end, overlap_start, overlap_end));
        }
        out
    }

    // r[verify unified.segment_options]
    #[test]
    fn segment_options_rejects_too_large_overlap() {
        let opts = SegmentOptions::new(NonZeroU32::new(100).unwrap());
        assert!(opts.with_overlap(99).is_ok());
        assert!(matches!(
            opts.with_overlap(100),
            Err(SegmentOptionsError::OverlapTooLarge { max_len: 100, overlap: 100 })
        ));
        assert!(matches!(opts.with_overlap(101), Err(SegmentOptionsError::OverlapTooLarge { .. })));
    }

    // r[verify unified.segment_struct]
    // r[verify unified.segment_overlap]
    #[test]
    fn segment_core_range_excludes_overlap() {
        let seg = Segment::new(tid(0), "c0".into(), p(100), p(199), 10, 5, p(999)).unwrap();
        assert_eq!(seg.len(), 100);
        let core = seg.core_range();
        assert_eq!(*core.start(), p(110));
        assert_eq!(*core.end(), p(194));
    }

    // r[verify unified.segment_struct]
    #[test]
    fn segment_new_rejects_start_after_end() {
        let err =
            Segment::new(tid(0), "c0".into(), p(200), p(100), 0, 0, p(999)).expect_err("invalid");
        assert!(matches!(err, SegmentInvariantError::StartAfterEnd { start: 200, end: 100 }));
    }

    // r[verify unified.segment_struct]
    #[test]
    fn segment_new_rejects_end_past_contig() {
        let err =
            Segment::new(tid(0), "c0".into(), p(0), p(500), 0, 0, p(100)).expect_err("invalid");
        assert!(matches!(
            err,
            SegmentInvariantError::EndPastContig { end: 500, contig_last_pos: 100 }
        ));
    }

    // r[verify unified.segment_struct]
    // r[verify unified.segment_overlap]
    #[test]
    fn segment_new_rejects_overlap_total_ge_len() {
        // Tile length is 100. overlap_start + overlap_end == 100 covers the
        // entire tile and leaves zero core, which would break neighbor-dedupe.
        let err =
            Segment::new(tid(0), "c0".into(), p(0), p(99), 50, 50, p(999)).expect_err("invalid");
        assert!(matches!(
            err,
            SegmentInvariantError::OverlapTooLarge { overlap_start: 50, overlap_end: 50, len: 100 }
        ));
        // overlap_start + overlap_end > len is also rejected.
        let err2 =
            Segment::new(tid(0), "c0".into(), p(0), p(99), 60, 50, p(999)).expect_err("invalid");
        assert!(matches!(err2, SegmentInvariantError::OverlapTooLarge { .. }));
        // Boundary: total == len - 1 is allowed (leaves 1 base of core).
        let ok = Segment::new(tid(0), "c0".into(), p(0), p(99), 50, 49, p(999)).unwrap();
        assert_eq!(ok.len(), 100);
        assert_eq!(*ok.core_range().start(), p(50));
        assert_eq!(*ok.core_range().end(), p(50));
    }

    // r[verify unified.segment_struct]
    /// `Segment::len()` is computed once at construction. Spot-check the
    /// stored value across a small spread of tile shapes.
    #[test]
    fn segment_len_matches_inclusive_span() {
        for (s, e, expected) in
            [(0u32, 0u32, 1u32), (0, 99, 100), (5, 14, 10), (1_000_000, 1_000_499, 500)]
        {
            let seg = Segment::new(tid(0), "c0".into(), p(s), p(e), 0, 0, p(2_000_000)).unwrap();
            assert_eq!(seg.len(), expected, "span [{s}..={e}]");
        }
    }

    #[test]
    fn segments_single_tile_when_range_fits() {
        let header = header_with_contigs(&[("chr1", 500)]);
        let opts = SegmentOptions::new(NonZeroU32::new(1000).unwrap());
        let ranges = ("chr1", p(0), p(99)).resolve_target(&header).unwrap();
        let segs: Vec<_> = Segments::new(ranges, opts).collect();
        assert_eq!(segs.len(), 1);
        let s = &segs[0];
        assert_eq!(s.start(), p(0));
        assert_eq!(s.end(), p(99));
        assert_eq!(s.overlap_start(), 0);
        assert_eq!(s.overlap_end(), 0);
        assert!(s.starts_at_contig_start());
        assert_eq!(s.contig_last_pos(), p(499));
    }

    #[test]
    fn segments_many_tiles_no_overlap() {
        let header = header_with_contigs(&[("chr1", 1_000_000)]);
        let opts = SegmentOptions::new(NonZeroU32::new(100).unwrap());
        let ranges = ("chr1", p(0), p(249)).resolve_target(&header).unwrap();
        let segs: Vec<_> = Segments::new(ranges, opts).collect();
        // Cores: [0..99], [100..199], [200..249] → 3 tiles.
        assert_eq!(segs.len(), 3);
        assert_eq!((segs[0].start(), segs[0].end()), (p(0), p(99)));
        assert_eq!((segs[1].start(), segs[1].end()), (p(100), p(199)));
        assert_eq!((segs[2].start(), segs[2].end()), (p(200), p(249)));
        // No overlaps requested.
        for s in &segs {
            assert_eq!(s.overlap_start(), 0);
            assert_eq!(s.overlap_end(), 0);
        }
    }

    #[test]
    fn segments_with_overlap() {
        let header = header_with_contigs(&[("chr1", 1000)]);
        let opts = SegmentOptions::new(NonZeroU32::new(100).unwrap()).with_overlap(10).unwrap();
        let ranges = ("chr1", p(0), p(249)).resolve_target(&header).unwrap();
        let segs: Vec<_> = Segments::new(ranges, opts).collect();

        // Cores: [0..99], [100..199], [200..249]
        // First tile: [0..109]   (overlap_start=0, overlap_end=10)
        // Mid tile:   [90..209]  (overlap_start=10, overlap_end=10)
        // Last tile:  [190..249] (overlap_start=10, overlap_end=0)
        assert_eq!(segs.len(), 3);
        assert_eq!(
            (segs[0].start(), segs[0].end(), segs[0].overlap_start(), segs[0].overlap_end()),
            (p(0), p(109), 0, 10)
        );
        assert_eq!(
            (segs[1].start(), segs[1].end(), segs[1].overlap_start(), segs[1].overlap_end()),
            (p(90), p(209), 10, 10)
        );
        assert_eq!(
            (segs[2].start(), segs[2].end(), segs[2].overlap_start(), segs[2].overlap_end()),
            (p(190), p(249), 10, 0)
        );

        // core_range() of all tiles must tile [0..249] exactly.
        let cores: Vec<_> = segs.iter().map(|s| s.core_range()).collect();
        assert_eq!(*cores[0].start(), p(0));
        assert_eq!(*cores[0].end(), p(99));
        assert_eq!(*cores[1].start(), p(100));
        assert_eq!(*cores[1].end(), p(199));
        assert_eq!(*cores[2].start(), p(200));
        assert_eq!(*cores[2].end(), p(249));
    }

    // r[verify unified.readers_segments]
    /// Regression: an iterator over a `ResolvedRange` whose `end ==
    /// Pos0::max_value()` must terminate. Before the fix, advancing past
    /// the last tile saturated `core_end + 1` back to `Pos0::max_value()`,
    /// failed the `next_core_start > range.end` check, and re-emitted the
    /// same single-base tile forever.
    ///
    /// The public `IntoSegmentTarget` impls cap `end` at `contig_last_pos
    /// <= i32::MAX - 1`, so this scenario is currently unreachable from
    /// outside the crate. The defence is still load-bearing: a future
    /// `IntoSegmentTarget` impl, an internal caller, or a fuzz harness
    /// could construct a `ResolvedRange` with `end == Pos0::max_value()`,
    /// and the iterator must not hang.
    #[test]
    fn segments_terminate_at_pos0_max() {
        let last = Pos0::max_value();
        // Build a ResolvedRange directly that ends at Pos0::max_value() and
        // is short enough that we expect exactly one tile.
        let near_max = Pos0::new(u32::try_from(last.as_u64()).unwrap() - 5).unwrap();
        let mut header_text = String::from("@HD\tVN:1.6\n");
        header_text.push_str("@SQ\tSN:fake\tLN:1000\n");
        let header = crate::bam::BamHeader::from_sam_text(&header_text).unwrap();
        let tid = 0u32.resolve_tid(&header).unwrap();
        let range = ResolvedRange {
            tid,
            contig: "fake".into(),
            start: near_max,
            end: last,
            contig_last_pos: last,
        };
        let opts = SegmentOptions::new(NonZeroU32::new(1_000_000).unwrap());
        let segs: Vec<_> = Segments::new(vec![range], opts).collect();
        assert_eq!(segs.len(), 1);
        assert_eq!(segs[0].end(), last);
        assert_eq!(segs[0].start(), near_max);
    }

    #[test]
    fn segments_whole_genome_skips_empty_contigs() {
        // Build header by hand so we can include a zero-length contig.
        let text = "@HD\tVN:1.6\n@SQ\tSN:a\tLN:50\n@SQ\tSN:empty\tLN:0\n@SQ\tSN:c\tLN:30\n";
        let header = crate::bam::BamHeader::from_sam_text(text).unwrap();
        let opts = SegmentOptions::new(NonZeroU32::new(100).unwrap());
        let ranges = ().resolve_target(&header).unwrap();
        let segs: Vec<_> = Segments::new(ranges, opts).collect();
        assert_eq!(segs.len(), 2);
        assert_eq!(segs[0].contig().as_str(), "a");
        assert_eq!(segs[1].contig().as_str(), "c");
    }

    // r[verify unified.into_segment_target]
    #[test]
    fn into_segment_target_empty_contig_errors_for_named() {
        let text = "@HD\tVN:1.6\n@SQ\tSN:empty\tLN:0\n";
        let header = crate::bam::BamHeader::from_sam_text(text).unwrap();
        let err = "empty".resolve_target(&header).unwrap_err();
        assert!(matches!(err, ReaderError::EmptyContig { .. }));
    }

    // r[verify unified.into_segment_target]
    #[test]
    fn into_segment_target_start_after_end_errors() {
        let header = header_with_contigs(&[("chr1", 1000)]);
        let err = ("chr1", p(100), p(50)).resolve_target(&header).unwrap_err();
        assert!(matches!(err, ReaderError::RegionStartAfterEnd { .. }));
    }

    // r[verify unified.into_segment_target]
    #[test]
    fn into_segment_target_end_past_contig_errors() {
        let header = header_with_contigs(&[("chr1", 100)]);
        let err = ("chr1", p(0), p(500)).resolve_target(&header).unwrap_err();
        assert!(matches!(err, ReaderError::RegionEndTooLarge { .. }));
    }

    proptest! {
        // r[verify unified.readers_segments]
        #[test]
        fn segments_match_oracle(
            range_start in 0u32..1_000_000,
            len in 1u32..1_000_000,
            max_len in 1u32..2_000,
            overlap in 0u32..2_000,
        ) {
            // Constrain combinations the iterator's contract requires.
            prop_assume!(overlap < max_len);
            let range_end = range_start.saturating_add(len - 1).min(1_999_999);

            let opts = SegmentOptions::new(NonZeroU32::new(max_len).unwrap())
                .with_overlap(overlap).unwrap();

            let header = header_with_contigs(&[("chr1", 2_000_000)]);
            let ranges = ("chr1", p(range_start), p(range_end)).resolve_target(&header).unwrap();
            let actual: Vec<(u32, u32, u32, u32)> = Segments::new(ranges, opts)
                .map(|s| (
                    u32::try_from(s.start().as_u64()).unwrap(),
                    u32::try_from(s.end().as_u64()).unwrap(),
                    s.overlap_start(),
                    s.overlap_end(),
                ))
                .collect();

            let expected = oracle_tiles(range_start, range_end, max_len, overlap);
            prop_assert_eq!(actual, expected);
        }

        // r[verify unified.segment_overlap]
        #[test]
        fn cores_tile_input_exactly(
            range_start in 0u32..100_000,
            len in 1u32..100_000,
            max_len in 1u32..500,
            overlap in 0u32..500,
        ) {
            prop_assume!(overlap < max_len);
            let range_end = range_start.saturating_add(len - 1).min(199_999);

            let opts = SegmentOptions::new(NonZeroU32::new(max_len).unwrap())
                .with_overlap(overlap).unwrap();
            let header = header_with_contigs(&[("chr1", 200_000)]);
            let ranges = ("chr1", p(range_start), p(range_end)).resolve_target(&header).unwrap();
            let segs: Vec<_> = Segments::new(ranges, opts).collect();

            // Each core_range must be contiguous with the next; first must
            // start at range_start; last must end at range_end.
            prop_assert!(!segs.is_empty());
            let first = &segs[0];
            let last = segs.last().unwrap();
            prop_assert_eq!(*first.core_range().start(), p(range_start));
            prop_assert_eq!(*last.core_range().end(), p(range_end));
            for w in segs.windows(2) {
                let a_end = w[0].core_range().end().as_u64();
                let b_start = w[1].core_range().start().as_u64();
                prop_assert_eq!(a_end + 1, b_start, "core ranges must be contiguous");
            }
            // No tile exceeds max_len; no tile is empty.
            for s in &segs {
                prop_assert!(!s.is_empty());
                prop_assert!(s.len() <= max_len + 2 * overlap);
            }
        }

        // r[verify unified.segment_overlap]
        #[test]
        fn first_and_last_overlaps_are_zero(
            // Vary range_start so the test isn't asserting a property about
            // contig boundaries; the range edge is what counts. The original
            // test pinned range_start=0 and built `total = n_tiles*max_len`,
            // making the inner-tile claim "overlap == overlap" tautological.
            range_start in 0u32..100_000,
            // `extra_len` is added to a deliberate non-multiple of max_len so
            // the last core may be short and the last tile must clip its
            // overlap_end to range_end without producing the requested overlap.
            n_full_cores in 1u32..20,
            extra_len in 0u32..200,
            max_len in 10u32..200,
            overlap in 0u32..9,
        ) {
            prop_assume!(overlap < max_len);
            // total bases = n_full_cores * max_len + extra_len, then -1 for inclusive.
            let total = u64::from(n_full_cores) * u64::from(max_len) + u64::from(extra_len);
            // Build the range somewhere in the middle of the contig.
            let range_end_u64 = u64::from(range_start).saturating_add(total).saturating_sub(1);
            let contig_len = range_end_u64.saturating_add(10);
            let contig_len_u32 = u32::try_from(contig_len).unwrap_or(u32::MAX);
            let range_end = u32::try_from(range_end_u64).unwrap_or(u32::MAX / 2);

            let opts = SegmentOptions::new(NonZeroU32::new(max_len).unwrap())
                .with_overlap(overlap).unwrap();
            let header = header_with_contigs(&[("chr1", contig_len_u32)]);
            let ranges =
                ("chr1", p(range_start), p(range_end)).resolve_target(&header).unwrap();
            let segs: Vec<_> = Segments::new(ranges, opts).collect();
            prop_assert!(!segs.is_empty());
            // Range-edge invariants: first segment of the *requested range*
            // has overlap_start == 0; last segment has overlap_end == 0.
            prop_assert_eq!(segs.first().unwrap().overlap_start(), 0);
            prop_assert_eq!(segs.last().unwrap().overlap_end(), 0);
            // First segment's tile_start equals the requested range_start;
            // last segment's tile_end equals the requested range_end. This
            // is a non-trivial check now that range_start != 0.
            prop_assert_eq!(segs.first().unwrap().start(), p(range_start));
            prop_assert_eq!(segs.last().unwrap().end(), p(range_end));
            // Every tile sits entirely inside the requested range.
            for s in &segs {
                prop_assert!(s.start() >= p(range_start));
                prop_assert!(s.end() <= p(range_end));
            }
        }
    }
}
