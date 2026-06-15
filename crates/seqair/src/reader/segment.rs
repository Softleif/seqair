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
use std::num::{NonZeroU32, NonZeroU64};

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

// `Segment::len()` always returns `>= 1` (the constructor rejects empty
// tiles) so a paired `is_empty()` would be a constant-`false` method. Skip
// the noise rather than expose it.
#[allow(
    clippy::len_without_is_empty,
    reason = "Segment is never empty by construction; an is_empty() method would be noise"
)]
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
        // The constructor enforces `overlap_start + overlap_end < len` and
        // `end <= i32::MAX` (Pos0 invariant), so:
        //   * `start + overlap_start <= end - overlap_end <= i32::MAX`
        //     (so both checked_*_offset calls succeed); and
        //   * `core_start <= core_end` (so the resulting range is
        //     non-empty).
        //
        // Using `expect` rather than a silent fallback means that if a
        // future change ever weakens those invariants the failure surfaces
        // as a panic instead of degrading silently to "the entire tile is
        // its own core" — which would double-count overlap regions in
        // every downstream caller.
        let core_start = self
            .start
            .checked_add_offset(seqair_types::Offset::new(i64::from(self.overlap_start)))
            .expect("constructor invariant: start + overlap_start <= end <= i32::MAX");
        let core_end = self
            .end
            .checked_sub_offset(seqair_types::Offset::new(i64::from(self.overlap_end)))
            .expect("constructor invariant: end - overlap_end >= start >= 0");
        debug_assert!(
            core_start <= core_end,
            "constructor invariant: overlap_start + overlap_end < len => core_start <= core_end"
        );
        core_start..=core_end
    }
}

// ── SegmentOptions ────────────────────────────────────────────────────────

// r[impl unified.segment_options]
/// Tile-size policy for [`Readers::segments`](super::Readers::segments).
///
/// `Default` returns 10 kb tiles with no overlap — a conservative choice
/// that keeps peak `RecordStore` memory low and works for most exploratory
/// or per-region use. For whole-genome scans, larger tiles (50–500 kb,
/// see [`SegmentOptions::new`]) amortise BAI/CRAI lookup overhead better.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct SegmentOptions {
    max_len: NonZeroU32,
    overlap: u32,
    max_bytes: Option<NonZeroU64>,
}

/// Default per-segment compressed-byte budget (256 MiB), matching the
/// `RegionBuf` bulk-load guard. The budget is in **compressed** bytes (what
/// the index can estimate); the decoded `RecordStore` is a few× larger.
const DEFAULT_MAX_BYTES: u64 = 256 * 1024 * 1024;

impl Default for SegmentOptions {
    /// 10 kb tiles, no overlap, 256 MiB byte budget. Suits exploratory and
    /// per-region work; pick a larger `max_len` via [`SegmentOptions::new`]
    /// for whole-genome scans.
    fn default() -> Self {
        // SAFETY: 10_000 is non-zero.
        let max_len = NonZeroU32::new(10_000).expect("non-zero literal");
        Self { max_len, overlap: 0, max_bytes: NonZeroU64::new(DEFAULT_MAX_BYTES) }
    }
}

impl SegmentOptions {
    /// Plain options: tile *cores* up to `max_len` bases, no overlap, with the
    /// default 256 MiB byte budget (see [`with_max_bytes`](Self::with_max_bytes)).
    ///
    /// `max_len` caps the **core** length of each tile
    /// (`Segment::core_range()`). With
    /// [`with_overlap(o)`](Self::with_overlap), an internal tile's full
    /// `[start, end]` is the core expanded by `o` bases on each side, so
    /// internal tiles can be up to `max_len + 2 * o` bases. Edge tiles
    /// clip their overlap to the requested range and are smaller.
    #[must_use]
    pub const fn new(max_len: NonZeroU32) -> Self {
        Self { max_len, overlap: 0, max_bytes: NonZeroU64::new(DEFAULT_MAX_BYTES) }
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
        Ok(Self { max_len: self.max_len, overlap, max_bytes: self.max_bytes })
    }

    /// Cap each segment's estimated **compressed** load to `max_bytes`.
    ///
    /// [`Readers::segments`](super::Readers::segments) consults the index and
    /// subdivides any tile whose estimated byte size exceeds this budget, so
    /// each emitted [`Segment`] stays loadable within it. A single index leaf
    /// bin (~16 kb) that alone exceeds the budget can't be split further and
    /// is emitted as-is.
    #[must_use]
    pub const fn with_max_bytes(self, max_bytes: NonZeroU64) -> Self {
        Self { max_len: self.max_len, overlap: self.overlap, max_bytes: Some(max_bytes) }
    }

    /// Disable the byte budget: segments are bounded by `max_len` (positions)
    /// only, regardless of how much data they cover.
    #[must_use]
    pub const fn without_byte_budget(self) -> Self {
        Self { max_len: self.max_len, overlap: self.overlap, max_bytes: None }
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

    /// Per-segment compressed-byte budget, or `None` if disabled.
    #[must_use]
    pub const fn max_bytes(&self) -> Option<NonZeroU64> {
        self.max_bytes
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
///
/// `Tid` carries no link back to the header it was resolved against (it's
/// a passthrough impl of `ResolveTid`), so we re-validate the tid against
/// the current header here. A stale tid surfaces as `TidOutOfRange` rather
/// than masquerading as `EmptyContig`.
fn whole_contig_with_name(
    header: &BamHeader,
    tid: Tid,
    name: SmolStr,
) -> Result<ResolvedRange, ReaderError> {
    let Some(len) = header.target_len(tid.as_u32()) else {
        let n_targets = target_count_u32(header)?;
        return Err(super::resolve::TidError::TidOutOfRange { tid: tid.as_u32(), n_targets }.into());
    };
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
    // Compute `n_targets` only when needed for the error path. The hot path
    // (tid in range) avoids a redundant `u32::try_from(usize)` per contig.
    let name = match header.target_name(tid.as_u32()) {
        Some(n) => n.to_owned(),
        None => {
            let n_targets = target_count_u32(header)?;
            return Err(
                super::resolve::TidError::TidOutOfRange { tid: tid.as_u32(), n_targets }.into()
            );
        }
    };
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
                return Err(ReaderError::RegionEndPastContig {
                    contig: range.contig.clone(),
                    end: end0.as_u64(),
                    contig_last_pos: range.contig_last_pos.as_u64(),
                });
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
            return Err(ReaderError::RegionEndPastContig {
                contig: range.contig.clone(),
                end: end.as_u64(),
                contig_last_pos: range.contig_last_pos.as_u64(),
            });
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
/// * Internal tile **cores** are exactly `max_len` bases long; the last
///   core of a range may be shorter. The full `[start, end]` of an
///   internal tile is its core expanded by `overlap` bases on each side
///   (so `len() == max_len + 2 * overlap`); first/last tiles of a range
///   clip their overlap to the requested range.
/// * `overlap_start == 0` for the first tile of each range; `overlap_end ==
///   0` for the last tile of each range; internal tiles carry the requested
///   overlap on both edges.
pub struct Segments {
    inner: SegmentsInner,
}

enum SegmentsInner {
    /// Lazy positional tiling — bounded by `max_len` only. Used directly when
    /// no byte budget applies (or `Segments` is built without index access).
    Planned(PlannedSegments),
    /// Byte-aware subdivision already applied by
    /// [`Readers::segments`](super::Readers::segments); replay the built list.
    Precomputed(std::vec::IntoIter<Segment>),
}

impl Segments {
    pub(crate) fn new(ranges: Vec<ResolvedRange>, opts: SegmentOptions) -> Self {
        Self { inner: SegmentsInner::Planned(PlannedSegments::new(ranges, opts)) }
    }

    /// Wrap a pre-built segment list (after byte-aware subdivision).
    pub(crate) fn precomputed(segments: Vec<Segment>) -> Self {
        Self { inner: SegmentsInner::Precomputed(segments.into_iter()) }
    }
}

impl Iterator for Segments {
    type Item = Segment;

    fn next(&mut self) -> Option<Segment> {
        match &mut self.inner {
            SegmentsInner::Planned(p) => p.next(),
            SegmentsInner::Precomputed(it) => it.next(),
        }
    }
}

/// Lazy positional tiler: walks `ResolvedRange`s and emits `max_len`-core
/// tiles. The byte budget in [`SegmentOptions`] is **not** applied here — that
/// happens in [`Readers::segments`](super::Readers::segments), which has index
/// access; this keeps the positional arithmetic (and its tests) index-free.
struct PlannedSegments {
    ranges: std::vec::IntoIter<ResolvedRange>,
    current: Option<ResolvedRange>,
    /// 0-based inclusive start of the *next core* to emit within the current
    /// range. Advanced by `max_len - 2*overlap` after each emit (less for the
    /// first tile of a range, where there is no left overlap to subtract).
    next_core_start: Pos0,
    opts: SegmentOptions,
}

impl PlannedSegments {
    fn new(ranges: Vec<ResolvedRange>, opts: SegmentOptions) -> Self {
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

// r[impl unified.segment_byte_budget]
/// Recursively bisect `seg` so that each emitted sub-segment's full
/// `[start, end]` is estimated by `estimate` to load at most `budget` bytes.
///
/// `overlap` is the per-edge overlap from [`SegmentOptions`]. Each sub-segment
/// is clamped to the parent segment's `[start, end]`, which reproduces the
/// positional iterator's overlap rules for free: an interior split boundary
/// gets the full `overlap` on each side, while a boundary that coincides with
/// the parent's edge inherits the parent's (already range-clamped) overlap.
///
/// A single-base core that still exceeds `budget` is emitted as-is — the index
/// can't address finer than its ~16 kb leaf bin, so callers must tolerate the
/// rare over-budget segment (it warns and loads in `RegionBuf`).
pub(crate) fn split_segment_by_bytes(
    seg: &Segment,
    overlap: u32,
    budget: u64,
    estimate: &mut dyn FnMut(u32, Pos0, Pos0) -> u64,
    out: &mut Vec<Segment>,
) {
    let ctx = SplitCtx {
        tid: seg.tid(),
        contig: seg.contig(),
        contig_last_pos: seg.contig_last_pos(),
        clamp_start: seg.start().as_u64(),
        clamp_end: seg.end().as_u64(),
        overlap: u64::from(overlap),
    };
    let core = seg.core_range();
    let total_end = core.end().as_u64();
    let mut core_start = core.start().as_u64();

    // Greedy forward growth: each segment grows its core to the largest end
    // still within budget, then the next starts after it. This is robust to the
    // index's ~16 kb quantization (bisecting at a midpoint that lands mid-bin
    // wouldn't lower the estimate), and it bounds the irreducible case: if even
    // a single position is over budget (a leaf bin bigger than the budget), the
    // estimate stays flat across that bin, so growth stops at the bin's end and
    // the bin becomes ONE over-budget segment — loaded once, not once per base.
    loop {
        let base = estimate_full(&ctx, core_start, core_start, estimate);
        let threshold = budget.max(base);
        let core_end = largest_core_end(&ctx, core_start, total_end, threshold, estimate);
        emit_segment(&ctx, core_start, core_end, out);
        if core_end >= total_end {
            break;
        }
        core_start = core_end.saturating_add(1);
    }
}

/// Binary-search the largest `core_end` in `[core_start, total_end]` whose
/// full-range estimate is `<= threshold`. Relies on the estimate being
/// monotonic non-decreasing in `core_end` (a larger range queries a superset of
/// index bins). `core_start` itself always qualifies, since the caller sets
/// `threshold >= estimate(core_start)`.
fn largest_core_end(
    ctx: &SplitCtx<'_>,
    core_start: u64,
    total_end: u64,
    threshold: u64,
    estimate: &mut dyn FnMut(u32, Pos0, Pos0) -> u64,
) -> u64 {
    let mut lo = core_start;
    let mut hi = total_end;
    while lo < hi {
        // Upper mid (in [lo+1, hi]) so the search always makes progress.
        let mid = lo.saturating_add(hi.saturating_sub(lo).saturating_add(1) / 2);
        if estimate_full(ctx, core_start, mid, estimate) <= threshold {
            lo = mid;
        } else {
            hi = mid.saturating_sub(1);
        }
    }
    lo
}

/// Invariant context for the byte-split helpers — fixed for one segment.
#[derive(Clone, Copy)]
struct SplitCtx<'a> {
    tid: Tid,
    contig: &'a SmolStr,
    contig_last_pos: Pos0,
    clamp_start: u64,
    clamp_end: u64,
    overlap: u64,
}

/// Full `[start, end]` of a core: core ± overlap, clamped to the parent extent.
fn full_range(ctx: &SplitCtx<'_>, core_start: u64, core_end: u64) -> (Pos0, Pos0) {
    let full_start = core_start.saturating_sub(ctx.overlap).max(ctx.clamp_start);
    let full_end = core_end.saturating_add(ctx.overlap).min(ctx.clamp_end);
    (
        pos0_from_u64(full_start).expect("clamped start within Pos0"),
        pos0_from_u64(full_end).expect("clamped end within Pos0"),
    )
}

fn estimate_full(
    ctx: &SplitCtx<'_>,
    core_start: u64,
    core_end: u64,
    estimate: &mut dyn FnMut(u32, Pos0, Pos0) -> u64,
) -> u64 {
    let (full_start, full_end) = full_range(ctx, core_start, core_end);
    estimate(ctx.tid.as_u32(), full_start, full_end)
}

fn emit_segment(ctx: &SplitCtx<'_>, core_start: u64, core_end: u64, out: &mut Vec<Segment>) {
    let (full_start, full_end) = full_range(ctx, core_start, core_end);
    let overlap_start =
        u32::try_from(core_start.saturating_sub(full_start.as_u64())).unwrap_or(u32::MAX);
    let overlap_end = u32::try_from(full_end.as_u64().saturating_sub(core_end)).unwrap_or(u32::MAX);
    let seg = Segment::new(
        ctx.tid,
        ctx.contig.clone(),
        full_start,
        full_end,
        overlap_start,
        overlap_end,
        ctx.contig_last_pos,
    )
    .expect("byte-split produces valid Segment invariants");
    out.push(seg);
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
    /// `whole_contig_with_name` must surface a stale (out-of-range) `Tid` as
    /// `TidOutOfRange`, not mask it as `EmptyContig`. The bug is unreachable
    /// from today's public API (the sealed `IntoSegmentTarget` impls always
    /// validate the tid by name first), but it would surface the moment a
    /// new internal caller — or a fuzz harness — passes a `Tid` resolved
    /// against a different header.
    #[test]
    fn whole_contig_with_name_rejects_stale_tid() {
        // Build header A with 5 contigs.
        let header_a = header_with_contigs(&[
            ("a0", 1000),
            ("a1", 1000),
            ("a2", 1000),
            ("a3", 1000),
            ("a4", 1000),
        ]);
        let stale_tid = 4u32.resolve_tid(&header_a).unwrap();

        // Build header B with only 2 contigs. The stale `Tid(4)` is out of
        // range against B even though it was a valid Tid against A.
        let header_b = header_with_contigs(&[("b0", 1000), ("b1", 1000)]);
        let err = whole_contig_with_name(&header_b, stale_tid, "b0".into()).unwrap_err();
        assert!(
            matches!(
                err,
                ReaderError::Tid {
                    source: super::super::resolve::TidError::TidOutOfRange { tid: 4, n_targets: 2 }
                }
            ),
            "expected TidOutOfRange, got {err:?}"
        );
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
        assert!(matches!(err, ReaderError::RegionEndPastContig { .. }));
    }

    // ── Byte-aware subdivision (split_segment_by_bytes) ──────────────────────

    /// Mock byte oracle: `bpb` compressed bytes per base of the *full* range.
    fn span_estimate(bpb: u64) -> impl FnMut(u32, Pos0, Pos0) -> u64 {
        move |_tid, s, e| (e.as_u64() - s.as_u64() + 1) * bpb
    }

    /// Mirror `Readers::segments`: positional tiling, then byte-aware split of
    /// each tile using `estimate`. Returns the final segment list.
    fn byte_aware(
        contig_len: u32,
        range: std::ops::RangeInclusive<u32>,
        opts: SegmentOptions,
        budget: u64,
        mut estimate: impl FnMut(u32, Pos0, Pos0) -> u64,
    ) -> Vec<Segment> {
        let header = header_with_contigs(&[("chr1", contig_len)]);
        let ranges = ("chr1", p(*range.start()), p(*range.end())).resolve_target(&header).unwrap();
        let mut out = Vec::new();
        for seg in Segments::new(ranges, opts) {
            split_segment_by_bytes(&seg, opts.overlap(), budget, &mut estimate, &mut out);
        }
        out
    }

    /// Assert the segments' cores partition `[start, end]` exactly (contiguous,
    /// no gap, no overlap).
    fn assert_cores_tile(segs: &[Segment], start: Pos0, end: Pos0) {
        assert!(!segs.is_empty(), "no segments produced");
        assert_eq!(*segs.first().unwrap().core_range().start(), start, "first core start");
        assert_eq!(*segs.last().unwrap().core_range().end(), end, "last core end");
        for w in segs.windows(2) {
            assert_eq!(
                w[0].core_range().end().as_u64() + 1,
                w[1].core_range().start().as_u64(),
                "cores must be contiguous"
            );
        }
    }

    // r[verify unified.segment_byte_budget]
    #[test]
    fn byte_split_passthrough_when_under_budget() {
        let opts = SegmentOptions::new(NonZeroU32::new(1000).unwrap());
        let segs = byte_aware(2000, 0..=99, opts, 1_000_000, span_estimate(10));
        assert_eq!(segs.len(), 1, "a small segment should not be split");
        assert_eq!((segs[0].start(), segs[0].end()), (p(0), p(99)));
        assert_eq!((segs[0].overlap_start(), segs[0].overlap_end()), (0, 0));
    }

    // r[verify unified.segment_byte_budget]
    #[test]
    fn byte_split_splits_until_under_budget() {
        // 1000 bases × 10 B = 10_000 B; budget 2000 B forces several splits.
        let budget = 2000;
        let opts = SegmentOptions::new(NonZeroU32::new(10_000).unwrap());
        let segs = byte_aware(2000, 0..=999, opts, budget, span_estimate(10));
        assert!(segs.len() > 1, "expected the segment to be split");
        for s in &segs {
            let bytes = (s.end().as_u64() - s.start().as_u64() + 1) * 10;
            assert!(
                bytes <= budget,
                "segment {:?} = {bytes} B exceeds {budget}",
                (s.start(), s.end())
            );
        }
        assert_cores_tile(&segs, p(0), p(999));
    }

    // r[verify unified.segment_byte_budget]
    // r[verify unified.segment_overlap]
    #[test]
    fn byte_split_preserves_overlap_and_tiling() {
        let opts = SegmentOptions::new(NonZeroU32::new(10_000).unwrap()).with_overlap(10).unwrap();
        let segs = byte_aware(5000, 0..=999, opts, 2000, span_estimate(10));
        assert!(segs.len() > 2, "need interior segments to check overlap");
        // Range edges have zero outward overlap; interior boundaries carry the
        // requested overlap on both sides.
        assert_eq!(segs.first().unwrap().overlap_start(), 0);
        assert_eq!(segs.last().unwrap().overlap_end(), 0);
        for s in &segs[1..segs.len() - 1] {
            assert_eq!(s.overlap_start(), 10, "interior overlap_start");
            assert_eq!(s.overlap_end(), 10, "interior overlap_end");
        }
        assert_cores_tile(&segs, p(0), p(999));
    }

    // r[verify unified.segment_byte_budget]
    /// The irreducible case: when bisecting never lowers the byte estimate (a
    /// single oversized index leaf bin reports the same size for any sub-range),
    /// the region must be emitted as ONE over-budget segment — not bisected into
    /// many pieces that would each reload the whole bin.
    #[test]
    fn byte_split_irreducible_region_emits_one_segment() {
        let opts = SegmentOptions::new(NonZeroU32::new(10_000).unwrap());
        // Constant estimate for every range → splitting can never help.
        let segs = byte_aware(2000, 100..=109, opts, 1, |_t, _s, _e| 1_000_000);
        assert_eq!(segs.len(), 1, "irreducible region must not be bisected");
        assert_eq!((segs[0].start(), segs[0].end()), (p(100), p(109)));
        assert_cores_tile(&segs, p(100), p(109));
    }

    // r[verify unified.segment_byte_budget]
    /// A reducible region containing one irreducible "leaf bin" must split the
    /// reducible part but keep the bin in a single segment (loaded once), never
    /// reloading it per sub-position.
    #[test]
    fn byte_split_irreducible_bin_loaded_once() {
        let opts = SegmentOptions::new(NonZeroU32::new(10_000).unwrap());
        // A 100-base "leaf bin" [200,299] reports 1 MB for ANY overlapping query;
        // everywhere else costs 10 B/base.
        let estimate = |_t: u32, s: Pos0, e: Pos0| {
            let (s, e) = (s.as_u64(), e.as_u64());
            if s <= 299 && e >= 200 { 1_000_000 } else { (e - s + 1) * 10 }
        };
        let segs = byte_aware(2000, 0..=399, opts, 5_000, estimate);

        assert_cores_tile(&segs, p(0), p(399));
        // Far fewer than the ~400 single-base pieces a naive bisect would make.
        assert!(segs.len() < 20, "got {} segments, expected a handful", segs.len());
        let bin_segments =
            segs.iter().filter(|s| s.start().as_u64() <= 299 && s.end().as_u64() >= 200).count();
        assert_eq!(bin_segments, 1, "the irreducible bin must sit in exactly one segment");
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
            // No tile is empty; tile length is bounded by core (≤ max_len)
            // plus overlap on each side (so ≤ max_len + 2*overlap).
            for s in &segs {
                prop_assert!(s.len() >= 1);
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

        // r[verify unified.segment_byte_budget]
        /// Byte-aware subdivision must (a) keep cores tiling the range exactly
        /// and (b) keep every segment under budget unless its core is a single
        /// base (the irreducible leaf-bin case).
        #[test]
        fn byte_split_stays_under_budget_and_tiles(
            range_start in 0u32..50_000,
            len in 1u32..2_000,
            bytes_per_base in 1u64..1_000,
            budget in 1u64..200_000,
            overlap in 0u32..50,
        ) {
            let range_end = range_start.saturating_add(len - 1).min(99_999);
            let max_len = NonZeroU32::new(10_000).unwrap();
            let opts = SegmentOptions::new(max_len).with_overlap(overlap).unwrap();
            let segs =
                byte_aware(100_000, range_start..=range_end, opts, budget, span_estimate(bytes_per_base));

            // Cores tile the requested range exactly.
            prop_assert_eq!(*segs.first().unwrap().core_range().start(), p(range_start));
            prop_assert_eq!(*segs.last().unwrap().core_range().end(), p(range_end));
            for w in segs.windows(2) {
                prop_assert_eq!(
                    w[0].core_range().end().as_u64() + 1,
                    w[1].core_range().start().as_u64()
                );
            }
            // Each segment fits the budget, or it can't be made cheaper: the
            // greedy floor is one core position expanded by overlap on each
            // side (the smallest range it can emit), so an over-budget segment
            // costs at most that floor.
            let floor = (1 + 2 * u64::from(overlap)) * bytes_per_base;
            for s in &segs {
                let bytes = (s.end().as_u64() - s.start().as_u64() + 1) * bytes_per_base;
                prop_assert!(
                    bytes <= budget.max(floor),
                    "segment {:?} = {} B exceeds max(budget {}, floor {})",
                    (s.start().as_u64(), s.end().as_u64()), bytes, budget, floor
                );
            }
        }
    }
}
