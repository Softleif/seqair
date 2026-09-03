//! Load-time depth cap for [`fetch_into_customized`](crate::reader::Readers).
//!
//! [`DepthCap`] wraps any [`CustomizeRecordStore`] and rejects reads in
//! [`filter_raw`](CustomizeRecordStore::filter_raw) — *before* any slab write
//! or base decode — once more than `cap` already-kept reads overlap the
//! incoming read's start position. This bounds the [`RecordStore`] to `cap`
//! reads per reference column regardless of how deep the underlying pileup is.
//!
//! Without it, a region fetch loads *every* overlapping read into the store
//! before any depth limit applies, so a collapsed-repeat decoy contig or an
//! EBV genome at tens of thousands × coverage materialises millions of reads
//! per segment. Because excess reads are dropped in `filter_raw`, the dropped
//! reads cost essentially nothing (no memcpy into slabs, no base decode, no
//! `compute`).
//!
//! The window is reset whenever records jump to a new `tid` or backwards in
//! coordinate (the start of a new fetch), so a [`DepthCap`] reused across
//! segments stays correct without an explicit per-fetch reset hook.
//!
//! Capping is greedy and streaming: at a saturated position the *earliest*
//! reads are kept and later-starting reads are dropped. This can truncate the
//! tail of a region whose depth stays above `cap` — coverage just past such a
//! run may dip below what a uniform random downsample would give. That only
//! happens where depth already exceeds `cap`, so any column the caller would
//! keep in full (depth `<= cap`) away from a saturated run is unaffected.

use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::num::NonZeroU32;

use seqair_types::Pos0;

use super::record_store::{CustomizeRecordStore, FilterRawFields, RecordStore, SlimRecord};

/// Per-column read-depth policy for store-filling fetches such as
/// [`Readers::pileup`](crate::reader::Readers::pileup).
///
/// There is no default on purpose: a region fetch loads every overlapping read
/// into memory, so callers must consciously choose between bounding that
/// (`PerColumn`) and accepting an unbounded load (`Unlimited`). `Unlimited` is
/// spelled out at the call site so it greps as a deliberate decision rather
/// than a silently-inherited default.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DepthLimit {
    /// Load every read in the region. Memory scales with coverage — a
    /// collapsed-repeat or viral contig at tens of thousands × can OOM.
    Unlimited,
    /// Keep at most this many reads overlapping any reference position,
    /// dropping the excess at load time (see [`DepthCap`]).
    PerColumn(NonZeroU32),
}

impl DepthLimit {
    /// The cap as an `Option`, `None` for [`Unlimited`](Self::Unlimited).
    pub fn per_column(self) -> Option<NonZeroU32> {
        match self {
            Self::Unlimited => None,
            Self::PerColumn(cap) => Some(cap),
        }
    }
}

/// Wraps a [`CustomizeRecordStore`] to cap fetched read depth at load time.
#[derive(Clone)]
pub struct DepthCap<E> {
    inner: E,
    cap: Option<NonZeroU32>,
    /// Last covered position of each currently-overlapping kept read
    /// (inclusive, min-heap via [`Reverse`]).
    active: BinaryHeap<Reverse<u64>>,
    /// Sentinel `i32::MIN` so the first record always resets the window.
    last_tid: i32,
    last_start: u64,
}

impl<E> DepthCap<E> {
    /// Wrap `inner` with capping disabled. Transparent until
    /// [`set_cap`](Self::set_cap) installs a budget.
    pub fn new(inner: E) -> Self {
        Self::with_cap(inner, None)
    }

    /// Wrap `inner` with an initial cap (`None` disables capping).
    pub fn with_cap(inner: E, cap: Option<NonZeroU32>) -> Self {
        Self { inner, cap, active: BinaryHeap::new(), last_tid: i32::MIN, last_start: 0 }
    }

    /// Set the per-column read cap. `None` disables capping.
    pub fn set_cap(&mut self, cap: Option<NonZeroU32>) {
        self.cap = cap;
    }

    pub fn cap(&self) -> Option<NonZeroU32> {
        self.cap
    }

    pub fn inner(&self) -> &E {
        &self.inner
    }

    pub fn inner_mut(&mut self) -> &mut E {
        &mut self.inner
    }

    pub fn into_inner(self) -> E {
        self.inner
    }
}

impl<E: CustomizeRecordStore> CustomizeRecordStore for DepthCap<E> {
    type Extra = E::Extra;

    fn filter_raw(&mut self, fields: &FilterRawFields<'_>) -> bool {
        // Delegate first: a read the inner filter rejects (flags, mapq, …) must
        // not consume a depth slot.
        if !self.inner.filter_raw(fields) {
            return false;
        }
        let Some(cap) = self.cap else {
            return true;
        };

        let start = fields.pos.as_u64();
        // A new fetch — new contig or any backwards coordinate jump — restarts
        // the window. Records within one fetch arrive in ascending order.
        if fields.tid != self.last_tid || start < self.last_start {
            self.active.clear();
        }
        self.last_tid = fields.tid;
        self.last_start = start;

        // r[impl record_store.depth_cap.window]
        // Last covered position, inclusive. `end_pos` is pre-computed for
        // SAM/CRAM; for BAM (`push_raw`) it is `None`, so derive the reference
        // span from the raw CIGAR. Guarantee a minimum 1bp span so a read
        // always occupies a slot.
        let end = fields.end_pos.map(Pos0::as_u64).unwrap_or_else(|| {
            let span = ref_span_from_cigar(fields.cigar).max(1);
            start.saturating_add(span.saturating_sub(1))
        });

        // Evict reads whose last covered position is behind `start`.
        while let Some(&Reverse(min_end)) = self.active.peek() {
            if min_end < start {
                self.active.pop();
            } else {
                break;
            }
        }

        // Reject once the column is already full. Comparing in `u32` keeps the
        // cast-free path; an `active` larger than `u32::MAX` is over any cap.
        if u32::try_from(self.active.len()).map_or(true, |n| n >= cap.get()) {
            return false;
        }
        self.active.push(Reverse(end));
        true
    }

    fn compute(&mut self, rec: &SlimRecord, store: &RecordStore<Self::Extra>) -> Self::Extra {
        self.inner.compute(rec, store)
    }

    fn filter(&mut self, rec: &SlimRecord, store: &RecordStore<Self::Extra>) -> bool {
        self.inner.filter(rec, store)
    }
}

/// Sum of reference-consuming CIGAR op lengths from raw BAM CIGAR bytes
/// (`n_cigar_op` little-endian `u32`s: `len << 4 | op`).
#[allow(clippy::indexing_slicing, reason = "chunks_exact(4) yields exactly 4 bytes")]
fn ref_span_from_cigar(cigar: &[u8]) -> u64 {
    let mut span = 0u64;
    for chunk in cigar.chunks_exact(4) {
        let val = u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]);
        // Reference-consuming ops: M(0), D(2), N(3), =(7), X(8).
        if matches!(val & 0xf, 0 | 2 | 3 | 7 | 8) {
            span = span.saturating_add(u64::from(val >> 4));
        }
    }
    span
}

#[cfg(test)]
#[allow(clippy::arithmetic_side_effects, reason = "test arithmetic on known small values")]
mod tests {
    use super::*;
    use crate::bam::record_store::Sequence;
    use seqair_types::BamFlags;

    /// One `<len>M` CIGAR op as raw BAM bytes.
    fn match_cigar(len: u32) -> Vec<u8> {
        (len << 4).to_le_bytes().to_vec()
    }

    /// The fields a `push_raw` (BAM) record presents: no `end_pos`, raw CIGAR.
    fn bam_fields<'a>(pos: Pos0, cigar: &'a [u8]) -> FilterRawFields<'a> {
        FilterRawFields {
            pos,
            end_pos: None,
            flags: BamFlags::empty(),
            mapq: 60,
            n_cigar_ops: 1,
            seq_len: 0,
            matching_bases: 0,
            indel_bases: 0,
            tid: 0,
            next_ref_id: -1,
            next_pos: -1,
            template_len: 0,
            qname: b"r",
            qual_bytes: &[],
            aux_bytes: &[],
            cigar,
            seq: Sequence::Bases(&[]),
        }
    }

    /// The same record as a `push_fields` (SAM/CRAM) record: `end_pos` supplied.
    fn field_fields<'a>(pos: Pos0, end_pos: Pos0, cigar: &'a [u8]) -> FilterRawFields<'a> {
        FilterRawFields { end_pos: Some(end_pos), ..bam_fields(pos, cigar) }
    }

    fn pos(p: u32) -> Pos0 {
        Pos0::new(p).expect("valid position")
    }

    // r[verify record_store.depth_cap.window]
    /// A kept read occupies a slot at every position it covers, up to and
    /// including its last one — and the boundary lands identically whether the
    /// span was supplied (SAM/CRAM) or derived from the CIGAR (BAM).
    #[test]
    fn the_cap_window_ends_at_the_last_covered_position() {
        let cigar = match_cigar(5); // covers pos ..= pos + 4
        let cap = NonZeroU32::new(1);

        // Second read starts on the first read's last base: still overlapping.
        for supplied in [false, true] {
            let mut capped = DepthCap::with_cap((), cap);
            let first = |p: u32| {
                if supplied {
                    field_fields(pos(p), pos(p + 4), &cigar)
                } else {
                    bam_fields(pos(p), &cigar)
                }
            };
            assert!(capped.filter_raw(&first(100)), "first read is always kept");
            assert!(
                !capped.filter_raw(&first(104)),
                "a read starting on the last covered base overlaps (supplied end_pos: {supplied})"
            );
        }

        // One position further along, the first read is out of the window.
        for supplied in [false, true] {
            let mut capped = DepthCap::with_cap((), cap);
            let build = |p: u32| {
                if supplied {
                    field_fields(pos(p), pos(p + 4), &cigar)
                } else {
                    bam_fields(pos(p), &cigar)
                }
            };
            assert!(capped.filter_raw(&build(100)));
            assert!(
                capped.filter_raw(&build(105)),
                "a read starting past the last covered base is kept (supplied end_pos: {supplied})"
            );
        }
    }

    // r[verify record_store.depth_cap.window]
    /// A CIGAR that consumes no reference still occupies exactly one position,
    /// so two such reads at the same start cannot both pass a cap of one.
    #[test]
    fn a_zero_span_read_occupies_one_position() {
        let cigar = match_cigar(0);
        let mut capped = DepthCap::with_cap((), NonZeroU32::new(1));
        assert!(capped.filter_raw(&bam_fields(pos(100), &cigar)));
        assert!(!capped.filter_raw(&bam_fields(pos(100), &cigar)));
        assert!(capped.filter_raw(&bam_fields(pos(101), &cigar)));
    }
}
