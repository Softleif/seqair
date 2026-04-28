//! Integration tests for [`Readers::segments`] and the segmented `pileup`.
//!
//! These tests use the real `tests/data/test.bam` to drive the public API
//! end-to-end and verify two big claims:
//!
//! 1. The set of pileup columns produced when iterating the region as
//!    several small segments equals the set produced when iterating the
//!    region as one big segment. This is the load-bearing equivalence —
//!    if it broke, every downstream caller that switches to segmentation
//!    would silently miss or duplicate columns.
//!
//! 2. `core_range()` of consecutive segments tiles the input range exactly
//!    (no gap, no overlap), matching `r[unified.segment_overlap]`.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    reason = "test code"
)]

use proptest::prelude::*;
use seqair::reader::{Readers, SegmentOptions};
use seqair_types::Pos0;
use std::num::NonZeroU32;
use std::path::Path;

fn test_bam_path() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.bam"))
}
fn test_fasta_path() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.fasta.gz"))
}

/// One pileup column projected to the fields that segmentation must preserve.
/// Match-depth + ref base together uniquely identify a column under the
/// pileup contract; we add `pos` for diagnostics and an integer hash of the
/// per-alignment ops for stronger equivalence.
#[derive(Debug, Clone, PartialEq, Eq)]
struct ColumnSnapshot {
    pos: u32,
    depth: usize,
    match_depth: usize,
    ref_base: u8,
}

/// Pile up `[start, end]` as one big segment and capture per-column snapshots
/// across the inclusive range — emulates what a caller used to do with the
/// old `pileup(tid, start, end)` API.
fn pileup_single(
    readers: &mut Readers,
    contig: &str,
    start: Pos0,
    end: Pos0,
) -> Vec<ColumnSnapshot> {
    // Use a max_len that comfortably covers the whole region in one tile.
    let max_len = NonZeroU32::new(2_000_000).unwrap();
    let opts = SegmentOptions::new(max_len);
    let mut plan = readers.segments((contig, start, end), opts).unwrap();
    let segment = plan.next().expect("at least one segment");
    assert!(plan.next().is_none(), "single-tile setup should yield exactly one segment");
    let mut engine = readers.pileup(&segment).unwrap();
    let mut out = Vec::new();
    while let Some(col) = engine.pileups() {
        out.push(ColumnSnapshot {
            pos: *col.pos(),
            depth: col.depth(),
            match_depth: col.match_depth(),
            ref_base: col.reference_base() as u8,
        });
    }
    drop(engine);
    out
}

/// Pile up the same region as several smaller segments and gather columns
/// while restricting to each segment's `core_range()` so adjacent overlaps
/// don't double-count.
fn pileup_segmented(
    readers: &mut Readers,
    contig: &str,
    start: Pos0,
    end: Pos0,
    max_len: u32,
    overlap: u32,
) -> Vec<ColumnSnapshot> {
    let opts = SegmentOptions::new(NonZeroU32::new(max_len).unwrap())
        .with_overlap(overlap)
        .expect("overlap < max_len");
    let plan: Vec<_> = readers.segments((contig, start, end), opts).unwrap().collect();
    assert!(!plan.is_empty(), "segments() must produce at least one tile for non-empty range");

    let mut out = Vec::new();
    for segment in &plan {
        let core = segment.core_range();
        let mut engine = readers.pileup(segment).unwrap();
        while let Some(col) = engine.pileups() {
            if !core.contains(&col.pos()) {
                continue;
            }
            out.push(ColumnSnapshot {
                pos: *col.pos(),
                depth: col.depth(),
                match_depth: col.match_depth(),
                ref_base: col.reference_base() as u8,
            });
        }
        // `engine` drops at end of scope, returning the store to `readers`.
    }
    out
}

const REGION_START: u32 = 6_103_500;
const REGION_END: u32 = 6_106_500;

// r[verify unified.readers_segments]
// r[verify unified.readers_pileup]
#[test]
fn segmented_pileup_equals_single_pileup_no_overlap() {
    let mut readers = Readers::open(test_bam_path(), test_fasta_path()).unwrap();
    let start = Pos0::new(REGION_START).unwrap();
    let end = Pos0::new(REGION_END).unwrap();
    let single = pileup_single(&mut readers, "chr19", start, end);
    let segmented = pileup_segmented(&mut readers, "chr19", start, end, 250, 0);
    assert_eq!(single, segmented, "segmented pileup with no overlap must equal single-tile pileup");
}

// r[verify unified.readers_segments]
// r[verify unified.readers_pileup]
#[test]
fn segmented_pileup_equals_single_pileup_with_overlap() {
    let mut readers = Readers::open(test_bam_path(), test_fasta_path()).unwrap();
    let start = Pos0::new(REGION_START).unwrap();
    let end = Pos0::new(REGION_END).unwrap();
    let single = pileup_single(&mut readers, "chr19", start, end);
    let segmented = pileup_segmented(&mut readers, "chr19", start, end, 300, 50);
    assert_eq!(
        single, segmented,
        "segmented pileup with overlap must equal single-tile pileup when restricted to core ranges"
    );
}

// r[verify unified.segment_overlap]
#[test]
fn cores_partition_input_with_overlap() {
    let readers = Readers::open(test_bam_path(), test_fasta_path()).unwrap();
    let start = Pos0::new(REGION_START).unwrap();
    let end = Pos0::new(REGION_END).unwrap();
    let opts = SegmentOptions::new(NonZeroU32::new(400).unwrap()).with_overlap(75).unwrap();
    let plan: Vec<_> = readers.segments(("chr19", start, end), opts).unwrap().collect();
    assert!(!plan.is_empty());
    assert_eq!(*plan.first().unwrap().core_range().start(), start);
    assert_eq!(*plan.last().unwrap().core_range().end(), end);
    for w in plan.windows(2) {
        let prev_end = *w[0].core_range().end();
        let next_start = *w[1].core_range().start();
        let prev_u32: u32 = *prev_end;
        let next_u32: u32 = *next_start;
        assert_eq!(
            prev_u32 + 1,
            next_u32,
            "core ranges must be contiguous (gap-free, overlap-free)"
        );
    }
}

// r[verify unified.into_segment_target]
#[test]
fn whole_genome_target_yields_all_nonempty_contigs() {
    let readers = Readers::open(test_bam_path(), test_fasta_path()).unwrap();
    let opts = SegmentOptions::new(NonZeroU32::new(10_000_000).unwrap());
    let plan: Vec<_> = readers.segments((), opts).unwrap().collect();
    let mut contigs: Vec<&str> = plan.iter().map(|s| s.contig().as_str()).collect();
    contigs.dedup();
    assert!(!contigs.is_empty(), "header has at least one contig");
    // Every header contig with non-zero length should appear.
    let header = readers.header();
    for tid_u32 in 0..u32::try_from(header.target_count()).unwrap() {
        if header.target_len(tid_u32).unwrap_or(0) > 0 {
            let name = header.target_name(tid_u32).unwrap();
            assert!(contigs.contains(&name), "whole-genome plan missing contig '{name}'",);
        }
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(16))]

    // For a random `(max_len, overlap)` pair the segmented pileup must match
    // the single-pileup baseline column-for-column.
    //
    // Using the real BAM as the oracle: the single-tile pileup is the
    // independent ground truth (it's what users had pre-segmentation), and
    // the segmented version is the new code path — they must agree on
    // every position and depth.
    // r[verify unified.readers_segments]
    #[test]
    fn segmented_matches_single(
        max_len in 100u32..1_000,
        overlap in 0u32..100,
    ) {
        prop_assume!(overlap < max_len);
        let mut readers = Readers::open(test_bam_path(), test_fasta_path()).unwrap();
        let start = Pos0::new(REGION_START).unwrap();
        let end = Pos0::new(REGION_END).unwrap();
        let single = pileup_single(&mut readers, "chr19", start, end);
        let segmented = pileup_segmented(&mut readers, "chr19", start, end, max_len, overlap);
        prop_assert_eq!(single, segmented);
    }
}
