//! Index statistics and per-fetch counts against `samtools idxstats` on the
//! shared test BAM, whose numbers are hard-coded below.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::indexing_slicing, reason = "test code")]

use seqair::bam::{CsiIndex, MappedCounts, ReferenceStats};
use seqair::reader::{DepthLimit, FetchCounts, IndexedReader, Readers, SegmentOptions};
use std::num::NonZeroU32;
use std::path::Path;

fn data(name: &str) -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/data").join(name)
}

/// `samtools idxstats tests/data/test.bam`: (contig, length, mapped, unmapped).
const IDXSTATS: [(&str, u32, u64, u64); 3] = [
    ("chr19", 61_431_566, 12_297, 70),
    ("2kb_3_Unmodified", 2_018, 2_146, 67),
    ("bacteriophage_lambda_CpG", 48_502, 10_500, 278),
];

/// Leaf bins with a chunk, per contig, as both `test.bam.bai` and the CSI
/// `samtools index -c` built for the same file record them.
const LEAVES: [u64; 3] = [3, 1, 3];

fn expected() -> Vec<ReferenceStats> {
    IDXSTATS
        .iter()
        .zip(LEAVES)
        .map(|(&(_, _, mapped, unmapped), covered_leaf_bins)| ReferenceStats {
            counts: Some(MappedCounts { mapped, unmapped }),
            covered_leaf_bins,
            leaf_bin_len: 16_384,
        })
        .collect()
}

// r[verify index.reference_stats]
#[test]
fn bai_reference_stats_match_samtools_idxstats() {
    let reader = IndexedReader::open(&data("test.bam")).unwrap();
    assert_eq!(reader.reference_stats().unwrap(), expected());
}

// r[verify index.reference_stats]
/// `test.bam.samtools.csi` is `samtools index -c tests/data/test.bam`, kept
/// under a name the BAM reader does not pick up in place of the `.bai`.
#[test]
fn csi_reference_stats_match_samtools_idxstats() {
    let index = CsiIndex::from_path(&data("test.bam.samtools.csi")).unwrap();
    assert_eq!(index.reference_stats(), expected());
}

// r[verify index.reference_stats]
// r[verify index.unplaced_unmapped]
#[test]
fn cram_has_no_reference_stats() {
    let readers = Readers::open(&data("test_v30.cram"), &data("test.fasta.gz")).unwrap();
    assert_eq!(readers.alignment().reference_stats(), None);
    assert_eq!(readers.alignment().unplaced_unmapped(), None);
}

// r[verify index.unplaced_unmapped]
/// idxstats' `*` line for test.bam: 4,000 reads with no reference.
#[test]
fn unplaced_unmapped_matches_samtools_idxstats() {
    let reader = IndexedReader::open(&data("test.bam")).unwrap();
    assert_eq!(reader.unplaced_unmapped(), Some(4_000));
    let csi = CsiIndex::from_path(&data("test.bam.samtools.csi")).unwrap();
    assert_eq!(csi.unplaced_unmapped(), Some(4_000));
}

/// The fetch counts of one whole-contig pileup per contig.
fn whole_contig_fetches(alignment: &str) -> Vec<FetchCounts> {
    let mut readers = Readers::open(&data(alignment), &data("test.fasta.gz")).unwrap();
    IDXSTATS
        .iter()
        .map(|&(contig, len, _, _)| {
            let opts = SegmentOptions::new(NonZeroU32::new(len).unwrap()).without_byte_budget();
            let segments: Vec<_> = readers.segments(contig, opts).unwrap().collect();
            assert_eq!(segments.len(), 1, "{contig} is one tile");
            readers.pileup(&segments[0], DepthLimit::Unlimited).run().unwrap().fetch_counts()
        })
        .collect()
}

// r[verify unified.fetch_counts+3]
// r[verify unified.pileup_fetch_counts]
/// A pileup over a whole contig fetches every placed record of it, so its
/// counts are the index's: mapped as idxstats has it, and the placed
/// unmapped ones on top.
#[test]
fn a_whole_contig_pileup_fetches_what_idxstats_counts() {
    let expected: Vec<FetchCounts> = IDXSTATS
        .iter()
        .map(|&(_, _, mapped, unmapped)| {
            let fetched = usize::try_from(mapped + unmapped).unwrap();
            FetchCounts {
                fetched,
                mapped: usize::try_from(mapped).unwrap(),
                started: fetched,
                kept: fetched,
            }
        })
        .collect();
    assert_eq!(whole_contig_fetches("test.bam"), expected, "BAM");
    assert_eq!(whole_contig_fetches("test_v30.cram"), expected, "CRAM");
}

// r[verify unified.pileup_fetch_counts]
/// A customize value that drops records leaves `fetched` and `mapped` alone.
#[test]
fn a_filter_changes_kept_only() {
    let mut readers = Readers::open_customized(
        &data("test.bam"),
        &data("test.fasta.gz"),
        seqair::bam::record_store::RejectUnmapped,
    )
    .unwrap();
    let (contig, len, mapped, unmapped) = IDXSTATS[1];
    let opts = SegmentOptions::new(NonZeroU32::new(len).unwrap()).without_byte_budget();
    let segment = readers.segments(contig, opts).unwrap().next().unwrap();
    let counts = readers.pileup(&segment, DepthLimit::Unlimited).run().unwrap().fetch_counts();
    let mapped = usize::try_from(mapped).unwrap();
    let fetched = mapped + usize::try_from(unmapped).unwrap();
    assert_eq!(counts, FetchCounts { fetched, mapped, started: fetched, kept: mapped });
}

// r[verify unified.fetch_counts_started]
/// Tiles of 300 bp over the 2,018 bp contig: reads crossing a boundary are
/// fetched by both tiles, but start in only one, so the `started` counts add
/// up to idxstats' total for the contig and the `fetched` counts exceed it.
#[test]
fn started_counts_tile_without_double_counting() {
    for alignment in ["test.bam", "test_v30.cram"] {
        let mut readers = Readers::open(&data(alignment), &data("test.fasta.gz")).unwrap();
        let (contig, _, mapped, unmapped) = IDXSTATS[1];
        let opts = SegmentOptions::new(NonZeroU32::new(300).unwrap()).without_byte_budget();
        let segments: Vec<_> = readers.segments(contig, opts).unwrap().collect();
        assert!(segments.len() > 1);
        let (mut started, mut fetched) = (0, 0);
        for segment in &segments {
            assert_eq!(segment.overlap_start() + segment.overlap_end(), 0, "tiles are disjoint");
            let guard = readers.pileup(segment, DepthLimit::Unlimited).run().unwrap();
            started += guard.fetch_counts().started;
            fetched += guard.fetch_counts().fetched;
        }
        assert_eq!(started, usize::try_from(mapped + unmapped).unwrap(), "{alignment}");
        assert!(fetched > started, "{alignment}");
    }
}
