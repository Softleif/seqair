//! Index statistics against `samtools idxstats` on the shared test BAM, whose
//! numbers are hard-coded below.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::indexing_slicing, reason = "test code")]

use seqair::bam::{CsiIndex, MappedCounts, ReferenceStats};
use seqair::reader::{IndexedReader, Readers};
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
#[test]
fn cram_has_no_reference_stats() {
    let readers = Readers::open(&data("test_v30.cram"), &data("test.fasta.gz")).unwrap();
    assert_eq!(readers.alignment().reference_stats(), None);
}
