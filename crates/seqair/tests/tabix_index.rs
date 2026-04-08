//! Tests for tabix (.tbi) index parsing.
//! Creates a bgzf-compressed SAM from the test BAM, indexes it with samtools/tabix,
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    reason = "test code"
)]

use seqair::bam::index::BamIndex;
use seqair_types::{Pos, Zero};
use std::path::Path;
use std::process::Command;

fn test_bam_path() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.bam"))
}

const CONTIGS: &[(&str, u64, u64)] = &[
    ("chr19", 6_103_076, 6_143_229),
    ("2kb_3_Unmodified", 1, 2_018),
    ("bacteriophage_lambda_CpG", 1, 48_502),
];

fn create_sam_gz_with_tabix(tmpdir: &tempfile::TempDir) -> std::path::PathBuf {
    let sam_gz = tmpdir.path().join("test.sam.gz");

    let status = Command::new("samtools")
        .args(["view", "-h", "--output-fmt", "SAM,level=6", "-o"])
        .arg(&sam_gz)
        .arg(test_bam_path())
        .status()
        .expect("samtools not found");
    assert!(status.success(), "samtools view failed");

    let status =
        Command::new("tabix").args(["-p", "sam"]).arg(&sam_gz).status().expect("tabix not found");
    assert!(status.success(), "tabix indexing failed");

    sam_gz
}

// r[verify tabix.magic]
// r[verify tabix.header]
// r[verify tabix.bai_reuse]
// r[verify tabix.index_data]
#[test]
fn tabix_parses_without_error() {
    let tmpdir = tempfile::tempdir().unwrap();
    let sam_gz = create_sam_gz_with_tabix(&tmpdir);
    let tbi_path = sam_gz.with_extension("gz.tbi");
    assert!(tbi_path.exists(), "tabix index should exist at {tbi_path:?}");

    let index = BamIndex::from_tabix_path(&tbi_path).expect("tabix parse");

    // Verify we can query — the index should have references
    let header = seqair::bam::BamHeader::from_bam_path(test_bam_path()).unwrap();
    let tid = header.tid("chr19").unwrap();
    let chunks = index.query(
        tid,
        Pos::<Zero>::new(6_105_700).unwrap(),
        Pos::<Zero>::new(6_105_800).unwrap(),
    );
    assert!(!chunks.is_empty(), "tabix query should return chunks for chr19");
}

// r[verify tabix.query]
// r[verify tabix.pseudo_bin]
#[test]
fn tabix_query_matches_bai_query() {
    let tmpdir = tempfile::tempdir().unwrap();
    let sam_gz = create_sam_gz_with_tabix(&tmpdir);
    let tbi_path = sam_gz.with_extension("gz.tbi");

    let tabix_index = BamIndex::from_tabix_path(&tbi_path).expect("tabix parse");
    let bai_path = test_bam_path().with_extension("bam.bai");
    let bai_index = BamIndex::from_path(&bai_path).expect("BAI parse");

    let header = seqair::bam::BamHeader::from_bam_path(test_bam_path()).unwrap();

    for &(contig, start, end) in CONTIGS {
        let tid = header.tid(contig).unwrap();
        let start_pos = Pos::<Zero>::try_from_u64(start).expect("start fits in u32");
        let end_pos = Pos::<Zero>::try_from_u64(end).expect("end fits in u32");
        let bai_chunks = bai_index.query(tid, start_pos, end_pos);
        let tabix_chunks = tabix_index.query(tid, start_pos, end_pos);

        // Both should return non-empty results (the test data has reads in all contigs)
        assert!(!bai_chunks.is_empty(), "{contig}: BAI query returned no chunks");
        assert!(!tabix_chunks.is_empty(), "{contig}: tabix query returned no chunks");

        // Chunk counts should match (same binning scheme, same data)
        assert_eq!(
            bai_chunks.len(),
            tabix_chunks.len(),
            "{contig}: BAI has {} chunks but tabix has {}",
            bai_chunks.len(),
            tabix_chunks.len(),
        );
    }
}

// r[verify tabix.compression]
#[test]
fn tabix_rejects_non_bgzf() {
    let tmpdir = tempfile::tempdir().unwrap();
    let fake_tbi = tmpdir.path().join("fake.tbi");
    std::fs::write(&fake_tbi, b"not a bgzf file at all").unwrap();

    let result = BamIndex::from_tabix_path(&fake_tbi);
    assert!(result.is_err());
}

// TODO (nice-to-have): Compare tabix query results against samtools CLI output.
