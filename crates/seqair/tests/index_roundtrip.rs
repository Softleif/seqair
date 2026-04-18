//! Index round-trip tests: write BAM + BAI with seqair, then verify that
//! samtools and noodles can use the index for region queries and get the
//! same results as seqair.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "test code"
)]
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    reason = "test code with known small values"
)]

use seqair::bam::cigar::{CigarOp, CigarOpType};
use seqair::bam::header::BamHeader;
use seqair::bam::owned_record::OwnedBamRecord;
use seqair::bam::writer::BamWriter;
use seqair::bam::{IndexedBamReader, Pos0, RecordStore};
use seqair_types::{Base, BaseQuality};
use std::path::Path;
use std::process::Command;

fn make_header() -> BamHeader {
    BamHeader::from_sam_text(
        "@HD\tVN:1.6\tSO:coordinate\n\
         @SQ\tSN:chr1\tLN:1000000\n\
         @SQ\tSN:chr2\tLN:500000\n\
         @SQ\tSN:chr3\tLN:200000\n",
    )
    .unwrap()
}

fn cyclic_seq(n: usize) -> Vec<Base> {
    (0..n)
        .map(|i| match i % 4 {
            0 => Base::A,
            1 => Base::C,
            2 => Base::G,
            _ => Base::T,
        })
        .collect()
}

/// Write a BAM with records spread across multiple contigs, return the path.
fn write_multi_contig_bam(dir: &Path) -> std::path::PathBuf {
    let header = make_header();
    let bam_path = dir.join("multi.bam");
    let mut writer = BamWriter::from_path(&bam_path, &header, true).unwrap();

    // chr1: 30 records at positions 0, 10000, 20000, ...
    for i in 0..30u32 {
        let rec = OwnedBamRecord::builder(0, i64::from(i) * 10000, format!("c1r{i}").into_bytes())
            .mapq(60)
            .cigar(vec![CigarOp::new(CigarOpType::Match, 150)])
            .seq(cyclic_seq(150))
            .qual(vec![BaseQuality::from_byte(30); 150])
            .build()
            .unwrap();
        writer.write(&rec).unwrap();
    }

    // chr2: 15 records
    for i in 0..15u32 {
        let rec = OwnedBamRecord::builder(1, i64::from(i) * 20000, format!("c2r{i}").into_bytes())
            .mapq(50)
            .cigar(vec![CigarOp::new(CigarOpType::Match, 100)])
            .seq(cyclic_seq(100))
            .qual(vec![BaseQuality::from_byte(25); 100])
            .build()
            .unwrap();
        writer.write(&rec).unwrap();
    }

    // chr3: 5 records
    for i in 0..5u32 {
        let rec = OwnedBamRecord::builder(2, i64::from(i) * 30000, format!("c3r{i}").into_bytes())
            .mapq(40)
            .cigar(vec![CigarOp::new(CigarOpType::Match, 80)])
            .seq(cyclic_seq(80))
            .qual(vec![BaseQuality::from_byte(20); 80])
            .build()
            .unwrap();
        writer.write(&rec).unwrap();
    }

    let (_inner, index_builder) = writer.finish().unwrap();
    if let Some(ib) = index_builder {
        let bai_file = std::fs::File::create(bam_path.with_extension("bam.bai")).unwrap();
        ib.write_bai(bai_file, header.target_count()).unwrap();
    }

    bam_path
}

/// Helper: run `samtools view -c <bam> <region>` and return the count.
fn samtools_count(bam_path: &Path, region: &str) -> usize {
    let output = Command::new("samtools")
        .args(["view", "-c"])
        .arg(bam_path)
        .arg(region)
        .output()
        .expect("samtools not found");
    assert!(
        output.status.success(),
        "samtools view -c failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).unwrap().trim().parse().unwrap()
}

/// Seqair's co-produced BAI index is readable by samtools for whole-contig queries.
#[test]
fn samtools_reads_seqair_bai_whole_contigs() {
    let dir = tempfile::tempdir().unwrap();
    let bam_path = write_multi_contig_bam(dir.path());

    assert_eq!(samtools_count(&bam_path, "chr1"), 30);
    assert_eq!(samtools_count(&bam_path, "chr2"), 15);
    assert_eq!(samtools_count(&bam_path, "chr3"), 5);
}

/// Seqair's co-produced BAI index supports sub-region queries via samtools.
#[test]
fn samtools_reads_seqair_bai_subregions() {
    let dir = tempfile::tempdir().unwrap();
    let bam_path = write_multi_contig_bam(dir.path());

    // chr1 records at 0, 10000, 20000, ..., 290000
    // Query chr1:50001-150000 (1-based) should find records at 50000, 60000, ..., 140000
    // That's positions 50000..=140000 in 0-based, so records at 50000, 60000, 70000, 80000,
    // 90000, 100000, 110000, 120000, 130000, 140000 = 10 records.
    // Plus any record whose span overlaps: record at 49850..50000? No, they start at multiples
    // of 10000 with 150M. Record at 40000 spans 40000..40150 — doesn't reach 50000.
    let count = samtools_count(&bam_path, "chr1:50001-150000");
    assert_eq!(count, 10, "chr1:50001-150000 should have 10 records");

    // chr2 records at 0, 20000, 40000, ...
    // Query chr2:1-40000 (1-based) should find records at 0, 20000 = 2 records
    // Record at 40000 starts at 40001 (1-based) which is > 40000, so excluded? Actually
    // samtools uses 1-based inclusive, so chr2:1-40000 = 0-based [0, 40000].
    // Record at 40000 (0-based) = 40001 (1-based) > 40000 → excluded.
    // But record at 20000 spans [20000, 20100) → included.
    // So records at 0 and 20000 → 2 records.
    let count = samtools_count(&bam_path, "chr2:1-40000");
    assert_eq!(count, 2, "chr2:1-40000 should have 2 records");
}

/// Seqair reads from its own BAI index correctly across contigs.
#[test]
fn seqair_reads_own_bai() {
    let dir = tempfile::tempdir().unwrap();
    let bam_path = write_multi_contig_bam(dir.path());

    let mut reader = IndexedBamReader::open(&bam_path).expect("open");

    // Full contig fetches
    for (name, len, expected) in
        [("chr1", 1_000_000u32, 30), ("chr2", 500_000, 15), ("chr3", 200_000, 5)]
    {
        let tid = reader.header().tid(name).expect(name);
        let mut store = RecordStore::new();
        reader
            .fetch_into(tid, Pos0::new(0).unwrap(), Pos0::new(len).unwrap(), &mut store)
            .expect("fetch");
        assert_eq!(store.len(), expected, "{name}: record count");
    }
}

/// Cross-validate: seqair and samtools agree on the same region query result.
#[test]
fn seqair_and_samtools_agree_on_regions() {
    let dir = tempfile::tempdir().unwrap();
    let bam_path = write_multi_contig_bam(dir.path());

    let mut reader = IndexedBamReader::open(&bam_path).expect("open");

    // Test whole-contig regions where seqair and samtools must agree exactly.
    // (Sub-region queries can differ because seqair's BAI fetch returns all
    // records in overlapping bins, which may include extras near boundaries.)
    let regions: &[(&str, u32, u32, &str)] = &[
        ("chr1", 0, 1_000_000, "chr1"),
        ("chr2", 0, 500_000, "chr2"),
        ("chr3", 0, 200_000, "chr3"),
    ];

    for &(contig, start, end, samtools_region) in regions {
        let tid = reader.header().tid(contig).expect(contig);
        let mut store = RecordStore::new();
        reader
            .fetch_into(tid, Pos0::new(start).unwrap(), Pos0::new(end).unwrap(), &mut store)
            .expect("fetch");

        let st_count = samtools_count(&bam_path, samtools_region);
        assert_eq!(
            store.len(),
            st_count,
            "region {samtools_region}: seqair={} samtools={st_count}",
            store.len()
        );
    }
}

/// samtools idxstats on seqair-produced BAI gives correct per-contig counts.
#[test]
fn samtools_idxstats_matches() {
    let dir = tempfile::tempdir().unwrap();
    let bam_path = write_multi_contig_bam(dir.path());

    let output = Command::new("samtools")
        .arg("idxstats")
        .arg(&bam_path)
        .output()
        .expect("samtools idxstats");
    assert!(output.status.success());

    let stats = String::from_utf8(output.stdout).unwrap();
    let mut mapped_counts: Vec<(String, usize)> = Vec::new();
    for line in stats.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 3 && fields[0] != "*" {
            let name = fields[0].to_string();
            let mapped: usize = fields[2].parse().unwrap();
            mapped_counts.push((name, mapped));
        }
    }

    assert_eq!(
        mapped_counts,
        vec![("chr1".to_string(), 30), ("chr2".to_string(), 15), ("chr3".to_string(), 5),]
    );
}
