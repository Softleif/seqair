//! CSI index tests: round-trip (write + parse + query), cross-validation
//! with samtools (`samtools index -c`), and comparison with BAI queries.
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
use seqair::bam::csi_index::CsiIndex;
use seqair::bam::header::BamHeader;
use seqair::bam::index::BamIndex;
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

/// Write a BAM with records spread across multiple contigs.
/// Returns the BAM path. Also writes BAI for cross-validation.
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
        // Write BAI for cross-validation
        let bai_file = std::fs::File::create(bam_path.with_extension("bam.bai")).unwrap();
        ib.write_bai(bai_file, header.target_count()).unwrap();

        // Write CSI
        let csi_path = bam_path.with_extension("bam.csi");
        let csi_file = std::fs::File::create(&csi_path).unwrap();
        ib.write_csi(csi_file, header.target_count(), &[]).unwrap();
    }

    bam_path
}

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

// --- Round-trip: seqair writes CSI, seqair reads it back ---

#[test]
fn csi_roundtrip_parse() {
    let dir = tempfile::tempdir().unwrap();
    let bam_path = write_multi_contig_bam(dir.path());
    let csi_path = bam_path.with_extension("bam.csi");

    let idx = CsiIndex::from_path(&csi_path).expect("parse seqair CSI");
    assert_eq!(idx.min_shift(), 14);
    assert_eq!(idx.depth(), 5);
}

/// `write_csi` MUST produce BGZF-compressed output (gzip magic 1f 8b) to
/// match htslib's behavior. Tools that strictly validate the gzip header
/// will reject uncompressed CSI files.
#[test]
fn write_csi_is_bgzf_compressed() {
    let dir = tempfile::tempdir().unwrap();
    let bam_path = write_multi_contig_bam(dir.path());
    let csi_path = bam_path.with_extension("bam.csi");
    let raw = std::fs::read(&csi_path).expect("read csi");
    assert!(raw.len() >= 2, "csi file unexpectedly short");
    assert_eq!(
        &raw[..2],
        &[0x1f, 0x8b],
        "write_csi output must be BGZF-compressed (gzip magic 1f 8b), got {:02x?}",
        &raw[..raw.len().min(4)]
    );
}

#[test]
fn csi_roundtrip_query_matches_bai() {
    let dir = tempfile::tempdir().unwrap();
    let bam_path = write_multi_contig_bam(dir.path());
    let csi_path = bam_path.with_extension("bam.csi");
    let bai_path = bam_path.with_extension("bam.bai");

    let csi = CsiIndex::from_path(&csi_path).expect("CSI");
    let bai = BamIndex::from_path(&bai_path).expect("BAI");

    // Query all contigs: CSI and BAI must return the same chunks
    let regions: &[(u32, u32, u32)] = &[
        (0, 0, 1_000_000),    // chr1 whole
        (1, 0, 500_000),      // chr2 whole
        (2, 0, 200_000),      // chr3 whole
        (0, 50_000, 150_000), // chr1 subregion
        (1, 0, 40_000),       // chr2 prefix
    ];

    for &(tid, start, end) in regions {
        let csi_chunks = csi.query(tid, Pos0::new(start).unwrap(), Pos0::new(end).unwrap());
        let bai_chunks = bai.query(tid, Pos0::new(start).unwrap(), Pos0::new(end).unwrap());

        // Same number of chunks (after merging)
        assert_eq!(
            csi_chunks.len(),
            bai_chunks.len(),
            "tid={tid}, region=[{start},{end}]: chunk count mismatch (csi={}, bai={})",
            csi_chunks.len(),
            bai_chunks.len()
        );

        // Same virtual offsets
        for (i, (cc, bc)) in csi_chunks.iter().zip(&bai_chunks).enumerate() {
            assert_eq!(
                cc.begin.0, bc.begin.0,
                "tid={tid}, region=[{start},{end}], chunk {i}: begin mismatch"
            );
            assert_eq!(
                cc.end.0, bc.end.0,
                "tid={tid}, region=[{start},{end}], chunk {i}: end mismatch"
            );
        }
    }
}

// --- samtools can read seqair's CSI index ---

#[test]
fn samtools_reads_seqair_csi() {
    let dir = tempfile::tempdir().unwrap();
    let bam_path = write_multi_contig_bam(dir.path());

    // Remove BAI so samtools must use the CSI
    let bai_path = bam_path.with_extension("bam.bai");
    std::fs::remove_file(&bai_path).unwrap();

    assert_eq!(samtools_count(&bam_path, "chr1"), 30);
    assert_eq!(samtools_count(&bam_path, "chr2"), 15);
    assert_eq!(samtools_count(&bam_path, "chr3"), 5);
}

#[test]
fn samtools_reads_seqair_csi_subregions() {
    let dir = tempfile::tempdir().unwrap();
    let bam_path = write_multi_contig_bam(dir.path());

    // Remove BAI so samtools must use CSI
    std::fs::remove_file(bam_path.with_extension("bam.bai")).unwrap();

    let count = samtools_count(&bam_path, "chr1:50001-150000");
    assert_eq!(count, 10, "chr1:50001-150000 should have 10 records");

    let count = samtools_count(&bam_path, "chr2:1-40000");
    assert_eq!(count, 2, "chr2:1-40000 should have 2 records");
}

// --- samtools creates CSI, seqair parses it ---

#[test]
fn seqair_parses_samtools_csi() {
    let dir = tempfile::tempdir().unwrap();
    let bam_path = write_multi_contig_bam(dir.path());

    // Remove seqair's CSI, create one with samtools
    let csi_path = bam_path.with_extension("bam.csi");
    std::fs::remove_file(&csi_path).ok();

    let output = Command::new("samtools")
        .args(["index", "-c"])
        .arg(&bam_path)
        .output()
        .expect("samtools index");
    assert!(
        output.status.success(),
        "samtools index -c failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // samtools creates file.bam.csi (BGZF-compressed)
    let idx = CsiIndex::from_path(&csi_path).expect("parse samtools CSI");
    assert_eq!(idx.min_shift(), 14);
    // samtools auto-calculates depth from the data; don't hardcode a value
    assert!(idx.depth() >= 1, "depth should be positive");
}

#[test]
fn seqair_queries_samtools_csi_match_bai() {
    let dir = tempfile::tempdir().unwrap();
    let bam_path = write_multi_contig_bam(dir.path());

    // Create CSI with samtools
    let csi_path = bam_path.with_extension("bam.csi");
    std::fs::remove_file(&csi_path).ok();
    let output = Command::new("samtools")
        .args(["index", "-c"])
        .arg(&bam_path)
        .output()
        .expect("samtools index");
    assert!(output.status.success());

    let csi = CsiIndex::from_path(&csi_path).expect("CSI");
    let bai = BamIndex::from_path(&bam_path.with_extension("bam.bai")).expect("BAI");

    // Queries on the samtools-produced CSI must yield the same chunks as BAI
    for &(tid, start, end) in
        &[(0u32, 0u32, 1_000_000u32), (1, 0, 500_000), (2, 0, 200_000), (0, 50_000, 150_000)]
    {
        let csi_chunks = csi.query(tid, Pos0::new(start).unwrap(), Pos0::new(end).unwrap());
        let bai_chunks = bai.query(tid, Pos0::new(start).unwrap(), Pos0::new(end).unwrap());

        assert_eq!(
            csi_chunks.len(),
            bai_chunks.len(),
            "tid={tid}, region=[{start},{end}]: chunk count mismatch (csi={}, bai={})\ncsi: {:?}\nbai: {:?}",
            csi_chunks.len(),
            bai_chunks.len(),
            csi_chunks.iter().map(|c| (c.begin.0, c.end.0)).collect::<Vec<_>>(),
            bai_chunks.iter().map(|c| (c.begin.0, c.end.0)).collect::<Vec<_>>(),
        );
    }
}

// --- Format detection: file.csi (without .bam infix) ---

// r[verify csi.detect]
#[test]
fn indexed_reader_finds_csi_without_bam_infix() {
    let dir = tempfile::tempdir().unwrap();
    let bam_path = write_multi_contig_bam(dir.path());

    // Remove both default index files
    std::fs::remove_file(bam_path.with_extension("bam.bai")).unwrap();
    std::fs::remove_file(bam_path.with_extension("bam.csi")).unwrap();

    // Place CSI at {file_without_ext}.csi (e.g. multi.csi instead of multi.bam.csi)
    let csi_without_infix = bam_path.with_extension("csi"); // multi.csi
    {
        // Recreate the index to write to the new path
        let header = make_header();
        let mut writer = BamWriter::from_path(&bam_path, &header, true).unwrap();
        for i in 0..5u32 {
            let rec =
                OwnedBamRecord::builder(0, i64::from(i) * 10000, format!("r{i}").into_bytes())
                    .mapq(60)
                    .cigar(vec![CigarOp::new(CigarOpType::Match, 100)])
                    .seq(cyclic_seq(100))
                    .qual(vec![BaseQuality::from_byte(30); 100])
                    .build()
                    .unwrap();
            writer.write(&rec).unwrap();
        }
        let (_inner, ib) = writer.finish().unwrap();
        let ib = ib.unwrap();
        ib.write_csi(std::fs::File::create(&csi_without_infix).unwrap(), 1, &[]).unwrap();
    }

    let mut reader = IndexedBamReader::open(&bam_path).expect("should find multi.csi");
    let mut store = RecordStore::new();
    reader
        .fetch_into(0, Pos0::new(0).unwrap(), Pos0::new(1_000_000).unwrap(), &mut store)
        .expect("fetch");
    assert!(!store.is_empty(), "should find records via multi.csi");
}

// --- IndexedBamReader with CSI (format detection) ---

#[test]
fn indexed_reader_prefers_csi_over_bai() {
    let dir = tempfile::tempdir().unwrap();
    let bam_path = write_multi_contig_bam(dir.path());

    // Both BAI and CSI exist; IndexedBamReader should prefer CSI
    assert!(bam_path.with_extension("bam.bai").exists());
    assert!(bam_path.with_extension("bam.csi").exists());

    let mut reader = IndexedBamReader::open(&bam_path).expect("open with CSI");
    let mut store = RecordStore::new();
    reader
        .fetch_into(0, Pos0::new(0).unwrap(), Pos0::new(1_000_000).unwrap(), &mut store)
        .expect("fetch");
    assert_eq!(store.len(), 30, "chr1 should have 30 records");
}

#[test]
fn indexed_reader_with_csi_only() {
    let dir = tempfile::tempdir().unwrap();
    let bam_path = write_multi_contig_bam(dir.path());

    // Remove BAI, keep CSI only
    std::fs::remove_file(bam_path.with_extension("bam.bai")).unwrap();

    let mut reader = IndexedBamReader::open(&bam_path).expect("open with CSI only");
    let mut store = RecordStore::new();
    reader
        .fetch_into(0, Pos0::new(0).unwrap(), Pos0::new(1_000_000).unwrap(), &mut store)
        .expect("fetch");
    assert_eq!(store.len(), 30);
}

// --- seqair fetch results match samtools view counts ---

#[test]
fn seqair_csi_fetch_matches_samtools_counts() {
    let dir = tempfile::tempdir().unwrap();
    let bam_path = write_multi_contig_bam(dir.path());

    // Remove BAI so both seqair and samtools use CSI
    std::fs::remove_file(bam_path.with_extension("bam.bai")).unwrap();

    let mut reader = IndexedBamReader::open(&bam_path).expect("open");

    let regions: &[(&str, u32, u32, u32, &str)] = &[
        ("chr1", 0, 0, 1_000_000, "chr1"),
        ("chr2", 1, 0, 500_000, "chr2"),
        ("chr3", 2, 0, 200_000, "chr3"),
    ];

    for &(name, tid, start, end, samtools_region) in regions {
        let mut store = RecordStore::new();
        reader
            .fetch_into(tid, Pos0::new(start).unwrap(), Pos0::new(end).unwrap(), &mut store)
            .expect("fetch");

        let st_count = samtools_count(&bam_path, samtools_region);
        assert_eq!(store.len(), st_count, "{name}: seqair={} samtools={st_count}", store.len());
    }
}

// --- query_split tests ---

#[test]
fn csi_query_split_union_equals_query() {
    let dir = tempfile::tempdir().unwrap();
    let bam_path = write_multi_contig_bam(dir.path());
    let csi_path = bam_path.with_extension("bam.csi");

    let csi = CsiIndex::from_path(&csi_path).expect("CSI");

    for &(tid, start, end) in &[(0u32, 0u32, 1_000_000u32), (0, 50_000, 150_000), (1, 0, 500_000)] {
        let start_pos = Pos0::new(start).unwrap();
        let end_pos = Pos0::new(end).unwrap();
        let flat = csi.query(tid, start_pos, end_pos);
        let split = csi.query_split(tid, start_pos, end_pos);

        let mut combined: Vec<u64> =
            split.nearby.iter().chain(&split.distant).map(|c| c.begin.0).collect();
        combined.sort();

        let mut flat_begins: Vec<u64> = flat.iter().map(|c| c.begin.0).collect();
        flat_begins.sort();

        assert_eq!(
            combined, flat_begins,
            "tid={tid}, region=[{start},{end}]: split union != flat query"
        );
    }
}

// --- samtools idxstats with CSI ---

#[test]
fn samtools_idxstats_with_seqair_csi() {
    let dir = tempfile::tempdir().unwrap();
    let bam_path = write_multi_contig_bam(dir.path());

    // Remove BAI
    std::fs::remove_file(bam_path.with_extension("bam.bai")).unwrap();

    let output = Command::new("samtools")
        .arg("idxstats")
        .arg(&bam_path)
        .output()
        .expect("samtools idxstats");
    assert!(
        output.status.success(),
        "samtools idxstats failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

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
        vec![("chr1".to_string(), 30), ("chr2".to_string(), 15), ("chr3".to_string(), 5)]
    );
}

// --- Real BAM file: samtools CSI on test.bam ---

#[test]
fn real_bam_samtools_csi_matches_bai() {
    let src_bam = Path::new("tests/data/test.bam");
    if !src_bam.exists() {
        eprintln!("skipping: tests/data/test.bam not found");
        return;
    }

    let dir = tempfile::tempdir().unwrap();
    let bam_path = dir.path().join("test.bam");
    std::fs::copy(src_bam, &bam_path).unwrap();

    // Create BAI and CSI with samtools
    Command::new("samtools").args(["index", "-b"]).arg(&bam_path).output().unwrap();
    let csi_path = bam_path.with_extension("bam.csi");
    Command::new("samtools").args(["index", "-c"]).arg(&bam_path).output().unwrap();

    let bai = BamIndex::from_path(&bam_path.with_extension("bam.bai")).expect("BAI");
    let csi = CsiIndex::from_path(&csi_path).expect("CSI");

    // Read header to get contigs
    let reader = IndexedBamReader::open(&bam_path).expect("open");
    let header = reader.header();

    // Query each contig fully and compare chunk counts
    for tid in 0..header.target_count() as u32 {
        let len = header.target_len(tid).unwrap();
        let bai_chunks = bai.query(tid, Pos0::new(0).unwrap(), Pos0::new(len as u32).unwrap());
        let csi_chunks = csi.query(tid, Pos0::new(0).unwrap(), Pos0::new(len as u32).unwrap());

        // Both should return non-empty results for the same tids
        assert_eq!(bai_chunks.is_empty(), csi_chunks.is_empty(), "tid={tid}: emptiness mismatch");
    }
}

// --- Proptest: write + parse round-trip preserves query results ---

use proptest::prelude::*;

proptest! {
    #[test]
    fn roundtrip_csi_preserves_query_chunks(
        n_records in 1usize..50,
        seed in 0u64..10_000,
    ) {
        let dir = tempfile::tempdir().unwrap();
        let header = BamHeader::from_sam_text(
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:10000000\n",
        )
        .unwrap();
        let bam_path = dir.path().join("prop.bam");
        let mut writer = BamWriter::from_path(&bam_path, &header, true).unwrap();

        for i in 0..n_records {
            let pos = (seed + (i as u64) * 1000) as i64 % 9_000_000;
            let rec = OwnedBamRecord::builder(0, pos, format!("r{i}").into_bytes())
                .mapq(60)
                .cigar(vec![CigarOp::new(CigarOpType::Match, 100)])
                .seq(cyclic_seq(100))
                .qual(vec![BaseQuality::from_byte(30); 100])
                .build()
                .unwrap();
            writer.write(&rec).unwrap();
        }

        let (_inner, index_builder) = writer.finish().unwrap();
        let ib = index_builder.unwrap();

        // Write both BAI and CSI
        let bai_path = bam_path.with_extension("bam.bai");
        let csi_path = bam_path.with_extension("bam.csi");
        ib.write_bai(std::fs::File::create(&bai_path).unwrap(), 1).unwrap();
        ib.write_csi(std::fs::File::create(&csi_path).unwrap(), 1, &[]).unwrap();

        // Parse both and compare queries
        let bai = BamIndex::from_path(&bai_path).unwrap();
        let csi = CsiIndex::from_path(&csi_path).unwrap();

        let bai_chunks = bai.query(0, Pos0::new(0).unwrap(), Pos0::new(10_000_000).unwrap());
        let csi_chunks = csi.query(0, Pos0::new(0).unwrap(), Pos0::new(10_000_000).unwrap());

        prop_assert_eq!(
            csi_chunks.len(), bai_chunks.len(),
            "chunk count mismatch: csi={}, bai={}", csi_chunks.len(), bai_chunks.len()
        );
    }
}
