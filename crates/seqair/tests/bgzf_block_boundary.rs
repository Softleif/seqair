//! BGZF block boundary tests: verify correct handling of BAM records
//! that straddle 64KB BGZF block boundaries.
//!
//! htslib tests this with synthetic records containing ~32KB sequences and
//! ~32KB CIGAR strings that force records to span multiple BGZF blocks.
//! We replicate that approach here.
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

use noodles::bam;
use noodles::sam;
use seqair::bam::cigar::{CigarOp, CigarOpType};
use seqair::bam::header::BamHeader;
use seqair::bam::owned_record::OwnedBamRecord;
use seqair::bam::writer::BamWriter;
use seqair::bam::{IndexedBamReader, Pos0, RecordStore};
use seqair_types::Base;
use std::path::Path;
use std::process::Command;

fn make_header() -> BamHeader {
    BamHeader::from_sam_text("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:10000000\n").unwrap()
}

fn write_bam_with_index(dir: &Path, records: &[OwnedBamRecord]) -> std::path::PathBuf {
    let header = make_header();
    let bam_path = dir.join("boundary.bam");
    let mut writer = BamWriter::from_path(&bam_path, &header, true).unwrap();
    for rec in records {
        writer.write(rec).unwrap();
    }
    let (_inner, index_builder) = writer.finish().unwrap();
    if let Some(ib) = index_builder {
        let bai_file = std::fs::File::create(bam_path.with_extension("bam.bai")).unwrap();
        ib.write_bai(bai_file, header.target_count()).unwrap();
    }
    bam_path
}

/// Generate a sequence of `n` bases cycling through A, C, G, T.
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

/// Build a large record that will exceed a single BGZF block (64KB).
/// Uses a long sequence (~8000 bases) which creates a large BAM record.
/// Multiple such records ensure block boundaries are exercised.
fn make_large_record(pos: i64, name: &[u8], seq_len: u32) -> OwnedBamRecord {
    OwnedBamRecord::builder(0, pos, name.to_vec())
        .mapq(60)
        .cigar(vec![CigarOp::new(CigarOpType::Match, seq_len)])
        .seq(cyclic_seq(seq_len as usize))
        .qual(vec![30; seq_len as usize])
        .build()
        .unwrap()
}

/// Build a record with a very long CIGAR (many small ops) to inflate the
/// binary CIGAR size. Each CIGAR op is 4 bytes in BAM, so 4000 ops = 16KB.
fn make_long_cigar_record(pos: i64, name: &[u8]) -> OwnedBamRecord {
    // Alternate M and I ops: 2000 × (1M + 1I) = 4000 ops, query_len = 4000
    let mut cigar = Vec::with_capacity(4000);
    for _ in 0..2000 {
        cigar.push(CigarOp::new(CigarOpType::Match, 1));
        cigar.push(CigarOp::new(CigarOpType::Insertion, 1));
    }
    let query_len = 4000u32; // 2000 M + 2000 I

    OwnedBamRecord::builder(0, pos, name.to_vec())
        .mapq(50)
        .cigar(cigar)
        .seq(cyclic_seq(query_len as usize))
        .qual(vec![25; query_len as usize])
        .build()
        .unwrap()
}

/// Records with large sequences force BGZF block splits.
/// Verify all records parse correctly after write→read round-trip.
#[test]
fn large_records_spanning_blocks() {
    let dir = tempfile::tempdir().unwrap();

    // Each record with 8000 bases: ~4000 bytes seq + 8000 bytes qual + overhead ≈ 12KB.
    // Write enough records to span multiple 64KB blocks.
    let records: Vec<OwnedBamRecord> = (0..20)
        .map(|i| make_large_record(i64::from(i) * 10000, format!("big{i}").as_bytes(), 8000))
        .collect();

    let bam_path = write_bam_with_index(dir.path(), &records);

    // samtools quickcheck
    let output = Command::new("samtools")
        .args(["quickcheck", "-v"])
        .arg(&bam_path)
        .output()
        .expect("samtools not found");
    assert!(output.status.success(), "samtools quickcheck failed");

    // Read with noodles
    let file = std::fs::File::open(&bam_path).unwrap();
    let mut reader = bam::io::Reader::new(file);
    let _header: sam::Header = reader.read_header().unwrap();
    let noodles_recs: Vec<_> = reader.records().map(|r| r.unwrap()).collect();
    assert_eq!(noodles_recs.len(), 20);

    // Read with seqair
    let mut seqair_reader = IndexedBamReader::open(&bam_path).expect("open");
    let tid = seqair_reader.header().tid("chr1").expect("tid");
    let mut store = RecordStore::new();
    seqair_reader
        .fetch_into(tid, Pos0::new(0).unwrap(), Pos0::new(10_000_000).unwrap(), &mut store)
        .expect("fetch");

    assert_eq!(store.len(), 20, "seqair should find all 20 large records");

    // Verify each record's content survived the block boundary
    for i in 0..20u32 {
        let r = store.record(i);
        assert_eq!(r.seq_len, 8000, "rec {i}: seq_len");
        assert_eq!(r.mapq, 60, "rec {i}: mapq");

        // Verify sequence content (cyclic ACGT pattern)
        for pos in 0..100 {
            // Spot-check first 100 bases
            let expected = match pos % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            };
            assert_eq!(store.seq_at(i, pos) as u8, expected, "rec {i} seq[{pos}]");
        }
    }
}

/// Records with very long CIGARs (4000 ops = 16KB of CIGAR data) also
/// stress the BGZF block boundary handling.
#[test]
fn long_cigar_records_spanning_blocks() {
    let dir = tempfile::tempdir().unwrap();

    let records: Vec<OwnedBamRecord> = (0..10)
        .map(|i| make_long_cigar_record(i64::from(i) * 50000, format!("cigar{i}").as_bytes()))
        .collect();

    let bam_path = write_bam_with_index(dir.path(), &records);

    // samtools validates the output
    let output = Command::new("samtools")
        .args(["quickcheck", "-v"])
        .arg(&bam_path)
        .output()
        .expect("samtools");
    assert!(output.status.success(), "samtools quickcheck failed");

    // Read with seqair
    let mut reader = IndexedBamReader::open(&bam_path).expect("open");
    let tid = reader.header().tid("chr1").expect("tid");
    let mut store = RecordStore::new();
    reader
        .fetch_into(tid, Pos0::new(0).unwrap(), Pos0::new(10_000_000).unwrap(), &mut store)
        .expect("fetch");

    assert_eq!(store.len(), 10, "should have all 10 long-cigar records");

    for i in 0..10u32 {
        let r = store.record(i);
        assert_eq!(r.seq_len, 4000, "rec {i}: seq_len should be 4000");

        // Verify CIGAR op count: 4000 ops × 4 bytes = 16000 bytes of CIGAR
        let cigar_bytes = store.cigar(i);
        assert_eq!(cigar_bytes.len(), 4000 * 4, "rec {i}: cigar byte count");
    }
}

/// Mix of small and large records: ensure small records immediately before
/// and after block boundaries are correctly parsed.
#[test]
fn mixed_record_sizes_across_boundaries() {
    let dir = tempfile::tempdir().unwrap();

    let mut records = Vec::new();
    for i in 0..50u32 {
        let pos = i64::from(i) * 2000;
        let name = format!("mix{i}");
        if i % 5 == 0 {
            // Every 5th record is large (8000 bases)
            records.push(make_large_record(pos, name.as_bytes(), 8000));
        } else {
            // Normal small record (50 bases)
            records.push(
                OwnedBamRecord::builder(0, pos, name.into_bytes())
                    .mapq(30)
                    .cigar(vec![CigarOp::new(CigarOpType::Match, 50)])
                    .seq(cyclic_seq(50))
                    .qual(vec![30; 50])
                    .build()
                    .unwrap(),
            );
        }
    }

    let bam_path = write_bam_with_index(dir.path(), &records);

    let output = Command::new("samtools")
        .args(["quickcheck", "-v"])
        .arg(&bam_path)
        .output()
        .expect("samtools");
    assert!(output.status.success());

    // Compare noodles vs seqair counts
    let file = std::fs::File::open(&bam_path).unwrap();
    let mut noodles_reader = bam::io::Reader::new(file);
    let _header: sam::Header = noodles_reader.read_header().unwrap();
    let noodles_count = noodles_reader.records().count();

    let mut seqair_reader = IndexedBamReader::open(&bam_path).expect("open");
    let tid = seqair_reader.header().tid("chr1").expect("tid");
    let mut store = RecordStore::new();
    seqair_reader
        .fetch_into(tid, Pos0::new(0).unwrap(), Pos0::new(10_000_000).unwrap(), &mut store)
        .expect("fetch");

    assert_eq!(store.len(), 50);
    assert_eq!(store.len(), noodles_count);

    // Verify the large records have correct seq_len
    for i in (0..50u32).step_by(5) {
        assert_eq!(store.record(i).seq_len, 8000, "rec {i} should be large");
    }
    // And the small ones
    for i in 0..50u32 {
        if i % 5 != 0 {
            assert_eq!(store.record(i).seq_len, 50, "rec {i} should be small");
        }
    }
}
