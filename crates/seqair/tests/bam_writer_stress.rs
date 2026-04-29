//! BAM writer stress tests: maximum field sizes, error poisoning integration,
//! placed-unmapped records, and samtools validation.
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

use seqair::bam::aux_data::AuxData;
use seqair::bam::cigar::{CigarOp, CigarOpType};
use seqair::bam::header::BamHeader;
use seqair::bam::owned_record::OwnedBamRecord;
use seqair::bam::writer::{BamWriteError, BamWriter};
use seqair::bam::{IndexedBamReader, Pos0, RecordStore};
use seqair_types::bam_flags::BamFlags;
use seqair_types::{Base, BaseQuality};
use std::path::Path;
use std::process::{Command, Stdio};

fn make_header() -> BamHeader {
    BamHeader::from_sam_text("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000000\n").unwrap()
}

fn write_bam(dir: &Path, records: &[OwnedBamRecord]) -> std::path::PathBuf {
    let header = make_header();
    let bam_path = dir.join("stress.bam");
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

fn samtools_quickcheck(bam_path: &Path) {
    let output = Command::new("samtools")
        .args(["quickcheck", "-v"])
        .arg(bam_path)
        .output()
        .expect("samtools not found");
    assert!(
        output.status.success(),
        "samtools quickcheck failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
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

// --- Maximum field sizes ---

/// Record with maximum allowed qname length (254 bytes).
#[test]
fn max_qname_length() {
    let dir = tempfile::tempdir().unwrap();
    let qname = vec![b'X'; 254];
    let records = vec![
        OwnedBamRecord::builder(0, Some(Pos0::new(100).unwrap()), qname)
            .mapq(60)
            .cigar(vec![CigarOp::new(CigarOpType::Match, 10)])
            .seq(cyclic_seq(10))
            .qual(vec![BaseQuality::from_byte(30); 10])
            .build()
            .unwrap(),
    ];
    let bam_path = write_bam(dir.path(), &records);
    samtools_quickcheck(&bam_path);

    let mut reader = IndexedBamReader::open(&bam_path).expect("open");
    let tid = reader.header().tid("chr1").expect("tid");
    let mut store = RecordStore::new();
    reader
        .fetch_into(tid, Pos0::new(0).unwrap(), Pos0::new(1_000_000).unwrap(), &mut store)
        .expect("fetch");
    assert_eq!(store.len(), 1);
    assert_eq!(store.qname(0).len(), 254);
}

/// Record with a very long sequence (50,000 bases).
#[test]
fn large_sequence() {
    let dir = tempfile::tempdir().unwrap();
    let seq_len = 50_000u32;
    let records = vec![
        OwnedBamRecord::builder(0, Some(Pos0::new(100).unwrap()), b"big".to_vec())
            .mapq(60)
            .cigar(vec![CigarOp::new(CigarOpType::Match, seq_len)])
            .seq(cyclic_seq(seq_len as usize))
            .qual(vec![BaseQuality::from_byte(30); seq_len as usize])
            .build()
            .unwrap(),
    ];
    let bam_path = write_bam(dir.path(), &records);
    samtools_quickcheck(&bam_path);

    let mut reader = IndexedBamReader::open(&bam_path).expect("open");
    let tid = reader.header().tid("chr1").expect("tid");
    let mut store = RecordStore::new();
    reader
        .fetch_into(tid, Pos0::new(0).unwrap(), Pos0::new(1_000_000).unwrap(), &mut store)
        .expect("fetch");
    assert_eq!(store.len(), 1);
    assert_eq!(store.record(0).seq_len, seq_len);
}

/// Record with many aux tags (fills the aux data buffer).
#[test]
fn many_aux_tags() {
    let dir = tempfile::tempdir().unwrap();
    let mut aux = AuxData::new();
    // Write 100 integer tags: X0..X9, Y0..Y9, etc.
    for prefix in b'A'..=b'J' {
        for digit in b'0'..=b'9' {
            aux.set_int([prefix, digit], i64::from(digit) * 100).unwrap();
        }
    }

    let records = vec![
        OwnedBamRecord::builder(0, Some(Pos0::new(100).unwrap()), b"auxtags".to_vec())
            .mapq(60)
            .cigar(vec![CigarOp::new(CigarOpType::Match, 10)])
            .seq(cyclic_seq(10))
            .qual(vec![BaseQuality::from_byte(30); 10])
            .aux(aux)
            .build()
            .unwrap(),
    ];
    let bam_path = write_bam(dir.path(), &records);
    samtools_quickcheck(&bam_path);

    // Verify samtools can see the tags
    let output =
        Command::new("samtools").args(["view"]).arg(&bam_path).output().expect("samtools view");
    let sam = String::from_utf8(output.stdout).unwrap();
    assert!(sam.contains("A0:"), "A0 tag missing");
    assert!(sam.contains("J9:"), "J9 tag missing");
}

// --- Index dispatch: mapped, placed-unmapped, fully-unmapped ---

/// Placed-unmapped reads (FLAG 0x4, `ref_id` >= 0) should be included in
/// the BAI index (in the unmapped count) but not in normal fetch results.
#[test]
fn placed_unmapped_records() {
    let dir = tempfile::tempdir().unwrap();
    let header = make_header();
    let bam_path = dir.path().join("placed.bam");
    let mut writer = BamWriter::from_path(&bam_path, &header, true).unwrap();

    // Mapped record
    let mapped = OwnedBamRecord::builder(0, Some(Pos0::new(100).unwrap()), b"mapped".to_vec())
        .mapq(60)
        .cigar(vec![CigarOp::new(CigarOpType::Match, 10)])
        .seq(cyclic_seq(10))
        .qual(vec![BaseQuality::from_byte(30); 10])
        .build()
        .unwrap();
    writer.write(&mapped).unwrap();

    // Placed-unmapped: has ref_id=0 and pos but flag 0x4
    let placed =
        OwnedBamRecord::builder(0, Some(Pos0::new(150).unwrap()), b"placed_unmap".to_vec())
            .flags(BamFlags::from(4u16)) // unmapped
            .mapq(0)
            .build()
            .unwrap();
    writer.write(&placed).unwrap();

    // Another mapped
    let mapped2 = OwnedBamRecord::builder(0, Some(Pos0::new(200).unwrap()), b"mapped2".to_vec())
        .mapq(55)
        .cigar(vec![CigarOp::new(CigarOpType::Match, 10)])
        .seq(cyclic_seq(10))
        .qual(vec![BaseQuality::from_byte(30); 10])
        .build()
        .unwrap();
    writer.write(&mapped2).unwrap();

    let (_inner, index_builder) = writer.finish().unwrap();
    if let Some(ib) = index_builder {
        let bai_file = std::fs::File::create(bam_path.with_extension("bam.bai")).unwrap();
        ib.write_bai(bai_file, header.target_count()).unwrap();
    }

    samtools_quickcheck(&bam_path);

    // samtools should see all 3 records
    let output =
        Command::new("samtools").args(["view", "-c"]).arg(&bam_path).output().expect("samtools");
    let total: usize = String::from_utf8(output.stdout).unwrap().trim().parse().unwrap();
    assert_eq!(total, 3, "samtools should see 3 records total");

    // samtools idxstats should show 2 mapped + 1 unmapped
    let output = Command::new("samtools")
        .arg("idxstats")
        .arg(&bam_path)
        .output()
        .expect("samtools idxstats");
    let stats = String::from_utf8(output.stdout).unwrap();
    let chr1_line = stats.lines().find(|l| l.starts_with("chr1")).unwrap();
    let fields: Vec<&str> = chr1_line.split('\t').collect();
    let n_mapped: usize = fields[2].parse().unwrap();
    let n_unmapped: usize = fields[3].parse().unwrap();
    assert_eq!(n_mapped, 2, "chr1 should have 2 mapped");
    assert_eq!(n_unmapped, 1, "chr1 should have 1 placed-unmapped");
}

/// Fully-unmapped records (`ref_id` == -1) should not be pushed to the index.
/// Verified via samtools: the BAM has both mapped and unmapped, but indexed
/// region queries only return the mapped record.
#[test]
fn fully_unmapped_not_indexed() {
    let dir = tempfile::tempdir().unwrap();

    // Use samtools to create a BAM with both mapped and unmapped records,
    // since writing ref_id=-1 records requires careful ordering.
    let sam_text = "@HD\tVN:1.6\tSO:coordinate\n\
                    @SQ\tSN:chr1\tLN:1000000\n\
                    mapped\t0\tchr1\t101\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n\
                    unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n";
    let sam_path = dir.path().join("mixed.sam");
    std::fs::write(&sam_path, sam_text).unwrap();

    let bam_path = dir.path().join("mixed.bam");
    let status = Command::new("samtools")
        .args(["sort", "-o"])
        .arg(&bam_path)
        .arg(&sam_path)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .expect("samtools sort");
    assert!(status.success());
    let status = Command::new("samtools")
        .arg("index")
        .arg(&bam_path)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .unwrap();
    assert!(status.success());

    // samtools sees both
    let output =
        Command::new("samtools").args(["view", "-c"]).arg(&bam_path).output().expect("samtools");
    let total: usize = String::from_utf8(output.stdout).unwrap().trim().parse().unwrap();
    assert_eq!(total, 2, "samtools should see 2 records total");

    // seqair fetch should only find the mapped record
    let mut reader = IndexedBamReader::open(&bam_path).expect("open");
    let tid = reader.header().tid("chr1").expect("tid");
    let mut store = RecordStore::new();
    reader
        .fetch_into(tid, Pos0::new(0).unwrap(), Pos0::new(1_000_000).unwrap(), &mut store)
        .expect("fetch");
    assert_eq!(store.len(), 1, "fetch should only return mapped records");
    assert_eq!(store.qname(0), b"mapped");
}

// --- Error poisoning integration ---

/// After a write error, subsequent writes return Poisoned, but the BAM
/// written up to that point is still readable.
#[test]
fn poisoned_writer_partial_output_is_readable() {
    let dir = tempfile::tempdir().unwrap();
    let header = make_header();
    let bam_path = dir.path().join("poison.bam");
    let mut writer = BamWriter::from_path(&bam_path, &header, false).unwrap();

    // Write two good records
    for i in 0..2u32 {
        let rec = OwnedBamRecord::builder(
            0,
            Some(Pos0::new(i * 100).unwrap()),
            format!("ok{i}").into_bytes(),
        )
        .mapq(60)
        .cigar(vec![CigarOp::new(CigarOpType::Match, 10)])
        .seq(cyclic_seq(10))
        .qual(vec![BaseQuality::from_byte(30); 10])
        .build()
        .unwrap();
        writer.write(&rec).unwrap();
    }

    // Force an error with a qname that's too long
    let bad_rec = OwnedBamRecord {
        ref_id: 0,
        pos: Some(Pos0::new(200).unwrap()),
        mapq: 0,
        flags: BamFlags::empty(),
        next_ref_id: -1,
        next_pos: None,
        template_len: 0,
        qname: vec![b'X'; 255],
        cigar: Vec::new(),
        seq: Vec::new(),
        qual: Vec::new(),
        aux: AuxData::new(),
    };
    assert!(writer.write(&bad_rec).is_err());

    // BamWriter is now poisoned
    let good_rec = OwnedBamRecord::builder(0, Some(Pos0::new(300).unwrap()), b"after".to_vec())
        .mapq(60)
        .cigar(vec![CigarOp::new(CigarOpType::Match, 5)])
        .seq(cyclic_seq(5))
        .qual(vec![BaseQuality::from_byte(30); 5])
        .build()
        .unwrap();
    let err = writer.write(&good_rec).unwrap_err();
    assert!(matches!(err, BamWriteError::Poisoned));

    // Finish flushes what was written before the error + writes EOF marker
    let _ = writer.finish();

    // The partial BAM should be valid and contain the 2 good records
    samtools_quickcheck(&bam_path);
    let output =
        Command::new("samtools").args(["view", "-c"]).arg(&bam_path).output().expect("samtools");
    let count: usize = String::from_utf8(output.stdout).unwrap().trim().parse().unwrap();
    assert_eq!(count, 2, "partial BAM should have 2 records");
}

// --- Dense records at same position ---

/// Many records at the same position: tests BAI bin handling.
#[test]
fn dense_records_same_position() {
    let dir = tempfile::tempdir().unwrap();
    let records: Vec<OwnedBamRecord> = (0..100)
        .map(|i| {
            OwnedBamRecord::builder(0, Some(Pos0::new(1000).unwrap()), format!("r{i}").into_bytes())
                .mapq(60)
                .cigar(vec![CigarOp::new(CigarOpType::Match, 50)])
                .seq(cyclic_seq(50))
                .qual(vec![BaseQuality::from_byte(30); 50])
                .build()
                .unwrap()
        })
        .collect();

    let bam_path = write_bam(dir.path(), &records);
    samtools_quickcheck(&bam_path);

    let mut reader = IndexedBamReader::open(&bam_path).expect("open");
    let tid = reader.header().tid("chr1").expect("tid");
    let mut store = RecordStore::new();
    reader
        .fetch_into(tid, Pos0::new(0).unwrap(), Pos0::new(1_000_000).unwrap(), &mut store)
        .expect("fetch");
    assert_eq!(store.len(), 100, "should find all 100 records at same pos");
}
