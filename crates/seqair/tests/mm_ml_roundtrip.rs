//! MM/ML tag round-trip tests: verify that base modification tags survive
//! BAM and CRAM read/write paths as raw aux data.
//!
//! MM (Z-type): encodes modification positions, e.g. "C+m,0,2,5;"
//! ML (B:C array): encodes modification probabilities, e.g. [200, 180, 230]
//!
//! seqair currently preserves these as raw aux bytes but does not parse
//! them into structured modification data. These tests validate the
//! pass-through correctness.
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

use seqair::bam::aux::AuxValue;
use seqair::bam::aux_data::AuxData;
use seqair::bam::cigar::{CigarOp, CigarOpType};
use seqair::bam::header::BamHeader;
use seqair::bam::owned_record::OwnedBamRecord;
use seqair::bam::writer::BamWriter;
use seqair::bam::{IndexedBamReader, Pos0, RecordStore};
use seqair_types::{Base, BaseQuality};
use std::path::Path;
use std::process::{Command, Stdio};

fn make_header() -> BamHeader {
    BamHeader::from_sam_text("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n").unwrap()
}

fn write_bam(dir: &Path, records: &[OwnedBamRecord]) -> std::path::PathBuf {
    let header = make_header();
    let bam_path = dir.join("test.bam");
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

/// Build a record with typical `CpG` methylation MM/ML tags.
fn make_methylation_record(
    pos: u32,
    name: &[u8],
    mm_str: &[u8],
    ml_probs: &[u8],
) -> OwnedBamRecord {
    let seq = vec![
        Base::A,
        Base::C,
        Base::G,
        Base::T,
        Base::A,
        Base::C,
        Base::G,
        Base::C,
        Base::G,
        Base::A,
        Base::C,
        Base::G,
        Base::T,
        Base::C,
        Base::G,
        Base::A,
        Base::C,
        Base::G,
        Base::T,
        Base::A,
    ];
    let mut aux = AuxData::new();
    aux.set_string(*b"MM", mm_str);
    aux.set_array_u8(*b"ML", ml_probs).unwrap();

    OwnedBamRecord::builder(0, Some(Pos0::new(pos).unwrap()), name.to_vec())
        .mapq(60)
        .cigar(vec![CigarOp::new(CigarOpType::Match, 20)])
        .seq(seq)
        .qual(vec![BaseQuality::from_byte(30); 20])
        .aux(aux)
        .build()
        .unwrap()
}

// ---- BAM round-trip ----

// r[verify base_mod.passthrough.bam]
/// MM/ML tags survive BAM write -> read round-trip via seqair.
#[test]
fn mm_ml_bam_roundtrip_seqair() {
    let dir = tempfile::tempdir().unwrap();

    let records = vec![
        make_methylation_record(100, b"read1", b"C+m,0,2,5;", &[200, 180, 230]),
        make_methylation_record(200, b"read2", b"C+m,1;C+h,0,3;", &[150, 100, 220, 190]),
    ];

    let bam_path = write_bam(dir.path(), &records);

    let mut reader = IndexedBamReader::open(&bam_path).expect("open BAM");
    let tid = reader.header().tid("chr1").expect("tid");
    let mut store = RecordStore::new();
    reader
        .fetch_into(tid, Pos0::new(0).unwrap(), Pos0::new(100_000).unwrap(), &mut store)
        .expect("fetch");

    assert_eq!(store.len(), 2);

    // Record 1: check MM tag
    let aux1 = AuxData::from_bytes(store.aux(0).to_vec());
    let mm1 = aux1.get(*b"MM").expect("MM tag missing on read1");
    match mm1 {
        AuxValue::String(s) => assert_eq!(s, b"C+m,0,2,5;"),
        other => panic!("MM should be String, got {other:?}"),
    }

    // Record 1: check ML tag
    let ml1 = aux1.get(*b"ML").expect("ML tag missing on read1");
    match ml1 {
        AuxValue::ArrayU8(bytes) => assert_eq!(bytes, &[200, 180, 230]),
        other => panic!("ML should be ArrayU8, got {other:?}"),
    }

    // Record 2: multi-modification MM string
    let aux2 = AuxData::from_bytes(store.aux(1).to_vec());
    let mm2 = aux2.get(*b"MM").expect("MM tag missing on read2");
    match mm2 {
        AuxValue::String(s) => assert_eq!(s, b"C+m,1;C+h,0,3;"),
        other => panic!("MM should be String, got {other:?}"),
    }
    let ml2 = aux2.get(*b"ML").expect("ML tag missing on read2");
    match ml2 {
        AuxValue::ArrayU8(bytes) => assert_eq!(bytes, &[150, 100, 220, 190]),
        other => panic!("ML should be ArrayU8, got {other:?}"),
    }
}

/// MM/ML tags survive BAM write -> samtools view -> check SAM text.
#[test]
fn mm_ml_bam_samtools_readable() {
    let dir = tempfile::tempdir().unwrap();

    let records = vec![make_methylation_record(100, b"read1", b"C+m,0,2,5;", &[200, 180, 230])];

    let bam_path = write_bam(dir.path(), &records);

    let output =
        Command::new("samtools").args(["view"]).arg(&bam_path).output().expect("samtools view");
    assert!(output.status.success());

    let sam_line = String::from_utf8(output.stdout).unwrap();
    assert!(
        sam_line.contains("MM:Z:C+m,0,2,5;"),
        "MM tag missing from samtools output: {sam_line}"
    );
    assert!(
        sam_line.contains("ML:B:C,200,180,230"),
        "ML tag missing from samtools output: {sam_line}"
    );
}

// ---- CRAM round-trip ----

// r[verify base_mod.passthrough.cram]
/// MM/ML tags survive BAM -> CRAM -> BAM round-trip via samtools + seqair.
#[test]
fn mm_ml_cram_roundtrip() {
    let dir = tempfile::tempdir().unwrap();

    let records = vec![
        make_methylation_record(100, b"read1", b"C+m,0,2,5;", &[200, 180, 230]),
        make_methylation_record(200, b"read2", b"C+h,1,3;", &[170, 210]),
    ];

    let bam_path = write_bam(dir.path(), &records);

    // BAM -> CRAM via samtools
    let cram_path = dir.path().join("test.cram");
    let ref_path = dir.path().join("ref.fa");
    // Create a minimal reference for CRAM encoding
    std::fs::write(&ref_path, ">chr1\n".to_string() + &"A".repeat(100_000) + "\n").unwrap();
    let status = Command::new("samtools")
        .arg("faidx")
        .arg(&ref_path)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .expect("samtools faidx");
    assert!(status.success());

    let status = Command::new("samtools")
        .args(["view", "-C", "-T"])
        .arg(&ref_path)
        .arg("-o")
        .arg(&cram_path)
        .arg(&bam_path)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .expect("samtools view -C");
    assert!(status.success());

    let status = Command::new("samtools")
        .arg("index")
        .arg(&cram_path)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .expect("samtools index");
    assert!(status.success());

    // Read CRAM with seqair
    let mut readers = seqair::reader::Readers::open(&cram_path, &ref_path).expect("open CRAM");
    let tid = readers.header().tid("chr1").expect("tid");
    let mut store = RecordStore::new();
    readers
        .fetch_into(tid, Pos0::new(0).unwrap(), Pos0::new(100_000).unwrap(), &mut store)
        .expect("fetch");

    assert_eq!(store.len(), 2);

    // Verify MM/ML survived the CRAM round-trip
    let aux1 = AuxData::from_bytes(store.aux(0).to_vec());
    let mm1 = aux1.get(*b"MM").expect("MM tag lost in CRAM round-trip");
    match mm1 {
        AuxValue::String(s) => assert_eq!(s, b"C+m,0,2,5;"),
        other => panic!("MM should be String after CRAM, got {other:?}"),
    }
    let ml1 = aux1.get(*b"ML").expect("ML tag lost in CRAM round-trip");
    match ml1 {
        AuxValue::ArrayU8(bytes) => assert_eq!(bytes, &[200, 180, 230]),
        other => panic!("ML should be ArrayU8 after CRAM, got {other:?}"),
    }

    let aux2 = AuxData::from_bytes(store.aux(1).to_vec());
    let mm2 = aux2.get(*b"MM").expect("MM tag lost on read2 in CRAM");
    match mm2 {
        AuxValue::String(s) => assert_eq!(s, b"C+h,1,3;"),
        other => panic!("MM should be String, got {other:?}"),
    }
}

// ---- Edge cases ----

/// Empty MM string (no modifications) with empty ML array.
#[test]
fn mm_ml_empty_values() {
    let dir = tempfile::tempdir().unwrap();

    let mut aux = AuxData::new();
    aux.set_string(*b"MM", b"");
    aux.set_array_u8(*b"ML", &[]).unwrap();

    let records = vec![
        OwnedBamRecord::builder(0, Some(Pos0::new(100).unwrap()), b"empty_mods".to_vec())
            .mapq(60)
            .cigar(vec![CigarOp::new(CigarOpType::Match, 10)])
            .seq(vec![
                Base::A,
                Base::C,
                Base::G,
                Base::T,
                Base::A,
                Base::C,
                Base::G,
                Base::T,
                Base::A,
                Base::C,
            ])
            .qual(vec![BaseQuality::from_byte(30); 10])
            .aux(aux)
            .build()
            .unwrap(),
    ];

    let bam_path = write_bam(dir.path(), &records);

    let mut reader = IndexedBamReader::open(&bam_path).expect("open BAM");
    let tid = reader.header().tid("chr1").expect("tid");
    let mut store = RecordStore::new();
    reader
        .fetch_into(tid, Pos0::new(0).unwrap(), Pos0::new(100_000).unwrap(), &mut store)
        .expect("fetch");

    let aux = AuxData::from_bytes(store.aux(0).to_vec());
    match aux.get(*b"MM").expect("MM tag missing") {
        AuxValue::String(s) => assert!(s.is_empty(), "empty MM should be empty string"),
        other => panic!("MM should be String, got {other:?}"),
    }
    match aux.get(*b"ML").expect("ML tag missing") {
        AuxValue::ArrayU8(bytes) => assert!(bytes.is_empty(), "empty ML should be empty array"),
        other => panic!("ML should be ArrayU8, got {other:?}"),
    }
}

/// Large MM/ML tags (simulating long-read nanopore data with many
/// modification calls) must survive round-trip without truncation.
#[test]
fn mm_ml_large_tags() {
    let dir = tempfile::tempdir().unwrap();

    // Build a long MM string: C+m with 500 position deltas
    let mm_str: String = {
        let mut s = String::from("C+m");
        for i in 0..500u32 {
            s.push(',');
            s.push_str(&(i % 10).to_string());
        }
        s.push(';');
        s
    };
    let ml_probs: Vec<u8> = (0..500).map(|i| (i % 256) as u8).collect();

    // Need a sequence long enough for the CIGAR
    let seq: Vec<Base> = (0..1000)
        .map(|i| match i % 4 {
            0 => Base::A,
            1 => Base::C,
            2 => Base::G,
            _ => Base::T,
        })
        .collect();

    let mut aux = AuxData::new();
    aux.set_string(*b"MM", mm_str.as_bytes());
    aux.set_array_u8(*b"ML", &ml_probs).unwrap();

    let records = vec![
        OwnedBamRecord::builder(0, Some(Pos0::new(100).unwrap()), b"long_read".to_vec())
            .mapq(60)
            .cigar(vec![CigarOp::new(CigarOpType::Match, 1000)])
            .seq(seq)
            .qual(vec![BaseQuality::from_byte(30); 1000])
            .aux(aux)
            .build()
            .unwrap(),
    ];

    let bam_path = write_bam(dir.path(), &records);

    let mut reader = IndexedBamReader::open(&bam_path).expect("open BAM");
    let tid = reader.header().tid("chr1").expect("tid");
    let mut store = RecordStore::new();
    reader
        .fetch_into(tid, Pos0::new(0).unwrap(), Pos0::new(100_000).unwrap(), &mut store)
        .expect("fetch");

    let aux = AuxData::from_bytes(store.aux(0).to_vec());

    let mm = aux.get(*b"MM").expect("MM tag missing");
    match mm {
        AuxValue::String(s) => {
            assert_eq!(s, mm_str.as_bytes(), "large MM string was truncated or corrupted");
        }
        other => panic!("MM should be String, got {other:?}"),
    }

    let ml = aux.get(*b"ML").expect("ML tag missing");
    match ml {
        AuxValue::ArrayU8(bytes) => {
            assert_eq!(bytes.len(), 500, "ML array length mismatch");
            assert_eq!(bytes, ml_probs.as_slice(), "large ML array was corrupted");
        }
        other => panic!("ML should be ArrayU8, got {other:?}"),
    }
}

// r[verify base_mod.passthrough.coexistence]
/// MM/ML tags coexist with other aux tags without interference.
#[test]
fn mm_ml_coexist_with_other_tags() {
    let dir = tempfile::tempdir().unwrap();

    let mut aux = AuxData::new();
    aux.set_int(*b"NM", 2).unwrap();
    aux.set_string(*b"RG", b"sample1");
    aux.set_string(*b"MM", b"C+m,0,5;");
    aux.set_array_u8(*b"ML", &[200, 180]).unwrap();
    aux.set_float(*b"XS", 0.95);
    aux.set_string(*b"MD", b"5A4^GT10");

    let records = vec![
        OwnedBamRecord::builder(0, Some(Pos0::new(100).unwrap()), b"multi_aux".to_vec())
            .mapq(60)
            .cigar(vec![CigarOp::new(CigarOpType::Match, 10)])
            .seq(vec![
                Base::A,
                Base::C,
                Base::G,
                Base::T,
                Base::A,
                Base::C,
                Base::G,
                Base::T,
                Base::A,
                Base::C,
            ])
            .qual(vec![BaseQuality::from_byte(30); 10])
            .aux(aux)
            .build()
            .unwrap(),
    ];

    let bam_path = write_bam(dir.path(), &records);

    let mut reader = IndexedBamReader::open(&bam_path).expect("open BAM");
    let tid = reader.header().tid("chr1").expect("tid");
    let mut store = RecordStore::new();
    reader
        .fetch_into(tid, Pos0::new(0).unwrap(), Pos0::new(100_000).unwrap(), &mut store)
        .expect("fetch");

    let aux = AuxData::from_bytes(store.aux(0).to_vec());

    // All tags must survive
    assert!(matches!(aux.get(*b"NM"), Some(AuxValue::U8(2))));
    assert!(matches!(aux.get(*b"RG"), Some(AuxValue::String(b"sample1"))));
    assert!(matches!(aux.get(*b"MM"), Some(AuxValue::String(b"C+m,0,5;"))));
    assert!(matches!(aux.get(*b"XS"), Some(AuxValue::Float(_))));
    assert!(matches!(aux.get(*b"MD"), Some(AuxValue::String(b"5A4^GT10"))));

    match aux.get(*b"ML").expect("ML missing") {
        AuxValue::ArrayU8(bytes) => assert_eq!(bytes, &[200, 180]),
        other => panic!("ML should be ArrayU8, got {other:?}"),
    }
}
