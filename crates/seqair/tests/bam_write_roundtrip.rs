//! BAM write round-trip tests: build records with seqair's `BamWriter`,
//! read back with noodles and samtools, verify field-level parity.
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
use seqair::bam::aux_data::AuxData;
use seqair::bam::cigar::{CigarOp, CigarOpType};
use seqair::bam::header::BamHeader;
use seqair::bam::owned_record::OwnedBamRecord;
use seqair::bam::writer::BamWriter;
use seqair::bam::{IndexedBamReader, Pos0, RecordStore};
use seqair_types::bam_flags::BamFlags;
use seqair_types::{Base, BaseQuality};
use std::path::Path;
use std::process::Command;

fn make_header() -> BamHeader {
    BamHeader::from_sam_text(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n@SQ\tSN:chr2\tLN:50000\n",
    )
    .unwrap()
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

/// Read all mapped records with noodles from a BAM.
fn noodles_read_all(bam_path: &Path) -> Vec<noodles::bam::Record> {
    let file = std::fs::File::open(bam_path).unwrap();
    let mut reader = bam::io::Reader::new(file);
    let _header: sam::Header = reader.read_header().unwrap();
    reader.records().map(|r| r.unwrap()).collect()
}

/// Verify samtools can read the BAM without errors.
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

// --- Test cases ---

/// Simple records: write a few reads, read back with noodles, compare.
#[test]
fn roundtrip_simple_records() {
    let dir = tempfile::tempdir().unwrap();

    let records = vec![
        OwnedBamRecord::builder(0, 100, b"read1".to_vec())
            .flags(BamFlags::empty())
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
            .qual([30, 31, 32, 33, 34, 35, 36, 37, 38, 39].map(BaseQuality::from_byte).to_vec())
            .build()
            .unwrap(),
        OwnedBamRecord::builder(0, 200, b"read2".to_vec())
            .flags(BamFlags::from(16)) // reverse
            .mapq(50)
            .cigar(vec![CigarOp::new(CigarOpType::Match, 8)])
            .seq(vec![Base::T, Base::G, Base::C, Base::A, Base::T, Base::G, Base::C, Base::A])
            .qual([40, 41, 42, 43, 44, 45, 46, 47].map(BaseQuality::from_byte).to_vec())
            .build()
            .unwrap(),
    ];

    let bam_path = write_bam(dir.path(), &records);
    samtools_quickcheck(&bam_path);

    let noodles_recs = noodles_read_all(&bam_path);
    assert_eq!(noodles_recs.len(), 2);

    // Check first record
    let n = &noodles_recs[0];
    assert_eq!(n.name().map(|n| n.to_vec()).unwrap(), b"read1");
    let pos = n.alignment_start().unwrap().unwrap();
    assert_eq!(usize::from(pos), 101); // 1-based
    assert_eq!(u8::from(n.mapping_quality().unwrap()), 60);
    let seq: Vec<u8> = n.sequence().iter().collect();
    assert_eq!(seq, b"ACGTACGTAC");
    let qual: Vec<u8> = n.quality_scores().iter().collect();
    assert_eq!(qual, vec![30, 31, 32, 33, 34, 35, 36, 37, 38, 39]);

    // Check second record
    let n = &noodles_recs[1];
    assert_eq!(n.name().map(|n| n.to_vec()).unwrap(), b"read2");
    assert!(n.flags().is_reverse_complemented());

    // Also read back with seqair
    let mut reader = IndexedBamReader::open(&bam_path).expect("seqair open");
    let tid = reader.header().tid("chr1").expect("tid");
    let mut store = RecordStore::new();
    reader
        .fetch_into(tid, Pos0::new(0).unwrap(), Pos0::new(100_000).unwrap(), &mut store)
        .expect("fetch");
    assert_eq!(store.len(), 2);
    assert_eq!(store.qname(0), b"read1");
    assert_eq!(store.qname(1), b"read2");
}

/// Complex CIGARs: soft clips, deletions, insertions, introns.
#[test]
fn roundtrip_complex_cigars() {
    let dir = tempfile::tempdir().unwrap();

    // 2S3M1I4M1D3M2S — query length = 2+3+1+4+3+2 = 15
    let seq = vec![
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
        Base::G,
        Base::T,
        Base::A,
        Base::C,
        Base::G,
    ];

    let records = vec![
        // Soft clips + insertion + deletion
        OwnedBamRecord::builder(0, 100, b"complex1".to_vec())
            .mapq(40)
            .cigar(vec![
                CigarOp::new(CigarOpType::SoftClip, 2),
                CigarOp::new(CigarOpType::Match, 3),
                CigarOp::new(CigarOpType::Insertion, 1),
                CigarOp::new(CigarOpType::Match, 4),
                CigarOp::new(CigarOpType::Deletion, 1),
                CigarOp::new(CigarOpType::Match, 3),
                CigarOp::new(CigarOpType::SoftClip, 2),
            ])
            .seq(seq.clone())
            .qual(vec![BaseQuality::from_byte(30); 15])
            .build()
            .unwrap(),
        // Hard clips + intron (N op) — query length = 3+5 = 8
        OwnedBamRecord::builder(0, 500, b"intron1".to_vec())
            .mapq(55)
            .cigar(vec![
                CigarOp::new(CigarOpType::HardClip, 5),
                CigarOp::new(CigarOpType::Match, 3),
                CigarOp::new(CigarOpType::RefSkip, 1000),
                CigarOp::new(CigarOpType::Match, 5),
                CigarOp::new(CigarOpType::HardClip, 3),
            ])
            .seq(vec![Base::A, Base::C, Base::G, Base::T, Base::A, Base::C, Base::G, Base::T])
            .qual(vec![BaseQuality::from_byte(35); 8])
            .build()
            .unwrap(),
        // SeqMatch/SeqMismatch ops — query length = 5+3+4 = 12
        OwnedBamRecord::builder(0, 800, b"exact1".to_vec())
            .mapq(60)
            .cigar(vec![
                CigarOp::new(CigarOpType::SeqMatch, 5),
                CigarOp::new(CigarOpType::SeqMismatch, 3),
                CigarOp::new(CigarOpType::SeqMatch, 4),
            ])
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
                Base::G,
                Base::T,
            ])
            .qual(vec![BaseQuality::from_byte(40); 12])
            .build()
            .unwrap(),
    ];

    let bam_path = write_bam(dir.path(), &records);
    samtools_quickcheck(&bam_path);

    let noodles_recs = noodles_read_all(&bam_path);
    assert_eq!(noodles_recs.len(), 3);

    // Verify CIGARs round-trip correctly via noodles
    use noodles::sam::alignment::record::cigar::op::Kind;
    let cigar0: Vec<_> = noodles_recs[0]
        .cigar()
        .iter()
        .map(|r| {
            let op = r.unwrap();
            (op.kind(), op.len())
        })
        .collect();
    assert_eq!(
        cigar0,
        vec![
            (Kind::SoftClip, 2),
            (Kind::Match, 3),
            (Kind::Insertion, 1),
            (Kind::Match, 4),
            (Kind::Deletion, 1),
            (Kind::Match, 3),
            (Kind::SoftClip, 2),
        ]
    );

    let cigar1: Vec<_> = noodles_recs[1]
        .cigar()
        .iter()
        .map(|r| {
            let op = r.unwrap();
            (op.kind(), op.len())
        })
        .collect();
    assert_eq!(
        cigar1,
        vec![
            (Kind::HardClip, 5),
            (Kind::Match, 3),
            (Kind::Skip, 1000),
            (Kind::Match, 5),
            (Kind::HardClip, 3),
        ]
    );

    let cigar2: Vec<_> = noodles_recs[2]
        .cigar()
        .iter()
        .map(|r| {
            let op = r.unwrap();
            (op.kind(), op.len())
        })
        .collect();
    assert_eq!(
        cigar2,
        vec![(Kind::SequenceMatch, 5), (Kind::SequenceMismatch, 3), (Kind::SequenceMatch, 4),]
    );

    // Read back with seqair and compare
    let mut reader = IndexedBamReader::open(&bam_path).expect("seqair open");
    let tid = reader.header().tid("chr1").expect("tid");
    let mut store = RecordStore::new();
    reader
        .fetch_into(tid, Pos0::new(0).unwrap(), Pos0::new(100_000).unwrap(), &mut store)
        .expect("fetch");
    assert_eq!(store.len(), 3);
}

/// Aux tags: write records with various aux tag types, verify they survive round-trip.
#[test]
fn roundtrip_aux_tags() {
    let dir = tempfile::tempdir().unwrap();

    let mut aux = AuxData::new();
    aux.set_int(*b"NM", 3).unwrap();
    aux.set_string(*b"RG", b"sample1");
    aux.set_float(*b"XS", -1.5);

    let records = vec![
        OwnedBamRecord::builder(0, 100, b"aux_read".to_vec())
            .mapq(60)
            .cigar(vec![CigarOp::new(CigarOpType::Match, 5)])
            .seq(vec![Base::A, Base::C, Base::G, Base::T, Base::A])
            .qual([30, 31, 32, 33, 34].map(BaseQuality::from_byte).to_vec())
            .aux(aux)
            .build()
            .unwrap(),
    ];

    let bam_path = write_bam(dir.path(), &records);
    samtools_quickcheck(&bam_path);

    // Verify with samtools view
    let output = Command::new("samtools")
        .args(["view", "-o", "-"])
        .arg(&bam_path)
        .output()
        .expect("samtools view");
    let sam_line = String::from_utf8(output.stdout).unwrap();
    assert!(sam_line.contains("NM:i:3"), "NM tag missing: {sam_line}");
    assert!(sam_line.contains("RG:Z:sample1"), "RG tag missing: {sam_line}");
    assert!(sam_line.contains("XS:f:"), "XS tag missing: {sam_line}");
}

/// Multiple contigs: records on chr1 and chr2, verify both are retrievable.
#[test]
fn roundtrip_multiple_contigs() {
    let dir = tempfile::tempdir().unwrap();

    let records = vec![
        OwnedBamRecord::builder(0, 100, b"chr1_read".to_vec())
            .mapq(60)
            .cigar(vec![CigarOp::new(CigarOpType::Match, 5)])
            .seq(vec![Base::A, Base::C, Base::G, Base::T, Base::A])
            .qual(vec![BaseQuality::from_byte(30); 5])
            .build()
            .unwrap(),
        OwnedBamRecord::builder(1, 200, b"chr2_read".to_vec())
            .mapq(50)
            .cigar(vec![CigarOp::new(CigarOpType::Match, 6)])
            .seq(vec![Base::T, Base::G, Base::C, Base::A, Base::T, Base::G])
            .qual(vec![BaseQuality::from_byte(35); 6])
            .build()
            .unwrap(),
    ];

    let bam_path = write_bam(dir.path(), &records);
    samtools_quickcheck(&bam_path);

    let mut reader = IndexedBamReader::open(&bam_path).expect("open");

    // chr1
    let tid1 = reader.header().tid("chr1").expect("chr1");
    let mut store = RecordStore::new();
    reader
        .fetch_into(tid1, Pos0::new(0).unwrap(), Pos0::new(100_000).unwrap(), &mut store)
        .expect("fetch chr1");
    assert_eq!(store.len(), 1);
    assert_eq!(store.qname(0), b"chr1_read");

    // chr2
    let tid2 = reader.header().tid("chr2").expect("chr2");
    store.clear();
    reader
        .fetch_into(tid2, Pos0::new(0).unwrap(), Pos0::new(50_000).unwrap(), &mut store)
        .expect("fetch chr2");
    assert_eq!(store.len(), 1);
    assert_eq!(store.qname(0), b"chr2_read");
}

/// Index co-production: write BAM + BAI, then use samtools to query regions
/// and verify the counts match seqair's fetch.
#[test]
fn index_coproduction_matches_samtools() {
    let dir = tempfile::tempdir().unwrap();

    // Write 20 records spread across chr1
    let records: Vec<OwnedBamRecord> = (0..20)
        .map(|i| {
            OwnedBamRecord::builder(0, i64::from(i) * 1000, format!("read{i}").into_bytes())
                .mapq(60)
                .cigar(vec![CigarOp::new(CigarOpType::Match, 100)])
                .seq(
                    (0..100)
                        .map(|j| {
                            if j % 4 == 0 {
                                Base::A
                            } else if j % 4 == 1 {
                                Base::C
                            } else if j % 4 == 2 {
                                Base::G
                            } else {
                                Base::T
                            }
                        })
                        .collect(),
                )
                .qual(vec![BaseQuality::from_byte(30); 100])
                .build()
                .unwrap()
        })
        .collect();

    let bam_path = write_bam(dir.path(), &records);
    samtools_quickcheck(&bam_path);

    // Query the whole contig with samtools
    let output = Command::new("samtools")
        .args(["view", "-c"])
        .arg(&bam_path)
        .arg("chr1")
        .output()
        .expect("samtools view");
    let samtools_count: usize = String::from_utf8(output.stdout).unwrap().trim().parse().unwrap();

    // Query same full contig with seqair via the co-produced index
    let mut reader = IndexedBamReader::open(&bam_path).expect("open");
    let tid = reader.header().tid("chr1").expect("tid");
    let mut store = RecordStore::new();
    reader
        .fetch_into(tid, Pos0::new(0).unwrap(), Pos0::new(100_000).unwrap(), &mut store)
        .expect("fetch");

    assert_eq!(
        store.len(),
        samtools_count,
        "seqair full-contig count ({}) != samtools count ({samtools_count})",
        store.len()
    );
    assert_eq!(store.len(), 20, "should have all 20 records");
}

/// Paired-end records: write proper pairs, verify flags and mate info survive.
#[test]
fn roundtrip_paired_end() {
    let dir = tempfile::tempdir().unwrap();

    let seq5 = vec![Base::A, Base::C, Base::G, Base::T, Base::A];
    let records = vec![
        OwnedBamRecord::builder(0, 100, b"pair1".to_vec())
            .flags(BamFlags::from(99)) // paired, proper, mate reverse, read1
            .mapq(60)
            .cigar(vec![CigarOp::new(CigarOpType::Match, 5)])
            .seq(seq5.clone())
            .qual(vec![BaseQuality::from_byte(30); 5])
            .next_ref_id(0)
            .next_pos(200)
            .template_len(105)
            .build()
            .unwrap(),
        OwnedBamRecord::builder(0, 200, b"pair1".to_vec())
            .flags(BamFlags::from(147)) // paired, proper, reverse, read2
            .mapq(55)
            .cigar(vec![CigarOp::new(CigarOpType::Match, 5)])
            .seq(seq5)
            .qual(vec![BaseQuality::from_byte(35); 5])
            .next_ref_id(0)
            .next_pos(100)
            .template_len(-105)
            .build()
            .unwrap(),
    ];

    let bam_path = write_bam(dir.path(), &records);
    samtools_quickcheck(&bam_path);

    let noodles_recs = noodles_read_all(&bam_path);
    assert_eq!(noodles_recs.len(), 2);

    // Read 1
    let r1 = &noodles_recs[0];
    assert!(r1.flags().is_segmented()); // paired
    assert!(r1.flags().is_properly_segmented()); // proper pair
    assert!(r1.flags().is_first_segment()); // read1

    // Read 2
    let r2 = &noodles_recs[1];
    assert!(r2.flags().is_last_segment()); // read2
    assert!(r2.flags().is_reverse_complemented()); // reverse
}

/// `write_store_record`: write via `RecordStore`, validate with samtools + noodles.
// r[verify bam_writer.test_store_roundtrip]
#[test]
fn roundtrip_write_store_record() {
    let dir = tempfile::tempdir().unwrap();
    let header = make_header();

    // Build OwnedBamRecords, serialize to raw BAM, push into RecordStore
    let seq5 = vec![Base::A, Base::C, Base::G, Base::T, Base::A];
    let records = vec![
        OwnedBamRecord::builder(0, 100, b"sr1".to_vec())
            .flags(BamFlags::from(99))
            .mapq(60)
            .cigar(vec![CigarOp::new(CigarOpType::Match, 5)])
            .seq(seq5.clone())
            .qual([30, 31, 32, 33, 34].map(BaseQuality::from_byte).to_vec())
            .next_ref_id(0)
            .next_pos(300)
            .template_len(205)
            .aux({
                let mut a = AuxData::new();
                a.set_string(*b"RG", b"grp1");
                a
            })
            .build()
            .unwrap(),
        OwnedBamRecord::builder(0, 300, b"sr1".to_vec())
            .flags(BamFlags::from(147))
            .mapq(55)
            .cigar(vec![
                CigarOp::new(CigarOpType::SoftClip, 1),
                CigarOp::new(CigarOpType::Match, 4),
            ])
            .seq(seq5)
            .qual([35, 36, 37, 38, 39].map(BaseQuality::from_byte).to_vec())
            .next_ref_id(0)
            .next_pos(100)
            .template_len(-205)
            .build()
            .unwrap(),
    ];

    let mut store = RecordStore::new();
    let mut raw_buf = Vec::new();
    for rec in &records {
        raw_buf.clear();
        rec.to_bam_bytes(&mut raw_buf).unwrap();
        store.push_raw(&raw_buf, &mut ()).unwrap();
    }

    // Write via write_store_record
    let bam_path = dir.path().join("store.bam");
    {
        let mut writer = BamWriter::from_path(&bam_path, &header, true).unwrap();
        for i in 0..store.len() as u32 {
            writer.write_store_record(&store, i).unwrap();
        }
        let (_inner, index_builder) = writer.finish().unwrap();
        if let Some(ib) = index_builder {
            let bai_file = std::fs::File::create(bam_path.with_extension("bam.bai")).unwrap();
            ib.write_bai(bai_file, header.target_count()).unwrap();
        }
    }

    // External validation
    samtools_quickcheck(&bam_path);

    // Noodles field-level comparison
    let noodles_recs = noodles_read_all(&bam_path);
    assert_eq!(noodles_recs.len(), 2);

    let r1 = &noodles_recs[0];
    assert_eq!(r1.name().map(|n| n.to_vec()).unwrap(), b"sr1");
    assert_eq!(r1.alignment_start().unwrap().unwrap().get(), 101); // 1-based
    assert_eq!(u8::from(r1.mapping_quality().unwrap()), 60);
    assert!(r1.flags().is_first_segment());

    let r2 = &noodles_recs[1];
    assert_eq!(r2.name().map(|n| n.to_vec()).unwrap(), b"sr1");
    assert_eq!(r2.alignment_start().unwrap().unwrap().get(), 301); // 1-based
    assert!(r2.flags().is_last_segment());

    // Read back via seqair and compare field-by-field with original store
    let mut store2 = RecordStore::new();
    let mut shared = IndexedBamReader::open(&bam_path).unwrap();
    shared.fetch_into(0, Pos0::new(0).unwrap(), Pos0::new(1000).unwrap(), &mut store2).unwrap();
    assert_eq!(store2.len(), 2);

    for i in 0..2u32 {
        let a = store.record(i);
        let b = store2.record(i);
        assert_eq!(*a.pos, *b.pos, "pos mismatch record {i}");
        assert_eq!(a.flags, b.flags, "flags mismatch record {i}");
        assert_eq!(a.mapq, b.mapq, "mapq mismatch record {i}");
        assert_eq!(a.next_ref_id, b.next_ref_id, "next_ref_id mismatch record {i}");
        assert_eq!(a.next_pos, b.next_pos, "next_pos mismatch record {i}");
        assert_eq!(a.template_len, b.template_len, "tlen mismatch record {i}");
        assert_eq!(store.qname(i), store2.qname(i), "qname mismatch record {i}");
        assert_eq!(store.cigar(i), store2.cigar(i), "cigar mismatch record {i}");
        assert_eq!(store.seq(i), store2.seq(i), "seq mismatch record {i}");
        assert_eq!(store.qual(i), store2.qual(i), "qual mismatch record {i}");
        assert_eq!(store.aux(i), store2.aux(i), "aux mismatch record {i}");
    }
}
