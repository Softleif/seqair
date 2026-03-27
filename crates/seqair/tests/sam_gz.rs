//! Tests for reading bgzf-compressed indexed SAM files (.sam.gz).
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::indexing_slicing)]
use seqair::bam::RecordStore;
use std::path::Path;
use std::process::Command;
use tempfile::TempDir;

fn test_bam_path() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.bam"))
}

const CONTIGS: &[(&str, u64, u64)] = &[
    ("chr19", 6_103_076, 6_143_229),
    ("2kb_3_Unmodified", 1, 2_018),
    ("bacteriophage_lambda_CpG", 1, 48_502),
];

fn create_sam_gz(dir: &Path) -> std::path::PathBuf {
    let sam_gz = dir.join("test.sam.gz");
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

struct RecordSnapshot {
    pos: i64,
    end_pos: i64,
    flags: u16,
    mapq: u8,
    qname: Vec<u8>,
    seq_len: u32,
}

fn snapshot_store(store: &RecordStore) -> Vec<RecordSnapshot> {
    (0..store.len() as u32)
        .map(|i| {
            let r = store.record(i);
            RecordSnapshot {
                pos: r.pos,
                end_pos: r.end_pos,
                flags: r.flags,
                mapq: r.mapq,
                qname: store.qname(i).to_vec(),
                seq_len: r.seq_len,
            }
        })
        .collect()
}

// r[verify sam.reader.open]
// r[verify sam.reader.fetch_into]
// r[verify sam.reader.overlap_filter+2]
// r[verify sam.perf.bulk_read]
#[test]
fn sam_gz_record_count_matches_bam() {
    let dir = TempDir::new().expect("create temp dir");
    let sam_gz = create_sam_gz(dir.path());

    let mut bam_reader = seqair::bam::IndexedBamReader::open(test_bam_path()).expect("bam");
    let mut sam_reader = seqair::sam::reader::IndexedSamReader::open(&sam_gz).expect("sam");

    for &(contig, start, end) in CONTIGS {
        let tid = bam_reader.header().tid(contig).expect("tid");

        let mut bam_store = RecordStore::new();
        bam_reader.fetch_into(tid, start, end, &mut bam_store).expect("bam fetch");

        let sam_tid = sam_reader.header().tid(contig).expect("sam tid");
        let mut sam_store = RecordStore::new();
        sam_reader.fetch_into(sam_tid, start, end, &mut sam_store).expect("sam fetch");

        assert_eq!(
            sam_store.len(),
            bam_store.len(),
            "{contig}: record count mismatch sam_gz={} bam={}",
            sam_store.len(),
            bam_store.len(),
        );
    }
}

// r[verify sam.record.parse]
// r[verify sam.record.coordinate_conversion]
// r[verify sam.record.cigar_parse]
// r[verify sam.record.seq_decode]
// r[verify sam.record.qual_decode]
#[test]
fn sam_gz_record_fields_match_bam() {
    let dir = TempDir::new().expect("create temp dir");
    let sam_gz = create_sam_gz(dir.path());

    let mut bam_reader = seqair::bam::IndexedBamReader::open(test_bam_path()).expect("bam");
    let mut sam_reader = seqair::sam::reader::IndexedSamReader::open(&sam_gz).expect("sam");

    for &(contig, start, end) in CONTIGS {
        let tid = bam_reader.header().tid(contig).expect("tid");

        let mut bam_store = RecordStore::new();
        bam_reader.fetch_into(tid, start, end, &mut bam_store).expect("bam fetch");

        let sam_tid = sam_reader.header().tid(contig).expect("sam tid");
        let mut sam_store = RecordStore::new();
        sam_reader.fetch_into(sam_tid, start, end, &mut sam_store).expect("sam fetch");

        let bam_recs = snapshot_store(&bam_store);
        let sam_recs = snapshot_store(&sam_store);

        assert_eq!(sam_recs.len(), bam_recs.len(), "{contig}: count mismatch");

        for (i, (sam, bam)) in sam_recs.iter().zip(&bam_recs).enumerate() {
            assert_eq!(sam.pos, bam.pos, "{contig} rec {i}: pos");
            assert_eq!(sam.end_pos, bam.end_pos, "{contig} rec {i}: end_pos");
            assert_eq!(sam.flags, bam.flags, "{contig} rec {i}: flags");
            assert_eq!(sam.mapq, bam.mapq, "{contig} rec {i}: mapq");
            assert_eq!(sam.qname, bam.qname, "{contig} rec {i}: qname");
            assert_eq!(sam.seq_len, bam.seq_len, "{contig} rec {i}: seq_len");
        }
    }
}

// r[verify sam.record.seq_decode]
// r[verify sam.record.qual_decode]
#[test]
fn sam_gz_sequence_and_quality_match_bam() {
    let dir = TempDir::new().expect("create temp dir");
    let sam_gz = create_sam_gz(dir.path());

    let mut bam_reader = seqair::bam::IndexedBamReader::open(test_bam_path()).expect("bam");
    let mut sam_reader = seqair::sam::reader::IndexedSamReader::open(&sam_gz).expect("sam");

    let contig = "chr19";
    let start = 6_105_700;
    let end = 6_105_800;

    let tid = bam_reader.header().tid(contig).expect("tid");
    let mut bam_store = RecordStore::new();
    bam_reader.fetch_into(tid, start, end, &mut bam_store).expect("bam fetch");

    let sam_tid = sam_reader.header().tid(contig).expect("sam tid");
    let mut sam_store = RecordStore::new();
    sam_reader.fetch_into(sam_tid, start, end, &mut sam_store).expect("sam fetch");

    assert_eq!(sam_store.len(), bam_store.len());

    for i in 0..bam_store.len() as u32 {
        assert_eq!(sam_store.seq(i), bam_store.seq(i), "rec {i}: seq mismatch");
        assert_eq!(sam_store.qual(i), bam_store.qual(i), "rec {i}: qual mismatch");
    }
}

// r[verify sam.header.parse]
// r[verify sam.header.preserve_full_text]
#[test]
fn sam_gz_header_matches_bam() {
    let dir = TempDir::new().expect("create temp dir");
    let sam_gz = create_sam_gz(dir.path());

    let bam_header = seqair::bam::BamHeader::from_bam_path(test_bam_path()).expect("bam header");
    let sam_reader = seqair::sam::reader::IndexedSamReader::open(&sam_gz).expect("sam");
    let sam_header = sam_reader.header();

    assert_eq!(bam_header.target_count(), sam_header.target_count());
    for tid in 0..bam_header.target_count() as u32 {
        assert_eq!(bam_header.target_name(tid), sam_header.target_name(tid));
        assert_eq!(bam_header.target_len(tid), sam_header.target_len(tid));
    }
}
