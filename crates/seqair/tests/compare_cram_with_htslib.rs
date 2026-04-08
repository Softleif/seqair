//! Full record-level comparison of CRAM decoding against htslib.
//! Fetches the same regions from BAM (via htslib) and CRAM (via seqair),
//! CRAM v3.0 and v3.1 to ensure codec correctness across versions.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::type_complexity,
    reason = "test code"
)]
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    reason = "test code with known small values"
)]
use rust_htslib::bam::{self, FetchDefinition, Read as _};
use seqair::bam::{Pos, RecordStore, Zero};
use seqair::reader::Readers;
use std::path::Path;

fn test_bam_path() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.bam"))
}

fn test_cram_v30_path() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test_v30.cram"))
}

fn test_cram_v31_path() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.cram"))
}

fn test_fasta_path() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.fasta.gz"))
}

const CRAM_VERSIONS: &[(&str, fn() -> &'static Path)] = &[
    ("v3.0", test_cram_v30_path as fn() -> &'static Path),
    ("v3.1", test_cram_v31_path as fn() -> &'static Path),
];

struct HtsRecord {
    pos: i64,
    end_pos: i64,
    flags: u16,
    mapq: u8,
    seq_len: usize,
}

fn fetch_hts_records(contig: &str, start: u64, end: u64) -> Vec<HtsRecord> {
    let mut reader = bam::IndexedReader::from_path(test_bam_path()).expect("htslib open");
    let tid = reader.header().tid(contig.as_bytes()).expect("tid");
    reader.fetch(FetchDefinition::Region(tid as i32, start as i64, end as i64)).expect("fetch");

    let mut records = Vec::new();
    let mut record = bam::Record::new();
    while reader.read(&mut record) == Some(Ok(())) {
        if record.flags() & 0x4 != 0 {
            continue;
        }
        record.cache_cigar();
        let end_pos = record.cigar_cached().unwrap().end_pos();
        records.push(HtsRecord {
            pos: record.pos(),
            end_pos,
            flags: record.flags(),
            mapq: record.mapq(),
            seq_len: record.seq_len(),
        });
    }
    records
}

fn fetch_hts_quals(contig: &str, start: u64, end: u64) -> Vec<Vec<u8>> {
    let mut reader = bam::IndexedReader::from_path(test_bam_path()).unwrap();
    let tid = reader.header().tid(contig.as_bytes()).unwrap();
    reader.fetch(FetchDefinition::Region(tid as i32, start as i64, end as i64)).unwrap();

    let mut quals = Vec::new();
    let mut record = bam::Record::new();
    while reader.read(&mut record) == Some(Ok(())) {
        if record.flags() & 0x4 != 0 {
            continue;
        }
        quals.push(record.qual().to_vec());
    }
    quals
}

fn fetch_hts_seqs(contig: &str, start: u64, end: u64) -> Vec<Vec<u8>> {
    let mut reader = bam::IndexedReader::from_path(test_bam_path()).unwrap();
    let tid = reader.header().tid(contig.as_bytes()).unwrap();
    reader.fetch(FetchDefinition::Region(tid as i32, start as i64, end as i64)).unwrap();

    let mut seqs = Vec::new();
    let mut record = bam::Record::new();
    while reader.read(&mut record) == Some(Ok(())) {
        if record.flags() & 0x4 != 0 {
            continue;
        }
        seqs.push((0..record.seq_len()).map(|i| record.seq()[i]).collect());
    }
    seqs
}

// r[verify cram.record.decode_order]
// r[verify cram.record.flags]
// r[verify cram.record.position]
// r[verify cram.record.cigar_reconstruction]
// r[verify cram.record.read_length]
// r[verify cram.record.read_name]
// r[verify cram.slice.header]
// r[verify cram.slice.multi_ref]
// r[verify cram.record.features]
// r[verify cram.scope.versions]
// r[verify cram.container.region_skip]
// r[verify unified.fetch_equivalence]
// r[verify region_buf.not_cram]
#[test]
fn cram_records_match_htslib_for_chr19() {
    let contig = "chr19";
    let start = 6_103_076u64;
    let end = 6_143_229u64;
    let hts_records = fetch_hts_records(contig, start, end);

    for &(version, cram_path_fn) in CRAM_VERSIONS {
        let mut readers = Readers::open(cram_path_fn(), test_fasta_path()).unwrap();
        let mut store = RecordStore::new();
        let tid = readers.header().tid(contig).unwrap();
        let cram_count = readers
            .fetch_into(
                tid,
                Pos::<Zero>::new(start as u32).unwrap(),
                Pos::<Zero>::new(end as u32).unwrap(),
                &mut store,
            )
            .unwrap();

        assert_eq!(
            cram_count,
            hts_records.len(),
            "{version}: record count mismatch: cram={cram_count}, htslib={}",
            hts_records.len()
        );

        for (i, hts) in hts_records.iter().enumerate() {
            let cram = store.record(i as u32);
            assert_eq!(cram.pos.as_i64(), hts.pos, "{version} rec {i}: pos");
            assert_eq!(cram.flags, hts.flags, "{version} rec {i}: flags");
            assert_eq!(cram.mapq, hts.mapq, "{version} rec {i}: mapq");
            assert_eq!(cram.seq_len as usize, hts.seq_len, "{version} rec {i}: seq_len");
            assert_eq!(cram.end_pos.as_i64(), hts.end_pos, "{version} rec {i}: end_pos");
        }
    }
}

// r[verify cram.record.quality]
// r[verify cram.codec.rans4x8]
#[test]
fn cram_quality_scores_match_htslib() {
    let contig = "chr19";
    let start = 6_105_700u64;
    let end = 6_105_800u64;
    let hts_quals = fetch_hts_quals(contig, start, end);

    for &(version, cram_path_fn) in CRAM_VERSIONS {
        let mut readers = Readers::open(cram_path_fn(), test_fasta_path()).unwrap();
        let mut store = RecordStore::new();
        let tid = readers.header().tid(contig).unwrap();
        readers
            .fetch_into(
                tid,
                Pos::<Zero>::new(start as u32).unwrap(),
                Pos::<Zero>::new(end as u32).unwrap(),
                &mut store,
            )
            .unwrap();

        assert_eq!(store.len(), hts_quals.len(), "{version}: record count");

        for (i, hts_qual) in hts_quals.iter().enumerate() {
            let cram_qual = store.qual(i as u32);
            assert_eq!(cram_qual, hts_qual.as_slice(), "{version} rec {i}: quality mismatch");
        }
    }
}

// r[verify cram.record.sequence]
// r[verify cram.record.read_group]
// r[verify cram.record.aux_tags]
// r[verify cram.record.rg_tag]
// r[verify cram.record.mate_detached]
// r[verify cram.record.mate_attached]
// r[verify cram.edge.position_overflow]
// r[verify cram.edge.unmapped_reads]
#[test]
fn cram_sequences_match_htslib() {
    let contig = "chr19";
    let start = 6_105_700u64;
    let end = 6_105_800u64;
    let hts_seqs = fetch_hts_seqs(contig, start, end);

    for &(version, cram_path_fn) in CRAM_VERSIONS {
        let mut readers = Readers::open(cram_path_fn(), test_fasta_path()).unwrap();
        let mut store = RecordStore::new();
        let tid = readers.header().tid(contig).unwrap();
        readers
            .fetch_into(
                tid,
                Pos::<Zero>::new(start as u32).unwrap(),
                Pos::<Zero>::new(end as u32).unwrap(),
                &mut store,
            )
            .unwrap();

        assert_eq!(store.len(), hts_seqs.len(), "{version}: record count");

        for (i, hts_seq) in hts_seqs.iter().enumerate() {
            let cram_bases: Vec<u8> = (0..store.record(i as u32).seq_len as usize)
                .map(|pos| {
                    let base = store.seq_at(i as u32, pos);
                    match base {
                        seqair_types::Base::A => b'A',
                        seqair_types::Base::C => b'C',
                        seqair_types::Base::G => b'G',
                        seqair_types::Base::T => b'T',
                        seqair_types::Base::Unknown => b'N',
                    }
                })
                .collect();

            assert_eq!(
                cram_bases,
                *hts_seq,
                "{version} rec {i}: sequence mismatch at pos {}",
                cram_bases
                    .iter()
                    .zip(hts_seq.iter())
                    .position(|(a, b)| a != b)
                    .unwrap_or(cram_bases.len())
            );
        }
    }
}

// r[verify unified.readers_fork]
// r[verify unified.fork_cram]
#[test]
fn cram_fork_produces_same_records() {
    for &(version, cram_path_fn) in CRAM_VERSIONS {
        let readers = Readers::open(cram_path_fn(), test_fasta_path()).unwrap();
        let mut fork1 = readers.fork().unwrap();
        let mut fork2 = readers.fork().unwrap();

        let tid = fork1.header().tid("chr19").unwrap();
        let start = 6_105_700u64;
        let end = 6_105_800u64;

        let mut store1 = RecordStore::new();
        let mut store2 = RecordStore::new();

        fork1
            .fetch_into(
                tid,
                Pos::<Zero>::new(start as u32).unwrap(),
                Pos::<Zero>::new(end as u32).unwrap(),
                &mut store1,
            )
            .unwrap();
        fork2
            .fetch_into(
                tid,
                Pos::<Zero>::new(start as u32).unwrap(),
                Pos::<Zero>::new(end as u32).unwrap(),
                &mut store2,
            )
            .unwrap();

        assert_eq!(store1.len(), store2.len(), "{version}: forked count mismatch");
        for i in 0..store1.len() as u32 {
            assert_eq!(store1.record(i).pos, store2.record(i).pos, "{version} rec {i}: pos");
            assert_eq!(store1.record(i).flags, store2.record(i).flags, "{version} rec {i}: flags");
        }
    }
}
