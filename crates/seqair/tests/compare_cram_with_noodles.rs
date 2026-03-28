//! Compares CRAM record decoding between seqair and noodles.
//! Reads all records sequentially from the same CRAM file using both
//! Tests run against both CRAM v3.0 and v3.1.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::type_complexity,
    dead_code
)]
use noodles::cram;
use noodles::fasta;
use noodles::sam;
use noodles::sam::alignment::record::Sequence as _;
use seqair::bam::RecordStore;
use seqair::reader::Readers;
use std::path::Path;

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

struct NoodlesRecord {
    pos: i64,
    flags: u16,
    mapq: u8,
    ref_id: Option<usize>,
    seq: Vec<u8>,
    qual: Vec<u8>,
    qname: Vec<u8>,
    seq_len: usize,
}

fn read_noodles_records(cram_path: &Path) -> (sam::Header, Vec<NoodlesRecord>) {
    let fasta_reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(test_fasta_path())
        .expect("open indexed FASTA");
    let repository =
        fasta::Repository::new(fasta::repository::adapters::IndexedReader::new(fasta_reader));

    let mut reader = cram::io::reader::Builder::default()
        .set_reference_sequence_repository(repository)
        .build_from_path(cram_path)
        .expect("open CRAM");

    let header = reader.read_header().expect("read header");

    let mut records = Vec::new();
    for result in reader.records(&header) {
        let record = result.expect("read record");

        let flags = record.flags();
        let flags_bits = u16::from(flags);

        if flags.is_unmapped() {
            continue;
        }

        let pos = record.alignment_start().map(|p| usize::from(p) as i64 - 1).unwrap_or(-1);
        let mapq = record.mapping_quality().map(u8::from).unwrap_or(255);
        let ref_id = record.reference_sequence_id();
        let seq: Vec<u8> = record.sequence().iter().collect();
        let qual: Vec<u8> = record.quality_scores().iter().collect();
        let qname: Vec<u8> = record.name().map(|n| n.to_vec()).unwrap_or_default();
        let seq_len = seq.len();

        records.push(NoodlesRecord {
            pos,
            flags: flags_bits,
            mapq,
            ref_id,
            seq,
            qual,
            qname,
            seq_len,
        });
    }

    (header, records)
}

// r[verify cram.record.decode_order]
// r[verify unified.fetch_equivalence]
#[test]
fn cram_chr19_count_matches_noodles() {
    let start = 6_105_000u64;
    let end = 6_140_000u64;

    for &(version, cram_path_fn) in CRAM_VERSIONS {
        let cram_path = cram_path_fn();
        let (_header, noodles_records) = read_noodles_records(cram_path);

        let noodles_in_range = noodles_records
            .iter()
            .filter(|r| r.ref_id == Some(0) && r.pos >= start as i64 && r.pos < end as i64)
            .count();

        let mut readers = Readers::open(cram_path, test_fasta_path()).unwrap();
        let chr19_tid = readers.header().tid("chr19").unwrap();
        let mut store = RecordStore::new();
        readers.fetch_into(chr19_tid, start, end, &mut store).unwrap();

        let our_in_range =
            (0..store.len() as u32).filter(|&i| store.record(i).pos >= start as i64).count();

        assert_eq!(
            our_in_range, noodles_in_range,
            "{version}: chr19 count ours={our_in_range} noodles={noodles_in_range}",
        );
    }
}

// r[verify cram.record.position]
// r[verify cram.record.flags]
#[test]
fn cram_chr19_records_match_noodles_field_by_field() {
    let start = 6_105_000i64;
    let end = 6_140_000i64;

    for &(version, cram_path_fn) in CRAM_VERSIONS {
        let cram_path = cram_path_fn();
        let (_header, noodles_records) = read_noodles_records(cram_path);

        let noodles_chr19: Vec<&NoodlesRecord> = noodles_records
            .iter()
            .filter(|r| r.ref_id == Some(0) && r.pos >= start && r.pos < end)
            .collect();

        let mut readers = Readers::open(cram_path, test_fasta_path()).unwrap();
        let chr19_tid = readers.header().tid("chr19").unwrap();
        let mut store = RecordStore::new();
        readers.fetch_into(chr19_tid, start as u64, end as u64, &mut store).unwrap();

        let our_records: Vec<u32> =
            (0..store.len() as u32).filter(|&i| store.record(i).pos >= start).collect();

        assert_eq!(our_records.len(), noodles_chr19.len(), "{version}: count mismatch");

        for (i, noodles_rec) in noodles_chr19.iter().enumerate() {
            let our_rec = store.record(our_records[i]);
            assert_eq!(our_rec.pos, noodles_rec.pos, "{version} rec {i}: pos");
            assert_eq!(our_rec.flags, noodles_rec.flags, "{version} rec {i}: flags");
            assert_eq!(our_rec.mapq, noodles_rec.mapq, "{version} rec {i}: mapq");
        }
    }
}

// r[verify cram.record.decode_order]
#[test]
fn cram_chr19_sequences_match_noodles() {
    let start = 6_105_000i64;
    let end = 6_140_000i64;

    for &(version, cram_path_fn) in CRAM_VERSIONS {
        let cram_path = cram_path_fn();
        let (_header, noodles_records) = read_noodles_records(cram_path);

        let noodles_chr19: Vec<&NoodlesRecord> = noodles_records
            .iter()
            .filter(|r| r.ref_id == Some(0) && r.pos >= start && r.pos < end)
            .collect();

        let mut readers = Readers::open(cram_path, test_fasta_path()).unwrap();
        let chr19_tid = readers.header().tid("chr19").unwrap();
        let mut store = RecordStore::new();
        readers.fetch_into(chr19_tid, start as u64, end as u64, &mut store).unwrap();

        let our_records: Vec<u32> =
            (0..store.len() as u32).filter(|&i| store.record(i).pos >= start).collect();

        assert_eq!(our_records.len(), noodles_chr19.len(), "{version}: count mismatch");

        for (i, noodles_rec) in noodles_chr19.iter().enumerate() {
            let idx = our_records[i];
            let our_bases = store.seq(idx);
            // Only compare positions where noodles emits a standard base (A/C/G/T).
            // Noodles preserves IUPAC codes; seqair normalizes ambiguous bases to Unknown.
            for (pos, (&our_base, &noodles_byte)) in
                our_bases.iter().zip(noodles_rec.seq.iter()).enumerate()
            {
                if matches!(noodles_byte, b'A' | b'C' | b'G' | b'T') {
                    let our_byte = match our_base {
                        seqair_types::Base::A => b'A',
                        seqair_types::Base::C => b'C',
                        seqair_types::Base::G => b'G',
                        seqair_types::Base::T => b'T',
                        seqair_types::Base::Unknown => b'N',
                    };
                    assert_eq!(
                        our_byte, noodles_byte,
                        "{version} rec {i} pos {pos}: sequence mismatch"
                    );
                }
            }
        }
    }
}

// r[verify cram.record.decode_order]
#[test]
fn cram_chr19_quality_scores_match_noodles() {
    let start = 6_105_000i64;
    let end = 6_140_000i64;

    for &(version, cram_path_fn) in CRAM_VERSIONS {
        let cram_path = cram_path_fn();
        let (_header, noodles_records) = read_noodles_records(cram_path);

        let noodles_chr19: Vec<&NoodlesRecord> = noodles_records
            .iter()
            .filter(|r| r.ref_id == Some(0) && r.pos >= start && r.pos < end)
            .collect();

        let mut readers = Readers::open(cram_path, test_fasta_path()).unwrap();
        let chr19_tid = readers.header().tid("chr19").unwrap();
        let mut store = RecordStore::new();
        readers.fetch_into(chr19_tid, start as u64, end as u64, &mut store).unwrap();

        let our_records: Vec<u32> =
            (0..store.len() as u32).filter(|&i| store.record(i).pos >= start).collect();

        assert_eq!(our_records.len(), noodles_chr19.len(), "{version}: count mismatch");

        for (i, noodles_rec) in noodles_chr19.iter().enumerate() {
            let idx = our_records[i];
            let our_qual = store.qual(idx);
            assert_eq!(
                our_qual,
                noodles_rec.qual.as_slice(),
                "{version} rec {i}: quality mismatch"
            );
        }
    }
}

// r[verify cram.record.decode_order]
#[test]
fn cram_chr19_qnames_match_noodles() {
    let start = 6_105_000i64;
    let end = 6_140_000i64;

    for &(version, cram_path_fn) in CRAM_VERSIONS {
        let cram_path = cram_path_fn();
        let (_header, noodles_records) = read_noodles_records(cram_path);

        let noodles_chr19: Vec<&NoodlesRecord> = noodles_records
            .iter()
            .filter(|r| r.ref_id == Some(0) && r.pos >= start && r.pos < end)
            .collect();

        let mut readers = Readers::open(cram_path, test_fasta_path()).unwrap();
        let chr19_tid = readers.header().tid("chr19").unwrap();
        let mut store = RecordStore::new();
        readers.fetch_into(chr19_tid, start as u64, end as u64, &mut store).unwrap();

        let our_records: Vec<u32> =
            (0..store.len() as u32).filter(|&i| store.record(i).pos >= start).collect();

        assert_eq!(our_records.len(), noodles_chr19.len(), "{version}: count mismatch");

        for (i, noodles_rec) in noodles_chr19.iter().enumerate() {
            let idx = our_records[i];
            let our_qname = store.qname(idx);
            assert_eq!(
                our_qname,
                noodles_rec.qname.as_slice(),
                "{version} rec {i}: qname mismatch"
            );
        }
    }
}
