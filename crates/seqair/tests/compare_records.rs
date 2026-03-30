//! Comparison tests: BAM record decoding.
//! Fetches the same region from a BAM file using both htslib and seqair,
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::indexing_slicing)]

use rust_htslib::bam::{self, Read as _};
use seqair::bam::{Pos, Zero};
use std::path::Path;

fn test_bam_path() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.bam"))
}

const TEST_REGION: &str = "chr19";
const TEST_START: u64 = 6_105_700;
const TEST_END: u64 = 6_105_800;

struct HtsRecord {
    pos: i64,
    end_pos: i64,
    flags: u16,
    mapq: u8,
    qname: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>,
    seq_len: usize,
}

fn fetch_htslib_records() -> Vec<HtsRecord> {
    let bam_path = test_bam_path();
    let mut reader = bam::IndexedReader::from_path(bam_path).expect("htslib open");
    let tid = reader.header().tid(TEST_REGION.as_bytes()).expect("tid");
    reader.fetch((tid, TEST_START as i64, TEST_END as i64)).expect("htslib fetch");

    let mut records = Vec::new();
    let mut record = bam::Record::new();
    while let Some(Ok(())) = reader.read(&mut record) {
        // Skip unmapped reads to match seqair behavior
        if record.flags() & 0x4 != 0 {
            continue;
        }
        let seq: Vec<u8> = (0..record.seq_len()).map(|i| record.seq()[i]).collect();
        record.cache_cigar();
        let end_pos = record.cigar_cached().unwrap().end_pos();

        records.push(HtsRecord {
            pos: record.pos(),
            end_pos,
            flags: record.flags(),
            mapq: record.mapq(),

            qname: record.qname().to_vec(),
            seq,
            qual: record.qual().to_vec(),
            seq_len: record.seq_len(),
        });
    }
    records
}

// r[verify bam.reader.fetch_into+2]
// r[verify bam.reader.overlap_filter]
// r[verify bam.reader.sorted_order+2]
// r[verify bam.reader.unmapped_skipped+2]
// r[verify arena.fixed_size_records]
// r[verify arena.no_per_record_alloc]
// r[verify arena.record_ref]
// r[verify arena.index_access]
// r[verify arena.iteration]
// r[verify bam.index.bai_magic]
// r[verify bam.index.bai_parse]
// r[verify bam.index.region_query]
// r[verify bam.index.bin_calculation]
// r[verify bgzf.seek]
// r[verify bgzf.virtual_offset]
// r[verify bgzf.eof]
// r[verify bgzf.read_partial]
#[test]
fn record_count_matches() {
    let hts_records = fetch_htslib_records();

    let bam_path = test_bam_path();
    let mut reader = seqair::bam::IndexedBamReader::open(bam_path).expect("seqair open");
    let mut store = seqair::bam::RecordStore::new();
    let tid = reader.header().tid(TEST_REGION).expect("tid");
    reader
        .fetch_into(
            tid,
            Pos::<Zero>::new(TEST_START as u32),
            Pos::<Zero>::new(TEST_END as u32),
            &mut store,
        )
        .expect("seqair fetch");

    assert_eq!(
        store.len(),
        hts_records.len(),
        "record count mismatch: seqair={} htslib={}",
        store.len(),
        hts_records.len()
    );
}

// r[verify bam.record.decode]
// r[verify bam.record.fields]
// r[verify bam.record.end_pos]
#[test]
fn record_fields_match() {
    let hts_records = fetch_htslib_records();

    let bam_path = test_bam_path();
    let mut reader = seqair::bam::IndexedBamReader::open(bam_path).expect("seqair open");
    let mut store = seqair::bam::RecordStore::new();
    let tid = reader.header().tid(TEST_REGION).expect("tid");
    reader
        .fetch_into(
            tid,
            Pos::<Zero>::new(TEST_START as u32),
            Pos::<Zero>::new(TEST_END as u32),
            &mut store,
        )
        .expect("seqair fetch");

    for (i, hts) in hts_records.iter().enumerate() {
        let rio = store.record(i as u32);

        assert_eq!(rio.pos.as_i64(), hts.pos, "pos mismatch at record {i}");
        assert_eq!(rio.flags, hts.flags, "flags mismatch at record {i}");
        assert_eq!(rio.mapq, hts.mapq, "mapq mismatch at record {i}");
        assert_eq!(store.qname(i as u32), hts.qname.as_slice(), "qname mismatch at record {i}");
        assert_eq!(store.qual(i as u32), hts.qual.as_slice(), "qual mismatch at record {i}");

        // htslib's CigarStringView::end_pos() returns pos + ref_consumed (exclusive)
        // Our end_pos = pos + ref_consumed - 1 (inclusive)
        assert_eq!(
            rio.end_pos.as_i64(),
            hts.end_pos - 1,
            "end_pos mismatch at record {i}: rio={} hts_exclusive={}",
            rio.end_pos.as_i64(),
            hts.end_pos
        );
    }
}

// r[verify bam.record.seq_4bit]
// r[verify bam.record.seq_at]
// r[verify bam.record.seq_at_simd+2]
#[test]
fn sequence_matches() {
    let hts_records = fetch_htslib_records();

    let bam_path = test_bam_path();
    let mut reader = seqair::bam::IndexedBamReader::open(bam_path).expect("seqair open");
    let mut store = seqair::bam::RecordStore::new();
    let tid = reader.header().tid(TEST_REGION).expect("tid");
    reader
        .fetch_into(
            tid,
            Pos::<Zero>::new(TEST_START as u32),
            Pos::<Zero>::new(TEST_END as u32),
            &mut store,
        )
        .expect("seqair fetch");

    for (i, hts) in hts_records.iter().enumerate() {
        let rio = store.record(i as u32);
        assert_eq!(rio.seq_len as usize, hts.seq_len, "seq_len mismatch at record {i}");

        for pos in 0..hts.seq_len {
            let base = store.seq_at(i as u32, pos);
            let hts_base = seqair_types::Base::from(hts.seq[pos]);
            assert_eq!(
                base, hts_base,
                "seq mismatch at record {i} pos {pos}: rio={base} hts={hts_base}",
            );
        }
    }
}

// r[verify bam.record.flag_reverse]
// r[verify bam.record.flag_first]
// r[verify bam.record.flag_second]
// r[verify bam.record.flag_unmapped]
#[test]
fn flag_helpers_match() {
    let hts_records = fetch_htslib_records();

    let bam_path = test_bam_path();
    let mut reader = seqair::bam::IndexedBamReader::open(bam_path).expect("seqair open");
    let mut store = seqair::bam::RecordStore::new();
    let tid = reader.header().tid(TEST_REGION).expect("tid");
    reader
        .fetch_into(
            tid,
            Pos::<Zero>::new(TEST_START as u32),
            Pos::<Zero>::new(TEST_END as u32),
            &mut store,
        )
        .expect("seqair fetch");

    for (i, hts) in hts_records.iter().enumerate() {
        let rio = store.record(i as u32);
        assert_eq!(
            rio.flags & 0x10 != 0,
            hts.flags & 0x10 != 0,
            "is_reverse mismatch at record {i}"
        );
        assert_eq!(rio.flags & 0x40 != 0, hts.flags & 0x40 != 0, "is_first mismatch at record {i}");
        assert_eq!(
            rio.flags & 0x80 != 0,
            hts.flags & 0x80 != 0,
            "is_second mismatch at record {i}"
        );
        assert_eq!(
            rio.flags & 0x4 != 0,
            hts.flags & 0x4 != 0,
            "is_unmapped mismatch at record {i}"
        );
    }
}

// r[verify bam.record.aux_parse]
// r[verify bam.record.raw_aux]
#[test]
fn aux_tags_accessible() {
    let bam_path = test_bam_path();
    let mut reader = seqair::bam::IndexedBamReader::open(bam_path).expect("seqair open");
    let mut store = seqair::bam::RecordStore::new();
    let tid = reader.header().tid(TEST_REGION).expect("tid");
    reader
        .fetch_into(
            tid,
            Pos::<Zero>::new(TEST_START as u32),
            Pos::<Zero>::new(TEST_END as u32),
            &mut store,
        )
        .expect("seqair fetch");

    let mut found_any_tag = false;
    for i in 0..store.len() {
        if !store.aux(i as u32).is_empty() {
            found_any_tag = true;
            break;
        }
    }
    assert!(found_any_tag, "expected at least one record with aux tags");
}
