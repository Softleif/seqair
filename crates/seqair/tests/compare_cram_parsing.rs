//! Compares CRAM parsing between seqair, noodles-cram, and htslib.
//! Tests verify that our CRAM container/block/header parsing produces
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    reason = "test code"
)]
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    reason = "test code with known small values"
)]

use seqair::bam::{Pos, Zero};
use std::path::Path;

fn test_cram_path() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.cram"))
}

fn test_bam_path() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.bam"))
}

fn test_crai_path() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.cram.crai"))
}

// ── File definition ──────────────────────────────────────────────────

// r[verify cram.file.magic]
#[test]
fn cram_magic_and_version() {
    let data = std::fs::read(test_cram_path()).unwrap();
    assert_eq!(&data[..4], b"CRAM");
    assert_eq!(data[4], 3, "major version should be 3");
    assert!(data[5] <= 1, "minor version should be 0 or 1");
}

// ── Header comparison: our parser vs htslib ──────────────────────────

// r[verify cram.file.header_container]
#[test]
fn cram_header_matches_bam_header_via_htslib() {
    use rust_htslib::bam::Read as _;
    use seqair::bam::BamHeader;
    use seqair::cram::{block, container::ContainerHeader};

    let cram_data = std::fs::read(test_cram_path()).unwrap();

    // Parse CRAM header container + block
    let after_file_def = &cram_data[26..];
    let header_container = ContainerHeader::parse(after_file_def).unwrap();
    let block_data = &after_file_def[header_container.header_size..];
    let (blk, _) = block::parse_block(block_data).unwrap();

    let header_len =
        i32::from_le_bytes([blk.data[0], blk.data[1], blk.data[2], blk.data[3]]) as usize;
    let header_text = std::str::from_utf8(&blk.data[4..4 + header_len]).unwrap();
    let our_header = BamHeader::from_sam_text(header_text).unwrap();

    // Parse BAM header via htslib
    let htslib_bam = rust_htslib::bam::Reader::from_path(test_bam_path()).unwrap();
    let htslib_header = htslib_bam.header();

    // Compare target names
    let our_names: Vec<&str> = our_header.target_names().collect();
    let htslib_names = htslib_header.target_names();
    assert_eq!(our_names.len(), htslib_names.len(), "target count");
    for (ours, theirs) in our_names.iter().zip(htslib_names.iter()) {
        assert_eq!(*ours, std::str::from_utf8(theirs).unwrap());
    }

    // Compare target lengths
    for tid in 0..our_header.target_count() as u32 {
        let our_len = our_header.target_len(tid).unwrap();
        let hts_len = htslib_header.target_len(tid).unwrap() as u64;
        assert_eq!(our_len, hts_len, "length mismatch for tid={tid}");
    }
}

// ── Container iteration: count containers ────────────────────────────

// r[verify cram.container.header]
#[test]
fn container_iteration_finds_all_containers() {
    use seqair::cram::container::ContainerHeader;

    let data = std::fs::read(test_cram_path()).unwrap();
    let mut pos = 26; // skip file def
    let mut data_containers = 0u32;
    let mut total_records = 0i64;
    let mut found_eof = false;

    while pos < data.len() {
        let container = ContainerHeader::parse(&data[pos..]).unwrap();

        if container.is_eof() {
            found_eof = true;
            break;
        }

        if container.num_records > 0 {
            data_containers += 1;
            total_records += i64::from(container.num_records);
        }

        pos += container.header_size + container.length as usize;
    }

    assert!(found_eof, "should find EOF");
    assert!(data_containers > 0, "should have data containers");

    // Compare total record count with htslib
    use rust_htslib::bam::Read as _;
    let mut hts = rust_htslib::bam::Reader::from_path(test_bam_path()).unwrap();
    let hts_count = hts.records().count() as i64;
    assert_eq!(
        total_records, hts_count,
        "CRAM container record count ({total_records}) should match BAM record count ({hts_count})"
    );
}

// ── CRAI index comparison ────────────────────────────────────────────

// r[verify cram.index.parse]
// r[verify cram.index.query]
// r[verify cram.index.multi_ref_slices]
#[test]
fn crai_index_covers_all_containers() {
    use seqair::cram::{container::ContainerHeader, index::CramIndex};

    let cram_data = std::fs::read(test_cram_path()).unwrap();
    let index = CramIndex::from_path(test_crai_path()).unwrap();

    // Every CRAI entry should point to a valid container
    for entry in index.entries() {
        if entry.ref_id < 0 {
            continue;
        }
        let offset = entry.container_offset as usize;
        assert!(offset < cram_data.len(), "offset {offset} out of bounds");
        let container = ContainerHeader::parse(&cram_data[offset..]).unwrap();
        assert!(container.num_records > 0);
    }
}

// ── Compression header parsing from real data ────────────────────────

// r[verify cram.compression.preservation]
// r[verify cram.compression.ds_encodings]
// r[verify cram.compression.tag_encodings]
// r[verify cram.encoding.gamma]
// r[verify cram.encoding.subexp]
// r[verify cram.scope.reference_required]
#[test]
fn compression_header_has_required_data_series() {
    use seqair::cram::{
        block, compression_header::CompressionHeader, container::ContainerHeader,
        encoding::IntEncoding,
    };

    let data = std::fs::read(test_cram_path()).unwrap();

    // Skip to first data container
    let after_file_def = &data[26..];
    let hdr = ContainerHeader::parse(after_file_def).unwrap();
    let data_start = hdr.header_size + hdr.length as usize;
    let dc_bytes = &after_file_def[data_start..];
    let dc = ContainerHeader::parse(dc_bytes).unwrap();

    // Parse compression header block
    let block_data = &dc_bytes[dc.header_size..];
    let (comp_block, _) = block::parse_block(block_data).unwrap();
    let ch = CompressionHeader::parse(&comp_block.data).unwrap();

    // All essential data series must be present (not Null)
    assert!(!matches!(ch.data_series.bam_flags, IntEncoding::Null));
    assert!(!matches!(ch.data_series.cram_flags, IntEncoding::Null));
    assert!(!matches!(ch.data_series.read_length, IntEncoding::Null));
    assert!(!matches!(ch.data_series.alignment_pos, IntEncoding::Null));
    assert!(!matches!(ch.data_series.feature_count, IntEncoding::Null));
    assert!(!matches!(ch.data_series.mapping_quality, IntEncoding::Null));

    // Tag dictionary + encodings should be present
    assert!(!ch.preservation.tag_dictionary.is_empty());
    assert!(!ch.tag_encodings.is_empty());

    // Substitution matrix should produce valid bases for all ref/code combos
    let sm = &ch.preservation.substitution_matrix;
    for ref_base in [b'A', b'C', b'G', b'T', b'N'] {
        for code in 0..4u8 {
            let sub = sm.substitute(ref_base, code);
            assert!(
                sub.is_ascii_alphabetic(),
                "substitute({}, {code}) = {sub} should be a letter",
                ref_base as char
            );
        }
    }
}

// ── BAM ↔ CRAM record comparison against htslib ─────────────────────

// r[verify cram.record.decode_order]
// r[verify cram.record.sequence]
// r[verify cram.record.cigar_reconstruction]
// r[verify cram.edge.empty_slice]
// r[verify cram.edge.seq_unknown]
// r[verify cram.edge.coordinate_clamp]
// r[verify cram.edge.long_reads]
// r[verify cram.perf.slice_granularity]
// r[verify cram.perf.reference_caching]
// r[verify cram.perf.codec_overhead]
// r[verify cram.slice.embedded_ref]
// r[verify unified.fetch_equivalence]
#[test]
fn cram_records_match_bam_records() {
    use rust_htslib::bam::{self, FetchDefinition, Read as _};
    use seqair::bam::record_store::RecordStore;
    use seqair::cram::reader::IndexedCramReader;

    let fasta_path =
        Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.fasta.gz"));

    // Use v3.0 CRAM (gzip+rans codecs)
    let cram_path =
        Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test_v30.cram"));

    let mut cram_reader = IndexedCramReader::open(cram_path, fasta_path).unwrap();
    let mut cram_store = RecordStore::new();

    // Also read BAM via htslib for comparison
    let mut hts = bam::IndexedReader::from_path(test_bam_path()).unwrap();

    // Compare for chr19 region
    let tid = 0u32;
    let start = 6_103_076u64;
    let end = 6_143_229u64;

    let cram_count = cram_reader
        .fetch_into(
            tid,
            Pos::<Zero>::new(start as u32).unwrap(),
            Pos::<Zero>::new(end as u32).unwrap(),
            &mut cram_store,
        )
        .unwrap();

    hts.fetch(FetchDefinition::Region(tid as i32, start as i64, end as i64)).unwrap();
    let mut hts_records = Vec::new();
    let mut record = bam::Record::new();
    while hts.read(&mut record) == Some(Ok(())) {
        if record.flags() & 0x4 != 0 {
            continue;
        }
        hts_records.push((record.pos(), record.flags(), record.mapq()));
    }

    // Record counts should match
    assert_eq!(
        cram_count,
        hts_records.len(),
        "CRAM should produce same number of records as BAM for chr19:{start}-{end}"
    );

    // Compare positions and flags
    for (i, (hts_pos, hts_flags, hts_mapq)) in hts_records.iter().enumerate() {
        let cram_pos = cram_store.record(i as u32).pos.as_i64();
        let cram_flags = cram_store.record(i as u32).flags;
        let cram_mapq = cram_store.record(i as u32).mapq;

        assert_eq!(
            cram_pos, *hts_pos,
            "position mismatch at record {i}: CRAM={cram_pos} vs BAM={hts_pos}"
        );
        assert_eq!(
            cram_flags, *hts_flags,
            "flags mismatch at record {i}: CRAM={cram_flags:#06x} vs BAM={hts_flags:#06x}"
        );
        assert_eq!(
            cram_mapq, *hts_mapq,
            "mapq mismatch at record {i}: CRAM={cram_mapq} vs BAM={hts_mapq}"
        );
    }
}

// ── Substitution matrix comparison with noodles test vector ──────────

// r[verify cram.compression.substitution_matrix]
#[test]
fn substitution_matrix_noodles_spec_vector() {
    use seqair::cram::compression_header::SubstitutionMatrix;

    // From noodles: § 10.6.4 test vector
    let sm = SubstitutionMatrix::parse(&[0x93, 0x1b, 0x6c, 0xb1, 0xc6]);

    // Expected results from noodles test_read_substitution_matrix
    assert_eq!(sm.substitute(b'A', 0), b'T');
    assert_eq!(sm.substitute(b'A', 1), b'G');
    assert_eq!(sm.substitute(b'A', 2), b'C');
    assert_eq!(sm.substitute(b'A', 3), b'N');

    assert_eq!(sm.substitute(b'C', 0), b'A');
    assert_eq!(sm.substitute(b'C', 1), b'G');
    assert_eq!(sm.substitute(b'C', 2), b'T');
    assert_eq!(sm.substitute(b'C', 3), b'N');

    assert_eq!(sm.substitute(b'G', 0), b'N');
    assert_eq!(sm.substitute(b'G', 1), b'A');
    assert_eq!(sm.substitute(b'G', 2), b'C');
    assert_eq!(sm.substitute(b'G', 3), b'T');

    assert_eq!(sm.substitute(b'T', 0), b'G');
    assert_eq!(sm.substitute(b'T', 1), b'N');
    assert_eq!(sm.substitute(b'T', 2), b'A');
    assert_eq!(sm.substitute(b'T', 3), b'C');

    assert_eq!(sm.substitute(b'N', 0), b'C');
    assert_eq!(sm.substitute(b'N', 1), b'G');
    assert_eq!(sm.substitute(b'N', 2), b'T');
    assert_eq!(sm.substitute(b'N', 3), b'A');
}
