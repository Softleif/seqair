//! Property-based tests for SAM parsing internals.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    reason = "test code"
)]
#![allow(clippy::arithmetic_side_effects, reason = "test code")]
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::format_collect,
    reason = "test code with known small values"
)]
mod helpers;
use hegel::prelude::*;
use helpers::ri;
use seqair::bam::Pos0;
use seqair_types::Base;

// We can't directly call the private parse functions, but we can test the
// public API with generated inputs. For CIGAR, we generate valid CIGAR strings,
// create SAM records, push through the reader, and verify round-trip properties.

/// Generate a valid CIGAR as a list of `(length, op_char)` parts plus the
/// assembled string.  Returning the parts lets callers compute reference-
/// consuming totals directly from the generator, without reimplementing
/// the parser.
#[hegel::composite]
fn arb_cigar_parts(tc: &TestCase) -> (Vec<(u32, char)>, String) {
    let parts: Vec<(u32, char)> = tc.draw_silent(
        gs::vecs(gs::tuples!(
            gs::integers::<u32>().min_value(1).max_value(199),
            gs::sampled_from(&['M', 'I', 'D', 'S', 'N', '=', 'X']),
        ))
        .min_size(1)
        .max_size(9),
    );
    let s = parts.iter().map(|(len, op)| format!("{len}{op}")).collect::<String>();
    (parts, s)
}

// r[verify sam.record.cigar_parse]
// r[verify sam.edge.long_cigar]
#[hegel::test(test_cases = 200)]
fn cigar_end_pos_invariants(tc: TestCase) {
    let (parts, cigar) = tc.draw(arb_cigar_parts());
    let pos = tc.draw(gs::integers::<i64>().min_value(100).max_value(9999));
    // Compute query length and ref-consuming total directly from the
    // generated parts — no reimplementation of the parser is needed.
    let query_consuming: u32 = parts
        .iter()
        .filter(|(_, op)| matches!(op, 'M' | 'I' | 'S' | '=' | 'X'))
        .map(|(len, _)| len)
        .sum();
    if query_consuming == 0 {
        return;
    }
    let ref_consuming: u32 = parts
        .iter()
        .filter(|(_, op)| matches!(op, 'M' | 'D' | 'N' | '=' | 'X'))
        .map(|(len, _)| len)
        .sum();

    let seq: String = (0..query_consuming).map(|_| 'A').collect();
    let qual: String = (0..query_consuming).map(|_| 'I').collect();
    let sam_text = format!(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000000\n\
         read1\t99\tchr1\t{pos_1}\t60\t{cigar}\t=\t200\t150\t{seq}\t{qual}\n",
        pos_1 = pos + 1, // SAM is 1-based
    );

    let dir = tempfile::tempdir().expect("tempdir");
    let sam_gz = create_indexed_sam(dir.path(), &sam_text);
    let mut reader =
        seqair::sam::reader::IndexedSamReader::open(&sam_gz).expect("open indexed SAM");
    let tid = reader.header().tid("chr1").unwrap();
    let mut store = seqair::bam::RecordStore::new();
    reader
        .fetch_into(tid, Pos0::new(1).unwrap(), Pos0::new(100_000_000).unwrap(), &mut store)
        .expect("fetch");

    assert!(!store.is_empty(), "should fetch at least 1 record");
    let rec = store.record(ri(0)).unwrap();

    // Invariant 1: alignment spans at least one base.
    assert!(rec.end_pos >= rec.pos, "end_pos must be >= pos");

    // Invariant 2: the ref span equals the sum of reference-consuming ops.
    // Special case: zero-refspan reads (pure insertion/soft-clip) have
    // end_pos == pos per r[pileup.zero_refspan_reads], so span is 1.
    let span = rec.end_pos.as_i64() - rec.pos.as_i64() + 1;
    let expected_span = if ref_consuming == 0 { 1 } else { i64::from(ref_consuming) };
    assert_eq!(
        span, expected_span,
        "cigar={}: ref span {} != expected {}",
        cigar, span, expected_span
    );

    // Invariant 3: seq_len is positive for a non-empty CIGAR.
    assert!(rec.seq_len > 0, "seq_len must be positive");
}

// r[verify sam.record.seq_decode]
//
// Extends coverage beyond ACGT to include all IUPAC ambiguity codes and N,
// verifying that the BAM spec mapping (M/R/W/S/Y/K/V/H/D/B/N → Unknown) is
// correctly applied end-to-end through the SAM→BAM encoding pipeline.
#[hegel::test(test_cases = 200)]
fn seq_decode_maps_iupac_to_unknown(tc: TestCase) {
    // ACGT, the two-fold and three-fold IUPAC ambiguity codes, and N.
    const IUPAC: &[u8] = b"ACGTMRWSYKVHDBN";
    let bases = tc.draw(gs::vecs(gs::sampled_from(IUPAC)).min_size(1).max_size(199));
    let seq_str: String = bases.iter().map(|&b| b as char).collect();
    let qual_str: String = bases.iter().map(|_| 'I').collect();
    let cigar = format!("{}M", bases.len());
    let sam_text = format!(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000000\n\
         read1\t99\tchr1\t100\t60\t{cigar}\t=\t200\t150\t{seq_str}\t{qual_str}\n"
    );

    let dir = tempfile::tempdir().expect("tempdir");
    let sam_gz = create_indexed_sam(dir.path(), &sam_text);
    let mut reader =
        seqair::sam::reader::IndexedSamReader::open(&sam_gz).expect("open indexed SAM");
    let tid = reader.header().tid("chr1").unwrap();
    let mut store = seqair::bam::RecordStore::new();
    reader
        .fetch_into(tid, Pos0::new(1).unwrap(), Pos0::new(100_000_000).unwrap(), &mut store)
        .expect("fetch");

    assert_eq!(store.len(), 1);
    let seq = store.record(ri(0)).unwrap().seq();
    assert_eq!(seq.len(), bases.len());
    for (i, (&input_byte, &actual)) in bases.iter().zip(seq.iter()).enumerate() {
        let expected_base = match input_byte {
            b'A' => Base::A,
            b'C' => Base::C,
            b'G' => Base::G,
            b'T' => Base::T,
            // All IUPAC ambiguity codes and N must decode to Unknown.
            _ => Base::Unknown,
        };
        assert_eq!(
            actual, expected_base,
            "pos {}: input '{}' expected {:?}",
            i, input_byte as char, expected_base
        );
    }
}

// r[verify sam.record.qual_decode]
#[hegel::test(test_cases = 200)]
fn qual_decode_roundtrip(tc: TestCase) {
    let quals = tc.draw(
        gs::vecs(
            // `*` is SAM's unknown-quality sentinel, so its Phred value is out.
            gs::integers::<u8>().max_value(93).filter(|&q| q + 33 != b'*'),
        )
        .min_size(1)
        .max_size(199),
    );
    let seq_str: String = quals.iter().map(|_| 'A').collect();
    let qual_str: String = quals.iter().map(|&q| (q + 33) as char).collect();
    let cigar = format!("{}M", quals.len());
    let sam_text = format!(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000000\n\
         read1\t99\tchr1\t100\t60\t{cigar}\t=\t200\t150\t{seq_str}\t{qual_str}\n"
    );

    let dir = tempfile::tempdir().expect("tempdir");
    let sam_gz = create_indexed_sam(dir.path(), &sam_text);
    let mut reader =
        seqair::sam::reader::IndexedSamReader::open(&sam_gz).expect("open indexed SAM");
    let tid = reader.header().tid("chr1").unwrap();
    let mut store = seqair::bam::RecordStore::new();
    reader
        .fetch_into(tid, Pos0::new(1).unwrap(), Pos0::new(100_000_000).unwrap(), &mut store)
        .expect("fetch");

    assert_eq!(store.len(), 1);
    let stored_qual = store.record(ri(0)).unwrap().qual();
    assert_eq!(stored_qual.len(), quals.len());
    for (i, (&expected, &actual)) in quals.iter().zip(stored_qual.iter()).enumerate() {
        assert_eq!(actual.as_byte(), expected, "pos {}", i);
    }
}

fn create_indexed_sam(dir: &std::path::Path, sam_text: &str) -> std::path::PathBuf {
    use std::process::{Command, Stdio};

    let sam_path = dir.join("test.sam");
    let sam_gz_path = dir.join("test.sam.gz");

    std::fs::write(&sam_path, sam_text).expect("write SAM");

    let status = Command::new("bgzip")
        .arg("-c")
        .arg(&sam_path)
        .stdout(std::fs::File::create(&sam_gz_path).expect("create .sam.gz"))
        .stderr(Stdio::null())
        .status()
        .expect("bgzip not found");
    assert!(status.success(), "bgzip failed");

    let status = Command::new("tabix")
        .args(["-p", "sam"])
        .arg(&sam_gz_path)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .expect("tabix not found");
    assert!(status.success(), "tabix failed");

    sam_gz_path
}
