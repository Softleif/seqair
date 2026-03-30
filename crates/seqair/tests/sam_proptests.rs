//! Property-based tests for SAM parsing internals.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::indexing_slicing)]
use proptest::prelude::*;
use seqair::bam::{Pos, Zero};
use seqair_types::Base;

// We can't directly call the private parse functions, but we can test the
// public API with generated inputs. For CIGAR, we generate valid CIGAR strings,
// create SAM records, push through the reader, and verify round-trip properties.

/// Generate a valid CIGAR as a list of `(length, op_char)` parts plus the
/// assembled string.  Returning the parts lets callers compute reference-
/// consuming totals directly from the generator, without reimplementing
/// the parser.
fn arb_cigar_parts() -> impl Strategy<Value = (Vec<(u32, char)>, String)> {
    prop::collection::vec(
        (1..200u32, prop::sample::select(vec!['M', 'I', 'D', 'S', 'N', '=', 'X'])),
        1..10,
    )
    .prop_map(|parts| {
        let s = parts.iter().map(|(len, op)| format!("{len}{op}")).collect::<String>();
        (parts, s)
    })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(200))]

    // r[verify sam.record.cigar_parse]
    // r[verify sam.edge.long_cigar]
    #[test]
    fn cigar_end_pos_invariants(
        (parts, cigar) in arb_cigar_parts(),
        pos in 100i64..10000,
    ) {
        // Compute query length and ref-consuming total directly from the
        // generated parts — no reimplementation of the parser is needed.
        let query_consuming: u32 = parts.iter()
            .filter(|(_, op)| matches!(op, 'M' | 'I' | 'S' | '=' | 'X'))
            .map(|(len, _)| len)
            .sum();
        if query_consuming == 0 {
            return Ok(());
        }
        let ref_consuming: u32 = parts.iter()
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

        let dir = tempfile::tempdir()?;
        let sam_gz = create_indexed_sam_proptest(dir.path(), &sam_text);
        let mut reader = seqair::sam::reader::IndexedSamReader::open(&sam_gz)?;
        let tid = reader.header().tid("chr1").unwrap();
        let mut store = seqair::bam::RecordStore::new();
        reader.fetch_into(tid, Pos::<Zero>::new(1), Pos::<Zero>::new(100_000_000), &mut store)?;

        prop_assert!(!store.is_empty(), "should fetch at least 1 record");
        let rec = store.record(0);

        // Invariant 1: alignment spans at least one base.
        prop_assert!(rec.end_pos >= rec.pos, "end_pos must be >= pos");

        // Invariant 2: the ref span equals the sum of reference-consuming ops.
        // Special case: zero-refspan reads (pure insertion/soft-clip) have
        // end_pos == pos per r[pileup.zero_refspan_reads], so span is 1.
        let span = rec.end_pos.as_i64() - rec.pos.as_i64() + 1;
        let expected_span = if ref_consuming == 0 { 1 } else { i64::from(ref_consuming) };
        prop_assert_eq!(
            span,
            expected_span,
            "cigar={}: ref span {} != expected {}",
            cigar, span, expected_span
        );

        // Invariant 3: seq_len is positive for a non-empty CIGAR.
        prop_assert!(rec.seq_len > 0, "seq_len must be positive");
    }

    // r[verify sam.record.seq_decode]
    //
    // Extends coverage beyond ACGT to include all IUPAC ambiguity codes and N,
    // verifying that the BAM spec mapping (M/R/W/S/Y/K/V/H/D/B/N → Unknown) is
    // correctly applied end-to-end through the SAM→BAM encoding pipeline.
    #[test]
    fn seq_decode_maps_iupac_to_unknown(
        bases in prop::collection::vec(
            prop::sample::select(vec![
                b'A', b'C', b'G', b'T',     // definite bases
                b'M', b'R', b'W', b'S', b'Y', b'K', // IUPAC ambiguity (2-fold)
                b'V', b'H', b'D', b'B',     // IUPAC ambiguity (3-fold)
                b'N',                        // fully unknown
            ]),
            1..200,
        ),
    ) {
        let seq_str: String = bases.iter().map(|&b| b as char).collect();
        let qual_str: String = bases.iter().map(|_| 'I').collect();
        let cigar = format!("{}M", bases.len());
        let sam_text = format!(
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000000\n\
             read1\t99\tchr1\t100\t60\t{cigar}\t=\t200\t150\t{seq_str}\t{qual_str}\n"
        );

        let dir = tempfile::tempdir()?;
        let sam_gz = create_indexed_sam_proptest(dir.path(), &sam_text);
        let mut reader = seqair::sam::reader::IndexedSamReader::open(&sam_gz)?;
        let tid = reader.header().tid("chr1").unwrap();
        let mut store = seqair::bam::RecordStore::new();
        reader.fetch_into(tid, Pos::<Zero>::new(1), Pos::<Zero>::new(100_000_000), &mut store)?;

        prop_assert_eq!(store.len(), 1);
        let seq = store.seq(0);
        prop_assert_eq!(seq.len(), bases.len());
        for (i, (&input_byte, &actual)) in bases.iter().zip(seq.iter()).enumerate() {
            let expected_base = match input_byte {
                b'A' => Base::A,
                b'C' => Base::C,
                b'G' => Base::G,
                b'T' => Base::T,
                // All IUPAC ambiguity codes and N must decode to Unknown.
                _ => Base::Unknown,
            };
            prop_assert_eq!(
                actual, expected_base,
                "pos {}: input '{}' expected {:?}",
                i, input_byte as char, expected_base
            );
        }
    }

    // r[verify sam.record.qual_decode]
    #[test]
    fn qual_decode_roundtrip(
        // Exclude quality 9: Phred+33 encodes it as ASCII '*' (42), which SAM
        // treats as "quality unavailable" when it is the entire QUAL field.
        quals in prop::collection::vec(
            (0u8..=93).prop_filter("not SAM unknown-qual sentinel", |&q| q + 33 != b'*'),
            1..200,
        ),
    ) {
        let seq_str: String = quals.iter().map(|_| 'A').collect();
        let qual_str: String = quals.iter().map(|&q| (q + 33) as char).collect();
        let cigar = format!("{}M", quals.len());
        let sam_text = format!(
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000000\n\
             read1\t99\tchr1\t100\t60\t{cigar}\t=\t200\t150\t{seq_str}\t{qual_str}\n"
        );

        let dir = tempfile::tempdir()?;
        let sam_gz = create_indexed_sam_proptest(dir.path(), &sam_text);
        let mut reader = seqair::sam::reader::IndexedSamReader::open(&sam_gz)?;
        let tid = reader.header().tid("chr1").unwrap();
        let mut store = seqair::bam::RecordStore::new();
        reader.fetch_into(tid, Pos::<Zero>::new(1), Pos::<Zero>::new(100_000_000), &mut store)?;

        prop_assert_eq!(store.len(), 1);
        let stored_qual = store.qual(0);
        prop_assert_eq!(stored_qual.len(), quals.len());
        for (i, (&expected, &actual)) in quals.iter().zip(stored_qual.iter()).enumerate() {
            prop_assert_eq!(actual, expected, "pos {}", i);
        }
    }
}

fn create_indexed_sam_proptest(dir: &std::path::Path, sam_text: &str) -> std::path::PathBuf {
    use std::process::Command;

    let sam_path = dir.join("test.sam");
    let sam_gz_path = dir.join("test.sam.gz");

    std::fs::write(&sam_path, sam_text).expect("write SAM");

    let status = Command::new("bgzip")
        .arg("-c")
        .arg(&sam_path)
        .stdout(std::fs::File::create(&sam_gz_path).expect("create .sam.gz"))
        .status()
        .expect("bgzip not found");
    assert!(status.success(), "bgzip failed");

    let status = Command::new("tabix")
        .args(["-p", "sam"])
        .arg(&sam_gz_path)
        .status()
        .expect("tabix not found");
    assert!(status.success(), "tabix failed");

    sam_gz_path
}
