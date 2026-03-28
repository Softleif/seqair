//! Compares seqair FASTA reader against noodles-fasta across the test data.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::indexing_slicing)]
use noodles::core::{Position, Region};
use noodles::fasta;
use proptest::prelude::*;
use seqair::fasta::IndexedFastaReader;
use std::path::Path;

fn test_fasta_path() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.fasta.gz"))
}

/// Sequences in the test FASTA with their lengths (from .fai).
const SEQUENCES: &[(&str, u64)] =
    &[("chr19", 61_431_566), ("2kb_3_Unmodified", 2018), ("bacteriophage_lambda_CpG", 48502)];

fn noodles_fetch(name: &str, start: u64, stop: u64) -> Vec<u8> {
    let mut reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(test_fasta_path())
        .expect("noodles fasta open");
    // seqair uses 0-based half-open [start, stop); noodles uses 1-based closed [start+1, stop]
    let region = Region::new(
        name,
        Position::try_from(start as usize + 1).unwrap()
            ..=Position::try_from(stop as usize).unwrap(),
    );
    let record = reader.query(&region).unwrap();
    record.sequence().as_ref().to_ascii_uppercase()
}

fn seqair_fetch(reader: &mut IndexedFastaReader, name: &str, start: u64, stop: u64) -> Vec<u8> {
    reader.fetch_seq(name, start, stop).expect("seqair fetch_seq")
}

// ---- Full sequence comparison for small sequences ----

// r[verify fasta.fetch.coordinates]
// r[verify fasta.fetch.newline_stripping]
#[test]
fn full_small_sequences_match_noodles() {
    let mut rio = IndexedFastaReader::open(test_fasta_path()).expect("rio open");

    for &(name, length) in SEQUENCES {
        if length > 100_000 {
            continue; // Skip large sequences for full comparison
        }

        let noodles = noodles_fetch(name, 0, length);
        let seq = seqair_fetch(&mut rio, name, 0, length);

        assert_eq!(
            noodles.len(),
            seq.len(),
            "length mismatch for {name}: noodles={}, seqair={}",
            noodles.len(),
            seq.len()
        );
        assert_eq!(noodles, seq, "sequence mismatch for {name} (full)");
    }
}

// ---- Windowed comparison across chr19 ----

// r[verify fasta.bgzf.decompress]
#[test]
fn chr19_windowed_comparison_noodles() {
    let mut rio = IndexedFastaReader::open(test_fasta_path()).expect("rio open");

    let windows: &[(u64, u64)] = &[
        (0, 100),
        (0, 1),
        (1000, 1100),
        (6_105_000, 6_106_000),
        (61_431_466, 61_431_566),
        (61_431_565, 61_431_566),
        (49, 52), // line boundary cross
        (100, 200),
        (30_000_000, 30_000_500),
    ];

    for &(start, stop) in windows {
        let noodles = noodles_fetch("chr19", start, stop);
        let seq = seqair_fetch(&mut rio, "chr19", start, stop);

        assert_eq!(
            noodles,
            seq,
            "mismatch for chr19:{start}-{stop} (len noodles={}, seqair={})",
            noodles.len(),
            seq.len()
        );
    }
}

// ---- Random region comparison via proptest ----

// r[verify fasta.bgzf.decompress]
// r[verify fasta.bgzf.sequential_read]
proptest! {
    #![proptest_config(ProptestConfig::with_cases(100))]

    #[test]
    fn random_regions_match_noodles(
        seq_idx in 0usize..3,
        start_frac in 0.0f64..0.99,
        len_frac in 0.001f64..0.01,
    ) {
        let (name, length) = SEQUENCES[seq_idx];
        let start = (start_frac * length as f64) as u64;
        let fetch_len = ((len_frac * length as f64) as u64).max(1).min(length - start);
        let stop = start + fetch_len;

        let mut rio = IndexedFastaReader::open(test_fasta_path()).expect("rio open");
        let seq = seqair_fetch(&mut rio, name, start, stop);
        let noodles = noodles_fetch(name, start, stop);

        prop_assert!(
            seq == noodles,
            "mismatch for {}:{}-{} (len seqair={}, noodles={})",
            name, start, stop, seq.len(), noodles.len()
        );
    }
}
