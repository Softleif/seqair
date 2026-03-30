//! Compares seqair FASTA reader against htslib faidx across the test data.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::indexing_slicing)]
use rust_htslib::faidx;
use seqair::fasta::{FastaError, IndexedFastaReader};
use seqair_types::{Pos, Zero};
use std::path::Path;
use tempfile::TempDir;

fn test_fasta_path() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.fasta.gz"))
}

/// Sequences in the test FASTA with their lengths (from .fai).
const SEQUENCES: &[(&str, u64)] =
    &[("chr19", 61_431_566), ("2kb_3_Unmodified", 2018), ("bacteriophage_lambda_CpG", 48502)];

fn htslib_fetch(name: &str, start: u64, stop: u64) -> Vec<u8> {
    let reader = faidx::Reader::from_path(test_fasta_path()).expect("htslib faidx open");
    let begin = start as usize;
    let end = stop.saturating_sub(1) as usize;
    reader.fetch_seq(name, begin, end).expect("htslib fetch_seq").to_ascii_uppercase()
}

fn fetch(reader: &mut IndexedFastaReader, name: &str, start: u64, stop: u64) -> Vec<u8> {
    let start = Pos::<Zero>::try_from_u64(start).expect("start fits in u32");
    let stop = Pos::<Zero>::try_from_u64(stop).expect("stop fits in u32");
    reader.fetch_seq(name, start, stop).expect("rio fetch_seq")
}

// ---- Full sequence comparison for small sequences ----

// r[verify fasta.bgzf.detect]
// r[verify fasta.bgzf.gzi_required]
#[test]
fn full_small_sequences_match_htslib() {
    let mut rio = IndexedFastaReader::open(test_fasta_path()).expect("rio open");

    for &(name, length) in SEQUENCES {
        if length > 100_000 {
            continue; // Skip large sequences for full comparison
        }

        let hts = htslib_fetch(name, 0, length);
        let seq = fetch(&mut rio, name, 0, length);

        assert_eq!(
            hts.len(),
            seq.len(),
            "length mismatch for {name}: htslib={}, rio={}",
            hts.len(),
            seq.len()
        );
        assert_eq!(hts, seq, "sequence mismatch for {name} (full)");
    }
}

// ---- Windowed comparison across chr19 ----

#[test]
fn chr19_windowed_comparison() {
    let mut rio = IndexedFastaReader::open(test_fasta_path()).expect("rio open");

    // Test various windows across chr19
    let windows: &[(u64, u64)] = &[
        (0, 100),                 // start of chromosome
        (0, 1),                   // single base at start
        (1000, 1100),             // arbitrary middle region
        (6_105_000, 6_106_000),   // region with BAM coverage
        (6_105_700, 6_105_800),   // typical segment fetch
        (61_431_466, 61_431_566), // last 100 bases
        (61_431_565, 61_431_566), // very last base
        (0, 50),                  // first partial line
        (49, 52),                 // cross first line boundary (50 bases/line)
        (100, 200),               // cross multiple lines
        (30_000_000, 30_000_500), // middle of chromosome
        (6_103_076, 6_143_229),   // full BAM coverage region
    ];

    for &(start, stop) in windows {
        let hts = htslib_fetch("chr19", start, stop);
        let seq = fetch(&mut rio, "chr19", start, stop);

        assert_eq!(
            hts,
            seq,
            "mismatch for chr19:{start}-{stop} (len hts={}, rio={})",
            hts.len(),
            seq.len()
        );
    }
}

// ---- Large fetch comparison ----

#[test]
fn large_fetch_matches_htslib() {
    let mut rio = IndexedFastaReader::open(test_fasta_path()).expect("rio open");

    // 100KB fetch — typical segment size
    let start = 6_100_000;
    let stop = 6_200_000;
    let hts = htslib_fetch("chr19", start, stop);
    let seq = fetch(&mut rio, "chr19", start, stop);

    assert_eq!(hts.len(), seq.len(), "length mismatch for 100KB fetch");
    assert_eq!(hts, seq, "100KB fetch mismatch");
}

// ---- Fork produces identical results ----

// r[verify fasta.fork.operation]
// r[verify fasta.fork.independence]
#[test]
fn forked_reader_matches_original() {
    let mut rio = IndexedFastaReader::open(test_fasta_path()).expect("rio open");
    let mut forked = rio.fork().expect("fork");

    for &(name, length) in SEQUENCES {
        let stop = length.min(1000);
        let orig = fetch(&mut rio, name, 0, stop);
        let fork_seq = fetch(&mut forked, name, 0, stop);
        assert_eq!(orig, fork_seq, "fork mismatch for {name}");
    }
}

// ---- Fork shares Arc ----

// r[verify fasta.fork.shared_state]
#[test]
fn fork_shares_arc() {
    let rio = IndexedFastaReader::open(test_fasta_path()).expect("rio open");
    let forked = rio.fork().expect("fork");
    assert!(rio.shares_index_with(&forked));
}

// ---- Buffer reuse ----

#[test]
fn fetch_into_matches_fetch() {
    let mut rio = IndexedFastaReader::open(test_fasta_path()).expect("rio open");
    let mut buf = Vec::new();

    for &(name, length) in SEQUENCES {
        let stop = Pos::<Zero>::try_from_u64(length.min(500)).expect("stop fits in u32");
        let zero = Pos::<Zero>::new(0).unwrap();
        let alloc = rio.fetch_seq(name, zero, stop).expect("fetch_seq");
        rio.fetch_seq_into(name, zero, stop, &mut buf).expect("fetch_seq_into");
        assert_eq!(alloc, buf, "fetch_seq vs fetch_seq_into mismatch for {name}");
    }
}

// ---- Error cases ----

// r[verify fasta.fetch.unknown_sequence]
// r[verify fasta.errors]
#[test]
fn unknown_sequence_error() {
    let mut rio = IndexedFastaReader::open(test_fasta_path()).expect("rio open");
    let err = rio
        .fetch_seq("nonexistent", Pos::<Zero>::new(0).unwrap(), Pos::<Zero>::new(100).unwrap())
        .unwrap_err();
    let msg = format!("{err}");
    assert!(msg.contains("nonexistent"), "error should name the sequence: {msg}");
    assert!(msg.contains("chr19"), "error should list available: {msg}");
}

// r[verify fasta.fetch.bounds_check]
// r[verify fasta.errors]
#[test]
fn out_of_bounds_error() {
    let mut rio = IndexedFastaReader::open(test_fasta_path()).expect("rio open");
    let err = rio
        .fetch_seq(
            "2kb_3_Unmodified",
            Pos::<Zero>::new(0).unwrap(),
            Pos::<Zero>::new(999999).unwrap(),
        )
        .unwrap_err();
    let msg = format!("{err}");
    assert!(msg.contains("out of bounds"), "error: {msg}");
}

#[test]
fn empty_range_error() {
    let mut rio = IndexedFastaReader::open(test_fasta_path()).expect("rio open");
    let err = rio
        .fetch_seq("chr19", Pos::<Zero>::new(100).unwrap(), Pos::<Zero>::new(100).unwrap())
        .unwrap_err();
    let msg = format!("{err}");
    assert!(msg.contains("out of bounds"), "error: {msg}");
}

// r[verify fasta.bgzf.gzi_missing]
// r[verify fasta.bgzf.gzi_required]
#[test]
fn missing_gzi_error() {
    // Copy the bgzf FASTA and its .fai into a tempdir, but omit the .gzi.
    // Opening should produce GziIndexNotFound with the expected path in the error.
    let dir = TempDir::new().unwrap();
    let src = test_fasta_path();
    let dst = dir.path().join("test.fasta.gz");
    std::fs::copy(src, &dst).expect("copy fasta");
    let fai_src = src.parent().unwrap().join("test.fasta.gz.fai");
    let fai_dst = dir.path().join("test.fasta.gz.fai");
    std::fs::copy(&fai_src, &fai_dst).expect("copy fai");
    // Intentionally do NOT copy .gzi

    let err = IndexedFastaReader::open(&dst).unwrap_err();
    match &err {
        FastaError::GziIndexNotFound { expected_path, .. } => {
            let path_str = expected_path.to_str().unwrap();
            assert!(path_str.ends_with(".gzi"), "expected .gzi path, got: {path_str}");
        }
        _ => panic!("expected GziIndexNotFound, got {err:?}"),
    }
}

// ---- Concurrent fork usage (simulates multi-threaded) ----

#[test]
fn concurrent_forks() {
    let rio = IndexedFastaReader::open(test_fasta_path()).expect("rio open");

    std::thread::scope(|s| {
        let handles: Vec<_> = (0..4)
            .map(|i| {
                let ref_ = &rio;
                s.spawn(move || {
                    let mut forked = ref_.fork().expect("fork");
                    let start = (i as u64) * 10000;
                    let stop = start + 10000;
                    let start_pos = Pos::<Zero>::try_from_u64(start).expect("start fits");
                    let stop_pos = Pos::<Zero>::try_from_u64(stop).expect("stop fits");
                    forked.fetch_seq("chr19", start_pos, stop_pos).expect("fetch")
                })
            })
            .collect();

        let results: Vec<_> = handles.into_iter().map(|h| h.join().unwrap()).collect();

        // Verify each result matches htslib
        for (i, seq) in results.iter().enumerate() {
            let start = (i as u64) * 10000;
            let stop = start + 10000;
            let hts = htslib_fetch("chr19", start, stop);
            assert_eq!(&hts, seq, "concurrent fork mismatch for window {i}");
        }
    });
}
