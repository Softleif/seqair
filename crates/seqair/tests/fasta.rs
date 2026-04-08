//! Property-based tests for the FASTA reader.
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
    reason = "test code with known small values"
)]
use proptest::prelude::*;
use rust_htslib::faidx;
use seqair::bam::{Pos, Zero};
use seqair::fasta::{FaiEntry, IndexedFastaReader};
use std::{io::Write, path::Path};
use tempfile::TempDir;

fn test_fasta_path() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.fasta.gz"))
}

// ---------------------------------------------------------------------------
// FAI byte offset properties
// ---------------------------------------------------------------------------

// r[verify fasta.index.offset_calculation]
proptest! {
    /// True end-to-end oracle: build an in-memory FASTA content string from
    /// generated parameters, then verify that `byte_offset(pos)` points to
    /// the correct character in the actual FASTA bytes.
    #[test]
    fn byte_offset_indexes_correct_character_in_fasta(
        linebases in 1u64..=80u64,
        extra in 1u64..=2u64,  // 1 = \n, 2 = \r\n
        bases in prop::collection::vec(
            prop::sample::select(vec![b'A', b'C', b'G', b'T']),
            1..=200usize,
        ),
    ) {
        let n_bases = bases.len() as u64;
        let linewidth = linebases + extra;

        // Build the FASTA sequence body: bases with newlines every `linebases`
        // characters. extra=1 means \n, extra=2 means \r\n.
        let mut content: Vec<u8> = Vec::new();
        for (i, &b) in bases.iter().enumerate() {
            content.push(b);
            let pos_in_line = (i as u64 + 1) % linebases;
            // After every `linebases` bases (but not after the very last base),
            // insert the line ending.
            if pos_in_line == 0 && (i as u64 + 1) < n_bases {
                if extra == 2 {
                    content.push(b'\r');
                }
                content.push(b'\n');
            }
        }

        // The FaiEntry has offset=0 (sequence starts at byte 0 of content).
        let entry = FaiEntry {
            name: "test".into(),
            length: n_bases,
            offset: 0,
            linebases,
            linewidth,
        };

        for pos in 0..n_bases {
            let byte_off = entry.byte_offset(pos);
            let actual_byte = content.get(byte_off as usize).copied();
            let expected_byte = bases.get(pos as usize).copied();
            prop_assert_eq!(
                actual_byte, expected_byte,
                "byte_offset({}) = {} indexes {:?} but expected {:?} \
                 (linebases={}, linewidth={}, content_len={})",
                pos, byte_off, actual_byte, expected_byte,
                linebases, linewidth, content.len()
            );
        }
    }
}

// ---------------------------------------------------------------------------
// GZI translate properties
// ---------------------------------------------------------------------------

fn make_gzi_data(entries: &[(u64, u64)]) -> Vec<u8> {
    let mut data = Vec::new();
    data.extend_from_slice(&(entries.len() as u64).to_le_bytes());
    for &(compressed, uncompressed) in entries {
        data.extend_from_slice(&compressed.to_le_bytes());
        data.extend_from_slice(&uncompressed.to_le_bytes());
    }
    data
}

// r[verify fasta.gzi.translate]
proptest! {
    /// Translating offsets within a block (≤ u16::MAX from block start) must
    /// produce the correct compressed_offset and within_block_offset.
    #[test]
    fn gzi_translate_same_block_valid(
        block_start_compressed in 0u64..1_000_000,
        block_start_uncompressed in 1u64..1_000_000,
        offset_a in 0u16..=u16::MAX,
        offset_b in 0u16..=u16::MAX,
    ) {
        let entries = vec![(block_start_compressed, block_start_uncompressed)];
        let data = make_gzi_data(&entries);
        let gzi = seqair::fasta::GziIndex::parse_test(&data);

        let target_a = block_start_uncompressed + u64::from(offset_a);
        let target_b = block_start_uncompressed + u64::from(offset_b);

        let loc_a = gzi.translate(target_a).unwrap();
        let loc_b = gzi.translate(target_b).unwrap();

        prop_assert_eq!(loc_a.compressed_offset, block_start_compressed);
        prop_assert_eq!(loc_b.compressed_offset, block_start_compressed);
        prop_assert_eq!(loc_a.within_block_offset, offset_a);
        prop_assert_eq!(loc_b.within_block_offset, offset_b);
    }

    /// For offsets before the first GZI entry, translation must map to
    /// compressed offset 0 with the uncompressed offset as within_block.
    #[test]
    fn gzi_translate_before_first_entry(
        block_start in 1000u64..1_000_000,
        target in 0u64..999,
    ) {
        let entries = vec![(block_start, 1000)];
        let data = make_gzi_data(&entries);
        let gzi = seqair::fasta::GziIndex::parse_test(&data);

        let loc = gzi.translate(target).unwrap();
        prop_assert_eq!(loc.compressed_offset, 0);
        prop_assert_eq!(loc.within_block_offset, target as u16);
    }

    /// Within-block offset exceeding u16::MAX must be rejected.
    #[test]
    fn gzi_translate_rejects_overflow(
        overflow in 1u64..100_000,
    ) {
        let data = make_gzi_data(&[]);
        let gzi = seqair::fasta::GziIndex::parse_test(&data);

        let target = u64::from(u16::MAX) + overflow;
        let result = gzi.translate(target);
        prop_assert!(result.is_err(),
            "translate({}) should fail but returned {:?}", target, result);
    }
}

// ---------------------------------------------------------------------------
// Plain FASTA roundtrip: generate, write, read back, verify
// ---------------------------------------------------------------------------

/// Strategy for generating a valid FASTA sequence with known content.
fn fasta_sequence_strategy() -> impl Strategy<Value = (Vec<u8>, u64)> {
    // linebases between 3 and 80, sequence length between 1 and 500
    (3u64..80, 1usize..500)
        .prop_flat_map(|(linebases, seq_len)| {
            let bases = prop::collection::vec(
                prop::sample::select(vec![b'A', b'C', b'G', b'T', b'N']),
                seq_len,
            );
            (bases, Just(linebases))
        })
        .prop_map(|(bases, linebases)| (bases, linebases))
}

// r[verify fasta.plain.read]
// r[verify fasta.fetch.newline_stripping]
// r[verify fasta.fetch.uppercase]
// r[verify fasta.fetch.coordinates]
proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    /// Write a random FASTA, read it back at a random position, verify content.
    #[test]
    fn plain_fasta_roundtrip(
        (bases, linebases) in fasta_sequence_strategy(),
        // We'll pick a random sub-range to fetch
        range_start_frac in 0.0f64..1.0,
        range_len_frac in 0.01f64..1.0,
    ) {
        let seq_len = bases.len() as u64;
        let start = (range_start_frac * seq_len as f64) as u64;
        let remaining = seq_len - start;
        let fetch_len = ((range_len_frac * remaining as f64) as u64).max(1).min(remaining);
        let stop = start + fetch_len;

        let dir = TempDir::new().unwrap();
        let fasta_path = dir.path().join("test.fa");
        let fai_path = dir.path().join("test.fa.fai");

        // Write FASTA with the given line length
        let linewidth = linebases + 1; // +1 for \n
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        writeln!(f, ">seq1").unwrap();
        for (i, &base) in bases.iter().enumerate() {
            f.write_all(&[base]).unwrap();
            if ((i + 1) as u64).is_multiple_of(linebases) && (i + 1) < bases.len() {
                f.write_all(b"\n").unwrap();
            }
        }
        f.write_all(b"\n").unwrap();

        // Compute FAI offset: ">seq1\n" = 6 bytes
        let offset = 6u64;
        let mut fai = std::fs::File::create(&fai_path).unwrap();
        writeln!(fai, "seq1\t{}\t{}\t{}\t{}", seq_len, offset, linebases, linewidth).unwrap();

        let mut reader = IndexedFastaReader::open(&fasta_path).unwrap();
        let start_pos = Pos::<Zero>::try_from_u64(start).unwrap();
        let stop_pos = Pos::<Zero>::try_from_u64(stop).unwrap();
        let fetched = reader.fetch_seq("seq1", start_pos, stop_pos).unwrap();

        let expected: Vec<u8> = bases[start as usize..stop as usize]
            .iter()
            .map(|b| b.to_ascii_uppercase())
            .collect();

        prop_assert_eq!(
            fetched, expected,
            "roundtrip failed: start={}, stop={}, linebases={}, seq_len={}",
            start, stop, linebases, seq_len
        );
    }
}

// ---------------------------------------------------------------------------
// Bgzip FASTA: random positions compared to htslib
// ---------------------------------------------------------------------------

const SEQUENCES: &[(&str, u64)] =
    &[("chr19", 61_431_566), ("2kb_3_Unmodified", 2018), ("bacteriophage_lambda_CpG", 48502)];

fn htslib_fetch(name: &str, start: u64, stop: u64) -> Vec<u8> {
    let reader = faidx::Reader::from_path(test_fasta_path()).expect("htslib faidx open");
    reader
        .fetch_seq(name, start as usize, (stop - 1) as usize)
        .expect("htslib fetch_seq")
        .to_ascii_uppercase()
}

// r[verify fasta.bgzf.decompress]
// r[verify fasta.bgzf.sequential_read]
// r[verify fasta.fork.equivalence]
proptest! {
    #![proptest_config(ProptestConfig::with_cases(100))]

    /// Fetch random regions from the bgzip test FASTA and compare to htslib.
    #[test]
    fn bgzf_random_regions_match_htslib(
        seq_idx in 0usize..3,
        start_frac in 0.0f64..0.99,
        len_frac in 0.001f64..0.01,
    ) {
        let (name, length) = SEQUENCES[seq_idx];
        let start = (start_frac * length as f64) as u64;
        let fetch_len = ((len_frac * length as f64) as u64).max(1).min(length - start);
        let stop = start + fetch_len;

        let mut rio = IndexedFastaReader::open(test_fasta_path()).expect("rio open");
        let start_pos = Pos::<Zero>::try_from_u64(start).unwrap();
        let stop_pos = Pos::<Zero>::try_from_u64(stop).unwrap();
        let seq = rio.fetch_seq(name, start_pos, stop_pos).expect("rio fetch");
        let hts_seq = htslib_fetch(name, start, stop);

        prop_assert!(
            seq == hts_seq,
            "mismatch for {}:{}-{} (len rio={}, hts={})",
            name, start, stop, seq.len(), hts_seq.len()
        );
    }

    /// Forked readers must produce identical results to the original at any position.
    #[test]
    fn fork_matches_original_at_random_positions(
        start_frac in 0.0f64..0.99,
        len_frac in 0.001f64..0.01,
    ) {
        let length = 61_431_566u64; // chr19
        let start = (start_frac * length as f64) as u64;
        let fetch_len = ((len_frac * length as f64) as u64).max(1).min(length - start);
        let stop = start + fetch_len;

        let mut rio = IndexedFastaReader::open(test_fasta_path()).expect("rio open");
        let mut forked = rio.fork().expect("fork");

        let start_pos = Pos::<Zero>::try_from_u64(start).unwrap();
        let stop_pos = Pos::<Zero>::try_from_u64(stop).unwrap();
        let orig = rio.fetch_seq("chr19", start_pos, stop_pos).expect("orig fetch");
        let fork_result = forked.fetch_seq("chr19", start_pos, stop_pos).expect("fork fetch");

        prop_assert_eq!(orig, fork_result,
            "fork mismatch at chr19:{}-{}", start, stop);
    }
}

// ---------------------------------------------------------------------------
// fetch_base_seq: buffer reuse across calls
// ---------------------------------------------------------------------------

// r[verify fasta.fetch.buffer_reuse]
#[test]
fn fetch_base_seq_reuses_buffer() {
    use seqair_types::Base;

    let mut readers = seqair::Readers::open(
        Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.bam")),
        test_fasta_path(),
    )
    .unwrap();

    // First call: buffer starts empty, gets allocated
    let seq1 = readers
        .fetch_base_seq(
            "bacteriophage_lambda_CpG",
            Pos::<Zero>::new(0).unwrap(),
            Pos::<Zero>::new(100).unwrap(),
        )
        .unwrap();
    assert_eq!(seq1.len(), 100);
    // Every element must be a valid Base
    for &b in seq1.iter() {
        assert!(matches!(b, Base::A | Base::C | Base::G | Base::T | Base::Unknown));
    }

    // Second call: should reuse the internal buffer (no way to observe capacity
    // directly, but we can verify correctness across calls)
    let seq2 = readers
        .fetch_base_seq(
            "bacteriophage_lambda_CpG",
            Pos::<Zero>::new(100).unwrap(),
            Pos::<Zero>::new(300).unwrap(),
        )
        .unwrap();
    assert_eq!(seq2.len(), 200);

    // Third call: same region as first, must produce identical result
    let seq3 = readers
        .fetch_base_seq(
            "bacteriophage_lambda_CpG",
            Pos::<Zero>::new(0).unwrap(),
            Pos::<Zero>::new(100).unwrap(),
        )
        .unwrap();
    assert_eq!(seq1, seq3, "repeated fetch must produce identical results");

    // Verify against raw fetch + manual conversion
    let mut raw_reader = IndexedFastaReader::open(test_fasta_path()).unwrap();
    let raw = raw_reader
        .fetch_seq(
            "bacteriophage_lambda_CpG",
            Pos::<Zero>::new(0).unwrap(),
            Pos::<Zero>::new(100).unwrap(),
        )
        .unwrap();
    let expected: Vec<Base> = Base::from_ascii_vec(raw);
    assert_eq!(&*seq1, &expected[..], "fetch_base_seq must match from_ascii_vec on raw fetch");
}

// ---------------------------------------------------------------------------
// Fetch into buffer: must produce same result as allocating fetch
// ---------------------------------------------------------------------------

// r[verify fasta.fetch.buffer_reuse]
proptest! {
    #![proptest_config(ProptestConfig::with_cases(50))]

    #[test]
    fn fetch_into_matches_fetch_at_random_positions(
        start_frac in 0.0f64..0.99,
        len_frac in 0.001f64..0.01,
    ) {
        let length = 48502u64; // bacteriophage_lambda_CpG
        let start = (start_frac * length as f64) as u64;
        let fetch_len = ((len_frac * length as f64) as u64).max(1).min(length - start);
        let stop = start + fetch_len;

        let mut rio = IndexedFastaReader::open(test_fasta_path()).expect("rio open");
        let mut buf = Vec::new();

        let start_pos = Pos::<Zero>::try_from_u64(start).unwrap();
        let stop_pos = Pos::<Zero>::try_from_u64(stop).unwrap();
        let alloc = rio.fetch_seq("bacteriophage_lambda_CpG", start_pos, stop_pos).expect("fetch_seq");
        rio.fetch_seq_into("bacteriophage_lambda_CpG", start_pos, stop_pos, &mut buf)
            .expect("fetch_seq_into");

        prop_assert_eq!(alloc.clone(), buf, "fetch vs fetch_into mismatch at {}-{}", start, stop);
        let hts = htslib_fetch("bacteriophage_lambda_CpG", start, stop);
        prop_assert_eq!(
            alloc, hts,
            "seqair vs htslib mismatch at {}-{}", start, stop
        );
    }
}
