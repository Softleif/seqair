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
use core::range::RangeInclusive;
use hegel::prelude::*;
use rust_htslib::faidx;
use seqair::bam::Pos0;
use seqair::fasta::{FaiEntry, IndexedFastaReader};
use std::{io::Write, path::Path};
use tempfile::TempDir;

fn test_fasta_path() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.fasta.gz"))
}

/// `start` inclusive, `stop` exclusive, the half-open convention every fetch
/// in this file is expressed in — converted to seqair's closed
/// `[start, stop - 1]` span here.
fn span(start: u64, stop: u64) -> RangeInclusive<Pos0> {
    RangeInclusive {
        start: Pos0::try_from(start).unwrap(),
        last: Pos0::try_from(stop - 1).unwrap(),
    }
}

// ---------------------------------------------------------------------------
// FAI byte offset properties
// ---------------------------------------------------------------------------

// r[verify fasta.index.offset_calculation]
/// True end-to-end oracle: build an in-memory FASTA content string from
/// generated parameters, then verify that `byte_offset(pos)` points to
/// the correct character in the actual FASTA bytes.
#[hegel::test]
fn byte_offset_indexes_correct_character_in_fasta(tc: TestCase) {
    let linebases = tc.draw(gs::integers::<u64>().min_value(1).max_value(80u64));
    let extra = tc.draw(gs::integers::<u64>().min_value(1).max_value(2u64));
    let bases = tc.draw(gs::vecs(gs::sampled_from(ACGT)).min_size(1).max_size(200));
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
        qual_offset: None,
    };

    for pos in 0..n_bases {
        let byte_off = entry.byte_offset(pos).expect("small test offsets never overflow u64");
        let actual_byte = content.get(byte_off as usize).copied();
        let expected_byte = bases.get(pos as usize).copied();
        assert_eq!(
            actual_byte,
            expected_byte,
            "byte_offset({}) = {} indexes {:?} but expected {:?} \
             (linebases={}, linewidth={}, content_len={})",
            pos,
            byte_off,
            actual_byte,
            expected_byte,
            linebases,
            linewidth,
            content.len()
        );
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
/// Translating offsets within a block (≤ `u16::MAX` from block start) must
/// produce the correct `compressed_offset` and `within_block_offset`.
#[hegel::test]
fn gzi_translate_same_block_valid(tc: TestCase) {
    let block_start_compressed = tc.draw(gs::integers::<u64>().max_value(999999));
    let block_start_uncompressed = tc.draw(gs::integers::<u64>().min_value(1).max_value(999999));
    let offset_a = tc.draw(gs::integers::<u16>());
    let offset_b = tc.draw(gs::integers::<u16>());
    let entries = vec![(block_start_compressed, block_start_uncompressed)];
    let data = make_gzi_data(&entries);
    let gzi = seqair::fasta::GziIndex::parse_test(&data);

    let target_a = block_start_uncompressed + u64::from(offset_a);
    let target_b = block_start_uncompressed + u64::from(offset_b);

    let loc_a = gzi.translate(target_a).unwrap();
    let loc_b = gzi.translate(target_b).unwrap();

    assert_eq!(loc_a.compressed_offset, block_start_compressed);
    assert_eq!(loc_b.compressed_offset, block_start_compressed);
    assert_eq!(loc_a.within_block_offset, offset_a);
    assert_eq!(loc_b.within_block_offset, offset_b);
}

/// For offsets before the first GZI entry, translation must map to
/// compressed offset 0 with the uncompressed offset as `within_block`.
#[hegel::test]
fn gzi_translate_before_first_entry(tc: TestCase) {
    let block_start = tc.draw(gs::integers::<u64>().min_value(1000).max_value(999999));
    let target = tc.draw(gs::integers::<u64>().max_value(998));
    let entries = vec![(block_start, 1000)];
    let data = make_gzi_data(&entries);
    let gzi = seqair::fasta::GziIndex::parse_test(&data);

    let loc = gzi.translate(target).unwrap();
    assert_eq!(loc.compressed_offset, 0);
    assert_eq!(loc.within_block_offset, target as u16);
}

/// Within-block offset exceeding `u16::MAX` must be rejected.
#[hegel::test]
fn gzi_translate_rejects_overflow(tc: TestCase) {
    let overflow = tc.draw(gs::integers::<u64>().min_value(1).max_value(99999));
    let data = make_gzi_data(&[]);
    let gzi = seqair::fasta::GziIndex::parse_test(&data);

    let target = u64::from(u16::MAX) + overflow;
    let err = gzi
        .translate(target)
        .expect_err("an offset past u16::MAX inside one block cannot be addressed");
    assert!(
        matches!(err, seqair::fasta::GziError::WithinBlockOverflow { .. }),
        "translate({target}) failed with {err:?}, not WithinBlockOverflow"
    );
}

// ---------------------------------------------------------------------------
// Plain FASTA roundtrip: generate, write, read back, verify
// ---------------------------------------------------------------------------

/// The bases a generated FASTA sequence is drawn from.
const ACGT: &[u8] = b"ACGT";
const ACGTN: &[u8] = b"ACGTN";

/// A valid FASTA sequence with known content, as `(bases, linebases)`.
#[hegel::composite]
fn fasta_sequence(tc: &TestCase) -> (Vec<u8>, u64) {
    // linebases between 3 and 79, sequence length between 1 and 499
    let linebases = tc.draw_silent(gs::integers::<u64>().min_value(3).max_value(79));
    let seq_len = tc.draw_silent(gs::integers::<usize>().min_value(1).max_value(499));
    let bases =
        tc.draw_silent(gs::vecs(gs::sampled_from(ACGTN)).min_size(seq_len).max_size(seq_len));
    (bases, linebases)
}

// r[verify fasta.plain.read]
// r[verify fasta.fetch.newline_stripping]
// r[verify fasta.fetch.uppercase]
// r[verify fasta.fetch.coordinates]

/// Write a random FASTA, read it back at a random position, verify content.
#[hegel::test(test_cases = 64)]
fn plain_fasta_roundtrip(tc: TestCase) {
    let (bases, linebases) = tc.draw(fasta_sequence());
    let range_start_frac = tc.draw(gs::floats::<f64>().min_value(0.0).max_value_exclusive(1.0));
    let range_len_frac = tc.draw(gs::floats::<f64>().min_value(0.01).max_value_exclusive(1.0));
    let seq_len = bases.len() as u64;
    let start = (range_start_frac * seq_len as f64) as u64;
    let remaining = seq_len - start;
    let fetch_len = ((range_len_frac * remaining as f64) as u64).max(1).min(remaining);
    let stop = start + fetch_len;

    let dir = TempDir::new().unwrap();
    let fasta_path = write_plain_fasta(dir.path(), &bases, linebases);

    let mut reader = IndexedFastaReader::open(&fasta_path).unwrap();
    let fetched = reader.fetch_seq("seq1", span(start, stop)).unwrap();

    let expected: Vec<u8> =
        bases[start as usize..stop as usize].iter().map(|b| b.to_ascii_uppercase()).collect();

    assert_eq!(
        fetched, expected,
        "roundtrip failed: start={}, stop={}, linebases={}, seq_len={}",
        start, stop, linebases, seq_len
    );
}

/// Write `bases` as a one-sequence FASTA named `seq1`, wrapped every
/// `linebases`, with its `.fai` beside it. Returns the FASTA path.
fn write_plain_fasta(dir: &Path, bases: &[u8], linebases: u64) -> std::path::PathBuf {
    let fasta_path = dir.join("test.fa");
    let fai_path = dir.join("test.fa.fai");

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

    // FAI offset: ">seq1\n" = 6 bytes
    let offset = 6u64;
    let mut fai = std::fs::File::create(&fai_path).unwrap();
    writeln!(fai, "seq1\t{}\t{}\t{}\t{}", bases.len(), offset, linebases, linewidth).unwrap();
    fasta_path
}

// r[verify fasta.fetch.coordinates]
// r[verify fasta.fetch.bounds_check]
// r[verify interval.span_type]
/// The closed span is the whole contract: a request is served iff
/// `start <= last < len`, and then it is exactly `bases[start..=last]`.
/// Spans are drawn around every edge — the first base, the last base, one
/// past it, and reversed — because that is where a `+ 1` would hide.
#[hegel::test(test_cases = 64)]
fn fetch_seq_serves_exactly_the_closed_span_inside_the_sequence(tc: TestCase) {
    let (bases, linebases) = tc.draw(fasta_sequence());
    let len = bases.len() as u64;
    // Either end lands anywhere in `0..=len + 2`, biased to the edges.
    let edge = |tc: &TestCase| -> u64 {
        match tc.draw(gs::integers::<u8>().max_value(2)) {
            0 => tc.draw(gs::integers::<u64>().max_value(2)),
            1 => tc.draw(gs::integers::<u64>().min_value(len.saturating_sub(2)).max_value(len + 2)),
            _ => tc.draw(gs::integers::<u64>().max_value(len + 2)),
        }
    };
    let (start, last) = (edge(&tc), edge(&tc));
    let request = RangeInclusive {
        start: Pos0::try_from(start).unwrap(),
        last: Pos0::try_from(last).unwrap(),
    };

    let dir = TempDir::new().unwrap();
    let mut reader =
        IndexedFastaReader::open(&write_plain_fasta(dir.path(), &bases, linebases)).unwrap();

    match reader.fetch_seq("seq1", request) {
        Ok(fetched) => {
            assert!(start <= last && last < len, "{start}..={last} served from {len} bases");
            let expected: Vec<u8> =
                bases[start as usize..=last as usize].iter().map(u8::to_ascii_uppercase).collect();
            assert_eq!(fetched, expected, "{start}..={last} of {len}");
            tc.event("served");
        }
        Err(seqair::fasta::FastaError::RegionOutOfBounds {
            start: reported_start,
            last: reported_last,
            seq_len,
            ..
        }) => {
            assert!(start > last || last >= len, "{start}..={last} refused from {len} bases");
            assert_eq!((reported_start, reported_last, seq_len), (start, last, len));
            tc.event(if start > last { "reversed" } else { "past the end" });
        }
        Err(other) => panic!("{start}..={last} of {len}: unexpected error {other}"),
    }
}

// r[verify fasta.index.terminator_bound]
/// A near-`u64::MAX` offset makes the byte-offset calculation overflow `u64`
/// for the end of the requested span. That is a corrupt index, so it must be
/// a typed error — not the panic the arithmetic would otherwise produce, and
/// not a wrapped value silently treated as valid.
#[test]
fn wrapped_byte_span_is_an_error_not_a_panic() {
    let dir = TempDir::new().unwrap();
    let fasta_path = dir.path().join("wrap.fa");
    let mut f = std::fs::File::create(&fasta_path).unwrap();
    writeln!(f, ">seq1").unwrap();
    writeln!(f, "ACGTACGTAC").unwrap();
    drop(f);

    // Line geometry is legal; only `offset` is hostile.
    let mut fai = std::fs::File::create(dir.path().join("wrap.fa.fai")).unwrap();
    writeln!(fai, "seq1\t1000\t{}\t50\t51", u64::MAX - 10).unwrap();
    drop(fai);

    let mut reader = IndexedFastaReader::open(&fasta_path).unwrap();
    let err = reader
        .fetch_seq("seq1", span(0, 100))
        .expect_err("an overflowing offset must not be served");
    assert!(
        matches!(err, seqair::fasta::FastaError::IndexOffsetOverflow { .. }),
        "expected IndexOffsetOverflow, got {err:?}"
    );
}

// r[verify fasta.index.terminator_bound]
/// Before `byte_offset` used checked arithmetic, a `u64::MAX`-adjacent offset
/// that overflowed for *both* endpoints of a span wrapped each endpoint down
/// by the same amount, so the span length came out correct even though both
/// endpoints landed at the wrong (small, in-bounds-looking) file location.
/// The reader would have served 11 bytes from near the start of the file
/// instead of erroring. Checked arithmetic must catch this before any read.
#[test]
fn wraparound_offset_is_rejected_not_served_from_wrong_location() {
    let dir = TempDir::new().unwrap();
    let fasta_path = dir.path().join("wrap2.fa");
    let mut f = std::fs::File::create(&fasta_path).unwrap();
    writeln!(f, ">seq1").unwrap();
    writeln!(f, "ACGTACGTACGTACGTACGTACGTACGTAC").unwrap();
    drop(f);

    // A single-line sequence (linebases > every requested position), so
    // `byte_offset(pos) == offset + pos` exactly. `offset` sits 5 below
    // `u64::MAX`, so `offset + 10` and `offset + 20` both overflow, wrapping
    // to 4 and 14 under naive wrapping arithmetic — a plausible, ordered,
    // in-bounds span that is nonetheless nowhere near the real data.
    let mut fai = std::fs::File::create(dir.path().join("wrap2.fa.fai")).unwrap();
    writeln!(fai, "seq1\t1000\t{}\t1000\t1000", u64::MAX - 5).unwrap();
    drop(fai);

    let mut reader = IndexedFastaReader::open(&fasta_path).unwrap();
    let err = reader
        .fetch_seq("seq1", span(10, 21))
        .expect_err("an overflowing offset must not be served");
    assert!(
        matches!(err, seqair::fasta::FastaError::IndexOffsetOverflow { .. }),
        "expected IndexOffsetOverflow, got {err:?}"
    );
}

// r[verify fasta.index.data_bound]
/// A `.fai` can claim any `length` for a tiny plain file. The span computed
/// from that lie must be rejected against the real file size before the read
/// buffer is sized, not discovered as a short read (or an OOM) afterward.
#[test]
fn plain_fasta_span_beyond_file_is_rejected_before_allocation() {
    let dir = TempDir::new().unwrap();
    let fasta_path = dir.path().join("short.fa");
    let mut f = std::fs::File::create(&fasta_path).unwrap();
    writeln!(f, ">seq1").unwrap();
    writeln!(f, "ACGT").unwrap();
    drop(f); // the file holds ~11 bytes total

    // The FAI claims a sequence far longer than the file backing it, with
    // line geometry that puts the requested range's byte offset well past
    // the real end of the file.
    let mut fai = std::fs::File::create(dir.path().join("short.fa.fai")).unwrap();
    writeln!(fai, "seq1\t2000000000\t6\t60\t61").unwrap();
    drop(fai);

    let mut reader = IndexedFastaReader::open(&fasta_path).unwrap();
    let err = reader
        .fetch_seq("seq1", span(1_999_999_900, 1_999_999_950))
        .expect_err("a span past the real file size must not be served");
    assert!(
        matches!(err, seqair::fasta::FastaError::SpanBeyondData { .. }),
        "expected SpanBeyondData, got {err:?}"
    );
}

// r[verify fasta.index.data_bound]
/// The BGZF/GZI counterpart of [`plain_fasta_span_beyond_file_is_rejected_before_allocation`]:
/// the data bound comes from the GZI (`last entry + one block`), not the
/// compressed file's own size, and a `.fai` claiming a far longer sequence
/// than that bound must still be rejected before the read buffer is sized.
#[test]
fn bgzf_fasta_span_beyond_gzi_bound_is_rejected_before_allocation() {
    let dir = TempDir::new().unwrap();
    let fasta_path = dir.path().join("short.fa.gz");

    {
        let file = std::fs::File::create(&fasta_path).unwrap();
        let mut writer = seqair::io::BgzfWriter::new(file);
        writer.write_all(b">seq1\nACGT\n").unwrap();
        writer.finish().unwrap();
    }

    // An empty GZI (no entries) means "one block, at most 64 KiB
    // uncompressed" — `FaiEntry::byte_offset`'s data-bound counterpart to the
    // plain-file test above.
    let gzi_data = make_gzi_data(&[]);
    std::fs::write(dir.path().join("short.fa.gz.gzi"), &gzi_data).unwrap();

    // The FAI claims a sequence long enough that this range's byte offset
    // (~101 KB) sits past the ~64 KiB GZI bound.
    let mut fai = std::fs::File::create(dir.path().join("short.fa.gz.fai")).unwrap();
    writeln!(fai, "seq1\t1000000\t6\t60\t61").unwrap();
    drop(fai);

    let mut reader = IndexedFastaReader::open(&fasta_path).unwrap();
    let err = reader
        .fetch_seq("seq1", span(100_000, 100_050))
        .expect_err("a span past the GZI-derived bound must not be served");
    assert!(
        matches!(err, seqair::fasta::FastaError::SpanBeyondData { .. }),
        "expected SpanBeyondData, got {err:?}"
    );
}

// r[verify fasta.index.terminator_bound]
/// The bound that keeps a fetch's allocation proportional to the bases asked
/// for. Without it this index would size a read buffer at ~97.7 PiB and abort
/// the process before touching the file.
#[test]
fn absurd_line_width_is_rejected_at_open() {
    let dir = TempDir::new().unwrap();
    let fasta_path = dir.path().join("evil.fa");
    let mut f = std::fs::File::create(&fasta_path).unwrap();
    writeln!(f, ">evil").unwrap();
    writeln!(f, "ACGT").unwrap();
    drop(f);

    let mut fai = std::fs::File::create(dir.path().join("evil.fa.fai")).unwrap();
    writeln!(fai, "evil\t1000000\t0\t1\t1099511627776").unwrap();
    drop(fai);

    let err = IndexedFastaReader::open(&fasta_path).expect_err("hostile index must be rejected");
    assert!(
        matches!(
            err,
            seqair::fasta::FastaError::Fai {
                source: seqair::fasta::FaiError::InvalidEntry {
                    kind: seqair::fasta::FaiEntryError::LineTerminatorTooLong { .. },
                    ..
                }
            }
        ),
        "expected LineTerminatorTooLong, got {err:?}"
    );
}

// r[verify fasta.plain.positional_read]

/// The plain-FASTA fetch is a positional read, so one reader must answer an
/// arbitrarily ordered stream of overlapping requests without any of them
/// perturbing the next. Each is checked against the generated bases, and the
/// last request is repeated at the end so a drifting file position would
/// show up as a disagreement with its own earlier answer.
#[hegel::test(test_cases = 48)]
fn out_of_order_fetches_on_one_reader_are_independent(tc: TestCase) {
    let (bases, linebases) = tc.draw(fasta_sequence());
    let frac = || gs::floats::<f64>().min_value(0.0).max_value_exclusive(1.0);
    let raw_ranges = tc.draw(gs::vecs(gs::tuples!(frac(), frac())).min_size(2).max_size(11));
    let seq_len = bases.len() as u64;
    let dir = TempDir::new().unwrap();
    let fasta_path = dir.path().join("plain.fa");
    let linewidth = linebases + 1;

    let mut f = std::fs::File::create(&fasta_path).unwrap();
    writeln!(f, ">seq1").unwrap();
    for chunk in bases.chunks(linebases as usize) {
        f.write_all(chunk).unwrap();
        f.write_all(b"\n").unwrap();
    }
    drop(f);
    let mut fai = std::fs::File::create(dir.path().join("plain.fa.fai")).unwrap();
    writeln!(fai, "seq1\t{}\t6\t{}\t{}", seq_len, linebases, linewidth).unwrap();
    drop(fai);

    // Map the generated fractions onto in-bounds, non-empty ranges.
    let ranges: Vec<(u64, u64)> = raw_ranges
        .iter()
        .map(|(a, b)| {
            let x = ((a * seq_len as f64) as u64).min(seq_len - 1);
            let y = ((b * seq_len as f64) as u64).min(seq_len - 1);
            let (lo, hi) = if x <= y { (x, y) } else { (y, x) };
            (lo, hi + 1)
        })
        .collect();

    let mut reader = IndexedFastaReader::open(&fasta_path).unwrap();
    let fetch = |reader: &mut IndexedFastaReader, start: u64, stop: u64| {
        reader.fetch_seq("seq1", span(start, stop)).unwrap()
    };

    let mut first_answer = None;
    for (i, &(start, stop)) in ranges.iter().enumerate() {
        let got = fetch(&mut reader, start, stop);
        let mut expected = bases[start as usize..stop as usize].to_vec();
        expected.make_ascii_uppercase();
        assert_eq!(&got, &expected, "request {} [{}, {})", i, start, stop);
        if i == 0 {
            first_answer = Some(got);
        }
    }

    // Re-issue the first request last: every intervening read must have
    // left the reader able to answer it identically.
    let (start, stop) = ranges[0];
    assert_eq!(
        fetch(&mut reader, start, stop),
        first_answer.unwrap(),
        "repeating the first request after {} others changed its answer",
        ranges.len() - 1
    );
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

/// Fetch random regions from the bgzip test FASTA and compare to htslib.
#[hegel::test(test_cases = 100)]
fn bgzf_random_regions_match_htslib(tc: TestCase) {
    let seq_idx = tc.draw(gs::integers::<usize>().max_value(2));
    let start_frac = tc.draw(gs::floats::<f64>().min_value(0.0).max_value_exclusive(0.99));
    let len_frac = tc.draw(gs::floats::<f64>().min_value(0.001).max_value_exclusive(0.01));
    let (name, length) = SEQUENCES[seq_idx];
    let start = (start_frac * length as f64) as u64;
    let fetch_len = ((len_frac * length as f64) as u64).max(1).min(length - start);
    let stop = start + fetch_len;

    let mut rio = IndexedFastaReader::open(test_fasta_path()).expect("rio open");
    let seq = rio.fetch_seq(name, span(start, stop)).expect("rio fetch");
    let hts_seq = htslib_fetch(name, start, stop);

    assert!(
        seq == hts_seq,
        "mismatch for {}:{}-{} (len rio={}, hts={})",
        name,
        start,
        stop,
        seq.len(),
        hts_seq.len()
    );
}

/// Forked readers must produce identical results to the original at any position.
#[hegel::test(test_cases = 100)]
fn fork_matches_original_at_random_positions(tc: TestCase) {
    let start_frac = tc.draw(gs::floats::<f64>().min_value(0.0).max_value_exclusive(0.99));
    let len_frac = tc.draw(gs::floats::<f64>().min_value(0.001).max_value_exclusive(0.01));
    let length = 61_431_566u64; // chr19
    let start = (start_frac * length as f64) as u64;
    let fetch_len = ((len_frac * length as f64) as u64).max(1).min(length - start);
    let stop = start + fetch_len;

    let mut rio = IndexedFastaReader::open(test_fasta_path()).expect("rio open");
    let mut forked = rio.fork().expect("fork");

    let orig = rio.fetch_seq("chr19", span(start, stop)).expect("orig fetch");
    let fork_result = forked.fetch_seq("chr19", span(start, stop)).expect("fork fetch");

    assert_eq!(orig, fork_result, "fork mismatch at chr19:{}-{}", start, stop);
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
    let seq1 = readers.fetch_base_seq("bacteriophage_lambda_CpG", span(0, 100)).unwrap();
    assert_eq!(seq1.len(), 100);
    // Every element must be a valid Base
    for &b in seq1.iter() {
        assert!(matches!(b, Base::A | Base::C | Base::G | Base::T | Base::Unknown));
    }

    // Second call: should reuse the internal buffer (no way to observe capacity
    // directly, but we can verify correctness across calls)
    let seq2 = readers.fetch_base_seq("bacteriophage_lambda_CpG", span(100, 300)).unwrap();
    assert_eq!(seq2.len(), 200);

    // Third call: same region as first, must produce identical result
    let seq3 = readers.fetch_base_seq("bacteriophage_lambda_CpG", span(0, 100)).unwrap();
    assert_eq!(seq1, seq3, "repeated fetch must produce identical results");

    // Verify against raw fetch + manual conversion
    let mut raw_reader = IndexedFastaReader::open(test_fasta_path()).unwrap();
    let raw = raw_reader.fetch_seq("bacteriophage_lambda_CpG", span(0, 100)).unwrap();
    let expected: Vec<Base> = Base::from_ascii_vec(raw);
    assert_eq!(&*seq1, &expected[..], "fetch_base_seq must match from_ascii_vec on raw fetch");
}

// ---------------------------------------------------------------------------
// Fetch into buffer: must produce same result as allocating fetch
// ---------------------------------------------------------------------------

// r[verify fasta.fetch.buffer_reuse]

#[hegel::test(test_cases = 50)]
fn fetch_into_matches_fetch_at_random_positions(tc: TestCase) {
    let start_frac = tc.draw(gs::floats::<f64>().min_value(0.0).max_value_exclusive(0.99));
    let len_frac = tc.draw(gs::floats::<f64>().min_value(0.001).max_value_exclusive(0.01));
    let length = 48502u64; // bacteriophage_lambda_CpG
    let start = (start_frac * length as f64) as u64;
    let fetch_len = ((len_frac * length as f64) as u64).max(1).min(length - start);
    let stop = start + fetch_len;

    let mut rio = IndexedFastaReader::open(test_fasta_path()).expect("rio open");
    let mut buf = Vec::new();

    let region = span(start, stop);
    let alloc = rio.fetch_seq("bacteriophage_lambda_CpG", region).expect("fetch_seq");
    rio.fetch_seq_into("bacteriophage_lambda_CpG", region, &mut buf).expect("fetch_seq_into");

    assert_eq!(alloc.clone(), buf, "fetch vs fetch_into mismatch at {}-{}", start, stop);
    let hts = htslib_fetch("bacteriophage_lambda_CpG", start, stop);
    assert_eq!(alloc, hts, "seqair vs htslib mismatch at {}-{}", start, stop);
}
