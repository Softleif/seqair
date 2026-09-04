//! Indexed access to FASTQ, cross-checked against htslib.
//!
//! A FASTQ index is not a distinct format: `samtools faidx` and `samtools
//! fqidx` both write `<file>.fai`, carrying the five FASTA columns plus a
//! sixth, `qual_offset`. Everything seqair needs for coordinate queries lives
//! in the first five, and the byte-span math is bounded by `length`, so the
//! reader never walks into the `+` separator or the quality block.
//!
//! That bound is what makes indexed access safe on input no line-oriented
//! parser survives. `@` and `+` are valid quality characters (Phred 31 and
//! 10), so a wrapped record's quality block may begin lines with either
//! sigil. Every Rust FASTQ crate surveyed as of 2026-09 gets this wrong:
//! `noodles-fastq`, `seq_io`, `needletail`, `paraseq` and `fastx` all assume
//! exactly four lines per record and error on wrapped input (noodles' own
//! `fastq::fs::index` cannot index the very files its `quality_scores_offset`
//! column describes), while `helicase` accepts it and mis-parses silently,
//! taking a `@`-prefixed quality line as the next record's header. htslib is
//! the only prior implementation that handles it, which is why it is the
//! oracle here.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    reason = "test code"
)]
#![allow(clippy::arithmetic_side_effects, reason = "test code")]
#![allow(clippy::unwrap_in_result, reason = "test helper propagates only the parse error")]
#![allow(clippy::cast_possible_truncation, reason = "test code with known small values")]

use proptest::prelude::*;
use rust_htslib::faidx;
use seqair::bam::Pos0;
use seqair::fasta::{FaiEntry, FaiEntryError, FaiError, FastaIndex, IndexedFastaReader};
use std::io::Write;
use std::path::{Path, PathBuf};
use tempfile::TempDir;

/// One generated FASTQ record: name, bases, and the width its lines wrap at.
#[derive(Debug, Clone)]
struct Rec {
    name: String,
    seq: Vec<u8>,
    qual: Vec<u8>,
    linebases: usize,
}

/// Write records as a wrapped FASTQ and build its index with htslib.
///
/// Returns the temp dir (kept alive by the caller) and the FASTQ path. The
/// `.fai` htslib writes is the real 6-column article, not a hand-rolled one.
fn write_and_index(recs: &[Rec]) -> (TempDir, PathBuf) {
    let dir = TempDir::new().unwrap();
    let path = dir.path().join("reads.fq");
    {
        let mut f = std::io::BufWriter::new(std::fs::File::create(&path).unwrap());
        for r in recs {
            writeln!(f, "@{}", r.name).unwrap();
            for chunk in r.seq.chunks(r.linebases) {
                f.write_all(chunk).unwrap();
                f.write_all(b"\n").unwrap();
            }
            f.write_all(b"+\n").unwrap();
            // The quality block must repeat the sequence block's line
            // geometry; that is what `samtools fqidx` assumes.
            for chunk in r.qual.chunks(r.linebases) {
                f.write_all(chunk).unwrap();
                f.write_all(b"\n").unwrap();
            }
        }
        f.flush().unwrap();
    }
    faidx::build(&path).expect("htslib must be able to index generated FASTQ");
    (dir, path)
}

/// Parse `.fai` text through the real file-based entry point.
fn parse_fai(contents: &str) -> Result<FastaIndex, FaiError> {
    let dir = TempDir::new().unwrap();
    let path = dir.path().join("x.fa.fai");
    std::fs::write(&path, contents).unwrap();
    FastaIndex::from_file(&path)
}

fn fai_path(path: &Path) -> PathBuf {
    let mut p = path.as_os_str().to_owned();
    p.push(".fai");
    PathBuf::from(p)
}

fn fai_text(path: &Path) -> String {
    std::fs::read_to_string(fai_path(path)).unwrap()
}

// ---------------------------------------------------------------------------
// Index parsing
// ---------------------------------------------------------------------------

// r[verify fastq.access.index_parse]
#[test]
fn six_column_index_parses_and_keeps_qual_offset() {
    let contents = "chr1\t20\t6\t10\t11\t30\nchr2\t10\t58\t10\t11\t71\n";
    let idx = parse_fai(contents).unwrap();

    let chr1 = idx.get("chr1").unwrap();
    assert_eq!(chr1.length, 20);
    assert_eq!(chr1.offset, 6);
    assert_eq!(chr1.qual_offset, Some(30));
    assert!(chr1.is_fastq());

    let chr2 = idx.get("chr2").unwrap();
    assert_eq!(chr2.qual_offset, Some(71));
}

// r[verify fastq.access.index_parse]
#[test]
fn five_column_index_still_parses_as_fasta() {
    let idx = parse_fai("chr1\t20\t6\t10\t11\n").unwrap();
    let chr1 = idx.get("chr1").unwrap();
    assert_eq!(chr1.qual_offset, None);
    assert!(!chr1.is_fastq());
}

// r[verify fastq.access.index_parse]
#[test]
fn seven_columns_is_still_rejected() {
    let err = parse_fai("chr1\t20\t6\t10\t11\t30\t99\n").unwrap_err();
    assert!(
        matches!(
            err,
            FaiError::InvalidEntry { kind: FaiEntryError::TooManyFields { found: 7 }, .. }
        ),
        "expected TooManyFields, got {err:?}"
    );
}

// r[verify fastq.access.index_parse]
#[test]
fn non_numeric_qual_offset_is_rejected() {
    let err = parse_fai("chr1\t20\t6\t10\t11\tNOPE\n").unwrap_err();
    assert!(
        matches!(
            err,
            FaiError::InvalidEntry {
                kind: FaiEntryError::InvalidField { field: "qual_offset", .. },
                ..
            }
        ),
        "expected InvalidField(qual_offset), got {err:?}"
    );
}

// ---------------------------------------------------------------------------
// The case that defeats every line-oriented parser
// ---------------------------------------------------------------------------

// r[verify fastq.access.indexed_seq]
// r[verify fastq.multiline.ambiguity]
#[test]
fn quality_lines_beginning_with_at_and_plus_do_not_confuse_indexed_access() {
    // chr1's quality block is `@@@@@@@@@@` then `++++++++++` — Phred 31 and 10,
    // both entirely legal, both indistinguishable from a header or separator to
    // a parser that scans for sigils.
    let recs = vec![
        Rec {
            name: "chr1".into(),
            seq: b"ACGTACGTACGTTTTTACGT".to_vec(),
            qual: b"@@@@@@@@@@++++++++++".to_vec(),
            linebases: 10,
        },
        Rec {
            name: "chr2".into(),
            seq: b"TTTTGGGGCC".to_vec(),
            qual: b"!!!!!!!!!!".to_vec(),
            linebases: 10,
        },
    ];
    let (_dir, path) = write_and_index(&recs);

    // htslib produced a 6-column index and seqair now accepts it directly.
    assert_eq!(fai_text(&path).lines().next().unwrap().split('\t').count(), 6);

    let mut reader = IndexedFastaReader::open(&path).unwrap();
    for r in &recs {
        let got = reader
            .fetch_seq(&r.name, Pos0::new(0).unwrap(), Pos0::new(r.seq.len() as u32).unwrap())
            .unwrap();
        assert_eq!(got, r.seq, "{} misread", r.name);
    }
}

// ---------------------------------------------------------------------------
// Differential properties against htslib
// ---------------------------------------------------------------------------

/// Bases include lowercase so the uppercasing rule is exercised, and `N` so
/// the generator is not limited to unambiguous calls.
fn base_strategy() -> impl Strategy<Value = u8> {
    prop::sample::select(vec![b'A', b'C', b'G', b'T', b'N', b'a', b'c', b'g', b't', b'n'])
}

/// The full Sanger quality range, `!` (33) through `~` (126) — `@` and `+`
/// included, which is the whole point.
fn qual_strategy() -> impl Strategy<Value = u8> {
    (33u8..=126u8).boxed()
}

fn rec_strategy(idx: usize) -> impl Strategy<Value = Rec> {
    (1usize..=200, 1usize..=60).prop_flat_map(move |(len, linebases)| {
        (
            prop::collection::vec(base_strategy(), len..=len),
            prop::collection::vec(qual_strategy(), len..=len),
            Just(linebases.min(len)),
        )
            .prop_map(move |(seq, qual, linebases)| Rec {
                name: format!("ctg{idx}"),
                seq,
                qual,
                linebases,
            })
    })
}

fn recs_strategy() -> impl Strategy<Value = Vec<Rec>> {
    (1usize..=3).prop_flat_map(|n| {
        let parts: Vec<_> = (0..n).map(rec_strategy).collect();
        parts
    })
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 64, ..ProptestConfig::default() })]

    // r[verify fastq.access.indexed_seq]
    /// Every subrange of every record must match what htslib returns from the
    /// same file and the same index — modulo case, which seqair normalises per
    /// `r[fasta.fetch.uppercase]` and htslib does not.
    #[test]
    fn indexed_fetch_matches_htslib(recs in recs_strategy()) {
        let (_dir, path) = write_and_index(&recs);
        let mut seqair_reader = IndexedFastaReader::open(&path).unwrap();
        let htslib_reader = faidx::Reader::from_path(&path).unwrap();

        for r in &recs {
            let len = r.seq.len();
            // Whole record, a prefix, a suffix, and an interior window.
            let ranges = [
                (0usize, len),
                (0, len.div_ceil(2)),
                (len / 2, len),
                (len / 3, (2 * len) / 3),
            ];
            for (start, stop) in ranges {
                if start >= stop {
                    continue;
                }
                let got = seqair_reader
                    .fetch_seq(
                        &r.name,
                        Pos0::new(start as u32).unwrap(),
                        Pos0::new(stop as u32).unwrap(),
                    )
                    .unwrap();

                let mut want = htslib_reader.fetch_seq(&r.name, start, stop - 1).unwrap();
                want.make_ascii_uppercase();
                prop_assert_eq!(&got, &want, "{} [{}, {})", r.name, start, stop);

                // Independent of htslib: it must also equal the generated bases.
                let mut expected = r.seq[start..stop].to_vec();
                expected.make_ascii_uppercase();
                prop_assert_eq!(&got, &expected, "{} [{}, {}) vs source", r.name, start, stop);
            }
        }
    }

    // r[verify fastq.access.index_parse]
    /// seqair's parse of htslib's own index must agree with htslib's geometry:
    /// same sequence names, lengths, and a `qual_offset` on every entry.
    #[test]
    fn parsed_index_matches_htslib_geometry(recs in recs_strategy()) {
        let (_dir, path) = write_and_index(&recs);
        let idx = FastaIndex::from_file(&fai_path(&path)).unwrap();
        let htslib_reader = faidx::Reader::from_path(&path).unwrap();

        prop_assert_eq!(idx.len(), recs.len());
        for r in &recs {
            let e = idx.get(&r.name).expect("name present");
            prop_assert_eq!(e.length, r.seq.len() as u64);
            prop_assert_eq!(e.length, htslib_reader.fetch_seq_len(&r.name));
            prop_assert!(e.is_fastq(), "FASTQ index entry must carry qual_offset");
        }
    }

    // r[verify fastq.access.indexed_qual]
    /// `qual_byte_offset` must land on the right quality character. The oracle
    /// is the raw file: read the byte at the computed offset and compare it to
    /// the quality value that was generated for that position.
    #[test]
    fn qual_byte_offset_locates_the_right_character(recs in recs_strategy()) {
        let (_dir, path) = write_and_index(&recs);
        let raw = std::fs::read(&path).unwrap();
        let idx = FastaIndex::from_file(&fai_path(&path)).unwrap();

        for r in &recs {
            let e = idx.get(&r.name).unwrap();
            for pos in 0..r.qual.len() {
                let off = e.qual_byte_offset(pos as u64).expect("FASTQ entry") as usize;
                prop_assert_eq!(
                    raw.get(off).copied(),
                    Some(r.qual[pos]),
                    "{} qual pos {}", r.name, pos
                );
            }
        }
    }
}

// r[verify fastq.access.indexed_qual]
#[test]
fn qual_byte_offset_is_none_for_a_fasta_index() {
    let e = FaiEntry {
        name: "s".into(),
        length: 100,
        offset: 10,
        linebases: 50,
        linewidth: 51,
        qual_offset: None,
    };
    assert_eq!(e.qual_byte_offset(0), None);
    assert!(!e.is_fastq());
}

// ---------------------------------------------------------------------------
// End to end: a FASTQ reference through `Readers::open` / `segments` / `pileup`
// ---------------------------------------------------------------------------

/// Re-emit a FASTA as a FASTQ with the same names, bases and line geometry,
/// filling the quality block with characters drawn from the full Sanger range
/// so `@` and `+` lead some lines.
fn fasta_to_fastq(fasta: &Path, out: &Path) {
    let text = std::fs::read_to_string(fasta).unwrap();
    let mut f = std::io::BufWriter::new(std::fs::File::create(out).unwrap());
    let mut widths: Vec<usize> = Vec::new();
    let mut next_qual = 0u8;
    let mut flush = |f: &mut std::io::BufWriter<std::fs::File>, widths: &mut Vec<usize>| {
        if widths.is_empty() {
            return;
        }
        f.write_all(b"+\n").unwrap();
        for w in widths.iter() {
            let line: Vec<u8> = (0..*w)
                .map(|_| {
                    next_qual = next_qual.wrapping_add(7) % 94;
                    33 + next_qual
                })
                .collect();
            f.write_all(&line).unwrap();
            f.write_all(b"\n").unwrap();
        }
        widths.clear();
    };
    for line in text.lines() {
        if let Some(name) = line.strip_prefix('>') {
            flush(&mut f, &mut widths);
            writeln!(f, "@{name}").unwrap();
        } else {
            widths.push(line.len());
            f.write_all(line.as_bytes()).unwrap();
            f.write_all(b"\n").unwrap();
        }
    }
    flush(&mut f, &mut widths);
    f.flush().unwrap();
}

// r[verify fastq.access.indexed_seq]
/// The request this whole thing came from: hand `Readers::open` a FASTQ
/// reference next to an alignment file and get identical pileup columns to the
/// FASTA it was converted from.
#[test]
fn readers_pileup_with_a_fastq_reference_matches_the_fasta() {
    use seqair::Readers;
    use seqair::bam::DepthLimit;
    use seqair::reader::SegmentOptions;

    let data = Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data"));
    let cram = data.join("unplaced_multiref.cram");
    let fasta = data.join("unplaced_multiref.fa");

    let dir = TempDir::new().unwrap();
    let fastq = dir.path().join("unplaced_multiref.fq");
    fasta_to_fastq(&fasta, &fastq);
    faidx::build(&fastq).expect("htslib must index the converted reference");
    assert_eq!(
        fai_text(&fastq).lines().next().unwrap().split('\t').count(),
        6,
        "converted reference must carry a 6-column FASTQ index"
    );

    let reference_bases = |reference: &Path| -> Vec<(u64, seqair_types::Base)> {
        let mut readers = Readers::open(&cram, reference).unwrap();
        let mut out = Vec::new();
        let segments: Vec<_> =
            readers.segments("chr1", SegmentOptions::default()).unwrap().collect();
        for segment in segments {
            let mut pileup = readers.pileup(&segment, DepthLimit::Unlimited).unwrap();
            while let Some(col) = pileup.pileups() {
                out.push((col.pos().as_u64(), col.reference_base()));
            }
        }
        out
    };

    let from_fasta = reference_bases(&fasta);
    let from_fastq = reference_bases(&fastq);

    assert!(!from_fasta.is_empty(), "expected pileup columns from the test CRAM");
    assert_eq!(from_fastq, from_fasta, "a FASTQ reference must give identical columns");
}
