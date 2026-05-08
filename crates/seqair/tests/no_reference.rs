//! Reading BAM/SAM/CRAM without a FASTA reference.
//!
//! Covers `r[unified.readers_open]` (the no-FASTA half of the polymorphic
//! signature), `r[cram.fasta.optional]`, and `r[pileup.reference_base.optional]`.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "test code"
)]
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    reason = "test code with known small values"
)]

use seqair::bam::{Pos0, RecordStore};
use seqair::cram::reader::CramError;
use seqair::reader::{ReaderError, Readers, SegmentOptions};
use seqair_types::Base;
use std::num::NonZeroU32;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn test_bam_path() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.bam"))
}

fn test_cram_v30_path() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test_v30.cram"))
}

fn test_fasta_path() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.fasta.gz"))
}

fn htslib_sam(name: &str) -> PathBuf {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/htslib/sam/")).join(name)
}

fn htslib_fasta(name: &str) -> PathBuf {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/htslib/fasta/")).join(name)
}

/// Convert a SAM to CRAM, optionally embedding the reference.
/// Mirrors `tests/cram_version_matrix.rs::sam_to_cram`.
fn sam_to_cram(
    dir: &Path,
    sam_path: &Path,
    fasta_path: &Path,
    version: &str,
    extra_opts: &[&str],
) -> PathBuf {
    let bam_path = dir.join("sorted.bam");
    let cram_path = dir.join(format!("test_{version}.cram"));

    let status = Command::new("samtools")
        .args(["sort", "-o"])
        .arg(&bam_path)
        .arg(sam_path)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .expect("samtools not found");
    assert!(status.success(), "samtools sort failed for {}", sam_path.display());

    let mut cmd = Command::new("samtools");
    cmd.args(["view", "-C"]).arg("--output-fmt-option").arg(format!("version={version}"));
    for opt in extra_opts {
        cmd.arg("--output-fmt-option").arg(opt);
    }
    cmd.arg("-T").arg(fasta_path).arg("-o").arg(&cram_path).arg(&bam_path);
    let status =
        cmd.stdout(Stdio::null()).stderr(Stdio::null()).status().expect("samtools view -C");
    assert!(status.success(), "samtools view -C failed");

    let status = Command::new("samtools")
        .arg("index")
        .arg(&cram_path)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .expect("samtools index");
    assert!(status.success(), "samtools index failed");

    cram_path
}

// r[verify unified.readers_open]
/// Opening a BAM without a FASTA via `Readers::open(path, None)` must succeed
/// and let `fetch_into` work.
#[test]
fn bam_open_without_reference_fetches_records() {
    let mut readers = Readers::open(test_bam_path(), None).expect("open BAM without FASTA");
    let tid = readers.header().tid("chr19").expect("test BAM has chr19");
    let mut store = RecordStore::new();
    let count = readers
        .fetch_into(tid, Pos0::new(6_100_000).unwrap(), Pos0::new(6_200_000).unwrap(), &mut store)
        .expect("fetch_into without FASTA");
    assert!(count > 0, "test BAM should yield records in chr19:6.1-6.2M");
    assert_eq!(store.len(), count);
}

// r[verify pileup.reference_base.optional]
/// Pileup without a FASTA must succeed; every column reports `Base::Unknown`
/// for the reference base, but per-read bases are still present.
#[test]
fn bam_pileup_without_reference_returns_unknown_ref_base() {
    let mut readers = Readers::open(test_bam_path(), None).expect("open BAM without FASTA");
    let opts = SegmentOptions::new(NonZeroU32::new(3_000).unwrap());
    let segments: Vec<_> = readers
        .segments(("chr19", Pos0::new(6_103_500).unwrap(), Pos0::new(6_106_500).unwrap()), opts)
        .expect("segments")
        .collect();
    assert!(!segments.is_empty(), "test BAM should yield at least one segment");

    let mut total_columns = 0usize;
    let mut total_alignments = 0usize;
    {
        let mut p = readers.pileup(&segments[0]).expect("pileup without FASTA");
        while let Some(col) = p.pileups() {
            assert_eq!(
                col.reference_base(),
                Base::Unknown,
                "pileup without FASTA must report reference_base = Unknown"
            );
            total_columns += 1;
            total_alignments += col.alignments().count();
        }
    }
    assert!(total_columns > 0, "expected at least one pileup column");
    assert!(total_alignments > 0, "pileup columns should still expose per-read alignments");
}

// r[verify cram.fasta.optional]
/// CRAM with embedded reference (RR=true but `embedded_reference >= 0` per
/// slice) MUST decode without an external FASTA. samtools writes such files
/// with `--output-fmt-option embed_ref=2`.
#[test]
fn cram_embedded_reference_no_external_fasta() {
    let dir = tempfile::tempdir().expect("tempdir");
    let cram = sam_to_cram(
        dir.path(),
        &htslib_sam("c1#clip.sam"),
        &htslib_fasta("c1.fa"),
        "3.0",
        &["embed_ref=2"],
    );

    let fasta = htslib_fasta("c1.fa");
    let with_fasta_count = {
        let mut readers = Readers::open(&cram, fasta.as_path()).expect("open with FASTA");
        let tid = readers.header().tid("c1").expect("tid");
        let mut store = RecordStore::new();
        readers
            .fetch_into(tid, Pos0::new(0).unwrap(), Pos0::new(10).unwrap(), &mut store)
            .expect("fetch with FASTA");
        store.len()
    };

    let mut readers =
        Readers::open(cram.as_path(), None).expect("open embed-ref CRAM without FASTA");
    let tid = readers.header().tid("c1").expect("tid");
    let mut store = RecordStore::new();
    readers
        .fetch_into(tid, Pos0::new(0).unwrap(), Pos0::new(10).unwrap(), &mut store)
        .expect("fetch embed-ref CRAM without FASTA");
    assert!(!store.is_empty(), "embed-ref CRAM should still yield records without FASTA");
    assert_eq!(
        store.len(),
        with_fasta_count,
        "embed-ref CRAM record count should match between with-FASTA and without-FASTA"
    );
}

// r[verify cram.edge.missing_reference]
/// CRAM with a normal external reference (RR=true, no embedded ref) MUST
/// surface a `MissingReference` error when no FASTA was supplied.
#[test]
fn cram_external_reference_required_errors_without_fasta() {
    let mut readers = Readers::open(test_cram_v30_path(), None)
        .expect("open external-ref CRAM without FASTA (open is fine, only fetch errors)");
    let tid = readers.header().tid("chr19").expect("tid");
    let mut store = RecordStore::new();
    let err = readers
        .fetch_into(tid, Pos0::new(6_100_000).unwrap(), Pos0::new(6_200_000).unwrap(), &mut store)
        .expect_err("fetch must fail when CRAM needs an external reference and we have none");
    match err {
        ReaderError::Cram { source: CramError::MissingReference { contig } } => {
            assert_eq!(contig.as_str(), "chr19");
        }
        other => panic!("expected CramError::MissingReference, got {other:?}"),
    }
}

// r[verify unified.readers_backward_compat]
/// `Readers::open` (with FASTA) must still work — backwards compatibility.
#[test]
fn readers_open_with_fasta_still_works() {
    let mut readers = Readers::open(test_bam_path(), test_fasta_path()).expect("open with FASTA");
    let tid = readers.header().tid("chr19").expect("tid");
    let mut store = RecordStore::new();
    let count = readers
        .fetch_into(tid, Pos0::new(6_100_000).unwrap(), Pos0::new(6_200_000).unwrap(), &mut store)
        .expect("fetch with FASTA");
    assert!(count > 0);
}

// r[verify unified.readers_fork]
/// Forking a `Readers` opened without a FASTA must also produce a fork
/// without a FASTA. A subsequent `pileup` on the fork must work the same way.
#[test]
fn fork_without_reference_keeps_no_reference() {
    let readers = Readers::open(test_bam_path(), None).expect("open without FASTA");
    let mut forked = readers.fork().expect("fork preserves no-reference");
    let opts = SegmentOptions::new(NonZeroU32::new(3_000).unwrap());
    let segments: Vec<_> = forked
        .segments(("chr19", Pos0::new(6_103_500).unwrap(), Pos0::new(6_106_500).unwrap()), opts)
        .expect("segments on forked reader")
        .collect();
    assert!(!segments.is_empty());
    let mut p = forked.pileup(&segments[0]).expect("pileup on forked reader without FASTA");
    while let Some(col) = p.pileups() {
        assert_eq!(col.reference_base(), Base::Unknown);
    }
}
