//! Smoke tests for the example binaries.
//!
//! Each test builds and runs the example with the repo's test data,
//! asserting that it exits successfully and produces some output.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    reason = "test code"
)]

use std::process::Command;

fn cargo_bin() -> String {
    std::env::var("CARGO").unwrap_or_else(|_| "cargo".to_string())
}

/// Resolve a test-data path to an absolute path rooted at the workspace root.
fn test_data(relative: &str) -> String {
    let workspace_root =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap().parent().unwrap();
    workspace_root.join(relative).to_str().unwrap().to_string()
}

fn run_example(name: &str, args: &[&str]) -> std::process::Output {
    let output = Command::new(cargo_bin())
        .args(["run", "--example", name, "--"])
        .args(args)
        .output()
        .unwrap_or_else(|e| panic!("failed to run example {name}: {e}"));
    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        panic!("example {name} failed (exit {:?}):\n{stderr}", output.status.code());
    }
    output
}

const REGION: &str = "chr19:6103076-6103200";

fn assert_mpileup_tsv(stdout: &str, expected_contig: &str) {
    assert!(!stdout.is_empty(), "mpileup should produce output");
    let first_line = stdout.lines().next().unwrap();
    let fields: Vec<&str> = first_line.split('\t').collect();
    assert_eq!(fields.len(), 6, "mpileup should produce 6 TSV columns, got: {first_line}");
    assert_eq!(fields[0], expected_contig);
}

#[test]
fn example_mpileup() {
    let bam = test_data("tests/data/test.bam");
    let fasta = test_data("tests/data/test.fasta.gz");
    let output = run_example("mpileup", &[&bam, "-f", &fasta, "-r", REGION]);
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert_mpileup_tsv(&stdout, "chr19");
}

/// `chr19` (no start, no end) must select the whole contig.
#[test]
fn example_mpileup_region_chrom_only() {
    let bam = test_data("tests/data/test.bam");
    let fasta = test_data("tests/data/test.fasta.gz");
    let output = run_example("mpileup", &[&bam, "-f", &fasta, "-r", "chr19"]);
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert_mpileup_tsv(&stdout, "chr19");
}

/// `chr19:6103076` (start only, no end) must run from the start position to
/// the end of the contig — the case the old example panicked on with
/// "no end pos given".
#[test]
fn example_mpileup_region_start_only() {
    let bam = test_data("tests/data/test.bam");
    let fasta = test_data("tests/data/test.fasta.gz");
    let output = run_example("mpileup", &[&bam, "-f", &fasta, "-r", "chr19:6103076"]);
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert_mpileup_tsv(&stdout, "chr19");
    // First emitted position should be at or after the requested start.
    let first = stdout.lines().next().unwrap();
    let pos1: u32 = first.split('\t').nth(1).expect("pos column").parse().expect("pos is integer");
    assert!(pos1 >= 6_103_076, "first emitted 1-based pos {pos1} < requested start 6103076");
}

/// No `-r` argument: scan every contig in the header.
#[test]
fn example_mpileup_no_region() {
    let bam = test_data("tests/data/test.bam");
    let fasta = test_data("tests/data/test.fasta.gz");
    let output = run_example("mpileup", &[&bam, "-f", &fasta]);
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(!stdout.is_empty(), "whole-genome mpileup should produce output");
    // Output spans more than one contig in this BAM (chr19 + lambda + 2kb_3).
    let contigs: std::collections::BTreeSet<&str> =
        stdout.lines().map(|l| l.split('\t').next().unwrap()).collect();
    assert!(
        contigs.len() >= 2,
        "whole-genome mpileup should emit ≥2 contigs in this BAM, saw {contigs:?}"
    );
}

/// CRAM input with the same reference must produce byte-for-byte identical
/// output to the BAM input over the same region.
#[test]
fn example_mpileup_cram_matches_bam() {
    let bam = test_data("tests/data/test.bam");
    let cram = test_data("tests/data/test.cram");
    let fasta = test_data("tests/data/test.fasta.gz");

    let bam_output = run_example("mpileup", &[&bam, "-f", &fasta, "-r", REGION]);
    let cram_output = run_example("mpileup", &[&cram, "-f", &fasta, "-r", REGION]);

    let bam_stdout = String::from_utf8_lossy(&bam_output.stdout);
    let cram_stdout = String::from_utf8_lossy(&cram_output.stdout);

    assert_mpileup_tsv(&bam_stdout, "chr19");
    assert_mpileup_tsv(&cram_stdout, "chr19");
    assert_eq!(
        bam_stdout, cram_stdout,
        "CRAM and BAM must produce identical mpileup output for the same region + reference"
    );
}

#[test]
fn example_pileup_extras() {
    let bam = test_data("tests/data/test.bam");
    let fasta = test_data("tests/data/test.fasta.gz");
    let output = run_example("pileup_extras", &[&bam, &fasta, "-r", REGION]);
    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("Loaded pileup"), "should print region summary to stderr");
    // Header + data lines
    let lines: Vec<&str> = stdout.lines().collect();
    assert!(lines.len() > 1, "should produce header + data lines");
    assert!(lines[0].starts_with("pos\t"), "first line should be header");
    // Data line has 5 columns: pos, depth, ref, read_groups, mean_aligned_frac
    let data_fields: Vec<&str> = lines[1].split('\t').collect();
    assert_eq!(data_fields.len(), 5, "data lines should have 5 columns, got: {}", lines[1]);
}

#[test]
fn example_simple_variant_caller() {
    let bam = test_data("tests/data/test.bam");
    let fasta = test_data("tests/data/test.fasta.gz");
    let output = run_example(
        "simple_variant_caller",
        &[&bam, &fasta, "-r", REGION, "--min-depth", "2", "--min-af", "0.01"],
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    // VCF output starts with header lines
    assert!(stdout.starts_with("##"), "VCF output should start with ## header");
    assert!(stdout.contains("#CHROM"), "VCF output should contain #CHROM header line");
}

#[test]
fn example_realignment() {
    let bam = test_data("tests/data/test.bam");
    let output = run_example("realignment", &[&bam, "-r", REGION]);
    let stdout = String::from_utf8_lossy(&output.stdout);
    // Realignment prints summary lines
    assert!(!stdout.is_empty(), "realignment should produce output");
}

#[test]
fn example_base_mods() {
    let bam = test_data("tests/data/test.bam");
    // base_mods needs a region with MM/ML tags — the test BAM has them
    // on the bacteriophage_lambda_CpG contig.
    let output = run_example("base_mods", &[&bam, "-r", "bacteriophage_lambda_CpG:1-500"]);
    // Just verify successful exit (handled by run_example).
    let _ = output;
}
