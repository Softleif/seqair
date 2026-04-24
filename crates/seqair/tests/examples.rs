//! Smoke tests for the example binaries.
//!
//! Each test builds and runs the example with the repo's test data,
//! asserting that it exits successfully and produces some output.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, reason = "test code")]

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

#[test]
fn example_mpileup() {
    let bam = test_data("tests/data/test.bam");
    let fasta = test_data("tests/data/test.fasta.gz");
    let output = run_example("mpileup", &[&bam, "-f", &fasta, "-r", REGION]);
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(!stdout.is_empty(), "mpileup should produce output");
    // Output is TSV with contig, pos, ref, depth, bases, quals
    let first_line = stdout.lines().next().unwrap();
    let fields: Vec<&str> = first_line.split('\t').collect();
    assert_eq!(fields.len(), 6, "mpileup should produce 6 TSV columns, got: {first_line}");
    assert_eq!(fields[0], "chr19");
}

#[test]
fn example_pileup_extras() {
    let bam = test_data("tests/data/test.bam");
    let fasta = test_data("tests/data/test.fasta.gz");
    let output = run_example("pileup_extras", &[&bam, &fasta, "-r", REGION]);
    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("Loaded"), "should print record count to stderr");
    // Header + data lines
    let lines: Vec<&str> = stdout.lines().collect();
    assert!(lines.len() > 1, "should produce header + data lines");
    assert!(lines[0].starts_with("pos\t"), "first line should be header");
    // Data line has 6 columns
    let data_fields: Vec<&str> = lines[1].split('\t').collect();
    assert_eq!(data_fields.len(), 6, "data lines should have 6 columns, got: {}", lines[1]);
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
