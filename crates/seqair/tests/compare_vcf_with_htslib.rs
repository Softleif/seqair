//! Validate seqair BCF output using bcftools (htslib reference implementation).
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects
)]
//!
//! Writes BCF with seqair, runs `bcftools view` to verify the file is valid
//! and fields are parsed correctly. This is the strongest validation available:
//! bcftools IS the reference BCF implementation.

use proptest::prelude::*;
use seqair::vcf::alleles::Alleles;
use seqair::vcf::bcf_writer::BcfWriter;
use seqair::vcf::header::{ContigDef, FormatDef, InfoDef, Number, ValueType, VcfHeader};
use seqair::vcf::record::{Genotype, SampleValue, VcfRecordBuilder};
use seqair_types::{Base, One, Pos, SmolStr};
use std::sync::Arc;

fn shared_header() -> Arc<VcfHeader> {
    Arc::new(
        VcfHeader::builder()
            .add_contig("chr1", ContigDef { length: Some(250_000_000) })
            .unwrap()
            .add_info(
                "DP",
                InfoDef {
                    number: Number::Count(1),
                    typ: ValueType::Integer,
                    description: SmolStr::from("Total Depth"),
                },
            )
            .unwrap()
            .add_format(
                "GT",
                FormatDef {
                    number: Number::Count(1),
                    typ: ValueType::String,
                    description: SmolStr::from("Genotype"),
                },
            )
            .unwrap()
            .add_format(
                "DP",
                FormatDef {
                    number: Number::Count(1),
                    typ: ValueType::Integer,
                    description: SmolStr::from("Read Depth"),
                },
            )
            .unwrap()
            .add_sample("sample1")
            .unwrap()
            .build()
            .unwrap(),
    )
}

/// Write a BCF file with seqair to a temp file.
#[allow(clippy::too_many_arguments)]
fn write_seqair_bcf(
    header: &Arc<VcfHeader>,
    pos: u32,
    ref_base: Base,
    alt_base: Base,
    qual: f32,
    depth: i32,
    gt_a0: u16,
    gt_a1: u16,
    phased: bool,
) -> tempfile::NamedTempFile {
    let tmp = tempfile::Builder::new().suffix(".bcf").tempfile().unwrap();
    let alleles = Alleles::snv(ref_base, alt_base).unwrap();
    let gt = if phased {
        Genotype::phased_diploid(gt_a0, gt_a1)
    } else {
        Genotype::unphased(gt_a0, gt_a1)
    };
    let record = VcfRecordBuilder::new("chr1", Pos::<One>::new(pos).unwrap(), alleles)
        .qual(qual)
        .filter_pass()
        .info_integer("DP", depth)
        .format_keys(&["GT", "DP"])
        .add_sample(vec![SampleValue::Genotype(gt), SampleValue::Integer(depth)])
        .build(header)
        .unwrap();

    let mut buf = Vec::new();
    {
        let mut writer = BcfWriter::new(&mut buf, header.clone(), false);
        writer.write_header().unwrap();
        writer.write_vcf_record(&record).unwrap();
        writer.finish().unwrap();
    }
    std::fs::write(tmp.path(), &buf).unwrap();
    tmp
}

/// Run `bcftools view` on a BCF file and return the VCF text output.
/// Returns None if bcftools is not available.
fn bcftools_view(path: &std::path::Path) -> Option<String> {
    let output = std::process::Command::new("bcftools")
        .args(["view", "-H"]) // -H = no header, just data lines
        .arg(path)
        .output()
        .ok()?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        panic!("bcftools view failed: {stderr}");
    }

    Some(String::from_utf8(output.stdout).unwrap())
}

/// Check if bcftools is available.
fn has_bcftools() -> bool {
    std::process::Command::new("bcftools").arg("--version").output().is_ok()
}

// ── Deterministic tests ────────────────────────────────────────────────

#[test]
fn bcftools_reads_seqair_simple_snv() {
    if !has_bcftools() {
        eprintln!("skipping: bcftools not found");
        return;
    }

    let header = shared_header();
    let tmp = write_seqair_bcf(&header, 12345, Base::A, Base::T, 30.0, 50, 0, 1, false);

    let vcf_text = bcftools_view(tmp.path()).unwrap();
    let fields: Vec<&str> = vcf_text.trim().split('\t').collect();

    assert_eq!(fields[0], "chr1", "CHROM");
    assert_eq!(fields[1], "12345", "POS");
    assert_eq!(fields[3], "A", "REF");
    assert_eq!(fields[4], "T", "ALT");
    assert_eq!(fields[5], "30", "QUAL");
    assert_eq!(fields[6], "PASS", "FILTER");
    assert!(fields[7].contains("DP=50"), "INFO should contain DP=50, got: {}", fields[7]);
}

#[test]
fn bcftools_reads_seqair_phased_gt() {
    if !has_bcftools() {
        eprintln!("skipping: bcftools not found");
        return;
    }

    let header = shared_header();
    let tmp = write_seqair_bcf(&header, 1000, Base::C, Base::G, 99.0, 100, 0, 1, true);

    let vcf_text = bcftools_view(tmp.path()).unwrap();
    let fields: Vec<&str> = vcf_text.trim().split('\t').collect();

    // GT should show phased separator
    assert!(fields[9].starts_with("0|1"), "expected phased GT 0|1, got: {}", fields[9]);
}

#[test]
fn bcftools_reads_seqair_hom_ref() {
    if !has_bcftools() {
        eprintln!("skipping: bcftools not found");
        return;
    }

    let header = shared_header();
    let tmp = write_seqair_bcf(&header, 500, Base::G, Base::A, 10.0, 5, 0, 0, false);

    let vcf_text = bcftools_view(tmp.path()).unwrap();
    let fields: Vec<&str> = vcf_text.trim().split('\t').collect();

    assert!(fields[9].starts_with("0/0"), "expected hom ref GT 0/0, got: {}", fields[9]);
}

// ── Proptest: bcftools validates random seqair BCF ─────────────────────

fn arb_base() -> impl Strategy<Value = Base> {
    prop_oneof![Just(Base::A), Just(Base::C), Just(Base::G), Just(Base::T)]
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(20))]

    /// Write random SNVs with seqair, validate with bcftools (the reference implementation).
    #[test]
    fn bcftools_validates_random_seqair_bcf(
        pos in 1u32..10_000_000,
        ref_base in arb_base(),
        alt_base in arb_base(),
        qual in 1.0f32..1000.0,
        depth in 1i32..10000,
        gt_a0 in 0u16..2,
        gt_a1 in 0u16..2,
        phased in proptest::bool::ANY,
    ) {
        prop_assume!(ref_base != alt_base);

        if !has_bcftools() {
            return Ok(());
        }

        let header = shared_header();
        let tmp = write_seqair_bcf(
            &header, pos, ref_base, alt_base, qual, depth, gt_a0, gt_a1, phased,
        );

        // bcftools must parse without error and produce correct POS
        let vcf_text = bcftools_view(tmp.path()).unwrap();
        let fields: Vec<&str> = vcf_text.trim().split('\t').collect();

        prop_assert_eq!(fields[0], "chr1");
        prop_assert_eq!(fields[1], &pos.to_string());
        prop_assert_eq!(fields[3], &format!("{}", ref_base.as_char()));
        prop_assert_eq!(fields[4], &format!("{}", alt_base.as_char()));
    }
}
