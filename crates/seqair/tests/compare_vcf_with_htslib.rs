//! Validate seqair BCF output using bcftools (htslib reference implementation).
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "test code"
)]
//!
//! Writes BCF with seqair, runs `bcftools view` to verify the file is valid
//! and fields are parsed correctly. This is the strongest validation available:
//! bcftools IS the reference BCF implementation.

use proptest::prelude::*;
use seqair::vcf::alleles::Alleles;
use seqair::vcf::header::{ContigDef, Number, ValueType};
use seqair::vcf::record::Genotype;
use seqair::vcf::record_encoder::{FormatFieldDef, Gt, InfoFieldDef, Scalar};
use seqair::vcf::{ContigId, FormatGt, FormatInt, InfoInt, OutputFormat, VcfHeader, Writer};
use seqair_types::{Base, Pos1};
use std::sync::Arc;

struct TestSetup {
    header: Arc<VcfHeader>,
    contig: ContigId,
    dp_info: InfoInt,
    gt_fmt: FormatGt,
    dp_fmt: FormatInt,
}

fn shared_setup() -> TestSetup {
    let mut builder = VcfHeader::builder();
    let contig = builder.register_contig("chr1", ContigDef { length: Some(250_000_000) }).unwrap();
    let mut builder = builder.infos();
    let dp_info = builder
        .register_info(&InfoFieldDef::new(
            "DP",
            Number::Count(1),
            ValueType::Integer,
            "Total Depth",
        ))
        .unwrap();
    let mut builder = builder.formats();
    let gt_fmt = builder
        .register_format(&FormatFieldDef::<Gt>::new(
            "GT",
            Number::Count(1),
            ValueType::String,
            "Genotype",
        ))
        .unwrap();
    let dp_fmt = builder
        .register_format(&FormatFieldDef::<Scalar<i32>>::new(
            "DP",
            Number::Count(1),
            ValueType::Integer,
            "Read Depth",
        ))
        .unwrap();
    let mut builder = builder.samples();
    builder.add_sample("sample1").unwrap();
    let header = Arc::new(builder.build().unwrap());
    TestSetup { header, contig, dp_info, gt_fmt, dp_fmt }
}

/// Write a BCF file with seqair to a temp file.
#[allow(clippy::too_many_arguments, reason = "test helper with many configurable fields")]
fn write_seqair_bcf(
    setup: &TestSetup,
    pos: u32,
    ref_base: Base,
    alt_base: Base,
    qual: f32,
    depth: i32,
    gt_a0: u16,
    gt_a1: u16,
    phased: bool,
) -> tempfile::NamedTempFile {
    let gt = if phased {
        Genotype::phased_diploid(gt_a0, gt_a1)
    } else {
        Genotype::unphased(gt_a0, gt_a1)
    };
    write_seqair_bcf_with_gt(setup, pos, ref_base, alt_base, qual, depth, gt)
}

/// As `write_seqair_bcf`, with the genotype supplied directly — for genotypes
/// the diploid constructors cannot express.
#[allow(clippy::too_many_arguments, reason = "test helper mirroring the record shape")]
fn write_seqair_bcf_with_gt(
    setup: &TestSetup,
    pos: u32,
    ref_base: Base,
    alt_base: Base,
    qual: f32,
    depth: i32,
    gt: Genotype,
) -> tempfile::NamedTempFile {
    let tmp = tempfile::Builder::new().suffix(".bcf").tempfile().unwrap();
    let alleles = Alleles::snv(ref_base, alt_base).unwrap();

    let mut buf = Vec::new();
    {
        let writer = Writer::new(&mut buf, OutputFormat::Bcf);
        let mut writer = writer.write_header(&setup.header).unwrap();
        let mut enc = writer
            .begin_record(&setup.contig, Pos1::new(pos).unwrap(), &alleles, Some(qual))
            .unwrap()
            .filter_pass();
        setup.dp_info.encode(&mut enc, depth);
        let mut enc = enc.begin_samples();
        setup.gt_fmt.encode(&mut enc, &[gt]).unwrap();
        setup.dp_fmt.encode(&mut enc, &[depth]).unwrap();
        enc.emit().unwrap();
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

    let setup = shared_setup();
    let tmp = write_seqair_bcf(&setup, 12345, Base::A, Base::T, 30.0, 50, 0, 1, false);

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

/// A partially phased genotype round-trips through htslib unchanged.
///
/// VCF 4.5 makes the first indicator "implicitly defined as `/` if any phasing
/// indicators are `/` and `|` otherwise", so this genotype's leading bit is
/// unphased — and htslib is entitled to *omit* the indicator it can infer,
/// which it does: `0|1/1`, not `/0|1/1`. Both spell the same genotype.
///
/// So this case cannot distinguish the first-allele phase bit;
/// `bcftools_reads_seqair_phased_gt` is what pins that, because a fully phased
/// genotype is where an unset leading bit becomes visible (`/0|1`). This test
/// guards the other half: that setting the bit correctly did not disturb a
/// genotype whose leading indicator must stay unphased.
#[test]
fn bcftools_reads_a_partially_phased_gt() {
    if !has_bcftools() {
        eprintln!("skipping: bcftools not found");
        return;
    }

    let setup = shared_setup();
    let gt = Genotype {
        alleles: [Some(0), Some(1), Some(1)].into_iter().collect(),
        phased: [true, false].into_iter().collect(),
    };
    let tmp = write_seqair_bcf_with_gt(&setup, 1000, Base::C, Base::G, 99.0, 100, gt);

    let vcf_text = bcftools_view(tmp.path()).unwrap();
    let fields: Vec<&str> = vcf_text.trim().split('\t').collect();
    assert!(
        fields[9].starts_with("0|1/1"),
        "mixed phasing must survive the round trip, got: {}",
        fields[9]
    );
}

#[test]
fn bcftools_reads_seqair_phased_gt() {
    if !has_bcftools() {
        eprintln!("skipping: bcftools not found");
        return;
    }

    let setup = shared_setup();
    let tmp = write_seqair_bcf(&setup, 1000, Base::C, Base::G, 99.0, 100, 0, 1, true);

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

    let setup = shared_setup();
    let tmp = write_seqair_bcf(&setup, 500, Base::G, Base::A, 10.0, 5, 0, 0, false);

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

        let setup = shared_setup();
        let tmp = write_seqair_bcf(
            &setup, pos, ref_base, alt_base, qual, depth, gt_a0, gt_a1, phased,
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
