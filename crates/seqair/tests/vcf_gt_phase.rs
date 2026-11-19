//! `GenotypeRef::is_phased` against the first allele's phase htslib infers
//! for a VCF 4.4+ genotype without a leading phase indicator
//! (`r[vcf_values.gt_phased]`), as bcftools stores it when it converts the
//! text to BCF; rust-htslib (htslib 1.19, which predates the inference)
//! reads the stored bit back. In an older file the same genotypes must give
//! the same answer; htslib's text parser leaves the first allele unphased
//! before 4.4 while its BCF check infers it, so it is no oracle there.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "test code"
)]

#[path = "../src/pinned.rs"]
mod pinned;

use hegel::prelude::*;
use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib::bcf::{self as hts, Read as _};
use seqair::vcf::{FormatGt, OnlySample, VcfFile};
use std::fmt::Write as _;
use std::io::Write as _;
use std::path::Path;
use std::process::Command;

const BCFTOOLS: &str = "/opt/homebrew/bin/bcftools";

fn bcftools_bin() -> Option<&'static str> {
    let bin = if Path::new(BCFTOOLS).exists() { BCFTOOLS } else { "bcftools" };
    Command::new(bin).arg("--version").output().ok().map(|_| bin)
}

/// A GT of ploidy 1–4 without a leading phase indicator: `0|1/.`, `1`, `.`.
fn genotype(tc: &TestCase) -> String {
    let ploidy = tc.draw(gs::integers::<usize>().min_value(1).max_value(4));
    let mut gt = String::new();
    for i in 0..ploidy {
        if i > 0 {
            gt.push(if tc.draw(gs::booleans()) { '|' } else { '/' });
        }
        match tc.draw(gs::integers::<u8>().max_value(2)) {
            2 => gt.push('.'),
            a => gt.push_str(&a.to_string()),
        }
    }
    gt
}

fn text(version: &str, gts: &[String]) -> String {
    let mut text = format!(
        "##fileformat=VCFv{version}\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         ##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
    );
    for (i, gt) in gts.iter().enumerate() {
        writeln!(text, "chr1\t{}\t.\tA\tC,G\t.\t.\t.\tGT\t{gt}", i + 1).unwrap();
    }
    text
}

/// `is_phased` of each record's genotype; `None` for a sample without one.
fn ours(version: &str, gts: &[String]) -> Vec<Option<bool>> {
    let file = VcfFile::from_bytes(text(version, gts).into_bytes(), None).unwrap();
    let key: FormatGt = file.format_key("GT").unwrap();
    let mut reader = file.reader::<OnlySample>().open().unwrap();
    let mut phased = Vec::new();
    while let Some(sample) = reader.next_record().unwrap() {
        phased.push(sample.get(&key).unwrap().map(|g| g.is_phased()));
    }
    phased
}

// r[verify vcf_values.gt_phased]
#[hegel::test(test_cases = pinned::cases(100))]
fn is_phased_is_htslibs_first_allele_phase(tc: TestCase) {
    let gts: Vec<String> = (0..4).map(|_| genotype(&tc)).collect();
    let ours = ours("4.4", &gts);
    assert_eq!(self::ours("4.3", &gts), ours, "the version does not matter");

    let Some(bin) = bcftools_bin() else { return };
    let mut vcf = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
    vcf.write_all(text("4.4", &gts).as_bytes()).unwrap();
    vcf.flush().unwrap();
    let bcf = tempfile::Builder::new().suffix(".bcf").tempfile().unwrap();
    let out = Command::new(bin)
        .args(["view", "-Ob", "-o"])
        .arg(bcf.path())
        .arg(vcf.path())
        .output()
        .unwrap();
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let mut htslib = hts::Reader::from_path(bcf.path()).unwrap();
    let theirs: Vec<bool> = htslib
        .records()
        .map(|rec| {
            let first = rec.unwrap().genotypes().unwrap().get(0)[0];
            matches!(first, GenotypeAllele::Phased(_) | GenotypeAllele::PhasedMissing)
        })
        .collect();
    assert_eq!(ours.len(), theirs.len());
    for ((gt, ours), theirs) in gts.iter().zip(&ours).zip(&theirs) {
        // A haploid `.` is no genotype for seqair; htslib keeps one missing allele.
        match ours {
            Some(ours) => assert_eq!(ours, theirs, "GT {gt}"),
            None => assert_eq!(gt, "."),
        }
    }
}
