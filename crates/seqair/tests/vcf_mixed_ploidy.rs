//! Samples of differing ploidy in one record (`r[record_encoder.mixed_ploidy]`),
//! read back by bcftools from VCF text and BCF alike: haploid next to diploid
//! next to triploid, with a `Number=G` field whose count follows each sample's
//! own ploidy.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::format_push_string,
    clippy::format_collect,
    reason = "test code"
)]

#[path = "../src/pinned.rs"]
mod pinned;

use hegel::prelude::*;
use seqair::vcf::record_encoder::{Arr, FormatFieldDef, Gt};
use seqair::vcf::{
    Alleles, ContigDef, ContigId, FormatGt, FormatInts, Genotype, Number, OutputFormat, VcfError,
    VcfHeader, Writer,
};
use seqair_types::{Base, Pos1, SmallVec};
use std::process::Command;

fn has_bcftools() -> bool {
    Command::new("bcftools").arg("--version").output().is_ok()
}

struct Setup {
    header: VcfHeader,
    chr1: ContigId,
    gt: FormatGt,
    pl: FormatInts,
}

fn setup(n_samples: usize) -> Setup {
    let mut builder = VcfHeader::builder();
    let chr1 = builder.register_contig("chr1", ContigDef::new(Some(1_000_000))).unwrap();
    let mut builder = builder.formats();
    let gt = builder.register_format(&FormatFieldDef::<Gt>::new("GT", "Genotype")).unwrap();
    let pl = builder
        .register_format(&FormatFieldDef::<Arr<i32>>::with_number(
            "PL",
            Number::Genotypes,
            "Likelihoods",
        ))
        .unwrap();
    let mut builder = builder.samples();
    for i in 0..n_samples {
        builder.add_sample(format!("S{i}")).unwrap();
    }
    Setup { header: builder.build().unwrap(), chr1, gt, pl }
}

/// A genotype of 0 to 3 alleles over `n_allele` alleles, any missing, each
/// separator phased or not.
#[hegel::composite]
fn genotype(tc: &TestCase, n_allele: u16) -> Genotype {
    let ploidy = tc.draw(gs::integers::<usize>().max_value(3));
    let alleles: SmallVec<Option<u16>, 2> = (0..ploidy)
        .map(|_| tc.draw(gs::optional(gs::integers::<u16>().max_value(n_allele - 1))))
        .collect();
    let phased = (1..ploidy).map(|_| tc.draw(gs::booleans())).collect();
    Genotype { alleles, phased }
}

/// The VCF text of a genotype, by the spec's grammar: no padding, `.` for
/// none at all.
fn gt_text(g: &Genotype) -> String {
    if g.alleles.is_empty() {
        return ".".to_owned();
    }
    let mut out = String::new();
    for (i, a) in g.alleles.iter().enumerate() {
        if i > 0 {
            out.push(if g.phased[i - 1] { '|' } else { '/' });
        }
        match a {
            Some(a) => out.push_str(&a.to_string()),
            None => out.push('.'),
        }
    }
    out
}

/// `binomial(n_allele + ploidy - 1, ploidy)`, by counting the multisets.
fn genotype_count(n_allele: usize, ploidy: usize) -> usize {
    fn count(n: usize, k: usize, from: usize) -> usize {
        if k == 0 { 1 } else { (from..n).map(|a| count(n, k - 1, a)).sum() }
    }
    count(n_allele, ploidy, 0)
}

fn write(setup: &Setup, format: OutputFormat, records: &[Vec<Genotype>]) -> Vec<u8> {
    let mut out = Vec::new();
    let mut writer = Writer::new(&mut out, format).write_header(&setup.header).unwrap();
    let alleles = Alleles::snv_multi(Base::A, &[Base::C, Base::G]).unwrap();
    for (i, gts) in records.iter().enumerate() {
        let pos = Pos1::new(u32::try_from(i).unwrap() * 10 + 1).unwrap();
        let mut enc = writer
            .begin_record(&setup.chr1, pos, &alleles, None)
            .unwrap()
            .filter_pass()
            .begin_samples();
        setup.gt.encode(&mut enc, gts).unwrap();
        let pls: Vec<Vec<i32>> = gts
            .iter()
            .map(|g| match g.alleles.len() {
                0 => Vec::new(),
                p => (0..genotype_count(3, p)).map(|v| i32::try_from(v).unwrap()).collect(),
            })
            .collect();
        let pls: Vec<&[i32]> = pls.iter().map(Vec::as_slice).collect();
        setup.pl.encode(&mut enc, &pls).unwrap();
        enc.emit().unwrap();
    }
    drop(writer.finish().unwrap());
    out
}

fn query(bytes: &[u8], suffix: &str) -> String {
    let mut file = tempfile::Builder::new().suffix(suffix).tempfile().unwrap();
    std::io::Write::write_all(&mut file, bytes).unwrap();
    let out = Command::new("bcftools")
        .args(["query", "-f", "[%GT:%PL\t]\n"])
        .arg(file.path())
        .output()
        .unwrap();
    assert!(out.status.success(), "bcftools: {}", String::from_utf8_lossy(&out.stderr));
    String::from_utf8(out.stdout).unwrap()
}

// r[verify record_encoder.mixed_ploidy]
// r[verify vcf_writer.genotype_serialization+2]
// r[verify bcf_writer.gt_encoding]
// r[verify record_encoder.count_check+3]
/// bcftools reads every sample at its own ploidy from both formats, with the
/// `Number=G` values its ploidy calls for.
#[hegel::test(test_cases = pinned::cases(30))]
fn bcftools_reads_each_sample_at_its_own_ploidy(tc: TestCase) {
    if !has_bcftools() {
        return;
    }
    let n_samples = tc.draw(gs::integers::<usize>().min_value(1).max_value(5));
    let n_records = tc.draw(gs::integers::<usize>().min_value(1).max_value(4));
    let records: Vec<Vec<Genotype>> =
        (0..n_records).map(|_| (0..n_samples).map(|_| tc.draw(genotype(3))).collect()).collect();
    let setup = setup(n_samples);

    let expected: String = records
        .iter()
        .map(|gts| {
            let mut line = String::new();
            for g in gts {
                let pl = match g.alleles.len() {
                    0 => ".".to_owned(),
                    p => (0..genotype_count(3, p))
                        .map(|v| v.to_string())
                        .collect::<Vec<_>>()
                        .join(","),
                };
                line.push_str(&format!("{}:{pl}\t", gt_text(g)));
            }
            line + "\n"
        })
        .collect();
    assert_eq!(query(&write(&setup, OutputFormat::Vcf, &records), ".vcf"), expected, "VCF");
    assert_eq!(query(&write(&setup, OutputFormat::Bcf, &records), ".bcf"), expected, "BCF");
}

// r[verify record_encoder.count_check+3]
/// A `Number=G` count is each sample's own: two values for a haploid sample
/// at a biallelic site, three for a diploid one, and a diploid count given
/// to the haploid sample names it.
#[test]
fn a_genotype_count_follows_the_samples_ploidy() {
    let setup = setup(2);
    let mut writer =
        Writer::new(Vec::new(), OutputFormat::Bcf).write_header(&setup.header).unwrap();
    let alleles = Alleles::snv(Base::A, Base::C).unwrap();
    let mut enc = writer
        .begin_record(&setup.chr1, Pos1::new(5).unwrap(), &alleles, None)
        .unwrap()
        .filter_pass()
        .begin_samples();
    setup.gt.encode(&mut enc, &[Genotype::haploid(1), Genotype::unphased(0, 1)]).unwrap();
    let err = setup.pl.encode(&mut enc, &[&[0, 1, 2], &[0, 1, 2]]).unwrap_err();
    assert!(matches!(err, VcfError::CountMismatch { expected: 2, found: 3, .. }), "{err:?}");
    setup.pl.encode(&mut enc, &[&[0, 1], &[0, 1, 2]]).unwrap();
    enc.emit().unwrap();
}
