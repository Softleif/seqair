//! A VCF text record's end against htslib's (`r[vcf_values.rlen]`): random
//! records whose span depends on REF, INFO END, SVLEN of symbolic SV
//! alleles and FORMAT LEN of gVCF blocks, converted to BCF by bcftools,
//! whose stored `rlen` (read by seqair's BCF reader) is the oracle for the
//! end seqair computes from the text.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "test code"
)]

#[path = "../src/pinned.rs"]
mod pinned;

use hegel::prelude::*;
use seqair::vcf::{Locus, VcfFile, VcfPos, VcfReadError};
use std::io::Write;
use std::process::Command;

const BCFTOOLS: &str = "/opt/homebrew/bin/bcftools";

fn bcftools() -> Option<&'static str> {
    let bin = if std::path::Path::new(BCFTOOLS).exists() { BCFTOOLS } else { "bcftools" };
    Command::new(bin).arg("--version").output().ok().map(|_| bin)
}

const HEADER: &str = "##fileformat=VCFv4.5\n\
##INFO=<ID=END,Number=1,Type=Integer,Description=\"end\">\n\
##INFO=<ID=SVLEN,Number=A,Type=Integer,Description=\"svlen\">\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"dp\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"gt\">\n\
##FORMAT=<ID=LEN,Number=1,Type=Integer,Description=\"len\">\n\
##contig=<ID=chr1,length=100000000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n";

const ALTS: [&str; 12] = [
    "A",
    "CT",
    "<DEL>",
    "<DUP>",
    "<INV>",
    "<CNV>",
    "<INS>",
    "<DEL:ME>",
    "<*>",
    "<NON_REF>",
    "*",
    "<DUP:TANDEM>",
];

/// An integer as a VCF value: small (so it lands near POS), large, negative
/// or missing.
fn value(tc: &TestCase, pos: u32) -> String {
    match tc.draw(gs::integers::<u8>().max_value(5)) {
        0 => ".".to_owned(),
        1 => tc.draw(gs::integers::<i32>().min_value(-500).max_value(500)).to_string(),
        2 => tc.draw(gs::integers::<i32>()).to_string(),
        _ => (i64::from(pos)
            + i64::from(tc.draw(gs::integers::<i32>().min_value(-20).max_value(2000))))
        .to_string(),
    }
}

/// One data line at `pos`.
fn line(tc: &TestCase, pos: u32) -> String {
    let ref_len = tc.draw(gs::integers::<usize>().min_value(1).max_value(6));
    let reference: String = "ACGTAC"[..ref_len].to_owned();
    // Up to 10 ALTs, so SVLEN's eighth element (htslib reads ALT 0's bit
    // for it) is reached.
    let n_alts = tc.draw(gs::integers::<usize>().max_value(10));
    let alts: Vec<&str> = (0..n_alts).map(|_| tc.draw(gs::sampled_from(ALTS.to_vec()))).collect();
    let alt = if alts.is_empty() { ".".to_owned() } else { alts.join(",") };
    let mut info = Vec::new();
    if tc.draw(gs::booleans()) {
        info.push(format!("DP={}", tc.draw(gs::integers::<u8>())));
    }
    if tc.draw(gs::booleans()) {
        let end = if tc.draw(gs::integers::<u8>().max_value(9)) == 0 {
            format!("{},{}", value(tc, pos), value(tc, pos))
        } else {
            value(tc, pos)
        };
        info.push(format!("END={end}"));
    }
    if tc.draw(gs::booleans()) {
        let n = tc.draw(gs::integers::<usize>().min_value(1).max_value(11));
        let values: Vec<String> = (0..n).map(|_| value(tc, 0)).collect();
        info.push(format!("SVLEN={}", values.join(",")));
    }
    let info = if info.is_empty() { ".".to_owned() } else { info.join(";") };
    let (format, samples) = if tc.draw(gs::booleans()) {
        ("GT:LEN", [format!("0/1:{}", value(tc, 0)), format!("0/0:{}", value(tc, 0))])
    } else {
        ("GT", ["0/1".to_owned(), "1/1".to_owned()])
    };
    format!(
        "chr1\t{pos}\t.\t{reference}\t{alt}\t.\t.\t{info}\t{format}\t{}\t{}",
        samples[0], samples[1]
    )
}

/// Each record's start and end, or the kind of error reading its locus failed with.
fn ends(file: &VcfFile) -> Vec<Result<(VcfPos, VcfPos), String>> {
    let mut reader = file.reader::<Locus>().open().unwrap();
    let mut out = Vec::new();
    loop {
        match reader.read() {
            Ok(Some(l)) => out.push(Ok((l.start(), l.end()))),
            Ok(None) => return out,
            // The kind only: BCF stores the span in 32 bits, the text's is wider.
            Err(VcfReadError::Record { source, .. }) => {
                out.push(Err(format!("{source:?}")
                    .split(' ')
                    .next()
                    .unwrap_or_default()
                    .to_owned()));
            }
            Err(e) => panic!("{e:?}"),
        }
    }
}

// r[verify vcf_values.rlen+3]
// r[verify vcf_reader.oracles]
/// Every generated record ends where the BCF bcftools writes from the same
/// text says it ends, whatever the version the file declares.
#[hegel::test(test_cases = pinned::cases(100))]
fn text_ends_are_the_rlen_bcftools_writes(tc: TestCase) {
    let Some(bin) = bcftools() else {
        eprintln!("bcftools not found, skipping");
        return;
    };
    let n = tc.draw(gs::integers::<usize>().min_value(1).max_value(20));
    let mut pos = 0u32;
    let version = tc.draw(gs::sampled_from(vec!["4.2", "4.3", "4.4", "4.5"]));
    let mut text = HEADER.replacen("4.5", version, 1);
    for _ in 0..n {
        pos += tc.draw(gs::integers::<u32>().min_value(1).max_value(5000));
        text.push_str(&line(&tc, pos));
        text.push('\n');
    }
    let mut vcf = tempfile::NamedTempFile::new().unwrap();
    vcf.write_all(text.as_bytes()).unwrap();
    let bcf = tempfile::NamedTempFile::new().unwrap();
    let out = Command::new(bin)
        .args(["view", "-Ob", "-o"])
        .arg(bcf.path())
        .arg(vcf.path())
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "bcftools: {}",
        String::from_utf8(out.stderr.clone()).unwrap_or_default()
    );
    let theirs = ends(&VcfFile::open(bcf.path()).unwrap());
    let ours = ends(&VcfFile::from_bytes(text.clone().into_bytes(), None).unwrap());
    assert_eq!(ours.len(), theirs.len());
    for ((o, t), l) in ours.iter().zip(&theirs).zip(text.lines().skip(HEADER.lines().count())) {
        assert_eq!(o, t, "{l}");
    }
}

// r[verify vcf_values.rlen+3]
/// htslib reads SVLEN's eighth element against ALT 0's bit (REF, never
/// set), so an eighth-ALT `<DEL>` does not extend the record, a seventh
/// does. Ends as `bcftools query -f '%END'` prints them.
#[test]
fn the_eighth_svlen_is_read_as_htslib_reads_it() {
    let text = "##fileformat=VCFv4.5\n\
        ##INFO=<ID=SVLEN,Number=A,Type=Integer,Description=\"s\">\n##contig=<ID=chr1>\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
        chr1\t100\t.\tA\tC,C,C,C,C,C,C,<DEL>\t.\t.\tSVLEN=.,.,.,.,.,.,.,500\n\
        chr1\t200\t.\tA\tC,C,C,C,C,C,<DEL>,C\t.\t.\tSVLEN=.,.,.,.,.,.,500,.\n";
    let ends: Vec<_> = ends(&VcfFile::from_bytes(text.as_bytes().to_vec(), None).unwrap())
        .into_iter()
        .map(|r| r.unwrap().1)
        .collect();
    let at = |p| VcfPos::At(seqair_types::Pos1::new(p).unwrap());
    assert_eq!(ends, [at(100), at(700)]);
}
