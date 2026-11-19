//! Record IDs (`r[record_encoder.record_id]`): written `;`-separated in VCF
//! text and BCF alike, as bcftools reads them back, and refused where VCF's
//! ID grammar refuses them.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::format_collect,
    reason = "test code"
)]

#[path = "../src/pinned.rs"]
mod pinned;

use hegel::prelude::*;
use seqair::vcf::record_encoder::{FormatFieldDef, Gt, InfoFieldDef, Scalar};
use seqair::vcf::{
    Alleles, ContigDef, ContigId, FormatGt, Genotype, InfoInt, OutputFormat, RecordIdProblem,
    VcfError, VcfHeader, Writer,
};
use seqair_types::{Base, Pos1};
use std::process::Command;

fn has_bcftools() -> bool {
    Command::new("bcftools").arg("--version").output().is_ok()
}

struct Setup {
    header: VcfHeader,
    chr1: ContigId,
    dp: InfoInt,
    gt: FormatGt,
}

fn setup() -> Setup {
    let mut builder = VcfHeader::builder();
    let chr1 = builder.register_contig("chr1", ContigDef::new(Some(1_000))).unwrap();
    let mut builder = builder.infos();
    let dp = builder.register_info(&InfoFieldDef::<Scalar<i32>>::new("DP", "Depth")).unwrap();
    let mut builder = builder.formats();
    let gt = builder.register_format(&FormatFieldDef::<Gt>::new("GT", "Genotype")).unwrap();
    let mut builder = builder.samples();
    builder.add_sample("S1").unwrap();
    Setup { header: builder.build().unwrap(), chr1, dp, gt }
}

/// Each record's ID lists, every list but the last set and then replaced.
fn write(setup: &Setup, format: OutputFormat, records: &[Vec<Vec<String>>]) -> Vec<u8> {
    let mut out = Vec::new();
    let mut writer = Writer::new(&mut out, format).write_header(&setup.header).unwrap();
    for (i, lists) in records.iter().enumerate() {
        let alleles = Alleles::complex("ACG", ["A"]);
        let pos = Pos1::new(u32::try_from(i).unwrap() + 1).unwrap();
        let mut enc = writer.begin_record(&setup.chr1, pos, &alleles, Some(5.0)).unwrap();
        for ids in lists {
            enc = enc.ids(ids.iter().map(String::as_str)).unwrap();
        }
        let mut enc = enc.filter_pass();
        setup.dp.encode(&mut enc, 7).unwrap();
        let mut enc = enc.begin_samples();
        setup.gt.encode(&mut enc, &[Genotype::unphased(0, 1)]).unwrap();
        enc.emit().unwrap();
    }
    drop(writer.finish().unwrap());
    out
}

fn query(bytes: &[u8], suffix: &str) -> String {
    let mut file = tempfile::Builder::new().suffix(suffix).tempfile().unwrap();
    std::io::Write::write_all(&mut file, bytes).unwrap();
    let out = Command::new("bcftools")
        .args(["query", "-f", "%ID %REF %ALT %QUAL %DP [%GT]\n"])
        .arg(file.path())
        .output()
        .unwrap();
    assert!(out.status.success(), "bcftools: {}", String::from_utf8_lossy(&out.stderr));
    String::from_utf8(out.stdout).unwrap()
}

/// Up to four distinct identifiers of the characters real IDs use, `.` excluded.
#[hegel::composite]
fn ids(tc: &TestCase) -> Vec<String> {
    let n = tc.draw(gs::integers::<usize>().max_value(4));
    let mut out: Vec<String> = Vec::new();
    for i in 0..n {
        let id: String = tc.draw(gs::from_regex("[A-Za-z0-9_:,|.é-]{0,40}").fullmatch(true));
        // Distinct by construction, and never a lone `.`.
        out.push(format!("{id}{i}"));
    }
    out
}

// r[verify record_encoder.record_id]
// r[verify record_encoder.typestate_transitions+3]
/// bcftools reads the ID column seqair writes, the last list set winning, and
/// the rest of the record unmoved by the replacement, from VCF text and BCF.
#[hegel::test(test_cases = pinned::cases(40))]
fn bcftools_reads_the_ids_as_written(tc: TestCase) {
    if !has_bcftools() {
        return;
    }
    let n = tc.draw(gs::integers::<usize>().min_value(1).max_value(4));
    let records: Vec<Vec<Vec<String>>> = (0..n)
        .map(|_| {
            let lists = tc.draw(gs::integers::<usize>().max_value(3));
            (0..lists).map(|_| tc.draw(ids())).collect()
        })
        .collect();
    let setup = setup();
    let expected: String = records
        .iter()
        .map(|lists| {
            let id = lists
                .last()
                .filter(|l| !l.is_empty())
                .map_or_else(|| ".".to_owned(), |l| l.join(";"));
            format!("{id} ACG A 5 7 0/1\n")
        })
        .collect();
    assert_eq!(query(&write(&setup, OutputFormat::Vcf, &records), ".vcf"), expected, "VCF");
    assert_eq!(query(&write(&setup, OutputFormat::Bcf, &records), ".bcf"), expected, "BCF");
}

// r[verify record_encoder.record_id]
/// Every way an identifier breaks VCF's grammar is refused, naming it.
#[test]
fn ids_vcf_forbids_are_refused() {
    let setup = setup();
    let alleles = Alleles::snv(Base::A, Base::C).unwrap();
    for format in [OutputFormat::Vcf, OutputFormat::Bcf] {
        let mut writer = Writer::new(Vec::new(), format).write_header(&setup.header).unwrap();
        let cases: [(&[&str], RecordIdProblem); 5] = [
            (&["rs1", ""], RecordIdProblem::Empty),
            (&["."], RecordIdProblem::Missing),
            (&["rs 1"], RecordIdProblem::Whitespace),
            (&["rs1\t"], RecordIdProblem::Whitespace),
            (&["rs1;rs2"], RecordIdProblem::Semicolon),
        ];
        for (ids, problem) in cases {
            let enc =
                writer.begin_record(&setup.chr1, Pos1::new(1).unwrap(), &alleles, None).unwrap();
            let err = enc.ids(ids.iter().copied()).err();
            assert!(
                matches!(&err, Some(VcfError::InvalidRecordId { problem: p, .. }) if *p == problem),
                "{ids:?}: {err:?}"
            );
        }
        let enc = writer.begin_record(&setup.chr1, Pos1::new(1).unwrap(), &alleles, None).unwrap();
        let err = enc.ids(["rs1", "rs2", "rs1"]).err();
        assert!(matches!(&err, Some(VcfError::DuplicateRecordId { id }) if id == "rs1"), "{err:?}");
        // The writer carries on: a refused record was never written.
        let enc = writer.begin_record(&setup.chr1, Pos1::new(2).unwrap(), &alleles, None).unwrap();
        enc.ids(["rs9"]).unwrap().filter_pass().emit().unwrap();
        drop(writer.finish().unwrap());
    }
}
