//! Tests for string FORMAT fields (`format_string`): one string per sample.
//! Covers VCF text output, BCF round-trip via bcftools, the fixed-width
//! NUL-padded BCF char vector, and an empty (missing) sample reading back as `.`.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    reason = "test code"
)]

use seqair::vcf::record_encoder::{FormatFieldDef, FormatGt, FormatString};
use seqair::vcf::{
    Alleles, ContigDef, Genotype, Number, OutputFormat, ValueType, VcfHeader, Writer,
};
use seqair_types::{Base, Pos1};
use std::sync::Arc;

struct Setup {
    header: Arc<VcfHeader>,
    contig: seqair::vcf::ContigId,
    gt_fmt: FormatGt,
    lbl_fmt: FormatString,
}

fn make_setup(samples: &[&str]) -> Setup {
    let mut builder = VcfHeader::builder();
    let contig = builder.register_contig("chr1", ContigDef { length: Some(1000) }).unwrap();
    let mut builder = builder.formats();
    let gt_fmt = builder
        .register_format(&FormatFieldDef::new(
            "GT",
            Number::Count(1),
            ValueType::String,
            "Genotype",
        ))
        .unwrap();
    let lbl_fmt = builder
        .register_format(&FormatFieldDef::new(
            "LBL",
            Number::Count(1),
            ValueType::String,
            "Per-sample label",
        ))
        .unwrap();
    let mut builder = builder.samples();
    for name in samples {
        builder.add_sample(*name).unwrap();
    }
    Setup { header: Arc::new(builder.build().unwrap()), contig, gt_fmt, lbl_fmt }
}

// r[verify record_encoder.format_methods]
#[test]
fn vcf_text_format_string() {
    let s = make_setup(&["a", "b", "c"]);
    let mut out = Vec::new();
    let writer = Writer::new(&mut out, OutputFormat::Vcf);
    let mut writer = writer.write_header(&s.header).unwrap();

    let alleles = Alleles::snv(Base::A, Base::T).unwrap();
    let enc = writer.begin_record(&s.contig, Pos1::new(10).unwrap(), &alleles, None).unwrap();
    let enc = enc.filter_pass();
    let mut enc = enc.begin_samples();
    let gts: [Genotype; 3] = std::array::from_fn(|_| Genotype::unphased(0, 1));
    s.gt_fmt.encode(&mut enc, &gts).unwrap();
    // Differing lengths, plus an empty sample → `.`.
    s.lbl_fmt.encode(&mut enc, &["pass", "lo", ""]).unwrap();
    enc.emit().unwrap();
    writer.finish().unwrap();

    let text = String::from_utf8(out).unwrap();
    let data_line = text.lines().find(|l| !l.starts_with('#')).unwrap();
    assert!(data_line.ends_with("0/1:pass\t0/1:lo\t0/1:."), "samples: {data_line}");
}

// r[verify record_encoder.format_methods]
#[test]
fn bcf_format_string_roundtrips_through_bcftools() {
    if std::process::Command::new("bcftools").arg("--version").output().is_err() {
        eprintln!("bcftools not available, skipping");
        return;
    }
    let s = make_setup(&["a", "b", "c"]);
    let mut out = Vec::new();
    let writer = Writer::new(&mut out, OutputFormat::Bcf);
    let mut writer = writer.write_header(&s.header).unwrap();

    let alleles = Alleles::snv(Base::A, Base::T).unwrap();
    let enc = writer.begin_record(&s.contig, Pos1::new(10).unwrap(), &alleles, None).unwrap();
    let enc = enc.filter_pass();
    let mut enc = enc.begin_samples();
    let gts: [Genotype; 3] = std::array::from_fn(|_| Genotype::unphased(0, 1));
    s.gt_fmt.encode(&mut enc, &gts).unwrap();
    s.lbl_fmt.encode(&mut enc, &["pass", "lo", ""]).unwrap();
    enc.emit().unwrap();
    writer.finish().unwrap();

    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("t.bcf");
    std::fs::write(&path, &out).unwrap();
    let q = std::process::Command::new("bcftools")
        .args(["query", "-f", "[%LBL\t]\n"])
        .arg(&path)
        .output()
        .unwrap();
    assert!(q.status.success(), "bcftools stderr: {}", String::from_utf8_lossy(&q.stderr));
    let lbl = String::from_utf8(q.stdout).unwrap();
    // Longest string is "pass" (width 4); "lo" NUL-padded to 4, empty sample → `.`.
    assert_eq!(lbl.trim(), "pass\tlo\t.", "bcftools read LBL back as: {lbl}");
}
