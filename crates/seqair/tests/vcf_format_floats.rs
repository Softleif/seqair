//! Tests for multi-valued float FORMAT fields (`format_floats`): a per-sample
//! float array (e.g. `Number=A`/`Number=M` FORMAT fields). Covers VCF text
//! output and BCF round-trip via bcftools.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    reason = "test code"
)]

use seqair::vcf::record_encoder::{FormatFieldDef, FormatFloats, FormatGt, InfoFieldDef, InfoInt};
use seqair::vcf::{
    Alleles, ContigDef, Genotype, Number, OutputFormat, ValueType, VcfHeader, Writer,
};
use seqair_types::{Base, Pos1};
use std::sync::Arc;

struct Setup {
    header: Arc<VcfHeader>,
    contig: seqair::vcf::ContigId,
    dp_info: InfoInt,
    gt_fmt: FormatGt,
    ml_fmt: FormatFloats,
}

fn make_setup() -> Setup {
    let mut builder = VcfHeader::builder();
    let contig = builder.register_contig("chr1", ContigDef { length: Some(1000) }).unwrap();
    let mut builder = builder.infos();
    let dp_info = builder
        .register_info(&InfoFieldDef::new("DP", Number::Count(1), ValueType::Integer, "Depth"))
        .unwrap();
    let mut builder = builder.formats();
    let gt_fmt = builder
        .register_format(&FormatFieldDef::new(
            "GT",
            Number::Count(1),
            ValueType::String,
            "Genotype",
        ))
        .unwrap();
    let ml_fmt = builder
        .register_format(&FormatFieldDef::new(
            "ML",
            Number::AlternateBases,
            ValueType::Float,
            "ML prediction per alt",
        ))
        .unwrap();
    let mut builder = builder.samples();
    builder.add_sample("s1").unwrap();
    Setup { header: Arc::new(builder.build().unwrap()), contig, dp_info, gt_fmt, ml_fmt }
}

#[test]
fn vcf_text_multi_value_format_float() {
    let s = make_setup();
    let mut out = Vec::new();
    let writer = Writer::new(&mut out, OutputFormat::Vcf);
    let mut writer = writer.write_header(&s.header).unwrap();

    // Two ALT alleles → ML has two values for the single sample.
    let alleles = Alleles::snv_multi(Base::A, &[Base::T, Base::C]).unwrap();
    let enc = writer.begin_record(&s.contig, Pos1::new(10).unwrap(), &alleles, Some(30.0)).unwrap();
    let mut enc = enc.filter_pass();
    s.dp_info.encode(&mut enc, 50);
    let mut enc = enc.begin_samples();
    s.gt_fmt.encode(&mut enc, &[Genotype::unphased(0, 1)]).unwrap();
    s.ml_fmt.encode(&mut enc, &[&[0.25f32, 0.75]]).unwrap();
    enc.emit().unwrap();
    writer.finish().unwrap();

    let text = String::from_utf8(out).unwrap();
    let data_line = text.lines().find(|l| !l.starts_with('#')).unwrap();
    // GT:ML  with ML = 0.25,0.75
    assert!(data_line.contains("\tGT:ML\t"), "FORMAT keys: {data_line}");
    assert!(data_line.ends_with("0/1:0.25,0.75"), "sample column: {data_line}");
}

#[test]
fn bcf_multi_value_format_float_roundtrips_through_bcftools() {
    if std::process::Command::new("bcftools").arg("--version").output().is_err() {
        eprintln!("bcftools not available, skipping");
        return;
    }
    let s = make_setup();
    let mut out = Vec::new();
    let writer = Writer::new(&mut out, OutputFormat::Bcf);
    let mut writer = writer.write_header(&s.header).unwrap();

    let alleles = Alleles::snv_multi(Base::A, &[Base::T, Base::C]).unwrap();
    let enc = writer.begin_record(&s.contig, Pos1::new(10).unwrap(), &alleles, Some(30.0)).unwrap();
    let mut enc = enc.filter_pass();
    s.dp_info.encode(&mut enc, 50);
    let mut enc = enc.begin_samples();
    s.gt_fmt.encode(&mut enc, &[Genotype::unphased(0, 1)]).unwrap();
    s.ml_fmt.encode(&mut enc, &[&[0.25f32, 0.75]]).unwrap();
    enc.emit().unwrap();
    writer.finish().unwrap();

    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("t.bcf");
    std::fs::write(&path, &out).unwrap();
    let q = std::process::Command::new("bcftools")
        .args(["query", "-f", "[%ML]\n"])
        .arg(&path)
        .output()
        .unwrap();
    assert!(q.status.success(), "bcftools stderr: {}", String::from_utf8_lossy(&q.stderr));
    let ml = String::from_utf8(q.stdout).unwrap();
    assert_eq!(ml.trim(), "0.25,0.75", "bcftools read ML back as: {ml}");
}
