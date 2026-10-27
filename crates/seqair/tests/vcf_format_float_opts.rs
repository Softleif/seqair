//! A float FORMAT array with per-element missing (`FormatKey<OptArr<f32>>`): `None` is the BCF
//! float missing value and `.` in VCF text.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    reason = "test code"
)]

use seqair::vcf::record_encoder::{FormatFieldDef, FormatFloatOpts, OptArr};
use seqair::vcf::{
    Alleles, ContigDef, Number, OutputFormat, ValueType, VcfEncodeError, VcfError, VcfHeader,
    Writer,
};
use seqair_types::{Base, Pos1};
use std::sync::Arc;

fn write(format: OutputFormat, values: &[Option<f32>]) -> (Vec<u8>, Result<(), VcfError>) {
    let mut builder = VcfHeader::builder();
    let contig = builder.register_contig("chr1", ContigDef { length: Some(1000) }).unwrap();
    let mut builder = builder.infos().formats();
    let ml: FormatFloatOpts = builder
        .register_format(&FormatFieldDef::<OptArr<f32>>::new(
            "ML",
            Number::AlternateBases,
            ValueType::Float,
            "ML per alt",
        ))
        .unwrap();
    let mut builder = builder.samples();
    builder.add_sample("s1").unwrap();
    let header = Arc::new(builder.build().unwrap());

    let mut out = Vec::new();
    let mut writer = Writer::new(&mut out, format).write_header(&header).unwrap();
    let alleles = Alleles::snv_multi(Base::A, &[Base::T, Base::C]).unwrap();
    let enc = writer.begin_record(&contig, Pos1::new(10).unwrap(), &alleles, None).unwrap();
    let mut enc = enc.filter_pass().begin_samples();
    let result = ml.encode_single_sample(&mut enc, values);
    if result.is_ok() {
        enc.emit().unwrap();
    }
    writer.finish().unwrap();
    (out, result)
}

// r[verify record_encoder.format_methods]
// r[verify vcf_writer.missing_dot]
#[test]
fn vcf_text_writes_dot_for_none() {
    let (out, result) = write(OutputFormat::Vcf, &[Some(0.5), None]);
    result.unwrap();
    let text = String::from_utf8(out).unwrap();
    let line = text.lines().find(|l| !l.starts_with('#')).unwrap();
    assert!(line.ends_with("\tML\t0.5,."), "{line}");
}

// r[verify bcf_writer.missing_sentinels]
#[test]
fn bcf_writes_float_missing_for_none() {
    if std::process::Command::new("bcftools").arg("--version").output().is_err() {
        eprintln!("bcftools not available, skipping");
        return;
    }
    let (out, result) = write(OutputFormat::Bcf, &[Some(0.5), None]);
    result.unwrap();
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("t.bcf");
    std::fs::write(&path, &out).unwrap();
    let q = std::process::Command::new("bcftools")
        .args(["query", "-f", "[%ML]\n"])
        .arg(&path)
        .output()
        .unwrap();
    assert!(q.status.success(), "bcftools stderr: {}", String::from_utf8_lossy(&q.stderr));
    assert_eq!(String::from_utf8(q.stdout).unwrap().trim(), "0.5,.");
}

// r[verify record_encoder.reserved_rejection]
#[test]
fn reserved_float_still_rejected_next_to_none() {
    let reserved = f32::from_bits(0x7F80_0001);
    let (_, result) = write(OutputFormat::Bcf, &[None, Some(reserved)]);
    assert!(matches!(result, Err(VcfError::Encode(VcfEncodeError::ReservedFloatValue { .. }))));
}
