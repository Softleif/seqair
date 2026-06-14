//! Tests for multi-valued integer FORMAT fields (`format_ints`): a per-sample
//! integer array (e.g. `AD`, `Number=R`). Covers VCF text output, BCF
//! round-trip via bcftools, end-of-vector padding for ragged/empty samples, and
//! the cross-sample integer type scan (BCF stores one int type for the whole
//! FORMAT column, so it must fit every value across every sample).
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    reason = "test code"
)]

use seqair::vcf::record_encoder::{FormatFieldDef, FormatGt, FormatInts};
use seqair::vcf::{
    Alleles, ContigDef, Genotype, Number, OutputFormat, ValueType, VcfHeader, Writer,
};
use seqair_types::{Base, Pos1};
use std::sync::Arc;

struct Setup {
    header: Arc<VcfHeader>,
    contig: seqair::vcf::ContigId,
    gt_fmt: FormatGt,
    ad_fmt: FormatInts,
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
    let ad_fmt = builder
        .register_format(&FormatFieldDef::new(
            "AD",
            Number::Unknown,
            ValueType::Integer,
            "Per-sample integer list",
        ))
        .unwrap();
    let mut builder = builder.samples();
    for name in samples {
        builder.add_sample(*name).unwrap();
    }
    Setup { header: Arc::new(builder.build().unwrap()), contig, gt_fmt, ad_fmt }
}

// r[verify record_encoder.format_methods]
#[test]
fn vcf_text_ragged_and_empty_format_ints() {
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
    // Sample a: two values, sample b: one value, sample c: empty.
    s.ad_fmt.encode(&mut enc, &[&[12i32, 3][..], &[40][..], &[][..]]).unwrap();
    enc.emit().unwrap();
    writer.finish().unwrap();

    let text = String::from_utf8(out).unwrap();
    let data_line = text.lines().find(|l| !l.starts_with('#')).unwrap();
    assert!(data_line.ends_with("0/1:12,3\t0/1:40\t0/1:."), "samples: {data_line}");
}

// r[verify record_encoder.format_methods]
#[test]
fn bcf_ragged_and_empty_format_ints_roundtrips_through_bcftools() {
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
    s.ad_fmt.encode(&mut enc, &[&[12i32, 3][..], &[40][..], &[][..]]).unwrap();
    enc.emit().unwrap();
    writer.finish().unwrap();

    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("t.bcf");
    std::fs::write(&path, &out).unwrap();
    let q = std::process::Command::new("bcftools")
        .args(["query", "-f", "[%AD\t]\n"])
        .arg(&path)
        .output()
        .unwrap();
    assert!(q.status.success(), "bcftools stderr: {}", String::from_utf8_lossy(&q.stderr));
    let ad = String::from_utf8(q.stdout).unwrap();
    assert_eq!(ad.trim(), "12,3\t40\t.", "bcftools read AD back as: {ad}");
}

// r[verify record_encoder.format_methods]
/// The BCF integer type for the whole column must fit the largest value across
/// *all* samples — not just the first. Sample a holds a value that fits INT8;
/// sample b needs INT32. A per-sample (first-element) type choice would corrupt
/// the column, which bcftools would read back wrong.
#[test]
fn bcf_format_ints_type_scan_spans_all_samples() {
    if std::process::Command::new("bcftools").arg("--version").output().is_err() {
        eprintln!("bcftools not available, skipping");
        return;
    }
    let s = make_setup(&["a", "b"]);
    let mut out = Vec::new();
    let writer = Writer::new(&mut out, OutputFormat::Bcf);
    let mut writer = writer.write_header(&s.header).unwrap();

    let alleles = Alleles::snv(Base::A, Base::T).unwrap();
    let enc = writer.begin_record(&s.contig, Pos1::new(10).unwrap(), &alleles, None).unwrap();
    let enc = enc.filter_pass();
    let mut enc = enc.begin_samples();
    let gts: [Genotype; 2] = std::array::from_fn(|_| Genotype::unphased(0, 1));
    s.gt_fmt.encode(&mut enc, &gts).unwrap();
    // 1 fits INT8; 40000 needs INT32 → column must be INT32 for both.
    s.ad_fmt.encode(&mut enc, &[&[1i32][..], &[40000][..]]).unwrap();
    enc.emit().unwrap();
    writer.finish().unwrap();

    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("t.bcf");
    std::fs::write(&path, &out).unwrap();
    let q = std::process::Command::new("bcftools")
        .args(["query", "-f", "[%AD\t]\n"])
        .arg(&path)
        .output()
        .unwrap();
    assert!(q.status.success(), "bcftools stderr: {}", String::from_utf8_lossy(&q.stderr));
    let ad = String::from_utf8(q.stdout).unwrap();
    assert_eq!(ad.trim(), "1\t40000", "bcftools read AD back as: {ad}");
}
