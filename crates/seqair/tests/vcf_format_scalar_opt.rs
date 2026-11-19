//! `FormatKey<Scalar<T>>::encode_opt` (`r[record_encoder.format_scalar_opt]`):
//! one value per sample, any of them missing, read back by bcftools from VCF
//! text and BCF, and in BCF stored as the missing value of the width the
//! column takes.
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
use seqair::vcf::record_encoder::{FormatFieldDef, Scalar};
use seqair::vcf::{
    Alleles, ContigDef, ContigId, FormatFloat, FormatInt, OutputFormat, VcfError, VcfHeader, Writer,
};
use seqair_types::{Base, Pos1};
use std::process::Command;

fn has_bcftools() -> bool {
    Command::new("bcftools").arg("--version").output().is_ok()
}

struct Setup {
    header: VcfHeader,
    chr1: ContigId,
    dp: FormatInt,
    vf: FormatFloat,
}

fn setup(n_samples: usize) -> Setup {
    let mut builder = VcfHeader::builder();
    let chr1 = builder.register_contig("chr1", ContigDef::new(Some(1_000))).unwrap();
    let mut builder = builder.formats();
    let dp = builder.register_format(&FormatFieldDef::<Scalar<i32>>::new("DP", "Depth")).unwrap();
    let vf =
        builder.register_format(&FormatFieldDef::<Scalar<f32>>::new("VF", "Fraction")).unwrap();
    let mut builder = builder.samples();
    for i in 0..n_samples {
        builder.add_sample(format!("S{i}")).unwrap();
    }
    Setup { header: builder.build().unwrap(), chr1, dp, vf }
}

fn write(setup: &Setup, format: OutputFormat, dp: &[Option<i32>], vf: &[Option<f32>]) -> Vec<u8> {
    let mut out = Vec::new();
    let mut writer = Writer::new(&mut out, format).write_header(&setup.header).unwrap();
    let alleles = Alleles::snv(Base::A, Base::C).unwrap();
    let mut enc = writer
        .begin_record(&setup.chr1, Pos1::new(10).unwrap(), &alleles, None)
        .unwrap()
        .filter_pass()
        .begin_samples();
    setup.dp.encode_opt(&mut enc, dp).unwrap();
    setup.vf.encode_opt(&mut enc, vf).unwrap();
    enc.emit().unwrap();
    drop(writer.finish().unwrap());
    out
}

fn query(bytes: &[u8], suffix: &str) -> String {
    let mut file = tempfile::Builder::new().suffix(suffix).tempfile().unwrap();
    std::io::Write::write_all(&mut file, bytes).unwrap();
    let out = Command::new("bcftools")
        .args(["query", "-f", "[%DP:%VF\t]\n"])
        .arg(file.path())
        .output()
        .unwrap();
    assert!(out.status.success(), "bcftools: {}", String::from_utf8_lossy(&out.stderr));
    String::from_utf8(out.stdout).unwrap()
}

/// Integers outside the int32 reserved band, at every width, or none.
#[hegel::composite]
fn int(tc: &TestCase) -> Option<i32> {
    tc.draw(gs::optional(hegel::one_of!(
        gs::integers::<i32>().min_value(-120).max_value(127),
        gs::integers::<i32>().min_value(-32_760).max_value(32_767),
        gs::integers::<i32>().min_value(i32::MIN + 8),
    )))
}

// r[verify record_encoder.format_scalar_opt]
// r[verify record_encoder.key_encode+2]
/// bcftools reads every sample's value, and `.` for every missing one, from
/// VCF text and BCF alike. (Floats are drawn from values `%g` prints exactly.)
#[hegel::test(test_cases = pinned::cases(40))]
fn bcftools_reads_missing_samples_as_missing(tc: TestCase) {
    if !has_bcftools() {
        return;
    }
    let n = tc.draw(gs::integers::<usize>().min_value(1).max_value(6));
    let dp: Vec<Option<i32>> = (0..n).map(|_| tc.draw(int())).collect();
    let vf: Vec<Option<f32>> = (0..n)
        .map(|_| tc.draw(gs::optional(gs::integers::<i16>())).map(|v| f32::from(v) / 4.0))
        .collect();
    let setup = setup(n);
    let text = |v: Option<String>| v.unwrap_or_else(|| ".".to_owned());
    let expected: String = dp
        .iter()
        .zip(&vf)
        .map(|(d, f)| {
            format!("{}:{}\t", text(d.map(|d| d.to_string())), text(f.map(|f| f.to_string())))
        })
        .collect::<String>()
        + "\n";
    assert_eq!(query(&write(&setup, OutputFormat::Vcf, &dp, &vf), ".vcf"), expected, "VCF");
    assert_eq!(query(&write(&setup, OutputFormat::Bcf, &dp, &vf), ".bcf"), expected, "BCF");
}

/// The bytes of the record's FORMAT block, its last bytes in the file.
fn bcf_tail(bytes: Vec<u8>, n: usize) -> Vec<u8> {
    let mut data = Vec::new();
    seqair::bam::bgzf::BgzfReader::from_reader(std::io::Cursor::new(bytes))
        .read_to_end(&mut data)
        .unwrap();
    data[data.len() - n..].to_vec()
}

// r[verify record_encoder.format_scalar_opt]
// r[verify bcf_writer.missing_sentinels]
/// A missing sample is the missing value of the width the column needs for
/// its concrete values — `0x8000` next to 300, which needs int16 — and int8's
/// `0x80` when no sample has a value; a missing float is `0x7F800001`.
#[test]
fn a_missing_sample_takes_the_columns_width() {
    let setup = setup(3);
    let dp_idx = u8::try_from(setup.header.string_map().get("DP").unwrap()).unwrap();
    let vf_idx = u8::try_from(setup.header.string_map().get("VF").unwrap()).unwrap();
    let bytes =
        write(&setup, OutputFormat::Bcf, &[Some(5), None, Some(300)], &[None, Some(0.5), None]);
    #[rustfmt::skip]
    let expected = [
        0x11, dp_idx, 0x12, 0x05, 0x00, 0x00, 0x80, 0x2c, 0x01,
        0x11, vf_idx, 0x15,
        0x01, 0x00, 0x80, 0x7f, 0x00, 0x00, 0x00, 0x3f, 0x01, 0x00, 0x80, 0x7f,
    ];
    assert_eq!(bcf_tail(bytes, expected.len()), expected);

    let bytes = write(&setup, OutputFormat::Bcf, &[None, None, None], &[Some(1.0), None, None]);
    let tail = bcf_tail(bytes, 6 + 15);
    assert_eq!(&tail[..6], [0x11, dp_idx, 0x11, 0x80, 0x80, 0x80], "int8 when all are missing");
}

// r[verify record_encoder.format_scalar_opt]
/// One value per sample, and a concrete value BCF reserves is refused.
#[test]
fn encode_opt_checks_count_and_reserved_values() {
    let setup = setup(2);
    let mut writer =
        Writer::new(Vec::new(), OutputFormat::Bcf).write_header(&setup.header).unwrap();
    let alleles = Alleles::snv(Base::A, Base::C).unwrap();
    let mut enc = writer
        .begin_record(&setup.chr1, Pos1::new(10).unwrap(), &alleles, None)
        .unwrap()
        .filter_pass()
        .begin_samples();
    let err = setup.dp.encode_opt(&mut enc, &[None]).unwrap_err();
    assert!(matches!(err, VcfError::SampleCountMismatch { expected: 2, got: 1 }), "{err:?}");
    let err = setup.dp.encode_opt(&mut enc, &[None, Some(i32::MIN)]).unwrap_err();
    assert!(matches!(err, VcfError::Encode(_)), "{err:?}");
    let err =
        setup.vf.encode_opt(&mut enc, &[Some(f32::from_bits(0x7F80_0001)), None]).unwrap_err();
    assert!(matches!(err, VcfError::Encode(_)), "{err:?}");
    setup.dp.encode_opt(&mut enc, &[None, Some(3)]).unwrap();
    enc.emit().unwrap();
}
