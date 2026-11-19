//! Arrays with missing elements in every section (`r[vcf_passthrough.opt_arrays]`):
//! INFO floats, FORMAT integers and FORMAT floats for every sample, written as
//! VCF text and BCF and read by bcftools against a line formatted by hand.
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
use seqair::vcf::record_encoder::{FormatFieldDef, InfoFieldDef, OptArr};
use seqair::vcf::{
    Alleles, ContigDef, FormatFloatOpts, FormatIntOpts, InfoFloatOpts, Number, OutputFormat,
    VcfError, VcfHeader, Writer,
};
use seqair_types::{Base, Pos1};
use std::process::Command;

/// Floats whose `%g` text is known exactly.
const FLOATS: [(f32, &str); 5] =
    [(0.5, "0.5"), (-2.0, "-2"), (1.25, "1.25"), (3.0, "3"), (100.0, "100")];

fn join<T>(values: &[Option<T>], text: impl Fn(&T) -> String) -> String {
    if values.is_empty() {
        return ".".to_owned();
    }
    let parts: Vec<String> =
        values.iter().map(|v| v.as_ref().map_or_else(|| ".".to_owned(), &text)).collect();
    parts.join(",")
}

struct Case {
    info: Vec<Option<usize>>,
    ints: Vec<Vec<Option<i32>>>,
    floats: Vec<Vec<Option<usize>>>,
}

fn write(case: &Case, format: OutputFormat) -> Vec<u8> {
    let mut b = VcfHeader::builder();
    let chr1 = b.register_contig("chr1", ContigDef::new(Some(1000))).unwrap();
    let mut b = b.infos();
    let fo: InfoFloatOpts = b
        .register_info(&InfoFieldDef::<OptArr<f32>>::with_number("FO", Number::Unknown, "f"))
        .unwrap();
    let mut b = b.formats();
    let io: FormatIntOpts = b
        .register_format(&FormatFieldDef::<OptArr<i32>>::with_number("IO", Number::Unknown, "i"))
        .unwrap();
    let ff: FormatFloatOpts = b
        .register_format(&FormatFieldDef::<OptArr<f32>>::with_number("FF", Number::Unknown, "f"))
        .unwrap();
    let mut b = b.samples();
    for i in 0..case.ints.len() {
        b.add_sample(format!("S{i}")).unwrap();
    }
    let header = b.build().unwrap();
    let mut out = Vec::new();
    let mut w = Writer::new(&mut out, format).write_header(&header).unwrap();
    let alleles = Alleles::snv(Base::A, Base::T).unwrap();
    let mut enc =
        w.begin_record(&chr1, Pos1::new(10).unwrap(), &alleles, None).unwrap().filter_pass();
    let info: Vec<Option<f32>> = case.info.iter().map(|v| v.map(|i| FLOATS[i].0)).collect();
    fo.encode(&mut enc, &info).unwrap();
    let mut enc = enc.begin_samples();
    let ints: Vec<&[Option<i32>]> = case.ints.iter().map(Vec::as_slice).collect();
    io.encode(&mut enc, &ints).unwrap();
    let floats: Vec<Vec<Option<f32>>> =
        case.floats.iter().map(|s| s.iter().map(|v| v.map(|i| FLOATS[i].0)).collect()).collect();
    let floats: Vec<&[Option<f32>]> = floats.iter().map(Vec::as_slice).collect();
    ff.encode(&mut enc, &floats).unwrap();
    enc.emit().unwrap();
    w.finish().unwrap();
    out
}

fn bcftools_line(bytes: &[u8], suffix: &str) -> String {
    let file = tempfile::Builder::new().suffix(suffix).tempfile().unwrap();
    std::fs::write(file.path(), bytes).unwrap();
    let out = Command::new("bcftools").args(["view", "-H"]).arg(file.path()).output().unwrap();
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    String::from_utf8(out.stdout).unwrap().trim_end().to_owned()
}

// r[verify vcf_passthrough.opt_arrays]
/// bcftools reads every element where it was written, `.` for `None`, in
/// VCF text and in BCF alike, whatever integer width the column needs.
#[hegel::test(test_cases = pinned::cases(40))]
fn arrays_with_missing_elements_read_back_in_bcftools(tc: TestCase) {
    if Command::new("bcftools").arg("--version").output().is_err() {
        return;
    }
    let n = tc.draw(gs::integers::<usize>().min_value(1).max_value(3));
    let float = || gs::optional(gs::integers::<usize>().max_value(FLOATS.len() - 1));
    let int = || {
        hegel::one_of!(
            gs::integers::<i32>().min_value(-200).max_value(200),
            gs::integers::<i32>().min_value(i32::MIN + 8)
        )
    };
    let case = Case {
        info: tc.draw(gs::vecs(float()).min_size(1).max_size(4)),
        ints: (0..n).map(|_| tc.draw(gs::vecs(gs::optional(int())).max_size(3))).collect(),
        floats: (0..n).map(|_| tc.draw(gs::vecs(float()).max_size(3))).collect(),
    };
    let float_text = |i: &usize| FLOATS[*i].1.to_owned();
    let samples: Vec<String> = (0..n)
        .map(|s| {
            let ints = join(&case.ints[s], i32::to_string);
            format!("{ints}:{}", join(&case.floats[s], float_text))
        })
        .collect();
    let expected = format!(
        "chr1\t10\t.\tA\tT\t.\tPASS\tFO={}\tIO:FF\t{}",
        join(&case.info, float_text),
        samples.join("\t")
    );
    assert_eq!(bcftools_line(&write(&case, OutputFormat::Vcf), ".vcf"), expected);
    assert_eq!(bcftools_line(&write(&case, OutputFormat::Bcf), ".bcf"), expected);
}

// r[verify vcf_passthrough.opt_arrays]
#[test]
fn a_sample_count_other_than_the_headers_is_refused() {
    let mut b = VcfHeader::builder();
    let chr1 = b.register_contig("chr1", ContigDef::new(Some(1000))).unwrap();
    let mut b = b.infos().formats();
    let io: FormatIntOpts = b
        .register_format(&FormatFieldDef::<OptArr<i32>>::with_number("IO", Number::Unknown, "i"))
        .unwrap();
    let mut b = b.samples();
    b.add_sample("S0").unwrap();
    let header = b.build().unwrap();
    let mut w = Writer::new(Vec::new(), OutputFormat::Vcf).write_header(&header).unwrap();
    let alleles = Alleles::snv(Base::A, Base::T).unwrap();
    let mut enc = w
        .begin_record(&chr1, Pos1::new(10).unwrap(), &alleles, None)
        .unwrap()
        .filter_pass()
        .begin_samples();
    let err = io.encode(&mut enc, &[&[Some(1)], &[None]]).unwrap_err();
    assert!(matches!(err, VcfError::SampleCountMismatch { expected: 1, got: 2 }), "{err:?}");
}
