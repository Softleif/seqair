//! Multi-sample FORMAT arrays, generated, read back by bcftools.
//!
//! BCF stores one integer width for a whole FORMAT column, so the encoder has
//! to scan every value of every sample before choosing it; ragged rows are
//! padded to the widest sample with an end-of-vector marker, and a missing
//! value uses a per-type sentinel rather than one universal marker. Each of
//! those has been a bug here before, and each is a rule that only bites on a
//! combination of samples — the fixtures pin three or four such combinations
//! by hand.
//!
//! The generator draws values that straddle the int8/int16/int32 boundaries
//! and rows of different lengths, and bcftools decodes what seqair wrote.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    reason = "test code"
)]

use hegel::prelude::*;
use seqair::vcf::record_encoder::{
    FormatFieldDef, FormatFloats, FormatGt, FormatInts, FormatString,
};
use seqair::vcf::{
    Alleles, ContigDef, ContigId, Genotype, Number, OutputFormat, ValueType, VcfHeader, Writer,
};
use seqair_types::{Base, Pos1};
use std::process::Command;
use std::sync::Arc;

const MAX_SAMPLES: usize = 4;

struct Setup {
    header: Arc<VcfHeader>,
    contig: ContigId,
    gt: FormatGt,
    ints: FormatInts,
    floats: FormatFloats,
    text: FormatString,
}

fn setup(n_samples: usize) -> Setup {
    let mut builder = VcfHeader::builder();
    let contig = builder.register_contig("chr1", ContigDef { length: Some(100_000) }).unwrap();
    let mut builder = builder.formats();
    let gt = builder
        .register_format(&FormatFieldDef::new("GT", Number::Count(1), ValueType::String, "GT"))
        .unwrap();
    let ints = builder
        .register_format(&FormatFieldDef::new(
            "AD",
            Number::Unknown,
            ValueType::Integer,
            "Per-sample integer list",
        ))
        .unwrap();
    let floats = builder
        .register_format(&FormatFieldDef::new(
            "VL",
            Number::Unknown,
            ValueType::Float,
            "Per-sample float list",
        ))
        .unwrap();
    let text = builder
        .register_format(&FormatFieldDef::new(
            "BC",
            Number::Count(1),
            ValueType::String,
            "Per-sample string",
        ))
        .unwrap();
    let mut builder = builder.samples();
    for i in 0..n_samples {
        builder.add_sample(format!("s{i}").as_str()).unwrap();
    }
    Setup { header: Arc::new(builder.build().unwrap()), contig, gt, ints, floats, text }
}

/// The widest BCF integer type a column needs. Drawn once per record: the
/// column's width is a maximum over every value of every sample, so drawing
/// each value independently over the whole `i32` range makes almost every
/// column an `int32` one and the narrow cases — where a missing value uses the
/// `int8` or `int16` sentinel — never get tested.
#[derive(Debug, Clone, Copy)]
enum Width {
    Int8,
    Int16,
    Int32,
}

/// A value inside `width`, usually smaller: samples in one column differ in
/// what they need, and the encoder has to take the maximum over all of them.
///
/// The ranges stop where `r[bcf_writer.smallest_int_type]` does: BCF reserves
/// the eight most negative values of each width for the missing and
/// end-of-vector sentinels, so `-2147483647` is not an `int32` — it is the
/// `int32` end-of-vector marker, and a reader handed it truncates the row.
/// Neither seqair nor htslib rejects one on the way in.
fn arb_int(width: Width) -> impl PrintableGenerator<i32> {
    let narrow = gs::integers::<i32>().min_value(-120).max_value(127);
    match width {
        Width::Int8 => hegel::one_of!(
            narrow,
            // The edges of the usable int8 range.
            gs::sampled_from(vec![-120i32, 127]),
        )
        .boxed_printable(),
        Width::Int16 => hegel::one_of!(
            narrow,
            gs::integers::<i32>().min_value(-32_760).max_value(32_767),
            gs::sampled_from(vec![-32_760i32, 32_767, 128, -121]),
        )
        .boxed_printable(),
        Width::Int32 => hegel::one_of!(
            narrow,
            gs::integers::<i32>().min_value(i32::MIN + 8).max_value(i32::MAX),
            gs::sampled_from(vec![32_768i32, -32_761, i32::MAX, i32::MIN + 8]),
        )
        .boxed_printable(),
    }
}

/// One record's worth of per-sample values.
#[derive(Debug, Clone)]
struct Record {
    pos: u32,
    /// Ragged on purpose: rows of different lengths, including empty ones.
    ints: Vec<Vec<i32>>,
    floats: Vec<Vec<f32>>,
    text: Vec<String>,
}

#[hegel::composite]
fn arb_records(tc: &TestCase, n_samples: usize) -> Vec<Record> {
    let n_records = tc.draw_silent(gs::integers::<usize>().min_value(1).max_value(5));
    let mut positions: Vec<u32> = (0..n_records)
        .map(|_| tc.draw_silent(gs::integers::<u32>().min_value(1).max_value(99_999)))
        .collect();
    positions.sort_unstable();

    positions
        .into_iter()
        .map(|pos| {
            let width =
                tc.draw_silent(gs::sampled_from(&[Width::Int8, Width::Int16, Width::Int32]));
            let row = |tc: &TestCase| tc.draw_silent(gs::vecs(arb_int(width)).max_size(4));
            let frow = |tc: &TestCase| {
                tc.draw_silent(
                    gs::vecs(gs::floats::<f32>().min_value(-1000.0).max_value(1000.0)).max_size(4),
                )
            };
            Record {
                pos,
                ints: (0..n_samples).map(|_| row(tc)).collect(),
                floats: (0..n_samples).map(|_| frow(tc)).collect(),
                text: (0..n_samples)
                    .map(|_| tc.draw_silent(gs::from_regex("[A-Za-z0-9]{0,6}")))
                    .collect(),
            }
        })
        .collect()
}

/// Write the records as BCF and return the bytes.
fn write_bcf(s: &Setup, records: &[Record]) -> Vec<u8> {
    let mut out = Vec::new();
    let writer = Writer::new(&mut out, OutputFormat::Bcf);
    let mut writer = writer.write_header(&s.header).unwrap();
    let n_samples = records.first().map_or(0, |r| r.ints.len());

    for rec in records {
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let enc = writer
            .begin_record(&s.contig, Pos1::new(rec.pos).unwrap(), &alleles, None)
            .unwrap()
            .filter_pass();
        let mut enc = enc.begin_samples();
        let gts: Vec<Genotype> = (0..n_samples).map(|_| Genotype::unphased(0, 1)).collect();
        s.gt.encode(&mut enc, &gts).unwrap();
        let ints: Vec<&[i32]> = rec.ints.iter().map(Vec::as_slice).collect();
        s.ints.encode(&mut enc, &ints).unwrap();
        let floats: Vec<&[f32]> = rec.floats.iter().map(Vec::as_slice).collect();
        s.floats.encode(&mut enc, &floats).unwrap();
        let text: Vec<&str> = rec.text.iter().map(String::as_str).collect();
        s.text.encode(&mut enc, &text).unwrap();
        enc.emit().unwrap();
    }
    writer.finish().unwrap();
    out
}

/// `bcftools query -f <fmt>` over the BCF, one line per record.
fn bcftools_query(bcf: &[u8], format: &str) -> Vec<String> {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("generated.bcf");
    std::fs::write(&path, bcf).unwrap();
    let out = Command::new("bcftools")
        .args(["query", "-f", format])
        .arg(&path)
        .output()
        .expect("bcftools not found");
    assert!(
        out.status.success(),
        "bcftools query failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap().lines().map(str::to_owned).collect()
}

/// How bcftools prints a row: values comma-separated, an empty row as `.`.
fn expect_ints(row: &[i32]) -> String {
    if row.is_empty() {
        ".".to_string()
    } else {
        row.iter().map(i32::to_string).collect::<Vec<_>>().join(",")
    }
}

// r[verify record_encoder.format_methods]
/// Every integer of every sample must come back from bcftools unchanged,
/// whatever widths the other samples in the same column needed.
#[hegel::test(test_cases = 40)]
fn format_int_arrays_round_trip_through_bcftools(tc: TestCase) {
    let n_samples = tc.draw(gs::integers::<usize>().min_value(1).max_value(MAX_SAMPLES));
    let records = tc.draw(arb_records(n_samples).print_as_debug());

    let s = setup(n_samples);
    let bcf = write_bcf(&s, &records);
    let lines = bcftools_query(&bcf, "[%AD\\t]\\n");

    assert_eq!(lines.len(), records.len(), "record count");
    for (line, rec) in lines.iter().zip(&records) {
        let want: Vec<String> = rec.ints.iter().map(|row| expect_ints(row)).collect();
        let got: Vec<&str> = line.trim_end().split('\t').collect();
        assert_eq!(got, want, "AD at pos {}", rec.pos);
    }

    // The column's width is the widest value anywhere in it; record when the
    // generator actually produced a column that needs more than one byte.
    let widest = records
        .iter()
        .flat_map(|r| r.ints.iter().flatten())
        .map(|v| v.unsigned_abs())
        .max()
        .unwrap_or(0);
    tc.event(match widest {
        0..=127 => "column fits int8",
        128..=32_767 => "column needs int16",
        _ => "column needs int32",
    });
    if records.iter().any(|r| r.ints.iter().any(Vec::is_empty)) {
        tc.event("has an empty row");
    }
    if records
        .iter()
        .any(|r| r.ints.iter().map(Vec::len).min() != r.ints.iter().map(Vec::len).max())
    {
        tc.event("has a ragged column");
    }
}

// r[verify record_encoder.format_methods]
/// The same for float rows, which are fixed-width but still padded per sample.
#[hegel::test(test_cases = 40)]
fn format_float_arrays_round_trip_through_bcftools(tc: TestCase) {
    let n_samples = tc.draw(gs::integers::<usize>().min_value(1).max_value(MAX_SAMPLES));
    let records = tc.draw(arb_records(n_samples).print_as_debug());

    let s = setup(n_samples);
    let bcf = write_bcf(&s, &records);
    let lines = bcftools_query(&bcf, "[%VL\\t]\\n");

    assert_eq!(lines.len(), records.len(), "record count");
    for (line, rec) in lines.iter().zip(&records) {
        let got: Vec<&str> = line.trim_end().split('\t').collect();
        assert_eq!(got.len(), rec.floats.len(), "sample count at pos {}", rec.pos);
        for (sample, (printed, row)) in got.iter().zip(&rec.floats).enumerate() {
            if row.is_empty() {
                assert_eq!(*printed, ".", "empty row at pos {} sample {sample}", rec.pos);
                continue;
            }
            let values: Vec<f32> = printed
                .split(',')
                .map(|v| v.parse().expect("bcftools prints parseable floats"))
                .collect();
            assert_eq!(values.len(), row.len(), "arity at pos {} sample {sample}", rec.pos);
            for (got, want) in values.iter().zip(row) {
                // BCF stores f32 exactly; bcftools prints six significant
                // digits, which is the round-trip the VCF spec asks for.
                let tolerance = want.abs() * 1e-5;
                assert!(
                    (got - want).abs() <= tolerance,
                    "pos {} sample {sample}: {got} != {want}",
                    rec.pos
                );
            }
        }
    }
}

// r[verify record_encoder.format_methods]
// r[verify vcf_writer.format_serialization]
/// seqair writes the sample columns twice over: as VCF text, and as BCF that
/// bcftools turns back into VCF text. The two paths share no encoding code, so
/// agreeing is evidence about both.
#[hegel::test(test_cases = 40)]
fn bcf_and_vcf_text_agree_on_sample_columns(tc: TestCase) {
    let n_samples = tc.draw(gs::integers::<usize>().min_value(1).max_value(MAX_SAMPLES));
    let records = tc.draw(arb_records(n_samples).print_as_debug());

    let s = setup(n_samples);

    let mut text = Vec::new();
    {
        let writer = Writer::new(&mut text, OutputFormat::Vcf);
        let mut writer = writer.write_header(&s.header).unwrap();
        for rec in &records {
            let alleles = Alleles::snv(Base::A, Base::T).unwrap();
            let enc = writer
                .begin_record(&s.contig, Pos1::new(rec.pos).unwrap(), &alleles, None)
                .unwrap()
                .filter_pass();
            let mut enc = enc.begin_samples();
            let gts: Vec<Genotype> = (0..n_samples).map(|_| Genotype::unphased(0, 1)).collect();
            s.gt.encode(&mut enc, &gts).unwrap();
            let ints: Vec<&[i32]> = rec.ints.iter().map(Vec::as_slice).collect();
            s.ints.encode(&mut enc, &ints).unwrap();
            let floats: Vec<&[f32]> = rec.floats.iter().map(Vec::as_slice).collect();
            s.floats.encode(&mut enc, &floats).unwrap();
            let strings: Vec<&str> = rec.text.iter().map(String::as_str).collect();
            s.text.encode(&mut enc, &strings).unwrap();
            enc.emit().unwrap();
        }
        writer.finish().unwrap();
    }
    let direct: Vec<String> = String::from_utf8(text)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').skip(8).collect::<Vec<_>>().join("\t"))
        .collect();

    let bcf = write_bcf(&s, &records);
    let via_bcf = bcftools_query(&bcf, "%LINE\\n")
        .into_iter()
        .map(|l| l.split('\t').skip(8).collect::<Vec<_>>().join("\t"))
        .collect::<Vec<_>>();

    assert_eq!(direct, via_bcf, "FORMAT and sample columns");
}
