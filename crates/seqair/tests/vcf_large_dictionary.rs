//! Headers with far more FILTER/INFO/FORMAT ids than the 64 the string
//! dictionary once allowed, written by seqair as VCF and BCF and read back by
//! bcftools and noodles.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "test code"
)]

use noodles::vcf::variant::io::Write as _;
use seqair::vcf::record_encoder::{FormatFieldDef, InfoFieldDef};
use seqair::vcf::{
    Alleles, ContigDef, ContigId, FormatInt, InfoInt, Number, OutputFormat, ValueType, VcfHeader,
    Writer,
};
use seqair_types::{Base, Pos1};
use std::io::Cursor;
use std::process::Command;

const N_INFO: usize = 320;
const N_FORMAT: usize = 80;
/// Written once with a wrong value first, so it is a duplicate within the record.
const DUP_INFO: usize = 300;
const DUP_FORMAT: usize = 70;

struct Setup {
    header: VcfHeader,
    contig: ContigId,
    infos: Vec<InfoInt>,
    formats: Vec<FormatInt>,
}

fn leak(name: String) -> &'static str {
    Box::leak(name.into_boxed_str())
}

fn setup() -> Setup {
    setup_with(N_INFO, N_FORMAT)
}

fn setup_with(n_info: usize, n_format: usize) -> Setup {
    let mut builder = VcfHeader::builder();
    let contig = builder.register_contig("chr1", ContigDef { length: Some(1_000) }).unwrap();
    let mut builder = builder.infos();
    let infos = (0..n_info)
        .map(|k| {
            let def =
                InfoFieldDef::new(leak(format!("I{k}")), Number::Count(1), ValueType::Integer, "i");
            builder.register_info(&def).unwrap()
        })
        .collect();
    let mut builder = builder.formats();
    // PASS and the INFO ids come first, so every FORMAT id is past 320.
    let formats = (0..n_format)
        .map(|k| {
            let def = FormatFieldDef::new(
                leak(format!("F{k}")),
                Number::Count(1),
                ValueType::Integer,
                "f",
            );
            builder.register_format(&def).unwrap()
        })
        .collect();
    let mut builder = builder.samples();
    builder.add_sample("S1").unwrap();
    builder.add_sample("S2").unwrap();
    Setup { header: builder.build().unwrap(), contig, infos, formats }
}

fn info_value(k: usize) -> i32 {
    i32::try_from(k * 7).unwrap()
}

fn format_value(k: usize, sample: usize) -> i32 {
    i32::try_from(k * 10 + sample).unwrap()
}

/// One parsed data line: INFO and FORMAT `(key, value)` pairs in column order,
/// FORMAT values as one string per sample.
type Line = (Vec<(String, String)>, Vec<(String, Vec<String>)>);

/// What the two records below must read back as. The duplicates keep the
/// position of their first write and the value of their last.
fn expected() -> Vec<Line> {
    let mut info_order = vec![DUP_INFO];
    info_order.extend((0..N_INFO).filter(|&k| k != DUP_INFO));
    let mut format_order = vec![DUP_FORMAT];
    format_order.extend((0..N_FORMAT).filter(|&k| k != DUP_FORMAT));
    let first = (
        info_order.iter().map(|&k| (format!("I{k}"), info_value(k).to_string())).collect(),
        format_order
            .iter()
            .map(|&k| (format!("F{k}"), (0..2).map(|s| format_value(k, s).to_string()).collect()))
            .collect(),
    );
    // The second record writes ids the first already used: the per-record
    // sets must have been cleared.
    let second = (
        vec![("I5".into(), "5".into()), ("I310".into(), "310".into())],
        vec![("F75".into(), vec!["1".into(), "2".into()])],
    );
    vec![first, second]
}

fn write(setup: &Setup, format: OutputFormat) -> Vec<u8> {
    let alleles = Alleles::snv(Base::A, Base::T).unwrap();
    let mut out = Vec::new();
    let mut writer = Writer::new(&mut out, format).write_header(&setup.header).unwrap();

    let mut enc = writer
        .begin_record(&setup.contig, Pos1::new(100).unwrap(), &alleles, None)
        .unwrap()
        .filter_pass();
    setup.infos[DUP_INFO].encode(&mut enc, -1).unwrap();
    for (k, key) in setup.infos.iter().enumerate() {
        key.encode(&mut enc, info_value(k)).unwrap();
    }
    let mut enc = enc.begin_samples();
    setup.formats[DUP_FORMAT].encode(&mut enc, &[-1, -1]).unwrap();
    for (k, key) in setup.formats.iter().enumerate() {
        key.encode(&mut enc, &[format_value(k, 0), format_value(k, 1)]).unwrap();
    }
    enc.emit().unwrap();

    let mut enc = writer
        .begin_record(&setup.contig, Pos1::new(200).unwrap(), &alleles, None)
        .unwrap()
        .filter_pass();
    setup.infos[5].encode(&mut enc, 5).unwrap();
    setup.infos[310].encode(&mut enc, 310).unwrap();
    let mut enc = enc.begin_samples();
    setup.formats[75].encode(&mut enc, &[1, 2]).unwrap();
    enc.emit().unwrap();

    writer.finish().unwrap();
    out
}

/// Parse the INFO, FORMAT and sample columns of every data line of VCF text.
fn parse(vcf: &str) -> Vec<Line> {
    vcf.lines()
        .filter(|l| !l.starts_with('#'))
        .map(|line| {
            let cols: Vec<&str> = line.split('\t').collect();
            assert_eq!(cols.len(), 11, "8 fixed columns, FORMAT and two samples: {line}");
            let info = (cols[7] != ".")
                .then(|| cols[7].split(';'))
                .into_iter()
                .flatten()
                .map(|kv| {
                    let (k, v) = kv.split_once('=').expect("every INFO field has a value");
                    (k.to_owned(), v.to_owned())
                })
                .collect();
            let samples: Vec<Vec<&str>> =
                cols[9..].iter().map(|s| s.split(':').collect()).collect();
            let format = cols[8]
                .split(':')
                .enumerate()
                .map(|(i, key)| (key.to_owned(), samples.iter().map(|s| s[i].to_owned()).collect()))
                .collect();
            (info, format)
        })
        .collect()
}

fn bcftools_view(bytes: &[u8], suffix: &str) -> String {
    let file = tempfile::Builder::new().suffix(suffix).tempfile().unwrap();
    std::fs::write(file.path(), bytes).unwrap();
    let out = Command::new("bcftools").arg("view").arg(file.path()).output().expect("bcftools");
    assert!(
        out.status.success(),
        "bcftools view failed: {}",
        String::from_utf8(out.stderr).unwrap()
    );
    String::from_utf8(out.stdout).unwrap()
}

/// noodles reads the file and writes it back out as VCF text.
fn noodles_view(bytes: &[u8], format: OutputFormat) -> String {
    let mut text = noodles::vcf::io::Writer::new(Vec::new());
    match format {
        OutputFormat::Bcf => {
            let mut reader = noodles::bcf::io::Reader::new(Cursor::new(bytes));
            let header = reader.read_header().unwrap();
            text.write_header(&header).unwrap();
            for record in reader.record_bufs(&header) {
                text.write_variant_record(&header, &record.unwrap()).unwrap();
            }
        }
        _ => {
            let mut reader = noodles::vcf::io::Reader::new(Cursor::new(bytes));
            let header = reader.read_header().unwrap();
            text.write_header(&header).unwrap();
            for record in reader.record_bufs(&header) {
                text.write_variant_record(&header, &record.unwrap()).unwrap();
            }
        }
    }
    String::from_utf8(text.into_inner()).unwrap()
}

// r[verify vcf_header.dict_capacity+2]
// r[verify record_encoder.info_dedup+2]
// r[verify record_encoder.format_dedup+2]
#[test]
fn hundreds_of_fields_read_back_through_bcftools() {
    let setup = setup();
    for (format, suffix) in [(OutputFormat::Vcf, ".vcf"), (OutputFormat::Bcf, ".bcf")] {
        let bytes = write(&setup, format);
        assert_eq!(parse(&bcftools_view(&bytes, suffix)), expected(), "{format:?}");
    }
}

// r[verify vcf_header.dict_capacity+2]
// r[verify record_encoder.info_dedup+2]
// r[verify record_encoder.format_dedup+2]
#[test]
fn hundreds_of_fields_read_back_through_noodles() {
    let setup = setup();
    for format in [OutputFormat::Vcf, OutputFormat::Bcf] {
        let bytes = write(&setup, format);
        assert_eq!(parse(&noodles_view(&bytes, format)), expected(), "{format:?}");
    }
}
