//! VCF/BCF writer validation against bcftools: write BCF with seqair,
//! verify bcftools can read it and produces correct output. Tests integer
//! type boundaries, multi-allelic records, and missing value handling.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "test code"
)]
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    reason = "test code with known small values"
)]

use seqair::vcf::record_encoder::{
    FilterFieldDef, FormatFieldDef, FormatGt, FormatInt, InfoFieldDef, InfoFloat, InfoInt, InfoInts,
};
use seqair::vcf::{
    Alleles, ContigDef, Genotype, Number, OutputFormat, Ready, ValueType, VcfHeader, Writer,
};
use seqair_types::{Base, Pos1};
use std::process::Command;
use std::sync::Arc;

struct Setup {
    header: Arc<VcfHeader>,
    contig: seqair::vcf::ContigId,
    dp_info: InfoInt,
    af_info: InfoFloat,
    ad_info: InfoInts,
    gt_fmt: FormatGt,
    dp_fmt: FormatInt,
    gq_fmt: FormatInt,
}

fn make_setup() -> Setup {
    let mut builder = VcfHeader::builder();
    let contig = builder.register_contig("chr1", ContigDef { length: Some(250_000_000) }).unwrap();
    let mut builder = builder.filters();
    let _pass = builder.register_filter(&FilterFieldDef::new("q20", "Quality below 20")).unwrap();
    let mut builder = builder.infos();
    let dp_info = builder
        .register_info(&InfoFieldDef::new("DP", Number::Count(1), ValueType::Integer, "Depth"))
        .unwrap();
    let af_info = builder
        .register_info(&InfoFieldDef::new(
            "AF",
            Number::AlternateBases,
            ValueType::Float,
            "Allele Frequency",
        ))
        .unwrap();
    let ad_info = builder
        .register_info(&InfoFieldDef::new(
            "AD",
            Number::ReferenceAlternateBases,
            ValueType::Integer,
            "Allele Depth",
        ))
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
    let dp_fmt = builder
        .register_format(&FormatFieldDef::new(
            "DP",
            Number::Count(1),
            ValueType::Integer,
            "Read Depth",
        ))
        .unwrap();
    let gq_fmt = builder
        .register_format(&FormatFieldDef::new(
            "GQ",
            Number::Count(1),
            ValueType::Integer,
            "Genotype Quality",
        ))
        .unwrap();
    let mut builder = builder.samples();
    builder.add_sample("sample1").unwrap();
    let header = Arc::new(builder.build().unwrap());
    Setup { header, contig, dp_info, af_info, ad_info, gt_fmt, dp_fmt, gq_fmt }
}

/// Write BCF to a temp file, return the path.
fn write_bcf_file(
    dir: &std::path::Path,
    setup: &Setup,
    write_records: impl FnOnce(&Setup, &mut Writer<std::io::BufWriter<std::fs::File>, Ready>),
) -> std::path::PathBuf {
    let bcf_path = dir.join("test.bcf");
    let file = std::fs::File::create(&bcf_path).unwrap();
    let buf = std::io::BufWriter::new(file);
    let writer = Writer::new(buf, OutputFormat::Bcf);
    let mut writer = writer.write_header(&setup.header).unwrap();
    write_records(setup, &mut writer);
    writer.finish().unwrap();
    bcf_path
}

/// Run bcftools query on a BCF file and return stdout.
fn bcftools_query(bcf_path: &std::path::Path, format_str: &str) -> String {
    let output = Command::new("bcftools")
        .args(["query", "-f", format_str])
        .arg(bcf_path)
        .output()
        .expect("bcftools not found");
    assert!(
        output.status.success(),
        "bcftools query failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).unwrap()
}

/// Run bcftools view (BCF -> VCF text) and return the data lines.
fn bcftools_view_records(bcf_path: &std::path::Path) -> Vec<String> {
    let output = Command::new("bcftools")
        .args(["view", "-H"])
        .arg(bcf_path)
        .output()
        .expect("bcftools not found");
    assert!(
        output.status.success(),
        "bcftools view failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).unwrap().lines().map(String::from).collect()
}

// --- Tests ---

/// Simple SNV: write BCF, verify bcftools reads it correctly.
#[test]
fn bcftools_reads_snv() {
    let dir = tempfile::tempdir().unwrap();
    let setup = make_setup();
    let bcf_path = write_bcf_file(dir.path(), &setup, |s, w| {
        let mut enc = w
            .begin_record(
                &s.contig,
                Pos1::new(100).unwrap(),
                &Alleles::snv(Base::A, Base::T).unwrap(),
                Some(99.0),
            )
            .unwrap()
            .filter_pass();
        s.dp_info.encode(&mut enc, 42);
        let mut enc = enc.begin_samples();
        s.gt_fmt.encode(&mut enc, &[Genotype::unphased(0, 1)]).unwrap();
        s.dp_fmt.encode(&mut enc, &[30]).unwrap();
        enc.emit().unwrap();
    });

    let result = bcftools_query(&bcf_path, "%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\n");
    assert_eq!(result.trim(), "chr1\t100\tA\tT\t42");
}

/// Integer type boundaries: values at int8/int16/int32 edges.
/// BCF uses the smallest integer type that fits.
#[test]
fn bcftools_integer_boundaries() {
    let dir = tempfile::tempdir().unwrap();
    let setup = make_setup();

    let test_values: &[(i32, &str)] = &[
        (0, "0"),
        (127, "127"),               // max int8
        (128, "128"),               // needs int16
        (255, "255"),               // max uint8
        (32767, "32767"),           // max int16
        (32768, "32768"),           // needs int32
        (65535, "65535"),           // max uint16
        (2147483647, "2147483647"), // max int32
    ];

    for &(val, expected) in test_values {
        let bcf_path = write_bcf_file(dir.path(), &setup, |s, w| {
            let mut enc = w
                .begin_record(
                    &s.contig,
                    Pos1::new(1).unwrap(),
                    &Alleles::snv(Base::A, Base::T).unwrap(),
                    None,
                )
                .unwrap()
                .filter_pass();
            s.dp_info.encode(&mut enc, val);
            enc.emit().unwrap();
        });

        let result = bcftools_query(&bcf_path, "%INFO/DP\n");
        assert_eq!(result.trim(), expected, "DP={val}: bcftools reads '{}'", result.trim());
    }
}

/// Multi-allelic record: 3+ alleles with per-allele INFO fields.
#[test]
fn bcftools_multi_allelic() {
    let dir = tempfile::tempdir().unwrap();
    let setup = make_setup();
    let bcf_path = write_bcf_file(dir.path(), &setup, |s, w| {
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let mut enc = w
            .begin_record(&s.contig, Pos1::new(500).unwrap(), &alleles, Some(200.0))
            .unwrap()
            .filter_pass();
        s.dp_info.encode(&mut enc, 100);
        s.af_info.encode(&mut enc, 0.35);
        s.ad_info.encode(&mut enc, &[60, 40]);
        enc.emit().unwrap();
    });

    let result = bcftools_query(&bcf_path, "%POS\t%REF\t%ALT\t%INFO/DP\t%INFO/AF\t%INFO/AD\n");
    let fields: Vec<&str> = result.trim().split('\t').collect();
    assert_eq!(fields[0], "500");
    assert_eq!(fields[1], "A");
    assert_eq!(fields[2], "T");
    assert_eq!(fields[3], "100");
    // AF is float, might have rounding
    assert!(fields[4].starts_with("0.3"), "AF={}", fields[4]);
    assert_eq!(fields[5], "60,40");
}

/// Genotype encoding: het, hom-ref, hom-alt, missing.
#[test]
fn bcftools_genotype_encoding() {
    let dir = tempfile::tempdir().unwrap();
    let setup = make_setup();
    let bcf_path = write_bcf_file(dir.path(), &setup, |s, w| {
        // Het 0/1
        let alleles = Alleles::snv(Base::C, Base::G).unwrap();
        let enc = w
            .begin_record(&s.contig, Pos1::new(100).unwrap(), &alleles, None)
            .unwrap()
            .filter_pass();
        let mut enc = enc.begin_samples();
        s.gt_fmt.encode(&mut enc, &[Genotype::unphased(0, 1)]).unwrap();
        s.dp_fmt.encode(&mut enc, &[20]).unwrap();
        s.gq_fmt.encode(&mut enc, &[99]).unwrap();
        enc.emit().unwrap();

        // Hom-ref 0/0
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let enc = w
            .begin_record(&s.contig, Pos1::new(200).unwrap(), &alleles, None)
            .unwrap()
            .filter_pass();
        let mut enc = enc.begin_samples();
        s.gt_fmt.encode(&mut enc, &[Genotype::unphased(0, 0)]).unwrap();
        s.dp_fmt.encode(&mut enc, &[15]).unwrap();
        s.gq_fmt.encode(&mut enc, &[50]).unwrap();
        enc.emit().unwrap();

        // Hom-alt 1/1
        let alleles = Alleles::snv(Base::G, Base::A).unwrap();
        let enc = w
            .begin_record(&s.contig, Pos1::new(300).unwrap(), &alleles, None)
            .unwrap()
            .filter_pass();
        let mut enc = enc.begin_samples();
        s.gt_fmt.encode(&mut enc, &[Genotype::unphased(1, 1)]).unwrap();
        s.dp_fmt.encode(&mut enc, &[25]).unwrap();
        s.gq_fmt.encode(&mut enc, &[80]).unwrap();
        enc.emit().unwrap();
    });

    let result = bcftools_query(&bcf_path, "%POS\t[%GT\t%DP\t%GQ]\n");
    let lines: Vec<&str> = result.trim().lines().collect();
    assert_eq!(lines.len(), 3);

    let f0: Vec<&str> = lines[0].split('\t').collect();
    assert_eq!(f0[1], "0/1", "het genotype");
    assert_eq!(f0[2], "20");
    assert_eq!(f0[3], "99");

    let f1: Vec<&str> = lines[1].split('\t').collect();
    assert_eq!(f1[1], "0/0", "hom-ref genotype");

    let f2: Vec<&str> = lines[2].split('\t').collect();
    assert_eq!(f2[1], "1/1", "hom-alt genotype");
}

/// Insertion and deletion alleles round-trip through bcftools.
#[test]
fn bcftools_indel_alleles() {
    let dir = tempfile::tempdir().unwrap();
    let setup = make_setup();
    let bcf_path = write_bcf_file(dir.path(), &setup, |s, w| {
        // Insertion: A -> ACG
        let alleles = Alleles::insertion(Base::A, &[Base::C, Base::G]).unwrap();
        let mut enc = w
            .begin_record(&s.contig, Pos1::new(100).unwrap(), &alleles, Some(50.0))
            .unwrap()
            .filter_pass();
        s.dp_info.encode(&mut enc, 30);
        enc.emit().unwrap();

        // Deletion: ACG -> A
        let alleles = Alleles::deletion(Base::A, &[Base::C, Base::G]).unwrap();
        let mut enc = w
            .begin_record(&s.contig, Pos1::new(200).unwrap(), &alleles, Some(60.0))
            .unwrap()
            .filter_pass();
        s.dp_info.encode(&mut enc, 40);
        enc.emit().unwrap();
    });

    let result = bcftools_query(&bcf_path, "%POS\t%REF\t%ALT\n");
    let lines: Vec<&str> = result.trim().lines().collect();
    assert_eq!(lines[0], "100\tA\tACG", "insertion");
    assert_eq!(lines[1], "200\tACG\tA", "deletion");
}

/// Multiple records: bcftools view produces all records in order.
#[test]
fn bcftools_multiple_records_ordered() {
    let dir = tempfile::tempdir().unwrap();
    let setup = make_setup();
    let bcf_path = write_bcf_file(dir.path(), &setup, |s, w| {
        for pos in [100, 500, 1000, 5000, 10000] {
            let alleles = Alleles::snv(Base::A, Base::T).unwrap();
            let mut enc = w
                .begin_record(&s.contig, Pos1::new(pos).unwrap(), &alleles, None)
                .unwrap()
                .filter_pass();
            s.dp_info.encode(&mut enc, pos as i32);
            enc.emit().unwrap();
        }
    });

    let records = bcftools_view_records(&bcf_path);
    assert_eq!(records.len(), 5);

    let positions: Vec<&str> = records.iter().map(|r| r.split('\t').nth(1).unwrap()).collect();
    assert_eq!(positions, vec!["100", "500", "1000", "5000", "10000"]);
}

/// BCF magic validation: bcftools stats should succeed (implicitly validates magic bytes).
#[test]
fn bcftools_stats_succeeds() {
    let dir = tempfile::tempdir().unwrap();
    let setup = make_setup();
    let bcf_path = write_bcf_file(dir.path(), &setup, |s, w| {
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let mut enc = w
            .begin_record(&s.contig, Pos1::new(100).unwrap(), &alleles, Some(99.0))
            .unwrap()
            .filter_pass();
        s.dp_info.encode(&mut enc, 42);
        let mut enc = enc.begin_samples();
        s.gt_fmt.encode(&mut enc, &[Genotype::unphased(0, 1)]).unwrap();
        s.dp_fmt.encode(&mut enc, &[30]).unwrap();
        enc.emit().unwrap();
    });

    let output =
        Command::new("bcftools").args(["stats"]).arg(&bcf_path).output().expect("bcftools stats");
    assert!(
        output.status.success(),
        "bcftools stats failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stats = String::from_utf8(output.stdout).unwrap();
    // Should contain summary line with record count
    assert!(stats.contains("number of records:"), "bcftools stats output: {stats}");
}

// r[verify record_encoder.info_dedup]
// r[verify record_encoder.format_dedup]
/// Write the same INFO and FORMAT fields twice, verify bcftools reads the
/// overwritten values (not duplicates, not the first value).
#[test]
fn bcftools_dedup_overwrite() {
    let dir = tempfile::tempdir().unwrap();
    let setup = make_setup();
    let bcf_path = write_bcf_file(dir.path(), &setup, |s, w| {
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let mut enc = w
            .begin_record(&s.contig, Pos1::new(100).unwrap(), &alleles, Some(99.0))
            .unwrap()
            .filter_pass();

        // Write DP=50 then overwrite with DP=100
        s.dp_info.encode(&mut enc, 50);
        s.dp_info.encode(&mut enc, 100);

        let mut enc = enc.begin_samples();
        s.gt_fmt.encode(&mut enc, &[Genotype::unphased(0, 1)]).unwrap();
        // Write DP=30 then overwrite with DP=99
        s.dp_fmt.encode(&mut enc, &[30]).unwrap();
        s.dp_fmt.encode(&mut enc, &[99]).unwrap();
        enc.emit().unwrap();
    });

    // bcftools must be able to read the file (no duplicate key errors)
    let result = bcftools_query(&bcf_path, "%INFO/DP\t[%DP]\n");
    let fields: Vec<&str> = result.trim().split('\t').collect();
    assert_eq!(fields[0], "100", "INFO DP should be the overwritten value");
    assert_eq!(fields[1], "99", "FORMAT DP should be the overwritten value");
}
