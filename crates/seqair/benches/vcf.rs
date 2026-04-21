//! Criterion benchmarks for VCF/BCF writing.
//!
//! Models a realistic germline WGS callset with GATK-style fields (DP, AN,
//! AC[A], AF[A], MQ, QD, FS), mixed variant types (SNV/indel/multi-allelic),
//! mixed genotypes, and occasional `LowQual` filter failures.
//!
//! Compares seqair against htslib and noodles across three output formats
//! (plain VCF text, BGZF-compressed VCF, BCF binary), two field counts
//! ("minimal" = 1 INFO + 1 FORMAT vs "full" = 7 INFO + 3 FORMAT), and two
//! record counts (1k and 10k) to isolate per-field and per-record overhead.
//!
//! **I/O destination**: all three libraries write to `/dev/null` so that
//! timings reflect pure encoding overhead without memory-allocation or
//! filesystem differences.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::indexing_slicing, reason = "benches")]
#![allow(clippy::arithmetic_side_effects, reason = "benches")]
#![allow(clippy::cast_possible_truncation, reason = "benches")]
#![allow(clippy::cast_possible_wrap, reason = "benches")]
#![allow(clippy::type_complexity, reason = "benches")]

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use noodles::vcf::variant::io::Write as _;
use std::fs::File;
use std::io::BufWriter;
use std::sync::Arc;

use seqair::vcf::record_encoder::{FilterFieldDef, FormatFieldDef, InfoFieldDef};
use seqair::vcf::{
    Alleles, ContigDef, ContigId, FilterId, FormatGt, FormatInt, Genotype, InfoFloat, InfoFloats,
    InfoInt, InfoInts, Number, OutputFormat, Ready, ValueType, VcfHeader, Writer,
};
use seqair_types::{Base, Pos1};

const N_SAMPLES: usize = 5;
const RECORD_COUNTS: [u32; 2] = [1_000, 10_000];

/// Allele bank: (REF, ALT1, optional ALT2). Shared by all three libraries.
static RAW_ALLELES: &[(&[u8], &[u8], Option<&[u8]>)] = &[
    (b"A", b"T", None),       // SNV
    (b"C", b"G", None),       // SNV
    (b"G", b"A", None),       // SNV
    (b"T", b"C", None),       // SNV
    (b"A", b"C", None),       // SNV
    (b"AC", b"A", None),      // deletion (1 bp)
    (b"TAC", b"T", None),     // deletion (2 bp)
    (b"G", b"GAT", None),     // insertion (2 bp)
    (b"A", b"T", Some(b"C")), // multi-allelic SNV
    (b"G", b"T", None),       // SNV
];

/// Buffered `/dev/null` sink — all three libraries write to this so timings
/// reflect encoding overhead, not syscall frequency or memory allocation.
/// htslib opens `/dev/null` by path (with its own C-level stdio buffering).
fn dev_null() -> BufWriter<File> {
    BufWriter::new(File::create("/dev/null").unwrap())
}

// ===========================================================================
// seqair
// ===========================================================================

struct SeqairSetup {
    header: Arc<VcfHeader>,
    contig: ContigId,
    lowqual_filter: FilterId,
    dp_info: InfoInt,
    an_info: InfoInt,
    ac_info: InfoInts,
    af_info: InfoFloats,
    mq_info: InfoFloat,
    qd_info: InfoFloat,
    fs_info: InfoFloat,
    gt_fmt: FormatGt,
    dp_fmt: FormatInt,
    gq_fmt: FormatInt,
    alleles_bank: Vec<Alleles>,
}

fn seqair_setup() -> SeqairSetup {
    let mut b = VcfHeader::builder();
    let contig = b.register_contig("chr1", ContigDef { length: Some(250_000_000) }).unwrap();

    let mut b = b.filters();
    let lowqual_filter = b.register_filter(&FilterFieldDef::new("LowQual", "Low quality")).unwrap();

    let mut b = b.infos();
    let dp_info = b
        .register_info(&InfoFieldDef::new(
            "DP",
            Number::Count(1),
            ValueType::Integer,
            "Total read depth",
        ))
        .unwrap();
    let an_info = b
        .register_info(&InfoFieldDef::new(
            "AN",
            Number::Count(1),
            ValueType::Integer,
            "Total allele count in genotypes",
        ))
        .unwrap();
    let ac_info = b
        .register_info(&InfoFieldDef::new(
            "AC",
            Number::AlternateBases,
            ValueType::Integer,
            "Allele count in genotypes, for each ALT",
        ))
        .unwrap();
    let af_info = b
        .register_info(&InfoFieldDef::new(
            "AF",
            Number::AlternateBases,
            ValueType::Float,
            "Allele frequency from GT counts, for each ALT",
        ))
        .unwrap();
    let mq_info = b
        .register_info(&InfoFieldDef::new(
            "MQ",
            Number::Count(1),
            ValueType::Float,
            "RMS mapping quality",
        ))
        .unwrap();
    let qd_info = b
        .register_info(&InfoFieldDef::new(
            "QD",
            Number::Count(1),
            ValueType::Float,
            "Variant confidence divided by unfiltered depth",
        ))
        .unwrap();
    let fs_info = b
        .register_info(&InfoFieldDef::new(
            "FS",
            Number::Count(1),
            ValueType::Float,
            "Phred-scaled p-value using Fisher's exact test for strand bias",
        ))
        .unwrap();

    let mut b = b.formats();
    let gt_fmt = b
        .register_format(&FormatFieldDef::new(
            "GT",
            Number::Count(1),
            ValueType::String,
            "Genotype",
        ))
        .unwrap();
    let dp_fmt = b
        .register_format(&FormatFieldDef::new(
            "DP",
            Number::Count(1),
            ValueType::Integer,
            "Approximate read depth",
        ))
        .unwrap();
    let gq_fmt = b
        .register_format(&FormatFieldDef::new(
            "GQ",
            Number::Count(1),
            ValueType::Integer,
            "Genotype quality",
        ))
        .unwrap();

    let mut b = b.samples();
    for s in 1..=N_SAMPLES {
        b.add_sample(format!("SAMPLE{s:02}")).unwrap();
    }
    let header = Arc::new(b.build().unwrap());

    let alleles_bank = vec![
        Alleles::snv(Base::A, Base::T).unwrap(),
        Alleles::snv(Base::C, Base::G).unwrap(),
        Alleles::snv(Base::G, Base::A).unwrap(),
        Alleles::snv(Base::T, Base::C).unwrap(),
        Alleles::snv(Base::A, Base::C).unwrap(),
        Alleles::deletion(Base::A, &[Base::C]).unwrap(),
        Alleles::deletion(Base::T, &[Base::A, Base::C]).unwrap(),
        Alleles::insertion(Base::G, &[Base::A, Base::T]).unwrap(),
        Alleles::snv_multi(Base::A, &[Base::T, Base::C]).unwrap(),
        Alleles::snv(Base::G, Base::T).unwrap(),
    ];

    SeqairSetup {
        header,
        contig,
        lowqual_filter,
        dp_info,
        an_info,
        ac_info,
        af_info,
        mq_info,
        qd_info,
        fs_info,
        gt_fmt,
        dp_fmt,
        gq_fmt,
        alleles_bank,
    }
}

/// Encode one VCF record using the seqair typestate encoder.
///
/// When `minimal` is true, only DP (INFO) and GT (FORMAT) are written —
/// the remaining 6 INFO and 2 FORMAT fields are skipped.
fn seqair_encode_record<W: std::io::Write>(
    writer: &mut Writer<W, Ready>,
    s: &SeqairSetup,
    i: u32,
    minimal: bool,
) {
    let alleles = &s.alleles_bank[i as usize % s.alleles_bank.len()];
    let n_alt = alleles.n_allele() - 1;
    let an = (N_SAMPLES * 2) as i32;
    let pos = Pos1::new(i * 237 + 1).unwrap();

    let enc =
        writer.begin_record(&s.contig, pos, alleles, Some(10.0 + (i % 50) as f32 * 2.0)).unwrap();
    let mut enc =
        if i % 10 == 7 { enc.filter_fail(&[&s.lowqual_filter]) } else { enc.filter_pass() };

    let dp = 30 + (i % 70) as i32;
    s.dp_info.encode(&mut enc, dp);

    if !minimal {
        s.an_info.encode(&mut enc, an);
        let ac: Vec<i32> = (0..n_alt).map(|k| 1 + (i as i32 + k as i32) % 9).collect();
        let af: Vec<f32> = ac.iter().map(|&a| a as f32 / an as f32).collect();
        s.ac_info.encode(&mut enc, &ac);
        s.af_info.encode(&mut enc, &af);
        s.mq_info.encode(&mut enc, 55.0 + (i % 20) as f32 * 0.5);
        s.qd_info.encode(&mut enc, 5.0 + (i % 25) as f32 * 0.8);
        s.fs_info.encode(&mut enc, (i % 30) as f32 * 0.5);
    }

    let mut enc = enc.begin_samples();
    let gts: Vec<Genotype> = (0..N_SAMPLES)
        .map(|si| match (si + i as usize) % 6 {
            0 => Genotype::unphased(0, 0),
            1 => Genotype::unphased(0, 1),
            2 => Genotype::unphased(1, 1),
            3 => Genotype::phased_diploid(0, 1),
            4 => Genotype::missing_diploid(),
            _ => Genotype::unphased(0, 1),
        })
        .collect();
    s.gt_fmt.encode(&mut enc, &gts).unwrap();

    if !minimal {
        let dp_vals: Vec<i32> = (0..N_SAMPLES as i32).map(|si| dp / 2 + si % 20).collect();
        let gq_vals: Vec<i32> = (0..N_SAMPLES as i32).map(|si| 30 + (i as i32 + si) % 60).collect();
        s.dp_fmt.encode(&mut enc, &dp_vals).unwrap();
        s.gq_fmt.encode(&mut enc, &gq_vals).unwrap();
    }

    enc.emit().unwrap();
}

// ===========================================================================
// htslib
// ===========================================================================

fn htslib_header() -> rust_htslib::bcf::Header {
    let mut h = rust_htslib::bcf::Header::new();
    h.push_record(b"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth\">");
    h.push_record(
        b"##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total allele count in genotypes\">",
    );
    h.push_record(
        b"##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes, for each ALT\">",
    );
    h.push_record(
        b"##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency from GT counts, for each ALT\">",
    );
    h.push_record(b"##INFO=<ID=MQ,Number=1,Type=Float,Description=\"RMS mapping quality\">");
    h.push_record(
        b"##INFO=<ID=QD,Number=1,Type=Float,Description=\"Variant confidence divided by unfiltered depth\">",
    );
    h.push_record(
        b"##INFO=<ID=FS,Number=1,Type=Float,Description=\"Phred-scaled p-value using Fisher's exact test for strand bias\">",
    );
    h.push_record(b"##FILTER=<ID=LowQual,Description=\"Low quality\">");
    h.push_record(b"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    h.push_record(b"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth\">");
    h.push_record(b"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">");
    h.push_record(b"##contig=<ID=chr1,length=250000000>");
    for s in 1..=N_SAMPLES {
        h.push_sample(format!("SAMPLE{s:02}").as_bytes());
    }
    h
}

/// Write `n_records` germline-style records with rust-htslib.
fn htslib_write_records(writer: &mut rust_htslib::bcf::Writer, n_records: u32, minimal: bool) {
    use rust_htslib::bcf::record::GenotypeAllele;

    let rid = writer.header().name2rid(b"chr1").unwrap();
    let pass_id = writer.header().name_to_id(cstr8::cstr8!("PASS")).unwrap();
    let lowqual_id = writer.header().name_to_id(cstr8::cstr8!("LowQual")).unwrap();
    let an = (N_SAMPLES * 2) as i32;

    for i in 0..i64::from(n_records) {
        let allele_idx = i as usize % RAW_ALLELES.len();
        let (ref_bytes, alt1_bytes, alt2_opt) = RAW_ALLELES[allele_idx];
        let n_alt = if alt2_opt.is_some() { 2usize } else { 1 };

        let mut record = writer.empty_record();
        record.set_rid(Some(rid));
        record.set_pos(i * 237);
        record.set_qual(10.0 + (i % 50) as f32 * 2.0);
        let alleles: Vec<&[u8]> = match alt2_opt {
            Some(alt2) => vec![ref_bytes, alt1_bytes, alt2],
            None => vec![ref_bytes, alt1_bytes],
        };
        record.set_alleles(&alleles).unwrap();
        if i % 10 == 7 {
            record.push_filter(&lowqual_id).unwrap();
        } else {
            record.push_filter(&pass_id).unwrap();
        }

        let dp = 30 + (i % 70) as i32;
        record.push_info_integer(cstr8::cstr8!("DP"), &[dp]).unwrap();

        if !minimal {
            let ac: Vec<i32> = (0..n_alt).map(|k| 1 + ((i + k as i64) % 9) as i32).collect();
            let af: Vec<f32> = ac.iter().map(|&a| a as f32 / an as f32).collect();
            record.push_info_integer(cstr8::cstr8!("AN"), &[an]).unwrap();
            record.push_info_integer(cstr8::cstr8!("AC"), &ac).unwrap();
            record.push_info_float(cstr8::cstr8!("AF"), &af).unwrap();
            record
                .push_info_float(cstr8::cstr8!("MQ"), &[55.0f32 + (i % 20) as f32 * 0.5])
                .unwrap();
            record.push_info_float(cstr8::cstr8!("QD"), &[5.0f32 + (i % 25) as f32 * 0.8]).unwrap();
            record.push_info_float(cstr8::cstr8!("FS"), &[(i % 30) as f32 * 0.5]).unwrap();
        }

        let gts: Vec<GenotypeAllele> = (0..N_SAMPLES)
            .flat_map(|s| match (s + i as usize) % 6 {
                0 => [GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(0)],
                1 => [GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)],
                2 => [GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(1)],
                3 => [GenotypeAllele::Phased(0), GenotypeAllele::Phased(1)],
                4 => [GenotypeAllele::UnphasedMissing, GenotypeAllele::UnphasedMissing],
                _ => [GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)],
            })
            .collect();
        record.push_genotypes(&gts).unwrap();

        if !minimal {
            let dp_vals: Vec<i32> = (0..N_SAMPLES as i32).map(|s| dp / 2 + s % 20).collect();
            let gq_vals: Vec<i32> =
                (0..N_SAMPLES as i32).map(|s| 30 + (i as i32 + s) % 60).collect();
            record.push_format_integer(cstr8::cstr8!("DP"), &dp_vals).unwrap();
            record.push_format_integer(cstr8::cstr8!("GQ"), &gq_vals).unwrap();
        }

        writer.write(&record).unwrap();
    }
}

// ===========================================================================
// noodles
// ===========================================================================

/// Build one noodles `RecordBuf` with the same data pattern as the other libs.
fn noodles_build_record(an: i32, i: u32, minimal: bool) -> noodles::vcf::variant::RecordBuf {
    use noodles::vcf;
    use noodles::vcf::variant::record_buf::info::field::Value as InfoValue;
    use noodles::vcf::variant::record_buf::info::field::value::Array as InfoArray;
    use noodles::vcf::variant::record_buf::samples::sample::Value as SampleValue;

    let allele_idx = i as usize % RAW_ALLELES.len();
    let (ref_bytes, alt1_bytes, alt2_opt) = RAW_ALLELES[allele_idx];
    let n_alt = if alt2_opt.is_some() { 2usize } else { 1 };
    let dp = 30 + (i % 70) as i32;
    let pos = noodles::core::Position::try_from(i as usize * 237 + 1).unwrap();

    let mut alt_bases = vec![std::str::from_utf8(alt1_bytes).unwrap().to_owned()];
    if let Some(alt2) = alt2_opt {
        alt_bases.push(std::str::from_utf8(alt2).unwrap().to_owned());
    }
    let filters = if i % 10 == 7 {
        [String::from("LowQual")].into_iter().collect()
    } else {
        vcf::variant::record_buf::Filters::pass()
    };

    let info: vcf::variant::record_buf::Info = if minimal {
        [("DP".parse().unwrap(), Some(InfoValue::Integer(dp)))].into_iter().collect()
    } else {
        let ac: Vec<Option<i32>> =
            (0..n_alt).map(|k| Some(1 + (i as i32 + k as i32) % 9)).collect();
        let af: Vec<Option<f32>> = ac.iter().map(|&a| a.map(|v| v as f32 / an as f32)).collect();
        [
            ("DP".parse().unwrap(), Some(InfoValue::Integer(dp))),
            ("AN".parse().unwrap(), Some(InfoValue::Integer(an))),
            ("AC".parse().unwrap(), Some(InfoValue::Array(InfoArray::Integer(ac)))),
            ("AF".parse().unwrap(), Some(InfoValue::Array(InfoArray::Float(af)))),
            ("MQ".parse().unwrap(), Some(InfoValue::Float(55.0 + (i % 20) as f32 * 0.5))),
            ("QD".parse().unwrap(), Some(InfoValue::Float(5.0 + (i % 25) as f32 * 0.8))),
            ("FS".parse().unwrap(), Some(InfoValue::Float((i % 30) as f32 * 0.5))),
        ]
        .into_iter()
        .collect()
    };

    let gt_str = |s: usize| match (s + i as usize) % 6 {
        0 => "0/0",
        1 => "0/1",
        2 => "1/1",
        3 => "0|1",
        4 => "./.",
        _ => "0/1",
    };

    let samples = if minimal {
        let format_keys = ["GT".parse().unwrap()].into_iter().collect();
        let rows = (0..N_SAMPLES)
            .map(|s| {
                [Some(SampleValue::Genotype(gt_str(s).parse().unwrap()))].into_iter().collect()
            })
            .collect();
        vcf::variant::record_buf::Samples::new(format_keys, rows)
    } else {
        let format_keys = ["GT".parse().unwrap(), "DP".parse().unwrap(), "GQ".parse().unwrap()]
            .into_iter()
            .collect();
        let rows = (0..N_SAMPLES)
            .map(|s| {
                [
                    Some(SampleValue::Genotype(gt_str(s).parse().unwrap())),
                    Some(SampleValue::Integer(dp / 2 + s as i32 % 20)),
                    Some(SampleValue::Integer(30 + (i as i32 + s as i32) % 60)),
                ]
                .into_iter()
                .collect()
            })
            .collect();
        vcf::variant::record_buf::Samples::new(format_keys, rows)
    };

    vcf::variant::RecordBuf::builder()
        .set_reference_sequence_name("chr1")
        .set_variant_start(pos)
        .set_reference_bases(std::str::from_utf8(ref_bytes).unwrap())
        .set_alternate_bases(vcf::variant::record_buf::AlternateBases::from(alt_bases))
        .set_quality_score(10.0 + (i % 50) as f32 * 2.0)
        .set_filters(filters)
        .set_info(info)
        .set_samples(samples)
        .build()
}

// ===========================================================================
// Benchmark groups
// ===========================================================================

fn write_vcf_text(c: &mut Criterion) {
    let setup = seqair_setup();
    let noodles_hdr: noodles::vcf::Header = setup.header.to_vcf_text().parse().unwrap();
    let htslib_hdr = htslib_header();
    let an = (N_SAMPLES * 2) as i32;

    for n_records in RECORD_COUNTS {
        let mut group = c.benchmark_group(format!("vcf_text_{}k", n_records / 1000));
        group.throughput(Throughput::Elements(u64::from(n_records)));

        for minimal in [true, false] {
            let label = if minimal { "minimal" } else { "full" };

            group.bench_function(BenchmarkId::new("seqair", label), |b| {
                b.iter(|| {
                    let sink = dev_null();
                    let writer = Writer::new(sink, OutputFormat::Vcf);
                    let mut writer = writer.write_header(&setup.header).unwrap();
                    for i in 0u32..n_records {
                        seqair_encode_record(&mut writer, &setup, i, minimal);
                    }
                    writer.finish().unwrap();
                });
            });

            // uncompressed=true → mode "w" → plain VCF text (not BGZF).
            group.bench_function(BenchmarkId::new("htslib", label), |b| {
                b.iter(|| {
                    let mut w = rust_htslib::bcf::Writer::from_path(
                        "/dev/null",
                        &htslib_hdr,
                        true,
                        rust_htslib::bcf::Format::Vcf,
                    )
                    .unwrap();
                    htslib_write_records(&mut w, n_records, minimal);
                });
            });

            group.bench_function(BenchmarkId::new("noodles", label), |b| {
                b.iter(|| {
                    let sink = dev_null();
                    let mut w = noodles::vcf::io::Writer::new(sink);
                    w.write_header(&noodles_hdr).unwrap();
                    for i in 0u32..n_records {
                        let record = noodles_build_record(an, i, minimal);
                        w.write_variant_record(&noodles_hdr, &record).unwrap();
                    }
                });
            });
        }

        group.finish();
    }
}

fn write_vcf_gz(c: &mut Criterion) {
    let setup = seqair_setup();
    let noodles_hdr: noodles::vcf::Header = setup.header.to_vcf_text().parse().unwrap();
    let htslib_hdr = htslib_header();
    let an = (N_SAMPLES * 2) as i32;

    for n_records in RECORD_COUNTS {
        let mut group = c.benchmark_group(format!("vcf_gz_{}k", n_records / 1000));
        group.throughput(Throughput::Elements(u64::from(n_records)));

        for minimal in [true, false] {
            let label = if minimal { "minimal" } else { "full" };

            group.bench_function(BenchmarkId::new("seqair", label), |b| {
                b.iter(|| {
                    let sink = dev_null();
                    let writer = Writer::new(sink, OutputFormat::VcfGz);
                    let mut writer = writer.write_header(&setup.header).unwrap();
                    for i in 0u32..n_records {
                        seqair_encode_record(&mut writer, &setup, i, minimal);
                    }
                    writer.finish().unwrap();
                });
            });

            // uncompressed=false → mode "wz" → BGZF-compressed VCF.
            group.bench_function(BenchmarkId::new("htslib", label), |b| {
                b.iter(|| {
                    let mut w = rust_htslib::bcf::Writer::from_path(
                        "/dev/null",
                        &htslib_hdr,
                        false,
                        rust_htslib::bcf::Format::Vcf,
                    )
                    .unwrap();
                    htslib_write_records(&mut w, n_records, minimal);
                });
            });

            group.bench_function(BenchmarkId::new("noodles", label), |b| {
                b.iter(|| {
                    let sink = dev_null();
                    let bgzf = noodles_bgzf::io::Writer::new(sink);
                    let mut w = noodles::vcf::io::Writer::new(bgzf);
                    w.write_header(&noodles_hdr).unwrap();
                    for i in 0u32..n_records {
                        let record = noodles_build_record(an, i, minimal);
                        w.write_variant_record(&noodles_hdr, &record).unwrap();
                    }
                    drop(w);
                });
            });
        }

        group.finish();
    }
}

fn write_bcf(c: &mut Criterion) {
    let setup = seqair_setup();
    let noodles_hdr: noodles::vcf::Header = setup.header.to_vcf_text().parse().unwrap();
    let htslib_hdr = htslib_header();
    let an = (N_SAMPLES * 2) as i32;

    for n_records in RECORD_COUNTS {
        let mut group = c.benchmark_group(format!("bcf_{}k", n_records / 1000));
        group.throughput(Throughput::Elements(u64::from(n_records)));

        for minimal in [true, false] {
            let label = if minimal { "minimal" } else { "full" };

            group.bench_function(BenchmarkId::new("seqair", label), |b| {
                b.iter(|| {
                    let sink = dev_null();
                    let writer = Writer::new(sink, OutputFormat::Bcf);
                    let mut writer = writer.write_header(&setup.header).unwrap();
                    for i in 0u32..n_records {
                        seqair_encode_record(&mut writer, &setup, i, minimal);
                    }
                    writer.finish().unwrap();
                });
            });

            // uncompressed=false → mode "wb" → compressed BCF.
            group.bench_function(BenchmarkId::new("htslib", label), |b| {
                b.iter(|| {
                    let mut w = rust_htslib::bcf::Writer::from_path(
                        "/dev/null",
                        &htslib_hdr,
                        false,
                        rust_htslib::bcf::Format::Bcf,
                    )
                    .unwrap();
                    htslib_write_records(&mut w, n_records, minimal);
                });
            });

            group.bench_function(BenchmarkId::new("noodles", label), |b| {
                b.iter(|| {
                    let sink = dev_null();
                    let mut w = noodles::bcf::io::Writer::new(sink);
                    w.write_header(&noodles_hdr).unwrap();
                    for i in 0u32..n_records {
                        let record = noodles_build_record(an, i, minimal);
                        w.write_variant_record(&noodles_hdr, &record).unwrap();
                    }
                    drop(w);
                });
            });
        }

        group.finish();
    }
}

criterion_group!(benches, write_vcf_text, write_vcf_gz, write_bcf);
criterion_main!(benches);
