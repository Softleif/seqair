//! Criterion benchmarks for VCF/BCF writing.
//!
//! Models a realistic germline WGS callset: 10 samples, GATK-style fields
//! (DP, AN, AC[A], AF[A], MQ, QD, FS), mixed variant types (SNV/indel/multi-
//! allelic), mixed genotypes (0/0, 0/1, 1/1, 0|1, ./.), and occasional
//! `LowQual` filter failures.
//!
//! Compares seqair against htslib and noodles across three output formats:
//! plain VCF text, BGZF-compressed VCF, and BCF binary.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::indexing_slicing, reason = "benches")]
#![allow(clippy::arithmetic_side_effects, reason = "benches")]
#![allow(clippy::cast_possible_truncation, reason = "benches")]
#![allow(clippy::cast_possible_wrap, reason = "benches")]
#![allow(clippy::type_complexity, reason = "benches")]

use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use std::hint::black_box;
use std::sync::Arc;

use seqair::vcf::record_encoder::{FilterFieldDef, FormatFieldDef, InfoFieldDef};
use seqair::vcf::{
    Alleles, ContigDef, ContigId, FilterId, FormatGt, FormatInt, Genotype, InfoFloat, InfoFloats,
    InfoInt, InfoInts, Number, OutputFormat, ValueType, VcfHeader, Writer,
};
use seqair_types::Base;

const N_RECORDS: u32 = 10_000;
const N_GERMLINE_SAMPLES: usize = 10;

/// Allele bank for htslib/noodles: (REF bytes, ALT1 bytes, optional ALT2 bytes).
/// Matches the seqair `GermlineSetup::alleles_bank` index-for-index.
static HTSLIB_ALLELES: &[(&[u8], &[u8], Option<&[u8]>)] = &[
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

struct GermlineSetup {
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

fn germline_setup() -> GermlineSetup {
    let mut builder = VcfHeader::builder();
    let contig = builder.register_contig("chr1", ContigDef { length: Some(250_000_000) }).unwrap();
    let lowqual_filter =
        builder.register_filter(&FilterFieldDef::new("LowQual", "Low quality")).unwrap();
    let dp_info = builder
        .register_info(&InfoFieldDef::new(
            "DP",
            Number::Count(1),
            ValueType::Integer,
            "Total read depth",
        ))
        .unwrap();
    let an_info = builder
        .register_info(&InfoFieldDef::new(
            "AN",
            Number::Count(1),
            ValueType::Integer,
            "Total allele count in genotypes",
        ))
        .unwrap();
    let ac_info = builder
        .register_info(&InfoFieldDef::new(
            "AC",
            Number::AlternateBases,
            ValueType::Integer,
            "Allele count in genotypes, for each ALT",
        ))
        .unwrap();
    let af_info = builder
        .register_info(&InfoFieldDef::new(
            "AF",
            Number::AlternateBases,
            ValueType::Float,
            "Allele frequency from GT counts, for each ALT",
        ))
        .unwrap();
    let mq_info = builder
        .register_info(&InfoFieldDef::new(
            "MQ",
            Number::Count(1),
            ValueType::Float,
            "RMS mapping quality",
        ))
        .unwrap();
    let qd_info = builder
        .register_info(&InfoFieldDef::new(
            "QD",
            Number::Count(1),
            ValueType::Float,
            "Variant confidence divided by unfiltered depth",
        ))
        .unwrap();
    let fs_info = builder
        .register_info(&InfoFieldDef::new(
            "FS",
            Number::Count(1),
            ValueType::Float,
            "Phred-scaled p-value using Fisher's exact test for strand bias",
        ))
        .unwrap();
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
            "Approximate read depth",
        ))
        .unwrap();
    let gq_fmt = builder
        .register_format(&FormatFieldDef::new(
            "GQ",
            Number::Count(1),
            ValueType::Integer,
            "Genotype quality",
        ))
        .unwrap();
    let mut builder = builder.add_sample("SAMPLE01").unwrap();
    for s in 2..=N_GERMLINE_SAMPLES {
        builder = builder.add_sample(format!("SAMPLE{s:02}")).unwrap();
    }
    let header = Arc::new(builder.build().unwrap());
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
    GermlineSetup {
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

fn htslib_germline_header() -> rust_htslib::bcf::Header {
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
    for s in 1..=N_GERMLINE_SAMPLES {
        h.push_sample(format!("SAMPLE{s:02}").as_bytes());
    }
    h
}

/// Write `N_RECORDS` germline-style records with rust-htslib to a temp file.
fn write_htslib_germline(format: rust_htslib::bcf::Format) -> tempfile::NamedTempFile {
    use rust_htslib::bcf;
    use rust_htslib::bcf::record::GenotypeAllele;

    let tmp = tempfile::NamedTempFile::new().unwrap();
    let hdr = htslib_germline_header();
    let mut writer = bcf::Writer::from_path(tmp.path(), &hdr, false, format).unwrap();
    let rid = writer.header().name2rid(b"chr1").unwrap();
    let pass_id = writer.header().name_to_id(cstr8::cstr8!("PASS")).unwrap();
    let lowqual_id = writer.header().name_to_id(cstr8::cstr8!("LowQual")).unwrap();
    let an = (N_GERMLINE_SAMPLES * 2) as i32;

    for i in 0..i64::from(N_RECORDS) {
        let allele_idx = i as usize % HTSLIB_ALLELES.len();
        let (ref_bytes, alt1_bytes, alt2_opt) = HTSLIB_ALLELES[allele_idx];
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
        let ac: Vec<i32> = (0..n_alt).map(|k| 1 + ((i + k as i64) % 9) as i32).collect();
        let af: Vec<f32> = ac.iter().map(|&a| a as f32 / an as f32).collect();
        record.push_info_integer(cstr8::cstr8!("DP"), &[dp]).unwrap();
        record.push_info_integer(cstr8::cstr8!("AN"), &[an]).unwrap();
        record.push_info_integer(cstr8::cstr8!("AC"), &ac).unwrap();
        record.push_info_float(cstr8::cstr8!("AF"), &af).unwrap();
        record.push_info_float(cstr8::cstr8!("MQ"), &[55.0f32 + (i % 20) as f32 * 0.5]).unwrap();
        record.push_info_float(cstr8::cstr8!("QD"), &[5.0f32 + (i % 25) as f32 * 0.8]).unwrap();
        record.push_info_float(cstr8::cstr8!("FS"), &[(i % 30) as f32 * 0.5]).unwrap();
        // Genotypes: N_GERMLINE_SAMPLES × 2 alleles flattened
        let gts: Vec<GenotypeAllele> = (0..N_GERMLINE_SAMPLES)
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
        let dp_vals: Vec<i32> = (0..N_GERMLINE_SAMPLES as i32).map(|s| dp / 2 + s % 20).collect();
        let gq_vals: Vec<i32> =
            (0..N_GERMLINE_SAMPLES as i32).map(|s| 30 + (i as i32 + s) % 60).collect();
        record.push_format_integer(cstr8::cstr8!("DP"), &dp_vals).unwrap();
        record.push_format_integer(cstr8::cstr8!("GQ"), &gq_vals).unwrap();
        writer.write(&record).unwrap();
    }
    drop(writer);
    tmp
}

// Shared inner loop for seqair: encode one record into the open writer.
// Called identically from all three format benchmarks.
macro_rules! seqair_encode_record {
    ($writer:expr, $setup:expr, $an:expr, $i:expr) => {{
        let GermlineSetup {
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
            ..
        } = $setup;
        let i: u32 = $i;
        let alleles = &alleles_bank[i as usize % alleles_bank.len()];
        let n_alt = alleles.n_allele() - 1;
        let pos = seqair_types::Pos1::new(i * 237 + 1).unwrap();
        let enc =
            $writer.begin_record(contig, pos, alleles, Some(10.0 + (i % 50) as f32 * 2.0)).unwrap();
        let mut enc =
            if i % 10 == 7 { enc.filter_fail(&[lowqual_filter]) } else { enc.filter_pass() };
        let dp = 30 + (i % 70) as i32;
        dp_info.encode(&mut enc, dp);
        an_info.encode(&mut enc, $an);
        let ac: Vec<i32> = (0..n_alt).map(|k| 1 + (i as i32 + k as i32) % 9).collect();
        let af: Vec<f32> = ac.iter().map(|&a| a as f32 / $an as f32).collect();
        ac_info.encode(&mut enc, &ac);
        af_info.encode(&mut enc, &af);
        mq_info.encode(&mut enc, 55.0 + (i % 20) as f32 * 0.5);
        qd_info.encode(&mut enc, 5.0 + (i % 25) as f32 * 0.8);
        fs_info.encode(&mut enc, (i % 30) as f32 * 0.5);
        let mut enc = enc.begin_samples(N_GERMLINE_SAMPLES as u32);
        let gts: Vec<Genotype> = (0..N_GERMLINE_SAMPLES)
            .map(|s| match (s + i as usize) % 6 {
                0 => Genotype::unphased(0, 0),
                1 => Genotype::unphased(0, 1),
                2 => Genotype::unphased(1, 1),
                3 => Genotype::phased_diploid(0, 1),
                4 => Genotype::missing_diploid(),
                _ => Genotype::unphased(0, 1),
            })
            .collect();
        let dp_vals: Vec<i32> = (0..N_GERMLINE_SAMPLES as i32).map(|s| dp / 2 + s % 20).collect();
        let gq_vals: Vec<i32> =
            (0..N_GERMLINE_SAMPLES as i32).map(|s| 30 + (i as i32 + s) % 60).collect();
        gt_fmt.encode(&mut enc, &gts);
        dp_fmt.encode(&mut enc, &dp_vals);
        gq_fmt.encode(&mut enc, &gq_vals);
        enc.emit().unwrap();
    }};
}

// Shared inner loop for noodles: build and write one RecordBuf.
macro_rules! noodles_write_record {
    ($writer:expr, $header:expr, $an:expr, $i:expr) => {{
        use noodles::vcf;
        use noodles::vcf::variant::io::Write as _;
        use noodles::vcf::variant::record_buf::info::field::Value as InfoValue;
        use noodles::vcf::variant::record_buf::info::field::value::Array as InfoArray;
        use noodles::vcf::variant::record_buf::samples::sample::Value as SampleValue;

        let i: u32 = $i;
        let an: i32 = $an;
        let allele_idx = i as usize % HTSLIB_ALLELES.len();
        let (ref_bytes, alt1_bytes, alt2_opt) = HTSLIB_ALLELES[allele_idx];
        let n_alt = if alt2_opt.is_some() { 2usize } else { 1 };
        let dp = 30 + (i % 70) as i32;
        let ac: Vec<Option<i32>> =
            (0..n_alt).map(|k| Some(1 + (i as i32 + k as i32) % 9)).collect();
        let af: Vec<Option<f32>> = ac.iter().map(|&a| a.map(|v| v as f32 / an as f32)).collect();
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
        let info = [
            ("DP".parse().unwrap(), Some(InfoValue::Integer(dp))),
            ("AN".parse().unwrap(), Some(InfoValue::Integer(an))),
            ("AC".parse().unwrap(), Some(InfoValue::Array(InfoArray::Integer(ac)))),
            ("AF".parse().unwrap(), Some(InfoValue::Array(InfoArray::Float(af)))),
            ("MQ".parse().unwrap(), Some(InfoValue::Float(55.0 + (i % 20) as f32 * 0.5))),
            ("QD".parse().unwrap(), Some(InfoValue::Float(5.0 + (i % 25) as f32 * 0.8))),
            ("FS".parse().unwrap(), Some(InfoValue::Float((i % 30) as f32 * 0.5))),
        ]
        .into_iter()
        .collect();
        let format_keys = ["GT".parse().unwrap(), "DP".parse().unwrap(), "GQ".parse().unwrap()]
            .into_iter()
            .collect();
        let sample_rows = (0..N_GERMLINE_SAMPLES)
            .map(|s| {
                let gt_str = match (s + i as usize) % 6 {
                    0 => "0/0",
                    1 => "0/1",
                    2 => "1/1",
                    3 => "0|1",
                    4 => "./.",
                    _ => "0/1",
                };
                [
                    Some(SampleValue::Genotype(gt_str.parse().unwrap())),
                    Some(SampleValue::Integer(dp / 2 + s as i32 % 20)),
                    Some(SampleValue::Integer(30 + (i as i32 + s as i32) % 60)),
                ]
                .into_iter()
                .collect()
            })
            .collect();
        let record = vcf::variant::RecordBuf::builder()
            .set_reference_sequence_name(std::str::from_utf8(ref_bytes).unwrap())
            .set_variant_start(pos)
            .set_reference_bases(std::str::from_utf8(ref_bytes).unwrap())
            .set_alternate_bases(vcf::variant::record_buf::AlternateBases::from(alt_bases))
            .set_quality_score(10.0 + (i % 50) as f32 * 2.0)
            .set_filters(filters)
            .set_info(info)
            .set_samples(vcf::variant::record_buf::Samples::new(format_keys, sample_rows))
            .build();
        $writer.write_variant_record($header, &record).unwrap();
    }};
}

// ---------------------------------------------------------------------------
// Group 5: Write VCF text (to memory buffer)
// seqair vs htslib vs noodles — plain uncompressed VCF, 10 samples
// ---------------------------------------------------------------------------

fn write_vcf_text(c: &mut Criterion) {
    let mut group = c.benchmark_group("write_vcf_text");
    group.throughput(Throughput::Elements(u64::from(N_RECORDS)));

    let setup = germline_setup();
    let an = (N_GERMLINE_SAMPLES * 2) as i32;

    // seqair: typestate encoder, zero-alloc, pre-resolved keys
    group.bench_function("seqair", |b| {
        b.iter(|| {
            let mut output = Vec::with_capacity(4_000_000);
            let writer = Writer::new(&mut output, OutputFormat::Vcf);
            let mut writer = writer.write_header(&setup.header).unwrap();
            for i in 1u32..=N_RECORDS {
                seqair_encode_record!(writer, &setup, an, i);
            }
            writer.finish().unwrap();
            black_box(output.len())
        });
    });

    // htslib: writes to temp file (htslib VCF writer requires a path)
    group.bench_function("htslib", |b| {
        b.iter(|| {
            let tmp = write_htslib_germline(rust_htslib::bcf::Format::Vcf);
            black_box(std::fs::metadata(tmp.path()).unwrap().len())
        });
    });

    // noodles: VCF text writer via RecordBuf
    group.bench_function("noodles", |b| {
        use noodles::vcf;
        let noodles_header: vcf::Header = setup.header.to_vcf_text().parse().unwrap();
        b.iter(|| {
            let mut output = Vec::with_capacity(4_000_000);
            let mut writer = vcf::io::Writer::new(&mut output);
            writer.write_header(&noodles_header).unwrap();
            for i in 1u32..=N_RECORDS {
                noodles_write_record!(writer, &noodles_header, an, i);
            }
            black_box(output.len())
        });
    });

    group.finish();
}

// ---------------------------------------------------------------------------
// Group 6: Write .vcf.gz (BGZF-compressed VCF), 10 samples
// seqair vs htslib — tests BGZF compression overhead on top of encoding
// ---------------------------------------------------------------------------

fn write_vcf_gz(c: &mut Criterion) {
    let mut group = c.benchmark_group("write_vcf_gz");
    group.throughput(Throughput::Elements(u64::from(N_RECORDS)));

    let setup = germline_setup();
    let an = (N_GERMLINE_SAMPLES * 2) as i32;

    // seqair: BGZF-compressed VCF to memory buffer
    group.bench_function("seqair", |b| {
        b.iter(|| {
            let mut output = Vec::with_capacity(2_000_000);
            let writer = Writer::new(&mut output, OutputFormat::VcfGz);
            let mut writer = writer.write_header(&setup.header).unwrap();
            for i in 1u32..=N_RECORDS {
                seqair_encode_record!(writer, &setup, an, i);
            }
            writer.finish().unwrap();
            black_box(output.len())
        });
    });

    // htslib: writes bgzf-compressed VCF to temp file
    group.bench_function("htslib", |b| {
        b.iter(|| {
            let tmp = tempfile::Builder::new().suffix(".vcf.gz").tempfile().unwrap();
            use rust_htslib::bcf;
            use rust_htslib::bcf::record::GenotypeAllele;
            let hdr = htslib_germline_header();
            let mut writer =
                bcf::Writer::from_path(tmp.path(), &hdr, false, bcf::Format::Vcf).unwrap();
            let rid = writer.header().name2rid(b"chr1").unwrap();
            let pass_id = writer.header().name_to_id(cstr8::cstr8!("PASS")).unwrap();
            let lowqual_id = writer.header().name_to_id(cstr8::cstr8!("LowQual")).unwrap();
            for i in 0..i64::from(N_RECORDS) {
                let allele_idx = i as usize % HTSLIB_ALLELES.len();
                let (ref_bytes, alt1_bytes, alt2_opt) = HTSLIB_ALLELES[allele_idx];
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
                let ac: Vec<i32> = (0..n_alt).map(|k| 1 + ((i + k as i64) % 9) as i32).collect();
                let af: Vec<f32> = ac.iter().map(|&a| a as f32 / an as f32).collect();
                record.push_info_integer(cstr8::cstr8!("DP"), &[dp]).unwrap();
                record.push_info_integer(cstr8::cstr8!("AN"), &[an]).unwrap();
                record.push_info_integer(cstr8::cstr8!("AC"), &ac).unwrap();
                record.push_info_float(cstr8::cstr8!("AF"), &af).unwrap();
                record
                    .push_info_float(cstr8::cstr8!("MQ"), &[55.0f32 + (i % 20) as f32 * 0.5])
                    .unwrap();
                record
                    .push_info_float(cstr8::cstr8!("QD"), &[5.0f32 + (i % 25) as f32 * 0.8])
                    .unwrap();
                record.push_info_float(cstr8::cstr8!("FS"), &[(i % 30) as f32 * 0.5]).unwrap();
                let gts: Vec<GenotypeAllele> = (0..N_GERMLINE_SAMPLES)
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
                let dp_vals: Vec<i32> =
                    (0..N_GERMLINE_SAMPLES as i32).map(|s| dp / 2 + s % 20).collect();
                let gq_vals: Vec<i32> =
                    (0..N_GERMLINE_SAMPLES as i32).map(|s| 30 + (i as i32 + s) % 60).collect();
                record.push_format_integer(cstr8::cstr8!("DP"), &dp_vals).unwrap();
                record.push_format_integer(cstr8::cstr8!("GQ"), &gq_vals).unwrap();
                writer.write(&record).unwrap();
            }
            drop(writer);
            black_box(std::fs::metadata(tmp.path()).unwrap().len())
        });
    });

    group.finish();
}

// ---------------------------------------------------------------------------
// Group 7: Write .bcf (binary), 10 samples
// seqair (typestate encoder, zero-alloc) vs htslib vs noodles
// ---------------------------------------------------------------------------

fn write_bcf(c: &mut Criterion) {
    let mut group = c.benchmark_group("write_bcf");
    group.throughput(Throughput::Elements(u64::from(N_RECORDS)));

    let setup = germline_setup();
    let an = (N_GERMLINE_SAMPLES * 2) as i32;

    // seqair: typestate encoder path (zero-alloc, pre-resolved dict indices)
    group.bench_function("seqair", |b| {
        b.iter(|| {
            let mut output = Vec::with_capacity(2_000_000);
            let writer = Writer::new(&mut output, OutputFormat::Bcf);
            let mut writer = writer.write_header(&setup.header).unwrap();
            for i in 1u32..=N_RECORDS {
                seqair_encode_record!(writer, &setup, an, i);
            }
            writer.finish().unwrap();
            black_box(output.len())
        });
    });

    // htslib
    group.bench_function("htslib", |b| {
        b.iter(|| {
            let tmp = write_htslib_germline(rust_htslib::bcf::Format::Bcf);
            black_box(std::fs::metadata(tmp.path()).unwrap().len())
        });
    });

    // noodles BCF writer
    group.bench_function("noodles", |b| {
        use noodles::bcf;
        use noodles::vcf;
        let noodles_header: vcf::Header = setup.header.to_vcf_text().parse().unwrap();
        b.iter(|| {
            let mut output = Vec::with_capacity(2_000_000);
            let mut writer = bcf::io::Writer::new(&mut output);
            writer.write_header(&noodles_header).unwrap();
            for i in 1u32..=N_RECORDS {
                noodles_write_record!(writer, &noodles_header, an, i);
            }
            drop(writer);
            black_box(output.len())
        });
    });

    group.finish();
}

criterion_group!(benches, write_vcf_text, write_vcf_gz, write_bcf);
criterion_main!(benches);
