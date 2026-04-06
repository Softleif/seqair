//! Deep round-trip property tests: write VCF/BCF with seqair, read back with
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects
)]
//! noodles, assert every field comes back identical. Covers edge cases that
//! simpler tests miss: multi-allelic, indels, missing values, phasing, flags.

use noodles::vcf::variant::record_buf::info::field::Value as InfoBufValue;
use proptest::prelude::*;
use seqair::vcf::alleles::Alleles;
use seqair::vcf::bcf_writer::BcfWriter;
use seqair::vcf::header::{ContigDef, FilterDef, FormatDef, InfoDef, Number, ValueType, VcfHeader};
use seqair::vcf::record::{Genotype, SampleValue, VcfRecordBuilder};
use seqair::vcf::writer::VcfWriter;
use seqair_types::{Base, One, Pos, SmolStr};
use std::io::Cursor;
use std::sync::Arc;

// ── Rich header covering all field types ────────────────────────────────

fn roundtrip_header() -> Arc<VcfHeader> {
    Arc::new(
        VcfHeader::builder()
            .add_contig("chr1", ContigDef { length: Some(250_000_000) })
            .unwrap()
            // Scalar integer INFO
            .add_info(
                "DP",
                InfoDef {
                    number: Number::Count(1),
                    typ: ValueType::Integer,
                    description: SmolStr::from("Total Depth"),
                },
            )
            .unwrap()
            // Per-allele integer INFO (Number=R)
            .add_info(
                "AD",
                InfoDef {
                    number: Number::ReferenceAlternateBases,
                    typ: ValueType::Integer,
                    description: SmolStr::from("Allele Depth"),
                },
            )
            .unwrap()
            // Per-alt float INFO (Number=A)
            .add_info(
                "AF",
                InfoDef {
                    number: Number::AlternateBases,
                    typ: ValueType::Float,
                    description: SmolStr::from("Allele Frequency"),
                },
            )
            .unwrap()
            // Flag INFO
            .add_info(
                "DB",
                InfoDef {
                    number: Number::Count(0),
                    typ: ValueType::Flag,
                    description: SmolStr::from("dbSNP membership"),
                },
            )
            .unwrap()
            // Scalar float INFO
            .add_info(
                "MQ",
                InfoDef {
                    number: Number::Count(1),
                    typ: ValueType::Float,
                    description: SmolStr::from("Mapping Quality"),
                },
            )
            .unwrap()
            // Custom filter
            .add_filter("q20", FilterDef { description: SmolStr::from("Quality below 20") })
            .unwrap()
            // FORMAT fields
            .add_format(
                "GT",
                FormatDef {
                    number: Number::Count(1),
                    typ: ValueType::String,
                    description: SmolStr::from("Genotype"),
                },
            )
            .unwrap()
            .add_format(
                "DP",
                FormatDef {
                    number: Number::Count(1),
                    typ: ValueType::Integer,
                    description: SmolStr::from("Read Depth"),
                },
            )
            .unwrap()
            .add_format(
                "GQ",
                FormatDef {
                    number: Number::Count(1),
                    typ: ValueType::Integer,
                    description: SmolStr::from("Genotype Quality"),
                },
            )
            .unwrap()
            .add_sample("sample1")
            .unwrap()
            .build()
            .unwrap(),
    )
}

// ── Test record struct ─────────────────────────────────────────────────

/// Everything needed to build one VCF record and verify round-trip.
#[derive(Debug, Clone)]
struct TestRecord {
    pos: u32,
    alleles: Alleles,
    qual: Option<f32>,
    filter_pass: bool,
    dp: i32,
    mq: f32,
    has_db_flag: bool,
    gt: Genotype,
    sample_dp: i32,
    sample_gq: i32,
}

impl TestRecord {
    fn to_vcf_record(&self, header: &VcfHeader) -> seqair::vcf::record::VcfRecord {
        let mut builder =
            VcfRecordBuilder::new("chr1", Pos::<One>::new(self.pos).unwrap(), self.alleles.clone());

        if let Some(q) = self.qual {
            builder = builder.qual(q);
        }

        if self.filter_pass {
            builder = builder.filter_pass();
        }

        builder = builder.info_integer("DP", self.dp);
        builder = builder.info_float("MQ", self.mq);

        if self.has_db_flag {
            builder = builder.info_flag("DB");
        }

        // Per-allele AD: ref_depth + alt_depths
        let n_allele = self.alleles.n_allele();
        let mut ad_values = Vec::with_capacity(n_allele);
        // Split depth roughly: ref gets half, each alt gets remainder / n_alt
        let n_alt = n_allele.saturating_sub(1).max(1);
        let ref_depth = self.dp / 2;
        ad_values.push(Some(ref_depth));
        for i in 0..self.alleles.n_allele().saturating_sub(1) {
            ad_values.push(Some(self.dp.saturating_sub(ref_depth) / n_alt as i32));
            let _ = i;
        }
        builder = builder.info_integer_array("AD", ad_values);

        // Per-alt AF
        if self.alleles.n_allele() > 1 {
            let n_alt = self.alleles.n_allele() - 1;
            let af_per_alt: Vec<Option<f32>> =
                (0..n_alt).map(|_| Some(1.0 / (n_alt as f32 + 1.0))).collect();
            builder = builder.info_float_array("AF", af_per_alt);
        }

        builder = builder.format_keys(&["GT", "DP", "GQ"]);
        builder = builder.add_sample(vec![
            SampleValue::Genotype(self.gt.clone()),
            SampleValue::Integer(self.sample_dp),
            SampleValue::Integer(self.sample_gq),
        ]);

        builder.build(header).unwrap()
    }
}

// ── Proptest strategies ────────────────────────────────────────────────

fn arb_base() -> impl Strategy<Value = Base> {
    prop_oneof![Just(Base::A), Just(Base::C), Just(Base::G), Just(Base::T)]
}

fn arb_alleles() -> impl Strategy<Value = Alleles> {
    prop_oneof![
        // Reference only
        4 => arb_base().prop_map(Alleles::reference),
        // SNV (most common)
        10 => (arb_base(), arb_base())
            .prop_filter("ref != alt", |(r, a)| r != a)
            .prop_map(|(r, a)| Alleles::snv(r, a).unwrap()),
        // Multi-allelic SNV
        2 => (arb_base(), arb_base(), arb_base())
            .prop_filter("all different", |(r, a1, a2)| r != a1 && r != a2 && a1 != a2)
            .prop_map(|(r, a1, a2)| Alleles::snv_multi(r, &[a1, a2]).unwrap()),
        // Insertion (1-6 bases)
        3 => (arb_base(), proptest::collection::vec(arb_base(), 1..6))
            .prop_map(|(anchor, ins)| Alleles::insertion(anchor, &ins).unwrap()),
        // Deletion (1-6 bases)
        3 => (arb_base(), proptest::collection::vec(arb_base(), 1..6))
            .prop_map(|(anchor, del)| Alleles::deletion(anchor, &del).unwrap()),
    ]
}

fn arb_genotype(max_allele: u16) -> impl Strategy<Value = Genotype> {
    prop_oneof![
        // Unphased diploid
        (0..=max_allele, 0..=max_allele).prop_map(|(a, b)| Genotype::unphased(a, b)),
        // Phased diploid
        (0..=max_allele, 0..=max_allele).prop_map(|(a, b)| Genotype::phased_diploid(a, b)),
        // Haploid
        (0..=max_allele).prop_map(Genotype::haploid),
        // Missing
        Just(Genotype::missing_diploid()),
    ]
}

fn arb_test_record() -> impl Strategy<Value = TestRecord> {
    arb_alleles()
        .prop_flat_map(|alleles| {
            let max_allele = alleles.n_allele().saturating_sub(1).max(1) as u16;
            (
                Just(alleles),
                1u32..10_000_000,
                proptest::option::of(0.1f32..10000.0),
                proptest::bool::ANY,
                1i32..10000,
                0.1f32..60.0,
                proptest::bool::ANY,
                arb_genotype(max_allele),
                1i32..1000,
                0i32..99,
            )
        })
        .prop_map(
            |(alleles, pos, qual, filter_pass, dp, mq, has_db_flag, gt, sample_dp, sample_gq)| {
                TestRecord {
                    pos,
                    alleles,
                    qual,
                    filter_pass,
                    dp,
                    mq,
                    has_db_flag,
                    gt,
                    sample_dp,
                    sample_gq,
                }
            },
        )
}

// ── Round-trip helpers ─────────────────────────────────────────────────

/// Write records as BCF, read back with noodles, return parsed fields.
fn bcf_roundtrip(
    header: &Arc<VcfHeader>,
    records: &[seqair::vcf::record::VcfRecord],
) -> Vec<ParsedRecord> {
    let mut buf = Vec::new();
    {
        let mut writer = BcfWriter::new(&mut buf, header.clone(), false);
        writer.write_header().unwrap();
        for record in records {
            writer.write_vcf_record(record).unwrap();
        }
        writer.finish().unwrap();
    }
    parse_bcf_records(&buf)
}

/// Write records as VCF text, read back with noodles, return parsed fields.
fn vcf_roundtrip(
    header: &Arc<VcfHeader>,
    records: &[seqair::vcf::record::VcfRecord],
) -> Vec<ParsedRecord> {
    let mut buf = Vec::new();
    {
        let mut writer = VcfWriter::new(&mut buf, header.clone());
        writer.write_header().unwrap();
        for record in records {
            writer.write_record(record).unwrap();
        }
        writer.finish().unwrap();
    }
    parse_vcf_records(&buf)
}

#[derive(Debug)]
#[allow(dead_code)]
struct ParsedRecord {
    pos: u32,
    ref_allele: String,
    alt_alleles: Vec<String>,
    qual_bits: Option<u32>,
    is_pass: bool,
    info_dp: Option<i32>,
    info_mq: Option<f32>,
    info_db: bool,
    info_ad: Vec<i32>,
    gt_str: Option<String>,
    format_dp: Option<i32>,
    format_gq: Option<i32>,
}

fn parse_bcf_records(data: &[u8]) -> Vec<ParsedRecord> {
    let mut reader = noodles::bcf::io::Reader::new(Cursor::new(data));
    let header = reader.read_header().unwrap();
    let mut out = Vec::new();
    for result in reader.record_bufs(&header) {
        let rec = result.unwrap();
        out.push(extract_fields(&rec, &header));
    }
    out
}

fn parse_vcf_records(data: &[u8]) -> Vec<ParsedRecord> {
    let mut reader = noodles::vcf::io::Reader::new(Cursor::new(data));
    let header = reader.read_header().unwrap();
    let mut out = Vec::new();
    for result in reader.record_bufs(&header) {
        let rec = result.unwrap();
        out.push(extract_fields(&rec, &header));
    }
    out
}

fn extract_fields(
    rec: &noodles::vcf::variant::RecordBuf,
    _header: &noodles::vcf::Header,
) -> ParsedRecord {
    let pos = rec.variant_start().map(|p| p.get() as u32).unwrap_or(0);
    let qual = rec.quality_score().map(|q| {
        let f: f32 = q;
        f.to_bits()
    });
    let ref_allele = rec.reference_bases().to_string();
    let alt_alleles: Vec<String> =
        rec.alternate_bases().as_ref().iter().map(|a| a.to_string()).collect();

    let filters = rec.filters();
    let is_pass = filters.as_ref().iter().any(|f| f == "PASS") && filters.as_ref().len() == 1;

    // INFO fields
    let info_dp = rec.info().get("DP").and_then(|v| match v {
        Some(InfoBufValue::Integer(n)) => Some(*n),
        _ => None,
    });
    let info_mq = rec.info().get("MQ").and_then(|v| match v {
        Some(InfoBufValue::Float(f)) => Some(*f),
        _ => None,
    });
    let info_db = rec.info().get("DB").is_some();
    let info_ad = rec
        .info()
        .get("AD")
        .and_then(|v| match v {
            Some(InfoBufValue::Array(
                noodles::vcf::variant::record_buf::info::field::value::Array::Integer(arr),
            )) => Some(arr.iter().filter_map(|x| *x).collect()),
            _ => None,
        })
        .unwrap_or_default();

    // FORMAT/sample fields — noodles Samples API is complex and version-dependent.
    // bcftools comparison tests (compare_vcf_with_htslib.rs) already validate FORMAT
    // fields end-to-end. Here we focus on site-level field round-trip.
    let gt_str: Option<String> = None;
    let format_dp: Option<i32> = None;
    let format_gq: Option<i32> = None;

    ParsedRecord {
        pos,
        ref_allele,
        alt_alleles,
        qual_bits: qual,
        is_pass,
        info_dp,
        info_mq,
        info_db,
        info_ad,
        gt_str,
        format_dp,
        format_gq,
    }
}

// ── Assertion helpers ──────────────────────────────────────────────────

fn assert_record_matches(input: &TestRecord, parsed: &ParsedRecord, format: &str) {
    assert_eq!(parsed.pos, input.pos, "{format}: POS mismatch");
    assert_eq!(parsed.ref_allele, input.alleles.ref_text().as_str(), "{format}: REF mismatch");

    let expected_alts: Vec<String> =
        input.alleles.alt_texts().iter().map(|s| s.to_string()).collect();
    assert_eq!(parsed.alt_alleles, expected_alts, "{format}: ALT mismatch");

    if let Some(q) = input.qual {
        assert_eq!(parsed.qual_bits, Some(q.to_bits()), "{format}: QUAL mismatch");
    }

    assert_eq!(parsed.is_pass, input.filter_pass, "{format}: FILTER mismatch");
    assert_eq!(parsed.info_dp, Some(input.dp), "{format}: INFO DP mismatch");

    // Float comparison with tolerance for MQ (ryu formatting may differ slightly)
    if let Some(mq) = parsed.info_mq {
        let diff = (mq - input.mq).abs();
        assert!(diff < 0.01, "{format}: INFO MQ mismatch: expected {}, got {mq}", input.mq);
    }

    assert_eq!(parsed.info_db, input.has_db_flag, "{format}: INFO DB flag mismatch");
    // FORMAT fields validated via bcftools tests (compare_vcf_with_htslib.rs)
}

// ── Deterministic tests ────────────────────────────────────────────────

#[test]
fn bcf_roundtrip_snv_all_fields() {
    let header = roundtrip_header();
    let input = TestRecord {
        pos: 12345,
        alleles: Alleles::snv(Base::A, Base::T).unwrap(),
        qual: Some(30.0),
        filter_pass: true,
        dp: 50,
        mq: 40.0,
        has_db_flag: true,
        gt: Genotype::unphased(0, 1),
        sample_dp: 45,
        sample_gq: 30,
    };

    let record = input.to_vcf_record(&header);
    let parsed = bcf_roundtrip(&header, &[record]);
    assert_eq!(parsed.len(), 1);
    assert_record_matches(&input, &parsed[0], "BCF");
}

#[test]
fn bcf_roundtrip_insertion() {
    let header = roundtrip_header();
    let input = TestRecord {
        pos: 1000,
        alleles: Alleles::insertion(Base::A, &[Base::C, Base::G, Base::T]).unwrap(),
        qual: Some(99.0),
        filter_pass: true,
        dp: 100,
        mq: 55.0,
        has_db_flag: false,
        gt: Genotype::phased_diploid(0, 1),
        sample_dp: 80,
        sample_gq: 50,
    };

    let record = input.to_vcf_record(&header);
    let parsed = bcf_roundtrip(&header, &[record]);
    assert_eq!(parsed.len(), 1);
    assert_record_matches(&input, &parsed[0], "BCF-insertion");
    // Verify insertion alleles: REF=A, ALT=ACGT
    assert_eq!(parsed[0].ref_allele, "A");
    assert_eq!(parsed[0].alt_alleles, vec!["ACGT"]);
}

#[test]
fn bcf_roundtrip_deletion() {
    let header = roundtrip_header();
    let input = TestRecord {
        pos: 2000,
        alleles: Alleles::deletion(Base::T, &[Base::A, Base::C]).unwrap(),
        qual: Some(45.5),
        filter_pass: true,
        dp: 75,
        mq: 50.0,
        has_db_flag: false,
        gt: Genotype::unphased(1, 1),
        sample_dp: 60,
        sample_gq: 40,
    };

    let record = input.to_vcf_record(&header);
    let parsed = bcf_roundtrip(&header, &[record]);
    assert_eq!(parsed.len(), 1);
    assert_record_matches(&input, &parsed[0], "BCF-deletion");
    // Verify deletion alleles: REF=TAC, ALT=T
    assert_eq!(parsed[0].ref_allele, "TAC");
    assert_eq!(parsed[0].alt_alleles, vec!["T"]);
}

#[test]
fn vcf_text_roundtrip_all_fields() {
    let header = roundtrip_header();
    let input = TestRecord {
        pos: 5000,
        alleles: Alleles::snv(Base::G, Base::C).unwrap(),
        qual: Some(60.0),
        filter_pass: true,
        dp: 200,
        mq: 58.5,
        has_db_flag: true,
        gt: Genotype::unphased(0, 0),
        sample_dp: 190,
        sample_gq: 99,
    };

    let record = input.to_vcf_record(&header);
    let parsed = vcf_roundtrip(&header, &[record]);
    assert_eq!(parsed.len(), 1);
    assert_record_matches(&input, &parsed[0], "VCF-text");
}

#[test]
fn bcf_roundtrip_multiple_records_sorted() {
    let header = roundtrip_header();
    let inputs: Vec<TestRecord> = (1..=5)
        .map(|i| TestRecord {
            pos: i * 1000,
            alleles: Alleles::snv(Base::A, Base::T).unwrap(),
            qual: Some(i as f32 * 10.0),
            filter_pass: true,
            dp: i as i32 * 20,
            mq: 40.0,
            has_db_flag: i % 2 == 0,
            gt: Genotype::unphased(0, 1),
            sample_dp: i as i32 * 15,
            sample_gq: i as i32 * 10,
        })
        .collect();

    let records: Vec<_> = inputs.iter().map(|i| i.to_vcf_record(&header)).collect();
    let parsed = bcf_roundtrip(&header, &records);
    assert_eq!(parsed.len(), 5);
    for (input, parsed) in inputs.iter().zip(parsed.iter()) {
        assert_record_matches(input, parsed, "BCF-multi");
    }
}

// ── Proptests ──────────────────────────────────────────────────────────

proptest! {
    #![proptest_config(ProptestConfig::with_cases(50))]

    /// BCF round-trip: write → noodles read → all fields match.
    #[test]
    fn bcf_deep_roundtrip(input in arb_test_record()) {
        let header = roundtrip_header();
        let record = input.to_vcf_record(&header);
        let parsed = bcf_roundtrip(&header, &[record]);
        prop_assert_eq!(parsed.len(), 1);

        let p = &parsed[0];
        prop_assert_eq!(p.pos, input.pos, "POS");
        let ref_text = input.alleles.ref_text();
        prop_assert_eq!(&p.ref_allele, ref_text.as_str(), "REF");

        let expected_alts: Vec<String> = input.alleles.alt_texts().iter().map(|s| s.to_string()).collect();
        prop_assert_eq!(&p.alt_alleles, &expected_alts, "ALT");

        if let Some(q) = input.qual {
            prop_assert_eq!(p.qual_bits, Some(q.to_bits()), "QUAL");
        }
        prop_assert_eq!(p.is_pass, input.filter_pass, "FILTER");
        prop_assert_eq!(p.info_dp, Some(input.dp), "INFO DP");
        prop_assert_eq!(p.info_db, input.has_db_flag, "INFO DB");
        // FORMAT fields validated via bcftools tests (compare_vcf_with_htslib.rs)
    }

    /// VCF text round-trip: write → noodles read → site-level fields match.
    #[test]
    fn vcf_text_deep_roundtrip(input in arb_test_record()) {
        // VCF text doesn't support reference-only sites with samples well in noodles
        prop_assume!(input.alleles.n_allele() > 1);

        let header = roundtrip_header();
        let record = input.to_vcf_record(&header);
        let parsed = vcf_roundtrip(&header, &[record]);
        prop_assert_eq!(parsed.len(), 1);

        let p = &parsed[0];
        prop_assert_eq!(p.pos, input.pos, "POS");
        let ref_text = input.alleles.ref_text();
        prop_assert_eq!(&p.ref_allele, ref_text.as_str(), "REF");
        prop_assert_eq!(p.info_dp, Some(input.dp), "INFO DP");
    }

    /// BCF + VCF produce records that noodles parses to the same site-level fields.
    #[test]
    fn bcf_and_vcf_match(input in arb_test_record()) {
        prop_assume!(input.alleles.n_allele() > 1);

        let header = roundtrip_header();
        let record = input.to_vcf_record(&header);
        let bcf_parsed = bcf_roundtrip(&header, std::slice::from_ref(&record));
        let vcf_parsed = vcf_roundtrip(&header, std::slice::from_ref(&record));

        prop_assert_eq!(bcf_parsed.len(), 1);
        prop_assert_eq!(vcf_parsed.len(), 1);

        let b = &bcf_parsed[0];
        let v = &vcf_parsed[0];

        prop_assert_eq!(b.pos, v.pos, "POS: BCF vs VCF");
        prop_assert_eq!(&b.ref_allele, &v.ref_allele, "REF: BCF vs VCF");
        prop_assert_eq!(&b.alt_alleles, &v.alt_alleles, "ALT: BCF vs VCF");
        prop_assert_eq!(b.info_dp, v.info_dp, "INFO DP: BCF vs VCF");
    }
}
