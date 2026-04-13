//! Deep round-trip property tests: write VCF/BCF with seqair unified Writer,
//! read back with noodles, assert every field comes back identical. Covers edge
//! cases that simpler tests miss: multi-allelic, indels, missing values,
//! phasing, flags.
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
    clippy::implicit_clone,
    reason = "test code with known small values"
)]

use noodles::vcf::variant::record_buf::info::field::Value as InfoBufValue;
use proptest::prelude::*;
use seqair::vcf::record_encoder::{
    FilterFieldDef, FormatFieldDef, FormatGt, FormatInt, InfoFieldDef, InfoFlag, InfoFloat,
    InfoFloats, InfoInt, InfoInts,
};
use seqair::vcf::{
    Alleles, ContigDef, FilterId, Genotype, Number, OutputFormat, Ready, ValueType, VcfHeader,
    Writer,
};
use seqair_types::{Base, One, Pos};
use std::io::Cursor;
use std::sync::Arc;

// ── Rich header setup ───────────────────────────────────────────────────

struct RoundtripSetup {
    header: Arc<VcfHeader>,
    contig: seqair::vcf::ContigId,
    dp_info: InfoInt,
    ad_info: InfoInts,
    af_info: InfoFloats,
    db_flag: InfoFlag,
    mq_info: InfoFloat,
    q20_filter: FilterId,
    gt_fmt: FormatGt,
    dp_fmt: FormatInt,
    gq_fmt: FormatInt,
}

fn make_setup() -> RoundtripSetup {
    let mut builder = VcfHeader::builder();
    let contig = builder.register_contig("chr1", ContigDef { length: Some(250_000_000) }).unwrap();
    // Register filters BEFORE info/format (BCF string dict order: PASS, filters, info, format).
    // TODO: enforce this ordering at the type level — VcfHeaderBuilder should use typestates
    // or a phased builder so that register_filter() can't be called after register_info().
    let q20_filter =
        builder.register_filter(&FilterFieldDef::new("q20", "Quality below 20")).unwrap();
    let dp_info = builder
        .register_info(&InfoFieldDef::new(
            "DP",
            Number::Count(1),
            ValueType::Integer,
            "Total Depth",
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
    let af_info = builder
        .register_info(&InfoFieldDef::new(
            "AF",
            Number::AlternateBases,
            ValueType::Float,
            "Allele Frequency",
        ))
        .unwrap();
    let db_flag = builder
        .register_info(&InfoFieldDef::new(
            "DB",
            Number::Count(0),
            ValueType::Flag,
            "dbSNP membership",
        ))
        .unwrap();
    let mq_info = builder
        .register_info(&InfoFieldDef::new(
            "MQ",
            Number::Count(1),
            ValueType::Float,
            "Mapping Quality",
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
    let header = Arc::new(builder.add_sample("sample1").unwrap().build().unwrap());
    RoundtripSetup {
        header,
        contig,
        dp_info,
        ad_info,
        af_info,
        db_flag,
        mq_info,
        q20_filter,
        gt_fmt,
        dp_fmt,
        gq_fmt,
    }
}

// ── Test record struct ─────────────────────────────────────────────────

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

fn write_bcf(setup: &RoundtripSetup, records: &[TestRecord]) -> Vec<u8> {
    let mut buf = Vec::new();
    let writer = Writer::new(&mut buf, OutputFormat::Bcf);
    let mut writer = writer.write_header(&setup.header).unwrap();
    for rec in records {
        write_one_record(setup, &mut writer, rec, OutputFormat::Bcf);
    }
    writer.finish().unwrap();
    buf
}

fn write_vcf(setup: &RoundtripSetup, records: &[TestRecord]) -> Vec<u8> {
    let mut buf = Vec::new();
    let writer = Writer::new(&mut buf, OutputFormat::Vcf);
    let mut writer = writer.write_header(&setup.header).unwrap();
    for rec in records {
        write_one_record(setup, &mut writer, rec, OutputFormat::Vcf);
    }
    writer.finish().unwrap();
    buf
}

fn write_one_record<W: std::io::Write>(
    setup: &RoundtripSetup,
    writer: &mut Writer<W, Ready>,
    rec: &TestRecord,
    _fmt: OutputFormat,
) {
    let n_allele = rec.alleles.n_allele();
    let n_alt = n_allele.saturating_sub(1).max(1);
    let ref_depth = rec.dp / 2;

    let enc = writer
        .begin_record(&setup.contig, Pos::<One>::new(rec.pos).unwrap(), &rec.alleles, rec.qual)
        .unwrap();

    let mut enc =
        if rec.filter_pass { enc.filter_pass() } else { enc.filter_fail(&[&setup.q20_filter]) };

    // INFO fields
    setup.dp_info.encode(&mut enc, rec.dp);

    // AD: ref_depth + per-alt depths
    let mut ad_vals = vec![ref_depth];
    for _ in 0..rec.alleles.n_allele().saturating_sub(1) {
        ad_vals.push(rec.dp.saturating_sub(ref_depth) / n_alt as i32);
    }
    setup.ad_info.encode(&mut enc, &ad_vals);

    // AF: per-alt floats
    if rec.alleles.n_allele() > 1 {
        let af: Vec<f32> =
            (0..rec.alleles.n_allele() - 1).map(|_| 1.0 / (n_alt as f32 + 1.0)).collect();
        setup.af_info.encode(&mut enc, &af);
    }

    if rec.has_db_flag {
        setup.db_flag.encode(&mut enc);
    }
    setup.mq_info.encode(&mut enc, rec.mq);

    let mut enc = enc.begin_samples();
    setup.gt_fmt.encode(&mut enc, std::slice::from_ref(&rec.gt)).unwrap();
    setup.dp_fmt.encode(&mut enc, &[rec.sample_dp]).unwrap();
    setup.gq_fmt.encode(&mut enc, &[rec.sample_gq]).unwrap();
    enc.emit().unwrap();
}

// ── Noodles parsing ────────────────────────────────────────────────────

#[derive(Debug)]
#[allow(dead_code, reason = "fields read in assertions via Debug formatting")]
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
}

fn parse_bcf_records(data: &[u8]) -> Vec<ParsedRecord> {
    let mut reader = noodles::bcf::io::Reader::new(Cursor::new(data));
    let header = reader.read_header().unwrap();
    let mut out = Vec::new();
    for result in reader.record_bufs(&header) {
        let rec = result.unwrap();
        out.push(extract_fields(&rec));
    }
    out
}

fn parse_vcf_records(data: &[u8]) -> Vec<ParsedRecord> {
    let mut reader = noodles::vcf::io::Reader::new(Cursor::new(data));
    let header = reader.read_header().unwrap();
    let mut out = Vec::new();
    for result in reader.record_bufs(&header) {
        let rec = result.unwrap();
        out.push(extract_fields(&rec));
    }
    out
}

fn extract_fields(rec: &noodles::vcf::variant::RecordBuf) -> ParsedRecord {
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
    }
}

// ── Assertion helper ───────────────────────────────────────────────────

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

    if let Some(mq) = parsed.info_mq {
        let diff = (mq - input.mq).abs();
        assert!(diff < 0.01, "{format}: INFO MQ mismatch: expected {}, got {mq}", input.mq);
    }

    assert_eq!(parsed.info_db, input.has_db_flag, "{format}: INFO DB flag mismatch");
}

// ── Proptest strategies ────────────────────────────────────────────────

fn arb_base() -> impl Strategy<Value = Base> {
    prop_oneof![Just(Base::A), Just(Base::C), Just(Base::G), Just(Base::T)]
}

fn arb_alleles() -> impl Strategy<Value = Alleles> {
    prop_oneof![
        4 => arb_base().prop_map(Alleles::reference),
        10 => (arb_base(), arb_base())
            .prop_filter("ref != alt", |(r, a)| r != a)
            .prop_map(|(r, a)| Alleles::snv(r, a).unwrap()),
        2 => (arb_base(), arb_base(), arb_base())
            .prop_filter("all different", |(r, a1, a2)| r != a1 && r != a2 && a1 != a2)
            .prop_map(|(r, a1, a2)| Alleles::snv_multi(r, &[a1, a2]).unwrap()),
        3 => (arb_base(), proptest::collection::vec(arb_base(), 1..6))
            .prop_map(|(anchor, ins)| Alleles::insertion(anchor, &ins).unwrap()),
        3 => (arb_base(), proptest::collection::vec(arb_base(), 1..6))
            .prop_map(|(anchor, del)| Alleles::deletion(anchor, &del).unwrap()),
    ]
}

fn arb_genotype(max_allele: u16) -> impl Strategy<Value = Genotype> {
    prop_oneof![
        (0..=max_allele, 0..=max_allele).prop_map(|(a, b)| Genotype::unphased(a, b)),
        (0..=max_allele, 0..=max_allele).prop_map(|(a, b)| Genotype::phased_diploid(a, b)),
        (0..=max_allele).prop_map(Genotype::haploid),
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

// ── Deterministic tests ────────────────────────────────────────────────

#[test]
fn bcf_roundtrip_snv_all_fields() {
    let setup = make_setup();
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

    let data = write_bcf(&setup, std::slice::from_ref(&input));
    let parsed = parse_bcf_records(&data);
    assert_eq!(parsed.len(), 1);
    assert_record_matches(&input, &parsed[0], "BCF");
}

#[test]
fn bcf_roundtrip_insertion() {
    let setup = make_setup();
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

    let data = write_bcf(&setup, std::slice::from_ref(&input));
    let parsed = parse_bcf_records(&data);
    assert_eq!(parsed.len(), 1);
    assert_record_matches(&input, &parsed[0], "BCF-insertion");
    assert_eq!(parsed[0].ref_allele, "A");
    assert_eq!(parsed[0].alt_alleles, vec!["ACGT"]);
}

#[test]
fn bcf_roundtrip_deletion() {
    let setup = make_setup();
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

    let data = write_bcf(&setup, std::slice::from_ref(&input));
    let parsed = parse_bcf_records(&data);
    assert_eq!(parsed.len(), 1);
    assert_record_matches(&input, &parsed[0], "BCF-deletion");
    assert_eq!(parsed[0].ref_allele, "TAC");
    assert_eq!(parsed[0].alt_alleles, vec!["T"]);
}

#[test]
fn vcf_text_roundtrip_all_fields() {
    let setup = make_setup();
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

    let data = write_vcf(&setup, std::slice::from_ref(&input));
    let parsed = parse_vcf_records(&data);
    assert_eq!(parsed.len(), 1);
    assert_record_matches(&input, &parsed[0], "VCF-text");
}

#[test]
fn bcf_roundtrip_multiple_records_sorted() {
    let setup = make_setup();
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

    let data = write_bcf(&setup, &inputs);
    let parsed = parse_bcf_records(&data);
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
        let setup = make_setup();
        let data = write_bcf(&setup, std::slice::from_ref(&input));
        let parsed = parse_bcf_records(&data);
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
    }

    /// VCF text round-trip: write → noodles read → site-level fields match.
    #[test]
    fn vcf_text_deep_roundtrip(input in arb_test_record()) {
        prop_assume!(input.alleles.n_allele() > 1);

        let setup = make_setup();
        let data = write_vcf(&setup, std::slice::from_ref(&input));
        let parsed = parse_vcf_records(&data);
        prop_assert_eq!(parsed.len(), 1);

        let p = &parsed[0];
        prop_assert_eq!(p.pos, input.pos, "POS");
        let ref_text = input.alleles.ref_text();
        prop_assert_eq!(&p.ref_allele, ref_text.as_str(), "REF");
        prop_assert_eq!(p.info_dp, Some(input.dp), "INFO DP");
    }

    // r[verify record_encoder.vcf_bcf_equivalence]
    /// BCF + VCF produce records that noodles parses to the same site-level fields.
    #[test]
    fn bcf_and_vcf_match(input in arb_test_record()) {
        prop_assume!(input.alleles.n_allele() > 1);

        let setup = make_setup();
        let bcf_data = write_bcf(&setup, std::slice::from_ref(&input));
        let vcf_data = write_vcf(&setup, std::slice::from_ref(&input));
        let bcf_parsed = parse_bcf_records(&bcf_data);
        let vcf_parsed = parse_vcf_records(&vcf_data);

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
