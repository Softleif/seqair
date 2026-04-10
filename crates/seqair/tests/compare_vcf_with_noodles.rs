//! Round-trip tests: write BCF/VCF with seqair, read back with noodles, compare.
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

use proptest::prelude::*;
use seqair::vcf::alleles::Alleles;
use seqair::vcf::header::{ContigDef, Number, ValueType};
use seqair::vcf::record::Genotype;
use seqair::vcf::record_encoder::{FormatFieldDef, Gt, InfoFieldDef, Scalar};
use seqair::vcf::{
    ContigId, FormatGt, FormatInt, InfoFlag, InfoFloat, InfoInt, InfoInts, OutputFormat, VcfHeader,
    Writer,
};
use seqair_types::{Base, One, Pos};
use std::io::Cursor;
use std::sync::Arc;

struct RichSetup {
    header: Arc<VcfHeader>,
    contig_chr1: ContigId,
    dp_info: InfoInt,
    gt_fmt: FormatGt,
    dp_fmt: FormatInt,
}

fn rich_setup() -> RichSetup {
    let mut builder = VcfHeader::builder();
    let contig_chr1 =
        builder.register_contig("chr1", ContigDef { length: Some(250_000_000) }).unwrap();
    builder.register_contig("chr2", ContigDef { length: Some(243_000_000) }).unwrap();
    let dp_info = builder
        .register_info(&InfoFieldDef::new(
            "DP",
            Number::Count(1),
            ValueType::Integer,
            "Total Depth",
        ))
        .unwrap();
    let _af_info: InfoFloat = builder
        .register_info(&InfoFieldDef::new(
            "AF",
            Number::AlternateBases,
            ValueType::Float,
            "Allele Frequency",
        ))
        .unwrap();
    let _ad_info: InfoInts = builder
        .register_info(&InfoFieldDef::new(
            "AD",
            Number::ReferenceAlternateBases,
            ValueType::Integer,
            "Allele Depth",
        ))
        .unwrap();
    let _db_flag: InfoFlag = builder
        .register_info(&InfoFieldDef::new(
            "DB",
            Number::Count(0),
            ValueType::Flag,
            "dbSNP membership",
        ))
        .unwrap();
    let gt_fmt = builder
        .register_format(&FormatFieldDef::<Gt>::new(
            "GT",
            Number::Count(1),
            ValueType::String,
            "Genotype",
        ))
        .unwrap();
    let dp_fmt = builder
        .register_format(&FormatFieldDef::<Scalar<i32>>::new(
            "DP",
            Number::Count(1),
            ValueType::Integer,
            "Read Depth",
        ))
        .unwrap();
    let _gq_fmt: FormatInt = builder
        .register_format(&FormatFieldDef::<Scalar<i32>>::new(
            "GQ",
            Number::Count(1),
            ValueType::Integer,
            "Genotype Quality",
        ))
        .unwrap();
    let header = Arc::new(builder.add_sample("sample1").unwrap().build().unwrap());
    RichSetup { header, contig_chr1, dp_info, gt_fmt, dp_fmt }
}

// ── VCF text round-trip via noodles ────────────────────────────────────

#[test]
fn vcf_text_readable_by_noodles() {
    let setup = rich_setup();
    let alleles = Alleles::snv(Base::A, Base::T).unwrap();

    let mut output = Vec::new();
    {
        let writer = Writer::new(&mut output, OutputFormat::Vcf);
        let mut writer = writer.write_header(&setup.header).unwrap();
        let mut enc = writer
            .begin_record(&setup.contig_chr1, Pos::<One>::new(12345).unwrap(), &alleles, Some(30.0))
            .unwrap()
            .filter_pass();
        setup.dp_info.encode(&mut enc, 50);
        let mut enc = enc.begin_samples(1);
        setup.gt_fmt.encode(&mut enc, &[Genotype::unphased(0, 1)]);
        setup.dp_fmt.encode(&mut enc, &[30]);
        enc.emit().unwrap();
        writer.finish().unwrap();
    }

    // Read back with noodles VCF reader
    let mut reader = noodles::vcf::io::Reader::new(Cursor::new(&output));
    let noodles_header = reader.read_header().unwrap();

    let mut records = Vec::new();
    for result in reader.record_bufs(&noodles_header) {
        records.push(result.unwrap());
    }

    assert_eq!(records.len(), 1, "expected 1 record");
    let rec = &records[0];

    // Verify position (noodles uses 1-based)
    assert_eq!(rec.variant_start().map(|p| p.get()), Some(12345usize), "position mismatch");
}

// ── BCF round-trip via noodles ─────────────────────────────────────────

#[test]
fn bcf_write_record_readable_by_noodles() {
    let setup = rich_setup();
    let alleles = Alleles::snv(Base::A, Base::T).unwrap();

    let mut bcf_output = Vec::new();
    {
        let writer = Writer::new(&mut bcf_output, OutputFormat::Bcf);
        let mut writer = writer.write_header(&setup.header).unwrap();
        let mut enc = writer
            .begin_record(&setup.contig_chr1, Pos::<One>::new(100).unwrap(), &alleles, Some(29.5))
            .unwrap()
            .filter_pass();
        setup.dp_info.encode(&mut enc, 42);
        let mut enc = enc.begin_samples(1);
        setup.gt_fmt.encode(&mut enc, &[Genotype::unphased(0, 1)]);
        setup.dp_fmt.encode(&mut enc, &[25]);
        enc.emit().unwrap();
        writer.finish().unwrap();
    }

    // noodles BCF reader
    let mut reader = noodles::bcf::io::Reader::new(Cursor::new(&bcf_output));
    let noodles_header = reader.read_header().unwrap();

    let mut records = Vec::new();
    for result in reader.record_bufs(&noodles_header) {
        records.push(result.unwrap());
    }

    assert_eq!(records.len(), 1, "expected 1 record from BCF");
}

// ── BCF encoder round-trip ─────────────────────────────────────────────

#[test]
fn bcf_encoder_readable_by_noodles() {
    let setup = rich_setup();
    let alleles = Alleles::snv(Base::C, Base::G).unwrap();

    let mut bcf_output = Vec::new();
    {
        let writer = Writer::new(&mut bcf_output, OutputFormat::Bcf);
        let mut writer = writer.write_header(&setup.header).unwrap();
        let mut enc = writer
            .begin_record(&setup.contig_chr1, Pos::<One>::new(500).unwrap(), &alleles, Some(45.0))
            .unwrap()
            .filter_pass();
        setup.dp_info.encode(&mut enc, 100);
        let mut enc = enc.begin_samples(1);
        setup.gt_fmt.encode(&mut enc, &[Genotype::phased_diploid(0, 1)]);
        setup.dp_fmt.encode(&mut enc, &[80]);
        enc.emit().unwrap();
        writer.finish().unwrap();
    }

    // Read back with noodles
    let mut reader = noodles::bcf::io::Reader::new(Cursor::new(&bcf_output));
    let noodles_header = reader.read_header().unwrap();

    let mut records = Vec::new();
    for result in reader.record_bufs(&noodles_header) {
        records.push(result.unwrap());
    }

    assert_eq!(records.len(), 1, "expected 1 record from encoder");
}

// ── Proptest strategies ────────────────────────────────────────────────

fn arb_base() -> impl Strategy<Value = Base> {
    prop_oneof![Just(Base::A), Just(Base::C), Just(Base::G), Just(Base::T),]
}

fn arb_alleles() -> impl Strategy<Value = Alleles> {
    prop_oneof![
        arb_base().prop_map(Alleles::reference),
        (arb_base(), arb_base())
            .prop_filter("ref != alt", |(r, a)| r != a)
            .prop_map(|(r, a)| Alleles::snv(r, a).unwrap()),
        (arb_base(), proptest::collection::vec(arb_base(), 1..4))
            .prop_map(|(anchor, ins)| Alleles::insertion(anchor, &ins).unwrap()),
        (arb_base(), proptest::collection::vec(arb_base(), 1..4))
            .prop_map(|(anchor, del)| Alleles::deletion(anchor, &del).unwrap()),
    ]
}

#[allow(dead_code, reason = "may be used in future test expansions")]
fn arb_genotype() -> impl Strategy<Value = Genotype> {
    prop_oneof![
        (0u16..4, 0u16..4).prop_map(|(a, b)| Genotype::unphased(a, b)),
        (0u16..4, 0u16..4).prop_map(|(a, b)| Genotype::phased_diploid(a, b)),
        (0u16..4).prop_map(Genotype::haploid),
        Just(Genotype::missing_diploid()),
    ]
}

// ── Proptest: BCF via writer round-trips through noodles ───────────────

proptest! {
    #![proptest_config(ProptestConfig::with_cases(50))]

    #[test]
    fn bcf_write_record_roundtrip(
        pos in 1u32..10_000_000,
        alleles in arb_alleles(),
        qual in proptest::option::of(0.1f32..1000.0),
        dp in 0i32..10000,
    ) {
        let setup = rich_setup();
        let pos = Pos::<One>::new(pos).unwrap();

        let mut bcf_output = Vec::new();
        {
            let writer = Writer::new(&mut bcf_output, OutputFormat::Bcf);
            let mut writer = writer.write_header(&setup.header).unwrap();
            let mut enc = writer
                .begin_record(&setup.contig_chr1, pos, &alleles, qual)
                .unwrap()
                .filter_pass();
            setup.dp_info.encode(&mut enc, dp);
            enc.emit().unwrap();
            writer.finish().unwrap();
        }

        // Must be parseable by noodles
        let mut reader = noodles::bcf::io::Reader::new(Cursor::new(&bcf_output));
        let noodles_header = reader.read_header().unwrap();
        let mut records = Vec::new();
        for result in reader.record_bufs(&noodles_header) {
            records.push(result.unwrap());
        }
        prop_assert_eq!(records.len(), 1, "noodles should parse exactly 1 record");
    }

    #[test]
    fn vcf_text_roundtrip_through_noodles(
        pos in 1u32..10_000_000,
        alleles in arb_alleles(),
        dp in 0i32..10000,
    ) {
        // Sites-only header (no samples) to avoid FORMAT column requirement
        let mut builder = VcfHeader::builder();
        let contig =
            builder.register_contig("chr1", ContigDef { length: Some(250_000_000) }).unwrap();
        let dp_info: InfoInt = builder
            .register_info(&InfoFieldDef::new(
                "DP",
                Number::Count(1),
                ValueType::Integer,
                "Depth",
            ))
            .unwrap();
        let header = Arc::new(builder.build().unwrap());

        let pos = Pos::<One>::new(pos).unwrap();

        let mut output = Vec::new();
        {
            let writer = Writer::new(&mut output, OutputFormat::Vcf);
            let mut writer = writer.write_header(&header).unwrap();
            let mut enc = writer
                .begin_record(&contig, pos, &alleles, None)
                .unwrap()
                .filter_pass();
            dp_info.encode(&mut enc, dp);
            enc.emit().unwrap();
            writer.finish().unwrap();
        }

        let mut reader = noodles::vcf::io::Reader::new(Cursor::new(&output));
        let noodles_header = reader.read_header().unwrap();
        let mut records = Vec::new();
        for result in reader.record_bufs(&noodles_header) {
            records.push(result.unwrap());
        }
        prop_assert_eq!(records.len(), 1, "noodles VCF reader should parse exactly 1 record");
    }
}
