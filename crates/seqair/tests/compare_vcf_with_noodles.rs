//! Round-trip tests: write BCF/VCF with seqair, read back with noodles, compare.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects
)]

use proptest::prelude::*;
use seqair::vcf::alleles::Alleles;
use seqair::vcf::bcf_writer::BcfWriter;
use seqair::vcf::encoder::{
    ContigHandle, FilterHandle, GtFormatHandle, ScalarFormatHandle, ScalarInfoHandle,
};
use seqair::vcf::header::{ContigDef, FormatDef, InfoDef, Number, ValueType, VcfHeader};
use seqair::vcf::record::{Genotype, SampleValue, VcfRecordBuilder};
use seqair::vcf::writer::VcfWriter;
use seqair_types::{Base, One, Pos, SmolStr};
use std::io::Cursor;
use std::marker::PhantomData;
use std::sync::Arc;

fn rich_header() -> Arc<VcfHeader> {
    Arc::new(
        VcfHeader::builder()
            .add_contig("chr1", ContigDef { length: Some(250_000_000) })
            .unwrap()
            .add_contig("chr2", ContigDef { length: Some(243_000_000) })
            .unwrap()
            .add_info(
                "DP",
                InfoDef {
                    number: Number::Count(1),
                    typ: ValueType::Integer,
                    description: SmolStr::from("Total Depth"),
                },
            )
            .unwrap()
            .add_info(
                "AF",
                InfoDef {
                    number: Number::AlternateBases,
                    typ: ValueType::Float,
                    description: SmolStr::from("Allele Frequency"),
                },
            )
            .unwrap()
            .add_info(
                "AD",
                InfoDef {
                    number: Number::ReferenceAlternateBases,
                    typ: ValueType::Integer,
                    description: SmolStr::from("Allele Depth"),
                },
            )
            .unwrap()
            .add_info(
                "DB",
                InfoDef {
                    number: Number::Count(0),
                    typ: ValueType::Flag,
                    description: SmolStr::from("dbSNP membership"),
                },
            )
            .unwrap()
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

// ── VCF text round-trip via noodles ────────────────────────────────────

#[test]
fn vcf_text_readable_by_noodles() {
    let header = rich_header();
    let record = VcfRecordBuilder::new(
        "chr1",
        Pos::<One>::new(12345).unwrap(),
        Alleles::snv(Base::A, Base::T).unwrap(),
    )
    .qual(30.0)
    .filter_pass()
    .info_integer("DP", 50)
    .format_keys(&["GT", "DP"])
    .add_sample(vec![SampleValue::Genotype(Genotype::unphased(0, 1)), SampleValue::Integer(30)])
    .build(&header)
    .unwrap();

    let mut output = Vec::new();
    {
        let mut writer = VcfWriter::new(&mut output, header);
        writer.write_header().unwrap();
        writer.write_record(&record).unwrap();
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
    let header = rich_header();
    let record = VcfRecordBuilder::new(
        "chr1",
        Pos::<One>::new(100).unwrap(),
        Alleles::snv(Base::A, Base::T).unwrap(),
    )
    .qual(29.5)
    .filter_pass()
    .info_integer("DP", 42)
    .format_keys(&["GT", "DP"])
    .add_sample(vec![SampleValue::Genotype(Genotype::unphased(0, 1)), SampleValue::Integer(25)])
    .build(&header)
    .unwrap();

    let mut bcf_output = Vec::new();
    {
        let mut writer = BcfWriter::new(&mut bcf_output, header, false);
        writer.write_header().unwrap();
        writer.write_vcf_record(&record).unwrap();
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
    let header = rich_header();
    let mut bcf_output = Vec::new();
    {
        let mut writer = BcfWriter::new(&mut bcf_output, header.clone(), false);
        writer.write_header().unwrap();

        // Use the encoder API directly
        let contig = ContigHandle(0); // chr1
        let dp_info = ScalarInfoHandle::<i32> {
            dict_idx: header.string_map().get("DP").unwrap() as u32,
            _marker: PhantomData,
        };
        let gt_fmt = GtFormatHandle { dict_idx: header.string_map().get("GT").unwrap() as u32 };
        let dp_fmt = ScalarFormatHandle::<i32> {
            dict_idx: header.string_map().get("DP").unwrap() as u32,
            _marker: PhantomData,
        };

        let alleles = Alleles::snv(Base::C, Base::G).unwrap();
        let mut enc = writer.record_encoder();
        alleles.begin_record(&mut enc, contig, Pos::<One>::new(500).unwrap(), Some(45.0));
        FilterHandle::PASS.encode(&mut enc);
        dp_info.encode(&mut enc, 100);
        enc.begin_samples(1);
        gt_fmt.encode(&mut enc, &Genotype::phased_diploid(0, 1));
        dp_fmt.encode(&mut enc, 80);
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

#[allow(dead_code)]
fn arb_genotype() -> impl Strategy<Value = Genotype> {
    prop_oneof![
        (0u16..4, 0u16..4).prop_map(|(a, b)| Genotype::unphased(a, b)),
        (0u16..4, 0u16..4).prop_map(|(a, b)| Genotype::phased_diploid(a, b)),
        (0u16..4).prop_map(Genotype::haploid),
        Just(Genotype::missing_diploid()),
    ]
}

// ── Proptest: BCF via write_record round-trips through noodles ─────────

proptest! {
    #![proptest_config(ProptestConfig::with_cases(50))]

    #[test]
    fn bcf_write_record_roundtrip(
        pos in 1u32..10_000_000,
        alleles in arb_alleles(),
        qual in proptest::option::of(0.1f32..1000.0),
        dp in 0i32..10000,
    ) {
        let header = rich_header();
        let pos = Pos::<One>::new(pos).unwrap();
        let record = VcfRecordBuilder::new("chr1", pos, alleles)
            .qual(qual.unwrap_or(0.0))
            .filter_pass()
            .info_integer("DP", dp)
            .build(&header)
            .unwrap();

        let mut bcf_output = Vec::new();
        {
            let mut writer = BcfWriter::new(&mut bcf_output, header, false);
            writer.write_header().unwrap();
            writer.write_vcf_record(&record).unwrap();
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
        let header = Arc::new(
            VcfHeader::builder()
                .add_contig("chr1", ContigDef { length: Some(250_000_000) }).unwrap()
                .add_info("DP", InfoDef {
                    number: Number::Count(1), typ: ValueType::Integer,
                    description: SmolStr::from("Depth"),
                }).unwrap()
                .build().unwrap(),
        );
        let pos = Pos::<One>::new(pos).unwrap();
        let record = VcfRecordBuilder::new("chr1", pos, alleles)
            .filter_pass()
            .info_integer("DP", dp)
            .build(&header)
            .unwrap();

        let mut output = Vec::new();
        {
            let mut writer = VcfWriter::new(&mut output, header);
            writer.write_header().unwrap();
            writer.write_record(&record).unwrap();
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
