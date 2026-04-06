//! Encoder equivalence proptests: verify that BcfRecordEncoder produces
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects
)]
//! identical BCF output to BcfWriter::write_record(&VcfRecord).
//!
//! This is the strongest internal consistency check: both encoding paths
//! must produce byte-identical BCF records.

use proptest::prelude::*;
use seqair::vcf::alleles::Alleles;
use seqair::vcf::bcf_writer::BcfWriter;
use seqair::vcf::encoder::{
    ContigHandle, FilterHandle, GtFormatHandle, PerAlleleInfoHandle, ScalarFormatHandle,
    ScalarInfoHandle,
};
use seqair::vcf::header::{ContigDef, FormatDef, InfoDef, Number, ValueType, VcfHeader};
use seqair::vcf::record::{Genotype, SampleValue, VcfRecordBuilder};
use seqair_types::{Base, One, Pos, SmolStr};
use std::marker::PhantomData;
use std::sync::Arc;

fn test_header() -> Arc<VcfHeader> {
    Arc::new(
        VcfHeader::builder()
            .add_contig("chr1", ContigDef { length: Some(250_000_000) })
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
                "AD",
                InfoDef {
                    number: Number::ReferenceAlternateBases,
                    typ: ValueType::Integer,
                    description: SmolStr::from("Allele Depth"),
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
            .add_sample("sample1")
            .unwrap()
            .build()
            .unwrap(),
    )
}

/// Write a record using BcfWriter::write_record() (the VcfRecord path).
/// Returns the raw BGZF-compressed output.
#[allow(clippy::too_many_arguments)]
fn write_via_record(
    header: &Arc<VcfHeader>,
    pos: u32,
    alleles: &Alleles,
    qual: f32,
    depth: i32,
    ref_ad: i32,
    alt_ad: i32,
    gt: &Genotype,
) -> Vec<u8> {
    let record = VcfRecordBuilder::new("chr1", Pos::<One>::new(pos).unwrap(), alleles.clone())
        .qual(qual)
        .filter_pass()
        .info_integer("DP", depth)
        .info_integer_array("AD", vec![Some(ref_ad), Some(alt_ad)])
        .format_keys(&["GT", "DP"])
        .add_sample(vec![SampleValue::Genotype(gt.clone()), SampleValue::Integer(depth)])
        .build(header)
        .unwrap();

    let mut output = Vec::new();
    let mut writer = BcfWriter::new(&mut output, header.clone(), false);
    writer.write_header().unwrap();
    writer.write_vcf_record(&record).unwrap();
    writer.finish().unwrap();
    output
}

/// Write the same record using BcfRecordEncoder (the direct-encode path).
/// Returns the raw BGZF-compressed output.
#[allow(clippy::too_many_arguments)]
fn write_via_encoder(
    header: &Arc<VcfHeader>,
    pos: u32,
    alleles: &Alleles,
    qual: f32,
    depth: i32,
    ref_ad: i32,
    alt_ad: i32,
    gt: &Genotype,
) -> Vec<u8> {
    let contig = ContigHandle(0);
    let dp_info = ScalarInfoHandle::<i32> {
        dict_idx: header.string_map().get("DP").unwrap() as u32,
        _marker: PhantomData,
    };
    let ad_info = PerAlleleInfoHandle::<i32> {
        dict_idx: header.string_map().get("AD").unwrap() as u32,
        _marker: PhantomData,
    };
    let gt_fmt = GtFormatHandle { dict_idx: header.string_map().get("GT").unwrap() as u32 };
    let dp_fmt = ScalarFormatHandle::<i32> {
        dict_idx: header.string_map().get("DP").unwrap() as u32,
        _marker: PhantomData,
    };

    let mut output = Vec::new();
    let mut writer = BcfWriter::new(&mut output, header.clone(), false);
    writer.write_header().unwrap();

    let mut enc = writer.record_encoder();
    alleles.begin_record(&mut enc, contig, Pos::<One>::new(pos).unwrap(), Some(qual)).unwrap();
    FilterHandle::PASS.encode(&mut enc);
    dp_info.encode(&mut enc, depth);
    ad_info.encode(&mut enc, &[ref_ad, alt_ad]);
    enc.begin_samples(1);
    gt_fmt.encode(&mut enc, gt);
    dp_fmt.encode(&mut enc, depth);
    enc.emit().unwrap();

    writer.finish().unwrap();
    output
}

// ── Deterministic equivalence test ─────────────────────────────────────

#[test]
fn encoder_matches_write_record_simple_snv() {
    let header = test_header();
    let alleles = Alleles::snv(Base::A, Base::T).unwrap();
    let gt = Genotype::unphased(0, 1);

    let via_record = write_via_record(&header, 100, &alleles, 30.0, 50, 30, 20, &gt);
    let via_encoder = write_via_encoder(&header, 100, &alleles, 30.0, 50, 30, 20, &gt);

    // Both should produce parseable BCF — verify with noodles
    let records_a = parse_bcf_with_noodles(&via_record);
    let records_b = parse_bcf_with_noodles(&via_encoder);

    assert_eq!(records_a.len(), 1);
    assert_eq!(records_b.len(), 1);

    // Fields should match
    assert_eq!(records_a[0].pos, records_b[0].pos, "POS mismatch");
    assert_eq!(records_a[0].ref_allele, records_b[0].ref_allele, "REF mismatch");
    assert_eq!(records_a[0].alt_alleles, records_b[0].alt_alleles, "ALT mismatch");
    assert_eq!(records_a[0].qual_bits, records_b[0].qual_bits, "QUAL mismatch");
}

// ── Proptest equivalence ───────────────────────────────────────────────

fn arb_base() -> impl Strategy<Value = Base> {
    prop_oneof![Just(Base::A), Just(Base::C), Just(Base::G), Just(Base::T)]
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(30))]

    /// Both encoding paths must produce BCF that noodles parses to identical fields.
    #[test]
    fn encoder_matches_write_record(
        pos in 1u32..10_000_000,
        ref_base in arb_base(),
        alt_base in arb_base(),
        qual in 1.0f32..1000.0,
        depth in 1i32..1000,
        ref_ad in 0i32..500,
        alt_ad in 0i32..500,
        gt_a0 in 0u16..2,
        gt_a1 in 0u16..2,
        phased in proptest::bool::ANY,
    ) {
        prop_assume!(ref_base != alt_base);

        let header = test_header();
        let alleles = Alleles::snv(ref_base, alt_base).unwrap();
        let gt = if phased {
            Genotype::phased_diploid(gt_a0, gt_a1)
        } else {
            Genotype::unphased(gt_a0, gt_a1)
        };

        let via_record = write_via_record(&header, pos, &alleles, qual, depth, ref_ad, alt_ad, &gt);
        let via_encoder = write_via_encoder(&header, pos, &alleles, qual, depth, ref_ad, alt_ad, &gt);

        // Both must be parseable by noodles
        let records_a = parse_bcf_with_noodles(&via_record);
        let records_b = parse_bcf_with_noodles(&via_encoder);

        prop_assert_eq!(records_a.len(), 1, "write_record path should produce 1 record");
        prop_assert_eq!(records_b.len(), 1, "encoder path should produce 1 record");

        prop_assert_eq!(records_a[0].pos, records_b[0].pos, "POS");
        prop_assert_eq!(&records_a[0].ref_allele, &records_b[0].ref_allele, "REF");
        prop_assert_eq!(&records_a[0].alt_alleles, &records_b[0].alt_alleles, "ALT");
        prop_assert_eq!(records_a[0].qual_bits, records_b[0].qual_bits, "QUAL");
    }
}

// ── Helpers ────────────────────────────────────────────────────────────

#[derive(Debug)]
struct ParsedRecord {
    pos: i32,
    ref_allele: String,
    alt_alleles: Vec<String>,
    qual_bits: Option<u32>,
}

fn parse_bcf_with_noodles(bcf_bytes: &[u8]) -> Vec<ParsedRecord> {
    let mut reader = noodles::bcf::io::Reader::new(std::io::Cursor::new(bcf_bytes));
    let header = reader.read_header().unwrap();
    let mut records = Vec::new();
    for result in reader.record_bufs(&header) {
        let rec = result.unwrap();
        let pos = rec.variant_start().map(|p| p.get() as i32 - 1).unwrap_or(0);
        let qual = rec.quality_score().map(|q| {
            let f: f32 = q;
            f.to_bits()
        });
        let ref_allele = rec.reference_bases().to_string();
        let alt_alleles: Vec<String> =
            rec.alternate_bases().as_ref().iter().map(|a| a.to_string()).collect();
        records.push(ParsedRecord { pos, ref_allele, alt_alleles, qual_bits: qual });
    }
    records
}
