//! Proptest round-trips for the unified Writer in BCF mode.
//!
//! Generates random records, writes them with seqair's unified `Writer`, and
//! verifies that noodles can parse the resulting BCF to identical site-level
//! fields.  This replaced the old two-path equivalence test (`BcfWriter` vs
//! `BcfRecordEncoder`) — both old paths have been removed; the unified Writer is
//! the only path now.
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

use hegel::prelude::*;
use seqair::vcf::record_encoder::{FormatFieldDef, FormatGt, FormatInt, InfoFieldDef, InfoInt};
use seqair::vcf::{
    Alleles, ContigDef, Genotype, Number, OutputFormat, ValueType, VcfHeader, Writer,
};
use seqair_types::{Base, Pos1};
use std::sync::Arc;

struct TestSetup {
    header: Arc<VcfHeader>,
    contig: seqair::vcf::ContigId,
    dp_info: InfoInt,
    gt_fmt: FormatGt,
    dp_fmt: FormatInt,
}

fn make_setup() -> TestSetup {
    let mut builder = VcfHeader::builder();
    let contig = builder.register_contig("chr1", ContigDef { length: Some(250_000_000) }).unwrap();
    let mut builder = builder.infos();
    let dp_info = builder
        .register_info(&InfoFieldDef::new(
            "DP",
            Number::Count(1),
            ValueType::Integer,
            "Total Depth",
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
    let mut builder = builder.samples();
    builder.add_sample("sample1").unwrap();
    let header = Arc::new(builder.build().unwrap());
    TestSetup { header, contig, dp_info, gt_fmt, dp_fmt }
}

fn write_bcf(
    setup: &TestSetup,
    pos: u32,
    alleles: &Alleles,
    qual: f32,
    depth: i32,
    gt: &Genotype,
) -> Vec<u8> {
    let mut output = Vec::new();
    let writer = Writer::new(&mut output, OutputFormat::Bcf);
    let mut writer = writer.write_header(&setup.header).unwrap();

    let mut enc = writer
        .begin_record(&setup.contig, Pos1::new(pos).unwrap(), alleles, Some(qual))
        .unwrap()
        .filter_pass();
    setup.dp_info.encode(&mut enc, depth);
    let mut enc = enc.begin_samples();
    setup.gt_fmt.encode(&mut enc, std::slice::from_ref(gt)).unwrap();
    setup.dp_fmt.encode(&mut enc, &[depth]).unwrap();
    enc.emit().unwrap();

    writer.finish().unwrap();
    output
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
        let alt_alleles: Vec<String> = rec.alternate_bases().as_ref().to_vec();
        records.push(ParsedRecord { pos, ref_allele, alt_alleles, qual_bits: qual });
    }
    records
}

/// The four concrete bases a REF or ALT allele can hold.
const ACGT: [Base; 4] = [Base::A, Base::C, Base::G, Base::T];

fn arb_base() -> impl PrintableGenerator<Base> {
    gs::sampled_from(&ACGT).print_as_debug()
}

// ── Deterministic test ─────────────────────────────────────────────────

// r[verify bcf_writer.coordinate_system]
#[test]
fn bcf_snv_parseable_by_noodles() {
    let setup = make_setup();
    let alleles = Alleles::snv(Base::A, Base::T).unwrap();
    let gt = Genotype::unphased(0, 1);

    let bcf_bytes = write_bcf(&setup, 100, &alleles, 30.0, 50, &gt);

    let records = parse_bcf_with_noodles(&bcf_bytes);
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].pos, 99, "POS should be 0-based 99");
    assert_eq!(records[0].ref_allele, "A");
    assert_eq!(records[0].alt_alleles, vec!["T"]);
    assert_eq!(records[0].qual_bits, Some(30.0f32.to_bits()));
}

// ── Proptest ───────────────────────────────────────────────────────────

/// Unified Writer BCF output must be parseable by noodles with correct site fields.
#[hegel::test(test_cases = 30)]
fn bcf_writer_noodles_roundtrip(tc: TestCase) {
    let pos = tc.draw(gs::integers::<u32>().min_value(1).max_value(9999999));
    let ref_base = tc.draw(arb_base());
    let alt_base = tc.draw(arb_base());
    let qual = tc.draw(gs::floats::<f32>().min_value(1.0).max_value_exclusive(1000.0));
    let depth = tc.draw(gs::integers::<i32>().min_value(1).max_value(999));
    let gt_a0 = tc.draw(gs::integers::<u16>().max_value(1));
    let gt_a1 = tc.draw(gs::integers::<u16>().max_value(1));
    let phased = tc.draw(gs::booleans());
    tc.assume(ref_base != alt_base);

    let setup = make_setup();
    let alleles = Alleles::snv(ref_base, alt_base).unwrap();
    let gt = if phased {
        Genotype::phased_diploid(gt_a0, gt_a1)
    } else {
        Genotype::unphased(gt_a0, gt_a1)
    };

    let bcf_bytes = write_bcf(&setup, pos, &alleles, qual, depth, &gt);
    let records = parse_bcf_with_noodles(&bcf_bytes);

    assert_eq!(records.len(), 1, "should produce exactly 1 record");

    let expected_pos = pos as i32 - 1;
    assert_eq!(records[0].pos, expected_pos, "POS (0-based)");
    assert_eq!(&records[0].ref_allele, ref_base.as_str(), "REF");
    assert_eq!(&records[0].alt_alleles, &vec![alt_base.as_str().to_string()], "ALT");
    assert_eq!(records[0].qual_bits, Some(qual.to_bits()), "QUAL");
}

// r[verify record_encoder.info_dedup]
// r[verify record_encoder.format_dedup]
/// Writing INFO/FORMAT fields twice must produce output identical to writing
/// only the final value. Tested for both VCF text and BCF binary.
#[hegel::test(test_cases = 30)]
fn dedup_matches_single_write(tc: TestCase) {
    let pos = tc.draw(gs::integers::<u32>().min_value(1).max_value(9999999));
    let ref_base = tc.draw(arb_base());
    let alt_base = tc.draw(arb_base());
    let first_depth = tc.draw(gs::integers::<i32>().min_value(1).max_value(999));
    let final_depth = tc.draw(gs::integers::<i32>().min_value(1).max_value(999));
    let gt_a0 = tc.draw(gs::integers::<u16>().max_value(1));
    let gt_a1 = tc.draw(gs::integers::<u16>().max_value(1));
    let format_idx = tc.draw(gs::integers::<usize>().max_value(2));
    tc.assume(ref_base != alt_base);
    let output_format = match format_idx {
        0 => OutputFormat::Bcf,
        1 => OutputFormat::Vcf,
        _ => OutputFormat::VcfGz,
    };
    let setup = make_setup();
    let alleles = Alleles::snv(ref_base, alt_base).unwrap();
    let gt = Genotype::unphased(gt_a0, gt_a1);

    // Write with dedup: write first_depth then overwrite with final_depth
    let dedup_output = {
        let mut output = Vec::new();
        let writer = Writer::new(&mut output, output_format);
        let mut writer = writer.write_header(&setup.header).unwrap();
        let mut enc = writer
            .begin_record(&setup.contig, Pos1::new(pos).unwrap(), &alleles, Some(30.0))
            .unwrap()
            .filter_pass();
        setup.dp_info.encode(&mut enc, first_depth);
        setup.dp_info.encode(&mut enc, final_depth); // overwrite
        let mut enc = enc.begin_samples();
        setup.gt_fmt.encode(&mut enc, std::slice::from_ref(&gt)).unwrap();
        setup.dp_fmt.encode(&mut enc, &[first_depth]).unwrap();
        setup.dp_fmt.encode(&mut enc, &[final_depth]).unwrap(); // overwrite
        enc.emit().unwrap();
        writer.finish().unwrap();
        output
    };

    // Write without dedup: write final_depth directly
    let clean_output = {
        let mut output = Vec::new();
        let writer = Writer::new(&mut output, output_format);
        let mut writer = writer.write_header(&setup.header).unwrap();
        let mut enc = writer
            .begin_record(&setup.contig, Pos1::new(pos).unwrap(), &alleles, Some(30.0))
            .unwrap()
            .filter_pass();
        setup.dp_info.encode(&mut enc, final_depth);
        let mut enc = enc.begin_samples();
        setup.gt_fmt.encode(&mut enc, std::slice::from_ref(&gt)).unwrap();
        setup.dp_fmt.encode(&mut enc, &[final_depth]).unwrap();
        enc.emit().unwrap();
        writer.finish().unwrap();
        output
    };

    assert_eq!(
        dedup_output, clean_output,
        "dedup output must be byte-identical to single-write output"
    );
}
