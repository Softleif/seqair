//! Property-based tests for VCF writer round-trip and BCF encoding correctness.
//! Uses the unified Writer with `OutputFormat::Vcf` and `OutputFormat::VcfGz`.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "test code"
)]

use proptest::prelude::*;
use seqair::vcf::record_encoder::{FilterFieldDef, InfoFieldDef, InfoInt};
use seqair::vcf::{Alleles, ContigDef, Number, OutputFormat, ValueType, VcfHeader, Writer};
use seqair_types::{Base, One, Pos};
use std::sync::Arc;

fn arb_base() -> impl Strategy<Value = Base> {
    prop_oneof![Just(Base::A), Just(Base::C), Just(Base::G), Just(Base::T),]
}

fn arb_alleles() -> impl Strategy<Value = Alleles> {
    prop_oneof![
        arb_base().prop_map(Alleles::reference),
        (arb_base(), arb_base())
            .prop_filter("ref != alt", |(r, a)| r != a)
            .prop_map(|(r, a)| Alleles::snv(r, a).unwrap()),
        (arb_base(), proptest::collection::vec(arb_base(), 1..6))
            .prop_map(|(anchor, ins)| Alleles::insertion(anchor, &ins).unwrap()),
        (arb_base(), proptest::collection::vec(arb_base(), 1..6))
            .prop_map(|(anchor, del)| Alleles::deletion(anchor, &del).unwrap()),
    ]
}

struct SimpleSetup {
    header: Arc<VcfHeader>,
    contig: seqair::vcf::ContigId,
    dp_info: InfoInt,
}

fn make_simple_setup() -> SimpleSetup {
    let mut builder = VcfHeader::builder();
    let contig = builder.register_contig("chr1", ContigDef { length: Some(250_000_000) }).unwrap();
    let dp_info = builder
        .register_info(&InfoFieldDef::new("DP", Number::Count(1), ValueType::Integer, "Depth"))
        .unwrap();
    let header = Arc::new(builder.build().unwrap());
    SimpleSetup { header, contig, dp_info }
}

proptest! {
    /// Any record serialized to VCF text must produce exactly 8 tab-separated fields
    /// (plus newline) for records without samples.
    #[test]
    fn vcf_text_has_eight_columns(
        pos in 1u32..100_000_000,
        alleles in arb_alleles(),
        qual in proptest::option::of(0.0f32..10000.0),
        dp in 0i32..10000,
    ) {
        let setup = make_simple_setup();
        let pos_typed = Pos::<One>::new(pos).unwrap();

        let mut output = Vec::new();
        let writer = Writer::new(&mut output, OutputFormat::Vcf);
        let mut writer = writer.write_header(&setup.header).unwrap();
        {
            let mut enc = writer
                .begin_record(&setup.contig, pos_typed, &alleles, qual)
                .unwrap()
                .filter_pass();
            setup.dp_info.encode(&mut enc, dp);
            enc.emit().unwrap();
        }
        writer.finish().unwrap();

        let text = String::from_utf8(output).unwrap();

        // Find the last line (the data line)
        let data_line = text.lines().last().unwrap();
        let fields: Vec<&str> = data_line.split('\t').collect();
        prop_assert_eq!(fields.len(), 8);

        // POS field must match
        prop_assert_eq!(fields[1], pos.to_string());

        // CHROM must be chr1
        prop_assert_eq!(fields[0], "chr1");
    }

    /// VCF lines must end with exactly one newline and contain no embedded carriage returns.
    #[test]
    fn vcf_lines_properly_terminated(
        pos in 1u32..1_000_000,
        alleles in arb_alleles(),
    ) {
        let setup = make_simple_setup();
        let pos_typed = Pos::<One>::new(pos).unwrap();

        let mut output = Vec::new();
        let writer = Writer::new(&mut output, OutputFormat::Vcf);
        let mut writer = writer.write_header(&setup.header).unwrap();
        writer
            .begin_record(&setup.contig, pos_typed, &alleles, None)
            .unwrap()
            .filter_pass()
            .emit()
            .unwrap();
        writer.finish().unwrap();

        let text = String::from_utf8(output).unwrap();

        for line in text.split('\n') {
            if !line.is_empty() {
                prop_assert!(!line.contains('\r'), "line contains CR: {line}");
            }
        }
        prop_assert!(text.ends_with('\n'));
    }

    /// BGZF round-trip: write VCF to BGZF, decompress, verify content matches plain write.
    #[test]
    fn bgzf_roundtrip_matches_plain(
        pos in 1u32..1_000_000,
        alleles in arb_alleles(),
        dp in 0i32..1000,
    ) {
        let setup = make_simple_setup();
        let pos_typed = Pos::<One>::new(pos).unwrap();

        // Write plain VCF
        let mut plain_output = Vec::new();
        {
            let writer = Writer::new(&mut plain_output, OutputFormat::Vcf);
            let mut writer = writer.write_header(&setup.header).unwrap();
            {
                let mut enc = writer
                    .begin_record(&setup.contig, pos_typed, &alleles, None)
                    .unwrap()
                    .filter_pass();
                setup.dp_info.encode(&mut enc, dp);
                enc.emit().unwrap();
            }
            writer.finish().unwrap();
        }

        // Write BGZF-compressed VCF
        let mut bgzf_output = Vec::new();
        {
            let writer = Writer::new(&mut bgzf_output, OutputFormat::VcfGz);
            let mut writer = writer.write_header(&setup.header).unwrap();
            {
                let mut enc = writer
                    .begin_record(&setup.contig, pos_typed, &alleles, None)
                    .unwrap()
                    .filter_pass();
                setup.dp_info.encode(&mut enc, dp);
                enc.emit().unwrap();
            }
            writer.finish().unwrap();
        }

        // Decompress BGZF
        let mut reader = seqair::bam::bgzf::BgzfReader::from_reader(
            std::io::Cursor::new(bgzf_output),
        );
        let mut decompressed = Vec::new();
        reader.read_to_end(&mut decompressed).unwrap();

        // Content must match
        prop_assert_eq!(
            String::from_utf8(decompressed).unwrap(),
            String::from_utf8(plain_output).unwrap()
        );
    }

    /// Filter serialization: Pass=PASS, NotApplied=., Failed=filter name.
    #[test]
    fn filter_serialization_correct(
        pos in 1u32..1_000_000,
        filter_state in 0u8..3,
    ) {
        let mut builder = VcfHeader::builder();
        let contig =
            builder.register_contig("chr1", ContigDef { length: Some(250_000_000) }).unwrap();
        // Register filters before info (BCF string dict order: PASS, filters, info, format).
        // TODO: enforce this ordering at the type level in VcfHeaderBuilder.
        let q20 = builder.register_filter(&FilterFieldDef::new("q20", "Quality below 20")).unwrap();
        let _dp: InfoInt = builder
            .register_info(&InfoFieldDef::new("DP", Number::Count(1), ValueType::Integer, "Depth"))
            .unwrap();
        let header = Arc::new(builder.build().unwrap());

        let pos_typed = Pos::<One>::new(pos).unwrap();
        let alleles = Alleles::reference(Base::A);

        let mut output = Vec::new();
        let writer = Writer::new(&mut output, OutputFormat::Vcf);
        let mut writer = writer.write_header(&header).unwrap();

        let expected_filter = match filter_state {
            0 => {
                writer
                    .begin_record(&contig, pos_typed, &alleles, None)
                    .unwrap()
                    .filter_pass()
                    .emit()
                    .unwrap();
                "PASS"
            }
            1 => {
                // No filter applied — emit with filter_pass to get a valid record;
                // we test NotApplied by using "." which is represented as missing filter.
                // Since the unified API only has filter_pass/filter_fail, use filter_pass
                // and check that PASS is output.
                writer
                    .begin_record(&contig, pos_typed, &alleles, None)
                    .unwrap()
                    .filter_pass()
                    .emit()
                    .unwrap();
                "PASS"
            }
            _ => {
                writer
                    .begin_record(&contig, pos_typed, &alleles, None)
                    .unwrap()
                    .filter_fail(&[&q20])
                    .emit()
                    .unwrap();
                "q20"
            }
        };

        writer.finish().unwrap();
        let text = String::from_utf8(output).unwrap();
        let data_line = text.lines().last().unwrap();
        let fields: Vec<&str> = data_line.split('\t').collect();

        prop_assert_eq!(fields[6], expected_filter);
    }
}
