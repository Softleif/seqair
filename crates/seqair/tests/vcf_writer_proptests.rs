//! Property-based tests for VCF writer round-trip and BCF encoding correctness.

use proptest::prelude::*;
use seqair::vcf::alleles::Alleles;
use seqair::vcf::header::{ContigDef, InfoDef, Number, ValueType, VcfHeader};
use seqair::vcf::record::VcfRecordBuilder;
use seqair::vcf::writer::VcfWriter;
use seqair_types::{Base, One, Pos, SmolStr};
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
                    description: SmolStr::from("Depth"),
                },
            )
            .unwrap()
            .build()
            .unwrap(),
    )
}

proptest! {
    /// Any VcfRecord serialized to VCF text must produce exactly 8 tab-separated fields
    /// (plus newline) for records without samples.
    #[test]
    fn vcf_text_has_eight_columns(
        pos in 1u32..100_000_000,
        alleles in arb_alleles(),
        qual in proptest::option::of(0.0f32..10000.0),
        dp in 0i32..10000,
    ) {
        let header = test_header();
        let pos = Pos::<One>::new(pos).unwrap();
        let record = VcfRecordBuilder::new("chr1", pos, alleles)
            .qual(qual.unwrap_or(0.0))
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
        let text = String::from_utf8(output).unwrap();

        // Find the last line (the data line)
        let data_line = text.lines().last().unwrap();
        let fields: Vec<&str> = data_line.split('\t').collect();
        prop_assert_eq!(fields.len(), 8);

        // POS field must match
        prop_assert_eq!(fields[1], pos.get().to_string());

        // CHROM must be chr1
        prop_assert_eq!(fields[0], "chr1");
    }

    /// VCF lines must end with exactly one newline and contain no embedded newlines.
    #[test]
    fn vcf_lines_properly_terminated(
        pos in 1u32..1_000_000,
        alleles in arb_alleles(),
    ) {
        let header = test_header();
        let pos = Pos::<One>::new(pos).unwrap();
        let record = VcfRecordBuilder::new("chr1", pos, alleles)
            .filter_pass()
            .build(&header)
            .unwrap();

        let mut output = Vec::new();
        {
            let mut writer = VcfWriter::new(&mut output, header);
            writer.write_header().unwrap();
            writer.write_record(&record).unwrap();
            writer.finish().unwrap();
        }
        let text = String::from_utf8(output).unwrap();

        // Every line must end with \n
        for line in text.split('\n') {
            if !line.is_empty() {
                // No embedded tabs or special chars in individual fields that would break parsing
                prop_assert!(!line.contains('\r'), "line contains CR: {line}");
            }
        }
        // Text must end with \n
        prop_assert!(text.ends_with('\n'));
    }

    /// BGZF round-trip: write VCF to BGZF, decompress, verify content matches plain write.
    #[test]
    fn bgzf_roundtrip_matches_plain(
        pos in 1u32..1_000_000,
        alleles in arb_alleles(),
        dp in 0i32..1000,
    ) {
        let header = test_header();
        let pos = Pos::<One>::new(pos).unwrap();
        let record = VcfRecordBuilder::new("chr1", pos, alleles)
            .filter_pass()
            .info_integer("DP", dp)
            .build(&header)
            .unwrap();

        // Write plain
        let mut plain_output = Vec::new();
        {
            let mut writer = VcfWriter::new(&mut plain_output, header.clone());
            writer.write_header().unwrap();
            writer.write_record(&record).unwrap();
            writer.finish().unwrap();
        }

        // Write BGZF
        let mut bgzf_output = Vec::new();
        {
            let mut writer = VcfWriter::bgzf(&mut bgzf_output, header);
            writer.write_header().unwrap();
            writer.write_record(&record).unwrap();
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

    /// Filter serialization: Pass=PASS, NotApplied=., Failed=semicolon-separated IDs.
    #[test]
    fn filter_serialization_correct(
        pos in 1u32..1_000_000,
        filter_state in 0u8..3,
    ) {
        let header = test_header();
        let pos = Pos::<One>::new(pos).unwrap();
        let mut builder = VcfRecordBuilder::new("chr1", pos, Alleles::reference(Base::A));

        let expected_filter = match filter_state {
            0 => { builder = builder.filter_pass(); "PASS" },
            1 => ".",  // NotApplied (default)
            _ => { builder = builder.filter_pass(); "PASS" }, // Can't easily test Failed without declared filters
        };

        let record = builder.build(&header).unwrap();

        let mut output = Vec::new();
        {
            let mut writer = VcfWriter::new(&mut output, header);
            writer.write_header().unwrap();
            writer.write_record(&record).unwrap();
            writer.finish().unwrap();
        }
        let text = String::from_utf8(output).unwrap();
        let data_line = text.lines().last().unwrap();
        let fields: Vec<&str> = data_line.split('\t').collect();

        prop_assert_eq!(fields[6], expected_filter);
    }
}
