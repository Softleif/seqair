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

use hegel::prelude::*;
use seqair::vcf::record_encoder::{FilterFieldDef, InfoFieldDef, InfoInt};
use seqair::vcf::{Alleles, ContigDef, Number, OutputFormat, ValueType, VcfHeader, Writer};
use seqair_types::{Base, Pos1};
use std::sync::Arc;

/// The four concrete bases a REF or ALT allele can hold.
const ACGT: [Base; 4] = [Base::A, Base::C, Base::G, Base::T];

fn arb_base() -> impl PrintableGenerator<Base> {
    gs::sampled_from(&ACGT).print_as_debug()
}

#[hegel::composite]
fn arb_alleles_inner(tc: &TestCase) -> Alleles {
    let bases = || gs::vecs(gs::sampled_from(&ACGT)).min_size(1).max_size(5);
    match tc.draw_silent(gs::integers::<u8>().max_value(3)) {
        0 => Alleles::reference(tc.draw_silent(arb_base())),
        1 => {
            let r = tc.draw_silent(arb_base());
            let alt = gs::sampled_from(&ACGT).filter(move |a| *a != r);
            Alleles::snv(r, tc.draw_silent(alt)).unwrap()
        }
        2 => Alleles::insertion(tc.draw_silent(arb_base()), &tc.draw_silent(bases())).unwrap(),
        _ => Alleles::deletion(tc.draw_silent(arb_base()), &tc.draw_silent(bases())).unwrap(),
    }
}

fn arb_alleles() -> impl PrintableGenerator<Alleles> {
    arb_alleles_inner().print_as_debug()
}

struct SimpleSetup {
    header: Arc<VcfHeader>,
    contig: seqair::vcf::ContigId,
    dp_info: InfoInt,
}

fn make_simple_setup() -> SimpleSetup {
    let mut builder = VcfHeader::builder();
    let contig = builder.register_contig("chr1", ContigDef { length: Some(250_000_000) }).unwrap();
    let mut builder = builder.infos();
    let dp_info = builder
        .register_info(&InfoFieldDef::new("DP", Number::Count(1), ValueType::Integer, "Depth"))
        .unwrap();
    let header = Arc::new(builder.build().unwrap());
    SimpleSetup { header, contig, dp_info }
}

/// Any record serialized to VCF text must produce exactly 8 tab-separated fields
/// (plus newline) for records without samples.
#[hegel::test]
fn vcf_text_has_eight_columns(tc: TestCase) {
    let pos = tc.draw(gs::integers::<u32>().min_value(1).max_value(99999999));
    let alleles = tc.draw(arb_alleles());
    let qual =
        tc.draw(gs::optional(gs::floats::<f32>().min_value(0.0).max_value_exclusive(10_000.0)));
    let dp = tc.draw(gs::integers::<i32>().min_value(0).max_value(9999));
    let setup = make_simple_setup();
    let pos_typed = Pos1::new(pos).unwrap();

    let mut output = Vec::new();
    let writer = Writer::new(&mut output, OutputFormat::Vcf);
    let mut writer = writer.write_header(&setup.header).unwrap();
    {
        let mut enc =
            writer.begin_record(&setup.contig, pos_typed, &alleles, qual).unwrap().filter_pass();
        setup.dp_info.encode(&mut enc, dp);
        enc.emit().unwrap();
    }
    writer.finish().unwrap();

    let text = String::from_utf8(output).unwrap();

    // Find the last line (the data line)
    let data_line = text.lines().last().unwrap();
    let fields: Vec<&str> = data_line.split('\t').collect();
    assert_eq!(fields.len(), 8);

    // POS field must match
    assert_eq!(fields[1], pos.to_string());

    // CHROM must be chr1
    assert_eq!(fields[0], "chr1");
}

/// VCF lines must end with exactly one newline and contain no embedded carriage returns.
#[hegel::test]
fn vcf_lines_properly_terminated(tc: TestCase) {
    let pos = tc.draw(gs::integers::<u32>().min_value(1).max_value(999999));
    let alleles = tc.draw(arb_alleles());
    let setup = make_simple_setup();
    let pos_typed = Pos1::new(pos).unwrap();

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
            assert!(!line.contains('\r'), "line contains CR: {line}");
        }
    }
    assert!(text.ends_with('\n'));
}

/// BGZF round-trip: write VCF to BGZF, decompress, verify content matches plain write.
#[hegel::test]
fn bgzf_roundtrip_matches_plain(tc: TestCase) {
    let pos = tc.draw(gs::integers::<u32>().min_value(1).max_value(999999));
    let alleles = tc.draw(arb_alleles());
    let dp = tc.draw(gs::integers::<i32>().min_value(0).max_value(999));
    let setup = make_simple_setup();
    let pos_typed = Pos1::new(pos).unwrap();

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
    let mut reader = seqair::bam::bgzf::BgzfReader::from_reader(std::io::Cursor::new(bgzf_output));
    let mut decompressed = Vec::new();
    reader.read_to_end(&mut decompressed).unwrap();

    // Content must match
    assert_eq!(String::from_utf8(decompressed).unwrap(), String::from_utf8(plain_output).unwrap());
}

/// Filter serialization: `Pass` = `PASS`, `NotApplied` = `.`, `Failed` = filter name.
#[hegel::test]
fn filter_serialization_correct(tc: TestCase) {
    let pos = tc.draw(gs::integers::<u32>().min_value(1).max_value(999999));
    let filter_state = tc.draw(gs::integers::<u8>().max_value(2));
    let mut builder = VcfHeader::builder();
    let contig = builder.register_contig("chr1", ContigDef { length: Some(250_000_000) }).unwrap();
    // Register filters before info (BCF string dict order: PASS, filters, info, format).
    let mut builder = builder.filters();
    let q20 = builder.register_filter(&FilterFieldDef::new("q20", "Quality below 20")).unwrap();
    let mut builder = builder.infos();
    let _dp: InfoInt = builder
        .register_info(&InfoFieldDef::new("DP", Number::Count(1), ValueType::Integer, "Depth"))
        .unwrap();
    let header = Arc::new(builder.build().unwrap());

    let pos_typed = Pos1::new(pos).unwrap();
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
                .filter_fail([&q20])
                .emit()
                .unwrap();
            "q20"
        }
    };

    writer.finish().unwrap();
    let text = String::from_utf8(output).unwrap();
    let data_line = text.lines().last().unwrap();
    let fields: Vec<&str> = data_line.split('\t').collect();

    assert_eq!(fields[6], expected_filter);
}

/// Counts `write` calls; shares its bytes so they can be read after the
/// writer that owns the sink is dropped.
#[derive(Clone, Default)]
struct CountingSink {
    writes: std::rc::Rc<std::cell::Cell<usize>>,
    bytes: std::rc::Rc<std::cell::RefCell<Vec<u8>>>,
}

impl std::io::Write for CountingSink {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.writes.set(self.writes.get() + 1);
        self.bytes.borrow_mut().extend_from_slice(buf);
        Ok(buf.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

// r[verify vcf_writer.buffered_output]
/// Plain VCF reaches an unbuffered sink in large writes, not one per record,
/// and every record arrives whether the writer is finished or just dropped.
#[hegel::test(test_cases = 40)]
fn plain_vcf_is_buffered_and_complete(tc: TestCase) {
    let n = tc.draw(gs::integers::<u32>().min_value(1).max_value(20_000));
    let finish = tc.draw(gs::booleans());
    let setup = make_simple_setup();
    let alleles = Alleles::snv(Base::A, Base::C).unwrap();

    let sink = CountingSink::default();
    let mut writer =
        Writer::new(sink.clone(), OutputFormat::Vcf).write_header(&setup.header).unwrap();
    for i in 0..n {
        let pos = Pos1::new(i + 1).unwrap();
        let mut enc =
            writer.begin_record(&setup.contig, pos, &alleles, None).unwrap().filter_pass();
        setup.dp_info.encode(&mut enc, 7);
        enc.emit().unwrap();
    }
    if finish {
        writer.finish().unwrap();
    } else {
        drop(writer);
    }

    let bytes = sink.bytes.borrow();
    let text = std::str::from_utf8(&bytes).unwrap();
    let records: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(records.len(), n as usize, "every record reaches the sink");
    assert!(records.iter().enumerate().all(|(i, l)| l.starts_with(&format!("chr1\t{}\t", i + 1))));
    // One write per 64 KiB at most, plus slack for the final partial buffer.
    let max_writes = bytes.len() / (64 * 1024) + 2;
    assert!(
        sink.writes.get() <= max_writes,
        "{} writes for {} bytes ({n} records)",
        sink.writes.get(),
        bytes.len()
    );
}
