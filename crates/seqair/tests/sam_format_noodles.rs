//! `format_sam_record` against noodles' SAM writer, field by field.
//!
//! noodles is an independent implementation, but it is not htslib, so the
//! comparison is over what both print canonically. The records come from the
//! generator's `NOODLES` profile (SAM-spec-conformant content: noodles refuses
//! to write the rest), and the remaining differences are taken apart here:
//!
//! - Floats (`f`, `B:f`): noodles prints the shortest decimal that reads back
//!   as the same `f32`, htslib (and seqair) six significant digits. Compared
//!   by value, to six significant digits.
//! - Long CIGARs: noodles also reads the real CIGAR out of `CG`, but keeps the
//!   `CG` tag in its output; htslib drops it. The first `CG` tag is removed
//!   from noodles' tags before comparing. (noodles also takes the CIGAR from
//!   `CG` for unplaced records and for any array subtype, where htslib does
//!   not; the generator only builds long CIGARs both read alike.)
//!
//! Everything else — the eleven mandatory fields and every other tag — must
//! be identical.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "test code"
)]

#[path = "../src/pinned.rs"]
mod pinned;
#[path = "support/raw_record_gen.rs"]
mod raw_record_gen;

use hegel::TestCase;
use hegel::generators as gs;
use noodles::sam::alignment::io::Write as _;
use raw_record_gen::{GenRecord, NOODLES, arb_record};
use seqair::bam::RawBamRecord;
use seqair::bam::writer::BamWriterBuilder;
use seqair::sam::format::format_sam_record;

fn noodles_lines(path: &std::path::Path) -> Vec<Vec<u8>> {
    let mut reader = noodles::bam::io::Reader::new(std::fs::File::open(path).unwrap());
    let header = reader.read_header().unwrap();
    let mut writer = noodles::sam::io::Writer::new(Vec::new());
    for record in reader.records() {
        writer.write_alignment_record(&header, &record.unwrap()).unwrap();
    }
    let text = writer.into_inner();
    let mut lines: Vec<Vec<u8>> = text.split(|&b| b == b'\n').map(<[u8]>::to_vec).collect();
    assert_eq!(lines.pop().as_deref(), Some(&b""[..]), "noodles ends every line with a newline");
    lines
}

fn parse_f32(text: &[u8]) -> f32 {
    std::str::from_utf8(text).unwrap().parse().unwrap()
}

/// `ours` is six significant digits of the value `theirs` spells exactly.
fn assert_six_digits(ours: &[u8], theirs: &[u8]) {
    let (a, b) = (f64::from(parse_f32(ours)), f64::from(parse_f32(theirs)));
    let tolerance = b.abs() * 5e-6;
    // A rounded f32 may also land a subnormal step away.
    let tolerance = tolerance.max(f64::from(f32::from_bits(1)));
    assert!((a - b).abs() <= tolerance, "{} vs {}", ours.escape_ascii(), theirs.escape_ascii());
}

fn assert_same_tag(ours: &[u8], theirs: &[u8]) {
    let show = || format!("{} vs {}", ours.escape_ascii(), theirs.escape_ascii());
    if ours.get(2..5) == Some(b":f:") {
        assert_eq!(ours[..5], theirs[..5], "{}", show());
        assert_six_digits(&ours[5..], &theirs[5..]);
    } else if ours.get(2..7) == Some(b":B:f,") {
        assert_eq!(ours[..6], theirs[..6], "{}", show());
        let a: Vec<_> = ours[7..].split(|&b| b == b',').collect();
        let b: Vec<_> = theirs[7..].split(|&b| b == b',').collect();
        assert_eq!(a.len(), b.len(), "{}", show());
        for (a, b) in a.iter().zip(&b) {
            assert_six_digits(a, b);
        }
    } else {
        assert_eq!(ours.escape_ascii().to_string(), theirs.escape_ascii().to_string());
    }
}

// r[verify sam.format]
// r[verify sam.format.fixed_fields]
// r[verify sam.format.cigar]
// r[verify sam.format.long_cigar]
// r[verify sam.format.seq_qual]
// r[verify sam.format.aux]
#[hegel::test(test_cases = pinned::cases(300))]
fn format_agrees_with_noodles(tc: TestCase) {
    let generated: Vec<GenRecord> = tc.draw(gs::vecs(arb_record(NOODLES)).min_size(1).max_size(6));
    let header = raw_record_gen::header();
    let records: Vec<RawBamRecord> = generated
        .iter()
        .map(|g| RawBamRecord::from_payload(g.payload(), &header).unwrap())
        .collect();
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("records.bam");
    let mut writer = BamWriterBuilder::to_path(&path, &header).build().unwrap();
    for record in &records {
        writer.write_raw(record).unwrap();
    }
    writer.finish().unwrap();

    let theirs = noodles_lines(&path);
    assert_eq!(theirs.len(), records.len());
    for (record, theirs) in records.iter().zip(&theirs) {
        let mut ours = Vec::new();
        format_sam_record(record, &header, &mut ours).unwrap();
        let ours: Vec<&[u8]> = ours.split(|&b| b == b'\t').collect();
        let mut theirs: Vec<&[u8]> = theirs.split(|&b| b == b'\t').collect();
        assert!(ours.len() >= 11 && theirs.len() >= 11);
        let mandatory =
            |f: &[&[u8]]| f[..11].iter().map(|f| f.escape_ascii().to_string()).collect::<Vec<_>>();
        assert_eq!(mandatory(&ours), mandatory(&theirs));
        if record.has_cg_placeholder() {
            let cg = theirs[11..]
                .iter()
                .position(|t| t.starts_with(b"CG:B:"))
                .expect("noodles keeps CG");
            theirs.remove(11 + cg);
        }
        assert_eq!(ours.len(), theirs.len(), "same tags");
        for (a, b) in ours[11..].iter().zip(&theirs[11..]) {
            assert_same_tag(a, b);
        }
    }
}
