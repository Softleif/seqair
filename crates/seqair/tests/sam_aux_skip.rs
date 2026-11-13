//! `SamRecordEncoder::with_aux_tags`: unrequested optional fields are
//! skipped unparsed (`r[sam.stream.fields.aux]`).
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "test code"
)]

#[path = "support/sam_lines.rs"]
mod sam_lines;

use hegel::TestCase;
use hegel::generators as gs;
use sam_lines::{Dialect, canonical_line, header_text};
use seqair::bam::{BamHeader, RawBamRecord};
use seqair::reader::AuxTags;
use seqair::sam::{SamRecordEncoder, SamRecordError};

fn encode(line: &str, tags: AuxTags, header: &BamHeader) -> Result<RawBamRecord, SamRecordError> {
    let mut record = RawBamRecord::default();
    SamRecordEncoder::new().with_aux_tags(tags).encode(line.as_bytes(), header, &mut record)?;
    Ok(record)
}

// r[verify sam.stream.fields.aux]
/// Against the full parse, with the independent oracle of filtering its
/// decoded tags: same tags, same order, every other byte equal.
#[hegel::test]
fn a_filtered_parse_is_the_full_parse_filtered(tc: TestCase) {
    let header = BamHeader::from_sam_text(&header_text()).unwrap();
    let line = tc.draw(canonical_line(Dialect { noodles: false }));
    let wanted: Vec<[u8; 2]> = line
        .fields
        .iter()
        .skip(11)
        .filter(|_| tc.draw_silent(gs::booleans()))
        .map(|f| [f.as_bytes()[0], f.as_bytes()[1]])
        .chain([*b"ZZ"])
        .collect();
    let tags = match tc.draw_silent(gs::integers::<u8>().max_value(2)) {
        0 => AuxTags::All,
        1 => AuxTags::None,
        _ => AuxTags::only(wanted),
    };
    let text = line.text();
    let full = encode(&text, AuxTags::All, &header).unwrap();
    let part = encode(&text, tags.clone(), &header).unwrap();

    let full_aux_len = full.aux().as_bytes().len();
    let full_payload = full.clone().into_payload();
    let expected: Vec<_> = full
        .aux()
        .iter_tags()
        .filter(|(t, _)| tags.contains(*t))
        .map(|(t, v)| (t, format!("{v:?}")))
        .collect();
    let got: Vec<_> = part.aux().iter_tags().map(|(t, v)| (t, format!("{v:?}"))).collect();
    assert_eq!(got, expected);

    let fixed = full_payload.len() - full_aux_len;
    assert_eq!(part.into_payload().get(..fixed), full_payload.get(..fixed), "non-aux bytes differ");
}

// r[verify sam.stream.fields.aux]
#[test]
fn a_malformed_unrequested_tag_does_not_fail_the_record() {
    let header = BamHeader::from_sam_text(&header_text()).unwrap();
    let base = "r\t0\tchr2\t10\t60\t4M\t*\t0\t0\tACGT\tIIII";
    let nm = || AuxTags::only([*b"NM"]);
    for bad in ["XX:q:1", "XX:i:notanumber", "XX:B:A,1", "XX:i:99999999999"] {
        assert!(encode(&format!("{base}\t{bad}\tNM:i:2"), AuxTags::All, &header).is_err(), "{bad}");
        let rec = encode(&format!("{base}\t{bad}\tNM:i:2"), nm(), &header).unwrap();
        assert_eq!(rec.aux().iter_tags().count(), 1, "{bad}");
        // The same fault in a requested tag still errors.
        let line = format!("{base}\t{}", bad.replace("XX", "NM"));
        assert!(encode(&line, nm(), &header).is_err(), "{bad}");
    }
    // Required fields are still validated.
    assert!(encode("r\t0\tchr2\t10\t60\t4M\t*\t0\t0\tACGT\tII", AuxTags::None, &header).is_err());
}
