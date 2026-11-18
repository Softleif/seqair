//! The reading guide's examples are the reader's acceptance target
//! (`crates/seqair/src/vcf/reading.md`, `r[vcf_fields.guide]`): each runs as
//! a doctest of `seqair::vcf::reader`, and none may stay `ignore`.
//!
//! Until the VCF text, region query and partition phases land, the examples
//! that need them stay `ignore`. This test pins exactly which: making an
//! example run means taking it off the list, and the list must end empty.
#![allow(clippy::unwrap_used, reason = "test code")]

/// (section heading, example index within the section) of every example
/// still marked `ignore`, and what it waits for.
const STILL_IGNORED: &[(&str, usize, &str)] = &[
    ("One sample: no sample handles", 0, "VCF text (truth.vcf.gz)"),
    ("Region queries", 0, "region queries"),
    ("Region queries", 1, "query_many"),
    ("Many threads", 0, "partitions"),
    ("Owning records: `read_into`, `records()`", 0, "VCF text (rastair.vcf.gz)"),
    ("Positions only", 0, "VCF text (truth.vcf.gz)"),
    (
        "Your own record type",
        0,
        "its last assertion needs a passing CpG SNV, which rastair.bcf does not have",
    ),
];

// r[verify vcf_fields.guide]
#[test]
fn every_guide_example_runs_but_the_listed_ones() {
    let guide = include_str!("../src/vcf/reading.md");
    let mut heading = "";
    let mut index = 0;
    let mut ignored = Vec::new();
    let mut examples = 0;
    for line in guide.lines() {
        if let Some(h) = line.strip_prefix("## ") {
            heading = h;
            index = 0;
        } else if let Some(info) = line.strip_prefix("```rust") {
            examples += 1;
            if info.split(',').any(|a| a.trim() == "ignore") {
                ignored.push((heading, index));
            }
            index += 1;
        }
    }
    let expected: Vec<(&str, usize)> = STILL_IGNORED.iter().map(|&(h, i, _)| (h, i)).collect();
    assert_eq!(ignored, expected, "update STILL_IGNORED: it must only ever shrink");
    assert!(examples > ignored.len(), "some examples run");
}
