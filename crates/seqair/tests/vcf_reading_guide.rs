//! The reading guide's examples are the reader's acceptance target
//! (`crates/seqair/src/vcf/reading.md`, `r[vcf_fields.guide]`): each runs as
//! a doctest of `seqair::vcf::reader`, and none may stay `ignore`.
//!
//! The list of examples allowed to stay `ignore` is empty and must stay so.
#![allow(clippy::unwrap_used, reason = "test code")]

/// (section heading, example index within the section) of every example
/// still marked `ignore`, and what it waits for.
const STILL_IGNORED: &[(&str, usize, &str)] = &[];

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
