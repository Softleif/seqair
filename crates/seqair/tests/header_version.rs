//! seqair writes one VCF version, and it is the one that defines what it emits.
//!
//! Covers `vcf_header.file_format`.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use seqair::vcf::{FormatFieldDef, Number, Scalar, ValueType, VcfHeader};

const fn format_def(name: &'static str, number: Number) -> FormatFieldDef<Scalar<f32>> {
    FormatFieldDef::new(name, number, ValueType::Float, "test field")
}

/// The declared version has to be one that defines `Number=M`, which is VCF
/// 4.5 and nothing earlier. The header used to default to 4.3 while offering
/// `M` — a header promising a grammar it then violates, with nothing
/// downstream able to tell.
///
/// There is no setter, so this is the only value there is; the test pins it
/// against the emitted text rather than against the constant alone, which is
/// what a reader actually sees.
#[test]
fn the_declared_version_defines_the_cardinalities_we_emit() {
    let mut builder = VcfHeader::builder().filters().infos().formats();
    builder.register_format(&format_def("M5mC", Number::BaseModification)).unwrap();
    let header = builder.samples().build().expect("valid header");

    assert_eq!(header.file_format(), "VCFv4.5");
    assert!(
        header.to_vcf_text().starts_with("##fileformat=VCFv4.5\n"),
        "the emitted header must declare it, not just the accessor"
    );
}
