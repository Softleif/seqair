//! seqair writes one VCF version, and it is the one that defines what it emits.
//!
//! Covers `vcf_header.file_format` and `vcf_header.file_format.single_version`.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use seqair::vcf::{Arr, FormatFieldDef, Number, VcfHeader, VcfVersion};

const fn format_def(name: &'static str, number: Number) -> FormatFieldDef<Arr<f32>> {
    FormatFieldDef::<Arr<f32>>::with_number(name, number, "test field")
}

/// The declared version has to be one that defines `Number=M`, which is VCF
/// 4.5 and nothing earlier. The header used to default to 4.3 while offering
/// `M` — a header promising a grammar it then violates, with nothing
/// downstream able to tell.
///
/// There is no setter, so this is the only value there is; the test pins it
/// against the emitted text rather than against the constant alone, which is
/// what a reader actually sees.
// r[verify vcf_header.file_format+2]
// r[verify vcf_header.file_format.single_version+4]
#[test]
fn the_declared_version_defines_the_cardinalities_we_emit() {
    let mut builder = VcfHeader::builder().filters().infos().formats();
    builder.register_format(&format_def("M5mC", Number::BaseModification)).unwrap();
    let header = builder.samples().build().expect("valid header");

    assert_eq!(header.version(), VcfVersion::V4_5);
    assert_eq!(VcfHeader::WRITE_VERSION, VcfVersion::V4_5);
    assert!(VcfVersion::V4_5 > VcfVersion::V4_4, "versions are ordered");
    assert_eq!(VcfVersion::V4_3.to_string(), "VCFv4.3");
    assert!(
        header.to_vcf_text().starts_with("##fileformat=VCFv4.5\n"),
        "the emitted header must declare it, not just the accessor"
    );
}
