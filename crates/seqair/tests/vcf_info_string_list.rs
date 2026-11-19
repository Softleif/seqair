//! A list of strings in one INFO field (`r[record_encoder.string_list]`):
//! each element percent-encoded, joined with `,`, the same text in VCF and
//! BCF, as bcftools reads it.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::indexing_slicing, reason = "test code")]

use seqair::vcf::record_encoder::InfoFieldDef;
use seqair::vcf::{Alleles, ContigDef, Number, OutputFormat, Str, VcfError, VcfHeader, Writer};
use seqair_types::{Base, Pos1};
use std::process::Command;

/// One record per list, written in `format`; `None` where encoding failed.
fn write(format: OutputFormat, lists: &[&[&str]]) -> (Vec<u8>, Vec<Option<VcfError>>) {
    let mut builder = VcfHeader::builder();
    let chr1 = builder.register_contig("chr1", ContigDef::new(Some(1000))).unwrap();
    let mut builder = builder.infos();
    let sources = builder
        .register_info(&InfoFieldDef::<Str>::with_number("SOURCES", Number::Unknown, "s"))
        .unwrap();
    let one = builder
        .register_info(&InfoFieldDef::<Str>::with_number("ONE", Number::Count(1), "o"))
        .unwrap();
    let header = builder.build().unwrap();
    let mut out = Vec::new();
    let mut errors = Vec::new();
    {
        let mut writer = Writer::new(&mut out, format).write_header(&header).unwrap();
        let alleles = Alleles::snv(Base::A, Base::C).unwrap();
        for (i, list) in lists.iter().enumerate() {
            let pos = Pos1::new(u32::try_from(i).unwrap().saturating_add(1)).unwrap();
            let mut enc = writer.begin_record(&chr1, pos, &alleles, None).unwrap().filter_pass();
            sources.encode_list(&mut enc, list).unwrap();
            errors.push(one.encode_list(&mut enc, list).err());
            enc.emit().unwrap();
        }
        writer.finish().unwrap();
    }
    (out, errors)
}

// r[verify record_encoder.string_list]
// r[verify bcf_writer.string_escape]
#[test]
fn a_string_list_is_its_elements_joined() {
    let lists: &[&[&str]] = &[&["vcf_0", "vcf_1"], &["a,b", "c;d"], &[], &["x"]];
    let (text, errors) = write(OutputFormat::Vcf, lists);
    let text = String::from_utf8(text).unwrap();
    let info: Vec<&str> =
        text.lines().filter(|l| !l.starts_with('#')).filter_map(|l| l.split('\t').nth(7)).collect();
    assert_eq!(
        info,
        ["SOURCES=vcf_0,vcf_1", "SOURCES=a%2Cb,c%3Bd", "SOURCES=.", "SOURCES=x;ONE=x"]
    );
    let mismatch = |e: &Option<VcfError>| matches!(e, Some(VcfError::CountMismatch { .. }));
    assert!(errors[..3].iter().all(mismatch), "Number=1 takes one element: {errors:?}");
    assert!(errors[3].is_none());

    let Ok(bcftools) = Command::new("bcftools").arg("--version").output() else { return };
    assert!(bcftools.status.success());
    let dir = tempfile::tempdir().unwrap();
    let query = |name: &str, bytes: &[u8]| {
        let path = dir.path().join(name);
        std::fs::write(&path, bytes).unwrap();
        let q = Command::new("bcftools").args(["query", "-f", "%SOURCES\n"]).arg(&path).output();
        String::from_utf8(q.unwrap().stdout).unwrap()
    };
    let (bcf, _) = write(OutputFormat::Bcf, lists);
    let from_text = query("t.vcf", text.as_bytes());
    assert_eq!(from_text, "vcf_0,vcf_1\na%2Cb,c%3Bd\n.\nx\n");
    assert_eq!(query("t.bcf", &bcf), from_text, "BCF stores the VCF text");
}
