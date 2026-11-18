//! `VcfHeader::parse` against htslib's reading of the same text: real headers
//! (written by bcftools, with and without `IDX=`), and every header the
//! builder can produce, round-tripped.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "test code"
)]

#[path = "../src/pinned.rs"]
mod pinned;

use hegel::TestCase;
use hegel::generators as gs;
use seqair::vcf::{
    Attribute, ContigDef, FilterDef, FormatDef, InfoDef, Number, ValueType, VcfHeader, VcfVersion,
};
use seqair_types::SmolStr;
use std::collections::BTreeMap;
use std::io::{Read, Write};
use std::path::Path;
use std::process::Command;

const FIXTURES: [&str; 5] = [
    "clinvar.vcf.header",
    "giab_hg002.vcf.header",
    "gnomad41_chr22.bcf.header",
    "kgp_chr22.bcf.header",
    "rastair_sv168.bcf.header",
];

fn fixture(name: &str) -> String {
    let path =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/vcf/headers").join(name);
    std::fs::read_to_string(path).unwrap()
}

type Dictionary = BTreeMap<String, usize>;

/// The string and contig dictionaries bcftools assigns to `text`: convert it to
/// BCF, where htslib writes `IDX=` on every definition line, and read those
/// back with plain string search — not with the parser under test.
fn bcftools_dictionaries(text: &str) -> (Dictionary, Dictionary) {
    let mut vcf = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
    vcf.write_all(text.as_bytes()).unwrap();
    vcf.flush().unwrap();
    let bcf = tempfile::Builder::new().suffix(".bcf").tempfile().unwrap();
    let out = Command::new("bcftools")
        .args(["view", "--no-version", "-Ob", "-o"])
        .arg(bcf.path())
        .arg(vcf.path())
        .output()
        .expect("bcftools");
    assert!(out.status.success(), "bcftools: {}", String::from_utf8_lossy(&out.stderr));

    let mut bytes = Vec::new();
    noodles_bgzf::io::Reader::new(std::fs::File::open(bcf.path()).unwrap())
        .read_to_end(&mut bytes)
        .unwrap();
    assert_eq!(bytes.get(..5), Some(&b"BCF\x02\x02"[..]));
    let l_text = u32::from_le_bytes(bytes[5..9].try_into().unwrap()) as usize;
    let l_text = std::str::from_utf8(&bytes[9..9 + l_text]).unwrap();

    let (mut fields, mut contigs) = (Dictionary::new(), Dictionary::new());
    for line in l_text.lines() {
        let dict = if ["##FILTER=<", "##INFO=<", "##FORMAT=<"].iter().any(|p| line.starts_with(p)) {
            &mut fields
        } else if line.starts_with("##contig=<") {
            &mut contigs
        } else {
            continue;
        };
        let id = line.split_once("<ID=").unwrap().1;
        let id = &id[..id.find([',', '>']).unwrap()];
        let idx = &line[line.rfind(",IDX=").unwrap() + 5..line.len() - 1];
        let previous = dict.insert(id.to_owned(), idx.parse().unwrap());
        assert!(previous.is_none_or(|p| p == idx.parse().unwrap()), "{id} twice in {line}");
    }
    (fields, contigs)
}

fn dictionaries(header: &VcfHeader) -> (Dictionary, Dictionary) {
    let collect = |map: &seqair::vcf::header::StringMap| {
        map.iter().map(|(idx, id)| (id.to_owned(), idx)).collect()
    };
    let contigs = header.contigs().map(|(c, _)| (c.name().to_string(), c.index())).collect();
    (collect(header.string_map()), contigs)
}

/// Everything a header holds, in order (an `IndexMap`'s `==` ignores order).
#[derive(Debug, PartialEq)]
struct Snapshot {
    version: VcfVersion,
    filters: Vec<(SmolStr, FilterDef)>,
    infos: Vec<(SmolStr, InfoDef)>,
    formats: Vec<(SmolStr, FormatDef)>,
    contigs: Vec<(SmolStr, ContigDef)>,
    samples: Vec<SmolStr>,
    other_lines: Vec<SmolStr>,
    dictionaries: (Dictionary, Dictionary),
}

fn snapshot(h: &VcfHeader) -> Snapshot {
    Snapshot {
        version: h.version(),
        filters: h.filters().clone().into_iter().collect(),
        infos: h.infos().clone().into_iter().collect(),
        formats: h.formats().clone().into_iter().collect(),
        contigs: h.contigs().map(|(c, def)| (c.name().clone(), def.clone())).collect(),
        samples: h.samples().to_vec(),
        other_lines: h.other_lines().to_vec(),
        dictionaries: dictionaries(h),
    }
}

// r[verify vcf_header.parse.round_trip]
// r[verify vcf_header.parse.dictionary]
/// Real headers get the dictionary htslib gives them: the BCF ones by their
/// `IDX=` (with the gaps a trimmed gnomAD header leaves, and 1000 Genomes'
/// FORMAT GT between INFO lines), the VCF ones by order of appearance.
#[test]
fn real_headers_get_the_dictionaries_bcftools_assigns() {
    for name in FIXTURES {
        let text = fixture(name);
        let header = VcfHeader::parse(&text).unwrap_or_else(|e| panic!("{name}: {e}"));
        assert_eq!(dictionaries(&header), bcftools_dictionaries(&text), "{name}");
    }
}

// r[verify vcf_header.parse.preserve]
// r[verify vcf_header.parse.round_trip]
/// Writing a parsed header back out loses nothing a re-read sees, and where
/// its dictionary is not the emission order, the `IDX=` it writes give
/// htslib the same dictionary too.
#[test]
fn real_headers_round_trip() {
    for name in FIXTURES {
        let header = VcfHeader::parse(&fixture(name)).unwrap();
        let text = header.to_vcf_text();
        let again = VcfHeader::parse(&text).unwrap_or_else(|e| panic!("{name}: {e}\n{text}"));
        assert_eq!(snapshot(&again), snapshot(&header), "{name}");
        assert_eq!(bcftools_dictionaries(&text), dictionaries(&header), "{name}");
        assert_eq!(again.to_vcf_text(), text, "{name}: a second round trip is a fixed point");

        let reordered = matches!(name, "gnomad41_chr22.bcf.header" | "kgp_chr22.bcf.header");
        assert_eq!(text.contains(",IDX="), reordered, "{name}: IDX only where needed");
    }
}

// r[verify vcf_header.parse.version]
// r[verify vcf_header.file_format.single_version+2]
#[test]
fn a_parsed_header_keeps_its_version_and_attributes() {
    let rastair = VcfHeader::parse(&fixture("rastair_sv168.bcf.header")).unwrap();
    assert_eq!(rastair.version(), VcfVersion::V4_3);
    assert!(rastair.to_vcf_text().starts_with("##fileformat=VCFv4.3\n"));
    assert_eq!(rastair.formats()["M5mC"].number, Number::BaseModification);
    assert!(rastair.other_lines().iter().any(|l| l.starts_with("rastairConfig={")));

    let giab = VcfHeader::parse(&fixture("giab_hg002.vcf.header")).unwrap();
    let (_, chr1) = giab.contigs().find(|(c, _)| c.name() == "chr1").unwrap();
    assert_eq!(chr1.length, Some(248_956_422));
    assert_eq!(
        chr1.attributes.first().map(|a| (a.key.as_str(), a.value.as_str())),
        Some(("assembly", "human_GRCh38_no_alt_analysis_set.fasta"))
    );
    assert_eq!(giab.samples(), ["HG002"]);

    let kgp = VcfHeader::parse(&fixture("kgp_chr22.bcf.header")).unwrap();
    assert_eq!(kgp.version(), VcfVersion::V4_2);
    assert_eq!(kgp.string_map().get("GT"), Some(4), "FORMAT GT sits between INFO lines");
    assert_eq!(kgp.samples().len(), 4);
}

// ── Builder headers, round-tripped ─────────────────────────────────────

const PRINTABLE: &str = r#"[ -~é]{0,16}"#;

#[hegel::composite]
fn attributes(tc: &TestCase) -> Vec<Attribute> {
    let n = tc.draw(gs::integers::<usize>().max_value(2));
    (0..n)
        .map(|i| {
            let quoted = tc.draw(gs::booleans());
            let value = if quoted {
                tc.draw(gs::from_regex(PRINTABLE).fullmatch(true))
            } else {
                tc.draw(gs::from_regex("[A-Za-z0-9_.:/]{1,10}").fullmatch(true))
            };
            Attribute { key: format!("x{i}").into(), value: value.into(), quoted }
        })
        .collect()
}

#[hegel::composite]
fn number(tc: &TestCase) -> Number {
    match tc.draw(gs::integers::<u8>().max_value(10)) {
        0 => Number::AlternateBases,
        1 => Number::ReferenceAlternateBases,
        2 => Number::Genotypes,
        3 => Number::BaseModification,
        4 => Number::Ploidy,
        5 => Number::LocalAlternateBases,
        6 => Number::LocalReferenceAlternateBases,
        7 => Number::LocalGenotypes,
        8 => Number::Unknown,
        _ => Number::Count(tc.draw(gs::integers::<u32>().max_value(1000))),
    }
}

/// A header in the builder's shape: ids made unique by an index suffix, each
/// declared in a random non-empty subset of FILTER/INFO/FORMAT, so INFO and
/// FORMAT (and FILTER) share ids.
#[hegel::composite]
fn builder_header(tc: &TestCase) -> VcfHeader {
    let types = [
        ValueType::Integer,
        ValueType::Float,
        ValueType::Flag,
        ValueType::Character,
        ValueType::String,
    ];
    let n_contigs = tc.draw(gs::integers::<usize>().max_value(6));
    let n_ids = tc.draw(gs::integers::<usize>().max_value(12));
    let n_samples = tc.draw(gs::integers::<usize>().max_value(5));
    let n_other = tc.draw(gs::integers::<usize>().max_value(3));

    let mut builder = VcfHeader::builder();
    for i in 0..n_other {
        let value = tc.draw(gs::from_regex("[!-;=-~][ -~]{0,15}").fullmatch(true));
        builder.add_other_line(format!("meta{i}={value}"));
    }
    for i in 0..n_contigs {
        let prefix = tc.draw(gs::from_regex("[A-Za-z0-9_.]{0,6}").fullmatch(true));
        let length = tc.draw(gs::optional(gs::integers::<u64>()));
        let mut def = ContigDef::new(length);
        def.attributes = tc.draw(attributes());
        builder.add_contig(format!("{prefix}c{i}"), def).unwrap();
    }
    let ids: Vec<(String, u8)> = (0..n_ids)
        .map(|i| {
            let prefix = tc.draw(gs::from_regex("[A-Za-z_][A-Za-z0-9_.]{0,5}").fullmatch(true));
            (format!("{prefix}{i}"), tc.draw(gs::integers::<u8>().min_value(1).max_value(7)))
        })
        .collect();
    let description = |tc: &TestCase| tc.draw(gs::from_regex(PRINTABLE).fullmatch(true));

    let mut builder = builder.filters();
    for (id, _) in ids.iter().filter(|(_, mask)| mask & 1 != 0) {
        let mut def = FilterDef::new(description(tc));
        def.attributes = tc.draw(attributes());
        builder.add_filter(id.as_str(), def).unwrap();
    }
    let mut builder = builder.infos();
    for (id, _) in ids.iter().filter(|(_, mask)| mask & 2 != 0) {
        let typ = tc.draw(gs::sampled_from(types.to_vec()));
        let number = if typ == ValueType::Flag { Number::Count(0) } else { tc.draw(number()) };
        let mut def = InfoDef::new(number, typ, description(tc));
        def.attributes = tc.draw(attributes());
        builder.add_info(id.as_str(), def).unwrap();
    }
    let mut builder = builder.formats();
    for (id, _) in ids.iter().filter(|(_, mask)| mask & 4 != 0) {
        let typ = tc.draw(gs::sampled_from(
            types.iter().copied().filter(|t| *t != ValueType::Flag).collect::<Vec<_>>(),
        ));
        let mut def = FormatDef::new(tc.draw(number()), typ, description(tc));
        def.attributes = tc.draw(attributes());
        builder.add_format(id.as_str(), def).unwrap();
    }
    let mut builder = builder.samples();
    for i in 0..n_samples {
        let prefix = tc.draw(gs::from_regex("[A-Za-z0-9_.-]{0,6}").fullmatch(true));
        builder.add_sample(format!("{prefix}s{i}")).unwrap();
    }
    builder.build().unwrap()
}

// r[verify vcf_header.parse.round_trip]
#[hegel::test]
fn every_builder_header_parses_back_to_itself(tc: TestCase) {
    let header = tc.draw(builder_header());
    let text = header.to_vcf_text();
    let parsed = VcfHeader::parse(&text).unwrap_or_else(|e| panic!("{e}\n{text}"));
    assert_eq!(snapshot(&parsed), snapshot(&header), "{text}");
    assert!(!text.contains(",IDX="), "a built header's dictionary is its emission order");
    assert_ne!(parsed.id(), header.id(), "a parsed header is a header of its own");
}

// r[verify vcf_header.parse.round_trip]
// r[verify vcf_header.parse.dictionary]
/// bcftools assigns a built header the dictionary seqair encodes records with.
#[hegel::test(test_cases = pinned::cases(30))]
fn bcftools_assigns_a_builder_header_the_same_dictionary(tc: TestCase) {
    let header = tc.draw(builder_header());
    assert_eq!(bcftools_dictionaries(&header.to_vcf_text()), dictionaries(&header));
}
