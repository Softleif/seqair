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
    let vcf = vcf_file(text);
    let bcf = tempfile::Builder::new().suffix(".bcf").tempfile().unwrap();
    let out = Command::new("bcftools")
        .args(["view", "--no-version", "-Ob", "-o"])
        .arg(bcf.path())
        .arg(vcf.path())
        .output()
        .expect("bcftools");
    assert!(out.status.success(), "bcftools: {}", String::from_utf8_lossy(&out.stderr));
    l_text_dictionaries(bcf.path())
}

/// `text` in a temporary `.vcf` file.
fn vcf_file(text: &str) -> tempfile::NamedTempFile {
    let mut vcf = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
    vcf.write_all(text.as_bytes()).unwrap();
    vcf.flush().unwrap();
    vcf
}

/// The dictionaries a BCF file's `l_text` declares, every definition line
/// carrying the `IDX=` htslib writes.
fn l_text_dictionaries(bcf: &Path) -> (Dictionary, Dictionary) {
    let mut bytes = Vec::new();
    noodles_bgzf::io::Reader::new(std::fs::File::open(bcf).unwrap())
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
        // Where the line had it, or last when htslib added it.
        let idx = &line[line.rfind(",IDX=").unwrap() + 5..];
        let idx = &idx[..idx.find(|c: char| !c.is_ascii_digit()).unwrap()];
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

// r[verify vcf_header.parse.round_trip+2]
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

// r[verify vcf_header.parse.preserve+2]
// r[verify vcf_header.parse.round_trip+2]
// r[verify vcf_header.serialization+2]
/// Writing a parsed header back out loses nothing a re-read sees, and where
/// its dictionary is not its line order, the `IDX=` its BCF text writes give
/// htslib the same dictionary too. Its VCF text has no `IDX`, and loses
/// nothing but the dictionary's holes.
#[test]
fn real_headers_round_trip() {
    for name in FIXTURES {
        let header = VcfHeader::parse(&fixture(name)).unwrap();
        let text = header.to_bcf_text();
        let again = VcfHeader::parse(&text).unwrap_or_else(|e| panic!("{name}: {e}\n{text}"));
        assert_eq!(snapshot(&again), snapshot(&header), "{name}");
        assert_eq!(bcftools_dictionaries(&text), dictionaries(&header), "{name}");
        assert_eq!(again.to_bcf_text(), text, "{name}: a second round trip is a fixed point");

        // kgp declares FORMAT GT between INFO lines: its line order is its dictionary.
        let gaps = name == "gnomad41_chr22.bcf.header";
        assert_eq!(text.contains(",IDX="), gaps, "{name}: IDX only where needed");

        let vcf = header.to_vcf_text();
        assert!(!vcf.contains("IDX="), "{name}");
        let again = VcfHeader::parse(&vcf).unwrap();
        let (mut ours, mut theirs) = (snapshot(&again), snapshot(&header));
        assert_eq!(ours.dictionaries, bcftools_dictionaries(&vcf), "{name}");
        (ours.dictionaries, theirs.dictionaries) = Default::default();
        assert_eq!(ours, theirs, "{name}: the VCF text loses no line");
        assert_eq!(again.to_vcf_text(), vcf, "{name}");
    }
}

// r[verify vcf_header.parse.version]
// r[verify vcf_header.file_format.single_version+4]
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

/// Printable text that does not end in a space: a quoted value's trailing
/// spaces are dropped on reading, by htslib and by `VcfHeader::parse`.
const PRINTABLE: &str = r#"([ -~é]{0,15}[!-~é])?"#;

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

// r[verify vcf_header.parse.round_trip+2]
#[hegel::test]
fn every_builder_header_parses_back_to_itself(tc: TestCase) {
    let header = tc.draw(builder_header());
    let text = header.to_vcf_text();
    let parsed = VcfHeader::parse(&text).unwrap_or_else(|e| panic!("{e}\n{text}"));
    assert_eq!(snapshot(&parsed), snapshot(&header), "{text}");
    assert_eq!(header.to_bcf_text(), text, "a built header's dictionary is its emission order");
    assert_ne!(parsed.id(), header.id(), "a parsed header is a header of its own");
}

// r[verify vcf_header.parse.round_trip+2]
// r[verify vcf_header.parse.dictionary]
/// bcftools assigns a built header the dictionary seqair encodes records with.
#[hegel::test(test_cases = pinned::cases(30))]
fn bcftools_assigns_a_builder_header_the_same_dictionary(tc: TestCase) {
    let header = tc.draw(builder_header());
    assert_eq!(bcftools_dictionaries(&header.to_bcf_text()), dictionaries(&header));
}

// ── Derived headers ────────────────────────────────────────────────────

// r[verify vcf_header.to_builder+2]
// r[verify vcf_header.brand+3]
/// `to_builder().build()` is the header again, under an id of its own: built
/// headers, and real ones whose `IDX=` leave gaps and reorder the dictionary.
#[hegel::test]
fn a_derived_header_is_its_source_but_for_the_id(tc: TestCase) {
    let header = tc.draw(builder_header());
    let derived = header.to_builder().build().unwrap();
    assert_eq!(snapshot(&derived), snapshot(&header));
    assert_eq!(derived.to_vcf_text(), header.to_vcf_text());
    assert_ne!(derived.id(), header.id());
}

// r[verify vcf_header.to_builder+2]
#[test]
fn a_derived_real_header_is_its_source_but_for_the_id() {
    for name in FIXTURES {
        let header = VcfHeader::parse(&fixture(name)).unwrap();
        let derived = header.to_builder().build().unwrap();
        assert_eq!(snapshot(&derived), snapshot(&header), "{name}");
        assert_eq!(derived.string_map(), header.string_map(), "{name}: holes included");
        assert_eq!(derived.to_vcf_text(), header.to_vcf_text(), "{name}");
        assert_ne!(derived.id(), header.id(), "{name}");
    }
}

// r[verify vcf_header.remove+2]
// r[verify vcf_header.to_builder+2]
/// Removing definitions leaves every other id at its index and the removed
/// ones' slots empty unless another category still declares them; new ids
/// go past every slot, and a removed id registered again gets its slot back,
/// as `bcftools annotate -x INFO/X -h` (with X's line) gives it. bcftools,
/// reading the written header, assigns the same dictionary.
#[hegel::test(test_cases = pinned::cases(40))]
fn removing_definitions_keeps_every_other_index(tc: TestCase) {
    let header = tc.draw(builder_header());
    let mut builder = header.to_builder();
    let mut gone = std::collections::BTreeSet::new();
    for (kind, ids) in [
        (0, header.filters().keys().cloned().collect::<Vec<_>>()),
        (1, header.infos().keys().cloned().collect()),
        (2, header.formats().keys().cloned().collect()),
    ] {
        for id in ids {
            if !tc.draw(gs::booleans()) {
                continue;
            }
            let removed = match kind {
                0 => builder.remove_filter(&id),
                1 => builder.remove_info(&id),
                _ => builder.remove_format(&id),
            };
            assert_eq!(removed, id != "PASS", "{id}");
            gone.insert((kind, id));
        }
    }
    let declared = |id: &str| {
        let still = |kind, present: bool| present && !gone.contains(&(kind, SmolStr::from(id)));
        still(0, header.filters().contains_key(id))
            || still(1, header.infos().contains_key(id))
            || still(2, header.formats().contains_key(id))
            || id == "PASS"
    };
    let (fields, contigs) = dictionaries(&header);
    let expected: Dictionary = fields.into_iter().filter(|(id, _)| declared(id)).collect();

    let slots = header.string_map().len();
    let mut builder = builder.filters();
    builder.add_filter("zzNewFilter", FilterDef::new("new")).unwrap();
    let mut builder = builder.infos();
    let def = InfoDef::new(Number::Count(1), ValueType::Integer, "new");
    builder.add_info("zzNewInfo", def.clone()).unwrap();
    let again = gone.iter().find(|(kind, _)| *kind == 1).map(|(_, id)| id.clone());
    if let Some(id) = &again {
        builder.add_info(id.as_str(), def).unwrap();
    }
    let derived = builder.build().unwrap();

    let mut expected = expected;
    expected.insert("zzNewFilter".to_owned(), slots);
    expected.insert("zzNewInfo".to_owned(), slots + 1);
    if let Some(id) = again {
        expected.insert(id.to_string(), header.string_map().get(&id).unwrap());
    }
    assert_eq!(dictionaries(&derived), (expected, contigs));
    let text = derived.to_bcf_text();
    assert_eq!(bcftools_dictionaries(&text), dictionaries(&derived), "{text}");
    let again = VcfHeader::parse(&text).unwrap();
    assert_eq!(snapshot(&again), snapshot(&derived), "{text}");
}

// r[verify vcf_header.remove+2]
#[test]
fn pass_stays_and_a_missing_id_is_reported() {
    let mut builder = VcfHeader::builder().build().unwrap().to_builder();
    assert!(!builder.remove_filter("PASS"));
    assert!(!builder.remove_info("DP"));
    assert!(!builder.remove_format("DP"));
    let header = builder.build().unwrap();
    assert_eq!(header.string_map().get("PASS"), Some(0));
    assert!(header.filters().contains_key("PASS"));
}

// r[verify vcf_header.set_samples]
#[test]
fn set_samples_replaces_them_or_leaves_them() {
    let mut builder = VcfHeader::builder().samples();
    builder.add_sample("a").unwrap();
    builder.set_samples(["c", "b"]).unwrap();
    let err = builder.set_samples(["x", "y", "x"]).unwrap_err();
    assert!(
        matches!(err, seqair::vcf::VcfHeaderError::DuplicateSample { ref name } if name == "x"),
        "{err:?}"
    );
    assert_eq!(builder.build().unwrap().samples(), ["c", "b"]);
}

// r[verify vcf_header.to_builder+2]
// r[verify vcf_header.remove+2]
/// `HiPhase`'s case: copy the input's header, swap its FORMAT PS for its own
/// PS and PF, add its `##HiPhase_*` lines, and write records with keys of
/// copied and new fields alike. bcftools reads them from VCF text and BCF.
#[test]
fn a_phaser_rewrites_ps_in_a_copied_header() {
    use seqair::vcf::record_encoder::{Arr, FormatFieldDef, Gt, Scalar, Str};
    use seqair::vcf::{Alleles, Genotype, OutputFormat, Writer};
    use seqair_types::{Base, Pos1};

    let input = VcfHeader::parse(&fixture("giab_hg002.vcf.header")).unwrap();
    let mut builder = input.to_builder();
    assert!(builder.remove_format("PS"));
    assert!(!builder.remove_format("PF"), "GIAB declares no PF");
    builder.add_other_line("HiPhase_version=1.4.5");
    builder.add_other_line("HiPhase_command=hiphase --vcf in.vcf.gz");
    let mut builder = builder.formats();
    let ps = builder
        .register_format(&FormatFieldDef::<Scalar<i32>>::new("PS", "Phase set identifier"))
        .unwrap();
    let pf = builder
        .register_format(&FormatFieldDef::<Str>::new("PF", Number::Count(1), "Phasing flag"))
        .unwrap();
    let header = builder.build().unwrap();
    assert!(header.other_lines().iter().any(|l| l == "HiPhase_version=1.4.5"));
    assert_eq!(header.version(), input.version());

    let gt = header.format_key::<Gt>("GT").unwrap();
    let ad = header.format_key::<Arr<i32>>("AD").unwrap();
    let chr1 = header.contig("chr1").unwrap();
    for (format, suffix) in [(OutputFormat::Vcf, ".vcf"), (OutputFormat::Bcf, ".bcf")] {
        let mut file = tempfile::Builder::new().suffix(suffix).tempfile().unwrap();
        let mut writer = Writer::new(file.as_file_mut(), format).write_header(&header).unwrap();
        for (pos, set) in [(1000, 1000), (1500, 1000)] {
            let alleles = Alleles::snv(Base::A, Base::G).unwrap();
            let mut enc = writer
                .begin_record(&chr1, Pos1::new(pos).unwrap(), &alleles, None)
                .unwrap()
                .filter_pass()
                .begin_samples();
            gt.encode(&mut enc, &[Genotype::phased_diploid(0, 1)]).unwrap();
            ad.encode(&mut enc, &[&[12, 9]]).unwrap();
            ps.encode(&mut enc, &[set]).unwrap();
            pf.encode(&mut enc, &["HP"]).unwrap();
            enc.emit().unwrap();
        }
        drop(writer.finish().unwrap());
        let out = Command::new("bcftools")
            .args(["query", "-f", "%POS [%GT %AD %PS %PF]\n"])
            .arg(file.path())
            .output()
            .unwrap();
        assert!(out.status.success(), "{suffix}: {}", String::from_utf8_lossy(&out.stderr));
        assert_eq!(
            String::from_utf8(out.stdout).unwrap(),
            "1000 0|1 12,9 1000 HP\n1500 0|1 12,9 1000 HP\n",
            "{suffix}"
        );
        let out = Command::new("bcftools").args(["view", "-h"]).arg(file.path()).output().unwrap();
        let text = String::from_utf8(out.stdout).unwrap();
        assert!(text.contains("##HiPhase_version=1.4.5"), "{suffix}");
        assert!(text.contains(
            "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set identifier\""
        ));
    }
}

// r[verify vcf_header.file_format.single_version+4]
/// A header derived from a VCF 4.2 or 4.3 file refuses a `Number` its
/// version does not define, naming both versions, and keeps its own version;
/// one from a 4.4 file takes `P` but not `M`; a built (4.5) header takes all.
#[test]
fn a_derived_header_refuses_a_number_newer_than_its_version() {
    use seqair::vcf::{VcfHeaderError, VcfVersion as V};
    let newer = [
        (Number::Ploidy, V::V4_4),
        (Number::BaseModification, V::V4_5),
        (Number::LocalAlternateBases, V::V4_5),
        (Number::LocalReferenceAlternateBases, V::V4_5),
        (Number::LocalGenotypes, V::V4_5),
    ];
    for (source, version) in
        [("kgp_chr22.bcf.header", V::V4_2), ("rastair_sv168.bcf.header", V::V4_3)]
    {
        let input = VcfHeader::parse(&fixture(source)).unwrap();
        assert_eq!(input.version(), version);
        for (number, needs) in newer {
            let mut builder = input.to_builder().infos();
            let err = builder
                .add_info("NEWI", InfoDef::new(number, ValueType::Integer, "x"))
                .unwrap_err();
            assert!(
                matches!(&err, VcfHeaderError::NumberNeedsNewerVersion { id, number: n, version: v, needs: w }
                    if id == "NEWI" && *n == number && *v == version && *w == needs),
                "{source} {number}: {err:?}"
            );
            let mut builder = builder.formats();
            let def = FormatDef::new(number, ValueType::Integer, "x");
            assert!(builder.add_format("NEWF", def).is_err(), "{source} {number}");
            let header = builder.build().unwrap();
            assert_eq!(header.version(), version, "the version is not raised");
            assert!(!header.infos().contains_key("NEWI") && !header.formats().contains_key("NEWF"));
            assert_eq!(
                header.string_map(),
                input.string_map(),
                "a refusal leaves the builder as it was"
            );
        }
        let mut builder = input.to_builder().formats();
        builder
            .add_format("NEWF", FormatDef::new(Number::Genotypes, ValueType::Integer, "x"))
            .unwrap();
    }

    let v44 =
        VcfHeader::parse("##fileformat=VCFv4.4\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            .unwrap();
    let mut builder = v44.to_builder().formats();
    builder.add_format("P1", FormatDef::new(Number::Ploidy, ValueType::Integer, "x")).unwrap();
    assert!(
        builder
            .add_format("M1", FormatDef::new(Number::BaseModification, ValueType::Integer, "x"))
            .is_err()
    );

    let mut built = VcfHeader::builder().formats();
    for (i, (number, _)) in newer.into_iter().enumerate() {
        built.add_format(format!("F{i}"), FormatDef::new(number, ValueType::Integer, "x")).unwrap();
    }
}

// ── Line order and form, against bcftools ──────────────────────────────

/// `bcftools <args> <path>`'s standard output.
fn bcftools(args: &[&str], path: &Path) -> String {
    let out = Command::new("bcftools").args(args).arg(path).output().expect("bcftools");
    assert!(out.status.success(), "bcftools {args:?}: {}", String::from_utf8_lossy(&out.stderr));
    String::from_utf8(out.stdout).unwrap()
}

// r[verify vcf_header.line_order]
// r[verify vcf_header.line_form]
// r[verify vcf_header.serialization+2]
/// A real file's header, read and written back, is what `bcftools view -h`
/// prints for it: same lines, same order, same attributes as written. The
/// fixtures include a BCF `l_text` with `IDX=` gaps, FORMAT GT between INFO
/// lines, and 1000 Genomes' `Type` before `Number`; the files, a GIAB truth
/// set with quoted contig attributes.
#[test]
fn real_headers_print_as_bcftools_prints_them() {
    for name in FIXTURES {
        let path =
            Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/vcf/headers").join(name);
        let ours = VcfHeader::parse(&fixture(name)).unwrap().to_vcf_text();
        assert_eq!(ours, bcftools(&["view", "-h", "--no-version"], &path), "{name}");
    }
    for name in ["rastair.bcf", "rastair.vcf.gz", "truth.vcf.gz", "panel.bcf"] {
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/vcf").join(name);
        let ours = seqair::vcf::VcfFile::open(&path).unwrap().header().to_vcf_text();
        assert_eq!(ours, bcftools(&["view", "-h", "--no-version"], &path), "{name}");
    }
}

fn shuffled<T>(tc: &TestCase, items: Vec<T>) -> Vec<T> {
    let mut keyed: Vec<(u32, T)> =
        items.into_iter().map(|i| (tc.draw(gs::integers::<u32>()), i)).collect();
    keyed.sort_by_key(|(k, _)| *k);
    keyed.into_iter().map(|(_, i)| i).collect()
}

/// `n` distinct indices from `from` up, with gaps, in random order: what a
/// BCF header's `IDX=` may hold.
fn indices(tc: &TestCase, n: usize, from: usize) -> Vec<usize> {
    let mut next = from;
    let increasing = (0..n)
        .map(|_| {
            next += tc.draw(gs::integers::<usize>().max_value(2));
            next += 1;
            next - 1
        })
        .collect();
    shuffled(tc, increasing)
}

/// A quoted value as a file writes it, with an escaped quote now and then.
fn quoted(tc: &TestCase) -> String {
    let text = tc.draw(gs::from_regex(r#"[ !#-\[\]-~]{0,10}"#).fullmatch(true));
    let escape = if tc.draw(gs::booleans()) { r#"\"x\\"# } else { "" };
    format!("\"{text}{escape}\"")
}

/// A header as files hold them: FILTER, INFO, FORMAT, contig and other lines
/// interleaved, attributes in any order after `ID`, values quoted or bare,
/// spaces after commas, `IDX=` with gaps (or none), a `PASS` line anywhere,
/// repeated other lines.
#[hegel::composite]
fn file_header(tc: &TestCase) -> String {
    let n_ids = tc.draw(gs::integers::<usize>().max_value(8));
    let ids: Vec<(String, u8)> = (0..n_ids)
        .map(|i| {
            let prefix = tc.draw(gs::from_regex("[A-Za-z_][A-Za-z0-9_]{0,3}").fullmatch(true));
            (format!("{prefix}{i}"), tc.draw(gs::integers::<u8>().min_value(1).max_value(7)))
        })
        .collect();
    let n_contigs = tc.draw(gs::integers::<usize>().max_value(4));
    let with_idx = tc.draw(gs::booleans());
    let field_idx = indices(tc, n_ids, 1);
    let contig_idx = indices(tc, n_contigs, 0);
    let idx = |i: usize| with_idx.then(|| format!("IDX={i}"));

    let extra = |tc: &TestCase| -> Vec<String> {
        let n = tc.draw(gs::integers::<usize>().max_value(2));
        (0..n)
            .map(|j| {
                let value = if tc.draw(gs::booleans()) {
                    quoted(tc)
                } else {
                    tc.draw(gs::from_regex("[A-Za-z0-9_.:/]{1,8}").fullmatch(true))
                };
                format!("x{j}={value}")
            })
            .collect()
    };
    // One line: ID first, the rest (IDX among them) in any order.
    let line = |tc: &TestCase, key: &str, id: &str, mut attrs: Vec<String>| {
        attrs.extend(extra(tc));
        let mut attrs = shuffled(tc, attrs);
        attrs.insert(0, format!("ID={id}"));
        let sep = if tc.draw(gs::booleans()) { ", " } else { "," };
        format!("##{key}=<{}>", attrs.join(sep))
    };

    let mut lines = Vec::new();
    for ((id, mask), &i) in ids.iter().zip(&field_idx) {
        if mask & 1 != 0 {
            let attrs = vec![format!("Description={}", quoted(tc))];
            lines.push(line(tc, "FILTER", id, attrs.into_iter().chain(idx(i)).collect()));
        }
        for (bit, key) in [(2, "INFO"), (4, "FORMAT")] {
            if mask & bit == 0 {
                continue;
            }
            let (number, typ) = tc.draw(gs::sampled_from(vec![
                ("1", "Integer"),
                ("A", "Float"),
                ("R", "Integer"),
                ("G", "Float"),
                (".", "String"),
                ("2", "Character"),
            ]));
            let attrs = vec![
                format!("Number={number}"),
                format!("Type={typ}"),
                format!("Description={}", quoted(tc)),
            ];
            lines.push(line(tc, key, id, attrs.into_iter().chain(idx(i)).collect()));
        }
    }
    for (c, &i) in contig_idx.iter().enumerate() {
        let mut attrs: Vec<String> = idx(i).into_iter().collect();
        match tc.draw(gs::integers::<u8>().max_value(2)) {
            0 => {}
            1 => attrs.push(format!("length={}", tc.draw(gs::integers::<u32>()))),
            _ => attrs.push(format!("length=\"{}\"", tc.draw(gs::integers::<u32>()))),
        }
        lines.push(line(tc, "contig", &format!("chr{c}"), attrs));
    }
    if tc.draw(gs::booleans()) {
        let attrs = vec![format!("Description={}", quoted(tc))];
        lines.push(line(tc, "FILTER", "PASS", attrs.into_iter().chain(idx(0)).collect()));
    }
    let n_other = tc.draw(gs::integers::<usize>().max_value(4));
    for i in 0..n_other {
        let other = match tc.draw(gs::integers::<u8>().max_value(3)) {
            0 => format!(
                "##meta{i}={}",
                tc.draw(gs::from_regex("[!-;=-~][ -~]{0,10}").fullmatch(true))
            ),
            1 => line(tc, "ALT", &format!("<DEL{i}>"), vec![format!("Description={}", quoted(tc))]),
            2 => format!("##META=<Type={}>", tc.draw(gs::from_regex("[a-z]{1,4}").fullmatch(true))),
            // A repeat: of a line drawn before, or of the version line.
            _ => match lines.last() {
                Some(last) if !last.starts_with("##FILTER") && !last.starts_with("##INFO") => {
                    last.clone()
                }
                _ => "##fileformat=VCFv4.2".to_owned(),
            },
        };
        lines.push(other);
    }
    let lines = shuffled(tc, lines);
    let n_samples = tc.draw(gs::integers::<usize>().max_value(2));
    let samples: Vec<String> = (0..n_samples).map(|i| format!("s{i}")).collect();
    let samples: String = std::iter::once("FORMAT".to_owned())
        .filter(|_| n_samples > 0)
        .chain(samples)
        .map(|column| "\t".to_owned() + &column)
        .collect();
    format!(
        "##fileformat=VCFv4.2\n{}\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO{samples}\n",
        lines.join("\n")
    )
}

// r[verify vcf_header.line_order]
// r[verify vcf_header.line_form]
// r[verify vcf_header.serialization+2]
// r[verify vcf_header.parse.round_trip+2]
/// A header as files hold it prints as `bcftools view -h` prints it, gets the
/// dictionary htslib gives it, and its BCF text gives htslib that dictionary
/// again.
#[hegel::test(test_cases = pinned::cases(40))]
fn file_headers_print_as_bcftools_prints_them(tc: TestCase) {
    let text = tc.draw(file_header());
    let header = VcfHeader::parse(&text).unwrap_or_else(|e| panic!("{e}\n{text}"));
    let vcf = vcf_file(&text);
    assert_eq!(
        header.to_vcf_text(),
        bcftools(&["view", "-h", "--no-version"], vcf.path()),
        "{text}"
    );
    assert_eq!(dictionaries(&header), bcftools_dictionaries(&text), "{text}");
    let bcf_text = header.to_bcf_text();
    assert_eq!(bcftools_dictionaries(&bcf_text), dictionaries(&header), "{bcf_text}");
}

// r[verify vcf_header.to_builder+2]
// r[verify vcf_header.remove+2]
// r[verify vcf_header.line_order]
// r[verify vcf_header.parse.preserve+2]
/// A header derived with removals and additions is what `bcftools annotate
/// -x <removed> -h <added>` makes of the same file: the same lines in the
/// same order (the added ones last, in order; a repeated other line
/// dropped), and the dictionary htslib gives it (a removed id registered
/// again keeps its index). Its BCF text gives htslib that dictionary.
#[hegel::test(test_cases = pinned::cases(40))]
fn derived_headers_are_what_bcftools_annotate_makes(tc: TestCase) {
    let text = tc.draw(file_header());
    let header = VcfHeader::parse(&text).unwrap();
    let mut removed = Vec::new();
    let mut builder = header.to_builder();
    for id in header.filters().keys().filter(|id| *id != "PASS") {
        if tc.draw(gs::booleans()) {
            assert!(builder.remove_filter(id));
            removed.push(format!("FILTER/{id}"));
        }
    }
    let mut readd = Vec::new();
    let infos = header.infos().keys().map(|id| ("INFO", id.clone()));
    for (key, id) in infos.chain(header.formats().keys().map(|id| ("FORMAT", id.clone()))) {
        if tc.draw(gs::booleans()) {
            assert!(if key == "INFO" {
                builder.remove_info(&id)
            } else {
                builder.remove_format(&id)
            });
            removed.push(format!("{key}/{id}"));
            if tc.draw(gs::booleans()) {
                readd.push((key, id));
            }
        }
    }

    let mut added: Vec<String> = Vec::new();
    for i in 0..tc.draw(gs::integers::<usize>().max_value(2)) {
        let line = match header.other_lines().first().filter(|_| tc.draw(gs::booleans())) {
            Some(again) => again.to_string(),
            None => format!("tool{i}=v{i}"),
        };
        builder.add_other_line(line.as_str());
        added.push(format!("##{line}"));
    }
    let mut builder = builder.filters();
    if tc.draw(gs::booleans()) {
        builder.add_filter("zzF", FilterDef::new("new filter")).unwrap();
        added.push("##FILTER=<ID=zzF,Description=\"new filter\">".to_owned());
    }
    let mut builder = builder.infos();
    let new_infos = tc.draw(gs::integers::<usize>().max_value(1));
    let infos = (0..new_infos).map(|i| SmolStr::from(format!("zzI{i}")));
    let again = readd.iter().filter(|(key, _)| *key == "INFO").map(|(_, id)| id.clone());
    for id in infos.chain(again) {
        builder
            .add_info(id.clone(), InfoDef::new(Number::Count(1), ValueType::Integer, "n"))
            .unwrap();
        added.push(format!("##INFO=<ID={id},Number=1,Type=Integer,Description=\"n\">"));
    }
    let mut builder = builder.formats();
    let again = readd.iter().filter(|(key, _)| *key == "FORMAT").map(|(_, id)| id.clone());
    for id in std::iter::once(SmolStr::from("zzG")).filter(|_| tc.draw(gs::booleans())).chain(again)
    {
        let def = FormatDef::new(Number::Unknown, ValueType::String, "n");
        builder.add_format(id.clone(), def).unwrap();
        added.push(format!("##FORMAT=<ID={id},Number=.,Type=String,Description=\"n\">"));
    }
    let derived = builder.build().unwrap();

    let vcf = vcf_file(&text);
    let hdr = vcf_file(&added.iter().map(|l| l.clone() + "\n").collect::<String>());
    let bcf = tempfile::Builder::new().suffix(".bcf").tempfile().unwrap();
    let mut args = vec!["annotate", "--no-version"];
    let removed = removed.join(",");
    if !removed.is_empty() {
        args.extend(["-x", &removed]);
    }
    let hdr_path = hdr.path().to_str().unwrap();
    if !added.is_empty() {
        args.extend(["-h", hdr_path]);
    }
    if args.len() == 2 {
        args = vec!["view", "--no-version"];
    }
    let theirs = bcftools(&[&args[..], &["-O", "v"]].concat(), vcf.path());
    let theirs: String =
        theirs.lines().take_while(|l| l.starts_with('#')).map(|l| l.to_owned() + "\n").collect();
    assert_eq!(derived.to_vcf_text(), theirs, "{text}\nremoved {removed}\nadded {added:?}");

    let bcf_path = bcf.path().to_str().unwrap();
    bcftools(&[&args[..], &["-O", "b", "-o", bcf_path]].concat(), vcf.path());
    assert_eq!(dictionaries(&derived), l_text_dictionaries(bcf.path()), "{text}\n{removed}");
    let bcf_text = derived.to_bcf_text();
    assert_eq!(bcftools_dictionaries(&bcf_text), dictionaries(&derived), "{bcf_text}");
}
