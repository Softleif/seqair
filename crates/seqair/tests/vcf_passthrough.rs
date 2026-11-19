//! Passthrough writing (`docs/spec/5-vcf-passthrough.md`): records read by the
//! reader, copied through a `FieldMap` into every output format, and edited.
//! The oracle is bcftools: `bcftools view` of a copy prints what it prints for
//! the input, and `bcftools annotate -x` removes what a removal removes.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::format_push_string,
    clippy::format_collect,
    reason = "test code building VCF text"
)]

#[path = "../src/pinned.rs"]
mod pinned;

use hegel::prelude::*;
use seqair::vcf::record_encoder::{FormatFieldDef, InfoFieldDef};
use seqair::vcf::{
    FieldMap, FieldMapError, FormatGt, FullRecord, Genotype, KeyKind, Number, OutputFormat, Scalar,
    Str, VcfError, VcfFile, VcfHeader, Writer,
};
use std::path::{Path, PathBuf};
use std::process::Command;

const FORMATS: [OutputFormat; 3] = [OutputFormat::Vcf, OutputFormat::VcfGz, OutputFormat::Bcf];

fn bcftools() -> bool {
    Command::new("bcftools").arg("--version").output().is_ok()
}

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/vcf").join(name)
}

fn suffix(format: OutputFormat) -> &'static str {
    match format {
        OutputFormat::Vcf => ".vcf",
        OutputFormat::VcfGz => ".vcf.gz",
        OutputFormat::Bcf => ".bcf",
    }
}

/// `bcftools <args> <path>`'s standard output; panics on failure.
fn run(args: &[&str], path: &Path) -> String {
    let out = Command::new("bcftools").args(args).arg(path).output().unwrap();
    assert!(out.status.success(), "bcftools {args:?}: {}", String::from_utf8_lossy(&out.stderr));
    String::from_utf8(out.stdout).unwrap()
}

/// The records `bcftools view` prints for `path`.
fn view(path: &Path) -> String {
    run(&["view", "-H"], path)
}

fn persist(bytes: &[u8], format: OutputFormat) -> tempfile::NamedTempFile {
    let file = tempfile::Builder::new().suffix(suffix(format)).tempfile().unwrap();
    std::fs::write(file.path(), bytes).unwrap();
    file
}

/// Copy every record of `file` into `format` under `map`, applying `edit`.
fn copy_with(
    file: &VcfFile,
    map: &FieldMap,
    format: OutputFormat,
    mut edit: impl FnMut(&FullRecord, &mut seqair::vcf::RecordEdit<'_>),
) -> Result<Vec<u8>, VcfError> {
    let mut out = Vec::new();
    let mut writer = Writer::new(&mut out, format).write_header(map.output())?;
    let mut reader = file.reader::<FullRecord>().open().map_err(VcfError::from)?;
    while let Some(record) = reader.next_record().map_err(VcfError::from)? {
        let mut copy = writer.copy_record_with(map, record)?;
        edit(record, &mut copy);
        copy.emit()?;
    }
    writer.finish()?;
    Ok(out)
}

/// Copy every record of `file` into `format` under its own header, through
/// the writer's implicit map.
fn copy(file: &VcfFile, format: OutputFormat) -> Vec<u8> {
    let mut out = Vec::new();
    let mut writer = Writer::new(&mut out, format).write_header(file.header()).unwrap();
    let mut reader = file.reader::<FullRecord>().open().unwrap();
    while let Some(record) = reader.next_record().unwrap() {
        writer.copy_record(record).unwrap().emit().unwrap();
    }
    writer.finish().unwrap();
    out
}

// r[verify vcf_passthrough.copy]
// r[verify vcf_passthrough.implicit_map]
// r[verify vcf_passthrough.copy.equivalence]
/// The test files, copied unchanged into every format, print in bcftools as
/// the originals do.
#[test]
fn copies_of_the_test_files_print_as_the_originals() {
    if !bcftools() {
        return;
    }
    for name in ["rastair.bcf", "rastair.vcf.gz", "truth.vcf.gz", "panel.bcf"] {
        let path = data(name);
        let file = VcfFile::open(&path).unwrap();
        let expected = view(&path);
        for format in FORMATS {
            let out = persist(&copy(&file, format), format);
            assert_eq!(view(out.path()), expected, "{name} copied to {format:?}");
        }
    }
}

// ── Generated VCF text ──────────────────────────────────────────────────

/// `(ID, Number, Type)` of the INFO fields a generated header declares.
const INFO: [(&str, &str, &str); 9] = [
    ("DP", "1", "Integer"),
    ("AF", "A", "Float"),
    ("AD", "R", "Integer"),
    ("IA", ".", "Integer"),
    ("FA", ".", "Float"),
    ("FL", "0", "Flag"),
    ("S1", "1", "String"),
    ("SS", ".", "String"),
    ("END", "1", "Integer"),
];

/// FORMAT fields; DP and AD share their ids with INFO's.
const FORMAT: [(&str, &str, &str); 6] = [
    ("GT", "1", "String"),
    ("DP", "1", "Integer"),
    ("AD", "R", "Integer"),
    ("FO", ".", "Float"),
    ("XS", ".", "String"),
    ("FT", "1", "String"),
];

const FLOATS: [&str; 7] = ["0.5", "-1", "1e-05", "3.14159", "100", "0", "2.5e+08"];

fn int(tc: &TestCase) -> String {
    let v = if tc.draw(gs::booleans()) {
        tc.draw(gs::integers::<i32>().min_value(-1000).max_value(1000))
    } else {
        // The int32 reserved band (MISSING, END_OF_VECTOR, …) is not data.
        tc.draw(gs::integers::<i32>().min_value(i32::MIN + 8))
    };
    v.to_string()
}

fn float(tc: &TestCase) -> String {
    tc.draw(gs::sampled_from(FLOATS.to_vec())).to_owned()
}

/// A string value: one word, or for a field that may hold a list, words
/// joined by `,` (which a copy must not percent-encode).
fn word(tc: &TestCase, list: bool) -> String {
    let n_words = if list { tc.draw(gs::integers::<usize>().min_value(1).max_value(3)) } else { 1 };
    let words: Vec<String> = (0..n_words)
        .map(|_| {
            let n = tc.draw(gs::integers::<usize>().min_value(1).max_value(4));
            (0..n).map(|_| tc.draw(gs::sampled_from(vec!['a', 'Z', '0', '_', '%']))).collect()
        })
        .collect();
    words.join(",")
}

/// A value of a field declared `number`, `typ` in a record of `n_allele`
/// alleles; `None` where the record cannot hold the field (`A` without ALT).
fn value(tc: &TestCase, number: &str, typ: &str, n_allele: usize) -> Option<String> {
    let count = match number {
        "A" => n_allele - 1,
        "R" => n_allele,
        "." => tc.draw(gs::integers::<usize>().min_value(1).max_value(3)),
        n => n.parse().unwrap(),
    };
    if count == 0 {
        return None;
    }
    if tc.draw(gs::integers::<u8>().max_value(5)) == 0 {
        return Some(".".to_owned());
    }
    let element = |tc: &TestCase| match typ {
        "Integer" | "Float" if tc.draw(gs::integers::<u8>().max_value(5)) == 0 => ".".to_owned(),
        "Integer" => int(tc),
        "Float" => float(tc),
        _ => word(tc, number == "."),
    };
    if typ == "String" {
        return Some(element(tc));
    }
    Some((0..count).map(|_| element(tc)).collect::<Vec<_>>().join(","))
}

fn genotype(tc: &TestCase, n_allele: usize) -> String {
    if tc.draw(gs::integers::<u8>().max_value(6)) == 0 {
        return ".".to_owned();
    }
    let ploidy = tc.draw(gs::integers::<usize>().min_value(1).max_value(3));
    let mut gt = String::new();
    for i in 0..ploidy {
        if i > 0 {
            gt.push(if tc.draw(gs::booleans()) { '|' } else { '/' });
        }
        match tc.draw(gs::optional(gs::integers::<usize>().max_value(n_allele - 1))) {
            Some(a) => gt.push_str(&a.to_string()),
            None => gt.push('.'),
        }
    }
    gt
}

fn shuffled<T: Clone>(tc: &TestCase, items: &[T]) -> Vec<T> {
    let mut keyed: Vec<(u32, T)> =
        items.iter().map(|i| (tc.draw(gs::integers::<u32>()), i.clone())).collect();
    keyed.sort_by_key(|(k, _)| *k);
    keyed.into_iter().map(|(_, i)| i).collect()
}

fn bases(tc: &TestCase) -> String {
    let n = tc.draw(gs::integers::<usize>().min_value(1).max_value(3));
    (0..n).map(|_| tc.draw(gs::sampled_from(vec!['A', 'C', 'G', 'T']))).collect()
}

/// A VCF file: a header declaring [`INFO`], [`FORMAT`], two FILTERs and two
/// contigs, and sorted records drawing on all of it.
fn vcf_text(tc: &TestCase) -> String {
    vcf_text_with(tc, false)
}

/// [`vcf_text`]; `gt`: every record with samples has GT.
fn vcf_text_with(tc: &TestCase, gt: bool) -> String {
    let version = tc.draw(gs::sampled_from(vec!["4.2", "4.3", "4.4", "4.5"]));
    let n_samples = tc.draw(gs::integers::<usize>().max_value(3));
    let mut text = format!("##fileformat=VCFv{version}\n");
    text.push_str("##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
    text.push_str("##FILTER=<ID=q10,Description=\"q\">\n##FILTER=<ID=lowDP,Description=\"d\">\n");
    for (id, n, t) in INFO {
        text.push_str(&format!("##INFO=<ID={id},Number={n},Type={t},Description=\"i\">\n"));
    }
    for (id, n, t) in FORMAT {
        text.push_str(&format!("##FORMAT=<ID={id},Number={n},Type={t},Description=\"f\">\n"));
    }
    text.push_str("##contig=<ID=chr1,length=100000>\n##contig=<ID=chr2,length=100000>\n");
    text.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
    if n_samples > 0 {
        text.push_str("\tFORMAT");
        for i in 0..n_samples {
            text.push_str(&format!("\tS{i}"));
        }
    }
    text.push('\n');
    let n_records = tc.draw(gs::integers::<usize>().min_value(1).max_value(6));
    let mut sites: Vec<(usize, u32)> = (0..n_records)
        .map(|_| {
            let contig = tc.draw(gs::integers::<usize>().max_value(1));
            (contig, tc.draw(gs::integers::<u32>().min_value(1).max_value(90_000)))
        })
        .collect();
    sites.sort_unstable();
    for (contig, pos) in sites {
        let reference = bases(tc);
        let n_alt = tc.draw(gs::integers::<usize>().max_value(2));
        let alts: Vec<String> = (0..n_alt)
            .map(|_| match tc.draw(gs::integers::<u8>().max_value(4)) {
                0 => "<DEL>".to_owned(),
                1 => "*".to_owned(),
                _ => bases(tc),
            })
            .collect();
        let n_allele = 1 + n_alt;
        let ids = match tc.draw(gs::integers::<u8>().max_value(2)) {
            0 => ".".to_owned(),
            1 => format!("rs{}", tc.draw(gs::integers::<u16>())),
            _ => {
                format!("rs{};id{}", tc.draw(gs::integers::<u16>()), tc.draw(gs::integers::<u16>()))
            }
        };
        let qual = if tc.draw(gs::booleans()) { float(tc) } else { ".".to_owned() };
        let filter = tc.draw(gs::sampled_from(vec!["PASS", ".", "q10", "q10;lowDP"]));
        let mut info = Vec::new();
        for (id, n, t) in shuffled(tc, &INFO) {
            if !tc.draw(gs::booleans()) {
                continue;
            }
            match (id, t) {
                (_, "Flag") => info.push(id.to_owned()),
                ("END", _) => {
                    let end = pos + tc.draw(gs::integers::<u32>().max_value(500));
                    info.push(format!("END={end}"));
                }
                _ => {
                    if let Some(v) = value(tc, n, t, n_allele) {
                        info.push(format!("{id}={v}"));
                    }
                }
            }
        }
        let info = if info.is_empty() { ".".to_owned() } else { info.join(";") };
        let alt = if alts.is_empty() { ".".to_owned() } else { alts.join(",") };
        let chrom = ["chr1", "chr2"][contig];
        let mut line =
            format!("{chrom}\t{pos}\t{ids}\t{reference}\t{alt}\t{qual}\t{filter}\t{info}");
        if n_samples > 0 {
            let mut keys: Vec<(&str, &str, &str)> = Vec::new();
            if gt || tc.draw(gs::booleans()) {
                keys.push(FORMAT[0]);
            }
            let rest: Vec<_> = FORMAT[1..]
                .iter()
                .copied()
                .filter(|(id, ..)| id != &"AD" || n_allele > 0)
                .collect();
            keys.extend(shuffled(tc, &rest).into_iter().filter(|_| tc.draw(gs::booleans())));
            if keys.is_empty() {
                keys.push(FORMAT[1]);
            }
            let names: Vec<&str> = keys.iter().map(|(id, ..)| *id).collect();
            line.push('\t');
            line.push_str(&names.join(":"));
            for _ in 0..n_samples {
                let mut values: Vec<String> = keys
                    .iter()
                    .map(|(id, n, t)| {
                        if *id == "GT" {
                            genotype(tc, n_allele)
                        } else {
                            value(tc, n, t, n_allele).unwrap_or_else(|| ".".to_owned())
                        }
                    })
                    .collect();
                // A sample may drop trailing fields.
                let keep = tc.draw(gs::integers::<usize>().min_value(1).max_value(values.len()));
                values.truncate(keep);
                line.push('\t');
                line.push_str(&values.join(":"));
            }
        }
        text.push_str(&line);
        text.push('\n');
    }
    text
}

/// What noodles reads from `path`: every record, fully decoded, as text.
/// htslib stores a FORMAT string no sample of the record holds as a NUL,
/// which noodles reads as `""` (or `[""]` for `Number=.`); seqair writes the missing `.`, read as `None`.
fn noodles_records(path: &Path) -> String {
    format!("{:?}", noodles_bufs(path))
        .replace("Some(Array(String([Some(\"\")])))", "None")
        .replace("Some(String(\"\"))", "None")
}

fn noodles_bufs(path: &Path) -> Vec<noodles::vcf::variant::RecordBuf> {
    let name = path.to_str().unwrap();
    if name.ends_with(".bcf") {
        let mut reader =
            noodles::bcf::io::reader::Builder::default().build_from_path(path).unwrap();
        let header = reader.read_header().unwrap();
        reader.record_bufs(&header).collect::<Result<_, _>>().unwrap()
    } else {
        let mut reader =
            noodles::vcf::io::reader::Builder::default().build_from_path(path).unwrap();
        let header = reader.read_header().unwrap();
        reader.record_bufs(&header).collect::<Result<_, _>>().unwrap()
    }
}

// r[verify vcf_passthrough.copy]
// r[verify vcf_passthrough.copy.encoder]
// r[verify vcf_passthrough.copy.values]
// r[verify vcf_passthrough.copy.equivalence]
/// Generated VCF text, and its BCF from bcftools, copied into every format:
/// bcftools prints every copy as it prints the input, and noodles reads each
/// copy as it reads bcftools' own file in that format.
#[hegel::test(test_cases = pinned::cases(30))]
fn every_copy_prints_in_bcftools_as_its_input(tc: TestCase) {
    if !bcftools() {
        return;
    }
    let text = vcf_text(&tc);
    let dir = tempfile::tempdir().unwrap();
    let vcf = dir.path().join("in.vcf");
    std::fs::write(&vcf, &text).unwrap();
    let expected = view(&vcf);
    // noodles percent-decodes VCF text strings and refuses a `%` that does
    // not start an escape, which htslib (and so a copy) keeps as written.
    let noodles = !text.contains('%');
    let mut reference = Vec::new();
    for format in FORMATS {
        let path = dir.path().join(format!("ref{}", suffix(format)));
        let flag = match format {
            OutputFormat::Vcf => "-Ov",
            OutputFormat::VcfGz => "-Oz",
            OutputFormat::Bcf => "-Ob",
        };
        run(&["view", "--no-version", flag, "-o", path.to_str().unwrap()], &vcf);
        reference.push(noodles.then(|| noodles_records(&path)));
    }
    let bcf = dir.path().join("refBcf.bcf");
    std::fs::rename(dir.path().join("ref.bcf"), &bcf).unwrap();
    for input in [&vcf, &bcf] {
        let file = VcfFile::open(input).unwrap();
        for (format, reference) in FORMATS.into_iter().zip(&reference) {
            let out = persist(&copy(&file, format), format);
            assert_eq!(view(out.path()), expected, "{input:?} to {format:?}\n{text}");
            if let Some(reference) = reference {
                let read = noodles_records(out.path());
                assert_eq!(&read, reference, "{input:?} to {format:?}\n{text}");
            }
        }
    }
}

// ── Field maps ──────────────────────────────────────────────────────────

const SMALL: &str = "##fileformat=VCFv4.2
##FILTER=<ID=q10,Description=\"q\">
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">
##INFO=<ID=AF,Number=A,Type=Float,Description=\"a\">
##INFO=<ID=S1,Number=1,Type=String,Description=\"s\">
##INFO=<ID=END,Number=1,Type=Integer,Description=\"e\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"d\">
##FORMAT=<ID=PS,Number=1,Type=String,Description=\"old phase set\">
##contig=<ID=chr1,length=1000>
##contig=<ID=chr2,length=1000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS0\tS1
chr1\t10\trs1\tA\tG,T\t30\tq10\tDP=5;AF=0.5,0.25;S1=a%b,c\tGT:DP:PS\t0/1:3:x\t1|2:.:y
chr2\t20\t.\tAC\t<DEL>\t.\tPASS\tEND=120;DP=7\tGT:DP\t0\t./.:4
";

fn small() -> VcfFile {
    VcfFile::from_bytes(SMALL.as_bytes().to_vec(), None).unwrap()
}

fn text(bytes: &[u8], format: OutputFormat) -> String {
    view(persist(bytes, format).path())
}

// r[verify vcf_passthrough.field_map.compatible]
#[test]
fn a_field_redefined_with_another_type_or_number_fails_the_map() {
    let file = small();
    // Not derived from the input's header: only the names relate the two.
    let unrelated = |h: VcfHeader| VcfHeader::parse(&h.to_vcf_text()).unwrap();
    let mut b = file.header().to_builder();
    assert!(b.remove_info("DP"));
    let mut b = b.infos();
    b.register_info(&InfoFieldDef::<Scalar<f32>>::new("DP", "depth as a float")).unwrap();
    let err = FieldMap::new(&file, &unrelated(b.build().unwrap())).unwrap_err();
    assert!(
        matches!(&err, FieldMapError::Incompatible { kind: KeyKind::Info, name, input_type, output_type, .. }
            if name == "DP" && *input_type == seqair::vcf::ValueType::Integer
                && *output_type == seqair::vcf::ValueType::Float),
        "{err:?}"
    );

    let mut b = file.header().to_builder();
    assert!(b.remove_format("PS"));
    let mut b = b.formats();
    b.register_format(&FormatFieldDef::<Str>::new("PS", Number::Unknown, "any count")).unwrap();
    let err = FieldMap::new(&file, &unrelated(b.build().unwrap())).unwrap_err();
    assert!(
        matches!(&err, FieldMapError::Incompatible { kind: KeyKind::Format, name, input_number, output_number, .. }
            if name == "PS" && *input_number == Number::Count(1) && *output_number == Number::Unknown),
        "{err:?}"
    );
}

// r[verify vcf_passthrough.field_map]
// r[verify vcf_passthrough.field_map.contigs]
// r[verify vcf_passthrough.field_map.samples]
#[test]
fn contigs_and_samples_are_checked_once() {
    let file = small();
    let mut b = VcfHeader::builder();
    b.register_contig("chr1", seqair::vcf::ContigDef::new(Some(1000))).unwrap();
    let mut b = b.samples();
    b.add_sample("S0").unwrap();
    b.add_sample("S1").unwrap();
    let err = FieldMap::new(&file, &b.build().unwrap()).unwrap_err();
    assert_eq!(err, FieldMapError::MissingContig { name: "chr2".into() });

    let mut b = file.header().to_builder().samples();
    b.set_samples(["S0", "S2"]).unwrap();
    let err = FieldMap::new(&file, &b.build().unwrap()).unwrap_err();
    assert_eq!(
        err,
        FieldMapError::SampleMismatch {
            index: 1,
            input: Some("S1".into()),
            output: Some("S2".into())
        }
    );
}

// r[verify vcf_passthrough.copy]
#[test]
fn a_record_or_writer_of_another_map_is_refused() {
    let file = small();
    let other = small();
    let map = FieldMap::new(&file, file.header()).unwrap();
    let mut out = Vec::new();
    let mut writer = Writer::new(&mut out, OutputFormat::Vcf).write_header(other.header()).unwrap();
    let mut reader = file.reader::<FullRecord>().open().unwrap();
    let record = reader.next_record().unwrap().unwrap();
    let err = writer.copy_record_with(&map, record).err().unwrap();
    assert!(matches!(err, VcfError::ForeignFieldMap), "{err:?}");

    let map = FieldMap::new(&other, other.header()).unwrap();
    let err = writer.copy_record_with(&map, record).err().unwrap();
    assert!(matches!(err, VcfError::ForeignRecord), "{err:?}");
    let err = writer.copy_record(record).err().unwrap();
    assert!(matches!(err, VcfError::ForeignRecord), "{err:?}");
}

// r[verify vcf_passthrough.copy]
#[test]
fn a_record_at_the_telomere_is_refused() {
    let text = SMALL.replace("chr1\t10\t", "chr1\t0\t");
    let file = VcfFile::from_bytes(text.into_bytes(), None).unwrap();
    let map = FieldMap::new(&file, file.header()).unwrap();
    let err = copy_with(&file, &map, OutputFormat::Vcf, |_, _| {}).unwrap_err();
    assert!(matches!(&err, VcfError::TelomereRecord { contig } if contig == "chr1"), "{err:?}");
}

// r[verify vcf_header.lineage]
// r[verify record_encoder.write_header+2]
// r[verify vcf_passthrough.field_map.removed]
// r[verify vcf_passthrough.implicit_map]
/// A field removed from a derived header and registered again is a new
/// field: the input's values of that name are not copied into it, even
/// where the definitions agree; edits write it.
#[test]
fn a_field_removed_and_registered_again_is_a_new_field() {
    if !bcftools() {
        return;
    }
    let file = small();
    let mut b = file.header().to_builder();
    assert!(b.remove_format("PS"));
    assert!(b.remove_info("S1"));
    let mut b = b.infos();
    b.register_info(&InfoFieldDef::<Str>::new("S1", Number::Count(1), "same as before")).unwrap();
    let mut b = b.formats();
    let ps = b.register_format(&FormatFieldDef::<Scalar<i32>>::new("PS", "Phase set")).unwrap();
    let output = b.build().unwrap();
    assert_eq!(output.parent(), Some(file.header().id()));
    for format in FORMATS {
        let mut out = Vec::new();
        let mut writer = Writer::new(&mut out, format).write_header(&output).unwrap();
        let mut reader = file.reader::<FullRecord>().open().unwrap();
        let mut first = true;
        while let Some(record) = reader.next_record().unwrap() {
            let mut copy = writer.copy_record(record).unwrap();
            if first {
                ps.encode_opt(&mut copy, &[Some(7), None]).unwrap();
            }
            first = false;
            copy.emit().unwrap();
        }
        writer.finish().unwrap();
        let lines = text_of(&out, format);
        let lines: Vec<&str> = lines.lines().collect();
        assert_eq!(
            lines[0],
            "chr1\t10\trs1\tA\tG,T\t30\tq10\tDP=5;AF=0.5,0.25\tGT:DP:PS\t0/1:3:7\t1|2:.:.",
            "{format:?}"
        );
        assert_eq!(lines[1], "chr2\t20\t.\tAC\t<DEL>\t.\tPASS\tEND=120;DP=7\tGT:DP\t0:.\t./.:4");
    }
}

// r[verify vcf_passthrough.copy.values]
/// A copied string keeps `%` and `,` as stored; only what would break the
/// line is percent-encoded, here a `;` that only BCF can hold in INFO.
#[test]
fn copied_strings_are_written_as_stored() {
    let file = small();
    let out = copy(&file, OutputFormat::Vcf);
    assert!(text(&out, OutputFormat::Vcf).contains("S1=a%b,c"));

    let header = file.header();
    let s1: seqair::vcf::InfoString = header.info_key("S1").unwrap();
    let chr1 = header.contig("chr1").unwrap();
    let mut bcf = Vec::new();
    let mut w = Writer::new(&mut bcf, OutputFormat::Bcf).write_header(header).unwrap();
    let alleles = seqair::vcf::Alleles::snv(seqair_types::Base::A, seqair_types::Base::G).unwrap();
    let mut enc = w
        .begin_record(&chr1, seqair_types::Pos1::new(5).unwrap(), &alleles, None)
        .unwrap()
        .no_filter();
    s1.encode(&mut enc, "x;y:z").unwrap();
    enc.begin_samples().emit().unwrap();
    w.finish().unwrap();
    let bcf = VcfFile::from_bytes(bcf, None).unwrap();
    let line = text(&copy(&bcf, OutputFormat::Vcf), OutputFormat::Vcf);
    assert!(line.contains("S1=x%3By:z"), "{line}");
}

// ── Removal and undeclared keys ─────────────────────────────────────────

// r[verify vcf_passthrough.field_map.removed]
// r[verify vcf_passthrough.edit.remove]
/// Fields removed from the output header, or from each copied record, leave
/// the copy as `bcftools annotate -x` leaves the input; a record whose
/// failed filters are all removed passes.
#[hegel::test(test_cases = pinned::cases(20))]
fn removals_match_bcftools_annotate(tc: TestCase) {
    if !bcftools() {
        return;
    }
    let text = vcf_text(&tc);
    let dir = tempfile::tempdir().unwrap();
    let vcf = dir.path().join("in.vcf");
    std::fs::write(&vcf, &text).unwrap();
    let infos: Vec<&str> =
        INFO.iter().map(|(id, ..)| *id).filter(|_| tc.draw(gs::booleans())).collect();
    let formats: Vec<&str> =
        FORMAT.iter().map(|(id, ..)| *id).filter(|_| tc.draw(gs::booleans())).collect();
    let drop_q10 = tc.draw(gs::booleans());
    let in_header = tc.draw(gs::booleans());
    let mut removed: Vec<String> = infos.iter().map(|id| format!("INFO/{id}")).collect();
    removed.extend(formats.iter().map(|id| format!("FORMAT/{id}")));
    if drop_q10 {
        removed.push("FILTER/q10".to_owned());
    }
    let expected = if removed.is_empty() {
        view(&vcf)
    } else {
        run(&["annotate", "-x", &removed.join(","), "-O", "v"], &vcf)
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| format!("{l}\n"))
            .collect()
    };
    let file = VcfFile::open(&vcf).unwrap();
    let mut b = file.header().to_builder();
    if drop_q10 {
        b.remove_filter("q10");
    }
    if in_header {
        for id in &infos {
            b.remove_info(id);
        }
        for id in &formats {
            b.remove_format(id);
        }
    }
    let output = b.build().unwrap();
    let map = FieldMap::new(&file, &output).unwrap();
    let info_keys: Vec<seqair::vcf::InfoAny> =
        infos.iter().filter(|_| !in_header).map(|id| output.info_key(id).unwrap()).collect();
    let format_keys: Vec<seqair::vcf::FormatAny> =
        formats.iter().filter(|_| !in_header).map(|id| output.format_key(id).unwrap()).collect();
    for format in FORMATS {
        let out = copy_with(&file, &map, format, |_, edit| {
            for key in &info_keys {
                edit.remove_info(key).unwrap();
            }
            for key in &format_keys {
                edit.remove_format(key).unwrap();
            }
        })
        .unwrap();
        assert_eq!(text_of(&out, format), expected, "{format:?}, removed {removed:?}\n{text}");
    }
}

fn text_of(bytes: &[u8], format: OutputFormat) -> String {
    text(bytes, format)
}

// r[verify vcf_writer.format_serialization+2]
/// A record without FORMAT fields in a file with samples keeps every column.
#[test]
fn a_record_without_format_fields_keeps_its_sample_columns() {
    if !bcftools() {
        return;
    }
    let file = small();
    let gt: FormatGt = file.header().format_key("GT").unwrap();
    let dp: seqair::vcf::FormatInt = file.header().format_key("DP").unwrap();
    let ps: seqair::vcf::FormatString = file.header().format_key("PS").unwrap();
    let map = FieldMap::new(&file, file.header()).unwrap();
    for format in FORMATS {
        let out = copy_with(&file, &map, format, |_, edit| {
            edit.remove_format(&gt).unwrap();
            edit.remove_format(&dp).unwrap();
            edit.remove_format(&ps).unwrap();
        })
        .unwrap();
        let lines = text_of(&out, format);
        assert!(lines.lines().all(|l| l.ends_with("\t.\t.\t.")), "{format:?}: {lines}");
    }
}

const UNDECLARED: &str = "##fileformat=VCFv4.2
##FILTER=<ID=q10,Description=\"q\">
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">
##contig=<ID=chr1,length=1000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS0
chr1\t10\t.\tA\tG\t.\tq10;foo\tDP=5;TRID=x1\tGT:XX\t0/1:abc
chr1\t20\t.\tA\tG\t.\tPASS\tDP=6\tGT\t1/1
";

// r[verify vcf_passthrough.undeclared]
/// Keys and filters no header declares fail the copy, or are left out, or,
/// declared in the output header, kept as bcftools keeps them.
#[test]
fn undeclared_keys_fail_are_dropped_or_declared() {
    if !bcftools() {
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let vcf = dir.path().join("in.vcf");
    std::fs::write(&vcf, UNDECLARED).unwrap();
    let file = VcfFile::open(&vcf).unwrap();

    let map = FieldMap::new(&file, file.header()).unwrap();
    let err = copy_with(&file, &map, OutputFormat::Vcf, |_, _| {}).unwrap_err();
    assert!(
        matches!(&err, VcfError::UndeclaredField { kind: KeyKind::Filter, name } if name == "foo"),
        "{err:?}"
    );

    // bcftools refuses to annotate keys its header lacks; the lines it would write:
    let dropped =
        "chr1\t10\t.\tA\tG\t.\tq10\tDP=5\tGT\t0/1\nchr1\t20\t.\tA\tG\t.\tPASS\tDP=6\tGT\t1/1\n";
    let map = FieldMap::new(&file, file.header()).unwrap().drop_undeclared();
    for format in FORMATS {
        let out = copy_with(&file, &map, format, |_, _| {}).unwrap();
        assert_eq!(text_of(&out, format), dropped, "{format:?}");
    }

    let mut b = file.header().to_builder().filters();
    b.register_filter(&seqair::vcf::FilterFieldDef::new("foo", "f")).unwrap();
    let mut b = b.infos();
    b.register_info(&InfoFieldDef::<Str>::new("TRID", Number::Count(1), "t")).unwrap();
    let mut b = b.formats();
    b.register_format(&FormatFieldDef::<Str>::new("XX", Number::Count(1), "x")).unwrap();
    let output = b.build().unwrap();
    let map = FieldMap::new(&file, &output).unwrap();
    for format in FORMATS {
        let out = copy_with(&file, &map, format, |_, _| {}).unwrap();
        assert_eq!(text_of(&out, format), view(&vcf), "{format:?}");
    }

    let mut b = file.header().to_builder().infos();
    b.register_info(&InfoFieldDef::<Scalar<i32>>::new("TRID", "not a string")).unwrap();
    let map = FieldMap::new(&file, &b.build().unwrap()).unwrap().drop_undeclared();
    let err = copy_with(&file, &map, OutputFormat::Vcf, |_, _| {}).unwrap_err();
    assert!(
        matches!(&err, VcfError::FieldMap(FieldMapError::Incompatible { name, .. }) if name == "TRID"),
        "{err:?}"
    );
}

// ── Edits ───────────────────────────────────────────────────────────────

// r[verify vcf_passthrough.edit]
// r[verify vcf_passthrough.edit.emit]
/// `HiPhase`'s rewrite: every sample's GT replaced (phased, its ploidy kept),
/// PS added as an optional integer, PF on some records, the input's String
/// PS removed from the header; everything else is the input's.
#[hegel::test(test_cases = pinned::cases(20))]
fn rewritten_genotypes_and_phase_sets_read_back_and_the_rest_is_copied(tc: TestCase) {
    if !bcftools() {
        return;
    }
    let text = vcf_text_with(&tc, true);
    let dir = tempfile::tempdir().unwrap();
    let vcf = dir.path().join("in.vcf");
    std::fs::write(&vcf, &text).unwrap();
    let file = VcfFile::open(&vcf).unwrap();
    if file.header().samples().is_empty() {
        return;
    }
    let mut b = file.header().to_builder();
    b.add_other_line("phaser=test");
    let mut b = b.formats();
    let ps = b.register_format(&FormatFieldDef::<Scalar<i32>>::new("PS", "Phase set")).unwrap();
    let pf: seqair::vcf::FormatString =
        b.register_format(&FormatFieldDef::<Str>::new("PF", Number::Count(1), "Flag")).unwrap();
    let output = b.build().unwrap();
    let in_gt: FormatGt = file.header().format_key("GT").unwrap();
    let out_gt: FormatGt = output.format_key("GT").unwrap();
    let map = FieldMap::new(&file, &output).unwrap();
    let format = tc.draw(gs::sampled_from(FORMATS.to_vec()));
    let mut expected = String::new();
    let out = copy_with(&file, &map, format, |record, edit| {
        let pos = record.locus().start().pos1().unwrap().as_u32();
        let gts: Vec<Genotype> = record
            .samples()
            .column(&in_gt)
            .unwrap()
            .unwrap()
            .map(|gt| {
                let n = gt.map_or(2, |g| g.ploidy());
                let mut g = Genotype::phased_diploid(1, 0);
                g.alleles.truncate(n);
                g.phased.truncate(n.saturating_sub(1));
                while g.alleles.len() < n {
                    g.alleles.push(Some(0));
                    g.phased.push(true);
                }
                g
            })
            .collect();
        out_gt.encode(edit, &gts).unwrap();
        let sets: Vec<Option<i32>> =
            (0..gts.len()).map(|i| (i % 2 == 0).then_some(i32::try_from(pos).unwrap())).collect();
        ps.encode_opt(edit, &sets).unwrap();
        let flagged = pos % 2 == 0;
        if flagged {
            let flags: Vec<&str> = (0..gts.len()).map(|_| "TR_OVERLAP").collect();
            pf.encode(edit, &flags).unwrap();
        }
        let samples: Vec<String> = gts
            .iter()
            .zip(&sets)
            .map(|(g, s)| {
                let gt: Vec<String> = g.alleles.iter().map(|a| a.unwrap().to_string()).collect();
                let ps = s.map_or_else(|| ".".to_owned(), |s| s.to_string());
                let pf = if flagged { "TR_OVERLAP" } else { "." };
                format!("{}\t{ps}\t{pf}", gt.join("|"))
            })
            .collect();
        expected.push_str(&samples.join("\t"));
        expected.push('\n');
    })
    .unwrap();
    let out = persist(&out, format);
    let query = run(&["query", "-f", "[%GT\\t%PS\\t%PF\\t]\\n"], out.path());
    let query: String = query.lines().map(|l| format!("{}\n", l.trim_end_matches('\t'))).collect();
    assert_eq!(query, expected, "{format:?}\n{text}");
    let rest = |path: &Path| {
        run(&["annotate", "-x", "FORMAT/GT,FORMAT/PS,FORMAT/PF", "-O", "v"], path)
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| format!("{l}\n"))
            .collect::<String>()
    };
    assert_eq!(rest(out.path()), rest(&vcf), "{format:?}\n{text}");
}

// r[verify vcf_passthrough.edit]
/// An encode replaces a copied field where it stands, appends a new one,
/// and INFO may follow FORMAT.
#[test]
fn edits_replace_in_place_in_any_order() {
    if !bcftools() {
        return;
    }
    let file = small();
    let header = file.header();
    let dp: seqair::vcf::InfoInt = header.info_key("DP").unwrap();
    let fdp: seqair::vcf::FormatInt = header.format_key("DP").unwrap();
    let s1: seqair::vcf::InfoString = header.info_key("S1").unwrap();
    let map = FieldMap::new(&file, header).unwrap();
    for format in FORMATS {
        let out = copy_with(&file, &map, format, |_, edit| {
            fdp.encode(edit, &[9, 8]).unwrap();
            dp.encode(edit, 42).unwrap();
            s1.encode(edit, "new").unwrap();
        })
        .unwrap();
        let lines = text_of(&out, format);
        let lines: Vec<&str> = lines.lines().collect();
        assert_eq!(
            lines[0],
            "chr1\t10\trs1\tA\tG,T\t30\tq10\tDP=42;AF=0.5,0.25;S1=new\tGT:DP:PS\t0/1:9:x\t1|2:8:y",
            "{format:?}"
        );
        assert_eq!(
            lines[1], "chr2\t20\t.\tAC\t<DEL>\t.\tPASS\tEND=120;DP=42;S1=new\tGT:DP\t0:9\t./.:8",
            "{format:?}"
        );
    }
}

// r[verify vcf_passthrough.edit.remove]
/// Removing END shortens the record to its REF, which bcftools reports as
/// the end of a BCF record without END.
#[test]
fn removing_end_takes_it_out_of_the_length() {
    if !bcftools() {
        return;
    }
    let file = small();
    let end: seqair::vcf::InfoInt = file.header().info_key("END").unwrap();
    let gt: FormatGt = file.header().format_key("GT").unwrap();
    let map = FieldMap::new(&file, file.header()).unwrap();
    let kept = copy_with(&file, &map, OutputFormat::Bcf, |_, _| {}).unwrap();
    let mut had_end = Vec::new();
    let removed = copy_with(&file, &map, OutputFormat::Bcf, |_, edit| {
        had_end.push(edit.remove_info(&end).unwrap());
        assert!(!edit.remove_info(&end).unwrap());
        assert!(edit.remove_format(&gt).unwrap());
        assert!(!edit.remove_format(&gt).unwrap());
    })
    .unwrap();
    assert_eq!(had_end, [false, true]);
    let ends = |bytes: &[u8]| {
        run(&["query", "-f", "%POS %END\\n"], persist(bytes, OutputFormat::Bcf).path())
    };
    assert_eq!(ends(&kept), "10 10\n20 120\n");
    assert_eq!(ends(&removed), "10 10\n20 21\n");
}
