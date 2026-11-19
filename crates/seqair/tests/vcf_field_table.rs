//! A VCF text file's field table (`r[vcf_reader.field_table]`): keys the
//! header leaves out resolve through a lenient key — reserved ones with
//! their VCF 4.5 definitions, others as Strings — and read what htslib
//! reads; plain keys resolve the header's fields only, in VCF and BCF alike.
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

use hegel::prelude::*;
use rust_htslib::bcf::{self as hts, Read as _};
use seqair::vcf::{
    DynValue, FormatInt, FullRecord, InfoFlag, InfoInt, InfoInts, InfoKey, InfoString, KeyError,
    KeyKind, Locus, Number, OutputFormat, Scalar, Str, ValueType, VcfError, VcfFile, VcfPos,
    Writer,
};
use std::io::Write as _;

const HEADER: &str = "##fileformat=VCFv4.3\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
##contig=<ID=chr1,length=100000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n";

fn text(lines: &[&str]) -> String {
    let mut text = HEADER.to_owned();
    for l in lines {
        text.push_str(l);
        text.push('\n');
    }
    text
}

fn file(lines: &[&str]) -> VcfFile {
    VcfFile::from_bytes(text(lines).into_bytes(), None).unwrap()
}

// r[verify vcf_values.file_keys+3]
// r[verify vcf_values.optional_keys+2]
// r[verify vcf_values.lenient_keys]
// r[verify vcf_reader.field_table]
#[test]
fn undeclared_keys_resolve_through_lenient_keys() {
    let file = file(&[
        "chr1\t5\t.\tA\t<DEL>\t.\t.\tSVTYPE=DEL;SVLEN=10;END=20;DB;TRID=t1;AN=4\tGT:GQ:ZZ\t0/1:30:q\t1/1:.:.",
    ]);
    // The header resolves only what it declares: it is what a writer writes.
    let err = file.header().info_key::<Str>("SVTYPE").unwrap_err();
    assert_eq!(err, KeyError::UnknownField { name: "SVTYPE".into(), kind: KeyKind::Info });
    // So does the file, unless asked to be lenient.
    assert_eq!(file.info_key::<Str>("SVTYPE").unwrap_err(), err);
    assert_eq!(file.info_key_opt::<Str>("TRID").unwrap(), None);
    assert_eq!(file.format_key_opt::<Scalar<i32>>("GQ").unwrap(), None);
    let lenient = file.lenient();

    // Reserved, undeclared: typed by VCF 4.5.
    let svtype: InfoString = lenient.info_key("SVTYPE").unwrap();
    let svlen: InfoInts = lenient.info_key("SVLEN").unwrap();
    let db: InfoFlag = lenient.info_key("DB").unwrap();
    let an: InfoInt = lenient.info_key("AN").unwrap();
    let gq: FormatInt = lenient.format_key("GQ").unwrap();
    assert_eq!(svlen.id().number(), Number::AlternateBases);
    assert!(matches!(
        lenient.info_key::<Scalar<f32>>("AN"),
        Err(KeyError::TypeMismatch { declared: ValueType::Integer, .. })
    ));
    assert!(matches!(
        lenient.info_key::<Scalar<i32>>("SVLEN"),
        Err(KeyError::NumberMismatch { .. })
    ));
    // Declared keys resolve as the header declares them.
    let dp: InfoInt = lenient.info_key("DP").unwrap();
    assert_eq!(dp, file.header().info_key::<Scalar<i32>>("DP").unwrap());
    assert_eq!(dp, file.info_key::<Scalar<i32>>("DP").unwrap());

    // Other undeclared keys: Strings, for a Str (or Any) key only.
    let trid: InfoString = lenient.info_key("TRID").unwrap();
    let zz: seqair::vcf::FormatString = lenient.format_key("ZZ").unwrap();
    assert_eq!(trid.id().number(), Number::Unknown);
    assert_eq!(
        lenient.info_key::<Scalar<i32>>("TRID").unwrap_err(),
        KeyError::UnknownField { name: "TRID".into(), kind: KeyKind::Info }
    );
    assert_eq!(lenient.info_key_opt::<Scalar<i32>>("TRID").unwrap(), None);
    assert!(lenient.info_key_opt::<Str>("NOPE").unwrap().is_some(), "a String key resolves");
    assert!(matches!(
        lenient.format_key_opt::<Scalar<f32>>("GQ"),
        Err(KeyError::TypeMismatch { .. })
    ));

    let rec = file.reader::<FullRecord>().open().unwrap().read().unwrap().unwrap();
    let info = rec.info();
    assert_eq!(info.get(&svtype).unwrap(), Some("DEL"));
    assert_eq!(info.get(&svlen).unwrap().map(|v| v.collect::<Vec<_>>()), Some(vec![10]));
    assert!(info.get(&db).unwrap());
    assert_eq!(info.get(&an).unwrap(), Some(4));
    assert_eq!(info.get(&trid).unwrap(), Some("t1"));
    assert_eq!(info.get(&dp).unwrap(), None);
    // r[verify vcf_values.access+2]
    // Present whatever its type, through `Any`.
    let any = |name| -> seqair::vcf::InfoAny { lenient.info_key(name).unwrap() };
    assert!(info.get(&any("TRID")).unwrap().is_some_and(|v| !v.is_missing()));
    assert!(matches!(info.get(&any("DB")).unwrap(), Some(DynValue::Flag)));
    assert!(matches!(info.get(&any("AN")).unwrap(), Some(DynValue::Int(_))));
    assert!(info.get(&any("NOPE")).unwrap().is_none());
    let s1 = file.header().sample("S1").unwrap();
    let s2 = file.header().sample("S2").unwrap();
    let samples = rec.samples();
    assert_eq!(samples.sample(s1).get(&gq).unwrap(), Some(30));
    assert_eq!(samples.sample(s2).get(&gq).unwrap(), None);
    assert_eq!(samples.sample(s1).get(&zz).unwrap(), Some("q"));
    // r[verify vcf_values.rlen+3]
    // An undeclared END is VCF's Integer END, as tabix reads it.
    assert_eq!(rec.locus().end(), VcfPos::At(seqair_types::Pos1::new(20).unwrap()));

    // Every field the record holds is listed with its table definition.
    let defs: Vec<(String, ValueType, Number)> = info
        .iter()
        .map(|f| f.map(|(name, def, _)| (name.to_owned(), def.typ, def.number)).unwrap())
        .collect();
    assert_eq!(
        defs,
        [
            ("SVTYPE".into(), ValueType::String, Number::Count(1)),
            ("SVLEN".into(), ValueType::Integer, Number::AlternateBases),
            ("END".into(), ValueType::Integer, Number::Count(1)),
            ("DB".into(), ValueType::Flag, Number::Count(0)),
            ("TRID".into(), ValueType::String, Number::Unknown),
            ("AN".into(), ValueType::Integer, Number::Count(1)),
        ]
    );
}

// r[verify vcf_reader.field_table]
// r[verify vcf_values.file_keys+3]
// r[verify vcf_values.lenient_keys]
/// One header as VCF text and as BCF: plain keys resolve its fields only,
/// lenient keys the same undeclared names to the same definitions, and a
/// BCF record, which names fields by dictionary index, never holds one.
#[test]
fn text_and_bcf_resolve_the_same_keys() {
    let header = file(&[]).header().clone();
    let files: Vec<VcfFile> = [OutputFormat::Vcf, OutputFormat::Bcf]
        .into_iter()
        .map(|format| {
            let (bytes, _) =
                Writer::new(Vec::new(), format).write_header(&header).unwrap().finish().unwrap();
            VcfFile::from_bytes(bytes, None).unwrap()
        })
        .collect();
    for file in &files {
        let unknown =
            |name: &str| KeyError::UnknownField { name: name.into(), kind: KeyKind::Info };
        assert_eq!(file.info_key::<Str>("SVTYPE").unwrap_err(), unknown("SVTYPE"));
        assert_eq!(file.info_key_opt::<Str>("TRID").unwrap(), None);
        assert_eq!(file.format_key_opt::<Scalar<i32>>("GQ").unwrap(), None);
        let lenient = file.lenient();
        let svlen: InfoInts = lenient.info_key("SVLEN").unwrap();
        assert_eq!(svlen.id().number(), Number::AlternateBases);
        assert!(lenient.info_key_opt::<Str>("TRID").unwrap().is_some());
        assert!(lenient.format_key_opt::<Scalar<i32>>("GQ").unwrap().is_some());
        assert_eq!(lenient.info_key_opt::<Scalar<i32>>("TRID").unwrap(), None);
    }
}

// r[verify vcf_values.header_brand+2]
/// A key the header does not declare reads its file, and every writer
/// refuses it, even one with the file's header: that header has no id for
/// the field.
#[test]
fn table_keys_read_their_file_and_no_writer_takes_them() {
    let file = file(&["chr1\t5\t.\tA\tC\t.\t.\tTRID=t1\t.\t.\t."]);
    let trid: InfoString = file.lenient().info_key("TRID").unwrap();
    let svtype: InfoKey<Str> = file.lenient().info_key("SVTYPE").unwrap();
    let mut out = Vec::new();
    let mut writer = Writer::new(&mut out, OutputFormat::Vcf).write_header(file.header()).unwrap();
    let chr1 = file.header().contig("chr1").unwrap();
    let alleles = seqair::vcf::Alleles::snv(seqair_types::Base::A, seqair_types::Base::C).unwrap();
    let pos = seqair_types::Pos1::new(5).unwrap();
    let mut enc = writer.begin_record(&chr1, pos, &alleles, None).unwrap().filter_pass();
    assert!(matches!(trid.encode(&mut enc, "x"), Err(VcfError::ForeignKey { .. })));
    assert!(matches!(svtype.encode(&mut enc, "DEL"), Err(VcfError::ForeignKey { .. })));
    // Another file's table keys are foreign to this one.
    let other = self::file(&[]);
    let foreign: InfoString = other.lenient().info_key("TRID").unwrap();
    let info = file.reader::<seqair::vcf::Info>().open().unwrap().read().unwrap().unwrap();
    assert_eq!(info.get(&trid).unwrap(), Some("t1"));
    assert!(matches!(info.get(&foreign), Err(seqair::vcf::VcfReadError::ForeignKey { .. })));
}

/// Keys resolved on one thread and records read on another agree on ids.
#[test]
fn the_table_is_shared_by_threads() {
    let lines: Vec<String> =
        (1..200).map(|i| format!("chr1\t{i}\t.\tA\tC\t.\t.\tK{}=v{i}\t.\t.\t.", i % 20)).collect();
    let lines: Vec<&str> = lines.iter().map(String::as_str).collect();
    let file = file(&lines);
    std::thread::scope(|s| {
        let readers: Vec<_> = (0..4)
            .map(|t| {
                let file = &file;
                s.spawn(move || {
                    // Each thread resolves the keys in another order.
                    let keys: std::collections::HashMap<u32, InfoString> = (0..20)
                        .map(|k| (k + t * 7) % 20)
                        .map(|k| (k, file.lenient().info_key(&format!("K{k}")).unwrap()))
                        .collect();
                    let mut n = 0;
                    let mut reader = file.reader::<(Locus, seqair::vcf::Info)>().open().unwrap();
                    while let Some((locus, info)) = reader.next_record().unwrap() {
                        let pos = locus.start().pos1().unwrap().as_u32();
                        let value = info.get(&keys[&(pos % 20)]).unwrap();
                        assert_eq!(value, Some(format!("v{pos}").as_str()));
                        n += 1;
                    }
                    n
                })
            })
            .collect();
        for r in readers {
            assert_eq!(r.join().unwrap(), 199);
        }
    });
}

// ── htslib as the oracle ────────────────────────────────────────────────

/// One random record with INFO and FORMAT keys a header leaves out: VCF
/// reserved ones (AN, AF, CIPOS, SVTYPE; GQ) and others (TRID, XX; ZZ).
fn line(tc: &TestCase, pos: u32) -> String {
    let int = |tc: &TestCase| -> String {
        if tc.draw(gs::integers::<u8>().max_value(4)) == 0 {
            ".".into()
        } else {
            tc.draw(gs::integers::<i32>().min_value(-1000).max_value(1000)).to_string()
        }
    };
    let float = |tc: &TestCase| -> String {
        tc.draw(gs::sampled_from(vec![".", "0.5", "0.25", "1.5", "3", "-2"])).to_owned()
    };
    let word = |tc: &TestCase| -> String {
        tc.draw(gs::sampled_from(vec![".", "t1", "DEL", "x_y", "7"])).to_owned()
    };
    let mut info = Vec::new();
    if tc.draw(gs::booleans()) {
        info.push(format!("DP={}", tc.draw(gs::integers::<u16>())));
    }
    if tc.draw(gs::booleans()) {
        info.push(format!("AN={}", int(tc)));
    }
    if tc.draw(gs::booleans()) {
        let n = tc.draw(gs::integers::<usize>().min_value(1).max_value(3));
        info.push(format!("AF={}", (0..n).map(|_| float(tc)).collect::<Vec<_>>().join(",")));
    }
    if tc.draw(gs::booleans()) {
        info.push(format!("CIPOS={},{}", int(tc), int(tc)));
    }
    if tc.draw(gs::booleans()) {
        info.push(format!("SVTYPE={}", word(tc)));
    }
    if tc.draw(gs::booleans()) {
        info.push(format!("TRID={}", word(tc)));
    }
    if tc.draw(gs::booleans()) {
        info.push(format!("XX={},{}", int(tc), word(tc)));
    }
    let info = if info.is_empty() { ".".to_owned() } else { info.join(";") };
    let sample = |tc: &TestCase| format!("0/1:{}:{}", int(tc), word(tc));
    format!("chr1\t{pos}\t.\tA\tC\t.\t.\t{info}\tGT:GQ:ZZ\t{}\t{}", sample(tc), sample(tc))
}

/// A value as VCF text writes it, `.` for missing.
fn render(value: &DynValue<'_>) -> String {
    let join = |parts: Vec<String>| if parts.is_empty() { ".".to_owned() } else { parts.join(",") };
    match value {
        DynValue::Int(v) => {
            join(v.clone().map(|e| e.map_or(".".into(), |e| e.to_string())).collect())
        }
        DynValue::Float(v) => {
            join(v.clone().map(|e| e.map_or(".".into(), |e| e.to_string())).collect())
        }
        DynValue::Str(s) => s.unwrap_or(".").to_owned(),
        DynValue::Flag => "flag".to_owned(),
        DynValue::Gt(_) => "gt".to_owned(),
    }
}

/// htslib's strings, joined as written.
fn joined(parts: &[&[u8]]) -> String {
    let parts: Vec<&str> = parts.iter().map(|p| std::str::from_utf8(p).unwrap()).collect();
    parts.join(",")
}

// r[verify vcf_reader.field_table]
// r[verify vcf_values.dynamic+2]
/// htslib reads a key its header leaves out as a String; seqair reads a
/// reserved one with VCF's type, every other one as a String. Rendered as
/// text, every value of every record agrees. (bcftools cannot convert such
/// a file: it refuses to write a record whose tag it had to add, so htslib
/// reading the text itself is the oracle.)
#[hegel::test(test_cases = pinned::cases(40))]
fn undeclared_keys_read_what_htslib_reads(tc: TestCase) {
    let n = tc.draw(gs::integers::<u32>().min_value(1).max_value(8));
    let lines: Vec<String> = (1..=n).map(|i| line(&tc, i * 10)).collect();
    let lines: Vec<&str> = lines.iter().map(String::as_str).collect();
    let mut tmp = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
    tmp.write_all(text(&lines).as_bytes()).unwrap();
    tmp.flush().unwrap();

    let file = VcfFile::open(tmp.path()).unwrap();
    let ours: Vec<FullRecord> =
        file.reader::<FullRecord>().open().unwrap().into_iter().collect::<Result<_, _>>().unwrap();
    let mut htslib = hts::Reader::from_path(tmp.path()).unwrap();
    let theirs: Vec<hts::Record> = htslib.records().collect::<Result<_, _>>().unwrap();
    assert_eq!(ours.len(), theirs.len());
    let samples: Vec<_> = file.header().sample_indices().collect();
    for (ours, theirs) in ours.iter().zip(&theirs) {
        for field in ours.info().iter() {
            let (name, _, value) = field.unwrap();
            if name == "DP" {
                continue; // declared: the other oracles cover it
            }
            let strings = theirs.info(name.as_bytes()).string().unwrap().expect("present");
            assert_eq!(render(&value), joined(&strings), "INFO {name}");
        }
        for (s, &sample) in samples.iter().enumerate() {
            for field in ours.samples().sample(sample).iter() {
                let (name, _, value) = field.unwrap();
                if name == "GT" {
                    continue;
                }
                let strings = theirs.format(name.as_bytes()).string().unwrap();
                assert_eq!(render(&value), joined(&[strings[s]]), "FORMAT {name} sample {s}");
            }
        }
    }
}
