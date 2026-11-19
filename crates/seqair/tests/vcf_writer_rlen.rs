//! The writer's reference length (`r[record_encoder.rlen]`): records whose
//! span INFO END, the SVLEN of symbolic SV alleles or FORMAT LEN of gVCF
//! blocks extend, written by seqair as VCF.gz and BCF. bcftools is the
//! oracle twice over: the `rlen` it writes converting seqair's VCF text to
//! BCF is the one seqair's BCF must store, and the index `bcftools index`
//! builds from seqair's output is the one seqair must co-produce.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::format_push_string,
    clippy::format_collect,
    reason = "test code"
)]

#[path = "support/index_parse.rs"]
mod index_parse;
#[path = "../src/pinned.rs"]
mod pinned;

use hegel::prelude::*;
use index_parse::{assert_same_index, parse_index};
use seqair::vcf::record_encoder::{Arr, FormatFieldDef, Gt, InfoFieldDef, OptArr, Scalar};
use seqair::vcf::{
    Alleles, ContigDef, ContigId, FormatGt, FormatInts, Genotype, IndexKind, InfoInt, InfoIntOpts,
    Locus, Number, OutputFormat, VcfFile, VcfHeader, Writer,
};
use std::path::{Path, PathBuf};
use std::process::Command;

fn has_bcftools() -> bool {
    Command::new("bcftools").arg("--version").output().is_ok()
}

const ALTS: [&str; 10] =
    ["A", "CT", "<DEL>", "<DUP>", "<INV>", "<CNV:TR>", "<INS>", "<*>", "<NON_REF>", "*"];

#[derive(Debug, Clone)]
struct Rec {
    pos: u32,
    reference: &'static str,
    alts: Vec<&'static str>,
    /// Every END encoded, in order: the record keeps the last.
    ends: Vec<i32>,
    svlen: Option<Vec<Option<i32>>>,
    len: Option<[Option<i32>; 2]>,
}

struct Setup {
    header: VcfHeader,
    chr1: ContigId,
    end: InfoInt,
    svlen: InfoIntOpts,
    gt: FormatGt,
    len: FormatInts,
}

/// The contig is 100 Mb so that `bcftools index -c`, which sizes the CSI
/// from it starting at depth 0, lands on seqair's depth 5.
fn setup() -> Setup {
    let mut builder = VcfHeader::builder();
    let chr1 = builder.register_contig("chr1", ContigDef::new(Some(100_000_000))).unwrap();
    let mut builder = builder.infos();
    let end = builder.register_info(&InfoFieldDef::<Scalar<i32>>::new("END", "End")).unwrap();
    let svlen = builder
        .register_info(&InfoFieldDef::<OptArr<i32>>::new("SVLEN", Number::AlternateBases, "Length"))
        .unwrap();
    let mut builder = builder.formats();
    let gt = builder.register_format(&FormatFieldDef::<Gt>::new("GT", "Genotype")).unwrap();
    let len = builder
        .register_format(&FormatFieldDef::<Arr<i32>>::new("LEN", Number::Count(1), "Block length"))
        .unwrap();
    let mut builder = builder.samples();
    builder.add_sample("S1").unwrap();
    builder.add_sample("S2").unwrap();
    Setup { header: builder.build().unwrap(), chr1, end, svlen, gt, len }
}

/// Near `pos` (where it decides the end), anywhere, or negative.
fn value(tc: &TestCase, pos: u32) -> i32 {
    let pos = i32::try_from(pos).unwrap();
    match tc.draw(gs::integers::<u8>().max_value(4)) {
        0 => tc.draw(gs::integers::<i32>().min_value(-5000).max_value(5000)),
        1 => tc.draw(gs::integers::<i32>().min_value(-1_000_000).max_value(1_000_000)),
        _ => pos + tc.draw(gs::integers::<i32>().min_value(-20).max_value(3000)),
    }
}

#[hegel::composite]
fn record(tc: &TestCase, pos: u32) -> Rec {
    let reference = &"ACGTAC"[..tc.draw(gs::integers::<usize>().min_value(1).max_value(6))];
    let n_alts = tc.draw(gs::integers::<usize>().max_value(4));
    let alts = (0..n_alts).map(|_| tc.draw(gs::sampled_from(ALTS.to_vec()))).collect::<Vec<_>>();
    let n_ends = tc.draw(gs::integers::<usize>().max_value(2));
    let ends = (0..n_ends).map(|_| value(tc, pos)).collect();
    let svlen = (n_alts > 0 && tc.draw(gs::booleans()))
        .then(|| (0..n_alts).map(|_| tc.draw(gs::booleans()).then(|| value(tc, 0))).collect());
    let len = tc
        .draw(gs::booleans())
        .then(|| [0, 1].map(|_| tc.draw(gs::booleans()).then(|| value(tc, 0))));
    Rec { pos, reference, alts, ends, svlen, len }
}

fn write(dir: &Path, setup: &Setup, records: &[Rec], format: OutputFormat) -> PathBuf {
    let path = dir.join(if format == OutputFormat::Bcf { "ours.bcf" } else { "ours.vcf.gz" });
    let mut writer = Writer::new(std::fs::File::create(&path).unwrap(), format)
        .index_kind(if format == OutputFormat::Bcf { IndexKind::Csi } else { IndexKind::Tbi })
        .write_header(&setup.header)
        .unwrap();
    for r in records {
        let alleles = Alleles::complex(r.reference, r.alts.iter().copied());
        let pos = seqair_types::Pos1::new(r.pos).unwrap();
        let mut enc = writer.begin_record(&setup.chr1, pos, &alleles, None).unwrap().filter_pass();
        for &end in &r.ends {
            setup.end.encode(&mut enc, end).unwrap();
        }
        if let Some(svlen) = &r.svlen {
            setup.svlen.encode(&mut enc, svlen).unwrap();
        }
        let mut enc = enc.begin_samples();
        setup.gt.encode(&mut enc, &[Genotype::unphased(0, 0), Genotype::unphased(0, 0)]).unwrap();
        if let Some(len) = r.len {
            let [a, b] = len.map(|v| v.map_or_else(Vec::new, |v| vec![v]));
            setup.len.encode(&mut enc, &[&a, &b]).unwrap();
        }
        enc.emit().unwrap();
    }
    let (_, index) = writer.finish().unwrap();
    let index = index.unwrap();
    index.write_to_path(format!("{}.{}", path.display(), index.suffix())).unwrap();
    path
}

fn bcftools(args: &[&str]) -> String {
    let out = Command::new("bcftools").args(args).output().unwrap();
    assert!(out.status.success(), "bcftools {args:?}: {}", String::from_utf8_lossy(&out.stderr));
    String::from_utf8(out.stdout).unwrap()
}

/// Every record's 1-based start and inclusive end, as a reader reads them.
fn spans(path: &Path) -> Vec<String> {
    let file = VcfFile::open(path).unwrap();
    let mut reader = file.reader::<Locus>().open().unwrap();
    let mut out = Vec::new();
    while let Some(l) = reader.read().unwrap() {
        out.push(format!("{:?}..{:?}", l.start(), l.end()));
    }
    out
}

fn decompress(path: &Path) -> Vec<u8> {
    let mut reader = seqair::bam::bgzf::BgzfReader::from_reader(std::fs::File::open(path).unwrap());
    let mut data = Vec::new();
    reader.read_to_end(&mut data).unwrap();
    data
}

/// `bcftools index` of `path` into `out`, as `-c` (CSI) or `-t` (TBI).
fn bcftools_index(path: &Path, out: &Path, flag: &str) -> Vec<u8> {
    bcftools(&["index", flag, "-o", out.to_str().unwrap(), path.to_str().unwrap()]);
    decompress(out)
}

// r[verify record_encoder.rlen]
// r[verify vcf_record.alleles_rlen+2]
/// seqair's BCF stores the `rlen` bcftools computes from seqair's VCF text
/// of the same records, and both co-produced indexes (TBI for the VCF.gz,
/// CSI for the BCF) hold every bin and chunk `bcftools index` builds.
#[hegel::test(test_cases = pinned::cases(40))]
fn rlen_and_index_are_what_bcftools_makes_of_them(tc: TestCase) {
    if !has_bcftools() {
        return;
    }
    let n = tc.draw(gs::integers::<usize>().min_value(1).max_value(15));
    let mut pos = 1000;
    let mut records = Vec::new();
    for _ in 0..n {
        pos += tc.draw(gs::integers::<u32>().max_value(5000));
        records.push(tc.draw(record(pos)));
    }
    let setup = setup();
    let dir = tempfile::tempdir().unwrap();
    let vcf = write(dir.path(), &setup, &records, OutputFormat::VcfGz);
    let bcf = write(dir.path(), &setup, &records, OutputFormat::Bcf);

    let theirs = dir.path().join("theirs.bcf");
    bcftools(&[
        "view",
        "--no-version",
        "-Ob",
        "-o",
        theirs.to_str().unwrap(),
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(spans(&bcf), spans(&theirs), "{records:#?}");

    let ours = parse_index(&decompress(&dir.path().join("ours.bcf.csi")), true);
    let csi = parse_index(&bcftools_index(&bcf, &dir.path().join("theirs.csi"), "-c"), true);
    assert_same_index(&ours, &csi, "BCF CSI");
    let ours = parse_index(&decompress(&dir.path().join("ours.vcf.gz.tbi")), false);
    let tbi = parse_index(&bcftools_index(&vcf, &dir.path().join("theirs.tbi"), "-t"), false);
    assert_same_index(&ours, &tbi, "VCF.gz TBI");
}

// r[verify record_encoder.rlen]
/// A deletion and a gVCF block are found by a region that overlaps only their
/// far end, through the index seqair wrote: the deletion by its SVLEN, the
/// block by END and by FORMAT LEN.
#[test]
fn bcftools_finds_an_sv_by_its_end() {
    if !has_bcftools() {
        return;
    }
    let setup = setup();
    let records = [
        Rec {
            pos: 1000,
            reference: "A",
            alts: vec!["<DEL>"],
            ends: vec![51_000],
            svlen: Some(vec![Some(-50_000)]),
            len: None,
        },
        Rec {
            pos: 2000,
            reference: "C",
            alts: vec!["<*>"],
            ends: vec![],
            svlen: None,
            len: Some([Some(70_000), None]),
        },
        Rec {
            pos: 3000,
            reference: "G",
            alts: vec!["<*>"],
            ends: vec![90_000],
            svlen: None,
            len: None,
        },
    ];
    let dir = tempfile::tempdir().unwrap();
    for format in [OutputFormat::VcfGz, OutputFormat::Bcf] {
        let path = write(dir.path(), &setup, &records, format);
        let found = |region: &str| -> Vec<String> {
            bcftools(&["view", "-H", "-r", region, path.to_str().unwrap()])
                .lines()
                .map(|l| l.split('\t').nth(1).unwrap().to_owned())
                .collect()
        };
        assert_eq!(found("chr1:50000-50001"), ["1000", "2000", "3000"], "{format:?}");
        assert_eq!(found("chr1:60000-60001"), ["2000", "3000"], "{format:?}");
        assert_eq!(found("chr1:80000-80001"), ["3000"], "{format:?}");
        assert!(found("chr1:95000-95001").is_empty(), "{format:?}");
    }
}
