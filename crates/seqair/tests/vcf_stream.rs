//! `VcfStream` (`r[vcf_reader.stream]`): VCF, BGZF VCF and BCF from input
//! that cannot seek read as the same file read through `VcfFile`.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, reason = "test code")]

use seqair::vcf::{
    FullRecord, Locus, OpenError, PlanError, Sample, Site, VcfFile, VcfReadError, VcfStream,
};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

const BCFTOOLS: &str = "/opt/homebrew/bin/bcftools";

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/vcf").join(name)
}

/// Every record of `reader`: where it is and its bytes as stored.
fn records(reader: seqair::vcf::Reader<FullRecord>) -> Vec<(String, Vec<u8>)> {
    reader
        .into_iter()
        .map(|r| {
            let r = r.unwrap();
            let l = r.locus();
            (format!("{}:{:?}-{:?}", l.contig().name(), l.start(), l.end()), r.raw_bytes().to_vec())
        })
        .collect()
}

// r[verify vcf_reader.stream]
/// Each fixture through a `VcfStream` (on and off a thread pool) reads as
/// through `VcfFile`.
#[test]
fn a_stream_reads_as_the_file() {
    for name in ["rastair.bcf", "rastair.vcf.gz", "truth.vcf.gz", "panel.bcf"] {
        let path = fixture(name);
        let want = records(VcfFile::open(&path).unwrap().reader::<FullRecord>().open().unwrap());
        assert!(!want.is_empty());
        let stream = VcfStream::new(std::fs::File::open(&path).unwrap());
        assert_eq!(records(stream.reader::<FullRecord>().open().unwrap()), want, "{name}");
        let two = NonZeroUsize::new(2).unwrap();
        let stream = VcfStream::new(std::fs::File::open(&path).unwrap());
        assert_eq!(
            records(stream.reader::<FullRecord>().threads(two).open().unwrap()),
            want,
            "{name} threads"
        );
    }
}

// r[verify vcf_reader.stream]
/// `bcftools view` piped into a `VcfStream`, as BCF, BGZF VCF and plain
/// VCF, reads as the file it printed.
#[test]
fn a_pipe_reads_as_the_file() {
    let bin = if Path::new(BCFTOOLS).exists() { BCFTOOLS } else { "bcftools" };
    if Command::new(bin).arg("--version").output().is_err() {
        eprintln!("bcftools not found, skipping");
        return;
    }
    let path = fixture("rastair.vcf.gz");
    let want: Vec<(String, String)> =
        sites(VcfFile::open(&path).unwrap().reader::<(Locus, Site)>().open().unwrap());
    for out in ["-Ob", "-Oz", "-Ov"] {
        let mut child = Command::new(bin)
            .args(["view", "--no-version", out, path.to_str().unwrap()])
            .stdout(Stdio::piped())
            .spawn()
            .unwrap();
        let stdout = child.stdout.take().unwrap();
        let got = sites(VcfStream::new(stdout).reader::<(Locus, Site)>().open().unwrap());
        assert!(child.wait().unwrap().success());
        assert_eq!(got, want, "bcftools view {out}");
    }
}

fn sites(reader: seqair::vcf::Reader<(Locus, Site)>) -> Vec<(String, String)> {
    reader
        .into_iter()
        .map(|r| {
            let (locus, site) = r.unwrap();
            (
                format!("{}:{:?}", locus.contig().name(), locus.start()),
                format!("{:?}", site.alleles().reference()),
            )
        })
        .collect()
}

/// A stream's options behave as a file reader's; gzip that is not BGZF is
/// refused.
#[test]
fn stream_options_and_errors() {
    let panel = || VcfStream::new(std::fs::File::open(fixture("panel.bcf")).unwrap());
    let mut reader = panel().reader::<Sample>().sample("NA12878").open().unwrap();
    assert_eq!(reader.next_record().unwrap().unwrap().name(), "NA12878");
    assert!(matches!(
        panel().reader::<Sample>().open(),
        Err(VcfReadError::Plan(PlanError::SampleNotSelected { n: 4 }))
    ));
    let mut gz = Vec::new();
    {
        use std::io::Write;
        let mut enc = flate2::write::GzEncoder::new(&mut gz, flate2::Compression::default());
        enc.write_all(b"##fileformat=VCFv4.2\n").unwrap();
        enc.finish().unwrap();
    }
    assert!(matches!(
        VcfStream::new(std::io::Cursor::new(gz)).reader::<Locus>().open(),
        Err(VcfReadError::Open { path: None, source: OpenError::NotBgzf })
    ));
    assert!(matches!(
        VcfStream::new(std::io::Cursor::new(b"hello".to_vec())).reader::<Locus>().open(),
        Err(VcfReadError::Open { path: None, source: OpenError::UnknownFormat })
    ));
    // A stream reader has no index to requery.
    let mut reader = panel().reader::<Locus>().open().unwrap();
    assert!(matches!(
        reader.requery(seqair::vcf::Region::contig("chr22")),
        Err(VcfReadError::NoIndex)
    ));
}

// r[verify vcf_values.reader_keys]
/// A stream reader resolves keys as its file would, and they read its
/// records as the file's keys read the file's.
#[test]
fn a_stream_readers_keys_read_its_records() {
    let path = fixture("rastair.bcf");
    let file = VcfFile::open(&path).unwrap();
    let dp: seqair::vcf::InfoInt = file.info_key("DP").unwrap();
    let want: Vec<Option<i32>> = file
        .reader::<seqair::vcf::Info>()
        .open()
        .unwrap()
        .into_iter()
        .map(|info| info.unwrap().get(&dp).unwrap())
        .collect();

    let mut reader = VcfStream::new(std::fs::File::open(&path).unwrap())
        .reader::<seqair::vcf::Info>()
        .open()
        .unwrap();
    let keys = reader.keys();
    let dp: seqair::vcf::InfoInt = keys.info_key("DP").unwrap();
    assert_eq!(keys.info_key_opt::<seqair::vcf::Str>("NOPE").unwrap(), None);
    assert!(keys.lenient().info_key_opt::<seqair::vcf::Str>("NOPE").unwrap().is_some());
    assert!(matches!(
        keys.info_key::<seqair::vcf::Str>("DP"),
        Err(seqair::vcf::KeyError::TypeMismatch { .. })
    ));
    let mut got = Vec::new();
    while let Some(info) = reader.next_record().unwrap() {
        got.push(info.get(&dp).unwrap());
    }
    assert!(got.iter().any(Option::is_some));
    assert_eq!(got, want);
}
