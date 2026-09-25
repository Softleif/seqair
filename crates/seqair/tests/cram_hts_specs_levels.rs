//! The hts-specs CRAM 3.1 conformance files: `level-1.cram` to `level-4.cram`
//! hold the same 20000 records, each level adding codecs — the name tokeniser
//! and rANS Nx16 transforms at level 2, bzip2 and fqzcomp at level 3, lzma and
//! the arithmetic coder (also inside tok3) at level 4. seqair must decode every
//! level to exactly the records htslib decodes.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "test code"
)]
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "test code with known small values"
)]

use rust_htslib::bam::{self, Read as _};
use seqair::bam::{Pos0, RecordStore};
use seqair::reader::Readers;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

const CONTIG: &str = "chrM";
const CONTIG_LEN: u32 = 16_571;

fn data_dir() -> PathBuf {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/hts-specs/cram")).to_path_buf()
}

fn level_path(level: u8) -> PathBuf {
    data_dir().join(format!("3.1/passed/level-{level}.cram"))
}

/// The files embed their reference; this FASTA only satisfies the reader's
/// signature and is never consulted for `chrM`.
fn fasta() -> PathBuf {
    data_dir().join("ce.fa")
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Decoded {
    qname: Vec<u8>,
    pos: i64,
    flags: u16,
    mapq: u8,
    cigar: String,
    seq: Vec<u8>,
    qual: Vec<u8>,
    aux: Vec<u8>,
}

/// Copy the file next to a fresh CRAI, since the vendored files ship without.
fn indexed_copy(dir: &Path, level: u8) -> PathBuf {
    let cram = dir.join(format!("level-{level}.cram"));
    std::fs::copy(level_path(level), &cram).expect("copy");
    let status = Command::new("samtools")
        .arg("index")
        .arg(&cram)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .expect("samtools not found");
    assert!(status.success(), "samtools index level-{level}");
    cram
}

fn read_seqair(cram: &Path) -> Vec<Decoded> {
    let mut readers = Readers::open(cram, &fasta()).expect("open");
    let tid = readers.header().tid(CONTIG).expect("contig");
    let mut store = RecordStore::new();
    readers
        .fetch_into(
            tid,
            (Pos0::new(0).unwrap()..=Pos0::new(CONTIG_LEN - 1).unwrap()).into(),
            &mut store,
        )
        .expect("fetch");
    store
        .indices()
        .map(|idx| {
            let rec = store.record(idx).unwrap();
            Decoded {
                qname: rec.qname().to_vec(),
                pos: rec.pos.as_i64(),
                flags: rec.flags.raw(),
                mapq: rec.mapq,
                cigar: rec.cigar().iter().map(ToString::to_string).collect(),
                seq: rec.seq().iter().map(|b| **b).collect(),
                qual: seqair_types::BaseQuality::slice_to_bytes(rec.qual()).to_vec(),
                aux: rec.aux().to_vec(),
            }
        })
        .collect()
}

/// htslib's view of every record, in file order — all are placed on `chrM`,
/// the 1178 unmapped ones next to their mates.
fn read_htslib(cram: &Path) -> Vec<Decoded> {
    let mut reader = bam::Reader::from_path(cram).expect("htslib open");
    // htslib regenerates MD and NM from the reference unless told not to;
    // seqair returns the stored tags only.
    // SAFETY: `htsfile()` is the reader's live, open `htsFile`, and
    // `CRAM_OPT_DECODE_MD` takes one `int` argument.
    let rc = unsafe {
        rust_htslib::htslib::hts_set_opt(
            reader.htsfile(),
            rust_htslib::htslib::hts_fmt_option_CRAM_OPT_DECODE_MD,
            0 as std::ffi::c_int,
        )
    };
    assert_eq!(rc, 0, "hts_set_opt(CRAM_OPT_DECODE_MD)");
    let mut out = Vec::new();
    let mut record = bam::Record::new();
    while let Some(r) = reader.read(&mut record) {
        r.expect("htslib read");
        let inner = record.inner();
        // SAFETY: htslib keeps `l_data` initialised bytes at `data` for the
        // record it just read, and `record` outlives this borrow.
        let data = unsafe { std::slice::from_raw_parts(inner.data, inner.l_data as usize) };
        let aux_start = usize::from(inner.core.l_qname)
            + record.cigar_len() * 4
            + record.seq_len().div_ceil(2)
            + record.seq_len();
        out.push(Decoded {
            qname: record.qname().to_vec(),
            pos: record.pos(),
            flags: record.flags(),
            mapq: record.mapq(),
            cigar: record.cigar().to_string(),
            seq: record.seq().as_bytes(),
            qual: record.qual().to_vec(),
            aux: data[aux_start..].to_vec(),
        });
    }
    out
}

/// The block methods `samtools cram-size -v` reports, so the test states
/// which codecs each level exercises instead of trusting the README.
fn block_methods(cram: &Path) -> Vec<String> {
    let out = Command::new("samtools")
        .args(["cram-size", "-v"])
        .arg(cram)
        .output()
        .expect("samtools cram-size");
    assert!(out.status.success(), "samtools cram-size");
    String::from_utf8(out.stdout)
        .expect("ASCII")
        .lines()
        .filter(|l| l.starts_with("BLOCK"))
        .filter_map(|l| l.split_whitespace().nth(5).map(str::to_owned))
        .collect()
}

/// The file's records as SAM text, from the samtools on `PATH`.
fn samtools_sam(cram: &Path) -> Vec<u8> {
    let out = Command::new("samtools")
        .args(["view", "--no-PG", "--input-fmt-option", "decode_md=0"])
        .arg(cram)
        .output()
        .expect("samtools view");
    assert!(out.status.success(), "samtools view {}", cram.display());
    out.stdout
}

/// Level `level` decodes to htslib's records. The rust-htslib we link bundles
/// htslib 1.19, which fails on `level-4.cram` (`BamTruncatedRecord`) while
/// current samtools reads it; so the records come from htslib's decode of
/// `level-1.cram`, and samtools confirms that level `level` holds exactly
/// those records.
fn assert_level_matches_htslib(level: u8) {
    let dir = tempfile::tempdir().expect("tempdir");
    let cram = indexed_copy(dir.path(), level);
    let baseline = if level == 1 { cram.clone() } else { indexed_copy(dir.path(), 1) };
    if level != 1 {
        assert!(
            samtools_sam(&cram) == samtools_sam(&baseline),
            "samtools: level-{level} and level-1 hold different records"
        );
    }

    let ours = read_seqair(&cram);
    let theirs = read_htslib(&baseline);
    assert_eq!(ours.len(), theirs.len(), "level-{level}: record count");
    assert_eq!(ours.len(), 20_000, "level-{level}: every record");
    for (i, (a, b)) in ours.iter().zip(&theirs).enumerate() {
        assert_eq!(a, b, "level-{level}: record {i} differs from htslib");
    }
}

// r[verify cram.codec.rans_nx16]
// r[verify cram.slice.embedded_ref]
// r[verify cram.edge.missing_reference+2]
// r[verify cram.record.cf_tag]
#[test]
fn level_1_matches_htslib() {
    assert_level_matches_htslib(1);
}

// r[verify cram.codec.tok3]
#[test]
fn level_2_matches_htslib() {
    let methods = block_methods(&level_path(2));
    assert!(methods.iter().any(|m| m == "tok3-rans"), "{methods:?}");
    assert_level_matches_htslib(2);
}

// r[verify cram.codec.fqzcomp]
// r[verify cram.codec.bzip2]
#[test]
fn level_3_matches_htslib() {
    let methods = block_methods(&level_path(3));
    assert!(methods.iter().any(|m| m == "fqzcomp"), "{methods:?}");
    assert_level_matches_htslib(3);
}

// r[verify cram.codec.arith]
// r[verify cram.codec.lzma]
#[test]
fn level_4_matches_htslib() {
    let methods = block_methods(&level_path(4));
    assert!(methods.iter().any(|m| m == "fqzcomp"), "{methods:?}");
    assert!(methods.iter().any(|m| m == "tok3-arith"), "{methods:?}");
    assert!(methods.iter().any(|m| m.starts_with("arith-")), "{methods:?}");
    assert!(methods.iter().any(|m| m == "lzma"), "{methods:?}");
    assert_level_matches_htslib(4);
}
