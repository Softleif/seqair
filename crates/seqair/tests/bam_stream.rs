//! `BamStreamReader`, `RawBamRecord` and `BamWriter::write_raw`, against
//! htslib and noodles.
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
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss,
    reason = "test code with known small values"
)]

#[path = "../src/pinned.rs"]
mod pinned;

use hegel::prelude::*;
use seqair::bam::RawBamRecord;
use seqair::bam::stream::{BamStreamError, BamStreamReader};
use seqair::io::BgzfWriter;
use std::fmt::Write as _;
use std::io::Cursor;

const REFS: [(&str, u32); 2] = [("chr1", 100_000), ("chr2", 50_000)];

fn header_text() -> String {
    let mut text = String::from("@HD\tVN:1.6\tSO:unknown\n");
    for (name, len) in REFS {
        writeln!(text, "@SQ\tSN:{name}\tLN:{len}").unwrap();
    }
    text
}

/// The uncompressed BAM header for [`REFS`].
fn header_bytes() -> Vec<u8> {
    let text = header_text();
    let mut out = b"BAM\x01".to_vec();
    out.extend_from_slice(&(text.len() as i32).to_le_bytes());
    out.extend_from_slice(text.as_bytes());
    out.extend_from_slice(&(REFS.len() as i32).to_le_bytes());
    for (name, len) in REFS {
        out.extend_from_slice(&(name.len() as i32 + 1).to_le_bytes());
        out.extend_from_slice(name.as_bytes());
        out.push(0);
        out.extend_from_slice(&(len as i32).to_le_bytes());
    }
    out
}

/// A BGZF-compressed BAM of the header and `payloads` (each framed with its
/// `block_size`), one block per record so a bad one cannot hide another.
fn bam_of_payloads(payloads: &[Vec<u8>]) -> Vec<u8> {
    let mut writer = BgzfWriter::new(Vec::new());
    writer.write_all(&header_bytes()).unwrap();
    for payload in payloads {
        writer.flush_if_needed(usize::MAX).unwrap();
        writer.write_all(&(payload.len() as i32).to_le_bytes()).unwrap();
        writer.write_all(payload).unwrap();
    }
    writer.finish().unwrap()
}

/// The fields of one BAM record, packed by hand so nothing normalises them.
#[derive(Debug, Clone)]
struct Wire {
    tid: i32,
    pos: i32,
    flag: u16,
    qname: Vec<u8>,
    cigar: Vec<u32>,
    seq_len: usize,
    aux: Vec<u8>,
}

impl Wire {
    fn payload(&self) -> Vec<u8> {
        let mut out = Vec::new();
        out.extend_from_slice(&self.tid.to_le_bytes());
        out.extend_from_slice(&self.pos.to_le_bytes());
        out.push(self.qname.len() as u8 + 1);
        out.push(60);
        out.extend_from_slice(&4680u16.to_le_bytes());
        out.extend_from_slice(&(self.cigar.len() as u16).to_le_bytes());
        out.extend_from_slice(&self.flag.to_le_bytes());
        out.extend_from_slice(&(self.seq_len as u32).to_le_bytes());
        out.extend_from_slice(&(-1i32).to_le_bytes());
        out.extend_from_slice(&(-1i32).to_le_bytes());
        out.extend_from_slice(&0i32.to_le_bytes());
        out.extend_from_slice(&self.qname);
        out.push(0);
        for op in &self.cigar {
            out.extend_from_slice(&op.to_le_bytes());
        }
        out.extend(std::iter::repeat_n(0x12, self.seq_len.div_ceil(2)));
        out.extend(std::iter::repeat_n(30, self.seq_len));
        out.extend_from_slice(&self.aux);
        out
    }
}

/// An htslib reader that exits non-zero on a bad record: htslib's own
/// `test/test_view` when `SEQAIR_HTSLIB_TEST_VIEW` names one (build it from
/// a current htslib checkout), else the samtools on `PATH`.
fn htslib_reader() -> (std::process::Command, bool) {
    match std::env::var_os("SEQAIR_HTSLIB_TEST_VIEW") {
        Some(path) => (std::process::Command::new(path), true),
        None => {
            let mut cmd = std::process::Command::new("samtools");
            cmd.args(["view", "-o", "/dev/null"]);
            (cmd, false)
        }
    }
}

/// Whether the samtools on `PATH` is htslib 1.25 or later, which checks a
/// `CG` CIGAR against the sequence (htslib `cc6b1f7c`); 1.24 does not.
fn samtools_checks_cg_cigars() -> bool {
    let out = std::process::Command::new("samtools").arg("--version").output().unwrap();
    let version = String::from_utf8(out.stdout).unwrap();
    let (major, minor) = version
        .lines()
        .find_map(|l| l.strip_prefix("Using htslib "))
        .and_then(|v| v.split_once('.'))
        .map(|(a, b)| (a.parse::<u32>().unwrap(), b.trim().parse::<u32>().unwrap_or(0)))
        .unwrap();
    (major, minor) >= (1, 25)
}

/// For each record, whether htslib's `bam_read1` accepts it. One file per
/// record: htslib stops at the first bad one.
fn htslib_accepts(payloads: &[Vec<u8>]) -> Vec<bool> {
    payloads
        .iter()
        .map(|payload| {
            let file = tempfile::NamedTempFile::new().unwrap();
            std::fs::write(file.path(), bam_of_payloads(std::slice::from_ref(payload))).unwrap();
            let (mut cmd, _) = htslib_reader();
            cmd.arg(file.path())
                .stdout(std::process::Stdio::null())
                .stderr(std::process::Stdio::null())
                .status()
                .expect("htslib reader not found")
                .success()
        })
        .collect()
}

/// For each record of `bam`, whether seqair's stream reader accepts it.
fn seqair_accepts(bam: &[u8]) -> Vec<bool> {
    let mut reader = BamStreamReader::new(Cursor::new(bam)).unwrap();
    let mut record = RawBamRecord::default();
    let mut accepted = Vec::new();
    loop {
        match reader.read_record(&mut record) {
            Ok(true) => accepted.push(true),
            Ok(false) => return accepted,
            Err(BamStreamError::InvalidRecord { .. }) => accepted.push(false),
            Err(e) => panic!("stream failed: {e:?}"),
        }
    }
}

fn op(len: u32, code: u32) -> u32 {
    len << 4 | code
}

/// A record that may carry the long-CIGAR placeholder: a first op
/// `<seq_len>S` (or not), a `CG` tag of any array type and length, mapped or
/// not, placed or not.
#[hegel::composite]
fn cg_record(tc: &TestCase) -> Wire {
    let seq_len = tc.draw(gs::integers::<usize>().max_value(6));
    let ops = || gs::integers::<u32>().max_value(8);
    let lens = || gs::integers::<u32>().max_value(7);
    let first_len = if tc.draw(gs::booleans()) { seq_len as u32 } else { tc.draw(lens()) };
    let first_code = if tc.draw(gs::booleans()) { 4 } else { tc.draw(ops()) };
    let mut cigar = vec![op(first_len, first_code)];
    for _ in 0..tc.draw(gs::integers::<usize>().max_value(2)) {
        cigar.push(op(tc.draw(lens()), tc.draw(ops())));
    }
    let mut aux = Vec::new();
    if tc.draw(gs::booleans()) {
        let sub = tc.draw(gs::sampled_from(vec![b'I', b'i', b'S', b'c']));
        let width = match sub {
            b'I' | b'i' => 4,
            b'S' => 2,
            _ => 1,
        };
        let n = tc.draw(gs::integers::<usize>().max_value(4));
        aux.extend_from_slice(b"CGB");
        aux.push(sub);
        aux.extend_from_slice(&(n as u32).to_le_bytes());
        for _ in 0..n {
            let value = op(tc.draw(lens()), tc.draw(ops()));
            aux.extend_from_slice(&value.to_le_bytes()[..width]);
        }
    }
    Wire {
        tid: tc.draw(gs::sampled_from(vec![-1, 0])),
        pos: tc.draw(gs::sampled_from(vec![-1, 100])),
        flag: tc.draw(gs::sampled_from(vec![0u16, 4, 16])),
        qname: b"r".to_vec(),
        cigar,
        seq_len,
        aux,
    }
}

// r[verify bam.raw_record.cigar_seq]
/// seqair accepts exactly the CIGAR/sequence combinations htslib's
/// `bam_read1` accepts, around the long-CIGAR placeholder: `CG:B:i` counts
/// as well as `CG:B:I`, a `CG` shorter than the stored CIGAR is no
/// placeholder, and a real CIGAR in `CG` must match the sequence even on an
/// unmapped record.
#[hegel::test(test_cases = pinned::cases(150))]
fn cigar_and_sequence_checks_agree_with_htslib(tc: TestCase) {
    let records: Vec<Wire> = tc.draw(gs::vecs(cg_record()).min_size(1).max_size(8));
    let payloads: Vec<Vec<u8>> = records.iter().map(Wire::payload).collect();
    let seqair = seqair_accepts(&bam_of_payloads(&payloads));
    let htslib = htslib_accepts(&payloads);
    let strict = htslib_reader().1 || samtools_checks_cg_cigars();
    for ((wire, seqair), htslib) in records.iter().zip(seqair).zip(htslib) {
        // Before 1.25, htslib took a CG CIGAR without checking it against the
        // sequence; seqair follows 1.25 (checked with its `test_view`).
        let older_htslib_skips_the_cg_check = !strict && htslib && !seqair && !wire.aux.is_empty();
        assert!(
            seqair == htslib || older_htslib_skips_the_cg_check,
            "{wire:#?}: seqair {seqair}, htslib {htslib}"
        );
    }
}

// r[verify bam.raw_record.cigar_seq]
#[test]
fn a_cg_b_i_placeholder_is_a_placeholder() {
    // 4S2M over 4 bases: no placeholder, the stored CIGAR consumes 6; with
    // CG:B:i holding 4M it is one, as htslib reads it.
    let mut aux = b"CGBi".to_vec();
    aux.extend_from_slice(&2u32.to_le_bytes());
    aux.extend_from_slice(&op(4, 0).to_le_bytes());
    aux.extend_from_slice(&op(5, 3).to_le_bytes());
    let wire = Wire {
        tid: 0,
        pos: 100,
        flag: 0,
        qname: b"r".to_vec(),
        cigar: vec![op(4, 4), op(2, 0)],
        seq_len: 4,
        aux,
    };
    let payloads = [wire.payload()];
    assert_eq!(htslib_accepts(&payloads), [true]);
    assert_eq!(seqair_accepts(&bam_of_payloads(&payloads)), [true]);
}

/// The records of `bam`, all of which must be valid.
fn stream_records(bam: &[u8]) -> Vec<RawBamRecord> {
    let mut reader = BamStreamReader::new(Cursor::new(bam)).unwrap();
    reader.records().collect::<Result<_, _>>().unwrap()
}

// r[verify bam.raw_record]
// r[verify bam.raw_record.cigar_seq]
/// Clearing the unmapped bit on a record whose CIGAR does not cover its
/// sequence would make it one htslib rejects; `set_flags` refuses, and FLAG
/// stays as it was.
#[test]
fn set_flags_keeps_the_cigar_check() {
    let unmapped = Wire {
        tid: 0,
        pos: 100,
        flag: 4,
        qname: b"r".to_vec(),
        cigar: vec![op(2, 0)],
        seq_len: 4,
        aux: Vec::new(),
    };
    let mapped = Wire { flag: 0, ..unmapped.clone() };
    assert_eq!(htslib_accepts(&[unmapped.payload(), mapped.payload()]), [true, false]);

    let mut record = stream_records(&bam_of_payloads(&[unmapped.payload()])).remove(0);
    let flags = record.flags();
    assert!(record.set_flags(seqair::bam::flags::BamFlags::empty()).is_err());
    assert_eq!(record.flags(), flags);
    assert_eq!(record.as_payload(), unmapped.payload());
    let reverse_unmapped = seqair::bam::flags::BamFlags::from(4 | 16);
    record.set_flags(reverse_unmapped).unwrap();
    assert_eq!(record.flags(), reverse_unmapped);
}

// r[verify bam.raw_record]
/// `decode_seq_into` keeps every 4-bit code; it used to decode through the
/// lossy `Base` table, turning `=` and every ambiguity code into `N`.
#[test]
fn decode_seq_into_keeps_every_iupac_code() {
    let wire = Wire {
        tid: -1,
        pos: -1,
        flag: 4,
        qname: b"r".to_vec(),
        cigar: Vec::new(),
        seq_len: 16,
        aux: Vec::new(),
    };
    let mut payload = wire.payload();
    let seq_at = 32 + 2;
    payload[seq_at..seq_at + 8].copy_from_slice(&[0x01, 0x23, 0x45, 0x67, 0x89, 0xAB, 0xCD, 0xEF]);
    let record = stream_records(&bam_of_payloads(&[payload])).remove(0);
    let mut seq = b">".to_vec();
    record.decode_seq_into(&mut seq);
    assert_eq!(seq, b">=ACMGRSVTWYHKDBN");
}
