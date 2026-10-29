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
    // Line by line: Debian's build flags in the output are not UTF-8.
    let (major, minor) = out
        .stdout
        .split(|&b| b == b'\n')
        .filter_map(|l| std::str::from_utf8(l).ok())
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
    assert!(matches!(
        record.set_flags(seqair::bam::flags::BamFlags::empty()),
        Err(seqair::bam::RawRecordEditError::Flags(_))
    ));
    assert_eq!(record.flags(), flags);
    assert_eq!(record.as_payload(), unmapped.payload());
    let reverse_unmapped = seqair::bam::flags::BamFlags::from(4 | 16);
    record.set_flags(reverse_unmapped).unwrap();
    assert_eq!(record.flags(), reverse_unmapped);
}

// ── generated SAM, written to BAM by htslib and by noodles ──────────────

use rust_htslib::bam::{self as htsbam, Read as _};
use seqair::bam::BgzfStream;
use seqair::bam::writer::BamWriterBuilder;
use std::path::{Path, PathBuf};

/// One SAM line's fields, kept as text.
#[derive(Debug, Clone)]
struct SamLine {
    qname: String,
    flag: u16,
    /// Index into [`REFS`]; `None` is `*`.
    rname: Option<usize>,
    /// 1-based; 0 for none.
    pos: u32,
    mapq: u8,
    cigar: String,
    rnext: String,
    pnext: u32,
    tlen: i32,
    seq: String,
    qual: String,
    aux: Vec<String>,
}

impl SamLine {
    fn text(&self) -> String {
        let rname = self.rname.map_or("*", |i| REFS[i].0);
        let mut line = format!(
            "{}\t{}\t{rname}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.qname,
            self.flag,
            self.pos,
            self.mapq,
            self.cigar,
            self.rnext,
            self.pnext,
            self.tlen,
            self.seq,
            self.qual
        );
        for aux in &self.aux {
            line.push('\t');
            line.push_str(aux);
        }
        line
    }

    /// Coordinate order, unplaced records last.
    fn sort_key(&self) -> (usize, u32) {
        (self.rname.unwrap_or(usize::MAX), self.pos)
    }
}

fn int_in(tc: &TestCase, lo: i64, hi: i64) -> i64 {
    tc.draw(gs::integers::<i64>().min_value(lo).max_value(hi))
}

/// One aux field of any SAM type, `B` arrays of every subtype included.
#[hegel::composite]
fn sam_aux(tc: &TestCase, tag: String) -> String {
    match tc.draw(gs::integers::<u8>().max_value(6)) {
        0 => format!("{tag}:A:{}", char::from(int_in(tc, 33, 126) as u8)),
        // Every width the writers choose between, both signs.
        1 => format!("{tag}:i:{}", int_in(tc, i64::from(i32::MIN), i64::from(u32::MAX))),
        2 => format!("{tag}:i:{}", int_in(tc, -300, 70_000)),
        3 => format!("{tag}:f:{}", tc.draw(gs::floats::<f32>().min_value(-1e6).max_value(1e6))),
        4 => format!("{tag}:Z:{}", tc.draw(gs::from_regex("[ !-~]{0,12}").fullmatch(true))),
        5 => {
            let n = tc.draw(gs::integers::<usize>().min_value(1).max_value(6));
            let mut hex = String::new();
            for _ in 0..n {
                write!(hex, "{:02X}", int_in(tc, 0, 255)).unwrap();
            }
            format!("{tag}:H:{hex}")
        }
        _ => {
            let (sub, lo, hi) = tc.draw(gs::sampled_from(vec![
                ('c', -128, 127),
                ('C', 0, 255),
                ('s', -32_768, 32_767),
                ('S', 0, 65_535),
                ('i', i64::from(i32::MIN), i64::from(i32::MAX)),
                ('I', 0, i64::from(u32::MAX)),
                ('f', -1000, 1000),
            ]));
            let n = tc.draw(gs::integers::<usize>().min_value(1).max_value(5));
            let mut out = format!("{tag}:B:{sub}");
            for _ in 0..n {
                write!(out, ",{}", int_in(tc, lo, hi)).unwrap();
            }
            out
        }
    }
}

/// A CIGAR valid by construction, and the query and reference lengths it
/// consumes.
#[hegel::composite]
fn sam_cigar(tc: &TestCase) -> (String, u32, u32) {
    let len = |tc: &TestCase, hi: i64| int_in(tc, 1, hi) as u32;
    let mut ops: Vec<(u32, char)> = Vec::new();
    if tc.draw(gs::booleans()) {
        ops.push((len(tc, 4), 'H'));
    }
    if tc.draw(gs::booleans()) {
        ops.push((len(tc, 8), 'S'));
    }
    let m = |tc: &TestCase| (len(tc, 30), tc.draw(gs::sampled_from(vec!['M', '=', 'X'])));
    ops.push(m(tc));
    for _ in 0..tc.draw(gs::integers::<usize>().max_value(2)) {
        let op = tc.draw(gs::sampled_from(vec!['I', 'D', 'N', 'P']));
        ops.push((if op == 'N' { len(tc, 300) } else { len(tc, 5) }, op));
        ops.push(m(tc));
    }
    if tc.draw(gs::booleans()) {
        ops.push((len(tc, 8), 'S'));
    }
    if tc.draw(gs::booleans()) {
        ops.push((len(tc, 4), 'H'));
    }
    let qlen = ops.iter().filter(|(_, o)| "MIS=X".contains(*o)).map(|(l, _)| l).sum();
    let rlen = ops.iter().filter(|(_, o)| "MDN=X".contains(*o)).map(|(l, _)| l).sum();
    let mut text = String::new();
    for (l, o) in &ops {
        write!(text, "{l}{o}").unwrap();
    }
    (text, qlen, rlen)
}

/// Every printable quality character, Phred 0 to 93.
static QUAL_CHARS: std::sync::LazyLock<String> = std::sync::LazyLock::new(|| ('!'..='~').collect());

/// IUPAC bases, `=` and lower case among them.
fn seq_of(tc: &TestCase, len: u32) -> String {
    tc.draw(
        gs::text().alphabet("=ACMGRSVTWYHKDBNacgtn").min_size(len as usize).max_size(len as usize),
    )
}

#[hegel::composite]
fn sam_line(tc: &TestCase, index: usize) -> SamLine {
    // 0: mapped; 1: placed but unmapped; 2: unplaced.
    let kind = tc.draw(gs::integers::<u8>().max_value(2));
    let mut flag = tc.draw(gs::integers::<u16>().max_value(0xFFF)) & !0x4;
    let (rname, cigar, qlen, rlen) = match kind {
        0 => {
            let (cigar, qlen, rlen) = tc.draw(sam_cigar());
            (Some(tc.draw(gs::integers::<usize>().max_value(1))), cigar, qlen, rlen)
        }
        k => {
            flag |= 0x4;
            let rname = (k == 1).then(|| tc.draw(gs::integers::<usize>().max_value(1)));
            (rname, "*".to_owned(), int_in(tc, 0, 40) as u32, 1)
        }
    };
    let pos = match rname {
        Some(r) => int_in(tc, 1, i64::from(REFS[r].1 - rlen)) as u32,
        None => 0,
    };
    let (seq, qual) = match tc.draw(gs::integers::<u8>().max_value(4)) {
        0 => ("*".to_owned(), "*".to_owned()),
        1 if qlen > 0 => (seq_of(tc, qlen), "*".to_owned()),
        _ if qlen > 0 => {
            let qual = tc.draw(
                gs::text().alphabet(&QUAL_CHARS).min_size(qlen as usize).max_size(qlen as usize),
            );
            (seq_of(tc, qlen), qual)
        }
        _ => ("*".to_owned(), "*".to_owned()),
    };
    let (rnext, pnext) = match tc.draw(gs::integers::<u8>().max_value(2)) {
        0 => ("*".to_owned(), 0),
        1 if rname.is_some() => ("=".to_owned(), int_in(tc, 1, 40_000) as u32),
        _ => {
            let r = tc.draw(gs::integers::<usize>().max_value(1));
            (REFS[r].0.to_owned(), int_in(tc, 1, 40_000) as u32)
        }
    };
    let n_aux = tc.draw(gs::integers::<usize>().max_value(4));
    let aux = (0..n_aux).map(|i| tc.draw(sam_aux(format!("X{i}")))).collect();
    SamLine {
        qname: format!(
            "{index}{}",
            tc.draw(gs::from_regex("[A-Za-z0-9_:.#-]{0,20}").fullmatch(true))
        ),
        flag,
        rname,
        pos,
        mapq: int_in(tc, 0, 255) as u8,
        cigar,
        rnext,
        pnext,
        tlen: int_in(tc, -1000, 1000) as i32,
        seq,
        qual,
        aux,
    }
}

#[hegel::composite]
fn sam_lines(tc: &TestCase) -> Vec<SamLine> {
    let n = tc.draw(gs::integers::<usize>().max_value(25));
    let mut lines: Vec<SamLine> = (0..n).map(|i| tc.draw(sam_line(i))).collect();
    lines.sort_by_key(SamLine::sort_key);
    lines
}

fn sam_text(lines: &[SamLine]) -> String {
    let mut text = header_text().replace("SO:unknown", "SO:coordinate");
    for line in lines {
        text.push_str(&line.text());
        text.push('\n');
    }
    text
}

/// The SAM written to BAM by htslib (through rust-htslib).
fn htslib_bam(dir: &Path, sam: &Path) -> PathBuf {
    let out = dir.join("htslib.bam");
    let mut reader = htsbam::Reader::from_path(sam).unwrap();
    let header = htsbam::Header::from_template(reader.header());
    let mut writer = htsbam::Writer::from_path(&out, &header, htsbam::Format::Bam).unwrap();
    for record in reader.records() {
        writer.write(&record.unwrap()).unwrap();
    }
    drop(writer);
    out
}

/// The SAM written to BAM by noodles.
fn noodles_bam(dir: &Path, sam: &Path) -> PathBuf {
    use noodles::sam::alignment::io::Write as _;
    let out = dir.join("noodles.bam");
    let mut reader =
        noodles::sam::io::Reader::new(std::io::BufReader::new(std::fs::File::open(sam).unwrap()));
    let header = reader.read_header().unwrap();
    let mut writer = noodles::bam::io::Writer::new(std::fs::File::create(&out).unwrap());
    writer.write_header(&header).unwrap();
    for record in reader.records() {
        writer.write_alignment_record(&header, &record.unwrap()).unwrap();
    }
    writer.try_finish().unwrap();
    out
}

/// Every record of the BAM at `path`, streamed by seqair.
fn seqair_stream(path: &Path, threads: usize) -> Vec<RawBamRecord> {
    let file = std::fs::File::open(path).unwrap();
    let builder = BgzfStream::builder(file);
    let stream = match u8::try_from(threads).ok().and_then(std::num::NonZeroU8::new) {
        Some(n) => builder.threads(n),
        None => builder,
    }
    .build();
    let mut reader = BamStreamReader::from_bgzf(stream).unwrap();
    reader.records().collect::<Result<_, _>>().unwrap()
}

/// What htslib holds for a record, as the BAM payload it was read from:
/// htslib pads the read name with NULs to a multiple of four in memory, so
/// that padding is taken out again.
fn htslib_payload(record: &htsbam::Record) -> Vec<u8> {
    let inner = record.inner();
    let core = &inner.core;
    // SAFETY: `data` points at `l_data` initialised bytes owned by `record`,
    // which outlives the slice.
    let data = unsafe { std::slice::from_raw_parts(inner.data, inner.l_data as usize) };
    let l_qname = usize::from(core.l_qname);
    let extranul = usize::from(core.l_extranul);
    let mut out = Vec::new();
    out.extend_from_slice(&core.tid.to_le_bytes());
    out.extend_from_slice(&(core.pos as i32).to_le_bytes());
    out.push((l_qname - extranul) as u8);
    out.push(core.qual);
    out.extend_from_slice(&core.bin.to_le_bytes());
    out.extend_from_slice(&(core.n_cigar as u16).to_le_bytes());
    out.extend_from_slice(&core.flag.to_le_bytes());
    out.extend_from_slice(&core.l_qseq.to_le_bytes());
    out.extend_from_slice(&core.mtid.to_le_bytes());
    out.extend_from_slice(&(core.mpos as i32).to_le_bytes());
    out.extend_from_slice(&(core.isize_ as i32).to_le_bytes());
    out.extend_from_slice(&data[..l_qname - extranul]);
    out.extend_from_slice(&data[l_qname..]);
    out
}

/// Each seqair record's accessors against htslib's decoding of the same
/// record, and its payload against htslib's bytes.
fn assert_agrees_with_htslib(path: &Path, ours: &[RawBamRecord]) {
    let mut reader = htsbam::Reader::from_path(path).unwrap();
    let theirs: Vec<htsbam::Record> = reader.records().map(Result::unwrap).collect();
    assert_eq!(ours.len(), theirs.len(), "record count");
    for (i, (ours, theirs)) in ours.iter().zip(&theirs).enumerate() {
        let what = format!("record {i} ({})", theirs.qname().escape_ascii());
        assert!(ours.as_payload() == htslib_payload(theirs), "{what}: payload");
        assert_eq!(ours.qname(), theirs.qname(), "{what}: qname");
        assert_eq!(ours.flags().raw(), theirs.flags(), "{what}: flags");
        assert_eq!(ours.tid().map_or(-1, |t| t.as_u32() as i32), theirs.tid(), "{what}: tid");
        assert_eq!(ours.pos().map_or(-1, |p| i64::from(p.as_i32())), theirs.pos(), "{what}: pos");
        assert_eq!(ours.mapq().as_byte(), theirs.mapq(), "{what}: mapq");
        assert_eq!(ours.bin(), theirs.bin(), "{what}: bin");
        let cigar: Vec<u32> = ours.cigar().map(|op| op.to_bam_u32()).collect();
        assert_eq!(cigar, theirs.raw_cigar(), "{what}: cigar");
        assert_eq!(ours.cigar_len(), theirs.raw_cigar().len(), "{what}: cigar_len");
        assert_eq!(ours.seq_len(), theirs.seq_len(), "{what}: seq_len");
        let mut seq = Vec::new();
        ours.decode_seq_into(&mut seq);
        assert_eq!(seq, theirs.seq().as_bytes(), "{what}: seq");
        let qual: Vec<u8> = ours
            .qual()
            .map(|q| seqair_types::BaseQuality::slice_to_bytes(q).to_vec())
            .unwrap_or_else(|| vec![0xFF; ours.seq_len()]);
        assert_eq!(qual, theirs.qual(), "{what}: qual");
        let mate = ours.mate_fields();
        assert_eq!(mate.tid.map_or(-1, |t| t.as_u32() as i32), theirs.mtid(), "{what}: mtid");
        assert_eq!(mate.pos.map_or(-1, |p| i64::from(p.as_i32())), theirs.mpos(), "{what}: mpos");
        assert_eq!(
            mate.template_len.map_or(0, |t| std::num::NonZeroI64::from(t).get()),
            theirs.insert_size(),
            "{what}: tlen"
        );
        assert!(ours.aux().iter_checked().all(|t| t.is_ok()), "{what}: aux");
        if ours.tid().is_some() && !ours.flags().is_unmapped() {
            let end = ours.end_pos().map(|p| i64::from(p.as_i32()) + 1);
            assert_eq!(end, Some(theirs.cigar().end_pos()), "{what}: end_pos");
        }
    }
}

/// Each seqair record against noodles' decoding of the same record.
fn assert_agrees_with_noodles(path: &Path, ours: &[RawBamRecord]) {
    let mut reader = noodles::bam::io::Reader::new(std::fs::File::open(path).unwrap());
    reader.read_header().unwrap();
    let theirs: Vec<noodles::bam::Record> = reader.records().map(Result::unwrap).collect();
    assert_eq!(ours.len(), theirs.len(), "record count");
    let id = |r: Option<std::io::Result<usize>>| r.map(Result::unwrap);
    let start = |r: Option<std::io::Result<noodles::core::Position>>| {
        r.map(|p| usize::from(p.unwrap()) - 1)
    };
    for (i, (ours, theirs)) in ours.iter().zip(&theirs).enumerate() {
        let what = format!("record {i}");
        let name: &[u8] = theirs.name().map_or(b"*", |n| n.as_ref());
        assert_eq!(ours.qname(), name, "{what}: name");
        assert_eq!(ours.flags().raw(), theirs.flags().bits(), "{what}: flags");
        assert_eq!(ours.tid().map(|t| t.as_usize()), id(theirs.reference_sequence_id()));
        assert_eq!(ours.pos().map(|p| p.as_u64() as usize), start(theirs.alignment_start()));
        assert_eq!(
            Some(ours.mapq().as_byte()).filter(|&q| q != 255),
            theirs.mapping_quality().map(u8::from),
            "{what}: mapq"
        );
        assert_eq!(ours.stored_cigar_bytes(), theirs.cigar().as_ref(), "{what}: cigar");
        assert_eq!(ours.packed_seq(), theirs.sequence().as_ref(), "{what}: seq");
        assert_eq!(ours.seq_len(), theirs.sequence().len(), "{what}: seq_len");
        // noodles shows missing qualities as none at all.
        let qual = ours.qual().map_or(&[][..], seqair_types::BaseQuality::slice_to_bytes);
        assert_eq!(qual, theirs.quality_scores().as_ref(), "{what}: qual");
        assert_eq!(ours.aux().as_bytes(), theirs.data().as_ref(), "{what}: aux");
        let mate = ours.mate_fields();
        assert_eq!(mate.tid.map(|t| t.as_usize()), id(theirs.mate_reference_sequence_id()));
        assert_eq!(mate.pos.map(|p| p.as_u64() as usize), start(theirs.mate_alignment_start()));
        assert_eq!(mate.template_len.map_or(0, |t| t.get()), theirs.template_length());
    }
}

// r[verify bam.stream]
// r[verify bam.stream.framing]
// r[verify bam.raw_record]
// r[verify bam.raw_record.validation]
// r[verify bam.raw_record.bin_verbatim]
// r[verify bam.record.unplaced_pos]
/// Random SAM — mapped, placed-unmapped and unplaced records, IUPAC bases,
/// every aux type, missing SEQ or QUAL — written to BAM by htslib and by
/// noodles, then streamed by seqair at a random thread count: every record
/// is there, in order, and agrees field by field with htslib's and noodles'
/// sequential readers on both files, its payload byte for byte with htslib's.
#[hegel::test(test_cases = pinned::cases(60))]
fn the_stream_agrees_with_htslib_and_noodles(tc: TestCase) {
    let lines = tc.draw(sam_lines());
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("in.sam");
    std::fs::write(&sam, sam_text(&lines)).unwrap();
    let threads = tc.draw(gs::integers::<usize>().max_value(3));
    for bam in [htslib_bam(dir.path(), &sam), noodles_bam(dir.path(), &sam)] {
        let ours = seqair_stream(&bam, threads);
        assert_eq!(ours.len(), lines.len(), "{}", bam.display());
        for (record, line) in ours.iter().zip(&lines) {
            assert_eq!(record.qname(), line.qname.as_bytes());
        }
        assert_agrees_with_htslib(&bam, &ours);
        assert_agrees_with_noodles(&bam, &ours);
    }
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

// ── write_raw ───────────────────────────────────────────────────────────

/// The header and every record of the BAM at `path`.
fn stream_file(path: &Path) -> (seqair::bam::BamHeader, Vec<RawBamRecord>) {
    let mut reader = BamStreamReader::new(std::fs::File::open(path).unwrap()).unwrap();
    let records = reader.records().collect::<Result<_, _>>().unwrap();
    (reader.into_header(), records)
}

fn run(cmd: &mut std::process::Command) {
    let out = cmd.output().unwrap();
    assert!(out.status.success(), "{cmd:?}: {}", String::from_utf8(out.stderr).unwrap());
}

// r[verify bam_writer.write_raw]
// r[verify bam.raw_record]
// r[verify index_builder.amend_last]
// r[verify index_builder.compress_bins]
// r[verify index_builder.no_coor]
/// Streaming an htslib-written, coordinate-sorted BAM and writing every record
/// back with `write_raw` reproduces each record's bytes exactly, and the BAI
/// the writer builds alongside holds the same bins, chunks and linear index as
/// the one `samtools index` builds for the output (bins compared by id).
#[hegel::test(test_cases = pinned::cases(40))]
fn write_raw_round_trips_and_indexes_like_samtools(tc: TestCase) {
    let mut lines = tc.draw(sam_lines());
    // Sometimes enough records for many BGZF blocks, so chunks cross block
    // boundaries and htslib's bin compression has work to do.
    let copies = tc.draw(gs::sampled_from(vec![1usize, 1, 50, 800]));
    lines = lines
        .iter()
        .flat_map(|l| {
            (0..copies).map(|i| SamLine { qname: format!("{}_{i}", l.qname), ..l.clone() })
        })
        .collect();
    lines.sort_by_key(SamLine::sort_key);
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("in.sam");
    std::fs::write(&sam, sam_text(&lines)).unwrap();
    let (header, records) = stream_file(&htslib_bam(dir.path(), &sam));

    let out = dir.path().join("out.bam");
    let threads = tc.draw(gs::integers::<usize>().max_value(2));
    let mut writer = BamWriterBuilder::to_path(&out, &header)
        .write_index(true)
        .compression_threads(threads)
        .build()
        .unwrap();
    for record in &records {
        writer.write_raw(record).unwrap();
    }
    let (_, index) = writer.finish().unwrap();
    let ours_bai = dir.path().join("ours.bai");
    index
        .unwrap()
        .write_bai(std::fs::File::create(&ours_bai).unwrap(), header.target_count())
        .unwrap();

    let (_, again) = stream_file(&out);
    assert_eq!(again.len(), records.len());
    for (i, (a, b)) in again.iter().zip(&records).enumerate() {
        assert!(a.as_payload() == b.as_payload(), "record {i} changed");
    }
    let samtools_bai = dir.path().join("samtools.bai");
    run(std::process::Command::new("samtools").arg("index").arg(&out).arg("-o").arg(&samtools_bai));
    assert_eq!(
        parse_bai(&std::fs::read(&ours_bai).unwrap()),
        parse_bai(&std::fs::read(&samtools_bai).unwrap()),
        "BAI differs from samtools index's"
    );
}

/// One reference of a BAI: bins (by id) with their chunks, and the linear
/// index.
type BaiRef = (std::collections::BTreeMap<u32, Vec<(u64, u64)>>, Vec<u64>);

/// A BAI's references and its `n_no_coor`, with bins keyed by id: seqair
/// writes them in ascending order, htslib in its hash table's order, and the
/// format does not care.
fn parse_bai(bytes: &[u8]) -> (Vec<BaiRef>, Option<u64>) {
    fn take<'a>(rest: &mut &'a [u8], n: usize) -> &'a [u8] {
        let (head, tail) = rest.split_at(n);
        *rest = tail;
        head
    }
    let mut rest = bytes;
    let mut take = |n: usize| take(&mut rest, n);
    assert_eq!(take(4), b"BAI\x01");
    let u32_ = |b: &[u8]| u32::from_le_bytes(b.try_into().unwrap());
    let u64_ = |b: &[u8]| u64::from_le_bytes(b.try_into().unwrap());
    let n_ref = u32_(take(4));
    let refs = (0..n_ref)
        .map(|_| {
            let n_bin = u32_(take(4));
            let bins = (0..n_bin)
                .map(|_| {
                    let bin = u32_(take(4));
                    let n_chunk = u32_(take(4));
                    let chunks = (0..n_chunk).map(|_| (u64_(take(8)), u64_(take(8)))).collect();
                    (bin, chunks)
                })
                .collect();
            let n_intv = u32_(take(4));
            (bins, (0..n_intv).map(|_| u64_(take(8))).collect())
        })
        .collect();
    let n_no_coor = match rest.len() {
        0 => None,
        8 => Some(u64_(rest)),
        n => panic!("{n} trailing bytes"),
    };
    (refs, n_no_coor)
}

/// A writer refuses a record whose reference its header lacks, before
/// writing anything.
#[test]
fn write_raw_rejects_a_tid_the_header_lacks() {
    let wire = Wire {
        tid: 1,
        pos: 10,
        flag: 0,
        qname: b"r".to_vec(),
        cigar: vec![op(4, 0)],
        seq_len: 4,
        aux: Vec::new(),
    };
    let record = stream_records(&bam_of_payloads(&[wire.payload()])).remove(0);
    let one_ref = seqair::bam::BamHeader::from_sam_text("@SQ\tSN:chr1\tLN:100000\n").unwrap();
    let mut writer = BamWriterBuilder::to_writer(Vec::new(), &one_ref).build().unwrap();
    assert!(matches!(
        writer.write_raw(&record),
        Err(seqair::bam::BamWriteError::TidOutOfRange { tid: 1, n_targets: 1 })
    ));
}

// ── write_store_record's bin ────────────────────────────────────────────

use seqair::bam::RecordStore;
use seqair::bam::reader::IndexedBamReader;

/// A mapped read of `len` reference bases ending `delta` bases after a
/// 16 kbp window boundary (`delta` = 0: its last base is the window's first).
#[hegel::composite]
fn read_at_a_window_edge(tc: &TestCase, index: usize) -> SamLine {
    let window = int_in(tc, 1, 5) as u32 * 16_384;
    let len = int_in(tc, 1, 40) as u32;
    let delta = int_in(tc, -2, 2);
    // 0-based last base = window + delta; SAM POS is 1-based.
    let pos = (i64::from(window) + delta - i64::from(len) + 2) as u32;
    let unmapped = tc.draw(gs::integers::<u8>().max_value(5)) == 0;
    SamLine {
        qname: format!("r{index}"),
        flag: if unmapped { 4 } else { 0 },
        rname: Some(0),
        pos,
        mapq: 60,
        cigar: if unmapped { "*".to_owned() } else { format!("{len}M") },
        rnext: "*".to_owned(),
        pnext: 0,
        tlen: 0,
        seq: "A".repeat(len as usize),
        qual: "I".repeat(len as usize),
        aux: Vec::new(),
    }
}

// r[verify bam_writer.store_record_bin]
// r[verify bam_writer.write_store_record]
/// `write_store_record` gives a read ending at a 16 kbp window edge the bin
/// htslib computes, and indexes it as `samtools index` does. It used to bin
/// with the inclusive end, so a read whose last base is a window's first
/// landed one level too deep.
#[hegel::test(test_cases = pinned::cases(40))]
fn write_store_record_bins_like_htslib(tc: TestCase) {
    let n = tc.draw(gs::integers::<usize>().min_value(1).max_value(12));
    let mut lines: Vec<SamLine> = (0..n).map(|i| tc.draw(read_at_a_window_edge(i))).collect();
    lines.sort_by_key(SamLine::sort_key);
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("in.sam");
    std::fs::write(&sam, sam_text(&lines)).unwrap();
    let input = htslib_bam(dir.path(), &sam);
    run(std::process::Command::new("samtools").arg("index").arg(&input));
    let (header, theirs) = stream_file(&input);

    let mut reader = IndexedBamReader::open(&input).unwrap();
    let tid = header.tid("chr1").unwrap();
    let mut store = RecordStore::new();
    let whole =
        (seqair::bam::Pos0::new(0).unwrap()..=seqair::bam::Pos0::new(99_999).unwrap()).into();
    reader.fetch_into(tid, whole, &mut store).unwrap();
    let out = dir.path().join("out.bam");
    let mut writer = BamWriterBuilder::to_path(&out, &header).write_index(true).build().unwrap();
    for idx in store.indices() {
        writer.write_store_record(&store, idx).unwrap();
    }
    let (_, index) = writer.finish().unwrap();
    let ours_bai = dir.path().join("ours.bai");
    index
        .unwrap()
        .write_bai(std::fs::File::create(&ours_bai).unwrap(), header.target_count())
        .unwrap();

    let (_, ours) = stream_file(&out);
    // The store keeps the mapped reads; htslib's file has the unmapped too.
    let mapped: Vec<&RawBamRecord> = theirs.iter().filter(|r| !r.flags().is_unmapped()).collect();
    assert_eq!(ours.len(), store.len());
    let by_name = |r: &&RawBamRecord| r.qname().to_vec();
    let mut mapped_sorted = mapped.clone();
    mapped_sorted.sort_by_key(by_name);
    let mut ours_sorted: Vec<&RawBamRecord> =
        ours.iter().filter(|r| !r.flags().is_unmapped()).collect();
    ours_sorted.sort_by_key(by_name);
    assert_eq!(ours_sorted.len(), mapped_sorted.len());
    for (a, b) in ours_sorted.iter().zip(&mapped_sorted) {
        assert_eq!(a.bin(), b.bin(), "{}: bin", a.qname().escape_ascii());
    }
    let samtools_bai = dir.path().join("samtools.bai");
    run(std::process::Command::new("samtools").arg("index").arg(&out).arg("-o").arg(&samtools_bai));
    assert_eq!(
        parse_bai(&std::fs::read(&ours_bai).unwrap()),
        parse_bai(&std::fs::read(&samtools_bai).unwrap()),
        "BAI differs from samtools index's"
    );
}

#[test]
fn a_read_whose_last_base_opens_a_window_is_binned_one_level_up() {
    // 10M at 0-based 16375 covers 16375..=16384: two 16 kbp windows, so
    // htslib's bin is the level-4 bin 585, not the leaf 4681.
    let line = SamLine {
        qname: "r".to_owned(),
        flag: 0,
        rname: Some(0),
        pos: 16_376,
        mapq: 60,
        cigar: "10M".to_owned(),
        rnext: "*".to_owned(),
        pnext: 0,
        tlen: 0,
        seq: "ACGTACGTAC".to_owned(),
        qual: "IIIIIIIIII".to_owned(),
        aux: Vec::new(),
    };
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("in.sam");
    std::fs::write(&sam, sam_text(&[line])).unwrap();
    let input = htslib_bam(dir.path(), &sam);
    run(std::process::Command::new("samtools").arg("index").arg(&input));
    let (header, theirs) = stream_file(&input);
    assert_eq!(theirs[0].bin(), 585);

    let mut reader = IndexedBamReader::open(&input).unwrap();
    let mut store = RecordStore::new();
    let whole =
        (seqair::bam::Pos0::new(0).unwrap()..=seqair::bam::Pos0::new(99_999).unwrap()).into();
    reader.fetch_into(header.tid("chr1").unwrap(), whole, &mut store).unwrap();
    let mut writer = BamWriterBuilder::to_writer(Vec::new(), &header).build().unwrap();
    writer.write_store_record(&store, store.indices().next().unwrap()).unwrap();
    let (bytes, _) = writer.finish().unwrap();
    let ours = stream_records(&bytes);
    assert_eq!(ours[0].bin(), 585);
}

// ── framing ─────────────────────────────────────────────────────────────

/// A BAM of the header and `body` (raw bytes after it), compressed in one go.
fn bam_with_body(body: &[u8]) -> Vec<u8> {
    let mut writer = BgzfWriter::new(Vec::new());
    writer.write_all(&header_bytes()).unwrap();
    writer.write_all(body).unwrap();
    writer.finish().unwrap()
}

/// The records of `bam` up to the first error, and that error.
fn read_until_error(bam: &[u8]) -> (Vec<RawBamRecord>, Option<BamStreamError>) {
    let mut reader = BamStreamReader::new(Cursor::new(bam)).unwrap();
    let mut records = Vec::new();
    loop {
        let mut record = RawBamRecord::default();
        match reader.read_record(&mut record) {
            Ok(true) => records.push(record),
            Ok(false) => return (records, None),
            Err(e) => {
                let mut again = RawBamRecord::default();
                assert!(
                    matches!(reader.read_record(&mut again), Err(BamStreamError::Poisoned)),
                    "not poisoned after {e:?}"
                );
                return (records, Some(e));
            }
        }
    }
}

fn framed(payload: &[u8]) -> Vec<u8> {
    let mut out = (payload.len() as i32).to_le_bytes().to_vec();
    out.extend_from_slice(payload);
    out
}

fn simple_record(index: usize) -> Vec<u8> {
    Wire {
        tid: 0,
        pos: 100 + index as i32,
        flag: 0,
        qname: format!("r{index}").into_bytes(),
        cigar: vec![op(4, 0)],
        seq_len: 4,
        aux: b"XAAQ".to_vec(),
    }
    .payload()
}

// r[verify bam.stream.framing]
/// A `block_size` below 32 or above 2 MiB fails with `RecordSizeOutOfRange`
/// after the records before it, and poisons the stream.
#[test]
fn a_block_size_out_of_range_is_an_error() {
    for size in [i32::MIN, -1, 0, 31, 2 * 1024 * 1024 + 1, i32::MAX] {
        let mut body = framed(&simple_record(0));
        body.extend_from_slice(&size.to_le_bytes());
        body.extend_from_slice(&[0; 64]);
        let (records, error) = read_until_error(&bam_with_body(&body));
        assert_eq!(records.len(), 1, "size {size}");
        assert!(
            matches!(error, Some(BamStreamError::RecordSizeOutOfRange { record: 2, size: s }) if s == size),
            "size {size}: {error:?}"
        );
    }
}

// r[verify bam.stream.framing]
// r[verify bgzf.stream.eof_marker]
/// Cutting the uncompressed records anywhere (and compressing what is left
/// properly): a cut between records is a clean end after the records before
/// it; a cut inside a `block_size` or a payload is `TruncatedRecord`, naming
/// the record it cut, after every whole record before it.
#[test]
fn a_record_cut_anywhere_is_truncated_unless_cut_between_records() {
    let records: Vec<Vec<u8>> = (0..4).map(simple_record).collect();
    let body: Vec<u8> = records.iter().flat_map(|r| framed(r)).collect();
    let mut ends = vec![0];
    for r in &records {
        ends.push(ends.last().unwrap() + 4 + r.len());
    }
    for cut in 0..=body.len() {
        let (got, error) = read_until_error(&bam_with_body(&body[..cut]));
        let whole = ends.iter().filter(|&&e| e <= cut).count() - 1;
        assert_eq!(got.len(), whole, "cut {cut}");
        for (g, r) in got.iter().zip(&records) {
            assert_eq!(g.as_payload(), r.as_slice());
        }
        if ends.contains(&cut) {
            assert!(error.is_none(), "cut {cut}: {error:?}");
        } else {
            assert!(
                matches!(error, Some(BamStreamError::TruncatedRecord { record }) if record == whole as u64 + 1),
                "cut {cut}: {error:?}"
            );
        }
    }
}

/// Ways to make a well-framed record invalid.
#[derive(Debug, Clone, Copy)]
enum Defect {
    TidOutOfRange,
    MateTidOutOfRange,
    PosBelowMinusOne,
    MatePosBelowMinusOne,
    NameNotTerminated,
    EmptyName,
    CigarSeqMismatch,
    SeqLongerThanPayload,
}

impl Defect {
    fn apply(self, mut payload: Vec<u8>) -> Vec<u8> {
        let set_i32 =
            |p: &mut Vec<u8>, at: usize, v: i32| p[at..at + 4].copy_from_slice(&v.to_le_bytes());
        match self {
            Self::TidOutOfRange => set_i32(&mut payload, 0, 2),
            Self::MateTidOutOfRange => set_i32(&mut payload, 20, 7),
            Self::PosBelowMinusOne => set_i32(&mut payload, 4, -2),
            Self::MatePosBelowMinusOne => set_i32(&mut payload, 24, -9),
            Self::NameNotTerminated => {
                let l = usize::from(payload[8]);
                payload[32 + l - 1] = b'x';
            }
            Self::EmptyName => {
                let l = usize::from(payload[8]);
                payload.drain(32..32 + l);
                payload[8] = 0;
            }
            Self::CigarSeqMismatch => {
                let l = usize::from(payload[8]);
                payload[32 + l..32 + l + 4].copy_from_slice(&op(3, 0).to_le_bytes());
            }
            Self::SeqLongerThanPayload => set_i32(&mut payload, 16, 400),
        }
        payload
    }
}

// r[verify bam.stream.framing]
// r[verify bam.raw_record.validation]
/// A well-framed but invalid record is `InvalidRecord`, numbered, and the
/// stream reads on: every valid record after it still comes out, byte for
/// byte, in order.
#[hegel::test]
fn an_invalid_record_leaves_the_stream_usable(tc: TestCase) {
    let n = tc.draw(gs::integers::<usize>().min_value(1).max_value(12));
    let defects = [
        Defect::TidOutOfRange,
        Defect::MateTidOutOfRange,
        Defect::PosBelowMinusOne,
        Defect::MatePosBelowMinusOne,
        Defect::NameNotTerminated,
        Defect::EmptyName,
        Defect::CigarSeqMismatch,
        Defect::SeqLongerThanPayload,
    ];
    let plan: Vec<Option<Defect>> = (0..n)
        .map(|_| tc.draw(gs::booleans()).then(|| tc.draw(gs::sampled_from(defects.to_vec()))))
        .collect();
    let payloads: Vec<Vec<u8>> = plan
        .iter()
        .enumerate()
        .map(|(i, d)| d.map_or_else(|| simple_record(i), |d| d.apply(simple_record(i))))
        .collect();
    let body: Vec<u8> = payloads.iter().flat_map(|p| framed(p)).collect();
    let bam = bam_with_body(&body);
    let mut reader = BamStreamReader::new(Cursor::new(&bam)).unwrap();
    let mut record = RawBamRecord::default();
    for (i, (defect, payload)) in plan.iter().zip(&payloads).enumerate() {
        match (reader.read_record(&mut record), defect) {
            (Ok(true), None) => assert_eq!(record.as_payload(), payload.as_slice()),
            (Err(BamStreamError::InvalidRecord { record: number, .. }), Some(_)) => {
                assert_eq!(number, i as u64 + 1);
                // A failed refill leaves the default record, not half of one.
                assert_eq!(record, RawBamRecord::default());
            }
            (got, _) => panic!("record {i} ({defect:?}): {got:?}"),
        }
    }
    assert!(!reader.read_record(&mut record).unwrap());
}

// ── region queries over empty blocks ────────────────────────────────────

/// The 28-byte BGZF EOF marker block, which is an empty block.
const EMPTY_BLOCK: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0, 0x42, 0x43, 0x02, 0, 0x1b, 0, 0x03, 0, 0,
    0, 0, 0, 0, 0, 0, 0,
];

/// One BGZF block holding `data`.
fn one_block(data: &[u8]) -> Vec<u8> {
    let mut writer = BgzfWriter::new(Vec::new());
    writer.write_all(data).unwrap();
    let mut bytes = writer.finish().unwrap();
    bytes.truncate(bytes.len() - EMPTY_BLOCK.len());
    bytes
}

/// Mapped reads' names in `samtools view <bam> <region>`.
fn samtools_region(bam: &Path, region: &str) -> Vec<Vec<u8>> {
    let out =
        std::process::Command::new("samtools").arg("view").arg(bam).arg(region).output().unwrap();
    assert!(out.status.success());
    let mut names: Vec<Vec<u8>> = out
        .stdout
        .split(|&b| b == b'\n')
        .filter(|l| !l.is_empty())
        .filter_map(|l| {
            let fields: Vec<&[u8]> = l.split(|&b| b == b'\t').collect();
            let flag: u16 = std::str::from_utf8(fields[1]).unwrap().parse().unwrap();
            (flag & 4 == 0).then(|| fields[0].to_vec())
        })
        .collect();
    names.sort();
    names
}

// r[verify region_buf.empty_block]
// r[verify bgzf.empty_block]
/// Empty blocks between the records of a BAM — inside the ranges a region
/// query reads — do not end the query early: seqair's indexed reader returns
/// the reads samtools returns for every region, on an index samtools built.
#[hegel::test(test_cases = pinned::cases(30))]
fn region_queries_read_past_empty_blocks(tc: TestCase) {
    let mut lines = tc.draw(sam_lines());
    lines.retain(|l| l.rname == Some(0) && l.flag & 4 == 0);
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("in.sam");
    std::fs::write(&sam, sam_text(&lines)).unwrap();
    let (header, records) = stream_file(&htslib_bam(dir.path(), &sam));

    // Records in blocks of random size, empty blocks between them.
    let mut bam = one_block(&header_bytes());
    let mut at = 0;
    while at < records.len() {
        let n = tc.draw(gs::integers::<usize>().min_value(1).max_value(4)).min(records.len() - at);
        let body: Vec<u8> =
            records[at..at + n].iter().flat_map(|r| framed(r.as_payload())).collect();
        for _ in 0..tc.draw(gs::integers::<usize>().max_value(2)) {
            bam.extend_from_slice(&EMPTY_BLOCK);
        }
        bam.extend_from_slice(&one_block(&body));
        at += n;
    }
    bam.extend_from_slice(&EMPTY_BLOCK);
    let path = dir.path().join("empty_blocks.bam");
    std::fs::write(&path, &bam).unwrap();
    run(std::process::Command::new("samtools").arg("index").arg(&path));

    let mut reader = IndexedBamReader::open(&path).unwrap();
    let tid = header.tid("chr1").unwrap();
    for _ in 0..4 {
        let start = int_in(&tc, 0, 99_000) as u32;
        let end = start + int_in(&tc, 0, 5_000) as u32;
        let mut store = RecordStore::new();
        let span =
            (seqair::bam::Pos0::new(start).unwrap()..=seqair::bam::Pos0::new(end).unwrap()).into();
        reader.fetch_into(tid, span, &mut store).unwrap();
        let mut ours: Vec<Vec<u8>> =
            store.indices().map(|i| store.record(i).unwrap().qname().to_vec()).collect();
        ours.sort();
        let region = format!("chr1:{}-{}", start + 1, end + 1);
        assert_eq!(ours, samtools_region(&path, &region), "{region}");
    }
}

// r[verify region_buf.truncated]
/// A BAM cut short anywhere, queried through the index of the whole file,
/// either returns exactly the reads samtools returns from the whole file, or
/// fails — never a silently shorter answer.
#[hegel::test(test_cases = pinned::cases(30))]
fn region_queries_on_a_truncated_bam_are_whole_or_errors(tc: TestCase) {
    let mut lines = tc.draw(sam_lines());
    lines.retain(|l| l.rname == Some(0) && l.flag & 4 == 0);
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("in.sam");
    std::fs::write(&sam, sam_text(&lines)).unwrap();
    let (header, records) = stream_file(&htslib_bam(dir.path(), &sam));

    // Small blocks, so a cut can fall inside, between or after the ones a
    // query needs.
    let mut bam = one_block(&header_bytes());
    let mut at = 0;
    while at < records.len() {
        let n = tc.draw(gs::integers::<usize>().min_value(1).max_value(3)).min(records.len() - at);
        let body: Vec<u8> =
            records[at..at + n].iter().flat_map(|r| framed(r.as_payload())).collect();
        bam.extend_from_slice(&one_block(&body));
        at += n;
    }
    bam.extend_from_slice(&EMPTY_BLOCK);
    let whole = dir.path().join("whole.bam");
    std::fs::write(&whole, &bam).unwrap();
    run(std::process::Command::new("samtools").arg("index").arg(&whole));

    let header_len = one_block(&header_bytes()).len();
    let cut = int_in(&tc, header_len as i64, bam.len() as i64) as usize;
    let path = dir.path().join("cut.bam");
    std::fs::write(&path, &bam[..cut]).unwrap();
    std::fs::copy(dir.path().join("whole.bam.bai"), dir.path().join("cut.bam.bai")).unwrap();

    let mut reader = IndexedBamReader::open(&path).unwrap();
    let tid = header.tid("chr1").unwrap();
    for _ in 0..4 {
        let start = int_in(&tc, 0, 99_000) as u32;
        let end = start + int_in(&tc, 0, 20_000) as u32;
        let mut store = RecordStore::new();
        let span =
            (seqair::bam::Pos0::new(start).unwrap()..=seqair::bam::Pos0::new(end).unwrap()).into();
        let region = format!("chr1:{}-{}", start + 1, end + 1);
        if reader.fetch_into(tid, span, &mut store).is_err() {
            // A failed query poisons nothing a fresh reader cannot redo.
            reader = IndexedBamReader::open(&path).unwrap();
            continue;
        }
        let mut ours: Vec<Vec<u8>> =
            store.indices().map(|i| store.record(i).unwrap().qname().to_vec()).collect();
        ours.sort();
        assert_eq!(ours, samtools_region(&whole, &region), "{region}, cut at {cut}");
    }
}

// r[verify bam.stream.payloads]
/// `read_payload_into` frames exactly what `read_record` reads — validating
/// its payloads later gives the same records, invalid ones included — and a
/// truncated stream leaves the caller's buffer as it was, then is poisoned.
#[hegel::test]
fn payloads_validate_later_into_the_same_records(tc: TestCase) {
    let n = tc.draw(gs::integers::<usize>().min_value(1).max_value(12));
    let payloads: Vec<Vec<u8>> = (0..n)
        .map(|i| {
            let record = simple_record(i);
            if tc.draw(gs::weighted_booleans(0.2)) {
                Defect::EmptyName.apply(record)
            } else {
                record
            }
        })
        .collect();
    let body: Vec<u8> = payloads.iter().flat_map(|p| framed(p)).collect();
    let bam = bam_with_body(&body);
    let header = BamStreamReader::new(Cursor::new(&bam)).unwrap().into_header();

    let mut by_record = BamStreamReader::new(Cursor::new(&bam)).unwrap();
    let mut by_payload = BamStreamReader::new(Cursor::new(&bam)).unwrap();
    let mut record = RawBamRecord::default();
    for payload in &payloads {
        let mut buf = b"prefix".to_vec();
        assert!(by_payload.read_payload_into(&mut buf).unwrap());
        assert_eq!(&buf[..6], b"prefix");
        assert_eq!(&buf[6..], payload.as_slice());
        let later = RawBamRecord::from_payload(buf[6..].to_vec(), &header);
        match by_record.read_record(&mut record) {
            Ok(true) => assert_eq!(later.unwrap(), record),
            Err(BamStreamError::InvalidRecord { .. }) => assert!(later.is_err()),
            other => panic!("{other:?}"),
        }
    }
    let mut buf = vec![1, 2, 3];
    assert!(!by_payload.read_payload_into(&mut buf).unwrap());
    assert_eq!(buf, [1, 2, 3]);

    // Cut inside the last record's payload.
    let mut cut = header_bytes();
    cut.extend_from_slice(&body[..body.len() - 1]);
    let mut writer = BgzfWriter::new(Vec::new());
    writer.write_all(&cut).unwrap();
    let truncated = writer.finish().unwrap();
    let mut reader = BamStreamReader::new(Cursor::new(&truncated)).unwrap();
    let mut buf = Vec::new();
    for _ in 1..payloads.len() {
        assert!(reader.read_payload_into(&mut buf).unwrap());
    }
    let before = buf.clone();
    assert!(matches!(
        reader.read_payload_into(&mut buf),
        Err(BamStreamError::TruncatedRecord { .. })
    ));
    assert_eq!(buf, before);
    assert!(matches!(reader.read_payload_into(&mut buf), Err(BamStreamError::Poisoned)));
}

// r[verify bam.raw_record.aux_lazy]
/// A record whose aux block is malformed streams and validates, as htslib
/// copies it; the malformation surfaces where the tags are interpreted.
#[test]
fn a_malformed_aux_block_is_read_and_reported_when_interpreted() {
    let mut payload = simple_record(0);
    payload.extend_from_slice(b"XXQ\x01"); // unknown type code `Q`
    let bam = bam_with_body(&framed(&payload));
    let mut reader = BamStreamReader::new(Cursor::new(&bam)).unwrap();
    let mut record = RawBamRecord::default();
    assert!(reader.read_record(&mut record).unwrap());
    assert_eq!(record.as_payload(), payload.as_slice());
    let first = record.aux().iter_checked().find_map(Result::err);
    assert!(
        matches!(first, Some(seqair::bam::aux::AuxLayoutError::UnknownType { code: b'Q', .. })),
        "{first:?}"
    );
}

// r[verify bam.raw_record.validation]
/// `from_payload` refuses a payload past the record size limit, so no
/// `RawBamRecord` exists that `write_raw` would then refuse.
#[test]
fn from_payload_enforces_the_record_size_limit() {
    let header = BamStreamReader::new(Cursor::new(&bam_with_body(&[]))).unwrap().into_header();
    let mut payload = simple_record(0);
    payload.resize(2 * 1024 * 1024 + 1, 0);
    assert!(matches!(
        RawBamRecord::from_payload(payload, &header),
        Err(seqair::bam::DecodeError::RecordTooLarge { .. })
    ));
}
