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
    let stream = BgzfStream::builder(file).threads(threads).build().unwrap();
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
        assert_eq!(ours.cigar_bytes(), theirs.cigar().as_ref(), "{what}: cigar");
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
/// Streaming an htslib-written, coordinate-sorted BAM and writing every record
/// back with `write_raw` reproduces each record's bytes exactly, and the BAI
/// the writer builds alongside is byte for byte the one `samtools index`
/// builds for the output.
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
