//! The whole read-then-write path, on records nobody wrote by hand.
//!
//! A generated SAM goes in one end: seqair's SAM reader parses it into a
//! `RecordStore`, seqair's `BamWriter` writes that store out as BAM, and
//! samtools turns the BAM back into SAM text. Every field of every record must
//! survive the trip.
//!
//! That is the SAM text parser, the record store's slabs, `write_store_record`
//! and the BAM binary layout in one comparison — and samtools decoding the
//! result means seqair cannot agree with itself and be wrong.
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
    reason = "test code with known small values"
)]

use hegel::prelude::*;
use seqair::bam::writer::BamWriterBuilder;
use seqair::bam::{Pos0, RecordStore};
use seqair::sam::reader::IndexedSamReader;
use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::Path;
use std::process::{Command, Stdio};

const CONTIG: &str = "chr1";
const CONTIG_LEN: u32 = 20_000;
/// Reads start inside this many bases of the contig start.
const START_WINDOW: u32 = 500;

// ── the generated SAM ───────────────────────────────────────────────────

/// One auxiliary tag, in the SAM spelling and in the form the comparison
/// normalizes to.
#[derive(Debug, Clone)]
enum Aux {
    Char(char),
    Int(i64),
    Float(f32),
    Text(String),
    IntArray(char, Vec<i64>),
}

impl Aux {
    fn sam(&self, tag: &str) -> String {
        match self {
            Self::Char(c) => format!("{tag}:A:{c}"),
            Self::Int(v) => format!("{tag}:i:{v}"),
            Self::Float(v) => format!("{tag}:f:{v}"),
            Self::Text(s) => format!("{tag}:Z:{s}"),
            Self::IntArray(sub, vs) => {
                let mut out = format!("{tag}:B:{sub}");
                for v in vs {
                    write!(out, ",{v}").expect("writing to a String cannot fail");
                }
                out
            }
        }
    }
}

#[derive(Debug, Clone)]
struct Read {
    pos: u32,
    flags: u16,
    mapq: u32,
    cigar: Vec<(u32, char)>,
    seq: String,
    qual: String,
    /// `-1` for "no mate position", which SAM spells `0`.
    next_pos: i64,
    tlen: i64,
    aux: Vec<(String, Aux)>,
}

impl Read {
    fn cigar_text(&self) -> String {
        let mut out = String::new();
        for (len, op) in &self.cigar {
            write!(out, "{len}{op}").expect("writing to a String cannot fail");
        }
        out
    }
}

fn bases(len: usize) -> impl PrintableGenerator<String> {
    gs::text().alphabet("ACGTN").min_size(len).max_size(len)
}

/// A CIGAR valid by construction: matches at both ends of the aligned body,
/// the indel and padding ops between them, soft clips outside those, hard
/// clips outermost — the order the SAM spec requires.
#[hegel::composite]
fn arb_cigar(tc: &TestCase) -> Vec<(u32, char)> {
    let n = |lo: u32, hi: u32| gs::integers::<u32>().min_value(lo).max_value(hi);
    let m = |tc: &TestCase| {
        // `=` and `X` are M's two halves; a parser that folds them together
        // would still round-trip M, so generate all three.
        (tc.draw_silent(n(1, 30)), tc.draw_silent(gs::sampled_from(&['M', '=', 'X'])))
    };

    let mut body = vec![m(tc)];
    for _ in 0..tc.draw_silent(gs::integers::<usize>().max_value(2)) {
        let op = tc.draw_silent(gs::sampled_from(&['I', 'D', 'N', 'P']));
        let len = match op {
            'N' => tc.draw_silent(n(1, 200)),
            'P' => tc.draw_silent(n(1, 3)),
            _ => tc.draw_silent(n(1, 6)),
        };
        body.push((len, op));
        body.push(m(tc));
    }

    let mut ops = Vec::new();
    if tc.draw_silent(gs::booleans()) {
        ops.push((tc.draw_silent(n(1, 4)), 'H'));
    }
    if tc.draw_silent(gs::booleans()) {
        ops.push((tc.draw_silent(n(1, 8)), 'S'));
    }
    ops.extend(body);
    if tc.draw_silent(gs::booleans()) {
        ops.push((tc.draw_silent(n(1, 8)), 'S'));
    }
    if tc.draw_silent(gs::booleans()) {
        ops.push((tc.draw_silent(n(1, 4)), 'H'));
    }
    ops
}

#[hegel::composite]
fn arb_aux(tc: &TestCase) -> Aux {
    match tc.draw_silent(gs::integers::<u8>().max_value(4)) {
        0 => Aux::Char(tc.draw_silent(gs::sampled_from(&['A', 'z', '0', '+', '~']))),
        // Spanning the widths BAM picks between, and both signs: the writer
        // chooses c/C/s/S/i/I, and samtools has to print `i` back whichever
        // it picked.
        1 => Aux::Int(tc.draw_silent(
            gs::integers::<i64>().min_value(i64::from(i32::MIN)).max_value(i64::from(u32::MAX)),
        )),
        2 => Aux::Float(tc.draw_silent(gs::floats::<f32>().min_value(-1e6).max_value(1e6))),
        3 => Aux::Text(tc.draw_silent(gs::from_regex("[A-Za-z0-9 _.-]{0,12}"))),
        _ => {
            let (sub, lo, hi) = match tc.draw_silent(gs::integers::<u8>().max_value(2)) {
                0 => ('c', -128i64, 127i64),
                1 => ('S', 0, 65_535),
                _ => ('i', i64::from(i32::MIN), i64::from(i32::MAX)),
            };
            let len = tc.draw_silent(gs::integers::<usize>().min_value(1).max_value(4));
            let vs = (0..len)
                .map(|_| tc.draw_silent(gs::integers::<i64>().min_value(lo).max_value(hi)))
                .collect();
            Aux::IntArray(sub, vs)
        }
    }
}

#[hegel::composite]
fn arb_read(tc: &TestCase) -> Read {
    let cigar = tc.draw_silent(arb_cigar());
    let span: u32 = cigar
        .iter()
        .filter(|(_, op)| matches!(op, 'M' | 'D' | 'N' | '=' | 'X'))
        .map(|(l, _)| l)
        .sum();
    let last_start = CONTIG_LEN.saturating_sub(span).clamp(1, START_WINDOW);
    let pos = tc.draw_silent(gs::integers::<u32>().min_value(1).max_value(last_start));

    let qlen: usize = cigar
        .iter()
        .filter(|(_, op)| matches!(op, 'M' | 'I' | 'S' | '=' | 'X'))
        .map(|(l, _)| *l as usize)
        .sum();

    // Paired flags, strand, and the ones a reader must carry through without
    // acting on: secondary, supplementary, duplicate, QC-fail.
    let flags = tc.draw_silent(gs::sampled_from(&[
        0u16, 16, 99, 147, 83, 163, 0x100, 0x800, 0x400, 0x200, 0x41, 0x81,
    ]));
    let paired = flags & 0x1 != 0;
    let next_pos = if paired {
        i64::from(tc.draw_silent(gs::integers::<u32>().min_value(1).max_value(START_WINDOW)))
    } else {
        0
    };
    let tlen = if paired {
        i64::from(tc.draw_silent(gs::integers::<i32>().min_value(-600).max_value(600)))
    } else {
        0
    };

    let n_aux = tc.draw_silent(gs::integers::<usize>().max_value(3));
    let aux = (0..n_aux)
        .map(|i| {
            // Distinct two-letter tags, so order and identity are both pinned.
            (format!("X{}", (b'a' + i as u8) as char), tc.draw_silent(arb_aux()))
        })
        .collect();

    Read {
        pos,
        flags,
        mapq: tc.draw_silent(gs::integers::<u32>().max_value(60)),
        cigar,
        seq: tc.draw_silent(bases(qlen)),
        qual: tc.draw_silent(
            gs::text()
                .alphabet("!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI")
                .min_size(qlen)
                .max_size(qlen),
        ),
        next_pos,
        tlen,
        aux,
    }
}

/// Coordinate-sorted, which is what `tabix -p sam` requires.
#[hegel::composite]
fn arb_reads(tc: &TestCase) -> Vec<Read> {
    let mut reads = tc.draw_silent(gs::vecs(arb_read()).min_size(1).max_size(8));
    reads.sort_by_key(|r| r.pos);
    reads
}

fn sam_text(reads: &[Read]) -> String {
    let mut sam = format!("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:{CONTIG}\tLN:{CONTIG_LEN}\n");
    for (i, r) in reads.iter().enumerate() {
        let rnext = if r.flags & 0x1 != 0 { "=" } else { "*" };
        write!(
            sam,
            "r{i}\t{}\t{CONTIG}\t{}\t{}\t{}\t{rnext}\t{}\t{}\t{}\t{}",
            r.flags,
            r.pos,
            r.mapq,
            r.cigar_text(),
            r.next_pos,
            r.tlen,
            r.seq,
            r.qual
        )
        .expect("writing to a String cannot fail");
        for (tag, aux) in &r.aux {
            sam.push('\t');
            sam.push_str(&aux.sam(tag));
        }
        sam.push('\n');
    }
    sam
}

// ── the trip ────────────────────────────────────────────────────────────

fn run(cmd: &mut Command, what: &str) {
    let status = cmd
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .unwrap_or_else(|e| panic!("{what}: {e}"));
    assert!(status.success(), "{what} failed");
}

/// bgzip and tabix the SAM, which is the form seqair's SAM reader takes.
fn indexed_sam(dir: &Path, sam: &str) -> std::path::PathBuf {
    let plain = dir.join("generated.sam");
    std::fs::write(&plain, sam).expect("write SAM");
    let gz = dir.join("generated.sam.gz");
    // Not through `run`: that redirects stdout, which is where bgzip writes.
    let status = Command::new("bgzip")
        .arg("-c")
        .arg(&plain)
        .stdout(std::fs::File::create(&gz).expect("create .sam.gz"))
        .stderr(Stdio::null())
        .status()
        .expect("bgzip not found");
    assert!(status.success(), "bgzip failed");
    run(Command::new("tabix").args(["-p", "sam"]).arg(&gz), "tabix");
    gz
}

/// SAM → seqair's reader → seqair's writer → BAM, then samtools back to SAM.
fn round_trip(dir: &Path, sam: &str) -> String {
    let gz = indexed_sam(dir, sam);

    let mut reader = IndexedSamReader::open(&gz).expect("open indexed SAM");
    let tid = reader.header().tid(CONTIG).expect("contig");
    let mut store = RecordStore::new();
    reader
        .fetch_into(
            tid,
            (Pos0::new(0).unwrap()..=Pos0::new(CONTIG_LEN - 1).unwrap()).into(),
            &mut store,
        )
        .expect("fetch");

    let bam = dir.join("generated.bam");
    {
        let header = reader.header();
        let mut writer = BamWriterBuilder::to_path(&bam, header).build().expect("build writer");
        for idx in store.indices() {
            writer.write_store_record(&store, idx).expect("write record");
        }
        writer.finish().expect("finish");
    }

    let out = Command::new("samtools")
        .args(["view", "-h"])
        .arg(&bam)
        .output()
        .expect("samtools not found");
    assert!(out.status.success(), "samtools view failed: {}", String::from_utf8_lossy(&out.stderr));
    String::from_utf8(out.stdout).expect("samtools emits UTF-8")
}

// ── comparison ──────────────────────────────────────────────────────────

/// A SAM alignment line, normalized so that the spellings a round-trip is
/// allowed to change do not count as differences.
#[derive(Debug, Clone, PartialEq)]
struct Line {
    qname: String,
    flags: u16,
    rname: String,
    pos: u32,
    mapq: u32,
    cigar: String,
    /// `=` resolved to the contig it means.
    rnext: String,
    pnext: i64,
    tlen: i64,
    seq: String,
    qual: String,
    /// Integer widths are the writer's choice (`r[unified.fetch_equivalence]`),
    /// so `c`/`C`/`s`/`S`/`i`/`I` all normalize to one integer, and floats are
    /// compared at the six significant digits SAM text carries.
    aux: BTreeMap<String, String>,
}

/// Six significant digits, which is all a SAM text float carries: samtools
/// prints aux floats with `%g`, so `1.015625` comes back as `1.01562`. The
/// BAM in the middle holds the exact `f32` — the loss is in the last hop, and
/// rendering both sides the same way compares what the channel can carry
/// rather than pretending it carries more.
fn sam_float(value: &str) -> String {
    let v: f32 = value.parse().expect("a SAM float field is a float");
    format!("{v:.5e}")
}

fn normalize_aux(tag: &str, type_char: char, value: &str) -> (String, String) {
    let normalized = match type_char {
        'f' => format!("f:{}", sam_float(value)),
        'B' => {
            let mut parts = value.splitn(2, ',');
            let sub = parts.next().unwrap_or("");
            let rest = parts.next().unwrap_or("");
            if sub == "f" {
                let rendered: Vec<String> = rest.split(',').map(sam_float).collect();
                format!("B:f:{}", rendered.join(","))
            } else {
                // The subtype is the writer's choice the same way a scalar
                // integer's is; only the values are data.
                format!("B:i:{rest}")
            }
        }
        other => format!("{other}:{value}"),
    };
    (tag.to_string(), normalized)
}

fn parse_lines(sam: &str) -> Vec<Line> {
    sam.lines()
        .filter(|l| !l.starts_with('@'))
        .map(|line| {
            let f: Vec<&str> = line.split('\t').collect();
            let rname = f[2].to_string();
            let rnext = match f[6] {
                "=" => rname.clone(),
                other => other.to_string(),
            };
            let aux = f[11..]
                .iter()
                .map(|field| {
                    let tag = &field[0..2];
                    let type_char = field.as_bytes()[3] as char;
                    normalize_aux(tag, type_char, &field[5..])
                })
                .collect();
            Line {
                qname: f[0].to_string(),
                flags: f[1].parse().expect("FLAG"),
                rname,
                pos: f[3].parse().expect("POS"),
                mapq: f[4].parse().expect("MAPQ"),
                cigar: f[5].to_string(),
                rnext,
                pnext: f[7].parse().expect("PNEXT"),
                tlen: f[8].parse().expect("TLEN"),
                seq: f[9].to_string(),
                qual: f[10].to_string(),
                aux,
            }
        })
        .collect()
}

// r[verify unified.fetch_equivalence]
// r[verify sam.record.aux_tags]
// r[verify bam_writer.write_store_record]
/// Every field of every record must come back unchanged from
/// SAM → seqair reader → seqair writer → BAM → samtools → SAM.
#[hegel::test(test_cases = 48)]
fn sam_survives_the_trip_through_seqairs_bam_writer(tc: TestCase) {
    let reads = tc.draw(arb_reads().print_as_debug());
    let original = sam_text(&reads);

    let dir = tempfile::tempdir().expect("tempdir");
    let round_tripped = round_trip(dir.path(), &original);

    let want = parse_lines(&original);
    let got = parse_lines(&round_tripped);

    assert_eq!(got.len(), want.len(), "record count\n--- in ---\n{original}");
    for (got, want) in got.iter().zip(&want) {
        assert_eq!(got, want, "record {}\n--- in ---\n{original}", want.qname);
    }

    tc.event_value("reads", reads.len() as f64);
    if reads.iter().any(|r| !r.aux.is_empty()) {
        tc.event("has aux tags");
    }
    if reads.iter().any(|r| r.flags & 0x1 != 0) {
        tc.event("has a paired read");
    }
    if reads.iter().any(|r| r.cigar.iter().any(|(_, op)| matches!(op, 'H' | 'P' | 'N'))) {
        tc.event("has an H, P or N op");
    }
}
