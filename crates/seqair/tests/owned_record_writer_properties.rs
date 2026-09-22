//! `OwnedBamRecord` built through its builder, written by `BamWriter`, read
//! back by samtools.
//!
//! This is the other half of the write path. `write_store_record` copies a
//! record the reader already decoded — slab bytes straight through — and
//! `sam_bam_roundtrip_properties` covers that. `write(&OwnedBamRecord)` is the
//! path a caller who *builds* a record takes: `AuxData`'s typed setters choose
//! BAM type codes and append bytes, `to_bam_bytes` lays out the binary record,
//! and nothing the reader produced is involved. Its tests are fixtures.
//!
//! samtools reads the result, and the values the generator chose are the
//! oracle for what it should say.
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
use seqair::bam::aux_data::AuxData;
use seqair::bam::cigar::{CigarOp, CigarOpType};
use seqair::bam::header::BamHeader;
use seqair::bam::owned_record::OwnedBamRecord;
use seqair::bam::writer::BamWriterBuilder;
use seqair::bam::{Pos0, RecordStore};
use seqair_types::{BamFlags, Base, BaseQuality};
use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::Path;
use std::process::Command;

const CONTIGS: [(&str, u32); 2] = [("chr1", 40_000), ("chr2", 10_000)];
/// Records start inside this many bases of a contig's start.
const START_WINDOW: u32 = 400;

fn header() -> BamHeader {
    let mut text = String::from("@HD\tVN:1.6\tSO:coordinate\n");
    for (name, len) in CONTIGS {
        writeln!(text, "@SQ\tSN:{name}\tLN:{len}").expect("writing to a String cannot fail");
    }
    BamHeader::from_sam_text(&text).unwrap()
}

// ── what the generator decides ──────────────────────────────────────────

/// One aux tag, as the setter to call and as the SAM text samtools must print.
#[derive(Debug, Clone)]
enum Aux {
    Char(u8),
    Int(i64),
    Float(f32),
    Text(Vec<u8>),
    Hex(Vec<u8>),
    ArrayI8(Vec<i8>),
    ArrayU8(Vec<u8>),
    ArrayI16(Vec<i16>),
    ArrayU16(Vec<u16>),
    ArrayI32(Vec<i32>),
    ArrayU32(Vec<u32>),
    ArrayF32(Vec<f32>),
}

impl Aux {
    fn write_into(&self, aux: &mut AuxData, tag: [u8; 2]) {
        match self {
            Self::Char(c) => aux.set_char(tag, *c).expect("printable ASCII"),
            Self::Int(v) => aux.set_int(tag, *v).expect("inside the i32/u32 union"),
            Self::Float(v) => aux.set_float(tag, *v),
            Self::Text(s) => aux.set_string(tag, s),
            Self::Hex(bytes) => {
                let mut hex = Vec::with_capacity(bytes.len() * 2);
                for b in bytes {
                    write!(&mut HexSink(&mut hex), "{b:02X}").expect("cannot fail");
                }
                aux.set_hex(tag, &hex);
            }
            Self::ArrayI8(v) => aux.set_array_i8(tag, v).expect("short array"),
            Self::ArrayU8(v) => aux.set_array_u8(tag, v).expect("short array"),
            Self::ArrayI16(v) => aux.set_array_i16(tag, v).expect("short array"),
            Self::ArrayU16(v) => aux.set_array_u16(tag, v).expect("short array"),
            Self::ArrayI32(v) => aux.set_array_i32(tag, v).expect("short array"),
            Self::ArrayU32(v) => aux.set_array_u32(tag, v).expect("short array"),
            Self::ArrayF32(v) => aux.set_array_f32(tag, v).expect("short array"),
        }
    }

    /// Which setter this exercises, for the statistics: a generator that
    /// stopped producing one of the twelve would otherwise do so silently.
    fn kind(&self) -> &'static str {
        match self {
            Self::Char(_) => "A",
            Self::Int(_) => "i",
            Self::Float(_) => "f",
            Self::Text(_) => "Z",
            Self::Hex(_) => "H",
            Self::ArrayI8(_) => "B:c",
            Self::ArrayU8(_) => "B:C",
            Self::ArrayI16(_) => "B:s",
            Self::ArrayU16(_) => "B:S",
            Self::ArrayI32(_) => "B:i",
            Self::ArrayU32(_) => "B:I",
            Self::ArrayF32(_) => "B:f",
        }
    }

    /// What the tag must read back as, normalized the way `parse_aux` below
    /// normalizes what samtools prints.
    fn expected(&self) -> String {
        let join = |vs: Vec<String>| format!("B:i:{}", vs.join(","));
        match self {
            Self::Char(c) => format!("A:{}", *c as char),
            Self::Int(v) => format!("i:{v}"),
            Self::Float(v) => format!("f:{}", sam_float_text(*v)),
            Self::Text(s) => format!("Z:{}", String::from_utf8_lossy(s)),
            Self::Hex(bytes) => {
                let mut hex = String::new();
                for b in bytes {
                    write!(hex, "{b:02X}").expect("cannot fail");
                }
                format!("H:{hex}")
            }
            Self::ArrayI8(v) => join(v.iter().map(i8::to_string).collect()),
            Self::ArrayU8(v) => join(v.iter().map(u8::to_string).collect()),
            Self::ArrayI16(v) => join(v.iter().map(i16::to_string).collect()),
            Self::ArrayU16(v) => join(v.iter().map(u16::to_string).collect()),
            Self::ArrayI32(v) => join(v.iter().map(i32::to_string).collect()),
            Self::ArrayU32(v) => join(v.iter().map(u32::to_string).collect()),
            Self::ArrayF32(v) => {
                format!(
                    "B:f:{}",
                    v.iter().map(|f| sam_float_text(*f)).collect::<Vec<_>>().join(",")
                )
            }
        }
    }
}

/// A `fmt::Write` that appends to a byte buffer.
struct HexSink<'a>(&'a mut Vec<u8>);
impl std::fmt::Write for HexSink<'_> {
    fn write_str(&mut self, s: &str) -> std::fmt::Result {
        self.0.extend_from_slice(s.as_bytes());
        Ok(())
    }
}

/// Six significant digits, which is all SAM text carries: samtools prints aux
/// floats with `%g`, so the exact `f32` in the BAM comes back rounded. Both
/// sides render the same way, comparing what the channel can carry.
fn sam_float_text(v: f32) -> String {
    format!("{v:.5e}")
}

#[derive(Debug, Clone)]
struct Record {
    tid: usize,
    pos: u32,
    flags: u16,
    mapq: u8,
    cigar: Vec<(u32, CigarOpType)>,
    seq: Vec<Base>,
    qual: Vec<u8>,
    next_ref_id: i32,
    next_pos: Option<u32>,
    template_len: i32,
    aux: Vec<([u8; 2], Aux)>,
}

impl Record {
    fn cigar_text(&self) -> String {
        if self.cigar.is_empty() {
            return "*".to_string();
        }
        let mut out = String::new();
        for (len, op) in &self.cigar {
            write!(out, "{len}{op}").expect("writing to a String cannot fail");
        }
        out
    }

    fn ref_span(&self) -> u32 {
        self.cigar
            .iter()
            .filter(|(_, op)| {
                matches!(
                    op,
                    CigarOpType::Match
                        | CigarOpType::Deletion
                        | CigarOpType::RefSkip
                        | CigarOpType::SeqMatch
                        | CigarOpType::SeqMismatch
                )
            })
            .map(|(l, _)| l)
            .sum()
    }
}

fn short<T: Clone + Send + Sync + hegel::PrettyPrintable + 'static>(
    element: impl PrintableGenerator<T> + Send + Sync + 'static,
) -> impl PrintableGenerator<Vec<T>> {
    gs::vecs(element).min_size(1).max_size(4)
}

#[hegel::composite]
fn arb_aux(tc: &TestCase) -> Aux {
    let ints = |lo: i32, hi: i32| gs::integers::<i32>().min_value(lo).max_value(hi);
    match tc.draw_silent(gs::integers::<u8>().max_value(11)) {
        // `set_char` takes printable ASCII only.
        0 => Aux::Char(tc.draw_silent(gs::integers::<u8>().min_value(b'!').max_value(b'~'))),
        // The union of i32 and u32, which is what `set_int` accepts and where
        // it picks between six BAM type codes.
        1 => Aux::Int(tc.draw_silent(
            gs::integers::<i64>().min_value(i64::from(i32::MIN)).max_value(i64::from(u32::MAX)),
        )),
        2 => Aux::Float(tc.draw_silent(gs::floats::<f32>().min_value(-1e6).max_value(1e6))),
        3 => Aux::Text(tc.draw_silent(gs::from_regex("[A-Za-z0-9 _.-]{0,10}")).into_bytes()),
        4 => Aux::Hex(tc.draw_silent(gs::binary().max_size(6))),
        5 => Aux::ArrayI8(tc.draw_silent(short(gs::integers::<i8>())).into_iter().collect()),
        6 => Aux::ArrayU8(tc.draw_silent(short(gs::integers::<u8>()))),
        7 => Aux::ArrayI16(tc.draw_silent(short(gs::integers::<i16>()))),
        8 => Aux::ArrayU16(tc.draw_silent(short(gs::integers::<u16>()))),
        9 => Aux::ArrayI32(tc.draw_silent(short(ints(i32::MIN, i32::MAX)))),
        10 => Aux::ArrayU32(tc.draw_silent(short(gs::integers::<u32>()))),
        _ => {
            Aux::ArrayF32(tc.draw_silent(short(gs::floats::<f32>().min_value(-1e6).max_value(1e6))))
        }
    }
}

#[hegel::composite]
fn arb_record(tc: &TestCase) -> Record {
    let n = |lo: u32, hi: u32| gs::integers::<u32>().min_value(lo).max_value(hi);
    let tid = tc.draw_silent(gs::integers::<usize>().max_value(CONTIGS.len() - 1));
    let contig_len = CONTIGS[tid].1;

    // Matches at both ends of the aligned body, indels and skips between,
    // clips outside — the order the BAM spec requires.
    let mut cigar: Vec<(u32, CigarOpType)> = Vec::new();
    if tc.draw_silent(gs::booleans()) {
        cigar.push((tc.draw_silent(n(1, 4)), CigarOpType::HardClip));
    }
    if tc.draw_silent(gs::booleans()) {
        cigar.push((tc.draw_silent(n(1, 6)), CigarOpType::SoftClip));
    }
    let matchish = |tc: &TestCase| {
        (
            tc.draw_silent(n(1, 25)),
            tc.draw_silent(gs::sampled_from(&[
                CigarOpType::Match,
                CigarOpType::SeqMatch,
                CigarOpType::SeqMismatch,
            ])),
        )
    };
    cigar.push(matchish(tc));
    for _ in 0..tc.draw_silent(gs::integers::<usize>().max_value(2)) {
        let (op, len) = match tc.draw_silent(gs::integers::<u8>().max_value(3)) {
            0 => (CigarOpType::Insertion, tc.draw_silent(n(1, 5))),
            1 => (CigarOpType::Deletion, tc.draw_silent(n(1, 6))),
            2 => (CigarOpType::RefSkip, tc.draw_silent(n(1, 300))),
            _ => (CigarOpType::Padding, tc.draw_silent(n(1, 3))),
        };
        cigar.push((len, op));
        cigar.push(matchish(tc));
    }
    if tc.draw_silent(gs::booleans()) {
        cigar.push((tc.draw_silent(n(1, 6)), CigarOpType::SoftClip));
    }
    if tc.draw_silent(gs::booleans()) {
        cigar.push((tc.draw_silent(n(1, 4)), CigarOpType::HardClip));
    }

    let span: u32 = cigar
        .iter()
        .filter(|(_, op)| {
            matches!(
                op,
                CigarOpType::Match
                    | CigarOpType::Deletion
                    | CigarOpType::RefSkip
                    | CigarOpType::SeqMatch
                    | CigarOpType::SeqMismatch
            )
        })
        .map(|(l, _)| l)
        .sum();
    let pos = tc.draw_silent(n(0, contig_len.saturating_sub(span).clamp(1, START_WINDOW) - 1));

    let qlen: usize = cigar
        .iter()
        .filter(|(_, op)| {
            matches!(
                op,
                CigarOpType::Match
                    | CigarOpType::Insertion
                    | CigarOpType::SoftClip
                    | CigarOpType::SeqMatch
                    | CigarOpType::SeqMismatch
            )
        })
        .map(|(l, _)| *l as usize)
        .sum();

    let flags =
        tc.draw_silent(gs::sampled_from(&[0u16, 16, 99, 147, 83, 163, 0x100, 0x800, 0x400, 0x200]));
    let paired = flags & 0x1 != 0;

    let n_aux = tc.draw_silent(gs::integers::<usize>().max_value(4));
    let aux = (0..n_aux).map(|i| ([b'X', b'a' + i as u8], tc.draw_silent(arb_aux()))).collect();

    Record {
        tid,
        pos,
        flags,
        mapq: tc.draw_silent(gs::integers::<u8>().max_value(60)),
        cigar,
        seq: (0..qlen)
            .map(|_| {
                tc.draw_silent(gs::sampled_from(&[
                    Base::A,
                    Base::C,
                    Base::G,
                    Base::T,
                    Base::Unknown,
                ]))
            })
            .collect(),
        qual: (0..qlen).map(|_| tc.draw_silent(gs::integers::<u8>().max_value(40))).collect(),
        next_ref_id: if paired { tid as i32 } else { -1 },
        next_pos: paired.then(|| tc.draw_silent(n(0, START_WINDOW))),
        template_len: if paired {
            tc.draw_silent(gs::integers::<i32>().min_value(-500).max_value(500))
        } else {
            0
        },
        aux,
    }
}

/// Coordinate-sorted, which `BamWriter`'s index co-production requires.
#[hegel::composite]
fn arb_records(tc: &TestCase) -> Vec<Record> {
    let mut records = tc.draw_silent(gs::vecs(arb_record()).min_size(1).max_size(8));
    records.sort_by_key(|r| (r.tid, r.pos));
    records
}

// ── the trip ────────────────────────────────────────────────────────────

fn write_bam(dir: &Path, records: &[Record]) -> std::path::PathBuf {
    let header = header();
    let path = dir.join("built.bam");
    let mut writer =
        BamWriterBuilder::to_path(&path, &header).write_index(true).build().expect("build writer");

    for (i, rec) in records.iter().enumerate() {
        let mut aux = AuxData::new();
        for (tag, value) in &rec.aux {
            value.write_into(&mut aux, *tag);
        }
        let built = OwnedBamRecord::builder(
            rec.tid as i32,
            Some(Pos0::new(rec.pos).unwrap()),
            format!("r{i}").into_bytes(),
        )
        .flags(BamFlags::from(rec.flags))
        .mapq(rec.mapq)
        .cigar(rec.cigar.iter().map(|&(len, op)| CigarOp::new(op, len)).collect())
        .seq(rec.seq.clone())
        .qual(rec.qual.iter().copied().map(BaseQuality::from_byte).collect())
        .next_ref_id(rec.next_ref_id)
        .next_pos(rec.next_pos.map(|p| Pos0::new(p).unwrap()))
        .template_len(rec.template_len)
        .aux(aux)
        .build()
        .expect("synthetic record must build");
        writer.write(&built).expect("write");
    }
    let (_inner, index) = writer.finish().expect("finish");
    let bai = std::fs::File::create(path.with_extension("bam.bai")).expect("create BAI");
    index.expect("write_index(true)").write_bai(bai, header.target_count()).expect("write BAI");
    path
}

/// One SAM alignment line as samtools prints it, normalized for the spellings
/// a BAM round-trip is allowed to choose.
#[derive(Debug, Clone, PartialEq)]
struct Line {
    qname: String,
    flags: u16,
    rname: String,
    pos: u32,
    mapq: u8,
    cigar: String,
    rnext: String,
    pnext: i64,
    tlen: i64,
    seq: String,
    qual: String,
    aux: BTreeMap<String, String>,
}

fn parse_aux(field: &str) -> (String, String) {
    let tag = field[0..2].to_string();
    let type_char = field.as_bytes()[3] as char;
    let value = &field[5..];
    let normalized = match type_char {
        // BAM picks the width; only the number is data.
        'c' | 'C' | 's' | 'S' | 'i' | 'I' => format!("i:{value}"),
        'f' => format!("f:{}", sam_float_text(value.parse().expect("SAM float"))),
        'B' => {
            let (sub, rest) = value.split_at(1);
            let rest = rest.trim_start_matches(',');
            if sub == "f" {
                let rendered: Vec<String> = rest
                    .split(',')
                    .map(|v| sam_float_text(v.parse().expect("SAM float")))
                    .collect();
                format!("B:f:{}", rendered.join(","))
            } else {
                format!("B:i:{rest}")
            }
        }
        other => format!("{other}:{value}"),
    };
    (tag, normalized)
}

fn samtools_lines(path: &Path) -> Vec<Line> {
    let out = Command::new("samtools").arg("view").arg(path).output().expect("samtools not found");
    assert!(out.status.success(), "samtools view failed: {}", String::from_utf8_lossy(&out.stderr));
    std::str::from_utf8(&out.stdout)
        .expect("samtools emits UTF-8")
        .lines()
        .filter(|l| !l.is_empty())
        .map(|line| {
            let f: Vec<&str> = line.split('\t').collect();
            let rname = f[2].to_string();
            Line {
                qname: f[0].to_string(),
                flags: f[1].parse().expect("FLAG"),
                rname,
                pos: f[3].parse().expect("POS"),
                mapq: f[4].parse().expect("MAPQ"),
                cigar: f[5].to_string(),
                rnext: f[6].to_string(),
                pnext: f[7].parse().expect("PNEXT"),
                tlen: f[8].parse().expect("TLEN"),
                seq: f[9].to_string(),
                qual: f[10].to_string(),
                aux: f[11..].iter().map(|field| parse_aux(field)).collect(),
            }
        })
        .collect()
}

/// What samtools must print for a record the generator described.
fn expected_line(i: usize, rec: &Record) -> Line {
    let paired = rec.flags & 0x1 != 0;
    Line {
        qname: format!("r{i}"),
        flags: rec.flags,
        rname: CONTIGS[rec.tid].0.to_string(),
        pos: rec.pos + 1, // SAM is 1-based
        mapq: rec.mapq,
        cigar: rec.cigar_text(),
        // `next_ref_id` equal to this record's is printed `=`.
        rnext: if paired { "=".to_string() } else { "*".to_string() },
        pnext: rec.next_pos.map_or(0, |p| i64::from(p) + 1),
        tlen: i64::from(rec.template_len),
        seq: rec.seq.iter().map(|b| **b as char).collect(),
        qual: rec.qual.iter().map(|q| (q + 33) as char).collect(),
        aux: rec
            .aux
            .iter()
            .map(|(tag, value)| (String::from_utf8(tag.to_vec()).expect("ASCII"), value.expected()))
            .collect(),
    }
}

/// A record as seqair's own reader gives it back.
#[derive(Debug, Clone, PartialEq)]
struct Decoded {
    pos: u32,
    flags: u16,
    mapq: u8,
    seq: Vec<Base>,
    qual: Vec<u8>,
    next_ref_id: i32,
    template_len: i32,
}

// r[verify bam.owned_record.aux_data]
// r[verify bam.owned_record.to_bam_bytes]
// r[verify bam_writer.write_record]
/// A record built through the builder and written by `BamWriter` must read
/// back through samtools as the values the generator chose.
#[hegel::test(test_cases = 48)]
fn built_records_read_back_through_samtools(tc: TestCase) {
    let records = tc.draw(arb_records().print_as_debug());

    let dir = tempfile::tempdir().expect("tempdir");
    let path = write_bam(dir.path(), &records);

    let got = samtools_lines(&path);
    assert_eq!(got.len(), records.len(), "record count");
    for (i, (got, rec)) in got.iter().zip(&records).enumerate() {
        assert_eq!(*got, expected_line(i, rec), "record r{i}");
    }

    tc.event_value("records", records.len() as f64);
    tc.event_value("ref span", records.iter().map(Record::ref_span).max().unwrap_or(0).into());
    if records.iter().any(|r| r.flags & 0x1 != 0) {
        tc.event("has a paired record");
    }
    for kind in records.iter().flat_map(|r| r.aux.iter().map(|(_, a)| a.kind())) {
        tc.event(kind);
    }
}

// r[verify bam.owned_record.to_bam_bytes]
// r[verify unified.fetch_equivalence]
/// And seqair must read its own output back the same way samtools does.
#[hegel::test(test_cases = 48)]
fn built_records_read_back_through_seqair(tc: TestCase) {
    let records = tc.draw(arb_records().print_as_debug());

    let dir = tempfile::tempdir().expect("tempdir");
    let path = write_bam(dir.path(), &records);

    // The records are sorted by `(tid, pos)`, so fetching the contigs in order
    // yields them in the order they were written.
    let mut reader = seqair::bam::IndexedBamReader::open(&path).expect("open");
    let mut all: Vec<Decoded> = Vec::new();
    for (name, len) in CONTIGS {
        let tid = reader.header().tid(name).expect("contig");
        let mut store = RecordStore::new();
        reader
            .fetch_into(tid, Pos0::new(0).unwrap(), Pos0::new(len - 1).unwrap(), &mut store)
            .expect("fetch");
        for idx in store.indices() {
            let r = store.record(idx).unwrap();
            all.push(Decoded {
                pos: r.pos.as_u32(),
                flags: r.flags.raw(),
                mapq: r.mapq,
                seq: r.seq().to_vec(),
                qual: BaseQuality::slice_to_bytes(r.qual()).to_vec(),
                next_ref_id: r.next_ref_id,
                template_len: r.template_len,
            });
        }
    }

    assert_eq!(all.len(), records.len(), "record count");
    for (got, rec) in all.iter().zip(&records) {
        assert_eq!(got.pos, rec.pos, "pos");
        assert_eq!(got.flags, rec.flags, "flags");
        assert_eq!(got.mapq, rec.mapq, "mapq");
        assert_eq!(got.seq, rec.seq, "seq");
        assert_eq!(got.qual, rec.qual, "qual");
        assert_eq!(got.next_ref_id, rec.next_ref_id, "next_ref_id");
        assert_eq!(got.template_len, rec.template_len, "template_len");
    }
}
