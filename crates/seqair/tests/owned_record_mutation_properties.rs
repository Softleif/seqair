//! `OwnedBamRecord` *mutated* after construction, then written and read back.
//!
//! `owned_record_writer_properties` builds a record and writes it once; every
//! field it checks was decided by the builder. This file covers the other
//! half of the owned record's contract — the one the TAPS rewrite pipeline
//! actually uses. A record is built, then pushed through a generated sequence
//! of mutations: SEQ and QUAL replaced at new lengths, flags rewritten, aux
//! tags set, re-set with a different type *and* a different encoded length,
//! removed, set again. A small in-test model (expected SEQ, QUAL, flags, and
//! a `BTreeMap` of tag to value) is the oracle; samtools, reading a snapshot
//! of the record after every single step, is the referee.
//!
//! The mutators are where the invariants are easy to break and hard to
//! notice. `set_seq` has to keep the 4-bit packing honest across an odd/even
//! length change; `set_qual` has to stay pinned to `seq.len()`; the aux
//! setters remove-and-re-append, so a tag written twice must leave exactly
//! one copy and no orphaned bytes behind.
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
use seqair::bam::Pos0;
use seqair::bam::aux_data::AuxData;
use seqair::bam::header::BamHeader;
use seqair::bam::owned_record::OwnedBamRecord;
use seqair::bam::writer::BamWriterBuilder;
use seqair_types::{BamFlags, Base, BaseQuality};
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use std::process::Command;

const CONTIG: &str = "chr1";
/// The mutation sequences never move the record, so a short contig is enough.
const SHORT_CONTIG_LEN: u32 = 100_000;
/// The three tags every aux mutation draws from. Keeping the pool small is
/// the point — collisions are what `aux_replace_semantics` is about.
const TAGS: [[u8; 2]; 3] = [*b"Xa", *b"Xb", *b"Xc"];
/// Flag words a real BAM carries: unpaired, reverse, placed-unmapped, the two
/// halves of a proper pair, secondary, duplicate.
const FLAGS: [u16; 8] = [0x0, 0x10, 0x4, 0x14, 0x63, 0x93, 0x100, 0x400];

fn header(contig_len: u32) -> BamHeader {
    let text = format!("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:{CONTIG}\tLN:{contig_len}\n");
    BamHeader::from_sam_text(&text).unwrap()
}

fn write_records(dir: &Path, header: &BamHeader, records: &[OwnedBamRecord]) -> PathBuf {
    let path = dir.join("mutated.bam");
    let mut writer = BamWriterBuilder::to_path(&path, header).build().expect("build writer");
    for rec in records {
        writer.write(rec).expect("write");
    }
    writer.finish().expect("finish");
    path
}

// ── aux values ──────────────────────────────────────────────────────────

/// An aux value, as the setter to call and as the SAM text samtools prints.
///
/// Five variants with five different encoded lengths: replacing an `A` with a
/// `B:C` array grows the tag, replacing a `Z` with an `A` shrinks it. That is
/// what makes `set_*` exercise the remove-then-append path rather than an
/// in-place overwrite.
#[derive(Debug, Clone, PartialEq)]
enum Aux {
    Char(u8),
    Int(i64),
    Float(f32),
    Text(Vec<u8>),
    ArrayU8(Vec<u8>),
}

impl Aux {
    fn write_into(&self, aux: &mut AuxData, tag: [u8; 2]) {
        match self {
            Self::Char(c) => aux.set_char(tag, *c).expect("printable ASCII"),
            Self::Int(v) => aux.set_int(tag, *v).expect("inside the i32/u32 union"),
            Self::Float(v) => aux.set_float(tag, *v),
            Self::Text(s) => aux.set_string(tag, s),
            Self::ArrayU8(v) => aux.set_array_u8(tag, v).expect("short array"),
        }
    }

    /// What the tag must read back as, normalized the way `parse_aux`
    /// normalizes what samtools prints.
    fn expected(&self) -> String {
        match self {
            Self::Char(c) => format!("A:{}", *c as char),
            Self::Int(v) => format!("i:{v}"),
            Self::Float(v) => format!("f:{}", sam_float_text(*v)),
            Self::Text(s) => format!("Z:{}", std::str::from_utf8(s).expect("ASCII")),
            Self::ArrayU8(v) => {
                format!("B:i:{}", v.iter().map(u8::to_string).collect::<Vec<_>>().join(","))
            }
        }
    }

    fn kind(&self) -> &'static str {
        match self {
            Self::Char(_) => "A",
            Self::Int(_) => "i",
            Self::Float(_) => "f",
            Self::Text(_) => "Z",
            Self::ArrayU8(_) => "B:C",
        }
    }
}

/// Six significant digits, which is all SAM text carries. Both sides render
/// the same way, comparing what the channel can carry.
fn sam_float_text(v: f32) -> String {
    format!("{v:.5e}")
}

#[hegel::composite]
fn arb_aux(tc: &TestCase) -> Aux {
    match tc.draw_silent(gs::integers::<u8>().max_value(4)) {
        // `set_char` takes printable ASCII only.
        0 => Aux::Char(tc.draw_silent(gs::integers::<u8>().min_value(b'!').max_value(b'~'))),
        // The union of i32 and u32, which is what `set_int` accepts.
        1 => Aux::Int(tc.draw_silent(
            gs::integers::<i64>().min_value(i64::from(i32::MIN)).max_value(i64::from(u32::MAX)),
        )),
        2 => Aux::Float(tc.draw_silent(gs::floats::<f32>().min_value(-1e6).max_value(1e6))),
        3 => Aux::Text(tc.draw_silent(gs::from_regex("[A-Za-z0-9_.-]{0,8}")).into_bytes()),
        _ => Aux::ArrayU8(tc.draw_silent(gs::vecs(gs::integers::<u8>()).max_size(4))),
    }
}

// ── the mutations ───────────────────────────────────────────────────────

/// One step of a mutation sequence.
///
/// `Resize` carries both halves of the pair because that is the only coherent
/// way to change a read's length: `set_seq` re-anchors what `set_qual` will
/// accept, so the two calls belong together.
#[derive(Debug, Clone)]
enum Op {
    Resize {
        seq: Vec<Base>,
        qual: Vec<u8>,
    },
    /// New quality values at the current length, whatever that has become.
    Requal {
        fill: u8,
    },
    /// Drop the quality array entirely — BAM's "unavailable" spelling.
    ClearQual,
    SetAux {
        tag: [u8; 2],
        value: Aux,
    },
    RemoveAux {
        tag: [u8; 2],
    },
    SetFlags {
        flags: u16,
    },
}

/// What the record must look like after the ops applied so far.
#[derive(Debug, Clone, Default)]
struct Model {
    seq: Vec<Base>,
    /// Empty means "no quality array", which BAM writes as 0xFF fill and
    /// samtools prints as `*`.
    qual: Vec<u8>,
    flags: u16,
    aux: BTreeMap<[u8; 2], Aux>,
}

fn base_of(tc: &TestCase) -> Base {
    tc.draw_silent(gs::sampled_from(&[Base::A, Base::C, Base::G, Base::T, Base::Unknown]))
}

#[hegel::composite]
fn arb_op(tc: &TestCase) -> Op {
    let tag =
        |tc: &TestCase| TAGS[tc.draw_silent(gs::integers::<usize>().max_value(TAGS.len() - 1))];
    match tc.draw_silent(gs::integers::<u8>().max_value(7)) {
        // Weighted towards aux work: two of the eight arms set a tag, one
        // removes one. Lengths 0..=9 straddle the 4-bit packing boundary in
        // both directions.
        0 | 1 => {
            let n = tc.draw_silent(gs::integers::<usize>().max_value(9));
            Op::Resize {
                seq: (0..n).map(|_| base_of(tc)).collect(),
                qual: (0..n).map(|_| tc.draw_silent(gs::integers::<u8>().max_value(40))).collect(),
            }
        }
        2 => Op::Requal { fill: tc.draw_silent(gs::integers::<u8>().max_value(40)) },
        3 => Op::ClearQual,
        4 | 5 => Op::SetAux { tag: tag(tc), value: tc.draw_silent(arb_aux()) },
        6 => Op::RemoveAux { tag: tag(tc) },
        _ => Op::SetFlags { flags: tc.draw_silent(gs::sampled_from(&FLAGS)) },
    }
}

#[hegel::composite]
fn arb_ops(tc: &TestCase) -> Vec<Op> {
    tc.draw_silent(gs::vecs(arb_op()).min_size(1).max_size(8))
}

/// Apply one op to both the record and the model. Every call here is one the
/// mutator's contract accepts; the rejected calls live in their own tests.
fn apply(op: &Op, rec: &mut OwnedBamRecord, model: &mut Model) {
    match op {
        Op::Resize { seq, qual } => {
            rec.set_seq(seq.clone()).expect("no CIGAR to contradict the new length");
            rec.set_qual(qual.iter().copied().map(BaseQuality::from_byte).collect())
                .expect("qual length matches the sequence just set");
            model.seq = seq.clone();
            model.qual = qual.clone();
        }
        Op::Requal { fill } => {
            let qual: Vec<u8> =
                (0..model.seq.len()).map(|i| ((usize::from(*fill) + i) % 41) as u8).collect();
            rec.set_qual(qual.iter().copied().map(BaseQuality::from_byte).collect())
                .expect("qual length matches the current sequence");
            model.qual = qual;
        }
        Op::ClearQual => {
            rec.set_qual(Vec::new()).expect("an empty qual is always allowed");
            model.qual.clear();
        }
        Op::SetAux { tag, value } => {
            value.write_into(&mut rec.aux, *tag);
            model.aux.insert(*tag, value.clone());
        }
        Op::RemoveAux { tag } => {
            rec.aux.remove(*tag);
            model.aux.remove(tag);
        }
        Op::SetFlags { flags } => {
            rec.set_flags(BamFlags::from(*flags));
            model.flags = *flags;
        }
    }
}

/// The starting record: placed, no CIGAR, no sequence, no tags. Everything
/// the properties look at is put there by the mutators.
fn seed_record() -> OwnedBamRecord {
    OwnedBamRecord::builder(0, Some(Pos0::new(1_000).unwrap()), b"seed".to_vec())
        .flags(BamFlags::from(0x4))
        .mapq(17)
        .build()
        .unwrap()
}

// ── samtools as referee ─────────────────────────────────────────────────

/// One SAM alignment line, normalized for the spellings a BAM round-trip is
/// allowed to choose.
#[derive(Debug, Clone, PartialEq)]
struct Line {
    qname: String,
    flags: u16,
    seq: String,
    qual: String,
    aux: BTreeMap<String, String>,
    /// How many aux fields the line actually carried. `aux` is a map, so a
    /// duplicated tag would silently collapse into one entry; this catches it.
    n_aux: usize,
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
            let (_sub, rest) = value.split_at(1);
            format!("B:i:{}", rest.trim_start_matches(','))
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
            Line {
                qname: f[0].to_string(),
                flags: f[1].parse().expect("FLAG"),
                seq: f[9].to_string(),
                qual: f[10].to_string(),
                aux: f[11..].iter().map(|field| parse_aux(field)).collect(),
                n_aux: f.len() - 11,
            }
        })
        .collect()
}

/// What samtools must print for a record in the state the model describes.
fn expected_line(step: usize, model: &Model) -> Line {
    Line {
        qname: format!("s{step}"),
        flags: model.flags,
        seq: if model.seq.is_empty() {
            "*".to_string()
        } else {
            model.seq.iter().map(|b| **b as char).collect()
        },
        // An absent quality array is written as 0xFF fill, which samtools
        // prints as `*` — the same spelling an absent sequence gets.
        qual: if model.qual.is_empty() {
            "*".to_string()
        } else {
            model.qual.iter().map(|q| (q + 33) as char).collect()
        },
        aux: model
            .aux
            .iter()
            .map(|(tag, value)| (String::from_utf8(tag.to_vec()).expect("ASCII"), value.expected()))
            .collect(),
        n_aux: model.aux.len(),
    }
}

// r[verify bam.owned_record.set_seq]
// r[verify bam.owned_record.set_qual]
// r[verify bam.owned_record.flag_methods]
// r[verify bam.owned_record.aux_uniqueness]
// r[verify bam.owned_record.aux_replace_semantics]
// r[verify bam.owned_record.to_bam_bytes]
/// A record driven through a sequence of mutations must, after *every* step,
/// serialize to BAM that samtools reads back as the model says.
#[hegel::test(test_cases = 48)]
fn mutation_sequences_agree_with_the_model_at_every_step(tc: TestCase) {
    let ops = tc.draw(arb_ops().print_as_debug());

    let mut rec = seed_record();
    let mut model = Model { flags: 0x4, ..Model::default() };
    // One snapshot per step: the whole history goes into a single BAM, so one
    // samtools call checks every intermediate state, not just the last.
    let mut snapshots = Vec::with_capacity(ops.len());
    let mut expected = Vec::with_capacity(ops.len());

    for (step, op) in ops.iter().enumerate() {
        apply(op, &mut rec, &mut model);

        // The record's own view has to agree before anyone serializes it.
        assert_eq!(rec.seq, model.seq, "seq after step {step}");
        assert_eq!(
            BaseQuality::slice_to_bytes(&rec.qual),
            model.qual.as_slice(),
            "qual after step {step}"
        );
        assert!(
            rec.qual.is_empty() || rec.qual.len() == rec.seq.len(),
            "a present qual array must be exactly as long as the sequence, step {step}"
        );
        assert_eq!(rec.flags.raw(), model.flags, "flags after step {step}");
        assert_eq!(rec.is_unmapped(), model.flags & 0x4 != 0, "is_unmapped after step {step}");
        assert_eq!(rec.is_reverse(), model.flags & 0x10 != 0, "is_reverse after step {step}");
        for tag in TAGS {
            assert_eq!(
                rec.aux.get(tag).is_some(),
                model.aux.contains_key(&tag),
                "tag {} presence after step {step}",
                String::from_utf8_lossy(&tag)
            );
        }
        assert_eq!(rec.aux.is_empty(), model.aux.is_empty(), "aux emptiness after step {step}");

        let mut snapshot = rec.clone();
        snapshot.qname = format!("s{step}").into_bytes();
        snapshots.push(snapshot);
        expected.push(expected_line(step, &model));
    }

    let dir = tempfile::tempdir().expect("tempdir");
    let path = write_records(dir.path(), &header(SHORT_CONTIG_LEN), &snapshots);

    let got = samtools_lines(&path);
    assert_eq!(got.len(), expected.len(), "record count");
    for (step, (got, want)) in got.iter().zip(&expected).enumerate() {
        assert_eq!(got, want, "step {step}");
    }

    tc.event_value("steps", ops.len() as f64);
    tc.event_value("final seq length", model.seq.len() as f64);
    for op in &ops {
        match op {
            Op::Resize { seq, .. } => {
                tc.event(if seq.len() % 2 == 0 { "resize to even" } else { "resize to odd" });
            }
            Op::Requal { .. } => tc.event("requal"),
            Op::ClearQual => tc.event("clear qual"),
            Op::SetAux { value, .. } => tc.event(value.kind()),
            Op::RemoveAux { .. } => tc.event("remove aux"),
            Op::SetFlags { .. } => tc.event("set flags"),
        }
    }
    let mut seen: BTreeSet<[u8; 2]> = BTreeSet::new();
    for op in &ops {
        if let Op::SetAux { tag, .. } = op
            && !seen.insert(*tag)
        {
            tc.event("a tag was set twice");
        }
    }
}
