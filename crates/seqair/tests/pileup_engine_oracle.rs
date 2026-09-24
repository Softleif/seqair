//! The pileup engine against a naive per-position oracle.
//!
//! The engine keeps per-read state across columns and walks each read's CIGAR
//! incrementally. The oracle keeps nothing: for every position it asks every
//! record, from scratch, what the stateless [`CigarMapping`] lookups say is
//! there. Any column the two disagree on — which records, in what order, what
//! op, which anchored indel, which per-read fields — is a bug in the engine's
//! bookkeeping.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    reason = "test code"
)]

use hegel::prelude::*;
use seqair::bam::{
    Pos0, RecordStore,
    cigar::{
        CigarMapping, CigarOp, CigarOpType, CigarPosInfo, calc_matches_indels, compute_end_pos,
    },
    pileup::{Indel, PileupEngine, PileupOp},
    record_store::PileupInput,
};
use seqair_types::{BamFlags, Base, BaseQuality, Offset};
use std::num::NonZeroU32;

const BASES: [Base; 4] = [Base::A, Base::C, Base::G, Base::T];
const OPS: [CigarOpType; 9] = [
    CigarOpType::Match,
    CigarOpType::Insertion,
    CigarOpType::Deletion,
    CigarOpType::RefSkip,
    CigarOpType::SoftClip,
    CigarOpType::HardClip,
    CigarOpType::Padding,
    CigarOpType::SeqMatch,
    CigarOpType::SeqMismatch,
];

#[derive(Clone, Debug)]
struct Read {
    pos: u32,
    ops: Vec<CigarOp>,
    unmapped: bool,
    reverse: bool,
    mapq: u8,
    seq: Vec<Base>,
    qual: Vec<u8>,
}

#[hegel::composite]
fn arb_read(tc: &TestCase) -> Read {
    let pos = tc.draw_silent(gs::integers::<u32>().max_value(80));
    let n_ops = tc.draw_silent(gs::integers::<usize>().min_value(1).max_value(9));
    let mut ops = Vec::with_capacity(n_ops);
    for _ in 0..n_ops {
        let op = tc.draw_silent(gs::sampled_from(&OPS));
        // Zero-length ops are legal BAM and are where off-by-one bugs live.
        let len = tc.draw_silent(gs::integers::<u32>().max_value(7));
        ops.push(CigarOp::new(op, len));
    }
    let qlen = seqair::bam::cigar::calc_query_len(&ops) as usize;
    let seq_star = tc.draw_silent(gs::integers::<u8>().max_value(9)) == 0;
    let (seq, qual) = if seq_star {
        (Vec::new(), Vec::new())
    } else {
        let seq = (0..qlen).map(|_| tc.draw_silent(gs::sampled_from(&BASES))).collect();
        let qual = (0..qlen).map(|_| tc.draw_silent(gs::integers::<u8>().max_value(60))).collect();
        (seq, qual)
    };
    Read {
        pos,
        ops,
        unmapped: tc.draw_silent(gs::integers::<u8>().max_value(15)) == 0,
        reverse: tc.draw_silent(gs::booleans()),
        mapq: tc.draw_silent(gs::integers::<u8>().max_value(60)),
        seq,
        qual,
    }
}

#[derive(Clone, Debug)]
struct Case {
    reads: Vec<Read>,
    /// Pairs `(i, i + 1)` for every `i` listed become mates.
    pair_starts: Vec<usize>,
    overhang: u32,
    max_depth: Option<u32>,
    region: (u32, u32),
}

#[hegel::composite]
fn arb_case(tc: &TestCase) -> Case {
    let mut reads = tc.draw_silent(gs::vecs(arb_read()).min_size(1).max_size(14));
    reads.sort_by_key(|r| r.pos);
    let mut pair_starts = Vec::new();
    let mut i = 0;
    while i + 1 < reads.len() {
        if tc.draw_silent(gs::booleans()) {
            pair_starts.push(i);
            i += 2;
        } else {
            i += 1;
        }
    }
    let overhang = tc.draw_silent(gs::sampled_from(&[0u32, 0, 1, 3]));
    let max_depth = tc.draw_silent(gs::sampled_from(&[None, None, Some(1u32), Some(3)]));
    let start = tc.draw_silent(gs::integers::<u32>().max_value(40));
    let last = tc.draw_silent(gs::integers::<u32>().min_value(start).max_value(140));
    Case { reads, pair_starts, overhang, max_depth, region: (start, last) }
}

fn build_store(case: &Case) -> PileupInput<()> {
    let mut store = RecordStore::<()>::new();
    let mut mate_of: Vec<Option<usize>> = vec![None; case.reads.len()];
    for &i in &case.pair_starts {
        mate_of[i] = Some(i + 1);
        mate_of[i + 1] = Some(i);
    }
    for (i, r) in case.reads.iter().enumerate() {
        let p = Pos0::new(r.pos).unwrap();
        let end_pos = if r.unmapped { p } else { compute_end_pos(p, &r.ops).unwrap() };
        let (matching, indels) = calc_matches_indels(&r.ops);
        let mut flags = 0u16;
        if r.unmapped {
            flags |= 0x4;
        }
        if r.reverse {
            flags |= 0x10;
        }
        let (name, next_ref_id, next_pos) = match mate_of[i] {
            Some(m) => {
                flags |= 0x1 | if m > i { 0x40 } else { 0x80 };
                (format!("t{}", i.min(m)), 0, i32::try_from(case.reads[m].pos).unwrap())
            }
            None => (format!("r{i}"), -1, -1),
        };
        store
            .push_fields(
                p,
                end_pos,
                BamFlags::from(flags),
                r.mapq,
                matching,
                indels,
                name.as_bytes(),
                &r.ops,
                &r.seq,
                &r.qual,
                &[],
                0,
                next_ref_id,
                next_pos,
                0,
                &mut (),
            )
            .unwrap();
    }
    store.prepare_for_pileup().input
}

/// One column entry, with every field a consumer can read off it.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Entry {
    record_idx: u32,
    op: PileupOp,
    indel_after: Indel,
    in_mate_overlap: bool,
    mate_idx: Option<u32>,
    mapq: u8,
    flags: u16,
    seq_len: u32,
    matching_bases: u32,
    indel_bases: u32,
}

fn run_engine(case: &Case) -> Vec<(u32, Vec<Entry>)> {
    let input = build_store(case);
    let region = Pos0::new(case.region.0).unwrap()..=Pos0::new(case.region.1).unwrap();
    let mut engine = PileupEngine::new(input, region.into());
    engine.set_soft_clip_overhang(case.overhang);
    if let Some(d) = case.max_depth {
        engine.set_max_depth(NonZeroU32::new(d).unwrap());
    }
    let mut out = Vec::new();
    while let Some(col) = engine.pileups() {
        let entries = col
            .alignments()
            .map(|v| Entry {
                record_idx: v.record_idx().get(),
                op: v.op,
                indel_after: v.indel_after(),
                in_mate_overlap: v.in_mate_overlap(),
                mate_idx: v.mate_idx().map(|m| m.get()),
                mapq: v.mapq,
                flags: v.flags.raw(),
                seq_len: v.seq_len,
                matching_bases: v.matching_bases,
                indel_bases: v.indel_bases,
            })
            .collect();
        out.push((col.pos().as_u32(), entries));
    }
    out
}

fn base_qual(seq: &[Base], qual: &[BaseQuality], qpos: u32) -> (Base, BaseQuality) {
    match (seq.get(qpos as usize), qual.get(qpos as usize)) {
        (Some(&b), Some(&q)) => (b, q),
        _ => (Base::Unknown, BaseQuality::UNAVAILABLE),
    }
}

/// The same columns, derived per position with no state carried between them.
fn run_oracle(case: &Case) -> Vec<(u32, Vec<Entry>)> {
    let input = build_store(case);
    let store = input.store();
    let oh = i64::from(case.overhang);
    let mut out = Vec::new();
    for p in case.region.0..=case.region.1 {
        let pos = Pos0::new(p).unwrap();
        let mut entries = Vec::new();
        for idx in store.indices() {
            let rec = store.record(idx).unwrap();
            if rec.flags.is_unmapped() {
                continue;
            }
            let (first, last) = (rec.pos.as_i64() - oh, rec.end_pos.as_i64() + oh);
            if i64::from(p) < first || i64::from(p) > last {
                continue;
            }
            let mapping = CigarMapping::new(rec.pos, rec.cigar()).unwrap();
            let (seq, qual) = (rec.seq(), rec.qual());
            let (op, indel_after) = match mapping.pos_info_at(pos) {
                Some(CigarPosInfo::Match { qpos }) => {
                    let (base, qual) = base_qual(seq, qual, qpos.get());
                    let indel = match mapping.deletion_after_at(pos) {
                        Some(n) => Indel::Deletion(n),
                        None => Indel::None,
                    };
                    (PileupOp::Match { qpos, base, qual }, indel)
                }
                Some(CigarPosInfo::Insertion { qpos, insert_len }) => {
                    let (base, qual) = base_qual(seq, qual, qpos.get());
                    (
                        PileupOp::Insertion { qpos, base, qual, insert_len },
                        Indel::Insertion(insert_len),
                    )
                }
                Some(CigarPosInfo::Deletion { del_len }) => {
                    (PileupOp::Deletion { del_len }, Indel::None)
                }
                Some(CigarPosInfo::ComplexIndel { del_len, insert_len, is_refskip }) => {
                    (PileupOp::ComplexIndel { del_len, insert_len, is_refskip }, Indel::None)
                }
                Some(CigarPosInfo::RefSkip) => (PileupOp::RefSkip, Indel::None),
                None => {
                    let Some(qpos) = mapping.soft_clip_qpos_at(pos, case.overhang, rec.seq_len)
                    else {
                        continue;
                    };
                    let (base, qual) = base_qual(seq, qual, qpos.get());
                    (PileupOp::SoftClip { qpos, base, qual }, Indel::None)
                }
            };
            let in_mate_overlap = rec.mate_overlap().is_some_and(|ov| {
                let o = Offset::new(oh);
                let start = ov.start.checked_sub_offset(o).unwrap_or(Pos0::ZERO);
                let end = ov.end.checked_add_offset(o).unwrap_or(Pos0::MAX);
                (start..end).contains(&pos)
            });
            entries.push(Entry {
                record_idx: idx.get(),
                op,
                indel_after,
                in_mate_overlap,
                mate_idx: rec.mate_idx().map(|m| m.get()),
                mapq: rec.mapq,
                flags: rec.flags.raw(),
                seq_len: rec.seq_len,
                matching_bases: rec.matching_bases,
                indel_bases: rec.indel_bases,
            });
        }
        if let Some(d) = case.max_depth {
            entries.truncate(d as usize);
        }
        if !entries.is_empty() {
            out.push((p, entries));
        }
    }
    out
}

// r[verify pileup.column_contents]
// r[verify pileup.column_record_order]
// r[verify pileup.mate_link_cache]
// r[verify pileup.soft_clip_overhang.emit]
// r[verify pileup_indel.complex_indel]
// r[verify pileup_indel.insertion_len]
/// Every column the engine yields is the column a stateless per-position
/// lookup builds, for arbitrary CIGARs (zero-length ops, padding, `D I`,
/// `N I`, hard and soft clips), `SEQ=*`, unmapped reads, linked mates,
/// soft-clip overhang, and depth caps.
#[hegel::test]
fn engine_matches_stateless_oracle(tc: TestCase) {
    let case = tc.draw(arb_case().print_as_debug());
    let ours = run_engine(&case);
    let oracle = run_oracle(&case);
    assert_eq!(ours, oracle);
    tc.event_value("columns", ours.len() as f64);
}
