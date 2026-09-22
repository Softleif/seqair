//! Tests for performance-related spec rules.
//! after refactoring) rather than measuring wall-clock time.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    reason = "test code"
)]
#![allow(clippy::arithmetic_side_effects, reason = "test code")]
mod helpers;
use helpers::ri;

use hegel::prelude::*;
use helpers::{cigar_op, cigar_ops, make_record, make_record_with_cigar};
use seqair::bam::cigar::{CigarMapping, CigarPosInfo};
use seqair::bam::{Pos0, RecordStore, pileup::PileupEngine};
use seqair_types::QPos;

// ---- perf.reuse_alignment_vec ----
// Verified indirectly: if the engine reuses the vec internally, output must
// still be correct. This test ensures correctness is maintained after the
// optimization.

// r[verify perf.reuse_alignment_vec+3]
#[hegel::test]
fn reused_vec_depth_matches_sweep_line(tc: TestCase) {
    let records = tc.draw(
        gs::vecs(gs::tuples!(
            gs::integers::<i32>().min_value(0).max_value(199),
            gs::integers::<u32>().min_value(20).max_value(80),
        ))
        .min_size(2)
        .max_size(15),
    );
    // Sort by position to match fetch_into's coordinate-sorted guarantee.
    let mut sorted_records = records.clone();
    sorted_records.sort_by_key(|&(offset, _)| offset);

    let mut arena = RecordStore::new();
    // Sweep-line: +1 at read start, -1 at read end+1. Prefix sum = depth.
    const MAX_POS: usize = 301;
    let mut delta = vec![0i32; MAX_POS + 1];
    for &(offset, len) in &sorted_records {
        let start = offset as usize;
        let end_excl = (offset as usize + len as usize).min(MAX_POS);
        if start < MAX_POS {
            delta[start] += 1;
        }
        if end_excl <= MAX_POS {
            delta[end_excl] -= 1;
        }
        arena.push_raw(&make_record(0, offset, 99, 60, len), &mut ()).unwrap();
    }

    let mut expected_depth = vec![0usize; MAX_POS];
    let mut running = 0i32;
    for (i, &d) in delta.iter().enumerate().take(MAX_POS) {
        running += d;
        expected_depth[i] = running.max(0) as usize;
    }

    let mut engine = PileupEngine::new(
        arena.prepare_for_pileup().input,
        Pos0::new(0).unwrap(),
        Pos0::new(300).unwrap(),
    );
    while let Some(col) = engine.pileups() {
        let pos = col.pos().as_usize();
        let exp = expected_depth.get(pos).copied().unwrap_or(0);
        assert_eq!(col.depth(), exp, "depth wrong at pos {}", col.pos());
    }
}

// ---- perf.no_sorted_indices ----
// The engine must produce correct output when records are already sorted
// (which they always are from fetch_into).

// r[verify perf.no_sorted_indices]
#[test]
fn already_sorted_records_produce_correct_pileup() {
    let mut arena = RecordStore::new();
    // Push in sorted order (as fetch_into guarantees)
    arena.push_raw(&make_record(0, 10, 99, 60, 30), &mut ()).unwrap();
    arena.push_raw(&make_record(0, 20, 99, 60, 30), &mut ()).unwrap();
    arena.push_raw(&make_record(0, 30, 99, 60, 30), &mut ()).unwrap();

    let mut engine = PileupEngine::new(
        arena.prepare_for_pileup().input,
        Pos0::new(0).unwrap(),
        Pos0::new(70).unwrap(),
    );
    let columns = helpers::collect_columns(&mut engine);

    // pos 10-19: 1 read
    assert_eq!(columns.iter().find(|c| c.pos() == Pos0::new(15).unwrap()).unwrap().depth(), 1);
    // pos 20-29: 2 reads
    assert_eq!(columns.iter().find(|c| c.pos() == Pos0::new(25).unwrap()).unwrap().depth(), 2);
    // pos 30-39: 3 reads
    assert_eq!(columns.iter().find(|c| c.pos() == Pos0::new(35).unwrap()).unwrap().depth(), 3);
    // pos 40-49: 2 reads (first ended)
    assert_eq!(columns.iter().find(|c| c.pos() == Pos0::new(45).unwrap()).unwrap().depth(), 2);
}

// ---- perf.avoid_redundant_arena_get ----
// Correctness preserved when arena.get() is called fewer times.

// r[verify perf.avoid_redundant_arena_get+2]
#[test]
fn single_arena_get_per_record_entry_still_correct() {
    use seqair::bam::record_store::{CustomizeRecordStore, SlimRecord};

    #[derive(Clone)]
    struct DropSecondary;
    impl CustomizeRecordStore for DropSecondary {
        type Extra = ();
        fn filter(&mut self, rec: &SlimRecord, _: &RecordStore<()>) -> bool {
            !rec.flags.is_secondary()
        }
        fn compute(&mut self, _: &SlimRecord, _: &RecordStore<()>) {}
    }

    let mut arena = RecordStore::new();
    arena.push_raw(&make_record(0, 0, 99, 60, 50), &mut DropSecondary).unwrap();
    arena.push_raw(&make_record(0, 10, 99 | 0x100, 40, 50), &mut DropSecondary).unwrap(); // secondary flag

    let mut engine = PileupEngine::new(
        arena.prepare_for_pileup().input,
        Pos0::new(0).unwrap(),
        Pos0::new(59).unwrap(),
    );
    let columns = helpers::collect_columns(&mut engine);

    // Only the mapq=60 read should remain in the store after push-time filtering.
    for col in &columns {
        if col.pos() < Pos0::new(10).unwrap() {
            assert_eq!(col.depth(), 1, "only high-mapq read at pos {}", col.pos());
        }
    }
    // At pos 10+, still only 1 because second read was filtered before push.
    let col = columns.iter().find(|c| c.pos() == Pos0::new(25).unwrap()).unwrap();
    assert_eq!(col.depth(), 1);
}

// ---- perf.cigar_no_to_vec ----
// After removing .to_vec(), CigarMapping must still produce correct qpos.

// r[verify perf.cigar_no_to_vec]
#[test]
fn cigar_index_from_arena_slab_correct() {
    let mut arena = RecordStore::new();
    // 30M 5D 20M — tests that cigar bytes from arena slab work without copying
    let raw = make_record_with_cigar(
        0,
        100,
        99,
        60,
        &[cigar_op(30, 0), cigar_op(5, 2), cigar_op(20, 0)],
        50,
    );
    arena.push_raw(&raw, &mut ()).unwrap();

    let mut engine = PileupEngine::new(
        arena.prepare_for_pileup().input,
        Pos0::new(100).unwrap(),
        Pos0::new(154).unwrap(),
    );
    let columns = helpers::collect_columns(&mut engine);

    // Before deletion: qpos = pos - 100
    let col = columns.iter().find(|c| c.pos() == Pos0::new(110).unwrap()).unwrap();
    assert_eq!(col.alignments().next().unwrap().qpos(), Some(QPos::new(10)));

    // Inside deletion: alignment present with Deletion op (no qpos)
    let del_col = columns.iter().find(|c| c.pos() == Pos0::new(132).unwrap()).unwrap();
    assert_eq!(del_col.depth(), 1);
    assert!(del_col.alignments().next().unwrap().is_del());

    // After deletion: qpos = pos - 100 - 5 (deletion consumes 5 ref but 0 query)
    let col = columns.iter().find(|c| c.pos() == Pos0::new(140).unwrap()).unwrap();
    assert_eq!(col.alignments().next().unwrap().qpos(), Some(QPos::new(35)));
}

// ---- perf.precompute_matches_indels ----

// r[verify perf.precompute_matches_indels]
#[test]
fn precomputed_matches_indels_accessible_from_record() {
    let mut arena = RecordStore::new();
    // 30M 5I 15M = 45 matches, 5 indels
    let raw = make_record_with_cigar(
        0,
        0,
        99,
        60,
        &[cigar_op(30, 0), cigar_op(5, 1), cigar_op(15, 0)],
        50,
    );
    arena.push_raw(&raw, &mut ()).unwrap();

    let r = arena.record(ri(0)).unwrap();
    assert_eq!(r.matching_bases, 45);
    assert_eq!(r.indel_bases, 5);
}

// r[verify perf.precompute_matches_indels]
#[hegel::test]
fn matches_indels_from_varied_cigar_ops(tc: TestCase) {
    // Generate a sequence of (op_type, length) pairs from {M=0, I=1, D=2, S=4}.
    // S is only valid at the ends, so we place it as leading/trailing clips.
    let clip = || gs::integers::<u32>().max_value(10);
    let leading_clip = tc.draw(clip());
    let trailing_clip = tc.draw(clip());
    let inner_ops = tc.draw(
        gs::vecs(gs::tuples!(
            gs::sampled_from(&[0u8, 1u8, 2u8]),
            gs::integers::<u32>().min_value(1).max_value(50),
        ))
        .min_size(1)
        .max_size(8),
    );
    // Build CIGAR: optional leading S, inner ops {M,I,D}, optional trailing S.
    // Per r[cigar.matches_indels]: matching_bases = sum of M(0) lengths only,
    // indel_bases = sum of I(1) + D(2) lengths.
    let mut ops: Vec<u32> = Vec::new();
    let mut seq_len = 0u32;

    if leading_clip > 0 {
        ops.push(cigar_op(leading_clip, 4)); // S
        seq_len += leading_clip;
    }

    // Ensure at least one M op so the record has a ref-consuming op.
    let mut has_m = false;
    let mut exp_matches = 0u32;
    let mut exp_indels = 0u32;

    for &(op_type, len) in &inner_ops {
        ops.push(cigar_op(len, op_type));
        match op_type {
            0 => {
                exp_matches += len;
                seq_len += len;
                has_m = true;
            }
            1 => {
                exp_indels += len;
                seq_len += len;
            }
            2 => {
                exp_indels += len;
            } // D consumes ref, not query
            _ => {}
        }
    }

    // Guarantee at least one M op so the record is usable.
    if !has_m {
        ops.push(cigar_op(1, 0)); // M
        seq_len += 1;
        exp_matches += 1;
    }

    if trailing_clip > 0 {
        ops.push(cigar_op(trailing_clip, 4)); // S
        seq_len += trailing_clip;
    }

    tc.assume(seq_len > 0);

    let raw = make_record_with_cigar(0, 0, 99, 60, &ops, seq_len);
    let mut arena = RecordStore::new();
    arena.push_raw(&raw, &mut ()).unwrap();
    let r = arena.record(ri(0)).unwrap();

    assert_eq!(r.matching_bases, exp_matches, "matching_bases mismatch for ops={:?}", inner_ops);
    assert_eq!(r.indel_bases, exp_indels, "indel_bases mismatch for ops={:?}", inner_ops);
}

// ---- perf.cigar_binary_search ----

// r[verify perf.cigar_binary_search]
#[test]
fn binary_search_correct_for_many_ops() {
    // Simulate RNA-seq: 30M 5000N 30M 3000N 40M (5 ops, triggers binary search)
    let ops = cigar_ops(&[
        cigar_op(30, 0),
        cigar_op(5000, 3),
        cigar_op(30, 0),
        cigar_op(3000, 3),
        cigar_op(40, 0),
    ]);
    let mapping = CigarMapping::new(Pos0::new(1000).unwrap(), &ops).unwrap();

    // First M block: 1000-1029
    assert_eq!(
        mapping.pos_info_at(Pos0::new(1000).unwrap()),
        Some(CigarPosInfo::Match { qpos: QPos::new(0) })
    );
    assert_eq!(
        mapping.pos_info_at(Pos0::new(1029).unwrap()),
        Some(CigarPosInfo::Match { qpos: QPos::new(29) })
    );

    // N skip: 1030-6029 → RefSkip
    assert_eq!(mapping.pos_info_at(Pos0::new(1030).unwrap()), Some(CigarPosInfo::RefSkip));
    assert_eq!(mapping.pos_info_at(Pos0::new(5000).unwrap()), Some(CigarPosInfo::RefSkip));

    // Second M block: 6030-6059
    assert_eq!(
        mapping.pos_info_at(Pos0::new(6030).unwrap()),
        Some(CigarPosInfo::Match { qpos: QPos::new(30) })
    );
    assert_eq!(
        mapping.pos_info_at(Pos0::new(6059).unwrap()),
        Some(CigarPosInfo::Match { qpos: QPos::new(59) })
    );

    // Second N skip: 6060-9059 → RefSkip
    assert_eq!(mapping.pos_info_at(Pos0::new(6060).unwrap()), Some(CigarPosInfo::RefSkip));

    // Third M block: 9060-9099
    assert_eq!(
        mapping.pos_info_at(Pos0::new(9060).unwrap()),
        Some(CigarPosInfo::Match { qpos: QPos::new(60) })
    );
    assert_eq!(
        mapping.pos_info_at(Pos0::new(9099).unwrap()),
        Some(CigarPosInfo::Match { qpos: QPos::new(99) })
    );
    assert_eq!(mapping.pos_info_at(Pos0::new(9100).unwrap()), None);
}

// r[verify perf.cigar_binary_search]
#[hegel::test]
fn binary_search_correct_for_rna_seq_pattern(tc: TestCase) {
    let n_exons = tc.draw(gs::integers::<usize>().min_value(3).max_value(8));
    let exon_len = tc.draw(gs::integers::<u32>().min_value(20).max_value(100));
    let intron_len = tc.draw(gs::integers::<u32>().min_value(100).max_value(5000));
    // Build alternating M/N CIGAR (RNA-seq pattern, triggers bsearch path for ≥5 ops)
    let mut ops = Vec::new();
    let mut total_query = 0u32;
    for i in 0..n_exons {
        ops.push(cigar_op(exon_len, 0)); // M
        total_query += exon_len;
        if i < n_exons - 1 {
            ops.push(cigar_op(intron_len, 3)); // N
        }
    }
    let typed = cigar_ops(&ops);
    let mapping = CigarMapping::new(Pos0::new(0).unwrap(), &typed).unwrap();

    // Compute total ref span
    let total_ref: u32 = ops
        .iter()
        .map(|&op| {
            let len = op >> 4;
            let op_type = (op & 0xF) as u8;
            if matches!(op_type, 0 | 2 | 3 | 7 | 8) { len } else { 0 }
        })
        .sum();

    // Check every position — qpos must be in range and monotonically increasing
    let mut last_qpos: Option<QPos> = None;
    for pos in 0..total_ref {
        match mapping.pos_info_at(Pos0::new(pos).unwrap()) {
            Some(CigarPosInfo::Match { qpos }) => {
                assert!(
                    qpos.as_usize() < total_query as usize,
                    "qpos {qpos} out of range at ref {pos}"
                );
                if let Some(prev) = last_qpos {
                    assert!(qpos > prev, "qpos not monotonic at ref {pos}");
                }
                last_qpos = Some(qpos);
            }
            Some(CigarPosInfo::RefSkip) => {}
            other => panic!("unexpected result at ref {pos}: {other:?}"),
        }
    }
    // Out of range
    assert_eq!(mapping.pos_info_at(Pos0::new(total_ref).unwrap()), None);
}

// ---- perf.arena_capacity_hint ----

// r[verify perf.arena_capacity_hint+2]
#[test]
fn arena_with_capacity_avoids_realloc() {
    // Pre-size based on estimated compressed bytes
    let mut arena = RecordStore::with_byte_hint(50_000);

    for i in 0..100 {
        arena.push_raw(&make_record(0, i * 10, 99, 60, 50), &mut ()).unwrap();
    }
    assert_eq!(arena.len(), 100);

    // All records accessible
    for (n, r) in arena.records().enumerate() {
        assert_eq!(r.pos, Pos0::new(u32::try_from(n).unwrap() * 10).unwrap());
    }
}

// ---- bam.reader.early_exit ----

/// A narrow query must not walk the rest of the contig.
///
/// The reader used to treat a record past the query end as a skip and keep
/// reading, so every query ran to the end of the last chunk the index handed
/// back — and a BAI query hands back the chunks of every ancestor bin, whose
/// byte ranges run far past the region. The observable consequence is the
/// record count the query touches: with the early exit it is the region's own
/// reads plus the tail of the block it starts in, without it the whole contig.
// r[verify bam.reader.early_exit]
#[test]
fn a_narrow_query_stops_at_the_first_record_past_its_end() {
    use seqair::bam::IndexedBamReader;

    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.bam");
    let mut reader = IndexedBamReader::open(std::path::Path::new(path)).unwrap();
    let tid = reader.header().tid("chr19").unwrap();

    // The test BAM's reads cover chr19:6,103,075–6,143,000 (0-based).
    let mut examined = |start: u32, end: u32| {
        let mut query =
            reader.query(tid, Pos0::new(start).unwrap(), Pos0::new(end).unwrap()).unwrap();
        query.for_each(|_| {}).unwrap();
        let c = query.counts();
        (c.fetched, c.fetched + c.skipped_out_of_range + c.skipped_tid)
    };

    let (whole_fetched, whole_examined) = examined(0, 6_143_000);
    let (narrow_fetched, narrow_examined) = examined(6_103_000, 6_103_999);

    assert!(narrow_fetched > 0, "the narrow region must contain reads at all");
    assert!(narrow_fetched < whole_fetched, "the narrow region must be a strict subset");
    // Everything examined and not kept started before the region — the reader
    // may still be inside the block the region begins in. Nothing past the end
    // should be reached at all: skip-and-continue examined 2529 records here
    // for the same 313 it kept, on a BAM of one small contig.
    assert!(
        narrow_examined < narrow_fetched * 2,
        "a 1 kb query examined {narrow_examined} records to keep {narrow_fetched} \
         (the contig has {whole_examined}) — it is reading past its end"
    );
}
