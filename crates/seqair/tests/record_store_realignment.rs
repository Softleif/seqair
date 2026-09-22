//! Tests for `RecordStore` alignment mutation (cigar slab, `set_alignment`, tid, mate fields).
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::too_many_arguments,
    clippy::type_complexity,
    reason = "test code"
)]
#![allow(
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    reason = "test code with known small values"
)]

mod helpers;
use hegel::prelude::*;
use helpers::ri;
use seqair::bam::{DecodeError, Pos0};
use seqair::bam::record_store::RecordStore;
use seqair_types::{BamFlags, Base, BaseQuality};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// BAM CIGAR op type constants.
const CIGAR_M: u8 = 0;
const CIGAR_I: u8 = 1;
const CIGAR_D: u8 = 2;
const CIGAR_S: u8 = 4;

/// Encode a CIGAR op as a packed BAM u32 (len << 4 | op).
fn pack_cigar_op(len: u32, op: u8) -> u32 {
    (len << 4) | u32::from(op)
}

/// Encode a list of (len, op) pairs as typed `CigarOp` slots.
fn pack_cigar(ops: &[(u32, u8)]) -> Vec<seqair::bam::CigarOp> {
    ops.iter()
        .map(|&(len, op)| seqair::bam::CigarOp::from_bam_u32(pack_cigar_op(len, op)))
        .collect()
}

/// Compute query-consuming length from a list of (len, op) parts.
/// Independent oracle — does NOT call any seqair code.
fn query_len_from_parts(ops: &[(u32, u8)]) -> u32 {
    ops.iter()
        .filter(|(_, op)| matches!(*op, CIGAR_M | CIGAR_I | CIGAR_S | 7 | 8)) // M/I/S/=/X
        .map(|(len, _)| len)
        .sum()
}

/// Compute reference-consuming length from a list of (len, op) parts.
/// Independent oracle — does NOT call any seqair code.
fn ref_len_from_parts(ops: &[(u32, u8)]) -> u32 {
    ops.iter()
        .filter(|(_, op)| matches!(*op, CIGAR_M | CIGAR_D | 3 | 7 | 8)) // M/D/N/=/X
        .map(|(len, _)| *len)
        .sum()
}

/// Compute expected `end_pos` from pos and cigar parts (independent oracle).
fn expected_end_pos(pos: u32, ops: &[(u32, u8)]) -> u32 {
    let rlen = ref_len_from_parts(ops);
    if rlen == 0 { pos } else { pos + rlen - 1 }
}

/// Compute `matching_bases` from parts (independent oracle).
fn expected_matching_bases(ops: &[(u32, u8)]) -> u32 {
    ops.iter().filter(|(_, op)| *op == CIGAR_M).map(|(len, _)| len).sum()
}

/// Compute `indel_bases` from parts (independent oracle).
fn expected_indel_bases(ops: &[(u32, u8)]) -> u32 {
    ops.iter().filter(|(_, op)| *op == CIGAR_I || *op == CIGAR_D).map(|(len, _)| len).sum()
}

/// Build a minimal raw BAM record for testing.
/// Same helper as in `record_store.rs` tests, extended with `next_pos` and `template_len`.
fn make_test_record(
    tid: i32,
    pos: i32,
    flags: u16,
    mapq: u8,
    qname: &[u8],
    cigar_ops: &[u32],
    packed_seq: &[u8],
    qual: &[u8],
    aux: &[u8],
    next_ref_id: i32,
    next_pos: i32,
    template_len: i32,
) -> Vec<u8> {
    let name_len = qname.len() + 1;
    let n_cigar_ops = cigar_ops.len();
    let cigar_bytes = n_cigar_ops * 4;
    let seq_len = qual.len() as u32;
    let seq_bytes = packed_seq.len();

    let total = 32 + name_len + cigar_bytes + seq_bytes + qual.len() + aux.len();
    let mut raw = vec![0u8; total];

    raw[0..4].copy_from_slice(&tid.to_le_bytes());
    raw[4..8].copy_from_slice(&pos.to_le_bytes());
    raw[8] = name_len as u8;
    raw[9] = mapq;
    // bytes 10..12: bin (unused here)
    raw[12..14].copy_from_slice(&(n_cigar_ops as u16).to_le_bytes());
    raw[14..16].copy_from_slice(&flags.to_le_bytes());
    raw[16..20].copy_from_slice(&seq_len.to_le_bytes());
    raw[20..24].copy_from_slice(&next_ref_id.to_le_bytes());
    raw[24..28].copy_from_slice(&next_pos.to_le_bytes());
    raw[28..32].copy_from_slice(&template_len.to_le_bytes());

    let mut off = 32;
    raw[off..off + qname.len()].copy_from_slice(qname);
    raw[off + qname.len()] = 0;
    off += name_len;

    for &op in cigar_ops {
        raw[off..off + 4].copy_from_slice(&op.to_le_bytes());
        off += 4;
    }

    raw[off..off + seq_bytes].copy_from_slice(packed_seq);
    off += seq_bytes;

    raw[off..off + qual.len()].copy_from_slice(qual);
    off += qual.len();

    raw[off..off + aux.len()].copy_from_slice(aux);

    raw
}

/// Shorthand: make a record with simple defaults for mate fields.
fn make_simple_record(
    tid: i32,
    pos: i32,
    cigar_ops: &[u32],
    seq_len: u32,
    qname: &[u8],
) -> Vec<u8> {
    let packed_seq: Vec<u8> = (0..seq_len.div_ceil(2)).map(|_| 0x12).collect(); // A/C pattern
    let qual: Vec<u8> = vec![30; seq_len as usize];
    make_test_record(tid, pos, 0x63, 60, qname, cigar_ops, &packed_seq, &qual, &[], 0, 200, 150)
}

// ---------------------------------------------------------------------------
// Unit tests: new fields (tid, next_pos, template_len)
// ---------------------------------------------------------------------------

// r[verify record_store.slim_record_fields]
#[test]
fn push_raw_preserves_tid() {
    let raw = make_simple_record(5, 100, &[pack_cigar_op(4, CIGAR_M)], 4, b"read1");
    let mut store = RecordStore::new();
    let idx = store.push_raw(&raw, &mut ()).unwrap().expect("kept");
    assert_eq!(store.record(idx).unwrap().tid, 5);
}

// r[verify record_store.slim_record_fields]
#[test]
fn push_raw_preserves_next_pos() {
    let raw = make_test_record(
        0,
        100,
        0x63,
        60,
        b"r1",
        &[pack_cigar_op(4, CIGAR_M)],
        &[0x12, 0x48],
        &[30; 4],
        &[],
        0,
        500,
        300,
    );
    let mut store = RecordStore::new();
    let idx = store.push_raw(&raw, &mut ()).unwrap().expect("kept");
    assert_eq!(store.record(idx).unwrap().next_pos, 500);
    assert_eq!(store.record(idx).unwrap().template_len, 300);
}

// r[verify record_store.slim_record_fields]
#[test]
fn push_raw_preserves_negative_template_len() {
    let raw = make_test_record(
        0,
        500,
        0x63,
        60,
        b"r1",
        &[pack_cigar_op(4, CIGAR_M)],
        &[0x12, 0x48],
        &[30; 4],
        &[],
        0,
        200,
        -300,
    );
    let mut store = RecordStore::new();
    let idx = store.push_raw(&raw, &mut ()).unwrap().expect("kept");
    assert_eq!(store.record(idx).unwrap().template_len, -300);
}

// ---------------------------------------------------------------------------
// Unit tests: cigar slab separation
// ---------------------------------------------------------------------------

#[test]
fn cigar_in_own_slab() {
    let cigar_packed = pack_cigar(&[(4, CIGAR_M)]);
    let raw = make_simple_record(0, 100, &[pack_cigar_op(4, CIGAR_M)], 4, b"read1");
    let mut store = RecordStore::new();
    let idx = store.push_raw(&raw, &mut ()).unwrap().expect("kept");

    // cigar accessor should return the correct packed bytes
    assert_eq!(store.record(idx).unwrap().cigar(), &cigar_packed);
    // qual and aux should be unaffected by cigar separation
    assert_eq!(store.record(idx).unwrap().qual().len(), 4);
    assert!(store.record(idx).unwrap().qual().iter().all(|q| q.get() == Some(30)));
}

// ---------------------------------------------------------------------------
// Unit tests: set_alignment
// ---------------------------------------------------------------------------

// r[verify record_store.set_alignment]
#[test]
fn set_alignment_updates_cigar_and_pos() {
    // Start with 4M at pos 100
    let raw = make_simple_record(0, 100, &[pack_cigar_op(4, CIGAR_M)], 4, b"read1");
    let mut store = RecordStore::new();
    let idx = store.push_raw(&raw, &mut ()).unwrap().expect("kept");

    // Realign to pos 95 with 2S2M (still query length 4)
    let new_cigar = pack_cigar(&[(2, CIGAR_S), (2, CIGAR_M)]);
    store.set_alignment(idx, Pos0::new(95).unwrap(), &new_cigar).unwrap();

    let rec = store.record(idx).unwrap();
    assert_eq!(rec.pos, Pos0::new(95).unwrap());
    assert_eq!(rec.end_pos, Pos0::new(96).unwrap()); // 95 + 2M - 1
    assert_eq!(rec.n_cigar_ops, 2);
    assert_eq!(rec.matching_bases, 2); // only 2M now
    assert_eq!(rec.indel_bases, 0);
    assert_eq!(store.record(idx).unwrap().cigar(), &new_cigar);
}

// r[verify record_store.set_alignment]
#[test]
fn set_alignment_preserves_other_fields() {
    let raw = make_test_record(
        3,
        100,
        0x63,
        60,
        b"read1",
        &[pack_cigar_op(4, CIGAR_M)],
        &[0x12, 0x48],
        &[30; 4],
        b"RGZgrp\0",
        0,
        500,
        300,
    );
    let mut store = RecordStore::new();
    let idx = store.push_raw(&raw, &mut ()).unwrap().expect("kept");

    let orig_seq: Vec<Base> = store.record(idx).unwrap().seq().to_vec();
    let orig_qual: Vec<BaseQuality> = store.record(idx).unwrap().qual().to_vec();
    let orig_aux: Vec<u8> = store.record(idx).unwrap().aux().to_vec();
    let orig_qname: Vec<u8> = store.record(idx).unwrap().qname().to_vec();

    let new_cigar = pack_cigar(&[(2, CIGAR_S), (2, CIGAR_M)]);
    store.set_alignment(idx, Pos0::new(95).unwrap(), &new_cigar).unwrap();

    // These must be unchanged
    assert_eq!(store.record(idx).unwrap().seq(), &orig_seq);
    assert_eq!(store.record(idx).unwrap().qual(), &orig_qual);
    assert_eq!(store.record(idx).unwrap().aux(), &orig_aux);
    assert_eq!(store.record(idx).unwrap().qname(), &orig_qname);
    assert_eq!(store.record(idx).unwrap().tid, 3);
    assert_eq!(store.record(idx).unwrap().next_pos, 500);
    assert_eq!(store.record(idx).unwrap().template_len, 300);
    assert_eq!(store.record(idx).unwrap().flags, BamFlags::from(0x63));
    assert_eq!(store.record(idx).unwrap().mapq, 60);
}

// r[verify record_store.set_alignment]
#[test]
fn set_alignment_does_not_corrupt_other_records() {
    let raw1 = make_simple_record(0, 100, &[pack_cigar_op(4, CIGAR_M)], 4, b"read1");
    let raw2 = make_simple_record(0, 200, &[pack_cigar_op(4, CIGAR_M)], 4, b"read2");
    let mut store = RecordStore::new();
    let idx1 = store.push_raw(&raw1, &mut ()).unwrap().expect("kept");
    let idx2 = store.push_raw(&raw2, &mut ()).unwrap().expect("kept");

    let orig_cigar2 = store.record(idx2).unwrap().cigar().to_vec();
    let orig_qual2 = store.record(idx2).unwrap().qual().to_vec();

    // Realign only record 1
    let new_cigar = pack_cigar(&[(2, CIGAR_S), (2, CIGAR_M)]);
    store.set_alignment(idx1, Pos0::new(95).unwrap(), &new_cigar).unwrap();

    // Record 2 must be untouched
    assert_eq!(store.record(idx2).unwrap().pos, Pos0::new(200).unwrap());
    assert_eq!(store.record(idx2).unwrap().cigar(), &orig_cigar2);
    assert_eq!(store.record(idx2).unwrap().qual(), &orig_qual2);
}

// r[verify record_store.set_alignment.validation]
#[test]
fn set_alignment_rejects_query_length_mismatch() {
    let raw = make_simple_record(0, 100, &[pack_cigar_op(4, CIGAR_M)], 4, b"read1");
    let mut store = RecordStore::new();
    let idx = store.push_raw(&raw, &mut ()).unwrap().expect("kept");

    // New cigar has query length 6 but record has seq_len 4
    let bad_cigar = pack_cigar(&[(6, CIGAR_M)]);
    let result = store.set_alignment(idx, Pos0::new(95).unwrap(), &bad_cigar);
    assert!(result.is_err(), "should reject query length mismatch");

    // Record should be unchanged
    assert_eq!(store.record(idx).unwrap().pos, Pos0::new(100).unwrap());
}

// r[verify record_store.set_alignment]
#[test]
fn set_alignment_with_indels() {
    // seq_len = 6: 3M1I2M (query = 3+1+2 = 6, ref = 3+2 = 5)
    let raw = make_simple_record(0, 100, &[pack_cigar_op(6, CIGAR_M)], 6, b"read1");
    let mut store = RecordStore::new();
    let idx = store.push_raw(&raw, &mut ()).unwrap().expect("kept");

    let new_cigar = pack_cigar(&[(3, CIGAR_M), (1, CIGAR_I), (2, CIGAR_M)]);
    store.set_alignment(idx, Pos0::new(100).unwrap(), &new_cigar).unwrap();

    let rec = store.record(idx).unwrap();
    assert_eq!(rec.end_pos, Pos0::new(104).unwrap()); // 100 + 5 - 1
    assert_eq!(rec.matching_bases, 5); // 3M + 2M
    assert_eq!(rec.indel_bases, 1); // 1I
    assert_eq!(rec.n_cigar_ops, 3);
}

// r[verify record_store.set_alignment]
#[test]
fn set_alignment_then_sort_restores_order() {
    let raw1 = make_simple_record(0, 100, &[pack_cigar_op(4, CIGAR_M)], 4, b"read1");
    let raw2 = make_simple_record(0, 200, &[pack_cigar_op(4, CIGAR_M)], 4, b"read2");
    let raw3 = make_simple_record(0, 300, &[pack_cigar_op(4, CIGAR_M)], 4, b"read3");
    let mut store = RecordStore::new();
    store.push_raw(&raw1, &mut ()).unwrap();
    store.push_raw(&raw2, &mut ()).unwrap();
    store.push_raw(&raw3, &mut ()).unwrap();

    // Move record at index 2 (pos=300) to pos=50
    let new_cigar = pack_cigar(&[(4, CIGAR_M)]);
    store.set_alignment(ri(2), Pos0::new(50).unwrap(), &new_cigar).unwrap();

    store.sort_by_pos();

    // After sort: pos 50, 100, 200
    assert_eq!(store.record(ri(0)).unwrap().pos, Pos0::new(50).unwrap());
    assert_eq!(store.record(ri(1)).unwrap().pos, Pos0::new(100).unwrap());
    assert_eq!(store.record(ri(2)).unwrap().pos, Pos0::new(200).unwrap());

    // Verify qnames survived the sort correctly
    assert_eq!(store.record(ri(0)).unwrap().qname(), b"read3"); // was at index 2
    assert_eq!(store.record(ri(1)).unwrap().qname(), b"read1");
    assert_eq!(store.record(ri(2)).unwrap().qname(), b"read2");
}

// r[verify record_store.set_alignment]
#[test]
fn set_alignment_multiple_times_same_record() {
    let raw = make_simple_record(0, 100, &[pack_cigar_op(4, CIGAR_M)], 4, b"read1");
    let mut store = RecordStore::new();
    let idx = store.push_raw(&raw, &mut ()).unwrap().expect("kept");

    // First realignment
    let cigar1 = pack_cigar(&[(2, CIGAR_S), (2, CIGAR_M)]);
    store.set_alignment(idx, Pos0::new(90).unwrap(), &cigar1).unwrap();
    assert_eq!(store.record(idx).unwrap().pos, Pos0::new(90).unwrap());

    // Second realignment — should also work (appends again, more dead cigar bytes)
    let cigar2 = pack_cigar(&[(1, CIGAR_S), (3, CIGAR_M)]);
    store.set_alignment(idx, Pos0::new(80).unwrap(), &cigar2).unwrap();
    assert_eq!(store.record(idx).unwrap().pos, Pos0::new(80).unwrap());
    assert_eq!(store.record(idx).unwrap().cigar(), &cigar2);
    assert_eq!(store.record(idx).unwrap().matching_bases, 3);
}

// ---------------------------------------------------------------------------
// Property-based tests
// ---------------------------------------------------------------------------

/// Generate a valid CIGAR as (len, `op_type`) parts with a guaranteed minimum
/// query length. Returns parts and the total query-consuming length.
#[hegel::composite]
fn arb_cigar_parts_with_query_len(tc: &TestCase) -> Vec<(u32, u8)> {
    // 1..8 ops, M/I/D/S only (keeping it simple for tests). Seven ops of at
    // most 49 bases cannot exceed the 500-base ceiling callers rely on, and a
    // short draw is topped up with one M rather than rejected — rejecting here
    // makes the smallest input the generator can produce enormous, which in
    // turn makes shrinking useless.
    let mut parts: Vec<(u32, u8)> = tc.draw_silent(
        gs::vecs(gs::tuples!(
            gs::integers::<u32>().min_value(1).max_value(49),
            gs::sampled_from(&[CIGAR_M, CIGAR_I, CIGAR_D, CIGAR_S]),
        ))
        .min_size(1)
        .max_size(7),
    );
    let qlen = query_len_from_parts(&parts);
    if qlen < MIN_QUERY_LEN {
        parts.push((MIN_QUERY_LEN - qlen, CIGAR_M));
    }
    parts
}

/// The minimum query length `set_alignment` callers here need.
const MIN_QUERY_LEN: u32 = 4;

/// Generate a pair of CIGARs with the same query length (for `set_alignment`).
/// The second cigar is constructed to have exactly the right query length
/// by splitting it into a random number of M/I/S ops that sum to the target.
#[hegel::composite]
fn arb_cigar_pair(tc: &TestCase) -> (Vec<(u32, u8)>, Vec<(u32, u8)>) {
    let original = tc.draw_silent(arb_cigar_parts_with_query_len());
    let qlen = query_len_from_parts(&original);
    // Generate 0..2 D ops (don't consume query) and 1..3 query-consuming splits
    let del_lens: Vec<u32> =
        tc.draw_silent(gs::vecs(gs::integers::<u32>().min_value(1).max_value(19)).max_size(2));
    let weights: Vec<u32> = tc.draw_silent(
        gs::vecs(gs::integers::<u32>().min_value(1).max_value(99)).min_size(1).max_size(3),
    );

    // Distribute qlen across the weights proportionally
    let total_weight: u32 = weights.iter().sum();
    let mut new_parts: Vec<(u32, u8)> = Vec::new();
    let mut remaining = qlen;
    let query_ops = [CIGAR_M, CIGAR_S]; // query-consuming ops
    for (i, &w) in weights.iter().enumerate() {
        let portion = if i == weights.len() - 1 {
            remaining // last one gets the remainder
        } else {
            let p = (u64::from(qlen) * u64::from(w) / u64::from(total_weight)) as u32;
            p.max(1).min(remaining.saturating_sub(
                (weights.len() - i - 1) as u32, // reserve 1 for each remaining
            ))
        };
        if portion > 0 {
            new_parts.push((portion, query_ops[i % query_ops.len()]));
            remaining -= portion;
        }
    }
    // Interleave D ops
    for &d in &del_lens {
        new_parts.push((d, CIGAR_D));
    }
    (original, new_parts)
}

// r[verify record_store.set_alignment]
#[hegel::test(test_cases = 200)]
fn set_alignment_end_pos_matches_oracle(tc: TestCase) {
    let (orig_parts, new_parts) = tc.draw(arb_cigar_pair().print_as_debug());
    let pos = tc.draw(gs::integers::<u32>().max_value(9999));
    let new_pos = tc.draw(gs::integers::<u32>().max_value(9999));
    let qlen = query_len_from_parts(&orig_parts);
    let orig_cigar_packed: Vec<u32> =
        orig_parts.iter().map(|&(len, op)| pack_cigar_op(len, op)).collect();

    let raw = make_simple_record(0, pos as i32, &orig_cigar_packed, qlen, b"read1");
    let mut store = RecordStore::new();
    let idx = store.push_raw(&raw, &mut ()).unwrap().expect("kept");

    let new_cigar_bytes = pack_cigar(&new_parts);
    let new_pos_typed = Pos0::new(new_pos).unwrap();
    store.set_alignment(idx, new_pos_typed, &new_cigar_bytes).unwrap();

    let rec = store.record(idx).unwrap();

    // Oracle: compute expected values independently
    let expected_end = expected_end_pos(new_pos, &new_parts);
    assert_eq!(
        rec.end_pos,
        Pos0::new(expected_end).unwrap(),
        "end_pos mismatch for cigar {:?} at pos {}",
        new_parts,
        new_pos
    );

    assert_eq!(rec.matching_bases, expected_matching_bases(&new_parts));
    assert_eq!(rec.indel_bases, expected_indel_bases(&new_parts));
    assert_eq!(rec.n_cigar_ops, new_parts.len() as u16);
}

// r[verify record_store.set_alignment]
#[hegel::test(test_cases = 200)]
fn set_alignment_preserves_seq_and_qual(tc: TestCase) {
    let (orig_parts, new_parts) = tc.draw(arb_cigar_pair().print_as_debug());
    let pos = tc.draw(gs::integers::<u32>().max_value(9999));
    let qlen = query_len_from_parts(&orig_parts);
    let orig_cigar_packed: Vec<u32> =
        orig_parts.iter().map(|&(len, op)| pack_cigar_op(len, op)).collect();

    let raw = make_simple_record(0, pos as i32, &orig_cigar_packed, qlen, b"testread");
    let mut store = RecordStore::new();
    let idx = store.push_raw(&raw, &mut ()).unwrap().expect("kept");

    let orig_seq: Vec<Base> = store.record(idx).unwrap().seq().to_vec();
    let orig_qual: Vec<BaseQuality> = store.record(idx).unwrap().qual().to_vec();
    let orig_qname: Vec<u8> = store.record(idx).unwrap().qname().to_vec();

    let new_cigar_bytes = pack_cigar(&new_parts);
    store.set_alignment(idx, Pos0::new(pos).unwrap(), &new_cigar_bytes).unwrap();

    assert_eq!(store.record(idx).unwrap().seq(), &orig_seq[..], "seq changed after set_alignment");
    assert_eq!(
        store.record(idx).unwrap().qual(),
        &orig_qual[..],
        "qual changed after set_alignment"
    );
    assert_eq!(
        store.record(idx).unwrap().qname(),
        &orig_qname[..],
        "qname changed after set_alignment"
    );
}

// r[verify record_store.set_alignment]
#[hegel::test(test_cases = 200)]
fn sort_after_set_alignment_is_position_ordered(tc: TestCase) {
    let new_positions =
        tc.draw(gs::vecs(gs::integers::<u32>().max_value(49999)).min_size(3).max_size(9));
    let mut store = RecordStore::new();
    for (i, &_) in new_positions.iter().enumerate() {
        let raw = make_simple_record(
            0,
            (i * 100) as i32,
            &[pack_cigar_op(4, CIGAR_M)],
            4,
            format!("r{i}").as_bytes(),
        );
        store.push_raw(&raw, &mut ()).unwrap();
    }

    // Realign each record to its new position
    let cigar = pack_cigar(&[(4, CIGAR_M)]);
    for (i, &new_pos) in new_positions.iter().enumerate() {
        store
            .set_alignment(ri(u32::try_from(i).unwrap()), Pos0::new(new_pos).unwrap(), &cigar)
            .unwrap();
    }

    store.sort_by_pos();

    // Verify sorted order
    for (prev, rec) in store.records().zip(store.records().skip(1)) {
        assert!(rec.pos >= prev.pos, "records not sorted after set_alignment + sort_by_pos");
    }
}

// r[verify record_store.set_alignment.validation]
#[hegel::test(test_cases = 200)]
fn set_alignment_rejects_wrong_query_len(tc: TestCase) {
    let orig_parts = tc.draw(arb_cigar_parts_with_query_len().print_as_debug());
    let delta = tc.draw(gs::integers::<u32>().min_value(1).max_value(9));
    let qlen = query_len_from_parts(&orig_parts);
    let orig_cigar_packed: Vec<u32> =
        orig_parts.iter().map(|&(len, op)| pack_cigar_op(len, op)).collect();

    let raw = make_simple_record(0, 100, &orig_cigar_packed, qlen, b"read1");
    let mut store = RecordStore::new();
    let idx = store.push_raw(&raw, &mut ()).unwrap().expect("kept");

    let before = {
        let rec = store.record(idx).unwrap();
        (rec.pos, rec.end_pos, rec.cigar().to_vec())
    };

    // Build a cigar with wrong query length (qlen + delta)
    let wrong_cigar = pack_cigar(&[(qlen + delta, CIGAR_M)]);
    let err = store
        .set_alignment(idx, Pos0::new(200).unwrap(), &wrong_cigar)
        .expect_err("a CIGAR that spends a different number of query bases must be refused");
    assert!(
        matches!(
            err,
            DecodeError::CigarQueryLenMismatch { cigar_query_len, seq_len }
                if cigar_query_len == qlen + delta && seq_len == qlen
        ),
        "rejected with {err:?}, which does not name the mismatch"
    );

    // A refused realignment must leave the record exactly as it was — the
    // position it was asked to move to is a different one on purpose.
    let rec = store.record(idx).unwrap();
    assert_eq!((rec.pos, rec.end_pos, rec.cigar().to_vec()), before, "partially applied");
}

// ---------------------------------------------------------------------------
// Integration test: realignment workflow with real BAM
// ---------------------------------------------------------------------------

// r[verify record_store.set_alignment]
#[test]
fn realignment_workflow_with_real_bam() {
    use seqair::bam::reader::IndexedBamReader;
    use std::path::Path;

    let bam_path = Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.bam"));
    let mut reader = IndexedBamReader::open(bam_path).expect("open BAM");
    let tid = reader.header().tid("chr19").expect("tid");

    let mut store = RecordStore::new();
    reader
        .fetch_into(tid, Pos0::new(6_105_700).unwrap(), Pos0::new(6_105_800).unwrap(), &mut store)
        .expect("fetch");

    assert!(!store.is_empty());

    // Verify tid is populated for all records
    for i in store.indices() {
        assert_eq!(store.record(i).unwrap().tid, tid as i32, "tid mismatch for record {i}");
    }

    // "Realign" the first record by shifting it 10bp left with same cigar
    let orig_cigar = store.record(ri(0)).unwrap().cigar().to_vec();
    let orig_seq: Vec<Base> = store.record(ri(0)).unwrap().seq().to_vec();
    let orig_qual: Vec<BaseQuality> = store.record(ri(0)).unwrap().qual().to_vec();
    let orig_pos = store.record(ri(0)).unwrap().pos;

    let new_pos = Pos0::new(orig_pos.as_u32().saturating_sub(10)).unwrap();
    store.set_alignment(ri(0), new_pos, &orig_cigar).unwrap();

    assert_eq!(store.record(ri(0)).unwrap().pos, new_pos);
    assert_eq!(store.record(ri(0)).unwrap().seq(), &orig_seq);
    assert_eq!(store.record(ri(0)).unwrap().qual(), &orig_qual);

    // Sort and verify order
    store.sort_by_pos();
    for (i, (prev, rec)) in store.records().zip(store.records().skip(1)).enumerate() {
        assert!(rec.pos >= prev.pos, "sort order broken at index {}", i + 1);
    }
}
