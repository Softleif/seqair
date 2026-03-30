//! Tests for Base-typed sequence decoding.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::indexing_slicing)]
use seqair::bam::record_store::RecordStore;
use seqair::bam::{Pos, Zero};
use seqair_types::Base;

// r[verify base_decode.table]
// r[verify base_decode.decode]
#[test]
fn decode_standard_bases() {
    // packed_seq: 0x12 = A(1)|C(2), 0x48 = G(4)|T(8) → ACGT
    let raw = make_record(&[0x12, 0x48], 4);
    let mut store = RecordStore::new();
    let idx = store.push_raw(&raw).unwrap();

    assert_eq!(store.seq_at(idx, 0), Base::A);
    assert_eq!(store.seq_at(idx, 1), Base::C);
    assert_eq!(store.seq_at(idx, 2), Base::G);
    assert_eq!(store.seq_at(idx, 3), Base::T);
}

// r[verify base_decode.table]
#[test]
fn decode_unknown_and_n() {
    // nibble 0 (=) and nibble 15 (N) both map to Base::Unknown
    // 0x0F = =(0)|N(15)
    let raw = make_record(&[0x0F], 2);
    let mut store = RecordStore::new();
    let idx = store.push_raw(&raw).unwrap();

    assert_eq!(store.seq_at(idx, 0), Base::Unknown);
    assert_eq!(store.seq_at(idx, 1), Base::Unknown);
}

// r[verify base_decode.table]
#[test]
fn decode_iupac_ambiguity_maps_to_unknown() {
    // nibble 3 = M (IUPAC A|C), nibble 5 = R (IUPAC A|G) → both Unknown
    // 0x35 = M(3)|R(5)
    let raw = make_record(&[0x35], 2);
    let mut store = RecordStore::new();
    let idx = store.push_raw(&raw).unwrap();

    assert_eq!(store.seq_at(idx, 0), Base::Unknown);
    assert_eq!(store.seq_at(idx, 1), Base::Unknown);
}

// r[verify base_decode.slab]
#[test]
fn bases_stored_in_separate_slab() {
    let raw = make_record(&[0x12, 0x48], 4);
    let mut store = RecordStore::new();
    let idx = store.push_raw(&raw).unwrap();

    // Access the full base slice — should be 4 Base values
    let bases = store.seq(idx);
    assert_eq!(bases.len(), 4);
    assert_eq!(bases, &[Base::A, Base::C, Base::G, Base::T]);
}

// r[verify base_decode.alignment]
#[test]
fn pileup_alignment_has_base_type() {
    use seqair::bam::{pileup::PileupEngine, record_store::RecordStore};

    let raw = make_record(&[0x12, 0x48], 4);
    let mut store = RecordStore::new();
    store.push_raw(&raw).unwrap();

    let mut engine = PileupEngine::new(store, Pos::<Zero>::new(100), Pos::<Zero>::new(103));
    engine.set_max_depth(1000);

    let col = engine.next().expect("should have a column");
    let aln = col.alignments().next().expect("should have an alignment");

    // aln.base() returns Option<Base> — this verifies the type via the convenience method
    let base: Base = aln.base().unwrap();
    assert!(matches!(base, Base::A | Base::C | Base::G | Base::T | Base::Unknown));
}

/// Build a minimal raw BAM record for testing sequence decoding.
fn make_record(packed_seq: &[u8], seq_len: u32) -> Vec<u8> {
    let name = b"r\0";
    let name_len = name.len();
    let n_cigar_ops = 1u16;
    let cigar_bytes = 4;
    let qual_len = seq_len as usize;

    let total = 32 + name_len + cigar_bytes + packed_seq.len() + qual_len;
    let mut raw = vec![0u8; total];

    // tid=0, pos=100
    raw[4..8].copy_from_slice(&100i32.to_le_bytes());
    raw[8] = name_len as u8;
    raw[9] = 60; // mapq
    raw[12..14].copy_from_slice(&n_cigar_ops.to_le_bytes());
    raw[14..16].copy_from_slice(&0x63u16.to_le_bytes()); // flags
    raw[16..20].copy_from_slice(&seq_len.to_le_bytes());

    let mut off = 32;
    raw[off..off + name_len].copy_from_slice(name);
    off += name_len;

    // CIGAR: seq_len M
    let cigar_op = seq_len << 4; // M
    raw[off..off + 4].copy_from_slice(&cigar_op.to_le_bytes());
    off += 4;

    raw[off..off + packed_seq.len()].copy_from_slice(packed_seq);
    off += packed_seq.len();

    for i in 0..qual_len {
        raw[off + i] = 30;
    }

    raw
}
