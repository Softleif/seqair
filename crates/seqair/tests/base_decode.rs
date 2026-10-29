//! Tests for how SEQ is stored (as `IupacBase`) and read by the pileup (as `Base`).
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    reason = "test code"
)]
#![allow(clippy::arithmetic_side_effects, reason = "test code")]
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    reason = "test code with known small values"
)]
mod helpers;
use seqair::bam::Pos0;
use seqair::bam::record_store::RecordStore;
use seqair_types::QPos;
use seqair_types::{Base, IupacBase};
use std::num::NonZeroU32;

// r[verify base_decode.slab+1]
// r[verify record_store.seq_lossless]
#[test]
fn decode_standard_bases() {
    // packed_seq: 0x12 = A(1)|C(2), 0x48 = G(4)|T(8) → ACGT
    let raw = make_record(&[0x12, 0x48], 4);
    let mut store = RecordStore::new();
    let idx = store.push_raw(&raw, helpers::header(), &mut ()).unwrap().expect("kept");

    assert_eq!(store.record(idx).unwrap().base_at(QPos::new(0)), IupacBase::A);
    assert_eq!(store.record(idx).unwrap().base_at(QPos::new(1)), IupacBase::C);
    assert_eq!(store.record(idx).unwrap().base_at(QPos::new(2)), IupacBase::G);
    assert_eq!(store.record(idx).unwrap().base_at(QPos::new(3)), IupacBase::T);
}

// r[verify base_decode.slab+1]
// r[verify record_store.seq_lossless]
#[test]
fn decode_keeps_eq_n_and_ambiguity_codes() {
    // 0x0F = =(0)|N(15), 0x35 = M(3)|R(5)
    let raw = make_record(&[0x0F, 0x35], 4);
    let mut store = RecordStore::new();
    let idx = store.push_raw(&raw, helpers::header(), &mut ()).unwrap().expect("kept");

    let rec = store.record(idx).unwrap();
    assert_eq!(rec.seq(), &[IupacBase::Eq, IupacBase::N, IupacBase::M, IupacBase::R]);
    assert_eq!(rec.base_at(QPos::new(4)), IupacBase::N, "past the read");
}

// r[verify base_decode.alignment+1]
// r[verify pileup.base_from_iupac]
#[test]
fn pileup_converts_codes_to_base() {
    use seqair::bam::{pileup::PileupEngine, record_store::RecordStore};

    // A, =, M, T
    let raw = make_record(&[0x10, 0x38], 4);
    let mut store = RecordStore::new();
    store.push_raw(&raw, helpers::header(), &mut ()).unwrap();

    let mut engine = PileupEngine::new(
        store.prepare_for_pileup().input,
        (Pos0::new(100).unwrap()..=Pos0::new(103).unwrap()).into(),
    );
    engine.set_max_depth(NonZeroU32::new(1000).unwrap());

    let mut bases = Vec::new();
    while let Some(col) = engine.pileups() {
        let aln = col.alignments().next().expect("should have an alignment");
        let base: Base = aln.base().unwrap();
        bases.push(base);
    }
    assert_eq!(bases, [Base::A, Base::Unknown, Base::Unknown, Base::T]);
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
