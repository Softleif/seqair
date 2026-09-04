#![no_main]

//! Fuzz the AlignedPairs layered iterator chain end-to-end.
//!
//! Drives random structured BAM records through the full chain
//! (`aligned_pairs_with_read` → `with_reference` → consume) plus the
//! `matches_only` and `nm`/`md` derivatives. The goal is to surface panics
//! from adversarial CIGAR/seq/qual combinations the unit + property tests
//! don't reach: degenerate ops, length mismatches, runaway saturating
//! arithmetic, etc.

use arbitrary::Arbitrary;
use libfuzzer_sys::fuzz_target;
use seqair::bam::aligned_pairs::{AlignedPair, AlignedPairs, MatchKind};
use seqair::bam::aligned_pairs_view::{AlignedPairWithRead, AlignedPairWithRef};
use seqair::bam::pileup::RefSeq;
use seqair::bam::record_store::RecordStore;
use seqair_types::{Base, Pos0};
use std::rc::Rc;

const OP_COUNT: u8 = 9;

const fn consumes_query(op: u8) -> bool {
    matches!(op, 0 | 1 | 4 | 7 | 8) // M, I, S, =, X
}

const fn consumes_ref(op: u8) -> bool {
    matches!(op, 0 | 2 | 3 | 7 | 8) // M, D, N, =, X
}

#[derive(Arbitrary, Debug, Clone)]
struct CigarOp {
    op_type: u8,
    length: u16,
}

impl CigarOp {
    fn bam_op(&self) -> u8 {
        self.op_type % OP_COUNT
    }

    fn length_clamped(&self) -> u32 {
        // Allow zero-length ops (the iterator should skip them).
        u32::from(self.length).min(100)
    }
}

#[derive(Arbitrary, Debug)]
struct FuzzInput {
    pos: u16,
    flag_bits: u8,
    cigar_ops: Vec<CigarOp>,
    seq_bases: Vec<u8>,
    quals: Vec<u8>,
    /// Reference window start offset relative to record pos, signed.
    /// Allows windows that start before, at, or after the record's pos.
    ref_window_offset: i16,
    /// Reference window length.
    ref_window_len: u16,
    ref_bases: Vec<u8>,
}

/// Build a raw BAM record matching `parse_header`'s layout. Same shape as
/// `fuzz_structured_pileup` but a much smaller subset — we don't need the
/// full pileup machinery here.
fn build_raw_record(input: &FuzzInput) -> Vec<u8> {
    let raw_ops: Vec<(u8, u32)> =
        input.cigar_ops.iter().take(20).map(|op| (op.bam_op(), op.length_clamped())).collect();

    let l_seq: u32 = raw_ops
        .iter()
        .filter(|(op, _)| consumes_query(*op))
        .map(|(_, len)| *len)
        .fold(0u32, |acc, v| acc.saturating_add(v));

    let n_cigar_ops = raw_ops.len() as u16;
    let qname = b"f\0";
    let l_read_name: u8 = qname.len() as u8;
    let flags: u16 = u16::from(input.flag_bits & 0xD1) & !0x0004;
    let pos_val = i32::from(input.pos);

    let mut raw =
        Vec::with_capacity(32 + qname.len() + (n_cigar_ops as usize) * 4 + (l_seq as usize));
    raw.extend_from_slice(&0i32.to_le_bytes());
    raw.extend_from_slice(&pos_val.to_le_bytes());
    raw.push(l_read_name);
    raw.push(0u8);
    raw.extend_from_slice(&0u16.to_le_bytes());
    raw.extend_from_slice(&n_cigar_ops.to_le_bytes());
    raw.extend_from_slice(&flags.to_le_bytes());
    raw.extend_from_slice(&l_seq.to_le_bytes());
    raw.extend_from_slice(&(-1i32).to_le_bytes());
    raw.extend_from_slice(&(-1i32).to_le_bytes());
    raw.extend_from_slice(&0i32.to_le_bytes());
    raw.extend_from_slice(qname);

    for (op, len) in &raw_ops {
        let packed: u32 = (len << 4) | u32::from(*op);
        raw.extend_from_slice(&packed.to_le_bytes());
    }

    let seq_len = l_seq as usize;
    let packed_len = seq_len.div_ceil(2);
    let mut packed_seq = vec![0u8; packed_len];
    let nibble_table: [u8; 6] = [1, 2, 4, 8, 15, 0];
    for i in 0..seq_len {
        let src = input.seq_bases.get(i).copied().unwrap_or(15);
        let nibble = nibble_table[(src as usize) % nibble_table.len()];
        if i % 2 == 0 {
            if let Some(b) = packed_seq.get_mut(i / 2) {
                *b = nibble << 4;
            }
        } else if let Some(b) = packed_seq.get_mut(i / 2) {
            *b |= nibble;
        }
    }
    raw.extend_from_slice(&packed_seq);

    for i in 0..seq_len {
        let q = input.quals.get(i).copied().unwrap_or(0xFF);
        raw.push(q);
    }

    raw
}

fn build_ref_seq(input: &FuzzInput) -> Option<RefSeq> {
    // Anchor the window at `pos + ref_window_offset`, but clamp into Pos0
    // range. Length capped at 1000 so the fuzzer doesn't allocate huge.
    let pos_i64 = i64::from(input.pos);
    let start_i64 = pos_i64.saturating_add(i64::from(input.ref_window_offset));
    if start_i64 < 0 || start_i64 > i64::from(i32::MAX) {
        return None;
    }
    let start = Pos0::new(start_i64 as u32)?;

    let len = (input.ref_window_len as usize).min(1000);
    let bases: Vec<Base> = (0..len)
        .map(|i| {
            let src = input.ref_bases.get(i).copied().unwrap_or(b'N');
            Base::from(src)
        })
        .collect();
    Some(RefSeq::new(Rc::from(bases), start))
}

// Drive every public AlignedPairs* path with the same record + reference.
// Each iterator must complete without panicking, regardless of input shape.
fuzz_target!(|input: FuzzInput| {
    let raw = build_raw_record(&input);

    let mut store = RecordStore::<()>::new();
    match store.push_raw(&raw, &mut ()) {
        Ok(Some(_)) => (),
        _ => return,
    };
    if store.is_empty() {
        return;
    }
    let rec = store.record(0);
    let cigar = match rec.cigar(&store) {
        Ok(c) => c,
        Err(_) => return,
    };

    // ── Layer 0: bare AlignedPairs (every option combination) ──
    {
        let it = AlignedPairs::new(rec.pos, cigar);
        consume_bare(it);
    }
    {
        let it = AlignedPairs::new(rec.pos, cigar).with_soft_clips();
        consume_bare(it);
    }
    {
        let it = AlignedPairs::new(rec.pos, cigar).full();
        consume_bare(it);
    }

    // ── Layer 0 derivative: matches_only ──
    {
        let count = AlignedPairs::new(rec.pos, cigar).matches_only().count();
        let _ = count;
    }

    // ── Layer 0: size_hint must agree with collected count ──
    {
        let it = AlignedPairs::new(rec.pos, cigar);
        let predicted = it.size_hint().1.unwrap_or(0);
        let actual = AlignedPairs::new(rec.pos, cigar).count();
        // Documented invariant: size_hint upper == actual count.
        assert_eq!(predicted, actual, "size_hint disagreed with collected count");
    }

    // ── Layer 1: with_read ──
    let with_read_result = rec.aligned_pairs_with_read(&store).ok();
    if let Some(it) = with_read_result.clone() {
        consume_with_read(it);
    }
    if let Some(it) = with_read_result.clone() {
        consume_with_read(it.with_soft_clips());
    }
    if let Some(it) = with_read_result.clone() {
        consume_with_read(it.full());
    }

    // ── Layer 1 derivative: matches_only via with_read ──
    if let Some(it) = with_read_result.clone() {
        let _ = it.matches_only().count();
    }

    // ── Layer 2: with_reference ──
    let ref_seq = match build_ref_seq(&input) {
        Some(rs) => rs,
        None => return,
    };
    if let Some(it) = with_read_result.clone() {
        let with_ref = it.with_reference(&ref_seq);
        consume_with_ref(with_ref);
    }

    // ── nm() / md() — exercised on a fresh chain to avoid consuming the cached one. ──
    if let Ok(it) = rec.aligned_pairs_with_read(&store) {
        let _nm = it.with_reference(&ref_seq).nm();
    }
    if let Ok(it) = rec.aligned_pairs_with_read(&store) {
        // md() may error on missing ref — that's fine, we just want no panic.
        let _ = it.with_reference(&ref_seq).md();
    }
});

fn consume_bare(it: AlignedPairs<'_>) {
    for pair in it {
        // Touch every field to ensure no destructure panic.
        match pair {
            AlignedPair::Match { qpos, rpos, kind } => {
                let _ = (qpos, rpos, kind);
                debug_assert!(matches!(
                    kind,
                    MatchKind::Match | MatchKind::SeqMatch | MatchKind::SeqMismatch
                ));
            }
            AlignedPair::Insertion { qpos, insert_len } => {
                let _ = (qpos, insert_len);
                debug_assert!(insert_len > 0, "zero-len summary leaked");
            }
            AlignedPair::Deletion { rpos, del_len } => {
                let _ = (rpos, del_len);
                debug_assert!(del_len > 0, "zero-len deletion summary leaked");
            }
            AlignedPair::RefSkip { rpos, skip_len } => {
                let _ = (rpos, skip_len);
                debug_assert!(skip_len > 0, "zero-len refskip summary leaked");
            }
            AlignedPair::SoftClip { qpos, len } => {
                let _ = (qpos, len);
                debug_assert!(len > 0, "zero-len softclip summary leaked");
            }
            AlignedPair::Padding { len } => {
                let _ = len;
                debug_assert!(len > 0, "zero-len padding summary leaked");
            }
            AlignedPair::Unknown { code, len } => {
                let _ = (code, len);
                debug_assert!(len > 0, "zero-len unknown summary leaked");
            }
        }
    }
}

fn consume_with_read<'cigar, 'read>(
    it: seqair::bam::aligned_pairs_view::AlignedPairsWithRead<'cigar, 'read>,
) {
    for ev in it {
        match ev {
            AlignedPairWithRead::Match { qpos, rpos, query, qual, kind } => {
                let _ = (qpos, rpos, query, qual, kind);
            }
            AlignedPairWithRead::Insertion { qpos, query, qual } => {
                let _ = qpos;
                debug_assert_eq!(query.len(), qual.len(), "insert query/qual lengths must match");
                debug_assert!(!query.is_empty(), "zero-len insertion leaked");
            }
            AlignedPairWithRead::Deletion { rpos, del_len } => {
                let _ = (rpos, del_len);
            }
            AlignedPairWithRead::RefSkip { rpos, skip_len } => {
                let _ = (rpos, skip_len);
            }
            AlignedPairWithRead::SoftClip { qpos, query, qual } => {
                let _ = qpos;
                debug_assert_eq!(query.len(), qual.len());
                debug_assert!(!query.is_empty());
            }
            AlignedPairWithRead::Padding { len } => {
                let _ = len;
            }
            AlignedPairWithRead::Unknown { code, len } => {
                let _ = (code, len);
            }
        }
    }
}

fn consume_with_ref<'cigar, 'read, 'ref_seq>(
    it: seqair::bam::aligned_pairs_view::AlignedPairsWithRef<'cigar, 'read, 'ref_seq>,
) {
    for ev in it {
        match ev {
            AlignedPairWithRef::Match { qpos, rpos, query, qual, kind, ref_base } => {
                let _ = (qpos, rpos, query, qual, kind, ref_base);
            }
            AlignedPairWithRef::Insertion { qpos, query, qual } => {
                let _ = (qpos, query, qual);
            }
            AlignedPairWithRef::Deletion { rpos, del_len, ref_bases } => {
                let _ = (rpos, del_len);
                if let Some(bases) = ref_bases {
                    debug_assert_eq!(
                        bases.len(),
                        del_len as usize,
                        "ref_bases len must match del_len when present"
                    );
                }
            }
            AlignedPairWithRef::RefSkip { rpos, skip_len } => {
                let _ = (rpos, skip_len);
            }
            AlignedPairWithRef::SoftClip { qpos, query, qual } => {
                let _ = (qpos, query, qual);
            }
            AlignedPairWithRef::Padding { len } => {
                let _ = len;
            }
            AlignedPairWithRef::Unknown { code, len } => {
                let _ = (code, len);
            }
        }
    }
}

#[allow(dead_code)]
fn _consumes_ref_used() {
    // Reference exists to silence dead_code; consumes_ref is part of the
    // intent doc even though build_raw_record doesn't currently use it.
    let _ = consumes_ref(0);
}
