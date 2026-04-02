//! Shared test helpers for building synthetic BAM records.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::indexing_slicing)]
#![allow(clippy::arithmetic_side_effects)]
#![allow(dead_code)]

pub fn cigar_op(len: u32, op: u8) -> u32 {
    (len << 4) | u32::from(op)
}

pub fn cigar_bytes(ops: &[u32]) -> Vec<u8> {
    ops.iter().flat_map(|op| op.to_le_bytes()).collect()
}

/// Build a synthetic BAM record with a simple N×M CIGAR.
pub fn make_record(tid: i32, pos: i32, flags: u16, mapq: u8, seq_len: u32) -> Vec<u8> {
    make_record_with_cigar(tid, pos, flags, mapq, &[cigar_op(seq_len, 0)], seq_len)
}

/// Build a synthetic BAM record with a custom CIGAR.
pub fn make_record_with_cigar(
    tid: i32,
    pos: i32,
    flags: u16,
    mapq: u8,
    cigar_ops: &[u32],
    seq_len: u32,
) -> Vec<u8> {
    let name = b"read\0";
    let name_len = name.len();
    let n_cigar_ops = u16::try_from(cigar_ops.len()).unwrap();
    let cigar_bytes_len = cigar_ops.len() * 4;
    let seq_bytes = (seq_len as usize).div_ceil(2);

    let total = 32 + name_len + cigar_bytes_len + seq_bytes + seq_len as usize;
    let mut raw = vec![0u8; total];

    raw[0..4].copy_from_slice(&tid.to_le_bytes());
    raw[4..8].copy_from_slice(&pos.to_le_bytes());
    raw[8] = name_len as u8;
    raw[9] = mapq;
    raw[10..12].copy_from_slice(&0u16.to_le_bytes());
    raw[12..14].copy_from_slice(&n_cigar_ops.to_le_bytes());
    raw[14..16].copy_from_slice(&flags.to_le_bytes());
    raw[16..20].copy_from_slice(&seq_len.to_le_bytes());
    raw[20..32].fill(0);
    raw[32..32 + name_len].copy_from_slice(name);

    let cigar_start = 32 + name_len;
    for (i, op) in cigar_ops.iter().enumerate() {
        raw[cigar_start + i * 4..cigar_start + (i + 1) * 4].copy_from_slice(&op.to_le_bytes());
    }

    let seq_start = cigar_start + cigar_bytes_len;
    for i in 0..seq_bytes {
        raw[seq_start + i] = 0x11; // A,A
    }

    let qual_start = seq_start + seq_bytes;
    for i in 0..seq_len as usize {
        raw[qual_start + i] = 30;
    }

    raw
}

/// Build a BAM record with a custom qname and packed sequence bytes.
/// `seq_packed` is the 4-bit packed sequence. `seq_len` is the number of bases.
pub fn make_named_record(
    qname: &[u8],
    tid: i32,
    pos: i32,
    flags: u16,
    mapq: u8,
    seq_len: u32,
    seq_packed: &[u8],
) -> Vec<u8> {
    let mut name_with_nul = qname.to_vec();
    name_with_nul.push(0);
    // Pad to 4-byte alignment as BAM spec requires
    while !name_with_nul.len().is_multiple_of(4) {
        name_with_nul.push(0);
    }
    let name_len = name_with_nul.len();
    let n_cigar_ops = 1u16;
    let cigar = cigar_op(seq_len, 0); // simple NM
    let seq_bytes = (seq_len as usize).div_ceil(2);

    let total = 32 + name_len + 4 + seq_bytes + seq_len as usize;
    let mut raw = vec![0u8; total];

    raw[0..4].copy_from_slice(&tid.to_le_bytes());
    raw[4..8].copy_from_slice(&pos.to_le_bytes());
    raw[8] = name_len as u8;
    raw[9] = mapq;
    raw[10..12].copy_from_slice(&0u16.to_le_bytes());
    raw[12..14].copy_from_slice(&n_cigar_ops.to_le_bytes());
    raw[14..16].copy_from_slice(&flags.to_le_bytes());
    raw[16..20].copy_from_slice(&seq_len.to_le_bytes());
    raw[20..32].fill(0);
    raw[32..32 + name_len].copy_from_slice(&name_with_nul);

    let cigar_start = 32 + name_len;
    raw[cigar_start..cigar_start + 4].copy_from_slice(&cigar.to_le_bytes());

    let seq_start = cigar_start + 4;
    let copy_len = seq_packed.len().min(seq_bytes);
    raw[seq_start..seq_start + copy_len].copy_from_slice(&seq_packed[..copy_len]);

    let qual_start = seq_start + seq_bytes;
    for i in 0..seq_len as usize {
        raw[qual_start + i] = 30;
    }

    raw
}

/// BAM 4-bit encoding for common bases
pub const BASE_A: u8 = 1;
pub const BASE_C: u8 = 2;
pub const BASE_G: u8 = 4;
pub const BASE_T: u8 = 8;

/// Pack two 4-bit bases into one byte
pub fn pack_bases(hi: u8, lo: u8) -> u8 {
    (hi << 4) | lo
}

// --- Proptest strategies for synthetic reads ---

use proptest::prelude::*;

/// A generated read with all the metadata needed to validate pileup behavior.
#[derive(Debug, Clone)]
pub struct SyntheticRead {
    pub pos: i32,
    pub flags: u16,
    pub mapq: u8,
    pub cigar_ops: Vec<(u32, u8)>, // (length, op_code)
    pub seq_len: u32,              // total query-consuming bases
    pub ref_span: u32,             // total ref-consuming bases
    pub raw: Vec<u8>,              // ready-to-push BAM bytes
}

impl SyntheticRead {
    /// Reference positions where this read has a qpos (M/=/X ops only).
    pub fn covered_ref_positions(&self) -> Vec<i64> {
        let mut positions = Vec::new();
        let mut ref_off = 0i64;
        let pos = i64::from(self.pos);
        for &(len, op) in &self.cigar_ops {
            let consumes_ref = matches!(op, 0 | 2 | 3 | 7 | 8); // M,D,N,=,X
            let consumes_query = matches!(op, 0 | 1 | 4 | 7 | 8); // M,I,S,=,X
            if consumes_ref && consumes_query {
                // M, =, X: produces a qpos at each ref position
                for j in 0..len {
                    positions.push(pos + ref_off + i64::from(j));
                }
            }
            if consumes_ref {
                ref_off += i64::from(len);
            }
        }
        positions
    }

    /// If the given reference position falls inside a D op, returns the D op's length.
    /// Otherwise returns None.
    pub fn del_len_at(&self, ref_pos: i64) -> Option<u32> {
        let mut ref_off = 0i64;
        let pos = i64::from(self.pos);
        for &(len, op) in &self.cigar_ops {
            let consumes_ref = matches!(op, 0 | 2 | 3 | 7 | 8);
            if consumes_ref {
                let ref_start = pos + ref_off;
                let ref_end = ref_start + i64::from(len);
                if ref_pos >= ref_start && ref_pos < ref_end {
                    return if op == 2 { Some(len) } else { None };
                }
                ref_off += i64::from(len);
            }
        }
        None
    }

    /// Expected qpos at a given reference position, or None.
    pub fn qpos_at(&self, ref_pos: i64) -> Option<usize> {
        let mut ref_off = 0i64;
        let mut query_off = 0usize;
        let pos = i64::from(self.pos);
        for &(len, op) in &self.cigar_ops {
            let consumes_ref = matches!(op, 0 | 2 | 3 | 7 | 8);
            let consumes_query = matches!(op, 0 | 1 | 4 | 7 | 8);
            if consumes_ref {
                let ref_start = pos + ref_off;
                let ref_end = ref_start + i64::from(len);
                if ref_pos >= ref_start && ref_pos < ref_end {
                    return if consumes_query {
                        Some(query_off + (ref_pos - ref_start) as usize)
                    } else {
                        None // deletion or ref-skip
                    };
                }
                ref_off += i64::from(len);
            }
            if consumes_query {
                query_off += len as usize;
            }
        }
        None
    }
}

/// Strategy for generating CIGAR patterns of various complexity classes.
fn arb_cigar_class() -> impl Strategy<Value = Vec<(u32, u8)>> {
    prop_oneof![
        // Simple: just M (exercises Linear fast path)
        (10u32..=200).prop_map(|len| vec![(len, 0u8)]),
        // Soft-clipped: S + M or S + M + S
        (1u32..=20, 10u32..=150, prop::option::of(1u32..=20)).prop_map(|(s1, m, s2)| {
            let mut ops = vec![(s1, 4u8), (m, 0u8)];
            if let Some(s) = s2 {
                ops.push((s, 4u8));
            }
            ops
        }),
        // Hard-clipped: H + M + H
        (1u32..=10, 10u32..=100, 1u32..=10)
            .prop_map(|(h1, m, h2)| { vec![(h1, 5u8), (m, 0u8), (h2, 5u8)] }),
        // Deletion: M + D + M
        (5u32..=50, 1u32..=20, 5u32..=50)
            .prop_map(|(m1, d, m2)| { vec![(m1, 0u8), (d, 2u8), (m2, 0u8)] }),
        // Deletion then insertion: M + D + I + M (orphan insertion, tests no_orphan_insertions rule)
        (5u32..=50, 1u32..=20, 1u32..=10, 5u32..=50)
            .prop_map(|(m1, d, i, m2)| { vec![(m1, 0u8), (d, 2u8), (i, 1u8), (m2, 0u8)] }),
        // Insertion: M + I + M
        (5u32..=50, 1u32..=20, 5u32..=50)
            .prop_map(|(m1, i, m2)| { vec![(m1, 0u8), (i, 1u8), (m2, 0u8)] }),
        // Intron/RefSkip: M + N + M (RNA-seq like)
        (10u32..=50, 100u32..=5000, 10u32..=50)
            .prop_map(|(m1, n, m2)| { vec![(m1, 0u8), (n, 3u8), (m2, 0u8)] }),
        // Complex: S + M + I + M + D + M + S (exercises Complex path)
        (1u32..=10, 5u32..=30, 1u32..=10, 5u32..=30, 1u32..=15, 5u32..=30, 1u32..=10).prop_map(
            |(s1, m1, i, m2, d, m3, s2)| {
                vec![(s1, 4u8), (m1, 0u8), (i, 1u8), (m2, 0u8), (d, 2u8), (m3, 0u8), (s2, 4u8)]
            }
        ),
        // Multi-exon RNA: M + N + M + N + M (many ops, exercises binary search)
        (10u32..=40, 100u32..=2000, 10u32..=40, 100u32..=2000, 10u32..=40).prop_map(
            |(m1, n1, m2, n2, m3)| { vec![(m1, 0u8), (n1, 3u8), (m2, 0u8), (n2, 3u8), (m3, 0u8)] }
        ),
        // =, X ops: explicit match/mismatch (tests SeqMatch/SeqMismatch)
        (5u32..=30, 5u32..=30, 5u32..=30)
            .prop_map(|(eq, x, eq2)| { vec![(eq, 7u8), (x, 8u8), (eq2, 7u8)] }),
        // Kitchen sink: H + S + M + I + M + D + M + N + M + S + H
        (
            1u32..=5,
            1u32..=10,
            5u32..=20,
            1u32..=5,
            5u32..=20,
            1u32..=10,
            5u32..=20,
            100u32..=500,
            5u32..=20,
            1u32..=10,
            1u32..=5
        )
            .prop_map(|(h1, s1, m1, i, m2, d, m3, n, m4, s2, h2)| {
                vec![
                    (h1, 5u8),
                    (s1, 4u8),
                    (m1, 0u8),
                    (i, 1u8),
                    (m2, 0u8),
                    (d, 2u8),
                    (m3, 0u8),
                    (n, 3u8),
                    (m4, 0u8),
                    (s2, 4u8),
                    (h2, 5u8),
                ]
            }),
    ]
}

/// Generate a single synthetic read with a specific CIGAR pattern.
pub fn arb_read() -> impl Strategy<Value = SyntheticRead> {
    (
        0i32..500,
        prop_oneof![Just(99u16), Just(163u16), Just(83u16), Just(147u16), Just(0u16),],
        (0u8..=60),
        arb_cigar_class(),
    )
        .prop_map(|(pos, flags, mapq, cigar_ops)| {
            let seq_len: u32 = cigar_ops
                .iter()
                .filter(|&&(_, op)| matches!(op, 0 | 1 | 4 | 7 | 8))
                .map(|&(len, _)| len)
                .sum();
            let ref_span: u32 = cigar_ops
                .iter()
                .filter(|&&(_, op)| matches!(op, 0 | 2 | 3 | 7 | 8))
                .map(|&(len, _)| len)
                .sum();
            let packed_ops: Vec<u32> =
                cigar_ops.iter().map(|&(len, op)| (len << 4) | u32::from(op)).collect();
            let raw = make_record_with_cigar(0, pos, flags, mapq, &packed_ops, seq_len);
            SyntheticRead { pos, flags, mapq, cigar_ops, seq_len, ref_span, raw }
        })
        .prop_filter("must have ref-consuming ops", |r| r.ref_span > 0)
}

/// Generate a sorted set of reads for a region.
pub fn arb_read_set(max_reads: usize) -> impl Strategy<Value = Vec<SyntheticRead>> {
    prop::collection::vec(arb_read(), 1..=max_reads).prop_map(|mut reads| {
        reads.sort_by_key(|r| r.pos);
        reads
    })
}
