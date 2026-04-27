#![no_main]

use arbitrary::Arbitrary;
use libfuzzer_sys::fuzz_target;
use seqair::bam::pileup::PileupEngine;
use seqair::bam::record_store::RecordStore;
use seqair_types::{Offset, Pos0};

// Valid BAM CIGAR op codes (0..=8, skipping 6=P which is rare and adds no coverage)
// M=0, I=1, D=2, N=3, S=4, H=5, P=6, ==7, X=8
const OP_COUNT: u8 = 9;

// Which ops consume the query sequence
const fn consumes_query(op: u8) -> bool {
    matches!(op, 0 | 1 | 4 | 7 | 8) // M, I, S, =, X
}

/// A single structured CIGAR operation.
#[derive(Arbitrary, Debug, Clone)]
struct CigarOp {
    /// Mapped to a valid BAM op via `% OP_COUNT`.
    op_type: u8,
    /// Operation length, clamped to 1..=200.
    length: u16,
}

impl CigarOp {
    fn bam_op(&self) -> u8 {
        self.op_type % OP_COUNT
    }

    fn length_clamped(&self) -> u32 {
        (self.length as u32).clamp(1, 200)
    }
}

/// A structurally valid fuzz BAM record.
#[derive(Arbitrary, Debug)]
struct FuzzRecord {
    /// 0-based alignment position (small for pileup overlap).
    pos: u16,
    mapq: u8,
    /// Lower flag bits for strand / paired-end variation.
    /// Bit 0x4 (unmapped) is always cleared so the record reaches the pileup.
    flag_bits: u8,
    /// Up to 8 CIGAR operations.
    cigar_ops: Vec<CigarOp>,
    /// Raw bases (any byte — will be packed 4-bit then decoded back to Base).
    seq_bases: Vec<u8>,
    /// Quality scores.
    quals: Vec<u8>,
}

/// Top-level fuzz input: a small set of records + a region.
#[derive(Arbitrary, Debug)]
struct FuzzPileupInput {
    records: Vec<FuzzRecord>,
    region_start: u16,
    region_len: u16,
}

/// Build a raw BAM record byte slice (everything AFTER the 4-byte block_size prefix).
///
/// The layout matches what `parse_header` / `push_raw` expects:
/// ```
/// [0..4]   ref_id (i32 LE)
/// [4..8]   pos    (i32 LE)
/// [8]      l_read_name (u8) — includes NUL terminator
/// [9]      mapq   (u8)
/// [10..12] bin    (u16 LE)
/// [12..14] n_cigar_op (u16 LE)
/// [14..16] flag   (u16 LE)
/// [16..20] l_seq  (u32 LE)
/// [20..24] next_ref_id (i32 LE)
/// [24..28] next_pos    (i32 LE)
/// [28..32] tlen        (i32 LE)
/// [32..]   qname (l_read_name bytes, NUL terminated)
///          cigar (n_cigar_op × 4 bytes)
///          seq   (ceil(l_seq/2) bytes, 4-bit packed)
///          qual  (l_seq bytes)
/// ```
fn build_raw_record(rec: &FuzzRecord) -> Vec<u8> {
    // --- CIGAR ---
    // Limit to 8 ops. Build the list of (op, len) pairs.
    let raw_ops: Vec<(u8, u32)> =
        rec.cigar_ops.iter().take(8).map(|op| (op.bam_op(), op.length_clamped())).collect();

    // l_seq = sum of query-consuming op lengths
    let l_seq: u32 = raw_ops
        .iter()
        .filter(|(op, _)| consumes_query(*op))
        .map(|(_, len)| *len)
        .fold(0u32, |acc, v| acc.saturating_add(v));

    let n_cigar_ops = raw_ops.len() as u16;

    // --- Fixed header ---
    let qname = b"fuzz\0"; // 5 bytes, includes NUL
    let l_read_name: u8 = qname.len() as u8;

    // Flags: bit 0x4 (unmapped) cleared; preserve strand/pair bits from flag_bits.
    // Only pass through the meaningful lower bits: paired(1), rev(0x10), first(0x40), second(0x80).
    let flags: u16 = u16::from(rec.flag_bits & 0xD1) & !0x0004;

    let pos_val = i32::from(rec.pos);

    let mut raw = Vec::with_capacity(32 + qname.len() + n_cigar_ops as usize * 4 + l_seq as usize);

    // ref_id = 0 (valid tid for mapped reads)
    raw.extend_from_slice(&0i32.to_le_bytes());
    // pos
    raw.extend_from_slice(&pos_val.to_le_bytes());
    // l_read_name, mapq
    raw.push(l_read_name);
    raw.push(rec.mapq);
    // bin (ignored)
    raw.extend_from_slice(&0u16.to_le_bytes());
    // n_cigar_op
    raw.extend_from_slice(&n_cigar_ops.to_le_bytes());
    // flag
    raw.extend_from_slice(&flags.to_le_bytes());
    // l_seq
    raw.extend_from_slice(&l_seq.to_le_bytes());
    // next_ref_id, next_pos, tlen
    raw.extend_from_slice(&(-1i32).to_le_bytes());
    raw.extend_from_slice(&(-1i32).to_le_bytes());
    raw.extend_from_slice(&0i32.to_le_bytes());

    // --- Variable-length: qname ---
    raw.extend_from_slice(qname);

    // --- Variable-length: CIGAR (packed u32 LE: len<<4 | op) ---
    for (op, len) in &raw_ops {
        let packed: u32 = (len << 4) | u32::from(*op);
        raw.extend_from_slice(&packed.to_le_bytes());
    }

    // --- Variable-length: 4-bit packed seq ---
    // Use provided seq_bases; pad with 0 (nibble 'N' in BAM = 15, but 0 = '=' which
    // decodes to Unknown too). We use nibble 0xF for padding (BAM 'N' = 15).
    let seq_len = l_seq as usize;
    let packed_len = seq_len.div_ceil(2);
    let mut packed_seq = vec![0u8; packed_len];
    for i in 0..seq_len {
        // BAM 4-bit encoding: A=1,C=2,G=4,T=8,N=15 (low nibble = first base in byte)
        // Map arbitrary byte to nibble via modulo over valid nibble values.
        let nibble_table: [u8; 6] = [1, 2, 4, 8, 15, 0]; // A,C,G,T,N,=
        let src = rec.seq_bases.get(i).copied().unwrap_or(15);
        let nibble = nibble_table[(src as usize) % nibble_table.len()];
        if i % 2 == 0 {
            // High nibble
            if let Some(b) = packed_seq.get_mut(i / 2) {
                *b = nibble << 4;
            }
        } else {
            // Low nibble
            if let Some(b) = packed_seq.get_mut(i / 2) {
                *b |= nibble;
            }
        }
    }
    raw.extend_from_slice(&packed_seq);

    // --- Variable-length: qual (l_seq bytes) ---
    for i in 0..seq_len {
        let q = rec.quals.get(i).copied().unwrap_or(0xFF);
        raw.push(q);
    }

    raw
}

fuzz_target!(|input: FuzzPileupInput| {
    const MAX_RECORDS: usize = 32;
    const MAX_REGION_LEN: u32 = 1000;

    let region_start = match Pos0::new(u32::from(input.region_start)) {
        Some(p) => p,
        None => return,
    };
    let region_len = u32::from(input.region_len).min(MAX_REGION_LEN);
    let region_end = match region_start.checked_add_offset(Offset::new(region_len as i64)) {
        Some(p) => p,
        None => return,
    };

    let mut store = RecordStore::new();

    for rec in input.records.iter().take(MAX_RECORDS) {
        let raw = build_raw_record(rec);
        // Errors from push_raw are expected (e.g. zero-length CIGAR → zero ref-span).
        // We intentionally ignore them — the goal is to reach pileup logic.
        let _ = store.push_raw(&raw, &mut ());
    }

    if store.is_empty() {
        return;
    }

    // Sort by position as the pileup engine requires.
    store.sort_by_pos();

    let mut engine = PileupEngine::new(store, region_start, region_end);
    engine.set_max_depth(200);

    let mut columns_seen: u32 = 0;
    while let Some(col) = engine.pileups() {
        let _pos = col.pos();
        let _depth = col.depth();
        let _mdepth = col.match_depth();
        let _refbase = col.reference_base();

        for aln in col.alignments() {
            let _idx = aln.record_idx();
            let _op = aln.op();
            let _qpos = aln.qpos();
            let _base = aln.base();
            let _qual = aln.qual();
            let _is_del = aln.is_del();
            let _is_refskip = aln.is_refskip();
            let _ilen = aln.insert_len();
            let _dlen = aln.del_len();
            let _mapq = aln.mapq;
            let _flags = aln.flags;
            let _strand = aln.strand;
            let _seq_len = aln.seq_len;
            let _matching = aln.matching_bases;
            let _indels = aln.indel_bases;
        }
        columns_seen = columns_seen.saturating_add(1);
        if columns_seen > MAX_REGION_LEN {
            break;
        }
    }
});
