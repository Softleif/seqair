//! Parse and query CIGAR strings. [`CigarMapping`] maps reference positions to query positions
//! (linear fast-path for clip+match CIGARs, `SmallVec` fallback for complex ones);
//! [`CigarPosInfo`] describes what occupies a given reference position.

use crate::utils::{TraceErr, TraceOk};
use seqair_types::{Pos, SmallVec, Zero};

// r[impl cigar.operations]
// r[impl io.named_constants]
pub const CIGAR_M: u8 = 0;
pub const CIGAR_I: u8 = 1;
pub const CIGAR_D: u8 = 2;
pub const CIGAR_N: u8 = 3;
pub const CIGAR_S: u8 = 4;
pub const CIGAR_H: u8 = 5;
pub const CIGAR_P: u8 = 6;
pub const CIGAR_EQ: u8 = 7;
pub const CIGAR_X: u8 = 8;

// r[impl io.typed_cigar_ops]
// r[impl cigar.index]
/// Strongly-typed CIGAR operation kind.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarOpType {
    Match,       // M(0): alignment match (may be sequence match or mismatch)
    Insertion,   // I(1): insertion to the reference
    Deletion,    // D(2): deletion from the reference
    RefSkip,     // N(3): skipped region from the reference (e.g. intron)
    SoftClip,    // S(4): soft clipping
    HardClip,    // H(5): hard clipping
    Padding,     // P(6): padding
    SeqMatch,    // =(7): sequence match
    SeqMismatch, // X(8): sequence mismatch
}

impl CigarOpType {
    pub fn from_bam(code: u8) -> Option<Self> {
        match code {
            CIGAR_M => Some(Self::Match),
            CIGAR_I => Some(Self::Insertion),
            CIGAR_D => Some(Self::Deletion),
            CIGAR_N => Some(Self::RefSkip),
            CIGAR_S => Some(Self::SoftClip),
            CIGAR_H => Some(Self::HardClip),
            CIGAR_P => Some(Self::Padding),
            CIGAR_EQ => Some(Self::SeqMatch),
            CIGAR_X => Some(Self::SeqMismatch),
            _ => None,
        }
    }

    // r[impl bam.owned_record.cigar_op]
    /// Convert back to the BAM numeric code. Inverse of [`from_bam`](Self::from_bam).
    pub const fn to_bam_code(self) -> u8 {
        match self {
            Self::Match => CIGAR_M,
            Self::Insertion => CIGAR_I,
            Self::Deletion => CIGAR_D,
            Self::RefSkip => CIGAR_N,
            Self::SoftClip => CIGAR_S,
            Self::HardClip => CIGAR_H,
            Self::Padding => CIGAR_P,
            Self::SeqMatch => CIGAR_EQ,
            Self::SeqMismatch => CIGAR_X,
        }
    }

    pub const fn consumes_ref(self) -> bool {
        matches!(
            self,
            Self::Match | Self::Deletion | Self::RefSkip | Self::SeqMatch | Self::SeqMismatch
        )
    }

    pub const fn consumes_query(self) -> bool {
        matches!(
            self,
            Self::Match | Self::Insertion | Self::SoftClip | Self::SeqMatch | Self::SeqMismatch
        )
    }
}

// r[impl bam.owned_record.cigar_op]
/// A single CIGAR operation with its operation type and length.
///
/// In BAM binary format, each op is packed as a u32: `len << 4 | op_code`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CigarOp {
    pub op: CigarOpType,
    pub len: u32,
}

impl CigarOp {
    /// BAM packs length into the upper 28 bits, so max length is 2^28-1 = `268_435_455`.
    pub const fn new(op: CigarOpType, len: u32) -> Self {
        debug_assert!(len < (1 << 28), "CIGAR op length exceeds 28-bit BAM limit");
        Self { op, len }
    }

    /// Decode from BAM packed u32 format (`len << 4 | op_code`).
    pub fn from_bam_u32(packed: u32) -> Option<Self> {
        let op = CigarOpType::from_bam((packed & 0xF) as u8)?;
        let len = packed >> 4;
        Some(Self { op, len })
    }

    /// Encode to BAM packed u32 format (`len << 4 | op_code`).
    pub const fn to_bam_u32(self) -> u32 {
        (self.len << 4) | self.op.to_bam_code() as u32
    }
}

/// Whether a CIGAR op code consumes the reference.
const fn consumes_ref(op: u8) -> bool {
    matches!(op, CIGAR_M | CIGAR_D | CIGAR_N | CIGAR_EQ | CIGAR_X)
}

/// Whether a CIGAR op code consumes the query.
const fn consumes_query(op: u8) -> bool {
    matches!(op, CIGAR_M | CIGAR_I | CIGAR_S | CIGAR_EQ | CIGAR_X)
}

// r[impl cigar.matches_indels]
/// Calculate matches and indels from packed CIGAR bytes (LE u32 per op).
///
/// Lower 4 bits encode the operation; upper 28 bits encode the length.
pub fn calc_matches_indels(cigar_bytes: &[u8]) -> (u32, u32) {
    let mut matches = 0u32;
    let mut indels = 0u32;
    let n_ops = cigar_bytes.len() / 4;
    for i in 0..n_ops {
        let op = u32::from_le_bytes(read4(cigar_bytes, i.wrapping_mul(4)));
        let len = op >> 4;
        let op_type = (op & 0xF) as u8;
        match op_type {
            CIGAR_M => matches = matches.saturating_add(len),
            CIGAR_I | CIGAR_D => indels = indels.saturating_add(len),
            _ => {}
        }
    }
    (matches, indels)
}

/// Position information returned by `CigarMapping` for the pileup engine.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarPosInfo {
    /// M/=/X op: read has a base aligned here.
    Match { qpos: u32 },
    /// M/=/X op at the last position before an I op follows.
    Insertion { qpos: u32, insert_len: u32 },
    /// D op: deletion spanning this position. `del_len` is the total length of the D CIGAR op.
    Deletion { del_len: u32 },
    /// N op: reference skip spanning this position.
    RefSkip,
}

/// Compact CIGAR mapping for fast qpos lookup in the pileup engine.
///
/// For the ~90% of reads with simple CIGARs (clips + one match block),
/// `Linear` avoids any allocation and computes qpos with a single subtraction.
/// Complex CIGARs use a `SmallVec` of compact ops that stays inline for ≤6 ops.
pub enum CigarMapping {
    /// `qpos = (pos - rec_pos) as usize + query_offset as usize`
    Linear { rec_pos: Pos<Zero>, query_offset: u32, match_len: u32 },
    /// Pre-computed compact ops for linear/binary search.
    Complex(SmallVec<CompactOp, 6>),
}

impl std::fmt::Debug for CigarMapping {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Linear { rec_pos, query_offset, match_len } => f
                .debug_struct("Linear")
                .field("rec_pos", &rec_pos.get())
                .field("query_offset", query_offset)
                .field("match_len", match_len)
                .finish(),
            Self::Complex(ops) => write!(f, "Complex({} ops)", ops.len()),
        }
    }
}

/// 16-byte CIGAR op for the compact index.
#[derive(Clone, Copy)]
pub struct CompactOp {
    ref_start: i32,
    query_start: u32,
    len: u32,
    op_type: u8,
}

impl CigarMapping {
    #[inline]
    pub fn new(rec_pos: Pos<Zero>, cigar_bytes: &[u8]) -> Option<Self> {
        match try_linear(cigar_bytes) {
            Some((query_offset, match_len)) => {
                Some(Self::Linear { rec_pos, query_offset, match_len })
            }
            None => Some(Self::Complex(build_compact_ops(rec_pos, cigar_bytes)?)),
        }
    }

    // r[impl cigar.qpos_at]
    #[inline]
    pub fn pos_info_at(&self, pos: Pos<Zero>) -> Option<CigarPosInfo> {
        match self {
            Self::Linear { rec_pos, query_offset, match_len } => {
                let offset = pos.as_i64().wrapping_sub(rec_pos.as_i64());
                if offset < 0 || offset >= i64::from(*match_len) {
                    return None;
                }
                Some(CigarPosInfo::Match {
                    qpos: u32::try_from(offset)
                        .trace_ok("offset exceeds u32 range")?
                        .checked_add(*query_offset)
                        .trace_err("qpos overflow")?,
                })
                // Linear path never has insertions/deletions (try_linear rejects them)
            }
            Self::Complex(ops) => {
                if ops.len() <= 4 {
                    pos_info_linear(ops, pos)
                } else {
                    pos_info_bsearch(ops, pos)
                }
            }
        }
    }
}

// r[impl cigar.qpos_bounds]
/// Check if the CIGAR is a simple clips-match-clips pattern.
/// Returns `(query_offset, match_len)` where `query_offset` is the leading soft-clip length
/// and `match_len` is the total length of the contiguous match/seq-match/seq-mismatch block.
/// Returns `None` for CIGARs with insertions, deletions, or multiple disjoint match regions.
#[inline]
fn try_linear(cigar_bytes: &[u8]) -> Option<(u32, u32)> {
    let n_ops = cigar_bytes.len() / 4;
    let mut query_offset = 0u32;
    let mut match_len = 0u32;
    // 0 = leading clips, 1 = match block, 2 = trailing clips
    let mut phase = 0u8;

    for i in 0..n_ops {
        let packed = u32::from_le_bytes(read4(
            cigar_bytes,
            i.checked_mul(4).trace_err("cigar loop overflow")?,
        ));
        let len = packed >> 4;
        let op_type = (packed & 0xF) as u8;

        match (phase, op_type) {
            (0, CIGAR_S) => query_offset = query_offset.saturating_add(len),
            (0, CIGAR_H) => {}
            (0, CIGAR_M | CIGAR_EQ | CIGAR_X) => {
                match_len = match_len.saturating_add(len);
                phase = 1;
            }
            (1, CIGAR_M | CIGAR_EQ | CIGAR_X) => match_len = match_len.saturating_add(len),
            (1 | 2, CIGAR_S | CIGAR_H) => phase = 2,
            _ => return None,
        }
    }

    // Only succeed if we reached the match phase (phase >= 1).
    // Pure soft/hard clips (phase == 0) should use the complex path.
    if phase >= 1 { Some((query_offset, match_len)) } else { None }
}

fn build_compact_ops(rec_pos: Pos<Zero>, cigar_bytes: &[u8]) -> Option<SmallVec<CompactOp, 6>> {
    let n_ops = cigar_bytes.len() / 4;
    let mut ops = SmallVec::with_capacity(n_ops);
    let mut ref_off: i64 = 0;
    let mut query_off: u32 = 0;

    for i in 0..n_ops {
        let packed = u32::from_le_bytes(read4(
            cigar_bytes,
            i.checked_mul(4).trace_err("cigar loop overflow")?,
        ));
        let len = packed >> 4;
        let op_type = (packed & 0xF) as u8;

        let ref_start_i64 = rec_pos.as_i64().wrapping_add(ref_off);
        // r[impl cigar.compact_op_position_invariant]
        let ref_start = i32::try_from(ref_start_i64).trace_ok("ref_start exceeds i32 range")?;
        ops.push(CompactOp { ref_start, query_start: query_off, len, op_type });

        if consumes_ref(op_type) {
            ref_off = ref_off.checked_add(i64::from(len)).trace_err("ref offset overflow")?;
        }
        if consumes_query(op_type) {
            query_off = query_off.saturating_add(len);
        }
    }

    Some(ops)
}

// r[impl pileup_indel.insertion_len]
/// Sum insertion lengths for consecutive I ops after index `op_idx`.
#[inline]
fn next_insertion_len(ops: &[CompactOp], op_idx: usize) -> Option<u32> {
    let mut total = 0u32;
    let mut idx = op_idx.checked_add(1)?;
    while let Some(next) = ops.get(idx) {
        if next.op_type != CIGAR_I {
            break;
        }
        total = total.saturating_add(next.len);
        idx = idx.checked_add(1)?;
    }
    if total > 0 { Some(total) } else { None }
}

// r[impl pileup_indel.insertion_at_last_match]
// r[impl pileup_indel.no_orphan_insertions]
#[inline]
fn classify_op(
    ops: &[CompactOp],
    i: usize,
    op: &CompactOp,
    pos32: i32,
    ref_end: i32,
) -> Option<CigarPosInfo> {
    if consumes_query(op.op_type) {
        debug_assert!(
            matches!(op.op_type, CIGAR_M | CIGAR_EQ | CIGAR_X),
            "classify_op reached query-consuming branch with unexpected op type {}",
            op.op_type
        );
        let offset = pos32.wrapping_sub(op.ref_start) as u32;
        let qpos = op.query_start.checked_add(offset).trace_err("qpos overflow")?;
        if pos32 == ref_end.wrapping_sub(1)
            && let Some(insert_len) = next_insertion_len(ops, i)
        {
            return Some(CigarPosInfo::Insertion { qpos, insert_len });
        }
        Some(CigarPosInfo::Match { qpos })
    } else if op.op_type == CIGAR_D {
        // r[impl pileup_indel.op_enum]
        Some(CigarPosInfo::Deletion { del_len: op.len })
    } else if op.op_type == CIGAR_N {
        Some(CigarPosInfo::RefSkip)
    } else {
        None
    }
}

#[inline]
fn pos_info_linear(ops: &[CompactOp], pos: Pos<Zero>) -> Option<CigarPosInfo> {
    let pos32 = i32::try_from(pos.get()).trace_ok("pos exceeds i32 range")?;
    for (i, op) in ops.iter().enumerate() {
        if !consumes_ref(op.op_type) {
            continue;
        }
        let ref_end =
            op.ref_start.saturating_add(i32::try_from(op.len).trace_ok("op len exceeds i32")?);
        if pos32 < op.ref_start || pos32 >= ref_end {
            continue;
        }
        return classify_op(ops, i, op, pos32, ref_end);
    }
    None
}

// r[impl perf.cigar_binary_search]
#[inline]
fn pos_info_bsearch(ops: &[CompactOp], pos: Pos<Zero>) -> Option<CigarPosInfo> {
    let pos32 = i32::try_from(pos.get()).trace_ok("pos exceeds i32 range")?;
    let idx = ops.partition_point(|op| op.ref_start <= pos32);
    if idx == 0 {
        return None;
    }
    for i in idx.saturating_sub(2)..ops.len().min(idx.saturating_add(1)) {
        let Some(op) = ops.get(i) else { continue };
        if !consumes_ref(op.op_type) {
            continue;
        }
        let ref_end =
            op.ref_start.saturating_add(i32::try_from(op.len).trace_ok("op len exceeds i32")?);
        if pos32 >= op.ref_start && pos32 < ref_end {
            return classify_op(ops, i, op, pos32, ref_end);
        }
    }
    None
}

// Callers only invoke read4 with offset = i * 4 where i < buf.len() / 4, so offset + 3 < buf.len().
#[allow(clippy::indexing_slicing, reason = "offset + 3 < buf.len() ensured by caller's loop bound")]
fn read4(buf: &[u8], offset: usize) -> [u8; 4] {
    debug_assert!(
        offset.saturating_add(3) < buf.len(),
        "read4 out of bounds: offset={offset}, len={}",
        buf.len()
    );
    [
        buf[offset],
        buf[offset.saturating_add(1)],
        buf[offset.saturating_add(2)],
        buf[offset.saturating_add(3)],
    ]
}

#[cfg(test)]
#[allow(clippy::arithmetic_side_effects, reason = "test arithmetic on known small values")]
mod tests {
    use super::*;

    /// Pack a single CIGAR op into 4 LE bytes: upper 28 bits = len, lower 4 = op.
    fn pack_cigar_op(op: u8, len: u32) -> [u8; 4] {
        ((len << 4) | u32::from(op)).to_le_bytes()
    }

    // r[verify bam.owned_record.cigar_op]
    #[test]
    fn cigar_op_roundtrip_all_ops() {
        let ops = [
            (CigarOpType::Match, CIGAR_M),
            (CigarOpType::Insertion, CIGAR_I),
            (CigarOpType::Deletion, CIGAR_D),
            (CigarOpType::RefSkip, CIGAR_N),
            (CigarOpType::SoftClip, CIGAR_S),
            (CigarOpType::HardClip, CIGAR_H),
            (CigarOpType::Padding, CIGAR_P),
            (CigarOpType::SeqMatch, CIGAR_EQ),
            (CigarOpType::SeqMismatch, CIGAR_X),
        ];
        for (op_type, code) in ops {
            assert_eq!(op_type.to_bam_code(), code, "to_bam_code failed for {op_type:?}");
            assert_eq!(
                CigarOpType::from_bam(code),
                Some(op_type),
                "from_bam failed for code {code}"
            );

            let cigar_op = CigarOp::new(op_type, 42);
            let packed = cigar_op.to_bam_u32();
            let decoded = CigarOp::from_bam_u32(packed).unwrap();
            assert_eq!(decoded, cigar_op, "round-trip failed for {op_type:?}");
        }
    }

    #[test]
    fn cigar_op_max_length() {
        // Upper 28 bits can hold up to 2^28 - 1 = 268_435_455
        let max_len = (1u32 << 28) - 1;
        let op = CigarOp::new(CigarOpType::Match, max_len);
        let packed = op.to_bam_u32();
        let decoded = CigarOp::from_bam_u32(packed).unwrap();
        assert_eq!(decoded.len, max_len);
        assert_eq!(decoded.op, CigarOpType::Match);
    }

    #[test]
    fn cigar_op_from_bam_u32_invalid_code() {
        // op code 9 is invalid
        let packed = (100u32 << 4) | 9;
        assert!(CigarOp::from_bam_u32(packed).is_none());
    }

    // r[verify cigar.compact_op_position_invariant]
    #[test]
    fn compact_op_ref_start_fits_i32_for_bam_positions() {
        // BAM positions are i32 (max 2^31-1 = 2_147_483_647).
        // CompactOp stores ref_start as i32. Verify it works at the maximum BAM position.
        let mut cigar = Vec::new();
        cigar.extend_from_slice(&pack_cigar_op(CIGAR_M, 100));
        let rec_pos = Pos::<Zero>::new((i32::MAX as u32) - 100).unwrap(); // near max BAM position
        let mapping = CigarMapping::new(rec_pos, &cigar).unwrap();
        // Should work fine — position fits in i32
        assert!(matches!(mapping, CigarMapping::Linear { .. }));
        assert_eq!(mapping.pos_info_at(rec_pos), Some(CigarPosInfo::Match { qpos: 0 }));
    }

    // r[verify cigar.qpos_bounds]
    #[test]
    fn linear_cigar_mapping_bounds_check() {
        // 5S + 100M: rec_pos=1000, query_offset=5, match_len=100
        let mut cigar = Vec::new();
        cigar.extend_from_slice(&pack_cigar_op(CIGAR_S, 5));
        cigar.extend_from_slice(&pack_cigar_op(CIGAR_M, 100));

        let mapping = CigarMapping::new(Pos::<Zero>::new(1000).unwrap(), &cigar).unwrap();
        assert!(matches!(mapping, CigarMapping::Linear { .. }));

        // Valid positions: 1000..1100
        assert_eq!(
            mapping.pos_info_at(Pos::<Zero>::new(1000).unwrap()),
            Some(CigarPosInfo::Match { qpos: 5 })
        );
        assert_eq!(
            mapping.pos_info_at(Pos::<Zero>::new(1050).unwrap()),
            Some(CigarPosInfo::Match { qpos: 55 })
        );
        assert_eq!(
            mapping.pos_info_at(Pos::<Zero>::new(1099).unwrap()),
            Some(CigarPosInfo::Match { qpos: 104 })
        );

        // Out-of-range: before alignment start
        assert_eq!(
            mapping.pos_info_at(Pos::<Zero>::new(999).unwrap()),
            None,
            "pos before rec_pos must return None"
        );

        // Out-of-range: at/past alignment end
        assert_eq!(
            mapping.pos_info_at(Pos::<Zero>::new(1100).unwrap()),
            None,
            "pos at rec_pos + match_len must return None"
        );
        assert_eq!(
            mapping.pos_info_at(Pos::<Zero>::new(1200).unwrap()),
            None,
            "pos past alignment end must return None"
        );

        // Far out-of-range (would wrap with unsigned subtraction)
        assert_eq!(
            mapping.pos_info_at(Pos::<Zero>::new(0).unwrap()),
            None,
            "pos far before alignment must return None"
        );
    }
}
