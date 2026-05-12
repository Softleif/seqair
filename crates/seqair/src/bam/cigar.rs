//! Parse and query CIGAR strings. [`CigarMapping`] maps reference positions to query positions
//! (linear fast-path for clip+match CIGARs, `SmallVec` fallback for complex ones);
//! [`CigarPosInfo`] describes what occupies a given reference position.

use crate::utils::{TraceErr, TraceOk};
use seqair_types::{Pos0, SmallVec};

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
/// Strongly-typed CIGAR operation kind. `Unknown` carries the raw 4-bit
/// op code from BAM for codes that aren't in the SAM spec — we tolerate
/// them on read so a single bad op doesn't fail the whole record, and
/// round-trip them verbatim on write.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarOpType {
    /// M(0): alignment match (may be sequence match or mismatch)
    Match,
    /// I(1): insertion to the reference
    Insertion,
    /// D(2): deletion from the reference
    Deletion,
    /// N(3): skipped region from the reference (e.g. intron)
    RefSkip,
    /// S(4): soft clipping
    SoftClip,
    /// H(5): hard clipping
    HardClip,
    /// P(6): padding
    Padding,
    /// =(7): sequence match
    SeqMatch,
    /// X(8): sequence mismatch
    SeqMismatch,
    /// 9..=15: reserved/unspecified code from BAM; the inner byte is the raw 4-bit op
    Unknown(u8),
}

impl CigarOpType {
    /// Decode a 4-bit BAM op code. Codes outside the spec become [`Self::Unknown`].
    pub const fn from_bam(code: u8) -> Self {
        match code {
            CIGAR_M => Self::Match,
            CIGAR_I => Self::Insertion,
            CIGAR_D => Self::Deletion,
            CIGAR_N => Self::RefSkip,
            CIGAR_S => Self::SoftClip,
            CIGAR_H => Self::HardClip,
            CIGAR_P => Self::Padding,
            CIGAR_EQ => Self::SeqMatch,
            CIGAR_X => Self::SeqMismatch,
            other => Self::Unknown(other),
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
            Self::Unknown(code) => code,
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
/// A single CIGAR operation, stored in the same packed u32 layout BAM uses
/// on disk: `len << 4 | op_code`. Decoding the op type is a single shift+mask;
/// re-encoding for BAM output is a no-op.
#[repr(transparent)]
#[derive(Clone, Copy, PartialEq, Eq, Hash, bytemuck::Pod, bytemuck::Zeroable)]
pub struct CigarOp(u32);

const _: () = assert!(size_of::<CigarOp>() == 4, "CigarOp must stay 4 bytes");

impl CigarOp {
    /// BAM packs length into the upper 28 bits, so max length is 2^28-1 = `268_435_455`.
    pub const fn new(op: CigarOpType, len: u32) -> Self {
        debug_assert!(len < (1 << 28), "CIGAR op length exceeds 28-bit BAM limit");
        Self((len << 4) | op.to_bam_code() as u32)
    }

    /// Wrap a raw BAM-packed u32 (`len << 4 | op_code`). Always succeeds —
    /// unknown 4-bit codes surface as [`CigarOpType::Unknown`] when decoded.
    pub const fn from_bam_u32(packed: u32) -> Self {
        Self(packed)
    }

    /// Encode to BAM packed u32 format (`len << 4 | op_code`).
    pub const fn to_bam_u32(self) -> u32 {
        self.0
    }

    /// Length of the operation (upper 28 bits).
    #[allow(clippy::len_without_is_empty, reason = "len is the CIGAR op span, not a collection")]
    pub const fn len(self) -> u32 {
        self.0 >> 4
    }

    /// Raw BAM op code (lower 4 bits).
    pub const fn op_code(self) -> u8 {
        (self.0 & 0xF) as u8
    }

    /// Decoded op type. Reserved op codes (9..=15) decode to [`CigarOpType::Unknown`].
    pub const fn op_type(self) -> CigarOpType {
        CigarOpType::from_bam(self.op_code())
    }

    pub const fn consumes_ref(self) -> bool {
        consumes_ref(self.op_code())
    }

    pub const fn consumes_query(self) -> bool {
        consumes_query(self.op_code())
    }

    /// Reinterpret a slice of ops as their on-the-wire BAM bytes
    /// (LE u32 per op). Safe via `CigarOp: Pod` — `repr(transparent)` over
    /// `u32` with every bit pattern a valid op.
    ///
    /// On big-endian hosts the bytes are still in native order, NOT BAM-on-disk
    /// order. Callers that write to disk on a BE host would need to byte-swap
    /// each op first; we don't, because seqair targets LE only.
    #[inline]
    pub fn ops_as_bytes(ops: &[Self]) -> &[u8] {
        bytemuck::cast_slice(ops)
    }

    /// Append BAM-on-disk CIGAR bytes (LE u32 per op) into a typed `Vec<CigarOp>`.
    ///
    /// Source bytes are unaligned; on LE hosts this is a single memcpy, on BE
    /// it would byte-swap per op (we don't compile that path — seqair is LE only).
    /// Caller must ensure `bytes.len()` is a multiple of 4 (true for valid BAM).
    #[inline]
    pub fn extend_from_bam_bytes(dst: &mut Vec<Self>, bytes: &[u8]) {
        debug_assert!(bytes.len().is_multiple_of(4), "BAM CIGAR bytes must be multiple of 4");
        #[cfg(target_endian = "little")]
        {
            let n_ops = bytes.len() / 4;
            let n_bytes = n_ops.checked_mul(4).expect("cigar byte length overflow");
            let new_len = dst.len().checked_add(n_ops).expect("cigar slab len overflow");
            dst.reserve(n_ops);
            // SAFETY: CigarOp is repr(transparent) over u32 (4 bytes) with a
            // compile-time size guard on line 107. Every bit pattern is a valid
            // CigarOp (op codes 0..=15 all decodable, length field is unconstrained
            // u28). On LE the in-memory u32 byte order matches BAM's LE on-disk
            // order, so a bytewise memcpy of n_ops * 4 source bytes produces n_ops
            // valid CigarOp values. dst has at least n_ops uninit slots after
            // reserve(); we copy n_bytes into it, then bump len.
            unsafe {
                let dst_ptr = dst.as_mut_ptr().add(dst.len()).cast::<u8>();
                std::ptr::copy_nonoverlapping(bytes.as_ptr(), dst_ptr, n_bytes);
                dst.set_len(new_len);
            }
        }
        #[cfg(not(target_endian = "little"))]
        {
            for chunk in bytes.chunks_exact(4) {
                let arr: [u8; 4] = chunk.try_into().expect("chunks_exact(4) yields 4 bytes");
                dst.push(CigarOp::from_bam_u32(u32::from_le_bytes(arr)));
            }
        }
    }

    /// View raw BAM CIGAR bytes as `&[CigarOp]` without copying.
    ///
    /// Returns `None` if the length is not a multiple of 4, the pointer is
    /// not 4-byte aligned, or the target is big-endian. On aligned LE input
    /// this is a zero-cost transmute. BAM CIGAR data is **not** guaranteed
    /// aligned — the read name length determines the start offset, and
    /// samtools produces unaligned CIGAR data when `l_read_name` is not a
    /// multiple of 4. Use [`extend_from_bam_bytes`] if alignment cannot be
    /// guaranteed (it's always correct, just requires a copy).
    #[inline]
    pub fn slice_from_bam_bytes(bytes: &[u8]) -> Option<&[Self]> {
        if bytes.is_empty() {
            return Some(&[]);
        }
        if !bytes.len().is_multiple_of(4) {
            return None;
        }
        #[cfg(target_endian = "little")]
        {
            // BAM data is packed sequentially: 32-byte fixed header, then
            // qname (l_read_name bytes, including NUL), then CIGAR. If
            // l_read_name is not a multiple of 4, the CIGAR pointer is
            // unaligned. from_raw_parts requires CigarOp alignment (4).
            #[allow(
                clippy::arithmetic_side_effects,
                clippy::manual_is_multiple_of,
                reason = "pointer-to-usize cast is defined on all supported targets; align_of is non-zero so % is safe"
            )]
            if bytes.as_ptr() as usize % align_of::<CigarOp>() != 0 {
                return None;
            }
            let n = bytes.len() / 4;
            // SAFETY:
            // - Length: `n` is `bytes.len() / 4`, checked to be a multiple of
            //   4 above. The pointer advance of `n * size_of::<CigarOp>()`
            //   (== n * 4) stays within the allocation.
            // - Alignment: checked via `is_aligned()` above.
            // - Validity: every `u32` bit pattern is a valid `CigarOp` — the
            //   op code field is 4 bits (0..=15, all decodable) and the length
            //   field is an unconstrained 28-bit integer.
            // - Endianness: guarded by `#[cfg(target_endian = "little")]`;
            //   on LE hosts BAM's on-disk LE `u32` byte order matches native
            //   in-memory layout.
            Some(unsafe { std::slice::from_raw_parts(bytes.as_ptr() as *const CigarOp, n) })
        }
        #[cfg(not(target_endian = "little"))]
        {
            None
        }
    }
}

impl std::fmt::Debug for CigarOp {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("CigarOp").field("op", &self.op_type()).field("len", &self.len()).finish()
    }
}

impl std::fmt::Display for CigarOpType {
    /// SAM 1.6 single-character op code. `Unknown` renders as `?`.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        use std::fmt::Write as _;
        let c = match self {
            Self::Match => 'M',
            Self::Insertion => 'I',
            Self::Deletion => 'D',
            Self::RefSkip => 'N',
            Self::SoftClip => 'S',
            Self::HardClip => 'H',
            Self::Padding => 'P',
            Self::SeqMatch => '=',
            Self::SeqMismatch => 'X',
            Self::Unknown(_) => '?',
        };
        f.write_char(c)
    }
}

impl std::fmt::Display for CigarOp {
    /// `{len}{op_char}`, e.g. `10M`, `2I`, `4=`.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}{}", self.len(), self.op_type())
    }
}

/// Display wrapper for a slice of [`CigarOp`]s, formatting as a SAM CIGAR
/// string (e.g. `"10M2I3D"`). Allocation-free — `format!("{}", CigarStr(ops))`
/// when you need an owned `String`.
pub struct CigarStr<'a>(pub &'a [CigarOp]);

impl std::fmt::Display for CigarStr<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for op in self.0 {
            write!(f, "{op}")?;
        }
        Ok(())
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

/// Sum of query-consuming op lengths (M, I, S, =, X).
pub fn calc_query_len(ops: &[CigarOp]) -> u32 {
    let mut qlen = 0u32;
    for op in ops {
        if op.consumes_query() {
            qlen = qlen.saturating_add(op.len());
        }
    }
    qlen
}

// r[impl cigar.matches_indels]
/// Sum of M op lengths and I/D op lengths.
pub fn calc_matches_indels(ops: &[CigarOp]) -> (u32, u32) {
    let mut matches = 0u32;
    let mut indels = 0u32;
    for op in ops {
        match op.op_code() {
            CIGAR_M => matches = matches.saturating_add(op.len()),
            CIGAR_I | CIGAR_D => indels = indels.saturating_add(op.len()),
            _ => {}
        }
    }
    (matches, indels)
}

// r[impl bam.record.end_pos]
// r[impl bam.record.zero_refspan]
/// 0-based exclusive end position from `pos` + reference-consuming op lengths.
pub fn compute_end_pos(pos: Pos0, ops: &[CigarOp]) -> Option<Pos0> {
    let mut ref_len: i64 = 0;
    for op in ops {
        if op.consumes_ref() {
            ref_len = ref_len.checked_add(i64::from(op.len()))?;
        }
    }
    if ref_len == 0 {
        Some(pos)
    } else {
        pos.checked_add_offset(seqair_types::Offset::new(ref_len.checked_sub(1)?))
    }
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
    /// D or N op at its last position, followed by an I op (e.g. `D I M`).
    /// `del_len` is the total D/N op length; `insert_len` is the total following insertion length.
    // r[impl pileup_indel.complex_indel]
    ComplexIndel { del_len: u32, insert_len: u32, is_refskip: bool },
    /// N op: reference skip spanning this position.
    RefSkip,
}

/// Compact CIGAR mapping for fast qpos lookup in the pileup engine.
///
/// Two variants:
///
/// **`Linear`**: For the ~96% of reads whose CIGAR is a single contiguous
/// match block with optional leading/trailing clips (`[HS]* [M=X]+ [HS]*`).
/// Computes `qpos` with a single subtraction. 16 bytes.
///
/// **`Complex`**: For reads with insertions, deletions, reference skips,
/// or multiple disjoint match blocks (~4% of reads). Stores a pre-computed
/// index of [`CompactOp`]s — 12 bytes each, 6 inline (72-byte buffer, no
/// heap allocation for 98.1% of complex CIGARs).
///
/// Real-world data from NA12878 chr12 Illumina WGS (100k reads):
/// 96.1% Linear, 3.9% Complex. Of complex reads: 79% have 3 ops (e.g. `M I M`),
/// 98.1% fit in ≤6 ops (inline).
pub enum CigarMapping {
    /// `qpos = (pos - rec_pos) + query_offset`. No insertion/deletion support —
    /// `try_linear` rejects CIGARs containing I, D, or N ops.
    Linear { rec_pos: Pos0, query_offset: u32, match_len: u32 },
    /// Pre-computed compact ops for position lookup via linear scan (≤4 ops)
    /// or binary search (>4 ops). The mapping is immutable — `pos_info_at`
    /// is a pure function of `pos`, allowing the active set to be reordered
    /// without invalidating state.
    Complex(SmallVec<CompactOp, 6>),
}

impl std::fmt::Debug for CigarMapping {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Linear { rec_pos, query_offset, match_len } => f
                .debug_struct("Linear")
                .field("rec_pos", &rec_pos)
                .field("query_offset", query_offset)
                .field("match_len", match_len)
                .finish(),
            Self::Complex(ops) => write!(f, "Complex({} ops)", ops.len()),
        }
    }
}

/// Compact CIGAR operation for the position-lookup index used by the pileup engine.
///
/// 12 bytes (down from 16). Each op records the absolute reference start, absolute
/// query start, operation length, and CIGAR op code. The length and op code are
/// packed into a single `u32` (`len << 4 | op_code`), matching BAM's on-disk
/// packing — extraction is a single shift or mask.
///
/// Stored inline in [`CigarMapping::Complex`] via `SmallVec<CompactOp, 6>`
/// (72-byte inline buffer, no heap allocation for 98.1% of complex CIGARs).
#[derive(Clone, Copy)]
pub struct CompactOp {
    /// 0-based reference position where this operation begins.
    /// Range: `0..=i32::MAX` (BAM's position limit).
    ref_start: u32,
    /// 0-based query position where this operation begins. Equal to the
    /// cumulative leading clip length plus all preceding query-consuming
    /// op lengths.
    query_start: u32,
    /// Packed: `len << 4 | op_code`. Same layout as BAM's on-disk CIGAR
    /// encoding. Upper 28 bits = operation length, lower 4 bits = CIGAR
    /// op code (0=M, 1=I, 2=D, 3=N, 4=S, 5=H, 6=P, 7==, 8=X).
    len_op: u32,
}

impl CompactOp {
    /// CIGAR operation length (upper 28 bits of the BAM-packed `u32`).
    #[inline]
    #[allow(
        clippy::len_without_is_empty,
        reason = "CIGAR ops don't have a meaningful 'empty' state"
    )]
    pub fn len(self) -> u32 {
        self.len_op >> 4
    }
    /// BAM CIGAR op code (lower 4 bits). See module-level constants
    /// `CIGAR_M` through `CIGAR_X`.
    #[inline]
    pub fn op_type(self) -> u8 {
        (self.len_op & 0xF) as u8
    }
}

impl CigarMapping {
    #[inline]
    pub fn new(rec_pos: Pos0, ops: &[CigarOp]) -> Option<Self> {
        match try_linear(ops) {
            Some((query_offset, match_len)) => {
                Some(Self::Linear { rec_pos, query_offset, match_len })
            }
            None => Some(Self::Complex(build_compact_ops(rec_pos, ops)?)),
        }
    }

    // r[impl cigar.qpos_at]
    #[inline]
    pub fn pos_info_at(&self, pos: Pos0) -> Option<CigarPosInfo> {
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
fn try_linear(ops: &[CigarOp]) -> Option<(u32, u32)> {
    let mut query_offset = 0u32;
    let mut match_len = 0u32;
    // 0 = leading clips, 1 = match block, 2 = trailing clips
    let mut phase = 0u8;

    for op in ops {
        let len = op.len();
        match (phase, op.op_code()) {
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

fn build_compact_ops(rec_pos: Pos0, ops: &[CigarOp]) -> Option<SmallVec<CompactOp, 6>> {
    let mut out = SmallVec::with_capacity(ops.len());
    let mut ref_off: i64 = 0;
    let mut query_off: u32 = 0;

    for op in ops {
        let len = op.len();
        let op_type = op.op_code();

        let ref_start_i64 = rec_pos.as_i64().wrapping_add(ref_off);
        // r[impl cigar.compact_op_position_invariant]
        let ref_start = u32::try_from(ref_start_i64).trace_ok("ref_start exceeds u32 range")?;
        let len_op = op.to_bam_u32();
        out.push(CompactOp { ref_start, query_start: query_off, len_op });

        if consumes_ref(op_type) {
            ref_off = ref_off.checked_add(i64::from(len)).trace_err("ref offset overflow")?;
        }
        if consumes_query(op_type) {
            query_off = query_off.saturating_add(len);
        }
    }

    Some(out)
}

// r[impl pileup_indel.insertion_len]
// r[impl pileup_indel.complex_indel]
/// Sum insertion lengths for consecutive I ops after index `op_idx`,
/// skipping P (padding) ops. htslib skips P when peeking ahead for
/// insertions (handles `M P I`, `D P I`, etc.).
#[inline]
fn next_insertion_len(ops: &[CompactOp], op_idx: usize) -> Option<u32> {
    let mut total = 0u32;
    let mut idx = op_idx.checked_add(1)?;
    while let Some(next) = ops.get(idx) {
        if next.op_type() == CIGAR_P {
            idx = idx.checked_add(1)?;
            continue;
        }
        if next.op_type() != CIGAR_I {
            break;
        }
        total = total.saturating_add(next.len());
        idx = idx.checked_add(1)?;
    }
    if total > 0 { Some(total) } else { None }
}

// r[impl pileup_indel.insertion_at_last_match]
// r[impl pileup_indel.complex_indel]
#[inline]
fn classify_op(
    ops: &[CompactOp],
    i: usize,
    op: &CompactOp,
    pos: u32,
    ref_end: u32,
) -> Option<CigarPosInfo> {
    if consumes_query(op.op_type()) {
        debug_assert!(
            matches!(op.op_type(), CIGAR_M | CIGAR_EQ | CIGAR_X),
            "classify_op reached query-consuming branch with unexpected op type {}",
            op.op_type()
        );
        let offset = pos.wrapping_sub(op.ref_start);
        let qpos = op.query_start.checked_add(offset).trace_err("qpos overflow")?;
        if pos == ref_end.wrapping_sub(1)
            && let Some(insert_len) = next_insertion_len(ops, i)
        {
            return Some(CigarPosInfo::Insertion { qpos, insert_len });
        }
        Some(CigarPosInfo::Match { qpos })
    } else if op.op_type() == CIGAR_D {
        if pos == ref_end.wrapping_sub(1)
            && let Some(insert_len) = next_insertion_len(ops, i)
        {
            return Some(CigarPosInfo::ComplexIndel {
                del_len: op.len(),
                insert_len,
                is_refskip: false,
            });
        }
        Some(CigarPosInfo::Deletion { del_len: op.len() })
    } else if op.op_type() == CIGAR_N {
        if pos == ref_end.wrapping_sub(1)
            && let Some(insert_len) = next_insertion_len(ops, i)
        {
            return Some(CigarPosInfo::ComplexIndel {
                del_len: op.len(),
                insert_len,
                is_refskip: true,
            });
        }
        Some(CigarPosInfo::RefSkip)
    } else {
        None
    }
}

#[inline]
fn pos_info_linear(ops: &[CompactOp], pos: Pos0) -> Option<CigarPosInfo> {
    let pos: u32 = *pos;
    for (i, op) in ops.iter().enumerate() {
        if !consumes_ref(op.op_type()) {
            continue;
        }
        let ref_end = op.ref_start.saturating_add(op.len());
        if pos < op.ref_start || pos >= ref_end {
            continue;
        }
        return classify_op(ops, i, op, pos, ref_end);
    }
    None
}

// r[impl perf.cigar_binary_search]
#[inline]
fn pos_info_bsearch(ops: &[CompactOp], pos: Pos0) -> Option<CigarPosInfo> {
    let pos: u32 = *pos;
    let idx = ops.partition_point(|op| op.ref_start <= pos);
    if idx == 0 {
        return None;
    }
    for i in idx.saturating_sub(2)..ops.len().min(idx.saturating_add(1)) {
        let Some(op) = ops.get(i) else { continue };
        if !consumes_ref(op.op_type()) {
            continue;
        }
        let ref_end = op.ref_start.saturating_add(op.len());
        if pos >= op.ref_start && pos < ref_end {
            return classify_op(ops, i, op, pos, ref_end);
        }
    }
    None
}

#[cfg(test)]
#[allow(clippy::arithmetic_side_effects, reason = "test arithmetic on known small values")]
mod tests {
    use super::*;
    use proptest::prelude::*;

    fn op(op_type: CigarOpType, len: u32) -> CigarOp {
        CigarOp::new(op_type, len)
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
            assert_eq!(CigarOpType::from_bam(code), op_type, "from_bam failed for code {code}");

            let cigar_op = CigarOp::new(op_type, 42);
            let packed = cigar_op.to_bam_u32();
            let decoded = CigarOp::from_bam_u32(packed);
            assert_eq!(decoded, cigar_op, "round-trip failed for {op_type:?}");
        }
    }

    #[test]
    fn cigar_op_max_length() {
        // Upper 28 bits can hold up to 2^28 - 1 = 268_435_455
        let max_len = (1u32 << 28) - 1;
        let op = CigarOp::new(CigarOpType::Match, max_len);
        let packed = op.to_bam_u32();
        let decoded = CigarOp::from_bam_u32(packed);
        assert_eq!(decoded.len(), max_len);
        assert_eq!(decoded.op_type(), CigarOpType::Match);
    }

    #[test]
    fn cigar_op_unknown_code_round_trips() {
        // Op codes 9..=15 are reserved; they decode to Unknown and round-trip verbatim.
        let packed = (100u32 << 4) | 9;
        let op = CigarOp::from_bam_u32(packed);
        assert_eq!(op.op_type(), CigarOpType::Unknown(9));
        assert_eq!(op.len(), 100);
        assert_eq!(op.to_bam_u32(), packed);
        // Unknown ops consume neither ref nor query.
        assert!(!op.consumes_ref());
        assert!(!op.consumes_query());
    }

    // r[verify cigar.compact_op_position_invariant]
    #[test]
    fn compact_op_ref_start_fits_i32_for_bam_positions() {
        // BAM positions are i32 (max 2^31-1 = 2_147_483_647).
        // CompactOp stores ref_start as i32. Verify it works at the maximum BAM position.
        let cigar = [op(CigarOpType::Match, 100)];
        let rec_pos = Pos0::new((i32::MAX as u32) - 100).unwrap(); // near max BAM position
        let mapping = CigarMapping::new(rec_pos, &cigar).unwrap();
        // Should work fine — position fits in i32
        assert!(matches!(mapping, CigarMapping::Linear { .. }));
        assert_eq!(mapping.pos_info_at(rec_pos), Some(CigarPosInfo::Match { qpos: 0 }));
    }

    #[test]
    fn cigar_op_type_display_matches_sam_spec() {
        assert_eq!(CigarOpType::Match.to_string(), "M");
        assert_eq!(CigarOpType::Insertion.to_string(), "I");
        assert_eq!(CigarOpType::Deletion.to_string(), "D");
        assert_eq!(CigarOpType::RefSkip.to_string(), "N");
        assert_eq!(CigarOpType::SoftClip.to_string(), "S");
        assert_eq!(CigarOpType::HardClip.to_string(), "H");
        assert_eq!(CigarOpType::Padding.to_string(), "P");
        assert_eq!(CigarOpType::SeqMatch.to_string(), "=");
        assert_eq!(CigarOpType::SeqMismatch.to_string(), "X");
        assert_eq!(CigarOpType::Unknown(9).to_string(), "?");
    }

    #[test]
    fn cigar_op_display_concatenates_len_and_op() {
        assert_eq!(op(CigarOpType::Match, 10).to_string(), "10M");
        assert_eq!(op(CigarOpType::Insertion, 2).to_string(), "2I");
        assert_eq!(op(CigarOpType::SeqMatch, 4).to_string(), "4=");
    }

    #[test]
    fn cigar_str_renders_full_cigar() {
        let ops = [
            op(CigarOpType::Match, 10),
            op(CigarOpType::Insertion, 2),
            op(CigarOpType::Deletion, 3),
        ];
        assert_eq!(CigarStr(&ops).to_string(), "10M2I3D");

        let clipped = [
            op(CigarOpType::SoftClip, 5),
            op(CigarOpType::Match, 50),
            op(CigarOpType::SoftClip, 3),
        ];
        assert_eq!(CigarStr(&clipped).to_string(), "5S50M3S");

        assert_eq!(CigarStr(&[]).to_string(), "");
    }

    // r[verify cigar.qpos_bounds]
    #[test]
    fn linear_cigar_mapping_bounds_check() {
        // 5S + 100M: rec_pos=1000, query_offset=5, match_len=100
        let cigar = [op(CigarOpType::SoftClip, 5), op(CigarOpType::Match, 100)];

        let mapping = CigarMapping::new(Pos0::new(1000).unwrap(), &cigar).unwrap();
        assert!(matches!(mapping, CigarMapping::Linear { .. }));

        // Valid positions: 1000..1100
        assert_eq!(
            mapping.pos_info_at(Pos0::new(1000).unwrap()),
            Some(CigarPosInfo::Match { qpos: 5 })
        );
        assert_eq!(
            mapping.pos_info_at(Pos0::new(1050).unwrap()),
            Some(CigarPosInfo::Match { qpos: 55 })
        );
        assert_eq!(
            mapping.pos_info_at(Pos0::new(1099).unwrap()),
            Some(CigarPosInfo::Match { qpos: 104 })
        );

        // Out-of-range: before alignment start
        assert_eq!(
            mapping.pos_info_at(Pos0::new(999).unwrap()),
            None,
            "pos before rec_pos must return None"
        );

        // Out-of-range: at/past alignment end
        assert_eq!(
            mapping.pos_info_at(Pos0::new(1100).unwrap()),
            None,
            "pos at rec_pos + match_len must return None"
        );
        assert_eq!(
            mapping.pos_info_at(Pos0::new(1200).unwrap()),
            None,
            "pos past alignment end must return None"
        );

        // Far out-of-range (would wrap with unsigned subtraction)
        assert_eq!(
            mapping.pos_info_at(Pos0::new(0).unwrap()),
            None,
            "pos far before alignment must return None"
        );
    }

    // r[verify cigar.slice_from_bam_bytes]
    proptest::proptest! {
        /// Round-trip: encode CigarOp → BAM LE bytes → slice_from_bam_bytes → same ops.
        #[test]
        fn roundtrip_slice_from_bam_bytes(
            ops in proptest::collection::vec(
                (0u32..=8u32, 1u32..1000u32),
                1..20,
            ).prop_map(|pairs| {
                pairs.into_iter()
                    .map(|(code, len)| CigarOp::from_bam_u32((len << 4) | code))
                    .collect::<Vec<_>>()
            })
        ) {
            let bytes: Vec<u8> = ops.iter()
                .flat_map(|op| op.to_bam_u32().to_le_bytes())
                .collect();
            // 8 bytes padding ensures 8-byte alignment (≥ 4-byte for CigarOp)
            let mut buf = vec![0u8; 8];
            buf.extend_from_slice(&bytes);
            let aligned = &buf[8..];

            let result = CigarOp::slice_from_bam_bytes(aligned);
            prop_assert!(result.is_some(), "aligned buffer must succeed");
            prop_assert_eq!(result.unwrap(), ops.as_slice());
        }
    }

    #[test]
    fn slice_from_bam_bytes_rejects_odd_length() {
        assert!(CigarOp::slice_from_bam_bytes(&[0u8; 7]).is_none());
    }

    #[test]
    fn slice_from_bam_bytes_empty_ok() {
        assert_eq!(CigarOp::slice_from_bam_bytes(&[]), Some(&[][..]));
    }

    // r[verify cigar.slice_from_bam_bytes]
    #[test]
    fn slice_from_bam_bytes_rejects_unaligned_input() {
        // BAM CIGAR data can be unaligned within the record: the read name
        // sits between the 32-byte fixed header and the CIGAR data. If
        // l_read_name is not a multiple of 4, the CIGAR bytes start at a
        // non-4-byte-aligned offset. Samtools produces such records.
        let cigar_op = CigarOp::new(CigarOpType::Match, 100);
        let packed = cigar_op.to_bam_u32();
        let mut buf = vec![0u8];
        buf.extend_from_slice(&packed.to_le_bytes());
        buf.extend_from_slice(&packed.to_le_bytes());
        buf.extend_from_slice(&packed.to_le_bytes());

        let mut any_unaligned = false;
        for off in 0..=3 {
            let sub = &buf[off..off + 4];
            let result = CigarOp::slice_from_bam_bytes(sub);
            if sub.as_ptr() as usize % 4 != 0 {
                assert!(result.is_none(), "unaligned slice at offset {off} must return None");
                any_unaligned = true;
            }
        }
        assert!(any_unaligned, "at least one test offset must be unaligned");
    }
}
