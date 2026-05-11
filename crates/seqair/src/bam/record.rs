//! Shared BAM-record wire-format primitives — header parse + encode, raw-byte
//! readers, end-position computation, and the [`DecodeError`] enum.
//!
//! There is intentionally **no** `BamRecord` struct here. The production decode
//! path is [`RecordStore::push_raw`](super::record_store::RecordStore::push_raw),
//! which streams raw BAM bytes straight into the slab buffers without
//! allocating a per-record owned struct. The functions in this module
//! (`parse_header`, `encode_fixed_header`, `compute_end_pos_from_raw`,
//! `read2`, `read4`) are the wire-format pieces shared between the read and
//! write paths.
//!
//! If you need an owned, mutable record (for writing or in-memory editing),
//! use [`OwnedBamRecord`](super::owned_record::OwnedBamRecord). If you need
//! to decode a region's worth of records for pileup or read-level work, push
//! them into a [`RecordStore`](super::record_store::RecordStore) and access
//! fields via [`SlimRecord`](super::record_store::SlimRecord).

use seqair_types::{BamFlags, Pos0};

/// Compute `end_pos` from raw BAM record bytes (before full decode).
///
/// Walks the packed CIGAR bytes in place — used on the hot pre-push scan
/// path before we know whether to keep the record, so we can't afford to
/// allocate a `Vec<CigarOp>` here.
pub fn compute_end_pos_from_raw(raw: &[u8]) -> Option<Pos0> {
    use super::cigar::{CIGAR_D, CIGAR_EQ, CIGAR_M, CIGAR_N, CIGAR_X};
    use seqair_types::Offset;

    let h = parse_header(raw).ok()?;
    // All bounds ≤ qual_end ≤ raw.len() checked by parse_header
    debug_assert!(h.cigar_end <= raw.len(), "cigar overrun: {} > {}", h.cigar_end, raw.len());
    #[allow(clippy::indexing_slicing, reason = "cigar_end ≤ raw.len() checked by parse_header")]
    let cigar_bytes = &raw[h.var_start..h.cigar_end];

    let mut ref_len: i64 = 0;
    for chunk in cigar_bytes.chunks_exact(4) {
        let arr: [u8; 4] = chunk.try_into().expect("chunks_exact(4) yields 4 bytes");
        let op = u32::from_le_bytes(arr);
        let op_len = i64::from(op >> 4);
        let op_type = (op & 0xF) as u8;
        match op_type {
            CIGAR_M | CIGAR_D | CIGAR_N | CIGAR_EQ | CIGAR_X => {
                ref_len = ref_len.checked_add(op_len)?;
            }
            _ => {}
        }
    }
    if ref_len == 0 {
        Some(h.pos)
    } else {
        h.pos.checked_add_offset(Offset::new(ref_len.checked_sub(1)?))
    }
}

// Callers only invoke read2/read4 after validating raw.len() >= 32 or equivalent.
#[allow(
    clippy::indexing_slicing,
    reason = "offset + 1 < raw.len() ensured by caller's length check"
)]
pub(crate) fn read2(buf: &[u8], offset: usize) -> [u8; 2] {
    debug_assert!(
        offset.saturating_add(1) < buf.len(),
        "read2 out of bounds: offset={offset}, len={}",
        buf.len()
    );
    [buf[offset], buf[offset.wrapping_add(1)]]
}

#[allow(
    clippy::indexing_slicing,
    reason = "offset + 3 < raw.len() ensured by caller's length check"
)]
pub(crate) fn read4(buf: &[u8], offset: usize) -> [u8; 4] {
    debug_assert!(
        offset.saturating_add(3) < buf.len(),
        "read4 out of bounds: offset={offset}, len={}",
        buf.len()
    );
    [
        buf[offset],
        buf[offset.wrapping_add(1)],
        buf[offset.wrapping_add(2)],
        buf[offset.wrapping_add(3)],
    ]
}

/// Parsed BAM fixed header fields and computed variable-length offsets.
///
/// Used by [`RecordStore::push_raw`](super::record_store::RecordStore::push_raw)
/// and [`compute_end_pos_from_raw`] to avoid duplicating the 36-byte header
/// parsing and checked offset arithmetic.
#[derive(Debug)]
pub(crate) struct ParsedHeader {
    pub tid: i32,
    pub pos: Pos0,
    pub mapq: u8,
    pub flags: BamFlags,
    pub n_cigar_ops: u16,
    pub seq_len: u32,
    pub next_ref_id: i32,
    pub next_pos: i32,
    pub template_len: i32,
    /// Start of variable-length data (32 + `name_len`).
    pub var_start: usize,
    /// End of CIGAR bytes.
    pub cigar_end: usize,
    /// End of packed sequence bytes.
    pub seq_end: usize,
    /// End of quality scores.
    pub qual_end: usize,
}

/// Parse the fixed 36-byte BAM header and compute checked offsets for
/// variable-length fields.
///
/// # Errors
/// Returns `DecodeError::TooShort` if `raw` is shorter than 32 bytes or
/// shorter than `qual_end`, and `DecodeError::OffsetOverflow` if any offset
/// arithmetic overflows.
// r[impl bam.record.checked_offsets]
// r[impl bam.record.fields]
pub(crate) fn parse_header(raw: &[u8]) -> Result<ParsedHeader, DecodeError> {
    if raw.len() < 32 {
        return Err(DecodeError::TooShort { len: raw.len() });
    }

    debug_assert!(raw.len() >= 32, "raw record too short for fixed fields: {}", raw.len());
    let tid = i32::from_le_bytes(read4(raw, 0));
    let pos_i32 = i32::from_le_bytes(read4(raw, 4));
    let pos =
        Pos0::try_from(pos_i32).map_err(|_| DecodeError::InvalidPosition { value: pos_i32 })?;
    // raw.len() >= 32, so raw[8] and raw[9] are in bounds
    #[allow(clippy::indexing_slicing, reason = "raw.len() >= 32 checked above")]
    let name_len = raw[8] as usize;
    #[allow(clippy::indexing_slicing, reason = "raw.len() >= 32 checked above")]
    let mapq = raw[9];
    let n_cigar_ops = u16::from_le_bytes(read2(raw, 12));
    let flags = BamFlags::from(u16::from_le_bytes(read2(raw, 14)));
    let seq_len = u32::from_le_bytes(read4(raw, 16));
    let next_ref_id = i32::from_le_bytes(read4(raw, 20));
    let next_pos = i32::from_le_bytes(read4(raw, 24));
    let template_len = i32::from_le_bytes(read4(raw, 28));

    let cigar_bytes = usize::from(n_cigar_ops) * 4;
    let seq_bytes = (seq_len as usize).div_ceil(2);

    let var_start = 32usize.checked_add(name_len).ok_or(DecodeError::OffsetOverflow)?;
    let cigar_end = var_start.checked_add(cigar_bytes).ok_or(DecodeError::OffsetOverflow)?;
    let seq_end = cigar_end.checked_add(seq_bytes).ok_or(DecodeError::OffsetOverflow)?;
    let qual_end = seq_end.checked_add(seq_len as usize).ok_or(DecodeError::OffsetOverflow)?;

    if raw.len() < qual_end {
        return Err(DecodeError::TooShort { len: raw.len() });
    }

    Ok(ParsedHeader {
        tid,
        pos,
        mapq,
        flags,
        n_cigar_ops,
        seq_len,
        next_ref_id,
        next_pos,
        template_len,
        var_start,
        cigar_end,
        seq_end,
        qual_end,
    })
}

/// Field values for the 32-byte BAM fixed record header, in BAM wire types.
///
/// Callers validate field-size limits (qname ≤ 254, cigar ≤ 65535, seq ≤
/// `i32::MAX`) before constructing this; the encoder only packs bytes.
/// `pos` and `next_pos` are already resolved to wire form (`-1` for
/// unmapped/unavailable per `[SAM1] §1.4`).
pub(crate) struct FixedHeaderFields {
    pub ref_id: i32,
    pub pos: i32,
    pub bin: u16,
    pub mapq: u8,
    /// `qname.len() + 1` (NUL terminator); caller validated qname ≤ 254 bytes.
    pub l_read_name: u8,
    pub flags: u16,
    pub n_cigar_op: u16,
    pub l_seq: i32,
    pub next_ref_id: i32,
    pub next_pos: i32,
    pub template_len: i32,
}

/// Pack the 11 fixed-record fields into the canonical 32-byte BAM record
/// header and append to `buf`. Inverse of [`parse_header`].
pub(crate) fn encode_fixed_header(buf: &mut Vec<u8>, f: &FixedHeaderFields) {
    let bin_mq_nl = (u32::from(f.bin) << 16) | (u32::from(f.mapq) << 8) | u32::from(f.l_read_name);
    let flag_nc = (u32::from(f.flags) << 16) | u32::from(f.n_cigar_op);
    buf.extend_from_slice(&f.ref_id.to_le_bytes());
    buf.extend_from_slice(&f.pos.to_le_bytes());
    buf.extend_from_slice(&bin_mq_nl.to_le_bytes());
    buf.extend_from_slice(&flag_nc.to_le_bytes());
    buf.extend_from_slice(&f.l_seq.to_le_bytes());
    buf.extend_from_slice(&f.next_ref_id.to_le_bytes());
    buf.extend_from_slice(&f.next_pos.to_le_bytes());
    buf.extend_from_slice(&f.template_len.to_le_bytes());
}

#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum DecodeError {
    #[error("BAM record too short: {len} bytes")]
    TooShort { len: usize },

    #[error("arithmetic overflow computing BAM record field offsets")]
    OffsetOverflow,

    #[error("slab offset exceeds u32::MAX")]
    SlabOverflow,

    #[error("invalid BAM position value {value}: negative positions are reserved")]
    InvalidPosition { value: i32 },

    #[error("CIGAR query length {cigar_query_len} does not match seq_len {seq_len}")]
    CigarQueryLenMismatch { cigar_query_len: u32, seq_len: u32 },

    #[error("CIGAR op count {count} exceeds u16::MAX (BAM n_cigar_op limit)")]
    CigarOpCountOverflow { count: usize },

    #[error("qual length {qual_len} does not match seq length {seq_len}")]
    QualLenMismatch { qual_len: usize, seq_len: usize },
}

#[cfg(test)]
#[allow(clippy::arithmetic_side_effects, reason = "test arithmetic on known small values")]
mod tests {
    use super::*;

    #[test]
    fn test_compute_end_pos() {
        use super::super::cigar::{CigarOp, CigarOpType, compute_end_pos};
        let ops = [CigarOp::new(CigarOpType::Match, 50)];
        assert_eq!(compute_end_pos(Pos0::new(100).unwrap(), &ops), Some(Pos0::new(149).unwrap()));
    }

    // r[verify bam.record.checked_offsets]
    #[test]
    fn parse_header_rejects_overflow_in_offset_calc() {
        // Craft a 32-byte record with fields that would overflow if added unchecked:
        // name_len = 255, n_cigar_ops = 65535, seq_len = u32::MAX
        let mut raw = [0u8; 32];
        raw[0..4].copy_from_slice(&0i32.to_le_bytes()); // tid
        raw[4..8].copy_from_slice(&0i32.to_le_bytes()); // pos
        raw[8] = 255; // name_len (l_read_name)
        raw[9] = 0; // mapq
        raw[12..14].copy_from_slice(&u16::MAX.to_le_bytes()); // n_cigar_ops
        raw[14..16].copy_from_slice(&0u16.to_le_bytes()); // flags
        raw[16..20].copy_from_slice(&u32::MAX.to_le_bytes()); // seq_len

        let result = parse_header(&raw);
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(
            matches!(err, DecodeError::OffsetOverflow | DecodeError::TooShort { .. }),
            "expected OffsetOverflow or TooShort, got {err:?}"
        );
    }
}
