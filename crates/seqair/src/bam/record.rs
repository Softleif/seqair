//! Decode a BAM record from raw bytes into a [`BamRecord`] with owned variable-length fields.
//! Used transiently during region loading; the pileup engine works from [`crate::bam::RecordStore`] instead.

use super::{
    flags::{
        BamFlags, FLAG_FIRST_IN_TEMPLATE, FLAG_REVERSE, FLAG_SECOND_IN_TEMPLATE, FLAG_UNMAPPED,
    },
    seq,
};
use seqair_types::{Offset, Pos, Zero};

/// A decoded BAM record with owned variable-length data.
///
/// Wrapped in `Rc` for cheap sharing between the arena and pileup columns.
/// All accessor methods live directly on the struct — no separate `RecordRef`.
// r[impl bam.record.fields]
// r[impl bam.record.decode]
#[derive(Debug, Clone)]
pub struct BamRecord {
    pub pos: Pos<Zero>,
    pub end_pos: Pos<Zero>,
    pub tid: i32,
    pub seq_len: u32,
    pub flags: u16,
    pub n_cigar_ops: u16,
    pub mapq: u8,
    // r[impl perf.precompute_matches_indels]
    pub matching_bases: u32,
    pub indel_bases: u32,
    pub qname: Box<[u8]>,
    pub cigar: Box<[u8]>,
    pub seq: Box<[u8]>,
    pub qual: Box<[u8]>,
    pub aux: Box<[u8]>,
}

impl BamRecord {
    /// Decode from raw BAM bytes (after the 4-byte block_size prefix).
    pub fn decode(raw: &[u8]) -> Result<Self, DecodeError> {
        let h = parse_header(raw)?;
        let seq_len_usize = h.seq_len as usize;

        // All slice bounds (32..var_start, var_start..cigar_end, etc.) are ≤ qual_end ≤ raw.len()
        debug_assert!(h.qual_end <= raw.len(), "qual_end overrun: {} > {}", h.qual_end, raw.len());
        #[allow(
            clippy::indexing_slicing,
            reason = "all bounds ≤ qual_end ≤ raw.len() checked by parse_header"
        )]
        let qname_raw = &raw[32..h.var_start];
        let qname_actual_len = qname_raw.iter().position(|&b| b == 0).unwrap_or(qname_raw.len());

        #[allow(
            clippy::indexing_slicing,
            reason = "all bounds ≤ qual_end ≤ raw.len() checked by parse_header"
        )]
        let cigar_slice = &raw[h.var_start..h.cigar_end];
        let end_pos = compute_end_pos(h.pos, cigar_slice);
        let (matching_bases, indel_bases) = super::cigar::calc_matches_indels(cigar_slice);

        #[allow(
            clippy::indexing_slicing,
            reason = "all bounds ≤ qual_end ≤ raw.len() checked by parse_header"
        )]
        Ok(BamRecord {
            pos: h.pos,
            end_pos,
            tid: h.tid,
            seq_len: h.seq_len,
            flags: h.flags,
            n_cigar_ops: h.n_cigar_ops,
            mapq: h.mapq,
            matching_bases,
            indel_bases,
            qname: qname_raw[..qname_actual_len].into(),
            cigar: cigar_slice.into(),
            // r[impl bam.record.seq_4bit]
            // r[impl bam.record.seq_at_simd+2]
            seq: seq::decode_seq(&raw[h.cigar_end..h.seq_end], seq_len_usize).into_boxed_slice(),
            qual: raw[h.seq_end..h.qual_end].into(),
            aux: raw[h.qual_end..].into(),
        })
    }

    // --- accessors ---

    pub fn bam_flags(&self) -> BamFlags {
        BamFlags::new(self.flags)
    }

    // r[impl bam.record.flag_reverse]
    pub fn is_reverse(&self) -> bool {
        self.flags & FLAG_REVERSE != 0
    }

    // r[impl bam.record.flag_first]
    pub fn is_first_in_template(&self) -> bool {
        self.flags & FLAG_FIRST_IN_TEMPLATE != 0
    }

    // r[impl bam.record.flag_second]
    pub fn is_second_in_template(&self) -> bool {
        self.flags & FLAG_SECOND_IN_TEMPLATE != 0
    }

    // r[impl bam.record.flag_unmapped]
    pub fn is_unmapped(&self) -> bool {
        self.flags & FLAG_UNMAPPED != 0
    }

    // r[impl bam.record.seq_at]
    pub fn seq_at(&self, pos: usize) -> u8 {
        self.seq.get(pos).copied().unwrap_or(b'N')
    }

    // r[impl bam.record.aux_parse]
    pub fn aux(&self, tag: &[u8; 2]) -> Option<super::aux::AuxValue<'_>> {
        super::aux::find_tag(&self.aux, *tag)
    }
}

/// Compute end_pos from raw BAM record bytes (before full decode).
pub fn compute_end_pos_from_raw(raw: &[u8]) -> Option<Pos<Zero>> {
    let h = parse_header(raw).ok()?;
    // All bounds ≤ qual_end ≤ raw.len() checked by parse_header
    debug_assert!(h.cigar_end <= raw.len(), "cigar overrun: {} > {}", h.cigar_end, raw.len());
    #[allow(clippy::indexing_slicing, reason = "cigar_end ≤ raw.len() checked by parse_header")]
    Some(compute_end_pos(h.pos, &raw[h.var_start..h.cigar_end]))
}

// r[impl bam.record.end_pos]
// r[impl bam.record.zero_refspan]
pub(crate) fn compute_end_pos(pos: Pos<Zero>, cigar_bytes: &[u8]) -> Pos<Zero> {
    use super::cigar::{CIGAR_D, CIGAR_EQ, CIGAR_M, CIGAR_N, CIGAR_X};

    let mut ref_len: i64 = 0;
    let n_ops = cigar_bytes.len() / 4;
    for i in 0..n_ops {
        let op = u32::from_le_bytes(read4(cigar_bytes, i * 4));
        let op_len = i64::from(op >> 4);
        let op_type = (op & 0xF) as u8;
        match op_type {
            CIGAR_M | CIGAR_D | CIGAR_N | CIGAR_EQ | CIGAR_X => ref_len += op_len,
            _ => {}
        }
    }
    if ref_len == 0 { pos } else { pos + Offset(ref_len - 1) }
}

// Callers only invoke read2/read4 after validating raw.len() >= 32 or equivalent.
#[allow(
    clippy::indexing_slicing,
    reason = "offset + 1 < raw.len() ensured by caller's length check"
)]
pub(crate) fn read2(buf: &[u8], offset: usize) -> [u8; 2] {
    debug_assert!(
        offset + 1 < buf.len(),
        "read2 out of bounds: offset={offset}, len={}",
        buf.len()
    );
    [buf[offset], buf[offset + 1]]
}

#[allow(
    clippy::indexing_slicing,
    reason = "offset + 3 < raw.len() ensured by caller's length check"
)]
pub(crate) fn read4(buf: &[u8], offset: usize) -> [u8; 4] {
    debug_assert!(
        offset + 3 < buf.len(),
        "read4 out of bounds: offset={offset}, len={}",
        buf.len()
    );
    [buf[offset], buf[offset + 1], buf[offset + 2], buf[offset + 3]]
}

/// Parsed BAM fixed header fields and computed variable-length offsets.
///
/// Shared between `BamRecord::decode` and `RecordStore::push_raw` to avoid
/// duplicating the 36-byte header parsing and checked offset arithmetic.
pub(crate) struct ParsedHeader {
    pub tid: i32,
    pub pos: Pos<Zero>,
    pub mapq: u8,
    pub flags: u16,
    pub n_cigar_ops: u16,
    pub seq_len: u32,
    /// Start of variable-length data (32 + name_len).
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
pub(crate) fn parse_header(raw: &[u8]) -> Result<ParsedHeader, DecodeError> {
    if raw.len() < 32 {
        return Err(DecodeError::TooShort { len: raw.len() });
    }

    debug_assert!(raw.len() >= 32, "raw record too short for fixed fields: {}", raw.len());
    let tid = i32::from_le_bytes(read4(raw, 0));
    let pos_i32 = i32::from_le_bytes(read4(raw, 4));
    let pos = Pos::<Zero>::try_from_i64(i64::from(pos_i32)).unwrap_or(Pos::<Zero>::new(0));
    // raw.len() >= 32, so raw[8] and raw[9] are in bounds
    #[allow(clippy::indexing_slicing, reason = "raw.len() >= 32 checked above")]
    let name_len = raw[8] as usize;
    #[allow(clippy::indexing_slicing, reason = "raw.len() >= 32 checked above")]
    let mapq = raw[9];
    let n_cigar_ops = u16::from_le_bytes(read2(raw, 12));
    let flags = u16::from_le_bytes(read2(raw, 14));
    let seq_len = u32::from_le_bytes(read4(raw, 16));

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
        var_start,
        cigar_end,
        seq_end,
        qual_end,
    })
}

#[derive(Debug, thiserror::Error)]
pub enum DecodeError {
    #[error("BAM record too short: {len} bytes")]
    TooShort { len: usize },

    #[error("arithmetic overflow computing BAM record field offsets")]
    OffsetOverflow,

    #[error("slab offset exceeds u32::MAX")]
    SlabOverflow,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compute_end_pos() {
        let op = 50u32 << 4;
        assert_eq!(
            compute_end_pos(Pos::<Zero>::new(100), &op.to_le_bytes()),
            Pos::<Zero>::new(149)
        );
    }

    // r[verify bam.record.checked_offsets]
    #[test]
    fn decode_rejects_overflow_in_offset_calc() {
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

        let result = BamRecord::decode(&raw);
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(
            matches!(err, DecodeError::OffsetOverflow | DecodeError::TooShort { .. }),
            "expected OffsetOverflow or TooShort, got {err:?}"
        );
    }

    #[test]
    fn test_decode_record() {
        let mut raw = [0u8; 64];
        raw[0..4].copy_from_slice(&0i32.to_le_bytes());
        raw[4..8].copy_from_slice(&100i32.to_le_bytes());
        raw[8] = 5;
        raw[9] = 60;
        raw[10..12].copy_from_slice(&0u16.to_le_bytes());
        raw[12..14].copy_from_slice(&1u16.to_le_bytes());
        raw[14..16].copy_from_slice(&99u16.to_le_bytes());
        raw[16..20].copy_from_slice(&4u32.to_le_bytes());
        raw[20..32].fill(0);
        raw[32..37].copy_from_slice(b"read\0");
        raw[37..41].copy_from_slice(&(4u32 << 4).to_le_bytes());
        raw[41] = 0x12;
        raw[42] = 0x48;
        raw[43..47].copy_from_slice(&[30, 30, 30, 30]);

        let rec = BamRecord::decode(&raw[..47]).unwrap();
        assert_eq!(rec.pos, Pos::<Zero>::new(100));
        assert_eq!(rec.end_pos, Pos::<Zero>::new(103));
        assert_eq!(rec.flags, 99);
        assert_eq!(rec.mapq, 60);
        assert_eq!(&*rec.qname, b"read");
        assert_eq!(rec.seq_at(0), b'A');
        assert_eq!(rec.seq_at(3), b'T');
        assert_eq!(rec.matching_bases, 4);
    }
}
