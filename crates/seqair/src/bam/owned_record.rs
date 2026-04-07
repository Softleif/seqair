//! Owned, mutable BAM record for writing and in-memory modification.
//!
//! [`OwnedBamRecord`] is the write-path counterpart to the read-path [`BamRecord`](super::record::BamRecord).
//! It stores all fields in owned, mutable collections and supports serialization to BAM binary format
//! via [`to_bam_bytes`](OwnedBamRecord::to_bam_bytes).

use super::aux_data::{AuxData, AuxDataError};
use super::cigar::CigarOp;
use super::flags::BamFlags;
use super::seq;
use crate::vcf::index_builder::reg2bin;
use seqair_types::Base;
use thiserror::Error;

/// Errors from owned record construction or serialization.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum OwnedRecordError {
    /// Query name exceeds 254 bytes (BAM `l_read_name` is u8, includes NUL).
    #[error("qname length {len} exceeds maximum of 254 bytes")]
    QnameTooLong { len: usize },

    /// CIGAR operation count exceeds 65535 (BAM `n_cigar_op` is u16).
    #[error("CIGAR op count {count} exceeds maximum of 65535")]
    CigarCountOverflow { count: usize },

    /// Sequence length exceeds i32::MAX (BAM `l_seq` is i32).
    #[error("sequence length {len} exceeds i32::MAX")]
    SeqLengthOverflow { len: usize },

    /// Position does not fit in BAM i32 field.
    #[error("pos {pos} does not fit in BAM i32 range [-1, {max}]")]
    PosOverflow { pos: i64, max: i32 },

    /// Next position does not fit in BAM i32 field.
    #[error("next_pos {next_pos} does not fit in BAM i32 range [-1, {max}]")]
    NextPosOverflow { next_pos: i64, max: i32 },

    /// Sequence and quality score lengths do not match.
    #[error("seq length {seq_len} != qual length {qual_len}")]
    SeqQualLengthMismatch { seq_len: usize, qual_len: usize },

    /// CIGAR query-consuming length does not match sequence length.
    #[error("CIGAR query length {cigar_query_len} != seq length {seq_len}")]
    CigarSeqLengthMismatch { cigar_query_len: u64, seq_len: usize },

    /// Serialized record exceeds size limit.
    #[error("record size {size} exceeds limit")]
    RecordTooLarge { size: usize },

    /// Auxiliary tag error.
    #[error("aux data error: {source}")]
    AuxData {
        #[from]
        source: AuxDataError,
    },
}

// r[impl bam.owned_record.fields]
/// An owned, mutable BAM record for writing and in-memory modification.
#[derive(Debug, Clone)]
pub struct OwnedBamRecord {
    pub ref_id: i32,
    pub pos: i64,
    pub mapq: u8,
    pub flags: u16,
    pub next_ref_id: i32,
    pub next_pos: i64,
    pub template_len: i32,
    pub qname: Vec<u8>,
    pub cigar: Vec<CigarOp>,
    pub seq: Vec<Base>,
    pub qual: Vec<u8>,
    pub aux: AuxData,
}

/// Builder for constructing [`OwnedBamRecord`] values.
pub struct OwnedBamRecordBuilder {
    ref_id: i32,
    pos: i64,
    mapq: u8,
    flags: u16,
    next_ref_id: i32,
    next_pos: i64,
    template_len: i32,
    qname: Vec<u8>,
    cigar: Vec<CigarOp>,
    seq: Vec<Base>,
    qual: Vec<u8>,
    aux: AuxData,
}

// r[impl bam.owned_record.builder]
impl OwnedBamRecordBuilder {
    pub fn flags(mut self, flags: u16) -> Self {
        self.flags = flags;
        self
    }
    pub fn mapq(mut self, mapq: u8) -> Self {
        self.mapq = mapq;
        self
    }
    pub fn cigar(mut self, cigar: Vec<CigarOp>) -> Self {
        self.cigar = cigar;
        self
    }
    pub fn seq(mut self, seq: Vec<Base>) -> Self {
        self.seq = seq;
        self
    }
    pub fn qual(mut self, qual: Vec<u8>) -> Self {
        self.qual = qual;
        self
    }
    pub fn next_ref_id(mut self, v: i32) -> Self {
        self.next_ref_id = v;
        self
    }
    pub fn next_pos(mut self, v: i64) -> Self {
        self.next_pos = v;
        self
    }
    pub fn template_len(mut self, v: i32) -> Self {
        self.template_len = v;
        self
    }
    pub fn aux(mut self, aux: AuxData) -> Self {
        self.aux = aux;
        self
    }

    // r[impl bam.owned_record.qname_limit]
    // r[impl bam.owned_record.cigar_count_limit]
    // r[impl bam.owned_record.seq_length_limit]
    pub fn build(self) -> Result<OwnedBamRecord, OwnedRecordError> {
        if self.qname.len() > 254 {
            return Err(OwnedRecordError::QnameTooLong { len: self.qname.len() });
        }
        if self.cigar.len() > 65535 {
            return Err(OwnedRecordError::CigarCountOverflow { count: self.cigar.len() });
        }
        if self.seq.len() > i32::MAX as usize {
            return Err(OwnedRecordError::SeqLengthOverflow { len: self.seq.len() });
        }
        if !self.qual.is_empty() && self.qual.len() != self.seq.len() {
            return Err(OwnedRecordError::SeqQualLengthMismatch {
                seq_len: self.seq.len(),
                qual_len: self.qual.len(),
            });
        }

        Ok(OwnedBamRecord {
            ref_id: self.ref_id,
            pos: self.pos,
            mapq: self.mapq,
            flags: self.flags,
            next_ref_id: self.next_ref_id,
            next_pos: self.next_pos,
            template_len: self.template_len,
            qname: self.qname,
            cigar: self.cigar,
            seq: self.seq,
            qual: self.qual,
            aux: self.aux,
        })
    }
}

impl OwnedBamRecord {
    /// Decode a complete owned record from raw BAM bytes (after the 4-byte block_size prefix).
    ///
    /// Unlike [`BamRecord::decode`](super::record::BamRecord::decode), this preserves ALL fields
    /// including mate info (next_ref_id, next_pos, template_len) that the read-path type drops.
    /// This is essential for BAM rewriting pipelines that need to modify records and write them back.
    pub fn from_raw_bam(raw: &[u8]) -> Result<Self, OwnedRecordError> {
        use super::record::parse_header;

        let h = parse_header(raw).map_err(|_| OwnedRecordError::QnameTooLong { len: 0 })?;

        // Extract mate fields from the 32-byte fixed header (offsets 20-31)
        // parse_header already validated raw.len() >= 32
        #[allow(clippy::indexing_slicing, reason = "raw.len() >= 32 validated by parse_header")]
        let next_ref_id = i32::from_le_bytes([raw[20], raw[21], raw[22], raw[23]]);
        #[allow(clippy::indexing_slicing, reason = "raw.len() >= 32 validated by parse_header")]
        let next_pos_i32 = i32::from_le_bytes([raw[24], raw[25], raw[26], raw[27]]);
        #[allow(clippy::indexing_slicing, reason = "raw.len() >= 32 validated by parse_header")]
        let template_len = i32::from_le_bytes([raw[28], raw[29], raw[30], raw[31]]);

        // Extract qname (strip trailing NUL)
        #[allow(clippy::indexing_slicing, reason = "bounds checked by parse_header")]
        let qname_raw = &raw[32..h.var_start];
        let qname_len = qname_raw.iter().position(|&b| b == 0).unwrap_or(qname_raw.len());
        let qname = qname_raw.get(..qname_len).unwrap_or(qname_raw).to_vec();

        // Convert packed CIGAR u32s to Vec<CigarOp>
        #[allow(clippy::indexing_slicing, reason = "bounds checked by parse_header")]
        let cigar_bytes = &raw[h.var_start..h.cigar_end];
        let mut cigar = Vec::with_capacity(h.n_cigar_ops as usize);
        for i in 0..h.n_cigar_ops as usize {
            let offset =
                i.checked_mul(4).ok_or(OwnedRecordError::CigarCountOverflow { count: i })?;
            let packed = u32::from_le_bytes([
                *cigar_bytes
                    .get(offset)
                    .ok_or(OwnedRecordError::CigarCountOverflow { count: i })?,
                *cigar_bytes
                    .get(offset.saturating_add(1))
                    .ok_or(OwnedRecordError::CigarCountOverflow { count: i })?,
                *cigar_bytes
                    .get(offset.saturating_add(2))
                    .ok_or(OwnedRecordError::CigarCountOverflow { count: i })?,
                *cigar_bytes
                    .get(offset.saturating_add(3))
                    .ok_or(OwnedRecordError::CigarCountOverflow { count: i })?,
            ]);
            let op = CigarOp::from_bam_u32(packed)
                .ok_or(OwnedRecordError::CigarCountOverflow { count: i })?;
            cigar.push(op);
        }

        // Decode 4-bit packed sequence to Base values
        #[allow(clippy::indexing_slicing, reason = "bounds checked by parse_header")]
        let seq_packed = &raw[h.cigar_end..h.seq_end];
        let seq_ascii = super::seq::decode_seq(seq_packed, h.seq_len as usize);
        let seq: Vec<Base> = seq_ascii.iter().map(|&b| Base::from(b)).collect();

        // Quality scores
        #[allow(clippy::indexing_slicing, reason = "bounds checked by parse_header")]
        let qual = raw[h.seq_end..h.qual_end].to_vec();

        // Auxiliary data (everything after qual)
        #[allow(clippy::indexing_slicing, reason = "qual_end <= raw.len()")]
        let aux = AuxData::from_bytes(raw[h.qual_end..].to_vec());

        Ok(OwnedBamRecord {
            ref_id: h.tid,
            pos: i64::from(h.pos.get() as i32),
            mapq: h.mapq,
            flags: h.flags,
            next_ref_id,
            next_pos: i64::from(next_pos_i32),
            template_len,
            qname,
            cigar,
            seq,
            qual,
            aux,
        })
    }

    /// Start building a record with the required fields.
    pub fn builder(ref_id: i32, pos: i64, qname: Vec<u8>) -> OwnedBamRecordBuilder {
        OwnedBamRecordBuilder {
            ref_id,
            pos,
            mapq: 0,
            flags: 0,
            next_ref_id: -1,
            next_pos: -1,
            template_len: 0,
            qname,
            cigar: Vec::new(),
            seq: Vec::new(),
            qual: Vec::new(),
            aux: AuxData::new(),
        }
    }

    /// Flag predicates via the [`BamFlags`] newtype.
    pub fn bam_flags(&self) -> BamFlags {
        BamFlags::new(self.flags)
    }

    // r[impl bam.owned_record.flag_methods]
    pub fn set_flags(&mut self, flags: u16) {
        self.flags = flags;
    }

    pub fn is_unmapped(&self) -> bool {
        self.bam_flags().is_unmapped()
    }

    pub fn is_reverse(&self) -> bool {
        self.bam_flags().is_reverse()
    }

    // r[impl bam.owned_record.end_pos]
    /// Compute the 0-based exclusive end position from pos + reference-consuming CIGAR ops.
    pub fn end_pos(&self) -> i64 {
        let ref_len: u64 =
            self.cigar.iter().filter(|op| op.op.consumes_ref()).map(|op| u64::from(op.len)).sum();
        self.pos.saturating_add(ref_len as i64)
    }

    // r[impl bam.owned_record.bin]
    /// Compute the BAI bin value from pos and end_pos (BAI scheme: min_shift=14, depth=5).
    pub fn bin(&self) -> u16 {
        let beg = if self.pos < 0 { 0u64 } else { self.pos as u64 };
        let end = if self.end_pos() <= self.pos {
            beg.saturating_add(1) // zero-span reads get beg..beg+1
        } else {
            self.end_pos() as u64
        };
        reg2bin(beg, end, 14, 5) as u16
    }

    // r[impl bam.owned_record.aligned_pairs]
    /// Iterator over (query_pos, ref_pos) alignment pairs.
    ///
    /// Yields one tuple per consumed base:
    /// - M/=/X: `(Some(qpos), Some(rpos))`
    /// - I: `(Some(qpos), None)`
    /// - D/N: `(None, Some(rpos))`
    /// - S/H/P: skipped entirely
    pub fn aligned_pairs(&self) -> AlignedPairsIter<'_> {
        AlignedPairsIter {
            cigar: &self.cigar,
            op_idx: 0,
            op_offset: 0,
            query_pos: 0,
            ref_pos: self.pos,
        }
    }

    /// Query-consuming length of the CIGAR.
    pub fn cigar_query_len(&self) -> u64 {
        self.cigar.iter().filter(|op| op.op.consumes_query()).map(|op| u64::from(op.len)).sum()
    }

    // r[impl bam.owned_record.set_alignment]
    /// Replace alignment (pos + cigar). Validates CIGAR query length == seq length for mapped reads.
    pub fn set_alignment(&mut self, pos: i64, cigar: Vec<CigarOp>) -> Result<(), OwnedRecordError> {
        if cigar.len() > 65535 {
            return Err(OwnedRecordError::CigarCountOverflow { count: cigar.len() });
        }
        // For mapped reads, validate CIGAR query length matches seq length
        if !self.is_unmapped() && !self.seq.is_empty() {
            let query_len: u64 =
                cigar.iter().filter(|op| op.op.consumes_query()).map(|op| u64::from(op.len)).sum();
            if query_len != self.seq.len() as u64 {
                return Err(OwnedRecordError::CigarSeqLengthMismatch {
                    cigar_query_len: query_len,
                    seq_len: self.seq.len(),
                });
            }
        }
        self.pos = pos;
        self.cigar = cigar;
        Ok(())
    }

    // r[impl bam.owned_record.set_seq]
    /// Replace the sequence. For mapped reads, validates length matches CIGAR query length.
    pub fn set_seq(&mut self, seq: Vec<Base>) -> Result<(), OwnedRecordError> {
        if !self.is_unmapped() && !self.cigar.is_empty() {
            let query_len = self.cigar_query_len();
            if query_len != seq.len() as u64 {
                return Err(OwnedRecordError::CigarSeqLengthMismatch {
                    cigar_query_len: query_len,
                    seq_len: seq.len(),
                });
            }
        }
        self.seq = seq;
        Ok(())
    }

    /// Mutable access to the sequence for in-place base modification.
    pub fn seq_mut(&mut self) -> &mut [Base] {
        &mut self.seq
    }

    // r[impl bam.owned_record.set_qual]
    /// Replace quality scores. Length must equal seq length.
    pub fn set_qual(&mut self, qual: Vec<u8>) -> Result<(), OwnedRecordError> {
        if !qual.is_empty() && qual.len() != self.seq.len() {
            return Err(OwnedRecordError::SeqQualLengthMismatch {
                seq_len: self.seq.len(),
                qual_len: qual.len(),
            });
        }
        self.qual = qual;
        Ok(())
    }

    // r[impl bam.owned_record.to_bam_bytes]
    // r[impl bam.owned_record.to_bam_field_overflow]
    // r[impl bam.owned_record.to_bam_bin_field]
    /// Serialize to BAM binary format by appending into `buf`.
    ///
    /// Does NOT include the 4-byte `block_size` prefix — that is the writer's responsibility.
    pub fn to_bam_bytes(&self, buf: &mut Vec<u8>) -> Result<(), OwnedRecordError> {
        // Validate field limits
        if self.qname.len() > 254 {
            return Err(OwnedRecordError::QnameTooLong { len: self.qname.len() });
        }
        if self.cigar.len() > 65535 {
            return Err(OwnedRecordError::CigarCountOverflow { count: self.cigar.len() });
        }
        if self.seq.len() > i32::MAX as usize {
            return Err(OwnedRecordError::SeqLengthOverflow { len: self.seq.len() });
        }
        if self.pos < -1 || self.pos > i64::from(i32::MAX) {
            return Err(OwnedRecordError::PosOverflow { pos: self.pos, max: i32::MAX });
        }
        if self.next_pos < -1 || self.next_pos > i64::from(i32::MAX) {
            return Err(OwnedRecordError::NextPosOverflow {
                next_pos: self.next_pos,
                max: i32::MAX,
            });
        }

        // Compute bin at serialization time
        let bin = self.bin();
        let l_read_name = (self.qname.len() as u32).saturating_add(1); // +1 for NUL
        let n_cigar_op = self.cigar.len() as u32;
        let l_seq = self.seq.len() as i32;

        // Pack bin_mq_nl: bin(16) | mapq(8) | l_read_name(8)
        let bin_mq_nl = (u32::from(bin) << 16) | (u32::from(self.mapq) << 8) | (l_read_name & 0xFF);

        // Pack flag_nc: flag(16) | n_cigar_op(16)
        let flag_nc = (u32::from(self.flags) << 16) | (n_cigar_op & 0xFFFF);

        // 32-byte fixed header
        buf.extend_from_slice(&self.ref_id.to_le_bytes());
        buf.extend_from_slice(&(self.pos as i32).to_le_bytes());
        buf.extend_from_slice(&bin_mq_nl.to_le_bytes());
        buf.extend_from_slice(&flag_nc.to_le_bytes());
        buf.extend_from_slice(&l_seq.to_le_bytes());
        buf.extend_from_slice(&self.next_ref_id.to_le_bytes());
        buf.extend_from_slice(&(self.next_pos as i32).to_le_bytes());
        buf.extend_from_slice(&self.template_len.to_le_bytes());

        // qname + NUL
        buf.extend_from_slice(&self.qname);
        buf.push(0);

        // CIGAR: each op as u32 LE
        for op in &self.cigar {
            buf.extend_from_slice(&op.to_bam_u32().to_le_bytes());
        }

        // Sequence: encode Base values to 4-bit packed
        if !self.seq.is_empty() {
            // Base is #[repr(u8)] with ASCII discriminants (A=0x41, C=0x43, G=0x47, T=0x54, N=0x4E).
            // The compile-time assert in seqair-types guarantees size_of::<Base>() == 1.
            // We collect to a byte vec rather than using unsafe transmute.
            let seq_bytes: Vec<u8> = self.seq.iter().map(|b| *b as u8).collect();
            let encoded = seq::encode_seq(&seq_bytes);
            buf.extend_from_slice(&encoded);
        }

        // Quality scores
        if self.qual.is_empty() {
            // Missing qual: fill with 0xFF per SAM spec
            buf.resize(buf.len().saturating_add(self.seq.len()), 0xFF);
        } else {
            buf.extend_from_slice(&self.qual);
        }

        // Auxiliary data
        buf.extend_from_slice(self.aux.as_bytes());

        Ok(())
    }
}

/// Iterator yielding `(Option<usize>, Option<i64>)` alignment pairs from CIGAR operations.
pub struct AlignedPairsIter<'a> {
    cigar: &'a [CigarOp],
    op_idx: usize,
    op_offset: u32,
    query_pos: usize,
    ref_pos: i64,
}

impl<'a> Iterator for AlignedPairsIter<'a> {
    type Item = (Option<usize>, Option<i64>);

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let op = self.cigar.get(self.op_idx)?;
            if self.op_offset >= op.len {
                self.op_idx = self.op_idx.checked_add(1)?;
                self.op_offset = 0;
                continue;
            }

            let consumes_q = op.op.consumes_query();
            let consumes_r = op.op.consumes_ref();

            // S/H/P: skip entirely
            if !consumes_q && !consumes_r {
                self.op_idx = self.op_idx.checked_add(1)?;
                self.op_offset = 0;
                continue;
            }

            // Soft clip: advance query but don't yield
            if matches!(op.op, super::cigar::CigarOpType::SoftClip) {
                self.query_pos = self.query_pos.saturating_add(1);
                self.op_offset = self.op_offset.saturating_add(1);
                continue;
            }

            let qpos = if consumes_q {
                let q = self.query_pos;
                self.query_pos = self.query_pos.saturating_add(1);
                Some(q)
            } else {
                None
            };

            let rpos = if consumes_r {
                let r = self.ref_pos;
                self.ref_pos = self.ref_pos.saturating_add(1);
                Some(r)
            } else {
                None
            };

            self.op_offset = self.op_offset.saturating_add(1);
            return Some((qpos, rpos));
        }
    }
}

#[cfg(test)]
#[allow(clippy::arithmetic_side_effects, reason = "test arithmetic on known small values")]
mod tests {
    use super::super::cigar::CigarOpType;
    use super::*;

    fn simple_record() -> OwnedBamRecord {
        OwnedBamRecord::builder(0, 100, b"read1".to_vec())
            .flags(0) // mapped
            .mapq(30)
            .cigar(vec![CigarOp::new(CigarOpType::Match, 5)])
            .seq(vec![Base::A, Base::C, Base::G, Base::T, Base::A])
            .qual(vec![30, 31, 32, 33, 34])
            .build()
            .unwrap()
    }

    // r[verify bam.owned_record.to_bam_bytes]
    // r[verify bam.owned_record.test_roundtrip_bytes]
    #[test]
    fn serialize_and_decode_roundtrip() {
        let rec = simple_record();
        let mut buf = Vec::new();
        rec.to_bam_bytes(&mut buf).unwrap();

        // Prepend block_size and decode with the existing BamRecord::decode
        let block_size = buf.len() as i32;
        let mut full = block_size.to_le_bytes().to_vec();
        full.extend_from_slice(&buf);

        let decoded = super::super::record::BamRecord::decode(&buf).unwrap();

        assert_eq!(decoded.tid, 0);
        assert_eq!(decoded.pos.get(), 100);
        assert_eq!(decoded.mapq, 30);
        assert_eq!(decoded.flags, 0);
        assert_eq!(decoded.seq_len, 5);
        assert_eq!(decoded.n_cigar_ops, 1);
        assert_eq!(&*decoded.qname, b"read1");
    }

    #[test]
    fn serialize_with_aux_tags() {
        let mut aux = AuxData::new();
        aux.set_int(*b"NM", 3).unwrap();
        aux.set_string(*b"RG", b"grp1");

        let rec = OwnedBamRecord::builder(0, 100, b"r1".to_vec())
            .cigar(vec![CigarOp::new(CigarOpType::Match, 3)])
            .seq(vec![Base::A, Base::C, Base::G])
            .qual(vec![30, 31, 32])
            .aux(aux)
            .build()
            .unwrap();

        let mut buf = Vec::new();
        rec.to_bam_bytes(&mut buf).unwrap();

        let decoded = super::super::record::BamRecord::decode(&buf).unwrap();
        let nm = super::super::aux::find_tag(&decoded.aux, *b"NM");
        assert_eq!(nm, Some(super::super::aux::AuxValue::U8(3)));

        let rg = super::super::aux::find_tag(&decoded.aux, *b"RG");
        assert_eq!(rg, Some(super::super::aux::AuxValue::String(b"grp1")));
    }

    // r[verify bam.owned_record.to_bam_bin_field]
    #[test]
    fn bin_computed_at_serialize_time() {
        let rec = simple_record();
        let bin = rec.bin();
        // 100..105 should be in a leaf bin
        assert!(bin > 0);
    }

    // r[verify bam.owned_record.end_pos]
    #[test]
    fn end_pos_simple() {
        let rec = simple_record();
        assert_eq!(rec.end_pos(), 105); // pos=100 + 5M
    }

    #[test]
    fn end_pos_with_insertion() {
        let rec = OwnedBamRecord::builder(0, 100, b"r".to_vec())
            .cigar(vec![
                CigarOp::new(CigarOpType::Match, 3),
                CigarOp::new(CigarOpType::Insertion, 2),
                CigarOp::new(CigarOpType::Match, 3),
            ])
            .seq(vec![Base::A; 8]) // 3M + 2I + 3M = 8 query bases
            .qual(vec![30; 8])
            .build()
            .unwrap();
        // Only M ops consume ref: 3 + 3 = 6
        assert_eq!(rec.end_pos(), 106);
    }

    // r[verify bam.owned_record.to_bam_field_overflow]
    #[test]
    fn pos_overflow_rejected() {
        let rec =
            OwnedBamRecord::builder(0, i64::from(i32::MAX) + 1, b"r".to_vec()).build().unwrap();
        let mut buf = Vec::new();
        assert!(rec.to_bam_bytes(&mut buf).is_err());
    }

    #[test]
    fn next_pos_overflow_rejected() {
        let rec = OwnedBamRecord::builder(0, 0, b"r".to_vec())
            .next_pos(i64::from(i32::MAX) + 1)
            .build()
            .unwrap();
        let mut buf = Vec::new();
        assert!(rec.to_bam_bytes(&mut buf).is_err());
    }

    // r[verify bam.owned_record.qname_limit]
    #[test]
    fn qname_too_long_rejected() {
        let long_name = vec![b'A'; 255];
        let result = OwnedBamRecord::builder(0, 0, long_name).build();
        assert!(result.is_err());
    }

    #[test]
    fn qname_at_limit_accepted() {
        let name = vec![b'A'; 254];
        let result = OwnedBamRecord::builder(0, 0, name).build();
        assert!(result.is_ok());
    }

    // r[verify bam.owned_record.set_alignment]
    // r[verify bam.owned_record.test_modification]
    #[test]
    fn set_alignment_validates_cigar_seq_mismatch() {
        let mut rec = simple_record(); // seq len = 5
        let result = rec.set_alignment(
            200,
            vec![CigarOp::new(CigarOpType::Match, 10)], // query len = 10 != 5
        );
        assert!(result.is_err());
    }

    #[test]
    fn set_alignment_succeeds_when_matching() {
        let mut rec = simple_record(); // seq len = 5
        let result = rec.set_alignment(200, vec![CigarOp::new(CigarOpType::Match, 5)]);
        assert!(result.is_ok());
        assert_eq!(rec.pos, 200);
        assert_eq!(rec.end_pos(), 205);
    }

    #[test]
    fn empty_qual_fills_with_0xff() {
        let rec = OwnedBamRecord::builder(0, 100, b"r".to_vec())
            .cigar(vec![CigarOp::new(CigarOpType::Match, 3)])
            .seq(vec![Base::A, Base::C, Base::G])
            // no qual set
            .build()
            .unwrap();

        let mut buf = Vec::new();
        rec.to_bam_bytes(&mut buf).unwrap();

        let decoded = super::super::record::BamRecord::decode(&buf).unwrap();
        assert!(decoded.qual.iter().all(|&q| q == 0xFF));
    }

    #[test]
    fn unmapped_record_serializes() {
        // Use ref_id=-1 and pos=-1 for a fully unmapped record.
        // The existing reader rejects pos=-1 (InvalidPosition), so we verify
        // serialization succeeds and the bytes have the right structure.
        let rec = OwnedBamRecord::builder(-1, -1, b"unmapped".to_vec())
            .flags(0x4)
            .seq(vec![Base::A, Base::C, Base::G])
            .build()
            .unwrap();

        let mut buf = Vec::new();
        rec.to_bam_bytes(&mut buf).unwrap();

        // Verify the fixed header bytes directly
        let ref_id = i32::from_le_bytes([buf[0], buf[1], buf[2], buf[3]]);
        let pos = i32::from_le_bytes([buf[4], buf[5], buf[6], buf[7]]);
        assert_eq!(ref_id, -1);
        assert_eq!(pos, -1);
    }

    #[test]
    fn placed_unmapped_record_roundtrips() {
        // Placed unmapped: has ref_id and pos, but flag 0x4 set
        let rec = OwnedBamRecord::builder(0, 500, b"placed".to_vec())
            .flags(0x4)
            .seq(vec![Base::A, Base::C, Base::G])
            .build()
            .unwrap();

        let mut buf = Vec::new();
        rec.to_bam_bytes(&mut buf).unwrap();

        let decoded = super::super::record::BamRecord::decode(&buf).unwrap();
        assert_eq!(decoded.n_cigar_ops, 0);
        assert_eq!(decoded.seq_len, 3);
        assert_eq!(decoded.flags & 0x4, 0x4);
    }

    // r[verify bam.owned_record.aligned_pairs]
    #[test]
    fn aligned_pairs_simple_match() {
        let rec = simple_record(); // 5M at pos 100
        let pairs: Vec<_> = rec.aligned_pairs().collect();
        assert_eq!(pairs.len(), 5);
        assert_eq!(pairs[0], (Some(0), Some(100)));
        assert_eq!(pairs[4], (Some(4), Some(104)));
    }

    #[test]
    fn aligned_pairs_with_insertion() {
        let rec = OwnedBamRecord::builder(0, 100, b"r".to_vec())
            .cigar(vec![
                CigarOp::new(CigarOpType::Match, 2),
                CigarOp::new(CigarOpType::Insertion, 1),
                CigarOp::new(CigarOpType::Match, 2),
            ])
            .seq(vec![Base::A; 5])
            .qual(vec![30; 5])
            .build()
            .unwrap();

        let pairs: Vec<_> = rec.aligned_pairs().collect();
        assert_eq!(pairs.len(), 5);
        assert_eq!(pairs[0], (Some(0), Some(100))); // M
        assert_eq!(pairs[1], (Some(1), Some(101))); // M
        assert_eq!(pairs[2], (Some(2), None)); // I — no ref pos
        assert_eq!(pairs[3], (Some(3), Some(102))); // M
        assert_eq!(pairs[4], (Some(4), Some(103))); // M
    }

    #[test]
    fn aligned_pairs_with_deletion() {
        let rec = OwnedBamRecord::builder(0, 100, b"r".to_vec())
            .cigar(vec![
                CigarOp::new(CigarOpType::Match, 2),
                CigarOp::new(CigarOpType::Deletion, 3),
                CigarOp::new(CigarOpType::Match, 2),
            ])
            .seq(vec![Base::A; 4])
            .qual(vec![30; 4])
            .build()
            .unwrap();

        let pairs: Vec<_> = rec.aligned_pairs().collect();
        assert_eq!(pairs.len(), 7); // 2M + 3D + 2M
        assert_eq!(pairs[0], (Some(0), Some(100)));
        assert_eq!(pairs[1], (Some(1), Some(101)));
        assert_eq!(pairs[2], (None, Some(102))); // D
        assert_eq!(pairs[3], (None, Some(103))); // D
        assert_eq!(pairs[4], (None, Some(104))); // D
        assert_eq!(pairs[5], (Some(2), Some(105)));
        assert_eq!(pairs[6], (Some(3), Some(106)));
    }

    #[test]
    fn aligned_pairs_with_soft_clip() {
        let rec = OwnedBamRecord::builder(0, 100, b"r".to_vec())
            .cigar(vec![
                CigarOp::new(CigarOpType::SoftClip, 2),
                CigarOp::new(CigarOpType::Match, 3),
            ])
            .seq(vec![Base::A; 5])
            .qual(vec![30; 5])
            .build()
            .unwrap();

        let pairs: Vec<_> = rec.aligned_pairs().collect();
        // Soft clips are skipped, only 3M pairs
        assert_eq!(pairs.len(), 3);
        assert_eq!(pairs[0], (Some(2), Some(100))); // qpos starts at 2 (after 2S)
        assert_eq!(pairs[2], (Some(4), Some(102)));
    }

    #[test]
    fn aligned_pairs_empty_cigar() {
        let rec = OwnedBamRecord::builder(-1, -1, b"r".to_vec()).flags(0x4).build().unwrap();
        let pairs: Vec<_> = rec.aligned_pairs().collect();
        assert!(pairs.is_empty());
    }

    // Test from_raw_bam round-trip
    #[test]
    fn from_raw_bam_preserves_all_fields() {
        let original = OwnedBamRecord::builder(1, 500, b"read1".to_vec())
            .flags(0x63)
            .mapq(42)
            .cigar(vec![CigarOp::new(CigarOpType::Match, 5)])
            .seq(vec![Base::A, Base::C, Base::G, Base::T, Base::A])
            .qual(vec![30, 31, 32, 33, 34])
            .next_ref_id(2)
            .next_pos(1000)
            .template_len(150)
            .build()
            .unwrap();

        let mut aux = AuxData::new();
        aux.set_string(*b"RG", b"grp1");
        aux.set_int(*b"NM", 3).unwrap();
        let mut original = original;
        original.aux = aux;

        // Serialize to BAM bytes
        let mut buf = Vec::new();
        original.to_bam_bytes(&mut buf).unwrap();

        // Reconstruct from raw bytes
        let reconstructed = OwnedBamRecord::from_raw_bam(&buf).unwrap();

        assert_eq!(reconstructed.ref_id, 1);
        assert_eq!(reconstructed.pos, 500);
        assert_eq!(reconstructed.flags, 0x63);
        assert_eq!(reconstructed.mapq, 42);
        assert_eq!(reconstructed.next_ref_id, 2);
        assert_eq!(reconstructed.next_pos, 1000);
        assert_eq!(reconstructed.template_len, 150);
        assert_eq!(reconstructed.qname, b"read1");
        assert_eq!(reconstructed.cigar.len(), 1);
        assert_eq!(reconstructed.cigar[0], CigarOp::new(CigarOpType::Match, 5));
        assert_eq!(reconstructed.seq.len(), 5);
        assert_eq!(reconstructed.qual, vec![30, 31, 32, 33, 34]);
        assert_eq!(reconstructed.aux.get(*b"NM"), Some(super::super::aux::AuxValue::U8(3)));
        assert_eq!(
            reconstructed.aux.get(*b"RG"),
            Some(super::super::aux::AuxValue::String(b"grp1"))
        );
    }

    // Test full round-trip: build → serialize → from_raw_bam → serialize → compare
    #[test]
    fn full_roundtrip_via_raw_bam() {
        let rec = simple_record();
        let mut buf1 = Vec::new();
        rec.to_bam_bytes(&mut buf1).unwrap();

        let reconstructed = OwnedBamRecord::from_raw_bam(&buf1).unwrap();
        let mut buf2 = Vec::new();
        reconstructed.to_bam_bytes(&mut buf2).unwrap();

        assert_eq!(buf1, buf2);
    }
}
