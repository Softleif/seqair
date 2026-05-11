//! Owned, mutable BAM record for writing and in-memory modification.
//!
//! [`OwnedBamRecord`] is the write-path entry point for the crate. Construct
//! one with [`OwnedBamRecord::builder`] (typed fields, validated at build
//! time), then hand it to [`BamWriter::write`](super::writer::BamWriter::write)
//! to serialize. [`to_bam_bytes`](OwnedBamRecord::to_bam_bytes) exposes the
//! same encoder for callers that need raw BAM bytes without going through a
//! writer.
//!
//! For the read path — bulk decode of a region's records for pileup or
//! per-read iteration — use [`RecordStore`](super::record_store::RecordStore)
//! instead. The store streams raw BAM bytes straight into slab buffers
//! without allocating per-record owned structs and is what every production
//! reader path uses. `OwnedBamRecord` exists separately because the write
//! path needs typed, mutable fields with builder ergonomics; pushing
//! pre-encoded bytes into slabs is the wrong shape for that.

use super::aux_data::{AuxData, AuxDataError};
use super::cigar::CigarOp;
use super::seq;
use crate::io::reg2bin;
use seqair_types::{BamFlags, Base, BaseQuality, Pos0};
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

    /// Sequence length exceeds `i32::MAX` (BAM `l_seq` is i32).
    #[error("sequence length {len} exceeds i32::MAX")]
    SeqLengthOverflow { len: usize },

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
    #[error("aux data error")]
    AuxData {
        #[from]
        source: AuxDataError,
    },
}

/// Encode `Option<Pos0>` to a BAM-wire `i32`. `None` → `-1`.
#[inline]
fn pos_to_bam_i32(p: Option<Pos0>) -> i32 {
    p.map_or(-1, Pos0::as_i32)
}

// r[impl bam.owned_record.fields]
/// An owned, mutable BAM record for writing and in-memory modification.
#[derive(Debug, Clone)]
pub struct OwnedBamRecord {
    /// Reference sequence index. `-1` is the canonical "unmapped" sentinel
    /// per [SAM1] §1.4 — it is preserved here as a raw `i32` (not wrapped in
    /// an `Option`) because reading code uses it as a header-table index and
    /// every check site already handles the -1 case explicitly.
    pub ref_id: i32,
    /// 0-based leftmost position. `None` for unmapped/unavailable (BAM
    /// wire `-1`); `Some(p)` for any in-range mapping. Construction via
    /// [`Pos0`] makes overflow unconstructable.
    pub pos: Option<Pos0>,
    pub mapq: u8,
    pub flags: BamFlags,
    /// Mate's reference sequence index. `-1` for unavailable (same rationale
    /// as `ref_id`).
    pub next_ref_id: i32,
    /// Mate's 0-based position. `None` for unavailable (BAM wire `-1`).
    pub next_pos: Option<Pos0>,
    pub template_len: i32,
    pub qname: Vec<u8>,
    pub cigar: Vec<CigarOp>,
    pub seq: Vec<Base>,
    // r[impl types.base_quality.field_type]
    pub qual: Vec<BaseQuality>,
    pub aux: AuxData,
}

/// Builder for constructing [`OwnedBamRecord`] values.
pub struct OwnedBamRecordBuilder {
    ref_id: i32,
    pos: Option<Pos0>,
    mapq: u8,
    flags: BamFlags,
    next_ref_id: i32,
    next_pos: Option<Pos0>,
    template_len: i32,
    qname: Vec<u8>,
    cigar: Vec<CigarOp>,
    seq: Vec<Base>,
    qual: Vec<BaseQuality>,
    aux: AuxData,
}

// r[impl bam.owned_record.builder]
impl OwnedBamRecordBuilder {
    pub fn flags(mut self, flags: BamFlags) -> Self {
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
    pub fn qual(mut self, qual: Vec<BaseQuality>) -> Self {
        self.qual = qual;
        self
    }
    pub fn next_ref_id(mut self, v: i32) -> Self {
        self.next_ref_id = v;
        self
    }
    pub fn next_pos(mut self, v: Option<Pos0>) -> Self {
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
    /// Start building a record with the required fields.
    ///
    /// Pass `None` for `pos` to mark the record as unmapped (BAM wire `-1`).
    pub fn builder(ref_id: i32, pos: Option<Pos0>, qname: Vec<u8>) -> OwnedBamRecordBuilder {
        OwnedBamRecordBuilder {
            ref_id,
            pos,
            mapq: 0,
            flags: BamFlags::empty(),
            next_ref_id: -1,
            next_pos: None,
            template_len: 0,
            qname,
            cigar: Vec::new(),
            seq: Vec::new(),
            qual: Vec::new(),
            aux: AuxData::new(),
        }
    }

    // r[impl bam.owned_record.flag_methods]
    pub fn set_flags(&mut self, flags: BamFlags) {
        self.flags = flags;
    }

    pub fn is_unmapped(&self) -> bool {
        self.flags.is_unmapped()
    }

    pub fn is_reverse(&self) -> bool {
        self.flags.is_reverse()
    }

    // r[impl bam.owned_record.end_pos]
    /// Compute the 0-based exclusive end position from pos + reference-consuming CIGAR ops.
    /// Returns `None` if the record is unmapped (`pos.is_none()`).
    pub fn end_pos(&self) -> Option<Pos0> {
        let pos = self.pos?;
        let ref_len: u64 =
            self.cigar.iter().filter(|op| op.consumes_ref()).map(|op| u64::from(op.len())).sum();
        let end = pos.as_u64().saturating_add(ref_len);
        let clamped = end.min(Pos0::max_value().as_u64());
        #[expect(
            clippy::cast_possible_truncation,
            reason = "clamped <= i32::MAX < u32::MAX; fits in u32"
        )]
        Pos0::new(clamped as u32)
    }

    // r[impl bam.owned_record.bin]
    /// Compute the BAI bin value from pos and `end_pos` (BAI scheme: `min_shift=14`, depth=5).
    /// For unmapped records (`pos.is_none()`), returns the unmapped bin (`reg2bin(0, 1, ...)`).
    pub fn bin(&self) -> u16 {
        let beg = self.pos.map(Pos0::as_u64).unwrap_or(0);
        let end = self.end_pos().map_or(beg.saturating_add(1), |p| {
            let v = p.as_u64();
            if v <= beg { beg.saturating_add(1) } else { v }
        });
        #[expect(
            clippy::cast_possible_truncation,
            reason = "BAI bin numbers with min_shift=14, depth=5 are bounded by 37449, fits in u16"
        )]
        let bin = reg2bin(beg, end, 14, 5) as u16;
        bin
    }

    /// Query-consuming length of the CIGAR.
    pub fn cigar_query_len(&self) -> u64 {
        self.cigar.iter().filter(|op| op.consumes_query()).map(|op| u64::from(op.len())).sum()
    }

    // r[impl bam.owned_record.set_alignment]
    /// Replace alignment (pos + cigar). Validates CIGAR query length == seq length for mapped reads.
    pub fn set_alignment(
        &mut self,
        pos: Option<Pos0>,
        cigar: Vec<CigarOp>,
    ) -> Result<(), OwnedRecordError> {
        if cigar.len() > 65535 {
            return Err(OwnedRecordError::CigarCountOverflow { count: cigar.len() });
        }
        // For mapped reads, validate CIGAR query length matches seq length
        if !self.is_unmapped() && !self.seq.is_empty() {
            let query_len: u64 =
                cigar.iter().filter(|op| op.consumes_query()).map(|op| u64::from(op.len())).sum();
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
    pub fn set_qual(&mut self, qual: Vec<BaseQuality>) -> Result<(), OwnedRecordError> {
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
    // r[impl bam.owned_record.to_bam_bin_field]
    /// Serialize to BAM binary format by appending into `buf`.
    ///
    /// Does NOT include the 4-byte `block_size` prefix — that is the writer's responsibility.
    pub fn to_bam_bytes(&self, buf: &mut Vec<u8>) -> Result<(), OwnedRecordError> {
        // Validate field limits. `pos`/`next_pos` are constructively bounded
        // by `Option<Pos0>` (Pos0 ≤ i32::MAX) so no overflow check is needed.
        let name_len = self.qname.len();
        if name_len > 254 {
            return Err(OwnedRecordError::QnameTooLong { len: name_len });
        }
        if self.cigar.len() > 65535 {
            return Err(OwnedRecordError::CigarCountOverflow { count: self.cigar.len() });
        }
        if self.seq.len() > i32::MAX as usize {
            return Err(OwnedRecordError::SeqLengthOverflow { len: self.seq.len() });
        }

        #[expect(
            clippy::cast_possible_truncation,
            reason = "qname.len() ≤ 254 (validated above); +1 fits in u8"
        )]
        let l_read_name = (self.qname.len() as u8).saturating_add(1); // +1 for NUL
        #[expect(
            clippy::cast_possible_truncation,
            reason = "cigar.len() ≤ 65535 (validated above); fits in u16"
        )]
        let n_cigar_op = self.cigar.len() as u16;
        #[expect(
            clippy::cast_possible_truncation,
            clippy::cast_possible_wrap,
            reason = "seq.len() ≤ i32::MAX (validated above); fits in i32"
        )]
        let l_seq = self.seq.len() as i32;

        super::record::encode_fixed_header(
            buf,
            &super::record::FixedHeaderFields {
                ref_id: self.ref_id,
                pos: pos_to_bam_i32(self.pos),
                bin: self.bin(),
                mapq: self.mapq,
                l_read_name,
                flags: self.flags.raw(),
                n_cigar_op,
                l_seq,
                next_ref_id: self.next_ref_id,
                next_pos: pos_to_bam_i32(self.next_pos),
                template_len: self.template_len,
            },
        );

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
            buf.extend_from_slice(BaseQuality::slice_to_bytes(&self.qual));
        }

        // Auxiliary data
        buf.extend_from_slice(self.aux.as_bytes());

        Ok(())
    }
}

#[cfg(test)]
#[allow(clippy::arithmetic_side_effects, reason = "test arithmetic on known small values")]
mod tests {
    use super::super::cigar::CigarOpType;
    use super::*;
    use seqair_types::Pos0;

    fn simple_record() -> OwnedBamRecord {
        OwnedBamRecord::builder(0, Some(Pos0::new(100).unwrap()), b"read1".to_vec())
            .flags(BamFlags::empty()) // mapped
            .mapq(30)
            .cigar(vec![CigarOp::new(CigarOpType::Match, 5)])
            .seq(vec![Base::A, Base::C, Base::G, Base::T, Base::A])
            .qual([30, 31, 32, 33, 34].map(BaseQuality::from_byte).to_vec())
            .build()
            .unwrap()
    }

    use super::super::test_util::decode_into_store;

    // r[verify bam.owned_record.to_bam_bytes]
    // r[verify bam.owned_record.test_roundtrip_bytes]
    #[test]
    fn serialize_and_decode_roundtrip() {
        let rec = simple_record();
        let mut buf = Vec::new();
        rec.to_bam_bytes(&mut buf).unwrap();

        let store = decode_into_store(&buf);
        let decoded = store.record(0);

        assert_eq!(decoded.tid, 0);
        assert_eq!(*decoded.pos, 100);
        assert_eq!(decoded.mapq, 30);
        assert_eq!(decoded.flags, BamFlags::empty());
        assert_eq!(decoded.seq_len, 5);
        assert_eq!(decoded.n_cigar_ops, 1);
        assert_eq!(store.qname(0), b"read1");
    }

    #[test]
    fn serialize_with_aux_tags() {
        let mut aux = AuxData::new();
        aux.set_int(*b"NM", 3).unwrap();
        aux.set_string(*b"RG", b"grp1");

        let rec = OwnedBamRecord::builder(0, Some(Pos0::new(100).unwrap()), b"r1".to_vec())
            .cigar(vec![CigarOp::new(CigarOpType::Match, 3)])
            .seq(vec![Base::A, Base::C, Base::G])
            .qual([30, 31, 32].map(BaseQuality::from_byte).to_vec())
            .aux(aux)
            .build()
            .unwrap();

        let mut buf = Vec::new();
        rec.to_bam_bytes(&mut buf).unwrap();

        let store = decode_into_store(&buf);
        let aux = store.aux(0);
        let nm = super::super::aux::find_tag(aux, *b"NM");
        assert_eq!(nm, Some(super::super::aux::AuxValue::U8(3)));

        let rg = super::super::aux::find_tag(aux, *b"RG");
        assert_eq!(rg, Some(super::super::aux::AuxValue::String(b"grp1")));
    }

    // r[verify bam.owned_record.to_bam_bin_field]
    #[test]
    fn bin_recomputed_after_modification() {
        let mut rec = simple_record(); // pos=100, 5M
        let bin_before = rec.bin();

        // Move to a completely different position — bin must change
        rec.set_alignment(
            Some(Pos0::new(10_000_000).unwrap()),
            vec![CigarOp::new(CigarOpType::Match, 5)],
        )
        .unwrap();
        let bin_after = rec.bin();
        assert_ne!(bin_before, bin_after, "bin must change when pos changes");

        // Serialize and verify the bin in the bytes matches the new value
        let mut buf = Vec::new();
        rec.to_bam_bytes(&mut buf).unwrap();
        // bin_mq_nl is at offset 8 in the 32-byte header (bytes 8-11)
        let bin_mq_nl = u32::from_le_bytes([buf[8], buf[9], buf[10], buf[11]]);
        let serialized_bin = (bin_mq_nl >> 16) as u16;
        assert_eq!(serialized_bin, bin_after, "serialized bin must match current bin()");
    }

    // r[verify bam.owned_record.end_pos]
    #[test]
    fn end_pos_simple() {
        let rec = simple_record();
        assert_eq!(rec.end_pos(), Some(Pos0::new(105).unwrap())); // pos=100 + 5M
    }

    #[test]
    fn end_pos_with_insertion() {
        let rec = OwnedBamRecord::builder(0, Some(Pos0::new(100).unwrap()), b"r".to_vec())
            .cigar(vec![
                CigarOp::new(CigarOpType::Match, 3),
                CigarOp::new(CigarOpType::Insertion, 2),
                CigarOp::new(CigarOpType::Match, 3),
            ])
            .seq(vec![Base::A; 8]) // 3M + 2I + 3M = 8 query bases
            .qual(vec![BaseQuality::from_byte(30); 8])
            .build()
            .unwrap();
        // Only M ops consume ref: 3 + 3 = 6
        assert_eq!(rec.end_pos(), Some(Pos0::new(106).unwrap()));
    }

    // r[verify bam.owned_record.qname_limit]
    #[test]
    fn qname_too_long_rejected() {
        let long_name = vec![b'A'; 255];
        let result = OwnedBamRecord::builder(0, None, long_name).build();
        assert!(result.is_err());
    }

    #[test]
    fn qname_at_limit_accepted() {
        let name = vec![b'A'; 254];
        let result = OwnedBamRecord::builder(0, None, name).build();
        assert!(result.is_ok());
    }

    // r[verify bam.owned_record.set_alignment]
    // r[verify bam.owned_record.test_modification]
    #[test]
    fn set_alignment_validates_cigar_seq_mismatch() {
        let mut rec = simple_record(); // seq len = 5
        let result = rec.set_alignment(
            Some(Pos0::new(200).unwrap()),
            vec![CigarOp::new(CigarOpType::Match, 10)], // query len = 10 != 5
        );
        assert!(result.is_err());
    }

    #[test]
    fn set_alignment_succeeds_when_matching() {
        let mut rec = simple_record(); // seq len = 5
        let result = rec.set_alignment(
            Some(Pos0::new(200).unwrap()),
            vec![CigarOp::new(CigarOpType::Match, 5)],
        );
        assert!(result.is_ok());
        assert_eq!(rec.pos, Some(Pos0::new(200).unwrap()));
        assert_eq!(rec.end_pos(), Some(Pos0::new(205).unwrap()));
    }

    #[test]
    fn empty_qual_fills_with_0xff() {
        let rec = OwnedBamRecord::builder(0, Some(Pos0::new(100).unwrap()), b"r".to_vec())
            .cigar(vec![CigarOp::new(CigarOpType::Match, 3)])
            .seq(vec![Base::A, Base::C, Base::G])
            // no qual set
            .build()
            .unwrap();

        let mut buf = Vec::new();
        rec.to_bam_bytes(&mut buf).unwrap();

        let store = decode_into_store(&buf);
        assert!(store.qual(0).iter().all(|&q| q == BaseQuality::UNAVAILABLE));
    }

    #[test]
    fn unmapped_record_serializes() {
        // Use ref_id=-1 and pos=None for a fully unmapped record.
        // The existing reader rejects pos=-1 (InvalidPosition), so we verify
        // serialization succeeds and the bytes have the right structure.
        let rec = OwnedBamRecord::builder(-1, None, b"unmapped".to_vec())
            .flags(BamFlags::from(0x4))
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
        let rec = OwnedBamRecord::builder(0, Some(Pos0::new(500).unwrap()), b"placed".to_vec())
            .flags(BamFlags::from(0x4))
            .seq(vec![Base::A, Base::C, Base::G])
            .build()
            .unwrap();

        let mut buf = Vec::new();
        rec.to_bam_bytes(&mut buf).unwrap();

        let store = decode_into_store(&buf);
        let decoded = store.record(0);
        assert_eq!(decoded.n_cigar_ops, 0);
        assert_eq!(decoded.seq_len, 3);
        assert_eq!(decoded.flags & BamFlags::from(0x4), BamFlags::from(0x4));
    }

    // r[verify cigar.aligned_pairs.owned_record]
    #[test]
    fn aligned_pairs_simple_match() {
        use super::super::aligned_pairs::{AlignedPair, MatchKind};
        let rec = simple_record(); // 5M at pos 100
        let pairs: Vec<_> = rec.aligned_pairs().unwrap().collect();
        assert_eq!(pairs.len(), 5);
        assert_eq!(
            pairs[0],
            AlignedPair::Match { qpos: 0, rpos: Pos0::new(100).unwrap(), kind: MatchKind::Match }
        );
        assert_eq!(
            pairs[4],
            AlignedPair::Match { qpos: 4, rpos: Pos0::new(104).unwrap(), kind: MatchKind::Match }
        );
    }

    #[test]
    fn aligned_pairs_with_insertion() {
        use super::super::aligned_pairs::{AlignedPair, MatchKind};
        let rec = OwnedBamRecord::builder(0, Some(Pos0::new(100).unwrap()), b"r".to_vec())
            .cigar(vec![
                CigarOp::new(CigarOpType::Match, 2),
                CigarOp::new(CigarOpType::Insertion, 1),
                CigarOp::new(CigarOpType::Match, 2),
            ])
            .seq(vec![Base::A; 5])
            .qual(vec![BaseQuality::from_byte(30); 5])
            .build()
            .unwrap();

        let pairs: Vec<_> = rec.aligned_pairs().unwrap().collect();
        let m = |q, r| AlignedPair::Match {
            qpos: q,
            rpos: Pos0::new(r).unwrap(),
            kind: MatchKind::Match,
        };
        assert_eq!(pairs.len(), 5);
        assert_eq!(pairs[0], m(0, 100));
        assert_eq!(pairs[1], m(1, 101));
        assert_eq!(pairs[2], AlignedPair::Insertion { qpos: 2, insert_len: 1 });
        assert_eq!(pairs[3], m(3, 102));
        assert_eq!(pairs[4], m(4, 103));
    }

    #[test]
    fn aligned_pairs_with_deletion() {
        use super::super::aligned_pairs::{AlignedPair, MatchKind};
        let rec = OwnedBamRecord::builder(0, Some(Pos0::new(100).unwrap()), b"r".to_vec())
            .cigar(vec![
                CigarOp::new(CigarOpType::Match, 2),
                CigarOp::new(CigarOpType::Deletion, 3),
                CigarOp::new(CigarOpType::Match, 2),
            ])
            .seq(vec![Base::A; 4])
            .qual(vec![BaseQuality::from_byte(30); 4])
            .build()
            .unwrap();

        let pairs: Vec<_> = rec.aligned_pairs().unwrap().collect();
        let m = |q, r| AlignedPair::Match {
            qpos: q,
            rpos: Pos0::new(r).unwrap(),
            kind: MatchKind::Match,
        };
        assert_eq!(pairs.len(), 5); // 2M + D(summary) + 2M = 5
        assert_eq!(pairs[0], m(0, 100));
        assert_eq!(pairs[1], m(1, 101));
        // Deletion summary: rpos=102, del_len=3
        assert_eq!(pairs[2], AlignedPair::Deletion { rpos: Pos0::new(102).unwrap(), del_len: 3 });
        assert_eq!(pairs[3], m(2, 105));
        assert_eq!(pairs[4], m(3, 106));
    }

    #[test]
    fn aligned_pairs_with_soft_clip() {
        use super::super::aligned_pairs::{AlignedPair, MatchKind};
        let rec = OwnedBamRecord::builder(0, Some(Pos0::new(100).unwrap()), b"r".to_vec())
            .cigar(vec![
                CigarOp::new(CigarOpType::SoftClip, 2),
                CigarOp::new(CigarOpType::Match, 3),
            ])
            .seq(vec![Base::A; 5])
            .qual(vec![BaseQuality::from_byte(30); 5])
            .build()
            .unwrap();

        let pairs: Vec<_> = rec.aligned_pairs().unwrap().collect();
        let m = |q, r| AlignedPair::Match {
            qpos: q,
            rpos: Pos0::new(r).unwrap(),
            kind: MatchKind::Match,
        };
        // Soft clips are skipped by default, only 3 Match pairs
        assert_eq!(pairs.len(), 3);
        assert_eq!(pairs[0], m(2, 100));
        assert_eq!(pairs[2], m(4, 102));
    }

    #[test]
    fn aligned_pairs_empty_cigar() {
        let rec = OwnedBamRecord::builder(-1, None, b"r".to_vec())
            .flags(BamFlags::from(0x4))
            .build()
            .unwrap();
        let pairs: Vec<_> = rec.aligned_pairs().unwrap().collect();
        assert!(pairs.is_empty());
    }

    // r[verify cigar.aligned_pairs.owned_record]
    #[test]
    fn aligned_pairs_rejects_unmapped_with_cigar() {
        // Degenerate but constructible: pos = None + non-empty CIGAR. Walking
        // would produce nonsense rpos values anchored at Pos0::ZERO, so we
        // refuse rather than silently corrupt the output.
        let rec = OwnedBamRecord::builder(-1, None, b"r".to_vec())
            .flags(BamFlags::from(0x4))
            .cigar(vec![CigarOp::new(CigarOpType::Match, 5)])
            .seq(vec![Base::A; 5])
            .qual(vec![BaseQuality::from_byte(30); 5])
            .build()
            .unwrap();
        let result = rec.aligned_pairs();
        assert!(matches!(
            result,
            Err(super::super::aligned_pairs::AlignedPairsError::UnmappedWithCigar { cigar_ops: 1 })
        ));
    }

    // r[verify cigar.aligned_pairs.with_read.owned_record]
    #[test]
    fn aligned_pairs_with_read_attaches_query_base() {
        // 3M record with seq [A, C, G] — verify the helper attaches per-pos
        // bases without going through the SlimRecord/RecordStore round-trip.
        use super::super::aligned_pairs_view::AlignedPairWithRead;
        let rec = OwnedBamRecord::builder(0, Some(Pos0::new(50).unwrap()), b"r".to_vec())
            .flags(BamFlags::empty())
            .cigar(vec![CigarOp::new(CigarOpType::Match, 3)])
            .seq(vec![Base::A, Base::C, Base::G])
            .qual([30, 31, 32].map(BaseQuality::from_byte).to_vec())
            .build()
            .unwrap();

        let events: Vec<_> = rec.aligned_pairs_with_read().unwrap().collect();
        assert_eq!(events.len(), 3);
        match events[0] {
            AlignedPairWithRead::Match { qpos: 0, query: Base::A, .. } => {}
            other => panic!("expected Match{{qpos=0, query=A}}, got {other:?}"),
        }
        match events[2] {
            AlignedPairWithRead::Match { qpos: 2, query: Base::G, .. } => {}
            other => panic!("expected Match{{qpos=2, query=G}}, got {other:?}"),
        }
    }

    // r[verify cigar.aligned_pairs.with_read.owned_record]
    #[test]
    fn aligned_pairs_with_read_propagates_unmapped_error() {
        // pos = None + non-empty CIGAR → UnmappedWithCigar surfaces through
        // the helper just like through the bare aligned_pairs().
        let rec = OwnedBamRecord::builder(-1, None, b"r".to_vec())
            .flags(BamFlags::from(0x4))
            .cigar(vec![CigarOp::new(CigarOpType::Match, 3)])
            .seq(vec![Base::A; 3])
            .qual(vec![BaseQuality::from_byte(30); 3])
            .build()
            .unwrap();
        let result = rec.aligned_pairs_with_read();
        assert!(matches!(
            result,
            Err(super::super::aligned_pairs::AlignedPairsError::UnmappedWithCigar { cigar_ops: 1 })
        ));
    }

    // r[verify cigar.aligned_pairs.owned_record]
    #[test]
    fn aligned_pairs_accepts_unmapped_with_empty_cigar() {
        // The "fully unmapped" case: pos = None AND cigar = empty. Iteration
        // yields nothing — no rpos anchoring needed.
        let rec = OwnedBamRecord::builder(-1, None, b"r".to_vec())
            .flags(BamFlags::from(0x4))
            .build()
            .unwrap();
        let pairs: Vec<_> = rec.aligned_pairs().unwrap().collect();
        assert!(pairs.is_empty());
    }
}
