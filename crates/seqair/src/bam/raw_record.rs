//! [`RawBamRecord`]: one BAM record kept in its wire form, validated once.
//!
//! The streaming readers (`r[bam.stream]`, `r[sam.stream]`) yield these, and
//! [`BamWriter::write_raw`](super::BamWriter::write_raw) writes one back byte
//! for byte. Nothing is decoded until an accessor asks for it, and nothing is
//! lost: unplaced records, IUPAC bases, every aux tag and its integer width,
//! and the stored `bin` all survive a read and a write unchanged.

use super::{
    BamHeader, CigarOp, DecodeError, MateFields,
    aux::{Aux, AuxValue, find_tag},
    cigar::CIGAR_S,
    record::{self, read2},
    seq,
};
use crate::reader::Tid;
use seqair_types::{BamFlags, Base, BaseQuality, MappingQuality, Pos0};
use std::fmt;

/// Fixed part of a BAM record, before the read name.
const FIXED_LEN: usize = 32;

/// A minimal valid record: unplaced, unmapped, name `*`, nothing else. What a
/// [`RawBamRecord`] holds after [`RawBamRecord::default`] or a failed refill.
const EMPTY_RECORD: [u8; FIXED_LEN + 2] = {
    let mut r = [0u8; FIXED_LEN + 2];
    let none = (-1i32).to_le_bytes();
    let mut i = 0;
    while i < 4 {
        r[i] = none[i]; // refID
        r[4 + i] = none[i]; // pos
        r[20 + i] = none[i]; // next refID
        r[24 + i] = none[i]; // next pos
        i += 1;
    }
    r[8] = 2; // l_read_name: "*\0"
    r[9] = 255; // MAPQ unavailable
    // bin 4680 = reg2bin(-1, 0), what htslib writes for unplaced records.
    let bin = 4680u16.to_le_bytes();
    r[10] = bin[0];
    r[11] = bin[1];
    r[14] = 4; // FLAG: unmapped
    r[FIXED_LEN] = b'*';
    r
};

// r[impl bam.raw_record]
// r[impl bam.raw_record.aux_lazy]
// r[impl bam.raw_record.bin_verbatim]
/// One BAM alignment record, kept in its wire form: the bytes after BAM's
/// `block_size` prefix.
///
/// # Invariants
///
/// Every value of this type has passed [`Self::from_payload`]'s checks
/// against some header (`r[bam.raw_record.validation]`): the variable-length
/// fields lie inside the payload, the read name is NUL-terminated, the
/// reference ids name targets of that header (or are `-1`), the positions are
/// `-1` or non-negative, and a mapped record's CIGAR covers its sequence. So
/// no accessor can fail or panic.
///
/// What is *not* checked, as htslib does not check it on read: the aux block
/// (htslib passes a malformed one through untouched and only fails when it
/// formats it; [`Aux::iter_checked`] and [`validate_aux`](super::aux::validate_aux)
/// report it), and the stored `bin`, which is kept as the file has it (htslib
/// recomputes it on read, so `samtools view -b` may rewrite it). A read name
/// without its NUL terminator is rejected, where htslib repairs it.
#[derive(Clone, PartialEq, Eq)]
pub struct RawBamRecord {
    payload: Vec<u8>,
    tid: Option<Tid>,
    pos: Option<Pos0>,
    mate: MateFields,
    /// End of the read name (NUL included) = start of the CIGAR.
    cigar_start: usize,
    cigar_end: usize,
    seq_end: usize,
    qual_end: usize,
}

/// Why [`RawBamRecord::refill`] failed: the caller's fill, or validation.
#[derive(Debug)]
pub(crate) enum Refill<E> {
    Fill(E),
    Invalid(DecodeError),
}

impl Default for RawBamRecord {
    /// An unplaced, unmapped record named `*` with no sequence.
    fn default() -> Self {
        Self {
            payload: EMPTY_RECORD.to_vec(),
            tid: None,
            pos: None,
            mate: MateFields::default(),
            cigar_start: EMPTY_RECORD.len(),
            cigar_end: EMPTY_RECORD.len(),
            seq_end: EMPTY_RECORD.len(),
            qual_end: EMPTY_RECORD.len(),
        }
    }
}

impl fmt::Debug for RawBamRecord {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("RawBamRecord")
            .field("qname", &String::from_utf8_lossy(self.qname()))
            .field("flags", &self.flags())
            .field("tid", &self.tid)
            .field("pos", &self.pos)
            .field("mapq", &self.mapq())
            .field("payload_len", &self.payload.len())
            .finish_non_exhaustive()
    }
}

impl RawBamRecord {
    // r[impl bam.raw_record.validation]
    /// Validate `payload` — one record's bytes, without the `block_size`
    /// prefix — against `header`, and keep it.
    ///
    /// # Errors
    /// A [`DecodeError`] naming the first check that failed.
    pub fn from_payload(payload: Vec<u8>, header: &BamHeader) -> Result<Self, DecodeError> {
        let mut record = Self { payload, ..Self::default() };
        record.revalidate(header)?;
        Ok(record)
    }

    /// Replace this record with the payload `fill` writes into the (cleared)
    /// buffer, reusing its allocation, and validate it against `header`. On
    /// error the record is [`Self::default`] again, so the invariants hold
    /// whatever happens. For callers outside the crate, reuse goes through
    /// [`Self::into_payload`] and [`Self::from_payload`].
    pub(crate) fn refill<E>(
        &mut self,
        header: &BamHeader,
        fill: impl FnOnce(&mut Vec<u8>) -> Result<(), E>,
    ) -> Result<(), Refill<E>> {
        self.payload.clear();
        if let Err(e) = fill(&mut self.payload) {
            self.reset();
            return Err(Refill::Fill(e));
        }
        self.revalidate(header).map_err(Refill::Invalid)
    }

    /// Check `self.payload` and fill in the cached layout. On error, reset.
    fn revalidate(&mut self, header: &BamHeader) -> Result<(), DecodeError> {
        let result = self.validate(header);
        if result.is_err() {
            self.reset();
        }
        result
    }

    fn validate(&mut self, header: &BamHeader) -> Result<(), DecodeError> {
        let h = record::parse_header(&self.payload)?;
        let refs = h.refs.resolve(header)?;
        let name = self.payload.get(FIXED_LEN..h.var_start).unwrap_or_default();
        match name.split_last() {
            Some((0, _)) => {}
            Some(_) => return Err(DecodeError::QnameNotTerminated),
            None => return Err(DecodeError::EmptyQname),
        }
        check_cigar_against_seq(&self.payload, &h, refs.tid.is_some() && h.pos.is_some())?;
        self.tid = refs.tid;
        self.pos = h.pos;
        self.mate = refs.mate_fields;
        self.cigar_start = h.var_start;
        self.cigar_end = h.cigar_end;
        self.seq_end = h.seq_end;
        self.qual_end = h.qual_end;
        Ok(())
    }

    /// Whether the stored CIGAR is the placeholder BAM writes for a record
    /// with more than 65535 operations — `<seq_len>S<ref_len>N` with the real
    /// CIGAR in a `CG:B:I` tag (`[SAM1] §4.2.2`), as htslib recognises it
    /// (`r[bam.raw_record.long_cigar]`).
    pub fn has_cg_placeholder(&self) -> bool {
        self.long_cigar().is_some()
    }

    /// The real CIGAR, in BAM's packed form, of a record that carries a long
    /// one in its `CG` tag (`r[bam.raw_record.long_cigar]`).
    pub(crate) fn long_cigar(&self) -> Option<&[u8]> {
        long_cigar(
            self.cigar_bytes(),
            self.seq_len(),
            self.tid.is_some() && self.pos.is_some(),
            self.aux(),
        )
    }

    fn reset(&mut self) {
        let mut payload = std::mem::take(&mut self.payload);
        payload.clear();
        payload.extend_from_slice(&EMPTY_RECORD);
        *self = Self { payload, ..Self::default() };
    }

    /// The record's bytes as BAM stores them, without `block_size`.
    pub fn as_payload(&self) -> &[u8] {
        &self.payload
    }

    /// Give up the payload bytes.
    pub fn into_payload(self) -> Vec<u8> {
        self.payload
    }

    /// The reference, or `None` for an unplaced record.
    pub fn tid(&self) -> Option<Tid> {
        self.tid
    }

    /// The 0-based leftmost position, or `None` (BAM's `-1`).
    pub fn pos(&self) -> Option<Pos0> {
        self.pos
    }

    /// The mate's reference and position, and the template length.
    pub fn mate_fields(&self) -> MateFields {
        self.mate
    }

    pub fn flags(&self) -> BamFlags {
        BamFlags::from(u16::from_le_bytes(read2(&self.payload, 14)))
    }

    pub fn mapq(&self) -> MappingQuality {
        MappingQuality::from_byte(self.payload.get(9).copied().unwrap_or(255))
    }

    /// The BAI bin as the file stores it, which is not checked against the
    /// alignment.
    pub fn bin(&self) -> u16 {
        u16::from_le_bytes(read2(&self.payload, 10))
    }

    /// The read name, without its NUL terminator, up to its first NUL (as
    /// htslib reads it).
    pub fn qname(&self) -> &[u8] {
        let name = self.payload.get(FIXED_LEN..self.cigar_start).unwrap_or_default();
        let len = crate::io::text_scan::find_byte(name, 0).unwrap_or(name.len());
        name.get(..len).unwrap_or_default()
    }

    /// The whole stored read name without its NUL terminator: every byte up
    /// to it, NULs inside the name included, as `sam_format1` prints it.
    pub(crate) fn qname_field(&self) -> &[u8] {
        let end = self.cigar_start.saturating_sub(1);
        self.payload.get(FIXED_LEN..end).unwrap_or_default()
    }

    /// The number of CIGAR operations stored (the placeholder ones for a
    /// record whose real CIGAR is in a `CG` tag).
    pub fn cigar_len(&self) -> usize {
        self.cigar_bytes().len() / 4
    }

    /// The stored CIGAR operations.
    pub fn cigar(&self) -> impl ExactSizeIterator<Item = CigarOp> + '_ {
        self.cigar_bytes()
            .as_chunks::<4>()
            .0
            .iter()
            .map(|op| CigarOp::from_bam_u32(u32::from_le_bytes(*op)))
    }

    /// The stored CIGAR in BAM's packed form (little-endian `u32` per op).
    pub fn cigar_bytes(&self) -> &[u8] {
        self.payload.get(self.cigar_start..self.cigar_end).unwrap_or_default()
    }

    /// The last reference position the alignment covers, by its CIGAR
    /// (`r[bam.record.end_pos]`); `None` for an unplaced record.
    pub fn end_pos(&self) -> Option<Pos0> {
        record::compute_end_pos_from_raw(&self.payload)
    }

    /// The number of bases in the sequence (0 when SEQ is `*`).
    pub fn seq_len(&self) -> usize {
        self.qual_end.saturating_sub(self.seq_end)
    }

    /// The sequence in BAM's 4-bit packed form, two bases per byte.
    pub fn packed_seq(&self) -> &[u8] {
        self.payload.get(self.cigar_end..self.seq_end).unwrap_or_default()
    }

    /// Append the sequence as ASCII, every 4-bit code kept
    /// (`=ACMGRSVTWYHKDBN`).
    pub fn decode_seq_into(&self, out: &mut Vec<u8>) {
        let start = out.len();
        out.resize(start.saturating_add(self.seq_len()), 0);
        if let Some(dst) = out.get_mut(start..) {
            seq::decode_iupac_into(self.packed_seq(), dst);
        }
    }

    /// The sequence as [`Base`]s: `A`/`C`/`G`/`T`, and everything else
    /// (ambiguity codes, `=`) [`Base::Unknown`].
    pub fn bases(&self) -> Vec<Base> {
        seq::decode_bases(self.packed_seq(), self.seq_len())
    }

    /// Base qualities, or `None` when the record has none (no sequence, or
    /// QUAL stored as `0xFF`, SAM's `*`).
    pub fn qual(&self) -> Option<&[BaseQuality]> {
        let qual = self.payload.get(self.seq_end..self.qual_end).unwrap_or_default();
        match qual.first() {
            None | Some(&0xFF) => None,
            Some(_) => Some(bytemuck::cast_slice(qual)),
        }
    }

    /// The aux tags, as stored. Not validated: use [`Aux::iter_checked`] to
    /// see a malformed tag as an error rather than the end of the tags.
    pub fn aux(&self) -> Aux<'_> {
        Aux::new(self.payload.get(self.qual_end..).unwrap_or_default())
    }

    /// Set FLAG in place.
    pub fn set_flags(&mut self, flags: BamFlags) {
        if let Some(dst) = self.payload.get_mut(14..16) {
            dst.copy_from_slice(&flags.raw().to_le_bytes());
        }
    }

    /// Set MAPQ in place.
    pub fn set_mapq(&mut self, mapq: MappingQuality) {
        if let Some(dst) = self.payload.get_mut(9) {
            *dst = mapq.as_byte();
        }
    }
}

// r[impl bam.raw_record.long_cigar]
/// See [`RawBamRecord::long_cigar`]. htslib's `bam_tag2cigar` test: the first
/// op is exactly `<seq_len>S`, the record is placed, and the first `CG` tag is
/// a `B:I`/`B:i` array of at least `n_cigar` and under 2^29 operations.
fn long_cigar<'a>(cigar: &[u8], seq_len: usize, placed: bool, aux: Aux<'a>) -> Option<&'a [u8]> {
    let first = u32::from_le_bytes(*cigar.as_chunks::<4>().0.first()?);
    let placeholder = u32::try_from(seq_len).ok()?.checked_mul(16)? | u32::from(CIGAR_S);
    if first != placeholder || !placed {
        return None;
    }
    let (Some(AuxValue::ArrayU32(ops)) | Some(AuxValue::ArrayI32(ops))) =
        find_tag(aux.as_bytes(), *b"CG")
    else {
        return None;
    };
    let n_ops = ops.len() / 4;
    (n_ops >= cigar.len() / 4 && n_ops < 1 << 29).then_some(ops)
}

/// The query length of packed CIGAR operations (`bam_cigar2qlen`).
fn query_len(cigar: &[u8]) -> u64 {
    cigar
        .as_chunks::<4>()
        .0
        .iter()
        .map(|op| CigarOp::from_bam_u32(u32::from_le_bytes(*op)))
        .filter(|op| op.consumes_query())
        .fold(0u64, |n, op| n.saturating_add(u64::from(op.len())))
}

// r[impl bam.raw_record.cigar_seq]
/// htslib's read-time sanity check (`bam_read1` and `bam_tag2cigar`): the
/// long CIGAR of a record carrying one, else the stored CIGAR of a mapped
/// record, must consume exactly the sequence.
fn check_cigar_against_seq(
    payload: &[u8],
    h: &record::ParsedHeader,
    placed: bool,
) -> Result<(), DecodeError> {
    let seq_len = h.seq_len;
    let stored = payload.get(h.var_start..h.cigar_end).unwrap_or_default();
    let aux = Aux::new(payload.get(h.qual_end..).unwrap_or_default());
    let checked = match long_cigar(stored, seq_len as usize, placed, aux) {
        Some(cg) => cg,
        None if h.n_cigar_ops != 0 && !h.flags.is_unmapped() => stored,
        None => return Ok(()),
    };
    let query_len = query_len(checked);
    if seq_len == 0 || query_len == u64::from(seq_len) {
        return Ok(());
    }
    Err(DecodeError::CigarQueryLenMismatch {
        cigar_query_len: u32::try_from(query_len).unwrap_or(u32::MAX),
        seq_len,
    })
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::arithmetic_side_effects, reason = "test code")]
mod tests {
    use super::*;

    /// A record payload on target 0 at `pos` (or unplaced), named `r`.
    pub(super) fn payload(
        placed: bool,
        flags: u16,
        cigar: &[u32],
        seq_codes: &[u8],
        aux: &[u8],
    ) -> Vec<u8> {
        let refid: i32 = if placed { 0 } else { -1 };
        let mut p = Vec::new();
        p.extend_from_slice(&refid.to_le_bytes());
        p.extend_from_slice(&refid.to_le_bytes()); // pos 0 or -1
        p.push(2); // l_read_name
        p.push(60);
        p.extend_from_slice(&4680u16.to_le_bytes());
        p.extend_from_slice(&u16::try_from(cigar.len()).unwrap().to_le_bytes());
        p.extend_from_slice(&flags.to_le_bytes());
        p.extend_from_slice(&u32::try_from(seq_codes.len()).unwrap().to_le_bytes());
        p.extend_from_slice(&(-1i32).to_le_bytes());
        p.extend_from_slice(&(-1i32).to_le_bytes());
        p.extend_from_slice(&0i32.to_le_bytes());
        p.extend_from_slice(b"r\0");
        for op in cigar {
            p.extend_from_slice(&op.to_le_bytes());
        }
        for pair in seq_codes.chunks(2) {
            p.push(pair[0] << 4 | pair.get(1).copied().unwrap_or(0));
        }
        p.extend(std::iter::repeat_n(30u8, seq_codes.len()));
        p.extend_from_slice(aux);
        p
    }

    fn header() -> BamHeader {
        BamHeader::from_sam_text("@SQ\tSN:chr1\tLN:1000\n").unwrap()
    }

    #[test]
    fn decode_seq_into_keeps_every_iupac_code() {
        let codes: Vec<u8> = (0..16).chain([3]).collect();
        let record =
            RawBamRecord::from_payload(payload(false, 4, &[], &codes, &[]), &header()).unwrap();
        let mut out = b">".to_vec();
        record.decode_seq_into(&mut out);
        assert_eq!(out, b">=ACMGRSVTWYHKDBNM");
    }

    fn op(len: u32, code: u32) -> u32 {
        len << 4 | code
    }

    fn cg(subtype: u8, ops: &[u32]) -> Vec<u8> {
        let mut aux = vec![b'C', b'G', b'B', subtype];
        aux.extend_from_slice(&u32::try_from(ops.len()).unwrap().to_le_bytes());
        for op in ops {
            aux.extend_from_slice(&op.to_le_bytes());
        }
        aux
    }

    fn record(
        placed: bool,
        flags: u16,
        cigar: &[u32],
        aux: &[u8],
    ) -> Result<RawBamRecord, DecodeError> {
        RawBamRecord::from_payload(payload(placed, flags, cigar, &[1, 2, 4, 8], aux), &header())
    }

    // r[verify bam.raw_record.long_cigar]
    #[test]
    fn long_cigar_follows_bam_tag2cigar() {
        let placeholder = [op(4, 4), op(3, 3)];
        let real = [op(2, 0), op(1, 1), op(1, 4)];
        // B:I and B:i both count.
        for subtype in *b"Ii" {
            let r = record(true, 0, &placeholder, &cg(subtype, &real)).unwrap();
            assert_eq!(r.long_cigar(), Some(&cg(subtype, &real)[8..]));
            assert!(r.has_cg_placeholder());
        }
        // Unplaced, other subtypes, too few ops, first op not exactly
        // `<seq_len>S`: no long CIGAR.
        let not_long = [
            record(false, 4, &placeholder, &cg(b'I', &real)),
            record(true, 0, &placeholder, &cg(b'S', &[1, 2, 3])),
            record(true, 0, &placeholder, &cg(b'I', &[op(4, 0)])),
            record(true, 0, &[op(4, 4)], &cg(b'I', &[])),
            record(true, 0, &[op(4, 0), op(3, 3)], &cg(b'I', &real)),
        ];
        for r in not_long {
            assert!(!r.unwrap().has_cg_placeholder());
        }
        // Only the first CG tag is looked at.
        let mut aux = cg(b'S', &[1, 2]);
        aux.extend(cg(b'I', &real));
        assert!(!record(true, 0, &placeholder, &aux).unwrap().has_cg_placeholder());
    }

    // r[verify bam.raw_record.cigar_seq]
    #[test]
    fn long_cigar_must_cover_the_sequence_even_when_unmapped() {
        let placeholder = [op(4, 4), op(3, 3)];
        let short = [op(2, 0), op(1, 2)];
        for flags in [0, 4] {
            let err = record(true, flags, &placeholder, &cg(b'I', &short)).unwrap_err();
            assert!(matches!(
                err,
                DecodeError::CigarQueryLenMismatch { cigar_query_len: 2, seq_len: 4 }
            ));
        }
        // The stored CIGAR of an unmapped record is not checked.
        assert!(record(true, 4, &[op(1, 0)], &[]).is_ok());
        assert!(record(true, 0, &[op(1, 0)], &[]).is_err());
    }
}
