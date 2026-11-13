//! [`RawBamRecord`]: one BAM record kept in its wire form, validated once.
//!
//! The streaming readers (`r[bam.stream]`, `r[sam.stream]`) yield these, and
//! [`BamWriter::write_raw`](super::BamWriter::write_raw) writes one back byte
//! for byte. Nothing is decoded until an accessor asks for it, and nothing is
//! lost: unplaced records, IUPAC bases, every aux tag and its integer width,
//! and the stored `bin` all survive a read and a write unchanged.

use super::{
    AuxData, AuxDataError, BamHeader, CigarOp, DecodeError, MateFields,
    aux::{Aux, AuxLayoutError, AuxValue, find_tag, validate_aux},
    cigar::CIGAR_S,
    record::{self, read2},
    seq,
};
use crate::reader::Tid;
use core::range::RangeInclusive;
use seqair_types::{BamFlags, Base, BaseQuality, IupacBase, MappingQuality, Pos0, QPos};
use std::fmt;
use thiserror::Error;

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
/// Two things are *not* checked, because htslib does not check them on read:
///
/// - The aux block. htslib passes a malformed one through untouched and fails
///   only when formatting it. [`Aux::iter_checked`] and [`validate_aux`]
///   report the problem.
/// - The stored `bin`, kept as the file has it. htslib recomputes it on read,
///   so `samtools view -b` may rewrite it.
///
/// One check is stricter than htslib: a read name without its NUL terminator
/// is rejected, where htslib repairs it.
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

/// Why an in-place edit of a [`RawBamRecord`] was refused. The record is
/// unchanged.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum RawRecordEditError {
    #[error("a sequence of {got} bases cannot replace one of {expected}")]
    SeqLenMismatch { expected: usize, got: usize },

    #[error("base {qpos} is past the end of a sequence of {len} bases")]
    BaseOutOfRange { qpos: QPos, len: usize },

    #[error("the aux edit would change the record's long CIGAR")]
    LongCigarChanged,

    #[error(
        "the edited record would be {len} bytes, over the {} byte limit",
        record::MAX_RECORD_SIZE
    )]
    RecordTooLarge { len: usize },

    /// The record's offsets do not describe its payload; `from_payload`
    /// rules this out, so it is a bug in seqair, not the caller's input.
    #[error("the record's field offsets are inconsistent with its payload (a seqair bug)")]
    InconsistentLayout,

    #[error("FLAG would contradict the record's CIGAR and sequence")]
    Flags(#[from] DecodeError),

    #[error("aux edit failed")]
    Aux(#[from] AuxDataError),

    #[error("the record's aux block is malformed")]
    MalformedAux(#[from] AuxLayoutError),
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
            .field("qname", &self.qname().escape_ascii().to_string())
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
        if self.payload.len() > record::MAX_RECORD_SIZE {
            return Err(DecodeError::RecordTooLarge { len: self.payload.len() });
        }
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
            self.stored_cigar_bytes(),
            self.seq_len(),
            self.tid.is_some() && self.pos.is_some(),
            self.aux(),
        )
    }

    /// Become [`Self::default`], keeping the payload's allocation.
    fn reset(&mut self) {
        self.payload.clear();
        self.payload.extend_from_slice(&EMPTY_RECORD);
        self.tid = None;
        self.pos = None;
        self.mate = MateFields::default();
        self.cigar_start = EMPTY_RECORD.len();
        self.cigar_end = EMPTY_RECORD.len();
        self.seq_end = EMPTY_RECORD.len();
        self.qual_end = EMPTY_RECORD.len();
    }

    /// The record's bytes as BAM stores them, without `block_size`.
    pub fn as_payload(&self) -> &[u8] {
        &self.payload
    }

    /// Give up the payload bytes.
    pub fn into_payload(self) -> Vec<u8> {
        self.payload
    }

    // r[impl bam.raw_record.take_payload]
    /// Move the payload into `out`, whose old contents are dropped, and
    /// become [`Self::default`] in `out`'s old allocation: a reader refilling
    /// this record reuses it, so handing records on costs neither a copy nor
    /// an allocation.
    pub fn take_payload_into(&mut self, out: &mut Vec<u8>) {
        std::mem::swap(&mut self.payload, out);
        self.reset();
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

    /// The number of CIGAR operations: of the real CIGAR, for a record that
    /// carries a long one in its `CG` tag (`r[bam.raw_record.long_cigar]`).
    pub fn cigar_len(&self) -> usize {
        self.cigar_packed().len() / 4
    }

    // r[impl bam.raw_record.cigar_access]
    /// The CIGAR operations — for a record with a long CIGAR, the real ones
    /// from its `CG` tag, as htslib's reader returns them. See
    /// [`Self::stored_cigar_bytes`] for the placeholder the record stores.
    pub fn cigar(&self) -> impl ExactSizeIterator<Item = CigarOp> + '_ {
        CigarOp::iter_bam_bytes(self.cigar_packed())
    }

    /// The real CIGAR in BAM's packed form: the long one, else the stored one.
    pub(crate) fn cigar_packed(&self) -> &[u8] {
        self.long_cigar().unwrap_or(self.stored_cigar_bytes())
    }

    /// The CIGAR field as stored, in BAM's packed form (little-endian `u32`
    /// per op) — for a long CIGAR, the `<seq_len>S<ref_len>N` placeholder.
    pub fn stored_cigar_bytes(&self) -> &[u8] {
        self.payload.get(self.cigar_start..self.cigar_end).unwrap_or_default()
    }

    /// The last reference position the alignment covers, by its CIGAR
    /// (`r[bam.record.end_pos]`); `None` for an unplaced record.
    pub fn end_pos(&self) -> Option<Pos0> {
        record::compute_end_pos_from_raw(&self.payload)
    }

    // r[impl interval.overlap_test]
    /// Whether the record overlaps `span` as region queries test it: an
    /// unmapped one sits at its POS whatever its CIGAR
    /// (`r[record_store.end_pos_htslib]`); an unplaced one overlaps nothing.
    pub(crate) fn overlaps(&self, span: RangeInclusive<Pos0>) -> bool {
        let Some(pos) = self.pos else { return false };
        let last = if self.flags().is_unmapped() { pos } else { self.end_pos().unwrap_or(pos) };
        pos <= span.last && last >= span.start
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
        self.iupac_bases().map(IupacBase::to_base).collect()
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

    // r[impl bam.raw_record.edit_qual]
    /// Base qualities to edit in place, or `None` exactly when [`Self::qual`]
    /// is. Writing [`BaseQuality::UNAVAILABLE`] to the first base makes the
    /// record one without qualities.
    pub fn qual_mut(&mut self) -> Option<&mut [BaseQuality]> {
        let qual = self.payload.get_mut(self.seq_end..self.qual_end).unwrap_or_default();
        match qual.first() {
            None | Some(&0xFF) => None,
            Some(_) => Some(bytemuck::cast_slice_mut(qual)),
        }
    }

    /// The aux tags, as stored. Not validated: use [`Aux::iter_checked`] to
    /// see a malformed tag as an error rather than the end of the tags.
    pub fn aux(&self) -> Aux<'_> {
        Aux::new(self.payload.get(self.qual_end..).unwrap_or_default())
    }

    /// Set FLAG in place.
    ///
    /// # Errors
    /// [`RawRecordEditError::Flags`] wrapping [`DecodeError::CigarQueryLenMismatch`], leaving FLAG as it was, if the
    /// new flags make the record one whose CIGAR must match its sequence and
    /// it does not (`r[bam.raw_record.cigar_seq]`) — clearing `0x4` on an
    /// unmapped record with a stale CIGAR, say. htslib rejects such a record.
    pub fn set_flags(&mut self, flags: BamFlags) -> Result<(), RawRecordEditError> {
        let old = self.flags();
        self.write_flags(flags);
        let checked = record::parse_header(&self.payload).and_then(|h| {
            check_cigar_against_seq(&self.payload, &h, self.tid.is_some() && self.pos.is_some())
        });
        if checked.is_err() {
            self.write_flags(old);
        }
        Ok(checked?)
    }

    fn write_flags(&mut self, flags: BamFlags) {
        if let Some(dst) = self.payload.get_mut(14..16) {
            dst.copy_from_slice(&flags.raw().to_le_bytes());
        }
    }

    // r[impl bam.raw_record.iupac_bases]
    /// The sequence, one [`IupacBase`] per base, every 4-bit code kept.
    pub fn iupac_bases(&self) -> impl ExactSizeIterator<Item = IupacBase> + '_ {
        let packed = self.packed_seq();
        (0..self.seq_len()).map(move |i| {
            let byte = packed.get(i / 2).copied().unwrap_or_default();
            IupacBase::from_nibble(if i.is_multiple_of(2) { byte >> 4 } else { byte })
        })
    }

    // r[impl bam.raw_record.base_at]
    /// The base at `qpos`, its 4-bit code kept; `None` past the end of the
    /// sequence.
    pub fn base_at(&self, qpos: QPos) -> Option<IupacBase> {
        let i = qpos.as_usize();
        if i >= self.seq_len() {
            return None;
        }
        let byte = *self.packed_seq().get(i / 2)?;
        Some(IupacBase::from_nibble(if i.is_multiple_of(2) { byte >> 4 } else { byte }))
    }

    // r[impl seq.iupac_bulk+3]
    /// Append the sequence to `out`, one [`IupacBase`] per base, every 4-bit
    /// code kept: [`Self::iupac_bases`] into a caller's buffer.
    pub fn iupac_bases_into(&self, out: &mut Vec<IupacBase>) {
        let start = out.len();
        out.resize(start.saturating_add(self.seq_len()), IupacBase::N);
        if let Some(dst) = out.get_mut(start..) {
            seq::decode_iupac_bases_into(self.packed_seq(), dst);
        }
    }

    // r[impl bam.raw_record.edit_seq]
    /// Replace the sequence with `seq`, in place.
    ///
    /// # Errors
    /// [`RawRecordEditError::SeqLenMismatch`] if `seq` is not as long as the
    /// current sequence.
    pub fn set_seq(&mut self, seq: &[IupacBase]) -> Result<(), RawRecordEditError> {
        let expected = self.seq_len();
        if seq.len() != expected {
            return Err(RawRecordEditError::SeqLenMismatch { expected, got: seq.len() });
        }
        let Some(dst) = self.payload.get_mut(self.cigar_end..self.seq_end) else {
            debug_assert!(false, "offsets within the payload");
            return Err(RawRecordEditError::InconsistentLayout);
        };
        // The unused low nibble of an odd-length sequence is not ours to touch.
        let pad = if expected % 2 == 1 { dst.last().map(|b| b & 0x0F) } else { None };
        // r[impl seq.iupac_bulk+3]
        seq::encode_iupac_bases_to(seq, dst);
        if let (Some(pad), Some(last)) = (pad, dst.last_mut()) {
            *last = (*last & 0xF0) | pad;
        }
        Ok(())
    }

    // r[impl bam.raw_record.edit_seq]
    /// Replace the base at `qpos` with `base`, in place.
    ///
    /// # Errors
    /// [`RawRecordEditError::BaseOutOfRange`] past the end of the sequence.
    pub fn set_base(&mut self, qpos: QPos, base: IupacBase) -> Result<(), RawRecordEditError> {
        let len = self.seq_len();
        let i = qpos.as_usize();
        if i >= len {
            return Err(RawRecordEditError::BaseOutOfRange { qpos, len });
        }
        let start = self.cigar_end;
        let Some(byte) = self.payload.get_mut(start.saturating_add(i / 2)) else {
            debug_assert!(false, "offsets within the payload");
            return Err(RawRecordEditError::InconsistentLayout);
        };
        let code = base.code();
        *byte =
            if i.is_multiple_of(2) { (code << 4) | (*byte & 0x0F) } else { (*byte & 0xF0) | code };
        Ok(())
    }

    // r[impl bam.raw_record.edit_aux]
    /// Edit the aux block: `edit` gets the tags as an [`AuxData`] and the
    /// record takes what it leaves. Tags it does not touch keep their bytes.
    /// Costs one copy of the aux block.
    ///
    /// # Errors
    /// [`RawRecordEditError::MalformedAux`] if the existing tags are not
    /// well-formed, `edit`'s own error, [`RawRecordEditError::RecordTooLarge`], or
    /// [`RawRecordEditError::LongCigarChanged`] if the real CIGAR would
    /// differ (a `CG` tag removed or rewritten).
    pub fn edit_aux(
        &mut self,
        edit: impl FnOnce(&mut AuxData) -> Result<(), AuxDataError>,
    ) -> Result<(), RawRecordEditError> {
        let current = self.payload.get(self.qual_end..).unwrap_or_default();
        validate_aux(current)?;
        let mut aux = AuxData::from_bytes(current.to_vec());
        edit(&mut aux)?;
        let len = self.qual_end.saturating_add(aux.as_bytes().len());
        if len > record::MAX_RECORD_SIZE {
            return Err(RawRecordEditError::RecordTooLarge { len });
        }
        let after = long_cigar(
            self.stored_cigar_bytes(),
            self.seq_len(),
            self.tid.is_some() && self.pos.is_some(),
            Aux::new(aux.as_bytes()),
        );
        if self.long_cigar() != after {
            return Err(RawRecordEditError::LongCigarChanged);
        }
        self.payload.truncate(self.qual_end);
        self.payload.extend_from_slice(aux.as_bytes());
        Ok(())
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
    let first = CigarOp::iter_bam_bytes(cigar).next()?.to_bam_u32();
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
    let query_len = super::cigar::query_len(CigarOp::iter_bam_bytes(checked));
    if seq_len == 0 || query_len == u64::from(seq_len) {
        return Ok(());
    }
    Err(DecodeError::CigarQueryLenMismatch { cigar_query_len: query_len, seq_len })
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

    // r[verify bam.raw_record.take_payload]
    #[test]
    fn take_payload_into_hands_over_the_bytes_and_keeps_the_allocation() {
        let bytes = payload(true, 0, &[4 << 4], &[1, 2, 4, 8], b"");
        let mut record = RawBamRecord::from_payload(bytes.clone(), &header()).unwrap();
        let mut out = Vec::with_capacity(4096);
        out.extend_from_slice(b"old");
        let spare = out.as_ptr();
        record.take_payload_into(&mut out);
        assert_eq!(out, bytes);
        assert_eq!(record, RawBamRecord::default());
        assert_eq!(record.as_payload().as_ptr(), spare);
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

    // r[verify bam.raw_record.edit_aux]
    #[test]
    fn edit_aux_refuses_a_malformed_block_without_panicking() {
        let truncated_ml: Vec<u8> = [b"MLBC".as_slice(), &100u32.to_le_bytes(), &[7]].concat();
        let shapes: [&[u8]; 3] = [b"XMZabc", &truncated_ml, b"NMi\x01"];
        for aux in shapes {
            let mut record =
                RawBamRecord::from_payload(payload(false, 4, &[], &[1, 2], aux), &header())
                    .unwrap();
            let before = record.as_payload().to_vec();
            let err = record
                .edit_aux(|a| {
                    a.remove(*b"XM");
                    a.remove(*b"ML");
                    a.remove(*b"NM");
                    Ok(())
                })
                .unwrap_err();
            assert!(matches!(err, RawRecordEditError::MalformedAux(_)), "{aux:?}: {err:?}");
            assert_eq!(record.as_payload(), before);
        }
    }

    // r[verify bam.raw_record.edit_aux]
    #[test]
    fn edit_aux_stops_at_the_record_size_limit() {
        let mut record = RawBamRecord::default();
        // 34 fixed bytes + tag name, type and NUL = 38 bytes besides the text.
        let fits = record::MAX_RECORD_SIZE - 38;
        record
            .edit_aux(|a| {
                a.set_string(*b"XM", &vec![b'a'; fits]).unwrap();
                Ok(())
            })
            .unwrap();
        assert_eq!(record.as_payload().len(), record::MAX_RECORD_SIZE);
        let before = record.as_payload().to_vec();
        let err = record
            .edit_aux(|a| {
                a.set_string(*b"XM", &vec![b'a'; fits + 1]).unwrap();
                Ok(())
            })
            .unwrap_err();
        assert!(
            matches!(err, RawRecordEditError::RecordTooLarge { len } if len == record::MAX_RECORD_SIZE + 1)
        );
        assert_eq!(record.as_payload(), before);
    }

    // r[verify bam.raw_record.base_at]
    // r[verify bam.raw_record.edit_seq]
    /// `base_at` reads the code stored at each position, `None` past the
    /// end, and gives back what `set_base` wrote.
    #[hegel::test]
    fn base_at_reads_each_stored_code(tc: hegel::TestCase) {
        use hegel::generators as gs;
        let len = tc.draw(gs::integers::<usize>().max_value(21));
        let codes: Vec<u8> =
            (0..len).map(|_| tc.draw_silent(gs::integers::<u8>().max_value(15))).collect();
        let mut record =
            RawBamRecord::from_payload(payload(false, 4, &[], &codes, &[]), &header()).unwrap();
        for i in 0..len + 3 {
            let qpos = QPos::new(u32::try_from(i).unwrap());
            assert_eq!(record.base_at(qpos), codes.get(i).map(|&c| IupacBase::from_nibble(c)));
        }
        if len > 0 {
            let at = QPos::new(
                u32::try_from(tc.draw(gs::integers::<usize>().max_value(len - 1))).unwrap(),
            );
            let base = IupacBase::from_nibble(tc.draw(gs::integers::<u8>().max_value(15)));
            record.set_base(at, base).unwrap();
            assert_eq!(record.base_at(at), Some(base));
        }
    }

    // r[verify bam.raw_record.iupac_bases]
    // r[verify seq.iupac_bulk+3]
    #[test]
    fn iupac_bases_agree_with_the_ascii_decode() {
        let codes: Vec<u8> = (0..16).chain([3]).collect();
        let record =
            RawBamRecord::from_payload(payload(false, 4, &[], &codes, &[]), &header()).unwrap();
        let bases: Vec<IupacBase> = record.iupac_bases().collect();
        assert_eq!(bases.len(), codes.len());
        let mut appended = vec![IupacBase::R];
        record.iupac_bases_into(&mut appended);
        assert_eq!(appended[0], IupacBase::R);
        assert_eq!(appended[1..], bases);
        assert_eq!(bases.iter().map(|b| b.code()).collect::<Vec<_>>(), codes);
        let mut ascii = Vec::new();
        record.decode_seq_into(&mut ascii);
        assert_eq!(IupacBase::to_ascii_slice(&bases), ascii);
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

    // r[verify cigar.aligned_pairs.raw_record]
    /// Walking the packed CIGAR yields what walking the decoded ops does, a
    /// long CIGAR's real ops included; an unplaced record with a CIGAR is
    /// refused.
    #[hegel::test]
    fn aligned_pairs_walk_the_packed_cigar_as_the_decoded_one(tc: hegel::TestCase) {
        use super::super::aligned_pairs::{AlignedPairs, AlignedPairsError};
        use hegel::generators as gs;
        let n = tc.draw(gs::integers::<usize>().max_value(6));
        let real: Vec<u32> = (0..n)
            .map(|_| {
                let code = tc.draw_silent(gs::sampled_from(vec![0u32, 1, 2, 3, 4, 7, 8]));
                op(tc.draw_silent(gs::integers::<u32>().max_value(4)), code)
            })
            .collect();
        let decode = |ops: &[u32]| -> Vec<CigarOp> {
            ops.iter().map(|&op| CigarOp::from_bam_u32(op)).collect()
        };
        let real_ops = decode(&real);
        let qlen: u32 = real_ops.iter().filter(|o| o.consumes_query()).map(|o| o.len()).sum();
        let rlen: u32 = real_ops.iter().filter(|o| o.consumes_ref()).map(|o| o.len()).sum();
        let placed = tc.draw(gs::booleans());
        let long = placed && qlen > 0 && tc.draw(gs::booleans());
        let (stored, aux) = if long {
            (vec![op(qlen, 4), op(rlen, 3)], cg(b'I', &real))
        } else {
            (real.clone(), vec![])
        };
        let codes = vec![1u8; qlen as usize];
        let flags = if placed { 0 } else { 4 };
        let record =
            RawBamRecord::from_payload(payload(placed, flags, &stored, &codes, &aux), &header())
                .unwrap();
        // htslib takes a CG tag as the real CIGAR only under its own rules
        // (`r[bam.raw_record.long_cigar]`), which a short real CIGAR fails.
        let walked = decode(if record.has_cg_placeholder() { &real } else { &stored });

        match record.aligned_pairs() {
            Ok(pairs) => {
                assert!(placed || walked.is_empty());
                let expected: Vec<_> = AlignedPairs::new(Pos0::ZERO, &walked).full().collect();
                assert_eq!(pairs.full().collect::<Vec<_>>(), expected);
            }
            Err(AlignedPairsError::UnmappedWithCigar { cigar_ops }) => {
                assert!(!placed && !walked.is_empty());
                assert_eq!(cigar_ops, walked.len());
            }
            Err(other) => panic!("unexpected error {other:?}"),
        }
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
