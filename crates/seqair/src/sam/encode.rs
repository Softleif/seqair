//! SAM text → BAM records, byte for byte as htslib encodes them.
//!
//! [`SamRecordEncoder`] turns one alignment line into a [`RawBamRecord`]
//! (`r[sam.encode]`); [`parse_cigar`] is its CIGAR parser, usable on its own
//! (for CIGARs in `SA`/`XA` tags, say). The optional-field parser here is the
//! one the indexed reader uses too (`r[sam.aux.parse]`).

use crate::bam::{
    BamHeader, CigarOp, DecodeError, MateFields, RawBamRecord,
    aux::tag_end,
    cigar::CigarOpType,
    raw_record::Refill,
    record::{FixedHeaderFields, MAX_RECORD_SIZE, encode_fixed_header},
    seq,
};
use crate::io::text_scan;
use crate::reader::Tid;
use seqair_types::{IupacBase, MappingQuality, Pos0};
use std::fmt;
use std::num::NonZeroI32;
use std::sync::atomic::{AtomicBool, Ordering};
use tracing::warn;

/// Shows bytes for an error message: printable ASCII as is, anything else as
/// `\xNN`, at most 32 input bytes.
struct Escaped<'a>(&'a [u8]);

impl fmt::Display for Escaped<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        const MAX: usize = 32;
        for &b in self.0.iter().take(MAX) {
            if b.is_ascii_graphic() || b == b' ' {
                fmt::Write::write_char(f, char::from(b))?;
            } else {
                write!(f, "\\x{b:02x}")?;
            }
        }
        if self.0.len() > MAX {
            f.write_str("...")?;
        }
        Ok(())
    }
}

/// Why a SAM alignment line could not be turned into a BAM record.
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum SamRecordError {
    #[error("empty line")]
    EmptyLine,

    #[error("expected >= 11 TAB-separated fields, got {found}")]
    TooFewFields { found: usize },

    /// A NUL byte, which neither SAM text nor a BAM string can hold.
    #[error("NUL byte at offset {offset}")]
    NulByte { offset: usize },

    #[error("QNAME is {len} bytes; BAM allows at most 254")]
    QnameTooLong { len: usize },

    /// Not `[0-9]+`, above 65535, or a leading zero (which htslib reads as
    /// octal or hex).
    #[error("invalid FLAG field: {}", Escaped(value))]
    InvalidFlag { value: Box<[u8]> },

    #[error("RNAME is not valid UTF-8: {}", Escaped(value))]
    InvalidRname { value: Box<[u8]> },

    /// RNAME names no target of the header.
    #[error("RNAME {} is not a reference of the header", Escaped(name))]
    UnknownReference { name: Box<[u8]> },

    /// Not `[0-9]+`, or above 2^31.
    #[error("invalid POS field: {}", Escaped(value))]
    InvalidPos { value: Box<[u8]> },

    #[error("invalid MAPQ field: {}", Escaped(value))]
    InvalidMapq { value: Box<[u8]> },

    /// PNEXT is not a valid 1-based position: not a number, negative, or out
    /// of range.
    #[error("invalid PNEXT field: {}", Escaped(value))]
    InvalidPnext { value: Box<[u8]> },

    /// TLEN is not an integer that fits in an `i32`.
    #[error("invalid TLEN field: {}", Escaped(value))]
    InvalidTlen { value: Box<[u8]> },

    #[error("empty CIGAR")]
    EmptyCigar,

    /// A CIGAR operation's length is missing, not a number, or 2^28 or more.
    #[error("invalid CIGAR operation length: {}", Escaped(value))]
    InvalidCigarLength { value: Box<[u8]> },

    #[error("unknown CIGAR operation: {}", Escaped(&[*op]))]
    UnknownCigarOp { op: u8 },

    /// The CIGAR ends in a length without an operation.
    #[error("CIGAR ends in a length without an operation: {}", Escaped(value))]
    MissingCigarOp { value: Box<[u8]> },

    #[error("CIGAR query length {cigar_query_len} does not match SEQ length {seq_len}")]
    CigarSeqLenMismatch { cigar_query_len: u64, seq_len: usize },

    /// A CIGAR of more than 65535 operations is stored as a placeholder whose
    /// `N` operation holds the reference length, which must be below 2^28.
    #[error(
        "CIGAR of more than 65535 operations covers {ref_len} reference bases, BAM allows < 2^28"
    )]
    CigarRefLenTooLong { ref_len: u64 },

    #[error("invalid SEQ byte {} at offset {offset}", Escaped(&[*byte]))]
    InvalidSeq { byte: u8, offset: usize },

    #[error("QUAL has {qual_len} characters but SEQ has {seq_len}")]
    QualLenMismatch { qual_len: usize, seq_len: usize },

    #[error("invalid QUAL byte {} at offset {offset}", Escaped(&[*byte]))]
    InvalidQual { byte: u8, offset: usize },

    /// An optional field that is not `TG:T:VALUE`.
    #[error("malformed optional field: {}", Escaped(field))]
    InvalidAuxField { field: Box<[u8]> },

    #[error("optional field {} has unknown type {}", Escaped(tag), Escaped(&[*type_code]))]
    UnknownAuxType { tag: [u8; 2], type_code: u8 },

    #[error("invalid value for optional field {}: {}", Escaped(tag), Escaped(value))]
    InvalidAuxValue { tag: [u8; 2], value: Box<[u8]> },

    #[error(
        "optional field {} value {} does not fit any BAM integer type (i8/u8/i16/u16/i32/u32)",
        Escaped(tag),
        Escaped(value)
    )]
    AuxIntOutOfRange { tag: [u8; 2], value: Box<[u8]> },

    #[error("aux B-array has {len} elements, exceeding u32::MAX")]
    AuxArrayTooLarge { len: usize },

    #[error("optional field {} has B-array subtype {}, not one of cCsSiIf", Escaped(tag), Escaped(&[*subtype]))]
    InvalidAuxSubtype { tag: [u8; 2], subtype: u8 },

    #[error("record would be {size} bytes, over the {MAX_RECORD_SIZE}-byte limit")]
    RecordTooLarge { size: usize },

    #[error("SAM header is not valid UTF-8")]
    HeaderNotUtf8 { source: std::string::FromUtf8Error },

    /// The encoded bytes failed `RawBamRecord` validation. The encoder checks
    /// everything that validation does first, so this means a bug.
    #[error("encoded record failed validation")]
    InvalidRecord { source: DecodeError },
}

/// Reusable SAM-line → BAM-record encoder (`r[sam.encode]`). Keep one per
/// thread and call [`Self::encode`] for every line: its scratch buffers and
/// the output record's buffer are reused, so steady-state encoding does not
/// allocate.
///
/// ```
/// use seqair::bam::{BamHeader, RawBamRecord};
/// use seqair::sam::SamRecordEncoder;
///
/// let header = BamHeader::from_sam_text("@SQ\tSN:chr1\tLN:1000\n")?;
/// let mut encoder = SamRecordEncoder::new();
/// let mut record = RawBamRecord::default();
/// encoder.encode(b"r1\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0", &header, &mut record)?;
/// assert_eq!(record.qname(), b"r1");
/// assert_eq!(record.seq_len(), 4);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Default)]
pub struct SamRecordEncoder {
    cigar: Vec<CigarOp>,
    aux: Vec<u8>,
}

impl SamRecordEncoder {
    pub fn new() -> Self {
        Self::default()
    }

    // r[impl sam.encode]
    /// Encode one SAM alignment line (no line terminator) against `header`
    /// into `out`, reusing its buffer.
    ///
    /// # Errors
    /// The first problem with the line. `out` is then the default record.
    pub fn encode(
        &mut self,
        line: &[u8],
        header: &BamHeader,
        out: &mut RawBamRecord,
    ) -> Result<(), SamRecordError> {
        out.refill(header, |payload| self.encode_payload(line, header, payload)).map_err(
            |e| match e {
                Refill::Fill(e) => e,
                Refill::Invalid(source) => SamRecordError::InvalidRecord { source },
            },
        )
    }

    // r[impl sam.encode.deviations]
    fn encode_payload(
        &mut self,
        line: &[u8],
        header: &BamHeader,
        payload: &mut Vec<u8>,
    ) -> Result<(), SamRecordError> {
        // r[impl sam.encode.fields]
        if line.is_empty() {
            return Err(SamRecordError::EmptyLine);
        }
        if let Some(offset) = text_scan::find_byte(line, 0) {
            return Err(SamRecordError::NulByte { offset });
        }
        let (fields, n_fields) = text_scan::split_tabs(line);
        let [qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, rest] = fields;
        if n_fields < 11 {
            return Err(SamRecordError::TooFewFields { found: n_fields });
        }
        if qname.len() > 254 {
            return Err(SamRecordError::QnameTooLong { len: qname.len() });
        }

        // r[impl sam.encode.flag]
        let mut flags = parse_flag(flag)?;

        // r[impl sam.encode.references+2]
        // r[impl sam.encode.fixups]
        // The order is htslib's: each fix-up sees the ones before it.
        let mut tid = match rname {
            b"*" => None,
            name => Some(
                resolve(header, name)
                    .ok_or_else(|| SamRecordError::UnknownReference { name: name.into() })?,
            ),
        };
        let pos = parse_position(pos, |value| SamRecordError::InvalidPos { value })?;
        if pos.is_none() {
            tid = None;
        }
        if tid.is_none() {
            flags |= UNMAPPED;
        }
        let mapq = parse_digits(mapq)
            .and_then(|v| u8::try_from(v).ok())
            .ok_or_else(|| SamRecordError::InvalidMapq { value: mapq.into() })?;
        if cigar == b"*" {
            self.cigar.clear();
            flags |= UNMAPPED;
        } else {
            parse_cigar(cigar, &mut self.cigar)?;
        }
        let mut mate_tid = match rnext {
            b"=" => tid,
            b"*" => None,
            name => resolve(header, name).or_else(|| {
                warn_unknown_mate_reference(name);
                None
            }),
        };
        let mate_pos = parse_position(pnext, |value| SamRecordError::InvalidPnext { value })?;
        if mate_pos.is_none() {
            mate_tid = None;
        }
        let template_len = parse_signed(tlen)
            .ok()
            .and_then(|(negative, magnitude)| signed_value(negative, magnitude))
            .and_then(|v| i32::try_from(v).ok())
            .ok_or_else(|| SamRecordError::InvalidTlen { value: tlen.into() })?;

        // r[impl sam.encode.seq]
        let seq: &[u8] = if seq == b"*" { b"" } else { seq };
        if let Some(offset) = seq.iter().position(|b| !(b'!'..=b'~').contains(b)) {
            let byte = seq.get(offset).copied().unwrap_or_default();
            return Err(SamRecordError::InvalidSeq { byte, offset });
        }
        // r[impl sam.encode.qual]
        let qual: Option<&[u8]> = if qual == b"*" { None } else { Some(qual) };
        if let Some(qual) = qual {
            if qual.len() != seq.len() {
                return Err(SamRecordError::QualLenMismatch {
                    qual_len: qual.len(),
                    seq_len: seq.len(),
                });
            }
            if let Some(offset) = qual.iter().position(|b| !(33..=160).contains(b)) {
                let byte = qual.get(offset).copied().unwrap_or_default();
                return Err(SamRecordError::InvalidQual { byte, offset });
            }
        }

        self.aux.clear();
        if n_fields > 11 {
            parse_aux_fields(rest, &mut self.aux)?;
        }

        // r[impl sam.encode.cg_tag]
        if tid.is_some() && pos.is_some() {
            self.take_cg_tag(seq.len());
        }

        // r[impl sam.encode.cigar_seq]
        if !seq.is_empty() && !self.cigar.is_empty() {
            let cigar_query_len = query_len(&self.cigar);
            if cigar_query_len != seq.len() as u64 {
                return Err(SamRecordError::CigarSeqLenMismatch {
                    cigar_query_len,
                    seq_len: seq.len(),
                });
            }
        }

        encode_alignment(
            payload,
            &AlignmentFields {
                tid,
                pos,
                flags,
                mapq: MappingQuality::from_byte(mapq),
                mate_fields: MateFields {
                    tid: mate_tid,
                    pos: mate_pos,
                    template_len: NonZeroI32::new(template_len),
                },
                qname,
                cigar: &self.cigar,
                seq: Seq::Ascii(seq),
                qual: match qual {
                    None => Qualities::Missing,
                    Some(qual) => Qualities::Phred33(qual),
                },
                aux: &self.aux,
            },
        )
        .map_err(|e| match e {
            RecordEncodeError::QnameTooLong { len } => SamRecordError::QnameTooLong { len },
            RecordEncodeError::CigarRefLenTooLong { ref_len } => {
                SamRecordError::CigarRefLenTooLong { ref_len }
            }
            RecordEncodeError::RecordTooLarge { size } => SamRecordError::RecordTooLarge { size },
        })
    }

    /// htslib's `bam_tag2cigar` on a placed record: a first CIGAR operation
    /// `<l_seq>S` and a first `CG` tag of type `B:I`/`B:i` holding at least as
    /// many (and fewer than 2^29) operations make that tag the CIGAR.
    fn take_cg_tag(&mut self, l_seq: usize) {
        let Some(first) = self.cigar.first() else { return };
        if first.op_type() != CigarOpType::SoftClip || first.len() as usize != l_seq {
            return;
        }
        let mut at = 0;
        while at < self.aux.len() {
            // Our own parser wrote these bytes, so every tag is well formed.
            let Ok(end) = tag_end(&self.aux, at) else { return };
            let tag = self.aux.get(at..end).unwrap_or_default();
            if let Some((b"CG", value)) = tag.split_first_chunk::<2>() {
                let Some((&[b'B', sub, n0, n1, n2, n3], ops)) = value.split_first_chunk::<6>()
                else {
                    return;
                };
                let count = u32::from_le_bytes([n0, n1, n2, n3]);
                if !matches!(sub, b'I' | b'i')
                    || (count as usize) < self.cigar.len()
                    || count >= 1 << 29
                {
                    return;
                }
                self.cigar.clear();
                CigarOp::extend_from_bam_bytes(&mut self.cigar, ops);
                self.aux.drain(at..end);
                return;
            }
            at = end;
        }
    }
}

/// Why [`encode_alignment`] could not encode a record: a field BAM cannot
/// hold.
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum RecordEncodeError {
    #[error("QNAME is {len} bytes; BAM allows at most 254")]
    QnameTooLong { len: usize },

    #[error(
        "CIGAR of more than 65535 operations covers {ref_len} reference bases, BAM allows < 2^28"
    )]
    CigarRefLenTooLong { ref_len: u64 },

    #[error("record would be {size} bytes, over the {MAX_RECORD_SIZE}-byte limit")]
    RecordTooLarge { size: usize },
}

/// A record's qualities as [`encode_alignment`] takes them.
pub(crate) enum Qualities<'a> {
    /// None: every byte `0xFF`.
    Missing,
    /// SAM's QUAL text, Phred+33.
    Phred33(&'a [u8]),
    /// BAM's bytes, copied as they are.
    Raw(&'a [u8]),
}

/// A record's SEQ as [`encode_alignment`] takes it.
#[derive(Debug, Clone, Copy)]
pub(crate) enum Seq<'a> {
    /// ASCII, `=ACMGRSVTWYHKDBN` and lower case, as SAM writes it.
    Ascii(&'a [u8]),
    /// BAM's 4-bit codes.
    Iupac(&'a [IupacBase]),
}

impl Seq<'_> {
    fn len(self) -> usize {
        match self {
            Self::Ascii(s) => s.len(),
            Self::Iupac(s) => s.len(),
        }
    }

    fn encode_into(self, payload: &mut Vec<u8>) {
        match self {
            Self::Ascii(s) => seq::encode_seq_into(s, payload),
            Self::Iupac(s) => payload.extend(s.chunks(2).map(|pair| {
                let code = |i: usize| pair.get(i).map_or(0, |b| b.code());
                (code(0) << 4) | code(1)
            })),
        }
    }
}

/// One alignment's fields, as htslib's `bam_set1` takes them.
pub(crate) struct AlignmentFields<'a> {
    pub tid: Option<Tid>,
    pub pos: Option<Pos0>,
    pub flags: u16,
    pub mapq: MappingQuality,
    pub mate_fields: MateFields,
    pub qname: &'a [u8],
    pub cigar: &'a [CigarOp],
    pub seq: Seq<'a>,
    /// As many as SEQ, unless missing.
    pub qual: Qualities<'a>,
    pub aux: &'a [u8],
}

// r[impl sam.encode.bin]
// r[impl sam.encode.long_cigar]
// r[impl sam.encode.limits]
/// Append one record's BAM payload (its bytes after `block_size`) to
/// `payload`, as htslib's `bam_set1` and `bam_write1` lay it out: the bin
/// from the position and the CIGAR's reference length, and a CIGAR of more
/// than 65535 operations as the `<l_seq>S<rlen>N` placeholder with the real
/// one in a trailing `CG:B:I` tag.
pub(crate) fn encode_alignment(
    payload: &mut Vec<u8>,
    f: &AlignmentFields<'_>,
) -> Result<(), RecordEncodeError> {
    let ref_len = ref_len(f.cigar);
    let span = if f.flags & UNMAPPED != 0 || ref_len == 0 { 1 } else { ref_len };
    let beg = f.pos.map_or(-1, |p| i64::from(p.as_i32()));
    let bin = reg2bin(beg, beg.saturating_add_unsigned(span));

    let long = f.cigar.len() > usize::from(u16::MAX);
    if long && ref_len >= 1 << 28 {
        return Err(RecordEncodeError::CigarRefLenTooLong { ref_len });
    }
    let n_cigar_op = if long { 2 } else { u16::try_from(f.cigar.len()).unwrap_or(u16::MAX) };

    let cigar_bytes = f.cigar.len().saturating_mul(4);
    let size = [
        FIXED_LEN,
        f.qname.len().saturating_add(1),
        if long { 8 } else { cigar_bytes },
        f.seq.len().div_ceil(2),
        f.seq.len(),
        f.aux.len(),
        if long { cigar_bytes.saturating_add(8) } else { 0 },
    ]
    .into_iter()
    .fold(0usize, usize::saturating_add);
    if size > MAX_RECORD_SIZE {
        return Err(RecordEncodeError::RecordTooLarge { size });
    }
    // Both fit: `size` (which includes them) is at most 2 MiB.
    let l_seq =
        i32::try_from(f.seq.len()).map_err(|_| RecordEncodeError::RecordTooLarge { size })?;
    let l_read_name = u8::try_from(f.qname.len().saturating_add(1))
        .map_err(|_| RecordEncodeError::QnameTooLong { len: f.qname.len() })?;

    payload.reserve(size);
    encode_fixed_header(
        payload,
        &FixedHeaderFields {
            tid: f.tid,
            pos: f.pos,
            bin,
            mapq: f.mapq,
            l_read_name,
            flags: f.flags,
            n_cigar_op,
            l_seq,
            mate_fields: f.mate_fields,
        },
    );
    payload.extend_from_slice(f.qname);
    payload.push(0);
    if long {
        // `bam_write1`'s placeholder: `<l_seq>S<ref_len>N`. Both lengths
        // are below 2^28 (checked above; `l_seq` by the size limit).
        let ref_len = u32::try_from(ref_len).unwrap_or(u32::MAX);
        for op in [
            CigarOp::new(CigarOpType::SoftClip, l_seq.cast_unsigned()),
            CigarOp::new(CigarOpType::RefSkip, ref_len),
        ] {
            payload.extend_from_slice(&op.to_bam_u32().to_le_bytes());
        }
    } else {
        payload.extend_from_slice(CigarOp::ops_as_bytes(f.cigar));
    }
    f.seq.encode_into(payload);
    match f.qual {
        Qualities::Missing => payload.resize(payload.len().saturating_add(f.seq.len()), 0xFF),
        Qualities::Phred33(qual) => payload.extend(qual.iter().map(|q| q.wrapping_sub(33))),
        Qualities::Raw(qual) => payload.extend_from_slice(qual),
    }
    payload.extend_from_slice(f.aux);
    if long {
        payload.extend_from_slice(b"CGBI");
        let n_ops =
            u32::try_from(f.cigar.len()).map_err(|_| RecordEncodeError::RecordTooLarge { size })?;
        payload.extend_from_slice(&n_ops.to_le_bytes());
        payload.extend_from_slice(CigarOp::ops_as_bytes(f.cigar));
    }
    Ok(())
}

/// Size of BAM's fixed record fields.
const FIXED_LEN: usize = 32;
/// FLAG `0x4`.
const UNMAPPED: u16 = 0x4;

fn resolve(header: &BamHeader, name: &[u8]) -> Option<Tid> {
    header.tid(std::str::from_utf8(name).ok()?)
}

// r[impl sam.encode.references+2]
/// htslib's warning for an RNEXT the header lacks, once per process: a file
/// written against another header has one on every paired record.
pub(crate) fn warn_unknown_mate_reference(name: &[u8]) {
    static WARNED: AtomicBool = AtomicBool::new(false);
    if !WARNED.swap(true, Ordering::Relaxed) {
        warn!(
            rnext = %name.escape_ascii(),
            "unknown mate reference name; treated as unmapped (further ones not logged)"
        );
    }
}

/// The CIGAR's query length (M, I, S, =, X), as `bam_cigar2qlen`.
fn query_len(ops: &[CigarOp]) -> u64 {
    ops.iter()
        .filter(|op| op.consumes_query())
        .map(|op| u64::from(op.len()))
        .fold(0, u64::saturating_add)
}

/// The CIGAR's reference length (M, D, N, =, X), as `bam_cigar2rlen`.
fn ref_len(ops: &[CigarOp]) -> u64 {
    ops.iter()
        .filter(|op| op.consumes_ref())
        .map(|op| u64::from(op.len()))
        .fold(0, u64::saturating_add)
}

/// htslib's `hts_reg2bin(beg, end, 14, 5)`, stored in BAM's 16-bit field.
fn reg2bin(beg: i64, end: i64) -> u16 {
    let last = end.saturating_sub(1);
    // (shift, first bin of the level), leaf level first.
    for (shift, offset) in [(14, 4681), (17, 585), (20, 73), (23, 9), (26, 1)] {
        if beg >> shift == last >> shift {
            let bin = (beg >> shift).saturating_add(offset);
            #[expect(
                clippy::cast_possible_truncation,
                clippy::cast_sign_loss,
                reason = "htslib stores the bin in a uint16_t; past 2^29 it wraps, and so must we"
            )]
            return bin as u16;
        }
    }
    0
}

// r[impl sam.encode.numbers]
/// `[0-9]+`, or `None` (empty, another byte, or past `u64`).
fn parse_digits(bytes: &[u8]) -> Option<u64> {
    if bytes.is_empty() {
        return None;
    }
    bytes.iter().try_fold(0u64, |n, &b| {
        let digit = b.checked_sub(b'0').filter(|d| *d < 10)?;
        n.checked_mul(10)?.checked_add(u64::from(digit))
    })
}

/// Why `[-+]?[0-9]+` did not parse.
enum SignedError {
    Malformed,
    Overflow,
}

/// `[-+]?[0-9]+` as (has a leading `-`, magnitude).
fn parse_signed(bytes: &[u8]) -> Result<(bool, u64), SignedError> {
    let (negative, digits) = match bytes.split_first() {
        Some((b'-', rest)) => (true, rest),
        Some((b'+', rest)) => (false, rest),
        _ => (false, bytes),
    };
    if digits.is_empty() || !digits.iter().all(u8::is_ascii_digit) {
        return Err(SignedError::Malformed);
    }
    parse_digits(digits).map(|m| (negative, m)).ok_or(SignedError::Overflow)
}

fn signed_value(negative: bool, magnitude: u64) -> Option<i64> {
    let m = i64::try_from(magnitude).ok()?;
    Some(if negative { m.checked_neg()? } else { m })
}

/// FLAG: `[0-9]+` without a leading zero, at most 65535.
fn parse_flag(field: &[u8]) -> Result<u16, SamRecordError> {
    let leading_zero = field.len() > 1 && field.first() == Some(&b'0');
    parse_digits(field)
        .filter(|_| !leading_zero)
        .and_then(|v| u16::try_from(v).ok())
        .ok_or_else(|| SamRecordError::InvalidFlag { value: field.into() })
}

/// POS/PNEXT: `[0-9]+` up to 2^31, as a 0-based position (`None` for 0).
fn parse_position(
    field: &[u8],
    invalid: impl FnOnce(Box<[u8]>) -> SamRecordError,
) -> Result<Option<Pos0>, SamRecordError> {
    match parse_digits(field) {
        Some(0) => Ok(None),
        Some(one_based) => u32::try_from(one_based.wrapping_sub(1))
            .ok()
            .and_then(Pos0::new)
            .map(Some)
            .ok_or_else(|| invalid(field.into())),
        None => Err(invalid(field.into())),
    }
}

// r[impl sam.encode.cigar]
/// Parse a SAM CIGAR string into `out`, which is cleared first. `*` gives no
/// operations. Lengths must be below 2^28 (BAM's 28-bit field); the
/// operators are `MIDNSHP=X` and htslib's `B`.
///
/// # Errors
/// [`SamRecordError::EmptyCigar`], [`SamRecordError::InvalidCigarLength`],
/// [`SamRecordError::UnknownCigarOp`] or [`SamRecordError::MissingCigarOp`].
pub fn parse_cigar(text: &[u8], out: &mut Vec<CigarOp>) -> Result<(), SamRecordError> {
    out.clear();
    if text == b"*" {
        return Ok(());
    }
    if text.is_empty() {
        return Err(SamRecordError::EmptyCigar);
    }
    let mut rest = text;
    while !rest.is_empty() {
        let digits = rest.iter().take_while(|b| b.is_ascii_digit()).count();
        let (len, tail) = rest.split_at(digits);
        let Some((&op, tail)) = tail.split_first() else {
            return Err(SamRecordError::MissingCigarOp { value: len.into() });
        };
        let code = match op {
            b'M' => 0,
            b'I' => 1,
            b'D' => 2,
            b'N' => 3,
            b'S' => 4,
            b'H' => 5,
            b'P' => 6,
            b'=' => 7,
            b'X' => 8,
            b'B' => 9,
            _ => return Err(SamRecordError::UnknownCigarOp { op }),
        };
        let len_value = parse_digits(len)
            .and_then(|n| u32::try_from(n).ok())
            .filter(|n| *n < 1 << 28)
            .ok_or_else(|| SamRecordError::InvalidCigarLength { value: len.into() })?;
        out.push(CigarOp::from_bam_u32((len_value << 4) | code));
        rest = tail;
    }
    Ok(())
}

// r[impl sam.aux.parse]
/// Append the BAM encoding of SAM optional fields — `text` is everything
/// after the TAB that ends QUAL — to `out`. On error `out` holds a partial
/// encoding the caller must discard.
pub(crate) fn parse_aux_fields(text: &[u8], out: &mut Vec<u8>) -> Result<(), SamRecordError> {
    let mut rest = text;
    loop {
        let (field, next) = match text_scan::find_byte(rest, b'\t') {
            Some(tab) => (rest.get(..tab).unwrap_or_default(), rest.get(tab.saturating_add(1)..)),
            None => (rest, None),
        };
        // One TAB at the very end of the line is accepted, as htslib does.
        if field.is_empty() && next.is_none() {
            return Ok(());
        }
        parse_aux_field(field, out)?;
        match next {
            Some(next) => rest = next,
            None => return Ok(()),
        }
    }
}

/// A tag byte as htslib (with a signed `char`) accepts it: `!` to DEL.
fn is_tag_byte(b: u8) -> bool {
    (0x21..=0x7F).contains(&b)
}

// r[impl sam.aux.types]
fn parse_aux_field(field: &[u8], out: &mut Vec<u8>) -> Result<(), SamRecordError> {
    let malformed = || SamRecordError::InvalidAuxField { field: field.into() };
    let Some((&[t0, t1, colon1, ty, colon2], value)) = field.split_first_chunk::<5>() else {
        return Err(malformed());
    };
    if colon1 != b':' || colon2 != b':' || !is_tag_byte(t0) || !is_tag_byte(t1) {
        return Err(malformed());
    }
    let tag = [t0, t1];
    let invalid = || SamRecordError::InvalidAuxValue { tag, value: value.into() };
    if value.is_empty() && !matches!(ty, b'Z' | b'H') {
        return Err(invalid());
    }
    out.extend_from_slice(&tag);
    match ty {
        b'A' | b'a' | b'c' | b'C' => {
            // htslib (signed `char`) refuses a first byte <= TAB or >= 0x80.
            let &[ch] = value else { return Err(invalid()) };
            if !(0x0A..0x80).contains(&ch) {
                return Err(invalid());
            }
            out.extend_from_slice(&[b'A', ch]);
        }
        b'i' | b'I' => push_int(tag, value, out)?,
        b'f' => {
            let v = parse_float(value).ok_or_else(invalid)?;
            out.push(b'f');
            out.extend_from_slice(&to_f32(v).to_le_bytes());
        }
        b'd' => {
            let v = parse_float(value).ok_or_else(invalid)?;
            out.push(b'd');
            out.extend_from_slice(&v.to_le_bytes());
        }
        b'Z' => {
            out.push(b'Z');
            out.extend_from_slice(value);
            out.push(0);
        }
        b'H' => {
            if !value.len().is_multiple_of(2) || !value.iter().all(u8::is_ascii_hexdigit) {
                return Err(invalid());
            }
            out.push(b'H');
            out.extend_from_slice(value);
            out.push(0);
        }
        b'B' => push_array(tag, value, out)?,
        type_code => return Err(SamRecordError::UnknownAuxType { tag, type_code }),
    }
    Ok(())
}

/// An `i` value, unsigned-first like htslib: a leading `-` picks the smallest
/// of `c`/`s`/`i`, anything else the smallest of `C`/`S`/`I`.
fn push_int(tag: [u8; 2], value: &[u8], out: &mut Vec<u8>) -> Result<(), SamRecordError> {
    let out_of_range = || SamRecordError::AuxIntOutOfRange { tag, value: value.into() };
    let (negative, magnitude) = parse_signed(value).map_err(|e| match e {
        SignedError::Malformed => SamRecordError::InvalidAuxValue { tag, value: value.into() },
        SignedError::Overflow => out_of_range(),
    })?;
    if negative {
        let v = signed_value(true, magnitude).ok_or_else(out_of_range)?;
        if let Ok(v) = i8::try_from(v) {
            out.extend_from_slice(&[b'c', v.cast_unsigned()]);
        } else if let Ok(v) = i16::try_from(v) {
            out.push(b's');
            out.extend_from_slice(&v.to_le_bytes());
        } else {
            let v = i32::try_from(v).map_err(|_| out_of_range())?;
            out.push(b'i');
            out.extend_from_slice(&v.to_le_bytes());
        }
    } else if let Ok(v) = u8::try_from(magnitude) {
        out.extend_from_slice(&[b'C', v]);
    } else if let Ok(v) = u16::try_from(magnitude) {
        out.push(b'S');
        out.extend_from_slice(&v.to_le_bytes());
    } else {
        let v = u32::try_from(magnitude).map_err(|_| out_of_range())?;
        out.push(b'I');
        out.extend_from_slice(&v.to_le_bytes());
    }
    Ok(())
}

/// A `B` array: subtype, then nothing or `,`-separated elements.
fn push_array(tag: [u8; 2], value: &[u8], out: &mut Vec<u8>) -> Result<(), SamRecordError> {
    let invalid = |v: &[u8]| SamRecordError::InvalidAuxValue { tag, value: v.into() };
    let Some((&subtype, rest)) = value.split_first() else { return Err(invalid(value)) };
    if !matches!(subtype, b'c' | b'C' | b's' | b'S' | b'i' | b'I' | b'f') {
        return Err(SamRecordError::InvalidAuxSubtype { tag, subtype });
    }
    out.extend_from_slice(&[b'B', subtype]);
    let count_at = out.len();
    out.extend_from_slice(&[0; 4]);
    let mut count = 0usize;
    if !rest.is_empty() {
        let elements = rest.strip_prefix(b",").ok_or_else(|| invalid(value))?;
        for element in elements.split(|&b| b == b',') {
            count = count.saturating_add(1);
            push_element(subtype, element, out).ok_or_else(|| invalid(element))?;
        }
    }
    let n = u32::try_from(count).map_err(|_| SamRecordError::AuxArrayTooLarge { len: count })?;
    if let Some(slot) = out.get_mut(count_at..count_at.saturating_add(4)) {
        slot.copy_from_slice(&n.to_le_bytes());
    }
    Ok(())
}

/// One `B` element through its subtype. An unsigned subtype refuses a
/// leading `-` even on `-0`, as htslib does.
fn push_element(subtype: u8, element: &[u8], out: &mut Vec<u8>) -> Option<()> {
    if subtype == b'f' {
        out.extend_from_slice(&to_f32(parse_float(element)?).to_le_bytes());
        return Some(());
    }
    let (negative, magnitude) = parse_signed(element).ok()?;
    let v = signed_value(negative, magnitude)?;
    match subtype {
        b'c' => out.push(i8::try_from(v).ok()?.cast_unsigned()),
        b's' => out.extend_from_slice(&i16::try_from(v).ok()?.to_le_bytes()),
        b'i' => out.extend_from_slice(&i32::try_from(v).ok()?.to_le_bytes()),
        _ if negative => return None,
        b'C' => out.push(u8::try_from(v).ok()?),
        b'S' => out.extend_from_slice(&u16::try_from(v).ok()?.to_le_bytes()),
        _ => out.extend_from_slice(&u32::try_from(v).ok()?.to_le_bytes()),
    }
    Some(())
}

/// A decimal float, or `inf`/`infinity`/`nan` (any case, optional sign) —
/// exactly what `f64::from_str` accepts. A NaN is the default quiet NaN with
/// the written sign, as `strtod` returns it.
fn parse_float(value: &[u8]) -> Option<f64> {
    let v: f64 = std::str::from_utf8(value).ok()?.parse().ok()?;
    if v.is_nan() {
        let sign = if value.first() == Some(&b'-') { 1u64 << 63 } else { 0 };
        return Some(f64::from_bits(0x7FF8_0000_0000_0000 | sign));
    }
    Some(v)
}

/// `(float) strtod(...)`: round the double to the nearest `f32`.
fn to_f32(v: f64) -> f32 {
    if v.is_nan() {
        // `as` leaves a NaN's bits unspecified; strtod's NaN rounds to this.
        let sign = u32::from(v.is_sign_negative()) << 31;
        return f32::from_bits(0x7FC0_0000 | sign);
    }
    #[expect(
        clippy::cast_possible_truncation,
        reason = "htslib rounds the double to float, as `as` does"
    )]
    let v = v as f32;
    v
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used, clippy::indexing_slicing, reason = "test code")]
mod tests {
    use super::*;

    fn aux(field: &str) -> Result<Vec<u8>, SamRecordError> {
        let mut out = Vec::new();
        parse_aux_fields(field.as_bytes(), &mut out)?;
        Ok(out)
    }

    // r[verify sam.encode.bin]
    /// Values from htslib's own `hts_reg2bin`, including the unplaced bin.
    #[test]
    fn reg2bin_matches_htslib() {
        assert_eq!(reg2bin(-1, 0), 4680);
        assert_eq!(reg2bin(0, 1), 4681);
        assert_eq!(reg2bin(16383, 16384), 4681);
        assert_eq!(reg2bin(16383, 16385), 585);
        assert_eq!(reg2bin(0, 1 << 29), 0);
        // Past 2^29 the leaf bin no longer fits 16 bits and wraps.
        let beg = 2_000_000_000i64;
        // (2e9 >> 14) + 4681 = 126751, stored mod 2^16.
        assert_eq!(reg2bin(beg, beg + 1), 61_215);
    }

    // r[verify sam.aux.types]
    #[test]
    fn integers_are_unsigned_first() {
        assert_eq!(aux("XX:i:0").unwrap(), b"XXC\x00");
        assert_eq!(aux("XX:i:127").unwrap(), b"XXC\x7f");
        assert_eq!(aux("XX:i:255").unwrap(), b"XXC\xff");
        assert_eq!(aux("XX:i:+256").unwrap(), b"XXS\x00\x01");
        assert_eq!(aux("XX:I:65536").unwrap(), b"XXI\x00\x00\x01\x00");
        assert_eq!(aux("XX:i:4294967295").unwrap(), b"XXI\xff\xff\xff\xff");
        assert_eq!(aux("XX:i:-0").unwrap(), b"XXc\x00");
        assert_eq!(aux("XX:i:-128").unwrap(), b"XXc\x80");
        assert_eq!(aux("XX:i:-129").unwrap(), b"XXs\x7f\xff");
        assert_eq!(aux("XX:i:-32769").unwrap(), b"XXi\xff\x7f\xff\xff");
        assert_eq!(aux("XX:i:-2147483648").unwrap(), b"XXi\x00\x00\x00\x80");
        for v in ["4294967296", "-2147483649", "99999999999999999999999"] {
            assert!(
                matches!(aux(&format!("XX:i:{v}")), Err(SamRecordError::AuxIntOutOfRange { .. })),
                "{v}"
            );
        }
    }

    // r[verify sam.aux.types]
    // r[verify sam.encode.deviations]
    /// What htslib skips or defaults, seqair refuses.
    #[test]
    fn junk_is_rejected() {
        for field in [
            "NM:i:5abc",
            "NM:i:",
            "NM:i:+",
            "NM:i:abc",
            "XA:A:ab",
            "XA:A:",
            "XF:f:1.5x",
            "XF:f: 1",
            "XF:f:0x1p3",
            "XF:f:",
            "XH:H:ABC",
            "XH:H:GG",
            "XB:B:c,",
            "XB:B:c,1,",
            "XB:B:c1",
            "XB:B:c,128",
            "XB:B:C,-1",
            "XB:B:C,-0",
            "XB:B:S,1 ",
        ] {
            assert!(
                matches!(aux(field), Err(SamRecordError::InvalidAuxValue { .. })),
                "{field}: {:?}",
                aux(field)
            );
        }
        for field in ["NMxiy5", "NM:i", "N", "\x01M:i:5", "NM:i:5\t\tXX:i:1"] {
            assert!(matches!(aux(field), Err(SamRecordError::InvalidAuxField { .. })), "{field:?}");
        }
        assert!(matches!(
            aux("XX:q:1"),
            Err(SamRecordError::UnknownAuxType { type_code: b'q', .. })
        ));
        assert!(matches!(
            aux("XX:B:A"),
            Err(SamRecordError::InvalidAuxSubtype { subtype: b'A', .. })
        ));
        assert!(matches!(
            aux("XX:B:d,1"),
            Err(SamRecordError::InvalidAuxSubtype { subtype: b'd', .. })
        ));
    }

    // r[verify sam.aux.parse]
    #[test]
    fn a_trailing_tab_is_accepted_an_inner_empty_field_is_not() {
        assert_eq!(aux("").unwrap(), b"");
        assert_eq!(aux("NM:i:1\t").unwrap(), b"NMC\x01");
        assert!(aux("\tNM:i:1").is_err());
    }

    // r[verify sam.aux.types]
    #[test]
    fn every_type_encodes() {
        assert_eq!(aux("XA:A:x\tXa:a:y\tXc:c:5\tXC:C:!").unwrap(), b"XAAxXaAyXcA5XCA!");
        assert_eq!(aux("XZ:Z:\tXY:Z:a b").unwrap(), b"XZZ\x00XYZa b\x00");
        assert_eq!(aux("XH:H:\tXG:H:1aF0").unwrap(), b"XHH\x00XGH1aF0\x00");
        assert_eq!(aux("XF:f:1.5").unwrap(), [b"XFf".as_slice(), &1.5f32.to_le_bytes()].concat());
        assert_eq!(
            aux("XD:d:-2.25e3").unwrap(),
            [b"XDd".as_slice(), &(-2250f64).to_le_bytes()].concat()
        );
        assert_eq!(aux("XF:f:-nan").unwrap(), b"XFf\x00\x00\xc0\xff");
        assert_eq!(
            aux("XF:f:INF").unwrap(),
            [b"XFf".as_slice(), &f32::INFINITY.to_le_bytes()].concat()
        );
        assert_eq!(aux("XB:B:c").unwrap(), b"XBBc\x00\x00\x00\x00");
        assert_eq!(aux("XB:B:c,-1,+2").unwrap(), b"XBBc\x02\x00\x00\x00\xff\x02");
        assert_eq!(aux("XB:B:S,65535").unwrap(), b"XBBS\x01\x00\x00\x00\xff\xff");
        assert_eq!(
            aux("XB:B:f,1,.5").unwrap(),
            [b"XBBf\x02\x00\x00\x00".as_slice(), &1f32.to_le_bytes(), &0.5f32.to_le_bytes()]
                .concat()
        );
    }

    // r[verify sam.encode.cigar]
    #[test]
    fn cigar_errors() {
        let mut ops = vec![CigarOp::from_bam_u32(0)];
        parse_cigar(b"*", &mut ops).unwrap();
        assert!(ops.is_empty());
        parse_cigar(b"0M268435455N3B", &mut ops).unwrap();
        assert_eq!(
            ops.iter().map(|o| o.to_bam_u32()).collect::<Vec<_>>(),
            [0, (268_435_455 << 4) | 3, (3 << 4) | 9]
        );
        assert!(matches!(parse_cigar(b"", &mut ops), Err(SamRecordError::EmptyCigar)));
        assert!(matches!(
            parse_cigar(b"M", &mut ops),
            Err(SamRecordError::InvalidCigarLength { .. })
        ));
        assert!(matches!(
            parse_cigar(b"+5M", &mut ops),
            Err(SamRecordError::UnknownCigarOp { op: b'+' })
        ));
        assert!(matches!(
            parse_cigar(b"268435456M", &mut ops),
            Err(SamRecordError::InvalidCigarLength { .. })
        ));
        assert!(matches!(
            parse_cigar(b"5M3", &mut ops),
            Err(SamRecordError::MissingCigarOp { .. })
        ));
        assert!(matches!(
            parse_cigar(b"5Z", &mut ops),
            Err(SamRecordError::UnknownCigarOp { op: b'Z' })
        ));
        assert!(matches!(
            parse_cigar(b"*5M", &mut ops),
            Err(SamRecordError::UnknownCigarOp { op: b'*' })
        ));
    }
}
