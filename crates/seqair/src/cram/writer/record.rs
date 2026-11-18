//! One record into a slice's data series: the checks that refuse it, then
//! its read features and tags.

use super::CramWriteError;
use super::container::{Series, SliceBuilder, SliceKind};
use super::header::CramHeader;
use super::mates::MateRecord;
use crate::bam::aux::AuxTag;
use crate::bam::aux_data::advance_past_value;
use crate::bam::cigar::{CigarOp, CigarOpType};
use crate::bam::mate::MateFields;
use crate::cram::compression_header::SubstitutionMatrix;
use crate::cram::slice::MdNm;
use crate::reader::Tid;
use seqair_types::{BamFlags, IupacBase, MappingQuality, Pos0};

/// Where a record sits in the sort order: its reference and position.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Placement {
    pub tid: Option<Tid>,
    pub pos: Option<Pos0>,
}

impl Placement {
    /// Unplaced records sort after every placed one, in any order: htslib
    /// writes their `POS` (often 1) but nothing sorts on it.
    fn key(self) -> (bool, Option<Tid>, Option<Pos0>) {
        (self.tid.is_none(), self.tid, self.tid.and(self.pos))
    }
}

/// A record's fields, from whichever input it came.
pub(super) struct RecordView<'a> {
    pub tid: Option<Tid>,
    pub pos: Option<Pos0>,
    pub flags: BamFlags,
    pub mapq: MappingQuality,
    pub mate: MateFields,
    pub qname: &'a [u8],
    pub cigar: &'a [CigarOp],
    pub seq: &'a [IupacBase],
    /// As many as SEQ; `None` when the record has none.
    pub qual: Option<&'a [u8]>,
    pub aux: &'a [u8],
}

/// Buffers a raw record's fields are unpacked into.
#[derive(Default)]
pub(super) struct InputScratch {
    pub cigar: Vec<CigarOp>,
    pub seq: Vec<IupacBase>,
}

#[derive(Default)]
pub(super) struct Scratch {
    pub input: InputScratch,
    md: MdNm,
    line: Vec<u8>,
}

/// What [`check`] found out about a record it accepts.
#[derive(Debug, Clone)]
pub(super) struct Checked {
    pub placement: Placement,
    /// `RL`.
    pub read_len: i32,
    /// 0-based exclusive end of the reference span (an unmapped record
    /// spans its position).
    pub end: i64,
    /// Whether the record has an `MD` tag, and an `NM` tag.
    has_md: bool,
    has_nm: bool,
    /// An unmapped record whose CIGAR or MAPQ is not stored.
    pub drops_alignment: bool,
}

/// One aux tag's bytes.
struct Tag<'a> {
    head: [u8; 3],
    value: &'a [u8],
    /// A `Z` or `H` value holds a TAB (before its NUL).
    tab: bool,
}

/// A `Z`/`H` value's length up to its NUL (`None` without one) and whether
/// a TAB comes before it, eight bytes at a time: one pass where finding the
/// NUL and then checking for a TAB took two.
#[allow(clippy::arithmetic_side_effects, reason = "SWAR on u64 lanes; i + 8 ≤ len")]
#[allow(clippy::cast_possible_truncation, reason = "trailing_zeros / 8 < 8")]
fn string_value(rest: &[u8]) -> (Option<usize>, bool) {
    const LOW: u64 = 0x0101_0101_0101_0101;
    const HIGH: u64 = 0x8080_8080_8080_8080;
    // The lowest set bit marks the first zero byte exactly; higher ones
    // may be false positives, only above a true zero.
    let zeros = |w: u64| w.wrapping_sub(LOW) & !w & HIGH;
    let mut tab = false;
    let (words, tail) = rest.as_chunks::<8>();
    for (k, word) in words.iter().enumerate() {
        let w = u64::from_le_bytes(*word);
        let (nul, tabs) = (zeros(w), zeros(w ^ (LOW * u64::from(b'\t'))));
        if nul != 0 {
            // A TAB marker below the NUL's is a true TAB.
            let below = nul.isolate_lowest_one() - 1;
            return (Some(k * 8 + (nul.trailing_zeros() / 8) as usize), tab || tabs & below != 0);
        }
        tab |= tabs != 0;
    }
    let start = words.len() * 8;
    match tail.iter().position(|&b| b == 0) {
        Some(i) => (Some(start + i), tab || tail.get(..i).is_some_and(|t| t.contains(&b'\t'))),
        None => (None, tab || tail.contains(&b'\t')),
    }
}

/// The tags of `aux`, which [`check`] has validated.
fn tags(aux: &[u8]) -> impl Iterator<Item = Result<Tag<'_>, CramWriteError>> {
    let mut pos = 0usize;
    std::iter::from_fn(move || {
        let head: [u8; 3] = aux.get(pos..pos.checked_add(3)?)?.try_into().ok()?;
        let start = pos.saturating_add(3);
        let (end, tab) = match head[2] {
            b'Z' | b'H' => match string_value(aux.get(start..).unwrap_or_default()) {
                (Some(len), tab) => (Some(start.saturating_add(len).saturating_add(1)), tab),
                (None, _) => (None, false),
            },
            ty => (advance_past_value(aux, start, ty), false),
        };
        let Some(end) = end.filter(|&e| e <= aux.len()) else {
            pos = aux.len();
            return Some(Err(CramWriteError::InvalidAux));
        };
        pos = end;
        Some(Ok(Tag { head, value: aux.get(start..end).unwrap_or_default(), tab }))
    })
}

// r[impl cram_writer.sort_order]
// r[impl cram_writer.record.unmapped]
// r[impl cram_writer.record.name]
/// Everything that can refuse a record, so that encoding it cannot fail.
pub(super) fn check(
    view: &RecordView<'_>,
    header: &CramHeader,
    last: Option<Placement>,
) -> Result<Checked, CramWriteError> {
    if view.qname.is_empty() {
        return Err(CramWriteError::EmptyReadName);
    }
    // `RN` is NUL-terminated.
    if view.qname.contains(&0) {
        return Err(CramWriteError::NulInReadName);
    }
    let n_targets = header.n_targets();
    for tid in [view.tid, view.mate.tid].into_iter().flatten() {
        if tid.as_usize() >= n_targets {
            return Err(CramWriteError::TidOutOfRange { tid: tid.as_u32(), n_targets });
        }
    }

    // `AP` and `NP` are 1-based i32s, which `Pos0::MAX` overflows.
    for pos in [view.pos, view.mate.pos].into_iter().flatten() {
        if pos == Pos0::MAX {
            return Err(CramWriteError::ValueTooLarge { what: "position", value: pos.as_u64() });
        }
    }
    let placement = Placement { tid: view.tid, pos: view.pos };
    let unmapped = view.flags.is_unmapped();
    if (view.tid.is_some() && view.pos.is_none()) || (!unmapped && view.tid.is_none()) {
        return Err(CramWriteError::InvalidPlacement { placement });
    }
    if let (Some(tid), Some(pos)) = (view.tid, view.pos) {
        let length = header.target_len(tid);
        if pos.as_u64() >= length {
            return Err(CramWriteError::PositionPastReferenceEnd { pos, length });
        }
    }
    if let Some(previous) = last
        && placement.key() < previous.key()
    {
        return Err(CramWriteError::Unsorted { previous, record: placement });
    }

    // r[impl cram_writer.record.unmapped]
    // CRAM has no CIGAR or MAPQ for an unmapped record; they are dropped.
    let drops_alignment = unmapped && (!view.cigar.is_empty() || view.mapq.as_byte() != 0);
    let mut query_len = 0u64;
    let mut ref_len = 0u64;
    for op in view.cigar.iter().filter(|_| !unmapped) {
        let len = u64::from(op.len());
        match op.op_type() {
            CigarOpType::Match | CigarOpType::SeqMatch | CigarOpType::SeqMismatch => {
                query_len = query_len.saturating_add(len);
                ref_len = ref_len.saturating_add(len);
            }
            CigarOpType::Insertion | CigarOpType::SoftClip => {
                query_len = query_len.saturating_add(len);
            }
            CigarOpType::Deletion | CigarOpType::RefSkip => ref_len = ref_len.saturating_add(len),
            CigarOpType::HardClip | CigarOpType::Padding => {}
            CigarOpType::Unknown(code) => return Err(CramWriteError::UnknownCigarOp { code }),
        }
    }
    if !view.seq.is_empty() && !unmapped && query_len != view.seq.len() as u64 {
        return Err(CramWriteError::CigarSeqLengthMismatch {
            cigar: query_len,
            seq: view.seq.len(),
        });
    }
    if let Some(qual) = view.qual
        && qual.len() != view.seq.len()
    {
        return Err(CramWriteError::QualLengthMismatch { qual: qual.len(), seq: view.seq.len() });
    }
    // r[impl cram_writer.record.position]
    let read_len = if view.seq.is_empty() && !unmapped { query_len } else { view.seq.len() as u64 };
    let read_len = i32::try_from(read_len)
        .map_err(|_| CramWriteError::ValueTooLarge { what: "read length", value: read_len })?;

    let (mut has_md, mut has_nm) = (false, false);
    for tag in tags(view.aux) {
        let Tag { head: [a, b, ty], tab, .. } = tag?;
        let name = AuxTag([a, b]);
        // The tag dictionary is NUL-separated; SAM allows no other names.
        if !a.is_ascii_alphabetic() || !b.is_ascii_alphanumeric() {
            return Err(CramWriteError::InvalidTagName { tag: name });
        }
        match ty {
            b'A' | b'c' | b'C' | b's' | b'S' | b'i' | b'I' | b'f' | b'B' => {}
            b'Z' | b'H' if tab => {
                return Err(CramWriteError::TabInStringTag { tag: name });
            }
            b'Z' | b'H' => {}
            _ => return Err(CramWriteError::UnsupportedTagType { tag: name, ty: char::from(ty) }),
        }
        if [a, b] == *b"cF" {
            return Err(CramWriteError::ReservedTag { tag: name });
        }
        has_md |= [a, b] == *b"MD";
        has_nm |= [a, b] == *b"NM";
    }

    let pos = view.pos.map_or(0, Pos0::as_i64);
    let span = if unmapped { 1 } else { i64::try_from(ref_len).unwrap_or(i64::MAX) };
    Ok(Checked {
        placement,
        read_len,
        end: pos.saturating_add(span),
        has_md,
        has_nm,
        drops_alignment,
    })
}

/// The uppercased bases a slice's records are stored against: the FASTA's
/// reference from its start, or the reference a slice embeds from its own.
#[derive(Debug, Clone, Copy)]
pub(super) struct SliceRef<'a> {
    pub bases: &'a [u8],
    /// 0-based reference position of `bases[0]`.
    pub start: usize,
    /// Built from the slice's reads (`r[cram_writer.embedded_reference]`).
    pub embedded: bool,
    /// Every byte of `bases` is an IUPAC code, so no record needs checking.
    pub iupac: bool,
}

/// A reference's uppercased bases, checked once for a byte that is no IUPAC
/// code, rather than once a record over its span.
#[derive(Debug, Default)]
pub(super) struct RefBases {
    pub bases: Vec<u8>,
    pub iupac: bool,
}

impl RefBases {
    pub(super) fn new(bases: Vec<u8>) -> Self {
        let iupac = all_iupac(&bases);
        Self { bases, iupac }
    }
}

/// Whether every byte is an upper-case IUPAC code.
pub(super) fn all_iupac(bases: &[u8]) -> bool {
    bases.iter().all(|&b| IS_IUPAC.get(usize::from(b)).copied().unwrap_or(false))
}

impl<'a> SliceRef<'a> {
    pub(super) fn fasta(reference: &'a RefBases) -> Self {
        Self { bases: &reference.bases, start: 0, embedded: false, iupac: reference.iupac }
    }

    /// `pos`'s index in `bases`; past them (never) for one before them.
    fn index(self, pos: usize) -> usize {
        pos.checked_sub(self.start).unwrap_or(usize::MAX)
    }

    fn get(self, pos: usize) -> Option<u8> {
        self.bases.get(self.index(pos)).copied()
    }

    fn range(self, from: usize, to: usize) -> Option<&'a [u8]> {
        self.bases.get(self.index(from)..self.index(to))
    }

    /// 0-based exclusive end of the bases.
    fn end(self) -> usize {
        self.start.saturating_add(self.bases.len())
    }
}

/// The record's data series, appended to `slice`. `reference` is the bases
/// the slice is stored against, if any.
pub(super) fn encode(
    slice: &mut SliceBuilder,
    kind: SliceKind,
    view: &RecordView<'_>,
    checked: &Checked,
    reference: Option<SliceRef<'_>>,
    header: &CramHeader,
    scratch: &mut Scratch,
) {
    let unmapped = view.flags.is_unmapped();
    let seq_unknown = view.seq.is_empty();
    let pos1 = view.pos.map_or(0, |p| p.as_i64().saturating_add(1));

    // r[impl cram_writer.record.order]
    // r[impl cram_writer.record.flags]
    slice.int(Series::Bf, i32::from(view.flags.raw()));
    slice.int(Series::Rl, checked.read_len);
    slice.begin_record(kind, pos1, checked.read_len);
    slice.reach(checked.end);
    slice.reference_id(view.tid, pos1, if unmapped { pos1 } else { pos1.max(checked.end) });

    let name_start = slice.len(Series::Rn);
    slice.bytes(Series::Rn, view.qname);
    let name = name_start..slice.len(Series::Rn);
    slice.byte(Series::Rn, 0);

    // r[impl cram_writer.record.mate]
    // `CF` and the mate series are written with the slice, attached or not.
    let pos0 = pos1.saturating_sub(1);
    slice.push_mate(MateRecord {
        name,
        flags: view.flags,
        tid: view.tid,
        pos: view.pos,
        end: if unmapped {
            pos0.saturating_add(1)
        } else {
            checked.end.max(pos0.saturating_add(1))
        },
        mate: view.mate,
        seq_unknown,
    });

    // `MD`/`NM` are generated on decode only for these.
    let regenerates = !unmapped && !seq_unknown && reference.is_some();
    // seqair's decode writes a reference byte that is no IUPAC code into
    // `MD` as `N`, htslib's as it is: only drop what both regenerate.
    let iupac_reference = || {
        reference.is_some_and(|r| {
            let from = view.pos.map_or(0, Pos0::as_usize);
            let to = usize::try_from(checked.end).unwrap_or(usize::MAX).min(r.end());
            r.iupac || all_iupac(r.range(from, to).unwrap_or_default())
        })
    };
    let track = regenerates && (checked.has_md || checked.has_nm) && iupac_reference();
    let md = &mut scratch.md;
    md.start();
    let mut deletion_past_end = false;

    if unmapped {
        for &b in view.seq {
            slice.byte(Series::Ba, b.to_ascii());
        }
    } else {
        let features = Features {
            slice: &mut *slice,
            view,
            reference,
            md: track.then_some(&mut *md),
            count: 0,
            prev_pos: 0,
            deletion_past_end: false,
        };
        let count;
        (count, deletion_past_end) = features.walk();
        slice.int(Series::Fn, count);
        // r[impl cram_writer.record.mapq]
        slice.int(Series::Mq, i32::from(view.mapq.as_byte()));
    }
    // r[impl cram_writer.record.quality]
    match view.qual {
        Some(qual) => slice.bytes(Series::Qs, qual),
        None => {
            for _ in 0..checked.read_len {
                slice.byte(Series::Qs, 0xFF);
            }
        }
    }

    // r[impl cram_writer.md_nm+2]
    // htslib's decoder fails a record whose deletion runs past its
    // reference's end when it generates either tag (`beyond_slice`).
    let track = track && !deletion_past_end;
    let (generated_md, generated_nm) = if track { md.finish() } else { (&[][..], 0) };
    // r[impl cram_writer.record.read_group]
    let mut read_group = -1;
    let line = &mut scratch.line;
    line.clear();
    for tag in tags(view.aux).map_while(Result::ok) {
        let [a, b, ty] = tag.head;
        let drop = match &[a, b] {
            b"RG" if ty == b'Z' && read_group < 0 => {
                match header.read_group(tag.value.strip_suffix(&[0]).unwrap_or(tag.value)) {
                    Some(rg) => {
                        read_group = rg;
                        true
                    }
                    None => false,
                }
            }
            b"MD" => track && ty == b'Z' && tag.value.strip_suffix(&[0]) == Some(generated_md),
            b"NM" => track && int_value(ty, tag.value) == Some(i64::from(generated_nm)),
            _ => false,
        };
        if !drop {
            line.extend_from_slice(&tag.head);
            slice.tag(tag.head, tag.value);
        }
    }
    // r[impl cram_writer.md_nm.past_reference_end]
    // r[impl cram_writer.embedded_reference.md_nm]
    let embedded = reference.is_some_and(|r| r.embedded);
    if regenerates && (deletion_past_end || embedded) {
        let flags = u8::from(!checked.has_md) | u8::from(!checked.has_nm) << 1;
        if flags != 0 {
            line.extend_from_slice(b"cFC");
            slice.tag(*b"cFC", &[flags]);
        }
    }
    slice.int(Series::Rg, read_group);
    let line_id = slice.tag_line(line);
    slice.int(Series::Tl, line_id);
}

/// An integer tag's value.
fn int_value(ty: u8, value: &[u8]) -> Option<i64> {
    Some(match (ty, value) {
        (b'c', &[v]) => i64::from(v.cast_signed()),
        (b'C', &[v]) => i64::from(v),
        (b's', &[a, b]) => i64::from(i16::from_le_bytes([a, b])),
        (b'S', &[a, b]) => i64::from(u16::from_le_bytes([a, b])),
        (b'i', &[a, b, c, d]) => i64::from(i32::from_le_bytes([a, b, c, d])),
        (b'I', &[a, b, c, d]) => i64::from(u32::from_le_bytes([a, b, c, d])),
        _ => return None,
    })
}

/// A mapped record's CIGAR walked against the reference into read
/// features, and the `MD`/`NM` a decode would generate from them.
struct Features<'a, 's> {
    slice: &'s mut SliceBuilder,
    view: &'a RecordView<'a>,
    reference: Option<SliceRef<'a>>,
    md: Option<&'s mut MdNm>,
    count: i32,
    /// 1-based read position of the previous feature.
    prev_pos: i32,
    /// Whether a deletion runs past the reference's end.
    deletion_past_end: bool,
}

impl Features<'_, '_> {
    /// Emit the features; returns their count, and whether a deletion runs
    /// past the reference's end.
    // r[impl cram_writer.features]
    fn walk(mut self) -> (i32, bool) {
        let view = self.view;
        let seq_unknown = view.seq.is_empty();
        let mut read = 0usize;
        let mut refpos = view.pos.map_or(0, Pos0::as_usize);
        // A zero-length operation is no feature; it reads back as nothing.
        for op in view.cigar.iter().filter(|op| op.len() > 0) {
            let len = op.len() as usize;
            match op.op_type() {
                CigarOpType::Match | CigarOpType::SeqMatch | CigarOpType::SeqMismatch => {
                    if !seq_unknown {
                        self.aligned(read, refpos, len);
                    }
                    read = read.saturating_add(len);
                    refpos = refpos.saturating_add(len);
                }
                CigarOpType::Insertion => {
                    if len == 1 {
                        self.feature(b'i', read);
                        let base = self.base_ascii(read);
                        self.slice.byte(Series::Ba, base);
                    } else {
                        self.feature(b'I', read);
                        for i in read..read.saturating_add(len) {
                            let base = self.base_ascii(i);
                            self.slice.byte(Series::In, base);
                        }
                        self.slice.byte(Series::In, 0);
                    }
                    if let Some(md) = self.md.as_deref_mut() {
                        md.insertion(len);
                    }
                    read = read.saturating_add(len);
                }
                CigarOpType::Deletion => {
                    self.feature(b'D', read);
                    self.slice.int(Series::Dl, op.len().cast_signed());
                    if let (Some(md), Some(r)) = (self.md.as_deref_mut(), self.reference) {
                        md.deletion(r.bases, r.index(refpos), len);
                    }
                    refpos = refpos.saturating_add(len);
                    self.deletion_past_end |= self.reference.is_some_and(|r| refpos > r.end());
                }
                CigarOpType::RefSkip => {
                    self.feature(b'N', read);
                    self.slice.int(Series::Rs, op.len().cast_signed());
                    refpos = refpos.saturating_add(len);
                }
                CigarOpType::SoftClip => {
                    self.feature(b'S', read);
                    for i in read..read.saturating_add(len) {
                        let base = self.base_ascii(i);
                        self.slice.byte(Series::Sc, base);
                    }
                    self.slice.byte(Series::Sc, 0);
                    read = read.saturating_add(len);
                }
                CigarOpType::HardClip => {
                    self.feature(b'H', read);
                    self.slice.int(Series::Hc, op.len().cast_signed());
                }
                CigarOpType::Padding => {
                    self.feature(b'P', read);
                    self.slice.int(Series::Pd, op.len().cast_signed());
                }
                // `check` refused these.
                CigarOpType::Unknown(_) => {}
            }
        }
        (self.count, self.deletion_past_end)
    }

    /// A feature at 0-based read position `read`.
    fn feature(&mut self, code: u8, read: usize) {
        let pos = i32::try_from(read.saturating_add(1)).unwrap_or(i32::MAX);
        self.slice.byte(Series::Fc, code);
        self.slice.int(Series::Fp, pos.saturating_sub(self.prev_pos));
        self.prev_pos = pos;
        self.count = self.count.saturating_add(1);
    }

    /// A read base as the features store it; `N` for SEQ `*`.
    fn base_ascii(&self, read: usize) -> u8 {
        self.view.seq.get(read).map_or(b'N', |b| b.to_ascii())
    }

    /// `len` aligned bases from read position `read` and reference position
    /// `refpos`: implicit where they match, substitutions where not; `b`
    /// features without a reference.
    fn aligned(&mut self, read: usize, refpos: usize, len: usize) {
        let view = self.view;
        let Some(reference) = self.reference else {
            // r[impl cram_writer.reference_free+2]
            self.feature(b'b', read);
            self.slice.int(Series::BbLen, i32::try_from(len).unwrap_or(i32::MAX));
            for i in read..read.saturating_add(len) {
                let base = self.base_ascii(i);
                self.slice.byte(Series::Bb, base);
            }
            return;
        };
        let mut run_start = refpos;
        let (read_end, ref_end) = (read.saturating_add(len), refpos.saturating_add(len));
        let mut i = 0;
        while i < len {
            // The run of matches up to the next mismatch, a vector at a time:
            // the same letter, as htslib compares them (a reference letter
            // that is no IUPAC code, or `U`, reads as `N` but is no match).
            let seq = view.seq.get(read.saturating_add(i)..read_end).unwrap_or_default();
            let g = refpos.saturating_add(i);
            let refs = reference.range(g, ref_end.min(reference.end())).unwrap_or_default();
            i = i.saturating_add(IupacBase::ascii_prefix_len(seq, refs));
            if i >= len {
                break;
            }
            let g = refpos.saturating_add(i);
            let r = read.saturating_add(i);
            i = i.saturating_add(1);
            let base = view.seq.get(r).copied().unwrap_or(IupacBase::N);
            let ref_base = reference.get(g);
            if ref_base == Some(base.to_ascii()) {
                continue;
            }
            self.matches(run_start, g);
            run_start = g.saturating_add(1);
            // r[impl cram_writer.substitution+2]
            // r[impl cram_writer.past_reference_end]
            match ref_base.and_then(|rb| substitution(rb, base)) {
                Some((row, alt)) => {
                    self.feature(b'X', r);
                    self.slice.substitution(row, alt);
                    if let Some(md) = self.md.as_deref_mut() {
                        md.substitution(ref_base);
                    }
                }
                None => {
                    self.feature(b'B', r);
                    self.slice.byte(Series::Ba, base.to_ascii());
                    let qual = view.qual.and_then(|q| q.get(r)).copied().unwrap_or(0xFF);
                    self.slice.byte(Series::Qs, qual);
                    if let Some(md) = self.md.as_deref_mut() {
                        md.read_base(ref_base);
                    }
                }
            }
        }
        self.matches(run_start, refpos.saturating_add(len));
    }

    /// The implicit matches over reference positions `from..to`.
    fn matches(&mut self, from: usize, to: usize) {
        if let (Some(md), Some(reference)) = (self.md.as_deref_mut(), self.reference)
            && to > from
        {
            md.matches(reference.range(from, to), to.saturating_sub(from));
        }
    }
}

/// Whether a byte is an upper-case IUPAC code, as a table: the check runs
/// on every reference base a record spans.
const IS_IUPAC: [bool; 256] = {
    let mut table = [false; 256];
    let mut i = 0;
    while i < IupacBase::ASCII.len() {
        #[allow(clippy::indexing_slicing, reason = "const: i < 16, a u8 indexes 256")]
        {
            table[IupacBase::ASCII[i] as usize] = true;
        }
        i += 1;
    }
    table
};

/// The matrix row of reference base `ref_base` (`ACGTN`), and the index of
/// `read` among that row's other bases in `ACGTN` order: `None` where a `B`
/// feature must store it. Only a literal `ACGTN` reference base takes a
/// substitution, so every decoder picks the same matrix row.
fn substitution(ref_base: u8, read: IupacBase) -> Option<(u8, u8)> {
    let row: u8 = match ref_base {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        b'N' => 4,
        _ => return None,
    };
    let bases = SubstitutionMatrix::READ_BASES.get(usize::from(row))?;
    let alt = bases.iter().position(|&b| b == read.to_ascii())?;
    Some((row, u8::try_from(alt).ok()?))
}

#[cfg(test)]
mod tests {
    #![allow(clippy::arithmetic_side_effects, clippy::indexing_slicing, reason = "test code")]
    use super::*;
    use hegel::TestCase;
    use hegel::generators as gs;

    // r[verify cram_writer.tags]
    /// The SWAR scan finds the first NUL and a TAB before it as a byte loop
    /// does, wherever they fall in or across the eight-byte words.
    #[hegel::test]
    fn string_values_end_at_the_first_nul(tc: TestCase) {
        let len = tc.draw(gs::integers::<usize>().max_value(40));
        // Mostly letters; NULs, TABs and high bytes where drawn.
        let bytes: Vec<u8> = (0..len)
            .map(|_| match tc.draw(gs::integers::<u8>().max_value(9)) {
                0 => 0,
                1 => b'\t',
                2 => tc.draw(gs::integers::<u8>()),
                _ => b'A' + tc.draw(gs::integers::<u8>().max_value(25)),
            })
            .collect();
        let nul = bytes.iter().position(|&b| b == 0);
        let tab = bytes[..nul.unwrap_or(bytes.len())].contains(&b'\t');
        let (got_nul, got_tab) = string_value(&bytes);
        assert_eq!(got_nul, nul, "{bytes:?}");
        // Without a NUL the value is refused whatever the TAB says.
        if nul.is_some() {
            assert_eq!(got_tab, tab, "{bytes:?}");
        }
    }
}
