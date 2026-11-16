//! BAM's escape for CIGARs of more than 65535 operations (`[SAM1] §4.2.2`).
//!
//! BAM stores the op count in 16 bits. A longer CIGAR is written as the
//! placeholder `<l_seq>S<ref_len>N` with the real ops in a trailing `CG:B:I`
//! tag. That is purely a matter of the BAM wire format: everywhere else in
//! seqair a CIGAR is its real ops. This module is the only place that knows
//! the placeholder — [`find`] and [`take`] undo it on read, as htslib's
//! `bam_tag2cigar` does, and [`WireCigar`] applies it on write, as
//! `bam_write1` does.

use super::aux::{find_tag_offset, tag_end};
use super::cigar::{self, CIGAR_S, CigarOp, CigarOpType, PackedCigar};
use super::record::DecodeError;
use std::borrow::Cow;
use std::ops::Range;
use thiserror::Error;

/// The bytes before a `CG` tag's ops: the tag, `B`, the subtype and the count.
pub(crate) const TAG_HEADER_LEN: usize = 8;

/// The most ops BAM's CIGAR field holds.
const MAX_INLINE_OPS: usize = u16::MAX as usize;

/// `bam_tag2cigar`'s bound on a `CG` array.
const MAX_CG_OPS: usize = 1 << 29;

/// The bound on the placeholder's lengths: a CIGAR op holds 28 bits.
const MAX_OP_LEN: u64 = 1 << 28;

/// A long CIGAR found in a record's tags.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct LongCigar<'a> {
    /// The real ops.
    pub ops: PackedCigar<'a>,
    /// The whole `CG` tag, as a range of the aux block it was found in.
    pub tag: Range<usize>,
}

/// Whether `stored` begins as BAM's placeholder does, for a placed record of
/// `seq_len` bases: the one condition of `bam_tag2cigar` that does not need
/// the tags, and the only cost a record without a long CIGAR pays.
pub(crate) fn starts_like_placeholder(
    stored: PackedCigar<'_>,
    seq_len: usize,
    placed: bool,
) -> bool {
    let Some(placeholder) = u32::try_from(seq_len).ok().and_then(|n| n.checked_mul(16)) else {
        return false;
    };
    placed && stored.first().is_some_and(|op| op.to_bam_u32() == placeholder | u32::from(CIGAR_S))
}

// r[impl bam.raw_record.long_cigar]
// r[impl bam.raw_record.cigar_seq]
/// htslib's `bam_tag2cigar`: a placed record whose stored CIGAR starts with
/// exactly `<seq_len>S`, and whose first `CG` tag is a `B:I`/`B:i` array of
/// at least as many and fewer than 2^29 ops, has those ops as its CIGAR.
/// They must account for the sequence, as htslib checks there, so no reader
/// can take a CIGAR from a tag without that check.
pub(crate) fn find<'a>(
    stored: PackedCigar<'_>,
    seq_len: usize,
    placed: bool,
    aux: &'a [u8],
) -> Result<Option<LongCigar<'a>>, DecodeError> {
    if !starts_like_placeholder(stored, seq_len, placed) {
        return Ok(None);
    }
    let Some(long) = cg_tag(stored.len(), aux) else { return Ok(None) };
    let cigar_query_len = cigar::query_len(long.ops);
    let seq_len = u32::try_from(seq_len).unwrap_or(u32::MAX);
    if seq_len != 0 && cigar_query_len != u64::from(seq_len) {
        return Err(DecodeError::CigarQueryLenMismatch { cigar_query_len, seq_len });
    }
    Ok(Some(long))
}

/// The first `CG` tag of `aux`, if it is an array `bam_tag2cigar` takes.
fn cg_tag(n_stored: usize, aux: &[u8]) -> Option<LongCigar<'_>> {
    let at = find_tag_offset(aux, *b"CG")?;
    let header = aux.get(at..at.checked_add(TAG_HEADER_LEN)?)?;
    let &[_, _, b'B', b'I' | b'i', n0, n1, n2, n3] = header else { return None };
    let count = usize::try_from(u32::from_le_bytes([n0, n1, n2, n3])).ok()?;
    if count < n_stored || count >= MAX_CG_OPS {
        return None;
    }
    let ops_start = at.checked_add(TAG_HEADER_LEN)?;
    let end = ops_start.checked_add(count.checked_mul(4)?)?;
    let ops = PackedCigar::whole_ops(aux.get(ops_start..end)?);
    Some(LongCigar { ops, tag: at..end })
}

impl LongCigar<'_> {
    // r[impl bam.raw_record.cg_hidden]
    /// The tags of `aux` — the block this was found in — once the CIGAR has
    /// left them: without any `CG` tag. `bam_tag2cigar` drops only the one
    /// it took; another one stays a tag that the next reader may take for
    /// the CIGAR (when the real one starts with `<seq_len>S`, or BAM stores
    /// it behind the placeholder again), so the record would not read back
    /// as itself. Borrowed in the common case, the tag last and alone.
    pub(crate) fn remaining_tags<'b>(&self, aux: &'b [u8]) -> Cow<'b, [u8]> {
        let before = aux.get(..self.tag.start).unwrap_or_default();
        let after = aux.get(self.tag.end..).unwrap_or_default();
        if after.is_empty() && find_tag_offset(before, *b"CG").is_none() {
            return before.into();
        }
        let mut kept = Vec::with_capacity(aux.len());
        extend_without_cg(aux, &mut kept);
        kept.into()
    }
}

/// [`find`] on a record being decoded into buffers (SAM, CRAM): make the
/// `CG` tag's ops the CIGAR and drop the tags, as
/// [`LongCigar::remaining_tags`] leaves them. Whether it did.
pub(crate) fn take(
    cigar: &mut Vec<CigarOp>,
    aux: &mut Vec<u8>,
    seq_len: usize,
    placed: bool,
) -> Result<bool, DecodeError> {
    let Some(long) = find(PackedCigar::from_ops(cigar), seq_len, placed, aux)? else {
        return Ok(false);
    };
    cigar.clear();
    long.ops.extend_into(cigar);
    let rest = long.remaining_tags(aux).into_owned();
    *aux = rest;
    Ok(true)
}

/// Drop every `CG` tag from `aux`: for a reader that decodes `CG` whatever
/// tags were asked for, because the CIGAR may be in one, and must not hand
/// one out that was not asked for and was not the CIGAR.
pub(crate) fn drop_cg_tags(aux: &mut Vec<u8>) {
    if find_tag_offset(aux, *b"CG").is_some() {
        let mut kept = Vec::with_capacity(aux.len());
        extend_without_cg(aux, &mut kept);
        *aux = kept;
    }
}

// r[impl sam.encode.long_cigar+1]
/// Why BAM cannot store a CIGAR of more than 65535 operations, which goes in
/// a `CG` tag behind the placeholder `<l_seq>S<ref_len>N`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Error)]
#[non_exhaustive]
pub enum LongCigarError {
    /// The placeholder's `N` op holds the reference length in 28 bits.
    #[error(
        "a CIGAR of more than 65535 operations covers {ref_len} reference bases; BAM allows < 2^28"
    )]
    RefLenTooLong { ref_len: u64 },

    /// The placeholder's `S` op holds the sequence length in 28 bits.
    #[error("a CIGAR of more than 65535 operations on {seq_len} bases; BAM allows < 2^28")]
    SeqTooLong { seq_len: usize },

    /// A reader takes the `CG` tag as the CIGAR only for a record with a
    /// reference and a position, so an unplaced record's would be lost.
    #[error("an unplaced record cannot carry a CIGAR of {n_ops} operations in BAM")]
    Unplaced { n_ops: usize },
}

/// A CIGAR as a BAM record stores it: inline, or behind the placeholder with
/// its ops in a trailing `CG:B:I` tag.
#[derive(Debug, Clone, Copy)]
pub(crate) enum WireCigar<'a> {
    /// In the CIGAR field. `looks_escaped` when it starts as the placeholder
    /// does, so that a reader looks for a `CG` tag to take in its place.
    Inline {
        ops: &'a [CigarOp],
        looks_escaped: bool,
    },
    Escaped {
        placeholder: [CigarOp; 2],
        ops: &'a [CigarOp],
    },
}

impl<'a> WireCigar<'a> {
    /// How BAM stores `ops`, the CIGAR of a record of `seq_len` bases,
    /// `placed` when it has a reference and a position.
    pub(crate) fn new(
        ops: &'a [CigarOp],
        seq_len: usize,
        placed: bool,
    ) -> Result<Self, LongCigarError> {
        if ops.len() <= MAX_INLINE_OPS {
            let looks_escaped =
                starts_like_placeholder(PackedCigar::from_ops(ops), seq_len, placed);
            return Ok(Self::Inline { ops, looks_escaped });
        }
        if !placed {
            return Err(LongCigarError::Unplaced { n_ops: ops.len() });
        }
        let ref_len = cigar::ref_len(ops);
        let n = u32::try_from(ref_len)
            .ok()
            .filter(|&n| u64::from(n) < MAX_OP_LEN)
            .ok_or(LongCigarError::RefLenTooLong { ref_len })?;
        let s = u32::try_from(seq_len)
            .ok()
            .filter(|&s| u64::from(s) < MAX_OP_LEN)
            .ok_or(LongCigarError::SeqTooLong { seq_len })?;
        let placeholder =
            [CigarOp::new(CigarOpType::SoftClip, s), CigarOp::new(CigarOpType::RefSkip, n)];
        Ok(Self::Escaped { placeholder, ops })
    }

    /// BAM's 16-bit `n_cigar_op`.
    pub(crate) fn n_cigar_op(&self) -> u16 {
        match self {
            // `new` keeps inline CIGARs within the field.
            Self::Inline { ops, .. } => u16::try_from(ops.len()).unwrap_or(u16::MAX),
            Self::Escaped { .. } => 2,
        }
    }

    /// The bytes the CIGAR field and the `CG` tag add to a record.
    pub(crate) fn wire_len(&self) -> usize {
        match self {
            Self::Inline { ops, .. } => ops.len().saturating_mul(4),
            Self::Escaped { ops, .. } => {
                ops.len().saturating_mul(4).saturating_add(TAG_HEADER_LEN).saturating_add(8)
            }
        }
    }

    /// Append the CIGAR field.
    pub(crate) fn write_cigar(&self, out: &mut Vec<u8>) {
        let ops: &[CigarOp] = match self {
            Self::Inline { ops, .. } => ops,
            Self::Escaped { placeholder, .. } => placeholder,
        };
        out.extend_from_slice(CigarOp::ops_as_bytes(ops));
    }

    /// Append the aux block `aux` and, for an escaped CIGAR, the `CG` tag
    /// after it. A `CG` tag in `aux` that a reader would take as the CIGAR
    /// is left out, with every other `CG` tag: it would stand in for the real
    /// ops. (htslib's `bam_write1` keeps it, and loses the CIGAR.)
    pub(crate) fn write_aux(&self, aux: &[u8], out: &mut Vec<u8>) {
        match self {
            Self::Inline { ops, looks_escaped: true } if cg_tag(ops.len(), aux).is_some() => {
                extend_without_cg(aux, out);
            }
            Self::Inline { .. } => out.extend_from_slice(aux),
            Self::Escaped { ops, .. } => {
                extend_without_cg(aux, out);
                out.extend_from_slice(b"CGBI");
                // At most a record's worth of ops (`MAX_RECORD_SIZE`), far
                // below 2^32.
                out.extend_from_slice(&u32::try_from(ops.len()).unwrap_or(u32::MAX).to_le_bytes());
                out.extend_from_slice(CigarOp::ops_as_bytes(ops));
            }
        }
    }
}

/// Append `aux` to `out` without its `CG` tags. A malformed tail is copied
/// as it is, so reading it still reports it.
pub(crate) fn extend_without_cg(aux: &[u8], out: &mut Vec<u8>) {
    let mut at = 0usize;
    while let Some(next) = find_tag_offset(aux.get(at..).unwrap_or_default(), *b"CG") {
        let start = at.saturating_add(next);
        out.extend_from_slice(aux.get(at..start).unwrap_or_default());
        let Ok(end) = tag_end(aux, start) else {
            out.extend_from_slice(aux.get(start..).unwrap_or_default());
            return;
        };
        at = end;
    }
    out.extend_from_slice(aux.get(at..).unwrap_or_default());
}

#[cfg(test)]
#[allow(
    clippy::unwrap_used,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    reason = "test code"
)]
mod tests {
    use super::*;
    use hegel::TestCase;
    use hegel::generators as gs;

    fn op(t: CigarOpType, len: u32) -> CigarOp {
        CigarOp::new(t, len)
    }

    fn cg(subtype: u8, ops: &[CigarOp]) -> Vec<u8> {
        let mut tag = vec![b'C', b'G', b'B', subtype];
        tag.extend_from_slice(&u32::try_from(ops.len()).unwrap().to_le_bytes());
        tag.extend_from_slice(CigarOp::ops_as_bytes(ops));
        tag
    }

    /// A CIGAR that may start with a soft clip, as a real one can.
    #[hegel::composite]
    fn cigar(tc: &TestCase) -> Vec<CigarOp> {
        let n = tc.draw(gs::integers::<usize>().min_value(1).max_value(8));
        (0..n)
            .map(|_| {
                let t = tc.draw(gs::sampled_from(&[
                    CigarOpType::Match,
                    CigarOpType::Insertion,
                    CigarOpType::Deletion,
                    CigarOpType::SoftClip,
                    CigarOpType::RefSkip,
                ]));
                op(t, tc.draw(gs::integers::<u32>().min_value(1).max_value(9)))
            })
            .collect()
    }

    /// One tag of a block around the `CG` tag: another tag, or a stale `CG`.
    #[derive(Debug, Clone)]
    enum Tag {
        Other(Vec<u8>),
        StaleCg(Vec<u8>),
    }

    #[hegel::composite]
    fn tags(tc: &TestCase) -> Vec<Tag> {
        let n = tc.draw(gs::integers::<usize>().max_value(3));
        (0..n)
            .map(|_| match tc.draw(gs::integers::<u8>().max_value(3)) {
                0 => Tag::Other(b"NMCx".to_vec()),
                1 => Tag::Other(b"RGZgrp\0".to_vec()),
                // An array whose bytes spell `CG`: no tag starts inside a value.
                2 => Tag::Other(b"XBBC\x02\x00\x00\x00CG".to_vec()),
                _ => Tag::StaleCg(cg(b'I', &[op(CigarOpType::Match, 7)])),
            })
            .collect()
    }

    fn bytes(tags: &[Tag], keep_cg: bool) -> Vec<u8> {
        tags.iter()
            .filter_map(|t| match t {
                Tag::Other(b) => Some(b.as_slice()),
                Tag::StaleCg(b) => keep_cg.then_some(b.as_slice()),
            })
            .flatten()
            .copied()
            .collect()
    }

    // r[verify bam.raw_record.long_cigar]
    // r[verify bam.raw_record.cg_hidden]
    /// `take` makes the first `CG` tag's ops the CIGAR exactly when
    /// `bam_tag2cigar` would, and then leaves the other tags, in order, byte
    /// for byte, and no `CG` tag; `find` reports the tag's exact range.
    #[hegel::test]
    fn take_moves_the_tag_into_the_cigar(tc: TestCase) {
        let real = tc.draw(cigar());
        let seq_len = usize::try_from(cigar::query_len(&real)).unwrap();
        let before = tc.draw(tags());
        let after = tc.draw(tags());
        let subtype = tc.draw(gs::sampled_from(b"Ii"));
        let placed = tc.draw(gs::booleans());
        let placeholder =
            vec![op(CigarOpType::SoftClip, seq_len as u32), op(CigarOpType::RefSkip, 1)];
        let tag = cg(subtype, &real);
        let aux = [bytes(&before, true), tag.clone(), bytes(&after, true)].concat();

        // The oracle: the first `CG` decides, and must hold at least as many
        // ops as the placeholder's two.
        let first_is_ours = !before.iter().any(|t| matches!(t, Tag::StaleCg(_)));
        let resolvable = placed && first_is_ours && real.len() >= placeholder.len();
        let found = find(PackedCigar::from_ops(&placeholder), seq_len, placed, &aux).unwrap();
        assert_eq!(found.is_some(), resolvable);
        if let Some(found) = &found {
            assert_eq!(found.tag, before_len(&before)..before_len(&before) + tag.len());
            assert_eq!(found.ops.to_vec(), real);
            assert_eq!(
                &*found.remaining_tags(&aux),
                [bytes(&before, false), bytes(&after, false)].concat()
            );
        }

        let (mut c, mut a) = (placeholder.clone(), aux.clone());
        assert_eq!(take(&mut c, &mut a, seq_len, placed).unwrap(), resolvable);
        if resolvable {
            assert_eq!(c, real);
            assert_eq!(a, [bytes(&before, false), bytes(&after, false)].concat());
        } else {
            assert_eq!((c, a), (placeholder, aux));
        }
    }

    fn before_len(tags: &[Tag]) -> usize {
        bytes(tags, true).len()
    }

    // r[verify bam.raw_record.long_cigar]
    // r[verify bam.raw_record.cigar_seq]
    /// Every condition of `bam_tag2cigar` is needed: an unplaced record, a
    /// first op other than `<seq_len>S`, another subtype, fewer ops than
    /// stored, or a truncated array leave the CIGAR alone; and the ops must
    /// account for the sequence, unless there is none.
    #[test]
    fn find_needs_every_condition() {
        let real = vec![
            op(CigarOpType::Match, 2),
            op(CigarOpType::Insertion, 1),
            op(CigarOpType::Match, 1),
        ];
        let stored = [op(CigarOpType::SoftClip, 4), op(CigarOpType::RefSkip, 3)];
        let p = PackedCigar::from_ops(&stored);
        let some = |r: Result<Option<LongCigar<'_>>, DecodeError>| r.unwrap().is_some();
        assert!(some(find(p, 4, true, &cg(b'I', &real))));
        assert!(some(find(p, 4, true, &cg(b'i', &real))));
        assert!(!some(find(p, 4, false, &cg(b'I', &real))));
        assert!(!some(find(p, 5, true, &cg(b'I', &real))));
        assert!(!some(find(p, 4, true, &cg(b'S', &real))));
        assert!(!some(find(p, 4, true, &cg(b'I', &real[..1]))));
        let mut short = cg(b'I', &real);
        short.pop();
        assert!(!some(find(p, 4, true, &short)));
        // Exactly as many ops as stored is enough.
        let two = [op(CigarOpType::Match, 3), op(CigarOpType::Insertion, 1)];
        assert!(some(find(p, 4, true, &cg(b'I', &two))));

        let over = vec![op(CigarOpType::Match, 3), op(CigarOpType::Insertion, 3)];
        assert!(matches!(
            find(p, 4, true, &cg(b'I', &over)),
            Err(DecodeError::CigarQueryLenMismatch { cigar_query_len: 6, seq_len: 4 })
        ));
        // No sequence (`*`): nothing to check against, as in htslib.
        let none = [op(CigarOpType::SoftClip, 0), op(CigarOpType::RefSkip, 3)];
        assert!(some(find(PackedCigar::from_ops(&none), 0, true, &cg(b'I', &over))));
    }

    // r[verify sam.encode.long_cigar+1]
    /// An escaped CIGAR writes the placeholder and one trailing `CG` tag, any
    /// stale `CG` dropped, and `find` reads it back; an inline one keeps a
    /// stale `CG` unless the CIGAR starts as a placeholder does, where a
    /// reader would take it for the CIGAR. Either way the record reads back
    /// with its own CIGAR.
    #[hegel::test]
    fn wire_cigars_read_back(tc: TestCase) {
        let long = tc.draw(gs::booleans());
        let n = if long {
            tc.draw(
                gs::integers::<usize>().min_value(MAX_INLINE_OPS - 1).max_value(MAX_INLINE_OPS + 2),
            )
        } else {
            tc.draw(gs::integers::<usize>().min_value(1).max_value(4))
        };
        let clipped = tc.draw(gs::booleans());
        let ops: Vec<CigarOp> = (0..n)
            .map(|i| match i {
                0 if clipped => op(CigarOpType::SoftClip, 1),
                _ if i % 2 == 0 => op(CigarOpType::Match, 1),
                _ => op(CigarOpType::Insertion, 1),
            })
            .collect();
        // A lone `1S` is as long as the sequence: the CIGAR starts as the
        // placeholder does.
        let seq_len = usize::try_from(cigar::query_len(&ops)).unwrap();
        let placed = tc.draw(gs::booleans());
        let stale = tc.draw(tags());
        let aux = bytes(&stale, true);

        let wire = match WireCigar::new(&ops, seq_len, placed) {
            Err(LongCigarError::Unplaced { n_ops }) => {
                assert!(!placed && n > MAX_INLINE_OPS && n_ops == n);
                return;
            }
            other => other.unwrap(),
        };
        let mut cigar_field = Vec::new();
        wire.write_cigar(&mut cigar_field);
        let mut tags = Vec::new();
        wire.write_aux(&aux, &mut tags);
        let stored = PackedCigar::from_bam_bytes(&cigar_field).unwrap();
        assert_eq!(stored.len(), usize::from(wire.n_cigar_op()));

        let read = find(stored, seq_len, placed, &tags).unwrap();
        let cigar = read.as_ref().map_or(stored, |l| l.ops);
        assert_eq!(cigar.to_vec(), ops, "the record reads back with its own CIGAR");
        assert_eq!(read.is_some(), n > MAX_INLINE_OPS);
        let tags_seen =
            read.as_ref().map_or(Cow::from(tags.as_slice()), |l| l.remaining_tags(&tags));
        let a_reader_would_look = placed && clipped && n == 1;
        let keeps_stale = n <= MAX_INLINE_OPS && !a_reader_would_look;
        assert_eq!(&*tags_seen, bytes(&stale, keeps_stale));
        assert_eq!(
            cigar_field.len() + tags.len() - bytes(&stale, keeps_stale).len(),
            wire.wire_len()
        );
    }

    #[test]
    fn a_placeholder_needs_28_bit_lengths() {
        let ops = vec![op(CigarOpType::Deletion, (1 << 28) - 1); MAX_INLINE_OPS + 1];
        assert!(matches!(WireCigar::new(&ops, 0, true), Err(LongCigarError::RefLenTooLong { .. })));
        let ops = vec![op(CigarOpType::Match, 1); MAX_INLINE_OPS + 1];
        assert!(matches!(
            WireCigar::new(&ops, 1 << 28, true),
            Err(LongCigarError::SeqTooLong { seq_len }) if seq_len == 1 << 28
        ));
    }
}
