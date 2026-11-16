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
use std::ops::Range;

/// The bytes before a `CG` tag's ops: the tag, `B`, the subtype and the count.
pub(crate) const TAG_HEADER_LEN: usize = 8;

/// The most ops BAM's CIGAR field holds.
const MAX_INLINE_OPS: usize = u16::MAX as usize;

/// `bam_tag2cigar`'s bound on a `CG` array.
const MAX_CG_OPS: usize = 1 << 29;

/// A long CIGAR found in a record's tags.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct LongCigar<'a> {
    /// The real ops.
    pub ops: PackedCigar<'a>,
    /// The whole `CG` tag, as a range of the aux block it was found in.
    pub tag: Range<usize>,
}

// r[impl bam.raw_record.long_cigar]
/// htslib's `bam_tag2cigar` test: a placed record whose stored CIGAR starts
/// with exactly `<seq_len>S`, and whose first `CG` tag is a `B:I`/`B:i` array
/// of at least as many and fewer than 2^29 ops. The one cheap check — the
/// first op — comes first, so a record without a long CIGAR costs a compare.
pub(crate) fn find<'a>(
    stored: PackedCigar<'_>,
    seq_len: usize,
    placed: bool,
    aux: &'a [u8],
) -> Option<LongCigar<'a>> {
    let placeholder = u32::try_from(seq_len).ok()?.checked_mul(16)? | u32::from(CIGAR_S);
    if stored.first()?.to_bam_u32() != placeholder || !placed {
        return None;
    }
    let at = find_tag_offset(aux, *b"CG")?;
    let header = aux.get(at..at.checked_add(TAG_HEADER_LEN)?)?;
    let &[_, _, b'B', b'I' | b'i', n0, n1, n2, n3] = header else { return None };
    let count = usize::try_from(u32::from_le_bytes([n0, n1, n2, n3])).ok()?;
    if count < stored.len() || count >= MAX_CG_OPS {
        return None;
    }
    let ops_start = at.checked_add(TAG_HEADER_LEN)?;
    let end = ops_start.checked_add(count.checked_mul(4)?)?;
    let ops = PackedCigar::whole_ops(aux.get(ops_start..end)?);
    Some(LongCigar { ops, tag: at..end })
}

impl LongCigar<'_> {
    /// The tags of `aux` — the block this was found in — once the CIGAR has
    /// left them: without this tag, and past 65535 ops, without any other
    /// `CG` tag either, since BAM then stores the CIGAR in a `CG` tag again
    /// and a reader would take another one in its place. Borrowed when the
    /// tag is the last one, which is where every writer puts it.
    pub(crate) fn remaining_tags<'b>(&self, aux: &'b [u8]) -> std::borrow::Cow<'b, [u8]> {
        if self.ops.len() > MAX_INLINE_OPS {
            let mut kept = Vec::with_capacity(aux.len());
            extend_without_cg(aux, &mut kept);
            return kept.into();
        }
        let (before, after) = (
            aux.get(..self.tag.start).unwrap_or_default(),
            aux.get(self.tag.end..).unwrap_or_default(),
        );
        if after.is_empty() { before.into() } else { [before, after].concat().into() }
    }
}

/// [`find`] on a record being decoded into buffers (SAM, CRAM): make the
/// `CG` tag's ops the CIGAR and drop the tag, as `bam_tag2cigar` does.
/// Whether it did.
pub(crate) fn take(
    cigar: &mut Vec<CigarOp>,
    aux: &mut Vec<u8>,
    seq_len: usize,
    placed: bool,
) -> bool {
    let Some(LongCigar { ops, tag }) = find(PackedCigar::from_ops(cigar), seq_len, placed, aux)
    else {
        return false;
    };
    cigar.clear();
    ops.extend_into(cigar);
    aux.drain(tag);
    true
}

/// The reference span of a CIGAR too long for BAM's escape, whose `N` op
/// holds it in 28 bits.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct RefLenTooLong {
    pub ref_len: u64,
}

// r[impl sam.encode.long_cigar+1]
/// A CIGAR as a BAM record stores it: inline, or behind the placeholder with
/// its ops in a trailing `CG:B:I` tag.
#[derive(Debug, Clone, Copy)]
pub(crate) enum WireCigar<'a> {
    Inline(&'a [CigarOp]),
    Escaped { placeholder: [CigarOp; 2], ops: &'a [CigarOp] },
}

impl<'a> WireCigar<'a> {
    /// How BAM stores `ops`, the CIGAR of a record of `seq_len` bases.
    pub(crate) fn new(ops: &'a [CigarOp], seq_len: u32) -> Result<Self, RefLenTooLong> {
        if ops.len() <= MAX_INLINE_OPS {
            return Ok(Self::Inline(ops));
        }
        let ref_len = cigar::ref_len(ops);
        let Ok(n) = u32::try_from(ref_len) else { return Err(RefLenTooLong { ref_len }) };
        if n >= 1 << 28 || seq_len >= 1 << 28 {
            return Err(RefLenTooLong { ref_len });
        }
        let placeholder =
            [CigarOp::new(CigarOpType::SoftClip, seq_len), CigarOp::new(CigarOpType::RefSkip, n)];
        Ok(Self::Escaped { placeholder, ops })
    }

    /// BAM's 16-bit `n_cigar_op`.
    pub(crate) fn n_cigar_op(&self) -> u16 {
        match self {
            // `new` keeps inline CIGARs within the field.
            Self::Inline(ops) => u16::try_from(ops.len()).unwrap_or(u16::MAX),
            Self::Escaped { .. } => 2,
        }
    }

    /// The bytes the CIGAR field and the `CG` tag add to a record.
    pub(crate) fn wire_len(&self) -> usize {
        match self {
            Self::Inline(ops) => ops.len().saturating_mul(4),
            Self::Escaped { ops, .. } => {
                ops.len().saturating_mul(4).saturating_add(TAG_HEADER_LEN).saturating_add(8)
            }
        }
    }

    /// Append the CIGAR field.
    pub(crate) fn write_cigar(&self, out: &mut Vec<u8>) {
        let ops: &[CigarOp] = match self {
            Self::Inline(ops) => ops,
            Self::Escaped { placeholder, .. } => placeholder,
        };
        out.extend_from_slice(CigarOp::ops_as_bytes(ops));
    }

    /// Append the aux block `aux` and, for an escaped CIGAR, the `CG` tag
    /// after it. A `CG` tag already in `aux` is left out then: a reader takes
    /// the first one, so a stale tag would win over the real ops.
    pub(crate) fn write_aux(&self, aux: &[u8], out: &mut Vec<u8>) {
        let Self::Escaped { ops, .. } = self else {
            out.extend_from_slice(aux);
            return;
        };
        extend_without_cg(aux, out);
        out.extend_from_slice(b"CGBI");
        // At most a record's worth of ops (`MAX_RECORD_SIZE`), far below 2^32.
        out.extend_from_slice(&u32::try_from(ops.len()).unwrap_or(u32::MAX).to_le_bytes());
        out.extend_from_slice(CigarOp::ops_as_bytes(ops));
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
    reason = "test code"
)]
mod tests {
    use super::*;
    use hegel::TestCase;
    use hegel::generators as gs;

    fn op(t: CigarOpType, len: u32) -> CigarOp {
        CigarOp::new(t, len)
    }

    fn cg_tag(subtype: u8, ops: &[CigarOp]) -> Vec<u8> {
        let mut tag = vec![b'C', b'G', b'B', subtype];
        tag.extend_from_slice(&u32::try_from(ops.len()).unwrap().to_le_bytes());
        tag.extend_from_slice(CigarOp::ops_as_bytes(ops));
        tag
    }

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

    /// Tags around a `CG` tag: none of them is `CG`.
    #[hegel::composite]
    fn other_tags(tc: &TestCase) -> Vec<u8> {
        let n = tc.draw(gs::integers::<usize>().max_value(3));
        let mut tags = Vec::new();
        for _ in 0..n {
            match tc.draw(gs::integers::<u8>().max_value(2)) {
                0 => tags.extend_from_slice(b"NMCx"),
                1 => tags.extend_from_slice(b"RGZgrp\0"),
                _ => {
                    // An array, so a `CG`-looking byte run inside a value
                    // never starts a tag.
                    tags.extend_from_slice(b"XBBC\x02\x00\x00\x00CG");
                }
            }
        }
        tags
    }

    // r[verify bam.raw_record.long_cigar]
    /// `take` makes the tag's ops the CIGAR and leaves the other tags, in
    /// order, byte for byte; `find` reports the tag's exact range.
    #[hegel::test]
    fn take_moves_the_tag_into_the_cigar(tc: TestCase) {
        let real = tc.draw(cigar());
        let seq_len = u32::try_from(cigar::query_len(&real)).unwrap();
        let before = tc.draw(other_tags());
        let after = tc.draw(other_tags());
        let subtype = tc.draw(gs::sampled_from(b"Ii"));
        let placeholder = vec![op(CigarOpType::SoftClip, seq_len), op(CigarOpType::RefSkip, 1)];
        let tag = cg_tag(subtype, &real);
        let aux = [before.as_slice(), &tag, &after].concat();

        let found = find(PackedCigar::from_ops(&placeholder), seq_len as usize, true, &aux);
        // Asking for at least as many ops as the placeholder has.
        let resolvable = real.len() >= placeholder.len();
        assert_eq!(found.is_some(), resolvable);
        if let Some(found) = found {
            assert_eq!(found.tag, before.len()..before.len() + tag.len());
            assert_eq!(found.ops.to_vec(), real);
        }

        let (mut c, mut a) = (placeholder.clone(), aux.clone());
        assert_eq!(take(&mut c, &mut a, seq_len as usize, true), resolvable);
        if resolvable {
            assert_eq!(c, real);
            assert_eq!(a, [before, after].concat());
        } else {
            assert_eq!((c, a), (placeholder, aux));
        }
    }

    // r[verify bam.raw_record.long_cigar]
    /// Every condition of `bam_tag2cigar` is needed: an unplaced record, a
    /// first op other than `<seq_len>S`, another subtype, or fewer ops than
    /// stored leave the CIGAR alone.
    #[test]
    fn find_needs_every_condition() {
        let real = vec![
            op(CigarOpType::Match, 2),
            op(CigarOpType::Insertion, 1),
            op(CigarOpType::Match, 1),
        ];
        let stored = [op(CigarOpType::SoftClip, 4), op(CigarOpType::RefSkip, 3)];
        let p = PackedCigar::from_ops(&stored);
        assert!(find(p, 4, true, &cg_tag(b'I', &real)).is_some());
        assert!(find(p, 4, false, &cg_tag(b'I', &real)).is_none());
        assert!(find(p, 5, true, &cg_tag(b'I', &real)).is_none());
        assert!(find(p, 4, true, &cg_tag(b'S', &real)).is_none());
        assert!(find(p, 4, true, &cg_tag(b'I', &real[..1])).is_none());
        // A truncated array is no tag.
        let mut short = cg_tag(b'I', &real);
        short.pop();
        assert!(find(p, 4, true, &short).is_none());
    }

    // r[verify sam.encode.long_cigar+1]
    /// An escaped CIGAR writes the placeholder and one trailing `CG` tag, a
    /// stale `CG` dropped, and reading it back with `find` gives the ops.
    #[hegel::test]
    fn escaped_cigars_read_back(tc: TestCase) {
        let n = tc.draw(
            gs::integers::<usize>().min_value(MAX_INLINE_OPS - 1).max_value(MAX_INLINE_OPS + 2),
        );
        let ops: Vec<CigarOp> = (0..n)
            .map(|i| op(if i % 2 == 0 { CigarOpType::Match } else { CigarOpType::Insertion }, 1))
            .collect();
        let seq_len = u32::try_from(cigar::query_len(&ops)).unwrap();
        let mut aux = tc.draw(other_tags());
        if tc.draw(gs::booleans()) {
            aux.extend_from_slice(&cg_tag(b'I', &[op(CigarOpType::Match, 7)]));
            aux.extend(tc.draw(other_tags()));
        }

        let wire = WireCigar::new(&ops, seq_len).unwrap();
        let mut cigar_field = Vec::new();
        wire.write_cigar(&mut cigar_field);
        let mut tags = Vec::new();
        wire.write_aux(&aux, &mut tags);
        let mut kept = Vec::new();
        extend_without_cg(&aux, &mut kept);
        let kept = if n > MAX_INLINE_OPS { kept } else { aux.clone() };
        assert_eq!(cigar_field.len() + tags.len(), wire.wire_len() + kept.len());

        let stored = PackedCigar::from_bam_bytes(&cigar_field).unwrap();
        assert_eq!(stored.len(), usize::from(wire.n_cigar_op()));
        match find(stored, seq_len as usize, true, &tags) {
            Some(found) => {
                assert!(n > MAX_INLINE_OPS);
                assert_eq!(found.ops.to_vec(), ops);
                assert_eq!(found.tag.end, tags.len(), "the CG tag is last");
                assert_eq!(&tags[..found.tag.start], kept.as_slice());
            }
            None => {
                assert!(n <= MAX_INLINE_OPS);
                assert_eq!(stored.to_vec(), ops);
                assert_eq!(tags, aux);
            }
        }
    }

    #[test]
    fn a_span_past_28_bits_cannot_be_escaped() {
        let ops = vec![op(CigarOpType::Deletion, (1 << 28) - 1); MAX_INLINE_OPS + 1];
        assert!(WireCigar::new(&ops, 0).is_err());
    }
}
