//! Slices on a reference the FASTA lacks: their records are held until the
//! slice is complete, then stored against a reference built from them, or
//! without one where that is smaller.

use super::record::{Checked, RecordView};
use crate::bam::aux::{AuxValue, find_tag};
use crate::bam::cigar::{CigarOp, CigarOpType};
use crate::bam::mate::MateFields;
use crate::reader::Tid;
use seqair_types::{BamFlags, IupacBase, MappingQuality, Pos0};

/// A record held until its slice is complete.
pub(super) struct Held {
    tid: Option<Tid>,
    pos: Option<Pos0>,
    flags: BamFlags,
    mapq: MappingQuality,
    mate: MateFields,
    qname: Vec<u8>,
    cigar: Vec<CigarOp>,
    seq: Vec<IupacBase>,
    qual: Option<Vec<u8>>,
    aux: Vec<u8>,
    pub checked: Checked,
}

impl Held {
    pub(super) fn new(view: &RecordView<'_>, checked: Checked) -> Self {
        Self {
            tid: view.tid,
            pos: view.pos,
            flags: view.flags,
            mapq: view.mapq,
            mate: view.mate,
            qname: view.qname.to_vec(),
            cigar: view.cigar.to_vec(),
            seq: view.seq.to_vec(),
            qual: view.qual.map(<[u8]>::to_vec),
            aux: view.aux.to_vec(),
            checked,
        }
    }

    pub(super) fn view(&self) -> RecordView<'_> {
        RecordView {
            tid: self.tid,
            pos: self.pos,
            flags: self.flags,
            mapq: self.mapq,
            mate: self.mate,
            qname: &self.qname,
            cigar: &self.cigar,
            seq: &self.seq,
            qual: self.qual.as_deref(),
            aux: &self.aux,
        }
    }

    /// Each `M`/`=`/`X` base with its 0-based reference position.
    fn aligned(&self) -> impl Iterator<Item = (usize, IupacBase)> + '_ {
        self.walk().filter_map(|(pos, base)| Some((pos, base?)))
    }

    /// Each reference position an `M`/`=`/`X` (with its read base) or a `D`
    /// covers, in order; nothing for a record without bases.
    fn walk(&self) -> impl Iterator<Item = (usize, Option<IupacBase>)> + '_ {
        let mapped = !self.flags.is_unmapped() && !self.seq.is_empty();
        let mut read = 0usize;
        let mut refpos = self.pos.map_or(0, Pos0::as_usize);
        self.cigar.iter().filter(move |_| mapped).flat_map(move |op| {
            let len = op.len() as usize;
            let (r, g) = (read, refpos);
            // Whether the operation consumes the read, the reference, and
            // whether its positions are walked.
            let (reads, refs, walked) = match op.op_type() {
                CigarOpType::Match | CigarOpType::SeqMatch | CigarOpType::SeqMismatch => {
                    (true, true, true)
                }
                CigarOpType::Insertion | CigarOpType::SoftClip => (true, false, false),
                CigarOpType::Deletion => (false, true, true),
                CigarOpType::RefSkip => (false, true, false),
                CigarOpType::HardClip | CigarOpType::Padding | CigarOpType::Unknown(_) => {
                    (false, false, false)
                }
            };
            if reads {
                read = read.saturating_add(len);
            }
            if refs {
                refpos = refpos.saturating_add(len);
            }
            let n = if walked { len } else { 0 };
            (0..n).map(move |i| {
                let base = if reads { self.seq.get(r.saturating_add(i)).copied() } else { None };
                (g.saturating_add(i), base)
            })
        })
    }

    /// The reference bases the record's `MD` gives, with their positions:
    /// a match run's are the read's, a mismatch's and a deletion's the
    /// letters. `None` without `MD`, or where it does not fit the CIGAR.
    fn md_bases(&self) -> Option<Vec<(usize, u8)>> {
        let AuxValue::String(md) = find_tag(&self.aux, *b"MD")? else { return None };
        let mut walk = self.walk();
        let mut out = Vec::new();
        let mut md = md.iter().copied().peekable();
        while let Some(c) = md.next() {
            if c.is_ascii_digit() {
                let mut n = usize::from(c.wrapping_sub(b'0'));
                while let Some(d) = md.next_if(u8::is_ascii_digit) {
                    n = n.checked_mul(10)?.checked_add(usize::from(d.wrapping_sub(b'0')))?;
                }
                for _ in 0..n {
                    let (pos, base) = walk.next()?;
                    out.push((pos, base?.to_ascii()));
                }
            } else if c == b'^' {
                while let Some(letter) = md.next_if(u8::is_ascii_alphabetic) {
                    let (pos, None) = walk.next()? else { return None };
                    out.push((pos, letter));
                }
            } else {
                let (pos, Some(_)) = walk.next()? else { return None };
                out.push((pos, c));
            }
        }
        walk.next().is_none().then_some(out)
    }
}

/// The 0-based start and exclusive end of `held`'s slice, its first record's
/// position to the furthest end, cut at `length`
/// (`r[cram_writer.slice_extent+2]`).
pub(super) fn extent(held: &[Held], length: u64) -> (usize, usize) {
    let start = held.first().and_then(|h| h.pos).map_or(0, Pos0::as_usize);
    let end = held.iter().map(|h| h.checked.end).max().unwrap_or(0);
    let end = usize::try_from(end.min(i64::try_from(length).unwrap_or(i64::MAX))).unwrap_or(0);
    (start, end.max(start))
}

/// The longest extent a slice embeds a reference for: its base counts take
/// 8 bytes a position while it is built (htslib's limit is near 22 Mi).
const MAX_EMBEDDED: usize = 1 << 24;

// r[impl cram_writer.embedded_reference]
/// The reference to embed for `held`, over `start..end`, as htslib's
/// `cram_generate_reference` builds it: where a record's `MD` gives a base,
/// that one (the last record's); elsewhere the most frequent of `ACGT`
/// among the aligned bases there, `N` where none is. `None` where the
/// aligned bases number under 1% of the extent (htslib's limit, which also
/// spares the counts of a nearly empty contig) or the extent is longer than
/// [`MAX_EMBEDDED`].
pub(super) fn consensus(held: &[Held], start: usize, end: usize) -> Option<Vec<u8>> {
    let len = end.saturating_sub(start);
    let aligned = held.iter().flat_map(Held::aligned).count();
    if aligned == 0 || len == 0 || len > MAX_EMBEDDED || aligned.saturating_mul(100) < len {
        return None;
    }
    let mut counts = vec![[0u16; 4]; len];
    for (pos, base) in held.iter().flat_map(Held::aligned) {
        let i = match base {
            IupacBase::A => 0,
            IupacBase::C => 1,
            IupacBase::G => 2,
            IupacBase::T => 3,
            _ => continue,
        };
        if let Some(c) = pos.checked_sub(start).and_then(|p| counts.get_mut(p))
            && let Some(n) = c.get_mut(i)
        {
            *n = n.saturating_add(1);
        }
    }
    let mut bases: Vec<u8> = counts
        .iter()
        .map(|c| {
            // The first of the most frequent: `max_by_key` keeps the last.
            let (best, &n) = c.iter().enumerate().rev().max_by_key(|&(_, n)| *n)?;
            (n > 0).then(|| b"ACGT".get(best).copied()).flatten()
        })
        .map(|b| b.unwrap_or(b'N'))
        .collect();
    for (pos, base) in held.iter().filter_map(Held::md_bases).flatten() {
        let base = base.to_ascii_uppercase();
        if let (b'A' | b'C' | b'G' | b'T' | b'N', Some(slot)) =
            (base, pos.checked_sub(start).and_then(|p| bases.get_mut(p)))
        {
            *slot = base;
        }
    }
    Some(bases)
}
