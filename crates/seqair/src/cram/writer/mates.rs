//! Mates attached within a slice: a pair is stored as a link from the first
//! record to the second, from which the decode rebuilds both records' mate
//! fields, where storing the fields themselves (detached) costs a few bytes
//! a record. A pair is attached only where every decoder rebuilds the
//! fields as the records have them.

use crate::bam::mate::MateFields;
use crate::cram::slice::ChainExtent;
use crate::cram::varint::write_itf8;
use crate::reader::Tid;
use rustc_hash::FxHashMap;
use seqair_types::{BamFlags, Pos0};
use std::ops::Range;

/// A record as the decode's mate resolution sees it, and the mate fields it
/// must give back.
#[derive(Debug, Clone)]
pub(super) struct MateRecord {
    /// The record's name in the slice's `RN` block.
    pub name: Range<usize>,
    pub flags: BamFlags,
    pub tid: Option<Tid>,
    /// `None` unplaced (`AP` 0).
    pub pos: Option<Pos0>,
    /// 0-based exclusive end: the reference span, at least one base.
    pub end: i64,
    pub mate: MateFields,
    /// SEQ is `*` (CRAM flag `0x8`).
    pub seq_unknown: bool,
}

/// The series the mates decide: `CF` and the mate series.
#[derive(Default)]
pub(super) struct MateSeries {
    pub cf: Vec<u8>,
    pub mf: Vec<u8>,
    pub ns: Vec<u8>,
    pub np: Vec<u8>,
    pub ts: Vec<u8>,
    pub nf: Vec<u8>,
}

const QUALITY_ARRAY: i32 = 0x1;
const DETACHED: i32 = 0x2;
const MATE_DOWNSTREAM: i32 = 0x4;
const NO_SEQ: i32 = 0x8;

/// How a record is stored.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Link {
    Detached,
    /// The first of a pair, its mate the record at this index.
    Downstream(usize),
    /// The second of a pair.
    Upstream,
}

// r[impl cram_writer.mates+2]
/// `CF` and the mate series of `records`, whose names are in `names` (the
/// `RN` block).
pub(super) fn series(records: &[MateRecord], names: &[u8]) -> MateSeries {
    // The first and second record of each name, and whether there are more.
    let mut by_name: FxHashMap<&[u8], (usize, Option<usize>, bool)> = FxHashMap::default();
    for (i, r) in records.iter().enumerate() {
        let name = names.get(r.name.clone()).unwrap_or_default();
        by_name
            .entry(name)
            .and_modify(|(_, second, many)| {
                *many |= second.is_some();
                *second = Some(i);
            })
            .or_insert((i, None, false));
    }
    let mut links = vec![Link::Detached; records.len()];
    for &(a, b, many) in by_name.values() {
        let (Some(b), false) = (b, many) else { continue };
        if let (Some(ra), Some(rb)) = (records.get(a), records.get(b))
            && rebuilds(ra, rb)
        {
            if let Some(link) = links.get_mut(a) {
                *link = Link::Downstream(b);
            }
            if let Some(link) = links.get_mut(b) {
                *link = Link::Upstream;
            }
        }
    }

    let mut out = MateSeries::default();
    for (i, (r, link)) in records.iter().zip(&links).enumerate() {
        let cf = QUALITY_ARRAY | if r.seq_unknown { NO_SEQ } else { 0 };
        match *link {
            Link::Downstream(b) => {
                write_itf8(&mut out.cf, cf | MATE_DOWNSTREAM);
                let between = b.saturating_sub(i).saturating_sub(1);
                write_itf8(&mut out.nf, i32::try_from(between).unwrap_or(i32::MAX));
            }
            Link::Upstream => write_itf8(&mut out.cf, cf),
            Link::Detached => {
                write_itf8(&mut out.cf, cf | DETACHED);
                // r[impl cram_writer.record.mate]
                // Only a paired record's mate bits: the decode adds 0x1 with
                // 0x20, which BF already holds.
                let mut mf = 0;
                if r.flags.is_paired() && r.flags.is_mate_reverse() {
                    mf |= 0x1;
                }
                if r.flags.is_mate_unmapped() {
                    mf |= 0x2;
                }
                write_itf8(&mut out.mf, mf);
                let tid = r.mate.tid.map_or(-1, |t| i32::try_from(t.as_u32()).unwrap_or(-1));
                write_itf8(&mut out.ns, tid);
                write_itf8(&mut out.np, r.mate.pos.map_or(0, |p| p.as_i32().saturating_add(1)));
                write_itf8(&mut out.ts, r.mate.template_len.map_or(0, i32::from));
            }
        }
    }
    out
}

/// Whether the decode rebuilds `a` and its downstream mate `b` as they are:
/// each one's FLAG holds the bits the decode adds, its RNEXT and PNEXT are
/// the other's reference and position, and its TLEN is the one the decode
/// computes (seqair's `resolve_mates` for a ring of two, htslib's
/// `cram_decode_slice_xref`). CRAM 3 has no other way to store a TLEN
/// (htslib's CRAM flag `0x10` is CRAM 4's), so a pair that would need one
/// is detached.
///
/// Both must be mapped on one reference: for a pair with an unmapped read,
/// or one across references (in a multi-reference slice), noodles computes
/// a TLEN where htslib and seqair give 0.
fn rebuilds(a: &MateRecord, b: &MateRecord) -> bool {
    // htslib's decode gives a pair across references TLEN 0, noodles one.
    if a.flags.is_unmapped() || b.flags.is_unmapped() || a.tid != b.tid {
        return false;
    }
    let start = |r: &MateRecord| r.pos.map_or(-1, Pos0::as_i64);
    let extent = ChainExtent {
        left: start(a).min(start(b)),
        right: a.end.max(b.end),
        left_cnt: u32::from(start(a) == start(a).min(start(b)))
            .saturating_add(u32::from(start(b) == start(a).min(start(b)))),
        right_cnt: u32::from(a.end == a.end.max(b.end))
            .saturating_add(u32::from(b.end == a.end.max(b.end))),
    };
    let first = extent.first_tlen(start(a), a.end, a.flags.is_first_in_template());
    let tlens = [first, first.saturating_neg()];
    [(a, b), (b, a)].into_iter().zip(tlens).all(|((r, mate), tlen)| {
        let mut flags = r.flags.with(0x1);
        if mate.flags.is_unmapped() {
            flags = flags.with(0x8);
        }
        if mate.flags.is_reverse() {
            flags = flags.with(0x20);
        }
        flags == r.flags
            && r.mate.tid == mate.tid
            && r.mate.pos == mate.pos
            && i64::from(r.mate.template_len.map_or(0, i32::from)) == tlen
    })
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::arithmetic_side_effects, reason = "test code")]
mod tests {
    use super::*;
    use crate::bam::test_util::tid;
    use std::num::NonZeroI32;

    /// A mapped 100-base read named `name` at `pos`, whose mate is at
    /// `mate_pos` with TLEN `tlen`.
    fn read(
        names: &mut Vec<u8>,
        name: &[u8],
        flags: u16,
        pos: u32,
        mate_pos: u32,
        tlen: i32,
    ) -> MateRecord {
        let start = names.len();
        names.extend_from_slice(name);
        let range = start..names.len();
        names.push(0);
        MateRecord {
            name: range,
            flags: BamFlags::from(flags),
            tid: Some(tid(0)),
            pos: Pos0::new(pos),
            end: i64::from(pos) + 100,
            mate: MateFields {
                tid: Some(tid(0)),
                pos: Pos0::new(mate_pos),
                template_len: NonZeroI32::new(tlen),
            },
            seq_unknown: false,
        }
    }

    // r[verify cram_writer.mates+2]
    /// A pair whose fields the decode rebuilds is attached: the first record
    /// downstream, nothing detached. A wrong PNEXT or a TLEN the decode
    /// would compute otherwise keeps both detached.
    #[test]
    fn pairs_are_attached_only_where_the_decode_rebuilds_them() {
        let mut names = Vec::new();
        // 99/147: read1 forward at 100, mate reverse at 400, ending at 500.
        let a = read(&mut names, b"p", 99, 100, 400, 400);
        let other = read(&mut names, b"x", 0, 200, 0, 0);
        let b = read(&mut names, b"p", 147, 400, 100, -400);
        let s = series(&[a.clone(), other.clone(), b.clone()], &names);
        assert_eq!(s.cf, [0x1 | 0x4, 0x1 | 0x2, 0x1]);
        assert_eq!(s.nf, [1]);
        assert_eq!(s.ts, [0], "only the detached record's TLEN");

        let detached = [0x1 | 0x2; 3];
        let mut wrong_pos = b.clone();
        wrong_pos.mate.pos = Pos0::new(101);
        let s = series(&[a.clone(), other.clone(), wrong_pos], &names);
        assert_eq!(s.cf, detached);
        assert!(s.nf.is_empty());

        let mut odd_tlen = b;
        odd_tlen.mate.template_len = NonZeroI32::new(-7);
        assert_eq!(series(&[a, other, odd_tlen], &names).cf, detached);
    }
}
