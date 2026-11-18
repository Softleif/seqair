//! Mates attached within a slice: a pair is stored as a link from the first
//! record to the second, from which the decode rebuilds both records' mate
//! fields, where storing the fields themselves (detached) costs a few bytes
//! a record. A pair is attached only where the rebuilt fields are the
//! record's own.

use crate::bam::mate::MateFields;
use crate::cram::varint::write_itf8;
use crate::reader::Tid;
use rustc_hash::FxHashMap;
use seqair_types::{BamFlags, Pos0};
use std::num::NonZeroI32;
use std::ops::Range;

/// A record as the decode's mate resolution sees it, and the mate fields it
/// must give back.
#[derive(Debug, Clone)]
pub(super) struct MateRecord {
    /// The record's name in the slice's `RN` block.
    pub name: Range<usize>,
    pub flags: BamFlags,
    pub tid: Option<Tid>,
    /// 0-based position as the decode has it: `AP - 1`, so -1 unplaced.
    pub pos: i64,
    /// 0-based exclusive end: the reference span, at least one base.
    pub end: i64,
    pub mate: MateFields,
    /// The CRAM flags other than the mate bits: `0x1`, and `0x8` for SEQ `*`.
    pub cram_flags: i32,
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

const DETACHED: i32 = 0x2;
const MATE_DOWNSTREAM: i32 = 0x4;
const EXPLICIT_TLEN: i32 = 0x10;

/// How one record of an attached pair is stored.
#[derive(Debug, Clone, Copy)]
struct Attached {
    /// The record's mate, for the first of a pair.
    downstream: Option<usize>,
    /// TLEN, where the decode would compute another.
    explicit_tlen: Option<i32>,
}

// r[impl cram_writer.mates]
/// `CF` and the mate series of `records`, whose names are in `names` (the
/// `RN` block).
pub(super) fn series(records: &[MateRecord], names: &[u8]) -> MateSeries {
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
    let mut attached: Vec<Option<Attached>> = vec![None; records.len()];
    for &(a, b, many) in by_name.values() {
        let (Some(b), false) = (b, many) else { continue };
        if let (Some(ra), Some(rb)) = (records.get(a), records.get(b))
            && let Some((tlen_a, tlen_b)) = rebuilt(ra, rb)
        {
            // The decode keeps a TLEN that does not fit as 0.
            let explicit = |want: Option<NonZeroI32>, got: i64| {
                let want = want.map_or(0, i32::from);
                (want != i32::try_from(got).unwrap_or(0)).then_some(want)
            };
            let tlen_a = explicit(ra.mate.template_len, tlen_a);
            let tlen_b = explicit(rb.mate.template_len, tlen_b);
            if let Some(slot) = attached.get_mut(a) {
                *slot = Some(Attached { downstream: Some(b), explicit_tlen: tlen_a });
            }
            if let Some(slot) = attached.get_mut(b) {
                *slot = Some(Attached { downstream: None, explicit_tlen: tlen_b });
            }
        }
    }

    let mut out = MateSeries::default();
    for (i, (r, link)) in records.iter().zip(&attached).enumerate() {
        match link {
            Some(Attached { downstream, explicit_tlen }) => {
                let mut cf = r.cram_flags;
                if explicit_tlen.is_some() {
                    cf |= EXPLICIT_TLEN;
                }
                if let Some(b) = downstream {
                    cf |= MATE_DOWNSTREAM;
                    write_itf8(&mut out.cf, cf);
                    let between = b.saturating_sub(i).saturating_sub(1);
                    write_itf8(&mut out.nf, i32::try_from(between).unwrap_or(i32::MAX));
                } else {
                    write_itf8(&mut out.cf, cf);
                }
                if let Some(tlen) = explicit_tlen {
                    write_itf8(&mut out.ts, *tlen);
                }
            }
            None => {
                write_itf8(&mut out.cf, r.cram_flags | DETACHED);
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

/// The TLENs the decode computes for `a` and its downstream mate `b`, if it
/// rebuilds both records' FLAG, RNEXT and PNEXT as they are; `None` if not.
/// This is the decode's `resolve_mates` and `close_chain` for a ring of two.
fn rebuilt(a: &MateRecord, b: &MateRecord) -> Option<(i64, i64)> {
    for (r, mate) in [(a, b), (b, a)] {
        let mut flags = r.flags.with(0x1);
        if mate.flags.is_unmapped() {
            flags = flags.with(0x8);
        }
        if mate.flags.is_reverse() {
            flags = flags.with(0x20);
        }
        if flags != r.flags || r.mate.tid != mate.tid || r.mate.pos != Pos0::try_from(mate.pos).ok()
        {
            return None;
        }
    }

    let aleft = a.pos.min(b.pos);
    let aright = a.end.max(b.end);
    let left_cnt = [a.pos, b.pos].iter().filter(|&&p| p == aleft).count();
    let right_cnt = [a.end, b.end].iter().filter(|&&e| e == aright).count();
    let span = aright.saturating_sub(aleft);
    let first = if a.tid != b.tid {
        0
    } else if a.pos == aleft && (a.end < aright || left_cnt <= 1) {
        span
    } else if a.pos == aleft && a.end == aright && left_cnt > 1 && right_cnt > 1 {
        if a.flags.is_first_in_template() { span } else { span.saturating_neg() }
    } else {
        span.saturating_neg()
    };
    let unmapped = a.flags.is_unmapped() || b.flags.is_unmapped();
    let tlen = |t: i64| if unmapped { 0 } else { t };
    Some((tlen(first), tlen(first.saturating_neg())))
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::arithmetic_side_effects, reason = "test code")]
mod tests {
    use super::*;
    use crate::bam::test_util::tid;

    /// A paired, mapped 100-base read named `name` at `pos`, whose mate is
    /// at `mate_pos` with TLEN `tlen`.
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
            pos: i64::from(pos),
            end: i64::from(pos) + 100,
            mate: MateFields {
                tid: Some(tid(0)),
                pos: Pos0::new(mate_pos),
                template_len: NonZeroI32::new(tlen),
            },
            cram_flags: 0x1,
        }
    }

    // r[verify cram_writer.mates]
    /// A pair whose fields the decode rebuilds is attached: the first record
    /// downstream, nothing detached; a wrong PNEXT keeps both detached; a
    /// TLEN the decode would compute otherwise is stored explicitly.
    #[test]
    fn pairs_are_attached_only_where_the_decode_rebuilds_them() {
        let mut names = Vec::new();
        // 99/147: read1 forward, mate reverse; TLEN 300 - 100 + ... = 400.
        let a = read(&mut names, b"p", 99, 100, 400, 400);
        let other = read(&mut names, b"x", 0, 200, 0, 0);
        let b = read(&mut names, b"p", 147, 400, 100, -400);
        let s = series(&[a.clone(), other.clone(), b.clone()], &names);
        assert_eq!(s.cf, [0x1 | 0x4, 0x1 | 0x2, 0x1]);
        assert_eq!(s.nf, [1]);
        assert_eq!(s.ts, [0], "only the detached record's TLEN");

        let mut wrong_pos = b.clone();
        wrong_pos.mate.pos = Pos0::new(101);
        let s = series(&[a.clone(), other.clone(), wrong_pos], &names);
        assert_eq!(s.cf, [0x1 | 0x2, 0x1 | 0x2, 0x1 | 0x2]);
        assert!(s.nf.is_empty());

        let mut odd_tlen = b;
        odd_tlen.mate.template_len = NonZeroI32::new(-7);
        let s = series(&[a, other, odd_tlen], &names);
        assert_eq!(s.cf, [0x1 | 0x4, 0x1 | 0x2, 0x1 | 0x10]);
        let mut ts = Vec::new();
        write_itf8(&mut ts, 0);
        write_itf8(&mut ts, -7);
        assert_eq!(s.ts, ts);
    }
}
