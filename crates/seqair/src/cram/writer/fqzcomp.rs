//! fqzcomp quality encoder (CRAM 3.1, method 7): the inverse of
//! [`crate::cram::fqzcomp`], as htscodecs' `fqz_compress`
//! (`fqzcomp_qual.c`) writes it for CRAM 3.x with one of its built-in
//! strategies.
//!
//! Parameters are picked as `fqz_pick_parameters` and `fqz_qual_stats` pick
//! them, in the same floating-point operations (so up to rounding near a
//! threshold, where a C build with fused multiply-adds can decide the other
//! way); the stream is coded as `compress_block_fqz2f` codes it. The context of each quality comes from
//! the decoder's own parameter block, parsed back from the bytes written, so
//! the two cannot disagree on it.

use super::rans_nx16::put_uint7;
use crate::cram::fqzcomp::{
    GFLAG_DO_REV, GFLAG_HAVE_STAB, PFLAG_DO_DEDUP, PFLAG_DO_SEL, PFLAG_FIXED_LEN, PFLAG_HAVE_DTAB,
    PFLAG_HAVE_PTAB, PFLAG_HAVE_QMAP, Params, QualModels, QualState, VERSION,
};
use crate::cram::range_coder::{AdaptiveModel, RangeEncoder};

// r[impl cram_writer.fqzcomp]

/// One record's share of a `QS` block.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct QualRecord {
    /// How many of the block's bytes are the record's.
    pub len: u32,
    /// BAM flag `0x10`: the qualities are coded reversed.
    pub reverse: bool,
    /// BAM flag `0x80`: a strategy may model second reads apart.
    pub read2: bool,
}

/// htscodecs' built-in strategies, as htslib's methods name them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Strategy {
    /// `FQZ`, htscodecs' basic options; htslib tries it from level 1.
    Fqz,
    /// `FQZ_b`, tuned on `HiSeq` 2000 data; from level 5.
    FqzB,
    /// `FQZ_c`, tuned on `MiSeq` data; from level 7.
    FqzC,
    /// `FQZ_d`, an adaptive order 1 for Ion Torrent; from level 7.
    FqzD,
}

impl Strategy {
    #[allow(dead_code, reason = "the writer's codec choice plumbs the strategies in")]
    pub(super) const ALL: [Self; 4] = [Self::Fqz, Self::FqzB, Self::FqzC, Self::FqzD];

    /// htscodecs' `strat_opts` row: qbits, qshift, pbits, pshift, dbits,
    /// dshift, qloc, sloc, ploc, dloc, `do_r2`, `do_qa`.
    const fn opts(self) -> [i32; 12] {
        match self {
            Self::Fqz => [10, 5, 4, -1, 2, 1, 0, 14, 10, 14, 0, -1],
            Self::FqzB => [8, 5, 7, 0, 0, 0, 0, 14, 8, 14, 1, -1],
            Self::FqzC => [12, 6, 2, 0, 2, 3, 0, 9, 12, 14, 0, 0],
            Self::FqzD => [12, 6, 0, 0, 0, 0, 0, 12, 0, 0, 0, 0],
        }
    }
}

/// `qualities`, the concatenated `QS` bytes of `records`, as an fqzcomp
/// stream picked and coded with `strategy`. `None` when the records' lengths
/// do not add up to the qualities, or the block is too large for the
/// format's sizes.
#[allow(dead_code, reason = "the writer's codec choice plumbs the encoder in")]
pub(super) fn encode(
    qualities: &[u8],
    records: &[QualRecord],
    strategy: Strategy,
) -> Option<Vec<u8>> {
    // An empty record cannot be coded (the decoder rejects length 0), and
    // htscodecs refuses one before the last non-empty record: those are left
    // out. Empty records after it are never coded, but count in htscodecs'
    // statistics (record count, fixed length, selectors), so they stay.
    let coded = records.iter().rposition(|r| r.len > 0).map_or(0, |i| i.saturating_add(1));
    let (body, tail) = records.split_at(coded);
    let records: Vec<QualRecord> =
        body.iter().copied().filter(|r| r.len > 0).chain(tail.iter().copied()).collect();
    let total = records.iter().try_fold(0usize, |sum, r| sum.checked_add(r.len as usize))?;
    if total != qualities.len() || i32::try_from(total).is_err() {
        return None;
    }

    let picked = pick_parameters(qualities, &records, strategy)?;
    let mut out = Vec::with_capacity(qualities.len() / 2);
    put_uint7(&mut out, qualities.len())?;
    let params_start = out.len();
    picked.store(&mut out)?;
    // The decoder wants ten bytes of parameters and stream, as htscodecs
    // does; the range coder's five final bytes always follow these.
    let mut written = out.get(params_start..)?.to_vec();
    written.extend_from_slice(&[0; 5]);
    let params = Params::parse(&mut written.as_slice()).ok()?;
    let first = params.first()?;

    // CRAM 3.x keeps qualities in alignment orientation; fqzcomp models
    // reverse-strand records in the order they were sequenced.
    let mut buf = qualities.to_vec();
    let mut at = 0usize;
    for r in &records {
        let end = at.checked_add(r.len as usize)?;
        if r.reverse {
            buf.get_mut(at..end)?.reverse();
        }
        at = end;
    }

    let mut rc = RangeEncoder::new(out);
    let mut qual = QualModels::new(u16::from(picked.max_sym).checked_add(1)?);
    let mut len_models: [AdaptiveModel<256>; 4] = std::array::from_fn(|_| AdaptiveModel::new(256));
    let mut rev = AdaptiveModel::<2>::new(2);
    let mut dup = AdaptiveModel::<2>::new(2);
    let mut sel_model = AdaptiveModel::<256>::new(usize::from(picked.max_sel).saturating_add(1));

    let mut first_len = true;
    let mut last_len: Option<u32> = None;
    let mut at = 0usize;
    for (r, &sel) in records.iter().zip(&picked.selectors) {
        if r.len == 0 {
            // The trailing empty records: nothing is left to code.
            break;
        }
        if picked.do_sel {
            sel_model.encode(&mut rc, u16::try_from(sel).ok()?)?;
        }
        if !picked.fixed_len || first_len {
            for (model, byte) in len_models.iter_mut().zip(r.len.to_le_bytes()) {
                model.encode(&mut rc, u16::from(byte))?;
            }
            first_len = false;
        }
        rev.encode(&mut rc, u16::from(r.reverse))?;

        let end = at.checked_add(r.len as usize)?;
        let record = buf.get(at..end)?;
        if picked.do_dedup {
            let previous = last_len.and_then(|l| buf.get(at.checked_sub(l as usize)?..at));
            if previous == Some(record) {
                dup.encode(&mut rc, 1)?;
                at = end;
                continue;
            }
            dup.encode(&mut rc, 0)?;
            last_len = Some(r.len);
        }

        let mut st = QualState::new(r.len);
        let mut ctx = picked.context;
        for &q in record {
            let sym = picked.qmap.get(usize::from(q)).copied().flatten()?;
            qual.encode(ctx, &mut rc, u16::from(sym))?;
            ctx = first.update(&mut st, sym, sel);
        }
        at = end;
    }
    Some(rc.finish())
}

/// The parameters htscodecs picks, and what the coding loop needs of them.
#[derive(Debug)]
struct Picked {
    gflags: u8,
    max_sel: u8,
    context: u16,
    pflags: u8,
    max_sym: u8,
    /// `qbits << 4 | qshift`, `qloc << 4 | sloc`, `ploc << 4 | dloc`.
    nibbles: [u8; 3],
    /// Quality to model symbol; `None` for a quality the block lacks.
    qmap: [Option<u8>; 256],
    store_qmap: bool,
    ptab: Option<[u32; PTAB_SIZE]>,
    dtab: Option<[u32; DTAB_SIZE]>,
    selectors: Vec<u32>,
    do_sel: bool,
    fixed_len: bool,
    do_dedup: bool,
}

const PTAB_SIZE: usize = 1024;
const DTAB_SIZE: usize = 256;

impl Picked {
    /// htscodecs' `fqz_store_parameters`.
    fn store(&self, out: &mut Vec<u8>) -> Option<()> {
        out.extend_from_slice(&[VERSION, self.gflags]);
        if self.gflags & GFLAG_HAVE_STAB != 0 {
            out.push(self.max_sel);
            // Every selector maps to the one parameter block.
            store_array(out, &[0; 256])?;
        }
        out.extend_from_slice(&self.context.to_le_bytes());
        out.extend_from_slice(&[self.pflags, self.max_sym]);
        out.extend_from_slice(&self.nibbles);
        if self.store_qmap {
            out.extend((0..=255u8).zip(&self.qmap).filter(|(_, m)| m.is_some()).map(|(q, _)| q));
        }
        if let Some(ptab) = &self.ptab {
            store_array(out, ptab)?;
        }
        if let Some(dtab) = &self.dtab {
            store_array(out, dtab)?;
        }
        Some(())
    }
}

/// htscodecs' `fqz_param` fields that steer the picking, as its `int`s.
#[derive(Debug, Clone, Copy)]
struct Shape {
    qbits: i32,
    qshift: i32,
    pbits: i32,
    pshift: i32,
    dbits: i32,
    dshift: i32,
    qloc: i32,
    sloc: i32,
    ploc: i32,
    dloc: i32,
    do_r2: i32,
    do_qa: i32,
}

impl Shape {
    const fn new(o: [i32; 12]) -> Self {
        let [qbits, qshift, pbits, pshift, dbits, dshift, qloc, sloc, ploc, dloc, do_r2, do_qa] = o;
        Self { qbits, qshift, pbits, pshift, dbits, dshift, qloc, sloc, ploc, dloc, do_r2, do_qa }
    }
}

/// htscodecs' `fqz_pick_parameters` for one of its strategies, CRAM 3.x.
#[allow(
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "htscodecs' int arithmetic on nibble-sized fields"
)]
fn pick_parameters(quals: &[u8], records: &[QualRecord], strategy: Strategy) -> Option<Picked> {
    // Approximately sqrt(delta).
    let mut dsqr: [i32; 64] = [
        0, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5,
        5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        7, 7, 7, 7,
    ];
    let mut pm = Shape::new(strategy.opts());
    let stats = qual_stats(quals, records, &mut pm);

    let store_qmap = stats.nsym <= 8 && stats.nsym * 2 < stats.max_sym;
    let first_len = records.first().map_or(0, |r| r.len);
    let fixed_len = records.iter().all(|r| r.len == first_len);

    if pm.pshift < 0 {
        let p = (f64::from(first_len) / f64::from(1u32 << pm.pbits)).ln() / 2f64.ln() + 0.5;
        // C's `MAX(0, p)` then the conversion to int, which truncates.
        pm.pshift = if 0.0 > p { 0 } else { p as i32 };
    }
    if stats.nsym <= 4 {
        // NovaSeq
        pm.qshift = 2;
        if quals.len() < 5_000_000 {
            pm.pbits = 2;
            pm.pshift = 5;
        }
    } else if stats.nsym <= 8 {
        // HiSeq X
        pm.qbits = pm.qbits.min(9);
        pm.qshift = 3;
        if quals.len() < 5_000_000 {
            pm.qbits = 6;
        }
    }
    if quals.len() < 300_000 {
        pm.qbits = pm.qshift;
        pm.dbits = 2;
    }

    let dmax = (1i32 << pm.dbits) - 1;
    for d in &mut dsqr {
        *d = (*d).min(dmax);
    }

    let mut qmap = [None; 256];
    let max_sym = if store_qmap {
        let mut next = 0u8;
        for (m, &count) in qmap.iter_mut().zip(&stats.qhist) {
            if count > 0 {
                *m = Some(next);
                next += 1;
            }
        }
        stats.nsym
    } else {
        for (q, m) in (0..=255u8).zip(&mut qmap) {
            *m = Some(q);
        }
        stats.max_sym
    };

    let ptab = (pm.pbits > 0).then(|| {
        let cap = (1u32 << pm.pbits) - 1;
        std::array::from_fn(|i| cap.min((i as u32).checked_shr(pm.pshift as u32).unwrap_or(0)))
    });
    let dtab = (pm.dbits > 0).then(|| {
        // Deltas past the table share its last entry.
        let last = dsqr.len() - 1;
        std::array::from_fn(|i| dsqr.get((i >> pm.dshift).min(last)).copied().unwrap_or(0) as u32)
    });

    // The selectors stats picked, or none if none came out above zero (an
    // htscodecs stream would code them with an uninitialised model).
    let max_sel = stats.selectors.iter().copied().max().unwrap_or(0);
    let do_sel = stats.do_sel && max_sel > 0;

    let flag = |on: bool, f: u8| if on { f } else { 0 };
    let pflags = flag(dtab.is_some(), PFLAG_HAVE_DTAB)
        | flag(ptab.is_some(), PFLAG_HAVE_PTAB)
        | flag(do_sel, PFLAG_DO_SEL)
        | flag(fixed_len, PFLAG_FIXED_LEN)
        | flag(stats.do_dedup, PFLAG_DO_DEDUP)
        | flag(store_qmap, PFLAG_HAVE_QMAP);

    let nibble = |hi: i32, lo: i32| -> Option<u8> {
        let hi = u8::try_from(hi).ok().filter(|&h| h < 16)?;
        let lo = u8::try_from(lo).ok().filter(|&l| l < 16)?;
        Some(hi << 4 | lo)
    };
    Some(Picked {
        gflags: GFLAG_DO_REV | flag(do_sel, GFLAG_HAVE_STAB),
        max_sel: u8::try_from(if do_sel { max_sel } else { 0 }).ok()?,
        context: 0,
        pflags,
        max_sym: u8::try_from(max_sym).ok()?,
        nibbles: [
            nibble(pm.qbits, pm.qshift)?,
            nibble(pm.qloc, pm.sloc)?,
            nibble(pm.ploc, pm.dloc)?,
        ],
        qmap,
        store_qmap,
        ptab,
        dtab,
        selectors: if do_sel { stats.selectors } else { vec![0; records.len()] },
        do_sel,
        fixed_len,
        do_dedup: stats.do_dedup,
    })
}

/// What `fqz_qual_stats` learns of a block.
#[derive(Debug)]
struct Stats {
    qhist: [u32; 256],
    /// Distinct qualities, and the largest.
    nsym: u32,
    max_sym: u32,
    do_dedup: bool,
    do_sel: bool,
    /// Per record, htscodecs' `flags >> 16`.
    selectors: Vec<u32>,
}

/// Histogram buckets by position from the record's end.
const NP: usize = 32;
/// Average-quality buckets, in tenths.
const AVG: usize = 2560;

/// htscodecs' `fast_log`: a linear fit to the logarithm off a double's bits.
#[allow(clippy::cast_precision_loss, reason = "htscodecs converts the same way")]
fn fast_log(a: f64) -> f64 {
    let x = i64::from_ne_bytes(a.to_ne_bytes());
    x.wrapping_sub(4_606_921_278_410_026_770) as f64 * 1.539_095_918_623_324e-16
}

/// htscodecs' `fqz_qual_stats` over every record (`one_param` -1): the
/// quality histogram, whether to dedup, and the average-quality and
/// READ1/READ2 selectors where they pay, which take bits from `pm`.
#[allow(
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss,
    clippy::cast_precision_loss,
    clippy::indexing_slicing,
    clippy::too_many_lines,
    reason = "a port of htscodecs' C: counts below the block size (< 2^31), \
              qualities (u8) index 256-entry tables, j % 32 the 32 buckets, \
              averages are clamped below 2560, and selectors are 0..=3"
)]
fn qual_stats(quals: &[u8], records: &[QualRecord], pm: &mut Shape) -> Stats {
    let n = records.len();
    let mut qhist = [0u32; 256];
    let mut qhistb = vec![[0u32; 256]; NP];
    let mut qhist1 = vec![[0u32; 256]; NP];
    let mut qhist2 = vec![[0u32; 256]; NP];
    let mut t1 = [0u64; NP];
    let mut t2 = [0u64; NP];
    let mut avg = vec![0u32; AVG];
    let mut avg_qual = vec![0u32; n];
    let has_r2 = records.iter().any(|r| r.read2);

    let mut dedup = 0usize;
    let mut last_len = 0usize;
    let mut i = 0usize;
    // htscodecs stops at the end of the qualities: trailing empty records
    // count as records, but get no average and no dedup count.
    let mut seen = 0usize;
    for (rec, r) in records.iter().enumerate() {
        let len = r.len as usize;
        if len == 0 {
            break;
        }
        seen += 1;
        debug_assert!(i >= last_len && i + len <= quals.len(), "records add up to the qualities");
        if i > 0 && len == last_len && quals[i - last_len..i] == quals[i..i + len] {
            dedup += 1;
        }
        last_len = len;
        let (qh, th) = if r.read2 { (&mut qhist2, &mut t2) } else { (&mut qhist1, &mut t1) };
        let mut tot = 0u32;
        for (k, &q) in quals[i..i + len].iter().enumerate() {
            let j = (len - k) & (NP - 1);
            let q = usize::from(q);
            tot = tot.wrapping_add(q as u32);
            qhist[q] += 1;
            qhistb[j][q] += 1;
            qh[j][q] += 1;
            th[j] += 1;
        }
        i += len;
        let tot = (f64::from(tot) * 10.0 / len as f64 + 0.5) as u32;
        avg_qual[rec] = tot;
        avg[(tot as usize).min(AVG - 1)] += 1;
    }
    let do_dedup = (seen + 1) / (dedup + 1) < 500;

    let nsym = qhist.iter().filter(|&&c| c > 0).count() as u32;
    let max_sym = qhist.iter().rposition(|&c| c > 0).unwrap_or(0) as u32;

    let mut selectors = vec![0u32; n];
    let mut max_sel = 0u32;
    let mut do_sel = false;
    let nr = n as i32;

    if pm.do_qa != 0 {
        // Turn avg[] from counts into each average's bucket, 0..=3, so the
        // buckets split the records at these fractions.
        let qf = if nsym > 8 { [0.2, 0.5, 0.8] } else { [0.05, 0.22, 0.60] };
        let mut total = 0i64;
        let mut a = 0usize;
        for (bucket, &f) in qf.iter().enumerate() {
            while a < AVG {
                total += i64::from(avg[a]);
                if total as f64 > f * f64::from(nr) {
                    break;
                }
                avg[a] = bucket as u32;
                a += 1;
            }
        }
        for v in avg.iter_mut().skip(a) {
            *v = 3;
        }

        // Entropy of the qualities by position, unsplit, split in two and
        // in four by average quality, over every eighth record.
        let mut qbin4 = vec![[[0i32; 256]; NP]; 4];
        let mut qbin2 = vec![[[0i32; 256]; NP]; 2];
        let mut qbin1 = vec![[0i32; 256]; NP];
        let mut qcnt4 = [[0i32; NP]; 4];
        let mut qcnt2 = [[0i32; NP]; 2];
        let mut qcnt1 = [0i32; NP];
        let mut i = 0usize;
        for (rec, r) in records.iter().enumerate() {
            let len = r.len as usize;
            if rec & 7 == 0 {
                let qb4 = avg[(avg_qual[rec] as usize).min(AVG - 1)] as usize;
                let qb2 = qb4 / 2;
                debug_assert!(qb4 < 4 && i + len <= quals.len(), "a bucket per quarter");
                for (k, &q) in quals[i..i + len].iter().enumerate() {
                    let x = (len - k) & (NP - 1);
                    let q = usize::from(q);
                    qbin4[qb4][x][q] += 1;
                    qcnt4[qb4][x] += 1;
                    qbin2[qb2][x][q] += 1;
                    qcnt2[qb2][x] += 1;
                    qbin1[x][q] += 1;
                    qcnt1[x] += 1;
                }
            }
            i += len;
        }

        let term = |count: i32, total: i32| {
            f64::from(count) * fast_log(f64::from(count) / f64::from(total))
        };
        let (mut e1, mut e2, mut e4) = (0.0f64, 0.0f64, 0.0f64);
        for j in 0..NP {
            for q in 0..256 {
                if qbin1[j][q] != 0 {
                    e1 += term(qbin1[j][q], qcnt1[j]);
                }
                for b in 0..2 {
                    if qbin2[b][j][q] != 0 {
                        e2 += term(qbin2[b][j][q], qcnt2[b][j]);
                    }
                }
                for b in 0..4 {
                    if qbin4[b][j][q] != 0 {
                        e4 += term(qbin4[b][j][q], qcnt4[b][j]);
                    }
                }
            }
        }
        let bytes = -2f64.ln() / 8.0;
        e1 /= bytes;
        e2 /= bytes;
        e4 /= bytes;

        // A selector takes bits from the rest of the context: htscodecs
        // asks it to save 2% (a guess, its comment says).
        let qm = if pm.do_qa > 0 { 1.0 } else { 0.98 };
        let (n4, n8) = (f64::from(nr / 4), f64::from(nr / 8));
        if (pm.do_qa == -1 || pm.do_qa >= 4) && e4 + n4 < e2 * qm + n8 && e4 + n4 < e1 * qm {
            for (s, &a) in selectors.iter_mut().zip(&avg_qual) {
                *s |= avg[(a as usize).min(AVG - 1)];
            }
            do_sel = true;
            max_sel = 3;
        } else if (pm.do_qa == -1 || pm.do_qa >= 2) && e2 + n8 < e1 * qm {
            for (s, &a) in selectors.iter_mut().zip(&avg_qual) {
                *s |= avg[(a as usize).min(AVG - 1)] >> 1;
            }
            do_sel = true;
            max_sel = 1;
        }

        if pm.do_qa == -1 {
            // Make room for the selector: qual, pos, delta in that order.
            if pm.pbits > 0 && pm.dbits > 0 {
                pm.sloc = pm.dloc - 1;
                pm.pbits -= 1;
                pm.dbits -= 1;
                pm.dloc += 1;
            } else if pm.dbits >= 2 {
                pm.sloc = pm.dloc;
                pm.dbits -= 2;
                pm.dloc += 2;
            } else if pm.qbits >= 2 {
                pm.qbits -= 2;
                pm.ploc -= 2;
                pm.sloc = 16 - 2 - pm.do_r2;
                if pm.qbits == 6 && pm.qshift == 5 {
                    pm.qbits -= 1;
                }
            }
            pm.do_qa = 4;
        }
    }

    // Does modelling READ1 and READ2 apart pay?
    if has_r2 || pm.do_r2 != 0 {
        let (mut e1, mut e2) = (0.0f64, 0.0f64);
        for j in 0..NP {
            if t1[j] == 0 || t2[j] == 0 {
                continue;
            }
            for q in 0..256 {
                let b = qhistb[j][q];
                if b == 0 {
                    continue;
                }
                e1 -= f64::from(b) * (f64::from(b) / (t1[j] + t2[j]) as f64).ln();
                if qhist1[j][q] != 0 {
                    let h = f64::from(qhist1[j][q]);
                    e2 -= h * (h / t1[j] as f64).ln();
                }
                if qhist2[j][q] != 0 {
                    let h = f64::from(qhist2[j][q]);
                    e2 -= h * (h / t2[j] as f64).ln();
                }
            }
        }
        e1 /= 2f64.ln() * 8.0;
        e2 /= 2f64.ln() * 8.0;
        let qm = if pm.do_r2 > 0 { 1.0 } else { 0.95 };
        if e2 + f64::from(8 + nr / 8) < e1 * qm {
            for (s, r) in selectors.iter_mut().zip(records) {
                *s = *s * 2 + u32::from(r.read2);
                max_sel = max_sel.max(*s);
            }
        }
    }

    if max_sel > 0 {
        do_sel = true;
    }
    Stats { qhist, nsym, max_sym, do_dedup, do_sel, selectors }
}

/// htscodecs' `store_array`: the run length of each value in turn (runs of
/// 255 or more continue in the next byte), then those bytes run-length coded
/// by following a repeated byte with the count of further copies. `array`
/// must be non-decreasing from 0.
#[allow(
    clippy::arithmetic_side_effects,
    reason = "run lengths and values are bounded by the array's length"
)]
fn store_array(out: &mut Vec<u8>, array: &[u32]) -> Option<()> {
    // A value below one already stored would never end its run.
    let sorted = array.is_sorted();
    debug_assert!(sorted, "fqzcomp tables are non-decreasing");
    if !sorted {
        return None;
    }
    let mut runs = Vec::new();
    let mut i = 0usize;
    let mut value = 0u32;
    while i < array.len() {
        let mut run = array.get(i..)?.iter().take_while(|&&v| v == value).count();
        i += run;
        loop {
            let part = run.min(255);
            runs.push(u8::try_from(part).ok()?);
            run -= part;
            if part != 255 {
                break;
            }
        }
        value += 1;
    }

    let mut last = None;
    let mut rest = runs.as_slice();
    while let Some((&r, after)) = rest.split_first() {
        out.push(r);
        rest = after;
        if last == Some(r) {
            let copies = rest.iter().take_while(|&&x| x == r).count();
            out.push(u8::try_from(copies).ok()?);
            rest = rest.get(copies..)?;
        } else {
            last = Some(r);
        }
    }
    Some(())
}

#[cfg(test)]
#[allow(
    clippy::unwrap_used,
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    clippy::indexing_slicing,
    reason = "test code"
)]
mod tests {
    use super::*;
    use crate::cram::fqzcomp::decode;
    use crate::cram::fqzcomp_reference::htscodecs_fqz::{
        FREAD2, FREVERSE, Records, arb_records, htscodecs_compress,
    };
    use hegel::prelude::*;

    // r[verify cram_writer.fqzcomp]

    fn qual_records(records: &Records) -> Vec<QualRecord> {
        records
            .lens
            .iter()
            .zip(&records.flags)
            .map(|(&len, &f)| QualRecord {
                len,
                reverse: f & FREVERSE != 0,
                read2: f & FREAD2 != 0,
            })
            .collect()
    }

    fn round_trip(records: &Records, strategy: Strategy) -> Vec<u8> {
        let stream = encode(&records.quals, &qual_records(records), strategy).unwrap();
        let decoded = decode(&stream).unwrap_or_else(|e| panic!("{strategy:?}: {e}"));
        assert!(decoded == records.quals, "{strategy:?}: qualities differ");
        stream
    }

    #[hegel::composite]
    fn arb_block(tc: &TestCase) -> Records {
        let max_alphabet = tc.draw_silent(gs::sampled_from(&[4, 8, 64, 256]));
        tc.draw(arb_records(max_alphabet))
    }

    /// htscodecs' `fqz_compress`, linked through hts-sys, writes the same
    /// bytes for every strategy (`vers` 3, as htslib passes for CRAM 3.1).
    #[hegel::test(test_cases = crate::pinned::cases(100))]
    fn matches_htscodecs_byte_for_byte(tc: TestCase) {
        let records = tc.draw(arb_block());
        for (strat, strategy) in (0..).zip(Strategy::ALL) {
            let ours = round_trip(&records, strategy);
            let theirs = htscodecs_compress(&records, 3, strat, None).unwrap();
            assert!(
                ours == theirs,
                "{strategy:?}: {} bytes against htscodecs' {}",
                ours.len(),
                theirs.len()
            );
        }
    }

    #[test]
    fn empty_records_are_left_out() {
        let quals = b"ABCDEFG".to_vec();
        let rec = |len| QualRecord { len, reverse: false, read2: false };
        let stream = encode(&quals, &[rec(0), rec(3), rec(0), rec(4)], Strategy::Fqz).unwrap();
        assert_eq!(decode(&stream).unwrap(), quals);
        let stream = encode(&[], &[rec(0)], Strategy::Fqz).unwrap();
        assert_eq!(decode(&stream).unwrap(), b"");
        let stream = encode(&[], &[], Strategy::FqzD).unwrap();
        assert_eq!(decode(&stream).unwrap(), b"");
    }

    /// htscodecs codes no empty record but takes trailing ones into its
    /// statistics: they clear the fixed length and can carry a selector.
    #[test]
    fn trailing_empty_records_match_htscodecs() {
        let quals: Vec<u8> = (0..60u8).map(|i| 30 + i % 7).collect();
        for (lens, flags) in [
            (vec![30, 30, 0], vec![0, 0, 0]),
            (vec![30, 30, 0, 0], vec![0, FREAD2, FREAD2 | FREVERSE, 0]),
            (vec![20, 40, 0], vec![FREVERSE, FREAD2, FREAD2]),
        ] {
            let records = Records { quals: quals.clone(), lens, flags, alphabet: Vec::new() };
            for (strat, strategy) in (0..).zip(Strategy::ALL) {
                let ours = round_trip(&records, strategy);
                let theirs = htscodecs_compress(&records, 3, strat, None).unwrap();
                assert!(ours == theirs, "{:?} {strategy:?}", records.lens);
            }
        }
    }

    #[test]
    fn lengths_must_add_up() {
        let rec = QualRecord { len: 3, reverse: false, read2: false };
        assert!(encode(b"ABCD", &[rec], Strategy::Fqz).is_none());
        assert!(encode(b"AB", &[rec], Strategy::Fqz).is_none());
    }

    /// A block of about `len` qualities over `alphabet` values: records of
    /// 100 to 250, every thirteenth a duplicate of the one before, whose
    /// average quality steps through `levels` levels, second reads `skew`
    /// lower.
    fn shaped_block(len: usize, alphabet: u8, levels: u32, skew: u8) -> Records {
        let mut r = Records {
            quals: Vec::new(),
            lens: Vec::new(),
            flags: Vec::new(),
            alphabet: Vec::new(),
        };
        let mut x = 0x2545_f491_u32;
        let mut next = move || {
            x ^= x << 13;
            x ^= x >> 17;
            x ^= x << 5;
            x
        };
        let mut k = 0u32;
        while r.quals.len() < len {
            if k % 13 == 12 {
                let (&l, &f) = (r.lens.last().unwrap(), r.flags.last().unwrap());
                let start = r.quals.len() - l as usize;
                r.quals.extend_from_within(start..);
                r.lens.push(l);
                r.flags.push(f);
                k += 1;
                continue;
            }
            let read2 = k % 2 == 1;
            let level = (k / 5 % levels) as u8;
            let base =
                (level * (alphabet / levels as u8)).saturating_sub(if read2 { skew } else { 0 });
            let l = 100 + next() % 151;
            for _ in 0..l {
                let q = (base + (next() % u32::from(alphabet / 4 + 1)) as u8).min(alphabet - 1);
                r.quals.push(30 + q);
            }
            r.lens.push(l);
            r.flags
                .push(if read2 { FREAD2 } else { 0 } | if next() & 64 != 0 { FREVERSE } else { 0 });
            k += 1;
        }
        r
    }

    /// Which of htscodecs' selectors paid: average-quality buckets (1, 2 or
    /// 4) and whether READ1 and READ2 were split, read off the selectors.
    fn selector_branch(records: &Records, strategy: Strategy) -> (u32, bool) {
        let recs = qual_records(records);
        let mut pm = Shape::new(strategy.opts());
        let stats = qual_stats(&records.quals, &recs, &mut pm);
        let split = stats.do_sel
            && stats.selectors.iter().zip(&recs).all(|(s, r)| s & 1 == u32::from(r.read2));
        let avg = stats.selectors.iter().map(|s| s >> u32::from(split)).max().unwrap_or(0) + 1;
        (avg, split)
    }

    /// Blocks past htscodecs' size thresholds (300 000 and 5 000 000
    /// qualities) with 4, 8 and 40 quality values, variable lengths and
    /// duplicates, shaped so every selector branch fires somewhere (none;
    /// two and four average-quality buckets, each with and without the
    /// READ1/READ2 split; the split alone): still htscodecs' bytes. The 5 M
    /// block, there only for its size threshold, runs one strategy.
    #[test]
    fn large_blocks_and_selectors_match_htscodecs() {
        // Per strategy (FQZ, FQZ_b, FQZ_c, FQZ_d): the average-quality
        // buckets and whether reads split. FQZ_c and FQZ_d never bucket.
        type Branches = [(u32, bool); 4];
        let cases: [(usize, u8, u32, u8, Branches); 5] = [
            (200_000, 4, 4, 0, [(1, false), (1, false), (1, false), (1, false)]),
            (320_000, 8, 4, 4, [(2, true), (2, true), (1, true), (1, true)]),
            (320_000, 40, 2, 0, [(4, false), (4, true), (1, false), (1, false)]),
            (320_000, 40, 1, 12, [(1, false), (1, false), (1, false), (1, false)]),
            (5_000_100, 8, 4, 4, [(2, true), (2, true), (1, true), (1, true)]),
        ];
        for (len, alphabet, levels, skew, branches) in cases {
            let records = shaped_block(len, alphabet, levels, skew);
            let strategies = if len > 5_000_000 { 1 } else { 4 };
            for ((strat, strategy), expected) in
                (0..).zip(Strategy::ALL).zip(branches).take(strategies)
            {
                let ours = round_trip(&records, strategy);
                let theirs = htscodecs_compress(&records, 3, strat, None).unwrap();
                assert!(ours == theirs, "{len} x {alphabet}, {strategy:?}");
                assert_eq!(
                    selector_branch(&records, strategy),
                    expected,
                    "{len} x {alphabet}, {levels} levels, skew {skew}, {strategy:?}"
                );
            }
        }
    }

    #[test]
    fn store_array_matches_the_spec_example() {
        // CRAMcodecs §6: A = {0,1,3,4,5,6,7,7,7,7} is stored as {1,1,0,0,1,1,2,4}.
        let mut out = Vec::new();
        store_array(&mut out, &[0, 1, 3, 4, 5, 6, 7, 7, 7, 7]).unwrap();
        assert_eq!(out, [1, 1, 0, 0, 1, 1, 2, 4]);
        out.clear();
        store_array(&mut out, &[0; 256]).unwrap();
        assert_eq!(out, [255, 1]);
    }
}
