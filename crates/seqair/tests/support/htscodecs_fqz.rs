//! htscodecs' fqzcomp encoder through FFI, and generators of records and
//! parameter sets for it — shared by `tests/htscodecs_fqzcomp_oracle.rs` and
//! the reference decoder's unit tests (`src/cram/fqzcomp_reference.rs`),
//! which include this file by path. It uses nothing from seqair.
//!
//! htscodecs is linked through `hts-sys` (the `rust-htslib` dev-dependency).
#![allow(
    dead_code,
    unreachable_pub,
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss,
    reason = "test support: each including crate uses part of it; values bounded by the generators"
)]

// Link hts-sys's static htslib, which carries htscodecs.
extern crate rust_htslib;

use std::ffi::{c_char, c_int, c_uint};

use hegel::TestCase;
use hegel::generators as gs;

/// `FQZ_FREVERSE`: the record's qualities are stored reversed (CRAM 3.1).
pub const FREVERSE: u32 = 16;
/// `FQZ_FREAD2`: second of a pair; the built-in strategies may split on it.
pub const FREAD2: u32 = 128;

// ── htscodecs FFI (fqzcomp_qual.h) ───────────────────────────────────

#[repr(C)]
pub struct FqzSlice {
    pub num_records: c_int,
    pub len: *mut u32,
    pub flags: *mut u32,
}

#[repr(C)]
#[derive(Clone)]
pub struct FqzParam {
    pub context: u16,
    pub pflags: c_uint,
    pub do_sel: c_uint,
    pub do_dedup: c_uint,
    pub store_qmap: c_uint,
    pub fixed_len: c_uint,
    pub use_qtab: u8,
    pub use_dtab: u8,
    pub use_ptab: u8,
    pub qbits: c_uint,
    pub qloc: c_uint,
    pub pbits: c_uint,
    pub ploc: c_uint,
    pub dbits: c_uint,
    pub dloc: c_uint,
    pub sbits: c_uint,
    pub sloc: c_uint,
    pub max_sym: c_int,
    pub nsym: c_int,
    pub max_sel: c_int,
    pub qmap: [c_uint; 256],
    pub qtab: [c_uint; 256],
    pub ptab: [c_uint; 1024],
    pub dtab: [c_uint; 256],
    pub qshift: c_int,
    pub pshift: c_int,
    pub dshift: c_int,
    pub sshift: c_int,
    pub qmask: c_uint,
    pub do_r2: c_int,
    pub do_qa: c_int,
}

#[repr(C)]
pub struct FqzGparams {
    pub vers: c_int,
    pub gflags: c_uint,
    pub nparam: c_int,
    pub max_sel: c_int,
    pub stab: [c_uint; 256],
    pub max_sym: c_int,
    pub p: *mut FqzParam,
}

// SAFETY: `p` points into a `Vec` owned by the `CustomParams` that hands the
// struct out by `&mut`, so only one thread ever reaches the blocks.
unsafe impl Send for FqzGparams {}

unsafe extern "C" {
    fn fqz_compress(
        vers: c_int,
        s: *mut FqzSlice,
        input: *mut c_char,
        in_size: usize,
        out_size: *mut usize,
        strat: c_int,
        gp: *mut FqzGparams,
    ) -> *mut c_char;
}

pub const GFLAG_MULTI_PARAM: c_uint = 1;
pub const GFLAG_HAVE_STAB: c_uint = 2;
pub const GFLAG_DO_REV: c_uint = 4;
pub const PFLAG_DO_DEDUP: c_uint = 2;
pub const PFLAG_DO_LEN: c_uint = 4;
pub const PFLAG_DO_SEL: c_uint = 8;
pub const PFLAG_HAVE_QMAP: c_uint = 16;
pub const PFLAG_HAVE_PTAB: c_uint = 32;
pub const PFLAG_HAVE_DTAB: c_uint = 64;
pub const PFLAG_HAVE_QTAB: c_uint = 128;
/// htscodecs' unset `qmap` entry.
pub const UNMAPPED: c_uint = i32::MAX as c_uint;

/// Qualities of one block, one entry per record in `lens` and `flags`.
#[derive(Debug, Clone)]
pub struct Records {
    pub quals: Vec<u8>,
    pub lens: Vec<u32>,
    pub flags: Vec<u32>,
    /// The distinct quality values, ascending.
    pub alphabet: Vec<u8>,
}

hegel::pretty_print_as_debug!(Records);

/// Compress with htscodecs. `vers` 3 is CRAM 3.1 (reversed records stored
/// in read orientation), 4 is CRAM 4. `gp` replaces the strategy's own
/// parameter choice; htscodecs shifts its tables in place, so it is used up.
pub fn htscodecs_compress(
    records: &Records,
    vers: c_int,
    strat: c_int,
    gp: Option<&mut FqzGparams>,
) -> Option<Vec<u8>> {
    // htscodecs keeps ~350 KB of statistics tables on the stack, more than
    // hegel's test threads have, so it gets a thread of its own.
    std::thread::scope(|scope| {
        std::thread::Builder::new()
            .stack_size(16 << 20)
            .spawn_scoped(scope, || compress_on_this_thread(records, vers, strat, gp))
            .unwrap()
            .join()
            .unwrap()
    })
}

fn compress_on_this_thread(
    records: &Records,
    vers: c_int,
    strat: c_int,
    gp: Option<&mut FqzGparams>,
) -> Option<Vec<u8>> {
    let mut input = records.quals.clone();
    let mut lens = records.lens.clone();
    let mut flags = records.flags.clone();
    let mut slice = FqzSlice {
        num_records: c_int::try_from(lens.len()).unwrap(),
        len: lens.as_mut_ptr(),
        flags: flags.as_mut_ptr(),
    };
    let mut out_size = 0usize;
    let gp = gp.map_or(std::ptr::null_mut(), std::ptr::from_mut);
    // SAFETY: `slice` points at `num_records` lengths and flags that live
    // until the call returns, the lengths sum to `input.len()`, and `gp` is
    // null or a parameter set whose `p` holds `nparam` blocks (see
    // `CustomParams::gparams`). htscodecs only mutates what it is handed.
    let out = unsafe {
        fqz_compress(
            vers,
            &raw mut slice,
            input.as_mut_ptr().cast::<c_char>(),
            input.len(),
            &raw mut out_size,
            strat,
            gp,
        )
    };
    if out.is_null() {
        return None;
    }
    // SAFETY: on success htscodecs returns a malloc'd buffer of `out_size` bytes.
    let bytes = unsafe { std::slice::from_raw_parts(out.cast::<u8>(), out_size) }.to_vec();
    // SAFETY: `out` came from htscodecs' malloc and is freed once.
    unsafe { libc::free(out.cast()) };
    Some(bytes)
}

// ── Generators ───────────────────────────────────────────────────────

/// Records with a quality alphabet of at most `max_alphabet` values,
/// correlated runs, duplicated neighbours, and random reverse / READ2 flags.
#[hegel::composite]
pub fn arb_records(tc: &TestCase, max_alphabet: usize) -> Records {
    let mut alphabet = tc.draw(gs::vecs(gs::integers::<u8>()).min_size(1).max_size(max_alphabet));
    alphabet.sort_unstable();
    alphabet.dedup();

    let n = tc.draw(gs::integers::<usize>().min_value(1).max_value(40));
    let draw_len = |tc: &TestCase| {
        // Some records past 1023 qualities, where the position context saturates.
        if tc.draw(gs::weighted_booleans(0.1)) {
            tc.draw(gs::integers::<u32>().min_value(1000).max_value(2100))
        } else {
            tc.draw(gs::integers::<u32>().min_value(1).max_value(300))
        }
    };
    let fixed = tc.draw(gs::booleans()).then(|| draw_len(tc));
    let dup_rate = tc.draw(gs::sampled_from(&[0.0, 0.2, 0.6]));
    // Each quality repeats the previous one when its byte is below `stick`.
    let stick = tc.draw(gs::integers::<u8>());

    let mut quals = Vec::new();
    let mut lens = Vec::with_capacity(n);
    let mut flags = Vec::with_capacity(n);
    for i in 0..n {
        let dup = i > 0 && tc.draw(gs::weighted_booleans(dup_rate));
        let len = match (dup, fixed) {
            (true, _) => lens[i - 1],
            (false, Some(len)) => len,
            (false, None) => draw_len(tc),
        };
        let rev = if tc.draw(gs::booleans()) { FREVERSE } else { 0 };
        let read2 = if tc.draw(gs::booleans()) { FREAD2 } else { 0 };
        flags.push(rev | read2);
        if dup {
            let start = quals.len() - len as usize;
            quals.extend_from_within(start..);
        } else {
            let bytes = tc.draw(gs::binary().min_size(len as usize).max_size(len as usize));
            let mut prev = alphabet[0];
            for (j, b) in bytes.into_iter().enumerate() {
                if j == 0 || b >= stick {
                    prev = alphabet[usize::from(b) % alphabet.len()];
                }
                quals.push(prev);
            }
        }
        lens.push(len);
    }
    Records { quals, lens, flags, alphabet }
}

/// A monotonically non-decreasing table of `N` entries with values in
/// `0..=max`, as htscodecs' `store_array` requires: a few steps up, of
/// random height (skipped values are empty runs), at random places.
pub fn arb_table<const N: usize>(tc: &TestCase, max: u32) -> [c_uint; N] {
    let mut steps: Vec<(usize, u32)> = tc.draw(
        gs::vecs(gs::tuples!(
            gs::integers::<usize>().max_value(N - 1),
            gs::integers::<u32>().min_value(1).max_value(40)
        ))
        .max_size(12),
    );
    steps.sort_unstable();
    let mut table = [0; N];
    let mut value = 0u32;
    let mut steps = steps.into_iter().peekable();
    for (i, entry) in table.iter_mut().enumerate() {
        while let Some(&(at, by)) = steps.peek()
            && at <= i
        {
            value = (value + by).min(max);
            steps.next();
        }
        *entry = value;
    }
    // htscodecs cannot read back a table whose last run is a multiple of
    // 255: `store_array` ends it with a 0 continuation byte that `read_array`
    // stops short of, having already counted the whole table. Its built-in
    // strategies never make one; lengthen such a run by one.
    let last_run = table.iter().rev().take_while(|&&v| v == value).count();
    if last_run % 255 == 0 {
        table[N - last_run - 1] = value;
    }
    table
}

fn nibble(tc: &TestCase) -> c_uint {
    tc.draw(gs::integers::<c_uint>().max_value(15))
}

/// A parameter set htscodecs can encode `records` with, plus the per-record
/// selector values it expects in the top 16 bits of each flag.
pub struct CustomParams {
    pub gp: FqzGparams,
    pub blocks: Vec<FqzParam>,
    pub selectors: Vec<u32>,
}

impl CustomParams {
    /// The global parameters, pointing at `blocks`; valid while `self` is.
    pub fn gparams(&mut self) -> &mut FqzGparams {
        self.gp.p = self.blocks.as_mut_ptr();
        &mut self.gp
    }
}

impl std::fmt::Debug for CustomParams {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("CustomParams")
            .field("gflags", &self.gp.gflags)
            .field("nparam", &self.gp.nparam)
            .field("max_sel", &self.gp.max_sel)
            .field("max_sym", &self.gp.max_sym)
            .field("stab", &&self.gp.stab[..])
            .field(
                "blocks",
                &self
                    .blocks
                    .iter()
                    .map(|b| {
                        (b.context, b.pflags, b.max_sym, [b.qbits, b.qloc, b.ploc, b.dloc, b.sloc])
                    })
                    .collect::<Vec<_>>(),
            )
            .field("selectors", &self.selectors)
            .finish()
    }
}

/// Generate a parameter set for `records`. Constraints come from what
/// htscodecs' encoder assumes of a hand-made set: the first block drives the
/// quality mapping and context for every record, a selector needs
/// `max_sel > 0`, and a fixed-length block needs equal lengths.
pub fn arb_params(tc: &TestCase, records: &Records) -> CustomParams {
    let nparam = tc.draw(gs::integers::<c_int>().min_value(1).max_value(3));
    let multi = nparam > 1 || tc.draw(gs::booleans());
    // htscodecs encodes a selector whenever there are several blocks, but
    // decodes one only if the first block has do_sel.
    let use_sel = multi || tc.draw(gs::booleans());
    let have_stab = (use_sel && nparam == 1) || tc.draw(gs::booleans());
    let do_rev = tc.draw(gs::booleans());

    let (max_sel, stab) = if have_stab {
        let max_sel = tc.draw(gs::integers::<c_int>().min_value(c_int::from(use_sel)).max_value(5));
        (max_sel, arb_table::<256>(tc, (nparam - 1) as u32))
    } else {
        // As the decoder derives it.
        (if nparam > 1 { nparam } else { 0 }, [0; 256])
    };

    let selectors = records
        .lens
        .iter()
        .map(|_| {
            if !use_sel {
                0
            } else if have_stab {
                tc.draw(gs::integers::<u32>().max_value(max_sel as u32))
            } else {
                tc.draw(gs::integers::<u32>().max_value((nparam - 1) as u32))
            }
        })
        .collect();

    let equal_lengths = records.lens.windows(2).all(|w| w[0] == w[1]);
    let top = *records.alphabet.last().unwrap();
    let blocks: Vec<FqzParam> = (0..nparam)
        .map(|i| {
            let store_qmap = records.alphabet.len() < 256 && tc.draw(gs::booleans());
            let mut qmap = [UNMAPPED; 256];
            let max_sym = if store_qmap {
                for (rank, &q) in records.alphabet.iter().enumerate() {
                    qmap[usize::from(q)] = rank as c_uint;
                }
                records.alphabet.len() as c_int
            } else {
                qmap = std::array::from_fn(|q| q as c_uint);
                // Room for every quality; more only makes the models larger.
                tc.draw(gs::integers::<c_int>().min_value(c_int::from(top)).max_value(255))
            };

            let qbits = nibble(tc);
            let use_qtab = tc.draw(gs::booleans());
            let qtab = if qbits > 0 && use_qtab {
                arb_table::<256>(tc, 500)
            } else {
                std::array::from_fn(|q| q as c_uint)
            };
            let use_ptab = tc.draw(gs::booleans());
            let ptab = if use_ptab { arb_table::<1024>(tc, 500) } else { [0; 1024] };
            let use_dtab = tc.draw(gs::booleans());
            let dtab = if use_dtab { arb_table::<256>(tc, 500) } else { [0; 256] };

            let do_sel = if i == 0 { use_sel } else { max_sel > 0 && tc.draw(gs::booleans()) };
            let fixed_len = equal_lengths && tc.draw(gs::booleans());
            let do_dedup = tc.draw(gs::booleans());
            let reserved = tc.draw(gs::booleans());
            let flag = |on: bool, f: c_uint| if on { f } else { 0 };
            let pflags = flag(use_qtab, PFLAG_HAVE_QTAB)
                | flag(use_dtab, PFLAG_HAVE_DTAB)
                | flag(use_ptab, PFLAG_HAVE_PTAB)
                | flag(do_sel, PFLAG_DO_SEL)
                | flag(fixed_len, PFLAG_DO_LEN)
                | flag(do_dedup, PFLAG_DO_DEDUP)
                | flag(store_qmap, PFLAG_HAVE_QMAP)
                | flag(reserved, 1);

            FqzParam {
                context: tc.draw(gs::integers::<u16>()),
                pflags,
                do_sel: c_uint::from(do_sel),
                do_dedup: c_uint::from(do_dedup),
                store_qmap: c_uint::from(store_qmap),
                fixed_len: c_uint::from(fixed_len),
                use_qtab: u8::from(use_qtab),
                use_dtab: u8::from(use_dtab),
                use_ptab: u8::from(use_ptab),
                qbits,
                qloc: nibble(tc),
                // Only "is there a table" matters to the encoder.
                pbits: c_uint::from(use_ptab),
                ploc: nibble(tc),
                dbits: c_uint::from(use_dtab),
                dloc: nibble(tc),
                sbits: 0,
                sloc: nibble(tc),
                max_sym,
                nsym: 0,
                max_sel: 0,
                qmap,
                qtab,
                ptab,
                dtab,
                qshift: nibble(tc) as c_int,
                pshift: 0,
                dshift: 0,
                sshift: 0,
                qmask: (1 << qbits) - 1,
                do_r2: 0,
                do_qa: 0,
            }
        })
        .collect();

    let gflags = if multi { GFLAG_MULTI_PARAM } else { 0 }
        | if have_stab { GFLAG_HAVE_STAB } else { 0 }
        | if do_rev { GFLAG_DO_REV } else { 0 };
    let gp = FqzGparams {
        vers: 5,
        gflags,
        nparam,
        max_sel,
        stab,
        max_sym: blocks.iter().map(|b| b.max_sym).max().unwrap(),
        p: std::ptr::null_mut(),
    };
    CustomParams { gp, blocks, selectors }
}

// ── Streams ──────────────────────────────────────────────────────────

/// Records, and htscodecs' stream of them with one of its built-in
/// strategies: 0..=3 are tuned ones, 4 its all-zero "custom" one. `vers` 3
/// (CRAM 3.1) stores reverse-strand records reversed, 4 does not. `None` if
/// htscodecs declines (its output buffer estimate is too small).
pub fn builtin_stream(tc: &TestCase) -> Option<(Records, Vec<u8>)> {
    let max_alphabet = tc.draw(gs::sampled_from(&[4, 8, 64, 256]));
    let records = tc.draw(arb_records(max_alphabet));
    let strat = tc.draw(gs::integers::<c_int>().min_value(0).max_value(4));
    let vers = tc.draw(gs::sampled_from(&[3, 4]));
    let stream = htscodecs_compress(&records, vers, strat, None)?;
    Some((records, stream))
}

/// Records, and htscodecs' stream of them with a generated parameter set.
pub fn custom_stream(tc: &TestCase) -> Option<(Records, Vec<u8>)> {
    let max_alphabet = tc.draw(gs::sampled_from(&[4, 40, 256]));
    let mut records = tc.draw(arb_records(max_alphabet));
    let mut params = arb_params(tc, &records);
    tc.note(&format!("{params:?}"));
    for (flag, sel) in records.flags.iter_mut().zip(&params.selectors) {
        *flag |= sel << 16;
    }
    let stream = htscodecs_compress(&records, 3, 0, Some(params.gparams()))?;
    Some((records, stream))
}
