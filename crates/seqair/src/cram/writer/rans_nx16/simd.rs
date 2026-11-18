//! The rANS Nx16 encoder's fast path (`r[cram_writer.rans_nx16.simd]`):
//! order-1 counts that do not serialize on repeated pairs, order-1 tables
//! for the alphabet's rows only, and SIMD kernels for the interleaved
//! states, as htscodecs' `rANS_static32x16pr_neon.c` codes them. Every
//! function here has a scalar counterpart in the parent module, the
//! reference it must match byte for byte.

use super::{
    Coder, LOWER_BOUND, ORDER_0_BITS, Symbol, context, normalize_into, order_1_header,
    order_1_shift, target, write_row,
};
use fearless_simd::{Level, dispatch, prelude::*, u8x16, u32x4, u64x2};
use fearless_simd_macros::simd;
use std::cell::RefCell;
use std::thread::LocalKey;

// r[impl cram_writer.rans_nx16.simd]

/// A symbol as the kernels code it: the reciprocal `m` in the low half, and
/// in the high half the start (bits 0–11), the frequency (12–24) and the
/// shift `s` (25–28). 0 for a symbol without frequency, never looked up.
type Enc = u64;

/// htscodecs' `RansEncSymbolInit` reciprocal, with the state doubled
/// instead of the shift lowered by one so that `f = 1` needs no bias:
/// `q = ⌊2x·m / 2^32⌋ >> s = ⌊x·m / 2^(31+s)⌋` with `s = ⌈log2 f⌉` and
/// `m = ⌈2^(31+s) / f⌉`, exact for every `x < 2^31` (Alverson; the error
/// `m·f − 2^(31+s) < f ≤ 2^s` times `x` stays below `2^(31+s)`).
#[allow(clippy::arithmetic_side_effects, reason = "1 ≤ f ≤ 4096, so s ≤ 12 and m < 2^32")]
#[allow(clippy::cast_possible_truncation, reason = "m < 2^32")]
pub(super) fn enc(sym: Symbol) -> Enc {
    let f = sym.freq;
    if f == 0 {
        return 0;
    }
    let s = 32 - (f - 1).leading_zeros();
    let m = ((1u64 << (31 + s)) + u64::from(f) - 1) / u64::from(f);
    debug_assert!(m <= u64::from(u32::MAX) && sym.start < 1 << 12 && f <= 1 << 12);
    u64::from(m as u32) | (u64::from(sym.start | f << 12 | s << 25) << 32)
}

/// The renormalization words, written backwards from the end of `buf` so
/// that `buf[at..]` holds them in decode order. Every symbol emits at most
/// one word, so a buffer of 16 bytes plus 2 a symbol never runs out, and
/// keeps `at ≥ 16` for the kernels' 16-byte stores below it.
struct Words<'a> {
    buf: &'a mut [u8],
    at: usize,
}

impl<'a> Words<'a> {
    /// Room for `symbols` symbols' words in `scratch`, whose old contents
    /// are overwritten before they are read; `None` if that overflows.
    fn new(scratch: &'a mut Vec<u8>, symbols: usize) -> Option<Self> {
        scratch.resize(symbols.checked_mul(2)?.checked_add(16)?, 0);
        Some(Self { at: scratch.len(), buf: scratch })
    }

    #[allow(clippy::indexing_slicing, reason = "at ≤ buf.len()")]
    fn decode_order(&self) -> &[u8] {
        &self.buf[self.at..]
    }
}

thread_local! {
    /// [`Words`]' buffer, kept between blocks: a fresh one per block is a
    /// page fault per 4 KiB of words.
    static WORDS: RefCell<Vec<u8>> = const { RefCell::new(Vec::new()) };
    /// The order-1 tables, kept between blocks for the same reason: zeroing
    /// their 1.5 MiB cost a 4 KiB block several times its coding.
    static ORDER_1: RefCell<Order1Scratch> = RefCell::new(Order1Scratch::default());
}

/// A word buffer larger than this (a block of 16 MiB) is dropped after use
/// rather than kept for the thread's next block.
const MAX_KEPT_WORDS: usize = 32 << 20;

/// Per-thread scratch: what to give back after use.
trait Scratch: Default {
    fn trim(&mut self) {}
}

impl Scratch for Vec<u8> {
    fn trim(&mut self) {
        if self.capacity() > MAX_KEPT_WORDS {
            *self = Vec::new();
        }
    }
}

/// `f` on this thread's scratch, or on a fresh one while the thread-local
/// is being destroyed or already in use (neither happens in the encoder's
/// own calls, but neither may panic).
fn with_scratch<T: Scratch, R>(
    key: &'static LocalKey<RefCell<T>>,
    f: impl FnOnce(&mut T) -> Option<R>,
) -> Option<R> {
    let mut f = Some(f);
    let kept = key.try_with(|cell| {
        let mut scratch = cell.try_borrow_mut().ok()?;
        let r = f.take().map(|f| f(&mut scratch));
        scratch.trim();
        r
    });
    match kept {
        Ok(Some(r)) => r,
        // `f` has not run.
        _ => f.take()?(&mut T::default()),
    }
}

/// Order-1 counts, in four tables so that a run of one pair is four chains
/// of increments rather than one, and the coding table, each indexed by the
/// context and the symbol as a big-endian u16.
struct Order1Scratch {
    counts: [Vec<u32>; 4],
    /// Whether `counts` may be non-zero: only after a block given up part
    /// way, since a block's tables take every count they read.
    dirty: bool,
    /// The alphabet's rows over the alphabet's columns, row after row.
    rows: Vec<u32>,
    freqs: Vec<u32>,
    /// Entries outside the current alphabet's rows and columns hold an
    /// earlier block's, which no position of this block reads.
    enc: Vec<Enc>,
}

impl Scratch for Order1Scratch {}

impl Default for Order1Scratch {
    fn default() -> Self {
        Self {
            counts: std::array::from_fn(|_| vec![0; 1 << 16]),
            dirty: false,
            rows: Vec::new(),
            freqs: Vec::new(),
            enc: vec![0; 1 << 16],
        }
    }
}

/// One state step on a scalar, with [`enc`]'s reciprocal.
#[inline(always)]
#[allow(clippy::indexing_slicing, reason = "at ≥ 2 while symbols remain (see `Words`)")]
#[allow(clippy::arithmetic_side_effects, reason = "as in `put`; q·(M − f) + x % f + start < 2^31")]
#[allow(clippy::cast_possible_truncation, reason = "the low word; the high half of a u64 product")]
fn put(x: u32, e: Enc, bits: u32, words: &mut Words<'_>) -> u32 {
    let info = (e >> 32) as u32;
    let f = (info >> 12) & 0x1FFF;
    // Branch-free: the word is stored either way, below the words so far,
    // and kept only if the state renormalizes.
    let renorm = x >> (31 - bits) >= f;
    let at = words.at;
    debug_assert!(at >= 2, "a word buffer holds 2 bytes a symbol");
    words.buf[at - 2..at].copy_from_slice(&(x as u16).to_le_bytes());
    words.at = at - 2 * usize::from(renorm);
    let x = if renorm { x >> 16 } else { x };
    let q = ((u64::from(x + x) * (e & 0xFFFF_FFFF)) >> 32) as u32 >> (info >> 25);
    x + (info & 0xFFF) + q * ((1 << bits) - f)
}

/// For each 4-lane renormalization mask, the byte shuffle that moves the
/// low u16 of each lane that renormalizes to the top of 16 bytes, in lane
/// order: stored ending at the write position, they are in decode order.
#[allow(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    reason = "const table build: m < 16, lane < 4, k ≤ 4"
)]
const EMIT4: [[u8; 16]; 16] = {
    let mut t = [[0x80u8; 16]; 16];
    let mut m = 0;
    while m < 16 {
        let k = (m as u32).count_ones() as usize;
        let mut slot = 16 - 2 * k;
        let mut lane = 0;
        while lane < 4 {
            if m & (1 << lane) != 0 {
                t[m][slot] = (4 * lane) as u8;
                t[m][slot + 1] = (4 * lane + 1) as u8;
                slot += 2;
            }
            lane += 1;
        }
        m += 1;
    }
    t
};

/// Bytes emitted for a 4-lane mask: 2 · popcount.
const EMITTED: [u8; 16] = [0, 2, 2, 4, 2, 4, 4, 6, 2, 4, 4, 6, 4, 6, 6, 8];

/// Four states' steps at once: [`put`] lane by lane, the words of the lanes
/// that renormalize stored in lane order below the ones before.
/// `latency_bound`: the caller has no other states to step meanwhile.
#[inline(always)]
#[allow(clippy::indexing_slicing, reason = "the mask of 4 lanes is < 16; at ≥ 16 (see `Words`)")]
#[allow(clippy::cast_possible_truncation, reason = "the mask of 4 lanes is < 16")]
#[allow(clippy::arithmetic_side_effects, reason = "vector lanes wrap; as in `put`")]
fn step4<S: Simd>(
    simd: S,
    x: u32x4<S>,
    e: &[Enc; 4],
    bits: u32,
    words: &mut Words<'_>,
    latency_bound: bool,
) -> u32x4<S> {
    let (m, info) = split(simd, e);
    let f = (info >> 12) & 0x1FFF;
    let renorm = (x >> (31 - bits)).simd_ge(f);
    let mask = (renorm.to_bitmask() & 0xF) as usize;
    let at = words.at;
    debug_assert!(at >= 16, "a word buffer starts with 16 spare bytes");
    let low: u8x16<S> = x.bitcast();
    low.swizzle_dyn(u8x16::from_slice(simd, &EMIT4[mask])).store_slice(&mut words.buf[at - 16..at]);
    words.at = at - usize::from(EMITTED[mask]);
    // With one vector of states the loop waits on each step: the quotients
    // of both candidate states are taken alongside the comparison, which
    // then only picks one (4-state order 1 +5% on an M4). With eight vectors
    // the loop is throughput-bound and the extra multiply costs more than it
    // hides (32-state order 0 −7% on an M4), as it does on Zen 2 even with
    // one vector (4-state order 0 −15%: two `pmuludq` and shuffles a step).
    let (x, q) = if latency_bound && cfg!(target_arch = "aarch64") {
        let x16 = x >> 16;
        (renorm.select(x16, x), renorm.select(quotient(x16, m, info), quotient(x, m, info)))
    } else {
        let x = renorm.select(x >> 16, x);
        (x, quotient(x, m, info))
    };
    x + (info & 0xFFF) + q * (u32x4::splat(simd, 1 << bits) - f)
}

/// Four symbols' reciprocals and their start, frequency and shift.
#[inline(always)]
fn split<S: Simd>(simd: S, e: &[Enc; 4]) -> (u32x4<S>, u32x4<S>) {
    let (lo, hi) = e.split_at(2);
    let a: u32x4<S> = u64x2::from_slice(simd, lo).bitcast();
    let b: u32x4<S> = u64x2::from_slice(simd, hi).bitcast();
    (a.unzip_low(b), a.unzip_high(b))
}

/// The kernels' quotients `⌊x / f⌋` of four states by four symbols'
/// frequencies, for tests.
#[cfg(test)]
#[inline(always)]
pub(super) fn quotients<S: Simd>(simd: S, x: [u32; 4], e: &[Enc; 4]) -> [u32; 4] {
    let (m, info) = split(simd, e);
    let mut q = [0; 4];
    quotient(u32x4::from_slice(simd, &x), m, info).store_slice(&mut q);
    q
}

/// `⌊x / f⌋` lane by lane, with [`enc`]'s reciprocal; `x < 2^31`.
#[inline(always)]
#[allow(clippy::arithmetic_side_effects, reason = "vector lanes wrap")]
fn quotient<S: Simd>(x: u32x4<S>, m: u32x4<S>, info: u32x4<S>) -> u32x4<S> {
    let (lo, hi) = (x + x).widening_mul(m);
    lo.bitcast::<u32x4<S>>().unzip_high(hi.bitcast::<u32x4<S>>()) >> (info >> 25)
}

/// The parent's `encode_order_0` coding loop: `n` states, the table already
/// in `out`; the states and words follow.
#[allow(clippy::arithmetic_side_effects, reason = "i % n ≤ i")]
#[allow(clippy::indexing_slicing, reason = "a u8 indexes a 256-entry table; l < n ≤ 32")]
pub(super) fn code_order_0(
    level: Level,
    data: &[u8],
    table: &[Symbol; 256],
    n: usize,
    out: &mut Vec<u8>,
) -> Option<()> {
    debug_assert!(n == 4 || n == 32);
    let enc_table = table.map(enc);
    let whole = data.len() - data.len() % n;
    let body = &data[..whole];
    let mut states = [LOWER_BOUND; 32];
    with_scratch(&WORDS, |scratch| {
        let mut words = Words::new(scratch, data.len())?;
        // The last, partial round of `n`: byte `i` with state `i % n`.
        debug_assert!(whole <= data.len() && data.len() - whole < n);
        for (l, &b) in data[whole..].iter().enumerate().rev() {
            states[l] = put(states[l], enc_table[usize::from(b)], ORDER_0_BITS, &mut words);
        }
        let groups = states.as_chunks_mut::<4>().0;
        if let Some(s) = groups.first_chunk_mut::<8>().filter(|_| n == 32) {
            dispatch!(level, simd => code_o0(simd, body.as_chunks::<32>().0, &enc_table, s, &mut words));
        } else if let Some(s) = groups.first_chunk_mut::<1>() {
            dispatch!(level, simd => code_o0(simd, body.as_chunks::<4>().0, &enc_table, s, &mut words));
        }
        finish(out, &states[..n], words.decode_order());
        Some(())
    })
}

/// Order 0 over whole rounds of `N` bytes, last round first, each state's
/// step in groups of four lanes, the highest first (the order the scalar
/// encoder pushes words in).
#[simd]
#[allow(clippy::indexing_slicing, reason = "u8s index a 256-entry table")]
fn code_o0<S: Simd, const N: usize, const M: usize>(
    simd: S,
    rounds: &[[u8; N]],
    enc_table: &[Enc; 256],
    states: &mut [[u32; 4]; M],
    words: &mut Words<'_>,
) {
    debug_assert_eq!(Some(N), M.checked_mul(4));
    // Local copies, so that the write position and the states stay in
    // registers: through memory, a store-to-load forward joins every
    // state's chain (a quarter of the 4-state loop's time on Zen 2).
    let mut local = Words { buf: &mut *words.buf, at: words.at };
    let mut xs = load_states(simd, states);
    for round in rounds.iter().rev() {
        for (x, bytes) in xs.iter_mut().zip(round.as_chunks::<4>().0).rev() {
            let e = [
                enc_table[usize::from(bytes[0])],
                enc_table[usize::from(bytes[1])],
                enc_table[usize::from(bytes[2])],
                enc_table[usize::from(bytes[3])],
            ];
            *x = step4(simd, *x, &e, ORDER_0_BITS, &mut local, M == 1);
        }
    }
    store_states(&xs, states);
    words.at = local.at;
}

#[inline(always)]
fn load_states<S: Simd, const M: usize>(simd: S, states: &[[u32; 4]; M]) -> [u32x4<S>; M] {
    let mut xs = [u32x4::splat(simd, 0); M];
    for (x, s) in xs.iter_mut().zip(states) {
        *x = u32x4::from_slice(simd, s);
    }
    xs
}

#[inline(always)]
fn store_states<S: Simd, const M: usize>(xs: &[u32x4<S>; M], states: &mut [[u32; 4]; M]) {
    for (x, s) in xs.iter().zip(states) {
        x.store_slice(s);
    }
}

/// The parent's `encode_order_1` with the fast path's counts and tables
/// and the SIMD coding loop. `data` holds at least 8 bytes.
#[allow(clippy::indexing_slicing, reason = "u8s and u16s index 256- and 65536-entry tables")]
#[allow(clippy::arithmetic_side_effects, reason = "j·chunk ≤ len; n ≥ 1")]
pub(super) fn encode_order_1(data: &[u8], n: usize, level: Level) -> Option<Vec<u8>> {
    let chunk = data.len() / n;
    with_scratch(&ORDER_1, |o1| {
        let (shift, table) = order_1_tables(o1, data, n)?;
        let enc_table: &[Enc; 1 << 16] = o1.enc.as_slice().try_into().ok()?;
        let mut out = order_1_header(shift, &table, Coder::Simd(level))?;
        let mut states = [LOWER_BOUND; 32];
        with_scratch(&WORDS, |scratch| {
            let mut words = Words::new(scratch, data.len())?;
            // The bytes after the last whole stretch, with the last state.
            for p in (chunk * n..data.len()).rev() {
                let ctx = context(data, chunk, n, p);
                let e = enc_table[usize::from(ctx) << 8 | usize::from(data[p])];
                states[n - 1] = put(states[n - 1], e, shift, &mut words);
            }
            if chunk > 0 {
                let body = &data[..chunk * n];
                let groups = states.as_chunks_mut::<4>().0;
                if let Some(s) = groups.first_chunk_mut::<8>().filter(|_| n == 32) {
                    dispatch!(level, simd => code_o1(simd, body, chunk, enc_table, shift, s, &mut words));
                } else if let Some(s) = groups.first_chunk_mut::<1>() {
                    dispatch!(level, simd => code_o1(simd, body, chunk, enc_table, shift, s, &mut words));
                }
            }
            finish(&mut out, &states[..n], words.decode_order());
            Some(())
        })?;
        Some(out)
    })
}

/// Rows of at most this many symbols are gathered from the four count
/// tables symbol by symbol; longer ones are summed whole, a vector at a time.
const GATHERED: usize = 32;

/// The parent's `order_1_tables` and the table writing of its
/// `encode_order_1`, over the alphabet's rows and columns only: the shift
/// and the tables as written, with `o1.enc` filled for every pair of the
/// alphabet.
#[allow(clippy::indexing_slicing, reason = "u8s and u16s index 256- and 65536-entry tables")]
#[allow(
    clippy::arithmetic_side_effects,
    reason = "counts ≤ data.len(); each stretch start was counted once as a pair; \
              shift ≥ log2 of a row total; a row's starts sum to its total ≤ 4096"
)]
#[allow(clippy::cast_possible_truncation, reason = "an index into a 256-symbol alphabet")]
fn order_1_tables(o1: &mut Order1Scratch, data: &[u8], n: usize) -> Option<(u32, Vec<u8>)> {
    if o1.dirty {
        for c in &mut o1.counts {
            c.fill(0);
        }
    }
    o1.dirty = true;
    let [c0, c1, c2, c3] = &mut o1.counts;
    let mut alphabet = [false; 256];
    alphabet[0] = true;
    for &b in data {
        alphabet[usize::from(b)] = true;
    }
    // The alphabet in order, and each symbol's place in it (0 is first).
    let syms: Vec<usize> = (0..256).filter(|&s| alphabet[s]).collect();
    let mut index = [0u8; 256];
    for (i, &s) in syms.iter().enumerate() {
        index[s] = i as u8;
    }
    let k = syms.len();

    // Every byte with the one before it (the first with context 0), the
    // pairs taking the four tables in turn.
    let (&first, rest) = data.split_first()?;
    c0[usize::from(first)] += 1;
    let mut prev = usize::from(first);
    let (quads, tail) = rest.as_chunks::<4>();
    for &[a, b, c, d] in quads {
        let (a, b, c, d) = (usize::from(a), usize::from(b), usize::from(c), usize::from(d));
        c0[prev << 8 | a] += 1;
        c1[a << 8 | b] += 1;
        c2[b << 8 | c] += 1;
        c3[c << 8 | d] += 1;
        prev = d;
    }
    for &a in tail {
        c0[prev << 8 | usize::from(a)] += 1;
        prev = usize::from(a);
    }
    // Each row's counts over the alphabet, taken out of the tables, which
    // are left zero: no pair outside the alphabet was counted.
    o1.rows.clear();
    o1.rows.resize(k * k, 0);
    for (row, &ctx) in o1.rows.chunks_exact_mut(k).zip(&syms) {
        let at = ctx << 8;
        let [r0, r1, r2, r3] = [
            &mut c0[at..][..256],
            &mut c1[at..][..256],
            &mut c2[at..][..256],
            &mut c3[at..][..256],
        ];
        if k <= GATHERED {
            for (r, &s) in row.iter_mut().zip(&syms) {
                *r = std::mem::take(&mut r0[s])
                    + std::mem::take(&mut r1[s])
                    + std::mem::take(&mut r2[s])
                    + std::mem::take(&mut r3[s]);
            }
        } else {
            let mut sum = [0u32; 256];
            for ((((s, a), b), c), d) in sum.iter_mut().zip(r0).zip(r1).zip(r2).zip(r3) {
                *s = std::mem::take(a) + std::mem::take(b) + std::mem::take(c) + std::mem::take(d);
            }
            for (r, &s) in row.iter_mut().zip(&syms) {
                *r = sum[s];
            }
        }
    }
    o1.dirty = false;
    let rows = &mut o1.rows;
    // Each stretch start but the first moved to context 0, which `context`
    // gives it; context 0's row is the first.
    let chunk = data.len() / n;
    if chunk > 0 {
        for p in (1..n).map(|j| j * chunk) {
            let (before, at) = (
                usize::from(index[usize::from(data[p - 1])]),
                usize::from(index[usize::from(data[p])]),
            );
            rows[before * k + at] -= 1;
            rows[at] += 1;
        }
    }

    let mut targets = Vec::with_capacity(k);
    for row in rows.chunks_exact_mut(k) {
        let (total, symbols) = row
            .iter()
            .fold((0u32, 0usize), |(t, s), &c| (t.saturating_add(c), s + usize::from(c > 0)));
        targets.push(if total == 0 {
            // A symbol that never precedes another: a row that only ever
            // decodes symbol 0, as the parent writes it.
            row[0] = 1;
            target(1, 1)
        } else {
            target(total, symbols)
        });
    }
    let shift = order_1_shift(targets.iter().copied());

    let mut table = Vec::new();
    super::super::rans::write_runs(&mut table, &alphabet, |_, _| {});
    o1.freqs.resize(k, 0);
    let freqs = &mut o1.freqs;
    for ((row, &t), &ctx) in rows.chunks_exact(k).zip(&targets).zip(&syms) {
        normalize_into(row, t, freqs);
        write_row(&mut table, freqs.iter().copied())?;
        let total: u32 = freqs.iter().sum();
        let scale = shift - total.trailing_zeros();
        let enc_row = &mut o1.enc[ctx << 8..][..256];
        let mut start = 0;
        for (&f, &s) in freqs.iter().zip(&syms) {
            let freq = f << scale;
            enc_row[s] = enc(Symbol { start, freq });
            start += freq;
        }
    }
    Some((shift, table))
}
/// Positions an order-1 block covers: each lane's bytes for them, and the
/// one before, are copied out once, so the steps index a fixed-size array.
const BLOCK: usize = 16;

/// Order 1 over `4·M` stretches of `chunk` bytes, position by position
/// from the last, the states in groups of four lanes, the highest first. A
/// stretch's first byte has context 0.
#[simd]
#[allow(
    clippy::indexing_slicing,
    reason = "j·chunk + i < 4·M·chunk = data.len(); i ≥ BLOCK in i - BLOCK; u16s index the table"
)]
#[allow(clippy::arithmetic_side_effects, reason = "j·chunk + i < data.len(); i ≥ 1 in i - 1")]
fn code_o1<S: Simd, const M: usize>(
    simd: S,
    data: &[u8],
    chunk: usize,
    enc_table: &[Enc; 1 << 16],
    bits: u32,
    states: &mut [[u32; 4]; M],
    words: &mut Words<'_>,
) {
    let mut local = Words { buf: &mut *words.buf, at: words.at };
    let mut xs = load_states(simd, states);
    // Whole blocks of positions i - BLOCK .. i, all with a byte before them
    // in their stretch: each lane's BLOCK + 1 bytes from i - BLOCK - 1.
    let mut i = chunk;
    while i > BLOCK {
        i -= BLOCK;
        let mut window = [[[0u8; BLOCK + 1]; 4]; M];
        for (j, w) in window.as_flattened_mut().iter_mut().enumerate() {
            let from = j * chunk + i - 1;
            debug_assert!(from + BLOCK < data.len(), "lane {j}'s block ends in its stretch");
            w.copy_from_slice(&data[from..=from + BLOCK]);
        }
        for k in (0..BLOCK).rev() {
            for (x, w) in xs.iter_mut().zip(&window).rev() {
                let e = [
                    enc_table[usize::from(u16::from_be_bytes([w[0][k], w[0][k + 1]]))],
                    enc_table[usize::from(u16::from_be_bytes([w[1][k], w[1][k + 1]]))],
                    enc_table[usize::from(u16::from_be_bytes([w[2][k], w[2][k + 1]]))],
                    enc_table[usize::from(u16::from_be_bytes([w[3][k], w[3][k + 1]]))],
                ];
                *x = step4(simd, *x, &e, bits, &mut local, M == 1);
            }
        }
    }
    // The first positions, the stretch starts among them.
    for i in (0..i).rev() {
        let mut j = 4 * M;
        for x in xs.iter_mut().rev() {
            j -= 4;
            let mut e = [0 as Enc; 4];
            for (l, e) in e.iter_mut().enumerate() {
                let p = (j + l) * chunk + i;
                let ctx = if i == 0 { 0 } else { data[p - 1] };
                *e = enc_table[usize::from(u16::from_be_bytes([ctx, data[p]]))];
            }
            *x = step4(simd, *x, &e, bits, &mut local, M == 1);
        }
    }
    store_states(&xs, states);
    words.at = local.at;
}

/// The states, then the words, already in decode order.
fn finish(out: &mut Vec<u8>, states: &[u32], words: &[u8]) {
    out.reserve(states.len().saturating_mul(4).saturating_add(words.len()));
    for s in states {
        out.extend_from_slice(&s.to_le_bytes());
    }
    out.extend_from_slice(words);
}

#[cfg(test)]
#[allow(clippy::unwrap_used, reason = "test code")]
mod tests {
    use super::*;

    /// A scratch already in use hands out a fresh one rather than
    /// panicking, and a word buffer past the cap is not kept.
    #[test]
    fn scratch_in_use_or_too_large_is_not_kept() {
        let nested = with_scratch(&WORDS, |outer| {
            outer.resize(MAX_KEPT_WORDS + 1, 0);
            with_scratch(&WORDS, |inner: &mut Vec<u8>| Some(inner.capacity()))
        });
        assert_eq!(nested, Some(0), "the nested call gets a fresh buffer");
        assert_eq!(WORDS.with_borrow(Vec::capacity), 0, "an oversized buffer is dropped");
    }
}
