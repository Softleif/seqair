//! rANS 4×8 encoder, order 0 and order 1: the inverse of
//! [`crate::cram::rans`], laid out as htscodecs' `rANS_static.c` writes it.

use crate::cram::varint::write_itf8;
use fearless_simd::Level;
use std::cell::RefCell;
use std::thread::LocalKey;

mod reference;
mod simd;

// r[impl cram_writer.rans4x8]

#[derive(Debug, Clone, Copy)]
pub(super) enum Order {
    Zero,
    One,
}

/// htscodecs normalizes to one less than 4096 (`fsum++`); its decoder
/// accepts 4095 or 4096.
const TOTAL: u32 = 4095;
const LOWER_BOUND: u32 = 1 << 23;

/// A symbol's slot range in a frequency table, and the reciprocal that
/// divides by its frequency.
#[derive(Debug, Clone, Copy, Default)]
struct Symbol {
    start: u32,
    freq: u32,
    /// `⌈2^(31+shift) / freq⌉`.
    rcp: u32,
    shift: u32,
}

impl Symbol {
    /// htscodecs' `RansEncSymbolInit` reciprocal, the state doubled rather
    /// than the shift lowered so `freq = 1` needs no bias:
    /// `⌊x / f⌋ = ⌊2x·m / 2^32⌋ >> s` with `s = ⌈log2 f⌉`,
    /// `m = ⌈2^(31+s) / f⌉`, exact for every `x < 2^31` (the error
    /// `m·f − 2^(31+s) < f ≤ 2^s`, times `x`, stays below `2^(31+s)`).
    #[allow(clippy::arithmetic_side_effects, reason = "1 ≤ freq ≤ 4095: s ≤ 12, m < 2^32")]
    #[allow(clippy::cast_possible_truncation, reason = "m < 2^32")]
    fn new(start: u32, freq: u32) -> Self {
        if freq == 0 {
            return Self { start, freq, rcp: 0, shift: 0 };
        }
        let shift = 32 - (freq - 1).leading_zeros();
        let rcp = ((1u64 << (31 + shift)) + u64::from(freq) - 1) / u64::from(freq);
        debug_assert!(rcp <= u64::from(u32::MAX));
        Self { start, freq, rcp: rcp as u32, shift }
    }

    /// `⌊x / freq⌋` for `x < 2^31`.
    #[inline]
    #[allow(clippy::cast_possible_truncation, reason = "the high half of a u64 product")]
    #[allow(clippy::arithmetic_side_effects, reason = "a u32 times a u32 fits a u64")]
    fn quotient(self, x: u32) -> u32 {
        debug_assert!(x < 1 << 31);
        ((u64::from(x.wrapping_shl(1)) * u64::from(self.rcp)) >> 32) as u32 >> self.shift
    }
}

type Table = [Symbol; 256];
/// A symbol as the fast coding loop reads it: the reciprocal in the low
/// half, and in the high half the start (bits 0–11), the frequency
/// (12–24) and the shift (25–28). 0 for a symbol without frequency, which
/// no position looks up.
type Enc = u64;

impl Symbol {
    #[allow(clippy::arithmetic_side_effects, reason = "start < 4096, freq ≤ 4095, shift ≤ 12")]
    fn pack(self) -> Enc {
        debug_assert!(self.start < 1 << 12 && self.freq < 1 << 13 && self.shift < 1 << 4);
        u64::from(self.rcp) | u64::from(self.start | self.freq << 12 | self.shift << 25) << 32
    }
}

// r[impl cram_writer.rans4x8.fast]
/// How [`encode_with`] counts and codes: the reference loops, or the fast
/// path. Both write the same bytes.
#[derive(Debug, Clone, Copy)]
enum Coder {
    Reference,
    /// The SIMD kernels at a level.
    Fast(Level),
}

/// `data` compressed, or `None` when it is empty, too large for the
/// format's 32-bit sizes, or too small for order 1.
pub(super) fn encode(data: &[u8], order: Order) -> Option<Vec<u8>> {
    encode_with(data, order, Coder::Fast(Level::new()))
}

/// [`encode`] for benchmarks and fuzzing: order 1 or 0, with the reference
/// loops the fast path is checked against, or with the fast path.
#[doc(hidden)]
pub fn encode_for_bench(data: &[u8], order_1: bool, reference: bool) -> Option<Vec<u8>> {
    let order = if order_1 { Order::One } else { Order::Zero };
    encode_with(data, order, if reference { Coder::Reference } else { Coder::Fast(Level::new()) })
}

fn encode_with(data: &[u8], order: Order, coder: Coder) -> Option<Vec<u8>> {
    // An empty table would read as symbol 0 followed by the states.
    if data.is_empty() {
        return None;
    }
    let raw_len = u32::try_from(data.len()).ok()?;
    match coder {
        Coder::Reference => {
            let mut table = Vec::new();
            let mut rev = Vec::with_capacity(data.len());
            let states = match order {
                Order::Zero => reference::encode_order_0(data, &mut table, &mut rev),
                Order::One => reference::encode_order_1(data, &mut table, &mut rev)?,
            };
            rev.reverse();
            assemble(order, raw_len, &table, states, &rev)
        }
        Coder::Fast(level) => match order {
            Order::Zero => fast_order_0(data, raw_len, level),
            Order::One => fast_order_1(data, raw_len, level),
        },
    }
}

/// The stream: order, sizes, table, the four states, and the
/// renormalization bytes in decode order.
fn assemble(
    order: Order,
    raw_len: u32,
    table: &[u8],
    states: [u32; 4],
    words: &[u8],
) -> Option<Vec<u8>> {
    let payload_len = table.len().saturating_add(16).saturating_add(words.len());
    let mut out = Vec::with_capacity(payload_len.saturating_add(9));
    out.push(match order {
        Order::Zero => 0,
        Order::One => 1,
    });
    out.extend_from_slice(&u32::try_from(payload_len).ok()?.to_le_bytes());
    out.extend_from_slice(&raw_len.to_le_bytes());
    out.extend_from_slice(table);
    for state in states {
        out.extend_from_slice(&state.to_le_bytes());
    }
    out.extend_from_slice(words);
    Some(out)
}

/// The renormalization bytes, written backwards from the end of `buf` so
/// that `buf[at..]` holds them in decode order. A symbol emits at most two
/// bytes, so 2 a symbol never run out, and 16 spare keep `at ≥ 16` for the
/// kernels' 16-byte stores below it.
struct Words<'a> {
    buf: &'a mut [u8],
    at: usize,
}

impl<'a> Words<'a> {
    /// Room for `symbols` symbols' bytes in `scratch`, whose old contents
    /// are overwritten before they are read; `None` if that overflows.
    fn new(scratch: &'a mut Vec<u8>, symbols: usize) -> Option<Self> {
        scratch.resize(symbols.checked_mul(2)?.checked_add(16)?, 0);
        Some(Self { at: scratch.len(), buf: scratch })
    }

    fn decode_order(&self) -> &[u8] {
        self.buf.get(self.at..).unwrap_or_default()
    }
}

thread_local! {
    /// [`Words`]' buffer, kept between blocks: a fresh one per block is a
    /// page fault per 4 KiB of bytes.
    static WORDS: RefCell<Vec<u8>> = const { RefCell::new(Vec::new()) };
    /// The order-1 counts and coding table (2.5 MiB), kept between blocks:
    /// zeroing them would cost a small block more than its coding.
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

/// The order-1 pair counts are split over this many tables, taking the
/// pairs in turn, so that a run of one pair (a quality run) is this many
/// chains of increments rather than one.
const PAIR_TABLES: usize = 8;

/// A pair table with 64 bytes of padding, so that one pair's count in two
/// tables is not 4 KiB apart: a store to one would hold up the load from the
/// other.
const PAIR_STRIDE: usize = (1 << 16) + 16;

/// Order-1 counts and the coding table, each indexed by the context and
/// the symbol as a big-endian u16 ([`pair`]).
struct Order1Scratch {
    /// [`PAIR_TABLES`] tables of [`PAIR_STRIDE`] counts, end to end.
    counts: Vec<u32>,
    /// The rows of `counts` the last block may have left non-zero.
    dirty: [bool; 256],
    /// Rows of contexts the current block lacks hold an earlier block's
    /// entries, which no position of this block reads.
    enc: Box<[Enc; 1 << 16]>,
}

impl Scratch for Order1Scratch {}

/// [`Order1Scratch::counts`] as fixed-size tables the counting loops index
/// without bounds checks; `None` never happens.
fn pair_tables(counts: &mut [u32]) -> Option<&mut [[u32; PAIR_STRIDE]; PAIR_TABLES]> {
    counts.as_chunks_mut::<PAIR_STRIDE>().0.try_into().ok()
}

impl Default for Order1Scratch {
    fn default() -> Self {
        // Zeroed by the allocator: no table is built on the stack first.
        let Ok(enc) = vec![0; 1 << 16].into_boxed_slice().try_into() else {
            unreachable!("a Vec of 2^16 elements converts to [Enc; 2^16]")
        };
        Self { counts: vec![0; PAIR_STRIDE * PAIR_TABLES], dirty: [false; 256], enc }
    }
}

/// Order-0 counts, in [`PAIR_TABLES`] tables taking the bytes in turn, so
/// that a run of one byte is that many chains of increments rather than one.
#[allow(clippy::indexing_slicing, reason = "u8s index 256-entry tables")]
fn count_order_0(data: &[u8]) -> [u32; 256] {
    let mut c = [[0u32; 256]; PAIR_TABLES];
    let (rounds, tail) = data.as_chunks::<PAIR_TABLES>();
    for round in rounds {
        for (c, &b) in c.iter_mut().zip(round) {
            c[usize::from(b)] = c[usize::from(b)].wrapping_add(1);
        }
    }
    for &b in tail {
        c[0][usize::from(b)] = c[0][usize::from(b)].wrapping_add(1);
    }
    let (sum, others) = c.split_at_mut(1);
    let mut sum = sum[0];
    for c in others {
        for (s, &n) in sum.iter_mut().zip(&*c) {
            // Each count is at most `data.len()` < 2^32.
            *s = s.wrapping_add(n);
        }
    }
    sum
}

/// A table's symbols packed for the coding loop.
fn pack_row(freq: &[u32; 256], row: &mut [Enc]) {
    let mut start = 0u32;
    for (e, &f) in row.iter_mut().zip(freq) {
        *e = Symbol::new(start, f).pack();
        start = start.saturating_add(f);
    }
}

fn fast_order_0(data: &[u8], raw_len: u32, level: Level) -> Option<Vec<u8>> {
    let freq = normalize(&count_order_0(data));
    let mut table = Vec::new();
    write_frequencies(&mut table, &freq);
    let mut enc = [0; 256];
    pack_row(&freq, &mut enc);
    with_scratch(&WORDS, |scratch| {
        let mut words = Words::new(scratch, data.len())?;
        let states = code_order_0(data, &enc, &mut words, level);
        assemble(Order::Zero, raw_len, &table, states, words.decode_order())
    })
}

fn fast_order_1(data: &[u8], raw_len: u32, level: Level) -> Option<Vec<u8>> {
    // htscodecs reads the byte before the remainder from the last quarter.
    if data.len() < 4 {
        return None;
    }
    with_scratch(&ORDER_1, |o1| {
        let table = order_1_tables(o1, data, level)?;
        let enc = &*o1.enc;
        with_scratch(&WORDS, |scratch| {
            let mut words = Words::new(scratch, data.len())?;
            let states = code_order_1(data, enc, &mut words, level);
            assemble(Order::One, raw_len, &table, states, words.decode_order())
        })
    })
}

/// The order-1 tables as written, with `o1.enc` filled for every context
/// of `data` (at least 4 bytes): every byte counted with the one before it,
/// then the quarter starts moved to context 0.
#[allow(clippy::indexing_slicing, reason = "u8s and u16s index 256- and 65536-entry tables")]
#[allow(
    clippy::arithmetic_side_effects,
    reason = "counts ≤ data.len() < 2^32; each quarter start was counted once as a pair"
)]
fn order_1_tables(o1: &mut Order1Scratch, data: &[u8], level: Level) -> Option<Vec<u8>> {
    let counts = pair_tables(&mut o1.counts)?;
    let dirty = std::mem::replace(&mut o1.dirty, [false; 256]);
    // Only the first table can hold counts: the merge below empties the
    // others' rows of every context it touched.
    for ctx in (0..256).filter(|&c| dirty[c]) {
        counts[0][ctx << 8..][..256].fill(0);
    }
    // Every context is a byte of `data`, or 0 at a quarter's start.
    let mut contexts = [false; 256];
    contexts[0] = true;
    simd::mark_present(level, data, &mut contexts);
    o1.dirty = contexts;

    let (&first, rest) = data.split_first()?;
    counts[0][pair(0, first)] += 1;
    let mut prev = first;
    let (rounds, tail) = rest.as_chunks::<PAIR_TABLES>();
    for round in rounds {
        // The pairs within the round as 16-bit windows of it read
        // big-endian: two operations a pair rather than six.
        let v = u64::from_be_bytes(*round);
        let (first, others) = counts.split_at_mut(1);
        first[0][pair(prev, round[0])] += 1;
        for (k, c) in (1..).zip(others) {
            #[allow(clippy::cast_possible_truncation, reason = "a 16-bit window")]
            let p = (v >> (56 - 8 * k)) as u16;
            c[usize::from(p)] += 1;
        }
        prev = round[PAIR_TABLES - 1];
    }
    for &b in tail {
        counts[0][pair(prev, b)] += 1;
        prev = b;
    }
    let (c0, others) = counts.split_at_mut(1);
    let c0 = &mut c0[0];
    for ctx in (0..256).filter(|&c| contexts[c]) {
        for c in &mut *others {
            for (a, b) in c0[ctx << 8..][..256].iter_mut().zip(&mut c[ctx << 8..][..256]) {
                *a += std::mem::take(b);
            }
        }
    }
    let quarter = data.len() / 4;
    for p in (1..4).map(|j| j * quarter) {
        c0[pair(data[p - 1], data[p])] -= 1;
        c0[pair(0, data[p])] += 1;
    }

    let rows = c0.as_chunks::<256>().0;
    let mut active = [false; 256];
    for ((on, row), &ctx) in active.iter_mut().zip(rows).zip(&contexts) {
        *on = ctx && row.iter().any(|&c| c > 0);
    }
    let mut table = Vec::new();
    let enc = &mut o1.enc;
    write_runs(&mut table, &active, |table, ctx| {
        if let Some(row) = rows.get(ctx) {
            let freq = normalize(row);
            write_frequencies(table, &freq);
            pack_row(&freq, &mut enc[ctx << 8..][..256]);
        }
    });
    Some(table)
}

/// One state's step: renormalize (0–2 bytes, `x` being below 2^31 and
/// `x_max = f·2^19` at least 2^19), then encode, dividing by the
/// reciprocal. Branch-free: both bytes are stored either way, below the
/// bytes so far, and kept only as far as the state renormalizes.
#[inline(always)]
#[allow(clippy::indexing_slicing, reason = "at ≥ 2 while symbols remain (see `Words`)")]
#[allow(
    clippy::arithmetic_side_effects,
    reason = "n ≤ 2 ≤ at; q·(4096 − f) + x % f + start < 2^31"
)]
#[allow(clippy::cast_possible_truncation, reason = "the low bytes; the high half of a u64 product")]
fn put(x: u32, e: Enc, words: &mut Words<'_>) -> u32 {
    let info = (e >> 32) as u32;
    let f = (info >> 12) & 0x1FFF;
    let x_max = f << 19;
    let n = u32::from(x >= x_max) + u32::from(x >> 8 >= x_max);
    let at = words.at;
    debug_assert!(at >= 2, "a word buffer holds 2 bytes a symbol");
    // Decode order: the second byte emitted (bits 8–15) first.
    words.buf[at - 2..at].copy_from_slice(&(x as u16).to_be_bytes());
    words.at = at - n as usize;
    let x = x >> (8 * n);
    let q = ((u64::from(x + x) * (e & 0xFFFF_FFFF)) >> 32) as u32 >> (info >> 25);
    x + (info & 0xFFF) + q * ((1 << 12) - f)
}

/// Order 0: byte `i` with state `i % 4`, backwards; the last, partial round
/// one state at a time, the whole rounds in the kernel.
#[allow(clippy::indexing_slicing, reason = "u8s index a 256-entry table; l < 4")]
fn code_order_0(data: &[u8], enc: &[Enc; 256], words: &mut Words<'_>, level: Level) -> [u32; 4] {
    let mut x = [LOWER_BOUND; 4];
    let (rounds, tail) = data.as_chunks::<4>();
    for (l, &b) in tail.iter().enumerate().rev() {
        x[l] = put(x[l], enc[usize::from(b)], words);
    }
    simd::code_order_0(level, rounds, enc, &mut x, words);
    x
}

/// The order-1 lookup of `b` after `ctx`.
fn pair(ctx: u8, b: u8) -> usize {
    usize::from(u16::from_be_bytes([ctx, b]))
}

/// Order 1 (`data` at least 4 bytes): position by position across the four
/// quarters, then the remainder with state 3, backwards. A quarter's first
/// byte has context 0; the remainder continues the last quarter. The
/// remainder and the quarters' first bytes one state at a time, the rest in
/// the kernel.
#[allow(clippy::indexing_slicing, reason = "u16s index a 65536-entry table; windows of 2")]
#[allow(clippy::arithmetic_side_effects, reason = "4·quarter ≥ 4")]
fn code_order_1(
    data: &[u8],
    enc: &[Enc; 1 << 16],
    words: &mut Words<'_>,
    level: Level,
) -> [u32; 4] {
    let mut x = [LOWER_BOUND; 4];
    let quarter = data.len() / 4;
    debug_assert!(quarter >= 1);
    for w in data[4 * quarter - 1..].windows(2).rev() {
        x[3] = put(x[3], enc[pair(w[0], w[1])], words);
    }
    let (q0, rest) = data.split_at(quarter);
    let (q1, rest) = rest.split_at(quarter);
    let (q2, rest) = rest.split_at(quarter);
    let q3 = &rest[..quarter];
    simd::code_order_1(level, [q0, q1, q2, q3], enc, &mut x, words);
    for (state, q) in x.iter_mut().zip([q0, q1, q2, q3]).rev() {
        *state = put(*state, enc[pair(0, q[0])], words);
    }
    x
}

/// Scale `counts` to frequencies summing to [`TOTAL`], each present symbol
/// at least 1.
fn normalize(counts: &[u32; 256]) -> [u32; 256] {
    let total: u64 = counts.iter().map(|&c| u64::from(c)).sum();
    let mut freq = [0u32; 256];
    for (f, &c) in freq.iter_mut().zip(counts) {
        if c > 0 {
            let scaled =
                u64::from(c).saturating_mul(u64::from(TOTAL)).checked_div(total).unwrap_or(0);
            *f = u32::try_from(scaled).unwrap_or(TOTAL).max(1);
        }
    }
    let mut sum: u32 = freq.iter().sum();
    // At most 256 symbols at 1 each, so the target is always reachable by
    // taking from the largest.
    while sum != TOTAL {
        let Some(max) = freq.iter_mut().max_by_key(|f| **f) else { break };
        if sum < TOTAL {
            *max = max.saturating_add(TOTAL.saturating_sub(sum));
            sum = TOTAL;
        } else {
            let take = sum.saturating_sub(TOTAL).min(max.saturating_sub(1));
            *max = max.saturating_sub(take);
            sum = sum.saturating_sub(take);
            if take == 0 {
                break;
            }
        }
    }
    freq
}

/// An order-0 frequency table: the present symbols, run-length coded as
/// htscodecs writes them, each with its frequency as ITF8; then 0.
fn write_frequencies(out: &mut Vec<u8>, freq: &[u32; 256]) {
    write_runs(out, &freq.map(|f| f > 0), |out, sym| {
        let f = freq.get(sym).copied().unwrap_or(0);
        write_itf8(out, i32::try_from(f).unwrap_or(i32::MAX));
    });
}

/// The alphabet coding both table levels share: a present symbol is written
/// unless a run covers it; one following a present symbol starts a run,
/// the count of present symbols right after it. `value` writes each one's
/// payload. The list ends with 0.
#[allow(clippy::cast_possible_truncation, reason = "symbols and run lengths are below 256")]
pub(super) fn write_runs(
    out: &mut Vec<u8>,
    present: &[bool; 256],
    mut value: impl FnMut(&mut Vec<u8>, usize),
) {
    let mut run = 0usize;
    for j in 0..256usize {
        if !present.get(j).copied().unwrap_or(false) {
            continue;
        }
        if run > 0 {
            run = run.saturating_sub(1);
        } else {
            out.push(j as u8);
            if j > 0 && present.get(j.wrapping_sub(1)).copied().unwrap_or(false) {
                run = present.iter().skip(j.saturating_add(1)).take_while(|&&p| p).count();
                out.push(run as u8);
            }
        }
        value(out, j);
    }
    out.push(0);
}
#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::arithmetic_side_effects, reason = "test code")]
mod tests {
    use super::*;
    use crate::cram::rans::decode;
    use hegel::prelude::*;

    /// A block: bytes of an alphabet of 1 to 256 symbols, uniform, in runs
    /// (repeated pairs, which the fast counts split eight ways), or skewed
    /// (one symbol and, one byte in 512, others: frequencies of 1 against
    /// one near 4095, so states climb high enough to emit two bytes), of
    /// any length, the short ones below 4 included.
    #[hegel::composite]
    fn arb_block(tc: &TestCase) -> Vec<u8> {
        let alphabet: u32 = tc.draw(gs::integers::<u32>().min_value(1).max_value(256));
        // 0: skewed.
        let max_run: u32 = tc.draw(gs::sampled_from(&[0, 1, 4, 64]));
        let len: usize = tc.draw(gs::integers::<usize>().max_value(20_000));
        let mut x = tc.draw(gs::integers::<u64>()) | 1;
        let mut next = || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x
        };
        let mut data = Vec::with_capacity(len);
        if max_run == 0 {
            #[allow(clippy::cast_possible_truncation, reason = "below alphabet ≤ 256")]
            let symbol = |r: u64| (r % u64::from(alphabet)) as u8;
            let dominant = symbol(next());
            data.extend((0..len).map(|_| {
                let r = next();
                if r % 512 == 0 { symbol(r >> 16) } else { dominant }
            }));
            return data;
        }
        while data.len() < len {
            let r = next();
            #[allow(clippy::cast_possible_truncation, reason = "below alphabet ≤ 256")]
            let b = (r % u64::from(alphabet)) as u8;
            let run = usize::try_from((r >> 32) % u64::from(max_run)).unwrap() + 1;
            data.extend(std::iter::repeat_n(b, run.min(len - data.len())));
        }
        data
    }

    // r[verify cram_writer.rans4x8.simd]
    /// The byte marking finds every value a block holds and no other,
    /// keeping the ones marked before.
    #[hegel::test]
    fn marking_finds_every_byte_value(tc: TestCase) {
        let data = tc.draw(arb_block());
        let preset: u8 = tc.draw(gs::integers::<u8>());
        let mut expected = [false; 256];
        expected[usize::from(preset)] = true;
        for &b in &data {
            expected[usize::from(b)] = true;
        }
        for level in crate::simd_levels::levels() {
            let mut present = [false; 256];
            present[usize::from(preset)] = true;
            simd::mark_present(level, &data, &mut present);
            assert_eq!(present, expected, "{level:?}");
        }
    }

    // r[verify cram_writer.rans4x8.fast]
    // r[verify cram_writer.rans4x8.simd]
    /// The fast path writes the reference loops' bytes at every SIMD level,
    /// and a block's leftovers in the thread's scratch change nothing in the
    /// next one.
    #[hegel::test]
    fn the_fast_path_writes_the_reference_bytes(tc: TestCase) {
        let first = tc.draw(arb_block());
        let second = tc.draw(arb_block());
        for data in [&first, &second, &first] {
            for order in [Order::Zero, Order::One] {
                let reference = encode_with(data, order, Coder::Reference);
                if let Some(coded) = &reference {
                    assert!(decode(coded).unwrap() == *data, "{order:?}");
                }
                for level in crate::simd_levels::levels() {
                    let fast = encode_with(data, order, Coder::Fast(level));
                    assert!(fast == reference, "{order:?} {level:?} of {} bytes", data.len());
                }
            }
        }
    }

    // r[verify cram_writer.rans4x8]
    #[hegel::test]
    fn every_input_round_trips_through_the_decoder(tc: TestCase) {
        // Small alphabets and long runs as well as uniform bytes: the
        // frequency tables' runs and the renormalization both get exercised.
        let alphabet: u8 = tc.draw(gs::integers::<u8>().min_value(1));
        let len: usize = tc.draw(gs::integers::<usize>().max_value(5000));
        let seed: u64 = tc.draw(gs::integers::<u64>());
        let mut x = seed | 1;
        let data: Vec<u8> = (0..len)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                #[allow(clippy::cast_possible_truncation, reason = "below alphabet ≤ 255")]
                let b = (x % u64::from(alphabet)) as u8;
                b
            })
            .collect();
        for order in [Order::Zero, Order::One] {
            let encoded = encode(&data, order);
            let expected_some = match order {
                Order::Zero => !data.is_empty(),
                Order::One => data.len() >= 4,
            };
            assert_eq!(encoded.is_some(), expected_some, "{order:?} of {} bytes", data.len());
            if let Some(encoded) = encoded {
                assert_eq!(decode(&encoded).unwrap(), data, "{order:?}");
            }
        }
    }

    /// Every table sums to 4095, htscodecs' total, with each present symbol
    /// at least 1; the decoder takes 4096 too, so a round trip cannot tell.
    #[hegel::test]
    fn frequencies_sum_to_htscodecs_total(tc: TestCase) {
        let mut counts = [0u32; 256];
        let n = tc.draw(gs::integers::<usize>().min_value(1).max_value(256));
        for c in counts.iter_mut().take(n) {
            *c = tc.draw(gs::integers::<u32>().max_value(1_000_000));
        }
        if counts.iter().all(|&c| c == 0) {
            counts[0] = 1;
        }
        let freq = normalize(&counts);
        assert_eq!(freq.iter().sum::<u32>(), TOTAL);
        for (&f, &c) in freq.iter().zip(&counts) {
            assert_eq!(f > 0, c > 0);
        }
    }

    // r[verify cram_writer.rans4x8]
    /// The reciprocal divides exactly: every frequency, at the states that
    /// stress it (0, the multiples of `f` and their neighbours, `2^31 − 1`)
    /// and at drawn ones.
    #[hegel::test]
    fn the_reciprocal_divides_exactly(tc: TestCase) {
        let f = tc.draw(gs::integers::<u32>().min_value(1).max_value(TOTAL));
        let sym = Symbol::new(0, f);
        let top = (1u32 << 31) - 1;
        let k = tc.draw(gs::integers::<u32>().max_value(top / f));
        let drawn = tc.draw(gs::integers::<u32>().max_value(top));
        for x in
            [0, 1, f - 1, f, k * f, (k * f).saturating_sub(1), (k * f + 1).min(top), top, drawn]
        {
            assert_eq!(sym.quotient(x), x / f, "x {x}, f {f}");
        }
    }

    #[test]
    // r[verify cram_writer.rans4x8.simd]
    /// Edge inputs round-trip, and the fast path writes the reference's
    /// bytes for them at every SIMD level.
    fn edge_inputs_round_trip() {
        let every_byte: Vec<u8> = (0..=255).cycle().take(10_000).collect();
        // One symbol at frequency 1 among 20 000 of another.
        let mut rare = vec![b'A'; 20_000];
        rare[19_999] = b'B';
        rare[7] = b'C';
        let inputs: [&[u8]; 8] = [
            b"",
            b"a",
            b"abcd",
            b"abcde",
            &[0; 1000],
            &[255, 0, 255, 0, 1, 254, 254, 254, 254],
            &every_byte,
            &rare,
        ];
        for data in inputs {
            for order in [Order::Zero, Order::One] {
                let reference = encode_with(data, order, Coder::Reference);
                if let Some(encoded) = &reference {
                    assert_eq!(decode(encoded).unwrap(), data, "{order:?} {data:?}");
                }
                for level in crate::simd_levels::levels() {
                    let fast = encode_with(data, order, Coder::Fast(level));
                    assert!(fast == reference, "{order:?} {level:?} of {} bytes", data.len());
                }
            }
        }
        for order in [Order::Zero, Order::One] {
            assert!(encode(&every_byte, order).is_some());
        }
    }

    /// A scratch already in use hands out a fresh one rather than
    /// panicking, and a word buffer past the cap is not kept. (The other
    /// fallback, a thread-local being destroyed, takes the same path.)
    #[test]
    fn scratch_in_use_or_too_large_is_not_kept() {
        let nested = with_scratch(&WORDS, |outer| {
            outer.resize(MAX_KEPT_WORDS + 1, 0);
            with_scratch(&WORDS, |inner: &mut Vec<u8>| Some(inner.capacity()))
        });
        assert_eq!(nested, Some(0), "the nested call gets a fresh buffer");
        assert_eq!(WORDS.with_borrow(Vec::capacity), 0, "an oversized buffer is dropped");
        // An encode while the order-1 scratch is held still writes the
        // reference bytes, on a fresh scratch.
        let data: Vec<u8> = (0..5000u32).map(|i| u8::try_from(i % 7).unwrap()).collect();
        let reference = encode_with(&data, Order::One, Coder::Reference);
        let nested = with_scratch(&ORDER_1, |_| Some(encode(&data, Order::One)));
        assert!(nested == Some(reference));
    }

    #[test]
    fn skewed_input_compresses() {
        let data: Vec<u8> = (0..10_000u32).map(|i| if i % 50 == 0 { b'B' } else { b'A' }).collect();
        let encoded = encode(&data, Order::Zero).unwrap();
        assert!(encoded.len() < data.len() / 4, "{} bytes", encoded.len());
    }
}
