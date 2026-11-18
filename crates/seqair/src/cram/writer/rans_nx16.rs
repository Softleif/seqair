//! rANS Nx16 encoder (CRAM 3.1, method 5): the inverse of
//! [`crate::cram::rans_nx16`], laid out as htscodecs'
//! `rans_compress_to_4x16` (`rANS_static4x16pr.c`) writes it.
#![cfg_attr(
    not(test),
    expect(dead_code, reason = "the CRAM 3.1 integration wires it into `Codecs::block`")
)]

use super::rans::Order;

// r[impl cram_writer.rans_nx16]

const FLAG_ORDER: u8 = 0x01;
const FLAG_N32: u8 = 0x04;
const FLAG_STRIPE: u8 = 0x08;
const FLAG_NO_SIZE: u8 = 0x10;
const FLAG_CAT: u8 = 0x20;
const FLAG_RLE: u8 = 0x40;
const FLAG_PACK: u8 = 0x80;

const LOWER_BOUND: u32 = 1 << 15;
const ORDER_0_BITS: u32 = 12;
/// htscodecs' default STRIPE width (`RANS_ORDER_STRIPE` with no `N`).
const STRIPES: usize = 4;
/// htscodecs drops STRIPE at or below this many bytes.
const STRIPE_MIN: usize = 20;
/// htscodecs drops order 1 below this many bytes (after PACK and RLE).
const ORDER_1_MIN: usize = 8;
/// htslib's `RANS_ORDER_SIMD_AUTO`: 32 states from this many bytes on.
const N32_MIN: usize = 50_000;

/// What to ask of [`encode`]: the stream's flags short of CAT, which
/// [`encode_cat`] writes. A transform that does not apply is dropped from
/// the stream's flags, as htscodecs drops it (`r[cram_writer.rans_nx16]`).
#[derive(Debug, Clone, Copy)]
pub(super) struct Nx16Flags {
    pub order: Order,
    /// 32 interleaved states instead of 4.
    pub n32: bool,
    /// Four byte-interleaved substreams, each coded on its own.
    pub stripe: bool,
    pub rle: bool,
    /// Two, four or eight symbols a byte, for up to 16 distinct symbols.
    pub pack: bool,
}

impl Nx16Flags {
    const fn new(order: Order, stripe: bool, rle: bool, pack: bool) -> Self {
        Self { order, n32: false, stripe, rle, pack }
    }
}

/// The combinations htslib's normal profile tries for a block of `len`
/// bytes (`cram_compress_slice` at level 5: `RANS_PR0`, `PR1`, `PR64`,
/// `PR9`, `PR128`, `PR193`), each with `RANS_ORDER_SIMD_AUTO`'s 32 states
/// from 50 000 bytes on, which htscodecs does not apply to STRIPE.
pub(super) fn normal_profile(len: usize) -> [Nx16Flags; 6] {
    let n32 = len >= N32_MIN;
    [
        Nx16Flags::new(Order::Zero, false, false, false),
        Nx16Flags::new(Order::One, false, false, false),
        Nx16Flags::new(Order::Zero, false, true, false),
        Nx16Flags::new(Order::One, true, false, false),
        Nx16Flags::new(Order::Zero, false, false, true),
        Nx16Flags::new(Order::One, false, true, true),
    ]
    .map(|f| Nx16Flags { n32: n32 && !f.stripe, ..f })
}

/// `data` as a rANS Nx16 stream, or `None` when it is too large for the
/// format's sizes (htscodecs refuses more than `i32::MAX` bytes).
pub(super) fn encode(data: &[u8], flags: Nx16Flags) -> Option<Vec<u8>> {
    check_len(data)?;
    if flags.stripe && data.len() > STRIPE_MIN {
        encode_stripe(data, flags)
    } else {
        encode_stream(data, flags, true)
    }
}

/// `data` stored as it is: a CAT stream, as a name tokeniser's tiny
/// streams use.
pub(super) fn encode_cat(data: &[u8]) -> Option<Vec<u8>> {
    check_len(data)?;
    let mut out = Vec::with_capacity(data.len().saturating_add(6));
    out.push(FLAG_CAT);
    put_uint7(&mut out, data.len())?;
    out.extend_from_slice(data);
    Some(out)
}

fn check_len(data: &[u8]) -> Option<()> {
    i32::try_from(data.len()).ok().map(|_| ())
}

/// One stream without STRIPE: PACK, then RLE, then the entropy coder, the
/// body stored as CAT when coding does not shrink it.
fn encode_stream(data: &[u8], flags: Nx16Flags, sized: bool) -> Option<Vec<u8>> {
    let n = if flags.n32 { 32 } else { 4 };
    let mut byte = if flags.n32 { FLAG_N32 } else { 0 };
    let mut out = vec![0];
    if sized {
        put_uint7(&mut out, data.len())?;
    } else {
        byte |= FLAG_NO_SIZE;
    }

    let packed = if flags.pack { pack(data) } else { None };
    let mut body = data;
    if let Some((meta, packed)) = &packed {
        byte |= FLAG_PACK;
        out.extend_from_slice(meta);
        put_uint7(&mut out, packed.len())?;
        body = packed;
    }

    let rle = if flags.rle { rle(body) } else { None };
    if let Some((meta, literals)) = &rle {
        byte |= FLAG_RLE;
        let meta_len = meta.len().checked_mul(2)?;
        match encode_order_0(meta, n) {
            Some(coded) if coded.len() < meta.len() => {
                put_uint7(&mut out, meta_len)?;
                put_uint7(&mut out, literals.len())?;
                put_uint7(&mut out, coded.len())?;
                out.extend_from_slice(&coded);
            }
            _ => {
                put_uint7(&mut out, meta_len.checked_add(1)?)?;
                put_uint7(&mut out, literals.len())?;
                out.extend_from_slice(meta);
            }
        }
        body = literals;
    }

    let order_1 = matches!(flags.order, Order::One) && body.len() >= ORDER_1_MIN;
    let coded = if order_1 { encode_order_1(body, n) } else { encode_order_0(body, n) };
    match coded {
        Some(coded) if coded.len() < body.len() => {
            if order_1 {
                byte |= FLAG_ORDER;
            }
            out.extend_from_slice(&coded);
        }
        _ => {
            byte |= FLAG_CAT;
            out.extend_from_slice(body);
        }
    }
    if let Some(first) = out.first_mut() {
        *first = byte;
    }
    Some(out)
}

/// STRIPE: the bytes dealt round-robin into four substreams, each written
/// without its size by whichever of the asked-for order, RLE, PACK and
/// order 0 is smallest, as htscodecs tries them.
fn encode_stripe(data: &[u8], flags: Nx16Flags) -> Option<Vec<u8>> {
    let order_bit = if matches!(flags.order, Order::One) { FLAG_ORDER } else { 0 };
    let n32_bit = if flags.n32 { FLAG_N32 } else { 0 };
    let mut out = vec![FLAG_STRIPE | order_bit | n32_bit];
    put_uint7(&mut out, data.len())?;
    out.push(u8::try_from(STRIPES).ok()?);

    let plain = Nx16Flags { stripe: false, rle: false, pack: false, ..flags };
    let candidates = [
        matches!(flags.order, Order::One).then_some(plain),
        flags.rle.then_some(Nx16Flags { order: Order::Zero, rle: true, ..plain }),
        flags.pack.then_some(Nx16Flags { order: Order::Zero, pack: true, ..plain }),
        Some(Nx16Flags { order: Order::Zero, ..plain }),
    ];
    let mut parts = Vec::with_capacity(STRIPES);
    for j in 0..STRIPES {
        let part: Vec<u8> = data.iter().skip(j).step_by(STRIPES).copied().collect();
        let mut best: Option<Vec<u8>> = None;
        for f in candidates.iter().flatten() {
            let coded = encode_stream(&part, *f, false)?;
            if best.as_ref().is_none_or(|b| coded.len() < b.len()) {
                best = Some(coded);
            }
        }
        let best = best?;
        put_uint7(&mut out, best.len())?;
        parts.push(best);
    }
    for part in parts {
        out.extend_from_slice(&part);
    }
    Some(out)
}

/// htscodecs' `hts_pack`: the symbol map (count, then the symbols in
/// order) and the packed bytes, each symbol's index into the map, lowest
/// bits first. `None` for more than 16 symbols.
#[allow(clippy::indexing_slicing, reason = "a u8 indexes a 256-entry table")]
#[allow(clippy::cast_possible_truncation, reason = "at most 16 symbols, so indices fit 4 bits")]
#[allow(clippy::arithmetic_side_effects, reason = "8 / per_byte with per_byte ∈ {2, 4, 8}")]
fn pack(data: &[u8]) -> Option<(Vec<u8>, Vec<u8>)> {
    let mut present = [false; 256];
    for &b in data {
        present[usize::from(b)] = true;
    }
    let mut code = [0u8; 256];
    let mut meta = vec![0];
    for (sym, _) in (0u8..=255).zip(present).filter(|(_, p)| *p) {
        code[usize::from(sym)] = (meta.len() - 1) as u8;
        meta.push(sym);
    }
    let symbols = meta.len() - 1;
    if symbols == 0 || symbols > 16 {
        return None;
    }
    if let Some(count) = meta.first_mut() {
        *count = symbols as u8;
    }
    let per_byte = match symbols {
        1 => return Some((meta, Vec::new())),
        2 => 8,
        3 | 4 => 4,
        _ => 2,
    };
    let bits = 8 / per_byte;
    let packed = data
        .chunks(per_byte)
        .map(|c| c.iter().rev().fold(0u8, |acc, &b| (acc << bits) | code[usize::from(b)]))
        .collect();
    Some((meta, packed))
}

/// htscodecs' `hts_rle_encode`: the metadata (the run symbols, as a count
/// with 256 written as 0, then the symbols, then each run's extra length
/// as uint7) and the literals. Run symbols are those that repeat more
/// often than they start (`rle_find_syms`). `None` when no symbol
/// qualifies, or when runs and literals together would not save 1%.
#[allow(clippy::indexing_slicing, reason = "a u8 indexes a 256-entry table")]
#[allow(
    clippy::cast_possible_truncation,
    reason = "256 symbols are written as 0, as htscodecs does"
)]
#[allow(clippy::arithmetic_side_effects, reason = "i + run + 1 ≤ data.len()")]
fn rle(data: &[u8]) -> Option<(Vec<u8>, Vec<u8>)> {
    let mut saved = [0i64; 256];
    let mut last = None;
    for &b in data {
        let s = &mut saved[usize::from(b)];
        if last == Some(b) {
            *s = s.saturating_add(1);
        } else {
            *s = s.saturating_sub(1);
            last = Some(b);
        }
    }
    let mut is_run = [false; 256];
    let mut meta = vec![0];
    for ((sym, &s), r) in (0u8..=255).zip(&saved).zip(&mut is_run) {
        if s > 0 {
            *r = true;
            meta.push(sym);
        }
    }
    let symbols = meta.len() - 1;
    if symbols == 0 {
        return None;
    }
    if let Some(count) = meta.first_mut() {
        *count = symbols as u8;
    }

    let mut literals = Vec::with_capacity(data.len());
    let mut i = 0;
    while let Some(&b) = data.get(i) {
        literals.push(b);
        i += 1;
        if is_run[usize::from(b)] {
            let run = data.get(i..).unwrap_or_default().iter().take_while(|&&c| c == b).count();
            put_uint7(&mut meta, run)?;
            i += run;
        }
    }
    let coded = literals.len().saturating_add(meta.len()).saturating_mul(100);
    (coded < data.len().saturating_mul(99)).then_some((meta, literals))
}

/// A symbol's slot range in a frequency table.
#[derive(Debug, Clone, Copy, Default)]
struct Symbol {
    start: u32,
    freq: u32,
}

/// Order 0, without the flags byte: the table, normalized to the next
/// power of two of the length (at most 4096) and written as that, which
/// the decoder scales up to 4096; the `n` states; the renormalization
/// words. `None` when `data` is empty (htscodecs writes no table then).
#[allow(clippy::indexing_slicing, reason = "a u8 indexes a 256-entry table; i % n < n")]
#[allow(clippy::arithmetic_side_effects, reason = "i % n with n ∈ {4, 32}")]
fn encode_order_0(data: &[u8], n: usize) -> Option<Vec<u8>> {
    if data.is_empty() {
        return None;
    }
    let mut counts = [0u32; 256];
    for &b in data {
        counts[usize::from(b)] = counts[usize::from(b)].saturating_add(1);
    }
    let total = u32::try_from(data.len().next_power_of_two()).ok()?.min(1 << ORDER_0_BITS);
    let freqs = normalize(&counts, total);
    let mut out = Vec::new();
    write_alphabet(&mut out, &freqs.map(|f| f > 0));
    for &f in freqs.iter().filter(|&&f| f > 0) {
        put_uint7(&mut out, f)?;
    }
    let table = symbols(&freqs, ORDER_0_BITS - total.trailing_zeros());

    let mut states = vec![LOWER_BOUND; n];
    let mut words = Vec::with_capacity(data.len() / 2);
    // The decoder reads byte `i` with state `i % n`; the encoder runs it
    // backwards.
    for (i, &b) in data.iter().enumerate().rev() {
        put(&mut states[i % n], &mut words, table[usize::from(b)], ORDER_0_BITS);
    }
    finish(&mut out, &states, &words);
    Some(out)
}

/// The order-1 tables as htscodecs' `encode_freq1` sizes them: each row
/// normalized to a power of two (the row total's, halved when the row has
/// under 64 symbols and the total is over 128, halved again over 1024),
/// written as that, and scaled up to `1 << shift` by the decoder.
struct Order1Tables {
    shift: u32,
    /// The alphabet: every symbol of the data, and 0, the context every
    /// state starts in. Each has a row.
    alphabet: [bool; 256],
    /// Each context's row as written.
    rows: Vec<[u32; 256]>,
}

/// Byte `p`'s context: the byte before it in its state's stretch, 0 at a
/// stretch's start. State `j` of `n` decodes `data[j·chunk..(j + 1)·chunk]`
/// with `chunk = len / n`, and the last state also the `len % n` bytes after.
#[allow(clippy::arithmetic_side_effects, reason = "p ≥ 1 in p - 1; chunk ≥ 1 in p % chunk")]
fn context(data: &[u8], chunk: usize, n: usize, p: usize) -> u8 {
    let starts_stretch = p == 0 || (p < chunk * n && p.is_multiple_of(chunk));
    if starts_stretch { 0 } else { data.get(p - 1).copied().unwrap_or(0) }
}

#[allow(clippy::indexing_slicing, reason = "u8s index 256-entry tables")]
fn order_1_tables(data: &[u8], n: usize) -> Order1Tables {
    let chunk = data.len().checked_div(n).unwrap_or(0);
    let mut counts = vec![[0u32; 256]; 256];
    let mut alphabet = [false; 256];
    alphabet[0] = true;
    for (p, &b) in data.iter().enumerate() {
        let c = &mut counts[usize::from(context(data, chunk, n, p))][usize::from(b)];
        *c = c.saturating_add(1);
        alphabet[usize::from(b)] = true;
    }
    // A symbol that never precedes another still has a row: one that only
    // ever decodes symbol 0, so both decoders take it as a complete table.
    for (row, _) in
        counts.iter_mut().zip(alphabet).filter(|(row, a)| *a && row.iter().all(|&c| c == 0))
    {
        row[0] = 1;
    }
    let targets: Vec<u32> = counts
        .iter()
        .map(|row| {
            let total: u32 = row.iter().fold(0, |a, &c| a.saturating_add(c));
            let symbols = row.iter().filter(|&&c| c > 0).count();
            let mut target = total.checked_next_power_of_two().unwrap_or(1 << 31);
            if symbols < 64 && target > 128 {
                target /= 2;
            }
            if target > 1024 {
                target /= 2;
            }
            target.min(1 << ORDER_0_BITS)
        })
        .collect();
    // htscodecs also takes 10 bits when its entropy estimate says 12 gain
    // under 1%; here only when no row needs more.
    let shift = if targets.iter().zip(alphabet).any(|(&t, a)| a && t > 1024) { 12 } else { 10 };
    let rows = counts
        .iter()
        .zip(&targets)
        .map(|(row, &t)| if row.iter().any(|&c| c > 0) { normalize(row, t) } else { [0; 256] })
        .collect();
    Order1Tables { shift, alphabet, rows }
}

/// Order 1, without the flags byte: the tables (order-0 coded when larger
/// than 1000 bytes and that saves more than 6), the `n` states, the words.
/// `data` holds at least [`ORDER_1_MIN`] bytes.
#[allow(clippy::indexing_slicing, reason = "u8s index 256-entry tables; j < n")]
#[allow(clippy::arithmetic_side_effects, reason = "j·chunk + i < n·chunk ≤ len")]
#[allow(clippy::cast_possible_truncation, reason = "shift ≤ 12 fits the high nibble")]
fn encode_order_1(data: &[u8], n: usize) -> Option<Vec<u8>> {
    let Order1Tables { shift, alphabet, rows } = order_1_tables(data, n);
    let mut table = Vec::new();
    write_alphabet(&mut table, &alphabet);
    let mut syms = vec![[Symbol::default(); 256]; 256];
    for ((row, on), s) in rows.iter().zip(alphabet).zip(syms.iter_mut()) {
        if !on {
            continue;
        }
        write_row(&mut table, &alphabet, row)?;
        let total: u32 = row.iter().sum();
        *s = symbols(row, shift - total.trailing_zeros());
    }

    let header = (shift as u8) << 4;
    let mut out = Vec::new();
    let coded = if table.len() > 1000 { encode_order_0(&table, 4) } else { None };
    match coded {
        Some(coded) if coded.len().saturating_add(6) < table.len().saturating_add(1) => {
            out.push(header | 1);
            put_uint7(&mut out, table.len())?;
            put_uint7(&mut out, coded.len())?;
            out.extend_from_slice(&coded);
        }
        _ => {
            out.push(header);
            out.extend_from_slice(&table);
        }
    }

    let chunk = data.len() / n;
    let mut states = vec![LOWER_BOUND; n];
    let mut words = Vec::with_capacity(data.len() / 2);
    let mut put_at = |state: usize, p: usize| {
        let sym = syms[usize::from(context(data, chunk, n, p))][usize::from(data[p])];
        put(&mut states[state], &mut words, sym, shift);
    };
    // Decode order: position by position across the n stretches, then the
    // remainder with the last state. Encode it backwards.
    for p in (chunk * n..data.len()).rev() {
        put_at(n - 1, p);
    }
    for i in (0..chunk).rev() {
        for j in (0..n).rev() {
            put_at(j, j * chunk + i);
        }
    }
    finish(&mut out, &states, &words);
    Some(out)
}

/// One order-1 row over the alphabet: each frequency as uint7, a run of
/// zeros as 0 and the count of zeros after the first (`encode_freq_d`).
fn write_row(out: &mut Vec<u8>, alphabet: &[bool; 256], row: &[u32; 256]) -> Option<()> {
    let mut zeros = 0u32;
    let flush = |out: &mut Vec<u8>, zeros: &mut u32| -> Option<()> {
        if *zeros > 0 {
            out.push(0);
            out.push(u8::try_from(zeros.saturating_sub(1)).ok()?);
            *zeros = 0;
        }
        Some(())
    };
    for (&f, _) in row.iter().zip(alphabet).filter(|(_, a)| **a) {
        if f == 0 {
            zeros = zeros.saturating_add(1);
        } else {
            flush(out, &mut zeros)?;
            put_uint7(out, f)?;
        }
    }
    flush(out, &mut zeros)
}

/// Scale `counts` to frequencies summing to `total`, each present symbol
/// at least 1, rounding to nearest and settling the difference on the
/// largest. `total` is at least the number of present symbols.
#[allow(
    clippy::arithmetic_side_effects,
    reason = "each subtraction is guarded by the comparison before it"
)]
fn normalize(counts: &[u32; 256], total: u32) -> [u32; 256] {
    let sum: u64 = counts.iter().map(|&c| u64::from(c)).sum();
    let mut freq = [0u32; 256];
    for (f, &c) in freq.iter_mut().zip(counts) {
        if c > 0 {
            let scaled = u64::from(c)
                .saturating_mul(u64::from(total))
                .saturating_add(sum / 2)
                .checked_div(sum)
                .unwrap_or(0);
            *f = u32::try_from(scaled).unwrap_or(total).max(1);
        }
    }
    let mut got: u32 = freq.iter().fold(0, |a, &f| a.saturating_add(f));
    while got != total {
        let Some(max) = freq.iter_mut().max_by_key(|f| **f) else { break };
        if got < total {
            *max = max.saturating_add(total - got);
            got = total;
        } else {
            let take = (got - total).min(max.saturating_sub(1));
            if take == 0 {
                break;
            }
            *max -= take;
            got -= take;
        }
    }
    debug_assert_eq!(got, total, "a total below the symbol count cannot be reached");
    freq
}

/// Each symbol's slot range once its frequency is scaled up by `scale` bits.
fn symbols(freqs: &[u32; 256], scale: u32) -> [Symbol; 256] {
    let mut table = [Symbol::default(); 256];
    let mut start = 0u32;
    for (s, &f) in table.iter_mut().zip(freqs) {
        let freq = f << scale;
        *s = Symbol { start, freq };
        start = start.saturating_add(freq);
    }
    table
}

/// `RansEncPutSymbol` for 16-bit renormalization: emit the low word if the
/// state would leave `[L, L << 16)`, then encode.
#[allow(clippy::arithmetic_side_effects, reason = "freq ≥ 1; the state stays below 2^31")]
#[allow(clippy::cast_possible_truncation, reason = "the low word")]
fn put(x: &mut u32, words: &mut Vec<u16>, sym: Symbol, bits: u32) {
    debug_assert!(sym.freq > 0, "encoding a symbol without frequency");
    let freq = sym.freq.max(1);
    if *x >= ((LOWER_BOUND >> bits) << 16).saturating_mul(freq) {
        words.push(*x as u16);
        *x >>= 16;
    }
    *x = ((*x / freq) << bits) + (*x % freq) + sym.start;
}

/// The states, then the words in decode order (they were emitted backwards).
fn finish(out: &mut Vec<u8>, states: &[u32], words: &[u16]) {
    out.reserve(states.len().saturating_mul(4).saturating_add(words.len().saturating_mul(2)));
    for s in states {
        out.extend_from_slice(&s.to_le_bytes());
    }
    for w in words.iter().rev() {
        out.extend_from_slice(&w.to_le_bytes());
    }
}

/// The present symbols, run-length coded as htscodecs' `encode_alphabet`
/// writes them: a symbol following a present one starts a run, the count
/// of present symbols right after it; then 0.
#[allow(clippy::cast_possible_truncation, reason = "symbols and run lengths are below 256")]
fn write_alphabet(out: &mut Vec<u8>, present: &[bool; 256]) {
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
    }
    out.push(0);
}

/// uint7: 7 bits a byte, most significant first, the high bit set on all
/// but the last.
#[allow(clippy::cast_possible_truncation, reason = "masked to 7 bits")]
#[allow(clippy::arithmetic_side_effects, reason = "g < 5")]
fn put_uint7(out: &mut Vec<u8>, v: impl TryInto<u32>) -> Option<()> {
    let v: u32 = v.try_into().ok()?;
    let groups = (0..5u32).rev().skip_while(|&g| g > 0 && v >> (7 * g) == 0);
    for g in groups {
        let bits = ((v >> (7 * g)) & 0x7f) as u8;
        out.push(if g == 0 { bits } else { bits | 0x80 });
    }
    Some(())
}

#[cfg(test)]
#[allow(
    clippy::unwrap_used,
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::indexing_slicing,
    reason = "test code"
)]
mod tests {
    use super::*;
    use crate::cram::rans_nx16::decode;
    use hegel::prelude::*;

    // r[verify cram_writer.rans_nx16]

    /// Every combination [`encode`] takes: order × states × STRIPE × RLE × PACK.
    fn every_flags() -> impl Iterator<Item = Nx16Flags> {
        (0u8..32).map(|m| Nx16Flags {
            order: if m & 1 == 0 { Order::Zero } else { Order::One },
            n32: m & 2 != 0,
            stripe: m & 4 != 0,
            rle: m & 8 != 0,
            pack: m & 16 != 0,
        })
    }

    fn flags_byte(f: Nx16Flags) -> u8 {
        [
            (matches!(f.order, Order::One), FLAG_ORDER),
            (f.n32, FLAG_N32),
            (f.stripe, FLAG_STRIPE),
            (f.rle, FLAG_RLE),
            (f.pack, FLAG_PACK),
        ]
        .iter()
        .filter(|(on, _)| *on)
        .fold(0, |byte, (_, bit)| byte | bit)
    }

    fn round_trip(data: &[u8], flags: Nx16Flags) -> Vec<u8> {
        let encoded = encode(data, flags).unwrap();
        let decoded = decode(&encoded, data.len())
            .unwrap_or_else(|e| panic!("{flags:?} of {} bytes: {e}", data.len()));
        assert!(decoded == data, "{flags:?} of {} bytes: {:#04x}", data.len(), encoded[0]);
        encoded
    }

    /// Bytes from a small, skewed alphabet in runs: what PACK and RLE need,
    /// and order 1's contexts with few successors. xorshift, so that hegel
    /// draws a few numbers rather than a vector it would keep short.
    fn draw_data(tc: &TestCase) -> Vec<u8> {
        let alphabet: u16 = tc.draw(gs::integers::<u16>().min_value(1).max_value(256));
        let base: u8 = tc.draw(gs::integers::<u8>());
        let skew: u32 = tc.draw(gs::integers::<u32>().max_value(8));
        let max_run: u32 = tc.draw(gs::integers::<u32>().min_value(1).max_value(64));
        let len: usize = tc.draw(gs::integers::<usize>().max_value(5000));
        let mut x = tc.draw(gs::integers::<u64>()) | 1;
        let mut next = move || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x
        };
        let mut data = Vec::with_capacity(len);
        while data.len() < len {
            // The least of `skew + 1` draws favours the low symbols.
            let sym = (0..=skew).map(|_| next() % u64::from(alphabet)).min().unwrap_or(0);
            let run = (next() % u64::from(max_run)) as usize + 1;
            let byte = base.wrapping_add(sym as u8);
            data.extend(std::iter::repeat_n(byte, run.min(len - data.len())));
        }
        data
    }

    #[hegel::test(test_cases = crate::pinned::cases(60))]
    fn every_input_round_trips_under_every_combination(tc: TestCase) {
        let data = draw_data(&tc);
        for flags in every_flags() {
            round_trip(&data, flags);
        }
    }

    #[hegel::test]
    fn every_order_1_row_is_a_complete_power_of_two_table(tc: TestCase) {
        let data = draw_data(&tc);
        let n = if tc.draw(gs::booleans()) { 32 } else { 4 };
        let Order1Tables { shift, alphabet, rows } = order_1_tables(&data, n);
        assert!(shift == 10 || shift == 12, "htscodecs builds tables for 10 and 12 bits only");
        for (sym, on) in alphabet.iter().enumerate() {
            assert_eq!(*on, sym == 0 || data.contains(&(sym as u8)), "symbol {sym}");
        }
        for (ctx, row) in rows.iter().enumerate().filter(|(c, _)| alphabet[*c]) {
            let total: u32 = row.iter().sum();
            assert!(total.is_power_of_two() && total <= 1 << shift, "context {ctx}: {total}");
            for (sym, &f) in row.iter().enumerate() {
                assert!(f == 0 || alphabet[sym], "context {ctx} codes {sym}, outside the alphabet");
            }
        }
        // Each byte's context: the byte before it, or 0 where a state's
        // stretch starts — counted independently of `context`.
        let chunk = data.len() / n;
        let mut starts: Vec<usize> = (0..n).map(|j| j * chunk).collect();
        starts.push(0);
        for (p, &b) in data.iter().enumerate() {
            let ctx = if starts.contains(&p) && p < n * chunk.max(1) { 0 } else { data[p - 1] };
            assert!(rows[usize::from(ctx)][usize::from(b)] > 0, "byte {p} in context {ctx}");
        }
    }

    #[hegel::test]
    fn normalized_frequencies_sum_to_the_total(tc: TestCase) {
        let mut counts = [0u32; 256];
        let n = tc.draw(gs::integers::<usize>().min_value(1).max_value(256));
        for c in counts.iter_mut().take(n) {
            *c = tc.draw(gs::integers::<u32>().max_value(1_000_000));
        }
        counts[n - 1] = counts[n - 1].max(1);
        let present = counts.iter().filter(|&&c| c > 0).count() as u32;
        let shift = tc.draw(gs::integers::<u32>().max_value(12));
        let total = (1u32 << shift).max(present.next_power_of_two());
        let freqs = normalize(&counts, total);
        assert_eq!(freqs.iter().sum::<u32>(), total);
        for (f, c) in freqs.iter().zip(&counts) {
            assert_eq!(*f > 0, *c > 0);
        }
    }

    #[test]
    fn edge_inputs_round_trip() {
        let every_byte: Vec<u8> = (0..=255).cycle().take(10_000).collect();
        let long_runs: Vec<u8> = (0..60_000u32).map(|i| (i / 997 % 3) as u8).collect();
        let inputs: [&[u8]; 10] = [
            b"",
            b"a",
            b"abcdefg",
            b"abcdefgh",
            b"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
            &[0; 1000],
            &[255, 0, 255, 0, 1, 254, 254, 254, 254],
            &every_byte,
            &long_runs,
            &[7; 70_000],
        ];
        for data in inputs {
            for flags in every_flags().chain(normal_profile(data.len())) {
                round_trip(data, flags);
            }
        }
    }

    #[test]
    fn cat_round_trips() {
        for data in [&b""[..], b"x", &[3; 300]] {
            let encoded = encode_cat(data).unwrap();
            assert_eq!(encoded[0], FLAG_CAT);
            assert_eq!(decode(&encoded, data.len()).unwrap(), data);
        }
    }

    #[test]
    fn transforms_that_do_not_apply_are_dropped() {
        let o1 = |stripe, rle, pack| Nx16Flags::new(Order::One, stripe, rle, pack);
        // 17 symbols: no PACK; no symbol repeats: no RLE.
        let distinct: Vec<u8> = (0..17).cycle().take(1000).collect();
        assert_eq!(round_trip(&distinct, o1(false, true, true))[0], FLAG_ORDER);
        // Under 8 bytes: order 0, which cannot beat 7 bytes, so CAT.
        assert_eq!(round_trip(b"abcabca", o1(false, false, false))[0], FLAG_CAT);
        // One symbol: PACK leaves nothing to code.
        assert_eq!(round_trip(&[9; 100], o1(false, false, true))[0], FLAG_PACK | FLAG_CAT);
        let runs: Vec<u8> = (0..4000u32).map(|i| (i / 40 % 2) as u8).collect();
        let both = FLAG_PACK | FLAG_RLE;
        assert_eq!(round_trip(&runs, o1(false, true, true))[0] & both, both);
        assert_eq!(round_trip(&distinct, o1(true, false, false))[0], FLAG_STRIPE | FLAG_ORDER);
        assert_eq!(round_trip(&distinct[..20], o1(true, false, false))[0] & FLAG_STRIPE, 0);
    }

    #[test]
    fn normal_profile_is_htslibs_level_5() {
        assert_eq!(normal_profile(1000).map(flags_byte), [0, 1, 64, 9, 128, 193]);
        assert_eq!(normal_profile(50_000).map(flags_byte), [4, 5, 68, 9, 132, 197]);
    }

    /// hts-specs' ransNx16 vectors are htscodecs' output under the flags
    /// their suffix names: ours under the same flags round-trips and comes
    /// within 1% of htscodecs' size.
    #[test]
    fn hts_specs_vectors_compress_like_htscodecs() {
        let dir =
            concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/hts-specs/cram/codecs/ransNx16");
        let mut seen = 0;
        for entry in std::fs::read_dir(dir).unwrap() {
            let path = entry.unwrap().path();
            let theirs = std::fs::read(&path).unwrap();
            let suffix: u8 = path.extension().unwrap().to_str().unwrap().parse().unwrap();
            let original = decode(&theirs, 0).unwrap();
            let flags = Nx16Flags {
                order: if suffix & FLAG_ORDER != 0 { Order::One } else { Order::Zero },
                n32: suffix & FLAG_N32 != 0,
                stripe: suffix & FLAG_STRIPE != 0,
                rle: suffix & FLAG_RLE != 0,
                pack: suffix & FLAG_PACK != 0,
            };
            assert_eq!(flags_byte(flags), suffix, "{}", path.display());
            let ours = round_trip(&original, flags);
            let ratio = ours.len() as f64 / theirs.len() as f64;
            assert!(
                ratio < 1.01,
                "{}: {} bytes against htscodecs' {}",
                path.display(),
                ours.len(),
                theirs.len()
            );
            seen += 1;
        }
        assert!(seen >= 31, "expected the hts-specs ransNx16 vectors, found {seen}");
    }

    /// htscodecs' own decoders as the oracle: `HTSCODECS_RANS4X16PR` names
    /// its `tests/rANS_static4x16pr_test.c` built as a program
    /// (`-d -r in out` decodes a stream). Each stream is decoded with the
    /// default decoders (SIMD where the CPU has it) and with `-c 0`, the
    /// scalar ones. Skipped when the variable is unset.
    #[hegel::test(test_cases = crate::pinned::cases(20))]
    fn htscodecs_decodes_every_combination(tc: TestCase) {
        htscodecs_decodes(&draw_data(&tc));
    }

    /// [`htscodecs_decodes_every_combination`] on blocks large enough for
    /// the SIMD decoders' main loops: real qualities and long runs.
    #[test]
    fn htscodecs_decodes_large_blocks() {
        let q8 = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../tests/data/hts-specs/cram/codecs/ransNx16/q8.0"
        );
        let qualities = decode(&std::fs::read(q8).unwrap(), 0).unwrap();
        let long_runs: Vec<u8> = (0..80_000u32).map(|i| (i / 997 % 5) as u8).collect();
        for data in [qualities, long_runs] {
            htscodecs_decodes(&data);
        }
    }

    fn htscodecs_decodes(data: &[u8]) {
        let Some(program) = std::env::var_os("HTSCODECS_RANS4X16PR") else {
            return;
        };
        let dir = tempfile::tempdir().unwrap();
        let (src, dst) = (dir.path().join("in"), dir.path().join("out"));
        for flags in every_flags() {
            std::fs::write(&src, round_trip(data, flags)).unwrap();
            for cpu in [None, Some("0")] {
                let mut cmd = std::process::Command::new(&program);
                if let Some(cpu) = cpu {
                    cmd.args(["-c", cpu]);
                }
                let status = cmd.arg("-d").arg("-r").arg(&src).arg(&dst).output().unwrap();
                assert!(status.status.success(), "{flags:?} cpu {cpu:?}: {status:?}");
                let decoded = std::fs::read(&dst).unwrap();
                assert!(decoded == data, "{flags:?} cpu {cpu:?} of {} bytes", data.len());
            }
        }
    }
}
