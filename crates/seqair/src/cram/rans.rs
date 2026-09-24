//! rANS 4x8 codec (CRAM compression method 4).
//!
//! 4-way interleaved asymmetric numeral systems with 8-bit renormalization.
//! Supports order-0 and order-1 modes.

// `clippy::unnecessary_lazy_evaluations` flags every `ok_or_else(||
// CramError::Truncated { ... })` here, because the variant's
// construction is "cheap". On the per-byte hot path it isn't: the
// eager `ok_or` form has the compiler build (and drop) a full-sized
// `CramError` on every successful read — the type's Drop has to
// dispatch on the discriminant because other variants own `PathBuf` /
// `std::io::Error` / `SmolStr`. samply showed
// `core::ptr::drop_in_place<CramError>` next to `decode_order_1` /
// `read_u8` for exactly this. Keep the lazy form here.
#![allow(
    clippy::unnecessary_lazy_evaluations,
    reason = "lazy form avoids per-call drop_in_place<CramError> on hot path"
)]

use super::reader::CramError;

const ALPHABET_SIZE: usize = 256;
const LOWER_BOUND: u32 = 1 << 23;

/// Reusable allocations for rANS 4x8 order-1 decoding.
///
/// Order-1 keeps two tables per context: the slot → symbol table (4 KiB)
/// and the symbol → step table (1 KiB, see [`step_entry`]). Allocating
/// them fresh for every block dominated CRAM decode profiles, so the
/// reader keeps one set; a block rebuilds only the rows it can reach.
pub(crate) struct Rans4x8Buf {
    freq: Box<[[u16; ALPHABET_SIZE]; ALPHABET_SIZE]>,
    sym_tables: Box<[SymRow; ALPHABET_SIZE]>,
    steps: Box<[StepRow; ALPHABET_SIZE]>,
}

type SymRow = [u8; 4096];
type StepRow = [u32; ALPHABET_SIZE];

/// A zeroed `Box<[T; N]>` straight from the allocator, never staged on the stack.
fn zeroed_box<T: Clone, const N: usize>(zero: T) -> Box<[T; N]> {
    let Ok(b) = vec![zero; N].into_boxed_slice().try_into() else {
        unreachable!("a Vec of N elements converts to [T; N]")
    };
    b
}

impl Rans4x8Buf {
    pub fn new() -> Self {
        Self {
            freq: zeroed_box([0u16; ALPHABET_SIZE]),
            sym_tables: zeroed_box([0u8; 4096]),
            steps: zeroed_box([0u32; ALPHABET_SIZE]),
        }
    }
}

impl Default for Rans4x8Buf {
    fn default() -> Self {
        Self::new()
    }
}

// r[impl cram.codec.rans4x8]
/// Decode a rANS 4x8 compressed block (allocating path — for callers
/// without a reusable buffer). Prefer `decode_with_buf` on the hot path.
pub fn decode(src: &[u8]) -> Result<Vec<u8>, CramError> {
    let mut buf = Rans4x8Buf::new();
    decode_with_buf(src, &mut buf)
}

/// Decode a rANS 4x8 compressed block, reusing the caller-owned
/// [`Rans4x8Buf`] for order-1 tables. Order-0 blocks don't touch the buf.
pub(crate) fn decode_with_buf(src: &[u8], buf: &mut Rans4x8Buf) -> Result<Vec<u8>, CramError> {
    let mut cur: &[u8] = src;

    // Header: order (u8), compressed_size (u32 LE), uncompressed_size (u32 LE).
    // The 9-byte header is followed by `compressed_size` payload bytes; htslib
    // rejects mismatch with `in_sz != in_size - 9` (rANS_static.c:245).
    // r[impl cram.codec.rans4x8_compressed_size_check]
    let order = read_u8(&mut cur).ok_or_else(|| CramError::Truncated { context: "rans header" })?;
    let compressed_size =
        read_u32_le(&mut cur).ok_or_else(|| CramError::Truncated { context: "rans header" })?;
    let uncompressed_size = read_u32_le(&mut cur)
        .ok_or_else(|| CramError::Truncated { context: "rans header" })?
        as usize;

    if cur.len() != compressed_size as usize {
        return Err(CramError::Rans4x8CompressedSizeMismatch {
            advertised: compressed_size,
            actual: cur.len(),
        });
    }

    super::reader::check_alloc_size(uncompressed_size, "rANS 4x8 output")?;
    let mut dst = vec![0u8; uncompressed_size];

    match order {
        0 => decode_order_0_fast(&mut cur, &mut dst)?,
        1 => decode_order_1_fast(&mut cur, &mut dst, buf)?,
        _ => return Err(CramError::InvalidRansOrder { order }),
    }

    Ok(dst)
}

/// The reference order-0 decoder the fast one is tested against.
#[cfg(test)]
#[allow(
    clippy::indexing_slicing,
    reason = "indices are bounded: f ≤ 4095 (12-bit mask), sym/i ≤ 255 (u8)"
)]
fn decode_order_0(src: &mut &[u8], dst: &mut [u8]) -> Result<(), CramError> {
    // CramError::Truncated is materialised lazily on the err path; the
    // hot inner loop never builds one because `renormalize` returns
    // `Option<()>`. See the helper-module comment for the cost story.
    let truncated = || CramError::Truncated { context: "rans order-0 truncated" };
    let freq = read_frequencies_0(src)?;
    let cum_freq = build_cumulative_frequencies(&freq);
    let sym_table = build_symbol_table(&cum_freq);
    let mut states = read_states(src).ok_or_else(truncated)?;

    for chunk in dst.chunks_mut(4) {
        for (d, state) in chunk.iter_mut().zip(states.iter_mut()) {
            let f = (*state & 0x0FFF) as u16;
            let sym = sym_table[f as usize];
            *d = sym;
            let i = sym as usize;
            // r[depends cram.codec.state_step_safety]
            *state = u32::from(freq[i])
                .wrapping_mul(*state >> 12)
                .wrapping_add(*state & 0x0FFF)
                .wrapping_sub(u32::from(cum_freq[i]));
            renormalize(state, src).ok_or_else(truncated)?;
        }
    }

    Ok(())
}

#[cfg(test)]
struct Tables {
    freq: Box<[[u16; ALPHABET_SIZE]; ALPHABET_SIZE]>,
    cum_freq: Box<[[u16; ALPHABET_SIZE]; ALPHABET_SIZE]>,
    sym_tables: Box<[[u8; 4096]; ALPHABET_SIZE]>,
}

/// The reference order-1 decoder the fast one is tested against: every
/// row built, inactive ones from all-zero frequencies.
#[cfg(test)]
#[allow(
    clippy::indexing_slicing,
    reason = "indices are bounded: ctx/sym ≤ 255 (u8), f ≤ 4095 (12-bit mask)"
)]
fn decode_order_1_buf(src: &mut &[u8], dst: &mut [u8]) -> Result<(), CramError> {
    let truncated = || CramError::Truncated { context: "rans order-1 truncated" };
    let mut freq = zeroed_box([0u16; ALPHABET_SIZE]);
    read_frequencies_1_into(src, &mut freq)?;
    let mut cum_freq = zeroed_box([0u16; ALPHABET_SIZE]);
    let mut sym_tables = zeroed_box([0u8; 4096]);
    for ((f, cum), table) in freq.iter().zip(cum_freq.iter_mut()).zip(sym_tables.iter_mut()) {
        *cum = build_cumulative_frequencies(f);
        *table = build_symbol_table(cum);
    }
    let buf = Tables { freq, cum_freq, sym_tables };

    let mut states = read_states(src).ok_or_else(truncated)?;
    let mut prev_syms = [0u8; 4];

    // Chunk-based interleaving: split output into 4 equal segments.
    // Each state decodes into its own segment sequentially.
    let chunk_size = dst.len() / 4;
    let bases: [usize; 4] = [0, chunk_size, chunk_size.wrapping_mul(2), chunk_size.wrapping_mul(3)];

    for pos in 0..chunk_size {
        for si in 0usize..4 {
            let mut state = states[si];
            let ctx = prev_syms[si] as usize;
            let f = (state & 0x0FFF) as u16;
            let sym = buf.sym_tables[ctx][f as usize];
            dst[bases[si].wrapping_add(pos)] = sym;
            let sym_idx = sym as usize;
            state = u32::from(buf.freq[ctx][sym_idx])
                .wrapping_mul(state >> 12)
                .wrapping_add(state & 0x0FFF)
                .wrapping_sub(u32::from(buf.cum_freq[ctx][sym_idx]));
            renormalize(&mut state, src).ok_or_else(truncated)?;
            states[si] = state;
            prev_syms[si] = sym;
        }
    }

    let remainder_start = chunk_size.wrapping_mul(4);
    for pos in remainder_start..dst.len() {
        let ctx = prev_syms[3] as usize;
        let f = (states[3] & 0x0FFF) as u16;
        let sym = buf.sym_tables[ctx][f as usize];
        dst[pos] = sym;
        let sym_idx = sym as usize;
        states[3] = u32::from(buf.freq[ctx][sym_idx])
            .wrapping_mul(states[3] >> 12)
            .wrapping_add(states[3] & 0x0FFF)
            .wrapping_sub(u32::from(buf.cum_freq[ctx][sym_idx]));
        if pos.checked_add(1).is_some_and(|next| next < dst.len()) {
            renormalize(&mut states[3], src).ok_or_else(truncated)?;
        }
        prev_syms[3] = sym;
    }

    Ok(())
}

// ── Fast path ────────────────────────────────────────────────────────
//
// htscodecs' layout: per context a slot → symbol table (u8) and a symbol →
// step table, 5 KiB a context where a packed slot table takes 16 KiB, so a
// quality stream's few dozen contexts stay in L1/L2. The step is the
// reference decoder's own arithmetic, so any table decodes exactly as there.
// What the fast loop adds is a branch-free renormalization of 0, 1 or 2
// bytes, all four states' taken from one big-endian load. That is exact
// when a step leaves x ≥ 2^11 (x ≥ 2^23 before, and the symbol's frequency
// ≥ 1); a round where some state falls below runs checked instead.

/// The unpacked decoder's step inputs for one symbol, `cum << 16 | freq`.
fn step_entry(freq: u16, cum: u16) -> u32 {
    (u32::from(cum) << 16) | u32::from(freq)
}

/// `freq · (x >> 12) + (x & 0xFFF) − cum`, the reference decoder's step.
#[inline(always)]
fn step(x: u32, entry: u32) -> u32 {
    (entry & 0xFFFF).wrapping_mul(x >> 12).wrapping_add(x & 0x0FFF).wrapping_sub(entry >> 16)
}

/// The rows a context's frequencies give the reference decoder: the scan's
/// slot → symbol table, and [`step_entry`] per symbol.
fn build_row(freq: &[u16; ALPHABET_SIZE], syms: &mut SymRow, steps: &mut StepRow) {
    let cum = build_cumulative_frequencies(freq);
    fill_symbol_table(freq, &cum, syms);
    for ((e, &f), &c) in steps.iter_mut().zip(freq).zip(&cum) {
        *e = step_entry(f, c);
    }
}

/// [`build_symbol_table`] without the per-slot scan when the frequencies
/// total at most 4096: each symbol then owns `cum..cum + freq`, and the
/// slots past the total fall through to 255, as the scan finds.
#[allow(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "cum + freq ≤ total ≤ 4096 is checked first"
)]
fn fill_symbol_table(freq: &[u16; ALPHABET_SIZE], cum: &[u16; ALPHABET_SIZE], syms: &mut SymRow) {
    let total: u32 = freq.iter().map(|&f| u32::from(f)).sum();
    if total > 4096 {
        *syms = build_symbol_table(cum);
        return;
    }
    let mut sym = 0u8;
    for (&f, &c) in freq.iter().zip(cum) {
        let c = usize::from(c);
        syms[c..c + usize::from(f)].fill(sym);
        sym = sym.wrapping_add(1);
    }
    syms[total as usize..].fill(255);
}

/// Renormalize four stepped states (each ≥ 2^11) from the next 8 bytes,
/// in state order; returns how many bytes they took.
#[inline(always)]
#[allow(
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    reason = "n ≤ 2 per state; the top 16 bits of a u64 fit a u32"
)]
fn renorm4(xs: &mut [u32; 4], next8: &[u8; 8]) -> usize {
    let mut w = u64::from_be_bytes(*next8);
    let mut used = 0;
    for x in xs {
        let n = u32::from(*x < LOWER_BOUND) + u32::from(*x < 1 << 15);
        let shift = 8 * n;
        *x = (*x << shift) | ((w >> 48) as u32 >> (16 - shift));
        w <<= shift;
        used += n as usize;
    }
    used
}

/// The next 8 stream bytes at `pos`, if there are 8.
#[inline(always)]
fn next8(bytes: &[u8], pos: usize) -> Option<&[u8; 8]> {
    bytes.get(pos..)?.first_chunk::<8>()
}

// ── Order 0: split slot tables ───────────────────────────────────────
//
// htscodecs' `rans_uncompress_O0` layout: per slot the symbol, its
// frequency and the slot's offset into the symbol's range, in three
// tables, so a step is `f·(x >> 12) + bias` with all three loads
// independent and no field extraction. Built from the reference decoder's
// scan, the tables give its step bit for bit, whatever the total.

/// Slots per order-0 table.
const SLOTS: usize = 4096;

/// The order-0 slot tables, and the first slot whose step may leave a
/// state below 2^11 (`SLOTS` if none).
struct Order0Tables {
    sym: [u8; SLOTS],
    freq: [u16; SLOTS],
    bias: [u16; SLOTS],
    fast_limit: usize,
}

impl Order0Tables {
    /// The reference decoder's slot → (symbol, frequency, `slot − cum`).
    /// Each symbol owns `cum..cum + freq`; the slots past the total go to
    /// symbol 255, as the scan finds, and so does every slot past a total
    /// over 4096, which takes the scan itself.
    #[allow(
        clippy::indexing_slicing,
        clippy::arithmetic_side_effects,
        clippy::cast_possible_truncation,
        reason = "cum + freq ≤ total ≤ SLOTS is checked first; slot − cum[sym] ∈ [0, 4096) \
                  because the scan only moves past a symbol whose range ends at or before the slot"
    )]
    fn new(freq: &[u16; ALPHABET_SIZE]) -> Self {
        let mut t = Self { sym: [0; SLOTS], freq: [0; SLOTS], bias: [0; SLOTS], fast_limit: SLOTS };
        let cum = build_cumulative_frequencies(freq);
        let total: u32 = freq.iter().map(|&f| u32::from(f)).sum();
        if total > SLOTS as u32 {
            t.sym = build_symbol_table(&cum);
            for (m, &s) in t.sym.iter().enumerate() {
                t.freq[m] = freq[usize::from(s)];
                t.bias[m] = (m as u16).wrapping_sub(cum[usize::from(s)]);
            }
            // Frequencies over 4096 can wrap the step: nothing is fast.
            t.fast_limit = 0;
            return t;
        }
        let mut s = 0u8;
        for (&f, &c) in freq.iter().zip(&cum) {
            let c = usize::from(c);
            let range = c..c + usize::from(f);
            t.sym[range.clone()].fill(s);
            t.freq[range.clone()].fill(f);
            for (b, y) in t.bias[range].iter_mut().zip(0u16..) {
                *b = y;
            }
            s = s.wrapping_add(1);
        }
        let total = total as usize;
        let (f255, c255) = (freq[255], usize::from(cum[255]));
        t.sym[total..].fill(255);
        t.freq[total..].fill(f255);
        for (m, b) in t.bias.iter_mut().enumerate().skip(total) {
            *b = (m - c255) as u16;
        }
        // A zero frequency leaves `bias < 2^12` alone, which may need more
        // than two renorm bytes.
        if f255 == 0 {
            t.fast_limit = total;
        }
        t
    }

    /// The reference decoder's step at slot `m = x & 0xFFF`.
    #[inline(always)]
    #[allow(clippy::indexing_slicing, reason = "m < SLOTS")]
    fn step_at(&self, x: u32, m: usize) -> u32 {
        u32::from(self.freq[m]).wrapping_mul(x >> 12).wrapping_add(u32::from(self.bias[m]))
    }

    /// The reference decoder's symbol and step at slot `x & 0xFFF`.
    #[inline(always)]
    #[allow(clippy::indexing_slicing, reason = "x & 0xFFF < SLOTS")]
    fn step(&self, x: u32) -> (u8, u32) {
        let m = (x & 0xFFF) as usize;
        (self.sym[m], self.step_at(x, m))
    }
}

/// htscodecs' `RansDecRenorm`: the first renorm byte taken branch-free (a
/// cmov and a carry into the cursor), the rare second one behind a branch.
/// Reads `w[k]`, `w[k + 1]` at most; the caller guarantees `x ≥ 2^7`, so
/// two bytes lift it to 2^23, and that `k + 2 ≤ 8` for the reads it makes
/// (the mask keeps the index in the window, and is exact under that bound).
#[inline(always)]
#[allow(clippy::indexing_slicing, reason = "k & 7 < 8")]
fn renorm_cmov(x: u32, w: &[u8; 8], k: &mut usize) -> u32 {
    let low = x < LOWER_BOUND;
    let y = (x << 8) | u32::from(w[*k & 7]);
    let x = core::hint::select_unpredictable(low, y, x);
    *k = k.wrapping_add(usize::from(low));
    if x < LOWER_BOUND {
        let x = (x << 8) | u32::from(w[*k & 7]);
        *k = k.wrapping_add(1);
        x
    } else {
        x
    }
}

// r[impl cram.codec.rans4x8_packed+2]
fn decode_order_0_fast(src: &mut &[u8], dst: &mut [u8]) -> Result<(), CramError> {
    let truncated = || CramError::Truncated { context: "rans order-0 truncated" };
    let freq = read_frequencies_0(src)?;
    let t = Order0Tables::new(&freq);
    let mut states = read_states(src).ok_or_else(truncated)?;
    let bytes = *src;
    let (chunks, remainder) = dst.as_chunks_mut::<4>();
    let (done, pos) = if states.iter().all(|&x| x >= LOWER_BOUND) {
        order_0_rounds(&t, bytes, &mut states, chunks)
    } else {
        (0, 0)
    };

    let mut rest = bytes.get(pos..).ok_or_else(truncated)?;
    let tail = chunks.get_mut(done..).unwrap_or_default().iter_mut();
    let tail = tail.flat_map(|c| c.iter_mut().zip(0usize..));
    for (d, j) in tail.chain(remainder.iter_mut().zip(0usize..)) {
        let x = states.get_mut(j).ok_or_else(truncated)?;
        let (sym, next) = t.step(*x);
        *d = sym;
        *x = next;
        renormalize(x, &mut rest).ok_or_else(truncated)?;
    }
    *src = rest;
    Ok(())
}

/// The fast order-0 rounds, four symbols each, from states all ≥ 2^23:
/// one 8-byte window per round (each state takes at most two bytes, a
/// step at a slot below `fast_limit` leaving x ≥ 2^11), so the stream needs
/// no per-state bounds check. Stops before a round that meets a slot at or
/// past `fast_limit` or lacks 8 bytes; returns the rounds done and the
/// bytes taken.
#[allow(clippy::indexing_slicing, clippy::cast_possible_truncation, reason = "x & 0xFFF < SLOTS")]
fn order_0_rounds(
    t: &Order0Tables,
    bytes: &[u8],
    states: &mut [u32; 4],
    chunks: &mut [[u8; 4]],
) -> (usize, usize) {
    let [mut x0, mut x1, mut x2, mut x3] = *states;
    let mut rest = bytes;
    let n = chunks.len();
    let mut out = chunks.iter_mut();
    let limit = t.fast_limit;
    while let Some(w) = rest.first_chunk::<8>() {
        let (m0, m1, m2, m3) = (
            (x0 & 0xFFF) as usize,
            (x1 & 0xFFF) as usize,
            (x2 & 0xFFF) as usize,
            (x3 & 0xFFF) as usize,
        );
        if m0 >= limit || m1 >= limit || m2 >= limit || m3 >= limit {
            break;
        }
        let Some(chunk) = out.next() else { break };
        *chunk = [t.sym[m0], t.sym[m1], t.sym[m2], t.sym[m3]];
        let mut k = 0usize;
        x0 = renorm_cmov(t.step_at(x0, m0), w, &mut k);
        x1 = renorm_cmov(t.step_at(x1, m1), w, &mut k);
        x2 = renorm_cmov(t.step_at(x2, m2), w, &mut k);
        x3 = renorm_cmov(t.step_at(x3, m3), w, &mut k);
        rest = rest.get(k..).unwrap_or_default();
    }
    *states = [x0, x1, x2, x3];
    (n.wrapping_sub(out.len()), bytes.len().wrapping_sub(rest.len()))
}

/// Builds the rows of every context the block can reach — 0, the symbols
/// its active contexts decode, and 255, which a slot past a row's total
/// and every inactive row decode to — inactive ones from all-zero
/// frequencies, as the reference decoder does. Other rows keep stale
/// contents that the stream never reads.
fn build_order_1_rows(active: &[bool; ALPHABET_SIZE], buf: &mut Rans4x8Buf) {
    let mut reachable = [false; ALPHABET_SIZE];
    let [first, .., last] = &mut reachable;
    (*first, *last) = (true, true);
    for (row, _) in buf.freq.iter().zip(active).filter(|(_, a)| **a) {
        for (r, &f) in reachable.iter_mut().zip(row) {
            *r |= f != 0;
        }
    }
    let rows = buf.freq.iter().zip(buf.sym_tables.iter_mut()).zip(buf.steps.iter_mut());
    for (((freq, syms), steps), (&a, &r)) in rows.zip(active.iter().zip(&reachable)) {
        if a {
            build_row(freq, syms, steps);
        } else if r {
            syms.fill(255);
            steps.fill(step_entry(0, 0));
        }
    }
}

// r[impl cram.codec.rans4x8_fast]
#[allow(
    clippy::indexing_slicing,
    clippy::cast_possible_truncation,
    clippy::arithmetic_side_effects,
    reason = "contexts are u8 < 256 rows, x & 0xFFF < 4096; i < chunk, the length of each quarter"
)]
fn decode_order_1_fast(
    src: &mut &[u8],
    dst: &mut [u8],
    buf: &mut Rans4x8Buf,
) -> Result<(), CramError> {
    let truncated = || CramError::Truncated { context: "rans order-1 truncated" };
    let active = read_frequencies_1_into(src, &mut buf.freq)?;
    build_order_1_rows(&active, buf);
    let (syms, steps) = (&*buf.sym_tables, &*buf.steps);
    let mut xs = read_states(src).ok_or_else(truncated)?;
    let bytes = *src;
    let mut pos = 0usize;
    let chunk = dst.len() / 4;
    let (q0, rest) = dst.split_at_mut(chunk);
    let (q1, rest) = rest.split_at_mut(chunk);
    let (q2, rest) = rest.split_at_mut(chunk);
    let (q3, leftover) = rest.split_at_mut(chunk);
    let mut quarters = [q0, q1, q2, q3];
    let mut ctx = [0u8; 4];
    let mut i = 0usize;
    // One checked step of state `j`.
    let checked = |x: &mut u32, c: &mut u8, rest: &mut &[u8], renorm: bool| {
        let row = usize::from(*c);
        *c = syms[row][(*x & 0xFFF) as usize];
        *x = step(*x, steps[row][usize::from(*c)]);
        if renorm {
            renormalize(x, rest).ok_or_else(truncated)?;
        }
        Ok::<u8, CramError>(*c)
    };

    if xs.iter().all(|&x| x >= LOWER_BOUND) {
        while i < chunk
            && let Some(w) = next8(bytes, pos)
        {
            let s: [u8; 4] =
                std::array::from_fn(|j| syms[usize::from(ctx[j])][(xs[j] & 0xFFF) as usize]);
            let mut next: [u32; 4] =
                std::array::from_fn(|j| step(xs[j], steps[usize::from(ctx[j])][usize::from(s[j])]));
            if next.iter().any(|&x| x < 1 << 11) {
                let mut rest = bytes.get(pos..).ok_or_else(truncated)?;
                for ((q, x), c) in quarters.iter_mut().zip(&mut xs).zip(&mut ctx) {
                    q[i] = checked(x, c, &mut rest, true)?;
                }
                pos = bytes.len().wrapping_sub(rest.len());
                i += 1;
                continue;
            }
            pos = pos.wrapping_add(renorm4(&mut next, w));
            for (q, &sym) in quarters.iter_mut().zip(&s) {
                q[i] = sym;
            }
            ctx = s;
            xs = next;
            i += 1;
        }
    }
    let mut rest = bytes.get(pos..).ok_or_else(truncated)?;
    for i in i..chunk {
        for ((q, x), c) in quarters.iter_mut().zip(&mut xs).zip(&mut ctx) {
            q[i] = checked(x, c, &mut rest, true)?;
        }
    }
    // The last state decodes the leftover bytes and skips the renorm
    // after the very last one, as the reference decoder does.
    let last = leftover.len().wrapping_sub(1);
    for (k, d) in leftover.iter_mut().enumerate() {
        *d = checked(&mut xs[3], &mut ctx[3], &mut rest, k != last)?;
    }
    *src = rest;
    Ok(())
}

/// Returns `None` when the source is exhausted mid-renormalize. The
/// caller materializes a `CramError::Truncated` from the `None` only
/// on the err path — see the helpers above for the size + drop story.
///
/// Mirrors htslib's `RansDecRenorm`: the common case for a well-formed
/// 4x8 stream is one renormalize byte per call (state was just above
/// `LOWER_BOUND >> 8` before the decode step). Two unrolled `if`-and-
/// read steps cover that and the slightly-rarer 2-byte case as
/// straight-line code, without a `while`-loop overhead per state
/// update. A small tail loop catches the pathological "state was
/// near 0" case (which only arises on malformed streams in practice
/// but is cheap to keep correct).
#[inline]
fn renormalize(state: &mut u32, src: &mut &[u8]) -> Option<()> {
    if *state >= LOWER_BOUND {
        return Some(());
    }
    *state = state.wrapping_shl(8) | u32::from(read_u8(src)?);
    if *state >= LOWER_BOUND {
        return Some(());
    }
    *state = state.wrapping_shl(8) | u32::from(read_u8(src)?);
    while *state < LOWER_BOUND {
        *state = state.wrapping_shl(8) | u32::from(read_u8(src)?);
    }
    Some(())
}

#[inline]
fn read_states(src: &mut &[u8]) -> Option<[u32; 4]> {
    let mut states = [0u32; 4];
    for s in &mut states {
        *s = read_u32_le(src)?;
    }
    Some(states)
}

// r[impl cram.codec.alphabet_run_bounded]
#[allow(clippy::indexing_slicing, reason = "sym is u8, so sym as usize ≤ 255 < ALPHABET_SIZE=256")]
fn read_frequencies_0(src: &mut &[u8]) -> Result<[u16; ALPHABET_SIZE], CramError> {
    let truncated = || CramError::Truncated { context: "rans 4x8 frequencies" };
    let mut freq = [0u16; ALPHABET_SIZE];
    let mut sym = read_u8(src).ok_or_else(truncated)?;
    let mut prev_sym = sym;

    loop {
        freq[sym as usize] = read_itf8_u16(src).ok_or_else(truncated)?;
        sym = read_u8(src).ok_or_else(truncated)?;

        if sym == 0 {
            break;
        }

        if sym == prev_sym.wrapping_add(1) {
            let run_len = read_u8(src).ok_or_else(truncated)?;
            // After the inner loop sym becomes start+run_len. htscodecs
            // rejects any run where the post-increment value wraps past
            // 255 (see `rANS_static.c::decode_freq` `if (j > 255) goto
            // cleanup`). Tolerating wraparound silently desyncs the source
            // stream — the next byte after the run would be misread as a
            // freq value.
            let end = u32::from(sym).wrapping_add(u32::from(run_len));
            if end > 255 {
                return Err(CramError::MalformedAlphabetRun { start: sym, len: run_len });
            }
            for _ in 0..run_len {
                freq[sym as usize] = read_itf8_u16(src).ok_or_else(truncated)?;
                sym = sym.wrapping_add(1);
            }
        }

        prev_sym = sym;
    }

    Ok(freq)
}

#[allow(clippy::indexing_slicing, reason = "ctx is u8, so ctx as usize ≤ 255 < ALPHABET_SIZE=256")]
/// Reads the order-1 tables into the rows of the contexts it names and
/// returns which those are; other rows keep what they had.
fn read_frequencies_1_into(
    src: &mut &[u8],
    freq: &mut [[u16; ALPHABET_SIZE]; ALPHABET_SIZE],
) -> Result<[bool; ALPHABET_SIZE], CramError> {
    let truncated = || CramError::Truncated { context: "rans 4x8 order-1 frequencies" };
    let mut active = [false; ALPHABET_SIZE];
    let mut ctx = read_u8(src).ok_or_else(truncated)?;
    let mut prev_ctx = ctx;

    loop {
        freq[ctx as usize] = read_frequencies_0(src)?;
        active[ctx as usize] = true;

        ctx = read_u8(src).ok_or_else(truncated)?;
        if ctx == 0 {
            break;
        }

        if ctx == prev_ctx.wrapping_add(1) {
            let run_len = read_u8(src).ok_or_else(truncated)?;
            // Same wraparound concern as read_frequencies_0 but for the
            // outer per-context dimension.
            let end = u32::from(ctx).wrapping_add(u32::from(run_len));
            if end > 255 {
                return Err(CramError::MalformedAlphabetRun { start: ctx, len: run_len });
            }
            for _ in 0..run_len {
                freq[ctx as usize] = read_frequencies_0(src)?;
                active[ctx as usize] = true;
                ctx = ctx.wrapping_add(1);
            }
        }

        prev_ctx = ctx;
    }

    Ok(active)
}

#[allow(
    clippy::indexing_slicing,
    reason = "i < ALPHABET_SIZE=256, sym ≤ 254 (guarded by < 255 check so sym+1 ≤ 255)"
)]
fn build_cumulative_frequencies(freq: &[u16; ALPHABET_SIZE]) -> [u16; ALPHABET_SIZE] {
    let mut cum = [0u16; ALPHABET_SIZE];
    let mut total = 0u16;
    for i in 0..ALPHABET_SIZE {
        cum[i] = total;
        total = total.wrapping_add(freq[i]);
    }
    cum
}

#[allow(
    clippy::indexing_slicing,
    reason = "sym < 255 (loop guard), so sym+1 ≤ 255 < ALPHABET_SIZE=256"
)]
fn build_symbol_table(cum_freq: &[u16; ALPHABET_SIZE]) -> [u8; 4096] {
    let mut table = [0u8; 4096];
    let mut sym = 0u8;

    for (f, g) in (0u16..).zip(&mut table) {
        while sym < 255
            && f >= cum_freq
                [(sym as usize).checked_add(1).expect("sym < 255 guarantees sym+1 ≤ 255")]
        {
            sym = sym.checked_add(1).expect("sym < 255 guarantees sym+1 ≤ 255");
        }
        *g = sym;
    }

    table
}

// Per-byte primitives (`read_u8`, `read_u32_le`) live in `super::codec_io`
// — see that module's docs for the `Option<T>` design rationale (avoiding
// `drop_in_place<CramError>` on the per-byte hot path).
use super::codec_io::{read_u8, read_u32_le};

#[inline]
fn read_itf8_u16(src: &mut &[u8]) -> Option<u16> {
    // ITF8 for frequency table values (they fit in u16)
    let val = super::varint::read_itf8_from(src)?;
    #[expect(
        clippy::cast_possible_truncation,
        reason = "rANS frequency table values are bounded by 4096 (12-bit), fits in u16"
    )]
    Some(val as u16)
}

#[cfg(test)]
#[allow(
    clippy::arithmetic_side_effects,
    clippy::needless_range_loop,
    reason = "test code: bounded values, range-by-index is clearer than enumerate here"
)]
mod tests {
    use super::*;
    use hegel::prelude::*;

    #[test]
    fn invalid_rans_order_returns_error() {
        // Build a minimal rANS 4x8 header with order=2 (invalid; only 0 and 1 are valid)
        let mut src = Vec::new();
        src.push(2u8); // order = 2 (invalid)
        src.extend_from_slice(&0u32.to_le_bytes()); // compressed_size
        src.extend_from_slice(&0u32.to_le_bytes()); // uncompressed_size

        let err = decode(&src).unwrap_err();
        assert!(matches!(err, CramError::InvalidRansOrder { order: 2 }));
    }

    #[test]
    fn invalid_rans_order_field_value() {
        let mut src = Vec::new();
        src.push(7u8); // order = 7 (invalid)
        src.extend_from_slice(&0u32.to_le_bytes());
        src.extend_from_slice(&0u32.to_le_bytes());

        let err = decode(&src).unwrap_err();
        assert!(matches!(err, CramError::InvalidRansOrder { order: 7 }));
    }

    // r[verify cram.codec.rans4x8_compressed_size_check]
    #[test]
    fn rans_compressed_size_mismatch_returns_error() {
        // Header advertises compressed_size = 100, but the actual remaining
        // payload after the 9-byte header is only ~16 bytes. htslib rejects
        // exactly this with `if (in_sz != in_size-9) return NULL`
        // (rANS_static.c:245). Catching the lie here instead of getting a
        // confusing downstream Truncated error.
        let mut src = Vec::new();
        src.push(0u8); // order = 0
        src.extend_from_slice(&100u32.to_le_bytes()); // compressed_size = 100 (lie)
        src.extend_from_slice(&7u32.to_le_bytes()); // uncompressed_size = 7
        // Only ~16 bytes of payload, far short of 100.
        src.extend_from_slice(&[0u8; 16]);

        let err = decode(&src).unwrap_err();
        assert!(
            matches!(err, CramError::Rans4x8CompressedSizeMismatch { advertised: 100, actual: 16 }),
            "expected Rans4x8CompressedSizeMismatch, got: {err:?}",
        );
    }

    // r[verify cram.codec.alphabet_run_bounded]
    #[test]
    fn read_frequencies_0_run_overflow_returns_error() {
        // sym=200, freq=ITF8(1)=0x01, then sym=201 (== prev+1) triggers a
        // run of length 100 — this would extend past sym 255. With ITF8(1)
        // payload bytes in between, the well-formed-prefix portion of the
        // stream is: [200, 0x01, 201, 100, ...freq values...]. The fix
        // must catch the bad `start + len` BEFORE consuming the run's
        // freq values.
        let mut src: &[u8] = &[200, 0x01, 201, 100];
        let err = read_frequencies_0(&mut src).unwrap_err();
        assert!(
            matches!(err, CramError::MalformedAlphabetRun { start: 201, len: 100 }),
            "expected MalformedAlphabetRun, got: {err:?}"
        );
    }

    #[test]
    fn read_frequencies_0_run_exactly_to_255_ok() {
        // Canonical htscodecs encoding for F[200..=255]=1: emit 200 (no
        // run), F[200]=1, then 201 with run-len 54 (covers F[201..254]),
        // then F[255] in the next outer iteration.
        // Stream: [200, 1, 201, 54, 1*54, 1, 0] = 60 bytes.
        let mut src: Vec<u8> = vec![200, 0x01, 201, 54];
        src.extend(std::iter::repeat_n(0x01u8, 54)); // F[201..254]
        src.push(0x01); // F[255]
        src.push(0); // terminator
        let mut cur: &[u8] = &src;
        let freq = read_frequencies_0(&mut cur).unwrap();
        for s in 200..=255usize {
            assert_eq!(freq[s], 1, "freq[{s}] should be 1");
        }
        assert_eq!(freq[199], 0, "freq[199] should be 0");
        assert_eq!(freq[0], 0, "freq[0] should be 0 (no wraparound)");
    }

    // Mirror of rans_nx16's bounds property. Only tests the rejection
    // path — constructing a fully valid 4x8 freq stream needs an extra
    // freq byte after the run for the trailing sym, which is fiddly to
    // get right at the boundary; the valid path is covered by the fixed
    // test above and integration round-trips.
    #[hegel::test]
    fn read_frequencies_0_invalid_run_rejected(tc: TestCase) {
        let start = tc.draw(gs::integers::<u8>().min_value(1));
        let len = tc.draw(gs::integers::<u8>());
        tc.assume(u32::from(start) + u32::from(len) > 255);
        let prev = start.checked_sub(1).expect("start >= 1");
        // Provide enough freq bytes that truncation can't fire before the
        // bounds check: prev + start are 2 reads; the run-len byte is 1
        // read; we don't need run payload because the check fires first.
        let stream = [prev, 0x01, start, len];
        let mut cur: &[u8] = &stream;
        let err = read_frequencies_0(&mut cur).unwrap_err();
        assert!(
            matches!(err, CramError::MalformedAlphabetRun { start: s, len: l }
                             if s == start && l == len),
            "expected MalformedAlphabetRun for start={start}, len={len}, got: {err:?}",
        );
    }

    // r[verify cram.codec.alphabet_run_bounded]
    #[test]
    fn read_frequencies_1_ctx_run_overflow_returns_error() {
        // Outer ctx loop has the same wraparound concern. Stream: ctx=200,
        // (entire sub-table for ctx=200), ctx=201 (== prev+1), run_len=100.
        // Bound check 201+100=301 > 255 must fire BEFORE the inner sub-table
        // reads. We use a minimal single-symbol sub-table for ctx=200:
        // [sym=0, freq=ITF8(4096), terminator=0].
        let src: Vec<u8> = vec![200, 0, 0x80, 0x20, 0, 201, 100];
        let mut freq = Box::new([[0u16; ALPHABET_SIZE]; ALPHABET_SIZE]);
        let mut cur: &[u8] = &src;
        let err = read_frequencies_1_into(&mut cur, &mut freq).unwrap_err();
        assert!(
            matches!(err, CramError::MalformedAlphabetRun { start: 201, len: 100 }),
            "expected MalformedAlphabetRun, got: {err:?}"
        );
    }

    // r[verify cram.codec.rans4x8]
    // r[verify cram.edge.rans_sym_overflow]
    #[test]
    fn decode_order_0_noodles_test_vector() {
        // Test vector from noodles rANS 4x8 tests: decodes to "noodles"
        let src = [
            0x00, 0x25, 0x00, 0x00, 0x00, 0x07, 0x00, 0x00, 0x00, 0x64, 0x82, 0x49, 0x65, 0x00,
            0x82, 0x49, 0x6c, 0x82, 0x49, 0x6e, 0x82, 0x49, 0x6f, 0x00, 0x84, 0x92, 0x73, 0x82,
            0x49, 0x00, 0xe2, 0x06, 0x83, 0x18, 0x74, 0x7b, 0x41, 0x0c, 0x2b, 0xa9, 0x41, 0x0c,
            0x25, 0x31, 0x80, 0x03,
        ];
        let result = decode(&src).unwrap();
        assert_eq!(result, b"noodles");
    }

    #[test]
    fn decode_order_1_noodles_test_vector() {
        // Test vector from noodles rANS 4x8 tests: decodes to "noodles"
        let src = [
            0x01, 0x3b, 0x00, 0x00, 0x00, 0x07, 0x00, 0x00, 0x00, 0x00, 0x64, 0x84, 0x00, 0x6e,
            0x84, 0x00, 0x6f, 0x00, 0x87, 0xff, 0x00, 0x64, 0x6c, 0x8f, 0xff, 0x00, 0x65, 0x00,
            0x73, 0x8f, 0xff, 0x00, 0x6c, 0x65, 0x8f, 0xff, 0x00, 0x6e, 0x6f, 0x8f, 0xff, 0x00,
            0x6f, 0x00, 0x64, 0x87, 0xff, 0x6f, 0x88, 0x00, 0x00, 0x00, 0x00, 0x04, 0x00, 0x02,
            0x02, 0x28, 0x00, 0x01, 0x02, 0x28, 0x00, 0x01, 0x02, 0x60, 0x00, 0x02,
        ];
        let result = decode(&src).unwrap();
        assert_eq!(result, b"noodles");
    }

    // ── Fast paths vs the reference decoders ────────────────────────────

    #[allow(clippy::cast_possible_truncation, reason = "v < 0x4000")]
    fn itf8_u16(v: u32) -> Vec<u8> {
        if v < 0x80 { vec![v as u8] } else { vec![0x80 | (v >> 8) as u8, (v & 0xFF) as u8] }
    }

    /// `syms` ascending, each followed by its payload, in the run-length
    /// form both frequency readers use.
    #[allow(clippy::cast_possible_truncation, reason = "a run is < 256 symbols")]
    fn encode_runs(syms: &[u8], payload: &dyn Fn(usize) -> Vec<u8>) -> Vec<u8> {
        let mut out = vec![syms[0]];
        out.extend(payload(0));
        let mut i = 1;
        while i < syms.len() {
            out.push(syms[i]);
            if syms[i] == syms[i - 1] + 1 {
                let mut k = i;
                while k + 1 < syms.len() && syms[k + 1] == syms[k] + 1 {
                    k += 1;
                }
                out.push((k - i) as u8);
                for m in i..=k {
                    out.extend(payload(m));
                }
                i = k + 1;
            } else {
                out.extend(payload(i));
                i += 1;
            }
        }
        out.push(0);
        out
    }

    /// `n` frequencies ≥ 1 summing to 4096, 4095 (htscodecs' 4x8 total), or
    /// something else that neither decoder packs.
    #[allow(clippy::cast_possible_truncation, reason = "n ≤ 256")]
    fn draw_freqs(tc: &TestCase, n: usize) -> Vec<u32> {
        let total = tc.draw(gs::sampled_from(&[4096u32, 4095, 4096, 4095, 3000, 4097]));
        let total = total.max(n as u32);
        let mut cuts: Vec<u32> = (0..n - 1)
            .map(|_| tc.draw(gs::integers::<u32>().min_value(1).max_value(total - 1)))
            .collect();
        cuts.extend([0, total]);
        cuts.sort_unstable();
        cuts.dedup();
        let mut f: Vec<u32> = cuts.windows(2).map(|w| w[1] - w[0]).collect();
        f.resize(n, 1);
        f
    }

    fn draw_syms(tc: &TestCase, with_zero: bool) -> Vec<u8> {
        let mut syms =
            tc.draw(gs::vecs(gs::integers::<u8>().min_value(1)).min_size(1).max_size(30));
        if with_zero {
            syms.push(0);
        }
        // Sometimes the top symbols, where the 4095 slot lands.
        if tc.draw(gs::booleans()) {
            syms.extend([254, 255]);
        }
        syms.sort_unstable();
        syms.dedup();
        syms
    }

    fn order0_table(tc: &TestCase, syms: &[u8]) -> Vec<u8> {
        let f = draw_freqs(tc, syms.len());
        encode_runs(syms, &|m| itf8_u16(f[m]))
    }

    /// States (above `LOWER_BOUND` or not) and payload, maybe cut short.
    fn states_and_payload(tc: &TestCase, len: usize) -> Vec<u8> {
        let low = tc.draw(gs::integers::<u8>().max_value(3)) == 0;
        let mut out = Vec::new();
        for _ in 0..4 {
            let x = if low {
                tc.draw(gs::integers::<u32>().max_value(1 << 24))
            } else {
                tc.draw(gs::integers::<u32>().min_value(LOWER_BOUND))
            };
            out.extend_from_slice(&x.to_le_bytes());
        }
        let mut payload = tc.draw(gs::binary().min_size(len + 16).max_size(len + 64));
        if tc.draw(gs::booleans()) {
            payload.truncate(tc.draw(gs::integers::<usize>().max_value(len + 16)));
        }
        out.extend(payload);
        out
    }

    fn assert_same(
        name: &str,
        got: (Result<(), CramError>, Vec<u8>, usize),
        want: &(Result<(), CramError>, Vec<u8>, usize),
    ) {
        assert_eq!(got.0.is_ok(), want.0.is_ok(), "{name}: {:?} vs {:?}", got.0, want.0);
        if got.0.is_ok() {
            assert_eq!(got.1, want.1, "{name}");
            assert_eq!(got.2, want.2, "{name} consumed differently");
        }
    }

    // r[verify cram.codec.rans4x8_packed+2]
    /// Order-0: the split-table decoder against the reference one.
    #[hegel::test]
    fn fast_order0_matches_reference(tc: TestCase) {
        let len = tc.draw(gs::integers::<usize>().max_value(400));
        let syms = draw_syms(&tc, tc.draw(gs::booleans()));
        let mut stream = order0_table(&tc, &syms);
        stream.extend(states_and_payload(&tc, len));
        let run = |f: fn(&mut &[u8], &mut [u8]) -> Result<(), CramError>| {
            let mut src: &[u8] = &stream;
            let mut dst = vec![0u8; len];
            let r = f(&mut src, &mut dst);
            (r, dst, src.len())
        };
        assert_same("fast", run(decode_order_0_fast), &run(decode_order_0));
    }

    // r[verify cram.codec.rans4x8_fast]
    /// Order-1: the split-table decoder against the reference one, and a
    /// reused buffer against a fresh one.
    #[hegel::test]
    fn fast_order1_matches_reference(tc: TestCase) {
        let mut reused = Rans4x8Buf::new();
        for _ in 0..2 {
            let len = tc.draw(gs::integers::<usize>().max_value(400));
            let ctxs = draw_syms(&tc, tc.draw(gs::integers::<u8>().max_value(7)) != 0);
            let rows: Vec<Vec<u8>> = ctxs
                .iter()
                .map(|_| {
                    let syms = draw_syms(&tc, false);
                    order0_table(&tc, &syms)
                })
                .collect();
            let mut stream = encode_runs(&ctxs, &|m| rows[m].clone());
            stream.extend(states_and_payload(&tc, len));
            type Decode = fn(&mut &[u8], &mut [u8], &mut Rans4x8Buf) -> Result<(), CramError>;
            let run = |f: Decode, buf: &mut Rans4x8Buf| {
                let mut src: &[u8] = &stream;
                let mut dst = vec![0u8; len];
                let r = f(&mut src, &mut dst, buf);
                (r, dst, src.len())
            };
            let want = run(|src, dst, _| decode_order_1_buf(src, dst), &mut Rans4x8Buf::new());
            assert_same("fast", run(decode_order_1_fast, &mut Rans4x8Buf::new()), &want);
            assert_same("fast, reused buffer", run(decode_order_1_fast, &mut reused), &want);
        }
    }
}
