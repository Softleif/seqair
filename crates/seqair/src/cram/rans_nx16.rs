// r[impl cram.codec.rans_nx16]
//! rANS Nx16 codec (CRAM compression method 5).
//!
//! N-way interleaved asymmetric numeral systems with 16-bit renormalization.
//! Supports order-0 and order-1 modes, plus STRIPE, RLE, and PACK transforms.

// See rans.rs for the rationale: lazy `ok_or_else(|| CramError::...)`
// is the deliberate hot-path choice because eager construction triggers
// a per-call `drop_in_place<CramError>` that shows up in profiles.
#![allow(
    clippy::unnecessary_lazy_evaluations,
    reason = "lazy form avoids per-call drop_in_place<CramError> on hot path"
)]

use super::codec_io::{self, Uint7Error};
use super::reader::CramError;
use fearless_simd::{Level, dispatch, prelude::*, u8x16, u32x4};
use fearless_simd_macros::simd;

/// Bridge `Uint7Error` (narrow, hot-path-friendly) to the rich `CramError`.
fn uint7_to_cram_error(e: Uint7Error) -> CramError {
    match e {
        Uint7Error::Truncated => CramError::Truncated { context: "rans_nx16 uint7" },
        Uint7Error::Overflow => CramError::Uint7Overflow,
    }
}

const ALPHABET_SIZE: usize = 256;

const FLAG_ORDER: u8 = 0x01;
const FLAG_N32: u8 = 0x04;
const FLAG_STRIPE: u8 = 0x08;
const FLAG_NO_SIZE: u8 = 0x10;
const FLAG_CAT: u8 = 0x20;
const FLAG_RLE: u8 = 0x40;
const FLAG_PACK: u8 = 0x80;

type Frequencies1 = Box<[[u32; ALPHABET_SIZE]; ALPHABET_SIZE]>;
type CumulativeFrequencies1 = Box<[[u32; ALPHABET_SIZE]; ALPHABET_SIZE]>;
/// 256 contexts x 4096-entry symbol-decode table. ~1 MiB.
type SymTables1 = Box<[[u8; 4096]; ALPHABET_SIZE]>;

/// Reusable allocations for rANS Nx16 order-1 decoding.
///
/// Allocated on the first order-1 block, so a buffer that only ever sees
/// order-0 streams costs nothing. The packed rows are 4 MiB of zeroed
/// (lazily mapped) memory of which a block touches only its active
/// contexts' rows; the unpacked tables serve malformed frequency tables only.
#[derive(Default)]
pub(crate) struct Nx16Order1Buf {
    tables: Option<Order1Tables>,
}

struct Order1Tables {
    frequencies: Frequencies1,
    packed: Box<[PackedRow; ALPHABET_SIZE]>,
    /// The unpacked path's tables, for frequency tables that do not pack.
    cumulative_frequencies: CumulativeFrequencies1,
    sym_tables: SymTables1,
    states: Vec<u32>,
    prev_syms: Vec<u8>,
}

/// A zeroed `Box<[T; N]>` straight from the allocator (`vec!` of a zero
/// value is a `calloc`), never staged on the stack.
fn zeroed_box<T: Clone, const N: usize>(zero: T) -> Box<[T; N]> {
    let Ok(b) = vec![zero; N].into_boxed_slice().try_into() else {
        unreachable!("a Vec of N elements converts to [T; N]")
    };
    b
}

impl Order1Tables {
    fn new() -> Self {
        Self {
            frequencies: zeroed_box([0; ALPHABET_SIZE]),
            packed: zeroed_box([0; ROW]),
            cumulative_frequencies: zeroed_box([0; ALPHABET_SIZE]),
            sym_tables: zeroed_box([0; 4096]),
            states: Vec::with_capacity(32),
            prev_syms: Vec::with_capacity(32),
        }
    }
}

impl Nx16Order1Buf {
    pub fn new() -> Self {
        Self::default()
    }
}

/// Decode a rANS Nx16 compressed block (allocating path — for callers
/// without a reusable buffer). Prefer `decode_with_buf` on the hot path.
pub fn decode(src: &[u8], uncompressed_size: usize) -> Result<Vec<u8>, CramError> {
    let mut buf = Nx16Order1Buf::new();
    decode_with_buf(src, uncompressed_size, &mut buf)
}

/// Decode a rANS Nx16 compressed block, reusing the caller-owned
/// [`Nx16Order1Buf`] for order-1 tables. Order-0 and transform paths
/// don't touch `buf` (the stripe path re-enters via `decode_with_buf`
/// for each sub-block).
pub(crate) fn decode_with_buf(
    src: &[u8],
    mut uncompressed_size: usize,
    buf: &mut Nx16Order1Buf,
) -> Result<Vec<u8>, CramError> {
    let mut cur: &[u8] = src;

    let flags =
        read_u8(&mut cur).ok_or_else(|| CramError::Truncated { context: "rans_nx16 flags" })?;
    let state_count = if flags & FLAG_N32 != 0 { 32 } else { 4 };

    if flags & FLAG_NO_SIZE == 0 {
        uncompressed_size = read_uint7(&mut cur)? as usize;
    }
    super::reader::check_alloc_size(uncompressed_size, "rANS Nx16 output")?;

    if flags & FLAG_STRIPE != 0 {
        return decode_stripe_with_buf(&mut cur, uncompressed_size, buf);
    }

    let bit_pack_ctx = if flags & FLAG_PACK != 0 {
        let (ctx, len) = read_bit_pack_context(&mut cur, uncompressed_size)?;
        uncompressed_size = len;
        Some(ctx)
    } else {
        None
    };

    let rle_ctx = if flags & FLAG_RLE != 0 {
        let (ctx, len) = read_rle_context(&mut cur, state_count, uncompressed_size)?;
        uncompressed_size = len;
        Some(ctx)
    } else {
        None
    };

    // r[impl io.fuzz.alloc_limits]
    // Re-check after PACK/RLE transforms may have updated uncompressed_size.
    super::reader::check_alloc_size(uncompressed_size, "rANS Nx16 output (post-transform)")?;
    let mut dst = vec![0u8; uncompressed_size];

    if flags & FLAG_CAT != 0 {
        let data = split_off(&mut cur, uncompressed_size)?;
        dst.copy_from_slice(data);
    } else if flags & FLAG_ORDER == 0 {
        decode_order_0(&mut cur, &mut dst, state_count)?;
    } else {
        decode_order_1_with_buf(&mut cur, &mut dst, state_count, buf)?;
    }

    if let Some(ctx) = rle_ctx {
        dst = apply_rle(&dst, &ctx)?;
    }

    if let Some(ctx) = bit_pack_ctx {
        dst = apply_bit_unpack(&dst, &ctx)?;
    }

    Ok(dst)
}

// ── Primitive readers ────────────────────────────────────────────────

// All per-byte primitives (`read_u8`, `read_u32_le`, `read_uint7`,
// `split_off`) live in `super::codec_io`; see that module's docs for
// the `Option<T>` design rationale (avoiding `drop_in_place<CramError>`
// on the per-byte hot path).
pub(crate) use super::codec_io::read_u16_le as read_u16_le_prv;
use super::codec_io::{read_u8, read_u32_le};

/// Local thin wrappers that produce a `CramError` with a tagged context
/// from the narrow `Option`/`Uint7Error` returns.
#[inline]
fn read_uint7(src: &mut &[u8]) -> Result<u32, CramError> {
    codec_io::read_uint7(src).map_err(uint7_to_cram_error)
}

#[inline]
fn split_off<'a>(src: &mut &'a [u8], len: usize) -> Result<&'a [u8], CramError> {
    codec_io::split_off(src, len)
        .ok_or_else(|| CramError::Truncated { context: "rans_nx16 split_off" })
}

fn read_states(src: &mut &[u8], state_count: usize) -> Result<Vec<u32>, CramError> {
    (0..state_count)
        .map(|_| {
            read_u32_le(src).ok_or_else(|| CramError::Truncated { context: "rans_nx16 state" })
        })
        .collect()
}

// ── Core rANS step functions ─────────────────────────────────────────

fn state_cumulative_frequency(s: u32, bits: u32) -> u32 {
    let mask = 1u32
        .checked_shl(bits)
        .and_then(|v| v.checked_sub(1))
        .expect("bits < 32 for valid rANS state");
    s & mask
}

/// Linear-scan oracle for the precomputed `build_symbol_table_nx16` table.
/// Kept for tests as an independent check; the decode hot loops now use
/// the precomputed tables directly (order-0: `build_symbol_table_nx16`;
/// order-1: `Nx16Order1Buf::sym_tables` populated by
/// `build_symbol_table_nx16_into`).
#[cfg(test)]
fn cumulative_frequencies_symbol(
    cumulative_frequencies: &[u32; ALPHABET_SIZE],
    frequency: u32,
) -> u8 {
    let mut sym: u8 = 0;
    while sym < 255
        && frequency >= *cumulative_frequencies.get(usize::from(sym) + 1).unwrap_or(&u32::MAX)
    {
        sym = sym.wrapping_add(1);
    }
    sym
}

// r[impl cram.codec.state_step_safety]
// The CRAM spec guarantees `g <= f * (s >> bits) + (s & mask)` for valid frequency
// tables, but malformed/fuzz input can violate it — `wrapping_sub` keeps the
// release path safe. No `debug_assert!`: the input is untrusted (CRAM data),
// so a debug assertion would panic on adversarial inputs in fuzz/debug builds.
pub(crate) fn state_step(s: u32, f: u32, g: u32, bits: u32) -> u32 {
    let result =
        f.wrapping_mul(s >> bits).wrapping_add(s & (1u32.wrapping_shl(bits).wrapping_sub(1)));
    // r[depends cram.codec.state_step_safety]
    result.wrapping_sub(g)
}

/// Returns `None` when the source is exhausted mid-renormalize. The
/// caller materializes a `CramError::Truncated` from the `None` only on
/// the err path — see the helper-module comment for the size + drop
/// story.
///
/// Two unrolled read steps cover the common 1- and 2-byte cases as
/// straight-line code (matching htslib's 16-bit renorm pattern). The
/// tail `while` only fires on malformed data where the state was
/// near zero. `wrapping_shl` + `|` replaces `checked_*` because
/// `s < 0x8000` guarantees `s << 16 ≤ 0x7FFF0000` (always fits in u32).
#[inline]
pub(crate) fn state_renormalize(mut s: u32, src: &mut &[u8]) -> Option<u32> {
    if s < (1 << 15) {
        s = s.wrapping_shl(16) | u32::from(read_u16_le_prv(src)?);
        if s < (1 << 15) {
            s = s.wrapping_shl(16) | u32::from(read_u16_le_prv(src)?);
            while s < (1 << 15) {
                s = s.wrapping_shl(16) | u32::from(read_u16_le_prv(src)?);
            }
        }
    }
    Some(s)
}

// ── Packed fast path ─────────────────────────────────────────────────
//
// htscodecs' layout: one `u32` per table slot holding everything a decode
// step needs, so a step is one load instead of three dependent ones
// (`sym_table[x]`, then `freq[sym]` and `cum[sym]`). A slot is
// `(f − 1) << 20 | bias << 8 | sym`, where `bias` is the slot's offset into
// its symbol's range, so the step is `f·(x >> bits) + bias`. Only tables
// whose frequencies sum to exactly `1 << bits` are packed: then every slot
// belongs to a symbol with f ≥ 1, `f − 1` and `bias` fit 12 bits, and the
// packed step equals `state_step` bit for bit. Any other table takes the
// unpacked scalar path, which keeps whatever it does on malformed tables.

const PACK_F_SHIFT: u32 = 20;
const PACK_BIAS_MASK: u32 = 0xFFF;
/// Slots per packed row: the table width at the largest allowed `bits` (12).
pub(crate) const ROW: usize = 1 << ORDER_0_BITS;
pub(crate) type PackedRow = [u32; ROW];

/// Pack one context's frequencies into `row[..1 << bits]`. False unless the
/// frequencies sum to exactly `1 << bits`; `row` is then partly written and
/// must not be used.
#[allow(clippy::indexing_slicing, reason = "x + f ≤ total ≤ ROW is checked before the slice")]
#[allow(clippy::arithmetic_side_effects, reason = "f ≥ 1 in the subtraction")]
pub(crate) fn pack_row(
    freqs: impl IntoIterator<Item = u32>,
    bits: u32,
    row: &mut PackedRow,
) -> bool {
    pack_row_to(freqs, 1 << bits, row)
}

/// [`pack_row`] for a table that must sum to exactly `total` (≤ `ROW`).
#[allow(clippy::indexing_slicing, reason = "x + f ≤ total ≤ ROW is checked before the slice")]
#[allow(clippy::arithmetic_side_effects, reason = "f ≥ 1 in the subtraction")]
pub(crate) fn pack_row_to(
    freqs: impl IntoIterator<Item = u32>,
    total: u32,
    row: &mut PackedRow,
) -> bool {
    let mut x = 0u32;
    for (sym, f) in (0u32..=255).zip(freqs) {
        if f == 0 {
            continue;
        }
        let Some(end) = x.checked_add(f).filter(|&end| end <= total) else {
            return false;
        };
        let base = ((f - 1) << PACK_F_SHIFT) | sym;
        for (y, slot) in (0u32..).zip(&mut row[x as usize..end as usize]) {
            *slot = base | (y << 8);
        }
        x = end;
    }
    x == total
}

/// The scalar packed step, for the scalar tails of the packed kernels.
#[inline(always)]
pub(crate) fn packed_step(x: u32, slot: u32, bits: u32) -> u32 {
    ((slot >> PACK_F_SHIFT).wrapping_add(1))
        .wrapping_mul(x >> bits)
        .wrapping_add((slot >> 8) & PACK_BIAS_MASK)
}

/// For each 4-bit "lane needs a renorm" mask, the byte shuffle that moves
/// the next u16s of the stream into the low halves of exactly those lanes,
/// in lane order. The other bytes are don't-care (0x80): what a swizzle
/// makes of an out-of-range index differs by level (zero on NEON and x86,
/// a wrapped index on the fallback), so the caller masks the low halves
/// and keeps the old state in lanes that do not renormalize.
#[allow(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "const table build, m < 16, lane < 4, k < 4"
)]
const RENORM_SHUFFLE: [[u8; 16]; 16] = {
    let mut t = [[0x80u8; 16]; 16];
    let mut m = 0;
    while m < 16 {
        let mut k = 0u8;
        let mut lane = 0;
        while lane < 4 {
            if m & (1 << lane) != 0 {
                t[m][lane * 4] = 2 * k;
                t[m][lane * 4 + 1] = 2 * k + 1;
                k += 1;
            }
            lane += 1;
        }
        m += 1;
    }
    t
};

/// Renormalize four states at once: every lane below `1 << 15` takes the
/// next u16 of the stream, in lane order, like four scalar
/// `state_renormalize` calls that each read at most once. Returns the new
/// states and the number of bytes consumed. `window` is the next 16 bytes.
#[inline(always)]
#[allow(clippy::indexing_slicing, reason = "the bitmask of 4 lanes is < 16")]
#[allow(clippy::cast_possible_truncation, reason = "the bitmask of 4 lanes is < 16")]
#[allow(clippy::arithmetic_side_effects, reason = "count_ones ≤ 4")]
fn renorm4<S: Simd>(simd: S, x: u32x4<S>, window: &[u8; 16]) -> (u32x4<S>, usize) {
    let low = x.simd_lt(u32x4::splat(simd, 1 << 15));
    let m = (low.to_bitmask() & 0xF) as usize;
    let raw = u8x16::from_slice(simd, window);
    let norm: u32x4<S> = raw.swizzle_dyn(u8x16::from_slice(simd, &RENORM_SHUFFLE[m])).bitcast();
    (low.select((x << 16) | (norm & 0xFFFF), x), 2 * m.count_ones() as usize)
}

/// Bytes the vector loop must have left before a step: at most 2 per lane,
/// plus the 16-byte window of the last group, which starts at most
/// `2·(n − 4)` bytes in.
#[allow(clippy::arithmetic_side_effects, reason = "n ≤ 32")]
const fn vector_headroom(n: usize) -> usize {
    2 * n + 16
}

/// Order-0 decode with `N` interleaved states (4 or 32) over a packed
/// table. The vector loop needs every state ≥ `1 << 15`, so that one
/// renorm read per step suffices, exactly as in the scalar decoder (a
/// packed slot has f ≥ 1, so a step leaves x ≥ `x >> 12` ≥ 8 and one u16
/// lifts it past `1 << 15`); near the end of the stream, or with a state
/// below that, it hands over to the scalar steps, which renormalize with
/// bounds checks.
#[simd]
#[allow(clippy::indexing_slicing, reason = "x & 0xFFF < ROW indexes the table")]
#[allow(clippy::cast_possible_truncation, reason = "the low byte of a slot is its symbol")]
#[allow(
    clippy::arithmetic_side_effects,
    reason = "vector lanes wrap, as state_step's do; pos advances by ≤ 2 per lane within the checked headroom"
)]
fn decode_o0_packed<S: Simd, const N: usize>(
    simd: S,
    src: &mut &[u8],
    dst: &mut [u8],
    table: &PackedRow,
    states: &mut [u32; N],
) -> Result<(), CramError> {
    debug_assert_eq!(N % 4, 0);
    let truncated = || CramError::Truncated { context: "rans_nx16 order-0 truncated" };
    let bytes: &[u8] = src;
    let mut pos = 0usize;
    let (chunks, remainder) = dst.as_chunks_mut::<N>();
    let mut chunks = chunks.iter_mut();

    if states.iter().all(|&x| x >= 1 << 15) {
        while bytes.len().saturating_sub(pos) >= vector_headroom(N) {
            let Some(chunk) = chunks.next() else { break };
            let mut slots = [0u32; N];
            for ((d, s), &x) in chunk.iter_mut().zip(&mut slots).zip(states.iter()) {
                *s = table[(x & 0xFFF) as usize];
                *d = *s as u8;
            }
            for (x4, s4) in states.as_chunks_mut::<4>().0.iter_mut().zip(slots.as_chunks::<4>().0) {
                let x = u32x4::from_slice(simd, x4);
                let s = u32x4::from_slice(simd, s4);
                let x =
                    ((s >> PACK_F_SHIFT) + 1) * (x >> ORDER_0_BITS) + ((s >> 8) & PACK_BIAS_MASK);
                let Some(window) = bytes.get(pos..).and_then(|w| w.first_chunk::<16>()) else {
                    return Err(truncated());
                };
                let (x, used) = renorm4(simd, x, window);
                pos += used;
                x.store_slice(x4);
            }
        }
    }

    let mut rest = bytes.get(pos..).ok_or_else(truncated)?;
    for chunk in chunks {
        for (d, x) in chunk.iter_mut().zip(states.iter_mut()) {
            let slot = table[(*x & 0xFFF) as usize];
            *d = slot as u8;
            *x = state_renormalize(packed_step(*x, slot, ORDER_0_BITS), &mut rest)
                .ok_or_else(truncated)?;
        }
    }
    for (d, x) in remainder.iter_mut().zip(states.iter_mut()) {
        let slot = table[(*x & 0xFFF) as usize];
        *d = slot as u8;
        *x = state_renormalize(packed_step(*x, slot, ORDER_0_BITS), &mut rest)
            .ok_or_else(truncated)?;
    }
    *src = rest;
    Ok(())
}

/// Order-1 decode with `N` interleaved states over packed per-context rows:
/// state `j` fills `dst[j·len/N ..][..len/N]`, its context is the symbol it
/// decoded last, and the last state also decodes the `len % N` leftover
/// bytes. Same vector/scalar split as [`decode_o0_packed`]. Only rows of
/// contexts reachable from context 0 are read, and the caller has packed
/// all of them.
#[simd]
#[allow(
    clippy::indexing_slicing,
    reason = "ctx: u8 < 256 rows, x & mask ≤ 0xFFF < ROW; j·chunk + i < N·chunk ≤ dst.len()"
)]
#[allow(clippy::cast_possible_truncation, reason = "the low byte of a slot is its symbol")]
#[allow(
    clippy::arithmetic_side_effects,
    reason = "vector lanes wrap, as state_step's do; j·chunk + i < dst.len(); pos within the checked headroom"
)]
fn decode_o1_packed<S: Simd, const N: usize>(
    simd: S,
    src: &mut &[u8],
    dst: &mut [u8],
    table: &[PackedRow; ALPHABET_SIZE],
    bits: u32,
    states: &mut [u32; N],
) -> Result<(), CramError> {
    debug_assert_eq!(N % 4, 0);
    let truncated = || CramError::Truncated { context: "rans_nx16 order-1 truncated" };
    let mask = ((1u32 << bits) - 1) & 0xFFF;
    let chunk = dst.len() / N;
    let bytes: &[u8] = src;
    let mut pos = 0usize;
    let mut ctx = [0u8; N];
    let mut i = 0usize;

    if states.iter().all(|&x| x >= 1 << 15) {
        while i < chunk && bytes.len().saturating_sub(pos) >= vector_headroom(N) {
            let mut slots = [0u32; N];
            for (j, ((s, c), &x)) in slots.iter_mut().zip(&mut ctx).zip(states.iter()).enumerate() {
                *s = table[usize::from(*c)][(x & mask) as usize];
                *c = *s as u8;
                dst[j * chunk + i] = *c;
            }
            for (x4, s4) in states.as_chunks_mut::<4>().0.iter_mut().zip(slots.as_chunks::<4>().0) {
                let x = u32x4::from_slice(simd, x4);
                let s = u32x4::from_slice(simd, s4);
                let x = ((s >> PACK_F_SHIFT) + 1) * (x >> bits) + ((s >> 8) & PACK_BIAS_MASK);
                let Some(window) = bytes.get(pos..).and_then(|w| w.first_chunk::<16>()) else {
                    return Err(truncated());
                };
                let (x, used) = renorm4(simd, x, window);
                pos += used;
                x.store_slice(x4);
            }
            i += 1;
        }
    }

    let mut rest = bytes.get(pos..).ok_or_else(truncated)?;
    for i in i..chunk {
        for (j, (c, x)) in ctx.iter_mut().zip(states.iter_mut()).enumerate() {
            let slot = table[usize::from(*c)][(*x & mask) as usize];
            *c = slot as u8;
            dst[j * chunk + i] = *c;
            *x = state_renormalize(packed_step(*x, slot, bits), &mut rest).ok_or_else(truncated)?;
        }
    }
    let (Some(c), Some(x)) = (ctx.last_mut(), states.last_mut()) else {
        return Err(truncated());
    };
    for d in &mut dst[chunk * N..] {
        let slot = table[usize::from(*c)][(*x & mask) as usize];
        *c = slot as u8;
        *d = *c;
        *x = state_renormalize(packed_step(*x, slot, bits), &mut rest).ok_or_else(truncated)?;
    }
    *src = rest;
    Ok(())
}

/// One branch-free renorm read: `x` takes the u16 at `pos` if it is below
/// `1 << 15`. The caller guarantees two readable bytes at `pos`.
#[inline(always)]
#[allow(clippy::arithmetic_side_effects, reason = "pos + 2 within the caller's headroom")]
fn renorm_branchless(x: u32, bytes: &[u8], pos: &mut usize) -> u32 {
    let y =
        bytes.get(*pos..).and_then(|w| w.first_chunk::<2>()).map_or(0, |w| u16::from_le_bytes(*w));
    let low = x < 1 << 15;
    *pos += 2 * usize::from(low);
    if low { (x << 16) | u32::from(y) } else { x }
}

/// [`decode_o0_packed`] without vectors: the same packed step and the same
/// one-read renorm, branch-free. For few states, where moving four lanes
/// in and out of a vector costs more than it saves.
#[allow(clippy::indexing_slicing, reason = "x & 0xFFF < ROW indexes the table")]
#[allow(clippy::cast_possible_truncation, reason = "the low byte of a slot is its symbol")]
fn decode_o0_packed_scalar<const N: usize>(
    src: &mut &[u8],
    dst: &mut [u8],
    table: &PackedRow,
    states: &mut [u32; N],
) -> Result<(), CramError> {
    let truncated = || CramError::Truncated { context: "rans_nx16 order-0 truncated" };
    let bytes: &[u8] = src;
    let mut pos = 0usize;
    let (chunks, remainder) = dst.as_chunks_mut::<N>();
    let mut chunks = chunks.iter_mut();

    if states.iter().all(|&x| x >= 1 << 15) {
        while bytes.len().saturating_sub(pos) >= N.saturating_mul(2) {
            let Some(chunk) = chunks.next() else { break };
            for (d, x) in chunk.iter_mut().zip(states.iter_mut()) {
                let slot = table[(*x & 0xFFF) as usize];
                *d = slot as u8;
                *x = renorm_branchless(packed_step(*x, slot, ORDER_0_BITS), bytes, &mut pos);
            }
        }
    }

    let mut rest = bytes.get(pos..).ok_or_else(truncated)?;
    for chunk in chunks {
        for (d, x) in chunk.iter_mut().zip(states.iter_mut()) {
            let slot = table[(*x & 0xFFF) as usize];
            *d = slot as u8;
            *x = state_renormalize(packed_step(*x, slot, ORDER_0_BITS), &mut rest)
                .ok_or_else(truncated)?;
        }
    }
    for (d, x) in remainder.iter_mut().zip(states.iter_mut()) {
        let slot = table[(*x & 0xFFF) as usize];
        *d = slot as u8;
        *x = state_renormalize(packed_step(*x, slot, ORDER_0_BITS), &mut rest)
            .ok_or_else(truncated)?;
    }
    *src = rest;
    Ok(())
}

/// [`decode_o1_packed`] without vectors; see [`decode_o0_packed_scalar`].
#[allow(
    clippy::indexing_slicing,
    reason = "ctx: u8 < 256 rows, x & mask ≤ 0xFFF < ROW; j·chunk + i < N·chunk ≤ dst.len()"
)]
#[allow(clippy::cast_possible_truncation, reason = "the low byte of a slot is its symbol")]
#[allow(clippy::arithmetic_side_effects, reason = "j·chunk + i < dst.len()")]
fn decode_o1_packed_scalar<const N: usize>(
    src: &mut &[u8],
    dst: &mut [u8],
    table: &[PackedRow; ALPHABET_SIZE],
    bits: u32,
    states: &mut [u32; N],
) -> Result<(), CramError> {
    let truncated = || CramError::Truncated { context: "rans_nx16 order-1 truncated" };
    let mask = ((1u32 << bits) - 1) & 0xFFF;
    let chunk = dst.len() / N;
    let bytes: &[u8] = src;
    let mut pos = 0usize;
    let mut ctx = [0u8; N];
    let mut i = 0usize;

    if states.iter().all(|&x| x >= 1 << 15) {
        while i < chunk && bytes.len().saturating_sub(pos) >= 2 * N {
            for (j, (c, x)) in ctx.iter_mut().zip(states.iter_mut()).enumerate() {
                let slot = table[usize::from(*c)][(*x & mask) as usize];
                *c = slot as u8;
                dst[j * chunk + i] = *c;
                *x = renorm_branchless(packed_step(*x, slot, bits), bytes, &mut pos);
            }
            i += 1;
        }
    }

    let mut rest = bytes.get(pos..).ok_or_else(truncated)?;
    for i in i..chunk {
        for (j, (c, x)) in ctx.iter_mut().zip(states.iter_mut()).enumerate() {
            let slot = table[usize::from(*c)][(*x & mask) as usize];
            *c = slot as u8;
            dst[j * chunk + i] = *c;
            *x = state_renormalize(packed_step(*x, slot, bits), &mut rest).ok_or_else(truncated)?;
        }
    }
    let (Some(c), Some(x)) = (ctx.last_mut(), states.last_mut()) else {
        return Err(truncated());
    };
    for d in &mut dst[chunk * N..] {
        let slot = table[usize::from(*c)][(*x & mask) as usize];
        *c = slot as u8;
        *d = *c;
        *x = state_renormalize(packed_step(*x, slot, bits), &mut rest).ok_or_else(truncated)?;
    }
    *src = rest;
    Ok(())
}

// ── Alphabet reading ─────────────────────────────────────────────────

// r[impl cram.codec.alphabet_run_bounded]
#[allow(
    clippy::indexing_slicing,
    reason = "sym is u8 so usize::from(sym) ≤ 255 < ALPHABET_SIZE=256"
)]
fn read_alphabet(src: &mut &[u8]) -> Result<[bool; ALPHABET_SIZE], CramError> {
    let truncated = || CramError::Truncated { context: "rans_nx16 alphabet" };
    let mut alphabet = [false; ALPHABET_SIZE];

    let mut sym = read_u8(src).ok_or_else(truncated)?;
    let mut prev_sym = sym;

    loop {
        alphabet[usize::from(sym)] = true;

        sym = read_u8(src).ok_or_else(truncated)?;

        if sym == 0 {
            break;
        }

        if sym == prev_sym.wrapping_add(1) {
            let len = read_u8(src).ok_or_else(truncated)?;
            // After the inner loop sym becomes start+len. htscodecs rejects
            // any run where the post-increment value wraps past 255 (see
            // `rANS_static16_int.h:decode_alphabet` `if (j > 255) return 0`).
            // Tolerating wraparound silently corrupts alphabet[0..] entries
            // and desyncs the source stream.
            let end = u32::from(sym).wrapping_add(u32::from(len));
            if end > 255 {
                return Err(CramError::MalformedAlphabetRun { start: sym, len });
            }
            for _ in 0..len {
                alphabet[usize::from(sym)] = true;
                sym = sym.wrapping_add(1);
            }
        }

        prev_sym = sym;
    }

    Ok(alphabet)
}

// ── Order-0 ──────────────────────────────────────────────────────────

pub(crate) const ORDER_0_BITS: u32 = 12;

fn decode_order_0(src: &mut &[u8], dst: &mut [u8], state_count: usize) -> Result<(), CramError> {
    decode_order_0_at(Level::new(), src, dst, state_count)
}

/// Frequency tables that pack take the packed kernels; the rest re-read the
/// table from the start through the unpacked paths.
// r[impl cram.codec.simd_dispatch+3]
fn decode_order_0_at(
    level: Level,
    src: &mut &[u8],
    dst: &mut [u8],
    state_count: usize,
) -> Result<(), CramError> {
    let mut cur = *src;
    let frequencies = read_frequencies_0(&mut cur)?;
    let mut table = [0u32; ROW];
    if pack_row(frequencies, ORDER_0_BITS, &mut table) {
        match state_count {
            32 => {
                let mut states = read_state_array::<32>(&mut cur)?;
                dispatch!(level, simd => decode_o0_packed(simd, &mut cur, dst, &table, &mut states))?;
                *src = cur;
                return Ok(());
            }
            // Four states: moving one vector's worth of lanes in and out
            // costs more than the vector step saves (measured).
            4 => {
                let mut states = read_state_array::<4>(&mut cur)?;
                decode_o0_packed_scalar(&mut cur, dst, &table, &mut states)?;
                *src = cur;
                return Ok(());
            }
            _ => {}
        }
    }
    if state_count == 32 {
        return decode_order_0_32state_at(level, src, dst);
    }
    decode_order_0_generic(src, dst, state_count)
}

fn read_state_array<const N: usize>(src: &mut &[u8]) -> Result<[u32; N], CramError> {
    let mut states = [0u32; N];
    for x in &mut states {
        *x = read_u32_le(src).ok_or_else(|| CramError::Truncated { context: "rans_nx16 state" })?;
    }
    Ok(states)
}

/// The 32-state order-0 decode, through the unpacked SIMD kernel.
#[cfg(test)]
fn decode_order_0_32state(src: &mut &[u8], dst: &mut [u8]) -> Result<(), CramError> {
    decode_order_0_32state_at(Level::new(), src, dst)
}

// r[impl io.simd_portable]
fn decode_order_0_32state_at(
    level: Level,
    src: &mut &[u8],
    dst: &mut [u8],
) -> Result<(), CramError> {
    let frequencies = read_frequencies_0(src)?;
    let cumulative_frequencies = build_cumulative_frequencies(&frequencies);
    let sym_table = build_symbol_table_nx16(&cumulative_frequencies);
    let Ok(mut states) = <[u32; 32]>::try_from(read_states(src, 32)?) else {
        return Err(CramError::Truncated { context: "rans_nx16 order-0 states" });
    };
    dispatch!(level, simd => decode_32state_simd(
        simd,
        src,
        dst,
        &frequencies,
        &cumulative_frequencies,
        &sym_table,
        &mut states,
    ))
}

/// Per 32-byte chunk: a scalar symbol lookup per state (a gather has no
/// portable spelling, and `sym_table` is 4 KiB), then the state step
/// `f·(x >> 12) + (x & 0xFFF) − g` on native-width vectors of states, then
/// the shared scalar `state_renormalize`, so this path cannot drift from the
/// scalar one. The tail cycles through the states like the scalar decoder
/// does: the stream interleaves bytes across all 32.
#[simd]
#[allow(
    clippy::indexing_slicing,
    reason = "sym ≤ 255 (u8) indexes the 256-entry tables, f < 4096 (12-bit mask) indexes sym_table"
)]
#[allow(
    clippy::arithmetic_side_effects,
    reason = "vector lanes wrap, as state_step's wrapping_* ops do"
)]
#[allow(
    clippy::chunks_exact_to_as_chunks,
    reason = "`S::u32s::LEN` depends on the generic `S`, so it cannot be a const argument"
)]
fn decode_32state_simd<S: Simd>(
    simd: S,
    src: &mut &[u8],
    dst: &mut [u8],
    frequencies: &[u32; ALPHABET_SIZE],
    cumulative_frequencies: &[u32; ALPHABET_SIZE],
    sym_table: &[u8; 4096],
    states: &mut [u32; 32],
) -> Result<(), CramError> {
    let lanes = S::u32s::LEN;
    debug_assert_eq!(32 % lanes, 0, "32 states split into whole vectors");
    let truncated = || CramError::Truncated { context: "rans_nx16 order-0 truncated" };
    let (chunks, remainder) = dst.as_chunks_mut::<32>();

    for chunk in chunks {
        let mut freq = [0u32; 32];
        let mut cum = [0u32; 32];
        for (((d, &x), f), g) in chunk.iter_mut().zip(states.iter()).zip(&mut freq).zip(&mut cum) {
            let sym = sym_table[(x & 0xFFF) as usize];
            *d = sym;
            *f = frequencies[usize::from(sym)];
            *g = cumulative_frequencies[usize::from(sym)];
        }

        for ((x, f), g) in states
            .chunks_exact_mut(lanes)
            .zip(freq.chunks_exact(lanes))
            .zip(cum.chunks_exact(lanes))
        {
            let v = S::u32s::from_slice(simd, x);
            let f = S::u32s::from_slice(simd, f);
            let g = S::u32s::from_slice(simd, g);
            (f * (v >> ORDER_0_BITS) + (v & 0xFFF) - g).store_slice(x);
        }

        for x in states.iter_mut() {
            *x = state_renormalize(*x, src).ok_or_else(truncated)?;
        }
    }

    for (d, x) in remainder.iter_mut().zip(states.iter_mut()) {
        let sym = sym_table[(*x & 0xFFF) as usize];
        *d = sym;
        let i = usize::from(sym);
        *x = state_step(*x, frequencies[i], cumulative_frequencies[i], ORDER_0_BITS);
        *x = state_renormalize(*x, src).ok_or_else(truncated)?;
    }

    Ok(())
}

fn decode_order_0_generic(
    src: &mut &[u8],
    dst: &mut [u8],
    state_count: usize,
) -> Result<(), CramError> {
    let truncated = || CramError::Truncated { context: "rans_nx16 order-0 truncated" };
    let frequencies = read_frequencies_0(src)?;
    let cumulative_frequencies = build_cumulative_frequencies(&frequencies);
    let sym_table = build_symbol_table_nx16(&cumulative_frequencies);
    let mut states = read_states(src, state_count)?;

    for chunk in dst.chunks_mut(states.len()) {
        for (d, state) in chunk.iter_mut().zip(states.iter_mut()) {
            let f = state_cumulative_frequency(*state, ORDER_0_BITS);
            debug_assert!((f as usize) < sym_table.len());
            #[allow(clippy::indexing_slicing, reason = "bounds checked by debug_assert above")]
            let sym = sym_table[f as usize];
            *d = sym;
            let i = usize::from(sym);
            *state = state_step(
                *state,
                *frequencies.get(i).unwrap_or(&0),
                *cumulative_frequencies.get(i).unwrap_or(&0),
                ORDER_0_BITS,
            );
            *state = state_renormalize(*state, src).ok_or_else(truncated)?;
        }
    }

    Ok(())
}

fn read_frequencies_0(src: &mut &[u8]) -> Result<[u32; ALPHABET_SIZE], CramError> {
    let alphabet = read_alphabet(src)?;
    let mut frequencies = [0u32; ALPHABET_SIZE];

    for (i, frequency) in alphabet.iter().zip(&mut frequencies) {
        if *i {
            *frequency = read_uint7(src)?;
        }
    }

    normalize_frequencies(&mut frequencies, ORDER_0_BITS)?;
    Ok(frequencies)
}

// r[impl cram.codec.normalize_checked]
fn normalize_frequencies(
    frequencies: &mut [u32; ALPHABET_SIZE],
    bits: u32,
) -> Result<(), CramError> {
    let sum: u32 = frequencies.iter().copied().try_fold(0u32, |acc, f| acc.checked_add(f)).ok_or(
        CramError::FrequencyNormalizationOverflow {
            sum: frequencies.iter().copied().map(u64::from).sum(),
        },
    )?;

    if sum == 0 || sum == (1 << bits) {
        return Ok(());
    }

    let mut running = sum;
    let mut shift = 0u32;
    while running < (1 << bits) {
        running = running.checked_mul(2).ok_or_else(|| {
            CramError::FrequencyNormalizationOverflow { sum: u64::from(sum) << shift }
        })?;
        shift = shift
            .checked_add(1)
            .ok_or_else(|| CramError::FrequencyNormalizationOverflow { sum: u64::from(sum) })?;
    }

    for f in frequencies {
        *f <<= shift;
    }

    Ok(())
}

fn build_cumulative_frequencies(frequencies: &[u32; ALPHABET_SIZE]) -> [u32; ALPHABET_SIZE] {
    let mut cumulative_frequencies = [0u32; ALPHABET_SIZE];
    let mut f = 0u32;

    for (i, cf) in cumulative_frequencies.iter_mut().enumerate().skip(1) {
        let prev = i.checked_sub(1).expect("skip(1) guarantees i ≥ 1");
        f = f.saturating_add(*frequencies.get(prev).unwrap_or(&0));
        *cf = f;
    }

    cumulative_frequencies
}

/// Pre-computed symbol table for order-0: `table[f]` returns the symbol
/// whose cumulative frequency range contains `f` (0 ≤ f < 4096).
/// Eliminates the per-state linear scan in `cumulative_frequencies_symbol`.
#[allow(
    clippy::indexing_slicing,
    reason = "sym ≤ 254 (loop guard), sym+1 ≤ 255 < ALPHABET_SIZE=256"
)]
fn build_symbol_table_nx16(cum: &[u32; ALPHABET_SIZE]) -> [u8; 4096] {
    let mut table = [0u8; 4096];
    build_symbol_table_nx16_into(cum, &mut table);
    table
}

/// In-place form of [`build_symbol_table_nx16`] for callers with a
/// reusable buffer (the order-1 path: 256 contexts × 4096 entries =
/// 1 MiB, repopulated per block).
#[allow(
    clippy::indexing_slicing,
    reason = "sym ≤ 254 (loop guard), sym+1 ≤ 255 < ALPHABET_SIZE=256"
)]
fn build_symbol_table_nx16_into(cum: &[u32; ALPHABET_SIZE], table: &mut [u8; 4096]) {
    let mut sym = 0u8;
    for (f, entry) in (0u32..4096).zip(table.iter_mut()) {
        while sym < 255 && f >= *cum.get(usize::from(sym).wrapping_add(1)).unwrap_or(&u32::MAX) {
            sym = sym.wrapping_add(1);
        }
        *entry = sym;
    }
}

fn decode_order_1_with_buf(
    src: &mut &[u8],
    dst: &mut [u8],
    state_count: usize,
    buf: &mut Nx16Order1Buf,
) -> Result<(), CramError> {
    decode_order_1_at(Level::new(), src, dst, state_count, buf)
}

/// Packs the rows of the active contexts — the only ones a stream can reach
/// when context 0 is among them, since every decoded symbol is in the
/// alphabet — and takes the packed kernels if all of them pack; otherwise
/// the unpacked scalar path, with the inactive rows zeroed so a reused
/// buffer's earlier blocks cannot leak into this one.
// r[impl cram.codec.simd_dispatch+3]
fn decode_order_1_at(
    level: Level,
    src: &mut &[u8],
    dst: &mut [u8],
    state_count: usize,
    buf: &mut Nx16Order1Buf,
) -> Result<(), CramError> {
    let buf = buf.tables.get_or_insert_with(Order1Tables::new);
    let (bits, alphabet) = read_frequencies_1(src, &mut buf.frequencies)?;

    let packs = alphabet[0]
        && alphabet
            .iter()
            .zip(buf.frequencies.iter())
            .zip(buf.packed.iter_mut())
            .filter(|((active, _), _)| **active)
            .all(|((_, freqs), row)| pack_row(freqs.iter().copied(), bits, row));
    if packs {
        match state_count {
            32 => {
                let mut states = read_state_array::<32>(src)?;
                return dispatch!(level, simd => decode_o1_packed(simd, src, dst, &buf.packed, bits, &mut states));
            }
            4 => {
                let mut states = read_state_array::<4>(src)?;
                return decode_o1_packed_scalar(src, dst, &buf.packed, bits, &mut states);
            }
            _ => {}
        }
    }

    // r[impl cram.codec.order1_table_reset]
    for (active, row) in alphabet.iter().zip(buf.frequencies.iter_mut()) {
        if !active {
            row.fill(0);
        }
    }
    decode_order_1_unpacked(src, dst, state_count, buf, bits)
}

/// The scalar order-1 decode over unpacked tables: the path for frequency
/// tables that do not pack, and the oracle the packed kernels are tested
/// against.
fn decode_order_1_unpacked(
    src: &mut &[u8],
    dst: &mut [u8],
    state_count: usize,
    buf: &mut Order1Tables,
    bits: u32,
) -> Result<(), CramError> {
    build_cumulative_frequencies_1_into(&buf.frequencies, &mut buf.cumulative_frequencies);
    // Pre-compute per-context symbol-decode tables. A row of zero
    // frequencies (every inactive context) decodes every slot to 255, which
    // is what the scan would compute, without the scan.
    for ((freqs, cum), table) in
        buf.frequencies.iter().zip(buf.cumulative_frequencies.iter()).zip(buf.sym_tables.iter_mut())
    {
        if freqs.iter().all(|&f| f == 0) {
            table.fill(255);
        } else {
            build_symbol_table_nx16_into(cum, table);
        }
    }

    let states = &mut buf.states;
    states.clear();
    (0..state_count)
        .try_for_each(|_| -> Option<()> {
            states.push(read_u32_le(src)?);
            Some(())
        })
        .ok_or_else(|| CramError::Truncated { context: "rans_nx16 order-1 states" })?;

    let prev_syms = &mut buf.prev_syms;
    prev_syms.clear();
    prev_syms.resize(state_count, 0);

    let chunk_size = dst
        .len()
        .checked_div(state_count)
        .ok_or_else(|| CramError::Truncated { context: "rans_nx16 order-1 zero state count" })?;

    // r[impl cram.codec.rans_nx16_bits_validation]
    #[allow(
        clippy::indexing_slicing,
        reason = "k/l ≤ 255, state/prev_sym indices from enumerate < state_count; \
                  f < 4096 (bits ≤ 12 validated by read_frequencies_1)"
    )]
    for i in 0..chunk_size {
        for (j, (state, prev_sym)) in states.iter_mut().zip(prev_syms.iter_mut()).enumerate() {
            let k = usize::from(*prev_sym);
            let f = state_cumulative_frequency(*state, bits);
            let sym = buf.sym_tables[k][f as usize];

            let out_idx =
                j.checked_mul(chunk_size).and_then(|v| v.checked_add(i)).ok_or_else(|| {
                    CramError::Truncated { context: "rans_nx16 order-1 index overflow" }
                })?;
            *dst.get_mut(out_idx)
                .ok_or_else(|| CramError::Truncated { context: "rans_nx16 order-1 output" })? = sym;

            let l = usize::from(sym);
            *state =
                state_step(*state, buf.frequencies[k][l], buf.cumulative_frequencies[k][l], bits);
            *state = state_renormalize(*state, src)
                .ok_or_else(|| CramError::Truncated { context: "rans_nx16 order-1 renormalize" })?;
            *prev_sym = sym;
        }
    }

    let last_chunk_start = chunk_size.checked_mul(state_count).ok_or_else(|| {
        CramError::Truncated { context: "rans_nx16 order-1 last chunk offset overflow" }
    })?;
    debug_assert!(last_chunk_start <= dst.len());
    #[allow(clippy::indexing_slicing, reason = "last_chunk_start bounded by checked_mul above")]
    let last_chunk = &mut dst[last_chunk_start..];
    if !last_chunk.is_empty() {
        let mut state = *states
            .last()
            .ok_or_else(|| CramError::Truncated { context: "rans_nx16 order-1 last state" })?;
        let mut prev_sym = *prev_syms
            .last()
            .ok_or_else(|| CramError::Truncated { context: "rans_nx16 order-1 last prev_sym" })?;

        for d in last_chunk {
            let k = usize::from(prev_sym);
            debug_assert!(k < ALPHABET_SIZE, "prev_sym {prev_sym} out of range");
            let f = state_cumulative_frequency(state, bits);
            #[allow(clippy::indexing_slicing, reason = "k < 256, f < 4096")]
            let sym = buf.sym_tables[k][f as usize];
            *d = sym;
            let l = usize::from(sym);
            debug_assert!(l < ALPHABET_SIZE, "sym {sym} out of range");
            #[allow(clippy::indexing_slicing, reason = "k < 256, l < 256")]
            {
                state = state_step(
                    state,
                    buf.frequencies[k][l],
                    buf.cumulative_frequencies[k][l],
                    bits,
                );
            }
            state = state_renormalize(state, src).ok_or_else(|| CramError::Truncated {
                context: "rans_nx16 order-1 last renormalize",
            })?;
            prev_sym = sym;
        }
    }

    Ok(())
}

fn build_cumulative_frequencies_1_into(
    frequencies: &Frequencies1,
    cf: &mut CumulativeFrequencies1,
) {
    for (f, g) in frequencies.iter().zip(cf.iter_mut()) {
        *g = build_cumulative_frequencies(f);
    }
}

/// Reads the order-1 frequency tables into the rows of the active contexts
/// (every other row is left as it was) and returns `bits` and the alphabet.
fn read_frequencies_1(
    src: &mut &[u8],
    frequencies: &mut Frequencies1,
) -> Result<(u32, [bool; ALPHABET_SIZE]), CramError> {
    let n =
        read_u8(src).ok_or_else(|| CramError::Truncated { context: "rans_nx16 freq1 header" })?;
    let bits = u32::from(n >> 4);
    // rANS Nx16 symbol tables are always 4096 entries (1 << 12).
    // bits > 12 would cause state_cumulative_frequency to produce
    // indices outside the table bounds.
    // r[impl cram.codec.rans_nx16_bits_validation]
    if bits > ORDER_0_BITS {
        return Err(CramError::InvalidRansNx16Bits { bits });
    }
    let is_compressed = (n & 0x01) != 0;

    if is_compressed {
        let uncompressed_size = read_uint7(src)? as usize;
        let compressed_size = read_uint7(src)? as usize;
        let mut compressed_data = split_off(src, compressed_size)?;
        // r[impl io.fuzz.alloc_limits]
        super::reader::check_alloc_size(uncompressed_size, "rans_nx16 freq1 compressed")?;
        let mut tmp = vec![0u8; uncompressed_size];
        decode_order_0(&mut compressed_data, &mut tmp, 4)?;
        let alphabet = read_frequencies_1_inner(&mut &tmp[..], frequencies, bits)?;
        Ok((bits, alphabet))
    } else {
        let alphabet = read_frequencies_1_inner(src, frequencies, bits)?;
        Ok((bits, alphabet))
    }
}

#[allow(
    clippy::indexing_slicing,
    reason = "ctx_idx/sym_idx come from enumerate() over [bool;256], so ≤ 255 < ALPHABET_SIZE=256"
)]
fn read_frequencies_1_inner(
    src: &mut &[u8],
    frequencies: &mut Frequencies1,
    bits: u32,
) -> Result<[bool; ALPHABET_SIZE], CramError> {
    let alphabet = read_alphabet(src)?;

    for (ctx_idx, ctx_active) in alphabet.iter().enumerate() {
        if !*ctx_active {
            continue;
        }

        let fs = &mut frequencies[ctx_idx];
        // r[impl cram.codec.order1_table_reset]
        // Symbols the stream skips (outside the alphabet, or in a run of
        // zeros) are zero, not whatever an earlier block left in a reused row.
        fs.fill(0);
        let mut sym_iter = alphabet.iter().enumerate().filter(|(_, b)| **b).peekable();

        while let Some((sym_idx, _)) = sym_iter.next() {
            let f = read_uint7(src)?;
            fs[sym_idx] = f;

            if f == 0 {
                let n = read_u8(src)
                    .ok_or_else(|| CramError::Truncated { context: "rans_nx16 freq1 inner run" })?
                    as usize;
                for _ in 0..n {
                    let _ = sym_iter.next();
                }
            }
        }

        normalize_frequencies(fs, bits)?;
    }

    Ok(alphabet)
}

// ── Stripe transform ─────────────────────────────────────────────────

fn decode_stripe_with_buf(
    src: &mut &[u8],
    uncompressed_size: usize,
    buf: &mut Nx16Order1Buf,
) -> Result<Vec<u8>, CramError> {
    let chunk_count = read_u8(src)
        .ok_or_else(|| CramError::Truncated { context: "rans_nx16 stripe chunk count" })?
        as usize;
    if chunk_count == 0 {
        return Err(CramError::RansStripeZeroChunks);
    }

    let compressed_sizes: Vec<usize> =
        (0..chunk_count).map(|_| read_uint7(src).map(|n| n as usize)).collect::<Result<_, _>>()?;

    let q = uncompressed_size
        .checked_div(chunk_count)
        .ok_or_else(|| CramError::RansStripeZeroChunks)?;
    let r = uncompressed_size
        .checked_rem(chunk_count)
        .ok_or_else(|| CramError::RansStripeZeroChunks)?;
    let uncompressed_sizes: Vec<usize> =
        (0..chunk_count).map(|i| if r > i { q.saturating_add(1) } else { q }).collect();

    let chunks: Vec<Vec<u8>> = compressed_sizes
        .iter()
        .zip(&uncompressed_sizes)
        .map(|(&cs, &us)| {
            let sub = split_off(src, cs)?;
            decode_with_buf(sub, us, buf)
        })
        .collect::<Result<_, _>>()?;

    let mut dst = vec![0u8; uncompressed_size];
    for (i, chunk) in chunks.iter().enumerate() {
        for (j, &s) in chunk.iter().enumerate() {
            let idx =
                j.checked_mul(chunk_count).and_then(|v| v.checked_add(i)).ok_or_else(|| {
                    CramError::Truncated { context: "rans_nx16 stripe index overflow" }
                })?;
            if let Some(d) = dst.get_mut(idx) {
                *d = s;
            }
        }
    }

    Ok(dst)
}

// ── Bit-pack transform ──────────────────────────────────────────────

struct BitPackContext {
    symbol_count: usize,
    mapping_table: Vec<u8>,
    uncompressed_size: usize,
}

fn read_bit_pack_context(
    src: &mut &[u8],
    uncompressed_size: usize,
) -> Result<(BitPackContext, usize), CramError> {
    let symbol_count = read_u8(src)
        .ok_or_else(|| CramError::Truncated { context: "rans_nx16 bit_pack symbol count" })?
        as usize;
    if symbol_count == 0 {
        return Err(CramError::RansBitPackZeroSymbols);
    }
    let mapping_table = split_off(src, symbol_count)?.to_vec();
    let packed_len = read_uint7(src)? as usize;

    Ok((BitPackContext { symbol_count, mapping_table, uncompressed_size }, packed_len))
}

fn apply_bit_unpack(src: &[u8], ctx: &BitPackContext) -> Result<Vec<u8>, CramError> {
    let mut dst = vec![0u8; ctx.uncompressed_size];

    match ctx.symbol_count {
        1 => {
            let sym = *ctx
                .mapping_table
                .first()
                .ok_or_else(|| CramError::Truncated { context: "bit_pack mapping" })?;
            dst.fill(sym);
        }
        2 => unpack(src, &ctx.mapping_table, 8, &mut dst),
        3..=4 => unpack(src, &ctx.mapping_table, 4, &mut dst),
        5..=16 => unpack(src, &ctx.mapping_table, 2, &mut dst),
        n => {
            return Err(CramError::RansBitPackTooManySymbols { symbol_count: n });
        }
    }

    Ok(dst)
}

fn unpack(src: &[u8], mapping_table: &[u8], chunk_size: usize, dst: &mut [u8]) {
    let bits = u8::BITS as usize;
    let shift = bits
        .checked_div(chunk_size)
        .expect("chunk_size is always 2, 4, or 8 from the match in apply_bit_unpack");
    let mask: u8 = (1u8 << shift).wrapping_sub(1);

    for (mut s, chunk) in src.iter().copied().zip(dst.chunks_mut(chunk_size)) {
        for d in chunk {
            let idx = usize::from(s & mask);
            *d = mapping_table.get(idx).copied().unwrap_or(0);
            s >>= shift;
        }
    }
}

// ── RLE transform ────────────────────────────────────────────────────

struct RleContext {
    rle_meta: Vec<u8>,
    output_len: usize,
}

fn read_rle_context(
    src: &mut &[u8],
    state_count: usize,
    uncompressed_size: usize,
) -> Result<(RleContext, usize), CramError> {
    let header = read_uint7(src)?;
    let context_size = (header >> 1) as usize;
    // r[impl io.fuzz.alloc_limits]
    super::reader::check_alloc_size(context_size, "rans_nx16 rle context")?;
    let is_compressed = (header & 0x01) == 0;

    let rle_encoded_len = read_uint7(src)? as usize;

    let rle_meta = if is_compressed {
        let compressed_size = read_uint7(src)? as usize;
        let mut buf = split_off(src, compressed_size)?;
        let mut tmp = vec![0u8; context_size];
        decode_order_0(&mut buf, &mut tmp, state_count)?;
        tmp
    } else {
        split_off(src, context_size)?.to_vec()
    };

    Ok((RleContext { rle_meta, output_len: uncompressed_size }, rle_encoded_len))
}

#[allow(
    clippy::indexing_slicing,
    reason = "sym is u8 so usize::from(sym) ≤ 255 < ALPHABET_SIZE=256"
)]
fn apply_rle(src: &[u8], ctx: &RleContext) -> Result<Vec<u8>, CramError> {
    let mut meta_cur: &[u8] = &ctx.rle_meta;

    let rle_alphabet = read_rle_alphabet(&mut meta_cur)?;

    let mut dst = vec![0u8; ctx.output_len];
    let mut dst_iter = dst.iter_mut();
    let mut src_iter = src.iter();

    while let Some(d) = dst_iter.next() {
        let &sym =
            src_iter.next().ok_or_else(|| CramError::Truncated { context: "rans_nx16 rle src" })?;
        *d = sym;

        if rle_alphabet[usize::from(sym)] {
            let len = read_uint7(&mut meta_cur)? as usize;
            for e in dst_iter.by_ref().take(len) {
                *e = sym;
            }
        }
    }

    Ok(dst)
}

#[allow(
    clippy::indexing_slicing,
    reason = "sym is u8 so usize::from(sym) ≤ 255 < ALPHABET_SIZE=256"
)]
fn read_rle_alphabet(src: &mut &[u8]) -> Result<[bool; ALPHABET_SIZE], CramError> {
    let truncated = || CramError::Truncated { context: "rans_nx16 rle alphabet" };
    let mut alphabet = [false; ALPHABET_SIZE];

    let n = read_u8(src).ok_or_else(truncated)? as usize;
    let symbol_count = if n == 0 { ALPHABET_SIZE } else { n };

    for _ in 0..symbol_count {
        let sym = read_u8(src).ok_or_else(truncated)?;
        alphabet[usize::from(sym)] = true;
    }

    Ok(alphabet)
}

#[cfg(test)]
#[allow(clippy::cast_possible_truncation, reason = "test code")]
#[allow(clippy::arithmetic_side_effects, reason = "test code")]
mod tests {
    use super::*;
    use hegel::prelude::*;

    // r[verify cram.codec.rans_nx16]

    #[test]
    fn rans_stripe_zero_chunks_returns_error() {
        // Build a byte stream with FLAG_STRIPE set and chunk_count = 0
        let src = vec![
            FLAG_STRIPE, // flags = STRIPE
            4u8,         // uncompressed_size = 4 (uint7)
            0u8,         // chunk_count = 0 → triggers RansStripeZeroChunks
        ];

        let err = decode(&src, 0).unwrap_err();
        assert!(matches!(err, CramError::RansStripeZeroChunks));
    }

    #[test]
    fn rans_bit_pack_zero_symbols_returns_error() {
        // Build a stream with FLAG_PACK and symbol_count = 0
        let src = vec![
            FLAG_PACK, // flags = PACK
            5u8,       // uncompressed_size = 5 (uint7)
            0u8,       // symbol_count = 0 → triggers RansBitPackZeroSymbols
        ];

        let err = decode(&src, 0).unwrap_err();
        assert!(matches!(err, CramError::RansBitPackZeroSymbols));
    }

    #[test]
    fn rans_bit_pack_too_many_symbols_returns_error() {
        // symbol_count = 17 (> 16) — triggers RansBitPackTooManySymbols.
        // After reading symbol_count, the code reads mapping_table (symbol_count bytes)
        // and packed_len (uint7). We need enough bytes to get past those reads.
        let symbol_count: u8 = 17;
        let packed_len: u8 = 0; // uint7 = 0
        let mut src = Vec::new();
        src.push(FLAG_PACK); // flags = PACK
        // FLAG_NO_SIZE not set → read uncompressed_size as uint7
        src.push(10u8); // uncompressed_size
        src.push(symbol_count); // symbol_count = 17
        // mapping_table: 17 bytes
        src.extend(std::iter::repeat_n(b'A', symbol_count as usize));
        src.push(packed_len); // packed_len uint7

        // The rest of the stream would be consumed for the inner rANS decode, but the
        // error fires before that in apply_bit_unpack. We need the CAT flag set
        // to bypass the rANS decode and reach apply_bit_unpack directly.
        // Rebuild with CAT | PACK so we skip the rANS step.
        let mut src2 = Vec::new();
        src2.push(FLAG_PACK | FLAG_CAT); // flags = PACK | CAT
        src2.push(10u8); // uncompressed_size
        src2.push(symbol_count); // symbol_count = 17
        src2.extend(std::iter::repeat_n(b'A', symbol_count as usize));
        src2.push(packed_len); // packed_len uint7
        // CAT data: 0 bytes (packed_len = 0)

        let err = decode(&src2, 0).unwrap_err();
        assert!(
            matches!(err, CramError::RansBitPackTooManySymbols { symbol_count: 17 }),
            "expected RansBitPackTooManySymbols, got: {err:?}"
        );
    }

    #[test]
    fn decode_order_0_noodles_test_vector() {
        let src = [
            0x00, // flags = {empty}
            0x07, // uncompressed len = 7
            0x64, 0x65, 0x00, 0x6c, 0x6e, 0x6f, 0x00, 0x73, 0x00, 0x01, 0x01, 0x01, 0x01, 0x03,
            0x01, 0x00, 0x26, 0x20, 0x00, 0x00, 0xb8, 0x0a, 0x00, 0x00, 0xd8, 0x0a, 0x00, 0x00,
            0x00, 0x04, 0x00,
        ];
        assert_eq!(decode(&src, 0).unwrap(), b"noodles");
    }

    #[test]
    fn decode_order_1_noodles_test_vector() {
        let src = [
            0x01, // flags = ORDER
            0x4d, // uncompressed len = 77
            0xa0, 0x00, 0x64, 0x65, 0x00, 0x6c, 0x6e, 0x6f, 0x00, 0x73, 0x00, 0x00, 0x00, 0x01,
            0x01, 0x00, 0x00, 0x01, 0x01, 0x00, 0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x01, 0x00,
            0x02, 0x00, 0x01, 0x0f, 0x00, 0x02, 0x01, 0x00, 0x01, 0x01, 0x0f, 0x00, 0x02, 0x00,
            0x03, 0x0f, 0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x02, 0x0f, 0x00, 0x00, 0x00,
            0x05, 0x10, 0x80, 0x72, 0x60, 0x00, 0x80, 0x8b, 0x5f, 0x00, 0xc0, 0xb0, 0x60, 0x00,
            0x40, 0x49, 0x39, 0x00,
        ];
        assert_eq!(
            decode(&src, 0).unwrap(),
            b"nnnnnnnnnnnnooooooooooooooooddddddddddddddllllllllllllllleeeeeeeeeessssssssss"
        );
    }

    #[test]
    fn decode_stripe_noodles_test_vector() {
        let src = [
            0x08, // flags = STRIPE
            0x07, // uncompressed len = 7
            0x04, 0x17, 0x17, 0x17, 0x15, 0x00, 0x02, 0x6c, 0x6e, 0x00, 0x01, 0x01, 0x00, 0x08,
            0x01, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00,
            0x00, 0x02, 0x65, 0x6f, 0x00, 0x01, 0x01, 0x00, 0x08, 0x01, 0x00, 0x00, 0x00, 0x01,
            0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x02, 0x6f, 0x73, 0x00,
            0x01, 0x01, 0x00, 0x00, 0x01, 0x00, 0x00, 0x08, 0x01, 0x00, 0x00, 0x80, 0x00, 0x00,
            0x00, 0x80, 0x00, 0x00, 0x00, 0x01, 0x64, 0x00, 0x01, 0x00, 0x80, 0x00, 0x00, 0x00,
            0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x02, 0x00,
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x22, 0x00, 0x81, 0x11, 0x01, 0x7f, 0x00,
        ];
        assert_eq!(decode(&src, 0).unwrap(), b"noodles");
    }

    #[test]
    fn decode_uncompressed_noodles_test_vector() {
        let src = [
            0x20, // flags = CAT
            0x07, // uncompressed len = 7
            0x6e, 0x6f, 0x6f, 0x64, 0x6c, 0x65, 0x73,
        ];
        assert_eq!(decode(&src, 0).unwrap(), b"noodles");
    }

    #[test]
    fn decode_rle_noodles_test_vector() {
        let src = [
            0x40, // flags = RLE
            0x0d, // uncompressed len = 13
            0x06, 0x06, 0x17, 0x01, 0x07, 0x6f, 0x00, 0x02, 0x01, 0x01, 0x00, 0x00, 0x01, 0x00,
            0x00, 0x0c, 0x02, 0x00, 0x00, 0x08, 0x02, 0x00, 0x00, 0x80, 0x00, 0x00, 0x64, 0x65,
            0x00, 0x6c, 0x6e, 0x6f, 0x00, 0x73, 0x00, 0x03, 0x01, 0x01, 0x01, 0x01, 0x01, 0x00,
            0x3a, 0x20, 0x00, 0x00, 0x7c, 0x20, 0x00, 0x00, 0x52, 0x01, 0x00, 0x00, 0x08, 0x04,
            0x00,
        ];
        assert_eq!(decode(&src, 0).unwrap(), b"noooooooodles");
    }

    #[test]
    fn decode_bit_packing_noodles_test_vector() {
        let src = [
            0x80, // flags = PACK
            0x07, // uncompressed len = 7
            0x06, 0x64, 0x65, 0x6c, 0x6e, 0x6f, 0x73, 0x04, 0x04, 0x05, 0x00, 0x12, 0x43, 0x00,
            0x01, 0x01, 0x01, 0x01, 0x00, 0x0c, 0x02, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x08,
            0x02, 0x00, 0x00, 0x04, 0x02, 0x00,
        ];
        assert_eq!(decode(&src, 0).unwrap(), b"noodles");
    }

    // r[verify cram.codec.alphabet_run_bounded]
    #[test]
    fn read_alphabet_run_overflow_returns_error() {
        // sym=250, then sym=251 (== prev+1) triggers a run of length 10. The run
        // would write alphabet[251..261], wrapping past 255 and corrupting
        // alphabet[0..5]. Must error instead.
        let mut src: &[u8] = &[250, 251, 10, 0];
        let err = read_alphabet(&mut src).unwrap_err();
        assert!(
            matches!(err, CramError::MalformedAlphabetRun { start: 251, len: 10 }),
            "expected MalformedAlphabetRun for run from 251 with len=10, got: {err:?}"
        );
    }

    // For all (start, len) where the first symbol is `start-1` (so the
    // run trigger fires), the decoder must accept iff `start + len <= 255`
    // and reject with `MalformedAlphabetRun` otherwise. Mirrors the
    // boundary in htscodecs `decode_alphabet`.
    #[hegel::test]
    fn read_alphabet_run_bounds_match_spec(tc: TestCase) {
        let start = tc.draw(gs::integers::<u8>().min_value(1));
        let len = tc.draw(gs::integers::<u8>());
        // Stream: [start-1, start, len, 0]. The first byte is needed
        // because the run-trigger requires `sym == prev_sym+1`.
        let prev = start.checked_sub(1).expect("start >= 1");
        let stream = [prev, start, len, 0];
        let mut cur: &[u8] = &stream;
        let result = read_alphabet(&mut cur);

        let end = u32::from(start) + u32::from(len);
        if end > 255 {
            assert!(
                matches!(result, Err(CramError::MalformedAlphabetRun { start: s, len: l })
                                 if s == start && l == len),
                "expected MalformedAlphabetRun for start={start}, len={len}, got: {result:?}",
            );
        } else {
            let alphabet = result.expect("valid bounds should not error");
            assert!(alphabet[usize::from(prev)], "alphabet[{prev}] must be set");
            if len > 0 {
                let last = start
                    .checked_add(len.checked_sub(1).expect("len>0"))
                    .expect("end ≤ 255 so add fits in u8");
                assert!(alphabet[usize::from(last)], "alphabet[{last}] must be set");
            }
        }
    }

    #[test]
    fn read_alphabet_run_exactly_to_255_ok() {
        // Canonical htscodecs encoding for alphabet {200..=255}: emit 200
        // (no run since alphabet[199] is unset), then 201 with run-len 54
        // (writes alphabet[201..254]; sym wraps to 255 at the end of the
        // inner loop). The next outer iteration writes alphabet[255] and
        // reads the terminator. Stream: [200, 201, 54, 0].
        let mut src: &[u8] = &[200, 201, 54, 0];
        let alphabet = read_alphabet(&mut src).unwrap();
        for s in 200..=255u32 {
            assert!(alphabet[s as usize], "alphabet[{s}] should be set");
        }
        assert!(!alphabet[199], "alphabet[199] should not be set");
        assert!(!alphabet[0], "alphabet[0] should not be set (no wraparound)");
    }

    // r[verify cram.codec.uint7_bounded]
    #[test]
    fn read_uint7_overflow_returns_error() {
        // 6 continuation bytes (all with high bit set) exceed the 5-iteration limit.
        let mut src: &[u8] = &[0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x00];
        let err = read_uint7(&mut src).unwrap_err();
        assert!(matches!(err, CramError::Uint7Overflow));
    }

    #[test]
    fn read_uint7_five_bytes_ok() {
        // Exactly 5 continuation bytes is the maximum valid sequence.
        let mut src: &[u8] = &[0x80, 0x80, 0x80, 0x80, 0x01];
        let result = read_uint7(&mut src);
        assert!(result.is_ok());
    }

    #[test]
    fn read_uint7_spec_vectors() {
        // Hard-coded byte sequences derived from htscodecs `var_put_u32`
        // (htslib/htscodecs/htscodecs/varint.h:206, BIG_END / MSB-first).
        // These are an independent oracle: a future refactor that flips
        // MSB↔LSB byte order in `read_uint7` would silently keep the
        // round-trip property passing (because the encoder is in the same
        // file), but would break against these fixed bytes.
        let cases: &[(u32, &[u8])] = &[
            (0, &[0x00]),
            (1, &[0x01]),
            (127, &[0x7F]),
            (128, &[0x81, 0x00]),
            (200, &[0x81, 0x48]),
            (16_383, &[0xFF, 0x7F]),
            (16_384, &[0x81, 0x80, 0x00]),
            (0x12345, &[0x84, 0xC6, 0x45]),
            ((1u32 << 28) - 1, &[0xFF, 0xFF, 0xFF, 0x7F]),
            (1u32 << 28, &[0x81, 0x80, 0x80, 0x80, 0x00]),
            (u32::MAX, &[0x8F, 0xFF, 0xFF, 0xFF, 0x7F]),
        ];
        for (val, encoded) in cases {
            let mut cur: &[u8] = encoded;
            let decoded = read_uint7(&mut cur)
                .unwrap_or_else(|e| panic!("decode of {encoded:02x?} failed: {e:?}"));
            assert_eq!(decoded, *val, "decoded value mismatch for {encoded:02x?}");
            assert!(
                cur.is_empty(),
                "decoder consumed wrong byte count for {val} encoded as {encoded:02x?}",
            );
        }
    }

    // r[verify cram.codec.rans_nx16]
    // Direct test of state_renormalize. Models the spec exactly:
    // - if state >= 1<<15, no bytes consumed, state unchanged.
    // - else, repeatedly shift-left-16 and OR a u16 LE from src
    //   until state >= 1<<15.
    #[hegel::test]
    fn state_renormalize_matches_spec(tc: TestCase) {
        let initial_state = tc.draw(gs::integers::<u32>());
        let bytes = tc.draw(gs::binary().max_size(16));
        let mut src: &[u8] = &bytes;
        let initial_src_len = src.len();
        let result = state_renormalize(initial_state, &mut src);

        if initial_state >= (1 << 15) {
            // Fast path: no bytes consumed, state unchanged.
            assert_eq!(result, Some(initial_state));
            assert_eq!(src.len(), initial_src_len);
            return;
        }

        // Slow path: build the expected result by replaying the spec.
        let mut expected_state = initial_state;
        let mut expected_src: &[u8] = &bytes;
        while expected_state < (1 << 15) {
            let Some((head, rest)) = expected_src.split_first_chunk::<2>() else {
                assert_eq!(result, None, "expected None when src exhausted mid-renorm");
                return;
            };
            let lo = u32::from(u16::from_le_bytes(*head));
            expected_state = expected_state.wrapping_shl(16) | lo;
            expected_src = rest;
        }
        assert_eq!(result, Some(expected_state));
        assert_eq!(src.len(), expected_src.len());
    }

    #[hegel::test]
    fn read_uint7_roundtrip(tc: TestCase) {
        let val = tc.draw(gs::integers::<u32>());
        let mut stream = Vec::with_capacity(256);
        encode_uint7_prv(&mut stream, val);
        let mut cur: &[u8] = &stream;
        let decoded = read_uint7(&mut cur).unwrap();
        assert_eq!(decoded, val);
        // Encoder and decoder must consume exactly the same number of bytes.
        assert!(cur.is_empty(), "undecoded trailing bytes");
    }

    #[hegel::test]
    fn read_uint7_roundtrip_max_continuation(tc: TestCase) {
        // Values requiring 5 continuation bytes (the maximum).
        let val = tc.draw(gs::integers::<u32>().min_value(1 << 28));
        let mut stream = Vec::with_capacity(256);
        encode_uint7_prv(&mut stream, val);
        // Must produce exactly 5 bytes (all continuation except last).
        assert_eq!(stream.len(), 5, "max-continuation encodes to 5 bytes");
        let mut cur: &[u8] = &stream;
        let decoded = read_uint7(&mut cur).unwrap();
        assert_eq!(decoded, val);
    }

    // r[verify cram.codec.normalize_checked]
    // r[verify cram.codec.rans_nx16_bits_validation]
    #[test]
    fn invalid_bits_returns_error() {
        // rANS Nx16 symbol tables are 4096 entries (1 << 12). If the
        // frequency-table header specifies bits != 12, state_cumulative_frequency
        // would produce indices outside the table bounds, causing a panic.
        // This validates that read_frequencies_1 rejects invalid bits early.
        //
        // Crafted input: flags=ORDER|NO_SIZE (order-1, 4 states, skip stream
        // size), freq header byte with bits=13 (0xD0), single-symbol alphabet
        // (sym=0), frequency=1, 4 states. Before the fix the decode loop
        // would panic with index 6007 in a 4096-element table.
        let stream: &[u8] = &[
            0x11, // flags: ORDER | NO_SIZE
            0xD0, // freq header: bits=13 (0xD >> 4), not compressed
            0x00, // alphabet: sym=0
            0x00, // alphabet: sentinel
            0x01, // freq[0] = 1 (uint7)
            // 4 u32 LE states; first state's lower 13 bits = 6007
            0x77, 0x17, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x00, 0x00,
        ];
        let result = decode(stream, 4);
        assert!(
            matches!(result, Err(CramError::InvalidRansNx16Bits { bits: 13 })),
            "expected InvalidRansNx16Bits, got: {result:?}"
        );
    }

    #[test]
    fn normalize_frequencies_overflow_returns_error() {
        // Fill frequencies so sum overflows u32.
        let mut frequencies = [0u32; ALPHABET_SIZE];
        // 256 * 0x01_000_000 = 0x1_0000_0000 which overflows u32.
        frequencies.fill(0x0100_0000);
        let result = normalize_frequencies(&mut frequencies, ORDER_0_BITS);
        assert!(
            matches!(result, Err(CramError::FrequencyNormalizationOverflow { .. })),
            "expected overflow error, got: {result:?}"
        );
    }

    #[test]
    fn build_symbol_table_matches_linear_scan() {
        // Verify the pre-computed 4096-entry table matches
        // cumulative_frequencies_symbol for all 4096 possible f values.
        let freq: [u32; 256] = {
            let mut f = [0u32; 256];
            // A non-trivial frequency distribution
            f[0] = 100;
            f[1] = 200;
            f[2] = 50;
            f[3] = 3646; // sum = 3996, leaving 100 for rest
            f[4] = 100;
            f
        };
        let cum = build_cumulative_frequencies(&freq);
        let table = build_symbol_table_nx16(&cum);

        for f_val in 0u32..4096 {
            let expected = cumulative_frequencies_symbol(&cum, f_val);
            let actual = table[f_val as usize];
            assert_eq!(
                actual, expected,
                "mismatch at f={f_val}: table gives {actual}, linear scan gives {expected}"
            );
        }
    }

    #[test]
    fn build_symbol_table_matches_linear_scan_uniform() {
        // Also test with uniformly distributed frequencies (every symbol same freq)
        let freq = [16u32; 256]; // 256 * 16 = 4096
        let cum = build_cumulative_frequencies(&freq);
        let table = build_symbol_table_nx16(&cum);

        for f_val in 0u32..4096 {
            let expected = cumulative_frequencies_symbol(&cum, f_val);
            let actual = table[f_val as usize];
            assert_eq!(actual, expected);
        }
    }

    // For arbitrary cumulative-frequency distributions summing to 4096,
    // the precomputed table must agree with the linear-scan oracle on
    // every f in [0, 4096).
    #[hegel::test]
    fn build_symbol_table_agrees_with_linear_scan(tc: TestCase) {
        let seeds =
            tc.draw(gs::vecs(gs::integers::<u32>().max_value(4096)).min_size(8).max_size(8));
        // Build a frequency table from the seed weights, normalized to sum 4096.
        let total: u64 = seeds.iter().map(|&s| u64::from(s)).sum();
        tc.assume(total > 0);
        let mut freq = [0u32; 256];
        for (i, &s) in seeds.iter().enumerate() {
            let scaled = (u64::from(s) * 4096 / total) as u32;
            freq[i % 256] = freq[i % 256].saturating_add(scaled);
        }
        // Normalize residual to symbol 0.
        let sum: u32 = freq.iter().sum();
        if sum < 4096 {
            freq[0] = freq[0].saturating_add(4096 - sum);
        }
        // Skip if normalization overflow — only valid distributions matter.
        tc.assume(freq.iter().sum::<u32>() == 4096);

        let cum = build_cumulative_frequencies(&freq);
        let table = build_symbol_table_nx16(&cum);

        // Verify the in-place form populates identically.
        let mut table_into = [0u8; 4096];
        build_symbol_table_nx16_into(&cum, &mut table_into);
        assert_eq!(table, table_into);

        for f_val in 0u32..4096 {
            let expected = cumulative_frequencies_symbol(&cum, f_val);
            let actual = table[f_val as usize];
            assert_eq!(actual, expected, "mismatch at f={f_val}");
        }
    }

    #[allow(clippy::cast_possible_truncation, reason = "len bounded by the generator")]
    #[test]
    fn order0_32state_and_generic_produce_same_output() {
        // Minimal 32-state order-0 stream: symbol 0 with freq=4096
        // (covers the full 12-bit range), 32 states initialized
        // to freq<<12, decodes to all-zero output. Validates that
        // the 32-state path and generic path agree on the same input.
        let mut stream = Vec::with_capacity(256);
        stream.push(0x04); // flags: N32
        stream.push(4); // uncompressed_size = 4 (uint7)
        stream.push(0); // alphabet: sym=0
        stream.push(0); // alphabet terminator
        stream.push(0xA0); // freq=4096 (uint7, 2 bytes)
        stream.push(0x20);
        for _ in 0..32 {
            stream.extend_from_slice(&0x01000000u32.to_le_bytes());
        }

        let result = decode(&stream, 0).unwrap();
        assert_eq!(result, &[0, 0, 0, 0]);

        // Direct comparison: 32-state vs generic
        let mut src1: &[u8] = &stream;
        let mut src2: &[u8] = &stream;
        let flags = read_u8(&mut src1).unwrap();
        let _ = read_u8(&mut src2).unwrap();
        assert_eq!(flags & FLAG_N32, FLAG_N32);
        let _ = read_uint7(&mut src1).unwrap();
        let _ = read_uint7(&mut src2).unwrap();
        let mut dst1 = vec![0u8; 4];
        let mut dst2 = vec![0u8; 4];
        decode_order_0_32state(&mut src1, &mut dst1).unwrap();
        decode_order_0_generic(&mut src2, &mut dst2, 32).unwrap();
        assert_eq!(dst1, dst2);
    }

    #[allow(
        clippy::cast_possible_truncation,
        clippy::arithmetic_side_effects,
        reason = "val masked to 7 bits; n bounded to ≤5 by loop over 32-bit value shifted by 7"
    )]
    fn encode_uint7_prv(stream: &mut Vec<u8>, val: u32) {
        if val == 0 {
            stream.push(0);
            return;
        }
        // Encode into temp buffer LSB-first, then emit in reverse (MSB-first).
        let mut tmp = [0u8; 5];
        let mut n = 0;
        let mut v = val;
        while v > 0 {
            tmp[n] = v as u8 & 0x7F;
            v >>= 7;
            n += 1;
        }
        // MSB-first: all bytes except the last get continuation bit set.
        while n > 1 {
            n -= 1;
            stream.push(tmp[n] | 0x80);
        }
        stream.push(tmp[0]);
    }

    #[test]
    fn neon_simd_handles_len_32() {
        let len = 32;
        let mut stream = Vec::with_capacity(256);
        stream.push(FLAG_N32);
        encode_uint7_prv(&mut stream, len as u32);
        stream.extend_from_slice(&[0, 0]);
        stream.push(0x80);
        stream.push(0x20);
        for _ in 0..32 {
            stream.extend_from_slice(&0x01000000u32.to_le_bytes());
        }

        // Scalar path (bypass SIMD dispatch)
        let mut cur: &[u8] = &stream;
        read_u8(&mut cur).unwrap();
        read_uint7(&mut cur).unwrap();
        let mut dst = vec![0u8; len];
        decode_order_0_generic(&mut cur, &mut dst, 32).unwrap();
        assert_eq!(dst, vec![0u8; 32], "scalar path failed");

        // Full decode (may use SIMD)
        let result = decode(&stream, 0).unwrap();
        assert_eq!(result, vec![0u8; 32], "SIMD path failed");
    }

    #[test]
    fn neon_fallback_handles_len_128_direct() {
        let len = 128;
        let mut stream = Vec::with_capacity(256);
        stream.push(FLAG_N32);
        encode_uint7_prv(&mut stream, len as u32);
        stream.extend_from_slice(&[0, 0]);
        stream.push(0x80);
        stream.push(0x20);
        for _ in 0..32 {
            stream.extend_from_slice(&0x01000000u32.to_le_bytes());
        }

        // Direct call: decode_order_0_32state bypasses decode()'s wrapper
        let mut cur: &[u8] = &stream;
        read_u8(&mut cur).unwrap();
        read_uint7(&mut cur).unwrap();
        let mut dst = vec![0u8; len];
        decode_order_0_32state(&mut cur, &mut dst).unwrap();
        assert!(dst.iter().all(|&b| b == 0), "decode_order_0_32state produced non-zero output");
    }

    #[hegel::test]
    fn simd_matches_scalar_order0_32state(tc: TestCase) {
        let len = tc.draw(gs::integers::<usize>().max_value(1023));
        // Stream decoding to `len` zero bytes. Symbol 0 has freq=4096
        // (covers the full 12-bit range), all other symbols freq=0.
        // With freq[sym0]=4096 the state is invariant — no renorm needed.
        let mut stream = Vec::with_capacity(256);
        stream.push(FLAG_N32);
        encode_uint7_prv(&mut stream, len as u32);
        // Alphabet: sym 0 only
        stream.extend_from_slice(&[0, 0]);
        // Frequency: 4096 for sym 0 (2-byte uint7), 0 for all others (single 0 byte each)
        stream.push(0x80);
        stream.push(0x20); // freq[0] = 4096 (uint7: 0x1000 → 0x80, 0x20)
        // 32 initial states
        for _ in 0..32 {
            stream.extend_from_slice(&0x0100_0000u32.to_le_bytes());
        }

        let simd_result = decode(&stream, 0).unwrap();
        assert_eq!(simd_result.len(), len);

        let mut cur: &[u8] = &stream;
        read_u8(&mut cur).unwrap();
        read_uint7(&mut cur).unwrap();
        let mut scalar_dst = vec![0u8; len];
        decode_order_0_generic(&mut cur, &mut scalar_dst, 32).unwrap();
        assert_eq!(simd_result, scalar_dst);
    }

    // r[verify cram.codec.simd_dispatch+3]
    // r[verify io.simd_portable]
    /// Every level against the generic scalar decoder, including a partial
    /// last chunk: two symbols at 2048 each, so each state's low 12 bits pick
    /// the symbol, arbitrary initial states and arbitrary renorm bytes, which
    /// may run out, in which case every level must fail exactly where the
    /// scalar decoder does.
    #[hegel::test]
    fn every_level_matches_scalar_order0_32state(tc: TestCase) {
        // At least one whole 32-byte chunk, so the vector step always runs,
        // plus a partial one.
        let chunks = tc.draw(gs::integers::<usize>().min_value(1).max_value(8));
        let len = chunks * 32 + tc.draw(gs::integers::<usize>().max_value(31));
        // Small enough that each state halves below 1 << 15 within a few
        // chunks, so renormalization — and running out of bytes — happens.
        let states =
            tc.draw(gs::vecs(gs::integers::<u32>().max_value(1 << 20)).min_size(32).max_size(32));
        // A state renormalizes at most once per 16 steps, so 256 bytes cover
        // every renorm of 9 chunks; a short cut covers truncation.
        let mut renorm_bytes = tc.draw(gs::binary().min_size(256).max_size(512));
        if tc.draw(gs::booleans()) {
            renorm_bytes.truncate(tc.draw(gs::integers::<usize>().max_value(64)));
        }
        // Alphabet {0, 1}: sym 0, sym 1, a zero-length run, the terminator.
        let mut stream = vec![0, 1, 0, 0];
        encode_uint7_prv(&mut stream, 2048);
        encode_uint7_prv(&mut stream, 2048);
        for s in &states {
            stream.extend_from_slice(&s.to_le_bytes());
        }
        stream.extend_from_slice(&renorm_bytes);

        let mut scalar_src: &[u8] = &stream;
        let mut scalar_dst = vec![0u8; len];
        let scalar = decode_order_0_generic(&mut scalar_src, &mut scalar_dst, 32);

        for level in crate::simd_levels::levels() {
            let mut src: &[u8] = &stream;
            let mut dst = vec![0u8; len];
            let simd = decode_order_0_32state_at(level, &mut src, &mut dst);
            assert_eq!(simd.is_ok(), scalar.is_ok(), "{level:?}: {simd:?} vs {scalar:?}");
            // On a truncation the partial output may differ: the kernel
            // writes a whole chunk before renormalizing, the scalar decoder
            // one byte at a time.
            if simd.is_ok() {
                assert_eq!(dst, scalar_dst, "{level:?}");
                assert_eq!(src.len(), scalar_src.len(), "{level:?} consumed differently");
            }
        }
    }

    // Exercises the renormalization path that the other two SIMD/scalar
    // properties miss. With initial state = 0x8001 (just above the 1<<15
    // renorm threshold) and freq[sym0]=2048, the first decode step
    // produces new_state = 2048 * 8 + 1 = 0x4001, below the threshold,
    // forcing a renorm read for every lane on every outer iteration.
    // The trailing renorm bytes vary across the generator's seed space.
    #[hegel::test]
    fn simd_matches_scalar_with_renorm(tc: TestCase) {
        let renorm_bytes = tc.draw(gs::binary().min_size(256).max_size(511));
        let len = tc.draw(gs::integers::<usize>().min_value(32).max_value(128));
        // Round len down to a multiple of 32 so the test exercises only
        // the full-chunk path (the remainder path is covered separately).
        let len = (len / 32).checked_mul(32).expect("len ≤ 128 → fits in usize");
        let mut stream = Vec::with_capacity(256);
        stream.push(FLAG_N32);
        encode_uint7_prv(&mut stream, len as u32);
        // Alphabet: sym 0 only.
        stream.extend_from_slice(&[0, 0]);
        // freq[0] = 2048 (uint7: 0x90, 0x00). With one symbol active
        // and sum=2048, normalize_frequencies will scale to 4096 by
        // doubling (shift=1), giving an effective freq[0] = 4096
        // — but we want renorm, so use TWO symbols with freq=2048 each.
        stream.clear();
        stream.push(FLAG_N32);
        encode_uint7_prv(&mut stream, len as u32);
        stream.push(0); // sym 0
        stream.push(1); // sym 1 (consecutive triggers run-compress, but len=0 is fine)
        stream.push(0); // run len = 0 (no run)
        stream.push(0); // terminator
        encode_uint7_prv(&mut stream, 2048); // freq[0]
        encode_uint7_prv(&mut stream, 2048); // freq[1]
        // 32 initial states all at 0x8001 — first step drops state to
        // 0x4001 (< 1<<15) and forces renorm on every lane.
        for _ in 0..32 {
            stream.extend_from_slice(&0x0000_8001u32.to_le_bytes());
        }
        // Renorm-fodder bytes: each renorm consumes 2 bytes (u16 LE).
        // Provide enough for several renorms per lane × 32 lanes.
        stream.extend_from_slice(&renorm_bytes);

        // Scalar path runs first as the oracle. Whatever scalar
        // produces (Ok+output OR Err), SIMD must match exactly. Now
        // that the SIMD dispatch propagates errors instead of falling
        // back to scalar, a SIMD-only failure surfaces as a test
        // failure rather than being silently masked.
        let mut cur: &[u8] = &stream;
        read_u8(&mut cur).unwrap();
        read_uint7(&mut cur).unwrap();
        let mut scalar_dst = vec![0u8; len];
        let scalar_result = decode_order_0_generic(&mut cur, &mut scalar_dst, 32);

        let simd_result = decode(&stream, 0);

        match (simd_result, scalar_result) {
            (Ok(simd_dst), Ok(())) => {
                assert_eq!(simd_dst, scalar_dst);
            }
            (Err(_), Err(_)) => {
                // Both paths agreed by erroring — fine.
            }
            (Ok(_), Err(scalar_err)) => {
                panic!("SIMD succeeded where scalar failed: {scalar_err:?}");
            }
            (Err(simd_err), Ok(())) => {
                panic!("SIMD failed ({simd_err:?}) where scalar succeeded — SIMD bug");
            }
        }
    }

    #[hegel::test]
    fn simd_remainder_uses_correct_lanes(tc: TestCase) {
        let len = tc.draw(gs::integers::<usize>().max_value(255));
        // Two symbols with equal frequencies so states diverge on
        // different s&0xFFF values. Initial states start large enough
        // that no renormalization is needed for ≤256 rounds.
        let mut stream = Vec::with_capacity(512);
        stream.push(FLAG_N32);
        encode_uint7_prv(&mut stream, len as u32);
        // Alphabet: sym 5 and sym 100. Non-consecutive to avoid the
        // run-compression path in read_alphabet (5→100 is not a run).
        stream.push(5);
        stream.push(100);
        stream.push(0); // terminator
        // freq[5]=2048, freq[100]=2048 (sum = 4096 → no normalization)
        encode_uint7_prv(&mut stream, 2048);
        encode_uint7_prv(&mut stream, 2048);
        // 32 initial states. Start at (4096 << 12) | (j * 80) so lanes
        // differ in their s&0xFFF, producing different decoded symbols.
        for j in 0u32..32 {
            let s = (4096 << 12) | (j * 80);
            stream.extend_from_slice(&s.to_le_bytes());
        }

        let mut cur_simd: &[u8] = &stream;
        read_u8(&mut cur_simd).unwrap();
        let _ = read_uint7(&mut cur_simd).unwrap();
        let mut dst_simd = vec![0u8; len];
        decode_order_0_32state(&mut cur_simd, &mut dst_simd).unwrap();

        let mut cur_gen: &[u8] = &stream;
        read_u8(&mut cur_gen).unwrap();
        let _ = read_uint7(&mut cur_gen).unwrap();
        let mut dst_gen = vec![0u8; len];
        decode_order_0_generic(&mut cur_gen, &mut dst_gen, 32).unwrap();

        assert_eq!(dst_simd, dst_gen, "SIMD and scalar diverge for len={len}");
    }

    // ── Packed kernels vs the unpacked scalar decoders ──────────────────

    /// The alphabet encoding `read_alphabet` reads: symbols ascending, a
    /// symbol one above its predecessor followed by the count of further
    /// consecutive ones, a 0 terminator.
    fn encode_alphabet(syms: &[u8]) -> Vec<u8> {
        let mut out = Vec::new();
        let mut i = 0;
        while i < syms.len() {
            out.push(syms[i]);
            if i + 1 < syms.len() && syms[i + 1] == syms[i] + 1 {
                out.push(syms[i + 1]);
                let mut k = i + 2;
                while k < syms.len() && syms[k] == syms[k - 1] + 1 && k - i - 2 < 255 {
                    k += 1;
                }
                // The run marks syms[i+1] .. syms[k-2]; syms[k-1] is marked
                // when the loop comes round.
                out.push((k - i - 2) as u8);
                i = k;
            } else {
                i += 1;
            }
        }
        out.push(0);
        out
    }

    /// Frequencies for `n` symbols: all ≥ 1 and summing to exactly `total`
    /// (so the table packs), unless `exact` is false, when one is nudged
    /// off so the sum misses (and the decoders take the unpacked path).
    fn draw_frequencies(tc: &TestCase, n: usize, total: u32, exact: bool) -> Vec<u32> {
        let mut cuts: Vec<u32> = (0..n - 1)
            .map(|_| tc.draw(gs::integers::<u32>().min_value(1).max_value(total - 1)))
            .collect();
        cuts.push(0);
        cuts.push(total);
        cuts.sort_unstable();
        cuts.dedup();
        let mut f: Vec<u32> = cuts.windows(2).map(|w| w[1] - w[0]).collect();
        f.resize(n, 1); // dedup may have merged cuts; a surplus 1 breaks the sum
        if !exact {
            f[0] += 1;
        }
        f
    }

    /// A draw of the stream after the frequency table: `n` states (mostly
    /// above the renorm threshold, sometimes not) and renorm bytes, maybe
    /// cut short.
    fn draw_states_and_payload(tc: &TestCase, n: usize, len: usize) -> Vec<u8> {
        let low_states = tc.draw(gs::booleans());
        let mut out = Vec::new();
        for _ in 0..n {
            let x = if low_states {
                tc.draw(gs::integers::<u32>().max_value(1 << 20))
            } else {
                tc.draw(gs::integers::<u32>().min_value(1 << 15).max_value(1 << 31))
            };
            out.extend_from_slice(&x.to_le_bytes());
        }
        // ~1 byte per symbol covers a decode; a short cut covers truncation.
        let mut payload = tc.draw(gs::binary().min_size(len + 64).max_size(len + 128));
        if tc.draw(gs::booleans()) {
            payload.truncate(tc.draw(gs::integers::<usize>().max_value(len + 64)));
        }
        out.extend_from_slice(&payload);
        out
    }

    fn draw_alphabet(tc: &TestCase, with_zero: bool) -> Vec<u8> {
        let mut syms = tc.draw(gs::vecs(gs::integers::<u8>()).min_size(1).max_size(40));
        if with_zero {
            syms.push(0);
        }
        syms.sort_unstable();
        syms.dedup();
        syms
    }

    // r[verify cram.codec.simd_dispatch+3]
    // r[verify io.simd_portable]
    /// Order-0, 4 and 32 states: the packed decoders — at every SIMD level,
    /// and the scalar packed kernel at 32 states too — against the unpacked
    /// scalar decoder, on exact and inexact tables, high and low initial
    /// states, whole and truncated streams.
    #[hegel::test]
    fn packed_order0_matches_unpacked(tc: TestCase) {
        let n = if tc.draw(gs::booleans()) { 32 } else { 4 };
        let len = tc.draw(gs::integers::<usize>().max_value(600));
        let syms = draw_alphabet(&tc, false);
        let freqs = draw_frequencies(&tc, syms.len(), 4096, tc.draw(gs::booleans()));
        let mut stream = encode_alphabet(&syms);
        for &f in &freqs {
            encode_uint7_prv(&mut stream, f);
        }
        stream.extend(draw_states_and_payload(&tc, n, len));

        let mut oracle_src: &[u8] = &stream;
        let mut oracle_dst = vec![0u8; len];
        let oracle = decode_order_0_generic(&mut oracle_src, &mut oracle_dst, n);

        let check = |name: &str, got: Result<(), CramError>, dst: &[u8], src: &[u8]| {
            assert_eq!(got.is_ok(), oracle.is_ok(), "{name}: {got:?} vs {oracle:?}");
            if got.is_ok() {
                assert_eq!(dst, oracle_dst, "{name}");
                assert_eq!(src.len(), oracle_src.len(), "{name} consumed differently");
            }
        };
        for level in crate::simd_levels::levels() {
            let mut src: &[u8] = &stream;
            let mut dst = vec![0u8; len];
            let got = decode_order_0_at(level, &mut src, &mut dst, n);
            check(&format!("{level:?}"), got, &dst, src);
        }
        if n == 32 {
            let mut src: &[u8] = &stream;
            let mut dst = vec![0u8; len];
            let frequencies = read_frequencies_0(&mut src).unwrap();
            let mut table = [0u32; ROW];
            if pack_row(frequencies, ORDER_0_BITS, &mut table) {
                let mut states = read_state_array::<32>(&mut src).unwrap();
                let got = decode_o0_packed_scalar(&mut src, &mut dst, &table, &mut states);
                check("scalar packed", got, &dst, src);
            }
        }
    }

    /// An order-1 frequency table in `read_frequencies_1`'s uncompressed
    /// form: every alphabet symbol is a context, and each context's row over
    /// the alphabet has some zero entries (each followed by a run length of
    /// 0) and sums to `1 << bits` unless `exact` is false for that row.
    fn encode_order1_table(tc: &TestCase, syms: &[u8], bits: u32, exact_rows: bool) -> Vec<u8> {
        let mut out = vec![(bits as u8) << 4];
        out.extend(encode_alphabet(syms));
        let broken_row = if exact_rows {
            usize::MAX
        } else {
            tc.draw(gs::integers::<usize>().max_value(syms.len() - 1))
        };
        for row in 0..syms.len() {
            let live: Vec<bool> =
                syms.iter().map(|_| tc.draw(gs::integers::<u8>().max_value(3)) != 0).collect();
            let live_count = live.iter().filter(|&&l| l).count().max(1);
            let freqs = draw_frequencies(tc, live_count, 1 << bits, row != broken_row);
            let mut next = freqs.iter();
            for (k, &l) in live.iter().enumerate() {
                // At least one live symbol per row: the first, if none drew.
                if l || (k == 0 && !live.iter().any(|&l| l)) {
                    encode_uint7_prv(&mut out, *next.next().unwrap());
                } else {
                    out.push(0); // f = 0
                    out.push(0); // skip no further symbols
                }
            }
        }
        out
    }

    fn decode_order_1_oracle(src: &mut &[u8], dst: &mut [u8], n: usize) -> Result<(), CramError> {
        let mut t = Order1Tables::new();
        let (bits, alphabet) = read_frequencies_1(src, &mut t.frequencies)?;
        for (active, row) in alphabet.iter().zip(t.frequencies.iter_mut()) {
            if !active {
                row.fill(0);
            }
        }
        decode_order_1_unpacked(src, dst, n, &mut t, bits)
    }

    // r[verify cram.codec.simd_dispatch+3]
    // r[verify io.simd_portable]
    /// Order-1, 4 and 32 states, 10 and 12 bits: the packed decoders at every
    /// SIMD level (and the scalar packed kernel at 32 states) against the
    /// unpacked scalar decoder, with tables that pack, rows that do not, and
    /// alphabets without context 0.
    #[hegel::test]
    fn packed_order1_matches_unpacked(tc: TestCase) {
        let n = if tc.draw(gs::booleans()) { 32 } else { 4 };
        let bits = if tc.draw(gs::booleans()) { 12 } else { 10 };
        let len = tc.draw(gs::integers::<usize>().max_value(600));
        let syms = draw_alphabet(&tc, tc.draw(gs::integers::<u8>().max_value(7)) != 0);
        let mut stream =
            encode_order1_table(&tc, &syms, bits, tc.draw(gs::integers::<u8>().max_value(7)) != 0);
        stream.extend(draw_states_and_payload(&tc, n, len));

        let mut oracle_src: &[u8] = &stream;
        let mut oracle_dst = vec![0u8; len];
        let oracle = decode_order_1_oracle(&mut oracle_src, &mut oracle_dst, n);

        let check = |name: &str, got: Result<(), CramError>, dst: &[u8], src: &[u8]| {
            assert_eq!(got.is_ok(), oracle.is_ok(), "{name}: {got:?} vs {oracle:?}");
            if got.is_ok() {
                assert_eq!(dst, oracle_dst, "{name}");
                assert_eq!(src.len(), oracle_src.len(), "{name} consumed differently");
            }
        };
        for level in crate::simd_levels::levels() {
            let mut src: &[u8] = &stream;
            let mut dst = vec![0u8; len];
            let got = decode_order_1_at(level, &mut src, &mut dst, n, &mut Nx16Order1Buf::new());
            check(&format!("{level:?}"), got, &dst, src);
        }
        if n == 32 {
            let mut t = Order1Tables::new();
            let mut src: &[u8] = &stream;
            let (bits, alphabet) = read_frequencies_1(&mut src, &mut t.frequencies).unwrap();
            let packs = alphabet[0]
                && alphabet
                    .iter()
                    .zip(t.frequencies.iter())
                    .zip(t.packed.iter_mut())
                    .filter(|((a, _), _)| **a)
                    .all(|((_, f), row)| pack_row(f.iter().copied(), bits, row));
            if packs {
                let mut dst = vec![0u8; len];
                let mut states = read_state_array::<32>(&mut src).unwrap();
                let got = decode_o1_packed_scalar(&mut src, &mut dst, &t.packed, bits, &mut states);
                check("scalar packed", got, &dst, src);
            }
        }
    }

    // r[verify cram.codec.order1_table_reset]
    /// A reused buffer decodes a block exactly as a fresh one does, whatever
    /// an earlier block with a different alphabet left in it.
    #[hegel::test]
    fn order1_buffer_reuse_matches_fresh(tc: TestCase) {
        let n = 4;
        let blocks: Vec<(Vec<u8>, usize)> = (0..2)
            .map(|_| {
                let len = tc.draw(gs::integers::<usize>().max_value(300));
                let syms = draw_alphabet(&tc, true);
                let mut stream = encode_order1_table(&tc, &syms, 12, tc.draw(gs::booleans()));
                stream.extend(draw_states_and_payload(&tc, n, len));
                (stream, len)
            })
            .collect();
        let mut reused = Nx16Order1Buf::new();
        for (stream, len) in &blocks {
            let mut src: &[u8] = stream;
            let mut got = vec![0u8; *len];
            let r1 = decode_order_1_with_buf(&mut src, &mut got, n, &mut reused);
            let mut src: &[u8] = stream;
            let mut want = vec![0u8; *len];
            let r2 = decode_order_1_with_buf(&mut src, &mut want, n, &mut Nx16Order1Buf::new());
            assert_eq!(r1.is_ok(), r2.is_ok());
            if r1.is_ok() {
                assert_eq!(got, want);
            }
        }
    }
}
