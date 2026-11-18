//! The rANS 4×8 encoder's coding loops as SIMD kernels
//! (`r[cram_writer.rans4x8.simd]`): the four interleaved states in one
//! vector, each step renormalizing every lane by 0, 1 or 2 bytes, placed in
//! decode order by a shuffle the two lane masks choose.

use super::{Enc, Words, pair};
use fearless_simd::{Level, dispatch, i32x4, prelude::*, u8x16, u32x4, u64x2};
use fearless_simd_macros::simd;

// r[impl cram_writer.rans4x8.simd]

/// For each pair of 4-lane masks (`r1` the lanes emitting a byte, `r2 ⊆ r1`
/// those emitting two, as `r1 | r2 << 4`), the byte shuffle that moves
/// their bytes to the top of 16 bytes in decode order: lane by lane from
/// lane 0, bits 8–15 before bits 0–7.
#[allow(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    reason = "const table build: m < 256, lane < 4, at most 8 bytes"
)]
const EMIT: [[u8; 16]; 256] = {
    let mut t = [[0x80u8; 16]; 256];
    let mut m: usize = 0;
    while m < 256 {
        let (r1, r2) = (m & 0xF, m >> 4);
        let k = ((m as u32).count_ones()) as usize;
        let mut slot = 16 - k;
        let mut lane = 0;
        while lane < 4 {
            if r2 & (1 << lane) != 0 {
                t[m][slot] = (4 * lane + 1) as u8;
                slot += 1;
            }
            if r1 & (1 << lane) != 0 {
                t[m][slot] = (4 * lane) as u8;
                slot += 1;
            }
            lane += 1;
        }
        m += 1;
    }
    t
};

/// The bytes a mask pair of [`EMIT`] emits.
#[allow(clippy::cast_possible_truncation, reason = "at most 8")]
const EMITTED: [u8; 256] = {
    let mut t = [0u8; 256];
    let mut m: usize = 0;
    while m < 256 {
        t[m] = (m as u32).count_ones() as u8;
        m += 1;
    }
    t
};

/// Four symbols' reciprocals and their start, frequency and shift.
#[inline(always)]
fn split<S: Simd>(simd: S, e: &[Enc; 4]) -> (u32x4<S>, u32x4<S>) {
    let (lo, hi) = e.split_at(2);
    let a: u32x4<S> = u64x2::from_slice(simd, lo).bitcast();
    let b: u32x4<S> = u64x2::from_slice(simd, hi).bitcast();
    (a.unzip_low(b), a.unzip_high(b))
}

/// The parent's `put`, lane by lane, the lanes' bytes stored in lane order
/// below the ones before.
#[inline(always)]
#[allow(clippy::indexing_slicing, reason = "a mask pair is < 256; at ≥ 16 (see `Words`)")]
#[allow(clippy::cast_possible_truncation, reason = "a 4-lane mask is < 16")]
#[allow(clippy::arithmetic_side_effects, reason = "vector lanes wrap; as in `put`")]
fn step4<S: Simd>(simd: S, x: u32x4<S>, e: &[Enc; 4], words: &mut Words<'_>) -> u32x4<S> {
    let (m, info) = split(simd, e);
    let f = (info >> 12) & 0x1FFF;
    // x ≥ f·2^19 as x > f·2^19 − 1, both below 2^31: a signed compare.
    let below: i32x4<S> = ((f << 19) - 1).bitcast();
    let x8 = x >> 8;
    let r1 = x.bitcast::<i32x4<S>>().simd_gt(below);
    let r2 = x8.bitcast::<i32x4<S>>().simd_gt(below);
    let mask = ((r1.to_bitmask() & 0xF) | (r2.to_bitmask() & 0xF) << 4) as usize;
    let at = words.at;
    debug_assert!(at >= 16, "a word buffer starts with 16 spare bytes");
    let low: u8x16<S> = x.bitcast();
    low.swizzle_dyn(u8x16::from_slice(simd, &EMIT[mask])).store_slice(&mut words.buf[at - 16..at]);
    words.at = at - usize::from(EMITTED[mask]);
    let s = info >> 25;
    let x16 = x >> 16;
    let (x, q) = if cfg!(target_arch = "aarch64") {
        let q = r2.select(quotient(x16, m, s), r1.select(quotient(x8, m, s), quotient(x, m, s)));
        (r2.select(x16, r1.select(x8, x)), q)
    } else {
        // Both shifted states are ready with the masks: two selects in a row.
        let x = r2.select(x16, r1.select(x8, x));
        (x, quotient(x, m, s))
    };
    x + (info & 0xFFF) + q * (u32x4::splat(simd, 1 << 12) - f)
}

/// `⌊x / f⌋` lane by lane, with the reciprocal `m` and shift `s`; `x < 2^31`.
#[inline(always)]
#[allow(clippy::arithmetic_side_effects, reason = "vector lanes wrap")]
fn quotient<S: Simd>(x: u32x4<S>, m: u32x4<S>, s: u32x4<S>) -> u32x4<S> {
    let (lo, hi) = (x + x).widening_mul(m);
    lo.bitcast::<u32x4<S>>().unzip_high(hi.bitcast::<u32x4<S>>()) >> s
}

/// Marks every byte value of `data` in `present`. A vector of 16 bytes is
/// tested against the values marked so far, as a 16 × 16 bit table looked
/// up by each byte's low nibble (two shuffles) and masked by its high
/// nibble (two more); only a vector holding a new value is marked byte by
/// byte. A store per byte took a cycle per byte on Zen 2.
pub(super) fn mark_present(level: Level, data: &[u8], present: &mut [bool; 256]) {
    dispatch!(level, simd => mark(simd, data, present));
}

/// For a high nibble `h`, its bit in the half of the table it is in: bit
/// `h` of the first half for `h < 8`, else 0 (`HIGH_BITS` the other half).
const LOW_BITS: [u8; 16] = [1, 2, 4, 8, 16, 32, 64, 128, 0, 0, 0, 0, 0, 0, 0, 0];
const HIGH_BITS: [u8; 16] = [0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 4, 8, 16, 32, 64, 128];

#[simd]
#[allow(clippy::indexing_slicing, reason = "u8s index a 256-entry table; nibbles index 16")]
#[allow(clippy::arithmetic_side_effects, reason = "a vector shift of u8 lanes")]
fn mark<S: Simd>(simd: S, data: &[u8], present: &mut [bool; 256]) {
    // `table[h / 8][l]` holds bit `h % 8` for value `h·16 + l`.
    let mut table = [[0u8; 16]; 2];
    for (b, _) in (0u8..=255).zip(present.iter()).filter(|(_, p)| **p) {
        table[usize::from(b >> 7)][usize::from(b & 0xF)] |= 1 << ((b >> 4) & 7);
    }
    let low_bits = u8x16::from_slice(simd, &LOW_BITS);
    let high_bits = u8x16::from_slice(simd, &HIGH_BITS);
    let nibble = u8x16::splat(simd, 0xF);
    let zero = u8x16::splat(simd, 0);
    let mut low_half = u8x16::from_slice(simd, &table[0]);
    let mut high_half = u8x16::from_slice(simd, &table[1]);
    let (chunks, tail) = data.as_chunks::<16>();
    for chunk in chunks {
        let v = u8x16::from_slice(simd, chunk);
        let (lo, hi) = (v & nibble, v >> 4);
        let hits = (low_half.swizzle_dyn(lo) & low_bits.swizzle_dyn(hi))
            | (high_half.swizzle_dyn(lo) & high_bits.swizzle_dyn(hi));
        if hits.simd_eq(zero).any_true() {
            for &b in chunk {
                present[usize::from(b)] = true;
                table[usize::from(b >> 7)][usize::from(b & 0xF)] |= 1 << ((b >> 4) & 7);
            }
            low_half = u8x16::from_slice(simd, &table[0]);
            high_half = u8x16::from_slice(simd, &table[1]);
        }
    }
    for &b in tail {
        present[usize::from(b)] = true;
    }
}

/// Order 0 over whole rounds of four bytes, last first.
pub(super) fn code_order_0(
    level: Level,
    rounds: &[[u8; 4]],
    enc: &[Enc; 256],
    states: &mut [u32; 4],
    words: &mut Words<'_>,
) {
    dispatch!(level, simd => code_o0(simd, rounds, enc, states, words));
}

#[simd]
#[allow(clippy::indexing_slicing, reason = "u8s index a 256-entry table")]
fn code_o0<S: Simd>(
    simd: S,
    rounds: &[[u8; 4]],
    enc: &[Enc; 256],
    states: &mut [u32; 4],
    words: &mut Words<'_>,
) {
    // A local copy, so that the write position stays in a register.
    let mut local = Words { buf: &mut *words.buf, at: words.at };
    let mut x = u32x4::from_slice(simd, states);
    for &[a, b, c, d] in rounds.iter().rev() {
        let e =
            [enc[usize::from(a)], enc[usize::from(b)], enc[usize::from(c)], enc[usize::from(d)]];
        x = step4(simd, x, &e, &mut local);
    }
    x.store_slice(states);
    words.at = local.at;
}

/// Order 1 over the four quarters' positions from 1 on, last first; each
/// slice's windows hold a position's byte and the one before.
pub(super) fn code_order_1(
    level: Level,
    quarters: [&[u8]; 4],
    enc: &[Enc; 1 << 16],
    states: &mut [u32; 4],
    words: &mut Words<'_>,
) {
    dispatch!(level, simd => code_o1(simd, quarters, enc, states, words));
}

/// Positions an order-1 block covers: each quarter's bytes for them, and the
/// one before, are copied out first: read in place, the quarters' bytes
/// trail the bytes written at a near-constant distance, which can alias
/// them modulo 4 KiB for a whole block (Zen 2: order 1 at 10 to 13 cycles
/// a byte by where the buffers landed, 8 to 10 with the copy).
const BLOCK: usize = 16;

#[simd]
#[allow(
    clippy::indexing_slicing,
    reason = "u16s index a 65536-entry table; windows of 2; k + 1 ≤ BLOCK"
)]
#[allow(clippy::arithmetic_side_effects, reason = "end > BLOCK ≥ 1 in the block loop")]
fn code_o1<S: Simd>(
    simd: S,
    quarters: [&[u8]; 4],
    enc: &[Enc; 1 << 16],
    states: &mut [u32; 4],
    words: &mut Words<'_>,
) {
    let mut local = Words { buf: &mut *words.buf, at: words.at };
    let mut x = u32x4::from_slice(simd, states);
    let [q0, q1, q2, q3] = quarters;
    // Whole blocks of positions `end - BLOCK .. end`, all from 1 on.
    let mut end = q0.len();
    while end > BLOCK {
        let start = end - BLOCK;
        let mut window = [[0u8; BLOCK + 1]; 4];
        for (w, q) in window.iter_mut().zip(quarters) {
            w.copy_from_slice(&q[start - 1..end]);
        }
        let [a, b, c, d] = &window;
        for k in (0..BLOCK).rev() {
            let e = [
                enc[pair(a[k], a[k + 1])],
                enc[pair(b[k], b[k + 1])],
                enc[pair(c[k], c[k + 1])],
                enc[pair(d[k], d[k + 1])],
            ];
            x = step4(simd, x, &e, &mut local);
        }
        end = start;
    }
    let [q0, q1, q2, q3] = [&q0[..end], &q1[..end], &q2[..end], &q3[..end]];
    let positions = q0.windows(2).zip(q1.windows(2)).zip(q2.windows(2)).zip(q3.windows(2));
    for (((a, b), c), d) in positions.rev() {
        let e = [
            enc[pair(a[0], a[1])],
            enc[pair(b[0], b[1])],
            enc[pair(c[0], c[1])],
            enc[pair(d[0], d[1])],
        ];
        x = step4(simd, x, &e, &mut local);
    }
    x.store_slice(states);
    words.at = local.at;
}
