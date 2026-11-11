//! Encode and decode BAM 4-bit packed sequences. Each direction is one
//! `fearless_simd` kernel, run at the best level the CPU has; every level,
//! the scalar fallback included, produces identical output.
//!
//! - **Decode** (`decode_nibbles`): each nibble through a 16-entry table with
//!   a byte shuffle, the high and low results interleaved back into read
//!   order. A table and its 256-entry pair table (two nibbles per byte, for
//!   the scalar tail) come together: `DECODE_BASE`/`DECODE_PAIR` (ASCII,
//!   every IUPAC code) and `IDENTITY`/`IDENTITY_PAIR` (`IupacBase` codes).
//!   `NibbleOrder` says which nibble comes first: BAM's high one, or the low
//!   one of CRAM's rANS PACK transform (`cram::rans_nx16`).
//! - **Encode**: `IupacBase` codes are already nibbles and need no lookup:
//!   a plain pair loop, which LLVM vectorises, beat a byte-shuffle kernel
//!   1.4x on NEON. ASCII text (`encode_seq_into`, the SAM encoder) is one
//!   kernel: even and odd bases unzipped apart, each byte's low nibble
//!   shuffled through the `ENCODE_BASE` row its high nibble picks.
//! - The `IupacBase` kernels fill their output through
//!   `IupacBase::as_codes_mut`, whose contract is that only codes (< 16) are
//!   left behind; the identity tables hold nothing else.

use fearless_simd::{Level, dispatch, prelude::*, u8x16};
use fearless_simd_macros::simd;
use seqair_types::IupacBase;

// r[impl seq.decode_scalar]
// r[impl seq.decode_pair_table]
// r[impl seq.decode_simd+2]
// r[impl seq.decode_dispatch+2]
// r[impl seq.encode_scalar+2]
// r[impl seq.simd_scalar_equivalence]
// r[impl io.platform_optimizations]

/// BAM 4-bit encoding → ASCII base lookup.
static DECODE_BASE: &[u8; 16] = &IupacBase::ASCII;

/// Pre-computed pair table: `DECODE_PAIR[byte]` → `[high_nibble_base, low_nibble_base]`.
#[allow(clippy::indexing_slicing, reason = "i < 256 = table.len(), nibbles < 16 = BASE.len()")]
static DECODE_PAIR: [[u8; 2]; 256] = {
    const BASE: [u8; 16] = IupacBase::ASCII;
    let mut table = [[0u8; 2]; 256];
    let mut i = 0;
    while i < 256 {
        table[i] = [BASE[i >> 4], BASE[i & 0xF]];
        i += 1;
    }
    table
};

/// ASCII → 4-bit encoding table: htslib's `seq_nt16_table`, so SAM text
/// encodes exactly as htslib encodes it. Unknown characters map to 15 (N).
static ENCODE_BASE: [u8; 256] = IupacBase::CODE_OF_ASCII;

/// Decode a 4-bit packed BAM sequence to ASCII bases (scalar fallback).
pub fn decode_seq_scalar(encoded: &[u8], len: usize) -> Vec<u8> {
    let required = len.div_ceil(2);
    if encoded.len() < required {
        return vec![b'N'; len];
    }
    let full_bytes = len / 2;
    let mut result = vec![0u8; len];

    #[allow(clippy::indexing_slicing, reason = "bounds ensured by zip + as_chunks")]
    for (chunk, &byte) in result.as_chunks_mut::<2>().0.iter_mut().zip(&encoded[..full_bytes]) {
        let pair = DECODE_PAIR[byte as usize];
        chunk[0] = pair[0];
        chunk[1] = pair[1];
    }

    if len % 2 == 1
        && let Some(byte) = encoded.get(full_bytes)
        && let Some(slot) = result.get_mut(len.checked_sub(1).expect("len - 1 underflow"))
    {
        // byte is u8 so byte < 256 = DECODE_PAIR.len(); [0] is always valid on [u8; 2]
        #[allow(clippy::indexing_slicing, reason = "byte < 256 = DECODE_PAIR.len()")]
        {
            *slot = DECODE_PAIR[*byte as usize][0];
        }
    }

    result
}

/// Decode a 4-bit packed BAM sequence using the best available implementation.
pub fn decode_seq(encoded: &[u8], len: usize) -> Vec<u8> {
    let mut out = vec![0u8; len];
    decode_iupac_into(encoded, &mut out);
    out
}

/// Decode `out.len()` bases of a 4-bit packed BAM sequence into `out` as
/// ASCII, every code kept (`=ACMGRSVTWYHKDBN`); all `N` when `encoded` is too
/// short.
pub(crate) fn decode_iupac_into(encoded: &[u8], out: &mut [u8]) {
    if !decode_nibbles(
        Level::new(),
        NibbleOrder::HighFirst,
        DECODE_BASE,
        &DECODE_PAIR,
        encoded,
        out,
    ) {
        out.fill(b'N');
    }
}

/// Which nibble of a packed byte comes first: BAM's high one, or the low
/// one of CRAM's rANS Nx16 PACK transform.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum NibbleOrder {
    HighFirst,
    LowFirst,
}

/// Decode `out.len()` nibbles of `encoded` through `lut`, the 16-entry nibble
/// table, and `pairs`, the same table two nibbles at a time (in `order`) for
/// the scalar tail. Returns `false`, having written nothing, when `encoded`
/// is too short.
// r[impl io.simd_portable]
pub(crate) fn decode_nibbles(
    level: Level,
    order: NibbleOrder,
    lut: &[u8; 16],
    pairs: &[[u8; 2]; 256],
    encoded: &[u8],
    out: &mut [u8],
) -> bool {
    let Some(packed) = encoded.get(..out.len().div_ceil(2)) else { return false };
    dispatch!(level, simd => decode_nibbles_simd(simd, order, lut, pairs, packed, out));
    true
}

/// One native-width vector of packed bytes per iteration — two vectors of
/// bases out. Each nibble is looked up with a byte shuffle against `lut`
/// (`pshufb` / `tbl`), the high-nibble and low-nibble results are
/// interleaved back into read order, and whatever does not fill a vector goes
/// through `pairs`. `packed.len()` is `out.len().div_ceil(2)`.
#[simd]
#[allow(
    clippy::chunks_exact_to_as_chunks,
    reason = "`S::u8s::LEN` depends on the generic `S`, so it cannot be a const argument"
)]
#[allow(clippy::indexing_slicing, reason = "a `u8` index is < 256 = pairs.len()")]
#[allow(
    clippy::arithmetic_side_effects,
    reason = "a lanewise shift by 4 of a `u8` cannot overflow"
)]
fn decode_nibbles_simd<S: Simd>(
    simd: S,
    order: NibbleOrder,
    lut: &[u8; 16],
    pairs: &[[u8; 2]; 256],
    packed: &[u8],
    out: &mut [u8],
) {
    let n = S::u8s::LEN;
    // `swizzle_dyn_within_blocks` shuffles within each 128-bit block, so the
    // table is repeated into every block of a wider vector.
    let table = S::u8s::block_splat(u8x16::from_slice(simd, lut));
    let (out_pairs, odd) = out.as_chunks_mut::<2>();
    let (whole, last) = packed.split_at(out_pairs.len().min(packed.len()));

    let mut vec_in = whole.chunks_exact(n);
    for (p, o) in (&mut vec_in).zip(out_pairs.chunks_exact_mut(n)) {
        decode_vector(simd, order, table, p, o);
    }

    // A ragged tail is one more vector overlapping the last whole one: the
    // output is a pure function of the input, and the input is never written,
    // so recomputing the overlap is harmless. A 150 bp read leaves 11 of its
    // 75 packed bytes to the tail, which a byte loop would pay for.
    let tail = vec_in.remainder();
    let overlap = whole.len().checked_sub(n).filter(|_| !tail.is_empty());
    match overlap.and_then(|t| Some((whole.get(t..)?, out_pairs.get_mut(t..)?))) {
        Some((p, o)) => decode_vector(simd, order, table, p, o),
        // Shorter than one vector: every byte is tail.
        None => {
            for (o, &b) in out_pairs.iter_mut().zip(tail) {
                *o = pairs[usize::from(b)];
            }
        }
    }
    if let [o] = odd
        && let Some(&b) = last.first()
    {
        *o = pairs[usize::from(b)][0];
    }
}

/// One vector of packed bytes through `table` into two vectors of bases, in
/// read order. A function, not a closure: a closure does not get the
/// enclosing `#[simd]` level's target features, so on x86 every vector op in
/// it becomes an out-of-line call.
#[inline(always)]
#[allow(
    clippy::arithmetic_side_effects,
    reason = "a lanewise shift by 4 of a `u8` cannot overflow"
)]
fn decode_vector<S: Simd>(
    simd: S,
    order: NibbleOrder,
    table: S::u8s,
    packed: &[u8],
    out: &mut [[u8; 2]],
) {
    let p = S::u8s::from_slice(simd, packed);
    let hi = table.swizzle_dyn_within_blocks(p >> 4);
    let lo = table.swizzle_dyn_within_blocks(p & 0x0F);
    let (first, second) = match order {
        NibbleOrder::HighFirst => hi.interleave(lo),
        NibbleOrder::LowFirst => lo.interleave(hi),
    };
    let (o_first, o_second) = out.as_flattened_mut().split_at_mut(S::u8s::LEN);
    first.store_slice(o_first);
    second.store_slice(o_second);
}

/// The rows of [`ENCODE_BASE`] the vector encoder looks up: ASCII `0x30` to
/// `0x7F`, one 16-entry row per high nibble. Every other byte is `N`, which
/// the assertion below checks against the table at compile time.
const ENCODE_ROW_FIRST: u8 = 0x3;
const ENCODE_ROWS: [[u8; 16]; 5] = {
    let mut rows = [[0u8; 16]; 5];
    let mut b = 0;
    while b < 256 {
        let row = b >> 4;
        let code = IupacBase::CODE_OF_ASCII[b];
        if row >= ENCODE_ROW_FIRST as usize && row < ENCODE_ROW_FIRST as usize + 5 {
            rows[row - ENCODE_ROW_FIRST as usize][b & 0xF] = code;
        } else {
            assert!(code == 15, "ENCODE_BASE outside 0x30..0x80 is N");
        }
        b += 1;
    }
    rows
};

// r[impl seq.encode_simd]
/// Encode ASCII bases to 4-bit packed BAM format, appending
/// `bases.len().div_ceil(2)` bytes to `buf`; an odd length leaves the last
/// low nibble 0.
pub fn encode_seq_into(bases: &[u8], buf: &mut Vec<u8>) {
    encode_ascii_at(Level::new(), bases, buf);
}

// r[impl io.simd_portable]
fn encode_ascii_at(level: Level, bases: &[u8], buf: &mut Vec<u8>) {
    let start = buf.len();
    buf.resize(start.saturating_add(bases.len().div_ceil(2)), 0);
    let out = buf.get_mut(start..).unwrap_or_default();
    dispatch!(level, simd => encode_ascii_simd(simd, bases, out));
}

/// Two native-width vectors of ASCII in, one of packed bytes out per
/// iteration: the even and odd bases are unzipped apart, each coded by
/// [`encode_codes`], and the even code shifted into the high nibble. A ragged
/// tail is one more step overlapping the last whole one, as in
/// [`decode_nibbles_simd`]; inputs shorter than one step go through
/// [`ENCODE_BASE`] a pair at a time.
#[simd]
#[allow(
    clippy::chunks_exact_to_as_chunks,
    reason = "`S::u8s::LEN` depends on the generic `S`, so it cannot be a const argument"
)]
#[allow(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "a u8 indexes a 256-entry table; a 4-bit code shifted by 4 fits a u8"
)]
fn encode_ascii_simd<S: Simd>(simd: S, bases: &[u8], out: &mut [u8]) {
    let n = S::u8s::LEN;
    let rows = ENCODE_ROWS.map(|row| S::u8s::block_splat(u8x16::from_slice(simd, &row)));
    let (pairs, odd) = bases.as_chunks::<2>();
    let pairs = pairs.as_flattened();
    let (whole, last) = out.split_at_mut(pairs.len() / 2);

    let mut vec_in = pairs.chunks_exact(2 * n);
    for (b, o) in (&mut vec_in).zip(whole.chunks_exact_mut(n)) {
        encode_vector(simd, &rows, b, o);
    }
    let tail = vec_in.remainder();
    let overlap = whole.len().checked_sub(n).filter(|_| !tail.is_empty());
    match overlap.and_then(|t| Some((pairs.get(2 * t..)?, whole.get_mut(t..)?))) {
        Some((b, o)) => encode_vector(simd, &rows, b, o),
        None => {
            for (o, &[hi, lo]) in whole.iter_mut().zip(pairs.as_chunks::<2>().0) {
                *o = (ENCODE_BASE[usize::from(hi)] << 4) | ENCODE_BASE[usize::from(lo)];
            }
        }
    }
    if let ([b], [o]) = (odd, last) {
        *o = ENCODE_BASE[usize::from(*b)] << 4;
    }
}

/// `2 × LEN` ASCII bases into `LEN` packed bytes. A function, not a closure,
/// for the reason given at [`decode_vector`].
#[inline(always)]
#[allow(clippy::arithmetic_side_effects, reason = "a lanewise shift by 4 of a code < 16 fits a u8")]
fn encode_vector<S: Simd>(simd: S, rows: &[S::u8s; 5], bases: &[u8], out: &mut [u8]) {
    let (a, b) = bases.split_at(S::u8s::LEN);
    let (a, b) = (S::u8s::from_slice(simd, a), S::u8s::from_slice(simd, b));
    let even = encode_codes(simd, rows, a.unzip_low(b));
    let odd = encode_codes(simd, rows, a.unzip_high(b));
    ((even << 4) | odd).store_slice(out);
}

/// Each byte's [`ENCODE_BASE`] code: its low nibble shuffled through the row
/// its high nibble picks, `N` where no row is picked.
#[inline(always)]
#[allow(clippy::arithmetic_side_effects, reason = "a lanewise shift of a u8 cannot overflow")]
fn encode_codes<S: Simd>(simd: S, rows: &[S::u8s; 5], v: S::u8s) -> S::u8s {
    let high = v >> 4;
    let low = v & 0x0F;
    let mut code = S::u8s::splat(simd, IupacBase::N.code());
    for (i, row) in (ENCODE_ROW_FIRST..).zip(rows) {
        code = high.simd_eq(i).select(row.swizzle_dyn_within_blocks(low), code);
    }
    code
}

/// The identity table: `IupacBase` code `i` is nibble `i`.
static IDENTITY: &[u8; 16] = &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];

#[allow(clippy::indexing_slicing, reason = "i < 256, nibbles < 16")]
#[allow(clippy::cast_possible_truncation, reason = "nibbles < 16")]
static IDENTITY_PAIR: [[u8; 2]; 256] = {
    let mut table = [[0u8; 2]; 256];
    let mut i = 0;
    while i < 256 {
        table[i] = [(i >> 4) as u8, (i & 0xF) as u8];
        i += 1;
    }
    table
};

// r[impl seq.iupac_bulk+3]
/// Encode `IupacBase`s to 4-bit packed BAM format, appending
/// `bases.len().div_ceil(2)` bytes to `buf`; an odd length leaves the last
/// low nibble 0. The codes are already nibbles, so this is a pair loop LLVM
/// vectorises (`ld2` on NEON), measured faster than the shuffle kernel.
#[allow(clippy::arithmetic_side_effects, reason = "a 4-bit code shifted by 4 fits a u8")]
pub fn encode_iupac_bases_into(bases: &[IupacBase], buf: &mut Vec<u8>) {
    let (pairs, odd) = bytemuck::cast_slice::<IupacBase, u8>(bases).as_chunks::<2>();
    buf.reserve(bases.len().div_ceil(2));
    buf.extend(pairs.iter().map(|&[hi, lo]| (hi << 4) | lo));
    buf.extend(odd.iter().map(|last| last << 4));
}

/// [`encode_iupac_bases_into`] over `out`, which is `bases.len().div_ceil(2)`
/// bytes.
// r[impl seq.iupac_bulk+3]
#[allow(clippy::arithmetic_side_effects, reason = "a 4-bit code shifted by 4 fits a u8")]
pub(crate) fn encode_iupac_bases_to(bases: &[IupacBase], out: &mut [u8]) {
    let (pairs, odd) = bytemuck::cast_slice::<IupacBase, u8>(bases).as_chunks::<2>();
    for (o, &[hi, lo]) in out.iter_mut().zip(pairs) {
        *o = (hi << 4) | lo;
    }
    if let ([last], Some(o)) = (odd, out.get_mut(pairs.len())) {
        *o = last << 4;
    }
}

// r[impl seq.iupac_bulk+3]
// r[impl bam.record.seq_4bit]
/// Decode `out.len()` bases of a 4-bit packed BAM sequence into `out`, every
/// code kept; all `N` when `encoded` is too short.
pub fn decode_iupac_bases_into(encoded: &[u8], out: &mut [IupacBase]) {
    // SAFETY: `decode_nibbles` writes only entries of its tables, and the
    // identity tables hold only codes 0..16.
    let codes = unsafe { IupacBase::as_codes_mut(out) };
    if !decode_nibbles(
        Level::new(),
        NibbleOrder::HighFirst,
        IDENTITY,
        &IDENTITY_PAIR,
        encoded,
        codes,
    ) {
        codes.fill(IupacBase::N.code());
    }
}

/// Encode ASCII bases to 4-bit packed BAM format.
pub fn encode_seq(bases: &[u8]) -> Vec<u8> {
    let n_bytes = bases.len().div_ceil(2);
    let mut encoded = vec![0u8; n_bytes];

    #[allow(clippy::indexing_slicing, reason = "bounds ensured by step_by loop")]
    for j in (0..bases.len()).step_by(2) {
        debug_assert!(j < bases.len(), "step_by loop: j={j} out of bounds len={}", bases.len());
        debug_assert!(
            j / 2 < n_bytes,
            "encoded index j/2={} out of bounds n_bytes={n_bytes}",
            j / 2
        );
        let hi = ENCODE_BASE[bases[j] as usize];
        let lo = if j.saturating_add(1) < bases.len() {
            ENCODE_BASE[bases[j.saturating_add(1)] as usize]
        } else {
            0
        };
        encoded[j / 2] = (hi << 4) | lo;
    }

    encoded
}

#[cfg(test)]
#[allow(
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    reason = "test code with known small values"
)]
mod tests {
    use super::*;
    use hegel::prelude::*;

    /// Nibble `i` of `packed` through `lut`, one nibble at a time: shares
    /// neither the pair table nor the vector loop with the kernel.
    fn nibble_oracle(lut: &[u8; 16], packed: &[u8], len: usize) -> Vec<u8> {
        (0..len)
            .map(|i| {
                let byte = packed[i / 2];
                let nibble = if i % 2 == 0 { byte >> 4 } else { byte & 0x0F };
                lut[usize::from(nibble)]
            })
            .collect()
    }

    // r[verify seq.decode_simd+2]
    // r[verify seq.decode_dispatch+2]
    // r[verify seq.simd_scalar_equivalence]
    // r[verify seq.iupac_bulk+3]
    // r[verify io.simd_portable]
    #[hegel::test]
    fn every_level_matches_nibble_oracle(tc: TestCase) {
        let packed = tc.draw(gs::binary().max_size(300));
        // Up to one nibble fewer than `packed` holds, so an odd `len` leaves
        // the low nibble of the last byte unread.
        let len = tc.draw(gs::integers::<usize>().max_value(packed.len() * 2));
        for (lut, pairs) in [(DECODE_BASE, &DECODE_PAIR), (IDENTITY, &IDENTITY_PAIR)] {
            let expected = nibble_oracle(lut, &packed, len);
            for level in crate::simd_levels::levels() {
                let mut out = vec![0u8; len];
                assert!(decode_nibbles(
                    level,
                    NibbleOrder::HighFirst,
                    lut,
                    pairs,
                    &packed,
                    &mut out
                ));
                assert_eq!(out, expected, "{level:?}, len={len}");
            }
        }
    }

    // r[verify io.simd_portable]
    /// Low nibble first (CRAM's PACK), through an arbitrary table.
    #[hegel::test]
    fn every_level_matches_nibble_oracle_low_first(tc: TestCase) {
        let packed = tc.draw(gs::binary().max_size(300));
        let len = tc.draw(gs::integers::<usize>().max_value(packed.len() * 2));
        let lut: [u8; 16] = std::array::from_fn(|_| tc.draw(gs::integers::<u8>()));
        let pairs: [[u8; 2]; 256] = std::array::from_fn(|b| [lut[b & 0x0F], lut[b >> 4]]);
        let expected: Vec<u8> = (0..len)
            .map(|i| {
                let byte = packed[i / 2];
                lut[usize::from(if i % 2 == 0 { byte & 0x0F } else { byte >> 4 })]
            })
            .collect();
        for level in crate::simd_levels::levels() {
            let mut out = vec![0u8; len];
            assert!(decode_nibbles(level, NibbleOrder::LowFirst, &lut, &pairs, &packed, &mut out));
            assert_eq!(out, expected, "{level:?}, len={len}");
        }
    }

    // r[verify seq.encode_simd]
    // r[verify io.simd_portable]
    #[hegel::test]
    fn every_level_encodes_as_the_table(tc: TestCase) {
        // Mostly bases, so long runs exercise the vector path, plus any byte.
        let bases = tc.draw(
            gs::vecs(gs::one_of(vec![
                gs::sampled_from(b"ACGTNacgtn=RYMKSWBDHVUu0123".to_vec()).boxed(),
                gs::integers::<u8>().boxed(),
            ]))
            .max_size(300),
        );
        let prefix = tc.draw(gs::binary().max_size(3));
        // One base at a time through the 256-entry table: shares neither the
        // pair loop nor the vector kernel.
        let mut expected = prefix.clone();
        expected.extend(bases.chunks(2).map(|pair| {
            let code = |i: usize| pair.get(i).map_or(0, |&b| ENCODE_BASE[usize::from(b)]);
            (code(0) << 4) | code(1)
        }));
        for level in crate::simd_levels::levels() {
            let mut out = prefix.clone();
            encode_ascii_at(level, &bases, &mut out);
            assert_eq!(out, expected, "{level:?}, len={}", bases.len());
        }
    }

    // r[verify seq.iupac_bulk+3]
    // r[verify io.simd_portable]
    /// Codes pack as the spec's high-nibble-first pairs, from any offset
    /// into the input, appending after what `buf` holds, and decode back at
    /// the best level (`every_level_matches_nibble_oracle` covers the rest).
    #[hegel::test]
    fn encode_iupac_matches_oracle(tc: TestCase) {
        let len = tc.draw(gs::integers::<usize>().max_value(300));
        let offset = tc.draw(gs::integers::<usize>().max_value(63));
        let buf: Vec<IupacBase> = (0..offset + len)
            .map(|_| IupacBase::from_nibble(tc.draw(gs::integers::<u8>())))
            .collect();
        let seq = &buf[offset..];
        let expected: Vec<u8> =
            seq.chunks(2).map(|p| (p[0].code() << 4) | p.get(1).map_or(0, |b| b.code())).collect();
        let mut out = vec![0xAB];
        encode_iupac_bases_into(seq, &mut out);
        assert_eq!(out[0], 0xAB);
        let out = out.split_off(1);
        assert_eq!(out, expected);
        let mut to = vec![0xFF; expected.len()];
        encode_iupac_bases_to(seq, &mut to);
        assert_eq!(to, expected);
        let mut decoded = vec![IupacBase::R; len];
        decode_iupac_bases_into(&out, &mut decoded);
        assert_eq!(decoded, seq);
    }

    // r[verify seq.iupac_bulk+3]
    #[test]
    fn short_input_decodes_to_iupac_n() {
        let mut out = [IupacBase::A; 3];
        decode_iupac_bases_into(&[0x12], &mut out);
        assert_eq!(out, [IupacBase::N; 3]);
    }

    #[test]
    fn short_input_decodes_to_n() {
        assert_eq!(decode_seq(&[0x12], 3), b"NNN");
    }
}
