//! Encode and decode BAM 4-bit packed sequences. [`decode_seq`] runs one
//! `fearless_simd` kernel at the best SIMD level the CPU has; every level,
//! the scalar fallback included, produces identical output.

use fearless_simd::{Level, dispatch, prelude::*, u8x16};
use fearless_simd_macros::simd;

// r[impl seq.decode_scalar]
// r[impl seq.decode_pair_table]
// r[impl seq.decode_simd+2]
// r[impl seq.decode_dispatch+2]
// r[impl seq.encode_scalar]
// r[impl seq.simd_scalar_equivalence]
// r[impl io.platform_optimizations]

/// BAM 4-bit encoding → ASCII base lookup.
static DECODE_BASE: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

/// Pre-computed pair table: `DECODE_PAIR[byte]` → `[high_nibble_base, low_nibble_base]`.
#[allow(clippy::indexing_slicing, reason = "i < 256 = table.len(), nibbles < 16 = BASE.len()")]
static DECODE_PAIR: [[u8; 2]; 256] = {
    const BASE: [u8; 16] = *b"=ACMGRSVTWYHKDBN";
    let mut table = [[0u8; 2]; 256];
    let mut i = 0;
    while i < 256 {
        table[i] = [BASE[i >> 4], BASE[i & 0xF]];
        i += 1;
    }
    table
};

/// ASCII → 4-bit encoding table. Unknown characters map to 15 (N).
#[allow(
    clippy::indexing_slicing,
    reason = "all indices are ASCII byte literals < 256 = table.len()"
)]
static ENCODE_BASE: [u8; 256] = {
    let mut table = [15u8; 256];
    table[b'=' as usize] = 0;
    table[b'A' as usize] = 1;
    table[b'a' as usize] = 1;
    table[b'C' as usize] = 2;
    table[b'c' as usize] = 2;
    table[b'M' as usize] = 3;
    table[b'G' as usize] = 4;
    table[b'g' as usize] = 4;
    table[b'R' as usize] = 5;
    table[b'S' as usize] = 6;
    table[b'V' as usize] = 7;
    table[b'T' as usize] = 8;
    table[b't' as usize] = 8;
    table[b'W' as usize] = 9;
    table[b'Y' as usize] = 10;
    table[b'H' as usize] = 11;
    table[b'K' as usize] = 12;
    table[b'D' as usize] = 13;
    table[b'B' as usize] = 14;
    table[b'N' as usize] = 15;
    table[b'n' as usize] = 15;
    table
};

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
    if !decode_nibbles(Level::new(), DECODE_BASE, &DECODE_PAIR, encoded, &mut out) {
        out.fill(b'N');
    }
    out
}

/// BAM 4-bit encoding → Base discriminant lookup.
/// Maps A(1), C(2), G(4), T(8) to their ASCII/Base values;
/// all other nibbles (=, IUPAC ambiguity, N) → Unknown(78).
// r[impl base_decode.table]
// r[impl base_decode.table_invariant]
// Indexed by the 4-bit BAM nibble value (0–15):
// 0:= 1:A 2:C 3:M 4:G 5:R 6:S 7:V 8:T 9:W 10:Y 11:H 12:K 13:D 14:B 15:N
#[allow(clippy::byte_char_slices, reason = "per-index comments explain the IUPAC nibble mapping")]
static DECODE_BASE_TYPED: &[u8; 16] = &[
    b'N', // 0: = → Unknown
    b'A', // 1: A
    b'C', // 2: C
    b'N', // 3: M → Unknown
    b'G', // 4: G
    b'N', // 5: R → Unknown
    b'N', // 6: S → Unknown
    b'N', // 7: V → Unknown
    b'T', // 8: T
    b'N', // 9: W → Unknown
    b'N', // 10: Y → Unknown
    b'N', // 11: H → Unknown
    b'N', // 12: K → Unknown
    b'N', // 13: D → Unknown
    b'N', // 14: B → Unknown
    b'N', // 15: N → Unknown
];

/// Pre-computed pair table for Base-typed decoding.
// r[depends base_decode.table_invariant]
#[allow(clippy::indexing_slicing, reason = "i < 256, nibbles < 16")]
static DECODE_PAIR_TYPED: [[u8; 2]; 256] = {
    const B: [u8; 16] = *b"NACNGNNNTNNNNNNN";
    let mut table = [[0u8; 2]; 256];
    let mut i = 0;
    while i < 256 {
        table[i] = [B[i >> 4], B[i & 0xF]];
        i += 1;
    }
    table
};

/// Decode 4-bit packed bases directly into a pre-allocated buffer.
///
/// `out` must have length >= `len`. Writes exactly `len` bytes.
/// All written bytes are valid `Base` discriminants per `DECODE_BASE_TYPED`.
// r[impl bam.record.seq_4bit]
pub fn decode_bases_into(encoded: &[u8], len: usize, out: &mut [u8]) {
    debug_assert!(out.len() >= len, "output buffer too small: {} < {}", out.len(), len);
    let Some(out) = out.get_mut(..len) else { return };
    if !decode_nibbles(Level::new(), DECODE_BASE_TYPED, &DECODE_PAIR_TYPED, encoded, out) {
        out.fill(b'N');
    }
}

/// Decode `out.len()` nibbles of `encoded` through `lut`, the 16-entry nibble
/// table, and `pairs`, the same table two nibbles at a time for the scalar
/// tail. Returns `false`, having written nothing, when `encoded` is too short.
// r[impl io.simd_portable]
fn decode_nibbles(
    level: Level,
    lut: &[u8; 16],
    pairs: &[[u8; 2]; 256],
    encoded: &[u8],
    out: &mut [u8],
) -> bool {
    let Some(packed) = encoded.get(..out.len().div_ceil(2)) else { return false };
    dispatch!(level, simd => decode_nibbles_simd(simd, lut, pairs, packed, out));
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
        decode_vector(simd, table, p, o);
    }

    // A ragged tail is one more vector overlapping the last whole one: the
    // output is a pure function of the input, and the input is never written,
    // so recomputing the overlap is harmless. A 150 bp read leaves 11 of its
    // 75 packed bytes to the tail, which a byte loop would pay for.
    let tail = vec_in.remainder();
    let overlap = whole.len().checked_sub(n).filter(|_| !tail.is_empty());
    match overlap.and_then(|t| Some((whole.get(t..)?, out_pairs.get_mut(t..)?))) {
        Some((p, o)) => decode_vector(simd, table, p, o),
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
fn decode_vector<S: Simd>(simd: S, table: S::u8s, packed: &[u8], out: &mut [[u8; 2]]) {
    let p = S::u8s::from_slice(simd, packed);
    let hi = table.swizzle_dyn_within_blocks(p >> 4);
    let lo = table.swizzle_dyn_within_blocks(p & 0x0F);
    let (first, second) = hi.interleave(lo);
    let (o_first, o_second) = out.as_flattened_mut().split_at_mut(S::u8s::LEN);
    first.store_slice(o_first);
    second.store_slice(o_second);
}

/// Decode a 4-bit packed BAM sequence directly into `Base` values.
// r[impl base_decode.decode+2]
pub fn decode_bases(encoded: &[u8], len: usize) -> Vec<seqair_types::Base> {
    let bytes = decode_bases_raw(encoded, len);
    // Safety: DECODE_BASE_TYPED/DECODE_PAIR_TYPED only produce valid Base
    // discriminants (A=65, C=67, G=71, T=84, Unknown=78) in every byte.
    // r[depends base_decode.table_invariant]
    unsafe { seqair_types::Base::vec_u8_into_vec_base(bytes) }
}

/// Raw byte decode using the Base-typed lookup table.
// r[depends seq.simd_scalar_equivalence]
fn decode_bases_raw(encoded: &[u8], len: usize) -> Vec<u8> {
    let mut out = vec![0u8; len];
    decode_bases_into(encoded, len, &mut out);
    out
}

#[cfg(test)]
fn decode_bases_scalar(encoded: &[u8], len: usize) -> Vec<u8> {
    let full_bytes = len / 2;
    let mut result = vec![0u8; len];

    #[allow(clippy::indexing_slicing, reason = "bounds ensured by zip + as_chunks")]
    for (chunk, &byte) in result.as_chunks_mut::<2>().0.iter_mut().zip(&encoded[..full_bytes]) {
        let pair = DECODE_PAIR_TYPED[byte as usize];
        chunk[0] = pair[0];
        chunk[1] = pair[1];
    }

    if len % 2 == 1
        && let Some(byte) = encoded.get(full_bytes)
        && let Some(slot) = result.get_mut(len.saturating_sub(1))
    {
        #[allow(clippy::indexing_slicing, reason = "byte < 256")]
        {
            *slot = DECODE_PAIR_TYPED[*byte as usize][0];
        }
    }

    result
}

/// Encode ASCII bases to 4-bit packed BAM format, appending into `buf`.
///
/// This avoids the intermediate `Vec` allocation of [`encode_seq`].
#[allow(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "bounds ensured by step_by loop and reserve"
)]
pub fn encode_seq_into(bases: &[u8], buf: &mut Vec<u8>) {
    let n_bytes = bases.len().div_ceil(2);
    let start = buf.len();
    buf.resize(start + n_bytes, 0);

    for j in (0..bases.len()).step_by(2) {
        let hi = ENCODE_BASE[bases[j] as usize];
        let lo = if j.saturating_add(1) < bases.len() {
            ENCODE_BASE[bases[j.saturating_add(1)] as usize]
        } else {
            0
        };
        buf[start + j / 2] = (hi << 4) | lo;
    }
}

/// `Base` discriminant low nibble → 4-bit BAM code. The five discriminants
/// (A=0x41, C=0x43, G=0x47, T=0x54, N=0x4E) have distinct low nibbles, so a
/// 16-entry table indexed by `byte & 0x0F` is exact for every `Base`; the
/// unused slots are never read for a valid `Base`.
static ENCODE_BASE_LO_NIBBLE: &[u8; 16] = &{
    let mut table = [15u8; 16];
    table[0x1] = 1; // A
    table[0x3] = 2; // C
    table[0x7] = 4; // G
    table[0x4] = 8; // T
    table[0xE] = 15; // N
    table
};

/// Encode `Base`s to 4-bit packed BAM format, appending `bases.len().div_ceil(2)`
/// bytes to `buf`. Byte-identical to [`encode_seq_into`] on the same bases.
// r[impl seq.encode_bases_simd]
pub fn encode_bases_into(bases: &[seqair_types::Base], buf: &mut Vec<u8>) {
    encode_bases_at(Level::new(), bytemuck::cast_slice(bases), buf);
}

/// `bases` MUST hold only `Base` discriminants: the low-nibble table is exact
/// for those five bytes and nothing else.
// r[impl io.simd_portable]
fn encode_bases_at(level: Level, bases: &[u8], buf: &mut Vec<u8>) {
    let start = buf.len();
    buf.resize(start.saturating_add(bases.len().div_ceil(2)), 0);
    if let Some(out) = buf.get_mut(start..) {
        dispatch!(level, simd => encode_bases_simd(simd, bases, out));
    }
}

/// Two native-width vectors of bases in, one vector of packed bytes out: each
/// base is looked up by its low nibble (`pshufb` / `tbl`), the codes are
/// deinterleaved into even and odd positions, and packed `even << 4 | odd`.
/// `out.len()` is `bases.len().div_ceil(2)`.
#[simd]
#[allow(
    clippy::chunks_exact_to_as_chunks,
    reason = "`S::u8s::LEN` depends on the generic `S`, so it cannot be a const argument"
)]
#[allow(
    clippy::arithmetic_side_effects,
    reason = "lanewise shift of a 4-bit code by 4 cannot overflow a u8; LEN * 2 is a small constant"
)]
fn encode_bases_simd<S: Simd>(simd: S, bases: &[u8], out: &mut [u8]) {
    let n = S::u8s::LEN;
    let table = S::u8s::block_splat(u8x16::from_slice(simd, ENCODE_BASE_LO_NIBBLE));
    let vec_in = bases.chunks_exact(2 * n);
    let tail_in = vec_in.remainder();
    for (b, o) in vec_in.zip(out.chunks_exact_mut(n)) {
        encode_vector(simd, table, b, o);
    }
    if tail_in.is_empty() {
        return;
    }
    // A ragged tail is one more vector overlapping the last whole one — the
    // output is a pure function of the input, which is never written. It
    // must start on a pair boundary, so an odd final base is left over for
    // the scalar encoder. Shorter than two vectors: all scalar.
    let pairs = bases.len() / 2;
    let overlap = pairs.checked_sub(n).map(|p| (p * 2, p));
    match overlap.and_then(|(b, o)| Some((bases.get(b..b + 2 * n)?, out.get_mut(o..o + n)?))) {
        Some((b, o)) => {
            encode_vector(simd, table, b, o);
            if bases.len() % 2 == 1 {
                let last = bases.len() - 1;
                encode_seq_scalar_into(
                    bases.get(last..).unwrap_or_default(),
                    out.get_mut(last / 2..).unwrap_or_default(),
                );
            }
        }
        None => {
            let done = bases.len() - tail_in.len();
            encode_seq_scalar_into(tail_in, out.get_mut(done / 2..).unwrap_or_default());
        }
    }
}

/// One step of [`encode_bases_simd`]: `2 * LEN` bases into `LEN` packed
/// bytes. A function, not a closure, so it keeps the level's target features.
#[inline(always)]
#[allow(clippy::arithmetic_side_effects, reason = "a 4-bit code shifted by 4 fits a u8")]
fn encode_vector<S: Simd>(simd: S, table: S::u8s, bases: &[u8], out: &mut [u8]) {
    let (b0, b1) = bases.split_at(S::u8s::LEN);
    let c0 = table.swizzle_dyn_within_blocks(S::u8s::from_slice(simd, b0) & 0x0F);
    let c1 = table.swizzle_dyn_within_blocks(S::u8s::from_slice(simd, b1) & 0x0F);
    let (even, odd) = c0.deinterleave(c1);
    ((even << 4) | odd).store_slice(out);
}

/// The 256-entry-table encoder over a pre-sized `out` of
/// `bases.len().div_ceil(2)` bytes: the tail of the vector kernel.
#[allow(clippy::indexing_slicing, reason = "a `u8` index is < 256 = ENCODE_BASE.len()")]
#[allow(clippy::arithmetic_side_effects, reason = "a 4-bit code shifted by 4 fits a u8")]
fn encode_seq_scalar_into(bases: &[u8], out: &mut [u8]) {
    let (pairs, odd) = bases.as_chunks::<2>();
    for (o, &[hi, lo]) in out.iter_mut().zip(pairs) {
        *o = (ENCODE_BASE[usize::from(hi)] << 4) | ENCODE_BASE[usize::from(lo)];
    }
    if let [last] = odd
        && let Some(o) = out.get_mut(pairs.len())
    {
        *o = ENCODE_BASE[usize::from(*last)] << 4;
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

    /// Valid `Base` discriminants, derived from the enum so the test stays correct
    /// if discriminant values ever change.
    const VALID_DISCRIMINANTS: [u8; 5] = [
        seqair_types::Base::A as u8,
        seqair_types::Base::C as u8,
        seqair_types::Base::G as u8,
        seqair_types::Base::T as u8,
        seqair_types::Base::Unknown as u8,
    ];

    fn is_valid_base_discriminant(byte: u8) -> bool {
        VALID_DISCRIMINANTS.contains(&byte)
    }

    // r[verify base_decode.table_invariant]
    #[test]
    fn decode_base_typed_entries_are_valid_base_discriminants() {
        for (i, &byte) in DECODE_BASE_TYPED.iter().enumerate() {
            assert!(
                is_valid_base_discriminant(byte),
                "DECODE_BASE_TYPED[{i}] = {byte} (0x{byte:02x}) is not a valid Base discriminant"
            );
        }
    }

    // r[verify base_decode.table_invariant]
    #[test]
    fn decode_pair_typed_entries_are_valid_base_discriminants() {
        for (i, pair) in DECODE_PAIR_TYPED.iter().enumerate() {
            assert!(
                is_valid_base_discriminant(pair[0]),
                "DECODE_PAIR_TYPED[{i}][0] = {} (0x{:02x}) is not a valid Base discriminant",
                pair[0],
                pair[0]
            );
            assert!(
                is_valid_base_discriminant(pair[1]),
                "DECODE_PAIR_TYPED[{i}][1] = {} (0x{:02x}) is not a valid Base discriminant",
                pair[1],
                pair[1]
            );
        }
    }

    // r[verify seq.simd_scalar_equivalence]
    #[test]
    fn decode_bases_raw_matches_scalar_for_all_byte_values() {
        // Build an input containing all 256 byte values
        let input: Vec<u8> = (0u16..=255).map(|b| b as u8).collect();
        // Decode with scalar path as the reference
        let expected = decode_bases_scalar(&input, 512);
        // Decode with the dispatched path (SIMD on supported platforms)
        let actual = decode_bases_raw(&input, 512);
        assert_eq!(actual, expected, "SIMD and scalar paths diverge for all-byte-values input");
    }

    // r[verify base_decode.decode+2]
    #[test]
    fn decode_bases_into_matches_decode_bases() {
        for seq_len in [0usize, 1, 2, 15, 16, 17, 31, 32, 33, 63, 64, 65, 128, 150] {
            let encoded_len = seq_len.div_ceil(2);
            let input: Vec<u8> = (0..encoded_len).map(|i| (i & 0xFF) as u8).collect();
            let expected = decode_bases(&input, seq_len);
            let mut out = vec![0u8; seq_len];
            decode_bases_into(&input, seq_len, &mut out);
            // Safety: decode_bases_into writes only valid Base discriminants (A=65, C=67, G=71, T=84, Unknown=78)
            // as guaranteed by the DECODE_BASE_TYPED table. Base is repr(u8), so reinterpreting Vec<u8> as Vec<Base> is sound.
            let actual = unsafe { seqair_types::Base::vec_u8_into_vec_base(out) };
            assert_eq!(actual, expected, "mismatch at seq_len={seq_len}");
        }
    }

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
    // r[verify io.simd_portable]
    #[hegel::test]
    fn every_level_matches_nibble_oracle(tc: TestCase) {
        let packed = tc.draw(gs::binary().max_size(300));
        // Up to one nibble fewer than `packed` holds, so an odd `len` leaves
        // the low nibble of the last byte unread.
        let len = tc.draw(gs::integers::<usize>().max_value(packed.len() * 2));
        for (lut, pairs) in [(DECODE_BASE, &DECODE_PAIR), (DECODE_BASE_TYPED, &DECODE_PAIR_TYPED)] {
            let expected = nibble_oracle(lut, &packed, len);
            for level in crate::simd_levels::levels() {
                let mut out = vec![0u8; len];
                assert!(decode_nibbles(level, lut, pairs, &packed, &mut out));
                assert_eq!(out, expected, "{level:?}, len={len}");
            }
        }
    }

    fn to_bases(codes: Vec<u8>) -> Vec<seqair_types::Base> {
        use seqair_types::Base;
        codes
            .into_iter()
            .map(|c| [Base::A, Base::C, Base::G, Base::T, Base::Unknown][usize::from(c)])
            .collect()
    }

    /// Independent oracle: the spec's nibble values per base, packed one pair
    /// at a time — shares neither table with the kernel.
    fn encode_oracle(bases: &[seqair_types::Base]) -> Vec<u8> {
        use seqair_types::Base;
        let code = |b: &Base| match b {
            Base::A => 1u8,
            Base::C => 2,
            Base::G => 4,
            Base::T => 8,
            Base::Unknown => 15,
        };
        bases.chunks(2).map(|p| (code(&p[0]) << 4) | p.get(1).map_or(0, code)).collect()
    }

    // r[verify seq.encode_bases_simd]
    // r[verify io.simd_portable]
    #[hegel::test]
    fn encode_bases_every_level_matches_oracle(tc: TestCase) {
        let seq = to_bases(tc.draw(gs::vecs(gs::integers::<u8>().max_value(4)).max_size(400)));
        let expected = encode_oracle(&seq);
        let bytes: &[u8] = bytemuck::cast_slice(&seq);
        assert_eq!(encode_seq(bytes), expected, "256-entry table");
        for level in crate::simd_levels::levels() {
            // A non-empty prefix proves the kernel appends rather than overwrites.
            let mut buf = vec![0xAB];
            encode_bases_at(level, bytes, &mut buf);
            assert_eq!(buf[0], 0xAB);
            assert_eq!(&buf[1..], expected.as_slice(), "{level:?}, len={}", seq.len());
        }
    }

    // r[verify seq.encode_bases_simd]
    #[test]
    fn encode_bases_boundary_lengths_round_trip() {
        use seqair_types::Base;
        let cycle = [Base::A, Base::C, Base::G, Base::T, Base::Unknown];
        for len in [0usize, 1, 2, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127, 128, 129, 150] {
            let seq: Vec<Base> = cycle.iter().copied().cycle().take(len).collect();
            let mut buf = Vec::new();
            encode_bases_into(&seq, &mut buf);
            assert_eq!(decode_bases(&buf, len), seq, "len={len}");
        }
    }

    #[test]
    fn short_input_decodes_to_n() {
        assert_eq!(decode_seq(&[0x12], 3), b"NNN");
        assert_eq!(decode_bases_raw(&[0x12], 3), b"NNN");
    }

    // r[verify seq.simd_scalar_equivalence]
    #[test]
    fn decode_bases_raw_matches_scalar_at_simd_boundaries() {
        for seq_len in [0usize, 1, 2, 15, 16, 17, 31, 32, 33, 63, 64, 65, 128] {
            let encoded_len = seq_len.div_ceil(2);
            let input: Vec<u8> = (0..encoded_len).map(|i| (i & 0xFF) as u8).collect();
            let expected = decode_bases_scalar(&input, seq_len);
            let actual = decode_bases_raw(&input, seq_len);
            assert_eq!(actual, expected, "mismatch at seq_len={seq_len}");
        }
    }
}
