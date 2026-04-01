//! Encode and decode BAM 4-bit packed sequences. [`decode_seq`] and [`encode_seq`] dispatch to
//! SSSE3 (x86_64), NEON (aarch64), or scalar paths; all produce identical output.

// r[impl seq.decode_scalar]
// r[impl seq.decode_pair_table]
// r[impl seq.decode_simd]
// r[impl seq.decode_dispatch]
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
    let required = (len + 1) / 2;
    if encoded.len() < required {
        return vec![b'N'; len];
    }
    let full_bytes = len / 2;
    let mut result = vec![0u8; len];

    #[allow(clippy::indexing_slicing, reason = "bounds ensured by zip + chunks_exact")]
    for (chunk, &byte) in result[..full_bytes * 2].chunks_exact_mut(2).zip(&encoded[..full_bytes]) {
        let pair = DECODE_PAIR[byte as usize];
        chunk[0] = pair[0];
        chunk[1] = pair[1];
    }

    if len % 2 == 1
        && let Some(byte) = encoded.get(full_bytes)
        && let Some(slot) = result.get_mut(len - 1)
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
    let required = (len + 1) / 2;
    if encoded.len() < required {
        return vec![b'N'; len];
    }

    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("ssse3") {
            // Safety: we just verified SSSE3 is available.
            return unsafe { decode_seq_ssse3(encoded, len) };
        }
    }

    #[cfg(target_arch = "aarch64")]
    {
        // Safety: NEON is always available on aarch64.
        return unsafe { decode_seq_neon(encoded, len) };
    }

    #[allow(unreachable_code)]
    decode_seq_scalar(encoded, len)
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
    const B: [u8; 16] = [
        b'N', b'A', b'C', b'N', b'G', b'N', b'N', b'N', b'T', b'N', b'N', b'N', b'N', b'N', b'N',
        b'N',
    ];
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
pub fn decode_bases_into(encoded: &[u8], len: usize, out: &mut [u8]) {
    debug_assert!(out.len() >= len, "output buffer too small: {} < {}", out.len(), len);
    let required = (len + 1) / 2;
    if encoded.len() < required || out.len() < len {
        out.get_mut(..len).map(|s| s.fill(b'N'));
        return;
    }
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("ssse3") {
            // Safety: SSSE3 verified.
            unsafe {
                decode_bases_into_ssse3(encoded, len, out);
            }
            return;
        }
    }

    #[cfg(target_arch = "aarch64")]
    {
        // Safety: NEON always available on aarch64.
        unsafe {
            decode_bases_into_neon(encoded, len, out);
        }
        return;
    }

    #[allow(unreachable_code)]
    decode_bases_into_scalar(encoded, len, out);
}

fn decode_bases_into_scalar(encoded: &[u8], len: usize, out: &mut [u8]) {
    let full_bytes = len / 2;

    #[allow(clippy::indexing_slicing, reason = "bounds ensured by zip + chunks_exact")]
    for (chunk, &byte) in out[..full_bytes * 2].chunks_exact_mut(2).zip(&encoded[..full_bytes]) {
        let pair = DECODE_PAIR_TYPED[byte as usize];
        chunk[0] = pair[0];
        chunk[1] = pair[1];
    }

    if len % 2 == 1
        && let Some(byte) = encoded.get(full_bytes)
        && let Some(slot) = out.get_mut(len - 1)
    {
        #[allow(clippy::indexing_slicing, reason = "byte < 256")]
        {
            *slot = DECODE_PAIR_TYPED[*byte as usize][0];
        }
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "ssse3")]
#[allow(
    clippy::indexing_slicing,
    reason = "scalar tail bounds guaranteed by loop invariants; debug_asserts verify"
)]
unsafe fn decode_bases_into_ssse3(encoded: &[u8], len: usize, out: &mut [u8]) {
    use std::arch::x86_64::*;

    let full_bytes = len / 2;

    // Safety: SSSE3 availability is guaranteed by #[target_feature(enable = "ssse3")].
    // Pointer offsets stay within bounds: i + 16 <= full_bytes <= encoded.len(),
    // and o + 32 <= len <= out.len() (guaranteed by debug_assert in decode_bases_into).
    let (lut, mask_lo) = unsafe {
        (_mm_loadu_si128(DECODE_BASE_TYPED.as_ptr() as *const __m128i), _mm_set1_epi8(0x0F))
    };

    let mut i = 0;
    let mut o = 0;

    while i + 16 <= full_bytes {
        // Safety: see above; i + 16 <= full_bytes and o + 32 <= len are loop invariants.
        unsafe {
            let packed = _mm_loadu_si128(encoded.as_ptr().add(i) as *const __m128i);
            let hi = _mm_and_si128(_mm_srli_epi16(packed, 4), mask_lo);
            let lo = _mm_and_si128(packed, mask_lo);
            let decoded_hi = _mm_shuffle_epi8(lut, hi);
            let decoded_lo = _mm_shuffle_epi8(lut, lo);
            let out_a = _mm_unpacklo_epi8(decoded_hi, decoded_lo);
            let out_b = _mm_unpackhi_epi8(decoded_hi, decoded_lo);
            _mm_storeu_si128(out.as_mut_ptr().add(o) as *mut __m128i, out_a);
            _mm_storeu_si128(out.as_mut_ptr().add(o + 16) as *mut __m128i, out_b);
        }
        i += 16;
        o += 32;
    }

    debug_assert!(i <= full_bytes, "SSSE3 base loop overshot: i={i}, full_bytes={full_bytes}");
    debug_assert!(o == i * 2, "cursor invariant broken: o={o}, i*2={}", i * 2);

    while i < full_bytes {
        let pair = DECODE_PAIR_TYPED[encoded[i] as usize];
        out[o] = pair[0];
        out[o + 1] = pair[1];
        i += 1;
        o += 2;
    }

    if len % 2 == 1 {
        out[o] = DECODE_PAIR_TYPED[encoded[i] as usize][0];
    }
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn decode_bases_into_neon(encoded: &[u8], len: usize, out: &mut [u8]) {
    use std::arch::aarch64::*;

    let full_bytes = len / 2;

    // Safety: DECODE_BASE_TYPED is a &[u8; 16] static; out has at least len bytes
    // (guaranteed by debug_assert in decode_bases_into).
    // Pointer offsets stay within the allocated ranges enforced by the loop bounds.
    unsafe {
        let lut = vld1q_u8(DECODE_BASE_TYPED.as_ptr());
        let mask_lo = vdupq_n_u8(0x0F);

        let mut i = 0;
        let mut o = 0;

        while i + 16 <= full_bytes {
            let packed = vld1q_u8(encoded.as_ptr().add(i));
            let hi = vshrq_n_u8(packed, 4);
            let lo = vandq_u8(packed, mask_lo);
            let decoded_hi = vqtbl1q_u8(lut, hi);
            let decoded_lo = vqtbl1q_u8(lut, lo);
            let out_a = vzip1q_u8(decoded_hi, decoded_lo);
            let out_b = vzip2q_u8(decoded_hi, decoded_lo);
            vst1q_u8(out.as_mut_ptr().add(o), out_a);
            vst1q_u8(out.as_mut_ptr().add(o + 16), out_b);
            i += 16;
            o += 32;
        }

        debug_assert!(i <= full_bytes, "NEON base loop overshot: i={i}, full_bytes={full_bytes}");
        debug_assert!(o == i * 2, "cursor invariant broken: o={o}, i*2={}", i * 2);

        #[allow(clippy::indexing_slicing, reason = "bounds ensured by loop")]
        {
            while i < full_bytes {
                let pair = DECODE_PAIR_TYPED[encoded[i] as usize];
                out[o] = pair[0];
                out[o + 1] = pair[1];
                i += 1;
                o += 2;
            }
            if len % 2 == 1 {
                out[o] = DECODE_PAIR_TYPED[encoded[i] as usize][0];
            }
        }
    }
}

/// Decode a 4-bit packed BAM sequence directly into `Base` values.
// r[impl base_decode.decode]
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
    let required = (len + 1) / 2;
    if encoded.len() < required {
        return vec![b'N'; len];
    }

    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("ssse3") {
            // Safety: SSSE3 verified.
            return unsafe { decode_bases_ssse3(encoded, len) };
        }
    }

    #[cfg(target_arch = "aarch64")]
    {
        // Safety: NEON always available on aarch64.
        return unsafe { decode_bases_neon(encoded, len) };
    }

    #[allow(unreachable_code)]
    decode_bases_scalar(encoded, len)
}

fn decode_bases_scalar(encoded: &[u8], len: usize) -> Vec<u8> {
    let full_bytes = len / 2;
    let mut result = vec![0u8; len];

    #[allow(clippy::indexing_slicing, reason = "bounds ensured by zip + chunks_exact")]
    for (chunk, &byte) in result[..full_bytes * 2].chunks_exact_mut(2).zip(&encoded[..full_bytes]) {
        let pair = DECODE_PAIR_TYPED[byte as usize];
        chunk[0] = pair[0];
        chunk[1] = pair[1];
    }

    if len % 2 == 1
        && let Some(byte) = encoded.get(full_bytes)
        && let Some(slot) = result.get_mut(len - 1)
    {
        #[allow(clippy::indexing_slicing, reason = "byte < 256")]
        {
            *slot = DECODE_PAIR_TYPED[*byte as usize][0];
        }
    }

    result
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "ssse3")]
#[allow(
    clippy::indexing_slicing,
    reason = "scalar tail bounds guaranteed by loop invariants; debug_asserts verify"
)]
unsafe fn decode_bases_ssse3(encoded: &[u8], len: usize) -> Vec<u8> {
    use std::arch::x86_64::*;

    let full_bytes = len / 2;
    let mut result = vec![0u8; len];

    // Safety: SSSE3 availability is guaranteed by #[target_feature(enable = "ssse3")].
    // Pointer offsets stay within bounds: i + 16 <= full_bytes <= encoded.len(),
    // and o + 32 <= len = result.len().
    let (lut, mask_lo) = unsafe {
        (_mm_loadu_si128(DECODE_BASE_TYPED.as_ptr() as *const __m128i), _mm_set1_epi8(0x0F))
    };

    let mut i = 0;
    let mut o = 0;

    while i + 16 <= full_bytes {
        // Safety: see above; i + 16 <= full_bytes and o + 32 <= len are loop invariants.
        unsafe {
            let packed = _mm_loadu_si128(encoded.as_ptr().add(i) as *const __m128i);
            let hi = _mm_and_si128(_mm_srli_epi16(packed, 4), mask_lo);
            let lo = _mm_and_si128(packed, mask_lo);
            let decoded_hi = _mm_shuffle_epi8(lut, hi);
            let decoded_lo = _mm_shuffle_epi8(lut, lo);
            let out_a = _mm_unpacklo_epi8(decoded_hi, decoded_lo);
            let out_b = _mm_unpackhi_epi8(decoded_hi, decoded_lo);
            _mm_storeu_si128(result.as_mut_ptr().add(o) as *mut __m128i, out_a);
            _mm_storeu_si128(result.as_mut_ptr().add(o + 16) as *mut __m128i, out_b);
        }
        i += 16;
        o += 32;
    }

    debug_assert!(i <= full_bytes, "SSSE3 base loop overshot: i={i}, full_bytes={full_bytes}");
    debug_assert!(o == i * 2, "cursor invariant broken: o={o}, i*2={}", i * 2);

    while i < full_bytes {
        let pair = DECODE_PAIR_TYPED[encoded[i] as usize];
        result[o] = pair[0];
        result[o + 1] = pair[1];
        i += 1;
        o += 2;
    }

    if len % 2 == 1 {
        result[o] = DECODE_PAIR_TYPED[encoded[i] as usize][0];
    }

    result
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn decode_bases_neon(encoded: &[u8], len: usize) -> Vec<u8> {
    use std::arch::aarch64::*;

    let full_bytes = len / 2;
    let mut result = vec![0u8; len];

    // Safety: DECODE_BASE_TYPED is a &[u8; 16] static; result is allocated above with len bytes.
    // Pointer offsets stay within the allocated ranges enforced by the loop bounds.
    unsafe {
        let lut = vld1q_u8(DECODE_BASE_TYPED.as_ptr());
        let mask_lo = vdupq_n_u8(0x0F);

        let mut i = 0;
        let mut o = 0;

        while i + 16 <= full_bytes {
            let packed = vld1q_u8(encoded.as_ptr().add(i));
            let hi = vshrq_n_u8(packed, 4);
            let lo = vandq_u8(packed, mask_lo);
            let decoded_hi = vqtbl1q_u8(lut, hi);
            let decoded_lo = vqtbl1q_u8(lut, lo);
            let out_a = vzip1q_u8(decoded_hi, decoded_lo);
            let out_b = vzip2q_u8(decoded_hi, decoded_lo);
            vst1q_u8(result.as_mut_ptr().add(o), out_a);
            vst1q_u8(result.as_mut_ptr().add(o + 16), out_b);
            i += 16;
            o += 32;
        }

        debug_assert!(i <= full_bytes, "NEON base loop overshot: i={i}, full_bytes={full_bytes}");
        debug_assert!(o == i * 2, "cursor invariant broken: o={o}, i*2={}", i * 2);

        #[allow(clippy::indexing_slicing, reason = "bounds ensured by loop")]
        {
            while i < full_bytes {
                let pair = DECODE_PAIR_TYPED[encoded[i] as usize];
                result[o] = pair[0];
                result[o + 1] = pair[1];
                i += 1;
                o += 2;
            }
            if len % 2 == 1 {
                result[o] = DECODE_PAIR_TYPED[encoded[i] as usize][0];
            }
        }
    }

    result
}

/// Encode ASCII bases to 4-bit packed BAM format.
pub fn encode_seq(bases: &[u8]) -> Vec<u8> {
    let n_bytes = bases.len().div_ceil(2);
    let mut encoded = vec![0u8; n_bytes];

    #[allow(clippy::indexing_slicing, reason = "bounds ensured by step_by loop")]
    for j in (0..bases.len()).step_by(2) {
        let hi = ENCODE_BASE[bases[j] as usize];
        let lo = if j + 1 < bases.len() { ENCODE_BASE[bases[j + 1] as usize] } else { 0 };
        encoded[j / 2] = (hi << 4) | lo;
    }

    encoded
}

// ---- SSSE3 (x86_64) ----

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "ssse3")]
#[allow(
    clippy::indexing_slicing,
    reason = "scalar tail bounds guaranteed by loop invariants; debug_asserts verify"
)]
unsafe fn decode_seq_ssse3(encoded: &[u8], len: usize) -> Vec<u8> {
    use std::arch::x86_64::*;

    let full_bytes = len / 2;
    let mut result = vec![0u8; len];

    // Safety: SSSE3 availability is guaranteed by #[target_feature(enable = "ssse3")].
    // Pointer offsets stay within bounds: i + 16 <= full_bytes <= encoded.len(),
    // and o + 32 <= len = result.len().
    let (lut, mask_lo) =
        unsafe { (_mm_loadu_si128(DECODE_BASE.as_ptr() as *const __m128i), _mm_set1_epi8(0x0F)) };

    let mut i = 0;
    let mut o = 0;

    // 16 packed bytes → 32 decoded bases per iteration
    while i + 16 <= full_bytes {
        // Safety: see above; i + 16 <= full_bytes and o + 32 <= len are loop invariants.
        unsafe {
            let packed = _mm_loadu_si128(encoded.as_ptr().add(i) as *const __m128i);

            let hi = _mm_and_si128(_mm_srli_epi16(packed, 4), mask_lo);
            let lo = _mm_and_si128(packed, mask_lo);

            let decoded_hi = _mm_shuffle_epi8(lut, hi);
            let decoded_lo = _mm_shuffle_epi8(lut, lo);

            let out_a = _mm_unpacklo_epi8(decoded_hi, decoded_lo);
            let out_b = _mm_unpackhi_epi8(decoded_hi, decoded_lo);

            _mm_storeu_si128(result.as_mut_ptr().add(o) as *mut __m128i, out_a);
            _mm_storeu_si128(result.as_mut_ptr().add(o + 16) as *mut __m128i, out_b);
        }
        i += 16;
        o += 32;
    }

    debug_assert!(i <= full_bytes, "SSSE3 seq loop overshot: i={i}, full_bytes={full_bytes}");
    debug_assert!(o == i * 2, "cursor invariant broken: o={o}, i*2={}", i * 2);

    // Scalar tail
    while i < full_bytes {
        let pair = DECODE_PAIR[encoded[i] as usize];
        result[o] = pair[0];
        result[o + 1] = pair[1];
        i += 1;
        o += 2;
    }

    if len % 2 == 1 {
        result[o] = DECODE_PAIR[encoded[i] as usize][0];
    }

    result
}

// ---- NEON (aarch64) ----

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn decode_seq_neon(encoded: &[u8], len: usize) -> Vec<u8> {
    use std::arch::aarch64::*;

    let full_bytes = len / 2;
    let mut result = vec![0u8; len];

    // Safety: DECODE_BASE is a &[u8; 16] static; result is allocated above with len bytes.
    // Pointer offsets stay within the allocated ranges enforced by the loop bounds.
    unsafe {
        let lut = vld1q_u8(DECODE_BASE.as_ptr());
        let mask_lo = vdupq_n_u8(0x0F);

        let mut i = 0;
        let mut o = 0;

        // 16 packed bytes → 32 decoded bases per iteration
        while i + 16 <= full_bytes {
            let packed = vld1q_u8(encoded.as_ptr().add(i));

            let hi = vshrq_n_u8(packed, 4);
            let lo = vandq_u8(packed, mask_lo);

            let decoded_hi = vqtbl1q_u8(lut, hi);
            let decoded_lo = vqtbl1q_u8(lut, lo);

            let out_a = vzip1q_u8(decoded_hi, decoded_lo);
            let out_b = vzip2q_u8(decoded_hi, decoded_lo);

            vst1q_u8(result.as_mut_ptr().add(o), out_a);
            vst1q_u8(result.as_mut_ptr().add(o + 16), out_b);

            i += 16;
            o += 32;
        }

        debug_assert!(i <= full_bytes, "NEON seq loop overshot: i={i}, full_bytes={full_bytes}");
        debug_assert!(o == i * 2, "cursor invariant broken: o={o}, i*2={}", i * 2);

        #[allow(
            clippy::indexing_slicing,
            reason = "i < full_bytes ≤ encoded.len(), o/o+1 < len = result.len()"
        )]
        {
            while i < full_bytes {
                let pair = DECODE_PAIR[encoded[i] as usize];
                result[o] = pair[0];
                result[o + 1] = pair[1];
                i += 1;
                o += 2;
            }

            if len % 2 == 1 {
                result[o] = DECODE_PAIR[encoded[i] as usize][0];
            }
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

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

    // r[verify base_decode.decode]
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
