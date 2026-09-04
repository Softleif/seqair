//! Vectorized newline-strip + ASCII-uppercase for fetched FASTA/FASTQ spans.
//!
//! Line terminators are sparse (one or two per line), so instead of a general
//! byte-compress the pass is split into three vectorized steps: a memchr-style
//! scan for the next `\n`/`\r`, a `memcpy` of the run between delimiters, and
//! an in-place uppercase of the copied run. Uppercasing is range-masked to
//! `a..=z` (0x61..=0x7A) — every other byte, including `>= 0x80`, passes
//! through unchanged, matching `u8::to_ascii_uppercase` exactly.

/// Strip `\n`/`\r` from `raw`, appending the rest ASCII-uppercased to `out`.
///
/// Byte-for-byte equivalent to the scalar loop this replaced: a byte equal to
/// `b'\n'` or `b'\r'` is dropped, anything else goes through
/// `to_ascii_uppercase`. Consecutive delimiters (`\r\n`) yield empty runs,
/// which the copy turns into the no-op the old loop's skip was.
pub(super) fn strip_newlines_uppercase(raw: &[u8], out: &mut Vec<u8>) {
    let mut search = 0usize;
    while let Some(hit) = find_newline(raw, search) {
        if let Some(run) = raw.get(search..hit) {
            append_uppercased(run, out);
        }
        search = hit.saturating_add(1);
    }
    if let Some(tail) = raw.get(search..) {
        append_uppercased(tail, out);
    }
}

/// Copy `run` onto `out` and uppercase the copy in place, so the vector
/// kernel writes each byte once and `raw` is never mutated.
fn append_uppercased(run: &[u8], out: &mut Vec<u8>) {
    let start = out.len();
    out.extend_from_slice(run);
    if let Some(appended) = out.get_mut(start..) {
        uppercase_in_place(appended);
    }
}

/// Index of the first `\n` or `\r` at or after `from`, if any.
fn find_newline(raw: &[u8], from: usize) -> Option<usize> {
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") {
            // Safety: AVX2 verified by feature detection.
            return unsafe { find_newline_avx2(raw, from) };
        }
        if is_x86_feature_detected!("ssse3") {
            // Safety: SSSE3 verified by feature detection.
            return unsafe { find_newline_ssse3(raw, from) };
        }
    }

    #[cfg(target_arch = "aarch64")]
    {
        // Safety: NEON is always available on aarch64.
        return unsafe { find_newline_neon(raw, from) };
    }

    #[cfg_attr(
        target_arch = "aarch64",
        expect(unreachable_code, reason = "NEON return above makes this dead on aarch64")
    )]
    find_newline_scalar(raw, from)
}

/// Uppercase `a..=z` in place, leaving every other byte unchanged.
fn uppercase_in_place(bytes: &mut [u8]) {
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") {
            // Safety: AVX2 verified by feature detection.
            unsafe { uppercase_avx2(bytes) };
            return;
        }
        if is_x86_feature_detected!("ssse3") {
            // Safety: SSSE3 verified by feature detection.
            unsafe { uppercase_ssse3(bytes) };
            return;
        }
    }

    #[cfg(target_arch = "aarch64")]
    {
        // Safety: NEON is always available on aarch64.
        unsafe { uppercase_neon(bytes) };
        return;
    }

    #[cfg_attr(
        target_arch = "aarch64",
        expect(unreachable_code, reason = "NEON return above makes this dead on aarch64")
    )]
    uppercase_scalar(bytes);
}

/// Scalar fallback: the exact predicate of the loop this replaced.
fn find_newline_scalar(raw: &[u8], from: usize) -> Option<usize> {
    let rel = raw.get(from..)?.iter().position(|&b| b == b'\n' || b == b'\r')?;
    Some(from.saturating_add(rel))
}

/// Scalar fallback: exact `u8::to_ascii_uppercase` semantics.
fn uppercase_scalar(bytes: &mut [u8]) {
    for b in bytes.iter_mut() {
        b.make_ascii_uppercase();
    }
}

/// AVX2 memchr2: 32 bytes per iteration, `movemask` + trailing-zeros locate
/// the first `\n`/`\r`.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn find_newline_avx2(raw: &[u8], from: usize) -> Option<usize> {
    use std::arch::x86_64::*;

    // Safety: AVX2 intrinsics require AVX2 availability (ensured by #[target_feature]).
    let v_lf = _mm256_set1_epi8(b'\n'.cast_signed());
    let v_cr = _mm256_set1_epi8(b'\r'.cast_signed());

    let len = raw.len();
    let ptr = raw.as_ptr();
    let mut i = from;

    while let Some(next) = i.checked_add(32)
        && next <= len
    {
        // Safety: `i + 32 <= len` keeps the load inside the slice.
        let mask = unsafe {
            let chunk = _mm256_loadu_si256(ptr.add(i) as *const __m256i);
            let hit =
                _mm256_or_si256(_mm256_cmpeq_epi8(chunk, v_lf), _mm256_cmpeq_epi8(chunk, v_cr));
            _mm256_movemask_epi8(hit)
        };
        if mask != 0 {
            // mask != 0 ⇒ trailing_zeros() < 32, so the hit is inside the chunk.
            return Some(i.saturating_add(mask.trailing_zeros() as usize));
        }
        i = next;
    }

    find_newline_scalar(raw, i)
}

/// SSSE3 memchr2: 16 bytes per iteration. Every instruction used is actually
/// SSE2-era; gating on SSSE3 matches the rest of the codebase's dispatch.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "ssse3")]
unsafe fn find_newline_ssse3(raw: &[u8], from: usize) -> Option<usize> {
    use std::arch::x86_64::*;

    // Safety: SSSE3 intrinsics require SSSE3 availability (ensured by #[target_feature]).
    let v_lf = _mm_set1_epi8(b'\n'.cast_signed());
    let v_cr = _mm_set1_epi8(b'\r'.cast_signed());

    let len = raw.len();
    let ptr = raw.as_ptr();
    let mut i = from;

    while let Some(next) = i.checked_add(16)
        && next <= len
    {
        // Safety: `i + 16 <= len` keeps the load inside the slice.
        let mask = unsafe {
            let chunk = _mm_loadu_si128(ptr.add(i) as *const __m128i);
            let hit = _mm_or_si128(_mm_cmpeq_epi8(chunk, v_lf), _mm_cmpeq_epi8(chunk, v_cr));
            _mm_movemask_epi8(hit)
        };
        if mask != 0 {
            // mask != 0 ⇒ trailing_zeros() < 16, so the hit is inside the chunk.
            return Some(i.saturating_add(mask.trailing_zeros() as usize));
        }
        i = next;
    }

    find_newline_scalar(raw, i)
}

/// NEON memchr2: 16 bytes per iteration. NEON has no `movemask`, so the
/// 0x00/0xFF compare lanes are narrowed one nibble per byte into a u64 —
/// bit `4k` set ⇔ byte `k` hit — and trailing-zeros recovers the index.
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn find_newline_neon(raw: &[u8], from: usize) -> Option<usize> {
    use std::arch::aarch64::*;

    let v_lf = vdupq_n_u8(b'\n');
    let v_cr = vdupq_n_u8(b'\r');

    let len = raw.len();
    let ptr = raw.as_ptr();
    let mut i = from;

    while let Some(next) = i.checked_add(16)
        && next <= len
    {
        // Safety: `i + 16 <= len` keeps the load inside the slice.
        let mask = unsafe {
            let chunk = vld1q_u8(ptr.add(i));
            let hit = vorrq_u8(vceqq_u8(chunk, v_lf), vceqq_u8(chunk, v_cr));
            let nibbles = vshrn_n_u16::<4>(vreinterpretq_u16_u8(hit));
            vget_lane_u64::<0>(vreinterpret_u64_u8(nibbles))
        };
        if mask != 0 {
            // Each byte became 4 bits, so the byte index is bit_index / 4.
            // mask != 0 ⇒ trailing_zeros() < 64 ⇒ the quotient is < 16.
            #[allow(
                clippy::arithmetic_side_effects,
                reason = "trailing_zeros of a nonzero u64 is < 64, so >> 2 cannot overflow"
            )]
            let offset = (mask.trailing_zeros() >> 2) as usize;
            return Some(i.saturating_add(offset));
        }
        i = next;
    }

    find_newline_scalar(raw, i)
}

/// AVX2 uppercase: 32 bytes per iteration, range-masked to `a..=z`.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn uppercase_avx2(bytes: &mut [u8]) {
    use std::arch::x86_64::*;

    // Safety: AVX2 intrinsics require AVX2 availability (ensured by #[target_feature]).
    let v_a = _mm256_set1_epi8(b'a'.cast_signed());
    let v_z = _mm256_set1_epi8(b'z'.cast_signed());
    let v_flip = _mm256_set1_epi8(0x20);

    let len = bytes.len();
    let ptr = bytes.as_mut_ptr();
    let mut i = 0usize;

    while let Some(next) = i.checked_add(32)
        && next <= len
    {
        // Safety: `i + 32 <= len` keeps the load and store inside the slice.
        unsafe {
            let chunk = _mm256_loadu_si256(ptr.add(i) as *const __m256i);
            // Unsigned range test: max(chunk,'a') == chunk ⇔ chunk >= 'a',
            // min(chunk,'z') == chunk ⇔ chunk <= 'z'. Only that range earns
            // the 0x20 flip — NOT the `& 0xDF` blanket uppercase, which would
            // also mangle `[`, `{`, `~` and bytes >= 0x80.
            let ge_a = _mm256_cmpeq_epi8(_mm256_max_epu8(chunk, v_a), chunk);
            let le_z = _mm256_cmpeq_epi8(_mm256_min_epu8(chunk, v_z), chunk);
            let flip = _mm256_and_si256(_mm256_and_si256(ge_a, le_z), v_flip);
            _mm256_storeu_si256(ptr.add(i) as *mut __m256i, _mm256_xor_si256(chunk, flip));
        }
        i = next;
    }

    if let Some(tail) = bytes.get_mut(i..) {
        uppercase_scalar(tail);
    }
}

/// SSSE3 uppercase: 16 bytes per iteration, range-masked to `a..=z`.
/// (SSE2-era instructions only; SSSE3 gating matches the codebase.)
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "ssse3")]
unsafe fn uppercase_ssse3(bytes: &mut [u8]) {
    use std::arch::x86_64::*;

    // Safety: SSSE3 intrinsics require SSSE3 availability (ensured by #[target_feature]).
    let v_a = _mm_set1_epi8(b'a'.cast_signed());
    let v_z = _mm_set1_epi8(b'z'.cast_signed());
    let v_flip = _mm_set1_epi8(0x20);

    let len = bytes.len();
    let ptr = bytes.as_mut_ptr();
    let mut i = 0usize;

    while let Some(next) = i.checked_add(16)
        && next <= len
    {
        // Safety: `i + 16 <= len` keeps the load and store inside the slice.
        unsafe {
            let chunk = _mm_loadu_si128(ptr.add(i) as *const __m128i);
            // Same unsigned range test as the AVX2 kernel above.
            let ge_a = _mm_cmpeq_epi8(_mm_max_epu8(chunk, v_a), chunk);
            let le_z = _mm_cmpeq_epi8(_mm_min_epu8(chunk, v_z), chunk);
            let flip = _mm_and_si128(_mm_and_si128(ge_a, le_z), v_flip);
            _mm_storeu_si128(ptr.add(i) as *mut __m128i, _mm_xor_si128(chunk, flip));
        }
        i = next;
    }

    if let Some(tail) = bytes.get_mut(i..) {
        uppercase_scalar(tail);
    }
}

/// NEON uppercase: 16 bytes per iteration, range-masked to `a..=z`.
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn uppercase_neon(bytes: &mut [u8]) {
    use std::arch::aarch64::*;

    let v_a = vdupq_n_u8(b'a');
    let v_z = vdupq_n_u8(b'z');
    let v_flip = vdupq_n_u8(0x20);

    let len = bytes.len();
    let ptr = bytes.as_mut_ptr();
    let mut i = 0usize;

    while let Some(next) = i.checked_add(16)
        && next <= len
    {
        // Safety: `i + 16 <= len` keeps the load and store inside the slice.
        unsafe {
            let chunk = vld1q_u8(ptr.add(i));
            // Unsigned range test, same shape as the x86 kernels: only
            // 0x61..=0x7A earns the 0x20 flip, all other bytes pass through.
            let is_lower = vandq_u8(vcgeq_u8(chunk, v_a), vcgeq_u8(v_z, chunk));
            vst1q_u8(ptr.add(i), veorq_u8(chunk, vandq_u8(is_lower, v_flip)));
        }
        i = next;
    }

    if let Some(tail) = bytes.get_mut(i..) {
        uppercase_scalar(tail);
    }
}

#[cfg(test)]
#[allow(clippy::arithmetic_side_effects, reason = "test arithmetic is not safety-critical")]
mod tests {
    use super::*;
    use proptest::prelude::*;

    /// Independent scalar oracle: the exact semantics of the loop this
    /// replaced, written as plainly as possible.
    fn strip_upper_reference(raw: &[u8]) -> Vec<u8> {
        let mut out = Vec::with_capacity(raw.len());
        for &b in raw {
            if b != b'\n' && b != b'\r' {
                out.push(b.to_ascii_uppercase());
            }
        }
        out
    }

    fn run(raw: &[u8]) -> Vec<u8> {
        let mut out = Vec::with_capacity(raw.len());
        strip_newlines_uppercase(raw, &mut out);
        out
    }

    // r[verify fasta.fetch.newline_stripping]
    // r[verify fasta.fetch.uppercase]
    #[test]
    fn targeted_vectors_match_reference() {
        let mut cases: Vec<Vec<u8>> = vec![
            Vec::new(),
            b"\n".to_vec(),
            b"\r".to_vec(),
            b"\r\n".to_vec(),
            b"\n\r".to_vec(),
            b"\nacgt".to_vec(), // delimiter at offset 0
            b"acgt\n".to_vec(), // delimiter as the final byte
            b"acgt".to_vec(),   // no delimiter at all
            b"acgtnryswkmbdhv".to_vec(),
            b"ACGTNRYSWKMBDHV".to_vec(),
            (0u8..=255).collect(),
            vec![b'\n'; 100],
            b"\r\n".repeat(50),
            b"ACGT\n".repeat(40), // realistic FASTA line structure
            vec![0x80, 0xFF, 0xFE, b'[', b'{', b'~', b'@', b'+', b'`'],
        ];
        // Straddle the 16/32-byte vector bodies and their scalar tails, with
        // the delimiter leading, trailing, and absent.
        for len in [1usize, 2, 15, 16, 17, 31, 32, 33, 63, 64, 65, 100, 1024] {
            cases.push(vec![b'g'; len]);
            let mut lead = vec![b'g'; len];
            if let Some(first) = lead.first_mut() {
                *first = b'\n';
            }
            cases.push(lead);
            let mut trail = vec![b'g'; len];
            if let Some(last) = trail.last_mut() {
                *last = b'\r';
            }
            cases.push(trail);
        }
        for case in &cases {
            assert_eq!(run(case), strip_upper_reference(case), "input: {case:?}");
        }
    }

    // r[verify fasta.fetch.newline_stripping]
    // r[verify fasta.fetch.uppercase]
    proptest! {
        /// Every byte value at every length, covering the vector bodies, their
        /// scalar tails, and delimiters at arbitrary positions: the dispatched
        /// SIMD path must equal the oracle bit-for-bit.
        #[test]
        fn simd_matches_scalar_reference(raw in proptest::collection::vec(any::<u8>(), 0..4096)) {
            prop_assert_eq!(run(&raw), strip_upper_reference(&raw));
        }
    }
}
