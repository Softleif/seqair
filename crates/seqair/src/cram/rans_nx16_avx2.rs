//! AVX2-accelerated inner loop for rANS Nx16 order-0 32-state decode.
#![allow(clippy::indexing_slicing, reason = "SIMD lane indices bounded by j<32 invariant")]

#[cfg(target_arch = "x86_64")]
use super::reader::CramError;

#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

/// AVX2 32-state order-0 inner loop.
///
/// # Safety
///
/// Caller must verify AVX2 via `is_x86_feature_detected!("avx2")`.
/// `states.len() == 32`, all table references must be valid.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
pub(crate) unsafe fn decode_32state_loop(
    src: &mut &[u8],
    dst: &mut [u8],
    frequencies: &[u32; 256],
    cumulative_frequencies: &[u32; 256],
    sym_table: &[u8; 4096],
    states: &mut [u32],
) -> Result<(), CramError> {
    let full_chunks = dst.len() / 32;
    let mask_12bit = _mm256_set1_epi32(0xFFF);

    for chunk_idx in 0..full_chunks {
        let start = chunk_idx.wrapping_mul(32);
        let chunk = dst
            .get_mut(start..start.wrapping_add(32))
            .ok_or(CramError::Truncated { context: "avx2 chunk range" })?;

        // Scalar symbol lookup
        for j in 0..32 {
            let f = states.get(j).unwrap_or(&0) & 0xFFF;
            let sym = sym_table
                .get(f as usize)
                .ok_or(CramError::Truncated { context: "avx2 sym_table" })?;
            *chunk.get_mut(j).unwrap_or(&mut 0) = *sym;
        }

        // SIMD state update (4 groups of 8)
        for j in (0..32).step_by(8) {
            // Safety: states has ≥32 elements, j+8 ≤ 32.
            let s = unsafe { _mm256_loadu_si256(states.as_ptr().add(j) as *const __m256i) };
            let hi = _mm256_srli_epi32(s, 12);
            let f = _mm256_and_si256(s, mask_12bit);

            let syms = [
                usize::from(chunk[j]),
                usize::from(chunk[j.wrapping_add(1)]),
                usize::from(chunk[j.wrapping_add(2)]),
                usize::from(chunk[j.wrapping_add(3)]),
                usize::from(chunk[j.wrapping_add(4)]),
                usize::from(chunk[j.wrapping_add(5)]),
                usize::from(chunk[j.wrapping_add(6)]),
                usize::from(chunk[j.wrapping_add(7)]),
            ];
            let freq_arr = [
                frequencies[syms[0]],
                frequencies[syms[1]],
                frequencies[syms[2]],
                frequencies[syms[3]],
                frequencies[syms[4]],
                frequencies[syms[5]],
                frequencies[syms[6]],
                frequencies[syms[7]],
            ];
            let cum_arr = [
                cumulative_frequencies[syms[0]],
                cumulative_frequencies[syms[1]],
                cumulative_frequencies[syms[2]],
                cumulative_frequencies[syms[3]],
                cumulative_frequencies[syms[4]],
                cumulative_frequencies[syms[5]],
                cumulative_frequencies[syms[6]],
                cumulative_frequencies[syms[7]],
            ];

            // Safety: freq_arr and cum_arr are stack arrays with 8 valid u32 elements.
            let (freqs, cums) = unsafe {
                (
                    _mm256_loadu_si256(freq_arr.as_ptr() as *const __m256i),
                    _mm256_loadu_si256(cum_arr.as_ptr() as *const __m256i),
                )
            };
            let prod = _mm256_mullo_epi32(freqs, hi);
            let sum = _mm256_add_epi32(prod, f);
            let new_s = _mm256_sub_epi32(sum, cums);
            // Safety: states has ≥32 elements, j+8 ≤ 32.
            unsafe { _mm256_storeu_si256(states.as_mut_ptr().add(j) as *mut __m256i, new_s) };
        }

        // Scalar renormalization
        for state in states.iter_mut() {
            if *state < (1 << 15) {
                let lo = u32::from(
                    super::rans_nx16::read_u16_le_prv(src)
                        .ok_or(CramError::Truncated { context: "avx2 renorm" })?,
                );
                *state = state.wrapping_shl(16) | lo;
                if *state < (1 << 15) {
                    let lo = u32::from(
                        super::rans_nx16::read_u16_le_prv(src)
                            .ok_or(CramError::Truncated { context: "avx2 renorm" })?,
                    );
                    *state = state.wrapping_shl(16) | lo;
                    while *state < (1 << 15) {
                        let lo = u32::from(
                            super::rans_nx16::read_u16_le_prv(src)
                                .ok_or(CramError::Truncated { context: "avx2 renorm" })?,
                        );
                        *state = state.wrapping_shl(16) | lo;
                    }
                }
            }
        }
    }

    // Scalar remainder
    let remainder_start = full_chunks.wrapping_mul(32);
    if let Some(remainder) = dst.get_mut(remainder_start..) {
        for d in remainder.iter_mut() {
            let f = states[0] & 0xFFF;
            let sym = sym_table
                .get(f as usize)
                .ok_or(CramError::Truncated { context: "avx2 remainder sym_table" })?;
            *d = *sym;
            let i = usize::from(*sym);
            states[0] = super::rans_nx16::state_step(
                states[0],
                frequencies[i],
                cumulative_frequencies[i],
                super::rans_nx16::ORDER_0_BITS,
            );
            states[0] = super::rans_nx16::state_renormalize(states[0], src)
                .ok_or(CramError::Truncated { context: "avx2 remainder renorm" })?;
        }
    }

    Ok(())
}
