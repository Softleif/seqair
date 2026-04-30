//! NEON-accelerated inner loop for rANS Nx16 order-0 32-state decode.
#![allow(clippy::indexing_slicing, reason = "SIMD lane indices bounded by j<32 invariant")]

#[cfg(target_arch = "aarch64")]
use super::reader::CramError;

#[cfg(target_arch = "aarch64")]
use std::arch::aarch64::*;

/// NEON 32-state order-0 inner loop.
///
/// # Safety
///
/// Caller guarantees NEON availability (always true on aarch64),
/// `states.len() == 32`, and all table references are valid.
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
pub(crate) unsafe fn decode_32state_loop(
    src: &mut &[u8],
    dst: &mut [u8],
    frequencies: &[u32; 256],
    cumulative_frequencies: &[u32; 256],
    sym_table: &[u8; 4096],
    states: &mut [u32],
) -> Result<(), CramError> {
    let full_chunks = dst.len() / 32;
    let mask_12bit = vdupq_n_u32(0xFFF);

    for chunk_idx in 0..full_chunks {
        let start = chunk_idx.wrapping_mul(32);
        let chunk = dst
            .get_mut(start..start.wrapping_add(32))
            .ok_or(CramError::Truncated { context: "neon chunk range" })?;

        // Scalar symbol lookup
        for j in 0..32 {
            let f = states.get(j).unwrap_or(&0) & 0xFFF;
            let sym = sym_table
                .get(f as usize)
                .ok_or(CramError::Truncated { context: "neon sym_table" })?;
            *chunk.get_mut(j).unwrap_or(&mut 0) = *sym;
        }

        // SIMD state update (8 groups of 4)
        for j in (0..32).step_by(4) {
            // Safety: states has ≥32 elements, j+4 ≤ 32.
            let s = unsafe { vld1q_u32(states.as_ptr().add(j)) };
            let hi = vshrq_n_u32(s, 12);
            let f = vandq_u32(s, mask_12bit);

            let syms = [
                usize::from(chunk[j]),
                usize::from(chunk[j.wrapping_add(1)]),
                usize::from(chunk[j.wrapping_add(2)]),
                usize::from(chunk[j.wrapping_add(3)]),
            ];
            let freq_arr = [
                frequencies[syms[0]],
                frequencies[syms[1]],
                frequencies[syms[2]],
                frequencies[syms[3]],
            ];
            let cum_arr = [
                cumulative_frequencies[syms[0]],
                cumulative_frequencies[syms[1]],
                cumulative_frequencies[syms[2]],
                cumulative_frequencies[syms[3]],
            ];

            // Safety: freq_arr and cum_arr are stack arrays with 4 valid u32 elements.
            let (freqs, cums) =
                unsafe { (vld1q_u32(freq_arr.as_ptr()), vld1q_u32(cum_arr.as_ptr())) };
            let new_s = vsubq_u32(vaddq_u32(vmulq_u32(freqs, hi), f), cums);
            // Safety: states has ≥32 elements, j+4 ≤ 32.
            unsafe { vst1q_u32(states.as_mut_ptr().add(j), new_s) };
        }

        // Scalar renormalization
        for state in states.iter_mut() {
            if *state < (1 << 15) {
                let lo = u32::from(
                    super::rans_nx16::read_u16_le_prv(src)
                        .ok_or(CramError::Truncated { context: "neon renorm" })?,
                );
                *state = state.wrapping_shl(16) | lo;
                if *state < (1 << 15) {
                    let lo = u32::from(
                        super::rans_nx16::read_u16_le_prv(src)
                            .ok_or(CramError::Truncated { context: "neon renorm" })?,
                    );
                    *state = state.wrapping_shl(16) | lo;
                    while *state < (1 << 15) {
                        let lo = u32::from(
                            super::rans_nx16::read_u16_le_prv(src)
                                .ok_or(CramError::Truncated { context: "neon renorm" })?,
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
                .ok_or(CramError::Truncated { context: "neon remainder sym_table" })?;
            *d = *sym;
            let i = usize::from(*sym);
            states[0] = super::rans_nx16::state_step(
                states[0],
                frequencies[i],
                cumulative_frequencies[i],
                super::rans_nx16::ORDER_0_BITS,
            );
            states[0] = super::rans_nx16::state_renormalize(states[0], src)
                .ok_or(CramError::Truncated { context: "neon remainder renorm" })?;
        }
    }

    Ok(())
}
