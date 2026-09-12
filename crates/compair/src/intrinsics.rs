//! The same eight lanes, written as target intrinsics rather than through
//! `wide`.
//!
//! This exists to answer one question by measurement: does the strip kernel's
//! remaining cost come from what `wide` emits, or from what the loop asks the
//! register allocator for? Every operation the [`Lane`] trait needs maps onto
//! one NEON or AVX2 instruction, so an intrinsics lane can differ from `wide`
//! only in the shapes `wide` has no way to spell -- the lane shift, the pair
//! load, and whatever the extra `#[repr]` buys the allocator. It is an
//! experiment, not a second production path: `strip_kernel` is generic over
//! the lane, so plugging this in changes no arithmetic and no order.
//!
//! Bit-parity is by construction for the same reason the `wide` lane's is --
//! every lane operation is the IEEE operation on the same values in the same
//! order -- and is asserted against the scalar kernel by the same property
//! test.

use crate::banded::LANE_MAX;

/// Eight `f32` as two 128-bit vectors, the NEON register pair `wide`'s
/// `f32x8` already is on this target -- but with the lane shift as `vextq`
/// and the window load as `vld1q` rather than as a pattern the compiler has
/// to recognise through an array.
#[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
#[derive(Clone, Copy)]
#[repr(C, align(16))]
pub(crate) struct Simd8 {
    lo: core::arch::aarch64::float32x4_t,
    hi: core::arch::aarch64::float32x4_t,
}

#[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
#[allow(
    clippy::undocumented_unsafe_blocks,
    reason = "every block is one NEON intrinsic on values this type owns; the only\
              safety condition any of them has is the target feature, which the\
              module's cfg is"
)]
mod aarch64 {
    use core::arch::aarch64::{
        float32x4_t, vaddq_f32, vaddvq_f32, vandq_u32, vbslq_f32, vceqq_f32, vcltq_f32, vdupq_n_f32,
        vextq_f32, vgetq_lane_f32, vld1q_f32, vmaxq_f32, vmaxvq_f32, vmulq_f32, vorrq_u32,
        vreinterpretq_f32_u32, vreinterpretq_u32_f32, vsetq_lane_f32, vst1q_f32,
    };

    use super::Simd8;
    use crate::banded::{Lane, Window};

    impl core::ops::Add for Simd8 {
        type Output = Self;
        #[inline(always)]
        fn add(self, rhs: Self) -> Self {
            unsafe { Self { lo: vaddq_f32(self.lo, rhs.lo), hi: vaddq_f32(self.hi, rhs.hi) } }
        }
    }

    impl core::ops::Mul for Simd8 {
        type Output = Self;
        #[inline(always)]
        fn mul(self, rhs: Self) -> Self {
            unsafe { Self { lo: vmulq_f32(self.lo, rhs.lo), hi: vmulq_f32(self.hi, rhs.hi) } }
        }
    }

    #[inline(always)]
    fn and(a: float32x4_t, b: float32x4_t) -> float32x4_t {
        unsafe {
            vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(a), vreinterpretq_u32_f32(b)))
        }
    }

    #[inline(always)]
    fn or(a: float32x4_t, b: float32x4_t) -> float32x4_t {
        unsafe {
            vreinterpretq_f32_u32(vorrq_u32(vreinterpretq_u32_f32(a), vreinterpretq_u32_f32(b)))
        }
    }

    impl Lane for Simd8 {
        const LANES: usize = 8;

        #[inline(always)]
        fn splat(value: f32) -> Self {
            unsafe { Self { lo: vdupq_n_f32(value), hi: vdupq_n_f32(value) } }
        }

        /// One `vld1q` per half. `wide` spells this as `f32x8::from([f32; 8])`,
        /// which the compiler has to see through to reach the same `ldp`.
        #[inline(always)]
        fn load(source: &Window) -> Self {
            let base = source.as_ptr();
            unsafe { Self { lo: vld1q_f32(base), hi: vld1q_f32(base.add(4)) } }
        }

        #[inline(always)]
        fn store(self, destination: &mut Window) {
            let base = destination.as_mut_ptr();
            unsafe {
                vst1q_f32(base, self.lo);
                vst1q_f32(base.add(4), self.hi);
            }
        }

        #[inline(always)]
        fn vmax(self, other: Self) -> Self {
            unsafe { Self { lo: vmaxq_f32(self.lo, other.lo), hi: vmaxq_f32(self.hi, other.hi) } }
        }

        #[inline(always)]
        fn horizontal_max(self) -> f32 {
            unsafe { vmaxvq_f32(vmaxq_f32(self.lo, self.hi)) }
        }

        #[inline(always)]
        fn equals(self, other: Self) -> Self {
            unsafe {
                Self {
                    lo: vreinterpretq_f32_u32(vceqq_f32(self.lo, other.lo)),
                    hi: vreinterpretq_f32_u32(vceqq_f32(self.hi, other.hi)),
                }
            }
        }

        #[inline(always)]
        fn either(self, other: Self) -> Self {
            Self { lo: or(self.lo, other.lo), hi: or(self.hi, other.hi) }
        }

        #[inline(always)]
        fn below(self, other: Self) -> Self {
            unsafe {
                Self {
                    lo: vreinterpretq_f32_u32(vcltq_f32(self.lo, other.lo)),
                    hi: vreinterpretq_f32_u32(vcltq_f32(self.hi, other.hi)),
                }
            }
        }

        #[inline(always)]
        fn offsets() -> Self {
            Self::load(&[0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
        }

        #[inline(always)]
        fn both(self, other: Self) -> Self {
            Self { lo: and(self.lo, other.lo), hi: and(self.hi, other.hi) }
        }

        #[inline(always)]
        fn select(self, if_true: Self, if_false: Self) -> Self {
            unsafe {
                Self {
                    lo: vbslq_f32(vreinterpretq_u32_f32(self.lo), if_true.lo, if_false.lo),
                    hi: vbslq_f32(vreinterpretq_u32_f32(self.hi), if_true.hi, if_false.hi),
                }
            }
        }

        /// Two `vextq` and one `vsetq_lane`: the upper half takes lane 3 of the
        /// lower half and the lower three of the upper, the lower half rotates
        /// and has `first` written into lane 0.
        #[inline(always)]
        fn shift_in(self, first: f32) -> Self {
            unsafe {
                Self {
                    lo: vsetq_lane_f32::<0>(first, vextq_f32::<3>(self.lo, self.lo)),
                    hi: vextq_f32::<3>(self.lo, self.hi),
                }
            }
        }

        #[inline(always)]
        fn last(self) -> f32 {
            unsafe { vgetq_lane_f32::<3>(self.hi) }
        }

        #[inline(always)]
        fn horizontal_sum(self) -> f32 {
            unsafe { vaddvq_f32(vaddq_f32(self.lo, self.hi)) }
        }
    }
}

/// Eight `f32` as one AVX2 register.
#[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
#[derive(Clone, Copy)]
#[repr(C, align(32))]
pub(crate) struct Simd8 {
    ymm: core::arch::x86_64::__m256,
}

#[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
#[allow(
    clippy::undocumented_unsafe_blocks,
    reason = "every block is one AVX2 intrinsic on values this type owns; the only\
              safety condition any of them has is the target feature, which the\
              module's cfg is"
)]
mod x86 {
    use core::arch::x86_64::{
        _CMP_EQ_OQ, _CMP_LT_OQ, _mm256_add_ps, _mm256_and_ps, _mm256_blendv_ps, _mm256_cmp_ps,
        _mm256_loadu_ps, _mm256_max_ps, _mm256_mul_ps, _mm256_or_ps,
        _mm256_permutevar8x32_ps, _mm256_set1_ps, _mm256_setr_epi32, _mm256_setr_ps,
        _mm256_storeu_ps,
    };

    use super::Simd8;
    use crate::banded::{Lane, Window};

    impl core::ops::Add for Simd8 {
        type Output = Self;
        #[inline(always)]
        fn add(self, rhs: Self) -> Self {
            unsafe { Self { ymm: _mm256_add_ps(self.ymm, rhs.ymm) } }
        }
    }

    impl core::ops::Mul for Simd8 {
        type Output = Self;
        #[inline(always)]
        fn mul(self, rhs: Self) -> Self {
            unsafe { Self { ymm: _mm256_mul_ps(self.ymm, rhs.ymm) } }
        }
    }

    impl Lane for Simd8 {
        const LANES: usize = 8;

        #[inline(always)]
        fn splat(value: f32) -> Self {
            unsafe { Self { ymm: _mm256_set1_ps(value) } }
        }

        #[inline(always)]
        fn load(source: &Window) -> Self {
            unsafe { Self { ymm: _mm256_loadu_ps(source.as_ptr()) } }
        }

        #[inline(always)]
        fn store(self, destination: &mut Window) {
            unsafe { _mm256_storeu_ps(destination.as_mut_ptr(), self.ymm) }
        }

        #[inline(always)]
        fn vmax(self, other: Self) -> Self {
            unsafe { Self { ymm: _mm256_max_ps(self.ymm, other.ymm) } }
        }

        #[inline(always)]
        fn horizontal_max(self) -> f32 {
            let [a, b, c, d, e, f, g, h] = self.to_array();
            a.max(b).max(c.max(d)).max(e.max(f).max(g.max(h)))
        }

        #[inline(always)]
        fn equals(self, other: Self) -> Self {
            unsafe { Self { ymm: _mm256_cmp_ps::<_CMP_EQ_OQ>(self.ymm, other.ymm) } }
        }

        #[inline(always)]
        fn either(self, other: Self) -> Self {
            unsafe { Self { ymm: _mm256_or_ps(self.ymm, other.ymm) } }
        }

        #[inline(always)]
        fn below(self, other: Self) -> Self {
            unsafe { Self { ymm: _mm256_cmp_ps::<_CMP_LT_OQ>(self.ymm, other.ymm) } }
        }

        #[inline(always)]
        fn offsets() -> Self {
            unsafe { Self { ymm: _mm256_setr_ps(0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0) } }
        }

        #[inline(always)]
        fn both(self, other: Self) -> Self {
            unsafe { Self { ymm: _mm256_and_ps(self.ymm, other.ymm) } }
        }

        #[inline(always)]
        fn select(self, if_true: Self, if_false: Self) -> Self {
            unsafe { Self { ymm: _mm256_blendv_ps(if_false.ymm, if_true.ymm, self.ymm) } }
        }

        /// One cross-lane permute plus a blend of the shifted-in scalar:
        /// AVX2 has no whole-register byte shift, so this is the shape the
        /// array permutation `wide` forces the compiler to recognise.
        #[inline(always)]
        fn shift_in(self, first: f32) -> Self {
            unsafe {
                let rotated = _mm256_permutevar8x32_ps(
                    self.ymm,
                    _mm256_setr_epi32(7, 0, 1, 2, 3, 4, 5, 6),
                );
                let head = _mm256_set1_ps(first);
                // Lane 0 from `head`, the rest from the rotation.
                let mask = _mm256_cmp_ps::<_CMP_LT_OQ>(
                    _mm256_setr_ps(0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0),
                    _mm256_set1_ps(1.0),
                );
                Self { ymm: _mm256_blendv_ps(rotated, head, mask) }
            }
        }

        #[inline(always)]
        fn last(self) -> f32 {
            self.to_array()[7]
        }

        #[inline(always)]
        fn horizontal_sum(self) -> f32 {
            self.to_array().iter().sum()
        }
    }

    impl Simd8 {
        #[inline(always)]
        fn to_array(self) -> [f32; 8] {
            let mut out = [0.0f32; 8];
            unsafe { _mm256_storeu_ps(out.as_mut_ptr(), self.ymm) };
            out
        }
    }
}

/// Whether this build has an intrinsics lane at all.
pub(crate) const AVAILABLE: bool = cfg!(any(
    all(target_arch = "aarch64", target_feature = "neon"),
    all(target_arch = "x86_64", target_feature = "avx2")
));

const _: () = assert!(LANE_MAX == 8, "the intrinsics lanes are eight wide");
