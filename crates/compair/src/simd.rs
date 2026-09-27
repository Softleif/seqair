//! The eight-lane [`Lane`]: `fearless_simd`'s `f32x8<S>`, at whatever level
//! the CPU has.
//!
//! The three banded kernels are generic over the lane, so this module is one
//! `Lane` impl and one `#[simd]` entry point per kernel, each reached through
//! `dispatch!`. `dispatch!` compiles a copy of the kernel per level the target
//! can have -- NEON on aarch64; SSE2, SSE4.2, AVX2 and AVX-512 on x86-64 -- in
//! a `#[target_feature]` context, and picks one per call from a cached CPU
//! check. A default x86-64 build therefore runs AVX2 where the CPU has it,
//! with no `unsafe` and no build flags.
//!
//! **The kernels and every helper they call are `#[inline(always)]`, and no
//! vector op sits in a closure.** `#[simd]` enables the level's target
//! features for the body of the function it annotates and nothing else: a
//! helper that LLVM leaves out of line, or a closure (which does not inherit
//! target features), is compiled for the baseline, and on x86-64 every vector
//! op in it becomes a call. That was measured, not guessed: with the strip
//! kernel merely `#[inline]`, its AVX2 instance was an out-of-line SSE2
//! function making a call per multiply, add and splat. It cannot happen on
//! aarch64, where NEON is the baseline, so only an x86-64 build shows it.
//!
//! Bit-parity with the scalar lane holds at every level for the reason it
//! holds at all: every lane operation is the same IEEE operation on the same
//! values in the same order, and there is no fused multiply-add on [`Lane`].

use fearless_simd::{Level, Select, Simd, SimdBase, SimdFrom, dispatch, f32x8, mask32x8};
use fearless_simd_macros::simd;

use crate::{
    banded::{Band, Lane, LaneMask, Plan, PlanView, Ring, Shape, Window, banded_kernel},
    batch::{BATCH, BatchBuffer, BatchPlan, batch_kernel},
    pairs::{PAIRS, PairsBuffer, PairsPlan, pairs_kernel},
    strips::{RowBuffer, strip_kernel},
    types::Log10Likelihood,
};

/// The one escape hatch: `masked` and `masked_out` as a plain `and` and
/// `andnot` on the AVX2 level, each a safe `fearless_simd::kernel!`.
///
/// The portable spelling, `mask.select(value, zero)`, is a `vblendvps`, which
/// reads only each lane's sign bit. LLVM folds the zero operand into an `and`
/// but cannot see through `vcmpps` that the mask is already whole, so it
/// sign-extends it first with a `vpcmpgtb` -- one per site, four per strip
/// step, five per batch step. On the 3950X that was 11% more cycles on the
/// strip kernel and 11% on the batch kernel (notes §11). NEON needs no hatch:
/// the portable spelling already lowers to `and`/`bic` there. Nor does
/// AVX-512, whose masks are `k` registers, where a `select` against zero is a
/// zero-masked move and already the best spelling.
///
/// Other hatches were measured and dropped: the hand-written lane's `vpermps`
/// lane shift, and a `vperm2f128` + `vpalignr` one, each fewer instructions
/// than the portable `shift_elements_right` and each 1-2% slower with loops
/// aligned (§11).
#[cfg(target_arch = "x86_64")]
mod x86 {
    use core::arch::x86_64::{__m256, __m256i};

    fearless_simd::kernel!(
        /// `value` where `mask` is set, as one `vandps`.
        #[inline(always)]
        pub(super) fn and(_avx2: Avx2, mask: __m256i, value: __m256) -> __m256 {
            use core::arch::x86_64::{_mm256_and_ps, _mm256_castsi256_ps};
            _mm256_and_ps(_mm256_castsi256_ps(mask), value)
        }
    );

    fearless_simd::kernel!(
        /// `value` where `mask` is clear, as one `vandnps`.
        #[inline(always)]
        pub(super) fn and_not(_avx2: Avx2, mask: __m256i, value: __m256) -> __m256 {
            use core::arch::x86_64::{_mm256_andnot_ps, _mm256_castsi256_ps};
            _mm256_andnot_ps(_mm256_castsi256_ps(mask), value)
        }
    );
}

impl<S: Simd> LaneMask for mask32x8<S> {
    #[inline(always)]
    fn either(self, other: Self) -> Self {
        self | other
    }

    #[inline(always)]
    fn both(self, other: Self) -> Self {
        self & other
    }
}

impl<S: Simd> Lane for f32x8<S> {
    const LANES: usize = 8;
    type Token = S;
    type Mask = mask32x8<S>;

    #[inline(always)]
    fn token(self) -> S {
        self.simd
    }

    #[inline(always)]
    fn splat(token: S, value: f32) -> Self {
        <Self as SimdBase<S>>::splat(token, value)
    }

    #[inline(always)]
    fn load(token: S, source: &Window) -> Self {
        <Self as SimdBase<S>>::load_array_ref(token, source)
    }

    #[inline(always)]
    fn store(self, destination: &mut Window) {
        self.store_array(destination);
    }

    #[inline(always)]
    fn vmax(self, other: Self) -> Self {
        // Not `max_precise`: nothing the kernels feed a maximum is a NaN, and
        // the precise form is a compare and a blend more on x86.
        self.max(other)
    }

    #[inline(always)]
    fn horizontal_max(self) -> f32 {
        self.reduce_max()
    }

    #[inline(always)]
    fn equals(self, other: Self) -> mask32x8<S> {
        self.simd_eq(other)
    }

    #[inline(always)]
    fn below(self, other: Self) -> mask32x8<S> {
        self.simd_lt(other)
    }

    #[inline(always)]
    fn offsets(token: S) -> Self {
        Self::simd_from(token, [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
    }

    #[inline(always)]
    fn select(mask: mask32x8<S>, if_true: Self, if_false: Self) -> Self {
        mask.select(if_true, if_false)
    }

    #[inline(always)]
    fn masked(mask: mask32x8<S>, value: Self) -> Self {
        #[cfg(target_arch = "x86_64")]
        if let Level::Avx2(avx2) = value.simd.level() {
            return Self::simd_from(value.simd, x86::and(avx2, mask.into(), value.into()));
        }
        mask.select(value, <Self as SimdBase<S>>::splat(value.simd, 0.0))
    }

    #[inline(always)]
    fn masked_out(mask: mask32x8<S>, value: Self) -> Self {
        #[cfg(target_arch = "x86_64")]
        if let Level::Avx2(avx2) = value.simd.level() {
            return Self::simd_from(value.simd, x86::and_not(avx2, mask.into(), value.into()));
        }
        mask.select(<Self as SimdBase<S>>::splat(value.simd, 0.0), value)
    }

    #[inline(always)]
    fn shift_in(self, first: f32) -> Self {
        self.shift_elements_right::<1>(first)
    }

    #[inline(always)]
    fn last(self) -> f32 {
        let [.., last] = self.to_array();
        last
    }

    #[inline(always)]
    fn horizontal_sum(self) -> f32 {
        self.reduce_sum()
    }
}

#[simd]
fn strip_kernel_simd<S: Simd>(
    simd: S,
    plan: PlanView<'_>,
    rows: &mut RowBuffer,
    shape: Shape,
    band: Band,
) -> Log10Likelihood {
    strip_kernel::<f32x8<S>>(simd, plan, rows, shape, band)
}

#[simd]
fn batch_kernel_simd<S: Simd>(
    simd: S,
    plan: &BatchPlan,
    buffer: &mut BatchBuffer,
    read_len: usize,
    band: Band,
) -> [Log10Likelihood; BATCH] {
    batch_kernel::<f32x8<S>>(simd, plan, buffer, read_len, band)
}

#[simd]
fn pairs_kernel_simd<S: Simd>(
    simd: S,
    plan: &PairsPlan,
    buffer: &mut PairsBuffer,
) -> [Log10Likelihood; PAIRS] {
    pairs_kernel::<f32x8<S>>(simd, plan, buffer)
}

#[simd]
fn banded_kernel_simd<S: Simd>(
    simd: S,
    plan: &Plan,
    ring: &mut Ring,
    shape: Shape,
    band: Band,
) -> Log10Likelihood {
    banded_kernel::<f32x8<S>>(simd, plan, ring, shape, band)
}

/// The strip kernel at `level`.
pub(crate) fn strip_kernel_at(
    level: Level,
    plan: PlanView<'_>,
    rows: &mut RowBuffer,
    shape: Shape,
    band: Band,
) -> Log10Likelihood {
    dispatch!(level, simd => strip_kernel_simd(simd, plan, rows, shape, band))
}

/// The batch kernel at `level`.
pub(crate) fn batch_kernel_at(
    level: Level,
    plan: &BatchPlan,
    buffer: &mut BatchBuffer,
    read_len: usize,
    band: Band,
) -> [Log10Likelihood; BATCH] {
    dispatch!(level, simd => batch_kernel_simd(simd, plan, buffer, read_len, band))
}

/// The pairs kernel at `level`.
pub(crate) fn pairs_kernel_at(
    level: Level,
    plan: &PairsPlan,
    buffer: &mut PairsBuffer,
) -> [Log10Likelihood; PAIRS] {
    dispatch!(level, simd => pairs_kernel_simd(simd, plan, buffer))
}

/// The diagonal kernel at `level`.
pub(crate) fn banded_kernel_at(
    level: Level,
    plan: &Plan,
    ring: &mut Ring,
    shape: Shape,
    band: Band,
) -> Log10Likelihood {
    dispatch!(level, simd => banded_kernel_simd(simd, plan, ring, shape, band))
}

/// The level `dispatch!` picks on this CPU, by name.
#[must_use]
pub fn simd_level() -> &'static str {
    let level = Level::new();
    if level.is_fallback() {
        return "fallback";
    }
    #[cfg(target_arch = "aarch64")]
    if level.as_neon().is_some() {
        return "neon";
    }
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if level.as_avx512().is_some() {
            return "avx512";
        }
        if level.as_avx2().is_some() {
            return "avx2";
        }
        if level.as_sse4_2().is_some() {
            return "sse4.2";
        }
        if level.as_sse2().is_some() {
            return "sse2";
        }
    }
    "other"
}

#[cfg(test)]
mod tests {
    use fearless_simd::{Level, Simd, SimdFrom, dispatch, f32x8, mask32x8};
    use proptest::prelude::*;

    use crate::banded::{Lane, LaneMask};

    /// Every level this CPU can run, the scalar `Fallback` included (the
    /// dev-dependency compiles it in on every target).
    fn levels() -> Vec<Level> {
        let best = Level::new();
        #[allow(unused_mut, reason = "only x86 has levels below the best one")]
        let mut levels = vec![Level::fallback(), best];
        #[cfg(target_arch = "x86_64")]
        {
            levels.extend(best.as_sse2().map(Level::Sse2));
            levels.extend(best.as_sse4_2().map(Level::Sse4_2));
            levels.extend(best.as_avx2().map(Level::Avx2));
        }
        levels
    }

    /// A lane mask with lane `l` set where `bits[l]` is.
    #[inline(always)]
    fn mask_of<S: Simd>(simd: S, bits: [bool; 8]) -> mask32x8<S> {
        f32x8::simd_from(simd, bits.map(|set| if set { 1.0 } else { 0.0 }))
            .equals(<f32x8<S> as Lane>::splat(simd, 1.0))
    }

    /// What the lane operations the strip kernel adds give on one level:
    /// `shift_in`, `last`, a mask `both` another as bits, and the sum of one
    /// live lane.
    #[inline(always)]
    fn lane_ops<S: Simd>(
        simd: S,
        values: [f32; 8],
        first: f32,
        left: [bool; 8],
        right: [bool; 8],
        single: [f32; 8],
    ) -> ([f32; 8], f32, [u32; 8], f32) {
        let vector = f32x8::simd_from(simd, values);
        let shifted: [f32; 8] = vector.shift_in(first).into();
        let both = mask_of(simd, left).both(mask_of(simd, right));
        let ones = <f32x8<S> as Lane>::splat(simd, f32::from_bits(u32::MAX));
        let anded: [f32; 8] = <f32x8<S> as Lane>::masked(both, ones).into();
        let sum = f32x8::simd_from(simd, single).horizontal_sum();
        (shifted, vector.last(), anded.map(f32::to_bits), sum)
    }

    proptest! {
        /// The lane operations with no scalar counterpart in the kernels,
        /// against their definitions, at every level.
        #[test]
        fn lane_ops_do_what_the_trait_says_at_every_level(
            values in proptest::array::uniform8(-1e30f32..1e30),
            first in -1e30f32..1e30,
            left in proptest::array::uniform8(any::<bool>()),
            right in proptest::array::uniform8(any::<bool>()),
            value in -1e30f32..1e30,
            lane in 0usize..8,
        ) {
            let mut single = [0.0f32; 8];
            *single.get_mut(lane).expect("lane in 0..8") = value;
            let [v0, v1, v2, v3, v4, v5, v6, v7] = values;
            for level in levels() {
                let (shifted, last, anded, sum) = dispatch!(
                    level, simd => lane_ops(simd, values, first, left, right, single)
                );
                prop_assert_eq!(shifted, [first, v0, v1, v2, v3, v4, v5, v6], "{:?}", level);
                prop_assert_eq!(last.to_bits(), v7.to_bits(), "{:?}", level);
                for index in 0..8 {
                    let want = if left[index] && right[index] { u32::MAX } else { 0 };
                    prop_assert_eq!(anded[index], want, "{:?} lane {}", level, index);
                }
                prop_assert_eq!(sum.to_bits(), value.to_bits(), "{:?}", level);
            }
        }
    }
}
