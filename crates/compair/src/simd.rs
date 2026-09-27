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

use fearless_simd::{
    Level, Select, Simd, SimdBase, SimdCombine, SimdFrom, dispatch, f32x4, f32x8, mask32x8,
};
use fearless_simd_macros::simd;

use crate::{
    banded::{Band, LANE_MAX, Lane, LaneMask, Plan, PlanView, Ring, Shape, Window, banded_kernel},
    batch::{BATCH, BatchBuffer, BatchView, batch_kernel},
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
        /// `first, value[0], .., value[6]`, as a `vinsertf128` and a
        /// `vpalignr`: `[first; 4] | value[0..4]` beside `value`, and each
        /// half shifted up by one lane across the pair.
        ///
        /// In `asm!` because LLVM sees through the intrinsics: both are
        /// generic shuffles to it, and it re-lowers the pair to whatever its
        /// throughput model prefers -- a `vpermps` and a blend, or a
        /// `vperm2f128` and two `vshufps`, both of which are longer on the
        /// strip kernel's carried chains.
        #[inline(always)]
        pub(super) fn shift_in(_avx2: Avx2, value: __m256, first: f32) -> __m256 {
            use core::arch::x86_64::_mm256_set1_ps;
            let first = _mm256_set1_ps(first);
            let shifted: __m256;
            // SAFETY: two AVX2 register-to-register instructions, in a
            // function the macro compiles with AVX2 enabled; no memory, no
            // stack, no flags.
            unsafe {
                core::arch::asm!(
                    "vinsertf128 {below}, {first}, {value:x}, 1",
                    "vpalignr {shifted}, {value}, {below}, 12",
                    first = in(ymm_reg) first,
                    value = in(ymm_reg) value,
                    below = out(ymm_reg) _,
                    shifted = lateout(ymm_reg) shifted,
                    options(pure, nomem, nostack, preserves_flags),
                );
            }
            shifted
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
        #[cfg(target_arch = "x86_64")]
        if let Level::Avx2(avx2) = self.simd.level() {
            return Self::simd_from(self.simd, x86::shift_in(avx2, self.into(), first));
        }
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

    /// Four 4 x 4 transposes of the windows' halves, each two rounds of
    /// `zip`s (`unpcklps`/`unpckhps` on x86, `zip1`/`zip2` on NEON), and the
    /// halves of lanes 0-3 and 4-7 put side by side. `f32x8`'s own `zip`s
    /// cross the 128-bit halves, which costs a `vperm2f128` per `zip` on
    /// AVX2; the halves never need to.
    #[inline(always)]
    fn transpose(simd: S, entries: [&Window; LANE_MAX]) -> [Self; LANE_MAX] {
        let [e0, e1, e2, e3, e4, e5, e6, e7] = entries;
        let (l0, h0) = halves(simd, e0);
        let (l1, h1) = halves(simd, e1);
        let (l2, h2) = halves(simd, e2);
        let (l3, h3) = halves(simd, e3);
        let (l4, h4) = halves(simd, e4);
        let (l5, h5) = halves(simd, e5);
        let (l6, h6) = halves(simd, e6);
        let (l7, h7) = halves(simd, e7);
        let [a0, a1, a2, a3] = transpose4(l0, l1, l2, l3);
        let [b0, b1, b2, b3] = transpose4(l4, l5, l6, l7);
        let [c0, c1, c2, c3] = transpose4(h0, h1, h2, h3);
        let [d0, d1, d2, d3] = transpose4(h4, h5, h6, h7);
        [
            a0.combine(b0),
            a1.combine(b1),
            a2.combine(b2),
            a3.combine(b3),
            c0.combine(d0),
            c1.combine(d1),
            c2.combine(d2),
            c3.combine(d3),
        ]
    }
}

/// A window as its two `f32x4` halves.
#[inline(always)]
fn halves<S: Simd>(simd: S, window: &Window) -> (f32x4<S>, f32x4<S>) {
    let [a, b, c, d, e, f, g, h] = *window;
    (f32x4::load_array(simd, [a, b, c, d]), f32x4::load_array(simd, [e, f, g, h]))
}

/// `[a[t], b[t], c[t], d[t]]` for each `t`.
#[inline(always)]
fn transpose4<S: Simd>(a: f32x4<S>, b: f32x4<S>, c: f32x4<S>, d: f32x4<S>) -> [f32x4<S>; 4] {
    let (ac_low, ac_high) = (a.zip_low(c), a.zip_high(c));
    let (bd_low, bd_high) = (b.zip_low(d), b.zip_high(d));
    [
        ac_low.zip_low(bd_low),
        ac_low.zip_high(bd_low),
        ac_high.zip_low(bd_high),
        ac_high.zip_high(bd_high),
    ]
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
    plan: BatchView<'_>,
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
    plan: BatchView<'_>,
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

    /// Eight windows through `transpose` on one level, as arrays.
    #[inline(always)]
    fn transposed<S: Simd>(simd: S, windows: &[[f32; 8]; 8]) -> [[f32; 8]; 8] {
        let [w0, w1, w2, w3, w4, w5, w6, w7] = windows;
        <f32x8<S> as Lane>::transpose(simd, [w0, w1, w2, w3, w4, w5, w6, w7]).map(Into::into)
    }

    proptest! {
        /// Lane `k` of `transpose`'s result `t` is window `k`'s element `t`,
        /// bit for bit, at every level.
        #[test]
        fn transpose_moves_element_t_of_window_k_to_lane_k_of_vector_t(
            windows in proptest::array::uniform8(proptest::array::uniform8(any::<u32>())),
        ) {
            let windows = windows.map(|window| window.map(f32::from_bits));
            for level in levels() {
                let got = dispatch!(level, simd => transposed(simd, &windows));
                for (t, vector) in got.iter().enumerate() {
                    for (k, lane) in vector.iter().enumerate() {
                        let want = windows.get(k).and_then(|window| window.get(t));
                        prop_assert_eq!(
                            Some(lane.to_bits()), want.map(|value| value.to_bits()),
                            "{:?}: vector {} lane {}", level, t, k
                        );
                    }
                }
            }
        }

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
