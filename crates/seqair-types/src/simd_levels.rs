//! The SIMD levels a kernel test runs at.

use fearless_simd::Level;

/// Every SIMD level this CPU can run, the scalar fallback included, so each
/// kernel is tested at all of them — not only at the one `Level::new()` picks.
// r[verify io.simd_portable]
pub(crate) fn levels() -> Vec<Level> {
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
