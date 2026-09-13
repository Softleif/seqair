//! compair against gkl on the "10s" benchmark, which is the only PairHMM
//! benchmark that is both small enough to ship and published enough to compare
//! against.
//!
//! **Every arm scores all 3,550 pairs per iteration**, so criterion's reported
//! time per iteration is directly the number the papers quote as a total. For
//! reference, on that dataset:
//!
//! | source | machine | time |
//! |---|---|---|
//! | GKL AVX-512 FP32 | Xeon Gold 6438 | 15.39 ms |
//! | GKL AVX-512 FP32 | EPYC Zen 4 | 23.95 ms |
//! | gpuPairHMM | H100 | 0.20 ms |
//! | gpuPairHMM | V100 | 0.37 ms |
//!
//! (gpuPairHMM arXiv:2411.11547 Table 2, re-reported by Endeavor arXiv:2606.25738.)
//!
//! Throughput is set to the **full matrix** cell count, 62,380,634, for every
//! arm including the banded ones. Criterion then prints `Kelem/s`, which is
//! CUPS on the convention both papers use — cells of the problem, not cells
//! anyone chose to visit. The banded arms visit far fewer than that and their
//! `Kelem/s` is correspondingly flattering as a *kernel* measure; that is the
//! honest way round, because what a caller buys is a score per pair, and the
//! band is how we get it cheaply rather than a smaller problem.
//! `tests/tenspeed.rs` prints the visited-cell count if the kernel rate is what
//! is wanted.
//!
//! Two comparisons are not like-for-like and should not be read as one:
//!
//! - gkl computes the **whole matrix**; compair's banded arms compute a
//!   46-column band. The `reference` arm is compair's full matrix, and that is
//!   the one to put next to gkl's for work done.
//! - gkl's `f32x8` is one read-haplotype pair per register swept in row
//!   strips, which is the same traversal as compair's `strips` arms. That pair
//!   is the closest thing to a controlled comparison in this file.
//!
//! **Both gkl implementations, in two runs.** gkl's `c-avx` feature *replaces*
//! its Rust AVX kernels with Intel's original C++ rather than adding them
//! beside them -- `forward_f32x8` is one `cfg_if` over the two, and the Rust
//! AVX module is `#[cfg(not(feature = "c-avx"))]`. So they cannot share a
//! binary, and the comparison is two runs over the same vendored C++ source:
//!
//! ```text
//! cargo bench -p compair --bench tenspeed
//! cargo bench -p compair --bench tenspeed --features gkl-cpp
//! ```
//!
//! The vector arms name themselves `gkl-rust/...` or `gkl-cpp/...` so the runs
//! do not overwrite each other, and the compair arms are identical in both, so
//! they double as a check that the machine was equally quiet for each. The
//! scalar arms stay `gkl-rust/...` in both runs, because `forward_f32x1` and
//! `forward_f64x1` are pure Rust and not feature-gated.
//!
//! Do not read numbers off a laptop, and do not read the gkl arms off Apple
//! Silicon at all: they compile out. This wants the Linux box.

use compair::{Band, StandardEmission, Workspace, align_full};
use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use std::hint::black_box;

#[path = "../tests/support/tenspeed.rs"]
mod tenspeed;

fn tenspeed(c: &mut Criterion) {
    let pairs = tenspeed::pairs().expect("10s.in and 10s.out parse");
    let matrix = tenspeed::matrix_cells(&pairs);
    let banded = tenspeed::band_cells(&pairs, |pair| Band::anchored(pair.offset));
    #[allow(clippy::cast_precision_loss, reason = "reporting a ratio of cell counts")]
    let ratio = banded as f64 / matrix as f64;
    println!(
        "10s: {} pairs, {matrix} matrix cells, {banded} band cells ({:.1}% of the matrix)",
        pairs.len(),
        ratio * 100.0
    );

    let standard = StandardEmission::default();
    let mut group = c.benchmark_group("10s");
    group.throughput(Throughput::Elements(matrix));
    group.sample_size(20);

    group.bench_function("compair/strips-simd/banded", |b| {
        let mut workspace = Workspace::new();
        b.iter(|| {
            pairs
                .iter()
                .map(|pair| {
                    workspace
                        .align_strips_simd(
                            black_box(&pair.haplotype),
                            black_box(&pair.read),
                            &standard,
                            Band::anchored(pair.offset),
                        )
                        .get()
                })
                .sum::<f64>()
        });
    });
    #[cfg(feature = "intrinsics")]
    {
        println!("intrinsics lane available: {}", compair::intrinsics_lane_available());
        group.bench_function("compair/strips-intrinsics/banded", |b| {
            let mut workspace = Workspace::new();
            b.iter(|| {
                pairs
                    .iter()
                    .map(|pair| {
                        workspace
                            .align_strips_intrinsics(
                                black_box(&pair.haplotype),
                                black_box(&pair.read),
                                &standard,
                                Band::anchored(pair.offset),
                            )
                            .get()
                    })
                    .sum::<f64>()
            });
        });
    }
    group.bench_function("compair/banded-simd/banded", |b| {
        let mut workspace = Workspace::new();
        b.iter(|| {
            pairs
                .iter()
                .map(|pair| {
                    workspace
                        .align_banded_simd(
                            black_box(&pair.haplotype),
                            black_box(&pair.read),
                            &standard,
                            Band::anchored(pair.offset),
                        )
                        .get()
                })
                .sum::<f64>()
        });
    });
    group.bench_function("compair/reference/f64-full", |b| {
        b.iter(|| {
            pairs
                .iter()
                .map(|pair| {
                    align_full(black_box(&pair.haplotype), black_box(&pair.read), &standard).get()
                })
                .sum::<f64>()
        });
    });

    #[cfg(target_arch = "x86_64")]
    {
        // gkl's `c-avx` feature replaces its Rust AVX kernels with Intel's
        // original C++ rather than adding them, so only one of the two is in
        // this binary. The vector arms say which; the scalar arms are always
        // gkl's Rust, since `forward_f32x1`/`f64x1` are not feature-gated.
        let vector = if cfg!(feature = "gkl-cpp") { "gkl-cpp" } else { "gkl-rust" };
        println!(
            "gkl vector kernels: {vector}. available: f32x8 {}, f64x4 {}, f32x16 {}, f64x8 {}",
            gkl::pairhmm::forward_f32x8().is_some(),
            gkl::pairhmm::forward_f64x4().is_some(),
            gkl::pairhmm::forward_f32x16().is_some(),
            gkl::pairhmm::forward_f64x8().is_some(),
        );

        macro_rules! gkl_arm {
            ($name:expr, $forward:expr) => {
                let forward = $forward;
                group.bench_function($name, |b| {
                    b.iter(|| {
                        pairs
                            .iter()
                            .map(|pair| {
                                let raw = &pair.raw;
                                f64::from(forward(
                                    black_box(&raw.haplotype),
                                    black_box(&raw.read),
                                    &raw.base_quals,
                                    &raw.insertion_quals,
                                    &raw.deletion_quals,
                                    &raw.gap_quals,
                                ))
                            })
                            .sum::<f64>()
                    });
                });
            };
        }

        if let Some(forward) = gkl::pairhmm::forward_f32x8() {
            gkl_arm!(format!("{vector}/f32x8/full"), forward);
        }
        if let Some(forward) = gkl::pairhmm::forward_f32x16() {
            gkl_arm!(format!("{vector}/f32x16/full"), forward);
        }
        if let Some(forward) = gkl::pairhmm::forward_f64x4() {
            gkl_arm!(format!("{vector}/f64x4/full"), forward);
        }
        // Always gkl's Rust, whatever the feature says: these are not gated.
        gkl_arm!("gkl-rust/f32x1/full", gkl::pairhmm::forward_f32x1());
        gkl_arm!("gkl-rust/f64x1/full", gkl::pairhmm::forward_f64x1());
        // What GATK actually runs: f32, redone in f64 when the score is low.
        gkl_arm!(format!("{vector}/adaptive/full"), gkl::pairhmm::forward);
    }

    group.finish();
}

criterion_group!(benches, tenspeed);
criterion_main!(benches);
