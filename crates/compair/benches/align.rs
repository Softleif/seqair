//! The C4 cost gate: one 150 bp read against a ~200 bp haplotype through the
//! default band, which the design budgets at ~20 us with 8-wide f32 SIMD and
//! ~60 us scalar.
//!
//! Do not read numbers off a laptop under load; this is meant for the Linux
//! box. The TAPS and standard arms should be indistinguishable: the emission
//! is folded into the plan before the loop runs, so a gap between them means
//! that hoisting has regressed. The `workspace` arms are the shape a caller
//! has -- one read against eight haplotypes, buffers kept between calls --
//! and the difference to the plain arms is the cost of allocating per call.

mod fixture;

use compair::{
    Betas, ConversionModel, StandardEmission, TapsEmission, Workspace, align_banded,
    align_banded_simd, align_full,
};
use criterion::{Criterion, criterion_group, criterion_main};
use fixture::{Fixture, HAPLOTYPES_PER_READ, fixture};
use std::hint::black_box;

fn align(c: &mut Criterion) {
    let Fixture { haplotypes, read, betas, band } = fixture();
    let haplotype = haplotypes.first().expect("the fixture has eight haplotypes");
    let taps = TapsEmission::new(ConversionModel::taps_default(), Betas::PerSite(&betas));
    let standard = StandardEmission::default();

    let mut group = c.benchmark_group("align/150x46");
    group.bench_function("scalar/standard", |b| {
        b.iter(|| align_banded(black_box(haplotype), black_box(&read), &standard, band));
    });
    group.bench_function("simd/standard", |b| {
        b.iter(|| align_banded_simd(black_box(haplotype), black_box(&read), &standard, band));
    });
    group.bench_function("scalar/taps", |b| {
        b.iter(|| align_banded(black_box(haplotype), black_box(&read), &taps, band));
    });
    group.bench_function("simd/taps", |b| {
        b.iter(|| align_banded_simd(black_box(haplotype), black_box(&read), &taps, band));
    });
    group.bench_function("simd/taps/workspace", |b| {
        let mut workspace = Workspace::new();
        b.iter(|| workspace.align_banded_simd(black_box(haplotype), black_box(&read), &taps, band));
    });
    group.bench_function("reference/standard", |b| {
        b.iter(|| align_full(black_box(haplotype), black_box(&read), &standard));
    });
    group.finish();

    let mut group = c.benchmark_group(format!("align/150x46/{HAPLOTYPES_PER_READ}-haplotypes"));
    group.bench_function("simd/taps", |b| {
        b.iter(|| {
            haplotypes
                .iter()
                .map(|haplotype| {
                    align_banded_simd(black_box(haplotype), black_box(&read), &taps, band).get()
                })
                .sum::<f64>()
        });
    });
    group.bench_function("simd/taps/workspace", |b| {
        let mut workspace = Workspace::new();
        b.iter(|| {
            haplotypes
                .iter()
                .map(|haplotype| {
                    workspace
                        .align_banded_simd(black_box(haplotype), black_box(&read), &taps, band)
                        .get()
                })
                .sum::<f64>()
        });
    });
    group.finish();
}

criterion_group!(benches, align);
criterion_main!(benches);
