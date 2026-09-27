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
//! The `strips` arms are the same band through the row-strip traversal.
//!
//! `shadow/64x3` is a variant caller's shadow scoring: 64 reads, each with its
//! own band, against three haplotypes, 192 alignments per iteration. The
//! per-alignment cost is the reported time over 192.

mod fixture;

use compair::{
    Band, Betas, ConversionModel, Haplotype, Pair, Read, StandardEmission, TapsEmission, Workspace,
    align_banded, align_banded_simd, align_full, align_strips, align_strips_simd,
};
use criterion::{Criterion, criterion_group, criterion_main};
use fixture::{Fixture, Shadow, fixture, shadow};
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
    group.bench_function("strips/standard", |b| {
        b.iter(|| align_strips(black_box(haplotype), black_box(&read), &standard, band));
    });
    group.bench_function("strips-simd/standard", |b| {
        b.iter(|| align_strips_simd(black_box(haplotype), black_box(&read), &standard, band));
    });
    group.bench_function("strips/taps", |b| {
        b.iter(|| align_strips(black_box(haplotype), black_box(&read), &taps, band));
    });
    group.bench_function("strips-simd/taps", |b| {
        b.iter(|| align_strips_simd(black_box(haplotype), black_box(&read), &taps, band));
    });
    group.bench_function("strips-simd/taps/workspace", |b| {
        let mut workspace = Workspace::new();
        b.iter(|| workspace.align_strips_simd(black_box(haplotype), black_box(&read), &taps, band));
    });
    group.bench_function("reference/standard", |b| {
        b.iter(|| align_full(black_box(haplotype), black_box(&read), &standard));
    });
    group.finish();

    // One read against `n` candidates, the shape the shadow path has: the
    // per-alignment cost is the reported time over `n`. The batch arm computes
    // eight lanes whatever `n` is, so this is where its fill shows.
    for count in [2usize, 3, 4, 8] {
        let group_haplotypes: Vec<&_> = haplotypes.iter().take(count).collect::<Vec<_>>();
        let mut group = c.benchmark_group(format!("align/150x46/{count}-haplotypes"));
        group.bench_function("strips-simd/taps/workspace", |b| {
            let mut workspace = Workspace::new();
            b.iter(|| {
                group_haplotypes
                    .iter()
                    .map(|haplotype| {
                        workspace
                            .align_strips_simd(black_box(haplotype), black_box(&read), &taps, band)
                            .get()
                    })
                    .sum::<f64>()
            });
        });
        group.bench_function("simd/taps/workspace", |b| {
            let mut workspace = Workspace::new();
            b.iter(|| {
                group_haplotypes
                    .iter()
                    .map(|haplotype| {
                        workspace
                            .align_banded_simd(black_box(haplotype), black_box(&read), &taps, band)
                            .get()
                    })
                    .sum::<f64>()
            });
        });
        group.bench_function("batch/taps/workspace", |b| {
            let mut workspace = Workspace::new();
            let mut out = Vec::with_capacity(count);
            b.iter(|| {
                workspace.align_batch(
                    black_box(&group_haplotypes),
                    black_box(&read),
                    &taps,
                    band,
                    &mut out,
                );
                out.iter().map(|s| s.get()).sum::<f64>()
            });
        });
        // What a caller actually calls: the batch kernel at this fill, the
        // strip kernel below `BATCH_BREAK_EVEN`. At eight haplotypes this
        // should sit on top of the arm above; the arms diverge at a partial
        // batch, which `examples/batchfill.rs` sweeps.
        group.bench_function("candidates/taps/workspace", |b| {
            let mut workspace = Workspace::new();
            let mut out = Vec::with_capacity(count);
            b.iter(|| {
                workspace.align_candidates(
                    black_box(&group_haplotypes),
                    black_box(&read),
                    &taps,
                    band,
                    &mut out,
                );
                out.iter().map(|s| s.get()).sum::<f64>()
            });
        });
        group.finish();
    }

    shadow_locus(c);
}

/// Every read of a locus against every haplotype: the strip kernel one pair at
/// a time (what rastair calls today), the pairs kernel through
/// `align_reads`, and the batch kernel's entry point read by read, which at
/// three haplotypes is the strip kernel again.
fn shadow_locus(c: &mut Criterion) {
    let Shadow { haplotypes, reads, betas } = shadow();
    let taps = TapsEmission::new(ConversionModel::taps_default(), Betas::PerSite(&betas));
    let refs: Vec<&Haplotype> = haplotypes.iter().collect();
    let reads: Vec<(&Read, Band)> = reads.iter().map(|(read, band)| (read, *band)).collect();
    let pairs: Vec<Pair<'_>> = reads
        .iter()
        .flat_map(|&(read, band)| refs.iter().map(move |&h| Pair::new(h, read, band)))
        .collect();

    let mut group = c.benchmark_group(format!("shadow/{}x{}", reads.len(), refs.len()));
    group.bench_function("strips-simd/taps/workspace", |b| {
        let mut workspace = Workspace::new();
        b.iter(|| {
            black_box(&pairs)
                .iter()
                .map(|pair| {
                    workspace.align_strips_simd(pair.haplotype, pair.read, &taps, pair.band).get()
                })
                .sum::<f64>()
        });
    });
    group.bench_function("candidates/taps/workspace", |b| {
        let mut workspace = Workspace::new();
        let mut out = Vec::with_capacity(refs.len());
        b.iter(|| {
            let mut sum = 0.0;
            for &(read, band) in black_box(&reads) {
                workspace.align_candidates(&refs, read, &taps, band, &mut out);
                sum += out.iter().map(|s| s.get()).sum::<f64>();
            }
            sum
        });
    });
    group.bench_function("reads/taps/workspace", |b| {
        let mut workspace = Workspace::new();
        let mut out = Vec::with_capacity(pairs.len());
        b.iter(|| {
            workspace.align_reads(black_box(&refs), black_box(&reads), &taps, &mut out);
            out.iter().map(|s| s.get()).sum::<f64>()
        });
    });
    group.finish();
}

criterion_group!(benches, align);
criterion_main!(benches);
