//! The bench fixture in a loop, so `samply` has something to sample.
//!
//! ```bash
//! cargo build --profile profiling -p compair --example profile_150x48
//! samply record --save-only --unstable-presymbolicate -o p.json --rate 4000 -- \
//!   ./target/profiling/examples/profile_150x48 simd-taps 200000
//! ```
//!
//! The first argument picks the kernel: `scalar-standard`, `scalar-taps`,
//! `simd-standard`, `simd-taps` (the default), `strips-standard`,
//! `strips-taps`, `strips-simd-standard`, `strips-simd-taps`,
//! `strips-intrinsics-standard`, `strips-intrinsics-taps`, or any of those
//! with a `-workspace` suffix to reuse buffers between calls. Every round
//! scores the read against the fixture's eight candidate haplotypes, as a
//! caller does.
//!
//! `batch-standard`, `batch-taps`, `candidates-standard` and `candidates-taps`
//! take the whole group in one call rather than one haplotype at a time, which
//! is the only way to sample the batch kernel -- the per-haplotype loop above
//! can never fill its lanes from the outside.
#![allow(clippy::print_stdout, reason = "this example exists to print a checksum")]

#[path = "../benches/fixture.rs"]
mod fixture;

use std::hint::black_box;

use compair::{
    Betas, ConversionModel, Haplotype, Log10Likelihood, StandardEmission, TapsEmission, Workspace,
};
use fixture::{Fixture, fixture};

fn main() {
    let mut args = std::env::args().skip(1);
    let which = args.next().unwrap_or_else(|| "simd-taps".to_owned());
    let rounds: usize = args.next().and_then(|n| n.parse().ok()).unwrap_or(200_000);

    let Fixture { haplotypes, read, betas, band, .. } = fixture();
    let taps = TapsEmission::new(ConversionModel::taps_default(), Betas::PerSite(&betas));
    let standard = StandardEmission::default();
    let (kernel, reuse) = match which.strip_suffix("-workspace") {
        Some(kernel) => (kernel, true),
        None => (which.as_str(), false),
    };

    let mut workspace = Workspace::new();
    let mut checksum = 0.0f64;

    if let Some(group_kernel) = kernel.strip_prefix("batch-").map(|e| ("batch", e)).or_else(|| {
        kernel.strip_prefix("candidates-").map(|emission| ("candidates", emission))
    }) {
        let (which_kernel, which_emission) = group_kernel;
        let refs: Vec<&Haplotype> = haplotypes.iter().collect();
        let mut out: Vec<Log10Likelihood> = Vec::new();
        for _ in 0..rounds {
            let mut fresh = Workspace::new();
            let workspace = if reuse { &mut workspace } else { &mut fresh };
            let refs = black_box(&refs);
            match (which_kernel, which_emission) {
                ("batch", "standard") => {
                    workspace.align_batch(refs, &read, &standard, band, &mut out);
                }
                ("batch", _) => workspace.align_batch(refs, &read, &taps, band, &mut out),
                (_, "standard") => {
                    workspace.align_candidates(refs, &read, &standard, band, &mut out);
                }
                _ => workspace.align_candidates(refs, &read, &taps, band, &mut out),
            }
            checksum += out.iter().map(|score| score.get()).sum::<f64>();
        }
        println!("{which}: {rounds} rounds, checksum {checksum}");
        return;
    }

    for _ in 0..rounds {
        for haplotype in &haplotypes {
            let haplotype = black_box(haplotype);
            let mut fresh = Workspace::new();
            let workspace = if reuse { &mut workspace } else { &mut fresh };
            let score = match kernel {
                "scalar-standard" => workspace.align_banded(haplotype, &read, &standard, band),
                "scalar-taps" => workspace.align_banded(haplotype, &read, &taps, band),
                "simd-standard" => workspace.align_banded_simd(haplotype, &read, &standard, band),
                "strips-standard" => workspace.align_strips(haplotype, &read, &standard, band),
                "strips-taps" => workspace.align_strips(haplotype, &read, &taps, band),
                "strips-simd-standard" => {
                    workspace.align_strips_simd(haplotype, &read, &standard, band)
                }
                "strips-simd-taps" => workspace.align_strips_simd(haplotype, &read, &taps, band),
                "strips-intrinsics-standard" => {
                    workspace.align_strips_intrinsics(haplotype, &read, &standard, band)
                }
                "strips-intrinsics-taps" => {
                    workspace.align_strips_intrinsics(haplotype, &read, &taps, band)
                }
                _ => workspace.align_banded_simd(haplotype, &read, &taps, band),
            };
            checksum += score.get();
        }
    }
    println!("{which}: {rounds} rounds, checksum {checksum}");
}
