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
//! `strips-taps`, `strips-simd-standard`, `strips-simd-taps`, or any of those
//! with a `-workspace` suffix to reuse buffers between calls. Every round
//! scores the read against the fixture's eight candidate haplotypes, as a
//! caller does.
#![allow(clippy::print_stdout, reason = "this example exists to print a checksum")]

#[path = "../benches/fixture.rs"]
mod fixture;

use std::hint::black_box;

use compair::{Betas, ConversionModel, StandardEmission, TapsEmission, Workspace};
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
