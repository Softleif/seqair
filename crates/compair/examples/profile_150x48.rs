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
//!
//! `batch-standard`, `batch-taps`, `candidates-standard` and `candidates-taps`
//! take the whole group in one call rather than one haplotype at a time, which
//! is the only way to sample the batch kernel -- the per-haplotype loop above
//! can never fill its lanes from the outside. `pairs-standard` and
//! `pairs-taps` put the same eight alignments through the pairs kernel, so the
//! two kernels race on identical work.
//!
//! `shadow-strips-{standard,taps}`, `shadow-reads-{standard,taps}` and
//! `shadow-pairs-{standard,taps}` score the shadow fixture instead -- 64 reads,
//! each through its own band, against three haplotypes -- one pair at a time
//! through the strip kernel, through `align_reads`, and through `align_pairs`
//! whatever the fill. A pass is 192 alignments, and the rounds are scaled so a
//! run is as many alignments as the other arms' (`rounds * 8 / 192` passes).
#![allow(clippy::print_stdout, reason = "this example exists to print a checksum")]

#[path = "../benches/fixture.rs"]
mod fixture;

use std::hint::black_box;

use compair::{
    Band, Betas, ConversionModel, Emission, Haplotype, Log10Likelihood, Pair, Read,
    StandardEmission, TapsEmission, Workspace,
};

use fixture::{Fixture, HAPLOTYPES_PER_READ, Shadow, fixture, shadow};

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

    if let Some(shadow_kernel) = kernel.strip_prefix("shadow-") {
        let Shadow { haplotypes, reads, betas } = shadow();
        let taps = TapsEmission::new(ConversionModel::taps_default(), Betas::PerSite(&betas));
        let refs: Vec<&Haplotype> = haplotypes.iter().collect();
        let reads: Vec<(&Read, Band)> = reads.iter().map(|(read, band)| (read, *band)).collect();
        let passes = (rounds * HAPLOTYPES_PER_READ / (refs.len() * reads.len())).max(1);
        let (which_kernel, which_emission) =
            shadow_kernel.split_once('-').unwrap_or(("reads", "taps"));
        for _ in 0..passes {
            let mut fresh = Workspace::new();
            let workspace = if reuse { &mut workspace } else { &mut fresh };
            checksum += if which_emission == "standard" {
                shadow_pass(workspace, which_kernel, black_box(&refs), black_box(&reads), &standard)
            } else {
                shadow_pass(workspace, which_kernel, black_box(&refs), black_box(&reads), &taps)
            };
        }
        println!("{which}: {passes} passes, checksum {checksum}");
        return;
    }

    if let Some(group_kernel) = kernel
        .strip_prefix("batch-")
        .map(|e| ("batch", e))
        .or_else(|| kernel.strip_prefix("candidates-").map(|emission| ("candidates", emission)))
        .or_else(|| kernel.strip_prefix("pairs-").map(|emission| ("pairs", emission)))
    {
        let (which_kernel, which_emission) = group_kernel;
        let refs: Vec<&Haplotype> = haplotypes.iter().collect();
        let pairs: Vec<Pair<'_>> = refs.iter().map(|h| Pair::new(h, &read, band)).collect();
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
                ("pairs", "standard") => {
                    workspace.align_pairs(black_box(&pairs), &standard, &mut out);
                }
                ("pairs", _) => workspace.align_pairs(black_box(&pairs), &taps, &mut out),
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
                _ => workspace.align_banded_simd(haplotype, &read, &taps, band),
            };
            checksum += score.get();
        }
    }
    println!("{which}: {rounds} rounds, checksum {checksum}");
}

/// One pass over the shadow fixture through the named kernel.
fn shadow_pass<E: Emission>(
    workspace: &mut Workspace,
    kernel: &str,
    haplotypes: &[&Haplotype],
    reads: &[(&Read, Band)],
    emission: &E,
) -> f64 {
    let mut out: Vec<Log10Likelihood> = Vec::new();
    match kernel {
        "strips" => {
            return reads
                .iter()
                .flat_map(|&(read, band)| haplotypes.iter().map(move |h| (*h, read, band)))
                .map(|(haplotype, read, band)| {
                    workspace.align_strips_simd(haplotype, read, emission, band).get()
                })
                .sum();
        }
        "pairs" => {
            let pairs: Vec<Pair<'_>> = reads
                .iter()
                .flat_map(|&(read, band)| haplotypes.iter().map(move |h| Pair::new(h, read, band)))
                .collect();
            workspace.align_pairs(&pairs, emission, &mut out);
        }
        _ => workspace.align_reads(haplotypes, reads, emission, &mut out),
    }
    out.iter().map(|score| score.get()).sum()
}
