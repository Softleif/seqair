//! Where the pairs kernel starts beating the strip kernel, as a function of how
//! many pairs a group holds -- the number `PAIRS_BREAK_EVEN` has to be.
//!
//! The pairs come from the shadow fixture in the order `align_reads` packs
//! them, read-major, so a group of `n` is `n` different alignments with up to
//! `n` different reads, strands and offsets. Both sides race through the same
//! SIMD level, as in `examples/batchfill.rs`, so the crossover is the cost of
//! an empty lane and not the gap between two lanes.
//!
//! `cargo run -p compair --release --example pairsfill`
#![allow(clippy::print_stdout, reason = "this example exists to print a table")]

#[path = "../benches/fixture.rs"]
mod fixture;

use std::time::Instant;

use compair::{Betas, ConversionModel, Emission, Log10Likelihood, Pair, TapsEmission, Workspace};

use fixture::{Shadow, shadow};

fn main() {
    let Shadow { haplotypes, reads, betas } = shadow();
    let taps = TapsEmission::new(ConversionModel::taps_default(), Betas::PerSite(&betas));
    let all: Vec<Pair<'_>> = reads
        .iter()
        .flat_map(|(read, band)| haplotypes.iter().map(move |h| Pair::new(h, read, *band)))
        .collect();
    println!("SIMD level: {}", compair::simd_level());
    println!("{:>5} {:>12} {:>12} {:>9}", "pairs", "strips us/aln", "pairs us/aln", "pairs/s");
    for n in [1usize, 2, 3, 4, 5, 6, 7, 8, 9, 12, 16, 24] {
        let group = all.get(..n).unwrap_or_default();
        let (strips, pairs) = race(group, &taps);
        #[allow(clippy::cast_precision_loss, reason = "small counts")]
        let per = n as f64;
        println!("{n:>5} {:>12.3} {:>12.3} {:>8.2}x", strips / per, pairs / per, strips / pairs);
    }
}

/// Microseconds for one pass of each kernel over the group, best of many.
fn race<E: Emission>(group: &[Pair<'_>], emission: &E) -> (f64, f64) {
    let mut workspace = Workspace::new();
    let mut out: Vec<Log10Likelihood> = Vec::new();
    let strips = time(|| {
        group
            .iter()
            .map(|pair| {
                workspace.align_strips_simd(pair.haplotype, pair.read, emission, pair.band).get()
            })
            .sum()
    });
    let pairs = time(|| {
        workspace.align_pairs(group, emission, &mut out);
        out.iter().map(|score| score.get()).sum()
    });
    (strips, pairs)
}

fn time(mut body: impl FnMut() -> f64) -> f64 {
    for _ in 0..200 {
        std::hint::black_box(body());
    }
    let mut best = f64::INFINITY;
    for _ in 0..400 {
        let start = Instant::now();
        let value = body();
        let elapsed = start.elapsed().as_secs_f64() * 1e6;
        std::hint::black_box(value);
        best = best.min(elapsed);
    }
    best
}
