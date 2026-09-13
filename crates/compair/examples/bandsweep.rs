//! What a band width buys and costs, on 3,550 real read-haplotype pairs.
//!
//! `Band::DEFAULT_WIDTH` is currently justified by how it fills SIMD lanes, not
//! by what it loses on real data. This measures the second half of that: for
//! each width, how many of the 10s benchmark's pairs keep the score
//! `align_full` gives, what the ones that don't give up, and what share of the
//! matrix the width visits. Then it names every pair the default width loses,
//! because on this dataset the answer is not "the tail of a distribution".
//!
//! `cargo run -p compair --release --example bandsweep`

use compair::{Band, StandardEmission, align_full, align_strips_simd};

#[path = "../tests/support/tenspeed.rs"]
mod tenspeed;

/// A loss below this is rounding between an `f32` band and an `f64` matrix.
const KEPT: f64 = 1e-3;

fn main() {
    let pairs = tenspeed::pairs().expect("10s.in parses");
    let standard = StandardEmission::default();
    let matrix = tenspeed::matrix_cells(&pairs);
    let full: Vec<f64> =
        pairs.iter().map(|pair| align_full(&pair.haplotype, &pair.read, &standard).get()).collect();

    println!("10s: {} pairs, {matrix} matrix cells\n", pairs.len());
    println!(
        "{:>6} {:>6} {:>6} {:>12} {:>12} {:>8}",
        "width", "kept", "lost", "worst loss", "mean loss", "cells"
    );
    for width in [14u32, 30, 46, 62, 94, 126, 190, 254] {
        let band = |pair: &tenspeed::Pair| {
            Band::new(width, pair.offset).unwrap_or(Band::anchored(pair.offset))
        };
        let (mut kept, mut worst, mut total) = (0usize, 0.0f64, 0.0f64);
        for (pair, reference) in pairs.iter().zip(&full) {
            let loss = reference
                - align_strips_simd(&pair.haplotype, &pair.read, &standard, band(pair)).get();
            if loss < KEPT {
                kept += 1;
            } else {
                worst = worst.max(loss);
                total += loss;
            }
        }
        let lost = pairs.len() - kept;
        let cells = tenspeed::band_cells(&pairs, band);
        #[allow(clippy::cast_precision_loss, reason = "ratios over a few thousand pairs")]
        let (mean, share) = (
            if lost == 0 { 0.0 } else { total / lost as f64 },
            100.0 * cells as f64 / matrix as f64,
        );
        println!("{width:>6} {kept:>6} {lost:>6} {worst:>12.3} {mean:>12.3} {share:>7.1}%");
    }

    println!("\npairs the default width ({}) loses:", Band::DEFAULT_WIDTH);
    println!(
        "{:>6} {:>6} {:>6} {:>8} {:>10} {:>10} {:>9}",
        "pair", "read", "hap", "offset", "full", "banded", "loss"
    );
    for ((index, pair), reference) in pairs.iter().enumerate().zip(&full) {
        let band = Band::anchored(pair.offset);
        let got = align_strips_simd(&pair.haplotype, &pair.read, &standard, band).get();
        if reference - got >= KEPT {
            println!(
                "{index:>6} {:>6} {:>6} {:>8} {reference:>10.3} {got:>10.3} {:>9.3}",
                pair.read.len(),
                pair.haplotype.len(),
                pair.offset,
                reference - got
            );
        }
    }
}
