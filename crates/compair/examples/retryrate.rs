//! How often GATK's production PairHMM has to redo a pair in `f64`.
//!
//! gkl's `forward` runs the `f32` kernel and, if the score falls below
//! `-28 - log10(2^120)`, throws it away and recomputes the whole pair in `f64`.
//! That is the cost compair's lazy power-of-two rescale is designed never to
//! pay. This counts how many of the 3,550 pairs would trip it.
//!
//! `cargo run -p compair --release --example retryrate`

use compair::{StandardEmission, align_full};

#[path = "../tests/support/tenspeed.rs"]
mod tenspeed;

fn main() {
    // gkl: `LOG10_INITIAL_CONSTANT_32 = (120f32).exp2().log10()`, and the f32
    // result is kept only when `result > -28 - that`.
    let threshold = -28.0 - f64::from(120.0f32.exp2().log10());
    let pairs = tenspeed::pairs().expect("10s.in parses");
    let standard = StandardEmission::default();
    let scores: Vec<f64> =
        pairs.iter().map(|pair| align_full(&pair.haplotype, &pair.read, &standard).get()).collect();
    let retries = scores.iter().filter(|score| **score <= threshold).count();
    let mut sorted = scores.clone();
    sorted.sort_by(f64::total_cmp);
    #[allow(clippy::cast_precision_loss, reason = "reporting a percentage")]
    let percent = 100.0 * retries as f64 / pairs.len() as f64;
    println!("gkl f32 -> f64 retry threshold: {threshold:.3} log10");
    println!("pairs below it: {retries} of {} ({percent:.2}%)", pairs.len());
    println!(
        "score distribution: min {:.2}, 1st pct {:.2}, median {:.2}, max {:.2}",
        sorted.first().copied().unwrap_or(f64::NAN),
        sorted.get(pairs.len() / 100).copied().unwrap_or(f64::NAN),
        sorted.get(pairs.len() / 2).copied().unwrap_or(f64::NAN),
        sorted.last().copied().unwrap_or(f64::NAN),
    );
}
