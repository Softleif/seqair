//! What a GPU that flushes subnormal intermediates does to the scores, on the
//! CPU: the shader's Rust transcription with subnormals kept against the same
//! with them flushed, over pairs from the tests' generators.
//!
//! The GPU tests show both measured GPUs return the flushed transcription bit
//! for bit, so this is the GPU's error without needing one.
//!
//! The pairs are hegel's, derandomized so a run is reproducible. Hegel leans
//! towards the edges of each range, so the counts below are over that mix
//! and not over uniformly random pairs; the worst cases are what it is for.
//!
//! `cargo run -p compair --release --features gpu --example gpu_flush [pairs]`

#[allow(dead_code, reason = "the example uses two of the generators")]
#[path = "../tests/support/mod.rs"]
mod support;

use compair::{
    Band, StandardEmission,
    gpu::{GpuPairs, Subnormals},
};
use hegel::TestCase;
use hegel::generators as gs;

const FLOORS: [f64; 7] = [-10.0, -20.0, -30.0, -40.0, -60.0, -80.0, f64::NEG_INFINITY];

#[derive(Default)]
struct Survey {
    pairs: usize,
    differ: usize,
    raised: f64,
    worst: [f64; FLOORS.len()],
    count: [usize; FLOORS.len()],
}

/// One pair: half unconstrained pairs under any width, half reads derived
/// from their haplotype under the default band.
#[hegel::standalone_function(test_cases = total, derandomize = true, database = None)]
fn survey(tc: TestCase, total: u64, survey: &mut Survey) {
    let (case, band) = if tc.draw(gs::booleans()) {
        let case = tc.draw(support::arbitrary_case());
        let width = tc.draw(gs::integers::<u32>().min_value(2).max_value(160));
        let band = Band::new(width, case.offset).expect("every width in 2..=160 is legal");
        (case, band)
    } else {
        let case = tc.draw(support::derived_case(6));
        let band = case.band();
        (case, band)
    };
    let emission = StandardEmission::default();
    let mut pairs = GpuPairs::new();
    let read = pairs.push_read(&case.read, &emission).expect("the read pushes");
    let haplotype = pairs
        .push_haplotype(&case.haplotype, case.read.strand(), &emission)
        .expect("the haplotype pushes");
    pairs.push_pair(read, haplotype, band).expect("the pair pushes");
    let score =
        |subnormals| pairs.emulate(subnormals).first().map_or(f64::NAN, |score| score.get());
    let (kept, flushed) = (score(Subnormals::Kept), score(Subnormals::Flushed));
    survey.pairs += 1;
    if kept.to_bits() != flushed.to_bits() {
        survey.differ += 1;
    }
    if !kept.is_finite() {
        return;
    }
    if flushed.is_finite() {
        survey.raised = survey.raised.max(flushed - kept);
    }
    let error = if flushed.is_finite() { (kept - flushed).abs() } else { f64::INFINITY };
    for ((floor, worst), count) in FLOORS.iter().zip(&mut survey.worst).zip(&mut survey.count) {
        if kept > *floor {
            *worst = worst.max(error);
            *count += 1;
        }
    }
}

fn main() {
    let total: u64 = std::env::args().nth(1).and_then(|n| n.parse().ok()).unwrap_or(200_000);
    let mut result = Survey::default();
    survey(total, &mut result);
    let Survey { pairs, differ, raised, worst, count } = result;
    println!(
        "{differ}/{pairs} pairs differ in any bit; flushing raised a score by at most {raised:.2e}"
    );
    for ((floor, worst), count) in FLOORS.iter().zip(worst).zip(count) {
        println!("score > {floor:>6}: {count:>6} pairs, max |Δlog10| {worst:.2e}");
    }
}
