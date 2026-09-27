//! What a GPU that flushes subnormal intermediates does to the scores, on the
//! CPU: the shader's Rust transcription with subnormals kept against the same
//! with them flushed, over random pairs from the test strategies.
//!
//! The GPU tests show both measured GPUs return the flushed transcription bit
//! for bit, so this is the GPU's error without needing one.
//!
//! `cargo run -p compair --release --features gpu --example gpu_flush [pairs]`

#[allow(dead_code, reason = "the example uses two of the strategies")]
#[path = "../tests/support/mod.rs"]
mod support;

use compair::{
    Band, StandardEmission,
    gpu::{GpuPairs, Subnormals},
};
use proptest::{
    strategy::{Strategy, ValueTree},
    test_runner::TestRunner,
};

const FLOORS: [f64; 7] = [-10.0, -20.0, -30.0, -40.0, -60.0, -80.0, f64::NEG_INFINITY];

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let total: usize = std::env::args().nth(1).and_then(|n| n.parse().ok()).unwrap_or(200_000);
    let mut runner = TestRunner::deterministic();
    let emission = StandardEmission::default();
    let mut worst = [0.0f64; FLOORS.len()];
    let mut count = [0usize; FLOORS.len()];
    let (mut differ, mut raised) = (0usize, 0.0f64);
    for index in 0..total {
        // Half unconstrained pairs under any width, half reads derived from
        // their haplotype under the default band.
        let (case, band) = if index % 2 == 0 {
            let (case, width) = (support::arbitrary_case(), 2u32..=160)
                .new_tree(&mut runner)
                .map_err(|reason| format!("{reason}"))?
                .current();
            let band = Band::new(width, case.offset)?;
            (case, band)
        } else {
            let case = support::derived_case(6)
                .new_tree(&mut runner)
                .map_err(|reason| format!("{reason}"))?
                .current();
            let band = case.band();
            (case, band)
        };
        let mut pairs = GpuPairs::new();
        let read = pairs.push_read(&case.read, &emission)?;
        let haplotype = pairs.push_haplotype(&case.haplotype, case.read.strand(), &emission)?;
        pairs.push_pair(read, haplotype, band)?;
        let score =
            |subnormals| pairs.emulate(subnormals).first().map_or(f64::NAN, |score| score.get());
        let (kept, flushed) = (score(Subnormals::Kept), score(Subnormals::Flushed));
        if kept.to_bits() != flushed.to_bits() {
            differ += 1;
        }
        if !kept.is_finite() {
            continue;
        }
        if flushed.is_finite() {
            raised = raised.max(flushed - kept);
        }
        let error = if flushed.is_finite() { (kept - flushed).abs() } else { f64::INFINITY };
        for ((floor, worst), count) in FLOORS.iter().zip(&mut worst).zip(&mut count) {
            if kept > *floor {
                *worst = worst.max(error);
                *count += 1;
            }
        }
    }
    println!(
        "{differ}/{total} pairs differ in any bit; flushing raised a score by at most {raised:.2e}"
    );
    for ((floor, worst), count) in FLOORS.iter().zip(worst).zip(count) {
        println!("score > {floor:>6}: {count:>6} pairs, max |Δlog10| {worst:.2e}");
    }
    Ok(())
}
