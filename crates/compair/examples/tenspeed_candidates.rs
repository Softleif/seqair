//! The 10s benchmark through `align_candidates` per read, and through one
//! `Workspace::candidates` per group, in a loop, for interleaved A/B runs.
//!
//! ```bash
//! cargo run -p compair --release --example tenspeed_candidates -- prepared 200
//! ```
//!
//! The first argument is `candidates` (the default), `prepared`, or
//! `align-reads` (every read of a group in one `Workspace::align_reads`), the second
//! the number of rounds over the whole dataset. Both score the same 3,550
//! pairs with `StandardEmission` through the default band and print the same
//! checksum: each score is added on its own, in the same order.
#![allow(clippy::print_stdout, reason = "this example exists to print a checksum")]

#[path = "../tests/support/tenspeed.rs"]
mod tenspeed;

use std::{hint::black_box, time::Instant};

use compair::{Band, Haplotype, Log10Likelihood, StandardEmission, Workspace};

fn main() {
    let mut args = std::env::args().skip(1);
    let which = args.next().unwrap_or_else(|| "candidates".to_owned());
    let rounds: usize = args.next().and_then(|n| n.parse().ok()).unwrap_or(200);
    let Ok(groups) = tenspeed::groups() else {
        println!("the 10s dataset did not parse");
        return;
    };
    let emission = StandardEmission::default();
    let mut workspace = Workspace::new();
    let mut out: Vec<Log10Likelihood> = Vec::new();
    let (mut checksum, mut pairs) = (0.0f64, 0usize);

    let start = Instant::now();
    for _ in 0..rounds {
        for group in &groups {
            let reads = group.reads.iter().zip(&group.offsets);
            if which == "align-reads" {
                let refs: Vec<&Haplotype> = black_box(&group.haplotypes).iter().collect();
                let reads: Vec<(&compair::Read, Band)> =
                    reads.map(|(read, offset)| (read, Band::anchored(*offset))).collect();
                workspace.align_reads(&refs, &reads, &emission, &mut out);
                for score in &out {
                    checksum += score.get();
                    pairs += 1;
                }
            } else if which == "prepared" {
                let mut candidates = workspace.candidates(black_box(&group.haplotypes), &emission);
                for (read, offset) in reads {
                    candidates.align(read, Band::anchored(*offset), &mut out);
                    for score in &out {
                        checksum += score.get();
                        pairs += 1;
                    }
                }
            } else {
                let refs: Vec<&Haplotype> = black_box(&group.haplotypes).iter().collect();
                for (read, offset) in reads {
                    workspace.align_candidates(
                        &refs,
                        read,
                        &emission,
                        Band::anchored(*offset),
                        &mut out,
                    );
                    for score in &out {
                        checksum += score.get();
                        pairs += 1;
                    }
                }
            }
        }
    }
    let elapsed = start.elapsed();
    #[allow(clippy::cast_precision_loss, reason = "a per-pair rate for a human to read")]
    let per = elapsed.as_secs_f64() * 1e9 / pairs.max(1) as f64;
    println!("{which}: {rounds} rounds, {pairs} items, {per:.1} ns/item, checksum {checksum}");
}
