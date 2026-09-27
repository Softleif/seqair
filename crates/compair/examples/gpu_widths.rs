//! The GPU kernel's cost against band width: the 10s pairs, cycled to a
//! saturating batch, at each width, per kernel variant.
//!
//! `cargo run -p compair --release --features gpu --example gpu_widths -- [variant...]`
//! (variants as in `gpu_tenspeed`).

use compair::{
    Band, StandardEmission, Workspace,
    gpu::{Contraction, GpuAligner, GpuContext, GpuPairs, KernelOptions, Style, Subnormals},
};
use std::time::Duration;

#[path = "../tests/support/tenspeed.rs"]
mod tenspeed;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let groups = tenspeed::groups().map_err(|error| format!("10s: {error:?}"))?;
    let emission = StandardEmission::default();
    let context = GpuContext::new()?;
    println!("adapter: {}", context.adapter_info().name);
    let batch: usize = std::env::var("BATCH").ok().and_then(|v| v.parse().ok()).unwrap_or(28400);
    let variants: Vec<String> = std::env::args().skip(1).collect();
    let mut workspace = Workspace::new();
    let widths: Vec<u32> = std::env::var("WIDTHS")
        .ok()
        .map(|v| v.split(',').filter_map(|w| w.parse().ok()).collect())
        .unwrap_or_else(|| vec![32, 46, 48, 64, 80, 96, 128]);
    let mut last_emulated = std::collections::BTreeMap::new();
    for &width in &widths {
        let mut pairs = GpuPairs::new();
        let mut oracle = Vec::new();
        'fill: loop {
            for g in &groups {
                let reads: Vec<_> = g
                    .reads
                    .iter()
                    .map(|r| pairs.push_read(r, &emission))
                    .collect::<Result<_, _>>()?;
                let haps: Vec<_> = g
                    .haplotypes
                    .iter()
                    .map(|h| pairs.push_haplotype(h, compair::Strand::OT, &emission))
                    .collect::<Result<_, _>>()?;
                for ((read, slot), &offset) in g.reads.iter().zip(&reads).zip(&g.offsets) {
                    let band = Band::new(width, offset)?;
                    for (hap, &hap_slot) in g.haplotypes.iter().zip(&haps) {
                        if oracle.len() < 3550 {
                            oracle.push(workspace.align_strips(hap, read, &emission, band));
                        }
                        pairs.push_pair(*slot, hap_slot, band)?;
                        if pairs.len() == batch {
                            break 'fill;
                        }
                    }
                }
            }
        }
        for (label, subnormals) in [("kept", Subnormals::Kept), ("flushed", Subnormals::Flushed)] {
            let emulated = pairs.emulate(subnormals);
            let identical = emulated
                .iter()
                .zip(&oracle)
                .filter(|(a, b)| a.get().to_bits() == b.get().to_bits())
                .count();
            println!(
                "width {width:>3} emulated, subnormals {label:<7} identical {identical}/{}",
                oracle.len()
            );
            last_emulated.insert(label, emulated);
        }
        for name in &variants {
            let mut options = KernelOptions { timestamps: true, ..KernelOptions::default() };
            for part in name.split('+') {
                match part {
                    "loop" => options.style = Some(Style::Loop),
                    "unrolled" => options.style = Some(Style::Unrolled),
                    "blocked" => options.contraction = Contraction::Blocked,
                    "allowed" => options.contraction = Contraction::Allowed,
                    _ => {}
                }
            }
            let mut aligner = GpuAligner::with_options(context.clone(), options);
            let scores = aligner.align(&pairs)?;
            let identical = scores
                .iter()
                .zip(&oracle)
                .filter(|(a, b)| a.get().to_bits() == b.get().to_bits())
                .count();
            for (label, emulated) in &last_emulated {
                let same = scores
                    .iter()
                    .zip(emulated)
                    .filter(|(a, b)| a.get().to_bits() == b.get().to_bits())
                    .count();
                println!(
                    "    {name} vs emulation with subnormals {label}: {same}/{} identical",
                    scores.len()
                );
            }
            let worst = scores
                .iter()
                .zip(&oracle)
                .filter(|(a, b)| a.get().is_finite() && b.get().is_finite())
                .map(|(a, b)| (a.get() - b.get()).abs())
                .fold(0.0, f64::max);
            let mut kernel: Vec<Duration> = Vec::new();
            for _ in 0..15 {
                aligner.align(&pairs)?;
                kernel.extend(aligner.last_kernel_time());
            }
            kernel.sort();
            let best = kernel.first().copied().unwrap_or_default();
            println!(
                "width {width:>3} {name:<18} kernel {:>8.1} µs  {:>6.3} µs/pair  identical {identical}/{}  max|Δ| {worst:.1e}",
                best.as_secs_f64() * 1e6,
                best.as_secs_f64() * 1e6 / batch as f64,
                oracle.len(),
            );
        }
    }
    Ok(())
}
