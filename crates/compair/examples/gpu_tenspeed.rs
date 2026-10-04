//! The GPU kernel on the 10s benchmark: parity against the strip kernel, and
//! time per launch against the CPU kernels at several batch sizes.
//!
//! ```text
//! cargo run -p compair --release --features gpu --example gpu_tenspeed -- [variant...]
//! ```
//!
//! A variant is `loop` or `unrolled`, optionally `+blocked` (contraction
//! blocked) and `@N` (workgroup size), e.g. `unrolled+blocked@32`. With none,
//! the default kernel. Every variant prints its parity line and its timings.

use compair::{
    Band, Haplotype, Log10Likelihood, StandardEmission, Workspace,
    gpu::{Contraction, GpuAligner, GpuContext, GpuPairs, KernelOptions, Style, Wait},
};
use std::time::{Duration, Instant};

#[path = "../tests/support/tenspeed.rs"]
mod tenspeed;

/// One pair of the benchmark, by group, read and haplotype index.
#[derive(Clone, Copy)]
struct Index {
    group: usize,
    read: usize,
    haplotype: usize,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let groups = tenspeed::groups().map_err(|error| format!("10s: {error:?}"))?;
    let index: Vec<Index> = groups
        .iter()
        .enumerate()
        .flat_map(|(group, g)| {
            (0..g.reads.len()).flat_map(move |read| {
                (0..g.haplotypes.len()).map(move |haplotype| Index { group, read, haplotype })
            })
        })
        .collect();
    let emission = StandardEmission::default();
    let band_of = |at: Index| {
        groups
            .get(at.group)
            .and_then(|g| g.offsets.get(at.read))
            .map_or(Band::anchored(0), |&offset| Band::anchored(offset))
    };
    let cells: u64 = index
        .iter()
        .map(|&at| {
            let g = &groups[at.group];
            band_cells(g.reads[at.read].len(), g.haplotypes[at.haplotype].len(), band_of(at))
        })
        .sum();
    println!("10s: {} pairs, {cells} band cells", index.len());

    let context = GpuContext::new()?;
    let info = context.adapter_info();
    println!(
        "adapter: {} ({:?}, {:?}), uma {}",
        info.name,
        info.backend,
        info.device_type,
        context.is_uma()
    );

    // The oracle, once.
    let mut workspace = Workspace::new();
    let oracle: Vec<Log10Likelihood> = index
        .iter()
        .map(|&at| {
            let g = &groups[at.group];
            workspace.align_strips(
                &g.haplotypes[at.haplotype],
                &g.reads[at.read],
                &emission,
                band_of(at),
            )
        })
        .collect();

    let fill = |count: usize, pairs: &mut GpuPairs| -> Result<(), compair::gpu::GpuError> {
        pairs.clear();
        let mut last: Option<(usize, usize, compair::gpu::ReadSlot)> = None;
        let mut haps: Vec<Option<compair::gpu::HaplotypeSlot>> = Vec::new();
        let mut hap_group = usize::MAX;
        for &at in index.iter().cycle().take(count) {
            let g = &groups[at.group];
            if hap_group != at.group {
                hap_group = at.group;
                haps.clear();
                haps.resize(g.haplotypes.len(), None);
            }
            let read = &g.reads[at.read];
            let read_slot = match last {
                Some((group, r, slot)) if group == at.group && r == at.read => slot,
                _ => {
                    let slot = pairs.push_read(read, &emission)?;
                    last = Some((at.group, at.read, slot));
                    slot
                }
            };
            let hap_slot = match haps[at.haplotype] {
                Some(slot) => slot,
                None => {
                    let slot = pairs.push_haplotype(
                        &g.haplotypes[at.haplotype],
                        read.strand(),
                        &emission,
                    )?;
                    haps[at.haplotype] = Some(slot);
                    slot
                }
            };
            pairs.push_pair(read_slot, hap_slot, band_of(at))?;
        }
        Ok(())
    };

    let mut variants: Vec<(String, KernelOptions)> =
        std::env::args().skip(1).map(|arg| (arg.clone(), parse_variant(&arg))).collect();
    if variants.is_empty() {
        variants.push(("default".into(), KernelOptions::default()));
    }

    for (label, wait) in [("spin", Wait::Spin), ("block", Wait::Block)] {
        let times = repeat(0, || context.round_trip(wait))?;
        let (min, median) = stats(&times);
        println!(
            "empty round trip ({label}): {:.1} / {:.1} µs (min / median)",
            us(min),
            us(median)
        );
    }

    let mut pairs = GpuPairs::new();
    for (name, options) in &variants {
        let mut aligner = GpuAligner::with_options(context.clone(), *options);
        fill(index.len(), &mut pairs)?;
        let started = Instant::now();
        // The kernel's own scores: this compares kernels, before any rescue.
        let scores: Vec<_> = aligner.align(&pairs)?.iter().map(|s| s.untrusted()).collect();
        println!("\n== {name}: first launch (compiles) {:.1} ms", ms(started.elapsed()));
        parity(&scores, &oracle);

        println!("  batch | fill µs | launch µs (min / median) | µs/pair | GCUPS (band)");
        for size in [1usize, 8, 64, 512, 3550, 14200, 56800] {
            let started = Instant::now();
            fill(size, &mut pairs)?;
            let fill_time = started.elapsed();
            let mut kernel = Vec::new();
            let times = repeat(size, || {
                aligner.align(&pairs)?;
                kernel.extend(aligner.last_kernel_time());
                Ok::<(), compair::gpu::GpuError>(())
            })?;
            let (min, median) = stats(&times);
            if !kernel.is_empty() {
                let (kmin, kmedian) = stats(&kernel);
                print!("  [kernel {:.1} / {:.1} µs] ", us(kmin), us(kmedian));
            }
            // The same launches back to back with the GPU kept busy by a
            // second in flight, to see what the idle clock costs.
            let mut spare = GpuAligner::with_options(context.clone(), *options);
            let busy = repeat(size, || {
                let first = aligner.submit(&pairs)?;
                let second = spare.submit(&pairs)?;
                first.collect()?;
                second.collect().map(|_| ())
            })?;
            let (_, busy_median) = stats(&busy);
            print!("  [2 in flight: {:.1} µs per launch] ", us(busy_median) / 2.0);
            let band: u64 = index
                .iter()
                .cycle()
                .take(size)
                .map(|&at| {
                    let g = &groups[at.group];
                    band_cells(
                        g.reads[at.read].len(),
                        g.haplotypes[at.haplotype].len(),
                        band_of(at),
                    )
                })
                .sum();
            println!(
                "  {size:>5} | {:>7.1} | {:>9.1} / {:>9.1} | {:>7.2} | {:.2}",
                us(fill_time),
                us(min),
                us(median),
                us(median) / size as f64,
                band as f64 / median.as_secs_f64() / 1e9,
            );
        }
    }

    println!("\n== CPU (this machine, one thread unless stated)");
    println!("  batch | strips-simd µs | µs/pair | candidates µs (whole groups)");
    for size in [1usize, 8, 64, 512, 3550] {
        let size = size.min(index.len());
        let times = repeat(size, || {
            let mut sum = 0.0;
            for &at in index.iter().take(size) {
                let g = &groups[at.group];
                sum += workspace
                    .align_strips_simd(
                        &g.haplotypes[at.haplotype],
                        &g.reads[at.read],
                        &emission,
                        band_of(at),
                    )
                    .get();
            }
            std::hint::black_box(sum);
            Ok::<(), compair::gpu::GpuError>(())
        })?;
        let (_, median) = stats(&times);
        println!("  {size:>5} | {:>14.1} | {:>7.2}", us(median), us(median) / size as f64);
    }
    let mut out = Vec::new();
    let times = repeat(3550, || {
        for g in &groups {
            let refs: Vec<&Haplotype> = g.haplotypes.iter().collect();
            for (read, &offset) in g.reads.iter().zip(&g.offsets) {
                workspace.align_candidates(
                    &refs,
                    read,
                    &emission,
                    Band::anchored(offset),
                    &mut out,
                );
            }
        }
        Ok::<(), compair::gpu::GpuError>(())
    })?;
    let (_, median) = stats(&times);
    println!("  all 3550 through align_candidates: {:.1} µs", us(median));
    let threads = std::thread::available_parallelism().map_or(1, usize::from);
    let times = repeat(3550, || {
        std::thread::scope(|scope| {
            for chunk in index.chunks(index.len().div_ceil(threads)) {
                let (groups, emission) = (&groups, &emission);
                scope.spawn(move || {
                    let mut workspace = Workspace::new();
                    let mut sum = 0.0;
                    for &at in chunk {
                        let g = &groups[at.group];
                        sum += workspace
                            .align_strips_simd(
                                &g.haplotypes[at.haplotype],
                                &g.reads[at.read],
                                emission,
                                band_of(at),
                            )
                            .get();
                    }
                    std::hint::black_box(sum);
                });
            }
        });
        Ok::<(), compair::gpu::GpuError>(())
    })?;
    let (_, median) = stats(&times);
    println!("  all 3550 through strips-simd on {threads} threads: {:.1} µs", us(median));
    Ok(())
}

fn parse_variant(arg: &str) -> KernelOptions {
    let mut options = KernelOptions::default();
    let (body, workgroup) = arg.split_once('@').unwrap_or((arg, ""));
    if let Ok(size) = workgroup.parse() {
        options.workgroup_size = size;
    }
    for part in body.split('+') {
        match part {
            "loop" => options.style = Some(Style::Loop),
            "unrolled" => options.style = Some(Style::Unrolled),
            "blocked" => options.contraction = Contraction::Blocked,
            "allowed" => options.contraction = Contraction::Allowed,
            "block" => options.wait = Wait::Block,
            "timed" => options.timestamps = true,
            _ => {}
        }
    }
    options
}

fn parity(scores: &[Log10Likelihood], oracle: &[Log10Likelihood]) {
    let mut identical = 0usize;
    let mut worst = 0.0f64;
    let mut finite_mismatch = 0usize;
    for (gpu, cpu) in scores.iter().zip(oracle) {
        let (g, c) = (gpu.get(), cpu.get());
        if g.to_bits() == c.to_bits() {
            identical += 1;
        } else if g.is_finite() && c.is_finite() {
            worst = worst.max((g - c).abs());
        } else {
            finite_mismatch += 1;
        }
    }
    println!(
        "  parity vs align_strips: {identical}/{} bit-identical ({:.1}%), max |Δlog10| {worst:.3e}, finite/impossible disagreements {finite_mismatch}",
        oracle.len(),
        100.0 * identical as f64 / oracle.len() as f64
    );
}

/// Enough runs of `run` to fill ~0.5 s, at least 5.
fn repeat<E>(_size: usize, mut run: impl FnMut() -> Result<(), E>) -> Result<Vec<Duration>, E> {
    run()?;
    let target = Duration::from_millis(500);
    let mut times = Vec::new();
    let started = Instant::now();
    while times.len() < 5 || (started.elapsed() < target && times.len() < 2000) {
        let t = Instant::now();
        run()?;
        times.push(t.elapsed());
    }
    Ok(times)
}

fn stats(times: &[Duration]) -> (Duration, Duration) {
    let mut sorted = times.to_vec();
    sorted.sort();
    (
        sorted.first().copied().unwrap_or_default(),
        sorted.get(sorted.len() / 2).copied().unwrap_or_default(),
    )
}

fn us(d: Duration) -> f64 {
    d.as_secs_f64() * 1e6
}

fn ms(d: Duration) -> f64 {
    d.as_secs_f64() * 1e3
}

fn band_cells(read: usize, haplotype: usize, band: Band) -> u64 {
    let half = i64::from(band.width() / 2);
    let offset = i64::from(band.offset());
    let (rows, columns) = (read as i64, haplotype as i64);
    (0..rows)
        .map(|row| {
            let low = (row + offset - half).max(0);
            let high = (row + offset + half).min(columns - 1);
            (high - low + 1).max(0) as u64
        })
        .sum()
}
