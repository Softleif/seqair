//! rastair's shape on the GPU: many windows of three haplotypes and 128 reads
//! each, every (read, haplotype) pair of a batch of windows in one launch,
//! against the CPU strip kernel on the same pairs.
//!
//! The windows are `shadow_shape`'s, copied so the two measure the same
//! pairs: a 250 bp reference with a 2 bp insertion and a 3 bp deletion after
//! its middle base, reads of 100 to 150 bp with
//! sequencing errors, TAPS conversions, per-base qualities and gap-open
//! qualities, bands of width 48 centred near where each read was cut.
//!
//! `cargo run -p compair --release --features gpu --example gpu_shadow -- [windows per launch...]`
#![allow(clippy::print_stdout, reason = "this example exists to print a table")]

use std::time::{Duration, Instant};

use compair::{
    Band, Base, BaseQuality, Betas, ConversionModel, Haplotype, Log10Likelihood, Probability, Read,
    Strand, TapsEmission, Workspace,
    gpu::{GpuAligner, GpuContext, GpuPairs, KernelOptions},
};
/// The reference window's length.
const WINDOW: usize = 250;
/// Where the alleles apply: after this base.
const ANCHOR: usize = 124;
/// Reads per window.
const READS: usize = 128;
/// How far a band's centre may sit from where its read was cut.
const JITTER: u64 = 10;
/// rastair's width for a 3 bp allele: `max(46, 2 * 3 + 32)` rounded up to 8.
const WIDTH: u32 = 48;

/// xorshift64, so the workload is the same on every machine.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    /// Uniform in `0..n`, `n > 0`.
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }

    /// `true` with probability `numerator / 10_000`.
    fn chance(&mut self, numerator: u64) -> bool {
        self.below(10_000) < numerator
    }

    fn base(&mut self) -> Base {
        Base::KNOWN.get(usize::try_from(self.below(4)).unwrap_or(0)).copied().unwrap_or(Base::A)
    }
}

/// One read as rastair holds it before it becomes a `compair::Read`.
struct ShadowRead {
    bases: Vec<Base>,
    quals: Vec<BaseQuality>,
    gap_open: Vec<BaseQuality>,
    gap_extend: Vec<BaseQuality>,
    offset: i32,
    strand: Strand,
}

struct Window {
    reference: Vec<Base>,
    alleles: Vec<Vec<Base>>,
    reads: Vec<ShadowRead>,
}

fn window(rng: &mut Rng) -> Window {
    let reference: Vec<Base> = (0..WINDOW).map(|_| rng.base()).collect();
    let (head, tail) = reference.split_at(ANCHOR + 1);
    let insertion: Vec<Base> =
        head.iter().copied().chain([rng.base(), rng.base()]).chain(tail.iter().copied()).collect();
    let deletion: Vec<Base> = head.iter().copied().chain(tail.iter().skip(3).copied()).collect();
    let sources = [reference.clone(), insertion.clone(), deletion.clone()];

    let reads = (0..READS)
        .map(|_| {
            let source =
                sources.get(usize::try_from(rng.below(3)).unwrap_or(0)).unwrap_or(&reference);
            let len = 100 + usize::try_from(rng.below(51)).unwrap_or(0);
            // Spanning the anchor with ten bases to spare on either side.
            let lowest = (ANCHOR + 10).saturating_sub(len);
            let highest = (ANCHOR - 10).min(source.len().saturating_sub(len));
            let span = u64::try_from(highest.saturating_sub(lowest) + 1).unwrap_or(1);
            let start = lowest + usize::try_from(rng.below(span)).unwrap_or(0);
            let strand = if rng.below(2) == 0 { Strand::OT } else { Strand::OB };
            let cut = source.get(start..start + len).unwrap_or_default();
            let bases: Vec<Base> = cut
                .iter()
                .enumerate()
                .map(|(index, &base)| {
                    let next = cut.get(index + 1).copied();
                    let previous = index.checked_sub(1).and_then(|i| cut.get(i)).copied();
                    // TAPS: a methylated CpG converts at ~half, anything
                    // else at the false-conversion rate.
                    let converted = match (strand, base) {
                        (Strand::OT, Base::C) => rng
                            .chance(if next == Some(Base::G) { 4_800 } else { 40 })
                            .then_some(Base::T),
                        (Strand::OB, Base::G) => rng
                            .chance(if previous == Some(Base::C) { 4_800 } else { 40 })
                            .then_some(Base::A),
                        _ => None,
                    };
                    let base = converted.unwrap_or(base);
                    if rng.chance(100) { rng.base() } else { base }
                })
                .collect();
            let quals = (0..len)
                .map(|_| {
                    BaseQuality::from_byte(if rng.chance(300) {
                        2
                    } else {
                        20 + u8::try_from(rng.below(21)).unwrap_or(0)
                    })
                })
                .collect();
            // Gap opens cheaper, and extends much cheaper, inside a run.
            let in_run = |index: usize| {
                let at = |i: usize| bases.get(i).copied();
                index >= 2 && at(index) == at(index - 1) && at(index) == at(index - 2)
            };
            let gap_open = (0..len)
                .map(|index| BaseQuality::from_byte(if in_run(index) { 30 } else { 45 }))
                .collect();
            let gap_extend = (0..len)
                .map(|index| BaseQuality::from_byte(if in_run(index) { 3 } else { 10 }))
                .collect();
            let jitter = i64::try_from(rng.below(2 * JITTER + 1)).unwrap_or(0)
                - i64::try_from(JITTER).unwrap_or(0);
            let offset = i32::try_from(i64::try_from(start).unwrap_or(0) + jitter).unwrap_or(0);
            ShadowRead { bases, quals, gap_open, gap_extend, offset, strand }
        })
        .collect();
    Window { reference, alleles: vec![insertion, deletion], reads }
}

/// A window's haplotypes and reads, built once.
struct Prebuilt {
    haplotypes: Vec<Haplotype>,
    reads: Vec<(Read, Band)>,
}

fn hmm_read(read: &ShadowRead) -> Option<Read> {
    Read::new(
        read.bases.clone(),
        &read.quals,
        &read.gap_open,
        &read.gap_open,
        &read.gap_extend,
        read.strand,
    )
    .ok()
}

fn fill(
    pairs: &mut GpuPairs,
    windows: &[Prebuilt],
    emission: &TapsEmission<'_>,
) -> Result<(), compair::gpu::GpuError> {
    pairs.clear();
    for window in windows {
        // Each haplotype once per strand: TAPS weights depend on it.
        let mut slots = Vec::with_capacity(2 * window.haplotypes.len());
        for strand in [Strand::OT, Strand::OB] {
            for haplotype in &window.haplotypes {
                slots.push((strand, pairs.push_haplotype(haplotype, strand, emission)?));
            }
        }
        for (read, band) in &window.reads {
            let slot = pairs.push_read(read, emission)?;
            for &(strand, haplotype) in &slots {
                if strand == read.strand() {
                    pairs.push_pair(slot, haplotype, *band)?;
                }
            }
        }
    }
    Ok(())
}

fn cpu(
    workspace: &mut Workspace,
    windows: &[Prebuilt],
    emission: &TapsEmission<'_>,
    out: &mut Vec<Log10Likelihood>,
) {
    out.clear();
    for window in windows {
        for (read, band) in &window.reads {
            for haplotype in &window.haplotypes {
                out.push(workspace.align_strips_simd(haplotype, read, emission, *band));
            }
        }
    }
}

fn best<E>(runs: usize, mut run: impl FnMut() -> Result<(), E>) -> Result<Duration, E> {
    run()?;
    let mut best = Duration::MAX;
    for _ in 0..runs {
        let started = Instant::now();
        run()?;
        best = best.min(started.elapsed());
    }
    Ok(best)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let launches: Vec<usize> = {
        let args: Vec<usize> = std::env::args().skip(1).filter_map(|a| a.parse().ok()).collect();
        if args.is_empty() { vec![1, 4, 16, 64, 256] } else { args }
    };
    let most = launches.iter().copied().max().unwrap_or(1);
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    let prebuilt: Vec<Prebuilt> = (0..most)
        .map(|_| {
            let w = window(&mut rng);
            let haplotypes = std::iter::once(&w.reference)
                .chain(&w.alleles)
                .map(|bases| Haplotype::new(bases.clone()))
                .collect();
            let reads = w
                .reads
                .iter()
                .filter_map(|r| Some((hmm_read(r)?, Band::new(WIDTH, r.offset).ok()?)))
                .collect();
            Prebuilt { haplotypes, reads }
        })
        .collect();
    let emission =
        TapsEmission::new(ConversionModel::taps_default(), Betas::Uniform(Probability::new(0.5)?));

    let context = GpuContext::new()?;
    println!("adapter: {} ({:?})", context.adapter_info().name, context.adapter_info().backend);
    let mut aligner = GpuAligner::with_options(
        context.clone(),
        KernelOptions { timestamps: true, ..KernelOptions::default() },
    );
    let mut workspace = Workspace::new();
    let mut pairs = GpuPairs::new();
    let mut strips = Vec::new();
    let threads = std::thread::available_parallelism().map_or(1, usize::from);
    println!(
        "windows | pairs | fill µs | launch µs | kernel µs | GPU µs/pair (launch) | CPU 1-thread µs/pair | CPU {threads}-thread µs/pair | identical | max |Δ|"
    );
    for &count in &launches {
        let windows = prebuilt.get(..count).unwrap_or_default();
        fill(&mut pairs, windows, &emission)?;
        let gpu = aligner.align(&pairs)?;
        cpu(&mut workspace, windows, &emission, &mut strips);
        let identical =
            gpu.iter().zip(&strips).filter(|(a, b)| a.get().to_bits() == b.get().to_bits()).count();
        let worst = gpu
            .iter()
            .zip(&strips)
            .filter(|(a, b)| a.get().is_finite() && b.get().is_finite())
            .map(|(a, b)| (a.get() - b.get()).abs())
            .fold(0.0, f64::max);

        let runs = (2000 / count).clamp(5, 200);
        let fill_time = best(runs, || fill(&mut pairs, windows, &emission))?;
        let mut kernel = Duration::MAX;
        let launch = best(runs, || {
            aligner.align(&pairs)?;
            kernel = kernel.min(aligner.last_kernel_time().unwrap_or(Duration::MAX));
            Ok::<(), compair::gpu::GpuError>(())
        })?;
        let single = best(runs.min(20), || {
            cpu(&mut workspace, windows, &emission, &mut strips);
            Ok::<(), ()>(())
        })
        .unwrap_or_default();
        let all = best(runs.min(20), || {
            std::thread::scope(|scope| {
                for chunk in windows.chunks(windows.len().div_ceil(threads)) {
                    let emission = &emission;
                    scope.spawn(move || {
                        let mut workspace = Workspace::new();
                        let mut out = Vec::new();
                        cpu(&mut workspace, chunk, emission, &mut out);
                        std::hint::black_box(out);
                    });
                }
            });
            Ok::<(), ()>(())
        })
        .unwrap_or_default();
        let mut submit = Duration::MAX;
        let mut collect = Duration::MAX;
        for _ in 0..runs.min(10) {
            let started = Instant::now();
            let handle = aligner.submit(&pairs)?;
            let submitted = started.elapsed();
            let started = Instant::now();
            handle.collect()?;
            submit = submit.min(submitted);
            collect = collect.min(started.elapsed());
        }
        println!(
            "    submit (upload + encode) {:.1} µs, collect (wait + readback) {:.1} µs, rows {:.1} MB",
            submit.as_secs_f64() * 1e6,
            collect.as_secs_f64() * 1e6,
            pairs.upload_bytes() as f64 / 1e6
        );
        let n = pairs.len() as f64;
        let us = |d: Duration| d.as_secs_f64() * 1e6;
        println!(
            "{count:>7} | {:>6} | {:>8.1} | {:>9.1} | {:>9.1} | {:>8.3} | {:>8.3} | {:>8.3} | {identical}/{} | {worst:.1e}",
            pairs.len(),
            us(fill_time),
            us(launch),
            if kernel == Duration::MAX { f64::NAN } else { us(kernel) },
            us(launch) / n,
            us(single) / n,
            us(all) / n,
            strips.len(),
        );
    }
    Ok(())
}
