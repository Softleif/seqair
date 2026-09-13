//! `align_batch` against the per-haplotype kernels on the 10s benchmark.
//!
//! The 10s dataset is grouped the way a caller is: each group is a set of reads
//! and the candidate haplotypes for that region. Group 5 alone is 110 reads x
//! **24** haplotypes, so batches fill completely -- which is the condition
//! `align_batch` needs to pay off, and the condition the 150x46 fixture's
//! 8-haplotype group only just meets.
//!
//! `cargo run -p compair --release --example tenspeed_batch`

use compair::{
    Band, Base, BaseQuality, Haplotype, Log10Likelihood, Read, StandardEmission, Strand, Workspace,
    align_strips_simd,
};
use std::time::Instant;

struct Group {
    reads: Vec<Read>,
    /// One band anchor per read, seeded against the first haplotype.
    offsets: Vec<i32>,
    haplotypes: Vec<Haplotype>,
}

fn main() {
    let groups = parse();
    let pairs: usize = groups.iter().map(|g| g.reads.len() * g.haplotypes.len()).sum();
    let cells: u64 = groups
        .iter()
        .flat_map(|g| {
            g.reads.iter().flat_map(move |r| {
                g.haplotypes.iter().map(move |h| r.len() as u64 * h.len() as u64)
            })
        })
        .sum();
    println!("10s: {} groups, {pairs} pairs, {cells} matrix cells", groups.len());
    for (index, group) in groups.iter().enumerate() {
        println!(
            "  group {index}: {:>3} reads x {:>2} haplotypes",
            group.reads.len(),
            group.haplotypes.len()
        );
    }

    let standard = StandardEmission::default();
    let mut workspace = Workspace::new();
    let mut out: Vec<Log10Likelihood> = Vec::new();

    // Correctness first: a batch must equal the strip kernel haplotype by
    // haplotype, or its timing means nothing.
    let mut worst = 0.0f64;
    for group in &groups {
        let refs: Vec<&Haplotype> = group.haplotypes.iter().collect();
        for (read, offset) in group.reads.iter().zip(&group.offsets) {
            let band = Band::anchored(*offset);
            out.clear();
            workspace.align_batch(&refs, read, &standard, band, &mut out);
            for (haplotype, got) in group.haplotypes.iter().zip(&out) {
                let reference = align_strips_simd(haplotype, read, &standard, band).get();
                worst = worst.max((got.get() - reference).abs());
            }
        }
    }
    println!("\nworst |batch - strips| over all {pairs} pairs: {worst:.3e}");

    let mut worst_dispatch = 0.0f64;
    for group in &groups {
        let refs: Vec<&Haplotype> = group.haplotypes.iter().collect();
        for (read, offset) in group.reads.iter().zip(&group.offsets) {
            let band = Band::anchored(*offset);
            workspace.align_candidates(&refs, read, &standard, band, &mut out);
            for (haplotype, got) in group.haplotypes.iter().zip(&out) {
                let reference = align_strips_simd(haplotype, read, &standard, band).get();
                worst_dispatch = worst_dispatch.max((got.get() - reference).abs());
            }
        }
    }
    println!("worst |dispatch - strips| over all {pairs} pairs: {worst_dispatch:.3e}");

    let repeats = 20;
    let strips = time(repeats, || {
        let mut sum = 0.0;
        for group in &groups {
            for (read, offset) in group.reads.iter().zip(&group.offsets) {
                let band = Band::anchored(*offset);
                for haplotype in &group.haplotypes {
                    sum += workspace.align_strips_simd(haplotype, read, &standard, band).get();
                }
            }
        }
        sum
    });
    let batch = time(repeats, || {
        let mut sum = 0.0;
        for group in &groups {
            let refs: Vec<&Haplotype> = group.haplotypes.iter().collect();
            for (read, offset) in group.reads.iter().zip(&group.offsets) {
                workspace.align_batch(&refs, read, &standard, Band::anchored(*offset), &mut out);
                sum += out.iter().map(|score| score.get()).sum::<f64>();
            }
        }
        sum
    });
    let dispatch = time(repeats, || {
        let mut sum = 0.0;
        for group in &groups {
            let refs: Vec<&Haplotype> = group.haplotypes.iter().collect();
            for (read, offset) in group.reads.iter().zip(&group.offsets) {
                workspace.align_candidates(
                    &refs,
                    read,
                    &standard,
                    Band::anchored(*offset),
                    &mut out,
                );
                sum += out.iter().map(|score| score.get()).sum::<f64>();
            }
        }
        sum
    });

    println!("\n{:<28} {:>10} {:>12}", "arm", "ms / pass", "us / pair");
    for (name, ms) in [
        ("strips-simd, per haplotype", strips),
        ("align_batch (forced)", batch),
        ("align_candidates (dispatch)", dispatch),
    ] {
        println!("{name:<28} {ms:>10.2} {:>12.3}", ms * 1000.0 / pairs as f64);
    }
    println!(
        "\nspeedup over strips: forced batch {:.2}x, dispatch {:.2}x",
        strips / batch,
        strips / dispatch
    );

    println!(
        "\nper group -- lane fill is the whole story:\n{:>6} {:>7} {:>6} {:>9} {:>8} {:>8} {:>8} {:>8}",
        "group", "reads", "haps", "mean len", "fill", "strips", "batch", "dispatch"
    );
    for (index, group) in groups.iter().enumerate() {
        let refs: Vec<&Haplotype> = group.haplotypes.iter().collect();
        let mean =
            group.reads.iter().map(Read::len).sum::<usize>() as f64 / group.reads.len() as f64;
        // Lanes actually used: the last chunk of a group is partial.
        let chunks = group.haplotypes.len().div_ceil(8);
        let fill = group.haplotypes.len() as f64 / (chunks * 8) as f64;
        let s = time(repeats, || {
            let mut sum = 0.0;
            for (read, offset) in group.reads.iter().zip(&group.offsets) {
                let band = Band::anchored(*offset);
                for haplotype in &group.haplotypes {
                    sum += workspace.align_strips_simd(haplotype, read, &standard, band).get();
                }
            }
            sum
        });
        let b = time(repeats, || {
            let mut sum = 0.0;
            for (read, offset) in group.reads.iter().zip(&group.offsets) {
                workspace.align_batch(&refs, read, &standard, Band::anchored(*offset), &mut out);
                sum += out.iter().map(|score| score.get()).sum::<f64>();
            }
            sum
        });
        let d = time(repeats, || {
            let mut sum = 0.0;
            for (read, offset) in group.reads.iter().zip(&group.offsets) {
                workspace.align_candidates(
                    &refs,
                    read,
                    &standard,
                    Band::anchored(*offset),
                    &mut out,
                );
                sum += out.iter().map(|score| score.get()).sum::<f64>();
            }
            sum
        });
        println!(
            "{index:>6} {:>7} {:>6} {mean:>9.0} {:>7.0}% {s:>8.2} {b:>8.2} {d:>8.2}   {:.2}x {:.2}x",
            group.reads.len(),
            group.haplotypes.len(),
            fill * 100.0,
            s / b,
            s / d
        );
    }
}

fn time(repeats: u32, mut body: impl FnMut() -> f64) -> f64 {
    let mut best = f64::INFINITY;
    for _ in 0..repeats {
        let start = Instant::now();
        let sum = body();
        let elapsed = start.elapsed().as_secs_f64() * 1000.0;
        std::hint::black_box(sum);
        best = best.min(elapsed);
    }
    best
}

fn parse() -> Vec<Group> {
    let input = include_str!("../tests/data/10s.in");
    let mut lines = input.lines();
    let mut groups = Vec::new();
    while let Some(header) = lines.next() {
        let mut counts = header.split_whitespace();
        let (Some(nr), Some(nh)) = (counts.next(), counts.next()) else { break };
        let (nr, nh): (usize, usize) = (nr.parse().unwrap_or(0), nh.parse().unwrap_or(0));
        let mut raw = Vec::with_capacity(nr);
        for _ in 0..nr {
            let line = lines.next().unwrap_or_default();
            let f: Vec<&str> = line.split_whitespace().collect();
            let q = |s: &str, min: u8| -> Vec<BaseQuality> {
                s.bytes().map(|b| BaseQuality::from_byte(b.saturating_sub(33).max(min))).collect()
            };
            let bases: Vec<Base> = f[0].bytes().map(Base::from).collect();
            raw.push(
                Read::new(bases, &q(f[1], 6), &q(f[2], 0), &q(f[3], 0), &q(f[4], 0), Strand::OT)
                    .expect("10s reads are well formed"),
            );
        }
        let haplotypes: Vec<Haplotype> = (0..nh)
            .map(|_| {
                Haplotype::from_ascii(
                    lines
                        .next()
                        .unwrap_or_default()
                        .split_whitespace()
                        .next()
                        .unwrap_or("")
                        .as_bytes(),
                )
            })
            .collect();
        let offsets =
            raw.iter().map(|read| haplotypes.first().map_or(0, |h| seed_offset(h, read))).collect();
        groups.push(Group { reads: raw, offsets, haplotypes });
    }
    groups
}

fn seed_offset(haplotype: &Haplotype, read: &Read) -> i32 {
    let (h, r) = (haplotype.len() as i64, read.len() as i64);
    let mut best = (0i64, -1i64);
    for offset in -(r - 1)..h {
        let mut hits = 0i64;
        for (index, base) in read.bases().iter().enumerate() {
            let column = offset + index as i64;
            if column >= 0 && haplotype.bases().get(column as usize) == Some(base) {
                hits += 1;
            }
        }
        if hits > best.1 {
            best = (offset, hits);
        }
    }
    i32::try_from(best.0).unwrap_or(0)
}
