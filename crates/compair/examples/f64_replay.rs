//! `align_banded_f64`, the rescue path, over real pairs dumped from rastair,
//! for interleaved A/B runs of the `f64` kernel.
//!
//! ```bash
//! cargo run -p compair --release --example f64_replay -- pairs.tsv below 20
//! ```
//!
//! The file has one pair per line, tab-separated: the strip kernel's raw
//! score, the score after `trusted`, band width, band offset, strand,
//! haplotype, read bases, the four quality tracks, the read's epsilons and
//! the haplotype's site weights. The second argument is `all` (every pair)
//! or `below` (only the pairs whose strip score is below `trusted`'s
//! floor, the rescue's real workload), the third the number of rounds.
#![allow(
    clippy::print_stdout,
    clippy::indexing_slicing,
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    reason = "a benchmark harness over a trusted dump"
)]

use std::{hint::black_box, time::Instant};

use compair::{
    Band, Base, BaseQuality, Emission, HapSite, Haplotype, Observation, Read, SiteWeights, Strand,
    align_banded_f64, align_banded_f64_rows,
};

/// The dumped emission: per-site weights and per-observation epsilons.
struct Replay {
    eps: Vec<f64>,
    sites: Vec<SiteWeights>,
}

impl Emission for Replay {
    fn site_weights(&self, site: HapSite, _: Strand) -> SiteWeights {
        self.sites[site.index.0 as usize]
    }

    fn epsilon(&self, observation: Observation) -> f64 {
        self.eps[observation.index.get() as usize]
    }
}

struct Pair {
    rescued: bool,
    band: Band,
    haplotype: Haplotype,
    read: Read,
    emission: Replay,
}

fn load(path: &str) -> Vec<Pair> {
    let text = std::fs::read_to_string(path).expect("the pairs file");
    text.lines()
        .map(|line| {
            let f: Vec<&str> = line.split('\t').collect();
            let quals = |s: &str| {
                s.split(',')
                    .filter(|x| !x.is_empty())
                    .map(|x| BaseQuality::from_byte(x.parse().unwrap()))
                    .collect::<Vec<_>>()
            };
            let strand = if f[4] == "OT" { Strand::OT } else { Strand::OB };
            let read = Read::new(
                f[6].bytes().map(Base::from).collect::<Vec<_>>(),
                &quals(f[7]),
                &quals(f[8]),
                &quals(f[9]),
                &quals(f[10]),
                strand,
            )
            .unwrap();
            let eps = f[11].split(',').filter(|x| !x.is_empty()).map(|x| x.parse().unwrap());
            let sites = f[12].split(',').filter(|x| !x.is_empty()).map(|x| {
                let b = x.as_bytes();
                SiteWeights {
                    base: Base::from(b[0]),
                    converted: if b[1] == b'-' { None } else { Some(Base::from(b[1])) },
                    rate: x[2..].parse().unwrap(),
                }
            });
            let raw: f64 = f[0].parse().unwrap();
            let width: f64 = f[2].parse().unwrap();
            let (h, r) = (f[5].len() as f64, f[6].len() as f64);
            // `reference::trust_floor` at the strip kernel's scale, 2^96.
            let columns = (2.0 * (width / 2.0).floor() + 1.0).min(h + 1.0);
            let floor = (3.0 * (r + 1.0) * columns).log10()
                - (125.0 + 96.0) * core::f64::consts::LOG10_2
                + 6.0;
            Pair {
                rescued: raw < floor,
                band: Band::new(f[2].parse().unwrap(), f[3].parse().unwrap()).unwrap(),
                haplotype: Haplotype::new(f[5].bytes().map(Base::from).collect::<Vec<_>>()),
                read,
                emission: Replay { eps: eps.collect(), sites: sites.collect() },
            }
        })
        .collect()
}

fn main() {
    let mut args = std::env::args().skip(1);
    let path = args.next().expect("a pairs file");
    let which = args.next().unwrap_or_else(|| "below".to_owned());
    let rounds: usize = args.next().and_then(|n| n.parse().ok()).unwrap_or(1);
    let kernel = args.next().unwrap_or_default();
    let mut pairs = load(&path);
    if which == "below" {
        pairs.retain(|pair| pair.rescued);
    }
    let (mut checksum, mut bits, mut count) = (0.0f64, 0u64, 0usize);
    let start = Instant::now();
    for _ in 0..rounds {
        for pair in &pairs {
            let (haplotype, read) = (black_box(&pair.haplotype), &pair.read);
            let score = if kernel == "rows" {
                align_banded_f64_rows(haplotype, read, &pair.emission, pair.band)
            } else {
                align_banded_f64(haplotype, read, &pair.emission, pair.band)
            }
            .get();
            if kernel == "check" {
                let rows = align_banded_f64_rows(haplotype, read, &pair.emission, pair.band).get();
                assert_eq!(score.to_bits(), rows.to_bits(), "{score} against {rows}");
            }
            checksum += score;
            bits = bits.rotate_left(5) ^ score.to_bits();
            count += 1;
        }
    }
    let elapsed = start.elapsed();
    #[allow(clippy::cast_precision_loss, reason = "a per-pair rate for a human to read")]
    let per = elapsed.as_secs_f64() * 1e9 / count.max(1) as f64;
    println!(
        "{which}: {} pairs x {rounds}, {per:.1} ns/pair, bits {bits:016x} checksum {checksum}",
        pairs.len()
    );
}
