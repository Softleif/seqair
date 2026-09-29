//! rastair's `score_window` in a loop: few haplotypes, many reads, a `Read`
//! built per read. The shape the crate is actually called in, as opposed to
//! the one-read-eight-haplotypes shape of `profile_150x48`.
//!
//! ```bash
//! cargo build --profile profiling -p compair --example shadow_shape
//! samply record --save-only --unstable-presymbolicate -o p.json --rate 4000 -- \
//!   ./target/profiling/examples/shadow_shape rastair 400
//! ```
//!
//! Each window is a 250 bp reference with a 2 bp insertion and a 3 bp
//! deletion applied after its middle base, so three haplotypes, and
//! [`READS`] reads of 100 to 150 bp that span the anchor, cut from one of the
//! three with sequencing errors, TAPS conversions at `CpG`s, per-base
//! qualities and per-base gap-open qualities that drop in homopolymers.
//! Each read's band is centred on where it was cut from, give or take
//! [`JITTER`] bases, at rastair's width for a 3 bp allele: 48.
//!
//! The first argument picks what a round does over every window:
//!
//! - `rastair` (the default) is `score_window`'s loop: the haplotypes built
//!   from the window's bases, then per read a `Read::new` from cloned bases
//!   and one `align_strips_simd` per haplotype;
//! - `prebuilt` scores the same pairs with every `Read` and `Haplotype` built
//!   once up front, so `rastair - prebuilt` is the per-read construction;
//! - `reads` only builds the `Read`s, which is that difference on its own;
//! - `candidates` is `rastair` through `Workspace::candidates`: the same
//!   `Read` per read, but its row tracks derived once for the three
//!   haplotypes and each haplotype's column tracks once per strand for the
//!   window, rather than both once per pair;
//! - `align-reads` builds every `Read` of a window first and scores them all
//!   in one `Workspace::align_reads`, the pairs kernel's entry point;
//! - `f64` scores the prebuilt pairs with `align_banded_f64`, the rescue
//!   path, on every pair.
//!
//! The second argument is the number of rounds.
#![allow(clippy::print_stdout, reason = "this example exists to print a checksum")]

use std::{hint::black_box, time::Instant};

use compair::{
    Band, Base, BaseQuality, Betas, ConversionModel, Haplotype, Log10Likelihood, Probability, Read,
    Strand, TapsEmission, Workspace, align_banded_f64,
};

/// The reference window's length.
const WINDOW: usize = 250;
/// Where the alleles apply: after this base.
const ANCHOR: usize = 124;
/// Distinct windows a round cycles through.
const WINDOWS: usize = 8;
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

fn main() {
    let mut args = std::env::args().skip(1);
    let which = args.next().unwrap_or_else(|| "rastair".to_owned());
    let rounds: usize = args.next().and_then(|n| n.parse().ok()).unwrap_or(200);

    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    let windows: Vec<Window> = (0..WINDOWS).map(|_| window(&mut rng)).collect();
    let emission = TapsEmission::new(
        ConversionModel::taps_default(),
        Betas::Uniform(Probability::new_panicky(0.5)),
    );
    let mut workspace = Workspace::new();
    let mut checksum = 0.0f64;
    let mut alignments = 0usize;
    let mut scores: Vec<Log10Likelihood> = Vec::new();

    let prebuilt: Vec<Prebuilt> = windows
        .iter()
        .map(|w| {
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

    let start = Instant::now();
    for _ in 0..rounds {
        match which.as_str() {
            "prebuilt" => {
                for Prebuilt { haplotypes, reads } in &prebuilt {
                    for (read, band) in reads {
                        for haplotype in haplotypes {
                            let score = workspace.align_strips_simd(
                                black_box(haplotype),
                                read,
                                &emission,
                                *band,
                            );
                            checksum += score.get();
                            alignments += 1;
                        }
                    }
                }
            }
            "candidates" => {
                for w in &windows {
                    let reference = black_box(&w.reference).clone();
                    let haplotypes: Vec<Haplotype> = std::iter::once(reference)
                        .chain(w.alleles.iter().cloned())
                        .map(Haplotype::new)
                        .collect();
                    let mut candidates = workspace.candidates(&haplotypes, &emission);
                    for read in &w.reads {
                        let Ok(band) = Band::new(WIDTH, read.offset) else { continue };
                        let Some(read) = hmm_read(read) else { continue };
                        candidates.align(&read, band, &mut scores);
                        // One addition per score, as the other arms do, so the
                        // checksums compare bit for bit.
                        for score in &scores {
                            checksum += score.get();
                            alignments += 1;
                        }
                    }
                }
            }
            "align-reads" => {
                for w in &windows {
                    let reference = black_box(&w.reference).clone();
                    let haplotypes: Vec<Haplotype> = std::iter::once(reference)
                        .chain(w.alleles.iter().cloned())
                        .map(Haplotype::new)
                        .collect();
                    let refs: Vec<&Haplotype> = haplotypes.iter().collect();
                    let reads: Vec<(Read, Band)> = w
                        .reads
                        .iter()
                        .filter_map(|r| Some((hmm_read(r)?, Band::new(WIDTH, r.offset).ok()?)))
                        .collect();
                    let reads: Vec<(&Read, Band)> =
                        reads.iter().map(|(read, band)| (read, *band)).collect();
                    workspace.align_reads(&refs, &reads, &emission, &mut scores);
                    for score in &scores {
                        checksum += score.get();
                        alignments += 1;
                    }
                }
            }
            "f64" => {
                for Prebuilt { haplotypes, reads } in &prebuilt {
                    for (read, band) in reads {
                        for haplotype in haplotypes {
                            let score =
                                align_banded_f64(black_box(haplotype), read, &emission, *band);
                            checksum += score.get();
                            alignments += 1;
                        }
                    }
                }
            }
            "reads" => {
                for w in &windows {
                    for read in &w.reads {
                        if let Some(read) = hmm_read(black_box(read)) {
                            checksum += f64::from(u32::try_from(read.len()).unwrap_or(0));
                            alignments += 1;
                        }
                    }
                }
            }
            _ => {
                for w in &windows {
                    // `score_window`: the reference as a `Vec<Base>`, then one
                    // `Haplotype::new` per haplotype from owned bases.
                    let reference = black_box(&w.reference).clone();
                    let haplotypes: Vec<Haplotype> = std::iter::once(reference)
                        .chain(w.alleles.iter().cloned())
                        .map(Haplotype::new)
                        .collect();
                    for read in &w.reads {
                        let Ok(band) = Band::new(WIDTH, read.offset) else { continue };
                        let Some(read) = hmm_read(read) else { continue };
                        for haplotype in &haplotypes {
                            let score =
                                workspace.align_strips_simd(haplotype, &read, &emission, band);
                            checksum += score.get();
                            alignments += 1;
                        }
                    }
                }
            }
        }
    }
    let elapsed = start.elapsed();
    #[allow(clippy::cast_precision_loss, reason = "a per-item rate for a human to read")]
    let per = elapsed.as_secs_f64() * 1e9 / alignments.max(1) as f64;
    println!("{which}: {rounds} rounds, {alignments} items, {per:.1} ns/item, checksum {checksum}");
}
