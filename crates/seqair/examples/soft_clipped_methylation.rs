#![allow(
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    clippy::collapsible_if,
    clippy::doc_markdown,
    clippy::indexing_slicing,
    clippy::print_stdout,
    reason = "example"
)]

use clap::Parser;
use seqair::{
    bam::{AlignedPairWithRead, RecordStore},
    reader::{Segment, SegmentOptions},
};
use seqair_types::{Base, Pos0};
use std::num::NonZeroU32;

/// Find CpG positions where one half is soft-clipped and the other is aligned.
///
/// Detects six patterns without consulting the reference — the read's own
/// sequence at the boundary tells the story:
///
/// * Leading CpG: soft-clip ends in C, first aligned base is G
/// * Leading TpG: soft-clip ends in T (methylated C), first aligned base is G
/// * Leading CpA: soft-clip ends in C, first aligned base is A
/// * Trailing CpG: last aligned base is C, soft-clip starts with G
/// * Trailing CpA: last aligned base is C, soft-clip starts with A
/// * Trailing TpG: last aligned base is T, soft-clip starts with G
///
/// Control (aligned CpG): consecutive M/=/X C→G at adjoining reference positions.
#[derive(Parser)]
struct Cli {
    /// Path to sorted and indexed BAM/CRAM file
    #[arg(value_hint = clap::ValueHint::FilePath)]
    input: std::path::PathBuf,

    /// Path to reference FASTA (required for CRAM)
    #[arg(short = 'r', long, value_hint = clap::ValueHint::FilePath)]
    reference: std::path::PathBuf,

    /// Region to process (e.g. "chr1:1000-2000"). Omit for whole-genome scan.
    #[arg(short, long)]
    region: Option<String>,
}

fn main() -> anyhow::Result<()> {
    use anyhow::Context;

    let args = Cli::parse();

    let mut readers = seqair::Readers::open(&args.input, &args.reference)
        .context("could not open alignment file")?;

    let opts =
        SegmentOptions::new(NonZeroU32::new(10_000).context("segment size must be non-zero")?);
    let plan: Vec<Segment> = if let Some(region) = args.region.as_ref() {
        readers.segments(region.as_str(), opts).context("could not plan region")?.collect()
    } else {
        readers.segments((), opts).context("could not plan whole-genome scan")?.collect()
    };

    eprintln!("Scanning {} segment(s)…", plan.len());

    let mut total = Output::default();

    for segment in &plan {
        let mut store = RecordStore::<()>::new();
        readers.fetch_into(segment.tid().as_u32(), segment.start(), segment.end(), &mut store)?;
        total += count_segment(&store);
    }

    print_results(&total);

    Ok(())
}

#[derive(Debug, Default)]
struct Output {
    leading_cs: Count,
    leading_ts: Count,
    leading_ca: Count,
    trailing_cs: Count,
    trailing_ts: Count,
    trailing_tg: Count,
}

#[derive(Debug, Default)]
struct Count {
    clipped: u64,
    not_clipped: u64,
}

impl std::ops::AddAssign for Output {
    fn add_assign(&mut self, rhs: Self) {
        self.leading_cs.clipped += rhs.leading_cs.clipped;
        self.leading_cs.not_clipped += rhs.leading_cs.not_clipped;
        self.leading_ts.clipped += rhs.leading_ts.clipped;
        self.leading_ts.not_clipped += rhs.leading_ts.not_clipped;
        self.leading_ca.clipped += rhs.leading_ca.clipped;
        self.leading_ca.not_clipped += rhs.leading_ca.not_clipped;
        self.trailing_cs.clipped += rhs.trailing_cs.clipped;
        self.trailing_cs.not_clipped += rhs.trailing_cs.not_clipped;
        self.trailing_ts.clipped += rhs.trailing_ts.clipped;
        self.trailing_ts.not_clipped += rhs.trailing_ts.not_clipped;
        self.trailing_tg.clipped += rhs.trailing_tg.clipped;
        self.trailing_tg.not_clipped += rhs.trailing_tg.not_clipped;
    }
}

fn count_segment(store: &RecordStore<()>) -> Output {
    let mut out = Output::default();

    for idx in 0..store.len() as u32 {
        let rec = store.record(idx);

        // Skip reads that can't contribute signal.
        if rec.flags.is_unmapped()
            || rec.flags.is_secondary()
            || rec.flags.is_supplementary()
            || rec.flags.is_failed_qc()
            || rec.flags.is_duplicate()
        {
            continue;
        }

        // State machine: track the last Match and the last SoftClip event.
        // Only Match → Match and soft-clip ↔ Match adjacencies matter;
        // any other event breaks the chain.
        let mut last_match: Option<(Base, Pos0)> = None;
        let mut last_softclip_end: Option<Base> = None;

        let walk = match rec.aligned_pairs_with_read(store) {
            Ok(w) => w.with_soft_clips(),
            Err(_) => continue,
        };

        for ev in walk {
            match ev {
                AlignedPairWithRead::Match { rpos, query, .. } => {
                    // ── Leading soft-clip case ──
                    // SoftClip(last=C|T) → Match(first-base=G|A)
                    if let Some(sc_base) = last_softclip_end.take() {
                        match (sc_base, query) {
                            (Base::C, Base::G) => out.leading_cs.clipped += 1,
                            (Base::T, Base::G) => out.leading_ts.clipped += 1,
                            (Base::C, Base::A) => out.leading_ca.clipped += 1,
                            _ => {}
                        }
                    }

                    // ── Consecutive Match: aligned CpG (control) ──
                    if let Some((prev_base, prev_rpos)) = last_match {
                        if prev_base == Base::C
                            && query == Base::G
                            && rpos.as_u64() == prev_rpos.as_u64() + 1
                        {
                            out.leading_cs.not_clipped += 1;
                            out.leading_ts.not_clipped += 1;
                            out.leading_ca.not_clipped += 1;
                            out.trailing_cs.not_clipped += 1;
                            out.trailing_ts.not_clipped += 1;
                            out.trailing_tg.not_clipped += 1;
                        }
                    }

                    last_match = Some((query, rpos));
                    last_softclip_end = None;
                }
                AlignedPairWithRead::SoftClip { query, .. } => {
                    // ── Trailing soft-clip case ──
                    // Match(last-base=C|T) → SoftClip(first=G|A)
                    if let Some((prev_base, _prev_rpos)) = last_match.take() {
                        match (prev_base, query.first().copied()) {
                            (Base::C, Some(Base::G)) => out.trailing_cs.clipped += 1,
                            (Base::C, Some(Base::A)) => out.trailing_ts.clipped += 1,
                            (Base::T, Some(Base::G)) => out.trailing_tg.clipped += 1,
                            _ => {}
                        }
                    }
                    last_softclip_end = query.last().copied();
                }
                _ => {
                    // Insertion, Deletion, RefSkip, Padding, Unknown —
                    // break adjacency on either side.
                    last_match = None;
                    last_softclip_end = None;
                }
            }
        }
    }

    out
}

fn print_results(total: &Output) {
    fn rate(c: &Count) -> f64 {
        let denom = c.clipped + c.not_clipped;
        if denom == 0 { 0.0 } else { c.clipped as f64 / denom as f64 }
    }

    println!(
        "leading  C→G  (soft-clip C,  aligned G): {:>8} / {:>8}  ({:.4})",
        total.leading_cs.clipped,
        total.leading_cs.not_clipped,
        rate(&total.leading_cs)
    );
    println!(
        "leading  T→G  (soft-clip T,  aligned G): {:>8} / {:>8}  ({:.4})",
        total.leading_ts.clipped,
        total.leading_ts.not_clipped,
        rate(&total.leading_ts)
    );
    println!(
        "trailing C→G  (aligned C,  soft-clip G): {:>8} / {:>8}  ({:.4})",
        total.trailing_cs.clipped,
        total.trailing_cs.not_clipped,
        rate(&total.trailing_cs)
    );
    println!(
        "trailing C→A  (aligned C,  soft-clip A): {:>8} / {:>8}  ({:.4})",
        total.trailing_ts.clipped,
        total.trailing_ts.not_clipped,
        rate(&total.trailing_ts)
    );
    println!(
        "leading  C→A  (soft-clip C,  aligned A): {:>8} / {:>8}  ({:.4})",
        total.leading_ca.clipped,
        total.leading_ca.not_clipped,
        rate(&total.leading_ca)
    );
    println!(
        "trailing T→G  (aligned T,  soft-clip G): {:>8} / {:>8}  ({:.4})",
        total.trailing_tg.clipped,
        total.trailing_tg.not_clipped,
        rate(&total.trailing_tg)
    );
}
