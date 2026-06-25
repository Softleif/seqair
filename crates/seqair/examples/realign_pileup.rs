#![allow(
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    clippy::indexing_slicing,
    clippy::print_stdout,
    clippy::unwrap_used,
    reason = "example"
)]

//! Realign-then-pileup in a single call via [`Readers::pileup_with`].
//!
//! This is the ergonomic counterpart to `realignment.rs`: instead of fetching
//! into a manual [`RecordStore`], rewriting alignments, and building a
//! [`PileupEngine`] by hand, [`Readers::pileup_with`] runs a caller-supplied
//! mutator on the freshly fetched store, re-sorts it, and hands back a pileup
//! guard — all while preserving the internal buffer reuse of
//! [`Readers::pileup`]. This is the exact hook a methylation/variant caller
//! uses to splice in local realignment (e.g. a POA consensus) just before
//! piling up.
//!
//! The "realignment rule" here is deliberately trivial — soft-clip the first
//! base of every read whose CIGAR starts with an `M` op of length ≥ 2, shifting
//! `pos` right by one — so the example stays focused on the *workflow*, not the
//! alignment math. It prints the per-column depth before and after so the
//! effect of the rewrite is visible.

use anyhow::Context;
use clap::Parser as _;
use seqair::{
    bam::{CigarOp, RecordStore, cigar::CigarOpType},
    reader::{DepthLimit, Readers, SegmentOptions},
};
use seqair_types::{Pos0, RegionString};
use std::{num::NonZeroU32, path::PathBuf};

/// seqair realign-then-pileup — demonstrates `Readers::pileup_with`.
#[derive(Debug, clap::Parser)]
struct Cli {
    /// BAM/CRAM file to read (must be indexed).
    input: PathBuf,

    /// FASTA reference (faidx-indexed).
    fasta: PathBuf,

    /// Region to process: `chr`, `chr:start`, or `chr:start-end` (1-based).
    #[clap(long, short)]
    region: RegionString,

    /// Print this many columns of the before/after depth profile.
    #[clap(long, default_value_t = 12)]
    show: usize,
}

/// Soft-clip the first aligned base of every read whose CIGAR starts with an
/// `M` op of length ≥ 2, shifting `pos` right by one. Returns the number of
/// records rewritten.
fn realign_leading_clip(store: &mut RecordStore<()>) -> usize {
    let mut plan: Vec<(u32, Pos0, Vec<CigarOp>)> = Vec::new();
    for idx in 0..store.len() as u32 {
        let rec = store.record(idx);
        if rec.flags.is_unmapped() {
            continue;
        }
        let Ok(cigar) = rec.cigar(store) else { continue };
        let Some(first) = cigar.first() else { continue };
        if first.op_type() != CigarOpType::Match || first.len() < 2 {
            continue;
        }
        let mut new = Vec::with_capacity(cigar.len() + 1);
        new.push(CigarOp::new(CigarOpType::SoftClip, 1));
        new.push(CigarOp::new(CigarOpType::Match, first.len() - 1));
        new.extend_from_slice(&cigar[1..]);
        let Some(new_pos) = Pos0::new(rec.pos.as_u64() as u32 + 1) else { continue };
        plan.push((idx, new_pos, new));
    }
    for (idx, pos, cig) in &plan {
        // Query length is preserved (one M became one S), so this never fails.
        store.set_alignment(*idx, *pos, cig).expect("query length preserved");
    }
    plan.len()
}

fn depth_profile(
    readers: &mut Readers,
    segment: &seqair::reader::Segment,
    realign: bool,
) -> anyhow::Result<Vec<(u64, usize)>> {
    let mut profile = Vec::new();
    let mut count = 0usize;
    if realign {
        let mut guard = readers
            .pileup_with(segment, DepthLimit::Unlimited, |store| {
                count = realign_leading_clip(store);
            })
            .context("pileup_with failed")?;
        while let Some(col) = guard.pileups() {
            profile.push((col.pos().as_u64(), col.depth()));
        }
        eprintln!("realigned {count} record(s)");
    } else {
        let mut guard = readers.pileup(segment, DepthLimit::Unlimited).context("pileup failed")?;
        while let Some(col) = guard.pileups() {
            profile.push((col.pos().as_u64(), col.depth()));
        }
    }
    Ok(profile)
}

fn main() -> anyhow::Result<()> {
    let args = Cli::parse();
    let mut readers = Readers::open(&args.input, &args.fasta).context("could not open inputs")?;

    let start = args.region.start.map_or(Pos0::ZERO, |p| p.to_zero_based());
    let end = match args.region.end {
        Some(p) => p.to_zero_based(),
        None => {
            let tid = readers
                .header()
                .tid(args.region.chromosome.as_str())
                .context("contig not in header")?;
            Pos0::new(readers.header().target_len(tid).context("no contig length")? as u32 - 1)
                .context("contig length out of range")?
        }
    };

    let opts = SegmentOptions::new(NonZeroU32::new(1_000_000).unwrap());
    let segment = readers
        .segments((args.region.chromosome.as_str(), start, end), opts)
        .context("planning segment failed")?
        .next()
        .context("region produced no segment")?;

    let before = depth_profile(&mut readers, &segment, false)?;
    let after = depth_profile(&mut readers, &segment, true)?;

    println!("pos\tdepth_before\tdepth_after");
    for (b, a) in before.iter().zip(after.iter()).take(args.show) {
        let marker = if b.1 == a.1 { "" } else { "  <-- changed" };
        println!("{}\t{}\t{}{marker}", b.0 + 1, b.1, a.1);
    }

    Ok(())
}
