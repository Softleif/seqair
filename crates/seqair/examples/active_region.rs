#![allow(
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    clippy::print_stdout,
    clippy::unwrap_used,
    reason = "example"
)]

//! Active regions: gather the reads spanning a window around a column of
//! interest, from inside the pileup loop, via
//! [`PileupColumn::records_overlapping`].
//!
//! A two-tier caller runs every column through a fast path and, at a column
//! that trips a trigger — here: any read carries an indel at the column —
//! re-solves a window of `--flank` bases either side against the reads that
//! span it. Those reads are already in the engine's store; the window query
//! hands back their indices without a second fetch, and
//! [`PileupColumn::find_record`] resolves any of them on the column.
//!
//! The hook passed to [`Readers::pileup_with`] sees the same reference the
//! columns will report. It covers the segment, not the reads, so the hook
//! here counts the reads reaching past either end — the ones whose outer
//! bases a reference-aware realignment could not score from this `RefSeq`.
//!
//! Columns inside a region already printed are skipped, so one line comes out
//! per cluster of triggers.
//!
//! [`PileupColumn::records_overlapping`]: seqair::bam::pileup::PileupColumn::records_overlapping
//! [`PileupColumn::find_record`]: seqair::bam::pileup::PileupColumn::find_record
//! [`Readers::pileup_with`]: seqair::reader::Readers::pileup_with

use anyhow::Context;
use clap::Parser as _;
use seqair::{
    bam::RecordIdx,
    bam::pileup::PileupOp,
    reader::{DepthLimit, Readers, SegmentOptions},
};
use seqair_types::{Pos0, RegionString};
use std::{num::NonZeroU32, path::PathBuf};

/// seqair active regions — demonstrates `PileupColumn::records_overlapping`.
#[derive(Debug, clap::Parser)]
struct Cli {
    /// BAM/CRAM file to read (must be indexed).
    input: PathBuf,

    /// FASTA reference (faidx-indexed).
    fasta: PathBuf,

    /// Region to process: `chr`, `chr:start`, or `chr:start-end` (1-based).
    #[clap(long, short)]
    region: RegionString,

    /// Bases either side of a triggering column.
    #[clap(long, default_value_t = 150)]
    flank: u32,

    /// Stop after this many regions.
    #[clap(long, default_value_t = 5)]
    show: usize,
}

fn main() -> anyhow::Result<()> {
    let args = Cli::parse();
    let mut readers = Readers::open(&args.input, &args.fasta).context("could not open inputs")?;

    let contig = args.region.chromosome.as_str();
    let start = args.region.start.map_or(Pos0::ZERO, |p| p.to_zero_based());
    let end = match args.region.end {
        Some(p) => p.to_zero_based(),
        None => {
            let tid = readers.header().tid(contig).context("contig not in header")?;
            Pos0::new(readers.header().target_len(tid).context("no contig length")? as u32 - 1)
                .context("contig length out of range")?
        }
    };
    let opts = SegmentOptions::new(NonZeroU32::new(1_000_000).unwrap());
    let segment = readers
        .segments((contig, start, end), opts)
        .context("planning segment failed")?
        .next()
        .context("region produced no segment")?;

    let mut past_segment = 0usize;
    let mut pileup = readers
        .pileup_with(&segment, DepthLimit::Unlimited, |store, ref_seq| {
            past_segment = store
                .records()
                .filter(|rec| !rec.flags.is_unmapped())
                .filter(|rec| {
                    ref_seq.try_base_at(rec.pos).is_none()
                        || ref_seq.try_base_at(rec.end_pos).is_none()
                })
                .count();
        })
        .context("pileup_with failed")?;
    eprintln!(
        "{past_segment} read(s) reach past the segment; the hook's reference does not hold their outer bases"
    );

    println!("region\ttrigger_pos\tdepth\tspanning\tfully_spanning");
    let mut skip_through: Option<Pos0> = None;
    let mut shown = 0usize;
    while let Some(col) = pileup.pileups() {
        if shown >= args.show {
            break;
        }
        if skip_through.is_some_and(|through| col.pos() <= through) {
            continue;
        }
        let triggered = col.alignments().any(|aln| {
            matches!(
                aln.op(),
                PileupOp::Insertion { .. }
                    | PileupOp::Deletion { .. }
                    | PileupOp::ComplexIndel { .. }
            )
        });
        if !triggered {
            continue;
        }

        let win_start = Pos0::new(col.pos().as_u32().saturating_sub(args.flank))
            .context("window start out of range")?;
        let win_end =
            Pos0::new(col.pos().as_u32() + args.flank).context("window end out of range")?;
        let spanning: Vec<RecordIdx> = col.records_overlapping(win_start, win_end).collect();
        let fully = spanning
            .iter()
            .filter(|&&idx| {
                let rec = col.store().record(idx);
                rec.pos <= win_start && rec.end_pos >= win_end
            })
            .count();
        println!(
            "{contig}:{}-{}\t{}\t{}\t{}\t{fully}",
            win_start.as_u64() + 1,
            win_end.as_u64() + 1,
            col.pos().as_u64() + 1,
            col.depth(),
            spanning.len(),
        );
        skip_through = Some(win_end);
        shown += 1;
    }

    Ok(())
}
