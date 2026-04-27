#![allow(
    clippy::arithmetic_side_effects,
    clippy::cast_lossless,
    clippy::cast_possible_truncation,
    clippy::indexing_slicing,
    clippy::print_stdout,
    reason = "example"
)]

//! Per-record extras in a pileup.
//!
//! Demonstrates attaching custom per-record data to a `RecordStore` and
//! accessing it during pileup iteration via `Readers<E>` + `aln.extra()`.
//!
//! The workflow:
//!
//! 1. Define a `Clone` struct implementing [`CustomizeRecordStore`] that
//!    computes one `Extra` value per record (and optionally filters records
//!    via `keep_record`).
//! 2. Open the readers with `Readers::open_customized(bam, fasta, customize)`.
//! 3. Resolve the region string against the header with `readers.resolve_region`.
//! 4. Call `readers.pileup(tid, start, end)` — this fetches records (running
//!    `keep_record` at push time), then `compute` once per kept record,
//!    fetches the reference sequence, and returns a `PileupEngine<E::Extra>`
//!    ready for iteration.
//! 5. Iterate with `while let Some(col) = engine.pileups()`. Each alignment
//!    view exposes `aln.extra()` alongside the usual fields.

use anyhow::Context;
use clap::Parser as _;
use seqair::bam::record_store::{CustomizeRecordStore, RecordStore, SlimRecord};
use seqair::reader::Readers;
use seqair_types::RegionString;
use seqair_types::SmolStr;
use std::path::PathBuf;

/// seqair pileup-extras — demonstrate per-record extras in pileup
///
/// Loads a BAM region, computes per-record metadata (read group, aligned
/// fraction), and prints a per-column summary showing how extras are accessed
/// during iteration.
#[derive(Debug, clap::Parser)]
struct Cli {
    /// BAM/SAM/CRAM file (must be indexed).
    input: PathBuf,

    /// Reference FASTA (must be indexed).
    reference: PathBuf,

    /// Region to display (e.g. "chr1:1000-2000").
    #[clap(long, short)]
    region: RegionString,

    /// Minimum mapping quality.
    #[clap(long, default_value_t = 20)]
    min_mapq: u8,
}

/// Per-record data computed once at load time, accessed per-column.
#[derive(Debug, Clone)]
struct ReadInfo {
    /// Read group extracted from aux tags (RG:Z:...), if present.
    read_group: Option<SmolStr>,
    /// Fraction of the read that is aligned (`matching_bases` / `seq_len`).
    aligned_fraction: f64,
}

/// Extras customizer — a zero-sized marker that produces `ReadInfo` per record
/// and (optionally) filters out unmapped/secondary reads at push time.
///
/// `Clone` so `Readers::fork` can duplicate it for multi-threaded use.
#[derive(Debug, Clone, Default)]
struct ReadInfoBuilder;

impl CustomizeRecordStore for ReadInfoBuilder {
    type Extra = ReadInfo;

    fn keep_record(&mut self, rec: &SlimRecord, _: &RecordStore<ReadInfo>) -> bool {
        // Drop unmapped and secondary alignments before they enter the store —
        // saves slab space and skips compute() for records the pileup ignores.
        !rec.flags.is_unmapped() && !rec.flags.is_secondary()
    }

    fn compute(&mut self, rec: &SlimRecord, store: &RecordStore<ReadInfo>) -> ReadInfo {
        let read_group: Option<SmolStr> =
            if let Ok(aux) = rec.aux(store) { aux.get("RG").ok() } else { None };
        let aligned_fraction =
            if rec.seq_len > 0 { rec.matching_bases as f64 / rec.seq_len as f64 } else { 0.0 };
        ReadInfo { read_group, aligned_fraction }
    }
}

fn main() -> anyhow::Result<()> {
    let args = Cli::parse();

    let mut readers = Readers::open_customized(&args.input, &args.reference, ReadInfoBuilder)
        .context("could not open BAM + FASTA")?;

    // Resolve region string against the header — fills in default start/end
    // when missing and validates the contig name.
    let (tid, start, end) = readers.resolve_region(&args.region)?;

    // Unmapped/secondary records are already dropped at push time by
    // `ReadInfoBuilder::keep_record`, so the pileup engine never sees them.
    let min_mapq = args.min_mapq;
    let mut engine = readers.pileup(tid, start, end).context("could not build pileup")?;

    eprintln!(
        "Loaded pileup for {region} ({start}..={end})",
        region = args.region.chromosome,
        start = *start,
        end = *end,
    );

    println!("pos\tdepth\tref\tread_groups\tmean_aligned_frac");

    while let Some(column) = engine.pileups() {
        let depth = column.depth();
        if depth == 0 {
            continue;
        }

        let mut rg_counts: Vec<(&str, u32)> = Vec::new();
        let mut aligned_sum = 0.0;
        let mut counted = 0u32;

        for aln in column.alignments() {
            if aln.mapq < min_mapq {
                continue;
            }

            // Access per-record extras via the alignment view.
            let info = aln.extra();

            aligned_sum += info.aligned_fraction;
            counted += 1;

            // Tally read groups.
            let rg = info.read_group.as_deref().unwrap_or(".");
            if let Some(entry) = rg_counts.iter_mut().find(|(name, _)| *name == rg) {
                entry.1 += 1;
            } else {
                rg_counts.push((rg, 1));
            }
        }

        if counted == 0 {
            continue;
        }

        // Format read group summary: "RG1:5,RG2:3" or ".:10"
        rg_counts.sort_by_key(|b| std::cmp::Reverse(b.1));
        let rg_summary: Vec<String> =
            rg_counts.iter().map(|(name, count)| format!("{name}:{count}")).collect();

        println!(
            "{pos}\t{depth}\t{ref_base}\t{rgs}\t{mean:.3}",
            pos = *column.pos() + 1,
            ref_base = column.reference_base() as u8 as char,
            rgs = rg_summary.join(","),
            mean = aligned_sum / counted as f64,
        );
    }

    readers.recover_store(&mut engine);

    Ok(())
}
