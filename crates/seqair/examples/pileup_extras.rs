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
//! 1. Define a `Clone` struct implementing [`RecordStoreExtras`] that computes
//!    one `Extra` value per record.
//! 2. Open the readers with `Readers::open_with_extras(bam, fasta, provider)`.
//! 3. Resolve the region string against the header with `readers.resolve_region`.
//! 4. Call `readers.pileup(tid, start, end)` — this fetches records, runs the
//!    extras provider once per record, fetches the reference sequence, and
//!    returns a `PileupEngine<E::Extra>` ready for iteration.
//! 5. Iterate with `while let Some(col) = engine.pileups()`. Each alignment
//!    view exposes `aln.extra()` alongside the usual fields.

use anyhow::Context;
use clap::Parser as _;
use seqair::bam::record_store::{RecordStore, RecordStoreExtras};
use seqair::reader::Readers;
use seqair_types::RegionString;
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
    read_group: Option<String>,
    /// Fraction of the read that is aligned (`matching_bases` / `seq_len`).
    aligned_fraction: f64,
}

/// Extras provider — just a zero-sized marker that produces `ReadInfo` per record.
///
/// `Clone` so `Readers::fork` can duplicate the provider for multi-threaded use.
#[derive(Debug, Clone, Default)]
struct ReadInfoBuilder;

impl RecordStoreExtras for ReadInfoBuilder {
    type Extra = ReadInfo;

    fn compute(&mut self, idx: u32, store: &RecordStore<()>) -> ReadInfo {
        let rec = store.record(idx);
        let read_group = extract_rg(store.aux(idx));
        let aligned_fraction =
            if rec.seq_len > 0 { rec.matching_bases as f64 / rec.seq_len as f64 } else { 0.0 };
        ReadInfo { read_group, aligned_fraction }
    }
}

fn main() -> anyhow::Result<()> {
    let args = Cli::parse();

    let mut readers =
        Readers::<ReadInfoBuilder>::open_with_extras(&args.input, &args.reference, ReadInfoBuilder)
            .context("could not open BAM + FASTA")?;

    // Resolve region string against the header — fills in default start/end
    // when missing and validates the contig name.
    let (tid, start, end) = readers.resolve_region(&args.region)?;

    let min_mapq = args.min_mapq;
    let mut engine = readers.pileup(tid, start, end).context("could not build pileup")?;
    engine.set_filter(move |flags, _aux| !flags.is_unmapped() && !flags.is_secondary());

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
        rg_counts.sort_by(|a, b| b.1.cmp(&a.1));
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

/// Extract the RG:Z tag from raw BAM aux bytes.
fn extract_rg(aux: &[u8]) -> Option<String> {
    // BAM aux format: 2-byte tag + 1-byte type + value.
    // RG is type 'Z' (NUL-terminated string).
    let mut i = 0;
    while i + 3 <= aux.len() {
        let tag = [aux[i], aux[i + 1]];
        let val_type = aux[i + 2];
        i += 3;
        if tag == *b"RG" && val_type == b'Z' {
            let start = i;
            while i < aux.len() && aux[i] != 0 {
                i += 1;
            }
            return String::from_utf8(aux[start..i].to_vec()).ok();
        }
        // Skip value based on type.
        match val_type {
            b'A' | b'c' | b'C' => i += 1,
            b's' | b'S' => i += 2,
            b'i' | b'I' | b'f' => i += 4,
            b'Z' | b'H' => {
                while i < aux.len() && aux[i] != 0 {
                    i += 1;
                }
                i += 1; // skip NUL
            }
            b'B' => {
                if i + 5 > aux.len() {
                    break;
                }
                let elem_type = aux[i];
                let count = u32::from_le_bytes([aux[i + 1], aux[i + 2], aux[i + 3], aux[i + 4]]);
                i += 5;
                let elem_size = match elem_type {
                    b'c' | b'C' => 1,
                    b's' | b'S' => 2,
                    b'i' | b'I' | b'f' => 4,
                    _ => break,
                };
                i += count as usize * elem_size;
            }
            _ => break,
        }
    }
    None
}
