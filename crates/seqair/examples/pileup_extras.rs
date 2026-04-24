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
//! accessing it during pileup iteration via `columns_with_store`.
//!
//! The workflow:
//!
//! 1. Load records into a `RecordStore<()>` (the default).
//! 2. Compute per-record extras with `store.with_extras(...)`, producing
//!    a `RecordStore<ReadInfo>`.
//! 3. Sort/dedup the typed store — `extras_idx` keeps the mapping intact.
//! 4. Build a `PileupEngine` from the typed store.
//! 5. Iterate with `columns_with_store()` to access extras during pileup.
//!
//! This avoids pre-extracting per-record data into a separate `Vec` before
//! iteration — the extras live alongside the record slabs and are accessed
//! by record index through the store.

use anyhow::Context;
use clap::Parser as _;
use seqair::bam::RecordStore;
use seqair::bam::pileup::{PileupEngine, RefSeq};
use seqair::reader::Readers;
use seqair_types::{Base, Pos0};
use std::path::PathBuf;
use std::rc::Rc;

/// seqair pileup-extras — demonstrate per-record extras in pileup
///
/// Loads a BAM region, computes per-record metadata (read group, aligned
/// fraction, overlap-dedup key), and prints a per-column summary showing
/// how extras are accessed during iteration.
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
#[derive(Debug)]
struct ReadInfo {
    /// Read group extracted from aux tags (RG:Z:...), if present.
    read_group: Option<String>,
    /// Fraction of the read that is aligned (`matching_bases` / `seq_len`).
    aligned_fraction: f64,
}

fn main() -> anyhow::Result<()> {
    let args = Cli::parse();

    let store = RecordStore::new().with_extras(|idx, s| {
        let rec = s.record(idx);
        let read_group = extract_rg(s.aux(idx));
        let aligned_fraction =
            if rec.seq_len > 0 { rec.matching_bases as f64 / rec.seq_len as f64 } else { 0.0 };
        ReadInfo { read_group, aligned_fraction }
    });
    // future feature: pre-filter records at load time, e.g. by mapq:
    // let store = store.with_pre_filter(|record| record.mapq >= args.min_mapq);

    // todo: add this to give readers custom store and also take on its generic type, e.g. `Readers<ReadInfo>` — would be more ergonomic than passing the store separately to pileup and columns_with_store
    let mut readers = Readers::open_with_store(&args.input, &args.reference, store)
        .context("could not open BAM + FASTA")?;

    // todo: add this method that resolves region strings directly using headers to fill optional start and end positions, e.g. `chr1` gets start=0 and end=chr1_len, `chr1:1000` gets end=chrom_len
    let (tid, start, end) = args.region.resolve(&readers);

    // todo: build proper pileup engine, should set reference seq, etc
    let mut engine: PileupEngine = readers.pileup(store, start, end);
    engine.set_filter(|flags, _aux| !flags.is_unmapped());

    println!("pos\tdepth\tref\tuniq_reads\tread_groups\tmean_aligned_frac");

    // todo: make this work as iterator
    for column in engine.iter() {
        let depth = column.depth();
        if depth == 0 {
            continue;
        }

        let mut rg_counts: Vec<(&str, u32)> = Vec::new();
        let mut aligned_sum = 0.0;
        let mut counted = 0u32;
        seen_qnames.clear();

        for aln in column.alignments() {
            if aln.mapq < min_mapq {
                continue;
            }

            // Access per-record extras
            let info: ReadInfo = aln.extra();

            // Use the pre-computed aligned fraction.
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
            "{pos}\t{depth}\t{ref_base}\t{uniq}\t{rgs}\t{mean:.3}",
            pos = *column.pos() + 1,
            ref_base = column.reference_base() as u8 as char,
            uniq = seen_qnames.len(),
            rgs = rg_summary.join(","),
            mean = aligned_sum / counted as f64,
        );
    }

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
