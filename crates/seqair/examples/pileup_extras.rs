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
    region: String,

    /// Minimum mapping quality.
    #[clap(long, default_value_t = 20)]
    min_mapq: u8,
}

/// Per-record data computed once at load time, accessed per-column.
#[derive(Debug)]
struct ReadInfo {
    /// Qname bytes, owned — avoids repeated slab lookups during overlap dedup.
    qname: Vec<u8>,
    /// Read group extracted from aux tags (RG:Z:...), if present.
    read_group: Option<String>,
    /// Fraction of the read that is aligned (`matching_bases` / `seq_len`).
    aligned_fraction: f64,
}

fn main() -> anyhow::Result<()> {
    let args = Cli::parse();

    let mut readers =
        Readers::open(&args.input, &args.reference).context("could not open BAM + FASTA")?;

    let (tid, start, end) = parse_region(&args.region, readers.header())?;
    let contig_name = readers.header().target_name(tid).context("unknown tid")?.to_owned();

    // ── Step 1: Load records into RecordStore<()> ─────────────────────
    let mut store = RecordStore::new();
    readers
        .fetch_into(tid, start, end, &mut store)
        .map_err(anyhow::Error::from)
        .context("could not fetch region")?;

    let record_count = store.len();
    eprintln!(
        "Loaded {record_count} records from {contig_name}:{start}-{end}",
        start = *start,
        end = *end
    );

    // ── Step 2: Compute per-record extras ─────────────────────────────
    //
    // The closure receives (record_index, &RecordStore<()>) and can access
    // any slab: record fields, qname, aux, sequence, quality.
    let store = store.with_extras(|idx, s| {
        let rec = s.record(idx);
        let qname = s.qname(idx).to_vec();
        let read_group = extract_rg(s.aux(idx));
        let aligned_fraction =
            if rec.seq_len > 0 { rec.matching_bases as f64 / rec.seq_len as f64 } else { 0.0 };
        ReadInfo { qname, read_group, aligned_fraction }
    });

    // ── Step 3: Sort/dedup on the typed store ─────────────────────────
    //
    // This is safe because each SlimRecord carries an extras_idx that
    // survives reordering. Before this change, sort_by_pos and dedup
    // were restricted to RecordStore<()>.
    let mut store = store;
    store.sort_by_pos();
    store.dedup();

    let after_dedup = store.len();
    if after_dedup < record_count {
        eprintln!("Deduped: {record_count} → {after_dedup} records");
    }

    // ── Step 4: Build the pileup engine with the typed store ──────────
    let ref_bases = readers
        .fasta_mut()
        .fetch_seq(&contig_name, start, end)
        .context("could not fetch reference")?;
    let ref_seq = RefSeq::new(Rc::from(Base::from_ascii_vec(ref_bases)), start);

    let min_mapq = args.min_mapq;
    let mut engine = PileupEngine::new(store, start, end);
    engine.set_reference_seq(ref_seq);
    engine.set_filter(move |flags, _aux| !flags.is_unmapped() && !flags.is_secondary());

    // ── Step 5: Iterate with columns_with_store ───────────────────────
    //
    // columns_with_store yields (PileupColumn, &RecordStore<ReadInfo>)
    // per position, giving access to extras without pre-extraction.
    println!("pos\tdepth\tref\tuniq_reads\tread_groups\tmean_aligned_frac");

    let mut cols = engine.columns_with_store();
    let mut seen_qnames: Vec<u32> = Vec::new();
    while let Some((column, store)) = cols.next_column() {
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

            // Access per-record extras via the store.
            let info = store.extra(aln.record_idx());

            // Use the pre-computed aligned fraction.
            aligned_sum += info.aligned_fraction;
            counted += 1;

            // Count distinct read names (overlap dedup use case).
            // In a real caller you'd use this to resolve mate-pair overlaps.
            let already_seen =
                seen_qnames.iter().any(|&prev_idx| store.extra(prev_idx).qname == info.qname);
            if !already_seen {
                seen_qnames.push(aln.record_idx());
            }

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

fn parse_region(
    region: &str,
    header: &seqair::bam::BamHeader,
) -> anyhow::Result<(u32, Pos0, Pos0)> {
    let region: seqair_types::RegionString = region.parse().context("invalid region")?;
    let tid = header
        .tid(&region.chromosome)
        .with_context(|| format!("contig '{}' not found", region.chromosome))?;
    let start = *region.start.context("no start")?;
    let end = *region.end.context("no end")?;
    let start = Pos0::new(start).context("invalid start")?;
    let end = Pos0::new(end).context("invalid end")?;
    Ok((tid, start, end))
}
