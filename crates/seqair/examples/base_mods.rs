#![allow(
    clippy::arithmetic_side_effects,
    clippy::cast_lossless,
    clippy::cast_possible_truncation,
    clippy::collapsible_if,
    clippy::doc_markdown,
    clippy::indexing_slicing,
    clippy::print_stdout,
    reason = "example"
)]

//! Extract and summarize base modifications (5mC methylation, etc.) from
//! a BAM file.
//!
//! Oxford Nanopore and PacBio sequencers encode base-modification calls in
//! the `MM` (modification positions) and `ML` (modification likelihoods)
//! auxiliary tags. This example shows how to:
//!
//! 1. Fetch records for a region with [`IndexedReader`] + [`RecordStore`].
//! 2. Parse `MM`/`ML` tags into a [`BaseModState`] per record.
//! 3. Query per-base modification calls with [`mod_at_qpos`].
//! 4. Aggregate per-position methylation frequency across all reads.
//!
//! For a typical bisulfite or nanopore BAM, this produces a BED-like summary
//! of CpG methylation levels at each reference position.

use anyhow::Context;
use clap::Parser as _;
use seqair::bam::aux::{AuxValue, find_tag};
use seqair::bam::{BaseModState, ModType, Pos0, RecordStore};
use seqair::reader::IndexedReader;
use seqair_types::Base;
use std::collections::BTreeMap;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

/// seqair base-mods — extract base-modification calls from a BAM.
#[derive(Debug, clap::Parser)]
struct Cli {
    /// BAM/SAM/CRAM file (must be indexed).
    input: PathBuf,

    /// Region to inspect (e.g. "chr1:1000-2000"). Required.
    #[clap(long, short)]
    region: String,

    /// Print per-read modification calls (verbose). Without this flag,
    /// only the per-position summary is printed.
    #[clap(long, short)]
    verbose: bool,

    /// Minimum modification probability (0–255 scale) to count as
    /// "modified" in the per-position summary.
    #[clap(long, default_value_t = 128)]
    min_prob: u8,

    /// Only report this modification code (e.g. 'm' for 5mC, 'h' for 5hmC).
    /// Defaults to all modification types.
    #[clap(long)]
    mod_code: Option<char>,

    /// Output file (defaults to stdout).
    #[clap(long, short)]
    output: Option<PathBuf>,
}

/// Accumulator for per-reference-position modification statistics.
struct PosSummary {
    modified: u32,
    unmodified: u32,
    canonical_base: Base,
    mod_code: u8,
}

impl PosSummary {
    fn total(&self) -> u32 {
        self.modified + self.unmodified
    }

    fn frequency(&self) -> f64 {
        if self.total() == 0 {
            return 0.0;
        }
        self.modified as f64 / self.total() as f64
    }
}

fn main() -> anyhow::Result<()> {
    let args = Cli::parse();

    let mut reader = IndexedReader::open(&args.input).context("could not open alignment file")?;
    let (tid, start, end) = parse_region(&args.region, reader.header())?;
    let contig_name = reader.header().target_name(tid).context("unknown tid")?.to_owned();

    let start_pos = Pos0::new(start).context("invalid start position")?;
    let end_pos = Pos0::new(end).context("invalid end position")?;

    let mut store = RecordStore::new();
    reader
        .fetch_into(tid, start_pos, end_pos, &mut store)
        .with_context(|| format!("fetch failed for {contig_name}:{start}-{end}"))?;

    let mut output: Box<dyn Write> = if let Some(ref path) = args.output {
        Box::new(BufWriter::new(
            std::fs::File::create(path).with_context(|| format!("could not create {path:?}"))?,
        ))
    } else {
        Box::new(BufWriter::new(std::io::stdout().lock()))
    };

    // Per-reference-position methylation accumulator.
    // Key: (ref_pos, mod_code) — so we can track different mod types separately.
    let mut pos_stats: BTreeMap<(u32, u8), PosSummary> = BTreeMap::new();

    let mut records_with_mods = 0u32;
    let mut total_calls = 0u32;

    for idx in 0..store.len() as u32 {
        let rec = store.record(idx);

        // Skip unmapped reads.
        if rec.flags.is_unmapped() || rec.tid < 0 {
            continue;
        }

        let aux = store.aux(idx);
        let seq = store.seq(idx);

        // Look up MM and ML aux tags. Both must be present for modification
        // calls. MM is a Z-type string; ML is a B:C (byte array).
        let Some(AuxValue::String(mm)) = find_tag(aux, *b"MM") else {
            continue;
        };
        let Some(AuxValue::ArrayU8(ml)) = find_tag(aux, *b"ML") else {
            continue;
        };

        // Parse MM/ML into resolved per-position calls. The parser handles
        // reverse-strand coordinate flipping internally.
        let state = match BaseModState::parse(mm, ml, seq, rec.flags.is_reverse()) {
            Ok(s) => s,
            Err(e) => {
                let qname = std::str::from_utf8(store.qname(idx)).unwrap_or("<non-utf8>");
                eprintln!("warning: skipping {qname}: {e}");
                continue;
            }
        };

        if state.is_empty() {
            continue;
        }
        records_with_mods += 1;

        if args.verbose {
            let qname = std::str::from_utf8(store.qname(idx)).unwrap_or("<non-utf8>");
            writeln!(
                output,
                "# {qname} (pos={}, strand={})",
                *rec.pos,
                strand_char(rec.flags.is_reverse())
            )?;
        }

        // Walk every query position looking for modification calls.
        for qpos in 0..seq.len() {
            let Some(mods) = state.mod_at_qpos(qpos) else {
                continue;
            };

            for m in mods {
                // Filter by modification code if requested.
                let code = mod_code_char(m.mod_type);
                if let Some(filter_code) = args.mod_code {
                    if code != filter_code {
                        continue;
                    }
                }
                total_calls += 1;

                if args.verbose {
                    writeln!(
                        output,
                        "  qpos={qpos:>5}  base={:<1}  mod={code}  prob={:>3}  strand={}",
                        m.canonical_base as u8 as char,
                        m.probability,
                        strand_char(m.strand == seqair::bam::ModStrand::Minus),
                    )?;
                }

                // Map query position back to reference coordinates via the
                // CIGAR, so we can accumulate per-position statistics.
                if let Some(ref_pos) = qpos_to_ref_pos(store.cigar(idx), *rec.pos, qpos as u32) {
                    if ref_pos >= start && ref_pos < end {
                        let code_byte = mod_code_byte(m.mod_type);
                        let entry = pos_stats.entry((ref_pos, code_byte)).or_insert(PosSummary {
                            modified: 0,
                            unmodified: 0,
                            canonical_base: m.canonical_base,
                            mod_code: code_byte,
                        });
                        if m.probability >= args.min_prob {
                            entry.modified += 1;
                        } else {
                            entry.unmodified += 1;
                        }
                    }
                }
            }
        }
    }

    // ── Per-position summary ───────────────────────────────────────────
    if !pos_stats.is_empty() {
        writeln!(output)?;
        writeln!(output, "# Per-position modification summary (min_prob={})", args.min_prob)?;
        writeln!(output, "# chrom\tpos(0-based)\tbase\tmod\tmodified\ttotal\tfrequency")?;
        for ((ref_pos, _code), summary) in &pos_stats {
            writeln!(
                output,
                "{contig_name}\t{ref_pos}\t{}\t{}\t{}\t{}\t{:.3}",
                summary.canonical_base as u8 as char,
                summary.mod_code as char,
                summary.modified,
                summary.total(),
                summary.frequency(),
            )?;
        }
    }

    output.flush()?;

    eprintln!(
        "{contig_name}:{start}-{end}: {} records, {} with modifications, {} total calls",
        store.len(),
        records_with_mods,
        total_calls,
    );

    Ok(())
}

/// Map a query position back to a reference position by walking CIGAR ops.
///
/// Returns `None` if qpos falls in an insertion or soft clip (no reference
/// coordinate).
fn qpos_to_ref_pos(
    cigar: &[seqair::bam::CigarOp],
    record_pos: u32,
    target_qpos: u32,
) -> Option<u32> {
    let mut rpos = record_pos;
    let mut qpos = 0u32;

    for op in cigar {
        let len = op.len();
        match op.op_code() {
            // M, =, X: consume both ref and query
            0 | 7 | 8 => {
                if target_qpos >= qpos && target_qpos < qpos + len {
                    return Some(rpos + (target_qpos - qpos));
                }
                rpos += len;
                qpos += len;
            }
            // I: consume query only — qpos in an insertion has no ref coordinate
            1 => {
                if target_qpos >= qpos && target_qpos < qpos + len {
                    return None;
                }
                qpos += len;
            }
            // D, N: consume ref only
            2 | 3 => {
                rpos += len;
            }
            // S: consume query only (soft clip)
            4 => {
                if target_qpos >= qpos && target_qpos < qpos + len {
                    return None;
                }
                qpos += len;
            }
            // H, P: consume neither
            5 | 6 => {}
            _ => {}
        }
    }
    None
}

fn strand_char(is_reverse: bool) -> char {
    if is_reverse { '-' } else { '+' }
}

fn mod_code_char(mod_type: ModType) -> char {
    match mod_type {
        ModType::Code(c) => c as char,
        ModType::ChEBI(id) => {
            // Common ChEBI ids
            match id {
                27551 => 'm', // 5mC
                76792 => 'h', // 5hmC
                _ => '?',
            }
        }
    }
}

fn mod_code_byte(mod_type: ModType) -> u8 {
    match mod_type {
        ModType::Code(c) => c,
        ModType::ChEBI(27551) => b'm',
        ModType::ChEBI(76792) => b'h',
        ModType::ChEBI(_) => b'?',
    }
}

fn parse_region(region: &str, header: &seqair::bam::BamHeader) -> anyhow::Result<(u32, u32, u32)> {
    if let Some((name, range)) = region.split_once(':') {
        let tid = header.tid(name).with_context(|| format!("contig '{name}' not found"))?;
        let (start_str, end_str) =
            range.split_once('-').with_context(|| format!("invalid range: {range}"))?;
        let start: u32 = start_str
            .replace(',', "")
            .parse()
            .with_context(|| format!("invalid start: {start_str}"))?;
        let end: u32 =
            end_str.replace(',', "").parse().with_context(|| format!("invalid end: {end_str}"))?;
        Ok((tid, start.saturating_sub(1), end))
    } else {
        let tid = header.tid(region).with_context(|| format!("contig '{region}' not found"))?;
        let len = header.target_len(tid).unwrap_or(0) as u32;
        Ok((tid, 0, len))
    }
}
