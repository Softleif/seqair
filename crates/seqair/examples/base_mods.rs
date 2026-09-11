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

//! Extract and summarise base modifications (5mC, 5hmC, ...) from a BAM
//! file.
//!
//! Oxford Nanopore and PacBio sequencers encode base-modification calls in
//! the `MM` (modification positions) and `ML` (modification likelihoods)
//! auxiliary tags. This example shows how to:
//!
//! 1. Fetch records for a region with [`IndexedReader`] + [`RecordStore`].
//! 2. Parse `MM`/`ML` tags into a [`BaseModState`] per record.
//! 3. Walk the record's CIGAR with [`SlimRecord::aligned_pairs`](seqair::bam::record_store::SlimRecord::aligned_pairs) and join
//!    each `Match { qpos, rpos }` event with `state.mod_at_qpos(qpos)` to
//!    project per-base modification calls onto the reference.
//! 4. Aggregate per-position methylation frequency across all reads.
//!
//! For a typical bisulfite or nanopore BAM, this produces a BED-like summary
//! of CpG methylation levels at each reference position.

use anyhow::Context;
use clap::Parser as _;
use seqair::bam::{AlignedPair, BaseModState, ModType, Pos0, RecordStore};
use seqair::reader::{IndexedReader, ResolveTid};
use seqair_types::{Base, RegionString};
use std::collections::BTreeMap;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

/// seqair base-mods — extract base-modification calls from a BAM.
#[derive(Debug, clap::Parser)]
struct Cli {
    /// BAM/SAM file (must be indexed). CRAM is not supported here because
    /// this example does not load a FASTA.
    input: PathBuf,

    /// Region to inspect. `chr`, `chr:start`, or `chr:start-end`
    /// (1-based, inclusive).
    #[clap(long, short)]
    region: RegionString,

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

    // Resolve the region against the BAM header.
    let tid = args.region.chromosome.as_str().resolve_tid(reader.header())?.as_u32();
    let contig_name = reader.header().target_name(tid).context("unknown tid")?.to_owned();
    let contig_len =
        reader.header().target_len(tid).context("contig has no length in header")? as u32;
    let start = match args.region.start {
        Some(p) => p.to_zero_based(),
        None => Pos0::new(0).expect("0 is valid"),
    };
    let end = match args.region.end {
        Some(p) => p.to_zero_based(),
        None => Pos0::new(contig_len).context("contig length out of i32 range")?,
    };

    let mut store = RecordStore::new();
    reader.fetch_into(tid, start, end, &mut store).with_context(|| {
        format!("fetch failed for {contig_name}:{}-{}", start.as_u32(), end.as_u32())
    })?;

    let mut output: Box<dyn Write> = if let Some(ref path) = args.output {
        Box::new(BufWriter::new(
            std::fs::File::create(path).with_context(|| format!("could not create {path:?}"))?,
        ))
    } else {
        Box::new(BufWriter::new(std::io::stdout().lock()))
    };

    // Per-reference-position methylation accumulator.
    // Key: (ref_pos, mod_code) — so we track different mod types separately.
    let mut pos_stats: BTreeMap<(u32, u8), PosSummary> = BTreeMap::new();

    let mut records_with_mods = 0u32;
    let mut total_calls = 0u32;

    for idx in store.indices() {
        let Some(rec) = store.record(idx) else { continue };

        // Skip unmapped reads.
        if rec.flags.is_unmapped() || rec.tid < 0 {
            continue;
        }

        // Pull MM/ML/seq/is_reverse out of the record in one call. Returns
        // `Ok(None)` when the record carries no MM tag (most reads in a
        // mixed BAM), surfaces malformed payloads as typed errors.
        let state = match BaseModState::from_record(rec.slim(), &store) {
            Ok(Some(s)) => s,
            Ok(None) => continue,
            Err(e) => {
                let qname = std::str::from_utf8(rec.qname()).unwrap_or("<non-utf8>");
                eprintln!("warning: skipping {qname}: {e}");
                continue;
            }
        };

        if state.is_empty() {
            continue;
        }
        records_with_mods += 1;

        if args.verbose {
            let qname = std::str::from_utf8(rec.qname()).unwrap_or("<non-utf8>");
            writeln!(
                output,
                "# {qname} (pos={}, strand={})",
                rec.pos.as_u32(),
                strand_char(rec.flags.is_reverse())
            )?;
        }

        // Walk CIGAR once. Each `Match { qpos, rpos }` event gives both the
        // query and reference position; lift the modifications at that qpos
        // straight onto the reference. `aligned_pairs` skips Insertion and
        // SoftClip qpos for us — exactly the positions that have no ref
        // coordinate.
        for pair in rec.aligned_pairs(&store)? {
            let AlignedPair::Match { qpos, rpos, .. } = pair else { continue };
            let Some(mods) = state.mod_at_qpos(qpos) else { continue };

            for m in mods {
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
                        "  qpos={qpos:>5}  base={}  mod={code}  prob={:>3}  strand={}",
                        m.canonical_base.as_char(),
                        m.probability,
                        strand_char(m.strand == seqair::bam::ModStrand::Minus),
                    )?;
                }

                let ref_pos = rpos.as_u32();
                if ref_pos < start.as_u32() || ref_pos >= end.as_u32() {
                    continue;
                }

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

    // ── Per-position summary ───────────────────────────────────────────
    if !pos_stats.is_empty() {
        writeln!(output)?;
        writeln!(output, "# Per-position modification summary (min_prob={})", args.min_prob)?;
        writeln!(output, "# chrom\tpos(0-based)\tbase\tmod\tmodified\ttotal\tfrequency")?;
        for ((ref_pos, _code), summary) in &pos_stats {
            writeln!(
                output,
                "{contig_name}\t{ref_pos}\t{}\t{}\t{}\t{}\t{:.3}",
                summary.canonical_base.as_char(),
                summary.mod_code as char,
                summary.modified,
                summary.total(),
                summary.frequency(),
            )?;
        }
    }

    output.flush()?;

    eprintln!(
        "{contig_name}:{}-{}: {} records, {} with modifications, {} total calls",
        start.as_u32(),
        end.as_u32(),
        store.len(),
        records_with_mods,
        total_calls,
    );

    Ok(())
}

fn strand_char(is_reverse: bool) -> char {
    if is_reverse { '-' } else { '+' }
}

fn mod_code_char(mod_type: ModType) -> char {
    match mod_type {
        ModType::Code(c) => c as char,
        ModType::ChEBI(id) => match id {
            27551 => 'm', // 5mC
            76792 => 'h', // 5hmC
            _ => '?',
        },
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
