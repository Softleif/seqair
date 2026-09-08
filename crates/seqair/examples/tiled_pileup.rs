//! Replay a variant caller's access pattern on a real BAM and time each phase.
//!
//! A caller like rastair does not query a chromosome once — it tiles it into
//! ~10 kb segments and issues an independent index query per segment. Anything
//! paid per query (BAI lookup, `RegionBuf` setup, re-inflating the block a tile
//! boundary lands in, sorting and mate-linking the store) is therefore paid
//! tens of thousands of times per chromosome, and a benchmark that reads one
//! large region cannot see any of it.
//!
//! ```text
//! cargo run --release --example tiled_pileup -- \
//!     NA12878_chr12.bam --region chr12 --tile 10000
//! ```
#![allow(
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    clippy::print_stdout,
    reason = "example"
)]

use anyhow::Context as _;
use clap::Parser as _;
use seqair::bam::{IndexedBamReader, RecordStore};
use seqair::bam::pileup::PileupEngine;
use seqair_types::{Offset, Pos0, RegionString};
use std::path::PathBuf;
use std::time::{Duration, Instant};

#[derive(Debug, clap::Parser)]
struct Cli {
    /// Indexed BAM to read.
    input: PathBuf,

    /// Region to tile (1-based, inclusive). Defaults to the whole first contig.
    #[clap(long, short)]
    region: Option<RegionString>,

    /// Tile width in bases — a caller's `--segment-max-length`.
    #[clap(long, default_value_t = 10_000)]
    tile: u32,

    /// Bases each tile shares with its neighbour, as a caller's overlap does.
    #[clap(long, default_value_t = 200)]
    overlap: u32,

    /// Stop after this many tiles (0 = no limit), to keep a probe run short.
    #[clap(long, default_value_t = 0)]
    max_tiles: u32,

    /// Only produce columns; don't walk their alignments. Splits the engine's
    /// own cost from what a consumer pays per emitted base.
    #[clap(long)]
    columns_only: bool,
}

#[derive(Default)]
struct Phases {
    fetch: Duration,
    prepare: Duration,
    columns: Duration,
}

fn main() -> anyhow::Result<()> {
    let args = Cli::parse();
    let mut reader =
        IndexedBamReader::open(&args.input).context("could not open the alignment file")?;

    let (contig, start, end) = resolve_region(&reader, args.region.as_ref())?;
    let tid = reader
        .header()
        .tid(&contig)
        .with_context(|| format!("{contig} is not in the header"))?;

    let mut store = RecordStore::<()>::new();
    let mut phases = Phases::default();
    let (mut tiles, mut records, mut cols, mut bases) = (0u64, 0u64, 0u64, 0u64);

    let wall = Instant::now();
    let mut tile_start = start;
    loop {
        if tile_start > end || (args.max_tiles > 0 && tiles >= u64::from(args.max_tiles)) {
            break;
        }
        let span = i64::from(args.tile).saturating_add(i64::from(args.overlap));
        let tile_end =
            tile_start.checked_add_offset(Offset::new(span)).unwrap_or(Pos0::max_value()).min(end);

        let t = Instant::now();
        records += reader.fetch_into(tid, tile_start, tile_end, &mut store)? as u64;
        phases.fetch += t.elapsed();

        let t = Instant::now();
        let input = store.prepare_for_pileup().input;
        phases.prepare += t.elapsed();

        let t = Instant::now();
        let mut engine = PileupEngine::new(input, tile_start, tile_end);
        while let Some(col) = engine.pileups() {
            cols += 1;
            if args.columns_only {
                bases += col.depth() as u64;
            } else {
                for aln in col.alignments() {
                    bases += u64::from(aln.base().is_some());
                }
            }
        }
        phases.columns += t.elapsed();
        store = engine.reclaim_allocation().unwrap_or_default();

        tiles += 1;
        let Some(next) = tile_start.checked_add_offset(Offset::new(i64::from(args.tile))) else {
            break;
        };
        tile_start = next;
    }
    let wall = wall.elapsed();

    println!("tiles     {tiles}");
    println!("records   {records}");
    println!("columns   {cols}");
    println!("bases     {bases}");
    println!("wall      {:>9.3?}", wall);
    println!("  fetch   {:>9.3?}  {:5.1}%", phases.fetch, pct(phases.fetch, wall));
    println!("  prepare {:>9.3?}  {:5.1}%", phases.prepare, pct(phases.prepare, wall));
    println!("  columns {:>9.3?}  {:5.1}%", phases.columns, pct(phases.columns, wall));
    Ok(())
}

fn pct(part: Duration, whole: Duration) -> f64 {
    if whole.is_zero() { 0.0 } else { part.as_secs_f64() / whole.as_secs_f64() * 100.0 }
}

fn resolve_region(
    reader: &IndexedBamReader<std::fs::File>,
    region: Option<&RegionString>,
) -> anyhow::Result<(String, Pos0, Pos0)> {
    let header = reader.header();
    let contig = match region {
        Some(r) => r.chromosome.to_string(),
        None => header.target_names().next().context("header has no contigs")?.to_string(),
    };
    let len = header
        .tid(&contig)
        .and_then(|tid| header.target_len(tid))
        .with_context(|| format!("{contig} is not in the header"))?;
    let last = Pos0::try_from(len.saturating_sub(1)).context("contig longer than Pos0 allows")?;
    let (start, end) = match region {
        Some(r) => (
            r.start.map_or(Pos0::ZERO, |p| Pos0::new(p.as_u32().saturating_sub(1)).unwrap_or(Pos0::ZERO)),
            r.end.map_or(last, |p| Pos0::new(p.as_u32().saturating_sub(1)).unwrap_or(last)).min(last),
        ),
        None => (Pos0::ZERO, last),
    };
    Ok((contig, start, end))
}
