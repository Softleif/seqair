#![allow(clippy::print_stdout, clippy::arithmetic_side_effects, reason = "example")]

//! Profiling probe: decode every record of a region into a `RecordStore`,
//! tile by tile, with no pileup — the BAM read path on its own.

use anyhow::Context;
use clap::Parser as _;
use seqair::bam::{IndexedBamReader, Pos0, RecordStore};
use std::path::PathBuf;

#[derive(Debug, clap::Parser)]
struct Cli {
    input: PathBuf,
    #[clap(long)]
    contig: String,
    #[clap(long)]
    start: u32,
    #[clap(long)]
    end: u32,
    #[clap(long, default_value_t = 100_000)]
    tile: u32,
}

fn main() -> anyhow::Result<()> {
    let args = Cli::parse();
    let mut reader = IndexedBamReader::open(&args.input)?;
    let tid = reader.header().tid(&args.contig).context("contig")?;
    let mut store = RecordStore::new();
    let (mut records, mut bases) = (0usize, 0usize);
    let mut start = args.start;
    while start < args.end {
        let last = start.saturating_add(args.tile - 1).min(args.end);
        store.clear();
        let span = (Pos0::new(start).context("pos")?..=Pos0::new(last).context("pos")?).into();
        reader.fetch_into(tid, span, &mut store)?;
        records += store.len();
        bases += store.records().map(|r| r.seq_len as usize).sum::<usize>();
        start = last + 1;
    }
    println!("records={records} bases={bases}");
    Ok(())
}
