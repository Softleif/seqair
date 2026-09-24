#![allow(
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    clippy::indexing_slicing,
    clippy::print_stdout,
    reason = "example"
)]

//! Profiling probe: a rastair-shaped pileup loop. Tiles a region, piles up
//! each tile against the reference, and per column counts bases by strand
//! behind a base-quality and MAPQ filter. `--repeat N` runs the whole region
//! N times so short regions still give samply enough samples.

use anyhow::Context;
use clap::Parser as _;
use seqair::Readers;
use seqair::reader::SegmentOptions;
use seqair_types::RegionString;
use std::num::NonZeroU32;
use std::path::PathBuf;

#[derive(Debug, clap::Parser)]
struct Cli {
    input: PathBuf,
    #[clap(long, short = 'f')]
    reference: PathBuf,
    #[clap(long, short)]
    region: RegionString,
    #[clap(long, default_value_t = 10_000)]
    tile: u32,
    #[clap(long, default_value_t = 1)]
    repeat: u32,
}

fn main() -> anyhow::Result<()> {
    let args = Cli::parse();
    let mut readers = Readers::open(&args.input, &args.reference).context("open")?;
    let opts = SegmentOptions::new(NonZeroU32::new(args.tile).context("tile")?);
    let plan: Vec<_> = readers.segments(&args.region, opts)?.collect();

    let mut counts = [[0u64; 2]; 4];
    let (mut columns, mut depth, mut mismatches) = (0u64, 0u64, 0u64);
    for _ in 0..args.repeat {
        for seg in &plan {
            let mut p = readers.pileup(seg, seqair::reader::DepthLimit::Unlimited).run()?;
            while let Some(col) = p.pileups() {
                columns += 1;
                depth += col.depth() as u64;
                let ref_base = col.reference_base();
                for aln in col.alignments() {
                    if aln.mapq < 10 {
                        continue;
                    }
                    let (Some(base), Some(q)) = (aln.base(), aln.qual()) else { continue };
                    if q.get().unwrap_or(0) < 20 {
                        continue;
                    }
                    if let Some(i) = base.known_index() {
                        counts[i][usize::from(aln.flags.is_reverse())] += 1;
                    }
                    if base != ref_base {
                        mismatches += 1;
                    }
                }
            }
        }
    }
    println!("columns={columns} depth={depth} mismatches={mismatches} counts={counts:?}");
    Ok(())
}
