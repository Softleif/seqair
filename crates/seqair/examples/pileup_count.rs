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
//!
//! `--stream` reads the whole input front to back instead (`-` for stdin),
//! with no index, through `StreamPileup`; `--tile` is then the window length
//! and `--region` is ignored.

use anyhow::Context;
use clap::Parser as _;
use seqair::Readers;
use seqair::bam::pileup::PileupColumn;
use seqair::fasta::IndexedFastaReader;
use seqair::reader::{AlignmentStream, DepthLimit, SegmentOptions, StreamPileup};
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
    #[clap(long)]
    stream: bool,
}

#[derive(Debug, Default)]
struct Counts {
    columns: u64,
    depth: u64,
    mismatches: u64,
    by_base: [[u64; 2]; 4],
}

impl Counts {
    fn add(&mut self, col: &PileupColumn<'_, ()>) {
        self.columns += 1;
        self.depth += col.depth() as u64;
        let ref_base = col.reference_base();
        for aln in col.alignments() {
            if aln.mapq.get().is_some_and(|q| q < 10) {
                continue;
            }
            let (Some(base), Some(q)) = (aln.base(), aln.qual()) else { continue };
            if q.get().unwrap_or(0) < 20 {
                continue;
            }
            if let Some(i) = base.known_index() {
                self.by_base[i][usize::from(aln.flags.is_reverse())] += 1;
            }
            if base != ref_base {
                self.mismatches += 1;
            }
        }
    }
}

fn main() -> anyhow::Result<()> {
    let args = Cli::parse();
    let tile = NonZeroU32::new(args.tile).context("tile")?;
    let mut counts = Counts::default();
    if args.stream {
        let input: Box<dyn std::io::Read + Send> = if args.input.as_os_str() == "-" {
            Box::new(std::io::stdin())
        } else {
            Box::new(std::fs::File::open(&args.input).context("open")?)
        };
        let stream = AlignmentStream::from_reader(input).context("open stream")?;
        let fasta = IndexedFastaReader::open(&args.reference).context("open reference")?;
        let mut p = StreamPileup::new(stream, (), DepthLimit::Unlimited)
            .with_window_len(tile)
            .with_reference(fasta);
        while let Some((_, col)) = p.pileups()? {
            counts.add(&col);
        }
    } else {
        let mut readers = Readers::open(&args.input, &args.reference).context("open")?;
        let plan: Vec<_> = readers.segments(&args.region, SegmentOptions::new(tile))?.collect();
        for _ in 0..args.repeat {
            for seg in &plan {
                let mut p = readers.pileup(seg, DepthLimit::Unlimited).run()?;
                while let Some(col) = p.pileups() {
                    counts.add(&col);
                }
            }
        }
    }
    let Counts { columns, depth, mismatches, by_base } = counts;
    println!("columns={columns} depth={depth} mismatches={mismatches} counts={by_base:?}");
    Ok(())
}
