#![allow(unused, reason = "example code")]

use clap::Parser as _;
use seqair::{
    Readers,
    bam::{
        RecordStore,
        record_store::{CustomizeRecordStore, SlimRecord},
    },
    reader::SegmentOptions,
};

#[derive(Debug, clap::Parser)]
struct Cli {
    input: std::path::PathBuf,     // BAM file
    reference: std::path::PathBuf, // FASTA file
    #[clap(long)]
    region: seqair_types::RegionString,
    #[clap(long, default_value_t = 20)]
    min_mapq: u8,
}

fn main() -> anyhow::Result<()> {
    let args = Cli::parse();

    let mut readers = Readers::open_customized(
        &args.input,
        &args.reference,
        ReadInfoBuilder { min_mapq: args.min_mapq },
    )?;
    let max_len = std::num::NonZeroU32::new(100_000).expect("non-zero literal");
    let opts = SegmentOptions::new(max_len);
    let plan: Vec<_> = readers.segments(&args.region, opts)?.collect();

    for segment in &plan {
        let mut engine = readers.pileup(segment)?;

        while let Some(column) = engine.pileups() {
            for aln in column.alignments() {
                let info = aln.extra();

                // ...
            }
        }
        // `engine` drops here; the RecordStore goes back into `readers`.
    }

    Ok(())
}

struct ReadInfo {
    read_group: Option<seqair_types::SmolStr>,
}

#[derive(Debug, Clone)]
struct ReadInfoBuilder {
    min_mapq: u8,
}

impl CustomizeRecordStore for ReadInfoBuilder {
    type Extra = ReadInfo;

    fn keep_record(&mut self, rec: &SlimRecord, _: &RecordStore<ReadInfo>) -> bool {
        !rec.flags.is_unmapped()
    }

    fn compute(&mut self, rec: &SlimRecord, store: &RecordStore<ReadInfo>) -> ReadInfo {
        ReadInfo { read_group: rec.aux(store).ok().and_then(|aux| aux.get("RG").ok()) }
    }
}
