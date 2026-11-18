//! Rewrite a coordinate-sorted BAM as CRAM, and as BAM for comparison, and
//! time both.
//!
//! ```text
//! cargo run --release --example cram_write -- in.bam ref.fa out.cram [threads]
//! ```
#![allow(clippy::print_stdout, reason = "example output")]

use seqair::bam::BamWriterBuilder;
use seqair::bam::stream::BamStreamReader;
use seqair::cram::writer::CramWriterBuilder;
use seqair::fasta::IndexedFastaReader;
use std::path::PathBuf;
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1).map(PathBuf::from);
    let (Some(bam), Some(fasta), Some(cram)) = (args.next(), args.next(), args.next()) else {
        return Err("usage: cram_write <in.bam> <ref.fa> <out.cram> [threads]".into());
    };
    let threads: usize = args.next().map_or(Ok(0), |t| t.to_string_lossy().parse())?;

    let mut reader = BamStreamReader::new(std::fs::File::open(&bam)?)?;
    let records = reader.records().collect::<Result<Vec<_>, _>>()?;
    let header = reader.into_header();

    let start = Instant::now();
    let mut writer = CramWriterBuilder::to_path(&cram, &header, IndexedFastaReader::open(&fasta)?)
        .compression_threads(threads)
        .build()?;
    let built = start.elapsed();
    for record in &records {
        writer.write_raw(record)?;
    }
    writer.finish()?;
    let cram_time = start.elapsed();

    let bam_out = cram.with_extension("rewritten.bam");
    let start = Instant::now();
    let mut writer =
        BamWriterBuilder::to_path(&bam_out, &header).compression_threads(threads).build()?;
    for record in &records {
        writer.write_raw(record)?;
    }
    writer.finish()?;
    let bam_time = start.elapsed();

    let size = |p: &PathBuf| std::fs::metadata(p).map(|m| m.len());
    println!(
        "{} records\nCRAM {:>12} bytes in {cram_time:?} (header and M5s {built:?})\nBAM  {:>12} bytes in {bam_time:?}",
        records.len(),
        size(&cram)?,
        size(&bam_out)?,
    );
    Ok(())
}
