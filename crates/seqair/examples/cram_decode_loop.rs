//! Profiling harness: repeatedly decode a CRAM region (not a benchmark).
//!
//! `cargo run --profile profiling --example cram_decode_loop -- <cram> <fasta> <chrom> <iters> [reopen]`
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::print_stdout, reason = "example")]
#![allow(clippy::indexing_slicing, clippy::arithmetic_side_effects, reason = "example")]
use seqair::bam::RecordStore;
use seqair::cram::reader::IndexedCramReader;
use seqair_types::Pos0;
use std::path::Path;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let cram = Path::new(&args[1]);
    let fasta = Path::new(&args[2]);
    let chrom = &args[3];
    let iters: usize = args[4].parse().unwrap();
    let reopen = args.get(5).is_some_and(|s| s == "reopen");
    let mut reader = IndexedCramReader::open(cram, fasta).unwrap();
    let tid = reader.header().tid(chrom).unwrap();
    let len = reader.header().target_len(tid).unwrap();
    let span = (Pos0::new(0).unwrap()..=Pos0::new(u32::try_from(len).unwrap()).unwrap()).into();
    let mut store = RecordStore::new();
    let t = std::time::Instant::now();
    let mut n = 0;
    for _ in 0..iters {
        if reopen {
            reader = IndexedCramReader::open(cram, fasta).unwrap();
        }
        store.clear();
        reader.fetch_into(tid, span, &mut store).unwrap();
        n += store.len();
    }
    println!("{n} records, {:?}/iter", t.elapsed() / u32::try_from(iters).unwrap());
}
