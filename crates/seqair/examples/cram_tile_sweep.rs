//! Sweep a contig in tiles on several threads, the way a tiled caller fetches
//! a CRAM, and report what the shared decoded-slice cache did.
//!
//! `cargo run --release --example cram_tile_sweep -- <cram> <fasta> <chrom> <tile_len> <overlap> <threads> <budget_mib|auto> [start end]`
//!
//! Each thread forks the reader and takes the next tile in contig order, so
//! neighbouring tiles are fetched at about the same time. A budget of 0 turns
//! the cache off: every tile decodes every slice it touches.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::print_stdout, reason = "example")]
#![allow(clippy::indexing_slicing, clippy::arithmetic_side_effects, reason = "example")]
use seqair::bam::RecordStore;
use seqair::cram::reader::IndexedCramReader;
use seqair::cram::slice_cache::SliceCacheBudget;
use seqair_types::Pos0;
use std::path::Path;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let cram = Path::new(&args[1]);
    let fasta = Path::new(&args[2]);
    let chrom = &args[3];
    let tile_len: u64 = args[4].parse().unwrap();
    let overlap: u64 = args[5].parse().unwrap();
    let threads: usize = args[6].parse().unwrap();
    let budget = match args[7].as_str() {
        "auto" => SliceCacheBudget::Auto,
        mib => SliceCacheBudget::Bytes(mib.parse::<usize>().unwrap() << 20),
    };

    let reader = IndexedCramReader::open(cram, fasta).unwrap();
    reader.set_slice_cache_budget(budget);
    let tid = reader.header().tid(chrom).unwrap();
    let len = reader.header().target_len(tid).unwrap();
    let start: u64 = args.get(8).map_or(0, |s| s.parse().unwrap());
    let end: u64 = args.get(9).map_or(len, |s| s.parse().unwrap()).min(len);
    let tiles: Vec<(u64, u64)> = (start..end)
        .step_by(usize::try_from(tile_len).unwrap())
        .map(|s| (s, (s + tile_len).min(end) - 1))
        .collect();

    let next = AtomicUsize::new(0);
    let records = AtomicU64::new(0);
    let t = std::time::Instant::now();
    std::thread::scope(|scope| {
        for _ in 0..threads {
            let mut reader = reader.fork().unwrap();
            let (next, records, tiles) = (&next, &records, &tiles);
            scope.spawn(move || {
                let mut store = RecordStore::new();
                while let Some(&(s, e)) = tiles.get(next.fetch_add(1, Ordering::Relaxed)) {
                    let lo = Pos0::try_from(s.saturating_sub(overlap)).unwrap();
                    let hi = Pos0::try_from((e + overlap).min(len - 1)).unwrap();
                    let n = reader.fetch_into(tid, (lo..=hi).into(), &mut store).unwrap();
                    records.fetch_add(n as u64, Ordering::Relaxed);
                }
            });
        }
    });
    let elapsed = t.elapsed();
    let stats = reader.slice_cache_stats();
    println!(
        "{chrom}:{start}-{end} tiles={} tile_len={tile_len} overlap={overlap} threads={threads} \
         {budget:?}: {:.2?}, {} records fetched, {} slices decoded, \
         largest slice {:.1} MiB, cache holds {} slices / {:.1} MiB",
        tiles.len(),
        elapsed,
        records.load(Ordering::Relaxed),
        stats.decoded,
        stats.largest_slice_bytes as f64 / f64::from(1 << 20),
        stats.cached_slices,
        stats.cached_bytes as f64 / f64::from(1 << 20),
    );
}
