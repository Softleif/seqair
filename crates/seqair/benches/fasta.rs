//! Criterion benchmarks for FASTA indexed fetch.
//!
//! Compares seqair against htslib and noodles.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "benches")]
#![allow(clippy::cast_possible_truncation, reason = "benches")]
#![allow(clippy::arithmetic_side_effects, reason = "benches")]

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use std::hint::black_box;

#[path = "support/data.rs"]
#[allow(dead_code, reason = "each bench target uses a subset")]
mod data;

const FASTA_PATH: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.fasta.gz");
const CHROM: &str = "chr19";

/// Region sizes swept by the steady-state benches. 1 kb is the case that
/// matters most: the pileup fetches one reference slice per segment, so
/// per-call overhead dominates there.
const SIZES: [u64; 4] = [1_000, 10_000, 100_000, 1_000_000];
const START: u64 = 6_100_000;

fn pos(v: u64) -> seqair::bam::Pos<seqair::bam::Zero> {
    seqair::bam::Pos::<seqair::bam::Zero>::new(v as u32).unwrap()
}

fn region(chrom: &str, start: u64, end: u64) -> noodles::core::Region {
    use noodles::core::{Position, Region};
    Region::new(
        chrom,
        Position::try_from(start as usize + 1).unwrap()..=Position::try_from(end as usize).unwrap(),
    )
}

// ---------------------------------------------------------------------------
// Group 1: open + fetch (the original bench; open cost is inside the loop)
// ---------------------------------------------------------------------------

fn fasta_open_and_fetch(c: &mut Criterion) {
    use noodles::fasta as nfasta;
    use rust_htslib::faidx;

    let mut group = c.benchmark_group("fasta_open_and_fetch");

    const FASTA_START: u64 = 6_100_000;
    const FASTA_END: u64 = 6_200_000;

    group.bench_function("seqair", |b| {
        b.iter(|| {
            let path = std::path::Path::new(FASTA_PATH);
            let mut reader = seqair::fasta::IndexedFastaReader::open(path).unwrap();
            let seq = reader.fetch_seq(CHROM, pos(FASTA_START), pos(FASTA_END)).unwrap();
            black_box(seq.len())
        });
    });

    group.bench_function("htslib", |b| {
        b.iter(|| {
            let reader = faidx::Reader::from_path(FASTA_PATH).unwrap();
            let seq =
                reader.fetch_seq(CHROM, FASTA_START as usize, FASTA_END as usize - 1).unwrap();
            black_box(seq.len())
        });
    });

    group.bench_function("noodles", |b| {
        b.iter(|| {
            let path = std::path::Path::new(FASTA_PATH);
            let mut reader =
                nfasta::io::indexed_reader::Builder::default().build_from_path(path).unwrap();
            let record = reader.query(&region(CHROM, FASTA_START, FASTA_END)).unwrap();
            black_box(record.sequence().len())
        });
    });

    group.finish();
}

// ---------------------------------------------------------------------------
// Group 2: open only
// ---------------------------------------------------------------------------

fn fasta_open(c: &mut Criterion) {
    use noodles::fasta as nfasta;
    use rust_htslib::faidx;

    let plain = data::plain_fasta();

    for (label, path) in [("bgzf", std::path::Path::new(FASTA_PATH)), ("plain", plain.as_path())] {
        let mut group = c.benchmark_group(format!("fasta_open/{label}"));

        group.bench_function("seqair", |b| {
            b.iter(|| black_box(seqair::fasta::IndexedFastaReader::open(path).unwrap()));
        });
        group.bench_function("htslib", |b| {
            b.iter(|| black_box(faidx::Reader::from_path(path).unwrap()));
        });
        group.bench_function("noodles", |b| {
            b.iter(|| {
                black_box(
                    nfasta::io::indexed_reader::Builder::default().build_from_path(path).unwrap(),
                )
            });
        });

        group.finish();
    }
}

// ---------------------------------------------------------------------------
// Group 3: steady-state fetch — readers opened once, sweeping region size
// ---------------------------------------------------------------------------

fn fasta_fetch_steady(c: &mut Criterion) {
    use noodles::fasta as nfasta;
    use rust_htslib::faidx;

    let plain = data::plain_fasta();

    for (label, path) in [("bgzf", std::path::Path::new(FASTA_PATH)), ("plain", plain.as_path())] {
        let mut group = c.benchmark_group(format!("fasta_fetch_steady/{label}"));

        let mut seqair_reader = seqair::fasta::IndexedFastaReader::open(path).unwrap();
        let htslib_reader = faidx::Reader::from_path(path).unwrap();
        let mut noodles_reader =
            nfasta::io::indexed_reader::Builder::default().build_from_path(path).unwrap();
        let mut reuse = Vec::new();

        for size in SIZES {
            let end = START + size;
            group.throughput(Throughput::Bytes(size));

            // Allocates a fresh Vec per call.
            group.bench_with_input(BenchmarkId::new("seqair_fetch_seq", size), &size, |b, _| {
                b.iter(|| {
                    black_box(seqair_reader.fetch_seq(CHROM, pos(START), pos(end)).unwrap().len())
                });
            });

            // Reuses the caller's buffer — seqair's zero-alloc path.
            group.bench_with_input(BenchmarkId::new("seqair_fetch_into", size), &size, |b, _| {
                b.iter(|| {
                    seqair_reader.fetch_seq_into(CHROM, pos(START), pos(end), &mut reuse).unwrap();
                    black_box(reuse.len())
                });
            });

            group.bench_with_input(BenchmarkId::new("htslib", size), &size, |b, _| {
                b.iter(|| {
                    black_box(
                        htslib_reader
                            .fetch_seq(CHROM, START as usize, end as usize - 1)
                            .unwrap()
                            .len(),
                    )
                });
            });

            // `query` returns an owned Record — allocates per call, no reuse API.
            group.bench_with_input(BenchmarkId::new("noodles", size), &size, |b, _| {
                b.iter(|| {
                    black_box(
                        noodles_reader.query(&region(CHROM, START, end)).unwrap().sequence().len(),
                    )
                });
            });
        }

        group.finish();
    }
}

// ---------------------------------------------------------------------------
// Group 4: name lookup cost on a many-contig reference
// ---------------------------------------------------------------------------

/// 3000 contigs, querying the last one. `noodles_fasta::fai::Index::query` is a
/// linear `iter().find()` with string comparison; seqair uses an `FxHashMap`.
/// A 1 kb fetch keeps I/O small so the lookup dominates.
fn fasta_name_lookup(c: &mut Criterion) {
    use noodles::fasta as nfasta;
    use rust_htslib::faidx;

    let path = data::many_contig_fasta();
    let last = "contig_2999";

    let mut group = c.benchmark_group("fasta_name_lookup/last_of_3000");
    group.throughput(Throughput::Bytes(1_000));

    let mut seqair_reader = seqair::fasta::IndexedFastaReader::open(&path).unwrap();
    let htslib_reader = faidx::Reader::from_path(&path).unwrap();
    let mut noodles_reader =
        nfasta::io::indexed_reader::Builder::default().build_from_path(&path).unwrap();
    let mut reuse = Vec::new();

    group.bench_function("seqair_fetch_into", |b| {
        b.iter(|| {
            seqair_reader.fetch_seq_into(last, pos(0), pos(1_000), &mut reuse).unwrap();
            black_box(reuse.len())
        });
    });
    group.bench_function("htslib", |b| {
        b.iter(|| black_box(htslib_reader.fetch_seq(last, 0, 999).unwrap().len()));
    });
    group.bench_function("noodles", |b| {
        b.iter(|| {
            black_box(noodles_reader.query(&region(last, 0, 1_000)).unwrap().sequence().len())
        });
    });

    group.finish();
}

criterion_group!(benches, fasta_open_and_fetch, fasta_open, fasta_fetch_steady, fasta_name_lookup);
criterion_main!(benches);
