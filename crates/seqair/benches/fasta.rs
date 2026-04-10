//! Criterion benchmarks for FASTA indexed fetch.
//!
//! Compares seqair against htslib and noodles.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "benches")]
#![allow(clippy::cast_possible_truncation, reason = "benches")]

use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

const FASTA_PATH: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.fasta.gz");
const CHROM: &str = "chr19";

// ---------------------------------------------------------------------------
// Group 3: FASTA fetch (100KB region)
// ---------------------------------------------------------------------------

fn fasta_fetch(c: &mut Criterion) {
    use noodles::core::{Position, Region};
    use noodles::fasta as nfasta;
    use rust_htslib::faidx;

    let mut group = c.benchmark_group("fasta_fetch");

    const FASTA_START: u64 = 6_100_000;
    const FASTA_END: u64 = 6_200_000;

    group.bench_function("seqair", |b| {
        b.iter(|| {
            let path = std::path::Path::new(FASTA_PATH);
            let mut reader = seqair::fasta::IndexedFastaReader::open(path).unwrap();
            let seq = reader
                .fetch_seq(
                    CHROM,
                    seqair::bam::Pos::<seqair::bam::Zero>::new(FASTA_START as u32).unwrap(),
                    seqair::bam::Pos::<seqair::bam::Zero>::new(FASTA_END as u32).unwrap(),
                )
                .unwrap();
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
            let region = Region::new(
                CHROM,
                Position::try_from(FASTA_START as usize + 1).unwrap()
                    ..=Position::try_from(FASTA_END as usize).unwrap(),
            );
            let record = reader.query(&region).unwrap();
            black_box(record.sequence().len())
        });
    });

    group.finish();
}

criterion_group!(benches, fasta_fetch);
criterion_main!(benches);
