//! Front-to-back reading: `BgzfStream` and `BamStreamReader` at 0, 1 and 4
//! inflation threads, against noodles' single- and multi-threaded readers.
//!
//! The input is `tests/data/test.bam`'s records sixteen times over (~25 MB),
//! built in memory once; it is also written to `CARGO_TARGET_TMPDIR` so the
//! numbers can be put next to `samtools view -c -@N` on the same file.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "benches")]
#![allow(clippy::arithmetic_side_effects, reason = "benches")]
#![allow(clippy::print_stdout, reason = "benches report where the input file is")]

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use seqair::bam::stream::BamStreamReader;
use seqair::bam::writer::BamWriterBuilder;
use seqair::bam::{BgzfStream, RawBamRecord};
use std::hint::black_box;
use std::io::{Cursor, Read};
use std::num::NonZero;
use std::sync::LazyLock;

const BAM_PATH: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.bam");
const COPIES: usize = 16;

/// The compressed input, and its uncompressed size.
static INPUT: LazyLock<(Vec<u8>, u64)> = LazyLock::new(|| {
    let mut reader = BamStreamReader::new(std::fs::File::open(BAM_PATH).unwrap()).unwrap();
    let records: Vec<RawBamRecord> = reader.records().collect::<Result<_, _>>().unwrap();
    let header = reader.into_header();
    let mut writer = BamWriterBuilder::to_writer(Vec::new(), &header).build().unwrap();
    for _ in 0..COPIES {
        for record in &records {
            writer.write_raw(record).unwrap();
        }
    }
    let (bytes, _) = writer.finish().unwrap();
    let path = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("stream_bench.bam");
    std::fs::write(&path, &bytes).unwrap();
    println!("stream bench input: {} ({} bytes)", path.display(), bytes.len());
    let mut plain = Vec::new();
    noodles_bgzf::io::Reader::new(bytes.as_slice()).read_to_end(&mut plain).unwrap();
    (bytes, plain.len() as u64)
});

fn bgzf_stream(c: &mut Criterion) {
    let (bytes, plain_len) = &*INPUT;
    let mut group = c.benchmark_group("bgzf_stream");
    group.throughput(Throughput::Bytes(*plain_len));
    group.sample_size(20);

    for threads in [0usize, 1, 4] {
        group.bench_function(BenchmarkId::new("seqair", threads), |b| {
            b.iter(|| {
                let mut stream =
                    BgzfStream::builder(Cursor::new(bytes)).threads(threads).build().unwrap();
                let mut total = 0usize;
                loop {
                    let n = stream.fill().unwrap().len();
                    if n == 0 {
                        break;
                    }
                    total += n;
                    stream.advance(n);
                }
                black_box(total)
            });
        });
    }
    group.bench_function(BenchmarkId::new("noodles", 0), |b| {
        b.iter(|| {
            let mut reader = noodles_bgzf::io::Reader::new(Cursor::new(bytes));
            black_box(std::io::copy(&mut reader, &mut std::io::sink()).unwrap())
        });
    });
    group.bench_function(BenchmarkId::new("noodles_mt", 4), |b| {
        b.iter(|| {
            let mut reader = noodles_bgzf::io::MultithreadedReader::with_worker_count(
                NonZero::new(4).unwrap(),
                Cursor::new(bytes.clone()),
            );
            black_box(std::io::copy(&mut reader, &mut std::io::sink()).unwrap())
        });
    });
    group.finish();
}

fn bam_stream(c: &mut Criterion) {
    let (bytes, plain_len) = &*INPUT;
    let mut group = c.benchmark_group("bam_stream");
    group.throughput(Throughput::Bytes(*plain_len));
    group.sample_size(20);

    for threads in [0usize, 1, 4] {
        group.bench_function(BenchmarkId::new("seqair", threads), |b| {
            b.iter(|| {
                let stream =
                    BgzfStream::builder(Cursor::new(bytes)).threads(threads).build().unwrap();
                let mut reader = BamStreamReader::from_bgzf(stream).unwrap();
                let mut record = RawBamRecord::default();
                let mut n = 0usize;
                while reader.read_record(&mut record).unwrap() {
                    n += 1;
                }
                black_box(n)
            });
        });
    }
    group.bench_function(BenchmarkId::new("noodles", 0), |b| {
        b.iter(|| {
            let mut reader = noodles::bam::io::Reader::new(Cursor::new(bytes));
            reader.read_header().unwrap();
            let mut record = noodles::bam::Record::default();
            let mut n = 0usize;
            while reader.read_record(&mut record).unwrap() != 0 {
                n += 1;
            }
            black_box(n)
        });
    });
    group.bench_function(BenchmarkId::new("noodles_mt", 4), |b| {
        b.iter(|| {
            let inner = noodles_bgzf::io::MultithreadedReader::with_worker_count(
                NonZero::new(4).unwrap(),
                Cursor::new(bytes.clone()),
            );
            let mut reader = noodles::bam::io::Reader::from(inner);
            reader.read_header().unwrap();
            let mut record = noodles::bam::Record::default();
            let mut n = 0usize;
            while reader.read_record(&mut record).unwrap() != 0 {
                n += 1;
            }
            black_box(n)
        });
    });
    group.finish();
}

criterion_group!(benches, bgzf_stream, bam_stream);
criterion_main!(benches);
