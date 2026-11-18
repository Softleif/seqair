//! Criterion benchmarks for the CRAM writer's rANS Nx16 encoder: the scalar
//! loops (the reference), the SIMD kernels, and htscodecs'
//! `rans_compress_to_4x16` (linked through `hts-sys`), each without
//! transforms, at order 0 and 1 with 4 and 32 states.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "benches")]
#![allow(clippy::cast_possible_truncation, clippy::indexing_slicing, reason = "benches")]

use core::range::RangeInclusive;
use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use seqair::cram::rans_nx16::decode;
use seqair::cram::writer::{rans_nx16_encode, rans_nx16_encode_methods};
use seqair_types::Pos0;
use std::hint::black_box;

#[path = "support/htscodecs.rs"]
#[allow(dead_code, reason = "each bench target uses a subset")]
mod htscodecs;

const BAM_PATH: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.bam");

/// The qualities and the bases (as ASCII) of test.bam's reads, in order:
/// what a slice's `QS` and `BA` blocks hold.
fn bam_blocks() -> (Vec<u8>, Vec<u8>) {
    let path = std::path::Path::new(BAM_PATH);
    let mut reader = seqair::bam::IndexedBamReader::open(path).unwrap();
    let tid = reader.header().tid("chr19").unwrap();
    let span = RangeInclusive { start: Pos0::new(0).unwrap(), last: Pos0::new(6_140_000).unwrap() };
    let mut store = seqair::bam::RecordStore::new();
    reader.fetch_into(tid, span, &mut store).unwrap();
    let (mut quals, mut bases) = (Vec::new(), Vec::new());
    for idx in store.indices() {
        let record = store.record(idx).unwrap();
        quals.extend(record.qual().iter().map(|q| q.get().unwrap_or(0xFF)));
        bases.extend(record.seq().iter().map(|b| b.to_ascii()));
    }
    (quals, bases)
}

/// Uniformly random bytes (xorshift64*).
fn random(len: usize) -> Vec<u8> {
    let mut x = 0x9E37_79B9_7F4A_7C15u64;
    (0..len)
        .map(|_| {
            x ^= x >> 12;
            x ^= x << 25;
            x ^= x >> 27;
            (x.wrapping_mul(0x2545_F491_4F6C_DD1D) >> 56) as u8
        })
        .collect()
}

/// Method bytes: order 0 and 1, 4 and 32 states (`X32` = 0x04).
const METHODS: [(u8, &str); 4] = [(0, "o0x4"), (4, "o0x32"), (1, "o1x4"), (5, "o1x32")];

fn rans_nx16_encode_bench(c: &mut Criterion) {
    let (quals, bases) = bam_blocks();
    let inputs: Vec<(String, Vec<u8>)> = vec![
        ("qual_4k".into(), quals[..4096].to_vec()),
        ("qual_64k".into(), quals[..65_536].to_vec()),
        ("qual_1m".into(), quals[..quals.len().min(1 << 20)].to_vec()),
        ("bases_64k".into(), bases[..65_536].to_vec()),
        ("random_64k".into(), random(65_536)),
    ];
    for (name, data) in &inputs {
        let mut group = c.benchmark_group(format!("rans_nx16_encode/{name}"));
        group.throughput(Throughput::Bytes(data.len() as u64));
        for (method, label) in METHODS {
            for scalar in [true, false] {
                let coded = rans_nx16_encode(data, method, scalar).unwrap();
                assert!(decode(&coded, data.len()).unwrap() == *data, "{name} {label}");
            }
            let theirs = htscodecs::rans_compress_4x16(data, method).unwrap();
            assert!(decode(theirs.as_slice(), data.len()).unwrap() == *data, "{name} {label}");

            group.bench_function(BenchmarkId::new("scalar", label), |b| {
                b.iter(|| rans_nx16_encode(black_box(data), method, true).unwrap());
            });
            group.bench_function(BenchmarkId::new("simd", label), |b| {
                b.iter(|| rans_nx16_encode(black_box(data), method, false).unwrap());
            });
            group.bench_function(BenchmarkId::new("htscodecs", label), |b| {
                b.iter(|| htscodecs::rans_compress_4x16(black_box(data), method).unwrap());
            });
        }
        // What the writer does with a block: every combination of htslib's
        // normal profile (order 0, 1, RLE, STRIPE, PACK, RLE + PACK; 32
        // states from 50 000 bytes on, but for STRIPE).
        let x32 = if data.len() >= 50_000 { 4 } else { 0 };
        let profile =
            [0, 1, 0x40, 0x09, 0x80, 0xC1].map(|m| if m & 0x08 == 0 { m | x32 } else { m });
        for (label, scalar) in [("scalar", true), ("simd", false)] {
            group.bench_function(BenchmarkId::new(label, "profile"), |b| {
                b.iter(|| {
                    profile.map(|m| rans_nx16_encode(black_box(data), m, scalar).unwrap().len())
                });
            });
        }
        // The same, with the work the combinations share done once, as the
        // writer's block does.
        group.bench_function(BenchmarkId::new("shared", "profile"), |b| {
            b.iter(|| rans_nx16_encode_methods(black_box(data), &profile).unwrap());
        });
        group.bench_function(BenchmarkId::new("htscodecs", "profile"), |b| {
            b.iter(|| {
                profile.map(|m| {
                    htscodecs::rans_compress_4x16(black_box(data), m).unwrap().as_slice().len()
                })
            });
        });
        group.finish();
    }
}

criterion_group!(benches, rans_nx16_encode_bench);
criterion_main!(benches);
