//! The SIMD kernels one at a time, on synthetic input, so a change to one of
//! them can be A/B'd without the I/O around it:
//!
//! ```text
//! git checkout <before> && cargo bench --bench simd_kernels -- --save-baseline before
//! git checkout <after>  && cargo bench --bench simd_kernels -- --baseline before
//! ```
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "benches")]
#![allow(clippy::cast_possible_truncation, reason = "benches")]
#![allow(clippy::arithmetic_side_effects, clippy::indexing_slicing, reason = "benches")]

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use seqair::bam::seq::{decode_bases_into, decode_seq};
use seqair::cram::rans_nx16;
use seqair_types::Base;
use std::hint::black_box;

/// A short read, a long read, and a bulk buffer: the first is dominated by
/// dispatch and the scalar tail, the last by the vector loop.
const LENGTHS: [usize; 3] = [150, 10_000, 1 << 20];

fn bytes(len: usize, mut seed: u64) -> Vec<u8> {
    (0..len)
        .map(|_| {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed as u8
        })
        .collect()
}

fn bam_seq_decode(c: &mut Criterion) {
    let mut group = c.benchmark_group("bam_seq_decode");
    for len in LENGTHS {
        let packed = bytes(len.div_ceil(2), 1);
        let mut out = vec![0u8; len];
        group.throughput(Throughput::Bytes(len as u64));
        group.bench_with_input(BenchmarkId::new("decode_seq", len), &len, |b, &len| {
            b.iter(|| decode_seq(black_box(&packed), len));
        });
        group.bench_with_input(BenchmarkId::new("decode_bases_into", len), &len, |b, &len| {
            b.iter(|| decode_bases_into(black_box(&packed), len, black_box(&mut out)));
        });
    }
    group.finish();
}

fn ascii_to_base(c: &mut Criterion) {
    let mut group = c.benchmark_group("ascii_to_base");
    for len in LENGTHS {
        // The kernel is branch-free, so rewriting already-converted bytes
        // costs the same as the first pass: no need to reset between runs.
        let mut buf: Vec<u8> =
            bytes(len, 2).iter().map(|b| b"ACGTacgtNn"[*b as usize % 10]).collect();
        group.throughput(Throughput::Bytes(len as u64));
        group.bench_with_input(BenchmarkId::from_parameter(len), &len, |b, _| {
            b.iter(|| Base::convert_ascii_in_place(black_box(&mut buf)));
        });
    }
    group.finish();
}

fn push_uint7(stream: &mut Vec<u8>, val: u32) {
    let mut groups = vec![(val & 0x7F) as u8];
    let mut v = val >> 7;
    while v > 0 {
        groups.push((v & 0x7F) as u8 | 0x80);
        v >>= 7;
    }
    stream.extend(groups.iter().rev());
}

/// A 32-state order-0 rANS Nx16 stream over two equiprobable symbols, with
/// random renormalization bytes to spare.
fn rans_stream(len: usize) -> Vec<u8> {
    let mut stream = vec![0x04]; // FLAG_N32, order 0
    push_uint7(&mut stream, len as u32);
    stream.extend_from_slice(&[0, 1, 0, 0]); // alphabet {0, 1}
    push_uint7(&mut stream, 2048);
    push_uint7(&mut stream, 2048);
    for s in bytes(32, 3) {
        stream.extend_from_slice(&(0x0080_0000u32 | u32::from(s)).to_le_bytes());
    }
    // One renorm per state per 16 steps, two bytes each, doubled for margin.
    stream.extend_from_slice(&bytes(len / 4 + 1024, 4));
    stream
}

fn rans_nx16_order0(c: &mut Criterion) {
    let mut group = c.benchmark_group("rans_nx16_order0_32state");
    for len in LENGTHS {
        let stream = rans_stream(len);
        assert_eq!(rans_nx16::decode(&stream, 0).unwrap().len(), len);
        group.throughput(Throughput::Bytes(len as u64));
        group.bench_with_input(BenchmarkId::from_parameter(len), &len, |b, _| {
            b.iter(|| rans_nx16::decode(black_box(&stream), 0).unwrap());
        });
    }
    group.finish();
}

criterion_group!(benches, bam_seq_decode, ascii_to_base, rans_nx16_order0);
criterion_main!(benches);
