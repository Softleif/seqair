//! Criterion benchmarks for the fqzcomp quality decoder (CRAM method 7) on the
//! hts-specs vectors, against htscodecs' `fqz_decompress` (linked through
//! `hts-sys`, the `rust-htslib` dev-dependency).
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "benches")]

use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;

#[path = "support/htscodecs.rs"]
#[allow(dead_code, reason = "each bench target uses a subset")]
mod htscodecs;

fn fqzcomp_vectors(c: &mut Criterion) {
    let dir = PathBuf::from(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../tests/data/hts-specs/cram/codecs/fqzcomp"
    ));
    let mut vectors: Vec<PathBuf> =
        std::fs::read_dir(&dir).unwrap().map(|e| e.unwrap().path()).collect();
    vectors.sort();

    for path in vectors {
        let name = path.file_name().unwrap().to_str().unwrap().to_owned();
        let src = std::fs::read(&path).unwrap();
        // Both decoders must agree before either is timed (the integration
        // tests check them against the originals).
        let ours = seqair::cram::fqzcomp::decode(&src).unwrap();
        let theirs = htscodecs::fqz_decompress(&src).expect("htscodecs rejects a vector");
        assert!(ours == theirs.as_slice(), "{name}: seqair and htscodecs decode differently");

        let mut group = c.benchmark_group(format!("fqzcomp_decode/{name}"));
        group.throughput(Throughput::Bytes(ours.len() as u64));
        group.sample_size(20);
        group.bench_function("seqair", |b| {
            b.iter(|| seqair::cram::fqzcomp::decode(black_box(&src)).unwrap());
        });
        group.bench_function("htscodecs", |b| {
            b.iter(|| htscodecs::fqz_decompress(black_box(&src)).unwrap());
        });
        group.finish();
    }
}

criterion_group!(benches, fqzcomp_vectors);
criterion_main!(benches);
