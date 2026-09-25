//! Criterion benchmarks for the fqzcomp quality decoder (CRAM method 7) on the
//! hts-specs vectors, against htscodecs' `fqz_decompress` (linked through
//! `hts-sys`, the `rust-htslib` dev-dependency).
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "benches")]

// Link hts-sys's static htslib, which carries htscodecs.
extern crate rust_htslib;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use std::ffi::{c_char, c_int};
use std::hint::black_box;
use std::path::PathBuf;

unsafe extern "C" {
    fn fqz_decompress(
        input: *mut c_char,
        in_size: usize,
        out_size: *mut usize,
        lengths: *mut c_int,
        nlengths: c_int,
    ) -> *mut c_char;
}

/// Decode with htscodecs; returns the output length.
fn htscodecs_decode(src: &[u8]) -> usize {
    let mut input = src.to_vec();
    let mut out_size = 0usize;
    // SAFETY: `input` is a valid buffer of `input.len()` bytes that htscodecs
    // only reads; no record lengths are requested.
    let out = unsafe {
        fqz_decompress(
            input.as_mut_ptr().cast::<c_char>(),
            input.len(),
            &raw mut out_size,
            std::ptr::null_mut(),
            0,
        )
    };
    assert!(!out.is_null(), "htscodecs rejects a vector");
    // SAFETY: `out` is htscodecs' malloc'd output, freed once.
    unsafe { libc::free(out.cast()) };
    out_size
}

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
        let out_len = seqair::cram::fqzcomp::decode(&src).unwrap().len();
        assert_eq!(htscodecs_decode(&src), out_len);

        let mut group = c.benchmark_group(format!("fqzcomp_decode/{name}"));
        group.throughput(Throughput::Bytes(out_len as u64));
        group.sample_size(20);
        group.bench_function("seqair", |b| {
            b.iter(|| black_box(seqair::cram::fqzcomp::decode(black_box(&src)).unwrap()));
        });
        group.bench_function("htscodecs", |b| b.iter(|| black_box(htscodecs_decode(&src))));
        group.finish();
    }
}

criterion_group!(benches, fqzcomp_vectors);
criterion_main!(benches);
