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
use seqair::bam::BaseModState;
use seqair::bam::seq::{
    decode_bases_into, decode_iupac_bases_into, decode_seq, decode_seq_scalar,
    encode_iupac_bases_into,
};
use seqair::cram::rans_nx16;
use seqair_types::{Base, IupacBase};
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

/// Read lengths for the `IupacBase` kernels: a short read, where dispatch and
/// the tail dominate, and long reads.
const READ_LENGTHS: [usize; 3] = [150, 10_000, 100_000];

/// Each `IupacBase` bulk conversion against the element-wise loop it
/// replaces (`scalar`) and the tightest table loop on plain bytes (`table`),
/// in one run so all three see the same machine state.
fn iupac(c: &mut Criterion) {
    let mut group = c.benchmark_group("iupac");
    for len in READ_LENGTHS {
        group.throughput(Throughput::Bytes(len as u64));
        let ascii: Vec<u8> = bytes(len, 5).iter().map(|b| b"ACGTNacgtnRY="[*b as usize % 13]).collect();
        let mut bases = IupacBase::from_ascii_slice(&ascii);
        let packed = bytes(len.div_ceil(2), 6);
        let mut codes_out: Vec<IupacBase> = Vec::with_capacity(len);
        let mut bytes_out: Vec<u8> = Vec::with_capacity(len);
        let mut codes: Vec<u8> = bases.iter().map(|b| b.code()).collect();
        let complement: [u8; 16] = std::array::from_fn(|i| IupacBase::ALL[i].complement().code());
        let pairs: [[u8; 2]; 256] = std::array::from_fn(|b| [(b >> 4) as u8, (b & 15) as u8]);
        let mut bench = |name: &str, f: &mut dyn FnMut()| {
            group.bench_function(BenchmarkId::new(name, len), |b| b.iter(&mut *f));
        };

        bench("from_ascii/scalar", &mut || {
            codes_out.clear();
            codes_out.extend(black_box(&ascii).iter().map(|&b| IupacBase::from_ascii(b)));
        });
        bench("from_ascii/table", &mut || {
            bytes_out.clear();
            bytes_out.extend(
                black_box(&ascii).iter().map(|&b| IupacBase::CODE_OF_ASCII[usize::from(b)]),
            );
        });
        bench("from_ascii/simd", &mut || {
            codes_out.clear();
            IupacBase::from_ascii_into(black_box(&ascii), &mut codes_out);
        });
        bench("to_ascii/scalar", &mut || {
            bytes_out.clear();
            bytes_out.extend(black_box(&bases).iter().map(|b| b.to_ascii()));
        });
        bench("to_ascii/table", &mut || {
            bytes_out.clear();
            bytes_out.extend(black_box(&codes).iter().map(|&c| IupacBase::ASCII[usize::from(c)]));
        });
        bench("to_ascii/simd", &mut || {
            bytes_out.clear();
            IupacBase::to_ascii_into(black_box(&bases), &mut bytes_out);
        });
        bench("complement/scalar", &mut || {
            for b in black_box(&mut bases).iter_mut() {
                *b = b.complement();
            }
        });
        bench("complement/table", &mut || {
            for c in black_box(&mut codes).iter_mut() {
                *c = complement[usize::from(*c)];
            }
        });
        bench("complement/simd", &mut || IupacBase::complement_in_place(black_box(&mut bases)));
        bench("reverse_complement/scalar", &mut || {
            let bases = black_box(&mut bases);
            bases.reverse();
            for b in bases.iter_mut() {
                *b = b.complement();
            }
        });
        bench("reverse_complement/table", &mut || {
            let codes = black_box(&mut codes);
            codes.reverse();
            for c in codes.iter_mut() {
                *c = complement[usize::from(*c)];
            }
        });
        bench("reverse_complement/simd", &mut || {
            IupacBase::reverse_complement_in_place(black_box(&mut bases));
        });
        // What `RawBamRecord::iupac_bases` does a base at a time.
        bench("unpack/scalar", &mut || {
            let packed = black_box(&packed);
            codes_out.clear();
            codes_out.extend((0..len).map(|i| {
                let byte = packed[i / 2];
                IupacBase::from_nibble(if i % 2 == 0 { byte >> 4 } else { byte })
            }));
        });
        bench("unpack/table", &mut || {
            bytes_out.clear();
            bytes_out.extend(
                black_box(&packed[..len / 2]).iter().flat_map(|&b| pairs[usize::from(b)]),
            );
        });
        bench("unpack/simd", &mut || {
            codes_out.clear();
            codes_out.resize(len, IupacBase::N);
            decode_iupac_bases_into(black_box(&packed), &mut codes_out);
        });
        bench("unpack_ascii/scalar", &mut || {
            black_box(decode_seq_scalar(black_box(&packed), len));
        });
        bench("unpack_ascii/simd", &mut || {
            black_box(decode_seq(black_box(&packed), len));
        });
        // The pair loop `RawBamRecord::set_seq` used.
        bench("pack/scalar", &mut || {
            bytes_out.clear();
            bytes_out.extend(black_box(&bases).chunks(2).map(|pair| {
                let hi = pair.first().map_or(0, |b| b.code());
                let lo = pair.get(1).map_or(0, |b| b.code());
                (hi << 4) | lo
            }));
        });
        bench("pack/table", &mut || {
            bytes_out.clear();
            let (whole, _) = black_box(&codes).as_chunks::<2>();
            bytes_out.extend(whole.iter().map(|&[hi, lo]| (hi << 4) | lo));
        });
        bench("pack/simd", &mut || {
            bytes_out.clear();
            encode_iupac_bases_into(black_box(&bases), &mut bytes_out);
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

/// A random read with an ONT-style `C+h?;C+m?;` pair over the Cs a caller
/// picked with probability `1 / every`: `every = 4` is about what a CpG-only
/// caller does on random sequence, `every = 20` on CpG-depleted human DNA.
fn modbam_record(len: usize, every: u8, seed: u64) -> (Vec<Base>, Vec<u8>, Vec<u8>) {
    let seq: Vec<Base> = bytes(len, seed)
        .iter()
        .map(|b| [Base::A, Base::C, Base::G, Base::T][*b as usize % 4])
        .collect();
    let pick = bytes(len, seed ^ 0x55);
    let mut deltas = Vec::new();
    let mut skipped = 0u32;
    for (b, p) in seq.iter().zip(&pick) {
        if *b != Base::C {
            continue;
        }
        if p % every == 0 {
            deltas.push(skipped);
            skipped = 0;
        } else {
            skipped += 1;
        }
    }
    let list = deltas.iter().fold(String::new(), |mut list, d| {
        use std::fmt::Write as _;
        write!(list, ",{d}").unwrap();
        list
    });
    let mm = format!("C+h?{list};C+m?{list};").into_bytes();
    let ml = bytes(deltas.len() * 2, seed ^ 0xAA);
    (seq, mm, ml)
}

fn base_mod_parse(c: &mut Criterion) {
    let mut group = c.benchmark_group("base_mod_parse");
    for len in [150usize, 10_000, 100_000] {
        for every in [4u8, 20] {
            let (seq, mm, ml) = modbam_record(len, every, 7);
            group.throughput(Throughput::Elements(len as u64));
            let id = format!("{len}/1in{every}");
            group.bench_with_input(BenchmarkId::new("forward", &id), &(), |b, _| {
                b.iter(|| {
                    BaseModState::parse(black_box(&mm), black_box(&ml), black_box(&seq), false)
                        .unwrap()
                });
            });
            // The same deltas resolved against the reverse complement: the
            // walk runs from the end of the stored read towards its start.
            let rc: Vec<Base> = seq.iter().rev().map(Base::inverse).collect();
            group.bench_with_input(BenchmarkId::new("reverse", &id), &(), |b, _| {
                b.iter(|| {
                    BaseModState::parse(black_box(&mm), black_box(&ml), black_box(&rc), true)
                        .unwrap()
                });
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bam_seq_decode, ascii_to_base, iupac, rans_nx16_order0, base_mod_parse);
criterion_main!(benches);
