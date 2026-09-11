//! The C4 cost gate: one 150 bp read against a ~200 bp haplotype through a
//! 48-column band, which the design budgets at ~20 us with 8-wide f32 SIMD and
//! ~60 us scalar.
//!
//! Do not read numbers off a laptop under load; this is meant for the Linux
//! box. The TAPS arm is there because the emission is a scalar call per cell
//! and is the part SIMD cannot help with, so the gap between the two arms is
//! the real question, not either one alone.

use compair::{
    Band, Base, BaseQuality, Betas, ConversionModel, Haplotype, Probability, Read,
    StandardEmission, Strand, TapsEmission, align_banded, align_banded_simd, align_full,
};
use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

const HAPLOTYPE_LEN: usize = 200;
const READ_LEN: usize = 150;
const READ_OFFSET: usize = 25;

/// A deterministic sequence with a realistic `CpG` density, so the TAPS arm
/// exercises the conversion rows rather than falling through everywhere.
fn sequence(len: usize) -> Vec<Base> {
    let mut state = 0x9e37_79b9_7f4a_7c15_u64;
    (0..len)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            match state % 4 {
                0 => Base::A,
                1 => Base::C,
                2 => Base::G,
                _ => Base::T,
            }
        })
        .collect()
}

fn fixture() -> (Haplotype, Read, Vec<Probability>, Band) {
    let bases = sequence(HAPLOTYPE_LEN);
    let read_bases: Vec<Base> =
        bases.get(READ_OFFSET..READ_OFFSET + READ_LEN).map(<[Base]>::to_vec).unwrap_or_default();
    let quals: Vec<BaseQuality> = (0..read_bases.len())
        .map(|index| BaseQuality::from_byte(25 + u8::try_from(index % 15).unwrap_or(0)))
        .collect();
    let read = Read::uniform(
        read_bases,
        &quals,
        BaseQuality::from_byte(45),
        BaseQuality::from_byte(45),
        BaseQuality::from_byte(10),
        Strand::OT,
    )
    .expect("the bench fixture must build");
    let betas = (0..HAPLOTYPE_LEN)
        .map(|index| {
            Probability::new(f64::from(u32::try_from(index % 11).unwrap_or(0)) / 10.0)
                .unwrap_or(Probability::ZERO)
        })
        .collect();
    let band = Band::anchored(i32::try_from(READ_OFFSET).unwrap_or(0));
    (Haplotype::new(bases), read, betas, band)
}

fn align(c: &mut Criterion) {
    let (haplotype, read, betas, band) = fixture();
    let taps = TapsEmission::new(ConversionModel::taps_default(), Betas::PerSite(&betas));

    let mut group = c.benchmark_group("align/150x48");
    group.bench_function("scalar/standard", |b| {
        b.iter(|| {
            align_banded(
                black_box(&haplotype),
                black_box(&read),
                &StandardEmission::default(),
                band,
            )
        });
    });
    group.bench_function("simd/standard", |b| {
        b.iter(|| {
            align_banded_simd(
                black_box(&haplotype),
                black_box(&read),
                &StandardEmission::default(),
                band,
            )
        });
    });
    group.bench_function("scalar/taps", |b| {
        b.iter(|| align_banded(black_box(&haplotype), black_box(&read), &taps, band));
    });
    group.bench_function("simd/taps", |b| {
        b.iter(|| align_banded_simd(black_box(&haplotype), black_box(&read), &taps, band));
    });
    group.bench_function("reference/standard", |b| {
        b.iter(|| {
            align_full(black_box(&haplotype), black_box(&read), &StandardEmission::default())
        });
    });
    group.finish();
}

criterion_group!(benches, align);
criterion_main!(benches);
