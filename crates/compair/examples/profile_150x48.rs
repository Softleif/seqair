//! The bench fixture in a loop, so `samply` has something to sample.
//!
//! ```bash
//! cargo build --profile profiling -p pair-hmm --example profile_150x48
//! samply record --save-only --unstable-presymbolicate -o p.json --rate 4000 -- \
//!   ./target/profiling/examples/profile_150x48 simd-taps 200000
//! ```
#![allow(clippy::print_stdout, reason = "this example exists to print a checksum")]

use std::hint::black_box;

use compair::{
    Band, Base, BaseQuality, Betas, ConversionModel, Haplotype, Probability, Read,
    StandardEmission, Strand, TapsEmission, align_banded, align_banded_simd,
};

const HAPLOTYPE_LEN: usize = 200;
const READ_LEN: usize = 150;
const READ_OFFSET: usize = 25;

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

fn main() {
    let mut args = std::env::args().skip(1);
    let which = args.next().unwrap_or_else(|| "simd-taps".to_owned());
    let rounds: usize = args.next().and_then(|n| n.parse().ok()).unwrap_or(200_000);

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
    .expect("the fixture must build");
    let betas: Vec<Probability> = (0..HAPLOTYPE_LEN)
        .map(|index| {
            Probability::new(f64::from(u32::try_from(index % 11).unwrap_or(0)) / 10.0)
                .unwrap_or(Probability::ZERO)
        })
        .collect();
    let band = Band::anchored(i32::try_from(READ_OFFSET).unwrap_or(0));
    let haplotype = Haplotype::new(bases);
    let taps = TapsEmission::new(ConversionModel::taps_default(), Betas::PerSite(&betas));
    let standard = StandardEmission::default();

    let mut checksum = 0.0f64;
    for _ in 0..rounds {
        let score = match which.as_str() {
            "scalar-standard" => align_banded(black_box(&haplotype), &read, &standard, band),
            "scalar-taps" => align_banded(black_box(&haplotype), &read, &taps, band),
            "simd-standard" => align_banded_simd(black_box(&haplotype), &read, &standard, band),
            _ => align_banded_simd(black_box(&haplotype), &read, &taps, band),
        };
        checksum += score.get();
    }
    println!("{which}: {rounds} rounds, checksum {checksum}");
}
