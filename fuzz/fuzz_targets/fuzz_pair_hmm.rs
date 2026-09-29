#![no_main]

//! Fuzz the `compair` pair-HMM end to end: arbitrary haplotype, read, quality
//! tracks, band and emission model through the reference, the diagonal
//! kernels and the strip kernels, with the invariants the crate documents as
//! the oracle.
//!
//! - the scalar and SIMD instances of each traversal are bit-identical, fresh
//!   or through a reused `Workspace`;
//! - every kernel's score is the `f64` recurrence over its band, at any
//!   score -- where `f32` cannot hold it, the entry point rescores in `f64`;
//! - no kernel ever returns a `NaN`, only a finite score or `IMPOSSIBLE`;
//! - a band never beats the unbanded reference, and widening a band never
//!   lowers the score;
//! - a band that contains the whole matrix is the reference to rounding.

use arbitrary::Arbitrary;
use compair::{
    BATCH, Band, Base, BaseQuality, Betas, ConversionModel, Emission, Haplotype, Log10Likelihood,
    Probability, Read, StandardEmission, Strand, TapsEmission, Workspace, align_banded,
    align_banded_f64, align_banded_simd, align_full, align_strips, align_strips_simd,
};
use libfuzzer_sys::fuzz_target;

/// Sequences a few hundred bases long are what the kernels are for; longer
/// input only slows the fuzzer down without reaching new code.
const MAX_LEN: usize = 300;

#[derive(Arbitrary, Debug)]
struct Input {
    haplotype: Vec<u8>,
    read: Vec<u8>,
    base_quals: Vec<u8>,
    insertion_quals: Vec<u8>,
    deletion_quals: Vec<u8>,
    gap_quals: Vec<u8>,
    bottom_strand: bool,
    width: u32,
    offset: i32,
    efficiency: u8,
    false_conversion: u8,
    floor: u8,
    betas: Vec<u8>,
}

fn quals(track: &[u8], len: usize) -> Vec<BaseQuality> {
    (0..len)
        .map(|index| {
            // Cycle a possibly short track over the read; Phred 0..=93 is the
            // range a FASTQ can carry, and `UNAVAILABLE` (255) is rejected by
            // `Read::new`, which is a path worth reaching too.
            let byte = track.get(index % track.len().max(1)).copied().unwrap_or(30);
            BaseQuality::from_byte(if byte == 255 { 255 } else { byte % 94 })
        })
        .collect()
}

fn probability(byte: u8) -> Probability {
    Probability::new(f64::from(byte) / 255.0).unwrap_or(Probability::ZERO)
}

fn check_parity(name: &str, scalar: Log10Likelihood, simd: Log10Likelihood) {
    assert_eq!(scalar.get().to_bits(), simd.get().to_bits(), "{name}: {scalar:?} vs {simd:?}");
    assert!(!scalar.get().is_nan(), "{name}: NaN");
}

/// Every invariant one emission model has to satisfy on one pair.
fn check<E: Emission>(
    name: &str,
    haplotype: &Haplotype,
    read: &Read,
    emission: &E,
    band: Band,
    workspace: &mut Workspace,
) {
    let full = align_full(haplotype, read, emission);
    let scalar = align_banded(haplotype, read, emission, band);
    let simd = align_banded_simd(haplotype, read, emission, band);
    let reused = workspace.align_banded_simd(haplotype, read, emission, band);
    check_parity(name, scalar, simd);
    check_parity(name, scalar, reused);
    let strips = align_strips(haplotype, read, emission, band);
    let strips_simd = align_strips_simd(haplotype, read, emission, band);
    let strips_reused = workspace.align_strips_simd(haplotype, read, emission, band);
    check_parity(name, strips, strips_simd);
    check_parity(name, strips, strips_reused);
    // Prepared candidates are the per-pair path: alone through the strip
    // kernel, and as a full batch through the batch kernel.
    let mut prepared = Vec::new();
    for group in [&[haplotype][..], &[haplotype; BATCH][..]] {
        workspace.candidates(group, emission).align(read, band, &mut prepared);
        assert_eq!(prepared.len(), group.len(), "{name}: one score per haplotype");
        for score in &prepared {
            check_parity(name, strips, *score);
        }
    }
    assert!(!full.get().is_nan(), "{name}: reference NaN");
    // Both traversals are the `f64` recurrence over the band, to rounding,
    // at any score: they agree on whether a path exists at all, and on its
    // score.
    let want = align_banded_f64(haplotype, read, emission, band).get();
    assert!(!want.is_nan(), "{name}: f64 NaN");
    for (kernel, got) in [("diagonals", scalar.get()), ("strips", strips.get())] {
        assert_eq!(got.is_finite(), want.is_finite(), "{name}: {kernel} {got} vs f64 {want}");
        if want.is_finite() {
            assert!(
                (got - want).abs() <= 1e-4 * (1.0 + want.abs()),
                "{name}: {kernel} {got} vs f64 {want}"
            );
        }
    }
    if scalar.get().is_finite() {
        assert!(
            scalar.get() <= full.get() + 1e-3 * (1.0 + full.get().abs()),
            "{name}: banded {scalar:?} beats the reference {full:?}"
        );
    }
    if strips.get().is_finite() {
        assert!(
            strips.get() <= full.get() + 1e-3 * (1.0 + full.get().abs()),
            "{name}: strips {strips:?} beats the reference {full:?}"
        );
    }
    // Widening can only add paths.
    if let Ok(wider) = Band::new(band.width().saturating_add(8), band.offset())
        && scalar.get().is_finite()
    {
        let b = align_banded_simd(haplotype, read, emission, wider).get();
        assert!(
            scalar.get() <= b + 1e-4 * (1.0 + b.abs()),
            "{name}: narrow {scalar:?} beats wider {b}"
        );
        let b = align_strips_simd(haplotype, read, emission, wider).get();
        assert!(
            strips.get() <= b + 1e-4 * (1.0 + b.abs()),
            "{name}: strips narrow {strips:?} beats wider {b}"
        );
    }
}

fuzz_target!(|input: Input| {
    let Input {
        haplotype,
        read,
        base_quals,
        insertion_quals,
        deletion_quals,
        gap_quals,
        bottom_strand,
        width,
        offset,
        efficiency,
        false_conversion,
        floor,
        betas,
    } = input;
    let haplotype =
        Haplotype::from_ascii(haplotype.get(..haplotype.len().min(MAX_LEN)).unwrap_or(&[]));
    let read_bases: Vec<Base> = read
        .get(..read.len().min(MAX_LEN))
        .unwrap_or(&[])
        .iter()
        .copied()
        .map(Base::from)
        .collect();
    let n = read_bases.len();
    let strand = if bottom_strand { Strand::OB } else { Strand::OT };
    let Ok(read) = Read::new(
        read_bases,
        &quals(&base_quals, n),
        &quals(&insertion_quals, n),
        &quals(&deletion_quals, n),
        &quals(&gap_quals, n),
        strand,
    ) else {
        return;
    };

    let Ok(band) = Band::new(width, offset) else {
        return;
    };
    let conversion = ConversionModel::new(probability(efficiency), probability(false_conversion));
    let betas: Vec<Probability> = betas.iter().copied().map(probability).collect();
    let standard = StandardEmission::default().with_artifact_floor(probability(floor));
    let taps = TapsEmission::new(conversion, Betas::PerSite(&betas))
        .with_artifact_floor(probability(floor));

    let mut workspace = Workspace::new();
    check("standard", &haplotype, &read, &standard, band, &mut workspace);
    check("taps", &haplotype, &read, &taps, band, &mut workspace);

    // A band holding the whole matrix is the reference, at any score.
    let span = 2 * (haplotype.len() + read.len() + 4);
    if let Ok(width) = u32::try_from(span)
        && let Ok(whole) = Band::new(width, 0)
    {
        let full = align_full(&haplotype, &read, &standard).get();
        let banded = align_banded_simd(&haplotype, &read, &standard, whole).get();
        let strips = align_strips_simd(&haplotype, &read, &standard, whole).get();
        for (kernel, got) in [("band", banded), ("strips", strips)] {
            assert_eq!(got.is_finite(), full.is_finite(), "whole-matrix {kernel} {got} vs {full}");
            if full.is_finite() {
                assert!(
                    (got - full).abs() < 1e-3 * (1.0 + full.abs()),
                    "whole-matrix {kernel} {got} against the reference {full}"
                );
            }
        }
    }
});
