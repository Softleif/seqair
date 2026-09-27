#![no_main]

//! Fuzz the `compair` pair-HMM end to end: arbitrary haplotype, read, quality
//! tracks, band and emission model through the reference, the diagonal
//! kernels and the strip kernels, with the invariants the crate documents as
//! the oracle.
//!
//! - the scalar and SIMD instances of each traversal are bit-identical, fresh
//!   or through a reused `Workspace`;
//! - the two traversals agree on whether a path exists, and on the score
//!   where `f32` holds it;
//! - no kernel ever returns a `NaN`, only a finite score or `IMPOSSIBLE`;
//! - a band never beats the unbanded reference, and widening a band never
//!   lowers the score;
//! - a band that contains the whole matrix is the reference to rounding.

use arbitrary::Arbitrary;
use compair::{
    BATCH, Band, Base, BaseQuality, Betas, ConversionModel, Emission, Haplotype, Log10Likelihood,
    Probability, Read, StandardEmission, Strand, TapsEmission, Workspace, align_banded,
    align_banded_simd, align_full, align_strips, align_strips_simd,
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
    // The strip kernel is the same recurrence over the same band, so the two
    // traversals agree on whether a path exists at all, and on the score
    // where `f32` holds it: they renormalise at different points and sum the
    // last row in different precisions, which is rounding.
    assert_eq!(
        scalar.get().is_finite(),
        strips.get().is_finite(),
        "{name}: diagonals {scalar:?} vs strips {strips:?}"
    );
    if scalar.get().is_finite() && scalar.get() > -30.0 {
        assert!(
            (scalar.get() - strips.get()).abs() <= 1e-4 * (1.0 + scalar.get().abs()),
            "{name}: diagonals {scalar:?} vs strips {strips:?}"
        );
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
    // Widening can only add paths -- where f32 holds the score. A cell the
    // kernel flushes is below `2^-126` of its diagonal's maximum, so below
    // ~1e-38 in absolute terms, and cannot move a total of 1e-30 or more;
    // below that, a wider band's diagonals span more dynamic range and can
    // flush mass a narrower band kept, which the crate documents.
    if let Ok(wider) = Band::new(band.width().saturating_add(8), band.offset())
        && scalar.get().is_finite()
    {
        let b = align_banded_simd(haplotype, read, emission, wider);
        if b.get() > -30.0 {
            assert!(
                scalar.get() <= b.get() + 1e-4 * (1.0 + b.get().abs()),
                "{name}: narrow {scalar:?} beats wider {b:?}"
            );
        }
        let b = align_strips_simd(haplotype, read, emission, wider);
        if b.get() > -30.0 {
            assert!(
                strips.get() <= b.get() + 1e-4 * (1.0 + b.get().abs()),
                "{name}: strips narrow {strips:?} beats wider {b:?}"
            );
        }
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

    // A band holding the whole matrix is the reference, where f32 can hold
    // the score at all. A diagonal of such a band crosses every read row, so
    // its cells span the alignment's whole dynamic range and everything more
    // than ~1e-38 below the diagonal's maximum flushes to zero -- down to
    // `IMPOSSIBLE` for a read that scores hundreds of log10 below its start
    // prior. The crate documents that; the comparison is only meaningful
    // where the score is within f32's reach.
    let span = 2 * (haplotype.len() + read.len() + 4);
    if let Ok(width) = u32::try_from(span)
        && let Ok(whole) = Band::new(width, 0)
    {
        let full = align_full(&haplotype, &read, &standard).get();
        let banded = align_banded_simd(&haplotype, &read, &standard, whole).get();
        if full > -30.0 {
            assert!(
                (banded - full).abs() < 1e-3 * (1.0 + full.abs()),
                "whole-matrix band {banded} against the reference {full}"
            );
        } else if banded.is_finite() {
            assert!(
                banded <= full + 1e-3 * (1.0 + full.abs()),
                "whole-matrix band {banded} over {full}"
            );
        }
        // The strip kernel's dynamic range is eight rows whatever the band,
        // so it is held to the reference wherever the reference is finite
        // and within f32's absolute range at all.
        let strips = align_strips_simd(&haplotype, &read, &standard, whole).get();
        if full > -30.0 {
            assert!(
                (strips - full).abs() < 1e-3 * (1.0 + full.abs()),
                "whole-matrix strips {strips} against the reference {full}"
            );
        } else if strips.is_finite() {
            assert!(
                strips <= full + 1e-3 * (1.0 + full.abs()),
                "whole-matrix strips {strips} over {full}"
            );
        }
    }
});
