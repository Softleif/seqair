#[allow(dead_code, reason = "each integration test uses a different part of this")]
mod support;

use compair::{Band, StandardEmission, align_banded, align_banded_simd, align_full};
use support::{gatk_vectors, seed_offset};

/// The gate for the reference implementation: GATK's own numbers, to 1e-4 in
/// log10 space.
#[test]
fn reference_reproduces_gatk() {
    let vectors = gatk_vectors().expect("the bundled test data parses");
    assert_eq!(vectors.len(), 104);
    let mut worst = 0.0f64;
    for (index, vector) in vectors.iter().enumerate() {
        let got = align_full(&vector.haplotype, &vector.read, &StandardEmission::default()).get();
        let difference = (got - vector.expected_log10).abs();
        assert!(
            difference < 1e-4,
            "vector {index}: got {got}, GATK says {}, difference {difference}",
            vector.expected_log10
        );
        worst = worst.max(difference);
    }
    assert!(worst < 1e-4, "worst difference {worst}");
}

/// The same vectors through the band. They are 101 bp reads inside 164 bp
/// haplotypes with the read a near-substring, so every optimal path fits a
/// 48-column band anchored on the best-matching offset; this asserts that and
/// prints the count, because a vector whose path left the band would be scored
/// lower rather than rejected.
#[test]
fn gatk_vectors_through_the_band() {
    let vectors = gatk_vectors().expect("the bundled test data parses");
    let mut inside = 0usize;
    let mut worst = 0.0f64;
    for vector in &vectors {
        let band = Band::anchored(seed_offset(&vector.haplotype, &vector.read));
        let banded =
            align_banded(&vector.haplotype, &vector.read, &StandardEmission::default(), band).get();
        let difference = (banded - vector.expected_log10).abs();
        if difference < 1e-3 {
            inside += 1;
            worst = worst.max(difference);
        }
    }
    assert_eq!(
        inside,
        vectors.len(),
        "only {inside} of {} vectors stay inside a {}-column band",
        vectors.len(),
        Band::DEFAULT_WIDTH
    );
    assert!(worst < 1e-3, "worst banded difference {worst}");
}

#[test]
fn simd_matches_scalar_on_every_gatk_vector() {
    let vectors = gatk_vectors().expect("the bundled test data parses");
    for (index, vector) in vectors.iter().enumerate() {
        let band = Band::anchored(seed_offset(&vector.haplotype, &vector.read));
        let scalar =
            align_banded(&vector.haplotype, &vector.read, &StandardEmission::default(), band);
        let simd =
            align_banded_simd(&vector.haplotype, &vector.read, &StandardEmission::default(), band);
        assert_eq!(
            scalar.get().to_bits(),
            simd.get().to_bits(),
            "vector {index}: scalar {scalar:?} simd {simd:?}"
        );
    }
}
