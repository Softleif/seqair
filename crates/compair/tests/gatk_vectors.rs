#[allow(dead_code, reason = "each integration test uses a different part of this")]
mod support;

use compair::{
    Band, BaseQuality, Read, StandardEmission, align_banded, align_banded_simd, align_full,
    align_strips, align_strips_simd,
};
use support::{gatk_vectors, seed_offset};

/// The gate for the reference implementation: GATK's own numbers, plus the
/// paths GATK does not allow, to 1e-5 in log10 space.
///
/// compair lets a read begin with an insertion (the crate docs, *Where
/// compair differs from GATK*), so its likelihood is GATK's plus every path
/// that inserts the read's first `k` bases and aligns the rest from a match.
/// Those paths sum, by their first state after the insertion, to
/// `match_to_insertion(1) * S_2` with `S_k = G_k + gap_continuation(k) *
/// S_(k + 1)`, `G_k` GATK's likelihood of the read from base `k` on. compair's
/// own likelihood of that suffix is `C_k = G_k + match_to_insertion(k) *
/// S_(k + 1)`, so `S_k = C_k + (gap_continuation(k) - match_to_insertion(k)) *
/// S_(k + 1)`, summed from the read's end. That checks the start's weight
/// against GATK rather than only bounding it: the residuals are the ones
/// compair had against GATK before it allowed the start (6.1e-6 at most).
/// GATK's convention alone moves every vector by 4e-5 to 1.1e-4.
#[test]
fn reference_reproduces_gatk() {
    let vectors = gatk_vectors().expect("the bundled test data parses");
    assert_eq!(vectors.len(), 104);
    let probability = |q: BaseQuality| 10f64.powf(-f64::from(q.as_byte()) / 10.0);
    let mut worst = 0.0f64;
    for (index, vector) in vectors.iter().enumerate() {
        let read = &vector.read;
        let score = |read: &Read| align_full(&vector.haplotype, read, &StandardEmission::default());
        let got = score(read).get();
        let mut suffixes = 0.0f64;
        for k in (1..read.len()).rev() {
            let suffix = Read::new(
                read.bases()[k..].to_vec(),
                &read.base_quals()[k..],
                &read.insertion_quals()[k..],
                &read.deletion_quals()[k..],
                &read.gap_quals()[k..],
                read.strand(),
            )
            .expect("a suffix of a valid read");
            let opens = probability(read.gap_quals()[k]) - probability(read.insertion_quals()[k]);
            suffixes = 10f64.powf(score(&suffix).get()) + opens * suffixes;
        }
        let begins_inserted = probability(read.insertion_quals()[0]) * suffixes;
        let want = (10f64.powf(vector.expected_log10) + begins_inserted).log10();
        let difference = (got - want).abs();
        assert!(
            difference < 1e-5,
            "vector {index}: got {got}, GATK says {} and the paths it does not allow make that \
             {want}, difference {difference}",
            vector.expected_log10
        );
        worst = worst.max(difference);
    }
    assert!(worst < 1e-5, "worst difference {worst}");
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
        let strips =
            align_strips_simd(&vector.haplotype, &vector.read, &StandardEmission::default(), band)
                .get();
        let difference = (banded - vector.expected_log10).abs();
        let strips_difference = (strips - vector.expected_log10).abs();
        if difference < 1e-3 && strips_difference < 1e-3 {
            inside += 1;
            worst = worst.max(difference).max(strips_difference);
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
        let scalar =
            align_strips(&vector.haplotype, &vector.read, &StandardEmission::default(), band);
        let simd =
            align_strips_simd(&vector.haplotype, &vector.read, &StandardEmission::default(), band);
        assert_eq!(
            scalar.get().to_bits(),
            simd.get().to_bits(),
            "vector {index}: strips scalar {scalar:?} simd {simd:?}"
        );
    }
}
