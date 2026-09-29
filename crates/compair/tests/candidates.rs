//! `Workspace::candidates` against the per-pair path it caches setup for.

#[path = "../src/pinned.rs"]
mod pinned;
#[allow(dead_code, reason = "each integration test uses a different part of this")]
mod support;

use compair::{
    Band, Base, BaseQuality, Betas, ConversionModel, Emission, Haplotype, Log10Likelihood,
    Probability, Read, StandardEmission, Strand, TapsEmission, Workspace, align_candidates,
    align_strips, align_strips_simd,
};
use hegel::TestCase;
use hegel::generators as gs;
use support::{any_conversion, any_probability, arbitrary_case, levels, mirror_strand};

fn bits(scores: &[Log10Likelihood]) -> Vec<u64> {
    scores.iter().map(|score| score.get().to_bits()).collect()
}

/// The same read on the other strand: the same bases and qualities, so only
/// the column tracks a strand selects can tell the two apart.
fn on_other_strand(read: &Read) -> Option<Read> {
    Read::new(
        read.bases().to_vec(),
        read.base_quals(),
        read.insertion_quals(),
        read.deletion_quals(),
        read.gap_quals(),
        mirror_strand(read.strand()),
    )
    .ok()
}

/// Every read against every haplotype through one `Candidates`, in the given
/// order, next to the free functions on a fresh workspace per pair.
fn check<E: Emission + Copy>(
    workspace: &mut Workspace,
    haplotypes: &[Haplotype],
    emission: E,
    reads: &[(&Read, Band)],
) {
    let mut candidates = workspace.candidates(haplotypes, emission);
    assert_eq!(candidates.len(), haplotypes.len());
    let mut out = vec![Log10Likelihood::new(1.0); 3];
    for &(read, band) in reads {
        let simd: Vec<Log10Likelihood> =
            haplotypes.iter().map(|h| align_strips_simd(h, read, &emission, band)).collect();
        let scalar: Vec<Log10Likelihood> =
            haplotypes.iter().map(|h| align_strips(h, read, &emission, band)).collect();
        candidates.align_strips_simd(read, band, &mut out);
        assert_eq!(bits(&out), bits(&simd), "eight lanes, {:?}", read.strand());
        candidates.align_strips(read, band, &mut out);
        assert_eq!(bits(&out), bits(&scalar), "one lane, {:?}", read.strand());
        for (level_name, level) in levels() {
            candidates.align_strips_simd_at(level, read, band, &mut out);
            assert_eq!(bits(&out), bits(&scalar), "{}, {:?}", level_name, read.strand());
            // The dispatching entry point, batch kernel and all: the batch
            // is bit-identical to the strip kernel, so the per-pair strip
            // scores are its oracle too.
            candidates.align_at(level, read, band, &mut out);
            assert_eq!(bits(&out), bits(&scalar), "align, {}, {:?}", level_name, read.strand());
        }
        let refs: Vec<&Haplotype> = haplotypes.iter().collect();
        let fresh = align_candidates(&refs, read, &emission, band);
        candidates.align(read, band, &mut out);
        assert_eq!(bits(&out), bits(&fresh), "align vs align_candidates");
    }
}

/// Scoring through prepared candidates is scoring one pair at a time, to
/// the bit: with both strands interleaved, so a haplotype's column tracks
/// are derived for one strand and then reused after the other's; with
/// bands that differ per read, which the whole-haplotype columns have to
/// serve alike; and with one workspace carried across emissions and
/// haplotype sets, so a second `candidates` call that inherited anything
/// from the first would show here.
#[hegel::test(test_cases = pinned::cases(64))]
fn candidates_score_what_one_pair_at_a_time_scores(tc: TestCase) {
    let case = tc.draw(arbitrary_case());
    let conversion = tc.draw(any_conversion());
    let uniform = tc.draw(any_probability());
    let count = tc.draw(gs::integers::<usize>().min_value(1).max_value(17));
    let shift = tc.draw(gs::integers::<i32>().min_value(-12).max_value(11));
    let width = tc.draw(gs::integers::<u32>().min_value(2).max_value(79));
    // `count` variants of the haplotype, each trimmed, edited or both
    // by its index, so a batch's lanes differ in length and in sequence;
    // up to seventeen, so the dispatch sees full batches, a short
    // remainder and a single straggler.
    let bases = case.haplotype.bases().to_vec();
    let haplotypes: Vec<Haplotype> = (0..count)
        .map(|k| {
            let mut variant =
                bases.get(..bases.len().saturating_sub(k % 4)).unwrap_or_default().to_vec();
            let at = (k * 7) % variant.len().max(1);
            if k % 3 == 1
                && let Some(base) = variant.get_mut(at)
            {
                *base = base.inverse();
            }
            Haplotype::new(variant)
        })
        .collect();
    let fewer = [Haplotype::new(bases.iter().rev().copied().collect::<Vec<_>>())];
    let other = on_other_strand(&case.read).expect("valid by construction");
    let wide = Band::new(width, case.offset.saturating_add(shift)).expect("a legal width");
    let reads =
        [(&case.read, case.band()), (&other, case.band()), (&case.read, wide), (&other, wide)];

    let taps = TapsEmission::new(conversion, Betas::PerSite(&case.betas));
    let uniform = TapsEmission::new(conversion, Betas::Uniform(uniform));
    let standard = StandardEmission::default();
    let mut workspace = Workspace::new();
    check(&mut workspace, &haplotypes, taps, &reads);
    check(&mut workspace, &haplotypes, standard, &reads);
    check(&mut workspace, &fewer, uniform, &reads);
    check(&mut workspace, &haplotypes, uniform, &reads);
}

/// A new candidate set derives its columns afresh, even over the same
/// haplotypes: the same workspace scored under `StandardEmission` and then
/// under `TapsEmission` gives each emission's own numbers. The fixture is
/// chosen so the two emissions disagree, which is what makes a stale column
/// visible at all; that is asserted first.
#[test]
fn a_new_candidate_set_forgets_the_last_one() {
    let haplotype = Haplotype::from_ascii(b"ACGTTACGGATCCGATTACGGCATCGATCGGACGTTAGC");
    let read_bases: Vec<Base> = haplotype.bases().get(4..34).unwrap_or_default().to_vec();
    let quals = vec![BaseQuality::from_byte(30); read_bases.len()];
    // Every `C` read as `T`, so the TAPS conversion rows are what score it.
    let converted: Vec<Base> =
        read_bases.iter().map(|&base| if base == Base::C { Base::T } else { base }).collect();
    let read = Read::uniform(
        converted,
        &quals,
        BaseQuality::from_byte(45),
        BaseQuality::from_byte(45),
        BaseQuality::from_byte(10),
        Strand::OT,
    )
    .expect("valid");
    let band = Band::anchored(4);
    let haplotypes = [haplotype.clone()];
    let standard = StandardEmission::default();
    let taps = TapsEmission::new(
        ConversionModel::taps_default(),
        Betas::Uniform(Probability::new(0.9).expect("in [0, 1]")),
    );
    let want_standard = bits(&[align_strips_simd(&haplotype, &read, &standard, band)]);
    let want_taps = bits(&[align_strips_simd(&haplotype, &read, &taps, band)]);
    assert_ne!(want_standard, want_taps, "the fixture has to tell the emissions apart");

    let mut workspace = Workspace::new();
    let mut out = Vec::new();
    for _ in 0..2 {
        workspace.candidates(&haplotypes, &standard).align_strips_simd(&read, band, &mut out);
        assert_eq!(bits(&out), want_standard);
        workspace.candidates(&haplotypes, &taps).align_strips_simd(&read, band, &mut out);
        assert_eq!(bits(&out), want_taps);
    }
    // A full batch, so `align` takes the batch kernel and its interleaved
    // columns are what has to be forgotten.
    let batch = vec![haplotype; compair::BATCH];
    for _ in 0..2 {
        workspace.candidates(&batch, &standard).align(&read, band, &mut out);
        assert_eq!(bits(&out), want_standard.repeat(compair::BATCH));
        workspace.candidates(&batch, &taps).align(&read, band, &mut out);
        assert_eq!(bits(&out), want_taps.repeat(compair::BATCH));
    }
}

/// The degenerate inputs the per-pair path scores as impossible are
/// impossible here too, one slot per haplotype, and `out` holds only this
/// call's scores.
#[test]
fn empty_haplotypes_and_missed_bands_are_impossible() {
    let haplotypes = [Haplotype::new(Vec::new()), Haplotype::from_ascii(b"ACGTACGTAC")];
    let read = Read::uniform(
        Haplotype::from_ascii(b"ACGTA").bases().to_vec(),
        &[BaseQuality::from_byte(30); 5],
        BaseQuality::from_byte(45),
        BaseQuality::from_byte(45),
        BaseQuality::from_byte(10),
        Strand::OB,
    )
    .expect("valid");
    let mut workspace = Workspace::new();
    let mut candidates = workspace.candidates(&haplotypes, StandardEmission::default());
    let mut out = vec![Log10Likelihood::new(0.0); 7];
    candidates.align_strips_simd(&read, Band::anchored(0), &mut out);
    assert_eq!(out.len(), 2);
    assert_eq!(out.first(), Some(&Log10Likelihood::IMPOSSIBLE));
    assert!(out.get(1).is_some_and(|score| score.get().is_finite()));
    // A band entirely past the haplotype's end.
    let missed = Band::new(4, 400).expect("a legal band");
    candidates.align_strips_simd(&read, missed, &mut out);
    assert_eq!(out, vec![Log10Likelihood::IMPOSSIBLE; 2]);
}
