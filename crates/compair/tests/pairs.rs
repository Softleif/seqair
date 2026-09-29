//! The pairs kernel's gate: every lane of a group is bit-identical to the
//! scalar strip kernel on that lane's own pair, whatever the other lanes hold.
//!
//! The lanes of a group differ in everything a pair has -- read, read length,
//! strand, haplotype, haplotype length and band offset -- so these tests build
//! groups that differ in all of them at once, and hold every lane to
//! `align_strips` alone. Parity with the strip kernel is what lets
//! `align_reads` route a short group through the strip kernel without changing
//! an answer, so the dispatch is tested at every group size too.

#[path = "../src/pinned.rs"]
mod pinned;
#[allow(dead_code, reason = "each integration test uses a different part of this")]
mod support;

use compair::{
    BATCH, Band, Base, BaseQuality, Betas, ConversionModel, Emission, Haplotype, Log10Likelihood,
    PAIRS, PAIRS_BREAK_EVEN, Pair, Probability, Read, StandardEmission, Strand, TapsEmission,
    Workspace, align_pairs, align_reads, align_strips, align_strips_simd,
};
use hegel::TestCase;
use hegel::generators as gs;
use support::{Case, any_conversion, any_probability, any_quality_case, length};

fn bits(scores: &[Log10Likelihood]) -> Vec<u64> {
    scores.iter().map(|score| score.get().to_bits()).collect()
}

/// Each pair alone through the scalar strip kernel: the oracle.
fn one_at_a_time<E: Emission>(pairs: &[Pair<'_>], emission: &E) -> Vec<u64> {
    pairs
        .iter()
        .map(|pair| align_strips(pair.haplotype, pair.read, emission, pair.band).get().to_bits())
        .collect()
}

/// The pairs kernel through every instance there is: the eight-lane one at
/// the best level and at every level the CPU has, the one-lane one, and the
/// free function -- all against the strip kernel, through one workspace so
/// nothing a call leaves behind can leak into the next.
fn check<E: Emission>(pairs: &[Pair<'_>], emission: &E, name: &str) {
    let want = one_at_a_time(pairs, emission);
    let mut workspace = Workspace::new();
    let mut out = Vec::new();
    workspace.align_pairs(pairs, emission, &mut out);
    assert_eq!(&want, &bits(&out), "{}: strips vs eight-lane pairs", name);
    workspace.align_pairs_scalar(pairs, emission, &mut out);
    assert_eq!(&want, &bits(&out), "{}: strips vs one-lane pairs", name);
    for (level_name, level) in support::levels() {
        workspace.align_pairs_at(level, pairs, emission, &mut out);
        assert_eq!(&want, &bits(&out), "{}: strips vs pairs at {}", name, level_name);
    }
    assert_eq!(&want, &bits(&align_pairs(pairs, emission)), "{}: free function", name);
}

/// A band width for a pair: mostly the default, sometimes narrower, so that
/// the groups `align_pairs` cuts at a change of width are exercised too.
#[hegel::composite]
fn any_width(tc: &TestCase) -> u32 {
    if tc.draw_silent(gs::weighted_booleans(0.25)) {
        tc.draw_silent(gs::sampled_from(&[16u32, 8]))
    } else {
        Band::DEFAULT_WIDTH
    }
}

// A case is up to nineteen pairs through every instance at every level
// and three emissions, so it is some forty alignments of each kernel; 128
// of them is the other parity blocks' budget in wall-clock.

/// Unrelated pairs in one call: every lane its own read, haplotype,
/// strand and offset, including bands that miss the alignment entirely,
/// and a count that leaves groups of every fill.
#[hegel::test(test_cases = pinned::cases(128))]
fn unrelated_pairs_are_bit_identical_to_one_alignment_at_a_time(tc: TestCase) {
    let count = tc.draw(length(1, 2 * PAIRS + 3));
    let cases: Vec<Case> = (0..count).map(|_| tc.draw(any_quality_case())).collect();
    let widths: Vec<u32> = tc.draw(gs::vecs(any_width()).min_size(count).max_size(count));
    let one_width = tc.draw(gs::booleans());
    let conversion = tc.draw(any_conversion());
    let uniform = tc.draw(any_probability());
    let bands: Vec<Band> = cases
        .iter()
        .zip(&widths)
        .map(|(case, &width)| {
            let width = if one_width { Band::DEFAULT_WIDTH } else { width };
            Band::new(width, case.offset).expect("the widths above are valid")
        })
        .collect();
    let pairs: Vec<Pair<'_>> = cases
        .iter()
        .zip(&bands)
        .map(|(case, &band)| Pair::new(&case.haplotype, &case.read, band))
        .collect();
    let betas = &cases.first().expect("at least one case").betas;
    check(&pairs, &StandardEmission::default(), "standard");
    check(&pairs, &TapsEmission::new(conversion, Betas::PerSite(betas)), "taps");
    check(&pairs, &TapsEmission::new(conversion, Betas::Uniform(uniform)), "uniform");
}

/// One read against eight haplotypes is the batch kernel's shape, and the
/// two kernels have to agree on it: same bits, lane for lane.
#[hegel::test(test_cases = pinned::cases(128))]
fn a_batch_shaped_group_scores_as_the_batch_kernel_does(tc: TestCase) {
    let case = tc.draw(any_quality_case());
    let conversion = tc.draw(any_conversion());
    let spread = tc.draw(gs::integers::<usize>().max_value(3));
    let band = case.band();
    let taps = TapsEmission::new(conversion, Betas::PerSite(&case.betas));
    let bases = case.haplotype.bases().to_vec();
    let batch: Vec<Haplotype> = (0..BATCH)
        .map(|k| Haplotype::new(bases[..bases.len().saturating_sub(k * spread)].to_vec()))
        .collect();
    let refs: Vec<&Haplotype> = batch.iter().collect();
    let pairs: Vec<Pair<'_>> =
        batch.iter().map(|haplotype| Pair::new(haplotype, &case.read, band)).collect();
    let mut workspace = Workspace::new();
    let mut batched = Vec::new();
    workspace.align_batch(&refs, &case.read, &taps, band, &mut batched);
    let mut paired = Vec::new();
    workspace.align_pairs(&pairs, &taps, &mut paired);
    assert_eq!(bits(&batched), bits(&paired));
}

/// A deterministic pseudo-random sequence.
fn sequence(len: usize, seed: u64) -> Vec<Base> {
    let mut state = seed;
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

/// A shadow-scoring locus, the shape `align_reads` exists for: a reference
/// window, an insertion and a deletion allele, and reads of 100 to 150 bases
/// from both strands, each placed at its own offset -- several of them a
/// multiple of eight long, the length where the last row is also a
/// renormalisation row.
struct Locus {
    haplotypes: Vec<Haplotype>,
    reads: Vec<(Read, Band)>,
    betas: Vec<Probability>,
}

fn locus(reads: usize) -> Locus {
    let reference = sequence(250, 0x9e37_79b9_7f4a_7c15);
    let mut insertion = reference.clone();
    insertion.splice(125..125, [Base::A, Base::C]);
    let mut deletion = reference.clone();
    deletion.drain(125..128);
    let haplotypes: Vec<Haplotype> =
        [reference.clone(), insertion, deletion].into_iter().map(Haplotype::new).collect();

    let reads = (0..reads)
        .map(|index| {
            let len = [150usize, 144, 128, 101, 136, 150, 120, 149][index % 8];
            let start = (index * 7) % (250 - len);
            let mut bases: Vec<Base> = reference[start..start + len].to_vec();
            if let Some(base) = bases.get_mut((index * 13) % len) {
                *base = base.inverse();
            }
            let quals: Vec<BaseQuality> = (0..len)
                .map(|at| BaseQuality::from_byte(20 + u8::try_from((at + index) % 21).unwrap_or(0)))
                .collect();
            let gaps: Vec<BaseQuality> = (0..len)
                .map(|at| BaseQuality::from_byte(if at % 17 == 0 { 30 } else { 45 }))
                .collect();
            let strand = if index % 3 == 0 { Strand::OB } else { Strand::OT };
            let read = Read::new(
                bases,
                &quals,
                &gaps,
                &gaps,
                &[BaseQuality::from_byte(10); 150][..len],
                strand,
            )
            .expect("the locus's reads build");
            // The seeded offset, give or take a few bases of anchor error.
            let jitter = i32::try_from(index % 11).unwrap_or(0) - 5;
            let offset = i32::try_from(start).unwrap_or(0) + jitter;
            (read, Band::new(48, offset).expect("48 is a valid width"))
        })
        .collect();
    let betas = (0..260)
        .map(|at| {
            Probability::new(f64::from(u32::try_from(at % 11).unwrap_or(0)) / 10.0)
                .unwrap_or(Probability::ZERO)
        })
        .collect();
    Locus { haplotypes, reads, betas }
}

/// Full-length reads through nineteen renormalisations, mixed lengths and
/// strands, every lane bit-identical to the strip kernel. The property tests'
/// reads are under 70 bases; this is the size the shadow path scores.
#[test]
fn a_shadow_locus_scores_as_the_strip_kernel_does() {
    let Locus { haplotypes, reads, betas } = locus(40);
    let taps = TapsEmission::new(ConversionModel::taps_default(), Betas::PerSite(&betas));
    let pairs: Vec<Pair<'_>> = reads
        .iter()
        .flat_map(|(read, band)| haplotypes.iter().map(move |h| Pair::new(h, read, *band)))
        .collect();
    let want = one_at_a_time(&pairs, &taps);
    assert!(want.iter().all(|&score| f64::from_bits(score).is_finite()), "every pair aligns");
    check(&pairs, &taps, "taps");
    check(&pairs, &StandardEmission::default(), "standard");
}

/// `align_reads` packs the reads-by-haplotypes product eight at a time and
/// sends a short last group through the strip kernel. Whichever it picks, the
/// scores have to be the strip kernel's, in read-major order -- at every
/// product size up to three full groups, and through a workspace that has
/// just scored something else.
#[test]
fn the_dispatch_never_changes_a_score() {
    let Locus { haplotypes, reads, betas } = locus(3 * PAIRS);
    let taps = TapsEmission::new(ConversionModel::taps_default(), Betas::PerSite(&betas));
    let mut workspace = Workspace::new();
    let mut out = Vec::new();
    for count in 1..=haplotypes.len() {
        let refs: Vec<&Haplotype> = haplotypes.iter().take(count).collect();
        for read_count in 0..=3 * PAIRS {
            let reads: Vec<(&Read, Band)> =
                reads.iter().take(read_count).map(|(read, band)| (read, *band)).collect();
            let want: Vec<u64> = reads
                .iter()
                .flat_map(|&(read, band)| {
                    refs.iter()
                        .map(move |h| align_strips_simd(h, read, &taps, band).get().to_bits())
                })
                .collect();
            workspace.align_reads(&refs, &reads, &taps, &mut out);
            assert_eq!(bits(&out), want, "{read_count} reads x {count} haplotypes, workspace");
            assert_eq!(
                bits(&align_reads(&refs, &reads, &taps)),
                want,
                "{read_count} reads x {count} haplotypes, free function"
            );
        }
    }
}

/// Every group size from one pair to a full group, through `align_pairs`:
/// the fill changes nothing but the cost, and the sizes on either side of
/// `PAIRS_BREAK_EVEN` go through different kernels in `align_reads`.
#[test]
fn every_fill_scores_the_same() {
    let Locus { haplotypes, reads, betas } = locus(PAIRS);
    let taps = TapsEmission::new(ConversionModel::taps_default(), Betas::PerSite(&betas));
    let all: Vec<Pair<'_>> = reads
        .iter()
        .zip(haplotypes.iter().cycle())
        .map(|((read, band), haplotype)| Pair::new(haplotype, read, *band))
        .collect();
    let want = one_at_a_time(&all, &taps);
    let mut workspace = Workspace::new();
    let mut out = Vec::new();
    for fill in 1..=PAIRS {
        let group = &all[..fill];
        workspace.align_pairs(group, &taps, &mut out);
        assert_eq!(bits(&out), want[..fill], "a group of {fill}");
    }
    const { assert!(PAIRS_BREAK_EVEN <= PAIRS) };
}

/// Both entry points replace what is in `out` rather than appending to it.
#[test]
fn a_reused_output_buffer_holds_only_this_call() {
    let Locus { haplotypes, reads, .. } = locus(5);
    let refs: Vec<&Haplotype> = haplotypes.iter().collect();
    let reads: Vec<(&Read, Band)> = reads.iter().map(|(read, band)| (read, *band)).collect();
    let pairs: Vec<Pair<'_>> =
        reads.iter().map(|&(read, band)| Pair::new(refs[0], read, band)).collect();
    let mut workspace = Workspace::new();
    let mut out = Vec::new();
    for _ in 0..3 {
        workspace.align_reads(&refs, &reads, &StandardEmission::default(), &mut out);
        assert_eq!(out.len(), refs.len() * reads.len(), "align_reads must not accumulate");
    }
    for _ in 0..3 {
        workspace.align_pairs(&pairs, &StandardEmission::default(), &mut out);
        assert_eq!(out.len(), pairs.len(), "align_pairs must not accumulate");
    }
}

/// An empty read or haplotype, a band that misses, and a lane whose offset is
/// wildly out of range: each scores impossible, and none of them moves the
/// scores of the live lanes beside it -- in particular the far offset must not
/// set the group's shifted columns.
#[test]
fn a_dead_lane_is_impossible_and_changes_nothing_beside_it() {
    let Locus { haplotypes, reads, .. } = locus(4);
    let empty_haplotype = Haplotype::new(Vec::new());
    let empty_read = Read::uniform(
        Vec::<Base>::new(),
        &[],
        BaseQuality::from_byte(45),
        BaseQuality::from_byte(45),
        BaseQuality::from_byte(10),
        Strand::OT,
    )
    .expect("an empty read builds");
    let (read, band) = (&reads[0].0, reads[0].1);
    let far = Band::new(48, i32::MAX).expect("48 is a valid width");
    let pairs = [
        Pair::new(&haplotypes[0], read, band),
        Pair::new(&empty_haplotype, read, band),
        Pair::new(&haplotypes[1], &empty_read, band),
        Pair::new(&haplotypes[2], read, far),
        Pair::new(&haplotypes[2], read, Band::new(48, i32::MIN).expect("valid width")),
        Pair::new(&haplotypes[1], &reads[1].0, reads[1].1),
        Pair::new(&haplotypes[2], &reads[2].0, reads[2].1),
        Pair::new(&haplotypes[0], &reads[3].0, reads[3].1),
    ];
    let standard = StandardEmission::default();
    let want = one_at_a_time(&pairs, &standard);
    for (lane, &score) in want.iter().enumerate().take(5).skip(1) {
        assert_eq!(f64::from_bits(score), Log10Likelihood::IMPOSSIBLE.get(), "lane {lane}");
    }
    let mut workspace = Workspace::new();
    let mut out = Vec::new();
    workspace.align_pairs(&pairs, &standard, &mut out);
    assert_eq!(bits(&out), want);
}

/// The pairs kernel derives a read's rows and a haplotype's columns once per
/// call and reuses them for every lane that holds the same read or
/// haplotype, keyed by address. So a workspace that scores one locus, drops
/// it, and scores a different one whose reads and haplotypes land at the same
/// addresses -- which the allocator makes likely here -- must derive them
/// afresh: every call is held to the strip kernel on its own pairs.
#[test]
fn a_new_call_never_reuses_what_an_earlier_call_derived() {
    let mut workspace = Workspace::new();
    let mut out = Vec::new();
    for seed in 0..12u64 {
        let reference = sequence(220, 0x9e37_79b9_7f4a_7c15 ^ (seed * 0x2545_f491));
        let haplotypes: Vec<Haplotype> =
            vec![Haplotype::new(reference.clone()), Haplotype::new(reference[3..].to_vec())];
        let strand = if seed % 2 == 0 { Strand::OT } else { Strand::OB };
        let reads: Vec<(Read, Band)> = (0..5)
            .map(|index| {
                let start = 10 + index * 9;
                let quals: Vec<BaseQuality> = (0..120)
                    .map(|at| {
                        BaseQuality::from_byte(15 + u8::try_from((at + index) % 25).unwrap_or(0))
                    })
                    .collect();
                let read = Read::uniform(
                    reference[start..start + 120].to_vec(),
                    &quals,
                    BaseQuality::from_byte(40 + u8::try_from(seed).unwrap_or(0)),
                    BaseQuality::from_byte(45),
                    BaseQuality::from_byte(10),
                    strand,
                )
                .expect("the reads build");
                let offset = i32::try_from(start).unwrap_or(0);
                (read, Band::new(48, offset).expect("48 is a valid width"))
            })
            .collect();
        let betas: Vec<Probability> = (0..220)
            .map(|at| {
                Probability::new(
                    f64::from(u32::try_from((at + 3 * seed as usize) % 11).unwrap_or(0)) / 10.0,
                )
                .unwrap_or(Probability::ZERO)
            })
            .collect();
        let taps = TapsEmission::new(ConversionModel::taps_default(), Betas::PerSite(&betas));
        let refs: Vec<&Haplotype> = haplotypes.iter().collect();
        let reads: Vec<(&Read, Band)> = reads.iter().map(|(read, band)| (read, *band)).collect();
        let want: Vec<u64> = reads
            .iter()
            .flat_map(|&(read, band)| {
                refs.iter().map(move |h| align_strips(h, read, &taps, band).get().to_bits())
            })
            .collect();
        workspace.align_reads(&refs, &reads, &taps, &mut out);
        assert_eq!(bits(&out), want, "locus {seed}");
    }
}

/// A group's lanes read their rows and columns from the call's caches until
/// its kernel has run, so a table one lane found in the cache must survive
/// the misses of the lanes after it. Both caches evict first-in first-out:
/// after sixteen reads (or sixty-four haplotypes) have filled one, a group
/// whose first lane hits the oldest entry and whose other seven lanes miss
/// evicts seven slots starting at that very one -- unless it is pinned.
#[test]
fn a_group_never_evicts_a_table_its_own_lanes_read() {
    let Locus { haplotypes, reads, betas } = locus(24);
    let taps = TapsEmission::new(ConversionModel::taps_default(), Betas::PerSite(&betas));
    let mut workspace = Workspace::new();
    let mut out = Vec::new();

    // Rows: reads 0-15 fill the sixteen slots, then read 0 leads reads 16-22.
    let order: Vec<usize> = (0..16).chain([0]).chain(16..23).collect();
    let reference = haplotypes.first().expect("the locus has haplotypes");
    let pairs: Vec<Pair<'_>> = order
        .iter()
        .map(|&index| {
            let (read, band) = reads.get(index).expect("24 reads");
            Pair::new(reference, read, *band)
        })
        .collect();
    workspace.align_pairs(&pairs, &taps, &mut out);
    assert_eq!(bits(&out), one_at_a_time(&pairs, &taps), "rows");

    // Columns: sixty-four haplotypes fill the slots for one strand, then the
    // first leads seven new ones.
    let shifted: Vec<Haplotype> =
        (0..71).map(|seed| Haplotype::new(sequence(220, 0x5851_f42d_4c95_7f2d ^ seed))).collect();
    let (read, band) = reads.first().expect("24 reads");
    let order: Vec<usize> = (0..64).chain([0]).chain(64..71).collect();
    let pairs: Vec<Pair<'_>> = order
        .iter()
        .map(|&index| Pair::new(shifted.get(index).expect("71 haplotypes"), read, *band))
        .collect();
    workspace.align_pairs(&pairs, &taps, &mut out);
    assert_eq!(bits(&out), one_at_a_time(&pairs, &taps), "columns");
}
