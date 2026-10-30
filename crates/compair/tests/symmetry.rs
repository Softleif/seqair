//! The model scores a read and a haplotype as it scores their reverse
//! complements, wherever the read's gap qualities are the same at every
//! base: the property GATK's conventions lack and compair restores (the crate
//! docs, *Where compair differs from GATK*). With the tandem-repeat gap model
//! the deletions keep it to a bound and the insertions do not; both are
//! pinned at the end.
//!
//! The emission's half of it, `C` on OT against `G` on OB, is exact and
//! tested in `tests/emission.rs`; this is the dynamic program's half, which
//! rests on the two ends of a read being scored alike. GATK's start would fail
//! it by up to `log10(1 + match_to_insertion / (eps / 3))` on any read whose
//! last base disagrees with the haplotype.

#[allow(dead_code, reason = "each integration test uses a different part of this")]
mod support;

use compair::{
    Band, Base, BaseQuality, Betas, Haplotype, PAIRS, PcrIndelModel, Read, StandardEmission,
    Strand, TapsEmission, Workspace, align_banded_f64_rows, align_full, align_strips_simd,
};
use hegel::TestCase;
use hegel::generators::{self as gs, Generator, PrintableGenerator};
use support::{Case, any_base, any_conversion, any_probability, any_strand, length, mirror};

/// A read cut out of a haplotype, up to its ends, with its first and last
/// bases substituted at random, and one gap quality set for every base: the
/// reads whose ends GATK scores differently.
#[hegel::composite]
fn end_case_inner(tc: &TestCase) -> Case {
    let len = tc.draw_silent(length(2, 80));
    let hap: Vec<Base> = (0..len).map(|_| tc.draw_silent(any_base())).collect();
    let start = tc.draw(gs::integers::<usize>().max_value(len - 1));
    let end = tc.draw(gs::integers::<usize>().min_value(start + 1).max_value(len));
    let mut bases = hap.get(start..end).unwrap_or_default().to_vec();
    for at in [0, bases.len() - 1] {
        if tc.draw(gs::booleans())
            && let Some(base) = bases.get_mut(at)
        {
            *base = tc.draw(any_base());
        }
    }
    let quals: Vec<BaseQuality> = (0..bases.len())
        .map(|_| BaseQuality::from_byte(tc.draw_silent(gs::integers::<u8>().max_value(60))))
        .collect();
    let [insertion, deletion, gap] =
        [0; 3].map(|_| BaseQuality::from_byte(tc.draw(gs::integers::<u8>().max_value(60))));
    let strand = tc.draw_silent(any_strand());
    let read = Read::uniform(bases, &quals, insertion, deletion, gap, strand)
        .expect("tracks as long as the read");
    let betas = (0..len).map(|_| tc.draw_silent(any_probability())).collect();
    let offset = i32::try_from(start).expect("a short haplotype");
    Case { haplotype: Haplotype::new(hap), read, offset, betas }
}

fn end_case() -> impl PrintableGenerator<Case> {
    end_case_inner().print_as_debug()
}

/// Scores that agree to `f64` rounding, impossible ones included.
fn assert_same(forward: f64, backward: f64, what: &str) {
    let same = forward == backward || (forward - backward).abs() <= 1e-9 * forward.abs().max(1.0);
    assert!(same, "{what}: {forward} forward, {backward} mirrored");
}

#[hegel::test]
fn the_reference_is_mirror_symmetric(tc: TestCase) {
    let case = tc.draw(end_case());
    let mirrored = mirror(&case).expect("the mirror of a valid read is valid");
    let standard = StandardEmission::default();
    assert_same(
        align_full(&case.haplotype, &case.read, &standard).get(),
        align_full(&mirrored.haplotype, &mirrored.read, &standard).get(),
        "plain emission",
    );
    let conversion = tc.draw(any_conversion());
    let forward = TapsEmission::new(conversion, Betas::PerSite(&case.betas));
    let backward = TapsEmission::new(conversion, Betas::PerSite(&mirrored.betas));
    assert_same(
        align_full(&case.haplotype, &case.read, &forward).get(),
        align_full(&mirrored.haplotype, &mirrored.read, &backward).get(),
        "conversion-aware emission",
    );
}

/// The band is mirrored with the pair: row `i` at column `j` is row `r + 1 -
/// i` at column `h + 1 - j`, so a band around `j = i + offset` becomes one
/// around `j = i + h - r - offset`, of the same width.
#[hegel::test]
fn the_banded_reference_is_mirror_symmetric(tc: TestCase) {
    let case = tc.draw(end_case());
    let mirrored = mirror(&case).expect("the mirror of a valid read is valid");
    let (h, r) = (case.haplotype.len(), case.read.len());
    let width = tc.draw(gs::integers::<u32>().min_value(2).max_value(64));
    let offset = case.offset + tc.draw(gs::integers::<i32>().min_value(-3).max_value(3));
    let mirrored_offset =
        i32::try_from(h).expect("short") - i32::try_from(r).expect("short") - offset;
    let (Ok(band), Ok(mirrored_band)) =
        (Band::new(width, offset), Band::new(width, mirrored_offset))
    else {
        tc.reject()
    };
    let standard = StandardEmission::default();
    assert_same(
        align_banded_f64_rows(&case.haplotype, &case.read, &standard, band).get(),
        align_banded_f64_rows(&mirrored.haplotype, &mirrored.read, &standard, mirrored_band).get(),
        "banded",
    );
}

/// A band around the case's offset, and its mirror.
fn mirrored_bands(tc: &TestCase, case: &Case) -> (Band, Band) {
    let (h, r) = (case.haplotype.len(), case.read.len());
    let width = tc.draw(gs::integers::<u32>().min_value(2).max_value(64));
    let offset = case.offset + tc.draw(gs::integers::<i32>().min_value(-3).max_value(3));
    let mirrored_offset =
        i32::try_from(h).expect("short") - i32::try_from(r).expect("short") - offset;
    match (Band::new(width, offset), Band::new(width, mirrored_offset)) {
        (Ok(band), Ok(mirrored_band)) => (band, mirrored_band),
        _ => tc.reject(),
    }
}

/// The kernels a caller runs, the strip kernel and the pairs kernel behind
/// `align_reads`, to their precision: a millionth of a log10 above their trust
/// floor, which is what every kernel promises of the `f64` recurrence.
#[hegel::test]
fn the_kernels_are_mirror_symmetric(tc: TestCase) {
    let case = tc.draw(end_case());
    let mirrored = mirror(&case).expect("the mirror of a valid read is valid");
    let (band, mirrored_band) = mirrored_bands(&tc, &case);
    let standard = StandardEmission::default();
    let close = |forward: f64, backward: f64, what: &str| {
        let same = forward == backward || (forward - backward).abs() <= 2e-6;
        assert!(same, "{what}: {forward} forward, {backward} mirrored");
    };
    close(
        align_strips_simd(&case.haplotype, &case.read, &standard, band).get(),
        align_strips_simd(&mirrored.haplotype, &mirrored.read, &standard, mirrored_band).get(),
        "strips",
    );
    let mut workspace = Workspace::new();
    let (mut forward, mut backward) = (Vec::new(), Vec::new());
    // A full group of pairs, so that the pairs kernel scores them, not the
    // strip kernel it hands a short group to.
    workspace.align_reads(
        &[&case.haplotype; PAIRS],
        &[(&case.read, band)],
        &standard,
        &mut forward,
    );
    workspace.align_reads(
        &[&mirrored.haplotype; PAIRS],
        &[(&mirrored.read, mirrored_band)],
        &standard,
        &mut backward,
    );
    for (forward, backward) in forward.iter().zip(&backward) {
        close(forward.get(), backward.get(), "align_reads");
    }
}

/// A haplotype of homopolymer runs, a read cut out of it with one base
/// substituted, and the run edges where a gap quality from the repeat model
/// changes along the read.
#[hegel::composite]
fn repeat_case_inner(tc: &TestCase) -> (Haplotype, Vec<Base>) {
    let len = tc.draw_silent(length(20, 60));
    let mut hap = Vec::with_capacity(len + 8);
    while hap.len() < len {
        let base = tc.draw(any_base());
        let run = tc.draw(gs::integers::<usize>().min_value(1).max_value(8));
        hap.extend(std::iter::repeat_n(base, run));
    }
    let start = tc.draw(gs::integers::<usize>().max_value(5));
    let end = hap.len() - tc.draw(gs::integers::<usize>().max_value(5));
    let mut bases = hap.get(start..end).unwrap_or_default().to_vec();
    let at = tc.draw(gs::integers::<usize>().max_value(bases.len() - 1));
    if let Some(base) = bases.get_mut(at) {
        *base = tc.draw(any_base());
    }
    (Haplotype::new(hap), bases)
}

fn repeat_case() -> impl PrintableGenerator<(Haplotype, Vec<Base>)> {
    repeat_case_inner().print_as_debug()
}

/// rastair's tandem-repeat gap model: GATK's, with the gap open floored at Q20.
fn pcr_model() -> PcrIndelModel {
    PcrIndelModel { floor: BaseQuality::from_byte(20), ..PcrIndelModel::default() }
}

fn reverse_complement(bases: &[Base]) -> Vec<Base> {
    bases.iter().rev().map(Base::inverse).collect()
}

/// GATK's `match_to_match` of each read base, `1 - (p_ins + p_del)` with both
/// raised to Q6, in natural log.
fn ln_match_to_match(read: &Read) -> Vec<f64> {
    let open = |q: &BaseQuality| 10f64.powf(-f64::from(q.as_byte().max(6)) / 10.0);
    read.insertion_quals()
        .iter()
        .zip(read.deletion_quals())
        .map(|(insertion, deletion)| (1.0 - open(insertion) - open(deletion)).ln())
        .collect()
}

/// With the deletion qualities of the repeat model and the insertion and
/// gap-continuation qualities the same at every base, a read scores as its
/// reverse complement does but for GATK's `match_to_match`, which is all the
/// bound below allows.
///
/// A deletion between two bases opens at the higher of their qualities
/// (`Read::with_pcr_indel_model`), so it costs the same in both orientations,
/// as every other transition and emission does. `match_to_match` does not:
/// the step from read base `i - 1` to `i` takes `1 - p_ins - p_del(i)`, the
/// deletion quality of the gap after base `i`, so in the mirror the same step
/// takes the gap before base `i - 1`. So no path's weight moves by more than
/// a factor of `e^B`, with `B` the sum over the steps of how far the two
/// `ln(match_to_match)` lie apart, and neither does the total. With GATK's
/// rule, a deletion charged to the base before it, the scores part by up to
/// two in log10 here, where a deletion beside a run competes with a mismatch.
#[hegel::test]
fn the_repeat_models_deletions_are_mirror_symmetric(tc: TestCase) {
    let (haplotype, bases) = tc.draw(repeat_case());
    let r = bases.len();
    let read = |bases: Vec<Base>| {
        let quals = vec![BaseQuality::from_byte(40); r];
        let gap = BaseQuality::from_byte(10);
        let deletions =
            Read::with_pcr_indel_model(bases.clone(), &quals, &pcr_model(), gap, Strand::OT)
                .expect("a valid read")
                .deletion_quals()
                .to_vec();
        let insertions = vec![BaseQuality::from_byte(44); r];
        Read::new(bases, &quals, &insertions, &deletions, &vec![gap; r], Strand::OT)
            .expect("a valid read")
    };
    let (forward, backward) = (read(bases.clone()), read(reverse_complement(&bases)));
    let (ln_forward, ln_backward) = (ln_match_to_match(&forward), ln_match_to_match(&backward));
    // The step into forward row `i` (0-based base `i`, `1..r`) is the step
    // into mirror base `r - i` (0-based).
    let bound: f64 = (1..r)
        .filter_map(|i| Some((ln_forward.get(i)? - ln_backward.get(r - i)?).abs()))
        .sum::<f64>()
        / core::f64::consts::LN_10;
    let standard = StandardEmission::default();
    let (one, other) = (
        align_full(&haplotype, &forward, &standard).get(),
        align_full(&haplotype.reverse_complement(), &backward, &standard).get(),
    );
    assert!(
        (one - other).abs() <= bound + 1e-9,
        "{one} forward, {other} mirrored, {} apart, more than the bound {bound}",
        (one - other).abs()
    );
}

/// A property the model does **not** have, pinned so that nobody adds it as
/// a test and watches it fail: with the repeat model's insertion qualities a
/// read does not score as its reverse complement.
///
/// A multi-base insertion is charged its first base's insertion quality. At
/// a read's end, the insertion that explains a stray base can open inside a
/// homopolymer, at the repeat's Q20, and run to the read's last base; at its
/// start, the mirror of that insertion has to open at the read's first base,
/// outside the run, at Q44. Here the read's second base is a mismatch next to
/// a run of `G`, and the two orientations part by more than one in log10. A
/// symmetric charge would have to look at both ends of an insertion, which
/// GATK's recurrence, charging each transition to one row, cannot.
#[test]
fn a_multi_base_insertion_is_charged_at_its_first_base() {
    let haplotype = Haplotype::from_ascii(b"TTTTGGGGGGGGGGGGTTTTT");
    let bases = Haplotype::from_ascii(b"GCGGGGGGGGG").bases().to_vec();
    let score = |haplotype: &Haplotype, bases: Vec<Base>| {
        let quals = vec![BaseQuality::from_byte(40); bases.len()];
        let gap = BaseQuality::from_byte(10);
        let read = Read::with_pcr_indel_model(bases, &quals, &pcr_model(), gap, Strand::OT)
            .expect("a valid read");
        align_full(haplotype, &read, &StandardEmission::default()).get()
    };
    let forward = score(&haplotype, bases.clone());
    let backward = score(&haplotype.reverse_complement(), reverse_complement(&bases));
    assert!(
        (forward - backward).abs() > 1.0,
        "the asymmetry is real and worth documenting: {forward} forward, {backward} mirrored"
    );
}
