//! The model scores a read and a haplotype as it scores their reverse
//! complements, wherever the read's gap qualities are the same at every
//! base: the property GATK's conventions lack and compair restores (the crate
//! docs, *Where compair differs from GATK*).
//!
//! The emission's half of it, `C` on OT against `G` on OB, is exact and
//! tested in `tests/emission.rs`; this is the dynamic program's half, which
//! rests on the two ends of a read being scored alike. GATK's start would fail
//! it by up to `log10(1 + match_to_insertion / (eps / 3))` on any read whose
//! last base disagrees with the haplotype.

#[allow(dead_code, reason = "each integration test uses a different part of this")]
mod support;

use compair::{
    Band, Base, BaseQuality, Betas, Haplotype, PAIRS, Read, StandardEmission, TapsEmission,
    Workspace, align_banded_f64_rows, align_full, align_strips_simd,
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
