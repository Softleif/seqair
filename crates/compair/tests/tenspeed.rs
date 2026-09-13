//! compair over the "10s" benchmark: 3,550 real read-haplotype pairs.
//!
//! This dataset has no trustworthy expected scores -- see
//! `tests/data/PROVENANCE.md` for why the `10s.out` that ships beside it is not
//! one -- so nothing here compares against a file. What it is good for is
//! *shape*: 3,550 pairs cut from real data, with reads from 10 to 247 bases,
//! haplotypes from 41 to 263, and `N` in the reads. The 104 bundled GATK
//! vectors have none of that spread, so these are the tests that exercise it:
//!
//! - the band's contract at scale. The band may lose a path and score lower; it
//!   may never score higher than the full matrix. This is the largest sample we
//!   have of that direction holding.
//! - what the band actually costs and saves on realistic inputs, reported
//!   rather than asserted, because the number is the output.
//! - on x86_64, agreement with gkl -- an independent Rust port of the kernel
//!   GATK calls through JNI, on exactly the bytes gkl's own harness would get.

use compair::{Band, StandardEmission, align_banded_simd, align_full, align_strips_simd};

#[path = "support/tenspeed.rs"]
mod tenspeed;

/// How far a banded score may sit *above* the full matrix before we call it a
/// broken contract rather than rounding. The kernels are `f32` and the reported
/// value is a log10, so a few ulp of the sum lands here.
const OVER_ESTIMATE_SLACK: f64 = 1e-3;

#[test]
fn the_dataset_is_the_one_the_papers_benchmark() {
    let pairs = tenspeed::pairs().expect("10s.in parses");
    assert_eq!(pairs.len(), 3550, "the 10s benchmark is 3,550 pairs");
    assert_eq!(
        tenspeed::matrix_cells(&pairs),
        62_380_634,
        "the full matrices are 62,380,634 cells"
    );
}

/// The band never scores above the full matrix, on every pair and both
/// traversals. This is the contract.
#[test]
fn the_band_never_over_estimates() {
    let pairs = tenspeed::pairs().expect("10s.in parses");
    let standard = StandardEmission::default();
    let mut worst_over = f64::NEG_INFINITY;
    for (index, pair) in pairs.iter().enumerate() {
        let band = Band::anchored(pair.offset);
        let full = align_full(&pair.haplotype, &pair.read, &standard).get();
        for (name, got) in [
            ("banded", align_banded_simd(&pair.haplotype, &pair.read, &standard, band).get()),
            ("strips", align_strips_simd(&pair.haplotype, &pair.read, &standard, band).get()),
        ] {
            let over = got - full;
            assert!(
                over < OVER_ESTIMATE_SLACK,
                "pair {index}: {name} scored {got}, above the full matrix's {full} by {over}"
            );
            worst_over = worst_over.max(over);
        }
    }
    println!("worst over-estimate across 3,550 pairs x 2 traversals: {worst_over:+.3e}");
}

/// What a 46-column band costs and saves on real inputs, against compair's own
/// full matrix. Reports the numbers; asserts only that the default band is not
/// wildly wrong for this shape of data.
#[test]
fn the_band_holds_the_optimal_path_on_most_pairs() {
    let pairs = tenspeed::pairs().expect("10s.in parses");
    let standard = StandardEmission::default();
    let mut inside = 0usize;
    let mut worst_loss = 0.0f64;
    let mut total_loss = 0.0f64;
    for pair in &pairs {
        let band = Band::anchored(pair.offset);
        let full = align_full(&pair.haplotype, &pair.read, &standard).get();
        let got = align_strips_simd(&pair.haplotype, &pair.read, &standard, band).get();
        let loss = full - got;
        if loss < 1e-3 {
            inside += 1;
        } else {
            worst_loss = worst_loss.max(loss);
            total_loss += loss;
        }
    }
    let matrix = tenspeed::matrix_cells(&pairs);
    let banded = tenspeed::band_cells(&pairs, |pair| Band::anchored(pair.offset));
    #[allow(clippy::cast_precision_loss, reason = "reporting a ratio of cell counts")]
    let saved = 100.0 - 100.0 * banded as f64 / matrix as f64;
    #[allow(clippy::cast_precision_loss, reason = "reporting a mean over a few thousand pairs")]
    let mean_loss =
        if inside == pairs.len() { 0.0 } else { total_loss / (pairs.len() - inside) as f64 };
    println!(
        "band {} wide: {inside}/{} pairs within 1e-3 of the full matrix; \
         worst loss {worst_loss:.3} log10, mean loss over the rest {mean_loss:.3}; \
         {banded} of {matrix} cells visited ({saved:.1}% saved)",
        Band::DEFAULT_WIDTH,
        pairs.len(),
    );
    assert!(
        inside * 100 >= pairs.len() * 70,
        "only {inside} of {} pairs stayed inside the band",
        pairs.len()
    );
}

/// An independent implementation of the same recurrence, on the same bytes.
/// gkl is a port of Intel's GKL, which is what GATK actually calls, so
/// agreement here is agreement with the production kernel rather than with our
/// own reading of the spec.
///
/// Run this twice to cover both implementations. gkl's `c-avx` feature
/// *replaces* its Rust AVX kernels with Intel's original C++, so the vector
/// arm below is whichever the build selected:
///
/// ```text
/// cargo test -p compair --test tenspeed --release -- --nocapture
/// cargo test -p compair --test tenspeed --release --features gkl-cpp -- --nocapture
/// ```
///
/// The scalar arm is gkl's Rust in both runs (`forward_f64x1` is not gated), so
/// the second run is what puts compair next to Intel's C++ directly.
#[cfg(target_arch = "x86_64")]
#[test]
fn the_reference_agrees_with_gkl() {
    let pairs = tenspeed::pairs().expect("10s.in parses");
    let vector = if cfg!(feature = "gkl-cpp") { "gkl-cpp" } else { "gkl-rust" };
    println!(
        "gkl vector kernels: {vector}. available: f32x8 {}, f64x4 {}, f32x16 {}, f64x8 {}",
        gkl::pairhmm::forward_f32x8().is_some(),
        gkl::pairhmm::forward_f64x4().is_some(),
        gkl::pairhmm::forward_f32x16().is_some(),
        gkl::pairhmm::forward_f64x8().is_some(),
    );

    let standard = StandardEmission::default();
    let ours: Vec<f64> =
        pairs.iter().map(|pair| align_full(&pair.haplotype, &pair.read, &standard).get()).collect();

    // `f64` kernels should agree with our `f64` matrix to rounding; the `f32`
    // ones carry their own error and only have to be close.
    let mut checked = 0usize;
    for (name, tolerance, scores) in [
        ("gkl-rust/f64x1", 1e-6, Some(run(&pairs, gkl::pairhmm::forward_f64x1()))),
        (
            &*format!("{vector}/f64x4"),
            1e-6,
            gkl::pairhmm::forward_f64x4().map(|forward| run(&pairs, forward)),
        ),
        (
            &*format!("{vector}/f32x8"),
            1e-3,
            gkl::pairhmm::forward_f32x8().map(|forward| run32(&pairs, forward)),
        ),
        (
            &*format!("{vector}/f32x16"),
            1e-3,
            gkl::pairhmm::forward_f32x16().map(|forward| run32(&pairs, forward)),
        ),
        (
            &*format!("{vector}/f64x8"),
            1e-6,
            gkl::pairhmm::forward_f64x8().map(|forward| run(&pairs, forward)),
        ),
    ] {
        let Some(theirs) = scores else {
            println!("{name:>18}: not available on this CPU");
            continue;
        };
        let mut worst = 0.0f64;
        let mut worst_at = 0usize;
        for (index, (ours, theirs)) in ours.iter().zip(&theirs).enumerate() {
            if (ours - theirs).abs() > worst {
                worst = (ours - theirs).abs();
                worst_at = index;
            }
        }
        println!("{name:>18}: worst disagreement over {} pairs {worst:.3e}", pairs.len());
        assert!(worst < tolerance, "{name} disagrees by {worst} at pair {worst_at}");
        checked += 1;
    }
    assert!(checked > 0, "no gkl kernel was available to check against");
}

#[cfg(target_arch = "x86_64")]
fn run(pairs: &[tenspeed::Pair], forward: gkl::pairhmm::Forward) -> Vec<f64> {
    pairs
        .iter()
        .map(|pair| {
            let raw = &pair.raw;
            forward(
                &raw.haplotype,
                &raw.read,
                &raw.base_quals,
                &raw.insertion_quals,
                &raw.deletion_quals,
                &raw.gap_quals,
            )
        })
        .collect()
}

#[cfg(target_arch = "x86_64")]
fn run32(pairs: &[tenspeed::Pair], forward: gkl::pairhmm::ForwardF32) -> Vec<f64> {
    pairs
        .iter()
        .map(|pair| {
            let raw = &pair.raw;
            f64::from(forward(
                &raw.haplotype,
                &raw.read,
                &raw.base_quals,
                &raw.insertion_quals,
                &raw.deletion_quals,
                &raw.gap_quals,
            ))
        })
        .collect()
}
