//! The GPU kernel's oracles.
//!
//! The kernel is a separate implementation of the batch kernel's arithmetic,
//! and the GPU compiler is free to do things to `f32` arithmetic a CPU
//! kernel's compiler is not, so its gates come in two layers:
//!
//! - **The algorithm**, with no GPU at all: the Rust transcription of the
//!   shader (`GpuPairs::emulate`) is bit-identical to the strip kernel. This is
//!   what says the band-relative row buffer, the masks and the renormalisation
//!   are right, and it runs wherever the tests do.
//! - **The device**: what a GPU returns is, pair for pair, the transcription
//!   with subnormal intermediates either kept or flushed -- the one thing the
//!   measured GPUs do differently from a CPU once contraction is blocked --
//!   and within `1e-4` in log10 of the strip kernel, like the GATK vectors.
//!   These skip, loudly, where wgpu finds no adapter.

#![cfg(feature = "gpu")]

#[allow(dead_code, reason = "each integration test uses a different part of this")]
mod support;

#[path = "support/tenspeed.rs"]
mod tenspeed;

use std::sync::{Arc, OnceLock};

use compair::{
    Band, Base, BaseQuality, Betas, Emission, Haplotype, Log10Likelihood, Probability, Read,
    StandardEmission, Strand, TapsEmission, Workspace, align_full, align_strips,
    gpu::{GpuAligner, GpuContext, GpuError, GpuPairs, Subnormals},
};
use proptest::prelude::*;
use support::{Case, any_conversion, any_probability, arbitrary_case, derived_case};

/// `|Δ log10|` the GPU may differ from the strip kernel by: the GATK vectors'
/// tolerance. Measured: 6e-7 at worst on the 10s pairs.
const TOLERANCE: f64 = 1e-4;

/// The score below which [`TOLERANCE`] is not promised; see
/// `the_gpu_is_the_strip_kernel_within_tolerance`.
const FLOOR: f64 = -60.0;

/// The shared device, or `None` where there is none; each GPU test then
/// reports that it skipped rather than passing silently.
fn context() -> Option<Arc<GpuContext>> {
    static CONTEXT: OnceLock<Option<Arc<GpuContext>>> = OnceLock::new();
    CONTEXT
        .get_or_init(|| match GpuContext::new() {
            Ok(context) => Some(context),
            Err(error) => {
                eprintln!("no GPU for the compair::gpu tests: {error}");
                None
            }
        })
        .clone()
}

/// Every case against its read under its own band width, one pair each.
fn plan<E: Emission>(cases: &[(Case, Band)], emission: &E) -> Result<GpuPairs, GpuError> {
    let mut pairs = GpuPairs::new();
    for (case, band) in cases {
        let read = pairs.push_read(&case.read, emission)?;
        let haplotype = pairs.push_haplotype(&case.haplotype, case.read.strand(), emission)?;
        pairs.push_pair(read, haplotype, *band)?;
    }
    Ok(pairs)
}

fn bits(scores: &[Log10Likelihood]) -> Vec<u64> {
    scores.iter().map(|score| score.get().to_bits()).collect()
}

/// A case with a band of its own width, from 2 up to past the unrolled
/// kernel's limit, so one launch spans several band classes.
fn banded_case() -> impl Strategy<Value = (Case, Band)> {
    (arbitrary_case(), 2u32..=160).prop_filter_map("a legal width", |(case, width)| {
        let band = Band::new(width, case.offset).ok()?;
        Some((case, band))
    })
}

/// The emissions the parity gates sweep, for one batch of cases.
fn for_each_emission(
    cases: &[(Case, Band)],
    conversion: compair::ConversionModel,
    uniform: Probability,
    mut check: impl FnMut(&str, &GpuPairs, Vec<Log10Likelihood>) -> Result<(), TestCaseError>,
) -> Result<(), TestCaseError> {
    let mut workspace = Workspace::new();
    let fail = |error: GpuError| TestCaseError::fail(format!("{error}"));
    {
        let emission = StandardEmission::default();
        let strips = cases
            .iter()
            .map(|(case, band)| {
                workspace.align_strips(&case.haplotype, &case.read, &emission, *band)
            })
            .collect();
        check("standard", &plan(cases, &emission).map_err(fail)?, strips)?;
    }
    {
        // Betas per site differ per case, so each case gets its own emission;
        // the plan is built pair by pair for the same reason.
        let mut pairs = GpuPairs::new();
        let mut strips = Vec::new();
        for (case, band) in cases {
            let emission = TapsEmission::new(conversion, Betas::PerSite(&case.betas));
            let read = pairs.push_read(&case.read, &emission).map_err(fail)?;
            let haplotype = pairs
                .push_haplotype(&case.haplotype, case.read.strand(), &emission)
                .map_err(fail)?;
            pairs.push_pair(read, haplotype, *band).map_err(fail)?;
            strips.push(workspace.align_strips(&case.haplotype, &case.read, &emission, *band));
        }
        check("taps", &pairs, strips)?;
    }
    {
        let emission = TapsEmission::new(conversion, Betas::Uniform(uniform));
        let strips = cases
            .iter()
            .map(|(case, band)| {
                workspace.align_strips(&case.haplotype, &case.read, &emission, *band)
            })
            .collect();
        check("uniform", &plan(cases, &emission).map_err(fail)?, strips)?;
    }
    Ok(())
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 512, ..ProptestConfig::default() })]

    /// The algorithm gate, CPU only: the shader's transcription with
    /// subnormals kept is the strip kernel, bit for bit, at every width and
    /// offset -- including bands that miss the alignment or the haplotype,
    /// `N` on both sides, and every emission.
    #[test]
    fn the_kernel_algorithm_is_the_strip_kernel(
        cases in prop::collection::vec(banded_case(), 1..6),
        conversion in any_conversion(),
        uniform in any_probability(),
    ) {
        for_each_emission(&cases, conversion, uniform, |name, pairs, strips| {
            prop_assert_eq!(bits(&pairs.emulate(Subnormals::Kept)), bits(&strips), "{}", name);
            Ok(())
        })?;
    }
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 64, ..ProptestConfig::default() })]

    /// The device gate: every score is the transcription's, with subnormal
    /// intermediates flushed (Metal, RADV) or kept (an IEEE-denormal GPU).
    /// Bit-exact, not a tolerance: this is what makes the kernel's numbers
    /// reproducible from one GPU to another, and it needs contraction
    /// blocked, which the default kernel does.
    #[test]
    fn the_gpu_runs_the_transcribed_kernel(
        cases in prop::collection::vec(banded_case(), 1..48),
        conversion in any_conversion(),
        uniform in any_probability(),
    ) {
        let Some(context) = context() else { return Ok(()) };
        let mut aligner = GpuAligner::new(context);
        for_each_emission(&cases, conversion, uniform, |name, pairs, _| {
            let gpu = bits(&aligner.align(pairs).map_err(|e| TestCaseError::fail(format!("{e}")))?);
            let flushed = bits(&pairs.emulate(Subnormals::Flushed));
            let kept = bits(&pairs.emulate(Subnormals::Kept));
            for (index, &score) in gpu.iter().enumerate() {
                prop_assert!(
                    Some(&score) == flushed.get(index) || Some(&score) == kept.get(index),
                    "{}: pair {} scored {} on the GPU, {:?} flushed, {:?} kept",
                    name, index, f64::from_bits(score),
                    flushed.get(index).map(|b| f64::from_bits(*b)),
                    kept.get(index).map(|b| f64::from_bits(*b)),
                );
            }
            Ok(())
        })?;
    }

    /// The tolerance oracle against the CPU's strip kernel.
    ///
    /// Within [`TOLERANCE`] above [`FLOOR`], never above it anywhere, and
    /// `IMPOSSIBLE` wherever the strip kernel is. Below the floor a GPU that
    /// flushes subnormal intermediates can lose more: a pair that bad has its
    /// surviving path so far under each row's maximum that the flushed
    /// products were the path (measured over 200,000 random pairs: nothing
    /// above -40 changes at all, 2.9e-5 at worst above -60, up to 0.76 at
    /// -85 and below). Flushing only ever removes mass, so the GPU errs the
    /// way a band does: low.
    #[test]
    fn the_gpu_is_the_strip_kernel_within_tolerance(
        cases in prop::collection::vec(banded_case(), 1..48),
        conversion in any_conversion(),
        uniform in any_probability(),
    ) {
        let Some(context) = context() else { return Ok(()) };
        let mut aligner = GpuAligner::new(context);
        for_each_emission(&cases, conversion, uniform, |name, pairs, strips| {
            let gpu = aligner.align(pairs).map_err(|e| TestCaseError::fail(format!("{e}")))?;
            for (index, (g, s)) in gpu.iter().zip(&strips).enumerate() {
                let (g, s) = (g.get(), s.get());
                if !s.is_finite() {
                    prop_assert_eq!(g, s, "{}: pair {}", name, index);
                    continue;
                }
                prop_assert!(!g.is_nan() && g <= s + TOLERANCE, "{}: pair {}: gpu {} above strips {}", name, index, g, s);
                if s > FLOOR {
                    prop_assert!((g - s).abs() <= TOLERANCE, "{}: pair {}: gpu {} strips {}", name, index, g, s);
                }
            }
            Ok(())
        })?;
    }

    /// The band contract on the GPU: a band never scores above the unbanded
    /// `f64` reference.
    #[test]
    fn the_gpu_band_never_beats_the_reference(cases in prop::collection::vec(derived_case(4), 1..32)) {
        let Some(context) = context() else { return Ok(()) };
        let emission = StandardEmission::default();
        let banded: Vec<(Case, Band)> = cases.iter().map(|case| (case.clone(), case.band())).collect();
        let pairs = plan(&banded, &emission).map_err(|e| TestCaseError::fail(format!("{e}")))?;
        let gpu = GpuAligner::new(context).align(&pairs).map_err(|e| TestCaseError::fail(format!("{e}")))?;
        for (case, score) in cases.iter().zip(&gpu) {
            let full = align_full(&case.haplotype, &case.read, &emission).get();
            prop_assert!(!score.get().is_nan());
            prop_assert!(score.get() <= full + TOLERANCE, "banded {} beats the reference {}", score.get(), full);
        }
    }
}

/// The 10s benchmark in one launch: every pair within tolerance of the strip
/// kernel. Prints how many are bit-identical, which the notes track.
#[test]
fn the_10s_pairs_are_the_strip_kernel_within_tolerance() -> Result<(), Box<dyn std::error::Error>> {
    let Some(context) = context() else { return Ok(()) };
    let groups = tenspeed::groups().map_err(|error| format!("10s: {error:?}"))?;
    let emission = StandardEmission::default();
    let mut pairs = GpuPairs::new();
    let mut strips = Vec::new();
    for group in &groups {
        let haplotypes = group
            .haplotypes
            .iter()
            .map(|h| pairs.push_haplotype(h, Strand::OT, &emission))
            .collect::<Result<Vec<_>, _>>()?;
        for (read, &offset) in group.reads.iter().zip(&group.offsets) {
            let slot = pairs.push_read(read, &emission)?;
            for (haplotype, &hap_slot) in group.haplotypes.iter().zip(&haplotypes) {
                pairs.push_pair(slot, hap_slot, Band::anchored(offset))?;
                strips.push(align_strips(haplotype, read, &emission, Band::anchored(offset)));
            }
        }
    }
    let gpu = GpuAligner::new(context).align(&pairs)?;
    assert_eq!(gpu.len(), strips.len());
    let identical = bits(&gpu).iter().zip(bits(&strips)).filter(|(a, b)| **a == *b).count();
    eprintln!("10s: {identical}/{} bit-identical to align_strips", strips.len());
    for (index, (g, s)) in gpu.iter().zip(&strips).enumerate() {
        assert!(
            (g.get() - s.get()).abs() <= TOLERANCE,
            "pair {index}: gpu {} strips {}",
            g.get(),
            s.get()
        );
    }
    Ok(())
}

fn read(bases: &[u8], strand: Strand) -> Read {
    let bases: Vec<Base> = Haplotype::from_ascii(bases).bases().to_vec();
    Read::uniform(
        bases.clone(),
        &vec![BaseQuality::from_byte(30); bases.len()],
        BaseQuality::from_byte(45),
        BaseQuality::from_byte(45),
        BaseQuality::from_byte(10),
        strand,
    )
    .expect("valid")
}

/// An aligner reused across launches of different sizes and bands scores
/// like a fresh one, and a pair scores the same alone as in a batch: nothing
/// from a previous launch or a neighbouring pair leaks in.
#[test]
fn a_reused_aligner_is_a_fresh_aligner() -> Result<(), GpuError> {
    let Some(context) = context() else { return Ok(()) };
    let long = Haplotype::from_ascii(b"ACGTTAGCATCGGATCCGATTACAGGCATTACGGATCCAGTACGGT");
    let short = Haplotype::from_ascii(b"TTAGCAT");
    let reads = [
        read(b"TAGCATCGGATCCGATTACAGG", Strand::OT),
        read(b"AGCAT", Strand::OT),
        read(b"GGATCCAGTAAAA", Strand::OT),
    ];
    let emission = StandardEmission::default();
    let mut all = GpuPairs::new();
    let mut singles = Vec::new();
    for haplotype in [&long, &short] {
        for read in &reads {
            for band in [
                Band::anchored(2),
                Band::new(12, 3).expect("legal"),
                Band::new(100, -4).expect("legal"),
            ] {
                let mut one = GpuPairs::new();
                let (r, h) = (
                    one.push_read(read, &emission)?,
                    one.push_haplotype(haplotype, Strand::OT, &emission)?,
                );
                one.push_pair(r, h, band)?;
                singles.push(one);
                let (r, h) = (
                    all.push_read(read, &emission)?,
                    all.push_haplotype(haplotype, Strand::OT, &emission)?,
                );
                all.push_pair(r, h, band)?;
            }
        }
    }
    let mut reused = GpuAligner::new(context.clone());
    let batch = reused.align(&all)?;
    for (index, single) in singles.iter().enumerate() {
        let alone = GpuAligner::new(context.clone()).align(single)?;
        let again = reused.align(single)?;
        assert_eq!(bits(&alone), bits(&again), "pair {index}: reused aligner");
        assert_eq!(
            bits(&alone),
            bits(batch.get(index..=index).unwrap_or_default()),
            "pair {index}: in a batch"
        );
    }
    assert_eq!(bits(&reused.align(&all)?), bits(&batch), "the batch again after the singles");
    Ok(())
}

/// A context on a device the caller created -- here with wgpu's default
/// limits and features, as another library sharing its device might --
/// scores what a context of the crate's own does, bit for bit.
#[test]
fn a_context_on_a_borrowed_device_scores_as_its_own() -> Result<(), Box<dyn std::error::Error>> {
    let Some(own) = context() else { return Ok(()) };
    let instance = wgpu::Instance::new(wgpu::InstanceDescriptor::new_without_display_handle());
    let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
        power_preference: wgpu::PowerPreference::HighPerformance,
        ..Default::default()
    }))?;
    let (device, queue) =
        pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor::default()))?;
    let borrowed = GpuContext::from_device(adapter.get_info(), device, queue);

    let groups = tenspeed::groups().map_err(|error| format!("10s: {error:?}"))?;
    let emission = StandardEmission::default();
    let mut pairs = GpuPairs::new();
    for group in &groups {
        let haplotypes = group
            .haplotypes
            .iter()
            .map(|h| pairs.push_haplotype(h, Strand::OT, &emission))
            .collect::<Result<Vec<_>, _>>()?;
        for (read, &offset) in group.reads.iter().zip(&group.offsets) {
            let slot = pairs.push_read(read, &emission)?;
            for &hap_slot in &haplotypes {
                pairs.push_pair(slot, hap_slot, Band::anchored(offset))?;
            }
        }
    }
    let theirs = GpuAligner::new(borrowed).align(&pairs)?;
    let ours = GpuAligner::new(own).align(&pairs)?;
    assert!(!ours.is_empty());
    assert_eq!(bits(&theirs), bits(&ours));
    Ok(())
}

/// No read, no haplotype, a band that misses the haplotype, and an empty
/// launch: `IMPOSSIBLE` or nothing, never a launch error.
#[test]
fn the_degenerate_pairs_are_impossible() -> Result<(), GpuError> {
    let emission = StandardEmission::default();
    let mut pairs = GpuPairs::new();
    let full = pairs.push_read(&read(b"ACGT", Strand::OT), &emission)?;
    let empty = pairs.push_read(&read(b"", Strand::OT), &emission)?;
    let haplotype =
        pairs.push_haplotype(&Haplotype::from_ascii(b"ACGTACGT"), Strand::OT, &emission)?;
    let nothing = pairs.push_haplotype(&Haplotype::new(Vec::new()), Strand::OT, &emission)?;
    pairs.push_pair(empty, haplotype, Band::anchored(0))?;
    pairs.push_pair(full, nothing, Band::anchored(0))?;
    pairs.push_pair(full, haplotype, Band::new(4, 40).expect("legal"))?;
    pairs.push_pair(full, haplotype, Band::new(4, -40).expect("legal"))?;
    let impossible = vec![Log10Likelihood::IMPOSSIBLE; 4];
    assert_eq!(pairs.emulate(Subnormals::Kept), impossible);
    if let Some(context) = context() {
        let mut aligner = GpuAligner::new(context);
        assert_eq!(aligner.align(&pairs)?, impossible);
        assert!(aligner.align(&GpuPairs::new())?.is_empty());
    }
    Ok(())
}

/// The plan refuses what would score the wrong thing: a read paired with a
/// haplotype pushed for the other strand, and a slot it did not hand out.
#[test]
fn the_plan_rejects_mismatched_pairs() -> Result<(), GpuError> {
    let emission = StandardEmission::default();
    let mut pairs = GpuPairs::new();
    let ot = pairs.push_read(&read(b"ACGT", Strand::OT), &emission)?;
    let ob_haplotype =
        pairs.push_haplotype(&Haplotype::from_ascii(b"ACGTACGT"), Strand::OB, &emission)?;
    assert!(matches!(
        pairs.push_pair(ot, ob_haplotype, Band::anchored(0)),
        Err(GpuError::StrandMismatch { read: Strand::OT, haplotype: Strand::OB })
    ));
    let mut other = GpuPairs::new();
    other.push_read(&read(b"ACGT", Strand::OT), &emission)?;
    let second = other.push_read(&read(b"ACGT", Strand::OT), &emission)?;
    let hap = other.push_haplotype(&Haplotype::from_ascii(b"ACGT"), Strand::OT, &emission)?;
    other.push_pair(second, hap, Band::anchored(0))?;
    assert!(matches!(
        pairs.push_pair(second, ob_haplotype, Band::anchored(0)),
        Err(GpuError::UnknownSlot)
    ));
    assert_eq!(pairs.len(), 0);
    Ok(())
}
