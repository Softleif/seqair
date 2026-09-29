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

#[path = "../src/pinned.rs"]
mod pinned;
#[allow(dead_code, reason = "each integration test uses a different part of this")]
mod support;

#[path = "support/tenspeed.rs"]
mod tenspeed;

use std::sync::{Arc, OnceLock};

use compair::{
    Band, Base, BaseQuality, Betas, Emission, Haplotype, Log10Likelihood, Probability, Read,
    StandardEmission, Strand, TapsEmission, Workspace, align_full, align_strips,
    gpu::{GpuAligner, GpuContext, GpuError, GpuPairs, Subnormals},
    trusted,
};
use hegel::TestCase;
use hegel::generators::{self as gs, Generator, PrintableGenerator};
use support::{Case, any_conversion, any_probability, arbitrary_case, derived_case, length};

/// `|Δ log10|` the GPU may differ from the strip kernel by: the GATK vectors'
/// tolerance. Measured: 6e-7 at worst on the 10s pairs.
const TOLERANCE: f64 = 1e-4;

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
#[hegel::composite]
fn banded_case_inner(tc: &TestCase) -> (Case, Band) {
    let case = tc.draw_silent(arbitrary_case());
    let width = tc.draw_silent(gs::integers::<u32>().min_value(2).max_value(160));
    let Ok(band) = Band::new(width, case.offset) else { tc.reject() };
    (case, band)
}

fn banded_case() -> impl PrintableGenerator<(Case, Band)> {
    banded_case_inner().print_as_debug()
}

/// A GPU launch's scores, each passed through [`trusted`] as a caller does.
type Rescue<'a> = dyn Fn(&[Log10Likelihood]) -> Vec<Log10Likelihood> + 'a;

/// Each score through [`trusted`] with its case's pair, as a GPU caller does.
fn rescued<'c>(
    scores: &[Log10Likelihood],
    cases: &'c [(Case, Band)],
    emission: impl Fn(&'c Case) -> Box<dyn Emission + 'c>,
) -> Vec<Log10Likelihood> {
    scores
        .iter()
        .zip(cases)
        .map(|(score, (case, band))| {
            trusted(*score, &case.haplotype, &case.read, &*emission(case), *band)
        })
        .collect()
}

/// The emissions the parity gates sweep, for one batch of cases: the plan,
/// the strip kernel's scores, and the rescue a GPU caller applies.
fn for_each_emission(
    cases: &[(Case, Band)],
    conversion: compair::ConversionModel,
    uniform: Probability,
    mut check: impl FnMut(&str, &GpuPairs, Vec<Log10Likelihood>, &Rescue<'_>),
) {
    let mut workspace = Workspace::new();
    {
        let emission = StandardEmission::default();
        let strips = cases
            .iter()
            .map(|(case, band)| {
                workspace.align_strips(&case.haplotype, &case.read, &emission, *band)
            })
            .collect();
        let rescue = |scores: &[Log10Likelihood]| {
            rescued(scores, cases, |_| Box::new(StandardEmission::default()))
        };
        check("standard", &plan(cases, &emission).expect("the plan builds"), strips, &rescue);
    }
    {
        // Betas per site differ per case, so each case gets its own emission;
        // the plan is built pair by pair for the same reason.
        let mut pairs = GpuPairs::new();
        let mut strips = Vec::new();
        for (case, band) in cases {
            let emission = TapsEmission::new(conversion, Betas::PerSite(&case.betas));
            let read = pairs.push_read(&case.read, &emission).expect("the read pushes");
            let haplotype = pairs
                .push_haplotype(&case.haplotype, case.read.strand(), &emission)
                .expect("the haplotype pushes");
            pairs.push_pair(read, haplotype, *band).expect("the pair pushes");
            strips.push(workspace.align_strips(&case.haplotype, &case.read, &emission, *band));
        }
        let rescue = |scores: &[Log10Likelihood]| {
            rescued(scores, cases, |case| {
                Box::new(TapsEmission::new(conversion, Betas::PerSite(&case.betas)))
            })
        };
        check("taps", &pairs, strips, &rescue);
    }
    {
        let emission = TapsEmission::new(conversion, Betas::Uniform(uniform));
        let strips = cases
            .iter()
            .map(|(case, band)| {
                workspace.align_strips(&case.haplotype, &case.read, &emission, *band)
            })
            .collect();
        let rescue = |scores: &[Log10Likelihood]| {
            rescued(scores, cases, |_| {
                Box::new(TapsEmission::new(conversion, Betas::Uniform(uniform)))
            })
        };
        check("uniform", &plan(cases, &emission).expect("the plan builds"), strips, &rescue);
    }
}

/// The algorithm gate, CPU only: the shader's transcription with
/// subnormals kept is the strip kernel, bit for bit, at every width and
/// offset -- including bands that miss the alignment or the haplotype,
/// `N` on both sides, and every emission.
#[hegel::test(test_cases = pinned::cases(512))]
fn the_kernel_algorithm_is_the_strip_kernel(tc: TestCase) {
    let count = tc.draw(length(1, 5));
    let cases: Vec<_> = (0..count).map(|_| tc.draw(banded_case())).collect();
    let conversion = tc.draw(any_conversion());
    let uniform = tc.draw(any_probability());
    for_each_emission(&cases, conversion, uniform, |name, pairs, strips, rescue| {
        assert_eq!(bits(&rescue(&pairs.emulate(Subnormals::Kept))), bits(&strips), "{}", name);
    });
}

/// Plans appended one after another are the plan built in one go: the
/// same bytes uploaded, the same scores (through the transcription, so no
/// GPU is needed), and each part's scores at the range `append` named.
#[hegel::test(test_cases = pinned::cases(256))]
fn appended_plans_are_one_plan(tc: TestCase) {
    let count = tc.draw(length(1, 23));
    let cases: Vec<_> = (0..count).map(|_| tc.draw(banded_case())).collect();
    let cut_count = tc.draw(length(0, 3));
    let cuts: Vec<usize> = (0..cut_count).map(|_| tc.draw(length(0, 23))).collect();
    let emission = StandardEmission::default();
    let whole = plan(&cases, &emission).expect("the plan builds");
    let mut cuts: Vec<usize> = cuts.into_iter().map(|cut| cut.min(cases.len())).collect();
    cuts.push(0);
    cuts.push(cases.len());
    cuts.sort_unstable();
    let mut appended = GpuPairs::new();
    let mut ranges = Vec::new();
    for window in cuts.windows(2) {
        let &[from, to] = window else { continue };
        let part =
            plan(cases.get(from..to).unwrap_or_default(), &emission).expect("the plan builds");
        let range = appended.append(&part).expect("the plans append");
        assert_eq!(range.clone(), from..to);
        ranges.push((range, bits(&part.emulate(Subnormals::Kept))));
    }
    assert_eq!(appended.len(), whole.len());
    assert_eq!(appended.upload_bytes(), whole.upload_bytes());
    let scores = bits(&appended.emulate(Subnormals::Kept));
    assert_eq!(&scores, &bits(&whole.emulate(Subnormals::Kept)));
    for (range, part) in ranges {
        assert_eq!(scores.get(range).unwrap_or_default(), part.as_slice());
    }
}

/// The device gate: every score is the transcription's, with subnormal
/// intermediates flushed (Metal, RADV) or kept (an IEEE-denormal GPU).
/// Bit-exact, not a tolerance: this is what makes the kernel's numbers
/// reproducible from one GPU to another, and it needs contraction
/// blocked, which the default kernel does.
#[hegel::test(test_cases = pinned::cases(64))]
fn the_gpu_runs_the_transcribed_kernel(tc: TestCase) {
    let count = tc.draw(length(1, 47));
    let cases: Vec<_> = (0..count).map(|_| tc.draw(banded_case())).collect();
    let conversion = tc.draw(any_conversion());
    let uniform = tc.draw(any_probability());
    let Some(context) = context() else { return };
    let mut aligner = GpuAligner::new(context);
    for_each_emission(&cases, conversion, uniform, |name, pairs, _, _| {
        let gpu = bits(&aligner.align(pairs).expect("the launch runs"));
        let flushed = bits(&pairs.emulate(Subnormals::Flushed));
        let kept = bits(&pairs.emulate(Subnormals::Kept));
        for (index, &score) in gpu.iter().enumerate() {
            assert!(
                Some(&score) == flushed.get(index) || Some(&score) == kept.get(index),
                "{}: pair {} scored {} on the GPU, {:?} flushed, {:?} kept",
                name,
                index,
                f64::from_bits(score),
                flushed.get(index).map(|b| f64::from_bits(*b)),
                kept.get(index).map(|b| f64::from_bits(*b)),
            );
        }
    });
}

/// The tolerance oracle against the CPU's strip kernel.
///
/// Within [`TOLERANCE`] everywhere once each score has been through
/// [`trusted`], and `IMPOSSIBLE` wherever the strip kernel is. A GPU that
/// flushes subnormal intermediates loses more than the CPU kernel only
/// where a pair scores so low that its surviving path sat far under each
/// row's maximum -- which is below the floor `trusted` rescores under.
#[hegel::test(test_cases = pinned::cases(64))]
fn the_gpu_is_the_strip_kernel_within_tolerance(tc: TestCase) {
    let count = tc.draw(length(1, 47));
    let cases: Vec<_> = (0..count).map(|_| tc.draw(banded_case())).collect();
    let conversion = tc.draw(any_conversion());
    let uniform = tc.draw(any_probability());
    let Some(context) = context() else { return };
    let mut aligner = GpuAligner::new(context);
    for_each_emission(&cases, conversion, uniform, |name, pairs, strips, rescue| {
        let gpu = rescue(&aligner.align(pairs).expect("the launch runs"));
        for (index, (g, s)) in gpu.iter().zip(&strips).enumerate() {
            let (g, s) = (g.get(), s.get());
            if !s.is_finite() {
                assert_eq!(g, s, "{}: pair {}", name, index);
                continue;
            }
            assert!(
                !g.is_nan() && g <= s + TOLERANCE,
                "{}: pair {}: gpu {} above strips {}",
                name,
                index,
                g,
                s
            );
            assert!((g - s).abs() <= TOLERANCE, "{}: pair {}: gpu {} strips {}", name, index, g, s);
        }
    });
}

/// The band contract on the GPU: a band never scores above the unbanded
/// `f64` reference.
#[hegel::test(test_cases = pinned::cases(64))]
fn the_gpu_band_never_beats_the_reference(tc: TestCase) {
    let count = tc.draw(length(1, 31));
    let cases: Vec<_> = (0..count).map(|_| tc.draw(derived_case(4))).collect();
    let Some(context) = context() else { return };
    let emission = StandardEmission::default();
    let banded: Vec<(Case, Band)> = cases.iter().map(|case| (case.clone(), case.band())).collect();
    let pairs = plan(&banded, &emission).expect("the launch runs");
    let gpu = GpuAligner::new(context).align(&pairs).expect("the launch runs");
    for (case, score) in cases.iter().zip(&gpu) {
        let full = align_full(&case.haplotype, &case.read, &emission).get();
        assert!(!score.get().is_nan());
        assert!(
            score.get() <= full + TOLERANCE,
            "banded {} beats the reference {}",
            score.get(),
            full
        );
    }
}

/// The 10s benchmark in one launch: every pair within tolerance of the strip
/// kernel. Prints how many are bit-identical.
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
