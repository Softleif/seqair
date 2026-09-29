//! An external oracle: rust-bio's `bio::stats::pairhmm::PairHMM`.
//!
//! rust-bio implements the textbook (Durbin) pair-HMM, not GATK's, so the two
//! only line up once the parameters are mapped. The mapping below is *exact* --
//! transcribing rust-bio's recurrence into linear-space `f64` reproduces
//! `align_full`'s number to the last bit, plus one extra term named below --
//! and it is worth writing down, because every part of it is a place the two
//! conventions could have been confused.
//!
//! # The mapping
//!
//! rust-bio's `x` is the **haplotype** and its `y` is the **read**. That
//! direction is forced: `free_start_gap_x` and `free_end_gap_x` are the only
//! free ends rust-bio offers, and GATK's free ends are both on the haplotype
//! (the read may start anywhere on it, and the total is taken at the last read
//! base over every haplotype column).
//!
//! With `x` the haplotype, rust-bio's two gap states land as:
//!
//! | rust-bio | consumes | compair | opened with | extended with | left with |
//! | --- | --- | --- | --- | --- | --- |
//! | `fx` ("gap in y") | `x` | deletion | `prob_gap_y` | `prob_gap_y_extend` | `1 - prob_gap_x_extend` |
//! | `fy` ("gap in x") | `y` | insertion | `prob_gap_x` | `prob_gap_x_extend` | `1 - prob_gap_y_extend` |
//!
//! Note the crossed names in the last column: rust-bio leaves a gap state
//! through the *other* state's extension probability. That is harmless here
//! only because GATK has a single gap-continuation penalty, so
//! `prob_gap_x_extend == prob_gap_y_extend == g` and both exits are `1 - g` --
//! which is exactly why this file compares only reads with **uniform**
//! insertion, deletion and gap-continuation qualities. Per-base *base*
//! qualities are fine; they live in the emission, which is per cell.
//!
//! So `prob_gap_x = p_ins`, `prob_gap_y = p_del`, both extensions are `g`, and
//! rust-bio's `prob_no_gap = 1 - (p_ins + p_del)` is GATK's clamped
//! `match_to_match` verbatim. The match emission is `1 - eps` on a match and
//! `eps / 3` on a mismatch, with an unknown base on either side matching at
//! `1 - eps`; the gap states emit with probability one.
//!
//! # The free start
//!
//! GATK spreads the start over the haplotype as `D[0][j] = 1 / h`, and every
//! path leaves it through one `indel_to_match = 1 - g` factor. rust-bio instead
//! parks the start in the *match* row (`fm[i][0] = prob_start_gap_x(i)`) and
//! every path leaves it through one `prob_no_gap` factor. Equal weights, then,
//! need
//!
//! ```text
//! prob_start_gap_x = (1 - g) / (h * (1 - p_ins - p_del))
//! ```
//!
//! and no post-hoc normalisation at all. Two wrinkles force the **dummy
//! haplotype base** this file prepends:
//!
//! - rust-bio seeds `fm[0][0] = 1` *before* the loop and then adds
//!   `prob_start_gap_x(0)` on top, so offset zero would be weighted
//!   `1 + start` while every other offset is weighted `start`. Giving offset
//!   zero a base that cannot be emitted (`prob_emit_xy(0, _) = 0`) discards
//!   that anomalous start and leaves exactly `h` starts of weight `start`.
//! - rust-bio clamps its result at probability one. Carrying the `1 / h`
//!   through the start keeps the total a true probability, clear of the clamp.
//!
//! # The free end, and the one term that does not map
//!
//! `free_end_gap_x` sums `fm`, `fx` and `fy` at the last `y` column over every
//! `x` row: the last read base, over every haplotype position. That is GATK's
//! termination -- except that GATK sums the match and insertion matrices only,
//! and rust-bio also adds `fx`, the deletion matrix. A read may therefore end
//! inside a deletion in rust-bio's total but not in compair's, so
//! `rust_bio - compair` is `log10(1 + S_D / (S_M + S_I))`, which the deletion
//! recurrence bounds rigorously:
//!
//! ```text
//! S_D * (1 - g) = p_del * (S_M - M_h) - g * D_h <= p_del * S_M
//! ```
//!
//! so the difference lies in `[0, log10(1 + p_del / (1 - g))]`. That is
//! [`deletion_matrix_bound`], asserted below, and at Q45 / Q10 it is 1.5e-5.
//!
//! # Why the agreement is not exact
//!
//! rust-bio's `ln_sum3_exp_approx` sorts its three arguments with
//!
//! ```text
//! if p1 < p2 { swap(p1, p2) }   // p1 = max(p1, p2)
//! if p1 > p0 { swap(p1, p0) }   // p0 = max of all three
//! if p0 - p1 > 10.0 { return p0 }
//! ```
//!
//! When the match term `p0` is the *smallest* of the three, the second swap
//! leaves `p1` holding it, so the shortcut compares the maximum against the
//! **minimum** and discards `p2` -- which is then the second largest, and can
//! be within a factor of one of the maximum. On reads that need a gap this
//! silently drops real probability mass: measured at up to 0.25 in log10 over
//! random pairs, but only 1.7e-5 over the GATK vectors, whose reads align
//! cleanly enough that the match term always dominates. (`fastexp`, rust-bio's
//! polynomial `exp`, contributes a further ~2e-7 and is not the limit.)
//!
//! The shortcut only ever *removes* mass, so the upper bound above stays
//! rigorous and is asserted everywhere; the lower bound is empirical, and tight
//! only where the reads align cleanly.

use bio::stats::{
    LogProb, Prob,
    pairhmm::{EmissionParameters, GapParameters, PairHMM, StartEndGapParameters, XYEmission},
};
use compair::{Base, BaseQuality, Haplotype, Read, StandardEmission, Strand, align_full};
use hegel::TestCase;
use hegel::generators::{self as gs, Generator, PrintableGenerator};

/// Floating-point slack on top of [`deletion_matrix_bound`], which is otherwise
/// exact. Measured worst overshoot over 20 000 random pairs: 1.1e-6.
const UPPER_SLACK: f64 = 1e-5;

/// How far *below* `align_full` rust-bio may land on the GATK vectors at
/// Q45 / Q45 / Q10, where its `ln_sum3_exp_approx` shortcut only ever discards
/// terms that really are negligible. Measured worst case over all 104: 1.7e-5.
const GATK_LOWER_SLACK: f64 = 1e-4;

/// The same, for random pairs whose read is an exact cut of the haplotype but
/// whose indel qualities go as low as Q10. Measured worst case over 20 000:
/// 3.2e-3.
const UNEDITED_LOWER_SLACK: f64 = 2e-2;

/// The same, for reads carrying indels, where the shortcut's mis-sorted
/// comparison can discard a near-maximal term. Measured worst case over 20 000
/// random pairs: 0.25.
const EDITED_LOWER_SLACK: f64 = 0.5;

const LN_10: f64 = core::f64::consts::LN_10;

const DATA: &str = include_str!("data/pairhmm-testdata.txt");

/// `10^(-Q/10)`.
fn error_probability(phred: u8) -> f64 {
    10f64.powf(-f64::from(phred) / 10.0)
}

/// The most, in log10, that rust-bio's inclusion of the deletion matrix in the
/// free end can raise its total above `align_full`'s.
fn deletion_matrix_bound(deletion_qual: u8, gap_qual: u8) -> f64 {
    let (p_del, g) = (error_probability(deletion_qual), error_probability(gap_qual));
    (1.0 + p_del / (1.0 - g)).log10()
}

/// `x` is the haplotype with one unemittable base prepended, `y` is the read.
struct Emissions<'a> {
    haplotype: &'a [Base],
    read: &'a [Base],
    epsilon: &'a [f64],
}

impl EmissionParameters for Emissions<'_> {
    fn prob_emit_xy(&self, i: usize, j: usize) -> XYEmission {
        let Some(&hap) = i.checked_sub(1).and_then(|index| self.haplotype.get(index)) else {
            // The dummy base. Nothing can be emitted against it, which is what
            // discards rust-bio's double-weighted start at offset zero.
            return XYEmission::Mismatch(LogProb::ln_zero());
        };
        let (Some(&base), Some(&eps)) = (self.read.get(j), self.epsilon.get(j)) else {
            return XYEmission::Mismatch(LogProb::ln_zero());
        };
        if hap == Base::Unknown || base == Base::Unknown || hap == base {
            XYEmission::Match(LogProb::from(Prob(1.0 - eps)))
        } else {
            XYEmission::Mismatch(LogProb::from(Prob(eps / 3.0)))
        }
    }

    /// GATK's gap states emit with probability one.
    fn prob_emit_x(&self, _i: usize) -> LogProb {
        LogProb::ln_one()
    }

    fn prob_emit_y(&self, _j: usize) -> LogProb {
        LogProb::ln_one()
    }

    fn len_x(&self) -> usize {
        self.haplotype.len().saturating_add(1)
    }

    fn len_y(&self) -> usize {
        self.read.len()
    }
}

struct Gaps {
    insertion: f64,
    deletion: f64,
    continuation: f64,
}

impl GapParameters for Gaps {
    /// Opens `fy`, the state that consumes a read base: `match_to_insertion`.
    fn prob_gap_x(&self) -> LogProb {
        LogProb::from(Prob(self.insertion))
    }

    /// Opens `fx`, the state that consumes a haplotype base: `match_to_deletion`.
    fn prob_gap_y(&self) -> LogProb {
        LogProb::from(Prob(self.deletion))
    }

    fn prob_gap_x_extend(&self) -> LogProb {
        LogProb::from(Prob(self.continuation))
    }

    fn prob_gap_y_extend(&self) -> LogProb {
        LogProb::from(Prob(self.continuation))
    }
}

struct Mode {
    start: LogProb,
}

impl StartEndGapParameters for Mode {
    fn prob_start_gap_x(&self, _i: usize) -> LogProb {
        self.start
    }

    fn free_start_gap_x(&self) -> bool {
        true
    }

    fn free_end_gap_x(&self) -> bool {
        true
    }
}

/// rust-bio's answer to the same problem, in log10 and already on
/// `align_full`'s scale: there is no normalisation left to apply.
fn rust_bio_log10(
    haplotype: &Haplotype,
    read: &Read,
    insertion_qual: u8,
    deletion_qual: u8,
    gap_qual: u8,
) -> f64 {
    let epsilon: Vec<f64> =
        read.base_quals().iter().map(|q| error_probability(q.as_byte())).collect();
    let gaps = Gaps {
        insertion: error_probability(insertion_qual),
        deletion: error_probability(deletion_qual),
        continuation: error_probability(gap_qual),
    };
    let match_to_match = 1.0 - (gaps.insertion + gaps.deletion);
    #[allow(
        clippy::cast_precision_loss,
        reason = "a haplotype long enough to lose precision here does not exist"
    )]
    let h = haplotype.len() as f64;
    let start = LogProb(((1.0 - gaps.continuation) / (h * match_to_match)).ln());
    let emissions =
        Emissions { haplotype: haplotype.bases(), read: read.bases(), epsilon: &epsilon };
    let mut hmm = PairHMM::new(&gaps);
    *hmm.prob_related(&emissions, &Mode { start }, None) / LN_10
}

/// The GATK vectors' haplotype, read and base-quality columns. The per-base
/// indel columns are dropped: rust-bio's gap parameters are not per base, so
/// the comparison is only meaningful with uniform ones.
fn gatk_pairs() -> Vec<(Haplotype, Vec<Base>, Vec<BaseQuality>)> {
    DATA.lines()
        .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
        .filter_map(|line| {
            let fields: Vec<&str> = line.split_whitespace().collect();
            let [haplotype, read, base_quals, ..] = fields[..] else { return None };
            Some((
                Haplotype::from_ascii(haplotype.as_bytes()),
                read.bytes().map(Base::from).collect(),
                base_quals
                    .bytes()
                    .map(|byte| BaseQuality::from_byte(byte.saturating_sub(33)))
                    .collect(),
            ))
        })
        .collect()
}

fn uniform_read(
    bases: Vec<Base>,
    base_quals: &[BaseQuality],
    insertion_qual: u8,
    deletion_qual: u8,
    gap_qual: u8,
) -> Option<Read> {
    Read::uniform(
        bases,
        base_quals,
        BaseQuality::from_byte(insertion_qual),
        BaseQuality::from_byte(deletion_qual),
        BaseQuality::from_byte(gap_qual),
        Strand::OT,
    )
    .ok()
}

/// The 104 GATK vectors, re-scored with uniform Q45 / Q45 / Q10 indel
/// qualities so that rust-bio's single gap-continuation parameter is faithful
/// to the read.
///
/// These reads align cleanly, so rust-bio's `ln_sum3_exp_approx` only ever
/// discards terms below `e^-10` of the maximum and the two agree to 1.7e-5.
#[test]
fn gatk_vectors_agree_with_rust_bio() {
    let (insertion_qual, deletion_qual, gap_qual) = (45u8, 45u8, 10u8);
    let bound = deletion_matrix_bound(deletion_qual, gap_qual);
    let pairs = gatk_pairs();
    assert_eq!(pairs.len(), 104, "the GATK vectors must all parse");
    for (index, (haplotype, bases, base_quals)) in pairs.into_iter().enumerate() {
        let read = uniform_read(bases, &base_quals, insertion_qual, deletion_qual, gap_qual)
            .expect("the vector's read must build");
        let mine = align_full(&haplotype, &read, &StandardEmission::default()).get();
        let theirs = rust_bio_log10(&haplotype, &read, insertion_qual, deletion_qual, gap_qual);
        let delta = theirs - mine;
        assert!(
            delta <= bound + UPPER_SLACK,
            "vector {index}: rust-bio {theirs} is {delta} above compair {mine}, more than the \
             deletion-matrix bound {bound}"
        );
        assert!(
            delta >= -GATK_LOWER_SLACK,
            "vector {index}: rust-bio {theirs} is {delta} below compair {mine}"
        );
    }
}

/// A haplotype, a read cut out of it and perturbed, and one uniform indel
/// quality set.
#[derive(Debug, Clone)]
struct Case {
    haplotype: Haplotype,
    read: Read,
    insertion_qual: u8,
    deletion_qual: u8,
    gap_qual: u8,
    edited: bool,
}

/// One of the four nucleotides, or one time in ten an `N`.
#[hegel::composite]
fn any_byte_or_n(tc: &TestCase) -> u8 {
    if tc.draw_silent(gs::weighted_booleans(0.1)) {
        b'N'
    } else {
        tc.draw_silent(gs::sampled_from(b"ACGT"))
    }
}

/// A haplotype, a read cut out of it at a random start, perturbed by at most
/// `max_edits` edits, and a random uniform indel quality set.
#[hegel::composite]
fn any_case_inner(tc: &TestCase, max_edits: usize) -> Case {
    // Lengths are drawn as integers, which hegel spreads over the range, and
    // not as a vector's size, which it keeps short.
    let len = tc.draw_silent(gs::integers::<usize>().min_value(20).max_value(119));
    let hap: Vec<u8> = tc.draw_silent(gs::vecs(any_byte_or_n()).min_size(len).max_size(len));
    let haplotype = Haplotype::from_ascii(&hap);
    let start = tc.draw_silent(gs::integers::<usize>().max_value(119));
    let start = start.min(haplotype.len().saturating_sub(10));
    let want = tc.draw_silent(gs::integers::<usize>().min_value(20).max_value(89));
    let end = start.saturating_add(want).min(haplotype.len());
    let Some(cut) = haplotype.bases().get(start..end) else { tc.reject() };
    let mut bases = cut.to_vec();
    let mut edited = false;
    for _ in 0..tc.draw_silent(gs::integers::<usize>().max_value(max_edits)) {
        if bases.len() < 10 {
            break;
        }
        let at = tc.draw_silent(gs::integers::<usize>().max_value(bases.len() - 1));
        edited = true;
        match tc.draw_silent(gs::integers::<u8>().max_value(2)) {
            0 => {
                let base = Base::from(tc.draw_silent(gs::sampled_from(b"ACGT")));
                if let Some(slot) = bases.get_mut(at) {
                    *slot = base;
                }
            }
            1 => bases.insert(at, Base::from(tc.draw_silent(gs::sampled_from(b"ACGT")))),
            _ => {
                let _removed = bases.remove(at);
            }
        }
    }
    tc.assume(bases.len() >= 10);
    let quals: Vec<u8> = tc.draw_silent(
        gs::vecs(gs::integers::<u8>().min_value(2).max_value(45))
            .min_size(bases.len())
            .max_size(bases.len()),
    );
    let base_quals: Vec<BaseQuality> = quals.into_iter().map(BaseQuality::from_byte).collect();
    // Q10 is the floor that keeps `p_ins + p_del < 1`, without which
    // rust-bio's `prob_no_gap` does not exist.
    let insertion_qual = tc.draw_silent(gs::integers::<u8>().min_value(10).max_value(45));
    let deletion_qual = tc.draw_silent(gs::integers::<u8>().min_value(10).max_value(45));
    let gap_qual = tc.draw_silent(gs::integers::<u8>().min_value(2).max_value(25));
    let Some(read) = uniform_read(bases, &base_quals, insertion_qual, deletion_qual, gap_qual)
    else {
        tc.reject()
    };
    Case { haplotype, read, insertion_qual, deletion_qual, gap_qual, edited }
}

fn any_case(max_edits: usize) -> impl PrintableGenerator<Case> {
    any_case_inner(max_edits).print_as_debug()
}

fn check(tc: &TestCase, case: &Case) {
    let Case { haplotype, read, insertion_qual, deletion_qual, gap_qual, edited } = case;
    let mine = align_full(haplotype, read, &StandardEmission::default()).get();
    let theirs = rust_bio_log10(haplotype, read, *insertion_qual, *deletion_qual, *gap_qual);
    tc.assume(mine.is_finite() && theirs.is_finite());
    let bound = deletion_matrix_bound(*deletion_qual, *gap_qual);
    let delta = theirs - mine;
    // Rigorous: the deletion matrix in rust-bio's free end is the only thing
    // that can put it above `align_full`, and the deletion recurrence bounds it.
    assert!(
        delta <= bound + UPPER_SLACK,
        "rust-bio {theirs} is {delta} above compair {mine}, more than {bound}"
    );
    let lower = if *edited { EDITED_LOWER_SLACK } else { UNEDITED_LOWER_SLACK };
    assert!(delta >= -lower, "rust-bio {theirs} is {delta} below compair {mine}");
}

/// Reads cut straight out of the haplotype. The match term dominates
/// almost every three-way sum rust-bio takes a shortcut on, so the two
/// agree to the deletion-matrix bound above and 2e-2 below.
#[hegel::test]
fn unedited_reads_agree_with_rust_bio(tc: TestCase) {
    let case = tc.draw(any_case(0));
    check(&tc, &case);
}

/// Reads carrying substitutions and short indels. The upper bound is still
/// rigorous; the lower one is loose because rust-bio's mis-sorted
/// three-way shortcut discards gap mass on exactly these reads.
#[hegel::test]
fn edited_reads_agree_with_rust_bio(tc: TestCase) {
    let case = tc.draw(any_case(3));
    check(&tc, &case);
}
