use compair::{
    Base, BaseQuality, Error, HapSite, Haplotype, Observation, Read, Strand, error_probability,
};

/// The strand a mirrored case reports. `compair::Strand` has no `mirror`:
/// `Unknown` has no mirror image, and a `Read` cannot carry it anyway.
pub fn mirror_strand(strand: Strand) -> Strand {
    match strand {
        Strand::OT => Strand::OB,
        Strand::OB => Strand::OT,
        Strand::Unknown => Strand::Unknown,
    }
}

/// One row of GATK's published pair-HMM test data.
pub struct Vector {
    pub haplotype: Haplotype,
    pub read: Read,
    pub expected_log10: f64,
}

const DATA: &str = include_str!("../data/pairhmm-testdata.txt");

fn quals(field: &str) -> Vec<BaseQuality> {
    field.bytes().map(|byte| BaseQuality::from_byte(byte.saturating_sub(33))).collect()
}

pub fn gatk_vectors() -> Result<Vec<Vector>, Error> {
    let mut out = Vec::new();
    for line in DATA.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        let [haplotype, read, base_q, ins_q, del_q, gcp, expected] = fields[..] else {
            continue;
        };
        out.push(Vector {
            haplotype: Haplotype::from_ascii(haplotype.as_bytes()),
            read: Read::new(
                bases(read),
                &quals(base_q),
                &quals(ins_q),
                &quals(del_q),
                &quals(gcp),
                Strand::OT,
            )?,
            expected_log10: expected.parse().unwrap_or(f64::NAN),
        });
    }
    Ok(out)
}

pub fn bases(sequence: &str) -> Vec<Base> {
    sequence.bytes().map(Base::from).collect()
}

/// The haplotype offset at which the read has the most matching bases.
///
/// Test scaffolding only: the band needs an anchor, and a caller in a variant
/// caller has one from the alignment. This stands in for it.
#[allow(
    clippy::cast_possible_wrap,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "test scaffolding over sequences of a few hundred bases"
)]
pub fn seed_offset(haplotype: &Haplotype, read: &Read) -> i32 {
    let (h, r) = (haplotype.len() as i64, read.len() as i64);
    let mut best = (0i64, -1i64);
    for offset in -(r - 1)..h {
        let mut hits = 0i64;
        for (index, base) in read.bases().iter().enumerate() {
            let column = offset + index as i64;
            if column >= 0 && haplotype.bases().get(column as usize) == Some(base) {
                hits += 1;
            }
        }
        if hits > best.1 {
            best = (offset, hits);
        }
    }
    i32::try_from(best.0).unwrap_or(0)
}

use compair::{Band, ConversionModel, Probability};
use hegel::TestCase;
use hegel::generators::{self as gs, Generator, PrintableGenerator};

/// `lo..=hi`, uniformly.
///
/// Every length below is drawn this way and the vector then built to it:
/// hegel spreads an integer over its range, but keeps a vector's own size
/// short, and the long reads are where the DP has room to go wrong.
pub fn length(lo: usize, hi: usize) -> impl PrintableGenerator<usize> {
    gs::integers::<usize>().min_value(lo).max_value(hi)
}

/// One of the four nucleotides.
pub fn any_base() -> impl PrintableGenerator<Base> {
    gs::sampled_from(&[Base::A, Base::C, Base::G, Base::T]).print_as_debug()
}

/// A nucleotide, or every so often an `N`: the kernels have a code for it and
/// the reference has a branch, and the two have to agree through the whole
/// DP, not only per cell.
#[hegel::composite]
pub fn any_base_or_n(tc: &TestCase) -> Base {
    if tc.draw_silent(gs::weighted_booleans(1.0 / 9.0)) {
        Base::Unknown
    } else {
        tc.draw_silent(any_base())
    }
}

pub fn any_strand() -> impl PrintableGenerator<Strand> {
    gs::sampled_from(&[Strand::OT, Strand::OB]).print_as_debug()
}

fn unit_interval() -> impl Generator<f64> {
    gs::floats::<f64>().min_value(0.0).max_value(1.0)
}

fn probability(value: f64) -> Probability {
    Probability::new(value).unwrap_or(Probability::ZERO)
}

pub fn any_probability() -> impl PrintableGenerator<Probability> {
    unit_interval().map(probability).print_as_debug()
}

pub fn any_conversion() -> impl PrintableGenerator<ConversionModel> {
    hegel::tuples!(unit_interval(), unit_interval())
        .map(|(c, f)| ConversionModel::new(probability(c), probability(f)))
        .print_as_debug()
}

/// A conversion model a real library could have: efficient, and rarely wrong
/// about an unmethylated base.
///
/// `any_conversion` includes `c = 0, f = 1` -- "every cytosine converts, none
/// of the methylated ones" -- under which a read is *unlikely* against the very
/// haplotype it was cut from, and alignments elsewhere become competitive. Any
/// property whose premise is "the optimal path is near the seed offset" needs
/// this one instead.
pub fn plausible_conversion() -> impl PrintableGenerator<ConversionModel> {
    hegel::tuples!(
        gs::floats::<f64>().min_value(0.5).max_value(1.0),
        gs::floats::<f64>().min_value(0.0).max_value(0.05),
    )
    .map(|(c, f)| {
        ConversionModel::new(
            Probability::new(c).unwrap_or(Probability::ONE),
            Probability::new(f).unwrap_or(Probability::ZERO),
        )
    })
    .print_as_debug()
}

/// A haplotype, a read derived from it, and the offset the read starts at.
#[derive(Debug, Clone)]
pub struct Case {
    pub haplotype: Haplotype,
    pub read: Read,
    pub offset: i32,
    pub betas: Vec<Probability>,
}

impl Case {
    pub fn band(&self) -> Band {
        Band::anchored(self.offset)
    }
}

/// `len` draws from `element`.
fn exactly<T>(tc: &TestCase, element: impl Generator<T>, len: usize) -> Vec<T> {
    tc.draw_silent(gs::vecs(element).min_size(len).max_size(len))
}

/// `len` qualities in `lo..=hi`.
fn draw_quals(tc: &TestCase, lo: u8, hi: u8, len: usize) -> Vec<BaseQuality> {
    exactly(tc, gs::integers::<u8>().min_value(lo).max_value(hi), len)
        .into_iter()
        .map(BaseQuality::from_byte)
        .collect()
}

/// A read cut out of a haplotype and then perturbed by at most `max_edits`
/// point substitutions and one-base indels, so the optimal alignment stays
/// within a few columns of the offset the case reports and a default band
/// contains it by construction.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    reason = "test scaffolding over sequences of a few hundred bases"
)]
#[hegel::composite]
fn derived_case_inner(tc: &TestCase, max_edits: usize) -> Case {
    let len = tc.draw_silent(length(80, 159));
    let hap = exactly(tc, any_base_or_n(), len);
    let start = tc.draw_silent(length(0, 24)).min(hap.len().saturating_sub(10));
    let want = tc.draw_silent(length(40, 89));
    let end = (start + want).min(hap.len());
    let Some(cut) = hap.get(start..end) else { tc.reject() };
    let mut bases = cut.to_vec();
    for _ in 0..tc.draw_silent(length(0, max_edits)) {
        if bases.len() < 8 {
            break;
        }
        let at = tc.draw_silent(length(0, bases.len() - 1));
        match tc.draw_silent(gs::integers::<u8>().max_value(2)) {
            0 => {
                let to = tc.draw_silent(any_base_or_n());
                if let Some(slot) = bases.get_mut(at) {
                    *slot = to;
                }
            }
            1 => bases.insert(at, tc.draw_silent(any_base_or_n())),
            _ => {
                let _removed = bases.remove(at);
            }
        }
    }
    tc.assume(bases.len() >= 8);
    let quals = draw_quals(tc, 2, 45, bases.len());
    let strand = tc.draw_silent(any_strand());
    let Ok(read) = Read::uniform(
        bases,
        &quals,
        BaseQuality::from_byte(45),
        BaseQuality::from_byte(45),
        BaseQuality::from_byte(10),
        strand,
    ) else {
        tc.reject()
    };
    let betas = exactly(tc, any_probability(), hap.len());
    Case { haplotype: Haplotype::new(hap), read, offset: start as i32, betas }
}

pub fn derived_case(max_edits: usize) -> impl PrintableGenerator<Case> {
    derived_case_inner(max_edits).print_as_debug()
}

/// An unconstrained pair, for the bit-parity check: the band is allowed to miss
/// the alignment entirely, because parity must hold there too.
///
/// Qualities run over every value a `Read` accepts, `Q0` to `Q254`: below
/// Q2 the emission caps `eps` at 3/4, below Q6 a gap-open quality counts as
/// Q6, and at the top a single mismatch costs `10^-25`, which is where the
/// `f32` kernels' range runs out and the `f64` rescue takes over.
pub fn arbitrary_case() -> impl PrintableGenerator<Case> {
    arbitrary_case_inner(0, 254).print_as_debug()
}

#[hegel::composite]
fn arbitrary_case_inner(tc: &TestCase, lowest: u8, highest: u8) -> Case {
    let len = tc.draw_silent(length(1, 89));
    let hap = exactly(tc, any_base_or_n(), len);
    let len = tc.draw_silent(length(1, 69));
    let bases = exactly(tc, any_base_or_n(), len);
    let quals = draw_quals(tc, lowest, highest, bases.len());
    let gaps = draw_quals(tc, lowest, highest, bases.len());
    let strand = tc.draw_silent(any_strand());
    let offset = tc.draw_silent(gs::integers::<i32>().min_value(-20).max_value(39));
    let Ok(read) = Read::new(bases, &quals, &gaps, &gaps, &gaps, strand) else { tc.reject() };
    let betas = exactly(tc, any_probability(), hap.len());
    Case { haplotype: Haplotype::new(hap), read, offset, betas }
}

/// [`arbitrary_case`] with the insertion, deletion and gap-continuation
/// qualities each the same at every base, as a caller with no per-base gap
/// model has: the reads whose cells stay probabilities, so the trust floor is
/// at its lowest and a flush it misses would show.
pub fn steady_case() -> impl PrintableGenerator<Case> {
    steady_case_inner().print_as_debug()
}

#[hegel::composite]
fn steady_case_inner(tc: &TestCase) -> Case {
    let len = tc.draw_silent(length(1, 89));
    let hap = exactly(tc, any_base_or_n(), len);
    let len = tc.draw_silent(length(1, 69));
    let bases = exactly(tc, any_base_or_n(), len);
    let quals = draw_quals(tc, 0, 254, bases.len());
    let [insertion, deletion, gap] = [0; 3].map(|_| {
        BaseQuality::from_byte(tc.draw_silent(gs::integers::<u8>().max_value(254)))
    });
    let strand = tc.draw_silent(any_strand());
    let offset = tc.draw_silent(gs::integers::<i32>().min_value(-20).max_value(39));
    let Ok(read) = Read::uniform(bases, &quals, insertion, deletion, gap, strand) else {
        tc.reject()
    };
    let betas = exactly(tc, any_probability(), hap.len());
    Case { haplotype: Haplotype::new(hap), read, offset, betas }
}

/// The mirror of a case: reverse-complement the haplotype and the read, swap
/// the strand, and reverse the per-site betas.
pub fn mirror(case: &Case) -> Option<Case> {
    let bases: Vec<Base> = case.read.bases().iter().rev().map(|base| base.inverse()).collect();
    Some(Case {
        haplotype: case.haplotype.reverse_complement(),
        read: Read::new(
            bases,
            &reversed(case.read.base_quals()),
            &reversed(case.read.insertion_quals()),
            &reversed(case.read.deletion_quals()),
            &reversed(case.read.gap_quals()),
            mirror_strand(case.read.strand()),
        )
        .ok()?,
        offset: case.offset,
        betas: case.betas.iter().rev().copied().collect(),
    })
}

fn reversed(quals: &[BaseQuality]) -> Vec<BaseQuality> {
    quals.iter().rev().copied().collect()
}

/// GATK's forward recurrence in `f64` over the whole matrix, with the prior
/// supplied per cell instead of by an `Emission`: the way to score a
/// composition the crate does not offer, since `MatchProbability` cannot be
/// overridden. Written out rather than derived from `align_full`, and without
/// its rescaling, so reads of a few hundred bases at most. `None` if the pair
/// cannot be scored at all.
#[allow(
    clippy::indexing_slicing,
    reason = "every index is in 0..=h and every row is allocated with h + 1 entries"
)]
pub fn forward_with(
    haplotype: &Haplotype,
    read: &Read,
    prior: impl Fn(HapSite, Observation) -> f64,
) -> Option<f64> {
    let (h, r) = (haplotype.len(), read.len());
    if h == 0 || r == 0 {
        return None;
    }
    let quality = |track: &[BaseQuality], index: usize| {
        error_probability(track.get(index).copied().unwrap_or(BaseQuality::from_byte(0)))
    };
    #[allow(clippy::cast_precision_loss, reason = "test haplotypes are short")]
    let init = 1.0 / h as f64;
    let mut prev_m = vec![0.0f64; h + 1];
    let mut prev_i = vec![0.0f64; h + 1];
    let mut prev_d = vec![init; h + 1];
    let (mut cur_m, mut cur_i, mut cur_d) =
        (vec![0.0f64; h + 1], vec![0.0f64; h + 1], vec![0.0f64; h + 1]);

    for i in 1..=r {
        let observation = read.observation(i - 1)?;
        let p_ins = quality(read.insertion_quals(), i - 1);
        let p_del = quality(read.deletion_quals(), i - 1);
        let gap = quality(read.gap_quals(), i - 1);
        let m2m = 1.0 - (p_ins + p_del).min(1.0);
        let i2m = 1.0 - gap;
        cur_m[0] = 0.0;
        cur_i[0] = 0.0;
        cur_d[0] = 0.0;
        for j in 1..=h {
            let site = haplotype.site(j - 1)?;
            let p = prior(site, observation);
            cur_m[j] = p * (prev_m[j - 1] * m2m + prev_i[j - 1] * i2m + prev_d[j - 1] * i2m);
            cur_i[j] = prev_m[j] * p_ins + prev_i[j] * gap;
            cur_d[j] = cur_m[j - 1] * p_del + cur_d[j - 1] * gap;
        }
        core::mem::swap(&mut prev_m, &mut cur_m);
        core::mem::swap(&mut prev_i, &mut cur_i);
        core::mem::swap(&mut prev_d, &mut cur_d);
    }
    let total: f64 = prev_m.iter().zip(prev_i.iter()).skip(1).map(|(m, i)| m + i).sum();
    (total > 0.0).then(|| total.log10())
}

/// Every SIMD level this CPU can run, the scalar `Fallback` included (the
/// dev-dependency compiles it in on every target), named for failure
/// messages. The kernels are tested at each of them, not only at the one
/// `Level::new()` picks: on x86-64 each is a separate instantiation.
pub fn levels() -> Vec<(&'static str, fearless_simd::Level)> {
    use fearless_simd::Level;
    let best = Level::new();
    #[allow(unused_mut, reason = "only x86 has levels below the best one")]
    let mut levels = vec![("fallback", Level::fallback()), ("best", best)];
    #[cfg(target_arch = "x86_64")]
    {
        levels.extend(best.as_sse2().map(|token| ("sse2", Level::Sse2(token))));
        levels.extend(best.as_sse4_2().map(|token| ("sse4.2", Level::Sse4_2(token))));
        levels.extend(best.as_avx2().map(|token| ("avx2", Level::Avx2(token))));
    }
    levels
}
