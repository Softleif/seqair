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
use proptest::prelude::*;

/// One of the four nucleotides.
pub fn any_base() -> impl Strategy<Value = Base> {
    prop_oneof![Just(Base::A), Just(Base::C), Just(Base::G), Just(Base::T)]
}

/// A nucleotide, or every so often an `N`: the kernels have a code for it and
/// the reference has a branch, and the two have to agree through the whole
/// DP, not only per cell.
pub fn any_base_or_n() -> impl Strategy<Value = Base> {
    prop_oneof![
        8 => any_base(),
        1 => Just(Base::Unknown),
    ]
}

pub fn any_strand() -> impl Strategy<Value = Strand> {
    prop_oneof![Just(Strand::OT), Just(Strand::OB)]
}

pub fn any_probability() -> impl Strategy<Value = Probability> {
    (0.0f64..=1.0).prop_map(|value| Probability::new(value).unwrap_or(Probability::ZERO))
}

pub fn any_conversion() -> impl Strategy<Value = ConversionModel> {
    (any_probability(), any_probability()).prop_map(|(c, f)| ConversionModel::new(c, f))
}

/// A conversion model a real library could have: efficient, and rarely wrong
/// about an unmethylated base.
///
/// `any_conversion` includes `c = 0, f = 1` -- "every cytosine converts, none
/// of the methylated ones" -- under which a read is *unlikely* against the very
/// haplotype it was cut from, and alignments elsewhere become competitive. Any
/// property whose premise is "the optimal path is near the seed offset" needs
/// this one instead.
pub fn plausible_conversion() -> impl Strategy<Value = ConversionModel> {
    (0.5f64..=1.0, 0.0f64..=0.05).prop_map(|(c, f)| {
        ConversionModel::new(
            Probability::new(c).unwrap_or(Probability::ONE),
            Probability::new(f).unwrap_or(Probability::ZERO),
        )
    })
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

/// One point substitution or short indel applied to a read drawn from the
/// haplotype.
#[derive(Debug, Clone, Copy)]
pub enum Edit {
    Substitute { at: usize, to: Base },
    Insert { at: usize, base: Base },
    Delete { at: usize },
}

pub fn any_edit() -> impl Strategy<Value = Edit> {
    prop_oneof![
        (0usize..4096, any_base_or_n()).prop_map(|(at, to)| Edit::Substitute { at, to }),
        (0usize..4096, any_base_or_n()).prop_map(|(at, base)| Edit::Insert { at, base }),
        (0usize..4096).prop_map(|at| Edit::Delete { at }),
    ]
}

/// A read cut out of a haplotype and then perturbed by at most `max_edits`
/// edits, so the optimal alignment stays within a few columns of the offset the
/// case reports and a default band contains it by construction.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    reason = "test scaffolding over sequences of a few hundred bases"
)]
pub fn derived_case(max_edits: usize) -> impl Strategy<Value = Case> {
    (
        proptest::collection::vec(any_base_or_n(), 80..160),
        0usize..25,
        40usize..90,
        proptest::collection::vec(any_edit(), 0..=max_edits),
        proptest::collection::vec(2u8..=45, 200),
        any_strand(),
        proptest::collection::vec(any_probability(), 200),
    )
        .prop_filter_map(
            "the read must survive its edits",
            |(hap, start, want, edits, quals, strand, betas)| {
                let start = start.min(hap.len().saturating_sub(10));
                let end = (start + want).min(hap.len());
                let mut bases: Vec<Base> = hap.get(start..end)?.to_vec();
                for edit in edits {
                    if bases.len() < 8 {
                        break;
                    }
                    match edit {
                        Edit::Substitute { at, to } => {
                            let at = at % bases.len();
                            *bases.get_mut(at)? = to;
                        }
                        Edit::Insert { at, base } => bases.insert(at % bases.len(), base),
                        Edit::Delete { at } => {
                            let _removed = bases.remove(at % bases.len());
                        }
                    }
                }
                if bases.len() < 8 {
                    return None;
                }
                let quals: Vec<BaseQuality> = bases
                    .iter()
                    .enumerate()
                    .map(|(index, _)| {
                        BaseQuality::from_byte(*quals.get(index % quals.len()).unwrap_or(&30))
                    })
                    .collect();
                let read = Read::uniform(
                    bases,
                    &quals,
                    BaseQuality::from_byte(45),
                    BaseQuality::from_byte(45),
                    BaseQuality::from_byte(10),
                    strand,
                )
                .ok()?;
                let betas = (0..hap.len())
                    .map(|index| *betas.get(index % betas.len()).unwrap_or(&Probability::ZERO))
                    .collect();
                Some(Case { haplotype: Haplotype::new(hap), read, offset: start as i32, betas })
            },
        )
}

/// An unconstrained pair, for the bit-parity check: the band is allowed to miss
/// the alignment entirely, because parity must hold there too.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    reason = "test scaffolding over sequences of a few hundred bases"
)]
pub fn arbitrary_case() -> impl Strategy<Value = Case> {
    (
        proptest::collection::vec(any_base_or_n(), 1..90),
        proptest::collection::vec(any_base_or_n(), 1..70),
        proptest::collection::vec(2u8..=45, 70),
        proptest::collection::vec(2u8..=45, 70),
        any_strand(),
        -20i32..40,
        proptest::collection::vec(any_probability(), 90),
    )
        .prop_filter_map(
            "the read must be constructible",
            |(hap, bases, quals, gaps, strand, offset, betas)| {
                let quals: Vec<BaseQuality> = bases
                    .iter()
                    .enumerate()
                    .map(|(index, _)| {
                        BaseQuality::from_byte(*quals.get(index % quals.len()).unwrap_or(&30))
                    })
                    .collect();
                let gaps: Vec<BaseQuality> = bases
                    .iter()
                    .enumerate()
                    .map(|(index, _)| {
                        BaseQuality::from_byte(*gaps.get(index % gaps.len()).unwrap_or(&10))
                    })
                    .collect();
                let read = Read::new(bases, &quals, &gaps, &gaps, &gaps, strand).ok()?;
                let betas = (0..hap.len())
                    .map(|index| *betas.get(index % betas.len()).unwrap_or(&Probability::ZERO))
                    .collect();
                Some(Case { haplotype: Haplotype::new(hap), read, offset, betas })
            },
        )
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
