#[allow(dead_code, reason = "each integration test uses a different part of this")]
mod support;

use compair::{
    Base, BaseQuality, Betas, ConversionModel, CpgRole, Emission, HapSite, Haplotype, MAX_EPSILON,
    MatchProbability, Observation, Probability, QPos, Read, StandardEmission, Strand, TapsEmission,
    error_probability,
};
use hegel::TestCase;
use hegel::generators::{self as gs, Generator};
use support::{
    any_base, any_conversion, any_probability, any_strand, derived_case, forward_with, mirror,
};

fn observe(base: Base, qual: BaseQuality, strand: Strand) -> Observation {
    Observation { index: QPos::new(0), base, error_probability: error_probability(qual), strand }
}

fn perfect() -> ConversionModel {
    ConversionModel::new(Probability::ONE, Probability::ZERO)
}

fn taps(betas: &[Probability], conversion: ConversionModel) -> TapsEmission<'_> {
    TapsEmission::new(conversion, Betas::PerSite(betas))
}

/// `CpG` context comes from the haplotype's own sequence, so an allele that
/// creates one turns the conversion rows on and an allele that destroys one
/// turns them off. That is the whole of de-novo `CpG` handling.
#[test]
fn cpg_context_is_read_from_the_haplotype() {
    let reference = Haplotype::from_ascii(b"AATGAA");
    let created = Haplotype::from_ascii(b"AACGAA");
    let destroyed = Haplotype::from_ascii(b"AACAAA");

    assert_eq!(
        reference.site(2).map(|site| site.cpg),
        Some(CpgRole::None),
        "a T is not a CpG cytosine"
    );
    assert_eq!(
        created.site(2).map(|site| site.cpg),
        Some(CpgRole::TopC),
        "a T>C allele creates the CpG"
    );
    assert_eq!(created.site(3).map(|site| site.cpg), Some(CpgRole::BottomG), "and its partner G");
    assert_eq!(
        destroyed.site(2).map(|site| site.cpg),
        Some(CpgRole::None),
        "a G>A allele on the partner destroys it"
    );
}

#[test]
fn beta_zero_makes_a_converted_base_a_mismatch() {
    let haplotype = Haplotype::from_ascii(b"AACGAA");
    let site = haplotype.site(2).expect("in range");
    let betas = vec![Probability::ZERO; 6];
    let emission = taps(&betas, perfect());
    let read = observe(Base::T, BaseQuality::from_byte(30), Strand::OT);

    assert!(
        (emission.match_probability(site, read)
            - StandardEmission::default().match_probability(site, read))
        .abs()
            < f64::EPSILON,
        "with f = 0 an unmethylated CpG scores a T exactly as a mismatch"
    );
}

#[test]
fn beta_one_makes_a_converted_base_a_near_match() {
    let haplotype = Haplotype::from_ascii(b"AACGAA");
    let site = haplotype.site(2).expect("in range");
    let betas = vec![Probability::ONE; 6];
    let emission = taps(&betas, perfect());
    let qual = BaseQuality::from_byte(30);
    let read = observe(Base::T, qual, Strand::OT);

    assert!(
        emission.match_probability(site, read) >= 1.0 - error_probability(qual),
        "with c = 1 a fully methylated CpG scores a T at least as well as a match"
    );
    assert!(
        emission.match_probability(site, observe(Base::C, qual, Strand::OT))
            < emission.match_probability(site, read),
        "and the unconverted C is now the unlikely observation"
    );
}

/// The leniency is granted only where the model says conversion happens: the
/// same `C` read on the bottom strand is a plain reference base.
#[test]
fn the_other_strand_sees_no_conversion() {
    let haplotype = Haplotype::from_ascii(b"AACGAA");
    let site = haplotype.site(2).expect("in range");
    let betas = vec![Probability::ONE; 6];
    let emission = taps(&betas, perfect());
    let read = observe(Base::T, BaseQuality::from_byte(30), Strand::OB);

    assert_eq!(
        emission.match_probability(site, read),
        StandardEmission::default().match_probability(site, read)
    );
}

/// Strand symmetry: reverse-complement the
/// haplotype and the read, swap the strand, reverse the per-site betas, and
/// every emission is the one it mirrors — bit for bit, not to a tolerance.
#[hegel::test]
fn emission_mirrors_between_strands(tc: TestCase) {
    let case = tc.draw(derived_case(3));
    let conversion = tc.draw(any_conversion());
    let mirrored = mirror(&case).expect("mirror must build");
    // No assumption about a CpG any more: since the conversion rows apply
    // at every cytosine, `C` on OT mirrors `G` on OB whether or not either
    // is in a CpG, and the CpG only decides which rate is used.
    let forward = taps(&case.betas, conversion);
    let backward = taps(&mirrored.betas, conversion);
    let (h, r) = (case.haplotype.len(), case.read.len());

    for j in 0..h {
        for i in 0..r {
            let site = case.haplotype.site(j).expect("site");
            let obs = case.read.observation(i).expect("obs");
            let mirrored_site = mirrored.haplotype.site(h - 1 - j).expect("site");
            let mirrored_obs = mirrored.read.observation(r - 1 - i).expect("obs");
            assert_eq!(
                forward.match_probability(site, obs).to_bits(),
                backward.match_probability(mirrored_site, mirrored_obs).to_bits(),
                "hap {} read {}: {:?} vs {:?}",
                j,
                i,
                site,
                mirrored_site
            );
        }
    }
}

/// Where the chemistry cannot act -- the haplotype base is not a cytosine
/// on the strand this read reports -- `TapsEmission` is `StandardEmission`
/// for every conversion model, not only the perfect one.
#[hegel::test]
fn taps_falls_through_where_the_chemistry_cannot_act(tc: TestCase) {
    let strand = tc.draw(any_strand());
    // A lone base has no partner, so it is never in CpG context; what
    // matters here is that it is also not a cytosine this read could see
    // converted.
    let haplotype_base = tc.draw(
        gs::sampled_from(match strand {
            Strand::OB => &[Base::A, Base::C, Base::T],
            _ => &[Base::A, Base::G, Base::T],
        })
        .print_as_debug(),
    );
    let read_base = tc.draw(any_base());
    let qual = tc.draw(gs::integers::<u8>().min_value(2).max_value(45));
    let beta = tc.draw(any_probability());
    let conversion = tc.draw(any_conversion());
    let haplotype = Haplotype::new(vec![haplotype_base]);
    let site = haplotype.site(0).expect("site");
    assert_eq!(site.cpg, CpgRole::None);
    let betas = vec![beta];
    let observation = observe(read_base, BaseQuality::from_byte(qual), strand);
    assert_eq!(
        taps(&betas, conversion).match_probability(site, observation).to_bits(),
        StandardEmission::default().match_probability(site, observation).to_bits()
    );
}

/// And where it can, the rows apply at the false-conversion rate whatever
/// the site's `beta`, because a cytosine outside a `CpG` is unmethylated.
#[hegel::test]
fn a_non_cpg_cytosine_reads_t_at_the_false_conversion_rate(tc: TestCase) {
    let read_base = tc.draw(any_base());
    let qual = tc.draw(gs::integers::<u8>().min_value(2).max_value(45));
    let strand = tc.draw(any_strand());
    let beta = tc.draw(any_probability());
    let conversion = tc.draw(any_conversion());
    let (haplotype_base, converted, unconverted) = match strand {
        Strand::OT => (Base::C, Base::T, Base::C),
        Strand::OB => (Base::G, Base::A, Base::G),
        Strand::Unknown => tc.reject(),
    };
    let haplotype = Haplotype::new(vec![haplotype_base]);
    let site = haplotype.site(0).expect("site");
    assert_eq!(site.cpg, CpgRole::None);
    let betas = vec![beta];
    let emission = taps(&betas, conversion);
    let observation = observe(read_base, BaseQuality::from_byte(qual), strand);
    let eps = error_probability(BaseQuality::from_byte(qual));
    let f = *conversion.false_conversion();

    let want = if read_base == converted {
        f * (1.0 - eps) + (1.0 - f) * eps / 3.0
    } else if read_base == unconverted {
        (1.0 - f) * (1.0 - eps) + f * eps / 3.0
    } else {
        eps / 3.0
    };
    assert!(
        (emission.match_probability(site, observation) - want).abs() < 1e-15,
        "{:?} over a non-CpG {:?} on {:?}: {} wanted {}",
        read_base,
        haplotype_base,
        strand,
        emission.match_probability(site, observation),
        want
    );
    // `beta` is ignored outside a CpG, so the answer does not move with it.
    let other = vec![Probability::ONE];
    assert_eq!(
        emission.match_probability(site, observation).to_bits(),
        taps(&other, conversion).match_probability(site, observation).to_bits()
    );
}

/// Every emission is a probability, whatever the conversion model.
#[hegel::test]
fn emissions_stay_in_range(tc: TestCase) {
    let case = tc.draw(derived_case(3));
    let conversion = tc.draw(any_conversion());
    let emission = taps(&case.betas, conversion);
    for j in 0..case.haplotype.len() {
        for i in 0..case.read.len() {
            let site = case.haplotype.site(j).expect("site");
            let obs = case.read.observation(i).expect("obs");
            let value = emission.match_probability(site, obs);
            assert!((0.0..=1.0).contains(&value), "emission {value} out of range");
        }
    }
}

/// The strand rule and the converted base, written out once more here so the
/// oracles below do not borrow them from the crate: TAPS converts a cytosine
/// on the strand the read reports, and a `G` in reference orientation is the
/// bottom strand's cytosine.
fn conversion_of(site: HapSite, strand: Strand) -> Option<(Base, bool)> {
    match (site.base, strand) {
        (Base::C, Strand::OT) => Some((Base::T, site.cpg == CpgRole::TopC)),
        (Base::G, Strand::OB) => Some((Base::A, site.cpg == CpgRole::BottomG)),
        _ => None,
    }
}

/// bsgenova's observation model (Feng & Gao 2024, "Bayesian probabilistic
/// model of bsgenova"), the oracle for how conversion and sequencing error
/// compose. A base passes through three stages -- an error before conversion,
/// the conversion, an error after it -- and the probability of the sequenced
/// base sums over the two latent bases in between:
///
/// ```text
/// P(s | g) = sum over x1, x2 of  P(s | x2) P(x2 | x1) P(x1 | g)
/// ```
///
/// Written here for TAPS. The pre-conversion stage is the identity: bsgenova's
/// `p1` models damage during sample preparation, which a pair-HMM leaves to
/// the haplotype. Conversion turns the strand's cytosine into its converted
/// base with probability `rate` and touches nothing else; bsgenova has it
/// deterministic given the methylation state, so `rate` is its `pm` folded
/// with the chemistry's two efficiencies. The post-conversion stage is the
/// sequencer's uniform error, `eps / 3` to each of the other three bases.
fn bsgenova(site: HapSite, strand: Strand, rate: f64, observed: Base, eps: f64) -> f64 {
    let mut after_conversion = [0.0f64; 4];
    let index = |base: Base| base.known_index().expect("known bases only");
    match conversion_of(site, strand) {
        Some((converted, _)) => {
            after_conversion[index(converted)] += rate;
            after_conversion[index(site.base)] += 1.0 - rate;
        }
        None => after_conversion[index(site.base)] = 1.0,
    }
    Base::KNOWN
        .iter()
        .zip(after_conversion)
        .map(|(latent, mass)| mass * if *latent == observed { 1.0 - eps } else { eps / 3.0 })
        .sum()
}

/// The chance a site at level `beta` reads as converted before sequencing,
/// computed from the raw efficiencies rather than through `ConversionModel`.
fn taps_rate(site: HapSite, strand: Strand, beta: f64, c: f64, f: f64) -> f64 {
    match conversion_of(site, strand) {
        Some((_, true)) => beta * c + (1.0 - beta) * f,
        Some((_, false)) => f,
        None => 0.0,
    }
}

/// A haplotype with every kind of site: a `CpG` (`C` at 1, `G` at 2), plain
/// bases, a non-`CpG` `C` at 5 and a non-`CpG` `G` at 8.
fn every_kind_of_site() -> Haplotype {
    Haplotype::from_ascii(b"ACGTACATG")
}

/// The `CpG` base the strand can see converted on [`every_kind_of_site`].
fn cpg_site_for(strand: Strand) -> Option<HapSite> {
    every_kind_of_site().site(if strand == Strand::OB { 2 } else { 1 })
}

/// Bis-SNP's row (Liu et al. 2012, eq. 5) for the strand that sees
/// conversion, in its own parameters: `beta` methylated, `alpha` the
/// under-conversion rate (an unmethylated cytosine bisulfite leaves alone),
/// `gamma` the over-conversion rate (a methylated one it converts anyway).
///
/// ```text
/// P(t | c) = (1 - eps) [beta gamma + (1 - beta)(1 - alpha)] + eps / 3
/// P(c | c) = (1 - eps) [beta (1 - gamma) + (1 - beta) alpha] + eps / 3
/// otherwise  eps / 3
/// ```
///
/// Note the bare `+ eps / 3`: the row sums to `1 + eps / 3`. That is the
/// composition the joint model's first draft copied and this crate no longer
/// uses; see `bis_snp_agrees_on_the_mixture_and_not_on_the_error_term`.
struct BisSnp {
    beta: f64,
    alpha: f64,
    gamma: f64,
}

impl BisSnp {
    fn row(&self, observed: Base, converted: Base, original: Base, eps: f64) -> f64 {
        let latent_t = self.beta * self.gamma + (1.0 - self.beta) * (1.0 - self.alpha);
        let latent_c = self.beta * (1.0 - self.gamma) + (1.0 - self.beta) * self.alpha;
        if observed == converted {
            (1.0 - eps) * latent_t + eps / 3.0
        } else if observed == original {
            (1.0 - eps) * latent_c + eps / 3.0
        } else {
            eps / 3.0
        }
    }
}

/// Every row is a probability distribution over the four bases -- at every
/// kind of site, on both strands, at every quality including `eps = 1`,
/// under every conversion model. The previous composition summed to
/// `1 + eps / 3` on the converting sites, and a row that sums to more than
/// one is free likelihood for every haplotype that has that site.
#[hegel::test]
fn every_row_sums_to_one(tc: TestCase) {
    let qual = tc.draw(gs::integers::<u8>().max_value(93));
    let strand = tc.draw(any_strand());
    let beta = tc.draw(any_probability());
    let conversion = tc.draw(any_conversion());
    let haplotype = every_kind_of_site();
    let emission = TapsEmission::new(conversion, Betas::Uniform(beta));
    for index in 0..haplotype.len() {
        let site = haplotype.site(index).expect("site");
        let row: f64 = Base::KNOWN
            .iter()
            .map(|base| {
                emission
                    .match_probability(site, observe(*base, BaseQuality::from_byte(qual), strand))
            })
            .sum();
        assert!(
            (row - 1.0).abs() < 1e-12,
            "site {} on {:?} at q{}: row sums to {}",
            index,
            strand,
            qual,
            row
        );
    }
}

/// A quality below Q2 claims `eps > 3/4`: that the base is *likelier* to be
/// any one of the other three than the one the sequencer called, which no
/// base call means. At Q0 it claims `eps = 1`, so a read base that matches
/// its haplotype scores exactly zero. Both are read as the most a quality can
/// say against its call, `eps = 3/4`: every base equally likely, whatever the
/// site, the strand, the chemistry or the artifact floor.
#[hegel::test]
fn a_quality_below_q2_carries_no_information(tc: TestCase) {
    let qual = tc.draw(gs::integers::<u8>().max_value(1));
    let strand = tc.draw(any_strand());
    let beta = tc.draw(any_probability());
    let conversion = tc.draw(any_conversion());
    let floor = tc.draw(any_probability());
    let haplotype = every_kind_of_site();
    let taps = TapsEmission::new(conversion, Betas::Uniform(beta)).with_artifact_floor(floor);
    let standard = StandardEmission::default().with_artifact_floor(floor);
    for index in 0..haplotype.len() {
        let site = haplotype.site(index).expect("site");
        for base in Base::KNOWN.into_iter().chain([Base::Unknown]) {
            let observation = observe(base, BaseQuality::from_byte(qual), strand);
            for (name, p) in [
                ("standard", standard.match_probability(site, observation)),
                ("taps", taps.match_probability(site, observation)),
            ] {
                assert!(
                    (p - 0.25).abs() <= f64::EPSILON,
                    "{name}: site {index} on {strand:?} at Q{qual} reads {base:?} at {p}"
                );
            }
        }
    }
}

/// So is an artifact floor above 3/4: the floor raises `eps` to at most the
/// point where the base carries no information, and never past it.
#[test]
fn an_artifact_floor_above_three_quarters_stops_at_three_quarters() {
    let haplotype = Haplotype::from_ascii(b"AACGAA");
    let emission = StandardEmission::default().with_artifact_floor(Probability::ONE);
    let site = haplotype.site(0).expect("A");
    for base in Base::KNOWN {
        let p =
            emission.match_probability(site, observe(base, BaseQuality::from_byte(40), Strand::OT));
        assert_eq!(p, 0.25, "{base:?}");
    }
}

/// The emission is bsgenova's three-stage marginalisation with the
/// pre-conversion stage removed: convert, then sequence, summing over the
/// base the site carried in between.
#[hegel::test]
fn the_emission_is_bsgenova_without_the_pre_conversion_stage(tc: TestCase) {
    let qual = tc.draw(gs::integers::<u8>().max_value(93));
    let strand = tc.draw(any_strand());
    let beta = tc.draw(any_probability());
    let conversion = tc.draw(any_conversion());
    let observed = tc.draw(any_base());
    let haplotype = every_kind_of_site();
    let emission = TapsEmission::new(conversion, Betas::Uniform(beta));
    let (c, f) = (*conversion.efficiency(), *conversion.false_conversion());
    let eps = error_probability(BaseQuality::from_byte(qual)).min(MAX_EPSILON);
    for index in 0..haplotype.len() {
        let site = haplotype.site(index).expect("site");
        let rate = taps_rate(site, strand, *beta, c, f);
        let want = bsgenova(site, strand, rate, observed, eps);
        let got = emission
            .match_probability(site, observe(observed, BaseQuality::from_byte(qual), strand));
        assert!(
            (got - want).abs() < 1e-15,
            "site {} ({:?}) read {:?} on {:?} at q{}: {} but bsgenova says {}",
            index,
            site.base,
            observed,
            strand,
            qual,
            got,
            want
        );
    }
}

/// Bis-SNP is the oracle for the *mixture* and a counter-example for the
/// *composition*. Its bracketed term is this crate's `rate` once the
/// chemistries are mapped onto each other -- bisulfite converts the
/// unmethylated cytosine and TAPS the methylated one, so TAPS's
/// efficiency is bisulfite's over-conversion (`c = gamma`) and its false
/// conversion is what bisulfite does to an unmethylated base
/// (`f = 1 - alpha`). Its row then differs from ours by exactly the
/// unweighted error term's excess, `eps / 3` times the latent weight, on
/// the two conversion outcomes and nowhere else: Bis-SNP's row sums to
/// `1 + eps / 3`.
#[hegel::test]
fn bis_snp_agrees_on_the_mixture_and_not_on_the_error_term(tc: TestCase) {
    let qual = tc.draw(gs::integers::<u8>().max_value(93));
    let strand = tc.draw(any_strand());
    let beta = tc.draw(any_probability());
    let conversion = tc.draw(any_conversion());
    let observed = tc.draw(any_base());
    let emission = TapsEmission::new(conversion, Betas::Uniform(beta));
    let bis_snp = BisSnp {
        beta: *beta,
        alpha: 1.0 - *conversion.false_conversion(),
        gamma: *conversion.efficiency(),
    };
    let eps = error_probability(BaseQuality::from_byte(qual)).min(MAX_EPSILON);
    let cpg = cpg_site_for(strand).expect("site");
    let (converted, _) = conversion_of(cpg, strand).expect("a CpG converts");

    // The mixture: Bis-SNP's P(latent t) is our rate.
    let weights = emission.site_weights(cpg, strand);
    let latent_t = bis_snp.beta * bis_snp.gamma + (1.0 - bis_snp.beta) * (1.0 - bis_snp.alpha);
    assert!(
        (weights.rate - latent_t).abs() < 1e-15,
        "rate {} vs Bis-SNP {}",
        weights.rate,
        latent_t
    );

    // The composition: differs by the bare error term, weighted by us and
    // not by Bis-SNP.
    let ours =
        emission.match_probability(cpg, observe(observed, BaseQuality::from_byte(qual), strand));
    let theirs = bis_snp.row(observed, converted, cpg.base, eps);
    let excess = weights.latent_weight(observed) * eps / 3.0;
    assert!(
        (theirs - ours - excess).abs() < 1e-15,
        "read {:?}: Bis-SNP {}, ours {}, excess should be {}",
        observed,
        theirs,
        ours,
        excess
    );
    let their_row: f64 =
        Base::KNOWN.iter().map(|base| bis_snp.row(*base, converted, cpg.base, eps)).sum();
    assert!((their_row - (1.0 + eps / 3.0)).abs() < 1e-12, "Bis-SNP's row sums to {}", their_row);
}

/// A conversion model with nothing to convert -- no false conversion, and
/// either no methylation or no efficiency -- is `StandardEmission` bit for
/// bit, at every cell and through the whole dynamic program. The previous
/// composition failed this: its `+ eps / 3` rode along on every cytosine
/// the strand could see, converting or not.
#[hegel::test]
fn a_model_with_nothing_to_convert_is_the_standard_emission(tc: TestCase) {
    let case = tc.draw(derived_case(3));
    let efficiency = tc.draw(any_probability());
    let beta = tc.draw(any_probability());
    let kill_beta = tc.draw(gs::booleans());
    let (efficiency, beta) =
        if kill_beta { (efficiency, Probability::ZERO) } else { (Probability::ZERO, beta) };
    let conversion = ConversionModel::new(efficiency, Probability::ZERO);
    let taps = TapsEmission::new(conversion, Betas::Uniform(beta));
    let standard = StandardEmission::default();
    for j in 0..case.haplotype.len() {
        for i in 0..case.read.len() {
            let site = case.haplotype.site(j).expect("site");
            let obs = case.read.observation(i).expect("obs");
            assert_eq!(
                taps.match_probability(site, obs).to_bits(),
                standard.match_probability(site, obs).to_bits(),
                "hap {} read {}",
                j,
                i
            );
        }
    }
    assert_eq!(
        compair::align_full(&case.haplotype, &case.read, &taps).get().to_bits(),
        compair::align_full(&case.haplotype, &case.read, &standard).get().to_bits()
    );
}

/// At the two ends of `beta` the site is a two-state chemistry: fully
/// methylated converts at the efficiency, unmethylated at the false
/// conversion rate, and the rows are the marginalisation at that one rate.
#[hegel::test]
fn the_beta_limits_are_the_two_state_rows(tc: TestCase) {
    let qual = tc.draw(gs::integers::<u8>().max_value(93));
    let strand = tc.draw(any_strand());
    let conversion = tc.draw(any_conversion());
    let observed = tc.draw(any_base());
    let eps = error_probability(BaseQuality::from_byte(qual)).min(MAX_EPSILON);
    let cpg = cpg_site_for(strand).expect("site");
    for (beta, rate) in [
        (Probability::ONE, *conversion.efficiency()),
        (Probability::ZERO, *conversion.false_conversion()),
    ] {
        let emission = TapsEmission::new(conversion, Betas::Uniform(beta));
        assert_eq!(emission.site_weights(cpg, strand).rate.to_bits(), rate.to_bits());
        let got = emission
            .match_probability(cpg, observe(observed, BaseQuality::from_byte(qual), strand));
        let want = bsgenova(cpg, strand, rate, observed, eps);
        assert!((got - want).abs() < 1e-15, "beta {:?}: {} vs {}", beta, got, want);
    }
}

/// The composition this crate used before, and Bis-SNP still does: the error
/// term added bare rather than weighted by the latent base. A function, not an
/// `Emission`, because the crate composes the halves itself and offers no
/// other composition -- `MatchProbability` cannot be overridden -- so scoring
/// it takes the test suite's own forward recurrence.
fn flat_composition(emission: &TapsEmission<'_>, site: HapSite, observation: Observation) -> f64 {
    let weights = emission.site_weights(site, observation.strand);
    let eps = emission.epsilon(observation);
    if observation.base == Base::Unknown || site.base == Base::Unknown {
        return 1.0 - eps;
    }
    match weights.converted {
        Some(converted) if observation.base == converted => (1.0 - eps) * weights.rate + eps / 3.0,
        Some(_) if observation.base == weights.base => {
            (1.0 - eps) * (1.0 - weights.rate) + eps / 3.0
        }
        _ if observation.base == weights.base => 1.0 - eps,
        _ => eps / 3.0,
    }
}

/// What the correction is worth, on the comparison a caller actually makes:
/// one read against two haplotypes that differ by a `CpG`-destroying `C>T`.
///
/// The flat composition inflates every converting column by up to `eps / 3`,
/// so it favours the haplotype *with* the cytosine -- the reference over the
/// `C>T` allele -- whatever the read says, and by more as quality falls. The
/// marginalised composition has no such term: with nothing to convert it is
/// `StandardEmission` exactly (pinned above), and under a real chemistry the
/// two haplotypes are compared on the conversion evidence alone.
#[test]
fn the_flat_composition_favours_the_cytosine_haplotype_by_more_at_low_quality() {
    let with_cpg = Haplotype::from_ascii(b"TTAGCATCGGATCCGATTACAGGCATTACGGATCCAGT");
    let destroyed = Haplotype::from_ascii(b"TTAGCATTGGATCCGATTACAGGCATTACGGATCCAGT");
    assert_eq!(with_cpg.site(7).map(|site| site.cpg), Some(CpgRole::TopC));
    assert_eq!(destroyed.site(7).map(|site| site.cpg), Some(CpgRole::None));
    // A read cut from the destroyed haplotype: it carries the T, so the C
    // haplotype can only explain it as a conversion.
    let read_bases = destroyed.bases().get(2..36).expect("in range").to_vec();
    let taps = TapsEmission::new(ConversionModel::taps_default(), Betas::Uniform(Probability::ONE));

    let mut gap_by_quality = Vec::new();
    for qual in [30u8, 10] {
        let read = Read::uniform(
            read_bases.clone(),
            &vec![BaseQuality::from_byte(qual); read_bases.len()],
            BaseQuality::from_byte(45),
            BaseQuality::from_byte(45),
            BaseQuality::from_byte(10),
            Strand::OT,
        )
        .expect("valid");
        let score = |haplotype: &Haplotype, prior: &dyn Fn(HapSite, Observation) -> f64| {
            forward_with(haplotype, &read, prior).expect("the pair is scorable")
        };
        let marginalised = |site, obs| taps.match_probability(site, obs);
        let flat = |site, obs| flat_composition(&taps, site, obs);

        // The suite's recurrence is the crate's, to rounding.
        for haplotype in [&with_cpg, &destroyed] {
            let crate_score = compair::align_full(haplotype, &read, &taps).get();
            assert!((score(haplotype, &marginalised) - crate_score).abs() < 1e-9);
        }

        let margin_marginalised =
            score(&with_cpg, &marginalised) - score(&destroyed, &marginalised);
        let margin_flat = score(&with_cpg, &flat) - score(&destroyed, &flat);
        assert!(
            margin_flat > margin_marginalised,
            "q{qual}: the flat composition should favour the C haplotype more: \
             {margin_flat} vs {margin_marginalised}"
        );
        gap_by_quality.push(margin_flat - margin_marginalised);
    }
    let (Some(at_q30), Some(at_q10)) = (gap_by_quality.first(), gap_by_quality.get(1)) else {
        panic!("two qualities");
    };
    // Every one of the read's cytosines on this strand carries the excess, not
    // only the CpG, so the bias is a property of the read's composition and
    // grows roughly a hundredfold from Q30 to Q10 as `eps` does.
    assert!(*at_q10 > 10.0 * at_q30, "Q10 gap {at_q10} should dwarf Q30 gap {at_q30}");
    assert!(*at_q10 > 0.01, "at Q10 the bias is not bookkeeping: {at_q10} log10");
}

/// `N` matches everything at full probability on either side, and the `CpG`
/// rows do not claim it: an `N` read over a methylated `CpG` is neither the
/// converted nor the unconverted base, so it falls through to the match row.
#[test]
fn n_matches_everything_and_never_takes_a_conversion_row() {
    let with_n = Haplotype::from_ascii(b"ACNGT");
    let n_site = with_n.site(2).expect("in range");
    let qual = BaseQuality::from_byte(30);
    let matched = 1.0 - error_probability(qual);
    for base in [Base::A, Base::C, Base::G, Base::T, Base::Unknown] {
        let got =
            StandardEmission::default().match_probability(n_site, observe(base, qual, Strand::OT));
        assert!((got - matched).abs() < 1e-12, "hap N against read {base:?} gave {got}");
    }
    let cpg = Haplotype::from_ascii(b"ACGT");
    let cpg_c = cpg.site(1).expect("in range");
    let got = StandardEmission::default()
        .match_probability(cpg_c, observe(Base::Unknown, qual, Strand::OT));
    assert!((got - matched).abs() < 1e-12, "hap C against read N gave {got}");
    let betas = vec![Probability::ONE; cpg.len()];
    let under_taps = taps(&betas, ConversionModel::taps_default())
        .match_probability(cpg_c, observe(Base::Unknown, qual, Strand::OT));
    assert!((under_taps - matched).abs() < 1e-12, "a CpG C against read N gave {under_taps}");
}

/// **A caller contract worth knowing.** `Betas::PerSite` is indexed by
/// haplotype *position*, and a `CpG` occupies two of them. The `C` reads
/// `betas[j]` and its partner `G` reads `betas[j + 1]`, so a caller that fills
/// only one of the pair makes the two strands of one `CpG` disagree about its
/// methylation. `Betas::Uniform` has no such trap.
#[test]
fn the_two_bases_of_one_cpg_read_their_own_betas() {
    let haplotype = Haplotype::from_ascii(b"AACGAA");
    let (cytosine, guanine) = (haplotype.site(2).expect("C"), haplotype.site(3).expect("G"));
    let mut betas = vec![Probability::ZERO; haplotype.len()];
    let Some(slot) = betas.get_mut(2) else { panic!("six entries") };
    *slot = Probability::ONE;

    let emission = taps(&betas, perfect());
    let qual = BaseQuality::from_byte(30);
    let converted_top = emission.match_probability(cytosine, observe(Base::T, qual, Strand::OT));
    let converted_bottom = emission.match_probability(guanine, observe(Base::A, qual, Strand::OB));
    assert!(converted_top > 0.9, "the C reads betas[2] = 1 and sees a T as a near-match");
    assert!(
        converted_bottom < 0.01,
        "its partner G reads betas[3] = 0 and sees an A as a mismatch: {converted_bottom}"
    );

    // With `Uniform` the pair cannot disagree.
    let uniform = TapsEmission::new(perfect(), Betas::Uniform(Probability::ONE));
    assert_eq!(
        uniform.match_probability(cytosine, observe(Base::T, qual, Strand::OT)).to_bits(),
        uniform.match_probability(guanine, observe(Base::A, qual, Strand::OB)).to_bits()
    );
}

/// A converted base only ever gets likelier as the site gets more
/// methylated, and its unconverted partner only ever gets unlikelier --
/// whenever the chemistry converts methylated bases more often than
/// unmethylated ones, which is what `efficiency >= false_conversion` says.
#[hegel::test]
fn a_converted_base_is_monotone_in_beta(tc: TestCase) {
    // Every quality: below Q2 `eps` is capped at 3/4 and every row is flat
    // in the latent weight, which is monotone too.
    let unit = || gs::floats::<f64>().min_value(0.0).max_value(1.0);
    let (a, b) = (tc.draw(unit()), tc.draw(unit()));
    let (efficiency, false_conversion) = if a >= b { (a, b) } else { (b, a) };
    let qual = tc.draw(gs::integers::<u8>().max_value(93));
    let levels = [tc.draw(unit()), tc.draw(unit())];
    let conversion = ConversionModel::new(
        Probability::new(efficiency).expect("c"),
        Probability::new(false_conversion).expect("f"),
    );
    let mut levels = levels;
    levels.sort_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));
    let (low, high) = (*levels.first().expect("two"), *levels.last().expect("two"));
    let haplotype = Haplotype::from_ascii(b"AACGAA");
    let site = haplotype.site(2).expect("site");
    let read = observe(Base::T, BaseQuality::from_byte(qual), Strand::OT);
    let unconverted = observe(Base::C, BaseQuality::from_byte(qual), Strand::OT);

    let at = |level: f64| {
        let betas = vec![Probability::new(level).expect("beta"); haplotype.len()];
        let emission = taps(&betas, conversion);
        (emission.match_probability(site, read), emission.match_probability(site, unconverted))
    };
    let (t_low, c_low) = at(low);
    let (t_high, c_high) = at(high);
    // Where `efficiency` and `false_conversion` are equal or nearly so, the
    // rows barely move with beta, and rounding alone can tip them either
    // way. Each row is a few products of probabilities, so rounding is
    // bounded by a few ulps of one, absolutely.
    let slack = 8.0 * f64::EPSILON;
    assert!(t_high >= t_low - slack, "T: {t_low} at beta {low} but {t_high} at beta {high}");
    assert!(c_high <= c_low + slack, "C: {c_low} at beta {low} but {c_high} at beta {high}");
}

/// And the same through the whole dynamic program: a read whose every `CpG`
/// cytosine is converted scores higher as the haplotype's methylation rises.
#[test]
fn the_dp_is_monotone_in_beta() {
    let haplotype = Haplotype::from_ascii(
        b"TCATTGGCTATCCTAACCCGACCCTAGGAGCGGTTGGCGTGTATGCCGTGAATTTTCTCATTTCCGCTAGACATAATCGTTCTGCCTATA",
    );
    let mut bases = haplotype.bases().get(20..60).expect("in range").to_vec();
    for (offset, base) in bases.iter_mut().enumerate() {
        if haplotype.site(20 + offset).map(|site| site.cpg) == Some(CpgRole::TopC) {
            *base = Base::T;
        }
    }
    let read = Read::uniform(
        bases,
        &[BaseQuality::from_byte(30); 40],
        BaseQuality::from_byte(45),
        BaseQuality::from_byte(45),
        BaseQuality::from_byte(10),
        Strand::OT,
    )
    .expect("valid");

    let mut previous = f64::NEG_INFINITY;
    for level in [0.0f64, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0] {
        let betas = vec![Probability::new(level).expect("in [0, 1]"); haplotype.len()];
        let score =
            compair::align_full(&haplotype, &read, &taps(&betas, ConversionModel::taps_default()))
                .get();
        assert!(score >= previous - 1e-12, "beta {level} scored {score} under {previous}");
        previous = score;
    }
}

/// The false-conversion row, at the numbers rastair measured, with the size of the
/// correction spelled out.
///
/// `f` was measured on the unmethylated pUC19 spike-in's cytosines, so it is
/// the chemistry acting on any unmethylated `C` -- nothing about false
/// conversion is `CpG`-specific. Scoring a non-`CpG` `C`-to-`T` as a plain
/// sequencing error puts it at `eps / 3`, which at Q40 is 120 times too
/// confident that the read carries a real `C>T` allele.
#[test]
fn a_non_cpg_c_to_t_is_a_hundred_times_likelier_than_a_sequencing_error() {
    let haplotype = Haplotype::from_ascii(b"AACAAA");
    let site = haplotype.site(2).expect("in range");
    assert_eq!(site.cpg, CpgRole::None, "no G after it, so no CpG");
    let betas = vec![Probability::ONE; haplotype.len()];
    let emission = taps(&betas, ConversionModel::taps_default());
    let qual = BaseQuality::from_byte(40);
    let eps = error_probability(qual);

    let converted = emission.match_probability(site, observe(Base::T, qual, Strand::OT));
    let as_error =
        StandardEmission::default().match_probability(site, observe(Base::T, qual, Strand::OT));
    assert!((converted - (0.004 * (1.0 - eps) + 0.996 * eps / 3.0)).abs() < 1e-15);
    assert!((as_error - eps / 3.0).abs() < 1e-15);
    assert!(
        converted / as_error > 100.0,
        "the f row is {}x the sequencing-error row, not ~1x",
        converted / as_error
    );

    // And it is still a long way from a match, which is the specificity the
    // symbol scoring lost: f is small, the row is not absent.
    let matched = emission.match_probability(site, observe(Base::C, qual, Strand::OT));
    assert!(
        matched / converted > 200.0,
        "a real C still beats a converted T by {}x",
        matched / converted
    );

    // On the other strand the chemistry cannot act and nothing changes.
    assert_eq!(
        emission.match_probability(site, observe(Base::T, qual, Strand::OB)).to_bits(),
        StandardEmission::default()
            .match_probability(site, observe(Base::T, qual, Strand::OB))
            .to_bits()
    );
}

/// The artifact floor is `eps = max(eps_from_qual, floor)` and nothing else.
///
/// Default zero, so it changes no number that existed before it; above the
/// quality's own `eps` it takes over, and below it is inert. rastair
/// measured mate-disagreement at ~0.01 against `eps / 3 = 3.3e-5` at Q40, which
/// is the gap it exists to close.
#[test]
fn the_artifact_floor_is_a_lower_bound_on_eps() {
    let haplotype = Haplotype::from_ascii(b"AACGAA");
    let plain = haplotype.site(0).expect("A");
    let cpg = haplotype.site(2).expect("C");
    let betas = vec![Probability::ONE; haplotype.len()];
    let floor = Probability::new(0.01).expect("in [0, 1]");

    assert_eq!(StandardEmission::default().artifact_floor(), Probability::ZERO);
    assert_eq!(
        TapsEmission::new(ConversionModel::taps_default(), Betas::PerSite(&betas)).artifact_floor(),
        Probability::ZERO
    );

    let floored_standard = StandardEmission::default().with_artifact_floor(floor);
    let floored_taps = taps(&betas, ConversionModel::taps_default()).with_artifact_floor(floor);

    // Q40: eps = 1e-4, well under the floor, so the floor decides.
    for (site, read) in [(plain, Base::A), (cpg, Base::T)] {
        let high = observe(read, BaseQuality::from_byte(40), Strand::OT);
        let at_floor = observe(read, BaseQuality::from_byte(20), Strand::OT); // eps = 0.01 exactly
        assert!(
            (floored_standard.match_probability(site, high)
                - StandardEmission::default().match_probability(site, at_floor))
            .abs()
                < 1e-15,
            "a floored Q40 scores as an unfloored Q20"
        );
        assert!(
            (floored_taps.match_probability(site, high)
                - taps(&betas, ConversionModel::taps_default()).match_probability(site, at_floor))
            .abs()
                < 1e-15,
            "and the same through the conversion rows"
        );
    }

    // Q10: eps = 0.1, above the floor, so the floor is inert.
    for base in Base::KNOWN {
        let low = observe(base, BaseQuality::from_byte(10), Strand::OT);
        assert_eq!(
            floored_standard.match_probability(plain, low).to_bits(),
            StandardEmission::default().match_probability(plain, low).to_bits()
        );
        assert_eq!(
            floored_taps.match_probability(cpg, low).to_bits(),
            taps(&betas, ConversionModel::taps_default()).match_probability(cpg, low).to_bits()
        );
    }
}
