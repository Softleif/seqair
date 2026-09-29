#![allow(dead_code, reason = "shared with the tests")]
#![allow(clippy::print_stdout, reason = "this example exists to print them")]
include!("../tests/support/mod.rs");

// `TapsEmission`'s rows beside the composition they replaced, and what the
// difference is worth through the dynamic program.
//
// `SiteWeights::probability` marginalises over the base the chemistry left:
//
// ```text
// P(converted) = rate * (1 - eps) + (1 - rate) * eps / 3
// P(original)  = (1 - rate) * (1 - eps) + rate * eps / 3
// ```
//
// The crate's first draft, following Bis-SNP's
// equation 5, added the error term bare:
//
// ```text
// P(converted) = (1 - eps) * rate + eps / 3
// P(original)  = (1 - eps) * (1 - rate) + eps / 3
// ```
//
// which sums to `1 + eps / 3` over the four bases. The excess is per column
// the strand can see converted, so it is free likelihood for every haplotype
// with a cytosine there, and the second table below scores one read against
// two haplotypes that differ by a `CpG`-destroying `C>T` to show the size of
// the bias. See `tests/emission.rs` for the pins.
//
// `cargo run -p compair --release --example emission_rows`

use compair::{Betas, CpgRole, Emission, MatchProbability, QPos, StandardEmission, TapsEmission};

/// The composition this crate no longer uses. A function rather than an
/// `Emission`: the crate offers exactly one composition (`MatchProbability`
/// cannot be overridden), so this is scored with the test suite's own
/// forward recurrence, `forward_with`.
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

fn main() {
    let haplotype = Haplotype::from_ascii(b"ACGT");
    let site = haplotype.site(1).expect("index 1 is the C of a CpG");

    for quality in [30u8, 20, 10] {
        let eps = compair::error_probability(BaseQuality::from_byte(quality));
        let observation = |base| Observation {
            index: QPos::new(0),
            base,
            error_probability: eps,
            strand: Strand::OT,
        };
        println!("\n=== Q{quality}: eps = {eps:.4}, eps/3 = {:.4} ===", eps / 3.0);
        println!(
            "standard emission, P(read C | hap C) = {:.9}",
            StandardEmission::default().match_probability(site, observation(Base::C))
        );
        println!(
            "{:>5} {:>7} {:>13} {:>13} {:>11} {:>11}",
            "beta", "rate", "P(T) marginal", "P(T) flat", "row marginal", "row flat"
        );
        for beta in [0.0, 0.25, 0.5, 0.75, 1.0] {
            let taps = TapsEmission::new(
                ConversionModel::taps_default(),
                Betas::Uniform(Probability::new(beta).expect("beta is in range")),
            );
            let rate = taps.site_weights(site, Strand::OT).rate;
            let row = |prior: &dyn Fn(HapSite, Observation) -> f64| -> f64 {
                Base::KNOWN.iter().map(|base| prior(site, observation(*base))).sum()
            };
            let marginal = |site, obs| taps.match_probability(site, obs);
            let flat = |site, obs| flat_composition(&taps, site, obs);
            println!(
                "{beta:>5.2} {rate:>7.4} {:>13.9} {:>13.9} {:>11.7} {:>11.7}",
                marginal(site, observation(Base::T)),
                flat(site, observation(Base::T)),
                row(&marginal),
                row(&flat),
            );
        }
    }

    // One read against two haplotypes that differ by a CpG-destroying C>T,
    // scored both ways: the read carries the T (cut from the variant), or the
    // C (cut from the reference).
    let with_cpg = Haplotype::from_ascii(b"TTAGCATCGGATCCGATTACAGGCATTACGGATCCAGT");
    let destroyed = Haplotype::from_ascii(b"TTAGCATTGGATCCGATTACAGGCATTACGGATCCAGT");
    assert_eq!(with_cpg.site(7).map(|site| site.cpg), Some(CpgRole::TopC));
    println!(
        "\n=== log10 P(read | C haplotype) - log10 P(read | T haplotype), beta = 1, TAPS default ===\n\
         {:>6} {:>10} {:>12} {:>12} {:>12}",
        "qual", "read has", "marginal", "flat", "flat excess"
    );
    for quality in [30u8, 20, 10] {
        for (source, label) in [(&destroyed, "T"), (&with_cpg, "C")] {
            let bases = source.bases().get(2..36).expect("in range").to_vec();
            let read = Read::uniform(
                bases.clone(),
                &vec![BaseQuality::from_byte(quality); bases.len()],
                BaseQuality::from_byte(45),
                BaseQuality::from_byte(45),
                BaseQuality::from_byte(10),
                Strand::OT,
            )
            .expect("valid");
            let taps = TapsEmission::new(
                ConversionModel::taps_default(),
                Betas::Uniform(Probability::ONE),
            );
            let margin = |prior: &dyn Fn(HapSite, Observation) -> f64| {
                forward_with(&with_cpg, &read, prior).expect("scorable")
                    - forward_with(&destroyed, &read, prior).expect("scorable")
            };
            let marginal = margin(&|site, obs| taps.match_probability(site, obs));
            let flat = margin(&|site, obs| flat_composition(&taps, site, obs));
            println!(
                "{quality:>6} {label:>10} {marginal:>12.5} {flat:>12.5} {:>12.5}",
                flat - marginal
            );
        }
    }
    println!(
        "\nWhichever base the read carries, the flat composition moves the margin toward\n\
         the C haplotype: the excess is not evidence, it is paid at every cytosine the\n\
         strand can see, in proportion to how well the read's base there is explained."
    );
}
