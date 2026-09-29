use crate::types::{CpgRole, HapSite, Observation};
use seqair_types::{Base, Probability, Strand};

/// The row-independent half of an emission: everything it needs about one
/// haplotype column, for a read on one strand.
///
/// A pair-HMM's emission is a function of `(site, read base, eps)`, and the
/// banded kernels evaluate it once per *cell* — the whole matrix. Splitting it
/// here is what lets them hoist the site half out of the inner loop and leave a
/// two-term `f32` select behind, and [`probability`] is the only place the
/// arithmetic lives, so the hoisted form cannot drift from the reference's.
///
/// [`probability`]: Self::probability
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SiteWeights {
    /// The haplotype base, which is also the base this site is observed as when
    /// the chemistry does *not* convert it.
    pub base: Base,
    /// What it is observed as when the chemistry does, if the read's strand can
    /// see a conversion here at all.
    pub converted: Option<Base>,
    /// The chance it does, before sequencing error.
    pub rate: f64,
}

impl SiteWeights {
    /// A site no chemistry acts on.
    #[must_use]
    pub const fn plain(base: Base) -> Self {
        Self { base, converted: None, rate: 0.0 }
    }

    /// `P(observed | this site)` at this error probability.
    ///
    /// The chemistry acts before the sequencer reads, so this marginalises
    /// over the base the site carried *after* conversion: the observation is
    /// either that base read correctly, or one of the other three misread as
    /// it.
    ///
    /// ```text
    /// P(observed) = w * (1 - eps) + (1 - w) * eps / 3
    ///           w = P(post-conversion base == observed), see `latent_weight`
    /// ```
    ///
    /// Every row sums to one over the four bases, and a site the chemistry
    /// does not touch gives `1 - eps` on a match and `eps / 3` on a mismatch
    /// exactly, so a conversion model with nothing to convert is
    /// [`StandardEmission`] bit for bit. An unknown base on either side matches
    /// at full probability, as GATK's pair-HMM does.
    ///
    /// **Not** the form Bis-SNP (Liu et al. 2012, eq. 5) and the joint model's
    /// first draft wrote, `(1 - eps) * w + eps / 3`, which adds the error term
    /// unweighted and so sums to `1 + eps / 3`: that scores a cytosine the
    /// chemistry leaves alone as *likelier* than any other matching base, by
    /// `eps / 3`, at every cytosine on the read's strand. `tests/emission.rs`
    /// pins both the marginalisation (against bsgenova's three-stage model,
    /// Feng & Gao 2024) and the size of the difference.
    #[inline]
    #[must_use]
    pub fn probability(self, observed: Base, eps: f64) -> f64 {
        let weight = self.latent_weight(observed);
        weight * (1.0 - eps) + (1.0 - weight) * (eps / 3.0)
    }

    /// The chance the base the chemistry left at this site is `observed`:
    /// the half of the emission that knows nothing about sequencing error.
    ///
    /// `rate` for the converted base, `1 - rate` for the site's own base, and
    /// one or zero where nothing converts; an unknown base on either side
    /// counts as the observed one.
    #[inline]
    #[must_use]
    pub fn latent_weight(self, observed: Base) -> f64 {
        if observed == Base::Unknown || self.base == Base::Unknown {
            return 1.0;
        }
        match self.converted {
            Some(converted) if observed == converted => self.rate,
            Some(_) if observed == self.base => 1.0 - self.rate,
            _ if observed == self.base => 1.0,
            _ => 0.0,
        }
    }
}

/// `P(observed read base | haplotype base)` for the pair-HMM's match state,
/// in two halves.
///
/// An implementation supplies what it knows about a haplotype column and what
/// it knows about a read base; [`SiteWeights::probability`] composes them,
/// and that composition is the only one there is. The reference DP asks for it
/// per cell through [`MatchProbability`]; the banded kernels call
/// [`site_weights`] once per haplotype column and [`epsilon`] once per read
/// row and compose the tracks lanewise. Nothing an implementation can write
/// makes the two disagree.
///
/// The result is used as a multiplicative weight and is never renormalised,
/// so a row that summed to more than one over the four bases would favour
/// every haplotype column it applied to, by the excess, in every alignment.
/// The composition keeps every row a distribution.
///
/// [`site_weights`]: Self::site_weights
/// [`epsilon`]: Self::epsilon
pub trait Emission {
    /// Everything about one haplotype column, for a read on this strand.
    fn site_weights(&self, site: HapSite, strand: Strand) -> SiteWeights;

    /// The error probability this model assigns one observation, which is the
    /// base quality's unless something floors it.
    fn epsilon(&self, observation: Observation) -> f64;
}

/// `P(observed | site)`: an [`Emission`]'s two halves composed by
/// [`SiteWeights::probability`].
///
/// Implemented for every `Emission` and for nothing else, so it cannot be
/// overridden. An earlier draft had this as a provided method of `Emission`,
/// which let an implementation override it -- and the banded kernels, which
/// never call it, would silently score something else.
pub trait MatchProbability {
    fn match_probability(&self, site: HapSite, observation: Observation) -> f64;
}

impl<E: Emission + ?Sized> MatchProbability for E {
    #[inline]
    fn match_probability(&self, site: HapSite, observation: Observation) -> f64 {
        self.site_weights(site, observation.strand)
            .probability(observation.base, self.epsilon(observation))
    }
}

impl<T: Emission + ?Sized> Emission for &T {
    #[inline]
    fn site_weights(&self, site: HapSite, strand: Strand) -> SiteWeights {
        (*self).site_weights(site, strand)
    }
    #[inline]
    fn epsilon(&self, observation: Observation) -> f64 {
        (*self).epsilon(observation)
    }
}

/// GATK's emission: `1 - eps` on a match, `eps / 3` on a mismatch, with an
/// unknown base on either side matching everything.
///
/// `Default` is GATK's exactly; see [`with_artifact_floor`] for the one knob.
///
/// [`with_artifact_floor`]: Self::with_artifact_floor
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct StandardEmission {
    artifact_floor: Probability,
}

impl StandardEmission {
    /// A lower bound on the per-base error probability, applied as
    /// `eps = max(eps_from_qual, floor)`.
    ///
    /// A base quality describes the sequencer's confidence and nothing else,
    /// so at Q40 it claims `eps = 1e-4` and `eps / 3 = 3.3e-5` for a
    /// substitution. rastair measured the rate at which the two mates of
    /// one template disagree about a base they both cover at **~0.01** — three
    /// orders of magnitude higher, because library preparation, alignment and
    /// the reference itself contribute errors the quality score knows nothing
    /// about. A model that believes the quality score at face value scores a
    /// single artifact read as a 1e-4 event and will pay real likelihood to
    /// explain it with a haplotype.
    ///
    /// The default is zero, so every number this crate produced before the knob
    /// existed is unchanged.
    #[must_use]
    pub const fn with_artifact_floor(mut self, artifact_floor: Probability) -> Self {
        self.artifact_floor = artifact_floor;
        self
    }

    pub const fn artifact_floor(self) -> Probability {
        self.artifact_floor
    }
}

/// `max(eps_from_qual, floor)`.
#[inline]
fn floored(observation: Observation, floor: Probability) -> f64 {
    observation.error_probability.max(*floor)
}

impl Emission for StandardEmission {
    #[inline]
    fn site_weights(&self, site: HapSite, _strand: Strand) -> SiteWeights {
        SiteWeights::plain(site.base)
    }
    #[inline]
    fn epsilon(&self, observation: Observation) -> f64 {
        floored(observation, self.artifact_floor)
    }
}

/// How reliably the chemistry converts, per library.
///
/// For TAPS, `efficiency` is the chance a methylated cytosine reads as `T` and
/// `false_conversion` the chance an unmethylated one does. Swapping the two
/// roles is what makes the same type describe bisulfite.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConversionModel {
    efficiency: Probability,
    false_conversion: Probability,
}

impl ConversionModel {
    #[must_use]
    pub const fn new(efficiency: Probability, false_conversion: Probability) -> Self {
        Self { efficiency, false_conversion }
    }

    /// The values measured from this project's spike-in controls: lambda pools
    /// to beta 0.964 and the unmethylated control to 0.0038.
    #[must_use]
    pub const fn taps_default() -> Self {
        Self {
            efficiency: Probability::new_panicky(0.96),
            false_conversion: Probability::new_panicky(0.004),
        }
    }

    pub const fn efficiency(self) -> Probability {
        self.efficiency
    }

    pub const fn false_conversion(self) -> Probability {
        self.false_conversion
    }

    /// `beta * c + (1 - beta) * f`, the chance a base at this level reads as
    /// converted before sequencing error is applied.
    #[inline]
    #[must_use]
    pub fn conversion_rate(self, beta: Probability) -> f64 {
        let beta = *beta;
        beta * *self.efficiency + (1.0 - beta) * *self.false_conversion
    }
}

/// Where a `TapsEmission` reads each `CpG`'s methylation level.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Betas<'a> {
    /// One level for every `CpG` on the haplotype.
    Uniform(Probability),
    /// One level per haplotype position. A position past the end is
    /// unmethylated, which keeps the emission total and is what a caller with
    /// no estimate for a site means.
    PerSite(&'a [Probability]),
}

impl Betas<'_> {
    #[inline]
    fn at(self, index: usize) -> Probability {
        match self {
            Self::Uniform(beta) => beta,
            Self::PerSite(betas) => betas.get(index).copied().unwrap_or(Probability::ZERO),
        }
    }
}

/// The conversion-aware emission of the joint model.
///
/// Conversion is something the chemistry does to a cytosine on the strand the
/// read reports, so the rows apply at **every** haplotype `C` read on the top
/// strand and every `G` read on the bottom strand (a `G` in reference
/// orientation is the bottom strand's `C`). What the `CpG` context changes is
/// the rate, not whether the rows apply:
///
/// ```text
/// P(T | C, OT, beta) = rate * (1 - eps) + (1 - rate) * eps / 3
/// P(C | C, OT, beta) = (1 - rate) * (1 - eps) + rate * eps / 3
///
///   in a CpG:  rate = beta * c + (1 - beta) * f
///   otherwise: rate = f
/// ```
///
/// `rate` is the mixture Bis-SNP writes for bisulfite with the chemistry's
/// roles swapped: bisulfite converts the *unmethylated* cytosine, so its
/// converted-`T` probability is `beta * gamma + (1 - beta) * (1 - alpha)` with
/// `gamma` the over-conversion and `alpha` the under-conversion rate, and TAPS
/// is the same expression at `c = gamma`, `f = 1 - alpha`. The composition
/// with sequencing error is not Bis-SNP's; see [`SiteWeights::probability`].
///
/// **Why a non-`CpG` `C` is not a plain mismatch.** `f` is the false-conversion
/// rate, and it was measured on the unmethylated pUC19 spike-in's cytosines --
/// the chemistry acting on an unmethylated `C`, with nothing `CpG`-specific
/// about it. At `f = 0.004` that is a hundred times `eps / 3` at Q40, so
/// scoring an isolated `C`-to-`T` outside a `CpG` as a sequencing error
/// overstates the evidence for a real `C>T` allele by two orders of magnitude.
/// So non-`CpG` `C` and `G` use the "no" rows too; earlier drafts of this
/// crate restricted the rows to `CpG`s, and that is now corrected.
///
/// Specificity is preserved by `f` being small, not by the rows being absent:
/// at `f = 0.004` a non-`CpG` `T` over a `C` still costs a factor of ~250
/// against a match, which is what keeps a real `C>T` visible. Every other
/// combination -- a `C` read on the bottom strand included -- falls through to
/// [`StandardEmission`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TapsEmission<'a> {
    conversion: ConversionModel,
    betas: Betas<'a>,
    artifact_floor: Probability,
}

impl<'a> TapsEmission<'a> {
    #[must_use]
    pub const fn new(conversion: ConversionModel, betas: Betas<'a>) -> Self {
        Self { conversion, betas, artifact_floor: Probability::ZERO }
    }

    /// See [`StandardEmission::with_artifact_floor`]; the floor is the same
    /// `eps = max(eps_from_qual, floor)` and reaches the conversion rows too.
    #[must_use]
    pub const fn with_artifact_floor(mut self, artifact_floor: Probability) -> Self {
        self.artifact_floor = artifact_floor;
        self
    }

    pub const fn artifact_floor(self) -> Probability {
        self.artifact_floor
    }
}

/// What this read would see if the chemistry converted this site, and whether
/// the site's own `CpG` context is the one the strand can report.
///
/// Anchored on the haplotype base, so every known base belongs to exactly one
/// strand and the two are mirror images.
/// `Strand::Unknown` is unreachable here -- [`crate::Read`] rejects it at
/// construction -- and falls through to the plain, non-converting rows.
#[inline]
fn converted_base(site: HapSite, strand: Strand) -> Option<(Base, bool)> {
    match (site.base, strand) {
        (Base::C, Strand::OT) => Some((Base::T, site.cpg == CpgRole::TopC)),
        (Base::G, Strand::OB) => Some((Base::A, site.cpg == CpgRole::BottomG)),
        _ => None,
    }
}

impl Emission for TapsEmission<'_> {
    #[inline]
    fn site_weights(&self, site: HapSite, strand: Strand) -> SiteWeights {
        let Some((converted, in_cpg)) = converted_base(site, strand) else {
            return SiteWeights::plain(site.base);
        };
        let beta = if in_cpg { self.betas.at(site.index.0 as usize) } else { Probability::ZERO };
        SiteWeights {
            base: site.base,
            converted: Some(converted),
            rate: self.conversion.conversion_rate(beta),
        }
    }

    #[inline]
    fn epsilon(&self, observation: Observation) -> f64 {
        floored(observation, self.artifact_floor)
    }
}
