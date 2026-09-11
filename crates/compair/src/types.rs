use seqair_types::{Base, BaseQuality, QPos, Strand};

/// `10^(-Q/10)` for a base quality, the error probability the pair-HMM works
/// in.
///
/// [`BaseQuality::UNAVAILABLE`] has no Phred value; the constructors reject it,
/// so this reads the wire byte and every caller has already checked.
#[must_use]
pub fn error_probability(qual: BaseQuality) -> f64 {
    10f64.powf(-f64::from(qual.as_byte()) / 10.0)
}

/// A likelihood in **log base 10**, the scale GATK's test vectors use.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Log10Likelihood(f64);

impl Log10Likelihood {
    /// The score of a pair no alignment can relate, `-inf`.
    pub const IMPOSSIBLE: Self = Self(f64::NEG_INFINITY);

    #[must_use]
    pub const fn new(value: f64) -> Self {
        Self(value)
    }

    #[must_use]
    pub const fn get(self) -> f64 {
        self.0
    }
}

/// Index of a base in a haplotype, zero-based.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct HapPos(pub u32);

/// What role a haplotype position plays in a `CpG`, read off that haplotype's
/// own sequence.
///
/// This is the whole of de-novo `CpG` handling: a haplotype carrying a `T>C`
/// allele that creates a `CpG` reports `TopC` there, and one carrying `C>A`
/// stops reporting it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CpgRole {
    /// Not a `CpG` cytosine or guanine.
    None,
    /// A `C` whose next base is `G`.
    TopC,
    /// A `G` whose previous base is `C`.
    BottomG,
}

impl CpgRole {
    #[must_use]
    pub const fn mirror(self) -> Self {
        match self {
            Self::None => Self::None,
            Self::TopC => Self::BottomG,
            Self::BottomG => Self::TopC,
        }
    }
}

/// One haplotype position, as the emission model sees it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct HapSite {
    pub index: HapPos,
    pub base: Base,
    pub cpg: CpgRole,
}

/// One read base, as the emission model sees it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Observation {
    pub index: QPos,
    pub base: Base,
    /// `10^(-Q/10)` for this base, precomputed once per read.
    pub error_probability: f64,
    pub strand: Strand,
}
