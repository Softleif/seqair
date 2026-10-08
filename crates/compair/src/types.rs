use std::sync::LazyLock;

use seqair_types::{Base, BaseQuality, Probability, QPos, Strand};

/// `10^(-Q/10)` for a base quality, the error probability the pair-HMM works
/// in.
///
/// [`BaseQuality::UNAVAILABLE`] has no Phred value; the constructors reject it,
/// so this reads the wire byte and every caller has already checked.
#[must_use]
pub fn error_probability(qual: BaseQuality) -> f64 {
    ErrorTable::get().probability(qual)
}

/// The `powf` every error probability is, in one place.
fn phred_to_error(byte: u8) -> f64 {
    10f64.powf(-f64::from(byte) / 10.0)
}

/// [`error_probability`] for every quality byte, computed once per process.
///
/// A quality is a byte, so there are 256 inputs, and a `Read` asks for four
/// per base: on glibc's `pow` that was ~7 us per 125 bp read, a fifth of
/// scoring it against three haplotypes. The entries are the same
/// `powf` calls, so a lookup is bit-identical to computing one.
pub(crate) struct ErrorTable([f64; 256]);

static ERROR_TABLE: LazyLock<ErrorTable> = LazyLock::new(|| {
    ErrorTable(core::array::from_fn(|index| phred_to_error(u8::try_from(index).unwrap_or(u8::MAX))))
});

impl ErrorTable {
    /// The table, built on first use. Take it once per read rather than once
    /// per base: every access to a `LazyLock` is an atomic load and a branch.
    pub(crate) fn get() -> &'static Self {
        &ERROR_TABLE
    }

    #[inline]
    pub(crate) fn probability(&self, qual: BaseQuality) -> f64 {
        let byte = qual.as_byte();
        // A `u8` indexes 256 entries, so the fallback is unreachable and the
        // compiler drops it; it is spelled out rather than indexed.
        self.0.get(usize::from(byte)).copied().unwrap_or_else(|| phred_to_error(byte))
    }
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

    /// The likelihood narrowed to `f32`, for storing many of them. Rounds to
    /// nearest; a score below `f32`'s range (around `-3.4e38`, far below any
    /// alignment's) becomes `-inf`, and `-inf` stays `-inf`.
    #[must_use]
    #[expect(clippy::cast_possible_truncation, reason = "the narrowing is the point")]
    pub const fn as_f32(self) -> f32 {
        self.0 as f32
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
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HapSite {
    pub index: HapPos,
    pub base: Base,
    pub cpg: CpgRole,
    /// The methylation level the haplotype carries for this position, if it
    /// carries any (see [`Haplotype::with_betas`](crate::Haplotype::with_betas)).
    pub beta: Option<Probability>,
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

#[cfg(test)]
mod tests {
    use super::{ErrorTable, Log10Likelihood, error_probability};

    #[test]
    fn as_f32_rounds_and_keeps_impossible_impossible() {
        assert_eq!(Log10Likelihood::new(-0.1).as_f32(), -0.1f32);
        assert_eq!(Log10Likelihood::IMPOSSIBLE.as_f32(), f32::NEG_INFINITY);
        assert_eq!(Log10Likelihood::new(-1e300).as_f32(), f32::NEG_INFINITY);
    }
    use seqair_types::BaseQuality;

    /// Every byte, against the formula written out again rather than through
    /// `phred_to_error`: the table is a cache of `powf`, bit for bit.
    #[test]
    fn the_table_is_powf_at_every_byte() {
        let table = ErrorTable::get();
        for byte in 0..=u8::MAX {
            let want = 10f64.powf(-f64::from(byte) / 10.0).to_bits();
            let qual = BaseQuality::from_byte(byte);
            assert_eq!(table.probability(qual).to_bits(), want, "Q{byte}");
            assert_eq!(error_probability(qual).to_bits(), want, "Q{byte}");
        }
    }
}
