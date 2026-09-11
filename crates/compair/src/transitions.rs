use crate::types::error_probability;
use seqair_types::BaseQuality;

/// The five transition probabilities GATK derives from one read base's
/// insertion, deletion and gap-continuation qualities.
///
/// Computed once per base when a [`Read`] is built, because every one of them
/// is a `powf` and the kernels score one read against several haplotypes.
///
/// [`Read`]: crate::Read
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub(crate) struct Transition {
    pub(crate) match_to_match: f64,
    pub(crate) match_to_insertion: f64,
    pub(crate) match_to_deletion: f64,
    pub(crate) indel_to_match: f64,
    pub(crate) gap_continuation: f64,
}

impl Transition {
    /// `match_to_match` is `1 - min(1, 10^(-i/10) + 10^(-d/10))`.
    ///
    /// Not `(1 - p_i)(1 - p_d)`: GATK sums the two error probabilities and
    /// clamps, which is what makes the three transitions out of a match sum
    /// to exactly one.
    pub(crate) fn from_qualities(
        insertion: BaseQuality,
        deletion: BaseQuality,
        gap: BaseQuality,
    ) -> Self {
        let match_to_insertion = error_probability(insertion);
        let match_to_deletion = error_probability(deletion);
        let gap_continuation = error_probability(gap);
        Self {
            match_to_match: 1.0 - (match_to_insertion + match_to_deletion).min(1.0),
            match_to_insertion,
            match_to_deletion,
            indel_to_match: 1.0 - gap_continuation,
            gap_continuation,
        }
    }
}
