use crate::{read::Read, types::error_probability};
use seqair_types::BaseQuality;

/// The five transition probabilities GATK derives from one read base's
/// insertion, deletion and gap-continuation qualities.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub(crate) struct Transition {
    pub(crate) match_to_match: f64,
    pub(crate) match_to_insertion: f64,
    pub(crate) match_to_deletion: f64,
    pub(crate) indel_to_match: f64,
    pub(crate) gap_continuation: f64,
}

/// `1 - min(1, 10^(-i/10) + 10^(-d/10))`.
///
/// Not `(1 - p_i)(1 - p_d)`: GATK sums the two error probabilities and clamps,
/// which is what makes the three transitions out of a match sum to exactly one.
fn match_to_match(insertion: BaseQuality, deletion: BaseQuality) -> f64 {
    1.0 - (error_probability(insertion) + error_probability(deletion)).min(1.0)
}

/// One entry per read row, with index 0 an unused zero row so that the DP can
/// address row `i` by the read base `i - 1` without a branch.
pub(crate) fn transitions(read: &Read) -> Vec<Transition> {
    let mut out = Vec::with_capacity(read.len() + 1);
    out.push(Transition::default());
    for index in 0..read.len() {
        let zero = BaseQuality::from_byte(0);
        let insertion = read.insertion_quals().get(index).copied().unwrap_or(zero);
        let deletion = read.deletion_quals().get(index).copied().unwrap_or(zero);
        let gap = read.gap_quals().get(index).copied().unwrap_or(zero);
        let gap_continuation = error_probability(gap);
        out.push(Transition {
            match_to_match: match_to_match(insertion, deletion),
            match_to_insertion: error_probability(insertion),
            match_to_deletion: error_probability(deletion),
            indel_to_match: 1.0 - gap_continuation,
            gap_continuation,
        });
    }
    out
}
