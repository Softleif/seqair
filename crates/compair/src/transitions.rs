use crate::types::ErrorTable;
use seqair_types::BaseQuality;

/// The five transition probabilities GATK derives from one read base's
/// insertion, deletion and gap-continuation qualities.
///
/// Computed once per base when a [`Read`] is built, because the kernels score
/// one read against several haplotypes. Each probability is a lookup in the
/// [`ErrorTable`], not a `powf`.
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

/// The lowest gap-open quality the model takes at face value: GATK's
/// `QualityUtils.MIN_USABLE_Q_SCORE`, to which its
/// `PairHMMLikelihoodCalculationEngine` raises every insertion and deletion
/// quality below it.
///
/// Below it the two gap-open probabilities can sum past one -- at Q0 each is
/// one -- and a match would hand on more than it holds. At Q6 the two sum to
/// at most a half, so `match_to_match` is never below a half either.
pub(crate) const MIN_GAP_OPEN_QUALITY: u8 = 6;

impl Transition {
    /// `match_to_match` is `1 - (10^(-i/10) + 10^(-d/10))`, with `i` and `d`
    /// raised to [`MIN_GAP_OPEN_QUALITY`].
    ///
    /// Not `(1 - p_i)(1 - p_d)`: GATK sums the two error probabilities, which
    /// is what makes the three transitions out of a match sum to exactly one.
    pub(crate) fn from_qualities(
        table: &ErrorTable,
        insertion: BaseQuality,
        deletion: BaseQuality,
        gap: BaseQuality,
    ) -> Self {
        let gap_open = |quality: BaseQuality| {
            table.probability(BaseQuality::from_byte(quality.as_byte().max(MIN_GAP_OPEN_QUALITY)))
        };
        let match_to_insertion = gap_open(insertion);
        let match_to_deletion = gap_open(deletion);
        let gap_continuation = table.probability(gap);
        Self {
            match_to_match: 1.0 - (match_to_insertion + match_to_deletion),
            match_to_insertion,
            match_to_deletion,
            indel_to_match: 1.0 - gap_continuation,
            gap_continuation,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::Transition;
    use crate::types::ErrorTable;
    use seqair_types::BaseQuality;

    /// Out of every state the transitions are a distribution, at every
    /// quality byte a read can carry. The kernels' precision bound rests on
    /// it: with no state handing on more mass than it holds, every cell is a
    /// probability, and what a flush loses is bounded by where it flushes.
    #[test]
    fn every_state_hands_on_exactly_what_it_holds() {
        let table = ErrorTable::get();
        let gap = BaseQuality::from_byte(10);
        for insertion in 0..=u8::MAX {
            for deletion in 0..=u8::MAX {
                let t = Transition::from_qualities(
                    table,
                    BaseQuality::from_byte(insertion),
                    BaseQuality::from_byte(deletion),
                    gap,
                );
                let out = [t.match_to_match, t.match_to_insertion, t.match_to_deletion];
                assert!(
                    out.iter().all(|p| (0.0..=1.0).contains(p)),
                    "Q{insertion}/Q{deletion}: {out:?}"
                );
                let sum: f64 = out.iter().sum();
                assert!((sum - 1.0).abs() <= 2.0 * f64::EPSILON, "Q{insertion}/Q{deletion}: {sum}");
            }
        }
        for byte in 0..=u8::MAX {
            let q = BaseQuality::from_byte(byte);
            let t = Transition::from_qualities(table, q, q, q);
            assert!((0.0..=1.0).contains(&t.gap_continuation), "Q{byte}");
            assert_eq!(t.indel_to_match + t.gap_continuation, 1.0, "Q{byte}");
        }
    }

    /// GATK's `PairHMMLikelihoodCalculationEngine` raises insertion and
    /// deletion qualities below `QualityUtils.MIN_USABLE_Q_SCORE` (6) to 6
    /// before its pair-HMM sees them; compair does the same, so a gap-open
    /// quality under 6 means what Q6 means.
    #[test]
    fn a_gap_open_quality_below_six_counts_as_six() {
        let table = ErrorTable::get();
        let six = BaseQuality::from_byte(6);
        let gap = BaseQuality::from_byte(10);
        for byte in 0..6 {
            let low = BaseQuality::from_byte(byte);
            for (insertion, deletion, floored_insertion, floored_deletion) in
                [(low, six, six, six), (six, low, six, six), (low, low, six, six)]
            {
                assert_eq!(
                    Transition::from_qualities(table, insertion, deletion, gap),
                    Transition::from_qualities(table, floored_insertion, floored_deletion, gap),
                    "Q{byte}"
                );
            }
        }
        // And nothing at or above 6 moves.
        let seven = BaseQuality::from_byte(7);
        let t = Transition::from_qualities(table, seven, six, gap);
        assert_eq!(t.match_to_insertion.to_bits(), 10f64.powf(-0.7).to_bits());
        assert_eq!(t.match_to_deletion.to_bits(), 10f64.powf(-0.6).to_bits());
    }
}
