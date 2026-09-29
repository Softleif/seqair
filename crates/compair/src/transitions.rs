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

/// How much more than it received GATK's recurrence can hand on, over a read.
///
/// The transitions out of one state sum to one, but the recurrence does not
/// take them from one row. A match cell in row `i` continues into row `i + 1`
/// with that row's `match_to_match` and `match_to_insertion`, and into a
/// deletion along its own row with row `i`'s `match_to_deletion`; a deletion
/// in row `i` continues along the row with row `i`'s `gap_continuation` and
/// down with row `i + 1`'s `indel_to_match`. Where those rows' qualities
/// differ a cell can hand on more than it holds, and a cell is not a
/// probability: 150 bases whose gap-continuation quality alternates between
/// Q254 and Q0 score `+37.7` against a homopolymer. So every bound that rests
/// on "a cell is at most one" is taken times these factors instead.
///
/// The argument. Leave the deletion cells out and let a match cell reach the
/// match cells of the next row through a run of `k >= 1` deletions directly,
/// with weight `match_to_deletion(i) * gap_continuation(i)^(k - 1) *
/// indel_to_match(i + 1)`. In a band no run is longer than [`RUN`] cells, so
/// the runs sum to at most `match_to_deletion(i) * indel_to_match(i + 1) *
/// run(i)` with `run(i) = min(RUN, 1 / indel_to_match(i))`. Out of a match
/// cell of row `i` that is at most `gamma(i) = max(1, match_to_match(i + 1) +
/// match_to_insertion(i + 1) + match_to_deletion(i) * indel_to_match(i + 1) *
/// run(i))`, out of an insertion cell exactly one, and every one of these
/// steps goes from row `i` to row `i + 1`. Divide each step out of row `i` by
/// `gamma(i)` and they are a sub-stochastic chain on the match and insertion
/// cells, whose paths visit a cell at most once. So, every emission being at
/// most one, the paths from one cell to another sum to at most the product
/// of `gamma` over the rows they cross, and so do the paths from one cell to
/// the last row's match and insertion cells, of which a path visits one. A
/// deletion cell of row `i` is a sum over its row's match cells to its left,
/// so it is at most `match_to_deletion(i) * run(i)` times the largest of
/// them, and it hands on at most `indel_to_match(i + 1) * run(i)` times what
/// a match cell of row `i + 1` does.
///
/// Where the qualities are constant along the read every factor is one and
/// every cell a probability. rastair's gap-open qualities step between Q43
/// and Q20 at tandem repeats, which costs about one percent per step.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Growth {
    /// `log10` of the product of `gamma` over every row crossing of the
    /// read, `1..r`: what the floors add.
    pub(crate) log10_total: f64,
    /// The largest product of `gamma` over the crossings out of rows `R +
    /// 1..R + 8`, for `R = 0, 8, 16, ..`: from a row the strip kernels
    /// renormalise to the next. The crossing out of `R` itself is not in it,
    /// since that bound starts from row `R`'s own cells.
    pub(crate) window: f64,
    /// The largest `match_to_deletion(i) * run(i)`: a deletion cell over the
    /// largest match cell of its row.
    pub(crate) spill: f64,
    /// At least one, and the largest `indel_to_match(i + 1) * run(i)`: what a
    /// deletion cell hands on, over what a match cell of the next row does.
    pub(crate) carry: f64,
}

/// The longest run of deletions inside any band: `Band::MAX_WIDTH + 1`
/// columns hold at most this many cells right of the match that opens it.
const RUN: f64 = 1024.0;
const _: () = assert!(crate::Band::MAX_WIDTH == 1024);

/// Rows between two renormalisations of the strip kernels.
const WINDOW: usize = 8;

/// `min(RUN, 1 / indel_to_match)`: what the deletion runs of a row sum to,
/// per unit of what opens them.
fn run(indel_to_match: f64) -> f64 {
    if indel_to_match * RUN > 1.0 { indel_to_match.recip() } else { RUN }
}

impl Growth {
    /// The factors for a read whose row `i` has `transitions[i - 1]`.
    pub(crate) fn of(transitions: &[Transition]) -> Self {
        let mut growth = Self { log10_total: 0.0, window: 1.0, spill: 0.0, carry: 1.0 };
        let (mut total, mut window) = (1.0f64, 1.0f64);
        let mut previous: Option<(&Transition, f64)> = None;
        for (row, t) in (1usize..).zip(transitions) {
            // A read whose gap-continuation quality holds still divides once.
            let run = match previous {
                Some((above, run)) if above.indel_to_match.to_bits() == t.indel_to_match.to_bits() => {
                    run
                }
                _ => run(t.indel_to_match),
            };
            growth.spill = growth.spill.max(t.match_to_deletion * run);
            if let Some((above, above_run)) = previous {
                // The crossing out of row `row - 1`.
                let carry = t.indel_to_match * above_run;
                growth.carry = growth.carry.max(carry);
                let gamma =
                    (t.match_to_match + t.match_to_insertion + above.match_to_deletion * carry)
                        .max(1.0);
                total *= gamma;
                if (row - 1) % WINDOW == 0 {
                    window = 1.0;
                } else {
                    window *= gamma;
                    growth.window = growth.window.max(window);
                }
            }
            previous = Some((t, run));
        }
        growth.log10_total = total.log10();
        growth
    }

    /// A bound on [`Growth::of`] for `rows` rows whose gap-continuation
    /// quality holds still at `indel_to_match`, and whose `match_to_deletion`
    /// stays within `low..=high`.
    ///
    /// With the gap continuation steady, `carry = indel_to_match * run` is at
    /// most one, so `gamma(i) <= 1 + max(0, match_to_deletion(i) -
    /// match_to_deletion(i + 1))`, at most `exp` of that. Over a stretch of
    /// crossings the product is at most `exp` of how far `match_to_deletion`
    /// falls along it, which is at most `high - low` in each of the at most
    /// `ceil(n / 2)` separate descents that `n` crossings hold.
    #[allow(clippy::cast_precision_loss, reason = "a count of rows, far below 2^52")]
    pub(crate) fn steady(rows: usize, high: f64, low: f64, indel_to_match: f64) -> Self {
        let fall = |crossings: usize| (high - low) * crossings.div_ceil(2) as f64;
        let (total, window) = (fall(rows.saturating_sub(1)), fall(WINDOW - 1));
        Self {
            log10_total: total * core::f64::consts::LOG10_E,
            window: total.min(window).exp(),
            spill: high * run(indel_to_match),
            carry: 1.0,
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
