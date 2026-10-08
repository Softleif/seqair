use crate::{
    error::Error,
    transitions::{Growth, MIN_GAP_OPEN_QUALITY, Transition},
    types::{ErrorTable, Observation},
};
use seqair_types::{Base, BaseQuality, QPos, Strand};

/// A read with the four per-base quality tracks GATK's pair-HMM needs.
///
/// Every `10^(-Q/10)` the kernels need -- the base error probabilities and
/// the five transition probabilities per base -- is computed here, once, so
/// that scoring the read against a haplotype does no transcendental
/// arithmetic at all.
///
/// Every quality but [`BaseQuality::UNAVAILABLE`] is accepted, and two are
/// read as less than they claim: insertion and deletion qualities below Q6
/// count as Q6, as GATK raises them before its pair-HMM, so a match never
/// hands on more than it holds; and a base quality below Q2 counts as
/// `eps = 3/4` ([`MAX_EPSILON`]), a base that carries no information. The
/// accessors return the qualities as given.
///
/// [`MAX_EPSILON`]: crate::MAX_EPSILON
#[derive(Debug, Clone, PartialEq)]
pub struct Read {
    bases: Box<[Base]>,
    /// The base, insertion, deletion and gap-continuation qualities, one
    /// track of `len()` after another: one allocation instead of four.
    quals: Box<[BaseQuality]>,
    error_probabilities: Box<[f64]>,
    transitions: Box<[Transition]>,
    /// [`Read::growth_bound`], taken once per read rather than once per pair.
    growth_bound: Option<Growth>,
    strand: Strand,
}

/// The tracks in `Read::quals`, in order, by the name errors give them.
const TRACKS: [&str; 4] = ["base_quals", "insertion_quals", "deletion_quals", "gap_quals"];

impl Read {
    pub fn new(
        bases: impl Into<Box<[Base]>>,
        base_quals: &[BaseQuality],
        insertion_quals: &[BaseQuality],
        deletion_quals: &[BaseQuality],
        gap_quals: &[BaseQuality],
        strand: Strand,
    ) -> Result<Self, Error> {
        if strand == Strand::Unknown {
            return Err(Error::UnknownStrand);
        }
        let bases = bases.into();
        let tracks = [base_quals, insertion_quals, deletion_quals, gap_quals];
        let expected = bases.len();
        for (field, track) in TRACKS.into_iter().zip(tracks) {
            let actual = track.len();
            if actual != expected {
                return Err(Error::ReadLengthMismatch { field, expected, actual });
            }
        }
        Self::assemble(bases, tracks.concat(), strand)
    }

    /// A read whose insertion, deletion and gap-continuation qualities are the
    /// same at every base, which is what a caller with no per-base indel model
    /// has.
    pub fn uniform(
        bases: impl Into<Box<[Base]>>,
        base_quals: &[BaseQuality],
        insertion_qual: BaseQuality,
        deletion_qual: BaseQuality,
        gap_qual: BaseQuality,
        strand: Strand,
    ) -> Result<Self, Error> {
        Self::with_indel_tracks(bases, base_quals, strand, |_, insertion, deletion, gap| {
            insertion.fill(insertion_qual);
            deletion.fill(deletion_qual);
            gap.fill(gap_qual);
        })
    }

    /// A read whose insertion and deletion qualities come from `model` (one
    /// gap-open track for both, from the tandem repeats in its own bases) and
    /// whose gap-continuation quality is `gap_qual` at every base.
    pub fn with_pcr_indel_model(
        bases: impl Into<Box<[Base]>>,
        base_quals: &[BaseQuality],
        model: &PcrIndelModel,
        gap_qual: BaseQuality,
        strand: Strand,
    ) -> Result<Self, Error> {
        Self::with_indel_tracks(bases, base_quals, strand, |bases, insertion, deletion, gap| {
            model.fill_gap_open(bases, insertion);
            deletion.copy_from_slice(insertion);
            gap.fill(gap_qual);
        })
    }

    /// A read whose three indel tracks `fill` writes straight into the read's
    /// own storage, each as long as the read.
    fn with_indel_tracks(
        bases: impl Into<Box<[Base]>>,
        base_quals: &[BaseQuality],
        strand: Strand,
        fill: impl FnOnce(&[Base], &mut [BaseQuality], &mut [BaseQuality], &mut [BaseQuality]),
    ) -> Result<Self, Error> {
        if strand == Strand::Unknown {
            return Err(Error::UnknownStrand);
        }
        let bases = bases.into();
        let n = bases.len();
        if base_quals.len() != n {
            return Err(Error::ReadLengthMismatch {
                field: "base_quals",
                expected: n,
                actual: base_quals.len(),
            });
        }
        let mut quals = vec![BaseQuality::UNAVAILABLE; n.saturating_mul(4)];
        let (base, indels) = quals.split_at_mut(n);
        base.copy_from_slice(base_quals);
        let (insertion, rest) = indels.split_at_mut(n);
        let (deletion, gap) = rest.split_at_mut(n);
        fill(&bases, insertion, deletion, gap);
        Self::assemble(bases, quals, strand)
    }

    /// The read over `quals`, its four tracks back to back and each as long
    /// as `bases`.
    fn assemble(
        bases: Box<[Base]>,
        quals: Vec<BaseQuality>,
        strand: Strand,
    ) -> Result<Self, Error> {
        let n = bases.len();
        debug_assert_eq!(quals.len(), n * 4, "four tracks of the read's length");
        for (field, track) in TRACKS.into_iter().zip(quals.chunks(n.max(1))) {
            if let Some(index) = track.iter().position(|q| *q == BaseQuality::UNAVAILABLE) {
                return Err(Error::MissingQuality { field, index });
            }
        }
        let track = |k: usize| quals.get(k * n..(k + 1) * n).unwrap_or_default();
        let (base_quals, insertion_quals, deletion_quals, gap_quals) =
            (track(0), track(1), track(2), track(3));
        let table = ErrorTable::get();
        let transitions: Box<[Transition]> = insertion_quals
            .iter()
            .zip(deletion_quals)
            .zip(gap_quals)
            .map(|((insertion, deletion), gap)| {
                Transition::from_qualities(table, *insertion, *deletion, *gap)
            })
            .collect();
        let error_probabilities = base_quals.iter().map(|q| table.probability(*q)).collect();
        let growth_bound = steady_growth(&transitions, deletion_quals, gap_quals);
        Ok(Self {
            bases,
            quals: quals.into_boxed_slice(),
            error_probabilities,
            transitions,
            growth_bound,
            strand,
        })
    }

    /// Quality track `k` of [`TRACKS`].
    fn track(&self, k: usize) -> &[BaseQuality] {
        let n = self.bases.len();
        self.quals.get(k * n..(k + 1) * n).unwrap_or_default()
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.bases.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.bases.is_empty()
    }

    pub fn bases(&self) -> &[Base] {
        &self.bases
    }

    pub fn strand(&self) -> Strand {
        self.strand
    }

    #[must_use]
    pub fn base_quals(&self) -> &[BaseQuality] {
        self.track(0)
    }

    #[must_use]
    pub fn insertion_quals(&self) -> &[BaseQuality] {
        self.track(1)
    }

    #[must_use]
    pub fn deletion_quals(&self) -> &[BaseQuality] {
        self.track(2)
    }

    #[must_use]
    pub fn gap_quals(&self) -> &[BaseQuality] {
        self.track(3)
    }

    /// The transitions out of read base `index`, which the DP applies on row
    /// `index + 1`.
    pub(crate) fn transition(&self, index: usize) -> Option<Transition> {
        self.transitions.get(index).copied()
    }

    /// What the recurrence can grow a value by over this read's rows: one
    /// pass over them, which is why it is not taken for every pair.
    pub(crate) fn growth(&self) -> Growth {
        Growth::of(&self.transitions)
    }

    /// A bound on [`Read::growth`] from three scans of quality bytes, where
    /// the gap-continuation quality is the same at every base (as it is for
    /// a caller with no per-base gap model); `None` where it is not.
    pub(crate) fn growth_bound(&self) -> Option<Growth> {
        self.growth_bound
    }

    #[must_use]
    pub fn observation(&self, index: usize) -> Option<Observation> {
        let base = *self.bases.get(index)?;
        let error_probability = *self.error_probabilities.get(index)?;
        let position = u32::try_from(index).ok()?;
        Some(Observation {
            index: QPos::new(position),
            base,
            error_probability,
            strand: self.strand,
        })
    }
}

/// [`Read::growth_bound`] for these rows.
fn steady_growth(
    transitions: &[Transition],
    deletion_quals: &[BaseQuality],
    gap_quals: &[BaseQuality],
) -> Option<Growth> {
    let (Some(first), Some(&gap)) = (transitions.first(), gap_quals.first()) else {
        return None;
    };
    // One pass, and no early exit, so that it is a vector loop.
    let gap = gap.as_byte();
    let (mut differ, mut low, mut high) = (0u8, u8::MAX, 0u8);
    for (deletion, continuation) in deletion_quals.iter().zip(gap_quals) {
        differ |= continuation.as_byte() ^ gap;
        low = low.min(deletion.as_byte());
        high = high.max(deletion.as_byte());
    }
    if differ != 0 {
        return None;
    }
    let table = ErrorTable::get();
    let open = |q: u8| table.probability(BaseQuality::from_byte(q.max(MIN_GAP_OPEN_QUALITY)));
    Some(Growth::steady(transitions.len(), open(low), open(high), first.indel_to_match))
}

/// GATK's PCR indel model: a gap opens more easily inside a tandem repeat,
/// where polymerase slippage makes stutter, and the more repetitions the
/// easier. At a base inside `r` repetitions of some unit the gap-open quality
/// is `gap_open - exp(r / (rate * ln 2)) + 1`, rounded, and kept between
/// `floor` and `gap_open`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PcrIndelModel {
    /// Gap-open quality outside any repeat.
    pub gap_open: BaseQuality,
    /// The lowest a repeat can bring the gap-open quality: the library's
    /// stutter rate in a long tract.
    pub floor: BaseQuality,
    /// How slowly the quality falls with repetitions; larger is slower.
    pub rate: f64,
    /// Longest tandem unit looked for.
    pub max_unit: usize,
    /// Repetitions beyond this count as this many.
    pub max_repeats: usize,
}

impl Default for PcrIndelModel {
    /// GATK's conservative setting: Q45 outside repeats, a Q10 floor, rate 3,
    /// units up to 8 bases, at most 20 repetitions counted.
    fn default() -> Self {
        Self {
            gap_open: BaseQuality::from_byte(45),
            floor: BaseQuality::from_byte(10),
            rate: 3.0,
            max_unit: 8,
            max_repeats: 20,
        }
    }
}

impl PcrIndelModel {
    /// The gap-open quality at a base inside `repeats` repetitions of a unit.
    #[must_use]
    #[expect(clippy::cast_possible_truncation, reason = "rounded and clamped into u8 first")]
    #[expect(clippy::cast_sign_loss, reason = "clamped to at least the floor, a u8")]
    pub fn gap_open_qual(&self, repeats: usize) -> BaseQuality {
        let open = f64::from(self.gap_open.as_byte());
        let floor = f64::from(self.floor.as_byte());
        let repeats = f64::from(u32::try_from(repeats.min(self.max_repeats)).unwrap_or(u32::MAX));
        let quality = open - (repeats / (self.rate * std::f64::consts::LN_2)).exp() + 1.0;
        // `max` then `min`, not `clamp`: a floor above the open quality must
        // not panic, and a NaN rate lands on the floor.
        BaseQuality::from_byte(quality.round().max(floor).min(open) as u8)
    }

    /// Write the gap-open quality of each base of `bases` into the matching
    /// slot of `track`. A base in no repeat counts as one repetition (Q44 by
    /// default, not Q45); one inside several runs takes the most repetitions
    /// any of them makes.
    pub fn fill_gap_open(&self, bases: &[Base], track: &mut [BaseQuality]) {
        track.fill(self.gap_open_qual(1));
        self.runs(bases, |run, repeats| {
            let quality = self.gap_open_qual(repeats);
            for slot in track.iter_mut().take(run.end).skip(run.start) {
                if quality.as_byte() < slot.as_byte() {
                    *slot = quality;
                }
            }
        });
    }

    /// Every tandem run in `bases` with its repetition count: for each unit
    /// length `u` from 1 to `max_unit`, a maximal stretch where each base
    /// equals the one `u` before it, together with the `u` bases it repeats.
    fn runs(&self, bases: &[Base], mut visit: impl FnMut(std::ops::Range<usize>, usize)) {
        for unit in 1..=self.max_unit {
            let repeats_back = |k: usize| k.checked_sub(unit).and_then(|j| bases.get(j));
            let mut i = unit;
            while i < bases.len() {
                if bases.get(i) != repeats_back(i) {
                    i += 1;
                    continue;
                }
                let mut end = i;
                while end + 1 < bases.len() && bases.get(end + 1) == repeats_back(end + 1) {
                    end += 1;
                }
                visit(i - unit..end + 1, (end - i + 1) / unit + 1);
                i = end + 1;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bases(ascii: &[u8]) -> Vec<Base> {
        ascii.iter().map(|&b| Base::from(b)).collect()
    }

    /// For every base, the most repetitions any run through it makes.
    fn repeats_per_base(model: &PcrIndelModel, bases: &[Base]) -> Vec<usize> {
        let mut best = vec![1; bases.len()];
        model.runs(bases, |run, repeats| {
            for slot in best.iter_mut().take(run.end).skip(run.start) {
                *slot = (*slot).max(repeats);
            }
        });
        best
    }

    #[test]
    fn repeats_are_counted_per_unit_through_every_base_of_the_run() {
        let r = repeats_per_base(&PcrIndelModel::default(), &bases(b"ACGTTTTTTGACACACACGT"));
        assert_eq!(r[..3], [1, 1, 1]);
        assert!(r[3..9].iter().all(|&x| x == 6), "the T run: {:?}", &r[3..9]);
        assert_eq!(r[9], 1);
        assert!(r[10..18].iter().all(|&x| x == 4), "the AC run: {:?}", &r[10..18]);
    }

    /// Gap-open falls with the repeat count, from Q45 outside any repeat to
    /// the floor in a long tract; rastair floors at Q20.
    #[test]
    fn the_gap_open_quality_falls_with_repeats_to_the_floor() {
        let model = PcrIndelModel { floor: BaseQuality::from_byte(20), ..PcrIndelModel::default() };
        let q = |repeats| model.gap_open_qual(repeats).as_byte();
        assert!(q(1) >= 44);
        assert_eq!(q(5), 35);
        assert_eq!(q(6), 28);
        assert_eq!(q(7), 20);
        assert_eq!(q(40), 20);
        // GATK's own floor is lower, so a long tract keeps falling past Q20.
        let gatk = PcrIndelModel::default();
        assert_eq!(gatk.gap_open_qual(7), BaseQuality::from_byte(17));
        assert_eq!(gatk.gap_open_qual(40), BaseQuality::from_byte(10));
    }

    /// The track is the quality of the most repetitions at each base, and a
    /// read built from the model is the read built from that track.
    #[test]
    fn a_pcr_model_read_is_the_read_of_its_track() {
        let model = PcrIndelModel { floor: BaseQuality::from_byte(20), ..PcrIndelModel::default() };
        let seq = bases(b"ACGTTTTTTGACACACACGTTTTTTTTTTTTA");
        let mut track = vec![BaseQuality::UNAVAILABLE; seq.len()];
        model.fill_gap_open(&seq, &mut track);
        let want: Vec<_> =
            repeats_per_base(&model, &seq).into_iter().map(|r| model.gap_open_qual(r)).collect();
        assert_eq!(track, want);

        let quals = vec![BaseQuality::from_byte(30); seq.len()];
        let extend = BaseQuality::from_byte(10);
        let got = Read::with_pcr_indel_model(seq.clone(), &quals, &model, extend, Strand::OB);
        let gaps = vec![extend; seq.len()];
        assert_eq!(got, Read::new(seq, &quals, &want, &want, &gaps, Strand::OB));
    }
}
