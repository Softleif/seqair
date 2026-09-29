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
    base_quals: Box<[BaseQuality]>,
    error_probabilities: Box<[f64]>,
    insertion_quals: Box<[BaseQuality]>,
    deletion_quals: Box<[BaseQuality]>,
    gap_quals: Box<[BaseQuality]>,
    transitions: Box<[Transition]>,
    /// [`Read::growth_bound`], taken once per read rather than once per pair.
    growth_bound: Option<Growth>,
    strand: Strand,
}

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
        let expected = bases.len();
        for (field, track) in [
            ("base_quals", base_quals),
            ("insertion_quals", insertion_quals),
            ("deletion_quals", deletion_quals),
            ("gap_quals", gap_quals),
        ] {
            let actual = track.len();
            if actual != expected {
                return Err(Error::ReadLengthMismatch { field, expected, actual });
            }
            if let Some(index) = track.iter().position(|q| *q == BaseQuality::UNAVAILABLE) {
                return Err(Error::MissingQuality { field, index });
            }
        }
        let table = ErrorTable::get();
        let transitions: Box<[Transition]> = insertion_quals
            .iter()
            .zip(deletion_quals)
            .zip(gap_quals)
            .map(|((insertion, deletion), gap)| {
                Transition::from_qualities(table, *insertion, *deletion, *gap)
            })
            .collect();
        Ok(Self {
            bases,
            base_quals: base_quals.into(),
            error_probabilities: base_quals.iter().map(|q| table.probability(*q)).collect(),
            insertion_quals: insertion_quals.into(),
            deletion_quals: deletion_quals.into(),
            gap_quals: gap_quals.into(),
            growth_bound: steady_growth(&transitions, deletion_quals, gap_quals),
            transitions,
            strand,
        })
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
        let bases = bases.into();
        let n = bases.len();
        Self::new(
            bases,
            base_quals,
            &vec![insertion_qual; n],
            &vec![deletion_qual; n],
            &vec![gap_qual; n],
            strand,
        )
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
        &self.base_quals
    }

    #[must_use]
    pub fn insertion_quals(&self) -> &[BaseQuality] {
        &self.insertion_quals
    }

    #[must_use]
    pub fn deletion_quals(&self) -> &[BaseQuality] {
        &self.deletion_quals
    }

    #[must_use]
    pub fn gap_quals(&self) -> &[BaseQuality] {
        &self.gap_quals
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
