use crate::{
    error::Error,
    types::{Observation, error_probability},
};
use seqair_types::{Base, BaseQuality, QPos, Strand};

/// A read with the four per-base quality tracks GATK's pair-HMM needs.
#[derive(Debug, Clone, PartialEq)]
pub struct Read {
    bases: Box<[Base]>,
    base_quals: Box<[BaseQuality]>,
    error_probabilities: Box<[f64]>,
    insertion_quals: Box<[BaseQuality]>,
    deletion_quals: Box<[BaseQuality]>,
    gap_quals: Box<[BaseQuality]>,
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
        Ok(Self {
            bases,
            base_quals: base_quals.into(),
            error_probabilities: base_quals.iter().copied().map(error_probability).collect(),
            insertion_quals: insertion_quals.into(),
            deletion_quals: deletion_quals.into(),
            gap_quals: gap_quals.into(),
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
