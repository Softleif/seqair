use crate::types::{CpgRole, HapPos, HapSite};
use seqair_types::Base;

/// A candidate sequence reads are scored against.
///
/// `CpG` roles are derived from this sequence alone, so two haplotypes over the
/// same locus disagree about where the `CpG`s are exactly when their alleles
/// create or destroy one.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Haplotype {
    bases: Box<[Base]>,
    roles: Box<[CpgRole]>,
}

impl Haplotype {
    #[must_use]
    pub fn new(bases: impl Into<Box<[Base]>>) -> Self {
        let bases = bases.into();
        let roles = bases
            .iter()
            .enumerate()
            .map(|(index, base)| match base {
                Base::C if bases.get(index + 1) == Some(&Base::G) => CpgRole::TopC,
                Base::G if index > 0 && bases.get(index - 1) == Some(&Base::C) => CpgRole::BottomG,
                _ => CpgRole::None,
            })
            .collect();
        Self { bases, roles }
    }

    /// Upper- and lowercase `ACGT` are the four nucleotides -- a soft-masked
    /// reference base is still that base -- and anything else, `N` and IUPAC
    /// codes included, becomes `Unknown`, which matches every base.
    #[must_use]
    pub fn from_ascii(sequence: &[u8]) -> Self {
        Self::new(sequence.iter().copied().map(Base::from).collect::<Vec<_>>())
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

    #[must_use]
    pub fn site(&self, index: usize) -> Option<HapSite> {
        let base = *self.bases.get(index)?;
        let cpg = *self.roles.get(index)?;
        let index = u32::try_from(index).ok()?;
        Some(HapSite { index: HapPos(index), base, cpg })
    }

    /// This haplotype with `inserted` after the base at `anchor`; `None` when
    /// `anchor` is past its end. `CpG` roles are recomputed, so an insertion
    /// can create one.
    #[must_use]
    pub fn with_insertion(&self, anchor: usize, inserted: &[Base]) -> Option<Self> {
        let (head, tail) = self.bases.split_at_checked(anchor.checked_add(1)?)?;
        Some(Self::new(concat(&[head, inserted, tail])))
    }

    /// This haplotype with the `len` bases after the base at `anchor` deleted;
    /// `None` when the deletion runs past its end. `CpG` roles are recomputed,
    /// so a deletion can create or destroy one.
    #[must_use]
    pub fn with_deletion(&self, anchor: usize, len: usize) -> Option<Self> {
        let (head, tail) = self.bases.split_at_checked(anchor.checked_add(1)?)?;
        Some(Self::new(concat(&[head, tail.get(len..)?])))
    }

    /// The reverse complement, with `CpG` roles recomputed from the new
    /// sequence rather than mirrored from the old one.
    #[must_use]
    pub fn reverse_complement(&self) -> Self {
        Self::new(self.bases.iter().rev().map(|base| base.inverse()).collect::<Vec<_>>())
    }
}

/// `parts` in one allocation of exactly their length, so the box `Haplotype`
/// keeps is that allocation rather than a copy of it.
fn concat(parts: &[&[Base]]) -> Vec<Base> {
    let mut bases = Vec::with_capacity(parts.iter().map(|part| part.len()).sum());
    for part in parts {
        bases.extend_from_slice(part);
    }
    bases
}
