use crate::{
    error::Error,
    types::{CpgRole, HapPos, HapSite},
};
use seqair_types::{Base, Probability};

/// A candidate sequence reads are scored against.
///
/// `CpG` roles are derived from this sequence alone, so two haplotypes over the
/// same locus disagree about where the `CpG`s are exactly when their alleles
/// create or destroy one.
///
/// A haplotype can carry a methylation level per position
/// ([`with_betas`](Self::with_betas)); an edit keeps the levels of the bases
/// it keeps, so a haplotype with an indel applied still knows its `CpG`s'
/// levels, which [`Betas::OfHaplotype`](crate::Betas::OfHaplotype) reads.
#[derive(Debug, Clone, PartialEq)]
pub struct Haplotype {
    bases: Box<[Base]>,
    roles: Box<[CpgRole]>,
    /// Empty when the haplotype carries no levels.
    betas: Box<[Option<Probability>]>,
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
        Self { bases, roles, betas: Box::default() }
    }

    /// This haplotype carrying a methylation level per position, `None` where
    /// the caller has none. An error unless there is one entry per base.
    pub fn with_betas(
        mut self,
        betas: impl Into<Box<[Option<Probability>]>>,
    ) -> Result<Self, Error> {
        let betas = betas.into();
        if betas.len() != self.bases.len() {
            return Err(Error::BetasLengthMismatch {
                expected: self.bases.len(),
                actual: betas.len(),
            });
        }
        self.betas = betas;
        Ok(self)
    }

    /// The methylation levels, one per base, or empty.
    pub fn betas(&self) -> &[Option<Probability>] {
        &self.betas
    }

    /// `Self::new` over `bases`, carrying `betas` when this haplotype carries any.
    fn edited(&self, bases: Vec<Base>, betas: Vec<Option<Probability>>) -> Self {
        let mut edited = Self::new(bases);
        if !self.betas.is_empty() {
            edited.betas = betas.into();
        }
        edited
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
        let beta = self.betas.get(index).copied().flatten();
        let index = u32::try_from(index).ok()?;
        Some(HapSite { index: HapPos(index), base, cpg, beta })
    }

    /// This haplotype with `inserted` after the base at `anchor`; `None` when
    /// `anchor` is past its end. `CpG` roles are recomputed, so an insertion
    /// can create one.
    #[must_use]
    pub fn with_insertion(&self, anchor: usize, inserted: &[Base]) -> Option<Self> {
        let at = anchor.checked_add(1)?;
        let (head, tail) = self.bases.split_at_checked(at)?;
        let betas = if self.betas.is_empty() {
            Vec::new()
        } else {
            let (before, after) = self.betas.split_at_checked(at)?;
            concat(&[before, &vec![None; inserted.len()], after])
        };
        Some(self.edited(concat(&[head, inserted, tail]), betas))
    }

    /// This haplotype with the `len` bases after the base at `anchor` deleted;
    /// `None` when the deletion runs past its end. `CpG` roles are recomputed,
    /// so a deletion can create or destroy one.
    #[must_use]
    pub fn with_deletion(&self, anchor: usize, len: usize) -> Option<Self> {
        let at = anchor.checked_add(1)?;
        let (head, tail) = self.bases.split_at_checked(at)?;
        let betas = if self.betas.is_empty() {
            Vec::new()
        } else {
            let (before, after) = self.betas.split_at_checked(at)?;
            concat(&[before, after.get(len..)?])
        };
        Some(self.edited(concat(&[head, tail.get(len..)?]), betas))
    }

    /// The reverse complement, with `CpG` roles recomputed from the new
    /// sequence rather than mirrored from the old one.
    #[must_use]
    pub fn reverse_complement(&self) -> Self {
        let bases = self.bases.iter().rev().map(|base| base.inverse()).collect::<Vec<_>>();
        self.edited(bases, self.betas.iter().rev().copied().collect())
    }
}

/// `parts` in one allocation of exactly their length, so the box `Haplotype`
/// keeps is that allocation rather than a copy of it.
fn concat<T: Copy>(parts: &[&[T]]) -> Vec<T> {
    let mut items = Vec::with_capacity(parts.iter().map(|part| part.len()).sum());
    for part in parts {
        items.extend_from_slice(part);
    }
    items
}

#[cfg(test)]
mod beta_tests {
    use super::*;

    fn p(x: f64) -> Option<Probability> {
        Some(Probability::new(x).expect("a probability"))
    }

    #[test]
    fn edits_keep_the_levels_of_the_bases_they_keep() {
        let hap = Haplotype::from_ascii(b"ACGTACG")
            .with_betas(vec![None, p(0.9), None, None, None, p(0.1), None])
            .expect("one level per base");
        assert_eq!(hap.site(1).and_then(|s| s.beta), p(0.9));
        let inserted = hap.with_insertion(2, &[Base::T, Base::T]).expect("fits");
        assert_eq!(inserted.betas().len(), 9);
        assert_eq!(inserted.site(1).and_then(|s| s.beta), p(0.9));
        assert_eq!(inserted.site(3).and_then(|s| s.beta), None);
        assert_eq!(inserted.site(7).and_then(|s| s.beta), p(0.1));
        let deleted = hap.with_deletion(1, 3).expect("fits");
        assert_eq!(deleted.betas().len(), 4);
        assert_eq!(deleted.site(2).and_then(|s| s.beta), p(0.1));
        assert!(Haplotype::from_ascii(b"ACG").with_betas(vec![None]).is_err());
        assert!(
            Haplotype::from_ascii(b"ACG").with_deletion(0, 1).expect("fits").betas().is_empty()
        );
    }
}
