//! VCF record types. Contains only the [`Genotype`] type used by
//! the record encoder for GT field encoding.

use seqair_types::SmallVec;

// r[impl vcf_record.genotype_encoding+2]
/// Genotype for a single sample.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Genotype {
    /// Allele indices (0=REF, 1+=ALT, None=missing).
    pub alleles: SmallVec<Option<u16>, 2>,
    /// Per-separator phasing (true=phased `|`, false=unphased `/`).
    /// Length = `alleles.len()` - 1 for diploid+.
    pub phased: SmallVec<bool, 2>,
}

impl Genotype {
    /// Unphased diploid genotype (e.g., 0/1).
    pub fn unphased(allele0: u16, allele1: u16) -> Self {
        Self {
            alleles: SmallVec::from([Some(allele0), Some(allele1)]),
            phased: SmallVec::from([false]),
        }
    }

    /// Phased diploid genotype (e.g., 0|1).
    pub fn phased_diploid(allele0: u16, allele1: u16) -> Self {
        Self {
            alleles: SmallVec::from([Some(allele0), Some(allele1)]),
            phased: SmallVec::from([true]),
        }
    }

    /// Haploid genotype (e.g., 0).
    pub fn haploid(allele: u16) -> Self {
        Self { alleles: SmallVec::from([Some(allele)]), phased: SmallVec::new() }
    }

    /// Missing genotype (./.  ).
    pub fn missing_diploid() -> Self {
        Self::missing(2)
    }

    /// A missing genotype of `ploidy` alleles, unphased: `.`, `./.`, `././.`.
    pub fn missing(ploidy: usize) -> Self {
        Self {
            alleles: std::iter::repeat_n(None, ploidy).collect(),
            phased: std::iter::repeat_n(false, ploidy.saturating_sub(1)).collect(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // r[verify vcf_record.genotype_encoding+2]
    #[test]
    fn genotype_construction() {
        let gt = Genotype::unphased(0, 1);
        assert_eq!(&gt.alleles, &[Some(0), Some(1)]);
        assert_eq!(&gt.phased, &[false]);

        let gt_phased = Genotype::phased_diploid(0, 1);
        assert_eq!(&gt_phased.phased, &[true]);

        let gt_haploid = Genotype::haploid(0);
        assert_eq!(gt_haploid.alleles.len(), 1);
        assert!(gt_haploid.phased.is_empty());

        let gt_missing = Genotype::missing_diploid();
        assert_eq!(&gt_missing.alleles, &[None, None]);
        assert_eq!(gt_missing.phased.len(), 1);
        for ploidy in 0..4 {
            let gt = Genotype::missing(ploidy);
            assert!(gt.alleles.len() == ploidy && gt.alleles.iter().all(Option::is_none));
            assert_eq!(gt.phased.len(), ploidy.saturating_sub(1), "one separator between alleles");
        }
    }
}
