//! VCF record types. Contains only the [`Genotype`] type used by
//! the record encoder for GT field encoding.

use seqair_types::SmallVec;

// r[impl vcf_record.genotype_encoding]
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
        Self { alleles: SmallVec::from([None, None]), phased: SmallVec::from([false]) }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // r[verify vcf_record.genotype_encoding]
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
    }
}
