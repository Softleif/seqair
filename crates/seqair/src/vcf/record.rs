//! VCF record model with builder pattern. Supports all VCF field types
//! and type-safe alleles.

use super::alleles::Alleles;
use super::error::VcfHeaderError;
use super::header::VcfHeader;
use indexmap::IndexMap;
use seqair_types::{One, Pos, SmallVec, SmolStr};

// r[impl vcf_record.filters]
/// Filter status for a VCF record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Filters {
    /// All filters passed.
    Pass,
    /// One or more filters failed.
    Failed(SmallVec<SmolStr, 2>),
    /// Filters have not been applied.
    NotApplied,
}

// r[impl vcf_record.info_fields]
// r[impl vcf_record.missing_values]
// r[impl vcf_record.smallvec]
/// Value for a single INFO field.
#[derive(Debug, Clone, PartialEq)]
pub enum InfoValue {
    Integer(i32),
    Float(f32),
    Flag,
    String(SmolStr),
    IntegerArray(SmallVec<Option<i32>, 4>),
    FloatArray(SmallVec<Option<f32>, 4>),
    StringArray(SmallVec<Option<SmolStr>, 4>),
}

// r[impl vcf_record.info_fields]
/// Collection of INFO key-value pairs. Duplicate keys are replaced (last value wins).
/// Iteration order matches insertion order.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct InfoFields {
    fields: IndexMap<SmolStr, InfoValue>,
}

impl InfoFields {
    pub fn new() -> Self {
        Self::default()
    }

    /// Insert or replace an INFO field. Duplicate keys are silently replaced.
    pub fn push(&mut self, key: SmolStr, value: InfoValue) {
        self.fields.insert(key, value);
    }

    pub fn is_empty(&self) -> bool {
        self.fields.is_empty()
    }

    pub fn iter(&self) -> impl Iterator<Item = (&SmolStr, &InfoValue)> {
        self.fields.iter()
    }

    pub fn len(&self) -> usize {
        self.fields.len()
    }
}

// r[impl vcf_record.genotype_encoding]
/// Genotype for a single sample.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Genotype {
    /// Allele indices (0=REF, 1+=ALT, None=missing).
    pub alleles: SmallVec<Option<u16>, 2>,
    /// Per-separator phasing (true=phased `|`, false=unphased `/`).
    /// Length = alleles.len() - 1 for diploid+.
    pub phased: SmallVec<bool, 2>,
}

impl Genotype {
    /// Unphased diploid genotype (e.g., 0/1).
    pub fn unphased(allele0: u16, allele1: u16) -> Self {
        use seqair_types::smallvec::smallvec;
        Self { alleles: smallvec![Some(allele0), Some(allele1)], phased: smallvec![false] }
    }

    /// Phased diploid genotype (e.g., 0|1).
    pub fn phased_diploid(allele0: u16, allele1: u16) -> Self {
        use seqair_types::smallvec::smallvec;
        Self { alleles: smallvec![Some(allele0), Some(allele1)], phased: smallvec![true] }
    }

    /// Haploid genotype (e.g., 0).
    pub fn haploid(allele: u16) -> Self {
        use seqair_types::smallvec::smallvec;
        Self { alleles: smallvec![Some(allele)], phased: smallvec![] }
    }

    /// Missing genotype (./.  ).
    pub fn missing_diploid() -> Self {
        use seqair_types::smallvec::smallvec;
        Self { alleles: smallvec![None, None], phased: smallvec![false] }
    }
}

/// Value for a single FORMAT field for one sample.
#[derive(Debug, Clone, PartialEq)]
pub enum SampleValue {
    Missing,
    Integer(i32),
    Float(f32),
    String(SmolStr),
    Genotype(Genotype),
    IntegerArray(SmallVec<Option<i32>, 4>),
    FloatArray(SmallVec<Option<f32>, 4>),
}

/// Per-sample FORMAT data.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct SampleFields {
    // r[impl vcf_record.format_gt_first]
    /// FORMAT keys in order (GT first if present).
    pub format_keys: SmallVec<SmolStr, 6>,
    /// Per-sample values. Outer = sample index, inner = format field index.
    pub values: Vec<SmallVec<SampleValue, 6>>,
}

// r[impl vcf_record.fields]
// r[impl vcf_record.pos_one_based]
/// A VCF record with type-safe alleles.
#[derive(Debug, Clone, PartialEq)]
pub struct VcfRecord {
    pub contig: SmolStr,
    pub pos: Pos<One>,
    pub id: Option<SmolStr>,
    pub alleles: Alleles,
    pub qual: Option<f32>,
    pub filters: Filters,
    pub info: InfoFields,
    pub samples: SampleFields,
}

// r[impl vcf_record.builder]
/// Builder for VcfRecord.
pub struct VcfRecordBuilder {
    contig: SmolStr,
    pos: Pos<One>,
    id: Option<SmolStr>,
    alleles: Alleles,
    qual: Option<f32>,
    filters: Filters,
    info: InfoFields,
    samples: SampleFields,
}

impl VcfRecordBuilder {
    /// Start building a record with required fields.
    pub fn new(contig: impl Into<SmolStr>, pos: Pos<One>, alleles: Alleles) -> Self {
        Self {
            contig: contig.into(),
            pos,
            alleles,
            id: None,
            qual: None,
            filters: Filters::NotApplied,
            info: InfoFields::new(),
            samples: SampleFields::default(),
        }
    }

    pub fn id(mut self, id: impl Into<SmolStr>) -> Self {
        self.id = Some(id.into());
        self
    }

    pub fn qual(mut self, qual: f32) -> Self {
        self.qual = Some(qual);
        self
    }

    pub fn filter_pass(mut self) -> Self {
        self.filters = Filters::Pass;
        self
    }

    pub fn filter_failed(mut self, filters: SmallVec<SmolStr, 2>) -> Self {
        self.filters = Filters::Failed(filters);
        self
    }

    pub fn info_integer(mut self, key: impl Into<SmolStr>, value: i32) -> Self {
        self.info.push(key.into(), InfoValue::Integer(value));
        self
    }

    pub fn info_float(mut self, key: impl Into<SmolStr>, value: f32) -> Self {
        self.info.push(key.into(), InfoValue::Float(value));
        self
    }

    pub fn info_flag(mut self, key: impl Into<SmolStr>) -> Self {
        self.info.push(key.into(), InfoValue::Flag);
        self
    }

    pub fn info_string(mut self, key: impl Into<SmolStr>, value: impl std::fmt::Display) -> Self {
        self.info.push(key.into(), InfoValue::String(SmolStr::from(value.to_string())));
        self
    }

    pub fn info_integer_array(
        mut self,
        key: impl Into<SmolStr>,
        values: impl Into<SmallVec<Option<i32>, 4>>,
    ) -> Self {
        self.info.push(key.into(), InfoValue::IntegerArray(values.into()));
        self
    }

    /// INFO integer array where all values are present (no missing).
    pub fn info_integers(mut self, key: impl Into<SmolStr>, values: &[i32]) -> Self {
        self.info
            .push(key.into(), InfoValue::IntegerArray(values.iter().map(|&v| Some(v)).collect()));
        self
    }

    pub fn info_float_array(
        mut self,
        key: impl Into<SmolStr>,
        values: impl Into<SmallVec<Option<f32>, 4>>,
    ) -> Self {
        self.info.push(key.into(), InfoValue::FloatArray(values.into()));
        self
    }

    /// INFO float array where all values are present (no missing).
    pub fn info_floats(mut self, key: impl Into<SmolStr>, values: &[f32]) -> Self {
        self.info
            .push(key.into(), InfoValue::FloatArray(values.iter().map(|&v| Some(v)).collect()));
        self
    }

    pub fn format_keys(mut self, keys: &[&str]) -> Self {
        self.samples.format_keys = keys.iter().map(|k| SmolStr::from(*k)).collect();
        self
    }

    pub fn add_sample(mut self, values: impl Into<SmallVec<SampleValue, 6>>) -> Self {
        self.samples.values.push(values.into());
        self
    }

    // r[impl vcf_record.builder]
    /// Build the record, validating against the header.
    #[must_use = "build() returns the record; ignoring it discards the variant data"]
    pub fn build(self, header: &VcfHeader) -> Result<VcfRecord, VcfHeaderError> {
        // r[impl vcf_record.fields]
        header.contig_id(&self.contig)?;

        // r[impl vcf_record.sample_count]
        let expected_samples = header.samples().len();
        let actual_samples = self.samples.values.len();
        if actual_samples != 0 && actual_samples != expected_samples {
            return Err(VcfHeaderError::SampleCountMismatch {
                expected: expected_samples,
                actual: actual_samples,
            });
        }

        // r[impl vcf_record.format_gt_first]
        if let Some(pos) = self.samples.format_keys.iter().position(|k| k == "GT")
            && pos != 0
        {
            return Err(VcfHeaderError::GtNotFirst { index: pos });
        }

        // r[impl vcf_record.info_fields]
        for (key, _) in self.info.iter() {
            if !header.infos().contains_key(key) {
                return Err(VcfHeaderError::MissingInfo { id: key.clone() });
            }
        }

        // r[impl vcf_record.format_gt_first]
        for key in &self.samples.format_keys {
            if !header.formats().contains_key(key) {
                return Err(VcfHeaderError::MissingFormat { id: key.clone() });
            }
        }

        Ok(VcfRecord {
            contig: self.contig,
            pos: self.pos,
            id: self.id,
            alleles: self.alleles,
            qual: self.qual,
            filters: self.filters,
            info: self.info,
            samples: self.samples,
        })
    }

    /// Build without header validation (for testing or when validation is deferred).
    pub fn build_unchecked(self) -> VcfRecord {
        VcfRecord {
            contig: self.contig,
            pos: self.pos,
            id: self.id,
            alleles: self.alleles,
            qual: self.qual,
            filters: self.filters,
            info: self.info,
            samples: self.samples,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::header::{ContigDef, InfoDef, Number, ValueType};
    use seqair_types::{Base, smallvec::smallvec};

    fn test_header() -> VcfHeader {
        VcfHeader::builder()
            .add_contig("chr1", ContigDef { length: Some(1000) })
            .unwrap()
            .add_info(
                "DP",
                InfoDef {
                    number: Number::Count(1),
                    typ: ValueType::Integer,
                    description: SmolStr::from("Depth"),
                },
            )
            .unwrap()
            .add_sample("S1")
            .unwrap()
            .build()
            .unwrap()
    }

    // r[verify vcf_record.builder]
    #[test]
    fn build_simple_snv() {
        let header = test_header();
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let record = VcfRecordBuilder::new("chr1", Pos::<One>::new(100).unwrap(), alleles)
            .qual(30.0)
            .filter_pass()
            .info_integer("DP", 50)
            .build(&header)
            .unwrap();

        assert_eq!(record.contig, "chr1");
        assert_eq!(record.pos.get(), 100);
        assert_eq!(record.alleles.ref_text(), "A");
        assert_eq!(record.qual, Some(30.0));
        assert_eq!(record.filters, Filters::Pass);
    }

    // r[verify vcf_record.builder]
    #[test]
    fn build_reference_site() {
        let header = test_header();
        let alleles = Alleles::reference(Base::G);
        let record = VcfRecordBuilder::new("chr1", Pos::<One>::new(1).unwrap(), alleles)
            .build(&header)
            .unwrap();

        assert_eq!(record.alleles.ref_text(), "G");
        assert!(record.alleles.alt_texts().is_empty());
        assert_eq!(record.filters, Filters::NotApplied);
    }

    // r[verify vcf_record.builder]
    #[test]
    fn build_rejects_unknown_contig() {
        let header = test_header();
        let alleles = Alleles::reference(Base::A);
        let result =
            VcfRecordBuilder::new("chrX", Pos::<One>::new(1).unwrap(), alleles).build(&header);
        assert!(result.is_err());
    }

    // r[verify vcf_record.filters]
    #[test]
    fn filter_states() {
        let pass = Filters::Pass;
        let failed = Filters::Failed(smallvec![SmolStr::from("q10"), SmolStr::from("dp20")]);
        let not_applied = Filters::NotApplied;

        assert_eq!(pass, Filters::Pass);
        assert!(matches!(failed, Filters::Failed(ref v) if v.len() == 2));
        assert_eq!(not_applied, Filters::NotApplied);
    }

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

    // r[verify vcf_record.info_fields]
    #[test]
    fn info_fields_dedup() {
        let record =
            VcfRecordBuilder::new("chr1", Pos::<One>::new(1).unwrap(), Alleles::reference(Base::A))
                .info_integer("DP", 50)
                .info_integer("DP", 100) // duplicate — last wins
                .build_unchecked();

        assert_eq!(record.info.len(), 1);
        let (key, val) = record.info.iter().next().unwrap();
        assert_eq!(key.as_str(), "DP");
        assert_eq!(*val, InfoValue::Integer(100));
    }

    // r[verify vcf_record.info_fields]
    #[test]
    fn info_fields_builder() {
        let record =
            VcfRecordBuilder::new("chr1", Pos::<One>::new(1).unwrap(), Alleles::reference(Base::A))
                .info_integer("DP", 50)
                .info_flag("DB")
                .info_float("AF", 0.5)
                .build_unchecked();

        assert_eq!(record.info.len(), 3);
    }

    // r[verify vcf_record.sample_count]
    #[test]
    fn rejects_sample_count_mismatch() {
        let header = test_header(); // has 1 sample
        let result =
            VcfRecordBuilder::new("chr1", Pos::<One>::new(1).unwrap(), Alleles::reference(Base::A))
                .format_keys(&["GT"])
                .add_sample(smallvec![SampleValue::Genotype(Genotype::unphased(0, 0))])
                .add_sample(smallvec![SampleValue::Genotype(Genotype::unphased(0, 0))])
                .build(&header);
        assert!(result.is_err());
    }

    // r[verify vcf_record.info_fields]
    #[test]
    fn rejects_undeclared_info_key() {
        let header = test_header();
        let result =
            VcfRecordBuilder::new("chr1", Pos::<One>::new(1).unwrap(), Alleles::reference(Base::A))
                .info_integer("BOGUS", 42)
                .build(&header);
        assert!(matches!(result, Err(VcfHeaderError::MissingInfo { .. })));
    }

    // r[verify vcf_record.sample_count]
    #[test]
    fn allows_zero_samples_sites_only() {
        let header = test_header(); // has 1 sample
        let result =
            VcfRecordBuilder::new("chr1", Pos::<One>::new(1).unwrap(), Alleles::reference(Base::A))
                .build(&header);
        assert!(result.is_ok());
    }

    // r[verify vcf_record.format_gt_first]
    #[test]
    fn rejects_gt_not_first() {
        let header = test_header();
        let result =
            VcfRecordBuilder::new("chr1", Pos::<One>::new(1).unwrap(), Alleles::reference(Base::A))
                .format_keys(&["DP", "GT"])
                .build(&header);
        assert!(result.is_err());
    }
}
