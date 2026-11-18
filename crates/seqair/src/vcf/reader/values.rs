//! The value views a reader hands out, the same for BCF and VCF text
//! (`r[vcf_values.access]`): [`OptValues`] yields each element as `Option<T>`,
//! [`Values`] promises none is missing. Both are lazy, allocation-free,
//! `ExactSizeIterator`s with `get(i)`.

use super::bcf_value::{BcfGt, BcfNumber, BcfVec};
use super::text_value::{TextGt, TextNumber, TextVec};
use crate::vcf::record::Genotype;
use seqair_types::SmallVec;

/// A number a VCF or BCF field holds: `i32` or `f32`.
pub(crate) trait Number: BcfNumber + TextNumber {}
impl Number for i32 {}
impl Number for f32 {}

/// A vector with an element missing where [`Values`] promised none.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
#[error("element {index} of a vector that must not have missing elements is missing")]
pub struct MissingElement {
    pub index: usize,
}

/// Where a vector's elements are stored.
#[derive(Debug, Clone)]
enum Repr<'a, T> {
    Bcf(BcfVec<'a, T>),
    Text(TextVec<'a, T>),
}

// r[impl vcf_values.access]
/// A number vector as the record stores it, whatever its width on disk: each
/// element `Some(value)`, or `None` where it is missing (`.`). An
/// `ExactSizeIterator`; `get(i)` reads any element not yet yielded.
#[derive(Clone)]
pub struct OptValues<'a, T>(Repr<'a, T>);

#[expect(private_bounds, reason = "crate-internal methods; the public ones are per type")]
impl<'a, T: Number> OptValues<'a, T> {
    pub(crate) fn bcf(v: BcfVec<'a, T>) -> Self {
        Self(Repr::Bcf(v))
    }

    pub(crate) fn text(v: TextVec<'a, T>) -> Self {
        Self(Repr::Text(v))
    }

    #[expect(
        clippy::option_option,
        reason = "in range, then the element as the iterator yields it"
    )]
    fn get_inner(&self, i: usize) -> Option<Option<T>> {
        match &self.0 {
            Repr::Bcf(v) => v.get(i),
            Repr::Text(v) => v.get(i),
        }
    }

    #[expect(clippy::option_option, reason = "the end, then the element")]
    #[inline]
    fn next_inner(&mut self) -> Option<Option<T>> {
        match &mut self.0 {
            Repr::Bcf(v) => v.next(),
            Repr::Text(v) => v.next(),
        }
    }

    fn len_inner(&self) -> usize {
        match &self.0 {
            Repr::Bcf(v) => v.len(),
            Repr::Text(v) => v.size_hint().0,
        }
    }

    /// No elements, or a lone missing one: how both formats spell a field
    /// that is `.` as a whole.
    pub(crate) fn is_missing(&self) -> bool {
        match self.len_inner() {
            0 => true,
            1 => matches!(self.get_inner(0), Some(None)),
            _ => false,
        }
    }

    // r[impl vcf_values.access]
    /// The same elements, checked to have none missing.
    pub(crate) fn all_present(self) -> Result<Values<'a, T>, MissingElement> {
        let mut rest = self.clone();
        let mut index = 0usize;
        while let Some(v) = rest.next_inner() {
            if v.is_none() {
                return Err(MissingElement { index });
            }
            index = index.saturating_add(1);
        }
        Ok(Values(self))
    }
}

/// A number vector with no missing element. An `ExactSizeIterator`; `get(i)`
/// reads any element not yet yielded.
#[derive(Clone)]
pub struct Values<'a, T>(OptValues<'a, T>);

// Per element type rather than over the crate's private `Number` bound, so
// no public impl names a private trait.
macro_rules! number_views {
    ($($t:ty),*) => {$(
        impl OptValues<'_, $t> {
            /// Element `i` of those not yet yielded: `None` past the end,
            /// `Some(None)` for a missing element.
            #[allow(
                clippy::option_option,
                reason = "in range, then the element as the iterator yields it"
            )]
            pub fn get(&self, i: usize) -> Option<Option<$t>> {
                self.get_inner(i)
            }
        }

        impl Iterator for OptValues<'_, $t> {
            type Item = Option<$t>;

            #[inline]
            fn next(&mut self) -> Option<Option<$t>> {
                self.next_inner()
            }

            fn size_hint(&self) -> (usize, Option<usize>) {
                let n = self.len_inner();
                (n, Some(n))
            }
        }

        impl ExactSizeIterator for OptValues<'_, $t> {}

        // r[impl vcf_reader.debug]
        /// The elements not yet yielded.
        impl std::fmt::Debug for OptValues<'_, $t> {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                f.debug_list().entries(self.clone()).finish()
            }
        }

        /// The elements not yet yielded.
        impl std::fmt::Debug for Values<'_, $t> {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                f.debug_list().entries(self.clone()).finish()
            }
        }

        impl Values<'_, $t> {
            /// Element `i` of those not yet yielded; `None` past the end.
            pub fn get(&self, i: usize) -> Option<$t> {
                self.0.get_inner(i).flatten()
            }
        }

        impl Iterator for Values<'_, $t> {
            type Item = $t;

            #[inline]
            fn next(&mut self) -> Option<$t> {
                let v = self.0.next_inner()?;
                debug_assert!(v.is_some(), "all_present checked every element");
                v
            }

            fn size_hint(&self) -> (usize, Option<usize>) {
                self.0.size_hint()
            }
        }

        impl ExactSizeIterator for Values<'_, $t> {}
    )*};
}

number_views!(i32, f32);

/// One allele of a genotype: its index (`None` for `.`) and whether it is
/// phased with the allele before it — for the first allele, the leading
/// phase indicator (`r[vcf_record.gt_first_phase+3]`), always unphased in a
/// file older than VCF 4.4.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GtAllele {
    pub allele: Option<u16>,
    pub phased: bool,
}

/// One sample's genotype, allele by allele.
#[derive(Debug, Clone)]
pub(crate) enum GtAlleles<'a> {
    Bcf(BcfGt<'a>),
    Text(TextGt<'a>),
}

impl GtAlleles<'_> {
    /// The writer's owned genotype. It has no slot for the first allele's
    /// phase: the writer derives that from the separators.
    pub(crate) fn into_genotype(self) -> Genotype {
        let mut alleles = SmallVec::new();
        let mut phased = SmallVec::new();
        for (i, a) in self.enumerate() {
            if i > 0 {
                phased.push(a.phased);
            }
            alleles.push(a.allele);
        }
        Genotype { alleles, phased }
    }
}

impl Iterator for GtAlleles<'_> {
    type Item = GtAllele;

    #[inline]
    fn next(&mut self) -> Option<GtAllele> {
        match self {
            Self::Bcf(gt) => gt.next(),
            Self::Text(gt) => gt.next(),
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        match self {
            Self::Bcf(gt) => gt.size_hint(),
            Self::Text(gt) => gt.size_hint(),
        }
    }
}

impl ExactSizeIterator for GtAlleles<'_> {}

// r[impl vcf_values.gt+2]
/// One sample's genotype as the record stores it. Re-iterable: every
/// [`alleles`](Self::alleles) call starts again from the first allele.
#[derive(Clone)]
pub struct GenotypeRef<'a>(pub(crate) GtAlleles<'a>);

// r[impl vcf_reader.debug]
/// As VCF 4.4 writes it: `0/1`, `0|1`, `|1` for a phased haploid allele.
impl std::fmt::Debug for GenotypeRef<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut gt = String::new();
        for (i, a) in self.alleles().enumerate() {
            if i > 0 || a.phased {
                gt.push(if a.phased { '|' } else { '/' });
            }
            match a.allele {
                Some(n) => gt.push_str(&n.to_string()),
                None => gt.push('.'),
            }
        }
        f.debug_tuple("GenotypeRef").field(&gt).finish()
    }
}

impl GenotypeRef<'_> {
    /// The alleles in order, each with its phase.
    pub fn alleles(&self) -> impl ExactSizeIterator<Item = GtAllele> + Clone + '_ {
        self.0.clone()
    }

    // r[impl vcf_values.gt+2]
    /// Allele `i`, from 0; `None` past the last.
    pub fn get(&self, i: usize) -> Option<GtAllele> {
        self.0.clone().nth(i)
    }

    /// How many alleles the sample has.
    pub fn ploidy(&self) -> usize {
        self.0.len()
    }

    /// The writer's owned [`Genotype`].
    pub fn to_genotype(&self) -> Genotype {
        self.0.clone().into_genotype()
    }
}
#[cfg(test)]
mod tests {
    #![allow(clippy::arithmetic_side_effects, reason = "test arithmetic on small values")]

    use super::*;
    use crate::vcf::bcf_encoding::{
        BCF_BT_FLOAT, FLOAT_MISSING, encode_type_byte, encode_typed_int_vec,
    };
    use crate::vcf::bcf_encoding::{BCF_BT_INT32, INT32_END_OF_VECTOR};
    use crate::vcf::reader::bcf_value::{FormatColumn, TypedValue};
    use hegel::prelude::*;

    fn bcf(bytes: &[u8]) -> OptValues<'_, i32> {
        OptValues::bcf(TypedValue::parse(bytes).unwrap().0.ints().unwrap())
    }

    // r[verify vcf_values.access]
    #[test]
    fn all_present_rejects_the_first_missing_element() {
        assert_eq!(
            bcf(&[0x31, 1, 0x80, 0x80]).all_present().unwrap_err(),
            MissingElement { index: 1 }
        );
        let v = bcf(&[0x31, 1, 2, 0x81]).all_present().unwrap();
        assert_eq!((v.len(), v.get(1), v.get(2)), (2, Some(2), None));
        assert_eq!(v.collect::<Vec<_>>(), [1, 2]);
    }

    #[test]
    fn a_lone_missing_element_or_nothing_is_missing() {
        assert!(bcf(&[0x01]).is_missing());
        assert!(bcf(&[0x11, 0x80]).is_missing());
        assert!(bcf(&[0x21, 0x80, 0x81]).is_missing());
        assert!(!bcf(&[0x21, 0x80, 0x80]).is_missing());
        assert!(!bcf(&[0x11, 0]).is_missing());
    }

    fn text<T: Number>(bytes: &[u8]) -> OptValues<'_, T> {
        OptValues::text(TextVec::new(bytes).unwrap())
    }

    // r[verify vcf_values.access]
    #[test]
    fn text_views_check_missing_elements_too() {
        assert!(text::<i32>(b".").is_missing());
        assert!(!text::<i32>(b".,.").is_missing());
        assert_eq!(text::<i32>(b"1,.,3").all_present().unwrap_err(), MissingElement { index: 1 });
        assert_eq!(text::<f32>(b"1,2.5").all_present().unwrap().collect::<Vec<_>>(), [1.0, 2.5]);
    }

    // r[verify vcf_values.bcf_sentinels]
    // r[verify vcf_values.text_numbers]
    /// One float vector as BCF stores it and as VCF prints it reads the same
    /// through both views (floats by value after the text round trip).
    #[hegel::test]
    fn bcf_and_text_float_vectors_agree(tc: TestCase) {
        let finite = gs::floats::<f32>().min_value(f32::MIN).max_value(f32::MAX);
        let values = tc.draw(gs::vecs(gs::optional(finite)).min_size(1).max_size(20));
        let mut bcf_bytes = Vec::new();
        encode_type_byte(&mut bcf_bytes, values.len(), BCF_BT_FLOAT);
        for v in &values {
            bcf_bytes.extend_from_slice(&v.map_or(FLOAT_MISSING, f32::to_bits).to_le_bytes());
        }
        let text_bytes: Vec<String> =
            values.iter().map(|v| v.map_or_else(|| ".".to_owned(), |v| v.to_string())).collect();
        let text_bytes = text_bytes.join(",");

        let from_bcf: OptValues<'_, f32> =
            OptValues::bcf(TypedValue::parse(&bcf_bytes).unwrap().0.floats().unwrap());
        let from_text = text::<f32>(text_bytes.as_bytes());
        assert_eq!(from_bcf.is_missing(), from_text.is_missing());
        assert_eq!(
            from_bcf.clone().all_present().map(Iterator::collect::<Vec<_>>),
            from_text.clone().all_present().map(Iterator::collect::<Vec<_>>)
        );
        assert_eq!(from_bcf.collect::<Vec<_>>(), values);
        assert_eq!(from_text.collect::<Vec<_>>(), values);
    }

    // r[verify vcf_values.bcf_sentinels]
    // r[verify vcf_values.text_numbers]
    #[hegel::test]
    fn bcf_and_text_int_vectors_agree(tc: TestCase) {
        let values = tc
            .draw(gs::vecs(gs::integers::<i32>().min_value(i32::MIN + 8)).min_size(1).max_size(20));
        let mut bcf_bytes = Vec::new();
        encode_typed_int_vec(&mut bcf_bytes, &values);
        let text_bytes: Vec<String> = values.iter().map(i32::to_string).collect();
        let text_bytes = text_bytes.join(",");
        let from_bcf: OptValues<'_, i32> =
            OptValues::bcf(TypedValue::parse(&bcf_bytes).unwrap().0.ints().unwrap());
        let from_text = text::<i32>(text_bytes.as_bytes());
        assert_eq!(from_bcf.all_present().unwrap().collect::<Vec<_>>(), values);
        assert_eq!(from_text.all_present().unwrap().collect::<Vec<_>>(), values);
    }

    fn a(allele: Option<u16>, phased: bool) -> GtAllele {
        GtAllele { allele, phased }
    }

    fn text_gt(gt: &str, read: bool) -> Vec<GtAllele> {
        GtAlleles::Text(TextGt::new(gt.as_bytes(), read).unwrap()).collect()
    }

    /// Each sample's alleles from a GT column: its descriptor and elements.
    fn bcf_gts(bytes: &[u8], n_samples: usize, read: bool) -> Vec<Vec<GtAllele>> {
        let (column, len) = FormatColumn::parse(bytes, n_samples).unwrap();
        assert_eq!(len, bytes.len());
        (0..n_samples)
            .map(|i| GtAlleles::Bcf(column.sample(i).unwrap().gt(read).unwrap()).collect())
            .collect()
    }

    // r[verify vcf_values.gt]
    // r[verify vcf_record.gt_first_phase+3]
    /// The examples of `[VCF45]` §6.3.3's GT table, and two VCF 4.4 prefixed
    /// ones, read from text and from BCF; before 4.4 the first phase is unphased.
    #[test]
    fn the_spec_genotype_examples() {
        let n = None;
        let rows: &[(&str, &[u8], &[GtAllele])] = &[
            ("0/1", &[0x21, 0x02, 0x04], &[a(Some(0), false), a(Some(1), false)]),
            ("0|1", &[0x21, 0x03, 0x05], &[a(Some(0), true), a(Some(1), true)]),
            ("./.", &[0x21, 0x00, 0x00], &[a(n, false), a(n, false)]),
            (".", &[0x11, 0x00], &[a(n, false)]),
            ("0", &[0x11, 0x03], &[a(Some(0), true)]),
            ("1", &[0x11, 0x05], &[a(Some(1), true)]),
            (
                "0/1/2",
                &[0x31, 0x02, 0x04, 0x06],
                &[a(Some(0), false), a(Some(1), false), a(Some(2), false)],
            ),
            (
                "0/1|2",
                &[0x31, 0x02, 0x04, 0x07],
                &[a(Some(0), false), a(Some(1), false), a(Some(2), true)],
            ),
            ("|0/1", &[0x21, 0x03, 0x04], &[a(Some(0), true), a(Some(1), false)]),
            ("/0|1", &[0x21, 0x02, 0x05], &[a(Some(0), false), a(Some(1), true)]),
        ];
        for &(text, bcf, expected) in rows {
            assert_eq!(text_gt(text, true), expected, "{text}");
            assert_eq!(bcf_gts(bcf, 1, true), [expected], "{text}");
            let mut before_44 = expected.to_vec();
            before_44[0].phased = false;
            assert_eq!(text_gt(text, false), before_44, "{text} before 4.4");
            assert_eq!(bcf_gts(bcf, 1, false), [before_44], "{text} before 4.4");
        }
        // "0 and 0/1": the haploid sample is padded with end-of-vector.
        assert_eq!(
            bcf_gts(&[0x21, 0x03, 0x81, 0x02, 0x04], 2, true),
            [vec![a(Some(0), true)], vec![a(Some(0), false), a(Some(1), false)]]
        );
        // htslib fills a sample without GT with the int missing value.
        assert_eq!(bcf_gts(&[0x21, 0x80, 0x81], 1, true), [vec![a(n, false)]]);
    }

    #[test]
    fn genotypes_convert_to_the_writers_genotype() {
        let gt = |text: &str| {
            GtAlleles::Text(TextGt::new(text.as_bytes(), true).unwrap()).into_genotype()
        };
        assert_eq!(gt("0|1"), Genotype::phased_diploid(0, 1));
        assert_eq!(gt("0/1"), Genotype::unphased(0, 1));
        assert_eq!(gt("1"), Genotype::haploid(1));
        assert_eq!(gt("./."), Genotype::missing_diploid());
        assert_eq!(gt("|0/1"), Genotype::unphased(0, 1));
    }

    /// A genotype as text renders it: alleles, the phase of each separator,
    /// and an explicit leading indicator or none.
    #[derive(Debug)]
    struct Drawn {
        alleles: Vec<Option<u16>>,
        separators: Vec<bool>,
        prefix: Option<bool>,
    }

    #[hegel::composite]
    fn genotype(tc: &TestCase) -> Drawn {
        let ploidy = tc.draw(gs::integers::<usize>().min_value(1).max_value(5));
        let allele = hegel::one_of!(gs::integers::<u16>().max_value(3), gs::integers::<u16>());
        Drawn {
            alleles: tc.draw(gs::vecs(gs::optional(allele)).min_size(ploidy).max_size(ploidy)),
            separators: tc.draw(gs::vecs(gs::booleans()).min_size(ploidy - 1).max_size(ploidy - 1)),
            prefix: tc.draw(gs::optional(gs::booleans())),
        }
    }

    impl Drawn {
        fn text(&self) -> String {
            let mut s = String::new();
            if let Some(p) = self.prefix {
                s.push(if p { '|' } else { '/' });
            }
            for (i, allele) in self.alleles.iter().enumerate() {
                if i > 0 {
                    s.push(if self.separators[i - 1] { '|' } else { '/' });
                }
                match allele {
                    Some(n) => s.push_str(&n.to_string()),
                    None => s.push('.'),
                }
            }
            s
        }

        /// `[VCF45]` §1.6.2: an omitted leading indicator is `/` if any
        /// separator is `/` and `|` otherwise; htslib reads a lone `.` as unphased.
        fn expected(&self, read: bool) -> Vec<GtAllele> {
            let implied = if self.alleles.len() == 1 {
                self.alleles[0].is_some()
            } else {
                self.separators.iter().all(|&p| p)
            };
            let first = read && self.prefix.unwrap_or(implied);
            let phases = std::iter::once(first).chain(self.separators.iter().copied());
            self.alleles.iter().zip(phases).map(|(&allele, phased)| a(allele, phased)).collect()
        }
    }

    // r[verify vcf_values.gt]
    // r[verify vcf_record.gt_first_phase+3]
    /// Text genotypes read as drawn, and BCF genotypes encoded per `[VCF45]`
    /// §6.3.3 (`(allele + 1) << 1 | phased`, int32, end-of-vector padding to
    /// the longest sample) read the same, in either version mode.
    #[hegel::test]
    fn text_and_bcf_genotypes_read_as_drawn(tc: TestCase) {
        let samples = tc.draw(gs::vecs(genotype()).min_size(1).max_size(4));
        let read = tc.draw(gs::booleans());
        let ploidy = samples.iter().map(|g| g.alleles.len()).max().unwrap_or(0);
        let mut bcf = Vec::new();
        encode_type_byte(&mut bcf, ploidy, BCF_BT_INT32);
        for g in &samples {
            assert_eq!(text_gt(&g.text(), read), g.expected(read), "{}", g.text());
            for e in g.expected(true) {
                let allele = e.allele.map_or(0, |n| i32::from(n) + 1);
                bcf.extend_from_slice(&((allele << 1) | i32::from(e.phased)).to_le_bytes());
            }
            for _ in g.alleles.len()..ploidy {
                bcf.extend_from_slice(&INT32_END_OF_VECTOR.to_le_bytes());
            }
        }
        let got = bcf_gts(&bcf, samples.len(), read);
        for (g, got) in samples.iter().zip(got) {
            assert_eq!(got, g.expected(read), "{}", g.text());
        }
    }
}
