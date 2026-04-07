//! Type-safe VCF allele representation. Each variant enforces structural
//! invariants and handles its own VCF REF/ALT serialization.

use super::error::AllelesError;
use seqair_types::{Base, SmallVec, SmolStr};

// r[impl vcf_record.alleles_typed]
/// Type-safe allele representation for VCF records.
///
/// Each variant enforces the structural rules of its VCF encoding and
/// provides serialization to REF/ALT text columns.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Alleles {
    /// Reference-only site (gVCF). REF=single base, no ALT.
    /// VCF: `A\t.`
    Reference { ref_base: Base },

    /// SNV(s) at a single position. REF=one base, each ALT=one base.
    /// Multi-allelic: REF=A, ALT=T,C
    Snv { ref_base: Base, alt_bases: SmallVec<Base, 2> },

    /// Insertion after anchor position.
    /// VCF encodes as REF=anchor, ALT=anchor+inserted.
    /// Position is at the anchor base (the base before the insertion).
    Insertion { anchor: Base, inserted: SmallVec<Base, 4> },

    /// Deletion starting after anchor.
    /// VCF encodes as REF=anchor+deleted, ALT=anchor.
    /// Position is at the anchor base.
    Deletion { anchor: Base, deleted: SmallVec<Base, 4> },

    /// Multi-allelic mixed site, complex variant (MNV), or symbolic allele.
    /// Escape hatch for anything that doesn't fit the above categories.
    Complex { ref_allele: SmolStr, alt_alleles: SmallVec<SmolStr, 2> },
}

impl Alleles {
    /// Reference-only site (gVCF).
    pub fn reference(ref_base: Base) -> Self {
        Self::Reference { ref_base }
    }

    /// Single SNV. Returns error if alt equals ref.
    pub fn snv(ref_base: Base, alt_base: Base) -> Result<Self, AllelesError> {
        if ref_base == alt_base {
            return Err(AllelesError::SnvAltEqualsRef);
        }
        Ok(Self::Snv {
            ref_base,
            alt_bases: {
                use seqair_types::smallvec::smallvec;
                smallvec![alt_base]
            },
        })
    }

    /// Multi-allelic SNV. Returns error if empty, any alt equals ref, or duplicates.
    pub fn snv_multi(ref_base: Base, alt_bases: &[Base]) -> Result<Self, AllelesError> {
        if alt_bases.is_empty() {
            return Err(AllelesError::SnvEmpty);
        }
        for alt in alt_bases {
            if *alt == ref_base {
                return Err(AllelesError::SnvAltEqualsRef);
            }
        }
        // Check duplicates (at most 4 bases, so O(n^2) is fine)
        for (i, a) in alt_bases.iter().enumerate() {
            for b in alt_bases.iter().skip(i.saturating_add(1)) {
                if a == b {
                    return Err(AllelesError::SnvDuplicateAlt);
                }
            }
        }
        Ok(Self::Snv { ref_base, alt_bases: SmallVec::from(alt_bases) })
    }

    /// Insertion after anchor base. `inserted` is the bases inserted after the anchor.
    pub fn insertion(anchor: Base, inserted: &[Base]) -> Result<Self, AllelesError> {
        if inserted.is_empty() {
            return Err(AllelesError::InsertionEmpty);
        }
        Ok(Self::Insertion { anchor, inserted: SmallVec::from(inserted) })
    }

    /// Deletion after anchor base. `deleted` is the bases being deleted.
    pub fn deletion(anchor: Base, deleted: &[Base]) -> Result<Self, AllelesError> {
        if deleted.is_empty() {
            return Err(AllelesError::DeletionEmpty);
        }
        Ok(Self::Deletion { anchor, deleted: SmallVec::from(deleted) })
    }

    /// Complex / symbolic / MNV (no structural validation).
    pub fn complex(ref_allele: SmolStr, alt_alleles: SmallVec<SmolStr, 2>) -> Self {
        Self::Complex { ref_allele, alt_alleles }
    }

    // r[impl vcf_record.alleles_serialization]
    /// VCF REF column value.
    pub fn ref_text(&self) -> SmolStr {
        match self {
            Self::Reference { ref_base } | Self::Snv { ref_base, .. } => {
                SmolStr::from(ref_base.as_str())
            }
            Self::Insertion { anchor, .. } => SmolStr::from(anchor.as_str()),
            Self::Deletion { anchor, deleted } => {
                let mut s = String::with_capacity(1usize.saturating_add(deleted.len()));
                s.push(anchor.as_char());
                for b in deleted {
                    s.push(b.as_char());
                }
                SmolStr::from(s)
            }
            Self::Complex { ref_allele, .. } => ref_allele.clone(),
        }
    }

    /// VCF ALT column entries. Empty for reference-only sites.
    pub fn alt_texts(&self) -> SmallVec<SmolStr, 2> {
        match self {
            Self::Reference { .. } => SmallVec::new(),
            Self::Snv { alt_bases, .. } => {
                alt_bases.iter().map(|b| SmolStr::from(b.as_str())).collect()
            }
            Self::Insertion { anchor, inserted } => {
                let mut s = String::with_capacity(1usize.saturating_add(inserted.len()));
                s.push(anchor.as_char());
                for b in inserted {
                    s.push(b.as_char());
                }
                {
                    use seqair_types::smallvec::smallvec;
                    smallvec![SmolStr::from(s)]
                }
            }
            Self::Deletion { anchor, .. } => {
                use seqair_types::smallvec::smallvec;
                smallvec![SmolStr::from(anchor.as_str())]
            }
            Self::Complex { alt_alleles, .. } => alt_alleles.clone(),
        }
    }

    /// Write the REF column directly into a byte buffer (zero-alloc for common cases).
    pub fn write_ref_into(&self, buf: &mut Vec<u8>) {
        match self {
            Self::Reference { ref_base } | Self::Snv { ref_base, .. } => {
                buf.push(ref_base.as_char() as u8);
            }
            Self::Insertion { anchor, .. } => {
                buf.push(anchor.as_char() as u8);
            }
            Self::Deletion { anchor, deleted } => {
                buf.push(anchor.as_char() as u8);
                for b in deleted {
                    buf.push(b.as_char() as u8);
                }
            }
            Self::Complex { ref_allele, .. } => {
                buf.extend_from_slice(ref_allele.as_bytes());
            }
        }
    }

    /// Write comma-separated ALT column directly into a byte buffer.
    /// Returns the number of alt alleles written.
    pub fn write_alts_into(&self, buf: &mut Vec<u8>) -> usize {
        match self {
            Self::Reference { .. } => 0,
            Self::Snv { alt_bases, .. } => {
                for (i, b) in alt_bases.iter().enumerate() {
                    if i > 0 {
                        buf.push(b',');
                    }
                    buf.push(b.as_char() as u8);
                }
                alt_bases.len()
            }
            Self::Insertion { anchor, inserted } => {
                buf.push(anchor.as_char() as u8);
                for b in inserted {
                    buf.push(b.as_char() as u8);
                }
                1
            }
            Self::Deletion { anchor, .. } => {
                buf.push(anchor.as_char() as u8);
                1
            }
            Self::Complex { alt_alleles, .. } => {
                for (i, alt) in alt_alleles.iter().enumerate() {
                    if i > 0 {
                        buf.push(b',');
                    }
                    buf.extend_from_slice(alt.as_bytes());
                }
                alt_alleles.len()
            }
        }
    }

    // r[impl vcf_record.alleles_rlen]
    // r[impl vcf_record.rlen]
    /// Reference length for BCF rlen field.
    pub fn rlen(&self) -> usize {
        match self {
            Self::Reference { .. } | Self::Snv { .. } | Self::Insertion { .. } => 1,
            Self::Deletion { deleted, .. } => 1usize.saturating_add(deleted.len()),
            Self::Complex { ref_allele, .. } => ref_allele.len(),
        }
    }

    /// Number of bytes `write_ref_into` writes (same as rlen for concrete alleles).
    pub fn ref_byte_len(&self) -> usize {
        self.rlen()
    }

    // r[impl vcf_record.allele_count]
    /// Number of alleles including REF (n_allele for BCF).
    pub fn n_allele(&self) -> usize {
        match self {
            Self::Reference { .. } => 1,
            Self::Snv { alt_bases, .. } => 1usize.saturating_add(alt_bases.len()),
            Self::Insertion { .. } | Self::Deletion { .. } => 2,
            Self::Complex { alt_alleles, .. } => 1usize.saturating_add(alt_alleles.len()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use seqair_types::smallvec::smallvec;

    // r[verify vcf_record.alleles_typed]
    #[test]
    fn reference_site() {
        let a = Alleles::reference(Base::A);
        assert_eq!(a.ref_text(), "A");
        assert!(a.alt_texts().is_empty());
        assert_eq!(a.rlen(), 1);
        assert_eq!(a.n_allele(), 1);
    }

    // r[verify vcf_record.alleles_typed]
    #[test]
    fn snv_basic() {
        let a = Alleles::snv(Base::A, Base::T).unwrap();
        assert_eq!(a.ref_text(), "A");
        assert_eq!(a.alt_texts().as_slice(), [SmolStr::from("T")]);
        assert_eq!(a.rlen(), 1);
        assert_eq!(a.n_allele(), 2);
    }

    // r[verify vcf_record.alleles_typed]
    #[test]
    fn snv_rejects_ref_equals_alt() {
        assert!(Alleles::snv(Base::A, Base::A).is_err());
    }

    // r[verify vcf_record.alleles_typed]
    #[test]
    fn snv_multi_allelic() {
        let a = Alleles::snv_multi(Base::A, &[Base::T, Base::C]).unwrap();
        assert_eq!(a.ref_text(), "A");
        let alts = a.alt_texts();
        assert_eq!(alts.len(), 2);
        assert_eq!(alts.as_slice(), [SmolStr::from("T"), SmolStr::from("C")]);
        assert_eq!(a.n_allele(), 3);
    }

    // r[verify vcf_record.alleles_typed]
    #[test]
    fn snv_multi_rejects_duplicates() {
        assert!(Alleles::snv_multi(Base::A, &[Base::T, Base::T]).is_err());
    }

    // r[verify vcf_record.alleles_typed]
    #[test]
    fn snv_multi_rejects_empty() {
        assert!(Alleles::snv_multi(Base::A, &[]).is_err());
    }

    // r[verify vcf_record.alleles_serialization]
    #[test]
    fn insertion_serialization() {
        // Insertion of CGT after anchor A → REF=A, ALT=ACGT
        let a = Alleles::insertion(Base::A, &[Base::C, Base::G, Base::T]).unwrap();
        assert_eq!(a.ref_text(), "A");
        assert_eq!(a.alt_texts().as_slice(), [SmolStr::from("ACGT")]);
        assert_eq!(a.rlen(), 1);
        assert_eq!(a.n_allele(), 2);
    }

    // r[verify vcf_record.alleles_typed]
    #[test]
    fn insertion_rejects_empty() {
        assert!(Alleles::insertion(Base::A, &[]).is_err());
    }

    // r[verify vcf_record.alleles_serialization]
    #[test]
    fn deletion_serialization() {
        // Deletion of CGT after anchor A → REF=ACGT, ALT=A
        let a = Alleles::deletion(Base::A, &[Base::C, Base::G, Base::T]).unwrap();
        assert_eq!(a.ref_text(), "ACGT");
        assert_eq!(a.alt_texts().as_slice(), [SmolStr::from("A")]);
        assert_eq!(a.rlen(), 4);
        assert_eq!(a.n_allele(), 2);
    }

    // r[verify vcf_record.alleles_typed]
    #[test]
    fn deletion_rejects_empty() {
        assert!(Alleles::deletion(Base::A, &[]).is_err());
    }

    // r[verify vcf_record.alleles_serialization]
    #[test]
    fn complex_alleles() {
        let a = Alleles::complex(
            SmolStr::from("ACG"),
            smallvec![SmolStr::from("TCC"), SmolStr::from("A")],
        );
        assert_eq!(a.ref_text(), "ACG");
        assert_eq!(a.alt_texts().len(), 2);
        assert_eq!(a.rlen(), 3);
        assert_eq!(a.n_allele(), 3);
    }

    // r[verify vcf_record.alleles_rlen]
    #[test]
    fn rlen_for_each_variant() {
        assert_eq!(Alleles::reference(Base::G).rlen(), 1);
        assert_eq!(Alleles::snv(Base::A, Base::T).unwrap().rlen(), 1);
        assert_eq!(Alleles::insertion(Base::A, &[Base::C, Base::G]).unwrap().rlen(), 1);
        assert_eq!(Alleles::deletion(Base::A, &[Base::C, Base::G]).unwrap().rlen(), 3);
        assert_eq!(Alleles::complex(SmolStr::from("ACGT"), smallvec![]).rlen(), 4);
    }

    // r[verify vcf_record.alleles_serialization]
    #[test]
    fn write_into_matches_text_methods() {
        // Verify the zero-alloc write_*_into methods produce identical output to ref_text/alt_texts
        let cases: Vec<Alleles> = vec![
            Alleles::reference(Base::A),
            Alleles::snv(Base::C, Base::T).unwrap(),
            Alleles::snv_multi(Base::G, &[Base::A, Base::T]).unwrap(),
            Alleles::insertion(Base::A, &[Base::C, Base::G, Base::T]).unwrap(),
            Alleles::deletion(Base::T, &[Base::A, Base::C]).unwrap(),
            Alleles::complex(
                SmolStr::from("ACG"),
                smallvec![SmolStr::from("T"), SmolStr::from("GG")],
            ),
        ];
        for alleles in &cases {
            let mut ref_buf = Vec::new();
            alleles.write_ref_into(&mut ref_buf);
            assert_eq!(
                String::from_utf8(ref_buf).unwrap(),
                alleles.ref_text().as_str(),
                "ref mismatch for {alleles:?}"
            );

            let mut alt_buf = Vec::new();
            let n = alleles.write_alts_into(&mut alt_buf);
            let alt_texts = alleles.alt_texts();
            assert_eq!(n, alt_texts.len(), "alt count mismatch for {alleles:?}");
            if !alt_texts.is_empty() {
                let expected: String =
                    alt_texts.iter().map(|s| s.as_str()).collect::<Vec<_>>().join(",");
                assert_eq!(
                    String::from_utf8(alt_buf).unwrap(),
                    expected,
                    "alt text mismatch for {alleles:?}"
                );
            }
        }
    }
}

#[cfg(test)]
mod proptests {
    use super::*;
    use proptest::prelude::*;

    fn arb_base() -> impl Strategy<Value = Base> {
        prop_oneof![Just(Base::A), Just(Base::C), Just(Base::G), Just(Base::T),]
    }

    fn arb_alleles() -> impl Strategy<Value = Alleles> {
        prop_oneof![
            // Reference
            arb_base().prop_map(Alleles::reference),
            // SNV — pick ref and alt that differ
            (arb_base(), arb_base())
                .prop_filter("ref != alt", |(r, a)| r != a)
                .prop_map(|(r, a)| Alleles::snv(r, a).unwrap()),
            // Insertion — 1-8 inserted bases
            (arb_base(), proptest::collection::vec(arb_base(), 1..8))
                .prop_map(|(anchor, ins)| Alleles::insertion(anchor, &ins).unwrap()),
            // Deletion — 1-8 deleted bases
            (arb_base(), proptest::collection::vec(arb_base(), 1..8))
                .prop_map(|(anchor, del)| Alleles::deletion(anchor, &del).unwrap()),
        ]
    }

    // r[verify vcf_record.alleles_serialization]
    proptest! {
        #[test]
        fn write_into_consistent_with_text_methods(alleles in arb_alleles()) {
            let mut ref_buf = Vec::new();
            alleles.write_ref_into(&mut ref_buf);
            let ref_str = String::from_utf8(ref_buf).unwrap();
            let ref_text = alleles.ref_text();
            prop_assert_eq!(ref_str.as_str(), ref_text.as_str());

            let mut alt_buf = Vec::new();
            let n = alleles.write_alts_into(&mut alt_buf);
            let alt_texts = alleles.alt_texts();
            prop_assert_eq!(n, alt_texts.len());
            if !alt_texts.is_empty() {
                let expected: String = alt_texts.iter().map(|s| s.as_str()).collect::<Vec<_>>().join(",");
                let alt_str = String::from_utf8(alt_buf).unwrap();
                prop_assert_eq!(alt_str, expected);
            }
        }

        // r[verify vcf_record.alleles_rlen]
        #[test]
        fn rlen_matches_ref_text_length(alleles in arb_alleles()) {
            let ref_text = alleles.ref_text();
            prop_assert_eq!(alleles.rlen(), ref_text.len());
        }

        // r[verify vcf_record.allele_count]
        #[test]
        fn n_allele_equals_one_plus_alts(alleles in arb_alleles()) {
            let alts = alleles.alt_texts();
            let expected = 1usize.saturating_add(alts.len());
            prop_assert_eq!(alleles.n_allele(), expected);
        }
    }
}
