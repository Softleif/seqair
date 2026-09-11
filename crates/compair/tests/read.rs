//! What `Read` and `Haplotype` accept, and what they say when they do not.

use compair::{Base, BaseQuality, Error, Haplotype, Read, Strand};
use proptest::prelude::*;

fn q(byte: u8) -> BaseQuality {
    BaseQuality::from_byte(byte)
}

/// Every quality track has to be as long as the read, and the error names the
/// track and both lengths so that a caller can see which one it got wrong.
#[test]
fn a_quality_track_of_the_wrong_length_names_itself() {
    let bases = vec![Base::A, Base::C, Base::G];
    let three = [q(30); 3];
    let two = [q(30); 2];
    let four = [q(30); 4];
    for (field, base_quals, insertion, deletion, gap, actual) in [
        ("base_quals", &two[..], &three[..], &three[..], &three[..], 2),
        ("insertion_quals", &three[..], &four[..], &three[..], &three[..], 4),
        ("deletion_quals", &three[..], &three[..], &two[..], &three[..], 2),
        ("gap_quals", &three[..], &three[..], &three[..], &[][..], 0),
    ] {
        let got = Read::new(bases.clone(), base_quals, insertion, deletion, gap, Strand::OT);
        assert_eq!(
            got.map(|_| ()),
            Err(Error::ReadLengthMismatch { field, expected: 3, actual }),
            "{field}"
        );
    }
    // `uniform` fills the three indel tracks itself, so only the base
    // qualities can be wrong.
    let got = Read::uniform(bases, &two, q(45), q(45), q(10), Strand::OT);
    assert_eq!(
        got.map(|_| ()),
        Err(Error::ReadLengthMismatch { field: "base_quals", expected: 3, actual: 2 })
    );
}

/// `BaseQuality::UNAVAILABLE` has no Phred value, so no track may carry it,
/// and the error says where it was.
#[test]
fn an_unavailable_quality_is_rejected_with_its_position() {
    let bases = vec![Base::A, Base::C, Base::G];
    let good = [q(30); 3];
    let mut bad = good;
    if let Some(slot) = bad.get_mut(1) {
        *slot = BaseQuality::UNAVAILABLE;
    }
    for (field, base_quals, insertion, deletion, gap) in [
        ("base_quals", &bad, &good, &good, &good),
        ("insertion_quals", &good, &bad, &good, &good),
        ("deletion_quals", &good, &good, &bad, &good),
        ("gap_quals", &good, &good, &good, &bad),
    ] {
        let got = Read::new(bases.clone(), base_quals, insertion, deletion, gap, Strand::OT);
        assert_eq!(got.map(|_| ()), Err(Error::MissingQuality { field, index: 1 }), "{field}");
    }
    let got = Read::uniform(bases, &good, BaseQuality::UNAVAILABLE, q(45), q(10), Strand::OT);
    assert_eq!(got.map(|_| ()), Err(Error::MissingQuality { field: "insertion_quals", index: 0 }));
}

/// The emission decides which strand can see a conversion, so a read has to
/// say which one it is on.
#[test]
fn an_unknown_strand_is_rejected() {
    let got = Read::uniform(vec![Base::A], &[q(30)], q(45), q(45), q(10), Strand::Unknown);
    assert_eq!(got.map(|_| ()), Err(Error::UnknownStrand));
    // Strand is checked first: it is the one thing no other argument can fix.
    let got = Read::new(vec![Base::A], &[], &[], &[], &[], Strand::Unknown);
    assert_eq!(got.map(|_| ()), Err(Error::UnknownStrand));
}

/// A soft-masked reference base is still that base: lowercase is not `N`.
#[test]
fn from_ascii_reads_lowercase_as_the_same_base() {
    let haplotype = Haplotype::from_ascii(b"acgtACGTnN.-RYK");
    assert_eq!(
        haplotype.bases(),
        &[
            Base::A,
            Base::C,
            Base::G,
            Base::T,
            Base::A,
            Base::C,
            Base::G,
            Base::T,
            Base::Unknown,
            Base::Unknown,
            Base::Unknown,
            Base::Unknown,
            Base::Unknown,
            Base::Unknown,
            Base::Unknown,
        ]
    );
    // And so a lowercase `cg` is a CpG.
    assert_eq!(
        Haplotype::from_ascii(b"acgt").site(1).map(|site| site.cpg),
        Haplotype::from_ascii(b"ACGT").site(1).map(|site| site.cpg),
    );
}

proptest! {
    /// Any byte string builds a haplotype, and the case of a letter never
    /// changes the base.
    #[test]
    fn from_ascii_is_case_insensitive(sequence in proptest::collection::vec(any::<u8>(), 0..64)) {
        let lower = Haplotype::from_ascii(&sequence.to_ascii_lowercase());
        let upper = Haplotype::from_ascii(&sequence.to_ascii_uppercase());
        prop_assert_eq!(lower.len(), sequence.len());
        prop_assert_eq!(lower, upper);
    }

    /// Every observation of a well-formed read carries its own error
    /// probability, and the transitions a `Read` precomputes are what the
    /// qualities say -- checked here with `powf` as the oracle, since a
    /// `Read` never calls it again after construction.
    #[test]
    fn a_read_precomputes_what_its_qualities_say(
        bases in proptest::collection::vec(0u8..=4, 1..40),
        quals in proptest::collection::vec(0u8..=93, 40),
    ) {
        let bases: Vec<Base> = bases.iter().map(|code| *Base::KNOWN.get(usize::from(*code)).unwrap_or(&Base::Unknown)).collect();
        let n = bases.len();
        let track = |offset: usize| -> Vec<BaseQuality> {
            (0..n).map(|index| q(*quals.get((index + offset) % quals.len()).unwrap_or(&30))).collect()
        };
        let base_quals = track(0);
        let read = Read::new(bases.clone(), &base_quals, &track(1), &track(2), &track(3), Strand::OB)
            .map_err(|_| TestCaseError::reject("valid by construction"))?;
        prop_assert_eq!(read.len(), n);
        prop_assert_eq!(read.strand(), Strand::OB);
        for (index, (base, qual)) in bases.iter().zip(&base_quals).enumerate() {
            let observation = read.observation(index).ok_or(TestCaseError::reject("in range"))?;
            prop_assert_eq!(observation.base, *base);
            prop_assert_eq!(observation.index.get(), u32::try_from(index).unwrap_or(0));
            let want = 10f64.powf(-f64::from(qual.as_byte()) / 10.0);
            prop_assert_eq!(observation.error_probability.to_bits(), want.to_bits());
        }
        prop_assert!(read.observation(n).is_none());
    }
}
