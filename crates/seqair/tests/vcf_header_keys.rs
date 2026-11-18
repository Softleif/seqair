//! Resolving typed keys, filters, contigs and samples from a header, checked
//! against its declarations and branded with its `HeaderId`.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use seqair::vcf::{
    Arr, ContigDef, Flag, FormatGt, FormatInt, Gt, InfoFieldDef, InfoFlag, InfoFloat, InfoInt,
    InfoInts, InfoString, KeyError, KeyKind, Number, OptArr, Scalar, Str, ValueType, VcfHeader,
};

const TEXT: &str = concat!(
    "##fileformat=VCFv4.2\n",
    "##FILTER=<ID=q10,Description=\"Quality below 10\">\n",
    "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n",
    "##INFO=<ID=BQ,Number=1,Type=Float,Description=\"Base quality\">\n",
    "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n",
    "##INFO=<ID=CPG,Number=0,Type=Flag,Description=\"CpG\">\n",
    "##INFO=<ID=DB,Number=1,Type=Flag,Description=\"Flag with a count\">\n",
    "##INFO=<ID=CH,Number=1,Type=Character,Description=\"A character\">\n",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
    "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n",
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Sample depth\">\n",
    "##contig=<ID=chr1,length=100>\n##contig=<ID=chr2>\n",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\tHG002\n",
);

fn header() -> VcfHeader {
    VcfHeader::parse(TEXT).unwrap()
}

// r[verify vcf_values.keys+2]
// r[verify vcf_values.field_kind]
#[test]
fn keys_resolve_when_the_marker_fits_the_declaration() {
    let h = header();
    let dp: InfoInt = h.info_key("DP").unwrap();
    assert_eq!((dp.id().name(), dp.id().dict_idx()), ("DP", 2));
    let _: InfoFloat = h.info_key("BQ").unwrap();
    let _: InfoInts = h.info_key("AC").unwrap();
    let _ = h.info_key::<OptArr<i32>>("AC").unwrap();
    let _ = h.info_key::<Arr<i32>>("DP").unwrap();
    let _: InfoFlag = h.info_key("CPG").unwrap();
    let _: InfoFlag = h.info_key("DB").unwrap(); // a Flag reads whatever its Number says
    let _: InfoString = h.info_key("CH").unwrap(); // Str reads Character
    let gt: FormatGt = h.format_key("GT").unwrap();
    let gq: FormatInt = h.format_key("GQ").unwrap();
    let fdp: FormatInt = h.format_key("DP").unwrap();
    assert_eq!(fdp.id().dict_idx(), dp.id().dict_idx(), "INFO and FORMAT DP share an index");
    assert_eq!((gt.id().dict_idx(), gq.id().dict_idx()), (8, 9));
}

// r[verify vcf_values.keys+2]
#[test]
fn a_key_that_does_not_fit_names_the_reason() {
    let h = header();
    let err = h.info_key::<Scalar<i32>>("BQ").unwrap_err();
    assert_eq!(
        err,
        KeyError::TypeMismatch {
            name: "BQ".into(),
            declared: ValueType::Float,
            requested: ValueType::Integer
        }
    );
    let err = h.info_key::<Flag>("DP").unwrap_err();
    assert!(matches!(err, KeyError::TypeMismatch { requested: ValueType::Flag, .. }), "{err:?}");
    let err = h.info_key::<Str>("DP").unwrap_err();
    assert!(matches!(err, KeyError::TypeMismatch { declared: ValueType::Integer, .. }), "{err:?}");
    let err = h.info_key::<Scalar<i32>>("AC").unwrap_err();
    assert_eq!(
        err,
        KeyError::NumberMismatch { name: "AC".into(), declared: Number::AlternateBases }
    );
    assert_eq!(h.format_key::<Gt>("GQ").unwrap_err(), KeyError::NotGenotype { name: "GQ".into() });
    let err = h.info_key::<Scalar<i32>>("NOPE").unwrap_err();
    assert_eq!(err, KeyError::UnknownField { name: "NOPE".into(), kind: KeyKind::Info });
    assert_eq!(
        h.format_key::<Scalar<f32>>("BQ").unwrap_err(),
        KeyError::UnknownField { name: "BQ".into(), kind: KeyKind::Format },
        "an INFO id is not a FORMAT field"
    );
    assert_eq!(h.filter("q20").unwrap_err(), KeyError::UnknownFilter { name: "q20".into() });
    assert_eq!(
        h.filter("DP").unwrap_err(),
        KeyError::UnknownFilter { name: "DP".into() },
        "an INFO id is not a FILTER"
    );
    assert_eq!(h.contig("chrM").unwrap_err(), KeyError::UnknownContig { name: "chrM".into() });
    assert_eq!(
        h.sample("NA24385").unwrap_err(),
        KeyError::UnknownSample { name: "NA24385".into() }
    );
}

// r[verify vcf_values.keys+2]
// r[verify vcf_fields.sample_idx]
#[test]
fn filters_contigs_and_samples_resolve() {
    let h = header();
    assert_eq!(h.filter("PASS").unwrap().dict_idx(), 0);
    assert_eq!(h.filter("q10").unwrap().dict_idx(), 1);
    let chr2 = h.contig("chr2").unwrap();
    assert_eq!((chr2.name().as_str(), chr2.index()), ("chr2", 1));
    let hg002 = h.sample("HG002").unwrap();
    assert_eq!(hg002.index(), 1);
    let all: Vec<_> = h.sample_indices().collect();
    assert_eq!(all.len(), 2);
    assert_eq!(all, [h.sample("NA12878").unwrap(), hg002]);
}

// r[verify vcf_header.brand+2]
// r[verify vcf_values.header_brand+2]
// r[verify vcf_fields.sample_idx]
/// Two parses of the same text are two headers: their keys, contigs and
/// samples differ, though name and index agree; a clone's are the same.
#[test]
fn keys_from_two_headers_differ() {
    let (a, b) = (header(), header());
    assert_ne!(a.id(), b.id());
    assert_eq!(a.clone().id(), a.id());
    let (ka, kb): (InfoInt, InfoInt) = (a.info_key("DP").unwrap(), b.info_key("DP").unwrap());
    assert_eq!(ka.id().dict_idx(), kb.id().dict_idx());
    assert_ne!(ka, kb);
    assert_eq!(ka, a.clone().info_key("DP").unwrap());
    assert_ne!(a.contig("chr1").unwrap(), b.contig("chr1").unwrap());
    assert_ne!(a.filter("q10").unwrap(), b.filter("q10").unwrap());
    assert_ne!(a.sample("HG002").unwrap(), b.sample("HG002").unwrap());
}

// r[verify vcf_values.keys+2]
/// A key the writer registered is the key the header resolves.
#[test]
fn a_registered_key_is_the_resolved_key() {
    let mut builder = VcfHeader::builder();
    let chr1 = builder.register_contig("chr1", ContigDef::new(Some(10))).unwrap();
    let mut builder = builder.infos();
    let dp = builder.register_info(&InfoFieldDef::<Scalar<i32>>::new("DP", "")).unwrap();
    let header = builder.build().unwrap();
    assert_eq!(header.info_key::<Scalar<i32>>("DP").unwrap(), dp);
    assert_eq!(header.contig("chr1").unwrap(), chr1);
}
