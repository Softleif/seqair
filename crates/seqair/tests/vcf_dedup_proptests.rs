//! Property-based tests for in-place duplicate-field overwriting.
//!
//! The encoder lets a field be written more than once per record and overwrites
//! the previous value **in place** (htslib `bcf_update_*` semantics): the field
//! keeps its original column position and takes its last-written value.
//!
//! The oracle here is deliberately independent of the dedup machinery: encoding
//! a sequence of writes that contains duplicates MUST be **byte-identical** to
//! encoding the de-duplicated sequence — each field once, in first-occurrence
//! order, with its last value — which drives only the no-duplicate code path.
//! Any divergence in value, ordering, separators, or BCF framing fails the test.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "test code"
)]

use std::sync::Arc;

use proptest::prelude::*;
use seqair::vcf::record_encoder::{Arr, Scalar, Str};
use seqair::vcf::{
    Alleles, ContigDef, ContigId, FormatFloats, FormatInt, InfoFieldDef, InfoFlag, InfoFloat,
    InfoInt, InfoInts, InfoString, Number, OutputFormat, ValueType, VcfHeader, Writer,
    record_encoder::{Flag, FormatFieldDef},
};
use seqair_types::{Base, Pos1};

// ── INFO ────────────────────────────────────────────────────────────────

#[derive(Clone, Debug)]
enum InfoOp {
    Dp(i32),
    Bq(f32),
    Ad(Vec<i32>),
    Sc(String),
    Flag,
}

impl InfoOp {
    /// Stable discriminant identifying the target field (for de-duplication).
    fn field(&self) -> u8 {
        match self {
            InfoOp::Dp(_) => 0,
            InfoOp::Bq(_) => 1,
            InfoOp::Ad(_) => 2,
            InfoOp::Sc(_) => 3,
            InfoOp::Flag => 4,
        }
    }
}

struct InfoKeys {
    dp: InfoInt,
    bq: InfoFloat,
    ad: InfoInts,
    sc: InfoString,
    flag: InfoFlag,
}

fn info_header() -> (Arc<VcfHeader>, ContigId, InfoKeys) {
    let mut b = VcfHeader::builder();
    let contig = b.register_contig("chr1", ContigDef { length: Some(1_000_000) }).unwrap();
    let mut b = b.infos();
    let dp = b
        .register_info(&InfoFieldDef::<Scalar<i32>>::new(
            "DP",
            Number::Count(1),
            ValueType::Integer,
            "d",
        ))
        .unwrap();
    let bq = b
        .register_info(&InfoFieldDef::<Scalar<f32>>::new(
            "BQ",
            Number::Count(1),
            ValueType::Float,
            "q",
        ))
        .unwrap();
    let ad = b
        .register_info(&InfoFieldDef::<Arr<i32>>::new(
            "AD",
            Number::Unknown,
            ValueType::Integer,
            "a",
        ))
        .unwrap();
    let sc = b
        .register_info(&InfoFieldDef::<Str>::new("SC", Number::Count(1), ValueType::String, "s"))
        .unwrap();
    let flag = b
        .register_info(&InfoFieldDef::<Flag>::new("FL", Number::Count(0), ValueType::Flag, "f"))
        .unwrap();
    let header = Arc::new(b.build().unwrap());
    (header, contig, InfoKeys { dp, bq, ad, sc, flag })
}

/// Reduce a sequence to "each field once, first-occurrence order, last value".
fn dedup<T: Clone>(ops: &[T], field: impl Fn(&T) -> u8) -> Vec<T> {
    let mut order: Vec<u8> = Vec::new();
    let mut last: Vec<(u8, T)> = Vec::new();
    for op in ops {
        let f = field(op);
        if !order.contains(&f) {
            order.push(f);
        }
        if let Some(slot) = last.iter_mut().find(|(k, _)| *k == f) {
            slot.1 = op.clone();
        } else {
            last.push((f, op.clone()));
        }
    }
    order.iter().map(|f| last.iter().find(|(k, _)| k == f).unwrap().1.clone()).collect()
}

fn encode_info(
    format: OutputFormat,
    keys: &InfoKeys,
    contig: &ContigId,
    header: &Arc<VcfHeader>,
    ops: &[InfoOp],
) -> Vec<u8> {
    let mut buf = Vec::new();
    let writer = Writer::new(&mut buf, format);
    let mut writer = writer.write_header(header).unwrap();
    let mut enc = writer
        .begin_record(
            contig,
            Pos1::new(42).unwrap(),
            &Alleles::snv(Base::A, Base::T).unwrap(),
            Some(30.0),
        )
        .unwrap()
        .filter_pass();
    for op in ops {
        match op {
            InfoOp::Dp(v) => keys.dp.encode(&mut enc, *v),
            InfoOp::Bq(v) => keys.bq.encode(&mut enc, *v),
            InfoOp::Ad(v) => keys.ad.encode(&mut enc, v),
            InfoOp::Sc(v) => keys.sc.encode(&mut enc, v),
            InfoOp::Flag => keys.flag.encode(&mut enc),
        }
    }
    enc.emit().unwrap();
    writer.finish().unwrap();
    buf
}

fn arb_info_op() -> impl Strategy<Value = InfoOp> {
    prop_oneof![
        (-5000i32..5000).prop_map(InfoOp::Dp),
        (-1000.0f32..1000.0).prop_map(InfoOp::Bq),
        proptest::collection::vec(-300i32..300, 1..4).prop_map(InfoOp::Ad),
        "[A-Za-z0-9]{1,6}".prop_map(InfoOp::Sc),
        Just(InfoOp::Flag),
    ]
}

// ── FORMAT (3 samples — exercises the per-sample VCF splice) ─────────────

#[derive(Clone, Debug)]
enum FmtOp {
    Dp([i32; 3]),
    Q([f32; 3]),
    Pl([Vec<f32>; 3]),
}

impl FmtOp {
    fn field(&self) -> u8 {
        match self {
            FmtOp::Dp(_) => 0,
            FmtOp::Q(_) => 1,
            FmtOp::Pl(_) => 2,
        }
    }
}

struct FmtKeys {
    dp: FormatInt,
    q: seqair::vcf::FormatFloat,
    pl: FormatFloats,
}

fn fmt_header() -> (Arc<VcfHeader>, ContigId, FmtKeys) {
    let mut b = VcfHeader::builder();
    let contig = b.register_contig("chr1", ContigDef { length: Some(1_000_000) }).unwrap();
    let mut b = b.formats();
    let dp = b
        .register_format(&FormatFieldDef::<Scalar<i32>>::new(
            "DP",
            Number::Count(1),
            ValueType::Integer,
            "d",
        ))
        .unwrap();
    let q = b
        .register_format(&FormatFieldDef::<Scalar<f32>>::new(
            "GQ",
            Number::Count(1),
            ValueType::Float,
            "q",
        ))
        .unwrap();
    let pl = b
        .register_format(&FormatFieldDef::<Arr<f32>>::new(
            "PL",
            Number::Unknown,
            ValueType::Float,
            "p",
        ))
        .unwrap();
    let mut b = b.samples();
    b.add_sample("s0").unwrap();
    b.add_sample("s1").unwrap();
    b.add_sample("s2").unwrap();
    let header = Arc::new(b.build().unwrap());
    (header, contig, FmtKeys { dp, q, pl })
}

fn encode_fmt(
    format: OutputFormat,
    keys: &FmtKeys,
    contig: &ContigId,
    header: &Arc<VcfHeader>,
    ops: &[FmtOp],
) -> Vec<u8> {
    let mut buf = Vec::new();
    let writer = Writer::new(&mut buf, format);
    let mut writer = writer.write_header(header).unwrap();
    let enc = writer
        .begin_record(
            contig,
            Pos1::new(42).unwrap(),
            &Alleles::snv(Base::A, Base::T).unwrap(),
            Some(30.0),
        )
        .unwrap()
        .filter_pass();
    let mut enc = enc.begin_samples();
    for op in ops {
        match op {
            FmtOp::Dp(v) => keys.dp.encode(&mut enc, v).unwrap(),
            FmtOp::Q(v) => keys.q.encode(&mut enc, v).unwrap(),
            FmtOp::Pl(v) => {
                let per_sample: [&[f32]; 3] = [&v[0], &v[1], &v[2]];
                keys.pl.encode(&mut enc, &per_sample).unwrap();
            }
        }
    }
    enc.emit().unwrap();
    writer.finish().unwrap();
    buf
}

fn arb_fmt_op() -> impl Strategy<Value = FmtOp> {
    prop_oneof![
        proptest::array::uniform3(-5000i32..5000).prop_map(FmtOp::Dp),
        proptest::array::uniform3(-1000.0f32..1000.0).prop_map(FmtOp::Q),
        proptest::array::uniform3(proptest::collection::vec(-100.0f32..100.0, 1..4))
            .prop_map(FmtOp::Pl),
    ]
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 400, ..ProptestConfig::default() })]

    /// INFO: a duplicate-laden write sequence encodes byte-for-byte identically
    /// to its de-duplicated form, for both BCF and VCF text.
    #[test]
    fn info_in_place_overwrite_matches_dedup_oracle(
        ops in proptest::collection::vec(arb_info_op(), 1..14)
    ) {
        let (header, contig, keys) = info_header();
        let oracle = dedup(&ops, InfoOp::field);
        for format in [OutputFormat::Bcf, OutputFormat::Vcf] {
            let with_dups = encode_info(format, &keys, &contig, &header, &ops);
            let once = encode_info(format, &keys, &contig, &header, &oracle);
            prop_assert_eq!(with_dups, once, "format={:?} ops={:?}", format, ops);
        }
    }

    /// FORMAT (3 samples): same property, exercising the per-sample colon splice.
    #[test]
    fn format_in_place_overwrite_matches_dedup_oracle(
        ops in proptest::collection::vec(arb_fmt_op(), 1..12)
    ) {
        let (header, contig, keys) = fmt_header();
        let oracle = dedup(&ops, FmtOp::field);
        for format in [OutputFormat::Bcf, OutputFormat::Vcf] {
            let with_dups = encode_fmt(format, &keys, &contig, &header, &ops);
            let once = encode_fmt(format, &keys, &contig, &header, &oracle);
            prop_assert_eq!(with_dups, once, "format={:?} ops={:?}", format, ops);
        }
    }
}
