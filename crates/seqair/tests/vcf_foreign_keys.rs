//! Every encode entry point refuses a handle resolved from another header
//! (`begin_record`'s contig check is tested in `vcf::header`).
//!
//! The two headers here are built the same way, so every foreign key has the
//! same dictionary index as its counterpart: without the check, each call below
//! would write a valid-looking field. Each test makes one foreign call, expects
//! `ForeignKey` naming the field, and checks the emitted record is byte for
//! byte the one written without that call, in BCF and in VCF text.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, reason = "test code")]

use seqair::vcf::{
    Alleles, Arr, ContigDef, ContigId, FilterFieldDef, FilterId, Filtered, Flag, FormatFieldDef,
    FormatFloat, FormatFloatOpts, FormatFloats, FormatGt, FormatInt, FormatInts, FormatString,
    Genotype, Gt, InfoFieldDef, InfoFlag, InfoFloat, InfoFloats, InfoInt, InfoIntOpts, InfoInts,
    InfoString, Number, OptArr, OutputFormat, RecordEncoder, Scalar, Str, VcfError, VcfHeader,
    WithSamples, Writer,
};
use seqair_types::{Base, Pos1};

struct Keys {
    header: VcfHeader,
    contig: ContigId,
    filter: FilterId,
    dp: InfoInt,
    bq: InfoFloat,
    ad: InfoInts,
    af: InfoFloats,
    db: InfoFlag,
    sc: InfoString,
    ao: InfoIntOpts,
    gt: FormatGt,
    gq: FormatInt,
    gf: FormatFloat,
    pl: FormatInts,
    gl: FormatFloats,
    fo: FormatFloatOpts,
    fs: FormatString,
}

fn keys() -> Keys {
    let mut b = VcfHeader::builder();
    let contig = b.register_contig("chr1", ContigDef::new(Some(1000))).unwrap();
    let mut b = b.filters();
    let filter = b.register_filter(&FilterFieldDef::new("q10", "low quality")).unwrap();
    let mut b = b.infos();
    let r = Number::ReferenceAlternateBases;
    let dp = b.register_info(&InfoFieldDef::<Scalar<i32>>::new("DP", "x")).unwrap();
    let bq = b.register_info(&InfoFieldDef::<Scalar<f32>>::new("BQ", "x")).unwrap();
    let ad = b.register_info(&InfoFieldDef::<Arr<i32>>::with_number("AD", r, "x")).unwrap();
    let af =
        b.register_info(&InfoFieldDef::<Arr<f32>>::with_number("AF", Number::AlternateBases, "x"));
    let db = b.register_info(&InfoFieldDef::<Flag>::new("DB", "x")).unwrap();
    let sc =
        b.register_info(&InfoFieldDef::<Str>::with_number("SC", Number::Count(1), "x")).unwrap();
    let ao = b.register_info(&InfoFieldDef::<OptArr<i32>>::with_number("AO", r, "x")).unwrap();
    let mut b = b.formats();
    let gt = b.register_format(&FormatFieldDef::<Gt>::new("GT", "x")).unwrap();
    let gq = b.register_format(&FormatFieldDef::<Scalar<i32>>::new("GQ", "x")).unwrap();
    let gf = b.register_format(&FormatFieldDef::<Scalar<f32>>::new("GF", "x")).unwrap();
    let pl =
        b.register_format(&FormatFieldDef::<Arr<i32>>::with_number("PL", Number::Genotypes, "x"));
    let gl =
        b.register_format(&FormatFieldDef::<Arr<f32>>::with_number("GL", Number::Genotypes, "x"));
    let fo = b.register_format(&FormatFieldDef::<OptArr<f32>>::with_number("FO", r, "x")).unwrap();
    let fs = b
        .register_format(&FormatFieldDef::<Str>::with_number("FS", Number::Count(1), "x"))
        .unwrap();
    let mut b = b.samples();
    b.add_sample("S1").unwrap();
    Keys {
        header: b.build().unwrap(),
        contig,
        filter,
        dp,
        bq,
        ad,
        af: af.unwrap(),
        db,
        sc,
        ao,
        gt,
        gq,
        gf,
        pl: pl.unwrap(),
        gl: gl.unwrap(),
        fo,
        fs,
    }
}

type InfoCall = fn(&Keys, &mut RecordEncoder<'_, Filtered>) -> Result<(), VcfError>;
type FormatCall = fn(&Keys, &mut RecordEncoder<'_, WithSamples>) -> Result<(), VcfError>;

enum Call {
    None,
    Info(InfoCall),
    Format(FormatCall),
}

const FORMATS: [OutputFormat; 2] = [OutputFormat::Bcf, OutputFormat::Vcf];

/// One record through `ours`, with `call` made with `other`'s keys; the call
/// must fail with `ForeignKey { name }`.
fn record(format: OutputFormat, ours: &Keys, other: &Keys, call: &Call, name: &str) -> Vec<u8> {
    let foreign = |result: Result<(), VcfError>| match result {
        Err(VcfError::ForeignKey { name: got }) => assert_eq!(got, name),
        other => panic!("{name}: expected ForeignKey, got {other:?}"),
    };
    let mut out = Vec::new();
    let mut writer = Writer::new(&mut out, format).write_header(&ours.header).unwrap();
    let alleles = Alleles::snv(Base::A, Base::T).unwrap();
    let mut enc = writer
        .begin_record(&ours.contig, Pos1::new(5).unwrap(), &alleles, None)
        .unwrap()
        .filter_pass();
    ours.dp.encode(&mut enc, 7).unwrap();
    if let Call::Info(f) = call {
        foreign(f(other, &mut enc));
    }
    let mut enc = enc.begin_samples();
    ours.gt.encode(&mut enc, &[Genotype::unphased(0, 1)]).unwrap();
    if let Call::Format(f) = call {
        foreign(f(other, &mut enc));
    }
    ours.gq.encode(&mut enc, &[30]).unwrap();
    enc.emit().unwrap();
    writer.finish().unwrap();
    out
}

fn assert_refused(name: &str, call: Call) {
    let (ours, other) = (keys(), keys());
    for format in FORMATS {
        let written = record(format, &ours, &other, &call, name);
        let clean = record(format, &ours, &other, &Call::None, name);
        assert_eq!(written, clean, "{name} ({format:?}): a refused call writes nothing");
    }
}

// r[verify record_encoder.foreign_keys]
#[test]
fn info_int_refuses_a_foreign_key() {
    assert_refused("DP", Call::Info(|k, e| k.dp.encode(e, 1)));
}

// r[verify record_encoder.foreign_keys]
#[test]
fn info_float_refuses_a_foreign_key() {
    assert_refused("BQ", Call::Info(|k, e| k.bq.encode(e, 1.5)));
}

// r[verify record_encoder.foreign_keys]
#[test]
fn info_ints_refuses_a_foreign_key() {
    assert_refused("AD", Call::Info(|k, e| k.ad.encode(e, &[1, 2])));
}

// r[verify record_encoder.foreign_keys]
#[test]
fn info_floats_refuses_a_foreign_key() {
    assert_refused("AF", Call::Info(|k, e| k.af.encode(e, &[0.5])));
}

// r[verify record_encoder.foreign_keys]
#[test]
fn info_flag_refuses_a_foreign_key() {
    assert_refused("DB", Call::Info(|k, e| k.db.encode(e)));
}

// r[verify record_encoder.foreign_keys]
#[test]
fn info_string_refuses_a_foreign_key() {
    assert_refused("SC", Call::Info(|k, e| k.sc.encode(e, "ACGT")));
}

// r[verify record_encoder.foreign_keys]
#[test]
fn info_int_opts_refuses_a_foreign_key() {
    assert_refused("AO", Call::Info(|k, e| k.ao.encode(e, &[Some(1), None])));
}

// r[verify record_encoder.foreign_keys]
#[test]
fn format_gt_refuses_a_foreign_key() {
    // A rewrite of the GT already present: refused, so the first one stays.
    assert_refused("GT", Call::Format(|k, e| k.gt.encode(e, &[Genotype::unphased(1, 1)])));
}

// r[verify record_encoder.foreign_keys]
#[test]
fn format_int_refuses_a_foreign_key() {
    assert_refused("GQ", Call::Format(|k, e| k.gq.encode(e, &[99])));
}

// r[verify record_encoder.foreign_keys]
#[test]
fn format_float_refuses_a_foreign_key() {
    assert_refused("GF", Call::Format(|k, e| k.gf.encode(e, &[0.25])));
}

// r[verify record_encoder.foreign_keys]
#[test]
fn format_ints_refuses_a_foreign_key() {
    assert_refused("PL", Call::Format(|k, e| k.pl.encode(e, &[&[0, 10, 20]])));
}

// r[verify record_encoder.foreign_keys]
#[test]
fn format_floats_refuses_a_foreign_key() {
    assert_refused("GL", Call::Format(|k, e| k.gl.encode(e, &[&[0.0, -1.0, -2.0]])));
}

// r[verify record_encoder.foreign_keys]
#[test]
fn format_float_opts_refuses_a_foreign_key() {
    assert_refused("FO", Call::Format(|k, e| k.fo.encode_single_sample(e, &[Some(0.5), None])));
}

// r[verify record_encoder.foreign_keys]
#[test]
fn format_string_refuses_a_foreign_key() {
    assert_refused("FS", Call::Format(|k, e| k.fs.encode(e, &["x"])));
}

// r[verify record_encoder.foreign_keys]
#[test]
fn filter_fail_refuses_a_foreign_filter_and_takes_pass_from_any_header() {
    let (ours, other) = (keys(), keys());
    let alleles = Alleles::snv(Base::A, Base::T).unwrap();
    for format in FORMATS {
        let mut out = Vec::new();
        let mut writer = Writer::new(&mut out, format).write_header(&ours.header).unwrap();
        let pos = Pos1::new(5).unwrap();
        let enc = writer.begin_record(&ours.contig, pos, &alleles, None).unwrap();
        match enc.filter_fail([&ours.filter, &other.filter]) {
            Err(VcfError::ForeignKey { name }) => assert_eq!(name, "q10"),
            Err(e) => panic!("{format:?}: expected ForeignKey, got {e:?}"),
            Ok(_) => panic!("{format:?}: a foreign filter was accepted"),
        }
        let enc = writer.begin_record(&ours.contig, pos, &alleles, None).unwrap();
        enc.filter_fail([&FilterId::PASS, &ours.filter]).unwrap().emit().unwrap();
        writer.finish().unwrap();
        // Only the second record was emitted.
        let mut clean = Vec::new();
        let mut writer = Writer::new(&mut clean, format).write_header(&ours.header).unwrap();
        let enc = writer.begin_record(&ours.contig, pos, &alleles, None).unwrap();
        enc.filter_fail([&FilterId::PASS, &ours.filter]).unwrap().emit().unwrap();
        writer.finish().unwrap();
        assert_eq!(out, clean, "{format:?}");
    }
}
