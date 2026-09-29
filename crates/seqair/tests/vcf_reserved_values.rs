//! The values BCF will not carry as data, and the ones it will.
//!
//! BCF2 reserves the eight most negative values of every integer width and the
//! float bit patterns `0x7F800001..=0x7F800007`. The first two of each band are
//! MISSING and `END_OF_VECTOR`. A reader that meets `END_OF_VECTOR` in the middle
//! of a row stops there, so writing one as data silently shortens the row —
//! see `reserved_end_of_vector_written_as_data_truncates_the_row`, which
//! demonstrates the loss on a file bcftools reads, by patching the byte the
//! encoder now refuses to write.
//!
//! The rule under test is the width one: a value is reserved only for the
//! width it is actually emitted at. `smallest_int_type` never selects int8 or
//! int16 for a value in those widths' bands, so `-128` is ordinary `i32` data
//! (written as an int16 `0xFF80`) while `i32::MIN` is not.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "test code"
)]

use hegel::prelude::*;
use seqair::vcf::record_encoder::{
    Arr, FormatEncoder, FormatFieldDef, FormatFloats, FormatGt, FormatInt, FormatInts, InfoEncoder,
    InfoFieldDef, InfoFloat, InfoFloats, InfoInt, InfoIntOpts, InfoInts, OptArr, Scalar,
};
use seqair::vcf::{
    Alleles, ContigDef, ContigId, Genotype, Number, OutputFormat, ReservedMarker, ValueType,
    VcfEncodeError, VcfError, VcfHeader, Writer,
};
use seqair_types::{Base, Pos1};
use std::process::Command;
use std::sync::Arc;

#[path = "../src/pinned.rs"]
mod pinned;

/// First reserved float bit pattern (MISSING); the band runs to `+ 6`.
const FLOAT_RESERVED_FIRST: u32 = 0x7F80_0001;
const FLOAT_RESERVED_LEN: u32 = 7;
/// The reserved integer band is the eight most negative values of the width.
const INT_RESERVED_LEN: i32 = 8;

struct Setup {
    header: Arc<VcfHeader>,
    contig: ContigId,
    info_int: InfoInt,
    info_ints: InfoInts,
    info_int_opts: InfoIntOpts,
    info_float: InfoFloat,
    info_floats: InfoFloats,
    gt: FormatGt,
    fmt_int: FormatInt,
    fmt_ints: FormatInts,
    fmt_float: seqair::vcf::FormatFloat,
    fmt_floats: FormatFloats,
}

fn setup() -> Setup {
    let mut builder = VcfHeader::builder();
    let contig = builder.register_contig("chr1", ContigDef { length: Some(100_000) }).unwrap();
    let mut builder = builder.infos();
    let info_int = builder
        .register_info(&InfoFieldDef::<Scalar<i32>>::new(
            "DP",
            Number::Count(1),
            ValueType::Integer,
            "Scalar integer",
        ))
        .unwrap();
    let info_ints = builder
        .register_info(&InfoFieldDef::<Arr<i32>>::new(
            "AD",
            Number::Unknown,
            ValueType::Integer,
            "Integer list",
        ))
        .unwrap();
    let info_int_opts = builder
        .register_info(&InfoFieldDef::<OptArr<i32>>::new(
            "OD",
            Number::Unknown,
            ValueType::Integer,
            "Integer list with missing elements",
        ))
        .unwrap();
    let info_float = builder
        .register_info(&InfoFieldDef::<Scalar<f32>>::new(
            "MQ",
            Number::Count(1),
            ValueType::Float,
            "Scalar float",
        ))
        .unwrap();
    let info_floats = builder
        .register_info(&InfoFieldDef::<Arr<f32>>::new(
            "AF",
            Number::Unknown,
            ValueType::Float,
            "Float list",
        ))
        .unwrap();
    let mut builder = builder.formats();
    let gt = builder
        .register_format(&FormatFieldDef::new("GT", Number::Count(1), ValueType::String, "GT"))
        .unwrap();
    let fmt_int = builder
        .register_format(&FormatFieldDef::<Scalar<i32>>::new(
            "FD",
            Number::Count(1),
            ValueType::Integer,
            "Per-sample integer",
        ))
        .unwrap();
    let fmt_ints = builder
        .register_format(&FormatFieldDef::<Arr<i32>>::new(
            "FA",
            Number::Unknown,
            ValueType::Integer,
            "Per-sample integer list",
        ))
        .unwrap();
    let fmt_float = builder
        .register_format(&FormatFieldDef::<Scalar<f32>>::new(
            "FQ",
            Number::Count(1),
            ValueType::Float,
            "Per-sample float",
        ))
        .unwrap();
    let fmt_floats = builder
        .register_format(&FormatFieldDef::<Arr<f32>>::new(
            "FL",
            Number::Unknown,
            ValueType::Float,
            "Per-sample float list",
        ))
        .unwrap();
    let mut builder = builder.samples();
    builder.add_sample("s0").unwrap();
    Setup {
        header: Arc::new(builder.build().unwrap()),
        contig,
        info_int,
        info_ints,
        info_int_opts,
        info_float,
        info_floats,
        gt,
        fmt_int,
        fmt_ints,
        fmt_float,
        fmt_floats,
    }
}

/// Every entry point that takes a caller-supplied number into a BCF numeric
/// field. `format_gt` is absent on purpose: its encoding is never negative.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Entry {
    InfoInt,
    InfoInts,
    InfoIntOpts,
    FormatInt,
    FormatInts,
    InfoFloat,
    InfoFloats,
    FormatFloat,
    FormatFloats,
    Qual,
}

const INT_ENTRIES: [Entry; 5] =
    [Entry::InfoInt, Entry::InfoInts, Entry::InfoIntOpts, Entry::FormatInt, Entry::FormatInts];
const FLOAT_ENTRIES: [Entry; 5] =
    [Entry::InfoFloat, Entry::InfoFloats, Entry::FormatFloat, Entry::FormatFloats, Entry::Qual];

impl Entry {
    fn is_info(self) -> bool {
        matches!(
            self,
            Entry::InfoInt
                | Entry::InfoInts
                | Entry::InfoIntOpts
                | Entry::InfoFloat
                | Entry::InfoFloats
        )
    }

    fn is_format(self) -> bool {
        matches!(
            self,
            Entry::FormatInt | Entry::FormatInts | Entry::FormatFloat | Entry::FormatFloats
        )
    }
}

/// Encode `int` or `float` through `entry` if it is an INFO entry point.
fn encode_info(
    s: &Setup,
    enc: &mut impl InfoEncoder,
    entry: Entry,
    int: i32,
    float: f32,
) -> Result<(), VcfError> {
    match entry {
        Entry::InfoInt => s.info_int.encode(enc, int),
        Entry::InfoInts => s.info_ints.encode(enc, &[int]),
        Entry::InfoIntOpts => s.info_int_opts.encode(enc, &[Some(int), None]),
        Entry::InfoFloat => s.info_float.encode(enc, float),
        Entry::InfoFloats => s.info_floats.encode(enc, &[float]),
        _ => Ok(()),
    }
}

/// Encode `int` or `float` through `entry` if it is a FORMAT entry point.
fn encode_format(
    s: &Setup,
    enc: &mut impl FormatEncoder,
    entry: Entry,
    int: i32,
    float: f32,
) -> Result<(), VcfError> {
    match entry {
        Entry::FormatInt => s.fmt_int.encode(enc, &[int]),
        Entry::FormatInts => s.fmt_ints.encode(enc, &[&[int][..]]),
        Entry::FormatFloat => s.fmt_float.encode(enc, &[float]),
        Entry::FormatFloats => s.fmt_floats.encode(enc, &[&[float][..]]),
        _ => Ok(()),
    }
}

/// What [`write_one`] does at its entry point.
#[derive(Debug, Clone, Copy)]
enum Call {
    /// Pass the value.
    Value,
    /// Make no call. For QUAL, which every record has, this is `Value`.
    Omitted,
    /// Pass this reserved pair, which must be rejected, then carry on as if no
    /// call had been made.
    Rejected(i32, f32),
}

/// Write one record carrying `int` or `float` through `entry`, in `format`.
///
/// Returns the encoder's error rather than panicking, so the reject-direction
/// properties can inspect it.
fn write_one(
    s: &Setup,
    format: OutputFormat,
    entry: Entry,
    int: i32,
    float: f32,
    call: Call,
) -> Result<Vec<u8>, VcfError> {
    let mut out = Vec::new();
    {
        let writer = Writer::new(&mut out, format);
        let mut writer = writer.write_header(&s.header)?;
        let alleles = Alleles::snv(Base::A, Base::T)?;
        let pos = Pos1::new(100).unwrap();
        let qual = |q: f32| if entry == Entry::Qual { Some(q) } else { Some(30.0) };
        if let Call::Rejected(_, bad) = call
            && entry == Entry::Qual
        {
            assert!(writer.begin_record(&s.contig, pos, &alleles, Some(bad)).is_err());
        }
        let mut enc = writer.begin_record(&s.contig, pos, &alleles, qual(float))?.filter_pass();
        match call {
            Call::Value => encode_info(s, &mut enc, entry, int, float)?,
            Call::Omitted => {}
            Call::Rejected(bad_int, bad_float) => assert_eq!(
                encode_info(s, &mut enc, entry, bad_int, bad_float).is_err(),
                entry.is_info()
            ),
        }
        let mut enc = enc.begin_samples();
        s.gt.encode(&mut enc, &[Genotype::unphased(0, 1)])?;
        match call {
            Call::Value => encode_format(s, &mut enc, entry, int, float)?,
            Call::Omitted => {}
            Call::Rejected(bad_int, bad_float) => assert_eq!(
                encode_format(s, &mut enc, entry, bad_int, bad_float).is_err(),
                entry.is_format()
            ),
        }
        enc.emit()?;
        writer.finish()?;
    }
    Ok(out)
}

/// The `bcftools query` expression that prints what `entry` wrote.
fn query_format(entry: Entry) -> &'static str {
    match entry {
        Entry::InfoInt => "%INFO/DP\n",
        Entry::InfoInts => "%INFO/AD\n",
        Entry::InfoIntOpts => "%INFO/OD\n",
        Entry::FormatInt => "[%FD\n]",
        Entry::FormatInts => "[%FA\n]",
        Entry::InfoFloat => "%INFO/MQ\n",
        Entry::InfoFloats => "%INFO/AF\n",
        Entry::FormatFloat => "[%FQ\n]",
        Entry::FormatFloats => "[%FL\n]",
        Entry::Qual => "%QUAL\n",
    }
}

fn bcftools_query(bcf: &[u8], format: &str) -> Vec<String> {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("generated.bcf");
    std::fs::write(&path, bcf).unwrap();
    let out = Command::new("bcftools")
        .args(["query", "-f", format])
        .arg(&path)
        .output()
        .expect("bcftools not found");
    assert!(
        out.status.success(),
        "bcftools query failed: {}",
        std::str::from_utf8(&out.stderr).unwrap()
    );
    String::from_utf8(out.stdout).unwrap().lines().map(str::to_owned).collect()
}

// ── Accept direction ───────────────────────────────────────────────────

// r[verify bcf_encoder.reserved_bands]
// r[verify bcf_encoder.reserved_width]
/// Every integer the encoder does accept must survive bcftools unchanged.
///
/// The domain is the whole of `i32` below the int32 reserved band, drawn
/// directly rather than filtered, so the smallest counterexample hegel can
/// report stays small. `i32::MIN + 8` is the first accepted value and is where
/// `r[bcf_writer.smallest_int_type]`'s int32 range starts.
#[hegel::test(test_cases = crate::pinned::cases(40))]
fn accepted_integers_round_trip_through_bcftools(tc: TestCase) {
    let value = tc.draw(gs::integers::<i32>().min_value(i32::MIN + INT_RESERVED_LEN));
    let entry = tc.draw(gs::sampled_from(&INT_ENTRIES).print_as_debug());

    let s = setup();
    let bcf =
        write_one(&s, OutputFormat::Bcf, entry, value, 30.0, Call::Value).expect("accepted value");
    let lines = bcftools_query(&bcf, query_format(entry));

    let want = match entry {
        // The trailing `None` is the missing marker, which bcftools prints `.`.
        Entry::InfoIntOpts => format!("{value},."),
        _ => value.to_string(),
    };
    assert_eq!(lines, vec![want], "{entry:?} value {value}");

    tc.event(match value {
        -120..=127 => "emitted as int8",
        -32_760..=32_767 => "emitted as int16",
        _ => "emitted as int32",
    });
}

// r[verify bcf_encoder.reserved_width]
/// The int8 and int16 reserved bands are ordinary `i32` data.
///
/// This is the half of the width rule that a naive "reject anything that looks
/// like a sentinel" check gets wrong. `-128` is the int8 MISSING *pattern*, but
/// `smallest_int_type` promotes it to int16 and writes `0xFF80`, so it must be
/// accepted and must read back as `-128`.
#[hegel::test(test_cases = crate::pinned::cases(40))]
fn the_int8_and_int16_reserved_bands_are_ordinary_integers(tc: TestCase) {
    let narrow = tc.draw(gs::booleans());
    let offset = tc.draw(gs::integers::<i32>().min_value(0).max_value(INT_RESERVED_LEN - 1));
    let value = if narrow { i32::from(i8::MIN) + offset } else { i32::from(i16::MIN) + offset };
    let entry = tc.draw(gs::sampled_from(&INT_ENTRIES).print_as_debug());

    let s = setup();
    let bcf = write_one(&s, OutputFormat::Bcf, entry, value, 30.0, Call::Value)
        .unwrap_or_else(|e| panic!("{value} must be ordinary data at its emitted width: {e}"));
    let lines = bcftools_query(&bcf, query_format(entry));
    let want = match entry {
        Entry::InfoIntOpts => format!("{value},."),
        _ => value.to_string(),
    };
    assert_eq!(lines, vec![want], "{entry:?} value {value}");
}

// r[verify bcf_encoder.reserved_bands]
/// Quiet NaN and both infinities are ordinary Float values per BCF2, so a check
/// that keyed on `is_nan()` — or on the exponent alone — would wrongly reject
/// them. `0x7F800000` is `+inf`, one below the reserved band's first member.
#[test]
fn quiet_nan_and_the_infinities_are_ordinary_floats() {
    let s = setup();
    for bits in [0x7FC0_0000u32, 0x7F80_0000, 0xFF80_0000] {
        let value = f32::from_bits(bits);
        for entry in FLOAT_ENTRIES {
            write_one(&s, OutputFormat::Bcf, entry, 1, value, Call::Value)
                .unwrap_or_else(|e| panic!("{bits:#010X} via {entry:?} must be accepted: {e}"));
        }
    }
}

// ── Reject direction ───────────────────────────────────────────────────

// r[verify bcf_encoder.reserved_rejected]
// r[verify record_encoder.reserved_rejection]
/// Every value of the int32 reserved band, through every integer entry point,
/// in every output format, must come back as `ReservedIntValue` naming that
/// band member — not as a panic, not as a silently written record.
#[hegel::test]
fn every_reserved_integer_is_rejected_at_every_entry_point(tc: TestCase) {
    let offset = tc.draw(gs::integers::<i32>().min_value(0).max_value(INT_RESERVED_LEN - 1));
    let entry = tc.draw(gs::sampled_from(&INT_ENTRIES).print_as_debug());
    let format = tc.draw(
        gs::sampled_from(&[OutputFormat::Vcf, OutputFormat::VcfGz, OutputFormat::Bcf])
            .print_as_debug(),
    );
    let value = i32::MIN + offset;

    let s = setup();
    let err = write_one(&s, format, entry, value, 30.0, Call::Value)
        .expect_err("a reserved integer must not be written");

    let want_marker = match offset {
        0 => ReservedMarker::Missing,
        1 => ReservedMarker::EndOfVector,
        _ => ReservedMarker::FutureUse,
    };
    match err {
        VcfError::Encode(VcfEncodeError::ReservedIntValue { value: got, marker, .. }) => {
            assert_eq!(got, value);
            assert_eq!(marker, want_marker);
        }
        other => panic!("{entry:?} {value}: expected ReservedIntValue, got {other:?}"),
    }
}

// r[verify bcf_encoder.reserved_rejected]
// r[verify record_encoder.reserved_rejection]
/// The float twin: every bit pattern of `0x7F800001..=0x7F800007`, through
/// every float entry point including QUAL.
#[hegel::test]
fn every_reserved_float_is_rejected_at_every_entry_point(tc: TestCase) {
    let offset = tc.draw(gs::integers::<u32>().max_value(FLOAT_RESERVED_LEN - 1));
    let entry = tc.draw(gs::sampled_from(&FLOAT_ENTRIES).print_as_debug());
    let format = tc.draw(
        gs::sampled_from(&[OutputFormat::Vcf, OutputFormat::VcfGz, OutputFormat::Bcf])
            .print_as_debug(),
    );
    let bits = FLOAT_RESERVED_FIRST + offset;

    let s = setup();
    let err = write_one(&s, format, entry, 1, f32::from_bits(bits), Call::Value)
        .expect_err("a reserved float must not be written");

    let want_marker = match offset {
        0 => ReservedMarker::Missing,
        1 => ReservedMarker::EndOfVector,
        _ => ReservedMarker::FutureUse,
    };
    match err {
        VcfError::Encode(VcfEncodeError::ReservedFloatValue { bits: got, marker, .. }) => {
            assert_eq!(got, bits);
            assert_eq!(marker, want_marker);
        }
        other => panic!("{entry:?} {bits:#010X}: expected ReservedFloatValue, got {other:?}"),
    }
}

// r[verify record_encoder.reserved_rejection]
/// A rejected call writes nothing. A caller that handles the error and
/// finishes the record without that field gets exactly the bytes it would have
/// got had it never made the call — in every format and at every entry point.
/// Re-encoding the field would prove nothing: a second call replaces the
/// first, bytes and all. For QUAL the writer must survive instead, and the
/// record is begun again with a valid QUAL.
#[hegel::test]
fn a_rejected_value_leaves_no_trace_in_the_output(tc: TestCase) {
    let entry = tc.draw(gs::sampled_from([INT_ENTRIES, FLOAT_ENTRIES].concat()).print_as_debug());
    let format = tc.draw(
        gs::sampled_from(&[OutputFormat::Vcf, OutputFormat::VcfGz, OutputFormat::Bcf])
            .print_as_debug(),
    );
    let int = tc.draw(gs::integers::<i32>().min_value(i32::MIN + INT_RESERVED_LEN));
    let float = tc.draw(gs::floats::<f32>().min_value(-1e6).max_value(1e6));
    let int_offset = tc.draw(gs::integers::<i32>().min_value(0).max_value(INT_RESERVED_LEN - 1));
    let float_offset = tc.draw(gs::integers::<u32>().max_value(FLOAT_RESERVED_LEN - 1));
    let reserved = (i32::MIN + int_offset, f32::from_bits(FLOAT_RESERVED_FIRST + float_offset));

    let s = setup();
    let without_the_call = write_one(&s, format, entry, int, float, Call::Omitted).unwrap();
    let after_rejection =
        write_one(&s, format, entry, int, float, Call::Rejected(reserved.0, reserved.1)).unwrap();
    assert!(
        without_the_call == after_rejection,
        "{entry:?} in {format:?}: the rejected call left bytes behind"
    );
}

// r[verify bcf_encoder.reserved_uniform]
/// Plain VCF text could render `i32::MIN` losslessly, but accepting it there
/// and refusing it in BCF would mean the same encoding calls produce different
/// records per `OutputFormat`. All three formats must agree, value for value.
#[hegel::test]
fn the_three_output_formats_accept_exactly_the_same_integers(tc: TestCase) {
    let value = tc.draw(hegel::one_of!(
        gs::integers::<i32>(),
        // Straddle the band: drawing over all of i32 reaches it once in 2^29.
        gs::integers::<i32>().min_value(i32::MIN).max_value(i32::MIN + 16),
    ));
    let entry = tc.draw(gs::sampled_from(&INT_ENTRIES).print_as_debug());

    let s = setup();
    let accepted: Vec<bool> = [OutputFormat::Vcf, OutputFormat::VcfGz, OutputFormat::Bcf]
        .into_iter()
        .map(|f| write_one(&s, f, entry, value, 30.0, Call::Value).is_ok())
        .collect();
    assert!(
        accepted.iter().all(|a| *a == accepted[0]),
        "{entry:?} value {value}: formats disagree — {accepted:?}"
    );
    assert_eq!(accepted[0], value > i32::MIN + INT_RESERVED_LEN - 1);
}

// ── Why any of this matters ────────────────────────────────────────────

// r[verify bcf_encoder.reserved_rejected]
/// The loss the rejection prevents, shown on a file bcftools reads.
///
/// The encoder will no longer emit `0x80000001` as data, so the file is built
/// with an ordinary value and that value's four bytes are patched afterwards —
/// which is exactly the file any writer without this check would produce. The
/// row goes in with three elements and comes out with one; nothing in the file
/// or in bcftools' exit status says anything was dropped.
#[test]
fn reserved_end_of_vector_written_as_data_truncates_the_row() {
    const VICTIM: i32 = 0x1234_5678;
    const END_OF_VECTOR: u32 = 0x8000_0001;

    let s = setup();
    let intact = write_one_with_ad(&s, &[1_000_000, VICTIM, 999]);
    assert_eq!(
        bcftools_query(&intact, "%INFO/AD\n"),
        vec![format!("1000000,{VICTIM},999")],
        "the row reads back whole before patching"
    );

    // Swap the middle element for END_OF_VECTOR inside the BGZF payload.
    let mut plain = bgzf_decompress(&intact);
    let needle = VICTIM.to_le_bytes();
    let hits: Vec<usize> =
        plain.windows(4).enumerate().filter(|(_, w)| *w == needle).map(|(i, _)| i).collect();
    assert_eq!(hits.len(), 1, "the sentinel value must occur exactly once");
    plain[hits[0]..hits[0] + 4].copy_from_slice(&END_OF_VECTOR.to_le_bytes());
    let patched = bgzf_compress(&plain);

    let got = bcftools_query(&patched, "%INFO/AD\n");
    assert_eq!(
        got,
        vec!["1000000".to_string()],
        "END_OF_VECTOR as data truncates the row at that element — silently"
    );
}

/// One record whose INFO/AD is `values`, as BCF bytes.
fn write_one_with_ad(s: &Setup, values: &[i32]) -> Vec<u8> {
    let mut out = Vec::new();
    {
        let writer = Writer::new(&mut out, OutputFormat::Bcf);
        let mut writer = writer.write_header(&s.header).unwrap();
        let alleles = Alleles::snv(Base::A, Base::T).unwrap();
        let mut enc = writer
            .begin_record(&s.contig, Pos1::new(100).unwrap(), &alleles, Some(30.0))
            .unwrap()
            .filter_pass();
        s.info_ints.encode(&mut enc, values).unwrap();
        let mut enc = enc.begin_samples();
        s.gt.encode(&mut enc, &[Genotype::unphased(0, 1)]).unwrap();
        enc.emit().unwrap();
        writer.finish().unwrap();
    }
    out
}

fn bgzf_decompress(compressed: &[u8]) -> Vec<u8> {
    use std::io::Read;
    let mut reader = bgzf::Reader::new(compressed);
    let mut out = Vec::new();
    reader.read_to_end(&mut out).unwrap();
    out
}

fn bgzf_compress(plain: &[u8]) -> Vec<u8> {
    use std::io::Write;
    let mut out = Vec::new();
    let mut writer = bgzf::Writer::new(&mut out, bgzf::CompressionLevel::new(6).unwrap());
    writer.write_all(plain).unwrap();
    writer.finish().unwrap();
    out
}
