//! `format_sam_record` against htslib's own `sam_format1`, byte for byte.
//!
//! Records are generated in wire form (`support/raw_record_gen.rs`), written
//! to a BAM file with `BamWriter::write_raw`, and read back by htslib, which
//! formats each one through `sam_format1` over FFI; seqair streams the same
//! file with `BamStreamReader` and formats each record. The two lines must be
//! identical. The float printers are pinned to their C originals the same way:
//! `format_g` to libc's `snprintf("%g")`, `format_kputd` to htslib's `kputd`.
//! `samtools view` and the htslib SAM test corpus are the end-to-end checks.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    unsafe_code,
    reason = "test code; FFI to htslib and libc"
)]

#[path = "../src/pinned.rs"]
mod pinned;
#[path = "support/raw_record_gen.rs"]
mod raw_record_gen;

use hegel::TestCase;
use hegel::generators as gs;
use raw_record_gen::{GenRecord, HTSLIB, arb_f32, arb_f64, arb_record};
use rust_htslib::bam::{self as htsbam, Read as _};
use rust_htslib::htslib;
use seqair::bam::aux::AuxLayoutError;
use seqair::bam::stream::BamStreamReader;
use seqair::bam::writer::BamWriterBuilder;
use seqair::bam::{BamHeader, RawBamRecord};
use seqair::sam::format::{
    ReferenceField, SamFormatError, format_g, format_kputd, format_sam_header, format_sam_record,
};
use std::path::{Path, PathBuf};
use std::process::Command;

// ---- float printers ------------------------------------------------------

fn libc_g(value: f64) -> Vec<u8> {
    let mut buf = [0u8; 64];
    // SAFETY: `buf` is writable for its length and the format takes one double.
    let n = unsafe { libc::snprintf(buf.as_mut_ptr().cast(), buf.len(), c"%g".as_ptr(), value) };
    buf[..usize::try_from(n).unwrap()].to_vec()
}

/// `%g` output with the trailing zeros of an exponent form's mantissa
/// removed, as the C standard requires. macOS's libc keeps them when the value
/// is an exact tie rounded down (`1000005` prints `1.00000e+06`, glibc `1e+06`).
fn strip_mantissa_zeros(g: &[u8]) -> Vec<u8> {
    let Some(e) = g.iter().position(|&b| b == b'e') else { return g.to_vec() };
    let (mantissa, exponent) = g.split_at(e);
    if !mantissa.contains(&b'.') {
        return g.to_vec();
    }
    let mut out = mantissa.to_vec();
    while out.last() == Some(&b'0') {
        out.pop();
    }
    if out.last() == Some(&b'.') {
        out.pop();
    }
    out.extend_from_slice(exponent);
    out
}

/// What `format_g` must print, given libc's output for `value`: the same,
/// except where macOS's libc departs from glibc and musl, which seqair
/// follows: a negative NaN prints `nan` there (`-nan` here), and the exact-tie
/// trailing zeros of [`strip_mantissa_zeros`].
fn expect_from_libc(value: f64, libc: Vec<u8>) -> Vec<u8> {
    if !cfg!(target_os = "macos") {
        return libc;
    }
    if value.is_nan() && value.is_sign_negative() {
        assert_eq!(libc, b"nan");
        return b"-nan".to_vec();
    }
    strip_mantissa_zeros(&libc)
}

/// A SAM line with macOS's `%g` trailing zeros removed from its `f`, `d` and
/// `B:f` values (see [`strip_mantissa_zeros`]); the line itself elsewhere.
/// Applied to both sides of a comparison, so a `Z` value that happens to look
/// like a float tag changes alike on both.
fn normalise_floats(line: &[u8]) -> Vec<u8> {
    if !cfg!(target_os = "macos") {
        return line.to_vec();
    }
    let fields: Vec<Vec<u8>> = line
        .split(|&b| b == b'\t')
        .map(|field| match field.get(2..) {
            Some([b':', b'f' | b'd', b':', ..]) => {
                [&field[..5], &strip_mantissa_zeros(&field[5..])[..]].concat()
            }
            Some([b':', b'B', b':', b'f', ..]) => {
                let mut out = field[..6].to_vec();
                for element in field[6..].split(|&b| b == b',').skip(1) {
                    out.push(b',');
                    out.extend(strip_mantissa_zeros(element));
                }
                out
            }
            _ => field.to_vec(),
        })
        .collect();
    fields.join(&b'\t')
}

fn kstring_bytes(ks: &mut htslib::kstring_t) -> Vec<u8> {
    let out = if ks.s.is_null() {
        Vec::new()
    } else {
        // SAFETY: htslib wrote `l` bytes at `s`.
        unsafe { std::slice::from_raw_parts(ks.s.cast::<u8>(), ks.l) }.to_vec()
    };
    // SAFETY: htslib allocates kstrings with malloc.
    unsafe { libc::free(ks.s.cast()) };
    ks.s = std::ptr::null_mut();
    out
}

fn htslib_kputd(value: f64) -> Vec<u8> {
    let mut ks = htslib::kstring_t { l: 0, m: 0, s: std::ptr::null_mut() };
    // SAFETY: an empty kstring is valid input; kputd grows it.
    let r = unsafe { htslib::kputd(value, &mut ks) };
    assert!(r >= 0, "kputd failed");
    kstring_bytes(&mut ks)
}

fn ours(f: fn(f64, &mut Vec<u8>), value: f64) -> Vec<u8> {
    // Something already in `out` must stay there.
    let mut out = b"x".to_vec();
    f(value, &mut out);
    assert_eq!(out[0], b'x');
    out.split_off(1)
}

// r[verify sam.format.float_g]
#[hegel::test]
fn format_g_matches_libc_for_doubles(tc: TestCase) {
    let value = arb_f64(&tc, true);
    let libc = expect_from_libc(value, libc_g(value));
    assert_eq!(
        ours(format_g, value).escape_ascii().to_string(),
        libc.escape_ascii().to_string(),
        "{value:e}"
    );
}

// r[verify sam.format.float_g]
#[hegel::test]
fn format_g_matches_libc_for_floats(tc: TestCase) {
    let value = f64::from(arb_f32(&tc, true));
    let libc = expect_from_libc(value, libc_g(value));
    assert_eq!(
        ours(format_g, value).escape_ascii().to_string(),
        libc.escape_ascii().to_string(),
        "{value:e}"
    );
}

// r[verify sam.format.float_kputd]
#[hegel::test]
fn format_kputd_matches_htslib(tc: TestCase) {
    // `B:f` elements are f32, but kputd takes any double.
    let value =
        if tc.draw(gs::booleans()) { f64::from(arb_f32(&tc, true)) } else { arb_f64(&tc, true) };
    let hts = expect_from_libc(value, htslib_kputd(value));
    assert_eq!(
        ours(format_kputd, value).escape_ascii().to_string(),
        hts.escape_ascii().to_string(),
        "{value:e}"
    );
}

// r[verify sam.format.float_g]
// r[verify sam.format.float_kputd]
#[test]
fn float_printers_on_known_values() {
    let g = |v| String::from_utf8(ours(format_g, v)).unwrap();
    let k = |v| String::from_utf8(ours(format_kputd, v)).unwrap();
    for (value, expected) in [
        (0.0, "0"),
        (-0.0, "-0"),
        (1.0, "1"),
        (0.5, "0.5"),
        (1e-5, "1e-05"),
        (0.0001, "0.0001"),
        (123_456.0, "123456"),
        (1_234_567.0, "1.23457e+06"),
        (1_234_565.0, "1.23456e+06"),
        (9.999_995, "10"),
        (1e100, "1e+100"),
        (f64::INFINITY, "inf"),
        (f64::NEG_INFINITY, "-inf"),
        (f64::NAN, "nan"),
        (-f64::NAN, "-nan"),
        (f64::from_bits(1), "4.94066e-324"),
    ] {
        assert_eq!(g(value), expected, "%g of {value:e}");
    }
    // kputd agrees with %g except where its scaled rounding differs.
    // On doubles near a tie the scaled product rounds the other way (never on
    // an f32, which is all `B:f` holds: every f32 in kputd's range agrees).
    assert_eq!((g(1.000_005), k(1.000_005)), ("1.00001".into(), "1".into()));
    assert_eq!((g(1000.015), k(1000.015)), ("1000.01".into(), "1000.02".into()));
    assert_eq!((g(1.000_045), k(1.000_045)), ("1.00005".into(), "1.00004".into()));
    assert_eq!(k(999_999.0), "999999");
    assert_eq!(k(999_999.5), "1e+06");
    assert_eq!(k(-0.0), "-0");
    assert_eq!(k(-f64::NAN), "-nan");
    assert_eq!(k(f64::NEG_INFINITY), "-inf");
}

// ---- records -------------------------------------------------------------

fn write_bam(path: &Path, header: &BamHeader, records: &[RawBamRecord]) {
    let mut writer = BamWriterBuilder::to_path(path, header).build().unwrap();
    for record in records {
        writer.write_raw(record).unwrap();
    }
    writer.finish().unwrap();
}

/// What htslib does with each record of `path`: its `sam_format1` line, or
/// `None` where reading or formatting failed.
fn htslib_lines(path: &Path) -> Vec<Option<Vec<u8>>> {
    let mut reader = htsbam::Reader::from_path(path).unwrap();
    let header = reader.header().inner_ptr();
    let mut record = htsbam::Record::new();
    let mut lines = Vec::new();
    loop {
        match reader.read(&mut record) {
            None => break,
            Some(Err(_)) => {
                lines.push(None);
                // A failed record read leaves htslib's stream misaligned.
                break;
            }
            Some(Ok(())) => {}
        }
        let mut ks = htslib::kstring_t { l: 0, m: 0, s: std::ptr::null_mut() };
        // SAFETY: header and record are live; sam_format1 fills the kstring.
        let r = unsafe { htslib::sam_format1(header, record.inner(), &mut ks) };
        let line = kstring_bytes(&mut ks);
        lines.push((r >= 0).then_some(line));
    }
    lines
}

/// Stream `path` with seqair and format every record.
fn seqair_lines(path: &Path) -> (BamHeader, Vec<Result<Vec<u8>, SamFormatError>>) {
    let mut reader = BamStreamReader::new(std::fs::File::open(path).unwrap()).unwrap();
    let header = BamHeader::from_template(reader.header());
    let mut record = RawBamRecord::default();
    let mut lines = Vec::new();
    while reader.read_record(&mut record).unwrap() {
        let mut line = Vec::new();
        lines.push(format_sam_record(&record, &header, &mut line).map(|()| line));
    }
    (header, lines)
}

fn show(line: &[u8]) -> String {
    line.escape_ascii().to_string()
}

fn temp_bam(dir: &tempfile::TempDir, name: &str) -> PathBuf {
    dir.path().join(name)
}

fn validated(header: &BamHeader, records: &[GenRecord]) -> Vec<RawBamRecord> {
    records
        .iter()
        .map(|g| {
            RawBamRecord::from_payload(g.payload(), header).expect("generated records are valid")
        })
        .collect()
}

// r[verify sam.format]
// r[verify sam.format.qname]
// r[verify sam.format.fixed_fields]
// r[verify sam.format.cigar]
// r[verify sam.format.long_cigar]
// r[verify sam.format.seq_qual]
// r[verify sam.format.aux]
#[hegel::test(test_cases = pinned::cases(300))]
fn format_matches_htslib_sam_format1(tc: TestCase) {
    let generated: Vec<GenRecord> = tc.draw(gs::vecs(arb_record(HTSLIB)).min_size(1).max_size(6));
    let header = raw_record_gen::header();
    let records = validated(&header, &generated);
    let dir = tempfile::tempdir().unwrap();
    let path = temp_bam(&dir, "records.bam");
    write_bam(&path, &header, &records);

    let hts = htslib_lines(&path);
    let (read_header, ours) = seqair_lines(&path);
    assert_eq!(hts.len(), records.len(), "htslib read every record");
    assert_eq!(ours.len(), records.len());
    for ((record, hts), ours) in records.iter().zip(&hts).zip(&ours) {
        let hts = hts.as_ref().expect("htslib formats every generated record");
        let ours = ours.as_ref().unwrap();
        assert_eq!(show(&normalise_floats(ours)), show(&normalise_floats(hts)));
        // The record seqair streamed back is the one it wrote.
        let mut direct = Vec::new();
        format_sam_record(record, &read_header, &mut direct).unwrap();
        assert_eq!(direct, *ours);
    }
}

fn samtools(args: &[&str]) -> Vec<u8> {
    let out = Command::new("samtools").args(args).output().expect("samtools is installed");
    assert!(out.status.success(), "samtools {args:?}: {}", out.stderr.escape_ascii());
    out.stdout
}

/// The header and records of `path` as seqair prints them, a newline after
/// each record, as `samtools view -h --no-PG` prints the file.
fn seqair_view(path: &Path) -> Vec<u8> {
    let (header, lines) = seqair_lines(path);
    let mut out = Vec::new();
    format_sam_header(&header, &mut out).unwrap();
    for line in lines {
        out.extend(line.unwrap());
        out.push(b'\n');
    }
    out
}

fn assert_same_text(ours: &[u8], theirs: &[u8]) {
    let lines = |text: &[u8]| -> Vec<String> {
        text.split(|&b| b == b'\n').map(|l| show(&normalise_floats(l))).collect()
    };
    let (ours, theirs) = (lines(ours), lines(theirs));
    assert_eq!(ours, theirs);
}

// r[verify sam.format]
// r[verify sam.format.header]
#[hegel::test(test_cases = pinned::cases(15))]
fn format_matches_samtools_view(tc: TestCase) {
    let generated: Vec<GenRecord> = tc.draw(gs::vecs(arb_record(HTSLIB)).min_size(1).max_size(10));
    let header = raw_record_gen::header();
    let records = validated(&header, &generated);
    let dir = tempfile::tempdir().unwrap();
    let path = temp_bam(&dir, "records.bam");
    write_bam(&path, &header, &records);
    let path_str = path.to_str().unwrap();
    assert_same_text(&seqair_view(&path), &samtools(&["view", "-h", "--no-PG", path_str]));
}

// r[verify sam.format]
// r[verify sam.format.header]
#[test]
fn htslib_sam_corpus_formats_like_samtools() {
    let corpus = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/htslib/sam");
    let dir = tempfile::tempdir().unwrap();
    let mut checked = 0;
    for entry in std::fs::read_dir(&corpus).unwrap() {
        let sam = entry.unwrap().path();
        if sam.extension().is_none_or(|e| e != "sam") {
            continue;
        }
        let bam = dir.path().join("corpus.bam");
        let converted = Command::new("samtools")
            .args(["view", "--no-PG", "-b", "-o"])
            .arg(&bam)
            .arg(&sam)
            .output()
            .unwrap();
        assert!(converted.status.success(), "samtools could not convert {}", sam.display());
        let expected = samtools(&["view", "-h", "--no-PG", bam.to_str().unwrap()]);
        let ours = seqair_view(&bam);
        assert_eq!(
            ours.split(|&b| b == b'\n').map(show).collect::<Vec<_>>(),
            expected.split(|&b| b == b'\n').map(show).collect::<Vec<_>>(),
            "{}",
            sam.display()
        );
        checked += 1;
    }
    assert!(checked >= 40, "only {checked} corpus files found");
}

// ---- header ----------------------------------------------------------

/// A BAM file with exactly this header text and these targets, and no records.
fn bam_with_raw_header(path: &Path, text: &[u8], targets: &[(&str, i32)]) {
    let mut raw = b"BAM\x01".to_vec();
    raw.extend_from_slice(&i32::try_from(text.len()).unwrap().to_le_bytes());
    raw.extend_from_slice(text);
    raw.extend_from_slice(&i32::try_from(targets.len()).unwrap().to_le_bytes());
    for (name, len) in targets {
        raw.extend_from_slice(&i32::try_from(name.len() + 1).unwrap().to_le_bytes());
        raw.extend_from_slice(name.as_bytes());
        raw.push(0);
        raw.extend_from_slice(&len.to_le_bytes());
    }
    let mut out = Vec::new();
    let mut writer = seqair::io::BgzfWriter::new(&mut out);
    writer.write_all(&raw).unwrap();
    writer.finish().unwrap();
    std::fs::write(path, out).unwrap();
}

// r[verify sam.format.header]
#[hegel::test(test_cases = pinned::cases(40))]
fn header_prints_like_samtools_view(tc: TestCase) {
    let targets: Vec<(&str, i32)> = [("chr1", 248_956_422), ("chrM", 16569), ("x", 1)]
        [..tc.draw(gs::integers::<usize>().max_value(3))]
        .to_vec();
    let mut lines: Vec<Vec<u8>> = Vec::new();
    if tc.draw(gs::booleans()) {
        lines.push(b"@HD\tVN:1.6\tSO:unsorted".to_vec());
    }
    let with_sq = tc.draw(gs::booleans());
    for (name, len) in &targets {
        if with_sq {
            lines.push(format!("@SQ\tSN:{name}\tLN:{len}").into_bytes());
        }
    }
    for _ in 0..tc.draw(gs::integers::<usize>().max_value(2)) {
        let text =
            tc.draw(gs::vecs(gs::integers::<u8>().min_value(b' ').max_value(b'~')).max_size(12));
        lines.push([&b"@CO\t"[..], &text].concat());
    }
    // A line htslib rejects: one not starting with `@`, or an empty one.
    let malformed = !lines.is_empty() && tc.draw(gs::weighted_booleans(0.2));
    if malformed {
        // Before another line: an empty last line would be the final newline.
        let at = tc.draw(gs::integers::<usize>().max_value(lines.len() - 1));
        let bad = tc.draw(gs::sampled_from(vec![b"".to_vec(), b"CO\tno at".to_vec()]));
        lines.insert(at, bad);
    }
    let mut text = lines.join(&b'\n');
    if !text.is_empty() && tc.draw(gs::booleans()) {
        text.push(b'\n');
    }
    // BAM writers pad the text with NULs.
    let padding = tc.draw(gs::integers::<usize>().max_value(3));
    text.extend(std::iter::repeat_n(0, padding));

    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("header.bam");
    bam_with_raw_header(&path, &text, &targets);
    let reader = BamStreamReader::new(std::fs::File::open(&path).unwrap()).unwrap();
    let mut ours = b"kept".to_vec();
    let result = format_sam_header(reader.header(), &mut ours);
    let samtools =
        Command::new("samtools").args(["view", "-h", "--no-PG"]).arg(&path).output().unwrap();
    if malformed {
        assert!(
            matches!(result, Err(SamFormatError::MalformedHeaderLine { .. })),
            "{result:?} for {}",
            text.escape_ascii()
        );
        assert_eq!(ours, b"kept");
        assert!(!samtools.status.success(), "samtools accepted {}", text.escape_ascii());
        return;
    }
    result.unwrap();
    assert!(samtools.status.success(), "{}", samtools.stderr.escape_ascii());
    // htslib copies the NUL padding into its output; seqair stops at the first NUL.
    let expected: Vec<u8> = samtools.stdout.iter().copied().filter(|&b| b != 0).collect();
    assert_eq!(show(&ours[4..]), show(&expected));
}

// ---- errors --------------------------------------------------------------

/// How to break a record's aux block.
#[derive(Debug, Clone, Copy)]
enum Damage {
    /// A tag with a type code that is not one.
    UnknownType,
    /// A `B` array with an element type that is not one (`A` included).
    UnknownSubtype,
    /// A `Z` or `H` string the block ends inside.
    Unterminated,
    /// A whole tag cut short.
    Truncated,
    /// 1–3 stray bytes after the last tag.
    Trailing,
}

/// Every byte but those in `codes`.
fn not_in(codes: &[u8]) -> Vec<u8> {
    (0..=255).filter(|c| !codes.contains(c)).collect()
}

/// Append a malformed tag to `aux`; returns how many bytes the malformed
/// part has (htslib stops quietly before fewer than 4).
fn damage(tc: &TestCase, damage: Damage, aux: &mut Vec<u8>) -> usize {
    let start = aux.len();
    match damage {
        Damage::UnknownType => {
            let code = tc.draw(gs::sampled_from(not_in(b"AcCsSiIfdZHB")));
            aux.extend_from_slice(&[b'X', b'X', code]);
            aux.extend(tc.draw(gs::vecs(gs::integers::<u8>()).min_size(1).max_size(8)));
        }
        Damage::UnknownSubtype => {
            let code = tc.draw(gs::sampled_from(not_in(b"cCsSiIf")));
            aux.extend_from_slice(&[b'X', b'X', b'B', code, 1, 0, 0, 0, 7, 7, 7, 7]);
        }
        Damage::Unterminated => {
            let kind = tc.draw(gs::sampled_from(b"ZH".to_vec()));
            aux.extend_from_slice(&[b'X', b'X', kind]);
            aux.extend(tc.draw(gs::vecs(gs::integers::<u8>().min_value(1)).max_size(8)));
        }
        Damage::Truncated => {
            let mut tag = Vec::new();
            raw_record_gen::arb_tag(tc, HTSLIB).encode(&mut tag);
            let keep = tc.draw(gs::integers::<usize>().min_value(1).max_value(tag.len() - 1));
            aux.extend_from_slice(&tag[..keep]);
        }
        Damage::Trailing => {
            aux.extend(tc.draw(gs::vecs(gs::integers::<u8>()).min_size(1).max_size(3)));
        }
    }
    aux.len() - start
}

// r[verify sam.format.aux_errors]
#[hegel::test(test_cases = pinned::cases(200))]
fn malformed_aux_is_a_typed_error(tc: TestCase) {
    let generated: GenRecord = tc.draw(arb_record(HTSLIB));
    let header = raw_record_gen::header();
    let intact = RawBamRecord::from_payload(generated.payload(), &header).unwrap();
    let kind = tc.draw(gs::sampled_from(vec![
        Damage::UnknownType,
        Damage::UnknownSubtype,
        Damage::Unterminated,
        Damage::Truncated,
        Damage::Trailing,
    ]));
    // The aux block ends the payload, so the damage goes after the last tag.
    let mut payload = generated.payload();
    let damaged_len = damage(&tc, kind, &mut payload);
    // Aux is not validated on read (`r[bam.raw_record.aux_lazy]`).
    let record = RawBamRecord::from_payload(payload, &header).unwrap();

    let mut out = b"kept".to_vec();
    let err = format_sam_record(&record, &header, &mut out).unwrap_err();
    assert!(matches!(err, SamFormatError::Aux(_)), "{err:?}");
    assert_eq!(out, b"kept", "a failed format leaves `out` as it was");

    let dir = tempfile::tempdir().unwrap();
    let path = temp_bam(&dir, "damaged.bam");
    write_bam(&path, &header, std::slice::from_ref(&record));
    let hts = htslib_lines(&path);
    match hts.as_slice() {
        [None] => {}
        [Some(line)] => {
            // htslib only gets through when the damage is under four bytes,
            // which it drops without a word: the line of the intact record.
            assert!(damaged_len < 4, "htslib formatted a {damaged_len}-byte malformed tag");
            let mut expected = Vec::new();
            format_sam_record(&intact, &header, &mut expected).unwrap();
            assert_eq!(show(&normalise_floats(line)), show(&normalise_floats(&expected)));
        }
        other => panic!("one record expected, htslib gave {other:?}"),
    }
}

#[test]
fn every_kind_of_malformed_aux_names_its_cause() {
    let header = raw_record_gen::header();
    let base = |aux: &[u8]| {
        let mut payload = GenRecord {
            refid: -1,
            pos: -1,
            mapq: 255,
            bin: 4680,
            flags: 4,
            next_refid: -1,
            next_pos: -1,
            tlen: 0,
            qname: b"r".to_vec(),
            cigar: Vec::new(),
            seq: Vec::new(),
            seq_pad: 0,
            qual: Vec::new(),
            tags: Vec::new(),
        }
        .payload();
        payload.extend_from_slice(aux);
        let record = RawBamRecord::from_payload(payload, &header).unwrap();
        format_sam_record(&record, &header, &mut Vec::new())
    };
    assert_eq!(base(b"NMi\x01\x00\x00\x00"), Ok(()));
    let err = |e| Err(SamFormatError::Aux(e));
    assert_eq!(base(b"XXq\x01\x00"), err(AuxLayoutError::UnknownType { offset: 0, code: b'q' }));
    assert_eq!(
        base(b"NMC\x01XXBA\x01\x00\x00\x00a"),
        err(AuxLayoutError::UnknownArraySubtype { offset: 4, code: b'A' })
    );
    assert_eq!(base(b"XXZabc"), err(AuxLayoutError::UnterminatedString { offset: 0 }));
    assert_eq!(base(b"XXi\x01\x00"), err(AuxLayoutError::Truncated { offset: 0 }));
    assert_eq!(
        base(b"XXB\x66\x02\x00\x00\x00\x00\x00\x00"),
        err(AuxLayoutError::Truncated { offset: 0 })
    );
    assert_eq!(base(b"NMC\x01XX"), err(AuxLayoutError::Truncated { offset: 4 }));
}

// r[verify sam.format.unknown_reference]
#[test]
fn a_reference_the_header_lacks_is_a_typed_error() {
    let header = raw_record_gen::header();
    let small = BamHeader::from_sam_text("@SQ\tSN:chr1\tLN:100\n").unwrap();
    let record = |refid: i32, next_refid: i32| {
        let g = GenRecord {
            refid,
            pos: if refid < 0 { -1 } else { 5 },
            mapq: 0,
            bin: 0,
            flags: 4,
            next_refid,
            next_pos: -1,
            tlen: 0,
            qname: b"r".to_vec(),
            cigar: Vec::new(),
            seq: Vec::new(),
            seq_pad: 0,
            qual: Vec::new(),
            tags: Vec::new(),
        };
        RawBamRecord::from_payload(g.payload(), &header).unwrap()
    };
    let fails = |r: &RawBamRecord| {
        let mut out = b"kept".to_vec();
        let e = format_sam_record(r, &small, &mut out).unwrap_err();
        assert_eq!(out, b"kept");
        e
    };
    let mut out = Vec::new();
    format_sam_record(&record(0, 0), &small, &mut out).unwrap();
    assert!(out.starts_with(b"r\t4\tchr1\t6\t0\t*\t=\t0\t0\t*\t*"), "{}", show(&out));
    assert!(matches!(
        fails(&record(2, -1)),
        SamFormatError::UnknownReference { field: ReferenceField::Rname, target_count: 1, .. }
    ));
    assert!(matches!(
        fails(&record(0, 1)),
        SamFormatError::UnknownReference { field: ReferenceField::Rnext, target_count: 1, .. }
    ));
    // RNEXT `=` needs no lookup, but RNAME does.
    assert!(matches!(
        fails(&record(1, 1)),
        SamFormatError::UnknownReference { field: ReferenceField::Rname, .. }
    ));
}
