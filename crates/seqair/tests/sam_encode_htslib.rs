//! `SamRecordEncoder` against htslib (samtools 1.24's `sam_parse1` +
//! `bam_write1`, and rust-htslib's bundled htslib) and noodles.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    reason = "test code"
)]

#[path = "../src/pinned.rs"]
mod pinned;
#[path = "support/sam_lines.rs"]
mod sam_lines;

use hegel::TestCase;
use hegel::generators as gs;
use sam_lines::{
    Dialect, Line, canonical_line, canonical_lines, header_text, rust_htslib_payloads, sam_text,
    samtools_payloads, seqair_payloads,
};
use seqair::bam::{BamHeader, CigarOp, CigarStr, RawBamRecord};
use seqair::sam::{SamRecordEncoder, SamRecordError, parse_cigar};

const HTSLIB: Dialect = Dialect { noodles: false };
const NOODLES: Dialect = Dialect { noodles: true };

fn header() -> BamHeader {
    BamHeader::from_sam_text(&header_text()).unwrap()
}

fn encode(line: &str) -> Result<Vec<u8>, SamRecordError> {
    let mut record = RawBamRecord::default();
    SamRecordEncoder::new().encode(line.as_bytes(), &header(), &mut record)?;
    Ok(record.into_payload())
}

/// Compare seqair's payloads with an oracle's, naming the first difference.
fn assert_same(ours: &[Vec<u8>], theirs: &[Vec<u8>], lines: &[Line], oracle: &str) {
    assert_eq!(ours.len(), theirs.len(), "{oracle}: record count");
    for (i, (a, b)) in ours.iter().zip(theirs).enumerate() {
        if a != b {
            let at = a.iter().zip(b).position(|(x, y)| x != y).unwrap_or(a.len().min(b.len()));
            panic!(
                "{oracle}: record {i} differs at byte {at} (ours {} bytes, theirs {})\nline: {}\nours:   {:02x?}\ntheirs: {:02x?}",
                a.len(),
                b.len(),
                lines[i].text().escape_debug(),
                &a[..a.len().min(at + 16)],
                &b[..b.len().min(at + 16)],
            );
        }
    }
}

fn has_uracil(line: &Line) -> bool {
    line.fields[9].contains(['U', 'u'])
}

// r[verify sam.encode]
// r[verify sam.encode.fixups]
// r[verify sam.encode.bin]
// r[verify sam.encode.seq]
// r[verify sam.encode.qual]
// r[verify sam.encode.references+2]
// r[verify sam.aux.types]
/// Canonical lines — every placement htslib fixes up, IUPAC in both cases,
/// `=`, `U`, missing SEQ/QUAL, every aux type and B subtype incl. empty
/// arrays, `=` RNEXT, MAPQ 255, POS up to 2^31 — encode to exactly the bytes
/// `samtools view -b` writes, and to rust-htslib's (bar `U`, which its
/// htslib 1.19 packs as `N`).
#[hegel::test(test_cases = crate::pinned::cases(150))]
fn canonical_lines_encode_to_htslibs_bytes(tc: TestCase) {
    let lines = tc.draw(canonical_lines(HTSLIB, 8));
    let sam = sam_text(&lines);
    let ours = seqair_payloads(sam.as_bytes()).unwrap_or_else(|e| panic!("seqair: {e:?}\n{sam}"));
    let samtools =
        samtools_payloads(sam.as_bytes()).unwrap_or_else(|e| panic!("samtools: {e}\n{sam}"));
    assert_same(&ours, &samtools, &lines, "samtools");

    let hts =
        rust_htslib_payloads(sam.as_bytes()).unwrap_or_else(|e| panic!("rust-htslib: {e}\n{sam}"));
    assert_eq!(hts.len(), ours.len());
    let (mut ours_cmp, mut hts_cmp, mut lines_cmp) = (vec![], vec![], vec![]);
    for ((o, h), l) in ours.iter().zip(&hts).zip(&lines) {
        if !has_uracil(l) {
            ours_cmp.push(o.clone());
            hts_cmp.push(h.clone());
            lines_cmp.push(l.clone());
        }
    }
    assert_same(&ours_cmp, &hts_cmp, &lines_cmp, "rust-htslib");
    tc.event_value("records", lines.len() as f64);
}

// r[verify sam.encode.long_cigar]
// r[verify sam.encode.cigar_seq]
/// More than 65535 CIGAR operations: the placeholder `<l_seq>S<rlen>N` and
/// `CG:B:I` after every other tag, exactly as `bam_write1` writes them.
#[hegel::test(test_cases = crate::pinned::cases(6))]
fn long_cigars_encode_to_htslibs_bytes(tc: TestCase) {
    let n_pairs = tc.draw(gs::integers::<usize>().min_value(32_768).max_value(33_500));
    let second = tc.draw(gs::sampled_from(&["I", "D", "N", "X"]));
    let cigar = format!("1M1{second}").repeat(n_pairs);
    let qlen = if matches!(second, "I" | "X") { 2 * n_pairs } else { n_pairs };
    let with_seq = tc.draw(gs::booleans());
    let (seq, qual) = if with_seq {
        ("A".repeat(qlen), "I".repeat(qlen))
    } else {
        ("*".to_string(), "*".to_string())
    };
    let flag = tc.draw(gs::sampled_from(&[0u16, 4, 16]));
    let line = Line {
        fields: [
            "long",
            &flag.to_string(),
            "chr2",
            "100",
            "60",
            &cigar,
            "*",
            "0",
            "0",
            &seq,
            &qual,
            "NM:i:3",
            "XZ:Z:after",
        ]
        .map(str::to_string)
        .to_vec(),
    };
    let lines = [line];
    let sam = sam_text(&lines);
    let ours = seqair_payloads(sam.as_bytes()).unwrap();
    assert_same(&ours, &samtools_payloads(sam.as_bytes()).unwrap(), &lines, "samtools");
    assert_same(&ours, &rust_htslib_payloads(sam.as_bytes()).unwrap(), &lines, "rust-htslib");
    let record = RawBamRecord::from_payload(ours[0].clone(), &header()).unwrap();
    assert_eq!(record.stored_cigar_bytes().len(), 2 * 4);
    // r[verify bam.raw_record.cigar_access]
    // The real CIGAR comes back from the `CG` tag wherever htslib's reader
    // would move it there.
    if record.has_cg_placeholder() {
        assert!(record.cigar_len() > 65535);
        assert_eq!(record.cigar().len(), record.cigar_len());
    } else {
        assert_eq!(record.cigar_len(), 2);
    }
    assert!(record.has_cg_placeholder() || flag == 4 || !with_seq);
}

// r[verify sam.encode.cg_tag]
/// A SAM that carries its CIGAR in a `CG:B:I` tag (placeholder CIGAR) gets
/// it moved back, as `bam_tag2cigar` does — and, when it is not longer than
/// 65535 operations, written inline without the tag.
#[hegel::test(test_cases = crate::pinned::cases(40))]
fn cg_tags_move_into_the_cigar_like_htslib(tc: TestCase) {
    let (real, qlen, rlen) = tc.draw(sam_lines::cigar(HTSLIB, 6));
    let mut ops = Vec::new();
    parse_cigar(real.as_bytes(), &mut ops).unwrap();
    let qlen = qlen as usize;
    let cg_type = tc.draw(gs::sampled_from(&["I", "i"]));
    let mut cg = String::new();
    for op in &ops {
        let v = op.to_bam_u32();
        cg.push(',');
        cg.push_str(&if cg_type == "I" { v.to_string() } else { v.cast_signed().to_string() });
    }
    // The placeholder may also be a CIGAR of its own that happens to start
    // with `<qlen>S`, with fewer ops than the tag.
    let placeholder =
        if tc.draw(gs::booleans()) { format!("{qlen}S{}N", rlen) } else { format!("{qlen}S") };
    let seq = if qlen == 0 { "*".to_string() } else { "A".repeat(qlen) };
    let before = tc.draw(gs::booleans());
    let mut fields: Vec<String> =
        ["cg", "0", "chr2", "10", "60", &placeholder, "*", "0", "0", &seq, "*"]
            .map(str::to_string)
            .to_vec();
    if before {
        fields.push("XA:i:1".to_string());
    }
    fields.push(format!("CG:B:{cg_type}{cg}"));
    fields.push("XB:Z:z".to_string());
    let lines = [Line { fields }];
    let sam = sam_text(&lines);
    let theirs = samtools_payloads(sam.as_bytes());
    let ours = seqair_payloads(sam.as_bytes());
    match (ours, theirs) {
        (Ok(ours), Ok(theirs)) => assert_same(&ours, &theirs, &lines, "samtools"),
        (Err(_), Err(_)) => tc.event("both reject"),
        (ours, theirs) => panic!("ours {ours:?}\ntheirs {theirs:?}\n{sam}"),
    }
}

/// Mutations of a canonical line, and which one was applied.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Mutation {
    ReplaceField(usize),
    AppendToField(usize),
    DropField(usize),
    ExtraTab,
    Truncate,
}

const JUNK: &str = "0123456789-+*=:,.;xXMSDIBNHPACGTUnZiIfdBcCsS \t";

#[hegel::composite]
fn mutated_line(tc: &TestCase) -> (Line, Mutation, String) {
    let line = tc.draw_silent(canonical_line(HTSLIB));
    let n = line.fields.len();
    let kind = tc.draw_silent(gs::integers::<u8>().max_value(4));
    let field = tc.draw_silent(gs::integers::<usize>().max_value(n - 1));
    let junk = tc.draw_silent(sam_lines::chars(JUNK, 0, 6));
    let mut fields = line.fields.clone();
    let mutation = match kind {
        0 => {
            fields[field] = junk;
            Mutation::ReplaceField(field)
        }
        1 => {
            fields[field].push_str(&junk);
            Mutation::AppendToField(field)
        }
        2 => {
            fields.remove(field);
            Mutation::DropField(field)
        }
        3 => {
            fields.insert(field, String::new());
            Mutation::ExtraTab
        }
        _ => Mutation::Truncate,
    };
    let mut text = Line { fields }.text();
    if mutation == Mutation::Truncate {
        let cut = tc.draw_silent(gs::integers::<usize>().max_value(text.len()));
        text.truncate(cut);
    }
    (line, mutation, text)
}

/// Whether htslib accepting a line that seqair rejects with `e` is one of
/// the documented deviations (`r[sam.encode.deviations]`), given the field
/// as written.
fn is_documented_deviation(e: &SamRecordError, fields: &[&str]) -> bool {
    let field = |i: usize| fields.get(i).copied().unwrap_or("");
    let empty_or_plus = |f: &str| f.is_empty() || f.starts_with('+');
    match e {
        SamRecordError::UnknownReference { .. } => true,
        SamRecordError::InvalidFlag { .. } => field(1).is_empty() || field(1).starts_with('0'),
        SamRecordError::InvalidPos { .. } => empty_or_plus(field(3)),
        SamRecordError::InvalidMapq { .. } => empty_or_plus(field(4)),
        SamRecordError::InvalidPnext { .. } => empty_or_plus(field(7)),
        SamRecordError::InvalidTlen { .. } => matches!(field(8), "" | "-" | "+"),
        SamRecordError::InvalidCigarLength { .. }
        | SamRecordError::UnknownCigarOp { .. }
        | SamRecordError::MissingCigarOp { .. }
        | SamRecordError::EmptyCigar => field(5).starts_with('*'),
        SamRecordError::InvalidSeq { .. } => field(9).bytes().any(|b| !(b'!'..=b'~').contains(&b)),
        SamRecordError::InvalidAuxField { .. }
        | SamRecordError::InvalidAuxValue { .. }
        | SamRecordError::InvalidAuxSubtype { .. } => fields.len() > 11,
        _ => false,
    }
}

// r[verify sam.encode]
// r[verify sam.encode.deviations]
// r[verify sam.encode.fields]
// r[verify sam.encode.numbers]
// r[verify sam.aux.parse]
/// Mutated lines: whatever htslib rejects, seqair rejects; whatever seqair
/// accepts, htslib accepts with the same bytes; and where only htslib
/// accepts, it is a documented deviation.
#[hegel::test(test_cases = crate::pinned::cases(300))]
fn malformed_lines_are_rejected_wherever_htslib_rejects_them(tc: TestCase) {
    let (_, mutation, text) = tc.draw(mutated_line());
    let sam = format!("{}{text}\n", header_text());
    let theirs = samtools_payloads(sam.as_bytes());
    let mut record = RawBamRecord::default();
    let ours = SamRecordEncoder::new().encode(text.as_bytes(), &header(), &mut record);
    match (ours, theirs) {
        (Ok(()), Ok(theirs)) => {
            assert_eq!(theirs.len(), 1, "{mutation:?} {text:?}");
            assert_eq!(record.as_payload(), theirs[0].as_slice(), "{mutation:?} {text:?}");
            tc.event("both accept");
        }
        (Err(_), Err(_)) => tc.event("both reject"),
        (Ok(()), Err(e)) => {
            panic!("seqair accepts what htslib rejects: {mutation:?} {text:?}\n{e}")
        }
        (Err(e), Ok(_)) => {
            let fields: Vec<&str> = text.split('\t').collect();
            assert!(
                is_documented_deviation(&e, &fields),
                "seqair rejects what htslib accepts, not a documented deviation: {mutation:?} {text:?}\n{e:?}"
            );
            tc.event("documented deviation");
        }
    }
}

// r[verify sam.encode.deviations]
// r[verify sam.encode.flag]
// r[verify sam.encode.numbers]
// r[verify sam.encode.references+2]
// r[verify sam.encode.cigar]
// r[verify sam.encode.seq]
// r[verify sam.aux.parse]
// r[verify sam.aux.types]
/// Each documented deviation: htslib accepts the line (losing or misreading
/// part of it), seqair rejects it with the named error.
#[test]
fn each_deviation_htslib_accepts_and_seqair_rejects() {
    type Check = fn(&SamRecordError) -> bool;
    let base = |f: &[(usize, &str)]| {
        let mut fields: Vec<String> =
            ["r", "0", "chr1", "10", "0", "4M", "*", "0", "0", "ACGT", "IIII"]
                .map(str::to_string)
                .to_vec();
        for &(i, v) in f {
            if i < fields.len() {
                fields[i] = v.to_string();
            } else {
                fields.push(v.to_string());
            }
        }
        fields.join("\t")
    };
    let cases: Vec<(String, Check)> = vec![
        (base(&[(1, "016")]), |e| matches!(e, SamRecordError::InvalidFlag { .. })),
        (base(&[(1, "0x10")]), |e| matches!(e, SamRecordError::InvalidFlag { .. })),
        (base(&[(1, "")]), |e| matches!(e, SamRecordError::InvalidFlag { .. })),
        (base(&[(2, "chrX")]), |e| matches!(e, SamRecordError::UnknownReference { .. })),
        (base(&[(3, "")]), |e| matches!(e, SamRecordError::InvalidPos { .. })),
        (base(&[(3, "+10")]), |e| matches!(e, SamRecordError::InvalidPos { .. })),
        (base(&[(4, "+1")]), |e| matches!(e, SamRecordError::InvalidMapq { .. })),
        (base(&[(7, "")]), |e| matches!(e, SamRecordError::InvalidPnext { .. })),
        (base(&[(8, "-")]), |e| matches!(e, SamRecordError::InvalidTlen { .. })),
        (base(&[(5, "*junk")]), |e| matches!(e, SamRecordError::UnknownCigarOp { op: b'*' })),
        (base(&[(9, "AC T")]), |e| matches!(e, SamRecordError::InvalidSeq { byte: b' ', .. })),
        (base(&[(11, "NMxiy5")]), |e| matches!(e, SamRecordError::InvalidAuxField { .. })),
        (base(&[(11, "NM:i:5abc")]), |e| matches!(e, SamRecordError::InvalidAuxValue { .. })),
        (base(&[(11, "NM:i:abc")]), |e| matches!(e, SamRecordError::InvalidAuxValue { .. })),
        (base(&[(11, "XA:A:ab")]), |e| matches!(e, SamRecordError::InvalidAuxValue { .. })),
        (base(&[(11, "XF:f: 1.5")]), |e| matches!(e, SamRecordError::InvalidAuxValue { .. })),
        (base(&[(11, "XF:f:0x1p3")]), |e| matches!(e, SamRecordError::InvalidAuxValue { .. })),
        (base(&[(11, "XH:H:ZZ")]), |e| matches!(e, SamRecordError::InvalidAuxValue { .. })),
        (base(&[(11, "XB:B:c,")]), |e| matches!(e, SamRecordError::InvalidAuxValue { .. })),
        (base(&[(11, "XB:B:c,200")]), |e| matches!(e, SamRecordError::InvalidAuxValue { .. })),
        (base(&[(11, "XB:B:C,-1")]), |e| matches!(e, SamRecordError::InvalidAuxValue { .. })),
        (base(&[(11, "XB:B:A")]), |e| {
            matches!(e, SamRecordError::InvalidAuxSubtype { subtype: b'A', .. })
        }),
    ];
    for (line, check) in cases {
        let sam = format!("{}{line}\n", header_text());
        assert!(samtools_payloads(sam.as_bytes()).is_ok(), "htslib should accept {line:?}");
        let err = encode(&line).expect_err(&line);
        assert!(check(&err), "{line:?}: {err:?}");
    }
}

// r[verify sam.encode.references+2]
/// An RNEXT the header lacks is no mate reference, PNEXT and FLAG kept: the
/// bytes `samtools view -b` writes.
#[test]
fn an_unknown_mate_reference_is_dropped_as_htslib_drops_it() {
    for line in [
        "r\t1\tchr1\t10\t0\t4M\tchrX\t5\t0\tACGT\tIIII",
        "r\t33\tchr1\t10\t0\t4M\tchrX\t0\t-7\tACGT\t*",
    ] {
        let sam = format!("{}{line}\n", header_text());
        let theirs = samtools_payloads(sam.as_bytes()).unwrap();
        assert_eq!(theirs.len(), 1);
        assert_eq!(encode(line).unwrap(), theirs[0], "{line:?}");
    }
}

// r[verify sam.encode.limits]
#[test]
fn a_record_over_the_size_limit_is_rejected() {
    let n = 1_500_000;
    let line = format!("big\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}", "A".repeat(n), "I".repeat(n));
    assert!(matches!(encode(&line), Err(SamRecordError::RecordTooLarge { .. })));
}

// r[verify sam.encode.fields]
#[test]
fn qname_limits_and_nul_bytes() {
    let line = |q: &str| format!("{q}\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*");
    let empty = encode(&line("")).unwrap();
    assert_eq!(empty[8], 1, "an empty QNAME is just its NUL");
    encode(&line(&"q".repeat(254))).unwrap();
    assert!(matches!(
        encode(&line(&"q".repeat(255))),
        Err(SamRecordError::QnameTooLong { len: 255 })
    ));
    assert!(matches!(encode("r\0\t4"), Err(SamRecordError::NulByte { offset: 1 })));
    assert!(matches!(encode(""), Err(SamRecordError::EmptyLine)));
    assert!(matches!(encode("r\t4\t*"), Err(SamRecordError::TooFewFields { found: 3 })));
}

// r[verify sam.encode.cigar_seq]
/// The CIGAR must cover SEQ even on an unmapped read, as `sam_parse1` checks.
#[test]
fn an_unmapped_read_still_needs_a_matching_cigar() {
    let line = "r\t4\tchr1\t10\t0\t3M\t*\t0\t0\tACGT\t*";
    let sam = format!("{}{line}\n", header_text());
    assert!(samtools_payloads(sam.as_bytes()).is_err());
    assert!(matches!(
        encode(line),
        Err(SamRecordError::CigarSeqLenMismatch { cigar_query_len: 3, seq_len: 4 })
    ));
}

/// The bytes of `sam_parse_cigar` (htslib 1.19.1, via rust-htslib's FFI),
/// or `None` where it fails.
fn htslib_parse_cigar(text: &str) -> Option<Vec<u32>> {
    use std::ffi::CString;
    let input = CString::new(text).ok()?;
    let mut buf: *mut u32 = std::ptr::null_mut();
    let mut mem: usize = 0;
    let mut end: *mut std::ffi::c_char = std::ptr::null_mut();
    // SAFETY: `input` is a NUL-terminated string that outlives the call;
    // htslib allocates `buf` with `realloc` and sets `mem` to its capacity.
    let n = unsafe {
        rust_htslib::htslib::sam_parse_cigar(
            input.as_ptr(),
            &raw mut end,
            &raw mut buf,
            &raw mut mem,
        )
    };
    let ops = if n >= 0 && !buf.is_null() {
        // SAFETY: htslib wrote `n` ops into `buf`, which holds `mem >= n`.
        Some(unsafe { std::slice::from_raw_parts(buf, n as usize) }.to_vec())
    } else if n == 0 {
        Some(vec![])
    } else {
        None
    };
    // SAFETY: `end` points into `input`, which is alive; consumed is the
    // distance htslib parsed.
    let consumed = if end.is_null() { 0 } else { unsafe { end.offset_from(input.as_ptr()) } };
    // SAFETY: `buf` came from htslib's `realloc` (or is null).
    unsafe { libc::free(buf.cast()) };
    // htslib stops at the first byte it cannot use; a whole CIGAR is one it
    // parsed to the end.
    (consumed as usize == text.len() && (n > 0 || text == "*")).then_some(ops?)
}

#[hegel::composite]
fn cigar_text(tc: &TestCase) -> String {
    if tc.draw_silent(gs::booleans()) {
        tc.draw_silent(sam_lines::cigar(HTSLIB, 8)).0
    } else {
        tc.draw_silent(sam_lines::chars("0123456789MIDNSHP=XB*+QZ", 0, 12))
    }
}

// r[verify sam.encode.cigar]
/// `parse_cigar` against htslib's `sam_parse_cigar`: both accept the same
/// strings, as the same ops.
#[hegel::test(test_cases = crate::pinned::cases(500))]
fn parse_cigar_matches_htslib(tc: TestCase) {
    let text = tc.draw(cigar_text());
    let mut ops = Vec::new();
    let ours = parse_cigar(text.as_bytes(), &mut ops)
        .map(|()| ops.iter().map(|o| o.to_bam_u32()).collect::<Vec<_>>());
    let theirs = htslib_parse_cigar(&text);
    match (ours, theirs) {
        (Ok(ours), Some(theirs)) => assert_eq!(ours, theirs, "{text:?}"),
        (Err(_), None) => {}
        (Ok(ours), None) => panic!("seqair accepts {text:?} as {ours:?}, htslib rejects it"),
        (Err(e), Some(_)) => panic!("seqair rejects {text:?} ({e:?}), htslib accepts it"),
    }
}

// r[verify sam.encode.cigar]
/// `parse_cigar` and `CigarStr` are inverse on `MIDNSHP=X` CIGARs.
#[hegel::test]
fn parse_cigar_round_trips_through_cigar_str(tc: TestCase) {
    let n = tc.draw(gs::integers::<usize>().max_value(20));
    let ops: Vec<CigarOp> = (0..n)
        .map(|_| {
            let len = tc.draw(gs::integers::<u32>().max_value((1 << 28) - 1));
            let code = tc.draw(gs::integers::<u32>().max_value(8));
            CigarOp::from_bam_u32((len << 4) | code)
        })
        .collect();
    let text = if ops.is_empty() { "*".to_string() } else { CigarStr(&ops).to_string() };
    let mut parsed = vec![CigarOp::from_bam_u32(0)];
    parse_cigar(text.as_bytes(), &mut parsed).unwrap();
    assert_eq!(parsed, ops, "{text}");
}

/// What a SAM sequence byte decodes to, written out by hand from the SAM
/// spec and htslib's table (not seqair's): IUPAC upper-cased, `U` is `T`,
/// colour-space digits are bases, anything else `N`.
fn expected_base(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'U' => b'T',
        b'0' => b'A',
        b'1' => b'C',
        b'2' => b'G',
        b'3' => b'T',
        c @ (b'=' | b'A' | b'C' | b'M' | b'G' | b'R' | b'S' | b'V' | b'T' | b'W' | b'Y' | b'H'
        | b'K' | b'D' | b'B' | b'N') => c,
        _ => b'N',
    }
}

// r[verify sam.encode]
// r[verify sam.encode.fixups]
// r[verify sam.encode.seq]
// r[verify sam.aux.types]
/// noodles' SAM reader agrees field by field with what seqair encoded, once
/// htslib's fix-ups (which noodles does not apply) are applied to its view.
#[hegel::test(test_cases = crate::pinned::cases(150))]
fn noodles_agrees_field_by_field(tc: TestCase) {
    use noodles::sam;
    use noodles::sam::alignment::record::data::field::Value;

    let lines = tc.draw(canonical_lines(NOODLES, 6));
    let text = sam_text(&lines);
    let ours = seqair_payloads(text.as_bytes()).unwrap();
    let our_header = header();
    let mut reader = sam::io::Reader::new(text.as_bytes());
    let n_header = reader.read_header().unwrap();
    for ((record, payload), line) in reader.records().zip(ours).zip(&lines) {
        let r = record.unwrap();
        let ctx = line.text();
        let ours = RawBamRecord::from_payload(payload, &our_header).unwrap();
        // noodles reads QNAME `*` as no name; BAM stores the `*`.
        assert_eq!(ours.qname(), r.name().map_or(&b"*"[..], <&[u8]>::from), "{ctx}");

        let n_flags = r.flags().unwrap().bits();
        let n_tid = r.reference_sequence_id(&n_header).transpose().unwrap();
        let n_pos = r.alignment_start().transpose().unwrap().map(|p| usize::from(p) - 1);
        let cigar_star = line.fields[5] == "*";
        // htslib's fix-ups, applied to noodles' literal reading.
        let tid = if n_pos.is_none() { None } else { n_tid };
        let unmapped = n_flags & 4 != 0 || tid.is_none() || cigar_star;
        // noodles keeps only the 12 defined FLAG bits.
        assert_eq!(
            ours.flags().raw() & 0xFFF,
            (n_flags & !4) | if unmapped { 4 } else { 0 },
            "{ctx}"
        );
        assert_eq!(ours.tid().map(|t| t.as_usize()), tid, "{ctx}");
        assert_eq!(ours.pos().map(|p| p.as_u64() as usize), n_pos, "{ctx}");
        let n_mapq = r.mapping_quality().transpose().unwrap().map_or(255, u8::from);
        assert_eq!(ours.mapq().as_byte(), n_mapq, "{ctx}");

        let n_cigar: Vec<(u8, usize)> = r
            .cigar()
            .iter()
            .map(|op| {
                let op = op.unwrap();
                (op.kind() as u8, op.len())
            })
            .collect();
        let our_cigar: Vec<(u8, usize)> =
            ours.cigar().map(|op| (op.op_code(), op.len() as usize)).collect();
        assert_eq!(our_cigar, n_cigar, "{ctx}");

        let mut seq = Vec::new();
        ours.decode_seq_into(&mut seq);
        let want: Vec<u8> = r.sequence().as_ref().iter().map(|&b| expected_base(b)).collect();
        assert_eq!(seq, want, "{ctx}");
        let n_qual: Vec<u8> = r.quality_scores().as_ref().iter().map(|q| q - 33).collect();
        let our_qual: Vec<u8> =
            ours.qual().map(|q| q.iter().map(|q| q.as_byte()).collect()).unwrap_or_default();
        assert_eq!(our_qual, n_qual, "{ctx}");

        let mate = ours.mate_fields();
        let n_mpos = r.mate_alignment_start().transpose().unwrap().map(|p| usize::from(p) - 1);
        let n_mtid = if n_mpos.is_none() {
            None
        } else if line.fields[6] == "=" {
            tid
        } else {
            r.mate_reference_sequence_id(&n_header).transpose().unwrap()
        };
        assert_eq!(mate.tid.map(|t| t.as_usize()), n_mtid, "{ctx}");
        assert_eq!(mate.pos.map(|p| p.as_u64() as usize), n_mpos, "{ctx}");
        assert_eq!(mate.template_len.map_or(0, |t| t.get()), r.template_length().unwrap(), "{ctx}");

        let our_aux: Vec<_> = ours.aux().iter_checked().map(Result::unwrap).collect();
        let n_aux: Vec<_> = r.data().iter().map(Result::unwrap).collect();
        assert_eq!(our_aux.len(), n_aux.len(), "{ctx}");
        for ((tag, ours), (n_tag, theirs)) in our_aux.iter().zip(&n_aux) {
            assert_eq!(tag, n_tag.as_ref(), "{ctx}");
            use seqair::bam::aux::AuxValue as A;
            let int = |v: &A<'_>| v.as_i64();
            let same = match (ours, theirs) {
                (A::Char(a), Value::Character(b)) => a == b,
                (o, Value::Int8(b)) => int(o) == Some(i64::from(*b)),
                (o, Value::UInt8(b)) => int(o) == Some(i64::from(*b)),
                (o, Value::Int16(b)) => int(o) == Some(i64::from(*b)),
                (o, Value::UInt16(b)) => int(o) == Some(i64::from(*b)),
                (o, Value::Int32(b)) => int(o) == Some(i64::from(*b)),
                (o, Value::UInt32(b)) => int(o) == Some(i64::from(*b)),
                (A::Float(a), Value::Float(b)) => a.to_bits() == b.to_bits(),
                (A::String(a), Value::String(b)) => *a == <&[u8]>::from(*b),
                (A::Hex(a), Value::Hex(b)) => *a == <&[u8]>::from(*b),
                (o, Value::Array(arr)) => our_array(o) == Some(noodles_array(arr)),
                _ => false,
            };
            assert!(same, "{ctx}: {tag:?} ours {ours:?} noodles {theirs:?}");
        }
    }
}

/// A B array as (subtype, elements as f64 bits).
type Elements = (char, Vec<u64>);

fn our_array(v: &seqair::bam::aux::AuxValue<'_>) -> Option<Elements> {
    use seqair::bam::aux::AuxValue as A;
    fn each<const N: usize>(b: &[u8], f: impl Fn([u8; N]) -> f64) -> Vec<u64> {
        b.as_chunks::<N>().0.iter().map(|c| f(*c).to_bits()).collect()
    }
    Some(match v {
        A::ArrayI8(b) => ('c', each::<1>(b, |c| f64::from(i8::from_le_bytes(c)))),
        A::ArrayU8(b) => ('C', each::<1>(b, |c| f64::from(c[0]))),
        A::ArrayI16(b) => ('s', each::<2>(b, |c| f64::from(i16::from_le_bytes(c)))),
        A::ArrayU16(b) => ('S', each::<2>(b, |c| f64::from(u16::from_le_bytes(c)))),
        A::ArrayI32(b) => ('i', each::<4>(b, |c| f64::from(i32::from_le_bytes(c)))),
        A::ArrayU32(b) => ('I', each::<4>(b, |c| f64::from(u32::from_le_bytes(c)))),
        A::ArrayFloat(b) => ('f', each::<4>(b, |c| f64::from(f32::from_le_bytes(c)))),
        _ => return None,
    })
}

fn noodles_array(arr: &noodles::sam::alignment::record::data::field::value::Array<'_>) -> Elements {
    use noodles::sam::alignment::record::data::field::value::Array;
    fn each<T: Into<f64>>(
        v: &(dyn noodles::sam::alignment::record::data::field::value::array::Values<'_, T> + '_),
    ) -> Vec<u64> {
        v.iter().map(|x| x.unwrap().into().to_bits()).collect()
    }
    match arr {
        Array::Int8(v) => ('c', each(v.as_ref())),
        Array::UInt8(v) => ('C', each(v.as_ref())),
        Array::Int16(v) => ('s', each(v.as_ref())),
        Array::UInt16(v) => ('S', each(v.as_ref())),
        Array::Int32(v) => ('i', each(v.as_ref())),
        Array::UInt32(v) => ('I', each(v.as_ref())),
        Array::Float(v) => ('f', each(v.as_ref())),
    }
}
