//! Canonical SAM line generators and the htslib/noodles oracles the SAM
//! encoder tests compare against. Included with `#[path]`.
#![allow(dead_code, reason = "each test binary uses a different subset")]
#![allow(
    clippy::expect_used,
    clippy::unwrap_in_result,
    clippy::format_push_string,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    reason = "test code"
)]

use hegel::TestCase;
use hegel::generators as gs;
use seqair::bam::{RawBamRecord, stream::BamStreamReader};
use seqair::sam::{SamStreamError, SamStreamReader};
use std::path::Path;
use std::process::Command;

/// Targets of [`header_text`]: a 2^31 − 1 bp contig, so positions reach the
/// range where the 16-bit `bin` wraps, and a name with `*` and `:`.
pub const TARGETS: [(&str, u64); 3] =
    [("chr1", 2_147_483_647), ("chr2", 5000), ("HLA-A*01:01", 3000)];

pub fn header_text() -> String {
    let mut h = String::from("@HD\tVN:1.6\n");
    for (name, len) in TARGETS {
        h.push_str(&format!("@SQ\tSN:{name}\tLN:{len}\n"));
    }
    h
}

/// What a generated line may contain. noodles is stricter than htslib on a
/// few things that are still valid SAM for htslib (`B` CIGAR op, `d` tags,
/// type-letter synonyms, explicit `+`), so its dialect leaves them out.
#[derive(Debug, Clone, Copy)]
pub struct Dialect {
    pub noodles: bool,
}

/// One SAM line as fields: the 11 mandatory ones, then the optional ones.
#[derive(Debug, Clone)]
pub struct Line {
    pub fields: Vec<String>,
}

impl Line {
    pub fn text(&self) -> String {
        self.fields.join("\t")
    }
}

/// `n` characters from `alphabet`.
#[hegel::composite]
pub fn chars(tc: &TestCase, alphabet: &'static str, min: usize, max: usize) -> String {
    let n = tc.draw_silent(gs::integers::<usize>().min_value(min).max_value(max));
    tc.draw_silent(gs::text().alphabet(alphabet).min_size(n).max_size(n))
}

const QNAME_CHARS: &str = "!\"#$%&'()*+,-./0123456789:;<=>?ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";
const PRINTABLE: &str = " !\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";
const QUAL_CHARS: &str = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";
/// IUPAC in both cases, `=`, `U`/`u`, `.`, colour-space digits and a few
/// letters htslib packs as `N`.
pub const SEQ_CHARS: &str = "ACGTNacgtn=MRWSYKVHDBmrwsykvhdbUu.0123XEJ";

/// An integer for an `i` tag: mostly the type boundaries.
#[hegel::composite]
fn aux_int(tc: &TestCase) -> i64 {
    if tc.draw_silent(gs::booleans()) {
        tc.draw_silent(gs::sampled_from(&[
            0i64,
            1,
            127,
            128,
            255,
            256,
            32767,
            32768,
            65535,
            65536,
            2_147_483_647,
            2_147_483_648,
            4_294_967_295,
            -1,
            -127,
            -128,
            -129,
            -32768,
            -32769,
            -2_147_483_648,
        ]))
    } else {
        tc.draw_silent(gs::integers::<i64>().min_value(-2_147_483_648).max_value(4_294_967_295))
    }
}

/// A float as text the way people and `%g` write it.
#[hegel::composite]
fn float_text(tc: &TestCase, d: Dialect, double: bool) -> String {
    let style = tc.draw_silent(gs::integers::<u8>().max_value(if d.noodles { 2 } else { 4 }));
    let v: f64 = if double {
        tc.draw_silent(gs::floats::<f64>().allow_nan(false).allow_infinity(false))
    } else {
        f64::from(tc.draw_silent(gs::floats::<f32>().allow_nan(false).allow_infinity(false)))
    };
    match style {
        0 => format!("{v}"),
        1 => format!("{v:e}"),
        2 => format!("{:.3}", v / 1e30),
        3 => tc
            .draw_silent(gs::sampled_from(&[
                "inf", "-inf", "nan", "-nan", "Infinity", "NaN", "+1.5", ".5", "5.", "1E-3",
            ]))
            .to_string(),
        _ => format!("{v:+e}"),
    }
}

/// One optional field with tag `tag`.
#[hegel::composite]
fn aux_field(tc: &TestCase, tag: String, d: Dialect) -> String {
    let kind = tc.draw_silent(gs::integers::<u8>().max_value(6));
    match kind {
        0 => {
            let ty = if d.noodles {
                "A"
            } else {
                tc.draw_silent(gs::sampled_from(&["A", "a", "c", "C"]))
            };
            format!(
                "{tag}:{ty}:{}",
                tc.draw_silent(chars(
                    "!#$%&'()*+,-./0123456789:;<=>?@ABCXYZ[\\]^_`abcxyz{|}~",
                    1,
                    1
                ))
            )
        }
        1 => {
            let v = tc.draw_silent(aux_int());
            let ty = if d.noodles { "i" } else { tc.draw_silent(gs::sampled_from(&["i", "I"])) };
            let plus = !d.noodles && v >= 0 && tc.draw_silent(gs::weighted_booleans(0.2));
            format!("{tag}:{ty}:{}{v}", if plus { "+" } else { "" })
        }
        2 => format!("{tag}:f:{}", tc.draw_silent(float_text(d, false))),
        3 if !d.noodles => format!("{tag}:d:{}", tc.draw_silent(float_text(d, true))),
        3 | 4 => format!("{tag}:Z:{}", tc.draw_silent(chars(PRINTABLE, 0, 20))),
        5 => {
            let n = tc.draw_silent(gs::integers::<usize>().max_value(6));
            let digits = if d.noodles { "0123456789ABCDEF" } else { "0123456789ABCDEFabcdef" };
            format!("{tag}:H:{}", tc.draw_silent(chars(digits, 2 * n, 2 * n)))
        }
        _ => {
            let sub = tc.draw_silent(gs::sampled_from(&['c', 'C', 's', 'S', 'i', 'I', 'f']));
            // noodles cannot read an empty array; htslib can.
            let n = tc.draw_silent(
                gs::integers::<usize>().min_value(usize::from(d.noodles)).max_value(6),
            );
            let mut s = format!("{tag}:B:{sub}");
            for _ in 0..n {
                let e = match sub {
                    'c' => tc.draw_silent(gs::integers::<i8>()).to_string(),
                    'C' => tc.draw_silent(gs::integers::<u8>()).to_string(),
                    's' => tc.draw_silent(gs::integers::<i16>()).to_string(),
                    'S' => tc.draw_silent(gs::integers::<u16>()).to_string(),
                    'i' => tc.draw_silent(gs::integers::<i32>()).to_string(),
                    'I' => tc.draw_silent(gs::integers::<u32>()).to_string(),
                    _ => tc.draw_silent(float_text(d, false)),
                };
                s.push(',');
                s.push_str(&e);
            }
            s
        }
    }
}

/// A CIGAR of `1..=max_ops` operations, as text plus its query and
/// reference lengths.
#[hegel::composite]
pub fn cigar(tc: &TestCase, d: Dialect, max_ops: usize) -> (String, u64, u64) {
    let n = tc.draw_silent(gs::integers::<usize>().min_value(1).max_value(max_ops));
    let ops: &[char] = if d.noodles {
        &['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X']
    } else {
        &['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', 'B']
    };
    let (mut text, mut qlen, mut rlen) = (String::new(), 0u64, 0u64);
    for _ in 0..n {
        let len =
            tc.draw_silent(gs::integers::<u64>().min_value(u64::from(d.noodles)).max_value(60));
        let op = tc.draw_silent(gs::sampled_from(ops));
        if matches!(op, 'M' | 'I' | 'S' | '=' | 'X') {
            qlen += len;
        }
        if matches!(op, 'M' | 'D' | 'N' | '=' | 'X') {
            rlen += len;
        }
        text.push_str(&format!("{len}{op}"));
    }
    (text, qlen, rlen)
}

/// Where a record sits, including the shapes htslib fixes up.
#[derive(Debug, Clone, Copy)]
enum Placement {
    /// RNAME `*`, POS 0.
    Unplaced,
    /// RNAME `*` but a POS.
    StarWithPos,
    Mapped,
    PlacedUnmapped,
    /// A reference but POS 0: htslib drops the reference.
    RefWithoutPos,
    /// Mapped but CIGAR `*`: htslib sets FLAG 0x4.
    MappedWithoutCigar,
}

/// A canonical line: valid SAM that htslib accepts unchanged in meaning.
#[hegel::composite]
pub fn canonical_line(tc: &TestCase, d: Dialect) -> Line {
    let qname = if tc.draw_silent(gs::weighted_booleans(0.03)) {
        tc.draw_silent(chars(QNAME_CHARS, 254, 254))
    } else {
        tc.draw_silent(chars(QNAME_CHARS, 1, 12))
    };
    let placement = tc.draw_silent(gs::sampled_from(&[
        Placement::Unplaced,
        Placement::StarWithPos,
        Placement::Mapped,
        Placement::Mapped,
        Placement::Mapped,
        Placement::PlacedUnmapped,
        Placement::RefWithoutPos,
        Placement::MappedWithoutCigar,
    ]));
    let mut flag = tc.draw_silent(gs::integers::<u16>());
    let target = tc.draw_silent(gs::integers::<usize>().max_value(TARGETS.len() - 1));
    let (name, len) = TARGETS[target];
    let random_pos = |tc: &TestCase| {
        // `len + 1` = 2^31 for chr1: the largest POS BAM can hold.
        let max = if d.noodles { len } else { len + 1 };
        if tc.draw_silent(gs::weighted_booleans(0.2)) {
            max
        } else {
            tc.draw_silent(gs::integers::<u64>().min_value(1).max_value(max))
        }
    };
    let (rname, pos) = match placement {
        Placement::Unplaced => ("*".to_string(), 0),
        Placement::StarWithPos => ("*".to_string(), random_pos(tc)),
        Placement::RefWithoutPos => (name.to_string(), 0),
        _ => (name.to_string(), random_pos(tc)),
    };
    match placement {
        Placement::Mapped | Placement::MappedWithoutCigar => flag &= !0x4,
        Placement::PlacedUnmapped => flag |= 0x4,
        _ => {}
    }
    let with_cigar = match placement {
        Placement::Mapped => true,
        Placement::MappedWithoutCigar => false,
        _ => tc.draw_silent(gs::booleans()),
    };
    let (cigar_text, qlen) = if with_cigar {
        let (text, qlen, _) = tc.draw_silent(cigar(d, 6));
        (text, Some(qlen))
    } else {
        ("*".to_string(), None)
    };
    let seq_len = match qlen {
        Some(0) => None,
        Some(q) => tc.draw_silent(gs::weighted_booleans(0.85)).then_some(q as usize),
        None => tc
            .draw_silent(gs::booleans())
            .then(|| tc.draw_silent(gs::integers::<usize>().min_value(1).max_value(40))),
    };
    let (seq, qual) = match seq_len {
        None => ("*".to_string(), "*".to_string()),
        Some(n) => {
            let seq = tc.draw_silent(chars(SEQ_CHARS, n, n));
            let qual = if tc.draw_silent(gs::weighted_booleans(0.3)) {
                "*".to_string()
            } else {
                tc.draw_silent(chars(QUAL_CHARS, n, n))
            };
            (seq, qual)
        }
    };
    let rnext_choice = tc.draw_silent(gs::integers::<usize>().max_value(TARGETS.len() + 1));
    let rnext = match rnext_choice {
        0 => "*".to_string(),
        1 => "=".to_string(),
        i => TARGETS[i - 2].0.to_string(),
    };
    let pnext = if tc.draw_silent(gs::weighted_booleans(0.2)) {
        0
    } else {
        tc.draw_silent(gs::integers::<u64>().min_value(1).max_value(if d.noodles {
            2_147_483_647
        } else {
            2_147_483_648
        }))
    };
    let tlen = tc.draw_silent(gs::integers::<i32>());
    let tlen = if !d.noodles && tlen >= 0 && tc.draw_silent(gs::weighted_booleans(0.1)) {
        format!("+{tlen}")
    } else {
        tlen.to_string()
    };
    let mapq = tc.draw_silent(gs::integers::<u8>());

    let mut fields = vec![
        qname,
        flag.to_string(),
        rname,
        pos.to_string(),
        mapq.to_string(),
        cigar_text,
        rnext,
        pnext.to_string(),
        tlen,
        seq,
        qual,
    ];
    let n_aux = tc.draw_silent(gs::integers::<usize>().max_value(6));
    let first =
        tc.draw_silent(chars("ABDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz", n_aux, n_aux));
    // Distinct by construction: the second byte is the field's index.
    for (i, c) in first.chars().enumerate() {
        let tag = format!("{c}{i}");
        fields.push(tc.draw_silent(aux_field(tag, d)));
    }
    Line { fields }
}

/// `n` canonical lines.
#[hegel::composite]
pub fn canonical_lines(tc: &TestCase, d: Dialect, max: usize) -> Vec<Line> {
    let n = tc.draw_silent(gs::integers::<usize>().min_value(1).max_value(max));
    (0..n).map(|_| tc.draw_silent(canonical_line(d))).collect()
}

pub fn sam_text(lines: &[Line]) -> String {
    let mut s = header_text();
    for l in lines {
        s.push_str(&l.text());
        s.push('\n');
    }
    s
}

/// The payloads of every record of a BAM, read by seqair's streaming reader.
pub fn bam_payloads(bam: &[u8]) -> Vec<Vec<u8>> {
    let mut reader = BamStreamReader::new(bam).expect("BAM header");
    reader.records().map(|r| r.expect("BAM record").into_payload()).collect()
}

/// `samtools view -b` (htslib's `sam_parse1` + `bam_write1`): the payloads,
/// or samtools' stderr if it rejects the input.
pub fn samtools_payloads(sam: &[u8]) -> Result<Vec<Vec<u8>>, String> {
    let dir = tempfile::tempdir().expect("tempdir");
    let input = dir.path().join("in.sam");
    let output = dir.path().join("out.bam");
    std::fs::write(&input, sam).expect("write SAM");
    let out = Command::new("samtools")
        .args(["view", "--no-PG", "-b", "-o"])
        .arg(&output)
        .arg(&input)
        .output()
        .expect("samtools not found — required for the SAM encoder oracle");
    if !out.status.success() {
        return Err(String::from_utf8_lossy(&out.stderr).into_owned());
    }
    Ok(bam_payloads(&std::fs::read(&output).expect("read BAM")))
}

/// The same through rust-htslib (its bundled htslib: 1.19.1, which predates
/// `U` → `T`).
pub fn rust_htslib_payloads(sam: &[u8]) -> Result<Vec<Vec<u8>>, String> {
    use rust_htslib::bam::{self, Read};
    let dir = tempfile::tempdir().expect("tempdir");
    let input = dir.path().join("in.sam");
    let output = dir.path().join("out.bam");
    std::fs::write(&input, sam).expect("write SAM");
    rust_htslib_convert(&input, &output).map_err(|e| e.to_string())?;
    fn rust_htslib_convert(input: &Path, output: &Path) -> Result<(), rust_htslib::errors::Error> {
        let mut reader = bam::Reader::from_path(input)?;
        let header = bam::Header::from_template(reader.header());
        let mut writer = bam::Writer::from_path(output, &header, bam::Format::Bam)?;
        for record in reader.records() {
            writer.write(&record?)?;
        }
        Ok(())
    }
    Ok(bam_payloads(&std::fs::read(&output).expect("read BAM")))
}

/// seqair: every record of `sam` through [`SamStreamReader`].
pub fn seqair_payloads(sam: &[u8]) -> Result<Vec<Vec<u8>>, SamStreamError> {
    let mut reader = SamStreamReader::new(sam)?;
    reader.records().map(|r| r.map(RawBamRecord::into_payload)).collect()
}
