//! Hegel generators for BAM record payloads, built directly in wire form.
//!
//! Every generated payload passes `RawBamRecord::from_payload` against
//! [`header`]: the generator constructs valid shapes instead of rejecting
//! invalid ones. It covers unplaced, placed-unmapped and mapped records, every
//! 4-bit base code, odd sequence lengths (with junk in the unused low nibble),
//! missing SEQ and QUAL, every aux type and `B` subtype (empty arrays
//! included), NaN/infinite/subnormal/huge floats, mate fields, every CIGAR op
//! code, and long CIGARs in `CG:B:I`/`CG:B:i` tags — real ones over 65535
//! operations included.
//!
//! Include it with `#[path = "support/raw_record_gen.rs"] mod raw_record_gen;`.
#![allow(
    dead_code,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss,
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    reason = "test generator; each including test uses a different part"
)]

use hegel::TestCase;
use hegel::generators as gs;
use seqair::bam::BamHeader;
use std::fmt::Write as _;

/// The targets of [`header`]: one long enough for any position BAM can hold.
pub const TARGETS: [(&str, u32); 3] = [("chr1", i32::MAX as u32), ("chr2", 1000), ("chrM", 16569)];

pub fn header_text() -> String {
    let mut text = String::from("@HD\tVN:1.6\tSO:unsorted\n");
    for (name, len) in TARGETS {
        writeln!(text, "@SQ\tSN:{name}\tLN:{len}").expect("writing to a String cannot fail");
    }
    text
}

pub fn header() -> BamHeader {
    BamHeader::from_sam_text(&header_text()).expect("the generator's header is valid")
}

/// Which corners of BAM a test can take.
#[derive(Debug, Clone, Copy)]
pub struct Profile {
    /// What BAM can hold and htslib prints, but the SAM spec (and noodles)
    /// rejects: any non-NUL byte in names and `A`/`Z`/`H` values, a NUL inside
    /// a read name, `@` in one, FLAG bits above 0x800, position `i32::MAX`, quality bytes above 93, reserved CIGAR op codes
    /// 9–15, an unmapped record's CIGAR that does not cover its sequence,
    /// NaN and infinite floats, and `d` (double) tags.
    pub exotic: bool,
    /// NaNs with their sign bit set, which glibc prints `-nan` and macOS's
    /// libc `nan`.
    pub negative_nan: bool,
}

/// Everything htslib reads; negative NaNs only where libc prints them as
/// seqair does.
pub const HTSLIB: Profile = Profile { exotic: true, negative_nan: !cfg!(target_os = "macos") };

/// What noodles reads and writes canonically.
pub const NOODLES: Profile = Profile { exotic: false, negative_nan: false };

/// One aux tag's type code and value bytes (everything after the key).
#[derive(Debug, Clone, PartialEq)]
pub struct GenTag {
    pub key: [u8; 2],
    pub value: Vec<u8>,
}

impl GenTag {
    pub fn encode(&self, out: &mut Vec<u8>) {
        out.extend_from_slice(&self.key);
        out.extend_from_slice(&self.value);
    }
}

/// A record before it becomes payload bytes, kept structured so a failing
/// case prints readably.
#[derive(Debug, Clone)]
pub struct GenRecord {
    pub refid: i32,
    pub pos: i32,
    pub mapq: u8,
    pub bin: u16,
    pub flags: u16,
    pub next_refid: i32,
    pub next_pos: i32,
    pub tlen: i32,
    /// Without the NUL terminator.
    pub qname: Vec<u8>,
    /// Packed `len << 4 | op`.
    pub cigar: Vec<u32>,
    /// One 4-bit code per base.
    pub seq: Vec<u8>,
    /// The unused low nibble of an odd-length sequence's last byte.
    pub seq_pad: u8,
    pub qual: Vec<u8>,
    pub tags: Vec<GenTag>,
}

impl GenRecord {
    /// The payload: the record's bytes without `block_size`.
    pub fn payload(&self) -> Vec<u8> {
        let mut p = Vec::new();
        p.extend_from_slice(&self.refid.to_le_bytes());
        p.extend_from_slice(&self.pos.to_le_bytes());
        p.push(u8::try_from(self.qname.len() + 1).expect("names are at most 254 bytes"));
        p.push(self.mapq);
        p.extend_from_slice(&self.bin.to_le_bytes());
        p.extend_from_slice(&u16::try_from(self.cigar.len()).expect("≤ 65535 ops").to_le_bytes());
        p.extend_from_slice(&self.flags.to_le_bytes());
        p.extend_from_slice(&(self.seq.len() as u32).to_le_bytes());
        p.extend_from_slice(&self.next_refid.to_le_bytes());
        p.extend_from_slice(&self.next_pos.to_le_bytes());
        p.extend_from_slice(&self.tlen.to_le_bytes());
        p.extend_from_slice(&self.qname);
        p.push(0);
        for op in &self.cigar {
            p.extend_from_slice(&op.to_le_bytes());
        }
        for pair in self.seq.chunks(2) {
            p.push(pair[0] << 4 | pair.get(1).copied().unwrap_or(self.seq_pad));
        }
        p.extend_from_slice(&self.qual);
        for tag in &self.tags {
            tag.encode(&mut p);
        }
        p
    }
}

const OP_M: u32 = 0;
const OP_I: u32 = 1;
const OP_D: u32 = 2;
const OP_N: u32 = 3;
const OP_S: u32 = 4;
const OP_H: u32 = 5;
const OP_P: u32 = 6;
const OP_EQ: u32 = 7;
const OP_X: u32 = 8;
const OP_B: u32 = 9;

pub fn op(len: u32, code: u32) -> u32 {
    len << 4 | code
}

fn n<T: hegel::generators::Integer>(lo: T, hi: T) -> gs::IntegerGenerator<T> {
    gs::integers::<T>().min_value(lo).max_value(hi)
}

/// An op length: mostly short, sometimes up to BAM's 28-bit limit.
fn op_len(tc: &TestCase) -> u32 {
    if tc.draw_silent(gs::weighted_booleans(0.05)) {
        tc.draw_silent(n(1, (1 << 28) - 1))
    } else {
        tc.draw_silent(n(1, 60))
    }
}

/// A run of ops that consume no query bases.
fn non_query_ops(tc: &TestCase, profile: Profile, out: &mut Vec<u32>) {
    let codes: &[u32] = if profile.exotic {
        &[OP_D, OP_N, OP_H, OP_P, OP_B, 10, 11, 12, 13, 14, 15]
    } else {
        &[OP_D, OP_N, OP_H, OP_P]
    };
    for _ in 0..tc.draw_silent(n(0usize, 2)) {
        let code = tc.draw_silent(gs::sampled_from(codes));
        out.push(op(op_len(tc), code));
    }
}

/// A CIGAR whose query-consuming ops sum to exactly `seq_len` (any ops at all
/// when `seq_len` is 0), never starting with `<seq_len>S` so it cannot be
/// mistaken for a long-CIGAR placeholder.
pub fn cigar_for(tc: &TestCase, profile: Profile, seq_len: u32) -> Vec<u32> {
    let mut ops = Vec::new();
    let query_ops = [OP_M, OP_I, OP_S, OP_EQ, OP_X];
    let mut left = seq_len;
    if seq_len == 0 {
        // htslib checks nothing here; zero-length query ops are legal too.
        for _ in 0..tc.draw_silent(n(0usize, 3)) {
            let code = tc.draw_silent(gs::sampled_from(&query_ops[..]));
            ops.push(op(0, code));
            non_query_ops(tc, profile, &mut ops);
        }
        return ops;
    }
    non_query_ops(tc, profile, &mut ops);
    while left > 0 {
        let len = tc.draw_silent(n(1, left));
        let code = tc.draw_silent(gs::sampled_from(&query_ops[..]));
        // `<seq_len>S` first would be a placeholder; make it a match.
        let code = if ops.is_empty() && code == OP_S && len == seq_len { OP_M } else { code };
        ops.push(op(len, code));
        left -= len;
        non_query_ops(tc, profile, &mut ops);
    }
    ops
}

/// A CIGAR of any ops, for unmapped records (htslib checks nothing).
fn any_cigar(tc: &TestCase, profile: Profile, seq_len: u32) -> Vec<u32> {
    let mut ops = Vec::new();
    let max_code = if profile.exotic { 15 } else { 8 };
    for _ in 0..tc.draw_silent(n(0usize, 4)) {
        let code = tc.draw_silent(n(0, max_code));
        let len = op_len(tc);
        // Keep clear of the placeholder shape.
        let len = if ops.is_empty() && code == OP_S && len == seq_len { len + 1 } else { len };
        ops.push(op(len, code));
    }
    ops
}

fn special_f64() -> Vec<f64> {
    vec![
        0.0,
        -0.0,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NAN,
        -f64::NAN,
        f64::MIN_POSITIVE,
        f64::from_bits(1),
        f64::MAX,
        f64::MIN,
        f64::EPSILON,
        1e-5,
        0.0001,
        0.000_099_999_99,
        0.001,
        0.009_999_999,
        0.099_999_95,
        0.999_999_5,
        1.0,
        9.999_995,
        99_999.95,
        99_999.97,
        123_456.5,
        999_998.5,
        999_999.0,
        999_999.4,
        999_999.5,
        1e6,
        1_234_565.0,
        1e15,
        1e16,
        1e100,
        1e-100,
        0.5,
        0.25,
        2.5e-5,
        -123.456,
        f64::from(f32::MAX),
        f64::from(f32::MIN_POSITIVE),
        f64::from(f32::from_bits(1)),
    ]
}

/// A double from every corner: raw bit patterns, the special list, and
/// decimals with seven significant digits (near-ties for six-digit rounding).
pub fn arb_f64(tc: &TestCase, negative_nan: bool) -> f64 {
    let v = match tc.draw_silent(n(0u8, 3)) {
        0 => f64::from_bits(tc.draw_silent(gs::integers::<u64>())),
        1 => tc.draw_silent(gs::sampled_from(special_f64())),
        2 => {
            let mantissa = tc.draw_silent(n(-10_000_000i64, 10_000_000));
            let exp = tc.draw_silent(n(-15i32, 15));
            mantissa as f64 * 10f64.powi(exp)
        }
        _ => tc.draw_silent(gs::floats::<f64>()),
    };
    if v.is_nan() && !negative_nan { f64::from_bits(v.to_bits() & !(1 << 63)) } else { v }
}

/// A float from the same corners, in `f32`.
pub fn arb_f32(tc: &TestCase, negative_nan: bool) -> f32 {
    let v = match tc.draw_silent(n(0u8, 3)) {
        0 => f32::from_bits(tc.draw_silent(gs::integers::<u32>())),
        1 => tc.draw_silent(gs::sampled_from(special_f64())) as f32,
        2 => {
            let mantissa = tc.draw_silent(n(-10_000_000i64, 10_000_000));
            let exp = tc.draw_silent(n(-15i32, 15));
            (mantissa as f64 * 10f64.powi(exp)) as f32
        }
        _ => tc.draw_silent(gs::floats::<f32>()),
    };
    if v.is_nan() && !negative_nan { f32::from_bits(v.to_bits() & !(1 << 31)) } else { v }
}

/// A float for a tag: finite unless `exotic`.
fn arb_f32_for(tc: &TestCase, profile: Profile) -> f32 {
    let v = arb_f32(tc, profile.negative_nan);
    if profile.exotic || v.is_finite() {
        v
    } else if v.is_nan() {
        0.0
    } else {
        f32::MAX.copysign(v)
    }
}

/// Bytes for a string value: printable ASCII (space included), or with
/// `exotic` any byte but NUL.
fn arb_text(tc: &TestCase, profile: Profile, max: usize) -> Vec<u8> {
    if profile.exotic && tc.draw_silent(gs::weighted_booleans(0.3)) {
        tc.draw_silent(gs::vecs(n(1u8, 255)).max_size(max))
    } else {
        tc.draw_silent(gs::vecs(n(b' ', b'~')).max_size(max))
    }
}

/// A tag key, never `CG` (which only [`arb_record`] places, on purpose).
fn arb_key(tc: &TestCase) -> [u8; 2] {
    let first = tc.draw_silent(gs::sampled_from(
        b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz".to_vec(),
    ));
    let second = tc.draw_silent(gs::sampled_from(
        b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789".to_vec(),
    ));
    if [first, second] == *b"CG" { *b"CH" } else { [first, second] }
}

/// A `B` array of `count` elements of `subtype`.
pub fn array_tag(key: [u8; 2], subtype: u8, elements: &[u8], count: u32) -> GenTag {
    let mut value = vec![b'B', subtype];
    value.extend_from_slice(&count.to_le_bytes());
    value.extend_from_slice(elements);
    GenTag { key, value }
}

pub fn arb_tag(tc: &TestCase, profile: Profile) -> GenTag {
    let key = arb_key(tc);
    let kinds: &[u8] = if profile.exotic { b"AcCsSiIfdZHB" } else { b"AcCsSiIfZHB" };
    let kind = tc.draw_silent(gs::sampled_from(kinds));
    let mut value = vec![kind];
    match kind {
        b'A' => value.push(if profile.exotic {
            tc.draw_silent(gs::integers::<u8>())
        } else {
            tc.draw_silent(n(b'!', b'~'))
        }),
        b'c' | b'C' => value.push(tc.draw_silent(gs::integers::<u8>())),
        b's' | b'S' => {
            value.extend_from_slice(&tc.draw_silent(gs::integers::<u16>()).to_le_bytes())
        }
        b'i' | b'I' => {
            value.extend_from_slice(&tc.draw_silent(gs::integers::<u32>()).to_le_bytes())
        }
        b'f' => value.extend_from_slice(&arb_f32_for(tc, profile).to_le_bytes()),
        b'd' => value.extend_from_slice(&arb_f64(tc, profile.negative_nan).to_le_bytes()),
        b'Z' => {
            value.extend(arb_text(tc, profile, 20));
            value.push(0);
        }
        b'H' => {
            if profile.exotic && tc.draw_silent(gs::weighted_booleans(0.2)) {
                value.extend(arb_text(tc, profile, 9));
            } else {
                let pairs = tc.draw_silent(n(0usize, 8));
                for _ in 0..pairs * 2 {
                    value.push(tc.draw_silent(gs::sampled_from(b"0123456789ABCDEF".to_vec())));
                }
            }
            value.push(0);
        }
        _ => {
            let subtype = tc.draw_silent(gs::sampled_from(b"cCsSiIf".to_vec()));
            let count = tc.draw_silent(n(0u32, 8));
            let mut elements = Vec::new();
            for _ in 0..count {
                match subtype {
                    b'c' | b'C' => elements.push(tc.draw_silent(gs::integers::<u8>())),
                    b's' | b'S' => elements
                        .extend_from_slice(&tc.draw_silent(gs::integers::<u16>()).to_le_bytes()),
                    b'i' | b'I' => elements
                        .extend_from_slice(&tc.draw_silent(gs::integers::<u32>()).to_le_bytes()),
                    _ => elements.extend_from_slice(&arb_f32_for(tc, profile).to_le_bytes()),
                }
            }
            return array_tag(key, subtype, &elements, count);
        }
    }
    GenTag { key, value }
}

fn arb_qname(tc: &TestCase, profile: Profile) -> Vec<u8> {
    let len = tc.draw_silent(n(0usize, 254));
    let short = tc.draw_silent(gs::weighted_booleans(0.9));
    let len = if short { len % 24 } else { len };
    let (lo, hi) = if profile.exotic && tc.draw_silent(gs::weighted_booleans(0.2)) {
        (1u8, 255u8)
    } else {
        (b'!', b'~')
    };
    let mut name: Vec<u8> = (0..len).map(|_| tc.draw_silent(n(lo, hi))).collect();
    if !profile.exotic {
        // SAM's QNAME is `[!-?A-~]`: no `@`.
        name.iter_mut().filter(|b| **b == b'@').for_each(|b| *b = b'A');
    }
    if profile.exotic && name.len() > 1 && tc.draw_silent(gs::weighted_booleans(0.05)) {
        let at = tc.draw_silent(n(0, name.len() - 1));
        name[at] = 0;
    }
    if name.is_empty() && !profile.exotic {
        // noodles reads an empty name as missing and prints `*`.
        name.push(b'r');
    }
    name
}

fn arb_seq(tc: &TestCase, len: usize) -> Vec<u8> {
    (0..len).map(|_| tc.draw_silent(n(0u8, 15))).collect()
}

fn arb_qual(tc: &TestCase, profile: Profile, len: usize) -> Vec<u8> {
    if len == 0 {
        return Vec::new();
    }
    if tc.draw_silent(gs::weighted_booleans(0.2)) {
        return vec![0xFF; len];
    }
    let hi = if profile.exotic && tc.draw_silent(gs::weighted_booleans(0.2)) { 254 } else { 93 };
    let mut qual: Vec<u8> = (0..len).map(|_| tc.draw_silent(n(0u8, hi))).collect();
    if profile.exotic && len > 1 && tc.draw_silent(gs::weighted_booleans(0.05)) {
        // Only the first byte decides QUAL `*`; a later 0xFF prints as a space.
        let at = tc.draw_silent(n(1, len - 1));
        qual[at] = 0xFF;
    }
    qual
}

/// A placed position: mostly small, sometimes anywhere BAM can hold (SAM's
/// 1-based POS stops one short of BAM's largest).
fn arb_pos(tc: &TestCase, profile: Profile) -> i32 {
    if tc.draw_silent(gs::weighted_booleans(0.1)) {
        tc.draw_silent(n(0, if profile.exotic { i32::MAX } else { i32::MAX - 1 }))
    } else {
        tc.draw_silent(n(0, 5000))
    }
}

/// How a record is placed and aligned.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Shape {
    Unplaced,
    PlacedUnmapped,
    Mapped,
    LongCigar,
}

/// A record that `RawBamRecord::from_payload` accepts against [`header`].
#[hegel::composite]
pub fn arb_record(tc: &TestCase, profile: Profile) -> GenRecord {
    let shape = tc.draw(gs::sampled_from(vec![
        Shape::Unplaced,
        Shape::PlacedUnmapped,
        Shape::Mapped,
        Shape::Mapped,
        Shape::LongCigar,
    ]));
    let targets = TARGETS.len() as i32;
    let placed = shape != Shape::Unplaced;
    let refid = if placed { tc.draw_silent(n(0, targets - 1)) } else { -1 };
    let pos = if placed { arb_pos(tc, profile) } else { -1 };
    let unmapped = match shape {
        Shape::Unplaced | Shape::PlacedUnmapped => true,
        Shape::Mapped => false,
        Shape::LongCigar => tc.draw_silent(gs::weighted_booleans(0.2)),
    };
    // FLAG bits above 0x800 are undefined; noodles drops them.
    let defined = if profile.exotic { u16::MAX } else { 0x0FFF };
    let flags = tc.draw_silent(gs::integers::<u16>()) & defined & !4 | if unmapped { 4 } else { 0 };

    // A real long CIGAR now and then: over 65535 ops, as BAM needs a CG tag for.
    let huge = shape == Shape::LongCigar && tc.draw_silent(gs::weighted_booleans(0.03));
    let seq_len = if huge {
        tc.draw_silent(n(32_768usize, 33_000))
    } else if tc.draw_silent(gs::weighted_booleans(0.15)) {
        0
    } else {
        tc.draw_silent(n(1usize, 40))
    };
    let seq = arb_seq(tc, seq_len);
    let qual = arb_qual(tc, profile, seq_len);
    let mut tags: Vec<GenTag> =
        (0..tc.draw_silent(n(0usize, 6))).map(|_| arb_tag(tc, profile)).collect();

    let cigar = match shape {
        Shape::Unplaced | Shape::PlacedUnmapped => {
            if tc.draw_silent(gs::booleans()) {
                Vec::new()
            } else if profile.exotic {
                any_cigar(tc, profile, seq_len as u32)
            } else {
                cigar_for(tc, profile, seq_len as u32)
            }
        }
        Shape::Mapped => cigar_for(tc, profile, seq_len as u32),
        Shape::LongCigar => {
            let real = if huge {
                // 1M1D per base: twice as many ops as bases.
                (0..seq_len).flat_map(|_| [op(1, OP_M), op(1, OP_D)]).collect()
            } else {
                let mut real = cigar_for(tc, profile, seq_len as u32);
                while real.len() < 2 {
                    real.push(op(op_len(tc), OP_D));
                }
                real
            };
            let rlen: u32 = real
                .iter()
                .filter(|&&o| matches!(o & 0xF, 0 | 2 | 3 | 7 | 8))
                .fold(0u32, |a, &o| a.saturating_add(o >> 4))
                .clamp(1, (1 << 28) - 1);
            let placeholder = vec![op(seq_len as u32, OP_S), op(rlen, OP_N)];
            let subtype = tc.draw_silent(gs::sampled_from(b"Ii".to_vec()));
            let bytes: Vec<u8> = real.iter().flat_map(|o| o.to_le_bytes()).collect();
            let cg = array_tag(*b"CG", subtype, &bytes, real.len() as u32);
            let at = tc.draw_silent(n(0, tags.len()));
            tags.insert(at, cg);
            // Sometimes a second CG after it, which stays a tag.
            if tc.draw_silent(gs::weighted_booleans(0.1)) {
                let at = tc.draw_silent(n(at + 1, tags.len()));
                tags.insert(at, array_tag(*b"CG", b'I', &op(3, OP_M).to_le_bytes(), 1));
            }
            placeholder
        }
    };

    let next_refid = match tc.draw_silent(n(0u8, 3)) {
        0 => -1,
        1 if placed => refid,
        _ => tc.draw_silent(n(0, targets - 1)),
    };
    let next_pos =
        if tc.draw_silent(gs::weighted_booleans(0.2)) { -1 } else { arb_pos(tc, profile) };
    GenRecord {
        refid,
        pos,
        mapq: tc.draw_silent(gs::integers::<u8>()),
        bin: tc.draw_silent(gs::integers::<u16>()),
        flags,
        next_refid,
        next_pos,
        tlen: tc.draw_silent(gs::integers::<i32>()),
        qname: arb_qname(tc, profile),
        cigar,
        seq,
        seq_pad: tc.draw_silent(n(0u8, 15)),
        qual,
        tags,
    }
}
