//! A reference decoder for tok3, the read-name tokeniser: the pseudocode of
//! `CRAMcodecs` §5 "Name tokenisation codec" (`DecodeNames`,
//! `DecodeTokenByteStreams`, `DecodeSingleName`, `LeftPadNumber`),
//! transliterated line by line, for testing [`super::tok3`] against.
//!
//! It exists to be obviously right, not fast, and shares no tok3 code with
//! the production decoder: the byte streams `B` are a map keyed by position
//! and type, the names `N` and tokens `T` plain vectors, every token a value.
//! Only the entropy decoders are shared — [`super::rans_nx16::decode`] and
//! [`super::arith::reference::decode`] (or whatever arithmetic decoder the
//! caller passes to [`decode_with`]), which are verified on their own.
//!
//! Where the pseudocode is silent, cannot run, or disagrees with htscodecs'
//! `tokenise_name3.c` (which writes the files) this follows htscodecs, each
//! place marked `htscodecs:`. Where htscodecs' behaviour on invalid input is
//! undefined or produces garbage, the input is rejected instead, each place
//! marked `Rejected:`. On everything htscodecs' encoder writes, all three
//! decode the same names. The differences from htscodecs' decoder, all on
//! blocks no encoder writes:
//! - the name count is the header's (as in `DecodeNames`); htscodecs decodes
//!   names until the position-0 type stream runs out, and rejects a count of 0;
//! - an unterminated STRING, a DELTA/DELTA0 whose previous token is not
//!   numeric (DIGITS / DIGITS0), and a token type other than CHAR, STRING,
//!   DIGITS, DIGITS0, DELTA, DELTA0, MATCH, NOP or END inside a name are
//!   rejected; htscodecs drops the string's last byte, adds the delta to
//!   whatever integer the token holds, and ends the name, respectively;
//! - a DIGITS0 whose length is above 9 or below its digit count is padded
//!   as `LeftPadNumber` says; htscodecs writes 9 digits at most and truncates
//!   the value's high digits;
//! - a stream header's type is `ttype & 63` and must name a token type
//!   (0..=12); htscodecs uses `ttype & 15`. A duplicated stream may come from
//!   any stream set before it; htscodecs also requires its (position, type)
//!   to sort before the target's and adds the type byte to the position
//!   unmasked;
//! - each stream is decoded from its `clen` bytes; htscodecs hands the
//!   entropy decoder the rest of the block.
//!
//! Beyond the format it applies the crate's resource policy, so it fails
//! exactly where the production decoder does: at most 10 million names
//! (htscodecs' limit), the codec output cap on the header's uncompressed
//! length, and no more output than that length plus 1 KiB (htscodecs'
//! margin: its encoder writes the exact length, noodles' one byte less).
//!
//! Every failure is `None`: the reference only has to say *whether* a block
//! decodes, and to what.

#![allow(
    clippy::arithmetic_side_effects,
    reason = "a transliteration of the pseudocode: counters bounded by the input and by the \
              checked output length, and `wrapping_add` where htscodecs' u32 arithmetic wraps"
)]

use std::collections::HashMap;

use super::reader::check_codec_output;

// Token types, `CRAMcodecs` §5's table.
const TYPE: u8 = 0;
const STRING: u8 = 1;
const CHAR: u8 = 2;
const DIGITS0: u8 = 3;
const DZLEN: u8 = 4;
const DUP: u8 = 5;
const DIGITS: u8 = 7;
const DELTA: u8 = 8;
const DELTA0: u8 = 9;
const MATCH: u8 = 10;
const NOP: u8 = 11;
const END: u8 = 12;

/// "A maximum of 128 tokens are permitted within any single read name";
/// htscodecs rejects a block with more token positions than that.
const MAX_TOKENS: usize = 128;

/// Not in the spec: htscodecs' limit on the name count.
const MAX_NAMES: usize = 10_000_000;

/// Not in the spec: how far past the header's uncompressed length the
/// output may run (htscodecs' margin).
const OUTPUT_SLACK: usize = 1024;

/// Decode a tok3 block, arithmetic-coded streams with
/// [`super::arith::reference::decode`]. `None` for an invalid block.
pub fn decode(src: &[u8]) -> Option<Vec<u8>> {
    decode_with(src, |data| super::arith::reference::decode(data, 0))
}

/// Decode a tok3 block, arithmetic-coded streams with `arith_decode`.
pub fn decode_with(src: &[u8], arith_decode: impl Fn(&[u8]) -> Option<Vec<u8>>) -> Option<Vec<u8>> {
    decode_names(&mut Input { bytes: src, pos: 0 }, &arith_decode)
}

/// The block being read.
struct Input<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl<'a> Input<'a> {
    fn read_uint8(&mut self) -> Option<u8> {
        let b = *self.bytes.get(self.pos)?;
        self.pos += 1;
        Some(b)
    }

    fn read_uint32(&mut self) -> Option<u32> {
        let bytes = self.read_data(4)?;
        Some(u32::from_le_bytes(bytes.try_into().ok()?))
    }

    /// `ReadUint7`: 7 bits a byte, most significant group first, the top
    /// bit set on all but the last byte; at most 5 bytes, keeping the low
    /// 32 bits (as the production decoder and the other codecs read it).
    fn read_uint7(&mut self) -> Option<usize> {
        let mut value: u64 = 0;
        for _ in 0..5 {
            let c = self.read_uint8()?;
            value = (value << 7) | u64::from(c & 127);
            if c < 128 {
                return usize::try_from(value & u64::from(u32::MAX)).ok();
            }
        }
        None
    }

    /// `ReadData(n)`: the next `n` bytes.
    fn read_data(&mut self, n: usize) -> Option<&'a [u8]> {
        let end = self.pos.checked_add(n)?;
        let data = self.bytes.get(self.pos..end)?;
        self.pos = end;
        Some(data)
    }

    fn eof(&self) -> bool {
        self.pos >= self.bytes.len()
    }
}

/// One byte stream `B_{pos,type}`, read from the front.
#[derive(Clone)]
struct ByteStream {
    bytes: Vec<u8>,
    pos: usize,
}

impl ByteStream {
    fn new(bytes: Vec<u8>) -> Self {
        Self { bytes, pos: 0 }
    }

    fn read_uint8(&mut self) -> Option<u8> {
        let b = *self.bytes.get(self.pos)?;
        self.pos += 1;
        Some(b)
    }

    fn read_uint32(&mut self) -> Option<u32> {
        let b = [self.read_uint8()?, self.read_uint8()?, self.read_uint8()?, self.read_uint8()?];
        Some(u32::from_le_bytes(b))
    }

    /// `ReadString`: bytes up to a NUL, which is consumed but not returned.
    /// Rejected: a string that runs off the end of the stream.
    fn read_string(&mut self) -> Option<Vec<u8>> {
        let mut s = Vec::new();
        loop {
            match self.read_uint8()? {
                0 => return Some(s),
                c => s.push(c),
            }
        }
    }
}

/// `B`, keyed by (position, type).
type Streams = HashMap<(usize, u8), ByteStream>;

/// Read from `B_{pos,type}`; a stream that was never set cannot be read.
fn stream(b: &mut Streams, pos: usize, ty: u8) -> Option<&mut ByteStream> {
    b.get_mut(&(pos, ty))
}

/// A token's value `T_{n,t}`, kept as a value rather than a string (the
/// spec's footnote allows either).
#[derive(Clone, Debug)]
enum Token {
    Char(u8),
    String(Vec<u8>),
    Digits(u32),
    /// A zero-padded number and the length of its string, which is what
    /// DELTA0 pads the next value to (`Length(T_{m,t})`).
    Digits0 {
        value: u32,
        len: usize,
    },
    /// NOP and END: the spec's empty string. MATCH may not refer to them.
    Empty,
}

impl Token {
    fn text(&self) -> Vec<u8> {
        match self {
            Self::Char(c) => vec![*c],
            Self::String(s) => s.clone(),
            Self::Digits(d) => d.to_string().into_bytes(),
            Self::Digits0 { value, len } => left_pad_number(*value, *len),
            Self::Empty => Vec::new(),
        }
    }

    /// A DIGITS0 value `d` padded to `l`, remembering the string's length.
    fn digits0(d: u32, l: usize) -> Self {
        Self::Digits0 { value: d, len: left_pad_number(d, l).len() }
    }
}

/// `LeftPadNumber`: `val` in base-10 digits, at least `len` bytes long with
/// leading zeros.
fn left_pad_number(val: u32, len: usize) -> Vec<u8> {
    let mut s = val.to_string().into_bytes();
    while s.len() < len {
        s.insert(0, b'0');
    }
    s
}

/// `DecodeNames`.
fn decode_names(
    input: &mut Input<'_>,
    arith_decode: &impl Fn(&[u8]) -> Option<Vec<u8>>,
) -> Option<Vec<u8>> {
    let ulen = usize::try_from(input.read_uint32()?).ok()?;
    let nnames = usize::try_from(input.read_uint32()?).ok()?;
    let use_arith = input.read_uint8()? != 0;

    // Not in the spec: the resource policy (see the module docs). Every name
    // adds at least its NUL, so more names than output bytes cannot decode.
    let max_output = ulen + OUTPUT_SLACK;
    if nnames > MAX_NAMES
        || check_codec_output(ulen, "tok3 reference").is_err()
        || nnames > max_output
    {
        return None;
    }

    let mut b = decode_token_byte_streams(input, use_arith, nnames, arith_decode)?;

    let mut names = Names { n: Vec::new(), t: Vec::new(), out: Vec::new(), max_output };
    for n in 0..nnames {
        names.decode_single_name(&mut b, n)?;
    }
    Some(names.out)
}

/// `DecodeTokenByteStreams`.
fn decode_token_byte_streams(
    input: &mut Input<'_>,
    use_arith: bool,
    nnames: usize,
    arith_decode: &impl Fn(&[u8]) -> Option<Vec<u8>>,
) -> Option<Streams> {
    let mut b = Streams::new();
    // `t ← -1`; `None` until the first new position.
    let mut t: Option<usize> = None;
    // htscodecs: `while (o < sz)`, where the pseudocode's repeat-until reads
    // one stream header even from an empty body.
    while !input.eof() {
        let ttype = input.read_uint8()?;
        let tok_new = ttype & 128 != 0;
        let tok_dup = ttype & 64 != 0;
        let ty = ttype & 63;
        // Rejected: a type that is not a token type (htscodecs masks with 15
        // and keeps streams for 13..=15 it never reads).
        if ty > END {
            return None;
        }
        if tok_new {
            let next = t.map_or(0, |t| t + 1);
            // htscodecs: at most MAX_TOKENS positions.
            if next >= MAX_TOKENS {
                return None;
            }
            t = Some(next);
            if ty != TYPE {
                // (type, MATCH, MATCH, ...), `nnames - 1` MATCHes.
                let mut types = vec![ty];
                types.resize(nnames.max(1), MATCH);
                b.insert((next, TYPE), ByteStream::new(types));
            }
        }
        // Rejected: a stream before the first position (`t = -1`).
        let t = t?;

        if tok_dup {
            let dup_pos = usize::from(input.read_uint8()?);
            let dup_type = input.read_uint8()?;
            // Rejected: copying a stream that was never set.
            let copy = b.get(&(dup_pos, dup_type))?.bytes.clone();
            b.insert((t, ty), ByteStream::new(copy));
        } else {
            let clen = input.read_uint7()?;
            let data = input.read_data(clen)?;
            // Every stream stores its uncompressed size, so none is passed.
            let decoded = if use_arith {
                arith_decode(data)?
            } else {
                super::rans_nx16::decode(data, 0).ok()?
            };
            b.insert((t, ty), ByteStream::new(decoded));
        }
    }
    Some(b)
}

/// The names decoded so far, `N`, their tokens, `T`, and the output.
struct Names {
    n: Vec<Vec<u8>>,
    /// `T[n][t - 1]` is `T_{n,t}`: positions start at 1.
    t: Vec<Vec<Token>>,
    out: Vec<u8>,
    max_output: usize,
}

impl Names {
    /// `T_{m,t}`, if name `m` has a token at `t`.
    fn token(&self, m: usize, t: usize) -> Option<&Token> {
        self.t.get(m)?.get(t - 1)
    }

    /// `DecodeSingleName`: decode name `n`, append it and its NUL to the
    /// output.
    fn decode_single_name(&mut self, b: &mut Streams, n: usize) -> Option<()> {
        let ty = stream(b, 0, TYPE)?.read_uint8()?;
        let dist = usize::try_from(stream(b, 0, ty)?.read_uint32()?).ok()?;
        // htscodecs: a distance past the first name is invalid.
        let m = n.checked_sub(dist)?;
        if ty == DUP {
            // htscodecs: a name cannot duplicate itself.
            if dist == 0 {
                return None;
            }
            let name = self.n.get(m)?.clone();
            let tokens = self.t.get(m)?.clone();
            return self.finish(name, tokens);
        }

        let mut name = Vec::new();
        let mut tokens: Vec<Token> = Vec::new();
        let mut t = 1;
        loop {
            let ty = stream(b, t, TYPE)?.read_uint8()?;
            let token = match ty {
                CHAR => Token::Char(stream(b, t, CHAR)?.read_uint8()?),
                STRING => Token::String(stream(b, t, STRING)?.read_string()?),
                DIGITS => Token::Digits(stream(b, t, DIGITS)?.read_uint32()?),
                DIGITS0 => {
                    let d = stream(b, t, DIGITS0)?.read_uint32()?;
                    let l = stream(b, t, DZLEN)?.read_uint8()?;
                    Token::digits0(d, usize::from(l))
                }
                DELTA => {
                    let delta = stream(b, t, DELTA)?.read_uint8()?;
                    // Rejected: a previous token that is not DIGITS (or a
                    // DELTA, which is stored as its DIGITS value).
                    // htscodecs: the sum wraps at 32 bits.
                    match self.token(m, t)? {
                        Token::Digits(v) => Token::Digits(v.wrapping_add(u32::from(delta))),
                        _ => return None,
                    }
                }
                DELTA0 => {
                    let delta = stream(b, t, DELTA0)?.read_uint8()?;
                    // Rejected: a previous token that is not DIGITS0 (or a
                    // DELTA0, stored as its DIGITS0 value).
                    // htscodecs: the sum wraps at 32 bits.
                    match self.token(m, t)? {
                        Token::Digits0 { value, len } => {
                            Token::digits0(value.wrapping_add(u32::from(delta)), *len)
                        }
                        _ => return None,
                    }
                }
                // htscodecs: MATCH may only refer to a CHAR, STRING, DIGITS
                // or DIGITS0 value (a matched MATCH or DELTA was stored as
                // one), not to NOP, END or a position the name does not
                // reach.
                MATCH => match self.token(m, t)? {
                    Token::Empty => return None,
                    token => token.clone(),
                },
                NOP | END => Token::Empty,
                // Rejected: TYPE, DZLEN, DUP, DIFF or an unknown type inside
                // a name. The pseudocode's `Else` makes them empty tokens;
                // htscodecs ends the name.
                _ => return None,
            };
            name.extend(token.text());
            tokens.push(token);
            // Not in the spec: stop as soon as the output would pass its
            // bound.
            if self.out.len() + name.len() + 1 > self.max_output {
                return None;
            }
            if ty == END {
                break;
            }
            t += 1;
        }
        self.finish(name, tokens)
    }

    /// Record `N_n` and `T_n`, and append the name and its NUL ("END: a nul
    /// byte is added to the name output buffer") to the output.
    fn finish(&mut self, name: Vec<u8>, tokens: Vec<Token>) -> Option<()> {
        if self.out.len() + name.len() + 1 > self.max_output {
            return None;
        }
        self.out.extend(&name);
        self.out.push(0);
        self.n.push(name);
        self.t.push(tokens);
        Some(())
    }
}

#[cfg(test)]
#[path = "../../tests/support/cram_raw_blocks.rs"]
#[allow(dead_code, reason = "the tests need only the tok3 blocks")]
mod cram_raw_blocks;

#[cfg(test)]
#[allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::cast_possible_truncation,
    reason = "test code; indices and sizes bounded by the generators"
)]
mod tests {
    use std::ffi::c_int;
    use std::io::Read;
    use std::path::PathBuf;

    use hegel::TestCase;
    use hegel::generators as gs;

    use super::cram_raw_blocks::{TOK3, raw_blocks};
    use super::left_pad_number;
    use crate::cram::arith::htscodecs_oracle::{htscodecs_tok3_decode, try_htscodecs_tokenise};
    use crate::cram::{arith, tok3};

    #[test]
    fn left_pad_number_pads_and_never_truncates() {
        assert_eq!(left_pad_number(123, 5), b"00123");
        assert_eq!(left_pad_number(123, 2), b"123");
        assert_eq!(left_pad_number(0, 0), b"0");
        assert_eq!(left_pad_number(7, 12), b"000000000007");
    }

    fn codecs_dir() -> PathBuf {
        PathBuf::from(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../tests/data/hts-specs/cram/codecs"
        ))
    }

    /// The gzipped original with each newline turned into tok3's NUL.
    fn original(stem: &str) -> Vec<u8> {
        let gz = std::fs::read(codecs_dir().join("gzip").join(format!("{stem}.gz"))).unwrap();
        let mut text = Vec::new();
        flate2::read::MultiGzDecoder::new(gz.as_slice()).read_to_end(&mut text).unwrap();
        text.iter().map(|&b| if b == b'\n' { 0 } else { b }).collect()
    }

    // r[verify cram.codec.tok3]
    // r[verify cram.codec.tok3_arith]
    #[test]
    fn reference_decodes_every_hts_specs_vector() {
        let mut n = 0;
        for entry in std::fs::read_dir(codecs_dir().join("tok3")).unwrap() {
            let path = entry.unwrap().path();
            let file_name = path.file_name().unwrap().to_str().unwrap().to_owned();
            let (stem, _) = file_name.rsplit_once('.').unwrap();
            let got = super::decode(&std::fs::read(&path).unwrap());
            assert!(got.as_deref() == Some(original(stem).as_slice()), "{file_name}");
            n += 1;
        }
        assert!(n >= 100, "expected the hts-specs tok3 vectors, found {n}");
    }

    /// The name blocks of the hts-specs CRAM 3.1 conformance files, tok3
    /// over rANS Nx16 (levels 2 and 3) and over the arithmetic coder
    /// (level 4), 20000 names each: the reference decodes them as htscodecs
    /// does, and to their declared size.
    // r[verify cram.codec.tok3]
    // r[verify cram.codec.tok3_arith]
    #[test]
    fn reference_decodes_the_conformance_files_names_as_htscodecs_does() {
        let mut n = 0;
        for level in 2..=4 {
            let path = codecs_dir().join(format!("../3.1/passed/level-{level}.cram"));
            for block in
                raw_blocks(&std::fs::read(path).unwrap()).iter().filter(|b| b.method == TOK3)
            {
                let ours = super::decode(&block.data).unwrap_or_else(|| panic!("level {level}"));
                assert_eq!(ours.len(), block.uncompressed_len, "level {level}");
                assert!(
                    Some(&ours) == htscodecs_tok3_decode(&block.data).as_ref(),
                    "level {level}"
                );
                n += 1;
            }
        }
        assert_eq!(n, 4, "the level-2 file has two name blocks, levels 3 and 4 one each");
    }

    /// The expected tok3 output: every name followed by a NUL.
    fn joined(names: &[Vec<u8>]) -> Vec<u8> {
        names.iter().flat_map(|n| n.iter().copied().chain([0])).collect()
    }

    /// The production decoder and the reference, over the same arithmetic
    /// decoder, agree on whether `block` decodes and to what — both with the
    /// production arithmetic decoder and with the reference one.
    // r[verify cram.codec.tok3.reference]
    fn assert_agree(block: &[u8]) {
        let production = tok3::decode(block).ok();
        let reference = super::decode_with(block, |s| arith::decode(s, 0).ok());
        assert_same("production arith", production.as_deref(), reference.as_deref());
        let production = tok3::decode_with_arith_reference(block).ok();
        let reference = super::decode(block);
        assert_same("reference arith", production.as_deref(), reference.as_deref());
    }

    fn assert_same(what: &str, production: Option<&[u8]>, reference: Option<&[u8]>) {
        let verdict = |o: Option<&[u8]>| {
            o.map_or_else(|| "fails".to_owned(), |o| format!("decodes {} bytes", o.len()))
        };
        assert!(
            production == reference,
            "{what}: production {} but reference {}",
            verdict(production),
            verdict(reference)
        );
    }

    /// A counter as a DIGITS or DIGITS0 field, zero-padded to `width`.
    fn padded(n: u64, width: usize) -> String {
        format!("{n:0width$}")
    }

    /// Sets of read names in the shapes tok3 meets: Illumina names sharing
    /// a prefix with numeric fields that creep forward, counters (padded
    /// or not) that step by small and large deltas, random printable
    /// strings, and a mix; with pairs, duplicates of earlier names, empty
    /// names, and sets of zero or one name.
    #[hegel::composite]
    fn arb_names(tc: &TestCase) -> Vec<Vec<u8>> {
        let count = tc.draw_silent(gs::integers::<usize>().max_value(400));
        let style = tc.draw_silent(gs::integers::<u8>().max_value(4));
        let paired = tc.draw_silent(gs::booleans());
        let width = tc.draw_silent(gs::integers::<usize>().min_value(1).max_value(12));
        let prefix = tc.draw_silent(gs::sampled_from(
            &["HWI-ST1234", "A00123:8:H7VJWDSXX", "SRR062634.", "m54006_171109", "r", ""][..],
        ));
        let mut counter = tc.draw_silent(gs::integers::<u64>().max_value(1 << 34));
        let (mut lane, mut tile, mut y) = (1u32, 1101u32, 2000u32);
        let mut names: Vec<Vec<u8>> = Vec::with_capacity(count);
        while names.len() < count {
            // A duplicate of a recent name.
            if !names.is_empty() && tc.draw_silent(gs::integers::<u8>().max_value(15)) == 0 {
                let back =
                    tc.draw_silent(gs::integers::<usize>().max_value(names.len().min(8) - 1));
                names.push(names[names.len() - 1 - back].clone());
                continue;
            }
            let step = tc.draw_silent(gs::sampled_from(
                &[0u64, 1, 1, 1, 2, 3, 17, 255, 256, 999, 65_536][..],
            ));
            counter = counter.saturating_add(step);
            let style =
                if style == 4 { tc.draw_silent(gs::integers::<u8>().max_value(3)) } else { style };
            let name = match style {
                0 => {
                    if tc.draw_silent(gs::integers::<u8>().max_value(50)) == 0 {
                        lane = lane % 8 + 1;
                    }
                    if tc.draw_silent(gs::integers::<u8>().max_value(20)) == 0 {
                        tile += 1;
                        y = 0;
                    }
                    let x = tc.draw_silent(gs::integers::<u32>().max_value(30_000));
                    y = y.saturating_add(u32::try_from(step % 5000).unwrap());
                    format!("{prefix}:{lane}:{tile}:{x}:{}", padded(u64::from(y), 5))
                }
                1 => format!("{prefix}{}", padded(counter, width)),
                2 => format!("{prefix}{counter}_{}", counter % 7),
                _ => tc
                    .draw_silent(
                        gs::vecs(gs::integers::<u8>().min_value(b' ').max_value(b'~')).max_size(40),
                    )
                    .into_iter()
                    .map(char::from)
                    .collect(),
            };
            if paired {
                for mate in ["/1", "/2"] {
                    if names.len() < count {
                        names.push(format!("{name}{mate}").into_bytes());
                    }
                }
            } else {
                names.push(name.into_bytes());
            }
        }
        names
    }

    /// htscodecs' tokenisation of `names` at a drawn level and entropy
    /// coder; `None` where htscodecs gives up (a name of more than 128
    /// tokens, or no names).
    fn encode(tc: &TestCase, names: &[Vec<u8>]) -> Option<Vec<u8>> {
        let level = tc.draw(gs::integers::<c_int>().min_value(1).max_value(9));
        let use_arith = tc.draw(gs::booleans());
        try_htscodecs_tokenise(names, level, use_arith)
    }

    // r[verify cram.codec.tok3]
    // r[verify cram.codec.tok3_arith]
    #[hegel::test(test_cases = 300)]
    fn both_decoders_give_back_htscodecs_names(tc: TestCase) {
        let names = tc.draw(arb_names());
        let Some(block) = encode(&tc, &names) else {
            // htscodecs only gives up on no names or on huge token counts,
            // which the random strings can reach.
            assert!(names.is_empty() || names.iter().any(|n| n.len() > 60));
            return;
        };
        let expected = joined(&names);
        assert!(super::decode(&block).as_deref() == Some(expected.as_slice()), "reference");
        assert!(tok3::decode(&block).unwrap() == expected, "production");
    }

    /// Blocks whose names end near the bound on the output (the declared
    /// length plus 1 KiB), where the decoders' writes take their slow paths:
    /// they decode exactly when the names fit, and alike.
    // r[verify cram.codec.tok3.limits]
    #[hegel::test(test_cases = 300)]
    fn decoders_agree_near_the_output_bound(tc: TestCase) {
        let names = tc.draw(arb_names());
        let Some(mut block) = encode(&tc, &names) else { return };
        let len = joined(&names).len();
        let short = tc.draw(gs::integers::<usize>().max_value(40));
        let ulen = (len + short).saturating_sub(1024 + 20);
        block[..4].copy_from_slice(&u32::try_from(ulen).unwrap().to_le_bytes());
        assert_agree(&block);
        assert_eq!(tok3::decode(&block).is_ok(), len <= ulen + 1024 && names.len() <= ulen + 1024);
    }

    /// Damaged blocks: whatever both decoders accept, they decode alike.
    #[hegel::test(test_cases = 500)]
    fn decoders_agree_on_damaged_blocks(tc: TestCase) {
        let names = tc.draw(arb_names());
        let Some(mut block) = encode(&tc, &names) else { return };
        match tc.draw(gs::integers::<u8>().max_value(2)) {
            0 => {
                let flips: Vec<(usize, u8)> = tc.draw(
                    gs::vecs(gs::tuples!(
                        gs::integers::<usize>(),
                        gs::integers::<u8>().min_value(1)
                    ))
                    .min_size(1)
                    .max_size(4),
                );
                for (at, x) in flips {
                    let len = block.len();
                    block[at % len] ^= x;
                }
            }
            1 => {
                let cut = tc.draw(gs::integers::<usize>().max_value(block.len()));
                block.truncate(cut);
            }
            _ => {
                let cut = tc.draw(gs::integers::<usize>().max_value(block.len()));
                block.truncate(cut);
                block.extend(tc.draw(gs::binary().max_size(256)));
            }
        }
        assert_agree(&block);
    }

    fn put_uint7(out: &mut Vec<u8>, v: usize) {
        let v = u32::try_from(v).unwrap();
        let groups = (0..5u32).rev().skip_while(|&g| g > 0 && v >> (7 * g) == 0);
        for g in groups {
            let bits = ((v >> (7 * g)) & 0x7f) as u8;
            out.push(if g == 0 { bits } else { bits | 0x80 });
        }
    }

    /// One raw (CAT) token stream: `ttype`, then the stream as the
    /// arithmetic coder and rANS Nx16 both store uncompressed data.
    fn put_stream(block: &mut Vec<u8>, ttype: u8, data: &[u8]) {
        let mut stream = vec![0x20];
        put_uint7(&mut stream, data.len());
        stream.extend(data);
        block.push(ttype);
        put_uint7(block, stream.len());
        block.extend(stream);
    }

    /// `n` draws of `values`.
    fn draw_n<T: Copy + Send + Sync + std::fmt::Debug + 'static>(
        tc: &TestCase,
        n: usize,
        values: &'static [T],
    ) -> Vec<T> {
        tc.draw_silent(gs::vecs(gs::sampled_from(values)).min_size(n).max_size(n))
    }

    /// Token streams written directly rather than by an encoder, each stored
    /// raw, so they reach what no encoder writes: DELTA past 2^32, over-long
    /// and too-short DIGITS0, MATCH of NOP, END or a missing token, DUP and
    /// DIFF distances out of range, strings without their NUL, invalid
    /// token types, elided and duplicated type streams, and names past the
    /// header's length. Mostly well-formed enough that names decode.
    #[hegel::composite]
    fn arb_raw_block(tc: &TestCase) -> Vec<u8> {
        let nnames = tc.draw_silent(gs::integers::<usize>().max_value(8));
        // Most blocks are well-formed but for their values; the rest break
        // a rule somewhere, which one drawn at each place below.
        let noisy = tc.draw_silent(gs::integers::<u8>().max_value(3)) == 0;
        let odd = |max: u8| noisy && tc.draw_silent(gs::integers::<u8>().max_value(max)) == 0;
        let ulen = if odd(3) { tc.draw_silent(gs::integers::<u32>().max_value(40)) } else { 1000 };
        let mut block = Vec::new();
        block.extend(ulen.to_le_bytes());
        block.extend(u32::try_from(nnames).unwrap().to_le_bytes());
        block.push(u8::from(tc.draw_silent(gs::booleans())));

        // Position 0: DIFF or DUP, and each kind's distances in the order
        // its names read them.
        let mut kinds = Vec::new();
        let (mut diff, mut dup) = (Vec::new(), Vec::new());
        for n in 0..u32::try_from(nnames).unwrap() {
            let kind = if n > 0 || odd(3) {
                tc.draw_silent(gs::sampled_from(&[6u8, 6, 6, 5][..]))
            } else {
                6
            };
            let d = if odd(3) {
                tc.draw_silent(gs::sampled_from(&[0, n + 1][..]))
            } else if kind == 5 || tc.draw_silent(gs::booleans()) {
                tc.draw_silent(gs::integers::<u32>().max_value(n)).max(n.min(1))
            } else {
                n.min(1)
            };
            kinds.push(kind);
            if kind == 5 { &mut dup } else { &mut diff }.extend(d.to_le_bytes());
        }
        put_stream(&mut block, 0x80, &kinds);
        put_stream(&mut block, 6, &diff);
        if !dup.is_empty() || odd(1) {
            put_stream(&mut block, 5, &dup);
        }

        let positions = tc.draw_silent(gs::integers::<u8>().min_value(1).max_value(5));
        for pos in 1..=positions {
            // The type stream: each position a column of one kind, with
            // MATCH (and DELTA for numbers) after the first name, and END
            // at the last position.
            let column: &'static [u8] = if pos == positions && !odd(4) {
                &[12]
            } else {
                tc.draw_silent(gs::sampled_from(
                    &[
                        &[1u8, 10, 10][..],
                        &[2, 10, 10],
                        &[7, 8, 8, 10],
                        &[3, 9, 9, 10],
                        &[11],
                        &[1, 2, 7, 11, 12],
                    ][..],
                ))
            };
            let types: Vec<u8> = (0..nnames)
                .map(|n| {
                    if odd(20) {
                        tc.draw_silent(gs::sampled_from(&[0u8, 4, 5, 6, 8, 9, 10, 13, 0x8a][..]))
                    } else if n == 0 || column.len() == 1 {
                        column[0]
                    } else {
                        tc.draw_silent(gs::sampled_from(column))
                    }
                })
                .collect();
            let first = types.first().copied().unwrap_or(12);
            let mut new = 0x80;
            match tc.draw_silent(gs::integers::<u8>().max_value(5)) {
                // Elided: all MATCH after the first, regenerated from the
                // first value stream's type.
                0 if first <= 12 && first != 0 => {}
                // A copy of an earlier position's type stream.
                1 if noisy => {
                    block.extend([0xC0, tc.draw_silent(gs::integers::<u8>().max_value(pos)), 0]);
                    new = 0;
                }
                _ => {
                    put_stream(&mut block, 0x80, &types);
                    new = 0;
                }
            }
            let mut value_types = vec![first];
            value_types.extend(types.iter().copied().filter(|&t| t != first));
            value_types.dedup();
            if value_types.contains(&3) {
                value_types.push(4);
            }
            for ty in value_types {
                if ty == 0 || ty > 12 {
                    continue;
                }
                let data: Vec<u8> = match ty {
                    1 => {
                        let mut s: Vec<u8> = (0..nnames)
                            .flat_map(|_| {
                                let mut s = draw_n(
                                    tc,
                                    tc.draw_silent(gs::integers::<usize>().max_value(3)),
                                    b"aB:_",
                                );
                                s.push(0);
                                s
                            })
                            .collect();
                        if odd(5) {
                            s.pop();
                        }
                        s
                    }
                    2 => draw_n(tc, nnames, &[b'a', b':', b'/', 0]),
                    3 | 7 => {
                        draw_n(tc, nnames, &[0u32, 1, 9, 10, 99, 123_456, 4_294_967_040, u32::MAX])
                            .into_iter()
                            .flat_map(u32::to_le_bytes)
                            .collect()
                    }
                    4 => draw_n(tc, nnames, &[0u8, 1, 2, 3, 5, 9, 10, 12, 20]),
                    8 | 9 => draw_n(tc, nnames, &[0u8, 1, 2, 200, 255]),
                    _ => continue,
                };
                // Now and then a stream copied from another instead.
                if tc.draw_silent(gs::integers::<u8>().max_value(15)) == 0 {
                    block.extend([
                        new | 0x40 | ty,
                        tc.draw_silent(gs::integers::<u8>().max_value(pos)),
                        tc.draw_silent(gs::integers::<u8>().max_value(12)),
                    ]);
                } else {
                    put_stream(&mut block, new | ty, &data);
                }
                new = 0;
            }
            if new != 0 {
                // No stream opened this position: an explicit type stream.
                put_stream(&mut block, 0x80, &types);
            }
        }
        // Now and then a damaged byte.
        if odd(3) {
            let at = tc.draw_silent(gs::integers::<usize>().max_value(block.len() - 1));
            block[at] ^= tc.draw_silent(gs::integers::<u8>().min_value(1));
        }
        block
    }

    #[hegel::test(test_cases = 2000)]
    fn decoders_agree_on_raw_blocks(tc: TestCase) {
        assert_agree(&tc.draw(arb_raw_block()));
    }

    #[hegel::test(test_cases = 500)]
    fn decoders_agree_on_arbitrary_bytes(tc: TestCase) {
        // Half the time a plausible header: a small length and name count.
        let mut block = Vec::new();
        if tc.draw(gs::booleans()) {
            block.extend(tc.draw(gs::integers::<u32>().max_value(1000)).to_le_bytes());
            block.extend(tc.draw(gs::integers::<u32>().max_value(20)).to_le_bytes());
            block.push(tc.draw(gs::integers::<u8>().max_value(1)));
        }
        block.extend(tc.draw(gs::binary().max_size(512)));
        assert_agree(&block);
    }
}
