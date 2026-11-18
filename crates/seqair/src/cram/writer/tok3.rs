//! tok3 name tokeniser encoder (CRAM 3.1 method 8): the inverse of
//! [`crate::cram::tok3`], tokenising as htscodecs' `tokenise_name3.c`
//! encoder does so its token streams are the ones htscodecs compresses.
//!
//! Each name is a `DUP` of an identical earlier name, or a `DIFF` against
//! one earlier name whose tokens it is compared with, position by position:
//! a token equal to that name's is a `MATCH`, a number a little above it a
//! `DELTA`. Every (position, token type) pair collects its bytes in a stream
//! of its own, and each stream is entropy-coded by [`compress_stream`].

// r[impl cram_writer.tok3]

use rustc_hash::FxHashMap;

// Token types, as the decoder numbers them.
const TYPE: u8 = 0;
const ALPHA: u8 = 1;
const CHAR: u8 = 2;
const DIGITS0: u8 = 3;
const DZLEN: u8 = 4;
const DUP: u8 = 5;
const DIFF: u8 = 6;
const DIGITS: u8 = 7;
const DELTA: u8 = 8;
const DELTA0: u8 = 9;
const MATCH: u8 = 10;
const END: u8 = 12;
const N_TYPES: usize = 13;

/// Token positions a block may have, the first holding `DUP`/`DIFF`
/// (htscodecs' `MAX_TOKENS`, which both decoders enforce).
const MAX_POSITIONS: usize = 128;
/// The last position a name's text may use: `END` takes the one after.
const LAST_TEXT_POSITION: usize = MAX_POSITIONS - 2;

/// A digit run is cut after this many digits, so its value fits a `u32`.
const MAX_DIGITS: usize = 9;

/// Stream header flags.
const NEW_POSITION: u8 = 0x80;
const COPY: u8 = 0x40;

/// The block header's `use_arith` byte: every stream below is rANS Nx16.
/// An arithmetic-coder variant would set it and code every stream with the
/// adaptive arithmetic coder instead, as one block uses one coder.
const USE_ARITH: u8 = 0;

/// The names of an `RN` block, each followed by a NUL, as a tok3 block; or
/// `None` for input the format cannot hold: empty (htscodecs' decoder
/// refuses a block of no names), without a NUL after its last name, or
/// past `u32` sizes.
pub(super) fn encode(names: &[u8]) -> Option<Vec<u8>> {
    let len = u32::try_from(names.len()).ok()?;
    let body = names.strip_suffix(&[0])?;
    let mut encoder = Encoder::default();
    for name in body.split(|&b| b == 0) {
        encoder.name(name)?;
    }
    let count = u32::try_from(encoder.names.len()).ok()?;

    let mut out = Vec::new();
    out.extend_from_slice(&len.to_le_bytes());
    out.extend_from_slice(&count.to_le_bytes());
    out.push(USE_ARITH);
    write_streams(&encoder.positions, &mut out)?;
    Some(out)
}

/// A token stream as an rANS Nx16 stream: the smallest of the flag bytes
/// htscodecs' `compress()` tries for its token type at the level htslib
/// asks for (3; `R[1]` in `tokenise_name3.c`). A STRIPE combination needs
/// a length divisible by four.
fn compress_stream(token_type: u8, data: &[u8]) -> Option<Vec<u8>> {
    let methods: &[u8] = match token_type {
        TYPE => &[192, 0],
        ALPHA => &[129, 1],
        DIGITS0 => &[136, 0],
        DUP | DIGITS => &[200],
        DIFF => &[136],
        DELTA0 => &[128],
        _ => &[0],
    };
    methods
        .iter()
        .filter(|&&m| m & 8 == 0 || data.len().is_multiple_of(4))
        .filter_map(|&m| super::rans_nx16::encode(data, super::rans_nx16::Nx16Flags::from_byte(m)))
        .min_by_key(Vec::len)
        .or_else(|| super::rans_nx16::encode_cat(data))
}

/// `v` as a uint7: 7 bits a byte, most significant first, the high bit
/// set on every byte but the last.
#[allow(clippy::cast_possible_truncation, reason = "masked to 7 bits")]
fn put_uint7(out: &mut Vec<u8>, v: u32) {
    for shift in [28, 21, 14, 7] {
        if v >> shift != 0 {
            out.push(((v >> shift) & 0x7f) as u8 | 0x80);
        }
    }
    out.push((v & 0x7f) as u8);
}

/// A token of an encoded name, kept for the names diffed against it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Token<'a> {
    Alpha(&'a [u8]),
    Char(u8),
    Digits(u32),
    /// A digit run stored with its length, leading zeros and all.
    Digits0 {
        value: u32,
        len: u8,
    },
}

/// One token position's streams, indexed by token type, and how its
/// numbers were stored so far.
#[derive(Default)]
struct Position {
    streams: [Vec<u8>; N_TYPES],
    deltas: u32,
    digits: u32,
}

impl Position {
    /// Append a token of type `ty`: its type to the type stream, `value` to
    /// its own stream.
    fn put(&mut self, ty: u8, value: &[u8]) -> Option<()> {
        self.streams.first_mut()?.push(ty);
        self.streams.get_mut(usize::from(ty))?.extend_from_slice(value);
        Some(())
    }

    /// Append `value` to the stream of type `ty` alone (`DZLEN`).
    fn put_value(&mut self, ty: u8, value: &[u8]) -> Option<()> {
        self.streams.get_mut(usize::from(ty))?.extend_from_slice(value);
        Some(())
    }
}

#[derive(Default)]
struct Encoder<'a> {
    positions: Vec<Position>,
    /// Every name's tokens, without its `END`.
    tokens: Vec<Token<'a>>,
    /// Per name, where its tokens are in `tokens`; a `DUP` shares its
    /// original's.
    names: Vec<(usize, usize)>,
    /// The most recent name of each text, for `DUP`.
    seen: FxHashMap<&'a [u8], usize>,
    /// The most recent name of each [`known_prefix`], to diff against.
    prefixes: FxHashMap<&'a [u8], usize>,
    /// The current name's tokens.
    scratch: Vec<Token<'a>>,
}

/// The position `t`, made if the block has none yet.
fn position(positions: &mut Vec<Position>, t: usize) -> Option<&mut Position> {
    while positions.len() <= t {
        positions.push(Position::default());
    }
    positions.get_mut(t)
}

impl<'a> Encoder<'a> {
    fn name(&mut self, name: &'a [u8]) -> Option<()> {
        let n = self.names.len();
        let prefix = known_prefix(name);
        let prefix_previous =
            prefix.and_then(|p| name.get(..p.len)).and_then(|p| self.prefixes.insert(p, n));

        // As htscodecs, an empty name is never a duplicate.
        if let Some(m) = self.seen.insert(name, n)
            && !name.is_empty()
        {
            let distance = u32::try_from(n.wrapping_sub(m)).ok()?;
            position(&mut self.positions, 0)?.put(DUP, &distance.to_le_bytes())?;
            let tokens = *self.names.get(m)?;
            self.names.push(tokens);
            return Some(());
        }

        // htscodecs diffs a name with a prefix against the last one sharing
        // it, or against nothing; any other against the name before.
        let previous = if prefix.is_some() { prefix_previous } else { n.checked_sub(1) };
        let distance = previous.map_or(0, |m| n.wrapping_sub(m));
        let distance = u32::try_from(distance).ok()?;
        position(&mut self.positions, 0)?.put(DIFF, &distance.to_le_bytes())?;
        let (first, count) = previous.and_then(|m| self.names.get(m)).copied().unwrap_or((0, 0));
        let previous = self.tokens.get(first..first.wrapping_add(count)).unwrap_or_default();

        self.scratch.clear();
        let mut rest = name;
        if let Some((head, tail)) =
            prefix.filter(|p| p.token).and_then(|p| name.split_at_checked(p.len))
        {
            self.scratch.push(Token::Alpha(head));
            rest = tail;
        }
        while !rest.is_empty() {
            let t = self.scratch.len().wrapping_add(1);
            let (token, tail) = if t == LAST_TEXT_POSITION {
                // htscodecs gives up on a name this long; its rest is one token.
                (Token::Alpha(rest), &[][..])
            } else {
                next_token(rest, previous.get(t.wrapping_sub(1)).copied())?
            };
            self.scratch.push(token);
            rest = tail;
        }

        for (i, &token) in self.scratch.iter().enumerate() {
            let position = position(&mut self.positions, i.wrapping_add(1))?;
            put_token(position, token, previous.get(i).copied())?;
        }
        position(&mut self.positions, self.scratch.len().wrapping_add(1))?.put(END, &[])?;

        let first = self.tokens.len();
        self.tokens.extend_from_slice(&self.scratch);
        self.names.push((first, self.scratch.len()));
        Some(())
    }
}

/// A name's constant prefix, as htscodecs guesses it from the platform's
/// naming scheme.
#[derive(Clone, Copy)]
struct Prefix {
    len: usize,
    /// Whether the prefix is the name's first token.
    token: bool,
}

/// htscodecs' `search_trie` platform checks, with its indices: a `PacBio`
/// movie name, an Ion Torrent name, an ONT read UUID, or an Illumina name,
/// whose prefix ends at the fourth colon from the end of its text before
/// the first space or control byte.
#[allow(clippy::arithmetic_side_effects, reason = "f is 0 or 1, plus small constants")]
fn known_prefix(name: &[u8]) -> Option<Prefix> {
    let d = name.strip_prefix(b"@").unwrap_or(name);
    let f = usize::from(name.first() == Some(&b'>'));
    let at = |i: usize| d.get(i).copied().unwrap_or(0);
    let hex = |b: u8| b.is_ascii_digit() || (b'a'..=b'f').contains(&b);
    if d.len() > 70 && at(f) == b'm' && at(7) == b'_' && at(f + 14) == b'_' && at(f + 61) == b'/' {
        return Some(Prefix { len: 60, token: false });
    }
    if d.len() == 17 && at(f + 5) == b':' && at(f + 11) == b':' {
        return Some(Prefix { len: 6, token: true });
    }
    if d.len() > 37
        && [8, 13, 18, 23].iter().all(|&i| at(f + i) == b'-')
        && hex(at(f))
        && hex(at(f + 35))
    {
        return Some(Prefix { len: 37, token: true });
    }
    let end = name.iter().position(|&b| b <= b' ').unwrap_or(name.len());
    let (colon, _) =
        name.get(..end)?.iter().enumerate().rev().filter(|&(_, &b)| b == b':').nth(3)?;
    Some(Prefix { len: colon + 1, token: true })
}

/// The token `rest` starts with, and what follows it; `previous` is the
/// token at this position of the name diffed against.
fn next_token<'a>(rest: &'a [u8], previous: Option<Token<'_>>) -> Option<(Token<'a>, &'a [u8])> {
    let &c = rest.first()?;
    if c.is_ascii_alphabetic() {
        let len = rest
            .iter()
            .position(|b| !(b.is_ascii_alphabetic() || b.is_ascii_punctuation()))
            .unwrap_or(rest.len());
        if len > 1 {
            let (head, tail) = rest.split_at_checked(len)?;
            return Some((Token::Alpha(head), tail));
        }
    } else if c.is_ascii_digit() {
        let len = rest.iter().take(MAX_DIGITS).take_while(|b| b.is_ascii_digit()).count();
        let (head, tail) = rest.split_at_checked(len)?;
        // At most nine digits: no overflow.
        let value = head
            .iter()
            .fold(0u32, |v, &d| v.wrapping_mul(10).wrapping_add(u32::from(d.wrapping_sub(b'0'))));
        let len = u8::try_from(len).ok()?;
        // A run as long as a zero-padded one before it is likely a
        // fixed-width column that happens to have no leading zero here.
        let fixed_width = matches!(previous, Some(Token::Digits0 { len: l, .. }) if l == len);
        let token = if c == b'0' || fixed_width {
            Token::Digits0 { value, len }
        } else {
            Token::Digits(value)
        };
        return Some((token, tail));
    }
    Some((Token::Char(c), rest.get(1..)?))
}

/// Append `token` to its position's streams, as a `MATCH` or `DELTA` of
/// `previous` where it can be.
fn put_token(position: &mut Position, token: Token<'_>, previous: Option<Token<'_>>) -> Option<()> {
    if previous == Some(token) {
        return position.put(MATCH, &[]);
    }
    match token {
        Token::Alpha(text) => {
            position.put(ALPHA, text)?;
            position.put_value(ALPHA, &[0])
        }
        Token::Char(c) => position.put(CHAR, &[c]),
        Token::Digits0 { value, len } => {
            if let Some(Token::Digits0 { value: p, len: l }) = previous
                && l == len
                && let Some(delta) = small_delta(p, value)
            {
                return position.put(DELTA0, &[delta]);
            }
            position.put_value(DZLEN, &[len])?;
            position.put(DIGITS0, &value.to_le_bytes())
        }
        Token::Digits(value) => {
            if let Some(Token::Digits(p)) = previous
                && let Some(delta) = small_delta(p, value)
                // htscodecs' heuristic: deltas while they keep up with plain values.
                && position.deltas.saturating_add(5) > position.digits
            {
                position.deltas = position.deltas.saturating_add(1);
                return position.put(DELTA, &[delta]);
            }
            if matches!(previous, Some(Token::Digits(_))) {
                position.digits = position.digits.saturating_add(1);
            }
            position.put(DIGITS, &value.to_le_bytes())
        }
    }
}

/// `value - previous` when it is in `1..256`.
fn small_delta(previous: u32, value: u32) -> Option<u8> {
    u8::try_from(value.checked_sub(previous)?).ok()
}

/// Every position's non-empty streams, each through [`compress_stream`],
/// in position and type order; the first of a position is flagged as such.
fn write_streams(positions: &[Position], out: &mut Vec<u8>) -> Option<()> {
    // Streams written so far, by position and type, to copy identical ones.
    let mut written: Vec<(u8, u8, Vec<u8>)> = Vec::new();
    for (p, position) in positions.iter().enumerate() {
        let p = u8::try_from(p).ok()?;
        let skip_types = regenerable_types(position);
        let mut new_position = NEW_POSITION;
        for (ty, data) in position.streams.iter().enumerate() {
            if data.is_empty() || (ty == usize::from(TYPE) && skip_types) {
                continue;
            }
            let ty = u8::try_from(ty).ok()?;
            let stream = compress_stream(ty, data)?;
            let flags = ty | new_position;
            new_position = 0;
            // As htscodecs: a copy costs 3 bytes, so only longer streams.
            if stream.len() > 4
                && let Some((from_p, from_ty, _)) = written.iter().find(|(_, _, s)| *s == stream)
            {
                out.extend_from_slice(&[flags | COPY, *from_p, *from_ty]);
            } else {
                out.push(flags);
                put_uint7(out, u32::try_from(stream.len()).ok()?);
                out.extend_from_slice(&stream);
                written.push((p, ty, stream));
            }
        }
    }
    Some(())
}

/// Whether the decoder can regenerate the position's type stream, which it
/// does from the type of the position's first stream: the type stream is
/// that type, then only `MATCH`.
fn regenerable_types(position: &Position) -> bool {
    let Some((&first, rest)) = position.streams.first().and_then(|s| s.split_first()) else {
        return false;
    };
    let first_stream = position.streams.iter().skip(1).position(|s| !s.is_empty());
    rest.iter().all(|&t| t == MATCH) && first_stream == Some(usize::from(first).wrapping_sub(1))
}

#[cfg(test)]
#[path = "../../../tests/support/cram_raw_blocks.rs"]
#[allow(
    dead_code,
    clippy::duplicate_mod,
    reason = "the tests need only the tok3 blocks; the decoder's tests include it privately too"
)]
mod cram_raw_blocks;

#[cfg(test)]
#[allow(
    clippy::unwrap_used,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    reason = "test code; sizes bounded by the generators"
)]
mod tests {
    use std::collections::BTreeMap;
    use std::io::Read;

    use hegel::TestCase;
    use hegel::generators as gs;

    use super::cram_raw_blocks::{TOK3, raw_blocks};
    use super::*;
    use crate::cram::arith::htscodecs_oracle::{htscodecs_tok3_decode, try_htscodecs_tokenise};
    use crate::cram::codec_io::{read_u8, read_uint7, split_off};
    use crate::cram::{arith, rans_nx16, tok3, tok3_reference};

    /// Every name followed by a NUL, as the `RN` block holds them.
    fn joined(names: &[Vec<u8>]) -> Vec<u8> {
        names.iter().flat_map(|n| n.iter().copied().chain([0])).collect()
    }

    /// `names` encoded, after checking that seqair's decoder, the reference
    /// decoder and htscodecs' all give them back.
    fn round_trip(names: &[u8]) -> Vec<u8> {
        let block = encode(names).unwrap();
        assert!(tok3::decode(&block).unwrap() == names, "production decoder");
        assert!(tok3_reference::decode(&block).as_deref() == Some(names), "reference decoder");
        assert!(htscodecs_tok3_decode(&block).as_deref() == Some(names), "htscodecs");
        block
    }

    /// A block's token streams, decoded, by position and type; copies
    /// resolved to the bytes they copy.
    fn streams(block: &[u8]) -> BTreeMap<(u8, u8), Vec<u8>> {
        let mut cur = block.get(9..).unwrap();
        let mut out = BTreeMap::new();
        let mut position = None::<u8>;
        while let Some(flags) = read_u8(&mut cur) {
            if flags & NEW_POSITION != 0 {
                position = Some(position.map_or(0, |p| p + 1));
            }
            let key = (position.unwrap(), flags & 0x3f);
            let data = if flags & COPY != 0 {
                let from = (read_u8(&mut cur).unwrap(), read_u8(&mut cur).unwrap());
                out.get(&from).cloned().unwrap()
            } else {
                let len = read_uint7(&mut cur).unwrap() as usize;
                let data = split_off(&mut cur, len).unwrap();
                if block[8] == 0 {
                    rans_nx16::decode(data, 0).unwrap()
                } else {
                    arith::decode(data, 0).unwrap()
                }
            };
            out.insert(key, data);
        }
        out
    }

    /// Illumina-style names as a coordinate-sorted slice holds them: a few
    /// run prefixes, lanes and tiles, coordinates in zero-padded or plain
    /// columns, both mates of some templates, and duplicates. Every name
    /// ends in a mate suffix, so none is a strict prefix of another, where
    /// htscodecs' trie picks a name to diff against that this encoder does
    /// not look for.
    #[hegel::composite]
    fn illumina_names(tc: &TestCase) -> Vec<Vec<u8>> {
        let count = tc.draw_silent(gs::integers::<usize>().min_value(1).max_value(300));
        let padded = tc.draw_silent(gs::booleans());
        let mut names: Vec<Vec<u8>> = Vec::with_capacity(count);
        let mut y = tc.draw_silent(gs::integers::<u32>().max_value(99_999));
        while names.len() < count {
            if !names.is_empty() && tc.draw_silent(gs::integers::<u8>().max_value(7)) == 0 {
                let back = tc.draw_silent(gs::integers::<usize>().max_value(names.len() - 1));
                names.push(names[names.len() - 1 - back].clone());
                continue;
            }
            let run = tc.draw_silent(gs::sampled_from(
                &["A00123:8:H7VJWDSXX", "HSQ1008_141", "VP2-06:112:H7LNDMCVY", "ERR174310.1"][..],
            ));
            let lane = tc.draw_silent(gs::integers::<u8>().min_value(1).max_value(2));
            let tile = tc.draw_silent(gs::sampled_from(&[1101u32, 1102, 2101, 2214][..]));
            let x = tc.draw_silent(gs::integers::<u32>().max_value(40_000));
            let step = tc.draw_silent(gs::sampled_from(&[0u32, 1, 3, 17, 255, 256, 4000][..]));
            y = (y + step) % 100_000;
            let y = if padded { format!("{y:05}") } else { y.to_string() };
            let mate = tc.draw_silent(gs::integers::<u8>().min_value(1).max_value(2));
            names.push(format!("{run}:{lane}:{tile}:{x}:{y}/{mate}").into_bytes());
        }
        names
    }

    /// Names in any shape: counters with leading zeros, long digit runs,
    /// many short tokens (past the 126 a name may have), any bytes but NUL,
    /// empty names; with duplicates.
    #[hegel::composite]
    fn any_names(tc: &TestCase) -> Vec<Vec<u8>> {
        let count = tc.draw_silent(gs::integers::<usize>().min_value(1).max_value(100));
        let mut names: Vec<Vec<u8>> = Vec::with_capacity(count);
        let mut counter = tc.draw_silent(gs::integers::<u64>().max_value(1 << 40));
        while names.len() < count {
            let name = match tc.draw_silent(gs::integers::<u8>().max_value(5)) {
                0 if !names.is_empty() => {
                    let back = tc.draw_silent(gs::integers::<usize>().max_value(names.len() - 1));
                    names[names.len() - 1 - back].clone()
                }
                0 | 1 => {
                    let width = tc.draw_silent(gs::integers::<usize>().max_value(14));
                    counter += u64::from(tc.draw_silent(gs::integers::<u16>()));
                    format!("SRR{counter:0width$}.{}", counter % 3).into_bytes()
                }
                2 => {
                    let x = tc.draw_silent(gs::integers::<u32>());
                    format!("r:{}:{x}:0{x}:{:012}:{}/1", x % 9, x / 7, u64::from(x) << 20)
                        .into_bytes()
                }
                3 => {
                    let tokens = tc.draw_silent(gs::integers::<usize>().max_value(140));
                    (0..tokens).flat_map(|i| [b'a', b'0' + (i % 10) as u8]).collect()
                }
                _ => {
                    let mut bytes = tc.draw_silent(gs::binary().max_size(300));
                    for b in &mut bytes {
                        *b = (*b).max(1);
                    }
                    bytes
                }
            };
            names.push(name);
        }
        names
    }

    // r[verify cram_writer.tok3]
    #[hegel::test]
    fn every_name_list_round_trips(tc: TestCase) {
        let names = tc.draw(any_names());
        round_trip(&joined(&names));
    }

    /// On names htscodecs takes, the token streams are htscodecs' own, so
    /// compressed alike the blocks are the same size.
    // r[verify cram_writer.tok3]
    #[hegel::test(test_cases = crate::pinned::cases(300))]
    fn illumina_names_tokenise_as_htscodecs_does(tc: TestCase) {
        let names = tc.draw(illumina_names());
        let ours = round_trip(&joined(&names));
        let theirs = try_htscodecs_tokenise(&names, 9, false).unwrap();
        assert_eq!(ours.get(..9), theirs.get(..9), "header");
        assert_eq!(streams(&ours), streams(&theirs));
    }

    /// htscodecs' name test sets (Illumina, Ion Torrent, `PacBio`, ...).
    fn htscodecs_name_set(stem: &str) -> Vec<Vec<u8>> {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/hts-specs/cram/codecs");
        let gz = std::fs::read(format!("{path}/gzip/{stem}.names.gz")).unwrap();
        let mut text = Vec::new();
        flate2::read::MultiGzDecoder::new(gz.as_slice()).read_to_end(&mut text).unwrap();
        text.split(|&b| b == b'\n').filter(|n| !n.is_empty()).map(<[u8]>::to_vec).collect()
    }

    /// Every one of htscodecs' name test sets, the platform prefixes
    /// included, tokenises into htscodecs' own streams.
    // r[verify cram_writer.tok3]
    #[test]
    fn htscodecs_name_sets_tokenise_as_htscodecs_does() {
        for stem in ["01", "02", "03", "05", "08", "09", "10", "20", "nv", "nv2", "rr"] {
            let names = htscodecs_name_set(stem);
            let ours = round_trip(&joined(&names));
            let theirs = try_htscodecs_tokenise(&names, 9, false).unwrap();
            assert!(streams(&ours) == streams(&theirs), "{stem}.names");
        }
    }

    // r[verify cram_writer.tok3]
    #[test]
    fn edge_inputs() {
        // No names, which htscodecs' decoder refuses, and a last name
        // without its NUL.
        assert!(encode(b"").is_none());
        assert!(encode(b"a\0b").is_none());
        for names in [
            &b"\0"[..],
            b"\0\0\0",
            b"a\0",
            b"a\0a\0a\0",
            b"x:1:2:3:4\0x:1:2:3:4\0x:1:2:3:5\0",
            b"007\x00008\x00010\x00999\x001000\0",
            b"\xff\xfe\x80\0\xff\xfe\x81\0",
            b"\n\t \0\n\t\0",
        ] {
            round_trip(names);
        }
        // More tokens than a block has positions: the rest of the name is one.
        let long: Vec<u8> = (0..300).flat_map(|i| [b'a', b'0' + (i % 10) as u8]).collect();
        round_trip(&joined(&[long.clone(), long[2..].to_vec(), long]));
    }

    /// One name over and over is a DIFF and DUPs, and every position with
    /// a token stream has its type stream regenerated by the decoder: only
    /// the first position and the one `END` alone fills keep theirs.
    #[test]
    fn duplicates_and_regenerated_type_streams() {
        let block = round_trip(&b"HSQ1008:5:1101:1454:3564\0".repeat(50));
        let streams = streams(&block);
        assert_eq!(streams[&(0, TYPE)], [&[DIFF][..], &[DUP; 49]].concat());
        let types: Vec<u8> =
            streams.keys().filter(|&&(_, ty)| ty == TYPE).map(|&(p, _)| p).collect();
        assert_eq!(types, [0, 9]);
        assert_eq!(streams[&(9, TYPE)], [END]);
    }

    /// The name blocks of the hts-specs CRAM 3.1 conformance files, 20000
    /// names each, which htslib wrote: this encoder's streams are the ones
    /// htscodecs stored there.
    // r[verify cram_writer.tok3]
    #[test]
    fn htslib_name_blocks_tokenise_as_htscodecs_does() {
        let mut n = 0;
        for level in 2..=4 {
            let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/hts-specs/cram/3.1");
            let file = std::fs::read(format!("{path}/passed/level-{level}.cram")).unwrap();
            for block in raw_blocks(&file).iter().filter(|b| b.method == TOK3) {
                let names = tok3::decode(&block.data).unwrap();
                let ours = round_trip(&names);
                assert!(streams(&ours) == streams(&block.data), "level {level}");
                n += 1;
            }
        }
        assert_eq!(n, 4);
    }
}
