// r[impl cram.codec.tok3]
//! tok3 / Name Tokenizer codec (CRAM compression method 8).
//!
//! A block holds up to 128 token positions, each with byte streams per
//! token type, coded with rANS Nx16 or the arithmetic coder. Every stream is
//! decoded up front; then each name is decoded straight into the output by
//! walking its positions' type streams. A name's tokens are kept as a kind,
//! a value and the place their text landed in the output, all names' in one
//! arena, so a later name's MATCH copies those bytes and DELTA adds to that
//! value. `docs/spec/2-cram-1-reader.md` has the rules, including where they
//! differ from the pseudocode and from htscodecs; the production decoder
//! must agree with [`super::tok3_reference`] on every input.

// See rans.rs: lazy `ok_or_else(|| CramError::...)` keeps error construction and its
// `drop_in_place<CramError>` off the per-record path.
#![allow(
    clippy::unnecessary_lazy_evaluations,
    reason = "lazy form avoids per-call drop_in_place<CramError> on hot path"
)]

use super::codec_io::{Uint7Error, read_u8, read_u32_le, read_uint7, split_off};
use super::rans_nx16::{self, Nx16Order1Buf};
use super::reader::{CramError, check_codec_output};

/// Bridge `Uint7Error` (narrow, hot-path-friendly) to the rich `CramError`.
fn uint7_to_cram_error(e: Uint7Error) -> CramError {
    match e {
        Uint7Error::Truncated => CramError::Truncated { context: "tok3 uint7" },
        Uint7Error::Overflow => CramError::Uint7Overflow,
    }
}

/// Maximum number of names allowed in a tok3 block (htscodecs' limit).
const TOK3_NAME_COUNT_LIMIT: usize = 10_000_000;

/// How far past the header's uncompressed length the output may run
/// (htscodecs' margin; noodles writes the length without the last NUL).
const OUTPUT_SLACK: usize = 1024;

/// Maximum number of token positions (htscodecs' `MAX_TOKENS`).
const MAX_POSITIONS: usize = 128;

// Token types.
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
const N_TYPES: usize = 13;

/// The fixed-size store short tokens are written with.
const CHUNK: usize = 16;

/// `v` (below 10^8) as 8 ASCII digits, zero-padded, most significant in the
/// lowest byte: the two 4-digit halves in 32-bit lanes, split into 2-digit
/// lanes of 16 bits, then into digits, each step one multiply by a
/// reciprocal (exact for these ranges, and no lane's product reaches the
/// next lane).
#[allow(clippy::arithmetic_side_effects, reason = "lane values below 10^4: no product overflows")]
#[inline]
fn eight_digits(v: u32) -> u64 {
    let halves = u64::from(v / 10_000) | (u64::from(v % 10_000) << 32);
    let hundreds = ((halves * 10_486) >> 20) & 0x0000_007F_0000_007F;
    let pairs = hundreds | ((halves - hundreds * 100) << 16);
    let tens = ((pairs * 103) >> 10) & 0x000F_000F_000F_000F;
    (tens | ((pairs - tens * 10) << 8)) + 0x3030_3030_3030_3030
}

/// `value` in decimal, left-padded with zeros to `width` (at most
/// [`CHUNK`]), at the start of a chunk; and its length. Bytes past the
/// length are unspecified.
#[allow(
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    reason = "len <= CHUNK; below 8 the shift is under 64; value / 10^8 is at most 42"
)]
#[inline(always)]
fn format_padded(value: u32, width: usize) -> ([u8; CHUNK], usize) {
    let digits = value.checked_ilog10().map_or(1, |l| l as usize + 1);
    let len = width.max(digits).min(CHUNK);
    let mut chunk = [b'0'; CHUNK];
    let low = eight_digits(value % 100_000_000);
    if len < 8 {
        // Then `value` has at most `len` digits: the last `len` of the eight.
        if let Some(dst) = chunk.first_chunk_mut::<8>() {
            *dst = (low >> (8 * (8 - len))).to_le_bytes();
        }
    } else {
        if let Some(dst) = chunk.get_mut(len - 8..len) {
            dst.copy_from_slice(&low.to_le_bytes());
        }
        // A ninth and tenth digit.
        let top = value / 100_000_000;
        if top >= 10 {
            if let Some(dst) = chunk.get_mut(len - 10..len - 8) {
                dst.copy_from_slice(&[b'0' + (top / 10) as u8, b'0' + (top % 10) as u8]);
            }
        } else if top > 0
            && let Some(d) = chunk.get_mut(len - 9)
        {
            *d = b'0' + top as u8;
        }
    }
    (chunk, len)
}

/// Each token type as a byte, so a regenerated type stream's first byte
/// can be a slice of it.
const TYPE_BYTES: [u8; N_TYPES] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];

/// Decode a tok3 compressed block.
///
/// Builds rANS Nx16 order-1 tables afresh for every block that needs them;
/// a [`Decoder`] keeps them.
pub fn decode(src: &[u8]) -> Result<Vec<u8>, CramError> {
    Decoder::new().decode(src)
}

/// A tok3 decoder that keeps its rANS Nx16 order-1 tables (about 5 MiB,
/// allocated on the first block that needs them) from one block to the
/// next, as the CRAM reader does. Building them per block costs more than
/// decoding a small block.
// r[impl cram.codec.tok3.table_reuse]
#[derive(Default)]
pub struct Decoder {
    rans_buf: Nx16Order1Buf,
}

impl Decoder {
    pub fn new() -> Self {
        Self::default()
    }

    /// Decode a tok3 compressed block.
    pub fn decode(&mut self, src: &[u8]) -> Result<Vec<u8>, CramError> {
        decode_with_buf(src, &mut self.rans_buf)
    }
}

/// [`decode`], reusing the caller's rANS Nx16 order-1 tables.
// r[impl cram.codec.tok3.table_reuse]
pub(crate) fn decode_with_buf(
    src: &[u8],
    rans_buf: &mut Nx16Order1Buf,
) -> Result<Vec<u8>, CramError> {
    // Every stream stores its length, so the size passed is unused.
    decode_with(src, |stream| super::arith::decode(stream, 0), rans_buf)
}

/// [`decode`] with the arithmetic coder's streams decoded by
/// [`super::arith::reference`], so tests and the fuzzer can compare it with
/// [`decode`].
#[cfg(any(test, feature = "fuzz"))]
#[doc(hidden)]
pub fn decode_with_arith_reference(src: &[u8]) -> Result<Vec<u8>, CramError> {
    decode_with(
        src,
        |stream| {
            super::arith::reference::decode(stream, 0)
                .ok_or(CramError::ArithCorruptData { context: "arith reference decoder" })
        },
        &mut Nx16Order1Buf::new(),
    )
}

/// Decode a tok3 block, decoding arith-coded streams with `arith` and
/// rANS Nx16 ones with `rans_buf`'s order-1 tables.
fn decode_with(
    src: &[u8],
    arith: impl Fn(&[u8]) -> Result<Vec<u8>, CramError>,
    rans_buf: &mut Nx16Order1Buf,
) -> Result<Vec<u8>, CramError> {
    let mut cur: &[u8] = src;
    let truncated = || CramError::Truncated { context: "tok3 header" };
    let ulen = read_u32_le(&mut cur).ok_or_else(truncated)? as usize;
    let name_count = read_u32_le(&mut cur).ok_or_else(truncated)? as usize;
    let use_arith = read_u8(&mut cur).ok_or_else(truncated)? != 0;

    // r[impl cram.tok3.name_count_limit]
    // r[impl cram.codec.tok3.limits]
    if name_count > TOK3_NAME_COUNT_LIMIT {
        return Err(CramError::Tok3NameCountExceedsLimit {
            count: name_count,
            limit: TOK3_NAME_COUNT_LIMIT,
        });
    }
    check_codec_output(ulen, "tok3 output")?;
    let max_output = ulen.saturating_add(OUTPUT_SLACK);
    // Every name takes at least its NUL.
    if name_count > max_output {
        return Err(CramError::Tok3NameCountExceedsLength { count: name_count, length: ulen });
    }

    let streams = Streams::read(&mut cur, use_arith.then_some(&arith), rans_buf)?;
    let mut names = Names::new(streams.cursors(name_count), max_output, name_count);
    for n in 0..name_count {
        names.decode_name(n)?;
    }
    Ok(names.finish())
}

// ── Token streams ────────────────────────────────────────────────────

/// Where a token stream's bytes come from.
#[derive(Clone, Copy)]
enum Source {
    /// `Streams::decoded[i]`.
    Decoded(usize),
    /// A regenerated TYPE stream: this type, then a MATCH for every other
    /// name.
    Regenerated(u8),
}

/// The decoded token streams, per position and type.
struct Streams {
    decoded: Vec<Vec<u8>>,
    positions: Vec<[Option<Source>; N_TYPES]>,
}

impl Streams {
    // r[impl cram.codec.tok3.streams]
    /// Read and decode every token stream; `arith` decodes them when the
    /// block uses the arithmetic coder, else they are rANS Nx16, decoded
    /// with `rans_buf`'s order-1 tables.
    fn read(
        src: &mut &[u8],
        arith: Option<&impl Fn(&[u8]) -> Result<Vec<u8>, CramError>>,
        rans_buf: &mut Nx16Order1Buf,
    ) -> Result<Self, CramError> {
        let mut streams = Self { decoded: Vec::new(), positions: Vec::new() };

        while let Some(ttype) = read_u8(src) {
            let ty = ttype & 0x3f;
            if usize::from(ty) >= N_TYPES {
                return Err(CramError::InvalidTok3TokenType { token_type: ttype });
            }
            if ttype & 0x80 != 0 {
                if streams.positions.len() >= MAX_POSITIONS {
                    return Err(CramError::Tok3TooManyPositions { limit: MAX_POSITIONS });
                }
                let mut position = [None; N_TYPES];
                if ty != TYPE
                    && let Some(types) = position.get_mut(usize::from(TYPE))
                {
                    *types = Some(Source::Regenerated(ty));
                }
                streams.positions.push(position);
            }
            if streams.positions.is_empty() {
                return Err(CramError::Truncated {
                    context: "tok3 stream before the first position",
                });
            }

            let source = if ttype & 0x40 != 0 {
                let truncated = || CramError::Truncated { context: "tok3 dup metadata" };
                let position = usize::from(read_u8(src).ok_or_else(truncated)?);
                let token_type = read_u8(src).ok_or_else(truncated)?;
                streams
                    .positions
                    .get(position)
                    .and_then(|p| p.get(usize::from(token_type)).copied().flatten())
                    .ok_or_else(|| CramError::Tok3DupStreamUnset { position, token_type })?
            } else {
                let compressed_size = read_uint7(src).map_err(uint7_to_cram_error)? as usize;
                let data = split_off(src, compressed_size)
                    .ok_or_else(|| CramError::Truncated { context: "tok3 compressed payload" })?;
                // r[impl cram.codec.tok3_arith]
                let decoded = match arith {
                    Some(arith) => arith(data)?,
                    None => rans_nx16::decode_with_buf(data, 0, rans_buf)?,
                };
                streams.decoded.push(decoded);
                Source::Decoded(streams.decoded.len().saturating_sub(1))
            };
            if let Some(slot) =
                streams.positions.last_mut().and_then(|p| p.get_mut(usize::from(ty)))
            {
                *slot = Some(source);
            }
        }
        Ok(streams)
    }

    /// A read cursor for every stream, unset ones empty.
    fn cursors(&self, name_count: usize) -> Vec<Position<'_>> {
        self.positions
            .iter()
            .map(|position| {
                Position(position.map(|source| match source {
                    None => Cursor::default(),
                    Some(Source::Decoded(i)) => {
                        Cursor { head: self.decoded.get(i).map_or(&[], Vec::as_slice), matches: 0 }
                    }
                    Some(Source::Regenerated(ty)) => Cursor {
                        head: TYPE_BYTES.get(usize::from(ty)..=usize::from(ty)).unwrap_or(&[]),
                        matches: name_count.saturating_sub(1),
                    },
                }))
            })
            .collect()
    }
}

/// One token position's streams, indexed by token type.
struct Position<'a>([Cursor<'a>; N_TYPES]);

impl<'a> Position<'a> {
    /// The stream of type `ty`; the constant types the decoder asks for
    /// always exist, so the check folds away.
    #[inline]
    fn stream(&mut self, ty: u8) -> Result<&mut Cursor<'a>, CramError> {
        self.0
            .get_mut(usize::from(ty))
            .ok_or_else(|| CramError::InvalidTok3TokenType { token_type: ty })
    }
}

/// Reads one token stream from the front: its bytes, then (for a
/// regenerated TYPE stream) `matches` MATCH bytes.
#[derive(Clone, Copy, Default)]
struct Cursor<'a> {
    head: &'a [u8],
    matches: usize,
}

impl<'a> Cursor<'a> {
    #[inline]
    fn u8(&mut self) -> Result<u8, CramError> {
        if let Some((&b, rest)) = self.head.split_first() {
            self.head = rest;
            Ok(b)
        } else {
            self.matches = self
                .matches
                .checked_sub(1)
                .ok_or_else(|| CramError::Truncated { context: "tok3 token stream" })?;
            Ok(MATCH)
        }
    }

    #[inline]
    fn u32(&mut self) -> Result<u32, CramError> {
        if let Some((bytes, rest)) = self.head.split_first_chunk::<4>() {
            self.head = rest;
            return Ok(u32::from_le_bytes(*bytes));
        }
        Ok(u32::from_le_bytes([self.u8()?, self.u8()?, self.u8()?, self.u8()?]))
    }

    /// Bytes up to a NUL, which is consumed. A regenerated stream's MATCH
    /// bytes hold no NUL, so the string must end in its bytes.
    #[inline]
    fn string(&mut self) -> Result<&'a [u8], CramError> {
        let nul = self
            .head
            .iter()
            .position(|&b| b == 0)
            .ok_or_else(|| CramError::Truncated { context: "tok3 unterminated string" })?;
        let (s, rest) = self.head.split_at(nul);
        self.head = rest.get(1..).unwrap_or_default();
        Ok(s)
    }
}

// ── Names ────────────────────────────────────────────────────────────

/// What a token holds, for a later name's MATCH, DELTA and DELTA0.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(u8)]
enum Kind {
    /// CHAR or STRING.
    Text = 1,
    /// DIGITS or DELTA.
    Digits = 2,
    /// DIGITS0 or DELTA0.
    Digits0 = 3,
    /// NOP or END.
    Empty = 4,
}

/// A decoded token: its kind, its numeric value, and where its text is in
/// the output.
#[derive(Clone, Copy)]
struct Token {
    kind: Kind,
    value: u32,
    start: u32,
    len: u32,
}

impl Token {
    const EMPTY: Self = Self { kind: Kind::Empty, value: 0, start: 0, len: 0 };
}

/// Token `t` (from 1) of a name whose tokens start at `first` and number
/// `count`.
#[inline]
fn token_at(tokens: &[Token], first: usize, count: usize, t: usize) -> Option<Token> {
    if t > count {
        return None;
    }
    tokens.get(first.wrapping_add(t).wrapping_sub(1)).copied()
}

/// A decoded name: its text in the output (without the NUL) and its tokens
/// in the arena.
#[derive(Clone, Copy)]
struct Name {
    start: u32,
    len: u32,
    first_token: u32,
    tokens: u32,
}

/// Decodes names into the output, keeping every name's tokens.
struct Names<'a> {
    positions: Vec<Position<'a>>,
    out: Vec<u8>,
    /// Bytes of `out` written; `out` is as long as the output may get, so a
    /// write past it is an error rather than a reallocation.
    len: usize,
    names: Vec<Name>,
    tokens: Vec<Token>,
}

/// A `usize` below the output length (at most `MAX_ALLOC_SIZE`) or the
/// token count as the `u32` it is stored as.
fn to_u32(v: usize) -> Result<u32, CramError> {
    u32::try_from(v).map_err(|_| CramError::Truncated { context: "tok3 offset overflow" })
}

impl<'a> Names<'a> {
    fn new(positions: Vec<Position<'a>>, max_output: usize, name_count: usize) -> Self {
        Self {
            positions,
            out: vec![0; max_output],
            len: 0,
            names: Vec::with_capacity(name_count),
            tokens: Vec::new(),
        }
    }

    fn finish(mut self) -> Vec<u8> {
        self.out.truncate(self.len);
        self.out
    }

    /// The next `n` output bytes.
    #[inline]
    fn reserve(&mut self, n: usize) -> Result<&mut [u8], CramError> {
        let end = self.len.checked_add(n);
        let limit = self.out.len();
        let dst = end
            .and_then(|end| self.out.get_mut(self.len..end))
            .ok_or_else(|| CramError::Tok3OutputOverflow { limit })?;
        self.len = self.len.wrapping_add(n);
        Ok(dst)
    }

    /// Write `bytes` as a text token.
    #[inline]
    fn put_text(&mut self, bytes: &[u8]) -> Result<Token, CramError> {
        let start = to_u32(self.len)?;
        self.reserve(bytes.len())?.copy_from_slice(bytes);
        Ok(Token { kind: Kind::Text, value: 0, start, len: to_u32(bytes.len())? })
    }

    /// Write `value` in decimal, left-padded with zeros to `width`.
    #[inline(always)]
    fn put_digits(&mut self, kind: Kind, value: u32, width: usize) -> Result<Token, CramError> {
        let start = self.len;
        if width <= CHUNK {
            let (chunk, len) = format_padded(value, width);
            self.put_chunk(&chunk, len)?;
        } else {
            // DZLEN up to 255: a rare, slow path.
            let (chunk, digits) = format_padded(value, 0);
            let dst = self.reserve(width)?;
            let (pad, rest) = dst.split_at_mut(width.saturating_sub(digits));
            pad.fill(b'0');
            rest.copy_from_slice(chunk.get(..digits).unwrap_or_default());
        }
        Ok(Token { kind, value, start: to_u32(start)?, len: to_u32(self.len.wrapping_sub(start))? })
    }

    /// Write the first `len` bytes of `chunk` (`len <= CHUNK`): all of it
    /// in one store when there is room (bytes past `len` are overwritten
    /// by what comes next, or cut off at the end).
    #[inline(always)]
    fn put_chunk(&mut self, chunk: &[u8; CHUNK], len: usize) -> Result<(), CramError> {
        let end = self.len.wrapping_add(len);
        if let Some(dst) = self.out.get_mut(self.len..).and_then(|d| d.first_chunk_mut::<CHUNK>()) {
            *dst = *chunk;
            if end > self.out.len() {
                return Err(CramError::Tok3OutputOverflow { limit: self.out.len() });
            }
        } else {
            let limit = self.out.len();
            let dst = self
                .out
                .get_mut(self.len..end)
                .ok_or_else(|| CramError::Tok3OutputOverflow { limit })?;
            dst.copy_from_slice(chunk.get(..len).unwrap_or_default());
        }
        self.len = end;
        Ok(())
    }

    /// Write again the text of `token`, from an earlier name.
    #[inline(always)]
    fn put_copy(&mut self, token: Token) -> Result<(), CramError> {
        let (start, len) = (token.start as usize, token.len as usize);
        // Short text: one 16-byte load and store. The load may run past
        // the written bytes; what it brings along lands past `len`.
        if len <= CHUNK
            && let Some(&chunk) = self.out.get(start..).and_then(|s| s.first_chunk::<CHUNK>())
        {
            return self.put_chunk(&chunk, len);
        }
        let (written, free) = self.out.split_at_mut(self.len);
        let src = written.get(start..start.wrapping_add(len));
        let dst = free.get_mut(..len);
        match (src, dst) {
            (Some(src), Some(dst)) => dst.copy_from_slice(src),
            _ => return Err(CramError::Tok3OutputOverflow { limit: self.out.len() }),
        }
        self.len = self.len.wrapping_add(len);
        Ok(())
    }

    // r[impl cram.codec.tok3.names]
    /// Decode name `n` and its NUL into the output.
    fn decode_name(&mut self, n: usize) -> Result<(), CramError> {
        let start = to_u32(self.len)?;
        let position0 = self
            .positions
            .first_mut()
            .ok_or_else(|| CramError::Truncated { context: "tok3 no token positions" })?;
        let ty = position0.stream(TYPE)?.u8()?;
        let distance = position0.stream(ty)?.u32()? as usize;
        let m = n
            .checked_sub(distance)
            .ok_or_else(|| CramError::Tok3DistanceExceedsIndex { distance, name_index: n })?;
        // `m == n` has no name yet: nothing to copy, match or add to.
        let previous = self.names.get(m).copied();

        if ty == DUP {
            let previous = previous.ok_or_else(|| CramError::Tok3DupRefOutOfRange { index: m })?;
            self.put_copy(Token { len: previous.len, start: previous.start, ..Token::EMPTY })?;
            self.reserve(1)?.fill(0);
            self.names.push(Name { start, ..previous });
            return Ok(());
        }

        let first_token = self.tokens.len();
        let (prev_first, prev_count) =
            previous.map_or((0, 0), |p| (p.first_token as usize, p.tokens as usize));
        let mut t = 1;
        loop {
            let position = self.positions.get_mut(t).ok_or_else(|| CramError::Truncated {
                context: "tok3 name past the last position",
            })?;
            let ty = position.stream(TYPE)?.u8()?;
            // The previous name's token at this position.
            let prev = || token_at(&self.tokens, prev_first, prev_count, t);
            let token = match ty {
                CHAR => {
                    let c = position.stream(CHAR)?.u8()?;
                    self.put_text(&[c])?
                }
                STRING => {
                    let s = position.stream(STRING)?.string()?;
                    self.put_text(s)?
                }
                DIGITS => {
                    let value = position.stream(DIGITS)?.u32()?;
                    self.put_digits(Kind::Digits, value, 0)?
                }
                DIGITS0 => {
                    let value = position.stream(DIGITS0)?.u32()?;
                    let width = position.stream(DZLEN)?.u8()?;
                    self.put_digits(Kind::Digits0, value, usize::from(width))?
                }
                DELTA => {
                    let delta = position.stream(DELTA)?.u8()?;
                    match prev() {
                        Some(p) if p.kind == Kind::Digits => self.put_digits(
                            Kind::Digits,
                            p.value.wrapping_add(u32::from(delta)),
                            0,
                        )?,
                        p => {
                            return Err(CramError::Tok3DeltaRequiresDigits {
                                found: p.map_or(0, |p| p.kind as u8),
                            });
                        }
                    }
                }
                DELTA0 => {
                    let delta = position.stream(DELTA0)?.u8()?;
                    match prev() {
                        Some(p) if p.kind == Kind::Digits0 => self.put_digits(
                            Kind::Digits0,
                            p.value.wrapping_add(u32::from(delta)),
                            p.len as usize,
                        )?,
                        p => {
                            return Err(CramError::Tok3Delta0RequiresPaddedDigits {
                                found: p.map_or(0, |p| p.kind as u8),
                            });
                        }
                    }
                }
                MATCH => match prev() {
                    Some(p) if p.kind != Kind::Empty => {
                        self.put_copy(p)?;
                        p
                    }
                    _ => return Err(CramError::Tok3MatchWithoutValue { position: t }),
                },
                NOP | END => Token::EMPTY,
                _ => return Err(CramError::InvalidTok3TokenType { token_type: ty }),
            };
            self.tokens.push(token);
            if ty == END {
                break;
            }
            t = t.wrapping_add(1);
        }

        let len = to_u32(self.len)?.wrapping_sub(start);
        self.reserve(1)?.fill(0);
        let tokens = to_u32(self.tokens.len().wrapping_sub(first_token))?;
        self.names.push(Name { start, len, first_token: to_u32(first_token)?, tokens });
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // r[verify cram.codec.tok3]

    /// One part of a hand-built block: a raw (CAT) stream, or a copy.
    enum Part<'a> {
        Stream(u8, &'a [u8]),
        Dup(u8, u8, u8),
    }

    /// A tok3 block over rANS Nx16 CAT streams.
    fn block(ulen: u32, names: u32, parts: &[Part<'_>]) -> Vec<u8> {
        let mut out = Vec::new();
        out.extend(ulen.to_le_bytes());
        out.extend(names.to_le_bytes());
        out.push(0);
        for part in parts {
            match *part {
                Part::Stream(ttype, data) => {
                    let mut stream = vec![0x20];
                    put_uint7(&mut stream, data.len());
                    stream.extend(data);
                    out.push(ttype);
                    put_uint7(&mut out, stream.len());
                    out.extend(stream);
                }
                Part::Dup(ttype, pos, ty) => out.extend([ttype | 0x40, pos, ty]),
            }
        }
        out
    }

    const NEW: u8 = 0x80;

    fn u32s(values: &[u32]) -> Vec<u8> {
        values.iter().flat_map(|v| v.to_le_bytes()).collect()
    }

    /// Two names, both DIFF against the one before (the first against
    /// itself), whose position 1 has the given type stream and value
    /// streams, and position 2 ends them.
    fn two_names(ulen: u32, types: &[u8], values: &[Part<'_>]) -> Vec<u8> {
        let dist = u32s(&[0, 1]);
        let mut parts = vec![
            Part::Stream(NEW | TYPE, &[6, 6]),
            Part::Stream(6, &dist),
            Part::Stream(NEW | TYPE, types),
        ];
        parts.extend(values.iter().map(|p| match *p {
            Part::Stream(t, d) => Part::Stream(t, d),
            Part::Dup(t, p, y) => Part::Dup(t, p, y),
        }));
        parts.push(Part::Stream(NEW | TYPE, &[END; 3]));
        block(ulen, 2, &parts)
    }

    /// The fixed-width digit writer is `format!`'s zero padding, at every
    /// width it takes.
    #[hegel::test]
    fn format_padded_matches_format(tc: hegel::TestCase) {
        use hegel::generators as gs;
        let value = if tc.draw(gs::booleans()) {
            tc.draw(gs::integers::<u32>())
        } else {
            tc.draw(gs::sampled_from(
                &[0, 9, 10, 99, 100, 999_999_999, 1_000_000_000, u32::MAX][..],
            ))
        };
        let width = tc.draw(gs::integers::<usize>().max_value(CHUNK));
        let (chunk, len) = format_padded(value, width);
        assert_eq!(&chunk[..len], format!("{value:0width$}").as_bytes());
    }

    /// Every digit count at every width, and a spread of values.
    #[test]
    fn format_padded_matches_format_at_the_edges() {
        let mut values: Vec<u32> = (0..=10_000).collect();
        values.extend((0..=9).flat_map(|k| {
            let p = 10u32.pow(k);
            [p - 1, p, p + 1, p.saturating_mul(5)]
        }));
        values.extend([u32::MAX, u32::MAX - 1, 4_000_000_000, 999_999_999, 1_000_000_000]);
        values.extend((0..100_000u32).map(|i| i.wrapping_mul(2_654_435_761)));
        for &value in &values {
            for width in 0..=CHUNK {
                let (chunk, len) = format_padded(value, width);
                assert_eq!(&chunk[..len], format!("{value:0width$}").as_bytes(), "{value} {width}");
            }
        }
    }

    // r[verify cram.codec.tok3.names]
    #[test]
    fn delta_wraps_at_32_bits() {
        let digits = u32s(&[u32::MAX]);
        let src = two_names(
            100,
            &[DIGITS, DELTA],
            &[Part::Stream(DIGITS, &digits), Part::Stream(DELTA, &[1])],
        );
        assert_eq!(decode(&src).unwrap(), b"4294967295\x000\0");
    }

    // r[verify cram.codec.tok3.names]
    #[test]
    fn digits0_pads_and_never_truncates() {
        let digits = u32s(&[123, 123]);
        let src = two_names(
            100,
            &[DIGITS0, DIGITS0],
            &[Part::Stream(DIGITS0, &digits), Part::Stream(DZLEN, &[5, 2])],
        );
        assert_eq!(decode(&src).unwrap(), b"00123\x00123\0");
        // DELTA0 pads to the previous token's length, itself at least its
        // digit count.
        let digits = u32s(&[98]);
        let src = two_names(
            100,
            &[DIGITS0, DELTA0],
            &[
                Part::Stream(DIGITS0, &digits),
                Part::Stream(DZLEN, &[1]),
                Part::Stream(DELTA0, &[3]),
            ],
        );
        assert_eq!(decode(&src).unwrap(), b"98\x00101\0");
        let digits = u32s(&[98]);
        let src = two_names(
            100,
            &[DIGITS0, DELTA0],
            &[
                Part::Stream(DIGITS0, &digits),
                Part::Stream(DZLEN, &[4]),
                Part::Stream(DELTA0, &[3]),
            ],
        );
        assert_eq!(decode(&src).unwrap(), b"0098\x000101\0");
    }

    // r[verify cram.codec.tok3.names]
    #[test]
    fn deltas_and_matches_need_the_right_previous_token() {
        let src =
            two_names(100, &[CHAR, DELTA], &[Part::Stream(CHAR, b"a"), Part::Stream(DELTA, &[1])]);
        assert!(matches!(decode(&src), Err(CramError::Tok3DeltaRequiresDigits { found: 1 })));
        let digits = u32s(&[5]);
        let src = two_names(
            100,
            &[DIGITS, DELTA0],
            &[Part::Stream(DIGITS, &digits), Part::Stream(DELTA0, &[1])],
        );
        assert!(matches!(
            decode(&src),
            Err(CramError::Tok3Delta0RequiresPaddedDigits { found: 2 })
        ));
        let src = two_names(100, &[NOP, MATCH], &[]);
        assert!(matches!(decode(&src), Err(CramError::Tok3MatchWithoutValue { position: 1 })));
        // The first name has no previous name to match.
        let src = two_names(100, &[MATCH, MATCH], &[]);
        assert!(matches!(decode(&src), Err(CramError::Tok3MatchWithoutValue { position: 1 })));
        // A MATCH of a DELTA copies the value it produced.
        let digits = u32s(&[7]);
        let src = block(
            100,
            3,
            &[
                Part::Stream(NEW | TYPE, &[6, 6, 6]),
                Part::Stream(6, &u32s(&[0, 1, 1])),
                Part::Stream(NEW | TYPE, &[DIGITS, DELTA, MATCH]),
                Part::Stream(DIGITS, &digits),
                Part::Stream(DELTA, &[2]),
                Part::Stream(NEW | TYPE, &[END; 3]),
            ],
        );
        assert_eq!(decode(&src).unwrap(), b"7\x009\x009\0");
    }

    // r[verify cram.codec.tok3.names]
    #[test]
    fn invalid_names_are_rejected() {
        // An unterminated string.
        let src = two_names(100, &[STRING, STRING], &[Part::Stream(STRING, b"ab\0cd")]);
        assert!(matches!(decode(&src), Err(CramError::Truncated { .. })));
        // A DZLEN type inside a name.
        let src = two_names(100, &[NOP, DZLEN], &[]);
        assert!(matches!(decode(&src), Err(CramError::InvalidTok3TokenType { token_type: 4 })));
        // A distance past the first name.
        let dist = u32s(&[1]);
        let src =
            block(100, 1, &[Part::Stream(NEW | 6, &dist), Part::Stream(NEW | TYPE, &[END; 3])]);
        assert!(matches!(
            decode(&src),
            Err(CramError::Tok3DistanceExceedsIndex { distance: 1, name_index: 0 })
        ));
        // A name duplicating itself.
        let dist = u32s(&[0]);
        let src =
            block(100, 1, &[Part::Stream(NEW | DUP, &dist), Part::Stream(NEW | TYPE, &[END; 3])]);
        assert!(matches!(decode(&src), Err(CramError::Tok3DupRefOutOfRange { index: 0 })));
    }

    // r[verify cram.codec.tok3.streams]
    #[test]
    fn stream_headers_are_checked() {
        let src = block(10, 1, &[Part::Stream(NEW | 13, &[])]);
        assert!(matches!(decode(&src), Err(CramError::InvalidTok3TokenType { token_type: 0x8D })));
        let src = block(10, 1, &[Part::Dup(NEW, 99, 0)]);
        assert!(matches!(
            decode(&src),
            Err(CramError::Tok3DupStreamUnset { position: 99, token_type: 0 })
        ));
        let src = block(10, 1, &[Part::Stream(TYPE, &[6])]);
        assert!(matches!(decode(&src), Err(CramError::Truncated { .. })));
        let parts: Vec<Part<'_>> = (0..129).map(|_| Part::Stream(NEW | NOP, &[])).collect();
        assert!(matches!(
            decode(&block(10, 1, &parts)),
            Err(CramError::Tok3TooManyPositions { limit: 128 })
        ));
    }

    /// A position whose first stream is not its type stream gets one
    /// regenerated: that type for the first name, MATCH for the rest; and a
    /// copied stream reads from its own start.
    // r[verify cram.codec.tok3.streams]
    #[test]
    fn regenerated_and_copied_streams() {
        let dist = u32s(&[0, 1, 1]);
        let src = block(
            100,
            3,
            &[
                Part::Stream(NEW | TYPE, &[6; 3]),
                Part::Stream(6, &dist),
                Part::Stream(NEW | STRING, b"read\0"),
                Part::Stream(NEW | CHAR, b"x"),
                Part::Dup(NEW | CHAR, 2, CHAR),
                Part::Stream(NEW | TYPE, &[END; 3]),
            ],
        );
        assert_eq!(decode(&src).unwrap(), b"readxx\0readxx\0readxx\0");
    }

    // r[verify cram.codec.tok3.limits]
    #[test]
    fn output_is_bounded_by_the_declared_length() {
        // 1 KiB past the declared length: one name of 1021 + 2 + 1 bytes
        // fits a declared length of 0, and of 1022 + 2 + 1 does not.
        let long = |n: usize| {
            let mut s = vec![b'a'; n];
            s.push(0);
            let dist = u32s(&[0]);
            block(
                0,
                1,
                &[
                    Part::Stream(NEW | TYPE, &[6]),
                    Part::Stream(6, &dist),
                    Part::Stream(NEW | STRING, &s),
                    Part::Stream(NEW | DIGITS0, &u32s(&[7])),
                    Part::Stream(DZLEN, &[2]),
                    Part::Stream(NEW | TYPE, &[END]),
                ],
            )
        };
        assert_eq!(decode(&long(1021)).unwrap().len(), 1024);
        assert!(matches!(decode(&long(1022)), Err(CramError::Tok3OutputOverflow { limit: 1024 })));
        let src = block(1, 1026, &[]);
        assert!(matches!(
            decode(&src),
            Err(CramError::Tok3NameCountExceedsLength { count: 1026, length: 1 })
        ));
        // No names decode to nothing, whatever the streams.
        assert_eq!(decode(&block(0, 0, &[])).unwrap(), b"");
    }

    // r[verify cram.codec.tok3_arith]
    #[test]
    fn decode_tok3_noodles_test_vector() {
        let src = [
            0x58, 0x00, 0x00, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x80, 0x15, 0x00, 0x03, 0x06,
            0x00, 0x04, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00,
            0x00, 0x80, 0x00, 0x00, 0x06, 0x18, 0x00, 0x0c, 0x00, 0x01, 0x00, 0x00, 0x0e, 0x02,
            0x00, 0x22, 0x25, 0x00, 0x00, 0xbc, 0x00, 0x00, 0x00, 0xbc, 0x00, 0x00, 0x00, 0xbc,
            0x00, 0x00, 0x80, 0x17, 0x00, 0x03, 0x01, 0x0a, 0x00, 0x01, 0x03, 0x00, 0x00, 0x02,
            0x00, 0x00, 0xac, 0x00, 0x00, 0x00, 0xac, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x01,
            0x1b, 0x00, 0x04, 0x00, 0x31, 0x37, 0x49, 0x00, 0x01, 0x01, 0x01, 0x01, 0x00, 0x0c,
            0x02, 0x00, 0x00, 0x04, 0x02, 0x00, 0x00, 0x08, 0x02, 0x00, 0x00, 0x00, 0x02, 0x00,
            0x80, 0x17, 0x00, 0x03, 0x02, 0x0a, 0x00, 0x01, 0x03, 0x00, 0x00, 0x02, 0x00, 0x00,
            0xac, 0x00, 0x00, 0x00, 0xac, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x02, 0x15, 0x00,
            0x01, 0x5f, 0x00, 0x01, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80,
            0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x80, 0x17, 0x00, 0x03, 0x03, 0x0a, 0x00, 0x01,
            0x03, 0x00, 0x00, 0x02, 0x00, 0x00, 0xac, 0x00, 0x00, 0x00, 0xac, 0x00, 0x00, 0x00,
            0x80, 0x00, 0x00, 0x03, 0x19, 0x00, 0x04, 0x00, 0x22, 0x3d, 0x00, 0x02, 0x01, 0x01,
            0x00, 0x0c, 0x02, 0x00, 0x00, 0x08, 0x02, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00,
            0x01, 0x00, 0x04, 0x15, 0x00, 0x01, 0x05, 0x00, 0x01, 0x00, 0x80, 0x00, 0x00, 0x00,
            0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x80, 0x17, 0x00,
            0x03, 0x02, 0x0a, 0x00, 0x01, 0x03, 0x00, 0x00, 0x02, 0x00, 0x00, 0xac, 0x00, 0x00,
            0x00, 0xac, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x02, 0x15, 0x00, 0x01, 0x3a, 0x00,
            0x01, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00,
            0x80, 0x00, 0x00, 0x80, 0x17, 0x00, 0x03, 0x07, 0x0a, 0x00, 0x01, 0x03, 0x00, 0x00,
            0x02, 0x00, 0x00, 0xac, 0x00, 0x00, 0x00, 0xac, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00,
            0x07, 0x17, 0x00, 0x04, 0x00, 0x02, 0x00, 0x03, 0x01, 0x00, 0x0c, 0x02, 0x00, 0x00,
            0xa8, 0x00, 0x00, 0x00, 0xa8, 0x00, 0x00, 0x00, 0xa8, 0x00, 0x00, 0x80, 0x17, 0x00,
            0x03, 0x02, 0x0a, 0x00, 0x01, 0x03, 0x00, 0x00, 0x02, 0x00, 0x00, 0xac, 0x00, 0x00,
            0x00, 0xac, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x02, 0x15, 0x00, 0x01, 0x3a, 0x00,
            0x01, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00,
            0x80, 0x00, 0x00, 0x80, 0x17, 0x00, 0x03, 0x07, 0x0a, 0x00, 0x03, 0x01, 0x00, 0xa8,
            0x00, 0x00, 0x00, 0x0c, 0x02, 0x00, 0x00, 0xa8, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00,
            0x07, 0x1a, 0x00, 0x08, 0x00, 0x7b, 0x7c, 0x00, 0x00, 0x06, 0x01, 0x01, 0x00, 0x7c,
            0x20, 0x00, 0x00, 0xe0, 0x00, 0x00, 0x00, 0xe0, 0x00, 0x00, 0x00, 0xe0, 0x00, 0x00,
            0x80, 0x17, 0x00, 0x03, 0x02, 0x0a, 0x00, 0x01, 0x03, 0x00, 0x00, 0x02, 0x00, 0x00,
            0xac, 0x00, 0x00, 0x00, 0xac, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x02, 0x15, 0x00,
            0x01, 0x3a, 0x00, 0x01, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80,
            0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x80, 0x15, 0x00, 0x03, 0x07, 0x00, 0x04, 0x00,
            0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00,
            0x00, 0x07, 0x22, 0x00, 0x0c, 0x00, 0x06, 0x2d, 0x64, 0x65, 0x00, 0xb2, 0xf0, 0x00,
            0x0a, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x00, 0xcd, 0x0b, 0x08, 0x00, 0xaf, 0x0e,
            0x08, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x02, 0x00, 0x80, 0x17, 0x00, 0x03, 0x02,
            0x0a, 0x00, 0x01, 0x03, 0x00, 0x00, 0x02, 0x00, 0x00, 0xac, 0x00, 0x00, 0x00, 0xac,
            0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x02, 0x15, 0x00, 0x01, 0x3a, 0x00, 0x01, 0x00,
            0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00,
            0x00, 0x80, 0x17, 0x00, 0x03, 0x03, 0x07, 0x00, 0x03, 0x01, 0x00, 0xa8, 0x00, 0x00,
            0x00, 0xa8, 0x00, 0x00, 0x00, 0x0c, 0x02, 0x00, 0x00, 0x80, 0x00, 0x00, 0x03, 0x1d,
            0x00, 0x08, 0x00, 0x06, 0x21, 0xa3, 0xe3, 0x00, 0x04, 0x01, 0x01, 0x01, 0x01, 0x00,
            0x6e, 0x20, 0x00, 0x00, 0x58, 0x20, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x02,
            0x00, 0x04, 0x15, 0x00, 0x02, 0x05, 0x00, 0x02, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80,
            0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x07, 0x19, 0x00, 0x04,
            0x00, 0x21, 0x3f, 0x00, 0x02, 0x01, 0x01, 0x00, 0x08, 0x02, 0x00, 0x00, 0x0c, 0x02,
            0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01, 0x00, 0x80, 0x17, 0x00, 0x03, 0x02,
            0x0a, 0x00, 0x01, 0x03, 0x00, 0x00, 0x02, 0x00, 0x00, 0xac, 0x00, 0x00, 0x00, 0xac,
            0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x02, 0x15, 0x00, 0x01, 0x23, 0x00, 0x01, 0x00,
            0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00,
            0x00, 0x80, 0x17, 0x00, 0x03, 0x07, 0x0a, 0x00, 0x01, 0x03, 0x00, 0x00, 0x02, 0x00,
            0x00, 0xac, 0x00, 0x00, 0x00, 0xac, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x07, 0x17,
            0x00, 0x04, 0x00, 0x09, 0x00, 0x03, 0x01, 0x00, 0x0c, 0x02, 0x00, 0x00, 0xa8, 0x00,
            0x00, 0x00, 0xa8, 0x00, 0x00, 0x00, 0xa8, 0x00, 0x00, 0x80, 0x15, 0x00, 0x03, 0x0c,
            0x00, 0x04, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00,
            0x00, 0x80, 0x00, 0x00,
        ];

        let result = decode(&src).unwrap();

        let expected = b"\
I17_08765:2:123:61541:01763#9\0\
I17_08765:2:123:1636:08611#9\0\
I17_08765:2:124:45613:16161#9\0\
";

        assert_eq!(result, expected);

        // The same block with every stream re-coded as an arith CAT stream
        // and `use_arith` set decodes to the same names.
        let arith = recode_streams_as_arith_cat(&src);
        assert_eq!(decode(&arith).unwrap(), expected);
    }

    #[allow(
        clippy::arithmetic_side_effects,
        clippy::cast_possible_truncation,
        reason = "g < 5, and the masked group fits a byte"
    )]
    fn put_uint7(out: &mut Vec<u8>, v: usize) {
        let v = u32::try_from(v).unwrap();
        let groups = (0..5u32).rev().skip_while(|&g| g > 0 && v >> (7 * g) == 0);
        for g in groups {
            let bits = ((v >> (7 * g)) & 0x7f) as u8;
            out.push(if g == 0 { bits } else { bits | 0x80 });
        }
    }

    /// Rewrite a rANS-coded tok3 block's streams as arith CAT streams
    /// (flags 0x20, the length, the bytes) and set `use_arith`.
    #[allow(clippy::indexing_slicing, reason = "test input is a known-valid block")]
    fn recode_streams_as_arith_cat(src: &[u8]) -> Vec<u8> {
        let (header, mut cur) = src.split_at(9);
        let mut out = header.to_vec();
        out[8] = 1;
        while let Some((&ttype, rest)) = cur.split_first() {
            out.push(ttype);
            cur = rest;
            if ttype & 0x40 != 0 {
                out.extend(&cur[..2]);
                cur = &cur[2..];
                continue;
            }
            let clen = read_uint7(&mut cur).unwrap() as usize;
            let data = super::super::rans_nx16::decode(&cur[..clen], 0).unwrap();
            cur = &cur[clen..];
            let mut stream = vec![0x20];
            put_uint7(&mut stream, data.len());
            stream.extend(data);
            put_uint7(&mut out, stream.len());
            out.extend(stream);
        }
        out
    }

    /// Arbitrary arith-coded blocks never panic, and tok3 over the
    /// production and the reference arithmetic decoder agree on whatever
    /// both accept.
    // r[verify cram.codec.tok3_arith]
    #[hegel::test]
    fn arbitrary_arith_blocks_never_panic(tc: hegel::TestCase) {
        use hegel::generators as gs;
        let mut src = tc.draw(gs::binary().min_size(9).max_size(1024));
        // A small name count, and `use_arith`.
        src.splice(4..9, [3, 0, 0, 0, 1]);
        if let (Ok(ours), Ok(theirs)) = (decode(&src), decode_with_arith_reference(&src)) {
            assert_eq!(ours, theirs);
        }
    }

    // r[verify cram.tok3.name_count_limit]
    #[test]
    fn tok3_name_count_exceeds_limit() {
        let mut src = Vec::new();
        src.extend_from_slice(&0u32.to_le_bytes()); // uncompressed_size
        src.extend_from_slice(&20_000_000u32.to_le_bytes()); // name_count > limit
        src.push(0u8); // method = 0 (rANS)

        let err = decode(&src).unwrap_err();
        assert!(
            matches!(err, CramError::Tok3NameCountExceedsLimit { count: 20_000_000, .. }),
            "expected Tok3NameCountExceedsLimit, got: {err:?}"
        );
    }
}
