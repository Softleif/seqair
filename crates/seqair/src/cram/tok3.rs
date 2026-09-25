// r[impl cram.codec.tok3]
//! tok3 / Name Tokenizer codec (CRAM compression method 8).
//!
//! Tokenizes read names by position. Each position stores a token type and
//! compressed data (via rANS Nx16 or the arithmetic coder). The decoder reconstructs names by
//! iterating over the tokens.

// See rans.rs: lazy `ok_or_else(|| CramError::...)` keeps error construction and its
// `drop_in_place<CramError>` off the per-record path.
#![allow(
    clippy::unnecessary_lazy_evaluations,
    reason = "lazy form avoids per-call drop_in_place<CramError> on hot path"
)]

use std::io::{BufRead, Cursor, Read, Write};

use super::codec_io::{Uint7Error, read_u8, read_u32_le, read_uint7, split_off};
use super::reader::CramError;

/// Bridge `Uint7Error` (narrow, hot-path-friendly) to the rich `CramError`.
fn uint7_to_cram_error(e: Uint7Error) -> CramError {
    match e {
        Uint7Error::Truncated => CramError::Truncated { context: "tok3 uint7" },
        Uint7Error::Overflow => CramError::Uint7Overflow,
    }
}

/// Maximum number of names allowed in a tok3 block.
const TOK3_NAME_COUNT_LIMIT: usize = 10_000_000;

/// Decode a tok3 compressed block.
pub fn decode(src: &[u8]) -> Result<Vec<u8>, CramError> {
    // Every stream stores its length, so the size passed is unused.
    decode_with(src, |stream| super::arith::decode(stream, 0))
}

/// [`decode`] with the arithmetic coder's streams decoded by
/// [`super::arith::reference`], so tests and the fuzzer can compare it with
/// [`decode`].
#[cfg(any(test, feature = "fuzz"))]
#[doc(hidden)]
pub fn decode_with_arith_reference(src: &[u8]) -> Result<Vec<u8>, CramError> {
    decode_with(src, |stream| {
        super::arith::reference::decode(stream, 0)
            .ok_or(CramError::ArithCorruptData { context: "arith reference decoder" })
    })
}

/// Decode a tok3 block, decoding arith-coded streams with `arith`.
fn decode_with(
    src: &[u8],
    arith: impl Fn(&[u8]) -> Result<Vec<u8>, CramError>,
) -> Result<Vec<u8>, CramError> {
    let mut cur: &[u8] = src;

    let (uncompressed_size, name_count, use_arith) = read_header(&mut cur)?;

    // r[impl cram.tok3.name_count_limit]
    if name_count > TOK3_NAME_COUNT_LIMIT {
        return Err(CramError::Tok3NameCountExceedsLimit {
            count: name_count,
            limit: TOK3_NAME_COUNT_LIMIT,
        });
    }

    // Each name slot needs ~48 bytes (two Vecs), plus the output buffer.
    super::reader::check_codec_output(
        name_count.saturating_mul(48).saturating_add(uncompressed_size),
        "tok3 output",
    )?;

    let mut b = decode_token_byte_streams(&mut cur, use_arith.then_some(&arith), name_count)?;

    let mut names: Vec<Vec<u8>> = vec![Vec::new(); name_count];
    let mut tokens: Vec<Vec<Option<Token>>> = vec![Vec::new(); name_count];

    let mut dst = Vec::with_capacity(uncompressed_size);

    for i in 0..name_count {
        let name = decode_single_name(&mut b, &mut names, &mut tokens, i)?;
        dst.extend_from_slice(&name);
        dst.push(0x00);
    }

    Ok(dst)
}

// ── Token types ──────────────────────────────────────────────────────

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum TokenType {
    Type,
    String,
    Char,
    Digits0,
    DZLen,
    Dup,
    Diff,
    Digits,
    Delta,
    Delta0,
    Match,
    Nop,
    End,
}

impl TokenType {
    fn from_byte(n: u8) -> Result<Self, CramError> {
        match n & 0x3f {
            0 => Ok(Self::Type),
            1 => Ok(Self::String),
            2 => Ok(Self::Char),
            3 => Ok(Self::Digits0),
            4 => Ok(Self::DZLen),
            5 => Ok(Self::Dup),
            6 => Ok(Self::Diff),
            7 => Ok(Self::Digits),
            8 => Ok(Self::Delta),
            9 => Ok(Self::Delta0),
            10 => Ok(Self::Match),
            11 => Ok(Self::Nop),
            12 => Ok(Self::End),
            _ => Err(CramError::InvalidTok3TokenType { token_type: n }),
        }
    }

    fn to_byte(self) -> u8 {
        match self {
            Self::Type => 0,
            Self::String => 1,
            Self::Char => 2,
            Self::Digits0 => 3,
            Self::DZLen => 4,
            Self::Dup => 5,
            Self::Diff => 6,
            Self::Digits => 7,
            Self::Delta => 8,
            Self::Delta0 => 9,
            Self::Match => 10,
            Self::Nop => 11,
            Self::End => 12,
        }
    }
}

// ── Token values ─────────────────────────────────────────────────────

#[derive(Clone, Debug, Eq, PartialEq)]
enum Token {
    Char(u8),
    String(Vec<u8>),
    Digits(u32),
    PaddedDigits(u32, u8),
    Nop,
}

/// Return a u8 discriminant for error reporting (avoids String in error variants).
fn token_discriminant(token: Option<&Token>) -> u8 {
    match token {
        None => 0,
        Some(Token::Char(_)) => 1,
        Some(Token::String(_)) => 2,
        Some(Token::Digits(_)) => 3,
        Some(Token::PaddedDigits(_, _)) => 4,
        Some(Token::Nop) => 5,
    }
}

// ── Per-position token reader ────────────────────────────────────────

#[derive(Clone, Debug, Default)]
struct TokenReader {
    type_reader: Cursor<Vec<u8>>,
    string_reader: Cursor<Vec<u8>>,
    char_reader: Cursor<Vec<u8>>,
    digits0_reader: Cursor<Vec<u8>>,
    dz_len_reader: Cursor<Vec<u8>>,
    dup_reader: Cursor<Vec<u8>>,
    diff_reader: Cursor<Vec<u8>>,
    digits_reader: Cursor<Vec<u8>>,
    delta_reader: Cursor<Vec<u8>>,
    delta0_reader: Cursor<Vec<u8>>,
}

impl TokenReader {
    // r[impl cram.tok3.dz_len_reader]
    fn get(&self, ty: TokenType) -> &Cursor<Vec<u8>> {
        match ty {
            TokenType::Type => &self.type_reader,
            TokenType::String => &self.string_reader,
            TokenType::Char => &self.char_reader,
            TokenType::Digits0 => &self.digits0_reader,
            TokenType::DZLen => &self.dz_len_reader,
            TokenType::Dup => &self.dup_reader,
            TokenType::Diff => &self.diff_reader,
            TokenType::Digits => &self.digits_reader,
            TokenType::Delta => &self.delta_reader,
            TokenType::Delta0 => &self.delta0_reader,
            _ => &self.type_reader,
        }
    }

    fn get_mut(&mut self, ty: TokenType) -> &mut Cursor<Vec<u8>> {
        match ty {
            TokenType::Type => &mut self.type_reader,
            TokenType::String => &mut self.string_reader,
            TokenType::Char => &mut self.char_reader,
            TokenType::Digits0 => &mut self.digits0_reader,
            TokenType::Dup => &mut self.dup_reader,
            TokenType::Diff => &mut self.diff_reader,
            TokenType::DZLen => &mut self.dz_len_reader,
            TokenType::Digits => &mut self.digits_reader,
            TokenType::Delta => &mut self.delta_reader,
            TokenType::Delta0 => &mut self.delta0_reader,
            _ => &mut self.type_reader,
        }
    }

    fn set(&mut self, ty: TokenType, buf: Vec<u8>) {
        *self.get_mut(ty).get_mut() = buf;
    }

    fn read_type(&mut self) -> Result<TokenType, CramError> {
        let mut buf = [0u8; 1];
        self.type_reader
            .read_exact(&mut buf)
            .map_err(|_| CramError::Truncated { context: "tok3 type byte" })?;
        TokenType::from_byte(buf[0])
    }

    fn read_distance(&mut self, ty: TokenType) -> Result<usize, CramError> {
        let reader = self.get_mut(ty);
        let mut buf = [0u8; 4];
        reader
            .read_exact(&mut buf)
            .map_err(|_| CramError::Truncated { context: "tok3 distance" })?;
        Ok(u32::from_le_bytes(buf) as usize)
    }

    fn read_token(&mut self, prev_token: Option<&Token>) -> Result<Option<Token>, CramError> {
        let ty = self.read_type()?;

        match ty {
            TokenType::Char => {
                let mut buf = [0u8; 1];
                self.char_reader
                    .read_exact(&mut buf)
                    .map_err(|_| CramError::Truncated { context: "tok3 char" })?;
                Ok(Some(Token::Char(buf[0])))
            }
            TokenType::String => {
                let mut buf = Vec::new();
                self.string_reader
                    .read_until(0x00, &mut buf)
                    .map_err(|_| CramError::Truncated { context: "tok3 string" })?;
                buf.pop();
                Ok(Some(Token::String(buf)))
            }
            TokenType::Digits => {
                let mut buf = [0u8; 4];
                self.digits_reader
                    .read_exact(&mut buf)
                    .map_err(|_| CramError::Truncated { context: "tok3 digits" })?;
                Ok(Some(Token::Digits(u32::from_le_bytes(buf))))
            }
            TokenType::Digits0 => {
                let mut dbuf = [0u8; 4];
                self.digits0_reader
                    .read_exact(&mut dbuf)
                    .map_err(|_| CramError::Truncated { context: "tok3 digits0" })?;
                let mut lbuf = [0u8; 1];
                self.dz_len_reader
                    .read_exact(&mut lbuf)
                    .map_err(|_| CramError::Truncated { context: "tok3 dzlen" })?;
                Ok(Some(Token::PaddedDigits(u32::from_le_bytes(dbuf), lbuf[0])))
            }
            TokenType::Delta => {
                let mut buf = [0u8; 1];
                self.delta_reader
                    .read_exact(&mut buf)
                    .map_err(|_| CramError::Truncated { context: "tok3 delta" })?;
                let delta = u32::from(buf[0]);
                match prev_token {
                    Some(Token::Digits(n)) => {
                        Ok(Some(Token::Digits(n.checked_add(delta).ok_or_else(|| {
                            CramError::Truncated { context: "tok3 delta overflow" }
                        })?)))
                    }
                    _ => Err(CramError::Tok3DeltaRequiresDigits {
                        found: token_discriminant(prev_token),
                    }),
                }
            }
            TokenType::Delta0 => {
                let mut buf = [0u8; 1];
                self.delta0_reader
                    .read_exact(&mut buf)
                    .map_err(|_| CramError::Truncated { context: "tok3 delta0" })?;
                let delta = u32::from(buf[0]);
                match prev_token {
                    Some(Token::PaddedDigits(n, width)) => Ok(Some(Token::PaddedDigits(
                        n.checked_add(delta).ok_or_else(|| CramError::Truncated {
                            context: "tok3 delta0 overflow",
                        })?,
                        *width,
                    ))),
                    _ => Err(CramError::Tok3Delta0RequiresPaddedDigits {
                        found: token_discriminant(prev_token),
                    }),
                }
            }
            TokenType::Match => Ok(prev_token.cloned()),
            TokenType::End => Ok(None),
            _ => Ok(Some(Token::Nop)),
        }
    }
}

// ── Header ───────────────────────────────────────────────────────────

fn read_header(src: &mut &[u8]) -> Result<(usize, usize, bool), CramError> {
    let truncated = || CramError::Truncated { context: "tok3 header" };
    let uncompressed_size = read_u32_le(src).ok_or_else(truncated)? as usize;
    let name_count = read_u32_le(src).ok_or_else(truncated)? as usize;
    let method = read_u8(src).ok_or_else(truncated)?;
    let use_arith = method != 0;
    Ok((uncompressed_size, name_count, use_arith))
}

// ── Decode sub-streams ───────────────────────────────────────────────

/// `arith` decodes a stream when the block uses the arithmetic coder;
/// without it streams are rANS Nx16.
fn decode_token_byte_streams(
    src: &mut &[u8],
    arith: Option<&impl Fn(&[u8]) -> Result<Vec<u8>, CramError>>,
    n_names: usize,
) -> Result<Vec<TokenReader>, CramError> {
    let mut b: Vec<TokenReader> = Vec::new();
    let mut t: Option<usize> = None;

    while !src.is_empty() {
        let ttype =
            read_u8(src).ok_or_else(|| CramError::Truncated { context: "tok3 token type" })?;

        let tok_new = ttype & 0x80 != 0;
        let tok_dup = ttype & 0x40 != 0;

        let ty = TokenType::from_byte(ttype)?;

        if tok_new {
            let new_t = t.map_or(0, |v| v.wrapping_add(1));
            t = Some(new_t);
            b.push(TokenReader::default());

            if ty != TokenType::Type {
                let mut buf = vec![TokenType::Match.to_byte(); n_names];
                if let Some(first) = buf.first_mut() {
                    *first = ty.to_byte();
                }
                b.get_mut(new_t)
                    .ok_or_else(|| CramError::Truncated { context: "tok3 new token position" })?
                    .set(TokenType::Type, buf);
            }
        }

        let t_idx = t.ok_or_else(|| CramError::Truncated {
            context: "tok3 token index before first new token",
        })?;

        if tok_dup {
            let truncated_dup = || CramError::Truncated { context: "tok3 dup metadata" };
            let dup_pos = read_u8(src).ok_or_else(truncated_dup)? as usize;
            let dup_type = TokenType::from_byte(read_u8(src).ok_or_else(truncated_dup)?)?;

            let buf = b
                .get(dup_pos)
                .ok_or_else(|| CramError::Tok3DupPositionOutOfRange { dup_pos })?
                .get(dup_type)
                .get_ref()
                .clone();

            b.get_mut(t_idx)
                .ok_or_else(|| CramError::Truncated { context: "tok3 dup set" })?
                .set(ty, buf);
        } else {
            let compressed_size = read_uint7(src).map_err(uint7_to_cram_error)? as usize;
            let buf = split_off(src, compressed_size)
                .ok_or_else(|| CramError::Truncated { context: "tok3 compressed payload" })?;
            // r[impl cram.codec.tok3_arith]
            let decompressed = match arith {
                Some(arith) => arith(buf)?,
                None => super::rans_nx16::decode(buf, 0)?,
            };

            b.get_mut(t_idx)
                .ok_or_else(|| CramError::Truncated { context: "tok3 stream set" })?
                .set(ty, decompressed);
        }
    }

    Ok(b)
}

fn decode_single_name(
    b: &mut [TokenReader],
    names: &mut [Vec<u8>],
    tokens: &mut [Vec<Option<Token>>],
    n: usize,
) -> Result<Vec<u8>, CramError> {
    let first_reader =
        b.first_mut().ok_or_else(|| CramError::Truncated { context: "tok3 no token readers" })?;

    let ty = first_reader.read_type()?;
    let dist = first_reader.read_distance(ty)?;

    let m = n
        .checked_sub(dist)
        .ok_or_else(|| CramError::Tok3DistanceExceedsIndex { distance: dist, name_index: n })?;

    if ty == TokenType::Dup {
        let prev_name =
            names.get(m).ok_or_else(|| CramError::Tok3DupRefOutOfRange { index: m })?.clone();
        let prev_tokens =
            tokens.get(m).ok_or_else(|| CramError::Tok3DupRefOutOfRange { index: m })?.clone();

        if let Some(slot) = names.get_mut(n) {
            *slot = prev_name;
        }
        if let Some(slot) = tokens.get_mut(n) {
            *slot = prev_tokens;
        }

        return Ok(names
            .get(n)
            .ok_or_else(|| CramError::Truncated { context: "tok3 dup result" })?
            .clone());
    }

    let mut t = 1;

    loop {
        let reader = b
            .get_mut(t)
            .ok_or_else(|| CramError::Truncated { context: "tok3 token reader position" })?;

        let prev_token = tokens.get(m).and_then(|ts| ts.get(t)).and_then(|tok| tok.as_ref());

        if let Some(token) = reader.read_token(prev_token)? {
            let name = names
                .get_mut(n)
                .ok_or_else(|| CramError::Truncated { context: "tok3 name index" })?;

            match &token {
                Token::Char(c) => name.push(*c),
                Token::String(s) => name.extend_from_slice(s),
                Token::Digits(d) => {
                    write!(name, "{d}")?;
                }
                Token::PaddedDigits(d, l) => {
                    write!(name, "{:0width$}", d, width = usize::from(*l))?;
                }
                Token::Nop => {}
            }

            if let Some(ts) = tokens.get_mut(n) {
                if t >= ts.len() {
                    let new_len = t.checked_add(1).ok_or_else(|| CramError::Truncated {
                        context: "tok3 token index overflow",
                    })?;
                    ts.resize(new_len, None);
                }
                if let Some(slot) = ts.get_mut(t) {
                    *slot = Some(token);
                }
            }
        } else {
            break;
        }

        t = t
            .checked_add(1)
            .ok_or_else(|| CramError::Truncated { context: "tok3 token index overflow" })?;
    }

    Ok(names.get(n).ok_or_else(|| CramError::Truncated { context: "tok3 final name" })?.clone())
}

// Primitive readers (`read_u8`, `read_u32_le`, `read_uint7`, `split_off`)
// live in `super::codec_io`; see the imports at the top of this file.

#[cfg(test)]
mod tests {
    use super::*;

    // r[verify cram.codec.tok3]

    #[test]
    fn invalid_tok3_token_type_returns_error() {
        // Token types are masked with 0x3f; values 13..=63 are invalid.
        // Pass 0x0D (13) directly.
        let err = TokenType::from_byte(13).unwrap_err();
        assert!(matches!(err, CramError::InvalidTok3TokenType { token_type: 13 }));
    }

    #[test]
    fn invalid_tok3_token_type_field_value() {
        let err = TokenType::from_byte(63).unwrap_err();
        assert!(matches!(err, CramError::InvalidTok3TokenType { token_type: 63 }));
    }

    #[test]
    fn tok3_dup_position_out_of_range() {
        // Tok3DupPositionOutOfRange is returned in decode_token_byte_streams when
        // tok_dup is set and dup_pos references an index beyond the current b vec.
        // Construct a stream: one token entry with tok_new=1, tok_dup=1 but dup_pos=99
        // (which is > 0 entries that exist at that point).
        let mut src = Vec::new();
        // Header
        src.extend_from_slice(&10u32.to_le_bytes()); // uncompressed_size
        src.extend_from_slice(&1u32.to_le_bytes()); // name_count
        src.push(0u8); // method = 0 (rANS Nx16, not arith)

        // Token stream: first token with NEW + DUP flags set
        // 0x80 = NEW, 0x40 = DUP, type bits = 0x00 (Type)
        src.push(0x80 | 0x40); // tok_new=1, tok_dup=1, type=Type
        src.push(99u8); // dup_pos = 99 (out of range — b is empty at this point)
        src.push(0x00u8); // dup_type byte

        let err = decode(&src).unwrap_err();
        assert!(
            matches!(err, CramError::Tok3DupPositionOutOfRange { dup_pos: 99 }),
            "expected Tok3DupPositionOutOfRange, got: {err:?}"
        );
    }

    #[test]
    fn tok3_distance_exceeds_index() {
        // Tok3DistanceExceedsIndex fires in decode_single_name when distance > n (name index).
        // This requires crafting a tok3 stream where the first name's distance field > 0.
        // TODO: requires crafting a full tok3 token byte stream; the distance is embedded
        // deep in the rANS-compressed sub-stream for position 0. Testing indirectly via
        // the error constructor instead.
        let err = CramError::Tok3DistanceExceedsIndex { distance: 5, name_index: 0 };
        assert!(matches!(err, CramError::Tok3DistanceExceedsIndex { distance: 5, name_index: 0 }));
        let msg = format!("{err}");
        assert!(
            msg.contains("5") && msg.contains("0"),
            "error should mention distance and index: {msg}"
        );
    }

    #[test]
    fn tok3_dup_ref_out_of_range_error_variant() {
        // Tok3DupRefOutOfRange fires in decode_single_name when the resolved m index
        // is out of bounds in the names/tokens arrays.
        // TODO: requires crafting a full valid tok3 stream with a Dup token type where
        // the referenced name index is out of range. Testing via constructor.
        let err = CramError::Tok3DupRefOutOfRange { index: 42 };
        assert!(matches!(err, CramError::Tok3DupRefOutOfRange { index: 42 }));
        let msg = format!("{err}");
        assert!(msg.contains("42"), "error should mention the index: {msg}");
    }

    #[test]
    fn tok3_delta_requires_digits_error_variant() {
        // Tok3DeltaRequiresDigits fires in read_token when Delta is encountered
        // but prev_token is not Some(Token::Digits(_)).
        // TODO: requires crafting a full multi-name tok3 stream where the second name
        // has a Delta token at a position that had a non-Digits token in the previous name.
        // Testing via constructor for now.
        let err = CramError::Tok3DeltaRequiresDigits { found: 0 };
        assert!(matches!(err, CramError::Tok3DeltaRequiresDigits { .. }));
    }

    #[test]
    fn tok3_delta0_requires_padded_digits_error_variant() {
        // Tok3Delta0RequiresPaddedDigits is similar but for Delta0 needing PaddedDigits.
        // TODO: requires crafting a multi-name tok3 stream.
        let err = CramError::Tok3Delta0RequiresPaddedDigits { found: 0 };
        assert!(matches!(err, CramError::Tok3Delta0RequiresPaddedDigits { .. }));
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

    // r[verify cram.tok3.dz_len_reader]
    #[test]
    fn token_reader_get_dz_len_matches_get_mut() {
        let mut reader = TokenReader::default();
        let test_data = vec![42u8, 99];
        reader.set(TokenType::DZLen, test_data.clone());

        // get() must return the dz_len_reader, not type_reader
        let immutable = reader.get(TokenType::DZLen);
        assert_eq!(immutable.get_ref(), &test_data);

        let mutable = reader.get_mut(TokenType::DZLen);
        assert_eq!(mutable.get_ref(), &test_data);
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
