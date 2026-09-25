// r[impl cram.codec.arith+2]
//! The CRAM 3.1 adaptive arithmetic coder (compression method 6, "arith" or
//! "range").
//!
//! An adaptive order-0 or order-1 model over the byte-wise range coder of
//! [`super::range_coder`], optionally with run lengths coded alongside, and
//! wrapped in the same transforms as rANS Nx16: STRIPE, PACK, CAT, plus EXT
//! (bzip2). Follows htscodecs' `arith_dynamic.c` where the spec's pseudocode
//! differs; `docs/spec/2-cram-1-reader.md` names each place.

use std::io::Read;

#[cfg(test)]
mod htscodecs_oracle;
#[cfg(any(test, feature = "fuzz"))]
#[doc(hidden)]
pub mod reference;

use super::codec_io::{self, Uint7Error, read_u8};
use super::range_coder::{AdaptiveModel, RangeDecoder};
use super::rans_nx16::{apply_bit_unpack, read_bit_pack_context};
use super::reader::{CramError, check_codec_output};

/// The ORDER field: order-1 when it is exactly 1 (htscodecs decodes the
/// reserved 2 and 3 as order-0).
const ORDER_MASK: u8 = 0x03;
const FLAG_EXT: u8 = 0x04;
const FLAG_STRIPE: u8 = 0x08;
const FLAG_NO_SIZE: u8 = 0x10;
const FLAG_CAT: u8 = 0x20;
const FLAG_RLE: u8 = 0x40;
const FLAG_PACK: u8 = 0x80;

/// How deep STRIPE streams may nest. htscodecs' encoder writes one level
/// (its substreams never stripe); the bound keeps a crafted stream from
/// recursing until the stack runs out.
const MAX_STRIPE_DEPTH: u8 = 4;

/// Symbols of a run-length model: parts of 0..=3, where 3 continues the run.
const RUN_SYMBOLS: usize = 4;
const RUN_CONTINUES: u16 = 3;
/// Run-length contexts: the literal (0..=255) for a run's first part, then
/// 256 for its second and 257 for every later one.
const RUN_CONTEXTS: usize = 258;
const RUN_CTX_SECOND: usize = 256;
const RUN_CTX_LATER: usize = 257;

/// A literal model covers every byte value.
type LiteralModel = AdaptiveModel<256>;
type RunModel = AdaptiveModel<RUN_SYMBOLS>;

/// Decode an arithmetic-coder stream.
///
/// `uncompressed_size` is used only by a stream that does not store its
/// length (the NOSZ flag); otherwise the stored length wins, and the result
/// always has the length the stream declares.
pub fn decode(src: &[u8], uncompressed_size: usize) -> Result<Vec<u8>, CramError> {
    decode_nested(src, uncompressed_size, 0)
}

fn read_uint7(src: &mut &[u8]) -> Result<usize, CramError> {
    match codec_io::read_uint7(src) {
        Ok(n) => Ok(n as usize),
        Err(Uint7Error::Truncated) => Err(CramError::Truncated { context: "arith uint7" }),
        Err(Uint7Error::Overflow) => Err(CramError::Uint7Overflow),
    }
}

// r[impl cram.codec.arith.wrapper]
fn decode_nested(src: &[u8], uncompressed_size: usize, depth: u8) -> Result<Vec<u8>, CramError> {
    let mut cur = src;
    let flags = read_u8(&mut cur).ok_or(CramError::Truncated { context: "arith flags" })?;
    if flags & FLAG_STRIPE != 0 {
        return decode_stripe(cur, depth);
    }

    let len = if flags & FLAG_NO_SIZE != 0 { uncompressed_size } else { read_uint7(&mut cur)? };
    // r[impl io.fuzz.alloc_limits]
    check_codec_output(len, "arith output")?;

    let pack = if flags & FLAG_PACK != 0 {
        let (ctx, packed_len) = read_bit_pack_context(&mut cur, len)?;
        if packed_len > len {
            return Err(CramError::ArithPackedLengthTooLarge { packed: packed_len, unpacked: len });
        }
        Some((ctx, packed_len))
    } else {
        None
    };

    let body_len = pack.as_ref().map_or(len, |&(_, packed_len)| packed_len);
    let body = decode_body(flags, cur, body_len)?;

    let out = match pack {
        // r[impl cram.codec.arith.pack]
        Some((ctx, _)) => {
            // `apply_bit_unpack` leaves what the packed bytes do not cover
            // zeroed; htscodecs' `hts_unpack` rejects such a stream.
            let needed = match ctx.symbol_count {
                1 => 0,
                2 => len.div_ceil(8),
                3..=4 => len.div_ceil(4),
                _ => len.div_ceil(2),
            };
            if body.len() < needed {
                return Err(CramError::ArithLengthMismatch {
                    expected: needed,
                    actual: body.len(),
                });
            }
            apply_bit_unpack(&body, &ctx)?
        }
        None => body,
    };

    if out.len() != len {
        return Err(CramError::ArithLengthMismatch { expected: len, actual: out.len() });
    }
    Ok(out)
}

/// Decode the body after the metadata to (at most) `len` bytes.
fn decode_body(flags: u8, body: &[u8], len: usize) -> Result<Vec<u8>, CramError> {
    // htscodecs decodes an empty body to nothing, whatever the flags say;
    // the length check after it decides whether that was valid.
    if body.is_empty() {
        return Ok(Vec::new());
    }
    if flags & FLAG_CAT != 0 {
        let data = body.get(..len).ok_or(CramError::Truncated { context: "arith CAT data" })?;
        return Ok(data.to_vec());
    }
    if flags & FLAG_EXT != 0 {
        return decode_ext(body, len);
    }
    if len == 0 {
        return Ok(Vec::new());
    }

    let order1 = flags & ORDER_MASK == 1;
    let rle = flags & FLAG_RLE != 0;
    match (rle, order1) {
        (false, false) => decode_order::<false>(body, len),
        (false, true) => decode_order::<true>(body, len),
        (true, false) => decode_rle::<false>(body, len),
        (true, true) => decode_rle::<true>(body, len),
    }
}

/// Read the `max_sym` byte (0 meaning 256) and start the range coder on the
/// rest.
fn start(body: &[u8]) -> Result<(usize, RangeDecoder<'_>), CramError> {
    let (&m, rest) = body.split_first().ok_or(CramError::Truncated { context: "arith max_sym" })?;
    let max_sym = if m == 0 { 256 } else { usize::from(m) };
    Ok((max_sym, RangeDecoder::new(rest)?))
}

/// The literal models: one for order-0, one per previous byte for order-1.
/// Every decoded byte is below `max_sym`, so contexts past it never occur.
fn literal_models<const ORDER1: bool>(max_sym: usize) -> Vec<LiteralModel> {
    let contexts = if ORDER1 { max_sym } else { 1 };
    vec![LiteralModel::new(max_sym); contexts]
}

fn corrupt(context: &'static str) -> impl FnOnce() -> CramError {
    move || CramError::ArithCorruptData { context }
}

/// A literal model's symbol as a byte; a `LiteralModel` has at most 256.
#[allow(clippy::cast_possible_truncation, reason = "LiteralModel symbols are < 256")]
fn byte(sym: u16) -> u8 {
    sym as u8
}

// r[impl cram.codec.arith.order]
fn decode_order<const ORDER1: bool>(body: &[u8], len: usize) -> Result<Vec<u8>, CramError> {
    let (max_sym, mut rc) = start(body)?;
    let mut models = literal_models::<ORDER1>(max_sym);
    let mut out = vec![0u8; len];
    let mut ctx = 0usize;
    for d in &mut out {
        let model = models.get_mut(ctx).ok_or_else(corrupt("arith literal context"))?;
        let sym = model.decode(&mut rc).ok_or_else(corrupt("arith literal"))?;
        *d = byte(sym);
        if ORDER1 {
            ctx = usize::from(sym);
        }
    }
    Ok(out)
}

// r[impl cram.codec.arith.rle]
fn decode_rle<const ORDER1: bool>(body: &[u8], len: usize) -> Result<Vec<u8>, CramError> {
    let (max_sym, mut rc) = start(body)?;
    let mut literals = literal_models::<ORDER1>(max_sym);
    let mut runs = vec![RunModel::new(RUN_SYMBOLS); RUN_CONTEXTS];
    let mut out = Vec::with_capacity(len);
    let mut ctx = 0usize;
    while out.len() < len {
        let model = literals.get_mut(ctx).ok_or_else(corrupt("arith literal context"))?;
        let sym = model.decode(&mut rc).ok_or_else(corrupt("arith literal"))?;
        let lit = byte(sym);
        out.push(lit);
        if ORDER1 {
            ctx = usize::from(sym);
        }

        let mut run = 0usize;
        let mut rctx = usize::from(sym);
        loop {
            let model = runs.get_mut(rctx).ok_or_else(corrupt("arith run context"))?;
            let part = model.decode(&mut rc).ok_or_else(corrupt("arith run length"))?;
            rctx = if rctx == usize::from(sym) { RUN_CTX_SECOND } else { RUN_CTX_LATER };
            run = run.saturating_add(usize::from(part));
            // htscodecs stops at the output length even mid-run.
            if part != RUN_CONTINUES || run >= len {
                break;
            }
        }
        let room = len.saturating_sub(out.len());
        out.resize(out.len().saturating_add(run.min(room)), lit);
    }
    Ok(out)
}

// r[impl cram.codec.arith.stripe]
fn decode_stripe(mut cur: &[u8], depth: u8) -> Result<Vec<u8>, CramError> {
    if depth >= MAX_STRIPE_DEPTH {
        return Err(CramError::ArithStripeTooDeep { limit: MAX_STRIPE_DEPTH });
    }
    // The length is stored whatever NOSZ says.
    let len = read_uint7(&mut cur)?;
    // r[impl io.fuzz.alloc_limits]
    check_codec_output(len, "arith stripe output")?;
    let n = usize::from(
        read_u8(&mut cur).ok_or(CramError::Truncated { context: "arith stripe count" })?,
    );
    let q = len.checked_div(n).ok_or(CramError::ArithStripeZeroStreams)?;
    let r = len.checked_rem(n).ok_or(CramError::ArithStripeZeroStreams)?;
    let clens = (0..n).map(|_| read_uint7(&mut cur)).collect::<Result<Vec<_>, _>>()?;

    let mut out = vec![0u8; len];
    for (j, &clen) in clens.iter().enumerate() {
        let sub = codec_io::split_off(&mut cur, clen)
            .ok_or(CramError::Truncated { context: "arith stripe substream" })?;
        let part_len = q.saturating_add(usize::from(j < r));
        let part = decode_nested(sub, part_len, depth.saturating_add(1))?;
        if part.len() != part_len {
            return Err(CramError::ArithLengthMismatch { expected: part_len, actual: part.len() });
        }
        for (d, &s) in out.iter_mut().skip(j).step_by(n).zip(&part) {
            *d = s;
        }
    }
    Ok(out)
}

// r[impl cram.codec.arith.ext]
fn decode_ext(body: &[u8], len: usize) -> Result<Vec<u8>, CramError> {
    if !body.starts_with(b"BZh") {
        return Err(CramError::ArithExtUnknownCodec);
    }
    let mut out = Vec::with_capacity(len);
    // One byte past `len` is enough to tell that the output does not fit.
    let limit = u64::try_from(len).unwrap_or(u64::MAX).saturating_add(1);
    bzip2::read::BzDecoder::new(body)
        .take(limit)
        .read_to_end(&mut out)
        .map_err(|source| CramError::Bzip2DecompressionFailed { source })?;
    if out.len() > len {
        return Err(CramError::ArithLengthMismatch { expected: len, actual: out.len() });
    }
    Ok(out)
}

#[cfg(test)]
#[allow(
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    clippy::indexing_slicing,
    reason = "test-only encoder and generators with bounded values"
)]
mod tests {
    //! Round trips through a test-only encoder: the spec's encoder side
    //! (`RangeEncode`, the RLE run split, `hts_pack`'s bit layout, the
    //! stripe transpose), written independently of the decoder.

    use super::super::range_coder::test_encoder::RangeEncoder;
    use super::*;
    use hegel::prelude::*;

    fn put_uint7(out: &mut Vec<u8>, v: usize) {
        let v = u32::try_from(v).unwrap();
        let groups = (0..5u32).rev().skip_while(|&g| g > 0 && v >> (7 * g) == 0);
        for g in groups {
            let bits = ((v >> (7 * g)) & 0x7f) as u8;
            out.push(if g == 0 { bits } else { bits | 0x80 });
        }
    }

    /// htscodecs writes `max + 1`, as a byte (256 wraps to 0).
    fn max_sym(data: &[u8]) -> u16 {
        data.iter().max().map_or(1, |&m| u16::from(m) + 1)
    }

    fn encode_order(data: &[u8], order1: bool) -> Vec<u8> {
        let m = max_sym(data);
        let mut rc = RangeEncoder::new();
        let mut models = vec![LiteralModel::new(usize::from(m)); 256];
        let mut last = 0usize;
        for &b in data {
            models[if order1 { last } else { 0 }].encode(&mut rc, u16::from(b));
            last = usize::from(b);
        }
        let mut out = vec![m as u8];
        out.extend(rc.finish());
        out
    }

    /// htscodecs' `arith_compress_O{0,1}_RLE`: a literal, then its extra
    /// copies in parts of at most 3 — context literal, then 256, then 257 —
    /// with a closing 0 after a final part of 3.
    fn encode_rle(data: &[u8], order1: bool) -> Vec<u8> {
        let m = max_sym(data);
        let mut rc = RangeEncoder::new();
        let mut lits = vec![LiteralModel::new(usize::from(m)); 256];
        let mut runs = vec![RunModel::new(4); 258];
        let mut last = 0usize;
        let mut i = 0;
        while i < data.len() {
            let lit = data[i];
            lits[if order1 { last } else { 0 }].encode(&mut rc, u16::from(lit));
            last = usize::from(lit);
            i += 1;
            let mut run = 0usize;
            while i < data.len() && data[i] == lit {
                run += 1;
                i += 1;
            }
            let mut rctx = usize::from(lit);
            loop {
                let part = run.min(3);
                runs[rctx].encode(&mut rc, part as u16);
                run -= part;
                rctx = if rctx == usize::from(lit) { 256 } else { 257 };
                if part == 3 && run == 0 {
                    runs[rctx].encode(&mut rc, 0);
                }
                if run == 0 {
                    break;
                }
            }
        }
        let mut out = vec![m as u8];
        out.extend(rc.finish());
        out
    }

    /// `hts_pack`: the sorted symbols as the map, then 8, 4 or 2 indices a
    /// byte, first in the lowest bits; one symbol packs to nothing.
    fn hts_pack(data: &[u8]) -> Option<(Vec<u8>, Vec<u8>)> {
        let mut seen = [false; 256];
        for &b in data {
            seen[usize::from(b)] = true;
        }
        let syms: Vec<u8> = (0..=255u8).filter(|&s| seen[usize::from(s)]).collect();
        if data.is_empty() || syms.len() > 16 {
            return None;
        }
        let mut meta = vec![syms.len() as u8];
        meta.extend(&syms);
        let index = |b: u8| syms.iter().position(|&s| s == b).unwrap() as u8;
        let per_byte = match syms.len() {
            1 => return Some((meta, Vec::new())),
            2 => 8,
            3..=4 => 4,
            _ => 2,
        };
        let bits = 8 / per_byte;
        let packed = data
            .chunks(per_byte)
            .map(|c| c.iter().enumerate().fold(0u8, |acc, (k, &b)| acc | (index(b) << (k * bits))))
            .collect();
        Some((meta, packed))
    }

    fn bzip2(data: &[u8]) -> Vec<u8> {
        let mut out = Vec::new();
        bzip2::read::BzEncoder::new(data, bzip2::Compression::fast())
            .read_to_end(&mut out)
            .unwrap();
        out
    }

    #[derive(Debug, Clone, Copy)]
    enum Coder {
        /// `order` is the whole ORDER field, reserved values included.
        Range {
            order: u8,
            rle: bool,
        },
        Cat,
        Ext,
    }

    #[derive(Debug, Clone)]
    enum Layout {
        Plain {
            coder: Coder,
            pack: bool,
            no_size: bool,
        },
        /// `junk` are other flag bits, which a STRIPE stream ignores.
        Stripe {
            junk: u8,
            parts: Vec<Layout>,
        },
    }

    fn encode(layout: &Layout, data: &[u8]) -> Vec<u8> {
        match layout {
            Layout::Plain { coder, pack, no_size } => {
                let mut flags = 0u8;
                let mut out = vec![0u8];
                if *no_size {
                    flags |= FLAG_NO_SIZE;
                } else {
                    put_uint7(&mut out, data.len());
                }
                let packed;
                let body: &[u8] = match hts_pack(data).filter(|_| *pack) {
                    Some((meta, p)) => {
                        flags |= FLAG_PACK;
                        out.extend(meta);
                        put_uint7(&mut out, p.len());
                        packed = p;
                        &packed
                    }
                    None => data,
                };
                match *coder {
                    Coder::Cat => {
                        flags |= FLAG_CAT;
                        out.extend(body);
                    }
                    Coder::Ext => {
                        flags |= FLAG_EXT;
                        out.extend(bzip2(body));
                    }
                    Coder::Range { order, rle } => {
                        flags |= order;
                        if rle {
                            flags |= FLAG_RLE;
                            out.extend(encode_rle(body, order == 1));
                        } else {
                            out.extend(encode_order(body, order == 1));
                        }
                    }
                }
                out[0] = flags;
                out
            }
            Layout::Stripe { junk, parts } => {
                let n = parts.len();
                let mut out = vec![FLAG_STRIPE | junk];
                put_uint7(&mut out, data.len());
                out.push(n as u8);
                let subs: Vec<Vec<u8>> = parts
                    .iter()
                    .enumerate()
                    .map(|(j, part)| {
                        let column: Vec<u8> = data.iter().skip(j).step_by(n).copied().collect();
                        encode(part, &column)
                    })
                    .collect();
                for s in &subs {
                    put_uint7(&mut out, s.len());
                }
                for s in subs {
                    out.extend(s);
                }
                out
            }
        }
    }

    /// Data in the shapes the coder meets: raw bytes, a small alphabet in
    /// runs, long runs of a couple of symbols, and little-endian u32s.
    #[hegel::composite]
    fn arb_data(tc: &TestCase) -> Vec<u8> {
        match tc.draw_silent(gs::integers::<u8>().max_value(3)) {
            0 => tc.draw_silent(gs::binary().max_size(2000)),
            1 => {
                let k = tc.draw_silent(gs::integers::<u8>().min_value(1).max_value(24));
                let base = tc.draw_silent(gs::integers::<u8>());
                let runs = tc.draw_silent(
                    gs::vecs(gs::tuples!(
                        gs::integers::<u8>().max_value(k - 1),
                        gs::integers::<usize>().max_value(12)
                    ))
                    .max_size(300),
                );
                runs.into_iter()
                    .flat_map(|(s, n)| std::iter::repeat_n(base.wrapping_add(s), n + 1))
                    .collect()
            }
            2 => {
                let runs = tc.draw_silent(
                    gs::vecs(gs::tuples!(
                        gs::integers::<u8>().max_value(1),
                        gs::integers::<usize>().max_value(3000)
                    ))
                    .max_size(8),
                );
                runs.into_iter().flat_map(|(s, n)| std::iter::repeat_n(b'A' + s, n)).collect()
            }
            _ => {
                let max = tc.draw_silent(gs::integers::<u32>().min_value(1));
                let vals =
                    tc.draw_silent(gs::vecs(gs::integers::<u32>().max_value(max)).max_size(600));
                vals.iter().flat_map(|v| v.to_le_bytes()).collect()
            }
        }
    }

    #[hegel::composite]
    fn arb_coder(tc: &TestCase) -> Coder {
        match tc.draw_silent(gs::integers::<u8>().max_value(5)) {
            0 => Coder::Cat,
            1 => Coder::Ext,
            _ => Coder::Range {
                order: tc.draw_silent(gs::integers::<u8>().max_value(3)),
                rle: tc.draw_silent(gs::booleans()),
            },
        }
    }

    #[hegel::composite]
    fn arb_plain(tc: &TestCase) -> Layout {
        Layout::Plain {
            coder: tc.draw_silent(arb_coder()),
            pack: tc.draw_silent(gs::booleans()),
            no_size: tc.draw_silent(gs::booleans()),
        }
    }

    /// A layout whose stripes nest at most `levels` deep.
    fn arb_layout(tc: &TestCase, levels: u8) -> Layout {
        if levels == 0 || !tc.draw_silent(gs::weighted_booleans(0.3)) {
            return tc.draw_silent(arb_plain());
        }
        let n = tc.draw_silent(gs::integers::<usize>().min_value(1).max_value(6));
        let junk = tc.draw_silent(gs::integers::<u8>()) & !FLAG_STRIPE;
        let parts = (0..n).map(|_| arb_layout(tc, levels - 1)).collect();
        Layout::Stripe { junk, parts }
    }

    // r[verify cram.codec.arith+2]
    // r[verify cram.codec.arith.wrapper]
    // r[verify cram.codec.arith.order]
    // r[verify cram.codec.arith.rle]
    // r[verify cram.codec.arith.pack]
    // r[verify cram.codec.arith.ext]
    // r[verify cram.codec.arith.stripe]
    #[hegel::test]
    fn round_trips_every_layout(tc: TestCase) {
        let data = tc.draw(arb_data());
        let layout = arb_layout(&tc, MAX_STRIPE_DEPTH);
        tc.note(&format!("layout: {layout:?}"));
        let stream = encode(&layout, &data);
        assert_eq!(decode(&stream, data.len()).unwrap(), data);
        assert_eq!(reference::decode(&stream, data.len()), Some(data), "reference decoder");
    }

    // r[verify cram.codec.arith.order]
    // r[verify cram.codec.arith.rle]
    #[hegel::test]
    fn range_coded_bodies_round_trip(tc: TestCase) {
        let data = tc.draw(arb_data());
        let order = tc.draw(gs::integers::<u8>().max_value(3));
        let rle = tc.draw(gs::booleans());
        let layout =
            Layout::Plain { coder: Coder::Range { order, rle }, pack: false, no_size: false };
        let stream = encode(&layout, &data);
        assert_eq!(decode(&stream, 0).unwrap(), data);
        assert_eq!(reference::decode(&stream, 0), Some(data), "reference decoder");
    }

    // r[verify cram.codec.arith.stripe]
    #[test]
    fn stripe_nesting_is_bounded() {
        let data = b"0123456789abcdefghijklmnopqrstuvwxyz".to_vec();
        let mut layout = Layout::Plain { coder: Coder::Cat, pack: false, no_size: true };
        for _ in 0..MAX_STRIPE_DEPTH {
            layout = Layout::Stripe { junk: 0, parts: vec![layout] };
        }
        assert_eq!(decode(&encode(&layout, &data), 0).unwrap(), data);
        let deeper = Layout::Stripe { junk: 0, parts: vec![layout] };
        assert!(matches!(
            decode(&encode(&deeper, &data), 0),
            Err(CramError::ArithStripeTooDeep { limit: MAX_STRIPE_DEPTH })
        ));
    }

    // r[verify cram.codec.arith.stripe]
    #[test]
    fn stripe_rejects_zero_streams_and_short_substreams() {
        assert!(matches!(decode(&[FLAG_STRIPE, 4, 0], 0), Err(CramError::ArithStripeZeroStreams)));
        // Two streams of 2 bytes each, but the second stores 1.
        let s = [FLAG_STRIPE, 4, 2, 4, 3, FLAG_CAT, 2, b'a', b'b', FLAG_CAT, 1, b'c'];
        assert!(matches!(
            decode(&s, 0),
            Err(CramError::ArithLengthMismatch { expected: 2, actual: 1 })
        ));
        // A compressed length past the end.
        let s = [FLAG_STRIPE, 2, 1, 9, FLAG_CAT, 2, b'a', b'b'];
        assert!(matches!(decode(&s, 0), Err(CramError::Truncated { .. })));
    }

    // r[verify cram.codec.arith.wrapper]
    #[test]
    fn no_size_takes_the_callers_length() {
        let s = [FLAG_CAT | FLAG_NO_SIZE, b'a', b'b', b'c'];
        assert_eq!(decode(&s, 3).unwrap(), b"abc");
        assert_eq!(decode(&s, 2).unwrap(), b"ab", "CAT copies the first bytes");
        // Without NOSZ the stored length wins over the caller's.
        assert_eq!(decode(&[FLAG_CAT, 2, b'a', b'b'], 7).unwrap(), b"ab");
        assert!(matches!(decode(&s, 4), Err(CramError::Truncated { .. })));
    }

    // r[verify cram.codec.arith.wrapper]
    #[test]
    fn an_empty_body_decodes_to_nothing() {
        assert_eq!(decode(&[0x00, 0], 0).unwrap(), b"");
        assert_eq!(decode(&[0x01 | FLAG_RLE, 0], 0).unwrap(), b"");
        assert!(matches!(
            decode(&[0x00, 5], 0),
            Err(CramError::ArithLengthMismatch { expected: 5, actual: 0 })
        ));
        // One symbol packs to nothing, so an empty body is its whole run.
        assert_eq!(decode(&[FLAG_PACK | FLAG_CAT, 4, 1, b'x', 0], 0).unwrap(), b"xxxx");
    }

    // r[verify cram.codec.arith.pack]
    #[test]
    fn pack_rejects_bad_lengths() {
        // A packed length past the uncompressed one.
        assert!(matches!(
            decode(&[FLAG_PACK | FLAG_CAT, 2, 2, b'a', b'b', 3, 0, 0, 0], 0),
            Err(CramError::ArithPackedLengthTooLarge { packed: 3, unpacked: 2 })
        ));
        // Nine 2-symbol values need two packed bytes; the body holds one.
        assert!(matches!(
            decode(&[FLAG_PACK | FLAG_CAT, 9, 2, b'a', b'b', 1, 0xff], 0),
            Err(CramError::ArithLengthMismatch { expected: 2, actual: 1 })
        ));
        assert_eq!(
            decode(&[FLAG_PACK | FLAG_CAT, 9, 2, b'a', b'b', 2, 0b0101_0110, 1], 0).unwrap(),
            b"abbababab"
        );
    }

    // r[verify cram.codec.arith.ext]
    #[test]
    fn ext_needs_bzip2_that_fits() {
        assert!(matches!(
            decode(&[FLAG_EXT, 3, 0x1f, 0x8b, 8, 0], 0),
            Err(CramError::ArithExtUnknownCodec)
        ));
        let mut s = vec![FLAG_EXT, 3];
        s.extend(bzip2(b"abcd"));
        assert!(matches!(
            decode(&s, 0),
            Err(CramError::ArithLengthMismatch { expected: 3, actual: 4 })
        ));
        let mut s = vec![FLAG_EXT, 4];
        s.extend(bzip2(b"abcd"));
        assert_eq!(decode(&s, 0).unwrap(), b"abcd");
    }

    // r[verify cram.codec.arith.order]
    #[test]
    fn truncated_range_coder_data_is_an_error() {
        let data = b"the quick brown fox jumps over the lazy dog".repeat(20);
        for rle in [false, true] {
            let layout = Layout::Plain {
                coder: Coder::Range { order: 1, rle },
                pack: false,
                no_size: false,
            };
            let stream = encode(&layout, &data);
            for cut in [4, 7, stream.len() / 2, stream.len() - 1] {
                let err = decode(&stream[..cut], 0).unwrap_err();
                assert!(
                    matches!(err, CramError::ArithCorruptData { .. } | CramError::Truncated { .. }),
                    "cut at {cut}: {err:?}"
                );
            }
        }
    }

    /// Where the production and reference decoders both accept a stream they
    /// agree on its bytes. Their verdicts can differ on invalid streams; the
    /// reference's module docs list where.
    fn assert_agree(stream: &[u8], size: usize) {
        if let (Ok(ours), Some(theirs)) = (decode(stream, size), reference::decode(stream, size)) {
            assert!(ours == theirs, "production and reference decoders disagree");
        }
    }

    // r[verify cram.codec.arith+2]
    #[hegel::test]
    fn arbitrary_input_never_panics(tc: TestCase) {
        let bytes = tc.draw(gs::binary().max_size(512));
        let size = tc.draw(gs::integers::<usize>().max_value(4096));
        assert_agree(&bytes, size);
    }

    /// Corrupting a valid stream reaches deeper paths than random bytes.
    // r[verify cram.codec.arith+2]
    #[hegel::test]
    fn corrupted_streams_never_panic(tc: TestCase) {
        let data = tc.draw(arb_data());
        let layout = arb_layout(&tc, 2);
        let mut stream = encode(&layout, &data);
        let flips = tc.draw(
            gs::vecs(gs::tuples!(gs::integers::<usize>(), gs::integers::<u8>().min_value(1)))
                .max_size(4),
        );
        for (at, x) in flips {
            let len = stream.len();
            stream[at % len] ^= x;
        }
        assert_agree(&stream, data.len());
    }
}
