//! Adaptive arithmetic coder encoder (CRAM 3.1, method 6): the inverse of
//! [`crate::cram::arith`], byte for byte as htscodecs' `arith_compress_to`
//! (`arith_dynamic.c`) writes it.

use super::rans::Order;
use super::rans_nx16::{pack, put_uint7};
use crate::cram::range_coder::{AdaptiveModel, RangeEncoder};

// r[impl cram_writer.arith]

const FLAG_ORDER: u8 = 0x01;
const FLAG_STRIPE: u8 = 0x08;
const FLAG_NO_SIZE: u8 = 0x10;
const FLAG_CAT: u8 = 0x20;
const FLAG_RLE: u8 = 0x40;
const FLAG_PACK: u8 = 0x80;

/// htscodecs' STRIPE width when the caller names none.
const STRIPES: usize = 4;
/// htscodecs drops STRIPE at or below this many bytes.
const STRIPE_MIN: usize = 20;
/// htscodecs codes order 1 as order 0 below this many bytes.
const ORDER_1_MIN: usize = 8;

/// Symbols of a run-length model: parts of 0..=3, where 3 continues the run.
const RUN_SYMBOLS: usize = 4;
const RUN_CONTINUES: usize = 3;
/// Run-length contexts: the literal for a run's first part, then
/// [`RUN_CTX_SECOND`] for its second and [`RUN_CTX_LATER`] for every later one.
const RUN_CONTEXTS: usize = 258;
const RUN_CTX_SECOND: usize = 256;
const RUN_CTX_LATER: usize = 257;

type LiteralModel = AdaptiveModel<256>;
type RunModel = AdaptiveModel<RUN_SYMBOLS>;

/// What to ask of [`encode`]. A transform that does not apply is dropped
/// from the stream's flags, as htscodecs drops it (`r[cram_writer.arith]`).
#[derive(Debug, Clone, Copy)]
pub(super) struct ArithFlags {
    pub order: Order,
    /// Four byte-interleaved substreams, each coded on its own.
    pub stripe: bool,
    /// Each literal's run length coded after it.
    pub rle: bool,
    /// Two, four or eight symbols a byte, for up to 16 distinct symbols.
    pub pack: bool,
}

impl ArithFlags {
    /// The combination a flag byte names, as htslib's `ARITH_PR*` method
    /// numbers do; its CAT and NOSZ bits are the encoder's to set.
    pub(super) const fn from_byte(byte: u8) -> Self {
        Self {
            order: if byte & FLAG_ORDER == 0 { Order::Zero } else { Order::One },
            stripe: byte & FLAG_STRIPE != 0,
            rle: byte & FLAG_RLE != 0,
            pack: byte & FLAG_PACK != 0,
        }
    }

    pub(super) const fn byte(self) -> u8 {
        let mut byte = 0;
        if matches!(self.order, Order::One) {
            byte |= FLAG_ORDER;
        }
        if self.stripe {
            byte |= FLAG_STRIPE;
        }
        if self.rle {
            byte |= FLAG_RLE;
        }
        if self.pack {
            byte |= FLAG_PACK;
        }
        byte
    }
}

/// The combinations htslib tries for a block above compression level 1
/// (`ARITH_PR0`, `PR64`, `PR128`, `PR192`, `PR1`, `PR9`, `PR129`, `PR193`).
#[allow(dead_code, reason = "the writer's codec choice plumbs the profile in")]
pub(super) const PROFILE: [ArithFlags; 8] = [
    ArithFlags::from_byte(0),
    ArithFlags::from_byte(FLAG_RLE),
    ArithFlags::from_byte(FLAG_PACK),
    ArithFlags::from_byte(FLAG_RLE | FLAG_PACK),
    ArithFlags::from_byte(FLAG_ORDER),
    ArithFlags::from_byte(FLAG_ORDER | FLAG_STRIPE),
    ArithFlags::from_byte(FLAG_ORDER | FLAG_PACK),
    ArithFlags::from_byte(FLAG_ORDER | FLAG_RLE | FLAG_PACK),
];

/// `data` as an arithmetic coder stream, or `None` when it is too large for
/// the format's sizes (htscodecs refuses more than `i32::MAX` bytes).
#[allow(dead_code, reason = "the writer's codec choice plumbs the encoder in")]
pub(super) fn encode(data: &[u8], flags: ArithFlags) -> Option<Vec<u8>> {
    i32::try_from(data.len()).ok()?;
    if flags.stripe && data.len() > STRIPE_MIN {
        encode_stripe(data, flags.byte())
    } else {
        encode_stream(data, flags.byte() & !FLAG_STRIPE, true)
    }
}

/// One stream without STRIPE, its flags in htscodecs' `order` byte: PACK,
/// then the range coder (with RLE), the body stored as CAT when coding does
/// not shrink it.
fn encode_stream(data: &[u8], mut byte: u8, sized: bool) -> Option<Vec<u8>> {
    let mut out = vec![0];
    if sized {
        put_uint7(&mut out, data.len())?;
    } else {
        byte |= FLAG_NO_SIZE;
    }

    let packed = if byte & FLAG_PACK != 0 { pack(data) } else { None };
    let mut body = data;
    match &packed {
        Some((meta, packed)) => {
            out.extend_from_slice(meta);
            put_uint7(&mut out, packed.len())?;
            body = packed;
        }
        None => byte &= !FLAG_PACK,
    }
    if body.is_empty() {
        byte &= !FLAG_RLE;
    }
    if body.len() < ORDER_1_MIN {
        byte &= !FLAG_ORDER;
    }

    let order_1 = byte & FLAG_ORDER != 0;
    let coded = if byte & FLAG_RLE != 0 {
        encode_rle(body, order_1)?
    } else {
        encode_order(body, order_1)?
    };
    if coded.len() < body.len() {
        out.extend_from_slice(&coded);
    } else {
        // htscodecs keeps the RLE bit, which a CAT body ignores.
        byte = (byte & !FLAG_ORDER) | FLAG_CAT;
        out.extend_from_slice(body);
    }
    if let Some(first) = out.first_mut() {
        *first = byte;
    }
    Some(out)
}

/// STRIPE: the bytes dealt round-robin into four substreams, each written
/// without its size by the smallest of htscodecs' candidates for its
/// position (the first of equals).
fn encode_stripe(data: &[u8], byte: u8) -> Option<Vec<u8>> {
    // htscodecs' candidates per substream (the last row for all past it);
    // order 1 only when the stream asks for it.
    const CANDIDATES: [&[u8]; 4] = [
        &[FLAG_ORDER, FLAG_RLE, 0],
        &[FLAG_ORDER, 0],
        &[FLAG_ORDER, FLAG_PACK],
        &[FLAG_ORDER, FLAG_PACK],
    ];
    let mut out = vec![byte];
    put_uint7(&mut out, data.len())?;
    out.push(u8::try_from(STRIPES).ok()?);
    let mut parts = Vec::with_capacity(STRIPES);
    for (j, candidates) in CANDIDATES.iter().enumerate() {
        let part: Vec<u8> = data.iter().skip(j).step_by(STRIPES).copied().collect();
        let mut best: Option<Vec<u8>> = None;
        for &c in candidates.iter().filter(|&&c| c & FLAG_ORDER == 0 || byte & FLAG_ORDER != 0) {
            let coded = encode_stream(&part, c, false)?;
            if best.as_ref().is_none_or(|b| coded.len() < b.len()) {
                best = Some(coded);
            }
        }
        let best = best?;
        put_uint7(&mut out, best.len())?;
        parts.push(best);
    }
    for part in parts {
        out.extend_from_slice(&part);
    }
    Some(out)
}

/// The symbol count every literal model covers: one past the largest byte,
/// written as a byte (256 as 0).
fn max_sym(data: &[u8]) -> usize {
    data.iter().max().map_or(1, |&m| usize::from(m).saturating_add(1))
}

/// The coded body's first byte, then the range coder starting after it.
#[allow(clippy::cast_possible_truncation, reason = "256 is written as 0, as htscodecs does")]
fn start(max_sym: usize) -> RangeEncoder {
    RangeEncoder::new(vec![max_sym as u8])
}

/// The literal models: one for order 0, one per previous byte for order 1.
/// Every byte is below `max_sym`, so contexts past it never occur.
fn literal_models(max_sym: usize, order_1: bool) -> Vec<LiteralModel> {
    vec![LiteralModel::new(max_sym); if order_1 { max_sym } else { 1 }]
}

/// htscodecs' `arith_compress_O0` / `arith_compress_O1`.
/// `None` never happens: every byte is below `max_sym`, so every context
/// has a model holding its symbol.
fn encode_order(data: &[u8], order_1: bool) -> Option<Vec<u8>> {
    let m = max_sym(data);
    let mut rc = start(m);
    let mut models = literal_models(m, order_1);
    let mut ctx = 0usize;
    for &b in data {
        models.get_mut(ctx)?.encode(&mut rc, u16::from(b))?;
        if order_1 {
            ctx = usize::from(b);
        }
    }
    Some(rc.finish())
}

/// htscodecs' `arith_compress_O0_RLE` / `arith_compress_O1_RLE`: a literal,
/// then its further copies in parts of at most 3 (the literal's context,
/// then [`RUN_CTX_SECOND`], then [`RUN_CTX_LATER`]), with a closing 0 after
/// a final part of 3. `None` never happens, as for [`encode_order`].
fn encode_rle(data: &[u8], order_1: bool) -> Option<Vec<u8>> {
    let m = max_sym(data);
    let mut rc = start(m);
    let mut literals = literal_models(m, order_1);
    let mut runs = vec![RunModel::new(RUN_SYMBOLS); RUN_CONTEXTS];
    let mut ctx = 0usize;
    let mut rest = data;
    while let Some((&lit, after)) = rest.split_first() {
        literals.get_mut(ctx)?.encode(&mut rc, u16::from(lit))?;
        if order_1 {
            ctx = usize::from(lit);
        }
        let mut run = after.iter().take_while(|&&b| b == lit).count();
        rest = after.get(run..).unwrap_or_default();

        let mut rctx = usize::from(lit);
        loop {
            let part = run.min(RUN_CONTINUES);
            run = run.saturating_sub(part);
            encode_run_part(&mut runs, rctx, part, &mut rc)?;
            rctx = if rctx == usize::from(lit) { RUN_CTX_SECOND } else { RUN_CTX_LATER };
            if part == RUN_CONTINUES && run == 0 {
                encode_run_part(&mut runs, rctx, 0, &mut rc)?;
            }
            if run == 0 {
                break;
            }
        }
    }
    Some(rc.finish())
}

fn encode_run_part(
    runs: &mut [RunModel],
    ctx: usize,
    part: usize,
    rc: &mut RangeEncoder,
) -> Option<()> {
    runs.get_mut(ctx)?.encode(rc, u16::try_from(part).ok()?)
}

#[cfg(test)]
#[allow(
    clippy::unwrap_used,
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    clippy::indexing_slicing,
    reason = "test code"
)]
mod tests {
    use super::*;
    use crate::cram::arith::decode;
    use crate::cram::arith::htscodecs_oracle::{arb_data, htscodecs_compress};
    use hegel::prelude::*;

    // r[verify cram_writer.arith]

    /// Every combination [`encode`] takes: order × STRIPE × RLE × PACK.
    fn every_flags() -> impl Iterator<Item = ArithFlags> {
        [0u8, 1, 8, 9, 64, 65, 72, 73, 128, 129, 136, 137, 192, 193, 200, 201]
            .into_iter()
            .map(ArithFlags::from_byte)
    }

    fn round_trip(data: &[u8], flags: ArithFlags) -> Vec<u8> {
        let stream = encode(data, flags).unwrap();
        let decoded = decode(&stream, 0).unwrap_or_else(|e| panic!("{flags:?}: {e}"));
        assert!(decoded == data, "{flags:?}: {} bytes decode to {}", data.len(), decoded.len());
        stream
    }

    #[hegel::test]
    fn every_combination_round_trips(tc: TestCase) {
        let data = tc.draw(arb_data());
        for flags in every_flags() {
            round_trip(&data, flags);
        }
    }

    /// htscodecs' `arith_compress_to`, linked through hts-sys, writes the
    /// same bytes for every combination.
    #[hegel::test(test_cases = crate::pinned::cases(200))]
    fn matches_htscodecs_byte_for_byte(tc: TestCase) {
        let data = tc.draw(arb_data());
        for flags in every_flags() {
            let ours = round_trip(&data, flags);
            let theirs = htscodecs_compress(&data, flags.byte().into());
            assert!(ours == theirs, "{flags:?} on {} bytes: {ours:?} vs {theirs:?}", data.len());
        }
    }

    #[test]
    fn edge_cases_match_htscodecs() {
        let inputs: [&[u8]; 8] = [
            b"",
            b"a",
            b"aaaaaaaa",
            b"abcdefg",
            &[0xFF; 30],
            &[0, 255, 0, 255, 0, 255, 0, 255, 0, 255, 0, 255, 0, 255, 0, 255, 0, 255, 0, 255, 0],
            &[7; 3000],
            b"the quick brown fox jumps over the lazy dog, again and again and again",
        ];
        for data in inputs {
            for flags in every_flags() {
                let ours = round_trip(data, flags);
                let theirs = htscodecs_compress(data, flags.byte().into());
                assert_eq!(ours, theirs, "{flags:?} on {data:?}");
            }
        }
    }

    /// Blocks long enough that models renormalise many times (a model's
    /// total passes `MAX_FREQ` after about 4000 codings): skewed bytes, a
    /// small alphabet in runs, and runs long enough for RLE's later contexts.
    #[test]
    fn large_inputs_match_htscodecs() {
        let mut x = 0x9e37_79b9_u32;
        let mut next = || {
            x ^= x << 13;
            x ^= x >> 17;
            x ^= x << 5;
            x
        };
        let skewed: Vec<u8> = (0..110_000).map(|_| (next() % 7 * (next() % 5)) as u8).collect();
        let runs: Vec<u8> = (0..20_000)
            .flat_map(|_| std::iter::repeat_n(b'A' + (next() % 4) as u8, 1 + next() as usize % 9))
            .collect();
        let flat: Vec<u8> = (0..100_000).map(|_| next() as u8).collect();
        for data in [skewed, runs, flat] {
            assert!(data.len() >= 100_000);
            for flags in every_flags() {
                let ours = round_trip(&data, flags);
                let theirs = htscodecs_compress(&data, flags.byte().into());
                assert!(ours == theirs, "{flags:?} on {} bytes", data.len());
            }
        }
    }

    #[test]
    fn profile_is_htslibs_methods() {
        let bytes: Vec<u8> = PROFILE.iter().map(|f| f.byte()).collect();
        assert_eq!(bytes, [0, 64, 128, 192, 1, 9, 129, 193]);
    }
}
