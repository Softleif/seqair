//! The byte-wise range coder and adaptive frequency model shared by the CRAM 3.1
//! arithmetic coder (method 6) and fqzcomp (method 7).
//!
//! Both follow htscodecs' `c_range_coder.h` and `c_simple_model.h` (Eugene
//! Shelwien's coder), which the spec's pseudocode (`CRAMcodecs` §4 "Range
//! coding", "Adaptive Modelling") describes.

use super::{codec_io::read_u8, reader::CramError};

/// Renormalise once the range drops below 2^24.
const TOP: u32 = 1 << 24;

/// A model's total frequency is halved once it exceeds this, so a `u16`
/// frequency never overflows (`2^16 - 17`).
pub(crate) const MAX_FREQ: u32 = (1 << 16) - 17;

/// What a decoded symbol adds to its frequency.
const STEP: u16 = 16;

/// The decoding side of the range coder.
// r[impl cram.codec.range_coder]
#[derive(Debug)]
pub(crate) struct RangeDecoder<'a> {
    src: &'a [u8],
    range: u32,
    code: u32,
}

impl<'a> RangeDecoder<'a> {
    /// Start decoding `src`, reading its first five bytes into `code`.
    ///
    /// The encoder's first byte is its initial (zero) cache byte, so shifting
    /// five bytes through a `u32` drops it, as in htscodecs' `RC_StartDecode`.
    pub(crate) fn new(mut src: &'a [u8]) -> Result<Self, CramError> {
        let mut code = 0u32;
        for _ in 0..5 {
            let b =
                read_u8(&mut src).ok_or(CramError::Truncated { context: "range coder start" })?;
            code = (code << 8) | u32::from(b);
        }
        Ok(Self { src, range: u32::MAX, code })
    }

    /// The bytes after what the decoder has consumed so far.
    #[cfg_attr(
        not(test),
        allow(dead_code, reason = "the arithmetic coder uses it; fqzcomp ignores trailing bytes")
    )]
    pub(crate) fn remaining(&self) -> &'a [u8] {
        self.src
    }

    /// Scale the range to `total` and return where `code` falls in it.
    ///
    /// A zero `total`, or one larger than the range, returns 0 without
    /// touching the range, as htscodecs' `RC_GetFreq` does; only corrupt
    /// input gets there, and the model rejects what follows.
    #[inline]
    pub(crate) fn get_freq(&mut self, total: u32) -> u32 {
        if total == 0 || self.range < total {
            return 0;
        }
        #[allow(
            clippy::arithmetic_side_effects,
            reason = "total > 0, and range >= total so the new range is >= 1"
        )]
        {
            self.range /= total;
            self.code / self.range
        }
    }

    /// Narrow the range to the symbol at `[cum, cum + freq)` of the total
    /// passed to the preceding [`get_freq`](Self::get_freq), and renormalise.
    ///
    /// Returns `None` when renormalising needs a byte past the end of the
    /// input, which a stream the encoder finished never does.
    #[inline]
    pub(crate) fn decode(&mut self, cum: u32, freq: u32) -> Option<()> {
        // Wrapping: on valid input `cum * range <= code` and `freq * range`
        // fits; corrupt input only produces wrong symbols, never a panic.
        self.code = self.code.wrapping_sub(cum.wrapping_mul(self.range));
        self.range = self.range.wrapping_mul(freq);
        while self.range < TOP {
            let b = read_u8(&mut self.src)?;
            self.code = (self.code << 8) | u32::from(b);
            self.range <<= 8;
        }
        Some(())
    }
}

/// One symbol and its frequency in an [`AdaptiveModel`].
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(crate) struct SymFreq {
    pub(crate) freq: u16,
    pub(crate) sym: u16,
}

/// An adaptive frequency model over the symbols `0..max_sym`, `max_sym <= N`.
///
/// Every symbol starts at frequency 1 and gains [`STEP`] each time it is
/// decoded. Entries are kept roughly most-frequent first — a decoded entry
/// swaps with its predecessor when it overtakes it — so the linear scan in
/// [`decode`](Self::decode) usually stops early. The scan order is part of the
/// format: the cumulative frequency of a symbol depends on it.
// r[impl cram.codec.adaptive_model]
#[derive(Debug, Clone)]
pub(crate) struct AdaptiveModel<const N: usize> {
    total: u32,
    /// Number of live entries in `syms`; the rest are never decoded.
    len: u16,
    syms: [SymFreq; N],
}

impl<const N: usize> AdaptiveModel<N> {
    /// A model over the symbols `0..max_sym`. `max_sym` is clamped to `N`.
    pub(crate) fn new(max_sym: usize) -> Self {
        const { assert!(N <= u16::MAX as usize) };
        let len = max_sym.min(N);
        let mut syms = [SymFreq::default(); N];
        for (i, s) in syms.iter_mut().enumerate().take(len) {
            #[allow(clippy::cast_possible_truncation, reason = "i < N <= u16::MAX")]
            {
                *s = SymFreq { freq: 1, sym: i as u16 };
            }
        }
        #[allow(clippy::cast_possible_truncation, reason = "len <= N <= u16::MAX")]
        Self { total: len as u32, len: len as u16, syms }
    }

    /// Decode one symbol from `rc` and update the model.
    ///
    /// Returns `None` if the input is corrupt (the coded value lies outside
    /// the model's total) or truncated.
    #[inline(always)]
    pub(crate) fn decode(&mut self, rc: &mut RangeDecoder<'_>) -> Option<u16> {
        let live = self.syms.get_mut(..usize::from(self.len))?;
        decode_symbol(rc, &mut self.total, live)
    }

    /// Halve every frequency, rounding up so none reaches zero.
    #[cfg(test)]
    fn renormalize(&mut self) {
        if let Some(live) = self.syms.get_mut(..usize::from(self.len)) {
            self.total = renormalize(live);
        }
    }
}

/// Decode one symbol of the adaptive model whose live entries are `syms` and
/// whose frequencies sum to `total`, and update the model — the body of
/// [`AdaptiveModel::decode`], for callers that store many models of a
/// run-time size in one arena (fqzcomp's 2^16 quality contexts).
///
/// Returns `None` if the input is corrupt (the coded value lies outside
/// `total`) or truncated.
#[inline]
pub(crate) fn decode_symbol(
    rc: &mut RangeDecoder<'_>,
    total: &mut u32,
    syms: &mut [SymFreq],
) -> Option<u16> {
    let target = rc.get_freq(*total);
    // All live frequencies sum to `total`, so a target inside it always
    // lands on a live entry.
    if target >= *total {
        return None;
    }
    let mut cum = 0u32;
    let mut i = 0usize;
    loop {
        let f = u32::from(syms.get(i)?.freq);
        #[allow(clippy::arithmetic_side_effects, reason = "cum <= total <= MAX_FREQ + STEP")]
        if cum + f > target {
            break;
        }
        #[allow(clippy::arithmetic_side_effects, reason = "as above; i < len")]
        {
            cum += f;
            i += 1;
        }
    }
    let entry = syms.get_mut(i)?;
    rc.decode(cum, u32::from(entry.freq))?;
    let sym = entry.sym;

    #[allow(
        clippy::arithmetic_side_effects,
        reason = "freq <= total <= MAX_FREQ < u16::MAX - STEP"
    )]
    {
        entry.freq += STEP;
        *total += u32::from(STEP);
    }
    if *total > MAX_FREQ {
        *total = renormalize(syms);
    }
    if let Some(prev) = i.checked_sub(1)
        && let (Some(&p), Some(&c)) = (syms.get(prev), syms.get(i))
        && c.freq > p.freq
    {
        syms.swap(prev, i);
    }
    Some(sym)
}

/// Halve every frequency, rounding up so none reaches zero; returns the new total.
fn renormalize(syms: &mut [SymFreq]) -> u32 {
    let mut total = 0u32;
    for s in syms {
        #[allow(
            clippy::arithmetic_side_effects,
            reason = "freq >> 1 <= freq; a model has far fewer than 2^16 entries, so the sum fits"
        )]
        {
            s.freq -= s.freq >> 1;
            total += u32::from(s.freq);
        }
    }
    total
}

#[cfg(test)]
#[allow(
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    clippy::indexing_slicing,
    reason = "test-only encoder and generators with bounded values"
)]
pub(crate) mod test_encoder {
    //! The encoding side, for round-trip tests only: `CRAMcodecs`
    //! `RangeEncode` / `RangeShiftLow` / `RangeEncodeEnd`, and the model's
    //! `encodeSymbol`.

    use super::{AdaptiveModel, MAX_FREQ, STEP, TOP};

    #[derive(Debug, Default)]
    pub(crate) struct RangeEncoder {
        low: u32,
        range: u32,
        ff_num: u32,
        cache: u32,
        carry: u32,
        pub(crate) out: Vec<u8>,
    }

    impl RangeEncoder {
        pub(crate) fn new() -> Self {
            Self { range: u32::MAX, ..Self::default() }
        }

        pub(crate) fn encode(&mut self, cum: u32, freq: u32, total: u32) {
            let old_low = self.low;
            self.range /= total;
            self.low = self.low.wrapping_add(cum * self.range);
            self.range *= freq;
            if self.low < old_low {
                self.carry = 1;
            }
            while self.range < TOP {
                self.range <<= 8;
                self.shift_low();
            }
        }

        fn shift_low(&mut self) {
            if self.low < 0xff00_0000 || self.carry != 0 {
                let (byte, fill) = if self.carry == 0 {
                    (self.cache as u8, 0xff)
                } else {
                    ((self.cache + 1) as u8, 0)
                };
                self.out.push(byte);
                for _ in 0..self.ff_num {
                    self.out.push(fill);
                }
                self.ff_num = 0;
                self.cache = self.low >> 24;
                self.carry = 0;
            } else {
                self.ff_num += 1;
            }
            self.low <<= 8;
        }

        pub(crate) fn finish(mut self) -> Vec<u8> {
            for _ in 0..5 {
                self.shift_low();
            }
            self.out
        }
    }

    impl<const N: usize> AdaptiveModel<N> {
        pub(crate) fn encode(&mut self, rc: &mut RangeEncoder, sym: u16) {
            let live = &mut self.syms[..usize::from(self.len)];
            let i = live.iter().position(|s| s.sym == sym).expect("symbol in model");
            let cum: u32 = live[..i].iter().map(|s| u32::from(s.freq)).sum();
            rc.encode(cum, u32::from(live[i].freq), self.total);
            live[i].freq += STEP;
            self.total += u32::from(STEP);
            if self.total > MAX_FREQ {
                self.renormalize();
            }
            if i > 0 && self.syms[i].freq > self.syms[i - 1].freq {
                self.syms.swap(i - 1, i);
            }
        }
    }
}

#[cfg(test)]
#[allow(
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    clippy::indexing_slicing,
    reason = "test code: bounded by generator ranges"
)]
mod tests {
    use super::test_encoder::RangeEncoder;
    use super::*;
    use hegel::prelude::*;

    // r[verify cram.codec.range_coder]
    // r[verify cram.codec.adaptive_model]
    #[hegel::test]
    fn model_round_trips_through_range_coder(tc: TestCase) {
        let max_sym = tc.draw(gs::integers::<usize>().min_value(1).max_value(256));
        // Skewed data: a small alphabet slice most of the time, so models
        // both renormalise (long runs of one symbol) and reorder.
        let hot = tc.draw(gs::integers::<usize>().min_value(1).max_value(max_sym));
        let data: Vec<u16> = tc.draw(
            gs::vecs(
                gs::integers::<usize>()
                    .min_value(0)
                    .max_value(max_sym - 1)
                    .map(move |s| (s % hot) as u16),
            )
            .max_size(20_000),
        );

        let mut enc = RangeEncoder::new();
        let mut m = AdaptiveModel::<256>::new(max_sym);
        for &s in &data {
            m.encode(&mut enc, s);
        }
        let bytes = enc.finish();

        let mut rc = RangeDecoder::new(&bytes).unwrap();
        let mut m = AdaptiveModel::<256>::new(max_sym);
        let decoded: Vec<u16> = (0..data.len()).map(|_| m.decode(&mut rc).unwrap()).collect();
        assert_eq!(decoded, data);
        assert!(rc.remaining().is_empty(), "decoder consumed exactly what the encoder wrote");
    }

    /// Interleaving several models over one coder is how every user of this
    /// module drives it; the contexts must stay independent.
    #[hegel::test]
    fn context_models_round_trip(tc: TestCase) {
        let data: Vec<u8> = tc.draw(gs::vecs(gs::integers::<u8>().max_value(3)).max_size(5_000));

        let mut enc = RangeEncoder::new();
        let mut ms = vec![AdaptiveModel::<4>::new(4); 4];
        let mut last = 0usize;
        for &s in &data {
            ms[last].encode(&mut enc, u16::from(s));
            last = usize::from(s);
        }
        let bytes = enc.finish();

        let mut rc = RangeDecoder::new(&bytes).unwrap();
        let mut ms = vec![AdaptiveModel::<4>::new(4); 4];
        let mut last = 0usize;
        for &s in &data {
            let d = ms[last].decode(&mut rc).unwrap();
            assert_eq!(d, u16::from(s));
            last = usize::from(d);
        }
    }

    #[test]
    fn model_renormalises_without_zero_frequencies() {
        let mut enc = RangeEncoder::new();
        let mut m = AdaptiveModel::<3>::new(3);
        // Enough of one symbol to cross MAX_FREQ several times.
        for _ in 0..20_000 {
            m.encode(&mut enc, 2);
        }
        assert!(m.total <= MAX_FREQ);
        assert!(m.syms.iter().all(|s| s.freq >= 1));
        assert_eq!(m.syms[0].sym, 2, "the hot symbol bubbled to the front");
        assert_eq!(m.total, m.syms.iter().map(|s| u32::from(s.freq)).sum::<u32>());
    }

    #[test]
    fn start_needs_five_bytes() {
        assert!(matches!(RangeDecoder::new(&[0; 4]), Err(CramError::Truncated { .. })));
        assert!(RangeDecoder::new(&[0; 5]).is_ok());
    }

    #[test]
    fn empty_model_rejects_every_input() {
        let bytes = [0u8; 16];
        let mut rc = RangeDecoder::new(&bytes).unwrap();
        let mut m = AdaptiveModel::<4>::new(0);
        assert_eq!(m.decode(&mut rc), None);
    }

    #[hegel::test]
    fn arbitrary_input_never_panics(tc: TestCase) {
        let bytes: Vec<u8> = tc.draw(gs::vecs(gs::integers::<u8>()).min_size(5).max_size(256));
        let max_sym = tc.draw(gs::integers::<usize>().max_value(256));
        let mut rc = RangeDecoder::new(&bytes).unwrap();
        let mut m = AdaptiveModel::<256>::new(max_sym);
        for _ in 0..1024 {
            if m.decode(&mut rc).is_none() {
                break;
            }
        }
    }

    /// Against the `CRAMcodecs` §4 pseudocode, transliterated in
    /// `fqzcomp_reference`: several models of different sizes interleaved
    /// over one coder give the same symbols and fail at the same point, on
    /// encoded streams and on noise.
    // r[verify cram.codec.range_coder]
    // r[verify cram.codec.adaptive_model]
    #[hegel::test]
    fn matches_the_spec_pseudocode(tc: TestCase) {
        use super::super::fqzcomp_reference as spec;

        let sizes: Vec<usize> = tc.draw(
            gs::vecs(gs::integers::<usize>().min_value(1).max_value(256)).min_size(1).max_size(4),
        );
        let n = tc.draw(gs::integers::<usize>().max_value(5_000));
        let bytes = if tc.draw(gs::booleans()) {
            // Skewed symbols, so models renormalise and reorder.
            let hot = tc.draw(gs::integers::<usize>().min_value(1).max_value(256));
            let picks = tc.draw(gs::binary().min_size(n).max_size(n));
            let mut enc = RangeEncoder::new();
            let mut ms: Vec<_> = sizes.iter().map(|&s| AdaptiveModel::<256>::new(s)).collect();
            for (i, &b) in picks.iter().enumerate() {
                let k = i % sizes.len();
                ms[k].encode(&mut enc, (usize::from(b) % hot % sizes[k]) as u16);
            }
            enc.finish()
        } else {
            tc.draw(gs::binary().max_size(64))
        };

        let (Ok(mut rc), Some(mut spec_rc)) =
            (RangeDecoder::new(&bytes), spec::RangeDecoder::new(&bytes))
        else {
            assert!(bytes.len() < 5, "both need exactly five bytes to start");
            return;
        };
        let mut ms: Vec<_> = sizes.iter().map(|&s| AdaptiveModel::<256>::new(s)).collect();
        let mut spec_ms: Vec<_> = sizes.iter().map(|&s| spec::Model::new(s)).collect();
        for i in 0..n.max(64) {
            let k = i % sizes.len();
            let got = ms[k].decode(&mut rc);
            assert_eq!(got, spec_ms[k].decode(&mut spec_rc), "symbol {i}");
            if got.is_none() {
                break;
            }
        }
    }
}
