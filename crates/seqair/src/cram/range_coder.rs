//! The byte-wise range coder and adaptive frequency model shared by the CRAM 3.1
//! arithmetic coder (method 6) and fqzcomp (method 7).
//!
//! Both follow htscodecs' `c_range_coder.h` and `c_simple_model.h` (Eugene
//! Shelwien's coder), which the spec's pseudocode (`CRAMcodecs` §4 "Range
//! coding", "Adaptive Modelling") describes.

use fearless_simd::{Simd, prelude::*, u8x16, u16x8};

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
    #[cfg(test)]
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
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, bytemuck::Pod, bytemuck::Zeroable)]
#[repr(C)]
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
    /// How deep recent symbols that missed the first [`SCALAR_PREFIX`]
    /// entries sat: an exponential moving average of their index, times 16.
    /// Picks the scalar or the vector search in
    /// [`decode_vectored`](Self::decode_vectored).
    depth: u16,
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
        Self { total: len as u32, len: len as u16, depth: 0, syms }
    }

    /// Decode one symbol from `rc` and update the model, scanning one
    /// entry at a time.
    ///
    /// Returns `None` if the input is corrupt (the coded value lies outside
    /// the model's total) or truncated.
    #[inline(always)]
    pub(crate) fn decode(&mut self, rc: &mut RangeDecoder<'_>) -> Option<u16> {
        let live = self.syms.get_mut(..usize::from(self.len))?;
        decode_symbol(rc, &mut self.total, live)
    }

    /// [`decode`](Self::decode), for models that may be wide and flat.
    ///
    /// The first [`SCALAR_PREFIX`] entries are scanned one at a time. Past
    /// them, a model whose symbols have lately sat [`DEEP`] in the scan order
    /// searches with vectors ([`find_vector`]); the rest keep scanning one
    /// at a time. Either way the entry found is the scan's.
    ///
    /// Returns `None` if the input is corrupt (the coded value lies outside
    /// the model's total) or truncated.
    // r[impl cram.codec.adaptive_model.simd_search]
    #[inline(always)]
    pub(crate) fn decode_vectored<S: Simd>(
        &mut self,
        simd: S,
        rc: &mut RangeDecoder<'_>,
    ) -> Option<u16> {
        let live = self.syms.get_mut(..usize::from(self.len))?;
        let target = rc.get_freq(self.total);
        if target >= self.total {
            return None;
        }
        let found = match find_prefix(live, target) {
            Ok(found) => found,
            Err(cum) => {
                let found = if N > BLOCK && self.depth >= DEEP {
                    find_vector(simd, live, cum, target)?
                } else {
                    find_scalar_from(live, SCALAR_PREFIX, cum, target)?
                };
                #[allow(
                    clippy::arithmetic_side_effects,
                    clippy::cast_possible_truncation,
                    reason = "depth <= 16 * 255 + 255 < 2^16 with the index capped at 255"
                )]
                {
                    self.depth = self.depth - (self.depth >> 4) + found.index.min(255) as u16;
                }
                found
            }
        };
        apply(rc, &mut self.total, live, found)
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
/// [`AdaptiveModel::decode`] with the scalar scan, for callers that store
/// many models of a run-time size in one arena (fqzcomp's 2^16 quality
/// contexts).
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
    let found = find_scalar(syms, target)?;
    apply(rc, total, syms, found)
}

/// Narrow the range to the entry `found` and update the model: add
/// [`STEP`] to the entry, renormalise past [`MAX_FREQ`], and swap the entry
/// forward if it overtook its predecessor. Returns its symbol.
#[inline(always)]
fn apply(
    rc: &mut RangeDecoder<'_>,
    total: &mut u32,
    syms: &mut [SymFreq],
    found: Found,
) -> Option<u16> {
    let i = found.index;
    rc.decode(found.cum, found.freq)?;
    let entry = syms.get_mut(i)?;
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

/// The model's scan, one entry at a time: the entry whose cumulative range
/// holds `target`. `None` if the frequencies sum to no more than `target`.
#[inline(always)]
fn find_scalar(syms: &[SymFreq], target: u32) -> Option<Found> {
    find_scalar_from(syms, 0, 0, target)
}

/// [`find_scalar`] from entry `index`, with `cum` the sum before it.
#[inline(always)]
fn find_scalar_from(
    syms: &[SymFreq],
    mut index: usize,
    mut cum: u32,
    target: u32,
) -> Option<Found> {
    loop {
        let freq = u32::from(syms.get(index)?.freq);
        #[allow(clippy::arithmetic_side_effects, reason = "cum <= total <= MAX_FREQ + STEP")]
        if cum + freq > target {
            return Some(Found { index, cum, freq });
        }
        #[allow(clippy::arithmetic_side_effects, reason = "as above; index < len")]
        {
            cum += freq;
            index += 1;
        }
    }
}

/// [`find_scalar`] over the first [`SCALAR_PREFIX`] entries: `Ok` if the
/// symbol is among them, else `Err` with their sum.
#[inline(always)]
fn find_prefix(syms: &[SymFreq], target: u32) -> Result<Found, u32> {
    let mut cum = 0u32;
    for (index, s) in syms.iter().enumerate().take(SCALAR_PREFIX) {
        let freq = u32::from(s.freq);
        #[allow(clippy::arithmetic_side_effects, reason = "cum <= total <= MAX_FREQ + STEP")]
        if cum + freq > target {
            return Ok(Found { index, cum, freq });
        }
        #[allow(clippy::arithmetic_side_effects, reason = "as above")]
        {
            cum += freq;
        }
    }
    Err(cum)
}

/// The [`AdaptiveModel::depth`] from which a model searches with vectors:
/// symbols that miss the scalar prefix sit 32 entries deep on average. Below
/// that the scalar scan wins on Apple M4 (measured over the hts-specs range
/// and tok3 vectors); wide flat models — packed qualities, the bytes of u32s,
/// name tokens — are far past it.
const DEEP: u16 = 16 * 32;

/// Entries checked one at a time before the vector search: on skewed data
/// the symbol is nearly always among the first few, where a predicted
/// branch beats the vector search's latency.
const SCALAR_PREFIX: usize = 4;

/// Entries per vector block of [`find_vector`].
const BLOCK: usize = 8;

/// Where the model's scan stops: the entry whose cumulative range holds the
/// coded value, the sum of the frequencies before it, and its frequency.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Found {
    index: usize,
    cum: u32,
    freq: u32,
}

/// [`find_scalar_from`] the entry after the first [`SCALAR_PREFIX`], whose
/// sum is `cum`, with vectors: the same entry, or `None`, as the scan.
/// [`find_in_blocks`] takes whole blocks, the entries past the last whole
/// block go one at a time.
#[inline(always)]
#[allow(
    clippy::arithmetic_side_effects,
    reason = "cumulative sums stay <= total <= MAX_FREQ + STEP, which fits u16 lanes and u32"
)]
fn find_vector<S: Simd>(simd: S, syms: &[SymFreq], cum: u32, target: u32) -> Option<Found> {
    let rest = syms.get(SCALAR_PREFIX..).unwrap_or_default();
    let (blocks, _) = rest.as_chunks::<BLOCK>();
    match find_in_blocks(simd, blocks, cum, target)? {
        Ok(found) => Some(Found { index: SCALAR_PREFIX + found.index, ..found }),
        Err(sum) => find_scalar_from(syms, SCALAR_PREFIX + blocks.len() * BLOCK, sum, target),
    }
}

/// The vector part of [`find_vector`], over blocks of eight entries: the
/// frequencies' running sums within the block (three shift-and-adds) plus
/// everything before it, compared against `target` all at once. The sums
/// only grow, so the lanes at or below `target` are the entries before the
/// symbol and their count is its index; the largest of those sums is its
/// `cum`, and the smallest sum past `target` is `cum + freq`. Blocks go in
/// pairs, with one exit test on the second: if the symbol is in the first,
/// no lane of the second is at or below `target`, so the same three
/// reductions over both blocks give the answer without a branch on which.
///
/// `Ok` with the index into `blocks`, or `Err` with the running sum after
/// the last block when `target` lies past them.
#[inline(always)]
#[allow(
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    reason = "running sums stay <= total <= MAX_FREQ + STEP < 2^16; a lane count fits usize"
)]
fn find_in_blocks<S: Simd>(
    simd: S,
    blocks: &[[SymFreq; BLOCK]],
    start: u32,
    target: u32,
) -> Option<Result<Found, u32>> {
    // Callers pass `target < total <= MAX_FREQ + STEP < 2^16`.
    let t = u16x8::splat(simd, u16::try_from(target).ok()?);
    let zero = u16x8::splat(simd, 0);
    let top = u16x8::splat(simd, u16::MAX);
    // Byte shuffle that broadcasts lane 7.
    let last =
        u8x16::from_slice(simd, &[14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15]);
    let mut carry = u16x8::splat(simd, u16::try_from(start).ok()?);
    let (pairs, rest) = blocks.as_chunks::<2>();
    for (b, [first, second]) in pairs.iter().enumerate() {
        let p0 = block_sums(simd, first, zero);
        let p1 = block_sums(simd, second, zero);
        let sums0 = p0 + carry;
        let carry1 = carry + p0.swizzle_dyn(last);
        let sums1 = p1 + carry1;
        let before1 = sums1.simd_le(t);
        let bits1 = before1.to_bitmask() & 0xFF;
        if bits1 != 0xFF {
            let before0 = sums0.simd_le(t);
            let bits0 = before0.to_bitmask() & 0xFF;
            let below = before0.select(sums0, carry).max(before1.select(sums1, carry));
            let above = before0.select(top, sums0).min(before1.select(top, sums1));
            let cum = below.reduce_max();
            let end = above.reduce_min();
            return Some(Ok(Found {
                index: 2 * b * BLOCK + (bits0 | (bits1 << BLOCK)).trailing_ones() as usize,
                cum: u32::from(cum),
                freq: u32::from(end - cum),
            }));
        }
        carry = carry1 + p1.swizzle_dyn(last);
    }
    if let Some(block) = rest.first() {
        let p = block_sums(simd, block, zero);
        let sums = p + carry;
        let before = sums.simd_le(t);
        let bits = before.to_bitmask() & 0xFF;
        if bits != 0xFF {
            let cum = before.select(sums, carry).reduce_max();
            let end = before.select(top, sums).reduce_min();
            return Some(Ok(Found {
                index: 2 * pairs.len() * BLOCK + bits.trailing_ones() as usize,
                cum: u32::from(cum),
                freq: u32::from(end - cum),
            }));
        }
        carry += p.swizzle_dyn(last);
    }
    // Every lane of `carry` holds the running sum.
    Some(Err(u32::from(carry.reduce_max())))
}

/// The running sums of one block's frequencies, lane `i` holding the sum of
/// entries `0..=i`.
#[inline(always)]
#[allow(clippy::arithmetic_side_effects, reason = "a model's frequencies sum to < 2^16")]
fn block_sums<S: Simd>(simd: S, block: &[SymFreq; BLOCK], zero: u16x8<S>) -> u16x8<S> {
    let raw: &[u16; 2 * BLOCK] = bytemuck::cast_ref(block);
    let (lo, hi) = raw.split_at(BLOCK);
    // `SymFreq` is `repr(C)` with `freq` first: the even halves.
    let f = u16x8::from_slice(simd, lo).unzip_low(u16x8::from_slice(simd, hi));
    let p = f + zero.slide::<7>(f);
    let p = p + zero.slide::<6>(p);
    p + zero.slide::<4>(p)
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
    use fearless_simd::{Level, dispatch};
    use hegel::prelude::*;

    /// Decode one symbol at `level`.
    fn decode_at<const N: usize>(
        level: Level,
        m: &mut AdaptiveModel<N>,
        rc: &mut RangeDecoder<'_>,
    ) -> Option<u16> {
        dispatch!(level, simd => m.decode_vectored(simd, rc))
    }

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

        for level in crate::simd_levels::levels() {
            let mut rc = RangeDecoder::new(&bytes).unwrap();
            let mut m = AdaptiveModel::<256>::new(max_sym);
            let decoded: Vec<u16> =
                (0..data.len()).map(|_| decode_at(level, &mut m, &mut rc).unwrap()).collect();
            assert_eq!(decoded, data, "{level:?}");
            assert!(rc.remaining().is_empty(), "decoder consumed exactly what the encoder wrote");
        }
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
            let d = decode_at(Level::new(), &mut ms[last], &mut rc).unwrap();
            assert_eq!(d, u16::from(s));
            last = usize::from(d);
        }
    }

    /// The model's scan as the spec states it, one entry at a time.
    fn scan(syms: &[SymFreq], target: u32) -> Option<Found> {
        let mut cum = 0u32;
        for (index, s) in syms.iter().enumerate() {
            if cum + u32::from(s.freq) > target {
                return Some(Found { index, cum, freq: u32::from(s.freq) });
            }
            cum += u32::from(s.freq);
        }
        None
    }

    /// Frequencies as a model holds them: at least 1 each, summing to at
    /// most `MAX_FREQ + STEP`; skewed or flat, so the symbol lands anywhere
    /// from the scalar prefix to the tail past the last whole block.
    #[hegel::composite]
    fn model_freqs(tc: &TestCase) -> Vec<SymFreq> {
        let n = tc.draw_silent(gs::integers::<usize>().min_value(1).max_value(256));
        let cap = tc.draw_silent(gs::sampled_from(vec![2u16, 40, 300, 4000]));
        let mut left = MAX_FREQ + u32::from(STEP) - n as u32;
        (0..n)
            .map(|i| {
                let extra = tc.draw_silent(gs::integers::<u16>().max_value(cap)).min(left as u16);
                left -= u32::from(extra);
                SymFreq { freq: 1 + extra, sym: i as u16 }
            })
            .collect()
    }

    // r[verify cram.codec.adaptive_model.simd_search]
    // r[verify io.simd_portable]
    #[hegel::test]
    fn symbol_search_matches_the_scan_at_every_level(tc: TestCase) {
        let syms = tc.draw(model_freqs().print_as_debug());
        let total: u32 = syms.iter().map(|s| u32::from(s.freq)).sum();
        // Past the total too: every level must find nothing there.
        let target = tc.draw(gs::integers::<u32>().max_value(total + 3));
        let expected = scan(&syms, target);
        for level in crate::simd_levels::levels() {
            let got = match find_prefix(&syms, target) {
                Ok(found) => Some(found),
                Err(cum) => dispatch!(level, simd => find_vector(simd, &syms, cum, target)),
            };
            assert_eq!(got, expected, "{level:?}, {} entries, target {target}", syms.len());
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
        assert_eq!(decode_at(Level::new(), &mut m, &mut rc), None);
    }

    #[hegel::test]
    fn arbitrary_input_never_panics(tc: TestCase) {
        let bytes: Vec<u8> = tc.draw(gs::vecs(gs::integers::<u8>()).min_size(5).max_size(256));
        let max_sym = tc.draw(gs::integers::<usize>().max_value(256));
        let mut rc = RangeDecoder::new(&bytes).unwrap();
        let mut m = AdaptiveModel::<256>::new(max_sym);
        for _ in 0..1024 {
            if decode_at(Level::new(), &mut m, &mut rc).is_none() {
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
    // r[verify cram.codec.adaptive_model.simd_search]
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

        let (Ok(_), Some(mut spec_rc)) =
            (RangeDecoder::new(&bytes), spec::RangeDecoder::new(&bytes))
        else {
            assert!(bytes.len() < 5, "both need exactly five bytes to start");
            return;
        };
        let mut spec_ms: Vec<_> = sizes.iter().map(|&s| spec::Model::new(s)).collect();
        let expected: Vec<Option<u16>> = (0..n.max(64))
            .map(|i| spec_ms[i % sizes.len()].decode(&mut spec_rc))
            .take_while(Option::is_some)
            .collect();
        for level in crate::simd_levels::levels() {
            let mut rc = RangeDecoder::new(&bytes).unwrap();
            let mut ms: Vec<_> = sizes.iter().map(|&s| AdaptiveModel::<256>::new(s)).collect();
            for i in 0..n.max(64) {
                let got = decode_at(level, &mut ms[i % sizes.len()], &mut rc);
                assert_eq!(got, expected.get(i).copied().flatten(), "{level:?} symbol {i}");
                if got.is_none() {
                    break;
                }
            }
        }
    }
}
