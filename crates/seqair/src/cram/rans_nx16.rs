// r[impl cram.codec.rans_nx16]
//! rANS Nx16 codec (CRAM compression method 5).
//!
//! N-way interleaved asymmetric numeral systems with 16-bit renormalization.
//! Supports order-0 and order-1 modes, plus STRIPE, RLE, and PACK transforms.

use super::reader::CramError;

const ALPHABET_SIZE: usize = 256;

const FLAG_ORDER: u8 = 0x01;
const FLAG_N32: u8 = 0x04;
const FLAG_STRIPE: u8 = 0x08;
const FLAG_NO_SIZE: u8 = 0x10;
const FLAG_CAT: u8 = 0x20;
const FLAG_RLE: u8 = 0x40;
const FLAG_PACK: u8 = 0x80;

/// Decode a rANS Nx16 compressed block.
pub fn decode(src: &[u8], mut uncompressed_size: usize) -> Result<Vec<u8>, CramError> {
    let mut cur: &[u8] = src;

    let flags = read_u8(&mut cur)?;
    let state_count = if flags & FLAG_N32 != 0 { 32 } else { 4 };

    if flags & FLAG_NO_SIZE == 0 {
        uncompressed_size = read_uint7(&mut cur)? as usize;
    }

    if flags & FLAG_STRIPE != 0 {
        return decode_stripe(&mut cur, uncompressed_size);
    }

    let bit_pack_ctx = if flags & FLAG_PACK != 0 {
        let (ctx, len) = read_bit_pack_context(&mut cur, uncompressed_size)?;
        uncompressed_size = len;
        Some(ctx)
    } else {
        None
    };

    let rle_ctx = if flags & FLAG_RLE != 0 {
        let (ctx, len) = read_rle_context(&mut cur, state_count, uncompressed_size)?;
        uncompressed_size = len;
        Some(ctx)
    } else {
        None
    };

    let mut dst = vec![0u8; uncompressed_size];

    if flags & FLAG_CAT != 0 {
        let data = split_off(&mut cur, uncompressed_size)?;
        dst.copy_from_slice(data);
    } else if flags & FLAG_ORDER == 0 {
        decode_order_0(&mut cur, &mut dst, state_count)?;
    } else {
        decode_order_1(&mut cur, &mut dst, state_count)?;
    }

    if let Some(ctx) = rle_ctx {
        dst = apply_rle(&dst, &ctx)?;
    }

    if let Some(ctx) = bit_pack_ctx {
        dst = apply_bit_unpack(&dst, &ctx)?;
    }

    Ok(dst)
}

// ── Primitive readers ────────────────────────────────────────────────

fn read_u8(src: &mut &[u8]) -> Result<u8, CramError> {
    let &b = src.first().ok_or(CramError::Truncated { context: "rans_nx16 u8" })?;
    *src = src.get(1..).ok_or(CramError::Truncated { context: "rans_nx16 u8" })?;
    Ok(b)
}

fn read_u16_le(src: &mut &[u8]) -> Result<u16, CramError> {
    let bytes: &[u8; 2] =
        src.first_chunk().ok_or(CramError::Truncated { context: "rans_nx16 u16" })?;
    let val = u16::from_le_bytes(*bytes);
    *src = src.get(2..).ok_or(CramError::Truncated { context: "rans_nx16 u16" })?;
    Ok(val)
}

fn read_u32_le(src: &mut &[u8]) -> Result<u32, CramError> {
    let bytes: &[u8; 4] =
        src.first_chunk().ok_or(CramError::Truncated { context: "rans_nx16 u32" })?;
    let val = u32::from_le_bytes(*bytes);
    *src = src.get(4..).ok_or(CramError::Truncated { context: "rans_nx16 u32" })?;
    Ok(val)
}

// r[impl cram.codec.uint7_bounded]
fn read_uint7(src: &mut &[u8]) -> Result<u32, CramError> {
    let mut n: u32 = 0;
    let mut count: u8 = 0;
    loop {
        if count >= 5 {
            return Err(CramError::Uint7Overflow);
        }
        count += 1;
        let b = read_u8(src)? as u32;
        n = (n << 7) | (b & 0x7f);
        if b & 0x80 == 0 {
            break;
        }
    }
    Ok(n)
}

fn split_off<'a>(src: &mut &'a [u8], len: usize) -> Result<&'a [u8], CramError> {
    let (head, rest) =
        src.split_at_checked(len).ok_or(CramError::Truncated { context: "rans_nx16 split_off" })?;
    *src = rest;
    Ok(head)
}

fn read_states(src: &mut &[u8], state_count: usize) -> Result<Vec<u32>, CramError> {
    (0..state_count).map(|_| read_u32_le(src)).collect()
}

// ── Core rANS step functions ─────────────────────────────────────────

fn state_cumulative_frequency(s: u32, bits: u32) -> u32 {
    s & ((1 << bits) - 1)
}

fn cumulative_frequencies_symbol(
    cumulative_frequencies: &[u32; ALPHABET_SIZE],
    frequency: u32,
) -> u8 {
    let mut sym: u8 = 0;
    while sym < 255
        && frequency >= *cumulative_frequencies.get(usize::from(sym) + 1).unwrap_or(&u32::MAX)
    {
        sym = sym.wrapping_add(1);
    }
    sym
}

// r[impl cram.codec.state_step_safety]
// SAFETY: spec guarantees g <= f * (s >> bits) + (s & mask) for valid frequency tables.
// Use wrapping_sub for release-mode robustness against malformed data.
fn state_step(s: u32, f: u32, g: u32, bits: u32) -> u32 {
    let result = f * (s >> bits) + (s & ((1 << bits) - 1));
    debug_assert!(result >= g, "state_step underflow: {result} < {g}");
    result.wrapping_sub(g)
}

fn state_renormalize(mut s: u32, src: &mut &[u8]) -> Result<u32, CramError> {
    if s < (1 << 15) {
        s = (s << 16) + u32::from(read_u16_le(src)?);
    }
    Ok(s)
}

// ── Alphabet reading ─────────────────────────────────────────────────

#[allow(
    clippy::indexing_slicing,
    reason = "sym is u8 so usize::from(sym) ≤ 255 < ALPHABET_SIZE=256"
)]
fn read_alphabet(src: &mut &[u8]) -> Result<[bool; ALPHABET_SIZE], CramError> {
    let mut alphabet = [false; ALPHABET_SIZE];

    let mut sym = read_u8(src)?;
    let mut prev_sym = sym;

    loop {
        alphabet[usize::from(sym)] = true;

        sym = read_u8(src)?;

        if sym == 0 {
            break;
        }

        if sym == prev_sym.wrapping_add(1) {
            let len = read_u8(src)?;
            for _ in 0..len {
                alphabet[usize::from(sym)] = true;
                sym = sym.wrapping_add(1);
            }
        }

        prev_sym = sym;
    }

    Ok(alphabet)
}

// ── Order-0 ──────────────────────────────────────────────────────────

const ORDER_0_BITS: u32 = 12;

fn decode_order_0(src: &mut &[u8], dst: &mut [u8], state_count: usize) -> Result<(), CramError> {
    let frequencies = read_frequencies_0(src)?;
    let cumulative_frequencies = build_cumulative_frequencies(&frequencies);
    let mut states = read_states(src, state_count)?;

    for chunk in dst.chunks_mut(states.len()) {
        for (d, state) in chunk.iter_mut().zip(states.iter_mut()) {
            let f = state_cumulative_frequency(*state, ORDER_0_BITS);
            let sym = cumulative_frequencies_symbol(&cumulative_frequencies, f);
            *d = sym;
            let i = usize::from(sym);
            *state = state_step(
                *state,
                *frequencies.get(i).unwrap_or(&0),
                *cumulative_frequencies.get(i).unwrap_or(&0),
                ORDER_0_BITS,
            );
            *state = state_renormalize(*state, src)?;
        }
    }

    Ok(())
}

fn read_frequencies_0(src: &mut &[u8]) -> Result<[u32; ALPHABET_SIZE], CramError> {
    let alphabet = read_alphabet(src)?;
    let mut frequencies = [0u32; ALPHABET_SIZE];

    for (i, frequency) in alphabet.iter().zip(&mut frequencies) {
        if *i {
            *frequency = read_uint7(src)?;
        }
    }

    normalize_frequencies(&mut frequencies, ORDER_0_BITS)?;
    Ok(frequencies)
}

// r[impl cram.codec.normalize_checked]
fn normalize_frequencies(
    frequencies: &mut [u32; ALPHABET_SIZE],
    bits: u32,
) -> Result<(), CramError> {
    let sum: u32 = frequencies.iter().copied().try_fold(0u32, |acc, f| acc.checked_add(f)).ok_or(
        CramError::FrequencyNormalizationOverflow {
            sum: frequencies.iter().copied().map(u64::from).sum(),
        },
    )?;

    if sum == 0 || sum == (1 << bits) {
        return Ok(());
    }

    let mut running = sum;
    let mut shift = 0u32;
    while running < (1 << bits) {
        running = running
            .checked_mul(2)
            .ok_or(CramError::FrequencyNormalizationOverflow { sum: u64::from(sum) << shift })?;
        shift += 1;
    }

    for f in frequencies {
        *f <<= shift;
    }

    Ok(())
}

fn build_cumulative_frequencies(frequencies: &[u32; ALPHABET_SIZE]) -> [u32; ALPHABET_SIZE] {
    let mut cumulative_frequencies = [0u32; ALPHABET_SIZE];
    let mut f = 0u32;

    for (i, cf) in cumulative_frequencies.iter_mut().enumerate().skip(1) {
        f += *frequencies.get(i - 1).unwrap_or(&0);
        *cf = f;
    }

    cumulative_frequencies
}

// ── Order-1 ──────────────────────────────────────────────────────────

type Frequencies1 = Box<[[u32; ALPHABET_SIZE]; ALPHABET_SIZE]>;
type CumulativeFrequencies1 = Box<[[u32; ALPHABET_SIZE]; ALPHABET_SIZE]>;

#[allow(
    clippy::indexing_slicing,
    reason = "k/l ≤ 255 (from u8 prev_sym/sym), arrays are [ALPHABET_SIZE=256]"
)]
fn decode_order_1(src: &mut &[u8], dst: &mut [u8], state_count: usize) -> Result<(), CramError> {
    let mut frequencies = Box::new([[0u32; ALPHABET_SIZE]; ALPHABET_SIZE]);
    let bits = read_frequencies_1(src, &mut frequencies)?;
    let cumulative_frequencies = build_cumulative_frequencies_1(&frequencies);

    let mut states = read_states(src, state_count)?;
    let mut prev_syms = vec![0u8; state_count];

    let chunk_size = dst.len() / state_count;

    for i in 0..chunk_size {
        for (j, (state, prev_sym)) in states.iter_mut().zip(&mut prev_syms).enumerate() {
            let k = usize::from(*prev_sym);
            let f = state_cumulative_frequency(*state, bits);
            let sym = cumulative_frequencies_symbol(
                cumulative_frequencies
                    .get(k)
                    .ok_or(CramError::Truncated { context: "rans_nx16 order-1 cumfreq" })?,
                f,
            );

            let out_idx = j * chunk_size + i;
            *dst.get_mut(out_idx)
                .ok_or(CramError::Truncated { context: "rans_nx16 order-1 output" })? = sym;

            let l = usize::from(sym);
            *state = state_step(*state, frequencies[k][l], cumulative_frequencies[k][l], bits);
            *state = state_renormalize(*state, src)?;
            *prev_sym = sym;
        }
    }

    let last_chunk = &mut dst[chunk_size * state_count..];
    if !last_chunk.is_empty() {
        let mut state = *states
            .last()
            .ok_or(CramError::Truncated { context: "rans_nx16 order-1 last state" })?;
        let mut prev_sym = *prev_syms
            .last()
            .ok_or(CramError::Truncated { context: "rans_nx16 order-1 last prev_sym" })?;

        for d in last_chunk {
            let k = usize::from(prev_sym);
            let f = state_cumulative_frequency(state, bits);
            let sym = cumulative_frequencies_symbol(&cumulative_frequencies[k], f);
            *d = sym;
            let l = usize::from(sym);
            state = state_step(state, frequencies[k][l], cumulative_frequencies[k][l], bits);
            state = state_renormalize(state, src)?;
            prev_sym = sym;
        }
    }

    Ok(())
}

fn read_frequencies_1(src: &mut &[u8], frequencies: &mut Frequencies1) -> Result<u32, CramError> {
    let n = read_u8(src)?;
    let bits = u32::from(n >> 4);
    let is_compressed = (n & 0x01) != 0;

    if is_compressed {
        let uncompressed_size = read_uint7(src)? as usize;
        let compressed_size = read_uint7(src)? as usize;
        let mut compressed_data = split_off(src, compressed_size)?;
        let mut tmp = vec![0u8; uncompressed_size];
        decode_order_0(&mut compressed_data, &mut tmp, 4)?;
        read_frequencies_1_inner(&mut &tmp[..], frequencies, bits)?;
    } else {
        read_frequencies_1_inner(src, frequencies, bits)?;
    }

    Ok(bits)
}

#[allow(
    clippy::indexing_slicing,
    reason = "ctx_idx/sym_idx come from enumerate() over [bool;256], so ≤ 255 < ALPHABET_SIZE=256"
)]
fn read_frequencies_1_inner(
    src: &mut &[u8],
    frequencies: &mut Frequencies1,
    bits: u32,
) -> Result<(), CramError> {
    let alphabet = read_alphabet(src)?;

    for (ctx_idx, ctx_active) in alphabet.iter().enumerate() {
        if !*ctx_active {
            continue;
        }

        let fs = &mut frequencies[ctx_idx];
        let mut sym_iter = alphabet.iter().enumerate().filter(|(_, b)| **b).peekable();

        while let Some((sym_idx, _)) = sym_iter.next() {
            let f = read_uint7(src)?;
            fs[sym_idx] = f;

            if f == 0 {
                let n = read_u8(src)? as usize;
                for _ in 0..n {
                    let _ = sym_iter.next();
                }
            }
        }

        normalize_frequencies(fs, bits)?;
    }

    Ok(())
}

fn build_cumulative_frequencies_1(frequencies: &Frequencies1) -> CumulativeFrequencies1 {
    let mut cf = Box::new([[0u32; ALPHABET_SIZE]; ALPHABET_SIZE]);
    for (f, g) in frequencies.iter().zip(cf.iter_mut()) {
        *g = build_cumulative_frequencies(f);
    }
    cf
}

// ── Stripe transform ─────────────────────────────────────────────────

fn decode_stripe(src: &mut &[u8], uncompressed_size: usize) -> Result<Vec<u8>, CramError> {
    let chunk_count = read_u8(src)? as usize;
    if chunk_count == 0 {
        return Err(CramError::RansStripeZeroChunks);
    }

    let compressed_sizes: Vec<usize> =
        (0..chunk_count).map(|_| read_uint7(src).map(|n| n as usize)).collect::<Result<_, _>>()?;

    let uncompressed_sizes: Vec<usize> = (0..chunk_count)
        .map(|i| {
            let (q, r) = (uncompressed_size / chunk_count, uncompressed_size % chunk_count);
            if r > i { q + 1 } else { q }
        })
        .collect();

    let chunks: Vec<Vec<u8>> = compressed_sizes
        .iter()
        .zip(&uncompressed_sizes)
        .map(|(&cs, &us)| {
            let buf = split_off(src, cs)?;
            decode(buf, us)
        })
        .collect::<Result<_, _>>()?;

    let mut dst = vec![0u8; uncompressed_size];
    for (i, chunk) in chunks.iter().enumerate() {
        for (j, &s) in chunk.iter().enumerate() {
            let idx = j * chunk_count + i;
            if let Some(d) = dst.get_mut(idx) {
                *d = s;
            }
        }
    }

    Ok(dst)
}

// ── Bit-pack transform ──────────────────────────────────────────────

struct BitPackContext {
    symbol_count: usize,
    mapping_table: Vec<u8>,
    uncompressed_size: usize,
}

fn read_bit_pack_context(
    src: &mut &[u8],
    uncompressed_size: usize,
) -> Result<(BitPackContext, usize), CramError> {
    let symbol_count = read_u8(src)? as usize;
    if symbol_count == 0 {
        return Err(CramError::RansBitPackZeroSymbols);
    }
    let mapping_table = split_off(src, symbol_count)?.to_vec();
    let packed_len = read_uint7(src)? as usize;

    Ok((BitPackContext { symbol_count, mapping_table, uncompressed_size }, packed_len))
}

fn apply_bit_unpack(src: &[u8], ctx: &BitPackContext) -> Result<Vec<u8>, CramError> {
    let mut dst = vec![0u8; ctx.uncompressed_size];

    match ctx.symbol_count {
        1 => {
            let sym = *ctx
                .mapping_table
                .first()
                .ok_or(CramError::Truncated { context: "bit_pack mapping" })?;
            dst.fill(sym);
        }
        2 => unpack(src, &ctx.mapping_table, 8, &mut dst),
        3..=4 => unpack(src, &ctx.mapping_table, 4, &mut dst),
        5..=16 => unpack(src, &ctx.mapping_table, 2, &mut dst),
        n => {
            return Err(CramError::RansBitPackTooManySymbols { symbol_count: n });
        }
    }

    Ok(dst)
}

fn unpack(src: &[u8], mapping_table: &[u8], chunk_size: usize, dst: &mut [u8]) {
    let bits = u8::BITS as usize;
    let shift = bits / chunk_size;
    let mask: u8 = (1u8 << shift).wrapping_sub(1);

    for (mut s, chunk) in src.iter().copied().zip(dst.chunks_mut(chunk_size)) {
        for d in chunk {
            let idx = usize::from(s & mask);
            *d = mapping_table.get(idx).copied().unwrap_or(0);
            s >>= shift;
        }
    }
}

// ── RLE transform ────────────────────────────────────────────────────

struct RleContext {
    rle_meta: Vec<u8>,
    output_len: usize,
}

fn read_rle_context(
    src: &mut &[u8],
    state_count: usize,
    uncompressed_size: usize,
) -> Result<(RleContext, usize), CramError> {
    let header = read_uint7(src)?;
    let context_size = (header >> 1) as usize;
    let is_compressed = (header & 0x01) == 0;

    let rle_encoded_len = read_uint7(src)? as usize;

    let rle_meta = if is_compressed {
        let compressed_size = read_uint7(src)? as usize;
        let mut buf = split_off(src, compressed_size)?;
        let mut tmp = vec![0u8; context_size];
        decode_order_0(&mut buf, &mut tmp, state_count)?;
        tmp
    } else {
        split_off(src, context_size)?.to_vec()
    };

    Ok((RleContext { rle_meta, output_len: uncompressed_size }, rle_encoded_len))
}

#[allow(
    clippy::indexing_slicing,
    reason = "sym is u8 so usize::from(sym) ≤ 255 < ALPHABET_SIZE=256"
)]
fn apply_rle(src: &[u8], ctx: &RleContext) -> Result<Vec<u8>, CramError> {
    let mut meta_cur: &[u8] = &ctx.rle_meta;

    let rle_alphabet = read_rle_alphabet(&mut meta_cur)?;

    let mut dst = vec![0u8; ctx.output_len];
    let mut dst_iter = dst.iter_mut();
    let mut src_iter = src.iter();

    while let Some(d) = dst_iter.next() {
        let &sym = src_iter.next().ok_or(CramError::Truncated { context: "rans_nx16 rle src" })?;
        *d = sym;

        if rle_alphabet[usize::from(sym)] {
            let len = read_uint7(&mut meta_cur)? as usize;
            for e in dst_iter.by_ref().take(len) {
                *e = sym;
            }
        }
    }

    Ok(dst)
}

#[allow(
    clippy::indexing_slicing,
    reason = "sym is u8 so usize::from(sym) ≤ 255 < ALPHABET_SIZE=256"
)]
fn read_rle_alphabet(src: &mut &[u8]) -> Result<[bool; ALPHABET_SIZE], CramError> {
    let mut alphabet = [false; ALPHABET_SIZE];

    let n = read_u8(src)? as usize;
    let symbol_count = if n == 0 { ALPHABET_SIZE } else { n };

    for _ in 0..symbol_count {
        let sym = read_u8(src)?;
        alphabet[usize::from(sym)] = true;
    }

    Ok(alphabet)
}

#[cfg(test)]
mod tests {
    use super::*;

    // r[verify cram.codec.rans_nx16]

    #[test]
    fn rans_stripe_zero_chunks_returns_error() {
        // Build a byte stream with FLAG_STRIPE set and chunk_count = 0
        let src = vec![
            FLAG_STRIPE, // flags = STRIPE
            4u8,         // uncompressed_size = 4 (uint7)
            0u8,         // chunk_count = 0 → triggers RansStripeZeroChunks
        ];

        let err = decode(&src, 0).unwrap_err();
        assert!(matches!(err, CramError::RansStripeZeroChunks));
    }

    #[test]
    fn rans_bit_pack_zero_symbols_returns_error() {
        // Build a stream with FLAG_PACK and symbol_count = 0
        let src = vec![
            FLAG_PACK, // flags = PACK
            5u8,       // uncompressed_size = 5 (uint7)
            0u8,       // symbol_count = 0 → triggers RansBitPackZeroSymbols
        ];

        let err = decode(&src, 0).unwrap_err();
        assert!(matches!(err, CramError::RansBitPackZeroSymbols));
    }

    #[test]
    fn rans_bit_pack_too_many_symbols_returns_error() {
        // symbol_count = 17 (> 16) — triggers RansBitPackTooManySymbols.
        // After reading symbol_count, the code reads mapping_table (symbol_count bytes)
        // and packed_len (uint7). We need enough bytes to get past those reads.
        let symbol_count: u8 = 17;
        let packed_len: u8 = 0; // uint7 = 0
        let mut src = Vec::new();
        src.push(FLAG_PACK); // flags = PACK
        // FLAG_NO_SIZE not set → read uncompressed_size as uint7
        src.push(10u8); // uncompressed_size
        src.push(symbol_count); // symbol_count = 17
        // mapping_table: 17 bytes
        src.extend(std::iter::repeat_n(b'A', symbol_count as usize));
        src.push(packed_len); // packed_len uint7

        // The rest of the stream would be consumed for the inner rANS decode, but the
        // error fires before that in apply_bit_unpack. We need the CAT flag set
        // to bypass the rANS decode and reach apply_bit_unpack directly.
        // Rebuild with CAT | PACK so we skip the rANS step.
        let mut src2 = Vec::new();
        src2.push(FLAG_PACK | FLAG_CAT); // flags = PACK | CAT
        src2.push(10u8); // uncompressed_size
        src2.push(symbol_count); // symbol_count = 17
        src2.extend(std::iter::repeat_n(b'A', symbol_count as usize));
        src2.push(packed_len); // packed_len uint7
        // CAT data: 0 bytes (packed_len = 0)

        let err = decode(&src2, 0).unwrap_err();
        assert!(
            matches!(err, CramError::RansBitPackTooManySymbols { symbol_count: 17 }),
            "expected RansBitPackTooManySymbols, got: {err:?}"
        );
    }

    #[test]
    fn decode_order_0_noodles_test_vector() {
        let src = [
            0x00, // flags = {empty}
            0x07, // uncompressed len = 7
            0x64, 0x65, 0x00, 0x6c, 0x6e, 0x6f, 0x00, 0x73, 0x00, 0x01, 0x01, 0x01, 0x01, 0x03,
            0x01, 0x00, 0x26, 0x20, 0x00, 0x00, 0xb8, 0x0a, 0x00, 0x00, 0xd8, 0x0a, 0x00, 0x00,
            0x00, 0x04, 0x00,
        ];
        assert_eq!(decode(&src, 0).unwrap(), b"noodles");
    }

    #[test]
    fn decode_order_1_noodles_test_vector() {
        let src = [
            0x01, // flags = ORDER
            0x4d, // uncompressed len = 77
            0xa0, 0x00, 0x64, 0x65, 0x00, 0x6c, 0x6e, 0x6f, 0x00, 0x73, 0x00, 0x00, 0x00, 0x01,
            0x01, 0x00, 0x00, 0x01, 0x01, 0x00, 0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x01, 0x00,
            0x02, 0x00, 0x01, 0x0f, 0x00, 0x02, 0x01, 0x00, 0x01, 0x01, 0x0f, 0x00, 0x02, 0x00,
            0x03, 0x0f, 0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x02, 0x0f, 0x00, 0x00, 0x00,
            0x05, 0x10, 0x80, 0x72, 0x60, 0x00, 0x80, 0x8b, 0x5f, 0x00, 0xc0, 0xb0, 0x60, 0x00,
            0x40, 0x49, 0x39, 0x00,
        ];
        assert_eq!(
            decode(&src, 0).unwrap(),
            b"nnnnnnnnnnnnooooooooooooooooddddddddddddddllllllllllllllleeeeeeeeeessssssssss"
        );
    }

    #[test]
    fn decode_stripe_noodles_test_vector() {
        let src = [
            0x08, // flags = STRIPE
            0x07, // uncompressed len = 7
            0x04, 0x17, 0x17, 0x17, 0x15, 0x00, 0x02, 0x6c, 0x6e, 0x00, 0x01, 0x01, 0x00, 0x08,
            0x01, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00,
            0x00, 0x02, 0x65, 0x6f, 0x00, 0x01, 0x01, 0x00, 0x08, 0x01, 0x00, 0x00, 0x00, 0x01,
            0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x02, 0x6f, 0x73, 0x00,
            0x01, 0x01, 0x00, 0x00, 0x01, 0x00, 0x00, 0x08, 0x01, 0x00, 0x00, 0x80, 0x00, 0x00,
            0x00, 0x80, 0x00, 0x00, 0x00, 0x01, 0x64, 0x00, 0x01, 0x00, 0x80, 0x00, 0x00, 0x00,
            0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x02, 0x00,
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x22, 0x00, 0x81, 0x11, 0x01, 0x7f, 0x00,
        ];
        assert_eq!(decode(&src, 0).unwrap(), b"noodles");
    }

    #[test]
    fn decode_uncompressed_noodles_test_vector() {
        let src = [
            0x20, // flags = CAT
            0x07, // uncompressed len = 7
            0x6e, 0x6f, 0x6f, 0x64, 0x6c, 0x65, 0x73,
        ];
        assert_eq!(decode(&src, 0).unwrap(), b"noodles");
    }

    #[test]
    fn decode_rle_noodles_test_vector() {
        let src = [
            0x40, // flags = RLE
            0x0d, // uncompressed len = 13
            0x06, 0x06, 0x17, 0x01, 0x07, 0x6f, 0x00, 0x02, 0x01, 0x01, 0x00, 0x00, 0x01, 0x00,
            0x00, 0x0c, 0x02, 0x00, 0x00, 0x08, 0x02, 0x00, 0x00, 0x80, 0x00, 0x00, 0x64, 0x65,
            0x00, 0x6c, 0x6e, 0x6f, 0x00, 0x73, 0x00, 0x03, 0x01, 0x01, 0x01, 0x01, 0x01, 0x00,
            0x3a, 0x20, 0x00, 0x00, 0x7c, 0x20, 0x00, 0x00, 0x52, 0x01, 0x00, 0x00, 0x08, 0x04,
            0x00,
        ];
        assert_eq!(decode(&src, 0).unwrap(), b"noooooooodles");
    }

    #[test]
    fn decode_bit_packing_noodles_test_vector() {
        let src = [
            0x80, // flags = PACK
            0x07, // uncompressed len = 7
            0x06, 0x64, 0x65, 0x6c, 0x6e, 0x6f, 0x73, 0x04, 0x04, 0x05, 0x00, 0x12, 0x43, 0x00,
            0x01, 0x01, 0x01, 0x01, 0x00, 0x0c, 0x02, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x08,
            0x02, 0x00, 0x00, 0x04, 0x02, 0x00,
        ];
        assert_eq!(decode(&src, 0).unwrap(), b"noodles");
    }

    // r[verify cram.codec.uint7_bounded]
    #[test]
    fn read_uint7_overflow_returns_error() {
        // 6 continuation bytes (all with high bit set) exceed the 5-iteration limit.
        let mut src: &[u8] = &[0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x00];
        let err = read_uint7(&mut src).unwrap_err();
        assert!(matches!(err, CramError::Uint7Overflow));
    }

    #[test]
    fn read_uint7_five_bytes_ok() {
        // Exactly 5 continuation bytes is the maximum valid sequence.
        let mut src: &[u8] = &[0x80, 0x80, 0x80, 0x80, 0x01];
        let result = read_uint7(&mut src);
        assert!(result.is_ok());
    }

    // r[verify cram.codec.normalize_checked]
    #[test]
    fn normalize_frequencies_overflow_returns_error() {
        // Fill frequencies so sum overflows u32.
        let mut frequencies = [0u32; ALPHABET_SIZE];
        // 256 * 0x01_000_000 = 0x1_0000_0000 which overflows u32.
        frequencies.fill(0x0100_0000);
        let result = normalize_frequencies(&mut frequencies, ORDER_0_BITS);
        assert!(
            matches!(result, Err(CramError::FrequencyNormalizationOverflow { .. })),
            "expected overflow error, got: {result:?}"
        );
    }
}
