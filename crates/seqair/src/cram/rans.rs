//! rANS 4x8 codec (CRAM compression method 4).
//!
//! 4-way interleaved asymmetric numeral systems with 8-bit renormalization.
//! Supports order-0 and order-1 modes.

// `clippy::unnecessary_lazy_evaluations` flags every `ok_or_else(||
// CramError::Truncated { ... })` here, because the variant's
// construction is "cheap". On the per-byte hot path it isn't: the
// eager `ok_or` form has the compiler build (and drop) a full-sized
// `CramError` on every successful read — the type's Drop has to
// dispatch on the discriminant because other variants own `PathBuf` /
// `std::io::Error` / `SmolStr`. samply showed
// `core::ptr::drop_in_place<CramError>` next to `decode_order_1` /
// `read_u8` for exactly this. Keep the lazy form here.
#![allow(
    clippy::unnecessary_lazy_evaluations,
    reason = "lazy form avoids per-call drop_in_place<CramError> on hot path"
)]

use super::reader::CramError;

const ALPHABET_SIZE: usize = 256;
const LOWER_BOUND: u32 = 1 << 23;

/// Reusable allocations for rANS 4x8 order-1 decoding.
///
/// Order-1 needs three large tables per block: frequency array (128 KB),
/// cumulative frequency array (128 KB), and symbol lookup table (1 MB).
/// Allocating them fresh for every block dominated CRAM decode profiles.
/// This struct holds one of each; zero it or overwrite it per block.
pub(crate) struct Rans4x8Buf {
    pub freq: Box<[[u16; ALPHABET_SIZE]; ALPHABET_SIZE]>,
    pub cum_freq: Box<[[u16; ALPHABET_SIZE]; ALPHABET_SIZE]>,
    pub sym_tables: Box<[[u8; 4096]; ALPHABET_SIZE]>,
}

impl Rans4x8Buf {
    pub fn new() -> Self {
        Self {
            freq: Box::new([[0u16; ALPHABET_SIZE]; ALPHABET_SIZE]),
            cum_freq: Box::new([[0u16; ALPHABET_SIZE]; ALPHABET_SIZE]),
            sym_tables: Box::new([[0u8; 4096]; ALPHABET_SIZE]),
        }
    }
}

impl Default for Rans4x8Buf {
    fn default() -> Self {
        Self::new()
    }
}

// r[impl cram.codec.rans4x8]
/// Decode a rANS 4x8 compressed block (allocating path — for callers
/// without a reusable buffer). Prefer [`decode_with_buf`] on the hot path.
pub fn decode(src: &[u8]) -> Result<Vec<u8>, CramError> {
    let mut buf = Rans4x8Buf::new();
    decode_with_buf(src, &mut buf)
}

/// Decode a rANS 4x8 compressed block, reusing the caller-owned
/// [`Rans4x8Buf`] for order-1 tables. Order-0 blocks don't touch the buf.
pub(crate) fn decode_with_buf(src: &[u8], buf: &mut Rans4x8Buf) -> Result<Vec<u8>, CramError> {
    let mut cur: &[u8] = src;

    // Header: order (u8), compressed_size (u32 LE), uncompressed_size (u32 LE)
    let order = read_u8(&mut cur).ok_or_else(|| CramError::Truncated { context: "rans header" })?;
    let _compressed_size =
        read_u32_le(&mut cur).ok_or_else(|| CramError::Truncated { context: "rans header" })?;
    let uncompressed_size = read_u32_le(&mut cur)
        .ok_or_else(|| CramError::Truncated { context: "rans header" })?
        as usize;

    super::reader::check_alloc_size(uncompressed_size, "rANS 4x8 output")?;
    let mut dst = vec![0u8; uncompressed_size];

    match order {
        0 => decode_order_0(&mut cur, &mut dst)?,
        1 => decode_order_1_buf(&mut cur, &mut dst, buf)?,
        _ => return Err(CramError::InvalidRansOrder { order }),
    }

    Ok(dst)
}

#[allow(
    clippy::indexing_slicing,
    reason = "indices are bounded: f ≤ 4095 (12-bit mask), sym/i ≤ 255 (u8)"
)]
fn decode_order_0(src: &mut &[u8], dst: &mut [u8]) -> Result<(), CramError> {
    // CramError::Truncated is materialised lazily on the err path; the
    // hot inner loop never builds one because `renormalize` returns
    // `Option<()>`. See the helper-module comment for the cost story.
    let truncated = || CramError::Truncated { context: "rans order-0 truncated" };
    let freq = read_frequencies_0(src).ok_or_else(truncated)?;
    let cum_freq = build_cumulative_frequencies(&freq);
    let sym_table = build_symbol_table(&cum_freq);
    let mut states = read_states(src).ok_or_else(truncated)?;

    for chunk in dst.chunks_mut(4) {
        for (d, state) in chunk.iter_mut().zip(states.iter_mut()) {
            let f = (*state & 0x0FFF) as u16;
            let sym = sym_table[f as usize];
            *d = sym;
            let i = sym as usize;
            // r[depends cram.codec.state_step_safety]
            *state = u32::from(freq[i])
                .wrapping_mul(*state >> 12)
                .wrapping_add(*state & 0x0FFF)
                .wrapping_sub(u32::from(cum_freq[i]));
            renormalize(state, src).ok_or_else(truncated)?;
        }
    }

    Ok(())
}

#[allow(
    clippy::indexing_slicing,
    reason = "indices are bounded: ctx/sym ≤ 255 (u8), f ≤ 4095 (12-bit mask)"
)]
fn decode_order_1_buf(
    src: &mut &[u8],
    dst: &mut [u8],
    buf: &mut Rans4x8Buf,
) -> Result<(), CramError> {
    // See decode_order_0 for the lazy-CramError pattern.
    let truncated = || CramError::Truncated { context: "rans order-1 truncated" };
    read_frequencies_1_into(src, &mut buf.freq).ok_or_else(truncated)?;

    for ctx in 0..ALPHABET_SIZE {
        buf.cum_freq[ctx] = build_cumulative_frequencies(&buf.freq[ctx]);
        buf.sym_tables[ctx] = build_symbol_table(&buf.cum_freq[ctx]);
    }

    let mut states = read_states(src).ok_or_else(truncated)?;
    let mut prev_syms = [0u8; 4];

    // Chunk-based interleaving: split output into 4 equal segments.
    // Each state decodes into its own segment sequentially.
    let chunk_size = dst.len() / 4;
    let bases: [usize; 4] = [0, chunk_size, chunk_size.wrapping_mul(2), chunk_size.wrapping_mul(3)];

    for pos in 0..chunk_size {
        for si in 0usize..4 {
            let mut state = states[si];
            let ctx = prev_syms[si] as usize;
            let f = (state & 0x0FFF) as u16;
            let sym = buf.sym_tables[ctx][f as usize];
            dst[bases[si].wrapping_add(pos)] = sym;
            let sym_idx = sym as usize;
            state = u32::from(buf.freq[ctx][sym_idx])
                .wrapping_mul(state >> 12)
                .wrapping_add(state & 0x0FFF)
                .wrapping_sub(u32::from(buf.cum_freq[ctx][sym_idx]));
            renormalize(&mut state, src).ok_or_else(truncated)?;
            states[si] = state;
            prev_syms[si] = sym;
        }
    }

    let remainder_start = chunk_size.wrapping_mul(4);
    for pos in remainder_start..dst.len() {
        let ctx = prev_syms[3] as usize;
        let f = (states[3] & 0x0FFF) as u16;
        let sym = buf.sym_tables[ctx][f as usize];
        dst[pos] = sym;
        let sym_idx = sym as usize;
        states[3] = u32::from(buf.freq[ctx][sym_idx])
            .wrapping_mul(states[3] >> 12)
            .wrapping_add(states[3] & 0x0FFF)
            .wrapping_sub(u32::from(buf.cum_freq[ctx][sym_idx]));
        if pos.checked_add(1).is_some_and(|next| next < dst.len()) {
            renormalize(&mut states[3], src).ok_or_else(truncated)?;
        }
        prev_syms[3] = sym;
    }

    Ok(())
}

/// Returns `None` when the source is exhausted mid-renormalize. The
/// caller materializes a `CramError::Truncated` from the `None` only
/// on the err path — see the helpers above for the size + drop story.
///
/// Mirrors htslib's `RansDecRenorm`: the common case for a well-formed
/// 4x8 stream is one renormalize byte per call (state was just above
/// `LOWER_BOUND >> 8` before the decode step). Two unrolled `if`-and-
/// read steps cover that and the slightly-rarer 2-byte case as
/// straight-line code, without a `while`-loop overhead per state
/// update. A small tail loop catches the pathological "state was
/// near 0" case (which only arises on malformed streams in practice
/// but is cheap to keep correct).
#[inline]
fn renormalize(state: &mut u32, src: &mut &[u8]) -> Option<()> {
    if *state >= LOWER_BOUND {
        return Some(());
    }
    *state = state.wrapping_shl(8) | u32::from(read_u8(src)?);
    if *state >= LOWER_BOUND {
        return Some(());
    }
    *state = state.wrapping_shl(8) | u32::from(read_u8(src)?);
    while *state < LOWER_BOUND {
        *state = state.wrapping_shl(8) | u32::from(read_u8(src)?);
    }
    Some(())
}

#[inline]
fn read_states(src: &mut &[u8]) -> Option<[u32; 4]> {
    let mut states = [0u32; 4];
    for s in &mut states {
        *s = read_u32_le(src)?;
    }
    Some(states)
}

#[allow(clippy::indexing_slicing, reason = "sym is u8, so sym as usize ≤ 255 < ALPHABET_SIZE=256")]
fn read_frequencies_0(src: &mut &[u8]) -> Option<[u16; ALPHABET_SIZE]> {
    let mut freq = [0u16; ALPHABET_SIZE];
    let mut sym = read_u8(src)?;
    let mut prev_sym = sym;

    loop {
        freq[sym as usize] = read_itf8_u16(src)?;
        sym = read_u8(src)?;

        if sym == 0 {
            break;
        }

        if sym == prev_sym.wrapping_add(1) {
            let run_len = read_u8(src)?;
            // r[impl cram.edge.rans_sym_overflow]
            // r[impl cram.codec.rans_sym_bounded+2]
            // sym can reach 255; break after writing freq[255] to avoid
            // wrapping to 0 and corrupting freq[0].
            for _ in 0..run_len {
                freq[sym as usize] = read_itf8_u16(src)?;
                if sym == 255 {
                    break;
                }
                sym = sym.wrapping_add(1);
            }
        }

        prev_sym = sym;
    }

    Some(freq)
}

#[allow(clippy::indexing_slicing, reason = "ctx is u8, so ctx as usize ≤ 255 < ALPHABET_SIZE=256")]
fn read_frequencies_1_into(
    src: &mut &[u8],
    freq: &mut [[u16; ALPHABET_SIZE]; ALPHABET_SIZE],
) -> Option<()> {
    let mut ctx = read_u8(src)?;
    let mut prev_ctx = ctx;

    loop {
        freq[ctx as usize] = read_frequencies_0(src)?;

        ctx = read_u8(src)?;
        if ctx == 0 {
            break;
        }

        if ctx == prev_ctx.wrapping_add(1) {
            let run_len = read_u8(src)?;
            for _ in 0..run_len {
                freq[ctx as usize] = read_frequencies_0(src)?;
                ctx = ctx.wrapping_add(1);
            }
        }

        prev_ctx = ctx;
    }

    Some(())
}

#[allow(
    clippy::indexing_slicing,
    reason = "i < ALPHABET_SIZE=256, sym ≤ 254 (guarded by < 255 check so sym+1 ≤ 255)"
)]
fn build_cumulative_frequencies(freq: &[u16; ALPHABET_SIZE]) -> [u16; ALPHABET_SIZE] {
    let mut cum = [0u16; ALPHABET_SIZE];
    let mut total = 0u16;
    for i in 0..ALPHABET_SIZE {
        cum[i] = total;
        total = total.wrapping_add(freq[i]);
    }
    cum
}

#[allow(
    clippy::indexing_slicing,
    reason = "sym < 255 (loop guard), so sym+1 ≤ 255 < ALPHABET_SIZE=256"
)]
fn build_symbol_table(cum_freq: &[u16; ALPHABET_SIZE]) -> [u8; 4096] {
    let mut table = [0u8; 4096];
    let mut sym = 0u8;

    for (f, g) in (0u16..).zip(&mut table) {
        while sym < 255
            && f >= cum_freq
                [(sym as usize).checked_add(1).expect("sym < 255 guarantees sym+1 ≤ 255")]
        {
            sym = sym.checked_add(1).expect("sym < 255 guarantees sym+1 ≤ 255");
        }
        *g = sym;
    }

    table
}

// These per-byte helpers return `Option<T>` rather than
// `Result<T, CramError>` for two compounding reasons:
//
// 1. Size. `CramError` is 80 bytes (sized to fit `Open { path: PathBuf,
//    source: std::io::Error }` and other heap-owning variants), so
//    `Result<u8, CramError>` is also 80 bytes — passed via memory/sret
//    in the non-inlined ABI. `Option<u8>` is 2 bytes, in a register.
// 2. Drop. `CramError` has a discriminant-dispatched Drop because some
//    variants own `PathBuf` / `SmolStr` / `io::Error`. samply flagged
//    `core::ptr::drop_in_place<CramError>` next to `read_u8` /
//    `renormalize` because the eager `ok_or(CramError::...)` form
//    constructed and then dropped a CramError on every successful read.
//
// Callers materialize a `CramError` only on the failure path, with
// `ok_or_else(|| CramError::Truncated { context: ... })` at the
// outermost function that knows the context.
//
// `split_first` / `split_first_chunk` (not `first` + `get(1..)`) is
// also load-bearing — the two-step form does two slice bounds checks
// and `slice::first` was at the top of the renormalize stack.
#[inline]
fn read_u8(src: &mut &[u8]) -> Option<u8> {
    let (&b, rest) = src.split_first()?;
    *src = rest;
    Some(b)
}

#[inline]
fn read_u32_le(src: &mut &[u8]) -> Option<u32> {
    let (head, rest) = src.split_first_chunk::<4>()?;
    *src = rest;
    Some(u32::from_le_bytes(*head))
}

#[inline]
fn read_itf8_u16(src: &mut &[u8]) -> Option<u16> {
    // ITF8 for frequency table values (they fit in u16)
    let val = super::varint::read_itf8_from(src)?;
    #[expect(
        clippy::cast_possible_truncation,
        reason = "rANS frequency table values are bounded by 4096 (12-bit), fits in u16"
    )]
    Some(val as u16)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn invalid_rans_order_returns_error() {
        // Build a minimal rANS 4x8 header with order=2 (invalid; only 0 and 1 are valid)
        let mut src = Vec::new();
        src.push(2u8); // order = 2 (invalid)
        src.extend_from_slice(&0u32.to_le_bytes()); // compressed_size
        src.extend_from_slice(&0u32.to_le_bytes()); // uncompressed_size

        let err = decode(&src).unwrap_err();
        assert!(matches!(err, CramError::InvalidRansOrder { order: 2 }));
    }

    #[test]
    fn invalid_rans_order_field_value() {
        let mut src = Vec::new();
        src.push(7u8); // order = 7 (invalid)
        src.extend_from_slice(&0u32.to_le_bytes());
        src.extend_from_slice(&0u32.to_le_bytes());

        let err = decode(&src).unwrap_err();
        assert!(matches!(err, CramError::InvalidRansOrder { order: 7 }));
    }

    // r[verify cram.codec.rans4x8]
    // r[verify cram.edge.rans_sym_overflow]
    #[test]
    fn decode_order_0_noodles_test_vector() {
        // Test vector from noodles rANS 4x8 tests: decodes to "noodles"
        let src = [
            0x00, 0x25, 0x00, 0x00, 0x00, 0x07, 0x00, 0x00, 0x00, 0x64, 0x82, 0x49, 0x65, 0x00,
            0x82, 0x49, 0x6c, 0x82, 0x49, 0x6e, 0x82, 0x49, 0x6f, 0x00, 0x84, 0x92, 0x73, 0x82,
            0x49, 0x00, 0xe2, 0x06, 0x83, 0x18, 0x74, 0x7b, 0x41, 0x0c, 0x2b, 0xa9, 0x41, 0x0c,
            0x25, 0x31, 0x80, 0x03,
        ];
        let result = decode(&src).unwrap();
        assert_eq!(result, b"noodles");
    }

    #[test]
    fn decode_order_1_noodles_test_vector() {
        // Test vector from noodles rANS 4x8 tests: decodes to "noodles"
        let src = [
            0x01, 0x3b, 0x00, 0x00, 0x00, 0x07, 0x00, 0x00, 0x00, 0x00, 0x64, 0x84, 0x00, 0x6e,
            0x84, 0x00, 0x6f, 0x00, 0x87, 0xff, 0x00, 0x64, 0x6c, 0x8f, 0xff, 0x00, 0x65, 0x00,
            0x73, 0x8f, 0xff, 0x00, 0x6c, 0x65, 0x8f, 0xff, 0x00, 0x6e, 0x6f, 0x8f, 0xff, 0x00,
            0x6f, 0x00, 0x64, 0x87, 0xff, 0x6f, 0x88, 0x00, 0x00, 0x00, 0x00, 0x04, 0x00, 0x02,
            0x02, 0x28, 0x00, 0x01, 0x02, 0x28, 0x00, 0x01, 0x02, 0x60, 0x00, 0x02,
        ];
        let result = decode(&src).unwrap();
        assert_eq!(result, b"noodles");
    }
}
