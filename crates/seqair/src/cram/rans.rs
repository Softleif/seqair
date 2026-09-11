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
/// without a reusable buffer). Prefer `decode_with_buf` on the hot path.
pub fn decode(src: &[u8]) -> Result<Vec<u8>, CramError> {
    let mut buf = Rans4x8Buf::new();
    decode_with_buf(src, &mut buf)
}

/// Decode a rANS 4x8 compressed block, reusing the caller-owned
/// [`Rans4x8Buf`] for order-1 tables. Order-0 blocks don't touch the buf.
pub(crate) fn decode_with_buf(src: &[u8], buf: &mut Rans4x8Buf) -> Result<Vec<u8>, CramError> {
    let mut cur: &[u8] = src;

    // Header: order (u8), compressed_size (u32 LE), uncompressed_size (u32 LE).
    // The 9-byte header is followed by `compressed_size` payload bytes; htslib
    // rejects mismatch with `in_sz != in_size - 9` (rANS_static.c:245).
    // r[impl cram.codec.rans4x8_compressed_size_check]
    let order = read_u8(&mut cur).ok_or_else(|| CramError::Truncated { context: "rans header" })?;
    let compressed_size =
        read_u32_le(&mut cur).ok_or_else(|| CramError::Truncated { context: "rans header" })?;
    let uncompressed_size = read_u32_le(&mut cur)
        .ok_or_else(|| CramError::Truncated { context: "rans header" })?
        as usize;

    if cur.len() != compressed_size as usize {
        return Err(CramError::Rans4x8CompressedSizeMismatch {
            advertised: compressed_size,
            actual: cur.len(),
        });
    }

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
    let freq = read_frequencies_0(src)?;
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
    read_frequencies_1_into(src, &mut buf.freq)?;

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

// r[impl cram.codec.alphabet_run_bounded]
#[allow(clippy::indexing_slicing, reason = "sym is u8, so sym as usize ≤ 255 < ALPHABET_SIZE=256")]
fn read_frequencies_0(src: &mut &[u8]) -> Result<[u16; ALPHABET_SIZE], CramError> {
    let truncated = || CramError::Truncated { context: "rans 4x8 frequencies" };
    let mut freq = [0u16; ALPHABET_SIZE];
    let mut sym = read_u8(src).ok_or_else(truncated)?;
    let mut prev_sym = sym;

    loop {
        freq[sym as usize] = read_itf8_u16(src).ok_or_else(truncated)?;
        sym = read_u8(src).ok_or_else(truncated)?;

        if sym == 0 {
            break;
        }

        if sym == prev_sym.wrapping_add(1) {
            let run_len = read_u8(src).ok_or_else(truncated)?;
            // After the inner loop sym becomes start+run_len. htscodecs
            // rejects any run where the post-increment value wraps past
            // 255 (see `rANS_static.c::decode_freq` `if (j > 255) goto
            // cleanup`). Tolerating wraparound silently desyncs the source
            // stream — the next byte after the run would be misread as a
            // freq value.
            let end = u32::from(sym).wrapping_add(u32::from(run_len));
            if end > 255 {
                return Err(CramError::MalformedAlphabetRun { start: sym, len: run_len });
            }
            for _ in 0..run_len {
                freq[sym as usize] = read_itf8_u16(src).ok_or_else(truncated)?;
                sym = sym.wrapping_add(1);
            }
        }

        prev_sym = sym;
    }

    Ok(freq)
}

#[allow(clippy::indexing_slicing, reason = "ctx is u8, so ctx as usize ≤ 255 < ALPHABET_SIZE=256")]
fn read_frequencies_1_into(
    src: &mut &[u8],
    freq: &mut [[u16; ALPHABET_SIZE]; ALPHABET_SIZE],
) -> Result<(), CramError> {
    let truncated = || CramError::Truncated { context: "rans 4x8 order-1 frequencies" };
    let mut ctx = read_u8(src).ok_or_else(truncated)?;
    let mut prev_ctx = ctx;

    loop {
        freq[ctx as usize] = read_frequencies_0(src)?;

        ctx = read_u8(src).ok_or_else(truncated)?;
        if ctx == 0 {
            break;
        }

        if ctx == prev_ctx.wrapping_add(1) {
            let run_len = read_u8(src).ok_or_else(truncated)?;
            // Same wraparound concern as read_frequencies_0 but for the
            // outer per-context dimension.
            let end = u32::from(ctx).wrapping_add(u32::from(run_len));
            if end > 255 {
                return Err(CramError::MalformedAlphabetRun { start: ctx, len: run_len });
            }
            for _ in 0..run_len {
                freq[ctx as usize] = read_frequencies_0(src)?;
                ctx = ctx.wrapping_add(1);
            }
        }

        prev_ctx = ctx;
    }

    Ok(())
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

// Per-byte primitives (`read_u8`, `read_u32_le`) live in `super::codec_io`
// — see that module's docs for the `Option<T>` design rationale (avoiding
// `drop_in_place<CramError>` on the per-byte hot path).
use super::codec_io::{read_u8, read_u32_le};

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
#[allow(
    clippy::arithmetic_side_effects,
    clippy::needless_range_loop,
    reason = "test code: bounded values, range-by-index is clearer than enumerate here"
)]
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

    // r[verify cram.codec.rans4x8_compressed_size_check]
    #[test]
    fn rans_compressed_size_mismatch_returns_error() {
        // Header advertises compressed_size = 100, but the actual remaining
        // payload after the 9-byte header is only ~16 bytes. htslib rejects
        // exactly this with `if (in_sz != in_size-9) return NULL`
        // (rANS_static.c:245). Catching the lie here instead of getting a
        // confusing downstream Truncated error.
        let mut src = Vec::new();
        src.push(0u8); // order = 0
        src.extend_from_slice(&100u32.to_le_bytes()); // compressed_size = 100 (lie)
        src.extend_from_slice(&7u32.to_le_bytes()); // uncompressed_size = 7
        // Only ~16 bytes of payload, far short of 100.
        src.extend_from_slice(&[0u8; 16]);

        let err = decode(&src).unwrap_err();
        assert!(
            matches!(err, CramError::Rans4x8CompressedSizeMismatch { advertised: 100, actual: 16 }),
            "expected Rans4x8CompressedSizeMismatch, got: {err:?}",
        );
    }

    // r[verify cram.codec.alphabet_run_bounded]
    #[test]
    fn read_frequencies_0_run_overflow_returns_error() {
        // sym=200, freq=ITF8(1)=0x01, then sym=201 (== prev+1) triggers a
        // run of length 100 — this would extend past sym 255. With ITF8(1)
        // payload bytes in between, the well-formed-prefix portion of the
        // stream is: [200, 0x01, 201, 100, ...freq values...]. The fix
        // must catch the bad `start + len` BEFORE consuming the run's
        // freq values.
        let mut src: &[u8] = &[200, 0x01, 201, 100];
        let err = read_frequencies_0(&mut src).unwrap_err();
        assert!(
            matches!(err, CramError::MalformedAlphabetRun { start: 201, len: 100 }),
            "expected MalformedAlphabetRun, got: {err:?}"
        );
    }

    #[test]
    fn read_frequencies_0_run_exactly_to_255_ok() {
        // Canonical htscodecs encoding for F[200..=255]=1: emit 200 (no
        // run), F[200]=1, then 201 with run-len 54 (covers F[201..254]),
        // then F[255] in the next outer iteration.
        // Stream: [200, 1, 201, 54, 1*54, 1, 0] = 60 bytes.
        let mut src: Vec<u8> = vec![200, 0x01, 201, 54];
        src.extend(std::iter::repeat_n(0x01u8, 54)); // F[201..254]
        src.push(0x01); // F[255]
        src.push(0); // terminator
        let mut cur: &[u8] = &src;
        let freq = read_frequencies_0(&mut cur).unwrap();
        for s in 200..=255usize {
            assert_eq!(freq[s], 1, "freq[{s}] should be 1");
        }
        assert_eq!(freq[199], 0, "freq[199] should be 0");
        assert_eq!(freq[0], 0, "freq[0] should be 0 (no wraparound)");
    }

    proptest::proptest! {
        // Mirror of rans_nx16's bounds proptest. Only tests the rejection
        // path — constructing a fully valid 4x8 freq stream needs an extra
        // freq byte after the run for the trailing sym, which is fiddly to
        // get right at the boundary; the valid path is covered by the fixed
        // test above and integration round-trips.
        #[test]
        fn read_frequencies_0_invalid_run_rejected(
            start in 1u8..=255,
            len in 0u8..=255,
        ) {
            proptest::prop_assume!(u32::from(start) + u32::from(len) > 255);
            let prev = start.checked_sub(1).expect("start >= 1");
            // Provide enough freq bytes that truncation can't fire before the
            // bounds check: prev + start are 2 reads; the run-len byte is 1
            // read; we don't need run payload because the check fires first.
            let stream = [prev, 0x01, start, len];
            let mut cur: &[u8] = &stream;
            let err = read_frequencies_0(&mut cur).unwrap_err();
            proptest::prop_assert!(
                matches!(err, CramError::MalformedAlphabetRun { start: s, len: l }
                             if s == start && l == len),
                "expected MalformedAlphabetRun for start={start}, len={len}, got: {err:?}",
            );
        }
    }

    // r[verify cram.codec.alphabet_run_bounded]
    #[test]
    fn read_frequencies_1_ctx_run_overflow_returns_error() {
        // Outer ctx loop has the same wraparound concern. Stream: ctx=200,
        // (entire sub-table for ctx=200), ctx=201 (== prev+1), run_len=100.
        // Bound check 201+100=301 > 255 must fire BEFORE the inner sub-table
        // reads. We use a minimal single-symbol sub-table for ctx=200:
        // [sym=0, freq=ITF8(4096), terminator=0].
        let src: Vec<u8> = vec![200, 0, 0x80, 0x20, 0, 201, 100];
        let mut freq = Box::new([[0u16; ALPHABET_SIZE]; ALPHABET_SIZE]);
        let mut cur: &[u8] = &src;
        let err = read_frequencies_1_into(&mut cur, &mut freq).unwrap_err();
        assert!(
            matches!(err, CramError::MalformedAlphabetRun { start: 201, len: 100 }),
            "expected MalformedAlphabetRun, got: {err:?}"
        );
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
