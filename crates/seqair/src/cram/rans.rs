//! rANS 4x8 codec (CRAM compression method 4).
//!
//! 4-way interleaved asymmetric numeral systems with 8-bit renormalization.
//! Supports order-0 and order-1 modes.

use super::reader::CramError;

const ALPHABET_SIZE: usize = 256;
const LOWER_BOUND: u32 = 1 << 23;

// r[impl cram.codec.rans4x8]
/// Decode a rANS 4x8 compressed block.
pub fn decode(src: &[u8]) -> Result<Vec<u8>, CramError> {
    let mut cur: &[u8] = src;

    // Header: order (u8), compressed_size (u32 LE), uncompressed_size (u32 LE)
    let order = read_u8(&mut cur)?;
    let _compressed_size = read_u32_le(&mut cur)?;
    let uncompressed_size = read_u32_le(&mut cur)? as usize;

    super::reader::check_alloc_size(uncompressed_size, "rANS 4x8 output")?;
    let mut dst = vec![0u8; uncompressed_size];

    match order {
        0 => decode_order_0(&mut cur, &mut dst)?,
        1 => decode_order_1(&mut cur, &mut dst)?,
        _ => return Err(CramError::InvalidRansOrder { order }),
    }

    Ok(dst)
}

#[allow(
    clippy::indexing_slicing,
    reason = "indices are bounded: f ≤ 4095 (12-bit mask), sym/i ≤ 255 (u8)"
)]
fn decode_order_0(src: &mut &[u8], dst: &mut [u8]) -> Result<(), CramError> {
    let freq = read_frequencies_0(src)?;
    let cum_freq = build_cumulative_frequencies(&freq);
    let sym_table = build_symbol_table(&cum_freq);
    let mut states = read_states(src)?;

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
            renormalize(state, src)?;
        }
    }

    Ok(())
}

#[allow(
    clippy::indexing_slicing,
    reason = "indices are bounded: ctx/sym ≤ 255 (u8), f ≤ 4095 (12-bit mask)"
)]
fn decode_order_1(src: &mut &[u8], dst: &mut [u8]) -> Result<(), CramError> {
    let freq = read_frequencies_1(src)?;

    let mut cum_freq = vec![[0u16; ALPHABET_SIZE]; ALPHABET_SIZE];
    let mut sym_tables: Box<[[u8; 4096]; ALPHABET_SIZE]> = Box::new([[0u8; 4096]; ALPHABET_SIZE]);
    for ctx in 0..ALPHABET_SIZE {
        cum_freq[ctx] = build_cumulative_frequencies(&freq[ctx]);
        sym_tables[ctx] = build_symbol_table(&cum_freq[ctx]);
    }

    let mut states = read_states(src)?;
    let mut prev_syms = [0u8; 4];

    // Chunk-based interleaving: split output into 4 equal segments.
    // Each state decodes into its own segment sequentially.
    let chunk_size = dst.len() / 4;

    for pos in 0..chunk_size {
        for si in 0usize..4 {
            let out_idx = si
                .checked_mul(chunk_size)
                .and_then(|base| base.checked_add(pos))
                .ok_or(CramError::Truncated { context: "rans order-1 out_idx overflow" })?;
            let ctx = prev_syms[si] as usize;
            let f = (states[si] & 0x0FFF) as u16;
            let sym = sym_tables[ctx][f as usize];
            dst[out_idx] = sym;
            let sym_idx = sym as usize;
            states[si] = u32::from(freq[ctx][sym_idx])
                .wrapping_mul(states[si] >> 12)
                .wrapping_add(states[si] & 0x0FFF)
                .wrapping_sub(u32::from(cum_freq[ctx][sym_idx]));
            renormalize(&mut states[si], src)?;
            prev_syms[si] = sym;
        }
    }

    // Remainder: state 3 handles trailing bytes
    let remainder_start = 4usize
        .checked_mul(chunk_size)
        .ok_or(CramError::Truncated { context: "rans order-1 remainder offset overflow" })?;
    for pos in remainder_start..dst.len() {
        let ctx = prev_syms[3] as usize;
        let f = (states[3] & 0x0FFF) as u16;
        let sym = sym_tables[ctx][f as usize];
        dst[pos] = sym;
        let sym_idx = sym as usize;
        states[3] = u32::from(freq[ctx][sym_idx])
            .wrapping_mul(states[3] >> 12)
            .wrapping_add(states[3] & 0x0FFF)
            .wrapping_sub(u32::from(cum_freq[ctx][sym_idx]));
        if pos.checked_add(1).is_some_and(|next| next < dst.len()) {
            renormalize(&mut states[3], src)?;
        }
        prev_syms[3] = sym;
    }

    Ok(())
}

fn renormalize(state: &mut u32, src: &mut &[u8]) -> Result<(), CramError> {
    while *state < LOWER_BOUND {
        let b = u32::from(read_u8(src)?);
        *state = (*state << 8) | b;
    }
    Ok(())
}

fn read_states(src: &mut &[u8]) -> Result<[u32; 4], CramError> {
    let mut states = [0u32; 4];
    for s in &mut states {
        *s = read_u32_le(src)?;
    }
    Ok(states)
}

#[allow(clippy::indexing_slicing, reason = "sym is u8, so sym as usize ≤ 255 < ALPHABET_SIZE=256")]
fn read_frequencies_0(src: &mut &[u8]) -> Result<[u16; ALPHABET_SIZE], CramError> {
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

    Ok(freq)
}

#[allow(clippy::indexing_slicing, reason = "ctx is u8, so ctx as usize ≤ 255 < ALPHABET_SIZE=256")]
fn read_frequencies_1(
    src: &mut &[u8],
) -> Result<Box<[[u16; ALPHABET_SIZE]; ALPHABET_SIZE]>, CramError> {
    let mut freq = Box::new([[0u16; ALPHABET_SIZE]; ALPHABET_SIZE]);

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

    Ok(freq)
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

fn read_u8(src: &mut &[u8]) -> Result<u8, CramError> {
    let &b = src.first().ok_or(CramError::Truncated { context: "rans u8" })?;
    *src = src.get(1..).ok_or(CramError::Truncated { context: "rans u8" })?;
    Ok(b)
}

fn read_u32_le(src: &mut &[u8]) -> Result<u32, CramError> {
    let bytes = src.get(..4).ok_or(CramError::Truncated { context: "rans u32" })?;
    let val = u32::from_le_bytes(
        bytes.try_into().map_err(|_| CramError::Truncated { context: "rans u32" })?,
    );
    *src = src.get(4..).ok_or(CramError::Truncated { context: "rans u32" })?;
    Ok(val)
}

fn read_itf8_u16(src: &mut &[u8]) -> Result<u16, CramError> {
    // ITF8 for frequency table values (they fit in u16)
    let val = super::varint::read_itf8_from(src)
        .ok_or(CramError::Truncated { context: "rans frequency itf8" })?;
    #[expect(
        clippy::cast_possible_truncation,
        reason = "rANS frequency table values are bounded by 4096 (12-bit), fits in u16"
    )]
    Ok(val as u16)
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
