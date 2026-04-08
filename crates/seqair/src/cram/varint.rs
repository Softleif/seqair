//! Encode and decode CRAM variable-length integers. [`decode_itf8`] and [`decode_ltf8`]
//! parse the ITF-8 and LTF-8 formats used in container/slice headers throughout the CRAM spec.

use std::io::Read;

// r[impl cram.itf8]
/// Decode an ITF8 value from a byte slice.
///
/// Returns `(value, bytes_consumed)`. The value is unsigned; use `as i32`
/// for fields that can be negative (`ref_seq_id`, `alignment_start`).
pub fn decode_itf8(buf: &[u8]) -> Option<(u32, usize)> {
    let &b0 = buf.first()?;

    if b0 & 0x80 == 0 {
        // 0xxxxxxx — 1 byte, 7 bits
        Some((u32::from(b0), 1))
    } else if b0 & 0x40 == 0 {
        // 10xxxxxx — 2 bytes, 14 bits
        let &b1 = buf.get(1)?;
        let val = (u32::from(b0 & 0x3F) << 8) | u32::from(b1);
        Some((val, 2))
    } else if b0 & 0x20 == 0 {
        // 110xxxxx — 3 bytes, 21 bits
        let &b1 = buf.get(1)?;
        let &b2 = buf.get(2)?;
        let val = (u32::from(b0 & 0x1F) << 16) | (u32::from(b1) << 8) | u32::from(b2);
        Some((val, 3))
    } else if b0 & 0x10 == 0 {
        // 1110xxxx — 4 bytes, 28 bits
        let &b1 = buf.get(1)?;
        let &b2 = buf.get(2)?;
        let &b3 = buf.get(3)?;
        let val = (u32::from(b0 & 0x0F) << 24)
            | (u32::from(b1) << 16)
            | (u32::from(b2) << 8)
            | u32::from(b3);
        Some((val, 4))
    } else {
        // 1111xxxx — 5 bytes, 32 bits
        let &b1 = buf.get(1)?;
        let &b2 = buf.get(2)?;
        let &b3 = buf.get(3)?;
        let &b4 = buf.get(4)?;
        let val = (u32::from(b0 & 0x0F) << 28)
            | (u32::from(b1) << 20)
            | (u32::from(b2) << 12)
            | (u32::from(b3) << 4)
            | u32::from(b4 & 0x0F);
        Some((val, 5))
    }
}

/// Read an ITF8 value from a `Read` source.
pub fn read_itf8<R: Read>(r: &mut R) -> std::io::Result<u32> {
    let mut b0 = [0u8; 1];
    r.read_exact(&mut b0)?;
    let b0 = b0[0];

    if b0 & 0x80 == 0 {
        Ok(u32::from(b0))
    } else if b0 & 0x40 == 0 {
        let mut buf = [0u8; 1];
        r.read_exact(&mut buf)?;
        Ok((u32::from(b0 & 0x3F) << 8) | u32::from(buf[0]))
    } else if b0 & 0x20 == 0 {
        let mut buf = [0u8; 2];
        r.read_exact(&mut buf)?;
        Ok((u32::from(b0 & 0x1F) << 16) | (u32::from(buf[0]) << 8) | u32::from(buf[1]))
    } else if b0 & 0x10 == 0 {
        let mut buf = [0u8; 3];
        r.read_exact(&mut buf)?;
        Ok((u32::from(b0 & 0x0F) << 24)
            | (u32::from(buf[0]) << 16)
            | (u32::from(buf[1]) << 8)
            | u32::from(buf[2]))
    } else {
        let mut buf = [0u8; 4];
        r.read_exact(&mut buf)?;
        Ok((u32::from(b0 & 0x0F) << 28)
            | (u32::from(buf[0]) << 20)
            | (u32::from(buf[1]) << 12)
            | (u32::from(buf[2]) << 4)
            | u32::from(buf[3] & 0x0F))
    }
}

/// Read an ITF8 value from a cursor (advancing it).
pub fn read_itf8_from(cursor: &mut &[u8]) -> Option<u32> {
    let (val, consumed) = decode_itf8(cursor)?;
    *cursor = cursor.get(consumed..)?;
    Some(val)
}

// r[impl cram.ltf8]
/// Decode an LTF8 value from a byte slice.
///
/// Returns `(value, bytes_consumed)`.
pub fn decode_ltf8(buf: &[u8]) -> Option<(u64, usize)> {
    let &b0 = buf.first()?;

    if b0 & 0x80 == 0 {
        Some((u64::from(b0), 1))
    } else if b0 & 0x40 == 0 {
        let &b1 = buf.get(1)?;
        Some(((u64::from(b0 & 0x3F) << 8) | u64::from(b1), 2))
    } else if b0 & 0x20 == 0 {
        let &b1 = buf.get(1)?;
        let &b2 = buf.get(2)?;
        Some(((u64::from(b0 & 0x1F) << 16) | (u64::from(b1) << 8) | u64::from(b2), 3))
    } else if b0 & 0x10 == 0 {
        let &b1 = buf.get(1)?;
        let &b2 = buf.get(2)?;
        let &b3 = buf.get(3)?;
        Some((
            (u64::from(b0 & 0x0F) << 24)
                | (u64::from(b1) << 16)
                | (u64::from(b2) << 8)
                | u64::from(b3),
            4,
        ))
    } else if b0 & 0x08 == 0 {
        let &b1 = buf.get(1)?;
        let &b2 = buf.get(2)?;
        let &b3 = buf.get(3)?;
        let &b4 = buf.get(4)?;
        Some((
            (u64::from(b0 & 0x07) << 32)
                | (u64::from(b1) << 24)
                | (u64::from(b2) << 16)
                | (u64::from(b3) << 8)
                | u64::from(b4),
            5,
        ))
    } else if b0 & 0x04 == 0 {
        let &b1 = buf.get(1)?;
        let &b2 = buf.get(2)?;
        let &b3 = buf.get(3)?;
        let &b4 = buf.get(4)?;
        let &b5 = buf.get(5)?;
        Some((
            (u64::from(b0 & 0x03) << 40)
                | (u64::from(b1) << 32)
                | (u64::from(b2) << 24)
                | (u64::from(b3) << 16)
                | (u64::from(b4) << 8)
                | u64::from(b5),
            6,
        ))
    } else if b0 & 0x02 == 0 {
        let &b1 = buf.get(1)?;
        let &b2 = buf.get(2)?;
        let &b3 = buf.get(3)?;
        let &b4 = buf.get(4)?;
        let &b5 = buf.get(5)?;
        let &b6 = buf.get(6)?;
        Some((
            (u64::from(b0 & 0x01) << 48)
                | (u64::from(b1) << 40)
                | (u64::from(b2) << 32)
                | (u64::from(b3) << 24)
                | (u64::from(b4) << 16)
                | (u64::from(b5) << 8)
                | u64::from(b6),
            7,
        ))
    } else if b0 & 0x01 == 0 {
        let &b1 = buf.get(1)?;
        let &b2 = buf.get(2)?;
        let &b3 = buf.get(3)?;
        let &b4 = buf.get(4)?;
        let &b5 = buf.get(5)?;
        let &b6 = buf.get(6)?;
        let &b7 = buf.get(7)?;
        Some((
            (u64::from(b1) << 48)
                | (u64::from(b2) << 40)
                | (u64::from(b3) << 32)
                | (u64::from(b4) << 24)
                | (u64::from(b5) << 16)
                | (u64::from(b6) << 8)
                | u64::from(b7),
            8,
        ))
    } else {
        // 0xFF prefix — full 8 continuation bytes
        let &b1 = buf.get(1)?;
        let &b2 = buf.get(2)?;
        let &b3 = buf.get(3)?;
        let &b4 = buf.get(4)?;
        let &b5 = buf.get(5)?;
        let &b6 = buf.get(6)?;
        let &b7 = buf.get(7)?;
        let &b8 = buf.get(8)?;
        Some((
            (u64::from(b1) << 56)
                | (u64::from(b2) << 48)
                | (u64::from(b3) << 40)
                | (u64::from(b4) << 32)
                | (u64::from(b5) << 24)
                | (u64::from(b6) << 16)
                | (u64::from(b7) << 8)
                | u64::from(b8),
            9,
        ))
    }
}

/// Read an LTF8 value from a `Read` source.
pub fn read_ltf8<R: Read>(r: &mut R) -> std::io::Result<u64> {
    let mut b0 = [0u8; 1];
    r.read_exact(&mut b0)?;
    let b0 = b0[0];

    let (n_extra, mask, shift) = match b0.leading_ones() {
        0 => return Ok(u64::from(b0)),
        1 => (1, 0x3Fu8, 8),
        2 => (2, 0x1Fu8, 16),
        3 => (3, 0x0Fu8, 24),
        4 => (4, 0x07u8, 32),
        5 => (5, 0x03u8, 40),
        6 => (6, 0x01u8, 48),
        7 => (7, 0x00u8, 48), // top byte is prefix-only
        _ => (8, 0x00u8, 56), // 0xFF: 8 continuation bytes
    };

    let mut extra = [0u8; 8];
    debug_assert!(n_extra <= 8, "n_extra out of bounds: {n_extra} > 8");
    #[allow(clippy::indexing_slicing, reason = "n_extra ≤ 8")]
    r.read_exact(&mut extra[..n_extra])?;

    let mut val = u64::from(b0 & mask) << shift;
    for (i, &byte) in extra.iter().take(n_extra).enumerate() {
        // i < n_extra ≤ 8, so n_extra-1-i ≥ 0 and the product ≤ 56
        let bit_shift = n_extra
            .checked_sub(1)
            .and_then(|v| v.checked_sub(i))
            .and_then(|v| v.checked_mul(8))
            .expect("n_extra ≤ 8 and i < n_extra, so (n_extra-1-i)*8 ≤ 56");
        val |= u64::from(byte) << bit_shift;
    }
    Ok(val)
}

/// Read an LTF8 value from a cursor (advancing it).
pub fn read_ltf8_from(cursor: &mut &[u8]) -> Option<u64> {
    let (val, consumed) = decode_ltf8(cursor)?;
    *cursor = cursor.get(consumed..)?;
    Some(val)
}

/// Encode an ITF8 value into a buffer.
///
/// Returns the number of bytes written (1-5).
/// Used only for testing roundtrips.
#[cfg(test)]
#[allow(
    clippy::cast_possible_truncation,
    reason = "intentional byte extraction via bitmask shifts in ITF8 encoding"
)]
pub(super) fn encode_itf8(val: u32, buf: &mut [u8; 5]) -> usize {
    if val < 0x80 {
        buf[0] = val as u8;
        1
    } else if val < 0x4000 {
        buf[0] = 0x80 | ((val >> 8) as u8 & 0x3F);
        buf[1] = val as u8;
        2
    } else if val < 0x20_0000 {
        buf[0] = 0xC0 | ((val >> 16) as u8 & 0x1F);
        buf[1] = (val >> 8) as u8;
        buf[2] = val as u8;
        3
    } else if val < 0x1000_0000 {
        buf[0] = 0xE0 | ((val >> 24) as u8 & 0x0F);
        buf[1] = (val >> 16) as u8;
        buf[2] = (val >> 8) as u8;
        buf[3] = val as u8;
        4
    } else {
        buf[0] = 0xF0 | ((val >> 28) as u8 & 0x0F);
        buf[1] = (val >> 20) as u8;
        buf[2] = (val >> 12) as u8;
        buf[3] = (val >> 4) as u8;
        buf[4] = (val & 0x0F) as u8;
        5
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    // ── ITF8 known values ──────────────────────────────────────────────

    // r[verify cram.itf8]
    #[test]
    fn itf8_single_byte() {
        // 0xxxxxxx: values 0..127
        assert_eq!(decode_itf8(&[0x00]), Some((0, 1)));
        assert_eq!(decode_itf8(&[0x01]), Some((1, 1)));
        assert_eq!(decode_itf8(&[0x7F]), Some((127, 1)));
    }

    #[test]
    fn itf8_two_bytes() {
        // 10xxxxxx: values 128..16383
        assert_eq!(decode_itf8(&[0x80, 0x80]), Some((128, 2)));
        assert_eq!(decode_itf8(&[0xBF, 0xFF]), Some((0x3FFF, 2)));
    }

    #[test]
    fn itf8_three_bytes() {
        // 110xxxxx: values up to 2^21-1
        assert_eq!(decode_itf8(&[0xC0, 0x40, 0x00]), Some((0x4000, 3)));
        assert_eq!(decode_itf8(&[0xDF, 0xFF, 0xFF]), Some((0x1FFFFF, 3)));
    }

    #[test]
    fn itf8_four_bytes() {
        // 1110xxxx: values up to 2^28-1
        assert_eq!(decode_itf8(&[0xE0, 0x20, 0x00, 0x00]), Some((0x20_0000, 4)));
        assert_eq!(decode_itf8(&[0xEF, 0xFF, 0xFF, 0xFF]), Some((0x0FFF_FFFF, 4)));
    }

    #[test]
    fn itf8_five_bytes() {
        // 1111xxxx: values up to 2^32-1
        assert_eq!(decode_itf8(&[0xF1, 0x00, 0x00, 0x00, 0x00]), Some((0x1000_0000, 5)));
        assert_eq!(decode_itf8(&[0xFF, 0xFF, 0xFF, 0xFF, 0x0F]), Some((0xFFFF_FFFF, 5)));
    }

    #[test]
    fn itf8_negative_as_i32() {
        // -1 in two's complement = 0xFFFF_FFFF
        let (val, _) = decode_itf8(&[0xFF, 0xFF, 0xFF, 0xFF, 0x0F]).unwrap();
        assert_eq!(val.cast_signed(), -1i32);

        // -2 = 0xFFFF_FFFE
        let (val, _) = decode_itf8(&[0xFF, 0xFF, 0xFF, 0xFF, 0x0E]).unwrap();
        assert_eq!(val.cast_signed(), -2i32);
    }

    #[test]
    fn itf8_empty_buffer() {
        assert_eq!(decode_itf8(&[]), None);
    }

    #[test]
    fn itf8_truncated() {
        // 2-byte encoding but only 1 byte available
        assert_eq!(decode_itf8(&[0x80]), None);
        // 5-byte encoding but only 3 bytes available
        assert_eq!(decode_itf8(&[0xF0, 0x00, 0x00]), None);
    }

    #[test]
    fn itf8_read_from_cursor() {
        let data = [0x80, 0x80, 0x05]; // ITF8(128) followed by ITF8(5)
        let mut cursor: &[u8] = &data;
        assert_eq!(read_itf8_from(&mut cursor), Some(128));
        assert_eq!(read_itf8_from(&mut cursor), Some(5));
        assert_eq!(read_itf8_from(&mut cursor), None); // exhausted
    }

    // ── ITF8 proptest roundtrip ────────────────────────────────────────

    proptest! {
        // Verify that encode_itf8 produces wire bytes that match the ITF8 spec bit patterns,
        // and that decode_itf8 recovers the original value.
        #[test]
        #[allow(clippy::cast_possible_truncation, reason = "val < 0x80 is checked before the cast")]
        fn itf8_wire_format_matches_spec(val: u32) {
            let mut buf = [0u8; 5];
            let n = encode_itf8(val, &mut buf);

            if val < 0x80 {
                prop_assert_eq!(n, 1);
                prop_assert_eq!(buf[0], val as u8, "1-byte: high bit must be 0");
                prop_assert_eq!(buf[0] & 0x80, 0);
            } else if val < 0x4000 {
                prop_assert_eq!(n, 2);
                prop_assert_eq!(buf[0] & 0xC0, 0x80, "2-byte: top 2 bits must be 10");
                let recovered = (u32::from(buf[0] & 0x3F) << 8) | u32::from(buf[1]);
                prop_assert_eq!(recovered, val);
            } else if val < 0x20_0000 {
                prop_assert_eq!(n, 3);
                prop_assert_eq!(buf[0] & 0xE0, 0xC0, "3-byte: top 3 bits must be 110");
                let recovered = (u32::from(buf[0] & 0x1F) << 16)
                    | (u32::from(buf[1]) << 8)
                    | u32::from(buf[2]);
                prop_assert_eq!(recovered, val);
            } else if val < 0x1000_0000 {
                prop_assert_eq!(n, 4);
                prop_assert_eq!(buf[0] & 0xF0, 0xE0, "4-byte: top 4 bits must be 1110");
                let recovered = (u32::from(buf[0] & 0x0F) << 24)
                    | (u32::from(buf[1]) << 16)
                    | (u32::from(buf[2]) << 8)
                    | u32::from(buf[3]);
                prop_assert_eq!(recovered, val);
            } else {
                prop_assert_eq!(n, 5);
                prop_assert_eq!(buf[0] & 0xF0, 0xF0, "5-byte: top 4 bits must be 1111");
                let recovered = (u32::from(buf[0] & 0x0F) << 28)
                    | (u32::from(buf[1]) << 20)
                    | (u32::from(buf[2]) << 12)
                    | (u32::from(buf[3]) << 4)
                    | u32::from(buf[4] & 0x0F);
                prop_assert_eq!(recovered, val);
            }

            // decode_itf8 must also recover the original value
            let (decoded, consumed) = decode_itf8(&buf).unwrap();
            prop_assert_eq!(decoded, val);
            prop_assert_eq!(consumed, n);
        }

        // All three decode paths must agree on the same encoded bytes.
        #[test]
        fn itf8_all_decoders_agree(val: u32) {
            let mut buf = [0u8; 5];
            let n = encode_itf8(val, &mut buf);
            let encoded = buf.get(..n).unwrap();

            let (slice_val, _) = decode_itf8(encoded).unwrap();
            let stream_val = read_itf8(&mut &*encoded).unwrap();
            let mut cursor: &[u8] = encoded;
            let cursor_val = read_itf8_from(&mut cursor).unwrap();

            prop_assert_eq!(slice_val, val);
            prop_assert_eq!(stream_val, val);
            prop_assert_eq!(cursor_val, val);
        }
    }

    // ── LTF8 known values ──────────────────────────────────────────────

    // r[verify cram.ltf8]
    #[test]
    fn ltf8_single_byte() {
        assert_eq!(decode_ltf8(&[0x00]), Some((0, 1)));
        assert_eq!(decode_ltf8(&[0x7F]), Some((127, 1)));
    }

    #[test]
    fn ltf8_two_bytes() {
        assert_eq!(decode_ltf8(&[0x80, 0x80]), Some((128, 2)));
    }

    #[test]
    fn ltf8_nine_bytes() {
        // 0xFF prefix, then 8 bytes for the full 64-bit value
        assert_eq!(
            decode_ltf8(&[0xFF, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01]),
            Some((1, 9))
        );
        assert_eq!(
            decode_ltf8(&[0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF]),
            Some((u64::MAX, 9))
        );
    }

    #[test]
    fn ltf8_empty_buffer() {
        assert_eq!(decode_ltf8(&[]), None);
    }

    #[test]
    fn ltf8_read_from_cursor() {
        let data = [0x80, 0x80, 0x05]; // LTF8(128) followed by LTF8(5)
        let mut cursor: &[u8] = &data;
        assert_eq!(read_ltf8_from(&mut cursor), Some(128));
        assert_eq!(read_ltf8_from(&mut cursor), Some(5));
        assert_eq!(read_ltf8_from(&mut cursor), None);
    }

    // ── LTF8 proptest ──────────────────────────────────────────────────

    proptest! {
        #[test]
        #[allow(clippy::cast_possible_truncation, reason = "val in 0..0x80 is checked by proptest range")]
        fn ltf8_small_values_match_itf8(val in 0u32..0x80) {
            // Values < 128 should encode identically in ITF8 and LTF8
            let (itf8_val, itf8_len) = decode_itf8(&[val as u8]).unwrap();
            let (ltf8_val, ltf8_len) = decode_ltf8(&[val as u8]).unwrap();
            prop_assert_eq!(u64::from(itf8_val), ltf8_val);
            prop_assert_eq!(itf8_len, ltf8_len);
        }

        // For values that fit in 1–4 ITF8 bytes (< 0x10000000), ITF8 and LTF8 use
        // the same prefix bit patterns, so both decoders must decode the same bytes
        // to the same value. decode_ltf8 and read_ltf8 must also agree with each other.
        // (Values >= 0x10000000 use a 5-byte ITF8 encoding whose wire format overlaps
        // with LTF8's 5-byte range but represents different values — not tested here.)
        #[test]
        fn ltf8_matches_itf8_for_small_values_and_decoders_agree(val in 0u32..0x1000_0000u32) {
            let mut itf8_buf = [0u8; 5];
            let n = encode_itf8(val, &mut itf8_buf);
            let encoded = itf8_buf.get(..n).unwrap();

            // ITF8 and LTF8 share the 1-4 byte prefix scheme; both must give the same value.
            let (itf8_val, _) = decode_itf8(encoded).unwrap();
            let (ltf8_val, ltf8_n) = decode_ltf8(encoded).unwrap();
            prop_assert_eq!(u64::from(itf8_val), ltf8_val);
            prop_assert_eq!(ltf8_n, n);

            // decode_ltf8 and read_ltf8 must agree on the same encoded bytes.
            let read_val = read_ltf8(&mut &*encoded).unwrap();
            prop_assert_eq!(ltf8_val, read_val);
        }
    }

    // ── Comparison with noodles ─────────────────────────────────────────

    #[test]
    fn itf8_matches_noodles_known_vectors() {
        // Test vectors from noodles-cram ITF8 tests
        // These are values that noodles specifically tests
        let cases: &[(u32, &[u8])] = &[
            (0, &[0x00]),
            (1, &[0x01]),
            (127, &[0x7F]),
            (128, &[0x80, 0x80]),
            (16383, &[0xBF, 0xFF]),
            (16384, &[0xC0, 0x40, 0x00]),
        ];
        for &(expected, bytes) in cases {
            let (val, _) = decode_itf8(bytes).unwrap();
            assert_eq!(val, expected, "mismatch for bytes {bytes:?}");
        }
    }
}
