//! Shared per-byte primitives for CRAM codec modules (`rans`, `rans_nx16`,
//! `tok3`). Each helper returns the smallest possible result type so the
//! per-byte hot path stays in registers.
//!
//! ## Why `Option<u8>` instead of `Result<u8, CramError>`
//!
//! The main [`CramError`](super::reader::CramError) is 80 bytes (sized to
//! fit variants like `Open { path: PathBuf, source: io::Error }`), so
//! `Result<u8, CramError>` is also 80 bytes — passed via memory/sret in
//! the non-inlined ABI. `Option<u8>` is 2 bytes, in a register. samply
//! flagged `core::ptr::drop_in_place<CramError>` next to `read_u8` /
//! `renormalize` because the eager `ok_or(CramError::...)` form
//! constructed and then dropped a `CramError` on every successful read.
//!
//! Callers materialize a `CramError` only on the failure path, with
//! `ok_or_else(|| CramError::Truncated { context: ... })` at the
//! outermost function that knows the context.
//!
//! `read_uint7` is different: it has two distinct failure modes
//! (truncation vs overflow) and benefits from a narrow typed error
//! ([`Uint7Error`]) — the enum is 1 byte so `Result<u32, Uint7Error>`
//! is the same size as `Option<u32>` and incurs no drop cost.

/// Failure modes for [`read_uint7`].
#[derive(Debug, thiserror::Error, PartialEq, Eq)]
pub enum Uint7Error {
    #[error("truncated uint7")]
    Truncated,
    #[error("uint7 overflow: more than 5 continuation bytes")]
    Overflow,
}

#[inline]
pub fn read_u8(src: &mut &[u8]) -> Option<u8> {
    let (&b, rest) = src.split_first()?;
    *src = rest;
    Some(b)
}

#[inline]
pub fn read_u16_le(src: &mut &[u8]) -> Option<u16> {
    let (head, rest) = src.split_first_chunk::<2>()?;
    *src = rest;
    Some(u16::from_le_bytes(*head))
}

#[inline]
pub fn read_u32_le(src: &mut &[u8]) -> Option<u32> {
    let (head, rest) = src.split_first_chunk::<4>()?;
    *src = rest;
    Some(u32::from_le_bytes(*head))
}

#[inline]
pub fn split_off<'a>(src: &mut &'a [u8], len: usize) -> Option<&'a [u8]> {
    let (head, rest) = src.split_at_checked(len)?;
    *src = rest;
    Some(head)
}

/// Decode a uint7 variable-length integer (MSB-first, max 5 bytes for u32).
///
/// `r[impl cram.codec.uint7_bounded]`
#[allow(clippy::arithmetic_side_effects, reason = "loop bounded to 5 iterations of <<7")]
pub fn read_uint7(src: &mut &[u8]) -> Result<u32, Uint7Error> {
    let mut n: u32 = 0;
    for _ in 0..5 {
        let b = u32::from(read_u8(src).ok_or(Uint7Error::Truncated)?);
        // ITF-8 / uint7 is MSB-first: shift accumulator left so each new
        // byte's 7 data bits become the new low bits.
        n = (n << 7) | (b & 0x7f);
        if b & 0x80 == 0 {
            return Ok(n);
        }
    }
    Err(Uint7Error::Overflow)
}

#[cfg(test)]
#[allow(
    clippy::arithmetic_side_effects,
    reason = "test code: bounded by generator ranges and fixed inputs"
)]
mod tests {
    use super::*;
    use hegel::prelude::*;

    // r[verify cram.codec.uint7_bounded]
    #[test]
    fn read_uint7_overflow_returns_error() {
        // 6 continuation bytes (all with high bit set) exceed the 5-iteration limit.
        let mut src: &[u8] = &[0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x00];
        assert_eq!(read_uint7(&mut src), Err(Uint7Error::Overflow));
    }

    #[test]
    fn read_uint7_truncated_returns_error() {
        let mut src: &[u8] = &[0x80, 0x80];
        assert_eq!(read_uint7(&mut src), Err(Uint7Error::Truncated));
    }

    #[test]
    fn read_uint7_five_bytes_ok() {
        let mut src: &[u8] = &[0x80, 0x80, 0x80, 0x80, 0x01];
        assert!(read_uint7(&mut src).is_ok());
    }

    #[test]
    fn read_uint7_spec_vectors() {
        // Hard-coded byte sequences derived from htscodecs `var_put_u32`
        // (htslib/htscodecs/htscodecs/varint.h:206, BIG_END / MSB-first).
        // Independent oracle — see commit 701c9ed2 for context.
        let cases: &[(u32, &[u8])] = &[
            (0, &[0x00]),
            (1, &[0x01]),
            (127, &[0x7F]),
            (128, &[0x81, 0x00]),
            (200, &[0x81, 0x48]),
            (16_383, &[0xFF, 0x7F]),
            (16_384, &[0x81, 0x80, 0x00]),
            (0x12345, &[0x84, 0xC6, 0x45]),
            ((1u32 << 28) - 1, &[0xFF, 0xFF, 0xFF, 0x7F]),
            (1u32 << 28, &[0x81, 0x80, 0x80, 0x80, 0x00]),
            (u32::MAX, &[0x8F, 0xFF, 0xFF, 0xFF, 0x7F]),
        ];
        for (val, encoded) in cases {
            let mut cur: &[u8] = encoded;
            let decoded = read_uint7(&mut cur)
                .unwrap_or_else(|e| panic!("decode of {encoded:02x?} failed: {e:?}"));
            assert_eq!(decoded, *val, "decoded value mismatch for {encoded:02x?}");
            assert!(cur.is_empty(), "wrong byte count for {val} ({encoded:02x?})");
        }
    }

    #[test]
    fn read_u8_u16_u32_basic() {
        let bytes = [0x12u8, 0x34, 0x56, 0x78, 0x9A, 0xBC, 0xDE];
        let mut cur: &[u8] = &bytes;
        assert_eq!(read_u8(&mut cur), Some(0x12));
        assert_eq!(read_u16_le(&mut cur), Some(0x5634));
        assert_eq!(read_u32_le(&mut cur), Some(0xDEBC9A78));
        assert_eq!(read_u8(&mut cur), None);
    }

    #[test]
    fn split_off_basic() {
        let bytes = [1u8, 2, 3, 4, 5];
        let mut cur: &[u8] = &bytes;
        let head = split_off(&mut cur, 3).unwrap();
        assert_eq!(head, &[1, 2, 3]);
        assert_eq!(cur, &[4, 5]);
        // Asking for more than available returns None and leaves cursor alone.
        let saved = cur;
        assert_eq!(split_off(&mut cur, 10), None);
        assert_eq!(cur, saved);
    }

    #[hegel::test]
    fn read_uint7_decode_matches_msb_first_assembly(tc: TestCase) {
        let val = tc.draw(gs::integers::<u32>());
        // MSB-first assembly: take the value's 5x 7-bit groups from
        // most significant down, set continuation bit on all but the
        // last non-leading byte, then verify decoder produces val.
        // Skip leading zero groups (canonical encoding).
        let mut groups: [u8; 5] = [
            ((val >> 28) & 0x7F) as u8,
            ((val >> 21) & 0x7F) as u8,
            ((val >> 14) & 0x7F) as u8,
            ((val >> 7) & 0x7F) as u8,
            (val & 0x7F) as u8,
        ];
        // Drop leading-zero groups, but keep at least one group.
        let start = groups.iter().position(|&g| g != 0).unwrap_or(4);
        let n = 5 - start;
        // Set continuation bit on all but the last byte.
        for g in groups.iter_mut().skip(start).take(n - 1) {
            *g |= 0x80;
        }
        let encoded = &groups[start..];
        let mut cur: &[u8] = encoded;
        assert_eq!(read_uint7(&mut cur).unwrap(), val);
        assert!(cur.is_empty());
    }
}
