//! Tests for BAM sequence encoding/decoding, including SIMD equivalence.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    reason = "test code"
)]
#![allow(clippy::arithmetic_side_effects, reason = "test code")]
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    reason = "test code with known small values"
)]
use proptest::prelude::*;
use seqair::bam::seq::{decode_seq, decode_seq_scalar, encode_seq};

// ---- seq.decode_scalar ----

// r[verify seq.decode_scalar]
#[test]
fn decode_all_16_bases() {
    // Each nibble 0-15 maps to =ACMGRSVTWYHKDBN
    let expected = b"=ACMGRSVTWYHKDBN";
    // Pack: 0x01 = (0,1) → (=,A), 0x23 = (2,3) → (C,M), etc.
    let encoded: Vec<u8> = (0..8).map(|i| ((i * 2) << 4) | (i * 2 + 1)).collect();
    let result = decode_seq_scalar(&encoded, 16);
    assert_eq!(&result, expected);
}

// r[verify seq.decode_scalar]
#[test]
fn decode_odd_length() {
    // 3 bases from 2 bytes: 0x12 = (A,C), 0x40 = (G,=) but only 3 bases
    let encoded = [0x12, 0x40];
    let result = decode_seq_scalar(&encoded, 3);
    assert_eq!(&result, b"ACG");
}

// r[verify seq.decode_scalar]
#[test]
fn decode_empty() {
    let result = decode_seq_scalar(&[], 0);
    assert!(result.is_empty());
}

// ---- seq.decode_pair_table ----

// r[verify seq.decode_pair_table]
#[test]
fn decode_pair_table_exhaustive() {
    // Every possible byte value should decode to two valid bases
    let valid_bases = b"=ACMGRSVTWYHKDBN";
    for byte_val in 0u8..=255 {
        let result = decode_seq_scalar(&[byte_val], 2);
        assert_eq!(result.len(), 2);
        assert!(
            valid_bases.contains(&result[0]) && valid_bases.contains(&result[1]),
            "byte 0x{byte_val:02X} decoded to invalid bases: {:?}",
            result
        );
    }
}

// ---- seq.decode_simd + seq.simd_scalar_equivalence ----

// r[verify seq.simd_scalar_equivalence]
// r[verify seq.decode_simd]
// r[verify seq.decode_dispatch]
// r[verify io.platform_optimizations]
proptest! {
    #[test]
    fn simd_matches_scalar_arbitrary(
        len in 0usize..=300,
        seed in 0u64..=u64::MAX,
    ) {
        let n_bytes = len.div_ceil(2);
        // Generate deterministic encoded bytes from seed
        let encoded: Vec<u8> = (0..n_bytes)
            .map(|i| ((seed.wrapping_mul(6364136223846793005).wrapping_add(i as u64)) >> 33) as u8)
            .collect();

        let scalar = decode_seq_scalar(&encoded, len);
        let dispatched = decode_seq(&encoded, len);

        prop_assert_eq!(&dispatched, &scalar,
            "SIMD/scalar mismatch for len={}", len);
    }
}

// r[verify seq.simd_scalar_equivalence]
#[test]
fn simd_boundary_lengths() {
    for len in [0usize, 1, 2, 15, 16, 30, 31, 32, 33, 34, 63, 64, 65, 127, 128, 129, 255, 256] {
        let n_bytes = (len + 1).div_ceil(2);
        let encoded: Vec<u8> = (0..n_bytes).map(|i| i as u8).collect();

        let scalar = decode_seq_scalar(&encoded, len);
        let dispatched = decode_seq(&encoded, len);

        assert_eq!(dispatched.len(), len, "wrong length for len={len}");
        assert_eq!(dispatched, scalar, "SIMD/scalar mismatch for len={len}");
    }
}

// r[verify seq.simd_scalar_equivalence]
#[test]
fn simd_uniform_patterns() {
    let len = 128;
    let n_bytes = len / 2;

    for byte_val in [0x00, 0x11, 0x24, 0x42, 0x88, 0xFF] {
        let encoded = vec![byte_val; n_bytes];

        let scalar = decode_seq_scalar(&encoded, len);
        let dispatched = decode_seq(&encoded, len);

        assert_eq!(dispatched, scalar, "mismatch for uniform byte 0x{byte_val:02X}");
    }
}

// ---- seq.encode_scalar ----

// r[verify seq.encode_scalar]
#[test]
fn encode_roundtrip() {
    let bases = b"ACGTACGTNN";
    let encoded = encode_seq(bases);
    let decoded = decode_seq(&encoded, bases.len());
    assert_eq!(&decoded, bases);
}

// r[verify seq.encode_scalar]
#[test]
fn encode_odd_length() {
    let bases = b"ACG";
    let encoded = encode_seq(bases);
    assert_eq!(encoded.len(), 2); // 3 bases → 2 bytes
    let decoded = decode_seq(&encoded, 3);
    assert_eq!(&decoded, b"ACG");
}

// r[verify seq.encode_scalar]
// Verifies the BAM spec (SAM1 §4.2.4) nibble encoding table:
// =ACMGRSVTWYHKDBN → nibble values 0–15, so A=1, C=2, G=4, T=8, N=15.
proptest! {
    #[test]
    fn encode_produces_spec_nibble_values(
        bases in prop::collection::vec(
            prop::sample::select(vec![b'A', b'C', b'G', b'T', b'N']),
            0..=200,
        ),
    ) {
        fn spec_nibble(b: u8) -> u8 {
            match b {
                b'A' => 1,
                b'C' => 2,
                b'G' => 4,
                b'T' => 8,
                b'N' => 15,
                _ => panic!("unexpected base {b}"),
            }
        }

        let encoded = encode_seq(&bases);

        // Verify each encoded byte has the correct nibble values per the BAM spec.
        for i in 0..bases.len() / 2 {
            let byte = encoded.get(i).copied().unwrap();
            let hi = byte >> 4;
            let lo = byte & 0x0F;
            prop_assert_eq!(
                hi,
                spec_nibble(*bases.get(i * 2).unwrap()),
                "high nibble mismatch at pair {}",
                i
            );
            prop_assert_eq!(
                lo,
                spec_nibble(*bases.get(i * 2 + 1).unwrap()),
                "low nibble mismatch at pair {}",
                i
            );
        }
        // For odd-length sequences the last byte's high nibble encodes the final base.
        if bases.len() % 2 == 1 {
            let last_byte = encoded.get(bases.len() / 2).copied().unwrap();
            prop_assert_eq!(
                last_byte >> 4,
                spec_nibble(*bases.last().unwrap()),
                "high nibble mismatch for final odd base"
            );
            prop_assert_eq!(last_byte & 0x0F, 0u8, "low nibble of padding byte must be 0");
        }

        // Secondary: roundtrip check.
        let decoded = decode_seq(&encoded, bases.len());
        prop_assert_eq!(&decoded, &bases);
    }
}
