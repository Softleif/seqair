#![no_main]

use libfuzzer_sys::fuzz_target;
use seqair::bam::aux::{self, AuxValue};
use seqair::bam::aux_data::AuxData;

/// Round-trip an `AuxValue` through `AuxData` and verify it comes back intact.
///
/// Integer types may select a wider BAM encoding (e.g., i8 42 → U8 42).
/// We accept `as_i64()` equivalence for integers. Floats compare bit patterns
/// since NaN ≠ NaN under `PartialEq`.
fn assert_value_roundtrip(original: &AuxValue<'_>, roundtripped: Option<AuxValue<'_>>) {
    let rt = roundtripped.expect("tag should survive round-trip");
    match (original, &rt) {
        (AuxValue::Float(a), AuxValue::Float(b)) => {
            assert_eq!(a.to_bits(), b.to_bits(), "float bit-pattern mismatch");
        }
        (AuxValue::Double(a), AuxValue::Double(b)) => {
            assert_eq!(a.to_bits(), b.to_bits(), "double bit-pattern mismatch");
        }
        (AuxValue::I8(a), _) => assert_eq!(rt.as_i64(), Some(i64::from(*a))),
        (AuxValue::U8(a), _) => assert_eq!(rt.as_i64(), Some(i64::from(*a))),
        (AuxValue::I16(a), _) => assert_eq!(rt.as_i64(), Some(i64::from(*a))),
        (AuxValue::U16(a), _) => assert_eq!(rt.as_i64(), Some(i64::from(*a))),
        (AuxValue::I32(a), _) => assert_eq!(rt.as_i64(), Some(i64::from(*a))),
        (AuxValue::U32(a), _) => assert_eq!(rt.as_i64(), Some(i64::from(*a))),
        _ => assert_eq!(original, &rt, "non-integer round-trip mismatch"),
    }
}

/// Decode a B-array byte slice (parser output) into a typed Vec for the
/// matching typed setter. The parser only emits length-aligned slices.
fn decode_i16(bytes: &[u8]) -> Vec<i16> {
    bytes.chunks_exact(2).map(|c| i16::from_le_bytes([c[0], c[1]])).collect()
}
fn decode_u16(bytes: &[u8]) -> Vec<u16> {
    bytes.chunks_exact(2).map(|c| u16::from_le_bytes([c[0], c[1]])).collect()
}
fn decode_i32(bytes: &[u8]) -> Vec<i32> {
    bytes.chunks_exact(4).map(|c| i32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect()
}
fn decode_u32(bytes: &[u8]) -> Vec<u32> {
    bytes.chunks_exact(4).map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect()
}
fn decode_f32(bytes: &[u8]) -> Vec<f32> {
    bytes.chunks_exact(4).map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect()
}

/// Set a single tag value into `AuxData`. Returns `false` if the writer
/// rejected it (e.g. `set_char` on a non-printable byte) so the caller can
/// drop the tag from the expected-roundtrip set.
fn set_one(aux: &mut AuxData, tag: [u8; 2], value: &AuxValue<'_>) -> bool {
    match value {
        AuxValue::Char(v) => aux.set_char(tag, *v).is_ok(),
        AuxValue::String(s) => {
            aux.set_string(tag, s);
            true
        }
        AuxValue::Hex(h) => {
            aux.set_hex(tag, h);
            true
        }
        AuxValue::I8(v) => aux.set_int(tag, i64::from(*v)).is_ok(),
        AuxValue::U8(v) => aux.set_int(tag, i64::from(*v)).is_ok(),
        AuxValue::I16(v) => aux.set_int(tag, i64::from(*v)).is_ok(),
        AuxValue::U16(v) => aux.set_int(tag, i64::from(*v)).is_ok(),
        AuxValue::I32(v) => aux.set_int(tag, i64::from(*v)).is_ok(),
        AuxValue::U32(v) => aux.set_int(tag, i64::from(*v)).is_ok(),
        AuxValue::Float(v) => {
            aux.set_float(tag, *v);
            true
        }
        AuxValue::Double(v) => {
            aux.set_double(tag, *v);
            true
        }
        AuxValue::ArrayI8(a) => {
            let typed: Vec<i8> = a.iter().map(|&b| b.cast_signed()).collect();
            aux.set_array_i8(tag, &typed).is_ok()
        }
        AuxValue::ArrayU8(a) => aux.set_array_u8(tag, a).is_ok(),
        AuxValue::ArrayI16(a) => aux.set_array_i16(tag, &decode_i16(a)).is_ok(),
        AuxValue::ArrayU16(a) => aux.set_array_u16(tag, &decode_u16(a)).is_ok(),
        AuxValue::ArrayI32(a) => aux.set_array_i32(tag, &decode_i32(a)).is_ok(),
        AuxValue::ArrayU32(a) => aux.set_array_u32(tag, &decode_u32(a)).is_ok(),
        AuxValue::ArrayFloat(a) => aux.set_array_f32(tag, &decode_f32(a)).is_ok(),
    }
}

fuzz_target!(|data: &[u8]| {
    // ── 1. Parse raw bytes (no panics) ──
    let tags: Vec<_> = aux::iter_tags(data).collect();

    // Also exercise find_tag for a few fixed tags
    let _ = aux::find_tag(data, *b"RG");
    let _ = aux::find_tag(data, *b"NM");
    let _ = aux::find_tag(data, *b"ZZ");

    // ── 2. Deduplicate: AuxData overwrites on set, so only the LAST
    //    occurrence of each tag name survives the round-trip.
    let mut last_value: rustc_hash::FxHashMap<[u8; 2], &AuxValue<'_>> =
        rustc_hash::FxHashMap::default();
    for (tag, value) in &tags {
        last_value.insert(*tag, value);
    }

    // ── 3. Build AuxData from deduplicated tags (last-wins) ──
    let mut aux = AuxData::new();
    let mut stored: Vec<[u8; 2]> = Vec::new();
    for (tag, value) in &last_value {
        if set_one(&mut aux, *tag, value) {
            stored.push(*tag);
        }
    }

    // ── 4. Verify round-trip for every stored tag ──
    for tag in &stored {
        let expected = last_value.get(tag).unwrap();
        assert_value_roundtrip(expected, aux.get(*tag));
    }

    // ── 5. Sequential removal → empty ──
    for tag in &stored {
        aux.remove(*tag);
    }
    assert!(aux.is_empty(), "all tags removed but aux not empty");
});
