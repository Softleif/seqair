//! `AuxData` and `RawBamRecord::edit_aux` on aux blocks that are not a clean
//! sequence of tags: well-formed tags cut at any byte, arrays whose count
//! runs past the data, strings without their NUL, and plain noise.
//!
//! A value running past the end of the block used to make `remove` drain a
//! range past the end and panic, and `edit_aux` reached it on any record.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "test code"
)]

#[path = "support/raw_record_gen.rs"]
mod raw_record_gen;

use hegel::TestCase;
use hegel::generators as gs;
use raw_record_gen::{GenTag, HTSLIB, arb_record, arb_tag, header};
use seqair::bam::aux::validate_aux;
use seqair::bam::{AuxData, RawBamRecord, RawRecordEditError};

/// An aux block and what the generator knows about it.
#[derive(Debug, Clone)]
struct Block {
    bytes: Vec<u8>,
    /// The tags `bytes` consists of, when it is exactly well-formed tags.
    tags: Option<Vec<GenTag>>,
    /// Whether the block is well-formed, when the generator can tell.
    valid: Option<bool>,
}

fn encode(tags: &[GenTag]) -> Vec<u8> {
    let mut out = Vec::new();
    for tag in tags {
        tag.encode(&mut out);
    }
    out
}

fn n<T: hegel::generators::Integer>(lo: T, hi: T) -> gs::IntegerGenerator<T> {
    gs::integers::<T>().min_value(lo).max_value(hi)
}

#[hegel::composite]
fn arb_block(tc: &TestCase) -> Block {
    let tags: Vec<GenTag> =
        (0..tc.draw_silent(n(0usize, 5))).map(|_| arb_tag(tc, HTSLIB)).collect();
    let whole = encode(&tags);
    match tc.draw_silent(n(0u8, 5)) {
        0 => Block { bytes: whole, tags: Some(tags), valid: Some(true) },
        // Cut anywhere: valid only on a tag boundary, and then it is the
        // tags before the cut.
        1 => {
            let cut = tc.draw_silent(n(0, whole.len()));
            let mut end = 0;
            let mut kept = Vec::new();
            for tag in &tags {
                if end + 2 + tag.value.len() > cut {
                    break;
                }
                end += 2 + tag.value.len();
                kept.push(tag.clone());
            }
            let on_boundary = end == cut;
            Block {
                bytes: whole[..cut].to_vec(),
                tags: on_boundary.then_some(kept),
                valid: Some(on_boundary),
            }
        }
        // A `B` array whose count needs more bytes than follow it.
        2 => {
            let (subtype, size) = tc.draw_silent(gs::sampled_from(vec![
                (b'c', 1u64),
                (b'C', 1),
                (b's', 2),
                (b'S', 2),
                (b'i', 4),
                (b'I', 4),
                (b'f', 4),
            ]));
            let elements = tc.draw_silent(gs::binary().max_size(12));
            let min_count = elements.len() as u64 / size + 1;
            let count = tc.draw_silent(n(min_count, u64::from(u32::MAX)));
            let mut bytes = whole;
            bytes.extend_from_slice(b"MLB");
            bytes.push(subtype);
            bytes.extend_from_slice(&u32::try_from(count).unwrap().to_le_bytes());
            bytes.extend_from_slice(&elements);
            Block { bytes, tags: None, valid: Some(false) }
        }
        // A `Z` or `H` string that runs to the end without a NUL.
        3 => {
            let mut bytes = whole;
            bytes.extend_from_slice(b"XZ");
            bytes.push(tc.draw_silent(gs::sampled_from(b"ZH".to_vec())));
            let len = tc.draw_silent(n(0usize, 10));
            bytes.extend((0..len).map(|_| tc.draw_silent(n(1u8, 255))));
            Block { bytes, tags: None, valid: Some(false) }
        }
        // Well-formed tags, then noise.
        4 => {
            let mut bytes = whole;
            bytes.extend(tc.draw_silent(gs::binary().min_size(1).max_size(16)));
            Block { bytes, tags: None, valid: None }
        }
        _ => Block { bytes: tc.draw_silent(gs::binary().max_size(48)), tags: None, valid: None },
    }
}

/// Keys worth asking for: every key in the block's first bytes and each
/// tag's, the ones the damage modes write, and one drawn.
fn probe_keys(tc: &TestCase, block: &Block) -> Vec<[u8; 2]> {
    let mut keys =
        vec![*b"ML", *b"XZ", *b"NM", tc.draw(gs::sampled_from(vec![*b"MM", *b"q9", *b"XM"]))];
    if let [a, b, ..] = block.bytes[..] {
        keys.push([a, b]);
    }
    if let Some(tags) = &block.tags {
        keys.extend(tags.iter().map(|t| t.key));
    }
    keys
}

/// `tags` with the first one named `key` taken out, as `AuxData::remove` does.
fn model_remove(tags: &[GenTag], key: [u8; 2]) -> Vec<GenTag> {
    let mut out = tags.to_vec();
    if let Some(i) = out.iter().position(|t| t.key == key) {
        out.remove(i);
    }
    out
}

// r[verify bam.owned_record.aux_data]
// r[verify bam.record.aux_checked]
/// On any bytes, `get`, `remove` and the setters return instead of
/// panicking; `remove` only ever cuts one span out; on well-formed tags it
/// takes out exactly the first tag of that name, and nothing for an absent
/// one. `validate_aux` agrees with the generator wherever that knows.
#[hegel::test]
fn aux_data_handles_any_bytes(tc: TestCase) {
    let block = tc.draw(arb_block());
    if let Some(valid) = block.valid {
        assert_eq!(validate_aux(&block.bytes).is_ok(), valid, "validate_aux");
    }
    for key in probe_keys(&tc, &block) {
        let aux = AuxData::from_bytes(block.bytes.clone());
        let found = aux.get(key).is_some();
        if let Some(tags) = &block.tags {
            assert_eq!(found, tags.iter().any(|t| t.key == key), "get {key:?}");
        }

        let mut removed = aux.clone();
        removed.remove(key);
        let after = removed.as_bytes();
        let kept = after.len();
        let prefix = block.bytes.iter().zip(after).take_while(|(a, b)| a == b).count();
        assert!(kept <= block.bytes.len(), "remove grew the block");
        assert_eq!(
            after[prefix..],
            block.bytes[block.bytes.len() - (kept - prefix)..],
            "remove cut more than one span"
        );
        if let Some(tags) = &block.tags {
            assert_eq!(after, encode(&model_remove(tags, key)), "remove {key:?}");
        }

        let mut set = aux.clone();
        set.set_string(key, b"abc");
        set.set_int(key, -7).unwrap();
        set.set_array_u8(key, &[1, 2, 3]).unwrap();
        set.set_float(key, 1.5);
        std::hint::black_box(set.get(key));
    }
}

// r[verify bam.raw_record.edit_aux]
/// `edit_aux` refuses a record exactly when `validate_aux` refuses its aux
/// block, as `MalformedAux` and with the payload untouched; otherwise the
/// edit lands as the tag model says.
#[hegel::test]
fn edit_aux_refuses_exactly_the_malformed_blocks(tc: TestCase) {
    let mut g = tc.draw(arb_record(HTSLIB));
    // Keep a long CIGAR's `CG` tag whole and in front, so the record stays
    // valid whatever is done to the rest.
    let cg: Vec<GenTag> = g.tags.drain(..).filter(|t| t.key == *b"CG").collect();
    let fixed = g.payload();
    let block = tc.draw(arb_block());
    let mut whole = encode(&cg);
    whole.extend_from_slice(&block.bytes);
    let mut payload = fixed.clone();
    payload.extend_from_slice(&whole);
    let mut record = RawBamRecord::from_payload(payload, &header()).expect("aux is not validated");
    let before = record.as_payload().to_vec();

    // Not `CG`: taking it out is the separate `LongCigarChanged` refusal.
    let keys: Vec<[u8; 2]> = probe_keys(&tc, &block).into_iter().filter(|k| k != b"CG").collect();
    let key = tc.draw(gs::sampled_from(keys));
    let result = record.edit_aux(|aux| {
        aux.remove(key);
        aux.set_string(*b"XR", b"CT");
        Ok(())
    });

    match validate_aux(&whole) {
        Err(_) => {
            assert!(matches!(result, Err(RawRecordEditError::MalformedAux(_))), "{result:?}");
            assert_eq!(record.as_payload(), before, "a refused edit changed the record");
        }
        Ok(()) => {
            result.unwrap();
            if let Some(tags) = &block.tags {
                let mut model = cg.clone();
                model.extend(tags.iter().cloned());
                let mut model = model_remove(&model_remove(&model, key), *b"XR");
                model.push(GenTag { key: *b"XR", value: b"ZCT\0".to_vec() });
                let mut want = fixed;
                want.extend(encode(&model));
                assert_eq!(record.as_payload(), want);
            }
        }
    }
}
