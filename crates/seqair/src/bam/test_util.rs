//! Shared test helpers for the BAM module.
//!
//! Test-only utilities that exercise the production decode path
//! (`RecordStore::push_raw`) so write-side tests can round-trip through the
//! same code real region loading runs.

#![cfg(test)]

use super::bgzf::BgzfReader;
use super::record_idx::RecordIdx;
use super::record_store::RecordStore;
use std::io::{Read, Seek};

/// The record index `n`. Tests address records by literal position; this is
/// the shorthand for the `Option` that `RecordIdx::new` returns.
#[allow(dead_code, reason = "used by sibling test modules; cfg gates may hide some call sites")]
#[track_caller]
pub(crate) fn ri(n: u32) -> RecordIdx {
    RecordIdx::new(n).expect("u32::MAX is not a record index")
}

/// Decode a single BAM record (sans the 4-byte `block_size` prefix) into a
/// fresh `RecordStore` via the production `push_raw` path.
///
/// Panics if the buffer is rejected, since these helpers exist for tests
/// that build the input themselves.
#[allow(
    dead_code,
    reason = "used by sibling test modules; rustc cfg gates may hide some call sites"
)]
pub(crate) fn decode_into_store(buf: &[u8]) -> RecordStore {
    let mut store = RecordStore::new();
    let idx = store
        .push_raw(buf, &mut ())
        .expect("push_raw returned an error")
        .expect("record was filtered out");
    assert_eq!(idx, RecordIdx::ZERO, "fresh store must produce index 0");
    store
}

/// Read one BAM record from a BGZF stream (`block_size` + bytes) and push it
/// into the given store via `push_raw`. Returns the assigned index.
#[allow(
    dead_code,
    reason = "used by sibling test modules; rustc cfg gates may hide some call sites"
)]
pub(crate) fn push_one_record_from_bgzf<R: Read + Seek>(
    reader: &mut BgzfReader<R>,
    store: &mut RecordStore,
) -> RecordIdx {
    let block_size: i32 = reader.read_i32().expect("read block_size");
    assert!(block_size > 0, "block_size must be positive, got {block_size}");
    let len = usize::try_from(block_size).expect("block_size fits in usize");
    let mut rec_data = vec![0u8; len];
    reader.read_exact_into(&mut rec_data).expect("read record bytes");
    store
        .push_raw(&rec_data, &mut ())
        .expect("push_raw returned an error")
        .expect("record was filtered out")
}
