#![no_main]

use libfuzzer_sys::fuzz_target;
use seqair::bam::RecordStore;

// Fuzz `RecordStore::push_raw` — the production decode path that every reader
// in the crate runs. There is intentionally no separate `BamRecord::decode`
// to fuzz; that struct was removed in favor of slab-only decode.
fuzz_target!(|data: &[u8]| {
    let mut store = RecordStore::default();
    let _ = store.push_raw(data, &mut ());
});
