#![no_main]

use libfuzzer_sys::fuzz_target;
use seqair::bam::record::BamRecord;
use seqair::bam::record_store::RecordStore;

fuzz_target!(|data: &[u8]| {
    // Fuzz the BAM record decoder
    let _ = BamRecord::decode(data);

    // Also fuzz RecordStore::push_raw which decodes and stores
    let mut store = RecordStore::default();
    let _ = store.push_raw(data);
});
