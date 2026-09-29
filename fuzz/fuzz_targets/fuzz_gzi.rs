#![no_main]

use libfuzzer_sys::fuzz_target;
use seqair::fasta::GziIndex;

fuzz_target!(|data: &[u8]| {
    let _ = GziIndex::from_bytes(data);
});
