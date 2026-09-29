#![no_main]

use libfuzzer_sys::fuzz_target;
use seqair::fasta::FastaIndex;

fuzz_target!(|data: &[u8]| {
    if let Ok(text) = std::str::from_utf8(data) {
        let _ = FastaIndex::from_contents(text);
    }
});
