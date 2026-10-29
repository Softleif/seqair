//! `CramStreamReader` and `CramContainerDecoder` on arbitrary bytes: reading
//! and decoding containers never panics, an error poisons the stream
//! (`r[cram.stream.framing+2]`), and a failed decode leaves its output as it
//! was (`r[cram.stream.decoder]`). Containers decode against the reference of
//! the `tests/data/cram-stream` fixtures, which make good seeds.
#![no_main]

use libfuzzer_sys::fuzz_target;
use seqair::cram::stream::{CramStreamError, CramStreamReader};
use std::path::Path;

const FASTA: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../tests/data/cram-stream/fidelity.fa");

fuzz_target!(|data: &[u8]| {
    let Ok(mut reader) = CramStreamReader::new(data) else { return };
    let mut decoder = reader.decoder(Path::new(FASTA)).expect("the fixture FASTA opens");
    let mut out = Vec::new();
    for _ in 0..10_000 {
        match reader.next_container() {
            Ok(Some(container)) => {
                if decoder.decode(&container, &mut out).is_err() {
                    assert!(out.is_empty(), "a failed decode changed its output");
                }
                for record in out.drain(..) {
                    std::hint::black_box((record.qname().len(), record.seq_len()));
                }
            }
            Ok(None) => {
                assert!(matches!(reader.next_container(), Ok(None)), "the end stays the end");
                break;
            }
            Err(_) => {
                assert!(
                    matches!(reader.next_container(), Err(CramStreamError::Poisoned)),
                    "not poisoned"
                );
                break;
            }
        }
    }
});
