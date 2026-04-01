#![no_main]

use arbitrary::Arbitrary;
use libfuzzer_sys::fuzz_target;
use seqair::cram::{rans, rans_nx16, tok3};

#[derive(Arbitrary, Debug)]
enum CodecInput {
    Rans4x8(Vec<u8>),
    RansNx16 { data: Vec<u8>, uncompressed_size: u16 },
    Tok3(Vec<u8>),
}

fuzz_target!(|input: CodecInput| {
    match input {
        CodecInput::Rans4x8(data) => {
            let _ = rans::decode(&data);
        }
        CodecInput::RansNx16 { data, uncompressed_size } => {
            let _ = rans_nx16::decode(&data, uncompressed_size as usize);
        }
        CodecInput::Tok3(data) => {
            let _ = tok3::decode(&data);
        }
    }
});
