#![no_main]

use arbitrary::Arbitrary;
use libfuzzer_sys::fuzz_target;
use seqair::cram::compression_header::CompressionHeader;
use seqair::cram::encoding::{ByteArrayEncoding, ByteEncoding, IntEncoding};
use seqair::cram::slice::SliceHeader;

#[derive(Arbitrary, Debug)]
enum HeaderInput {
    CompressionHeader(Vec<u8>),
    SliceHeader(Vec<u8>),
    IntEncoding(Vec<u8>),
    ByteEncoding(Vec<u8>),
    ByteArrayEncoding(Vec<u8>),
}

fuzz_target!(|input: HeaderInput| {
    match input {
        HeaderInput::CompressionHeader(data) => {
            let _ = CompressionHeader::parse(&data);
        }
        HeaderInput::SliceHeader(data) => {
            let _ = SliceHeader::parse(&data);
        }
        HeaderInput::IntEncoding(data) => {
            let mut cursor: &[u8] = &data;
            let _ = IntEncoding::parse(&mut cursor);
        }
        HeaderInput::ByteEncoding(data) => {
            let mut cursor: &[u8] = &data;
            let _ = ByteEncoding::parse(&mut cursor);
        }
        HeaderInput::ByteArrayEncoding(data) => {
            let mut cursor: &[u8] = &data;
            let _ = ByteArrayEncoding::parse(&mut cursor);
        }
    }
});
