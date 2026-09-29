#![no_main]

use arbitrary::Arbitrary;
use libfuzzer_sys::fuzz_target;
use seqair::cram::{
    bitstream::BitReader,
    encoding::{ByteEncoding, DecodeContext, ExternalCursor, HuffmanTable},
};
use seqair_types::SmallVec;

#[derive(Arbitrary, Debug)]
struct EncodingInput {
    alphabet: Vec<i32>,
    bit_lengths: Vec<u32>,
    core_data: Vec<u8>,
    external_data: Vec<u8>,
    external_ops: Vec<ExternalOp>,
    beta_offset: i32,
    beta_bits: u8,
    beta_count: u8,
}

#[derive(Arbitrary, Debug)]
enum ExternalOp {
    ReadByte,
    ReadItf8,
    ReadBytesUntil(u8),
    ReadBytes(u8),
    Remaining,
}

fuzz_target!(|input: EncodingInput| {
    // Fuzz HuffmanTable construction and decoding
    if let Ok(table) = HuffmanTable::new(&input.alphabet, &input.bit_lengths) {
        let mut reader = BitReader::new(&input.core_data);
        let _ = table.decode(&mut reader);
        // Try decoding multiple symbols
        let mut reader2 = BitReader::new(&input.core_data);
        for _ in 0..8 {
            if table.decode(&mut reader2).is_none() {
                break;
            }
        }
    }

    // Fuzz BETA byte decoding from the core bit stream; widths above 32 are
    // rejected at parse time, so construct only what parsing would accept.
    let beta =
        ByteEncoding::Beta { offset: input.beta_offset, bits: u32::from(input.beta_bits % 33) };
    let mut ctx = DecodeContext::new(&input.core_data, SmallVec::new());
    let _ = beta.decode(&mut ctx);
    let mut buf = Vec::new();
    let _ = beta.decode_n_into(&mut ctx, usize::from(input.beta_count), &mut buf);

    // Fuzz ExternalCursor operations
    let mut cursor = ExternalCursor::new(input.external_data);
    let mut buf = Vec::new();
    for op in &input.external_ops {
        match op {
            ExternalOp::ReadByte => {
                let _ = cursor.read_byte();
            }
            ExternalOp::ReadItf8 => {
                let _ = cursor.read_itf8();
            }
            ExternalOp::ReadBytesUntil(stop) => {
                buf.clear();
                let _ = cursor.read_bytes_until_into(*stop, &mut buf);
            }
            ExternalOp::ReadBytes(n) => {
                // Clamp to avoid huge allocations
                let n = (*n as usize) % 64;
                buf.clear();
                let _ = cursor.read_bytes_into(n, &mut buf);
            }
            ExternalOp::Remaining => {
                let _ = cursor.remaining();
            }
        }
    }
});
