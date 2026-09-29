#![no_main]

use arbitrary::Arbitrary;
use libfuzzer_sys::fuzz_target;
use seqair::cram::bitstream::BitReader;

#[derive(Arbitrary, Debug)]
struct BitstreamInput {
    data: Vec<u8>,
    ops: Vec<BitstreamOp>,
}

#[derive(Arbitrary, Debug)]
enum BitstreamOp {
    ReadBit,
    ReadBits(u8),
    ReadBitsI32(u8),
    RemainingBits,
}

fuzz_target!(|input: BitstreamInput| {
    let mut reader = BitReader::new(&input.data);
    for op in &input.ops {
        match op {
            BitstreamOp::ReadBit => {
                let _ = reader.read_bit();
            }
            BitstreamOp::ReadBits(n) => {
                // Clamp to 0..=32 as documented
                let bits = (*n as u32) % 33;
                let _ = reader.read_bits(bits);
            }
            BitstreamOp::ReadBitsI32(n) => {
                let bits = (*n as u32) % 33;
                let _ = reader.read_bits_i32(bits);
            }
            BitstreamOp::RemainingBits => {
                let _ = reader.remaining_bits();
            }
        }
    }
});
