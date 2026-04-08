//! Read individual bits from a byte slice. [`BitReader`] reads MSB-first and is used
//! to decode Huffman-coded data series from CRAM core data blocks.

// r[impl cram.bitstream]
/// MSB-first bit reader over a byte slice.
///
/// Used for reading the core data block in a CRAM slice.
/// Bits are read from the most significant bit first within each byte.
pub struct BitReader<'a> {
    data: &'a [u8],
    byte_pos: usize,
    bit_pos: u8, // 0 = MSB (0x80), 7 = LSB (0x01)
}

impl<'a> BitReader<'a> {
    pub fn new(data: &'a [u8]) -> Self {
        Self { data, byte_pos: 0, bit_pos: 0 }
    }

    /// Read a single bit, returning 0 or 1.
    pub fn read_bit(&mut self) -> Option<u8> {
        let byte = self.data.get(self.byte_pos)?;
        let bit = (byte >> (7u8.checked_sub(self.bit_pos)?)) & 1;
        self.bit_pos = self.bit_pos.checked_add(1)?;
        if self.bit_pos == 8 {
            self.bit_pos = 0;
            self.byte_pos = self.byte_pos.checked_add(1)?;
        }
        Some(bit)
    }

    /// Read `n` bits as a u32 (MSB-first). Max 32 bits.
    ///
    /// A zero-bit read returns 0 (used for single-symbol Huffman codes).
    pub fn read_bits(&mut self, n: u32) -> Option<u32> {
        if n == 0 {
            return Some(0);
        }
        let mut val = 0u32;
        for _ in 0..n {
            val = (val << 1) | u32::from(self.read_bit()?);
        }
        Some(val)
    }

    /// Read `n` bits as an i32 (MSB-first, unsigned interpretation).
    pub fn read_bits_i32(&mut self, n: u32) -> Option<i32> {
        self.read_bits(n).map(|v| v.cast_signed())
    }

    /// Returns the number of bits remaining.
    pub fn remaining_bits(&self) -> usize {
        let total_bits = self.data.len().saturating_mul(8);
        let consumed = self.byte_pos.saturating_mul(8).saturating_add(self.bit_pos as usize);
        total_bits.saturating_sub(consumed)
    }
}

#[cfg(test)]
#[allow(
    clippy::arithmetic_side_effects,
    reason = "test-only arithmetic on bounded constants and small test values"
)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    // r[verify cram.bitstream]
    #[test]
    fn read_single_bits() {
        let data = [0b1100_1010];
        let mut reader = BitReader::new(&data);
        assert_eq!(reader.read_bit(), Some(1)); // MSB
        assert_eq!(reader.read_bit(), Some(1));
        assert_eq!(reader.read_bit(), Some(0));
        assert_eq!(reader.read_bit(), Some(0));
        assert_eq!(reader.read_bit(), Some(1));
        assert_eq!(reader.read_bit(), Some(0));
        assert_eq!(reader.read_bit(), Some(1));
        assert_eq!(reader.read_bit(), Some(0)); // LSB
        assert_eq!(reader.read_bit(), None); // exhausted
    }

    #[test]
    fn read_bits_within_byte() {
        let data = [0b1100_1010];
        let mut reader = BitReader::new(&data);
        assert_eq!(reader.read_bits(4), Some(0b1100));
        assert_eq!(reader.read_bits(4), Some(0b1010));
    }

    #[test]
    fn read_bits_cross_byte_boundary() {
        let data = [0b1100_1111, 0b0100_0000];
        let mut reader = BitReader::new(&data);
        assert_eq!(reader.read_bits(4), Some(0b1100));
        assert_eq!(reader.read_bits(6), Some(0b11_1101));
        assert_eq!(reader.read_bits(6), Some(0b00_0000));
    }

    #[test]
    fn zero_bit_read() {
        let data = [0xFF];
        let mut reader = BitReader::new(&data);
        assert_eq!(reader.read_bits(0), Some(0));
        assert_eq!(reader.remaining_bits(), 8); // no bits consumed
    }

    #[test]
    fn read_full_u32() {
        let data = [0x12, 0x34, 0x56, 0x78];
        let mut reader = BitReader::new(&data);
        assert_eq!(reader.read_bits(32), Some(0x12345678));
    }

    #[test]
    fn empty_reader() {
        let data: &[u8] = &[];
        let mut reader = BitReader::new(data);
        assert_eq!(reader.read_bit(), None);
        assert_eq!(reader.read_bits(1), None);
        assert_eq!(reader.remaining_bits(), 0);
    }

    #[test]
    fn remaining_bits_tracks_correctly() {
        let data = [0xFF, 0xFF];
        let mut reader = BitReader::new(&data);
        assert_eq!(reader.remaining_bits(), 16);
        reader.read_bits(3);
        assert_eq!(reader.remaining_bits(), 13);
        reader.read_bits(8);
        assert_eq!(reader.remaining_bits(), 5);
    }

    #[test]
    fn read_bits_insufficient_data() {
        let data = [0xFF];
        let mut reader = BitReader::new(&data);
        assert_eq!(reader.read_bits(9), None); // only 8 bits available
    }

    proptest! {
        // Read each bit individually and check it matches (byte >> (7 - i)) & 1,
        // verifying that bits are read MSB-first.
        #[test]
        fn individual_bits_match_msb_order(byte: u8) {
            let data = [byte];
            let mut reader = BitReader::new(&data);
            for i in 0u8..8 {
                let expected = (byte >> (7 - i)) & 1;
                let got = reader.read_bit().unwrap();
                prop_assert_eq!(got, expected, "bit {} of byte {:#010b}", i, byte);
            }
            prop_assert_eq!(reader.read_bit(), None);
        }

        // Read 16 individual bits from a 2-byte buffer and verify the MSB-first bit pattern.
        #[test]
        fn u16_individual_bits_match_msb_order(hi: u8, lo: u8) {
            let data = [hi, lo];
            let mut reader = BitReader::new(&data);
            for i in 0u8..8 {
                let expected = (hi >> (7 - i)) & 1;
                let got = reader.read_bit().unwrap();
                prop_assert_eq!(got, expected, "hi bit {} of {:#010b}", i, hi);
            }
            for i in 0u8..8 {
                let expected = (lo >> (7 - i)) & 1;
                let got = reader.read_bit().unwrap();
                prop_assert_eq!(got, expected, "lo bit {} of {:#010b}", i, lo);
            }
            prop_assert_eq!(reader.read_bit(), None);
        }

        #[test]
        fn split_read_equals_combined(byte: u8, split in 0u32..=8) {
            let data = [byte];
            let mut reader1 = BitReader::new(&data);
            let mut reader2 = BitReader::new(&data);

            let combined = reader1.read_bits(8).unwrap();

            let hi = reader2.read_bits(split).unwrap();
            let lo = reader2.read_bits(8 - split).unwrap();
            let reassembled = (hi << (8 - split)) | lo;

            prop_assert_eq!(combined, reassembled, "split={}", split);
        }
    }
}
