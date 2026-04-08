//! Decode CRAM data series values. [`IntEncoding`], [`ByteEncoding`], and [`ByteArrayEncoding`]
//! cover all encoding kinds defined by the spec (external, Huffman, Beta, subexp, byte array
//! length/stop), and are driven by [`CompressionHeader`](super::compression_header::CompressionHeader).

use super::{bitstream::BitReader, reader::CramError, varint};
use rustc_hash::FxHashMap;

/// Integer encoding (for data series like BF, CF, RL, AP, etc.)
#[derive(Debug, Clone)]
pub enum IntEncoding {
    Null,
    External { content_id: i32 },
    Huffman(HuffmanTable),
    Beta { offset: i32, bits: u32 },
    Subexp { offset: i32, k: u32 },
    Gamma { offset: i32 },
}

/// Byte encoding (for single-byte data series like BA, QS, FC, BS)
#[derive(Debug, Clone)]
pub enum ByteEncoding {
    Null,
    External { content_id: i32 },
    Huffman(HuffmanTable),
}

/// Byte array encoding (for data series like RN, IN, SC, BB, QQ)
#[derive(Debug, Clone)]
pub enum ByteArrayEncoding {
    Null,
    External { content_id: i32 },
    ByteArrayLen { len_encoding: Box<IntEncoding>, val_encoding: Box<ByteEncoding> },
    ByteArrayStop { stop_byte: u8, content_id: i32 },
}

/// Precomputed canonical Huffman decode table.
#[derive(Debug, Clone)]
pub struct HuffmanTable {
    /// (symbol, `bit_length`) pairs sorted by (`bit_length`, symbol) for canonical assignment.
    symbols: Vec<(i32, u32)>,
    /// For single-symbol codes (`bit_length=0`), the symbol value.
    single_symbol: Option<i32>,
}

impl HuffmanTable {
    pub fn new(alphabet: &[i32], bit_lengths: &[u32]) -> Result<Self, CramError> {
        if alphabet.len() != bit_lengths.len() {
            return Err(CramError::HuffmanSizeMismatch {
                alphabet_size: alphabet.len(),
                bit_lengths_size: bit_lengths.len(),
            });
        }

        if alphabet.is_empty() {
            return Ok(Self { symbols: Vec::new(), single_symbol: None });
        }

        // Check for single-symbol code (bit_length=0)
        if let Some(&single) = alphabet.first().filter(|_| alphabet.len() == 1) {
            return Ok(Self { symbols: vec![(single, 0)], single_symbol: Some(single) });
        }

        // Build (symbol, bit_length) pairs and sort canonically
        let mut pairs: Vec<(i32, u32)> =
            alphabet.iter().zip(bit_lengths.iter()).map(|(&sym, &bl)| (sym, bl)).collect();
        pairs.sort_by(|a, b| a.1.cmp(&b.1).then(a.0.cmp(&b.0)));

        Ok(Self { symbols: pairs, single_symbol: None })
    }

    /// Decode a single symbol from the bit stream.
    ///
    /// Uses canonical Huffman decoding: read one bit at a time, tracking the
    /// current code value against the canonical code for each bit length level.
    pub fn decode(&self, reader: &mut BitReader<'_>) -> Option<i32> {
        if let Some(sym) = self.single_symbol {
            return Some(sym);
        }

        if self.symbols.is_empty() {
            return None;
        }

        // Pre-assign canonical codes: sorted by (bit_len, symbol),
        // codes start at 0 and increment, left-shifting when length increases.
        let mut codes: Vec<u32> = Vec::with_capacity(self.symbols.len());
        let mut code = 0u32;
        let mut prev_len = self.symbols.first()?.1;

        for &(_sym, bit_len) in &self.symbols {
            if bit_len > prev_len {
                let shift = bit_len.checked_sub(prev_len)?;
                if shift >= 32 {
                    return None;
                }
                code <<= shift;
                prev_len = bit_len;
            }
            codes.push(code);
            code = code.wrapping_add(1);
        }

        // Read bits incrementally, checking codes at each length level
        let mut value = 0u32;
        let mut bits_read = 0u32;

        let max_len = self.symbols.last()?.1;
        for target_len in 1..=max_len {
            while bits_read < target_len {
                value = (value << 1) | u32::from(reader.read_bit()?);
                bits_read = bits_read.checked_add(1)?;
            }
            for (i, &(sym, bl)) in self.symbols.iter().enumerate() {
                if bl == target_len && *codes.get(i)? == value {
                    return Some(sym);
                }
            }
        }

        None
    }
}

/// Cursor into an external data block.
#[derive(Debug)]
pub struct ExternalCursor {
    data: Vec<u8>,
    pos: usize,
}

impl ExternalCursor {
    pub fn new(data: Vec<u8>) -> Self {
        Self { data, pos: 0 }
    }

    pub fn read_byte(&mut self) -> Option<u8> {
        let b = self.data.get(self.pos).copied()?;
        self.pos = self.pos.checked_add(1)?;
        Some(b)
    }

    pub fn read_itf8(&mut self) -> Option<u32> {
        let remaining = self.data.get(self.pos..)?;
        let (val, n) = varint::decode_itf8(remaining)?;
        self.pos = self.pos.checked_add(n)?;
        Some(val)
    }

    // TODO(perf): read_bytes_until and read_bytes both .to_vec() from contiguous data.
    // Three options to avoid per-call allocation:
    //
    // 1. Quick win (no API change): for ByteArrayLen with External val_encoding,
    //    add read_bytes_slice() that returns &[u8] and memcpy once, instead of
    //    N push() calls. Same for read_bytes_until → read_bytes_until_slice.
    //    NLL allows returning &[u8] from &mut self (borrows data, not mutability).
    //
    // 2. decode_into(ctx, buf: &mut Vec<u8>): callers pass a reusable scratch buffer.
    //    Still copies, but buffer capacity stabilizes after first region — zero allocs
    //    in steady state. Requires changing call sites in slice.rs.
    //
    // 3. True zero-copy: change pos to Cell<usize>, all cursor methods take &self,
    //    get_external returns &ExternalCursor, decode returns Cow<'_, [u8]>.
    //    Problem: read_name is stored across many decode() calls, so Cow borrows ctx
    //    and blocks subsequent decodes. Would need split borrows or a different caller
    //    pattern (immediate .to_vec() for read_name, Cow for short-lived results).

    pub fn read_bytes_until(&mut self, stop: u8) -> Option<Vec<u8>> {
        let start = self.pos;
        while self.pos < self.data.len() {
            if *self.data.get(self.pos)? == stop {
                let result = self.data.get(start..self.pos)?.to_vec();
                self.pos = self.pos.checked_add(1)?; // skip stop byte
                return Some(result);
            }
            self.pos = self.pos.checked_add(1)?;
        }
        None
    }

    pub fn read_bytes(&mut self, n: usize) -> Option<Vec<u8>> {
        let end = self.pos.checked_add(n)?;
        let result = self.data.get(self.pos..end)?.to_vec();
        self.pos = end;
        Some(result)
    }

    pub fn into_data(self) -> Vec<u8> {
        self.data
    }

    pub fn remaining(&self) -> usize {
        self.data.len().saturating_sub(self.pos)
    }
}

/// Decode context holding core bit reader and external block cursors.
pub struct DecodeContext<'a> {
    pub core: BitReader<'a>,
    pub external: FxHashMap<i32, ExternalCursor>,
}

impl<'a> DecodeContext<'a> {
    pub fn new(core_data: &'a [u8], external_blocks: FxHashMap<i32, ExternalCursor>) -> Self {
        Self { core: BitReader::new(core_data), external: external_blocks }
    }

    fn get_external(&mut self, content_id: i32) -> Result<&mut ExternalCursor, CramError> {
        self.external.get_mut(&content_id).ok_or(CramError::ExternalBlockNotFound { content_id })
    }
}

impl IntEncoding {
    /// Decode an integer value from the decode context.
    pub fn decode(&self, ctx: &mut DecodeContext<'_>) -> Result<i32, CramError> {
        match self {
            // r[impl cram.encoding.null]
            Self::Null => Ok(0),
            // r[impl cram.encoding.external]
            Self::External { content_id } => {
                let cursor = ctx.get_external(*content_id)?;
                let val =
                    cursor.read_itf8().ok_or(CramError::Truncated { context: "external int" })?;
                Ok(val.cast_signed())
            }
            // r[impl cram.encoding.huffman]
            Self::Huffman(table) => {
                table.decode(&mut ctx.core).ok_or(CramError::Truncated { context: "huffman int" })
            }
            // r[impl cram.encoding.beta]
            Self::Beta { offset, bits } => {
                let raw = ctx
                    .core
                    .read_bits(*bits)
                    .ok_or(CramError::Truncated { context: "beta int" })?;
                raw.cast_signed()
                    .checked_sub(*offset)
                    .ok_or(CramError::Truncated { context: "beta int offset overflow" })
            }
            // r[impl cram.encoding.subexp]
            Self::Subexp { offset, k } => {
                let val = decode_subexp(&mut ctx.core, *k)
                    .ok_or(CramError::Truncated { context: "subexp int" })?;
                val.checked_sub(*offset)
                    .ok_or(CramError::Truncated { context: "subexp int offset overflow" })
            }
            // r[impl cram.encoding.gamma]
            Self::Gamma { offset } => {
                let val = decode_gamma(&mut ctx.core)
                    .ok_or(CramError::Truncated { context: "gamma int" })?;
                val.checked_sub(*offset)
                    .ok_or(CramError::Truncated { context: "gamma int offset overflow" })
            }
        }
    }

    /// Parse an encoding descriptor from a byte cursor.
    pub fn parse(cursor: &mut &[u8]) -> Result<Self, CramError> {
        let encoding_id = varint::read_itf8_from(cursor)
            .ok_or(CramError::Truncated { context: "encoding id" })?
            .cast_signed();
        let param_len = varint::read_itf8_from(cursor)
            .ok_or(CramError::Truncated { context: "encoding param length" })?
            as usize;

        let params =
            cursor.get(..param_len).ok_or(CramError::Truncated { context: "encoding params" })?;
        let mut pcur: &[u8] = params;

        let result = match encoding_id {
            0 => Ok(Self::Null),
            1 => {
                let content_id = varint::read_itf8_from(&mut pcur)
                    .ok_or(CramError::Truncated { context: "external content_id" })?
                    .cast_signed();
                Ok(Self::External { content_id })
            }
            3 => {
                let (alphabet, bit_lengths) = parse_huffman_params(&mut pcur)?;
                let table = HuffmanTable::new(&alphabet, &bit_lengths)?;
                Ok(Self::Huffman(table))
            }
            6 => {
                let offset = varint::read_itf8_from(&mut pcur)
                    .ok_or(CramError::Truncated { context: "beta offset" })?
                    .cast_signed();
                let bits = varint::read_itf8_from(&mut pcur)
                    .ok_or(CramError::Truncated { context: "beta bits" })?;
                Ok(Self::Beta { offset, bits })
            }
            7 => {
                let offset = varint::read_itf8_from(&mut pcur)
                    .ok_or(CramError::Truncated { context: "subexp offset" })?
                    .cast_signed();
                let k = varint::read_itf8_from(&mut pcur)
                    .ok_or(CramError::Truncated { context: "subexp k" })?;
                Ok(Self::Subexp { offset, k })
            }
            9 => {
                let offset = varint::read_itf8_from(&mut pcur)
                    .ok_or(CramError::Truncated { context: "gamma offset" })?
                    .cast_signed();
                Ok(Self::Gamma { offset })
            }
            _ => Err(CramError::UnsupportedEncoding { encoding_id }),
        };

        *cursor = cursor
            .get(param_len..)
            .ok_or(CramError::Truncated { context: "advance past encoding params" })?;
        result
    }
}

impl ByteEncoding {
    pub fn decode(&self, ctx: &mut DecodeContext<'_>) -> Result<u8, CramError> {
        match self {
            Self::Null => Ok(0),
            Self::External { content_id } => {
                let cursor = ctx.get_external(*content_id)?;
                cursor.read_byte().ok_or(CramError::Truncated { context: "external byte" })
            }
            Self::Huffman(table) => {
                let val = table
                    .decode(&mut ctx.core)
                    .ok_or(CramError::Truncated { context: "huffman byte" })?;
                #[expect(
                    clippy::cast_possible_truncation,
                    clippy::cast_sign_loss,
                    reason = "Huffman byte encoding returns i32 but values are 0..=255 for byte streams"
                )]
                Ok(val as u8)
            }
        }
    }

    pub fn parse(cursor: &mut &[u8]) -> Result<Self, CramError> {
        let encoding_id = varint::read_itf8_from(cursor)
            .ok_or(CramError::Truncated { context: "encoding id" })?
            .cast_signed();
        let param_len = varint::read_itf8_from(cursor)
            .ok_or(CramError::Truncated { context: "encoding param length" })?
            as usize;

        let params =
            cursor.get(..param_len).ok_or(CramError::Truncated { context: "encoding params" })?;
        let mut pcur: &[u8] = params;

        let result = match encoding_id {
            0 => Ok(Self::Null),
            1 => {
                let content_id = varint::read_itf8_from(&mut pcur)
                    .ok_or(CramError::Truncated { context: "external content_id" })?
                    .cast_signed();
                Ok(Self::External { content_id })
            }
            3 => {
                let (alphabet, bit_lengths) = parse_huffman_params(&mut pcur)?;
                let table = HuffmanTable::new(&alphabet, &bit_lengths)?;
                Ok(Self::Huffman(table))
            }
            _ => Err(CramError::UnsupportedEncoding { encoding_id }),
        };

        *cursor = cursor
            .get(param_len..)
            .ok_or(CramError::Truncated { context: "advance past encoding params" })?;
        result
    }
}

impl ByteArrayEncoding {
    pub fn decode(&self, ctx: &mut DecodeContext<'_>) -> Result<Vec<u8>, CramError> {
        match self {
            Self::Null => Ok(Vec::new()),
            // r[impl cram.encoding.external]
            Self::External { content_id } => {
                // For byte arrays with external encoding, the entire array is in the external block.
                // The length is determined by the caller or context.
                // This is typically used with a known length from elsewhere.
                Err(CramError::ExternalByteArrayNeedsLength { content_id: *content_id })
            }
            // r[impl cram.encoding.byte_array_len]
            Self::ByteArrayLen { len_encoding, val_encoding } => {
                let len_i32 = len_encoding.decode(ctx)?;
                let len = usize::try_from(len_i32)
                    .map_err(|_| super::reader::CramError::InvalidLength { value: len_i32 })?;
                super::reader::check_alloc_size(len, "byte array length")?;
                let mut result = Vec::with_capacity(len);
                for _ in 0..len {
                    result.push(val_encoding.decode(ctx)?);
                }
                Ok(result)
            }
            // r[impl cram.encoding.byte_array_stop]
            Self::ByteArrayStop { stop_byte, content_id } => {
                let cursor = ctx.get_external(*content_id)?;
                cursor
                    .read_bytes_until(*stop_byte)
                    .ok_or(CramError::Truncated { context: "byte array stop" })
            }
        }
    }

    pub fn parse(cursor: &mut &[u8]) -> Result<Self, CramError> {
        let encoding_id = varint::read_itf8_from(cursor)
            .ok_or(CramError::Truncated { context: "encoding id" })?
            .cast_signed();
        let param_len = varint::read_itf8_from(cursor)
            .ok_or(CramError::Truncated { context: "encoding param length" })?
            as usize;

        let params =
            cursor.get(..param_len).ok_or(CramError::Truncated { context: "encoding params" })?;
        let mut pcur: &[u8] = params;

        let result = match encoding_id {
            0 => Ok(Self::Null),
            1 => {
                let content_id = varint::read_itf8_from(&mut pcur)
                    .ok_or(CramError::Truncated { context: "external content_id" })?
                    .cast_signed();
                Ok(Self::External { content_id })
            }
            4 => {
                let len_enc = IntEncoding::parse(&mut pcur)?;
                let val_enc = ByteEncoding::parse(&mut pcur)?;
                Ok(Self::ByteArrayLen {
                    len_encoding: Box::new(len_enc),
                    val_encoding: Box::new(val_enc),
                })
            }
            5 => {
                let stop = *pcur
                    .first()
                    .ok_or(CramError::Truncated { context: "byte array stop byte" })?;
                pcur = pcur.get(1..).ok_or(CramError::Truncated { context: "byte array stop" })?;
                let content_id = varint::read_itf8_from(&mut pcur)
                    .ok_or(CramError::Truncated { context: "byte array stop content_id" })?
                    .cast_signed();
                Ok(Self::ByteArrayStop { stop_byte: stop, content_id })
            }
            _ => Err(CramError::UnsupportedEncoding { encoding_id }),
        };

        *cursor = cursor
            .get(param_len..)
            .ok_or(CramError::Truncated { context: "advance past encoding params" })?;
        result
    }
}

fn parse_huffman_params(cursor: &mut &[u8]) -> Result<(Vec<i32>, Vec<u32>), CramError> {
    let alpha_count = varint::read_itf8_from(cursor)
        .ok_or(CramError::Truncated { context: "huffman alphabet count" })?;
    let alpha_count_usize = alpha_count as usize;
    super::reader::check_alloc_size(alpha_count_usize.saturating_mul(4), "huffman alphabet")?;
    let mut alphabet = Vec::with_capacity(alpha_count_usize);
    for _ in 0..alpha_count {
        let sym = varint::read_itf8_from(cursor)
            .ok_or(CramError::Truncated { context: "huffman alphabet symbol" })?;
        alphabet.push(sym.cast_signed());
    }

    let bl_count = varint::read_itf8_from(cursor)
        .ok_or(CramError::Truncated { context: "huffman bit length count" })?;
    let bl_count_usize = bl_count as usize;
    super::reader::check_alloc_size(bl_count_usize.saturating_mul(4), "huffman bit lengths")?;
    let mut bit_lengths = Vec::with_capacity(bl_count_usize);
    for _ in 0..bl_count {
        let bl = varint::read_itf8_from(cursor)
            .ok_or(CramError::Truncated { context: "huffman bit length" })?;
        bit_lengths.push(bl);
    }

    Ok((alphabet, bit_lengths))
}

/// Decode Elias gamma code from a bit stream.
fn decode_gamma(reader: &mut BitReader<'_>) -> Option<i32> {
    let mut n = 0u32;
    while reader.read_bit()? == 0 {
        n = n.checked_add(1)?;
    }
    if n == 0 {
        return Some(0);
    }
    let val = reader.read_bits(n)?;
    let combined = (1u32.checked_shl(n)? | val).cast_signed();
    combined.checked_sub(1)
}

/// Decode sub-exponential code from a bit stream.
fn decode_subexp(reader: &mut BitReader<'_>, k: u32) -> Option<i32> {
    let mut n = 0u32;
    while reader.read_bit()? == 1 {
        n = n.checked_add(1)?;
    }
    if n == 0 {
        let val = reader.read_bits(k)?;
        return Some(val.cast_signed());
    }
    let bits = n.checked_add(k)?.checked_sub(1)?;
    let val = reader.read_bits(bits)?;
    let base = 1u32.checked_shl(bits)?.checked_sub(1u32.checked_shl(k)?)?;
    Some(base.checked_add(val)?.cast_signed())
}

#[cfg(test)]
#[allow(
    clippy::arithmetic_side_effects,
    reason = "test-only arithmetic on bounded/proptest-generated values"
)]
mod tests {
    use super::*;

    // r[verify cram.encoding.huffman]
    #[test]
    fn huffman_single_symbol() {
        let table = HuffmanTable::new(&[42], &[0]).unwrap();
        let data = [0u8; 1]; // no bits needed
        let mut reader = BitReader::new(&data);
        assert_eq!(table.decode(&mut reader), Some(42));
        // Should not consume any bits
        assert_eq!(reader.remaining_bits(), 8);
    }

    #[test]
    fn huffman_two_symbols() {
        // A=0 (1 bit: "0"), B=1 (1 bit: "1")
        let table = HuffmanTable::new(&[0, 1], &[1, 1]).unwrap();
        // Data: 0b10 = symbol 1 then symbol 0
        let data = [0b1000_0000];
        let mut reader = BitReader::new(&data);
        assert_eq!(table.decode(&mut reader), Some(1));
        assert_eq!(table.decode(&mut reader), Some(0));
    }

    #[test]
    fn huffman_three_symbols() {
        // Symbols sorted by (bit_len, symbol): A(1bit), B(2bit), C(2bit)
        // Canonical codes: A=0, B=10, C=11
        let table = HuffmanTable::new(&[65, 66, 67], &[1, 2, 2]).unwrap();
        // Data: A(0) B(10) C(11) = 0b0_10_11_000 = 0x58
        let data = [0b0101_1000];
        let mut reader = BitReader::new(&data);
        assert_eq!(table.decode(&mut reader), Some(65)); // A
        assert_eq!(table.decode(&mut reader), Some(66)); // B
        assert_eq!(table.decode(&mut reader), Some(67)); // C
    }

    #[test]
    fn huffman_empty() {
        let table = HuffmanTable::new(&[], &[]).unwrap();
        let data = [0u8; 1];
        let mut reader = BitReader::new(&data);
        assert_eq!(table.decode(&mut reader), None);
    }

    #[test]
    fn external_cursor_read_byte() {
        let mut cursor = ExternalCursor::new(vec![0x41, 0x42, 0x43]);
        assert_eq!(cursor.read_byte(), Some(0x41));
        assert_eq!(cursor.read_byte(), Some(0x42));
        assert_eq!(cursor.read_byte(), Some(0x43));
        assert_eq!(cursor.read_byte(), None);
    }

    #[test]
    fn external_cursor_read_itf8() {
        let mut cursor = ExternalCursor::new(vec![0x80, 0x80, 0x05]);
        assert_eq!(cursor.read_itf8(), Some(128));
        assert_eq!(cursor.read_itf8(), Some(5));
        assert_eq!(cursor.read_itf8(), None);
    }

    #[test]
    fn external_cursor_read_bytes_until() {
        let mut cursor = ExternalCursor::new(b"hello\x00world\x00".to_vec());
        assert_eq!(cursor.read_bytes_until(0), Some(b"hello".to_vec()));
        assert_eq!(cursor.read_bytes_until(0), Some(b"world".to_vec()));
        assert_eq!(cursor.read_bytes_until(0), None);
    }

    // r[verify cram.encoding.beta]
    #[test]
    fn beta_encoding_decode() {
        let data = [0b10110000]; // bits: 1011 = 11
        let mut ctx = DecodeContext::new(&data, FxHashMap::default());
        let enc = IntEncoding::Beta { offset: 0, bits: 4 };
        assert_eq!(enc.decode(&mut ctx).unwrap(), 11);
    }

    #[test]
    fn beta_encoding_with_offset() {
        let data = [0b10110000]; // bits: 1011 = 11, minus offset 5 = 6
        let mut ctx = DecodeContext::new(&data, FxHashMap::default());
        let enc = IntEncoding::Beta { offset: 5, bits: 4 };
        assert_eq!(enc.decode(&mut ctx).unwrap(), 6);
    }

    // r[verify cram.encoding.external]
    #[test]
    fn external_int_decode() {
        let mut external = FxHashMap::default();
        external.insert(7, ExternalCursor::new(vec![42])); // ITF8(42)
        let mut ctx = DecodeContext::new(&[], external);
        let enc = IntEncoding::External { content_id: 7 };
        assert_eq!(enc.decode(&mut ctx).unwrap(), 42);
    }

    #[test]
    fn external_byte_decode() {
        let mut external = FxHashMap::default();
        external.insert(3, ExternalCursor::new(vec![0xAB]));
        let mut ctx = DecodeContext::new(&[], external);
        let enc = ByteEncoding::External { content_id: 3 };
        assert_eq!(enc.decode(&mut ctx).unwrap(), 0xAB);
    }

    // r[verify cram.encoding.byte_array_stop]
    #[test]
    fn byte_array_stop_decode() {
        let mut external = FxHashMap::default();
        external.insert(5, ExternalCursor::new(b"test\x00".to_vec()));
        let mut ctx = DecodeContext::new(&[], external);
        let enc = ByteArrayEncoding::ByteArrayStop { stop_byte: 0, content_id: 5 };
        assert_eq!(enc.decode(&mut ctx).unwrap(), b"test");
    }

    // r[verify cram.encoding.byte_array_len]
    #[test]
    fn byte_array_len_decode() {
        let mut external = FxHashMap::default();
        // Length comes from core Huffman (single symbol = 3)
        // Values come from external block
        external.insert(1, ExternalCursor::new(vec![0x41, 0x42, 0x43]));
        let mut ctx = DecodeContext::new(&[], external);

        let enc = ByteArrayEncoding::ByteArrayLen {
            len_encoding: Box::new(IntEncoding::Huffman(HuffmanTable::new(&[3], &[0]).unwrap())),
            val_encoding: Box::new(ByteEncoding::External { content_id: 1 }),
        };
        assert_eq!(enc.decode(&mut ctx).unwrap(), b"ABC");
    }

    // r[verify cram.encoding.null]
    #[test]
    fn null_encodings_return_defaults() {
        let mut ctx = DecodeContext::new(&[], FxHashMap::default());
        assert_eq!(IntEncoding::Null.decode(&mut ctx).unwrap(), 0);
        assert_eq!(ByteEncoding::Null.decode(&mut ctx).unwrap(), 0);
        assert_eq!(ByteArrayEncoding::Null.decode(&mut ctx).unwrap(), Vec::<u8>::new());
    }

    #[test]
    fn huffman_size_mismatch_returns_error() {
        // alphabet has 3 elements, bit_lengths has 2 — should return HuffmanSizeMismatch
        let err = HuffmanTable::new(&[1, 2, 3], &[1, 1]).unwrap_err();
        assert!(matches!(
            err,
            CramError::HuffmanSizeMismatch { alphabet_size: 3, bit_lengths_size: 2 }
        ));
    }

    #[test]
    fn huffman_size_mismatch_reversed() {
        let err = HuffmanTable::new(&[1], &[1, 2]).unwrap_err();
        assert!(matches!(
            err,
            CramError::HuffmanSizeMismatch { alphabet_size: 1, bit_lengths_size: 2 }
        ));
    }

    #[test]
    fn external_byte_array_needs_length_returned() {
        // ByteArrayEncoding::External always returns ExternalByteArrayNeedsLength
        let enc = ByteArrayEncoding::External { content_id: 42 };
        let mut ctx = DecodeContext::new(&[], FxHashMap::default());
        let err = enc.decode(&mut ctx).unwrap_err();
        assert!(matches!(err, CramError::ExternalByteArrayNeedsLength { content_id: 42 }));
    }

    proptest::proptest! {
        // Verify Beta decoding with varying bit widths and offsets: decoded value must
        // equal the raw bit pattern minus the offset, as per the CRAM spec.
        #[test]
        fn beta_decode_with_varying_params(
            bits in 1u32..=16,
            offset in -100i32..=100,
            val in 0u32..=(u32::MAX),
        ) {
            // Clamp val to the valid range for this bit width
            let max_val = (1u32 << bits) - 1;
            let val = val % (max_val + 1);

            // Pack `val` MSB-first into 3 bytes (enough for up to 16 bits)
            let shift = 24u32.saturating_sub(bits);
            let packed = val << shift;
            let data = [(packed >> 16) as u8, (packed >> 8) as u8, packed as u8];

            let mut ctx = DecodeContext::new(&data, FxHashMap::default());
            let enc = IntEncoding::Beta { offset, bits };
            let decoded = enc.decode(&mut ctx).unwrap();
            proptest::prop_assert_eq!(decoded, val as i32 - offset);
        }

        #[test]
        fn external_byte_roundtrip(bytes in proptest::collection::vec(proptest::prelude::any::<u8>(), 1..32)) {
            let mut external = FxHashMap::default();
            external.insert(0, ExternalCursor::new(bytes.clone()));
            let mut ctx = DecodeContext::new(&[], external);
            let enc = ByteEncoding::External { content_id: 0 };
            for &expected in &bytes {
                let got = enc.decode(&mut ctx).unwrap();
                proptest::prop_assert_eq!(got, expected);
            }
        }
    }
}
