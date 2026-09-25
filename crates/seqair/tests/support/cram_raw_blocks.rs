//! Every block of a CRAM file as stored: its compression method and
//! compressed bytes, for tests and benchmarks that time or check one codec
//! on real data. Walks the containers after the 26-byte file definition —
//! the header container's blocks included — and reads each block's header
//! without decompressing it.
//!
//! Included via `#[path]`; panics on a malformed file, which only the
//! vendored test files are fed to.
//!
//! Self-contained (it parses ITF-8 itself) so that the library's own unit
//! tests can include it too.
#![allow(
    clippy::expect_used,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "test support for known-good files"
)]

/// Block compression methods, `[CRAM3]` §8.
pub const RANS_NX16: u8 = 5;
pub const ARITH: u8 = 6;
pub const FQZCOMP: u8 = 7;
pub const TOK3: u8 = 8;

/// One block as stored.
pub struct RawBlock {
    pub method: u8,
    /// The compressed payload, as handed to the codec.
    pub data: Vec<u8>,
    /// The size the block header declares for the decoded data.
    pub uncompressed_len: usize,
}

/// Reads the container and block headers' integers.
struct Cursor<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl Cursor<'_> {
    fn u8(&mut self) -> u8 {
        let b = self.bytes[self.pos];
        self.pos += 1;
        b
    }

    /// ITF-8: the count of leading 1 bits in the first byte is the count of
    /// bytes that follow; the fifth byte adds only its low 4 bits.
    fn itf8(&mut self) -> u32 {
        let first = self.u8();
        let extra = first.leading_ones().min(4);
        let mask = if extra == 4 { 0x0F } else { 0xFF >> (extra + 1) };
        let mut v = u32::from(first) & mask;
        for i in 0..extra {
            let b = u32::from(self.u8());
            v = if i == 3 { (v << 4) | (b & 0x0F) } else { (v << 8) | b };
        }
        v
    }

    /// LTF-8, skipped: 1 + the leading 1 bits of the first byte, at most 9.
    fn skip_ltf8(&mut self) {
        let first = self.u8();
        self.pos += usize::try_from(first.leading_ones()).expect("at most 8");
    }

    fn usize(&mut self) -> usize {
        usize::try_from(self.itf8()).expect("u32 fits usize")
    }
}

/// Every block of the CRAM file `file`, in file order.
pub fn raw_blocks(file: &[u8]) -> Vec<RawBlock> {
    let mut blocks = Vec::new();
    let mut cur = Cursor { bytes: file, pos: 26 };
    while cur.pos < file.len() {
        let length = u32::from_le_bytes(file[cur.pos..cur.pos + 4].try_into().expect("4 bytes"));
        cur.pos += 4;
        // Reference ID, start, span, record count.
        for _ in 0..4 {
            cur.itf8();
        }
        // Record counter, bases.
        cur.skip_ltf8();
        cur.skip_ltf8();
        let n_blocks = cur.itf8();
        for _ in 0..cur.itf8() {
            cur.itf8();
        }
        cur.pos += 4; // CRC32
        let body_end = cur.pos + usize::try_from(length).expect("u32 fits usize");
        // The header container may be padded past its blocks.
        for _ in 0..n_blocks {
            let method = cur.u8();
            let _content_type = cur.u8();
            let _content_id = cur.itf8();
            let compressed = cur.usize();
            let uncompressed_len = cur.usize();
            let data = file[cur.pos..cur.pos + compressed].to_vec();
            cur.pos += compressed + 4; // the payload, then its CRC32
            blocks.push(RawBlock { method, data, uncompressed_len });
        }
        cur.pos = body_end;
    }
    blocks
}
