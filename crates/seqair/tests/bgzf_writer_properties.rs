//! What `BgzfWriter` promises about the offsets it hands out, and about the
//! bytes it writes, on write sequences nobody wrote by hand.
//!
//! `r[bgzf.writer.virtual_offset]` is the rule a co-produced index rests on:
//! the offset taken before a record is written must name that record's first
//! byte. It has been wrong once — a write that filled a block exactly reported
//! the block's *last* byte, and every index offset taken there pointed into the
//! middle of a record. The test that caught it sweeps 96 record sizes of one
//! fixed shape through the BCF writer. This asks the BGZF writer directly, for
//! any sequence of writes, and deliberately keeps landing on the boundary that
//! broke it.
//!
//! noodles reads the result, so the offsets are checked by an implementation
//! that had no part in producing them.
//!
//! The bytes and the compression level are drawn too. A block of data that
//! does not compress is *stored*, a little larger than it went in, and the
//! writer once filled blocks to 64 KiB — which stored came to more than a BGZF
//! block may hold, so every level-0 file and every block of noise failed to
//! write (`r[bgzf.writer.block_size]`). Text-like bytes at the default level
//! never get near that, which is all these properties used to write.
//!
//! The case counts are pinned. Each case deflates up to 40 chunks, some of them
//! larger than a block, which is dear enough that the suite default would make
//! this file most of the wall time; and the generator reaches the boundary
//! these properties exist for in half the cases, so a hundred and fifty of them
//! land on it far more often than the 96-size sweep they supersede. For a
//! deeper run, `HEGEL_TEST_CASES=N` overrides the per-test count.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    reason = "test code"
)]

use hegel::prelude::*;
use noodles_bgzf::VirtualPosition;
use seqair::io::BgzfWriter;
use std::collections::BTreeMap;
use std::io::{Cursor, Read as _};

/// The uncompressed payload of a full BGZF block: htslib's `BGZF_BLOCK_SIZE`,
/// which `r[bgzf.writer.block_size]` adopts. The generator aims writes at this
/// boundary, so it has to be the writer's, not the format's 64 KiB ceiling.
const BLOCK: usize = 0xff00;

/// What the bytes of a case look like to DEFLATE.
#[derive(Debug, Clone, Copy)]
enum Content {
    /// A byte ramp: compresses to almost nothing at any level above 0.
    Ramp,
    /// xorshift output: does not compress at all, so every block is stored.
    Noise,
    /// Each write picks one, so blocks are part noise.
    Mixed,
}

/// A libdeflate level, 0 (store only) through 12.
fn arb_level() -> impl PrintableGenerator<i32> {
    gs::integers::<i32>().min_value(0).max_value(12)
}

/// A run of writes, each one a record as far as the writer is concerned.
///
/// Sizes are drawn against the block that is currently filling, not blindly:
/// one write in four is exactly the number of bytes left in the block, which
/// is the case `r[bgzf.writer.virtual_offset]` used to get wrong and which a
/// uniform size distribution reaches about once in 65 536 tries.
#[hegel::composite]
fn arb_writes(tc: &TestCase) -> Vec<Vec<u8>> {
    let content =
        tc.draw_silent(gs::sampled_from(&[Content::Ramp, Content::Noise, Content::Mixed]));
    let n = tc.draw_silent(gs::integers::<usize>().min_value(1).max_value(40));
    let mut writes = Vec::with_capacity(n);
    let mut fill = 0usize;

    for _ in 0..n {
        let remaining = BLOCK - fill;
        let len = if tc.draw_silent(gs::weighted_booleans(0.25)) {
            // Exactly finish the block.
            remaining
        } else if tc.draw_silent(gs::weighted_booleans(0.15)) {
            // Land one byte either side of the boundary, the neighbours of the
            // case above.
            let off = tc.draw_silent(gs::integers::<usize>().min_value(1).max_value(2));
            remaining.saturating_sub(off).max(1)
        } else if tc.draw_silent(gs::weighted_booleans(0.1)) {
            // A record larger than a whole block, which must be split across
            // blocks rather than refused.
            tc.draw_silent(gs::integers::<usize>().min_value(BLOCK).max_value(BLOCK * 2))
        } else {
            tc.draw_silent(gs::integers::<usize>().min_value(1).max_value(4_000))
        };

        // Bytes that differ between writes, so landing on the wrong record is
        // visible rather than a coincidence of identical content.
        let seed = tc.draw_silent(gs::integers::<u64>());
        let noise = match content {
            Content::Ramp => false,
            Content::Noise => true,
            Content::Mixed => tc.draw_silent(gs::booleans()),
        };
        writes.push(if noise { noise_bytes(seed, len) } else { ramp_bytes(seed as u8, len) });

        fill = if len >= remaining { (fill + len) % BLOCK } else { fill + len };
    }
    writes
}

fn ramp_bytes(seed: u8, len: usize) -> Vec<u8> {
    (0..len).map(|i| seed.wrapping_add(i as u8)).collect()
}

fn noise_bytes(seed: u64, len: usize) -> Vec<u8> {
    // xorshift64 has no zero state; any other seed gives a full-period stream.
    let mut x = seed | 1;
    (0..len)
        .map(|_| {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            (x >> 32) as u8
        })
        .collect()
}

/// Write each chunk, taking the offset the way `BamWriter` does: ask whether
/// the block has room first, then record where the write will land.
fn write_all(writes: &[Vec<u8>], level: i32, flush_first: bool) -> (Vec<u8>, Vec<u64>) {
    let mut out = Vec::new();
    let mut writer = BgzfWriter::with_compression_level(&mut out, level);
    let mut offsets = Vec::with_capacity(writes.len());

    for chunk in writes {
        if flush_first {
            writer.flush_if_needed(chunk.len()).expect("flush");
        }
        offsets.push(writer.virtual_offset().0);
        writer.write_all(chunk).expect("write");
    }
    writer.finish().expect("finish");
    (out, offsets)
}

// r[verify bgzf.writer.virtual_offset]
/// Seeking a reader to the offset the writer reported before a write must land
/// on that write's first byte — for every write, whatever sizes preceded it.
#[hegel::test(test_cases = 150)]
fn a_reported_offset_names_the_write_that_followed_it(tc: TestCase) {
    let writes = tc.draw(arb_writes().print_as_debug());
    let level = tc.draw(arb_level());
    let flush_first = tc.draw(gs::booleans());
    let (bgzf, offsets) = write_all(&writes, level, flush_first);

    let mut reader = noodles_bgzf::io::Reader::new(Cursor::new(&bgzf));
    for (i, (offset, chunk)) in offsets.iter().zip(&writes).enumerate() {
        reader.seek(VirtualPosition::from(*offset)).expect("seek");
        let mut got = vec![0u8; chunk.len()];
        reader.read_exact(&mut got).expect("read");
        assert_eq!(
            got,
            *chunk,
            "write {i} of {}: offset {offset:#x} (block {}, within {}) does not name it",
            writes.len(),
            offset >> 16,
            offset & 0xFFFF,
        );
    }

    // The case the rule exists for. `flush_if_needed` moves a write that would
    // not fit to the next block, so only the un-flushed run reaches it.
    if !flush_first && offsets.iter().any(|o| o & 0xFFFF == 0 && *o != 0) {
        tc.event("a write began exactly on a block boundary");
    }
    tc.event_value("blocks", blocks(&bgzf).len() as f64);
}

// r[verify bgzf.writer.virtual_offset]
/// Decoding the offset by hand — find its block in the file, add its
/// within-block part to that block's uncompressed start — must give the byte
/// position of the write. The failure mode that corrupted indexes produced a
/// *valid-looking* offset one byte short of the block's end, which a seek
/// happily accepts; this is the arithmetic that says which byte it really is.
#[hegel::test(test_cases = 150)]
fn an_offset_decodes_to_the_writes_own_byte_position(tc: TestCase) {
    let writes = tc.draw(arb_writes().print_as_debug());
    let level = tc.draw(arb_level());
    let flush_first = tc.draw(gs::booleans());
    let (bgzf, offsets) = write_all(&writes, level, flush_first);

    // Where each block starts, in the file and in the uncompressed stream.
    // Read from the gzip trailers, so nothing here comes from the writer's
    // own bookkeeping.
    let mut uncompressed_start: BTreeMap<u64, u64> = BTreeMap::new();
    let mut at = 0u64;
    for block in blocks(&bgzf) {
        uncompressed_start.insert(block.file_offset, at);
        at += block.uncompressed_len as u64;
    }

    let mut want = 0u64;
    for (i, (offset, chunk)) in offsets.iter().zip(&writes).enumerate() {
        let block_offset = offset >> 16;
        let within = offset & 0xFFFF;
        let start = *uncompressed_start.get(&block_offset).unwrap_or_else(|| {
            panic!(
                "write {i}: offset names block {block_offset:#x}, which is not a block in the file"
            )
        });
        assert_eq!(
            start + within,
            want,
            "write {i}: offset {offset:#x} decodes to byte {}, but the write starts at {want}",
            start + within,
        );
        want += chunk.len() as u64;
    }
}

// r[verify bgzf.writer.buffer]
// r[verify bgzf.writer.block_size]
/// Whatever the writes were, at whatever level, the writer accepts them and
/// the bytes come back out of an independent decompressor in one piece.
#[hegel::test(test_cases = 150)]
fn the_stream_decompresses_to_what_was_written(tc: TestCase) {
    let writes = tc.draw(arb_writes().print_as_debug());
    let level = tc.draw(arb_level());
    let flush_first = tc.draw(gs::booleans());
    let (bgzf, _) = write_all(&writes, level, flush_first);

    let expected: Vec<u8> = writes.concat();
    let mut got = Vec::with_capacity(expected.len());
    noodles_bgzf::io::Reader::new(Cursor::new(&bgzf))
        .read_to_end(&mut got)
        .expect("noodles must decompress it");

    assert_eq!(got.len(), expected.len(), "decompressed length");
    assert_eq!(got, expected, "decompressed bytes");

    // Every block stays within htslib's block size, and `blocks` walking the
    // whole file by BSIZE means every block's framing is sound.
    let found = blocks(&bgzf);
    for block in &found {
        assert!(
            block.uncompressed_len <= BLOCK,
            "a block holds {} uncompressed bytes, over htslib's {BLOCK}",
            block.uncompressed_len
        );
    }
    let covered: usize = found.iter().map(|b| b.uncompressed_len).sum();
    assert_eq!(covered, expected.len(), "blocks found by BSIZE cover the stream");

    let full = found.iter().filter(|b| b.uncompressed_len == BLOCK).count();
    if level == 0 && full > 0 {
        tc.event("a full block stored at level 0");
    }
}

/// Where each BGZF block starts in the file, and how many uncompressed bytes
/// it holds — read from the BSIZE subfield and the gzip trailer, so none of it
/// comes from the writer's own bookkeeping.
struct Block {
    file_offset: u64,
    uncompressed_len: usize,
}

fn blocks(bgzf: &[u8]) -> Vec<Block> {
    let mut found = Vec::new();
    let mut at = 0usize;
    while at + 18 <= bgzf.len() {
        // BSIZE is the last two bytes of the BC subfield, which BGZF requires
        // to be the first extra subfield; it holds the block's size minus one.
        let bsize = u16::from_le_bytes([bgzf[at + 16], bgzf[at + 17]]) as usize + 1;
        let end = at + bsize;
        if end > bgzf.len() || bsize < 26 {
            break;
        }
        let isize_at = end - 4;
        let uncompressed_len = u32::from_le_bytes([
            bgzf[isize_at],
            bgzf[isize_at + 1],
            bgzf[isize_at + 2],
            bgzf[isize_at + 3],
        ]) as usize;
        found.push(Block { file_offset: at as u64, uncompressed_len });
        at = end;
    }
    found
}
