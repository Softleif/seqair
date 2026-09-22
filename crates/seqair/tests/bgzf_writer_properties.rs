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
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "test code"
)]

use hegel::prelude::*;
use noodles_bgzf::VirtualPosition;
use seqair::io::BgzfWriter;
use std::collections::BTreeMap;
use std::io::{Cursor, Read as _};

/// The uncompressed payload of one BGZF block.
const BLOCK: usize = 65_536;

/// A run of writes, each one a record as far as the writer is concerned.
///
/// Sizes are drawn against the block that is currently filling, not blindly:
/// one write in four is exactly the number of bytes left in the block, which
/// is the case `r[bgzf.writer.virtual_offset]` used to get wrong and which a
/// uniform size distribution reaches about once in 65 536 tries.
#[hegel::composite]
fn arb_writes(tc: &TestCase) -> Vec<Vec<u8>> {
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
        let seed = tc.draw_silent(gs::integers::<u8>());
        writes.push(
            (0..len).map(|i| seed.wrapping_add(u8::try_from(i % 256).unwrap_or(0))).collect(),
        );

        fill = if len >= remaining { (fill + len) % BLOCK } else { fill + len };
    }
    writes
}

/// Write each chunk, taking the offset the way `BamWriter` does: ask whether
/// the block has room first, then record where the write will land.
fn write_all(writes: &[Vec<u8>], flush_first: bool) -> (Vec<u8>, Vec<u64>) {
    let mut out = Vec::new();
    let mut writer = BgzfWriter::new(&mut out);
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
#[hegel::test]
fn a_reported_offset_names_the_write_that_followed_it(tc: TestCase) {
    let writes = tc.draw(arb_writes().print_as_debug());
    let flush_first = tc.draw(gs::booleans());
    let (bgzf, offsets) = write_all(&writes, flush_first);

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
#[hegel::test]
fn an_offset_decodes_to_the_writes_own_byte_position(tc: TestCase) {
    let writes = tc.draw(arb_writes().print_as_debug());
    let flush_first = tc.draw(gs::booleans());
    let (bgzf, offsets) = write_all(&writes, flush_first);

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
/// Whatever the writes were, the bytes come back out of an independent
/// decompressor in one piece.
#[hegel::test]
fn the_stream_decompresses_to_what_was_written(tc: TestCase) {
    let writes = tc.draw(arb_writes().print_as_debug());
    let flush_first = tc.draw(gs::booleans());
    let (bgzf, _) = write_all(&writes, flush_first);

    let expected: Vec<u8> = writes.concat();
    let mut got = Vec::with_capacity(expected.len());
    noodles_bgzf::io::Reader::new(Cursor::new(&bgzf))
        .read_to_end(&mut got)
        .expect("noodles must decompress it");

    assert_eq!(got.len(), expected.len(), "decompressed length");
    assert_eq!(got, expected, "decompressed bytes");

    // Every block must stay inside BGZF's 64 KiB uncompressed limit, or htslib
    // refuses the file however well it round-trips here.
    for block in blocks(&bgzf) {
        assert!(
            block.uncompressed_len <= BLOCK,
            "a block holds {} uncompressed bytes, over the 64 KiB limit",
            block.uncompressed_len
        );
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
