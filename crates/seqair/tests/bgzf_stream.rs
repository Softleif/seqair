//! `BgzfStream` (front-to-back BGZF decoding) and the block framing it shares
//! with `BgzfReader` and `RegionBuf`.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "test code"
)]
#![allow(clippy::cast_possible_truncation, reason = "test code with known small values")]

use seqair::bam::BgzfError;
use seqair::bam::bgzf::{BgzfReader, BgzfStream, EofMarker, VirtualOffset};
use seqair::bam::index::Chunk;
use seqair::bam::region_buf::RegionBuf;
use seqair::io::BgzfWriter;
use std::io::Cursor;

/// The 28-byte BGZF EOF marker block (`[SAM1] §4.1.2`).
const EOF_BLOCK: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0, 0x42, 0x43, 0x02, 0, 0x1b, 0, 0x03, 0, 0,
    0, 0, 0, 0, 0, 0, 0,
];

/// One BGZF block holding `data` (at most 65280 bytes), as seqair writes it.
fn block(data: &[u8]) -> Vec<u8> {
    let mut writer = BgzfWriter::new(Vec::new());
    writer.write_all(data).unwrap();
    let mut bytes = writer.finish().unwrap();
    assert!(bytes.ends_with(&EOF_BLOCK));
    bytes.truncate(bytes.len() - EOF_BLOCK.len());
    bytes
}

/// Every byte of `stream` through the typed API, or the first error.
fn stream_all(bytes: &[u8], threads: usize, eof: EofMarker) -> Result<Vec<u8>, BgzfError> {
    let mut stream =
        BgzfStream::builder(Cursor::new(bytes)).threads(threads).eof_marker(eof).build()?;
    let mut out = Vec::new();
    loop {
        let chunk = stream.fill()?;
        if chunk.is_empty() {
            return Ok(out);
        }
        let n = chunk.len();
        out.extend_from_slice(chunk);
        stream.advance(n);
    }
}

fn reader_all(bytes: &[u8]) -> Result<Vec<u8>, BgzfError> {
    let mut reader = BgzfReader::from_reader(Cursor::new(bytes.to_vec()));
    let mut out = Vec::new();
    reader.read_to_end(&mut out)?;
    Ok(out)
}

/// Every byte `RegionBuf` reads from one range spanning the whole input.
fn region_buf_all(bytes: &[u8]) -> Result<Vec<u8>, BgzfError> {
    let mut cursor = Cursor::new(bytes.to_vec());
    let chunk =
        Chunk { begin: VirtualOffset::new(0, 0), end: VirtualOffset::new(bytes.len() as u64, 0) };
    let mut buf = RegionBuf::new(&mut cursor, &[chunk])?;
    let mut out = Vec::new();
    loop {
        let (_, data) = buf.fill_buf()?;
        if data.is_empty() {
            return Ok(out);
        }
        let n = data.len();
        out.extend_from_slice(data);
        buf.consume(n);
    }
}

// r[verify bgzf.empty_block]
// r[verify region_buf.empty_block]
/// A corrupt empty block — the EOF marker or one mid-file — fails like any
/// other corrupt block: htslib inflates it and checks its CRC, so a reader
/// that skipped it unchecked accepted files htslib rejects.
#[test]
fn a_corrupt_empty_block_is_an_error() {
    let data = block(b"hello, world");
    // Payload (18, 19), CRC (20..24) and ISIZE (24..28) of the empty block.
    for at in 18..28 {
        for (name, file) in [
            ("EOF marker", [data.clone(), EOF_BLOCK.to_vec()].concat()),
            ("mid-file", [EOF_BLOCK.to_vec(), data.clone(), EOF_BLOCK.to_vec()].concat()),
        ] {
            let corrupt_at = if name == "mid-file" { at } else { data.len() + at };
            let mut file = file;
            file[corrupt_at] ^= 0x01;
            for threads in [0, 2] {
                assert!(
                    stream_all(&file, threads, EofMarker::Optional).is_err(),
                    "{name}: stream (threads {threads}) took byte {at} flipped"
                );
            }
            assert!(reader_all(&file).is_err(), "{name}: BgzfReader took byte {at} flipped");
            assert!(region_buf_all(&file).is_err(), "{name}: RegionBuf took byte {at} flipped");
        }
    }
}

// r[verify bgzf.decompression+2]
/// A block whose ISIZE promises more bytes than its payload inflates to is
/// an error in every reader — before, `RegionBuf` truncated to what it got,
/// and the CRC (over the real bytes) still matched.
#[test]
fn an_isize_larger_than_the_payload_is_an_error() {
    let mut data = block(b"hello, world");
    let isize_at = data.len() - 4;
    data[isize_at] += 1;
    let file = [data, EOF_BLOCK.to_vec()].concat();
    assert!(matches!(
        stream_all(&file, 0, EofMarker::Required),
        Err(BgzfError::SizeMismatch { .. })
    ));
    assert!(matches!(reader_all(&file), Err(BgzfError::SizeMismatch { .. })));
    assert!(matches!(region_buf_all(&file), Err(BgzfError::SizeMismatch { .. })));
}
