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

/// What a stream yields: every byte up to the end or the first error, and
/// that error.
#[derive(Debug)]
struct Outcome {
    data: Vec<u8>,
    error: Option<BgzfError>,
}

/// Read `bytes` to the end through the typed `fill`/`advance` API.
fn stream_outcome(bytes: &[u8], threads: usize, eof: EofMarker) -> Outcome {
    let mut stream =
        BgzfStream::builder(Cursor::new(bytes)).threads(threads).eof_marker(eof).build().unwrap();
    let mut data = Vec::new();
    loop {
        match stream.fill() {
            Ok([]) => return Outcome { data, error: None },
            Ok(chunk) => {
                let n = chunk.len();
                data.extend_from_slice(chunk);
                stream.advance(n);
            }
            Err(e) => {
                assert!(
                    matches!(stream.fill(), Err(BgzfError::Poisoned)),
                    "not poisoned after {e:?}"
                );
                return Outcome { data, error: Some(e) };
            }
        }
    }
}

/// Every byte of `bytes` through the typed API, or the first error.
fn stream_all(bytes: &[u8], threads: usize, eof: EofMarker) -> Result<Vec<u8>, BgzfError> {
    let Outcome { data, error } = stream_outcome(bytes, threads, eof);
    error.map_or(Ok(data), Err)
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

// r[verify bgzf.stream.threads]
/// A block that fails to inflate is reported after every byte before it —
/// including the blocks before it in the same batch of eight, which the
/// stream used to drop with the batch.
#[test]
fn a_bad_block_is_reported_after_every_byte_before_it() {
    let blocks: Vec<Vec<u8>> = (0..20u8).map(|i| vec![i; 100 + usize::from(i)]).collect();
    for bad in [0, 1, 5, 7, 8, 13, 19] {
        let mut file = Vec::new();
        for (i, data) in blocks.iter().enumerate() {
            let mut b = block(data);
            if i == bad {
                let crc_at = b.len() - 8;
                b[crc_at] ^= 0x01;
            }
            file.extend_from_slice(&b);
        }
        file.extend_from_slice(&EOF_BLOCK);
        let before: Vec<u8> = blocks[..bad].concat();
        for threads in [0, 1, 3] {
            let outcome = stream_outcome(&file, threads, EofMarker::Required);
            assert_eq!(outcome.data, before, "bad block {bad}, threads {threads}");
            assert!(
                matches!(outcome.error, Some(BgzfError::ChecksumMismatch { .. })),
                "bad block {bad}, threads {threads}: {:?}",
                outcome.error
            );
        }
    }
}

// ── generated streams ───────────────────────────────────────────────────

#[path = "../src/pinned.rs"]
mod pinned;

use hegel::prelude::*;
use std::io::{BufRead, Read};

/// `len` bytes from `seed`: compressible runs or noise, without drawing every
/// byte through hegel (which keeps vectors short).
fn bytes_from_seed(seed: u64, len: usize, noisy: bool) -> Vec<u8> {
    let mut x = seed | 1;
    (0..len)
        .map(|i| {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            if noisy { x as u8 } else { b"ACGT\n"[(x as usize ^ (i / 7)) % 5] }
        })
        .collect()
}

/// One piece of a stream: a block with data, or an empty block.
#[derive(Debug, Clone)]
enum Piece {
    /// Data, and whether noodles (rather than seqair) compresses it.
    Data(Vec<u8>, bool),
    Empty,
}

/// A BGZF stream assembled block by block, so every block boundary is known.
#[derive(Debug, Clone)]
struct Layout {
    pieces: Vec<Piece>,
    /// Whether an EOF marker block closes the stream.
    eof_marker: bool,
}

impl Layout {
    /// The stream's bytes and the offset where each piece ends.
    fn encode(&self) -> (Vec<u8>, Vec<usize>) {
        let mut out = Vec::new();
        let mut ends = Vec::new();
        for piece in &self.pieces {
            match piece {
                Piece::Data(data, false) => out.extend_from_slice(&block(data)),
                Piece::Data(data, true) => out.extend_from_slice(&noodles_block(data)),
                Piece::Empty => out.extend_from_slice(&EOF_BLOCK),
            }
            ends.push(out.len());
        }
        if self.eof_marker {
            out.extend_from_slice(&EOF_BLOCK);
            ends.push(out.len());
        }
        (out, ends)
    }

    fn data(&self) -> Vec<u8> {
        self.data_before(self.pieces.len())
    }

    /// The data of the first `n` pieces.
    fn data_before(&self, n: usize) -> Vec<u8> {
        self.pieces
            .iter()
            .take(n)
            .filter_map(|p| match p {
                Piece::Data(data, _) => Some(data.as_slice()),
                Piece::Empty => None,
            })
            .flatten()
            .copied()
            .collect()
    }

    /// Whether the last block is empty — an EOF marker, as a reader sees it.
    fn ends_with_empty_block(&self) -> bool {
        self.eof_marker || matches!(self.pieces.last(), Some(Piece::Empty))
    }
}

/// One BGZF block holding `data`, as noodles writes it.
fn noodles_block(data: &[u8]) -> Vec<u8> {
    use std::io::Write;
    let mut writer = noodles::bgzf::io::Writer::new(Vec::new());
    writer.write_all(data).unwrap();
    let mut bytes = writer.finish().unwrap();
    assert!(bytes.ends_with(&EOF_BLOCK));
    bytes.truncate(bytes.len() - EOF_BLOCK.len());
    bytes
}

#[hegel::composite]
fn piece(tc: &TestCase, max_len: usize) -> Piece {
    if tc.draw(gs::integers::<u8>().max_value(5)) == 0 {
        return Piece::Empty;
    }
    // Mostly small blocks, so a stream spans several eight-block batches;
    // now and then a full one.
    let len = if tc.draw(gs::integers::<u8>().max_value(15)) == 0 {
        max_len
    } else {
        tc.draw(gs::integers::<usize>().min_value(1).max_value(max_len.min(2000)))
    };
    let data = bytes_from_seed(tc.draw(gs::integers::<u64>()), len, tc.draw(gs::booleans()));
    Piece::Data(data, tc.draw(gs::booleans()))
}

#[hegel::composite]
fn layout(tc: &TestCase, max_pieces: usize, max_len: usize) -> Layout {
    Layout {
        pieces: tc.draw(gs::vecs(piece(max_len)).max_size(max_pieces)),
        eof_marker: tc.draw(gs::booleans()),
    }
}

/// How a test reads a stream.
#[derive(Debug, Clone, Copy)]
enum Access {
    /// `fill`/`advance`, taking at most this many bytes at a time.
    FillAdvance(usize),
    /// `read_exact_into` in chunks of this size, then `read_up_to`.
    ReadExact(usize),
    /// `io::Read::read` into a buffer of this size.
    IoRead(usize),
    /// `io::BufRead::read_until` this delimiter.
    ReadUntil(u8),
}

#[hegel::composite]
fn access(tc: &TestCase) -> Access {
    let size = tc.draw(gs::sampled_from(vec![1usize, 3, 4, 17, 4096, 65536, 200_000]));
    match tc.draw(gs::integers::<u8>().max_value(3)) {
        0 => Access::FillAdvance(size),
        1 => Access::ReadExact(size),
        2 => Access::IoRead(size),
        _ => Access::ReadUntil(tc.draw(gs::sampled_from(vec![b'\n', 0u8, b'A']))),
    }
}

/// The `BgzfError` inside an `io::Error` from the `Read`/`BufRead` impls.
fn bgzf_error(e: std::io::Error) -> BgzfError {
    assert_eq!(e.kind(), std::io::ErrorKind::InvalidData);
    *e.into_inner().unwrap().downcast::<BgzfError>().unwrap()
}

/// Read `bytes` to the end, or to the first error, the way `access` says.
fn read_with(bytes: &[u8], threads: usize, eof: EofMarker, access: Access) -> Outcome {
    let mut stream =
        BgzfStream::builder(Cursor::new(bytes)).threads(threads).eof_marker(eof).build().unwrap();
    let mut data = Vec::new();
    let error = loop {
        let step: Result<bool, BgzfError> = match access {
            Access::FillAdvance(max) => match stream.fill() {
                Ok(chunk) => {
                    let n = chunk.len().min(max);
                    data.extend_from_slice(&chunk[..n]);
                    stream.advance(n);
                    Ok(n > 0)
                }
                Err(e) => Err(e),
            },
            Access::ReadExact(size) => {
                let mut buf = vec![0u8; size];
                match stream.read_exact_into(&mut buf) {
                    Ok(()) => {
                        data.extend_from_slice(&buf);
                        Ok(true)
                    }
                    // What was left is lost; read the tail of the stream again
                    // with `read_up_to`, which is what a caller would do.
                    Err(BgzfError::UnexpectedEof) => break None,
                    Err(e) => Err(e),
                }
            }
            Access::IoRead(size) => {
                let mut buf = vec![0u8; size];
                stream.read(&mut buf).map_err(bgzf_error).map(|n| {
                    data.extend_from_slice(&buf[..n]);
                    n > 0
                })
            }
            Access::ReadUntil(delim) => {
                // On error, `line` keeps what was read before it.
                let mut line = Vec::new();
                let read = stream.read_until(delim, &mut line);
                data.extend_from_slice(&line);
                read.map_err(bgzf_error).map(|n| n > 0)
            }
        };
        match step {
            Ok(true) => {}
            Ok(false) => break None,
            Err(e) => break Some(e),
        }
    };
    if let Some(e) = &error {
        let mut buf = [0u8; 8];
        assert!(
            matches!(stream.read_up_to(&mut buf), Err(BgzfError::Poisoned)),
            "not poisoned after {e:?}"
        );
    }
    Outcome { data, error }
}

/// The error's variant and message: equal for equal errors (`BgzfError` is not
/// `PartialEq`; the I/O errors it can carry are not either).
fn error_key(e: &BgzfError) -> String {
    format!("{e:?}")
}

// r[verify bgzf.stream]
// r[verify bgzf.stream.api]
// r[verify bgzf.stream.threads]
// r[verify bgzf.stream.eof_marker]
// r[verify bgzf.empty_block]
// r[verify bgzf.eof+2]
/// Whatever the block layout (seqair's and noodles' blocks, empty blocks
/// anywhere, with or without an EOF marker), thread count and way of
/// reading, the stream yields exactly what noodles' reader and flate2's
/// multi-member gzip decoder yield — and, with `EofMarker::Required`, fails
/// with `MissingEofMarker` after all of it when the marker is missing.
#[hegel::test(test_cases = pinned::cases(200))]
fn the_stream_reads_what_noodles_and_flate2_read(tc: TestCase) {
    let layout = tc.draw(layout(40, 65280));
    let (bytes, _) = layout.encode();
    let expected = layout.data();

    let mut noodles_out = Vec::new();
    noodles::bgzf::io::Reader::new(Cursor::new(&bytes)).read_to_end(&mut noodles_out).unwrap();
    assert_eq!(noodles_out, expected, "noodles");
    if !bytes.is_empty() {
        // flate2 wants at least one gzip member.
        let mut flate2_out = Vec::new();
        flate2::read::MultiGzDecoder::new(Cursor::new(&bytes))
            .read_to_end(&mut flate2_out)
            .unwrap();
        assert_eq!(flate2_out, expected, "flate2");
    }

    let threads = tc.draw(gs::integers::<usize>().max_value(4));
    let access = tc.draw(access());
    let eof = tc.draw(gs::sampled_from(vec![EofMarker::Required, EofMarker::Optional]));
    let outcome = read_with(&bytes, threads, eof, access);
    if let Access::ReadExact(_) = access {
        // The last partial chunk is dropped by `read_exact_into`.
        assert!(expected.starts_with(&outcome.data), "{access:?}");
        assert!(expected.len() - outcome.data.len() < 200_000);
    } else {
        assert_eq!(outcome.data.len(), expected.len(), "{access:?}");
        assert!(outcome.data == expected, "{access:?}: data differs");
    }
    let missing = eof == EofMarker::Required && !layout.ends_with_empty_block();
    match (&outcome.error, access) {
        (None, _) => assert!(!missing, "missing EOF marker accepted"),
        (Some(BgzfError::MissingEofMarker), _) => assert!(missing),
        (Some(e), _) => panic!("unexpected error {e:?}"),
    }
    assert_eq!(reader_all(&bytes).unwrap(), expected, "BgzfReader");
}

// r[verify bgzf.truncated_header]
// r[verify bgzf.stream.eof_marker]
// r[verify bgzf.stream.poison]
// r[verify bgzf.stream.threads]
/// Cutting a stream at every byte offset: a cut inside a block is an error
/// (never a short successful read) reported after every whole block before
/// it; a cut at a block boundary is the end with `EofMarker::Optional`, and
/// with `Required` too if the last whole block is empty, else
/// `MissingEofMarker` after all the data. Identical for every thread count,
/// and the seekable `BgzfReader` agrees (it does not ask for the marker).
#[hegel::test(test_cases = pinned::cases(40))]
fn a_stream_cut_anywhere_fails_unless_cut_at_a_block_boundary(tc: TestCase) {
    let layout = tc.draw(layout(6, 300));
    let (bytes, ends) = layout.encode();
    for cut in 0..=bytes.len() {
        let input = &bytes[..cut];
        // Pieces wholly before the cut; the EOF marker counts as the last.
        let whole = ends.iter().take_while(|&&end| end <= cut).count();
        let at_boundary = cut == 0 || ends.contains(&cut);
        let data_before = layout.data_before(whole.min(layout.pieces.len()));
        let last_whole_is_empty = match whole.checked_sub(1) {
            None => false,
            Some(i) if i == layout.pieces.len() => true,
            Some(i) => matches!(layout.pieces[i], Piece::Empty),
        };
        for eof in [EofMarker::Optional, EofMarker::Required] {
            let outcomes: Vec<Outcome> =
                [0, 1, 3].iter().map(|&threads| stream_outcome(input, threads, eof)).collect();
            let first = &outcomes[0];
            for other in &outcomes[1..] {
                assert_eq!(other.data, first.data, "cut {cut}: data differs by thread count");
                assert_eq!(
                    other.error.as_ref().map(error_key),
                    first.error.as_ref().map(error_key),
                    "cut {cut}: error differs by thread count"
                );
            }
            assert_eq!(first.data, data_before, "cut {cut} {eof:?}: data before the end");
            match (&first.error, at_boundary, eof) {
                (None, true, EofMarker::Optional) => {}
                (None, true, EofMarker::Required) => assert!(last_whole_is_empty, "cut {cut}"),
                (Some(BgzfError::MissingEofMarker), true, EofMarker::Required) => {
                    assert!(!last_whole_is_empty, "cut {cut}");
                }
                (Some(BgzfError::TruncatedBlock), false, _) => {}
                (outcome, _, _) => {
                    panic!("cut {cut} (boundary: {at_boundary}) {eof:?}: {outcome:?}")
                }
            }
        }
        match reader_all(input) {
            Ok(data) => {
                assert!(at_boundary, "BgzfReader took a cut at {cut}");
                assert_eq!(data, data_before);
            }
            Err(e) => assert!(!at_boundary, "BgzfReader failed a cut at a boundary: {e:?}"),
        }
    }
}

/// The byte ranges of a block that a reader checks, relative to its start:
/// magic, XLEN, the BC subfield, the payload, CRC and ISIZE. (MTIME, XFL and
/// OS are informational: htslib and seqair ignore them.)
fn checked_bytes(block_len: usize) -> Vec<usize> {
    [0..4, 10..18, 18..block_len].into_iter().flatten().collect()
}

// r[verify bgzf.magic]
// r[verify bgzf.bsize]
// r[verify bgzf.crc32]
// r[verify bgzf.decompression+2]
// r[verify bgzf.stream.poison]
/// Flipping any bit of a checked byte of any block — magic, XLEN, BC
/// subfield, BSIZE, CRC, ISIZE, the EOF marker's included — is an error for
/// the stream at every thread count (the same error), for `BgzfReader` and
/// for `RegionBuf`, and the stream is poisoned after it. A flipped payload
/// bit is an error exactly when flate2 rejects the stream too: the bits after
/// a DEFLATE stream's last code are padding, which no decoder can check.
#[hegel::test(test_cases = pinned::cases(300))]
fn corrupting_a_checked_byte_is_an_error(tc: TestCase) {
    let mut layout = tc.draw(layout(12, 3000));
    layout.eof_marker = true;
    let original = layout.data();
    let (mut bytes, ends) = layout.encode();
    let which = tc.draw(gs::integers::<usize>().max_value(ends.len() - 1));
    let start = which.checked_sub(1).map_or(0, |i| ends[i]);
    let len = ends[which] - start;
    let offsets = checked_bytes(len);
    let within = offsets[tc.draw(gs::integers::<usize>().max_value(offsets.len() - 1))];
    let at = start + within;
    bytes[at] ^= 1 << tc.draw(gs::integers::<u32>().max_value(7));

    let in_payload = (18..len - 8).contains(&within);
    let must_fail = !in_payload || {
        let mut out = Vec::new();
        let flate2 = flate2::read::MultiGzDecoder::new(Cursor::new(&bytes)).read_to_end(&mut out);
        flate2.is_err() || out != original
    };
    let keys: Vec<Option<String>> = [0, 2, 4]
        .iter()
        .map(|&threads| {
            let outcome = stream_outcome(&bytes, threads, EofMarker::Optional);
            assert_eq!(
                outcome.error.is_some(),
                must_fail,
                "byte {at} (threads {threads}): {:?}",
                outcome.error
            );
            if outcome.error.is_none() {
                assert!(outcome.data == original);
            }
            outcome.error.as_ref().map(error_key)
        })
        .collect();
    assert!(keys.windows(2).all(|w| w[0] == w[1]), "errors differ by thread count: {keys:?}");
    assert_eq!(reader_all(&bytes).is_err(), must_fail, "BgzfReader, byte {at}");
    // RegionBuf takes a block running past the end of its range for the end
    // of the range, so a BSIZE grown past the end of the file reads as a
    // clean end there.
    let bsize = usize::from(u16::from_le_bytes([bytes[start + 16], bytes[start + 17]]));
    let runs_past_the_end = (16..18).contains(&within) && start + bsize + 1 > bytes.len();
    if !runs_past_the_end {
        assert_eq!(region_buf_all(&bytes).is_err(), must_fail, "RegionBuf, byte {at}");
    }
}

// r[verify bgzf.empty_block]
// r[verify bgzf.eof+2]
/// The seekable `BgzfReader` reads past empty blocks mid-file exactly as
/// htslib's BGZF reader (through rust-htslib) does.
#[hegel::test(test_cases = pinned::cases(100))]
fn bgzf_reader_reads_what_htslib_reads(tc: TestCase) {
    let layout = tc.draw(layout(30, 65280));
    let (bytes, _) = layout.encode();
    let file = tempfile::NamedTempFile::new().unwrap();
    std::fs::write(file.path(), &bytes).unwrap();
    let mut htslib = Vec::new();
    rust_htslib::bgzf::Reader::from_path(file.path()).unwrap().read_to_end(&mut htslib).unwrap();
    assert!(htslib == layout.data(), "htslib disagrees with the layout");
    assert!(reader_all(&bytes).unwrap() == htslib, "BgzfReader disagrees with htslib");
}
