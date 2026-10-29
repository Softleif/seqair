//! `AlignmentStream`: one SAM, written as plain text, bgzip, gzip and BAM,
//! detected from its first bytes and read back as the same records.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "test code"
)]

#[path = "../src/pinned.rs"]
mod pinned;
#[path = "support/sam_lines.rs"]
mod sam_lines;

use hegel::TestCase;
use hegel::generators as gs;
use sam_lines::{Dialect, canonical_lines, header_text, sam_text, seqair_payloads};
use seqair::bam::{BamWriterBuilder, RawBamRecord, bgzf::EofMarker};
use seqair::reader::{AlignmentStream, AlignmentStreamError, StreamFormat};
use seqair::sam::SamStreamReader;
use std::fmt::Write as _;
use std::io::{Read, Write};

/// A reader that hands out at most `n` bytes per `read`, as a pipe may.
struct Trickle<R> {
    inner: R,
    n: usize,
}

impl<R: Read> Read for Trickle<R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let n = buf.len().min(self.n);
        self.inner.read(&mut buf[..n])
    }
}

fn bgzip(data: &[u8]) -> Vec<u8> {
    let mut w = noodles::bgzf::io::Writer::new(Vec::new());
    w.write_all(data).unwrap();
    w.finish().unwrap()
}

fn gzip(data: &[u8]) -> Vec<u8> {
    let mut w = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    w.write_all(data).unwrap();
    w.finish().unwrap()
}

/// The SAM as BAM, through seqair's own SAM reader and raw writer.
fn bam(sam: &[u8]) -> Vec<u8> {
    let mut reader = SamStreamReader::new(sam).unwrap();
    let header = reader.header().clone();
    let mut writer = BamWriterBuilder::to_writer(Vec::new(), &header).build().unwrap();
    let mut record = RawBamRecord::default();
    while reader.read_record(&mut record).unwrap() {
        writer.write_raw(&record).unwrap();
    }
    writer.finish().unwrap().0
}

fn stream(
    data: Vec<u8>,
    chunk: usize,
    threads: usize,
) -> Result<(StreamFormat, String, Vec<Vec<u8>>), AlignmentStreamError> {
    let input = Trickle { inner: std::io::Cursor::new(data), n: chunk };
    let builder = AlignmentStream::builder(input);
    let mut stream = match u8::try_from(threads).ok().and_then(std::num::NonZeroU8::new) {
        Some(n) => builder.threads(n),
        None => builder,
    }
    .build()?;
    let header = stream.header().header_text().to_owned();
    let mut records = Vec::new();
    let mut record = RawBamRecord::default();
    while stream.read_record(&mut record)? {
        records.push(record.as_payload().to_vec());
    }
    Ok((stream.format(), header, records))
}

// r[verify unified.stream.detect+2]
// r[verify unified.stream.records]
/// Every wrapping of the same SAM is detected as what it is, however the
/// bytes arrive, and yields the same header and records.
#[hegel::test(test_cases = crate::pinned::cases(100))]
fn every_format_is_detected_and_reads_the_same_records(tc: TestCase) {
    let lines = tc.draw(canonical_lines(Dialect { noodles: false }, 5));
    let sam = sam_text(&lines).into_bytes();
    let want = seqair_payloads(&sam).unwrap();
    let chunk = tc.draw(gs::integers::<usize>().min_value(1).max_value(64));
    let threads = tc.draw(gs::integers::<usize>().max_value(3));
    for (format, data) in [
        (StreamFormat::Sam, sam.clone()),
        (StreamFormat::BgzfSam, bgzip(&sam)),
        (StreamFormat::GzipSam, gzip(&sam)),
        (StreamFormat::Bam, bam(&sam)),
    ] {
        let (found, header, records) = stream(data, chunk, threads).unwrap();
        assert_eq!(found, format);
        assert_eq!(header, header_text(), "{format:?}");
        assert_eq!(records, want, "{format:?}");
    }
}

// r[verify unified.stream.detect+2]
#[test]
fn empty_input_is_an_empty_sam() {
    let (format, header, records) = stream(Vec::new(), 1, 0).unwrap();
    assert_eq!(format, StreamFormat::Sam);
    assert!(header.is_empty());
    assert!(records.is_empty());
}

// r[verify unified.stream.detect+2]
#[test]
fn unsupported_formats_are_named() {
    let sam = format!("{}r\t4\t*\t0\t0\t*\t*\t0\t0\tA\tI\n", header_text()).into_bytes();
    let bam = bam(&sam);
    let mut raw_bam = Vec::new();
    seqair::bam::bgzf::BgzfStream::new(&bam[..]).read_to_end(&mut raw_bam).unwrap();

    let cram = b"CRAM\x03\x00rest-of-a-cram-file".to_vec();
    assert!(matches!(stream(cram, 3, 0), Err(AlignmentStreamError::CramNeedsReference)));
    assert!(matches!(stream(gzip(&raw_bam), 1, 0), Err(AlignmentStreamError::GzipBamNotSupported)));
    assert!(matches!(
        stream(raw_bam, 2, 0),
        Err(AlignmentStreamError::UncompressedBamNotSupported)
    ));
}

// r[verify unified.stream.detect+2]
/// A BAM cut at a block boundary ends in `MissingEofMarker` by default, and
/// reads to the end with `EofMarker::Optional`.
#[test]
fn eof_marker_option_reaches_the_bgzf_layer() {
    let mut body = header_text();
    for i in 0..20 {
        writeln!(body, "r{i}\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII").unwrap();
    }
    let mut data = bam(body.as_bytes());
    data.truncate(data.len() - 28); // the EOF marker block
    let read = |marker| -> Result<usize, AlignmentStreamError> {
        let mut s = AlignmentStream::builder(std::io::Cursor::new(data.clone()))
            .eof_marker(marker)
            .build()?;
        let mut n = 0;
        let mut r = RawBamRecord::default();
        while s.read_record(&mut r)? {
            n += 1;
        }
        Ok(n)
    };
    assert_eq!(read(EofMarker::Optional).unwrap(), 20);
    let err = read(EofMarker::Required).unwrap_err();
    assert!(
        matches!(&err, AlignmentStreamError::Bam { source } if std::error::Error::source(source)
            .and_then(|s| s.downcast_ref::<seqair::bam::BgzfError>())
            .is_some_and(|e| matches!(e, seqair::bam::BgzfError::MissingEofMarker))),
        "{err:?}"
    );
}
