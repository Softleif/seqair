//! `SamStreamReader`: htslib's own SAM corpus, compressed inputs, line
//! endings, line numbers and failure modes.
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
use sam_lines::{
    Dialect, canonical_lines, header_text, sam_text, samtools_payloads, seqair_payloads,
};
use seqair::bam::{BamHeader, BgzfStream, RawBamRecord};
use seqair::sam::{SamRecordEncoder, SamRecordError, SamStreamError, SamStreamReader};
use std::io::{BufRead, BufReader, Read, Write};
use std::path::Path;

fn corpus_dir() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/htslib/sam"))
}

// r[verify sam.stream]
// r[verify sam.encode]
/// Every file of htslib's own SAM test corpus: where `samtools view -b`
/// converts it, seqair yields the same records byte for byte; where samtools
/// fails, seqair fails too.
#[test]
fn htslibs_sam_corpus_encodes_to_htslibs_bytes() {
    let mut files: Vec<_> = std::fs::read_dir(corpus_dir())
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.extension().is_some_and(|e| e == "sam"))
        .collect();
    files.sort();
    assert!(files.len() >= 40, "corpus missing: {}", files.len());
    let (mut agreed, mut rejected) = (0, 0);
    for path in &files {
        let sam = std::fs::read(path).unwrap();
        let name = path.file_name().unwrap().to_string_lossy();
        match (samtools_payloads(&sam), seqair_payloads(&sam)) {
            (Ok(theirs), Ok(ours)) => {
                assert_eq!(ours.len(), theirs.len(), "{name}: record count");
                for (i, (a, b)) in ours.iter().zip(&theirs).enumerate() {
                    assert_eq!(a, b, "{name}: record {i}");
                }
                agreed += 1;
            }
            (Err(_), Err(_)) => rejected += 1,
            (Ok(_), Err(e)) => panic!("{name}: htslib reads it, seqair fails: {e:?}"),
            (Err(e), Ok(_)) => panic!("{name}: htslib fails, seqair reads it: {e}"),
        }
    }
    assert!(agreed > 30, "only {agreed} files agreed ({rejected} rejected by both)");
}

/// The records of `text`, each line encoded on its own — what the stream
/// must reproduce whatever wraps the text.
fn line_by_line(text: &str) -> Vec<Vec<u8>> {
    let header = BamHeader::from_sam_text(&header_text()).unwrap();
    let mut encoder = SamRecordEncoder::new();
    text.lines()
        .filter(|l| !l.starts_with('@'))
        .map(|l| {
            let mut r = RawBamRecord::default();
            encoder.encode(l.as_bytes(), &header, &mut r).unwrap();
            r.into_payload()
        })
        .collect()
}

fn read_all<R: BufRead>(reader: R) -> Result<(String, Vec<Vec<u8>>), SamStreamError> {
    let mut s = SamStreamReader::new(reader)?;
    let records =
        s.records().map(|r| r.map(RawBamRecord::into_payload)).collect::<Result<_, _>>()?;
    Ok((s.header().header_text().to_string(), records))
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

/// A reader that hands out at most `n` bytes per `read`, so `BufReader`
/// over it splits lines anywhere.
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

// r[verify sam.stream]
// r[verify sam.stream.header]
// r[verify sam.stream.lines]
/// Plain, bgzipped and gzipped text, with `\n` or `\r\n` endings, with or
/// without a final newline, fed a few bytes per read: always the header
/// text with `\n` endings and the records each line encodes to.
#[hegel::test(test_cases = crate::pinned::cases(200))]
fn every_wrapping_of_the_text_gives_the_same_records(tc: TestCase) {
    let lines = tc.draw(canonical_lines(Dialect { noodles: false }, 5));
    let text = sam_text(&lines);
    let want = line_by_line(&text);
    let crlf = tc.draw(gs::booleans());
    let mut input = if crlf { text.replace('\n', "\r\n") } else { text.clone() }.into_bytes();
    if tc.draw(gs::booleans()) {
        input.truncate(input.len() - if crlf { 2 } else { 1 });
    }
    let chunk = tc.draw(gs::integers::<usize>().min_value(1).max_value(64));
    let capacity = tc.draw(gs::integers::<usize>().min_value(1).max_value(64));
    let wrapped: Box<dyn BufRead> = match tc.draw(gs::integers::<u8>().max_value(2)) {
        0 => Box::new(BufReader::with_capacity(capacity, Trickle { inner: &input[..], n: chunk })),
        1 => Box::new(BgzfStream::new(Trickle {
            inner: std::io::Cursor::new(bgzip(&input)),
            n: chunk,
        })),
        _ => Box::new(BufReader::with_capacity(
            capacity,
            flate2::read::MultiGzDecoder::new(Trickle {
                inner: std::io::Cursor::new(gzip(&input)),
                n: chunk,
            }),
        )),
    };
    let (header, records) = read_all(wrapped).unwrap();
    assert_eq!(header, header_text());
    assert_eq!(records, want);
}

// r[verify sam.stream.errors]
/// A truncated gzip or BGZF stream is an error, never a short clean end,
/// and the stream stays failed.
#[test]
fn truncated_compressed_input_errors() {
    let mut body = header_text().into_bytes();
    for i in 0..20_000 {
        body.extend_from_slice(format!("r{i}\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n").as_bytes());
    }
    for (name, compressed) in [("gzip", gzip(&body)), ("bgzf", bgzip(&body))] {
        for cut in [compressed.len() - 1, compressed.len() - 8, compressed.len() / 2] {
            let data = &compressed[..cut];
            let result = if name == "gzip" {
                read_all(BufReader::new(flate2::read::MultiGzDecoder::new(data)))
            } else {
                read_all(BgzfStream::new(data))
            };
            assert!(
                matches!(result, Err(SamStreamError::Io { .. })),
                "{name} cut at {cut}: {:?}",
                result.map(|(_, r)| r.len())
            );
        }
        let data = &compressed[..compressed.len() / 2];
        let mut reader = if name == "gzip" {
            Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(data))) as Box<dyn BufRead>
        } else {
            Box::new(BgzfStream::new(data))
        };
        let mut s = SamStreamReader::new(&mut reader).unwrap();
        let mut r = RawBamRecord::default();
        let err = loop {
            match s.read_record(&mut r) {
                Ok(true) => {}
                Ok(false) => panic!("{name}: truncated input ended cleanly"),
                Err(e) => break e,
            }
        };
        assert!(matches!(err, SamStreamError::Io { .. }), "{err:?}");
        assert!(matches!(s.read_record(&mut r), Err(SamStreamError::Poisoned)));
    }
}

// r[verify sam.stream.lines]
// r[verify sam.stream.errors]
/// A blank line fails as it does in htslib, with its 1-based line number
/// (header lines counted); the next line still reads.
#[test]
fn a_blank_line_is_an_error_and_the_stream_goes_on() {
    let text = format!(
        "{}a\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n\nb\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n",
        header_text()
    );
    assert!(samtools_payloads(text.as_bytes()).is_err());
    let n_header = header_text().lines().count() as u64;
    let mut s = SamStreamReader::new(text.as_bytes()).unwrap();
    let mut r = RawBamRecord::default();
    assert!(s.read_record(&mut r).unwrap());
    assert_eq!(r.qname(), b"a");
    let err = s.read_record(&mut r).unwrap_err();
    assert!(
        matches!(err, SamStreamError::InvalidRecord { line, source: SamRecordError::EmptyLine } if line == n_header + 2),
        "{err:?}"
    );
    assert!(s.read_record(&mut r).unwrap());
    assert_eq!(r.qname(), b"b");
    assert!(!s.read_record(&mut r).unwrap());
    assert!(!s.read_record(&mut r).unwrap(), "the end stays the end");
}

// r[verify sam.stream.header]
/// No header, or one without `@SQ`: unplaced records read, a record naming
/// a reference fails with its line number.
#[test]
fn a_sam_without_sq_lines_streams_its_unplaced_reads() {
    for header in ["", "@HD\tVN:1.6\n@CO\tno references\n"] {
        let text = format!(
            "{header}u\t4\t*\t0\t0\t*\t*\t0\t0\tACGU\t*\nm\t0\tchr1\t5\t0\t4M\t*\t0\t0\tACGT\t*\n"
        );
        let mut s = SamStreamReader::new(text.as_bytes()).unwrap();
        assert_eq!(s.header().header_text(), header);
        assert_eq!(s.header().target_count(), 0);
        let mut r = RawBamRecord::default();
        assert!(s.read_record(&mut r).unwrap());
        let mut seq = Vec::new();
        r.decode_seq_into(&mut seq);
        assert_eq!(seq, b"ACGT");
        let line = header.lines().count() as u64 + 2;
        assert!(matches!(
            s.read_record(&mut r),
            Err(SamStreamError::InvalidRecord { line: l, source: SamRecordError::UnknownReference { .. } }) if l == line
        ));
    }
}

// r[verify sam.stream.header]
#[test]
fn a_header_line_that_is_not_utf8_names_its_line() {
    let text = b"@HD\tVN:1.6\n@CO\t\xff\n";
    let err = SamStreamReader::new(&text[..]).unwrap_err();
    assert!(matches!(err, SamStreamError::HeaderNotUtf8 { line: 2, .. }), "{err:?}");
}

// r[verify sam.stream.lines]
// r[verify sam.stream.errors]
/// Input without newlines is cut off at 16 MiB instead of filling memory,
/// and ends the stream.
#[test]
fn an_endless_line_is_cut_off() {
    let endless = BufReader::new(std::io::repeat(b'A').take(17 * 1024 * 1024));
    let mut s = SamStreamReader::new(endless).unwrap();
    let mut r = RawBamRecord::default();
    assert!(matches!(s.read_record(&mut r), Err(SamStreamError::LineTooLong { line: 1 })));
    assert!(matches!(s.read_record(&mut r), Err(SamStreamError::Poisoned)));
}
