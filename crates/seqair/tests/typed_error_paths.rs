//! Which error comes back — not merely that one does.
//!
//! The fuzz targets prove that malformed input never panics. They say nothing
//! about *which* typed variant a given failure produces, so an error enum can
//! quietly collapse into a single catch-all without a test going red. These
//! tests pin the variant, the typed fields it carries, and — where the enums
//! nest — the wrapping chain, so a rewrapping that loses the cause is caught.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::unwrap_in_result,
    reason = "test code"
)]
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    reason = "test code with known small values"
)]

use hegel::prelude::*;
use seqair::bam::aux::AuxValue;
use seqair::bam::aux_data::{AuxData, AuxDataError};
use seqair::bam::cigar::{CigarOp, CigarOpType};
use seqair::bam::header::BamHeader;
use seqair::bam::owned_record::{OwnedBamRecord, OwnedRecordError};
use seqair::bam::writer::{BamWriteError, BamWriter, BamWriterBuilder};
use seqair::bam::{
    BamError, BamHeaderError, BgzfError, IndexedBamReader, Pos0, RecordIdx, RecordStore,
};
use seqair::io::IndexError;
use seqair_types::bam_flags::BamFlags;
use seqair_types::{Base, BaseQuality};
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::sync::{Arc, Mutex, OnceLock};

// ── A BAM that spans several BGZF blocks ────────────────────────────────

/// A valid, coordinate-sorted BAM with its BAI, plus the file offsets of every
/// BGZF block in it.
///
/// The records carry pseudo-random sequence and quality so the blocks do not
/// compress down to one — the reader-propagation tests need a header block and
/// several record blocks that can be corrupted independently.
struct Corpus {
    bam: Vec<u8>,
    bai: Vec<u8>,
    blocks: Vec<Block>,
}

/// One BGZF block: `start` is its file offset, `len` its total on-disk size
/// (`BSIZE + 1`). The 18-byte header is at `start`, the deflate payload at
/// `start + 18`, and the 8-byte CRC32/ISIZE footer at `start + len - 8`.
#[derive(Debug, Clone, Copy)]
struct Block {
    start: usize,
    len: usize,
}

impl Block {
    /// The first deflate byte: it carries the block's type and final-block bit,
    /// so flipping it reliably breaks the stream rather than one literal.
    fn payload_head(self) -> usize {
        self.start + 18
    }
    fn crc_byte(self) -> usize {
        self.start + self.len - 8
    }
    fn isize_byte(self) -> usize {
        self.start + self.len - 1
    }
}

/// Walk the BGZF block chain. Every block seqair writes carries the standard
/// 6-byte `BC` extra subfield, so BSIZE is the u16 at `start + 16`.
fn bgzf_blocks(bam: &[u8]) -> Vec<Block> {
    let mut blocks = Vec::new();
    let mut at = 0usize;
    while at + 18 <= bam.len() {
        assert_eq!(&bam[at..at + 2], &[0x1f, 0x8b], "BGZF magic at block start");
        assert_eq!(&bam[at + 12..at + 14], b"BC", "standard BSIZE subfield");
        let bsize = u16::from_le_bytes([bam[at + 16], bam[at + 17]]) as usize;
        let len = bsize + 1;
        assert!(len > 26, "block must hold header + footer");
        blocks.push(Block { start: at, len });
        at += len;
    }
    assert_eq!(at, bam.len(), "block chain must tile the file exactly");
    blocks
}

fn corpus() -> &'static Corpus {
    static CORPUS: OnceLock<Corpus> = OnceLock::new();
    CORPUS.get_or_init(|| {
        let header =
            BamHeader::from_sam_text("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:2000000\n")
                .unwrap();
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("corpus.bam");
        let mut writer =
            BamWriterBuilder::to_path(&path, &header).write_index(true).build().unwrap();

        // A deterministic LCG: the bytes must not compress away, but the corpus
        // must be identical on every run so the offsets in the sweep below mean
        // the same thing every time.
        let mut state: u32 = 0x1234_5678;
        let mut next = move || {
            state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
            state
        };
        for i in 0..1500u32 {
            let seq: Vec<Base> = (0..100)
                .map(|_| match (next() >> 28) & 3 {
                    0 => Base::A,
                    1 => Base::C,
                    2 => Base::G,
                    _ => Base::T,
                })
                .collect();
            let qual: Vec<BaseQuality> =
                (0..100).map(|_| BaseQuality::from_byte((next() % 40) as u8)).collect();
            let record =
                OwnedBamRecord::builder(0, Some(Pos0::new(i * 100 + 1).unwrap()), read_name(i))
                    .flags(BamFlags::empty())
                    .mapq(60)
                    .cigar(vec![CigarOp::new(CigarOpType::Match, 100)])
                    .seq(seq)
                    .qual(qual)
                    .build()
                    .unwrap();
            writer.write(&record).unwrap();
        }
        let (inner, index) = writer.finish().unwrap();
        drop(inner);
        let mut bai = Vec::new();
        index.unwrap().write_bai(&mut bai, header.target_count()).unwrap();
        let bam = std::fs::read(&path).unwrap();
        let blocks = bgzf_blocks(&bam);
        assert!(blocks.len() >= 4, "corpus must span several BGZF blocks, got {}", blocks.len());
        Corpus { bam, bai, blocks }
    })
}

fn read_name(i: u32) -> Vec<u8> {
    format!("read{i}").into_bytes()
}

/// Open `bam` (with the corpus's intact BAI beside it) and query the whole
/// contig, returning whichever [`BamError`] surfaces first.
fn open_and_fetch(bam: &[u8]) -> Result<usize, BamError> {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("t.bam");
    std::fs::write(&path, bam).unwrap();
    std::fs::write(dir.path().join("t.bam.bai"), &corpus().bai).unwrap();
    let mut reader = IndexedBamReader::open(&path)?;
    let mut store = RecordStore::default();
    reader.fetch_into(0, Pos0::new(0).unwrap(), Pos0::new(1_999_999).unwrap(), &mut store)
}

/// The `BgzfError` at the bottom of a `BamError`, whichever way it was wrapped,
/// together with how it got there. `None` for a `BamError` that is not a
/// rewrapped BGZF failure at all.
fn bgzf_cause(err: &BamError) -> Option<(&'static str, &BgzfError)> {
    match err {
        BamError::Bgzf { source } => Some(("BamError::Bgzf", source)),
        BamError::Header { source: BamHeaderError::Bgzf { source } } => {
            Some(("BamError::Header/BamHeaderError::Bgzf", source))
        }
        _ => None,
    }
}

/// A stable name for a `BgzfError` variant, for tallying which ones a sweep
/// actually reached.
fn bgzf_variant(err: &BgzfError) -> &'static str {
    match err {
        BgzfError::InvalidMagic => "InvalidMagic",
        BgzfError::MissingBsize => "MissingBsize",
        BgzfError::BlockSizeTooSmall { .. } => "BlockSizeTooSmall",
        BgzfError::TruncatedBlock => "TruncatedBlock",
        BgzfError::DecompressionFailed { .. } => "DecompressionFailed",
        BgzfError::ChecksumMismatch { .. } => "ChecksumMismatch",
        BgzfError::UnexpectedEof => "UnexpectedEof",
        BgzfError::UncompressedSizeTooLarge { .. } => "UncompressedSizeTooLarge",
        other => panic!("unexpected BgzfError from a corrupted BAM: {other:?}"),
    }
}

/// The corruptions this file applies, each aimed at one decode step in
/// `read_block`, and the `BgzfError` each one is expected to raise.
#[derive(Debug, Clone, Copy)]
enum Corruption {
    /// Flip the gzip magic — rejected before anything else is read.
    Magic,
    /// Rename the `BC` extra subfield so BSIZE cannot be found.
    BsizeSubfieldId,
    /// Shrink BSIZE below the 18-byte header.
    BsizeTooSmall,
    /// Flip a byte in the deflate payload.
    Payload,
    /// Flip a byte of the stored CRC32.
    Crc,
    /// Set ISIZE past the 64 KiB uncompressed block ceiling.
    Isize,
}

impl Corruption {
    const ALL: [Self; 6] = [
        Self::Magic,
        Self::BsizeSubfieldId,
        Self::BsizeTooSmall,
        Self::Payload,
        Self::Crc,
        Self::Isize,
    ];

    /// The `BgzfError` variant name this corruption must produce.
    ///
    /// A payload flip is the one that can land either way: libdeflate usually
    /// rejects the stream outright (`DecompressionFailed`), but a flip that
    /// still decodes leaves the CRC to catch it.
    fn expected(self) -> &'static [&'static str] {
        match self {
            Self::Magic => &["InvalidMagic"],
            Self::BsizeSubfieldId => &["MissingBsize"],
            Self::BsizeTooSmall => &["BlockSizeTooSmall"],
            Self::Payload => &["DecompressionFailed", "ChecksumMismatch"],
            Self::Crc => &["ChecksumMismatch"],
            Self::Isize => &["UncompressedSizeTooLarge"],
        }
    }

    fn apply(self, bam: &mut [u8], block: Block) {
        match self {
            Self::Magic => bam[block.start] ^= 0xff,
            Self::BsizeSubfieldId => bam[block.start + 12] = b'X',
            Self::BsizeTooSmall => {
                bam[block.start + 16] = 5;
                bam[block.start + 17] = 0;
            }
            Self::Payload => bam[block.payload_head()] ^= 0xff,
            Self::Crc => bam[block.crc_byte()] ^= 0xff,
            Self::Isize => bam[block.isize_byte()] = 0xff,
        }
    }
}

// r[verify bam.reader.propagate_errors]
/// Every distinct BGZF decode failure must reach the caller as its own
/// variant, wrapped but not flattened.
///
/// The sweep corrupts the header block (surfacing through
/// `IndexedBamReader::open`) and a record block (surfacing through
/// `fetch_into`), and tallies which `BgzfError` variants came back. It then
/// fails if any enumerated variant was never reached — a test that accepts
/// "any of these six" and only ever sees one proves nothing.
#[test]
fn every_bgzf_decode_failure_surfaces_as_its_own_variant() {
    let corpus = corpus();
    let header_block = corpus.blocks[0];
    // The last block is the 28-byte BGZF EOF marker; the one before it is a
    // record block far enough in that `open` has already succeeded.
    let record_block = corpus.blocks[corpus.blocks.len() - 2];

    let mut seen: BTreeMap<&'static str, usize> = BTreeMap::new();

    for corruption in Corruption::ALL {
        for (label, block, wrapper) in [
            ("header block", header_block, "BamError::Header/BamHeaderError::Bgzf"),
            ("record block", record_block, "BamError::Bgzf"),
        ] {
            let mut bam = corpus.bam.clone();
            corruption.apply(&mut bam, block);
            let err = open_and_fetch(&bam)
                .expect_err(&format!("{corruption:?} in the {label} must be rejected"));
            let (chain, cause) = bgzf_cause(&err).unwrap_or_else(|| {
                panic!("{corruption:?} in the {label} lost its BGZF cause: {err:?}")
            });
            assert_eq!(
                chain, wrapper,
                "{corruption:?} in the {label} took the wrong wrapping chain: {err:?}"
            );
            let variant = bgzf_variant(cause);
            assert!(
                corruption.expected().contains(&variant),
                "{corruption:?} in the {label} gave {variant}, expected one of {:?} ({err:?})",
                corruption.expected(),
            );
            *seen.entry(variant).or_default() += 1;
        }
    }

    // Truncation: cutting inside the first block leaves `open` with either no
    // bytes at all or a block header it cannot complete.
    for cut in [0, 1, 17, header_block.len / 2, header_block.len - 1] {
        let err = open_and_fetch(&corpus.bam[..cut])
            .expect_err("a BAM truncated inside its first BGZF block must be rejected");
        let (chain, cause) = bgzf_cause(&err)
            .unwrap_or_else(|| panic!("truncation at {cut} lost its BGZF cause: {err:?}"));
        assert_eq!(chain, "BamError::Header/BamHeaderError::Bgzf", "truncation at {cut}: {err:?}");
        let variant = bgzf_variant(cause);
        assert!(
            matches!(variant, "UnexpectedEof" | "TruncatedBlock"),
            "truncation at {cut} gave {variant} ({err:?})",
        );
        *seen.entry(variant).or_default() += 1;
    }

    let expected: &[&str] = &[
        "BlockSizeTooSmall",
        "ChecksumMismatch",
        "DecompressionFailed",
        "InvalidMagic",
        "MissingBsize",
        "TruncatedBlock",
        "UncompressedSizeTooLarge",
        "UnexpectedEof",
    ];
    let reached: Vec<&str> = seen.keys().copied().collect();
    assert_eq!(reached, expected, "the sweep must reach every enumerated variant exactly");
}

// r[verify bam.reader.propagate_errors]
/// A flipped CRC32 must come back as `ChecksumMismatch` carrying *both* checksums,
/// and `expected` must be the block's real stored CRC — a variant built with the
/// wrong data reads the same as a correct one until the fields are checked.
#[test]
fn a_flipped_checksum_reports_the_stored_crc_and_the_computed_one() {
    let corpus = corpus();
    let block = corpus.blocks[corpus.blocks.len() - 2];
    let stored =
        u32::from_le_bytes(corpus.bam[block.crc_byte()..block.crc_byte() + 4].try_into().unwrap());

    let mut bam = corpus.bam.clone();
    bam[block.crc_byte()] ^= 0xff;
    let corrupted =
        u32::from_le_bytes(bam[block.crc_byte()..block.crc_byte() + 4].try_into().unwrap());

    let err = open_and_fetch(&bam).expect_err("a flipped CRC32 must be rejected");
    let Some(("BamError::Bgzf", BgzfError::ChecksumMismatch { expected, found })) =
        bgzf_cause(&err)
    else {
        panic!("expected BamError::Bgzf(ChecksumMismatch), got {err:?}");
    };
    // The reader reports the checksum it read from the file (the corrupted one)
    // against the checksum it computed over the payload (the original).
    assert_eq!(*expected, corrupted, "`expected` must be the CRC32 as stored in the file");
    assert_eq!(*found, stored, "`found` must be the CRC32 computed over the decoded payload");
    assert_ne!(expected, found, "the two checksums must differ, or nothing was detected");
}

// r[verify bam.reader.propagate_errors]
/// The same six corruptions, applied at a drawn block and read back through
/// whichever entry point owns that block: the variant set is closed, and the
/// BGZF cause is never dropped on the way out.
#[hegel::test(test_cases = 48)]
fn a_corrupted_block_always_arrives_with_its_bgzf_cause_intact(tc: TestCase) {
    let corpus = corpus();
    // Block 0 holds the header; the last block is the EOF marker, whose 28
    // bytes carry no payload to corrupt.
    let idx = tc.draw(gs::integers::<usize>().max_value(corpus.blocks.len() - 2));
    let block = corpus.blocks[idx];
    let which = tc.draw(gs::integers::<usize>().max_value(Corruption::ALL.len() - 1));
    let corruption = Corruption::ALL[which];
    tc.event(if idx == 0 { "header block" } else { "record block" });
    tc.event(format!("{corruption:?}"));

    let mut bam = corpus.bam.clone();
    corruption.apply(&mut bam, block);

    let err = open_and_fetch(&bam)
        .expect_err("a BAM with a corrupted BGZF block must not read back clean");
    let Some((chain, cause)) = bgzf_cause(&err) else {
        panic!("{corruption:?} at block {idx} lost its BGZF cause: {err:?}");
    };
    // Which wrapper applies is decided by where the block sits, not by the
    // corruption: block 0 is consumed by the header parse, the rest by the query.
    let expected_chain =
        if idx == 0 { "BamError::Header/BamHeaderError::Bgzf" } else { "BamError::Bgzf" };
    assert_eq!(chain, expected_chain, "{corruption:?} at block {idx}: {err:?}");
    let variant = bgzf_variant(cause);
    assert!(
        corruption.expected().contains(&variant),
        "{corruption:?} at block {idx} gave {variant}, expected one of {:?} ({err:?})",
        corruption.expected(),
    );
}

// ── BamWriter: which error, and whether the stream was touched ──────────

/// A `Write` sink that accepts `limit` bytes and then fails, so the writer's
/// poisoning can be driven by a real I/O failure rather than a validation one.
#[derive(Debug, Clone)]
struct FailAfter {
    state: Arc<Mutex<SinkState>>,
    limit: usize,
}

#[derive(Debug, Default)]
struct SinkState {
    written: usize,
    failed: bool,
}

impl FailAfter {
    fn new(limit: usize) -> Self {
        Self { state: Arc::new(Mutex::new(SinkState::default())), limit }
    }
    fn written(&self) -> usize {
        self.state.lock().unwrap().written
    }
}

impl Write for FailAfter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let mut state = self.state.lock().unwrap();
        if state.written + buf.len() > self.limit {
            state.failed = true;
            return Err(io::Error::from(io::ErrorKind::BrokenPipe));
        }
        state.written += buf.len();
        Ok(buf.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

fn writer_header() -> BamHeader {
    BamHeader::from_sam_text("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n").unwrap()
}

/// A minimal mapped record. `len` bases of `A` with matching qualities.
fn mapped_record(ref_id: i32, pos: u32, name: &str, len: usize) -> OwnedBamRecord {
    OwnedBamRecord::builder(ref_id, Some(Pos0::new(pos).unwrap()), name.as_bytes().to_vec())
        .flags(BamFlags::empty())
        .mapq(60)
        .cigar(vec![CigarOp::new(CigarOpType::Match, len as u32)])
        .seq(vec![Base::A; len])
        .qual(vec![BaseQuality::from_byte(30); len])
        .build()
        .unwrap()
}

/// Bytes of a BAM that got a header and nothing else — the baseline a write
/// that was rejected before reaching BGZF must leave behind.
fn header_only_bytes(dir: &std::path::Path) -> Vec<u8> {
    let header = writer_header();
    let path = dir.join("control.bam");
    let writer = BamWriterBuilder::to_path(&path, &header).write_index(true).build().unwrap();
    writer.finish().unwrap();
    std::fs::read(&path).unwrap()
}

/// Run `attempt` against a fresh indexed path-target writer and return the
/// error it produced together with the file that was left behind.
fn rejected_write(
    dir: &std::path::Path,
    attempt: impl FnOnce(&mut BamWriter<BufWriter<File>>) -> BamWriteError,
) -> (BamWriteError, Vec<u8>) {
    let header = writer_header();
    let path = dir.join("subject.bam");
    let mut writer = BamWriterBuilder::to_path(&path, &header).write_index(true).build().unwrap();
    let err = attempt(&mut writer);
    writer.finish().unwrap();
    (err, std::fs::read(&path).unwrap())
}

// r[verify bam_writer.error_type]
// r[verify bam_writer.index_record_dispatch]
// r[verify bam_writer.validate_before_write]
/// A mapped record with `ref_id == -1` is structurally unindexable. The
/// validation runs *before* the BGZF write, so the rejection must leave a file
/// byte-identical to one that never saw the record — the part a test that only
/// checks the error variant would miss.
#[test]
fn a_mapped_record_without_a_reference_is_rejected_before_the_stream_is_touched() {
    let dir = tempfile::tempdir().unwrap();
    let control = header_only_bytes(dir.path());

    let (err, actual) = rejected_write(dir.path(), |writer| {
        let record = mapped_record(-1, 100, "orphan", 10);
        writer.write(&record).expect_err("a mapped record with ref_id == -1 must be refused")
    });

    assert!(
        matches!(err, BamWriteError::MappedWithoutReference),
        "expected MappedWithoutReference, got {err:?}",
    );
    assert_eq!(actual, control, "the rejected record must not have reached the BGZF stream");
}

// r[verify bam_writer.error_type]
// r[verify bam_writer.record_size_limit]
// r[verify bam_writer.validate_before_write]
/// The 2 MiB record limit is checked against the serialized bytes, before the
/// stream. `size` must be the real serialized length, not a guess.
#[test]
fn an_oversized_record_reports_its_serialized_size_and_leaves_the_stream_alone() {
    let dir = tempfile::tempdir().unwrap();
    let control = header_only_bytes(dir.path());

    let record = mapped_record(0, 100, "huge", 3_000_000);
    let mut serialized = Vec::new();
    record.to_bam_bytes(&mut serialized).unwrap();
    let serialized_len = serialized.len();
    assert!(serialized_len > 2 * 1024 * 1024, "the fixture must exceed the limit");

    let (err, actual) = rejected_write(dir.path(), |writer| {
        writer.write(&record).expect_err("a record past the 2 MiB limit must be refused")
    });

    let BamWriteError::RecordTooLarge { size } = err else {
        panic!("expected RecordTooLarge, got {err:?}");
    };
    assert_eq!(size, serialized_len, "`size` must be the serialized record length");
    assert_eq!(actual, control, "the oversized record must not have reached the BGZF stream");
}

// r[verify bam_writer.error_type]
// r[verify bam_writer.validate_before_write]
// r[verify bam.owned_record.seq_qual_length_at_serialization]
/// Serialization failures arrive wrapped as `BamWriteError::Record`, carrying
/// the `OwnedRecordError` that caused them — and the record never reaches BGZF,
/// because `to_bam_bytes` validates before the writer emits anything.
#[test]
fn a_record_whose_qual_no_longer_matches_its_seq_is_refused_with_both_lengths() {
    let dir = tempfile::tempdir().unwrap();
    let control = header_only_bytes(dir.path());

    // An unmapped record has no CIGAR to contradict, so `set_seq` accepts a
    // shorter sequence and leaves the ten quality scores behind it.
    let mut record = OwnedBamRecord::builder(0, Some(Pos0::new(100).unwrap()), b"shrunk".to_vec())
        .flags(BamFlags::from(4))
        .seq(vec![Base::A; 10])
        .qual(vec![BaseQuality::from_byte(30); 10])
        .build()
        .unwrap();
    record.set_seq(vec![Base::C; 4]).expect("no CIGAR to contradict the new length");

    let (err, actual) = rejected_write(dir.path(), |writer| {
        writer.write(&record).expect_err("seq/qual must agree at serialization time")
    });

    let BamWriteError::Record {
        source: OwnedRecordError::SeqQualLengthMismatch { seq_len, qual_len },
    } = err
    else {
        panic!("expected Record(SeqQualLengthMismatch), got {err:?}");
    };
    assert_eq!((seq_len, qual_len), (4, 10), "both lengths must be reported as they are");
    assert_eq!(actual, control, "a record that failed to serialize must not reach the stream");
}

// r[verify bam_writer.error_type]
// r[verify bam_writer.index_sort_order]
/// Out-of-order input fails in the index builder, which runs *after* the BGZF
/// write — so unlike the validations above, this record is already in the
/// stream when the error comes back. Pinning both halves of that asymmetry is
/// what keeps `r[bam_writer.index_record_dispatch]`'s ordering honest.
#[test]
fn an_out_of_order_record_fails_in_the_index_after_it_was_written() {
    let dir = tempfile::tempdir().unwrap();
    let header = writer_header();
    let path = dir.path().join("unsorted.bam");
    let mut writer = BamWriterBuilder::to_path(&path, &header).write_index(true).build().unwrap();
    writer.write(&mapped_record(0, 5_000, "first", 10)).unwrap();
    let err = writer
        .write(&mapped_record(0, 1_000, "backwards", 10))
        .expect_err("a record before the previous one must be refused");
    writer.finish().unwrap();
    let with_both = std::fs::read(&path).unwrap();

    let BamWriteError::Index { source: IndexError::UnsortedInput { tid, pos } } = err else {
        panic!("expected Index(UnsortedInput), got {err:?}");
    };
    assert_eq!((tid, pos), (0, 1_000), "the offending record's tid and position");

    // The same writer without the second record produces a strictly shorter
    // file: the rejected record really did go into the stream first.
    let control_path = dir.path().join("sorted.bam");
    let mut control =
        BamWriterBuilder::to_path(&control_path, &header).write_index(true).build().unwrap();
    control.write(&mapped_record(0, 5_000, "first", 10)).unwrap();
    control.finish().unwrap();
    let with_one = std::fs::read(&control_path).unwrap();
    assert_ne!(with_both, with_one, "the unsorted record was written before the index refused it");
}

// r[verify bam_writer.error_type]
// r[verify record_store.record_idx.resolution]
/// `write_store_record` with an index the store does not hold reports the index
/// it was handed, not a generic failure.
#[test]
fn writing_a_record_the_store_does_not_hold_names_the_index() {
    let dir = tempfile::tempdir().unwrap();
    let header = writer_header();
    let path = dir.path().join("empty-store.bam");
    let mut writer = BamWriterBuilder::to_path(&path, &header).write_index(true).build().unwrap();
    let store = RecordStore::default();
    let idx = RecordIdx::new(7).unwrap();

    let err = writer.write_store_record(&store, idx).expect_err("an empty store holds no record 7");
    let BamWriteError::NoSuchRecord { idx: reported } = err else {
        panic!("expected NoSuchRecord, got {err:?}");
    };
    assert_eq!(reported, idx, "the error must name the index it was given");
}

// r[verify bam_writer.error_poisoning]
/// A real I/O failure poisons the writer: the first error carries the BGZF
/// write failure, and every later write returns `Poisoned` without reaching
/// the sink at all.
#[test]
fn an_io_failure_poisons_the_writer_and_later_writes_never_reach_the_sink() {
    let header = writer_header();
    // Enough room for the header block, not enough for a flush of record data.
    let sink = FailAfter::new(4096);
    let mut writer = BamWriterBuilder::to_writer(sink.clone(), &header).build().unwrap();

    let mut first_error = None;
    for i in 0..4000u32 {
        let record = mapped_record(0, i * 10 + 1, "read", 100);
        if let Err(e) = writer.write(&record) {
            first_error = Some(e);
            break;
        }
    }
    let first_error = first_error.expect("the sink must run out of room");
    assert!(
        matches!(
            first_error,
            BamWriteError::Bgzf { source: BgzfError::WriteFailed { .. } }
                | BamWriteError::Io { .. }
        ),
        "the first failure must carry the I/O cause, got {first_error:?}",
    );

    let after_failure = sink.written();
    for i in 0..5u32 {
        let record = mapped_record(0, 900_000 + i, "later", 10);
        let err = writer.write(&record).expect_err("a poisoned writer accepts nothing");
        assert!(matches!(err, BamWriteError::Poisoned), "expected Poisoned, got {err:?}");
    }
    let store = RecordStore::default();
    let err = writer
        .write_store_record(&store, RecordIdx::new(0).unwrap())
        .expect_err("a poisoned writer accepts nothing");
    assert!(
        matches!(err, BamWriteError::Poisoned),
        "write_store_record must report Poisoned before NoSuchRecord, got {err:?}",
    );
    assert_eq!(sink.written(), after_failure, "a poisoned writer must not write to the sink");
}

/// The ways a `write()` can fail that this file drives from outside the crate.
#[derive(Debug, Clone, Copy)]
enum WriteFailure {
    MappedWithoutReference,
    RecordTooLarge,
    SeqQualMismatch,
    Unsorted,
}

impl WriteFailure {
    const ALL: [Self; 4] =
        [Self::MappedWithoutReference, Self::RecordTooLarge, Self::SeqQualMismatch, Self::Unsorted];

    fn record(self) -> OwnedBamRecord {
        match self {
            Self::MappedWithoutReference => mapped_record(-1, 100, "orphan", 10),
            // Large enough to pass the 2 MiB limit, small enough to build fast.
            Self::RecordTooLarge => mapped_record(0, 100, "huge", 3_000_000),
            Self::SeqQualMismatch => {
                let mut record =
                    OwnedBamRecord::builder(0, Some(Pos0::new(100).unwrap()), b"shrunk".to_vec())
                        .flags(BamFlags::from(4))
                        .seq(vec![Base::A; 10])
                        .qual(vec![BaseQuality::from_byte(30); 10])
                        .build()
                        .unwrap();
                record.set_seq(vec![Base::C; 4]).unwrap();
                record
            }
            // Position 1 is behind whatever the property already wrote.
            Self::Unsorted => mapped_record(0, 1, "backwards", 10),
        }
    }

    fn matches(self, err: &BamWriteError) -> bool {
        match self {
            Self::MappedWithoutReference => matches!(err, BamWriteError::MappedWithoutReference),
            Self::RecordTooLarge => matches!(err, BamWriteError::RecordTooLarge { .. }),
            Self::SeqQualMismatch => matches!(
                err,
                BamWriteError::Record { source: OwnedRecordError::SeqQualLengthMismatch { .. } }
            ),
            Self::Unsorted => {
                matches!(err, BamWriteError::Index { source: IndexError::UnsortedInput { .. } })
            }
        }
    }
}

// r[verify bam_writer.error_poisoning]
/// Whichever way the first write fails, and however many records preceded it,
/// every subsequent write returns `Poisoned` — never the original error again,
/// and never a success.
#[hegel::test(test_cases = 24)]
fn the_first_failure_wins_and_everything_after_it_is_poisoned(tc: TestCase) {
    let good = tc.draw(gs::integers::<u32>().min_value(1).max_value(6));
    let which = tc.draw(gs::integers::<usize>().max_value(WriteFailure::ALL.len() - 1));
    let failure = WriteFailure::ALL[which];
    let later = tc.draw(gs::integers::<u32>().min_value(1).max_value(4));
    tc.event(format!("{failure:?}"));

    let dir = tempfile::tempdir().unwrap();
    let header = writer_header();
    let path = dir.path().join("poison.bam");
    let mut writer = BamWriterBuilder::to_path(&path, &header).write_index(true).build().unwrap();
    for i in 0..good {
        writer.write(&mapped_record(0, 1_000 + i * 100, "ok", 10)).expect("sorted and valid");
    }

    let err = writer.write(&failure.record()).expect_err("the failing record must be refused");
    assert!(failure.matches(&err), "{failure:?} produced the wrong variant: {err:?}");

    for i in 0..later {
        // A record that would be perfectly acceptable on a healthy writer.
        let record = mapped_record(0, 50_000 + i * 100, "after", 10);
        let err = writer.write(&record).expect_err("the writer is poisoned");
        assert!(
            matches!(err, BamWriteError::Poisoned),
            "write {i} after {failure:?} gave {err:?}, expected Poisoned",
        );
    }
}

// ── AuxData: which error, and whether the block moved ───────────────────

/// An aux block with a few tags already in it, so a setter that wrote before it
/// validated would have something behind it to corrupt.
fn populated_aux() -> AuxData {
    let mut aux = AuxData::new();
    aux.set_string(*b"RG", b"group-1");
    aux.set_int(*b"NM", 3).unwrap();
    aux.set_float(*b"XF", 1.5);
    aux.set_int(*b"AS", -120).unwrap();
    aux
}

/// The two ends of what a BAM aux integer can spell: `c`/`s`/`i` reach down to
/// `i32::MIN`, `C`/`S`/`I` up to `u32::MAX`, and nothing outside that union has
/// a type byte ([SAM1] §4.2.5).
const AUX_INT_MAX: i64 = u32::MAX as i64;
const AUX_INT_MIN: i64 = i32::MIN as i64;

// r[verify bam.owned_record.aux_int_encoding]
// r[verify bam.owned_record.failed_mutation_is_inert]
/// One step past either end of the aux integer range is refused, the error
/// names the value it was handed, and the block is byte-identical afterwards.
#[test]
fn set_int_past_either_end_of_the_aux_range_names_the_value_and_changes_nothing() {
    for (label, value) in [("above u32::MAX", AUX_INT_MAX + 1), ("below i32::MIN", AUX_INT_MIN - 1)]
    {
        let mut aux = populated_aux();
        let before = aux.as_bytes().to_vec();

        let err = aux
            .set_int(*b"ZZ", value)
            .expect_err("a value with no BAM integer type must be refused");
        let AuxDataError::IntegerOutOfRange { value: reported } = err else {
            panic!("{label}: expected IntegerOutOfRange, got {err:?}");
        };
        assert_eq!(reported, value, "{label}: the error must carry the rejected value");
        assert_eq!(aux.as_bytes(), before, "{label}: the aux block must be untouched");
        assert_eq!(aux.get(*b"ZZ"), None, "{label}: no partial tag may be left behind");
    }
}

// r[verify bam.owned_record.aux_int_encoding]
/// The last accepted value at each end really is accepted — without this the
/// test above would still pass if `set_int` rejected the whole range.
#[test]
fn set_int_accepts_both_ends_of_the_aux_range() {
    let mut aux = AuxData::new();
    aux.set_int(*b"HI", AUX_INT_MAX).expect("u32::MAX has type I");
    aux.set_int(*b"LO", AUX_INT_MIN).expect("i32::MIN has type i");
    assert_eq!(aux.get(*b"HI"), Some(AuxValue::U32(u32::MAX)));
    assert_eq!(aux.get(*b"LO"), Some(AuxValue::I32(i32::MIN)));
}

// r[verify bam.owned_record.aux_data]
// r[verify bam.owned_record.failed_mutation_is_inert]
/// `set_char` is the other fallible setter that validates ahead of writing: a
/// byte outside the SAM `A`-type grammar `[!-~]` is refused by value, and the
/// block does not move.
#[test]
fn set_char_outside_the_printable_ascii_grammar_names_the_byte_and_changes_nothing() {
    for value in [0x00u8, 0x20, 0x7f, 0xff] {
        let mut aux = populated_aux();
        let before = aux.as_bytes().to_vec();

        let err = aux.set_char(*b"ZC", value).expect_err("A-type takes printable ASCII only");
        let AuxDataError::InvalidCharByte { value: reported } = err else {
            panic!("byte {value:#04x}: expected InvalidCharByte, got {err:?}");
        };
        assert_eq!(reported, value, "the error must carry the rejected byte");
        assert_eq!(aux.as_bytes(), before, "byte {value:#04x}: the aux block must be untouched");
    }
}

// r[verify bam.owned_record.aux_int_encoding]
// r[verify bam.owned_record.failed_mutation_is_inert]
/// Across the whole i64 domain, `set_int` splits exactly at the aux integer
/// range: inside it the tag round-trips, outside it the call reports the value
/// and leaves the block byte-identical. Testing both directions from one
/// generator is what rules out a setter that simply refuses everything.
#[hegel::test]
fn set_int_is_exactly_as_permissive_as_the_bam_integer_types(tc: TestCase) {
    let value = tc.draw(gs::integers::<i64>());
    let mut aux = populated_aux();
    let before = aux.as_bytes().to_vec();

    let in_range = (AUX_INT_MIN..=AUX_INT_MAX).contains(&value);
    tc.event(if in_range { "in range" } else { "out of range" });

    match aux.set_int(*b"ZZ", value) {
        Ok(()) => {
            assert!(in_range, "{value} is outside the BAM integer types but was accepted");
            let round_tripped = match aux.get(*b"ZZ") {
                Some(AuxValue::U8(v)) => i64::from(v),
                Some(AuxValue::U16(v)) => i64::from(v),
                Some(AuxValue::U32(v)) => i64::from(v),
                Some(AuxValue::I8(v)) => i64::from(v),
                Some(AuxValue::I16(v)) => i64::from(v),
                Some(AuxValue::I32(v)) => i64::from(v),
                other => panic!("set_int stored a non-integer tag: {other:?}"),
            };
            assert_eq!(round_tripped, value, "the stored tag must read back as the value set");
        }
        Err(err) => {
            assert!(!in_range, "{value} fits a BAM integer type but was refused: {err:?}");
            let AuxDataError::IntegerOutOfRange { value: reported } = err else {
                panic!("expected IntegerOutOfRange, got {err:?}");
            };
            assert_eq!(reported, value, "the error must carry the rejected value");
            assert_eq!(aux.as_bytes(), before, "a refused set_int must not move the block");
        }
    }
}
