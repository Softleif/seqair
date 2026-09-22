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
use seqair::bam::cigar::{CigarOp, CigarOpType};
use seqair::bam::header::BamHeader;
use seqair::bam::owned_record::OwnedBamRecord;
use seqair::bam::writer::BamWriterBuilder;
use seqair::bam::{BamError, BamHeaderError, BgzfError, IndexedBamReader, Pos0, RecordStore};
use seqair_types::bam_flags::BamFlags;
use seqair_types::{Base, BaseQuality};
use std::collections::BTreeMap;
use std::sync::OnceLock;

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
