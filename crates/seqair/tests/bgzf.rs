//! Tests for BGZF virtual offset representation and ordering.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::indexing_slicing)]
use proptest::prelude::*;
use seqair::bam::{BamError, BamHeaderError, BgzfError, bgzf::VirtualOffset};

// r[verify bgzf.virtual_offset]
// Verifies the BAI spec (SAM1 §4.1): upper 48 bits = compressed block offset,
// lower 16 bits = uncompressed offset within block.
proptest! {
    #[test]
    fn virtual_offset_bit_layout_matches_spec(
        block_offset in 0u64..=(1u64 << 48) - 1,
        within_block in 0u16..=u16::MAX,
    ) {
        let vo = VirtualOffset::new(block_offset, within_block);
        let expected_raw = (block_offset << 16) | u64::from(within_block);
        prop_assert_eq!(vo.0, expected_raw, "raw u64 must match BAI spec bit layout");
        prop_assert_eq!(vo.block_offset(), block_offset);
        prop_assert_eq!(vo.within_block(), within_block);
    }
}

// r[verify bgzf.virtual_offset]
// Verifies that VirtualOffsets constructed from raw BAI u64 values sort in the
// same order as their raw u64 representations, and that block_offset is the
// primary sort key per the BAI spec.
proptest! {
    #[test]
    fn virtual_offset_ordering_consistent_with_raw_u64(
        raw1 in 0u64..=u64::MAX,
        raw2 in 0u64..=u64::MAX,
    ) {
        let vo1 = VirtualOffset(raw1);
        let vo2 = VirtualOffset(raw2);
        // Ordering must be identical to the raw u64 ordering (since the struct
        // is a transparent u64 newtype and block_offset occupies the high bits).
        prop_assert_eq!(vo1.cmp(&vo2), raw1.cmp(&raw2));
        // When vo1 < vo2, block_offset of vo1 must be <= block_offset of vo2,
        // confirming block_offset is the primary sort key.
        if vo1 < vo2 {
            prop_assert!(vo1.block_offset() <= vo2.block_offset());
        }
    }
}

// r[verify io.errors]
#[test]
fn error_types_implement_std_error_and_capture_context() {
    use std::error::Error;

    let err = BamHeaderError::InvalidMagic;
    let _: &dyn Error = &err;
    assert!(format!("{err}").contains("BAM"));

    let err = BamError::ContigNotFound { name: seqair_types::SmolStr::new("chr99") };
    assert!(format!("{err}").contains("chr99"));

    let err = BamError::RegionOutOfBounds {
        contig: seqair_types::SmolStr::new("chr1"),
        start: 100,
        end: 200,
        contig_len: 50,
    };
    let msg = format!("{err}");
    assert!(msg.contains("chr1") && msg.contains("100") && msg.contains("50"));
}

// r[verify bgzf.crc32]
// r[verify bgzf.fast_header]
// r[verify bgzf.resize_uninit]
// r[verify bgzf.block_offset_tracking]
#[test]
fn read_bam_verifies_crc32_and_tracks_offsets() {
    use seqair::bam::reader::IndexedBamReader;
    use seqair::bam::record_store::RecordStore;
    use std::path::Path;

    let bam_path = Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.bam"));
    let mut reader = IndexedBamReader::open(bam_path).expect("open BAM");
    let tid = reader.header().tid("chr19").expect("tid");

    let mut store = RecordStore::new();
    use seqair::bam::{Pos, Zero};
    let count = reader
        .fetch_into(
            tid,
            Pos::<Zero>::new(6_105_700).unwrap(),
            Pos::<Zero>::new(6_105_800).unwrap(),
            &mut store,
        )
        .expect("CRC32 verification should pass on valid BAM data");

    // If CRC32 verification failed, fetch_into would have returned an error.
    // If block_offset tracking broke, seeking via the index would fail.
    // If resize_uninit corrupted data, record decoding would fail.
    // If fast_header parsing failed, the block wouldn't decompress.
    assert!(count > 0, "should have read some records");
}

// r[verify bgzf.crc32]
#[test]
fn checksum_mismatch_error_formats_correctly() {
    let err = BgzfError::ChecksumMismatch { expected: 0xDEADBEEF, found: 0xCAFEBABE };
    let msg = format!("{err}");
    assert!(msg.contains("0xdeadbeef"), "should contain expected CRC: {msg}");
    assert!(msg.contains("0xcafebabe"), "should contain found CRC: {msg}");
}

// r[verify io.named_constants]
#[test]
fn bam_flag_constants_match_sam_spec() {
    use seqair::bam::flags::*;
    assert_eq!(FLAG_PAIRED, 0x1);
    assert_eq!(FLAG_PROPER_PAIR, 0x2);
    assert_eq!(FLAG_UNMAPPED, 0x4);
    assert_eq!(FLAG_MATE_UNMAPPED, 0x8);
    assert_eq!(FLAG_REVERSE, 0x10);
    assert_eq!(FLAG_MATE_REVERSE, 0x20);
    assert_eq!(FLAG_FIRST_IN_TEMPLATE, 0x40);
    assert_eq!(FLAG_SECOND_IN_TEMPLATE, 0x80);
    assert_eq!(FLAG_SECONDARY, 0x100);
    assert_eq!(FLAG_FAILED_QC, 0x200);
    assert_eq!(FLAG_DUPLICATE, 0x400);
    assert_eq!(FLAG_SUPPLEMENTARY, 0x800);
}

// r[verify io.typed_flags]
#[test]
fn bam_flags_newtype_predicate_methods() {
    use seqair::bam::flags::BamFlags;

    // 0x53 = paired(0x1) + proper_pair(0x2) + reverse(0x10) + first_in_template(0x40)
    let flags = BamFlags::new(0x53);
    assert!(flags.is_paired());
    assert!(flags.is_proper_pair());
    assert!(flags.is_reverse());
    assert!(flags.is_first_in_template());
    assert!(!flags.is_unmapped());
    assert!(!flags.is_secondary());
    assert!(!flags.is_supplementary());

    let unmapped = BamFlags::new(0x4);
    assert!(unmapped.is_unmapped());
    assert!(!unmapped.is_reverse());

    assert_eq!(flags.raw(), 0x53);
}
