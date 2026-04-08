//! Cross-validation: seqair BGZF virtual offsets and decompression against noodles-bgzf.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    reason = "test code"
)]
use noodles_bgzf::VirtualPosition as NoodlesVP;
use proptest::prelude::*;
use seqair::bam::bgzf::{BgzfReader, VirtualOffset};

fn test_bam_path() -> &'static std::path::Path {
    std::path::Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.bam"))
}

// r[verify bgzf.virtual_offset]
// Both crates encode the BAI spec bit layout identically: upper 48 bits = compressed block
// offset, lower 16 bits = uncompressed within-block offset. Assert raw u64 and field accessors
// agree across both implementations.
proptest! {
    #[test]
    fn virtual_offset_encoding_matches_noodles(
        block_offset in 0u64..=(1u64 << 48) - 1,
        within_block in 0u16..=u16::MAX,
    ) {
        let vo = VirtualOffset::new(block_offset, within_block);
        // noodles VirtualPosition::new returns Option; the range constraint guarantees Some.
        let noodles_vp = NoodlesVP::new(block_offset, within_block)
            .expect("block_offset is within 48-bit range");

        prop_assert_eq!(
            vo.0,
            u64::from(noodles_vp),
            "raw u64 must agree: seqair={:#x}, noodles={:#x}",
            vo.0,
            u64::from(noodles_vp),
        );
        prop_assert_eq!(
            vo.block_offset(),
            noodles_vp.compressed(),
            "block_offset / compressed must agree",
        );
        prop_assert_eq!(
            vo.within_block(),
            noodles_vp.uncompressed(),
            "within_block / uncompressed must agree",
        );
    }
}

// r[verify bgzf.virtual_offset]
// Both crates must produce the same Ord result for arbitrary raw u64 virtual offsets.
proptest! {
    #[test]
    fn virtual_offset_ordering_matches_noodles(
        raw1 in 0u64..=u64::MAX,
        raw2 in 0u64..=u64::MAX,
    ) {
        let vo1 = VirtualOffset(raw1);
        let vo2 = VirtualOffset(raw2);
        let nv1 = NoodlesVP::from(raw1);
        let nv2 = NoodlesVP::from(raw2);

        prop_assert_eq!(
            vo1.cmp(&vo2),
            nv1.cmp(&nv2),
            "ordering must agree for raw1={:#x}, raw2={:#x}",
            raw1, raw2,
        );
    }
}

// r[verify bgzf.decompression]
// r[verify bgzf.crc32]
// Both seqair's BgzfReader and noodles' bgzf::io::Reader must decompress the same leading
// bytes from a real BAM file. BGZF decompression is deterministic, so byte-for-byte identity
// is required.
#[test]
fn bgzf_decompression_matches_noodles() {
    use std::io::Read;

    let path = test_bam_path();
    const TARGET: usize = 4096;

    // Decompress the first TARGET bytes with seqair's BgzfReader.
    let mut seqair_reader = BgzfReader::open(path).expect("open BgzfReader");
    let mut seqair_buf = Vec::with_capacity(TARGET);
    let mut tmp = [0u8; TARGET];
    while seqair_buf.len() < TARGET {
        let n = seqair_reader.read_up_to(&mut tmp).expect("seqair read_up_to");
        if n == 0 {
            break;
        }
        seqair_buf.extend_from_slice(&tmp[..n]);
    }
    let read_len = seqair_buf.len().min(TARGET);
    let seqair_buf = &seqair_buf[..read_len];

    // Decompress the same bytes with noodles' bgzf::io::Reader.
    let file = std::fs::File::open(path).expect("open BAM for noodles");
    let mut noodles_reader = noodles_bgzf::io::Reader::new(file);
    let mut noodles_buf = vec![0u8; read_len];
    noodles_reader.read_exact(&mut noodles_buf).expect("noodles read_exact");

    assert_eq!(
        seqair_buf,
        noodles_buf.as_slice(),
        "decompressed bytes must be identical between seqair and noodles",
    );
}
