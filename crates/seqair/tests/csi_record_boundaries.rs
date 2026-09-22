//! Every virtual offset a co-produced CSI hands a reader must be the first byte
//! of a record.
//!
//! A dense BCF is written through the unified `Writer`, its CSI is queried the
//! way a reader queries it, and each chunk start is checked against the set of
//! record starts recovered by decoding the file sequentially. The offsets are
//! only ever as good as the BGZF writer's `virtual_offset()`, so this is the
//! end-to-end form of `r[bgzf.writer.virtual_offset]`: before that rule was
//! implemented correctly, a record ending exactly on a 64 KiB block boundary
//! made the writer report the block's last byte, and htslib rejected the file
//! with "shared section malformed or too short" for any region query that
//! entered the affected bin.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "test code"
)]
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    reason = "test code with known small values"
)]

use seqair::bam::bgzf::VirtualOffset;
use seqair::bam::{CsiIndex, Pos0};
use seqair::vcf::record_encoder::{InfoFieldDef, InfoString, Str};
use seqair::vcf::{Alleles, ContigDef, Number, OutputFormat, ValueType, VcfHeader, Writer};
use seqair_types::{Base, Pos1};
use std::collections::BTreeSet;
use std::sync::Arc;

const CONTIG_LEN: u32 = 20_000_000;
/// Enough records to fill several BGZF blocks.
const N_RECORDS: u32 = 900;
/// One record per 16 kb bin, so *every* record boundary is also a bin boundary
/// and becomes the next bin's chunk start. A file where several records share a
/// bin hides the bug: the bad offset only reaches the index when the record that
/// ends on the block boundary is also the last one of its bin.
const STRIDE: u32 = 16_384;

struct Written {
    bcf: Vec<u8>,
    csi: Vec<u8>,
}

/// Write a BCF whose records each carry a `pad`-byte INFO string, so the sweep
/// below walks record sizes past every residue mod 64 KiB.
fn write_bcf(pad: usize) -> Written {
    let mut builder = VcfHeader::builder();
    let contig =
        builder.register_contig("chr1", ContigDef { length: Some(CONTIG_LEN.into()) }).unwrap();
    let mut infos = builder.infos();
    let sc: InfoString = infos
        .register_info(&InfoFieldDef::<Str>::new(
            "SC",
            Number::Count(1),
            ValueType::String,
            "padding",
        ))
        .unwrap();
    let header = Arc::new(infos.build().unwrap());

    let padding = "A".repeat(pad);
    let alleles = Alleles::snv(Base::C, Base::T).unwrap();

    let mut bcf = Vec::new();
    let mut writer = Writer::new(&mut bcf, OutputFormat::Bcf).write_header(&header).unwrap();
    for i in 0..N_RECORDS {
        let pos = Pos1::new(i * STRIDE + 1).unwrap();
        let mut enc =
            writer.begin_record(&contig, pos, &alleles, Some(30.0)).unwrap().filter_pass();
        sc.encode(&mut enc, &padding);
        enc.emit().unwrap();
    }
    let (_inner, index) = writer.finish().unwrap();

    let mut csi = Vec::new();
    index.expect("BCF output is indexed").write(&mut csi).unwrap();
    Written { bcf, csi }
}

/// One BGZF block: where it starts in the file, and where its first byte sits
/// in the uncompressed stream.
struct Block {
    compressed_offset: u64,
    uncompressed_start: usize,
    uncompressed_len: usize,
}

fn blocks_of(data: &[u8]) -> Vec<Block> {
    let mut blocks = Vec::new();
    let mut pos = 0usize;
    let mut uncompressed = 0usize;
    while pos + 18 <= data.len() {
        let bsize = u16::from_le_bytes([data[pos + 16], data[pos + 17]]) as usize + 1;
        let len =
            u32::from_le_bytes(data[pos + bsize - 4..pos + bsize].try_into().unwrap()) as usize;
        blocks.push(Block {
            compressed_offset: pos as u64,
            uncompressed_start: uncompressed,
            uncompressed_len: len,
        });
        uncompressed += len;
        pos += bsize;
    }
    blocks
}

/// Resolve a virtual offset to a position in the uncompressed stream, rejecting
/// one that names a block the file does not have or a byte past the block's end.
fn resolve(blocks: &[Block], voff: VirtualOffset) -> Option<usize> {
    let block = blocks.iter().find(|b| b.compressed_offset == voff.block_offset())?;
    let within = usize::from(voff.within_block());
    // == len is allowed only for the EOF block's zero length, which no chunk names.
    (within <= block.uncompressed_len).then_some(block.uncompressed_start + within)
}

/// Uncompressed byte offset of every record, by decoding the BCF sequentially.
fn record_starts(bcf: &[u8]) -> BTreeSet<usize> {
    let mut plain = Vec::new();
    let mut reader = seqair::bam::bgzf::BgzfReader::from_reader(std::io::Cursor::new(bcf));
    reader.read_to_end(&mut plain).unwrap();

    assert_eq!(&plain[..3], b"BCF", "not a BCF stream");
    let l_text = u32::from_le_bytes(plain[5..9].try_into().unwrap()) as usize;
    let mut pos = 9 + l_text;

    let mut starts = BTreeSet::new();
    while pos + 8 <= plain.len() {
        starts.insert(pos);
        let l_shared = u32::from_le_bytes(plain[pos..pos + 4].try_into().unwrap()) as usize;
        let l_indiv = u32::from_le_bytes(plain[pos + 4..pos + 8].try_into().unwrap()) as usize;
        pos += 8 + l_shared + l_indiv;
    }
    // r[verify bcf_writer.record_layout]
    // Walking the body purely by `l_shared` + `l_indiv` has to land exactly on
    // the end of the stream: any record whose two length prefixes disagree with
    // the bytes that follow desynchronises the walk and this fails.
    assert_eq!(pos, plain.len(), "BCF body did not decode to a whole number of records");
    assert_eq!(starts.len(), N_RECORDS as usize, "unexpected record count");
    starts
}

/// Check every offset the index would hand a reader for every bin-sized region.
/// Returns the number of chunk starts checked.
fn assert_chunk_starts_are_records(w: &Written, pad: usize) -> usize {
    let blocks = blocks_of(&w.bcf);
    let starts = record_starts(&w.bcf);
    let dir = tempfile::tempdir().unwrap();
    let csi_path = dir.path().join("sweep.bcf.csi");
    std::fs::write(&csi_path, &w.csi).unwrap();
    let index = CsiIndex::from_path(&csi_path).unwrap();

    let mut checked = 0usize;
    let last = N_RECORDS * STRIDE;
    let mut region_start = 0u32;
    while region_start <= last {
        let region_end = region_start + 16_383;
        for chunk in
            index.query(0, Pos0::new(region_start).unwrap(), Pos0::new(region_end).unwrap())
        {
            let resolved = resolve(&blocks, chunk.begin).unwrap_or_else(|| {
                panic!(
                    "pad {pad}: chunk start {:?} (block {}, within {}) is not a position in the file",
                    chunk.begin,
                    chunk.begin.block_offset(),
                    chunk.begin.within_block()
                )
            });
            assert!(
                starts.contains(&resolved),
                "pad {pad}: chunk start {:?} resolves to uncompressed byte {resolved}, \
                 which is inside a record, not at one — a reader seeking there reads garbage",
                chunk.begin
            );
            checked += 1;
        }
        region_start += 16_384;
    }
    assert!(checked > 0, "pad {pad}: the sweep queried no chunks at all");
    checked
}

// r[verify bgzf.writer.virtual_offset]
// r[verify csi.write_loffset]
#[test]
fn csi_offsets_always_name_a_record_start() {
    // Whether a record lands exactly on a block boundary is a question of record
    // size against 65536, so one size proves nothing: sweep them. The assertion
    // at the end is what keeps the sweep honest if the record shape changes.
    let mut saw_exactly_full_block = false;
    for pad in 0..96 {
        let w = write_bcf(pad);
        assert_chunk_starts_are_records(&w, pad);
        saw_exactly_full_block |= blocks_of(&w.bcf).iter().any(|b| b.uncompressed_len == 65_536);
    }
    // r[verify bcf_writer.bgzf_blocks]
    // Also pins the flush threshold: a writer that flushed early (say every
    // 1 KiB) would still produce a readable file, but no block would ever reach
    // 64 KiB and this assertion would fail.
    assert!(
        saw_exactly_full_block,
        "no file in the sweep had a block filled to exactly 64 KiB — the case this test exists \
         for was never exercised; widen the sweep or change the record shape"
    );
}

/// The first record size in the sweep whose file has a block filled to exactly
/// 64 KiB — the shape that used to corrupt the index.
fn block_aligned_pad() -> usize {
    (0..96)
        .find(|&pad| blocks_of(&write_bcf(pad).bcf).iter().any(|b| b.uncompressed_len == 65_536))
        .expect("no record size in 0..96 fills a block exactly")
}

// r[verify bgzf.writer.virtual_offset]
/// The symptom as it was reported: `bcftools view -R` over a dense BCF failed
/// with "shared section malformed or too short" while a sequential read of the
/// same file was fine, because the index — not the data — was wrong.
#[test]
fn bcftools_can_region_query_a_block_aligned_bcf() {
    if std::process::Command::new("bcftools").arg("--version").output().is_err() {
        eprintln!("skipping: bcftools not found");
        return;
    }

    let w = write_bcf(block_aligned_pad());
    let dir = tempfile::tempdir().unwrap();
    let bcf_path = dir.path().join("aligned.bcf");
    std::fs::write(&bcf_path, &w.bcf).unwrap();
    std::fs::write(dir.path().join("aligned.bcf.csi"), &w.csi).unwrap();

    // One BED interval per bin, which is how the failure was first seen.
    let bed_path = dir.path().join("regions.bed");
    let mut bed = String::new();
    for i in 0..N_RECORDS {
        use std::fmt::Write as _;
        writeln!(bed, "chr1\t{}\t{}", i * STRIDE, (i + 1) * STRIDE).unwrap();
    }
    std::fs::write(&bed_path, bed).unwrap();

    let out = std::process::Command::new("bcftools")
        .args(["view", "-H", "-R"])
        .arg(&bed_path)
        .arg(&bcf_path)
        .output()
        .expect("bcftools should run");

    assert!(
        out.status.success(),
        "bcftools rejected seqair's index: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert_eq!(
        String::from_utf8_lossy(&out.stdout).lines().count(),
        N_RECORDS as usize,
        "region query lost records"
    );
}
