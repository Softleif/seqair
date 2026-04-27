//! Full Readers::pileup pipeline: BAM + BAI + FASTA + FAI + GZI → pileup with reference.
//!
//! This is the highest-level fuzz target — it exercises the exact same code path
//! that rastair uses in production: open readers, fetch a region, iterate pileup
//! columns with reference bases.
//!
//! Seeded with real test files concatenated with a split header.
#![no_main]

use libfuzzer_sys::fuzz_target;
use seqair::reader::FuzzReaders;
use seqair_types::{Offset, Pos0};

fuzz_target!(|data: &[u8]| {
    // Layout: [bam_len:4][bai_len:4][fai_len:2][gzi_len:2][bam][bai][fasta_gz][fai_text][gzi]
    if data.len() < 12 {
        return;
    }

    let bam_len = u32::from_le_bytes([data[0], data[1], data[2], data[3]]) as usize;
    let bai_len = u32::from_le_bytes([data[4], data[5], data[6], data[7]]) as usize;
    let fai_len = u16::from_le_bytes([data[8], data[9]]) as usize;
    let gzi_len = u16::from_le_bytes([data[10], data[11]]) as usize;

    let header_size = 12;
    let total = header_size + bam_len + bai_len + fai_len + gzi_len;
    if data.len() < total || bam_len < 4 || bai_len < 8 {
        return;
    }

    let mut off = header_size;
    let bam_data = &data[off..off + bam_len];
    off += bam_len;
    let bai_data = &data[off..off + bai_len];
    off += bai_len;
    // Everything between fai and gzi is fasta_gz
    let fasta_end = data.len() - fai_len - gzi_len;
    let fasta_gz_data = &data[off..fasta_end.max(off)];
    let fai_text_bytes = &data[fasta_end.max(off)
        ..fasta_end.max(off) + fai_len.min(data.len().saturating_sub(fasta_end.max(off)))];
    let gzi_data = &data[data.len().saturating_sub(gzi_len)..];

    let fai_text = match std::str::from_utf8(fai_text_bytes) {
        Ok(s) => s,
        Err(_) => return,
    };

    let mut readers = match FuzzReaders::from_bam_bytes(
        bam_data.to_vec(),
        bai_data,
        fasta_gz_data.to_vec(),
        fai_text,
        gzi_data,
    ) {
        Ok(r) => r,
        Err(_) => return,
    };

    let target_count = readers.header().target_count();
    if target_count == 0 {
        return;
    }

    let start = Pos0::new(0).unwrap();
    let end = match start.checked_add_offset(Offset::new(5_000)) {
        Some(p) => p,
        None => return,
    };

    let mut engine = match readers.pileup(0, start, end) {
        Ok(e) => e,
        Err(_) => return,
    };

    // Iterate pileup — the production hot path
    let mut col_count: usize = 0;
    while let Some(col) = engine.pileups() {
        let _pos = col.pos();
        let _depth = col.depth();
        let _mdepth = col.match_depth();
        let _refbase = col.reference_base();
        for aln in col.alignments() {
            let _op = aln.op();
            let _base = aln.base();
            let _qual = aln.qual();
            let _qpos = aln.qpos();
        }
        col_count = col_count.saturating_add(1);
        if col_count >= 1000 {
            break;
        }
    }

    readers.recover_store(&mut engine);
});
