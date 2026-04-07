//! Full-stack BAM reader fuzz: BGZF stream → BAM header → records → pileup.
//!
//! Input is raw BGZF-compressed BAM bytes (same format as a .bam file).
//! Seeded with a real BAM file for high initial coverage.
#![no_main]

use libfuzzer_sys::fuzz_target;
use seqair::bam::bgzf::BgzfReader;
use seqair::bam::header::BamHeader;
use seqair::bam::pileup::PileupEngine;
use seqair::bam::record_store::RecordStore;
use seqair_types::{Offset, Pos, Zero};

fuzz_target!(|data: &[u8]| {
    if data.len() < 4 || data.len() > 256 * 1024 {
        return;
    }

    // Open BGZF cursor and parse BAM header
    let mut bgzf = BgzfReader::from_cursor(data.to_vec());
    let header = match BamHeader::parse(&mut bgzf) {
        Ok(h) => h,
        Err(_) => return,
    };

    let target_count = header.target_count();
    for i in 0..target_count.min(100) {
        let _ = header.target_name(i as u32);
        let _ = header.target_len(i as u32);
    }

    // Read records through the BGZF layer
    let mut store = RecordStore::new();
    let mut buf = [0u8; 4];
    for _ in 0..128 {
        if bgzf.read_exact_into(&mut buf).is_err() {
            break;
        }
        let block_size = i32::from_le_bytes(buf) as usize;
        if !(32..=100_000).contains(&block_size) {
            break;
        }

        let mut record_buf = vec![0u8; block_size];
        if bgzf.read_exact_into(&mut record_buf).is_err() {
            break;
        }

        let _ = store.push_raw(&record_buf);
    }

    if store.is_empty() {
        return;
    }

    // Pileup over first contig's records
    let region_start = Pos::<Zero>::new(0).unwrap();
    let region_end = match region_start.checked_add_offset(Offset::new(1000)) {
        Some(p) => p,
        None => return,
    };

    let mut engine = PileupEngine::new(store, region_start, region_end);
    engine.set_max_depth(50);
    for col in engine.by_ref().take(1000) {
        let _depth = col.depth();
        let _mdepth = col.match_depth();
        let _refbase = col.reference_base();
        for aln in col.alignments() {
            let _op = aln.op();
            let _base = aln.base();
            let _qual = aln.qual();
        }
    }
});
