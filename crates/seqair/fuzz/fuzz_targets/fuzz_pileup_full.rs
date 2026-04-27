#![no_main]

use arbitrary::Arbitrary;
use libfuzzer_sys::fuzz_target;
use seqair::bam::pileup::PileupEngine;
use seqair::bam::record_store::RecordStore;
use seqair_types::Pos0;

/// Fuzz input: multiple raw BAM records + a region to pileup over.
#[derive(Arbitrary, Debug)]
struct PileupInput {
    /// Raw BAM record bytes (multiple concatenated with length prefixes).
    records: Vec<Vec<u8>>,
    /// Region start (clamped to reasonable range).
    region_start: u32,
    /// Region length (clamped to prevent long iteration).
    region_len: u16,
}

fuzz_target!(|input: PileupInput| {
    // Limit records and region to prevent timeouts
    let max_records = 32;
    let max_region_len = 1000;

    let region_start = match Pos0::new(input.region_start) {
        Some(p) => p,
        None => return,
    };
    let region_len = (input.region_len as u32).min(max_region_len);
    let region_end =
        match region_start.checked_add_offset(seqair_types::Offset::new(region_len as i64)) {
            Some(p) => p,
            None => return,
        };

    let mut store = RecordStore::new();

    // Push up to max_records into the store
    for raw in input.records.iter().take(max_records) {
        let _ = store.push_raw(raw, &mut ());
    }

    if store.is_empty() {
        return;
    }

    // Create pileup engine and iterate all columns
    let mut engine = PileupEngine::new(store, region_start, region_end);
    engine.set_max_depth(100);

    let mut columns = 0u32;
    while let Some(col) = engine.pileups() {
        // Exercise all accessors on the column
        let _pos = col.pos();
        let _depth = col.depth();
        let _mdepth = col.match_depth();
        let _refbase = col.reference_base();

        for aln in col.alignments() {
            let _idx = aln.record_idx();
            let _op = aln.op();
            let _qpos = aln.qpos();
            let _base = aln.base();
            let _qual = aln.qual();
            let _del = aln.is_del();
            let _skip = aln.is_refskip();
            let _ilen = aln.insert_len();
            let _dlen = aln.del_len();
        }
        columns = columns.saturating_add(1);
        if columns > max_region_len {
            break;
        }
    }
});
