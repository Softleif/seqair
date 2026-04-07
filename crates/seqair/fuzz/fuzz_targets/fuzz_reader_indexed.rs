//! Full-stack fuzz through FuzzReaders: BAM/SAM/CRAM + index + FASTA → fetch → pileup.
//!
//! Uses a simple hand-rolled binary format (not Arbitrary) for zero-waste seed usage.
//! First byte selects BAM (0), SAM-BGZF (1), CRAM (2), or plain SAM (3), then
//! fixed-size length headers split the remaining bytes into alignment data, index,
//! FASTA, FAI, and GZI.
#![no_main]

use libfuzzer_sys::fuzz_target;
use seqair::{
    bam::{pileup::PileupEngine, record_store::RecordStore},
    reader::FuzzReaders,
};
use seqair_fuzz::indexed_reader::{Format, Input};
use seqair_types::{Offset, Pos, Zero};

fuzz_target!(|data: &[u8]| {
    let Some(input) = Input::parse(data) else {
        return;
    };
    let Ok(fai_str) = std::str::from_utf8(input.fai) else {
        return;
    };

    let readers = match input.format {
        Format::Bam => FuzzReaders::from_bam_bytes(
            input.data1.to_vec(),
            input.data2,
            input.fasta_gz.to_vec(),
            fai_str,
            input.gzi,
        ),
        Format::Sam => FuzzReaders::from_sam_bytes(
            input.data1.to_vec(),
            input.data2,
            input.fasta_gz.to_vec(),
            fai_str,
            input.gzi,
        ),
        Format::Cram => FuzzReaders::from_cram_bytes(
            input.data1.to_vec(),
            input.data2,
            input.fasta_gz.to_vec(),
            fai_str,
            input.gzi,
        ),
        Format::PlainSam => FuzzReaders::from_plain_sam_bytes(
            input.data1.to_vec(),
            input.fasta_gz.to_vec(),
            fai_str,
            input.gzi,
        ),
    };

    let Ok(mut readers) = readers else {
        return;
    };

    let header = readers.header();
    let target_count = header.target_count();
    if target_count == 0 {
        return;
    }
    for i in 0..target_count.min(50) {
        let _ = header.target_name(i as u32);
        let _ = header.target_len(i as u32);
    }

    let start = Pos::<Zero>::new(0).unwrap();
    let Some(end) = start.checked_add_offset(Offset::new(10_000)) else {
        return;
    };

    let mut store = RecordStore::new();
    let _ = readers.fetch_into(0, start, end, &mut store);

    if store.is_empty() {
        return;
    }

    let mut engine = PileupEngine::new(store, start, end);
    engine.set_max_depth(50);
    for col in engine.by_ref().take(1000) {
        let _depth = col.depth();
        let _mdepth = col.match_depth();
        for aln in col.alignments() {
            let _op = aln.op();
            let _base = aln.base();
            let _qual = aln.qual();
        }
    }
});
