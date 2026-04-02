//! Full-stack fuzz through IndexedReader<Cursor>: BAM+BAI bytes → IndexedReader →
//! header access → fetch_into → pileup.
//!
//! Exercises the unified reader enum with cursor-backed I/O — the same dispatch
//! path production code uses, just without file I/O.
#![no_main]

use arbitrary::Arbitrary;
use libfuzzer_sys::fuzz_target;
use seqair::{
    bam::{pileup::PileupEngine, record_store::RecordStore},
    reader::FuzzReaders,
};
use seqair_types::{Offset, Pos, Zero};

#[derive(Arbitrary, Debug)]
struct Input {
    readers: ReadersInput,
    alignment: FastaInput,
}

#[derive(Arbitrary, Debug)]
enum ReadersInput {
    Bam { bam: Vec<u8>, bai: Vec<u8> },
    Sam { sam: Vec<u8>, sai: Vec<u8> },
    Cram { cram: Vec<u8>, crai: Vec<u8> },
}

#[derive(Arbitrary, Debug)]
struct FastaInput {
    fasta_gz: Vec<u8>,
    fai: String,
    gzi: Vec<u8>,
}

fuzz_target!(|data: Input| {
    run(data);
});

fn run(data: Input) {
    let FastaInput { fasta_gz, fai, gzi } = data.alignment;
    let reader = match data.readers {
        ReadersInput::Bam { bam, bai } => {
            FuzzReaders::from_bam_bytes(bam, &bai, fasta_gz, &fai, &gzi)
        }
        ReadersInput::Sam { sam, sai } => return,
        ReadersInput::Cram { cram, crai } => {
            FuzzReaders::from_cram_bytes(cram, &crai, fasta_gz, &fai, &gzi)
        }
    };
    let Ok(mut reader) = reader else {
        return;
    };

    // Exercise header through the IndexedReader dispatch
    let header = reader.header();
    let target_count = header.target_count();
    if target_count == 0 {
        return;
    }
    for i in 0..target_count.min(50) {
        let _ = header.target_name(i as u32);
        let _ = header.target_len(i as u32);
    }

    // fetch_into through IndexedReader dispatch
    let start = Pos::<Zero>::new(0).unwrap();
    let Some(end) = start.checked_add_offset(Offset::new(10_000)) else {
        return;
    };

    let mut store = RecordStore::new();
    let _ = reader.fetch_into(0, start, end, &mut store);

    if store.len() == 0 {
        return;
    }

    // Pileup
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
}
