//! Criterion benchmarks for CRAM reading, decoding, and pileup.
//!
//! Compares seqair against htslib and noodles.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "benches")]
#![allow(clippy::arithmetic_side_effects, reason = "benches")]
#![allow(clippy::cast_possible_truncation, reason = "benches")]
#![allow(clippy::cast_possible_wrap, reason = "benches")]

use criterion::{Criterion, criterion_group, criterion_main};
use seqair_types::{Base, Pos0};
use std::hint::black_box;

const CRAM_PATH: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.cram");
const CRAI_PATH: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.cram.crai");
const FASTA_PATH: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.fasta.gz");

const CHROM: &str = "chr19";
const START: Pos0 = Pos0::new(0).unwrap();
const END: Pos0 = Pos0::new(6_140_000).unwrap();

// ---------------------------------------------------------------------------
// Group 1: CRAM indexed record decode
// Indexed fetch + decode for a region. Matches bam_record_decode.
// ---------------------------------------------------------------------------

fn cram_record_decode(c: &mut Criterion) {
    use noodles::core::{Position, Region};
    use noodles::cram as ncram;
    use noodles::fasta as nfasta;
    use noodles::sam;
    use rust_htslib::bam::{self, FetchDefinition, Read as _};

    let mut group = c.benchmark_group("cram_record_decode");

    // seqair: IndexedCramReader + RecordStore
    group.bench_function("seqair", |b| {
        b.iter(|| {
            let mut reader = seqair::cram::reader::IndexedCramReader::open(
                std::path::Path::new(CRAM_PATH),
                std::path::Path::new(FASTA_PATH),
            )
            .unwrap();
            let mut store = seqair::bam::RecordStore::new();
            let tid = reader.header().tid(CHROM).unwrap();
            reader.fetch_into(tid, START, END, &mut store).unwrap();
            black_box(store.len())
        });
    });

    // htslib: IndexedReader fetch + read loop
    group.bench_function("htslib", |b| {
        b.iter(|| {
            let mut reader = bam::IndexedReader::from_path(CRAM_PATH).unwrap();
            reader.set_reference(FASTA_PATH).unwrap();
            let tid = reader.header().tid(CHROM.as_bytes()).unwrap();
            reader
                .fetch(FetchDefinition::Region(tid as i32, START.as_i64(), END.as_i64()))
                .unwrap();
            let mut record = bam::Record::new();
            let mut count = 0usize;
            while reader.read(&mut record) == Some(Ok(())) {
                count += 1;
            }
            black_box(count)
        });
    });

    // noodles: CRAM Reader + query
    group.bench_function("noodles", |b| {
        b.iter(|| {
            let fasta_repo = {
                let fasta_reader = nfasta::io::indexed_reader::Builder::default()
                    .build_from_path(FASTA_PATH)
                    .unwrap();
                nfasta::repository::Repository::new(
                    nfasta::repository::adapters::IndexedReader::new(fasta_reader),
                )
            };
            let mut reader = ncram::io::reader::Builder::default()
                .set_reference_sequence_repository(fasta_repo)
                .build_from_path(CRAM_PATH)
                .unwrap();
            let header: sam::Header = reader.read_header().unwrap();
            let crai_index = ncram::crai::fs::read(CRAI_PATH).unwrap();
            let region = Region::new(
                CHROM,
                Position::try_from(1_usize).unwrap()
                    ..=Position::try_from(usize::try_from(END).unwrap()).unwrap(),
            );
            let query = reader.query(&header, &crai_index, &region).unwrap();
            let mut count = 0usize;
            for result in query {
                let _ = result.unwrap();
                count += 1;
            }
            black_box(count)
        });
    });

    group.finish();
}

// ---------------------------------------------------------------------------
// Group 2: CRAM pileup e2e
// End-to-end: open → fetch → decode → pileup → iterate all columns.
// Matches pileup_e2e in bam.rs.
// ---------------------------------------------------------------------------

struct Counter {
    a: u32,
    c: u32,
    g: u32,
    t: u32,
    n: u32,
}

impl Counter {
    fn new() -> Self {
        Self { a: 0, c: 0, g: 0, t: 0, n: 0 }
    }

    fn count(&mut self, base: Base) {
        match base {
            Base::A => self.a += 1,
            Base::C => self.c += 1,
            Base::G => self.g += 1,
            Base::T => self.t += 1,
            Base::Unknown => self.n += 1,
        }
    }
}

fn cram_pileup_e2e(c: &mut Criterion) {
    use rust_htslib::bam::{self, FetchDefinition, Read as _};

    let mut group = c.benchmark_group("cram_pileup_e2e");

    // seqair: fetch_into → PileupEngine → iterate columns
    group.bench_function("seqair", |b| {
        b.iter(|| {
            let mut reader = seqair::cram::reader::IndexedCramReader::open(
                std::path::Path::new(CRAM_PATH),
                std::path::Path::new(FASTA_PATH),
            )
            .unwrap();
            let mut store = seqair::bam::RecordStore::new();
            let tid = reader.header().tid(CHROM).unwrap();
            reader.fetch_into(tid, START, END, &mut store).unwrap();

            let mut engine = seqair::bam::PileupEngine::new(store, START, END);
            let mut total_depth: u64 = 0;
            let mut columns: u64 = 0;
            let mut counter = Counter::new();
            while let Some(col) = engine.pileups() {
                total_depth += col.depth() as u64;
                columns += 1;
                col.alignments().for_each(|aln| counter.count(aln.base().unwrap_or_default()));
            }
            black_box((columns, total_depth, counter))
        });
    });

    // htslib: IndexedReader → fetch → pileup() → iterate columns
    group.bench_function("htslib", |b| {
        b.iter(|| {
            let mut reader = bam::IndexedReader::from_path(CRAM_PATH).unwrap();
            reader.set_reference(FASTA_PATH).unwrap();
            let tid = reader.header().tid(CHROM.as_bytes()).unwrap();
            reader
                .fetch(FetchDefinition::Region(tid as i32, START.as_i64(), END.as_i64()))
                .unwrap();

            let mut total_depth: u64 = 0;
            let mut columns: u64 = 0;
            let mut counter = Counter::new();

            for p in reader.pileup() {
                let p = p.unwrap();
                let pos = u64::from(p.pos());
                if !(START.as_u64()..=END.as_u64()).contains(&pos) {
                    continue;
                }
                total_depth += u64::from(p.depth());
                columns += 1;
                p.alignments().for_each(|aln| {
                    let Some(qpos) = aln.qpos() else {
                        return;
                    };
                    let record = aln.record();
                    counter.count(Base::from(record.seq()[qpos]));
                });
            }
            black_box((columns, total_depth, counter))
        });
    });

    group.finish();
}

// ---------------------------------------------------------------------------
// Group 3: CRAM sequential full-file decode
// Read all records without index filtering. Tests pure decode throughput.
// seqair uses indexed fetch over the full chromosome; htslib and noodles use
// their non-indexed readers.
// ---------------------------------------------------------------------------

fn cram_full_decode(c: &mut Criterion) {
    use noodles::core::Region;
    use noodles::cram as ncram;
    use noodles::fasta as nfasta;
    use noodles::sam;
    use rust_htslib::bam::{self, FetchDefinition, Read as _};

    let mut group = c.benchmark_group("cram_full_decode");

    // seqair: IndexedCramReader + fetch_into for full chr19
    group.bench_function("seqair", |b| {
        b.iter(|| {
            let mut reader = seqair::cram::reader::IndexedCramReader::open(
                std::path::Path::new(CRAM_PATH),
                std::path::Path::new(FASTA_PATH),
            )
            .unwrap();
            let tid = reader.header().tid(CHROM).unwrap();
            let chr_len = reader.header().target_len(tid).unwrap();
            let full_start = Pos0::new(0).unwrap();
            let full_end = Pos0::new(u32::try_from(chr_len).unwrap()).unwrap();
            let mut store = seqair::bam::RecordStore::new();
            reader.fetch_into(tid, full_start, full_end, &mut store).unwrap();
            black_box(store.len())
        });
    });

    // htslib: IndexedReader + full chr19 range
    group.bench_function("htslib", |b| {
        b.iter(|| {
            let mut reader = bam::IndexedReader::from_path(CRAM_PATH).unwrap();
            reader.set_reference(FASTA_PATH).unwrap();
            let tid = reader.header().tid(CHROM.as_bytes()).unwrap();
            // Fetch full chr19: 0 to target length
            let chr_len = reader.header().target_len(tid).unwrap();
            reader.fetch(FetchDefinition::Region(tid as i32, 0, chr_len as i64)).unwrap();
            let mut record = bam::Record::new();
            let mut count = 0usize;
            while reader.read(&mut record) == Some(Ok(())) {
                count += 1;
            }
            black_box(count)
        });
    });

    // noodles: indexed query + full chr19 range
    group.bench_function("noodles", |b| {
        b.iter(|| {
            let fasta_repo = {
                let fasta_reader = nfasta::io::indexed_reader::Builder::default()
                    .build_from_path(FASTA_PATH)
                    .unwrap();
                nfasta::repository::Repository::new(
                    nfasta::repository::adapters::IndexedReader::new(fasta_reader),
                )
            };
            let mut reader = ncram::io::reader::Builder::default()
                .set_reference_sequence_repository(fasta_repo)
                .build_from_path(CRAM_PATH)
                .unwrap();
            let header: sam::Header = reader.read_header().unwrap();
            let crai_index = ncram::crai::fs::read(CRAI_PATH).unwrap();
            let region = Region::new(CHROM, ..);
            let query = reader.query(&header, &crai_index, &region).unwrap();
            let mut count = 0usize;
            for result in query {
                let _ = result.unwrap();
                count += 1;
            }
            black_box(count)
        });
    });

    group.finish();
}

criterion_group!(benches, cram_record_decode, cram_pileup_e2e, cram_full_decode,);
criterion_main!(benches);
