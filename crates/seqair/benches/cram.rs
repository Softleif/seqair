//! Criterion benchmarks for CRAM reading, decoding, and pileup.
//!
//! Compares seqair against htslib and noodles.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "benches")]
#![allow(clippy::arithmetic_side_effects, reason = "benches")]
#![allow(clippy::cast_possible_truncation, reason = "benches")]
#![allow(clippy::cast_possible_wrap, reason = "benches")]

use core::range::RangeInclusive;
use criterion::{Criterion, criterion_group, criterion_main};
use seqair_types::{Base, Pos0};
use std::hint::black_box;

const CRAM_PATH: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.cram");
const CRAI_PATH: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.cram.crai");
const FASTA_PATH: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.fasta.gz");

const CHROM: &str = "chr19";
const START: Pos0 = Pos0::new(0).unwrap();
const END: Pos0 = Pos0::new(6_140_000).unwrap();
/// The closed region every `fetch_into` below asks for.
const SPAN: RangeInclusive<Pos0> = RangeInclusive { start: START, last: END };

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
            reader.fetch_into(tid, SPAN, &mut store).unwrap();
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
                Position::try_from(1_usize).unwrap()..=Position::try_from(END.as_usize()).unwrap(),
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
            reader.fetch_into(tid, SPAN, &mut store).unwrap();

            let mut engine = seqair::bam::PileupEngine::new(store.prepare_for_pileup().input, SPAN);
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
            reader.fetch_into(tid, (full_start..=full_end).into(), &mut store).unwrap();
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

// ---------------------------------------------------------------------------
// Group 4: the arithmetic coder (method 6)
// hts-specs range vectors decoded by seqair and by htscodecs'
// `arith_uncompress_to` (in the htslib rust-htslib links); throughput is
// decoded bytes.
// ---------------------------------------------------------------------------

#[path = "support/htscodecs.rs"]
#[allow(dead_code, reason = "each bench target uses a subset")]
mod htscodecs;

const RANGE_VECTORS: &str =
    concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/hts-specs/cram/codecs/range/");

fn cram_arith_decode(c: &mut Criterion) {
    use criterion::Throughput;

    let mut group = c.benchmark_group("cram_arith_decode");
    // Order-0/1, RLE, PACK+RLE and STRIPE, over binned, unbinned and
    // long-read qualities and u32s.
    for name in ["q4.1", "q8.193", "q40+dir.0", "q40+dir.1", "q40+dir.65", "qvar.1", "u32.9"] {
        let src = std::fs::read(format!("{RANGE_VECTORS}{name}")).unwrap();
        // Both decoders must agree before either is timed (the integration
        // tests check them against the originals).
        let ours = seqair::cram::arith::decode(&src, 0).unwrap();
        let theirs = htscodecs::arith_uncompress(&src).expect("htscodecs rejects a vector");
        assert!(ours == theirs.as_slice(), "{name}: seqair and htscodecs decode differently");
        group.throughput(Throughput::Bytes(ours.len() as u64));

        group.bench_function(format!("seqair/{name}"), |b| {
            b.iter(|| seqair::cram::arith::decode(black_box(&src), 0).unwrap());
        });
        group.bench_function(format!("htscodecs/{name}"), |b| {
            b.iter(|| htscodecs::arith_uncompress(black_box(&src)).unwrap());
        });
    }
    group.finish();
}

// ---------------------------------------------------------------------------
// Group 5: the name tokeniser (method 8)
// The name blocks of the hts-specs CRAM 3.1 conformance files (20000 names:
// tok3 over rANS Nx16 at levels 2 and 3, over the arithmetic coder at level
// 4) and a few hts-specs tok3 vectors (1000 names each; `.1`/`.9` over rANS
// Nx16, `.11`/`.19` over the arithmetic coder), decoded by seqair and by
// htscodecs' `tok3_decode_names`; throughput is decoded bytes.
// ---------------------------------------------------------------------------

#[path = "../tests/support/cram_raw_blocks.rs"]
#[allow(dead_code, reason = "only the tok3 blocks are timed")]
mod cram_raw_blocks;

const TOK3_VECTORS: &str =
    concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/hts-specs/cram/codecs/tok3/");
const LEVELS: &str =
    concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/hts-specs/cram/3.1/passed/");

fn cram_tok3_decode(c: &mut Criterion) {
    use criterion::Throughput;

    let mut inputs: Vec<(String, Vec<Vec<u8>>)> = (2..=4)
        .map(|level| {
            let file = std::fs::read(format!("{LEVELS}level-{level}.cram")).unwrap();
            let blocks = cram_raw_blocks::raw_blocks(&file)
                .into_iter()
                .filter(|b| b.method == cram_raw_blocks::TOK3)
                .map(|b| b.data)
                .collect();
            (format!("level-{level}"), blocks)
        })
        .collect();
    for name in
        ["01.names.1", "01.names.9", "01.names.11", "01.names.19", "nv.names.9", "nv.names.19"]
    {
        inputs
            .push((name.to_owned(), vec![std::fs::read(format!("{TOK3_VECTORS}{name}")).unwrap()]));
    }

    let mut group = c.benchmark_group("cram_tok3_decode");
    for (name, blocks) in &inputs {
        // Both decoders must agree before either is timed (the tests check
        // them against the originals and against htslib).
        let mut decoded = 0;
        for block in blocks {
            let ours = seqair::cram::tok3::decode(block).unwrap();
            let theirs = htscodecs::tok3_decode_names(block).expect("htscodecs rejects a block");
            assert!(ours == theirs.as_slice(), "{name}: seqair and htscodecs decode differently");
            decoded += ours.len();
        }
        group.throughput(Throughput::Bytes(decoded as u64));

        group.bench_function(format!("seqair/{name}"), |b| {
            b.iter(|| {
                for block in blocks {
                    black_box(seqair::cram::tok3::decode(black_box(block)).unwrap());
                }
            });
        });
        group.bench_function(format!("htscodecs/{name}"), |b| {
            b.iter(|| {
                for block in blocks {
                    black_box(htscodecs::tok3_decode_names(black_box(block)).unwrap());
                }
            });
        });
    }
    group.finish();
}

criterion_group!(
    benches,
    cram_record_decode,
    cram_pileup_e2e,
    cram_full_decode,
    cram_arith_decode,
    cram_tok3_decode
);
criterion_main!(benches);
