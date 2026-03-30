//! Criterion benchmarks comparing seqair against htslib, noodles, and the bgzf crate.
//!
//! All libraries use libdeflate for BGZF decompression (fair comparison).
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::indexing_slicing)]

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use std::hint::black_box;
use std::io::Read;

const BAM_PATH: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.bam");
const FASTA_PATH: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.fasta.gz");

const CHROM: &str = "chr19";
const START: u64 = 6_105_000;
const END: u64 = 6_140_000;

fn start_pos() -> seqair::bam::Pos<seqair::bam::Zero> {
    seqair::bam::Pos::<seqair::bam::Zero>::new(START as u32)
}
fn end_pos() -> seqair::bam::Pos<seqair::bam::Zero> {
    seqair::bam::Pos::<seqair::bam::Zero>::new(END as u32)
}

// ---------------------------------------------------------------------------
// Group 1: BGZF decompression throughput
// Raw decompression comparison — all using libdeflate.
// ---------------------------------------------------------------------------

fn bgzf_decompress(c: &mut Criterion) {
    use bgzf::Reader as BgzfCrateReader;
    use noodles_bgzf as nbgzf;

    let mut group = c.benchmark_group("bgzf_decompress");

    for size in [65536usize, 1_048_576] {
        group.throughput(Throughput::Bytes(size as u64));

        // seqair: BgzfReader with read_up_to
        group.bench_with_input(BenchmarkId::new("seqair", size), &size, |b, &size| {
            b.iter(|| {
                let path = std::path::Path::new(BAM_PATH);
                let mut reader = seqair::bam::bgzf::BgzfReader::open(path).unwrap();
                let mut buf = vec![0u8; size];
                let mut total = 0usize;
                while total < size {
                    let n = reader.read_up_to(&mut buf[total..]).unwrap();
                    if n == 0 {
                        break;
                    }
                    total += n;
                }
                black_box(total)
            });
        });

        // noodles-bgzf (with libdeflate feature)
        group.bench_with_input(BenchmarkId::new("noodles", size), &size, |b, &size| {
            b.iter(|| {
                let file = std::fs::File::open(BAM_PATH).unwrap();
                let mut reader = nbgzf::io::Reader::new(file);
                let mut buf = vec![0u8; size];
                let mut total = 0usize;
                while total < size {
                    match reader.read(&mut buf[total..]) {
                        Ok(0) => break,
                        Ok(n) => total += n,
                        Err(_) => break,
                    }
                }
                black_box(total)
            });
        });

        // bgzf crate (uses libdeflate internally)
        group.bench_with_input(BenchmarkId::new("bgzf_crate", size), &size, |b, &size| {
            b.iter(|| {
                let file = std::fs::File::open(BAM_PATH).unwrap();
                let mut reader = BgzfCrateReader::new(file);
                let mut buf = vec![0u8; size];
                let mut total = 0usize;
                while total < size {
                    match reader.read(&mut buf[total..]) {
                        Ok(0) => break,
                        Ok(n) => total += n,
                        Err(_) => break,
                    }
                }
                black_box(total)
            });
        });
    }

    group.finish();
}

// ---------------------------------------------------------------------------
// Group 2: BAM record decode
// Indexed fetch + decode for a region. "htslib+base_decode" adds the cost of
// converting sequences to Vec<Base>, matching what seqair does upfront in its slab.
// ---------------------------------------------------------------------------

fn bam_record_decode(c: &mut Criterion) {
    use noodles::bam as nbam;
    use noodles::sam;
    use rust_htslib::bam::{self, FetchDefinition, Read as _};

    let mut group = c.benchmark_group("bam_record_decode");

    // seqair: IndexedBamReader + RecordStore (includes pre-decoded bases)
    group.bench_function("seqair", |b| {
        b.iter(|| {
            let path = std::path::Path::new(BAM_PATH);
            let mut reader = seqair::bam::IndexedBamReader::open(path).unwrap();
            let mut store = seqair::bam::RecordStore::new();
            let tid = reader.header().tid(CHROM).unwrap();
            reader.fetch_into(tid, start_pos(), end_pos(), &mut store).unwrap();
            black_box(store.len())
        });
    });

    // htslib: IndexedReader fetch + read loop (raw decode, no base conversion)
    group.bench_function("htslib", |b| {
        b.iter(|| {
            let mut reader = bam::IndexedReader::from_path(BAM_PATH).unwrap();
            reader.fetch(FetchDefinition::Region(0, START as i64, END as i64)).unwrap();
            let mut record = bam::Record::new();
            let mut count = 0usize;
            while let Some(Ok(())) = reader.read(&mut record) {
                count += 1;
            }
            black_box(count)
        });
    });

    // htslib + base decode: adds the per-record cost of decoding to Vec<Base>,
    // matching the work seqair does upfront in its slab during fetch_into.
    group.bench_function("htslib+base_decode", |b| {
        b.iter(|| {
            let mut reader = bam::IndexedReader::from_path(BAM_PATH).unwrap();
            reader.fetch(FetchDefinition::Region(0, START as i64, END as i64)).unwrap();
            let mut record = bam::Record::new();
            let mut total_bases = 0usize;
            while let Some(Ok(())) = reader.read(&mut record) {
                let bases: Vec<seqair_types::Base> = (0..record.seq_len())
                    .map(|i| seqair_types::Base::from(record.seq()[i]))
                    .collect();
                total_bases += bases.len();
            }
            black_box(total_bases)
        });
    });

    // noodles: sequential read + manual filter (no indexed reader without noodles-csi)
    group.bench_function("noodles", |b| {
        b.iter(|| {
            let file = std::fs::File::open(BAM_PATH).unwrap();
            let mut reader = nbam::io::Reader::new(file);
            let header: sam::Header = reader.read_header().unwrap();
            let contig_idx = header.reference_sequences().get_index_of(CHROM.as_bytes()).unwrap();

            let mut count = 0usize;
            for result in reader.records() {
                let record = result.unwrap();
                if record.flags().is_unmapped() {
                    continue;
                }
                let ref_seq_id = match record.reference_sequence_id().and_then(|r| r.ok()) {
                    Some(id) => id,
                    None => continue,
                };
                if ref_seq_id != contig_idx {
                    continue;
                }
                let aln_start = match record.alignment_start().and_then(|r| r.ok()) {
                    Some(pos) => usize::from(pos) as i64 - 1,
                    None => continue,
                };
                if aln_start >= START as i64 && aln_start < END as i64 {
                    count += 1;
                }
            }
            black_box(count)
        });
    });

    group.finish();
}

// ---------------------------------------------------------------------------
// Group 3: FASTA fetch (100KB region)
// ---------------------------------------------------------------------------

fn fasta_fetch(c: &mut Criterion) {
    use noodles::core::{Position, Region};
    use noodles::fasta as nfasta;
    use rust_htslib::faidx;

    let mut group = c.benchmark_group("fasta_fetch");

    const FASTA_START: u64 = 6_100_000;
    const FASTA_END: u64 = 6_200_000;

    group.bench_function("seqair", |b| {
        b.iter(|| {
            let path = std::path::Path::new(FASTA_PATH);
            let mut reader = seqair::fasta::IndexedFastaReader::open(path).unwrap();
            let seq = reader.fetch_seq(CHROM, FASTA_START, FASTA_END).unwrap();
            black_box(seq.len())
        });
    });

    group.bench_function("htslib", |b| {
        b.iter(|| {
            let reader = faidx::Reader::from_path(FASTA_PATH).unwrap();
            let seq =
                reader.fetch_seq(CHROM, FASTA_START as usize, FASTA_END as usize - 1).unwrap();
            black_box(seq.len())
        });
    });

    group.bench_function("noodles", |b| {
        b.iter(|| {
            let path = std::path::Path::new(FASTA_PATH);
            let mut reader =
                nfasta::io::indexed_reader::Builder::default().build_from_path(path).unwrap();
            let region = Region::new(
                CHROM,
                Position::try_from(FASTA_START as usize + 1).unwrap()
                    ..=Position::try_from(FASTA_END as usize).unwrap(),
            );
            let record = reader.query(&region).unwrap();
            black_box(record.sequence().len())
        });
    });

    group.finish();
}

// ---------------------------------------------------------------------------
// Group 4: End-to-end pileup
// The whole pipeline: open → fetch → decode → pileup iteration.
// Shows how all design choices compound.
// ---------------------------------------------------------------------------

fn pileup_e2e(c: &mut Criterion) {
    use noodles::bam as nbam;
    use noodles::sam;
    use noodles_util::alignment::iter::Depth;
    use rust_htslib::bam::{self, FetchDefinition, Read as _};

    let mut group = c.benchmark_group("pileup_e2e");

    // seqair: fetch_into → PileupEngine → iterate all columns
    group.bench_function("seqair", |b| {
        b.iter(|| {
            let path = std::path::Path::new(BAM_PATH);
            let mut reader = seqair::bam::IndexedBamReader::open(path).unwrap();
            let mut store = seqair::bam::RecordStore::new();
            let tid = reader.header().tid(CHROM).unwrap();
            reader.fetch_into(tid, start_pos(), end_pos(), &mut store).unwrap();

            let engine = seqair::bam::PileupEngine::new(store, start_pos(), end_pos());
            let mut total_depth: u64 = 0;
            let mut columns: u64 = 0;
            for col in engine {
                total_depth += col.depth() as u64;
                columns += 1;
            }
            black_box((columns, total_depth))
        });
    });

    // htslib: IndexedReader → pileup() → iterate columns
    group.bench_function("htslib", |b| {
        b.iter(|| {
            let mut reader = bam::IndexedReader::from_path(BAM_PATH).unwrap();
            let tid = reader.header().tid(CHROM.as_bytes()).unwrap();
            reader.fetch(FetchDefinition::Region(tid as i32, START as i64, END as i64)).unwrap();

            let mut total_depth: u64 = 0;
            let mut columns: u64 = 0;
            for p in reader.pileup() {
                let p = p.unwrap();
                let pos = p.pos() as u64;
                if !(START..=END).contains(&pos) {
                    continue;
                }
                total_depth += p.depth() as u64;
                columns += 1;
            }
            black_box((columns, total_depth))
        });
    });

    // noodles: sequential BAM read → noodles-util Pileup iterator
    group.bench_function("noodles", |b| {
        b.iter(|| {
            let file = std::fs::File::open(BAM_PATH).unwrap();
            let mut reader = nbam::io::Reader::new(file);
            let header: sam::Header = reader.read_header().unwrap();
            let contig_idx = header.reference_sequences().get_index_of(CHROM.as_bytes()).unwrap();

            // Collect records in the region as Box<dyn Record>
            let records: Vec<std::io::Result<Box<dyn sam::alignment::Record>>> = reader
                .records()
                .filter_map(|result| {
                    let record = result.ok()?;
                    if record.flags().is_unmapped() {
                        return None;
                    }
                    let ref_seq_id = record.reference_sequence_id().and_then(|r| r.ok())?;
                    if ref_seq_id != contig_idx {
                        return None;
                    }
                    let aln_start = record.alignment_start().and_then(|r| r.ok())?;
                    let pos = usize::from(aln_start) as u64 - 1;
                    if (START..END).contains(&pos) {
                        Some(Ok(Box::new(record) as Box<dyn sam::alignment::Record>))
                    } else {
                        None
                    }
                })
                .collect();

            let pileup = Depth::new(&header, records.into_iter());
            let mut total_depth: u64 = 0;
            let mut columns: u64 = 0;
            for result in pileup {
                let (_pos, depth): (noodles::core::Position, u64) = result.unwrap();
                total_depth += depth;
                columns += 1;
            }
            black_box((columns, total_depth))
        });
    });

    group.finish();
}

criterion_group!(benches, bgzf_decompress, bam_record_decode, fasta_fetch, pileup_e2e);
criterion_main!(benches);
