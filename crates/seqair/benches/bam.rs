//! Criterion benchmarks for BAM/BGZF reading, decoding, writing, and pileup.
//!
//! Compares seqair against htslib, noodles, and the bgzf crate.
//! All libraries use libdeflate for BGZF decompression (fair comparison).
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::indexing_slicing, reason = "benches")]
#![allow(clippy::arithmetic_side_effects, reason = "benches")]
#![allow(clippy::cast_possible_truncation, reason = "benches")]
#![allow(clippy::cast_possible_wrap, reason = "benches")]

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use seqair_types::{Base, Pos0};
use std::hint::black_box;
use std::io::Read;

const BAM_PATH: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.bam");

const CHROM: &str = "chr19";
const START: Pos0 = Pos0::new(0).unwrap();
const END: Pos0 = Pos0::new(6_140_000).unwrap();

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
            reader.fetch_into(tid, START, END, &mut store).unwrap();
            black_box(store.len())
        });
    });

    // htslib: IndexedReader fetch + read loop (raw decode, no base conversion)
    group.bench_function("htslib", |b| {
        b.iter(|| {
            let mut reader = bam::IndexedReader::from_path(BAM_PATH).unwrap();
            reader.fetch(FetchDefinition::Region(0, START.as_i64(), END.as_i64())).unwrap();
            let mut record = bam::Record::new();
            let mut count = 0usize;
            while reader.read(&mut record) == Some(Ok(())) {
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
            reader.fetch(FetchDefinition::Region(0, START.as_i64(), END.as_i64())).unwrap();
            let mut record = bam::Record::new();
            let mut total_bases = 0usize;
            while reader.read(&mut record) == Some(Ok(())) {
                let bases: Vec<Base> =
                    (0..record.seq_len()).map(|i| Base::from(record.seq()[i])).collect();
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
                let Some(ref_seq_id) = record.reference_sequence_id().and_then(|r| r.ok()) else {
                    continue;
                };
                if ref_seq_id != contig_idx {
                    continue;
                }
                let aln_start = match record.alignment_start().and_then(|r| r.ok()) {
                    Some(pos) => usize::from(pos) as i64 - 1,
                    None => continue,
                };
                if aln_start >= START.as_i64() && aln_start < END.as_i64() {
                    count += 1;
                }
            }
            black_box(count)
        });
    });

    group.finish();
}

// ---------------------------------------------------------------------------
// Group 2b: BAM read → write roundtrip
// Read all records from a region, write them to a new BAM in memory.
// Measures the full read-modify-write pipeline cost.
// ---------------------------------------------------------------------------

fn bam_roundtrip(c: &mut Criterion) {
    use noodles::bam as nbam;
    use noodles::sam;
    use noodles::sam::alignment::io::Write as _;
    use rust_htslib::bam::{self, FetchDefinition, Read as _};

    let mut group = c.benchmark_group("bam_roundtrip");

    // --- seqair: IndexedBamReader → OwnedBamRecord → BamWriter ---
    group.bench_function("seqair", |b| {
        b.iter(|| {
            use seqair::bam::BamHeader;
            use seqair::bam::owned_record::OwnedBamRecord;
            use seqair::bam::writer::BamWriter;
            use seqair::io::BgzfWriter;

            let path = std::path::Path::new(BAM_PATH);
            let mut reader = seqair::bam::IndexedBamReader::open(path).unwrap();
            let header = BamHeader::from_template(reader.header());
            let tid = reader.header().tid(CHROM).unwrap();
            let mut store = seqair::bam::RecordStore::new();
            reader.fetch_into(tid, START, END, &mut store).unwrap();

            let mut output = Vec::with_capacity(2_000_000);
            let bgzf = BgzfWriter::new(&mut output);
            let mut writer = BamWriter::new_inner(bgzf, &header, false).unwrap();

            for i in 0..store.len() {
                let idx = i as u32;
                let slim = store.record(idx);
                let cigar_bytes = store.cigar(idx);
                let mut cigar = Vec::with_capacity(slim.n_cigar_ops as usize);
                for j in 0..slim.n_cigar_ops as usize {
                    let off = j * 4;
                    let packed = u32::from_le_bytes([
                        cigar_bytes[off],
                        cigar_bytes[off + 1],
                        cigar_bytes[off + 2],
                        cigar_bytes[off + 3],
                    ]);
                    cigar.push(seqair::bam::CigarOp::from_bam_u32(packed).unwrap());
                }

                let rec = OwnedBamRecord {
                    ref_id: tid as i32,
                    pos: slim.pos.as_i64(),
                    mapq: slim.mapq,
                    flags: slim.flags,
                    next_ref_id: -1,
                    next_pos: -1,
                    template_len: 0,
                    qname: store.qname(idx).to_vec(),
                    cigar,
                    seq: store.seq(idx).to_vec(),
                    qual: store.qual(idx).to_vec(),
                    aux: seqair::bam::AuxData::from_bytes(store.aux(idx).to_vec()),
                };
                writer.write(&rec).unwrap();
            }
            writer.finish().unwrap();
            black_box(output.len())
        });
    });

    // --- htslib: IndexedReader → Record (reused) → Writer ---
    group.bench_function("htslib", |b| {
        b.iter(|| {
            let mut reader = bam::IndexedReader::from_path(BAM_PATH).unwrap();
            reader.fetch(FetchDefinition::Region(0, START.as_i64(), END.as_i64())).unwrap();
            let hdr = bam::Header::from_template(reader.header());

            // Write to temp file (htslib Writer requires a path or fd, not a Vec)
            let tmp = tempfile::NamedTempFile::new().unwrap();
            let mut writer = bam::Writer::from_path(tmp.path(), &hdr, bam::Format::Bam).unwrap();

            let mut record = bam::Record::new();
            while reader.read(&mut record) == Some(Ok(())) {
                writer.write(&record).unwrap();
            }
            drop(writer);
            black_box(std::fs::metadata(tmp.path()).unwrap().len())
        });
    });

    // --- noodles: Reader → Record → Writer ---
    group.bench_function("noodles", |b| {
        b.iter(|| {
            let file = std::fs::File::open(BAM_PATH).unwrap();
            let mut reader = nbam::io::Reader::new(file);
            let header: sam::Header = reader.read_header().unwrap();
            let contig_idx = header.reference_sequences().get_index_of(CHROM.as_bytes()).unwrap();

            let mut records = Vec::new();
            for result in reader.records() {
                let record = result.unwrap();
                if record.flags().is_unmapped() {
                    continue;
                }
                let Some(ref_seq_id) = record.reference_sequence_id().and_then(|r| r.ok()) else {
                    continue;
                };
                if ref_seq_id != contig_idx {
                    continue;
                }
                let aln_start = match record.alignment_start().and_then(|r| r.ok()) {
                    Some(pos) => usize::from(pos) as i64 - 1,
                    None => continue,
                };
                if aln_start >= START.as_i64() && aln_start < END.as_i64() {
                    records.push(record);
                }
            }

            let mut output = Vec::with_capacity(2_000_000);
            {
                let mut writer = nbam::io::Writer::new(&mut output);
                writer.write_header(&header).unwrap();
                for record in &records {
                    writer.write_alignment_record(&header, record).unwrap();
                }
            }
            black_box(output.len())
        });
    });

    group.finish();
}

// ---------------------------------------------------------------------------
// Group 4: End-to-end pileup
// The whole pipeline: open → fetch → decode → pileup iteration.
// Shows how all design choices compound.
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

/// Construct pileup, count depth and all bases
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

    // htslib: IndexedReader → pileup() → iterate columns
    group.bench_function("htslib", |b| {
        b.iter(|| {
            let mut reader = bam::IndexedReader::from_path(BAM_PATH).unwrap();
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
                    if let Some(b) = aln.record_view().seq().get(p.pos() as usize).map(Base::from) {
                        counter.count(b)
                    }
                });
            }
            black_box((columns, total_depth, counter))
        });
    });

    // noodles: sequential BAM read → noodles-util Pileup iterator
    group.bench_function("noodles-no-counter", |b| {
        b.iter(|| {
            let file = std::fs::File::open(BAM_PATH).unwrap();
            let mut reader = nbam::io::Reader::new(file);
            let header: sam::Header = reader.read_header().unwrap();
            let contig_idx = header.reference_sequences().get_index_of(CHROM.as_bytes()).unwrap();

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
                    if (START.as_u64()..END.as_u64()).contains(&pos) {
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
                let Ok((_pos, depth)) = result else {
                    continue;
                };
                total_depth += depth;
                columns += 1;
            }
            black_box((columns, total_depth))
        });
    });

    group.finish();
}

criterion_group!(benches, bgzf_decompress, bam_record_decode, bam_roundtrip, pileup_e2e);
criterion_main!(benches);
