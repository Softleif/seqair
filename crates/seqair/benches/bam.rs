//! Criterion benchmarks for BAM/BGZF reading, decoding, writing, and pileup.
//!
//! Compares seqair against htslib, noodles, and the bgzf crate.
//! All libraries use libdeflate for BGZF decompression (fair comparison).
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::indexing_slicing, reason = "benches")]
#![allow(clippy::arithmetic_side_effects, reason = "benches")]
#![allow(clippy::cast_possible_truncation, reason = "benches")]
#![allow(clippy::cast_possible_wrap, reason = "benches")]

use core::range::RangeInclusive;
use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use seqair_types::{Base, Pos0};
use std::hint::black_box;
use std::io::Read;

#[path = "support/data.rs"]
#[allow(dead_code, reason = "each bench target uses a subset")]
mod data;

const BAM_PATH: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.bam");

const CHROM: &str = "chr19";
const START: Pos0 = Pos0::new(0).unwrap();
const END: Pos0 = Pos0::new(6_140_000).unwrap();
/// The closed region every `fetch_into` below asks for.
const SPAN: RangeInclusive<Pos0> = RangeInclusive { start: START, last: END };

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
            reader.fetch_into(tid, SPAN, &mut store).unwrap();
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
            use seqair::bam::writer::BamWriterBuilder;

            let path = std::path::Path::new(BAM_PATH);
            let mut reader = seqair::bam::IndexedBamReader::open(path).unwrap();
            let header = BamHeader::from_template(reader.header());
            let tid = reader.header().tid(CHROM).unwrap();
            let mut store = seqair::bam::RecordStore::new();
            reader.fetch_into(tid, SPAN, &mut store).unwrap();

            let mut output = Vec::with_capacity(2_000_000);
            let mut writer = BamWriterBuilder::to_writer(&mut output, &header).build().unwrap();

            for idx in store.indices() {
                let slim = store.record(idx).unwrap();
                let cigar = store.record(idx).unwrap().cigar().to_vec();

                let rec = OwnedBamRecord {
                    ref_id: tid as i32,
                    pos: Some(slim.pos),
                    mapq: slim.mapq,
                    flags: slim.flags,
                    next_ref_id: -1,
                    next_pos: None,
                    template_len: 0,
                    qname: store.record(idx).unwrap().qname().to_vec(),
                    cigar,
                    seq: store.record(idx).unwrap().seq().to_vec(),
                    qual: store.record(idx).unwrap().qual().to_vec(),
                    aux: seqair::bam::AuxData::from_bytes(
                        store.record(idx).unwrap().aux().to_vec(),
                    ),
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
/// The slab-heavy shape: qname + seq + aux for every alignment of every
/// column, which is what a caller doing per-read work actually does — and what
/// `pileup_e2e` deliberately does not, since it only reads the column entry.
///
/// The two variants differ only in how often the record is resolved: once per
/// accessor, or once per alignment through `AlignmentView::record`. That
/// difference is the whole point of the handle, so it is measured rather than
/// asserted.
fn pileup_slab_access(c: &mut Criterion) {
    let mut group = c.benchmark_group("pileup_slab_access");

    let load = || {
        let path = std::path::Path::new(BAM_PATH);
        let mut reader = seqair::bam::IndexedBamReader::open(path).unwrap();
        let mut store = seqair::bam::RecordStore::new();
        let tid = reader.header().tid(CHROM).unwrap();
        reader.fetch_into(tid, SPAN, &mut store).unwrap();
        store
    };

    group.bench_function("per_accessor", |b| {
        b.iter(|| {
            let store = load();
            let mut engine = seqair::bam::PileupEngine::new(store.prepare_for_pileup().input, SPAN);
            let mut bytes: u64 = 0;
            while let Some(col) = engine.pileups() {
                for aln in col.alignments() {
                    bytes +=
                        aln.qname().len() as u64 + aln.seq().len() as u64 + aln.aux().len() as u64;
                }
            }
            black_box(bytes)
        });
    });

    group.bench_function("via_record", |b| {
        b.iter(|| {
            let store = load();
            let mut engine = seqair::bam::PileupEngine::new(store.prepare_for_pileup().input, SPAN);
            let mut bytes: u64 = 0;
            while let Some(col) = engine.pileups() {
                for aln in col.alignments() {
                    let rec = aln.record();
                    bytes +=
                        rec.qname().len() as u64 + rec.seq().len() as u64 + rec.aux().len() as u64;
                }
            }
            black_box(bytes)
        });
    });

    group.finish();
}

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
                    let Some(qpos) = aln.qpos() else { return };
                    let record = aln.record();
                    counter.count(Base::from(record.seq()[qpos]))
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

// ---------------------------------------------------------------------------
// Group 5: AlignedPairs walking (per-record CIGAR walk)
//
// Compares seqair's AlignedPairs (and the layered with_read variant) against
// rust-htslib's aligned_pairs_full (per-base) extension. Both walk every
// record in a region and visit every (qpos, rpos) pair.
//
// Records are pre-loaded once outside the timer so we measure only the walk,
// not the fetch overhead. (The fetch path is benchmarked separately in
// `bam_record_decode`.)
// ---------------------------------------------------------------------------

fn aligned_pairs_walk(c: &mut Criterion) {
    use rust_htslib::bam::{self, FetchDefinition, Read as _, ext::BamRecordExtensions};

    let mut group = c.benchmark_group("aligned_pairs_walk");

    // ── seqair: pre-load the RecordStore once ──
    let seqair_store = {
        let path = std::path::Path::new(BAM_PATH);
        let mut reader = seqair::bam::IndexedBamReader::open(path).unwrap();
        let mut store = seqair::bam::RecordStore::new();
        let tid = reader.header().tid(CHROM).unwrap();
        reader.fetch_into(tid, SPAN, &mut store).unwrap();
        store
    };

    // ── htslib: pre-load all records into a Vec once ──
    let htslib_records: Vec<bam::Record> = {
        let mut reader = bam::IndexedReader::from_path(BAM_PATH).unwrap();
        let tid = reader.header().tid(CHROM.as_bytes()).unwrap();
        reader.fetch(FetchDefinition::Region(tid as i32, START.as_i64(), END.as_i64())).unwrap();
        let mut out = Vec::new();
        let mut rec = bam::Record::new();
        while reader.read(&mut rec) == Some(Ok(())) {
            out.push(rec.clone());
        }
        out
    };

    // seqair bare AlignedPairs: positions only (no seq/qual lookup).
    group.bench_function("seqair_bare", |b| {
        b.iter(|| {
            let mut total_pairs: u64 = 0;
            for idx in seqair_store.indices() {
                let rec = seqair_store.record(idx).unwrap();
                for pair in rec.aligned_pairs(&seqair_store).unwrap() {
                    let _ = black_box(pair);
                    total_pairs += 1;
                }
            }
            black_box(total_pairs)
        });
    });

    // seqair with_read: pairs + per-position read base/qual.
    group.bench_function("seqair_with_read", |b| {
        b.iter(|| {
            let mut total_pairs: u64 = 0;
            for idx in seqair_store.indices() {
                let rec = seqair_store.record(idx).unwrap();
                for ev in rec.aligned_pairs_with_read(&seqair_store).unwrap() {
                    let _ = black_box(ev);
                    total_pairs += 1;
                }
            }
            black_box(total_pairs)
        });
    });

    // seqair matches_only: positions filtered to Match-only.
    group.bench_function("seqair_matches_only", |b| {
        b.iter(|| {
            let mut total: u64 = 0;
            for idx in seqair_store.indices() {
                let rec = seqair_store.record(idx).unwrap();
                for m in rec.aligned_pairs(&seqair_store).unwrap().matches_only() {
                    let _ = black_box(m);
                    total += 1;
                }
            }
            black_box(total)
        });
    });

    // rust-htslib aligned_pairs_full (the closest analogue — per-base events).
    group.bench_function("htslib_aligned_pairs_full", |b| {
        b.iter(|| {
            let mut total: u64 = 0;
            for record in &htslib_records {
                for pair in record.aligned_pairs_full() {
                    let _ = black_box(pair);
                    total += 1;
                }
            }
            black_box(total)
        });
    });

    // rust-htslib aligned_pairs (matches-only — filtered to M positions).
    group.bench_function("htslib_aligned_pairs", |b| {
        b.iter(|| {
            let mut total: u64 = 0;
            for record in &htslib_records {
                for pair in record.aligned_pairs() {
                    let _ = black_box(pair);
                    total += 1;
                }
            }
            black_box(total)
        });
    });

    group.finish();
}

// ---------------------------------------------------------------------------
// Group 5: End-to-end pileup *with reference fetch*
// `pileup_e2e` feeds the engine no reference at all, so the FASTA/FASTQ fetch
// path (newline-strip + uppercase) never shows up in it. This group times the
// full `Readers::pileup` loop — per-segment reference fetch included — once
// against a BGZF FASTA and once against a FASTQ view of the same reference.
// ---------------------------------------------------------------------------

fn pileup_with_reference(c: &mut Criterion) {
    use rust_htslib::bam::{self, FetchDefinition, Read as _};
    use rust_htslib::faidx;
    use seqair::reader::{DepthLimit, Readers, SegmentOptions};
    use std::num::NonZeroU32;

    let bam = std::path::Path::new(BAM_PATH);
    let fasta = data::plain_fasta();
    let fastq = data::chr19_fastq();

    // Populate reads: chr19 reads span ~6.10–6.14 Mb.
    const P_START: Pos0 = Pos0::new(6_100_000).unwrap();
    const P_END: Pos0 = Pos0::new(6_140_000).unwrap();

    let mut group = c.benchmark_group("pileup_with_reference");
    group.sample_size(30);

    let opts = SegmentOptions::new(NonZeroU32::new(100_000).unwrap());

    // seqair: Readers::pileup fetches the reference slice per segment, then
    // piles up reads against it. The closure keeps the borrow story simple:
    // the guard borrows `readers`, so it must be drained and dropped before
    // the next `segments()` call.
    let bench_seqair = |ref_path: &std::path::Path| {
        let mut readers = Readers::open(bam, ref_path).unwrap();
        let segments: Vec<_> =
            readers.segments((CHROM, (P_START..=P_END).into()), opts).unwrap().collect();
        let mut total_depth: u64 = 0;
        let mut mismatches: u64 = 0;
        for seg in &segments {
            let mut p = readers.pileup(seg, DepthLimit::Unlimited).run().unwrap();
            while let Some(col) = p.pileups() {
                total_depth += col.depth() as u64;
                let ref_base = col.reference_base();
                for aln in col.alignments() {
                    if let Some(base) = aln.base()
                        && base != ref_base
                    {
                        mismatches += 1;
                    }
                }
            }
        }
        (total_depth, mismatches)
    };

    // Cross-check once: FASTA and FASTQ views must agree, and both must see
    // real reads (otherwise the pileup loop silently benches nothing).
    let (expected_depth, expected_mm) = bench_seqair(&fasta);
    assert!(expected_depth > 0, "test BAM should produce pileup depth in this region");
    assert_eq!(
        bench_seqair(&fastq),
        (expected_depth, expected_mm),
        "FASTA and FASTQ references disagree"
    );

    group.bench_function("seqair_fasta", |b| {
        b.iter(|| black_box(bench_seqair(&fasta)));
    });
    group.bench_function("seqair_fastq", |b| {
        b.iter(|| black_box(bench_seqair(&fastq)));
    });

    // htslib: reads pileup + faidx slice per 100 kb tile, uppercased to match
    // seqair's `r[fasta.fetch.uppercase]`. Tiles are cached: a naive
    // fetch-per-column runs ~340x slower than seqair (3.1 s vs 9 ms), which
    // measures faidx syscall overhead, not anything a real caller would do.
    group.bench_function("htslib_fasta", |b| {
        b.iter(|| {
            let mut reader = bam::IndexedReader::from_path(bam).unwrap();
            let ref_reader = faidx::Reader::from_path(&fasta).unwrap();
            let tid = reader.header().tid(CHROM.as_bytes()).unwrap();
            reader
                .fetch(FetchDefinition::Region(tid as i32, P_START.as_i64(), P_END.as_i64()))
                .unwrap();

            let mut total_depth: u64 = 0;
            let mut mismatches: u64 = 0;
            let mut tile = u64::MAX; // sentinel: no tile fetched yet
            let mut ref_seq: Vec<u8> = Vec::new();
            let mut tile_start = 0u64;
            for p in reader.pileup() {
                let p = p.unwrap();
                let pos = u64::from(p.pos());
                if !(P_START.as_u64()..=P_END.as_u64()).contains(&pos) {
                    continue;
                }
                total_depth += u64::from(p.depth());
                // The tile index mirrors seqair's 100 kb segments; fetch_seq
                // is inclusive-end, hence the -1.
                let t = pos / 100_000;
                if t != tile {
                    tile = t;
                    tile_start = t * 100_000;
                    let e = ((t + 1) * 100_000).min(P_END.as_u64() + 1);
                    ref_seq =
                        ref_reader.fetch_seq(CHROM, tile_start as usize, e as usize - 1).unwrap();
                    ref_seq.make_ascii_uppercase();
                }
                let ref_base = ref_seq[(pos - tile_start) as usize];
                for aln in p.alignments() {
                    let Some(qpos) = aln.qpos() else { continue };
                    let record = aln.record();
                    if record.seq()[qpos] != ref_base {
                        mismatches += 1;
                    }
                }
            }
            black_box((total_depth, mismatches))
        });
    });

    group.finish();
}

// ---------------------------------------------------------------------------
// Group 7: tiled region queries — the shape a variant caller actually uses
//
// The groups above issue one query for a whole span. rastair does the
// opposite: it tiles a contig into 10 kb segments and issues an independent
// index query per segment, so everything that is paid *per query* — the BAI
// lookup, `RegionBuf` setup, re-decompressing the block a tile boundary falls
// in, sorting and mate-linking the store — is paid ~13 k times on a
// chromosome instead of once.
//
// Varying the tile size over one fixed span turns that overhead into a slope:
// the same reads, columns and bases come out of every tile size, so anything
// that grows as tiles shrink is per-query cost and nothing else.
// ---------------------------------------------------------------------------

fn pileup_tiled(c: &mut Criterion) {
    use seqair::reader::{DepthLimit, Readers, SegmentOptions};
    use std::num::NonZeroU32;

    let bam = std::path::Path::new(BAM_PATH);
    let fasta = data::plain_fasta();

    // The span the test BAM's reads actually cover.
    const T_START: Pos0 = Pos0::new(6_100_000).unwrap();
    const T_END: Pos0 = Pos0::new(6_140_000).unwrap();
    /// rastair's default `--segment-overlap`.
    const OVERLAP: u32 = 200;
    /// rastair's default `--max-coverage`, so the `DepthCap` wrapper the
    /// per-column limit installs is in the measured path.
    const DEPTH_CAP: u32 = 8_000;

    let mut group = c.benchmark_group("pileup_tiled");
    group.sample_size(20);

    let run = |readers: &mut Readers, opts: SegmentOptions| {
        let segments: Vec<_> =
            readers.segments((CHROM, (T_START..=T_END).into()), opts).unwrap().collect();
        let cap = DepthLimit::PerColumn(NonZeroU32::new(DEPTH_CAP).unwrap());
        let mut columns: u64 = 0;
        let mut bases: u64 = 0;
        for seg in &segments {
            let mut p = readers.pileup(seg, cap).run().unwrap();
            while let Some(col) = p.pileups() {
                // Tiles overlap, so count each column once — the totals have to
                // be identical across tile sizes for the comparison to mean
                // anything.
                if !seg.core_span().contains(&col.pos()) {
                    continue;
                }
                columns += 1;
                for aln in col.alignments() {
                    bases += u64::from(aln.base().is_some());
                }
            }
        }
        (columns, bases)
    };

    // One tile over the whole span is the floor: the same work with the
    // per-query overhead paid once.
    let mut readers = Readers::open(bam, &fasta).unwrap();
    let whole = SegmentOptions::new(NonZeroU32::new(u32::MAX).unwrap());
    let expected = run(&mut readers, whole);
    assert!(expected.0 > 0, "test BAM should produce columns in this region");

    for tile in [1_000u32, 10_000, 100_000] {
        let opts =
            SegmentOptions::new(NonZeroU32::new(tile).unwrap()).with_overlap(OVERLAP).unwrap();
        assert_eq!(
            run(&mut readers, opts),
            expected,
            "tile size {tile} changed what the pileup sees"
        );
        group.bench_with_input(BenchmarkId::new("seqair", tile), &opts, |b, &opts| {
            b.iter(|| black_box(run(&mut readers, opts)));
        });
    }

    group.bench_function("seqair_single_query", |b| {
        b.iter(|| black_box(run(&mut readers, whole)));
    });

    group.finish();
}

// ---------------------------------------------------------------------------
// Group: the pileup engine alone, on synthetic reads
// No I/O in the timed loop: the store is built in the batch setup, so this
// isolates column construction. `spliced` is RNA-seq shaped — every read is
// split by one or two `N` introns — which is where a stateless per-column
// CIGAR lookup costs the most.
// ---------------------------------------------------------------------------

/// Deterministic xorshift, so every run piles up the same reads.
struct XorShift(u64);

impl XorShift {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn below(&mut self, n: u64) -> u32 {
        (self.next() % n) as u32
    }
}

#[derive(Clone, Copy)]
enum ReadShape {
    /// 150 bp reads; one in 25 carries a 1–3 bp insertion or deletion.
    Dna,
    /// Two or three exons of 150 aligned bases in total, split by introns of
    /// 100–5000 bp.
    Spliced,
}

fn synthetic_input(shape: ReadShape) -> seqair::bam::record_store::PileupInput<()> {
    use seqair::bam::cigar::{CigarOp, CigarOpType, calc_matches_indels, compute_end_pos};
    const N_READS: u32 = 200_000;
    const CONTIG: u32 = 1_000_000;
    let mut rng = XorShift(0x9E37_79B9_7F4A_7C15);
    let mut starts: Vec<u32> = (0..N_READS).map(|_| rng.below(u64::from(CONTIG))).collect();
    starts.sort_unstable();
    let bases = [Base::A, Base::C, Base::G, Base::T];
    let mut store = seqair::bam::RecordStore::<()>::new();
    let m = |len| CigarOp::new(CigarOpType::Match, len);
    for (i, &start) in starts.iter().enumerate() {
        let ops: Vec<CigarOp> = match shape {
            ReadShape::Dna => match rng.below(50) {
                0 => vec![m(70), CigarOp::new(CigarOpType::Insertion, 1 + rng.below(3)), m(78)],
                1 => vec![m(70), CigarOp::new(CigarOpType::Deletion, 1 + rng.below(3)), m(80)],
                _ => vec![m(150)],
            },
            ReadShape::Spliced => {
                let intron =
                    |rng: &mut XorShift| CigarOp::new(CigarOpType::RefSkip, 100 + rng.below(4900));
                if rng.below(3) == 0 {
                    vec![m(40), intron(&mut rng), m(70), intron(&mut rng), m(40)]
                } else {
                    let a = 20 + rng.below(110);
                    vec![m(a), intron(&mut rng), m(150 - a)]
                }
            }
        };
        let qlen = seqair::bam::cigar::calc_query_len(&ops) as usize;
        let seq: Vec<Base> = (0..qlen).map(|_| bases[rng.below(4) as usize]).collect();
        let qual: Vec<u8> = (0..qlen).map(|_| 20 + rng.below(20) as u8).collect();
        let pos = Pos0::new(start).unwrap();
        let (matching, indels) = calc_matches_indels(&ops);
        store
            .push_fields(
                pos,
                compute_end_pos(pos, &ops).unwrap(),
                seqair_types::BamFlags::from(if i % 2 == 0 { 0 } else { 16 }),
                60,
                matching,
                indels,
                format!("r{i}").as_bytes(),
                &ops,
                &seq,
                &qual,
                &[],
                0,
                -1,
                -1,
                0,
                &mut (),
            )
            .unwrap();
    }
    store.prepare_for_pileup().input
}

fn pileup_synthetic(c: &mut Criterion) {
    let mut group = c.benchmark_group("pileup_synthetic");
    // The capped arm is an amplicon-style deep column under `max_depth`: the
    // spliced columns hold ~1000 entries, of which it keeps 50.
    for (name, shape, cap) in [
        ("dna", ReadShape::Dna, None),
        ("spliced", ReadShape::Spliced, None),
        ("spliced_max_depth_50", ReadShape::Spliced, std::num::NonZeroU32::new(50)),
    ] {
        group.bench_function(name, |b| {
            b.iter_batched(
                || synthetic_input(shape),
                |input| {
                    let span = RangeInclusive {
                        start: Pos0::new(0).unwrap(),
                        last: Pos0::new(1_010_000).unwrap(),
                    };
                    let mut engine = seqair::bam::PileupEngine::new(input, span);
                    if let Some(cap) = cap {
                        engine.set_max_depth(cap);
                    }
                    let mut counter = Counter::new();
                    let mut mapq: u64 = 0;
                    while let Some(col) = engine.pileups() {
                        for aln in col.alignments() {
                            mapq += u64::from(aln.mapq);
                            counter.count(aln.base().unwrap_or_default());
                        }
                    }
                    black_box((counter, mapq))
                },
                criterion::BatchSize::LargeInput,
            );
        });
    }
    group.finish();
}

criterion_group!(
    benches,
    bgzf_decompress,
    bam_record_decode,
    bam_roundtrip,
    pileup_e2e,
    pileup_slab_access,
    aligned_pairs_walk,
    pileup_with_reference,
    pileup_tiled,
    pileup_synthetic,
);
criterion_main!(benches);
