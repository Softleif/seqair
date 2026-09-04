//! Criterion benchmarks for FASTQ access.
//!
//! Three independent questions are measured here.
//!
//! **Indexed range fetch** — the access pattern `Readers::open` / `pileup`
//! needs. `samtools faidx` indexes FASTQ into a 6-column FAI (the extra column
//! is `qual_offset`); there is no separate "FQI" format. Of the Rust ecosystem:
//!
//! - `noodles-fastq` can *build* and *parse* that index (`fai::Record` exposes
//!   `quality_scores_offset`) but has no `IndexedReader`, no `query` and no
//!   `seek` — its `Reader` offers only `read_record`/`records`. It therefore
//!   cannot be benchmarked as a FASTQ contender at all; the capability is
//!   absent, not slow. `noodles-fasta`'s `IndexedReader` *is* benchmarked,
//!   pointed at the same 5-column view seqair gets, because the byte-span math
//!   is format-independent.
//! - `rust-htslib` exposes only `fai_load` (FASTA), not `fai_load_format` with
//!   `FAI_FASTQ`. htslib's FASTA index parser reads the first four numeric
//!   fields with `sscanf` and ignores a trailing fifth, so pointing it at a
//!   FASTQ's 6-column FAI still works; `htslib_on_fastq` measures that.
//! - `seqair` is measured through a 5-column view of the same index (see
//!   `data::five_column_view`), which is exactly the index a 6-column-aware
//!   parser would build.
//!
//! **Sequential iteration** — throughput for a streaming reader, comparing
//! `noodles-fastq`, `needletail` and `seq_io` over plain, gzip and BGZF inputs,
//! both counting bases only and touching every field.
//!
//! **The layer underneath** — raw BGZF/gzip decompression throughput, which
//! turns out to dominate every compressed measurement above.
//!
//! # Fairness
//!
//! What each contender is and is not charged for, since the numbers are only
//! worth anything if that is explicit.
//!
//! * **Buffering is matched.** `seq_io` and `needletail` each own a 64 KiB
//!   `buffer_redux` reader, so they are handed the raw source — wrapping them
//!   in a `BufReader` would stack a second buffer in front and charge them an
//!   extra copy of every byte. `noodles` takes `BufRead` by construction, so it
//!   gets a `BufReader` at that same 64 KiB.
//! * **seqair is charged for work the others skip.** `fetch_seq*` uppercases
//!   per `r[fasta.fetch.uppercase]`; htslib and noodles return the soft-mask
//!   case verbatim. The comparison is therefore against seqair, not for it.
//! * **Allocation is shown both ways.** `seqair_fetch_into` reuses a caller
//!   buffer, which neither of the others offers; `seqair_fetch_seq` allocates
//!   per call like they do, and both are reported.
//! * **Query construction counts only where it is real.** `noodles::query`
//!   needs an allocated `Region`; seqair and htslib take `&str`. Where the
//!   region is loop-invariant it is hoisted out of the timed section; in
//!   `fastq_pileup_walk`, where it changes every step, it is not.
//! * **The gap to htslib is structural, not binding overhead.**
//!   `rust_htslib`'s `fetch_seq` adopts htslib's `malloc`'d buffer rather than
//!   copying it, so the binding adds one `CString` per call and little else.
//!   The difference is in `fai_retrieve`, which issues one `bgzf_read` *per
//!   line* of the requested span and goes through the BGZF layer even for a
//!   plain file; seqair issues a single `pread` over the whole span and strips
//!   newlines in one pass. Expect the gap to widen as `linebases` shrinks.
//! * **noodles decodes eagerly**, copying name/sequence/quality into an owned
//!   record, where `seq_io` and `needletail` hand out borrowed slices. That is
//!   a design difference rather than an inefficiency, which is why
//!   `fastq_sequential_full` touches every field of every record.
//! * **`noodles_bgzf_mt` pays thread startup inside each iteration**, because
//!   the reader consumes the file and cannot be hoisted. On a ~70 MB input that
//!   is a real part of the cost, but it is not a steady-state number.
//! * **Everything runs warm.** Criterion re-reads the same files thousands of
//!   times, so the page cache is hot for every contender. These measure CPU and
//!   syscall cost, not cold-start I/O.
//! * **Inputs avoid the uniform-random trap.** See `data::reads_fastq` — reads
//!   are sampled from a synthetic genome with realistic error and 4-bin
//!   Illumina quality, giving a ~4.3x gzip ratio rather than the ~2x that
//!   uniform-random data yields and that would distort every compressed
//!   measurement. `fastq_indexed_fetch_real` runs against chr19 re-emitted as
//!   FASTQ, with its real base composition, soft-masking and centromeric `N`
//!   runs.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "benches")]
#![allow(clippy::cast_possible_truncation, reason = "benches")]
#![allow(clippy::arithmetic_side_effects, reason = "benches")]
#![allow(clippy::print_stdout, reason = "benches report to the console")]

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use std::fs::File;
use std::hint::black_box;
use std::io::{BufReader, Read};
use std::path::Path;

#[path = "support/data.rs"]
#[allow(dead_code, reason = "each bench target uses a subset")]
mod data;

const CTG: &str = "ctg_a";
const START: u64 = 5_000_000;
const SIZES: [u64; 3] = [1_000, 10_000, 100_000];

fn pos(v: u64) -> seqair::bam::Pos<seqair::bam::Zero> {
    seqair::bam::Pos::<seqair::bam::Zero>::new(v as u32).unwrap()
}

fn region(chrom: &str, start: u64, end: u64) -> noodles::core::Region {
    use noodles::core::{Position, Region};
    Region::new(
        chrom,
        Position::try_from(start as usize + 1).unwrap()..=Position::try_from(end as usize).unwrap(),
    )
}

// ---------------------------------------------------------------------------
// Group 1: indexed range fetch on FASTQ
// ---------------------------------------------------------------------------

/// Sweep region sizes on a FASTQ carrying three long "contigs".
fn fastq_indexed_fetch(c: &mut Criterion) {
    use rust_htslib::faidx;

    let real = data::ref_fastq();
    let view = data::ref_fastq_seqair_view();

    let mut group = c.benchmark_group("fastq_indexed_fetch");

    let mut seqair_reader = seqair::fasta::IndexedFastaReader::open(&view).unwrap();
    let htslib_reader = faidx::Reader::from_path(&real);
    let mut noodles_reader =
        noodles::fasta::io::indexed_reader::Builder::default().build_from_path(&view).unwrap();
    let mut reuse = Vec::new();

    // Correctness cross-check before timing: all three paths must agree.
    // seqair uppercases per `r[fasta.fetch.uppercase]`; htslib preserves the
    // soft-mask case, so normalise before comparing.
    seqair_reader.fetch_seq_into(CTG, pos(START), pos(START + 1_000), &mut reuse).unwrap();
    if let Ok(ref h) = htslib_reader {
        let mut hs = h.fetch_seq(CTG, START as usize, (START + 1_000) as usize - 1).unwrap();
        hs.make_ascii_uppercase();
        assert_eq!(reuse, hs, "seqair and htslib disagree on FASTQ subsequence");
    }
    let nr = noodles_reader.query(&region(CTG, START, START + 1_000)).unwrap();
    assert_eq!(reuse.as_slice(), nr.sequence().as_ref(), "seqair and noodles disagree");

    for size in SIZES {
        let end = START + size;
        group.throughput(Throughput::Bytes(size));

        group.bench_with_input(BenchmarkId::new("seqair_fetch_into", size), &size, |b, _| {
            b.iter(|| {
                seqair_reader.fetch_seq_into(CTG, pos(START), pos(end), &mut reuse).unwrap();
                black_box(reuse.len())
            });
        });

        // seqair's allocating API, for a like-for-like comparison with the
        // other two (neither offers a caller-supplied buffer).
        group.bench_with_input(BenchmarkId::new("seqair_fetch_seq", size), &size, |b, _| {
            b.iter(|| black_box(seqair_reader.fetch_seq(CTG, pos(START), pos(end)).unwrap().len()));
        });

        if let Ok(ref h) = htslib_reader {
            group.bench_with_input(BenchmarkId::new("htslib_on_fastq", size), &size, |b, _| {
                b.iter(|| {
                    black_box(h.fetch_seq(CTG, START as usize, end as usize - 1).unwrap().len())
                });
            });
        }

        // The region is constant across iterations here, so building it is
        // hoisted out of the timed loop: `Region::new` allocates, and neither
        // seqair nor htslib constructs a query object at all. (In
        // `fastq_pileup_walk` the region changes every step, so there building
        // it is genuinely part of the work and stays inside.)
        let query_region = region(CTG, START, end);
        group.bench_with_input(BenchmarkId::new("noodles_fasta_on_fastq", size), &size, |b, _| {
            b.iter(|| black_box(noodles_reader.query(&query_region).unwrap().sequence().len()));
        });
    }

    group.finish();
}

/// The same fetch on real reference content — chr19 re-emitted as FASTQ,
/// 50 bases/line, including the centromeric `N` runs.
fn fastq_indexed_fetch_real(c: &mut Criterion) {
    use rust_htslib::faidx;

    let real = data::chr19_fastq();
    let view = data::five_column_view(&real, "chr19_seqair.fq");

    let mut group = c.benchmark_group("fastq_indexed_fetch_real");

    let mut seqair_reader = seqair::fasta::IndexedFastaReader::open(&view).unwrap();
    let htslib_reader = faidx::Reader::from_path(&real);
    let mut reuse = Vec::new();

    const CHR: &str = "chr19";
    const CHR_START: u64 = 6_100_000;

    // chr19 carries soft-masked (lowercase) repeat regions; seqair uppercases
    // them per `r[fasta.fetch.uppercase]` while htslib returns them verbatim.
    seqair_reader.fetch_seq_into(CHR, pos(CHR_START), pos(CHR_START + 1_000), &mut reuse).unwrap();
    if let Ok(ref h) = htslib_reader {
        let mut hs =
            h.fetch_seq(CHR, CHR_START as usize, (CHR_START + 1_000) as usize - 1).unwrap();
        hs.make_ascii_uppercase();
        assert_eq!(reuse, hs, "seqair and htslib disagree on real FASTQ subsequence");
    }

    for size in SIZES {
        let end = CHR_START + size;
        group.throughput(Throughput::Bytes(size));

        group.bench_with_input(BenchmarkId::new("seqair_fetch_into", size), &size, |b, _| {
            b.iter(|| {
                seqair_reader.fetch_seq_into(CHR, pos(CHR_START), pos(end), &mut reuse).unwrap();
                black_box(reuse.len())
            });
        });

        if let Ok(ref h) = htslib_reader {
            group.bench_with_input(BenchmarkId::new("htslib_on_fastq", size), &size, |b, _| {
                b.iter(|| {
                    black_box(h.fetch_seq(CHR, CHR_START as usize, end as usize - 1).unwrap().len())
                });
            });
        }
    }

    group.finish();
}

/// The access pattern `Readers::pileup` actually generates: many small,
/// consecutive reference slices walking a contig. Per-call overhead — seek,
/// allocation, name lookup — dominates here, not bulk copy speed.
fn fastq_pileup_walk(c: &mut Criterion) {
    use rust_htslib::faidx;

    const WINDOW: u64 = 500;
    const STEPS: u64 = 2_000;

    let real = data::ref_fastq();
    let view = data::ref_fastq_seqair_view();

    let mut group = c.benchmark_group("fastq_pileup_walk");
    group.throughput(Throughput::Bytes(WINDOW * STEPS));
    group.sample_size(30);

    let mut seqair_reader = seqair::fasta::IndexedFastaReader::open(&view).unwrap();
    let htslib_reader = faidx::Reader::from_path(&real);
    let mut noodles_reader =
        noodles::fasta::io::indexed_reader::Builder::default().build_from_path(&view).unwrap();
    let mut reuse = Vec::new();

    group.bench_function("seqair_fetch_into", |b| {
        b.iter(|| {
            let mut total = 0usize;
            for i in 0..STEPS {
                let s = i * WINDOW;
                seqair_reader.fetch_seq_into(CTG, pos(s), pos(s + WINDOW), &mut reuse).unwrap();
                total += reuse.len();
            }
            black_box(total)
        });
    });

    if let Ok(ref h) = htslib_reader {
        group.bench_function("htslib_on_fastq", |b| {
            b.iter(|| {
                let mut total = 0usize;
                for i in 0..STEPS {
                    let s = (i * WINDOW) as usize;
                    total += h.fetch_seq(CTG, s, s + WINDOW as usize - 1).unwrap().len();
                }
                black_box(total)
            });
        });
    }

    group.bench_function("noodles_fasta_on_fastq", |b| {
        b.iter(|| {
            let mut total = 0usize;
            for i in 0..STEPS {
                let s = i * WINDOW;
                total +=
                    noodles_reader.query(&region(CTG, s, s + WINDOW)).unwrap().sequence().len();
            }
            black_box(total)
        });
    });

    group.finish();
}

// ---------------------------------------------------------------------------
// Group 2: sequential iteration throughput
// ---------------------------------------------------------------------------

fn count_noodles<R: std::io::BufRead>(inner: R) -> usize {
    use noodles::fastq;
    let mut reader = fastq::io::Reader::new(inner);
    let mut record = fastq::Record::default();
    let mut n = 0;
    while reader.read_record(&mut record).unwrap() != 0 {
        n += record.sequence().len();
    }
    n
}

/// Touch every field, not just the sequence length. `noodles` decodes eagerly
/// into an owned record, so this should cost it nothing extra; `seq_io` and
/// `needletail` hand out borrowed slices, so this is where any laziness they
/// exploit has to be paid for.
fn count_noodles_full<R: std::io::BufRead>(inner: R) -> usize {
    use noodles::fastq;
    let mut reader = fastq::io::Reader::new(inner);
    let mut record = fastq::Record::default();
    let mut n = 0;
    while reader.read_record(&mut record).unwrap() != 0 {
        n += record.name().len() + record.sequence().len() + record.quality_scores().len();
    }
    n
}

fn count_seq_io<R: Read>(inner: R) -> usize {
    use seq_io::fastq::Record as _;
    let mut reader = seq_io::fastq::Reader::new(inner);
    let mut n = 0;
    while let Some(rec) = reader.next() {
        n += rec.unwrap().seq().len();
    }
    n
}

fn count_seq_io_full<R: Read>(inner: R) -> usize {
    use seq_io::fastq::Record as _;
    let mut reader = seq_io::fastq::Reader::new(inner);
    let mut n = 0;
    while let Some(rec) = reader.next() {
        let rec = rec.unwrap();
        n += rec.head().len() + rec.seq().len() + rec.qual().len();
    }
    n
}

/// `seq_io`'s batched API: fill a `RecordSet` per call instead of one record at
/// a time. This is the shape a parallel reader would use.
fn count_seq_io_recordset<R: Read>(inner: R) -> usize {
    use seq_io::fastq::Record as _;
    let mut reader = seq_io::fastq::Reader::new(inner);
    let mut rset = seq_io::fastq::RecordSet::default();
    let mut n = 0;
    while let Some(res) = reader.read_record_set(&mut rset) {
        res.unwrap();
        for rec in &rset {
            n += rec.seq().len();
        }
    }
    n
}

fn count_needletail(path: &Path) -> usize {
    let mut reader = needletail::parse_fastx_file(path).unwrap();
    let mut n = 0;
    while let Some(rec) = reader.next() {
        n += rec.unwrap().seq().len();
    }
    n
}

/// `needletail` fed the same already-decompressed stream the others get, so the
/// compressed groups compare parsers rather than each crate's plumbing. Its own
/// `parse_fastx_file` re-sniffs the magic bytes and stacks two `Chain` adapters
/// in front of the decoder, which is measured separately.
fn count_needletail_reader<R: std::io::Read + Send>(inner: R) -> usize {
    let mut reader = needletail::parse_fastx_reader(inner).unwrap();
    let mut n = 0;
    while let Some(rec) = reader.next() {
        n += rec.unwrap().seq().len();
    }
    n
}

fn count_needletail_full(path: &Path) -> usize {
    let mut reader = needletail::parse_fastx_file(path).unwrap();
    let mut n = 0;
    while let Some(rec) = reader.next() {
        let rec = rec.unwrap();
        n += rec.id().len() + rec.seq().len() + rec.qual().map_or(0, |q| q.len());
    }
    n
}

/// `noodles` requires `BufRead`, so it gets a 64 KiB `BufReader` — the same
/// capacity `seq_io` and `needletail` use for the buffers they own.
fn open_plain(path: &Path) -> BufReader<File> {
    BufReader::with_capacity(1 << 16, File::open(path).unwrap())
}

fn open_gz(path: &Path) -> BufReader<flate2::read::MultiGzDecoder<File>> {
    BufReader::with_capacity(1 << 16, flate2::read::MultiGzDecoder::new(File::open(path).unwrap()))
}

/// `seq_io` and `needletail` take `Read` and wrap it in their *own* 64 KiB
/// `buffer_redux` reader. Handing them a `BufReader` would stack a second
/// buffer in front of that and charge them an extra copy of every byte, so
/// they get the raw source.
fn open_raw(path: &Path) -> File {
    File::open(path).unwrap()
}

fn open_gz_raw(path: &Path) -> flate2::read::MultiGzDecoder<File> {
    flate2::read::MultiGzDecoder::new(File::open(path).unwrap())
}

fn fastq_sequential(c: &mut Criterion) {
    let (plain, gz, bgz) = data::reads_fastq();
    let plain_len = std::fs::metadata(&plain).unwrap().len();

    for (label, path) in
        [("plain", plain.as_path()), ("gzip", gz.as_path()), ("bgzf", bgz.as_path())]
    {
        let mut group = c.benchmark_group(format!("fastq_sequential/{label}"));
        // Throughput is reported against the uncompressed size in every case,
        // so the three inputs are directly comparable.
        group.throughput(Throughput::Bytes(plain_len));
        group.sample_size(20);
        let compressed = label != "plain";

        group.bench_function("noodles", |b| {
            b.iter(|| {
                black_box(if compressed {
                    count_noodles(open_gz(path))
                } else {
                    count_noodles(open_plain(path))
                })
            });
        });

        group.bench_function("seq_io", |b| {
            b.iter(|| {
                black_box(if compressed {
                    count_seq_io(open_gz_raw(path))
                } else {
                    count_seq_io(open_raw(path))
                })
            });
        });

        group.bench_function("seq_io_recordset", |b| {
            b.iter(|| {
                black_box(if compressed {
                    count_seq_io_recordset(open_gz_raw(path))
                } else {
                    count_seq_io_recordset(open_raw(path))
                })
            });
        });

        group.bench_function("needletail", |b| {
            b.iter(|| {
                black_box(if compressed {
                    count_needletail_reader(open_gz_raw(path))
                } else {
                    count_needletail_reader(open_raw(path))
                })
            });
        });

        // needletail's convenience entry point, which opens the path and sniffs
        // the compression itself. Reported separately because it is measuring
        // that plumbing, not the parser the row above isolates.
        group.bench_function("needletail_parse_fastx_file", |b| {
            b.iter(|| black_box(count_needletail(path)));
        });

        group.finish();
    }
}

/// Same parsers, plain input only, with every field touched. Isolates parsing
/// cost from decompression cost.
fn fastq_sequential_full(c: &mut Criterion) {
    let (plain, _, _) = data::reads_fastq();
    let plain_len = std::fs::metadata(&plain).unwrap().len();

    let mut group = c.benchmark_group("fastq_sequential_full");
    group.throughput(Throughput::Bytes(plain_len));
    group.sample_size(20);

    group.bench_function("noodles", |b| {
        b.iter(|| black_box(count_noodles_full(open_plain(&plain))));
    });
    group.bench_function("seq_io", |b| {
        b.iter(|| black_box(count_seq_io_full(open_raw(&plain))));
    });
    group.bench_function("needletail", |b| {
        b.iter(|| black_box(count_needletail_full(&plain)));
    });

    group.finish();
}

// ---------------------------------------------------------------------------
// Group 3: the layer underneath — raw decompression throughput
// ---------------------------------------------------------------------------

/// How much of the compressed-input numbers is decompression rather than
/// parsing? This reads the same BGZF file to EOF without parsing anything.
fn fastq_decompress(c: &mut Criterion) {
    let (plain, gz, bgz) = data::reads_fastq();
    let plain_len = std::fs::metadata(&plain).unwrap().len();

    let mut group = c.benchmark_group("fastq_decompress");
    group.throughput(Throughput::Bytes(plain_len));
    group.sample_size(20);

    let mut sink = vec![0u8; 1 << 16];

    group.bench_function("flate2_gzip", |b| {
        b.iter(|| {
            let mut r = flate2::read::MultiGzDecoder::new(File::open(&gz).unwrap());
            let mut n = 0usize;
            loop {
                let k = r.read(&mut sink).unwrap();
                if k == 0 {
                    break;
                }
                n += k;
            }
            black_box(n)
        });
    });

    // `flate2_gzip` reads the decoder directly; the parsers all put a buffer in
    // front of it. Both shapes are measured so the parse groups have a control
    // that matches how they actually consume the stream.
    group.bench_function("flate2_gzip_bufread", |b| {
        b.iter(|| {
            use std::io::BufRead as _;
            let mut r = open_gz(&gz);
            let mut n = 0usize;
            loop {
                let k = r.fill_buf().unwrap().len();
                if k == 0 {
                    break;
                }
                r.consume(k);
                n += k;
            }
            black_box(n)
        });
    });

    group.bench_function("flate2_on_bgzf", |b| {
        b.iter(|| {
            let mut r = flate2::read::MultiGzDecoder::new(File::open(&bgz).unwrap());
            let mut n = 0usize;
            loop {
                let k = r.read(&mut sink).unwrap();
                if k == 0 {
                    break;
                }
                n += k;
            }
            black_box(n)
        });
    });

    // Same decompressor, but consumed exactly the way the parsers consume it
    // (`BufRead::fill_buf`/`consume` through `open_gz`). Pins down whether the
    // read shape, rather than the codec, explains the gap to the parse group.
    group.bench_function("flate2_on_bgzf_bufread", |b| {
        b.iter(|| {
            use std::io::BufRead as _;
            let mut r = open_gz(&bgz);
            let mut n = 0usize;
            loop {
                let k = r.fill_buf().unwrap().len();
                if k == 0 {
                    break;
                }
                r.consume(k);
                n += k;
            }
            black_box(n)
        });
    });

    group.bench_function("seqair_bgzf", |b| {
        b.iter(|| {
            let mut r = seqair::bam::bgzf::BgzfReader::open(&bgz).unwrap();
            let mut n = 0usize;
            loop {
                let k = r.read_up_to(&mut sink).unwrap();
                if k == 0 {
                    break;
                }
                n += k;
            }
            black_box(n)
        });
    });

    group.bench_function("noodles_bgzf", |b| {
        b.iter(|| {
            let mut r = noodles_bgzf::io::Reader::new(File::open(&bgz).unwrap());
            let mut n = 0usize;
            loop {
                let k = r.read(&mut sink).unwrap();
                if k == 0 {
                    break;
                }
                n += k;
            }
            black_box(n)
        });
    });

    // Spawning the worker pool is done in the (untimed) setup step: creating a
    // reader per iteration would measure thread startup on a ~17 MB input
    // rather than steady-state throughput, and the pool cannot be reused
    // because the reader owns the file.
    group.bench_function("noodles_bgzf_mt", |b| {
        b.iter_batched(
            || noodles_bgzf::io::MultithreadedReader::new(File::open(&bgz).unwrap()),
            |mut r| {
                let mut local = vec![0u8; 1 << 16];
                let mut n = 0usize;
                loop {
                    let k = r.read(&mut local).unwrap();
                    if k == 0 {
                        break;
                    }
                    n += k;
                }
                black_box(n)
            },
            criterion::BatchSize::PerIteration,
        );
    });

    group.finish();
}

// ---------------------------------------------------------------------------
// Group 4: building the index
// ---------------------------------------------------------------------------

/// `noodles-fastq` can build the 6-column FAI even though it cannot query it.
/// This is the cost a caller pays once, if they have to index a FASTQ that
/// arrived without one.
///
/// Measured on the *unwrapped* read-like FASTQ, because that is the only kind
/// noodles can index: pointed at the line-wrapped `ref_fastq` — precisely the
/// shape a reference FASTQ has, and the shape `quality_scores_offset` exists to
/// describe — `fs::index` fails with `invalid name prefix`, having taken the
/// record's second sequence line for the next record's header.
fn fastq_index_build(c: &mut Criterion) {
    let (src, _, _) = data::reads_fastq();
    let len = std::fs::metadata(&src).unwrap().len();

    // Documents the limitation above; `samtools faidx` indexes this file fine.
    // Reported, not asserted, so a future noodles release that gains wrapped
    // support does not turn this bench into a panic.
    if noodles::fastq::fs::index(data::ref_fastq()).is_ok() {
        println!("note: noodles-fastq now indexes line-wrapped FASTQ; revisit the module docs");
    }

    let mut group = c.benchmark_group("fastq_index_build");
    group.throughput(Throughput::Bytes(len));
    group.sample_size(10);

    group.bench_function("noodles_fastq_fs_index", |b| {
        b.iter(|| black_box(noodles::fastq::fs::index(&src).unwrap().len()));
    });

    group.finish();
}

criterion_group!(
    benches,
    fastq_indexed_fetch,
    fastq_indexed_fetch_real,
    fastq_pileup_walk,
    fastq_sequential,
    fastq_sequential_full,
    fastq_decompress,
    fastq_index_build
);
criterion_main!(benches);
