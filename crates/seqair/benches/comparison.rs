//! Criterion benchmarks comparing seqair against htslib, noodles, and the bgzf crate.
//!
//! All libraries use libdeflate for BGZF decompression (fair comparison).
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::indexing_slicing, reason = "benches")]
#![allow(clippy::arithmetic_side_effects, reason = "benches")]
#![allow(clippy::cast_possible_truncation, reason = "benches")]
#![allow(clippy::cast_possible_wrap, reason = "benches")]

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use std::hint::black_box;
use std::io::Read;

const BAM_PATH: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.bam");
const FASTA_PATH: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.fasta.gz");

const CHROM: &str = "chr19";
const START: u64 = 6_105_000;
const END: u64 = 6_140_000;

fn start_pos() -> seqair::bam::Pos<seqair::bam::Zero> {
    seqair::bam::Pos::<seqair::bam::Zero>::new(START as u32).unwrap()
}
fn end_pos() -> seqair::bam::Pos<seqair::bam::Zero> {
    seqair::bam::Pos::<seqair::bam::Zero>::new(END as u32).unwrap()
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
            reader.fetch(FetchDefinition::Region(0, START as i64, END as i64)).unwrap();
            let mut record = bam::Record::new();
            let mut total_bases = 0usize;
            while reader.read(&mut record) == Some(Ok(())) {
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
            use seqair::bam::bgzf_writer::BgzfWriter;
            use seqair::bam::owned_record::OwnedBamRecord;
            use seqair::bam::writer::BamWriter;

            let path = std::path::Path::new(BAM_PATH);
            let mut reader = seqair::bam::IndexedBamReader::open(path).unwrap();
            let header = BamHeader::from_template(reader.header());
            let tid = reader.header().tid(CHROM).unwrap();
            let mut store = seqair::bam::RecordStore::new();
            reader.fetch_into(tid, start_pos(), end_pos(), &mut store).unwrap();

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
                    pos: i64::from(slim.pos.get()),
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
            // Read
            let mut reader = bam::IndexedReader::from_path(BAM_PATH).unwrap();
            reader.fetch(FetchDefinition::Region(0, START as i64, END as i64)).unwrap();
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
            // Read (sequential scan — noodles BAM indexed access requires extra setup)
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
                if aln_start >= START as i64 && aln_start < END as i64 {
                    records.push(record);
                }
            }

            // Write — noodles BAM writer handles header + ref seqs + records
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
            let seq = reader
                .fetch_seq(
                    CHROM,
                    seqair::bam::Pos::<seqair::bam::Zero>::new(FASTA_START as u32).unwrap(),
                    seqair::bam::Pos::<seqair::bam::Zero>::new(FASTA_END as u32).unwrap(),
                )
                .unwrap();
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
                let pos = u64::from(p.pos());
                if !(START..=END).contains(&pos) {
                    continue;
                }
                total_depth += u64::from(p.depth());
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

// ---------------------------------------------------------------------------
// Shared VCF writing helpers for benchmarks
//
// Models a realistic germline WGS callset: 10 samples, GATK-style fields
// (DP, AN, AC[A], AF[A], MQ, QD, FS), mixed variant types (SNV/indel/multi-
// allelic), mixed genotypes, and occasional LowQual filter failures.
// ---------------------------------------------------------------------------

use seqair::vcf::record_encoder::{FilterFieldDef, FormatFieldDef, InfoFieldDef};
use seqair::vcf::{
    Alleles, ContigDef, ContigId, FilterId, FormatGt, FormatInt, Genotype, InfoFloat, InfoFloats,
    InfoInt, InfoInts, Number, OutputFormat, ValueType, VcfHeader, Writer,
};
use seqair_types::Base;
use std::sync::Arc;

const N_RECORDS: u32 = 10_000;
const N_GERMLINE_SAMPLES: usize = 10;

/// Allele bank for htslib/noodles: (REF bytes, ALT1 bytes, optional ALT2 bytes).
/// Matches the seqair `GermlineSetup::alleles_bank` index-for-index.
static HTSLIB_ALLELES: &[(&[u8], &[u8], Option<&[u8]>)] = &[
    (b"A", b"T", None),       // SNV
    (b"C", b"G", None),       // SNV
    (b"G", b"A", None),       // SNV
    (b"T", b"C", None),       // SNV
    (b"A", b"C", None),       // SNV
    (b"AC", b"A", None),      // deletion (1 bp)
    (b"TAC", b"T", None),     // deletion (2 bp)
    (b"G", b"GAT", None),     // insertion (2 bp)
    (b"A", b"T", Some(b"C")), // multi-allelic SNV
    (b"G", b"T", None),       // SNV
];

struct GermlineSetup {
    header: Arc<VcfHeader>,
    contig: ContigId,
    lowqual_filter: FilterId,
    dp_info: InfoInt,
    an_info: InfoInt,
    ac_info: InfoInts,
    af_info: InfoFloats,
    mq_info: InfoFloat,
    qd_info: InfoFloat,
    fs_info: InfoFloat,
    gt_fmt: FormatGt,
    dp_fmt: FormatInt,
    gq_fmt: FormatInt,
    alleles_bank: Vec<Alleles>,
}

fn germline_setup() -> GermlineSetup {
    let mut builder = VcfHeader::builder();
    let contig = builder.register_contig("chr1", ContigDef { length: Some(250_000_000) }).unwrap();
    let lowqual_filter =
        builder.register_filter(&FilterFieldDef::new("LowQual", "Low quality")).unwrap();
    let dp_info = builder
        .register_info(&InfoFieldDef::new(
            "DP",
            Number::Count(1),
            ValueType::Integer,
            "Total read depth",
        ))
        .unwrap();
    let an_info = builder
        .register_info(&InfoFieldDef::new(
            "AN",
            Number::Count(1),
            ValueType::Integer,
            "Total allele count in genotypes",
        ))
        .unwrap();
    let ac_info = builder
        .register_info(&InfoFieldDef::new(
            "AC",
            Number::AlternateBases,
            ValueType::Integer,
            "Allele count in genotypes, for each ALT",
        ))
        .unwrap();
    let af_info = builder
        .register_info(&InfoFieldDef::new(
            "AF",
            Number::AlternateBases,
            ValueType::Float,
            "Allele frequency from GT counts, for each ALT",
        ))
        .unwrap();
    let mq_info = builder
        .register_info(&InfoFieldDef::new(
            "MQ",
            Number::Count(1),
            ValueType::Float,
            "RMS mapping quality",
        ))
        .unwrap();
    let qd_info = builder
        .register_info(&InfoFieldDef::new(
            "QD",
            Number::Count(1),
            ValueType::Float,
            "Variant confidence divided by unfiltered depth",
        ))
        .unwrap();
    let fs_info = builder
        .register_info(&InfoFieldDef::new(
            "FS",
            Number::Count(1),
            ValueType::Float,
            "Phred-scaled p-value using Fisher's exact test for strand bias",
        ))
        .unwrap();
    let gt_fmt = builder
        .register_format(&FormatFieldDef::new(
            "GT",
            Number::Count(1),
            ValueType::String,
            "Genotype",
        ))
        .unwrap();
    let dp_fmt = builder
        .register_format(&FormatFieldDef::new(
            "DP",
            Number::Count(1),
            ValueType::Integer,
            "Approximate read depth",
        ))
        .unwrap();
    let gq_fmt = builder
        .register_format(&FormatFieldDef::new(
            "GQ",
            Number::Count(1),
            ValueType::Integer,
            "Genotype quality",
        ))
        .unwrap();
    let mut builder = builder.add_sample("SAMPLE01").unwrap();
    for s in 2..=N_GERMLINE_SAMPLES {
        builder = builder.add_sample(&format!("SAMPLE{s:02}")).unwrap();
    }
    let header = Arc::new(builder.build().unwrap());
    let alleles_bank = vec![
        Alleles::snv(Base::A, Base::T).unwrap(),
        Alleles::snv(Base::C, Base::G).unwrap(),
        Alleles::snv(Base::G, Base::A).unwrap(),
        Alleles::snv(Base::T, Base::C).unwrap(),
        Alleles::snv(Base::A, Base::C).unwrap(),
        Alleles::deletion(Base::A, &[Base::C]).unwrap(),
        Alleles::deletion(Base::T, &[Base::A, Base::C]).unwrap(),
        Alleles::insertion(Base::G, &[Base::A, Base::T]).unwrap(),
        Alleles::snv_multi(Base::A, &[Base::T, Base::C]).unwrap(),
        Alleles::snv(Base::G, Base::T).unwrap(),
    ];
    GermlineSetup {
        header,
        contig,
        lowqual_filter,
        dp_info,
        an_info,
        ac_info,
        af_info,
        mq_info,
        qd_info,
        fs_info,
        gt_fmt,
        dp_fmt,
        gq_fmt,
        alleles_bank,
    }
}

fn htslib_germline_header() -> rust_htslib::bcf::Header {
    let mut h = rust_htslib::bcf::Header::new();
    h.push_record(b"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth\">");
    h.push_record(
        b"##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total allele count in genotypes\">",
    );
    h.push_record(
        b"##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes, for each ALT\">",
    );
    h.push_record(
        b"##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency from GT counts, for each ALT\">",
    );
    h.push_record(b"##INFO=<ID=MQ,Number=1,Type=Float,Description=\"RMS mapping quality\">");
    h.push_record(
        b"##INFO=<ID=QD,Number=1,Type=Float,Description=\"Variant confidence divided by unfiltered depth\">",
    );
    h.push_record(
        b"##INFO=<ID=FS,Number=1,Type=Float,Description=\"Phred-scaled p-value using Fisher's exact test for strand bias\">",
    );
    h.push_record(b"##FILTER=<ID=LowQual,Description=\"Low quality\">");
    h.push_record(b"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    h.push_record(b"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth\">");
    h.push_record(b"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">");
    h.push_record(b"##contig=<ID=chr1,length=250000000>");
    for s in 1..=N_GERMLINE_SAMPLES {
        h.push_sample(format!("SAMPLE{s:02}").as_bytes());
    }
    h
}

/// Write N_RECORDS germline-style records with rust-htslib to a temp file.
fn write_htslib_germline(format: rust_htslib::bcf::Format) -> tempfile::NamedTempFile {
    use rust_htslib::bcf;
    use rust_htslib::bcf::record::GenotypeAllele;

    let tmp = tempfile::NamedTempFile::new().unwrap();
    let hdr = htslib_germline_header();
    let mut writer = bcf::Writer::from_path(tmp.path(), &hdr, false, format).unwrap();
    let rid = writer.header().name2rid(b"chr1").unwrap();
    let pass_id = writer.header().name_to_id(cstr8::cstr8!("PASS")).unwrap();
    let lowqual_id = writer.header().name_to_id(cstr8::cstr8!("LowQual")).unwrap();
    let an = (N_GERMLINE_SAMPLES * 2) as i32;

    for i in 0..i64::from(N_RECORDS) {
        let allele_idx = i as usize % HTSLIB_ALLELES.len();
        let (ref_bytes, alt1_bytes, alt2_opt) = HTSLIB_ALLELES[allele_idx];
        let n_alt = if alt2_opt.is_some() { 2usize } else { 1 };

        let mut record = writer.empty_record();
        record.set_rid(Some(rid));
        record.set_pos(i * 237);
        record.set_qual(10.0 + (i % 50) as f32 * 2.0);
        let alleles: Vec<&[u8]> = match alt2_opt {
            Some(alt2) => vec![ref_bytes, alt1_bytes, alt2],
            None => vec![ref_bytes, alt1_bytes],
        };
        record.set_alleles(&alleles).unwrap();
        if i % 10 == 7 {
            record.push_filter(&lowqual_id).unwrap();
        } else {
            record.push_filter(&pass_id).unwrap();
        }
        let dp = 30 + (i % 70) as i32;
        let ac: Vec<i32> = (0..n_alt).map(|k| 1 + ((i + k as i64) % 9) as i32).collect();
        let af: Vec<f32> = ac.iter().map(|&a| a as f32 / an as f32).collect();
        record.push_info_integer(cstr8::cstr8!("DP"), &[dp]).unwrap();
        record.push_info_integer(cstr8::cstr8!("AN"), &[an]).unwrap();
        record.push_info_integer(cstr8::cstr8!("AC"), &ac).unwrap();
        record.push_info_float(cstr8::cstr8!("AF"), &af).unwrap();
        record.push_info_float(cstr8::cstr8!("MQ"), &[55.0f32 + (i % 20) as f32 * 0.5]).unwrap();
        record.push_info_float(cstr8::cstr8!("QD"), &[5.0f32 + (i % 25) as f32 * 0.8]).unwrap();
        record.push_info_float(cstr8::cstr8!("FS"), &[(i % 30) as f32 * 0.5]).unwrap();
        // Genotypes: N_GERMLINE_SAMPLES × 2 alleles flattened
        let gts: Vec<GenotypeAllele> = (0..N_GERMLINE_SAMPLES)
            .flat_map(|s| match (s + i as usize) % 6 {
                0 => [GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(0)],
                1 => [GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)],
                2 => [GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(1)],
                3 => [GenotypeAllele::Phased(0), GenotypeAllele::Phased(1)],
                4 => [GenotypeAllele::UnphasedMissing, GenotypeAllele::UnphasedMissing],
                _ => [GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)],
            })
            .collect();
        record.push_genotypes(&gts).unwrap();
        let dp_vals: Vec<i32> = (0..N_GERMLINE_SAMPLES as i32).map(|s| dp / 2 + s % 20).collect();
        let gq_vals: Vec<i32> =
            (0..N_GERMLINE_SAMPLES as i32).map(|s| 30 + (i as i32 + s) % 60).collect();
        record.push_format_integer(cstr8::cstr8!("DP"), &dp_vals).unwrap();
        record.push_format_integer(cstr8::cstr8!("GQ"), &gq_vals).unwrap();
        writer.write(&record).unwrap();
    }
    drop(writer);
    tmp
}

// ---------------------------------------------------------------------------
// Group 5: Write VCF text (to memory buffer)
// seqair vs htslib vs noodles — plain uncompressed VCF, 10 samples
// ---------------------------------------------------------------------------

fn write_vcf_text(c: &mut Criterion) {
    let mut group = c.benchmark_group("write_vcf_text");
    group.throughput(Throughput::Elements(u64::from(N_RECORDS)));

    let setup = germline_setup();

    // seqair: typestate encoder, zero-alloc, pre-resolved keys
    group.bench_function("seqair", |b| {
        let GermlineSetup {
            header,
            contig,
            lowqual_filter,
            dp_info,
            an_info,
            ac_info,
            af_info,
            mq_info,
            qd_info,
            fs_info,
            gt_fmt,
            dp_fmt,
            gq_fmt,
            alleles_bank,
        } = &setup;
        let an = (N_GERMLINE_SAMPLES * 2) as i32;
        b.iter(|| {
            let mut output = Vec::with_capacity(4_000_000);
            let writer = Writer::new(&mut output, OutputFormat::Vcf);
            let mut writer = writer.write_header(header).unwrap();
            for i in 1u32..=N_RECORDS {
                let alleles = &alleles_bank[i as usize % alleles_bank.len()];
                let n_alt = alleles.n_allele() - 1;
                let pos = seqair_types::Pos1::new(i * 237 + 1).unwrap();
                let qual = Some(10.0 + (i % 50) as f32 * 2.0);
                let enc = writer.begin_record(contig, pos, alleles, qual).unwrap();
                let mut enc = if i % 10 == 7 {
                    enc.filter_fail(&[lowqual_filter])
                } else {
                    enc.filter_pass()
                };
                let dp = 30 + (i % 70) as i32;
                dp_info.encode(&mut enc, dp);
                an_info.encode(&mut enc, an);
                let ac: Vec<i32> = (0..n_alt).map(|k| 1 + (i as i32 + k as i32) % 9).collect();
                let af: Vec<f32> = ac.iter().map(|&a| a as f32 / an as f32).collect();
                ac_info.encode(&mut enc, &ac);
                af_info.encode(&mut enc, &af);
                mq_info.encode(&mut enc, 55.0 + (i % 20) as f32 * 0.5);
                qd_info.encode(&mut enc, 5.0 + (i % 25) as f32 * 0.8);
                fs_info.encode(&mut enc, (i % 30) as f32 * 0.5);
                let mut enc = enc.begin_samples(N_GERMLINE_SAMPLES as u32);
                for s in 0..N_GERMLINE_SAMPLES {
                    let gt = match (s + i as usize) % 6 {
                        0 => Genotype::unphased(0, 0),
                        1 => Genotype::unphased(0, 1),
                        2 => Genotype::unphased(1, 1),
                        3 => Genotype::phased_diploid(0, 1),
                        4 => Genotype::missing_diploid(),
                        _ => Genotype::unphased(0, 1),
                    };
                    gt_fmt.encode(&mut enc, &gt);
                    dp_fmt.encode(&mut enc, dp / 2 + s as i32 % 20);
                    gq_fmt.encode(&mut enc, 30 + (i as i32 + s as i32) % 60);
                }
                enc.emit().unwrap();
            }
            writer.finish().unwrap();
            black_box(output.len())
        });
    });

    // htslib: writes to temp file (htslib VCF writer requires a path)
    group.bench_function("htslib", |b| {
        b.iter(|| {
            let tmp = write_htslib_germline(rust_htslib::bcf::Format::Vcf);
            black_box(std::fs::metadata(tmp.path()).unwrap().len())
        });
    });

    // noodles: VCF text writer via RecordBuf
    group.bench_function("noodles", |b| {
        use noodles::vcf;
        use noodles::vcf::variant::io::Write as _;
        use noodles::vcf::variant::record_buf::info::field::Value as InfoValue;
        use noodles::vcf::variant::record_buf::info::field::value::Array as InfoArray;
        use noodles::vcf::variant::record_buf::samples::sample::Value as SampleValue;

        let noodles_header: vcf::Header = setup.header.to_vcf_text().parse().unwrap();
        let an = (N_GERMLINE_SAMPLES * 2) as i32;

        b.iter(|| {
            let mut output = Vec::with_capacity(4_000_000);
            let mut writer = vcf::io::Writer::new(&mut output);
            writer.write_header(&noodles_header).unwrap();

            for i in 1u32..=N_RECORDS {
                let allele_idx = i as usize % HTSLIB_ALLELES.len();
                let (ref_bytes, alt1_bytes, alt2_opt) = HTSLIB_ALLELES[allele_idx];
                let n_alt = if alt2_opt.is_some() { 2usize } else { 1 };
                let dp = 30 + (i % 70) as i32;
                let ac: Vec<Option<i32>> =
                    (0..n_alt).map(|k| Some(1 + (i as i32 + k as i32) % 9)).collect();
                let af: Vec<Option<f32>> =
                    ac.iter().map(|&a| a.map(|v| v as f32 / an as f32)).collect();

                let pos = noodles::core::Position::try_from(i as usize * 237 + 1).unwrap();
                let mut alt_bases = vec![std::str::from_utf8(alt1_bytes).unwrap().to_owned()];
                if let Some(alt2) = alt2_opt {
                    alt_bases.push(std::str::from_utf8(alt2).unwrap().to_owned());
                }
                let filters = if i % 10 == 7 {
                    [String::from("LowQual")].into_iter().collect()
                } else {
                    vcf::variant::record_buf::Filters::pass()
                };
                let info = [
                    ("DP".parse().unwrap(), Some(InfoValue::Integer(dp))),
                    ("AN".parse().unwrap(), Some(InfoValue::Integer(an))),
                    ("AC".parse().unwrap(), Some(InfoValue::Array(InfoArray::Integer(ac)))),
                    ("AF".parse().unwrap(), Some(InfoValue::Array(InfoArray::Float(af)))),
                    ("MQ".parse().unwrap(), Some(InfoValue::Float(55.0 + (i % 20) as f32 * 0.5))),
                    ("QD".parse().unwrap(), Some(InfoValue::Float(5.0 + (i % 25) as f32 * 0.8))),
                    ("FS".parse().unwrap(), Some(InfoValue::Float((i % 30) as f32 * 0.5))),
                ]
                .into_iter()
                .collect();
                let format_keys =
                    ["GT".parse().unwrap(), "DP".parse().unwrap(), "GQ".parse().unwrap()]
                        .into_iter()
                        .collect();
                let sample_rows = (0..N_GERMLINE_SAMPLES)
                    .map(|s| {
                        let gt_str = match (s + i as usize) % 6 {
                            0 => "0/0",
                            1 => "0/1",
                            2 => "1/1",
                            3 => "0|1",
                            4 => "./.",
                            _ => "0/1",
                        };
                        [
                            Some(SampleValue::Genotype(gt_str.parse().unwrap())),
                            Some(SampleValue::Integer(dp / 2 + s as i32 % 20)),
                            Some(SampleValue::Integer(30 + (i as i32 + s as i32) % 60)),
                        ]
                        .into_iter()
                        .collect()
                    })
                    .collect();
                let record = vcf::variant::RecordBuf::builder()
                    .set_reference_sequence_name(std::str::from_utf8(ref_bytes).unwrap())
                    .set_variant_start(pos)
                    .set_reference_bases(std::str::from_utf8(ref_bytes).unwrap())
                    .set_alternate_bases(vcf::variant::record_buf::AlternateBases::from(alt_bases))
                    .set_quality_score(10.0 + (i % 50) as f32 * 2.0)
                    .set_filters(filters)
                    .set_info(info)
                    .set_samples(vcf::variant::record_buf::Samples::new(format_keys, sample_rows))
                    .build();
                writer.write_variant_record(&noodles_header, &record).unwrap();
            }
            black_box(output.len())
        });
    });

    group.finish();
}

// ---------------------------------------------------------------------------
// Group 6: Write .vcf.gz (BGZF-compressed VCF), 10 samples
// seqair vs htslib — tests BGZF compression overhead on top of encoding
// ---------------------------------------------------------------------------

fn write_vcf_gz(c: &mut Criterion) {
    let mut group = c.benchmark_group("write_vcf_gz");
    group.throughput(Throughput::Elements(u64::from(N_RECORDS)));

    let setup = germline_setup();

    // seqair: BGZF-compressed VCF to memory buffer
    group.bench_function("seqair", |b| {
        let GermlineSetup {
            header,
            contig,
            lowqual_filter,
            dp_info,
            an_info,
            ac_info,
            af_info,
            mq_info,
            qd_info,
            fs_info,
            gt_fmt,
            dp_fmt,
            gq_fmt,
            alleles_bank,
        } = &setup;
        let an = (N_GERMLINE_SAMPLES * 2) as i32;
        b.iter(|| {
            let mut output = Vec::with_capacity(2_000_000);
            let writer = Writer::new(&mut output, OutputFormat::VcfGz);
            let mut writer = writer.write_header(header).unwrap();
            for i in 1u32..=N_RECORDS {
                let alleles = &alleles_bank[i as usize % alleles_bank.len()];
                let n_alt = alleles.n_allele() - 1;
                let pos = seqair_types::Pos1::new(i * 237 + 1).unwrap();
                let enc = writer
                    .begin_record(contig, pos, alleles, Some(10.0 + (i % 50) as f32 * 2.0))
                    .unwrap();
                let mut enc = if i % 10 == 7 {
                    enc.filter_fail(&[lowqual_filter])
                } else {
                    enc.filter_pass()
                };
                let dp = 30 + (i % 70) as i32;
                dp_info.encode(&mut enc, dp);
                an_info.encode(&mut enc, an);
                let ac: Vec<i32> = (0..n_alt).map(|k| 1 + (i as i32 + k as i32) % 9).collect();
                let af: Vec<f32> = ac.iter().map(|&a| a as f32 / an as f32).collect();
                ac_info.encode(&mut enc, &ac);
                af_info.encode(&mut enc, &af);
                mq_info.encode(&mut enc, 55.0 + (i % 20) as f32 * 0.5);
                qd_info.encode(&mut enc, 5.0 + (i % 25) as f32 * 0.8);
                fs_info.encode(&mut enc, (i % 30) as f32 * 0.5);
                let mut enc = enc.begin_samples(N_GERMLINE_SAMPLES as u32);
                for s in 0..N_GERMLINE_SAMPLES {
                    let gt = match (s + i as usize) % 6 {
                        0 => Genotype::unphased(0, 0),
                        1 => Genotype::unphased(0, 1),
                        2 => Genotype::unphased(1, 1),
                        3 => Genotype::phased_diploid(0, 1),
                        4 => Genotype::missing_diploid(),
                        _ => Genotype::unphased(0, 1),
                    };
                    gt_fmt.encode(&mut enc, &gt);
                    dp_fmt.encode(&mut enc, dp / 2 + s as i32 % 20);
                    gq_fmt.encode(&mut enc, 30 + (i as i32 + s as i32) % 60);
                }
                enc.emit().unwrap();
            }
            writer.finish().unwrap();
            black_box(output.len())
        });
    });

    // htslib: writes bgzf-compressed VCF to temp file
    group.bench_function("htslib", |b| {
        b.iter(|| {
            let tmp = tempfile::Builder::new().suffix(".vcf.gz").tempfile().unwrap();
            use rust_htslib::bcf;
            use rust_htslib::bcf::record::GenotypeAllele;
            let hdr = htslib_germline_header();
            let mut writer =
                bcf::Writer::from_path(tmp.path(), &hdr, false, bcf::Format::Vcf).unwrap();
            let rid = writer.header().name2rid(b"chr1").unwrap();
            let pass_id = writer.header().name_to_id(cstr8::cstr8!("PASS")).unwrap();
            let lowqual_id = writer.header().name_to_id(cstr8::cstr8!("LowQual")).unwrap();
            let an = (N_GERMLINE_SAMPLES * 2) as i32;
            for i in 0..i64::from(N_RECORDS) {
                let allele_idx = i as usize % HTSLIB_ALLELES.len();
                let (ref_bytes, alt1_bytes, alt2_opt) = HTSLIB_ALLELES[allele_idx];
                let n_alt = if alt2_opt.is_some() { 2usize } else { 1 };
                let mut record = writer.empty_record();
                record.set_rid(Some(rid));
                record.set_pos(i * 237);
                record.set_qual(10.0 + (i % 50) as f32 * 2.0);
                let alleles: Vec<&[u8]> = match alt2_opt {
                    Some(alt2) => vec![ref_bytes, alt1_bytes, alt2],
                    None => vec![ref_bytes, alt1_bytes],
                };
                record.set_alleles(&alleles).unwrap();
                if i % 10 == 7 {
                    record.push_filter(&lowqual_id).unwrap();
                } else {
                    record.push_filter(&pass_id).unwrap();
                }
                let dp = 30 + (i % 70) as i32;
                let ac: Vec<i32> = (0..n_alt).map(|k| 1 + ((i + k as i64) % 9) as i32).collect();
                let af: Vec<f32> = ac.iter().map(|&a| a as f32 / an as f32).collect();
                record.push_info_integer(cstr8::cstr8!("DP"), &[dp]).unwrap();
                record.push_info_integer(cstr8::cstr8!("AN"), &[an]).unwrap();
                record.push_info_integer(cstr8::cstr8!("AC"), &ac).unwrap();
                record.push_info_float(cstr8::cstr8!("AF"), &af).unwrap();
                record
                    .push_info_float(cstr8::cstr8!("MQ"), &[55.0f32 + (i % 20) as f32 * 0.5])
                    .unwrap();
                record
                    .push_info_float(cstr8::cstr8!("QD"), &[5.0f32 + (i % 25) as f32 * 0.8])
                    .unwrap();
                record.push_info_float(cstr8::cstr8!("FS"), &[(i % 30) as f32 * 0.5]).unwrap();
                let gts: Vec<GenotypeAllele> = (0..N_GERMLINE_SAMPLES)
                    .flat_map(|s| match (s + i as usize) % 6 {
                        0 => [GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(0)],
                        1 => [GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)],
                        2 => [GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(1)],
                        3 => [GenotypeAllele::Phased(0), GenotypeAllele::Phased(1)],
                        4 => [GenotypeAllele::UnphasedMissing, GenotypeAllele::UnphasedMissing],
                        _ => [GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)],
                    })
                    .collect();
                record.push_genotypes(&gts).unwrap();
                let dp_vals: Vec<i32> =
                    (0..N_GERMLINE_SAMPLES as i32).map(|s| dp / 2 + s % 20).collect();
                let gq_vals: Vec<i32> =
                    (0..N_GERMLINE_SAMPLES as i32).map(|s| 30 + (i as i32 + s) % 60).collect();
                record.push_format_integer(cstr8::cstr8!("DP"), &dp_vals).unwrap();
                record.push_format_integer(cstr8::cstr8!("GQ"), &gq_vals).unwrap();
                writer.write(&record).unwrap();
            }
            drop(writer);
            black_box(std::fs::metadata(tmp.path()).unwrap().len())
        });
    });

    group.finish();
}

// ---------------------------------------------------------------------------
// Group 7: Write .bcf (binary), 10 samples
// seqair (typestate encoder, zero-alloc) vs htslib vs noodles
// ---------------------------------------------------------------------------

fn write_bcf(c: &mut Criterion) {
    let mut group = c.benchmark_group("write_bcf");
    group.throughput(Throughput::Elements(u64::from(N_RECORDS)));

    let setup = germline_setup();

    // seqair: typestate encoder path (zero-alloc, pre-resolved dict indices)
    group.bench_function("seqair", |b| {
        let GermlineSetup {
            header,
            contig,
            lowqual_filter,
            dp_info,
            an_info,
            ac_info,
            af_info,
            mq_info,
            qd_info,
            fs_info,
            gt_fmt,
            dp_fmt,
            gq_fmt,
            alleles_bank,
        } = &setup;
        let an = (N_GERMLINE_SAMPLES * 2) as i32;
        b.iter(|| {
            let mut output = Vec::with_capacity(2_000_000);
            let writer = Writer::new(&mut output, OutputFormat::Bcf);
            let mut writer = writer.write_header(header).unwrap();
            for i in 1u32..=N_RECORDS {
                let alleles = &alleles_bank[i as usize % alleles_bank.len()];
                let n_alt = alleles.n_allele() - 1;
                let pos = seqair_types::Pos1::new(i * 237 + 1).unwrap();
                let enc = writer
                    .begin_record(contig, pos, alleles, Some(10.0 + (i % 50) as f32 * 2.0))
                    .unwrap();
                let mut enc = if i % 10 == 7 {
                    enc.filter_fail(&[lowqual_filter])
                } else {
                    enc.filter_pass()
                };
                let dp = 30 + (i % 70) as i32;
                dp_info.encode(&mut enc, dp);
                an_info.encode(&mut enc, an);
                let ac: Vec<i32> = (0..n_alt).map(|k| 1 + (i as i32 + k as i32) % 9).collect();
                let af: Vec<f32> = ac.iter().map(|&a| a as f32 / an as f32).collect();
                ac_info.encode(&mut enc, &ac);
                af_info.encode(&mut enc, &af);
                mq_info.encode(&mut enc, 55.0 + (i % 20) as f32 * 0.5);
                qd_info.encode(&mut enc, 5.0 + (i % 25) as f32 * 0.8);
                fs_info.encode(&mut enc, (i % 30) as f32 * 0.5);
                let mut enc = enc.begin_samples(N_GERMLINE_SAMPLES as u32);
                for s in 0..N_GERMLINE_SAMPLES {
                    let gt = match (s + i as usize) % 6 {
                        0 => Genotype::unphased(0, 0),
                        1 => Genotype::unphased(0, 1),
                        2 => Genotype::unphased(1, 1),
                        3 => Genotype::phased_diploid(0, 1),
                        4 => Genotype::missing_diploid(),
                        _ => Genotype::unphased(0, 1),
                    };
                    gt_fmt.encode(&mut enc, &gt);
                    dp_fmt.encode(&mut enc, dp / 2 + s as i32 % 20);
                    gq_fmt.encode(&mut enc, 30 + (i as i32 + s as i32) % 60);
                }
                enc.emit().unwrap();
            }
            writer.finish().unwrap();
            black_box(output.len())
        });
    });

    // htslib
    group.bench_function("htslib", |b| {
        b.iter(|| {
            let tmp = write_htslib_germline(rust_htslib::bcf::Format::Bcf);
            black_box(std::fs::metadata(tmp.path()).unwrap().len())
        });
    });

    // noodles BCF writer
    group.bench_function("noodles", |b| {
        use noodles::bcf;
        use noodles::vcf;
        use noodles::vcf::variant::io::Write as _;
        use noodles::vcf::variant::record_buf::info::field::Value as InfoValue;
        use noodles::vcf::variant::record_buf::info::field::value::Array as InfoArray;
        use noodles::vcf::variant::record_buf::samples::sample::Value as SampleValue;

        let noodles_header: vcf::Header = setup.header.to_vcf_text().parse().unwrap();
        let an = (N_GERMLINE_SAMPLES * 2) as i32;

        b.iter(|| {
            let mut output = Vec::with_capacity(2_000_000);
            let mut writer = bcf::io::Writer::new(&mut output);
            writer.write_header(&noodles_header).unwrap();

            for i in 1u32..=N_RECORDS {
                let allele_idx = i as usize % HTSLIB_ALLELES.len();
                let (ref_bytes, alt1_bytes, alt2_opt) = HTSLIB_ALLELES[allele_idx];
                let n_alt = if alt2_opt.is_some() { 2usize } else { 1 };
                let dp = 30 + (i % 70) as i32;
                let ac: Vec<Option<i32>> =
                    (0..n_alt).map(|k| Some(1 + (i as i32 + k as i32) % 9)).collect();
                let af: Vec<Option<f32>> =
                    ac.iter().map(|&a| a.map(|v| v as f32 / an as f32)).collect();

                let pos = noodles::core::Position::try_from(i as usize * 237 + 1).unwrap();
                let mut alt_bases = vec![std::str::from_utf8(alt1_bytes).unwrap().to_owned()];
                if let Some(alt2) = alt2_opt {
                    alt_bases.push(std::str::from_utf8(alt2).unwrap().to_owned());
                }
                let filters = if i % 10 == 7 {
                    [String::from("LowQual")].into_iter().collect()
                } else {
                    vcf::variant::record_buf::Filters::pass()
                };
                let info = [
                    ("DP".parse().unwrap(), Some(InfoValue::Integer(dp))),
                    ("AN".parse().unwrap(), Some(InfoValue::Integer(an))),
                    ("AC".parse().unwrap(), Some(InfoValue::Array(InfoArray::Integer(ac)))),
                    ("AF".parse().unwrap(), Some(InfoValue::Array(InfoArray::Float(af)))),
                    ("MQ".parse().unwrap(), Some(InfoValue::Float(55.0 + (i % 20) as f32 * 0.5))),
                    ("QD".parse().unwrap(), Some(InfoValue::Float(5.0 + (i % 25) as f32 * 0.8))),
                    ("FS".parse().unwrap(), Some(InfoValue::Float((i % 30) as f32 * 0.5))),
                ]
                .into_iter()
                .collect();
                let format_keys =
                    ["GT".parse().unwrap(), "DP".parse().unwrap(), "GQ".parse().unwrap()]
                        .into_iter()
                        .collect();
                let sample_rows = (0..N_GERMLINE_SAMPLES)
                    .map(|s| {
                        let gt_str = match (s + i as usize) % 6 {
                            0 => "0/0",
                            1 => "0/1",
                            2 => "1/1",
                            3 => "0|1",
                            4 => "./.",
                            _ => "0/1",
                        };
                        [
                            Some(SampleValue::Genotype(gt_str.parse().unwrap())),
                            Some(SampleValue::Integer(dp / 2 + s as i32 % 20)),
                            Some(SampleValue::Integer(30 + (i as i32 + s as i32) % 60)),
                        ]
                        .into_iter()
                        .collect()
                    })
                    .collect();
                let record = vcf::variant::RecordBuf::builder()
                    .set_reference_sequence_name(std::str::from_utf8(ref_bytes).unwrap())
                    .set_variant_start(pos)
                    .set_reference_bases(std::str::from_utf8(ref_bytes).unwrap())
                    .set_alternate_bases(vcf::variant::record_buf::AlternateBases::from(alt_bases))
                    .set_quality_score(10.0 + (i % 50) as f32 * 2.0)
                    .set_filters(filters)
                    .set_info(info)
                    .set_samples(vcf::variant::record_buf::Samples::new(format_keys, sample_rows))
                    .build();
                writer.write_variant_record(&noodles_header, &record).unwrap();
            }
            drop(writer);
            black_box(output.len())
        });
    });

    group.finish();
}

criterion_group!(
    benches,
    bgzf_decompress,
    bam_record_decode,
    bam_roundtrip,
    fasta_fetch,
    pileup_e2e,
    write_vcf_text,
    write_vcf_gz,
    write_bcf,
);
criterion_main!(benches);
