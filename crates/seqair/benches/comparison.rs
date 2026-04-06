//! Criterion benchmarks comparing seqair against htslib, noodles, and the bgzf crate.
//!
//! All libraries use libdeflate for BGZF decompression (fair comparison).
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::indexing_slicing)]
#![allow(clippy::arithmetic_side_effects)]

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

// ---------------------------------------------------------------------------
// Shared VCF writing helpers for benchmarks
// ---------------------------------------------------------------------------

use seqair::vcf::alleles::Alleles;
use seqair::vcf::bcf_writer::BcfWriter;
use seqair::vcf::encoder::*;
use seqair::vcf::header::*;
use seqair::vcf::record::{Genotype, SampleValue, VcfRecordBuilder};
use seqair::vcf::writer::VcfWriter;
use seqair_types::{Base, SmolStr};
use std::marker::PhantomData;
use std::sync::Arc;

const N_RECORDS: u32 = 10_000;

fn vcf_header() -> Arc<VcfHeader> {
    Arc::new(
        VcfHeader::builder()
            .add_contig("chr1", ContigDef { length: Some(250_000_000) })
            .unwrap()
            .add_info(
                "DP",
                InfoDef {
                    number: Number::Count(1),
                    typ: ValueType::Integer,
                    description: SmolStr::from("Depth"),
                },
            )
            .unwrap()
            .add_format(
                "GT",
                FormatDef {
                    number: Number::Count(1),
                    typ: ValueType::String,
                    description: SmolStr::from("Genotype"),
                },
            )
            .unwrap()
            .add_format(
                "DP",
                FormatDef {
                    number: Number::Count(1),
                    typ: ValueType::Integer,
                    description: SmolStr::from("Read Depth"),
                },
            )
            .unwrap()
            .add_sample("sample1")
            .unwrap()
            .build()
            .unwrap(),
    )
}

fn htslib_header() -> rust_htslib::bcf::Header {
    let mut h = rust_htslib::bcf::Header::new();
    h.push_record(b"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">");
    h.push_record(b"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    h.push_record(b"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">");
    h.push_record(b"##contig=<ID=chr1,length=250000000>");
    h.push_sample(b"sample1");
    h
}

/// Write N_RECORDS with rust-htslib to a temp file, return the file.
fn write_htslib(format: rust_htslib::bcf::Format) -> tempfile::NamedTempFile {
    use rust_htslib::bcf;

    let tmp = tempfile::NamedTempFile::new().unwrap();
    let hdr = htslib_header();
    let mut writer = bcf::Writer::from_path(tmp.path(), &hdr, false, format).unwrap();
    let rid = writer.header().name2rid(b"chr1").unwrap();
    let pass_id = writer.header().name_to_id(cstr8::cstr8!("PASS")).unwrap();

    for i in 0..N_RECORDS as i64 {
        let mut record = writer.empty_record();
        record.set_rid(Some(rid));
        record.set_pos(i);
        record.set_qual(30.0);
        record.set_alleles(&[b"A", b"T"]).unwrap();
        record.push_filter(&pass_id).unwrap();
        record.push_info_integer(cstr8::cstr8!("DP"), &[50]).unwrap();
        record
            .push_genotypes(&[
                bcf::record::GenotypeAllele::Unphased(0),
                bcf::record::GenotypeAllele::Unphased(1),
            ])
            .unwrap();
        record.push_format_integer(cstr8::cstr8!("DP"), &[30]).unwrap();
        writer.write(&record).unwrap();
    }
    drop(writer);
    tmp
}

// ---------------------------------------------------------------------------
// Group 5: Write VCF text (to memory buffer)
// seqair vs htslib vs noodles — plain uncompressed VCF
// ---------------------------------------------------------------------------

fn write_vcf_text(c: &mut Criterion) {
    let mut group = c.benchmark_group("write_vcf_text");
    group.throughput(Throughput::Elements(N_RECORDS as u64));

    let header = vcf_header();

    // seqair
    group.bench_function("seqair", |b| {
        b.iter(|| {
            let mut output = Vec::with_capacity(1_000_000);
            let mut writer = VcfWriter::new(&mut output, header.clone());
            writer.write_header().unwrap();
            for i in 1..=N_RECORDS {
                let alleles = Alleles::snv(Base::A, Base::T).unwrap();
                let record =
                    VcfRecordBuilder::new("chr1", seqair_types::Pos1::new(i).unwrap(), alleles)
                        .qual(30.0)
                        .filter_pass()
                        .info_integer("DP", 50)
                        .format_keys(&["GT", "DP"])
                        .add_sample({
                            use seqair_types::smallvec::smallvec;
                            smallvec![
                                SampleValue::Genotype(Genotype::unphased(0, 1)),
                                SampleValue::Integer(30),
                            ]
                        })
                        .build(&header)
                        .unwrap();
                writer.write_record(&record).unwrap();
            }
            writer.finish().unwrap();
            black_box(output.len())
        });
    });

    // htslib (writes to temp file, VCF format)
    group.bench_function("htslib", |b| {
        b.iter(|| {
            let tmp = write_htslib(rust_htslib::bcf::Format::Vcf);
            black_box(std::fs::metadata(tmp.path()).unwrap().len())
        });
    });

    // noodles: VCF text writer
    group.bench_function("noodles", |b| {
        use noodles::vcf;
        use noodles::vcf::variant::io::Write as _;
        use noodles::vcf::variant::record_buf::info::field::Value as InfoValue;
        use noodles::vcf::variant::record_buf::samples::sample::Value as SampleValue;

        let noodles_header: vcf::Header = header.to_vcf_text().parse().unwrap();

        b.iter(|| {
            let mut output = Vec::with_capacity(1_000_000);
            let mut writer = vcf::io::Writer::new(&mut output);
            writer.write_header(&noodles_header).unwrap();

            for i in 1..=N_RECORDS {
                let pos = noodles::core::Position::try_from(i as usize).unwrap();
                let record = vcf::variant::RecordBuf::builder()
                    .set_reference_sequence_name("chr1")
                    .set_variant_start(pos)
                    .set_reference_bases("A")
                    .set_alternate_bases(vcf::variant::record_buf::AlternateBases::from(vec![
                        "T".to_string(),
                    ]))
                    .set_quality_score(30.0)
                    .set_filters(vcf::variant::record_buf::Filters::pass())
                    .set_info(
                        [("DP".parse().unwrap(), Some(InfoValue::Integer(50)))]
                            .into_iter()
                            .collect(),
                    )
                    .set_samples(vcf::variant::record_buf::Samples::new(
                        ["GT".parse().unwrap(), "DP".parse().unwrap()].into_iter().collect(),
                        vec![
                            [
                                Some(SampleValue::Genotype("0/1".parse().unwrap())),
                                Some(SampleValue::Integer(30)),
                            ]
                            .into_iter()
                            .collect(),
                        ],
                    ))
                    .build();
                writer.write_variant_record(&noodles_header, &record).unwrap();
            }
            black_box(output.len())
        });
    });

    group.finish();
}

// ---------------------------------------------------------------------------
// Group 6: Write .vcf.gz (BGZF-compressed VCF)
// seqair vs htslib — noodles doesn't have a trivial bgzf VCF writer
// ---------------------------------------------------------------------------

fn write_vcf_gz(c: &mut Criterion) {
    let mut group = c.benchmark_group("write_vcf_gz");
    group.throughput(Throughput::Elements(N_RECORDS as u64));

    let header = vcf_header();

    // seqair (BGZF-compressed VCF to memory buffer)
    group.bench_function("seqair", |b| {
        b.iter(|| {
            let mut output = Vec::with_capacity(1_000_000);
            let mut writer = VcfWriter::bgzf(&mut output, header.clone());
            writer.write_header().unwrap();
            for i in 1..=N_RECORDS {
                let alleles = Alleles::snv(Base::A, Base::T).unwrap();
                let record =
                    VcfRecordBuilder::new("chr1", seqair_types::Pos1::new(i).unwrap(), alleles)
                        .qual(30.0)
                        .filter_pass()
                        .info_integer("DP", 50)
                        .format_keys(&["GT", "DP"])
                        .add_sample({
                            use seqair_types::smallvec::smallvec;
                            smallvec![
                                SampleValue::Genotype(Genotype::unphased(0, 1)),
                                SampleValue::Integer(30),
                            ]
                        })
                        .build(&header)
                        .unwrap();
                writer.write_record(&record).unwrap();
            }
            writer.finish().unwrap();
            black_box(output.len())
        });
    });

    // htslib (writes compressed VCF to temp file)
    group.bench_function("htslib", |b| {
        b.iter(|| {
            // htslib VCF writer with compression enabled
            let tmp = tempfile::Builder::new().suffix(".vcf.gz").tempfile().unwrap();
            let hdr = htslib_header();
            let mut writer = rust_htslib::bcf::Writer::from_path(
                tmp.path(),
                &hdr,
                false,
                rust_htslib::bcf::Format::Vcf,
            )
            .unwrap();
            let rid = writer.header().name2rid(b"chr1").unwrap();
            let pass_id = writer.header().name_to_id(cstr8::cstr8!("PASS")).unwrap();
            for i in 0..N_RECORDS as i64 {
                let mut record = writer.empty_record();
                record.set_rid(Some(rid));
                record.set_pos(i);
                record.set_qual(30.0);
                record.set_alleles(&[b"A", b"T"]).unwrap();
                record.push_filter(&pass_id).unwrap();
                record.push_info_integer(cstr8::cstr8!("DP"), &[50]).unwrap();
                record
                    .push_genotypes(&[
                        rust_htslib::bcf::record::GenotypeAllele::Unphased(0),
                        rust_htslib::bcf::record::GenotypeAllele::Unphased(1),
                    ])
                    .unwrap();
                record.push_format_integer(cstr8::cstr8!("DP"), &[30]).unwrap();
                writer.write(&record).unwrap();
            }
            drop(writer);
            black_box(std::fs::metadata(tmp.path()).unwrap().len())
        });
    });

    group.finish();
}

// ---------------------------------------------------------------------------
// Group 7: Write .bcf (binary)
// seqair (VcfRecord path, encoder path) vs htslib
// ---------------------------------------------------------------------------

fn write_bcf(c: &mut Criterion) {
    let mut group = c.benchmark_group("write_bcf");
    group.throughput(Throughput::Elements(N_RECORDS as u64));

    let header = vcf_header();

    // seqair: VcfRecord path
    group.bench_function("seqair", |b| {
        b.iter(|| {
            let mut output = Vec::with_capacity(1_000_000);
            let mut writer = BcfWriter::new(&mut output, header.clone(), false);
            writer.write_header().unwrap();
            for i in 1..=N_RECORDS {
                let alleles = Alleles::snv(Base::A, Base::T).unwrap();
                let record =
                    VcfRecordBuilder::new("chr1", seqair_types::Pos1::new(i).unwrap(), alleles)
                        .qual(30.0)
                        .filter_pass()
                        .info_integer("DP", 50)
                        .format_keys(&["GT", "DP"])
                        .add_sample({
                            use seqair_types::smallvec::smallvec;
                            smallvec![
                                SampleValue::Genotype(Genotype::unphased(0, 1)),
                                SampleValue::Integer(30),
                            ]
                        })
                        .build(&header)
                        .unwrap();
                writer.write_vcf_record(&record).unwrap();
            }
            writer.finish().unwrap();
            black_box(output.len())
        });
    });

    // seqair: direct encoder path (zero-alloc, pre-resolved handles)
    group.bench_function("seqair_encoder", |b| {
        let contig = ContigHandle(0);
        let dp_info = ScalarInfoHandle::<i32> {
            dict_idx: header.string_map().get("DP").unwrap() as u32,
            _marker: PhantomData,
        };
        let gt_fmt = GtFormatHandle { dict_idx: header.string_map().get("GT").unwrap() as u32 };
        let dp_fmt = ScalarFormatHandle::<i32> {
            dict_idx: header.string_map().get("DP").unwrap() as u32,
            _marker: PhantomData,
        };

        b.iter(|| {
            let mut output = Vec::with_capacity(1_000_000);
            let mut writer = BcfWriter::new(&mut output, header.clone(), false);
            writer.write_header().unwrap();

            let alleles = Alleles::snv(Base::A, Base::T).unwrap();
            for i in 1..=N_RECORDS {
                let pos = seqair_types::Pos1::new(i).unwrap();
                let mut enc = writer.record_encoder();
                alleles.begin_record(&mut enc, contig, pos, Some(30.0)).unwrap();
                FilterHandle::PASS.encode(&mut enc);
                dp_info.encode(&mut enc, 50);
                enc.begin_samples(1);
                gt_fmt.encode(&mut enc, &Genotype::unphased(0, 1));
                dp_fmt.encode(&mut enc, 30);
                enc.emit().unwrap();
            }
            writer.finish().unwrap();
            black_box(output.len())
        });
    });

    // htslib
    group.bench_function("htslib", |b| {
        b.iter(|| {
            let tmp = write_htslib(rust_htslib::bcf::Format::Bcf);
            black_box(std::fs::metadata(tmp.path()).unwrap().len())
        });
    });

    // noodles BCF writer
    group.bench_function("noodles", |b| {
        use noodles::bcf;
        use noodles::vcf;
        use noodles::vcf::variant::io::Write as _;
        use noodles::vcf::variant::record_buf::info::field::Value as InfoValue;
        use noodles::vcf::variant::record_buf::samples::sample::Value as SampleValue;

        let noodles_header: vcf::Header = header.to_vcf_text().parse().unwrap();

        b.iter(|| {
            let mut output = Vec::with_capacity(1_000_000);
            let mut writer = bcf::io::Writer::new(&mut output);
            writer.write_header(&noodles_header).unwrap();

            for i in 1..=N_RECORDS {
                let pos = noodles::core::Position::try_from(i as usize).unwrap();
                let record = vcf::variant::RecordBuf::builder()
                    .set_reference_sequence_name("chr1")
                    .set_variant_start(pos)
                    .set_reference_bases("A")
                    .set_alternate_bases(vcf::variant::record_buf::AlternateBases::from(vec![
                        "T".to_string(),
                    ]))
                    .set_quality_score(30.0)
                    .set_filters(vcf::variant::record_buf::Filters::pass())
                    .set_info(
                        [("DP".parse().unwrap(), Some(InfoValue::Integer(50)))]
                            .into_iter()
                            .collect(),
                    )
                    .set_samples(vcf::variant::record_buf::Samples::new(
                        ["GT".parse().unwrap(), "DP".parse().unwrap()].into_iter().collect(),
                        vec![
                            [
                                Some(SampleValue::Genotype("0/1".parse().unwrap())),
                                Some(SampleValue::Integer(30)),
                            ]
                            .into_iter()
                            .collect(),
                        ],
                    ))
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
    fasta_fetch,
    pileup_e2e,
    write_vcf_text,
    write_vcf_gz,
    write_bcf,
);
criterion_main!(benches);
