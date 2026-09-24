//! Parallel BGZF compression is an implementation detail: the BAM and
//! VCF/BCF writers must produce the same bytes — data file and co-produced
//! index — with `compression_threads(n)` as without.
//!
//! The single-threaded writer is the oracle. Its own correctness is covered by
//! the samtools/bcftools/noodles round trips; here the property is equality,
//! over record streams long enough to span many BGZF blocks.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    reason = "test code"
)]

use hegel::prelude::*;
use seqair::bam::aux_data::AuxData;
use seqair::bam::cigar::{CigarOp, CigarOpType};
use seqair::bam::header::BamHeader;
use seqair::bam::owned_record::OwnedBamRecord;
use seqair::bam::writer::BamWriterBuilder;
use seqair::vcf::record_encoder::{FormatFieldDef, InfoFieldDef};
use seqair::vcf::{
    Alleles, ContigDef, FormatInt, Genotype, InfoInts, Number, OutputFormat, ValueType, VcfHeader,
    Writer,
};
use seqair::vcf::{FormatGt, InfoFloat};
use seqair_types::{BamFlags, Base, BaseQuality, Pos0, Pos1};

const ACGT: [Base; 4] = [Base::A, Base::C, Base::G, Base::T];

/// A tiny xorshift so record content is cheap to produce in bulk: hegel picks
/// the stream's shape and the seed, not every base.
struct Noise(u64);

impl Noise {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
}

// ── BAM ─────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
struct BamStream {
    /// Per record: position step from the previous record, read length, and
    /// whether it is a placed-unmapped record.
    shape: Vec<(u32, usize, bool)>,
    /// Index of the first record on the second contig.
    switch_at: usize,
    seed: u64,
}

#[hegel::composite]
fn arb_bam_stream_inner(tc: &TestCase) -> BamStream {
    let n = tc.draw(gs::integers::<usize>().max_value(2_500));
    let shape = (0..n)
        .map(|_| {
            (
                tc.draw_silent(gs::integers::<u32>().max_value(200)),
                tc.draw_silent(gs::integers::<usize>().min_value(1).max_value(400)),
                tc.draw_silent(gs::integers::<u8>().max_value(30)) == 0,
            )
        })
        .collect();
    let switch_at = tc.draw(gs::integers::<usize>().max_value(n));
    let seed = tc.draw(gs::integers::<u64>().min_value(1));
    BamStream { shape, switch_at, seed }
}

fn arb_bam_stream() -> impl PrintableGenerator<BamStream> {
    arb_bam_stream_inner().print_as_debug()
}

fn bam_header() -> BamHeader {
    BamHeader::from_sam_text(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:10000000\n@SQ\tSN:chr2\tLN:10000000\n",
    )
    .unwrap()
}

fn bam_records(stream: &BamStream) -> Vec<OwnedBamRecord> {
    let mut noise = Noise(stream.seed);
    let mut pos = 0u32;
    stream
        .shape
        .iter()
        .enumerate()
        .map(|(i, &(step, len, unmapped))| {
            let tid = i32::from(i >= stream.switch_at);
            if i == stream.switch_at {
                pos = 0;
            }
            pos += step;
            let seq: Vec<Base> = (0..len).map(|_| ACGT[(noise.next() % 4) as usize]).collect();
            let qual: Vec<BaseQuality> =
                (0..len).map(|_| BaseQuality::from_byte(20 + (noise.next() % 20) as u8)).collect();
            let mut aux = AuxData::new();
            aux.set_string(*b"MM", b"C+m,0,3,1;");
            let ml: Vec<u8> = (0..(len % 7)).map(|_| noise.next() as u8).collect();
            aux.set_array_u8(*b"ML", &ml).unwrap();
            let flags = if unmapped { BamFlags::from(0x4) } else { BamFlags::empty() };
            let cigar = if unmapped {
                Vec::new()
            } else {
                vec![CigarOp::new(CigarOpType::Match, len as u32)]
            };
            OwnedBamRecord::builder(tid, Some(Pos0::new(pos).unwrap()), format!("r{i}").into())
                .flags(flags)
                .mapq(60)
                .cigar(cigar)
                .seq(seq)
                .qual(qual)
                .aux(aux)
                .build()
                .unwrap()
        })
        .collect()
}

/// BAM and BAI bytes from writing `records` with `threads` compression threads.
fn write_bam(records: &[OwnedBamRecord], level: i32, threads: usize) -> (Vec<u8>, Vec<u8>) {
    let header = bam_header();
    let file = tempfile::NamedTempFile::new().unwrap();
    let mut writer = BamWriterBuilder::to_path(file.path(), &header)
        .compression_level(level)
        .compression_threads(threads)
        .write_index(true)
        .build()
        .unwrap();
    for record in records {
        writer.write(record).unwrap();
    }
    let (_, index) = writer.finish().unwrap();
    let mut bai = Vec::new();
    index.unwrap().write_bai(&mut bai, header.target_count()).unwrap();
    (std::fs::read(file.path()).unwrap(), bai)
}

// r[verify bam_writer.multithreaded_compression]
// r[verify bgzf.writer.parallel.identical_output]
// r[verify bgzf.writer.parallel.index_offsets]
#[hegel::test(test_cases = 30)]
fn parallel_bam_matches_serial(tc: TestCase) {
    let stream = tc.draw(arb_bam_stream());
    let threads = tc.draw(gs::integers::<usize>().min_value(1).max_value(4));
    let level = tc.draw(gs::sampled_from(&[0, 1, 6]));
    let records = bam_records(&stream);

    let (bam, bai) = write_bam(&records, level, 0);
    let (par_bam, par_bai) = write_bam(&records, level, threads);
    assert!(par_bam == bam, "BAM bytes differ with {threads} threads at level {level}");
    assert!(par_bai == bai, "BAI bytes differ with {threads} threads at level {level}");
}

// ── VCF.gz / BCF ────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
struct VcfStream {
    /// Number of records; their content comes from `seed`.
    records: usize,
    seed: u64,
}

#[hegel::composite]
fn arb_vcf_stream_inner(tc: &TestCase) -> VcfStream {
    // Up to ~15 BGZF blocks of BCF.
    let records = tc.draw(gs::integers::<usize>().max_value(40_000));
    VcfStream { records, seed: tc.draw(gs::integers::<u64>().min_value(1)) }
}

fn arb_vcf_stream() -> impl PrintableGenerator<VcfStream> {
    arb_vcf_stream_inner().print_as_debug()
}

struct VcfSetup {
    header: VcfHeader,
    contig: seqair::vcf::ContigId,
    ad: InfoInts,
    mq: InfoFloat,
    gt: FormatGt,
    dp: FormatInt,
}

fn vcf_setup() -> VcfSetup {
    let mut b = VcfHeader::builder();
    let contig = b.register_contig("chr1", ContigDef { length: Some(250_000_000) }).unwrap();
    let mut b = b.infos();
    let ad = b
        .register_info(&InfoFieldDef::new(
            "AD",
            Number::ReferenceAlternateBases,
            ValueType::Integer,
            "Depth per allele",
        ))
        .unwrap();
    let mq = b
        .register_info(&InfoFieldDef::new("MQ", Number::Count(1), ValueType::Float, "RMS MAPQ"))
        .unwrap();
    let mut b = b.formats();
    let gt = b
        .register_format(&FormatFieldDef::new("GT", Number::Count(1), ValueType::String, "GT"))
        .unwrap();
    let dp = b
        .register_format(&FormatFieldDef::new("DP", Number::Count(1), ValueType::Integer, "DP"))
        .unwrap();
    let mut b = b.samples();
    b.add_sample("sample").unwrap();
    VcfSetup { header: b.build().unwrap(), contig, ad, mq, gt, dp }
}

/// Data and CSI bytes from writing the stream with `threads` compression threads.
fn write_vcf(
    setup: &VcfSetup,
    stream: &VcfStream,
    format: OutputFormat,
    threads: usize,
) -> (Vec<u8>, Vec<u8>) {
    let mut noise = Noise(stream.seed);
    let writer = Writer::new(Vec::new(), format).compression_threads(threads).unwrap();
    let mut writer = writer.write_header(&setup.header).unwrap();
    let mut pos = 0u32;
    for _ in 0..stream.records {
        pos += 1 + (noise.next() % 300) as u32;
        let r = ACGT[(noise.next() % 4) as usize];
        let a = ACGT[(noise.next() % 4) as usize];
        let alleles = if a == r { Alleles::reference(r) } else { Alleles::snv(r, a).unwrap() };
        let n_allele = alleles.n_allele();
        let qual = Some((noise.next() % 1000) as f32 / 7.0);
        let enc = writer.begin_record(&setup.contig, Pos1::new(pos).unwrap(), &alleles, qual);
        let mut enc = enc.unwrap().filter_pass();
        let ad: Vec<i32> = (0..n_allele).map(|_| (noise.next() % 400) as i32).collect();
        setup.ad.encode(&mut enc, &ad);
        setup.mq.encode(&mut enc, (noise.next() % 6000) as f32 / 100.0);
        let mut enc = enc.begin_samples();
        setup.gt.encode(&mut enc, &[Genotype::unphased(0, u16::from(n_allele > 1))]).unwrap();
        setup.dp.encode(&mut enc, &[(noise.next() % 40_000) as i32]).unwrap();
        enc.emit().unwrap();
    }
    let (data, index) = writer.finish().unwrap();
    let mut csi = Vec::new();
    index.unwrap().write(&mut csi).unwrap();
    (data, csi)
}

// r[verify record_encoder.compression_threads]
// r[verify bgzf.writer.parallel.identical_output]
// r[verify bgzf.writer.parallel.index_offsets]
#[hegel::test(test_cases = 30)]
fn parallel_bcf_and_vcf_gz_match_serial(tc: TestCase) {
    let stream = tc.draw(arb_vcf_stream());
    let threads = tc.draw(gs::integers::<usize>().min_value(1).max_value(4));
    let format =
        tc.draw(gs::sampled_from(&[OutputFormat::Bcf, OutputFormat::VcfGz]).print_as_debug());
    let setup = vcf_setup();

    let (data, csi) = write_vcf(&setup, &stream, format, 0);
    let (par_data, par_csi) = write_vcf(&setup, &stream, format, threads);
    assert!(par_data == data, "{format:?} bytes differ with {threads} threads");
    assert!(par_csi == csi, "{format:?} index bytes differ with {threads} threads");
}

// r[verify record_encoder.compression_threads]
/// Plain VCF is not compressed; asking for threads changes nothing.
#[test]
fn plain_vcf_ignores_compression_threads() {
    let setup = vcf_setup();
    let stream = VcfStream { records: 50, seed: 7 };
    let mut outputs = Vec::new();
    for threads in [0, 3] {
        let writer = Writer::new(Vec::new(), OutputFormat::Vcf).compression_threads(threads);
        let mut writer = writer.unwrap().write_header(&setup.header).unwrap();
        for i in 0..stream.records {
            let alleles = Alleles::reference(Base::A);
            let pos = Pos1::new(u32::try_from(i).unwrap() + 1).unwrap();
            let enc = writer.begin_record(&setup.contig, pos, &alleles, None).unwrap();
            enc.filter_pass().emit().unwrap();
        }
        outputs.push(writer.finish().unwrap().0);
    }
    assert_eq!(outputs[0], outputs[1]);
}
