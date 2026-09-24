//! Profiling harness for text parsing and write paths (not a bench).
#![allow(
    clippy::panic,
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::print_stdout,
    clippy::type_complexity,
    reason = "profiling harness"
)]
use core::range::RangeInclusive;
use seqair::vcf::record_encoder::{FilterFieldDef, FormatFieldDef, InfoFieldDef};
use seqair::vcf::{
    Alleles, ContigDef, ContigId, FilterId, FormatGt, FormatInt, Genotype, InfoFloat, InfoFloats,
    InfoInt, InfoInts, Number, OutputFormat, Ready, ValueType, VcfHeader, Writer,
};
use seqair_types::{Base, Pos0, Pos1};
use std::fs::File;
use std::hint::black_box;
use std::io::BufWriter;
use std::sync::Arc;
const N_SAMPLES: usize = 5;
fn dev_null() -> BufWriter<File> {
    BufWriter::new(File::create("/dev/null").unwrap())
}

struct SeqairSetup {
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

fn seqair_setup() -> SeqairSetup {
    let mut b = VcfHeader::builder();
    let contig = b.register_contig("chr1", ContigDef { length: Some(250_000_000) }).unwrap();

    let mut b = b.filters();
    let lowqual_filter = b.register_filter(&FilterFieldDef::new("LowQual", "Low quality")).unwrap();

    let mut b = b.infos();
    let dp_info = b
        .register_info(&InfoFieldDef::new(
            "DP",
            Number::Count(1),
            ValueType::Integer,
            "Total read depth",
        ))
        .unwrap();
    let an_info = b
        .register_info(&InfoFieldDef::new(
            "AN",
            Number::Count(1),
            ValueType::Integer,
            "Total allele count in genotypes",
        ))
        .unwrap();
    let ac_info = b
        .register_info(&InfoFieldDef::new(
            "AC",
            Number::AlternateBases,
            ValueType::Integer,
            "Allele count in genotypes, for each ALT",
        ))
        .unwrap();
    let af_info = b
        .register_info(&InfoFieldDef::new(
            "AF",
            Number::AlternateBases,
            ValueType::Float,
            "Allele frequency from GT counts, for each ALT",
        ))
        .unwrap();
    let mq_info = b
        .register_info(&InfoFieldDef::new(
            "MQ",
            Number::Count(1),
            ValueType::Float,
            "RMS mapping quality",
        ))
        .unwrap();
    let qd_info = b
        .register_info(&InfoFieldDef::new(
            "QD",
            Number::Count(1),
            ValueType::Float,
            "Variant confidence divided by unfiltered depth",
        ))
        .unwrap();
    let fs_info = b
        .register_info(&InfoFieldDef::new(
            "FS",
            Number::Count(1),
            ValueType::Float,
            "Phred-scaled p-value using Fisher's exact test for strand bias",
        ))
        .unwrap();

    let mut b = b.formats();
    let gt_fmt = b
        .register_format(&FormatFieldDef::new(
            "GT",
            Number::Count(1),
            ValueType::String,
            "Genotype",
        ))
        .unwrap();
    let dp_fmt = b
        .register_format(&FormatFieldDef::new(
            "DP",
            Number::Count(1),
            ValueType::Integer,
            "Approximate read depth",
        ))
        .unwrap();
    let gq_fmt = b
        .register_format(&FormatFieldDef::new(
            "GQ",
            Number::Count(1),
            ValueType::Integer,
            "Genotype quality",
        ))
        .unwrap();

    let mut b = b.samples();
    for s in 1..=N_SAMPLES {
        b.add_sample(format!("SAMPLE{s:02}")).unwrap();
    }
    let header = Arc::new(b.build().unwrap());

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

    SeqairSetup {
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

/// Encode one VCF record using the seqair typestate encoder.
///
/// When `minimal` is true, only DP (INFO) and GT (FORMAT) are written —
/// the remaining 6 INFO and 2 FORMAT fields are skipped.
fn seqair_encode_record<W: std::io::Write>(
    writer: &mut Writer<W, Ready>,
    s: &SeqairSetup,
    i: u32,
    minimal: bool,
) {
    let alleles = &s.alleles_bank[i as usize % s.alleles_bank.len()];
    let n_alt = alleles.n_allele() - 1;
    let an = (N_SAMPLES * 2) as i32;
    let pos = Pos1::new(i * 237 + 1).unwrap();

    let enc =
        writer.begin_record(&s.contig, pos, alleles, Some(10.0 + (i % 50) as f32 * 2.0)).unwrap();
    let mut enc =
        if i % 10 == 7 { enc.filter_fail([&s.lowqual_filter]) } else { enc.filter_pass() };

    let dp = 30 + (i % 70) as i32;
    s.dp_info.encode(&mut enc, dp);

    if !minimal {
        s.an_info.encode(&mut enc, an);
        let ac: Vec<i32> = (0..n_alt).map(|k| 1 + (i as i32 + k as i32) % 9).collect();
        let af: Vec<f32> = ac.iter().map(|&a| a as f32 / an as f32).collect();
        s.ac_info.encode(&mut enc, &ac);
        s.af_info.encode(&mut enc, &af);
        s.mq_info.encode(&mut enc, 55.0 + (i % 20) as f32 * 0.5);
        s.qd_info.encode(&mut enc, 5.0 + (i % 25) as f32 * 0.8);
        s.fs_info.encode(&mut enc, (i % 30) as f32 * 0.5);
    }

    let mut enc = enc.begin_samples();
    let gts: Vec<Genotype> = (0..N_SAMPLES)
        .map(|si| match (si + i as usize) % 6 {
            0 => Genotype::unphased(0, 0),
            1 => Genotype::unphased(0, 1),
            2 => Genotype::unphased(1, 1),
            3 => Genotype::phased_diploid(0, 1),
            4 => Genotype::missing_diploid(),
            _ => Genotype::unphased(0, 1),
        })
        .collect();
    s.gt_fmt.encode(&mut enc, &gts).unwrap();

    if !minimal {
        let dp_vals: Vec<i32> = (0..N_SAMPLES as i32).map(|si| dp / 2 + si % 20).collect();
        let gq_vals: Vec<i32> = (0..N_SAMPLES as i32).map(|si| 30 + (i as i32 + si) % 60).collect();
        s.dp_fmt.encode(&mut enc, &dp_vals).unwrap();
        s.gq_fmt.encode(&mut enc, &gq_vals).unwrap();
    }

    enc.emit().unwrap();
}

fn vcf(fmt: OutputFormat, iters: u32, n: u32) {
    let s = seqair_setup();
    for _ in 0..iters {
        let w = Writer::new(dev_null(), fmt);
        let mut w = w.write_header(&s.header).unwrap();
        for i in 0..n {
            seqair_encode_record(&mut w, &s, i, false);
        }
        w.finish().unwrap();
    }
}

fn bam(iters: u32) {
    use seqair::bam::BamHeader;
    use seqair::bam::owned_record::OwnedBamRecord;
    use seqair::bam::writer::BamWriterBuilder;
    let path =
        std::path::Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.bam"));
    let mut reader = seqair::bam::IndexedBamReader::open(path).unwrap();
    let header = BamHeader::from_template(reader.header());
    let tid = reader.header().tid("chr19").unwrap();
    let mut store = seqair::bam::RecordStore::new();
    let span =
        RangeInclusive { start: Pos0::new(0).unwrap(), last: Pos0::new(61_000_000).unwrap() };
    reader.fetch_into(tid, span, &mut store).unwrap();
    let recs: Vec<OwnedBamRecord> = store
        .indices()
        .map(|idx| {
            let r = store.record(idx).unwrap();
            OwnedBamRecord {
                ref_id: tid as i32,
                pos: Some(r.pos),
                mapq: r.mapq,
                flags: r.flags,
                next_ref_id: -1,
                next_pos: None,
                template_len: 0,
                qname: r.qname().to_vec(),
                cigar: r.cigar().to_vec(),
                seq: r.seq().to_vec(),
                qual: r.qual().to_vec(),
                aux: seqair::bam::AuxData::from_bytes(r.aux().to_vec()),
            }
        })
        .collect();
    println!("{} records", recs.len());
    let out_buf = std::env::temp_dir().join("prof_out.bam");
    let out = out_buf.as_path();
    for _ in 0..iters {
        let mut writer = BamWriterBuilder::to_path(out, &header).write_index(true).build().unwrap();
        for rec in &recs {
            writer.write(rec).unwrap();
        }
        black_box(writer.finish().unwrap());
    }
}

fn bam_serialize_only(iters: u32) {
    // only to_bam_bytes, no compression
    use seqair::bam::owned_record::OwnedBamRecord;
    let path =
        std::path::Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.bam"));
    let mut reader = seqair::bam::IndexedBamReader::open(path).unwrap();
    let tid = reader.header().tid("chr19").unwrap();
    let mut store = seqair::bam::RecordStore::new();
    let span =
        RangeInclusive { start: Pos0::new(0).unwrap(), last: Pos0::new(61_000_000).unwrap() };
    reader.fetch_into(tid, span, &mut store).unwrap();
    let recs: Vec<OwnedBamRecord> = store
        .indices()
        .map(|idx| {
            let r = store.record(idx).unwrap();
            OwnedBamRecord {
                ref_id: tid as i32,
                pos: Some(r.pos),
                mapq: r.mapq,
                flags: r.flags,
                next_ref_id: -1,
                next_pos: None,
                template_len: 0,
                qname: r.qname().to_vec(),
                cigar: r.cigar().to_vec(),
                seq: r.seq().to_vec(),
                qual: r.qual().to_vec(),
                aux: seqair::bam::AuxData::from_bytes(r.aux().to_vec()),
            }
        })
        .collect();
    let mut buf = Vec::with_capacity(1 << 16);
    let mut total = 0usize;
    for _ in 0..iters {
        for rec in &recs {
            buf.clear();
            rec.to_bam_bytes(&mut buf).unwrap();
            total += buf.len();
        }
    }
    println!("{total}");
}

fn sam(iters: u32, path: &str) {
    let mut r = seqair::sam::reader::IndexedSamReader::open(std::path::Path::new(path)).unwrap();
    let tid = r.header().tid("chr19").unwrap();
    let mut store = seqair::bam::RecordStore::new();
    let span =
        RangeInclusive { start: Pos0::new(0).unwrap(), last: Pos0::new(61_000_000).unwrap() };
    let mut n = 0;
    for _ in 0..iters {
        n += r.fetch_into(tid, span, &mut store).unwrap();
    }
    println!("{n}");
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let iters: u32 = args.get(2).map_or(10, |s| s.parse().unwrap());
    match args[1].as_str() {
        "vcf" => vcf(OutputFormat::Vcf, iters, 10_000),
        "vcfgz" => vcf(OutputFormat::VcfGz, iters, 10_000),
        "bcf" => vcf(OutputFormat::Bcf, iters, 10_000),
        "bam" => bam(iters),
        "bamser" => bam_serialize_only(iters),
        "sam" => sam(iters, &args[3]),
        _ => panic!("mode"),
    }
}
