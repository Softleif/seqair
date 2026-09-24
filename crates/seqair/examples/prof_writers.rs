//! Profiling harness for the write paths, shaped like rastair's use of them
//! (not a bench).
//!
//! Every mode loads its input into memory first and then times only the
//! writing, `iters` times, printing wall time and process CPU time
//! (user + sys, all threads) per iteration.
//!
//! ```text
//! # BCF/VCF: replay the records of a rastair BCF (default field set) through
//! # seqair's encoder. `fmt` is bcf | vcf | vcfgz.
//! prof_writers bcf-replay <in.bcf> <fmt> <iters> [out] [compression threads]
//!
//! # BAM rewrite: preload a region, add MM/ML-like tags, write.
//! # level -1: serialize only (to_bam_bytes, no BGZF)
//! prof_writers bam-seqair <in.bam> <contig> <start> <end> <level> <index 0|1> <iters> [out] [threads]
//! prof_writers bam-htslib <in.bam> <contig> <start> <end> <level> <threads> <iters> [out]
//! ```
#![allow(
    clippy::panic,
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss,
    clippy::cast_precision_loss,
    clippy::print_stdout,
    clippy::too_many_lines,
    clippy::type_complexity,
    reason = "profiling harness"
)]

use core::range::RangeInclusive;
use seqair::vcf::record_encoder::{Arr, FilterFieldDef, FormatFieldDef, InfoFieldDef, Scalar, Str};
use seqair::vcf::{
    Alleles, ContigDef, ContigId, FilterId, FormatFloat, FormatFloats, FormatGt, FormatInt,
    FormatInts, Genotype, InfoFlag, InfoFloat, InfoFloats, InfoInt, InfoInts, InfoString, Number,
    OutputFormat, Ready, ValueType, VcfHeader, Writer,
};
use seqair_types::{Base, Pos0, Pos1, SmallVec, SmolStr};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::Instant;

// ── timing ──────────────────────────────────────────────────────────────

fn cpu_seconds() -> f64 {
    // SAFETY: getrusage writes into the zeroed struct we own.
    let mut ru: libc::rusage = unsafe { std::mem::zeroed() };
    // SAFETY: valid pointer to an rusage struct.
    unsafe { libc::getrusage(libc::RUSAGE_SELF, &raw mut ru) };
    // `tv_usec` is i32 on macOS and i64 on Linux.
    #[allow(clippy::cast_lossless, clippy::unnecessary_cast, reason = "platform-dependent width")]
    let tv = |t: libc::timeval| t.tv_sec as f64 + t.tv_usec as f64 * 1e-6;
    tv(ru.ru_utime) + tv(ru.ru_stime)
}

/// User-space cycles and instructions of this process and every thread it
/// starts afterwards (Linux `perf_event_open` with `inherit`); `None` elsewhere
/// or when perf events are not permitted.
#[cfg(target_os = "linux")]
struct Counters([i32; 2]);

#[cfg(target_os = "linux")]
impl Counters {
    fn open() -> Option<Self> {
        let mut fds = [-1; 2];
        for (config, fd) in fds.iter_mut().enumerate() {
            // perf_event_attr: type, size, config, ..., flags at byte 40.
            let mut attr = [0u64; 16];
            attr[0] = 128 << 32; // type = PERF_TYPE_HARDWARE (0), size = 128
            attr[1] = config as u64; // 0 = cycles, 1 = instructions
            attr[5] = (1 << 1) | (1 << 5) | (1 << 6); // inherit, exclude_kernel, exclude_hv
            // SAFETY: attr is a zero-padded perf_event_attr of the declared size.
            let r =
                unsafe { libc::syscall(libc::SYS_perf_event_open, attr.as_ptr(), 0, -1, -1, 0) };
            if r < 0 {
                return None;
            }
            *fd = r as i32;
        }
        Some(Self(fds))
    }

    fn read(&self) -> [u64; 2] {
        self.0.map(|fd| {
            let mut v = 0u64;
            // SAFETY: reading 8 bytes into a u64 from a perf event fd.
            unsafe { libc::read(fd, (&raw mut v).cast(), 8) };
            v
        })
    }
}

fn timed(label: &str, iters: u32, mut f: impl FnMut()) {
    let mut walls = Vec::new();
    let mut cpus = Vec::new();
    #[cfg(target_os = "linux")]
    let counters = Counters::open();
    #[cfg(target_os = "linux")]
    let mut cycles = Vec::new();
    for _ in 0..iters {
        #[cfg(target_os = "linux")]
        let k0 = counters.as_ref().map(Counters::read);
        let c0 = cpu_seconds();
        let t0 = Instant::now();
        f();
        walls.push(t0.elapsed().as_secs_f64());
        cpus.push(cpu_seconds() - c0);
        #[cfg(target_os = "linux")]
        if let (Some(k0), Some(c)) = (k0, counters.as_ref()) {
            let k1 = c.read();
            cycles.push((k1[0] - k0[0], k1[1] - k0[1]));
        }
    }
    let min = |v: &[f64]| v.iter().copied().fold(f64::INFINITY, f64::min);
    println!(
        "{label}: wall min {:.3} s  cpu min {:.3} s  (walls {:?})",
        min(&walls),
        min(&cpus),
        walls.iter().map(|w| (w * 1000.0).round() / 1000.0).collect::<Vec<_>>()
    );
    #[cfg(target_os = "linux")]
    if let Some(&(c, i)) = cycles.iter().min() {
        println!(
            "{label}: user cycles min {:.3} G  instructions {:.3} G",
            c as f64 / 1e9,
            i as f64 / 1e9
        );
    }
}

fn out_path(arg: Option<&String>, default: &str) -> PathBuf {
    arg.map_or_else(|| std::env::temp_dir().join(default), PathBuf::from)
}

// ── BCF replay ──────────────────────────────────────────────────────────

/// rastair's schema (`src/vcf/schema.rs`), every field registered, same
/// handle kinds, same header order.
struct Schema {
    ad: InfoInts,
    bq: InfoFloat,
    dp: InfoInt,
    mq: InfoFloat,
    mq0: InfoInt,
    ns: InfoInt,
    as_sb_ot: InfoInts,
    as_sb_ob: InfoInts,
    sc5: InfoString,
    af: InfoFloats,
    abq: InfoFloats,
    amq: InfoFloats,
    as_ss_bq_ot: InfoFloats,
    as_ss_bq_ob: InfoFloats,
    as_ss_mq_ot: InfoFloats,
    as_ss_mq_ob: InfoFloats,
    pir: InfoFloats,
    ent100: InfoFloat,
    nab: InfoFloats,
    noi: InfoFloats,
    m5mc_strands: InfoInts,
    cpg: InfoFlag,
    cpgnovo: InfoFlag,
    gt: FormatGt,
    gl: FormatFloat,
    gc: FormatFloat,
    fdp: FormatInt,
    m5mc: FormatFloats,
    dpm5mc: FormatInts,
    adm5mc: FormatInts,
    ml: FormatFloats,
    filters: Vec<(SmolStr, FilterId)>,
    contigs: Vec<(SmolStr, ContigId)>,
}

const FILTERS: &[(&str, &str)] = &[
    ("lowDp", "Low read depth"),
    ("dnCpG_lowDp", "Low read depth for de-novo CpG candidate"),
    ("dnCpG_bq", "Low base quality for de-novo CpG candidate"),
    ("dnCpG_mapq", "Low mapping quality for de-novo CpG candidate"),
    ("dnCpG_vaf", "Low variant allele frequency for de-novo CpG candidate"),
    ("dnCpG_adj", "Included as adjacent position for de-novo CpG candidate"),
    ("m_vaf", "Low variant allele frequency for methylation candidate"),
    ("m_bq_ratio", "Low quality ratio for methylation candidate"),
    ("m_pos", "Alt allele evidence from read edges for methylation candidate"),
    ("m_highDp", "Excessive coverage for methylation candidate"),
    ("pre_ml", "Low amount of usable evidence, skipping ML"),
    ("low_ml_score", "Machine Learning module prediction below threshold"),
    ("indel_strand", "Indel allele supported on only one bisulfite strand"),
];

fn build_schema(contigs: &[(String, u64)]) -> (VcfHeader, Schema) {
    let mut b = VcfHeader::builder();
    let contig_ids = contigs
        .iter()
        .map(|(name, len)| {
            let id = b.register_contig(name.as_str(), ContigDef { length: Some(*len) }).unwrap();
            (SmolStr::new(name), id)
        })
        .collect();
    let mut b = b.filters();
    let filters = FILTERS
        .iter()
        .map(|(n, d)| (SmolStr::new(n), b.register_filter(&FilterFieldDef::new(n, d)).unwrap()))
        .collect();
    let mut b = b.infos();
    let r = Number::ReferenceAlternateBases;
    let one = Number::Count(1);
    let (int, float) = (ValueType::Integer, ValueType::Float);
    macro_rules! info {
        ($kind:ty, $name:expr, $num:expr, $ty:expr) => {
            b.register_info(&InfoFieldDef::<$kind>::new($name, $num, $ty, "d")).unwrap()
        };
    }
    let ad = info!(Arr<i32>, "AD", r, int);
    let bq = info!(Scalar<f32>, "BQ", one, float);
    let dp = info!(Scalar<i32>, "DP", one, int);
    let mq = info!(Scalar<f32>, "MQ", one, float);
    let mq0 = info!(Scalar<i32>, "MQ0", one, int);
    let ns = info!(Scalar<i32>, "NS", one, int);
    let as_sb_ot = info!(Arr<i32>, "AS_SB_OT", r, int);
    let as_sb_ob = info!(Arr<i32>, "AS_SB_OB", r, int);
    let sc5 = info!(Str, "SC5", one, ValueType::String);
    let af = info!(Arr<f32>, "AF", Number::AlternateBases, float);
    let abq = info!(Arr<f32>, "ABQ", r, float);
    let amq = info!(Arr<f32>, "AMQ", r, float);
    let as_ss_bq_ot = info!(Arr<f32>, "AS_SS_BQ_OT", r, float);
    let as_ss_bq_ob = info!(Arr<f32>, "AS_SS_BQ_OB", r, float);
    let as_ss_mq_ot = info!(Arr<f32>, "AS_SS_MQ_OT", r, float);
    let as_ss_mq_ob = info!(Arr<f32>, "AS_SS_MQ_OB", r, float);
    let pir = info!(Arr<f32>, "PIR", r, float);
    let ent100 = info!(Scalar<f32>, "ENT100", one, float);
    let nab = info!(Arr<f32>, "NAB", r, float);
    let noi = info!(Arr<f32>, "NOI", r, float);
    let m5mc_strands = info!(Arr<i32>, "M5mC_Strands", Number::Count(4), int);
    let cpg = info!(seqair::vcf::Flag, "CPG", Number::Count(0), ValueType::Flag);
    let cpgnovo = info!(seqair::vcf::Flag, "CPGnovo", Number::Count(0), ValueType::Flag);
    let mut b = b.formats();
    macro_rules! format {
        ($kind:ty, $name:expr, $num:expr, $ty:expr) => {
            b.register_format(&FormatFieldDef::<$kind>::new($name, $num, $ty, "d")).unwrap()
        };
    }
    let gt = format!(seqair::vcf::Gt, "GT", one, ValueType::String);
    let gl = format!(Scalar<f32>, "GL", Number::Genotypes, float);
    let gc = format!(Scalar<f32>, "GC", Number::Genotypes, float);
    let fdp = format!(Scalar<i32>, "DP", one, int);
    let m5mc = format!(Arr<f32>, "M5mC", Number::Unknown, float);
    let dpm5mc = format!(Arr<i32>, "DPM5mC", Number::Unknown, int);
    let adm5mc = format!(Arr<i32>, "ADM5mC", Number::Unknown, int);
    let ml = format!(Arr<f32>, "ML", Number::AlternateBases, float);
    let mut b = b.samples();
    b.add_sample("sample").unwrap();
    let header = b.build().unwrap();
    let schema = Schema {
        ad,
        bq,
        dp,
        mq,
        mq0,
        ns,
        as_sb_ot,
        as_sb_ob,
        sc5,
        af,
        abq,
        amq,
        as_ss_bq_ot,
        as_ss_bq_ob,
        as_ss_mq_ot,
        as_ss_mq_ob,
        pir,
        ent100,
        nab,
        noi,
        m5mc_strands,
        cpg,
        cpgnovo,
        gt,
        gl,
        gc,
        fdp,
        m5mc,
        dpm5mc,
        adm5mc,
        ml,
        filters,
        contigs: contig_ids,
    };
    (header, schema)
}

#[derive(Clone)]
enum Val {
    Int(i32),
    Ints(SmallVec<i32, 4>),
    Float(f32),
    Floats(SmallVec<f32, 3>),
    Flag,
    Str(SmolStr),
    Gt(Genotype),
}

/// A decoded record, in the order rastair encodes its fields.
struct Rec {
    contig: u16,
    pos: Pos1,
    alleles: Alleles,
    qual: Option<f32>,
    filters: SmallVec<u8, 2>,
    info: SmallVec<(u8, Val), 8>,
    format: SmallVec<(u8, Val), 8>,
}

const INFO_IDS: &[&str] = &[
    "AD",
    "BQ",
    "DP",
    "MQ",
    "MQ0",
    "NS",
    "AS_SB_OT",
    "AS_SB_OB",
    "SC5",
    "AF",
    "ABQ",
    "AMQ",
    "AS_SS_BQ_OT",
    "AS_SS_BQ_OB",
    "AS_SS_MQ_OT",
    "AS_SS_MQ_OB",
    "PIR",
    "ENT100",
    "NAB",
    "NOI",
    "M5mC_Strands",
    "CPG",
    "CPGnovo",
];
const FORMAT_IDS: &[&str] = &["GT", "GL", "GC", "DP", "M5mC", "DPM5mC", "ADM5mC", "ML"];

fn base(c: u8) -> Base {
    Base::from(c)
}

fn alleles_of(r: &str, alts: &[String]) -> Alleles {
    let rb = r.as_bytes();
    if alts.is_empty() {
        return Alleles::reference(base(rb[0]));
    }
    if rb.len() == 1 && alts.iter().all(|a| a.len() == 1) {
        let bases: SmallVec<Base, 2> = alts.iter().map(|a| base(a.as_bytes()[0])).collect();
        return Alleles::snv_multi(base(rb[0]), &bases).unwrap();
    }
    if alts.len() == 1 {
        let a = alts[0].as_bytes();
        if rb.len() == 1 && a.len() > 1 && a[0] == rb[0] {
            let ins: Vec<Base> = a[1..].iter().map(|&c| base(c)).collect();
            return Alleles::insertion(base(rb[0]), &ins).unwrap();
        }
        if a.len() == 1 && rb.len() > 1 && a[0] == rb[0] {
            let del: Vec<Base> = rb[1..].iter().map(|&c| base(c)).collect();
            return Alleles::deletion(base(rb[0]), &del).unwrap();
        }
    }
    Alleles::complex(SmolStr::new(r), alts.iter().map(SmolStr::new).collect())
}

fn load_bcf(path: &Path) -> (Vec<(String, u64)>, Vec<Rec>) {
    use noodles::vcf::variant::record_buf::info::field::Value as IV;
    use noodles::vcf::variant::record_buf::info::field::value::Array as IA;
    use noodles::vcf::variant::record_buf::samples::sample::Value as SV;
    use noodles::vcf::variant::record_buf::samples::sample::value::Array as SA;

    let mut reader = noodles::bcf::io::reader::Builder::default().build_from_path(path).unwrap();
    let header = reader.read_header().unwrap();
    let contigs: Vec<(String, u64)> = header
        .contigs()
        .iter()
        .map(|(name, c)| (name.clone(), c.length().unwrap_or(0) as u64))
        .collect();
    let contig_index: std::collections::HashMap<String, u16> =
        contigs.iter().enumerate().map(|(i, (n, _))| (n.clone(), i as u16)).collect();
    let filter_index: std::collections::HashMap<&str, u8> =
        FILTERS.iter().enumerate().map(|(i, (n, _))| (*n, i as u8)).collect();

    let mut out = Vec::new();
    for result in reader.record_bufs(&header) {
        let rec = result.unwrap();
        let contig = contig_index[rec.reference_sequence_name()];
        let pos = Pos1::new(rec.variant_start().unwrap().get() as u32).unwrap();
        let alts: Vec<String> = rec.alternate_bases().as_ref().to_vec();
        let alleles = alleles_of(rec.reference_bases(), &alts);
        let filters = rec
            .filters()
            .as_ref()
            .iter()
            .filter(|f| *f != "PASS")
            .map(|f| filter_index[f.as_str()])
            .collect();
        let mut info = SmallVec::new();
        for (key, value) in rec.info().as_ref() {
            let idx = INFO_IDS.iter().position(|k| k == key).unwrap() as u8;
            let v = match value {
                Some(IV::Integer(n)) => Val::Int(*n),
                Some(IV::Float(f)) => Val::Float(*f),
                Some(IV::Flag) => Val::Flag,
                Some(IV::String(s)) => Val::Str(SmolStr::new(s)),
                Some(IV::Array(IA::Integer(a))) => {
                    Val::Ints(a.iter().map(|x| x.unwrap_or(i32::MIN)).collect())
                }
                Some(IV::Array(IA::Float(a))) => {
                    Val::Floats(a.iter().map(|x| x.unwrap_or(f32::NAN)).collect())
                }
                other => panic!("unexpected INFO {key}: {other:?}"),
            };
            info.push((idx, v));
        }
        let mut format = SmallVec::new();
        let samples = rec.samples();
        let keys: Vec<&str> = samples.keys().as_ref().iter().map(String::as_str).collect();
        if let Some(sample) = samples.values().next() {
            for (key, value) in keys.iter().zip(sample.values()) {
                let idx = FORMAT_IDS.iter().position(|k| k == key).unwrap() as u8;
                let Some(value) = value else { continue };
                let scalar = matches!(*key, "GL" | "GC" | "DP");
                let v = match value {
                    SV::Genotype(g) => {
                        let a: Vec<Option<usize>> =
                            g.as_ref().iter().map(|a| a.position()).collect();
                        Val::Gt(match a.as_slice() {
                            [Some(x), Some(y)] => Genotype::unphased(*x as u16, *y as u16),
                            _ => Genotype::missing_diploid(),
                        })
                    }
                    SV::Integer(n) if scalar => Val::Int(*n),
                    SV::Integer(n) => Val::Ints(SmallVec::from_iter([*n])),
                    SV::Float(f) if scalar => Val::Float(*f),
                    SV::Float(f) => Val::Floats(SmallVec::from_iter([*f])),
                    SV::Array(SA::Float(a)) if scalar => Val::Float(a[0].unwrap_or(f32::NAN)),
                    SV::Array(SA::Float(a)) => {
                        Val::Floats(a.iter().map(|x| x.unwrap_or(f32::NAN)).collect())
                    }
                    SV::Array(SA::Integer(a)) if scalar => Val::Int(a[0].unwrap_or(i32::MIN)),
                    SV::Array(SA::Integer(a)) => {
                        Val::Ints(a.iter().map(|x| x.unwrap_or(i32::MIN)).collect())
                    }
                    other => panic!("unexpected FORMAT {key}: {other:?}"),
                };
                format.push((idx, v));
            }
        }
        out.push(Rec { contig, pos, alleles, qual: rec.quality_score(), filters, info, format });
    }
    (contigs, out)
}

fn encode_rec<W: Write>(w: &mut Writer<W, Ready>, s: &Schema, r: &Rec) {
    let contig = &s.contigs[r.contig as usize].1;
    let enc = w.begin_record(contig, r.pos, &r.alleles, r.qual).unwrap();
    let mut enc = if r.filters.is_empty() {
        enc.filter_pass()
    } else {
        enc.filter_fail(r.filters.iter().map(|&i| &s.filters[i as usize].1))
    };
    for (idx, v) in &r.info {
        let e = &mut enc;
        match (*idx, v) {
            (0, Val::Ints(a)) => s.ad.encode(e, a),
            (1, Val::Float(f)) => s.bq.encode(e, *f),
            (2, Val::Int(n)) => s.dp.encode(e, *n),
            (3, Val::Float(f)) => s.mq.encode(e, *f),
            (4, Val::Int(n)) => s.mq0.encode(e, *n),
            (5, Val::Int(n)) => s.ns.encode(e, *n),
            (6, Val::Ints(a)) => s.as_sb_ot.encode(e, a),
            (7, Val::Ints(a)) => s.as_sb_ob.encode(e, a),
            (8, Val::Str(x)) => s.sc5.encode(e, x.as_str()),
            (9, Val::Floats(a)) => s.af.encode(e, a),
            (10, Val::Floats(a)) => s.abq.encode(e, a),
            (11, Val::Floats(a)) => s.amq.encode(e, a),
            (12, Val::Floats(a)) => s.as_ss_bq_ot.encode(e, a),
            (13, Val::Floats(a)) => s.as_ss_bq_ob.encode(e, a),
            (14, Val::Floats(a)) => s.as_ss_mq_ot.encode(e, a),
            (15, Val::Floats(a)) => s.as_ss_mq_ob.encode(e, a),
            (16, Val::Floats(a)) => s.pir.encode(e, a),
            (17, Val::Float(f)) => s.ent100.encode(e, *f),
            (18, Val::Floats(a)) => s.nab.encode(e, a),
            (19, Val::Floats(a)) => s.noi.encode(e, a),
            (20, Val::Ints(a)) => s.m5mc_strands.encode(e, a),
            (21, Val::Flag) => s.cpg.encode(e),
            (22, Val::Flag) => s.cpgnovo.encode(e),
            _ => panic!("INFO kind mismatch for {}", INFO_IDS[*idx as usize]),
        }
    }
    let mut enc = enc.begin_samples();
    for (idx, v) in &r.format {
        let e = &mut enc;
        match (*idx, v) {
            (0, Val::Gt(g)) => s.gt.encode(e, std::slice::from_ref(g)),
            (1, Val::Float(f)) => s.gl.encode(e, &[*f]),
            (2, Val::Float(f)) => s.gc.encode(e, &[*f]),
            (3, Val::Int(n)) => s.fdp.encode(e, &[*n]),
            (4, Val::Floats(a)) => s.m5mc.encode_single_sample(e, a),
            (5, Val::Ints(a)) => s.dpm5mc.encode_single_sample(e, a),
            (6, Val::Ints(a)) => s.adm5mc.encode_single_sample(e, a),
            (7, Val::Floats(a)) => s.ml.encode_single_sample(e, a),
            _ => panic!("FORMAT kind mismatch for {}", FORMAT_IDS[*idx as usize]),
        }
        .unwrap();
    }
    enc.emit().unwrap();
}

fn bcf_replay(args: &[String]) {
    let path = Path::new(&args[2]);
    let fmt = match args[3].as_str() {
        "bcf" => OutputFormat::Bcf,
        "vcf" => OutputFormat::Vcf,
        "vcfgz" => OutputFormat::VcfGz,
        other => panic!("format {other}"),
    };
    let iters: u32 = args[4].parse().unwrap();
    let out = out_path(args.get(5), "prof_writers_out.bcf");
    let t = Instant::now();
    let (contigs, recs) = load_bcf(path);
    println!("loaded {} records in {:.1} s", recs.len(), t.elapsed().as_secs_f64());
    let (header, schema) = build_schema(&contigs);
    let threads: usize = args.get(6).map_or(0, |s| s.parse().unwrap());
    timed("bcf-replay", iters, || write_vcf(&out, fmt, threads, &header, &schema, &recs));
}

#[inline(never)]
fn write_vcf(
    out: &Path,
    fmt: OutputFormat,
    threads: usize,
    header: &VcfHeader,
    schema: &Schema,
    recs: &[Rec],
) {
    // rastair writes through a `Box<dyn Write + Send>` over an unbuffered file.
    let file: Box<dyn Write + Send> = Box::new(std::fs::File::create(out).unwrap());
    let w = Writer::new(file, fmt).compression_threads(threads).unwrap();
    let mut w = w.write_header(header).unwrap();
    for r in recs {
        encode_rec(&mut w, schema, r);
    }
    let (_, index) = w.finish().unwrap();
    if let Some(index) = index {
        let mut csi = Vec::new();
        index.write(&mut csi).unwrap();
        let mut path = out.as_os_str().to_owned();
        path.push(".csi");
        std::fs::write(path, csi).unwrap();
    }
}

// ── BAM rewrite ─────────────────────────────────────────────────────────

/// rastair-like `MM`/`ML`: one modification call per `CpG` C on the stored strand.
fn mm_ml(seq: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let mut mm = b"C+m".to_vec();
    let mut ml = Vec::new();
    let mut skipped = 0u32;
    let mut itoa = itoa::Buffer::new();
    for (i, w) in seq.windows(2).enumerate() {
        if w[0] != b'C' {
            continue;
        }
        if w[1] == b'G' {
            mm.push(b',');
            mm.extend_from_slice(itoa.format(skipped).as_bytes());
            ml.push(((i * 37) % 256) as u8);
            skipped = 0;
        } else {
            skipped += 1;
        }
    }
    mm.push(b';');
    (mm, ml)
}

fn parse_region(args: &[String]) -> (String, u32, u32) {
    (args[3].clone(), args[4].parse().unwrap(), args[5].parse().unwrap())
}

fn bam_seqair(args: &[String]) {
    use seqair::bam::{AuxData, BamHeader, OwnedBamRecord, RecordStore};
    let path = Path::new(&args[2]);
    let (contig, start, end) = parse_region(args);
    let level: i32 = args[6].parse().unwrap();
    let index = args[7] == "1";
    let iters: u32 = args[8].parse().unwrap();
    let out = out_path(args.get(9), "prof_writers_out.bam");

    let mut reader = seqair::bam::IndexedBamReader::open(path).unwrap();
    let header = BamHeader::from_template(reader.header());
    let tid = reader.header().tid(&contig).unwrap();
    let mut store = RecordStore::new();
    let span = RangeInclusive { start: Pos0::new(start).unwrap(), last: Pos0::new(end).unwrap() };
    reader.fetch_into(tid, span, &mut store).unwrap();
    let recs: Vec<OwnedBamRecord> = store
        .indices()
        .map(|idx| {
            let r = store.record(idx).unwrap();
            let seq = r.seq().to_vec();
            let ascii: Vec<u8> = seq.iter().map(|b| b.as_char() as u8).collect();
            let mut aux = AuxData::from_bytes(r.aux().to_vec());
            let (mm, ml) = mm_ml(&ascii);
            if !ml.is_empty() {
                aux.set_string(*b"MM", &mm);
                aux.set_array_u8(*b"ML", &ml).unwrap();
            }
            OwnedBamRecord {
                ref_id: r.tid,
                pos: Some(r.pos),
                mapq: r.mapq,
                flags: r.flags,
                next_ref_id: r.next_ref_id,
                next_pos: u32::try_from(r.next_pos).ok().and_then(Pos0::new),
                template_len: r.template_len,
                qname: r.qname().to_vec(),
                cigar: r.cigar().to_vec(),
                seq,
                qual: r.qual().to_vec(),
                aux,
            }
        })
        .collect();
    println!("{} records", recs.len());
    if level < 0 {
        // Serialization only: `to_bam_bytes` without BGZF or index.
        timed("bam-serialize", iters, || serialize_bam(&recs));
        return;
    }
    let threads: usize = args.get(10).map_or(0, |s| s.parse().unwrap());
    timed("bam-seqair", iters, || write_bam_seqair(&out, &header, level, index, threads, &recs));
}

#[inline(never)]
fn serialize_bam(recs: &[seqair::bam::OwnedBamRecord]) {
    let mut buf = Vec::with_capacity(1 << 16);
    let mut total = 0usize;
    for rec in recs {
        buf.clear();
        rec.to_bam_bytes(&mut buf).unwrap();
        total += buf.len();
    }
    std::hint::black_box(total);
}

#[inline(never)]
fn write_bam_seqair(
    out: &Path,
    header: &seqair::bam::BamHeader,
    level: i32,
    index: bool,
    threads: usize,
    recs: &[seqair::bam::OwnedBamRecord],
) {
    let mut writer = seqair::bam::BamWriterBuilder::to_path(out, header)
        .compression_level(level)
        .compression_threads(threads)
        .write_index(index)
        .build()
        .unwrap();
    for rec in recs {
        writer.write(rec).unwrap();
    }
    let (_, idx) = writer.finish().unwrap();
    if let Some(idx) = idx {
        let mut bai = Vec::new();
        idx.write_bai(&mut bai, header.target_count()).unwrap();
        std::fs::write(out.with_extension("bam.bai"), bai).unwrap();
    }
}

fn bam_htslib(args: &[String]) {
    use rust_htslib::bam::record::{Aux, AuxArray};
    use rust_htslib::bam::{self, Read};
    let path = Path::new(&args[2]);
    let (contig, start, end) = parse_region(args);
    let level: u32 = args[6].parse().unwrap();
    let threads: usize = args[7].parse().unwrap();
    let iters: u32 = args[8].parse().unwrap();
    let out = out_path(args.get(9), "prof_writers_out_hts.bam");

    let mut reader = bam::IndexedReader::from_path(path).unwrap();
    reader.fetch((contig.as_str(), i64::from(start), i64::from(end) + 1)).unwrap();
    let header = bam::Header::from_template(reader.header());
    let mut recs = Vec::new();
    let mut rec = bam::Record::new();
    while let Some(r) = reader.read(&mut rec) {
        r.unwrap();
        let (mm, ml) = mm_ml(&rec.seq().as_bytes());
        if !ml.is_empty() {
            rec.push_aux(b"MM", Aux::String(std::str::from_utf8(&mm).unwrap())).unwrap();
            rec.push_aux(b"ML", Aux::ArrayU8(AuxArray::from(&ml))).unwrap();
        }
        recs.push(rec.clone());
    }
    println!("{} records", recs.len());
    timed("bam-htslib", iters, || {
        let mut w = bam::Writer::from_path(&out, &header, bam::Format::Bam).unwrap();
        w.set_compression_level(bam::CompressionLevel::Level(level)).unwrap();
        if threads > 1 {
            w.set_threads(threads).unwrap();
        }
        for r in &recs {
            w.write(r).unwrap();
        }
        drop(w);
    });
}

// ── libdeflate level sweep ──────────────────────────────────────────────

/// Compress an uncompressed stream in BGZF-sized blocks at every libdeflate
/// level; print compressed size and single-thread CPU time.
fn levels(args: &[String]) {
    let data = std::fs::read(&args[2]).unwrap();
    let block: usize = args.get(3).map_or(65536, |s| s.parse().unwrap());
    let mut out = vec![0u8; 2 * block];
    println!("{} bytes, {block}-byte blocks", data.len());
    let max: i32 = args.get(4).map_or(9, |s| s.parse().unwrap());
    for level in 0..=max {
        let mut c = libdeflater::Compressor::new(libdeflater::CompressionLvl::new(level).unwrap());
        let mut total = 0usize;
        let c0 = cpu_seconds();
        for chunk in data.chunks(block) {
            total += c.deflate_compress(chunk, &mut out).unwrap() + 26;
        }
        let cpu = cpu_seconds() - c0;
        println!(
            "level {level:2}: {:>11} bytes  ratio {:.3}  cpu {:.3} s  {:.0} MB/s",
            total,
            data.len() as f64 / total as f64,
            cpu,
            data.len() as f64 / cpu / 1e6
        );
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    match args[1].as_str() {
        "bcf-replay" => bcf_replay(&args),
        "bam-seqair" => bam_seqair(&args),
        "bam-htslib" => bam_htslib(&args),
        "levels" => levels(&args),
        other => panic!("mode {other}"),
    }
}
