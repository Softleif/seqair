#![allow(
    clippy::arithmetic_side_effects,
    clippy::cast_lossless,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::doc_markdown,
    clippy::indexing_slicing,
    clippy::print_stdout,
    reason = "example"
)]

//! Toy single-sample SNV caller.
//!
//! Reads a BAM file with a reference FASTA, iterates the pileup in a region,
//! and emits a VCF record whenever the most-frequent non-reference base
//! exceeds a minimum allele frequency and depth.
//!
//! This is *not* a production variant caller — it has no base-quality
//! weighting, strand-bias filtering, or statistical model. Its purpose is
//! to demonstrate the full read → pileup → VCF-write pipeline:
//!
//! 1. Open alignment + reference with [`Readers`].
//! 2. Iterate columns with [`PileupEngine`].
//! 3. Count bases using [`Base::known_index`].
//! 4. Build a VCF header from the BAM header with [`VcfHeader::from_bam_header`].
//! 5. Encode records through the typestate chain:
//!    `begin_record` → `filter_pass` → `begin_samples` → `emit`.

use anyhow::Context;
use clap::Parser as _;
use seqair::bam::pileup::PileupOp;
use seqair::reader::Readers;
use seqair::vcf::{
    Alleles, FormatGt, FormatInt, Genotype, InfoInt, Number, OutputFormat, ValueType,
    VcfHeaderBuilder, Writer,
    record_encoder::{Arr, FormatFieldDef, Gt, InfoFieldDef, Scalar},
};
use seqair_types::{Base, Pos0, Pos1};
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::sync::Arc;

/// seqair simple-variant-caller — a toy SNV caller.
#[derive(Debug, clap::Parser)]
struct Cli {
    /// BAM/SAM/CRAM file (must be indexed).
    input: PathBuf,

    /// Reference FASTA (must be indexed).
    reference: PathBuf,

    /// Region to call (e.g. "chr1:1000-2000"). Required.
    #[clap(long, short)]
    region: String,

    /// Sample name to use in the VCF header.
    #[clap(long, default_value = "SAMPLE")]
    sample: String,

    /// Minimum depth at a position to consider calling.
    #[clap(long, default_value_t = 5)]
    min_depth: u32,

    /// Minimum alt allele frequency (0.0–1.0) to emit a call.
    #[clap(long, default_value_t = 0.1)]
    min_af: f64,

    /// Minimum mapping quality — reads below this are ignored.
    #[clap(long, default_value_t = 20)]
    min_mapq: u8,

    /// Output file. Defaults to stdout. Extension determines format
    /// (.vcf, .vcf.gz, .bcf).
    #[clap(long, short)]
    output: Option<PathBuf>,
}

/// Per-column base counts for A, C, G, T.
struct BaseCounts {
    counts: [i32; 4],
    total: i32,
}

impl BaseCounts {
    fn new() -> Self {
        Self { counts: [0; 4], total: 0 }
    }

    fn add(&mut self, base: Base) {
        if let Some(idx) = base.known_index() {
            self.counts[idx] += 1;
            self.total += 1;
        }
    }

    /// Returns (alt_base, alt_count) for the most frequent non-reference base,
    /// or `None` if every observed base matches the reference.
    fn best_alt(&self, ref_base: Base) -> Option<(Base, i32)> {
        const BASES: [Base; 4] = [Base::A, Base::C, Base::G, Base::T];
        let ref_idx = ref_base.known_index()?;

        let mut best: Option<(Base, i32)> = None;
        for (i, &count) in self.counts.iter().enumerate() {
            if i == ref_idx || count == 0 {
                continue;
            }
            match best {
                Some((_, prev)) if count <= prev => {}
                _ => best = Some((BASES[i], count)),
            }
        }
        best
    }
}

fn main() -> anyhow::Result<()> {
    let args = Cli::parse();

    // ── Open readers ───────────────────────────────────────────────────
    let mut readers =
        Readers::open(&args.input, &args.reference).context("could not open BAM + FASTA")?;

    let (tid, start, end) = parse_region(&args.region, readers.header())?;
    let contig_name = readers.header().target_name(tid).context("unknown tid")?.to_owned();

    // ── Build VCF header ───────────────────────────────────────────────
    // from_bam_header copies all @SQ lines so contig names/lengths match.
    let from_bam = VcfHeaderBuilder::from_bam_header(readers.header())
        .context("could not build VCF header from BAM")?;
    let mut builder = from_bam.builder.infos();

    // INFO fields: DP (total depth), AF (alt allele frequency)
    let dp_info: InfoInt = builder.register_info(&InfoFieldDef::<Scalar<i32>>::new(
        "DP",
        Number::Count(1),
        ValueType::Integer,
        "Total depth",
    ))?;
    let af_info: seqair::vcf::InfoFloats =
        builder.register_info(&InfoFieldDef::<Arr<f32>>::new(
            "AF",
            Number::AlternateBases,
            ValueType::Float,
            "Allele frequency",
        ))?;

    // FORMAT fields: GT (genotype), DP (sample depth), AD (allelic depths)
    let mut builder = builder.formats();
    let gt_fmt: FormatGt = builder.register_format(&FormatFieldDef::<Gt>::new(
        "GT",
        Number::Count(1),
        ValueType::String,
        "Genotype",
    ))?;
    let dp_fmt: FormatInt = builder.register_format(&FormatFieldDef::<Scalar<i32>>::new(
        "DP",
        Number::Count(1),
        ValueType::Integer,
        "Sample depth",
    ))?;

    let mut builder = builder.samples();
    builder.add_sample(&args.sample)?;
    let header = Arc::new(builder.build()?);

    // ── Set up VCF writer ──────────────────────────────────────────────
    let output_format = args
        .output
        .as_ref()
        .map(|p| OutputFormat::from_path(p))
        .transpose()?
        .unwrap_or(OutputFormat::Vcf);

    let mut output: Box<dyn Write> = if let Some(ref path) = args.output {
        Box::new(BufWriter::new(
            std::fs::File::create(path).with_context(|| format!("could not create {path:?}"))?,
        ))
    } else {
        Box::new(BufWriter::new(std::io::stdout().lock()))
    };

    let vcf_writer = Writer::new(&mut output, output_format);
    let mut vcf_writer = vcf_writer.write_header(&header)?;

    // ── Pileup and call variants ───────────────────────────────────────
    let start_pos = Pos0::new(start).context("invalid start position")?;
    let end_pos = Pos0::new(end).context("invalid end position")?;
    let min_mapq = args.min_mapq;

    let mut pileup = readers.pileup(tid, start_pos, end_pos)?;

    // Skip reads that are unmapped, secondary, supplementary, QC-failed,
    // duplicate, or below the minimum mapping quality.
    pileup.set_filter(move |flags, _cigar| {
        !flags.is_unmapped()
            && !flags.is_secondary()
            && !flags.is_supplementary()
            && !flags.is_failed_qc()
            && !flags.is_duplicate()
    });

    let mut n_calls = 0u32;

    while let Some(column) = pileup.pileups() {
        let ref_base = column.reference_base();
        // Skip positions where the reference is unknown (N).
        if ref_base == Base::Unknown {
            continue;
        }

        // Count bases, filtering by mapping quality.
        let mut bc = BaseCounts::new();
        for aln in column.alignments() {
            if aln.mapq < min_mapq {
                continue;
            }
            if let PileupOp::Match { base, .. } | PileupOp::Insertion { base, .. } = &aln.op {
                bc.add(*base);
            }
        }

        // Apply depth threshold.
        if bc.total < args.min_depth as i32 {
            continue;
        }

        // Find the best alt allele.
        let Some((alt_base, alt_count)) = bc.best_alt(ref_base) else {
            continue;
        };
        let af = alt_count as f64 / bc.total as f64;
        if af < args.min_af {
            continue;
        }

        // Convert 0-based pileup position to 1-based VCF position.
        let pos1 = Pos1::new(*column.pos() + 1).context("position overflow")?;

        // Build alleles — seqair enforces ref != alt at construction time.
        let alleles = Alleles::snv(ref_base, alt_base)?;

        // Encode through the typestate chain:
        //   begin_record → filter_pass (INFO) → begin_samples (FORMAT) → emit
        let enc = vcf_writer.begin_record(&from_bam.contigs[tid as usize], pos1, &alleles, None)?;
        let mut enc = enc.filter_pass();
        dp_info.encode(&mut enc, bc.total);
        af_info.encode(&mut enc, &[af as f32]);
        let mut enc = enc.begin_samples();

        // Simple genotype: 0/1 (het) if AF < 0.8, 1/1 (hom-alt) otherwise.
        let ref_count = bc.counts[ref_base.known_index().unwrap_or(0)];
        let gt = if af >= 0.8 { Genotype::unphased(1, 1) } else { Genotype::unphased(0, 1) };
        gt_fmt.encode(&mut enc, &[gt])?;
        dp_fmt.encode(&mut enc, &[bc.total])?;

        // AD: ref depth, alt depth — but we only registered DP and GT here,
        // so we skip AD to keep the example focused.
        _ = ref_count;

        enc.emit()?;
        n_calls += 1;
    }

    vcf_writer.finish()?;
    output.flush()?;
    readers.recover_store(&mut pileup);

    eprintln!(
        "called {n_calls} SNV(s) in {contig_name}:{start}-{end} (min_depth={}, min_af={:.2})",
        args.min_depth, args.min_af,
    );

    Ok(())
}

fn parse_region(region: &str, header: &seqair::bam::BamHeader) -> anyhow::Result<(u32, u32, u32)> {
    if let Some((name, range)) = region.split_once(':') {
        let tid = header.tid(name).with_context(|| format!("contig '{name}' not found"))?;
        let (start_str, end_str) =
            range.split_once('-').with_context(|| format!("invalid range: {range}"))?;
        let start: u32 = start_str
            .replace(',', "")
            .parse()
            .with_context(|| format!("invalid start: {start_str}"))?;
        let end: u32 =
            end_str.replace(',', "").parse().with_context(|| format!("invalid end: {end_str}"))?;
        Ok((tid, start.saturating_sub(1), end))
    } else {
        let tid = header.tid(region).with_context(|| format!("contig '{region}' not found"))?;
        let len = header.target_len(tid).unwrap_or(0) as u32;
        Ok((tid, 0, len))
    }
}
