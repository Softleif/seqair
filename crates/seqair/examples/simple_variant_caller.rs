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
use seqair::reader::{Readers, ResolveTid, SegmentOptions};
use seqair::vcf::{
    Alleles, FormatGt, FormatInt, Genotype, InfoInt, Number, OutputFormat, ValueType,
    VcfHeaderBuilder, Writer,
    record_encoder::{Arr, FormatFieldDef, Gt, InfoFieldDef, Scalar},
};
use seqair_types::{Base, Pos1, RegionString};
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
    region: RegionString,

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

    // ── Open readers with a push-time filter ───────────────────────────
    // Skip reads that are unmapped, secondary, supplementary, QC-failed, or
    // duplicate before they enter the store. Mapping quality is checked per
    // alignment in the pileup loop below since `min_mapq` is dynamic.
    use seqair::bam::record_store::{CustomizeRecordStore, RecordStore, SlimRecord};
    #[derive(Clone)]
    struct DropUseless;
    impl CustomizeRecordStore for DropUseless {
        type Extra = ();
        fn filter(&mut self, rec: &SlimRecord, _: &RecordStore<()>) -> bool {
            !rec.flags.is_unmapped()
                && !rec.flags.is_secondary()
                && !rec.flags.is_supplementary()
                && !rec.flags.is_failed_qc()
                && !rec.flags.is_duplicate()
        }
        fn compute(&mut self, _: &SlimRecord, _: &RecordStore<()>) {}
    }

    let mut readers =
        Readers::<DropUseless>::open_customized(&args.input, &args.reference, DropUseless)
            .context("could not open BAM + FASTA")?;

    // ── Build VCF header ───────────────────────────────────────────────
    // from_bam_header copies all @SQ lines so contig names/lengths match.
    let from_bam = VcfHeaderBuilder::from_bam_header(readers.header())
        .context("could not build VCF header from BAM")?;

    // All segments of a single-contig region share one ContigId, so resolve
    // it once via the new `FromBamHeader::contig(tid)` helper before we
    // partial-move `from_bam.builder` into the typestate chain below.
    let region_tid = args.region.chromosome.as_str().resolve_tid(readers.header())?;
    let region_contig =
        from_bam.contig(region_tid).context("contig not registered in VCF header")?.clone();

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
    let min_mapq = args.min_mapq;

    // Tile the requested region into 100 kb segments. This example processes
    // them sequentially; a real caller would parallelize across forks. The
    // parsed RegionString is fed straight to `segments()` — no manual region
    // parsing needed.
    let max_len = std::num::NonZeroU32::new(100_000).expect("non-zero literal");
    let opts = SegmentOptions::new(max_len);
    let plan: Vec<_> = readers.segments(&args.region, opts)?.collect();
    let mut n_calls = 0u32;

    for segment in &plan {
        let mut pileup = readers.pileup(segment)?;

        while let Some(column) = pileup.pileups() {
            // Restrict emission to the segment's "owned" core range so that
            // adjacent overlapping segments don't double-call the same site.
            let core = segment.core_range();
            if !core.contains(&column.pos()) {
                continue;
            }
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
            let enc = vcf_writer.begin_record(&region_contig, pos1, &alleles, None)?;
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

        // `pileup` drops at end of scope, returning the store to `readers`.
    }

    vcf_writer.finish()?;
    output.flush()?;

    eprintln!(
        "called {n_calls} SNV(s) in {region} (min_depth={}, min_af={:.2})",
        args.min_depth,
        args.min_af,
        region = args.region,
    );

    Ok(())
}
