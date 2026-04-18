//! Empirical calibration for `RecordStore::with_byte_hint`.
//!
//! Measures actual per-record slab sizes and compression ratios across a
//! user-supplied corpus of BAM files, then reports what the constants in
//! `with_byte_hint` *should* be for that workload.
//!
//! This is `#[ignore]`d so it stays out of CI — the BAMs are user-provided
//! and paths vary. Run with:
//!
//! ```text
//! SEQAIR_CALIBRATE_BAMS=/path/a.bam,/path/b.bam \
//!   cargo test -p seqair --test calibrate_byte_hint -- --ignored --nocapture
//! ```
//!
//! Set `SEQAIR_CALIBRATE_REGION_SIZE` (bp, default `100_000`) and
//! `SEQAIR_CALIBRATE_REGIONS_PER_BAM` (default 4) to tune sampling.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::print_stdout,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::arithmetic_side_effects,
    reason = "calibration harness, not production code"
)]

use std::path::{Path, PathBuf};

use seqair::bam::Pos0;
use seqair::bam::reader::IndexedBamReader;
use seqair::bam::record_store::RecordStore;

/// Measurements for one region load.
#[derive(Debug, Clone, Copy)]
struct Sample {
    /// Compressed bytes the reader would pass to `with_byte_hint` for this region.
    compressed_bytes: usize,
    records: usize,
    names_bytes: usize,
    bases_bytes: usize,
    cigar_bytes: usize,
    qual_bytes: usize,
    aux_bytes: usize,
}

impl Sample {
    fn uncompressed_total(&self) -> usize {
        self.names_bytes
            + self.bases_bytes
            + self.cigar_bytes
            + self.qual_bytes
            + self.aux_bytes
            // fixed-size record struct (SlimRecord is 64 bytes)
            + self.records * 64
    }
}

/// Running per-BAM aggregate.
#[derive(Default)]
struct Aggregate {
    label: String,
    samples: Vec<Sample>,
}

impl Aggregate {
    fn total_records(&self) -> usize {
        self.samples.iter().map(|s| s.records).sum()
    }
    fn total_compressed(&self) -> usize {
        self.samples.iter().map(|s| s.compressed_bytes).sum()
    }
    fn total_names(&self) -> usize {
        self.samples.iter().map(|s| s.names_bytes).sum()
    }
    fn total_bases(&self) -> usize {
        self.samples.iter().map(|s| s.bases_bytes).sum()
    }
    fn total_cigar(&self) -> usize {
        self.samples.iter().map(|s| s.cigar_bytes).sum()
    }
    fn total_qual(&self) -> usize {
        self.samples.iter().map(|s| s.qual_bytes).sum()
    }
    fn total_aux(&self) -> usize {
        self.samples.iter().map(|s| s.aux_bytes).sum()
    }
    fn total_uncompressed(&self) -> usize {
        self.samples.iter().map(Sample::uncompressed_total).sum()
    }
}

fn env_usize(key: &str, default: usize) -> usize {
    std::env::var(key).ok().and_then(|v| v.parse().ok()).unwrap_or(default)
}

fn bam_paths() -> Vec<PathBuf> {
    match std::env::var("SEQAIR_CALIBRATE_BAMS") {
        Ok(s) => s.split(',').filter(|p| !p.is_empty()).map(PathBuf::from).collect(),
        Err(_) => Vec::new(),
    }
}

/// Sample `regions_per_bam` windows of `region_size` bp each, spread across
/// the first reference that has indexable data.
fn sample_bam(path: &Path, region_size: u64, regions_per_bam: usize) -> Aggregate {
    let mut agg = Aggregate { label: path.display().to_string(), samples: Vec::new() };

    let mut reader = match IndexedBamReader::open(path) {
        Ok(r) => r,
        Err(e) => {
            println!("  ! failed to open {}: {e}", path.display());
            return agg;
        }
    };

    // Pick the first reference that produces non-empty query chunks at a
    // mid-chromosome position. For whole-genome BAMs this is usually chr1;
    // for chr-restricted BAMs (e.g. chr12-only) it's whichever chromosome
    // has actual alignments.
    let mut picked: Option<(u32, String, u64)> = None;
    for tid in 0..reader.header().target_count() as u32 {
        let name = reader.header().target_name(tid).map(str::to_owned).unwrap_or_default();
        let len = reader.header().target_len(tid).unwrap_or(0);
        if len < region_size * (regions_per_bam as u64 + 2) {
            continue;
        }
        let probe_start = Pos0::new(u32::try_from(len / 4).expect("probe_start fits u32")).unwrap();
        let probe_end =
            Pos0::new(u32::try_from(len / 4 + region_size).expect("probe_end fits u32")).unwrap();
        if !reader.shared().index().query(tid, probe_start, probe_end).is_empty() {
            picked = Some((tid, name, len));
            break;
        }
    }

    let Some((tid, name, ref_len)) = picked else {
        println!("  ! no reference with indexable data found");
        return agg;
    };

    println!("  using ref tid={tid} ({name}, len={ref_len})");

    // Spread sample windows evenly across [ref_len * 0.2, ref_len * 0.8] so
    // we avoid telomeres/centromeres that can be sparsely covered.
    let span_start = ref_len / 5;
    let span_end = ref_len * 4 / 5;
    let step =
        if regions_per_bam > 1 { (span_end - span_start) / (regions_per_bam as u64) } else { 0 };

    let mut store = RecordStore::new();
    for i in 0..regions_per_bam {
        let start_bp = span_start + step * (i as u64);
        let end_bp = start_bp + region_size;
        let start = Pos0::new(u32::try_from(start_bp).expect("start fits u32")).unwrap();
        let end = Pos0::new(u32::try_from(end_bp).expect("end fits u32")).unwrap();

        // Estimate compressed bytes by summing (merged) chunk byte spans.
        let chunks = reader.shared().index().query(tid, start, end);
        if chunks.is_empty() {
            continue;
        }
        let compressed_bytes: usize = {
            // Merge overlapping/adjacent ranges by file offset — BAI chunks can overlap.
            let mut ranges: Vec<(u64, u64)> = chunks
                .iter()
                .map(|c| (c.begin.block_offset(), c.end.block_offset().max(c.begin.block_offset())))
                .collect();
            ranges.sort_unstable();
            let mut merged: Vec<(u64, u64)> = Vec::with_capacity(ranges.len());
            for r in ranges {
                match merged.last_mut() {
                    Some(last) if r.0 <= last.1 => last.1 = last.1.max(r.1),
                    _ => merged.push(r),
                }
            }
            merged.iter().map(|r| (r.1 - r.0) as usize).sum()
        };

        store.clear();
        let n = match reader.fetch_into(tid, start, end, &mut store) {
            Ok(n) => n,
            Err(e) => {
                println!("    [skip] fetch {start_bp}..{end_bp}: {e}");
                continue;
            }
        };
        if n == 0 {
            continue;
        }

        // Walk each record and sum per-slab bytes.
        let mut names_bytes = 0usize;
        let mut bases_bytes = 0usize;
        let mut cigar_bytes = 0usize;
        let mut qual_bytes = 0usize;
        let mut aux_bytes = 0usize;
        for i in 0..store.len() as u32 {
            names_bytes += store.qname(i).len();
            bases_bytes += store.seq(i).len();
            cigar_bytes += store.cigar(i).len();
            qual_bytes += store.qual(i).len();
            aux_bytes += store.aux(i).len();
        }

        let sample = Sample {
            compressed_bytes,
            records: n,
            names_bytes,
            bases_bytes,
            cigar_bytes,
            qual_bytes,
            aux_bytes,
        };
        println!(
            "    {start_bp:>12}..{end_bp:<12}  recs={n:>6}  comp={compressed_bytes:>10}  \
             uncomp={:>10}  ratio={:.2}x  ~{}B/rec",
            sample.uncompressed_total(),
            sample.uncompressed_total() as f64 / compressed_bytes.max(1) as f64,
            sample.uncompressed_total() / n.max(1),
        );
        agg.samples.push(sample);
    }

    agg
}

fn report(aggs: &[Aggregate]) {
    println!();
    println!("=== per-BAM summary ===");
    for a in aggs {
        if a.samples.is_empty() {
            println!("  {} — no samples", a.label);
            continue;
        }
        let recs = a.total_records() as f64;
        let comp = a.total_compressed() as f64;
        let uncomp = a.total_uncompressed() as f64;
        println!(
            "  {}\n    samples={}  records={}  comp={}  uncomp={}  ratio={:.2}x",
            a.label,
            a.samples.len(),
            a.total_records(),
            a.total_compressed(),
            a.total_uncompressed(),
            uncomp / comp,
        );
        println!(
            "    per-record: total={:>6.1}B  names={:>5.1}  bases={:>5.1}  cigar={:>5.1}  \
             qual={:>5.1}  aux={:>6.1}",
            uncomp / recs,
            a.total_names() as f64 / recs,
            a.total_bases() as f64 / recs,
            a.total_cigar() as f64 / recs,
            a.total_qual() as f64 / recs,
            a.total_aux() as f64 / recs,
        );
        println!(
            "    records / compressed byte: {:.4} (i.e. {:.1} compressed bytes per record)",
            recs / comp,
            comp / recs,
        );
    }

    // Pooled recommendation: sum across all BAMs.
    let tot_recs: f64 = aggs.iter().map(|a| a.total_records() as f64).sum();
    let tot_comp: f64 = aggs.iter().map(|a| a.total_compressed() as f64).sum();
    let tot_uncomp: f64 = aggs.iter().map(|a| a.total_uncompressed() as f64).sum();
    if tot_recs == 0.0 {
        return;
    }
    let tot_names: f64 = aggs.iter().map(|a| a.total_names() as f64).sum();
    let tot_bases: f64 = aggs.iter().map(|a| a.total_bases() as f64).sum();
    let tot_cigar: f64 = aggs.iter().map(|a| a.total_cigar() as f64).sum();
    let tot_qual: f64 = aggs.iter().map(|a| a.total_qual() as f64).sum();
    let tot_aux: f64 = aggs.iter().map(|a| a.total_aux() as f64).sum();

    println!();
    println!("=== pooled recommendation ===");
    println!("  compression ratio:        {:.2}x", tot_uncomp / tot_comp);
    println!("  bytes per record (total): {:.1}", tot_uncomp / tot_recs);
    println!("  names per record:         {:.1}", tot_names / tot_recs);
    println!("  bases per record:         {:.1}", tot_bases / tot_recs);
    println!("  cigar per record:         {:.1}", tot_cigar / tot_recs);
    println!("  qual  per record:         {:.1}", tot_qual / tot_recs);
    println!("  aux   per record:         {:.1}", tot_aux / tot_recs);
    println!();
    println!("Current `with_byte_hint` assumes: 3.00x compression, 400 B/rec,");
    println!("names=25, cigar=20, bases=150, qual=150, aux=residual.");
}

// Run manually — not part of CI. See module doc-comment.
#[test]
#[ignore = "requires SEQAIR_CALIBRATE_BAMS env var pointing at user-owned BAMs"]
fn calibrate_byte_hint() {
    let paths = bam_paths();
    if paths.is_empty() {
        println!(
            "SEQAIR_CALIBRATE_BAMS not set. Example:\n  \
             SEQAIR_CALIBRATE_BAMS=/path/one.bam,/path/two.bam \\\n  \
             cargo test -p seqair --test calibrate_byte_hint -- --ignored --nocapture"
        );
        return;
    }

    let region_size = env_usize("SEQAIR_CALIBRATE_REGION_SIZE", 100_000) as u64;
    let regions_per_bam = env_usize("SEQAIR_CALIBRATE_REGIONS_PER_BAM", 4);

    println!(
        "calibrating against {} BAM(s), {} regions × {} bp each",
        paths.len(),
        regions_per_bam,
        region_size,
    );

    let mut aggs = Vec::with_capacity(paths.len());
    for path in &paths {
        println!();
        println!("== {} ==", path.display());
        aggs.push(sample_bam(path, region_size, regions_per_bam));
    }

    report(&aggs);
}
