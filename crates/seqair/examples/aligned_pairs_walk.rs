#![allow(
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    clippy::doc_markdown,
    clippy::indexing_slicing,
    clippy::print_stdout,
    reason = "example"
)]

//! Per-record AlignedPairs walking — three layers, three calling tasks.
//!
//! This example demonstrates the [`AlignedPairs`](seqair::bam::AlignedPairs)
//! iterator and its two layered adapters. Each adapter pays only for the
//! data it attaches, so callers compose the view they need:
//!
//! 1. **Bare walk** — [`SlimRecord::aligned_pairs(store)`] yields events
//!    with positions and op kinds only. Cheapest; useful for things like
//!    counting indel events per record.
//! 2. **With read** — `.with_read(seq, qual)` (or the one-shot
//!    [`SlimRecord::aligned_pairs_with_read`]) attaches the read's bases
//!    and quality scores. `Match` carries `query: Base` + `qual`; `Insertion`
//!    and `SoftClip` carry pre-sliced runs (no off-by-one risk for callers).
//! 3. **With reference** — `.with_reference(&ref_seq)` further attaches
//!    [`Base`] from the loaded reference window. `Match` gains
//!    `ref_base: Option<Base>` (`None` outside the loaded window);
//!    `Deletion` gains `ref_bases: Option<&[Base]>`.
//!
//! Three example calling tasks share the same walk:
//!
//! - **Methylation (TAPS)**: count C→T conversions where the reference base
//!   is `C` and the read base is `T` with a Phred quality threshold.
//! - **SNV evidence**: count per-record mismatches where the read base
//!   differs from the reference base.
//! - **Indel evidence**: count insertions and deletions per record, with
//!   their lengths.
//!
//! ## Run
//!
//! ```bash
//! cargo run --example aligned_pairs_walk -- \
//!     --input sample.bam --reference ref.fa --region chr1:1000-2000
//! ```
//!
//! The BAM and FASTA must be indexed.

use anyhow::Context;
use clap::Parser as _;
use seqair::bam::record_store::{CustomizeRecordStore, RecordStore, SlimRecord};
use seqair::bam::{AlignedPair, AlignedPairWithRef};
use seqair::bam::{MatchKind, RefSeq};
use seqair::reader::{Readers, ResolveTid};
use seqair_types::{Base, RegionString};
use std::path::PathBuf;
use std::rc::Rc;

#[derive(Debug, clap::Parser)]
struct Cli {
    /// BAM/SAM/CRAM file (must be indexed).
    #[clap(long, short)]
    input: PathBuf,

    /// Reference FASTA (must be indexed).
    #[clap(long, short)]
    reference: PathBuf,

    /// Region to scan (e.g. "chr1:1000-2000").
    #[clap(long, short)]
    region: RegionString,

    /// Phred quality threshold for methylation evidence.
    #[clap(long, default_value_t = 20)]
    min_qual: u8,
}

/// A push-time filter that drops records that won't produce useful
/// evidence — same defensive set as the variant caller example.
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

#[derive(Default, Debug)]
struct Counts {
    records: u32,
    methylation_calls: u32,
    snv_evidence: u32,
    insertions: u32,
    deletions: u32,
    inserted_bases: u32,
    deleted_bases: u32,
}

fn main() -> anyhow::Result<()> {
    let args = Cli::parse();

    let mut readers =
        Readers::<DropUseless>::open_customized(&args.input, &args.reference, DropUseless)
            .context("could not open BAM + FASTA")?;

    // Resolve the region against the BAM header.
    let tid = args.region.chromosome.as_str().resolve_tid(readers.header())?;
    let start = args
        .region
        .start
        .context("region must have a start (e.g. chr1:1000-2000)")?
        .to_zero_based();
    let end =
        args.region.end.context("region must have an end (e.g. chr1:1000-2000)")?.to_zero_based();

    // Fetch records into a freshly-built RecordStore. (For repeated calls
    // you'd reuse one — see Readers::pileup which retains capacity via the
    // PileupGuard. This example keeps the data flow visible.)
    let mut store = RecordStore::<()>::new();
    readers.fetch_into(tid.as_u32(), start, end, &mut store)?;

    // Load the reference window once for the whole region. RefSeq wraps an
    // Rc<[Base]>; the layered iterator borrows it for the duration of the walk.
    let ref_bases: Rc<[Base]> =
        readers.fetch_base_seq(args.region.chromosome.as_str(), start, end)?;
    let ref_seq = RefSeq::new(ref_bases, start);

    let mut counts = Counts::default();

    // Iterate every record in the fetched window.
    for idx in 0..store.len() as u32 {
        counts.records += 1;
        let rec = store.record(idx);

        // ── Layer 1: indel evidence (no read or ref data needed) ───────────
        for pair in rec.aligned_pairs(&store)? {
            match pair {
                AlignedPair::Insertion { insert_len, .. } => {
                    counts.insertions += 1;
                    counts.inserted_bases += insert_len;
                }
                AlignedPair::Deletion { del_len, .. } => {
                    counts.deletions += 1;
                    counts.deleted_bases += del_len;
                }
                _ => {}
            }
        }

        // ── Layers 2+3: methylation + SNV evidence (need ref + qual) ──────
        // One walk, two callers:
        //   - methylation: ref C + query T with qual >= threshold
        //   - SNV: any mismatch where ref_base is known
        // We use AlignedPairWithRef::Match's full destructure to express the
        // methylation predicate in the pattern itself — the compiler enforces
        // exhaustiveness across all variants.
        for ev in rec.aligned_pairs_with_read(&store)?.with_reference(&ref_seq) {
            if let AlignedPairWithRef::Match {
                kind: _, // methylation cares about M/=/X equally
                query,
                ref_base: Some(rb),
                qual,
                ..
            } = ev
            {
                // Methylation (TAPS C→T conversion).
                if rb == Base::C && query == Base::T && qual.get().unwrap_or(0) >= args.min_qual {
                    counts.methylation_calls += 1;
                }
                // Per-read SNV evidence — query disagrees with the reference.
                // (Real callers gate this on quality, MAPQ, strand bias, etc.)
                if query != rb && query != Base::Unknown && rb != Base::Unknown {
                    counts.snv_evidence += 1;
                }
            }
        }
    }

    println!("region: {}", args.region);
    println!("records walked:        {}", counts.records);
    println!(
        "methylation calls:     {} (ref C + query T, q≥{})",
        counts.methylation_calls, args.min_qual
    );
    println!("SNV evidence (matches): {}", counts.snv_evidence);
    println!("insertions: {} ({} bases)", counts.insertions, counts.inserted_bases);
    println!("deletions:  {} ({} bases)", counts.deletions, counts.deleted_bases);

    Ok(())
}

// ── Bonus: bare-walk indel-only counter ───────────────────────────────────
// Inline standalone function showing the cheapest possible layer — no read
// data, no reference, just CIGAR events. Useful when you only need event
// counts (for QC, sanity checks, or coarse filtering).
#[allow(dead_code, reason = "showcase function")]
fn count_indels_per_record(
    rec: &SlimRecord,
    store: &RecordStore<()>,
) -> anyhow::Result<(u32, u32)> {
    let mut ins = 0u32;
    let mut del = 0u32;
    for pair in rec.aligned_pairs(store)? {
        match pair {
            AlignedPair::Insertion { .. } => ins += 1,
            AlignedPair::Deletion { .. } => del += 1,
            _ => {}
        }
    }
    Ok((ins, del))
}

// ── Bonus: indel-evidence walker that exposes the inserted bases ───────────
// `with_read` gives Insertion variants pre-sliced bases and qual. This is
// what an indel-aware variant caller would use to record the inserted
// sequence as part of the variant's ALT allele.
#[allow(dead_code, reason = "showcase function")]
fn collect_insertions(
    rec: &SlimRecord,
    store: &RecordStore<()>,
) -> anyhow::Result<Vec<(u32, Vec<Base>)>> {
    use seqair::bam::AlignedPairWithRead;
    let mut out = Vec::new();
    for ev in rec.aligned_pairs_with_read(store)? {
        if let AlignedPairWithRead::Insertion { qpos, query, .. } = ev {
            out.push((qpos, query.to_vec()));
        }
    }
    Ok(out)
}

// ── Bonus: MatchKind dispatch ──────────────────────────────────────────────
// `=` and `X` CIGARs let an aligner record per-position match/mismatch
// information without forcing the caller to consult the reference. If a
// caller wants to count pre-known mismatches without a FASTA, `MatchKind`
// makes that a one-line filter.
#[allow(dead_code, reason = "showcase function")]
fn count_explicit_mismatches(rec: &SlimRecord, store: &RecordStore<()>) -> anyhow::Result<u32> {
    let mut count = 0u32;
    for pair in rec.aligned_pairs(store)? {
        if let AlignedPair::Match { kind: MatchKind::SeqMismatch, .. } = pair {
            count += 1;
        }
    }
    Ok(count)
}
