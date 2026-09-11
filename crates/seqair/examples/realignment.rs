#![allow(
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    clippy::indexing_slicing,
    clippy::print_stdout,
    reason = "example"
)]

//! Toy local-realignment demo.
//!
//! Reads a BAM/SAM region into a [`RecordStore`], then "realigns" every
//! eligible record in place via [`RecordStore::set_alignment`]. The
//! realignment is deliberately fake — it does not improve any alignment.
//! Its only job is to showcase the realignment workflow:
//!
//! 1. Fetch records for a region.
//! 2. Propose a new (pos, cigar) per record. The query length must stay the
//!    same; everything else (seq, qual, qname, aux, flags, mapq) is preserved.
//! 3. Apply via `set_alignment`, then re-sort with `sort_by_pos` so the store
//!    is ready for pileup / writing again.
//!
//! The "realignment rule" here: for each record whose CIGAR starts with an
//! `M` (match) op of length ≥ `min_match`, convert the first `clip` bases
//! of that leading match into a soft clip and shift `pos` right by `clip`.
//! Example: `10M` at pos 100 → `2S8M` at pos 102 (when `clip = 2`).

use anyhow::Context;
use clap::Parser as _;
use seqair::bam::RecordIdx;
use seqair::{
    bam::{BamWriterBuilder, CigarOp, CigarStr, Pos0, RecordStore, cigar::CigarOpType},
    reader::{IndexedReader, ResolveTid},
};
use seqair_types::{RegionString, SmolStr, smol_str::ToSmolStr};
use std::path::PathBuf;

/// seqair realignment — a toy local-realignment example.
#[derive(Debug, clap::Parser)]
struct Cli {
    /// BAM/SAM file to read (must be indexed).
    input: PathBuf,

    /// Region to process. `chr`, `chr:start`, or `chr:start-end`
    /// (1-based, inclusive). Omit for the first contig in the header.
    #[clap(long, short)]
    region: Option<RegionString>,

    /// Soft-clip this many leading match bases per eligible record.
    #[clap(long, default_value_t = 2)]
    clip: u32,

    /// Only realign reads whose leading M op is at least this long.
    #[clap(long, default_value_t = 4)]
    min_match: u32,

    /// Print a before/after line for up to this many modified records.
    #[clap(long, default_value_t = 5)]
    show: usize,

    /// Write realigned records to this BAM file.
    #[clap(long, short)]
    output: Option<PathBuf>,
}

fn main() -> anyhow::Result<()> {
    let args = Cli::parse();

    let mut reader = IndexedReader::open(&args.input).context("could not open alignment file")?;

    let (tid, start, end, contig) = resolve_region(args.region.as_ref(), reader.header())?;

    let mut store = RecordStore::new();
    reader.fetch_into(tid, start, end, &mut store).with_context(|| {
        format!("fetch failed for {contig}:{}-{}", start.as_u32(), end.as_u32())
    })?;

    println!("fetched {} record(s) from {contig}:{}-{}", store.len(), start.as_u32(), end.as_u32());

    // Pass 1: inspect records and plan changes without mutating the store.
    // `set_alignment` needs `&mut self`, so we collect the work first.
    let mut plan: Vec<(RecordIdx, Pos0, Vec<CigarOp>)> = Vec::new();
    for idx in store.indices() {
        let Some(rec) = store.record(idx) else { continue };
        // Skip unmapped reads — no alignment to update.
        if rec.flags.is_unmapped() || rec.tid < 0 {
            continue;
        }
        let Some((new_pos, new_cigar)) =
            propose_realignment(rec.cigar(), rec.pos, args.clip, args.min_match)
        else {
            continue;
        };
        plan.push((idx, new_pos, new_cigar));
    }

    // Pass 2: apply every planned change.
    let total_planned = plan.len();
    let mut applied = 0usize;
    for (i, (idx, new_pos, new_cigar)) in plan.iter().enumerate() {
        let before = snapshot(&store, *idx);
        store
            .set_alignment(*idx, *new_pos, new_cigar)
            .with_context(|| format!("set_alignment failed for record {idx}"))?;
        applied += 1;

        if i < args.show
            && let (Some(before), Some(after)) = (before.as_ref(), snapshot(&store, *idx))
        {
            println!(
                "  [{idx}] {:?}  {} @ {}  ->  {} @ {}",
                before.qname,
                before.cigar_str,
                before.pos.as_u32(),
                after.cigar_str,
                after.pos.as_u32(),
            );
        }
    }

    // After mutating positions, the store is no longer position-sorted.
    // Any further pileup or writing step needs a sort pass first.
    store.sort_by_pos();

    let skipped = total_planned.saturating_sub(applied);
    println!("realigned {applied} record(s); skipped {skipped}; store now sorted and ready");

    if let Some(ref out_path) = args.output {
        write_store(&store, reader.header(), out_path)?;
        println!("wrote {} record(s) to {}", store.len(), out_path.display());
    }

    Ok(())
}

/// Resolve a `RegionString` (or "first contig" fallback) against the BAM
/// header into a `(tid, start, end_exclusive, contig_name)` tuple.
fn resolve_region(
    region: Option<&RegionString>,
    header: &seqair::bam::BamHeader,
) -> anyhow::Result<(u32, Pos0, Pos0, String)> {
    match region {
        Some(r) => {
            let tid = r.chromosome.as_str().resolve_tid(header)?.as_u32();
            let len = header.target_len(tid).context("contig has no length")? as u32;
            let start = match r.start {
                Some(p) => p.to_zero_based(),
                None => Pos0::new(0).expect("0 is valid"),
            };
            let end = match r.end {
                Some(p) => p.to_zero_based(),
                None => Pos0::new(len).context("contig length out of range")?,
            };
            Ok((tid, start, end, r.chromosome.to_string()))
        }
        None => {
            let tid = 0u32;
            let len = header.target_len(tid).context("empty header")? as u32;
            let contig = header.target_name(tid).context("unknown tid")?.to_owned();
            let start = Pos0::new(0).expect("0 is valid");
            let end = Pos0::new(len).context("contig length out of range")?;
            Ok((tid, start, end, contig))
        }
    }
}

/// The proposed realignment: turn the first `clip` bases of the leading `M`
/// op into a soft clip and shift `pos` right by `clip`. Returns `None` when
/// the record does not match the rule (empty CIGAR, leading op not M, or M
/// too short).
fn propose_realignment(
    cigar: &[CigarOp],
    pos: Pos0,
    clip: u32,
    min_match: u32,
) -> Option<(Pos0, Vec<CigarOp>)> {
    let first = *cigar.first()?;
    if clip == 0 || !matches!(first.op_type(), CigarOpType::Match) {
        return None;
    }
    let first_len = first.len();
    if first_len < min_match || first_len <= clip {
        return None;
    }

    let new_match_len = first_len - clip;
    let new_pos = Pos0::new(pos.as_u32().checked_add(clip)?)?;

    // Rebuild: [clip]S + [new_match_len]M + cigar[1..]
    let mut out = Vec::with_capacity(cigar.len() + 1);
    out.push(CigarOp::new(CigarOpType::SoftClip, clip));
    out.push(CigarOp::new(CigarOpType::Match, new_match_len));
    out.extend_from_slice(&cigar[1..]);
    Some((new_pos, out))
}

struct Snapshot {
    qname: String,
    pos: Pos0,
    cigar_str: SmolStr,
}

fn snapshot(store: &RecordStore, idx: RecordIdx) -> Option<Snapshot> {
    let rec = store.record(idx)?;
    let qname = std::str::from_utf8(rec.qname()).unwrap_or("<non-utf8>").to_owned();
    Some(Snapshot { qname, pos: rec.pos, cigar_str: CigarStr(rec.cigar()).to_smolstr() })
}

fn write_store(
    store: &RecordStore,
    header: &seqair::bam::BamHeader,
    path: &std::path::Path,
) -> anyhow::Result<()> {
    let mut writer =
        BamWriterBuilder::to_path(path, header).build().context("could not create BAM writer")?;

    for i in store.indices() {
        writer
            .write_store_record(store, i)
            .with_context(|| format!("could not write record {i}"))?;
    }

    writer.finish().context("could not finalize BAM")?;
    Ok(())
}
