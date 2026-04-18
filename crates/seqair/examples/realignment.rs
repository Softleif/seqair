#![allow(
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    clippy::indexing_slicing,
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
//! `NM` (match) op of length ≥ `min_match`, convert the first `clip` bases
//! of that leading match into a soft clip and shift `pos` right by `clip`.
//! Example: `10M` at pos 100 → `2S8M` at pos 102 (when `clip = 2`).

use anyhow::Context;
use clap::Parser as _;
use seqair::bam::{Pos0, RecordStore};
use seqair::reader::IndexedReader;
use std::path::PathBuf;

/// seqair realignment — a toy local-realignment example.
#[derive(Debug, clap::Parser)]
struct Cli {
    /// BAM/SAM file to read (must be indexed).
    input: PathBuf,

    /// Region to process (e.g. "chr1:1000-2000"). Omit for the first contig.
    #[clap(long, short)]
    region: Option<String>,

    /// Soft-clip this many leading match bases per eligible record.
    #[clap(long, default_value_t = 2)]
    clip: u32,

    /// Only realign reads whose leading M op is at least this long.
    #[clap(long, default_value_t = 4)]
    min_match: u32,

    /// Print a before/after line for up to this many modified records.
    #[clap(long, default_value_t = 5)]
    show: usize,
}

fn main() -> anyhow::Result<()> {
    let args = Cli::parse();

    let mut reader = IndexedReader::open(&args.input).context("could not open alignment file")?;

    let (tid, start, end) = if let Some(ref r) = args.region {
        parse_region(r, reader.header())?
    } else {
        let tid = 0u32;
        let len = reader.header().target_len(tid).context("empty header")? as u32;
        (tid, 0, len)
    };
    let contig = reader.header().target_name(tid).context("unknown tid")?.to_owned();

    let start_pos = Pos0::new(start).context("invalid start position")?;
    let end_pos = Pos0::new(end).context("invalid end position")?;

    let mut store = RecordStore::new();
    reader
        .fetch_into(tid, start_pos, end_pos, &mut store)
        .with_context(|| format!("fetch failed for {contig}:{start}-{end}"))?;

    println!("fetched {} record(s) from {contig}:{start}-{end}", store.len());

    // Pass 1: inspect records and plan changes without mutating the store.
    // `set_alignment` needs `&mut self`, so we collect the work first.
    let mut plan: Vec<(u32, u32, Vec<u8>)> = Vec::new();
    for idx in 0..store.len() as u32 {
        let rec = store.record(idx);
        // Skip unmapped reads — no alignment to update.
        if rec.flags.is_unmapped() || rec.tid < 0 {
            continue;
        }
        let Some((new_pos, new_cigar)) =
            propose_realignment(store.cigar(idx), *rec.pos, args.clip, args.min_match)
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
        let Some(new_pos_typed) = Pos0::new(*new_pos) else {
            continue;
        };
        store
            .set_alignment(*idx, new_pos_typed, new_cigar)
            .with_context(|| format!("set_alignment failed for record {idx}"))?;
        applied += 1;

        if i < args.show {
            let after = snapshot(&store, *idx);
            println!(
                "  [{}] {:?}  {} @ {}  ->  {} @ {}",
                idx, before.qname, before.cigar_str, before.pos, after.cigar_str, after.pos,
            );
        }
    }

    // After mutating positions, the store is no longer position-sorted.
    // Any further pileup or writing step needs a sort pass first.
    store.sort_by_pos();

    let skipped = total_planned.saturating_sub(applied);
    println!("realigned {applied} record(s); skipped {skipped}; store now sorted and ready");

    Ok(())
}

/// The proposed realignment: turn the first `clip` bases of the leading `M`
/// op into a soft clip and shift `pos` right by `clip`. Returns `None` when
/// the record does not match the rule (empty CIGAR, leading op not M, or M
/// too short).
fn propose_realignment(
    cigar: &[u8],
    pos: u32,
    clip: u32,
    min_match: u32,
) -> Option<(u32, Vec<u8>)> {
    if clip == 0 || cigar.len() < 4 {
        return None;
    }

    let first = u32::from_le_bytes([cigar[0], cigar[1], cigar[2], cigar[3]]);
    let first_op = first & 0xF;
    let first_len = first >> 4;

    // Only rewrite a plain leading M op; leave soft-clipped/insertion-led
    // reads alone so this stays a narrow, local edit.
    const CIGAR_M: u32 = 0;
    const CIGAR_S: u32 = 4;
    if first_op != CIGAR_M || first_len < min_match || first_len <= clip {
        return None;
    }

    let new_match_len = first_len - clip;
    let new_pos = pos.checked_add(clip)?;

    // Rebuild: [clip]S + [new_match_len]M + cigar[4..]
    let mut out = Vec::with_capacity(cigar.len() + 4);
    out.extend_from_slice(&pack_op(clip, CIGAR_S).to_le_bytes());
    out.extend_from_slice(&pack_op(new_match_len, CIGAR_M).to_le_bytes());
    out.extend_from_slice(&cigar[4..]);
    Some((new_pos, out))
}

fn pack_op(len: u32, op: u32) -> u32 {
    (len << 4) | op
}

struct Snapshot {
    qname: String,
    pos: u32,
    cigar_str: String,
}

fn snapshot(store: &RecordStore, idx: u32) -> Snapshot {
    let rec = store.record(idx);
    let qname = String::from_utf8_lossy(store.qname(idx)).into_owned();
    Snapshot { qname, pos: *rec.pos, cigar_str: fmt_cigar(store.cigar(idx)) }
}

fn fmt_cigar(packed: &[u8]) -> String {
    const OPS: &[u8] = b"MIDNSHP=X";
    let mut s = String::new();
    for chunk in packed.chunks_exact(4) {
        let w = u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]);
        let op = (w & 0xF) as usize;
        let len = w >> 4;
        use std::fmt::Write as _;
        let _ = write!(s, "{len}");
        s.push(*OPS.get(op).unwrap_or(&b'?') as char);
    }
    s
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
