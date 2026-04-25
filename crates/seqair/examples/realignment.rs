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
use seqair::bam::{BamWriter, Pos0, RecordStore};
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

    /// Write realigned records to this BAM file.
    #[clap(long, short)]
    output: Option<PathBuf>,
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
    let mut plan: Vec<(u32, u32, Vec<seqair::bam::CigarOp>)> = Vec::new();
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

    if let Some(ref out_path) = args.output {
        write_store(&store, reader.header(), out_path)?;
        println!("wrote {} record(s) to {}", store.len(), out_path.display());
    }

    Ok(())
}

/// The proposed realignment: turn the first `clip` bases of the leading `M`
/// op into a soft clip and shift `pos` right by `clip`. Returns `None` when
/// the record does not match the rule (empty CIGAR, leading op not M, or M
/// too short).
fn propose_realignment(
    cigar: &[seqair::bam::CigarOp],
    pos: u32,
    clip: u32,
    min_match: u32,
) -> Option<(u32, Vec<seqair::bam::CigarOp>)> {
    use seqair::bam::CigarOp;
    use seqair::bam::cigar::CigarOpType;

    let first = *cigar.first()?;
    if clip == 0 || !matches!(first.op_type(), CigarOpType::Match) {
        return None;
    }
    let first_len = first.len();
    if first_len < min_match || first_len <= clip {
        return None;
    }

    let new_match_len = first_len - clip;
    let new_pos = pos.checked_add(clip)?;

    // Rebuild: [clip]S + [new_match_len]M + cigar[1..]
    let mut out = Vec::with_capacity(cigar.len() + 1);
    out.push(CigarOp::new(CigarOpType::SoftClip, clip));
    out.push(CigarOp::new(CigarOpType::Match, new_match_len));
    out.extend_from_slice(&cigar[1..]);
    Some((new_pos, out))
}

struct Snapshot {
    qname: String,
    pos: u32,
    cigar_str: String,
}

fn snapshot(store: &RecordStore, idx: u32) -> Snapshot {
    let rec = store.record(idx);
    let qname = std::str::from_utf8(store.qname(idx)).unwrap_or("<non-utf8>").to_owned();
    Snapshot { qname, pos: *rec.pos, cigar_str: fmt_cigar(store.cigar(idx)) }
}

fn fmt_cigar(ops: &[seqair::bam::CigarOp]) -> String {
    use seqair::bam::cigar::CigarOpType;
    let mut s = String::new();
    for op in ops {
        use std::fmt::Write as _;
        let _ = write!(s, "{}", op.len());
        let c = match op.op_type() {
            CigarOpType::Match => 'M',
            CigarOpType::Insertion => 'I',
            CigarOpType::Deletion => 'D',
            CigarOpType::RefSkip => 'N',
            CigarOpType::SoftClip => 'S',
            CigarOpType::HardClip => 'H',
            CigarOpType::Padding => 'P',
            CigarOpType::SeqMatch => '=',
            CigarOpType::SeqMismatch => 'X',
            CigarOpType::Unknown(_) => '?',
        };
        s.push(c);
    }
    s
}

fn write_store(
    store: &RecordStore,
    header: &seqair::bam::BamHeader,
    path: &std::path::Path,
) -> anyhow::Result<()> {
    let mut writer =
        BamWriter::from_path(path, header, false).context("could not create BAM writer")?;

    for i in 0..store.len() as u32 {
        writer
            .write_store_record(store, i)
            .with_context(|| format!("could not write record {i}"))?;
    }

    writer.finish().context("could not finalize BAM")?;
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
