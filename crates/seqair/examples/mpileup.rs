#![allow(
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    clippy::indexing_slicing,
    reason = "example"
)]

use anyhow::Context;
use clap::Parser as _;
use seqair::bam::pileup::{PileupEngine, PileupOp, RefSeq};
use seqair::reader::IndexedReader;
use seqair::{
    bam::{Pos0, RecordStore},
    fasta::IndexedFastaReader,
};
use seqair_types::{Base, RegionString, SmolStr};
use std::collections::HashMap;
use std::{
    io::{BufWriter, Write},
    path::PathBuf,
    rc::Rc,
};

/// seqair mpileup — a simple pileup viewer
///
/// Usage example replicating core features of `samtools mpileup`.
/// Outputs one line per covered reference position with depth and
/// per-read base/indel information.
///
/// By default, output matches htslib's `bam_plp_auto` semantics (includes
/// trailing D/N positions in depth). Use `--samtools-compat` to match
/// `samtools mpileup -B` output (excludes reads at trailing D/N positions).
#[derive(Debug, clap::Parser)]
struct Cli {
    /// BAM/SAM/CRAM file to read
    input: PathBuf,

    /// Reference FASTA (required for CRAM, optional for BAM/SAM).
    /// Enables reference base display in column 3.
    #[clap(long, short = 'f')]
    reference: Option<PathBuf>,

    /// Region to display (e.g. "chr1:1000-2000"). Omit for whole genome.
    #[clap(long, short)]
    region: Option<RegionString>,

    /// Output file (defaults to stdout)
    #[clap(long, short)]
    out: Option<PathBuf>,

    /// Match `samtools mpileup -B` output: exclude reads at trailing
    /// deletion/ref-skip positions (positions past the last M/=/X op).
    #[clap(long)]
    samtools_compat: bool,
}

fn main() -> anyhow::Result<()> {
    let args = Cli::parse();

    let mut reader = IndexedReader::open(&args.input).context("could not open alignment file")?;
    let mut fasta = args
        .reference
        .as_ref()
        .map(|p| IndexedFastaReader::open(p))
        .transpose()
        .context("could not open reference FASTA")?;

    let mut output: Box<dyn Write> = if let Some(out) = &args.out {
        let file =
            std::fs::File::create(out).with_context(|| format!("could not create {out:?}"))?;
        Box::new(BufWriter::new(file))
    } else {
        Box::new(BufWriter::new(std::io::stdout().lock()))
    };

    let regions = if let Some(region_str) = args.region.as_ref() {
        let name = &region_str.chromosome;
        let tid =
            reader.header().tid(name).with_context(|| format!("contig '{name}' not found"))?;
        vec![(
            tid,
            *region_str.start.context("no start pos given")?,
            *region_str.end.context("no end pos given")?,
        )]
    } else {
        (0..reader.header().target_count())
            .map(|i| {
                let tid = i as u32;
                let len = reader.header().target_len(tid).unwrap_or(0);
                (tid, 0u32, len as u32)
            })
            .collect()
    };

    // Window size for streaming pileup.  1 MiB of genomic positions keeps
    // the RegionBuf well under 256 MiB even at very high coverage, while
    // being large enough that index-query overhead is negligible.
    const WINDOW: u32 = 1_000_000;

    // Pre-compute all windows so the header borrow is released before the
    // fetch loop (fetch_into needs &mut reader).
    let mut windows: Vec<(u32, SmolStr, u32, u32)> = Vec::new();
    for (tid, region_start, region_end) in &regions {
        let contig_name: SmolStr = reader.header().target_name(*tid).context("unknown tid")?.into();
        let mut ws = *region_start;
        while ws < *region_end {
            let we = (*region_end).min(ws.saturating_add(WINDOW));
            windows.push((*tid, contig_name.clone(), ws, we));
            ws = we;
        }
    }

    let mut store = Some(RecordStore::new());
    // Cache: record_idx → last query-consuming ref position (exclusive, 0-based).
    // Only used in samtools-compat mode.
    let mut query_end_cache: HashMap<u32, u32> = HashMap::new();

    let mut bases = String::new();
    let mut quals = String::new();

    for (tid, contig_name, win_start, win_end) in &windows {
        let start_pos = Pos0::new(*win_start).context("invalid start position")?;
        let end_pos = Pos0::new(*win_end).context("invalid end position")?;

        let mut s = store.take().unwrap_or_default();
        reader
            .fetch_into(*tid, start_pos, end_pos, &mut s)
            .with_context(|| format!("fetch failed for {contig_name}:{win_start}-{win_end}"))?;

        // Pre-compute query-end positions for samtools-compat mode.
        if args.samtools_compat {
            query_end_cache.clear();
            for i in 0..s.len() as u32 {
                let rec = s.record(i);
                let qend = query_end_pos(&s, i, *rec.pos);
                query_end_cache.insert(i, qend);
            }
        }

        let ref_seq = if let Some(ref mut fasta) = fasta {
            let raw = fasta
                .fetch_seq(contig_name, start_pos, end_pos)
                .with_context(|| format!("could not fetch reference for {contig_name}"))?;
            let bases = Base::from_ascii_vec(raw);
            Some(RefSeq::new(Rc::from(bases), start_pos))
        } else {
            None
        };

        let mut engine = PileupEngine::new(s, start_pos, end_pos);
        if let Some(rs) = ref_seq {
            engine.set_reference_seq(rs);
        }

        while let Some(column) = engine.pileups() {
            let pos = column.pos();
            let pos1 = *pos + 1;
            let ref_base = column.reference_base();
            let store = column.store();

            let mut depth = 0u32;
            bases.clear();
            quals.clear();

            for aln in column.alignments() {
                // In samtools-compat mode, skip alignments at trailing D/N positions.
                if args.samtools_compat
                    && aln.qpos().is_none()
                    && let Some(&qend) = query_end_cache.get(&aln.record_idx())
                    && *pos >= qend
                {
                    continue;
                }

                depth += 1;
                let rec = store.record(aln.record_idx());
                let is_reverse = aln.flags.is_reverse();

                // Read-start marker
                if pos == rec.pos {
                    bases.push('^');
                    #[allow(
                        clippy::cast_possible_truncation,
                        reason = "mapq is u8, +33 fits in u8"
                    )]
                    bases.push((aln.mapq.saturating_add(33)) as char);
                }

                format_alignment(&mut bases, &aln.op, ref_base, is_reverse, || {
                    read_inserted_bases(store, aln.record_idx(), &aln.op)
                });

                // Read-end marker
                let read_end = if args.samtools_compat {
                    query_end_cache.get(&aln.record_idx()).copied().unwrap_or(*rec.end_pos)
                } else {
                    *rec.end_pos
                };
                if *pos + 1 >= read_end {
                    bases.push('$');
                }

                match aln.qual().and_then(|q| q.get()) {
                    Some(q) => quals.push((q.saturating_add(33)) as char),
                    None => quals.push('~'),
                }
            }

            if depth > 0 {
                writeln!(
                    output,
                    "{contig_name}\t{pos1}\t{}\t{depth}\t{bases}\t{quals}",
                    ref_base as u8 as char
                )?;
            }
        }

        store = engine.take_store();
        query_end_cache.clear();
    }

    output.flush()?;
    Ok(())
}

/// Compute the "samtools end" of a read: the full reference span minus any
/// trailing D/N ops that have no query-consuming ops (M/I/S/=/X) after them.
///
/// For CIGAR `7M2D2I`: the 2I after the D means the D is NOT trailing → full rlen.
/// For CIGAR `2D7M2D`: the trailing 2D has nothing after → strip it.
fn query_end_pos(store: &RecordStore, record_idx: u32, record_pos: u32) -> u32 {
    let cigar = store.cigar(record_idx);
    let n_ops = cigar.len() / 4;

    // Walk backwards to find trailing ref-consuming, non-query-consuming ops
    // (D=2, N=3) that have no query-consuming ops after them.
    let mut trailing_ref_len = 0u32;
    for i in (0..n_ops).rev() {
        let packed = u32::from_le_bytes([
            cigar[i * 4],
            cigar[i * 4 + 1],
            cigar[i * 4 + 2],
            cigar[i * 4 + 3],
        ]);
        let op = packed & 0xF;
        let len = packed >> 4;

        match op {
            2 | 3 => trailing_ref_len += len, // D, N: trailing ref-only ops
            5 | 6 => {}                       // H, P: skip (no ref, no query)
            _ => break, // M, I, S, =, X: query-consuming op found → stop stripping
        }
    }

    // Full rlen
    let mut rlen = 0u32;
    for chunk in cigar.chunks_exact(4) {
        let packed = u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]);
        let op = packed & 0xF;
        let len = packed >> 4;
        if matches!(op, 0 | 2 | 3 | 7 | 8) {
            rlen += len;
        }
    }

    record_pos + rlen - trailing_ref_len
}

/// Read the inserted bases from a record's sequence.
fn read_inserted_bases(store: &RecordStore, record_idx: u32, op: &PileupOp) -> Option<Vec<u8>> {
    let (qpos, insert_len) = match op {
        PileupOp::Insertion { qpos, insert_len, .. } if *insert_len > 0 => (*qpos, *insert_len),
        _ => return None,
    };
    let start = qpos as usize + 1;
    let seq = store.seq(record_idx);
    let end = (start + insert_len as usize).min(seq.len());
    Some(seq.get(start..end)?.iter().map(|b| *b as u8).collect())
}

/// Format a single alignment in samtools mpileup style.
fn format_alignment(
    buf: &mut String,
    op: &PileupOp,
    ref_base: Base,
    is_reverse: bool,
    get_ins_bases: impl FnOnce() -> Option<Vec<u8>>,
) {
    use std::fmt::Write;

    fn push_base(buf: &mut String, base: Base, ref_base: Base, is_reverse: bool) {
        if base == ref_base && ref_base != Base::Unknown {
            buf.push(if is_reverse { ',' } else { '.' });
        } else {
            let ch = base as u8 as char;
            if is_reverse {
                buf.push(ch.to_ascii_lowercase());
            } else {
                buf.push(ch.to_ascii_uppercase());
            }
        }
    }

    match op {
        PileupOp::Match { base, .. } => push_base(buf, *base, ref_base, is_reverse),
        PileupOp::Insertion { base, insert_len, .. } => {
            push_base(buf, *base, ref_base, is_reverse);
            let _ = write!(buf, "+{insert_len}");
            if let Some(bases) = get_ins_bases() {
                for b in &bases {
                    let ch = *b as char;
                    buf.push(if is_reverse { ch.to_ascii_lowercase() } else { ch });
                }
            }
        }
        PileupOp::Deletion { .. } => buf.push('*'),
        PileupOp::ComplexIndel { insert_len, .. } => {
            let _ = write!(buf, "*+{insert_len}");
        }
        PileupOp::RefSkip => buf.push(if is_reverse { '<' } else { '>' }),
    }
}
