#![allow(
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    clippy::indexing_slicing,
    reason = "example"
)]

use anyhow::Context;
use clap::Parser as _;
use seqair::Readers;
use seqair::bam::RecordStore;
use seqair::bam::pileup::PileupOp;
use seqair::reader::{Segment, SegmentOptions};
use seqair_types::{Base, RegionString};
use std::collections::HashMap;
use std::num::NonZeroU32;
use std::{
    io::{BufWriter, Write},
    path::PathBuf,
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
    /// BAM/SAM/CRAM file to read.
    input: PathBuf,

    /// Reference FASTA. Required for CRAM, used to render the reference
    /// base column for BAM/SAM. Not optional under the unified `Readers`
    /// API.
    #[clap(long, short = 'f')]
    reference: PathBuf,

    /// Region to display. Accepts `chr`, `chr:start`, or `chr:start-end`
    /// (1-based, inclusive). Omit to scan every contig in the header.
    #[clap(long, short)]
    region: Option<RegionString>,

    /// Output file (defaults to stdout).
    #[clap(long, short)]
    out: Option<PathBuf>,

    /// Match `samtools mpileup -B` output: exclude reads at trailing
    /// deletion/ref-skip positions (positions past the last M/=/X op).
    #[clap(long)]
    samtools_compat: bool,
}

/// Tile size for the segmenter. 1 Mbp keeps peak `RecordStore` memory
/// well-bounded while keeping per-tile fetch overhead negligible.
const MAX_TILE_LEN: u32 = 1_000_000;

fn main() -> anyhow::Result<()> {
    let args = Cli::parse();

    let mut readers =
        Readers::open(&args.input, &args.reference).context("could not open alignment file")?;

    let mut output: Box<dyn Write> = if let Some(out) = &args.out {
        let file =
            std::fs::File::create(out).with_context(|| format!("could not create {out:?}"))?;
        Box::new(BufWriter::new(file))
    } else {
        Box::new(BufWriter::new(std::io::stdout().lock()))
    };

    let max_len = NonZeroU32::new(MAX_TILE_LEN).expect("non-zero literal");
    let opts = SegmentOptions::new(max_len);

    // Plan segments: either a parsed region (whole contig / start-only /
    // start-end), or `()` for a whole-genome scan. The iterator borrows
    // `&readers` only, so we collect to release the header borrow before
    // the &mut self pileup loop.
    let plan: Vec<Segment> = if let Some(region) = args.region.as_ref() {
        readers.segments(region, opts).context("could not plan region")?.collect()
    } else {
        readers.segments((), opts).context("could not plan whole-genome scan")?.collect()
    };

    // Cache: record_idx → last query-consuming ref position (exclusive,
    // 0-based). Only used in samtools-compat mode. Cleared per-segment
    // because record_idx is only unique within a single fetch.
    let mut query_end_cache: HashMap<u32, u32> = HashMap::new();

    let mut bases = String::new();
    let mut quals = String::new();

    for segment in &plan {
        let mut engine = readers
            .pileup(segment)
            .with_context(|| format!("pileup failed for {}", segment.contig()))?;

        // Pre-compute query-end positions for samtools-compat mode using
        // the engine's store.
        if args.samtools_compat {
            query_end_cache.clear();
            let store = engine.store();
            for i in 0..store.len() as u32 {
                let rec = store.record(i);
                let qend = query_end_pos(store, i, *rec.pos);
                query_end_cache.insert(i, qend);
            }
        }

        let contig_name = segment.contig().clone();

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
        // `engine` drops here; the RecordStore goes back into `readers`.
    }

    output.flush()?;
    Ok(())
}

/// Compute the "samtools end" of a read: the full reference span minus any
/// trailing D/N ops that have no query-consuming ops (M/I/S/=/X) after them.
///
/// For CIGAR `7M2D2I`: the 2I after the D means the D is NOT trailing → full rlen.
/// For CIGAR `2D7M2D`: the trailing 2D has nothing after → strip it.
fn query_end_pos<U>(store: &RecordStore<U>, record_idx: u32, record_pos: u32) -> u32 {
    let cigar = store.cigar(record_idx);

    // Walk backwards to find trailing ref-consuming, non-query-consuming ops
    // (D=2, N=3) that have no query-consuming ops after them.
    let mut trailing_ref_len = 0u32;
    for op in cigar.iter().rev() {
        match op.op_code() {
            2 | 3 => trailing_ref_len += op.len(), // D, N: trailing ref-only ops
            5 | 6 => {}                            // H, P: skip (no ref, no query)
            _ => break, // M, I, S, =, X: query-consuming op found → stop stripping
        }
    }

    // Full rlen
    let mut rlen = 0u32;
    for op in cigar {
        if matches!(op.op_code(), 0 | 2 | 3 | 7 | 8) {
            rlen += op.len();
        }
    }

    record_pos + rlen - trailing_ref_len
}

/// Read the inserted bases from a record's sequence.
fn read_inserted_bases<U>(
    store: &RecordStore<U>,
    record_idx: u32,
    op: &PileupOp,
) -> Option<Vec<u8>> {
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
