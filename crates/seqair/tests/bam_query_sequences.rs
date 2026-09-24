//! Many queries on one reader give what each would give on its own.
//!
//! `IndexedBamReader` keeps decompressed BGZF blocks between queries
//! (`r[region_buf.block_cache]`) and starts a query where the previous one
//! first found a record (`r[bam.reader.resume]`), so a query's answer could
//! in principle depend on the queries before it. These tests run random query sequences —
//! repeats, small steps, backward jumps, switches between references — on one
//! long-lived reader and check every answer against a truth that has no
//! history: the generated reads themselves, or one whole-reference scan per
//! reference on a fresh reader.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "test code"
)]
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    reason = "test code with known small values"
)]

use hegel::prelude::*;
use seqair::bam::cigar::{CigarOp, CigarOpType};
use seqair::bam::header::BamHeader;
use seqair::bam::owned_record::OwnedBamRecord;
use seqair::bam::writer::BamWriterBuilder;
use seqair::bam::{IndexedBamReader, Pos0, RecordStore};
use seqair_types::{BamFlags, Base, BaseQuality};
use std::fmt::Write as _;
use std::path::Path;
use std::sync::OnceLock;

/// One query: reference index, first and last position (0-based inclusive).
type Query = (usize, u32, u32);

/// The next query, drawn relative to the previous one: the same window again,
/// a step either way, or anywhere on any reference.
fn draw_query(tc: &TestCase, prev: Option<Query>, contigs: &[(u32, u32)]) -> Query {
    let len = match tc.draw_silent(gs::integers::<u8>().max_value(3)) {
        0 => 0,
        1 => tc.draw_silent(gs::integers::<u32>().max_value(200)),
        2 => tc.draw_silent(gs::integers::<u32>().max_value(3_000)),
        _ => tc.draw_silent(gs::integers::<u32>().max_value(40_000)),
    };
    let (tid, start) = match (prev, tc.draw_silent(gs::integers::<u8>().max_value(3))) {
        (Some((tid, start, last)), 0) => return (tid, start, last),
        (Some((tid, start, _)), 1 | 2) => {
            let step = tc.draw_silent(gs::integers::<i32>().min_value(-20_000).max_value(20_000));
            let (lo, hi) = contigs[tid];
            (tid, start.saturating_add_signed(step).clamp(lo, hi))
        }
        _ => {
            let tid = tc.draw_silent(gs::integers::<usize>().max_value(contigs.len() - 1));
            let (lo, hi) = contigs[tid];
            (tid, tc.draw_silent(gs::integers::<u32>().min_value(lo).max_value(hi)))
        }
    };
    let hi = contigs[tid].1;
    (tid, start, start.saturating_add(len).min(hi))
}

// ── Generated BAM, brute-force truth ────────────────────────────────────

const GEN_CONTIGS: [(&str, u32); 3] = [("chr1", 200_000), ("chr2", 30_000), ("chr3", 5_000)];

#[derive(Debug, Clone)]
struct Read {
    tid: usize,
    pos: u32,
    cigar: Vec<(u32, CigarOpType)>,
    unmapped: bool,
}

impl Read {
    fn last(&self) -> u32 {
        if self.unmapped {
            return self.pos;
        }
        let span: u32 = self
            .cigar
            .iter()
            .filter(|(_, op)| {
                matches!(op, CigarOpType::Match | CigarOpType::Deletion | CigarOpType::RefSkip)
            })
            .map(|(len, _)| len)
            .sum();
        self.pos + span.saturating_sub(1)
    }

    fn seq_len(&self) -> u32 {
        self.cigar
            .iter()
            .filter(|(_, op)| {
                matches!(op, CigarOpType::Match | CigarOpType::Insertion | CigarOpType::SoftClip)
            })
            .map(|(len, _)| len)
            .sum::<u32>()
            .max(1)
    }
}

/// Enough reads to fill many BGZF blocks: they come from a drawn seed rather
/// than a draw per read, which would make every test case thousands of draws.
/// Coverage is uneven (dense patches, gaps), and a few reads carry long skips
/// or deletions so that the linear index pulls some queries far back.
fn generate_reads(tc: &TestCase) -> Vec<Read> {
    let seed = tc.draw(gs::integers::<u64>());
    let counts: Vec<usize> =
        GEN_CONTIGS.iter().map(|_| tc.draw(gs::integers::<usize>().max_value(2_500))).collect();
    let mut state = seed | 1;
    let mut next = move |bound: u32| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state % u64::from(bound.max(1))) as u32
    };

    let mut reads = Vec::new();
    for (tid, (&(_, contig_len), &n)) in GEN_CONTIGS.iter().zip(&counts).enumerate() {
        for _ in 0..n {
            // Half the reads crowd into a patch at the contig's start.
            let window = if next(2) == 0 { contig_len / 20 } else { contig_len - 600 };
            let pos = next(window);
            if next(20) == 0 {
                reads.push(Read { tid, pos, cigar: Vec::new(), unmapped: true });
                continue;
            }
            let mut cigar = vec![(20 + next(80), CigarOpType::Match)];
            match next(40) {
                0 => {
                    cigar.push((1 + next(20_000), CigarOpType::RefSkip));
                    cigar.push((20 + next(50), CigarOpType::Match));
                }
                1 => {
                    cigar.push((1 + next(300), CigarOpType::Deletion));
                    cigar.push((20 + next(50), CigarOpType::Match));
                }
                2 => {
                    cigar.push((1 + next(10), CigarOpType::Insertion));
                    cigar.push((20 + next(50), CigarOpType::Match));
                }
                _ => {}
            }
            let read = Read { tid, pos, cigar, unmapped: false };
            if read.last() < contig_len {
                reads.push(read);
            }
        }
    }
    reads.sort_by_key(|r| (r.tid, r.pos));
    reads
}

fn write_bam(dir: &Path, reads: &[Read]) -> std::path::PathBuf {
    let mut text = String::from("@HD\tVN:1.6\tSO:coordinate\n");
    for (name, len) in GEN_CONTIGS {
        writeln!(text, "@SQ\tSN:{name}\tLN:{len}").unwrap();
    }
    let header = BamHeader::from_sam_text(&text).unwrap();
    let bam_path = dir.join("generated.bam");
    let mut writer =
        BamWriterBuilder::to_path(&bam_path, &header).write_index(true).build().unwrap();
    for (i, read) in reads.iter().enumerate() {
        let n = read.seq_len() as usize;
        let flags = if read.unmapped { BamFlags::from(0x4u16) } else { BamFlags::empty() };
        let rec = OwnedBamRecord::builder(
            read.tid as i32,
            Some(Pos0::new(read.pos).unwrap()),
            format!("r{i}").into_bytes(),
        )
        .flags(flags)
        .mapq(60)
        .cigar(read.cigar.iter().map(|&(len, op)| CigarOp::new(op, len)).collect())
        .seq((0..n).map(|j| [Base::A, Base::C, Base::G, Base::T][(i + j) % 4]).collect())
        .qual((0..n).map(|j| BaseQuality::from_byte(((i * 7 + j) % 40) as u8)).collect())
        .build()
        .unwrap();
        writer.write(&rec).unwrap();
    }
    let (_inner, index_builder) = writer.finish().unwrap();
    let bai = std::fs::File::create(bam_path.with_extension("bam.bai")).unwrap();
    index_builder.unwrap().write_bai(bai, header.target_count()).unwrap();
    bam_path
}

/// Read indices named by the records in `store`, in store (file) order.
fn indices(store: &RecordStore) -> Vec<usize> {
    store
        .records()
        .map(|r| {
            let name = std::str::from_utf8(r.qname(store).unwrap()).unwrap();
            name.trim_start_matches('r').parse().unwrap()
        })
        .collect()
}

// r[verify region_buf.block_cache]
// r[verify bam.reader.resume]
// r[verify bam.reader.overlap_filter]
// r[verify bam.reader.sorted_order+2]
/// Every query of a random sequence on one reader returns exactly the
/// generated reads that overlap it, in file order.
#[hegel::test(test_cases = 64)]
fn query_sequence_matches_generated_reads(tc: TestCase) {
    let reads = generate_reads(&tc);
    let dir = tempfile::tempdir().unwrap();
    let bam_path = write_bam(dir.path(), &reads);
    let mut reader = IndexedBamReader::open(&bam_path).unwrap();
    let mut store = RecordStore::new();

    let contigs: Vec<(u32, u32)> = GEN_CONTIGS.iter().map(|&(_, len)| (0, len - 1)).collect();
    let n_queries = tc.draw(gs::integers::<usize>().min_value(1).max_value(40));
    let mut prev = None;
    for _ in 0..n_queries {
        let query @ (tid, start, last) = draw_query(&tc, prev, &contigs);
        prev = Some(query);
        let span = (Pos0::new(start).unwrap()..=Pos0::new(last).unwrap()).into();
        reader.fetch_into(tid as u32, span, &mut store).unwrap();

        let expected: Vec<usize> = reads
            .iter()
            .enumerate()
            .filter(|(_, r)| r.tid == tid && r.pos <= last && r.last() >= start)
            .map(|(i, _)| i)
            .collect();
        assert_eq!(indices(&store), expected, "query {query:?}");
    }
}

// ── tests/data/test.bam, whole-reference truth ──────────────────────────

fn test_bam_path() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.bam"))
}

/// Where each reference's reads are, 0-based inclusive.
const TEST_CONTIGS: [(&str, u32, u32); 3] = [
    ("chr19", 6_103_000, 6_143_300),
    ("2kb_3_Unmodified", 0, 2_017),
    ("bacteriophage_lambda_CpG", 0, 48_501),
];

/// What identifies a record here: name, position, last position, flags.
type Rec = (Vec<u8>, u32, u32, u16);

fn records(store: &RecordStore) -> Vec<Rec> {
    store
        .records()
        .map(|r| {
            (r.qname(store).unwrap().to_vec(), r.pos.as_u32(), r.end_pos.as_u32(), r.flags.raw())
        })
        .collect()
}

/// Each reference read once, whole, by a reader that has done nothing else.
fn whole_references() -> &'static [Vec<Rec>] {
    static CELL: OnceLock<Vec<Vec<Rec>>> = OnceLock::new();
    CELL.get_or_init(|| {
        TEST_CONTIGS
            .iter()
            .map(|&(name, _, _)| {
                let mut reader = IndexedBamReader::open(test_bam_path()).unwrap();
                let tid = reader.header().tid(name).unwrap();
                let len = reader.header().target_len(tid).unwrap() as u32;
                let span = (Pos0::new(0).unwrap()..=Pos0::new(len - 1).unwrap()).into();
                let mut store = RecordStore::new();
                reader.fetch_into(tid, span, &mut store).unwrap();
                records(&store)
            })
            .collect()
    })
}

// r[verify region_buf.block_cache]
// r[verify bam.reader.resume]
/// Every query of a random sequence on one reader of real data returns what
/// filtering a fresh whole-reference read for the window returns, in the
/// same order.
#[hegel::test(test_cases = 64)]
fn query_sequence_matches_whole_reference_scan(tc: TestCase) {
    let whole = whole_references();
    let mut reader = IndexedBamReader::open(test_bam_path()).unwrap();
    let mut store = RecordStore::new();

    let contigs: Vec<(u32, u32)> = TEST_CONTIGS.iter().map(|&(_, lo, hi)| (lo, hi)).collect();
    let n_queries = tc.draw(gs::integers::<usize>().min_value(1).max_value(40));
    let mut prev = None;
    for _ in 0..n_queries {
        let query @ (idx, start, last) = draw_query(&tc, prev, &contigs);
        prev = Some(query);
        let tid = reader.header().tid(TEST_CONTIGS[idx].0).unwrap();
        let span = (Pos0::new(start).unwrap()..=Pos0::new(last).unwrap()).into();
        reader.fetch_into(tid, span, &mut store).unwrap();

        let expected: Vec<Rec> =
            whole[idx].iter().filter(|r| r.1 <= last && r.2 >= start).cloned().collect();
        assert_eq!(records(&store), expected, "query {query:?}");
    }
}
