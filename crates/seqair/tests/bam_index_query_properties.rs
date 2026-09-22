//! The write → index → query round-trip, on records nobody wrote by hand.
//!
//! seqair writes the BAM and co-produces the BAI, then three parties are asked
//! which records overlap a region: seqair's own indexed reader, `samtools
//! view` through the same index, and a brute-force scan of the records that
//! were generated in the first place. All three must give the same answer.
//!
//! That one comparison covers `OwnedBamRecord` serialization, `BamWriter`,
//! `IndexBuilder`'s BAI, BGZF framing, the BAI query, `RegionBuf` and the
//! reader's overlap filter — and the brute-force scan means a bug shared by
//! seqair's writer and reader cannot hide behind htslib agreeing with it.
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
use std::process::Command;

/// Three contigs. `chr1` spans several 16 kb BAI leaf bins, so a query there
/// has to combine bins rather than read one; the others are short enough that
/// generated positions collide constantly.
const CONTIGS: [(&str, u32); 3] = [("chr1", 60_000), ("chr2", 8_000), ("chr3", 2_000)];

fn header() -> BamHeader {
    let mut text = String::from("@HD\tVN:1.6\tSO:coordinate\n");
    for (name, len) in CONTIGS {
        writeln!(text, "@SQ\tSN:{name}\tLN:{len}").expect("writing to a String cannot fail");
    }
    BamHeader::from_sam_text(&text).unwrap()
}

/// A record as the generator describes it, before it becomes BAM bytes.
#[derive(Debug, Clone)]
struct Read {
    tid: u32,
    pos: u32,
    /// `(len, op_type)`; empty for a placed-unmapped read.
    cigar: Vec<(u32, CigarOpType)>,
    /// Placed but unmapped: it has a position and is indexed there, but no
    /// alignment, so it covers exactly one reference base.
    unmapped: bool,
    seq_len: u32,
}

impl Read {
    /// Last reference position the record covers, inclusive.
    fn end_pos(&self) -> u32 {
        if self.unmapped {
            return self.pos;
        }
        let span: u32 = self
            .cigar
            .iter()
            .filter(|(_, op)| {
                matches!(
                    op,
                    CigarOpType::Match
                        | CigarOpType::Deletion
                        | CigarOpType::RefSkip
                        | CigarOpType::SeqMatch
                        | CigarOpType::SeqMismatch
                )
            })
            .map(|(len, _)| len)
            .sum();
        self.pos + span.saturating_sub(1)
    }

    fn overlaps(&self, start: u32, end: u32) -> bool {
        self.pos <= end && self.end_pos() >= start
    }
}

/// A start position on `contig_len`, drawn inside one of two clusters at
/// opposite ends of the contig. Uniform positions over 60 kb leave reads too
/// sparse to overlap each other or a query window; clustering keeps them dense
/// while still putting them in different BAI leaf bins.
fn arb_pos(tc: &TestCase, contig_len: u32, span: u32) -> u32 {
    let last = contig_len.saturating_sub(span).max(1) - 1;
    let cluster = tc.draw_silent(gs::integers::<u32>().max_value(1));
    let base = cluster * (contig_len / 2);
    let jitter = tc.draw_silent(gs::integers::<u32>().max_value(600));
    base.saturating_add(jitter).min(last)
}

#[hegel::composite]
fn arb_read(tc: &TestCase) -> Read {
    let n = |lo: u32, hi: u32| gs::integers::<u32>().min_value(lo).max_value(hi);
    let tid = tc.draw_silent(gs::integers::<u32>().max_value(CONTIGS.len() as u32 - 1));
    let contig_len = CONTIGS[tid as usize].1;

    // One in ten reads is placed but unmapped — the case with its own index
    // dispatch (`push_unmapped`) and its own end_pos rule.
    if !tc.draw_silent(gs::weighted_booleans(0.9)) {
        let pos = arb_pos(tc, contig_len, 1);
        let seq_len = tc.draw_silent(n(1, 40));
        return Read { tid, pos, cigar: Vec::new(), unmapped: true, seq_len };
    }

    // Matches at both ends, indels and skips between them, clips outside.
    let mut cigar: Vec<(u32, CigarOpType)> = Vec::new();
    if tc.draw_silent(gs::booleans()) {
        cigar.push((tc.draw_silent(n(1, 5)), CigarOpType::SoftClip));
    }
    cigar.push((tc.draw_silent(n(1, 30)), CigarOpType::Match));
    for _ in 0..tc.draw_silent(gs::integers::<usize>().max_value(2)) {
        let (op, len) = match tc.draw_silent(gs::integers::<u8>().max_value(2)) {
            0 => (CigarOpType::Insertion, tc.draw_silent(n(1, 5))),
            1 => (CigarOpType::Deletion, tc.draw_silent(n(1, 8))),
            // A spliced read: its span is its intron, which is exactly what
            // makes a window query need `reach` rather than a `pos` bound.
            _ => (CigarOpType::RefSkip, tc.draw_silent(n(1, 400))),
        };
        cigar.push((len, op));
        cigar.push((tc.draw_silent(n(1, 30)), CigarOpType::Match));
    }
    if tc.draw_silent(gs::booleans()) {
        cigar.push((tc.draw_silent(n(1, 5)), CigarOpType::SoftClip));
    }

    let span: u32 = cigar
        .iter()
        .filter(|(_, op)| {
            matches!(op, CigarOpType::Match | CigarOpType::Deletion | CigarOpType::RefSkip)
        })
        .map(|(len, _)| len)
        .sum();
    let pos = arb_pos(tc, contig_len, span);
    let seq_len = cigar
        .iter()
        .filter(|(_, op)| {
            matches!(op, CigarOpType::Match | CigarOpType::Insertion | CigarOpType::SoftClip)
        })
        .map(|(len, _)| len)
        .sum();

    Read { tid, pos, cigar, unmapped: false, seq_len }
}

/// A coordinate-sorted set of reads, which is what `BamWriter`'s index
/// co-production requires.
#[hegel::composite]
fn arb_reads(tc: &TestCase) -> Vec<Read> {
    let mut reads = tc.draw_silent(gs::vecs(arb_read()).min_size(6).max_size(40));
    reads.sort_by_key(|r| (r.tid, r.pos));
    reads
}

/// Write the reads as BAM with a co-produced BAI. Names are `r{index}`, so a
/// query's answer is a set of indices into `reads`.
fn write_bam(dir: &Path, reads: &[Read]) -> std::path::PathBuf {
    let header = header();
    let bam_path = dir.join("generated.bam");
    let mut writer =
        BamWriterBuilder::to_path(&bam_path, &header).write_index(true).build().unwrap();

    for (i, read) in reads.iter().enumerate() {
        let seq = vec![Base::A; read.seq_len as usize];
        let qual = vec![BaseQuality::from_byte(30); read.seq_len as usize];
        let cigar: Vec<CigarOp> =
            read.cigar.iter().map(|&(len, op)| CigarOp::new(op, len)).collect();
        let flags = if read.unmapped { BamFlags::from(0x4u16) } else { BamFlags::empty() };
        let rec = OwnedBamRecord::builder(
            read.tid as i32,
            Some(Pos0::new(read.pos).unwrap()),
            format!("r{i}").into_bytes(),
        )
        .flags(flags)
        .mapq(60)
        .cigar(cigar)
        .seq(seq)
        .qual(qual)
        .build()
        .unwrap();
        writer.write(&rec).unwrap();
    }

    let (_inner, index_builder) = writer.finish().unwrap();
    let bai = std::fs::File::create(bam_path.with_extension("bam.bai")).unwrap();
    index_builder.expect("write_index(true)").write_bai(bai, header.target_count()).unwrap();
    bam_path
}

/// `samtools view <bam> <region>`, as the set of read indices it returns.
fn samtools_query(bam_path: &Path, region: &str) -> Vec<usize> {
    let output = Command::new("samtools")
        .arg("view")
        .arg(bam_path)
        .arg(region)
        .output()
        .expect("samtools not found");
    assert!(
        output.status.success(),
        "samtools view {region} failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let text = std::str::from_utf8(&output.stdout).expect("samtools emits UTF-8");
    let mut ids: Vec<usize> = text
        .lines()
        .filter(|l| !l.is_empty())
        .map(|line| {
            let qname = line.split('\t').next().expect("a SAM line has fields");
            qname.trim_start_matches('r').parse().expect("generated qnames are r<index>")
        })
        .collect();
    ids.sort_unstable();
    ids
}

fn seqair_query(bam_path: &Path, tid_name: &str, start: u32, end: u32) -> Vec<usize> {
    let mut reader = IndexedBamReader::open(bam_path).expect("open");
    let tid = reader.header().tid(tid_name).expect("contig");
    let mut store = RecordStore::new();
    let span = (Pos0::new(start).unwrap()..=Pos0::new(end).unwrap()).into();
    reader.fetch_into(tid, span, &mut store).expect("fetch");
    let mut ids: Vec<usize> = store
        .indices()
        .map(|idx| {
            let qname = store.record(idx).unwrap().qname().to_vec();
            String::from_utf8(qname)
                .expect("ASCII")
                .trim_start_matches('r')
                .parse()
                .expect("generated qnames are r<index>")
        })
        .collect();
    ids.sort_unstable();
    ids
}

/// A query window as `(contig, start, end)`, both ends 0-based inclusive.
///
/// Half the windows are anchored a base or two off a generated read's own
/// start or end, which is where an off-by-one in the overlap filter or in a
/// BAI bin boundary hides. A uniform window over a 60 kb contig almost never
/// touches a read, and a suite of empty-vs-empty comparisons proves nothing.
#[hegel::composite]
fn arb_region(tc: &TestCase, reads: Vec<Read>) -> (usize, u32, u32) {
    // Weighted towards anchored: an absolute window over a 60 kb contig is
    // usually empty, and a suite of empty-vs-empty comparisons proves nothing.
    let anchored = tc.draw_silent(gs::weighted_booleans(0.7));
    let len = tc.draw_silent(gs::integers::<u32>().max_value(500));

    if anchored && !reads.is_empty() {
        let which = tc.draw_silent(gs::integers::<usize>().max_value(reads.len() - 1));
        let read = &reads[which];
        let contig_len = CONTIGS[read.tid as usize].1;
        let at_end = tc.draw_silent(gs::booleans());
        let boundary = if at_end { read.end_pos() } else { read.pos };
        let offset = tc.draw_silent(gs::integers::<i32>().min_value(-2).max_value(2));
        let start = boundary.saturating_add_signed(offset).min(contig_len - 1);
        return (read.tid as usize, start, start.saturating_add(len).min(contig_len - 1));
    }

    let contig_idx = tc.draw_silent(gs::integers::<usize>().max_value(CONTIGS.len() - 1));
    let contig_len = CONTIGS[contig_idx].1;
    let start = tc.draw_silent(gs::integers::<u32>().max_value(contig_len - 1));
    (contig_idx, start, start.saturating_add(len).min(contig_len - 1))
}

// r[verify bam.reader.overlap_filter]
// r[verify interval.span_type]
/// A span whose `last` is before its `start` names no positions, so it names
/// no records — quietly, the same as a window that happens to be empty. The
/// window is anchored on real reads so the bins the index would visit are not
/// empty; only the span is.
#[hegel::test(test_cases = 24)]
fn a_reversed_span_names_no_records(tc: TestCase) {
    let reads = tc.draw(arb_reads().print_as_debug());
    let (contig_idx, start, end) = tc.draw(arb_region(reads.clone()).print_as_debug());
    let (contig, _) = CONTIGS[contig_idx];
    // Cross the ends by at least one base; `start == end` is a real window.
    let gap = tc.draw(gs::integers::<u32>().min_value(1).max_value(3));
    let (start, last) = (end.saturating_add(gap), start);

    let dir = tempfile::tempdir().expect("tempdir");
    let bam_path = write_bam(dir.path(), &reads);
    assert_eq!(seqair_query(&bam_path, contig, start, last), Vec::<usize>::new());
}

// r[verify bam.reader.overlap_filter]
// r[verify index_builder.bai_format]
// r[verify index_builder.bai_all_refs]
// r[verify index_builder.binning]
/// seqair's reader, samtools through seqair's index, and the generated truth
/// must name the same records for the same region.
#[hegel::test(test_cases = 48)]
fn region_query_agrees_with_samtools_and_brute_force(tc: TestCase) {
    let reads = tc.draw(arb_reads().print_as_debug());
    let (contig_idx, start, end) = tc.draw(arb_region(reads.clone()).print_as_debug());
    let (contig, _) = CONTIGS[contig_idx];

    let dir = tempfile::tempdir().expect("tempdir");
    let bam_path = write_bam(dir.path(), &reads);

    let expected: Vec<usize> = reads
        .iter()
        .enumerate()
        // A placed-unmapped read is returned too: it has a position, it is
        // indexed there, and `r[bam.reader.overlap_filter]` is about positions,
        // not about having an alignment. It covers exactly `pos`.
        .filter(|(_, r)| r.tid as usize == contig_idx && r.overlaps(start, end))
        .map(|(i, _)| i)
        .collect();

    // samtools takes 1-based inclusive coordinates.
    let region = format!("{contig}:{}-{}", start + 1, end + 1);
    let from_samtools = samtools_query(&bam_path, &region);
    let from_seqair = seqair_query(&bam_path, contig, start, end);

    assert_eq!(from_seqair, expected, "seqair vs the generated reads, {region}");
    assert_eq!(from_samtools, expected, "samtools vs the generated reads, {region}");

    tc.event_value("hits", expected.len() as f64);
    if !expected.is_empty() {
        tc.event("non-empty region");
    }
    if reads.iter().any(|r| r.unmapped) {
        tc.event("has a placed-unmapped read");
    }
}
