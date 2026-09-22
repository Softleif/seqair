//! CRAM and BAM must decode to the same records, on data nobody wrote by hand.
//!
//! CRAM is the most involved decoder in the crate — sequences are stored as
//! edits against a reference, spread across data series that each have their
//! own codec — and its tests all run on the same handful of checked-in files.
//! Here samtools writes the *same* generated reads as both BAM and CRAM, and
//! seqair reads both: the BAM path is the oracle, and the generated reads are
//! the oracle for that in turn, so a misreading shared by neither can hide.
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
use seqair::bam::{Pos0, RecordStore};
use seqair::reader::Readers;
use std::fmt::Write as _;
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

const CONTIG: &str = "ref1";
const CONTIG_LEN: u32 = 1_000;
/// Reads start inside this many bases of the contig start, so they overlap.
const START_WINDOW: u32 = 300;
/// The Sanger quality range a generated QUAL string is drawn from.
const QUAL_CHARS: &str = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI";

fn bases(len: usize) -> impl PrintableGenerator<String> {
    gs::text().alphabet("ACGT").min_size(len).max_size(len)
}

/// One generated read, in the form the SAM line and the assertions both need.
#[derive(Debug, Clone)]
struct Read {
    pos: u32,
    cigar: Vec<(u32, char)>,
    seq: String,
    qual: String,
    flags: u32,
    mapq: u32,
}

impl Read {
    fn cigar_text(&self) -> String {
        let mut out = String::new();
        for (len, op) in &self.cigar {
            write!(out, "{len}{op}").expect("writing to a String cannot fail");
        }
        out
    }
}

fn run(cmd: &mut Command, what: &str) {
    let status = cmd
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .unwrap_or_else(|e| panic!("{what}: {e}"));
    assert!(status.success(), "{what} failed");
}

/// Write the reference and index it, so samtools can compute the M5 tags the
/// CRAM header carries and seqair can verify them on read.
fn write_reference(dir: &Path, bases: &str) -> PathBuf {
    let fasta = dir.join("ref.fa");
    let mut f = std::fs::File::create(&fasta).expect("create FASTA");
    writeln!(f, ">{CONTIG}").expect("write");
    for chunk in bases.as_bytes().chunks(60) {
        f.write_all(chunk).expect("write");
        f.write_all(b"\n").expect("write");
    }
    drop(f);
    run(Command::new("samtools").arg("faidx").arg(&fasta), "samtools faidx");
    fasta
}

/// Sort the SAM into a BAM, then transcode that BAM to CRAM against the same
/// reference. Both are indexed. Returns `(bam, cram)`.
fn write_bam_and_cram(dir: &Path, sam: &str, fasta: &Path) -> (PathBuf, PathBuf) {
    let sam_path = dir.join("generated.sam");
    std::fs::File::create(&sam_path)
        .expect("create SAM")
        .write_all(sam.as_bytes())
        .expect("write SAM");

    let bam = dir.join("generated.bam");
    run(Command::new("samtools").args(["sort", "-o"]).arg(&bam).arg(&sam_path), "samtools sort");
    run(Command::new("samtools").arg("index").arg(&bam), "samtools index (bam)");

    let cram = dir.join("generated.cram");
    run(
        Command::new("samtools")
            .args(["view", "-C", "-T"])
            .arg(fasta)
            .arg("-o")
            .arg(&cram)
            .arg(&bam),
        "samtools view -C",
    );
    run(Command::new("samtools").arg("index").arg(&cram), "samtools index (cram)");

    (bam, cram)
}

/// A decoded record, in the fields both formats must agree on.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Decoded {
    qname: Vec<u8>,
    pos: i64,
    flags: u16,
    mapq: u8,
    cigar: String,
    seq: Vec<u8>,
    /// Phred values, not the ASCII SAM writes.
    qual: Vec<u8>,
}

/// Every record seqair reads from `path`.
fn read_all(path: &Path, fasta: &Path) -> Vec<Decoded> {
    let mut readers = Readers::open(path, fasta).expect("open");
    let tid = readers.header().tid(CONTIG).expect("contig");
    let mut store = RecordStore::new();
    readers
        .fetch_into(
            tid,
            (Pos0::new(0).unwrap()..=Pos0::new(CONTIG_LEN - 1).unwrap()).into(),
            &mut store,
        )
        .expect("fetch");

    store
        .indices()
        .map(|idx| {
            let rec = store.record(idx).unwrap();
            let cigar = rec.cigar().iter().map(ToString::to_string).collect::<String>();
            let seq: Vec<u8> = rec.seq().iter().map(|b| **b).collect();
            let qual = seqair_types::BaseQuality::slice_to_bytes(rec.qual()).to_vec();
            Decoded {
                qname: rec.qname().to_vec(),
                pos: rec.pos.as_i64(),
                flags: rec.flags.raw(),
                mapq: rec.mapq,
                cigar,
                seq,
                qual,
            }
        })
        .collect()
}

/// A reference sequence: the bases the reads are drawn against, so most
/// aligned bases match and CRAM stores them as matches rather than as a
/// substitution at every position.
#[hegel::composite]
fn arb_reference(tc: &TestCase) -> String {
    tc.draw_silent(bases(CONTIG_LEN as usize))
}

/// A CIGAR valid by construction: matches at both ends of the aligned body,
/// insertions and deletions between them, soft clips outside. No `N` or `P` —
/// this test is about sequence reconstruction, and the pileup parity test
/// already covers those ops.
#[hegel::composite]
fn arb_cigar(tc: &TestCase) -> Vec<(u32, char)> {
    let n = |lo: u32, hi: u32| gs::integers::<u32>().min_value(lo).max_value(hi);
    let mut body = vec![(tc.draw_silent(n(5, 40)), 'M')];
    for _ in 0..tc.draw_silent(gs::integers::<usize>().max_value(2)) {
        let op = tc.draw_silent(gs::sampled_from(&['I', 'D']));
        body.push((tc.draw_silent(n(1, 6)), op));
        body.push((tc.draw_silent(n(5, 40)), 'M'));
    }

    let mut ops = Vec::new();
    if tc.draw_silent(gs::booleans()) {
        ops.push((tc.draw_silent(n(1, 8)), 'S'));
    }
    ops.extend(body);
    if tc.draw_silent(gs::booleans()) {
        ops.push((tc.draw_silent(n(1, 8)), 'S'));
    }
    ops
}

/// A read whose aligned bases start out equal to the reference, with a handful
/// of positions overwritten. CRAM stores an aligned base as a match or as a
/// substitution depending on that comparison, so a read that differs
/// everywhere would only ever exercise the substitution path — and one that
/// differs nowhere only the match path.
fn arb_read(tc: &TestCase, reference: &str) -> Read {
    let ops = tc.draw_silent(arb_cigar());
    let span: u32 = ops.iter().filter(|(_, op)| matches!(op, 'M' | 'D')).map(|(l, _)| l).sum();
    let last_start = CONTIG_LEN.saturating_sub(span).clamp(1, START_WINDOW);
    let pos = tc.draw_silent(gs::integers::<u32>().min_value(1).max_value(last_start));

    let ref_bytes = reference.as_bytes();
    let mut seq: Vec<u8> = Vec::new();
    let mut ref_cursor = (pos - 1) as usize; // 0-based
    for &(len, op) in &ops {
        match op {
            'M' => {
                seq.extend_from_slice(&ref_bytes[ref_cursor..ref_cursor + len as usize]);
                ref_cursor += len as usize;
            }
            'D' => ref_cursor += len as usize,
            'I' | 'S' => seq.extend_from_slice(tc.draw_silent(bases(len as usize)).as_bytes()),
            _ => unreachable!("arb_cigar emits only M/I/D/S"),
        }
    }

    for _ in 0..tc.draw_silent(gs::integers::<usize>().max_value(4)) {
        let at = tc.draw_silent(gs::integers::<usize>().max_value(seq.len() - 1));
        seq[at] = tc.draw_silent(gs::sampled_from(b"ACGT"));
    }

    let qual =
        tc.draw_silent(gs::text().alphabet(QUAL_CHARS).min_size(seq.len()).max_size(seq.len()));
    let flags = tc.draw_silent(gs::sampled_from(&[0u32, 16]));
    let mapq = tc.draw_silent(gs::integers::<u32>().max_value(60));

    Read { pos, cigar: ops, seq: String::from_utf8(seq).expect("ACGT is ASCII"), qual, flags, mapq }
}

/// A reference plus the reads aligned to it.
#[hegel::composite]
fn arb_alignment(tc: &TestCase) -> (String, Vec<Read>) {
    let reference = tc.draw_silent(arb_reference());
    let n_reads = tc.draw_silent(gs::integers::<usize>().min_value(1).max_value(8));
    let reads = (0..n_reads).map(|_| arb_read(tc, &reference)).collect();
    (reference, reads)
}

fn sam_text(reads: &[Read]) -> String {
    let mut sam = format!("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:{CONTIG}\tLN:{CONTIG_LEN}\n");
    for (i, r) in reads.iter().enumerate() {
        writeln!(
            sam,
            "r{i}\t{}\t{CONTIG}\t{}\t{}\t{}\t*\t0\t0\t{}\t{}",
            r.flags,
            r.pos,
            r.mapq,
            r.cigar_text(),
            r.seq,
            r.qual
        )
        .expect("writing to a String cannot fail");
    }
    sam
}

// r[verify unified.fetch_equivalence]
// r[verify cram.record.sequence]
/// The same reads, written by samtools as BAM and as CRAM, must come back
/// from seqair identical in every field — and equal to what was generated.
#[hegel::test(test_cases = 32)]
fn cram_and_bam_decode_to_the_same_records(tc: TestCase) {
    let (reference, reads) = tc.draw(arb_alignment().print_as_debug());

    let dir = tempfile::tempdir().expect("tempdir");
    let fasta = write_reference(dir.path(), &reference);
    let (bam, cram) = write_bam_and_cram(dir.path(), &sam_text(&reads), &fasta);

    let from_bam = read_all(&bam, &fasta);
    let from_cram = read_all(&cram, &fasta);

    assert_eq!(from_cram.len(), reads.len(), "CRAM record count");
    assert_eq!(from_cram, from_bam, "CRAM and BAM disagree");

    // And both agree with the reads that were generated. `samtools sort`
    // reorders them, so match on the qname the generator assigned.
    for got in &from_cram {
        let name = std::str::from_utf8(&got.qname).expect("ASCII");
        let i: usize = name.trim_start_matches('r').parse().expect("qnames are r<index>");
        let want = &reads[i];
        assert_eq!(got.pos, i64::from(want.pos) - 1, "{name}: pos");
        assert_eq!(u32::from(got.flags), want.flags, "{name}: flags");
        assert_eq!(u32::from(got.mapq), want.mapq, "{name}: mapq");
        assert_eq!(got.cigar, want.cigar_text(), "{name}: cigar");
        assert_eq!(got.seq, want.seq.as_bytes(), "{name}: seq");
        // SAM writes QUAL as `phred + 33`; seqair stores the Phred value.
        let want_qual: Vec<u8> = want.qual.bytes().map(|b| b - 33).collect();
        assert_eq!(got.qual, want_qual, "{name}: qual");
    }

    tc.event_value("reads", reads.len() as f64);
    if reads.iter().any(|r| r.cigar.iter().any(|(_, op)| *op == 'D')) {
        tc.event("has a deletion");
    }
    if reads.iter().any(|r| r.cigar.iter().any(|(_, op)| *op == 'I')) {
        tc.event("has an insertion");
    }
    if reads.iter().any(|r| r.cigar.iter().any(|(_, op)| *op == 'S')) {
        tc.event("has a soft clip");
    }
}
