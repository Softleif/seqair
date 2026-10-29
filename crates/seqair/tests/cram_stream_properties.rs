//! The CRAM stream decoder and the indexed reader on generated samples
//! written by samtools under random writer options, against `samtools view`.
//!
//! The samples hold the records CRAM stores least verbatim: mapped records
//! with SEQ and QUAL `*` (bwa's secondaries), placed-unmapped records and
//! unplaced ones. Bases are ACGT and read names are kept (`RN`), so the
//! stream's known gaps with htslib (`cram_stream.rs`'s `KNOWN_MISMATCHES`)
//! cannot show up, and the comparison is byte for byte.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    reason = "test code"
)]

#[path = "../src/pinned.rs"]
mod pinned;

use hegel::TestCase;
use hegel::generators as gs;
use seqair::bam::stream::BamStreamReader;
use seqair::bam::{RawBamRecord, RecordStore};
use seqair::cram::reader::IndexedCramReader;
use seqair::cram::stream::CramStreamReader;
use seqair::sam::format::format_sam_record;
use seqair_types::Pos0;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::process::Command;

const CONTIG_LEN: usize = 300;
const QUAL_CHARS: &str = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI";
/// `reg2bin(-1, 0)`: the bin the SAM spec (§4.2.1) gives a record with no position.
const UNPLACED_BIN: u16 = 4680;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kind {
    Mapped,
    /// Mapped, flag 0x100, SEQ and QUAL `*`.
    SecondaryNoSeq,
    PlacedUnmapped,
    Unplaced,
}

#[derive(Debug, Clone)]
struct Line {
    kind: Kind,
    flag: u16,
    contig: Option<usize>,
    /// 1-based; 0 when unplaced.
    pos: usize,
    mapq: u8,
    cigar: String,
    seq: String,
    qual: String,
}

#[derive(Debug, Clone)]
struct Sample {
    contigs: Vec<String>,
    lines: Vec<Line>,
}

fn n(lo: usize, hi: usize) -> gs::IntegerGenerator<usize> {
    gs::integers::<usize>().min_value(lo).max_value(hi)
}

fn text(tc: &TestCase, alphabet: &str, len: usize) -> String {
    tc.draw_silent(gs::text().alphabet(alphabet).min_size(len).max_size(len))
}

/// A mapped line: soft clips around a match, maybe split by one indel, its
/// aligned bases the reference's with a few substitutions.
fn mapped(tc: &TestCase, contig: usize, reference: &str) -> Line {
    let mut ops: Vec<(usize, char)> = Vec::new();
    if tc.draw_silent(gs::booleans()) {
        ops.push((tc.draw_silent(n(1, 5)), 'S'));
    }
    ops.push((tc.draw_silent(n(5, 30)), 'M'));
    if tc.draw_silent(gs::booleans()) {
        ops.push((tc.draw_silent(n(1, 4)), tc.draw_silent(gs::sampled_from(vec!['I', 'D']))));
        ops.push((tc.draw_silent(n(5, 30)), 'M'));
    }
    if tc.draw_silent(gs::booleans()) {
        ops.push((tc.draw_silent(n(1, 5)), 'S'));
    }
    let span: usize = ops.iter().filter(|(_, o)| matches!(o, 'M' | 'D')).map(|(l, _)| l).sum();
    let pos = tc.draw_silent(n(1, CONTIG_LEN - span + 1));
    let mut seq = String::new();
    let mut at = pos - 1;
    for &(len, op) in &ops {
        match op {
            'M' => {
                seq.push_str(&reference[at..at + len]);
                at += len;
            }
            'D' => at += len,
            _ => seq.push_str(&text(tc, "ACGT", len)),
        }
    }
    let mut seq = seq.into_bytes();
    for _ in 0..tc.draw_silent(n(0, 3)) {
        let i = tc.draw_silent(n(0, seq.len() - 1));
        seq[i] = tc.draw_silent(gs::sampled_from(b"ACGT".to_vec()));
    }
    let mut cigar = String::new();
    for (len, op) in &ops {
        write!(cigar, "{len}{op}").unwrap();
    }
    let qual = text(tc, QUAL_CHARS, seq.len());
    Line {
        kind: Kind::Mapped,
        flag: tc.draw_silent(gs::sampled_from(vec![0u16, 16])),
        contig: Some(contig),
        pos,
        mapq: tc.draw_silent(gs::integers::<u8>().max_value(60)),
        cigar,
        seq: String::from_utf8(seq).unwrap(),
        qual,
    }
}

fn line(tc: &TestCase, kind: Kind, contigs: &[String]) -> Line {
    let contig = tc.draw_silent(n(0, contigs.len() - 1));
    match kind {
        Kind::Mapped => mapped(tc, contig, &contigs[contig]),
        Kind::SecondaryNoSeq => {
            let primary = mapped(tc, contig, &contigs[contig]);
            Line {
                kind,
                flag: primary.flag | 0x100,
                seq: "*".to_owned(),
                qual: "*".to_owned(),
                ..primary
            }
        }
        Kind::PlacedUnmapped | Kind::Unplaced => {
            let len = tc.draw_silent(n(10, 30));
            let placed = kind == Kind::PlacedUnmapped;
            Line {
                kind,
                flag: 4,
                contig: placed.then_some(contig),
                pos: if placed { tc.draw_silent(n(1, CONTIG_LEN)) } else { 0 },
                mapq: 0,
                cigar: "*".to_owned(),
                seq: text(tc, "ACGT", len),
                qual: text(tc, QUAL_CHARS, len),
            }
        }
    }
}

/// At least `min` lines of each kind in `required`, and some of any kind.
fn arb_sample(tc: &TestCase, required: &[Kind]) -> Sample {
    let contigs: Vec<String> =
        (0..tc.draw_silent(n(1, 3))).map(|_| text(tc, "ACGT", CONTIG_LEN)).collect();
    let mut kinds: Vec<Kind> = required.to_vec();
    for _ in 0..tc.draw_silent(n(0, 12)) {
        kinds.push(tc.draw_silent(gs::sampled_from(vec![
            Kind::Mapped,
            Kind::Mapped,
            Kind::SecondaryNoSeq,
            Kind::PlacedUnmapped,
            Kind::Unplaced,
        ])));
    }
    let lines = kinds.into_iter().map(|k| line(tc, k, &contigs)).collect();
    Sample { contigs, lines }
}

#[hegel::composite]
fn arb_secondary_sample(tc: &TestCase) -> Sample {
    arb_sample(tc, &[Kind::SecondaryNoSeq, Kind::Mapped])
}

#[hegel::composite]
fn arb_unplaced_sample(tc: &TestCase) -> Sample {
    arb_sample(tc, &[Kind::Unplaced, Kind::PlacedUnmapped, Kind::Mapped])
}

/// `samtools view -C` options: version, slicing, multi-reference slices and
/// an embedded reference.
#[hegel::composite]
fn arb_opts(tc: &TestCase) -> Vec<String> {
    let mut opts = vec![
        format!("version={}", tc.draw(gs::sampled_from(vec!["3.0", "3.1"]))),
        format!("seqs_per_slice={}", tc.draw(gs::sampled_from(vec![1, 2, 5, 10_000]))),
        format!("slices_per_container={}", tc.draw(gs::sampled_from(vec![1, 2]))),
    ];
    if let Some(multi) = tc.draw(gs::sampled_from(vec![None, Some(0), Some(1)])) {
        opts.push(format!("multi_seq_per_slice={multi}"));
    }
    if tc.draw(gs::booleans()) {
        opts.push("embed_ref=1".to_owned());
    }
    opts
}

fn run(cmd: &mut Command) -> Vec<u8> {
    let out = cmd.output().expect("samtools on PATH");
    assert!(out.status.success(), "{cmd:?}: {}", String::from_utf8_lossy(&out.stderr));
    out.stdout
}

/// The sample as a coordinate-sorted, indexed CRAM under `opts`, and its
/// reference.
fn write_cram(dir: &Path, sample: &Sample, opts: &[String]) -> (PathBuf, PathBuf) {
    let fasta = dir.join("ref.fa");
    let mut fa = String::new();
    for (i, contig) in sample.contigs.iter().enumerate() {
        writeln!(fa, ">c{i}\n{contig}").unwrap();
    }
    std::fs::write(&fasta, fa).unwrap();
    run(Command::new("samtools").arg("faidx").arg(&fasta));

    let mut sam = String::from("@HD\tVN:1.6\tSO:unsorted\n");
    for i in 0..sample.contigs.len() {
        writeln!(sam, "@SQ\tSN:c{i}\tLN:{CONTIG_LEN}").unwrap();
    }
    for (i, l) in sample.lines.iter().enumerate() {
        let rname = l.contig.map_or("*".to_owned(), |c| format!("c{c}"));
        writeln!(
            sam,
            "r{i}\t{}\t{rname}\t{}\t{}\t{}\t*\t0\t0\t{}\t{}",
            l.flag, l.pos, l.mapq, l.cigar, l.seq, l.qual
        )
        .unwrap();
    }
    let sam_path = dir.join("in.sam");
    std::fs::write(&sam_path, sam).unwrap();
    let bam = dir.join("sorted.bam");
    run(Command::new("samtools").args(["sort", "--no-PG", "-o"]).arg(&bam).arg(&sam_path));

    let cram = dir.join("out.cram");
    let mut cmd = Command::new("samtools");
    cmd.args(["view", "--no-PG", "-C", "-T"]).arg(&fasta);
    for opt in opts {
        cmd.arg("--output-fmt-option").arg(opt);
    }
    run(cmd.arg("-o").arg(&cram).arg(&bam));
    run(Command::new("samtools").arg("index").arg(&cram));
    (cram, fasta)
}

/// Every record of the CRAM, streamed.
fn stream(cram: &Path, fasta: &Path) -> (seqair::bam::BamHeader, Vec<RawBamRecord>) {
    let mut reader = CramStreamReader::new(std::fs::File::open(cram).unwrap()).unwrap();
    let mut decoder = reader.decoder(fasta).unwrap();
    let mut records = Vec::new();
    while let Some(container) = reader.next_container().unwrap() {
        decoder.decode(&container, &mut records).unwrap();
    }
    (reader.into_header(), records)
}

fn samtools_view(cram: &Path, fasta: &Path, args: &[&str]) -> Vec<u8> {
    run(Command::new("samtools")
        .args(["view", "--no-PG", "--reference"])
        .arg(fasta)
        .args(args)
        .arg(cram))
}

/// samtools' records as BAM payloads, through `samtools view -b`.
fn samtools_payloads(cram: &Path, fasta: &Path) -> Vec<Vec<u8>> {
    let bam = samtools_view(cram, fasta, &["-b"]);
    let mut reader = BamStreamReader::new(&bam[..]).unwrap();
    let mut out = Vec::new();
    let mut payload = Vec::new();
    while reader.read_payload_into(&mut payload).unwrap() {
        out.push(std::mem::take(&mut payload));
    }
    out
}

fn sam_line(record: &RawBamRecord, header: &seqair::bam::BamHeader) -> String {
    let mut out = Vec::new();
    format_sam_record(record, header, &mut out).unwrap();
    String::from_utf8(out).unwrap()
}

fn kind_of(sample: &Sample, qname: &[u8]) -> Kind {
    let i: usize = std::str::from_utf8(qname).unwrap().trim_start_matches('r').parse().unwrap();
    sample.lines[i].kind
}

// r[verify cram.edge.seq_unknown+3]
// r[verify cram.stream.fidelity+2]
/// A mapped record stored with SEQ and QUAL `*` comes back with neither,
/// from the stream and from the indexed reader, under any writer layout;
/// samtools prints it as the stream formats it.
#[hegel::test(test_cases = pinned::cases(40))]
fn mapped_reads_without_seq_come_back_without_seq(tc: TestCase) {
    let sample = tc.draw(arb_secondary_sample());
    let opts = tc.draw(arb_opts());
    let dir = tempfile::tempdir().unwrap();
    let (cram, fasta) = write_cram(dir.path(), &sample, &opts);

    let (header, records) = stream(&cram, &fasta);
    let theirs = String::from_utf8(samtools_view(&cram, &fasta, &[])).unwrap();
    let theirs: Vec<&str> = theirs.lines().collect();
    assert_eq!(records.len(), theirs.len(), "{opts:?}: record count");
    let mut seen = 0;
    for (record, want) in records.iter().zip(&theirs) {
        if kind_of(&sample, record.qname()) == Kind::SecondaryNoSeq {
            seen += 1;
            assert_eq!(record.seq_len(), 0, "{opts:?}: {want}");
            assert!(record.qual().is_none_or(<[_]>::is_empty), "{opts:?}: {want}");
            assert_eq!(&sam_line(record, &header), want, "{opts:?}");
        }
    }
    assert!(seen > 0, "the sample holds a secondary without SEQ");

    let mut indexed = IndexedCramReader::open(&cram, &fasta).unwrap();
    let mut store = RecordStore::new();
    let mut fetched = 0;
    for c in 0..sample.contigs.len() {
        let tid = indexed.header().tid(&format!("c{c}")).unwrap();
        indexed.fetch_into(tid, (Pos0::ZERO..=Pos0::MAX).into(), &mut store).unwrap();
        for i in store.indices() {
            let r = store.record(i).unwrap();
            if kind_of(&sample, r.qname()) == Kind::SecondaryNoSeq {
                fetched += 1;
                assert_eq!(r.seq_len, 0, "{opts:?}: indexed {:?}", r.qname());
                assert!(r.qual().is_empty(), "{opts:?}: indexed {:?}", r.qname());
            }
        }
    }
    assert_eq!(fetched, seen, "{opts:?}: the indexed reader returns every secondary");
}

// r[verify cram.stream.all_records]
// r[verify cram.stream.raw_records]
// r[verify cram.stream.fidelity+2]
// r[verify cram.edge.unmapped_reads+2]
/// Streaming keeps every record — unplaced and placed-unmapped ones
/// included, in whatever slices the writer put them — as many as
/// `samtools view -c` counts; an unplaced one has no reference, no position
/// and bin 4680; and every record is byte for byte what `samtools view -b`
/// writes and formats as `samtools view` prints it.
#[hegel::test(test_cases = pinned::cases(40))]
fn streamed_unplaced_reads_are_kept_as_samtools_reads_them(tc: TestCase) {
    let sample = tc.draw(arb_unplaced_sample());
    let opts = tc.draw(arb_opts());
    let dir = tempfile::tempdir().unwrap();
    let (cram, fasta) = write_cram(dir.path(), &sample, &opts);

    let (header, records) = stream(&cram, &fasta);
    let count = String::from_utf8(samtools_view(&cram, &fasta, &["-c"])).unwrap();
    assert_eq!(records.len(), count.trim().parse::<usize>().unwrap(), "{opts:?}: samtools -c");
    assert_eq!(records.len(), sample.lines.len(), "{opts:?}: generated");

    let mut unplaced = 0;
    for record in &records {
        if kind_of(&sample, record.qname()) == Kind::Unplaced {
            unplaced += 1;
            let what = format!("{opts:?}: {:?}", record.qname().escape_ascii().to_string());
            assert_eq!(record.tid(), None, "{what}");
            assert_eq!(record.pos(), None, "{what}");
            assert_eq!(record.bin(), UNPLACED_BIN, "{what}");
        }
    }
    assert!(unplaced > 0, "the sample holds an unplaced read");

    let payloads = samtools_payloads(&cram, &fasta);
    let theirs = String::from_utf8(samtools_view(&cram, &fasta, &[])).unwrap();
    let theirs: Vec<&str> = theirs.lines().collect();
    assert_eq!(payloads.len(), records.len(), "{opts:?}");
    for ((record, payload), line) in records.iter().zip(&payloads).zip(&theirs) {
        assert_eq!(&sam_line(record, &header), line, "{opts:?}");
        assert!(record.as_payload() == payload.as_slice(), "{opts:?}: payload of {line}");
    }
}
