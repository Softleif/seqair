//! The streaming CRAM decoder against samtools on generated files: for any
//! generated sample and any CRAM writer options, the streamed records are
//! the records `samtools view -b` writes and format as `samtools view` prints
//! them (`r[cram.stream.fidelity+2]`).
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "test code"
)]

#[path = "../src/pinned.rs"]
mod pinned;
use hegel::prelude::*;
use seqair::bam::RawBamRecord;
use seqair::bam::stream::BamStreamReader;
use seqair::cram::stream::CramStreamReader;
use seqair::sam::format::format_sam_record;
use std::fmt::Write as _;
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

const CONTIG_LEN: usize = 200;
const QUAL_CHARS: &str = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI";

#[derive(Debug, Clone)]
struct Sample {
    contigs: Vec<String>,
    /// SAM alignment lines.
    lines: Vec<String>,
    /// Whether the header's last `@RG` is `UNKNOWN`, which htslib drops
    /// from records (`r[cram.record.rg_tag+2]`).
    unknown_rg: bool,
}

fn text(alphabet: &str, len: usize) -> impl PrintableGenerator<String> {
    gs::text().alphabet(alphabet).min_size(len).max_size(len)
}

fn int(lo: usize, hi: usize) -> impl PrintableGenerator<usize> {
    gs::integers::<usize>().min_value(lo).max_value(hi)
}

/// Aux tags in an order samtools keeps: stored MD/NM that may or may not
/// equal what htslib recomputes, a read group in the header or not.
fn arb_tags(tc: &TestCase) -> String {
    let mut tags = Vec::new();
    if tc.draw_silent(gs::booleans()) {
        tags.push(format!("NM:i:{}", tc.draw_silent(int(0, 3))));
    }
    if tc.draw_silent(gs::booleans()) {
        tags.push(format!("MD:Z:{}", tc.draw_silent(int(0, 40))));
    }
    if tc.draw_silent(gs::booleans()) {
        tags.push(format!("AS:i:{}", tc.draw_silent(int(0, 300))));
    }
    if let Some(rg) =
        tc.draw_silent(gs::sampled_from(&[None, Some("grpA"), Some("grpB"), Some("x"), Some("UNKNOWN")]))
    {
        tags.push(format!("RG:Z:{rg}"));
    }
    if tc.draw_silent(gs::booleans()) {
        tags.reverse();
    }
    tags.iter().fold(String::new(), |mut out, t| {
        out.push('\t');
        out.push_str(t);
        out
    })
}

/// A mapped read: soft clips, matches, at most one insertion, deletion or
/// reference skip, the aligned bases copied from the reference with a few
/// substituted.
fn arb_mapped(tc: &TestCase, contig: usize, reference: &str) -> String {
    let mut ops = Vec::new();
    if tc.draw_silent(gs::booleans()) {
        ops.push((tc.draw_silent(int(1, 5)), 'S'));
    }
    ops.push((tc.draw_silent(int(5, 25)), 'M'));
    if let Some(op) = tc.draw_silent(gs::sampled_from(&[None, Some('I'), Some('D'), Some('N')])) {
        let len = if op == 'N' { tc.draw_silent(int(1, 60)) } else { tc.draw_silent(int(1, 4)) };
        ops.push((len, op));
        ops.push((tc.draw_silent(int(5, 25)), 'M'));
    }
    if tc.draw_silent(gs::booleans()) {
        ops.push((tc.draw_silent(int(1, 5)), 'S'));
    }
    let span: usize =
        ops.iter().filter(|(_, op)| matches!(op, 'M' | 'D' | 'N')).map(|(l, _)| l).sum();
    let pos = tc.draw_silent(int(1, CONTIG_LEN - span + 1));

    let mut seq = Vec::new();
    let mut at = pos - 1;
    for &(len, op) in &ops {
        match op {
            'M' => {
                seq.extend_from_slice(&reference.as_bytes()[at..at + len]);
                at += len;
            }
            'D' | 'N' => at += len,
            _ => seq.extend_from_slice(tc.draw_silent(text("ACGT", len)).as_bytes()),
        }
    }
    for _ in 0..tc.draw_silent(int(0, 3)) {
        let i = tc.draw_silent(int(0, seq.len() - 1));
        // IUPAC codes, `=` and lower case are stored verbatim as `B`/`b`.
        seq[i] = tc.draw_silent(gs::sampled_from(b"ACGTNRY=acgtn"));
    }
    let mut cigar = String::new();
    for (l, op) in &ops {
        write!(cigar, "{l}{op}").unwrap();
    }
    let flags = tc.draw_silent(gs::sampled_from(&[0u32, 16, 256]));
    // bwa writes secondary records with SEQ and QUAL `*`.
    let (seq, qual) = if flags == 256 && tc.draw_silent(gs::booleans()) {
        ("*".to_owned(), "*".to_owned())
    } else {
        let qual = if tc.draw_silent(gs::booleans()) {
            tc.draw_silent(text(QUAL_CHARS, seq.len()))
        } else {
            "*".to_owned()
        };
        (String::from_utf8(seq).unwrap(), qual)
    };
    let mapq = tc.draw_silent(int(0, 60));
    format!("{flags}\tc{contig}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{seq}\t{qual}{}", arb_tags(tc))
}

#[hegel::composite]
fn arb_sample(tc: &TestCase) -> Sample {
    // Reference `N`s, alone and in runs, are mismatches in a generated MD;
    // other IUPAC codes are matches, copied into SEQ as themselves.
    let contigs: Vec<String> = (0..tc.draw_silent(int(1, 3)))
        .map(|_| tc.draw_silent(text("ACGTACGTACGTNRYMacgt", CONTIG_LEN)))
        .collect();
    let mut lines = Vec::new();
    for (i, reference) in contigs.iter().enumerate() {
        for _ in 0..tc.draw_silent(int(0, 5)) {
            lines.push(format!("r{}\t{}", lines.len(), arb_mapped(tc, i, reference)));
        }
    }
    for _ in 0..tc.draw_silent(int(0, 3)) {
        let pair = arb_pair(tc, &contigs, &format!("p{}", lines.len()));
        lines.extend(pair);
    }
    for _ in 0..tc.draw_silent(int(0, 2)) {
        let len = tc.draw_silent(int(5, 20));
        let seq = tc.draw_silent(text("ACGTN", len));
        let tags = arb_tags(tc);
        lines.push(format!("u{}\t4\t*\t0\t0\t*\t*\t0\t0\t{seq}\t*{tags}", lines.len()));
    }
    Sample { contigs, lines, unknown_rg: tc.draw_silent(gs::booleans()) }
}

/// The reference span of a SAM CIGAR.
fn ref_span(cigar: &str) -> usize {
    let mut span = 0;
    let mut len = 0;
    for c in cigar.chars() {
        match c.to_digit(10) {
            Some(d) => len = len * 10 + d as usize,
            None => {
                if matches!(c, 'M' | 'D' | 'N') {
                    span += len;
                }
                len = 0;
            }
        }
    }
    span
}

/// A read pair: both mapped on one reference, on two, or one mate unmapped
/// at the other's position; TLEN as samtools computes it, which lets the
/// writer attach the mates, or 0. In a `lossy_names` file attached mates
/// lose their names (`r[cram.record.read_name+2]`).
fn arb_pair(tc: &TestCase, contigs: &[String], name: &str) -> [String; 2] {
    let c1 = tc.draw_silent(int(0, contigs.len() - 1));
    let first = arb_mapped(tc, c1, &contigs[c1]);
    let f1: Vec<&str> = first.split('\t').collect();
    let (pos1, span1) = (f1[2].parse::<usize>().unwrap(), ref_span(f1[4]).max(1));
    let rev1 = tc.draw_silent(gs::booleans());
    let rev2 = tc.draw_silent(gs::booleans());
    let proper = tc.draw_silent(gs::booleans());
    let (r1, r2) = if tc.draw_silent(gs::booleans()) { (0x40, 0x80) } else { (0x80, 0x40) };
    let flags = |own_rev: bool, mate_rev: bool, read: u32| {
        0x1 | read
            | if proper { 0x2 } else { 0 }
            | if own_rev { 0x10 } else { 0 }
            | if mate_rev { 0x20 } else { 0 }
    };
    let rest1 = f1[3..5].join("\t");
    let tail1 = f1[8..].join("\t");
    match tc.draw_silent(int(0, 2)) {
        // the mate unmapped, placed at the first's position
        0 => {
            let len = tc.draw_silent(int(5, 20));
            let seq = tc.draw_silent(text("ACGT", len));
            [
                format!(
                    "{name}\t{}\tc{c1}\t{pos1}\t{rest1}\t=\t{pos1}\t0\t{tail1}",
                    flags(rev1, rev2, r1) | 0x8
                ),
                format!(
                    "{name}\t{}\tc{c1}\t{pos1}\t0\t*\t=\t{pos1}\t0\t{seq}\t*",
                    flags(rev2, rev1, r2) | 0x4
                ),
            ]
        }
        other => {
            let c2 = if other == 1 { c1 } else { tc.draw_silent(int(0, contigs.len() - 1)) };
            let second = arb_mapped(tc, c2, &contigs[c2]);
            let f2: Vec<&str> = second.split('\t').collect();
            let (pos2, span2) = (f2[2].parse::<usize>().unwrap(), ref_span(f2[4]).max(1));
            let (rnext1, rnext2) = if c1 == c2 {
                ("=".to_owned(), "=".to_owned())
            } else {
                (format!("c{c2}"), format!("c{c1}"))
            };
            let tlen = if c1 == c2 && tc.draw_silent(gs::booleans()) {
                let left = pos1.min(pos2);
                let right = (pos1 + span1).max(pos2 + span2);
                let t = i64::try_from(right - left).unwrap();
                if pos1 <= pos2 { t } else { -t }
            } else {
                0
            };
            [
                format!(
                    "{name}\t{}\tc{c1}\t{pos1}\t{rest1}\t{rnext1}\t{pos2}\t{tlen}\t{tail1}",
                    flags(rev1, rev2, r1)
                ),
                format!(
                    "{name}\t{}\tc{c2}\t{pos2}\t{}\t{rnext2}\t{pos1}\t{}\t{}",
                    flags(rev2, rev1, r2),
                    f2[3..5].join("\t"),
                    -tlen,
                    f2[8..].join("\t")
                ),
            ]
        }
    }
}

#[derive(Debug, Clone, Copy)]
enum Reference {
    External,
    Embedded,
    Consensus,
    None,
}

/// How the input is ordered before it becomes CRAM.
#[derive(Debug, Clone, Copy)]
enum Order {
    Coordinate,
    /// `samtools sort -n`.
    Queryname,
    /// As generated, `SO:unsorted`.
    Unsorted,
}

#[derive(Debug, Clone)]
struct Opts {
    order: Order,
    version: &'static str,
    reference: Reference,
    seqs_per_slice: u32,
    slices_per_container: u32,
    multi_seq: Option<bool>,
    lossy_names: bool,
}

impl Opts {
    fn args(&self) -> Vec<String> {
        let mut args = vec![
            format!("version={}", self.version),
            format!("seqs_per_slice={}", self.seqs_per_slice),
            format!("slices_per_container={}", self.slices_per_container),
        ];
        match self.reference {
            Reference::External => {}
            Reference::Embedded => args.push("embed_ref=1".into()),
            Reference::Consensus => args.push("embed_ref=2".into()),
            Reference::None => args.push("no_ref=1".into()),
        }
        if self.lossy_names {
            args.push("lossy_names=1".into());
        }
        if let Some(m) = self.multi_seq {
            args.push(format!("multi_seq_per_slice={}", u8::from(m)));
        }
        args
    }
}

#[hegel::composite]
fn arb_opts(tc: &TestCase) -> Opts {
    let reference = tc.draw_silent(gs::sampled_from(&[
        Reference::External,
        Reference::Embedded,
        Reference::Consensus,
        Reference::None,
    ]));
    // samtools builds a consensus from sorted data only.
    let orders: &[Order] = if matches!(reference, Reference::Consensus) {
        &[Order::Coordinate]
    } else {
        &[Order::Coordinate, Order::Queryname, Order::Unsorted]
    };
    Opts {
        order: tc.draw_silent(gs::sampled_from(orders)),
        version: tc.draw_silent(gs::sampled_from(&["3.0", "3.1"])),
        reference,
        seqs_per_slice: tc.draw_silent(gs::sampled_from(&[1u32, 2, 5, 10_000])),
        slices_per_container: tc.draw_silent(gs::sampled_from(&[1u32, 2])),
        multi_seq: tc.draw_silent(gs::sampled_from(&[None, Some(false), Some(true)])),
        lossy_names: tc.draw_silent(gs::booleans()),
    }
}

fn run(cmd: &mut Command) -> Vec<u8> {
    let out = cmd.stderr(Stdio::piped()).output().expect("samtools on PATH");
    assert!(out.status.success(), "{cmd:?}: {}", String::from_utf8_lossy(&out.stderr));
    out.stdout
}

/// The CRAM, its reference, and how many records it holds.
fn write_files(dir: &Path, sample: &Sample, opts: &Opts) -> (PathBuf, PathBuf, usize) {
    let fasta = dir.join("ref.fa");
    let mut fa = String::new();
    let mut sam = String::from("@HD\tVN:1.6\tSO:unsorted\n");
    for (i, contig) in sample.contigs.iter().enumerate() {
        writeln!(fa, ">c{i}\n{contig}").unwrap();
        writeln!(sam, "@SQ\tSN:c{i}\tLN:{CONTIG_LEN}").unwrap();
    }
    sam.push_str("@RG\tID:grpA\tSM:x\n@RG\tID:grpB\tSM:x\n");
    if sample.unknown_rg {
        sam.push_str("@RG\tID:UNKNOWN\tSM:x\n");
    }
    let mut written = 0;
    for line in &sample.lines {
        // samtools 1.25 fails to build a consensus for a slice of mapped
        // SEQ `*` reads ("Failed to build reference") and writes a broken file.
        let fields: Vec<&str> = line.split('\t').collect();
        let mapped_seq_star = fields[1].parse::<u32>().unwrap() & 4 == 0 && fields[9] == "*";
        if matches!(opts.reference, Reference::Consensus) && mapped_seq_star {
            continue;
        }
        writeln!(sam, "{line}").unwrap();
        written += 1;
    }
    std::fs::write(&fasta, fa).unwrap();
    run(Command::new("samtools").arg("faidx").arg(&fasta));
    let sam_path = dir.join("in.sam");
    std::fs::File::create(&sam_path).unwrap().write_all(sam.as_bytes()).unwrap();
    let bam = dir.join("in.bam");
    let mut sort = Command::new("samtools");
    match opts.order {
        Order::Coordinate => sort.args(["sort", "--no-PG"]),
        Order::Queryname => sort.args(["sort", "--no-PG", "-n"]),
        Order::Unsorted => sort.args(["view", "--no-PG", "-b"]),
    };
    run(sort.arg("-o").arg(&bam).arg(&sam_path));
    let cram = dir.join("generated.cram");
    let mut cmd = Command::new("samtools");
    cmd.args(["view", "--no-PG", "-C"]);
    for opt in opts.args() {
        cmd.arg("--output-fmt-option").arg(opt);
    }
    run(cmd.arg("-T").arg(&fasta).arg("-o").arg(&cram).arg(&bam));
    (cram, fasta, written)
}

fn streamed(cram: &Path, fasta: &Path, decode_md: bool) -> Vec<RawBamRecord> {
    let mut reader = CramStreamReader::new(std::fs::File::open(cram).unwrap()).unwrap();
    let prefix = cram.file_name().unwrap().to_str().unwrap();
    let mut decoder = reader.decoder(fasta).unwrap().decode_md(decode_md).name_prefix(prefix);
    let mut records = Vec::new();
    while let Some(container) = reader.next_container().unwrap() {
        decoder.decode(&container, &mut records).unwrap();
    }
    records
}

// r[verify cram.stream.fidelity+2]
// r[verify cram.stream.md_nm]
// r[verify cram.stream.md_nm.option]
// r[verify cram.stream.iupac]
// r[verify cram.stream+2]
// r[verify cram.record.read_name+2]
// r[verify cram.record.mate_tlen_reconstruction+2]
// r[verify cram.record.rg_tag+2]
/// For any generated sample and CRAM writer options, the streamed records
/// equal `samtools view -b`'s, byte for byte, and format as `samtools view`
/// prints them, with MD/NM generation on (samtools' default) or off.
#[hegel::test(test_cases = pinned::cases(40))]
fn streamed_records_are_what_samtools_writes(tc: TestCase) {
    let sample = tc.draw(arb_sample().print_as_debug());
    let opts = tc.draw(arb_opts().print_as_debug());
    let decode_md = tc.draw(gs::booleans());
    assert_streams_as_samtools(&sample, &opts, decode_md);
}

/// An unplaced read with SEQ and QUAL `*` has no bases, so a slice of only
/// such reads has no base or quality block at all.
const UNPLACED_NO_SEQ: &str = "4\t*\t0\t0\t*\t*\t0\t0\t*\t*";

#[hegel::composite]
fn arb_unplaced_no_seq_sample(tc: &TestCase) -> Sample {
    let reference = tc.draw_silent(text("ACGT", CONTIG_LEN));
    let mut lines: Vec<String> = (0..tc.draw_silent(int(0, 4)))
        .map(|i| format!("r{i}\t{}", arb_mapped(tc, 0, &reference)))
        .collect();
    for _ in 0..tc.draw_silent(int(1, 8)) {
        lines.push(format!("u{}\t{UNPLACED_NO_SEQ}", lines.len()));
    }
    Sample { contigs: vec![reference], lines, unknown_rg: false }
}

// r[verify cram.stream.all_records]
// r[verify cram.encoding.external_empty]
/// Unplaced reads with SEQ and QUAL `*`, alone in their slice or after
/// mapped reads, stream as samtools reads them.
#[hegel::test(test_cases = pinned::cases(24))]
fn unplaced_reads_without_seq_stream_as_samtools_reads_them(tc: TestCase) {
    let sample = tc.draw(arb_unplaced_no_seq_sample().print_as_debug());
    let opts = tc.draw(arb_opts().print_as_debug());
    assert_streams_as_samtools(&sample, &opts, true);
}

/// The streamed records of `sample` written under `opts` against samtools';
/// one per record written.
fn assert_streams_as_samtools(sample: &Sample, opts: &Opts, decode_md: bool) {
    let dir = tempfile::tempdir().unwrap();
    let (cram, fasta, written) = write_files(dir.path(), sample, opts);

    let ours = streamed(&cram, &fasta, decode_md);
    let header = CramStreamReader::new(std::fs::File::open(&cram).unwrap()).unwrap().into_header();
    let md_args: &[&str] = if decode_md { &[] } else { &["--input-fmt-option", "decode_md=0"] };
    let view = |extra: &[&str]| {
        run(Command::new("samtools")
            .args(["view", "--no-PG"])
            .args(extra)
            .args(md_args)
            .arg("--reference")
            .arg(&fasta)
            .arg(&cram))
    };

    let lines = String::from_utf8(view(&[])).unwrap();
    let theirs: Vec<&str> = lines.lines().collect();
    let ours_sam: Vec<String> = ours
        .iter()
        .map(|r| {
            let mut out = Vec::new();
            format_sam_record(r, &header, &mut out).unwrap();
            String::from_utf8(out).unwrap()
        })
        .collect();
    assert_eq!(ours_sam, theirs);

    let bam = view(&["-b"]);
    let mut reader = BamStreamReader::new(&bam[..]).unwrap();
    for (i, record) in ours.iter().enumerate() {
        let mut payload = Vec::new();
        assert!(reader.read_payload_into(&mut payload).unwrap());
        assert_eq!(record.as_payload(), payload.as_slice(), "record {i}: {}", theirs[i]);
    }
    let mut rest = Vec::new();
    assert!(!reader.read_payload_into(&mut rest).unwrap());
    assert_eq!(ours.len(), written);
}
