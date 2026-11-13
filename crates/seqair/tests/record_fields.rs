//! The field opt-in API (`RecordFields`) on every format and entry point.
//!
//! The same generated reads become SAM (plain and tabix-indexed bgzip), BAM
//! and CRAM. Each reader's own `RecordFields::all()` decode is the oracle for
//! what a narrower request must return (that full decode is checked against
//! samtools elsewhere); what was not requested must be refused identically
//! everywhere. CRAM adds MD/NM where the file lacks them, so formats are only
//! compared to each other with those two tags set aside.
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

#[path = "../src/pinned.rs"]
mod pinned;
use hegel::prelude::*;
use seqair::bam::aux::{AuxValue, GetAuxError};
use seqair::bam::{Pos0, RecordStore};
use seqair::reader::{
    AlignmentStream, AuxTags, Field, IndexedReader, PartialRecord, RecordFields, Skipped,
};
use std::fmt::Write as _;
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

const CONTIG: &str = "ref1";
const CONTIG_LEN: u32 = 400;
/// Tags the generator may put on a read, plus one nobody has.
const TAGS: [[u8; 2]; 6] = [*b"NM", *b"MC", *b"RG", *b"XA", *b"XF", *b"XC"];
const ABSENT: [u8; 2] = *b"ZZ";

fn run(cmd: &mut Command, what: &str) {
    let status = cmd
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .unwrap_or_else(|e| panic!("{what}: {e}"));
    assert!(status.success(), "{what} failed");
}

// ── generated reads ─────────────────────────────────────────────────────

#[derive(Debug, Clone)]
struct Read {
    /// 1-based; 0 for an unplaced read.
    pos: u32,
    flags: u32,
    cigar: String,
    seq: String,
    /// `*` or one character per base.
    qual: String,
    /// SAM text of each optional field, in the order written.
    tags: Vec<String>,
}

fn bases(len: usize) -> impl PrintableGenerator<String> {
    gs::text().alphabet("ACGT").min_size(len).max_size(len)
}

#[hegel::composite]
fn arb_read(tc: &TestCase) -> Read {
    let n = |lo: u32, hi: u32| gs::integers::<u32>().min_value(lo).max_value(hi);
    // Mapped reads dominate; the other two kinds are the awkward ones.
    let kind = tc.draw_silent(gs::sampled_from(&[0u8, 0, 0, 1, 2]));
    let len = tc.draw_silent(n(8, 40)) as usize;
    let seq = tc.draw_silent(bases(len));
    let qual = if tc.draw_silent(gs::sampled_from(&[true, true, false])) {
        tc.draw_silent(gs::text().alphabet("5?I#+").min_size(len).max_size(len))
    } else {
        "*".to_owned()
    };
    let (pos, flags, cigar) = match kind {
        0 => (
            tc.draw_silent(n(1, CONTIG_LEN - 80)),
            tc.draw_silent(gs::sampled_from(&[0u32, 16])),
            format!("{len}M"),
        ),
        1 => (tc.draw_silent(n(1, CONTIG_LEN - 80)), 4, "*".to_owned()),
        _ => (0, 4, "*".to_owned()),
    };

    let mut tags = Vec::new();
    for tag in TAGS {
        if !tc.draw_silent(gs::booleans()) {
            continue;
        }
        let name = std::str::from_utf8(&tag).expect("ASCII");
        tags.push(match &tag {
            b"NM" => format!("NM:i:{}", tc.draw_silent(n(0, 70_000))),
            b"MC" => format!("MC:Z:{}M", tc.draw_silent(n(1, 99))),
            b"RG" => format!("RG:Z:g{}", tc.draw_silent(n(1, 3))),
            b"XA" => format!("XA:B:C,{},{}", tc.draw_silent(n(0, 255)), tc.draw_silent(n(0, 255))),
            b"XF" => format!("XF:f:{}.25", tc.draw_silent(n(0, 9))),
            _ => format!("{name}:A:{}", tc.draw_silent(gs::sampled_from(&['a', 'Q']))),
        });
    }
    Read { pos, flags, cigar, seq, qual, tags }
}

/// A reference plus reads aligned to it; mapped reads copy the reference so
/// CRAM stores them as matches.
#[hegel::composite]
fn arb_input(tc: &TestCase) -> (String, Vec<Read>) {
    let reference = tc.draw_silent(bases(CONTIG_LEN as usize));
    let n_reads = tc.draw_silent(gs::integers::<usize>().min_value(1).max_value(8));
    let mut reads: Vec<Read> = (0..n_reads).map(|_| tc.draw_silent(arb_read())).collect();
    for r in &mut reads {
        if r.flags & 4 == 0 {
            let at = r.pos as usize - 1;
            r.seq = reference[at..at + r.seq.len()].to_owned();
        }
    }
    (reference, reads)
}

fn sam_text(reads: &[Read]) -> String {
    let mut sam = format!(
        "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:{CONTIG}\tLN:{CONTIG_LEN}\n\
         @RG\tID:g1\n@RG\tID:g2\n@RG\tID:g3\n"
    );
    for (i, r) in reads.iter().enumerate() {
        let rname = if r.pos == 0 { "*" } else { CONTIG };
        write!(
            sam,
            "r{i}\t{}\t{rname}\t{}\t{}\t{}\t*\t0\t0\t{}\t{}",
            r.flags,
            r.pos,
            if r.flags & 4 == 0 { 30 } else { 0 },
            r.cigar,
            r.seq,
            r.qual
        )
        .expect("String write");
        for t in &r.tags {
            write!(sam, "\t{t}").expect("String write");
        }
        sam.push('\n');
    }
    sam
}

struct Files {
    _dir: tempfile::TempDir,
    fasta: PathBuf,
    /// Plain SAM (stream only), tabix-indexed bgzip SAM, BAM, CRAM.
    sam: PathBuf,
    sam_gz: PathBuf,
    bam: PathBuf,
    cram: PathBuf,
}

fn write_files(reference: &str, reads: &[Read]) -> Files {
    let dir = tempfile::tempdir().expect("tempdir");
    let p = |name: &str| dir.path().join(name);

    let fasta = p("ref.fa");
    let mut f = std::fs::File::create(&fasta).expect("create FASTA");
    writeln!(f, ">{CONTIG}\n{reference}").expect("write");
    drop(f);
    run(Command::new("samtools").arg("faidx").arg(&fasta), "samtools faidx");

    let raw = p("raw.sam");
    std::fs::write(&raw, sam_text(reads)).expect("write SAM");
    let bam = p("g.bam");
    run(Command::new("samtools").args(["sort", "-o"]).arg(&bam).arg(&raw), "sort");
    run(Command::new("samtools").arg("index").arg(&bam), "index bam");

    let cram = p("g.cram");
    run(
        Command::new("samtools")
            .args(["view", "-C", "-T"])
            .arg(&fasta)
            .arg("-o")
            .arg(&cram)
            .arg(&bam),
        "view -C",
    );
    run(Command::new("samtools").arg("index").arg(&cram), "index cram");

    let sam = p("g.sam");
    run(Command::new("samtools").args(["view", "-h", "-o"]).arg(&sam).arg(&bam), "view sam");
    let sam_gz = p("g.sam.gz");
    run(
        Command::new("samtools")
            .args(["view", "-h", "--output-fmt", "SAM,level=6", "-o"])
            .arg(&sam_gz)
            .arg(&bam),
        "view sam.gz",
    );
    run(Command::new("tabix").args(["-p", "sam"]).arg(&sam_gz), "tabix");

    Files { _dir: dir, fasta, sam, sam_gz, bam, cram }
}

// ── what a record shows ─────────────────────────────────────────────────

/// An aux value with integer widths erased: formats may pick another width.
fn norm(v: &AuxValue<'_>) -> String {
    match *v {
        AuxValue::I8(x) => format!("int {x}"),
        AuxValue::U8(x) => format!("int {x}"),
        AuxValue::I16(x) => format!("int {x}"),
        AuxValue::U16(x) => format!("int {x}"),
        AuxValue::I32(x) => format!("int {x}"),
        AuxValue::U32(x) => format!("int {x}"),
        ref other => format!("{other:?}"),
    }
}

type Tags = Vec<([u8; 2], String)>;

fn sorted(mut t: Tags) -> Tags {
    t.sort();
    t
}

/// Everything a `PartialRecord` exposes, each gated field as a `Result`.
#[derive(Debug, Clone, PartialEq)]
struct Snap {
    core: String,
    seq_len: usize,
    qname: Result<Vec<u8>, Skipped>,
    seq: Result<Vec<u8>, Skipped>,
    qual: Result<Option<Vec<u8>>, Skipped>,
    aux: Tags,
    /// Per probed tag: refused as `NotRequested`?
    not_requested: Vec<([u8; 2], bool)>,
    complete: bool,
}

fn snap(rec: PartialRecord) -> Snap {
    let core = format!(
        "{:?} {:?} {:?} {:?} {} {:?} {:?}",
        rec.tid(),
        rec.pos(),
        rec.flags(),
        rec.mapq(),
        rec.cigar().map(|c| c.to_string()).collect::<String>(),
        rec.mate_fields(),
        rec.end_pos(),
    );
    let aux = sorted(rec.aux().iter_tags().map(|(t, v)| (t, norm(&v))).collect());
    let not_requested = TAGS
        .iter()
        .chain([&ABSENT])
        .map(|&t| {
            let refused = matches!(rec.aux().get::<i64>(t), Err(GetAuxError::NotRequested { .. }));
            (t, refused)
        })
        .collect();
    let qname = rec.qname().map(<[u8]>::to_vec);
    let seq = rec.seq().map(|s| s.iupac_bases().map(|b| b.to_ascii()).collect());
    let qual = rec.qual().map(|q| q.map(|q| seqair_types::BaseQuality::slice_to_bytes(q).to_vec()));
    let seq_len = rec.seq_len();
    let fields_were_all = rec.fields().is_all();
    let complete = match rec.into_complete() {
        Ok(_) => true,
        Err(back) => {
            assert!(!back.fields().is_all());
            false
        }
    };
    assert_eq!(complete, fields_were_all);
    Snap { core, seq_len, qname, seq, qual, aux, not_requested, complete }
}

/// The same in the shape `fetch_into` stores: no gated fields.
#[derive(Debug, Clone, PartialEq)]
struct Stored {
    core: String,
    qname: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>,
    aux: Tags,
}

#[derive(Clone, Copy, Debug)]
enum Fmt {
    Sam,
    SamGz,
    Bam,
    Cram,
}

fn path(files: &Files, f: Fmt) -> &Path {
    match f {
        Fmt::Sam => &files.sam,
        Fmt::SamGz => &files.sam_gz,
        Fmt::Bam => &files.bam,
        Fmt::Cram => &files.cram,
    }
}

fn via_stream(files: &Files, f: Fmt, fields: &RecordFields) -> Vec<Snap> {
    let mut s = AlignmentStream::open_path(path(files, f))
        .expect("open")
        .reference(&files.fasta)
        .fields(fields.clone())
        .build()
        .expect("build");
    s.records().map(|r| snap(r.expect("record"))).collect()
}

fn span() -> std::range::RangeInclusive<Pos0> {
    (Pos0::new(0).unwrap()..=Pos0::new(CONTIG_LEN - 1).unwrap()).into()
}

fn open(files: &Files, f: Fmt, fields: &RecordFields) -> seqair::reader::PartialIndexedReader {
    IndexedReader::builder(path(files, f))
        .reference(&files.fasta)
        .fields(fields.clone())
        .open()
        .expect("open indexed")
}

fn via_query(files: &Files, f: Fmt, fields: &RecordFields) -> Vec<Snap> {
    let mut r = open(files, f, fields);
    let tid = r.header().tid(CONTIG).expect("contig");
    r.query_raw(tid, span()).expect("query").map(|x| snap(x.expect("record"))).collect()
}

fn via_fetch(files: &Files, f: Fmt, fields: &RecordFields) -> Vec<Stored> {
    let mut r = open(files, f, fields);
    let tid = r.header().tid(CONTIG).expect("contig");
    let mut store = RecordStore::new();
    r.fetch_into(tid, span(), &mut store).expect("fetch");
    store
        .indices()
        .map(|i| {
            let rec = store.record(i).unwrap();
            Stored {
                core: format!("{:?} {:?} {}", rec.pos, rec.flags, rec.cigar().len()),
                qname: rec.qname().to_vec(),
                seq: rec.seq().iter().map(|b| b.to_ascii()).collect(),
                qual: seqair_types::BaseQuality::slice_to_bytes(rec.qual()).to_vec(),
                aux: sorted(
                    seqair::bam::aux::iter_tags(rec.aux()).map(|(t, v)| (t, norm(&v))).collect(),
                ),
            }
        })
        .collect()
}

// ── the checks ──────────────────────────────────────────────────────────

fn keep(tags: Tags, wanted: &AuxTags) -> Tags {
    tags.into_iter().filter(|(t, _)| wanted.contains(*t)).collect()
}

/// `got` (read with `fields`) against `full` (the same reader, all fields).
fn check_partial(got: &[Snap], full: &[Snap], fields: &RecordFields, what: &str) {
    assert_eq!(got.len(), full.len(), "{what}: record count");
    for (i, (g, f)) in got.iter().zip(full).enumerate() {
        let at = format!("{what} #{i}");
        assert!(f.complete, "{at}: all() must give a complete record");
        assert_eq!(g.complete, fields.is_all(), "{at}: into_complete");
        assert_eq!(g.core, f.core, "{at}: core");
        assert_eq!(g.seq_len, f.seq_len, "{at}: seq_len");
        assert_eq!(
            g.qname,
            if fields.names() { f.qname.clone() } else { Err(Skipped { field: Field::Name }) },
            "{at}: qname"
        );
        assert_eq!(
            g.seq,
            if fields.sequence() { f.seq.clone() } else { Err(Skipped { field: Field::Seq }) },
            "{at}: seq"
        );
        assert_eq!(
            g.qual,
            if fields.quality() { f.qual.clone() } else { Err(Skipped { field: Field::Qual }) },
            "{at}: qual"
        );
        assert_eq!(g.aux, keep(f.aux.clone(), fields.aux()), "{at}: aux tags");
        for (tag, refused) in &g.not_requested {
            assert_eq!(*refused, !fields.aux().contains(*tag), "{at}: NotRequested for {tag:?}");
        }
    }
}

/// What formats must agree on: MD/NM are set aside (CRAM may add them).
fn cross(s: &Snap) -> Snap {
    let mut s = s.clone();
    s.aux.retain(|(t, _)| t != b"MD" && t != b"NM");
    s.not_requested.retain(|(t, _)| t != b"NM");
    s
}

fn fields_from(names: bool, seq: bool, qual: bool, aux: AuxTags) -> RecordFields {
    let mut f = RecordFields::core().with_aux(aux);
    if names {
        f = f.with_names();
    }
    if seq {
        f = f.with_sequence();
    }
    if qual {
        f = f.with_quality();
    }
    f
}

#[hegel::composite]
fn arb_fields(tc: &TestCase) -> RecordFields {
    let aux = match tc.draw_silent(gs::integers::<u8>().max_value(2)) {
        0 => AuxTags::All,
        1 => AuxTags::None,
        _ => {
            let mut tags = Vec::new();
            for t in TAGS.iter().chain([&ABSENT]) {
                if tc.draw_silent(gs::booleans()) {
                    tags.push(*t);
                }
            }
            AuxTags::only(tags)
        }
    };
    let (n, s, q) = (
        tc.draw_silent(gs::booleans()),
        tc.draw_silent(gs::booleans()),
        tc.draw_silent(gs::booleans()),
    );
    // `all()` is also reachable by drawing every bit and AuxTags::All.
    fields_from(n, s, q, aux)
}

// r[verify unified.record_fields]
// r[verify unified.partial_record]
// r[verify unified.partial_record.skipped]
// r[verify unified.partial_record.l_seq]
// r[verify unified.stream.fields]
// r[verify unified.record_fields.indexed]
// r[verify unified.record_fields.fetch_into]
/// Whatever is requested, on every format and entry point: requested fields
/// equal the full decode, the rest are refused alike, and formats agree.
#[hegel::test(test_cases = crate::pinned::cases(24))]
fn requested_fields_match_the_full_decode_and_the_rest_are_refused(tc: TestCase) {
    let (reference, reads) = tc.draw(arb_input().print_as_debug());
    let fields = tc.draw(arb_fields().print_as_debug());
    let files = write_files(&reference, &reads);
    let all = RecordFields::all();
    let formats = [Fmt::Sam, Fmt::SamGz, Fmt::Bam, Fmt::Cram];

    // Stream: every format, including the unplaced reads.
    let mut streams = Vec::new();
    for f in formats {
        let full = via_stream(&files, f, &all);
        let got = via_stream(&files, f, &fields);
        check_partial(&got, &full, &fields, &format!("stream {f:?}"));
        streams.push(got);
    }
    cross_check(&streams, "stream");

    // Indexed query_raw and fetch_into: formats with an index.
    let indexed = [Fmt::SamGz, Fmt::Bam, Fmt::Cram];
    let mut queries = Vec::new();
    let mut stored = Vec::new();
    for f in indexed {
        let full = via_query(&files, f, &all);
        let got = via_query(&files, f, &fields);
        check_partial(&got, &full, &fields, &format!("query_raw {f:?}"));
        queries.push(got);

        let full = via_fetch(&files, f, &all);
        let got = via_fetch(&files, f, &fields);
        assert_eq!(got.len(), full.len(), "fetch_into {f:?}: count");
        for (i, (g, w)) in got.iter().zip(&full).enumerate() {
            let want = Stored { aux: keep(w.aux.clone(), fields.aux()), ..w.clone() };
            assert_eq!(g, &want, "fetch_into {f:?} #{i}");
        }
        stored.push(got);
    }
    cross_check(&queries, "query_raw");
    for s in &stored[1..] {
        assert_eq!(s.len(), stored[0].len(), "fetch_into count across formats");
        for (a, b) in stored[0].iter().zip(s) {
            let strip = |x: &Stored| Stored {
                aux: x.aux.iter().filter(|(t, _)| t != b"MD" && t != b"NM").cloned().collect(),
                ..x.clone()
            };
            assert_eq!(strip(a).aux, strip(b).aux, "fetch_into aux across formats");
            assert_eq!((&a.qname, &a.seq, &a.qual), (&b.qname, &b.seq, &b.qual));
        }
    }

    if fields.is_all() {
        tc.event("all fields");
    }
    if matches!(fields.aux(), AuxTags::Only(_)) {
        tc.event("aux subset");
    }
}

fn cross_check(per_format: &[Vec<Snap>], what: &str) {
    let first = &per_format[0];
    for other in &per_format[1..] {
        assert_eq!(other.len(), first.len(), "{what}: count across formats");
        for (i, (a, b)) in first.iter().zip(other).enumerate() {
            let (a, b) = (cross(a), cross(b));
            // Formats differ in how they spell an absent QUAL; only the
            // presence of the request is compared there.
            assert_eq!(
                (a.qname, a.seq, a.aux, a.not_requested, a.seq_len, a.complete),
                (b.qname, b.seq, b.aux, b.not_requested, b.seq_len, b.complete),
                "{what} #{i} across formats"
            );
            assert_eq!(a.qual.is_err(), b.qual.is_err(), "{what} #{i}: qual refusal");
        }
    }
}

/// A fixed request on a fixed file, so a failure here is readable.
#[test]
fn core_refuses_everything_else_on_every_format() {
    let reads = vec![Read {
        pos: 5,
        flags: 0,
        cigar: "10M".into(),
        seq: "ACGTACGTAC".into(),
        qual: "IIIIIIIIII".into(),
        tags: vec!["NM:i:1".into(), "MC:Z:10M".into()],
    }];
    let files = write_files(&"ACGT".repeat(100), &reads);
    for f in [Fmt::Sam, Fmt::SamGz, Fmt::Bam, Fmt::Cram] {
        let got = via_stream(&files, f, &RecordFields::core());
        assert_eq!(got.len(), 1);
        let g = &got[0];
        assert_eq!(g.qname, Err(Skipped { field: Field::Name }));
        assert_eq!(g.seq, Err(Skipped { field: Field::Seq }));
        assert_eq!(g.qual, Err(Skipped { field: Field::Qual }));
        assert!(g.aux.is_empty());
        assert_eq!(g.seq_len, 10);
        assert!(!g.complete);
    }
}
