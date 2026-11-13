//! CRAM decodes that skip fields (`r[cram.stream.fields]`): whatever a decode
//! leaves out, every field it was asked for is what the full decode gives, on
//! the stream, `query_raw` and `fetch_into`; a decode that needs no reference
//! runs without one. The oracle is the same reader's `RecordFields::all()`
//! decode, which other tests hold against samtools.
//!
//! samtools writes each series and tag to a block of its own; htsjdk shares
//! blocks between series, which the vendored htsjdk files cover.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "test code"
)]
#![allow(clippy::cast_possible_truncation, reason = "test code with known small values")]

#[path = "../src/pinned.rs"]
mod pinned;
use hegel::prelude::*;
use seqair::bam::{Pos0, RawBamRecord, RecordStore};
use seqair::cram::stream::CramStreamReader;
use seqair::reader::{
    AlignmentStream, AlignmentStreamError, AuxTags, FormatDetectionError, IndexedReader,
    PartialRecord, ReaderError, RecordFields, Tid,
};
use std::fmt::Write as _;
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn data(name: &str) -> PathBuf {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data")).join(name)
}

// ── comparing a partial record with the full one ────────────────────────

type Tags = Vec<([u8; 2], String)>;

fn tags<'a>(iter: impl Iterator<Item = ([u8; 2], seqair::bam::aux::AuxValue<'a>)>) -> Tags {
    iter.map(|(t, v)| (t, format!("{v:?}"))).collect()
}

fn requested(full: Tags, fields: &RecordFields) -> Tags {
    full.into_iter().filter(|(t, _)| fields.aux().contains(*t)).collect()
}

/// `got`, read with `fields`, against `full`, the same record read whole.
fn check(got: &PartialRecord, full: &RawBamRecord, fields: &RecordFields, at: &str) {
    assert_eq!(got.tid(), full.tid(), "{at}: tid");
    assert_eq!(got.pos(), full.pos(), "{at}: pos");
    assert_eq!(got.flags(), full.flags(), "{at}: flags");
    assert_eq!(got.mapq(), full.mapq(), "{at}: mapq");
    assert_eq!(got.mate_fields(), full.mate_fields(), "{at}: mate");
    assert_eq!(got.cigar().collect::<Vec<_>>(), full.cigar().collect::<Vec<_>>(), "{at}: cigar");
    assert_eq!(got.end_pos(), full.end_pos(), "{at}: end_pos");
    assert_eq!(got.seq_len(), full.seq_len(), "{at}: seq_len");
    assert_eq!(got.qname().ok(), fields.names().then(|| full.qname()), "{at}: qname");
    assert_eq!(
        got.seq().ok().map(|s| s.packed()),
        fields.sequence().then(|| full.packed_seq()),
        "{at}: seq"
    );
    assert_eq!(got.qual().ok(), fields.quality().then(|| full.qual()), "{at}: qual");
    assert_eq!(
        tags(got.aux().iter_tags()),
        requested(tags(full.aux().iter_tags()), fields),
        "{at}: tags"
    );
}

/// Every record of the stream of `cram`, read whole or with `fields`.
fn stream(
    cram: &Path,
    fasta: Option<&Path>,
    decode_md: bool,
    fields: Option<&RecordFields>,
) -> Result<Vec<PartialRecord>, seqair::cram::stream::CramStreamError> {
    let mut reader = CramStreamReader::new(
        std::fs::File::open(cram)
            .map_err(|source| seqair::cram::stream::CramStreamError::Io { source })?,
    )?;
    let decoder = match fasta {
        Some(fasta) => reader.decoder(fasta)?,
        None => reader.decoder_without_reference(),
    }
    .decode_md(decode_md);
    let mut decoder = decoder.fields(fields.cloned().unwrap_or_default());
    let mut out = Vec::new();
    while let Some(container) = reader.next_container()? {
        decoder.decode(&container, &mut out)?;
    }
    Ok(out)
}

fn whole(records: Vec<PartialRecord>) -> Vec<RawBamRecord> {
    records.into_iter().map(|r| r.into_complete().unwrap()).collect()
}

fn check_all(got: &[PartialRecord], full: &[RawBamRecord], fields: &RecordFields, what: &str) {
    assert_eq!(got.len(), full.len(), "{what}: records");
    for (i, (g, f)) in got.iter().zip(full).enumerate() {
        check(g, f, fields, &format!("{what} #{i}"));
    }
}

/// Whether a reader of `fields` needs the reference: for the sequence, or
/// to generate MD/NM (the stream and `query_raw` do by default).
fn needs_reference(fields: &RecordFields) -> bool {
    fields.sequence() || fields.aux().contains(*b"MD") || fields.aux().contains(*b"NM")
}

const PROBE_TAGS: [[u8; 2]; 9] =
    [*b"NM", *b"MD", *b"RG", *b"XA", *b"XF", *b"XB", *b"MC", *b"AS", *b"ZZ"];

#[hegel::composite]
fn arb_fields(tc: &TestCase) -> RecordFields {
    let mut f = RecordFields::core();
    if tc.draw(gs::booleans()) {
        f = f.with_names();
    }
    if tc.draw(gs::booleans()) {
        f = f.with_sequence();
    }
    if tc.draw(gs::booleans()) {
        f = f.with_quality();
    }
    let aux = match tc.draw(gs::integers::<u8>().max_value(2)) {
        0 => AuxTags::All,
        1 => AuxTags::None,
        _ => {
            let picks = tc.draw(gs::vecs(gs::sampled_from(&PROBE_TAGS)).max_size(4));
            AuxTags::only(picks)
        }
    };
    f.with_aux(aux)
}

// ── htsjdk's shared blocks ──────────────────────────────────────────────

// r[verify cram.stream.fields]
// r[verify cram.stream.fields.dropped_blocks]
// r[verify cram.stream.fields.seq_off]
// r[verify cram.stream.fields.qual_off]
// r[verify cram.stream.fields.md_nm]
/// htsjdk keeps `BA`, `IN` and `SC` in one block, `NP` and `TS` in another,
/// and the tags of one type in one: a decode that skips some of them must
/// still read the rest where the full decode does.
#[hegel::test]
fn shared_blocks_keep_the_kept_series_in_step(tc: TestCase) {
    let fields = tc.draw(arb_fields());
    let decode_md = tc.draw(gs::booleans());
    for (cram, fasta) in [
        (data("cram-stream/ce_5b_java.cram"), data("hts-specs/cram/ce.fa")),
        (data("cram-stream/auxf_values_java.cram"), data("cram-stream/auxf.fa")),
    ] {
        let full = whole(stream(&cram, Some(&fasta), decode_md, None).unwrap());
        assert!(!full.is_empty());
        let got = stream(&cram, Some(&fasta), decode_md, Some(&fields)).unwrap();
        check_all(&got, &full, &fields, &cram.display().to_string());
        if !fields.sequence() && !(decode_md && needs_reference(&fields)) {
            let got = stream(&cram, None, decode_md, Some(&fields)).unwrap();
            check_all(&got, &full, &fields, "without a reference");
        }
    }
}

// ── generated files ─────────────────────────────────────────────────────

const CONTIGS: [&str; 2] = ["c1", "c2"];
const CONTIG_LEN: u32 = 400;

fn run(cmd: &mut Command, what: &str) {
    let status = cmd
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .unwrap_or_else(|e| panic!("{what}: {e}"));
    assert!(status.success(), "{what} failed");
}

fn bases(len: usize) -> impl PrintableGenerator<String> {
    gs::text().alphabet("ACGT").min_size(len).max_size(len)
}

/// One SAM line's fields after the name.
#[derive(Debug, Clone)]
struct Read {
    name: String,
    flags: u32,
    /// Index into `CONTIGS`; `None` for an unplaced read.
    contig: Option<usize>,
    pos: u32,
    cigar: String,
    seq: String,
    qual: String,
    mate: Option<(usize, u32)>,
    tags: Vec<String>,
}

/// A mapped read over `reference` from 1-based `pos`: its CIGAR, and the
/// reference's bases with a few changed, some to `N` or an IUPAC code
/// (stored as `B` features), clipped and inserted bases random.
fn mapped(tc: &TestCase, reference: &str, pos: u32) -> (String, String) {
    let n = |lo: u32, hi: u32| gs::integers::<u32>().min_value(lo).max_value(hi);
    let mut ops: Vec<(u32, char)> = Vec::new();
    if tc.draw_silent(gs::booleans()) {
        ops.push((tc.draw_silent(n(1, 5)), 'H'));
    }
    if tc.draw_silent(gs::booleans()) {
        ops.push((tc.draw_silent(n(1, 6)), 'S'));
    }
    ops.push((tc.draw_silent(n(4, 30)), 'M'));
    for _ in 0..tc.draw_silent(gs::integers::<usize>().max_value(2)) {
        ops.push((tc.draw_silent(n(1, 5)), tc.draw_silent(gs::sampled_from(&['I', 'D', 'N']))));
        ops.push((tc.draw_silent(n(4, 30)), 'M'));
    }
    if tc.draw_silent(gs::booleans()) {
        ops.push((tc.draw_silent(n(1, 6)), 'S'));
    }
    let span: u32 = ops.iter().filter(|(_, op)| matches!(op, 'M' | 'D' | 'N')).map(|o| o.0).sum();
    let start = (pos as usize - 1).min((CONTIG_LEN - span) as usize);
    let mut seq = Vec::new();
    let mut at = start;
    for &(len, op) in &ops {
        let len = len as usize;
        match op {
            'M' => {
                seq.extend_from_slice(&reference.as_bytes()[at..at + len]);
                at += len;
            }
            'D' | 'N' => at += len,
            'I' | 'S' => seq.extend_from_slice(tc.draw_silent(bases(len)).as_bytes()),
            _ => {}
        }
    }
    for _ in 0..tc.draw_silent(gs::integers::<usize>().max_value(3)) {
        let i = tc.draw_silent(gs::integers::<usize>().max_value(seq.len() - 1));
        seq[i] = tc.draw_silent(gs::sampled_from(b"ACGTNR"));
    }
    let cigar = ops.iter().fold(String::new(), |mut s, (l, o)| {
        write!(s, "{l}{o}").unwrap();
        s
    });
    (cigar, String::from_utf8(seq).unwrap())
}

#[hegel::composite]
fn arb_tags(tc: &TestCase) -> Vec<String> {
    let n = |lo: u32, hi: u32| gs::integers::<u32>().min_value(lo).max_value(hi);
    let mut tags = Vec::new();
    for tag in PROBE_TAGS {
        if !tc.draw_silent(gs::booleans()) {
            continue;
        }
        tags.push(match &tag {
            b"NM" => format!("NM:i:{}", tc.draw_silent(n(0, 300))),
            b"MD" => format!("MD:Z:{}", tc.draw_silent(n(0, 99))),
            b"RG" => format!("RG:Z:g{}", tc.draw_silent(n(1, 2))),
            b"XA" => format!("XA:Z:{}", tc.draw_silent(bases(3))),
            b"XF" => format!("XF:f:{}.5", tc.draw_silent(n(0, 9))),
            b"XB" => format!("XB:B:c,{},-{}", tc.draw_silent(n(0, 99)), tc.draw_silent(n(0, 99))),
            b"MC" => format!("MC:Z:{}M", tc.draw_silent(n(1, 99))),
            b"AS" => format!("AS:i:-{}", tc.draw_silent(n(0, 70_000))),
            _ => continue,
        });
    }
    tags
}

#[hegel::composite]
fn arb_input(tc: &TestCase) -> (Vec<String>, Vec<Read>) {
    let n = |lo: u32, hi: u32| gs::integers::<u32>().min_value(lo).max_value(hi);
    let references: Vec<String> =
        CONTIGS.iter().map(|_| tc.draw_silent(bases(CONTIG_LEN as usize))).collect();
    let mut reads = Vec::new();
    for t in 0..tc.draw_silent(gs::integers::<usize>().min_value(1).max_value(6)) {
        let pair = tc.draw_silent(gs::booleans());
        let mut template = Vec::new();
        for end in 0..if pair { 2 } else { 1 } {
            let contig = tc.draw_silent(gs::sampled_from(&[Some(0), Some(0), Some(1), None]));
            let pos = if contig.is_some() { tc.draw_silent(n(1, CONTIG_LEN - 60)) } else { 0 };
            let unmapped =
                contig.is_none() || tc.draw_silent(gs::integers::<u8>().max_value(5)) == 0;
            let (cigar, seq) = match contig {
                Some(c) if !unmapped => mapped(tc, &references[c], pos),
                _ => {
                    let len = tc.draw_silent(n(5, 30)) as usize;
                    ("*".to_owned(), tc.draw_silent(bases(len)))
                }
            };
            let qual = match tc.draw_silent(gs::integers::<u8>().max_value(3)) {
                0 => "*".to_owned(),
                _ => tc.draw_silent(
                    gs::text().alphabet("#5?I").min_size(seq.len()).max_size(seq.len()),
                ),
            };
            let seq_star = !unmapped && tc.draw_silent(gs::integers::<u8>().max_value(7)) == 0;
            let mut flags =
                if unmapped { 4 } else { tc.draw_silent(gs::sampled_from(&[0u32, 16])) };
            if pair {
                flags |= 1 | if end == 0 { 0x40 } else { 0x80 };
            }
            template.push(Read {
                name: format!("t{t}"),
                flags,
                contig,
                pos,
                cigar,
                seq: if seq_star { "*".to_owned() } else { seq },
                qual: if seq_star { "*".to_owned() } else { qual },
                mate: None,
                tags: tc.draw_silent(arb_tags()),
            });
        }
        if let [a, b] = &mut template[..]
            && let (Some(ca), Some(cb)) = (a.contig, b.contig)
        {
            a.mate = Some((cb, b.pos));
            b.mate = Some((ca, a.pos));
        }
        reads.extend(template);
    }
    (references, reads)
}

fn sam_text(reads: &[Read]) -> String {
    let mut sam = "@HD\tVN:1.6\tSO:unsorted\n".to_owned();
    for c in CONTIGS {
        writeln!(sam, "@SQ\tSN:{c}\tLN:{CONTIG_LEN}").unwrap();
    }
    sam.push_str("@RG\tID:g1\n@RG\tID:g2\n");
    for r in reads {
        let rname = r.contig.map_or("*", |c| CONTIGS[c]);
        let (rnext, pnext) = match r.mate {
            Some((c, p)) if Some(c) == r.contig => ("=", p),
            Some((c, p)) => (CONTIGS[c], p),
            None => ("*", 0),
        };
        let mapq = if r.flags & 4 == 0 { 30 } else { 0 };
        write!(
            sam,
            "{}\t{}\t{rname}\t{}\t{mapq}\t{}\t{rnext}\t{pnext}\t0\t{}\t{}",
            r.name, r.flags, r.pos, r.cigar, r.seq, r.qual
        )
        .unwrap();
        for t in &r.tags {
            write!(sam, "\t{t}").unwrap();
        }
        sam.push('\n');
    }
    sam
}

/// samtools options that change how the CRAM is laid out.
const CRAM_OPTIONS: [&[&str]; 7] = [
    &[],
    &["version=3.0"],
    &["seqs_per_slice=3", "slices_per_container=2"],
    &["multi_seq_per_slice=1", "seqs_per_slice=8"],
    &["lossy_names=1"],
    &["no_ref=1"],
    &["store_md=1", "store_nm=1"],
];

struct Files {
    _dir: tempfile::TempDir,
    fasta: PathBuf,
    cram: PathBuf,
}

fn write_files(references: &[String], reads: &[Read], options: &[&str]) -> Files {
    let dir = tempfile::tempdir().unwrap();
    let p = |name: &str| dir.path().join(name);
    let fasta = p("ref.fa");
    let mut f = std::fs::File::create(&fasta).unwrap();
    for (c, bases) in CONTIGS.iter().zip(references) {
        writeln!(f, ">{c}\n{bases}").unwrap();
    }
    drop(f);
    run(Command::new("samtools").arg("faidx").arg(&fasta), "samtools faidx");
    let sam = p("raw.sam");
    std::fs::write(&sam, sam_text(reads)).unwrap();
    let bam = p("g.bam");
    run(Command::new("samtools").args(["sort", "-o"]).arg(&bam).arg(&sam), "samtools sort");
    let cram = p("g.cram");
    let mut view = Command::new("samtools");
    view.args(["view", "--no-PG", "-C", "-T"]).arg(&fasta);
    for o in options {
        view.args(["--output-fmt-option", o]);
    }
    run(view.arg("-o").arg(&cram).arg(&bam), "samtools view -C");
    run(Command::new("samtools").arg("index").arg(&cram), "samtools index");
    Files { _dir: dir, fasta, cram }
}

fn whole_contig(tid: Tid) -> (Tid, std::range::RangeInclusive<Pos0>) {
    (tid, (Pos0::new(0).unwrap()..=Pos0::new(CONTIG_LEN - 1).unwrap()).into())
}

/// The store's records as tuples, tags cut down to `aux`.
fn stored(store: &RecordStore, aux: &AuxTags) -> Vec<String> {
    store
        .indices()
        .map(|i| {
            let r = store.record(i).unwrap();
            let tags: Tags = tags(seqair::bam::aux::iter_tags(r.aux()))
                .into_iter()
                .filter(|(t, _)| aux.contains(*t))
                .collect();
            format!(
                "{:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {tags:?}",
                r.pos,
                r.end_pos,
                r.flags,
                r.mapq,
                r.qname(),
                r.cigar(),
                r.seq(),
                r.qual()
            )
        })
        .collect()
}

// r[verify cram.stream.fields]
// r[verify cram.stream.fields.dropped_blocks]
// r[verify cram.stream.fields.seq_off]
// r[verify cram.stream.fields.qual_off]
// r[verify cram.stream.fields.md_nm]
// r[verify cram.stream.fields.reference_optional]
// r[verify unified.record_fields.indexed]
// r[verify unified.record_fields.fetch_into]
// r[verify unified.stream.detect+4]
/// For any fields, every requested field of the stream's, `query_raw`'s and
/// `fetch_into`'s records is the full decode's, over reads with mates,
/// unmapped and unplaced reads, SEQ `*`, QUAL `*`, `B` features, clips,
/// skips and tags, in files laid out every way samtools can; and a reader
/// opens without a reference exactly when its fields need none.
#[hegel::test(test_cases = crate::pinned::cases(32))]
fn requested_fields_are_the_full_decodes(tc: TestCase) {
    let (references, reads) = tc.draw(arb_input().print_as_debug());
    let options = tc.draw(gs::sampled_from(&CRAM_OPTIONS));
    let fields = tc.draw(arb_fields());
    let decode_md = tc.draw(gs::booleans());
    let files = write_files(&references, &reads, options);
    let (cram, fasta) = (files.cram.as_path(), files.fasta.as_path());
    let needs = needs_reference(&fields);

    // The stream.
    let full = whole(stream(cram, Some(fasta), decode_md, None).unwrap());
    assert_eq!(full.len(), reads.len(), "stream records");
    let got = stream(cram, Some(fasta), decode_md, Some(&fields)).unwrap();
    check_all(&got, &full, &fields, "stream");
    if !fields.sequence() && !(decode_md && needs) {
        let got = stream(cram, None, decode_md, Some(&fields)).unwrap();
        check_all(&got, &full, &fields, "stream without a reference");
    }

    // `AlignmentStream`, which generates MD/NM, needs a reference exactly
    // when the fields do.
    let opened =
        AlignmentStream::builder(std::fs::File::open(cram).unwrap()).fields(fields.clone()).build();
    match opened {
        Err(AlignmentStreamError::CramNeedsReference) => assert!(needs),
        Err(e) => panic!("{e}"),
        Ok(mut s) => {
            assert!(!needs);
            let mut got = Vec::new();
            let mut r = PartialRecord::default();
            while s.read_record(&mut r).unwrap() {
                got.push(r.clone());
            }
            let full = whole(stream(cram, Some(fasta), true, None).unwrap());
            check_all(&got, &full, &fields, "AlignmentStream without a reference");
        }
    }

    // `query_raw` and `fetch_into`, with and without a reference.
    let mut whole_reader = IndexedReader::builder(cram).reference(fasta).open().unwrap();
    let mut partial =
        IndexedReader::builder(cram).reference(fasta).fields(fields.clone()).open().unwrap();
    let mut bare = match IndexedReader::builder(cram).fields(fields.clone()).open() {
        Err(ReaderError::Format { source: FormatDetectionError::CramRequiresFasta }) => {
            assert!(needs);
            None
        }
        Err(e) => panic!("{e}"),
        Ok(r) => {
            assert!(!needs);
            Some(r)
        }
    };
    for contig in CONTIGS {
        let tid = seqair::reader::ResolveTid::resolve_tid(&contig, whole_reader.header()).unwrap();
        let (tid, span) = whole_contig(tid);
        let full: Vec<RawBamRecord> =
            whole_reader.query_raw(tid, span).unwrap().map(Result::unwrap).collect();
        let got: Vec<PartialRecord> =
            partial.query_raw(tid, span).unwrap().map(Result::unwrap).collect();
        check_all(&got, &full, &fields, "query_raw");
        if let Some(bare) = &mut bare {
            let got: Vec<PartialRecord> =
                bare.query_raw(tid, span).unwrap().map(Result::unwrap).collect();
            check_all(&got, &full, &fields, "query_raw without a reference");
        }

        let mut want = RecordStore::new();
        whole_reader.fetch_into(tid, span, &mut want).unwrap();
        let mut store = RecordStore::new();
        partial.fetch_into(tid, span, &mut store).unwrap();
        assert_eq!(stored(&store, &AuxTags::All), stored(&want, fields.aux()), "fetch_into");
    }

    tc.event_value("reads", reads.len() as f64);
    if options.contains(&"multi_seq_per_slice=1") {
        tc.event("multi-reference slices");
    }
}
