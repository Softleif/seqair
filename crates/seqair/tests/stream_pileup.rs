//! The streaming pileup against the indexed one: for generated
//! coordinate-sorted inputs, `StreamPileup` over the unindexed stream yields
//! the columns `Readers::pileup` yields over one whole-contig segment of the
//! same file, indexed by samtools. Plus the sort check, and depth against
//! `samtools mpileup` on the repo's test BAM.
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
use seqair::Readers;
use seqair::bam::pileup::{Indel, PileupColumn, PileupOp};
use seqair::bam::record_store::{CustomizeRecordStore, FilterRawFields, RecordStore, SlimRecord};
use seqair::fasta::IndexedFastaReader;
use seqair::reader::{
    AlignmentStream, DepthLimit, SegmentOptions, StreamPileup, StreamPileupError, WindowBudget,
};
use seqair_types::Base;
use std::fmt::Write as _;
use std::num::NonZeroU32;
use std::path::Path;
use std::process::Command;

/// Contigs of the generated inputs: long enough for several windows, short
/// enough that a long deletion or intron runs off the end of one.
const CONTIGS: [(&str, u32); 3] = [("chrA", 3000), ("chrB", 1200), ("chrC", 500)];

/// One SAM line and where it sorts.
#[derive(Debug, Clone)]
struct Line {
    /// `(tid, pos)`, `None` for an unplaced record.
    key: Option<(usize, u32)>,
    text: String,
}

/// A CIGAR and the query length it consumes. Mostly short ops; sometimes a
/// deletion or intron long enough to cross many windows and the contig end;
/// sometimes `D` then `I`, the complex-indel case.
#[hegel::composite]
fn cigar(tc: &TestCase) -> (String, u32) {
    let mut cigar = String::new();
    let mut qlen = 0;
    let mut push = |cigar: &mut String, len: u32, op: char, query: bool| {
        write!(cigar, "{len}{op}").unwrap();
        if query {
            qlen += len;
        }
    };
    let clip = |tc: &TestCase| tc.draw_silent(gs::integers::<u32>().max_value(4));
    let lead = clip(tc);
    if lead > 0 {
        push(&mut cigar, lead, 'S', true);
    }
    let matches = |tc: &TestCase| {
        let op = tc.draw_silent(gs::sampled_from(vec!['M', 'M', 'M', '=', 'X']));
        (tc.draw_silent(gs::integers::<u32>().min_value(1).max_value(40)), op)
    };
    let (len, op) = matches(tc);
    push(&mut cigar, len, op, true);
    for _ in 0..tc.draw_silent(gs::integers::<u32>().max_value(3)) {
        match tc.draw_silent(gs::integers::<u8>().max_value(9)) {
            0..=2 => push(
                &mut cigar,
                tc.draw_silent(gs::integers().min_value(1).max_value(5)),
                'I',
                true,
            ),
            3..=5 => push(
                &mut cigar,
                tc.draw_silent(gs::integers().min_value(1).max_value(12)),
                'D',
                false,
            ),
            6 => push(
                &mut cigar,
                tc.draw_silent(gs::integers().min_value(1).max_value(60)),
                'N',
                false,
            ),
            7 => {
                push(
                    &mut cigar,
                    tc.draw_silent(gs::integers().min_value(1).max_value(4)),
                    'D',
                    false,
                );
                push(
                    &mut cigar,
                    tc.draw_silent(gs::integers().min_value(1).max_value(4)),
                    'I',
                    true,
                );
            }
            _ => {
                let op = tc.draw_silent(gs::sampled_from(vec!['D', 'N']));
                push(
                    &mut cigar,
                    tc.draw_silent(gs::integers().min_value(200).max_value(2500)),
                    op,
                    false,
                );
            }
        }
        let (len, op) = matches(tc);
        push(&mut cigar, len, op, true);
    }
    let trail = clip(tc);
    if trail > 0 {
        push(&mut cigar, trail, 'S', true);
    }
    (cigar, qlen)
}

fn bases(tc: &TestCase, len: u32) -> String {
    tc.draw_silent(gs::text().alphabet("ACGTN").min_size(len as usize).max_size(len as usize))
}

fn quals(tc: &TestCase, len: u32) -> String {
    if tc.draw_silent(gs::weighted_booleans(0.1)) {
        return "*".into();
    }
    tc.draw_silent(gs::text().alphabet("#+5?I").min_size(len as usize).max_size(len as usize))
}

/// A mapped read's fields after QNAME and FLAG: everything from RNAME on,
/// with the mate fields given.
fn mapped_fields(tc: &TestCase, tid: usize, pos: u32, mate: &str) -> String {
    let (cigar, qlen) = tc.draw_silent(cigar());
    let (seq, qual) = if tc.draw_silent(gs::weighted_booleans(0.05)) {
        ("*".to_owned(), "*".to_owned())
    } else {
        (bases(tc, qlen), quals(tc, qlen))
    };
    let mapq = tc.draw_silent(gs::integers::<u8>().max_value(60));
    format!("{}\t{}\t{mapq}\t{cigar}\t{mate}\t0\t{seq}\t{qual}", CONTIGS[tid].0, pos + 1)
}

/// Every record of a coordinate-sorted file: singles, overlapping and
/// distant pairs, pairs with a placed-unmapped mate, secondary alignments,
/// and unplaced reads at the end.
#[hegel::composite]
fn sorted_lines(tc: &TestCase) -> Vec<Line> {
    let mut lines = Vec::new();
    for (tid, &(_, len)) in CONTIGS.iter().enumerate() {
        let templates = tc.draw_silent(gs::integers::<usize>().max_value(14));
        for t in 0..templates {
            let qname = format!("t{tid}_{t}");
            let pos = tc.draw_silent(gs::integers::<u32>().max_value(len - 1));
            let reverse = if tc.draw_silent(gs::booleans()) { 0x10 } else { 0 };
            match tc.draw_silent(gs::integers::<u8>().max_value(5)) {
                0 | 1 => {
                    let secondary =
                        if tc.draw_silent(gs::weighted_booleans(0.15)) { 0x100 } else { 0 };
                    let fields = mapped_fields(tc, tid, pos, "*\t0");
                    lines.push(Line {
                        key: Some((tid, pos)),
                        text: format!("{qname}\t{}\t{fields}", reverse | secondary),
                    });
                }
                2..=4 => {
                    // Mostly near enough to overlap; sometimes far apart.
                    let offset = if tc.draw_silent(gs::weighted_booleans(0.8)) {
                        tc.draw_silent(gs::integers::<u32>().max_value(80))
                    } else {
                        tc.draw_silent(gs::integers::<u32>().max_value(2000))
                    };
                    let mate_pos = (pos + offset).min(len - 1);
                    let mate_reverse = if reverse == 0 { 0x10 } else { 0 };
                    let first = mapped_fields(tc, tid, pos, &format!("=\t{}", mate_pos + 1));
                    let second = mapped_fields(tc, tid, mate_pos, &format!("=\t{}", pos + 1));
                    lines.push(Line {
                        key: Some((tid, pos)),
                        text: format!("{qname}\t{}\t{first}", 0x41 | reverse | mate_reverse << 1),
                    });
                    lines.push(Line {
                        key: Some((tid, mate_pos)),
                        text: format!("{qname}\t{}\t{second}", 0x81 | mate_reverse | reverse << 1),
                    });
                }
                _ => {
                    // A mapped read whose mate is unmapped and placed beside it.
                    let mapped = mapped_fields(tc, tid, pos, &format!("=\t{}", pos + 1));
                    lines.push(Line {
                        key: Some((tid, pos)),
                        text: format!("{qname}\t{}\t{mapped}", 0x49 | reverse),
                    });
                    let seq = bases(tc, 8);
                    lines.push(Line {
                        key: Some((tid, pos)),
                        text: format!(
                            "{qname}\t{}\t{}\t{}\t0\t*\t=\t{}\t0\t{seq}\t*",
                            0x85 | reverse << 1,
                            CONTIGS[tid].0,
                            pos + 1,
                            pos + 1
                        ),
                    });
                }
            }
        }
    }
    // Stable: records at one position keep the order they were made in.
    lines.sort_by_key(|l| l.key);
    for u in 0..tc.draw_silent(gs::integers::<usize>().max_value(3)) {
        let seq = bases(tc, 6);
        lines.push(Line { key: None, text: format!("u{u}\t4\t*\t0\t0\t*\t*\t0\t0\t{seq}\t*") });
    }
    lines
}

fn sam_text(lines: &[Line]) -> String {
    let mut text = String::from("@HD\tVN:1.6\tSO:coordinate\n");
    for (name, len) in CONTIGS {
        writeln!(text, "@SQ\tSN:{name}\tLN:{len}").unwrap();
    }
    for line in lines {
        text.push_str(&line.text);
        text.push('\n');
    }
    text
}

/// A reference for [`CONTIGS`] from a seed, so it costs one draw.
fn fasta_text(seed: u64) -> String {
    let mut state = seed | 1;
    let mut text = String::new();
    for (name, len) in CONTIGS {
        writeln!(text, ">{name}").unwrap();
        for i in 0..len {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            text.push(b"ACGTN"[(state % 5) as usize] as char);
            if i % 60 == 59 {
                text.push('\n');
            }
        }
        text.push('\n');
    }
    text
}

fn samtools(args: &[&str], dir: &Path) {
    let out = Command::new("samtools").args(args).current_dir(dir).output().expect("samtools");
    assert!(out.status.success(), "samtools {args:?}: {}", String::from_utf8_lossy(&out.stderr));
}

/// One alignment of a column, by what it is rather than by store index.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Entry {
    qname: Vec<u8>,
    flags: u16,
    pos: u32,
    op: PileupOp,
    indel_after: Indel,
    in_mate_overlap: bool,
    /// The mate's flags, when the mate is in this column.
    mate: Option<u16>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Col {
    tid: u32,
    pos: u32,
    depth: usize,
    reference: Base,
    matches_reference: bool,
    entries: Vec<Entry>,
}

fn col(tid: u32, column: &PileupColumn<'_, ()>) -> Col {
    let entries = column
        .alignments()
        .map(|a| Entry {
            qname: a.qname().to_vec(),
            flags: a.flags.raw(),
            pos: a.record().pos.as_u32(),
            op: *a.op(),
            indel_after: a.indel_after(),
            in_mate_overlap: a.in_mate_overlap(),
            mate: column.mate_of(&a).map(|m| m.flags.raw()),
        })
        .collect();
    Col {
        tid,
        pos: column.pos().as_u32(),
        depth: column.depth(),
        reference: column.reference_base(),
        matches_reference: column.matches_reference(),
        entries,
    }
}

/// The oracle: one whole-contig segment per contig through the indexed path.
fn indexed_columns(bam: &Path, fasta: &Path, depth: DepthLimit) -> Vec<Col> {
    let mut readers = Readers::open(bam, fasta).unwrap();
    let opts = SegmentOptions::new(NonZeroU32::new(1 << 30).unwrap());
    let segments: Vec<_> = readers.segments((), opts).unwrap().collect();
    let mut cols = Vec::new();
    for segment in &segments {
        assert!(segment.starts_at_contig_start() && segment.ends_at_contig_end());
        let mut pileup = readers.pileup(segment, depth).run().unwrap();
        while let Some(column) = pileup.pileups() {
            cols.push(col(segment.tid().as_u32(), &column));
        }
    }
    cols
}

fn streamed_columns<E: CustomizeRecordStore<Extra = ()>>(mut pileup: StreamPileup<E>) -> Vec<Col> {
    let mut cols = Vec::new();
    while let Some((tid, column)) = pileup.pileups().unwrap() {
        cols.push(col(tid.as_u32(), &column));
    }
    cols
}

// r[verify pileup_stream.columns]
// r[verify pileup_stream.window+2]
// r[verify pileup_stream.carry]
// r[verify pileup_stream.window_len]
// r[verify pileup_stream.budget]
// r[verify pileup_stream.records+1]
// r[verify pileup_stream.mates+1]
// r[verify pileup_stream.depth_limit+1]
// r[verify pileup_stream.reference]
// r[verify pileup_stream.lending]
/// Every column, and the set of positions, matches the indexed pileup of the
/// same file — whatever the window length, including windows far shorter
/// than the reads, whatever the window budget, down to one record, from BAM or
/// SAM text, with or without a depth cap.
#[hegel::test(test_cases = pinned::cases(60))]
fn streaming_columns_equal_indexed_columns(tc: TestCase) {
    let lines = tc.draw(sorted_lines());
    let seed = tc.draw(gs::integers::<u64>());
    let window = if tc.draw(gs::booleans()) {
        tc.draw(gs::integers::<u32>().min_value(1).max_value(300))
    } else {
        seqair::reader::DEFAULT_WINDOW_LEN.get()
    };
    let depth = match tc.draw(gs::integers::<u32>().max_value(4)) {
        0 => DepthLimit::PerColumn(
            NonZeroU32::new(tc.draw(gs::integers().min_value(1).max_value(4))).unwrap(),
        ),
        _ => DepthLimit::Unlimited,
    };
    let budget = if tc.draw(gs::booleans()) {
        let small = |tc: &TestCase, max| {
            NonZeroU32::new(tc.draw(gs::integers::<u32>().min_value(1).max_value(max))).unwrap()
        };
        WindowBudget { records: small(&tc, 8), slab_len: small(&tc, 400) }
    } else {
        WindowBudget::default()
    };
    let from_sam = tc.draw(gs::booleans());

    let dir = tempfile::tempdir().unwrap();
    let sam = sam_text(&lines);
    std::fs::write(dir.path().join("in.sam"), &sam).unwrap();
    std::fs::write(dir.path().join("ref.fa"), fasta_text(seed)).unwrap();
    samtools(&["view", "-b", "-o", "in.bam", "in.sam"], dir.path());
    samtools(&["index", "in.bam"], dir.path());
    samtools(&["faidx", "ref.fa"], dir.path());
    let bam = dir.path().join("in.bam");
    let fasta = dir.path().join("ref.fa");

    let expected = indexed_columns(&bam, &fasta, depth);

    let stream = if from_sam {
        AlignmentStream::from_reader(std::io::Cursor::new(sam.into_bytes())).unwrap()
    } else {
        AlignmentStream::from_reader(std::fs::File::open(&bam).unwrap()).unwrap()
    };
    let pileup = StreamPileup::new(stream, (), depth)
        .with_window_len(NonZeroU32::new(window).unwrap())
        .with_window_budget(budget)
        .with_reference(IndexedFastaReader::open(&fasta).unwrap());
    let actual = streamed_columns(pileup);

    let positions = |cols: &[Col]| cols.iter().map(|c| (c.tid, c.pos)).collect::<Vec<_>>();
    assert_eq!(positions(&actual), positions(&expected), "column positions");
    for (a, e) in actual.iter().zip(&expected) {
        assert_eq!(a, e, "column {}:{}", a.tid, a.pos);
    }
}

fn is_mapped(line: &Line) -> bool {
    let flag: u16 = line.text.split('\t').nth(1).unwrap().parse().unwrap();
    flag & 0x4 == 0
}

// r[verify pileup_stream.sorted+1]
// r[verify pileup_stream.errors]
/// A mapped record moved ahead of a mapped one that sorts before it fails
/// with `Unsorted` at exactly the record after it, and the pileup stays failed.
#[hegel::test(test_cases = pinned::cases(100))]
fn unsorted_input_is_a_typed_error(tc: TestCase) {
    let mut lines = tc.draw(sorted_lines());
    // Bracket the placed records with the smallest and the largest key, so
    // the largest can always be moved in front of a smaller mapped one.
    let last_c = CONTIGS[2].1 - 1;
    lines.insert(
        0,
        Line { key: Some((0, 0)), text: "first\t0\tchrA\t1\t0\t1M\t*\t0\t0\tA\t*".into() },
    );
    let placed = lines.iter().take_while(|l| l.key.is_some()).count();
    let max = Line {
        key: Some((2, last_c)),
        text: format!("last\t0\tchrC\t{}\t0\t1M\t*\t0\t0\tA\t*", last_c + 1),
    };
    lines.insert(placed, max.clone());
    let smaller: Vec<usize> =
        (0..=placed).filter(|&i| lines[i].key < max.key && is_mapped(&lines[i])).collect();
    let at = tc.draw(gs::sampled_from(smaller));
    lines.remove(placed);
    lines.insert(at, max);
    let window = tc.draw(gs::integers::<u32>().min_value(1).max_value(100));

    let stream =
        AlignmentStream::from_reader(std::io::Cursor::new(sam_text(&lines).into_bytes())).unwrap();
    let mut pileup = StreamPileup::new(stream, (), DepthLimit::Unlimited)
        .with_window_len(NonZeroU32::new(window).unwrap());
    let err = loop {
        match pileup.pileups() {
            Ok(Some(_)) => {}
            Ok(None) => panic!("an unsorted input ended without an error"),
            Err(e) => break e,
        }
    };
    let StreamPileupError::Unsorted { record, tid, previous_tid, previous_pos, .. } = err else {
        panic!("expected Unsorted, got {err:?}");
    };
    // 1-based: the moved record is number `at + 1`, the one after it fails.
    assert_eq!(record, at as u64 + 2);
    assert_eq!(previous_tid.as_u32(), 2);
    assert_eq!(previous_pos.map(|p| p.as_u32()), Some(last_c));
    assert!(tid <= previous_tid);
    assert!(matches!(pileup.pileups(), Err(StreamPileupError::Poisoned)));
}

/// Keeps every record but secondary alignments, like a read callback that
/// filters on flags.
#[derive(Clone)]
struct RejectSecondary;

impl CustomizeRecordStore for RejectSecondary {
    type Extra = ();
    fn filter_raw(&mut self, fields: &FilterRawFields<'_>) -> bool {
        !fields.flags.is_secondary()
    }
    fn compute(&mut self, _: &SlimRecord, _: &RecordStore<()>) {}
}

// r[verify pileup_stream.sorted+1]
// r[verify pileup_stream.records+1]
/// What `bam_plp_push` does not check, seqair does not either: unmapped
/// records and records the customize value rejects may sit anywhere, out of
/// order, and change no column.
#[hegel::test(test_cases = pinned::cases(100))]
fn unchecked_records_may_be_out_of_order(tc: TestCase) {
    let lines = tc.draw(sorted_lines());
    let window = tc.draw(gs::integers::<u32>().min_value(1).max_value(100));
    let mut with_strays = lines.clone();
    for i in 0..tc.draw(gs::integers::<usize>().min_value(1).max_value(6)) {
        let tid = tc.draw(gs::integers::<usize>().max_value(CONTIGS.len() - 1));
        let (name, len) = CONTIGS[tid];
        let pos = tc.draw(gs::integers::<u32>().min_value(1).max_value(len));
        let text = if tc.draw(gs::booleans()) {
            format!("s{i}\t4\t{name}\t{pos}\t0\t*\t*\t0\t0\tACGT\t*")
        } else {
            format!("s{i}\t256\t{name}\t{pos}\t0\t4M\t*\t0\t0\tACGT\t*")
        };
        let at = tc.draw(gs::integers::<usize>().max_value(with_strays.len()));
        with_strays.insert(at, Line { key: None, text });
    }
    let columns = |lines: &[Line]| {
        let stream =
            AlignmentStream::from_reader(std::io::Cursor::new(sam_text(lines).into_bytes()))
                .unwrap();
        streamed_columns(
            StreamPileup::new(stream, RejectSecondary, DepthLimit::Unlimited)
                .with_window_len(NonZeroU32::new(window).unwrap()),
        )
    };
    assert_eq!(columns(&with_strays), columns(&lines));
}

// r[verify pileup_stream.budget]
/// The budget does close windows: with ten disjoint reads and a budget of two
/// records, or of a bases slab two reads long, no column's store holds more
/// than two records — the third read each time closes the window instead.
#[test]
fn the_budget_bounds_the_window_store() {
    let lines: Vec<Line> = (0..10)
        .map(|i| Line {
            key: None,
            text: format!("r{i}\t0\tchrA\t{}\t30\t10M\t*\t0\t0\tACGTACGTAC\t*", 1 + 20 * i),
        })
        .collect();
    let max_store = |budget: WindowBudget| {
        let stream =
            AlignmentStream::from_reader(std::io::Cursor::new(sam_text(&lines).into_bytes()))
                .unwrap();
        let mut pileup =
            StreamPileup::new(stream, (), DepthLimit::Unlimited).with_window_budget(budget);
        let mut max = 0;
        while let Some((_, column)) = pileup.pileups().unwrap() {
            max = max.max(column.store().len());
        }
        max
    };
    let n = |v| NonZeroU32::new(v).unwrap();
    assert_eq!(max_store(WindowBudget::default()), 10);
    assert_eq!(max_store(WindowBudget { records: n(2), slab_len: n(1 << 20) }), 2);
    assert_eq!(max_store(WindowBudget { records: n(1 << 20), slab_len: n(25) }), 2);
}

// r[verify pileup_stream.budget]
/// A pile deeper than the budget at one position cannot be split, so it grows
/// past the budget and yields the columns an unbudgeted window yields —
/// whether the pile opens its window or starts inside one.
#[test]
fn a_pile_deeper_than_the_budget_grows_past_it() {
    for lead in [0, 3] {
        let mut lines: Vec<Line> = (0..lead)
            .map(|i| Line {
                key: None,
                text: format!("lead{i}\t0\tchrA\t10\t30\t20M\t*\t0\t0\t*\t*"),
            })
            .collect();
        lines.extend((0..40).map(|i| Line {
            key: None,
            text: format!("pile{i}\t0\tchrA\t20\t30\t{}M\t*\t0\t0\t*\t*", 5 + i % 7),
        }));
        let columns = |budget: WindowBudget| {
            let stream =
                AlignmentStream::from_reader(std::io::Cursor::new(sam_text(&lines).into_bytes()))
                    .unwrap();
            let pileup =
                StreamPileup::new(stream, (), DepthLimit::Unlimited).with_window_budget(budget);
            streamed_columns(pileup)
        };
        let three = NonZeroU32::new(3).unwrap();
        let tiny = WindowBudget { records: three, slab_len: NonZeroU32::new(1 << 20).unwrap() };
        let expected = columns(WindowBudget::default());
        assert_eq!(expected.iter().map(|c| c.depth).max(), Some(40 + lead));
        assert_eq!(columns(tiny), expected, "lead {lead}");
    }
}

/// `(contig, 1-based position, depth)` per line of `samtools mpileup`.
fn mpileup_depths(bam: &Path) -> Vec<(String, u32, usize)> {
    // No BAQ, no quality or flag filters, orphans kept, no overlap
    // adjustment, no depth cap: every alignment htslib's pileup has.
    let out = Command::new("samtools")
        .args(["mpileup", "-B", "-Q", "0", "-q", "0", "-A", "-x", "-d", "0", "--ff", "0"])
        .arg(bam)
        .output()
        .expect("samtools");
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    String::from_utf8(out.stdout)
        .unwrap()
        .lines()
        .map(|line| {
            let f: Vec<&str> = line.split('\t').collect();
            (f[0].to_owned(), f[1].parse().unwrap(), f[3].parse().unwrap())
        })
        .filter(|&(_, _, depth)| depth > 0)
        .collect()
}

/// On the repo's test BAM (three contigs and a tail of unplaced reads) the
/// streaming pileup reports the depth `samtools mpileup` reports, at every
/// position it reports one.
#[test]
fn depth_matches_samtools_mpileup() {
    let bam = Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data/test.bam"));
    let expected = mpileup_depths(bam);

    let stream = AlignmentStream::from_reader(std::fs::File::open(bam).unwrap()).unwrap();
    let names: Vec<String> = stream.header().target_names().map(str::to_owned).collect();
    let mut pileup = StreamPileup::new(stream, (), DepthLimit::Unlimited);
    let mut actual = Vec::new();
    while let Some((tid, column)) = pileup.pileups().unwrap() {
        let name = names[tid.as_usize()].clone();
        actual.push((name, column.pos().as_u32() + 1, column.depth()));
    }
    assert!(!expected.is_empty());
    assert_eq!(actual, expected);
}
