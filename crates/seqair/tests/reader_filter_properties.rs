//! Reader-level filtering: `fetch_into_customized` and the customize value
//! attached by `Readers::open_customized`.
//!
//! This is the surface rastair actually integrates against — it never pushes
//! into a `RecordStore` by hand, it hands the reader a `CustomizeRecordStore`
//! and takes what comes back — and it had no property. The record-store tests
//! pin the slab bookkeeping of a single rejected push; what they cannot see is
//! the reader: whether the filter is consulted once per record the region query
//! produced, whether a rejection removes exactly that record from the store,
//! and whether the surviving records are still the ones an unfiltered fetch
//! would have given, in the same order, with the same bytes.
//!
//! So the oracle here is the *unfiltered* fetch of the same region. A filtered
//! fetch must equal it, restricted to the records the predicate accepts —
//! field by field, not merely in count, because a rollback that truncates a
//! slab too far leaves the right number of records holding the wrong bases.
//!
//! The predicate comes from a small algebra (flag mask, MAPQ floor, reference
//! span, qname digest; all-of / any-of; negation) so that "accepts everything"
//! and "accepts nothing" are the empty conjunction and the empty disjunction
//! rather than special cases bolted on. The clauses are evaluated twice from
//! two different views of the same record: `filter_raw` sees the raw BAM header
//! and CIGAR bytes before any slab is touched, while the oracle sees the
//! decoded `SlimRecord` the reader produced. A field the reader mis-slices on
//! its way into `FilterRawFields` makes the two disagree.
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
use seqair::bam::record_store::{CustomizeRecordStore, FilterRawFields, RecordStore, SlimRecord};
use seqair::bam::writer::BamWriterBuilder;
use seqair::bam::{IndexedBamReader, Pos0};
use seqair::reader::{FetchCounts, Readers};
use seqair_types::{BamFlags, Base, BaseQuality};
use std::fmt::Write as _;
use std::path::{Path, PathBuf};

/// One contig, long enough to span several 16 kb BAI leaf bins. The region
/// query's own correctness is pinned by `bam_index_query_properties`; here the
/// region is only the supply of records, so every fetch takes the whole contig
/// and gets as many of them as possible.
const CONTIG: (&str, u32) = ("chr1", 40_000);

fn header() -> BamHeader {
    let (name, len) = CONTIG;
    let mut text = String::from("@HD\tVN:1.6\tSO:coordinate\n");
    writeln!(text, "@SQ\tSN:{name}\tLN:{len}").expect("writing to a String cannot fail");
    BamHeader::from_sam_text(&text).unwrap()
}

// ── the reads ───────────────────────────────────────────────────────────

/// A record as the generator describes it, before it becomes BAM bytes.
#[derive(Debug, Clone)]
struct Read {
    pos: u32,
    flags: u16,
    mapq: u8,
    /// `(len, op)`; empty for a placed-unmapped read.
    cigar: Vec<(u32, CigarOpType)>,
    seq_len: u32,
    /// An `NM` tag on some records, so the aux slab is non-empty and a
    /// rollback that forgets to truncate it cannot hide behind an empty one.
    nm: Option<u8>,
}

/// Positions cluster at both ends of the contig: uniform starts over 40 kb
/// leave the reads too sparse to share a bin, and a filter that is never asked
/// about two records in the same chunk proves little.
fn arb_pos(tc: &TestCase, span: u32) -> u32 {
    let (_, contig_len) = CONTIG;
    let last = contig_len.saturating_sub(span).max(1) - 1;
    let cluster = tc.draw_silent(gs::integers::<u32>().max_value(1));
    let jitter = tc.draw_silent(gs::integers::<u32>().max_value(800));
    (cluster * (contig_len / 2)).saturating_add(jitter).min(last)
}

#[hegel::composite]
fn arb_read(tc: &TestCase) -> Read {
    let n = |lo: u32, hi: u32| gs::integers::<u32>().min_value(lo).max_value(hi);
    // A spread that makes every clause split the set: paired/unpaired, both
    // strands, secondary and supplementary, and the unmapped bit — which
    // `r[bam.reader.unmapped_skipped]` says reaches the customize layer like
    // any other record rather than being dropped by the reader.
    let flags = tc.draw_silent(gs::sampled_from(&[0u16, 16, 99, 147, 83, 163, 0x100, 0x800]));
    let mapq = tc.draw_silent(gs::integers::<u8>().max_value(60));
    let nm = tc.draw_silent(gs::optional(gs::integers::<u8>().max_value(9)));

    // One in ten is placed but unmapped: it has a position and is indexed
    // there, but no alignment, so its reference span is the degenerate one.
    if !tc.draw_silent(gs::weighted_booleans(0.9)) {
        let seq_len = tc.draw_silent(n(1, 40));
        let pos = arb_pos(tc, 1);
        return Read { pos, flags: flags | 0x4, mapq, cigar: Vec::new(), seq_len, nm };
    }

    let mut cigar: Vec<(u32, CigarOpType)> = Vec::new();
    if tc.draw_silent(gs::booleans()) {
        cigar.push((tc.draw_silent(n(1, 5)), CigarOpType::SoftClip));
    }
    cigar.push((tc.draw_silent(n(1, 30)), CigarOpType::Match));
    for _ in 0..tc.draw_silent(gs::integers::<usize>().max_value(2)) {
        let (op, len) = match tc.draw_silent(gs::integers::<u8>().max_value(2)) {
            0 => (CigarOpType::Insertion, tc.draw_silent(n(1, 5))),
            1 => (CigarOpType::Deletion, tc.draw_silent(n(1, 8))),
            // A spliced read: a long span for very few query bases, which is
            // what makes the span clause disagree with a seq-length clause.
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
    let seq_len = cigar
        .iter()
        .filter(|(_, op)| {
            matches!(op, CigarOpType::Match | CigarOpType::Insertion | CigarOpType::SoftClip)
        })
        .map(|(len, _)| len)
        .sum();

    Read { pos: arb_pos(tc, span), flags, mapq, cigar, seq_len, nm }
}

/// Coordinate-sorted, which `BamWriter`'s index co-production requires.
#[hegel::composite]
fn arb_reads(tc: &TestCase) -> Vec<Read> {
    let mut reads = tc.draw_silent(gs::vecs(arb_read()).min_size(4).max_size(30));
    reads.sort_by_key(|r| r.pos);
    reads
}

/// The generated reads as a BAM with a co-produced BAI, plus a FASTA beside it
/// so the same directory can also be opened through `Readers`.
fn write_inputs(dir: &Path, reads: &[Read]) -> (PathBuf, PathBuf) {
    let header = header();
    let bam_path = dir.join("generated.bam");
    let mut writer =
        BamWriterBuilder::to_path(&bam_path, &header).write_index(true).build().unwrap();

    for (i, read) in reads.iter().enumerate() {
        // Distinct bases per record, so a slab that is truncated to the wrong
        // length shows up as the wrong sequence rather than the same one.
        let base = [Base::A, Base::C, Base::G, Base::T][i % 4];
        let mut aux = seqair::bam::aux_data::AuxData::new();
        if let Some(nm) = read.nm {
            aux.set_int(*b"NM", i64::from(nm)).expect("a small non-negative int");
        }
        let rec = OwnedBamRecord::builder(
            0,
            Some(Pos0::new(read.pos).unwrap()),
            format!("r{i}").into_bytes(),
        )
        .flags(BamFlags::from(read.flags))
        .mapq(read.mapq)
        .cigar(read.cigar.iter().map(|&(len, op)| CigarOp::new(op, len)).collect())
        .seq(vec![base; read.seq_len as usize])
        .qual(vec![BaseQuality::from_byte(30); read.seq_len as usize])
        .aux(aux)
        .build()
        .unwrap();
        writer.write(&rec).unwrap();
    }

    let (_inner, index_builder) = writer.finish().unwrap();
    let bai = std::fs::File::create(bam_path.with_extension("bam.bai")).unwrap();
    index_builder.expect("write_index(true)").write_bai(bai, header.target_count()).unwrap();

    (bam_path, write_fasta(dir))
}

/// A plain FASTA and `.fai` for the contig. `Readers` opens a reference even
/// for BAM, where nothing reads it back.
fn write_fasta(dir: &Path) -> PathBuf {
    let (name, len) = CONTIG;
    const LINE: u32 = 60;
    let mut text = format!(">{name}\n");
    for start in (0..len).step_by(LINE as usize) {
        let this_line = LINE.min(len - start);
        text.push_str(&"A".repeat(this_line as usize));
        text.push('\n');
    }
    let fasta_path = dir.join("reference.fa");
    std::fs::write(&fasta_path, &text).unwrap();
    std::fs::write(
        dir.join("reference.fa.fai"),
        format!("{name}\t{len}\t{}\t{LINE}\t{}\n", name.len() + 2, LINE + 1),
    )
    .unwrap();
    fasta_path
}

// ── the filter algebra ──────────────────────────────────────────────────

/// One clause. Each reads a different part of the record, so a conjunction of
/// several cuts the set along independent axes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Clause {
    /// Every bit of the mask is set.
    FlagsSet(u16),
    /// Every bit of the mask is clear.
    FlagsClear(u16),
    /// MAPQ at or above the floor.
    MinMapq(u8),
    /// Reference span (`end_pos - pos + 1`) at or above the floor.
    MinRefSpan(u32),
    /// The read's name digests to `rest` modulo `modulus` — a predicate with no
    /// relation to any field the reader computes, so it cannot accidentally
    /// agree with one.
    QnameDigest { modulus: u8, rest: u8 },
}

impl Clause {
    fn accepts(self, flags: u16, mapq: u8, span: u32, qname: &[u8]) -> bool {
        match self {
            Self::FlagsSet(mask) => flags & mask == mask,
            Self::FlagsClear(mask) => flags & mask == 0,
            Self::MinMapq(floor) => mapq >= floor,
            Self::MinRefSpan(floor) => span >= floor,
            Self::QnameDigest { modulus, rest } => {
                // Both position- and length-sensitive, so appending the
                // trailing NUL of a BAM name that was handed over unstripped
                // changes the digest under every modulus this draws. A plain
                // byte sum would not — `\0` adds nothing — and a `filter_raw`
                // shown the padded name would agree with the oracle by
                // accident.
                let digest = qname
                    .iter()
                    .fold(0u32, |acc, &b| acc.wrapping_mul(31).wrapping_add(u32::from(b)))
                    .wrapping_add(qname.len() as u32);
                (digest % u32::from(modulus)) as u8 == rest
            }
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Mode {
    All,
    Any,
}

/// A predicate over records. The empty `All` accepts everything and the empty
/// `Any` accepts nothing, so the two degenerate filters the API promises fall
/// out of the algebra instead of being written by hand.
#[derive(Debug, Clone, PartialEq)]
struct Predicate {
    mode: Mode,
    negate: bool,
    clauses: Vec<Clause>,
}

impl Predicate {
    fn accepts(&self, flags: u16, mapq: u8, span: u32, qname: &[u8]) -> bool {
        let held = match self.mode {
            Mode::All => self.clauses.iter().all(|c| c.accepts(flags, mapq, span, qname)),
            Mode::Any => self.clauses.iter().any(|c| c.accepts(flags, mapq, span, qname)),
        };
        held != self.negate
    }

    fn always() -> Self {
        Self { mode: Mode::All, negate: false, clauses: Vec::new() }
    }

    fn never() -> Self {
        Self { mode: Mode::Any, negate: false, clauses: Vec::new() }
    }

    /// The conjunction of two all-of predicates, which is again an all-of
    /// predicate.
    fn and(a: &[Clause], b: &[Clause]) -> Self {
        Self { mode: Mode::All, negate: false, clauses: [a, b].concat() }
    }
}

#[hegel::composite]
fn arb_clause(tc: &TestCase) -> Clause {
    // Bits that actually vary across the generated flags, plus 0x4, which is
    // the one the reader is documented not to act on by itself.
    let mask = || gs::sampled_from(&[0x1u16, 0x4, 0x10, 0x40, 0x80, 0x100, 0x800]);
    match tc.draw_silent(gs::integers::<u8>().max_value(4)) {
        0 => Clause::FlagsSet(tc.draw_silent(mask())),
        1 => Clause::FlagsClear(tc.draw_silent(mask())),
        2 => Clause::MinMapq(tc.draw_silent(gs::integers::<u8>().max_value(60))),
        3 => Clause::MinRefSpan(tc.draw_silent(gs::integers::<u32>().min_value(1).max_value(60))),
        _ => {
            let modulus = tc.draw_silent(gs::integers::<u8>().min_value(2).max_value(4));
            Clause::QnameDigest {
                modulus,
                rest: tc.draw_silent(gs::integers::<u8>().max_value(modulus - 1)),
            }
        }
    }
}

/// Up to three clauses: enough for the modes and the negation to matter,
/// few enough that a shrunk counterexample is readable.
fn arb_clauses() -> impl PrintableGenerator<Vec<Clause>> {
    gs::vecs(arb_clause()).max_size(3).print_as_debug()
}

#[hegel::composite]
fn arb_predicate(tc: &TestCase) -> Predicate {
    Predicate {
        mode: if tc.draw_silent(gs::booleans()) { Mode::All } else { Mode::Any },
        negate: tc.draw_silent(gs::booleans()),
        clauses: tc.draw_silent(arb_clauses()),
    }
}

// ── the customize value ─────────────────────────────────────────────────

/// Which of the two rejection hooks the predicate is wired to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Gate {
    /// `filter_raw`, before any slab is touched.
    Raw,
    /// `filter`, after the record is written and with a rollback.
    Late,
}

#[derive(Clone)]
struct ClauseFilter {
    pred: Predicate,
    gate: Gate,
}

impl ClauseFilter {
    fn raw(pred: Predicate) -> Self {
        Self { pred, gate: Gate::Raw }
    }

    fn late(pred: Predicate) -> Self {
        Self { pred, gate: Gate::Late }
    }
}

/// Reference span from raw BAM CIGAR bytes — the span the reader will later
/// report as `end_pos - pos + 1`, derived here from the other end.
///
/// A record whose CIGAR consumes no reference still covers exactly its own
/// position (`r[interval.end_pos_inclusive]`), so the floor is one.
fn ref_span_from_cigar_bytes(cigar: &[u8]) -> u32 {
    let span: u32 = cigar
        .as_chunks::<4>()
        .0
        .iter()
        .map(|op| u32::from_le_bytes(*op))
        .filter(|op| matches!(op & 0xf, 0 | 2 | 3 | 7 | 8))
        .map(|op| op >> 4)
        .sum();
    span.max(1)
}

impl CustomizeRecordStore for ClauseFilter {
    type Extra = ();

    fn filter_raw(&mut self, fields: &FilterRawFields<'_>) -> bool {
        if self.gate != Gate::Raw {
            return true;
        }
        self.pred.accepts(
            fields.flags.raw(),
            fields.mapq,
            ref_span_from_cigar_bytes(fields.cigar),
            fields.qname,
        )
    }

    fn filter(&mut self, rec: &SlimRecord, store: &RecordStore<()>) -> bool {
        if self.gate != Gate::Late {
            return true;
        }
        let qname = rec.qname(store).expect("a just-pushed record's name is in the slab");
        self.pred.accepts(rec.flags.raw(), rec.mapq, span_of(rec.pos, rec.end_pos), qname)
    }

    fn compute(&mut self, _rec: &SlimRecord, _store: &RecordStore<()>) {}
}

fn span_of(pos: Pos0, end_pos: Pos0) -> u32 {
    (end_pos.as_u64() - pos.as_u64() + 1) as u32
}

// ── comparing stores ────────────────────────────────────────────────────

/// Everything a `RecordStore` will say about one record: the fixed fields and
/// every slab it owns. Comparing only positions, or only the count, would pass
/// for a rollback that truncated the bases slab one base short.
#[derive(Debug, Clone, PartialEq)]
struct Snapshot {
    pos: u64,
    end_pos: u64,
    flags: u16,
    mapq: u8,
    n_cigar_ops: u16,
    seq_len: u32,
    matching_bases: u32,
    indel_bases: u32,
    tid: i32,
    next_ref_id: i32,
    next_pos: i32,
    template_len: i32,
    qname: Vec<u8>,
    seq: Vec<Base>,
    qual: Vec<u8>,
    cigar: Vec<CigarOp>,
    aux: Vec<u8>,
}

impl Snapshot {
    fn accepts(&self, pred: &Predicate) -> bool {
        let span = (self.end_pos - self.pos + 1) as u32;
        pred.accepts(self.flags, self.mapq, span, &self.qname)
    }
}

fn snapshot(store: &RecordStore) -> Vec<Snapshot> {
    store
        .indices()
        .map(|idx| {
            let rec = store.record(idx).unwrap();
            Snapshot {
                pos: rec.pos.as_u64(),
                end_pos: rec.end_pos.as_u64(),
                flags: rec.flags.raw(),
                mapq: rec.mapq,
                n_cigar_ops: rec.n_cigar_ops,
                seq_len: rec.seq_len,
                matching_bases: rec.matching_bases,
                indel_bases: rec.indel_bases,
                tid: rec.tid,
                next_ref_id: rec.next_ref_id,
                next_pos: rec.next_pos,
                template_len: rec.template_len,
                qname: rec.qname().to_vec(),
                seq: rec.seq().to_vec(),
                qual: BaseQuality::slice_to_bytes(rec.qual()).to_vec(),
                cigar: rec.cigar().to_vec(),
                aux: rec.aux().to_vec(),
            }
        })
        .collect()
}

/// The whole contig, which is every record in the file.
fn whole_contig() -> (Pos0, Pos0) {
    (Pos0::new(0u32).unwrap(), Pos0::new(CONTIG.1 - 1).unwrap())
}

/// Fetch with a customize value, returning both what the store holds and what
/// the reader counted.
fn fetch_filtered(bam: &Path, pred: Predicate, gate: Gate) -> (Vec<Snapshot>, FetchCounts) {
    let mut reader = IndexedBamReader::open(bam).expect("open");
    let tid = reader.header().tid(CONTIG.0).expect("contig");
    let (start, end) = whole_contig();
    let mut store = RecordStore::new();
    let counts = reader
        .fetch_into_customized(tid, start, end, &mut store, &mut ClauseFilter { pred, gate })
        .expect("fetch_into_customized");
    (snapshot(&store), counts)
}

/// Fetch through the plain `fetch_into`, which passes `&mut ()`.
fn fetch_plain(bam: &Path) -> (Vec<Snapshot>, usize) {
    let mut reader = IndexedBamReader::open(bam).expect("open");
    let tid = reader.header().tid(CONTIG.0).expect("contig");
    let (start, end) = whole_contig();
    let mut store = RecordStore::new();
    let kept = reader.fetch_into(tid, start, end, &mut store).expect("fetch_into");
    (snapshot(&store), kept)
}

/// Fetch through `Readers`, which carries the customize value from open time
/// rather than taking it per call.
fn fetch_through_readers(bam: &Path, fasta: &Path, pred: Predicate) -> Vec<Snapshot> {
    let mut readers =
        Readers::open_customized(bam, fasta, ClauseFilter::raw(pred)).expect("open_customized");
    let tid = readers.header().tid(CONTIG.0).expect("contig");
    let (start, end) = whole_contig();
    let mut store = RecordStore::new();
    readers.fetch_into(tid, start, end, &mut store).expect("Readers::fetch_into");
    snapshot(&store)
}

/// Report how hard the predicate cut, so a suite that only ever generated
/// accept-everything or accept-nothing predicates would say so.
fn report_selectivity(tc: &TestCase, kept: usize, total: usize) {
    tc.event_value("fetched", total as f64);
    tc.event_value("kept", kept as f64);
    if kept > 0 && kept < total {
        tc.event("the filter split the region");
    } else if kept == 0 && total > 0 {
        tc.event("the filter rejected everything");
    } else if total > 0 {
        tc.event("the filter accepted everything");
    }
}

// ── properties ──────────────────────────────────────────────────────────

// r[verify unified.fetch_into_customized]
// r[verify unified.readers_open_customized]
// r[verify record_store.pre_filter.rollback]
// r[verify record_store.filter_raw]
// r[verify record_store.filter_raw_fields]
/// A filtered fetch is the unfiltered fetch restricted to the records the
/// predicate accepts — same records, same order, same bytes — and the same
/// predicate attached at open time through `Readers::open_customized` gives
/// the same store.
#[hegel::test]
fn a_filtered_fetch_is_the_unfiltered_fetch_restricted(tc: TestCase) {
    let reads = tc.draw(arb_reads().print_as_debug());
    let pred = tc.draw(arb_predicate().print_as_debug());

    let dir = tempfile::tempdir().expect("tempdir");
    let (bam, fasta) = write_inputs(dir.path(), &reads);

    let (all, plain_kept) = fetch_plain(&bam);
    let expected: Vec<Snapshot> = all.iter().filter(|s| s.accepts(&pred)).cloned().collect();

    let (got, counts) = fetch_filtered(&bam, pred.clone(), Gate::Raw);
    assert_eq!(got, expected, "filtered fetch vs the unfiltered fetch restricted");

    // `fetched` counts what the region query produced and must not move with
    // the filter; `kept` counts what survived it.
    assert_eq!(counts.fetched, plain_kept, "fetched must be filter-independent");
    assert_eq!(counts.kept, expected.len(), "kept must be the accepted count");

    assert_eq!(
        fetch_through_readers(&bam, &fasta, pred),
        expected,
        "the customize value attached at open time must filter the same way"
    );

    report_selectivity(&tc, expected.len(), all.len());
    if reads.iter().any(|r| r.flags & 0x4 != 0) {
        tc.event("has a placed-unmapped read");
    }
}

// r[verify unified.fetch_into_customized]
// r[verify record_store.pre_filter.rollback]
// r[verify record_store.filter_raw]
/// The two rejection hooks are a performance choice, not a semantic one: the
/// same decision taken in `filter_raw` (before any slab is written) and in
/// `filter` (after, with a rollback) must leave the reader's store holding the
/// same records with the same bytes.
#[hegel::test]
fn both_filter_gates_leave_the_same_store(tc: TestCase) {
    let reads = tc.draw(arb_reads().print_as_debug());
    let pred = tc.draw(arb_predicate().print_as_debug());

    let dir = tempfile::tempdir().expect("tempdir");
    let (bam, _fasta) = write_inputs(dir.path(), &reads);

    let (early, early_counts) = fetch_filtered(&bam, pred.clone(), Gate::Raw);
    let (late, late_counts) = fetch_filtered(&bam, pred, Gate::Late);

    assert_eq!(early, late, "rejecting early vs rejecting late");
    assert_eq!(early_counts, late_counts, "counts");

    report_selectivity(&tc, early.len(), early_counts.fetched);
}

// r[verify unified.fetch_into_customized]
// r[verify bam.reader.unmapped_skipped]
/// An always-true filter must be invisible: the store it leaves is what plain
/// `fetch_into` leaves, and every fetched record is a kept one.
///
/// This is also the one property with an absolute oracle rather than a
/// relative one: the generated reads all fit inside the contig, so the store
/// must hold every one of them — including the placed-unmapped ones, which
/// the reader is not allowed to drop on its own account.
#[hegel::test]
fn an_always_true_filter_matches_plain_fetch_into(tc: TestCase) {
    let reads = tc.draw(arb_reads().print_as_debug());

    let dir = tempfile::tempdir().expect("tempdir");
    let (bam, fasta) = write_inputs(dir.path(), &reads);

    let (plain, plain_kept) = fetch_plain(&bam);
    let (customized, counts) = fetch_filtered(&bam, Predicate::always(), Gate::Raw);

    assert_eq!(customized, plain, "an always-true filter vs fetch_into");
    assert_eq!(counts.kept, counts.fetched, "nothing was rejected");
    assert_eq!(counts.kept, plain_kept, "fetch_into's return is the kept count");
    assert_eq!(
        fetch_through_readers(&bam, &fasta, Predicate::always()),
        plain,
        "through Readers too"
    );

    assert_eq!(plain.len(), reads.len(), "every generated read fits inside the contig");
    assert_eq!(
        plain.iter().filter(|s| s.flags & 0x4 != 0).count(),
        reads.iter().filter(|r| r.flags & 0x4 != 0).count(),
        "a placed-unmapped read reaches the store; the reader is not a second gate"
    );
    tc.event_value("records", plain.len() as f64);
}

// r[verify unified.fetch_into_customized]
// r[verify record_store.pre_filter.rollback]
/// An always-false filter leaves an empty store that is still a working store:
/// nothing to iterate, nothing to pile up, and the next fetch into it refills
/// it exactly as a fresh store would be filled.
///
/// That last part is the one worth having, and it is why the rejection goes
/// through the *late* gate: every record is written into the slabs and then
/// rolled back, so the store that comes out is empty only because a few dozen
/// rollbacks each undid exactly their own push. The only place a caller can
/// observe one that undid too little is the *next* fetch —
/// `fetch_into_customized` reuses the store it is given, which is the whole
/// point of `Readers` holding one.
#[hegel::test]
fn an_always_false_filter_empties_the_store_without_breaking_it(tc: TestCase) {
    let reads = tc.draw(arb_reads().print_as_debug());

    let dir = tempfile::tempdir().expect("tempdir");
    let (bam, _fasta) = write_inputs(dir.path(), &reads);

    let mut reader = IndexedBamReader::open(&bam).expect("open");
    let tid = reader.header().tid(CONTIG.0).expect("contig");
    let (start, end) = whole_contig();
    let mut store = RecordStore::new();

    let counts = reader
        .fetch_into_customized(
            tid,
            start,
            end,
            &mut store,
            &mut ClauseFilter::late(Predicate::never()),
        )
        .expect("fetch_into_customized");

    assert!(counts.fetched > 0, "the region holds records");
    assert_eq!(counts.kept, 0, "an always-false filter keeps nothing");
    assert_eq!(store.len(), 0);
    assert!(store.is_empty());
    assert_eq!(store.indices().count(), 0);
    assert_eq!(store.records().len(), 0);

    // Refilling the rolled-back store must give what a fresh store gives.
    reader.fetch_into(tid, start, end, &mut store).expect("fetch_into");
    let (fresh, _) = fetch_plain(&bam);
    assert_eq!(snapshot(&store), fresh, "a reused store must refill like a fresh one");

    // The derived structures are consistent with an empty store, not merely
    // zero-length: preparing one yields no mate links and no columns.
    // `prepare_for_pileup` consumes the store, so this runs on its own.
    let mut emptied = RecordStore::new();
    reader
        .fetch_into_customized(
            tid,
            start,
            end,
            &mut emptied,
            &mut ClauseFilter::late(Predicate::never()),
        )
        .expect("fetch_into_customized");
    let prepared = emptied.prepare_for_pileup();
    assert_eq!(prepared.stats, seqair::bam::record_store::MateLinkStats::default());
    let mut engine = seqair::bam::pileup::PileupEngine::new(
        prepared.input,
        Pos0::ZERO,
        Pos0::new(99u32).unwrap(),
    );
    assert!(engine.pileups().is_none(), "an empty store has no columns");

    tc.event_value("rolled back", counts.fetched as f64);
}

// r[verify unified.fetch_into_customized]
/// Filtering composes the way the algebra says it does: fetching with `a ∧ b`
/// gives exactly the records that both `a` and `b` accept on their own, and
/// running the same fetch twice gives the same answer both times.
#[hegel::test]
fn a_conjunction_fetches_the_intersection_of_its_parts(tc: TestCase) {
    let reads = tc.draw(arb_reads().print_as_debug());
    let a = tc.draw(arb_clauses());
    let b = tc.draw(arb_clauses());

    let dir = tempfile::tempdir().expect("tempdir");
    let (bam, _fasta) = write_inputs(dir.path(), &reads);

    let only_a = Predicate { mode: Mode::All, negate: false, clauses: a.clone() };
    let only_b = Predicate { mode: Mode::All, negate: false, clauses: b.clone() };

    let (from_a, _) = fetch_filtered(&bam, only_a.clone(), Gate::Raw);
    let (from_b, _) = fetch_filtered(&bam, only_b, Gate::Raw);
    let (from_both, _) = fetch_filtered(&bam, Predicate::and(&a, &b), Gate::Raw);

    // `from_a` is already in store order, so intersecting by membership
    // preserves that order.
    let expected: Vec<Snapshot> = from_a.iter().filter(|s| from_b.contains(s)).cloned().collect();
    assert_eq!(from_both, expected, "a ∧ b vs the intersection of a and b");

    // The same fetch, twice: `fetch_into_customized` clears the store it is
    // given, so a second run must not see the first one.
    let (again, _) = fetch_filtered(&bam, only_a, Gate::Raw);
    assert_eq!(again, from_a, "repeating a filtered fetch");

    tc.event_value("a", from_a.len() as f64);
    tc.event_value("b", from_b.len() as f64);
    tc.event_value("a and b", from_both.len() as f64);
    if from_both.len() < from_a.len() && from_both.len() < from_b.len() {
        tc.event("the conjunction cut further than either part");
    }
}
